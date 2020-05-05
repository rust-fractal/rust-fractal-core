use crate::util::{ComplexArbitrary, to_fixed_exp};
use float_extended::complex_extended::ComplexExtended;
use crate::math::reference_extended::ReferenceExtended;
use rug::{Float, Assign};
use float_extended::float_extended::FloatExtended;

pub struct SeriesApproximationExtended {
    pub current_iteration: usize,
    maximum_iteration: usize,
    delta_pixel: FloatExtended,
    z: ComplexArbitrary,
    c: ComplexArbitrary,
    pub order: usize,
    coefficients: Vec<ComplexExtended>,
    current_probes: Vec<ComplexExtended>,
    original_probes: Vec<ComplexExtended>,
    approximation_probes: Vec<Vec<ComplexExtended>>,
    approximation_probes_derivative: Vec<Vec<ComplexExtended>>,
    delta_top_left: ComplexExtended
}

impl SeriesApproximationExtended {
    pub fn new(c: ComplexArbitrary, order: usize, maximum_iteration: usize, delta_pixel: FloatExtended, delta_top_left: ComplexExtended) -> Self {
        assert!(order >= 1);

        let mut coefficients = vec![ComplexExtended::new2(0.0, 0.0, 0); order as usize + 1];

        let temp = to_fixed_exp(&c);
        coefficients[0] = ComplexExtended::new(temp.0, temp.1);
        coefficients[1] = ComplexExtended::new2(1.0, 0.0, 0);

        // The current iteration is set to 1 as we set z = c
        SeriesApproximationExtended {
            current_iteration: 1,
            maximum_iteration,
            delta_pixel,
            z: c.clone(),
            c,
            order,
            coefficients,
            current_probes: Vec::new(),
            original_probes: Vec::new(),
            approximation_probes: Vec::new(),
            approximation_probes_derivative: Vec::new(),
            delta_top_left
        }
    }

    pub fn step(&mut self) -> bool {
        let z_next = self.z.clone().square() + &self.c;

        let mut next_coefficients = vec![ComplexExtended::new2(0.0, 0.0, 0); self.order as usize + 1];

        let temp = to_fixed_exp(&z_next);
        next_coefficients[0] = ComplexExtended::new(temp.0, temp.1);
        next_coefficients[1] = self.coefficients[0] * self.coefficients[1] * 2.0 + ComplexExtended::new2(1.0, 0.0, 0);

        // Calculate the new coefficents
        for k in 2..=self.order {
            let mut sum = ComplexExtended::new2(0.0, 0.0, 0);

            for j in 0..=((k - 1) / 2) {
                sum += self.coefficients[j] * self.coefficients[k - j];
            }
            sum *= 2.0;

            // If even, we include the mid term as well
            if k % 2 == 0 {
                sum += self.coefficients[k / 2] * self.coefficients[k / 2];
            }

            sum.reduce();
            next_coefficients[k] = sum;
        }

        let temp = to_fixed_exp(&self.z);
        let temp2 = ComplexExtended::new(temp.0, temp.1 + 1);

        // this section takes about half of the time
        for i in 0..self.original_probes.len() {
            // step the probe points using perturbation
            self.current_probes[i] = temp2 * self.current_probes[i] + self.current_probes[i] * self.current_probes[i] + self.original_probes[i];
            self.current_probes[i].reduce();

            // get the new approximations
            let mut series_probe = ComplexExtended::new2(0.0, 0.0, 0);
            let mut derivative_probe = ComplexExtended::new2(0.0, 0.0, 0);

            for k in 1..=self.order {
                series_probe += next_coefficients[k] * self.approximation_probes[i][k - 1];
                derivative_probe += next_coefficients[k] * self.approximation_probes_derivative[i][k - 1] * k as f64;

                if k % 16 == 0 {
                    series_probe.reduce();
                    derivative_probe.reduce();
                }
            };

            let relative_error = (self.current_probes[i] - series_probe).norm();
            let mut derivative = derivative_probe.norm();

            // Check to make sure that the derivative is greater than or equal to 1
            if derivative.to_float() < 1.0 {
                derivative = FloatExtended::new(1.0, 0);
            }

            // Check that the error over the derivative is less than the pixel spacing
            if relative_error / derivative > self.delta_pixel {
                return false;
            }
        }

        // If the approximation is valid, continue
        self.current_iteration += 1;
        self.z.assign(&z_next);
        self.coefficients = next_coefficients;
        true
    }

    pub fn run(&mut self) {
        self.add_probe(ComplexExtended::new2(self.delta_top_left.mantissa.re, self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        self.add_probe(ComplexExtended::new2(self.delta_top_left.mantissa.re, -self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        self.add_probe(ComplexExtended::new2(-self.delta_top_left.mantissa.re, self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        self.add_probe(ComplexExtended::new2(-self.delta_top_left.mantissa.re, -self.delta_top_left.mantissa.im, self.delta_top_left.exponent));

        // Can be changed later into a better loop - this function could also return some more information
        while self.step() && self.current_iteration < self.maximum_iteration {
            continue;
        }
    }

    // Get the current reference, and the current number of iterations done
    pub fn get_reference(&self, reference_delta: ComplexExtended) -> ReferenceExtended {
        let mut reference_c = self.c.clone();
        let temp = Float::with_val(self.c.real().prec(), reference_delta.exponent).exp2();
        let temp2 = Float::with_val(self.c.real().prec(), reference_delta.mantissa.re);
        let temp3 = Float::with_val(self.c.real().prec(), reference_delta.mantissa.im);

        *reference_c.mut_real() += &temp2 * &temp;
        *reference_c.mut_imag() += &temp3 * &temp;

        let mut reference_z = self.z.clone();
        let temp4 = self.evaluate(reference_delta);
        let temp = Float::with_val(self.c.real().prec(), temp4.exponent).exp2();
        let temp2 = Float::with_val(self.c.real().prec(), temp4.mantissa.re);
        let temp3 = Float::with_val(self.c.real().prec(), temp4.mantissa.im);

        *reference_z.mut_real() += &temp2 * &temp;
        *reference_z.mut_imag() += &temp3 * &temp;

        ReferenceExtended::new(reference_z, reference_c, self.current_iteration, self.maximum_iteration)
    }

    pub fn evaluate(&self, point_delta: ComplexExtended) -> ComplexExtended {
        let mut approximation = self.coefficients[1] * point_delta;
        let mut original_point_n = point_delta * point_delta;

        for k in 2..=self.order {
            approximation += self.coefficients[k] * original_point_n;
            original_point_n *= point_delta;

            if k % 16 == 0 {
                approximation.reduce();
                original_point_n.reduce();
            }
        };

        approximation
    }

    // pub fn evaluate_derivative(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
    //     let mut original_point_derivative_n = ComplexExtended::new(1.0, 0, 0.0, 0);
    //     let mut approximation_derivative = ComplexExtended::new(0.0, 0, 0.0, 0);
    //
    //     for k in 1..=self.order {
    //         approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
    //         original_point_derivative_n *= ComplexExtended::new(point_delta.re, 0, point_delta.im, 0);
    //     };
    //
    //     approximation_derivative.to_float()
    // }

    pub fn add_probe(&mut self, delta_probe: ComplexExtended) {
        // here we will need to check to make sure we are still at the first iteration, or use perturbation to go forward
        self.original_probes.push(delta_probe);
        self.current_probes.push(delta_probe);

        let mut delta_probe_n = Vec::with_capacity(self.order + 1);
        let mut delta_probe_n_derivative = Vec::with_capacity(self.order + 1);

        // The first element will be 1, in order for the derivative to be calculated
        delta_probe_n.push(delta_probe);
        delta_probe_n_derivative.push(ComplexExtended::new2(1.0, 0.0, 0));

        for i in 1..=self.order {
            let mut delta_probe1 = delta_probe_n[i - 1] * delta_probe;
            let mut delta_probe2 = delta_probe_n_derivative[i - 1] * delta_probe;
            delta_probe1.reduce();
            delta_probe2.reduce();
            delta_probe_n.push(delta_probe1);
            delta_probe_n_derivative.push(delta_probe2);
        }

        self.approximation_probes.push(delta_probe_n);
        self.approximation_probes_derivative.push(delta_probe_n_derivative);
    }
}