use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed, to_fixed_exp};
use float_extended::complex_extended::ComplexExtended;
use crate::math::reference::Reference;
use rug::Float;
use float_extended::float_extended::FloatExtended;

pub struct SeriesApproximation2 {
    pub current_iteration: usize,
    maximum_iteration: usize,
    delta_pixel: f64,
    z: ComplexArbitrary,
    c: ComplexArbitrary,
    pub order: usize,
    coefficients: Vec<ComplexExtended>,
    current_probes: Vec<ComplexFixed<f64>>,
    original_probes: Vec<ComplexFixed<f64>>
}

impl SeriesApproximation2 {
    pub fn new(c: ComplexArbitrary, order: usize, maximum_iteration: usize, delta_pixel: f64, delta_top_left: ComplexExtended) -> Self {
        assert!(order >= 1);

        let mut coefficients = vec![ComplexExtended::new2(0.0, 0.0, 0); order as usize + 1];

        let temp = to_fixed_exp(&c);
        coefficients[0] = ComplexExtended::new(temp.0, temp.1);
        coefficients[1] = ComplexExtended::new2(1.0, 0.0, 0);

        let probes = vec![
            ComplexExtended::new2(delta_top_left.mantissa.re, delta_top_left.mantissa.im, delta_top_left.exponent),
            ComplexExtended::new2(delta_top_left.mantissa.re, delta_top_left.mantissa.im * -1.0, delta_top_left.exponent),
            ComplexExtended::new2(delta_top_left.mantissa.re * -1.0, delta_top_left.mantissa.im, delta_top_left.exponent),
            ComplexExtended::new2(delta_top_left.mantissa.re * -1.0, delta_top_left.mantissa.im * -1.0, delta_top_left.exponent),
        ];

        // The current iteration is set to 1 as we set z = c
        SeriesApproximation2 {
            current_iteration: 1,
            maximum_iteration,
            delta_pixel,
            z: c.clone(),
            c,
            order,
            coefficients,
            current_probes: probes.clone(),
            original_probes: probes
        }
    }

    pub fn step(&mut self) -> bool {
        let z_next = self.z.clone().square() + self.c.clone();
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

            next_coefficients[k] = sum;
        }

        for i in 0..self.original_probes.len() {
            // step the probe points using perturbation
            let temp = to_fixed_exp(&self.z);
            self.current_probes[i] = temp.0 * 2.0f64.powi(temp.1) * self.current_probes[i] * 2.0 + self.current_probes[i] * self.current_probes[i] + self.original_probes[i];

            // get the new approximations
            let mut original_probe_n = ComplexExtended::new(self.original_probes[i], 0);
            let mut series_probe = ComplexExtended::new2(0.0, 0.0, 0);

            for k in 1..=self.order {
                series_probe += next_coefficients[k] * original_probe_n;
                original_probe_n *= ComplexExtended::new(self.original_probes[i], 0);
            };

            let mut original_probe_derivative_n = ComplexExtended::new2(1.0, 0.0, 0);
            let mut derivative_probe = ComplexExtended::new2(0.0, 0.0, 0);

            for k in 1..=self.order {
                derivative_probe += next_coefficients[k] * original_probe_derivative_n * k as f64;
                original_probe_derivative_n *= ComplexExtended::new(self.original_probes[i], 0);
            };

            let mut relative_error = (ComplexExtended::new(self.current_probes[i], 0) - series_probe).norm().to_float();
            let mut derivative = derivative_probe.norm().to_float();

            // Check to make sure that the derivative is greater than or equal to 1
            if derivative < 1.0 {
                derivative = 1.0;
            }

            // Check that the error over the derivative is less than the pixel spacing
            if relative_error / derivative > self.delta_pixel {
                return false;
            }
        }

        // If the approximation is valid, continue
        self.current_iteration += 1;
        self.z = z_next;
        self.coefficients = next_coefficients;
        true
    }

    pub fn run(&mut self) {
        // Can be changed later into a better loop - this function could also return some more information
        while self.step() && self.current_iteration < self.maximum_iteration {
            continue;
        }
    }

    // Get the current reference, and the current number of iterations done
    pub fn get_reference(&self) -> Reference {
        Reference::new(self.z.clone(), self.c.clone(), self.current_iteration, self.maximum_iteration)
    }

    pub fn evaluate(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
        let mut original_point_n = ComplexExtended::new(point_delta, 0);
        let mut approximation = ComplexExtended::new2(0.0, 0.0, 0);

        for k in 1..=self.order {
            approximation += self.coefficients[k] * original_point_n;
            original_point_n *= ComplexExtended::new(point_delta, 0);
        };

        approximation.to_float()
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
}