use crate::util::{ComplexArbitrary, to_extended};
use crate::util::complex_extended::ComplexExtended;
use crate::math::reference::Reference;
use rug::Float;
use crate::util::float_extended::FloatExtended;
use std::mem::swap;

pub struct SeriesApproximation {
    maximum_iteration: usize,
    delta_pixel_square: FloatExtended,
    pub order: usize,
    coefficients: Vec<Vec<ComplexExtended>>,
    original_probes: Vec<ComplexExtended>,
    current_probes: Vec<ComplexExtended>,
    approximation_probes: Vec<Vec<ComplexExtended>>,
    approximation_probes_derivative: Vec<Vec<ComplexExtended>>,
    delta_top_left: ComplexExtended,
    valid_coefficients: Vec<ComplexExtended>,
    pub valid_iteration: usize,
    center_reference: &Reference,
}

impl SeriesApproximation {
    pub fn new(center_reference: &Reference, order: usize, maximum_iteration: usize, delta_pixel_square: FloatExtended, delta_top_left: ComplexExtended) -> Self {
        let mut coefficients = vec![vec![ComplexExtended::new2(0.0, 0.0, 0); order as usize + 1]; 1];

        coefficients[0][0] = to_extended(&center_reference.c);
        coefficients[0][1] = ComplexExtended::new2(1.0, 0.0, 0);

        // The current iteration is set to 1 as we set z = c
        SeriesApproximation {
            maximum_iteration,
            delta_pixel_square,
            order,
            coefficients,
            original_probes: Vec::new(),
            current_probes: Vec::new(),
            approximation_probes: Vec::new(),
            approximation_probes_derivative: Vec::new(),
            delta_top_left,
            valid_coefficients: Vec::new(),
            valid_iteration: 1,
            center_reference
        }
    }

    pub fn generate_approximation(&mut self) {
        let add_value = ComplexExtended::new2(1.0, 0.0, 0);

        // Can be changed later into a better loop - this function could also return some more information
        for i in 1..self.maximum_iteration {
            let coefficients = &self.coefficients[i - 1];
            let mut next_coefficients = vec![ComplexExtended::new2(0.0, 0.0, 0); self.order as usize + 1];

            // This is checking if the approximation can step forward so takes the next iteration
            next_coefficients[0] = self.center_reference.approximation_data[i];
            next_coefficients[1] = coefficients[0] * coefficients[1] * 2.0 + add_value;
            next_coefficients[0].reduce();
            next_coefficients[1].reduce();

            // Calculate the new coefficents
            for k in 2..=self.order {
                let mut sum = coefficients[0] * coefficients[k];
                sum.reduce();

                for j in 1..=((k - 1) / 2) {
                    sum.reduce();
                    sum += coefficients[j] * coefficients[k - j];
                }
                sum *= 2.0;

                // If even, we include the mid term as well
                if k % 2 == 0 {
                    sum += coefficients[k / 2] * coefficients[k / 2];
                }

                sum.reduce();
                next_coefficients[k] = sum;
            }

            self.coefficients.push(next_coefficients);
        }
    }

    pub fn check_approximation(&mut self) {
        self.add_probe(ComplexExtended::new2(self.delta_top_left.mantissa.re, self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        self.add_probe(ComplexExtended::new2(self.delta_top_left.mantissa.re, -self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        self.add_probe(ComplexExtended::new2(-self.delta_top_left.mantissa.re, self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        self.add_probe(ComplexExtended::new2(-self.delta_top_left.mantissa.re, -self.delta_top_left.mantissa.im, self.delta_top_left.exponent));

        // Middle edges
        // self.add_probe(ComplexExtended::new2(self.delta_top_left.mantissa.re, 0.0, self.delta_top_left.exponent));
        // self.add_probe(ComplexExtended::new2(-self.delta_top_left.mantissa.re, 0.0, self.delta_top_left.exponent));
        // self.add_probe(ComplexExtended::new2(0.0, self.delta_top_left.mantissa.im, self.delta_top_left.exponent));
        // self.add_probe(ComplexExtended::new2(0.0, -self.delta_top_left.mantissa.im, self.delta_top_left.exponent));

        self.valid_iteration = 1;

        while self.valid_iteration < self.maximum_iteration {
            let coefficients = &self.coefficients[self.valid_iteration - 1];
            let next_coefficients = &self.coefficients[self.valid_iteration];

            // this section takes about half of the time
            for i in 0..self.original_probes.len() {
                // step the probe points using perturbation
                self.current_probes[i] = self.current_probes[i] * (coefficients[0] * 2.0 + self.current_probes[i]);
                self.current_probes[i] += self.original_probes[i];

                // TODO this probably does not need to be done every iteration
                self.current_probes[i].reduce();

                // get the new approximations
                let mut series_probe = next_coefficients[1] * self.approximation_probes[i][0];
                let mut derivative_probe = next_coefficients[1] * self.approximation_probes_derivative[i][0];

                for k in 2..=self.order {
                    series_probe += next_coefficients[k] * self.approximation_probes[i][k - 1];
                    derivative_probe += next_coefficients[k] * self.approximation_probes_derivative[i][k - 1];
                };

                let relative_error = (self.current_probes[i] - series_probe).norm_square();
                let mut derivative = derivative_probe.norm_square();

                // Check to make sure that the derivative is greater than or equal to 1
                if derivative.to_float() < 1.0 {
                    derivative.mantissa = 1.0;
                    derivative.exponent = 0;
                }

                // Check that the error over the derivative is less than the pixel spacing
                if relative_error / derivative > self.delta_pixel_square {
                    self.valid_coefficients = coefficients.clone();
                    return;
                }
            }

            // If the approximation is valid, continue
            self.valid_iteration += 1;
        }
    }

    pub fn add_probe(&mut self, delta_probe: ComplexExtended) {
        // here we will need to check to make sure we are still at the first iteration, or use perturbation to go forward
        self.original_probes.push(delta_probe);
        self.current_probes.push(delta_probe);

        let mut current_value = delta_probe;

        let mut delta_n = Vec::with_capacity(self.order + 1);
        let mut delta_derivative_n = Vec::with_capacity(self.order + 1);

        // The first element will be 1, in order for the derivative to be calculated
        delta_n.push(current_value);
        delta_derivative_n.push(ComplexExtended::new2(1.0, 0.0, 0));

        for i in 1..=self.order {
            delta_derivative_n.push(current_value * (i + 1) as f64);
            current_value *= delta_probe;
            delta_n.push(current_value);
        }

        self.approximation_probes.push(delta_n);
        self.approximation_probes_derivative.push(delta_derivative_n);
    }

    // Get the current reference, and the current number of iterations done
    pub fn get_reference(&self, reference_delta: ComplexExtended) -> Reference {
        let mut reference_c = self.center_reference.c.clone();
        let temp = Float::with_val(self.center_reference.c.real().prec(), reference_delta.exponent).exp2();
        let temp2 = Float::with_val(self.center_reference.c.real().prec(), reference_delta.mantissa.re);
        let temp3 = Float::with_val(self.center_reference.c.real().prec(), reference_delta.mantissa.im);

        *reference_c.mut_real() += &temp2 * &temp;
        *reference_c.mut_imag() += &temp3 * &temp;

        // let mut reference_z = self.center_reference.approximation_data[self.valid_iteration].clone();
        let mut reference_z = self.center_reference.approximation_data[self.valid_iteration]

        let temp4 = self.evaluate(reference_delta);
        let temp = Float::with_val(self.c.real().prec(), temp4.exponent).exp2();
        let temp2 = Float::with_val(self.c.real().prec(), temp4.mantissa.re);
        let temp3 = Float::with_val(self.c.real().prec(), temp4.mantissa.im);

        *reference_z.mut_real() += &temp2 * &temp;
        *reference_z.mut_imag() += &temp3 * &temp;

        

        Reference::new(reference_z, reference_c, self.current_iteration, self.maximum_iteration)
    }

    pub fn evaluate(&self, point_delta: ComplexExtended) -> ComplexExtended {
        // 1907 ms packing opus 4K
        // Horner's rule
        let mut approximation = self.valid_coefficients[self.order];

        for k in (1..=(self.order - 1)).rev() {
            approximation *= point_delta;
            approximation += self.valid_coefficients[k];
        }

        approximation *= point_delta;
        approximation.reduce();
        approximation
    }

    // pub fn evaluate_derivative(&self, point_delta: ComplexExtended) -> FloatExtended {
    //     let mut original_point_derivative_n = ComplexExtended::new(1.0, 0, 0.0, 0);
    //     let mut approximation_derivative = ComplexExtended::new(0.0, 0, 0.0, 0);
    
    //     for k in 1..=self.order {
    //         approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
    //         original_point_derivative_n *= ComplexExtended::new(point_delta.re, 0, point_delta.im, 0);
    //     };
    
    //     approximation_derivative.to_float()


    //     let mut approximation_derivative = self.coefficients[self.order] * self.order as f64;

    //     for k in (1..=(self.order - 1)).rev() {
    //         approximation *= point_delta;
    //         approximation += self.coefficients[k] * k as f64;
    //     }

    //     approximation *= point_delta;
    //     approximation.reduce();
    //     approximation


    // }
}