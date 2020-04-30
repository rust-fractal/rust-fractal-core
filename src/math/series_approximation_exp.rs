use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed_exp};
use crate::math::reference::Reference;
use std::ops::{Mul, Add};
use float_extended::float_extended::FloatExtended;
use float_extended::util::FloatExp;

pub struct SeriesApproximationExp {
    pub current_iteration: usize,
    maximum_iteration: usize,
    delta_pixel: ComplexExtended,
    delta_top_left: ComplexFixed<f64>,
    delta_top_left_scale: i32,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub order: usize,
    coefficients: Vec<ComplexExtended>,
    original_probes: Vec<ComplexExtended>,
    perturbation_probes: Vec<ComplexExtended>,
    approximation_probes: Vec<Vec<ComplexExtended>>,
    approximation_probes_derivative: Vec<Vec<ComplexExtended>>,
    delta_maximum: f64
}

#[derive(Copy, Clone)]
pub struct ComplexExtended {
    mantissa: ComplexFixed<f64>,
    exponent: i32,
}

impl ComplexExtended {
    pub fn new(mantissa: ComplexFixed<f64>, exponent: i32) -> Self {
        ComplexExtended {
            mantissa,
            exponent,
        }
    }

    pub fn norm(&self) -> (f64, i32) {
        (self.mantissa.norm(), self.exponent)
    }

    pub fn rescale(&mut self) {
        let (temp_mantissa, added_exponent) = self.mantissa.re.frexp();
        self.mantissa.re = temp_mantissa;
        self.mantissa.im = self.mantissa.im.ldexp(-added_exponent);
        self.exponent += added_exponent;
    }
}

impl Mul for ComplexExtended {
    type Output = ComplexExtended;

    fn mul(self, rhs: Self) -> Self::Output {
        ComplexExtended {
            mantissa: self.mantissa * rhs.mantissa,
            exponent: self.exponent + rhs.exponent
        }
    }
}

impl Mul<f64> for ComplexExtended {
    type Output = ComplexExtended;

    fn mul(self, rhs: f64) -> Self::Output {
        ComplexExtended {
            mantissa: self.mantissa * rhs,
            exponent: self.exponent
        }
    }
}

impl Add for ComplexExtended {
    type Output = ComplexExtended;

    fn add(self, rhs: Self) -> Self::Output {
        if self.mantissa.norm_sqr() == 0.0 {
            rhs
        } else if rhs.mantissa.norm_sqr() == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == rhs.exponent {
                (self.mantissa + rhs.mantissa, self.exponent)
            } else if self.exponent > rhs.exponent {
                (self.mantissa + rhs.mantissa / 2.0f64.powi(self.exponent - rhs.exponent), self.exponent)
            } else {
                (rhs.mantissa + self.mantissa / 2.0f64.powi(rhs.exponent - self.exponent), rhs.exponent)
            };
            ComplexExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl SeriesApproximationExp {
    pub fn new(c: ComplexArbitrary, order: usize, maximum_iteration: usize, delta_pixel: ComplexExtended, delta_top_left: ComplexFixed<f64>, delta_top_left_scale: i32) -> Self {
        assert!(order >= 1);
        let delta_maximum = delta_top_left.norm();

        // To avoid scaling issues we can put in delta / delta_maximum or something
        let mut coefficients = vec![ComplexExtended::new(ComplexFixed::new(0.0, 0.0), 0); order as usize + 1];

        let temp = to_fixed_exp(&c);

        coefficients[0] = ComplexExtended::new(temp.0, temp.1);

        // Set 1, 1, to the maximum delta
        coefficients[1] = ComplexExtended::new(ComplexFixed::new(1.0, 0.0), 0);

        let mut temp2 = delta_pixel.clone();
        temp2.rescale();

        // The current iteration is set to 1 as we set z = c
        SeriesApproximationExp {
            current_iteration: 1,
            maximum_iteration,
            delta_pixel: temp2,
            delta_top_left,
            delta_top_left_scale,
            z: c.clone(),
            c,
            order,
            coefficients,
            original_probes: Vec::new(),
            perturbation_probes: Vec::new(),
            approximation_probes: Vec::new(),
            approximation_probes_derivative: Vec::new(),
            delta_maximum,
        }
    }

    pub fn step(&mut self) -> bool {
        let z_next = self.z.clone().square() + self.c.clone();
        let mut next_coefficients = vec![ComplexExtended::new(ComplexFixed::new(0.0, 0.0), 0); self.order as usize + 1];

        // Store the x_n in the first element of the array to simplfy things
        let temp = to_fixed_exp(&z_next);
        next_coefficients[0] = ComplexExtended::new(temp.0, temp.1);

        next_coefficients[1] = self.coefficients[0] * self.coefficients[1] * 2.0 + ComplexExtended::new(ComplexFixed::new(1.0, 0.0), 0);

        // Calculate the new coefficients
        for k in 2..=self.order {
            let mut sum = ComplexExtended::new(ComplexFixed::new(0.0, 0.0), 0);

            for j in 0..=((k - 1) / 2) {
                sum = sum + self.coefficients[j] * self.coefficients[k - j];
                sum.rescale();
            }
            sum = sum * 2.0;

            // If even, we include the mid term as well
            if k % 2 == 0 {
                sum = sum + self.coefficients[k / 2] * self.coefficients[k / 2];
                sum.rescale();
            }

            next_coefficients[k] = sum;
        }

        // Step the probe points using perturbation

        // TODO - here we could also use the optimised perturbation loop
        for i in 0..self.perturbation_probes.len() {
            self.perturbation_probes[i] = next_coefficients[0] * self.perturbation_probes[i] * 2.0 + self.perturbation_probes[i] * self.perturbation_probes[i] + self.original_probes[i];
            self.perturbation_probes[i].rescale();
        }

        for i in 0..self.approximation_probes.len() {
            // Get the new approximations
            let mut series_probe = ComplexExtended::new(ComplexFixed::new(0.0, 0.0), 0);
            let mut derivative_probe = ComplexExtended::new(ComplexFixed::new(0.0, 0.0), 0);

            for k in 1..=self.order {
                series_probe = series_probe + next_coefficients[k] * self.approximation_probes[i][k - 1];
                derivative_probe = derivative_probe + next_coefficients[k] * self.approximation_probes_derivative[i][k - 1] * (k + 1) as f64;
                series_probe.rescale();
                derivative_probe.rescale();
            };

            let relative_error = (self.perturbation_probes[i] + series_probe * -1.0).norm();
            let mut derivative = derivative_probe.norm();

            // Check to make sure that the derivative is greater than or equal to 1
            if derivative.0 + 2.0f64.powi(derivative.1) < 1.0 {
                derivative = (1.0, 0);
            }

            // Check that the error over the derivative is less than the pixel spacing

            let mut temp5 = (relative_error.0 / derivative.0, relative_error.1 - derivative.1);

            let (temp_mantissa, added_exponent) = temp5.0.frexp();
            temp5.0 = temp_mantissa;
            temp5.1 += added_exponent;

            // if temp5.1 > self.delta_pixel.exponent || (temp5.1 == self.delta_pixel.exponent && temp5.0 > self.delta_pixel.mantissa.norm()) {
            //     return false;
            // }
        }

        // If the approximation is valid, continue
        self.current_iteration += 1;
        self.z = z_next;
        self.coefficients = next_coefficients;
        self.current_iteration <= 10000
    }

    pub fn run(&mut self) {
        // Add the probe points to the series approximation
        self.add_probe(ComplexExtended {
            mantissa: ComplexFixed::new(self.delta_top_left.re, self.delta_top_left.im),
            exponent: self.delta_top_left_scale
        });
        self.add_probe(ComplexExtended {
            mantissa: ComplexFixed::new(self.delta_top_left.re, -self.delta_top_left.im),
            exponent: self.delta_top_left_scale
        });
        self.add_probe(ComplexExtended {
            mantissa: ComplexFixed::new(-self.delta_top_left.re, self.delta_top_left.im),
            exponent: self.delta_top_left_scale
        });
        self.add_probe(ComplexExtended {
            mantissa: ComplexFixed::new(-self.delta_top_left.re, -self.delta_top_left.im),
            exponent: self.delta_top_left_scale
        });

        // Can be changed later into a better loop - this function could also return some more information
        while self.step() && self.current_iteration < self.maximum_iteration {
            continue;
        }
    }

    pub fn add_probe(&mut self, delta_probe: ComplexExtended) {
        // here we will need to check to make sure we are still at the first iteration, or use perturbation to go forward
        self.original_probes.push(delta_probe);
        self.perturbation_probes.push(delta_probe);

        let mut delta_probe_n = Vec::with_capacity(self.order + 1);
        let mut delta_probe_n_derivative = Vec::with_capacity(self.order + 1);

        // The first element will be 1, in order for the derivative to be calculated
        delta_probe_n.push(delta_probe);
        delta_probe_n_derivative.push(ComplexExtended::new(ComplexFixed::new(1.0, 0.0), 0));

        for i in 1..=self.order {
            delta_probe_n.push(delta_probe_n[i - 1] * delta_probe);
            delta_probe_n_derivative.push(delta_probe_n_derivative[i - 1] * delta_probe);
            delta_probe_n.last_mut().unwrap().rescale();
            delta_probe_n_derivative.last_mut().unwrap().rescale();
        }

        self.approximation_probes.push(delta_probe_n);
        self.approximation_probes_derivative.push(delta_probe_n_derivative);
    }

    // Get the current reference, and the current number of iterations done
    pub fn get_reference(&self, reference_delta: ComplexExtended) -> Reference {
        // TODO - fix this for values that are not 0
        let mut reference_c = self.c.clone();
        *reference_c.mut_real() = reference_c.real().clone() + reference_delta.mantissa.re;
        *reference_c.mut_imag() = reference_c.imag().clone() + reference_delta.mantissa.im;

        let mut reference_z = self.z.clone();
        let temp = self.evaluate(reference_delta);
        *reference_z.mut_real() = reference_z.real().clone() + temp.mantissa.re;
        *reference_z.mut_imag() = reference_z.imag().clone() + temp.mantissa.im;

        Reference::new(reference_z, reference_c, self.current_iteration, self.maximum_iteration)
    }

    pub fn evaluate(&self, point_delta: ComplexExtended) -> ComplexExtended {
        // could remove this divide

        let mut original_point_n = point_delta;
        let mut approximation = ComplexExtended::new(ComplexFixed::new(0.0, 0.0), 0);

        for k in 1..=self.order {
            // for this, we need to take the approximation as floatexp
            approximation = approximation + self.coefficients[k] * original_point_n;
            approximation.rescale();
            original_point_n = original_point_n * point_delta;
            original_point_n.rescale();
        };

        approximation
    }

    // pub fn evaluate_derivative(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
    //     let scaled_delta = point_delta / self.delta_maximum;
    //
    //     let mut original_point_derivative_n = ComplexFixed::new(1.0, 0.0) / self.delta_maximum;
    //     let mut approximation_derivative = ComplexFixed::new(0.0, 0.0);
    //
    //     for k in 1..=self.order {
    //         approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
    //         original_point_derivative_n *= scaled_delta;
    //     };
    //
    //     approximation_derivative
    // }
}