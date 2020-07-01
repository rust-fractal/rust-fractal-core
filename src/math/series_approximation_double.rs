use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed};
use crate::math::reference_double::ReferenceDouble;
use std::mem::swap;

pub struct SeriesApproximationDouble {
    pub current_iteration: usize,
    maximum_iteration: usize,
    delta_pixel: f64,
    delta_top_left: ComplexFixed<f64>,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub order: usize,
    coefficients: Vec<ComplexFixed<f64>>,
    next_coefficients: Vec<ComplexFixed<f64>>,
    original_probes: Vec<ComplexFixed<f64>>,
    perturbation_probes: Vec<ComplexFixed<f64>>,
    approximation_probes: Vec<Vec<ComplexFixed<f64>>>,
    approximation_probes_derivative: Vec<Vec<ComplexFixed<f64>>>,
    delta_maximum: f64
}

impl SeriesApproximationDouble {
    pub fn new(c: ComplexArbitrary, order: usize, maximum_iteration: usize, delta_pixel: f64, delta_top_left: ComplexFixed<f64>) -> Self {
        let delta_maximum = delta_top_left.norm();

        // To avoid scaling issues we can put in delta / delta_maximum or something
        let mut coefficients = vec![ComplexFixed::new(0.0, 0.0); order as usize + 1];
        coefficients[0] = to_fixed(&c);

        // Set 1, 1, to the maximum delta
        coefficients[1] = ComplexFixed::new(delta_maximum, 0.0);

        // The current iteration is set to 1 as we set z = c
        SeriesApproximationDouble {
            current_iteration: 1,
            maximum_iteration,
            delta_pixel,
            delta_top_left,
            z: c.clone(),
            c,
            order,
            coefficients: coefficients.clone(),
            next_coefficients: coefficients,
            original_probes: Vec::new(),
            perturbation_probes: Vec::new(),
            approximation_probes: Vec::new(),
            approximation_probes_derivative: Vec::new(),
            delta_maximum,
        }
    }

    pub fn step(&mut self) -> bool {
        // This section takes 2328 ms / ~6000ms
        self.z.square_mut();
        self.z += &self.c;

        // Store the x_n in the first element of the array to simplfy things
        self.next_coefficients[0] = to_fixed(&self.z);
        self.next_coefficients[1] = 2.0 * self.coefficients[0] * self.coefficients[1] + self.delta_maximum;

        // Calculate the new coefficients
        for k in 2..=self.order {
            let mut sum = ComplexFixed::new(0.0, 0.0);

            // This adds all opposite pair multiples
            for j in 0..=((k - 1) / 2) {
                sum += self.coefficients[j] * self.coefficients[k - j];
            }
            sum *= 2.0;

            // If even, we include the mid term as well
            if k % 2 == 0 {
                sum += self.coefficients[k / 2] * self.coefficients[k / 2];
            }

            self.next_coefficients[k] = sum;
        }

        // This section is about 2000 ms

        for i in 0..self.perturbation_probes.len() {
            self.perturbation_probes[i] = 2.0 * self.coefficients[0] * self.perturbation_probes[i] + self.perturbation_probes[i] * self.perturbation_probes[i] + self.original_probes[i];
        }

        for i in 0..self.approximation_probes.len() {
            // Get the new approximations
            let mut series_probe = ComplexFixed::new(0.0, 0.0);
            let mut derivative_probe = ComplexFixed::new(0.0, 0.0);

            for k in 1..=self.order {
                series_probe += self.next_coefficients[k] * self.approximation_probes[i][k - 1];
                derivative_probe += (k + 1) as f64 * self.next_coefficients[k] * self.approximation_probes_derivative[i][k - 1];
            };

            let relative_error = (self.perturbation_probes[i] - series_probe).norm();
            let mut derivative = derivative_probe.norm();

            // Check to make sure that the derivative is greater than or equal to 1
            if derivative < 1.0 {
                derivative = 1.0;
            }

            // Check that the error over the derivative is less than the pixel spacing
            // TODO maybe roll back one iteration when the best has been found (might be more stable)
            if relative_error / derivative > self.delta_pixel {
                self.z -= &self.c;
                self.z.sqrt_mut();
                return false;
            }
        }

        // If the approximation is valid, continue
        self.current_iteration += 1;

        // Swap the two coefficients buffers
        swap(&mut self.coefficients, &mut self.next_coefficients);
        true
    }

    pub fn run(&mut self) {
        // Add the probe points to the series approximation
        self.add_probe(1.5 * ComplexFixed::new(self.delta_top_left.re, self.delta_top_left.im));
        self.add_probe(1.5 * ComplexFixed::new(self.delta_top_left.re, -self.delta_top_left.im));
        self.add_probe(1.5 * ComplexFixed::new(-self.delta_top_left.re, self.delta_top_left.im));
        self.add_probe(1.5 * ComplexFixed::new(-self.delta_top_left.re, -self.delta_top_left.im));

        // Can be changed later into a better loop - this function could also return some more information
        while self.step() && self.current_iteration < self.maximum_iteration {
            continue;
        }
    }

    pub fn add_probe(&mut self, delta_probe: ComplexFixed<f64>) {
        // here we will need to check to make sure we are still at the first iteration, or use perturbation to go forward
        self.original_probes.push(delta_probe);
        self.perturbation_probes.push(delta_probe);

        let delta_probe_scaled = delta_probe / self.delta_maximum;
        let mut delta_probe_n = Vec::with_capacity(self.order + 1);
        let mut delta_probe_n_derivative = Vec::with_capacity(self.order + 1);

        // The first element will be 1, in order for the derivative to be calculated
        delta_probe_n.push(delta_probe_scaled);
        delta_probe_n_derivative.push(ComplexFixed::new(1.0, 0.0) / self.delta_maximum);

        for i in 1..=self.order {
            delta_probe_n.push(delta_probe_n[i - 1] * delta_probe_scaled);
            delta_probe_n_derivative.push(delta_probe_n_derivative[i - 1] * delta_probe_scaled);
        }

        self.approximation_probes.push(delta_probe_n);
        self.approximation_probes_derivative.push(delta_probe_n_derivative);
    }

    // Get the current reference, and the current number of iterations done
    pub fn get_reference(&self, reference_delta: ComplexFixed<f64>) -> ReferenceDouble {
        let mut reference_c = self.c.clone();
        *reference_c.mut_real() = reference_c.real().clone() + reference_delta.re;
        *reference_c.mut_imag() = reference_c.imag().clone() + reference_delta.im;

        let mut reference_z = self.z.clone();
        let temp = self.evaluate(reference_delta);
        *reference_z.mut_real() = reference_z.real().clone() + temp.re;
        *reference_z.mut_imag() = reference_z.imag().clone() + temp.im;

        ReferenceDouble::new(reference_z, reference_c, self.current_iteration, self.maximum_iteration)
    }

    pub fn evaluate(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
        // could remove this divide
        let scaled_delta = point_delta / self.delta_maximum;

        let mut original_point_n = scaled_delta;
        let mut approximation = ComplexFixed::new(0.0, 0.0);

        for k in 1..=self.order {
            approximation += self.coefficients[k] * original_point_n;
            original_point_n *= scaled_delta
        };

        approximation
    }

    pub fn evaluate_derivative(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
        let scaled_delta = point_delta / self.delta_maximum;

        let mut original_point_derivative_n = ComplexFixed::new(1.0, 0.0) / self.delta_maximum;
        let mut approximation_derivative = ComplexFixed::new(0.0, 0.0);

        for k in 1..=self.order {
            approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
            original_point_derivative_n *= scaled_delta;
        };

        approximation_derivative
    }
}