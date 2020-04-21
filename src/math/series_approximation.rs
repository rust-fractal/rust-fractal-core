use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed};
use float_extended::complex_extended::ComplexExtended;
use crate::math::reference::Reference;
use rug::Float;
use float_extended::float_extended::FloatExtended;

pub struct SeriesApproximation {
    pub current_iteration: usize,
    maximum_iteration: usize,
    delta_pixel: f64,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub order: usize,
    coefficients: Vec<ComplexFixed<f64>>,
    current_probes: Vec<ComplexFixed<f64>>,
    original_probes: Vec<ComplexFixed<f64>>,
    bmax: f64
}

impl SeriesApproximation {
    pub fn new(c: ComplexArbitrary, order: usize, maximum_iteration: usize, delta_pixel: f64, delta_top_left: ComplexFixed<f64>) -> Self {
        assert!(order >= 1);

        // To avoid scaling issues we can put in b/bmax or something
        let mut coefficents = vec![ComplexFixed::new(0.0, 0.0); order as usize + 1];
        coefficents[0] = ComplexFixed::new(c.real().to_f64(), c.imag().to_f64());

        // Set 1, 1, to bmax
        coefficents[1] = ComplexFixed::new(delta_top_left.norm(), 0.0);

        // I think these should already be scaled
        let probes = vec![
            ComplexFixed::new(delta_top_left.re, delta_top_left.im),
            ComplexFixed::new(delta_top_left.re, delta_top_left.im * -1.0),
            ComplexFixed::new(delta_top_left.re * -1.0, delta_top_left.im),
            ComplexFixed::new(delta_top_left.re * -1.0, delta_top_left.im * -1.0),
        ];

        // The current iteration is set to 1 as we set z = c
        SeriesApproximation {
            current_iteration: 1,
            maximum_iteration,
            delta_pixel,
            z: c.clone(),
            c,
            order,
            coefficients: coefficents,
            current_probes: probes.clone(),
            original_probes: probes,
            bmax: delta_top_left.norm(),
        }
    }

    pub fn step(&mut self) -> bool {
        let z_next = self.z.clone().square() + self.c.clone();
        let mut next_coefficients = vec![ComplexFixed::new(0.0, 0.0); self.order as usize + 1];

        // Store the x_n in the first element of the array to simplfy things
        next_coefficients[0] = ComplexFixed::new(z_next.real().to_f64(), z_next.imag().to_f64());
        next_coefficients[1] = 2.0 * self.coefficients[0] * self.coefficients[1] + self.bmax;

        // Calculate the new coefficients
        for k in 2..=self.order {
            let mut sum = ComplexFixed::new(0.0, 0.0);

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
            self.current_probes[i] = 2.0 * ComplexFixed::new(self.z.real().to_f64(), self.z.imag().to_f64()) * self.current_probes[i] + self.current_probes[i] * self.current_probes[i] + self.original_probes[i];
            // self.derivative_probes[i] = 2.0 * ComplexFixed::new(self.z.real().to_f64(), self.z.imag().to_f64()) * self.derivative_probes[i] + 1.0;

            let scaled_delta = self.original_probes[i] / self.bmax;

            // get the new approximations
            let mut original_probe_n = scaled_delta;
            let mut series_probe = ComplexFixed::new(0.0, 0.0);
            let mut original_probe_derivative_n = ComplexFixed::new(1.0, 0.0) / self.bmax;
            let mut derivative_probe = ComplexFixed::new(0.0, 0.0);

            for k in 1..=self.order {
                series_probe += next_coefficients[k] * original_probe_n;
                original_probe_n *= scaled_delta;

                derivative_probe += k as f64 * next_coefficients[k] * original_probe_derivative_n;
                original_probe_derivative_n *= scaled_delta;
            };

            let mut relative_error = (self.current_probes[i] - series_probe).norm();
            let mut derivative = derivative_probe.norm();

            // Check to make sure that the derivative is greater than or equal to 1
            if derivative < 1.0 {
                derivative = 1.0;
            }

            // Check that the error over the derivative is less than the pixel spacing
            if relative_error / derivative > self.delta_pixel{
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
        let scaled_delta = point_delta / self.bmax;

        let mut original_point_n = scaled_delta;
        let mut approximation = ComplexFixed::new(0.0, 0.0);

        for k in 1..=self.order {
            approximation += self.coefficients[k] * original_point_n;
            original_point_n *= scaled_delta
        };

        approximation
    }

    pub fn evaluate_derivative(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
        let scaled_delta = point_delta / self.bmax;

        let mut original_point_derivative_n = ComplexFixed::new(1.0, 0.0) / self.bmax;
        let mut approximation_derivative = ComplexFixed::new(0.0, 0.0);

        for k in 1..=self.order {
            approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
            original_point_derivative_n *= scaled_delta;
        };

        approximation_derivative
    }
}