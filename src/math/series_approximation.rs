use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed};
use float_extended::complex_extended::ComplexExtended;
use crate::math::reference::Reference;
use rug::Float;
use float_extended::float_extended::FloatExtended;

pub struct SeriesApproximation {
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

impl SeriesApproximation {
    pub fn new(c: ComplexArbitrary, order: usize, maximum_iteration: usize, delta_pixel: f64, delta_top_left: ComplexFixed<f64>) -> Self {
        assert!(order >= 1);

        let mut coefficents = vec![ComplexExtended::new(0.0, 0, 0.0, 0); order as usize + 1];
        coefficents[0] = ComplexExtended::new(c.real().to_f64(), 0, c.imag().to_f64(), 0);
        coefficents[1] = ComplexExtended::new(1.0, 0, 0.0, 0);

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
            original_probes: probes
        }
    }

    pub fn step(&mut self) -> bool {
        let z_next = self.z.clone().square() + self.c.clone();
        let mut next_coefficents = vec![ComplexExtended::new(0.0, 0, 0.0, 0); self.order as usize + 1];
        next_coefficents[0] = ComplexExtended::new(z_next.real().to_f64(), 0, z_next.imag().to_f64(), 0);
        next_coefficents[1] = 2.0 * self.coefficients[0] * self.coefficients[1] + 1.0;

        // Calculate the new coefficents
        for k in 2..=self.order {
            let mut sum = ComplexExtended::new(0.0, 0, 0.0, 0);

            for j in 0..=((k - 1) / 2) {
                sum += self.coefficients[j] * self.coefficients[k - j];
            }
            sum *= 2.0;

            // If even, we include the mid term as well
            if k % 2 == 0 {
                sum += self.coefficients[k / 2] * self.coefficients[k / 2];
            }

            next_coefficents[k] = sum;
        }

        for i in 0..self.original_probes.len() {
            // step the probe points using perturbation
            self.current_probes[i] = 2.0 * ComplexFixed::new(self.z.real().to_f64(), self.z.imag().to_f64()) * self.current_probes[i] + self.current_probes[i] * self.current_probes[i] + self.original_probes[i];

            // get the new approximations
            let mut original_probe_n = ComplexExtended::new(self.original_probes[i].re, 0, self.original_probes[i].im, 0);
            let mut series_probe = ComplexExtended::new(0.0, 0, 0.0, 0);

            for k in 1..=self.order {
                series_probe += next_coefficents[k] * original_probe_n;
                original_probe_n *= ComplexExtended::new(self.original_probes[i].re, 0, self.original_probes[i].im, 0);
            };

            let mut original_probe_derivative_n = ComplexExtended::new(1.0, 0, 0.0, 0);
            let mut derivative_probe = ComplexExtended::new(0.0, 0, 0.0, 0);

            for k in 1..=self.order {
                derivative_probe += k as f64 * next_coefficents[k] * original_probe_derivative_n;
                original_probe_derivative_n *= ComplexExtended::new(self.original_probes[i].re, 0, self.original_probes[i].im, 0);
            };

            let mut relative_error = (ComplexExtended::new(self.current_probes[i].re, 0, self.current_probes[i].im, 0) - series_probe).norm().to_float();
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
        self.coefficients = next_coefficents;
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
        let mut original_point_n = ComplexExtended::new(point_delta.re, 0, point_delta.im, 0);
        let mut approximation = ComplexExtended::new(0.0, 0, 0.0, 0);

        for k in 1..=self.order {
            approximation += self.coefficients[k] * original_point_n;
            original_point_n *= ComplexExtended::new(point_delta.re, 0, point_delta.im, 0);
        };

        approximation.to_float()
    }

    pub fn evaluate_derivative(&self, point_delta: ComplexFixed<f64>) -> ComplexFixed<f64> {
        let mut original_point_derivative_n = ComplexExtended::new(1.0, 0, 0.0, 0);
        let mut approximation_derivative = ComplexExtended::new(0.0, 0, 0.0, 0);

        for k in 1..=self.order {
            approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
            original_point_derivative_n *= ComplexExtended::new(point_delta.re, 0, point_delta.im, 0);
        };

        approximation_derivative.to_float()
    }
}