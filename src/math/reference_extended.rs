use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed, to_fixed_exp};

pub struct ReferenceExtended {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub z_reference: Vec<(ComplexFixed<f64>, i32)>,
    pub z_tolerance: Vec<f64>
}

impl ReferenceExtended {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize) -> ReferenceExtended {
        let z_fixed = (to_fixed(&z), 0);

        // Stored in premultiplied form
        let z_n = if z_fixed.0.re <= 1e-128 && z_fixed.0.re >= -1e-128 {
            let temp = to_fixed_exp(&z);
            (temp.0 * 2.0, temp.1)
        } else {
            (z_fixed.0 * 2.0, 0)
        };

        // 1e-6 is the threshold for pauldelbrot's criterion
        ReferenceExtended {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            z_reference: vec![z_n],
            z_tolerance: vec![1e-6 * z_fixed.0.norm_sqr()]
        }
    }

    pub fn step(&mut self) -> bool {
        self.z = self.z.clone().square() + &self.c;
        self.current_iteration += 1;

        let z_fixed = (to_fixed(&self.z), 0);

        // Stored in premultiplied form
        let z_n = if z_fixed.0.re <= 1e-128 && z_fixed.0.re >= -1e-128 {
            let temp = to_fixed_exp(&self.z);
            (temp.0 * 2.0, temp.1)
        } else {
            (z_fixed.0 * 2.0, 0)
        };

        self.z_reference.push(z_n);
        self.z_tolerance.push(1e-6 * z_fixed.0.norm_sqr());

        // If the value is not small we do the escape check, otherwise it has not escaped
        // as we do the check for 65536 on the perturbation, we need this to be more than that squared
        z_fixed.0.norm_sqr() <= 1e256
    }


    pub fn run(&mut self) -> bool {
        while self.current_iteration < self.maximum_iteration {
            if !self.step() {
                break;
            }
        };
        self.current_iteration == self.maximum_iteration
    }
}

