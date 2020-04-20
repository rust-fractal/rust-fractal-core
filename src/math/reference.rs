use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed};

pub struct Reference {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub z_reference: Vec<ComplexFixed<f64>>,
    pub z_tolerance: Vec<f64>
}

impl Reference {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize) -> Reference {
        let z_fixed = to_fixed(&z);

        // 1e-6 is the threshold for pauldelbrot's criterion
        Reference {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            z_reference: vec![z_fixed],
            z_tolerance: vec![1e-6 * z_fixed.norm_sqr()]
        }
    }

    pub fn step(&mut self) -> bool {
        self.z = self.z.clone().square() + self.c.clone();
        self.current_iteration += 1;
        let z_fixed = to_fixed(&self.z);
        self.z_reference.push(z_fixed);
        self.z_tolerance.push(1e-6 * z_fixed.norm_sqr());
        z_fixed.norm_sqr() <= 4.0
    }

    pub fn run(&mut self) -> bool {
        while self.current_iteration <= self.maximum_iteration {
            if !self.step() {
                break;
            }
        };
        self.current_iteration == self.maximum_iteration
    }
}

