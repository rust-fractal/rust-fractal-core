use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed};

pub struct Reference {
    pub start_iteration: i32,
    pub current_iteration: i32,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub z_reference: Vec<ComplexFixed<f64>>,
    pub z_tolerance: Vec<f64>
}

impl Reference {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, start_iteration: i32) -> Reference {
        let z_fixed = to_fixed(&z);

        Reference {
            start_iteration,
            current_iteration: start_iteration,
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

    pub fn iterate(&mut self, max_iteration: i32) -> bool {
        while self.current_iteration <= max_iteration {
            if !self.step() {
                break;
            }
        };
        self.current_iteration == max_iteration
    }
}

