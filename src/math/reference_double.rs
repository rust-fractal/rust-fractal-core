use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed};

pub struct ReferenceDouble {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub z_reference: Vec<ComplexFixed<f64>>,
    pub z_tolerance: Vec<f64>
}

impl ReferenceDouble {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize) -> ReferenceDouble {
        let z_fixed = to_fixed(&z);

        // 1e-6 is the threshold for pauldelbrot's criterion
        ReferenceDouble {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            z_reference: vec![z_fixed],
            z_tolerance: vec![1e-6 * z_fixed.norm_sqr()]
        }
    }

    pub fn run(&mut self) -> bool {
        // Loop until the reference point escapes or the specified maximum is reached
        while self.current_iteration < self.maximum_iteration {
            // Square and add c
            self.z.square_mut();
            self.z += &self.c;
            self.current_iteration += 1;
            let z_fixed = to_fixed(&self.z);
            self.z_reference.push(z_fixed);
            self.z_tolerance.push(1e-6 * z_fixed.norm_sqr());
            if z_fixed.norm_sqr() >= 1e256 {
                break;
            };
        };
        self.current_iteration == self.maximum_iteration
    }
}