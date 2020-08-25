use crate::util::{ComplexArbitrary, ComplexFixed, ComplexExtended, to_fixed, to_extended};

pub struct Reference {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub data: Vec<ReferenceIteration>
}

pub struct ReferenceIteration {
    pub z_fixed: ComplexFixed<f64>,
    pub z_extended: Option<ComplexExtended>,
    pub z_tolerance: f64,
}

impl Reference {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize) -> Reference {
        Reference {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            data: Vec::with_capacity(1000)
        }
    }

    pub fn step(&mut self) -> bool {
        self.z.square_mut();
        self.z += &self.c;
        self.current_iteration += 1;

        let z_fixed = to_fixed(&self.z);
        let z_tolerance = 1e-6 * z_fixed.norm_sqr();

        // This is if we need to use the extended precision for the reference
        if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
        // if true {
            println!("found slow at: {}", self.current_iteration);
            let mut temp = to_extended(&self.z);
            temp.reduce();
            self.data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: Some(temp),
                    z_tolerance
                }
            )
        } else {
            self.data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: None,
                    z_tolerance
                }
            )
        }

        // If the value is not small we do the escape check, otherwise it has not escaped
        // as we do the check for 65536 on the perturbation, we need this to be more than that squared
        z_fixed.norm_sqr() <= 1e256
    }


    pub fn run(&mut self) -> bool {
        let z_fixed = to_fixed(&self.z);
        let z_tolerance = 1e-6 * z_fixed.norm_sqr();

        // This is if we need to use the extended precision for the reference
        if z_fixed.re.abs() < 1e300 && z_fixed.im.abs() < 1e300 {
            self.data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: Some(to_extended(&self.z)),
                    z_tolerance
                }
            )
        } else {
            self.data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: None,
                    z_tolerance
                }
            )
        }

        while self.current_iteration < self.maximum_iteration {
            if !self.step() {
                break;
            }
        };

        self.current_iteration == self.maximum_iteration
    }
}

