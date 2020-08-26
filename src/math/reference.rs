use crate::util::{ComplexArbitrary, ComplexFixed, ComplexExtended, to_fixed, to_extended};

pub struct Reference {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub perturbation_data: Vec<ReferenceIteration>,
    pub approximation_data: Vec<ComplexExtended>,
    // This is for every 100th iteration, when we do glitch correction the new references will be spawed from these values
    // Storing every iteration is memory intensive.
    pub high_precision_data: Vec<ComplexArbitrary>
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
            perturbation_data: Vec::with_capacity(1000),
            approximation_data: Vec::with_capacity(1000),
            high_precision_data: Vec::with_capacity(1000)
        }
    }

    pub fn step(&mut self) -> bool {
        self.z.square_mut();
        self.z += &self.c;
        self.current_iteration += 1;

        let z_fixed = to_fixed(&self.z);
        let z_tolerance = 1e-6 * z_fixed.norm_sqr();

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        // This is if we need to use the extended precision for the reference
        if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
            // println!("found slow at: {}", self.current_iteration);
            self.perturbation_data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: Some(z_extended),
                    z_tolerance
                }
            )
        } else {
            self.perturbation_data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: None,
                    z_tolerance
                }
            )
        }

        self.approximation_data.push(z_extended);

        // If the value is not small we do the escape check, otherwise it has not escaped
        // as we do the check for 65536 on the perturbation, we need this to be more than that squared
        z_fixed.norm_sqr() <= 1e256
    }


    pub fn run(&mut self) {
        let z_fixed = to_fixed(&self.z);
        let z_tolerance = 1e-6 * z_fixed.norm_sqr();

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        // This is if we need to use the extended precision for the reference
        if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
            self.perturbation_data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: Some(z_extended),
                    z_tolerance
                }
            )
        } else {
            self.perturbation_data.push(
                ReferenceIteration {
                    z_fixed,
                    z_extended: None,
                    z_tolerance
                }
            )
        }

        self.approximation_data.push(z_extended);

        if self.start_iteration == 1 {
            // Add the first value for iterations 0-100
            self.high_precision_data.push(self.z.clone());

            while self.current_iteration < self.maximum_iteration {
                if !self.step() {
                    break;
                };
                if self.current_iteration % 100 == 1 {
                    self.high_precision_data.push(self.z.clone());
                }
            };
        } else {
            while self.current_iteration < self.maximum_iteration {
                if !self.step() {
                    break;
                }
            };
        }
    }
}

