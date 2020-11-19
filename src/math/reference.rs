use crate::util::{ComplexArbitrary, ComplexFixed, ComplexExtended, FloatExtended, to_fixed, to_extended};

#[derive(Clone)]
pub struct Reference {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub reference_data: Vec<ReferenceIteration>,
    // This is for every 100th iteration, when we do glitch correction the new references will be spawed from these values
    // Storing every iteration is memory intensive.
    pub zoom: FloatExtended,
    pub data_storage_interval: usize,
    pub high_precision_data: Vec<ComplexArbitrary>,
    pub glitch_tolerance: f64,
}

#[derive(Clone)]
pub struct ReferenceIteration {
    pub z_fixed: ComplexFixed<f64>,
    pub z_extended: ComplexExtended,
    pub z_tolerance: f64,
    pub extended_precision_required: bool
}

impl Reference {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize, data_storage_interval: usize, glitch_tolerance: f64, zoom: FloatExtended) -> Reference {
        Reference {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            reference_data: Vec::with_capacity(1000),
            zoom,
            data_storage_interval,
            high_precision_data: Vec::with_capacity(1000),
            glitch_tolerance
        }
    }

    pub fn step(&mut self) -> bool {
        self.z.square_mut();
        self.z += &self.c;
        self.current_iteration += 1;

        let z_fixed = to_fixed(&self.z);
        let z_tolerance = self.glitch_tolerance * z_fixed.norm_sqr();

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        // This is if we need to use the extended precision for the reference
        let extended_precision_required = z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300;

        self.reference_data.push(
            ReferenceIteration {
                z_fixed,
                z_extended,
                z_tolerance,
                extended_precision_required,
            }
        );

        // If the value is not small we do the escape check, otherwise it has not escaped
        // as we do the check for 65536 on the perturbation, we need this to be more than that squared
        z_fixed.norm_sqr() <= 1e256
    }


    pub fn run(&mut self) {
        let z_fixed = to_fixed(&self.z);
        let z_tolerance = self.glitch_tolerance * z_fixed.norm_sqr();

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        let extended_precision_required = z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300;

        // This is if we need to use the extended precision for the reference
        // We could use only the complex extended and check if the exponent is zero
        self.reference_data.push(
            ReferenceIteration {
                z_fixed,
                z_extended,
                z_tolerance,
                extended_precision_required,
            }
        );

        while self.current_iteration < self.maximum_iteration {
            if self.start_iteration == 1 {
                if self.data_storage_interval == 1 || self.current_iteration % self.data_storage_interval == 1 {
                    self.high_precision_data.push(self.z.clone());
                }
            }

            if !self.step() {
                break;
            };
        }
    }
}

