use crate::util::{ComplexArbitrary, ComplexFixed, ComplexExtended, FloatExtended, to_fixed, to_extended};

use atomic_counter::{AtomicCounter, RelaxedCounter};

use std::sync::Arc;

#[derive(Clone)]
pub struct Reference {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub reference_data: Vec<ReferenceIteration>,
    pub extended_iterations: Vec<usize>,
    pub reference_data_extended: Vec<ComplexExtended>,
    // This is for every 100th iteration, when we do glitch correction the new references will be spawed from these values
    // Storing every iteration is memory intensive.
    pub zoom: FloatExtended,
    pub data_storage_interval: usize,
    pub high_precision_data: Vec<ComplexArbitrary>,
    pub glitch_tolerance: f64,
}

#[derive(Clone)]
pub struct ReferenceIteration {
    pub z: ComplexFixed<f64>,
    pub tolerance: f64,
}

impl Reference {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize, data_storage_interval: usize, glitch_tolerance: f64, zoom: FloatExtended) -> Reference {
        Reference {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            reference_data: Vec::new(),
            extended_iterations: Vec::new(),
            reference_data_extended: Vec::new(),
            zoom,
            data_storage_interval,
            high_precision_data: Vec::new(),
            glitch_tolerance
        }
    }

    pub fn step(&mut self) -> bool {
        self.z.square_mut();
        self.z += &self.c;
        self.current_iteration += 1;

        let z_fixed = to_fixed(&self.z);
        let tolerance = self.glitch_tolerance * z_fixed.norm_sqr();

        // This is if we need to use the extended precision for the reference
        if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
            // this is stored without the offset
            self.extended_iterations.push(self.current_iteration);
        }

        // We pack these together as they are always accessed together
        self.reference_data.push(
            ReferenceIteration {
                z: z_fixed,
                tolerance,
            }
        );

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        self.reference_data_extended.push(z_extended);

        // If the value is not small we do the escape check, otherwise it has not escaped
        // as we do the check for 65536 on the perturbation, we need this to be more than that squared
        z_fixed.norm_sqr() <= 1e256
    }


    pub fn run(&mut self, reference_counter: &Arc<RelaxedCounter>, reference_maximum_iteration_counter: &Arc<RelaxedCounter>, stop_flag: &Arc<RelaxedCounter>) {
        let z_fixed = to_fixed(&self.z);
        let tolerance = self.glitch_tolerance * z_fixed.norm_sqr();

        // This is if we need to use the extended precision for the reference
        if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
            self.extended_iterations.push(self.current_iteration);
        }

        // We pack these together as they are always accessed together
        // The first iteration is z_1=c = iteration 1 is index 0
        // access with iteration - start_iteration
        self.reference_data.push(
            ReferenceIteration {
                z: z_fixed,
                tolerance,
            }
        );

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        self.reference_data_extended.push(z_extended);

        while self.current_iteration < self.maximum_iteration {
            if self.start_iteration == 1 && (self.data_storage_interval == 1 || self.current_iteration % self.data_storage_interval == 1 ) {
                self.high_precision_data.push(self.z.clone());
            }

            if stop_flag.get() >= 1 {
                return
            };

            reference_counter.inc();

            if !self.step() {
                reference_maximum_iteration_counter.add(usize::max_value() - self.maximum_iteration + reference_counter.get() + 1);

                break;
            };
        }

        // println!("{:?}", self.extended_iterations);
    }
}

