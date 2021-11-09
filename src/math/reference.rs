use rug::ops::PowAssign;

use crate::util::{ComplexArbitrary, ComplexFixed, ComplexExtended, FloatExtended, to_fixed, to_extended, FloatArbitrary};
use std::sync::{Arc, atomic::{AtomicUsize, AtomicBool, Ordering}};

#[derive(Clone)]
pub struct Reference {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub reference_data: Vec<ComplexFixed<f64>>,
    pub extended_iterations: Vec<usize>,
    pub reference_data_extended: Vec<ComplexExtended>,
    // This is for every 100th iteration, when we do glitch correction the new references will be spawed from these values
    // Storing every iteration is memory intensive.
    pub zoom: FloatExtended,
    pub data_storage_interval: usize,
    pub high_precision_data: Vec<ComplexArbitrary>,
    pub glitch_tolerance: f64,
}

// #[derive(Clone)]
// pub struct ReferenceIteration {
//     pub z: ComplexFixed<f64>
// }

impl Reference {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize, data_storage_interval: usize, glitch_tolerance: f64, zoom: FloatExtended) -> Reference {
        let zero = ComplexArbitrary::with_val(
            c.prec().0 as u32,
            ComplexArbitrary::parse("(0.0,0.0)").expect("provided location not valid"));

        Reference {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z: zero,
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

    pub fn run<const FRACTAL_TYPE: usize, const FRACTAL_POWER: usize>(&mut self, reference_counter: &Arc<AtomicUsize>, reference_maximum_iteration_counter: &Arc<AtomicUsize>, stop_flag: &Arc<AtomicBool>) {
        let z_fixed = to_fixed(&self.z);
        // let tolerance = self.glitch_tolerance * z_fixed.norm_sqr();

        // This is if we need to use the extended precision for the reference
        // if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
        //     self.extended_iterations.push(self.current_iteration);
        // }

        // We pack these together as they are always accessed together
        // The first iteration is z_1=c = iteration 1 is index 0
        // access with iteration - start_iteration
        self.reference_data.push(z_fixed);

        let mut z_extended = to_extended(&self.z);
        z_extended.reduce();

        self.reference_data_extended.push(z_extended);

        while self.current_iteration < self.maximum_iteration {
            if self.data_storage_interval == 1 || self.current_iteration % self.data_storage_interval == 1 {
                self.high_precision_data.push(self.z.clone());
            }

            if stop_flag.load(Ordering::SeqCst) {
                return;
            };

            reference_counter.fetch_add(1, Ordering::SeqCst);

            match FRACTAL_TYPE {
                1 => {
                    // Burning ship
                    self.z.mut_real().abs_mut();
                    self.z.mut_imag().abs_mut();

                    self.z.square_mut();
                    self.z += &self.c;
                }
                _ => {
                    match FRACTAL_POWER {
                        2 => {
                            // Power 2 mandelbrot
                            self.z.square_mut();
                            self.z += &self.c;
                        },
                        3 => {
                            // Power 3 mandelbrot
                            self.z *= self.z.clone().square();
                            self.z += &self.c;
                        },
                        _ => {
                            // Power N mandelbrot
                            self.z.pow_assign(FRACTAL_POWER as i64);
                            self.z += &self.c;
                        }
                    }
                }
            }

            self.current_iteration += 1;
    
            let z_fixed = to_fixed(&self.z);
            // let tolerance = self.glitch_tolerance * z_fixed.norm_sqr();
    
            // This is if we need to use the extended precision for the reference
            if z_fixed.re.abs() < 1e-300 && z_fixed.im.abs() < 1e-300 {
                // this is stored without the offset
                self.extended_iterations.push(self.current_iteration);
            }
    
            self.reference_data.push(z_fixed);
    
            let mut z_extended = to_extended(&self.z);
            z_extended.reduce();
    
            self.reference_data_extended.push(z_extended);
    
            // If the value is not small we do the escape check, otherwise it has not escaped
            // as we do the check for 65536 on the perturbation, we need this to be more than that squared
            if z_fixed.norm_sqr() >= 1e256 {
                break;
            }
        }

        reference_maximum_iteration_counter.store(self.current_iteration, Ordering::SeqCst);

        println!("{:?}", self.reference_data[0]);
        println!("{:?}", self.reference_data[1]);
        println!("{:?}", self.reference_data[2]);


        println!("{:?}", self.extended_iterations);
    }

    // This gets a reference that stores the high precision data every iteration
    pub fn get_central_glitch_resolving_reference(&self, iteration: usize) -> Reference {
        let iteration_reference = self.data_storage_interval * ((iteration - self.start_iteration) / self.data_storage_interval) + self.start_iteration;

        let reference_c = self.c.clone();
        let reference_z = self.high_precision_data[(iteration - self.start_iteration) / self.data_storage_interval].clone();

        Reference::new(reference_z, reference_c, iteration_reference, self.maximum_iteration, 1, self.glitch_tolerance, self.zoom)
    }

    // This is for use when getting new references others with full data
    pub fn get_glitch_resolving_reference(&self, iteration: usize, reference_delta: ComplexExtended, current_delta: ComplexExtended) -> Reference {
        assert!(self.data_storage_interval == 1);

        let precision = self.c.real().prec();

        let mut reference_c = self.c.clone();

        let temp = FloatArbitrary::with_val(precision, reference_delta.exponent).exp2();
        let temp2 = FloatArbitrary::with_val(precision, reference_delta.mantissa.re);
        let temp3 = FloatArbitrary::with_val(precision, reference_delta.mantissa.im);

        *reference_c.mut_real() += &temp2 * &temp;
        *reference_c.mut_imag() += &temp3 * &temp;

        let mut reference_z = self.high_precision_data[iteration - self.start_iteration].clone();

        let temp = FloatArbitrary::with_val(precision, current_delta.exponent).exp2();
        let temp2 = FloatArbitrary::with_val(precision, current_delta.mantissa.re);
        let temp3 = FloatArbitrary::with_val(precision, current_delta.mantissa.im);

        *reference_z.mut_real() += &temp2 * &temp;
        *reference_z.mut_imag() += &temp3 * &temp;

        Reference::new(reference_z, reference_c, iteration, self.maximum_iteration, 1, self.glitch_tolerance, self.zoom)
    }
}

