use crate::util::PixelData;

use rayon::prelude::*;
use crate::math::reference::Reference;

pub struct Perturbation {}

impl Perturbation {
    pub fn iterate(pixel_data: &mut Vec<PixelData>, reference: &Reference, reference_current_iteration: usize) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                for pixel in pixel_data {
                    let mut scaled_iterations = 0;
                    let mut scaled_scale_factor_1 = 2.0f64.powi(pixel.delta_current.exponent);
                    let mut scaled_delta_reference = 2.0f64.powi(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    while pixel.iteration < reference_current_iteration {
                        let delta_current_float = scaled_scale_factor_1 * pixel.delta_current.mantissa;

                        if pixel.delta_current.exponent > -250 {
                            let z_norm = (reference.data[pixel.iteration - reference.start_iteration].z_fixed + delta_current_float).norm_sqr();

                            if z_norm < reference.data[pixel.iteration - reference.start_iteration].z_tolerance {
                                pixel.glitched = true;
                                pixel.delta_current.reduce();
                                break;
                            }
    
                            if z_norm > 1e16 {
                                pixel.escaped = true;
                                pixel.delta_current.reduce();
                                break;
                            }
                        }

                        match reference.data[pixel.iteration - reference.start_iteration].z_extended {
                            // If the reference is small, use the slow extended method
                            Some(z_extended) => {
                                // do the slow 
                                pixel.delta_current *= z_extended * 2.0 + pixel.delta_current;
                                pixel.delta_current += pixel.delta_reference;

                                // reset the scaled counter
                                pixel.delta_current.reduce();

                                scaled_scale_factor_1 = 2.0f64.powi(pixel.delta_current.exponent);
                                scaled_delta_reference = 2.0f64.powi(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                                scaled_iterations = 0;
                            },
                            // If the reference is not small, use the usual method
                            None => {
                                // pixel.delta_current.mantissa = 2.0 * reference.data[pixel.iteration - reference.start_iteration].z_fixed * pixel.delta_current.mantissa + temp * pixel.delta_current.mantissa + scaled_delta_reference;

                                pixel.delta_current.mantissa *= 2.0 * reference.data[pixel.iteration - reference.start_iteration].z_fixed + delta_current_float;
                                pixel.delta_current.mantissa += scaled_delta_reference;

                                scaled_iterations += 1;

                                // check the counter, if it is > 250, do a normalisation
                                if scaled_iterations > 250 {
                                    pixel.delta_current.reduce();

                                    scaled_scale_factor_1 = 2.0f64.powi(pixel.delta_current.exponent);
                                    scaled_delta_reference = 2.0f64.powi(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                                    scaled_iterations = 0;
                                }
                            }
                        }

                        pixel.iteration += 1;
                    }
                }
            });
    }
}



