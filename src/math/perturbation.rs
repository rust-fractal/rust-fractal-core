use crate::util::PixelData;
use crate::util::FloatExp;

use rayon::prelude::*;
use crate::math::reference::Reference;

pub struct Perturbation {}

impl Perturbation {
    pub fn iterate(pixel_data: &mut Vec<PixelData>, reference: &Reference, reference_current_iteration: usize) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                for pixel in pixel_data {
                    let mut scaled_iterations = 0;

                    // maybe have a scale factor for the reference delta

                    let mut scaled_scale_factor_1 = 2.0f64.powi(pixel.delta_current.exponent);
                    let mut scaled_delta_reference = 2.0f64.powi(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    while pixel.iteration < reference_current_iteration {
                        let z_norm = (reference.data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.to_float()).norm_sqr();

                        if z_norm < reference.data[pixel.iteration - reference.start_iteration].z_tolerance {
                            pixel.glitched = true;
                            break;
                        }

                        if z_norm > 1e16 {
                            pixel.escaped = true;
                            break;
                        }

                        match reference.data[pixel.iteration - reference.start_iteration].z_extended {
                            // If the reference is small, use the slow extended method
                            Some(z_extended) => {
                                // do the slow 
                                pixel.delta_current = pixel.delta_current * z_extended * 2.0 + pixel.delta_current * pixel.delta_current + pixel.delta_reference;

                                // reset the scaled counter
                                pixel.delta_current.reduce();

                                scaled_scale_factor_1 = 2.0f64.powi(pixel.delta_current.exponent);
                                scaled_delta_reference = 2.0f64.powi(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                                scaled_iterations = 0;
                            },
                            // If the reference is not small, use the usual method
                            None => {
                                // possible to check if the exponent is below -500, if so dont bother with z^2
                                pixel.delta_current.mantissa = 2.0 * reference.data[pixel.iteration - reference.start_iteration].z_fixed * pixel.delta_current.mantissa + scaled_scale_factor_1 * pixel.delta_current.mantissa * pixel.delta_current.mantissa + scaled_delta_reference;

                                // pixel.delta_current.mantissa = if pixel.delta_current.exponent < -500 {
                                //     2.0 * reference.data[pixel.iteration - reference.start_iteration].z_fixed * pixel.delta_current.mantissa + scaled_delta_reference
                                // } else {
                                //     2.0 * reference.data[pixel.iteration - reference.start_iteration].z_fixed * pixel.delta_current.mantissa + scaled_scale_factor_1 * pixel.delta_current.mantissa * pixel.delta_current.mantissa + scaled_delta_reference
                                // };

                                scaled_iterations += 1;

                                // check the counter, if it is > 500, do a normalisation
                                if scaled_iterations > 500 {
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



