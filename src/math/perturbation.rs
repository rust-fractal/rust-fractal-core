use crate::util::{PixelData, FloatExp};

use rayon::prelude::*;
use crate::math::reference::Reference;

pub struct Perturbation {}

impl Perturbation {
    pub fn iterate(pixel_data: &mut Vec<PixelData>, reference: &Reference) {
        pixel_data.par_chunks_mut(8)
            .for_each(|pixel_data| {
                for pixel in pixel_data {
                    let mut scaled_iterations = 0;
                    let mut scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                    let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    while pixel.iteration < reference.current_iteration {
                        // 2 multiplications and 2 adds
                        let z = reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                        if pixel.delta_current.exponent > -500 {
                            // 2 multiplications and one add
                            let z_norm = z.norm_sqr();

                            if z_norm < reference.reference_data[pixel.iteration - reference.start_iteration].z_tolerance {
                                pixel.glitched = true;
                                break;
                            }
    
                            if z_norm > 1e16 {
                                pixel.escaped = true;
                                pixel.delta_current.mantissa = pixel.delta_current.to_float();
                                pixel.delta_current.exponent = 0;
                                break;
                            }
                        }

                        if reference.reference_data[pixel.iteration - reference.start_iteration].extended_precision_required {
                            // If the reference is small, use the slow extended method
                            pixel.delta_current *= reference.reference_data[pixel.iteration - reference.start_iteration].z_extended * 2.0 + pixel.delta_current;
                            pixel.delta_current += pixel.delta_reference;

                            // reset the scaled counter
                            pixel.delta_current.reduce();

                            scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                            scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                            scaled_iterations = 0;
                        } else {
                            // If the reference is not small, use the usual method

                            // pixel.delta_current.mantissa = 2.0 * reference.data[pixel.iteration - reference.start_iteration].z_fixed * pixel.delta_current.mantissa + temp * pixel.delta_current.mantissa + scaled_delta_reference;

                            // 4 multiplications and 2 additions
                            pixel.delta_current.mantissa *= z + reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed;
                            // 2 additions
                            pixel.delta_current.mantissa += scaled_delta_reference;

                            scaled_iterations += 1;

                            // check the counter, if it is > 250, do a normalisation
                            if scaled_iterations > 250 {
                                pixel.delta_current.reduce();

                                scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                                scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                                scaled_iterations = 0;
                            }
                        }

                        pixel.iteration += 1;
                    }
                }
            });
    }
}



