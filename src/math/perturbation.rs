use crate::util::{PixelData, FloatExp};

use rayon::prelude::*;
use crate::math::reference::Reference;

use crate::util::ComplexExtended;

use atomic_counter::{AtomicCounter, RelaxedCounter};

use std::sync::Arc;

pub struct Perturbation {}

impl Perturbation {
    pub fn iterate_normal(pixel_data: &mut [PixelData], reference: &Reference, pixels_complete: &Arc<RelaxedCounter>) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                let mut new_pixels_complete = 0;

                for pixel in pixel_data {
                    let mut scaled_iterations = 0;
                    let mut scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                    let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    let val1 = pixel.iteration - reference.start_iteration;
                    let val2 = reference.current_iteration - reference.start_iteration;
                    let reference_slice = &reference.reference_data[val1..=val2];

                    for additional_iterations in 0..=(reference.current_iteration - pixel.iteration) {
                        let reference_data = &reference_slice[additional_iterations];

                        // 2 multiplications and 2 adds
                        let z = reference_data.z_fixed + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                        if pixel.delta_current.exponent > -500 {
                            // 2 multiplications and one add
                            let z_norm = z.norm_sqr();

                            if z_norm < reference_data.z_tolerance {
                                pixel.iteration += additional_iterations;
                                pixel.glitched = true;
                                break;
                            }
    
                            if z_norm > 1e16 {
                                pixel.iteration += additional_iterations;
                                pixel.escaped = true;
                                pixel.delta_current.mantissa = pixel.delta_current.to_float();
                                pixel.delta_current.exponent = 0;
                                new_pixels_complete += 1;
                                break;
                            }
                        }

                        if reference_data.extended_precision_required {
                            // If the reference is small, use the slow extended method
                            pixel.delta_current *= reference_data.z_extended * 2.0 + pixel.delta_current;
                            pixel.delta_current += pixel.delta_reference;

                            // reset the scaled counter
                            pixel.delta_current.reduce();

                            scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                            scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                            scaled_iterations = 0;
                        } else {
                            // If the reference is not small, use the usual method

                            // 4 multiplications and 2 additions
                            pixel.delta_current.mantissa *= z + reference_data.z_fixed;
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
                    }

                    if !pixel.escaped && !pixel.glitched {
                        pixel.iteration = reference.current_iteration;
                        new_pixels_complete += 1;
                    }
                }

                pixels_complete.add(new_pixels_complete);
            });
    }

    pub fn iterate_normal_plus_derivative(pixel_data: &mut [PixelData], reference: &Reference, pixels_complete: &Arc<RelaxedCounter>) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                let mut new_pixels_complete = 0;

                for pixel in pixel_data {
                    let mut scaled_iterations = 0;
                    let mut scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                    let mut scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);
                    let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    let val1 = pixel.iteration - reference.start_iteration;
                    let val2 = reference.current_iteration - reference.start_iteration;
                    let reference_slice = &reference.reference_data[val1..=val2];

                    for additional_iterations in 0..=(reference.current_iteration - pixel.iteration) {
                        let reference_data = &reference_slice[additional_iterations];

                        // 2 multiplications and 2 adds
                        let z = reference_data.z_fixed + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                        if pixel.delta_current.exponent > -500 {
                            // 2 multiplications and one add
                            let z_norm = z.norm_sqr();

                            if z_norm < reference_data.z_tolerance {
                                pixel.iteration += additional_iterations;
                                pixel.glitched = true;
                                break;
                            }
    
                            if z_norm > 1e16 {
                                pixel.iteration += additional_iterations;
                                pixel.escaped = true;
                                pixel.delta_current.mantissa = pixel.delta_current.to_float();
                                pixel.delta_current.exponent = 0;
                                new_pixels_complete += 1;
                                break;
                            }
                        }

                        if reference_data.extended_precision_required {
                            // If the reference is small, use the slow extended method
                            pixel.derivative_current *= (reference_data.z_extended + pixel.delta_current) * 2.0;
                            pixel.derivative_current += ComplexExtended::new2(1.0, 0.0, 0);

                            pixel.delta_current *= reference_data.z_extended * 2.0 + pixel.delta_current;
                            pixel.delta_current += pixel.delta_reference;

                            // reset the scaled counter
                            pixel.delta_current.reduce();
                            pixel.derivative_current.reduce();

                            scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                            scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);
                            scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                            scaled_iterations = 0;
                        } else {
                            // If the reference is not small, use the usual method

                            // 4 multiplications and 2 additions
                            pixel.delta_current.mantissa *= z + reference_data.z_fixed;
                            // 2 additions
                            pixel.delta_current.mantissa += scaled_delta_reference;

                            pixel.derivative_current.mantissa *= 2.0 * z;
                            pixel.derivative_current.mantissa += scaled_scale_factor_2;

                            scaled_iterations += 1;

                            // check the counter, if it is > 250, do a normalisation
                            if scaled_iterations > 250 {
                                pixel.delta_current.reduce();
                                pixel.derivative_current.reduce();

                                scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                                scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);
                                scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                                scaled_iterations = 0;
                            }
                        }
                    }

                    pixel.derivative_current.reduce();

                    if !pixel.escaped && !pixel.glitched {
                        pixel.iteration = reference.current_iteration;
                        new_pixels_complete += 1;
                    }
                }

                pixels_complete.add(new_pixels_complete);
            });
    }
}



