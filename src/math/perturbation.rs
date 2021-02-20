use crate::util::{PixelData, FloatExp};

use rayon::prelude::*;
use crate::math::reference::Reference;

use crate::util::ComplexExtended;

use atomic_counter::{AtomicCounter, RelaxedCounter};

use std::{sync::Arc};

pub struct Perturbation {}

impl Perturbation {
    pub fn iterate_normal(pixel_data: &mut [PixelData], reference: &Reference, pixels_complete: &Arc<RelaxedCounter>, stop_flag: &Arc<RelaxedCounter>) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                if stop_flag.get() >= 1 {
                    return;
                }

                let mut new_pixels_complete = 0;

                for pixel in pixel_data {
                    // let mut scaled_iterations = 0;
                    let mut scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                    let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    let val1 = pixel.iteration - reference.start_iteration;
                    let val2 = reference.current_iteration - reference.start_iteration;
                    let reference_slice = &reference.reference_data[val1..=val2];

                    // get the first iteration that needs to be done in extended precision
                    let (first_extended_index, &first_extended_iteration) = reference.extended_iterations
                        .iter()
                        .enumerate()
                        .find(|&(_, &value)| value >= pixel.iteration)
                        .unwrap_or((0, &0xFFFFFFFF));

                    let mut first_extended_index = first_extended_index;
                    let mut additional_extended_iteration = first_extended_iteration - pixel.iteration;
                    let mut additional_iterations = 0;

                    // loop in validated chunks of iterations
                    loop {
                        let mut next_iteration_batch = 250;
                        let mut need_extended_iteration = false;

                        // check if we need to stop because of max iterations (within 250 of max)
                        if reference.current_iteration - pixel.iteration + additional_iterations < 250 {
                            next_iteration_batch = reference.current_iteration - pixel.iteration + additional_iterations
                        };

                        if additional_extended_iteration - additional_iterations <= 250 {
                            let temp = additional_extended_iteration - additional_iterations - 1;

                            if temp < next_iteration_batch {
                                next_iteration_batch = temp;
                                need_extended_iteration = true;
                            }

                            if first_extended_index < reference.extended_iterations.len() - 1 {
                                first_extended_index += 1;
                                additional_extended_iteration = reference.extended_iterations[first_extended_index] - pixel.iteration;
                            }
                        };

                        let need_escape_check = pixel.delta_current.exponent > -500;

                        if need_escape_check {
                            for _ in 0..next_iteration_batch {
                                let reference_data = &reference_slice[additional_iterations];

                                // 2 multiplications and 2 adds
                                let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;
        
                                // 2 multiplications and one add
                                let z_norm = z.norm_sqr();
    
                                if z_norm < reference_data.tolerance {
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

                                // 4 multiplications and 2 additions
                                pixel.delta_current.mantissa *= z + reference_data.z;
                                // 2 additions
                                pixel.delta_current.mantissa += scaled_delta_reference;

                                additional_iterations += 1;
                            }

                            if pixel.glitched || pixel.escaped {
                                break;
                            }
                        } else {
                            // assuming the point will not escape in the these iterations
                            for _ in 0..next_iteration_batch {
                                let reference_data = &reference_slice[additional_iterations];

                                // 2 multiplications and 2 adds
                                let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                                // 4 multiplications and 2 additions
                                pixel.delta_current.mantissa *= z + reference_data.z;
                                // 2 additions
                                pixel.delta_current.mantissa += scaled_delta_reference;
                            }

                            additional_iterations += next_iteration_batch
                        }

                        if pixel.iteration + additional_iterations > reference.current_iteration {
                            pixel.iteration = reference.current_iteration;
                            new_pixels_complete += 1;
                            break;
                        }

                        if need_extended_iteration {
                            let reference_data = &reference_slice[additional_iterations];

                            // 2 multiplications and 2 adds
                            let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                            // TODO this could be optimised out as well
                            if need_escape_check {
                                // 2 multiplications and one add
                                let z_norm = z.norm_sqr();
    
                                if z_norm < reference_data.tolerance {
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

                            pixel.delta_current *= reference.reference_data_extended[val1 + additional_iterations] * 2.0 + pixel.delta_current;
                            pixel.delta_current += pixel.delta_reference;

                            additional_iterations += 1;
                        }

                        pixel.delta_current.reduce();

                        scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                        scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;
                    }
                }

                pixels_complete.add(new_pixels_complete);
            });
    }

    pub fn iterate_normal_plus_derivative(pixel_data: &mut [PixelData], reference: &Reference, pixels_complete: &Arc<RelaxedCounter>, stop_flag: &Arc<RelaxedCounter>) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                if stop_flag.get() >= 1 {
                    return;
                }

                let mut new_pixels_complete = 0;

                for pixel in pixel_data {
                    // let mut scaled_iterations = 0;
                    let mut scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                    let mut scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);
                    let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                    let val1 = pixel.iteration - reference.start_iteration;
                    let val2 = reference.current_iteration - reference.start_iteration;
                    let reference_slice = &reference.reference_data[val1..=val2];

                    // get the first iteration that needs to be done in extended precision
                    let (first_extended_index, &first_extended_iteration) = reference.extended_iterations
                        .iter()
                        .enumerate()
                        .find(|&(_, &value)| value >= pixel.iteration)
                        .unwrap_or((0, &0xFFFFFFFF));

                    let mut first_extended_index = first_extended_index;
                    let mut additional_extended_iteration = first_extended_iteration - pixel.iteration;
                    let mut additional_iterations = 0;

                    // loop in validated chunks of iterations
                    loop {
                        let mut next_iteration_batch = 250;
                        let mut need_extended_iteration = false;

                        // check if we need to stop because of max iterations (within 250 of max)
                        if reference.current_iteration - pixel.iteration + additional_iterations < 250 {
                            next_iteration_batch = reference.current_iteration - pixel.iteration + additional_iterations
                        };

                        if additional_extended_iteration - additional_iterations <= 250 {
                            let temp = additional_extended_iteration - additional_iterations - 1;

                            if temp < next_iteration_batch {
                                next_iteration_batch = temp;
                                need_extended_iteration = true;
                            }

                            if first_extended_index < reference.extended_iterations.len() - 1 {
                                first_extended_index += 1;
                                additional_extended_iteration = reference.extended_iterations[first_extended_index] - pixel.iteration;
                            }
                        };

                        let need_escape_check = pixel.delta_current.exponent > -500;

                        if need_escape_check {
                            for _ in 0..next_iteration_batch {
                                let reference_data = &reference_slice[additional_iterations];

                                // 2 multiplications and 2 adds
                                let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;
        
                                // 2 multiplications and one add
                                let z_norm = z.norm_sqr();
    
                                if z_norm < reference_data.tolerance {
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

                                // 4 multiplications and 2 additions
                                pixel.delta_current.mantissa *= z + reference_data.z;
                                // 2 additions
                                pixel.delta_current.mantissa += scaled_delta_reference;

                                pixel.derivative_current.mantissa *= 2.0 * z;
                                pixel.derivative_current.mantissa += scaled_scale_factor_2;

                                additional_iterations += 1;
                            }

                            if pixel.glitched || pixel.escaped {
                                pixel.derivative_current.reduce();
                                break;
                            }
                        } else {
                            // assuming the point will not escape in the these iterations
                            for _ in 0..next_iteration_batch {
                                let reference_data = &reference_slice[additional_iterations];

                                // 2 multiplications and 2 adds
                                let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                                // 4 multiplications and 2 additions
                                pixel.delta_current.mantissa *= z + reference_data.z;
                                // 2 additions
                                pixel.delta_current.mantissa += scaled_delta_reference;

                                pixel.derivative_current.mantissa *= 2.0 * z;
                                pixel.derivative_current.mantissa += scaled_scale_factor_2;
                            }

                            additional_iterations += next_iteration_batch
                        }

                        // check if the pixel escapes
                        if pixel.iteration + additional_iterations > reference.current_iteration {
                            pixel.iteration = reference.current_iteration;
                            new_pixels_complete += 1;
                            break;
                        }

                        if need_extended_iteration {
                            let reference_data = &reference_slice[additional_iterations];

                            // 2 multiplications and 2 adds
                            let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                            // TODO this could be optimised out as well
                            if need_escape_check {
                                // 2 multiplications and one add
                                let z_norm = z.norm_sqr();
    
                                if z_norm < reference_data.tolerance {
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

                            pixel.delta_current *= reference.reference_data_extended[val1 + additional_iterations] * 2.0 + pixel.delta_current;
                            pixel.delta_current += pixel.delta_reference;

                            pixel.derivative_current *= (reference.reference_data_extended[val1 + additional_iterations] + pixel.delta_current) * 2.0;
                            pixel.derivative_current += ComplexExtended::new2(1.0, 0.0, 0);

                            additional_iterations += 1;
                        }

                        pixel.delta_current.reduce();
                        pixel.derivative_current.reduce();

                        scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                        scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);
                        scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;
                    }
                }

                pixels_complete.add(new_pixels_complete);
            });
    }
}



