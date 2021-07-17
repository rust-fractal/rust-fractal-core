

use crate::util::{FloatExp, FloatExtended, PixelData, data_export::DataExport};

use rayon::prelude::*;
use crate::math::reference::Reference;
use crate::util::{ComplexExtended, ComplexFixed};

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use parking_lot::Mutex;

use crate::math::SeriesApproximation;

const ESCAPE_RADIUS: f64 = 1e16;

pub struct Perturbation {}

impl Perturbation {
    #[inline(always)]
    fn perturb_function<const DATA_TYPE: usize, const FRACTAL_TYPE: usize, const FRACTAL_POWER: usize>(
        delta_current_mantissa: &mut ComplexFixed<f64>, 
        derivative_current_mantissa: &mut ComplexFixed<f64>,
        z: ComplexFixed<f64>, 
        delta_reference: ComplexFixed<f64>, 
        scale_factor_1: f64,
        scale_factor_2: f64,
        pascal: &Vec<f64>) {

        match DATA_TYPE {
            1 | 3 => {
                match FRACTAL_TYPE {
                    _ => {
                        match FRACTAL_POWER {
                            2 => {
                                let temp = scale_factor_1 * *delta_current_mantissa;
                                let temp2 = 2.0 * z;
        
                                *derivative_current_mantissa *= temp2 + 2.0 * temp;
                                *derivative_current_mantissa += scale_factor_2;
        
                                *delta_current_mantissa *= temp2 + temp;
                                *delta_current_mantissa += delta_reference;
                            },
                            3 => {
                                let temp = scale_factor_1 * *delta_current_mantissa;
                                let temp2 = 3.0 * z;
                                let temp3 = z + temp;
        
                                *derivative_current_mantissa *= temp3 * temp3 * 3.0;
                                *derivative_current_mantissa += scale_factor_2;
        
                                *delta_current_mantissa *= temp2 * temp3 + temp * temp;
                                *delta_current_mantissa += delta_reference;
                            },
                            // This should be a generic implementation for mandelbrot powers > 3
                            _ => {
                                let temp = scale_factor_1 * *delta_current_mantissa;

                                *derivative_current_mantissa *= FRACTAL_POWER as f64 * (z + temp).powi(FRACTAL_POWER as i32 - 1);
                                *derivative_current_mantissa += scale_factor_2;

                                // for 3rd power 3Z^2 + 3Z * z + z * z
                                // for 4th power 4Z^3 + 6Z^2 * z + 4Z * z^2 + z^3
                                let mut sum = pascal[1] * z + temp;
                                let mut z_p = z;

                                for i in 2..FRACTAL_POWER { 
                                    sum *= temp; 
                                    z_p *= z;
                                    sum += pascal[i] * z_p;
                                }

                                *delta_current_mantissa *= sum;
                                *delta_current_mantissa += delta_reference;
                            }
                        }
                    }
                }
            },
            _ => {
                match FRACTAL_TYPE {
                    _ => {
                        match FRACTAL_POWER {
                            2 => {
                                *delta_current_mantissa *= 2.0 * z + scale_factor_1 * *delta_current_mantissa;
                                *delta_current_mantissa += delta_reference;
                            },
                            3 => {
                                let temp = scale_factor_1 * *delta_current_mantissa;
                                let temp2 = 3.0 * z;

                                *delta_current_mantissa *= temp2 * z + temp * (temp2 + temp);
                                *delta_current_mantissa += delta_reference;
                            }
                            // This should be a generic implementation for mandelbrot powers > 3
                            _ => {
                                let temp = scale_factor_1 * *delta_current_mantissa;

                                // for 3rd power 3Z^2 + 3Z * z + z * z
                                // for 4th power 4Z^3 + 6Z^2 * z + 4Z * z^2 + z^3
                                let mut sum = pascal[1] * z + temp;
                                let mut z_p = z;

                                for i in 2..FRACTAL_POWER { 
                                    sum *= temp; 
                                    z_p *= z;
                                    sum += pascal[i] * z_p;
                                }

                                *delta_current_mantissa *= sum;
                                *delta_current_mantissa += delta_reference;
                            }
                        }
                        
                    }
                }
            }
        }
    }

    pub fn iterate<const DATA_TYPE: usize, const FRACTAL_TYPE: usize, const FRACTAL_POWER: usize>(
        pixel_data: &mut [PixelData], 
        reference: &Reference, 
        pixels_complete: &Arc<AtomicUsize>, 
        stop_flag: &Arc<AtomicBool>, 
        data_export: Arc<Mutex<DataExport>>, 
        delta_pixel: FloatExtended, 
        scale: usize, 
        chunk_size: usize, 
        series_approximation: &SeriesApproximation, 
        initial: bool,
        pascal: &Vec<f64>) {
        pixel_data.par_chunks_mut(chunk_size)
        .for_each(|pixel_data| {
            // Record the number of new pixels that have been completed
            let mut new_pixels_complete = 0;
            let mut pixel_index = 0;

            // Go through each pixel in the packet
            for pixel in pixel_data.iter_mut() {
                // Check if the stop flag has been hit
                if stop_flag.load(Ordering::SeqCst) {
                    break;
                };

                if initial {
                    pixel.delta_current = series_approximation.evaluate(pixel.delta_reference, pixel.iteration);

                    if DATA_TYPE == 1 || DATA_TYPE == 3 {
                        pixel.derivative_current = series_approximation.evaluate_derivative(pixel.delta_reference, pixel.iteration);
                    }
                }

                pixel_index += 1;

                // Variable to record the number of additional iterations
                let mut additional_iterations = 0;

                // Scaled factors and reference values for the scaled double implementation
                let mut scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                let mut scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);

                let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                // Get the reference slice that is worked on
                let val1 = pixel.iteration - reference.start_iteration;
                let val2 = reference.current_iteration - reference.start_iteration;
                let val3 = reference.current_iteration - pixel.iteration;
                let reference_slice = &reference.reference_data[val1..=val2];

                // Get the number of iterations to the first extended iteration
                let (mut extended_index, &first_extended_iteration) = reference.extended_iterations
                    .iter()
                    .enumerate()
                    .find(|&(_, &value)| value >= pixel.iteration)
                    .unwrap_or((0, &0xFFFFFFFF));

                // Number of iterations to the next extended iteration
                let mut next_extended_iteration = first_extended_iteration - pixel.iteration;

                // Start running the core iteration loop
                'outer: loop {
                    // Number of iterations remaining
                    let iterations_remaining = val3 - additional_iterations;
                    let mut next_iteration_batch = iterations_remaining.min(10);

                    // Check if we need to to a new extended iteration
                    let need_extended_iteration = if next_extended_iteration < next_iteration_batch {
                        next_iteration_batch = next_extended_iteration;
                        true
                    } else {
                        false
                    };

                    let need_escape_check = pixel.delta_current.exponent > -500;
                    let reference_batch = &reference_slice[additional_iterations..(additional_iterations + next_iteration_batch)];

                    // If we need to check for escaping in these iterations
                    if need_escape_check {
                        for (i, reference_data) in reference_batch.iter().enumerate() {
                            // This is Z + z
                            let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;
                            let z_norm = z.norm_sqr();

                            if z_norm < reference_data.tolerance {
                                pixel.iteration += additional_iterations + i;
                                pixel.glitched = true;

                                break 'outer;
                            }

                            if z_norm > ESCAPE_RADIUS {
                                pixel.iteration += additional_iterations + i;
                                pixel.z_norm = z_norm;
                                pixel.delta_current.mantissa = pixel.delta_current.to_float();
                                pixel.delta_current.exponent = 0;

                                new_pixels_complete += 1;
                                break 'outer;
                            }

                            if DATA_TYPE == 2 || DATA_TYPE == 3 {
                                pixel.stripe_iteration += 1;
                                pixel.stripe_iteration %= 4;

                                pixel.stripe_storage[pixel.stripe_iteration] = z;
                            }

                            Perturbation::perturb_function::<DATA_TYPE, FRACTAL_TYPE, FRACTAL_POWER>(
                                &mut pixel.delta_current.mantissa,
                                &mut pixel.derivative_current.mantissa,
                                reference_data.z,
                                scaled_delta_reference,
                                scaled_scale_factor_1,
                                scaled_scale_factor_2,
                                pascal
                            )
                        }
                    } else {
                        for reference_data in reference_batch.iter() {
                            Perturbation::perturb_function::<DATA_TYPE, FRACTAL_TYPE, FRACTAL_POWER>(
                                &mut pixel.delta_current.mantissa,
                                &mut pixel.derivative_current.mantissa,
                                reference_data.z,
                                scaled_delta_reference,
                                scaled_scale_factor_1,
                                scaled_scale_factor_2,
                                pascal
                            )
                        }
                    }

                    // If we have hit the iteration limit
                    if iterations_remaining == next_iteration_batch {
                        if (pixel.iteration + additional_iterations + next_iteration_batch) < reference.maximum_iteration {
                            pixel.glitched = true;
                            pixel.iteration = reference.current_iteration;
                        } else {
                            pixel.iteration = reference.maximum_iteration;
                        }

                        new_pixels_complete += 1;
                        break;
                    }

                    additional_iterations += next_iteration_batch;
                    next_extended_iteration -= next_iteration_batch;

                    if need_extended_iteration {
                        let reference_data = &reference_slice[additional_iterations];
                        let z = reference_data.z + scaled_scale_factor_1 * pixel.delta_current.mantissa;

                        if need_escape_check {
                            let z_norm = z.norm_sqr();

                            if z_norm < reference_data.tolerance {
                                pixel.iteration += additional_iterations;
                                pixel.glitched = true;

                                break;
                            }

                            if z_norm > ESCAPE_RADIUS {
                                pixel.iteration += additional_iterations;
                                pixel.z_norm = z_norm;
                                pixel.delta_current.mantissa = pixel.delta_current.to_float();
                                pixel.delta_current.exponent = 0;
                                
                                new_pixels_complete += 1;
                                break;
                            }

                            if DATA_TYPE == 2 || DATA_TYPE == 3 {
                                pixel.stripe_iteration += 1;
                                pixel.stripe_iteration %= 4;

                                pixel.stripe_storage[pixel.stripe_iteration] = z;
                            }
                        }

                        if DATA_TYPE == 2 || DATA_TYPE == 3 {
                            pixel.derivative_current *= (reference.reference_data_extended[val1 + additional_iterations] + pixel.delta_current) * 2.0;
                            pixel.derivative_current += ComplexExtended::new2(1.0, 0.0, 0);
                        }

                        pixel.delta_current *= reference.reference_data_extended[val1 + additional_iterations] * 2.0 + pixel.delta_current;
                        pixel.delta_current += pixel.delta_reference;
                        
                        additional_iterations += 1;

                        extended_index += 1;
                        next_extended_iteration = if extended_index < reference.extended_iterations.len() {
                            reference.extended_iterations[extended_index] - additional_iterations - pixel.iteration 
                        } else {
                            0xFFFFFFFF
                        };
                    }

                    pixel.delta_current.reduce();

                    if DATA_TYPE == 1 || DATA_TYPE == 3 {
                        pixel.derivative_current.reduce();
                        scaled_scale_factor_2 = 1.0f64.ldexp(-pixel.derivative_current.exponent);
                    }

                    scaled_scale_factor_1 = 1.0f64.ldexp(pixel.delta_current.exponent);
                    scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;
                }
            }

            data_export.lock().export_pixels(&pixel_data[0..pixel_index], reference, delta_pixel, scale);
            pixels_complete.fetch_add(new_pixels_complete, Ordering::Relaxed);
        });
    }
}