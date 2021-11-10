

use crate::util::{FloatExp, FloatExtended, PixelData, data_export::DataExport};

use rayon::prelude::*;
use crate::math::reference::Reference;
use crate::util::{ComplexExtended, ComplexFixed, diff_abs};

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use parking_lot::Mutex;

use crate::math::SeriesApproximation;

const ESCAPE_RADIUS: f64 = 1e16;

pub struct Perturbation {}

impl Perturbation {
    #[inline(always)]
    fn perturb_function<const DATA_TYPE: usize, const FRACTAL_TYPE: usize>(
        delta_current_mantissa: &mut ComplexFixed<f64>, 
        jacobian: &mut [ComplexExtended; 2],
        z: ComplexFixed<f64>, 
        delta_reference: ComplexFixed<f64>, 
        scale_factor_1: f64,
        scale_factor_2: f64,
        pascal: &Vec<f64>,
        fractal_power: usize) {

        match DATA_TYPE {
            1 | 3 => {
                match FRACTAL_TYPE {
                    1 => {
                        let jacobian_a_copy = jacobian[0].mantissa.clone();
                        let jacobian_b_copy = jacobian[1].mantissa.clone();

                        let loc = z + scale_factor_1 * *delta_current_mantissa;

                        let sign = 2.0 * (loc.re * loc.im).signum();

                        // TODO needs to be scaled
                        jacobian[0].mantissa.re = 2.0 * (loc.re * jacobian_a_copy.re - loc.im * jacobian_b_copy.re) + scale_factor_2;
                        jacobian[0].mantissa.im = 2.0 * (loc.re * jacobian_a_copy.im - loc.im * jacobian_b_copy.im);
                        jacobian[1].mantissa.re = sign * (loc.re * jacobian_b_copy.re + loc.im * jacobian_a_copy.re);
                        jacobian[1].mantissa.im = sign * (loc.re * jacobian_b_copy.im + loc.im * jacobian_a_copy.im) + scale_factor_2;

                        let temp_re = delta_current_mantissa.re;

                        delta_current_mantissa.re = (2.0 * z.re + temp_re * scale_factor_1) * temp_re - (2.0 * z.im + delta_current_mantissa.im * scale_factor_1) * delta_current_mantissa.im;
                        delta_current_mantissa.im = 2.0 * diff_abs(z.re * z.im / scale_factor_1, z.re * delta_current_mantissa.im + temp_re * (z.im + delta_current_mantissa.im * scale_factor_1));

                        *delta_current_mantissa += delta_reference;
                    }
                    _ => {
                        match fractal_power {
                            2 => {
                                let temp = scale_factor_1 * *delta_current_mantissa;
                                let temp2 = 2.0 * z;
        
                                jacobian[0].mantissa *= temp2 + 2.0 * temp;
                                jacobian[0].mantissa += scale_factor_2;
        
                                *delta_current_mantissa *= temp2 + temp;
                                *delta_current_mantissa += delta_reference;
                            },
                            3 => {
                                let temp = scale_factor_1 * *delta_current_mantissa;
                                let temp2 = 3.0 * z;
                                let temp3 = z + temp;
        
                                jacobian[0].mantissa *= temp3 * temp3 * 3.0;
                                jacobian[0].mantissa += scale_factor_2;
        
                                *delta_current_mantissa *= temp2 * temp3 + temp * temp;
                                *delta_current_mantissa += delta_reference;
                            },
                            // This should be a generic implementation for mandelbrot powers > 3
                            _ => {
                                let temp = scale_factor_1 * *delta_current_mantissa;

                                jacobian[0].mantissa *= fractal_power as f64 * (z + temp).powi(fractal_power as i32 - 1);
                                jacobian[0].mantissa += scale_factor_2;

                                // for 3rd power 3Z^2 + 3Z * z + z * z
                                // for 4th power 4Z^3 + 6Z^2 * z + 4Z * z^2 + z^3
                                let mut sum = pascal[1] * z + temp;
                                let mut z_p = z;

                                for i in 2..fractal_power { 
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
                    1 => {
                        let temp_re = delta_current_mantissa.re;

                        delta_current_mantissa.re = (2.0 * z.re + temp_re * scale_factor_1) * temp_re - (2.0 * z.im + delta_current_mantissa.im * scale_factor_1) * delta_current_mantissa.im;
                        delta_current_mantissa.im = 2.0 * diff_abs(z.re * z.im / scale_factor_1, z.re * delta_current_mantissa.im + temp_re * (z.im + delta_current_mantissa.im * scale_factor_1));

                        *delta_current_mantissa += delta_reference;
                    }
                    _ => {
                        match fractal_power {
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

                                for i in 2..fractal_power { 
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

    #[inline(always)]
    fn perturb_function_extended<const DATA_TYPE: usize, const FRACTAL_TYPE: usize>(
        delta_current: &mut ComplexExtended, 
        jacobian: &mut [ComplexExtended; 2],
        z: ComplexExtended, 
        delta_reference: ComplexExtended, 
        pascal: &Vec<f64>,
        fractal_power: usize) {

        match DATA_TYPE {
            1 | 3 => {
                match FRACTAL_TYPE {
                    1 => {
                        unimplemented!()
                    }
                    _ => {
                        match fractal_power {
                            2 => {
                                jacobian[0] *= (z + *delta_current) * 2.0;
                                jacobian[0] += ComplexExtended::new2(1.0, 0.0, 0);

                                *delta_current *= z * 2.0 + *delta_current;
                                *delta_current += delta_reference;
                            },
                            // This should be a generic implementation for mandelbrot powers > 3
                            _ => {
                                jacobian[0] *= (z + *delta_current).powi(fractal_power as i32 - 1) * fractal_power as f64;
                                jacobian[0] += ComplexExtended::new2(1.0, 0.0, 0);

                                let mut sum = z * pascal[1] + *delta_current;
                                let mut z_p = z;

                                for i in 2..fractal_power { 
                                    sum *= *delta_current; 
                                    z_p *= z;
                                    sum += z_p * pascal[i];
                                }

                                *delta_current *= sum;
                                *delta_current += delta_reference;
                            }
                        }
                    }
                }
            },
            _ => {
                match FRACTAL_TYPE {
                    1 => {
                        unimplemented!()
                    }
                    _ => {
                        match fractal_power {
                            2 => {
                                *delta_current *= z * 2.0 + *delta_current;
                                *delta_current += delta_reference;
                            },
                            _ => {
                                let mut sum = z * pascal[1] + *delta_current;
                                let mut z_p = z;

                                for i in 2..fractal_power { 
                                    sum *= *delta_current; 
                                    z_p *= z;
                                    sum += z_p * pascal[i];
                                }

                                *delta_current *= sum;
                                *delta_current += delta_reference;
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
        _initial: bool,
        pascal: &Vec<f64>) {

        let iterations_before_check = 400 / FRACTAL_POWER;

        pixel_data.par_chunks_mut(chunk_size)
        .for_each(|pixel_data| {
            // Record the number of new pixels that have been completed
            let mut new_pixels_complete = 0;
            let mut pixel_index = 0;

            // Iterate through each pixel in the packet of pixels
            for pixel in pixel_data.iter_mut() {
                // Check if the stop flag has been hit and STOP
                if stop_flag.load(Ordering::SeqCst) {
                    break;
                };

                // Evaluate series approximation
                pixel.delta_current = series_approximation.evaluate(pixel.delta_reference, pixel.iteration);

                if DATA_TYPE == 1 || DATA_TYPE == 3 {
                    pixel.jacobian_current[0] = series_approximation.evaluate_derivative(pixel.delta_reference, pixel.iteration);
                    pixel.jacobian_current[1].scale_to_exponent(pixel.jacobian_current[0].exponent);
                }

                pixel_index += 1;

                // Get index into reference array
                let mut reference_index = pixel.iteration;

                // Scaled factors and reference values for the scaled double implementation
                let mut scale_factor_delta = 1.0f64.ldexp(pixel.delta_current.exponent);
                let mut scale_factor_derivative = 1.0f64.ldexp(-pixel.jacobian_current[0].exponent);
                let mut scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;

                // Get the number of iterations to the first extended iteration
                let (mut extended_index, &next_extended_iteration) = reference.extended_iterations
                    .iter()
                    .enumerate()
                    .find(|&(_, &value)| value >= pixel.iteration)
                    .unwrap_or((0, &0xFFFFFFFF));

                let mut next_extended_iteration = next_extended_iteration;

                let first_extended_iteration = *reference.extended_iterations.get(0).unwrap_or(&0xFFFFFFFF);

                // CORE ITERATION LOOP
                'outer: loop {
                    let iterations_remaining = reference.maximum_iteration - pixel.iteration;
                    let next_iteration_batch = iterations_remaining
                        .min(iterations_before_check)
                        .min(next_extended_iteration - reference_index);

                    // If we should be doing escape checks
                    if pixel.delta_current.exponent > -500 {
                        // for loop to avoid bounds checks


                        let mut i = 0;

                        while i < next_iteration_batch && reference_index < next_extended_iteration {
                            let mut reference_z = reference.reference_data[reference_index];

                            let z = reference_z + scale_factor_delta * pixel.delta_current.mantissa;
                            let z_norm = z.norm_sqr();

                            // Check for escape
                            if z_norm > ESCAPE_RADIUS {
                                pixel.iteration += i;
                                pixel.reference_iteration = reference_index;
                                pixel.z_norm = z_norm;
                                pixel.delta_current.mantissa = pixel.delta_current.to_float();
                                pixel.delta_current.exponent = 0;

                                new_pixels_complete += 1;
                                break 'outer;
                            }

                            // Add iterations to stripe storage
                            if DATA_TYPE == 2 || DATA_TYPE == 3 {
                                pixel.stripe_iteration += 1;
                                pixel.stripe_iteration %= 4;

                                pixel.stripe_storage[pixel.stripe_iteration] = z;
                            }

                            // Check - could be optimised
                            if z_norm < (scale_factor_delta * pixel.delta_current.mantissa).norm_sqr() || reference_index == reference.current_iteration {
                                pixel.delta_current.mantissa = z / scale_factor_delta;

                                reference_index = 0;
                                reference_z = reference.reference_data[0];

                                extended_index = 0;
                                next_extended_iteration = first_extended_iteration;
                            }

                            Perturbation::perturb_function::<DATA_TYPE, FRACTAL_TYPE>(
                                &mut pixel.delta_current.mantissa,
                                &mut pixel.jacobian_current,
                                reference_z,
                                scaled_delta_reference,
                                scale_factor_delta,
                                scale_factor_derivative,
                                pascal,
                                FRACTAL_POWER
                            );

                            i += 1;
                            reference_index += 1;
                        };

                        pixel.iteration += i;
                    } else {
                        for reference_z in reference.reference_data[reference_index..(reference_index + next_iteration_batch)].iter() {
                            Perturbation::perturb_function::<DATA_TYPE, FRACTAL_TYPE>(
                                &mut pixel.delta_current.mantissa,
                                &mut pixel.jacobian_current,
                                *reference_z,
                                scaled_delta_reference,
                                scale_factor_delta,
                                scale_factor_derivative,
                                pascal,
                                FRACTAL_POWER
                            )
                        }

                        reference_index += next_iteration_batch;
                        pixel.iteration += next_iteration_batch;
                    }

                    if pixel.iteration >= reference.maximum_iteration {
                        // println!("hit3 {} {}", pixel.iteration, reference.maximum_iteration);
                        pixel.iteration = reference.maximum_iteration;
                        pixel.reference_iteration = reference_index;
                        new_pixels_complete += 1;
                        break;
                    }

                    if reference_index == next_extended_iteration {
                        let reference_z = reference.reference_data[reference_index];

                        let z = reference_z + scale_factor_delta * pixel.delta_current.mantissa;
                        let z_norm = z.norm_sqr();

                        // Check for escape
                        if z_norm > ESCAPE_RADIUS {
                            pixel.iteration += 1;
                            pixel.z_norm = z_norm;
                            pixel.delta_current.mantissa = pixel.delta_current.to_float();
                            pixel.delta_current.exponent = 0;

                            new_pixels_complete += 1;
                            break 'outer;
                        }

                        // Add iterations to stripe storage
                        if DATA_TYPE == 2 || DATA_TYPE == 3 {
                            pixel.stripe_iteration += 1;
                            pixel.stripe_iteration %= 4;

                            pixel.stripe_storage[pixel.stripe_iteration] = z;
                        }

                        // Check - could be optimised
                        if z_norm < (scale_factor_delta * pixel.delta_current.mantissa).norm_sqr() || reference_index == reference.current_iteration {
                            pixel.delta_current.mantissa = z / scale_factor_delta;

                            reference_index = 0;

                            extended_index = 0;
                            next_extended_iteration = first_extended_iteration;

                            Perturbation::perturb_function::<DATA_TYPE, FRACTAL_TYPE>(
                                &mut pixel.delta_current.mantissa,
                                &mut pixel.jacobian_current,
                                reference.reference_data[0],
                                scaled_delta_reference,
                                scale_factor_delta,
                                scale_factor_derivative,
                                pascal,
                                FRACTAL_POWER
                            );
                        } else {
                            Perturbation::perturb_function_extended::<DATA_TYPE, FRACTAL_TYPE>(
                                &mut pixel.delta_current,
                                &mut pixel.jacobian_current,
                                reference.reference_data_extended[reference_index],
                                pixel.delta_reference,
                                pascal,
                                FRACTAL_POWER
                            );

                            extended_index += 1;
                            reference_index += 1;

                            next_extended_iteration = *reference.extended_iterations.get(extended_index).unwrap_or(&0xFFFFFFFF);
                        }

                        pixel.iteration += 1;
                    }

                    pixel.delta_current.reduce();

                    if DATA_TYPE == 1 || DATA_TYPE == 3 {
                        pixel.jacobian_current[0].reduce();
                        pixel.jacobian_current[1].scale_to_exponent(pixel.jacobian_current[0].exponent);
                        scale_factor_derivative = 1.0f64.ldexp(-pixel.jacobian_current[0].exponent);
                    }

                    scale_factor_delta = 1.0f64.ldexp(pixel.delta_current.exponent);
                    scaled_delta_reference = 1.0f64.ldexp(pixel.delta_reference.exponent - pixel.delta_current.exponent) * pixel.delta_reference.mantissa;
                }
            }

            data_export.lock().export_pixels::<DATA_TYPE, FRACTAL_TYPE, FRACTAL_POWER>(&pixel_data[0..pixel_index], reference, delta_pixel, scale);
            pixels_complete.fetch_add(new_pixels_complete, Ordering::Relaxed);
        });
    }
}