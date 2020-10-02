use crate::util::{to_extended, ComplexArbitrary, ComplexFixed};
use crate::util::complex_extended::ComplexExtended;
use crate::math::reference::Reference;
use rug::Float;
use crate::util::float_extended::FloatExtended;
use rayon::prelude::*;

use std::cmp::max;

pub struct SeriesApproximation {
    pub maximum_iteration: usize,
    pub delta_pixel_square: FloatExtended,
    pub order: usize,
    coefficients: Vec<Vec<ComplexExtended>>,
    probe_start: Vec<ComplexExtended>,
    approximation_probes: Vec<Vec<ComplexExtended>>,
    approximation_probes_derivative: Vec<Vec<ComplexExtended>>,
    pub min_valid_iteration: usize,
    pub valid_iterations: Vec<usize>,
    pub valid_interpolation: Vec<usize>,
    pub probe_sampling: usize,
    pub experimental: bool,
    pub valid_iteration_frame_multiplier: f32,
    pub valid_iteration_probe_multiplier: f32,
    pub high_precision_data_interval: usize,
}

impl SeriesApproximation {
    pub fn new_central(c: &ComplexArbitrary, 
        order: usize, 
        maximum_iteration: usize, 
        delta_pixel_square: FloatExtended, 
        probe_sampling: usize, 
        experimental: bool, 
        valid_iteration_frame_multiplier: f32, 
        valid_iteration_probe_multiplier: f32,
        high_precision_data_interval: usize) -> Self {

        let mut coefficients = vec![vec![ComplexExtended::new2(0.0, 0.0, 0); order as usize + 1]; 1];

        coefficients[0][0] = to_extended(&c);
        coefficients[0][1] = ComplexExtended::new2(1.0, 0.0, 0);

        // The current iteration is set to 1 as we set z = c
        SeriesApproximation {
            maximum_iteration,
            delta_pixel_square,
            order,
            coefficients,
            probe_start: Vec::new(),
            approximation_probes: Vec::new(),
            approximation_probes_derivative: Vec::new(),
            min_valid_iteration: 1,
            valid_iterations: Vec::new(),
            valid_interpolation: Vec::new(),
            probe_sampling,
            experimental,
            valid_iteration_frame_multiplier,
            valid_iteration_probe_multiplier,
            high_precision_data_interval,
        }
    }

    pub fn generate_approximation(&mut self, center_reference: &Reference) {
        // Reset the coefficients
        self.coefficients = vec![vec![ComplexExtended::new2(0.0, 0.0, 0); self.order as usize + 1]; 1];

        self.coefficients[0][0] = to_extended(&center_reference.c);
        self.coefficients[0][1] = ComplexExtended::new2(1.0, 0.0, 0);

        let add_value = ComplexExtended::new2(1.0, 0.0, 0);

        // Can be changed later into a better loop - this function could also return some more information
        for i in 1..self.maximum_iteration {
            let coefficients = &self.coefficients[i - 1];
            let mut next_coefficients = vec![ComplexExtended::new2(0.0, 0.0, 0); self.order as usize + 1];

            // This is checking if the approximation can step forward so takes the next iteration
            next_coefficients[0] = center_reference.reference_data[i].z_extended;
            next_coefficients[1] = coefficients[0] * coefficients[1] * 2.0 + add_value;
            next_coefficients[0].reduce();
            next_coefficients[1].reduce();

            // Calculate the new coefficents
            for k in 2..=self.order {
                let mut sum = coefficients[0] * coefficients[k];

                for j in 1..=((k - 1) / 2) {
                    sum += coefficients[j] * coefficients[k - j];
                }
                sum *= 2.0;

                // If even, we include the mid term as well
                if k % 2 == 0 {
                    sum += coefficients[k / 2] * coefficients[k / 2];
                }

                sum.reduce();
                next_coefficients[k] = sum;
            }

            self.coefficients.push(next_coefficients);
        }
    }

    pub fn check_approximation(&mut self, 
        delta_top_left_mantissa: ComplexFixed<f64>, 
        delta_top_left_exponent: i32, 
        cos_rotate: f64, 
        sin_rotate: f64, 
        delta_pixel: f64,
        image_width: usize,
        image_height: usize) {
        // Delete the previous probes and calculate new ones
        self.probe_start = Vec::new();
        self.approximation_probes = Vec::new();
        self.approximation_probes_derivative = Vec::new();

        // Probes are stored in row first order
        for j in 0..self.probe_sampling {
            for i in 0..self.probe_sampling {

                let pos_x = image_width as f64 * (i as f64 / (self.probe_sampling as f64 - 1.0));
                let pos_y = image_height as f64 * (j as f64 / (self.probe_sampling as f64 - 1.0));

                // This could be changed to account for jittering if needed
                let real = pos_x * delta_pixel * cos_rotate - pos_y * delta_pixel * sin_rotate + delta_top_left_mantissa.re;
                let imag = pos_x * delta_pixel * sin_rotate + pos_y * delta_pixel * cos_rotate + delta_top_left_mantissa.im;

                self.add_probe(ComplexExtended::new2(
                    cos_rotate * real - sin_rotate * imag, 
                    sin_rotate * real + cos_rotate * imag, delta_top_left_exponent));
            }
        }

        if self.experimental {
            // iterate one probe and find its value

            // iterate one probe find value
            // multiply by 0.99 or something
            // if some probes instantly have too much error, iterate them again (decrease by 0.99) and so forth until no probes remain
            // these probes only have to be iterated until they first break tolerance with an lower tolerance beforehand

            // maybe just do all at once
            // start with 0.98 and step beckwards, also maybe like has to be 100 iterations after the approximation that are valid

            // If we have no idea what the min value might be we calculate one:
            // if self.min_valid_iteration == 1 {
            if true {
                // We check using the top left probe
                let i = 0;

                let test_val = max((self.min_valid_iteration as f32 * self.valid_iteration_frame_multiplier) as usize, 1000);

                let mut first_valid_iterations = if self.min_valid_iteration > test_val {
                    self.min_valid_iteration - test_val
                } else {
                    1
                };

                let mut probe = self.evaluate(self.probe_start[i], first_valid_iterations);
                
                while first_valid_iterations < self.maximum_iteration {
                    let coefficients = &self.coefficients[first_valid_iterations - 1];
                    let next_coefficients = &self.coefficients[first_valid_iterations];

                    // step the probe points using perturbation
                    probe = probe * (coefficients[0] * 2.0 + probe);
                    probe += self.probe_start[i];

                    // This is not done on every iteration
                    if first_valid_iterations % 250 == 0 {
                        probe.reduce();
                    }

                    // get the new approximations
                    let mut series_probe = next_coefficients[1] * self.approximation_probes[i][0];
                    let mut derivative_probe = next_coefficients[1] * self.approximation_probes_derivative[i][0];

                    for k in 2..=self.order {
                        series_probe += next_coefficients[k] * self.approximation_probes[i][k - 1];
                        derivative_probe += next_coefficients[k] * self.approximation_probes_derivative[i][k - 1];
                    };

                    let relative_error = (probe - series_probe).norm_square();
                    let mut derivative = derivative_probe.norm_square();

                    // Check to make sure that the derivative is greater than or equal to 1
                    if derivative.to_float() < 1.0 {
                        derivative.mantissa = 1.0;
                        derivative.exponent = 0;
                    }

                    // The first element is reduced, the second might need to be reduced a little more
                    // Check that the error over the derivative is less than the pixel spacing
                    if relative_error / derivative > self.delta_pixel_square {
                        break;
                    }

                    first_valid_iterations += 1;
                }

                self.min_valid_iteration = first_valid_iterations;
            }

            let test_val = max((self.min_valid_iteration as f32 * self.valid_iteration_probe_multiplier) as usize, 1000);

            let mut current_probe_check_value = if self.min_valid_iteration > test_val {
                self.min_valid_iteration - test_val
            } else {
                1
            };

            let mut next_probe_check_value = if current_probe_check_value > test_val {
                current_probe_check_value - test_val
            } else {
                1
            };

            // This is the array that will be iterated
            let mut valid_iterations = vec![current_probe_check_value; self.probe_sampling * self.probe_sampling];

            // preallocate array for the values
            // TODO there needs to be checks here to make sure that we are not too close to 0

            loop {
                valid_iterations.iter_mut().enumerate()
                    .for_each(|(i, probe_iteration_level)| {
                        // let mut probe_iteration_level = valid_iteration;

                        // check if the probe has already found its max skip
                        if *probe_iteration_level == current_probe_check_value {
                            let mut probe = self.evaluate(self.probe_start[i], *probe_iteration_level);

                            while *probe_iteration_level < self.maximum_iteration {
                                let coefficients = &self.coefficients[*probe_iteration_level - 1];
                                let next_coefficients = &self.coefficients[*probe_iteration_level];

                                // step the probe points using perturbation
                                probe = probe * (coefficients[0] * 2.0 + probe);
                                probe += self.probe_start[i];

                                // This is not done on every iteration
                                if *probe_iteration_level % 250 == 0 {
                                    probe.reduce();
                                }

                                // get the new approximations
                                let mut series_probe = next_coefficients[1] * self.approximation_probes[i][0];
                                let mut derivative_probe = next_coefficients[1] * self.approximation_probes_derivative[i][0];

                                for k in 2..=self.order {
                                    series_probe += next_coefficients[k] * self.approximation_probes[i][k - 1];
                                    derivative_probe += next_coefficients[k] * self.approximation_probes_derivative[i][k - 1];
                                };

                                let relative_error = (probe - series_probe).norm_square();
                                let mut derivative = derivative_probe.norm_square();

                                // Check to make sure that the derivative is greater than or equal to 1
                                if derivative.to_float() < 1.0 {
                                    derivative.mantissa = 1.0;
                                    derivative.exponent = 0;
                                }

                                // Check that the error over the derivative is less than the pixel spacing
                                if relative_error / derivative > self.delta_pixel_square {
                                    // set the current value to the lower one so it is checked in the next iteration
                                    if *probe_iteration_level <= (current_probe_check_value + 10) {
                                        *probe_iteration_level = next_probe_check_value;
                                    };
                                    
                                    break;
                                }

                                *probe_iteration_level += 1;
                            }
                        };
                    });

                // we have now iterated all the values, we need to update those which skipped too quickly

                self.min_valid_iteration = valid_iterations.iter().min().unwrap().clone();

                // this would indicate that no more of the probes are bad
                if self.min_valid_iteration != next_probe_check_value  || self.min_valid_iteration == 1 {
                    // println!("{:?}, {}, {}", valid_iterations, self.min_valid_iteration, current_probe_check_value);
                    break;
                } else {
                    current_probe_check_value = next_probe_check_value;

                    let test_val = max((self.min_valid_iteration as f32 * self.valid_iteration_probe_multiplier) as usize, 1000);
        
                    next_probe_check_value = if current_probe_check_value > test_val {
                        current_probe_check_value - test_val
                    } else {
                        1
                    };
                    // println!("{:?}, {}, {}, {}", valid_iterations, self.min_valid_iteration, current_probe_check_value, next_probe_check_value);
                }
            }

            self.valid_iterations = valid_iterations;

            // Also, here we do the interpolation and set up the array
            self.valid_interpolation = Vec::new();

            for j in 0..(self.probe_sampling - 1) {
                for i in 0..(self.probe_sampling - 1) {
                    // this is the index into the main array
                    let index = j * self.probe_sampling + i;

                    let min_interpolation = [self.valid_iterations[index], 
                        self.valid_iterations[index + 1], 
                        self.valid_iterations[index + self.probe_sampling], 
                        self.valid_iterations[index + self.probe_sampling + 1]].iter().min().unwrap().clone();

                    self.valid_interpolation.push(min_interpolation);
                }
            }

        } else {
            // Possible how to add the roots as probes

            // lets say we check the first 1000 iterations for roots with periods 0-1000

            // Each iteration, we NR the current SA polynomial??

            // Maybe only needs to be done once for the entire zoom sequence??

            // How do we know which probes to use??

            // Maybe some kind of distance check, where the probes are chosen that are in the image, and are within some tolerance
            // ~10 pixel spacing from others. 

            // Possible loop once the roots have been calculated

            // Get all of the root locations, and don't use any of them that are not in the image
            // There might still be quite a few roots, so remove those which are close to each other

            self.valid_iterations = (0..self.probe_start.len()).into_par_iter()
                .map(|i| {
                    let mut valid_iterations = 1;
                    let mut probe = self.probe_start[i];

                    while valid_iterations < self.maximum_iteration {
                        let coefficients = &self.coefficients[valid_iterations - 1];
                        let next_coefficients = &self.coefficients[valid_iterations];

                        // step the probe points using perturbation
                        probe = probe * (coefficients[0] * 2.0 + probe);
                        probe += self.probe_start[i];

                        // This is not done on every iteration
                        if valid_iterations % 250 == 0 {
                            probe.reduce();
                        }

                        // get the new approximations
                        let mut series_probe = next_coefficients[1] * self.approximation_probes[i][0];
                        let mut derivative_probe = next_coefficients[1] * self.approximation_probes_derivative[i][0];

                        for k in 2..=self.order {
                            series_probe += next_coefficients[k] * self.approximation_probes[i][k - 1];
                            derivative_probe += next_coefficients[k] * self.approximation_probes_derivative[i][k - 1];
                        };

                        let relative_error = (probe - series_probe).norm_square();
                        let mut derivative = derivative_probe.norm_square();

                        // Check to make sure that the derivative is greater than or equal to 1
                        if derivative.to_float() < 1.0 {
                            derivative.mantissa = 1.0;
                            derivative.exponent = 0;
                        }

                        // Check that the error over the derivative is less than the pixel spacing
                        if relative_error / derivative > self.delta_pixel_square {
                            break;
                        }

                        valid_iterations += 1;
                        
                    }

                valid_iterations
            }).collect::<Vec<usize>>();

            self.min_valid_iteration = self.valid_iterations.iter().min().unwrap().clone();
        }
    }

    pub fn add_probe(&mut self, delta_probe: ComplexExtended) {
        // here we will need to check to make sure we are still at the first iteration, or use perturbation to go forward
        self.probe_start.push(delta_probe);

        let mut current_value = delta_probe;

        let mut delta_n = Vec::with_capacity(self.order + 1);
        let mut delta_derivative_n = Vec::with_capacity(self.order + 1);

        // The first element will be 1, in order for the derivative to be calculated
        delta_n.push(current_value);
        delta_derivative_n.push(ComplexExtended::new2(1.0, 0.0, 0));

        for i in 1..=self.order {
            delta_derivative_n.push(current_value * (i + 1) as f64);
            current_value *= delta_probe;
            delta_n.push(current_value);
        }

        self.approximation_probes.push(delta_n);
        self.approximation_probes_derivative.push(delta_derivative_n);
    }

    // Get the current reference, and the current number of iterations done
    pub fn get_reference(&self, reference_delta: ComplexExtended, center_reference: &Reference) -> Reference {
        let precision = center_reference.c.real().prec();
        let iteration_reference = self.high_precision_data_interval * ((self.min_valid_iteration - 1) / self.high_precision_data_interval) + 1;

        let mut reference_c = center_reference.c.clone();
        let temp = Float::with_val(precision, reference_delta.exponent).exp2();
        let temp2 = Float::with_val(precision, reference_delta.mantissa.re);
        let temp3 = Float::with_val(precision, reference_delta.mantissa.im);

        *reference_c.mut_real() += &temp2 * &temp;
        *reference_c.mut_imag() += &temp3 * &temp;

        // let mut reference_z = self.center_reference.approximation_data[self.valid_iteration].clone();
        let mut reference_z = center_reference.high_precision_data[(self.min_valid_iteration - 1) / self.high_precision_data_interval].clone();

        let temp4 = self.evaluate(reference_delta, iteration_reference);
        let temp = Float::with_val(precision, temp4.exponent).exp2();
        let temp2 = Float::with_val(precision, temp4.mantissa.re);
        let temp3 = Float::with_val(precision, temp4.mantissa.im);

        *reference_z.mut_real() += &temp2 * &temp;
        *reference_z.mut_imag() += &temp3 * &temp;

        Reference::new(reference_z, reference_c, iteration_reference, center_reference.maximum_iteration, self.high_precision_data_interval, center_reference.glitch_tolerance)
    }

    pub fn evaluate(&self, point_delta: ComplexExtended, iteration: usize) -> ComplexExtended {
        // This could be improved to use the iteration option better

        let new_coefficients = &self.coefficients[iteration - 1];
        // Horner's rule
        let mut approximation = new_coefficients[self.order];

        for k in (1..=(self.order - 1)).rev() {
            approximation *= point_delta;
            approximation += new_coefficients[k];
        }

        approximation *= point_delta;
        approximation.reduce();
        approximation
    }

    // pub fn evaluate_derivative(&self, point_delta: ComplexExtended) -> FloatExtended {
    //     let mut original_point_derivative_n = ComplexExtended::new(1.0, 0, 0.0, 0);
    //     let mut approximation_derivative = ComplexExtended::new(0.0, 0, 0.0, 0);
    
    //     for k in 1..=self.order {
    //         approximation_derivative += k as f64 * self.coefficients[k] * original_point_derivative_n;
    //         original_point_derivative_n *= ComplexExtended::new(point_delta.re, 0, point_delta.im, 0);
    //     };
    
    //     approximation_derivative.to_float()


    //     let mut approximation_derivative = self.coefficients[self.order] * self.order as f64;

    //     for k in (1..=(self.order - 1)).rev() {
    //         approximation *= point_delta;
    //         approximation += self.coefficients[k] * k as f64;
    //     }

    //     approximation *= point_delta;
    //     approximation.reduce();
    //     approximation


    // }
}