use rayon::prelude::*;
use packed_simd::*;
use std::time::Instant;
use rand::seq::SliceRandom;

use crate::util::{ComplexFixed, ComplexArbitrary};
use crate::util::complex_vector::ComplexVector;
use crate::util::point::Point;
use crate::colouring::ColourMethod;

// TODO add a colouring method enum
pub struct Renderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: f64,
    maximum_iterations: usize,
    reference_points: usize,
    origin: ComplexArbitrary,
    reference: ComplexArbitrary,
    x_n: Vec<ComplexFixed<f64>>,
    x_n_2: Vec<ComplexFixed<f64>>,
    tolerance_check: Vec<f64>,
    reference_delta: ComplexFixed<f64>,
    resolution: f64,
    glitch_tolerance: f64,
    display_glitches: bool,
    colouring_method: ColourMethod
}

impl Renderer {
    pub fn new(image_width: usize, image_height: usize, initial_zoom: f64, maximum_iterations: usize, center_re: &str, center_complex: &str, precision: u32, glitch_tolerance: f64, display_glitches: bool) -> Self {
        let location = ComplexArbitrary::with_val(
            precision,
            ComplexArbitrary::parse("(".to_owned() + center_re + "," + center_complex + ")").expect("Location is not valid!"));

        let aspect = image_width as f64 / image_height as f64;
        let image_width = image_width;
        let image_height = image_height;
        let zoom = initial_zoom;

        // The height is kept to be the correct size (in terms of the zoom) and the width is scaled to counter for the aspect ratio
        // reference delta can be changes, but may need to be updated
        Renderer {
            image_width,
            image_height,
            aspect,
            zoom,
            maximum_iterations,
            reference_points: 0,
            origin: location.clone(),
            reference: location,
            x_n: Vec::new(),
            x_n_2: Vec::new(),
            tolerance_check: Vec::new(),
            reference_delta: ComplexFixed::new(0.0, 0.0),
            resolution: (-2.0 * (4.0 / image_height as f64 - 2.0) / zoom) / image_height as f64,
            glitch_tolerance,
            display_glitches,
            colouring_method: ColourMethod::Histogram
        }
    }

    pub fn render(&mut self) {
        // All of these buffers have the potential to be filled with all pixels
        let mut points_remaining = Vec::with_capacity(self.image_width * self.image_height);
        let mut points_glitched = Vec::with_capacity(self.image_width * self.image_height);
        let mut points_complete = Vec::with_capacity(self.image_width * self.image_height);

        let top_left_delta = ComplexFixed::new((4.0 / self.image_width as f64 - 2.0) / self.zoom * self.aspect, (4.0 / self.image_height as f64 - 2.0) / self.zoom);

        // Fill the first buffer with the correct points and the correct offset from the origin point
        for image_y in 0..self.image_height {
            for image_x in 0..self.image_width {
                points_remaining.push(Point::new(image_x, image_y, self.image_width, self.resolution, top_left_delta));
            }
        }

        // Start solving the points by iteratively moving around the reference point
        while points_remaining.len() as f64 > (self.glitch_tolerance * (self.image_width * self.image_height) as f64) {
            if self.reference_points != 0 {
                let temp = points_remaining.choose(&mut rand::thread_rng()).unwrap().delta;
                self.reference_delta = temp;
                *self.reference.mut_real() = self.origin.real().clone() + self.reference_delta.re;
                *self.reference.mut_imag() = self.origin.imag().clone() + self.reference_delta.im;
            }

            self.reference_points += 1;
            self.calculate_reference();
            self.calculate_perturbations_parallel(&points_remaining, &mut points_glitched, &mut points_complete);

            println!("{:<12}{:>7}", "Reference:", self.reference_points);
            println!("{:<12}{:>7}", "Glitched:", points_glitched.len());

            points_remaining = points_glitched.clone();
            points_glitched.clear();
        }

        // Check if there are any glitched points remaining and add them to the complete points as algorithm has terminated
        if points_remaining.len() > 0 {
            points_complete.append(&mut points_remaining);
        }

        let start = Instant::now();
        let mut image = vec![0u8; self.image_width * self.image_height * 3];
        self.colouring_method.run(&points_complete, &mut image, self.maximum_iterations, self.display_glitches);

        println!("{:<10}{:>6} ms", "Colouring", start.elapsed().as_millis());
        let start = Instant::now();

        image::save_buffer("output.png", &image, self.image_width as u32, self.image_height as u32, image::RGB(8)).unwrap();

        println!("{:<10}{:>6} ms", "Saving", start.elapsed().as_millis());
    }

    pub fn calculate_reference(&mut self) {
        let start = Instant::now();
        let glitch_tolerance = 1e-3;

        // Clear both buffers for new values
        self.x_n.clear();
        self.x_n_2.clear();
        self.tolerance_check.clear();

        let mut z = self.reference.clone();

        for _ in 0..=self.maximum_iterations {
            self.x_n.push(ComplexFixed::new(z.real().to_f64(), z.imag().to_f64()));
            self.x_n_2.push(ComplexFixed::new(2.0, 0.0) * self.x_n.last().unwrap().clone());
            self.tolerance_check.push(glitch_tolerance * self.x_n.last().unwrap().norm_sqr());

            z = z.square() + &self.reference
        }
        println!("{:<10}{:>6} ms", "Reference", start.elapsed().as_millis());
    }

    pub fn calculate_perturbations_parallel(&mut self, remaining_points: &Vec<Point>, glitched_points: &mut Vec<Point>, finished_points: &mut Vec<Point>) {
        let start = Instant::now();
        let values = remaining_points.into_par_iter()
                         .map(|point| {
                             let delta_0 = point.delta - self.reference_delta;
                             let mut delta_n = delta_0;
                             let mut iteration = 0;
                             let mut glitched = false;
                             let mut z_norm = 0.0;

                             while z_norm < 256.0 && iteration < self.maximum_iterations {
                                 let temp = delta_n;
                                 delta_n += self.x_n_2[iteration];
                                 delta_n *= temp;
                                 delta_n += delta_0;

                                 iteration += 1;
                                 z_norm = (self.x_n[iteration] + delta_n).norm_sqr();

                                 if z_norm < self.tolerance_check[iteration] {
                                     glitched = true;
                                     break;
                                 }
                             }

                             if iteration >= self.maximum_iterations || glitched {
                                 (iteration as f64, glitched)
                             } else {
                                 (iteration as f64 + 2.0 - (z_norm.log2() / 2.0).log2() as f64, glitched)
                             }
                         })
                         .collect::<Vec<(f64, bool)>>();

        for index in 0..remaining_points.len() {
            // check to see if a point is glitched
            if values[index].1 {
                glitched_points.push(Point {
                    delta: remaining_points[index].delta,
                    index: remaining_points[index].index,
                    iterations: values[index].0,
                    glitched: true
                });
            } else {
                finished_points.push(Point {
                    delta: remaining_points[index].delta,
                    index: remaining_points[index].index,
                    iterations: values[index].0,
                    glitched: false
                });
            }
        }

        println!("{:<10}{:>6} ms", "Iterating", start.elapsed().as_millis());
    }

    pub fn calculate_perturbations_parallel_vectorised(&mut self, remaining_points: &Vec<Point>, glitched_points: &mut Vec<Point>, finished_points: &mut Vec<Point>) {
        // here we pack the points into vectorised points
        // possibly only works on image sizes that are multiples of the vector lanes
        // investigate adding this earlier so that the vector of points never needs to be constructed

        type FloatVector = f64x4;
        type MaskVector = m64x4;
        let lanes = FloatVector::lanes();

        let start = Instant::now();
        let points_vectors = (0..remaining_points.len())
            .step_by(lanes)
            .map(|i| {
                let delta_re = (0..lanes)
                    .map(|j| {
                        if i + j >= remaining_points.len() {
                            100.0
                        } else {
                            remaining_points[i + j].delta.re - self.reference_delta.re
                        }
                    })
                    .collect::<Vec<f64>>();

                let delta_im = (0..lanes)
                    .map(|j| {
                        if i + j >= remaining_points.len() {
                            100.0
                        } else {
                            remaining_points[i + j].delta.im - self.reference_delta.im
                        }
                    })
                    .collect::<Vec<f64>>();

                ComplexVector::<FloatVector>::new(delta_re.as_slice(), delta_im.as_slice())
            })
            .collect::<Vec<ComplexVector<FloatVector>>>();

        println!("{:<10}{:>6} ms", "Packing", start.elapsed().as_millis());
        let start = Instant::now();

        let values = points_vectors.into_par_iter()
                         .map(|point_vector| {
                             let mut delta_n = point_vector;
                             let mut iterations = FloatVector::splat(0.0);
                             let mut glitched = MaskVector::splat(false);
                             let mut z_norm = FloatVector::splat(0.0);
                             let mut mask = MaskVector::splat(false);

                             for iteration in 0..self.maximum_iterations {
                                 let temp = delta_n;
                                 delta_n += ComplexVector::<FloatVector>::splat(self.x_n_2[iteration]);
                                 delta_n *= temp;
                                 delta_n += point_vector;

                                 // this line takes up more than 50% of the time
                                 let new_z_norm = (ComplexVector::<FloatVector>::splat(self.x_n[iteration + 1]) + delta_n).norm_sqr();
                                 z_norm = mask.select(z_norm, new_z_norm);
                                 mask = z_norm.ge(FloatVector::splat(256.0));
                                 // here keep all glitched points. We keep iterating even if all points are glitched as they might be used in the final image

                                 // approximately 200ms / 2800ms
                                 glitched = glitched.select(MaskVector::splat(true), z_norm.le(FloatVector::splat(self.tolerance_check[iteration + 1])));

                                 if (mask | glitched).all() {
                                     break;
                                 }

                                 iterations += mask.select(FloatVector::splat(0.0), FloatVector::splat(1.0));
                             }

                             let nu = FloatVector::splat(2.0) - (z_norm.ln() / (FloatVector::splat(2.0) * FloatVector::LN_2)).ln() / FloatVector::LN_2;
                             let out = iterations + mask.select(nu, FloatVector::splat(0.0));

                             let mut test = Vec::new();

                             for i in 0..lanes {
                                 test.push((out.extract(i) as f64, glitched.extract(i)))
                             }
                             test
                         })
                         .flatten()
                         .collect::<Vec<(f64, bool)>>();

        println!("{:<10}{:>6} ms", "Iterating", start.elapsed().as_millis());

        for index in 0..remaining_points.len() {
            // check to see if a point is glitched
            if values[index].1 {
                glitched_points.push(Point {
                    delta: remaining_points[index].delta,
                    index: remaining_points[index].index,
                    iterations: values[index].0,
                    glitched: true
                });
            } else {
                finished_points.push(Point {
                    delta: remaining_points[index].delta,
                    index: remaining_points[index].index,
                    iterations: values[index].0,
                    glitched: false
                });
            }
        }

        // this takes ~2800ms at f64x8
        // 2600ms f64x4
        // 1660ms f32x4
        // 1560ms f32x8
        // 1780ms f32x16
        // this is at 2000x2000 location 0, 0 and zoom 10
    }
}