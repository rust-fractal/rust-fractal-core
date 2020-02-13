use rug::Complex;
use rayon::prelude::*;
use packed_simd::*;
use crate::complex::{Complexf32, ComplexVector};
use std::time::Instant;
use rand::seq::SliceRandom;

pub type FloatVector = f32x8;

#[derive(Copy, Clone)]
pub struct Point {
    pub delta: Complexf32,
    pub index: usize,
    pub iterations: f32,
    pub glitched: bool
}

pub struct PointVector {
    pub delta: ComplexVector,
    pub index: Vec<usize>,
    pub glitched: bool
}

impl Point {
    // start_delta is the delta from reference of point (0, 0) in the image
    pub fn new(image_x: usize, image_y: usize, image_width: usize, resolution: f32, top_left_delta: Complexf32) -> Self {
        Point {
            delta: Complexf32::new(image_x as f32 * resolution + top_left_delta.real, image_y as f32 * resolution + top_left_delta.imaginary),
            index: image_y * image_width + image_x,
            iterations: 0.0,
            glitched: false
        }
    }
}

// TODO add a colouring method enum
pub struct Renderer {
    image_width: usize,
    image_height: usize,
    maximum_iterations: usize,
    reference_points: usize,
    origin: Complex,
    reference: Complex,
    x_n: Vec<Complexf32>,
    x_n_2: Vec<Complexf32>,
    tolerance_check: Vec<f32>,
    top_left_delta: Complexf32,
    reference_delta: Complexf32,
    resolution: f32,
}

impl Renderer {
    pub fn new(image_width: usize, image_height: usize, zoom: f32, maximum_iterations: usize, center_real: &str, center_complex: &str, precision: u32) -> Self {
        let location = Complex::with_val(
            precision,
            Complex::parse("(".to_owned() + center_real + "," + center_complex + ")").expect("Location is not valid!"));

        let aspect = image_width as f32 / image_height as f32;
        let image_width = image_width;
        let image_height = image_height;

        // The height is kept to be the correct size (in terms of the zoom) and the width is scaled to counter for the aspect ratio
        // reference delta can be changes, but may need to be updated
        Renderer {
            image_width,
            image_height,
            maximum_iterations,
            reference_points: 0,
            origin: location.clone(),
            reference: location,
            x_n: Vec::new(),
            x_n_2: Vec::new(),
            tolerance_check: Vec::new(),
            top_left_delta: Complexf32::new((4.0 / image_width as f32 - 2.0) / zoom * aspect, (4.0 / image_height as f32 - 2.0) / zoom),
            reference_delta: Complexf32::new(0.0, 0.0),
            resolution: (-2.0 * (4.0 / image_height as f32 - 2.0) / zoom) / image_height as f32,
        }
    }

    pub fn render(&mut self) {
        let mut remaining_points = Vec::with_capacity(self.image_width * self.image_height);
        let mut glitched_points = Vec::with_capacity(self.image_width * self.image_height);
        let mut finished_points = Vec::with_capacity(self.image_width * self.image_height);

//        let mut glitched_points = Vec::new();

        for image_y in 0..self.image_height {
            for image_x in 0..self.image_width {
                remaining_points.push(Point::new(image_x, image_y, self.image_width, self.resolution, self.top_left_delta));
            }
        }

        println!("Top left delta: {}, {}", self.top_left_delta.real, self.top_left_delta.imaginary);

        // This is meant to include the glitch tolerance as well
        while remaining_points.len() > (0.001 * self.image_width as f32 * self.image_height as f32) as usize {
            println!("Reference: \n{}", self.reference);
            println!("Reference Delta: {} {}", self.reference_delta.real, self.reference_delta.imaginary);

            self.reference_points += 1;
            self.calculate_reference();
            self.calculate_perturbations_parallel(&remaining_points, &mut glitched_points, &mut finished_points);

            println!("glitched {}", glitched_points.len());
            println!("total {}", remaining_points.len());


            remaining_points = glitched_points.clone();
            glitched_points.clear();

            // move to new reference location

            if remaining_points.len() > 0 {
                let temp = remaining_points.choose(&mut rand::thread_rng()).unwrap().delta;

                // point delta - reference delta

                self.reference_delta = temp;
                *self.reference.mut_real() = self.origin.real().clone() + self.reference_delta.real;
                *self.reference.mut_imag() = self.origin.imag().clone() + self.reference_delta.imaginary;
            }

        }

        // there are glitched points left
        if remaining_points.len() > 0 {
            println!("glitched points in image: {}", remaining_points.len());
            finished_points.append(&mut remaining_points);
        }

        // go through points - complex, index

        // here the calculate_perturbations function returns only iterations and glitched for remaining_points
        let start = Instant::now();
        let mut iteration_counts = vec![0usize; self.maximum_iterations + 1];

        for point in &finished_points {
            iteration_counts[point.iterations.floor() as usize] += 1
        }

        for i in 1..iteration_counts.len() {
            iteration_counts[i] += iteration_counts[i - 1];
        }

//        let total = iteration_counts[self.maximum_iterations - 1];

        let mut colours = Vec::new();

        for i in 0..8192 {
            let value = i as f32 / 8192 as f32;

            let red;
            let green;
            let blue;

            if value < 0.16 {
                let factor = (value - 0.0) / (0.16 - 0.0);

                red = 0.0 + factor * (32.0 - 0.0);
                green = 7.0 + factor * (107.0 - 7.0);
                blue = 100.0 + factor * (203.0 - 100.0);
            } else if value < 0.42 {
                let factor = (value - 0.16) / (0.42 - 0.16);

                red = 32.0 + factor * (237.0 - 32.0);
                green = 107.0 + factor * (255.0 - 107.0);
                blue = 203.0 + factor * (255.0 - 203.0);
            } else if value < 0.6425 {
                let factor = (value - 0.42) / (0.6425 - 0.42);

                red = 237.0 + factor * (255.0 - 237.0);
                green = 255.0 + factor * (170.0 - 255.0);
                blue = 255.0 + factor * (0.0 - 255.0);
            } else if value < 0.8575 {
                let factor = (value - 0.6425) / (0.8575 - 0.6425);

                red = 255.0 + factor * (0.0 - 255.0);
                green = 170.0 + factor * (2.0 - 170.0);
                blue = 0.0 + factor * (0.0 - 0.0);
            } else {
                let factor = (value - 0.8575) / (1.0 - 0.8575);

                red = 0.0 + factor * (0.0 - 0.0);
                green = 2.0 + factor * (7.0 - 2.0);
                blue = 0.0 + factor * (100.0 - 0.0);
            }

            colours.push((red, green, blue))
        }

        let mut image = vec![0u8; self.image_width * self.image_height * 3];

        let show_glitches = true;

        for point in finished_points {
            let index = point.index;

            if point.glitched && show_glitches {
                image[3 * index] = 255u8;
                image[3 * index + 1] = 0u8;
                image[3 * index + 2] = 0u8;
            } else if point.iterations.floor() >= self.maximum_iterations as f32 {
                image[3 * index] = 0u8;
                image[3 * index + 1] = 0u8;
                image[3 * index + 2] = 0u8;
            } else {
//                let v1 = iteration_counts[iterations[i].0.floor() as usize] as f32 / total as f32;
//                let v2 = iteration_counts[iterations[i].0.floor() as usize + 1] as f32 / total as f32;

                // the hue is used to smooth the histogram bins. The hue is in the range 0.0-1.0
//                let hue = (v1 + (v2 - v1) * iterations[i].0.fract() as f32) * 8192.0;
                let hue = ((point.iterations + 1.0).sqrt() * 1600.0) % 8192.0;

                let colour = colours[(hue.floor() as usize) % 8192];
                let colour2 = colours[(hue.floor() as usize + 1) % 8192];

                let red = (colour.0 + ((colour2.0 - colour.0) * hue.fract())) as u8;
                let green = (colour.1 + ((colour2.1 - colour.1) * hue.fract())) as u8;
                let blue = (colour.2 + ((colour2.2 - colour.2) * hue.fract())) as u8;

                image[3 * index] = red;
                image[3 * index + 1] = green;
                image[3 * index + 2] = blue;
            }
        }
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
            self.x_n.push(Complexf32::new(z.real().to_f32(), z.imag().to_f32()));
            self.x_n_2.push(2.0 * self.x_n.last().unwrap().clone());
            self.tolerance_check.push(glitch_tolerance * self.x_n.last().unwrap().norm());

            z = z.square() + &self.reference
        }
        println!("{:<10}{:>6} ms", "Reference", start.elapsed().as_millis());
    }

    pub fn calculate_perturbations_single(&mut self, points: &Vec<Point>) -> Vec<(f32, bool)> {
        let test = points.into_iter()
            .map(|point| {
                let mut delta_n = point.delta;
                let mut iteration = 0;
                let mut glitched = false;
                let mut z_norm = 0.0;

                while z_norm < 256.0 && iteration < self.maximum_iterations {
                    delta_n = self.x_n_2[iteration] * delta_n + delta_n * delta_n + point.delta;

                    iteration += 1;
                    z_norm = (self.x_n[iteration] + delta_n).norm();

//                     glitch detection
                    if z_norm < self.tolerance_check[iteration] {
                        glitched = true;
                        break;
                    }
                }

                if iteration >= self.maximum_iterations || glitched {
                    (iteration as f32, glitched)
                } else {
                    (iteration as f32 + 1.0 - (z_norm.log2() / 2.0).log2(), glitched)
                }
            })
            .collect::<Vec<(f32, bool)>>();

        test
    }

    pub fn calculate_perturbations_parallel(&mut self, remaining_points: &Vec<Point>, glitched_points: &mut Vec<Point>, finished_points: &mut Vec<Point>) {
        let start = Instant::now();
        let values = remaining_points.into_par_iter()
                         .map(|point| {
                             let mut delta_n = point.delta - self.reference_delta;
                             let mut iteration = 0;
                             let mut glitched = false;
                             let mut z_norm = 0.0;

                             while z_norm < 256.0 && iteration < self.maximum_iterations {
                                 delta_n = self.x_n_2[iteration] * delta_n + delta_n * delta_n + (point.delta - self.reference_delta);

                                 iteration += 1;
                                 z_norm = (self.x_n[iteration] + delta_n).norm();

//                     glitch detection
                                 if z_norm < self.tolerance_check[iteration] {
                                     glitched = true;
                                     break;
                                 }
                             }

                             if iteration >= self.maximum_iterations || glitched {
                                 (iteration as f32, glitched)
                             } else {
                                 (iteration as f32 + 2.0 - (z_norm.log2() / 2.0).log2(), glitched)
                             }
                         })
                         .collect::<Vec<(f32, bool)>>();

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
        // TODO investigate adding this earlier so that the vector of points never needs to be constructed

        let start = Instant::now();
        let points_vectors = (0..remaining_points.len())
            .step_by(8)
            .map(|i| {
                let delta_real = (0..8)
                    .map(|j| {
                        if i + j >= remaining_points.len() {
                            100.0
                        } else {
                            remaining_points[i + j].delta.real - self.reference_delta.real
                        }
                    })
                    .collect::<Vec<f32>>();

                let delta_imag = (0..8)
                    .map(|j| {
                        if i + j >= remaining_points.len() {
                            100.0
                        } else {
                            remaining_points[i + j].delta.imaginary - self.reference_delta.imaginary
                        }
                    })
                    .collect::<Vec<f32>>();

                ComplexVector::new(
                    FloatVector::from_slice_unaligned(delta_real.as_slice()),
                    FloatVector::from_slice_unaligned(delta_imag.as_slice()))
            })
            .collect::<Vec<ComplexVector>>();

        println!("{:<10}{:>6} ms", "Packing", start.elapsed().as_millis());
        let start = Instant::now();

        let values = points_vectors.into_par_iter()
                         .map(|point_vector| {
                             let mut delta_n = point_vector;
                             let mut iterations = f32x8::splat(0.0);
                             let mut glitched = m32x8::splat(false);
                             let mut z_norm = f32x8::splat(0.0);
                             let mut mask = m32x8::splat(true);

                             for iteration in 0..self.maximum_iterations {
                                 let temp = delta_n;
                                 delta_n += ComplexVector::splat(self.x_n_2[iteration]);
                                 delta_n *= temp;
                                 delta_n += point_vector;

                                 // this line takes up more than 50% of the time
                                 let new_z_norm = (ComplexVector::splat(self.x_n[iteration + 1]) + delta_n).norm();
                                 z_norm = mask.select(new_z_norm, z_norm);
                                 mask = z_norm.le(f32x8::splat(256.0));
                                 // here keep all glitched points. We keep iterating even if all points are glitched as they might be used in the final image

                                 // approximately 200ms / 2800ms
                                 glitched = glitched.select(m32x8::splat(true), z_norm.le(f32x8::splat(self.tolerance_check[iteration])));

                                 if mask.none() {
                                     break;
                                 }

                                 iterations += mask.select(f32x8::splat(1.0), f32x8::splat(0.0))
                             }

                             let nu = f32x8::splat(2.0) - (z_norm.ln() / (f32x8::splat(2.0) * f32x8::LN_2)).ln() / f32x8::LN_2;
                             let out = iterations + mask.select(f32x8::splat(0.0), nu);

                             let mut test = Vec::new();

                             for i in 0..8 {
                                 test.push((out.extract(i) as f32, glitched.extract(i)))
                             }
                             test
                         })
                         .flatten()
                         .collect::<Vec<(f32, bool)>>();

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