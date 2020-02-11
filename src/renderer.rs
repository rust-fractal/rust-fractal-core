use rug::{Complex, Float};
use rayon::prelude::*;
use packed_simd::*;
use std::f64::consts::LN_2;
use crate::complex::{ComplexF64, ComplexVector};
use std::time::Instant;

pub type FloatVector = f64x8;

pub struct Point {
    pub delta: ComplexF64,
    pub index: usize,
    pub glitched: bool
}

pub struct PointVector {
    pub delta: ComplexVector,
    pub index: Vec<usize>,
    pub glitched: bool
}

impl Point {
    // start_delta is the delta from reference of point (0, 0) in the image
    pub fn new(image_x: usize, image_y: usize, image_width: usize, resolution: f64, start_delta: ComplexF64) -> Self {
        Point {
            delta: ComplexF64::new(image_x as f64 * resolution + start_delta.real, image_y as f64 * resolution + start_delta.imaginary),
            index: image_y * image_width + image_x,
            glitched: false
        }
    }
}

pub struct Renderer {
    image_width: usize,
    image_height: usize,
    zoom: f64,
    maximum_iterations: usize,
    reference_points: usize,
    reference: Complex,
    delta: ComplexF64,
    x_n: Vec<ComplexF64>,
    x_n_2: Vec<ComplexF64>,
    tolerance_check: Vec<f64>,
    start_delta: ComplexF64,
    resolution: f64
}

impl Renderer {
    pub fn new(image_width: usize, image_height: usize, zoom: f64, maximum_iterations: usize, center_real: &str, center_complex: &str, precision: u32) -> Self {
        let location_string = "(".to_owned() + center_real + "," + center_complex + ")";
        let location = Complex::with_val(precision, Complex::parse(location_string).unwrap());

        Renderer {
            image_width,
            image_height,
            zoom,
            maximum_iterations,
            reference_points: 0,
            reference: location,
            delta: ComplexF64::new(0.0, 0.0),
            x_n: Vec::new(),
            x_n_2: Vec::new(),
            tolerance_check: Vec::new(),
            start_delta: ComplexF64::new((4.0 / image_width as f64 - 2.0) / zoom, (4.0 / image_height as f64 - 2.0) / zoom),
            resolution: (-2.0 * (4.0 / image_width as f64 - 2.0) / zoom) / image_width as f64
        }
    }

    pub fn render(&mut self) {
        let mut image = vec![0u8; self.image_width * self.image_height * 3];

        let mut remaining_points = Vec::with_capacity(self.image_width * self.image_height);
//        let mut glitched_points = Vec::new();

        for image_y in 0..self.image_height {
            for image_x in 0..self.image_width {
                remaining_points.push(Point::new(image_x, image_y, self.image_width, self.resolution, self.start_delta));
            }
        }

        // This is meant to include the glitch tolerance as well
//        while remaining_points.len() > 1000 {
            self.reference_points += 1;
            self.calculate_reference();
            let iterations = self.calculate_perturbations_parallel_vectorised(&mut remaining_points);


//        }

        // go through points - complex, index

        // here the calculate_perturbations function returns only iterations and glitched for remaining_points

        let mut iteration_counts = vec![0usize; self.maximum_iterations + 1];

        for iteration in &iterations {
            iteration_counts[iteration.0.floor() as usize] += 1
        }

        for i in 1..iteration_counts.len() {
            iteration_counts[i] += iteration_counts[i - 1];
        }

        let mut total = iteration_counts[self.maximum_iterations - 1];

        for i in 0..remaining_points.len() {
//            if iterations[i].1 {
//                image[3 * i] = 255u8;
//                image[3 * i + 1] = 0u8;
//                image[3 * i + 2] = 0u8;
            /*} else */if iterations[i].0.floor() >= self.maximum_iterations as f64 {
                image[3 * i] = 0u8;
                image[3 * i + 1] = 0u8;
                image[3 * i + 2] = 0u8;
            } else {
                let test = iteration_counts[iterations[i].0.floor() as usize] as f64 / total as f64;
                let test2 = iteration_counts[iterations[i].0.floor() as usize + 1] as f64 / total as f64;

                let hue = test + (test2 - test) * iterations[i].0.fract() as f64;

                let test = (hue * 255.99) as u8;
                image[3 * i] = test;
                image[3 * i + 1] = test;
                image[3 * i + 2] = test;
            }
        }


        image::save_buffer("output.png", &image, self.image_width as u32, self.image_height as u32, image::RGB(8)).unwrap();
    }

    pub fn calculate_reference(&mut self) {
        println!("Location: \n{}", self.reference);
        let start = Instant::now();
        let glitch_tolerance = 1e-3;

        // Clear both buffers for new values
        self.x_n.clear();
        self.x_n_2.clear();

        let mut z = self.reference.clone();

        for iteration in 0..=self.maximum_iterations {
            self.x_n.push(ComplexF64::new(z.real().to_f64(), z.imag().to_f64()));
            self.x_n_2.push(2.0 * self.x_n.last().unwrap().clone());
            self.tolerance_check.push(glitch_tolerance * self.x_n.last().unwrap().norm());

            z = z.square() + &self.reference
        }
        println!("Reference: {} ms", start.elapsed().as_millis());
    }

    pub fn calculate_perturbations_single(&mut self, points: &mut Vec<Point>) -> Vec<(f64, bool)> {
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
                    (iteration as f64, glitched)
                } else {
                    (iteration as f64 + 1.0 - (z_norm.log2() / 2.0).log2(), glitched)
                }
            })
            .collect::<Vec<(f64, bool)>>();

        test
    }

    pub fn calculate_perturbations_parallel(&mut self, points: &mut Vec<Point>) -> Vec<(f64, bool)> {
        let start = Instant::now();
        let test = points.into_par_iter()
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
                                 (iteration as f64, glitched)
                             } else {
                                 (iteration as f64 + 1.0 - (z_norm.log2() / 2.0).log2(), glitched)
                             }
                         })
                         .collect::<Vec<(f64, bool)>>();

        println!("Iterating: {} ms", start.elapsed().as_millis());
        test
    }

    pub fn calculate_perturbations_parallel_vectorised(&mut self, points: &mut Vec<Point>) -> Vec<(f64, bool)> {
        // here we pack the points into vectorised points
        // possibly only works on image sizes that are multiples of the vector lanes
        // TODO investigate adding this earlier so that the vector of points never needs to be constructed

        let start = Instant::now();
        let points_vectors = (0..points.len())
            .step_by(8)
            .map(|i| {

                let delta_real = (0..8)
                    .map(|j| {
                        points[i + j].delta.real
                    })
                    .collect::<Vec<f64>>();

                let delta_imag = (0..8)
                    .map(|j| {
                        points[i + j].delta.imaginary
                    })
                    .collect::<Vec<f64>>();

                ComplexVector::new(
                    FloatVector::from_slice_unaligned(delta_real.as_slice()),
                    FloatVector::from_slice_unaligned(delta_imag.as_slice()))
            })
            .collect::<Vec<ComplexVector>>();

        println!("Packing: {} ms", start.elapsed().as_millis());
        let start = Instant::now();

        let values = points_vectors.into_par_iter()
                         .map(|point_vector| {
                             let mut delta_n = point_vector;
                             let mut iterations = f64x8::splat(0.0);
                             let mut glitched = m64x8::splat(false);
                             let mut z_norm = f64x8::splat(0.0);
                             let mut mask = m64x8::splat(true);

                             for iteration in 0..self.maximum_iterations {
                                 let temp = delta_n;
                                 delta_n += ComplexVector::splat(self.x_n_2[iteration]);
                                 delta_n *= temp;
                                 delta_n += point_vector;

                                 z_norm = mask.select((ComplexVector::splat(self.x_n[iteration + 1]) + delta_n).norm(), z_norm);
                                 mask = z_norm.le(f64x8::splat(256.0));
                                 // here keep all glitched points. We keep iterating even if all points are glitched as they might be used in the final image
                                 glitched = glitched.select(m64x8::splat(true), z_norm.le(f64x8::splat(self.tolerance_check[iteration + 1])));

                                 if mask.none() {
                                     break;
                                 }

                                 iterations += mask.select(f64x8::splat(1.0), f64x8::splat(0.0))
                             }

                             let nu = f64x8::splat(1.0) - (z_norm.ln() / f64x8::splat(2.0 * LN_2)).ln() / f64x8::splat(LN_2);
                             let out = iterations + mask.select(f64x8::splat(0.0), nu);

                             let mut test = Vec::new();

                             for i in 0..8 {
                                 test.push((out.extract(i), glitched.extract(i)))
                             }
                             test
                         })
                         .flatten()
                         .collect::<Vec<(f64, bool)>>();
        println!("Iterating: {} ms", start.elapsed().as_millis());
        values
    }
}