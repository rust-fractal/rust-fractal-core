use rug::Complex;
use rayon::prelude::*;
use packed_simd::*;
use std::f64::consts::LN_2;
use crate::complex::ComplexF64;

pub struct Point {
    pub delta: ComplexF64,
    pub index: usize,
    pub glitched: bool
}

impl Point {
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
            let iterations = self.calculate_perturbations(&mut remaining_points);


//        }

        // go through points - complex, index

        let mut iteration_counts = vec![0usize; self.maximum_iterations + 1];

        for iteration in &iterations {
            iteration_counts[iteration.0.floor() as usize] += 1
        }

        for i in 1..iteration_counts.len() {
            iteration_counts[i] += iteration_counts[i - 1];
        }

        let mut total = iteration_counts[self.maximum_iterations - 1];

        for i in 0..remaining_points.len() {
            if iterations[i].0.floor() >= self.maximum_iterations as f64 {
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
        println!("Calculating reference: {}", self.reference);
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
    }

    pub fn calculate_perturbations(&mut self, points: &mut Vec<Point>) -> Vec<(f64, bool)> {
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

        test
    }
}