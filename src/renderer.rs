use rug::Complex;
use rayon::prelude::*;
use packed_simd::*;
use std::f64::consts::LN_2;

pub struct Point {
    pub image_x: usize,
    pub image_y: usize,
    pub iterations: f64,
    pub glitched: bool
}

impl Point {
    pub fn new(image_x: usize, image_y: usize) -> Self {
        Point {
            image_x,
            image_y,
            iterations: 0.0,
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
    delta: (f64, f64),
    x_n: Vec<(f64, f64)>,
    x_n_2: Vec<(f64, f64)>,
    tolerance_check: Vec<f64>,
    start_delta: (f64, f64),
    resolution: f64
}

impl Renderer {
    pub fn new(image_width: usize, image_height: usize, zoom: f64, maximum_iterations: usize, center_real: &str, center_complex: &str, precision: u32) -> Self {
        let location_string = "(".to_owned() + center_real + "," + center_complex + ")";
        let location = Complex::with_val(precision,
        Complex::parse(location_string).unwrap()
        );

        Renderer {
            image_width,
            image_height,
            zoom,
            maximum_iterations,
            reference_points: 0,
            reference: location,
            delta: (0.0, 0.0),
            x_n: Vec::new(),
            x_n_2: Vec::new(),
            tolerance_check: Vec::new(),
            start_delta: ((4.0 / image_width as f64 - 2.0) / zoom, (4.0 / image_height as f64 - 2.0) / zoom),
            resolution: (-2.0 * (4.0 / image_width as f64 - 2.0) / zoom) / image_width as f64
        }
    }

    pub fn render(&mut self) {
        let mut image = vec![0u8; self.image_width * self.image_height * 3];
        let mut iterations = vec![0f64; self.image_width * self.image_height];

        let mut remaining_points = Vec::with_capacity(self.image_width * self.image_height);
//        let mut glitched_points = Vec::new();

        for image_y in 0..self.image_height {
            for image_x in 0..self.image_width {
                remaining_points.push(Point::new(image_x, image_y));
            }
        }

        // This is meant to include the glitch tolerance as well
//        while remaining_points.len() > 1000 {
            self.reference_points += 1;
            self.calculate_reference();
            self.calculate_perturbations(&mut remaining_points);


//        }

        let mut iteration_counts = vec![0usize; self.maximum_iterations + 1];

        for point in &remaining_points {
            iteration_counts[point.iterations.floor() as usize] += 1
        }

        for i in 1..iteration_counts.len() {
            iteration_counts[i] += iteration_counts[i - 1];
        }

        let mut total = iteration_counts[self.maximum_iterations - 1];

        for i in 0..remaining_points.len() {
            if remaining_points[i].iterations.floor() >= self.maximum_iterations as f64 {
                image[3 * i] = 0u8;
                image[3 * i + 1] = 0u8;
                image[3 * i + 2] = 0u8;
            } else {
                let test = iteration_counts[remaining_points[i].iterations.floor() as usize] as f64 / total as f64;
                let test2 = iteration_counts[remaining_points[i].iterations.floor() as usize + 1] as f64 / total as f64;

                let hue = test + (test2 - test) * remaining_points[i].iterations.fract() as f64;

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
            self.x_n.push((z.real().to_f64(), z.imag().to_f64()));
            self.x_n_2.push((2.0 * z.real().to_f64(), 2.0 * z.imag().to_f64()));
            self.tolerance_check.push(glitch_tolerance * (z.real().to_f64().powi(2) + z.imag().to_f64().powi(2)));

            z = z.square() + &self.reference
        }
    }

    pub fn calculate_perturbations(&mut self, points: &mut Vec<Point>) {
        let block_size = f64x8::lanes();
        let width_in_blocks = self.image_width / block_size;

        points.as_mut_slice()
            .chunks_mut(block_size)
            .for_each(|test| {
                // test should be a mutable slice of points
                // here we get the x coordinate of each element in the slice

                let mut test2 = (0..block_size)
                    .map(|j| {
                        test[j].image_x as f64 * self.resolution + self.start_delta.0
                    }).collect::<Vec<f64>>();

                let delta_0_real = f64x8::from_slice_unaligned(test2.as_slice());

                // repeat for the y coordinates

                let mut test2 = (0..block_size)
                    .map(|j| {
                        test[j].image_y as f64 * self.resolution + self.start_delta.1
                    }).collect::<Vec<f64>>();

                let delta_0_imag = f64x8::from_slice_unaligned(test2.as_slice());


                let mut count = f64x8::splat(0.0);
                let mut norm = f64x8::splat(0.0);
                let mut mask = m64x8::splat(true);

                let mut delta_n_real = delta_0_real;
                let mut delta_n_imag = delta_0_imag;

                let mut iteration = 0;

                while iteration < self.maximum_iterations {
                    let delta_n_real_square = delta_n_real.powf(f64x8::splat(2.0));
                    let delta_n_imag_square = delta_n_imag.powf(f64x8::splat(2.0));
                    let delta_n_cross = f64x8::splat(2.0) * delta_n_real * delta_n_imag;

                    delta_n_real = f64x8::splat(self.x_n_2[iteration].0) * delta_n_real - f64x8::splat(self.x_n_2[iteration].1) * delta_n_imag + delta_n_real_square - delta_n_imag_square + delta_0_real;
                    delta_n_imag = f64x8::splat(self.x_n_2[iteration].1) * delta_n_real + f64x8::splat(self.x_n_2[iteration].0) * delta_n_imag + delta_n_cross + delta_0_imag;

                    norm = mask.select(
                        (self.x_n[iteration].0 + delta_n_real).powf(f64x8::splat(2.0)) + (self.x_n[iteration].1 + delta_n_imag).powf(f64x8::splat(2.0)),
                        norm);

                    iteration += 1;

                    mask = norm.le(f64x8::splat(256.0));
                    count += mask.select(f64x8::splat(1.0), f64x8::splat(0.0));

                    if mask.none() {
                        break;
                    }
                }

                let nu = f64x8::splat(2.0) - (norm.ln() / f64x8::splat(2.0 * LN_2)).ln() / f64x8::splat(LN_2);

                let out = count + mask.select(f64x8::splat(0.0), nu);

                for i in 0..test.len() {
                    test[i].iterations = out.extract(i)
                }
            });

//        for (i, element) in temp.into_iter().enumerate() {
//            element.write_to_slice_unaligned(&mut iterations[(block_size * i)..(block_size * (i + 1))]);
//        }
    }
}