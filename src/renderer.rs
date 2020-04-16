use std::time::Instant;
use rand::seq::SliceRandom;
use pbr::ProgressBar;
use std::io::Stdout;
use mantexp::mantexp_complex::MantExpComplex;

use crate::util::{ComplexFixed, ComplexArbitrary};
use crate::util::point::Point;
use crate::colouring::ColourMethod;
use mantexp::mantexp::MantExp;

enum DeltaType {
    Float32,
    Float64
}

pub struct ImageRenderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: f64,
    pub maximum_iterations: usize,
    reference_points: usize,
    origin: ComplexArbitrary,
    reference: ComplexArbitrary,
    pub x_n_f64: Vec<ComplexFixed<f64>>,
    pub tolerance_check_f64: Vec<f64>,
    pub reference_delta_f64: ComplexFixed<f64>,
    glitch_tolerance: f32,
    display_glitches: bool,
    colouring_method: ColourMethod,
    progress: ProgressBar<Stdout>,
    delta_type: DeltaType,
    precision: usize,
    num_coefficients: usize,
    coefficients: Vec<Vec<MantExpComplex>>,

}

impl ImageRenderer {
    pub fn new(image_width: usize, image_height: usize, initial_zoom: f64, maximum_iterations: usize, center_re: &str, center_complex: &str, precision: usize, glitch_tolerance: f32, display_glitches: bool) -> Self {
        let location = ComplexArbitrary::with_val(
            precision as u32,
            ComplexArbitrary::parse("(".to_owned() + center_re + "," + center_complex + ")").expect("Location is not valid!"));

        let aspect = image_width as f64 / image_height as f64;
        let image_width = image_width;
        let image_height = image_height;
        let zoom = initial_zoom;

        // Set just below the limit to allow all functions to work right
        // Currently for some reason the f64 is faster than the f32 version
        let delta_type = if zoom < 1e-2 {
            println!("Delta: f32");
            DeltaType::Float32
        } else {
            println!("Delta: f64");
            DeltaType::Float64
        };

        // The height is kept to be the correct size (in terms of the zoom) and the width is scaled to counter for the aspect ratio
        // reference delta can be changes, but may need to be updated
        ImageRenderer {
            image_width,
            image_height,
            aspect,
            zoom,
            maximum_iterations,
            reference_points: 0,
            origin: location.clone(),
            reference: location,
            x_n_f64: Vec::new(),
            tolerance_check_f64: Vec::new(),
            reference_delta_f64: ComplexFixed::new(0.0, 0.0),
            glitch_tolerance,
            display_glitches,
            colouring_method: ColourMethod::Iteration,
            progress: ProgressBar::new((image_width * image_height) as u64),
            delta_type,
            precision,
            num_coefficients: 3,
            coefficients: Vec::new()
        }
    }

    pub fn render(&mut self) {
        self.progress.set(0);
        let mut points_remaining_f64 = Vec::with_capacity(self.image_width * self.image_height);
        let mut points_complete_f64 = Vec::with_capacity(self.image_width * self.image_height);

        let mut image = vec![0u8; self.image_width * self.image_height * 3];

        let resolution =  (-2.0 * (4.0 / self.image_height as f64 - 2.0) / self.zoom) / self.image_height as f64;
        let top_left_delta = ComplexFixed::new(
            (4.0 / self.image_width as f64 - 2.0) / self.zoom as f64 * self.aspect as f64,
            (4.0 / self.image_height as f64 - 2.0) / self.zoom as f64);

        for image_y in 0..self.image_height {
            for image_x in 0..self.image_width {
                points_remaining_f64.push(Point::<ComplexFixed<f64>>::new(image_x, image_y, self.image_width, resolution, top_left_delta));
            }
        }

        let mut reference_time = 0;
        let mut iteration_time = 0;

        // Start solving the points by iteratively moving around the reference point
        while (points_remaining_f64.len()) as f32 > (self.glitch_tolerance * (self.image_width * self.image_height) as f32) {
            if self.reference_points != 0 {
                // At the moment this just randomly chooses a point in the glitched points.
                // A slightly better method might be to group all of the glitched points and then choose the center of the largest group
                // and continue with this - this may be more effective at reducing the glitch count fast.
                let temp = points_remaining_f64.choose(&mut rand::thread_rng()).unwrap().delta;
                self.reference_delta_f64 = temp;
                *self.reference.mut_real() = self.origin.real().clone() + self.reference_delta_f64.re;
                *self.reference.mut_imag() = self.origin.imag().clone() + self.reference_delta_f64.im;
            }

            self.reference_points += 1;
            let start = Instant::now();
            self.calculate_reference();
            reference_time += start.elapsed().as_millis();

            let mut temp_iterations = 0;

            if self.reference_points == 1 && false {
                let skipped_iterations = self.calculate_series(top_left_delta);
                println!("Series approximation skipped: {}", skipped_iterations);
                temp_iterations = skipped_iterations;

                if skipped_iterations != 0 {
                    for point in &mut points_remaining_f64 {
                        let mut d0_to_the = Vec::with_capacity(self.num_coefficients + 1);
                        d0_to_the.push(MantExpComplex::new(
                            MantExp::new(point.delta.re, 1),
                            MantExp::new(point.delta.im, 1)));

                        d0_to_the[0].reduce();

                        for i in 1..self.num_coefficients {
                            d0_to_the.push(d0_to_the[i - 1].clone() * d0_to_the[0].clone());
                            d0_to_the[i].reduce();
                            // there is a reduce here?
                        }

                        let mut delta_sub_n = MantExpComplex::new(
                            MantExp::new(0.0, 0),
                            MantExp::new(0.0, 0));

                        for i in 0..self.num_coefficients {
                            delta_sub_n = delta_sub_n + (self.coefficients[i][skipped_iterations] * d0_to_the[i]);
                        }
                        point.delta = delta_sub_n.to_complex();
                    }
                }
            }

            // go through all points and set the delta to the new delta by series approximation

            let start = Instant::now();
            self.calculate_perturbations_f64_vectorised(&mut points_remaining_f64, &mut points_complete_f64, temp_iterations);
            iteration_time += start.elapsed().as_millis();
            self.progress.set(points_complete_f64.len() as u64);
        }

        self.progress.finish();
        println!("\n{:<10}{:>6} ms", "Reference:", reference_time);
        println!("{:<10}{:>6} ms", "Iteration:", iteration_time);
        println!("{:<12}{:>7}", "References:", self.reference_points);
        println!("{:<12}{:>7}", "Glitched:", points_remaining_f64.len());

        // Check if there are any glitched points remaining and add them to the complete points as algorithm has terminated
        if points_remaining_f64.len() > 0 {
            points_complete_f64.append(&mut points_remaining_f64);
        }

        let start = Instant::now();
        self.colouring_method.run(&points_complete_f64, &mut image, self.maximum_iterations, self.display_glitches);
        println!("{:<10}{:>6} ms", "Colouring:", start.elapsed().as_millis());

        let start = Instant::now();

        image::save_buffer("output.png", &image, self.image_width as u32, self.image_height as u32, image::RGB(8)).unwrap();

        println!("{:<10}{:>6} ms", "Saving:", start.elapsed().as_millis());
    }

    pub fn calculate_reference(&mut self) {
        let glitch_tolerance = 1e-6;

        // Clear both buffers for new values
        self.x_n_f64.clear();
        self.tolerance_check_f64.clear();

        let mut z = self.reference.clone();

        for _ in 0..=self.maximum_iterations {
            self.x_n_f64.push(ComplexFixed::new(z.real().to_f64(), z.imag().to_f64()));
            self.tolerance_check_f64.push(glitch_tolerance as f64 * self.x_n_f64.last().unwrap().norm_sqr());

            z = z.square() + &self.reference
        }
    }

    pub fn calculate_series(&mut self, top_left_delta: ComplexFixed<f64>) -> usize {
        if self.num_coefficients < 2 {
            return 0;
        }

        if self.num_coefficients > 5 {
            self.num_coefficients = 5;
        }

        let mut d0_to_the = Vec::with_capacity(self.num_coefficients + 1);
        d0_to_the.push(MantExpComplex::new(
            MantExp::new(top_left_delta.re, 1),
            MantExp::new(top_left_delta.im, 1)));

        d0_to_the[0].reduce();

        for i in 1..self.num_coefficients {
            d0_to_the.push(d0_to_the[i - 1].clone() * d0_to_the[0].clone());
            d0_to_the[i].reduce();
            // there is a reduce here?
        }

        self.coefficients = vec!(vec![MantExpComplex::new(
            MantExp::new(0.0, 1),
            MantExp::new(0.0, 1)); self.x_n_f64.len()]; self.num_coefficients);

        self.coefficients[0][0] = MantExpComplex::new(
            MantExp::new(1.0, 1),
            MantExp::new(0.0, 1));

        let tol = 2.0f64.powf(-64.0);

        for i in 1..self.x_n_f64.len() {
            if self.num_coefficients >= 1 {
                // temp here need to override f64
                self.coefficients[0][i] = (self.coefficients[0][i - 1].clone() * self.x_n_f64[i - 1] * 2.0) + MantExpComplex::new(MantExp::new(1.0, 1), MantExp::new(0.0, 1));
            }
            if self.num_coefficients >= 2 {
                self.coefficients[1][i] = (self.coefficients[1][i - 1].clone() * self.x_n_f64[i - 1] * 2.0) + self.coefficients[0][i - 1].clone() * self.coefficients[0][i - 1].clone();
            }
            if self.num_coefficients >= 3 {
                self.coefficients[2][i] = (self.coefficients[2][i - 1].clone() * self.x_n_f64[i - 1] * 2.0) + (self.coefficients[1][i - 1].clone() * self.coefficients[0][i - 1].clone() * 2.0);
            }
            if self.num_coefficients >= 4 {
                self.coefficients[3][i] = (self.coefficients[3][i - 1].clone() * self.x_n_f64[i - 1] * 2.0) + (self.coefficients[0][i - 1].clone() * self.coefficients[2][i - 1].clone() * 2.0) + self.coefficients[1][i - 1].clone() * self.coefficients[1][i - 1].clone();
            }
            if self.num_coefficients >= 5 {
                self.coefficients[4][i] = (self.coefficients[4][i - 1].clone() * self.x_n_f64[i - 1] * 2.0) + (self.coefficients[0][i - 1].clone() * self.coefficients[3][i - 1].clone() * 2.0) + (self.coefficients[1][i - 1].clone() * self.coefficients[2][i - 1].clone() * 2.0);
            }

            for j in 0..self.num_coefficients {
                self.coefficients[j][i].reduce();
            }

            if ((self.coefficients[self.num_coefficients - 2][i] * d0_to_the[self.num_coefficients - 2]).norm_mantexp() * tol).compare_lt(&mut (self.coefficients[self.num_coefficients - 1][i] * d0_to_the[self.num_coefficients - 1]).norm_mantexp()) {
                return if i <= 3 {
                    for j in 0..self.num_coefficients {
                        self.coefficients[j].truncate(0);
                    }
                    0
                } else {
                    for j in 0..self.num_coefficients {
                        self.coefficients[j].truncate(i);
                    }
                    i - 3
                }
            }
        }

        for j in 0..self.num_coefficients {
            self.coefficients[j].truncate(self.x_n_f64.len());
        }

        self.x_n_f64.len() - 1
    }
}