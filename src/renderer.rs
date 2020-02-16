use rayon::prelude::*;
use packed_simd::*;
use std::time::Instant;
use rand::seq::SliceRandom;
use pbr::ProgressBar;
use std::io::Stdout;

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
    pub maximum_iterations: usize,
    reference_points: usize,
    origin: ComplexArbitrary,
    reference: ComplexArbitrary,
    pub x_n: Vec<ComplexFixed<f64>>,
    pub x_n_2: Vec<ComplexFixed<f64>>,
    pub tolerance_check: Vec<f64>,
    pub reference_delta: ComplexFixed<f64>,
    resolution: f64,
    glitch_tolerance: f64,
    display_glitches: bool,
    colouring_method: ColourMethod,
    progress: ProgressBar<Stdout>
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
            colouring_method: ColourMethod::Iteration,
            progress: ProgressBar::new((image_width * image_height) as u64)
        }
    }

    pub fn render(&mut self) {
        self.progress.set(0);
        // All of these buffers have the potential to be filled with all pixels
        let mut points_remaining = Vec::with_capacity(self.image_width * self.image_height);
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
                // At the moment this just randomly chooses a point in the glitched points.
                // A slightly better method might be to group all of the glitched points and then choose the center of the largest group
                // and continue with this - this may be more effective at reducing the glitch count fast.
                let temp = points_remaining.choose(&mut rand::thread_rng()).unwrap().delta;
                self.reference_delta = temp;
                *self.reference.mut_real() = self.origin.real().clone() + self.reference_delta.re;
                *self.reference.mut_imag() = self.origin.imag().clone() + self.reference_delta.im;
            }

            self.reference_points += 1;
            self.calculate_reference();
            self.calculate_perturbations_f64(&mut points_remaining, &mut points_complete);
            self.progress.set(points_complete.len() as u64);
        }

        self.progress.finish();
        println!("\n{:<12}{:>7}", "References:", self.reference_points);
        println!("{:<12}{:>7}", "Glitched:", points_remaining.len());

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
        let glitch_tolerance = 1e-6;

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
//        println!("{:<10}{:>6} ms", "Reference", start.elapsed().as_millis());
    }
}