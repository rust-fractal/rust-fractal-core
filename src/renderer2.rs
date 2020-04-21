use crate::util::image::Image;
use crate::math::perturbation::Perturbation;
use crate::util::{ComplexArbitrary, ComplexFixed};
use crate::math::series_approximation::SeriesApproximation;
use crate::math::reference::Reference;

use pbr::ProgressBar;
use std::time::Instant;
use rand::seq::SliceRandom;

pub struct FractalRenderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: f64,
    center_location: ComplexArbitrary,
    reference_location: ComplexArbitrary,
    maximum_iteration: usize,
    approximation_order: usize,
    glitch_tolerance: f64,
    image: Image,
}

impl FractalRenderer {
    pub fn new(image_width: usize,
               image_height: usize,
               initial_zoom: f64,
               maximum_iteration: usize,
               center_real: &str,
               center_imag: &str,
               precision: usize,
               glitch_tolerance: f64,
               display_glitches: bool,
               approximation_order: usize) -> Self {

        let aspect = image_width as f64 / image_height as f64;
        let image_width = image_width;
        let image_height = image_height;
        let zoom = initial_zoom;

        let center_location = ComplexArbitrary::with_val(
            precision as u32,
            ComplexArbitrary::parse("(".to_owned() + center_real + "," + center_imag + ")").expect("Location is not valid!"));

        FractalRenderer {
            image_width,
            image_height,
            aspect,
            zoom,
            center_location: center_location.clone(),
            reference_location: center_location,
            maximum_iteration,
            approximation_order,
            glitch_tolerance,
            image: Image::new(image_width, image_height, display_glitches)
        }
    }

    pub fn render(&mut self) {
        let delta_pixel =  (-2.0 * (4.0 / self.image_height as f64 - 2.0) / self.zoom) / self.image_height as f64;
        let delta_top_left = ComplexFixed::new(
            (4.0 / self.image_width as f64 - 2.0) / self.zoom as f64 * self.aspect as f64,
            (4.0 / self.image_height as f64 - 2.0) / self.zoom as f64);

        let time = Instant::now();

        // We run the series approximation using the center point as a reference
        let mut series_approximation = SeriesApproximation::new(
            self.center_location.clone(),
            self.approximation_order,
            self.maximum_iteration,
            delta_pixel,
            delta_top_left
        );
        series_approximation.run();

        println!("{:<14}{:>6} ms", "Approximation", time.elapsed().as_millis());
        println!("{:<14}{:>6} (order {})", "Skipped", series_approximation.current_iteration, series_approximation.order);

        let time = Instant::now();
        let mut reference = series_approximation.get_reference();
        reference.run();
        println!("{:<14}{:>6} ms", "Reference", time.elapsed().as_millis());

        // Now we need to do the perturbation

        let time = Instant::now();
        let mut image_x = Vec::new();
        let mut image_y = Vec::new();
        let mut iteration = Vec::new();
        let mut delta_reference = Vec::new();
        let mut delta_current = Vec::new();
        let mut derivative_current = Vec::new();

        for i in 0..self.image_width {
            for j in 0..self.image_height {
                image_x.push(i);
                image_y.push(j);
                iteration.push(reference.start_iteration);
                let element = ComplexFixed::new(i as f64 * delta_pixel + delta_top_left.re, j as f64 * delta_pixel + delta_top_left.im);
                delta_reference.push(element);
                delta_current.push(series_approximation.evaluate(element));
                derivative_current.push(series_approximation.evaluate_derivative(element));
            }
        }

        println!("{:<14}{:>6} ms", "Packing", time.elapsed().as_millis());

        let time = Instant::now();
        let mut perturbation = Perturbation::new(image_x, image_y, iteration, delta_reference, delta_current, derivative_current);
        perturbation.iterate(&reference, reference.current_iteration);
        println!("{:<14}{:>6} ms", "Iteration", time.elapsed().as_millis());

        self.image.plot_image(&perturbation, delta_pixel);

        for i in 0..1 {
            // could sort here

            // Get all of the new glitched pixels
            let mut glitched_image_x = Vec::new();
            let mut glitched_image_y = Vec::new();
            let mut glitched_iteration = Vec::new();
            let mut glitched_delta_reference = Vec::new();
            let mut glitched_delta_current = Vec::new();
            let mut glitched_derivative_current = Vec::new();

            for i in 0..perturbation.image_x.len() {
                if perturbation.glitched[i] {
                    // reset back to the starting position
                    glitched_image_x.push(perturbation.image_x[i]);
                    glitched_image_y.push(perturbation.image_y[i]);
                    glitched_iteration.push(reference.start_iteration);
                    glitched_delta_reference.push(perturbation.delta_reference[i]);
                    glitched_delta_current.push(series_approximation.evaluate(perturbation.delta_reference[i]));
                    glitched_derivative_current.push(series_approximation.evaluate_derivative(perturbation.delta_reference[i]));
                }
            }

            println!("{:<14}{:>6}", "Glitched", glitched_image_x.len());

            // Use the series approximation to get the current reference location
            let temp = glitched_delta_reference.choose(&mut rand::thread_rng()).unwrap().clone();
            let temp2 = series_approximation.evaluate(temp);
            let temp3 = series_approximation.evaluate_derivative(temp);
            *self.reference_location.mut_real() = self.center_location.real().clone() + temp.re;
            *self.reference_location.mut_imag() = self.center_location.imag().clone() + temp.im;

            let mut reference_location_current = series_approximation.c.clone();
            *reference_location_current.mut_real() = reference_location_current.real().clone() + temp2.re;
            *reference_location_current.mut_imag() = reference_location_current.imag().clone() + temp2.im;

            // Run the new reference using the series approximation skipped values
            let mut reference2 = Reference::new(reference_location_current.clone(), reference_location_current, reference.start_iteration, self.maximum_iteration);
            reference2.run();

            for i in 0..glitched_image_x.len() {
                glitched_delta_reference[i] = glitched_delta_reference[i] - temp;
                // glitched_delta_current[i] = glitched_delta_current[i] - temp2;
            }

            let mut perturbation2 = Perturbation::new(glitched_image_x, glitched_image_y, glitched_iteration, glitched_delta_reference, glitched_delta_current, glitched_derivative_current);
            perturbation2.iterate(&reference2, reference2.current_iteration);

            // Resolving glitches
            // while points_remaining_f64.len() as f64 > 100.0 * self.glitch_tolerance * (self.image_width * self.image_height) {
            //
            // }
            self.image.plot_image(&perturbation2, delta_pixel);
        }



        let time = Instant::now();
        self.image.save();
        println!("{:<14}{:>6} ms", "Saving", time.elapsed().as_millis());
    }
}