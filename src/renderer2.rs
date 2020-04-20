use crate::util::image::Image;
use crate::math::perturbation::Perturbation;
use crate::util::{ComplexArbitrary, ComplexFixed};
use crate::math::series_approximation::SeriesApproximation;
use crate::math::reference::Reference;

use pbr::ProgressBar;
use std::time::Instant;

pub struct FractalRenderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: f64,
    center_location: ComplexArbitrary,
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
            center_location,
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

        let time = Instant::now();
        self.image.plot_image(&perturbation, delta_pixel);
        self.image.save();
        println!("{:<14}{:>6} ms", "Saving", time.elapsed().as_millis());
    }
}