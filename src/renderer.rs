use crate::util::image::Image;
use crate::util::{ComplexArbitrary, ComplexFixed, PixelData, PixelDataDouble};
use crate::math::series_approximation_extended::SeriesApproximationExtended;

use std::time::Instant;
use rand::seq::SliceRandom;
use std::cmp::max;
use crate::util::float_extended::FloatExtended;
use std::f64::consts::LOG2_10;
use crate::util::complex_extended::ComplexExtended;
use itertools::Itertools;
use rayon::prelude::*;
use crate::math::series_approximation_double::SeriesApproximationDouble;
use crate::math::perturbation::Perturbation;
use crate::math::perturbation_double::PerturbationDouble;
use crate::util::colouring_double::ColouringDouble;
use crate::util::colouring_extended::ColouringExtended;

pub struct FractalRenderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: FloatExtended,
    center_location: ComplexArbitrary,
    maximum_iteration: usize,
    approximation_order: usize,
    glitch_tolerance: f64,
    image: Image,
}

impl FractalRenderer {
    pub fn new(image_width: usize,
               image_height: usize,
               initial_zoom: &str,
               maximum_iteration: usize,
               center_real: &str,
               center_imag: &str,
               glitch_tolerance: f64,
               display_glitches: bool,
               approximation_order: usize) -> Self {

        let aspect = image_width as f64 / image_height as f64;
        let image_width = image_width;
        let image_height = image_height;
        let temp: Vec<&str> = initial_zoom.split('E').collect();
        let zoom = FloatExtended::new(temp[0].parse::<f64>().unwrap() * 2.0_f64.powf((temp[1].parse::<f64>().unwrap() * LOG2_10).fract()), (temp[1].parse::<f64>().unwrap() * LOG2_10).floor() as i32);

        let delta_pixel =  (-2.0 * (4.0 / image_height as f64 - 2.0) / zoom) / image_height as f64;
        let radius = delta_pixel * image_width as f64;
        let precision = max(64, -radius.exponent + 64);

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
        let delta_pixel =  (-2.0 * (4.0 / self.image_height as f64 - 2.0) / self.zoom.mantissa) / self.image_height as f64;

        // this should be the delta relative to the image, without the big zoom factor applied.
        let delta_top_left = ComplexFixed::new((4.0 / self.image_width as f64 - 2.0) / self.zoom.mantissa * self.aspect as f64, (4.0 / self.image_height as f64 - 2.0) / self.zoom.mantissa);

        let time = Instant::now();

        println!("Rendering...");

        let mut series_approximation = SeriesApproximationExtended::new(
            self.center_location.clone(),
            self.approximation_order,
            self.maximum_iteration,
            FloatExtended::new(delta_pixel, -self.zoom.exponent),
            ComplexExtended::new(delta_top_left, -self.zoom.exponent),
        );

        series_approximation.run();

        println!("{:<14}{:>6} ms", "Approximation", time.elapsed().as_millis());
        println!("{:<16}{:>6} (order {})", "Skipped", series_approximation.current_iteration, series_approximation.order);

        let time = Instant::now();

        let mut reference = series_approximation.get_reference(ComplexExtended::new2(0.0, 0.0, 0));
        reference.run();

        println!("{:<14}{:>6} ms (precision {}, iterations {})", "Reference", time.elapsed().as_millis(), self.center_location.prec().0, reference.current_iteration);

        let time = Instant::now();

        let indices = (0..self.image_width).cartesian_product(0..self.image_height);

        let mut pixel_data = indices.into_iter().par_bridge()
                            .map(|(i, j)| {
                                let element = ComplexFixed::new(i as f64 * delta_pixel + delta_top_left.re, j as f64 * delta_pixel + delta_top_left.im);
                                let point_delta = ComplexExtended::new(element, -self.zoom.exponent);
                                let new_delta = series_approximation.evaluate(point_delta);

                                PixelData {
                                    image_x: i,
                                    image_y: j,
                                    iteration: reference.start_iteration,
                                    delta_reference: point_delta,
                                    delta_current: new_delta,
                                    derivative_current: ComplexFixed::new(1.0, 0.0),
                                    glitched: false,
                                    escaped: false
                                }
                            }).collect::<Vec<PixelData>>();

        println!("{:<14}{:>6} ms", "Packing", time.elapsed().as_millis());

        let time = Instant::now();
        Perturbation::iterate(&mut pixel_data, &reference, reference.current_iteration);
        println!("{:<14}{:>6} ms", "Iteration", time.elapsed().as_millis());

        let time = Instant::now();
        ColouringExtended::Iteration.run(&pixel_data, &mut self.image, self.maximum_iteration, delta_pixel);
        println!("{:<14}{:>6} ms", "Coloring", time.elapsed().as_millis());

        // Remove all non-glitched points from the remaining points
        pixel_data.retain(|packet| {
            packet.glitched
        });

        while pixel_data.len() as f64 > 0.01 * self.glitch_tolerance * (self.image_width * self.image_height) as f64 {
            // delta_c is the difference from the next reference from the previous one
            let delta_c = pixel_data.choose(&mut rand::thread_rng()).unwrap().clone();
            let element = ComplexFixed::new(delta_c.image_x as f64 * delta_pixel + delta_top_left.re, delta_c.image_y as f64 * delta_pixel + delta_top_left.im);

            let reference_wrt_sa = ComplexExtended::new(element, -self.zoom.exponent);

            let delta_z = series_approximation.evaluate(reference_wrt_sa);

            let mut r = series_approximation.get_reference(reference_wrt_sa);
            r.run();

            // this can be made faster, without having to do the series approximation again
            // this is done by storing more data in pixeldata2
            pixel_data.chunks_mut(1)
                        .for_each(|pixel_data| {
                            for data in pixel_data {
                                let element = ComplexFixed::new(data.image_x as f64 * delta_pixel + delta_top_left.re, data.image_y as f64 * delta_pixel + delta_top_left.im);
                                let point_delta = ComplexExtended::new(element, -self.zoom.exponent);

                                data.iteration = reference.start_iteration;
                                data.glitched = false;
                                data.escaped = false;
                                data.delta_current = series_approximation.evaluate(point_delta) - delta_z;
                                data.delta_reference = point_delta - reference_wrt_sa;
                                // might not need the evaluate here as if we store it separately, there is no need
                                data.derivative_current = ComplexFixed::new(1.0, 0.0);
                            }
                        });

            Perturbation::iterate(&mut pixel_data, &r, r.current_iteration);

            ColouringExtended::Iteration.run(&pixel_data, &mut self.image, self.maximum_iteration, delta_pixel);

            // Remove all non-glitched points from the remaining points
            pixel_data.retain(|packet| {
                packet.glitched
            });
        }

        println!("{:<14}{:>6} ms (remaining {})", "Fixing", time.elapsed().as_millis(), pixel_data.len());
        
        let time = Instant::now();
        self.image.save();
        println!("{:<14}{:>6} ms", "Saving", time.elapsed().as_millis());
        return;
    }
}