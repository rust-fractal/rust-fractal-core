use crate::util::{data_export::*, ComplexFixed, ComplexArbitrary, PixelData, complex_extended::ComplexExtended, float_extended::FloatExtended, string_to_extended, extended_to_string};
use crate::math::{SeriesApproximation, Perturbation, Reference};

use std::time::Instant;
use std::cmp::{min, max};

use rand::seq::SliceRandom;
use rayon::prelude::*;
use config::Config;

pub struct FractalRenderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: FloatExtended,
    center_location: ComplexArbitrary,
    maximum_iteration: usize,
    approximation_order: usize,
    glitch_tolerance: f64,
    data_export: DataExport,
    start_render_time: Instant,
    remaining_frames: usize,
    zoom_scale_factor: f64,
    center_reference: Reference,
}

impl FractalRenderer {
    pub fn new(settings: Config) -> Self {
        // Print out the status information
        println!("{:<6}| {:<14}| {:<14}| {:<14}| {:<14}| {:<14}| {:<14}| {:<14}| {:<14}| {:<14}| {:<14}", "Frame", "Zoom", "Approx [ms]", "Skipped [it]", "Reference [ms]", "Packing [ms]", "Iteration [ms]", "Colouring [ms]", "Correct [ms]", "Saving [ms]", "TOTAL [ms]");

        let image_width = settings.get_int("image_width").unwrap_or(1000) as usize;
        let image_height = settings.get_int("image_height").unwrap_or(1000) as usize;
        let maximum_iteration = settings.get_int("iterations").unwrap_or(1000) as usize;
        let initial_zoom = settings.get_str("zoom").unwrap_or("1E0".to_string());
        let center_real = settings.get_str("real").unwrap_or("-0.75".to_string());
        let center_imag = settings.get_str("imag").unwrap_or("0.0".to_string());
        let approximation_order = settings.get_int("approximation_order").unwrap_or(0) as usize;
        let glitch_tolerance = settings.get_float("glitch_tolerance").unwrap_or(0.01);
        let remaining_frames = settings.get_int("frames").unwrap_or(1) as usize;
        let zoom_scale_factor = settings.get_float("zoom_scale").unwrap_or(2.0);
        let display_glitches = false;

        let aspect = image_width as f64 / image_height as f64;
        let zoom = string_to_extended(&initial_zoom);

        let delta_pixel =  (-2.0 * (4.0 / image_height as f64 - 2.0) / zoom) / image_height as f64;

        let radius = delta_pixel * image_width as f64;
        let precision = max(64, -radius.exponent + 64);

        let center_location = ComplexArbitrary::with_val(
            precision as u32,
            ComplexArbitrary::parse("(".to_owned() + &center_real + "," + &center_imag + ")").expect("Location is not valid!"));

        let auto_approximation = if approximation_order == 0 {
            let auto = (((image_width * image_height) as f64).log(1e6).powf(6.619) * 16.0f64) as usize;
            min(max(auto, 3), 64)
        } else {
            approximation_order
        };

        FractalRenderer {
            image_width,
            image_height,
            aspect,
            zoom,
            center_location: center_location.clone(),
            maximum_iteration,
            approximation_order: auto_approximation,
            glitch_tolerance,
            data_export: DataExport::new(image_width, image_height, display_glitches, DataType::BOTH),
            start_render_time: Instant::now(),
            remaining_frames,
            zoom_scale_factor,
            center_reference: Reference::new(center_location.clone(), center_location, 1, maximum_iteration)
        }
    }

    pub fn render_frame(&mut self, index: usize, filename: String) {
        print!("{:<6}", index);
        print!("| {:<14}", extended_to_string(self.zoom));
        let frame_time = Instant::now();
        let time = Instant::now();
        let delta_pixel =  (-2.0 * (4.0 / self.image_height as f64 - 2.0) / self.zoom.mantissa) / self.image_height as f64;

        // this should be the delta relative to the image, without the big zoom factor applied.
        let delta_top_left = ComplexFixed::new((4.0 / self.image_width as f64 - 2.0) / self.zoom.mantissa * self.aspect as f64, (4.0 / self.image_height as f64 - 2.0) / self.zoom.mantissa);

        let delta_pixel_extended = FloatExtended::new(delta_pixel, -self.zoom.exponent);

        // Series approximation currently has some overskipping issues
        // this can be resolved by root finding and adding new probe points
        let mut series_approximation = SeriesApproximation::new(
            self.center_location.clone(),
            self.approximation_order,
            self.maximum_iteration,
            delta_pixel_extended * delta_pixel_extended,
            ComplexExtended::new(delta_top_left, -self.zoom.exponent),
        );

        series_approximation.run();

        print!("| {:<14}", time.elapsed().as_millis());
        print!("| {:<14}", series_approximation.current_iteration);

        let time = Instant::now();

        let mut reference = series_approximation.get_reference(ComplexExtended::new2(0.0, 0.0, 0));
        reference.run();

        print!("| {:<14}", time.elapsed().as_millis());

        let time = Instant::now();

        let mut pixel_data = (0..(self.image_width * self.image_height)).into_par_iter()
            .map(|index| {
                let i = index % self.image_width;
                let j = index / self.image_width;
                let element = ComplexFixed::new(i as f64 * delta_pixel + delta_top_left.re, j as f64 * delta_pixel + delta_top_left.im);
                let point_delta = ComplexExtended::new(element, -self.zoom.exponent);
                let new_delta = series_approximation.evaluate(point_delta);

                PixelData {
                    image_x: i,
                    image_y: j,
                    iteration: reference.start_iteration,
                    delta_centre: point_delta,
                    delta_reference: point_delta,
                    delta_start: new_delta,
                    delta_current: new_delta,
                    derivative_current: ComplexFixed::new(1.0, 0.0),
                    glitched: false,
                    escaped: false
                }
            }).collect::<Vec<PixelData>>();

        print!("| {:<14}", time.elapsed().as_millis());

        let time = Instant::now();
        Perturbation::iterate(&mut pixel_data, &reference, reference.current_iteration);
        print!("| {:<14}", time.elapsed().as_millis());

        let time = Instant::now();
        self.data_export.export_pixels(&pixel_data, self.maximum_iteration, &reference);
        print!("| {:<14}", time.elapsed().as_millis());

        let time = Instant::now();

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
            pixel_data.par_iter_mut()
                .for_each(|pixel| {
                    pixel.iteration = reference.start_iteration;
                    pixel.glitched = false;
                    pixel.delta_current = pixel.delta_start - delta_z;
                    pixel.delta_reference = pixel.delta_centre - reference_wrt_sa;
                        // might not need the evaluate here as if we store it separately, there is no need
                        // data.derivative_current = ComplexFixed::new(1.0, 0.0);
                });

            Perturbation::iterate(&mut pixel_data, &r, r.current_iteration);

            self.data_export.export_pixels(&pixel_data, self.maximum_iteration, &r);

            // Remove all non-glitched points from the remaining points
            pixel_data.retain(|packet| {
                packet.glitched
            });
        }

        print!("| {:<14}", time.elapsed().as_millis());
        
        let time = Instant::now();
        self.data_export.save(&filename);
        print!("| {:<14}", time.elapsed().as_millis());
        println!("| {:<8} {:<8}", frame_time.elapsed().as_millis(), self.start_render_time.elapsed().as_millis());
    }

    pub fn render(&mut self) {
        let mut count = 0;
        while self.remaining_frames > 0 && self.zoom.to_float() > 0.5 {
            self.render_frame(count, format!("output/keyframe_{:08}", count));
            self.zoom.mantissa /= self.zoom_scale_factor;
            self.zoom.reduce();
            self.remaining_frames -= 1;
            count += 1;
        }
    }
}