use crate::util::{data_export::*, ComplexFixed, ComplexArbitrary, PixelData, complex_extended::ComplexExtended, float_extended::FloatExtended, string_to_extended, extended_to_string};
use crate::math::{SeriesApproximation, Perturbation, Reference};

use std::time::Instant;
use std::cmp::{min, max};
use std::io::Write;

use rand::seq::SliceRandom;
use rayon::prelude::*;
use config::Config;

pub struct FractalRenderer {
    image_width: usize,
    image_height: usize,
    aspect: f64,
    zoom: FloatExtended,
    maximum_iteration: usize,
    glitch_tolerance: f64,
    data_export: DataExport,
    start_render_time: Instant,
    remaining_frames: usize,
    zoom_scale_factor: f64,
    center_reference: Reference,
    series_approximation: SeriesApproximation
}

impl FractalRenderer {
    pub fn new(settings: Config) -> Self {
        let image_width = settings.get_int("image_width").unwrap_or(1000) as usize;
        let image_height = settings.get_int("image_height").unwrap_or(1000) as usize;
        let maximum_iteration = settings.get_int("iterations").unwrap_or(1000) as usize;
        let initial_zoom = settings.get_str("zoom").unwrap_or(String::from("1E0")).to_ascii_uppercase();
        let center_real = settings.get_str("real").unwrap_or(String::from("-0.75"));
        let center_imag = settings.get_str("imag").unwrap_or(String::from("0.0"));
        let approximation_order = settings.get_int("approximation_order").unwrap_or(0) as usize;
        let glitch_tolerance = settings.get_float("glitch_tolerance").unwrap_or(0.01);
        let remaining_frames = settings.get_int("frames").unwrap_or(1) as usize;
        let zoom_scale_factor = settings.get_float("zoom_scale").unwrap_or(2.0);
        let data_type = match settings.get_str("export").unwrap_or(String::from("COLOUR")).to_ascii_uppercase().as_ref() {
            "RAW" => DataType::RAW,
            "COLOUR" => DataType::COLOUR,
            "BOTH" => DataType::BOTH,
            _ => DataType::COLOUR
        };
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

        let reference = Reference::new(center_location.clone(), center_location.clone(), 1, maximum_iteration);
        let series_approximation = SeriesApproximation::new_central(&center_location, 
            auto_approximation, 
            maximum_iteration, 
            FloatExtended::new(0.0, 0), 
            ComplexExtended::new2(0.0, 0.0, 0));

        FractalRenderer {
            image_width,
            image_height,
            aspect,
            zoom,
            maximum_iteration,
            glitch_tolerance,
            data_export: DataExport::new(image_width, image_height, display_glitches, data_type),
            start_render_time: Instant::now(),
            remaining_frames,
            zoom_scale_factor,
            center_reference: reference,
            series_approximation,
        }
    }

    pub fn render_frame(&mut self, index: usize, filename: String) {
        print!("{:<6}", index);
        print!("| {:<15}", extended_to_string(self.zoom));
        std::io::stdout().flush().unwrap();
        let frame_time = Instant::now();
        let approximation_time = Instant::now();

        // If we are on the first frame we need to run the central reference calculation
        if index == 0 {
            self.center_reference.run();
            self.series_approximation.maximum_iteration = self.center_reference.current_iteration;
            self.series_approximation.generate_approximation(&self.center_reference);
        }
        
        let delta_pixel =  (-2.0 * (4.0 / self.image_height as f64 - 2.0) / self.zoom.mantissa) / self.image_height as f64;
        // This should be the delta relative to the image, without the big zoom factor applied.
        let delta_top_left = ComplexFixed::new((4.0 / self.image_width as f64 - 2.0) / self.zoom.mantissa * self.aspect as f64, (4.0 / self.image_height as f64 - 2.0) / self.zoom.mantissa);
        let delta_pixel_extended = FloatExtended::new(delta_pixel, -self.zoom.exponent);

        self.series_approximation.delta_pixel_square = delta_pixel_extended * delta_pixel_extended;
        self.series_approximation.delta_top_left = ComplexExtended::new(delta_top_left, -self.zoom.exponent);
        self.series_approximation.check_approximation();

        print!("| {:<15}", approximation_time.elapsed().as_millis());
        print!("| {:<15}", self.series_approximation.valid_iteration);
        std::io::stdout().flush().unwrap();

        let packing_time = Instant::now();

        let mut pixel_data = (0..(self.image_width * self.image_height)).into_par_iter()
            .map(|index| {
                let i = index % self.image_width;
                let j = index / self.image_width;
                let element = ComplexFixed::new(i as f64 * delta_pixel + delta_top_left.re, j as f64 * delta_pixel + delta_top_left.im);
                let point_delta = ComplexExtended::new(element, -self.zoom.exponent);

                let new_delta = self.series_approximation.evaluate(point_delta, None);

                PixelData {
                    image_x: i,
                    image_y: j,
                    iteration: self.series_approximation.valid_iteration,
                    delta_centre: point_delta,
                    delta_reference: point_delta,
                    delta_current: new_delta,
                    derivative_current: ComplexFixed::new(1.0, 0.0),
                    glitched: false,
                    escaped: false
                }
            }).collect::<Vec<PixelData>>();

        print!("| {:<15}", packing_time.elapsed().as_millis());
        std::io::stdout().flush().unwrap();

        let iteration_time = Instant::now();

        // This one has no offset because it is not a glitch resolving reference
        Perturbation::iterate(&mut pixel_data, &self.center_reference);
        self.data_export.export_pixels(&pixel_data, self.maximum_iteration, &self.center_reference);
        print!("| {:<15}", iteration_time.elapsed().as_millis());
        std::io::stdout().flush().unwrap();


        let correction_time = Instant::now();

        // Remove all non-glitched points from the remaining points
        pixel_data.retain(|packet| {
            packet.glitched
        });

        while pixel_data.len() as f64 > 0.01 * self.glitch_tolerance * (self.image_width * self.image_height) as f64 {
            // delta_c is the difference from the next reference from the previous one
            let glitch_reference_pixel = pixel_data.choose(&mut rand::thread_rng()).unwrap().clone();

            let mut glitch_reference = self.series_approximation.get_reference(glitch_reference_pixel.delta_centre, &self.center_reference);
            glitch_reference.run();

            let delta_current_reference = self.series_approximation.evaluate(glitch_reference_pixel.delta_centre, Some(glitch_reference.start_iteration));

            pixel_data.par_iter_mut()
                .for_each(|pixel| {
                    pixel.iteration = glitch_reference.start_iteration;
                    pixel.glitched = false;
                    pixel.delta_current = self.series_approximation.evaluate( pixel.delta_centre, Some(glitch_reference.start_iteration)) - delta_current_reference;
                    pixel.delta_reference = pixel.delta_centre - glitch_reference_pixel.delta_centre;
                });

            Perturbation::iterate(&mut pixel_data, &glitch_reference);
            self.data_export.export_pixels(&pixel_data, self.maximum_iteration, &glitch_reference);

            // Remove all non-glitched points from the remaining points
            pixel_data.retain(|packet| {
                packet.glitched
            });
        }

        print!("| {:<15}", correction_time.elapsed().as_millis());
        std::io::stdout().flush().unwrap();

        
        let saving_time = Instant::now();
        self.data_export.save(&filename);
        print!("| {:<15}", saving_time.elapsed().as_millis());
        println!("| {:<15}| {:<15}", frame_time.elapsed().as_millis(), self.start_render_time.elapsed().as_millis());
        std::io::stdout().flush().unwrap();
    }

    pub fn render(&mut self) {
        // Print out the status information
        println!("{:<6}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}", "Frame", "Zoom", "Approx [ms]", "Skipped [it]", "Packing [ms]", "Iteration [ms]", "Correct [ms]", "Saving [ms]", "Frame [ms]", "TOTAL [ms]");

        let mut count = 0;
        while self.remaining_frames > 0 && self.zoom.to_float() > 0.5 {
            self.render_frame(count, format!("output/output_{:08}", count));
            self.zoom.mantissa /= self.zoom_scale_factor;
            self.zoom.reduce();
            self.remaining_frames -= 1;
            count += 1;
        }
    }
}