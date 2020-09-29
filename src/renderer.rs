use crate::util::{data_export::*, ComplexFixed, ComplexArbitrary, PixelData, complex_extended::ComplexExtended, float_extended::FloatExtended, string_to_extended, extended_to_string_short, extended_to_string_long, get_approximation_terms, get_delta_top_left};
use crate::math::{SeriesApproximation, Perturbation, Reference};

use std::time::Instant;
use std::io::Write;
use std::cmp::max;

use rand::seq::SliceRandom;
use rayon::prelude::*;
use config::Config;

pub struct FractalRenderer {
    image_width: usize,
    image_height: usize,
    rotate: f64,
    zoom: FloatExtended,
    auto_adjust_iterations: bool,
    maximum_iteration: usize,
    glitch_tolerance: f64,
    data_export: DataExport,
    start_render_time: Instant,
    remaining_frames: usize,
    frame_offset: usize,
    zoom_scale_factor: f64,
    center_reference: Reference,
    series_approximation: SeriesApproximation,
    render_indices: Vec<usize>,
    remove_centre: bool
}

impl FractalRenderer {
    pub fn new(settings: Config) -> Self {
        let image_width = settings.get_int("image_width").unwrap_or(1000) as usize;
        let image_height = settings.get_int("image_height").unwrap_or(1000) as usize;
        let rotate = settings.get_float("rotate").unwrap_or(0.0).to_radians();
        let maximum_iteration = settings.get_int("iterations").unwrap_or(1000) as usize;
        let initial_zoom = settings.get_str("zoom").unwrap_or(String::from("1E0")).to_ascii_uppercase();
        let center_real = settings.get_str("real").unwrap_or(String::from("-0.75"));
        let center_imag = settings.get_str("imag").unwrap_or(String::from("0.0"));
        let approximation_order = settings.get_int("approximation_order").unwrap_or(0) as usize;
        let glitch_tolerance = settings.get_float("glitch_tolerance").unwrap_or(0.001);
        let remaining_frames = settings.get_int("frames").unwrap_or(1) as usize;
        let frame_offset = settings.get_int("frame_offset").unwrap_or(0) as usize;
        let zoom_scale_factor = settings.get_float("zoom_scale").unwrap_or(2.0);
        let display_glitches = settings.get_bool("display_glitches").unwrap_or(false);
        let auto_adjust_iterations = settings.get_bool("auto_adjust_iterations").unwrap_or(true);
        let experimental = settings.get_bool("experimental").unwrap_or(false);
        let probe_sampling = settings.get_int("probe_sampling").unwrap_or(3) as usize;
        let remove_centre = settings.get_bool("remove_centre").unwrap_or(true);
        let data_type = match settings.get_str("export").unwrap_or(String::from("COLOUR")).to_ascii_uppercase().as_ref() {
            "RAW" | "EXR" => DataType::RAW,
            "COLOUR" | "PNG" => DataType::COLOUR,
            "KFB" => DataType::KFB,
            "BOTH" => DataType::BOTH,
            _ => DataType::COLOUR
        };

        let mut zoom = string_to_extended(&initial_zoom);
        let delta_pixel =  (-2.0 * (4.0 / image_height as f64 - 2.0) / zoom) / image_height as f64;
        let radius = delta_pixel * image_width as f64;
        let precision = max(64, -radius.exponent + 64);
        let center_location = ComplexArbitrary::with_val(
            precision as u32,
            ComplexArbitrary::parse("(".to_owned() + &center_real + "," + &center_imag + ")").expect("Location is not valid!"));
        let auto_approximation = get_approximation_terms(approximation_order, image_width, image_height);

        let reference = Reference::new(center_location.clone(), center_location.clone(), 1, maximum_iteration);
        let series_approximation = SeriesApproximation::new_central(&center_location, 
            auto_approximation, 
            maximum_iteration, 
            FloatExtended::new(0.0, 0), 
            ComplexExtended::new2(0.0, 0.0, 0),
            probe_sampling,
            experimental);
        let render_indices = (0..(image_width * image_height)).collect::<Vec<usize>>();

        // Change the zoom level to the correct one for the frame offset
        for _ in 0..frame_offset {
            zoom.mantissa /= zoom_scale_factor;
            zoom.reduce();
        }

        FractalRenderer {
            image_width,
            image_height,
            rotate,
            zoom,
            auto_adjust_iterations,
            maximum_iteration,
            glitch_tolerance,
            data_export: DataExport::new(image_width, image_height, display_glitches, data_type),
            start_render_time: Instant::now(),
            remaining_frames,
            frame_offset,
            zoom_scale_factor,
            center_reference: reference,
            series_approximation,
            render_indices,
            remove_centre,
        }
    }

    pub fn render_frame(&mut self, frame_index: usize, filename: String) {
        print!("{:<6}", frame_index + self.frame_offset);
        print!("| {:<15}", extended_to_string_short(self.zoom));
        std::io::stdout().flush().unwrap();
        let frame_time = Instant::now();
        let approximation_time = Instant::now();

        if frame_index == 0 {
            self.center_reference.run();
            self.series_approximation.maximum_iteration = self.center_reference.current_iteration;
            self.series_approximation.generate_approximation(&self.center_reference);
        }
        
        let cos_rotate = self.rotate.cos();
        let sin_rotate = self.rotate.sin();
        let delta_pixel = 4.0 / ((self.image_height - 1) as f64 * self.zoom.mantissa);
        let delta_top_left = get_delta_top_left(delta_pixel, self.image_width, self.image_height, cos_rotate, sin_rotate);

        let delta_pixel_extended = FloatExtended::new(delta_pixel, -self.zoom.exponent);

        self.series_approximation.delta_pixel_square = delta_pixel_extended * delta_pixel_extended;

        // Used for placing the probe points
        self.series_approximation.delta_top_left = ComplexExtended::new(delta_top_left, -self.zoom.exponent);
        self.series_approximation.check_approximation();

        print!("| {:<15}", approximation_time.elapsed().as_millis());
        print!("| {:<15}", self.series_approximation.valid_iteration);
        print!("| {:<6}", self.series_approximation.order);
        print!("| {:<15}", self.maximum_iteration);
        std::io::stdout().flush().unwrap();

        let packing_time = Instant::now();

        if (frame_index + self.frame_offset) != 0 && self.remove_centre {
            // This will remove the central pixels
            self.data_export.clear_buffers();

            let image_width = self.image_width;
            let image_height = self.image_height;
            let temp = 0.5 - 0.5 / self.zoom_scale_factor;

            // Set up new render indices
            self.render_indices.retain(|index| {
                let i = index % image_width;
                let j = index / image_width;

                // Add one to avoid rescaling artifacts
                let val1 = (image_width as f64 * temp).ceil() as usize;
                let val2 = (image_height as f64 * temp).ceil() as usize;

                i <= val1 || i >= image_width - val1 || j <= val2 || j >= image_height - val2
            });

            // The centre has already been removed
            self.remove_centre = false;
        }

        let mut pixel_data = (&self.render_indices).into_par_iter()
            .map(|index| {
                let i = index % self.image_width;
                let j = index / self.image_width;

                // This could be changed to account for jittering if needed
                let element = ComplexFixed::new(
                    i as f64 * delta_pixel * cos_rotate - j as f64 * delta_pixel * sin_rotate + delta_top_left.re, 
                    i as f64 * delta_pixel * sin_rotate + j as f64 * delta_pixel * cos_rotate + delta_top_left.im
                );

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
        self.data_export.save(&filename, self.maximum_iteration, self.series_approximation.order, &extended_to_string_long(self.zoom));
        print!("| {:<15}", saving_time.elapsed().as_millis());
        println!("| {:<15}| {:<15}", frame_time.elapsed().as_millis(), self.start_render_time.elapsed().as_millis());
        std::io::stdout().flush().unwrap();
    }

    pub fn render(&mut self) {
        // Print out the status information
        println!("{:<6}| {:<15}| {:<15}| {:<15}| {:<6}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}| {:<15}", "Frame", "Zoom", "Approx [ms]", "Skipped [it]", "Order", "Maximum [it]", "Packing [ms]", "Iteration [ms]", "Correct [ms]", "Saving [ms]", "Frame [ms]", "TOTAL [ms]");

        let mut count = 0;

        while self.remaining_frames > 0 && self.zoom.to_float() > 0.5 {
            self.render_frame(count, 
                format!("output/{:08}_{}", count + self.frame_offset, 
                extended_to_string_short(self.zoom)));

            self.zoom.mantissa /= self.zoom_scale_factor;
            self.zoom.reduce();

            if self.zoom.to_float() < 1e10 {
                // SA has some problems with precision with lots of terms at lot zoom levels
                if self.series_approximation.order > 8 {
                    // Overwrite the series approximation order
                    self.series_approximation.order = 8;
                    self.series_approximation.maximum_iteration = self.center_reference.current_iteration;
                    self.series_approximation.generate_approximation(&self.center_reference);
                }

                // Logic in here to automatically adjust the maximum number of iterations
                // This is done arbitrarily and could be done in a config file if required
                if self.auto_adjust_iterations && self.maximum_iteration > 10000 {
                    let new_iteration_value = 10000;

                    if self.center_reference.current_iteration >= 10000 {
                        self.center_reference.current_iteration = new_iteration_value;
                    };

                    self.center_reference.maximum_iteration = new_iteration_value;
                    self.maximum_iteration = new_iteration_value;
                }
            }
            
            self.remaining_frames -= 1;
            count += 1;
        }
    }
}