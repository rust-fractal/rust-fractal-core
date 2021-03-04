use crate::util::{data_export::*, ComplexFixed, ComplexArbitrary, PixelData, complex_extended::ComplexExtended, float_extended::FloatExtended, string_to_extended, extended_to_string_short, extended_to_string_long, get_approximation_terms, get_delta_top_left, generate_default_palette, ProgressCounters};
use crate::math::{SeriesApproximation, Perturbation, Reference};

use std::{cmp::Ordering, time::{Duration, Instant}};
use std::io::Write;
use std::cmp::{min, max};

use rand::{seq::SliceRandom, thread_rng};
use rand::Rng;

use rayon::prelude::*;
use config::Config;

use std::thread;
use std::sync::{Arc, mpsc};
use parking_lot::Mutex;

use atomic_counter::{AtomicCounter, RelaxedCounter};

pub struct FractalRenderer {
    pub image_width: usize,
    pub image_height: usize,
    pub total_pixels: usize,
    pub rotate: f64,
    pub zoom: FloatExtended,
    pub auto_adjust_iterations: bool,
    pub maximum_iteration: usize,
    pub glitch_percentage: f64,
    pub data_export: Arc<Mutex<DataExport>>,
    start_render_time: Instant,
    pub remaining_frames: usize,
    frame_offset: usize,
    pub zoom_scale_factor: f64,
    pub center_reference: Reference,
    pub series_approximation: SeriesApproximation,
    render_indices: Vec<usize>,
    pub remove_centre: bool,
    pub analytic_derivative: bool,
    pub jitter: bool,
    pub jitter_factor: f64,
    pub experimental: bool,
    show_output: bool,
    pub progress: ProgressCounters,
    pub render_time: u128
}

impl FractalRenderer {
    pub fn new(settings: Config) -> Self {
        let image_width = settings.get_int("image_width").unwrap_or(1000) as usize;
        let image_height = settings.get_int("image_height").unwrap_or(1000) as usize;
        let rotate = settings.get_float("rotate").unwrap_or(0.0).to_radians();
        let maximum_iteration = settings.get_int("iterations").unwrap_or(1000) as usize;
        let initial_zoom = settings.get_str("zoom").unwrap_or_else(|_| String::from("1E0")).to_ascii_uppercase();
        let center_real = settings.get_str("real").unwrap_or_else(|_| String::from("-0.75"));
        let center_imag = settings.get_str("imag").unwrap_or_else(|_| String::from("0.0"));
        let approximation_order = settings.get_int("approximation_order").unwrap_or(0) as usize;
        let glitch_percentage = settings.get_float("glitch_percentage").unwrap_or(0.001);
        let remaining_frames = settings.get_int("frames").unwrap_or(1) as usize;
        let frame_offset = settings.get_int("frame_offset").unwrap_or(0) as usize;
        let zoom_scale_factor = settings.get_float("zoom_scale").unwrap_or(2.0);
        let display_glitches = settings.get_bool("display_glitches").unwrap_or(false);
        let auto_adjust_iterations = settings.get_bool("auto_adjust_iterations").unwrap_or(true);
        let experimental = settings.get_bool("experimental").unwrap_or(false);
        let probe_sampling = settings.get_int("probe_sampling").unwrap_or(3) as usize;
        let remove_centre = settings.get_bool("remove_centre").unwrap_or(true);
        let iteration_division = settings.get_float("iteration_division").unwrap_or(0.1) as f32;
        let palette_offset = settings.get_float("palette_offset").unwrap_or(0.0) as f32;
        let valid_iteration_frame_multiplier = settings.get_float("valid_iteration_frame_multiplier").unwrap_or(0.25) as f32;
        let valid_iteration_probe_multiplier = settings.get_float("valid_iteration_probe_multiplier").unwrap_or(0.02) as f32;
        let glitch_tolerance = settings.get_float("glitch_tolerance").unwrap_or(1.4e-6) as f64;
        let data_storage_interval = settings.get_int("data_storage_interval").unwrap_or(10) as usize;
        let analytic_derivative = settings.get_bool("analytic_derivative").unwrap_or(false);
        let jitter = settings.get_bool("jitter").unwrap_or(false);
        let jitter_factor = settings.get_float("jitter_factor").unwrap_or(0.2);
        let show_output = settings.get_bool("show_output").unwrap_or(true);
        
        let data_type = match settings.get_str("export").unwrap_or_else(|_| String::from("COLOUR")).to_ascii_uppercase().as_ref() {
            "NONE" => DataType::None,
            "GUI" => DataType::Gui,
            "RAW" | "EXR" => DataType::Raw,
            "COLOUR" | "COLOR" | "PNG" => DataType::Color,
            "KFB" => DataType::Kfb,
            "BOTH" => DataType::Both,
            _ => DataType::Color
        };

        let palette = match data_type {
            DataType::Raw | DataType::None => {
                Vec::new()
            },
            DataType::Kfb | DataType::Color | DataType::Both | DataType::Gui => {
                if let Ok(colour_values) = settings.get_array("palette") {
                    colour_values.chunks_exact(3).map(|value| {
                        // We assume the palette is in BGR rather than RGB
                        (value[2].clone().into_int().unwrap() as u8, 
                            value[1].clone().into_int().unwrap() as u8, 
                            value[0].clone().into_int().unwrap() as u8)
                    }).collect::<Vec<(u8, u8, u8)>>()
                } else {
                    generate_default_palette()
                }
            }
        };

        let mut zoom = string_to_extended(&initial_zoom);
        let delta_pixel =  (-2.0 * (4.0 / image_height as f64 - 2.0) / zoom) / image_height as f64;
        let radius = delta_pixel * image_width as f64;
        let precision = max(64, -radius.exponent + 64);

        let center_location = ComplexArbitrary::with_val(
            precision as u32,
            ComplexArbitrary::parse("(".to_owned() + &center_real + "," + &center_imag + ")").expect("provided location not valid"));
        let auto_approximation = get_approximation_terms(approximation_order, image_width, image_height);

        let reference = Reference::new(center_location.clone(), 
            center_location, 
            1, 
            maximum_iteration, 
            data_storage_interval,
            glitch_tolerance,
            zoom);

        let series_approximation = SeriesApproximation::new_central(auto_approximation, 
            maximum_iteration, 
            FloatExtended::new(0.0, 0), 
            probe_sampling,
            experimental,
            valid_iteration_frame_multiplier,
            valid_iteration_probe_multiplier,
            data_storage_interval);

        let render_indices = FractalRenderer::generate_render_indices(image_width, image_height);

        // Change the zoom level to the correct one for the frame offset
        for _ in 0..frame_offset {
            zoom.mantissa /= zoom_scale_factor;
            zoom.reduce();
        }

        FractalRenderer {
            image_width,
            image_height,
            total_pixels: image_width * image_height,
            rotate,
            zoom,
            auto_adjust_iterations,
            maximum_iteration,
            glitch_percentage,
            data_export: Arc::new(Mutex::new(DataExport::new(image_width, image_height, display_glitches, data_type, palette, iteration_division, palette_offset, analytic_derivative))),
            start_render_time: Instant::now(),
            remaining_frames,
            frame_offset,
            zoom_scale_factor,
            center_reference: reference,
            series_approximation,
            render_indices,
            remove_centre,
            analytic_derivative,
            jitter,
            jitter_factor,
            experimental,
            show_output,
            progress: ProgressCounters::new(maximum_iteration),
            render_time: 0
        }
    }

    pub fn render_frame(&mut self, frame_index: usize, filename: String, stop_flag: Option<Arc<RelaxedCounter>>) {
        if self.show_output {
            print!("{:<6}", frame_index + self.frame_offset);
            print!("| {:<15}", extended_to_string_short(self.zoom));
            std::io::stdout().flush().unwrap();
        };

        let frame_time = Instant::now();
        let approximation_time = Instant::now();

        let stop_flag = match stop_flag {
            Some(stop_flag) => stop_flag,
            None => Arc::new(RelaxedCounter::new(0))
        };

        let (tx, rx) = mpsc::channel();

        if self.show_output {
            let thread_counter_1 = self.progress.reference.clone();
            let thread_counter_2 = self.progress.series_approximation.clone();
            let thread_counter_3 = self.progress.reference_maximum.clone();
            let thread_counter_4 = self.progress.series_validation.clone();
            print!("|               ");

            thread::spawn(move || {
                loop {
                    match rx.try_recv() {
                        Ok(_) => {
                            break;
                        },
                        Err(_) => {
                            // 40% weighting to first reference, 40% to SA calculation, 20% to SA checking
                            let mut percentage_complete = 0.0;

                            // println!("{} {} {}", thread_counter_1.get(), thread_counter_2.get(), thread_counter_3.get());

                            percentage_complete += 45.0 * thread_counter_1.get() as f64 / thread_counter_3.get() as f64;
                            percentage_complete += 45.0 * thread_counter_2.get() as f64 / thread_counter_3.get() as f64;
                            percentage_complete += 10.0 * thread_counter_4.get() as f64 / 2.0;

                            print!("\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08{:^14}", format!("{:.2}%", percentage_complete));
                            std::io::stdout().flush().unwrap();
                        }
                    };

                    thread::sleep(Duration::from_millis(100));
                };
            });
        };

        if frame_index == 0 {
            self.data_export.lock().maximum_iteration = self.maximum_iteration;

            self.center_reference.run(&self.progress.reference, &self.progress.reference_maximum, &stop_flag);

            if stop_flag.get() >= 1 {
                self.progress.reset();
                self.render_time = frame_time.elapsed().as_millis();
                return;
            };

            self.series_approximation.maximum_iteration = self.center_reference.current_iteration;
            self.series_approximation.generate_approximation(&self.center_reference, &self.progress.series_approximation, &stop_flag);
        } else {
            let mut export = self.data_export.lock();

            // If the image width/height changes intraframe (GUI) we need to regenerate some things
            if export.image_width != self.image_width || export.image_height != self.image_height {
                self.render_indices = FractalRenderer::generate_render_indices(self.image_width, self.image_height);

                export.image_width = self.image_width;
                export.image_height = self.image_height;
            }
            
            export.clear_buffers();

            drop(export);

            // Check to see if the series approximation order has changed intraframe
            if self.series_approximation.order != (self.series_approximation.coefficients[0].len() - 1) {
                // TODO make it so that the value is set back to zero, rather than remade
                // self.progress.reset_series_approximation();
                self.series_approximation.min_valid_iteration = 1;

                self.series_approximation.generate_approximation(&self.center_reference, &self.progress.series_approximation, &stop_flag);
            }
        }

        if stop_flag.get() >= 1 {
            self.progress.reset();
            self.render_time = frame_time.elapsed().as_millis();
            return;
        };
        
        let cos_rotate = self.rotate.cos();
        let sin_rotate = self.rotate.sin();
        let delta_pixel = 4.0 / ((self.image_height - 1) as f64 * self.zoom.mantissa);
        let delta_top_left = get_delta_top_left(delta_pixel, self.image_width, self.image_height, cos_rotate, sin_rotate);
        let delta_pixel_extended = FloatExtended::new(delta_pixel, -self.zoom.exponent);

        let minimum_dimension = min(self.image_width, self.image_height);

        self.series_approximation.delta_pixel_square = if minimum_dimension < 1000 {
            let fixed_delta_pixel_extended = FloatExtended::new(4.0 / (999.0 * self.zoom.mantissa), -self.zoom.exponent);
            
            fixed_delta_pixel_extended * fixed_delta_pixel_extended
        } else {
            delta_pixel_extended * delta_pixel_extended
        };

        // Used for placing the probe points
        self.series_approximation.check_approximation(
            delta_top_left, 
            -self.zoom.exponent, 
            cos_rotate, 
            sin_rotate, 
            delta_pixel,
            self.image_width,
            self.image_height,
            &self.center_reference,
            &self.progress.series_validation);

        // -1 because we already know that 1 iteration can be skipped (take z = c)
        self.progress.min_series_approximation.add(self.series_approximation.min_valid_iteration - 1);
        self.progress.max_series_approximation.add(self.series_approximation.max_valid_iteration - 1);

        if stop_flag.get() >= 1 {
            self.progress.reset();
            self.render_time = frame_time.elapsed().as_millis();
            return;
        };

        tx.send(()).unwrap();

        if self.show_output {
            print!("\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08{:<15}", approximation_time.elapsed().as_millis());
            print!("| {:<15}", self.series_approximation.min_valid_iteration);
            print!("| {:<6}", self.series_approximation.order);
            print!("| {:<15}", self.maximum_iteration);
            std::io::stdout().flush().unwrap();
        };

        let packing_time = Instant::now();

        if !self.remove_centre && self.data_export.lock().centre_removed {
            self.render_indices = FractalRenderer::generate_render_indices(self.image_width, self.image_height);

            self.data_export.lock().centre_removed = false;
        }

        // If the remove_centre flag is set, and either it is not the first frame or gui mode is enabled
        if (self.remove_centre && ((frame_index + self.frame_offset) != 0 || self.data_export.lock().data_type == DataType::Gui)) && !self.data_export.lock().centre_removed {
            // This will remove the central pixels
            let image_width = self.image_width;
            let image_height = self.image_height;
            let temp = 0.5 - 0.5 / self.zoom_scale_factor;

            // Add one to avoid rescaling artifacts
            let val1 = (image_width as f64 * temp).ceil() as usize;
            let val2 = (image_height as f64 * temp).ceil() as usize;

            // Set up new render indices
            self.render_indices.retain(|index| {
                let i = index % image_width;
                let j = index / image_width;

                i <= val1 || i >= image_width - val1 || j <= val2 || j >= image_height - val2
            });

            // The centre has already been removed
            self.data_export.lock().centre_removed = true;
            self.total_pixels = self.render_indices.len();
        }

        let mut pixel_data = (&self.render_indices).into_par_iter()
            .map(|index| {
                let image_x = index % self.image_width;
                let image_y = index / self.image_width;

                let mut i = image_x as f64;
                let mut j = image_y as f64;

                let chosen_iteration = if self.experimental {
                    let test1 = ((self.series_approximation.probe_sampling - 1) as f64 * i / self.image_width as f64).floor() as usize;
                    let test2 = ((self.series_approximation.probe_sampling - 1) as f64 * j / self.image_height as f64).floor() as usize;

                    let index = test2 * (self.series_approximation.probe_sampling - 1) + test1;

                    self.series_approximation.valid_interpolation[index]
                } else {
                    self.series_approximation.min_valid_iteration
                };

                if self.jitter {
                    let mut rng = rand::thread_rng();

                    i += rng.gen_range(-self.jitter_factor, self.jitter_factor);
                    j += rng.gen_range(-self.jitter_factor, self.jitter_factor);
                }

                // This could be changed to account for jittering if needed
                let element = ComplexFixed::new(
                    i * delta_pixel * cos_rotate - j * delta_pixel * sin_rotate + delta_top_left.re, 
                    i * delta_pixel * sin_rotate + j * delta_pixel * cos_rotate + delta_top_left.im
                );

                let point_delta = ComplexExtended::new(element, -self.zoom.exponent);
                let new_delta = self.series_approximation.evaluate(point_delta, chosen_iteration);

                let derivative = if self.analytic_derivative {
                    self.series_approximation.evaluate_derivative(point_delta, chosen_iteration)
                } else {
                    ComplexExtended::new2(1.0, 0.0, 0)
                };

                PixelData {
                    image_x,
                    image_y,
                    iteration: chosen_iteration,
                    delta_centre: point_delta,
                    delta_reference: point_delta,
                    delta_current: new_delta,
                    derivative_current: derivative,
                    glitched: false,
                    escaped: false,
                }
            }).collect::<Vec<PixelData>>();

        if self.show_output {
            print!("| {:<15}", packing_time.elapsed().as_millis());
            std::io::stdout().flush().unwrap();
        };

        if stop_flag.get() >= 1 {
            self.progress.reset();
            self.render_time = frame_time.elapsed().as_millis();
            return;
        };
        
        let iteration_time = Instant::now();

        let total_pixels = self.render_indices.len() as f64;

        let (tx, rx) = mpsc::channel();

        if self.show_output {
            let thread_counter = Arc::clone(&self.progress.iteration);
            print!("|               ");

            thread::spawn(move || {
                loop {
                    thread::sleep(Duration::from_millis(100));
                    match rx.try_recv() {
                        Ok(_) => {
                            break;
                        },
                        Err(_) => {
                            print!("\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08{:^14}", format!("{:.2}%", 100.0 * thread_counter.get() as f64 / total_pixels));
                            std::io::stdout().flush().unwrap();
                        }
                    };
                };
            });
        };

        let values = [16usize, 8, 4, 2, 1];
        let mut previous_value = 0;
        let number_pixels = self.render_indices.len();

        for value in values.iter() {
            // println!("{}", value);
            let end_value = number_pixels / (value * value);
            let chunk_size = max((end_value - previous_value) / 2048, 4);
            // println!("chunk size init: {}", chunk_size);

            // This one has no offset because it is not a glitch resolving reference
            if self.analytic_derivative {
                Perturbation::iterate_normal_plus_derivative(&mut pixel_data[previous_value..end_value], &self.center_reference, &self.progress.iteration, &stop_flag, self.data_export.clone(), delta_pixel_extended, *value, chunk_size);
            } else {
                Perturbation::iterate_normal(&mut pixel_data[previous_value..end_value], &self.center_reference, &self.progress.iteration, &stop_flag, self.data_export.clone(), delta_pixel_extended, *value, chunk_size);
            };

            previous_value = end_value;
        }
        

        tx.send(()).unwrap();

        if stop_flag.get() >= 1 {
            self.progress.reset();
            self.render_time = frame_time.elapsed().as_millis();
            return;
        };

        // self.data_export.export_pixels(&pixel_data, &self.center_reference, delta_pixel_extended);

        if self.show_output {
            print!("\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08{:<15}", iteration_time.elapsed().as_millis());
            std::io::stdout().flush().unwrap();
        };
        
        let correction_time = Instant::now();
        let mut correction_references = 1;

        // Remove all non-glitched points from the remaining points
        pixel_data.retain(|packet| {
            packet.glitched
        });

        let glitched_pixels = pixel_data.len() as f64;
        self.progress.glitched_maximum.add(pixel_data.len());

        let complete_pixels = total_pixels - glitched_pixels;

        let (tx, rx) = mpsc::channel();

        if self.show_output {
            let thread_counter = Arc::clone(&self.progress.iteration);
            print!("|               ");

            thread::spawn(move || {
                loop {
                    thread::sleep(Duration::from_millis(100));
                    match rx.try_recv() {
                        Ok(_) => {
                            break;
                        },
                        Err(_) => {
                            print!("\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08{:^14}", format!("{:.2}%", 100.0 * (thread_counter.get() as f64 - complete_pixels) / glitched_pixels));
                            std::io::stdout().flush().unwrap();
                        }
                    };
                };
            });
        };

        while pixel_data.len() as f64 > 0.01 * self.glitch_percentage * self.total_pixels as f64 {
            // delta_c is the difference from the next reference from the previous one

            let glitch_reference_pixel = pixel_data.choose(&mut rand::thread_rng()).unwrap().clone();
            // let glitch_reference_pixel = pixel_data[0].clone();

            let mut glitch_reference = self.series_approximation.get_reference(glitch_reference_pixel.delta_centre, &self.center_reference);

            correction_references += 1;
            glitch_reference.run(&Arc::new(RelaxedCounter::new(0)), &Arc::new(RelaxedCounter::new(0)), &stop_flag);

            if stop_flag.get() >= 1 {
                self.render_time = frame_time.elapsed().as_millis();
                self.progress.reset();
                return;
            };

            let delta_current_reference = self.series_approximation.evaluate(glitch_reference_pixel.delta_centre, glitch_reference.start_iteration);

            if self.analytic_derivative {
                pixel_data.par_iter_mut()
                    .for_each(|pixel| {
                        pixel.iteration = glitch_reference.start_iteration;
                        pixel.glitched = false;
                        pixel.delta_current = self.series_approximation.evaluate( pixel.delta_centre, glitch_reference.start_iteration) - delta_current_reference;
                        pixel.delta_reference = pixel.delta_centre - glitch_reference_pixel.delta_centre;
                        pixel.derivative_current = self.series_approximation.evaluate_derivative(pixel.delta_centre, glitch_reference.start_iteration);
                });
            } else {
                pixel_data.par_iter_mut()
                    .for_each(|pixel| {
                        pixel.iteration = glitch_reference.start_iteration;
                        pixel.glitched = false;
                        pixel.delta_current = self.series_approximation.evaluate( pixel.delta_centre, glitch_reference.start_iteration) - delta_current_reference;
                        pixel.delta_reference = pixel.delta_centre - glitch_reference_pixel.delta_centre;
                });
            }
            
            let chunk_size = max(pixel_data.len() / 2048, 4);
            // println!("chunk size: {}", chunk_size);

            if self.analytic_derivative {
                Perturbation::iterate_normal_plus_derivative(&mut pixel_data, &glitch_reference, &self.progress.iteration, &stop_flag, self.data_export.clone(), delta_pixel_extended, 1, chunk_size);
            } else {
                Perturbation::iterate_normal(&mut pixel_data, &glitch_reference, &self.progress.iteration, &stop_flag, self.data_export.clone(), delta_pixel_extended, 1, chunk_size);
            };

            // self.data_export.export_pixels(&pixel_data, &glitch_reference, delta_pixel_extended);

            // Remove all non-glitched points from the remaining points
            // TODO maybe there is a function to split into two vectors, we only need to export the non glitched pixels
            pixel_data.retain(|packet| {
                packet.glitched
            });
        };

        tx.send(()).unwrap();

        // Possibly here correct the small glitches
        // self.data_export.interpolate_glitches(&pixel_data);
        // self.data_export.export_pixels(&pixel_data, &glitch_reference, delta_pixel_extended);

        if self.show_output {
            print!("\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08\x08{:<15}", correction_time.elapsed().as_millis());
            print!("| {:<6}", correction_references);
            std::io::stdout().flush().unwrap();
        };
        
        let saving_time = Instant::now();
        self.data_export.lock().save(&filename, self.series_approximation.order, &extended_to_string_long(self.zoom));

        self.render_time = frame_time.elapsed().as_millis();

        self.progress.reset();

        if self.show_output {
            print!("| {:<15}", saving_time.elapsed().as_millis());
            println!("| {:<15}| {:<15}", frame_time.elapsed().as_millis(), self.start_render_time.elapsed().as_millis());
            std::io::stdout().flush().unwrap();
        }
    }

    // Returns true if the maximum iterations has been increased
    pub fn adjust_iterations(&mut self) -> bool {
        if self.auto_adjust_iterations {
            if 4 * self.series_approximation.max_valid_iteration > self.maximum_iteration {
                self.maximum_iteration *= 3;
                self.maximum_iteration /= 2;
                return true;
            }

            if 8 * self.series_approximation.max_valid_iteration < self.maximum_iteration {
                self.maximum_iteration *= 3;
                self.maximum_iteration /= 4;

                if self.maximum_iteration < 1000 {
                    self.maximum_iteration = 1000;
                }

                self.data_export.lock().maximum_iteration = self.maximum_iteration;
                
                if self.center_reference.current_iteration > self.maximum_iteration {
                    self.center_reference.current_iteration = self.maximum_iteration;
                }

                self.center_reference.maximum_iteration = self.maximum_iteration;
            }
        }
        false
    }

    pub fn render(&mut self) {
        // Print out the status information
        if self.show_output {
            println!("{:<6}| {:<15}| {:<15}| {:<15}| {:<6}| {:<15}| {:<15}| {:<15}| {:<15}| {:<6}| {:<15}| {:<15}| {:<15}", "Frame", "Zoom", "Approx [ms]", "Skipped [it]", "Order", "Maximum [it]", "Packing [ms]", "Iteration [ms]", "Correct [ms]", "Ref", "Saving [ms]", "Frame [ms]", "TOTAL [ms]");
        };

        let mut count = 0;

        while self.remaining_frames > 0 && self.zoom.to_float() > 0.5 {
            self.render_frame(count, 
                format!("output/{:08}_{}", count + self.frame_offset, 
                extended_to_string_short(self.zoom)), None);

            self.zoom.mantissa /= self.zoom_scale_factor;
            self.zoom.reduce();

            if self.zoom.to_float() < 1e10 {
                // Set these to start from the beginning
                self.series_approximation.valid_iteration_frame_multiplier = 1.0;
                self.series_approximation.valid_iteration_probe_multiplier = 1.0;

                // SA has some problems with precision with lots of terms at lot zoom levels
                if self.series_approximation.order > 8 {
                    // Overwrite the series approximation order
                    self.series_approximation.order = 8;
                    self.series_approximation.maximum_iteration = self.center_reference.current_iteration;
                    self.series_approximation.generate_approximation(&self.center_reference, &self.progress.series_approximation, &Arc::new(RelaxedCounter::new(0)));
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
            } else if self.series_approximation.min_valid_iteration < 1000 && self.series_approximation.order > 16 {
                    self.series_approximation.order = 16;
                    self.series_approximation.generate_approximation(&self.center_reference, &self.progress.series_approximation, &Arc::new(RelaxedCounter::new(0)));
            } else if self.series_approximation.min_valid_iteration < 10000 && self.series_approximation.order > 32 {
                self.series_approximation.order = 32;
                self.series_approximation.generate_approximation(&self.center_reference, &self.progress.series_approximation, &Arc::new(RelaxedCounter::new(0)));
            }
            
            self.remaining_frames -= 1;
            count += 1;
        }
    }

    pub fn generate_render_indices(image_width: usize, image_height: usize) -> Vec<usize> {
        // let time = Instant::now();

        let chooser = 3;
        let mut indices = Vec::with_capacity(image_width * image_height);

        match chooser {
            1 => {
                indices = (0..(image_width * image_height)).collect::<Vec<usize>>();
                indices.shuffle(&mut thread_rng());
            }
            2 => {
                indices = (0..(image_width * image_height)).collect::<Vec<usize>>();
                indices.sort_by(|a, b| {
                    let a_i = a % image_width;
                    let a_j = a / image_width;

                    let b_i = b % image_width;
                    let b_j = b / image_width;

                    let a_dist = (a_i - image_width / 2).pow(2) + (a_j - image_height / 2).pow(2);
                    let b_dist = (b_i - image_width / 2).pow(2) + (b_j - image_height / 2).pow(2);

                    if a_dist <= b_dist {
                        Ordering::Less
                    } else {
                        Ordering::Greater
                    }
                });
            }
            3 => {
                // Could order each subsection with circular ordering
                let values = [16, 8, 4, 2, 1];

                for (n, value) in values.iter().enumerate() {
                    for i in (0..image_height).step_by(*value) {
                        for j in (0..image_width).step_by(*value) {
                            if n == 0 || i & (values[n - 1] - 1) != 0 || j & (values[n - 1] - 1) != 0 {
                                indices.push(i * image_width + j);
                            }
                        }
                    }
                }
            }
            _ => {}
        }

        // println!("generate indices took {}ms", time.elapsed().as_millis());

        indices
    }

    pub fn regenerate_from_settings(&mut self, settings: Config) {
        self.image_width = settings.get_int("image_width").unwrap_or(1000) as usize;
        self.image_height = settings.get_int("image_height").unwrap_or(1000) as usize;
        self.rotate = settings.get_float("rotate").unwrap_or(0.0).to_radians();
        self.maximum_iteration = settings.get_int("iterations").unwrap_or(1000) as usize;
        let initial_zoom = settings.get_str("zoom").unwrap_or_else(|_| String::from("1E0")).to_ascii_uppercase();
        let center_real = settings.get_str("real").unwrap_or_else(|_| String::from("-0.75"));
        let center_imag = settings.get_str("imag").unwrap_or_else(|_| String::from("0.0"));
        let approximation_order = settings.get_int("approximation_order").unwrap_or(0) as usize;
        self.glitch_percentage = settings.get_float("glitch_percentage").unwrap_or(0.001);
        self.remaining_frames = settings.get_int("frames").unwrap_or(1) as usize;
        self.frame_offset = settings.get_int("frame_offset").unwrap_or(0) as usize;
        self.zoom_scale_factor = settings.get_float("zoom_scale").unwrap_or(2.0);
        self.data_export.lock().display_glitches = settings.get_bool("display_glitches").unwrap_or(false);
        self.auto_adjust_iterations = settings.get_bool("auto_adjust_iterations").unwrap_or(true);
        self.experimental = settings.get_bool("experimental").unwrap_or(false);
        let probe_sampling = settings.get_int("probe_sampling").unwrap_or(3) as usize;
        self.remove_centre = settings.get_bool("remove_centre").unwrap_or(true);
        self.data_export.lock().iteration_division = settings.get_float("iteration_division").unwrap_or(0.1) as f32;
        self.data_export.lock().palette_offset = settings.get_float("palette_offset").unwrap_or(0.0) as f32;
        let valid_iteration_frame_multiplier = settings.get_float("valid_iteration_frame_multiplier").unwrap_or(0.25) as f32;
        let valid_iteration_probe_multiplier = settings.get_float("valid_iteration_probe_multiplier").unwrap_or(0.02) as f32;
        let glitch_tolerance = settings.get_float("glitch_tolerance").unwrap_or(1.4e-6) as f64;
        let data_storage_interval = settings.get_int("data_storage_interval").unwrap_or(10) as usize;
        self.analytic_derivative = settings.get_bool("analytic_derivative").unwrap_or(false);
        self.jitter = settings.get_bool("jitter").unwrap_or(false);
        self.jitter_factor = settings.get_float("jitter_factor").unwrap_or(0.2);
        self.show_output = settings.get_bool("show_output").unwrap_or(true);
        
        let data_type = match settings.get_str("export").unwrap_or_else(|_| String::from("COLOUR")).to_ascii_uppercase().as_ref() {
            "NONE" => DataType::None,
            "GUI" => DataType::Gui,
            "RAW" | "EXR" => DataType::Raw,
            "COLOUR" | "COLOR" | "PNG" => DataType::Color,
            "KFB" => DataType::Kfb,
            "BOTH" => DataType::Both,
            _ => DataType::Color
        };

        self.data_export.lock().palette = match data_type {
            DataType::Raw | DataType::None => {
                Vec::new()
            },
            DataType::Kfb | DataType::Color | DataType::Both | DataType::Gui => {
                if let Ok(colour_values) = settings.get_array("palette") {
                    colour_values.chunks_exact(3).map(|value| {
                        // We assume the palette is in BGR rather than RGB
                        (value[2].clone().into_int().unwrap() as u8, 
                            value[1].clone().into_int().unwrap() as u8, 
                            value[0].clone().into_int().unwrap() as u8)
                    }).collect::<Vec<(u8, u8, u8)>>()
                } else {
                    generate_default_palette()
                }
            }
        };

        let mut zoom = string_to_extended(&initial_zoom);
        let delta_pixel =  (-2.0 * (4.0 / self.image_height as f64 - 2.0) / zoom) / self.image_height as f64;
        let radius = delta_pixel * self.image_width as f64;
        let precision = max(64, -radius.exponent + 64);

        let center_location = ComplexArbitrary::with_val(
            precision as u32,
            ComplexArbitrary::parse("(".to_owned() + &center_real + "," + &center_imag + ")").expect("provided location not valid"));
        let auto_approximation = get_approximation_terms(approximation_order, self.image_width, self.image_height);

        self.center_reference = Reference::new(center_location.clone(), 
            center_location, 
            1, 
            self.maximum_iteration, 
            data_storage_interval,
            glitch_tolerance,
            zoom);

        self.series_approximation = SeriesApproximation::new_central(auto_approximation, 
            self.maximum_iteration, 
            FloatExtended::new(0.0, 0), 
            probe_sampling,
            self.experimental,
            valid_iteration_frame_multiplier,
            valid_iteration_probe_multiplier,
            data_storage_interval);

        self.total_pixels = self.image_width * self.image_height;

        if self.remove_centre {
            let temp = 1.0 / self.zoom_scale_factor;

            // Add one to avoid rescaling artifacts
            let val1 = (self.image_width as f64 * temp).ceil() as usize - 1;
            let val2 = (self.image_height as f64 * temp).ceil() as usize - 1;
    
            self.total_pixels -= val1 * val2;
        }

        let mut data_export = self.data_export.lock();

        if self.image_width != data_export.image_width || self.image_height != data_export.image_height {
            self.render_indices = FractalRenderer::generate_render_indices(self.image_width, self.image_height);
            data_export.centre_removed = false;
        }

        // Change the zoom level to the correct one for the frame offset
        for _ in 0..self.frame_offset {
            zoom.mantissa /= self.zoom_scale_factor;
            zoom.reduce();
        }

        self.zoom = zoom;

        self.start_render_time = Instant::now();
        self.progress = ProgressCounters::new(self.maximum_iteration);

        data_export.image_width = self.image_width;
        data_export.image_height = self.image_height;
        data_export.analytic_derivative = self.analytic_derivative;

        data_export.clear_buffers();
    }
}