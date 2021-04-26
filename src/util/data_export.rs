use crate::util::{PixelData, FloatExtended, ComplexFixed, FractalType};
use crate::math::Reference;

use std::{collections::HashMap, f64::consts::LN_2};
// use std::cmp::{min, max};
use std::f32::consts::{FRAC_PI_4, TAU};

use exr::{prelude::simple_image};
use colorgrad::{Color, CustomGradient, Interpolation, BlendMode};

// This is 1e16f32.ln().log2() + 1.0
const ESCAPE_RADIUS_LN_LOG2_P1: f32 = 5.203254472696 + 1.0;
// const ESCAPE_RADIUS_LN_LOG3_P1: f32 = 3.282888062227 + 1.0;

#[derive(PartialEq, Clone, Copy)]
pub enum ExportType {
    Color,
    Raw,
    Both,
    Gui
}

#[derive(PartialEq, Clone, Copy)]
pub enum ColoringType {
    SmoothIteration,
    StepIteration,
    Distance,
    DistanceStripe,
    Stripe
}

#[derive(PartialEq, Clone, Copy)]
pub enum DataType { 
    Iteration,
    Distance,
    Stripe,
    DistanceStripe,
    AtomDomain
}

pub struct LightingParameters {
    pub diffuse: [f32; 4],
    pub specular: [f32; 4],
    pub shininess: i32,
    pub opacity: [f32; 2],
    pub ambient: f32
}

impl LightingParameters {
    pub fn new(direction: f32, azimuth: f32, opacity: f32, ambient: f32, diffuse: f32, specular: f32, shininess: i32) -> Self { 
        let phi_half = FRAC_PI_4 + azimuth / 2.0;

        LightingParameters {
            diffuse: [direction.cos() * azimuth.cos(), direction.sin() * azimuth.cos(), azimuth.sin(), diffuse],
            specular: [direction.cos() * phi_half.sin(), direction.sin() * phi_half.sin(), phi_half.cos(), specular],
            shininess,
            opacity: [opacity, (1.0 - opacity) / 2.0],
            ambient
        } 
    }
}

pub struct DataExport {
    pub image_width: usize,
    pub image_height: usize,
    pub buffer: Vec<u8>,
    pub iterations: Vec<u32>,
    pub smooth: Vec<f32>,
    pub stripe: Vec<f32>,
    pub distance_x: Vec<f32>,
    pub distance_y: Vec<f32>,
    pub glitched: Vec<bool>,
    pub palette_buffer: Vec<Color>,
    pub palette_interpolated_buffer: Vec<Color>,
    pub palette_cyclic: bool,
    pub display_glitches: bool,
    pub palette_iteration_span: f32,
    pub palette_offset: f32,
    pub distance_transition: f32,
    pub centre_removed: bool,
    pub data_type: DataType,
    pub coloring_type: ColoringType,
    pub maximum_iteration: usize,
    pub fractal_type: FractalType,
    pub export_type: ExportType,
    pub lighting_parameters: LightingParameters,
    pub lighting: bool,
}

impl DataExport {
    pub fn new(image_width: usize, 
        image_height: usize, 
        display_glitches: bool, 
        palette_buffer: Vec<Color>, 
        palette_interpolated_buffer: Vec<Color>, 
        palette_cyclic: bool,
        palette_iteration_span: f32, 
        palette_offset: f32, 
        distance_transition: f32, 
        lighting: bool,
        coloring_type: ColoringType,
        data_type: DataType,
        fractal_type: FractalType, 
        export_type: ExportType) -> Self {

        // TODO maybe make the distance estimate arrays empty until required
        DataExport {
            image_width,
            image_height,
            buffer: vec![0u8; image_width * image_height * 3],
            iterations: vec![0u32; image_width * image_height],
            smooth: vec![0.0f32; image_width * image_height],
            stripe: vec![0.0f32; image_width * image_height],
            distance_x: vec![0.0f32; image_width * image_height],
            distance_y: vec![0.0f32; image_width * image_height],
            glitched: vec![false; image_width * image_height],
            palette_buffer,
            palette_interpolated_buffer,
            palette_cyclic,
            display_glitches,
            palette_iteration_span,
            palette_offset,
            distance_transition,
            centre_removed: false,
            data_type,
            coloring_type,
            maximum_iteration: 0,
            fractal_type,
            export_type,
            lighting_parameters: LightingParameters::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0),
            lighting
        }
    }

    #[inline]
    pub fn export_pixels(&mut self, pixel_data: &[PixelData], reference: &Reference, delta_pixel: FloatExtended, scale: usize) {
        for pixel in pixel_data {
            let new_scale = if self.export_type == ExportType::Gui {
                scale
            } else {
                1
            };

            self.glitched[pixel.index] = pixel.glitched;

            if pixel.glitched && self.display_glitches {
                self.set_with_scale(pixel.index, [255, 0, 0], new_scale);
                continue;
            }

            self.iterations[pixel.index] = pixel.iteration as u32;

            if pixel.iteration >= self.maximum_iteration {
                self.set_with_scale(pixel.index, [0, 0, 0], new_scale);
                continue;
            }

            let z = reference.reference_data[pixel.iteration - reference.start_iteration].z + pixel.delta_current.mantissa;
            let z_norm = z.norm_sqr();

            // self.smooth[pixel.index] = match self.fractal_type {
            //     FractalType::Mandelbrot2 => {
            //         ESCAPE_RADIUS_LN_LOG2_P1 - (z_norm.ln() as f32).log2()
            //     }
            //     FractalType::Mandelbrot3 => {
            //         ESCAPE_RADIUS_LN_LOG3_P1 - (z_norm.ln() as f32).log(3.0)
            //     }
            // };

            self.smooth[pixel.index] = ESCAPE_RADIUS_LN_LOG2_P1 - (z_norm.ln() as f32).log2();

            if self.data_type == DataType::Stripe || self.data_type == DataType::DistanceStripe {
                self.stripe[pixel.index] = (pixel.stripe.0 + pixel.stripe.1 * self.smooth[pixel.index] + pixel.stripe.2 * (1.0 - self.smooth[pixel.index])) / 4.0;
            }

            if self.data_type == DataType::Distance || self.data_type == DataType::DistanceStripe {
                // This calculates the distance in terms of pixels
                let temp1 = reference.reference_data_extended[pixel.iteration - reference.start_iteration] + pixel.delta_current;
                let temp2 = temp1.norm();

                let temp3 = 2.0f64.powi(temp1.exponent - temp2.exponent) / temp2.mantissa;

                let norm_z_x = temp1.mantissa.re * temp3;
                let norm_z_y = temp1.mantissa.im * temp3;

                let jxa = FloatExtended::new(pixel.derivative_current.mantissa.re, pixel.derivative_current.exponent);
                let jya = FloatExtended::new(pixel.derivative_current.mantissa.im, pixel.derivative_current.exponent);

                let scaled_jxa = (jxa * delta_pixel).to_float();
                let scaled_jya = (jya * delta_pixel).to_float();

                let num = (temp2 * (temp2.mantissa.ln() + temp2.exponent as f64 * LN_2)).to_float();

                let den_1 = norm_z_x * scaled_jxa + norm_z_y * scaled_jya;
                let den_2 = norm_z_x * -1.0 * scaled_jya + norm_z_y * scaled_jxa;

                let output = num / ComplexFixed::new(den_1, den_2);

                self.distance_x[pixel.index] = output.re as f32;
                self.distance_y[pixel.index] = output.im as f32;
            };

            self.colour_index(pixel.index, new_scale)
        }
    }

    pub fn save(&mut self, filename: &str, approximation_order: usize, zoom: &str) {
        match self.export_type {
            ExportType::Color => {
                self.save_colour(filename);
            },
            ExportType::Raw => {
                self.save_raw(filename, approximation_order, zoom);
            },
            ExportType::Both => {
                self.save_colour(filename);
                self.save_raw(filename, approximation_order, zoom);
            }
            _ => {},
        }
    }

    pub fn save_colour(&mut self, filename: &str) {
        // Extension is specified
        if let Some(extension) = filename.split_terminator('.').last() {
            match extension {
                "jpg" | "jpeg" | "png" => {
                    image::save_buffer(
                        filename.to_owned(), 
                        &self.buffer, 
                        self.image_width as u32, 
                        self.image_height as u32, 
                        image::ColorType::Rgb8).unwrap();

                    return;
                }
                _ => {}
            }
        }

        image::save_buffer(
            filename.to_owned() + ".png", 
            &self.buffer, 
            self.image_width as u32, 
            self.image_height as u32, 
            image::ColorType::Rgb8).unwrap();
    }

    pub fn save_raw(&mut self, filename: &str, approximation_order: usize, zoom: &str) {
        let iterations = simple_image::Channel::non_color_data(simple_image::Text::from("N").unwrap(), simple_image::Samples::U32(self.iterations.clone()));
        let smooth = simple_image::Channel::non_color_data(simple_image::Text::from("NF").unwrap(), simple_image::Samples::F32(self.smooth.clone()));

        let channels = if self.data_type == DataType::Distance {
            let distance_x = simple_image::Channel::non_color_data(simple_image::Text::from("DEX").unwrap(), simple_image::Samples::F32(self.distance_x.clone()));
            let distance_y = simple_image::Channel::non_color_data(simple_image::Text::from("DEY").unwrap(), simple_image::Samples::F32(self.distance_y.clone()));

            smallvec::smallvec![iterations, smooth, distance_x, distance_y]
        } else {
            smallvec::smallvec![iterations, smooth]
        };

        let mut layer = simple_image::Layer::new(simple_image::Text::from("fractal_data").unwrap(), (self.image_width, self.image_height), channels)
            .with_compression(simple_image::Compression::PXR24)
            .with_block_format(None, simple_image::attribute::LineOrder::Increasing);   

        let mut attributes = HashMap::new();
        attributes.insert(simple_image::Text::from("IterationsBias").unwrap(), exr::meta::attribute::AttributeValue::I32(0));
        attributes.insert(simple_image::Text::from("Iterations").unwrap(), exr::meta::attribute::AttributeValue::I32(self.maximum_iteration as i32));
        attributes.insert(simple_image::Text::from("Zoom").unwrap(), exr::meta::attribute::AttributeValue::Text(simple_image::Text::from(zoom).unwrap()));
        attributes.insert(simple_image::Text::from("approximation_order").unwrap(), exr::meta::attribute::AttributeValue::I32(approximation_order as i32));

        layer.attributes = exr::meta::header::LayerAttributes::new(simple_image::Text::from("fractal_data").unwrap());
        layer.attributes.custom = attributes;

        let image = simple_image::Image::new_from_single_layer(layer);

        image.write_to_file(filename.to_owned() + ".exr", simple_image::write_options::high()).unwrap();
    }

    pub fn clear_buffers(&mut self) {
        self.buffer = vec![0u8; self.image_width * self.image_height * 3];
        self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
        self.smooth = vec![0.0f32; self.image_width * self.image_height];
        self.stripe = vec![0.0f32; self.image_width * self.image_height];
        self.distance_x = vec![0.0f32; self.image_width * self.image_height];
        self.distance_y = vec![0.0f32; self.image_width * self.image_height];
        self.glitched = vec![false; self.image_width * self.image_height];
    }

    pub fn regenerate(&mut self) {
        for i in 0..self.iterations.len() {
            if self.glitched[i] && self.display_glitches {
                self.set_with_scale(i, [255, 0, 0], 1);
                continue;
            }

            if self.iterations[i] >= self.maximum_iteration as u32 {
                self.set_with_scale(i, [0, 0, 0], 1);
                continue;
            }

            self.colour_index(i, 1);
        }
    }

    // pub fn interpolate_glitches(&mut self, pixel_data: &[PixelData]) {
    //     if !self.display_glitches {
    //         for pixel in pixel_data {
    //             let k = pixel.image_y * self.image_width + pixel.image_x;

    //             let k_up = (max(1, pixel.image_y) - 1) * self.image_width + pixel.image_x;
    //             let k_down = (min(self.image_height - 2, pixel.image_y) + 1) * self.image_width + pixel.image_x;

    //             let k_left = pixel.image_y * self.image_width + max(1, pixel.image_x) - 1;
    //             let k_right = pixel.image_y * self.image_width + min(self.image_width - 2, pixel.image_x) + 1;

    //             self.iterations[k] = (self.iterations[k_up] + self.iterations[k_down] + self.iterations[k_left] + self.iterations[k_right]) / 4;
    //             self.smooth[k] = (self.smooth[k_up] + self.smooth[k_down] + self.smooth[k_left] + self.smooth[k_right]) / 4.0;

    //             if self.data_type == DataType::Distance {
    //                 self.distance_x[k] = (self.distance_x[k_up] + self.distance_x[k_down] + self.distance_x[k_left] + self.distance_x[k_right]) / 4.0;
    //                 self.distance_y[k] = (self.distance_y[k_up] + self.distance_y[k_down] + self.distance_y[k_left] + self.distance_y[k_right]) / 4.0;
    //             }

    //             self.colour_index(k, 1);
    //         }
    //     }
    // }

    #[inline]
    pub fn colour_index(&mut self, k: usize, scale: usize) {
        let rgb: [u8; 3] = match self.coloring_type {
            ColoringType::Distance => {
                // TODO have some alternative algorithms
                // Calculate distance estimate in terms of pixels
                let length_pixel = (self.distance_x[k].powi(2) + self.distance_y[k].powi(2)).sqrt();

                // Scale so transition will happen about 50 pixels
                let length_scaled = (length_pixel / self.distance_transition).max(0.0);

                // colouring algorithm based on 'rainbow_fringe' by claude
                let angle = self.distance_y[k].atan2(self.distance_x[k]);
    
                let mut hue = angle / TAU;
                hue -= hue.floor();
    
                let saturation = (1.0 / (1.0 + length_scaled)).clamp(0.0, 1.0);
                let value = length_scaled.clamp(0.0, 1.0);
    
                let color = Color::from_hsv(hue as f64 * 360.0, saturation as f64, value as f64);
                let (r, g, b, _) = color.rgba_u8();
    
                [r, g, b]
            },
            ColoringType::SmoothIteration | ColoringType::StepIteration => {
                let mut floating_iteration = self.iterations[k] as f32 / self.palette_iteration_span;
                
                if self.coloring_type == ColoringType::SmoothIteration {
                    floating_iteration += self.smooth[k] / self.palette_iteration_span
                };
                
                let temp = self.palette_interpolated_buffer.len() as f32 * (floating_iteration + self.palette_offset).fract();
    
                let pos1 = temp.floor() as usize;
                let pos2 = (pos1 + 1) % self.palette_interpolated_buffer.len();
    
                let frac = temp.fract() as f64;
    
                let (r, g, b, _) = self.palette_interpolated_buffer[pos1].interpolate_rgb(&self.palette_interpolated_buffer[pos2], frac).rgba_u8();
    
                [r, g, b]
            },
            ColoringType::Stripe => {
                let mut floating_iteration = self.iterations[k] as f32 / self.palette_iteration_span;
                
                floating_iteration += self.smooth[k] / self.palette_iteration_span;
                
                let temp = self.palette_interpolated_buffer.len() as f32 * (floating_iteration + self.palette_offset).fract();
    
                let pos1 = temp.floor() as usize;
                let pos2 = (pos1 + 1) % self.palette_interpolated_buffer.len();
    
                let frac = temp.fract() as f64;
    
                let (r, g, b, _) = self.palette_interpolated_buffer[pos1].interpolate_rgb(&self.palette_interpolated_buffer[pos2], frac).rgba();

                let bright = self.stripe[k] as f64;

                let (r, g, b) = if bright < 0.5 {
                    ((2.0 * r * bright), (2.0 * g * bright), (2.0 * b * bright))
                } else {
                    (1.0 - 2.0 * (1.0 - r) * (1.0 - bright), 1.0 - 2.0 * (1.0 - g) * (1.0 - bright), 1.0 - 2.0 * (1.0 - b) * (1.0 - bright))
                };

                let (r, g, b, _) = Color::from_rgb(r, g, b).rgba_u8();

                [r, g, b]
            },
            ColoringType::DistanceStripe => {
                // TODO, it could be possible to have some kind of continous distance estimate - which is based on the zoom level as well; so that keyframes can be added together
                let floating_iteration = (self.iterations[k] as f32 + self.smooth[k]) / self.palette_iteration_span;
                
                let temp = self.palette_interpolated_buffer.len() as f32 * (floating_iteration + self.palette_offset).fract();
    
                let pos1 = temp.floor() as usize;
                let pos2 = (pos1 + 1) % self.palette_interpolated_buffer.len();
    
                let frac = temp.fract() as f64;
    
                let (r, g, b, _) = self.palette_interpolated_buffer[pos1].interpolate_rgb(&self.palette_interpolated_buffer[pos2], frac).rgba();

                // Blinn-phong from GPU mandelbrot
                let mut normal = ComplexFixed::new(self.distance_x[k], self.distance_y[k]);
                normal /= normal.norm();

                let bright = if self.lighting {
                    // This is diffuse lighting
                    let light_diffuse = (normal.re * self.lighting_parameters.diffuse[0] + normal.im * self.lighting_parameters.diffuse[1] + self.lighting_parameters.diffuse[2]) / (1.0 + self.lighting_parameters.diffuse[2]);

                    // This is specular lighting
                    let light_specular = ((normal.re * self.lighting_parameters.specular[0] + normal.im * self.lighting_parameters.specular[1] + self.lighting_parameters.specular[2]) / (1.0 + self.lighting_parameters.specular[2])).powi(self.lighting_parameters.shininess);

                    // Ambient + diffuse + specular
                    let bright = self.lighting_parameters.ambient + self.lighting_parameters.diffuse[3] * light_diffuse + self.lighting_parameters.specular[3] * light_specular;

                    // Add intensity
                    bright * self.lighting_parameters.opacity[0] - self.lighting_parameters.opacity[1]
                } else {
                    0.5
                };

                // Calculate distance estimate in terms of pixels
                let length_pixel = (self.distance_x[k].powi(2) + self.distance_y[k].powi(2)).sqrt();

                // Scale so transition will happen about 50 pixels
                let length_scaled = (length_pixel / self.distance_transition).max(0.0);

                // Apply sigmoid function for non-linear transform
                let value = if length_scaled < 1.0 {
                    1.0 - 1.0 / (1.0 + (-10.0 * (length_scaled - 0.5)).exp())
                } else {
                    0.0
                };

                let temp = if self.stripe[k] < 0.5 {
                    2.0 * self.stripe[k] * bright
                } else {
                    1.0 - 2.0 * (1.0 - bright) * (1.0 - self.stripe[k])
                };

                let bright = (temp * (1.0 - value) + bright * value) as f64;

                let (r, g, b) = if bright < 0.5 {
                    ((2.0 * r * bright), (2.0 * g * bright), (2.0 * b * bright))
                } else {
                    (1.0 - 2.0 * (1.0 - r) * (1.0 - bright), 1.0 - 2.0 * (1.0 - g) * (1.0 - bright), 1.0 - 2.0 * (1.0 - b) * (1.0 - bright))
                };

                let (r, g, b, _) = Color::from_rgb(r, g, b).rgba_u8();
    
                [r, g, b]
            }
        };

        self.set_with_scale(k, rgb, scale)
    }

    #[inline]
    pub fn change_palette(&mut self, palette: Option<Vec<(u8, u8, u8)>>, palette_iteration_span: f32, palette_offset: f32, distance_transition: f32, cyclic: bool, lighting: bool) {
        let mut new_palette = false;
        
        if let Some(palette) = palette {
            self.palette_buffer = palette.iter().map(|value| {
                Color::from_rgb_u8(value.0, value.1, value.2)
            }).collect::<Vec<Color>>();

            if self.palette_buffer[0] != *self.palette_buffer.last().unwrap() {
                self.palette_buffer.push(self.palette_buffer[0].clone());
            };

            new_palette = true;
        };

        if new_palette || (cyclic != self.palette_cyclic) {
            let mut number_colors = self.palette_buffer.len();

            if !cyclic {
                number_colors -= 1;
            }

            let palette_generator = CustomGradient::new()
                .colors(&self.palette_buffer[0..number_colors])
                .interpolation(Interpolation::CatmullRom)
                .mode(BlendMode::Oklab)
                .build().unwrap();

            self.palette_interpolated_buffer = palette_generator.colors(number_colors * 64);
        };

        self.palette_iteration_span = palette_iteration_span;
        self.palette_offset = palette_offset;
        self.palette_cyclic = cyclic;
        self.lighting = lighting;
        self.distance_transition = distance_transition;
    }

    #[inline]
    pub fn change_lighting(&mut self, direction: f32, azimuth: f32, opacity: f32, ambient: f32, diffuse: f32, specular: f32, shininess: i32) {
        self.lighting_parameters = LightingParameters::new(direction, azimuth, opacity, ambient, diffuse, specular, shininess);
    }

    #[inline]
    pub fn set_with_scale(&mut self, index: usize, value: [u8; 3], scale: usize) {
        if scale > 1 {
            let image_x = index % self.image_width;
            let image_y = index / self.image_width;

            let horizontal = if image_x + scale < self.image_width {
                scale
            } else {
                self.image_width - image_x
            };

            let vertical = if image_y + scale < self.image_height {
                scale
            } else {
                self.image_height - image_y
            };

            if self.data_type == DataType::DistanceStripe {
                for i in image_x..(image_x + horizontal) {
                    for j in image_y..(image_y + vertical) {
                        let scale_index = j * self.image_width + i;

                        self.iterations[scale_index] = self.iterations[index];
                        self.smooth[scale_index] = self.smooth[index];
                        self.glitched[scale_index] = self.glitched[index];
                        self.distance_x[scale_index] = self.distance_x[index];
                        self.distance_y[scale_index] = self.distance_y[index];
                        self.stripe[scale_index] = self.stripe[index];
    
                        self.buffer[3 * (scale_index)] = value[0];
                        self.buffer[3 * (scale_index) + 1] = value[1];
                        self.buffer[3 * (scale_index) + 2] = value[2];
                    }
                }
            } else if self.data_type == DataType::Distance {
                for i in image_x..(image_x + horizontal) {
                    for j in image_y..(image_y + vertical) {
                        let scale_index = j * self.image_width + i;

                        self.iterations[scale_index] = self.iterations[index];
                        self.smooth[scale_index] = self.smooth[index];
                        self.glitched[scale_index] = self.glitched[index];
                        self.distance_x[scale_index] = self.distance_x[index];
                        self.distance_y[scale_index] = self.distance_y[index];
    
                        self.buffer[3 * (scale_index)] = value[0];
                        self.buffer[3 * (scale_index) + 1] = value[1];
                        self.buffer[3 * (scale_index) + 2] = value[2];
                    }
                }
            } else if self.data_type == DataType::Stripe {
                for i in image_x..(image_x + horizontal) {
                    for j in image_y..(image_y + vertical) {
                        let scale_index = j * self.image_width + i;

                        self.iterations[scale_index] = self.iterations[index];
                        self.smooth[scale_index] = self.smooth[index];
                        self.glitched[scale_index] = self.glitched[index];
                        self.stripe[scale_index] = self.stripe[index];
    
                        self.buffer[3 * (scale_index)] = value[0];
                        self.buffer[3 * (scale_index) + 1] = value[1];
                        self.buffer[3 * (scale_index) + 2] = value[2];
                    }
                }
            } else {
                for i in image_x..(image_x + horizontal) {
                    for j in image_y..(image_y + vertical) {
                        let scale_index = j * self.image_width + i;

                        self.iterations[scale_index] = self.iterations[index];
                        self.smooth[scale_index] = self.smooth[index];
                        self.glitched[scale_index] = self.glitched[index];
    
                        self.buffer[3 * (scale_index)] = value[0];
                        self.buffer[3 * (scale_index) + 1] = value[1];
                        self.buffer[3 * (scale_index) + 2] = value[2];
                    }
                }
            }
        } else {
            self.buffer[3 * index] = value[0];
            self.buffer[3 * index + 1] = value[1];
            self.buffer[3 * index + 2] = value[2];
        }
    }
}