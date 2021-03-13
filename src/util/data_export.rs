use crate::util::{PixelData, FloatExtended, ComplexFixed};
use crate::math::Reference;

use std::collections::HashMap;
use std::cmp::{min, max};
use std::f32::consts::TAU;

use exr::{prelude::simple_image};
use colorgrad::{Color, CustomGradient, Interpolation, BlendMode, Gradient};

use super::FractalType;

#[derive(PartialEq, Clone, Copy)]
pub enum DataType {
    None,
    Gui,
    Color,
    Raw,
    Both
}

pub struct DataExport {
    pub image_width: usize,
    pub image_height: usize,
    pub buffer: Vec<u8>,
    pub iterations: Vec<u32>,
    pub smooth: Vec<f32>,
    pub distance_x: Vec<f32>,
    pub distance_y: Vec<f32>,
    pub glitched: Vec<bool>,
    pub palette_generator: Gradient,
    pub palette_buffer: Vec<Color>,
    pub display_glitches: bool,
    pub iteration_division: f32,
    pub palette_offset: f32,
    pub centre_removed: bool,
    pub data_type: DataType,
    pub analytic_derivative: bool,
    pub maximum_iteration: usize,
    pub fractal_type: FractalType
}

impl DataExport {
    pub fn new(image_width: usize, image_height: usize, display_glitches: bool, data_type: DataType, palette_generator: Gradient, palette_buffer: Vec<Color>, iteration_division: f32, palette_offset: f32, analytic_derivative: bool, fractal_type: FractalType) -> Self {
        let mut buffer = Vec::new();
        let mut smooth = Vec::new();
        let mut distance_x = Vec::new();
        let mut distance_y = Vec::new();
        let mut glitched = Vec::new();

        match data_type {
            DataType::None => {},
            DataType::Color => {
                buffer = vec![0u8; image_width * image_height * 3];
            },
            DataType::Gui => {
                buffer = vec![0u8; image_width * image_height * 3];
                smooth = vec![0.0f32; image_width * image_height];
                distance_x = vec![0.0f32; image_width * image_height];
                distance_y = vec![0.0f32; image_width * image_height];
                glitched = vec![false; image_width * image_height];
            }
            DataType::Raw => {
                smooth = vec![0.0f32; image_width * image_height];
                distance_x = vec![0.0f32; image_width * image_height];
                distance_y = vec![0.0f32; image_width * image_height];
            },
            DataType::Both => {
                buffer = vec![0u8; image_width * image_height * 3];
                smooth = vec![0.0f32; image_width * image_height];
                distance_x = vec![0.0f32; image_width * image_height];
                distance_y = vec![0.0f32; image_width * image_height];
            }
        }

        DataExport {
            image_width,
            image_height,
            buffer,
            iterations: vec![0u32; image_width * image_height],
            smooth,
            distance_x,
            distance_y,
            glitched,
            palette_generator,
            palette_buffer,
            display_glitches,
            iteration_division,
            palette_offset,
            centre_removed: false,
            data_type,
            analytic_derivative,
            maximum_iteration: 0,
            fractal_type
        }
    }

    pub fn export_pixels(&mut self, pixel_data: &[PixelData], reference: &Reference, delta_pixel: FloatExtended, scale: usize) {
        let escape_radius_ln = 1e16f32.ln();

        for pixel in pixel_data {
            match self.data_type {
                DataType::None => {},
                DataType::Color => {
                    let k = (pixel.image_y * self.image_width + pixel.image_x) * 3;
    
                    if pixel.glitched && self.display_glitches {
                        self.buffer[k] = 255;
                        self.buffer[k + 1] = 0;
                        self.buffer[k + 2] = 0;
                    } else if pixel.iteration >= self.maximum_iteration {
                        self.buffer[k] = 0;
                        self.buffer[k + 1] = 0;
                        self.buffer[k + 2] = 0;
                    } else {
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
    
                        let (r, g, b) = if self.analytic_derivative && pixel.iteration < self.maximum_iteration {
                            // let temp = pixel.delta_current.norm();
                            let temp = (reference.reference_data_extended[pixel.iteration - reference.start_iteration] + pixel.delta_current).norm();
    
                            let de = 2.0 * temp * (temp.mantissa.ln() + temp.exponent as f64 * 2.0f64.ln()) / pixel.derivative_current.norm();
    
                            let out = (255.0 * (de / delta_pixel).to_float().tanh()) as u8;
    
                            (out, out, out)
                        } else {
                            let temp = self.palette_buffer.len() as f32 * ((pixel.iteration as f32 + smooth) / self.iteration_division + self.palette_offset).fract();

                            let pos1 = temp.floor() as usize;
                            let pos2 = if pos1 == (self.palette_buffer.len() - 1) {
                                0
                            } else {
                                pos1 + 1
                            };
    
                            let frac = temp.fract() as f64;

                            let out = self.palette_buffer[pos1].interpolate_rgb(&self.palette_buffer[pos2], frac).rgba_u8();
                            (out.0, out.1, out.2)
                        };

                        self.buffer[3 * k] = r;
                        self.buffer[3 * k + 1] = g;
                        self.buffer[3 * k + 2] = b;
                    };
                },
                DataType::Gui => {
                    let k = pixel.image_y * self.image_width + pixel.image_x;
    
                    if pixel.glitched {
                        self.glitched[k] = true;
                    } else {
                        self.glitched[k] = false;
                    };
    
                    if pixel.glitched && self.display_glitches {    
                        self.set_rgb_with_scale(k, [255, 0, 0], scale, pixel.image_x, pixel.image_y);
                    } else {
                        self.iterations[k] = pixel.iteration as u32;
    
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z + pixel.delta_current.mantissa).norm_sqr();

                        self.smooth[k] = match self.fractal_type {
                            FractalType::Mandelbrot2 => {
                                1.0 - (z_norm.ln() as f32 / escape_radius_ln).log2()
                            }
                            FractalType::Mandelbrot3 => {
                                1.0 - (z_norm.ln() as f32 / escape_radius_ln).log(3.0)
                            }
                        };

                        if self.analytic_derivative {
                            let temp1 = reference.reference_data_extended[pixel.iteration - reference.start_iteration] + pixel.delta_current;
                            let temp2 = temp1.norm();
    
                            let norm_z_x = FloatExtended::new(
                                temp1.mantissa.re / temp2.mantissa,
                                temp1.exponent - temp2.exponent
                            ).to_float();
    
                            let norm_z_y = FloatExtended::new(
                                temp1.mantissa.im / temp2.mantissa,
                                temp1.exponent - temp2.exponent
                            ).to_float();
    
                            let jxa = FloatExtended::new(pixel.derivative_current.mantissa.re, pixel.derivative_current.exponent);
                            let jya = FloatExtended::new(pixel.derivative_current.mantissa.im, pixel.derivative_current.exponent);
    
                            let scaled_jxa = (jxa * delta_pixel).to_float();
                            let scaled_jya = (jya * delta_pixel).to_float();
    
                            let num = (temp2 * (temp2.mantissa.ln() + temp2.exponent as f64 * 2.0f64.ln())).to_float();
    
                            let den_1 = norm_z_x * scaled_jxa + norm_z_y * scaled_jya;
                            let den_2 = norm_z_x * -1.0 * scaled_jya + norm_z_y * scaled_jxa;
    
                            let output = num / ComplexFixed::new(den_1, den_2);
    
                            self.distance_x[k] = output.re as f32;
                            self.distance_y[k] = output.im as f32;
                        };
    
                        self.colour_index(k, scale, pixel.image_x, pixel.image_y)
                    };
                },
                DataType::Raw => {
                    let k = pixel.image_y * self.image_width + pixel.image_x;
    
                    self.iterations[k] = if pixel.glitched {
                        0x00000000
                    } else if pixel.iteration >= self.maximum_iteration {
                        0xFFFFFFFF
                    } else {
                        pixel.iteration as u32
                    };
    
                    let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth[k] = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
    
                    if self.analytic_derivative && pixel.iteration < self.maximum_iteration {
                        let temp1 = reference.reference_data_extended[pixel.iteration - reference.start_iteration] + pixel.delta_current;
                        let temp2 = temp1.norm();
    
                        let norm_z_x = FloatExtended::new(
                            temp1.mantissa.re / temp2.mantissa,
                            temp1.exponent - temp2.exponent
                        ).to_float();
    
                        let norm_z_y = FloatExtended::new(
                            temp1.mantissa.im / temp2.mantissa,
                            temp1.exponent - temp2.exponent
                        ).to_float();
    
                        let jxa = FloatExtended::new(pixel.derivative_current.mantissa.re, pixel.derivative_current.exponent);
                        let jya = FloatExtended::new(pixel.derivative_current.mantissa.im, pixel.derivative_current.exponent);
    
                        let scaled_jxa = (jxa * delta_pixel).to_float();
                        let scaled_jya = (jya * delta_pixel).to_float();
    
                        let num = (temp2 * (temp2.mantissa.ln() + temp2.exponent as f64 * 2.0f64.ln())).to_float();
    
                        let den_1 = norm_z_x * scaled_jxa + norm_z_y * scaled_jya;
                        let den_2 = norm_z_x * -1.0 * scaled_jya + norm_z_y * scaled_jxa;
    
                        let output = num / ComplexFixed::new(den_1, den_2);
    
                        self.distance_x[k] = output.re as f32;
                        self.distance_y[k] = output.im as f32;
                    }
                },
                DataType::Both => {
                    let k = pixel.image_y * self.image_width + pixel.image_x;
    
                    if pixel.glitched && self.display_glitches {
                        self.buffer[3 * k] = 255;
                        self.buffer[3 * k + 1] = 0;
                        self.buffer[3 * k + 2] = 0;
                        self.iterations[k] = 0x00000000;
                    } else if pixel.iteration >= self.maximum_iteration {
                        self.buffer[3 * k] = 0;
                        self.buffer[3 * k + 1] = 0;
                        self.buffer[3 * k + 2] = 0;
                        self.iterations[k] = 0xFFFFFFFF;
                    } else {
                        self.iterations[k] = pixel.iteration as u32;
    
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
    
                        self.smooth[k] = smooth;
    
                        if self.analytic_derivative {
                            let temp1 = reference.reference_data_extended[pixel.iteration - reference.start_iteration] + pixel.delta_current;
                            let temp2 = temp1.norm();
    
                            let norm_z_x = FloatExtended::new(
                                temp1.mantissa.re / temp2.mantissa,
                                temp1.exponent - temp2.exponent
                            ).to_float();
    
                            let norm_z_y = FloatExtended::new(
                                temp1.mantissa.im / temp2.mantissa,
                                temp1.exponent - temp2.exponent
                            ).to_float();
    
                            let jxa = FloatExtended::new(pixel.derivative_current.mantissa.re, pixel.derivative_current.exponent);
                            let jya = FloatExtended::new(pixel.derivative_current.mantissa.im, pixel.derivative_current.exponent);
    
                            let scaled_jxa = (jxa * delta_pixel).to_float();
                            let scaled_jya = (jya * delta_pixel).to_float();
    
                            let num = (temp2 * (temp2.mantissa.ln() + temp2.exponent as f64 * 2.0f64.ln())).to_float();
    
                            let den_1 = norm_z_x * scaled_jxa + norm_z_y * scaled_jya;
                            let den_2 = norm_z_x * -1.0 * scaled_jya + norm_z_y * scaled_jxa;
    
                            let output = num / ComplexFixed::new(den_1, den_2);
    
                            self.distance_x[k] = output.re as f32;
                            self.distance_y[k] = output.im as f32;
                        };
    
                        self.colour_index(k, 1, 0, 0)
                    };
                }
            }
        }
    }

    pub fn save(&mut self, filename: &str, approximation_order: usize, zoom: &str) {
        match self.data_type {
            DataType::None | DataType::Gui => {},
            DataType::Color => {
                self.save_colour(filename);
            },
            DataType::Raw => {
                self.save_raw(filename, approximation_order, zoom);
            },
            DataType::Both => {
                self.save_colour(filename);
                self.save_raw(filename, approximation_order, zoom);
            }
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

        let channels = if self.analytic_derivative {
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
        match self.data_type {
            DataType::None => {},
            DataType::Color => {
                self.buffer = vec![0u8; self.image_width * self.image_height * 3];
            }
            DataType::Gui => {
                self.buffer = vec![0u8; self.image_width * self.image_height * 3];
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
                self.distance_x = vec![0.0f32; self.image_width * self.image_height];
                self.distance_y = vec![0.0f32; self.image_width * self.image_height];
                self.glitched = vec![false; self.image_width * self.image_height];
            }
            DataType::Raw => {
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
                self.distance_x = vec![0.0f32; self.image_width * self.image_height];
                self.distance_y = vec![0.0f32; self.image_width * self.image_height];
            },
            DataType::Both => {
                self.buffer = vec![0u8; self.image_width * self.image_height * 3];
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
                self.distance_x = vec![0.0f32; self.image_width * self.image_height];
                self.distance_y = vec![0.0f32; self.image_width * self.image_height];
            }
        }
    }

    pub fn regenerate(&mut self) {
        if self.data_type == DataType::Gui {
            for i in 0..self.iterations.len() {
                self.colour_index(i, 1, 0, 0);
            }
        }
    }

    pub fn interpolate_glitches(&mut self, pixel_data: &[PixelData]) {
        if !self.display_glitches && self.data_type == DataType::Gui {
            for pixel in pixel_data {
                let k = pixel.image_y * self.image_width + pixel.image_x;

                let k_up = (max(1, pixel.image_y) - 1) * self.image_width + pixel.image_x;
                let k_down = (min(self.image_height - 2, pixel.image_y) + 1) * self.image_width + pixel.image_x;

                let k_left = pixel.image_y * self.image_width + max(1, pixel.image_x) - 1;
                let k_right = pixel.image_y * self.image_width + min(self.image_width - 2, pixel.image_x) + 1;

                self.iterations[k] = (self.iterations[k_up] + self.iterations[k_down] + self.iterations[k_left] + self.iterations[k_right]) / 4;
                self.smooth[k] = (self.smooth[k_up] + self.smooth[k_down] + self.smooth[k_left] + self.smooth[k_right]) / 4.0;

                if self.analytic_derivative {
                    self.distance_x[k] = (self.distance_x[k_up] + self.distance_x[k_down] + self.distance_x[k_left] + self.distance_x[k_right]) / 4.0;
                    self.distance_y[k] = (self.distance_y[k_up] + self.distance_y[k_down] + self.distance_y[k_left] + self.distance_y[k_right]) / 4.0;
                }

                self.colour_index(k, 1, 0, 0);
            }
        }
    }

    #[inline]
    pub fn colour_index(&mut self, i: usize, scale: usize, image_x: usize, image_y: usize) {
        if self.glitched[i] && self.display_glitches {
            self.set_rgb_with_scale(i, [255, 0, 0], scale, image_x, image_y)
        } else if self.iterations[i] >= self.maximum_iteration as u32 {
            self.set_rgb_with_scale(i, [0, 0, 0], scale, image_x, image_y)
        } else if self.analytic_derivative {
            let length = (self.distance_x[i].powi(2) + self.distance_y[i].powi(2)).sqrt();
            
            // colouring algorithm based on 'rainbow_fringe' by claude
            let angle = self.distance_y[i].atan2(self.distance_x[i]);

            let mut hue = angle / TAU;
            hue -= hue.floor();

            let saturation = (1.0 / (1.0 + length)).clamp(0.0, 1.0);
            let value = length.clamp(0.0, 1.0);

            let color = Color::from_hsv(hue as f64 * 360.0, saturation as f64, value as f64);
            let (r, g, b, _) = color.rgba_u8();

            self.set_rgb_with_scale(i, [r, g, b], scale, image_x, image_y)
        } else {
            let temp = self.palette_buffer.len() as f32 * ((self.iterations[i] as f32 + self.smooth[i]) / self.iteration_division + self.palette_offset).fract();

            let pos1 = temp.floor() as usize;
            let pos2 = if pos1 >= (self.palette_buffer.len() - 1) {
                0
            } else {
                pos1 + 1
            };

            let frac = temp.fract() as f64;

            let (r, g, b, _) = self.palette_buffer[pos1].interpolate_rgb(&self.palette_buffer[pos2], frac).rgba_u8();

            self.set_rgb_with_scale(i, [r, g, b], scale, image_x, image_y)
        }
    }

    #[inline]
    pub fn change_palette(&mut self, palette: Option<Vec<(u8, u8, u8)>>, iteration_division: f32, palette_offset: f32) {
        if let Some(palette) = palette {
            let color = palette.iter().map(|value| {
                // We assume the palette is in BGR rather than RGB
                Color::from_rgb_u8(value.0, value.1, value.2)
            }).collect::<Vec<Color>>();

            let palette_generator = CustomGradient::new()
                .colors(&color)
                .interpolation(Interpolation::CatmullRom)
                .mode(BlendMode::Oklab)
                .build().unwrap();

            self.palette_buffer = palette_generator.colors(color.len() * 64);
            self.palette_generator = palette_generator;
        };

        self.iteration_division = iteration_division;
        self.palette_offset = palette_offset;
    }

    pub fn set_rgb_with_scale(&mut self, index: usize, value: [u8; 3], scale: usize, image_x: usize, image_y: usize) {
        if scale > 1 {
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

            for i in image_x..(image_x + horizontal) {
                for j in image_y..(image_y + vertical) {
                    self.buffer[3 * (j * self.image_width + i)] = value[0];
                    self.buffer[3 * (j * self.image_width + i) + 1] = value[1];
                    self.buffer[3 * (j * self.image_width + i) + 2] = value[2];
                }
            }
        } else {
            self.buffer[3 * index] = value[0];
            self.buffer[3 * index + 1] = value[1];
            self.buffer[3 * index + 2] = value[2];
        }
    }
}