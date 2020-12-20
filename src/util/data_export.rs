use crate::util::PixelData;
use crate::math::Reference;
use crate::util::FloatExtended;
use crate::util::ComplexFixed;

use std::io::prelude::*;
use std::fs::File;
use std::slice;
use std::collections::HashMap;
use std::cmp::{min, max};
use std::f32::consts::TAU;

use exr::prelude::simple_image;
use color_space::{Rgb, Hsv};

pub enum DataType {
    NONE,
    GUI,
    COLOUR,
    RAW,
    KFB,
    BOTH
}

pub struct DataExport {
    pub image_width: usize,
    pub image_height: usize,
    pub rgb: Vec<u8>,
    pub palette: Vec<(u8, u8, u8)>,
    pub iterations: Vec<u32>,
    pub smooth: Vec<f32>,
    pub distance_x: Vec<f32>,
    pub distance_y: Vec<f32>,
    pub display_glitches: bool,
    pub iteration_division: f32,
    pub palette_offset: f32,
    palette_length: usize,
    scaled_offset: f32,
    data_type: DataType,
    pub analytic_derivative: bool,
    pub maximum_iteration: usize
}

impl DataExport {
    pub fn new(image_width: usize, image_height: usize, display_glitches: bool, data_type: DataType, palette: Vec<(u8, u8, u8)>, iteration_division: f32, palette_offset: f32, analytic_derivative: bool) -> Self {
        let mut rgb = Vec::new();
        let mut smooth = Vec::new();
        let mut distance_x = Vec::new();
        let mut distance_y = Vec::new();

        match data_type {
            DataType::NONE => {},
            DataType::COLOUR => {
                rgb = vec![0u8; image_width * image_height * 3];
            },
            DataType::GUI => {
                rgb = vec![0u8; image_width * image_height * 3];
                smooth = vec![0.0f32; image_width * image_height];
                distance_x = vec![0.0f32; image_width * image_height];
                distance_y = vec![0.0f32; image_width * image_height];
            }
            DataType::RAW => {
                smooth = vec![0.0f32; image_width * image_height];
                distance_x = vec![0.0f32; image_width * image_height];
                distance_y = vec![0.0f32; image_width * image_height];
            },
            DataType::KFB => {
                smooth = vec![0.0f32; image_width * image_height];
            },
            DataType::BOTH => {
                rgb = vec![0u8; image_width * image_height * 3];
                smooth = vec![0.0f32; image_width * image_height];
                distance_x = vec![0.0f32; image_width * image_height];
                distance_y = vec![0.0f32; image_width * image_height];
            }
        }

        let palette_length = palette.len();

        DataExport {
            image_width,
            image_height,
            rgb,
            palette,
            iterations: vec![0u32; image_width * image_height],
            smooth,
            distance_x,
            distance_y,
            display_glitches,
            iteration_division,
            palette_offset,
            palette_length,
            scaled_offset: palette_offset * palette_length as f32,
            data_type,
            analytic_derivative,
            maximum_iteration: 0
        }
    }

    pub fn export_pixels(&mut self, pixel_data: &[PixelData], reference: &Reference, delta_pixel: FloatExtended) {
        let escape_radius_ln = 1e16f32.ln();

        match self.data_type {
            DataType::NONE => {},
            DataType::COLOUR => {
                for pixel in pixel_data {
                    let k = (pixel.image_y * self.image_width + pixel.image_x) * 3;

                    if pixel.glitched && self.display_glitches {
                        self.rgb[k] = 255;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                    } else if pixel.iteration >= self.maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                    } else {
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        let temp = ((pixel.iteration as f32 + smooth) / self.iteration_division) + self.scaled_offset;

                        let temp2 = temp.floor() as usize % self.palette_length;
                        let temp3 = (temp as usize + 1) % self.palette_length;
                        let temp4 = temp.fract();

                        let colour1 = self.palette[temp2];
                        let colour2 = self.palette[temp3];

                        self.rgb[k] = (colour1.0 as f32 + temp4 * (colour2.0 as f32 - colour1.0 as f32)) as u8; 
                        self.rgb[k + 1] = (colour1.1 as f32 + temp4 * (colour2.1 as f32 - colour1.1 as f32)) as u8; 
                        self.rgb[k + 2] = (colour1.2 as f32 + temp4 * (colour2.2 as f32 - colour1.2 as f32)) as u8;

                        if self.analytic_derivative && pixel.iteration < self.maximum_iteration {
                            // let temp = pixel.delta_current.norm();
                            let temp = (reference.reference_data[pixel.iteration - reference.start_iteration].z_extended + pixel.delta_current).norm();

                            let de = 2.0 * temp * (temp.mantissa.ln() + temp.exponent as f64 * 2.0f64.ln()) / pixel.derivative_current.norm();

                            let out = (255.0 * (de / delta_pixel).to_float().tanh()) as u8;

                            self.rgb[k] = out; 
                            self.rgb[k + 1] = out; 
                            self.rgb[k + 2] = out;
                        };
                    };
                };
            },
            DataType::GUI => {
                for pixel in pixel_data {
                    let k = pixel.image_y * self.image_width + pixel.image_x;

                    if pixel.glitched && self.display_glitches {
                        self.rgb[3 * k] = 255;
                        self.rgb[3 * k + 1] = 0;
                        self.rgb[3 * k + 2] = 0;
                    } else {
                        self.iterations[k] = pixel.iteration as u32;

                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        self.smooth[k] = smooth;

                        if self.analytic_derivative {
                            let temp1 = reference.reference_data[pixel.iteration - reference.start_iteration].z_extended + pixel.delta_current;
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

                        self.colour_index(k)
                    };
                };
            },
            DataType::RAW => {
                for pixel in pixel_data {
                    let k = pixel.image_y * self.image_width + pixel.image_x;

                    self.iterations[k] = if pixel.glitched {
                        0x00000000
                    } else if pixel.iteration >= self.maximum_iteration {
                        0xFFFFFFFF
                    } else {
                        pixel.iteration as u32
                    };

                    let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth[k] = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                    if self.analytic_derivative && pixel.iteration < self.maximum_iteration {
                        let temp1 = reference.reference_data[pixel.iteration - reference.start_iteration].z_extended + pixel.delta_current;
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
                }
            },
            DataType::KFB => {
                for pixel in pixel_data {
                    let k = pixel.image_x * self.image_height + pixel.image_y;

                    self.iterations[k] = if pixel.glitched {
                        0x00000000
                    } else if pixel.iteration >= self.maximum_iteration {
                        0xFFFFFFFF
                    } else {
                        pixel.iteration as u32
                    };

                    let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth[k] = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
                }
            },
            DataType::BOTH => {
                for pixel in pixel_data {
                    let k = pixel.image_y * self.image_width + pixel.image_x;

                    if pixel.glitched && self.display_glitches {
                        self.rgb[3 * k] = 255;
                        self.rgb[3 * k + 1] = 0;
                        self.rgb[3 * k + 2] = 0;
                        self.iterations[k] = 0x00000000;
                    } else if pixel.iteration >= self.maximum_iteration {
                        self.rgb[3 * k] = 0;
                        self.rgb[3 * k + 1] = 0;
                        self.rgb[3 * k + 2] = 0;
                        self.iterations[k] = 0xFFFFFFFF;
                    } else {
                        self.iterations[k] = pixel.iteration as u32;

                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        self.smooth[k] = smooth;

                        if self.analytic_derivative {
                            let temp1 = reference.reference_data[pixel.iteration - reference.start_iteration].z_extended + pixel.delta_current;
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

                        self.colour_index(k)
                    };
                };
            }
        }
    }

    pub fn save(&mut self, filename: &str, approximation_order: usize, zoom: &str) {
        match self.data_type {
            DataType::NONE | DataType::GUI => {},
            DataType::COLOUR => {
                self.save_colour(filename);
            },
            DataType::RAW => {
                self.save_raw(filename, approximation_order, zoom);
            },
            DataType::KFB => {
                self.save_kfb(filename);
            }
            DataType::BOTH => {
                self.save_colour(filename);
                self.save_raw(filename, approximation_order, zoom);
            }
        }
    }

    pub fn save_colour(&mut self, filename: &str) {
        // Extension is specified
        match filename.split_terminator(".").last() {
            Some(extension) => {
                match extension {
                    "jpg" | "jpeg" | "png" => {
                        image::save_buffer(
                            filename.to_owned(), 
                            &self.rgb, 
                            self.image_width as u32, 
                            self.image_height as u32, 
                            image::ColorType::Rgb8).unwrap();

                        return;
                    }
                    _ => {}
                }
                
            }
            _ => {}
        }

        image::save_buffer(
            filename.to_owned() + ".png", 
            &self.rgb, 
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

    pub fn save_kfb(&mut self, filename: &str) {
        let mut file = File::create(filename.to_owned() + ".kfb").unwrap();

        file.write_all(b"KFB").unwrap();

        let test1 = [self.image_width as u32, self.image_height as u32];

        // iteration division??
        let test3 = [1u32, self.palette_length as u32];

        // Maxmimum iteration
        let test6 = [self.maximum_iteration as u32];

        file.write_all(unsafe {
            slice::from_raw_parts(test1.as_ptr() as *const u8, 4)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(self.iterations.as_ptr() as *const u8, self.iterations.len() * 4)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(test3.as_ptr() as *const u8, 4)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(self.palette.as_ptr() as *const u8, 3 * self.palette_length)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(test6.as_ptr() as *const u8, 4)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(self.smooth.as_ptr() as *const u8, self.iterations.len() * 4)
        }).unwrap();

    }

    pub fn clear_buffers(&mut self) {
        match self.data_type {
            DataType::NONE => {},
            DataType::COLOUR => {
                self.rgb = vec![0u8; self.image_width * self.image_height * 3];
            }
            DataType::GUI => {
                self.rgb = vec![0u8; self.image_width * self.image_height * 3];
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
                self.distance_x = vec![0.0f32; self.image_width * self.image_height];
                self.distance_y = vec![0.0f32; self.image_width * self.image_height];
            }
            DataType::RAW => {
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
                self.distance_x = vec![0.0f32; self.image_width * self.image_height];
                self.distance_y = vec![0.0f32; self.image_width * self.image_height];
            },
            DataType::KFB => {
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
            },
            DataType::BOTH => {
                self.rgb = vec![0u8; self.image_width * self.image_height * 3];
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth = vec![0.0f32; self.image_width * self.image_height];
                self.distance_x = vec![0.0f32; self.image_width * self.image_height];
                self.distance_y = vec![0.0f32; self.image_width * self.image_height];
            }
        }
    }

    pub fn regenerate(&mut self) {
        match self.data_type {
            DataType::GUI => {
                for i in 0..self.iterations.len() {
                    self.colour_index(i);
                }
            },
            _ => {}
        }
    }

    pub fn interpolate_glitches(&mut self, pixel_data: &[PixelData]) {
        if !self.display_glitches {
            match self.data_type {
                DataType::GUI => {
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

                        self.colour_index(k);
                    }
                },
                _ => {}
            }
        }
    }

    pub fn colour_index(&mut self, i: usize) {
        if self.iterations[i] >= self.maximum_iteration as u32 {
            self.rgb[3 * i] = 0u8; 
            self.rgb[3 * i + 1] = 0u8; 
            self.rgb[3 * i + 2] = 0u8;
        } else if self.analytic_derivative {
            let length = (self.distance_x[i].powi(2) + self.distance_y[i].powi(2)).sqrt();
            
            // colouring algorithm based on 'rainbow_fringe' by claude
            let angle = self.distance_y[i].atan2(self.distance_x[i]);

            let mut hue = angle / TAU;
            hue -= hue.floor();

            let saturation = (1.0 / (1.0 + length)).max(0.0).min(1.0);
            let value = length.max(0.0).min(1.0);

            let hsv = Hsv::new(hue as f64 * 360.0, saturation as f64, value as f64);
            let rgb = Rgb::from(hsv);

            self.rgb[3 * i] = rgb.r as u8; 
            self.rgb[3 * i + 1] = rgb.g as u8; 
            self.rgb[3 * i + 2] = rgb.b as u8;

            // default colouring algorithm
            // let out = (255.0 * length) as u8;

            // self.rgb[3 * i] = out as u8; 
            // self.rgb[3 * i + 1] = out as u8; 
            // self.rgb[3 * i + 2] = out as u8;
        } else {
            let temp = ((self.iterations[i] as f32 + self.smooth[i]) / self.iteration_division) + self.scaled_offset;

            let temp2 = temp.floor() as usize % self.palette_length;
            let temp3 = (temp as usize + 1) % self.palette_length;
            let temp4 = temp.fract();

            let colour1 = self.palette[temp2];
            let colour2 = self.palette[temp3];

            self.rgb[3 * i] = (colour1.0 as f32 + temp4 * (colour2.0 as f32 - colour1.0 as f32)) as u8; 
            self.rgb[3 * i + 1] = (colour1.1 as f32 + temp4 * (colour2.1 as f32 - colour1.1 as f32)) as u8; 
            self.rgb[3 * i + 2] = (colour1.2 as f32 + temp4 * (colour2.2 as f32 - colour1.2 as f32)) as u8;
        }
    }

    pub fn change_palette(&mut self, palette: Option<Vec<(u8, u8, u8)>>, iteration_division: f32, palette_offset: f32) {
        match palette {
            Some(palette) => {
                self.palette = palette;
                self.palette_length = self.palette.len();
            },
            None => {}
        };

        self.iteration_division = iteration_division;
        self.palette_offset = palette_offset;
        self.scaled_offset = palette_offset * self.palette_length as f32;
    }
}