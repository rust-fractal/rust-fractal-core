use crate::util::PixelData;
use crate::math::Reference;
use crate::util::FloatExtended;
use crate::util::ComplexFixed;

use std::io::prelude::*;
use std::fs::File;
use std::slice;
use std::collections::HashMap;

use exr::prelude::simple_image;
use half::f16;

pub enum DataType {
    NONE,
    GUI,
    COLOUR,
    RAW,
    KFB,
    BOTH
}

pub struct DataExport {
    image_width: usize,
    image_height: usize,
    pub rgb: Vec<u8>,
    palette: Vec<(u8, u8, u8)>,
    pub iterations: Vec<u32>,
    pub smooth_f16: Vec<f16>,
    pub smooth_f32: Vec<f32>,
    pub distance_x: Vec<f16>,
    pub distance_y: Vec<f16>,
    pub display_glitches: bool,
    iteration_division: f32,
    data_type: DataType,
    analytic_derivative: bool
}

impl DataExport {
    pub fn new(image_width: usize, image_height: usize, display_glitches: bool, data_type: DataType, palette: Vec<(u8, u8, u8)>, iteration_division: f32, analytic_derivative: bool) -> Self {
        let mut rgb = Vec::new();
        let mut smooth_f16 = Vec::new();
        let mut smooth_f32 = Vec::new();
        let mut distance_x = Vec::new();
        let mut distance_y = Vec::new();

        match data_type {
            DataType::NONE => {},
            DataType::COLOUR | DataType::GUI => {
                rgb = vec![0u8; image_width * image_height * 3];
            },
            DataType::RAW => {
                smooth_f16 = vec![f16::ZERO; image_width * image_height];
                distance_x = vec![f16::ZERO; image_width * image_height];
                distance_y = vec![f16::ZERO; image_width * image_height];
            },
            DataType::KFB => {
                smooth_f32 = vec![0.0f32; image_width * image_height];
            },
            DataType::BOTH => {
                rgb = vec![0u8; image_width * image_height * 3];
                smooth_f16 = vec![f16::ZERO; image_width * image_height];
            }
        }

        DataExport {
            image_width,
            image_height,
            rgb,
            palette,
            iterations: vec![0u32; image_width * image_height],
            smooth_f16,
            smooth_f32,
            distance_x,
            distance_y,
            display_glitches,
            iteration_division,
            data_type,
            analytic_derivative
        }
    }

    pub fn export_pixels(&mut self, pixel_data: &[PixelData], maximum_iteration: usize, reference: &Reference, delta_pixel: FloatExtended) {
        let escape_radius_ln = 1e16f32.ln();

        match self.data_type {
            DataType::NONE => {},
            DataType::COLOUR | DataType::GUI => {
                for pixel in pixel_data {
                    let k = (pixel.image_y * self.image_width + pixel.image_x) * 3;

                    if pixel.glitched && self.display_glitches {
                        self.rgb[k] = 255;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                    } else if pixel.iteration >= maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                    } else {
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        let temp = (pixel.iteration as f32 + smooth) / self.iteration_division;

                        let temp2 = temp.floor() as usize % self.palette.len();
                        let temp3 = (temp as usize + 1) % self.palette.len();
                        let temp4 = temp.fract();

                        let colour1 = self.palette[temp2];
                        let colour2 = self.palette[temp3];

                        self.rgb[k] = (colour1.0 as f32 + temp4 * (colour2.0 as f32 - colour1.0 as f32)) as u8; 
                        self.rgb[k + 1] = (colour1.1 as f32 + temp4 * (colour2.1 as f32 - colour1.1 as f32)) as u8; 
                        self.rgb[k + 2] = (colour1.2 as f32 + temp4 * (colour2.2 as f32 - colour1.2 as f32)) as u8;

                        if self.analytic_derivative && pixel.iteration < maximum_iteration {
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
            DataType::RAW => {
                for pixel in pixel_data {
                    let k = pixel.image_y * self.image_width + pixel.image_x;

                    self.iterations[k] = if pixel.glitched {
                        0x00000000
                    } else if pixel.iteration >= maximum_iteration {
                        0xFFFFFFFF
                    } else {
                        pixel.iteration as u32
                    };

                    let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth_f16[k] = f16::from_f32(1.0 - (z_norm.ln() / escape_radius_ln).log2());

                    if self.analytic_derivative && pixel.iteration < maximum_iteration {
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

                        self.distance_x[k] = f16::from_f64(output.re);
                        self.distance_y[k] = f16::from_f64(output.im);
                    }
                }
            },
            DataType::KFB => {
                for pixel in pixel_data {
                    let k = pixel.image_x * self.image_height + pixel.image_y;

                    self.iterations[k] = if pixel.glitched {
                        0x00000000
                    } else if pixel.iteration >= maximum_iteration {
                        0xFFFFFFFF
                    } else {
                        pixel.iteration as u32
                    };

                    let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth_f32[k] = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
                }
            },
            DataType::BOTH => {
                for pixel in pixel_data {
                    let k = (pixel.image_y * self.image_width + pixel.image_x) * 3;

                    if pixel.glitched && self.display_glitches {
                        self.rgb[k] = 255;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                        self.iterations[k / 3] = 0x00000000;
                    } else if pixel.iteration >= maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                        self.iterations[k / 3] = 0xFFFFFFFF;
                    } else {
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        let temp = (pixel.iteration as f32 + smooth) / self.iteration_division;

                        let temp2 = temp.floor() as usize % self.palette.len();
                        let temp3 = (temp as usize + 1) % self.palette.len();
                        let temp4 = temp.fract();

                        let colour1 = self.palette[temp2];
                        let colour2 = self.palette[temp3];

                        self.rgb[k] = (colour1.0 as f32 + temp4 * (colour2.0 as f32 - colour1.0 as f32)) as u8; 
                        self.rgb[k + 1] = (colour1.1 as f32 + temp4 * (colour2.1 as f32 - colour1.1 as f32)) as u8; 
                        self.rgb[k + 2] = (colour1.2 as f32 + temp4 * (colour2.2 as f32 - colour1.2 as f32)) as u8;
                        self.iterations[k / 3] = pixel.iteration as u32;
                        self.smooth_f16[k / 3] = f16::from_f32(smooth)
                    };
                }
            }
        }
    }

    pub fn save(&mut self, filename: &str, maximum_iteration: usize, approximation_order: usize, zoom: &str) {
        match self.data_type {
            DataType::NONE | DataType::GUI => {},
            DataType::COLOUR => {
                self.save_colour(filename);
            },
            DataType::RAW => {
                self.save_raw(filename, maximum_iteration, approximation_order, zoom);
            },
            DataType::KFB => {
                self.save_kfb(filename, maximum_iteration);
            }
            DataType::BOTH => {
                self.save_colour(filename);
                self.save_raw(filename, maximum_iteration, approximation_order, zoom);
            }
        }
    }

    fn save_colour(&mut self, filename: &str) {
        image::save_buffer(
            filename.to_owned() + ".png", 
            &self.rgb, 
            self.image_width as u32, 
            self.image_height as u32, 
            image::ColorType::Rgb8).unwrap();
    }

    fn save_raw(&mut self, filename: &str, maximum_iteration: usize, approximation_order: usize, zoom: &str) {
        let iterations = simple_image::Channel::non_color_data(simple_image::Text::from("N").unwrap(), simple_image::Samples::U32(self.iterations.clone()));
        let smooth = simple_image::Channel::non_color_data(simple_image::Text::from("NF").unwrap(), simple_image::Samples::F16(self.smooth_f16.clone()));

        let channels = if self.analytic_derivative {
            let distance_x = simple_image::Channel::non_color_data(simple_image::Text::from("DEX").unwrap(), simple_image::Samples::F16(self.distance_x.clone()));
            let distance_y = simple_image::Channel::non_color_data(simple_image::Text::from("DEY").unwrap(), simple_image::Samples::F16(self.distance_y.clone()));

            smallvec::smallvec![iterations, smooth, distance_x, distance_y]
        } else {
            smallvec::smallvec![iterations, smooth]
        };

        let mut layer = simple_image::Layer::new(simple_image::Text::from("fractal_data").unwrap(), (self.image_width, self.image_height), channels)
            .with_compression(simple_image::Compression::PXR24)
            .with_block_format(None, simple_image::attribute::LineOrder::Increasing);   

        let mut attributes = HashMap::new();
        attributes.insert(simple_image::Text::from("IterationsBias").unwrap(), exr::meta::attribute::AttributeValue::I32(0));
        attributes.insert(simple_image::Text::from("Iterations").unwrap(), exr::meta::attribute::AttributeValue::I32(maximum_iteration as i32));
        attributes.insert(simple_image::Text::from("Zoom").unwrap(), exr::meta::attribute::AttributeValue::Text(simple_image::Text::from(zoom).unwrap()));
        attributes.insert(simple_image::Text::from("approximation_order").unwrap(), exr::meta::attribute::AttributeValue::I32(approximation_order as i32));

        layer.attributes = exr::meta::header::LayerAttributes::new(simple_image::Text::from("fractal_data").unwrap());
        layer.attributes.custom = attributes;

        let image = simple_image::Image::new_from_single_layer(layer);

        image.write_to_file(filename.to_owned() + ".exr", simple_image::write_options::high()).unwrap();
    }

    fn save_kfb(&mut self, filename: &str, maximum_iteration: usize) {
        let mut file = File::create(filename.to_owned() + ".kfb").unwrap();

        file.write_all(b"KFB").unwrap();

        let test1 = [self.image_width as u32, self.image_height as u32];

        // iteration division??
        let test3 = [1u32, self.palette.len() as u32];

        // Maxmimum iteration
        let test6 = [maximum_iteration as u32];

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
            slice::from_raw_parts(self.palette.as_ptr() as *const u8, 3 * self.palette.len())
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(test6.as_ptr() as *const u8, 4)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(self.smooth_f32.as_ptr() as *const u8, self.iterations.len() * 4)
        }).unwrap();

    }

    pub fn clear_buffers(&mut self) {
        match self.data_type {
            DataType::NONE => {},
            DataType::COLOUR | DataType::GUI => {
                self.rgb = vec![0u8; self.image_width * self.image_height * 3];
            }
            DataType::RAW => {
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth_f16 = vec![f16::ZERO; self.image_width * self.image_height];
            },
            DataType::KFB => {
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth_f32 = vec![0.0f32; self.image_width * self.image_height];
            },
            DataType::BOTH => {
                self.rgb = vec![0u8; self.image_width * self.image_height * 3];
                self.iterations = vec![0xFFFFFFFF; self.image_width * self.image_height];
                self.smooth_f16 = vec![f16::ZERO; self.image_width * self.image_height];
            }
        }
    }
}