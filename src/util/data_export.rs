use crate::util::PixelData;
use crate::math::Reference;

use std::io::prelude::*;
use std::fs::File;
use std::slice;
use std::collections::HashMap;

use exr::prelude::simple_image;
use half::f16;

pub enum DataType {
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
    pub display_glitches: bool,
    data_type: DataType,
}

impl DataExport {
    pub fn new(image_width: usize, image_height: usize, display_glitches: bool, data_type: DataType) -> Self {
        match data_type {
            DataType::COLOUR => {
                DataExport {
                    image_width,
                    image_height,
                    rgb: vec![0u8; image_width * image_height * 3],
                    palette: DataExport::generate_default_palette(),
                    iterations: vec![0u32; image_width * image_height],
                    smooth_f16: Vec::new(),
                    smooth_f32: Vec::new(),
                    display_glitches,
                    data_type
                }
            },
            DataType::RAW => {
                DataExport {
                    image_width,
                    image_height,
                    rgb: Vec::new(),
                    palette: Vec::new(),
                    iterations: vec![0u32; image_width * image_height],
                    smooth_f16: vec![f16::ZERO; image_width * image_height],
                    smooth_f32: Vec::new(),
                    display_glitches,
                    data_type
                }
            },
            DataType::KFB => {
                DataExport {
                    image_width,
                    image_height,
                    rgb: Vec::new(),
                    palette: Vec::new(),
                    iterations: vec![0u32; image_width * image_height],
                    smooth_f16: Vec::new(),
                    smooth_f32: vec![0.0f32; image_width * image_height],
                    display_glitches,
                    data_type
                }
            },
            DataType::BOTH => {
                DataExport {
                    image_width,
                    image_height,
                    rgb: vec![0u8; image_width * image_height * 3],
                    palette: DataExport::generate_default_palette(),
                    iterations: vec![0u32; image_width * image_height],
                    smooth_f16: vec![f16::ZERO; image_width * image_height],
                    smooth_f32: Vec::new(),
                    display_glitches,
                    data_type
                }
            }
        }
    }

    pub fn export_pixels(&mut self, pixel_data: &[PixelData], maximum_iteration: usize, reference: &Reference) {
        let escape_radius_ln = 1e16f32.ln();

        match self.data_type {
            DataType::COLOUR => {
                for pixel in pixel_data {
                    let k = (pixel.image_y * self.image_width + pixel.image_x) * 3;

                    if pixel.glitched {
                        if self.display_glitches {
                            self.rgb[k] = 255;
                            self.rgb[k + 1] = 0;
                            self.rgb[k + 2] = 0;
                        } else {
                            // Pixel is glitched so the smooth won't do much

                            let colour = self.palette[10 * pixel.iteration % 1024];
                            self.rgb[k] = colour.0;
                            self.rgb[k + 1] = colour.1;
                            self.rgb[k + 2] = colour.2;
                        }
                    } else if pixel.iteration >= maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                    } else {
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;

                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        let colour = self.palette[(10.0 * (pixel.iteration as f32 + smooth)) as usize % 1024];
                        self.rgb[k] = colour.0;
                        self.rgb[k + 1] = colour.1;
                        self.rgb[k + 2] = colour.2;
                    };

                    self.iterations[k / 3] = pixel.iteration as u32;
                }
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

                    if pixel.glitched {
                        if self.display_glitches {
                            self.rgb[k] = 255;
                            self.rgb[k + 1] = 0;
                            self.rgb[k + 2] = 0;
                            self.iterations[k / 3] = 0x00000000
                        } else {
                            let colour = self.palette[10 * pixel.iteration % 1024];
                            self.rgb[k] = colour.0;
                            self.rgb[k + 1] = colour.1;
                            self.rgb[k + 2] = colour.2;
                            self.iterations[k / 3] = 0x00000000
                        }
                    } else if pixel.iteration >= maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                        self.iterations[k / 3] = 0xFFFFFFFF;
                    } else {
                        let z_norm = (reference.reference_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;

                        let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                        let colour = self.palette[(10.0 * (pixel.iteration as f32 + smooth)) as usize % 1024];
                        self.rgb[k] = colour.0;
                        self.rgb[k + 1] = colour.1;
                        self.rgb[k + 2] = colour.2;
                        self.iterations[k / 3] = pixel.iteration as u32;
                        self.smooth_f16[k / 3] = f16::from_f32(smooth)
                    };
                }
            }
        }
    }

    pub fn save(&mut self, filename: &str, maximum_iteration: usize, approximation_order: usize, zoom: &str) {
        match self.data_type {
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

        let mut layer = simple_image::Layer::new(simple_image::Text::from("fractal_data").unwrap(), (self.image_width, self.image_height), smallvec::smallvec![iterations, smooth])
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

        // Colours in colourmap
        let test5 = DataExport::generate_default_palette();

        // iteration division??
        let test3 = [1u32, test5.len() as u32];

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
            slice::from_raw_parts(test5.as_ptr() as *const u8, 3 * test5.len())
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(test6.as_ptr() as *const u8, 4)
        }).unwrap();

        file.write_all(unsafe {
            slice::from_raw_parts(self.smooth_f32.as_ptr() as *const u8, self.iterations.len() * 4)
        }).unwrap();

    }

    fn generate_default_palette() -> Vec<(u8, u8, u8)> {
        let mut colours = Vec::with_capacity(1024);

        for i in 0..1024 {
            let value = i as f32 / 1024.0;

            let red;
            let green;
            let blue;

            if value < 0.16 {
                let factor = (value - 0.0) / (0.16 - 0.0);

                red = 0.0 + factor * (32.0 - 0.0);
                green = 7.0 + factor * (107.0 - 7.0);
                blue = 100.0 + factor * (203.0 - 100.0);
            } else if value < 0.42 {
                let factor = (value - 0.16) / (0.42 - 0.16);

                red = 32.0 + factor * (237.0 - 32.0);
                green = 107.0 + factor * (255.0 - 107.0);
                blue = 203.0 + factor * (255.0 - 203.0);
            } else if value < 0.6425 {
                let factor = (value - 0.42) / (0.6425 - 0.42);

                red = 237.0 + factor * (255.0 - 237.0);
                green = 255.0 + factor * (170.0 - 255.0);
                blue = 255.0 + factor * (0.0 - 255.0);
            } else if value < 0.8575 {
                let factor = (value - 0.6425) / (0.8575 - 0.6425);

                red = 255.0 + factor * (0.0 - 255.0);
                green = 170.0 + factor * (2.0 - 170.0);
                blue = 0.0;
            } else {
                let factor = (value - 0.8575) / (1.0 - 0.8575);

                red = 0.0;
                green = 2.0 + factor * (7.0 - 2.0);
                blue = 0.0 + factor * (100.0 - 0.0);
            }

            colours.push((red as u8, green as u8, blue as u8))
        }

        colours
    }

    pub fn clear_buffers(&mut self) {
        match self.data_type {
            DataType::COLOUR => {
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