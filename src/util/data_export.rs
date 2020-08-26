use crate::util::PixelData;
use crate::math::Reference;

use exr::prelude::simple_image;

pub enum DataType {
    COLOUR,
    RAW,
    BOTH
}

pub struct DataExport {
    image_width: usize,
    image_height: usize,
    rgb: Vec<u8>,
    palette: Vec<(f32, f32, f32)>,
    iterations: Vec<u32>,
    smooth: Vec<f32>,
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
                    palette: DataExport::generate_colour_palette(),
                    iterations: Vec::new(),
                    smooth: Vec::new(),
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
                    smooth: vec![0.0f32; image_width * image_height],
                    display_glitches,
                    data_type
                }
            },
            DataType::BOTH => {
                DataExport {
                    image_width,
                    image_height,
                    rgb: vec![0u8; image_width * image_height * 3],
                    palette: DataExport::generate_colour_palette(),
                    iterations: vec![0u32; image_width * image_height],
                    smooth: vec![0.0f32; image_width * image_height],
                    display_glitches,
                    data_type
                }
            }
        }
    }

    pub fn export_pixels(&mut self, pixel_data: &Vec<PixelData>, maximum_iteration: usize, reference: &Reference) {
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
                            let colour = self.palette[10 * pixel.iteration % 1024];
                            self.rgb[k] = colour.0 as u8;
                            self.rgb[k + 1] = colour.1 as u8;
                            self.rgb[k + 2] = colour.2 as u8;
                        }
                    } else if pixel.iteration >= maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                    } else {
                        let colour = self.palette[10 * pixel.iteration % 1024];
                        self.rgb[k] = colour.0 as u8;
                        self.rgb[k + 1] = colour.1 as u8;
                        self.rgb[k + 2] = colour.2 as u8;
                    };
                }
            },
            DataType::RAW => {
                let escape_radius_ln = 1e16f32.ln();

                for pixel in pixel_data {
                    let k = pixel.image_y * self.image_width + pixel.image_x;

                    self.iterations[k] = if pixel.glitched {
                        0x00000000
                    } else if pixel.iteration >= maximum_iteration {
                        0xFFFFFFFF
                    } else {
                        pixel.iteration as u32
                    };

                    let z_norm = (reference.perturbation_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth[k] = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
                }
            },
            DataType::BOTH => {
                let escape_radius_ln = 1e16f32.ln();

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
                            self.rgb[k] = colour.0 as u8;
                            self.rgb[k + 1] = colour.1 as u8;
                            self.rgb[k + 2] = colour.2 as u8;
                            self.iterations[k / 3] = 0x00000000
                        }
                    } else if pixel.iteration >= maximum_iteration {
                        self.rgb[k] = 0;
                        self.rgb[k + 1] = 0;
                        self.rgb[k + 2] = 0;
                        self.iterations[k / 3] = 0xFFFFFFFF;
                    } else {
                        let colour = self.palette[10 * pixel.iteration % 1024];
                        self.rgb[k] = colour.0 as u8;
                        self.rgb[k + 1] = colour.1 as u8;
                        self.rgb[k + 2] = colour.2 as u8;
                        self.iterations[k / 3] = pixel.iteration as u32;
                    };

                    let z_norm = (reference.perturbation_data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    self.smooth[k / 3] = 1.0 - (z_norm.ln() / escape_radius_ln).log2();
                }
            }
        }
    }

    pub fn save(&mut self, filename: &String) {
        match self.data_type {
            DataType::COLOUR => {
                self.save_colour(filename);
            },
            DataType::RAW => {
                self.save_raw(filename);
            },
            DataType::BOTH => {
                self.save_colour(filename);
                self.save_raw(filename);
            }
        }
    }

    fn save_colour(&mut self, filename: &String) {
        image::save_buffer(filename.to_owned() + ".jpg", &self.rgb, self.image_height as u32, self.image_height as u32, image::ColorType::Rgb8).unwrap();
    }

    fn save_raw(&mut self, filename: &String) {
        let iterations = simple_image::Channel::non_color_data(simple_image::Text::from("N").unwrap(), simple_image::Samples::U32(self.iterations.clone()));
        let smooth = simple_image::Channel::non_color_data(simple_image::Text::from("NF").unwrap(), simple_image::Samples::F32(self.smooth.clone()));

        let layer = simple_image::Layer::new(simple_image::Text::from("fractal_data").unwrap(), (self.image_width, self.image_height), smallvec::smallvec![iterations, smooth])
            .with_compression(simple_image::Compression::PXR24)
            .with_block_format(None, simple_image::attribute::LineOrder::Increasing);

        let image = simple_image::Image::new_from_single_layer(layer);

        image.write_to_file(filename.to_owned() + ".exr", simple_image::write_options::high()).unwrap();
    }

    fn generate_colour_palette() -> Vec<(f32, f32, f32)> {
        let mut colours = Vec::with_capacity(1024);

        for i in 0..1024 {
            let value = i as f32 / 1024 as f32;

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
                blue = 0.0 + factor * (0.0 - 0.0);
            } else {
                let factor = (value - 0.8575) / (1.0 - 0.8575);

                red = 0.0 + factor * (0.0 - 0.0);
                green = 2.0 + factor * (7.0 - 2.0);
                blue = 0.0 + factor * (100.0 - 0.0);
            }

            colours.push((red, green, blue))
        }

        colours
    }
}