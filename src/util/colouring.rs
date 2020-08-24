use crate::util::image::Image;
use crate::util::PixelData;
use crate::math::Reference;

pub enum DataExport {
    PNG,
    EXR,
    BOTH
}

impl DataExport {
    pub fn run(&self, pixel_data: &Vec<PixelData>, image: &mut Image, maximum_iteration: usize, reference: &Reference) {
        // Palette is temporarily here

        match self {
            DataExport::PNG => {
                let colours = self.generate_colour_palette();

                for pixel in pixel_data {
                    let (r, g, b) = if pixel.glitched && image.display_glitches {
                        (255, 0, 0)
                    } else if pixel.iteration >= maximum_iteration {
                        (0, 0, 0)
                    } else {
                        // 0.1656
                        let hue = (40 * pixel.iteration) % 1024;

                        let colour = colours[hue];

                        (colour.0 as u8, colour.1 as u8, colour.2 as u8)
                    };

                    image.plot_jpg(pixel.image_x, pixel.image_y, r, g, b);
                }
            },
            DataExport::EXR => {
                let escape_radius_ln = 1e16f32.ln();

                for pixel in pixel_data {
                    let z_norm = (reference.data[pixel.iteration - reference.start_iteration].z_fixed + pixel.delta_current.mantissa).norm_sqr() as f32;
                    let smooth = 1.0 - (z_norm.ln() / escape_radius_ln).log2();

                    image.plot_exr(pixel.image_x, pixel.image_y, pixel.iteration as u32, smooth);
                }
            },
            DataExport::BOTH => {
                return;
            }
        }
    }

    fn generate_colour_palette(&self) -> Vec<(f32, f32, f32)> {
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