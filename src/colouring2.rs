use crate::util::image::Image;
use crate::util::{PixelData, PixelData2};

pub enum ColourMethod2 {
    Iteration,
    IterationSquareRoot,
    Histogram,
    Distance
}

impl ColourMethod2 {
    pub fn run(&self, pixel_data: &Vec<PixelData2>, image: &mut Image, maximum_iteration: usize, delta_pixel: f64) {
        // Palette is temporarily here
        let mut colours = Vec::new();

        for i in 0..8192 {
            let value = i as f32 / 8192 as f32;

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

        match self {
            ColourMethod2::Iteration => {
                // No smooth colouring at the moment

                for pixel in pixel_data {
                    let (r, g, b) = if pixel.glitched && image.display_glitches {
                        (255, 0, 0)
                    } else if pixel.iteration >= maximum_iteration {
                        (0, 0, 0)
                    } else {
                        // 0.1656
                        let hue = (100.0 * pixel.iteration as f64) as usize % 8192;

                        let colour = colours[hue];

                        (colour.0 as u8, colour.1 as u8, colour.2 as u8)
                    };

                    image.plot(pixel.image_x, pixel.image_y, r, g, b);
                }
            },
            ColourMethod2::IterationSquareRoot => {
                for pixel in pixel_data {
                    let (r, g, b) = if pixel.glitched && image.display_glitches {
                        (255, 0, 0)
                    } else if pixel.iteration >= maximum_iteration {
                        (0, 0, 0)
                    } else {
                        // 0.1656
                        let hue = (7.0 * pixel.iteration as f64).sqrt() as usize % 8192;

                        let colour = colours[hue];

                        (colour.0 as u8, colour.1 as u8, colour.2 as u8)
                    };

                    image.plot(pixel.image_x, pixel.image_y, r, g, b);
                }
            },
            ColourMethod2::Histogram => {
                let mut iteration_counts = vec![0usize; maximum_iteration + 2];

                for pixel in pixel_data {
                    iteration_counts[pixel.iteration as usize] += 1
                }

                for i in 1..iteration_counts.len() {
                    iteration_counts[i] += iteration_counts[i - 1];
                }

                // Don't count the pixels that are inside the set
                let total = iteration_counts[maximum_iteration - 1];


                for pixel in pixel_data {
                    let (r, g, b) = if pixel.glitched && image.display_glitches {
                        (255, 0, 0)
                    } else if pixel.iteration >= maximum_iteration {
                        (0, 0, 0)
                    } else {
                        let v1 = iteration_counts[pixel.iteration] as f32 / total as f32;

                        // the hue is used to smooth the histogram bins. The hue is in the range 0.0-1.0
                        let hue = (v1 * 8192.0) as usize;

                        let colour = colours[hue % 8192];

                        (colour.0 as u8, colour.1 as u8, colour.2 as u8)
                    };

                    image.plot(pixel.image_x, pixel.image_y, r, g, b);
                }
            },
            ColourMethod2::Distance => {
                // At the moment distance has a white-black gradient
                for pixel in pixel_data {
                    let (r, g, b) = if pixel.glitched && image.display_glitches {
                        (255, 0, 0)
                    } else {
                        if pixel.escaped {
                            println!("{}", pixel.derivative_current.norm());
                            let de = 2.0 * pixel.delta_current.norm() * pixel.delta_current.norm().ln() / pixel.derivative_current.norm();
                            let out = (255.0 * (de / delta_pixel).tanh()) as u8;
                            (out, out, out)
                        } else {
                            (0, 0, 0)
                        }
                    };

                    image.plot(pixel.image_x, pixel.image_y, r, g, b);
                }
            }
        }
    }
}