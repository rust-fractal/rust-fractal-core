use crate::util::image::Image;
use crate::util::PixelData;

pub enum ColourMethod {
    Iteration,
    // IterationSquareRoot,
    // Histogram
}

impl ColourMethod {
    pub fn run(&self, pixels: &Vec<PixelData>, image: &mut Image, maximum_iterations: usize) {
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
            ColourMethod::Iteration => {
                for pixel in pixels {
                    let (red, green, blue) = if pixel.glitched && image.display_glitches {
                        (255, 0, 0)
                    } else if pixel.iteration >= maximum_iterations {
                        (0, 0, 0)
                    } else {
                        // 0.1656
                        let hue = (0.1656 * pixel.iteration as f64) as usize % 8192;

                        let colour = colours[hue];

                        (colour.0 as u8, colour.1 as u8, colour.2 as u8)
                    };

                    image.plot(pixel.image_x, pixel.image_y, red, green, blue);
                }
            },
            // ColourMethod::IterationSquareRoot => {
            //     for point in points {
            //         let index = point.index;
            //
            //         if point.glitched && display_glitches {
            //             image[3 * index] = 255u8;
            //             image[3 * index + 1] = 0u8;
            //             image[3 * index + 2] = 0u8;
            //         } else if point.iterations >= maximum_iterations {
            //             image[3 * index] = 0u8;
            //             image[3 * index + 1] = 0u8;
            //             image[3 * index + 2] = 0u8;
            //         } else {
            //             let hue = 1600.0 * (point.iterations as f32 + point.smooth).sqrt() % 8192.0;
            //
            //             let colour = colours[(hue.floor() as usize) % 8192];
            //             let colour2 = colours[(hue.floor() as usize + 1) % 8192];
            //
            //             let red = (colour.0 + ((colour2.0 - colour.0) * hue.fract())) as u8;
            //             let green = (colour.1 + ((colour2.1 - colour.1) * hue.fract())) as u8;
            //             let blue = (colour.2 + ((colour2.2 - colour.2) * hue.fract())) as u8;
            //
            //             image[3 * index] = red;
            //             image[3 * index + 1] = green;
            //             image[3 * index + 2] = blue;
            //         }
            //     }
            // },
            // ColourMethod::Histogram => {
            //     let mut iteration_counts = vec![0usize; maximum_iterations + 1];
            //
            //     for point in points {
            //         iteration_counts[point.iterations as usize] += 1
            //     }
            //
            //     for i in 1..iteration_counts.len() {
            //         iteration_counts[i] += iteration_counts[i - 1];
            //     }
            //
            //     let total = iteration_counts[maximum_iterations - 1];
            //
            //     for point in points {
            //         let index = point.index;
            //
            //         if point.glitched && display_glitches {
            //             image[3 * index] = 255u8;
            //             image[3 * index + 1] = 0u8;
            //             image[3 * index + 2] = 0u8;
            //         } else if point.iterations >= maximum_iterations {
            //             image[3 * index] = 0u8;
            //             image[3 * index + 1] = 0u8;
            //             image[3 * index + 2] = 0u8;
            //         } else {
            //             let factor = if point.smooth == std::f32::NAN || point.iterations as f32 + point.smooth < 0.0 {
            //                 point.iterations as f32
            //             } else {
            //                 point.iterations as f32 + point.smooth
            //             };
            //
            //             let v1 = iteration_counts[factor as usize] as f32 / total as f32;
            //             let v2 = iteration_counts[factor as usize + 1] as f32 / total as f32;
            //
            //             // the hue is used to smooth the histogram bins. The hue is in the range 0.0-1.0
            //             let hue = (v1 + (v2 - v1) * factor.fract()) * 8192.0;
            //
            //             let colour = colours[hue.floor() as usize % 8192];
            //             let colour2 = colours[(hue.floor() as usize + 1) % 8192];
            //
            //             let red = (colour.0 + ((colour2.0 - colour.0) * factor.fract())) as u8;
            //             let green = (colour.1 + ((colour2.1 - colour.1) * factor.fract())) as u8;
            //             let blue = (colour.2 + ((colour2.2 - colour.2) * factor.fract())) as u8;
            //
            //             image[3 * index] = red;
            //             image[3 * index + 1] = green;
            //             image[3 * index + 2] = blue;
            //         }
            //     }
            // }
        }
    }
}