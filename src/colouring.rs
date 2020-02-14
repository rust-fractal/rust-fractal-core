use crate::util::point::Point;

pub enum ColourMethod {
    Iteration,
    IterationSquareRoot,
    Histogram
}

impl ColourMethod {
    pub fn run(&self, points: &Vec<Point>, image: &mut Vec<u8>, maximum_iterations: usize, display_glitches: bool) {
        // Palette is temporarily here
        let mut colours = Vec::new();

        for i in 0..8192 {
            let value = i as f64 / 8192 as f64;

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
                for point in points {
                    let index = point.index;

                    if point.glitched && display_glitches {
                        image[3 * index] = 255u8;
                        image[3 * index + 1] = 0u8;
                        image[3 * index + 2] = 0u8;
                    } else if point.iterations.floor() >= maximum_iterations as f64 {
                        image[3 * index] = 0u8;
                        image[3 * index + 1] = 0u8;
                        image[3 * index + 2] = 0u8;
                    } else {
                        let hue = 200.0 * point.iterations % 8192.0;

                        let colour = colours[(hue.floor() as usize) % 8192];
                        let colour2 = colours[(hue.floor() as usize + 1) % 8192];

                        let red = (colour.0 + ((colour2.0 - colour.0) * hue.fract())) as u8;
                        let green = (colour.1 + ((colour2.1 - colour.1) * hue.fract())) as u8;
                        let blue = (colour.2 + ((colour2.2 - colour.2) * hue.fract())) as u8;

                        image[3 * index] = red;
                        image[3 * index + 1] = green;
                        image[3 * index + 2] = blue;
                    }
                }
            },
            ColourMethod::IterationSquareRoot => {
                for point in points {
                    let index = point.index;

                    if point.glitched && display_glitches {
                        image[3 * index] = 255u8;
                        image[3 * index + 1] = 0u8;
                        image[3 * index + 2] = 0u8;
                    } else if point.iterations.floor() >= maximum_iterations as f64 {
                        image[3 * index] = 0u8;
                        image[3 * index + 1] = 0u8;
                        image[3 * index + 2] = 0u8;
                    } else {
                        let hue = 1600.0 * point.iterations.sqrt() % 8192.0;

                        let colour = colours[(hue.floor() as usize) % 8192];
                        let colour2 = colours[(hue.floor() as usize + 1) % 8192];

                        let red = (colour.0 + ((colour2.0 - colour.0) * hue.fract())) as u8;
                        let green = (colour.1 + ((colour2.1 - colour.1) * hue.fract())) as u8;
                        let blue = (colour.2 + ((colour2.2 - colour.2) * hue.fract())) as u8;

                        image[3 * index] = red;
                        image[3 * index + 1] = green;
                        image[3 * index + 2] = blue;
                    }
                }
            },
            ColourMethod::Histogram => {
                let mut iteration_counts = vec![0usize; maximum_iterations + 1];

                for point in points {
                    iteration_counts[point.iterations.floor() as usize] += 1
                }

                for i in 1..iteration_counts.len() {
                    iteration_counts[i] += iteration_counts[i - 1];
                }

                let total = iteration_counts[maximum_iterations - 1];

                for point in points {
                    let index = point.index;

                    if point.glitched && display_glitches {
                        image[3 * index] = 255u8;
                        image[3 * index + 1] = 0u8;
                        image[3 * index + 2] = 0u8;
                    } else if point.iterations.floor() >= maximum_iterations as f64 {
                        image[3 * index] = 0u8;
                        image[3 * index + 1] = 0u8;
                        image[3 * index + 2] = 0u8;
                    } else {
                        let v1 = iteration_counts[point.iterations.floor() as usize] as f64 / total as f64;
                        let v2 = iteration_counts[point.iterations.floor() as usize + 1] as f64 / total as f64;

                        // the hue is used to smooth the histogram bins. The hue is in the range 0.0-1.0
                        let hue = (v1 + (v2 - v1) * point.iterations.fract() as f64) * 8192.0;

                        let colour = colours[(hue.floor() as usize) % 8192];
                        let colour2 = colours[(hue.floor() as usize + 1) % 8192];

                        let red = (colour.0 + ((colour2.0 - colour.0) * hue.fract())) as u8;
                        let green = (colour.1 + ((colour2.1 - colour.1) * hue.fract())) as u8;
                        let blue = (colour.2 + ((colour2.2 - colour.2) * hue.fract())) as u8;

                        image[3 * index] = red;
                        image[3 * index + 1] = green;
                        image[3 * index + 2] = blue;
                    }
                }
            }
        }
    }
}