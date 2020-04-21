use crate::math::perturbation::Perturbation;
use crate::util::PixelData;

pub struct Image {
    width: usize,
    height: usize,
    rgb: Vec<u8>,
    pub display_glitches: bool
}

impl Image {
    pub fn new(width: usize, height: usize, display_glitches: bool) -> Self {
        Image {
            width,
            height,
            rgb: vec![0u8; width * height * 3],
            display_glitches
        }
    }

    pub fn plot(&mut self, x: usize, y: usize, r: u8, g: u8, b: u8) {
        let k = (y * self.width + x) * 3;
        self.rgb[k] = r;
        self.rgb[k + 1] = g;
        self.rgb[k + 2] = b;
    }

    pub fn save(&mut self) {
        image::save_buffer("output.png", &self.rgb, self.width as u32, self.height as u32, image::RGB(8)).unwrap();
    }

    pub fn plot_image(&mut self, pixel_data: &Vec<PixelData>, delta_pixel: f64) {
        for packet in pixel_data {
            let (r, g, b) = if packet.glitched && self.display_glitches {
                (255, 0, 0)
            } else {
                if packet.escaped {
                    let de = 2.0 * packet.delta_current.norm() * packet.delta_current.norm().ln() / packet.derivative_current.norm();
                    let out = (255.0 * (de / delta_pixel).tanh()) as u8;
                    (out, out, out)
                } else {
                    (0, 0, 0)
                }
            };

            self.plot(packet.image_x, packet.image_y, r, g, b);
        }
    }
}