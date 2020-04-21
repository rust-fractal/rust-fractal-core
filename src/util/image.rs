use crate::math::perturbation::Perturbation;

pub struct Image {
    width: usize,
    height: usize,
    rgb: Vec<u8>,
    display_glitches: bool
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

    pub fn plot_image(&mut self, perturbation: &Perturbation, delta_pixel: f64) {
        for i in 0..perturbation.image_x.len() {
            let (r, g, b) = if perturbation.glitched[i] && self.display_glitches {
                (255, 0, 0)
            } else {
                if perturbation.escaped[i] {
                    let de = 2.0 * perturbation.delta_current[i].norm() * perturbation.delta_current[i].norm().ln() / perturbation.derivative_current[i].norm();
                    let out = (255.0 * (de / delta_pixel).tanh()) as u8;
                    (out, out, out)
                } else {
                    (0, 0, 0)
                }
            };

            self.plot(perturbation.image_x[i], perturbation.image_y[i], r, g, b);
        }
    }
}