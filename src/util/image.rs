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
}