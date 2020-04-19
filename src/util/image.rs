pub struct Image {
    width: i32,
    height: i32,
    rgb: Vec<u8>
}

impl Image {
    pub fn new(width: i32, height: i32) -> Self {
        Image {
            width,
            height,
            rgb: vec![0u8; (width * height * 3) as usize]
        }
    }

    pub fn plot(&mut self, x: i32, y: i32, r: u8, g: u8, b: u8) {
        let k = ((y * self.width + x) * 3) as usize;
        self.rgb[k] = r;
        self.rgb[k + 1] = g;
        self.rgb[k + 2] = b;
    }

    pub fn save(&mut self) {
        image::save_buffer("output.png", &self.rgb, self.width as u32, self.height as u32, image::RGB(8)).unwrap();
    }
}