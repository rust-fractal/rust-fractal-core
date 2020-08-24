use exr::prelude::simple_image;

// TODO possibly have option to change between EXR and PNG saving
pub struct Image {
    width: usize,
    height: usize,
    rgb: Vec<u8>,
    iterations: Vec<u32>,
    smooth: Vec<f32>,
    pub display_glitches: bool
}

impl Image {
    pub fn new(width: usize, height: usize, display_glitches: bool) -> Self {
        Image {
            width,
            height,
            rgb: vec![0u8; width * height * 3],
            iterations: vec![0u32; width * height],
            smooth: vec![0.0f32; width * height],
            display_glitches
        }
    }

    pub fn plot_jpg(&mut self, x: usize, y: usize, r: u8, g: u8, b: u8) {
        let k = (y * self.width + x) * 3;
        self.rgb[k] = r;
        self.rgb[k + 1] = g;
        self.rgb[k + 2] = b;
    }

    pub fn plot_exr(&mut self, x: usize, y: usize, iterations: u32, smooth: f32) {
        let k = y * self.width + x;
        self.iterations[k] = iterations;
        self.smooth[k] = smooth;
    }

    pub fn save_jpg(&mut self, filename: &String) {
        image::save_buffer(filename, &self.rgb, self.width as u32, self.height as u32, image::ColorType::Rgb8).unwrap();
    }

    pub fn save_exr(&mut self, filename: &String) {
        let iterations = simple_image::Channel::non_color_data(simple_image::Text::from("N").unwrap(), simple_image::Samples::U32(self.iterations.clone()));
        let smooth = simple_image::Channel::non_color_data(simple_image::Text::from("NF").unwrap(), simple_image::Samples::F32(self.smooth.clone()));

        let layer = simple_image::Layer::new(simple_image::Text::from("fractal_data").unwrap(), (self.width, self.height), smallvec::smallvec![iterations, smooth])
            .with_compression(simple_image::Compression::PXR24)
            .with_block_format(None, simple_image::attribute::LineOrder::Increasing);

        let image = simple_image::Image::new_from_single_layer(layer);

        image.write_to_file(filename, simple_image::write_options::high()).unwrap();
    }
}