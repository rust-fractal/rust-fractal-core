use exr::prelude::simple_image::*;
use rayon::prelude::*;
use config::Config;

use std::fs;
use std::time::Instant;

use crate::util::generate_default_palette;

pub struct RecolourExr {
    palette: Vec<(u8, u8, u8)>,
    files: Vec<String>,
    iteration_division: f32
}

impl RecolourExr {
    pub fn new(settings: Config) -> Self {
        let palette = if let Ok(values) = settings.get_array("palette") {
            let palette = values.chunks_exact(3).map(|value| {
                // We assume the palette is in BGR rather than RGB
                (value[2].clone().into_int().unwrap() as u8, 
                    value[1].clone().into_int().unwrap() as u8, 
                    value[0].clone().into_int().unwrap() as u8)
            }).collect::<Vec<(u8, u8, u8)>>();

            palette
        } else {
            generate_default_palette()
        };
        let iteration_division = settings.get_float("iteration_division").unwrap_or(0.1) as f32;

        let paths = fs::read_dir("./output/").unwrap();
        let mut exr_files = Vec::new();
    
        for path in paths {
            let name = path.unwrap()
                .path()
                .to_owned()
                .to_str()
                .unwrap()
                .to_string();
            
            if name.contains(".exr") {
                exr_files.push(name)
            }
        };

        RecolourExr {
            palette,
            files: exr_files,
            iteration_division
        }
    }

    pub fn colour(&self) {
        let colouring_time = Instant::now();

        (&self.files).into_par_iter()
        .for_each(|exr_file| {
            let raw_data = Image::read_from_file(&exr_file, read_options::high()).unwrap();

            let mut iterations = Vec::new();
            let mut smooth = Vec::new();

            for layer in &raw_data.layers {
                for channel in &layer.channels {
                    match &channel.samples {
                        Samples::F16(f16_vec) => {
                            smooth = f16_vec.clone();
                        },
                        Samples::F32(_) => {},
                        Samples::U32(u32_vec) => {
                            iterations = u32_vec.clone();
                        }
                    };
                }
            }

            let file_name = exr_file.split(".exr").collect::<Vec<_>>()[0];
            let dimensions = raw_data.attributes.display_window.size;

            println!("{}.exr", file_name);

            let mut rgb_buffer = vec![0u8; iterations.len() * 3];
            
            for i in 0..iterations.len() {
                if iterations[i] == 0xFFFFFFFF {
                    rgb_buffer[3 * i] = 0u8;
                    rgb_buffer[3 * i + 1] = 0u8;
                    rgb_buffer[3 * i + 2] = 0u8;
                } else {
                    let temp = (iterations[i] as f32 + smooth[i].to_f32()) / self.iteration_division;

                    let temp2 = temp.floor() as usize % self.palette.len();
                    let temp3 = (temp as usize + 1) % self.palette.len();
                    let temp4 = temp.fract();

                    let colour1 = self.palette[temp2];
                    let colour2 = self.palette[temp3];

                    rgb_buffer[3 * i] = (colour1.0 as f32 + temp4 * (colour2.0 as f32 - colour1.0 as f32)) as u8; 
                    rgb_buffer[3 * i + 1] = (colour1.1 as f32 + temp4 * (colour2.1 as f32 - colour1.1 as f32)) as u8; 
                    rgb_buffer[3 * i + 2] = (colour1.2 as f32 + temp4 * (colour2.2 as f32 - colour1.2 as f32)) as u8;
                }
            }

            image::save_buffer(file_name.to_owned() + ".png", &rgb_buffer, dimensions.x() as u32, dimensions.y() as u32, image::ColorType::Rgb8).unwrap();
        });

        println!("Recolouring {} images took {} ms.", self.files.len(), colouring_time.elapsed().as_millis());
    }
}
