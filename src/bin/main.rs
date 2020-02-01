use std::time::Instant;
use rayon::prelude::*;

fn main() {
    let width: usize = 1920;
    let height: usize = 1080;

    let x_range = (-2.5f32, 1.0f32);
    let y_range = (-1.0, 1.0);

    let x_res = (x_range.1 - x_range.0) / width as f32;
    let y_res = (y_range.1 - y_range.0) / height as f32;

    let max_iterations = 1024;


    let time = Instant::now();
    let pixels = (0..height)
        .into_par_iter()
        .map(|y_pixel| {
            let mut iterations_buffer = Vec::with_capacity(width);

            for x_pixel in 0..width {
                let x0 = x_pixel as f32 * x_res + x_range.0;
                let y0 = y_pixel as f32 * y_res - y_range.1;
                let mut x = x0;
                let mut y = y0;
                let mut xx = 0.0;
                let mut yy = 0.0;
                let mut xy = 0.0;
                let mut ab = 0.0;
                let mut iterations = 0;

                while ab < 4.0 && iterations < max_iterations {
                    xx = x * x;
                    yy = y * y;
                    xy = 2.0 * x * y;
                    iterations += 1;
                    ab = xx + yy;
                    x = xx - yy + x0;
                    y = xy + y0;
                }
                iterations_buffer.push((iterations as f32 / max_iterations as f32 * 255.99) as u8);
            }

            iterations_buffer
        })
        .flatten()
        .collect::<Vec<u8>>();

    println!("Rendering took {}ms", time.elapsed().as_millis());

    image::save_buffer("output.png", &pixels, width as u32, height as u32, image::Gray(8)).unwrap();
}