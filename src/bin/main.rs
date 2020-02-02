use std::time::Instant;
use rayon::prelude::*;
use packed_simd::*;
use std::f32::consts::LN_2;

fn main() {
    let width: usize = 1920 * 16;
    let height: usize = 1080 * 16;

    let aspect = height as f32 / width as f32;

    let center = (-0.743314, 0.131612);
    let zoom = 80.0;

    let x_range = (center.0 - 0.5 / zoom, center.0 + 0.5 / zoom);
    let y_range = (center.1 - 0.5 / zoom * aspect, center.1 + 0.5 / zoom * aspect);

    let x_res = (x_range.1 - x_range.0) / width as f32;
    let y_res = (y_range.1 - y_range.0) / height as f32;

    let max_iterations = 8192;

    let block_size = f32x16::lanes();
    let width_in_blocks = width / block_size;

    let mut iterations = vec![0.0f32; width * height];

    let time = Instant::now();
    let pixels = (0..height)
        .into_par_iter()
        .map(|y_pixel| {
            let xs = unsafe {
                let mut buf: Vec<f32x16> = vec![f32x16::splat(0.0); width_in_blocks];

                std::slice::from_raw_parts_mut(buf.as_mut_ptr() as *mut f32, width)
                    .iter_mut()
                    .enumerate()
                    .for_each(|(j, x)| {
                        *x = j as f32 * x_res + x_range.0
                    });
                buf
            };

            let y0 = f32x16::splat(y_pixel as f32 * y_res - y_range.1);

            let mut out = vec![f32x16::splat(0.0); width_in_blocks];

            out.iter_mut()
                .enumerate()
                .for_each(|(i, element)| {
                    // this is each chuck of the buffer
                    let x0 = xs[i];
                    let mut count = f32x16::splat(0.0);
                    let mut sum = f32x16::splat(0.0);
                    let mut mask = m32x16::splat(true);

                    let mut x = x0;
                    let mut y = y0;

                    for _ in 0..max_iterations {
                        let xx = x * x;
                        let yy = y * y;
                        let xy = f32x16::splat(2.0) * x * y;
                        sum = mask.select(xx + yy, sum);
                        x = xx - yy + x0;
                        y = xy + y0;

                        mask = sum.le(f32x16::splat(10000000.0));

                        if mask.none() {
                            break;
                        }

                        count += mask.select(f32x16::splat(1.0), f32x16::splat(0.0));
                    }

                    let nu = f32x16::splat(1.0) - (sum.ln() / f32x16::splat(2.0 * LN_2)).ln() / f32x16::splat(LN_2);

                    *element = count + mask.select(f32x16::splat(0.0), nu);
                });

            out
        })
        .flatten()
        .collect::<Vec<f32x16>>();

    for (i, element) in pixels.into_iter().enumerate() {
        element.write_to_slice_unaligned(&mut iterations[(16 * i)..(16 * (i + 1))]);
    }

    println!("Calculating: {}ms", time.elapsed().as_millis());
    let time = Instant::now();

    // 2nd pass
    let mut iteration_counts = vec![0usize; max_iterations + 1];

    for iteration in &iterations {
        iteration_counts[iteration.floor() as usize] += 1
    }

    // 3rd pass

    // modify here to convert to cumulative
    for i in 1..iteration_counts.len() {
        iteration_counts[i] += iteration_counts[i - 1];
    }

    let mut total = iteration_counts[max_iterations - 1];

    let mut colors = vec![0u8; width * height * 3];

    for i in 0..(width * height) {
        if iterations[i].floor() >= max_iterations as f32 {
            colors[3 * i] = 0u8;
            colors[3 * i + 1] = 0u8;
            colors[3 * i + 2] = 0u8;
        } else {
            let test = iteration_counts[iterations[i].floor() as usize] as f64 / total as f64;
            let test2 = iteration_counts[iterations[i].floor() as usize + 1] as f64 / total as f64;

            let hue = test + (test2 - test) * iterations[i].fract() as f64;

            let test = (hue * 255.99) as u8;
            colors[3 * i] = test;
            colors[3 * i + 1] = test;
            colors[3 * i + 2] = test;
        }
    }

    println!("Colouring: {}ms", time.elapsed().as_millis());
    let time = Instant::now();

    image::save_buffer("output.png", &colors, width as u32, height as u32, image::RGB(8)).unwrap();
    println!("Saving: {}ms", time.elapsed().as_millis());
}