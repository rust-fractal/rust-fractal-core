use std::time::Instant;
use rayon::prelude::*;
use packed_simd::*;

fn main() {
    let width: usize = 1920;
    let height: usize = 1080;

    let aspect = height as f32 / width as f32;

    let center = (-0.745314, 0.131612);
    let zoom = 80.0;

    let x_range = (center.0 - 0.5 / zoom, center.0 + 0.5 / zoom);
    let y_range = (center.1 - 0.5 / zoom * aspect, center.1 + 0.5 / zoom * aspect);

    let x_res = (x_range.1 - x_range.0) / width as f32;
    let y_res = (y_range.1 - y_range.0) / height as f32;

    let max_iterations = 1024;

    let block_size = f32x16::lanes();
    let width_in_blocks = width / block_size;

    let mut iterations = vec![0u32; width * height];

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

            let mut out = vec![u32x16::splat(0); width_in_blocks];

            out.iter_mut()
                .enumerate()
                .for_each(|(i, element)| {
                    // this is each chuck of the buffer
                    let x0 = xs[i];
                    let mut x = x0;
                    let mut y = y0;
                    let mut count = u32x16::splat(0);

                    for _ in 0..max_iterations {
                        let xx = x * x;
                        let yy = y * y;
                        let xy = f32x16::splat(2.0) * x * y;
                        let sum = xx + yy;
                        x = xx - yy + x0;
                        y = xy + y0;
                        let mask = sum.le(f32x16::splat(4.0));

                        if mask.none() {
                            break;
                        }

                        count += mask.select(u32x16::splat(1), u32x16::splat(0));
                    }
                    *element = count;
                });

            out
        })
        .flatten()
        .collect::<Vec<u32x16>>();

    for (i, element) in pixels.into_iter().enumerate() {
        element.write_to_slice_unaligned(&mut iterations[(16 * i)..(16 * (i + 1))]);
    }

    let colors = iterations
        .into_iter()
        .map(|item| {
            ((item as f32 / max_iterations as f32) * 255.99) as u8
        })
        .collect::<Vec<u8>>();

    println!("Rendering: {}ms", time.elapsed().as_millis());
    let time = Instant::now();
    image::save_buffer("output.png", &colors, width as u32, height as u32, image::Gray(8)).unwrap();
    println!("Saving: {}ms", time.elapsed().as_millis());
}