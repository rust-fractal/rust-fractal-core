use std::time::Instant;
use rayon::prelude::*;
use packed_simd::*;

fn main() {
    let width: usize = 1920;
    let height: usize = 1080;

    let x_range = (-2.5, 1.0);
    let y_range = (-1.0, 1.0);

    let x_res = (x_range.1 - x_range.0) / width as f64;
    let y_res = (y_range.1 - y_range.0) / height as f64;

    let max_iterations = 1024;

    let block_size = f64x8::lanes();
    let width_in_blocks = width / block_size;

    let time = Instant::now();
    let pixels = (0..height)
        .into_par_iter()
        .map(|y_pixel| {
            let xs = unsafe {
                let mut buf: Vec<f64x8> = vec![f64x8::splat(0.0); width_in_blocks];

                std::slice::from_raw_parts_mut(buf.as_mut_ptr() as *mut f64, width)
                    .iter_mut()
                    .enumerate()
                    .for_each(|(j, x)| {
                        *x = j as f64 * x_res + x_range.0
                    });
                buf
            };

            let y0 = f64x8::splat(y_pixel as f64 * y_res - y_range.1);

            let mut out = vec![u64x8::splat(0); width_in_blocks];

            out.iter_mut()
                .enumerate()
                .for_each(|(i, element)| {
                    // this is each chuck of the buffer
                    let x0 = xs[i];
                    let mut x = x0;
                    let mut y = y0;
                    let mut count = u64x8::splat(0);

                    for _ in 0..max_iterations {
                        let xx = x * x;
                        let yy = y * y;
                        let xy = f64x8::splat(2.0) * x * y;
                        let sum = xx + yy;
                        x = xx - yy + x0;
                        y = xy + y0;
                        let mask = sum.le(f64x8::splat(4.0));

                        if mask.none() {
                            break;
                        }

                        count += mask.select(u64x8::splat(1), u64x8::splat(0));
                    }
                    *element = count;
                });

            out
        })
        .flatten()
        .collect::<Vec<u64x8>>();

    println!("Rendering took {}ms", time.elapsed().as_millis());
    let time = Instant::now();



    let mut buffer = Vec::new();

    for element in pixels {
        let mut part = vec![0u64; 8];
        element.write_to_slice_unaligned(&mut part);
        for element in part {
            buffer.push((element as f64 / max_iterations as f64 * 255.99) as u8)
        }
    }

    image::save_buffer("output.png", &buffer, width as u32, height as u32, image::Gray(8)).unwrap();
    println!("Saving took {}ms", time.elapsed().as_millis());
}