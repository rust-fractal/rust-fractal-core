use std::time::Instant;
use rust_fractal::renderer::Renderer;

fn main() {
    // use f32 to zoom = 1x10-308, then f64 to 3.65x10-4951 possibly?

    let width: usize = 1000;
    let height: usize = 1000;

    let center = (
        "-0.4032037107646445210993574275561294",
        "-0.5943378885148266682959301608267696");
    let zoom = 2.1474836480000001E9;

    let mut renderer = Renderer::new(
        width,
        height,
        zoom,
        2500,
        center.0,
        center.1,
        250,
    );

    println!("Mandelbrot Renderer");

    let time = Instant::now();
    renderer.render();
    println!("{:<10}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}