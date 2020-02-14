use std::time::Instant;
use rust_fractal::renderer::Renderer;

fn main() {
    // f32 max is 1x10-32, use f64 to zoom = 1x10-308, then f128 to 3.65x10-4951 possibly?

    let center = (
        "-0.6453593405513405298032904276942",
        "-0.4761981646219887966000672580610");
    let zoom = 2.0971520000000017E6;

    let mut renderer = Renderer::new(
        1600,
        1600,
        zoom,
        1000,
        center.0,
        center.1,
        250,
        0.001,
        false
    );

    println!("Mandelbrot Renderer");

    let time = Instant::now();
    renderer.render();
    println!("{:<10}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}