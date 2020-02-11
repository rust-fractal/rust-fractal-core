use std::time::Instant;
use rust_fractal::renderer::Renderer;

fn main() {
    // use f32 to zoom = 1x10-308, then f64 to 3.65x10-4951 possibly?

    let width: usize = 6000;
    let height: usize = 2000;

    let center = (
        "-0.8029700432880514986022283644336760",
        "0.1783402118497273069259942314259117");
    let zoom = 3E20;

    let mut renderer = Renderer::new(
        width,
        height,
        zoom,
        10000,
        center.0,
        center.1,
        250,
        1
    );

    println!("Mandelbrot Renderer");

    let time = Instant::now();
    renderer.render();
    println!("{:<10}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}