use std::time::Instant;
use rust_fractal::renderer::Renderer;

fn main() {
    // f32 max is 1x10-32, use f64 to zoom = 1x10-308, then f128 to 3.65x10-4951 possibly?

    let center = (
        "-0.802970043288051498602228364433675973121848792122673962954317",
        "0.178340211849727306925994231425911682945348875011614003515315");
    let zoom = 8.044158107639076e27;

    let mut renderer = Renderer::new(
        2000,
        2000,
        zoom,
        1000000,
        center.0,
        center.1,
        320,
        0.2,
        true
    );

    println!("Mandelbrot Renderer");

    let time = Instant::now();
    renderer.render();
    println!("{:<10}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}