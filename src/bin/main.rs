use std::time::Instant;
use rust_fractal::renderer::Renderer;

fn main() {
    // f32 max is 1x10-32, use f64 to zoom = 1x10-308, then f128 to 3.65x10-4951 possibly?

    let center = (
        "-1.99996619445037030418434688506350579675531241540724851511761922944801584242342684381376129778868913812287046406560949864353810575744772166485672496092803920095332",
        "+0.00000000000000000000000000000000030013824367909383240724973039775924987346831190773335270174257280120474975614823581185647299288414075519224186504978181625478529");
    let zoom = 2.3620330788506154104770818136626E157;

//    let center = (
//        "0.0",
//        "0.0");
//    let zoom = 1.0;

    let mut renderer = Renderer::new(
        1000,
        400,
        zoom,
        100000,
        center.0,
        center.1,
        1000,
        0.000,
        false
    );

    println!("Mandelbrot Renderer");

    let time = Instant::now();
    renderer.render();
    println!("{:<10}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}