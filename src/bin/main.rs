use std::time::Instant;
use rust_fractal::renderer::Renderer;

fn main() {
    let width: usize = 2000;
    let height: usize = 2000;

    let center = ("-0.802970043288051498602228364433675973121848792122673962954317", "0.178340211849727306925994231425911682945348875011614003515315");
    let zoom = 8.044158107639076e27;

    let mut renderer = Renderer::new(
        width,
        height,
        zoom,
        1000000,
        center.0,
        center.1,
        320
    );

    let time = Instant::now();
    renderer.render();
    println!("Rendering took: {} ms", time.elapsed().as_millis());

    // simple multithreading takes 13-15 sec f64
    // single threading takes 84 sec f64
}