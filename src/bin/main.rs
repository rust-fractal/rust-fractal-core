use std::time::Instant;

use rust_fractal::renderer::Renderer;


fn main() {
    let width: usize = 640;
    let height: usize = 640;

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

    renderer.render();
}