use std::time::Instant;
use rust_fractal::renderer::FractalRenderer;
use std::fs;
use std::io::{stdout, stdin, Write};

fn main() {
    println!("Mandelbrot Renderer");

    // TODO on some images, we may need to check the imaginary part of the floatexp to make sure that is not too large.

    let mut s = String::new();
    print!("File to render: ");
    let _ = stdout().flush();
    stdin().read_line(&mut s).expect("User did not enter a correct string.");

    s.pop();

    if !s.ends_with(".kfr") {
        s.push_str(".kfr")
    };

    println!("Reading: {}", s);

    let data = fs::read_to_string(s).expect("Unable to read file.");

    let mut center_re = "-0.75";
    let mut center_im  = "0.0";
    let mut zoom = "1E0";
    let mut iterations = "1000";

    for line in data.lines() {
        let mut parts = line.split_whitespace();
        // println!("{}", line);
        let temp = parts.next().unwrap();

        match temp {
            "Re:" => {
                center_re = parts.next().unwrap();
            },
            "Im:" => {
                center_im = parts.next().unwrap();
            }
            "Zoom:" => {
                zoom = parts.next().unwrap();
            },
            "Iterations:" => {
                iterations = parts.next().unwrap();
            }
            _ => {}
        }
    }

    println!("Zoom: {}", zoom);

    let mut renderer = FractalRenderer::new(
        1000,
        1000,
        zoom,
        iterations.parse::<usize>().unwrap(),
        center_re,
        center_im,
        0.01,
        false,
        32
    );

    let time = Instant::now();
    renderer.render();
    println!("{:<14}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}