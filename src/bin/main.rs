use std::time::Instant;
use rust_fractal::renderer::FractalRenderer;
use float_extended::float_extended::FloatExtended;
use std::{io, fs};
use std::io::{stdout, stdin, Write, BufReader, BufRead};
use std::fs::File;

fn main() {
    println!("Mandelbrot Renderer");

    let mut s = String::new();
    print!("File to render: ");
    let _ = stdout().flush();
    stdin().read_line(&mut s).expect("User did not enter a correct string.");

    s.pop();

    if !s.ends_with(".kfr") {
        s.push_str(".kfr")
    };

    println!("Rendering {}...", s);

    let data = fs::read_to_string(s).expect("Unable to read file.");

    let mut center_re = "-0.75";
    let mut center_im  = "0.0";
    let mut zoom = "1E0";
    let mut iterations = "1000";

    for line in data.lines() {
        let mut parts = line.split_whitespace();
        // println!("{}", line);
        let mut temp = parts.next().unwrap();

        match temp {
            "Re:" => {
                let mut temp2 = parts.next().unwrap();
                center_re = temp2;
            },
            "Im:" => {
                let mut temp2 = parts.next().unwrap();
                center_im = temp2;
            }
            "Zoom:" => {
                let mut temp2 = parts.next().unwrap();
                zoom = temp2;
            },
            "Iterations:" => {
                let mut temp2 = parts.next().unwrap();
                iterations = temp2;
            }
            _ => {}
        }
    }

    let mut renderer = FractalRenderer::new(
        1000,
        1000,
        zoom,
        iterations.parse::<usize>().unwrap(),
        center_re,
        center_im,
        0.001,
        true,
        16
    );

    let time = Instant::now();
    renderer.render();
    println!("{:<14}{:>6} ms", "TOTAL", time.elapsed().as_millis());
}