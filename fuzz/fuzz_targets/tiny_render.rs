#![no_main]

use fractal_core::{
    ComplexCoordinate, FractalFormula, RenderRequest, RenderSettings, Viewport, Zoom,
    compute_frame,
};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if data.len() < 12 {
        return;
    }
    let width = u32::from(data[0]) % 16 + 1;
    let height = u32::from(data[1]) % 16 + 1;
    let iterations = u32::from(data[2]) % 128 + 1;
    let power = u32::from(data[3]) % 15 + 2;
    let zoom = f64::from(u16::from_le_bytes([data[4], data[5]])) / 128.0 + 0.01;
    let real = f64::from(i16::from_le_bytes([data[6], data[7]])) / 4096.0;
    let imag = f64::from(i16::from_le_bytes([data[8], data[9]])) / 4096.0;
    let rotation = f64::from(data[10]) / 16.0;
    let request = RenderRequest {
        generation: 1,
        viewport: Viewport {
            center: ComplexCoordinate::new(real, imag),
            zoom: Zoom::new(zoom).unwrap(),
            rotation_radians: rotation,
            width,
            height,
        },
        settings: RenderSettings {
            fractal: if data[11] & 1 == 0 {
                FractalFormula::Mandelbrot { power }
            } else {
                FractalFormula::BurningShip { power }
            },
            iteration_limit: iterations,
            ..RenderSettings::default()
        },
    };
    let _ = compute_frame(&request);
});
