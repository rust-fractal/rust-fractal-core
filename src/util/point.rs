use crate::util::ComplexFixed;

#[derive(Copy, Clone)]
pub struct Point<T> {
    pub delta: T,
    pub index: usize,
    pub iterations: usize,
    pub smooth: f32,
    pub glitched: bool
}

impl Point<ComplexFixed<f32>> {
    // start_delta is the delta from reference of point (0, 0) in the image
    pub fn new(image_x: usize, image_y: usize, image_width: usize, resolution: f32, top_left_delta: ComplexFixed<f32>) -> Self {
        Point {
            delta: ComplexFixed::new(image_x as f32 * resolution + top_left_delta.re, image_y as f32 * resolution + top_left_delta.im),
            index: image_y * image_width + image_x,
            iterations: 0,
            smooth: 0.0,
            glitched: false
        }
    }
}

impl Point<ComplexFixed<f64>> {
    // start_delta is the delta from reference of point (0, 0) in the image
    pub fn new(image_x: usize, image_y: usize, image_width: usize, resolution: f64, top_left_delta: ComplexFixed<f64>) -> Self {
        Point {
            delta: ComplexFixed::new(image_x as f64 * resolution + top_left_delta.re, image_y as f64 * resolution + top_left_delta.im),
            index: image_y * image_width + image_x,
            iterations: 0,
            smooth: 0.0,
            glitched: false
        }
    }
}