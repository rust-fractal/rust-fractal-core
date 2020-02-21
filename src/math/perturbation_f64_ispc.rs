use crate::renderer::ImageRenderer;
use crate::util::point::Point;
use crate::util::ComplexFixed;

ispc_module!(mandelbrot);

impl ImageRenderer {
    pub fn calculate_perturbations_f64_ispc(&self, points_remaining: &mut Vec<Point<ComplexFixed<f64>>>, points_complete: &mut Vec<Point<ComplexFixed<f64>>>) {
        // 11438 ms 1000x1000 100000 1000

        let mut counts = vec![0; points_remaining.len()];
        let mut smooth = vec![0.0f32; points_remaining.len()];
        let mut glitched = vec![false; points_remaining.len()];

        let delta_real = points_remaining.into_iter().map(|point| {
            point.delta.re
        }).collect::<Vec<f64>>();

        let delta_imag = points_remaining.into_iter().map(|point| {
            point.delta.im
        }).collect::<Vec<f64>>();

        let reference_real = (&self.x_n_f64).into_iter().map(|value| {
            value.re
        }).collect::<Vec<f64>>();

        let reference_imag = (&self.x_n_f64).into_iter().map(|value| {
            value.im
        }).collect::<Vec<f64>>();

        unsafe {
            // about 1000 for the task size seems to work
            mandelbrot::mandelbrot_perturbation(1000,
                                                (points_remaining.len() / 1000) as i32,
                                                self.maximum_iterations as i32,
                                                delta_real.as_ptr(),
                                                delta_imag.as_ptr(),
                                                counts.as_mut_ptr(),
                                                smooth.as_mut_ptr(),
                                                glitched.as_mut_ptr(),
                                                reference_real.as_ptr(),
                                                reference_imag.as_ptr());
        }

        for i in 0..points_remaining.len() {
            points_remaining[i].iterations = counts[i] as usize;
            points_remaining[i].smooth = smooth[i];
            points_remaining[i].glitched = glitched[i];
        }

        for i in 0..points_remaining.len() {
            // check to see if a point is glitched
            if !glitched[i] {
                points_complete.push(points_remaining[i].clone());
            }
        }

        // Remove all non-glitched points from the remaining points
        points_remaining.retain(|point| {
            point.glitched
        });
    }
}