use crate::renderer::ImageRenderer;
use crate::util::point::Point;
use crate::util::ComplexFixed;

ispc_module!(mandelbrot);

impl ImageRenderer {
    pub fn calculate_perturbations_f64_ispc(&self, points_remaining: &mut Vec<Point<ComplexFixed<f64>>>, points_complete: &mut Vec<Point<ComplexFixed<f64>>>) {
        // 11438 ms 1000x1000 100000 1000
        // 17053 ms 1000x1000 100000 1000
        // TODO make it so that non-filled tasks are filled with dummy values (e.g. ones which will instantly escape)

        let mut counts = vec![0; points_remaining.len()];
        let mut smooth = vec![0.0f32; points_remaining.len()];
        let mut glitched = vec![0; points_remaining.len()];

        let delta_real = points_remaining.into_iter().map(|point| {
            point.delta.re - self.reference_delta_f64.re
        }).collect::<Vec<f64>>();

        let delta_imag = points_remaining.into_iter().map(|point| {
            point.delta.im - self.reference_delta_f64.im
        }).collect::<Vec<f64>>();

        let reference_real = (&self.x_n_f64).into_iter().map(|value| {
            value.re
        }).collect::<Vec<f64>>();

        let reference_imag = (&self.x_n_f64).into_iter().map(|value| {
            value.im
        }).collect::<Vec<f64>>();

        unsafe {
            // about 1000 for the task size seems to work
            mandelbrot::mandelbrot_perturbation_f64(64,
                (points_remaining.len() / 64) as i32,
                (points_remaining.len() % 64) as i32,
                self.maximum_iterations as i32,
                delta_real.as_ptr(),
                delta_imag.as_ptr(),
                counts.as_mut_ptr(),
                smooth.as_mut_ptr(),
                glitched.as_mut_ptr(),
                reference_real.as_ptr(),
                reference_imag.as_ptr(),
                self.tolerance_check_f64.as_ptr())
        }

        for i in 0..points_remaining.len() {
            points_remaining[i].iterations = counts[i] as usize;
            points_remaining[i].smooth = smooth[i];
            points_remaining[i].glitched = glitched[i] != 0;
        }

        for i in 0..points_remaining.len() {
            // check to see if a point is glitched
            if !points_remaining[i].glitched {
                points_complete.push(points_remaining[i].clone());
            }
        }

        // Remove all non-glitched points from the remaining points
        points_remaining.retain(|point| {
            point.glitched
        });
    }
}