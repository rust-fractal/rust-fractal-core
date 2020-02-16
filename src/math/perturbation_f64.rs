use crate::util::ComplexFixed;
use crate::renderer::Renderer;
use crate::util::point::Point;
use rayon::prelude::*;

impl Renderer {
    // This function will run the perturbation algorithm on all of the deltas in the locations,
    pub fn calculate_perturbations_f64(&self, points_remaining: &mut Vec<Point>, points_complete: &mut Vec<Point>) {
        *points_remaining = points_remaining.into_par_iter()
                                     .map(|point| {
                                         let delta_0 = point.delta - self.reference_delta;
                                         let mut delta_n = delta_0;
                                         let mut iteration = 0;
                                         let mut glitched = false;
                                         let mut z_norm = 0.0;

                                         while z_norm < 256.0 && iteration < self.maximum_iterations {
                                             let temp = delta_n;
                                             delta_n += self.x_n_2[iteration];
                                             delta_n *= temp;
                                             delta_n += delta_0;

                                             iteration += 1;
                                             z_norm = (self.x_n[iteration] + delta_n).norm_sqr();

                                             if z_norm < self.tolerance_check[iteration] {
                                                 glitched = true;
                                                 break;
                                             }
                                         }

                                         Point {
                                             delta: point.delta,
                                             index: point.index,
                                             iterations: iteration,
                                             smooth: 1.0 - (z_norm.log2() / 2.0).log2() as f32,
                                             glitched
                                         }
                                     })
                                     .collect::<Vec<Point>>();

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
