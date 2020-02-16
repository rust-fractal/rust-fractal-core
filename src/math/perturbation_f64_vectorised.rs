use crate::renderer::ImageRenderer;
use crate::util::point::Point;
use rayon::prelude::*;
use crate::util::complex_vector::ComplexVector;
use packed_simd::*;
use crate::util::ComplexFixed;

impl ImageRenderer {
    pub fn calculate_perturbations_f64_vectorised(&self, points_remaining: &mut Vec<Point<ComplexFixed<f64>>>, points_complete: &mut Vec<Point<ComplexFixed<f64>>>) {
        // here we pack the points into vectorised points
        // possibly only works on image sizes that are multiples of the vector lanes
        // investigate adding this earlier so that the vector of points never needs to be constructed

        // approximate speedup is 1.4x
        let lanes = f64x4::lanes();

        *points_remaining = points_remaining.chunks(lanes)
                                     .map(|chunk| {
                                         Vec::from(chunk)
                                     })
                                     .collect::<Vec<Vec<Point<ComplexFixed<f64>>>>>()
                                     .into_par_iter()
                                     .map(|chunk| {
                                         let mut points_re = vec![100.0; 4];
                                         let mut points_im = vec![100.0; 4];

                                         for i in 0..chunk.len() {
                                             points_re[i] = chunk[i].delta.re - self.reference_delta_f64.re;
                                             points_im[i] = chunk[i].delta.im - self.reference_delta_f64.im;
                                         }

                                         let delta_0 = ComplexVector::<f64x4>::new(points_re.as_slice(), points_im.as_slice());
                                         let mut delta_n = delta_0;
                                         let mut iterations = u64x4::splat(0);
                                         let mut z_norm = f64x4::splat(0.0);
                                         let mut glitched = m64x4::splat(false);
                                         let mut escaped = m64x4::splat(false);

                                         for iteration in 0..self.maximum_iterations {
                                             let temp = delta_n;
                                             delta_n += ComplexVector::<f64x4 >::splat( self.x_n_2_f64[iteration]);
                                             delta_n *= temp;
                                             delta_n += delta_0;

                                             let new_z_norm = (ComplexVector::< f64x4 >::splat( self.x_n_f64[iteration + 1]) + delta_n).norm_sqr();
                                             z_norm = escaped.select(z_norm, new_z_norm);
                                             escaped = z_norm.ge(f64x4::splat(256.0));

                                             glitched = glitched.select(m64x4::splat(true), z_norm.le(f64x4::splat( self.tolerance_check_f64[iteration + 1])));

                                             if (escaped | glitched).all() {
                                                 break;
                                             }

                                             iterations += escaped.select(u64x4::splat(0), u64x4::splat(1));
                                         }

                                         let nu = f64x4::splat(1.0) - (z_norm.ln() / (f64x4::splat(2.0) * f64x4::LN_2)).ln() / f64x4::LN_2;

                                         // if a pixel has escaped it is fine to run the smooth colouring algorithm
                                         let smooth = escaped.select(nu, f64x4::splat(0.0));

                                         let mut test = Vec::new();

                                         for i in 0..chunk.len() {
                                             test.push(Point {
                                                 delta: chunk[i].delta,
                                                 index: chunk[i].index,
                                                 iterations: iterations.extract(i) as usize,
                                                 smooth: smooth.extract(i) as f32,
                                                 glitched: glitched.extract(i)
                                             });
                                         }
                                         test
                                     })
                                     .flatten()
                                     .collect::<Vec<Point<ComplexFixed<f64>>>>();

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

