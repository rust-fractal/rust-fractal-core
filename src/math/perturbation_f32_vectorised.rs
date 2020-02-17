use crate::renderer::ImageRenderer;
use crate::util::point::Point;
use rayon::prelude::*;
use crate::util::complex_vector::ComplexVector;
use packed_simd::*;
use crate::util::ComplexFixed;

impl ImageRenderer {
    pub fn calculate_perturbations_f32_vectorised(&self, points_remaining: &mut Vec<Point<ComplexFixed<f32>>>, points_complete: &mut Vec<Point<ComplexFixed<f32>>>) {
        // here we pack the points into vectorised points
        // possibly only works on image sizes that are multiples of the vector lanes
        // investigate adding this earlier so that the vector of points never needs to be constructed

        // approximate speedup is 1.4x
        let lanes = f32x8::lanes();

        *points_remaining = points_remaining.par_chunks(lanes)
            .map(|chunk| {
                let mut points_re = vec![100.0; lanes];
                let mut points_im = vec![100.0; lanes];

                for i in 0..chunk.len() {
                    points_re[i] = chunk[i].delta.re - self.reference_delta_f32.re;
                    points_im[i] = chunk[i].delta.im - self.reference_delta_f32.im;
                }

                let delta_0 = ComplexVector::<f32x8>::new(points_re.as_slice(), points_im.as_slice());
                let mut delta_n = delta_0;
                let mut iterations = u32x8::splat(0);
                let mut z_norm = f32x8::splat(0.0);
                let mut glitched = m32x8::splat(false);
                let mut escaped = m32x8::splat(false);

                for iteration in 0..self.maximum_iterations {
                    let temp = delta_n;
                    delta_n += ComplexVector::<f32x8>::splat( self.x_n_2_f32[iteration]);
                    delta_n *= temp;
                    delta_n += delta_0;

                    let temp2 = ComplexVector::<f32x8>::splat( self.x_n_f32[iteration + 1]) + delta_n;
                    let new_z_norm = temp2.norm_sqr();
                    z_norm = escaped.select(z_norm, new_z_norm);
                    escaped = z_norm.ge(f32x8::splat(256.0));


                    glitched = glitched.select(m32x8::splat(true), z_norm.lt(f32x8::splat(self.tolerance_check_f32[iteration + 1])));

                    escaped |= glitched;

                    if escaped.all() {
                        break;
                    }

                    iterations += escaped.select(u32x8::splat(0), u32x8::splat(1));
                }

                let nu = f32x8::splat(1.0) - (z_norm.ln() / (f32x8::splat(2.0) * f32x8::LN_2)).ln() / f32x8::LN_2;

                let mut test = Vec::new();

                for i in 0..chunk.len() {
                    test.push(Point {
                        delta: chunk[i].delta,
                        index: chunk[i].index,
                        iterations: iterations.extract(i) as usize,
                        smooth: nu.extract(i),
                        glitched: glitched.extract(i)
                    });
                }
                test
            })
            .flatten()
            .collect::<Vec<Point<ComplexFixed<f32>>>>();

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

