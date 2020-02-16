use crate::util::ComplexFixed;
use crate::renderer::Renderer;
use crate::util::point::Point;
use rayon::prelude::*;
use crate::util::complex_vector::ComplexVector;
use packed_simd::*;

impl Renderer {
    pub fn calculate_perturbations_f64_vectorised(&self, points_remaining: &mut Vec<Point>, points_complete: &mut Vec<Point>) {
        // here we pack the points into vectorised points
// possibly only works on image sizes that are multiples of the vector lanes
// investigate adding this earlier so that the vector of points never needs to be constructed

        type FloatVector = f64x4;
        type MaskVector = m64x4;
        let lanes = FloatVector::lanes();

//        let points_vectors = (0..remaining_points.len())
//            .step_by(lanes)
//            .map( | i| {
//                let delta_re = (0..lanes)
//                    .map( | j | {
//                        if i + j >= remaining_points.len() {
//                            100.0
//                        } else {
//                            remaining_points[i + j].delta.re - self.reference_delta.re
//                        }
//                    })
//                    .collect::< Vec < f64 > >();
//
//                let delta_im = (0..lanes)
//                    .map( | j | {
//                        if i + j >= remaining_points.len() {
//                            100.0
//                        } else {
//                            remaining_points[i + j].delta.im - self.reference_delta.im
//                        }
//                    })
//                    .collect::< Vec < f64 > >();
//
//                ComplexVector::< FloatVector >::new(delta_re.as_slice(), delta_im.as_slice())
//            })
//            .collect::<Vec<ComplexVector<FloatVector>>>();

        // TODO make sure divisible by 4
        let values = points_remaining.into_par_iter()
                                     .chunks(lanes)
                                     .map(|chunk| {
                                         let points_re = chunk.into_iter()
                                           .map(|point| {
                                               point.delta.re - self.reference_delta.re
                                           })
                                           .collect::<Vec<f64>>();

                                       let points_im = chunk.into_iter()
                                           .map(|point| {
                                               point.delta.re - self.reference_delta.re
                                           })
                                           .collect::<Vec<f64>>();

                                       let mut delta_0 = ComplexVector::<FloatVector>::new(points_re.as_slice(), points_im.as_slice());
                                       let mut delta_n = delta_0;
                                       let mut iterations = FloatVector::splat(0.0);
                                       let mut z_norm = FloatVector::splat(0.0);
                                       let mut glitched = MaskVector::splat(false);
                                       let mut escaped = MaskVector::splat(false);

                                       for iteration in 0..self.maximum_iterations {
                                           let temp = delta_n;
                                           delta_n += ComplexVector::<FloatVector >::splat( self.x_n_2[iteration]);
                                           delta_n *= temp;
                                           delta_n += delta_0;

                                           let new_z_norm = (ComplexVector::< FloatVector >::splat( self.x_n[iteration + 1]) + delta_n).norm_sqr();
                                           z_norm = escaped.select(z_norm, new_z_norm);
                                           escaped = z_norm.ge(FloatVector::splat(256.0));

                                           glitched = glitched.select(MaskVector::splat(true), z_norm.le(FloatVector::splat( self.tolerance_check[iteration + 1])));

                                           if (escaped | glitched).all() {
                                               break;
                                           }

                                           iterations += escaped.select(FloatVector::splat(0.0), FloatVector::splat(1.0));
                                       }

                                       let nu = FloatVector::splat(2.0) - (z_norm.ln() / (FloatVector::splat(2.0) * FloatVector::LN_2)).ln() / FloatVector::LN_2;

                                       // if a pixel has escaped it is fine to run the smooth colouring algorithm
                                       let out = iterations + escaped.select(nu, FloatVector::splat(0.0));

                                       let mut test = Vec::new();

                                       for i in 0..lanes {
                                           test.push((out.extract(i) as f64, glitched.extract(i)))
                                       }
                                       test
                                   })
                                   .flatten()
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

