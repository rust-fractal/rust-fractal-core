use crate::math::reference::Reference;
use crate::util::ComplexFixed;
use crate::util::complex_vector::ComplexVector;

use rayon::prelude::*;
use packed_simd::*;

pub struct Perturbation {
    pub image_x: Vec<usize>,
    pub image_y: Vec<usize>,
    pub iteration: Vec<usize>,
    pub delta_reference: Vec<ComplexFixed<f64>>,
    pub delta_current: Vec<ComplexFixed<f64>>,
    pub derivative_current: Vec<ComplexFixed<f64>>,
    pub glitched: Vec<bool>,
    pub escaped: Vec<bool>,
}

// AoSoA layout
struct PerturbationData {
    iteration: u64x4,
    delta_reference: ComplexVector<f64x4>,
    delta_current: ComplexVector<f64x4>,
    derivative_current: ComplexVector<f64x4>,
    glitched: m64x4,
    escaped: m64x4,
}

impl Perturbation {
    // This sets up a new perturbation method. It will pack the points that are provided into the correct vectors
    pub fn new(image_x: Vec<usize>,
               image_y: Vec<usize>,
               iteration: Vec<usize>,
               delta_reference: Vec<ComplexFixed<f64>>,
               delta_current: Vec<ComplexFixed<f64>>,
               derivative_current: Vec<ComplexFixed<f64>>) -> Self {
        let length = image_x.len();
        Perturbation {
            image_x,
            image_y,
            iteration,
            delta_reference,
            delta_current,
            derivative_current,
            glitched: vec![false; length],
            escaped: vec![false; length],
        }
    }

    pub fn iterate(&mut self, reference: &Reference, maximum_iteration: usize) {
        // let mut packets = Vec::new();
        //
        // // TODO non 4x image
        // for i in (0..image_x.len()).step_by(4) {
        //     packets.push(
        //         PerturbationData {
        //             iteration: u64x4::splat(self.iteration[0]),
        //             delta_reference: ComplexVector::new2(&self.delta_reference[i..(i + 4)].as_slice()),
        //             delta_current: ComplexVector::new2(&self.delta_current[i..(i + 4)].as_slice()),
        //             derivative_current: ComplexVector::new2(&self.derivative_current[i..(i + 4)].as_slice()),
        //             glitched: m64x4::splat(false),
        //             escaped: m64x4::splat(false)
        //         }
        //     );
        // }
        //
        // *packets = packets.into_par_iter()
        //     .map(|packet| {
        //         // temporary - the number of iterations should be the same in all cases after the series approximation
        //         for iteration in self.iteration[0]..maximum_iteration {
        //             let temp = packet.delta_current;
        //             packet.delta_current +=
        //
        //
        //         }
        //
        //
        //
        //     }).collect::<Vec<PerturbationData>>();


        // TODO convert this into vectorised code
        for i in 0..self.image_x.len() {
            while self.iteration[i] < maximum_iteration {
                // This uses the difference between the starting iteration of the reference - can be used to skip some
                let z = self.delta_current[i] + reference.z_reference[(self.iteration[i] - reference.start_iteration) as usize];
                let z_norm = z.norm_sqr();

                if z_norm < reference.z_tolerance[(self.iteration[i] - reference.start_iteration) as usize] {
                    self.glitched[i] = true;
                    self.delta_current[i] = z;
                    break;
                }

                if z_norm > 65536.0 {
                    self.escaped[i] = true;
                    self.delta_current[i] = z;
                    break;
                }

                self.derivative_current[i] = 2.0 * z * self.derivative_current[i] + 1.0;
                self.delta_current[i] = 2.0 * reference.z_reference[(self.iteration[i] - reference.start_iteration) as usize] * self.delta_current[i] + self.delta_current[i] * self.delta_current[i] + self.delta_reference[i];
                self.iteration[i] += 1;
            }
        };
    }
}

