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
        let mut packets = Vec::new();

        // TODO non 4x image
        // for i in (0..self.image_x.len()).step_by(4) {
        //     packets.push(
        //         PerturbationData {
        //             iteration: u64x4::splat(self.iteration[0] as u64),
        //             delta_reference: ComplexVector::new2(&[self.delta_reference[i], self.delta_reference[i + 1], self.delta_reference[i + 2], self.delta_reference[i + 3]]),
        //             delta_current: ComplexVector::new2(&[self.delta_current[i], self.delta_current[i + 1], self.delta_current[i + 2], self.delta_current[i + 3]]),
        //             derivative_current: ComplexVector::new2(&[self.derivative_current[i], self.derivative_current[i + 1], self.derivative_current[i + 2], self.derivative_current[i + 3]]),
        //             glitched: m64x4::splat(false),
        //             escaped: m64x4::splat(false)
        //         }
        //     );
        // }

        // packets.par_chunks_mut(1)
        //     .for_each(|packet| {
        //         // temporary - the number of iterations should be the same in all cases after the series approximation
        //         for iteration in self.iteration[0]..maximum_iteration {
        //             let z = packet[0].delta_current + ComplexVector::<f64x4>::splat(reference.z_reference[(iteration - reference.start_iteration)]);
        //             let z_norm = z.norm_sqr();
        //
        //             packet[0].escaped = packet[0].escaped.select(m64x4::splat(true), z_norm.ge(f64x4::splat(65536.0)));
        //             packet[0].glitched = packet[0].glitched.select(m64x4::splat(true), z_norm.le(f64x4::splat(reference.z_tolerance[(iteration - reference.start_iteration)])));
        //
        //             packet[0].escaped |= packet[0].glitched;
        //
        //             if packet[0].escaped.all() {
        //                 break;
        //             }
        //
        //             packet[0].derivative_current = z * packet[0].derivative_current * f64x4::splat(2.0) + ComplexVector::<f64x4>::splat(ComplexFixed::new(1.0, 0.0));
        //             packet[0].delta_current = ComplexVector::<f64x4>::splat(2.0 * reference.z_reference[(iteration - reference.start_iteration)]) * packet[0].delta_current + packet[0].delta_current * packet[0].delta_current + packet[0].delta_reference;
        //             packet[0].iteration += packet[0].escaped.select(u64x4::splat(0), u64x4::splat(1));
        //         }
        //     });
        //
        // for i in (0..self.image_x.len()).step_by(4) {
        //     for j in 0..4 {
        //         self.iteration[i + j] = packets[i / 4].iteration.extract(j) as usize;
        //         self.escaped[i + j] = packets[i / 4].escaped.extract(j);
        //         self.glitched[i + j] = packets[i / 4].glitched.extract(j);
        //         self.derivative_current[i + j] = ComplexFixed::new(packets[i / 4].derivative_current.re.extract(j), packets[i / 4].derivative_current.im.extract(j));
        //         self.delta_current[i + j] = ComplexFixed::new(packets[i / 4].delta_current.re.extract(j), packets[i / 4].delta_current.im.extract(j));
        //     }
        // }




        // TODO convert this into vectorised code
        // for i in 0..self.image_x.len() {
        //     while self.iteration[i] < maximum_iteration {
        //         // This uses the difference between the starting iteration of the reference - can be used to skip some
        //         let z = self.delta_current[i] + reference.z_reference[(self.iteration[i] - reference.start_iteration) as usize];
        //         let z_norm = z.norm_sqr();
        //
        //         if z_norm < reference.z_tolerance[(self.iteration[i] - reference.start_iteration) as usize] {
        //             self.glitched[i] = true;
        //             self.delta_current[i] = z;
        //             break;
        //         }
        //
        //         if z_norm > 65536.0 {
        //             self.escaped[i] = true;
        //             self.delta_current[i] = z;
        //             break;
        //         }
        //
        //         self.derivative_current[i] = 2.0 * z * self.derivative_current[i] + 1.0;
        //         self.delta_current[i] = 2.0 * reference.z_reference[(self.iteration[i] - reference.start_iteration) as usize] * self.delta_current[i] + self.delta_current[i] * self.delta_current[i] + self.delta_reference[i];
        //         self.iteration[i] += 1;
        //     }
        // };
    }
}

