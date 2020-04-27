use crate::math::reference::Reference;
use crate::util::PixelData;

use rayon::prelude::*;

pub struct Perturbation {}

impl Perturbation {
    pub fn iterate(pixel_data: &mut Vec<PixelData>, reference: &Reference, maximum_iteration: usize) {
        pixel_data.par_chunks_mut(1)
            .for_each(|pixel_data| {
                for packet in pixel_data {
                    // TODO investigate setting up a floatexp-type thing which only rescales every 500 or so iterations

                    // normal
                    while packet.iteration < maximum_iteration {
                        // This uses the difference between the starting iteration of the reference - can be used to skip some
                        let z = packet.delta_current + reference.z_reference[packet.iteration - reference.start_iteration];
                        let z_norm = z.norm_sqr();

                        if z_norm < reference.z_tolerance[packet.iteration - reference.start_iteration] {
                            packet.glitched = true;
                            packet.delta_current = z;
                            break;
                        }

                        if z_norm > 65536.0 {
                            packet.escaped = true;
                            packet.delta_current = z;
                            break;
                        }

                        packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;
                        packet.delta_current = 2.0 * reference.z_reference[packet.iteration - reference.start_iteration] * packet.delta_current + packet.delta_current * packet.delta_current + packet.delta_reference;
                        packet.iteration += 1;
                    }
                }
            });
    }
}



