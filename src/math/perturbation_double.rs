use crate::math::reference_double::ReferenceDouble;
use crate::util::PixelDataDouble;
use rayon::prelude::*;

pub struct PerturbationDouble {}

impl PerturbationDouble {
    pub fn iterate(pixel_data: &mut Vec<PixelDataDouble>, reference: &ReferenceDouble, reference_current_iteration: usize) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                for packet in pixel_data {
                    // Normal loop for single elements
                    while packet.iteration < reference_current_iteration {
                        // This uses the difference between the starting iteration of the reference - can be used to skip some

                        let z = packet.delta_current + reference.z_reference[packet.iteration - reference.start_iteration];
                        let z_norm = z.norm_sqr();

                        if z_norm < reference.z_tolerance[packet.iteration - reference.start_iteration] {
                            packet.glitched = true;
                            packet.delta_current = z;
                            break;
                        }

                        if z_norm > 1e16 {
                            packet.escaped = true;
                            packet.delta_current = z;
                            break;
                        }

                        packet.delta_current = 2.0 * reference.z_reference[packet.iteration - reference.start_iteration] * packet.delta_current + packet.delta_current * packet.delta_current + packet.delta_reference;

                        packet.derivative_current *= 2.0 * z;
                        packet.derivative_current += 1.0;

                        packet.iteration += 1;
                    }
                }
            });
    }
}



