use crate::util::PixelDataExtended;
use float_extended::util::FloatExp;

use rayon::prelude::*;
use crate::math::reference_extended::ReferenceExtended;
use std::cmp::max;

pub struct PerturbationExtended {}

impl PerturbationExtended {
    pub fn iterate(pixel_data: &mut Vec<PixelDataExtended>, reference: &ReferenceExtended, reference_current_iteration: usize) {
        pixel_data.par_chunks_mut(4)
            .for_each(|pixel_data| {
                for packet in pixel_data {
                    let original_reference = (packet.delta_reference, packet.p_initial);

                    // normalise so that the current d0 is same as current dn
                    packet.delta_reference *= 1.0.ldexp(packet.p_initial - packet.p_current);
                    packet.p_initial = packet.p_current;
                    let mut temp = 1.0.ldexp(packet.p_current);

                    // normal but with added exponent
                    while packet.iteration < reference_current_iteration {
                        // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;
                        let z_reference = reference.z_reference[packet.iteration - reference.start_iteration];

                        if packet.p_current > -500 {
                            let mut z = 0.5 * z_reference.0;

                            if z_reference.1 < 0 {
                                z *= 1.0.ldexp(z_reference.1);
                            }

                            z += packet.delta_current * temp;
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

                            // We have already done this check, we only need to check if the current value is in floatexp form
                            if z_reference.1 < 0 {
                                // we rescale back to the initial delta level

                                // this is the new p value we will be using
                                let new_p = max(packet.p_current + z_reference.1, original_reference.1);

                                packet.delta_reference = original_reference.0 * 1.0.ldexp(original_reference.1 - new_p);

                                let temp9 = z_reference.0 * 1.0.ldexp(packet.p_current + z_reference.1 - new_p) + packet.delta_current * 1.0.ldexp(2 * packet.p_current - new_p);

                                packet.delta_current *= temp9;
                                packet.delta_current += packet.delta_reference;
                                packet.p_current = new_p;
                                temp = 1.0.ldexp(packet.p_current);
                            } else {
                                // for the delta squared bit, as we want it in terms of the current p level, the powers are 2^(2 * pinit) / 2^(pinit)
                                let temp2 = packet.delta_current * temp + z_reference.0;
                                packet.delta_current *= temp2;
                                packet.delta_current += packet.delta_reference;
                            }
                        } else {
                            // We have already done this check, we only need to check if the current value is in floatexp form
                            if z_reference.1 < 0 {
                                // we rescale back to the initial delta level

                                // this is the new p value we will be using
                                let new_p = max(packet.p_current + z_reference.1, original_reference.1);

                                packet.delta_reference = original_reference.0 * 1.0.ldexp(original_reference.1 - new_p);

                                let temp9 = z_reference.0 * 1.0.ldexp(packet.p_current + z_reference.1 - new_p) + packet.delta_current * 1.0.ldexp(2 * packet.p_current - new_p);

                                packet.delta_current *= temp9;
                                packet.delta_current += packet.delta_reference;
                                packet.p_current = new_p;
                                temp = 1.0.ldexp(packet.p_current);
                            } else {
                                packet.delta_current *= z_reference.0;
                                packet.delta_current += packet.delta_reference;
                            }
                        }

                        packet.iteration += 1;

                        // 2x time with 100 vs 1
                        if packet.iteration % 400 == 0 {
                            // here we should do the rescale, based on the real value
                            let (temp_mantissa, added_exponent) = packet.delta_current.re.frexp();
                            packet.delta_current.re = temp_mantissa;

                            packet.delta_current.im = packet.delta_current.im.ldexp(-added_exponent);
                            packet.p_current += added_exponent;
                            packet.delta_reference.re = packet.delta_reference.re.ldexp(-added_exponent);
                            packet.delta_reference.im = packet.delta_reference.im.ldexp(-added_exponent);
                            temp = 1.0.ldexp(packet.p_current);
                            // packet.derivative_current /= 1.0.ldexp(added_exponent);
                        }
                    }
                }
            });
    }
}



