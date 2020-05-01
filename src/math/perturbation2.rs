use crate::math::reference::Reference;
use crate::util::{PixelData, PixelData2};
use float_extended::util::FloatExp;

use rayon::prelude::*;
use float_extended::float_extended::FloatExtended;
use crate::math::reference2::Reference2;
use std::cmp::max;

pub struct Perturbation2 {}

impl Perturbation2 {
    pub fn iterate(pixel_data: &mut Vec<PixelData2>, reference: &Reference, maximum_iteration: usize) {
        pixel_data.par_chunks_mut(1)
            .for_each(|pixel_data| {
                for packet in pixel_data {
                    let original_reference = (packet.delta_reference, packet.p_initial);

                    // normalise so that the current d0 is same as current dn
                    packet.delta_reference *= 2.0f64.powi(packet.p_initial - packet.p_current);
                    packet.p_initial = packet.p_current;

                    // normal but with added exponent
                    while packet.iteration < maximum_iteration {
                        // This uses the difference between the starting iteration of the reference - can be used to skip some
                        // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;

                        let z_reference = reference.z_reference[packet.iteration - reference.start_iteration];

                        if packet.p_current > -800 {
                            // The -1 in the power here is because the reference is already pre-multiplied by 2
                            let z = z_reference.0 * 2.0f64.powi(z_reference.1 - 1) + packet.delta_current * 2.0f64.powi(packet.p_current);
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

                            // We have already done this check, we only need to check if the current value is in floatexp form
                            if z_reference.1 < 0 {
                                // this is the new p value we will be using
                                let new_p = max(packet.p_current + z_reference.1, original_reference.1);

                                packet.delta_reference = original_reference.0 * 2.0f64.powi(original_reference.1 - new_p);

                                // for the delta squared bit, as we want it in terms of the new p level, the powers are 2^(2 * pinit) / 2^(new_p)
                                packet.delta_current = packet.delta_current * z_reference.0 * 2.0f64.powi(packet.p_current + z_reference.1 - new_p) + packet.delta_current * packet.delta_current * 2.0f64.powi(2 * packet.p_current - new_p) + packet.delta_reference;
                                packet.p_current = new_p;
                            } else {
                                // for the delta squared bit, as we want it in terms of the current p level, the powers are 2^(2 * pinit) / 2^(pinit)
                                packet.delta_current = z_reference.0 * packet.delta_current + packet.delta_current * packet.delta_current * 2.0f64.powi(packet.p_current) + packet.delta_reference;
                            }
                        } else {
                            // here, we have that the exponent gets very low - need to do a relative thing

                            // We have already done this check, we only need to check if the current value is in floatexp form
                            if z_reference.1 < 0 {
                                // we rescale back to the initial delta level

                                // this is the new p value we will be using
                                let new_p = max(packet.p_current + z_reference.1, original_reference.1);

                                packet.delta_reference = original_reference.0 * 2.0f64.powi(original_reference.1 - new_p);

                                packet.delta_current = packet.delta_current * z_reference.0 * 2.0f64.powi(packet.p_current + z_reference.1 - new_p) + packet.delta_current * packet.delta_current * 2.0f64.powi(2 * packet.p_current - new_p) + packet.delta_reference;
                                packet.p_current = new_p;
                            } else {
                                packet.delta_current = z_reference.0 * packet.delta_current + packet.delta_reference;
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

                            // packet.derivative_current /= 2.0f64.powi(added_exponent);
                        }
                    }
                }
            });
    }
}



