use crate::math::reference::Reference;
use crate::util::{PixelData, PixelData2};
use float_extended::util::FloatExp;

use rayon::prelude::*;
use float_extended::float_extended::FloatExtended;

pub struct Perturbation2 {}

impl Perturbation2 {
    pub fn iterate(pixel_data: &mut Vec<PixelData2>, reference: &Reference, maximum_iteration: usize) {
        pixel_data.chunks_mut(1)
            .for_each(|pixel_data| {
                for packet in pixel_data {
                    // TODO investigate setting up a floatexp-type thing which only rescales every 500 or so iterations

                    let temp_ref = packet.delta_reference;
                    let temp_p = packet.p_initial;

                    // -500 is temporary (needs to be basically 300 * 3.13 or something)
                    while packet.p_initial < -700 && packet.iteration < maximum_iteration {
                        // let z = packet.delta_current + reference.z_reference[packet.iteration - reference.start_iteration];

                        // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;
                        // we need to do a check here to see if the z reference is 0
                        let temp = reference.z_reference[packet.iteration - reference.start_iteration];

                        if temp.norm_sqr() == 0.0 {
                            packet.delta_current = temp_ref;
                            packet.p_initial = temp_p;
                            packet.delta_reference = temp_ref;
                        } else {
                            packet.delta_current = 2.0 * temp * packet.delta_current + packet.delta_reference;
                        }
                        packet.iteration += 1;

                        if packet.iteration % 400 == 0 {
                            // here we should do the rescale, based on the real value
                            let (temp_mantissa, added_exponent) = packet.delta_current.re.frexp();
                            packet.delta_current.re = temp_mantissa;
                            // packet.delta_current.im *= 2.0f64.powi(-added_exponent);
                            // packet.p_initial += added_exponent;
                            // packet.delta_reference *= 2.0f64.powi(-added_exponent);


                            packet.delta_current.im = packet.delta_current.im.ldexp(-added_exponent);
                            packet.p_initial += added_exponent;
                            packet.delta_reference.re = packet.delta_reference.re.ldexp(-added_exponent);
                            packet.delta_reference.im = packet.delta_reference.im.ldexp(-added_exponent);


                            // packet.derivative_current /= 2.0f64.powi(added_exponent);
                        }
                    }

                    // ok, so now we should have the delta sized enough to fit in double
                    packet.delta_current.re = packet.delta_current.re.ldexp(packet.p_initial);
                    packet.delta_current.im = packet.delta_current.im.ldexp(packet.p_initial);
                    packet.delta_reference.re = packet.delta_reference.re.ldexp(packet.p_initial);
                    packet.delta_reference.im = packet.delta_reference.im.ldexp(packet.p_initial);
                    // packet.derivative_current *= 2.0f64.powi(packet.p_initial);

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

                        // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;
                        packet.delta_current = 2.0 * reference.z_reference[packet.iteration - reference.start_iteration] * packet.delta_current + packet.delta_current * packet.delta_current + packet.delta_reference;
                        packet.iteration += 1;
                    }
                }
            });
    }
}



