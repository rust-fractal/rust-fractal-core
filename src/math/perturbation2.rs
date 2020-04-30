use crate::math::reference::Reference;
use crate::util::{PixelData, PixelData2};
use float_extended::util::FloatExp;

use rayon::prelude::*;
use float_extended::float_extended::FloatExtended;
use crate::math::reference2::Reference2;
use std::f64::{NAN, INFINITY};
use std::cmp::max;

pub struct Perturbation2 {}

impl Perturbation2 {
    pub fn iterate(pixel_data: &mut Vec<PixelData2>, reference: &Reference, maximum_iteration: usize) {
        pixel_data.par_chunks_mut(1)
            .for_each(|pixel_data| {
                for packet in pixel_data {
                    // TODO investigate setting up a floatexp-type thing which only rescales every 500 or so iterations
                    // -500 is temporary (needs to be basically 300 * 3.13 or something)

                    let temp2 = packet.delta_reference;
                    let temp3 = packet.p_initial;

                    // normal but with added exponent
                    while packet.iteration < maximum_iteration {
                        // This uses the difference between the starting iteration of the reference - can be used to skip some
                        // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;

                        let temp = reference.z_reference[packet.iteration - reference.start_iteration];

                        if packet.p_initial > -800 {
                            let z = temp + packet.delta_current * 2.0f64.powi(packet.p_initial);
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

                            // packet.p_initial += temp.1;

                            // could be removed
                            // packet.delta_reference *= 2.0f64.powi(-temp.1);

                            if temp.re <= 1e-30 && temp.re >= -1e-30 {
                                let temp4 = reference.z_reference_extended[packet.iteration - reference.start_iteration];

                                // this is the new p value we will be using
                                let new_p = max(packet.p_initial + temp4.1, temp3);

                                packet.delta_reference = temp2 * 2.0f64.powi(temp3 - new_p);


                                // for the delta squared bit, as we want it in terms of the new p level, the powers are 2^(2 * pinit) / 2^(new_p)
                                packet.delta_current = 2.0 * packet.delta_current * temp4.0 * 2.0f64.powi(packet.p_initial + temp4.1 - new_p) + packet.delta_current * packet.delta_current * 2.0f64.powi(2 * packet.p_initial - new_p) + packet.delta_reference;
                                packet.p_initial = new_p;
                            } else {
                                // for the delta squared bit, as we want it in terms of the current p level, the powers are 2^(2 * pinit) / 2^(pinit)
                                packet.delta_current = 2.0 * temp * packet.delta_current + packet.delta_current * packet.delta_current * 2.0f64.powi(packet.p_initial) + packet.delta_reference;
                            }
                        } else {
                            // here, we have that the exponent gets very low - need to do a relative thing

                            // increment the exponent value
                            // packet.p_initial += temp.1;

                            // could be removed
                            // packet.delta_reference *= 2.0f64.powi(-temp.1);

                            // maybe check for near zero rather than 0
                            if temp.re <= 1e-30 && temp.re >= -1e-30 {
                                // we need to do a rescale here
                                // if we use doubles everywhere except near zero, where we need the extra precision

                                // we rescale back to the initial delta level

                                let temp4 = reference.z_reference_extended[packet.iteration - reference.start_iteration];

                                // this is the new p value we will be using
                                let new_p = max(packet.p_initial + temp4.1, temp3);

                                packet.delta_reference = temp2 * 2.0f64.powi(temp3 - new_p);
                                packet.delta_current = 2.0 * packet.delta_current * temp4.0 * 2.0f64.powi(packet.p_initial + temp4.1 - new_p) + packet.delta_reference;
                                packet.p_initial = new_p;
                            } else {
                                packet.delta_current = 2.0 * temp * packet.delta_current + packet.delta_reference;
                            }


                        }

                        packet.iteration += 1;

                        // 2x time with 100 vs 1
                        if packet.iteration % 100 == 0 {
                            // here we should do the rescale, based on the real value
                            let (temp_mantissa, added_exponent) = packet.delta_current.re.frexp();
                            packet.delta_current.re = temp_mantissa;

                            packet.delta_current.im = packet.delta_current.im.ldexp(-added_exponent);
                            packet.p_initial += added_exponent;
                            packet.delta_reference.re = packet.delta_reference.re.ldexp(-added_exponent);
                            packet.delta_reference.im = packet.delta_reference.im.ldexp(-added_exponent);


                            // packet.derivative_current /= 2.0f64.powi(added_exponent);
                        }
                    }

                    println!("{}, {}", packet.iteration, packet.delta_current);


                    // while packet.p_initial < -800 && packet.iteration < maximum_iteration {
                    //     // let z = packet.delta_current + reference.z_reference[packet.iteration - reference.start_iteration];
                    //
                    //     // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;
                    //     // we need to do a check here to see if the z reference is 0
                    //     let mut temp = reference.z_reference[packet.iteration - reference.start_iteration];
                    //
                    //     packet.delta_current = 2.0 * temp * packet.delta_current + packet.delta_reference;
                    //     packet.iteration += 1;
                    //
                    //     if packet.iteration % 400 == 0 {
                    //         // here we should do the rescale, based on the real value
                    //         let (temp_mantissa, added_exponent) = packet.delta_current.re.frexp();
                    //         packet.delta_current.re = temp_mantissa;
                    //         // packet.delta_current.im *= 2.0f64.powi(-added_exponent);
                    //         // packet.p_initial += added_exponent;
                    //         // packet.delta_reference *= 2.0f64.powi(-added_exponent);
                    //
                    //
                    //         packet.delta_current.im = packet.delta_current.im.ldexp(-added_exponent);
                    //         packet.p_initial += added_exponent;
                    //         packet.delta_reference.re = packet.delta_reference.re.ldexp(-added_exponent);
                    //         packet.delta_reference.im = packet.delta_reference.im.ldexp(-added_exponent);
                    //
                    //
                    //         // packet.derivative_current /= 2.0f64.powi(added_exponent);
                    //     }
                    // }
                    //
                    // // what if we use this as an acceleration for normal floatexp; we first do the "shortcut"
                    // // iterations then do the normal
                    //
                    // // ok, so now we should have the delta sized enough to fit in double
                    // // packet.delta_current.re = packet.delta_current.re.ldexp(packet.p_initial);
                    // // packet.delta_current.im = packet.delta_current.im.ldexp(packet.p_initial);
                    // // packet.delta_reference.re = packet.delta_reference.re.ldexp(packet.p_initial);
                    // // packet.delta_reference.im = packet.delta_reference.im.ldexp(packet.p_initial);
                    // // packet.derivative_current *= 2.0f64.powi(packet.p_initial);
                    //
                    // // normal but with added exponent
                    // while packet.iteration < maximum_iteration {
                    //     // This uses the difference between the starting iteration of the reference - can be used to skip some
                    //     let z = reference.z_reference[packet.iteration - reference.start_iteration] + packet.delta_current * 2.0f64.powi(packet.p_initial);
                    //     let z_norm = z.norm_sqr();
                    //
                    //     if z_norm < reference.z_tolerance[packet.iteration - reference.start_iteration] {
                    //         packet.glitched = true;
                    //         packet.delta_current = z;
                    //         break;
                    //     }
                    //
                    //     if z_norm > 65536.0 {
                    //         packet.escaped = true;
                    //         packet.delta_current = z;
                    //         break;
                    //     }
                    //
                    //     // packet.derivative_current = 2.0 * z * packet.derivative_current + 1.0;
                    //
                    //     packet.delta_current = 2.0 * reference.z_reference[packet.iteration - reference.start_iteration] * packet.delta_current + packet.delta_current * packet.delta_current * 2.0f64.powi(packet.p_initial) + packet.delta_reference;
                    //
                    //     packet.iteration += 1;
                    //
                    //     // 2x time with 100 vs 1
                    //     if packet.iteration % 400 == 0 {
                    //         // here we should do the rescale, based on the real value
                    //         let (temp_mantissa, added_exponent) = packet.delta_current.re.frexp();
                    //         packet.delta_current.re = temp_mantissa;
                    //         // packet.delta_current.im *= 2.0f64.powi(-added_exponent);
                    //         // packet.p_initial += added_exponent;
                    //         // packet.delta_reference *= 2.0f64.powi(-added_exponent);
                    //
                    //
                    //         packet.delta_current.im = packet.delta_current.im.ldexp(-added_exponent);
                    //         packet.p_initial += added_exponent;
                    //         packet.delta_reference.re = packet.delta_reference.re.ldexp(-added_exponent);
                    //         packet.delta_reference.im = packet.delta_reference.im.ldexp(-added_exponent);
                    //
                    //
                    //         // packet.derivative_current /= 2.0f64.powi(added_exponent);
                    //     }
                    // }

                }
            });
    }
}



