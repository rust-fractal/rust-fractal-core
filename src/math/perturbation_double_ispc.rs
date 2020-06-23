use crate::math::reference_double::ReferenceDouble;
use crate::util::PixelDataDouble;
use std::time::Instant;

ispc_module!(mandelbrot);
pub struct PerturbationDoubleISPC {}

impl PerturbationDoubleISPC {
    pub fn iterate(pixel_data: &mut Vec<PixelDataDouble>, reference: &ReferenceDouble, reference_current_iteration: usize) {
        
        let length = pixel_data.len();

        let mut iterations = vec![reference.start_iteration as i32; length];
        let mut glitched = vec![0; length];
        let mut escaped = vec![0; length];

        let delta_reference_re = pixel_data.into_iter().map({
            |data| {
                data.delta_reference.re
            }
        }).collect::<Vec<f64>>();

        let delta_reference_im = pixel_data.into_iter().map({
            |data| {
                data.delta_reference.im
            }
        }).collect::<Vec<f64>>();

        let mut delta_current_re = pixel_data.into_iter().map({
            |data| {
                data.delta_current.re
            }
        }).collect::<Vec<f64>>();

        let mut delta_current_im = pixel_data.into_iter().map({
            |data| {
                data.delta_current.im
            }
        }).collect::<Vec<f64>>();

        let reference_re = reference.z_reference.iter().map({
            |data| {
                data.re
            }
        }).collect::<Vec<f64>>();

        let reference_im = reference.z_reference.iter().map({
            |data| {
                data.im
            }
        }).collect::<Vec<f64>>();


        // let time = Instant::now();

        unsafe {
            mandelbrot::perturbation_double_ispc(
                64,
                (length / 64) as i32,
                iterations.as_mut_ptr(),
                delta_reference_re.as_ptr(), 
                delta_reference_im.as_ptr(), 
                delta_current_re.as_mut_ptr(), 
                delta_current_im.as_mut_ptr(), 
                glitched.as_mut_ptr(), 
                escaped.as_mut_ptr(), 
                reference_current_iteration as i32, 
                reference.start_iteration as i32,
                reference_re.as_ptr(),
                reference_im.as_ptr(),
                reference.z_tolerance.as_ptr()
            )
        }

        // println!("{:<14}{:>6} ms", "Actual Iteration", time.elapsed().as_millis());

        // println!("{}", iterations[0]);

        for i in 0..length {
            pixel_data[i].iteration = iterations[i] as usize;
            pixel_data[i].delta_current.re = delta_current_re[i];
            pixel_data[i].delta_current.im = delta_current_im[i];
            pixel_data[i].glitched = glitched[i] != 0;
            pixel_data[i].escaped = escaped[i] != 0;
        }
    }
}



