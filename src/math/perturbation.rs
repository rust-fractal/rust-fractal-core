use crate::math::reference_arbitrary::Reference;
use crate::util::ComplexFixed;

pub struct Perturbation {
    pub image_x: Vec<i32>,
    pub image_y: Vec<i32>,
    pub iteration: Vec<i32>,
    pub dc: Vec<ComplexFixed<f64>>,
    pub dz: Vec<ComplexFixed<f64>>,
    pub dzdc: Vec<ComplexFixed<f64>>,
    pub glitched: Vec<bool>,
    pub escaped: Vec<bool>
}

impl Perturbation {
    pub fn new(image_x: Vec<i32>, image_y: Vec<i32>, iteration: Vec<i32>, dc: Vec<ComplexFixed<f64>>, dz: Vec<ComplexFixed<f64>>, dzdc: Vec<ComplexFixed<f64>>) -> Self {
        let length = image_x.len();
        Perturbation {
            image_x,
            image_y,
            iteration,
            dc,
            dz,
            dzdc,
            glitched: vec![false; length],
            escaped: vec![false; length],
        }
    }

    pub fn iterate(&mut self, reference: &Reference, max_iterations: i32) {
        // TODO convert this into vectorised code
        for i in 0..self.image_x.len() {
            while self.iteration[i] < max_iterations {
                // This uses the difference between the starting iteration of the reference - can be used to skip some
                let z = self.dz[i] + reference.z_reference[(self.iteration[i] - reference.start_iteration) as usize];
                let z_norm = z.norm_sqr();

                if z_norm < reference.z_tolerance[(self.iteration[i] - reference.start_iteration) as usize] {
                    self.glitched[i] = true;
                    self.dz[i] = z;
                    break;
                }

                if z_norm > 65536.0 {
                    self.escaped[i] = true;
                    self.dz[i] = z;
                    break;
                }

                self.dzdc[i] = 2.0 * z * self.dzdc[i] + 1.0;
                self.dz[i] = 2.0 * reference.z_reference[(self.iteration[i] - reference.start_iteration) as usize] * self.dz[i] + self.dz[i] * self.dz[i] + self.dc[i];
                self.iteration[i] += 1;
            }
        }
    }
}

