use crate::util::{ComplexArbitrary, ComplexFixed, to_fixed_exp};

pub struct Reference2 {
    pub start_iteration: usize,
    pub current_iteration: usize,
    pub maximum_iteration: usize,
    pub z: ComplexArbitrary,
    pub c: ComplexArbitrary,
    pub z_reference: Vec<(ComplexFixed<f64>, i32)>,
    pub z_tolerance: Vec<f64>
}

impl Reference2 {
    pub fn new(z: ComplexArbitrary, c: ComplexArbitrary, current_iteration: usize, maximum_iteration: usize) -> Reference2 {
        let z_fixed = to_fixed_exp(&z);

        // 1e-6 is the threshold for pauldelbrot's criterion
        Reference2 {
            start_iteration: current_iteration,
            current_iteration,
            maximum_iteration,
            z,
            c,
            z_reference: vec![z_fixed],
            z_tolerance: vec![1e-6 * z_fixed.0.norm_sqr() * 2.0f64.powi(2 * z_fixed.1)]
        }
    }

    pub fn step(&mut self) -> bool {
        self.z = self.z.clone().square() + self.c.clone();
        self.current_iteration += 1;
        let z_fixed = to_fixed_exp(&self.z);

        self.z_reference.push(z_fixed);
        self.z_tolerance.push(1e-6 * z_fixed.0.norm_sqr() * 2.0f64.powi(2 * z_fixed.1));
        // println!("{}", z_fixed);

        (z_fixed.0.norm_sqr() * 2.0f64.powi(2 * z_fixed.1)) <= 4.0
    }


    pub fn run(&mut self) -> bool {
        while self.current_iteration <= self.maximum_iteration {
            if !self.step() {
                break;
            }
        };

        let mut temp = 0;

        for element in &self.z_reference {
            if element.1 < temp {
                temp = element.1;
            }
        }

        println!("element: {}", temp);

        self.current_iteration == self.maximum_iteration
    }
}

