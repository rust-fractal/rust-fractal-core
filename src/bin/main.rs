use std::ops::{Add, Mul, MulAssign, AddAssign};

#[derive(Copy, Clone)]
struct Complex {
    pub real: f64,
    pub complex: f64
}

impl Complex {
    pub fn new(real: f64, complex: f64) -> Self {
        Complex {
            real,
            complex
        }
    }

    pub fn length(&self) -> f64 {
        (self.real * self.real + self.complex * self.complex).sqrt()
    }

    pub fn length_squared(&self) -> f64 {
        self.real * self.real + self.complex * self.complex
    }

    pub fn square(&mut self) {
        self.real = self.real.powf(2.0);
        self.complex = self.complex.powf(2.0);
    }
}

impl Add<Complex> for Complex {
    type Output = Complex;

    fn add(self, other: Complex) -> Complex {
        Complex {
            real: self.real + other.real,
            complex: self.complex + other.complex
        }
    }
}

impl Mul<Complex> for Complex {
    type Output = Complex;

    fn mul(self, other: Complex) -> Complex {
        Complex {
            real: self.real * other.real - self.complex * other.complex,
            complex: self.real * other.complex + self.complex * other.real
        }
    }
}

impl AddAssign<Complex> for Complex {
    fn add_assign(&mut self, other: Complex) {
        self.real += other.real;
        self.complex += other.complex
    }
}

impl Mul<Complex> for f64 {
    type Output = Complex;

    fn mul(self, other: Complex) -> Complex {
        Complex {
            real: self * other.real,
            complex: self * other.complex
        }
    }
}

struct Sample {
    z: Complex,
    c: Complex,
    iterations: usize,
    max_iterations: usize
}

impl Sample {
    pub fn new(c: Complex, max_iterations: usize) -> Self {
        Sample {
            z: c,
            c,
            iterations: 0,
            max_iterations
        }
    }
}

impl Iterator for Sample {
    type Item = Complex;

    fn next(&mut self) -> Option<Complex> {
        if self.z.length_squared() > 4.0 || self.iterations > self.max_iterations {
            None
        } else {
            self.z.square();
            self.z += self.c;
            self.iterations += 1;
            Some(self.z)
        }
    }
}

fn main() {
    let width: usize = 50;
    let height: usize = 25;

    let top_left = Complex::new(-1.0, 1.0);
    let bottom_right = Complex::new(1.0, -1.0);

    let mut iterations_buffer = Vec::new();

    for j in 0..height {
        for i in 0..width {
            let c = top_left + (i as f64 / width as f64) * Complex::new(2.0, 0.0) + (j as f64 / height as f64) * Complex::new(0.0, 2.0);
            let iterations = Sample::new(c, 255).count();

            iterations_buffer.push(iterations as u8);
        }
    }

    image::save_buffer("output.png", &iterations_buffer, width as u32, height as u32, image::Gray(8)).unwrap();
}