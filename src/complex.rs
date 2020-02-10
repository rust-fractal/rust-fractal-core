use rug::Float;

#[derive(Copy, Clone)]
pub struct ComplexF64 {
    pub real: f64,
    pub imaginary: f64
}

impl ComplexF64 {
    pub fn new(real: f64, imaginary: f64) -> Self {
        ComplexF64 {
            real,
            imaginary
        }
    }

    pub fn norm(&self) -> f64 {
        self.real * self.real + self.imaginary * self.imaginary
    }

    pub fn square(&self) -> ComplexF64 {
        ComplexF64 {
            real: self.real * self.real - self.imaginary * self.imaginary,
            imaginary: 2.0 * self.real * self.imaginary
        }
    }
}

impl std::ops::Add<ComplexF64> for ComplexF64 {
    type Output = ComplexF64;

    fn add(self, other: ComplexF64) -> ComplexF64 {
        ComplexF64 {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary
        }
    }
}

impl std::ops::Mul<ComplexF64> for ComplexF64 {
    type Output = ComplexF64;

    fn mul(self, other: ComplexF64) -> ComplexF64 {
        ComplexF64 {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real
        }
    }
}

impl std::ops::Mul<f64> for ComplexF64 {
    type Output = ComplexF64;

    fn mul(self, other: f64) -> ComplexF64 {
        ComplexF64 {
            real: self.real * other,
            imaginary: self.imaginary * other
        }
    }
}

impl std::ops::Mul<ComplexF64> for f64 {
    type Output = ComplexF64;

    fn mul(self, other: ComplexF64) -> ComplexF64 {
        ComplexF64 {
            real: self * other.real,
            imaginary: self * other.imaginary
        }
    }
}

impl std::ops::Sub<ComplexF64> for ComplexF64 {
    type Output = ComplexF64;

    fn sub(self, other: ComplexF64) -> ComplexF64 {
        ComplexF64 {
            real: self.real - other.real,
            imaginary: self.imaginary - other.imaginary
        }
    }
}

impl std::ops::MulAssign<ComplexF64> for ComplexF64 {
    fn mul_assign(&mut self, other: ComplexF64) {
        *self = ComplexF64 {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real
        }
    }
}

impl std::ops::AddAssign<ComplexF64> for ComplexF64 {
    fn add_assign(&mut self, other: ComplexF64) {
        *self = ComplexF64 {
            real: self.real * other.real,
            imaginary: self.imaginary * other.imaginary
        }
    }
}