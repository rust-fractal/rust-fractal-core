use rug::Float;
use crate::renderer::FloatVector;

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

#[derive(Copy, Clone)]
pub struct ComplexVector {
    pub real: FloatVector,
    pub imaginary: FloatVector
}

impl ComplexVector {
    pub fn new(real: FloatVector, imaginary: FloatVector) -> Self {
        ComplexVector {
            real,
            imaginary
        }
    }

    pub fn splat(value: ComplexF64) -> Self {
        ComplexVector {
            real: FloatVector::splat(value.real),
            imaginary: FloatVector::splat(value.imaginary)
        }
    }

    pub fn norm(&self) -> FloatVector {
        self.real * self.real + self.imaginary * self.imaginary
    }

    pub fn square(&self) -> ComplexVector {
        ComplexVector {
            real: self.real * self.real - self.imaginary * self.imaginary,
            imaginary: 2.0 * self.real * self.imaginary
        }
    }
}

impl std::ops::Add<ComplexVector> for ComplexVector {
    type Output = ComplexVector;

    fn add(self, other: ComplexVector) -> ComplexVector {
        ComplexVector {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary
        }
    }
}

impl std::ops::Mul<ComplexVector> for ComplexVector {
    type Output = ComplexVector;

    fn mul(self, other: ComplexVector) -> ComplexVector {
        ComplexVector {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real
        }
    }
}

impl std::ops::Mul<f64> for ComplexVector {
    type Output = ComplexVector;

    fn mul(self, other: f64) -> ComplexVector {
        ComplexVector {
            real: self.real * other,
            imaginary: self.imaginary * other
        }
    }
}

impl std::ops::Mul<ComplexVector> for f64 {
    type Output = ComplexVector;

    fn mul(self, other: ComplexVector) -> ComplexVector {
        ComplexVector {
            real: self * other.real,
            imaginary: self * other.imaginary
        }
    }
}

impl std::ops::Sub<ComplexVector> for ComplexVector {
    type Output = ComplexVector;

    fn sub(self, other: ComplexVector) -> ComplexVector {
        ComplexVector {
            real: self.real - other.real,
            imaginary: self.imaginary - other.imaginary
        }
    }
}

impl std::ops::MulAssign<ComplexVector> for ComplexVector {
    fn mul_assign(&mut self, other: ComplexVector) {
        *self = ComplexVector {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real
        }
    }
}

impl std::ops::AddAssign<ComplexVector> for ComplexVector {
    fn add_assign(&mut self, other: ComplexVector) {
        *self = ComplexVector {
            real: self.real * other.real,
            imaginary: self.imaginary * other.imaginary
        }
    }
}