use crate::renderer::FloatVector;

#[derive(Copy, Clone)]
pub struct Complexf32 {
    pub real: f32,
    pub imaginary: f32
}

impl Complexf32 {
    pub fn new(real: f32, imaginary: f32) -> Self {
        Complexf32 {
            real,
            imaginary
        }
    }

    pub fn norm(&self) -> f32 {
        self.real * self.real + self.imaginary * self.imaginary
    }

    pub fn square(&self) -> Complexf32 {
        Complexf32 {
            real: self.real * self.real - self.imaginary * self.imaginary,
            imaginary: 2.0 * self.real * self.imaginary
        }
    }
}

impl std::ops::Add<Complexf32> for Complexf32 {
    type Output = Complexf32;

    fn add(self, other: Complexf32) -> Complexf32 {
        Complexf32 {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary
        }
    }
}

impl std::ops::Mul<Complexf32> for Complexf32 {
    type Output = Complexf32;

    fn mul(self, other: Complexf32) -> Complexf32 {
        Complexf32 {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real
        }
    }
}

impl std::ops::Mul<f32> for Complexf32 {
    type Output = Complexf32;

    fn mul(self, other: f32) -> Complexf32 {
        Complexf32 {
            real: self.real * other,
            imaginary: self.imaginary * other
        }
    }
}

impl std::ops::Mul<Complexf32> for f32 {
    type Output = Complexf32;

    fn mul(self, other: Complexf32) -> Complexf32 {
        Complexf32 {
            real: self * other.real,
            imaginary: self * other.imaginary
        }
    }
}

impl std::ops::Sub<Complexf32> for Complexf32 {
    type Output = Complexf32;

    fn sub(self, other: Complexf32) -> Complexf32 {
        Complexf32 {
            real: self.real - other.real,
            imaginary: self.imaginary - other.imaginary
        }
    }
}

impl std::ops::MulAssign<Complexf32> for Complexf32 {
    fn mul_assign(&mut self, other: Complexf32) {
        *self = Complexf32 {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real
        }
    }
}

impl std::ops::AddAssign<Complexf32> for Complexf32 {
    fn add_assign(&mut self, other: Complexf32) {
        *self = Complexf32 {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary
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

    pub fn splat(value: Complexf32) -> Self {
        ComplexVector {
            real: FloatVector::splat(value.real),
            imaginary: FloatVector::splat(value.imaginary)
        }
    }

    pub fn norm(&mut self) -> FloatVector {
        self.real * self.real + self.imaginary * self.imaginary
    }

    pub fn square(&self) -> ComplexVector {
        ComplexVector {
            real: self.real * self.real - self.imaginary * self.imaginary,
            imaginary: FloatVector::splat(2.0) * self.real * self.imaginary
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

impl std::ops::Mul<f32> for ComplexVector {
    type Output = ComplexVector;

    fn mul(self, other: f32) -> ComplexVector {
        ComplexVector {
            real: self.real * other,
            imaginary: self.imaginary * other
        }
    }
}

impl std::ops::Mul<ComplexVector> for f32 {
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
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary
        }
    }
}