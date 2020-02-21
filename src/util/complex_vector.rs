use crate::util::ComplexFixed;

#[derive(Copy, Clone)]
pub struct ComplexVector<T> {
    pub re: T,
    pub im: T
}

impl ComplexVector<f32x8> {
    #[inline]
    pub fn new(re: &[f32], im: &[f32]) -> Self {
        ComplexVector {
            re: f32x8::from_slice_unaligned(re),
            im: f32x8::from_slice_unaligned(im)
        }
    }

    #[inline]
    pub fn splat(value: ComplexFixed<f32>) -> Self {
        ComplexVector {
            re: f32x8::splat(value.re),
            im: f32x8::splat(value.im)
        }
    }

    #[inline]
    pub fn norm_sqr(&self) -> f32x8 {
        self.re * self.re + self.im * self.im
    }
}

impl ComplexVector<f64x4> {
    #[inline]
    pub fn new(re: &[f64], im: &[f64]) -> Self {
        ComplexVector {
            re: f64x4::from_slice_unaligned(re),
            im: f64x4::from_slice_unaligned(im)
        }
    }

    #[inline]
    pub fn splat(value: ComplexFixed<f64>) -> Self {
        ComplexVector {
            re: f64x4::splat(value.re),
            im: f64x4::splat(value.im)
        }
    }

    #[inline]
    pub fn splat2(value: ComplexFixed<f32>) -> Self {
        ComplexVector {
            re: f64x4::splat(value.re as f64),
            im: f64x4::splat(value.im as f64)
        }
    }

    #[inline]
    pub fn norm_sqr(&self) -> f64x4 {
        self.re * self.re + self.im * self.im
    }
}

impl std::ops::Add<ComplexVector<f32x8>> for ComplexVector<f32x8> {
    type Output = ComplexVector<f32x8>;

    #[inline]
    fn add(self, other: ComplexVector<f32x8>) -> ComplexVector<f32x8> {
        ComplexVector {
            re: self.re + other.re,
            im: self.im + other.im
        }
    }
}

impl std::ops::Add<ComplexVector<f64x4>> for ComplexVector<f64x4> {
    type Output = ComplexVector<f64x4>;

    #[inline]
    fn add(self, other: ComplexVector<f64x4>) -> ComplexVector<f64x4> {
        ComplexVector {
            re: self.re + other.re,
            im: self.im + other.im
        }
    }
}

impl std::ops::Sub<ComplexVector<f32x8>> for ComplexVector<f32x8> {
    type Output = ComplexVector<f32x8>;

    #[inline]
    fn sub(self, other: ComplexVector<f32x8>) -> ComplexVector<f32x8> {
        ComplexVector {
            re: self.re - other.re,
            im: self.im - other.im
        }
    }
}

impl std::ops::Sub<ComplexVector<f64x4>> for ComplexVector<f64x4> {
    type Output = ComplexVector<f64x4>;

    #[inline]
    fn sub(self, other: ComplexVector<f64x4>) -> ComplexVector<f64x4> {
        ComplexVector {
            re: self.re - other.re,
            im: self.im - other.im
        }
    }
}

impl std::ops::Mul<ComplexVector<f32x8>> for ComplexVector<f32x8> {
    type Output = ComplexVector<f32x8>;

    #[inline]
    fn mul(self, other: ComplexVector<f32x8>) -> ComplexVector<f32x8> {
        ComplexVector {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re
        }
    }
}

impl std::ops::Mul<ComplexVector<f64x4>> for ComplexVector<f64x4> {
    type Output = ComplexVector<f64x4>;

    #[inline]
    fn mul(self, other: ComplexVector<f64x4>) -> ComplexVector<f64x4> {
        ComplexVector {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re
        }
    }
}

impl std::ops::AddAssign<ComplexVector<f32x8>> for ComplexVector<f32x8> {
    #[inline]
    fn add_assign(&mut self, other: ComplexVector<f32x8>) {
        *self = ComplexVector {
            re: self.re + other.re,
            im: self.im + other.im
        }
    }
}

impl std::ops::AddAssign<ComplexVector<f64x4>> for ComplexVector<f64x4> {
    #[inline]
    fn add_assign(&mut self, other: ComplexVector<f64x4>) {
        *self = ComplexVector {
            re: self.re + other.re,
            im: self.im + other.im
        }
    }
}

impl std::ops::MulAssign<ComplexVector<f32x8>> for ComplexVector<f32x8> {
    #[inline]
    fn mul_assign(&mut self, other: ComplexVector<f32x8>) {
        *self = ComplexVector {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re
        }
    }
}

impl std::ops::MulAssign<ComplexVector<f64x4>> for ComplexVector<f64x4> {
    #[inline]
    fn mul_assign(&mut self, other: ComplexVector<f64x4>) {
        *self = ComplexVector {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re
        }
    }
}