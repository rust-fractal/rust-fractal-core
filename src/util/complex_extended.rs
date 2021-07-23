use crate::util::float_extended::FloatExtended;
use std::ops::{Mul, Add, AddAssign, MulAssign, Sub, DivAssign, SubAssign, Div};
use num_complex::Complex;
use std::fmt::{Display, Result, Formatter};
use crate::util::FloatExp;

#[derive(Debug, Copy, Clone)]
pub struct ComplexExtended {
    pub mantissa: Complex<f64>,
    pub exponent: i32,
}

impl ComplexExtended {
    #[inline]
    pub fn new(mantissa: Complex<f64>, exponent: i32) -> Self {
        ComplexExtended {
            mantissa,
            exponent
        }
    }

    #[inline]
    pub fn new2(re: f64, im: f64, exponent: i32) -> Self {
        ComplexExtended {
            mantissa: Complex::<f64>::new(re, im),
            exponent
        }
    }

    #[inline]
    pub fn norm(&self) -> FloatExtended {
        FloatExtended {
            mantissa: self.mantissa.norm(),
            exponent: self.exponent
        }
    }

    #[inline]
    pub fn norm_square(&self) -> FloatExtended {
        FloatExtended {
            mantissa: self.mantissa.norm_sqr(),
            exponent: self.exponent * 2
        }
    }

    #[inline]
    pub fn to_float(&self) -> Complex<f64> {
        self.mantissa * 1.0f64.ldexp(self.exponent)
    }

    #[inline]
    pub fn powi(&self, exp: i32) -> ComplexExtended {
        ComplexExtended {
            mantissa: self.mantissa.powi(exp),
            exponent: self.exponent * exp
        }
    }

    #[inline]
    pub fn reduce(&mut self) {
        if self.mantissa.re == 0.0 && self.mantissa.im == 0.0 {
            self.exponent = 0
        } else if self.mantissa.re * self.mantissa.re > self.mantissa.im * self.mantissa.im {
            let (temp_mantissa, added_exponent) = self.mantissa.re.frexp();
            self.mantissa.re = temp_mantissa;
            self.mantissa.im = self.mantissa.im.ldexp(-added_exponent);
            self.exponent += added_exponent;
        } else {
            let (temp_mantissa, added_exponent) = self.mantissa.im.frexp();
            self.mantissa.im = temp_mantissa;
            self.mantissa.re = self.mantissa.re.ldexp(-added_exponent);
            self.exponent += added_exponent;
        }
    }

    #[inline]
    pub fn scale_to_exponent(&mut self, exponent: i32) {
        let added_exponent = exponent - self.exponent;
        self.mantissa.re = self.mantissa.re.ldexp(added_exponent);
        self.mantissa.im = self.mantissa.im.ldexp(added_exponent);
        self.exponent = exponent;
    }
}

impl AddAssign for ComplexExtended {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        if self.exponent > other.exponent {
            self.mantissa += other.mantissa * 1.0f64.ldexp(other.exponent - self.exponent);
        } else {
            self.mantissa *= 1.0f64.ldexp(self.exponent - other.exponent);
            self.mantissa += other.mantissa;
            self.exponent = other.exponent;
        }
    }
}

impl SubAssign for ComplexExtended {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        if self.exponent > other.exponent {
            self.mantissa -= other.mantissa * 1.0f64.ldexp(other.exponent - self.exponent);
        } else {
            self.mantissa *= 1.0f64.ldexp(self.exponent - other.exponent);
            self.mantissa -= other.mantissa;
            self.exponent = other.exponent;
        }
    }
}

impl MulAssign<ComplexExtended> for ComplexExtended {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        self.mantissa *= other.mantissa;
        self.exponent += other.exponent;
    }
}

impl MulAssign<f64> for ComplexExtended {
    #[inline]
    fn mul_assign(&mut self, other: f64) {
        self.mantissa *= other;
    }
}

impl DivAssign for ComplexExtended {
    #[inline]
    fn div_assign(&mut self, other: Self) {
        self.mantissa /= other.mantissa;
        self.exponent -= other.exponent;
    }
}

impl Add<ComplexExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        if self.exponent > other.exponent {
            ComplexExtended::new(self.mantissa + other.mantissa * 1.0f64.ldexp(other.exponent - self.exponent), self.exponent)
        } else {
            ComplexExtended::new(other.mantissa + self.mantissa * 1.0f64.ldexp(self.exponent - other.exponent), other.exponent)
        }
    }
}

impl Sub<ComplexExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        if self.exponent > other.exponent {
            ComplexExtended::new(self.mantissa - other.mantissa * 1.0f64.ldexp(other.exponent - self.exponent), self.exponent)
        } else {
            ComplexExtended::new(self.mantissa * 1.0f64.ldexp(self.exponent - other.exponent) - other.mantissa, other.exponent)
        }
    }
}

impl Mul<ComplexExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        ComplexExtended::new(
            self.mantissa * other.mantissa,
            self.exponent + other.exponent
        )
    }
}

impl Div<ComplexExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn div(self, other: Self) -> Self::Output {
        ComplexExtended::new(
            self.mantissa / other.mantissa,
            self.exponent - other.exponent
        )
    }
}


impl Mul<FloatExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn mul(self, other: FloatExtended) -> Self::Output {
        ComplexExtended::new(
            self.mantissa * other.mantissa,
            self.exponent + other.exponent
        )
    }
}

impl Mul<f64> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn mul(self, other: f64) -> Self::Output {
        ComplexExtended::new(
            self.mantissa * other,
            self.exponent
        )
    }
}

impl Div<f64> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn div(self, other: f64) -> Self::Output {
        ComplexExtended::new(
            self.mantissa / other,
            self.exponent
        )
    }
}

impl Display for ComplexExtended {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}*2^{}", self.mantissa, self.exponent)
    }
}