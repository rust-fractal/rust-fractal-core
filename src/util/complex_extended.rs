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

// TODO we can make a struct with shared exponent
impl ComplexExtended {
    #[inline]
    pub fn new(mantissa: Complex<f64>, exponent: i32) -> Self {
        let mut temp = ComplexExtended {
            mantissa,
            exponent
        };

        temp.reduce();
        temp
    }

    #[inline]
    pub fn new2(re: f64, im: f64, exponent: i32) -> Self {
        let mut temp = ComplexExtended {
            mantissa: Complex::<f64>::new(re, im),
            exponent
        };

        temp.reduce();
        temp
    }

    #[inline]
    pub fn norm(&self) -> FloatExtended {
        self.norm_square().sqrt()
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
        self.mantissa * 2.0f64.powi(self.exponent)
    }

    #[inline]
    pub fn reduce(&mut self) {
        if self.mantissa.re.abs() > self.mantissa.im.abs() {
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
}

impl AddAssign for ComplexExtended {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        if self.mantissa.re == 0.0 && self.mantissa.im == 0.0 {
            self.mantissa = other.mantissa;
            self.exponent = other.exponent;
        } else if other.mantissa.re == 0.0 && other.mantissa.im == 0.0 {
            return;
        } else {
            if self.exponent == other.exponent {
                self.mantissa += other.mantissa;
            } else if self.exponent > other.exponent {
                self.mantissa += other.mantissa / 2.0f64.powi(self.exponent - other.exponent);
            } else {
                self.mantissa /= 2.0f64.powi(other.exponent - self.exponent);
                self.exponent = other.exponent;
                self.mantissa += other.mantissa;
            }
        }
    }
}

impl SubAssign for ComplexExtended {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        if self.exponent == other.exponent {
            self.mantissa -= other.mantissa;
        } else if self.exponent > other.exponent {
            self.mantissa -= other.mantissa / 2.0f64.powi(self.exponent - other.exponent);
        } else {
            self.mantissa /= 2.0f64.powi(other.exponent - self.exponent);
            self.exponent = other.exponent;
            self.mantissa -= other.mantissa;
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
        self.reduce();
    }
}

impl Add<ComplexExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        if self.mantissa.re == 0.0 && self.mantissa.im == 0.0 {
            other
        } else if other.mantissa.re == 0.0 && other.mantissa.im == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == other.exponent {
                (self.mantissa + other.mantissa, self.exponent)
            } else if self.exponent > other.exponent {
                (self.mantissa + other.mantissa / 2.0f64.powi(self.exponent - other.exponent), self.exponent)
            } else {
                (other.mantissa + self.mantissa / 2.0f64.powi(other.exponent - self.exponent), other.exponent)
            };
            ComplexExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Sub<ComplexExtended> for ComplexExtended {
    type Output = ComplexExtended;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        if self.mantissa.re == 0.0 && self.mantissa.im == 0.0 {
            ComplexExtended::new(-other.mantissa, other.exponent)
        } else if other.mantissa.re == 0.0 && other.mantissa.im == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == other.exponent {
                (self.mantissa - other.mantissa, self.exponent)
            } else if self.exponent > other.exponent {
                (self.mantissa - other.mantissa / 2.0f64.powi(self.exponent - other.exponent), self.exponent)
            } else {
                (-1.0 * other.mantissa + self.mantissa / 2.0f64.powi(other.exponent - self.exponent), other.exponent)
            };
            ComplexExtended::new(new_mantissa, new_exponent)
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

// impl Add<f64> for ComplexExtended {
//     type Output = ComplexExtended;
//
//     fn add(self, other: f64) -> Self::Output {
//         if self.mantissa.re == 0.0 && self.mantissa.im == 0.0 {
//             ComplexExtended::new2(other, 0.0, 0)
//         } else if other == 0.0 {
//             self
//         } else {
//             let (new_mantissa, new_exponent) = if self.exponent == temp.exponent {
//                 (self.mantissa + temp.mantissa, self.exponent)
//             } else if self.exponent > temp.exponent {
//                 (self.mantissa + temp.mantissa / 2.0f64.powi(self.exponent - temp.exponent), self.exponent)
//             } else {
//                 (temp.mantissa + self.mantissa / 2.0f64.powi(temp.exponent - self.exponent), temp.exponent)
//             };
//             ComplexExtended::new(new_mantissa, new_exponent)
//         }
//     }
// }
//
// impl Sub<f64> for ComplexExtended {
//     type Output = ComplexExtended;
//
//     fn sub(self, other: f64) -> Self::Output {
//         if self.mantissa.re == 0.0 && self.mantissa.im == 0.0 {
//             ComplexExtended::new2(-other, 0.0, 0)
//         } else if other == 0.0 {
//             self
//         } else {
//             let (new_mantissa, new_exponent) = if self.exponent == other.exponent {
//                 (self.mantissa - other.mantissa, self.exponent)
//             } else if self.exponent > other.exponent {
//                 (self.mantissa - other.mantissa / 2.0f64.powi(self.exponent - other.exponent), self.exponent)
//             } else {
//                 (-1.0 * other.mantissa + self.mantissa / 2.0f64.powi(other.exponent - self.exponent), other.exponent)
//             };
//             ComplexExtended::new(new_mantissa, new_exponent)
//         }
//     }
// }

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