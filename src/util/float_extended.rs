use crate::util::FloatExp;
use std::ops::{AddAssign, Mul, Add, Sub, SubAssign, MulAssign, DivAssign, Div};
use std::cmp::Ordering;
use std::fmt::{Display, Formatter, Result};

#[derive(Debug, Copy, Clone)]
pub struct FloatExtended {
    pub mantissa: f64,
    pub exponent: i32,
}

impl FloatExtended {
    #[inline]
    pub fn new(mantissa: f64, exponent: i32) -> Self {
        let mut output = FloatExtended {
            mantissa,
            exponent
        };
        output.reduce();
        output
    }

    #[inline]
    pub fn reduce(&mut self) {
        // Check to make sure that the mantissa is outside of the new bounds
        if self.mantissa < 0.5 || self.mantissa > 1.0 {
            let (temp_mantissa, added_exponent) = self.mantissa.frexp();
            self.mantissa = temp_mantissa;
            self.exponent += added_exponent;
        }
    }

    #[inline]
    pub fn sqrt(&self) -> FloatExtended {
        let (new_mantissa, new_exponent) = if self.exponent % 2 == 0 {
            (self.mantissa.sqrt(), self.exponent / 2)
        } else {
            ((2.0 * self.mantissa).sqrt(), (self.exponent - 1) / 2)
        };
        FloatExtended::new(new_mantissa, new_exponent)
    }

    #[inline]
    pub fn to_float(&self) -> f64 {
        self.mantissa.ldexp(self.exponent)
    }
}

impl PartialEq for FloatExtended {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.mantissa == other.mantissa && self.exponent == other.exponent
    }
}

impl PartialOrd for FloatExtended {
    // TODO sort out NaN values
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // This is horrible
        if self.mantissa == 0.0 && other.mantissa < 0.0 {
            Some(Ordering::Greater)
        } else if self.mantissa == 0.0 && other.mantissa > 0.0 {
            Some(Ordering::Less)
        } else if other.mantissa == 0.0 && self.mantissa < 0.0 {
            Some(Ordering::Less)
        } else if other.mantissa == 0.0 && self.mantissa > 0.0 {
            Some(Ordering::Greater)
        } else if self.exponent < other.exponent {
            Some(Ordering::Less)
        } else if self.exponent > other.exponent {
            Some(Ordering::Greater)
        } else {
            if self.mantissa == other.mantissa {
                Some(Ordering::Equal)
            } else if self.mantissa < other.mantissa {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            }
        }
    }
}

impl AddAssign for FloatExtended {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        if self.exponent == other.exponent {
            self.mantissa += other.mantissa;
        } else if self.exponent > other.exponent {
            self.mantissa += other.mantissa / 2.0f64.powi(self.exponent - other.exponent);
        } else {
            self.mantissa /= 2.0f64.powi(other.exponent - self.exponent);
            self.exponent = other.exponent;
            self.mantissa += other.mantissa;
        }
        self.reduce();
    }
}

impl SubAssign for FloatExtended {
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
        self.reduce();
    }
}

impl MulAssign<FloatExtended> for FloatExtended {
    #[inline]
    fn mul_assign(&mut self, other: Self) {
        self.mantissa *= other.mantissa;
        self.exponent += other.exponent;
        self.reduce();
    }
}

impl MulAssign<f64> for FloatExtended {
    #[inline]
    fn mul_assign(&mut self, other: f64) {
        self.mantissa *= other;
        self.reduce();
    }
}

impl DivAssign for FloatExtended {
    #[inline]
    fn div_assign(&mut self, other: Self) {
        self.mantissa /= other.mantissa;
        self.exponent -= other.exponent;
        self.reduce();
    }
}

impl Add<FloatExtended> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        if self.mantissa == 0.0 {
            other
        } else if other.mantissa == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == other.exponent {
                (self.mantissa + other.mantissa, self.exponent)
            } else if self.exponent > other.exponent {
                (self.mantissa + other.mantissa / 2.0f64.powi(self.exponent - other.exponent), self.exponent)
            } else {
                (other.mantissa + self.mantissa / 2.0f64.powi(other.exponent - self.exponent), other.exponent)
            };
            FloatExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Sub<FloatExtended> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        if self.mantissa == 0.0 {
            -1.0 * other
        } else if other.mantissa == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == other.exponent {
                (self.mantissa - other.mantissa, self.exponent)
            } else if self.exponent > other.exponent {
                (self.mantissa - other.mantissa / 2.0f64.powi(self.exponent - other.exponent), self.exponent)
            } else {
                (-1.0 * other.mantissa + self.mantissa / 2.0f64.powi(other.exponent - self.exponent), other.exponent)
            };
            FloatExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Mul<FloatExtended> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn mul(self, other: Self) -> Self::Output {
        FloatExtended::new(
            self.mantissa * other.mantissa,
            self.exponent + other.exponent
        )
    }
}

impl Div<FloatExtended> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn div(self, other: Self) -> Self::Output {
        FloatExtended::new(
            self.mantissa / other.mantissa,
            self.exponent - other.exponent
        )
    }
}

impl Add<f64> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn add(self, other: f64) -> Self::Output {
        let temp = FloatExtended::new(other, 0);

        if self.mantissa == 0.0 {
            temp
        } else if temp.mantissa == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == temp.exponent {
                (self.mantissa + temp.mantissa, self.exponent)
            } else if self.exponent > temp.exponent {
                (self.mantissa + temp.mantissa / 2.0f64.powi(self.exponent - temp.exponent), self.exponent)
            } else {
                (temp.mantissa + self.mantissa / 2.0f64.powi(temp.exponent - self.exponent), temp.exponent)
            };
            FloatExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Sub<f64> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn sub(self, other: f64) -> Self::Output {
        let temp = FloatExtended::new(other, 0);

        if self.mantissa == 0.0 {
            -1.0 * temp
        } else if temp.mantissa == 0.0 {
            self
        } else {
            let (new_mantissa, new_exponent) = if self.exponent == temp.exponent {
                (self.mantissa - temp.mantissa, self.exponent)
            } else if self.exponent > temp.exponent {
                (self.mantissa - temp.mantissa / 2.0f64.powi(self.exponent - temp.exponent), self.exponent)
            } else {
                (-1.0 * temp.mantissa + self.mantissa / 2.0f64.powi(temp.exponent - self.exponent), temp.exponent)
            };
            FloatExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Mul<f64> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn mul(self, other: f64) -> Self::Output {
        let temp = FloatExtended::new(other, 0);

        FloatExtended::new(
            self.mantissa * temp.mantissa,
            self.exponent + temp.exponent
        )
    }
}

impl Div<f64> for FloatExtended {
    type Output = FloatExtended;

    #[inline]
    fn div(self, other: f64) -> Self::Output {
        let temp = FloatExtended::new(other, 0);

        FloatExtended::new(
            self.mantissa / temp.mantissa,
            self.exponent - temp.exponent
        )
    }
}

impl Add<FloatExtended> for f64 {
    type Output = FloatExtended;

    #[inline]
    fn add(self, other: FloatExtended) -> Self::Output {
        let temp = FloatExtended::new(self, 0);

        if temp.mantissa == 0.0 {
            other
        } else if other.mantissa == 0.0 {
            temp
        } else {
            let (new_mantissa, new_exponent) = if temp.exponent == other.exponent {
                (temp.mantissa + other.mantissa, temp.exponent)
            } else if temp.exponent > other.exponent {
                (temp.mantissa + other.mantissa / 2.0f64.powi(temp.exponent - other.exponent), temp.exponent)
            } else {
                (other.mantissa + temp.mantissa / 2.0f64.powi(other.exponent - temp.exponent), other.exponent)
            };
            FloatExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Sub<FloatExtended> for f64 {
    type Output = FloatExtended;

    #[inline]
    fn sub(self, other: FloatExtended) -> Self::Output {
        let temp = FloatExtended::new(self, 0);

        if temp.mantissa == 0.0 {
            -1.0 * other
        } else if other.mantissa == 0.0 {
            temp
        } else {
            let (new_mantissa, new_exponent) = if temp.exponent == other.exponent {
                (temp.mantissa - other.mantissa, temp.exponent)
            } else if temp.exponent > other.exponent {
                (temp.mantissa - other.mantissa / 2.0f64.powi(temp.exponent - other.exponent), temp.exponent)
            } else {
                (-1.0 * other.mantissa + temp.mantissa / 2.0f64.powi(other.exponent - temp.exponent), other.exponent)
            };
            FloatExtended::new(new_mantissa, new_exponent)
        }
    }
}

impl Mul<FloatExtended> for f64 {
    type Output = FloatExtended;

    #[inline]
    fn mul(self, other: FloatExtended) -> Self::Output {
        let temp = FloatExtended::new(self, 0);

        FloatExtended::new(
            temp.mantissa * other.mantissa,
            temp.exponent + other.exponent
        )
    }
}

impl Div<FloatExtended> for f64 {
    type Output = FloatExtended;

    #[inline]
    fn div(self, other: FloatExtended) -> Self::Output {
        let temp = FloatExtended::new(self, 0);

        FloatExtended::new(
            temp.mantissa / other.mantissa,
            temp.exponent - other.exponent
        )
    }
}

impl Display for FloatExtended {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "{}*2^{}", self.mantissa, self.exponent)
    }
}