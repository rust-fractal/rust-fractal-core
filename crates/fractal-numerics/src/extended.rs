use core::{cmp::Ordering, fmt, str::FromStr};

/// A finite binary scientific value represented as `mantissa * 2^exponent`.
///
/// The mantissa is normalized to `[0.5, 1)` for non-zero values. This avoids
/// platform C floating-point decomposition APIs while retaining the useful
/// exponent range needed by deep-zoom coordinates and distance estimates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extended {
    mantissa: f64,
    exponent: i32,
}

impl Extended {
    pub const ZERO: Self = Self { mantissa: 0.0, exponent: 0 };
    pub const FOUR: Self = Self { mantissa: 0.5, exponent: 3 };

    pub fn from_f64(value: f64) -> Option<Self> {
        if !value.is_finite() {
            return None;
        }
        if value == 0.0 {
            return Some(Self::ZERO);
        }

        let sign = value.is_sign_negative();
        let absolute = value.abs();
        let bits = absolute.to_bits();
        let exponent_bits = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        if exponent_bits == 0 {
            return Self::from_f64(absolute * 2_f64.powi(54)).map(|mut result| {
                result.exponent -= 54;
                if sign {
                    result.mantissa = -result.mantissa;
                }
                result
            });
        }

        let mantissa_bits = (1022_u64 << 52) | fraction;
        let mut mantissa = f64::from_bits(mantissa_bits);
        if sign {
            mantissa = -mantissa;
        }
        Some(Self { mantissa, exponent: exponent_bits - 1022 })
    }

    pub const fn from_parts(mantissa: f64, exponent: i32) -> Option<Self> {
        if mantissa == 0.0 {
            return Some(Self::ZERO);
        }
        if !mantissa.is_finite() || mantissa.abs() < 0.5 || mantissa.abs() >= 1.0 {
            return None;
        }
        Some(Self { mantissa, exponent })
    }

    pub const fn mantissa(self) -> f64 {
        self.mantissa
    }

    pub const fn exponent(self) -> i32 {
        self.exponent
    }

    pub fn to_f64(self) -> f64 {
        self.mantissa * 2_f64.powi(self.exponent)
    }

    pub fn checked_mul(self, other: Self) -> Option<Self> {
        if self.mantissa == 0.0 || other.mantissa == 0.0 {
            return Some(Self::ZERO);
        }
        let mantissa = self.mantissa * other.mantissa;
        let exponent = self.exponent.checked_add(other.exponent)?;
        normalize(mantissa, exponent)
    }

    pub fn checked_add(self, other: Self) -> Option<Self> {
        if self.mantissa == 0.0 {
            return Some(other);
        }
        if other.mantissa == 0.0 {
            return Some(self);
        }

        let (larger, smaller) =
            if self.exponent >= other.exponent { (self, other) } else { (other, self) };
        let difference = larger.exponent.checked_sub(smaller.exponent)?;
        if difference > 1074 {
            return Some(larger);
        }
        let mantissa = larger.mantissa + smaller.mantissa * 2_f64.powi(-difference);
        normalize(mantissa, larger.exponent)
    }

    pub fn checked_sub(self, other: Self) -> Option<Self> {
        Self { mantissa: -other.mantissa, exponent: other.exponent }.checked_add(self)
    }

    pub fn checked_div(self, other: Self) -> Option<Self> {
        if other.mantissa == 0.0 {
            return None;
        }
        let mantissa = self.mantissa / other.mantissa;
        let exponent = self.exponent.checked_sub(other.exponent)?;
        normalize(mantissa, exponent)
    }

    pub fn checked_pow(self, exponent: u32) -> Option<Self> {
        let mut result = Self::from_f64(1.0)?;
        let mut base = self;
        let mut exponent = exponent;
        while exponent != 0 {
            if exponent & 1 == 1 {
                result = result.checked_mul(base)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                base = base.checked_mul(base)?;
            }
        }
        Some(result)
    }

    pub fn abs(self) -> Self {
        Self { mantissa: self.mantissa.abs(), exponent: self.exponent }
    }

    /// Converts a finite decimal/scientific value into the compact binary
    /// representation. Only the mantissa is rounded to f64; the exponent is
    /// retained separately, so values such as `1E8000` do not underflow.
    pub fn from_decimal(value: &str) -> Result<Self, &'static str> {
        let value = value.trim();
        let (significand, exponent) =
            value.split_once(['e', 'E']).map_or((value, 0_i64), |(mantissa, exponent)| {
                (mantissa, exponent.parse::<i64>().unwrap_or(i64::MIN))
            });
        if exponent == i64::MIN {
            return Err("invalid decimal exponent");
        }
        let negative = significand.starts_with('-');
        let unsigned = significand.trim_start_matches(['+', '-']);
        let fractional_digits = unsigned.split_once('.').map_or(0, |(_, fraction)| fraction.len());
        let digits = unsigned.replace('.', "");
        let digits = digits.trim_start_matches('0');
        if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
            return Err("invalid decimal significand");
        }
        let prefix_length = digits.len().min(17);
        let prefix =
            digits[..prefix_length].parse::<f64>().map_err(|_| "invalid decimal significand")?;
        let remaining_digits = digits.len() - prefix_length;
        let decimal_exponent = exponent
            .checked_sub(i64::try_from(fractional_digits).map_err(|_| "decimal too long")?)
            .and_then(|exponent| exponent.checked_sub(i64::try_from(remaining_digits).ok()?))
            .and_then(|exponent| exponent.checked_add(i64::try_from(prefix_length).ok()? - 1))
            .ok_or("decimal exponent overflow")?;
        let logarithm = prefix.log2() + (decimal_exponent as f64) * 10_f64.log2();
        if !logarithm.is_finite() {
            return Err("decimal exponent is outside the extended range");
        }
        let exponent = logarithm.floor();
        let exponent = i32::try_from(exponent as i64).map_err(|_| "extended exponent overflow")?;
        let mut mantissa = 2_f64.powf(logarithm - f64::from(exponent));
        if negative {
            mantissa = -mantissa;
        }
        Self::from_parts(
            mantissa * 0.5,
            exponent.checked_add(1).ok_or("extended exponent overflow")?,
        )
        .ok_or("decimal value is not finite")
    }
}

fn normalize(mantissa: f64, exponent: i32) -> Option<Extended> {
    if mantissa == 0.0 {
        return Some(Extended::ZERO);
    }
    if !mantissa.is_finite() {
        return None;
    }
    let bits = mantissa.abs().to_bits();
    let exponent_bits = ((bits >> 52) & 0x7ff) as i32;
    if exponent_bits == 0 {
        return Extended::from_f64(mantissa * 2_f64.powi(54)).and_then(|value| {
            Some(Extended {
                mantissa: value.mantissa,
                exponent: value.exponent.checked_add(exponent)?,
            })
        });
    }
    let scale = exponent_bits - 1022;
    let normalized = mantissa * 2_f64.powi(-scale);
    Some(Extended { mantissa: normalized, exponent: exponent.checked_add(scale)? })
}

/// A complex value whose real and imaginary components retain independent
/// binary exponents. This is the compact representation used by perturbation;
/// it never combines a deep pixel offset with its absolute reference center.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexExtended {
    real: Extended,
    imag: Extended,
}

impl ComplexExtended {
    pub const ZERO: Self = Self { real: Extended::ZERO, imag: Extended::ZERO };

    pub const fn from_components(real: Extended, imag: Extended) -> Self {
        Self { real, imag }
    }

    pub fn from_f64(real: f64, imag: f64) -> Option<Self> {
        Some(Self { real: Extended::from_f64(real)?, imag: Extended::from_f64(imag)? })
    }

    pub fn from_decimal(real: &str, imag: &str) -> Result<Self, &'static str> {
        Ok(Self { real: Extended::from_decimal(real)?, imag: Extended::from_decimal(imag)? })
    }

    pub const fn real(self) -> Extended {
        self.real
    }

    pub const fn imag(self) -> Extended {
        self.imag
    }

    #[allow(clippy::should_implement_trait)]
    pub fn add(self, other: Self) -> Option<Self> {
        Some(Self {
            real: self.real.checked_add(other.real)?,
            imag: self.imag.checked_add(other.imag)?,
        })
    }

    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, other: Self) -> Option<Self> {
        Some(Self {
            real: self.real.checked_sub(other.real)?,
            imag: self.imag.checked_sub(other.imag)?,
        })
    }

    pub fn multiply(self, other: Self) -> Option<Self> {
        let real =
            self.real.checked_mul(other.real)?.checked_sub(self.imag.checked_mul(other.imag)?)?;
        let imag =
            self.real.checked_mul(other.imag)?.checked_add(self.imag.checked_mul(other.real)?)?;
        Some(Self { real, imag })
    }

    pub fn scale(self, factor: f64) -> Option<Self> {
        Some(Self {
            real: self.real.checked_mul(Extended::from_f64(factor)?)?,
            imag: self.imag.checked_mul(Extended::from_f64(factor)?)?,
        })
    }

    pub fn norm_squared(self) -> Option<Extended> {
        self.real.checked_mul(self.real)?.checked_add(self.imag.checked_mul(self.imag)?)
    }

    pub fn abs_components(self) -> Self {
        Self { real: self.real.abs(), imag: self.imag.abs() }
    }

    pub fn to_f64(self) -> (f64, f64) {
        (self.real.to_f64(), self.imag.to_f64())
    }
}

impl fmt::Display for Extended {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} * 2^{}", self.mantissa, self.exponent)
    }
}

impl FromStr for Extended {
    type Err = &'static str;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let (mantissa, exponent) = value.split_once("* 2^").ok_or("invalid extended format")?;
        let mantissa = mantissa.trim().parse::<f64>().map_err(|_| "invalid mantissa")?;
        let exponent = exponent.trim().parse::<i32>().map_err(|_| "invalid exponent")?;
        Self::from_parts(mantissa, exponent).ok_or("value is not normalized")
    }
}

impl PartialOrd for Extended {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.mantissa == 0.0 && other.mantissa == 0.0 {
            return Some(Ordering::Equal);
        }
        if self.mantissa.is_sign_negative() != other.mantissa.is_sign_negative() {
            return self.mantissa.partial_cmp(&other.mantissa);
        }
        let magnitude = self
            .exponent
            .cmp(&other.exponent)
            .then_with(|| self.mantissa.abs().partial_cmp(&other.mantissa.abs()).unwrap());
        if self.mantissa.is_sign_negative() { Some(magnitude.reverse()) } else { Some(magnitude) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn round_trip_preserves_normal_floats() {
        for value in [1.0, -0.5, 3.125, 1e-200, 1e200] {
            let extended = Extended::from_f64(value).unwrap();
            assert_eq!(extended.to_f64(), value);
        }
    }

    #[test]
    fn multiplication_tracks_exponent_separately() {
        let left = Extended::from_f64(1e200).unwrap();
        let right = Extended::from_f64(1e200).unwrap();
        let product = left.checked_mul(right).unwrap();
        assert!(product.exponent() > left.exponent());
        assert!(product.to_f64().is_infinite());
    }

    #[test]
    fn large_exponent_addition_keeps_the_dominant_value() {
        let large = Extended::from_parts(0.5, 2_000).unwrap();
        let small = Extended::from_parts(0.5, 0).unwrap();
        assert_eq!(large.checked_add(small), Some(large));
    }

    #[test]
    fn subtraction_division_and_powers_preserve_values() {
        let value = Extended::from_f64(3.0).unwrap();
        assert!(
            (value.checked_sub(Extended::from_f64(1.0).unwrap()).unwrap().to_f64() - 2.0).abs()
                < 1e-12
        );
        assert!(
            (value.checked_div(Extended::from_f64(2.0).unwrap()).unwrap().to_f64() - 1.5).abs()
                < 1e-12
        );
        assert!((value.checked_pow(3).unwrap().to_f64() - 27.0).abs() < 1e-12);
    }

    #[test]
    fn ordering_and_format_round_trip_are_stable() {
        let value = Extended::from_f64(-3.125).unwrap();
        let formatted = value.to_string();
        assert_eq!(formatted.parse::<Extended>().unwrap(), value);
        assert!(Extended::from_f64(-4.0).unwrap() < Extended::from_f64(2.0).unwrap());
    }

    proptest! {
        #[test]
        fn finite_floats_round_trip(value in -1.0e200_f64..1.0e200) {
            prop_assume!(value.is_finite());
            let extended = Extended::from_f64(value).unwrap();
            prop_assert_eq!(extended.to_f64(), value);
            prop_assert!(extended.mantissa().abs() >= 0.5 || value == 0.0);
            prop_assert!(extended.mantissa().abs() < 1.0 || value == 0.0);
        }
    }
}
