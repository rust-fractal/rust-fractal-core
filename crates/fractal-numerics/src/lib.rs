//! Numerical building blocks kept separate from the render domain and runtime.

mod error;
mod extended;

pub use error::MalachiteError;
pub use extended::Extended;

#[cfg(feature = "malachite")]
pub mod malachite {
    //! Arbitrary-precision operations backed exclusively by Malachite.

    use super::MalachiteError;
    pub use malachite::Float;
    pub use malachite::base::{
        num::conversion::traits::RoundingFrom, rounding_modes::RoundingMode,
    };

    #[derive(Clone, Debug)]
    pub struct MalachiteComplex {
        real: Float,
        imag: Float,
        precision_bits: u32,
    }

    impl MalachiteComplex {
        pub fn from_decimal(
            real: &str,
            imag: &str,
            precision_bits: u32,
        ) -> Result<Self, MalachiteError> {
            if precision_bits == 0 {
                return Err(MalachiteError::InvalidPrecision);
            }

            let precision = u64::from(precision_bits);
            let real = Float::from_sci_string_prec(real, precision)
                .map(|(value, _)| value)
                .ok_or_else(|| MalachiteError::Native(format!("invalid real value: {real}")))?;
            let imag =
                Float::from_sci_string_prec(imag, precision).map(|(value, _)| value).ok_or_else(
                    || MalachiteError::Native(format!("invalid imaginary value: {imag}")),
                )?;
            if !real.is_finite() || !imag.is_finite() {
                return Err(MalachiteError::Native(
                    "Malachite coordinates must be finite".to_owned(),
                ));
            }
            Ok(Self { real, imag, precision_bits })
        }

        pub fn from_f64(real: f64, imag: f64, precision_bits: u32) -> Result<Self, MalachiteError> {
            if precision_bits == 0 || !real.is_finite() || !imag.is_finite() {
                return Err(MalachiteError::InvalidPrecision);
            }
            let precision = u64::from(precision_bits);
            let (real, _) = Float::from_primitive_float_prec(real, precision);
            let (imag, _) = Float::from_primitive_float_prec(imag, precision);
            Ok(Self { real, imag, precision_bits })
        }

        pub fn zero(precision_bits: u32) -> Result<Self, MalachiteError> {
            Self::from_f64(0.0, 0.0, precision_bits)
        }

        pub const fn precision_bits(&self) -> u32 {
            self.precision_bits
        }

        pub fn real(&self) -> &Float {
            &self.real
        }

        pub fn imag(&self) -> &Float {
            &self.imag
        }

        pub fn is_finite(&self) -> bool {
            self.real.is_finite() && self.imag.is_finite()
        }

        pub fn square_add(&mut self, constant: &Self) -> Result<(), MalachiteError> {
            ensure_same_precision(self, constant)?;

            let real_square = &self.real * &self.real;
            let imag_square = &self.imag * &self.imag;
            let cross = &self.real * &self.imag;
            self.real = real_square - imag_square;
            self.real += &constant.real;
            self.imag = cross.clone() + cross;
            self.imag += &constant.imag;
            Ok(())
        }

        pub fn add(&self, other: &Self) -> Result<Self, MalachiteError> {
            ensure_same_precision(self, other)?;
            Ok(Self {
                real: &self.real + &other.real,
                imag: &self.imag + &other.imag,
                precision_bits: self.precision_bits,
            })
        }

        pub fn sub(&self, other: &Self) -> Result<Self, MalachiteError> {
            ensure_same_precision(self, other)?;
            Ok(Self {
                real: &self.real - &other.real,
                imag: &self.imag - &other.imag,
                precision_bits: self.precision_bits,
            })
        }

        pub fn multiply(&self, other: &Self) -> Result<Self, MalachiteError> {
            ensure_same_precision(self, other)?;
            Ok(Self {
                real: &self.real * &other.real - &self.imag * &other.imag,
                imag: &self.real * &other.imag + &self.imag * &other.real,
                precision_bits: self.precision_bits,
            })
        }

        pub fn divide(&self, other: &Self) -> Result<Self, MalachiteError> {
            ensure_same_precision(self, other)?;
            let denominator = &other.real * &other.real + &other.imag * &other.imag;
            if denominator.is_zero() {
                return Err(MalachiteError::Native("complex division by zero".to_owned()));
            }
            Ok(Self {
                real: (&self.real * &other.real + &self.imag * &other.imag) / &denominator,
                imag: (&self.imag * &other.real - &self.real * &other.imag) / &denominator,
                precision_bits: self.precision_bits,
            })
        }

        pub fn powu(&self, exponent: u32) -> Result<Self, MalachiteError> {
            let mut result = Self::from_f64(1.0, 0.0, self.precision_bits)?;
            let mut base = self.clone();
            let mut exponent = exponent;
            while exponent != 0 {
                if exponent & 1 == 1 {
                    result = result.multiply(&base)?;
                }
                exponent >>= 1;
                if exponent != 0 {
                    base = base.multiply(&base)?;
                }
            }
            Ok(result)
        }

        pub fn norm_squared(&self) -> Float {
            &self.real * &self.real + &self.imag * &self.imag
        }

        pub fn scale(&self, factor: &Float) -> Self {
            Self {
                real: &self.real * factor,
                imag: &self.imag * factor,
                precision_bits: self.precision_bits,
            }
        }

        pub fn is_escaped(&self) -> bool {
            self.norm_squared() > 4u32
        }

        pub fn to_extended(&self) -> (Option<super::Extended>, Option<super::Extended>) {
            (float_to_extended(&self.real), float_to_extended(&self.imag))
        }

        pub fn to_f64(&self) -> (f64, f64) {
            (
                f64::rounding_from(&self.real, RoundingMode::Nearest).0,
                f64::rounding_from(&self.imag, RoundingMode::Nearest).0,
            )
        }
    }

    fn ensure_same_precision(
        left: &MalachiteComplex,
        right: &MalachiteComplex,
    ) -> Result<(), MalachiteError> {
        if left.precision_bits == right.precision_bits {
            Ok(())
        } else {
            Err(MalachiteError::Native("Malachite values use different precisions".to_owned()))
        }
    }

    fn float_to_extended(value: &Float) -> Option<super::Extended> {
        let (mantissa, exponent, _) =
            value.sci_mantissa_and_exponent_round::<f64>(RoundingMode::Nearest)?;
        super::Extended::from_parts(mantissa * 0.5, exponent.checked_add(1)?)
    }
}
