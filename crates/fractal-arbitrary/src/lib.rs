//! Malachite-backed arbitrary-precision reference orbits.
//!
//! This crate intentionally does not contain a second numeric implementation.
//! All arbitrary-precision values and operations come from Malachite in
//! `fractal-numerics`; this crate owns only the
//! fractal-specific recurrence, checkpoint storage, and statistics.

use fractal_numerics::{
    MalachiteError,
    malachite::{Float, MalachiteComplex, RoundingFrom, RoundingMode},
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReferenceError {
    InvalidIterationLimit,
    InvalidZoom,
    Malachite(MalachiteError),
}

impl core::fmt::Display for ReferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidIterationLimit => f.write_str("iteration limit must be non-zero"),
            Self::InvalidZoom => f.write_str("reference zoom does not have a valid exponent"),
            Self::Malachite(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ReferenceError {}

impl From<MalachiteError> for ReferenceError {
    fn from(error: MalachiteError) -> Self {
        Self::Malachite(error)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceStatistics {
    pub precision_bits: u32,
    pub requested_iterations: u32,
    pub generated_iterations: u32,
    pub escaped_at: Option<u32>,
}

#[derive(Clone, Debug)]
pub struct MalachiteReference {
    pub center: MalachiteComplex,
    pub points: Vec<MalachiteComplex>,
    pub statistics: ReferenceStatistics,
}

impl MalachiteReference {
    /// Converts a scientific-notation location zoom into a conservative
    /// Malachite precision. The raw zoom remains a string throughout.
    pub fn precision_bits_for_zoom(zoom: &str) -> Result<u32, ReferenceError> {
        let exponent_start = zoom.bytes().position(|byte| byte == b'e' || byte == b'E');
        let exponent = exponent_start
            .map(|index| zoom.get(index + 1..).ok_or(ReferenceError::InvalidZoom))
            .transpose()?
            .map(|value| value.parse::<u64>().map_err(|_| ReferenceError::InvalidZoom))
            .transpose()?
            .unwrap_or(0);
        let precision = 256_u64
            .checked_add(exponent.saturating_mul(3322).div_ceil(1000))
            .ok_or(ReferenceError::InvalidZoom)?;
        u32::try_from(precision).map_err(|_| ReferenceError::InvalidZoom)
    }

    pub fn generate_for_zoom(
        real: &str,
        imag: &str,
        zoom: &str,
        iteration_limit: u32,
    ) -> Result<Self, ReferenceError> {
        let precision_bits = Self::precision_bits_for_zoom(zoom)?;
        Self::generate(real, imag, precision_bits, iteration_limit)
    }

    pub fn generate(
        real: &str,
        imag: &str,
        precision_bits: u32,
        iteration_limit: u32,
    ) -> Result<Self, ReferenceError> {
        if iteration_limit == 0 {
            return Err(ReferenceError::InvalidIterationLimit);
        }
        let center = MalachiteComplex::from_decimal(real, imag, precision_bits)?;
        let mut value = MalachiteComplex::zero(precision_bits)?;
        let mut points = Vec::with_capacity(iteration_limit as usize + 1);
        points.push(value.clone());
        let mut escaped_at = None;
        for iteration in 1..=iteration_limit {
            value.square_add(&center)?;
            points.push(value.clone());
            if value.is_escaped() {
                escaped_at = Some(iteration);
                break;
            }
        }
        Ok(Self {
            center,
            statistics: ReferenceStatistics {
                precision_bits,
                requested_iterations: iteration_limit,
                generated_iterations: points.len().saturating_sub(1) as u32,
                escaped_at,
            },
            points,
        })
    }

    pub fn point(&self, iteration: u32) -> Option<&MalachiteComplex> {
        self.points.get(iteration as usize)
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    pub fn maximum_norm_squared_f64(&self) -> f64 {
        self.points
            .iter()
            .map(|point| {
                let (real, imag) = point.to_f64();
                real.mul_add(real, imag * imag)
            })
            .fold(0.0, f64::max)
    }
}

pub fn compare_to_f64(reference: &MalachiteReference, expected: &[(f64, f64)]) -> Option<f64> {
    if expected.len() != reference.points.len() {
        return None;
    }
    let mut maximum_error: f64 = 0.0;
    for (point, &(expected_real, expected_imag)) in reference.points.iter().zip(expected) {
        let (real, imag) = point.to_f64();
        maximum_error = maximum_error.max((real - expected_real).abs());
        maximum_error = maximum_error.max((imag - expected_imag).abs());
    }
    Some(maximum_error)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RootError {
    InvalidPeriod,
    InvalidIterationLimit,
    UnsupportedFormula,
    Malachite(MalachiteError),
    DidNotConverge,
}

impl core::fmt::Display for RootError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidPeriod => f.write_str("root period must be non-zero"),
            Self::InvalidIterationLimit => f.write_str("root iteration limit must be non-zero"),
            Self::UnsupportedFormula => {
                f.write_str("Malachite root finding supports Mandelbrot powers only")
            }
            Self::Malachite(error) => error.fmt(f),
            Self::DidNotConverge => f.write_str("Malachite root finding did not converge"),
        }
    }
}

impl std::error::Error for RootError {}

impl From<MalachiteError> for RootError {
    fn from(error: MalachiteError) -> Self {
        Self::Malachite(error)
    }
}

#[derive(Clone, Debug)]
pub struct RootResult {
    pub center: MalachiteComplex,
    pub period: u32,
    pub iterations: u32,
    pub residual: f64,
}

/// Evaluates `f^period(0)` and its derivative with respect to `c` using
/// Malachite arithmetic. The returned residual is the value whose root is the
/// period nucleus.
pub fn period_residual(
    center: &MalachiteComplex,
    period: u32,
    power: u32,
) -> Result<(MalachiteComplex, MalachiteComplex), RootError> {
    if period == 0 {
        return Err(RootError::InvalidPeriod);
    }
    if power < 2 {
        return Err(RootError::UnsupportedFormula);
    }
    let precision = center.precision_bits();
    let one = MalachiteComplex::from_f64(1.0, 0.0, precision)?;
    let power_value = MalachiteComplex::from_f64(power as f64, 0.0, precision)?;
    let mut value = MalachiteComplex::zero(precision)?;
    let mut derivative = MalachiteComplex::zero(precision)?;
    for _ in 0..period {
        let powered = value.powu(power)?;
        let derivative_power = value.powu(power - 1)?;
        derivative = derivative_power.multiply(&derivative)?.multiply(&power_value)?.add(&one)?;
        value = powered.add(center)?;
    }
    Ok((value, derivative))
}

/// Finds a Mandelbrot nucleus with Newton iteration at the requested
/// Malachite precision. `tolerance` is parsed as a decimal at that same
/// precision, so termination does not depend on an `f64` conversion.
pub fn find_nucleus(
    initial: MalachiteComplex,
    period: u32,
    power: u32,
    maximum_iterations: u32,
    tolerance: &str,
) -> Result<RootResult, RootError> {
    if maximum_iterations == 0 {
        return Err(RootError::InvalidIterationLimit);
    }
    let precision = initial.precision_bits();
    let tolerance = Float::from_sci_string_prec(tolerance, u64::from(precision))
        .map(|(value, _)| value)
        .ok_or_else(|| {
            RootError::Malachite(MalachiteError::Native("invalid root tolerance".to_owned()))
        })?;
    if tolerance <= 0u32 {
        return Err(RootError::Malachite(MalachiteError::Native(
            "root tolerance must be positive".to_owned(),
        )));
    }
    let tolerance_squared = &tolerance * &tolerance;
    let mut center = initial;
    for iteration in 0..maximum_iterations {
        let (residual, derivative) = period_residual(&center, period, power)?;
        let residual_squared = residual.norm_squared();
        if residual_squared <= tolerance_squared {
            return Ok(RootResult {
                residual: f64::rounding_from(&residual_squared, RoundingMode::Nearest).0.sqrt(),
                center,
                period,
                iterations: iteration,
            });
        }
        center = center.sub(&residual.divide(&derivative)?)?;
    }
    Err(RootError::DidNotConverge)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn malachite_reference_matches_a_quadratic_orbit() {
        let reference = MalachiteReference::generate("-0.75", "0", 256, 32).unwrap();
        assert_eq!(reference.statistics.generated_iterations, 32);
        let mut expected = vec![(0.0, 0.0)];
        let mut z = (0.0, 0.0);
        for _ in 0..32 {
            z = (z.0 * z.0 - z.1 * z.1 - 0.75, 2.0 * z.0 * z.1);
            expected.push(z);
        }
        assert!(compare_to_f64(&reference, &expected).unwrap() < 1e-12);
    }

    #[test]
    fn arbitrary_reference_uses_decimal_coordinates_without_f64_parsing() {
        let reference = MalachiteReference::generate(
            "-1.999966194450370304184346885063505796755312415407248515117619229448015842423",
            "-0.0000000000000000000000000000000003001382436790938324072497303977594683",
            1024,
            64,
        )
        .unwrap();
        assert_eq!(reference.statistics.precision_bits, 1024);
        assert!(reference.points.len() > 1);
    }

    #[test]
    fn location_zoom_selects_malachite_precision_from_the_raw_exponent() {
        assert_eq!(MalachiteReference::precision_bits_for_zoom("1E8000"), Ok(26_832));
        let reference = MalachiteReference::generate_for_zoom("-0.75", "0", "1E8000", 8).unwrap();
        assert_eq!(reference.statistics.precision_bits, 26_832);
    }

    #[test]
    fn escaped_reference_stops_at_the_first_arbitrary_norm_crossing() {
        let reference = MalachiteReference::generate("1", "0", 128, 32).unwrap();
        assert_eq!(reference.statistics.escaped_at, Some(3));
        assert_eq!(reference.statistics.generated_iterations, 3);
    }

    #[test]
    fn malachite_period_one_nucleus_is_zero() {
        let initial = MalachiteComplex::from_f64(0.1, 0.1, 256).unwrap();
        let result = find_nucleus(initial, 1, 2, 16, "1e-60").unwrap();
        let (real, imag) = result.center.to_f64();
        assert!(real.abs() < 1e-50);
        assert!(imag.abs() < 1e-50);
        assert!(result.residual < 1e-50);
    }

    #[test]
    fn malachite_period_two_nucleus_is_minus_one() {
        let initial = MalachiteComplex::from_f64(-0.8, 0.1, 256).unwrap();
        let result = find_nucleus(initial, 2, 2, 32, "1e-50").unwrap();
        let (real, imag) = result.center.to_f64();
        assert!((real + 1.0).abs() < 1e-40);
        assert!(imag.abs() < 1e-40);
        assert!(result.residual < 1e-40);
    }
}
