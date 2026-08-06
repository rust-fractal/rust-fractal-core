use crate::{ComplexCoordinate, FractalFormula, RenderError, reference};

#[derive(Clone, Debug, PartialEq)]
pub enum RootError {
    InvalidPeriod,
    SingularDerivative,
    DidNotConverge { residual: f64 },
    InvalidRender(RenderError),
}

impl core::fmt::Display for RootError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidPeriod => f.write_str("root period must be non-zero"),
            Self::SingularDerivative => f.write_str("root derivative is singular"),
            Self::DidNotConverge { residual } => {
                write!(f, "root finding did not converge: residual {residual}")
            }
            Self::InvalidRender(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for RootError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RootResult {
    pub center: ComplexCoordinate,
    pub period: u32,
    pub iterations: u32,
    pub residual: f64,
}

pub fn period_residual(
    center: ComplexCoordinate,
    formula: FractalFormula,
    period: u32,
) -> Result<(ComplexCoordinate, ComplexCoordinate), RootError> {
    if period == 0 {
        return Err(RootError::InvalidPeriod);
    }
    formula.validate().map_err(RootError::InvalidRender)?;
    let power = match formula {
        FractalFormula::Mandelbrot { power } | FractalFormula::BurningShip { power } => power,
    };
    let mut z = ComplexCoordinate::new(0.0, 0.0);
    let mut derivative = ComplexCoordinate::new(0.0, 0.0);
    for _ in 0..period {
        let powered = reference::complex_power(z.real, z.imag, power);
        let derivative_power = reference::complex_power(z.real, z.imag, power - 1);
        derivative = derivative
            .multiply(ComplexCoordinate::new(
                derivative_power.0 * power as f64,
                derivative_power.1 * power as f64,
            ))
            .add(ComplexCoordinate::new(1.0, 0.0));
        z = ComplexCoordinate::new(powered.0 + center.real, powered.1 + center.imag);
    }
    Ok((z, derivative))
}

pub fn find_nucleus(
    initial: ComplexCoordinate,
    formula: FractalFormula,
    period: u32,
    maximum_iterations: u32,
    tolerance: f64,
) -> Result<RootResult, RootError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(RootError::InvalidRender(RenderError::InvalidGlitchTolerance { tolerance }));
    }
    let mut center = initial;
    for iteration in 1..=maximum_iterations {
        let (residual, derivative) = period_residual(center, formula, period)?;
        let residual_norm = residual.norm();
        if residual_norm <= tolerance {
            return Ok(RootResult {
                center,
                period,
                iterations: iteration - 1,
                residual: residual_norm,
            });
        }
        center = center.sub(residual.divide(derivative).ok_or(RootError::SingularDerivative)?);
    }
    let (residual, _) = period_residual(center, formula, period)?;
    Err(RootError::DidNotConverge { residual: residual.norm() })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn period_one_mandelbrot_nucleus_is_zero() {
        let result =
            find_nucleus(ComplexCoordinate::new(0.1, 0.1), FractalFormula::default(), 1, 20, 1e-12)
                .unwrap();
        assert!(result.center.norm() < 1e-10);
        assert!(result.residual < 1e-12);
    }

    #[test]
    fn period_two_mandelbrot_nucleus_is_minus_one() {
        let result = find_nucleus(
            ComplexCoordinate::new(-0.8, 0.1),
            FractalFormula::default(),
            2,
            30,
            1e-10,
        )
        .unwrap();
        assert!((result.center.real + 1.0).abs() < 1e-8);
        assert!(result.center.imag.abs() < 1e-8);
    }
}
