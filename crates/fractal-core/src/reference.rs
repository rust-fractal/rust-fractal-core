use crate::{ComplexCoordinate, FractalFormula, RenderError};

#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceStatistics {
    pub requested_iterations: u32,
    pub generated_iterations: u32,
    pub escaped_at: Option<u32>,
    pub maximum_norm_squared: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceOrbit {
    pub center: ComplexCoordinate,
    pub formula: FractalFormula,
    pub points: Vec<ComplexCoordinate>,
    pub statistics: ReferenceStatistics,
}

impl ReferenceOrbit {
    pub fn point(&self, iteration: u32) -> Option<ComplexCoordinate> {
        self.points.get(iteration as usize).copied()
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

pub fn generate_reference(
    center: ComplexCoordinate,
    formula: FractalFormula,
    iteration_limit: u32,
) -> Result<ReferenceOrbit, RenderError> {
    formula.validate()?;
    if iteration_limit == 0 {
        return Err(RenderError::InvalidIterationLimit);
    }
    if !center.is_finite() {
        return Err(RenderError::NonFiniteCoordinate { real: center.real, imag: center.imag });
    }

    let mut points = Vec::with_capacity(iteration_limit as usize + 1);
    let mut z = ComplexCoordinate::new(0.0, 0.0);
    points.push(z);
    let mut escaped_at = None;
    let mut maximum_norm_squared: f64 = 0.0;
    for iteration in 1..=iteration_limit {
        z = step(z, center, formula);
        points.push(z);
        let norm_squared = z.norm_squared();
        maximum_norm_squared = maximum_norm_squared.max(norm_squared);
        if escaped_at.is_none() && norm_squared > 4.0 {
            escaped_at = Some(iteration);
            break;
        }
    }
    let generated_iterations = points.len().saturating_sub(1) as u32;
    Ok(ReferenceOrbit {
        center,
        formula,
        points,
        statistics: ReferenceStatistics {
            requested_iterations: iteration_limit,
            generated_iterations,
            escaped_at,
            maximum_norm_squared,
        },
    })
}

pub(crate) fn step(
    z: ComplexCoordinate,
    c: ComplexCoordinate,
    formula: FractalFormula,
) -> ComplexCoordinate {
    let (mut real, mut imag) = (z.real, z.imag);
    if matches!(formula, FractalFormula::BurningShip { .. }) {
        real = real.abs();
        imag = imag.abs();
    }
    let power = match formula {
        FractalFormula::Mandelbrot { power } | FractalFormula::BurningShip { power } => power,
    };
    let powered = complex_power(real, imag, power);
    ComplexCoordinate::new(powered.0 + c.real, powered.1 + c.imag)
}

pub(crate) fn complex_power(real: f64, imag: f64, power: u32) -> (f64, f64) {
    let mut result = (1.0, 0.0);
    let mut base = (real, imag);
    let mut exponent = power;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = complex_multiply(result, base);
        }
        exponent >>= 1;
        if exponent != 0 {
            base = complex_multiply(base, base);
        }
    }
    result
}

pub(crate) fn complex_multiply(left: (f64, f64), right: (f64, f64)) -> (f64, f64) {
    (left.0 * right.0 - left.1 * right.1, left.0 * right.1 + left.1 * right.0)
}

impl ComplexCoordinate {
    pub(crate) fn norm_squared(self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }

    pub(crate) fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }

    pub(crate) fn add(self, other: Self) -> Self {
        Self::new(self.real + other.real, self.imag + other.imag)
    }

    pub(crate) fn sub(self, other: Self) -> Self {
        Self::new(self.real - other.real, self.imag - other.imag)
    }

    pub(crate) fn scale(self, factor: f64) -> Self {
        Self::new(self.real * factor, self.imag * factor)
    }

    pub(crate) fn multiply(self, other: Self) -> Self {
        let product = complex_multiply((self.real, self.imag), (other.real, other.imag));
        Self::new(product.0, product.1)
    }

    pub(crate) fn divide(self, other: Self) -> Option<Self> {
        let denominator = other.norm_squared();
        if denominator == 0.0 || !denominator.is_finite() {
            return None;
        }
        Some(Self::new(
            (self.real * other.real + self.imag * other.imag) / denominator,
            (self.imag * other.real - self.real * other.imag) / denominator,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbit_has_explicit_zero_based_iteration_indexing() {
        let orbit =
            generate_reference(ComplexCoordinate::new(0.0, 0.0), FractalFormula::default(), 8)
                .unwrap();
        assert_eq!(orbit.point(0), Some(ComplexCoordinate::new(0.0, 0.0)));
        assert_eq!(orbit.point(8), Some(ComplexCoordinate::new(0.0, 0.0)));
        assert_eq!(orbit.statistics.generated_iterations, 8);
    }

    #[test]
    fn exterior_orbit_records_escape_iteration() {
        let orbit =
            generate_reference(ComplexCoordinate::new(1.0, 0.0), FractalFormula::default(), 32)
                .unwrap();
        assert_eq!(orbit.statistics.escaped_at, Some(3));
        assert_eq!(orbit.statistics.generated_iterations, 3);
    }
}
