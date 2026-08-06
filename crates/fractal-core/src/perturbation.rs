use crate::{
    ComplexCoordinate, FractalFormula, ReferenceOrbit, RenderError, SeriesApproximation, reference,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlitchDetector {
    tolerance: f64,
}

impl GlitchDetector {
    pub fn new(tolerance: f64) -> Result<Self, RenderError> {
        if tolerance.is_finite() && tolerance > 0.0 {
            Ok(Self { tolerance })
        } else {
            Err(RenderError::InvalidGlitchTolerance { tolerance })
        }
    }

    pub fn is_glitch(&self, reference: ComplexCoordinate, delta: ComplexCoordinate) -> bool {
        let reference_norm = reference.norm_squared().max(f64::MIN_POSITIVE);
        delta.norm_squared() > self.tolerance * reference_norm
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlitchResult {
    pub detected: bool,
    pub maximum_delta_norm: f64,
    pub first_iteration: Option<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PerturbationResult {
    pub iterations: u32,
    pub escaped: bool,
    pub smooth_iteration: f64,
    pub distance: f64,
    pub glitch: GlitchResult,
}

pub fn perturbate(
    reference_orbit: &ReferenceOrbit,
    pixel: ComplexCoordinate,
    iteration_limit: u32,
    glitch_detector: Option<GlitchDetector>,
) -> Result<PerturbationResult, RenderError> {
    perturbate_with_series(reference_orbit, pixel, iteration_limit, glitch_detector, None)
}

pub fn perturbate_with_series(
    reference_orbit: &ReferenceOrbit,
    pixel: ComplexCoordinate,
    iteration_limit: u32,
    glitch_detector: Option<GlitchDetector>,
    series: Option<&SeriesApproximation>,
) -> Result<PerturbationResult, RenderError> {
    if iteration_limit == 0 {
        return Err(RenderError::InvalidIterationLimit);
    }
    if !pixel.is_finite() {
        return Err(RenderError::NonFiniteCoordinate { real: pixel.real, imag: pixel.imag });
    }
    let delta_c = pixel.sub(reference_orbit.center);
    let mut delta = ComplexCoordinate::new(0.0, 0.0);
    let mut derivative = ComplexCoordinate::new(0.0, 0.0);
    let mut maximum_delta_norm: f64 = 0.0;
    let mut first_glitch = None;
    let mut z = ComplexCoordinate::new(0.0, 0.0);
    let power = match reference_orbit.formula {
        FractalFormula::Mandelbrot { power } | FractalFormula::BurningShip { power } => power,
    };
    let mut iteration = 0;
    if let Some(series) = series
        && let Some(skip) = series.largest_valid_iteration(iteration_limit, delta_c, 1e-3)
    {
        delta = series.evaluate(skip, delta_c).ok_or(RenderError::InvalidRawFrame)?;
        if let Some(reference_z) = reference_orbit.point(skip) {
            z = reference_z.add(delta);
            iteration = skip;
        }
    }
    while iteration < iteration_limit {
        let reference_z = reference_orbit.point(iteration).unwrap_or(z);
        let reference_next = reference_orbit.point(iteration + 1);
        if matches!(reference_orbit.formula, FractalFormula::BurningShip { .. }) {
            z = reference::step(z, pixel, reference_orbit.formula);
            delta = z.sub(reference_next.unwrap_or(z));
        } else if let Some(reference_next) = reference_next {
            delta = polynomial_delta(reference_z, delta, delta_c, power);
            z = reference_next.add(delta);
        } else {
            z = reference::step(z, pixel, reference_orbit.formula);
        }

        if first_glitch.is_none()
            && glitch_detector.is_some_and(|detector| detector.is_glitch(reference_z, delta))
        {
            first_glitch = Some(iteration);
        }
        maximum_delta_norm = maximum_delta_norm.max(delta.norm());
        iteration += 1;
        if z.norm_squared() > 4.0 || !z.is_finite() {
            break;
        }
        if matches!(reference_orbit.formula, FractalFormula::Mandelbrot { .. }) {
            let derivative_power = reference::complex_power(
                reference_z.real,
                reference_z.imag,
                power.saturating_sub(1),
            );
            derivative = derivative
                .multiply(ComplexCoordinate::new(
                    derivative_power.0 * power as f64,
                    derivative_power.1 * power as f64,
                ))
                .add(ComplexCoordinate::new(1.0, 0.0));
        }
    }
    let escaped = z.norm_squared() > 4.0;
    let smooth_iteration = if escaped && z.norm() > 1.0 {
        iteration as f64 + 1.0 - z.norm().ln().ln() / (power as f64).ln()
    } else {
        iteration as f64
    };
    let derivative_norm = derivative.norm();
    let distance = if escaped && derivative_norm > 0.0 {
        0.5 * z.norm().ln() * z.norm() / derivative_norm
    } else {
        0.0
    };
    Ok(PerturbationResult {
        iterations: iteration,
        escaped,
        smooth_iteration,
        distance,
        glitch: GlitchResult {
            detected: first_glitch.is_some(),
            maximum_delta_norm,
            first_iteration: first_glitch,
        },
    })
}

fn polynomial_delta(
    reference_z: ComplexCoordinate,
    delta: ComplexCoordinate,
    delta_c: ComplexCoordinate,
    power: u32,
) -> ComplexCoordinate {
    let mut result = delta_c;
    let mut delta_power = ComplexCoordinate::new(1.0, 0.0);
    for order in 1..=power {
        delta_power = delta_power.multiply(delta);
        let coefficient = binomial(power, order) as f64;
        let reference_term =
            reference::complex_power(reference_z.real, reference_z.imag, power - order);
        result = result.add(
            ComplexCoordinate::new(reference_term.0, reference_term.1)
                .multiply(delta_power)
                .scale(coefficient),
        );
    }
    result
}

fn binomial(n: u32, k: u32) -> u64 {
    let mut result = 1_u64;
    for index in 1..=k {
        result = result * u64::from(n + 1 - index) / u64::from(index);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate_reference;

    #[test]
    fn perturbation_matches_the_central_reference() {
        let center = ComplexCoordinate::new(-0.75, 0.0);
        let reference = generate_reference(center, FractalFormula::default(), 100).unwrap();
        let result =
            perturbate(&reference, center, 100, Some(GlitchDetector::new(1e-6).unwrap())).unwrap();
        assert_eq!(result.iterations, 100);
        assert!(!result.glitch.detected);
    }

    #[test]
    fn glitch_detector_uses_relative_delta_size() {
        let detector = GlitchDetector::new(1e-3).unwrap();
        assert!(
            !detector
                .is_glitch(ComplexCoordinate::new(1.0, 0.0), ComplexCoordinate::new(0.01, 0.0))
        );
        assert!(
            detector.is_glitch(ComplexCoordinate::new(1.0, 0.0), ComplexCoordinate::new(0.1, 0.0))
        );
    }
}
