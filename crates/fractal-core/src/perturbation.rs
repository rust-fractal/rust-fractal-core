use crate::{
    ComplexCoordinate, FractalFormula, ReferenceOrbit, RenderError, SeriesApproximation, reference,
};
use fractal_numerics::{ComplexExtended, Extended};

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
        let reference_norm = reference.norm_squared();
        if reference_norm == 0.0 {
            return false;
        }
        delta.norm_squared() > self.tolerance * reference_norm
    }

    pub(crate) fn is_extended_glitch(
        &self,
        reference: ComplexExtended,
        delta: ComplexExtended,
    ) -> Result<bool, RenderError> {
        let reference_norm = reference.norm_squared().ok_or(RenderError::NumericOverflow)?;
        if reference_norm == Extended::ZERO {
            return Ok(false);
        }
        let delta_norm = delta.norm_squared().ok_or(RenderError::NumericOverflow)?;
        let tolerance = Extended::from_f64(self.tolerance).ok_or(RenderError::NumericOverflow)?;
        Ok(delta_norm
            > reference_norm.checked_mul(tolerance).ok_or(RenderError::NumericOverflow)?)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlitchResult {
    pub detected: bool,
    pub corrected: bool,
    pub rebase_count: u32,
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
    perturbate_with_series_cancel(
        reference_orbit,
        pixel,
        iteration_limit,
        glitch_detector,
        series,
        || false,
    )
}

pub(crate) fn perturbate_with_series_cancel<F>(
    reference_orbit: &ReferenceOrbit,
    pixel: ComplexCoordinate,
    iteration_limit: u32,
    glitch_detector: Option<GlitchDetector>,
    series: Option<&SeriesApproximation>,
    mut cancelled: F,
) -> Result<PerturbationResult, RenderError>
where
    F: FnMut() -> bool,
{
    if iteration_limit == 0 {
        return Err(RenderError::InvalidIterationLimit);
    }
    if !pixel.is_finite() {
        return Err(RenderError::NonFiniteCoordinate { real: pixel.real, imag: pixel.imag });
    }
    let delta_c = pixel.clone().sub(reference_orbit.center.clone());
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
        && let Some(skip) = series.largest_valid_iteration(iteration_limit, delta_c.clone(), 1e-3)
    {
        delta = series.evaluate(skip, delta_c.clone()).ok_or(RenderError::InvalidRawFrame)?;
        if let Some(reference_z) = reference_orbit.point(skip) {
            z = reference_z.add(delta.clone());
            iteration = skip;
        }
    }
    while iteration < iteration_limit {
        if iteration.is_multiple_of(128) && cancelled() {
            return Err(RenderError::Cancelled);
        }
        let reference_z = reference_orbit.point(iteration).unwrap_or(z.clone());
        let reference_next = reference_orbit.point(iteration + 1);
        if matches!(reference_orbit.formula, FractalFormula::BurningShip { .. }) {
            z = reference::step(z, pixel.clone(), reference_orbit.formula);
            delta = z.clone().sub(reference_next.clone().unwrap_or(z.clone()));
        } else if let Some(reference_next) = reference_next.clone() {
            delta = polynomial_delta(reference_z.clone(), delta, delta_c.clone(), power);
            z = reference_next.add(delta.clone());
        } else {
            z = reference::step(z, pixel.clone(), reference_orbit.formula);
        }

        let reference_for_glitch = reference_next.unwrap_or(reference_z.clone());
        if first_glitch.is_none()
            && glitch_detector
                .is_some_and(|detector| detector.is_glitch(reference_for_glitch, delta.clone()))
        {
            first_glitch = Some(iteration);
        }
        maximum_delta_norm = maximum_delta_norm.max(delta.clone().norm());
        iteration += 1;
        if z.clone().norm_squared() > 4.0 || !z.is_finite() {
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
    let escaped = z.clone().norm_squared() > 4.0;
    let z_norm = z.clone().norm();
    let smooth_iteration = if escaped && z_norm > 1.0 {
        iteration as f64 + 1.0 - z_norm.ln().ln() / (power as f64).ln()
    } else {
        iteration as f64
    };
    let derivative_norm = derivative.norm();
    let distance = if escaped && derivative_norm > 0.0 {
        0.5 * z_norm.ln() * z_norm / derivative_norm
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
            corrected: false,
            rebase_count: 0,
            maximum_delta_norm,
            first_iteration: first_glitch,
        },
    })
}

/// Deep perturbation kernel. The reference center and each pixel offset are
/// kept as separate compact extended values, so no absolute deep coordinate is
/// ever rounded through f64.
pub fn perturbate_extended(
    reference_orbit: &ReferenceOrbit,
    delta_c: ComplexExtended,
    iteration_limit: u32,
    glitch_detector: Option<GlitchDetector>,
) -> Result<PerturbationResult, RenderError> {
    perturbate_extended_with_cancel(
        reference_orbit,
        delta_c,
        iteration_limit,
        glitch_detector,
        None,
        || false,
    )
}

pub(crate) fn perturbate_extended_with_cancel<F>(
    reference_orbit: &ReferenceOrbit,
    delta_c: ComplexExtended,
    iteration_limit: u32,
    glitch_detector: Option<GlitchDetector>,
    series: Option<&SeriesApproximation>,
    mut cancelled: F,
) -> Result<PerturbationResult, RenderError>
where
    F: FnMut() -> bool,
{
    if iteration_limit == 0 {
        return Err(RenderError::InvalidIterationLimit);
    }
    if reference_orbit.statistics.generated_iterations < iteration_limit {
        return Err(RenderError::ReferenceGeneration(
            "deep perturbation requires a complete reference orbit".to_owned(),
        ));
    }
    let power = match reference_orbit.formula {
        FractalFormula::Mandelbrot { power } | FractalFormula::BurningShip { power } => power,
    };
    let mut delta = ComplexExtended::ZERO;
    let mut z = ComplexExtended::ZERO;
    let mut derivative = ComplexExtended::ZERO;
    let mut maximum_delta_norm: f64 = 0.0;
    let mut first_glitch = None;
    let mut rebased = false;
    let mut rebase_count: u32 = 0;
    let pixel_c =
        reference_orbit.center_compact.add(delta_c).ok_or(RenderError::NumericOverflow)?;
    let mut iteration = 0;
    if let Some(series) = series
        && let Some(skip) = series.largest_valid_iteration_extended(iteration_limit, delta_c, 1e-3)
    {
        delta = series.evaluate_extended(skip, delta_c).ok_or(RenderError::NumericOverflow)?;
        z = reference_orbit
            .compact_point(skip)
            .ok_or(RenderError::ReferenceGeneration("reference point is missing".to_owned()))?
            .add(delta)
            .ok_or(RenderError::NumericOverflow)?;
        iteration = skip;
    }

    while iteration < iteration_limit {
        if iteration.is_multiple_of(128) && cancelled() {
            return Err(RenderError::Cancelled);
        }
        let reference_z = reference_orbit
            .compact_point(iteration)
            .ok_or(RenderError::ReferenceGeneration("reference point is missing".to_owned()))?;
        let reference_next = reference_orbit
            .compact_point(iteration + 1)
            .ok_or(RenderError::ReferenceGeneration("reference point is missing".to_owned()))?;
        match reference_orbit.formula {
            FractalFormula::Mandelbrot { power: 2 } => {
                if rebased {
                    z = z
                        .multiply(z)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(pixel_c)
                        .ok_or(RenderError::NumericOverflow)?;
                    delta = z.sub(reference_next).ok_or(RenderError::NumericOverflow)?;
                } else {
                    let two_reference =
                        reference_z.scale(2.0).ok_or(RenderError::NumericOverflow)?;
                    delta = delta
                        .multiply(two_reference)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(delta_c)
                        .ok_or(RenderError::NumericOverflow)?;
                    z = reference_next.add(delta).ok_or(RenderError::NumericOverflow)?;
                    derivative = derivative
                        .multiply(two_reference)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(
                            ComplexExtended::from_f64(1.0, 0.0)
                                .ok_or(RenderError::NumericOverflow)?,
                        )
                        .ok_or(RenderError::NumericOverflow)?;
                }
            }
            FractalFormula::Mandelbrot { power: 3 } => {
                if rebased {
                    z = complex_power_extended(z, 3)?
                        .add(pixel_c)
                        .ok_or(RenderError::NumericOverflow)?;
                    delta = z.sub(reference_next).ok_or(RenderError::NumericOverflow)?;
                } else {
                    let delta_squared =
                        delta.multiply(delta).ok_or(RenderError::NumericOverflow)?;
                    let reference_squared =
                        reference_z.multiply(reference_z).ok_or(RenderError::NumericOverflow)?;
                    let term = reference_squared
                        .multiply(delta)
                        .ok_or(RenderError::NumericOverflow)?
                        .scale(3.0)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(
                            reference_z
                                .multiply(delta_squared)
                                .ok_or(RenderError::NumericOverflow)?
                                .scale(3.0)
                                .ok_or(RenderError::NumericOverflow)?,
                        )
                        .ok_or(RenderError::NumericOverflow)?
                        .add(delta_squared.multiply(delta).ok_or(RenderError::NumericOverflow)?)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(delta_c)
                        .ok_or(RenderError::NumericOverflow)?;
                    delta = term;
                    z = reference_next.add(delta).ok_or(RenderError::NumericOverflow)?;
                    let derivative_factor =
                        reference_squared.scale(3.0).ok_or(RenderError::NumericOverflow)?;
                    derivative = derivative
                        .multiply(derivative_factor)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(
                            ComplexExtended::from_f64(1.0, 0.0)
                                .ok_or(RenderError::NumericOverflow)?,
                        )
                        .ok_or(RenderError::NumericOverflow)?;
                }
            }
            FractalFormula::Mandelbrot { power } => {
                if rebased {
                    z = complex_power_extended(z, power)?
                        .add(pixel_c)
                        .ok_or(RenderError::NumericOverflow)?;
                    delta = z.sub(reference_next).ok_or(RenderError::NumericOverflow)?;
                } else {
                    delta = polynomial_delta_extended(reference_z, delta, delta_c, power)?;
                    z = reference_next.add(delta).ok_or(RenderError::NumericOverflow)?;
                    let factor = complex_power_extended(reference_z, power - 1)?
                        .scale(f64::from(power))
                        .ok_or(RenderError::NumericOverflow)?;
                    derivative = derivative
                        .multiply(factor)
                        .ok_or(RenderError::NumericOverflow)?
                        .add(
                            ComplexExtended::from_f64(1.0, 0.0)
                                .ok_or(RenderError::NumericOverflow)?,
                        )
                        .ok_or(RenderError::NumericOverflow)?;
                }
            }
            FractalFormula::BurningShip { power } => {
                let absolute = z.abs_components();
                z = complex_power_extended(absolute, power)?
                    .add(pixel_c)
                    .ok_or(RenderError::NumericOverflow)?;
                delta = z.sub(reference_next).ok_or(RenderError::NumericOverflow)?;
            }
        }

        let glitch = if let Some(detector) = glitch_detector {
            detector.is_extended_glitch(reference_next, delta)?
        } else {
            false
        };
        if first_glitch.is_none() && glitch {
            first_glitch = Some(iteration);
            rebased = true;
            rebase_count = rebase_count.saturating_add(1);
        }
        maximum_delta_norm = maximum_delta_norm
            .max(delta.norm_squared().ok_or(RenderError::NumericOverflow)?.to_f64().sqrt());
        iteration += 1;
        if z.norm_squared().ok_or(RenderError::NumericOverflow)? > Extended::FOUR {
            break;
        }
    }

    let escaped = z.norm_squared().ok_or(RenderError::NumericOverflow)? > Extended::FOUR;
    let z_norm = z.norm_squared().ok_or(RenderError::NumericOverflow)?.to_f64().sqrt();
    let smooth_iteration = if escaped && z_norm > 1.0 {
        f64::from(iteration) + 1.0 - z_norm.ln().ln() / f64::from(power).ln()
    } else {
        f64::from(iteration)
    };
    let derivative_norm =
        derivative.norm_squared().map(Extended::to_f64).unwrap_or(f64::INFINITY).sqrt();
    let distance = if escaped && derivative_norm > 0.0 && derivative_norm.is_finite() {
        0.5 * z_norm.ln() * z_norm / derivative_norm
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
            corrected: rebase_count != 0,
            rebase_count,
            maximum_delta_norm,
            first_iteration: first_glitch,
        },
    })
}

fn polynomial_delta_extended(
    reference_z: ComplexExtended,
    delta: ComplexExtended,
    delta_c: ComplexExtended,
    power: u32,
) -> Result<ComplexExtended, RenderError> {
    let mut result = delta_c;
    let mut delta_power =
        ComplexExtended::from_f64(1.0, 0.0).ok_or(RenderError::NumericOverflow)?;
    for order in 1..=power {
        delta_power = delta_power.multiply(delta).ok_or(RenderError::NumericOverflow)?;
        let reference_power = complex_power_extended(reference_z, power - order)?;
        result = result
            .add(
                reference_power
                    .multiply(delta_power)
                    .ok_or(RenderError::NumericOverflow)?
                    .scale(binomial(power, order) as f64)
                    .ok_or(RenderError::NumericOverflow)?,
            )
            .ok_or(RenderError::NumericOverflow)?;
    }
    Ok(result)
}

fn complex_power_extended(
    value: ComplexExtended,
    power: u32,
) -> Result<ComplexExtended, RenderError> {
    let mut result = ComplexExtended::from_f64(1.0, 0.0).ok_or(RenderError::NumericOverflow)?;
    let mut base = value;
    let mut power = power;
    while power != 0 {
        if power & 1 == 1 {
            result = result.multiply(base).ok_or(RenderError::NumericOverflow)?;
        }
        power >>= 1;
        if power != 0 {
            base = base.multiply(base).ok_or(RenderError::NumericOverflow)?;
        }
    }
    Ok(result)
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
        delta_power = delta_power.multiply(delta.clone());
        let coefficient = binomial(power, order) as f64;
        let reference_term =
            reference::complex_power(reference_z.real, reference_z.imag, power - order);
        result = result.add(
            ComplexCoordinate::new(reference_term.0, reference_term.1)
                .multiply(delta_power.clone())
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
        let reference = generate_reference(center.clone(), FractalFormula::default(), 100).unwrap();
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

    #[test]
    fn deep_glitch_is_rebased_in_extended_arithmetic() {
        let reference =
            generate_reference(ComplexCoordinate::new(-0.75, 0.0), FractalFormula::default(), 32)
                .unwrap();
        let result = perturbate_extended(
            &reference,
            ComplexExtended::from_f64(0.01, 0.0).unwrap(),
            32,
            Some(GlitchDetector::new(1e-300).unwrap()),
        )
        .unwrap();
        assert!(result.glitch.detected);
        assert!(result.glitch.corrected);
        assert!(result.glitch.rebase_count >= 1);
    }
}
