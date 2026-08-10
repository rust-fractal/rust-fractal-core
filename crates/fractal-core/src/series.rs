use crate::{ComplexCoordinate, FractalFormula, ReferenceOrbit, RenderError};
use fractal_numerics::{ComplexExtended, Extended};

#[derive(Clone, Debug, PartialEq)]
pub struct SeriesCoefficients {
    pub iteration: u32,
    pub order: u32,
    pub terms: Vec<ComplexCoordinate>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SeriesApproximation {
    pub order: u32,
    pub storage_interval: u32,
    pub states: Vec<SeriesCoefficients>,
}

pub fn generate_series(
    reference: &ReferenceOrbit,
    order: u32,
) -> Result<SeriesApproximation, RenderError> {
    if order == 0 || order > 64 {
        return Err(RenderError::InvalidApproximationOrder { order });
    }
    if !matches!(reference.formula, FractalFormula::Mandelbrot { power: 2 }) {
        return Err(RenderError::InvalidApproximationOrder { order });
    }

    let mut terms = vec![ComplexCoordinate::new(0.0, 0.0); order as usize + 1];
    let mut next = terms.clone();
    let storage_interval = if reference.statistics.generated_iterations > 256 { 32 } else { 1 };
    let mut states = Vec::with_capacity(
        reference.statistics.generated_iterations.div_ceil(storage_interval) as usize + 1,
    );
    states.push(SeriesCoefficients { iteration: 0, order, terms: terms.clone() });
    for iteration in 0..reference.statistics.generated_iterations {
        let z = reference.point(iteration).ok_or(RenderError::InvalidRawFrame)?;
        next.fill(ComplexCoordinate::new(0.0, 0.0));
        for index in 1..terms.len() {
            next[index] = terms[index].clone().scale(2.0).multiply(z.clone());
        }
        next[1] = next[1].clone().add(ComplexCoordinate::new(1.0, 0.0));
        for left in 1..terms.len() {
            for right in 1..terms.len() - left {
                let index = left + right;
                next[index] =
                    next[index].clone().add(terms[left].clone().multiply(terms[right].clone()));
            }
        }
        core::mem::swap(&mut terms, &mut next);
        let next_iteration = iteration + 1;
        if next_iteration.is_multiple_of(storage_interval)
            || next_iteration == reference.statistics.generated_iterations
        {
            states.push(SeriesCoefficients {
                iteration: next_iteration,
                order,
                terms: terms.clone(),
            });
        }
    }
    Ok(SeriesApproximation { order, storage_interval, states })
}

impl SeriesApproximation {
    pub fn evaluate(
        &self,
        iteration: u32,
        delta_c: ComplexCoordinate,
    ) -> Option<ComplexCoordinate> {
        let state = self.states.iter().find(|state| state.iteration == iteration)?;
        let mut result = ComplexCoordinate::new(0.0, 0.0);
        let mut power = ComplexCoordinate::new(1.0, 0.0);
        for term in state.terms.iter().skip(1) {
            power = power.multiply(delta_c.clone());
            result = result.add(term.clone().multiply(power.clone()));
        }
        Some(result)
    }

    pub fn is_valid(&self, iteration: u32, delta_c: ComplexCoordinate, threshold: f64) -> bool {
        self.evaluate(iteration, delta_c)
            .is_some_and(|delta| delta.norm_squared() <= threshold * threshold)
    }

    pub fn evaluate_extended(
        &self,
        iteration: u32,
        delta_c: ComplexExtended,
    ) -> Option<ComplexExtended> {
        let state = self.states.iter().find(|state| state.iteration == iteration)?;
        let mut result = ComplexExtended::ZERO;
        let mut power = ComplexExtended::from_f64(1.0, 0.0)?;
        for term in state.terms.iter().skip(1) {
            power = power.multiply(delta_c)?;
            let coefficient = ComplexExtended::from_f64(term.real, term.imag)?;
            result = result.add(coefficient.multiply(power)?)?;
        }
        Some(result)
    }

    pub fn largest_valid_iteration_extended(
        &self,
        maximum_iteration: u32,
        delta_c: ComplexExtended,
        threshold: f64,
    ) -> Option<u32> {
        let threshold = Extended::from_f64(threshold)?;
        let threshold_squared = threshold.checked_mul(threshold)?;
        self.states
            .iter()
            .rev()
            .filter(|state| state.iteration <= maximum_iteration && state.iteration != 0)
            .find(|state| {
                self.evaluate_extended(state.iteration, delta_c)
                    .and_then(|delta| delta.norm_squared())
                    .is_some_and(|norm| norm <= threshold_squared)
            })
            .map(|state| state.iteration)
    }

    pub fn largest_valid_iteration(
        &self,
        maximum_iteration: u32,
        delta_c: ComplexCoordinate,
        threshold: f64,
    ) -> Option<u32> {
        self.states
            .iter()
            .rev()
            .filter(|state| state.iteration <= maximum_iteration && state.iteration != 0)
            .find(|state| self.is_valid(state.iteration, delta_c.clone(), threshold))
            .map(|state| state.iteration)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate_reference;

    #[test]
    fn quadratic_series_matches_direct_delta_for_small_offsets() {
        let reference =
            generate_reference(ComplexCoordinate::new(-0.75, 0.0), FractalFormula::default(), 4)
                .unwrap();
        let series = generate_series(&reference, 4).unwrap();
        let delta_c = ComplexCoordinate::new(1e-8, -2e-8);
        let estimate = series.evaluate(2, delta_c.clone()).unwrap();
        assert!(estimate.norm() < 1e-5);
        assert!(series.is_valid(2, delta_c, 1e-4));
    }
}
