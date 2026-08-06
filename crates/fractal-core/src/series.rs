use crate::{ComplexCoordinate, FractalFormula, ReferenceOrbit, RenderError};

#[derive(Clone, Debug, PartialEq)]
pub struct SeriesCoefficients {
    pub iteration: u32,
    pub order: u32,
    pub terms: Vec<ComplexCoordinate>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SeriesApproximation {
    pub order: u32,
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
    let mut states = Vec::with_capacity(reference.points.len());
    states.push(SeriesCoefficients { iteration: 0, order, terms: terms.clone() });
    for iteration in 0..reference.statistics.generated_iterations {
        let z = reference.point(iteration).ok_or(RenderError::InvalidRawFrame)?;
        next.fill(ComplexCoordinate::new(0.0, 0.0));
        for index in 1..terms.len() {
            next[index] = terms[index].scale(2.0).multiply(z);
        }
        next[1] = next[1].add(ComplexCoordinate::new(1.0, 0.0));
        for left in 1..terms.len() {
            for right in 1..terms.len() - left {
                let index = left + right;
                next[index] = next[index].add(terms[left].multiply(terms[right]));
            }
        }
        core::mem::swap(&mut terms, &mut next);
        states.push(SeriesCoefficients { iteration: iteration + 1, order, terms: terms.clone() });
    }
    Ok(SeriesApproximation { order, states })
}

impl SeriesApproximation {
    pub fn evaluate(
        &self,
        iteration: u32,
        delta_c: ComplexCoordinate,
    ) -> Option<ComplexCoordinate> {
        let state = self.states.get(iteration as usize)?;
        let mut result = ComplexCoordinate::new(0.0, 0.0);
        let mut power = ComplexCoordinate::new(1.0, 0.0);
        for term in state.terms.iter().skip(1) {
            power = power.multiply(delta_c);
            result = result.add(term.multiply(power));
        }
        Some(result)
    }

    pub fn is_valid(&self, iteration: u32, delta_c: ComplexCoordinate, threshold: f64) -> bool {
        self.evaluate(iteration, delta_c)
            .is_some_and(|delta| delta.norm_squared() <= threshold * threshold)
    }

    pub fn largest_valid_iteration(
        &self,
        maximum_iteration: u32,
        delta_c: ComplexCoordinate,
        threshold: f64,
    ) -> Option<u32> {
        let maximum_iteration = maximum_iteration.min(self.states.len().saturating_sub(1) as u32);
        (1..=maximum_iteration)
            .rev()
            .find(|&iteration| self.is_valid(iteration, delta_c, threshold))
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
        let estimate = series.evaluate(2, delta_c).unwrap();
        assert!(estimate.norm() < 1e-5);
        assert!(series.is_valid(2, delta_c, 1e-4));
    }
}
