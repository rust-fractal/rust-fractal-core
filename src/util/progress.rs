use atomic_counter::RelaxedCounter;

use std::sync::Arc;

pub struct ProgressCounters {
    pub reference: Arc<RelaxedCounter>,
    pub reference_maximum: Arc<RelaxedCounter>,
    pub series_approximation: Arc<RelaxedCounter>,
    pub series_validation: Arc<RelaxedCounter>,
    pub iteration: Arc<RelaxedCounter>,
    pub glitched_maximum: Arc<RelaxedCounter>,
}

impl ProgressCounters {
    pub fn new(maximum_iteration: usize) -> ProgressCounters {
        ProgressCounters {
            reference: Arc::new(RelaxedCounter::new(0)),
            reference_maximum: Arc::new(RelaxedCounter::new(maximum_iteration - 1)),
            series_approximation: Arc::new(RelaxedCounter::new(0)),
            series_validation: Arc::new(RelaxedCounter::new(0)),
            iteration: Arc::new(RelaxedCounter::new(0)),
            glitched_maximum: Arc::new(RelaxedCounter::new(0)),
        }
    }

    // TODO just set these to zero rather than reset
    // Reset without the series approximation changed
    pub fn reset(&mut self) {
        self.series_validation = Arc::new(RelaxedCounter::new(0));
        self.iteration = Arc::new(RelaxedCounter::new(0));
        self.glitched_maximum = Arc::new(RelaxedCounter::new(0));
    }

    // TODO just set these to zero rather than reset
    // Reset with the series approximation changed
    pub fn reset_series_approximation(&mut self) {
        self.series_approximation = Arc::new(RelaxedCounter::new(0));
    }
}