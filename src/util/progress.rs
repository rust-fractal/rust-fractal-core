use atomic_counter::RelaxedCounter;

use std::sync::Arc;

pub struct ProgressCounters {
    pub reference: Arc<RelaxedCounter>,
    pub reference_maximum: Arc<RelaxedCounter>,
    pub series_approximation: Arc<RelaxedCounter>,
    pub series_validation: Arc<RelaxedCounter>,
    pub iteration: Arc<RelaxedCounter>
}

impl ProgressCounters {
    pub fn new(maximum_iteration: usize) -> ProgressCounters {
        ProgressCounters {
            reference: Arc::new(RelaxedCounter::new(0)),
            reference_maximum: Arc::new(RelaxedCounter::new(maximum_iteration)),
            series_approximation: Arc::new(RelaxedCounter::new(0)),
            series_validation: Arc::new(RelaxedCounter::new(0)),
            iteration: Arc::new(RelaxedCounter::new(0)),
        }
    }

    pub fn reset(&mut self, maximum_iteration: usize) {
        self.reference = Arc::new(RelaxedCounter::new(0));
        self.reference_maximum = Arc::new(RelaxedCounter::new(maximum_iteration));
        self.series_approximation = Arc::new(RelaxedCounter::new(0));
        self.series_validation = Arc::new(RelaxedCounter::new(0));
        self.iteration = Arc::new(RelaxedCounter::new(0));
    }
}