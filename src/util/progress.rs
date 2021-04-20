use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};

pub struct ProgressCounters {
    pub reference: Arc<AtomicUsize>,
    pub reference_maximum: Arc<AtomicUsize>,
    pub reference_count: Arc<AtomicUsize>,
    pub series_approximation: Arc<AtomicUsize>,
    pub min_series_approximation: Arc<AtomicUsize>,
    pub max_series_approximation: Arc<AtomicUsize>,
    pub series_validation: Arc<AtomicUsize>,
    pub iteration: Arc<AtomicUsize>,
    pub glitched_maximum: Arc<AtomicUsize>,
}

impl ProgressCounters {
    pub fn new(maximum_iteration: usize) -> ProgressCounters {
        ProgressCounters {
            reference: Arc::new(AtomicUsize::new(1)),
            reference_maximum: Arc::new(AtomicUsize::new(maximum_iteration - 1)),
            reference_count: Arc::new(AtomicUsize::new(1)),
            series_approximation: Arc::new(AtomicUsize::new(0)),
            min_series_approximation: Arc::new(AtomicUsize::new(1)),
            max_series_approximation: Arc::new(AtomicUsize::new(1)),
            series_validation: Arc::new(AtomicUsize::new(0)),
            iteration: Arc::new(AtomicUsize::new(0)),
            glitched_maximum: Arc::new(AtomicUsize::new(0))
        }
    }

    // TODO just set these to zero rather than reset
    // Reset without the series approximation changed
    pub fn reset(&mut self) {
        self.min_series_approximation.store(1, Ordering::SeqCst);
        self.max_series_approximation.store(1, Ordering::SeqCst);
        self.series_validation.store(0, Ordering::SeqCst);
        self.iteration.store(0, Ordering::SeqCst);
        self.glitched_maximum.store(0, Ordering::SeqCst);
        self.reference_count.store(1, Ordering::SeqCst);
    }

    // TODO just set these to zero rather than reset
    // Reset with the series approximation changed
    pub fn reset_series_approximation(&mut self) {
        self.series_approximation.store(0, Ordering::SeqCst);
    }

    pub fn reset_all(&mut self, maximum_iteration: usize) {
        self.min_series_approximation.store(1, Ordering::SeqCst);
        self.max_series_approximation.store(1, Ordering::SeqCst);
        self.series_validation.store(0, Ordering::SeqCst);
        self.iteration.store(0, Ordering::SeqCst);
        self.glitched_maximum.store(0, Ordering::SeqCst);
        self.series_approximation.store(0, Ordering::SeqCst);
        self.reference.store(1, Ordering::SeqCst);
        self.reference_maximum.store(maximum_iteration - 1, Ordering::SeqCst);
        self.reference_count.store(1, Ordering::SeqCst);
    }
}