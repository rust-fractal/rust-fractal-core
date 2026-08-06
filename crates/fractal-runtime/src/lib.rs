//! Native scheduling around the deterministic [`fractal_core`] computation.
//!
//! The core computes one tile at a time and owns no threads. This crate owns
//! tile partitioning, cancellation, and progressive event delivery.

use fractal_core::{
    RenderError, RenderEvent, RenderProgress, RenderRequest, RenderResult, RenderSession,
    RenderStarted, RenderTile, Tile,
};
use std::sync::{
    Arc,
    atomic::{AtomicBool, AtomicU64, Ordering},
    mpsc::{Receiver, sync_channel},
};
use std::thread;

#[derive(Clone, Default)]
pub struct CancellationToken {
    generation: u64,
    cancelled: Arc<AtomicBool>,
}

impl CancellationToken {
    pub fn for_generation(generation: u64) -> Self {
        Self { generation, cancelled: Arc::new(AtomicBool::new(false)) }
    }

    pub const fn generation(&self) -> u64 {
        self.generation
    }

    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Release);
    }

    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Acquire)
    }
}

#[derive(Debug)]
pub struct RenderExecutor {
    tile_width: u32,
    tile_height: u32,
    progress_every_tiles: usize,
    #[cfg(feature = "native")]
    pool: Option<Arc<rayon::ThreadPool>>,
}

impl RenderExecutor {
    pub const fn new(tile_width: u32, tile_height: u32) -> Self {
        Self {
            tile_width,
            tile_height,
            progress_every_tiles: 1,
            #[cfg(feature = "native")]
            pool: None,
        }
    }

    pub const fn tile_size(&self) -> (u32, u32) {
        (self.tile_width, self.tile_height)
    }

    pub fn with_workers(
        tile_width: u32,
        tile_height: u32,
        workers: usize,
    ) -> Result<Self, RuntimeError> {
        if workers == 0 {
            return Err(RuntimeError::InvalidWorkerCount);
        }
        let mut executor = Self::new(tile_width, tile_height);
        #[cfg(feature = "native")]
        {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .map_err(|error| RuntimeError::ThreadPool(error.to_string()))?;
            executor.pool = Some(Arc::new(pool));
        }
        Ok(executor)
    }

    pub const fn progress_every_tiles(mut self, tiles: usize) -> Self {
        self.progress_every_tiles = if tiles == 0 { 1 } else { tiles };
        self
    }

    pub fn render(
        &self,
        request: &RenderRequest,
        cancellation: &CancellationToken,
    ) -> Vec<RenderEvent> {
        if let Err(error) = request.validate() {
            return vec![RenderEvent::Failed(error)];
        }
        if self.tile_width == 0 || self.tile_height == 0 {
            return vec![RenderEvent::Failed(RenderError::InvalidTile)];
        }
        if cancellation.is_cancelled() || cancellation.generation() != request.generation {
            return vec![RenderEvent::Cancelled { generation: request.generation }];
        }

        let session = match RenderSession::new(*request) {
            Ok(session) => session,
            Err(error) => return vec![RenderEvent::Failed(error)],
        };
        let tiles = tiles_for(request, self.tile_width, self.tile_height);
        let mut events = vec![RenderEvent::Started(RenderStarted {
            generation: request.generation,
            width: request.viewport.width,
            height: request.viewport.height,
        })];

        let tile_results = self.compute_tiles(&session, &tiles, cancellation);
        let mut result = empty_result(request);
        let total_pixels = u64::from(request.viewport.width) * u64::from(request.viewport.height);
        let mut completed_pixels = 0;

        for tile_result in tile_results {
            if cancellation.is_cancelled() {
                events.push(RenderEvent::Cancelled { generation: request.generation });
                return events;
            }

            let tile = match tile_result {
                Ok(tile) => tile,
                Err(RenderError::Cancelled) => {
                    events.push(RenderEvent::Cancelled { generation: request.generation });
                    return events;
                }
                Err(error) => {
                    events.push(RenderEvent::Failed(error));
                    return events;
                }
            };
            merge_tile(&mut result, &tile);
            completed_pixels += u64::from(tile.tile.width) * u64::from(tile.tile.height);
            events.push(RenderEvent::TileReady(tile));
            let tile_number =
                events.iter().filter(|event| matches!(event, RenderEvent::TileReady(_)));
            if tile_number.count().is_multiple_of(self.progress_every_tiles)
                || completed_pixels == total_pixels
            {
                events.push(RenderEvent::Progress(RenderProgress {
                    generation: request.generation,
                    completed_pixels,
                    total_pixels,
                }));
            }
        }

        events.push(RenderEvent::Finished(result));
        events
    }

    pub fn render_bounded(
        &self,
        request: RenderRequest,
        cancellation: CancellationToken,
        capacity: usize,
    ) -> Result<RenderEventStream, RuntimeError> {
        if capacity == 0 {
            return Err(RuntimeError::InvalidChannelCapacity);
        }
        let (sender, receiver) = sync_channel(capacity);
        let executor = self.clone();
        thread::spawn(move || {
            for event in executor.render(&request, &cancellation) {
                if sender.send(event).is_err() {
                    break;
                }
            }
        });
        Ok(RenderEventStream { receiver })
    }

    #[cfg(feature = "native")]
    fn compute_tiles(
        &self,
        session: &RenderSession,
        tiles: &[Tile],
        cancellation: &CancellationToken,
    ) -> Vec<Result<RenderTile, RenderError>> {
        use rayon::prelude::*;
        let compute = || {
            tiles
                .par_iter()
                .map(|&tile| session.compute_tile_with_cancel(tile, || cancellation.is_cancelled()))
                .collect()
        };
        self.pool.as_ref().map_or_else(compute, |pool| pool.install(compute))
    }

    #[cfg(not(feature = "native"))]
    fn compute_tiles(
        &self,
        session: &RenderSession,
        tiles: &[Tile],
        cancellation: &CancellationToken,
    ) -> Vec<Result<RenderTile, RenderError>> {
        tiles
            .iter()
            .copied()
            .map(|tile| session.compute_tile_with_cancel(tile, || cancellation.is_cancelled()))
            .collect()
    }
}

impl Default for RenderExecutor {
    fn default() -> Self {
        Self::new(64, 64)
    }
}

impl Clone for RenderExecutor {
    fn clone(&self) -> Self {
        Self {
            tile_width: self.tile_width,
            tile_height: self.tile_height,
            progress_every_tiles: self.progress_every_tiles,
            #[cfg(feature = "native")]
            pool: self.pool.clone(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RuntimeError {
    InvalidWorkerCount,
    InvalidChannelCapacity,
    ThreadPool(String),
}

impl std::fmt::Display for RuntimeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidWorkerCount => f.write_str("worker count must be non-zero"),
            Self::InvalidChannelCapacity => f.write_str("channel capacity must be non-zero"),
            Self::ThreadPool(error) => write!(f, "thread pool could not be created: {error}"),
        }
    }
}

impl std::error::Error for RuntimeError {}

pub struct RenderEventStream {
    receiver: Receiver<RenderEvent>,
}

impl Iterator for RenderEventStream {
    type Item = RenderEvent;

    fn next(&mut self) -> Option<Self::Item> {
        self.receiver.recv().ok()
    }
}

pub struct GenerationGate {
    current: AtomicU64,
}

impl Default for GenerationGate {
    fn default() -> Self {
        Self { current: AtomicU64::new(0) }
    }
}

impl GenerationGate {
    pub fn begin(&self, generation: u64) {
        self.current.fetch_max(generation, Ordering::AcqRel);
    }

    pub fn accepts(&self, generation: u64) -> bool {
        self.current.load(Ordering::Acquire) == generation
    }
}

fn tiles_for(request: &RenderRequest, tile_width: u32, tile_height: u32) -> Vec<Tile> {
    let mut tiles = Vec::new();
    let mut y = 0;
    while y < request.viewport.height {
        let mut x = 0;
        while x < request.viewport.width {
            tiles.push(Tile {
                x,
                y,
                width: tile_width.min(request.viewport.width - x),
                height: tile_height.min(request.viewport.height - y),
            });
            x += tile_width;
        }
        y += tile_height;
    }
    tiles
}

fn empty_result(request: &RenderRequest) -> RenderResult {
    let len = (request.viewport.width * request.viewport.height) as usize;
    RenderResult {
        generation: request.generation,
        width: request.viewport.width,
        height: request.viewport.height,
        iteration_limit: request.settings.iteration_limit,
        iterations: vec![0; len],
        smooth_iteration: request.settings.output.smooth_iteration.then(|| vec![0.0; len]),
        distance: request.settings.output.distance.then(|| vec![0.0; len]),
        stripes: request.settings.output.stripes.then(|| vec![0.0; len]),
        validity: request.settings.output.validity.then(|| vec![0; len]),
    }
}

fn merge_tile(result: &mut RenderResult, tile: &RenderTile) {
    for local_y in 0..tile.tile.height {
        let destination_start = ((tile.tile.y + local_y) * result.width + tile.tile.x) as usize;
        let source_start = (local_y * tile.tile.width) as usize;
        let len = tile.tile.width as usize;
        result.iterations[destination_start..destination_start + len]
            .copy_from_slice(&tile.iterations[source_start..source_start + len]);
        if let (Some(destination), Some(source)) =
            (result.smooth_iteration.as_mut(), tile.smooth_iteration.as_ref())
        {
            destination[destination_start..destination_start + len]
                .copy_from_slice(&source[source_start..source_start + len]);
        }
        if let (Some(destination), Some(source)) =
            (result.distance.as_mut(), tile.distance.as_ref())
        {
            destination[destination_start..destination_start + len]
                .copy_from_slice(&source[source_start..source_start + len]);
        }
        if let (Some(destination), Some(source)) = (result.stripes.as_mut(), tile.stripes.as_ref())
        {
            destination[destination_start..destination_start + len]
                .copy_from_slice(&source[source_start..source_start + len]);
        }
        if let (Some(destination), Some(source)) =
            (result.validity.as_mut(), tile.validity.as_ref())
        {
            destination[destination_start..destination_start + len]
                .copy_from_slice(&source[source_start..source_start + len]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fractal_core::{
        ComplexCoordinate, FractalFormula, OutputChannels, RenderSettings, Viewport, Zoom,
    };

    fn request() -> RenderRequest {
        RenderRequest {
            generation: 7,
            viewport: Viewport {
                center: ComplexCoordinate::new(-0.75, 0.0),
                zoom: Zoom::unit(),
                rotation_radians: 0.0,
                width: 9,
                height: 7,
            },
            settings: RenderSettings {
                fractal: FractalFormula::default(),
                output: OutputChannels {
                    smooth_iteration: true,
                    distance: true,
                    ..OutputChannels::default()
                },
                ..RenderSettings::default()
            },
        }
    }

    #[test]
    fn events_are_in_deterministic_tile_order() {
        let executor = RenderExecutor::new(4, 3);
        let token = CancellationToken::for_generation(7);
        let first = executor.render(&request(), &token);
        let second = executor.render(&request(), &token);
        assert_eq!(first, second);
        assert!(matches!(first.first(), Some(RenderEvent::Started(_))));
        assert!(matches!(first.last(), Some(RenderEvent::Finished(_))));
    }

    #[test]
    fn worker_counts_produce_equal_results() {
        let one = RenderExecutor::with_workers(4, 4, 1).unwrap();
        let many = RenderExecutor::with_workers(4, 4, 4).unwrap();
        let one_events = one.render(&request(), &CancellationToken::for_generation(7));
        let many_events = many.render(&request(), &CancellationToken::for_generation(7));
        assert_eq!(one_events, many_events);
    }

    #[test]
    fn cancellation_is_reported_before_work_starts() {
        let executor = RenderExecutor::default();
        let token = CancellationToken::for_generation(7);
        token.cancel();
        assert_eq!(
            executor.render(&request(), &token),
            vec![RenderEvent::Cancelled { generation: 7 }]
        );
    }

    #[test]
    fn stale_generations_are_rejected() {
        let gate = GenerationGate::default();
        gate.begin(2);
        assert!(!gate.accepts(1));
        assert!(gate.accepts(2));
    }

    #[test]
    fn bounded_stream_delivers_a_complete_render() {
        let executor = RenderExecutor::new(4, 4);
        let stream =
            executor.render_bounded(request(), CancellationToken::for_generation(7), 2).unwrap();
        let events = stream.collect::<Vec<_>>();
        assert!(matches!(events.first(), Some(RenderEvent::Started(_))));
        assert!(matches!(events.last(), Some(RenderEvent::Finished(_))));
    }
}
