//! Native scheduling around the deterministic [`fractal_core`] computation.
//!
//! The core computes one tile at a time and owns no threads. This crate owns
//! tile partitioning, cancellation, and progressive event delivery.

use fractal_core::{
    RawFrameView, RenderError, RenderEvent, RenderProgress, RenderRequest, RenderResult,
    RenderSession, RenderStarted, RenderTile, Tile,
};
use std::collections::BTreeMap;
use std::sync::{
    Arc,
    atomic::{AtomicBool, AtomicU64, Ordering},
    mpsc::{Receiver, SyncSender, sync_channel},
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
        self.cancelled.store(true, Ordering::Relaxed);
    }

    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
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
        let capacity = tiles_for(request, self.tile_width.max(1), self.tile_height.max(1))
            .len()
            .saturating_mul(3)
            .saturating_add(4)
            .max(1);
        let (sender, receiver) = sync_channel(capacity);
        self.render_to_sender(request.clone(), cancellation.clone(), capacity, true, &sender);
        drop(sender);
        receiver.into_iter().collect()
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
            executor.render_to_sender(request, cancellation, capacity, false, &sender);
        });
        Ok(RenderEventStream { receiver })
    }

    /// Computes a complete frame into one preallocated set of channels. Row
    /// views are split safely and filled in parallel, avoiding owned tiles and
    /// the progressive merge copy.
    #[cfg(feature = "native")]
    pub fn compute_frame_full(
        &self,
        request: &RenderRequest,
        cancellation: &CancellationToken,
    ) -> Result<RenderResult, RenderError> {
        request.validate()?;
        if self.tile_width == 0 || self.tile_height == 0 {
            return Err(RenderError::InvalidTile);
        }
        if cancellation.is_cancelled() || cancellation.generation() != request.generation {
            return Err(RenderError::Cancelled);
        }
        let session = Arc::new(RenderSession::new(request.clone())?);
        let mut result = empty_result(request);
        let view = result.view_mut();
        let compute = || compute_frame_view(&session, view, 0, cancellation);
        if let Some(pool) = &self.pool {
            pool.install(compute)?;
        } else {
            compute()?;
        }
        Ok(result)
    }

    fn render_to_sender(
        &self,
        request: RenderRequest,
        cancellation: CancellationToken,
        capacity: usize,
        ordered: bool,
        sender: &SyncSender<RenderEvent>,
    ) {
        if let Err(error) = request.validate() {
            let _ = sender.send(RenderEvent::Failed(error));
            return;
        }
        if self.tile_width == 0 || self.tile_height == 0 {
            let _ = sender.send(RenderEvent::Failed(RenderError::InvalidTile));
            return;
        }
        if cancellation.is_cancelled() || cancellation.generation() != request.generation {
            let _ = sender.send(RenderEvent::Cancelled { generation: request.generation });
            return;
        }
        let session = match RenderSession::new(request.clone()) {
            Ok(session) => Arc::new(session),
            Err(error) => {
                let _ = sender.send(RenderEvent::Failed(error));
                return;
            }
        };
        let tiles = tiles_for(&request, self.tile_width, self.tile_height);
        let _ = sender.send(RenderEvent::Started(RenderStarted {
            generation: request.generation,
            width: request.viewport.width,
            height: request.viewport.height,
        }));
        let total_pixels = u64::from(request.viewport.width) * u64::from(request.viewport.height);
        let (tile_sender, tile_receiver) =
            sync_channel::<(usize, Result<RenderTile, RenderError>)>(capacity);
        let producer_session = session;
        let producer_cancel = cancellation.clone();
        #[cfg(feature = "native")]
        let producer_pool = self.pool.clone();
        let producer = thread::spawn(move || {
            #[cfg(feature = "native")]
            {
                use rayon::prelude::*;
                let compute = || {
                    tiles.into_par_iter().enumerate().for_each_with(
                        tile_sender,
                        |sender, (index, tile)| {
                            let result = producer_session
                                .compute_tile_with_cancel(tile, || producer_cancel.is_cancelled());
                            let _ = sender.send((index, result));
                        },
                    );
                };
                if let Some(pool) = producer_pool {
                    pool.install(compute);
                } else {
                    compute();
                }
            }
            #[cfg(not(feature = "native"))]
            {
                for (index, tile) in tiles.into_iter().enumerate() {
                    let result = producer_session
                        .compute_tile_with_cancel(tile, || producer_cancel.is_cancelled());
                    if tile_sender.send((index, result)).is_err() {
                        break;
                    }
                }
            }
        });

        let mut result = empty_result(&request);
        let mut completed_pixels = 0;
        let mut completed_tiles = 0;
        let mut next_index = 0;
        let mut pending = BTreeMap::new();
        let mut stopped = false;
        let mut terminal_sent = false;
        while let Ok((index, tile_result)) = tile_receiver.recv() {
            if ordered {
                pending.insert(index, tile_result);
                while let Some(tile_result) = pending.remove(&next_index) {
                    match deliver_tile(
                        sender,
                        &mut result,
                        tile_result,
                        request.generation,
                        total_pixels,
                        &mut completed_pixels,
                        &mut completed_tiles,
                        self.progress_every_tiles,
                    ) {
                        Delivery::Continue => {}
                        Delivery::Terminal => {
                            stopped = true;
                            terminal_sent = true;
                            break;
                        }
                        Delivery::Closed => {
                            stopped = true;
                            break;
                        }
                    }
                    next_index += 1;
                }
            } else {
                match deliver_tile(
                    sender,
                    &mut result,
                    tile_result,
                    request.generation,
                    total_pixels,
                    &mut completed_pixels,
                    &mut completed_tiles,
                    self.progress_every_tiles,
                ) {
                    Delivery::Continue => {}
                    Delivery::Terminal => {
                        stopped = true;
                        terminal_sent = true;
                    }
                    Delivery::Closed => stopped = true,
                }
            }
            if stopped {
                break;
            }
        }
        drop(tile_receiver);
        let _ = producer.join();
        if !terminal_sent && (stopped || cancellation.is_cancelled()) {
            let _ = sender.send(RenderEvent::Cancelled { generation: request.generation });
        } else if completed_pixels == total_pixels {
            let _ = sender.send(RenderEvent::Finished(result));
        }
    }
}

#[cfg(feature = "native")]
fn compute_frame_view(
    session: &RenderSession,
    view: RawFrameView<'_>,
    first_row: u32,
    cancellation: &CancellationToken,
) -> Result<(), RenderError> {
    if cancellation.is_cancelled() {
        return Err(RenderError::Cancelled);
    }
    if view.height() <= 32 {
        let mut view = view;
        return session
            .compute_view_with_cancel(&mut view, first_row, || cancellation.is_cancelled());
    }
    let split = view.height() / 2;
    let (top, bottom) = view.split_rows(split)?;
    let (top_result, bottom_result) = rayon::join(
        || compute_frame_view(session, top, first_row, cancellation),
        || compute_frame_view(session, bottom, first_row + split, cancellation),
    );
    top_result.and(bottom_result)
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

enum Delivery {
    Continue,
    Terminal,
    Closed,
}

#[allow(clippy::too_many_arguments)]
fn deliver_tile(
    sender: &SyncSender<RenderEvent>,
    result: &mut RenderResult,
    tile_result: Result<RenderTile, RenderError>,
    generation: u64,
    total_pixels: u64,
    completed_pixels: &mut u64,
    completed_tiles: &mut usize,
    progress_every_tiles: usize,
) -> Delivery {
    let tile = match tile_result {
        Ok(tile) => tile,
        Err(RenderError::Cancelled) => {
            return if sender.send(RenderEvent::Cancelled { generation }).is_ok() {
                Delivery::Terminal
            } else {
                Delivery::Closed
            };
        }
        Err(error) => {
            return if sender.send(RenderEvent::Failed(error)).is_ok() {
                Delivery::Terminal
            } else {
                Delivery::Closed
            };
        }
    };
    merge_tile(result, &tile);
    *completed_pixels += u64::from(tile.tile.width) * u64::from(tile.tile.height);
    *completed_tiles += 1;
    if sender.send(RenderEvent::TileReady(tile)).is_err() {
        return Delivery::Closed;
    }
    if (completed_tiles.is_multiple_of(progress_every_tiles) || *completed_pixels == total_pixels)
        && sender
            .send(RenderEvent::Progress(RenderProgress {
                generation,
                completed_pixels: *completed_pixels,
                total_pixels,
            }))
            .is_err()
    {
        return Delivery::Closed;
    }
    Delivery::Continue
}

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

    #[cfg(feature = "native")]
    #[test]
    fn full_frame_path_matches_progressive_result() {
        let executor = RenderExecutor::with_workers(4, 4, 3).unwrap();
        let token = CancellationToken::for_generation(7);
        let full = executor.compute_frame_full(&request(), &token).unwrap();
        let progressive = executor.render(&request(), &token);
        let finished = progressive
            .into_iter()
            .find_map(|event| match event {
                RenderEvent::Finished(result) => Some(result),
                _ => None,
            })
            .unwrap();
        assert_eq!(full, finished);
    }
}
