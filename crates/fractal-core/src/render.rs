use crate::perturbation::{perturbate_extended_with_cancel, perturbate_with_series_cancel};
use crate::{
    ComplexCoordinate, FractalFormula, GlitchDetector, ReferenceOrbit, RenderError, RenderRequest,
    SeriesApproximation, generate_reference, generate_reference_exact, generate_series,
};
use fractal_numerics::{ComplexExtended, Extended};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Tile {
    pub x: u32,
    pub y: u32,
    pub width: u32,
    pub height: u32,
}

pub type PixelRect = Tile;

impl Tile {
    pub fn validate(self, request: &RenderRequest) -> Result<(), RenderError> {
        let x_end = self.x.checked_add(self.width).ok_or(RenderError::InvalidTile)?;
        let y_end = self.y.checked_add(self.height).ok_or(RenderError::InvalidTile)?;
        if self.width == 0
            || self.height == 0
            || x_end > request.viewport.width
            || y_end > request.viewport.height
        {
            Err(RenderError::InvalidTile)
        } else {
            Ok(())
        }
    }

    fn pixel_count(self) -> Result<usize, RenderError> {
        (self.width as usize)
            .checked_mul(self.height as usize)
            .ok_or(RenderError::InvalidDimensions { width: self.width, height: self.height })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RenderStarted {
    pub generation: u64,
    pub width: u32,
    pub height: u32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RenderProgress {
    pub generation: u64,
    pub completed_pixels: u64,
    pub total_pixels: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RawTile {
    pub generation: u64,
    pub tile: Tile,
    pub iteration_limit: u32,
    pub iterations: Vec<u32>,
    pub smooth_iteration: Option<Vec<f64>>,
    pub distance: Option<Vec<f64>>,
    pub stripes: Option<Vec<f64>>,
    pub validity: Option<Vec<u8>>,
}

pub type RenderTile = RawTile;

#[derive(Clone, Debug, PartialEq)]
pub struct RawFrame {
    pub generation: u64,
    pub width: u32,
    pub height: u32,
    pub iteration_limit: u32,
    pub iterations: Vec<u32>,
    pub smooth_iteration: Option<Vec<f64>>,
    pub distance: Option<Vec<f64>>,
    pub stripes: Option<Vec<f64>>,
    pub validity: Option<Vec<u8>>,
}

pub type RenderResult = RawFrame;

/// A safe disjoint mutable view over a set of complete frame rows. Native
/// schedulers can split this view between workers without tile allocations or
/// tile-to-frame copies.
pub struct RawFrameView<'a> {
    width: u32,
    height: u32,
    iterations: &'a mut [u32],
    smooth_iteration: Option<&'a mut [f64]>,
    distance: Option<&'a mut [f64]>,
    stripes: Option<&'a mut [f64]>,
    validity: Option<&'a mut [u8]>,
}

impl RawFrame {
    pub fn view_mut(&mut self) -> RawFrameView<'_> {
        RawFrameView {
            width: self.width,
            height: self.height,
            iterations: &mut self.iterations,
            smooth_iteration: self.smooth_iteration.as_deref_mut(),
            distance: self.distance.as_deref_mut(),
            stripes: self.stripes.as_deref_mut(),
            validity: self.validity.as_deref_mut(),
        }
    }
}

impl<'a> RawFrameView<'a> {
    pub fn height(&self) -> u32 {
        self.height
    }

    pub fn split_rows(self, rows: u32) -> Result<(Self, Self), RenderError> {
        if rows == 0 || rows >= self.height {
            return Err(RenderError::InvalidRawFrame);
        }
        let split_at = rows as usize * self.width as usize;
        let (iterations_top, iterations_bottom) = self.iterations.split_at_mut(split_at);
        let split_optional =
            |values: Option<&'a mut [f64]>| values.map(|values| values.split_at_mut(split_at));
        let (smooth_top, smooth_bottom) = split_optional(self.smooth_iteration)
            .map_or((None, None), |(top, bottom)| (Some(top), Some(bottom)));
        let (distance_top, distance_bottom) = split_optional(self.distance)
            .map_or((None, None), |(top, bottom)| (Some(top), Some(bottom)));
        let (stripes_top, stripes_bottom) = split_optional(self.stripes)
            .map_or((None, None), |(top, bottom)| (Some(top), Some(bottom)));
        let (validity_top, validity_bottom) = self
            .validity
            .map(|values| {
                let (top, bottom) = values.split_at_mut(split_at);
                (Some(top), Some(bottom))
            })
            .unwrap_or((None, None));
        Ok((
            Self {
                width: self.width,
                height: rows,
                iterations: iterations_top,
                smooth_iteration: smooth_top,
                distance: distance_top,
                stripes: stripes_top,
                validity: validity_top,
            },
            Self {
                width: self.width,
                height: self.height - rows,
                iterations: iterations_bottom,
                smooth_iteration: smooth_bottom,
                distance: distance_bottom,
                stripes: stripes_bottom,
                validity: validity_bottom,
            },
        ))
    }
}

#[derive(Clone, Debug)]
struct ShallowTransform {
    top_left: ComplexCoordinate,
    x_step: ComplexCoordinate,
    y_step: ComplexCoordinate,
}

#[derive(Clone, Copy, Debug)]
struct DeepTransform {
    top_left: ComplexExtended,
    x_step: ComplexExtended,
    y_step: ComplexExtended,
}

#[derive(Clone, Debug)]
enum ViewportTransform {
    Shallow(ShallowTransform),
    Deep(DeepTransform),
}

#[derive(Clone, Debug)]
enum PixelLocation {
    Shallow(ComplexCoordinate),
    Deep(ComplexExtended),
}

impl ViewportTransform {
    fn new(request: &RenderRequest, deep: bool) -> Result<Self, RenderError> {
        let viewport = &request.viewport;
        let aspect = f64::from(viewport.width) / f64::from(viewport.height);
        let (sin, cos) = viewport.rotation_radians.sin_cos();
        if !deep {
            let scale = 4.0 / (f64::from(viewport.height) * viewport.zoom.scale());
            let x_step = ComplexCoordinate::new(scale * aspect * cos, scale * aspect * sin);
            let y_step = ComplexCoordinate::new(-scale * sin, scale * cos);
            let x = -f64::from(viewport.width) / 2.0;
            let y = -f64::from(viewport.height) / 2.0;
            let top_left = ComplexCoordinate::new(
                viewport.center.real + x * x_step.real + y * y_step.real,
                viewport.center.imag + x * x_step.imag + y * y_step.imag,
            );
            return Ok(Self::Shallow(ShallowTransform { top_left, x_step, y_step }));
        }

        let zoom = if let Some(decimal) = viewport.zoom.decimal() {
            Extended::from_decimal(decimal).map_err(|_| RenderError::InvalidDecimalZoom)?
        } else {
            Extended::from_f64(viewport.zoom.scale())
                .ok_or(RenderError::InvalidZoom { scale: viewport.zoom.scale() })?
        };
        let base_scale = Extended::from_f64(4.0 / f64::from(viewport.height))
            .ok_or(RenderError::NumericOverflow)?
            .checked_div(zoom)
            .ok_or(RenderError::NumericOverflow)?;
        let horizontal_scale = base_scale
            .checked_mul(Extended::from_f64(aspect).ok_or(RenderError::NumericOverflow)?)
            .ok_or(RenderError::NumericOverflow)?;
        let x_step = ComplexExtended::from_components(
            horizontal_scale
                .checked_mul(Extended::from_f64(cos).ok_or(RenderError::NumericOverflow)?)
                .ok_or(RenderError::NumericOverflow)?,
            horizontal_scale
                .checked_mul(Extended::from_f64(sin).ok_or(RenderError::NumericOverflow)?)
                .ok_or(RenderError::NumericOverflow)?,
        );
        let y_step = ComplexExtended::from_components(
            base_scale
                .checked_mul(Extended::from_f64(-sin).ok_or(RenderError::NumericOverflow)?)
                .ok_or(RenderError::NumericOverflow)?,
            base_scale
                .checked_mul(Extended::from_f64(cos).ok_or(RenderError::NumericOverflow)?)
                .ok_or(RenderError::NumericOverflow)?,
        );
        let top_left = x_step
            .scale(-f64::from(viewport.width) / 2.0)
            .ok_or(RenderError::NumericOverflow)?
            .add(
                y_step
                    .scale(-f64::from(viewport.height) / 2.0)
                    .ok_or(RenderError::NumericOverflow)?,
            )
            .ok_or(RenderError::NumericOverflow)?;
        Ok(Self::Deep(DeepTransform { top_left, x_step, y_step }))
    }

    fn location(&self, x: f64, y: f64) -> Result<PixelLocation, RenderError> {
        match self {
            Self::Shallow(transform) => Ok(PixelLocation::Shallow(ComplexCoordinate::new(
                transform.top_left.real + x * transform.x_step.real + y * transform.y_step.real,
                transform.top_left.imag + x * transform.x_step.imag + y * transform.y_step.imag,
            ))),
            Self::Deep(transform) => {
                let first = transform
                    .top_left
                    .add(transform.x_step.scale(x).ok_or(RenderError::NumericOverflow)?)
                    .ok_or(RenderError::NumericOverflow)?;
                let value = first
                    .add(transform.y_step.scale(y).ok_or(RenderError::NumericOverflow)?)
                    .ok_or(RenderError::NumericOverflow)?;
                Ok(PixelLocation::Deep(value))
            }
        }
    }
}

/// Immutable deterministic state shared by all tiles in one render.
///
/// Reference and series generation happen once when the session is created,
/// rather than once for every tile. The session owns no threads, channels, or
/// UI handles and can therefore be shared by a native executor.
pub struct RenderSession {
    request: RenderRequest,
    reference: Option<ReferenceOrbit>,
    series: Option<SeriesApproximation>,
    glitch_detector: Option<GlitchDetector>,
    transform: ViewportTransform,
    wants_smooth: bool,
    wants_distance: bool,
}

impl RenderSession {
    pub fn new(request: RenderRequest) -> Result<Self, RenderError> {
        request.validate()?;
        let deep = requires_deep_precision(&request);
        let reference = approximation_reference(&request, deep)?;
        let series = request
            .settings
            .approximation
            .enabled
            .then(|| reference.as_ref().map(|reference| approximation_series(&request, reference)))
            .flatten()
            .transpose()?;
        let transform = ViewportTransform::new(&request, deep)?;
        let glitch_detector = request
            .settings
            .glitch
            .enabled
            .then(|| GlitchDetector::new(request.settings.glitch.tolerance).ok())
            .flatten();
        Ok(Self {
            request: request.clone(),
            reference,
            series,
            glitch_detector,
            transform,
            wants_smooth: request.settings.output.smooth_iteration
                || request.settings.output.stripes,
            wants_distance: request.settings.output.distance,
        })
    }

    pub const fn request(&self) -> &RenderRequest {
        &self.request
    }

    /// Returns the pixel displacement from the exact reference center. Deep
    /// callers can inspect this value without converting it to an absolute
    /// f64 coordinate.
    pub fn pixel_offset(&self, x: f64, y: f64) -> Result<ComplexExtended, RenderError> {
        match self.transform.location(x, y)? {
            PixelLocation::Deep(offset) => Ok(offset),
            PixelLocation::Shallow(coordinate) => ComplexExtended::from_f64(
                coordinate.real - self.request.viewport.center.real,
                coordinate.imag - self.request.viewport.center.imag,
            )
            .ok_or(RenderError::NumericOverflow),
        }
    }

    pub fn compute_frame(&self) -> Result<RenderResult, RenderError> {
        let tile = Tile {
            x: 0,
            y: 0,
            width: self.request.viewport.width,
            height: self.request.viewport.height,
        };
        let tile_result = self.compute_tile(tile)?;
        Ok(RawFrame {
            generation: self.request.generation,
            width: self.request.viewport.width,
            height: self.request.viewport.height,
            iteration_limit: self.request.settings.iteration_limit,
            iterations: tile_result.iterations,
            smooth_iteration: tile_result.smooth_iteration,
            distance: tile_result.distance,
            stripes: tile_result.stripes,
            validity: tile_result.validity,
        })
    }

    pub fn compute_tile(&self, tile: Tile) -> Result<RenderTile, RenderError> {
        self.compute_tile_with_cancel(tile, || false)
    }

    /// Fills a disjoint set of frame rows directly. The caller owns the view,
    /// so parallel callers cannot alias output channels.
    pub fn compute_view_with_cancel<F>(
        &self,
        view: &mut RawFrameView<'_>,
        first_row: u32,
        mut cancelled: F,
    ) -> Result<(), RenderError>
    where
        F: FnMut() -> bool,
    {
        let sample_context = SampleContext {
            request: &self.request,
            reference: self.reference.as_ref(),
            series: self.series.as_ref(),
            glitch_detector: self.glitch_detector,
            transform: &self.transform,
            want_smooth: self.wants_smooth,
            want_distance: self.wants_distance,
        };
        for local_y in 0..view.height {
            for x in 0..view.width {
                if cancelled() {
                    return Err(RenderError::Cancelled);
                }
                let sample = sample_pixel(&sample_context, x, first_row + local_y, &mut cancelled)?;
                let index = (local_y * view.width + x) as usize;
                view.iterations[index] = sample.iterations;
                if let Some(channel) = &mut view.smooth_iteration {
                    channel[index] = sample.smooth_iteration;
                }
                if let Some(channel) = &mut view.distance {
                    channel[index] = sample.distance;
                }
                if let Some(channel) = &mut view.stripes {
                    channel[index] = sample.stripes;
                }
                if let Some(channel) = &mut view.validity {
                    channel[index] = u8::from(sample.valid);
                }
            }
        }
        Ok(())
    }

    /// Computes one tile while allowing the owning executor to stop between pixels.
    pub fn compute_tile_with_cancel<F>(
        &self,
        tile: Tile,
        mut cancelled: F,
    ) -> Result<RenderTile, RenderError>
    where
        F: FnMut() -> bool,
    {
        tile.validate(&self.request)?;

        let capacity = tile.pixel_count()?;
        let mut iterations = Vec::with_capacity(capacity);
        let mut smooth_iteration =
            self.request.settings.output.smooth_iteration.then(|| Vec::with_capacity(capacity));
        let mut distance =
            self.request.settings.output.distance.then(|| Vec::with_capacity(capacity));
        let mut stripes =
            self.request.settings.output.stripes.then(|| Vec::with_capacity(capacity));
        let mut validity =
            self.request.settings.output.validity.then(|| Vec::with_capacity(capacity));
        let sample_context = SampleContext {
            request: &self.request,
            reference: self.reference.as_ref(),
            series: self.series.as_ref(),
            glitch_detector: self.glitch_detector,
            transform: &self.transform,
            want_smooth: self.wants_smooth,
            want_distance: self.wants_distance,
        };
        for y in tile.y..tile.y + tile.height {
            for x in tile.x..tile.x + tile.width {
                if cancelled() {
                    return Err(RenderError::Cancelled);
                }
                let sample = sample_pixel(&sample_context, x, y, &mut cancelled)?;
                iterations.push(sample.iterations);
                if let Some(channel) = &mut smooth_iteration {
                    channel.push(sample.smooth_iteration);
                }
                if let Some(channel) = &mut distance {
                    channel.push(sample.distance);
                }
                if let Some(channel) = &mut stripes {
                    channel.push(sample.stripes);
                }
                if let Some(channel) = &mut validity {
                    channel.push(u8::from(sample.valid));
                }
            }
        }

        Ok(RawTile {
            generation: self.request.generation,
            tile,
            iteration_limit: self.request.settings.iteration_limit,
            iterations,
            smooth_iteration,
            distance,
            stripes,
            validity,
        })
    }
}

impl RawFrame {
    pub fn pixel_count(&self) -> Result<usize, RenderError> {
        (self.width as usize)
            .checked_mul(self.height as usize)
            .ok_or(RenderError::InvalidDimensions { width: self.width, height: self.height })
    }

    pub fn validate(&self) -> Result<(), RenderError> {
        if self.width == 0 || self.height == 0 || self.iteration_limit == 0 {
            return Err(RenderError::InvalidRawFrame);
        }
        let pixel_count = self.pixel_count()?;
        let valid = |length: usize| length == pixel_count;
        if !valid(self.iterations.len())
            || self.smooth_iteration.as_ref().is_some_and(|v| !valid(v.len()))
            || self.distance.as_ref().is_some_and(|v| !valid(v.len()))
            || self.stripes.as_ref().is_some_and(|v| !valid(v.len()))
            || self.validity.as_ref().is_some_and(|v| !valid(v.len()))
        {
            return Err(RenderError::InvalidRawFrame);
        }
        Ok(())
    }

    /// Stable, dependency-free hash for compact raw regression fixtures.
    pub fn raw_hash(&self) -> Result<u64, RenderError> {
        self.validate()?;
        let mut hash = 0xcbf29ce484222325_u64;
        let mut update = |bytes: &[u8]| {
            for byte in bytes {
                hash ^= u64::from(*byte);
                hash = hash.wrapping_mul(0x100000001b3);
            }
        };
        update(&self.generation.to_le_bytes());
        update(&self.width.to_le_bytes());
        update(&self.height.to_le_bytes());
        update(&self.iteration_limit.to_le_bytes());
        for iteration in &self.iterations {
            update(&iteration.to_le_bytes());
        }
        for channel in
            [self.smooth_iteration.as_deref(), self.distance.as_deref(), self.stripes.as_deref()]
                .into_iter()
                .flatten()
        {
            for value in channel {
                update(&value.to_bits().to_le_bytes());
            }
        }
        if let Some(channel) = self.validity.as_deref() {
            update(channel);
        }
        Ok(hash)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RenderEvent {
    Started(RenderStarted),
    Progress(RenderProgress),
    TileReady(RenderTile),
    Finished(RenderResult),
    Cancelled { generation: u64 },
    Failed(RenderError),
}

pub fn compute_frame(request: &RenderRequest) -> Result<RenderResult, RenderError> {
    RenderSession::new(request.clone())?.compute_frame()
}

pub fn compute_tile(request: &RenderRequest, tile: Tile) -> Result<RenderTile, RenderError> {
    RenderSession::new(request.clone())?.compute_tile(tile)
}

fn requires_deep_precision(request: &RenderRequest) -> bool {
    request
        .settings
        .precision
        .required_bits(request.viewport.zoom.clone())
        .is_ok_and(|bits| bits > 128)
}

fn approximation_reference(
    request: &RenderRequest,
    deep: bool,
) -> Result<Option<ReferenceOrbit>, RenderError> {
    if !deep && !request.settings.approximation.enabled {
        return Ok(None);
    }
    let reference = if deep {
        let precision = request.settings.precision.required_bits(request.viewport.zoom.clone())?;
        generate_reference_exact(
            &request.viewport.center,
            request.settings.fractal,
            precision,
            request.settings.iteration_limit,
        )?
    } else {
        let reference = generate_reference(
            request.viewport.center.clone(),
            request.settings.fractal,
            request.settings.iteration_limit,
        )?;
        if reference.statistics.generated_iterations != request.settings.iteration_limit {
            return Ok(None);
        }
        reference
    };
    if deep && reference.statistics.generated_iterations != request.settings.iteration_limit {
        return Err(RenderError::ReferenceGeneration(
            "deep reference escaped before the requested iteration limit".to_owned(),
        ));
    }
    Ok(Some(reference))
}

fn approximation_series(
    request: &RenderRequest,
    reference: &ReferenceOrbit,
) -> Result<SeriesApproximation, RenderError> {
    generate_series(reference, request.settings.approximation.order)
}

struct SampleContext<'a> {
    request: &'a RenderRequest,
    reference: Option<&'a ReferenceOrbit>,
    series: Option<&'a SeriesApproximation>,
    glitch_detector: Option<GlitchDetector>,
    transform: &'a ViewportTransform,
    want_smooth: bool,
    want_distance: bool,
}

fn sample_pixel<F>(
    context: &SampleContext<'_>,
    x: u32,
    y: u32,
    cancelled: &mut F,
) -> Result<EscapeSample, RenderError>
where
    F: FnMut() -> bool,
{
    let request = context.request;
    let samples = request.settings.sampling.samples_per_pixel;
    if samples == 1 {
        return approximated_sample(context, f64::from(x) + 0.5, f64::from(y) + 0.5, cancelled);
    }

    let grid = (f64::from(samples).sqrt().ceil() as u32).max(1);
    let mut iteration_sum = 0.0;
    let mut total = EscapeSample {
        iterations: 0,
        smooth_iteration: 0.0,
        distance: 0.0,
        stripes: 0.0,
        valid: true,
    };
    for sample in 0..samples {
        let sample_x = sample % grid;
        let sample_y = sample / grid;
        let coordinate_x = f64::from(x) + (f64::from(sample_x) + 0.5) / f64::from(grid);
        let coordinate_y = f64::from(y) + (f64::from(sample_y) + 0.5) / f64::from(grid);
        let value = approximated_sample(context, coordinate_x, coordinate_y, cancelled)?;
        iteration_sum += f64::from(value.iterations);
        total.smooth_iteration += value.smooth_iteration;
        total.distance += value.distance;
        total.stripes += value.stripes;
        total.valid &= value.valid;
    }
    let divisor = f64::from(samples);
    total.iterations = (iteration_sum / divisor).round() as u32;
    total.smooth_iteration /= divisor;
    total.distance /= divisor;
    total.stripes /= divisor;
    Ok(total)
}

fn approximated_sample<F>(
    context: &SampleContext<'_>,
    x: f64,
    y: f64,
    cancelled: &mut F,
) -> Result<EscapeSample, RenderError>
where
    F: FnMut() -> bool,
{
    let request = context.request;
    let reference = context.reference;
    let series = context.series;
    let glitch_detector = context.glitch_detector;
    let want_smooth = context.want_smooth;
    let want_distance = context.want_distance;
    let location = context.transform.location(x, y)?;
    if let PixelLocation::Deep(delta_c) = location {
        let reference = reference.ok_or_else(|| {
            RenderError::ReferenceGeneration("deep render has no reference orbit".to_owned())
        })?;
        let result = perturbate_extended_with_cancel(
            reference,
            delta_c,
            request.settings.iteration_limit,
            glitch_detector,
            context.series,
            &mut *cancelled,
        )?;
        return Ok(EscapeSample {
            iterations: result.iterations,
            smooth_iteration: if want_smooth {
                result.smooth_iteration
            } else {
                result.iterations as f64
            },
            distance: if want_distance { result.distance } else { 0.0 },
            stripes: if want_smooth {
                0.5 + 0.5 * (result.smooth_iteration * core::f64::consts::TAU).sin()
            } else {
                0.0
            },
            valid: !result.glitch.detected || result.glitch.corrected,
        });
    }
    let PixelLocation::Shallow(coordinate) = location else {
        return Err(RenderError::NumericOverflow);
    };
    if let Some(reference) = reference
        && let Ok(result) = perturbate_with_series_cancel(
            reference,
            coordinate.clone(),
            request.settings.iteration_limit,
            glitch_detector,
            series,
            &mut *cancelled,
        )
    {
        if result.glitch.detected {
            let mut fallback = escape_sample_with_cancel(
                coordinate,
                request.settings.fractal,
                request.settings.iteration_limit,
                want_smooth,
                want_distance,
                &mut *cancelled,
            )?;
            fallback.valid = false;
            return Ok(fallback);
        }
        return Ok(EscapeSample {
            iterations: result.iterations,
            smooth_iteration: if want_smooth {
                result.smooth_iteration
            } else {
                result.iterations as f64
            },
            distance: if want_distance { result.distance } else { 0.0 },
            stripes: if want_smooth {
                0.5 + 0.5 * (result.smooth_iteration * core::f64::consts::TAU).sin()
            } else {
                0.0
            },
            valid: true,
        });
    }
    escape_sample_with_cancel(
        coordinate,
        request.settings.fractal,
        request.settings.iteration_limit,
        want_smooth,
        want_distance,
        &mut *cancelled,
    )
}

pub fn pixel_coordinate(request: &RenderRequest, x: u32, y: u32) -> ComplexCoordinate {
    pixel_coordinate_at(request, f64::from(x) + 0.5, f64::from(y) + 0.5)
}

fn pixel_coordinate_at(request: &RenderRequest, x: f64, y: f64) -> ComplexCoordinate {
    let viewport = &request.viewport;
    let aspect = viewport.width as f64 / viewport.height as f64;
    let scale = 4.0 / (viewport.height as f64 * viewport.zoom.scale());
    let real = (x - viewport.width as f64 / 2.0) * scale * aspect;
    let imag = (y - viewport.height as f64 / 2.0) * scale;
    let (sin, cos) = viewport.rotation_radians.sin_cos();

    ComplexCoordinate::new(
        viewport.center.real + real * cos - imag * sin,
        viewport.center.imag + real * sin + imag * cos,
    )
}

#[derive(Clone, Copy, Debug)]
struct EscapeSample {
    iterations: u32,
    smooth_iteration: f64,
    distance: f64,
    stripes: f64,
    valid: bool,
}

fn escape_sample_with_cancel<F>(
    coordinate: ComplexCoordinate,
    formula: FractalFormula,
    limit: u32,
    want_smooth_iteration: bool,
    want_distance: bool,
    mut cancelled: F,
) -> Result<EscapeSample, RenderError>
where
    F: FnMut() -> bool,
{
    let mut z = ComplexCoordinate::new(0.0, 0.0);
    let mut derivative = ComplexCoordinate::new(0.0, 0.0);
    let power = match formula {
        FractalFormula::Mandelbrot { power } | FractalFormula::BurningShip { power } => power,
    };
    let mut iteration = 0;
    while iteration < limit && z.real * z.real + z.imag * z.imag <= 4.0 {
        if iteration.is_multiple_of(128) && cancelled() {
            return Err(RenderError::Cancelled);
        }
        let (mut real, mut imag) = (z.real, z.imag);
        if matches!(formula, FractalFormula::BurningShip { .. }) {
            real = real.abs();
            imag = imag.abs();
        }

        let powered = complex_power(real, imag, power);
        if want_distance_iteration(formula, want_distance) {
            let derivative_power = complex_power(real, imag, power.saturating_sub(1));
            let scaled = (derivative_power.0 * power as f64, derivative_power.1 * power as f64);
            let next_derivative = complex_multiply(scaled, (derivative.real, derivative.imag));
            derivative = ComplexCoordinate::new(next_derivative.0 + 1.0, next_derivative.1);
        }
        z = ComplexCoordinate::new(powered.0 + coordinate.real, powered.1 + coordinate.imag);
        iteration += 1;
    }

    let magnitude = (z.real * z.real + z.imag * z.imag).sqrt();
    let escaped = iteration < limit;
    let smooth_iteration = if want_smooth_iteration && escaped && magnitude > 1.0 {
        iteration as f64 + 1.0 - magnitude.ln().ln() / (power as f64).ln()
    } else {
        iteration as f64
    };
    let derivative_magnitude =
        (derivative.real * derivative.real + derivative.imag * derivative.imag).sqrt();
    let distance = if want_distance && escaped && derivative_magnitude > 0.0 {
        0.5 * magnitude.ln() * magnitude / derivative_magnitude
    } else {
        0.0
    };
    let stripes = if want_smooth_iteration {
        0.5 + 0.5 * (smooth_iteration * core::f64::consts::TAU).sin()
    } else {
        0.0
    };

    Ok(EscapeSample {
        iterations: iteration,
        smooth_iteration,
        distance,
        stripes,
        valid: magnitude.is_finite() && z.real.is_finite() && z.imag.is_finite(),
    })
}

fn want_distance_iteration(formula: FractalFormula, want_distance: bool) -> bool {
    want_distance && matches!(formula, FractalFormula::Mandelbrot { .. })
}

fn complex_power(real: f64, imag: f64, power: u32) -> (f64, f64) {
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

fn complex_multiply(left: (f64, f64), right: (f64, f64)) -> (f64, f64) {
    (left.0 * right.0 - left.1 * right.1, left.0 * right.1 + left.1 * right.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ApproximationSettings, GlitchSettings, OutputChannels, PrecisionPolicy, RenderSettings,
        SamplingSettings, Viewport, Zoom,
    };

    fn request() -> RenderRequest {
        RenderRequest {
            generation: 42,
            viewport: Viewport {
                center: ComplexCoordinate::new(-0.75, 0.0),
                zoom: Zoom::unit(),
                rotation_radians: 0.0,
                width: 16,
                height: 12,
            },
            settings: RenderSettings {
                fractal: FractalFormula::Mandelbrot { power: 2 },
                iteration_limit: 100,
                approximation: ApproximationSettings::default(),
                glitch: GlitchSettings::default(),
                sampling: SamplingSettings::default(),
                output: OutputChannels::default(),
                precision: PrecisionPolicy::default(),
            },
        }
    }

    #[test]
    fn frame_computation_is_deterministic() {
        let request = request();
        assert_eq!(compute_frame(&request), compute_frame(&request));
    }

    #[test]
    fn tiles_match_their_frame_region() {
        let request = request();
        let frame = compute_frame(&request).unwrap();
        let tile = compute_tile(&request, Tile { x: 4, y: 3, width: 5, height: 4 }).unwrap();

        for local_y in 0..tile.tile.height {
            for local_x in 0..tile.tile.width {
                let tile_index = (local_y * tile.tile.width + local_x) as usize;
                let frame_index = ((tile.tile.y + local_y) * request.viewport.width
                    + tile.tile.x
                    + local_x) as usize;
                assert_eq!(tile.iterations[tile_index], frame.iterations[frame_index]);
            }
        }
    }

    #[test]
    fn requested_channels_are_allocated_and_aligned() {
        let mut request = request();
        request.settings.output = OutputChannels {
            smooth_iteration: true,
            distance: true,
            stripes: true,
            validity: true,
            ..OutputChannels::default()
        };
        let frame = compute_frame(&request).unwrap();
        assert_eq!(frame.iterations.len(), 192);
        assert_eq!(frame.smooth_iteration.as_ref().unwrap().len(), 192);
        assert_eq!(frame.distance.as_ref().unwrap().len(), 192);
        assert_eq!(frame.stripes.as_ref().unwrap().len(), 192);
        assert_eq!(frame.validity.as_ref().unwrap().len(), 192);
    }

    #[test]
    fn invalid_tiles_are_rejected() {
        let request = request();
        assert_eq!(
            compute_tile(&request, Tile { x: 15, y: 0, width: 2, height: 1 }),
            Err(RenderError::InvalidTile)
        );
    }

    #[test]
    fn tile_computation_honours_cancellation_between_pixels() {
        let session = RenderSession::new(request()).unwrap();
        assert_eq!(
            session.compute_tile_with_cancel(Tile { x: 0, y: 0, width: 4, height: 4 }, || true,),
            Err(RenderError::Cancelled)
        );
    }

    #[test]
    fn deep_requests_keep_adjacent_pixel_offsets_distinct() {
        let request = RenderRequest {
            generation: 1,
            viewport: Viewport {
                center: ComplexCoordinate::from_decimal("-0.743643887037151", "0.13182590420533")
                    .unwrap(),
                zoom: Zoom::from_decimal("1E8000").unwrap(),
                rotation_radians: 0.0,
                width: 3,
                height: 1,
            },
            settings: RenderSettings {
                iteration_limit: 32,
                output: OutputChannels { validity: true, ..OutputChannels::default() },
                ..RenderSettings::default()
            },
        };
        let session = RenderSession::new(request).unwrap();
        let first = session.pixel_offset(0.5, 0.5).unwrap();
        let second = session.pixel_offset(1.5, 0.5).unwrap();
        assert_ne!(first, second);
        let frame = session.compute_frame().unwrap();
        assert_eq!(frame.iterations.len(), 3);
        assert!(frame.validity.as_ref().unwrap().iter().all(|&valid| valid == 1));
    }

    #[test]
    fn exact_moderate_zoom_matches_the_shallow_reference() {
        let mut exact = request();
        exact.viewport.center = ComplexCoordinate::from_decimal("-0.75", "0").unwrap();
        exact.viewport.zoom = Zoom::from_decimal("1E2").unwrap();
        exact.viewport.width = 8;
        exact.viewport.height = 6;
        exact.settings.iteration_limit = 80;
        let shallow = RenderRequest {
            viewport: Viewport {
                center: ComplexCoordinate::new(-0.75, 0.0),
                zoom: Zoom::new(100.0).unwrap(),
                ..exact.viewport.clone()
            },
            ..exact.clone()
        };
        let exact_frame = compute_frame(&exact).unwrap();
        let shallow_frame = compute_frame(&shallow).unwrap();
        assert_eq!(exact_frame.iterations, shallow_frame.iterations);
    }

    #[test]
    fn deep_specialized_formulas_produce_frames() {
        for formula in
            [FractalFormula::Mandelbrot { power: 3 }, FractalFormula::BurningShip { power: 2 }]
        {
            let request = RenderRequest {
                generation: 1,
                viewport: Viewport {
                    center: ComplexCoordinate::from_decimal("0", "0").unwrap(),
                    zoom: Zoom::from_decimal("1E50").unwrap(),
                    rotation_radians: 0.0,
                    width: 3,
                    height: 2,
                },
                settings: RenderSettings {
                    fractal: formula,
                    iteration_limit: 24,
                    ..RenderSettings::default()
                },
            };
            let frame = compute_frame(&request).unwrap();
            assert_eq!(frame.iterations.len(), 6);
        }
    }
}
