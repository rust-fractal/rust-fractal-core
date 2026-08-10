//! Deterministic, UI-independent fractal computation.
//!
//! This crate intentionally has no configuration-file, image, filesystem, thread,
//! or native numerical-backend dependency. It is the portable boundary that a
//! native executor and a future frontend can consume.

mod colour;
mod domain;
mod perturbation;
mod reference;
mod render;
mod root;
mod series;

pub use colour::{ColourError, ColourPipeline, ColourSettings, ImageFrame, Palette};
pub use domain::{
    ApproximationSettings, ComplexCoordinate, FractalFormula, GlitchSettings, OutputChannels,
    PrecisionPolicy, RenderError, RenderRequest, RenderSettings, SamplingSettings, Viewport, Zoom,
};
pub use fractal_numerics::{ComplexExtended, Extended};
pub use perturbation::{
    GlitchDetector, GlitchResult, PerturbationResult, perturbate, perturbate_extended,
    perturbate_with_series,
};
pub use reference::{
    ReferenceCheckpoint, ReferenceOrbit, ReferenceStatistics, generate_reference,
    generate_reference_exact,
};
pub use render::{
    PixelRect, RawFrame, RawFrameView, RawTile, RenderEvent, RenderProgress, RenderResult,
    RenderSession, RenderStarted, RenderTile, Tile, compute_frame, compute_tile, pixel_coordinate,
};
pub use root::{RootError, RootResult, find_nucleus, period_residual};
pub use series::{SeriesApproximation, SeriesCoefficients, generate_series};
