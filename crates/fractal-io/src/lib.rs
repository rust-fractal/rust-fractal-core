//! File-format boundaries for the native application layer.

use fractal_core::{
    ApproximationSettings, ComplexCoordinate, FractalFormula, GlitchSettings, OutputChannels,
    RenderSettings, SamplingSettings, Viewport, Zoom,
};
#[cfg(any(feature = "png-export", feature = "jpeg-export", feature = "exr-export"))]
use fractal_core::{ImageFrame, RenderResult};
use fractal_core::{RenderError, RenderRequest};
#[cfg(any(feature = "png-export", feature = "jpeg-export", feature = "exr-export"))]
use std::path::Path;
use std::{error::Error, fmt};

#[derive(Debug)]
pub enum IoError {
    InvalidRequest(RenderError),
    Json(serde_json::Error),
    TomlDeserialize(toml::de::Error),
    TomlSerialize(toml::ser::Error),
    InvalidLegacyNumber(String),
    #[cfg(any(feature = "png-export", feature = "jpeg-export"))]
    Image(image::ImageError),
    #[cfg(feature = "exr-export")]
    Exr(exr::error::Error),
    BufferSize,
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRequest(error) => write!(f, "invalid render request: {error}"),
            Self::Json(error) => write!(f, "settings JSON error: {error}"),
            Self::TomlDeserialize(error) => write!(f, "settings TOML error: {error}"),
            Self::TomlSerialize(error) => write!(f, "settings TOML serialization error: {error}"),
            Self::InvalidLegacyNumber(value) => {
                write!(f, "legacy numeric value is invalid: {value}")
            }
            #[cfg(feature = "png-export")]
            Self::Image(error) => write!(f, "image export error: {error}"),
            #[cfg(feature = "exr-export")]
            Self::Exr(error) => write!(f, "EXR export error: {error}"),
            Self::BufferSize => f.write_str("render buffer size does not match its dimensions"),
        }
    }
}

impl Error for IoError {}

pub fn request_to_json(request: &RenderRequest) -> Result<String, IoError> {
    request.validate().map_err(IoError::InvalidRequest)?;
    serde_json::to_string_pretty(request).map_err(IoError::Json)
}

pub fn request_from_json(json: &str) -> Result<RenderRequest, IoError> {
    let request: RenderRequest = serde_json::from_str(json).map_err(IoError::Json)?;
    request.validate().map_err(IoError::InvalidRequest)?;
    Ok(request)
}

pub fn request_to_toml(request: &RenderRequest) -> Result<String, IoError> {
    request.validate().map_err(IoError::InvalidRequest)?;
    toml::to_string_pretty(request).map_err(IoError::TomlSerialize)
}

pub fn request_from_toml(toml_text: &str) -> Result<RenderRequest, IoError> {
    let request: RenderRequest = toml::from_str(toml_text).map_err(IoError::TomlDeserialize)?;
    request.validate().map_err(IoError::InvalidRequest)?;
    Ok(request)
}

#[derive(Clone, Debug, PartialEq, Eq, serde::Deserialize)]
pub struct LegacyLocation {
    pub real: String,
    pub imag: String,
    pub zoom: String,
    pub iterations: u32,
}

impl LegacyLocation {
    /// Returns the positive scientific-notation exponent, when the zoom has one.
    ///
    /// This intentionally inspects the original string rather than parsing it as
    /// `f64`, so deep locations remain available to an arbitrary-precision caller.
    pub fn zoom_exponent(&self) -> Option<u32> {
        let index = self.zoom.bytes().position(|byte| byte == b'e' || byte == b'E')?;
        self.zoom.get(index + 1..)?.parse::<u32>().ok()
    }

    /// Estimates the Malachite precision needed for this decimal zoom.
    pub fn required_precision_bits(&self, guard_bits: u32) -> Result<u32, IoError> {
        let exponent = u64::from(self.zoom_exponent().unwrap_or(0));
        let bits = 256_u64
            .checked_add(exponent.saturating_mul(3322).div_ceil(1000))
            .and_then(|bits| bits.checked_add(u64::from(guard_bits)))
            .ok_or_else(|| IoError::InvalidLegacyNumber("zoom precision overflow".to_owned()))?;
        u32::try_from(bits)
            .map_err(|_| IoError::InvalidLegacyNumber("zoom precision exceeds u32".to_owned()))
    }
}

#[derive(Debug, Default, serde::Deserialize)]
struct LegacyOptions {
    image_width: Option<u32>,
    image_height: Option<u32>,
    rotate: Option<f64>,
    approximation_order: Option<u32>,
    series_approximation_enabled: Option<bool>,
    glitch_tolerance: Option<f64>,
    distance_color: Option<bool>,
    stripe_scale: Option<f64>,
    display_glitches: Option<bool>,
}

pub fn request_from_legacy_toml(
    options_toml: &str,
    location_toml: &str,
    generation: u64,
) -> Result<RenderRequest, IoError> {
    let options: LegacyOptions = toml::from_str(options_toml).map_err(IoError::TomlDeserialize)?;
    let location = legacy_location_from_toml(location_toml)?;
    let center = ComplexCoordinate::from_decimal(&location.real, &location.imag)
        .map_err(IoError::InvalidRequest)?;
    let zoom = Zoom::from_decimal(&location.zoom).map_err(IoError::InvalidRequest)?;
    let viewport = Viewport {
        center,
        zoom,
        rotation_radians: options.rotate.unwrap_or(0.0).to_radians(),
        width: options.image_width.unwrap_or(3840),
        height: options.image_height.unwrap_or(2160),
    };
    let output = OutputChannels {
        distance: options.distance_color.unwrap_or(false),
        smooth_iteration: options.stripe_scale.is_some(),
        stripes: options.stripe_scale.is_some(),
        validity: options.display_glitches.unwrap_or(false),
        ..OutputChannels::default()
    };
    let settings = RenderSettings {
        fractal: FractalFormula::default(),
        iteration_limit: location.iterations,
        approximation: ApproximationSettings {
            enabled: options.series_approximation_enabled.unwrap_or(false),
            order: options.approximation_order.unwrap_or(0),
        },
        glitch: GlitchSettings {
            enabled: true,
            tolerance: options.glitch_tolerance.unwrap_or(1.4e-6),
        },
        sampling: SamplingSettings::default(),
        output,
        precision: Default::default(),
    };
    let request = RenderRequest { generation, viewport, settings };
    request.validate().map_err(IoError::InvalidRequest)?;
    Ok(request)
}

/// Parses a historical location while preserving all numeric fields verbatim.
///
/// Use this function for arbitrary-precision reference calculations. The
/// compatibility request converter above remains intentionally `f64`-based.
pub fn legacy_location_from_toml(location_toml: &str) -> Result<LegacyLocation, IoError> {
    toml::from_str(location_toml).map_err(IoError::TomlDeserialize)
}

#[derive(Debug, serde::Deserialize)]
struct LegacyPalette {
    palette: Vec<u8>,
    palette_iteration_span: Option<f64>,
}

pub fn palette_from_toml(palette_toml: &str) -> Result<(fractal_core::Palette, f64), IoError> {
    let palette: LegacyPalette = toml::from_str(palette_toml).map_err(IoError::TomlDeserialize)?;
    if !palette.palette.len().is_multiple_of(3) || palette.palette.is_empty() {
        return Err(IoError::InvalidLegacyNumber("palette length".to_owned()));
    }
    let stops = palette
        .palette
        .chunks_exact(3)
        .map(|rgb| [rgb[0] as f32 / 255.0, rgb[1] as f32 / 255.0, rgb[2] as f32 / 255.0, 1.0])
        .collect();
    let result = fractal_core::Palette::new(stops, true)
        .map_err(|_| IoError::InvalidLegacyNumber("palette values".to_owned()))?;
    Ok((result, palette.palette_iteration_span.unwrap_or(100.0)))
}

#[cfg(feature = "png-export")]
pub fn write_png(path: impl AsRef<Path>, result: &RenderResult) -> Result<(), IoError> {
    let expected =
        (result.width as usize).checked_mul(result.height as usize).ok_or(IoError::BufferSize)?;
    if result.iterations.len() != expected {
        return Err(IoError::BufferSize);
    }
    let maximum = result.iterations.iter().copied().max().unwrap_or(1).max(1);
    let pixels = result
        .iterations
        .iter()
        .map(|&iteration| ((u64::from(iteration) * 255) / u64::from(maximum)) as u8)
        .collect::<Vec<_>>();
    let image =
        image::ImageBuffer::<image::Luma<u8>, _>::from_vec(result.width, result.height, pixels)
            .ok_or(IoError::BufferSize)?;
    image.save(path).map_err(IoError::Image)
}

#[cfg(feature = "png-export")]
pub fn write_rgba_png(path: impl AsRef<Path>, image_frame: &ImageFrame) -> Result<(), IoError> {
    let expected = (image_frame.width as usize)
        .checked_mul(image_frame.height as usize)
        .and_then(|pixels| pixels.checked_mul(4))
        .ok_or(IoError::BufferSize)?;
    if image_frame.rgba8.len() != expected {
        return Err(IoError::BufferSize);
    }
    let image = image::ImageBuffer::<image::Rgba<u8>, _>::from_vec(
        image_frame.width,
        image_frame.height,
        image_frame.rgba8.clone(),
    )
    .ok_or(IoError::BufferSize)?;
    image.save(path).map_err(IoError::Image)
}

#[cfg(feature = "jpeg-export")]
pub fn write_jpeg(path: impl AsRef<Path>, image_frame: &ImageFrame) -> Result<(), IoError> {
    let expected = (image_frame.width as usize)
        .checked_mul(image_frame.height as usize)
        .and_then(|pixels| pixels.checked_mul(4))
        .ok_or(IoError::BufferSize)?;
    if image_frame.rgba8.len() != expected {
        return Err(IoError::BufferSize);
    }
    let image = image::ImageBuffer::<image::Rgba<u8>, _>::from_vec(
        image_frame.width,
        image_frame.height,
        image_frame.rgba8.clone(),
    )
    .ok_or(IoError::BufferSize)?;
    image::DynamicImage::ImageRgba8(image)
        .save_with_format(path, image::ImageFormat::Jpeg)
        .map_err(IoError::Image)
}

#[cfg(feature = "exr-export")]
pub fn write_exr(path: impl AsRef<Path>, result: &RenderResult) -> Result<(), IoError> {
    let expected =
        (result.width as usize).checked_mul(result.height as usize).ok_or(IoError::BufferSize)?;
    if result.iterations.len() != expected {
        return Err(IoError::BufferSize);
    }
    let maximum = result.iterations.iter().copied().max().unwrap_or(1).max(1) as f32;
    exr::prelude::write_rgb_file(path, result.width as usize, result.height as usize, |x, y| {
        let value = result.iterations[y * result.width as usize + x] as f32 / maximum;
        (value, value, value)
    })
    .map_err(IoError::Exr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fractal_core::{ComplexCoordinate, Viewport, Zoom};

    #[test]
    fn settings_round_trip_through_json() {
        let request = RenderRequest {
            generation: 9,
            viewport: Viewport {
                center: ComplexCoordinate::new(-0.75, 0.1),
                zoom: Zoom::unit(),
                rotation_radians: 0.0,
                width: 32,
                height: 24,
            },
            settings: Default::default(),
        };
        let restored = request_from_json(&request_to_json(&request).unwrap()).unwrap();
        assert_eq!(restored, request);
    }

    #[test]
    fn historical_split_toml_converts_to_a_valid_request() {
        let options = "image_width = 32\nimage_height = 24\nrotate = 90\nseries_approximation_enabled = true\napproximation_order = 8\n";
        let location = "real = \"-0.75\"\nimag = \"0.0\"\nzoom = \"1E2\"\niterations = 200\n";
        let request = request_from_legacy_toml(options, location, 12).unwrap();
        assert_eq!(request.viewport.width, 32);
        assert!((request.viewport.rotation_radians - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
        assert_eq!(request.settings.iteration_limit, 200);
        assert_eq!(request.settings.approximation.order, 8);
    }

    #[test]
    fn deep_legacy_location_preserves_decimal_strings_for_arbitrary_precision() {
        let location = legacy_location_from_toml(
            "real = \"-0.743643887037151\"\nimag = \"0.13182590420533\"\nzoom = \"1E8000\"\niterations = 5000\n",
        )
        .unwrap();
        assert_eq!(location.zoom, "1E8000");
        assert_eq!(location.zoom_exponent(), Some(8000));
        assert!(location.required_precision_bits(128).unwrap() > 26_000);
    }

    #[test]
    fn historical_palette_is_decoded_to_rgba_stops() {
        let (palette, span) =
            palette_from_toml("palette = [0, 0, 0, 255, 128, 64]\npalette_iteration_span = 50.0\n")
                .unwrap();
        assert_eq!(palette.stops().len(), 2);
        assert_eq!(span, 50.0);
    }
}
