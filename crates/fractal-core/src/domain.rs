use core::fmt;

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexCoordinate {
    pub real: f64,
    pub imag: f64,
    /// Original decimal values, when this coordinate came from an exact
    /// location. The f64 fields remain useful as a cheap shallow hint; deep
    /// rendering never reconstructs the reference from those hints.
    #[cfg_attr(feature = "serde", serde(default, skip_serializing_if = "Option::is_none"))]
    pub real_decimal: Option<Box<str>>,
    #[cfg_attr(feature = "serde", serde(default, skip_serializing_if = "Option::is_none"))]
    pub imag_decimal: Option<Box<str>>,
}

impl ComplexCoordinate {
    pub const fn new(real: f64, imag: f64) -> Self {
        Self { real, imag, real_decimal: None, imag_decimal: None }
    }

    /// Creates a coordinate without converting its decimal components through
    /// f64. Conversion is retained only as a shallow-path hint.
    pub fn from_decimal(real: &str, imag: &str) -> Result<Self, RenderError> {
        let real = real.trim();
        let imag = imag.trim();
        let real_hint =
            decimal_f64_hint(real).ok_or(RenderError::InvalidDecimalCoordinate { axis: "real" })?;
        let imag_hint =
            decimal_f64_hint(imag).ok_or(RenderError::InvalidDecimalCoordinate { axis: "imag" })?;
        Ok(Self {
            real: real_hint,
            imag: imag_hint,
            real_decimal: Some(real.into()),
            imag_decimal: Some(imag.into()),
        })
    }

    pub fn decimal_components(&self) -> Option<(&str, &str)> {
        Some((self.real_decimal.as_deref()?, self.imag_decimal.as_deref()?))
    }

    pub fn has_exact_components(&self) -> bool {
        self.decimal_components().is_some()
    }

    pub(crate) fn is_finite(&self) -> bool {
        self.real.is_finite() && self.imag.is_finite()
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Zoom {
    scale: f64,
    #[cfg_attr(feature = "serde", serde(default, skip_serializing_if = "Option::is_none"))]
    decimal: Option<Box<str>>,
}

impl Zoom {
    pub fn new(scale: f64) -> Result<Self, RenderError> {
        if scale.is_finite() && scale > 0.0 {
            Ok(Self { scale, decimal: None })
        } else {
            Err(RenderError::InvalidZoom { scale })
        }
    }

    /// Creates a zoom from an exact decimal/scientific value. Values such as
    /// `1E8000` are intentionally not parsed as f64.
    pub fn from_decimal(value: &str) -> Result<Self, RenderError> {
        let value = value.trim();
        let scale = decimal_f64_hint(value)
            .filter(|scale| *scale > 0.0)
            .unwrap_or_else(|| if decimal_is_positive(value) { f64::MAX } else { f64::NAN });
        if !decimal_is_positive(value) {
            return Err(RenderError::InvalidDecimalZoom);
        }
        Ok(Self { scale, decimal: Some(value.into()) })
    }

    pub const fn unit() -> Self {
        Self { scale: 1.0, decimal: None }
    }

    pub fn scale(&self) -> f64 {
        self.scale
    }

    pub fn decimal(&self) -> Option<&str> {
        self.decimal.as_deref()
    }

    pub fn has_exact_scale(&self) -> bool {
        self.decimal.is_some()
    }

    pub(crate) fn binary_scale_bits(&self) -> Result<u32, RenderError> {
        let logarithm = if let Some(decimal) = self.decimal() {
            decimal_log2(decimal).ok_or(RenderError::InvalidDecimalZoom)?
        } else {
            self.scale.log2()
        };
        if !logarithm.is_finite() || logarithm > f64::from(u32::MAX) {
            return Err(RenderError::InvalidZoom { scale: self.scale });
        }
        Ok(logarithm.ceil().max(0.0) as u32)
    }
}

impl Default for Zoom {
    fn default() -> Self {
        Self::unit()
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct Viewport {
    pub center: ComplexCoordinate,
    pub zoom: Zoom,
    pub rotation_radians: f64,
    pub width: u32,
    pub height: u32,
}

impl Viewport {
    pub fn validate(&self) -> Result<(), RenderError> {
        if !self.center.is_finite() && !self.center.has_exact_components() {
            return Err(RenderError::NonFiniteCoordinate {
                real: self.center.real,
                imag: self.center.imag,
            });
        }
        if !self.rotation_radians.is_finite() {
            return Err(RenderError::NonFiniteRotation { rotation_radians: self.rotation_radians });
        }
        if self.width == 0 || self.height == 0 {
            return Err(RenderError::InvalidDimensions { width: self.width, height: self.height });
        }
        Ok(())
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FractalFormula {
    Mandelbrot { power: u32 },
    BurningShip { power: u32 },
}

impl Default for FractalFormula {
    fn default() -> Self {
        Self::Mandelbrot { power: 2 }
    }
}

impl FractalFormula {
    pub fn validate(self) -> Result<(), RenderError> {
        let power = match self {
            Self::Mandelbrot { power } | Self::BurningShip { power } => power,
        };
        if (2..=16).contains(&power) { Ok(()) } else { Err(RenderError::InvalidPower { power }) }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ApproximationSettings {
    pub enabled: bool,
    pub order: u32,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlitchSettings {
    pub enabled: bool,
    pub tolerance: f64,
}

impl Default for GlitchSettings {
    fn default() -> Self {
        Self { enabled: true, tolerance: 1.4e-6 }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SamplingSettings {
    pub samples_per_pixel: u32,
}

impl Default for SamplingSettings {
    fn default() -> Self {
        Self { samples_per_pixel: 1 }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OutputChannels {
    pub iterations: bool,
    pub smooth_iteration: bool,
    pub distance: bool,
    pub stripes: bool,
    pub validity: bool,
}

impl Default for OutputChannels {
    fn default() -> Self {
        Self {
            iterations: true,
            smooth_iteration: false,
            distance: false,
            stripes: false,
            validity: false,
        }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PrecisionPolicy {
    pub minimum_bits: u32,
    pub guard_bits: u32,
    pub max_bits: Option<u32>,
}

impl Default for PrecisionPolicy {
    fn default() -> Self {
        Self { minimum_bits: 64, guard_bits: 32, max_bits: None }
    }
}

impl PrecisionPolicy {
    pub fn required_bits(self, zoom: Zoom) -> Result<u32, RenderError> {
        if self.minimum_bits == 0 || self.guard_bits == 0 {
            return Err(RenderError::InvalidPrecisionPolicy);
        }

        if let Some(max_bits) = self.max_bits
            && max_bits < self.minimum_bits
        {
            return Err(RenderError::InvalidPrecisionPolicy);
        }

        let zoom_bits = zoom.binary_scale_bits()?;
        let required_bits = self
            .minimum_bits
            .checked_add(self.guard_bits)
            .and_then(|bits| bits.checked_add(zoom_bits))
            .ok_or(RenderError::InvalidPrecisionPolicy)?;

        if let Some(max_bits) = self.max_bits
            && required_bits > max_bits
        {
            return Err(RenderError::PrecisionLimitExceeded { required_bits, max_bits });
        }

        Ok(required_bits)
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderSettings {
    pub fractal: FractalFormula,
    pub iteration_limit: u32,
    pub approximation: ApproximationSettings,
    pub glitch: GlitchSettings,
    pub sampling: SamplingSettings,
    pub output: OutputChannels,
    pub precision: PrecisionPolicy,
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            fractal: FractalFormula::default(),
            iteration_limit: 1_000,
            approximation: ApproximationSettings::default(),
            glitch: GlitchSettings::default(),
            sampling: SamplingSettings::default(),
            output: OutputChannels::default(),
            precision: PrecisionPolicy::default(),
        }
    }
}

impl RenderSettings {
    pub fn validate(self) -> Result<(), RenderError> {
        self.fractal.validate()?;
        if self.iteration_limit == 0 {
            return Err(RenderError::InvalidIterationLimit);
        }
        if self.approximation.order > 64
            || (self.approximation.enabled && self.approximation.order == 0)
        {
            return Err(RenderError::InvalidApproximationOrder { order: self.approximation.order });
        }
        if !self.glitch.tolerance.is_finite() || self.glitch.tolerance <= 0.0 {
            return Err(RenderError::InvalidGlitchTolerance { tolerance: self.glitch.tolerance });
        }
        if self.sampling.samples_per_pixel == 0 {
            return Err(RenderError::InvalidSampling {
                samples_per_pixel: self.sampling.samples_per_pixel,
            });
        }
        self.precision.required_bits(Zoom::unit())?;
        Ok(())
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct RenderRequest {
    pub generation: u64,
    pub viewport: Viewport,
    pub settings: RenderSettings,
}

impl RenderRequest {
    pub fn validate(&self) -> Result<(), RenderError> {
        self.viewport.validate()?;
        self.settings.validate()?;
        self.settings.precision.required_bits(self.viewport.zoom.clone())?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RenderError {
    InvalidDimensions { width: u32, height: u32 },
    InvalidZoom { scale: f64 },
    InvalidPower { power: u32 },
    InvalidIterationLimit,
    InvalidApproximationOrder { order: u32 },
    InvalidGlitchTolerance { tolerance: f64 },
    InvalidSampling { samples_per_pixel: u32 },
    InvalidPrecisionPolicy,
    PrecisionLimitExceeded { required_bits: u32, max_bits: u32 },
    InvalidColourSettings,
    InvalidRawFrame,
    Cancelled,
    NonFiniteCoordinate { real: f64, imag: f64 },
    InvalidDecimalCoordinate { axis: &'static str },
    InvalidDecimalZoom,
    ReferenceGeneration(String),
    NumericOverflow,
    NonFiniteRotation { rotation_radians: f64 },
    InvalidTile,
}

impl fmt::Display for RenderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimensions { width, height } => {
                write!(f, "render dimensions must be non-zero: {width}x{height}")
            }
            Self::InvalidZoom { scale } => write!(f, "zoom must be finite and positive: {scale}"),
            Self::InvalidPower { power } => {
                write!(f, "fractal power must be between 2 and 16: {power}")
            }
            Self::InvalidIterationLimit => write!(f, "iteration limit must be non-zero"),
            Self::InvalidApproximationOrder { order } => {
                write!(f, "approximation order is too large: {order}")
            }
            Self::InvalidGlitchTolerance { tolerance } => {
                write!(f, "glitch tolerance must be finite and positive: {tolerance}")
            }
            Self::InvalidSampling { samples_per_pixel } => {
                write!(f, "samples per pixel must be non-zero: {samples_per_pixel}")
            }
            Self::InvalidPrecisionPolicy => write!(f, "precision policy is inconsistent"),
            Self::PrecisionLimitExceeded { required_bits, max_bits } => write!(
                f,
                "requested zoom requires {required_bits} precision bits, exceeding the {max_bits}-bit limit"
            ),
            Self::InvalidColourSettings => write!(f, "colour settings are invalid"),
            Self::InvalidRawFrame => write!(f, "raw frame channels do not match its dimensions"),
            Self::Cancelled => write!(f, "render was cancelled"),
            Self::NonFiniteCoordinate { real, imag } => {
                write!(f, "viewport center is not finite: ({real}, {imag})")
            }
            Self::InvalidDecimalCoordinate { axis } => {
                write!(f, "decimal coordinate component is invalid: {axis}")
            }
            Self::InvalidDecimalZoom => f.write_str("decimal zoom is invalid or non-positive"),
            Self::ReferenceGeneration(error) => write!(f, "reference generation failed: {error}"),
            Self::NumericOverflow => f.write_str("extended numerical value overflowed"),
            Self::NonFiniteRotation { rotation_radians } => {
                write!(f, "viewport rotation is not finite: {rotation_radians}")
            }
            Self::InvalidTile => write!(f, "tile lies outside the viewport"),
        }
    }
}

fn decimal_is_positive(value: &str) -> bool {
    decimal_log2(value).is_some_and(|logarithm| logarithm.is_finite())
}

fn decimal_f64_hint(value: &str) -> Option<f64> {
    let value = value.trim();
    if value.is_empty()
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'+' | b'-' | b'.' | b'e' | b'E'))
        || !value.bytes().any(|byte| byte.is_ascii_digit())
    {
        return None;
    }
    let parsed = value.parse::<f64>().ok()?;
    if parsed.is_nan() {
        return None;
    }
    if parsed.is_finite() {
        Some(parsed)
    } else if parsed.is_sign_positive() {
        Some(f64::MAX)
    } else {
        Some(-f64::MAX)
    }
}

fn decimal_log2(value: &str) -> Option<f64> {
    let value = value.trim();
    if value.is_empty() {
        return None;
    }
    let (significand, exponent) =
        value.split_once(['e', 'E']).map_or((value, 0_i64), |(mantissa, exponent)| {
            (mantissa, exponent.parse::<i64>().unwrap_or(i64::MIN))
        });
    if exponent == i64::MIN {
        return None;
    }
    let mantissa = significand.parse::<f64>().ok()?;
    if !mantissa.is_finite() || mantissa == 0.0 {
        return if mantissa == 0.0 { None } else { Some(f64::INFINITY) };
    }
    let magnitude = mantissa.abs().log2() + (exponent as f64) * 10_f64.log2();
    (mantissa.is_sign_positive() && magnitude.is_finite()).then_some(magnitude)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_requests_with_typed_errors() {
        let request = RenderRequest {
            generation: 1,
            viewport: Viewport {
                center: ComplexCoordinate::new(0.0, 0.0),
                zoom: Zoom::new(0.0).unwrap_or_default(),
                rotation_radians: 0.0,
                width: 0,
                height: 100,
            },
            settings: RenderSettings::default(),
        };

        assert_eq!(
            request.validate(),
            Err(RenderError::InvalidDimensions { width: 0, height: 100 })
        );
    }

    #[test]
    fn precision_maximum_cannot_be_below_minimum() {
        let mut settings = RenderSettings::default();
        settings.precision.max_bits = Some(32);
        assert_eq!(settings.validate(), Err(RenderError::InvalidPrecisionPolicy));
    }

    #[test]
    fn required_precision_grows_monotonically_with_zoom() {
        let policy = PrecisionPolicy::default();
        let low = policy.required_bits(Zoom::new(1.0).unwrap()).unwrap();
        let high = policy.required_bits(Zoom::new(1_000_000.0).unwrap()).unwrap();
        assert!(high > low);
    }

    #[test]
    fn precision_limit_is_reported_as_a_typed_error() {
        let policy = PrecisionPolicy { max_bits: Some(100), ..PrecisionPolicy::default() };
        assert_eq!(
            policy.required_bits(Zoom::new(1_000_000.0).unwrap()),
            Err(RenderError::PrecisionLimitExceeded { required_bits: 116, max_bits: 100 })
        );
    }

    #[test]
    fn decimal_inputs_reject_special_float_strings_but_keep_deep_values() {
        assert_eq!(
            ComplexCoordinate::from_decimal("NaN", "0"),
            Err(RenderError::InvalidDecimalCoordinate { axis: "real" })
        );
        assert_eq!(Zoom::from_decimal("inf"), Err(RenderError::InvalidDecimalZoom));
        assert!(Zoom::from_decimal("1E8000").unwrap().has_exact_scale());
    }
}
