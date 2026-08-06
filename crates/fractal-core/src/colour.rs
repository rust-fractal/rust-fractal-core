use crate::{RawFrame, RenderError};
use core::{fmt, slice};

#[derive(Clone, Debug, PartialEq)]
pub struct Palette {
    stops: Box<[[f32; 4]]>,
    cyclic: bool,
}

impl Palette {
    pub fn new(stops: Vec<[f32; 4]>, cyclic: bool) -> Result<Self, ColourError> {
        if stops.is_empty() || stops.iter().flatten().any(|value| !value.is_finite()) {
            return Err(ColourError::InvalidPalette);
        }
        Ok(Self { stops: stops.into_boxed_slice(), cyclic })
    }

    pub fn grayscale() -> Self {
        Self::new(vec![[0.0, 0.0, 0.0, 1.0], [1.0, 1.0, 1.0, 1.0]], false)
            .expect("built-in grayscale palette is valid")
    }

    pub fn stops(&self) -> &[[f32; 4]] {
        &self.stops
    }

    pub fn cyclic(&self) -> bool {
        self.cyclic
    }

    fn sample(&self, value: f64) -> [f32; 4] {
        if self.stops.len() == 1 {
            return self.stops[0];
        }
        let value = if self.cyclic { value.rem_euclid(1.0) } else { value.clamp(0.0, 1.0) };
        let (position, upper_index) = if self.cyclic {
            (value * self.stops.len() as f64, self.stops.len())
        } else {
            (value * (self.stops.len() - 1) as f64, self.stops.len() - 1)
        };
        let lower = position.floor() as usize;
        let upper =
            if self.cyclic { (lower + 1) % upper_index } else { (lower + 1).min(upper_index) };
        let fraction = (position - lower as f64) as f32;
        let mut result = [0.0; 4];
        for (channel, output) in result.iter_mut().enumerate() {
            *output = self.stops[lower][channel]
                + (self.stops[upper][channel] - self.stops[lower][channel]) * fraction;
        }
        result
    }
}

impl Default for Palette {
    fn default() -> Self {
        Self::grayscale()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ColourSettings {
    pub scale: f64,
    pub offset: f64,
    pub distance_strength: f32,
    pub interior: [u8; 4],
    pub invalid: [u8; 4],
}

impl Default for ColourSettings {
    fn default() -> Self {
        Self {
            scale: 1.0,
            offset: 0.0,
            distance_strength: 0.0,
            interior: [0, 0, 0, 255],
            invalid: [255, 0, 255, 255],
        }
    }
}

impl ColourSettings {
    fn validate(self) -> Result<(), ColourError> {
        if !self.scale.is_finite() || !self.offset.is_finite() || self.scale == 0.0 {
            return Err(ColourError::InvalidSettings);
        }
        if !self.distance_strength.is_finite() || self.distance_strength < 0.0 {
            return Err(ColourError::InvalidSettings);
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImageFrame {
    pub generation: u64,
    pub width: u32,
    pub height: u32,
    pub rgba8: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ColourError {
    InvalidPalette,
    InvalidSettings,
    InvalidRawFrame,
}

impl fmt::Display for ColourError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidPalette => "palette must contain finite colour stops",
            Self::InvalidSettings => "colour settings are invalid",
            Self::InvalidRawFrame => "raw frame channels do not match its dimensions",
        })
    }
}

impl From<ColourError> for RenderError {
    fn from(error: ColourError) -> Self {
        match error {
            ColourError::InvalidPalette | ColourError::InvalidSettings => {
                RenderError::InvalidColourSettings
            }
            ColourError::InvalidRawFrame => RenderError::InvalidRawFrame,
        }
    }
}

pub struct ColourPipeline {
    palette: Palette,
    settings: ColourSettings,
}

impl ColourPipeline {
    pub fn new(palette: Palette, settings: ColourSettings) -> Result<Self, ColourError> {
        settings.validate()?;
        Ok(Self { palette, settings })
    }

    pub fn palette(&self) -> &Palette {
        &self.palette
    }

    pub fn settings(&self) -> ColourSettings {
        self.settings
    }

    pub fn map(&self, frame: &RawFrame) -> Result<ImageFrame, ColourError> {
        frame.validate().map_err(|_| ColourError::InvalidRawFrame)?;
        let pixel_count = (frame.width as usize)
            .checked_mul(frame.height as usize)
            .ok_or(ColourError::InvalidRawFrame)?;
        let mut rgba8 = Vec::with_capacity(pixel_count * 4);
        let smooth = frame.smooth_iteration.as_deref();
        let distance = frame.distance.as_deref();
        let validity = frame.validity.as_deref();
        for index in 0..pixel_count {
            let pixel = if validity.is_some_and(|channel| channel[index] == 0) {
                self.settings.invalid
            } else if frame.iterations[index] >= frame.iteration_limit {
                self.settings.interior
            } else {
                let value = smooth.map_or(frame.iterations[index] as f64, |channel| channel[index]);
                let mut colour =
                    self.palette.sample(value * self.settings.scale + self.settings.offset);
                if let Some(channel) = distance {
                    let light = (1.0 - channel[index] as f32 * self.settings.distance_strength)
                        .clamp(0.0, 1.0);
                    for component in &mut colour[..3] {
                        *component *= light;
                    }
                }
                colour.map(|component| (component.clamp(0.0, 1.0) * 255.0).round() as u8)
            };
            rgba8.extend_from_slice(&pixel);
        }
        Ok(ImageFrame {
            generation: frame.generation,
            width: frame.width,
            height: frame.height,
            rgba8,
        })
    }
}

impl<'a> IntoIterator for &'a Palette {
    type Item = &'a [f32; 4];
    type IntoIter = slice::Iter<'a, [f32; 4]>;

    fn into_iter(self) -> Self::IntoIter {
        self.stops.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ComplexCoordinate, FractalFormula, OutputChannels, RenderSettings, Viewport, Zoom,
        compute_frame,
    };

    fn frame() -> RawFrame {
        let request = crate::RenderRequest {
            generation: 4,
            viewport: Viewport {
                center: ComplexCoordinate::new(0.0, 0.0),
                zoom: Zoom::unit(),
                rotation_radians: 0.0,
                width: 4,
                height: 4,
            },
            settings: RenderSettings {
                fractal: FractalFormula::default(),
                iteration_limit: 8,
                output: OutputChannels {
                    smooth_iteration: true,
                    distance: true,
                    ..OutputChannels::default()
                },
                ..RenderSettings::default()
            },
        };
        compute_frame(&request).unwrap()
    }

    #[test]
    fn palette_wraps_only_when_requested() {
        let cyclic = Palette::new(vec![[0.0, 0.0, 0.0, 1.0], [1.0, 1.0, 1.0, 1.0]], true).unwrap();
        let pipeline = ColourPipeline::new(cyclic, ColourSettings::default()).unwrap();
        let image = pipeline.map(&frame()).unwrap();
        assert_eq!(image.rgba8.len(), 64);
    }

    #[test]
    fn cyclic_palette_interpolates_from_last_stop_to_first() {
        let palette = Palette::new(vec![[1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0]], true).unwrap();
        let pipeline = ColourPipeline::new(
            palette,
            ColourSettings { scale: 0.75, ..ColourSettings::default() },
        )
        .unwrap();
        let image = pipeline.map(&frame()).unwrap();
        assert!(image.rgba8.chunks_exact(4).any(|pixel| pixel[0] > 0 && pixel[2] > 0));
    }

    #[test]
    fn interior_pixels_are_black_by_default() {
        let pipeline = ColourPipeline::new(Palette::default(), ColourSettings::default()).unwrap();
        let image = pipeline.map(&frame()).unwrap();
        assert!(image.rgba8.chunks_exact(4).any(|pixel| pixel == [0, 0, 0, 255]));
    }
}
