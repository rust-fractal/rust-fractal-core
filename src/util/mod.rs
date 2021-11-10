use std::f64::consts::{LOG2_10, LOG10_2};

pub mod data_export;
pub mod float_extended;
pub mod complex_extended;
pub mod recolour_exr;
pub mod progress;

use colorgrad::{CustomGradient, Interpolation, Color, BlendMode};
pub use complex_extended::ComplexExtended;
use config::Config;
pub use float_extended::FloatExtended;
pub use recolour_exr::RecolourExr;
pub use progress::ProgressCounters;

pub type ComplexFixed<T> = num_complex::Complex<T>;
pub type ComplexArbitrary = rug::Complex;
pub type FloatArbitrary = rug::Float;

use std::os::raw::{c_double, c_int};

use self::data_export::{ColoringType, DataType};

extern "C" {
    fn frexp(x: c_double, exp: *mut c_int) -> c_double;
    fn ldexp(x: c_double, exp: c_int) -> c_double;
}

pub trait FloatExp: Sized {
    fn frexp(self) -> (Self, i32);
    fn ldexp(self, exp: i32) -> Self;
}

impl FloatExp for f64 {
    #[inline]
    fn frexp(self) -> (Self, i32) {
        let mut exp: c_int = 0;
        let res = unsafe { frexp(self, &mut exp) };
        (res, exp)
    }

    #[inline]
    fn ldexp(self, exp: i32) -> Self {
        unsafe { ldexp(self, exp) }
    }
}

#[inline]
pub fn to_fixed(value: &ComplexArbitrary) -> ComplexFixed<f64> {
    let re = value.real().to_f64();
    let im = value.imag().to_f64();

    ComplexFixed::new(re, im)
}

pub fn to_extended(value: &ComplexArbitrary) -> ComplexExtended {
    let (mut re, p1) = value.real().to_f64_exp();
    let (mut im, p2) = value.imag().to_f64_exp();

    let exponent = if p1 < p2 {
        re *= 2f64.powi(p1 - p2);
        p2
    } else {
        im *= 2f64.powi(p2 - p1);
        p1
    };

    ComplexExtended::new2(re, im, exponent)
}

pub fn string_to_extended(string: &str) -> FloatExtended {
    let temp: Vec<&str> = string.split(&['E', 'e'][..]).collect();

    let first = temp[0].parse::<f64>().unwrap();

    let second = if temp[1].is_empty() {
        0.0
    } else {
        temp[1].parse::<f64>().unwrap() * LOG2_10
    };

    if second < 0.0 {
        FloatExtended::new(first * 2.0f64.powf(1.0 + second.fract()), second.floor() as i32)
    } else {
        FloatExtended::new(first * 2.0f64.powf(second.fract()), second.floor() as i32)
    }
}

pub fn extended_to_string_short(value: FloatExtended) -> String {
    let first = value.mantissa;
    let second = value.exponent as f64 * LOG10_2;

    if second < 0.0 {
        format!("{:.2}E{}", first * 10.0f64.powf(1.0 + second.fract()), second.floor() as i32)
    } else {
        format!("{:.2}E{}", first * 10.0f64.powf(second.fract()), second.floor() as i32)
    }
}

pub fn extended_to_string_long(value: FloatExtended) -> String {
    let first = value.mantissa;
    let second = value.exponent as f64 * LOG10_2;

    if second < 0.0 {
        format!("{}E{}", first * 10.0f64.powf(1.0 + second.fract()), second.floor() as i32)
    } else {
        format!("{}E{}", first * 10.0f64.powf(second.fract()), second.floor() as i32)
    }
}

pub fn linear_interpolation_between_zoom(zoom1: FloatExtended, zoom2: FloatExtended, factor: f64) -> FloatExtended {
    let temp = (1.0 - factor) * (zoom1.mantissa.log2() + zoom1.exponent as f64) + factor * (zoom2.mantissa.log2() + zoom2.exponent as f64);

    FloatExtended::new(2.0f64.powf(temp.fract()), temp.floor() as i32)
}

pub fn generate_pascal_coefficients(row: usize) -> Vec<f64> {
    let mut current_row = vec![1.0; row + 1];

    for i in 2..row {
        let previous_row = current_row.clone();

        for j in 1..i {
            current_row[j] = previous_row[j - 1] + previous_row[j]
        }
    }

    current_row
}

#[inline]
pub fn diff_abs(a: f64, b: f64) -> f64 {
    match (a >= 0.0, a + b >= 0.0) {
        (true, true) => b,
        (true, _) => -2.0 * a - b,
        (_, true) => 2.0 * a + b,
        (_, _) => -b
    }
}

pub fn get_delta_top_left(delta_pixel: f64, image_width: usize, image_height: usize, cos_rotate: f64, sin_rotate: f64) -> ComplexFixed<f64> {
    let aspect = image_width as f64 / image_height as f64;

    let temp_real = -0.5 * (image_height - 1) as f64 * delta_pixel * aspect as f64;
    let temp_imag = -0.5 * (image_height - 1) as f64 * delta_pixel;

    ComplexFixed::new(
        temp_real * cos_rotate - temp_imag * sin_rotate, 
        temp_real * sin_rotate + temp_imag * cos_rotate)
}

pub fn get_approximation_terms(approximation_order: usize, image_width: usize, image_height: usize) -> usize {
    if approximation_order == 0 {
        let auto = (((image_width * image_height) as f64).log(1e6).powf(6.619) * 16.0f64) as usize;
        auto.clamp(3, 64)
    } else {
        approximation_order
    }
}

pub fn generate_default_palette() -> (Vec<Color>, Vec<Color>) {
    let palette_generator = CustomGradient::new()
        .colors(
            &[Color::from_rgb_u8(0, 2, 0),
            Color::from_rgb_u8(0, 7, 100), 
            Color::from_rgb_u8(32, 107, 203), 
            Color::from_rgb_u8(237, 255, 255),
            Color::from_rgb_u8(255, 170, 0),
            Color::from_rgb_u8(0, 2, 0)])
            .domain(&[0.0, 0.1425, 0.3025, 0.5625, 0.785, 1.0])
            .interpolation(Interpolation::CatmullRom).mode(BlendMode::Oklab)
            .build().unwrap();

    (palette_generator.colors(6), palette_generator.colors(6 * 64))
}

pub fn get_fractal_type_from_settings(settings: &Config) -> FractalType {
    let fractal_power = settings.get_int("fractal_power").unwrap_or(2) as usize;

    match settings.get("fractal_type").unwrap_or_else(|_| String::from("mandelbrot")).to_ascii_uppercase().as_ref() {
        "MANDELBROT" => FractalType::Mandelbrot(fractal_power),
        "BURNINGSHIP" => FractalType::BurningShip(fractal_power),
        _ => FractalType::Mandelbrot(fractal_power)
    }
}

pub fn get_data_coloring_type_from_settings(settings: &Config) -> (ColoringType, DataType) {
    let coloring_type = match settings.get("coloring_type").unwrap_or_else(|_| String::from("smooth_iteration")).to_ascii_uppercase().as_ref() {
        "SMOOTH_ITERATION" | "SMOOTH" => ColoringType::SmoothIteration,
        "STEP_ITERATION" | "STEP" => ColoringType::StepIteration,
        "DISTANCE" => ColoringType::Distance,
        "STRIPE" => ColoringType::Stripe,
        "DISTANCE_STRIPE" => ColoringType::DistanceStripe,
        _ => ColoringType::SmoothIteration
    };

    let data_type = match coloring_type {
        ColoringType::SmoothIteration | ColoringType::StepIteration => DataType::Iteration,
        ColoringType::Stripe => DataType::Stripe,
        ColoringType::DistanceStripe => DataType::DistanceStripe,
        _ => DataType::Distance
    };

    (coloring_type, data_type)
}

#[derive(Clone)]
pub struct PixelData {
    pub index: usize,
    pub iteration: usize,
    pub reference_iteration: usize,
    pub delta_reference: ComplexExtended,
    pub delta_current: ComplexExtended,
    pub jacobian_current: [ComplexExtended; 2],
    pub z_norm: f64,
    pub stripe_storage: [ComplexFixed<f64>; 4],
    pub stripe_iteration: usize,
}

#[derive(Copy, Clone, PartialEq)]
pub enum FractalType {
    Mandelbrot(usize),
    BurningShip(usize)
}