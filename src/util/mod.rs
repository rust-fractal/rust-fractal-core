use std::f64::consts::{LOG2_10, LOG10_2};

pub mod data_export;
pub mod float_extended;
pub mod complex_extended;

pub use complex_extended::ComplexExtended;
pub use float_extended::FloatExtended;

pub type ComplexFixed<T> = num_complex::Complex<T>;
pub type ComplexArbitrary = rug::Complex;

use std::os::raw::{c_double, c_int};

extern "C" {
    fn frexp(x: c_double, exp: *mut c_int) -> c_double;
    fn ldexp(x: c_double, exp: c_int) -> c_double;
}

pub trait FloatExp: Sized {
    fn frexp(self) -> (Self, i32);
    fn ldexp(self, exp: i32) -> Self;
}

impl FloatExp for f64 {
    fn frexp(self) -> (Self, i32) {
        let mut exp: c_int = 0;
        let res = unsafe { frexp(self, &mut exp) };
        (res, exp)
    }

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

pub fn string_to_extended(string: &String) -> FloatExtended {
    // Split on E as the exponent
    let temp: Vec<&str> = string.split('E').collect();

    let first = temp[0].parse::<f64>().unwrap();
    let second = temp[1].parse::<f64>().unwrap() * LOG2_10;

    FloatExtended::new(first * 2.0f64.powf(second.fract()), second.floor() as i32)
}

pub fn extended_to_string(value: FloatExtended) -> String {
    let first = value.mantissa;
    let second = value.exponent as f64 * LOG10_2;

    format!("{:.2}E{}", first * 10.0f64.powf(second.fract()), second.floor() as i32)
}

#[derive(Clone)]
pub struct PixelData {
    pub image_x: usize,
    pub image_y: usize,
    pub iteration: usize,
    pub delta_centre: ComplexExtended,
    pub delta_reference: ComplexExtended,
    pub delta_current: ComplexExtended,
    pub derivative_current: ComplexFixed<f64>,
    pub glitched: bool,
    pub escaped: bool,
}