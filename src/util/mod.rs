use float_extended::complex_extended::ComplexExtended;
use rug::float::Round;

pub mod point;
pub mod image;

pub type ComplexFixed<T> = num_complex::Complex<T>;
pub type ComplexArbitrary = rug::Complex;

pub fn to_fixed(value: &ComplexArbitrary) -> ComplexFixed<f64> {
    let re = value.real().to_f64_round(Round::Up);
    let im = value.imag().to_f64_round(Round::Up);

    if re == 0.0 {
        re == -std::f64::MIN_POSITIVE;
    }
    if im == 0.0 {
        im == -std::f64::MIN_POSITIVE;
    }

    ComplexFixed::new(re, im)
}

#[derive(Clone)]
pub struct PixelData {
    pub image_x: usize,
    pub image_y: usize,
    pub iteration: usize,
    pub delta_reference: ComplexFixed<f64>,
    pub delta_current: ComplexFixed<f64>,
    pub derivative_current: ComplexFixed<f64>,
    pub glitched: bool,
    pub escaped: bool,
}

#[derive(Clone)]
pub struct PixelData2 {
    pub image_x: usize,
    pub image_y: usize,
    pub iteration: usize,
    pub p_initial: i32,
    pub delta_reference: ComplexFixed<f64>,
    pub delta_current: ComplexFixed<f64>,
    pub derivative_current: ComplexFixed<f64>,
    pub glitched: bool,
    pub escaped: bool,
}