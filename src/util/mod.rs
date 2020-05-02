use float_extended::complex_extended::ComplexExtended;
use rug::float::Round;

pub mod point;
pub mod image;

pub type ComplexFixed<T> = num_complex::Complex<T>;
pub type ComplexArbitrary = rug::Complex;

pub fn to_fixed(value: &ComplexArbitrary) -> ComplexFixed<f64> {
    let re = value.real().to_f64();
    let im = value.imag().to_f64();

    ComplexFixed::new(re, im)
}

pub fn to_fixed_exp(value: &ComplexArbitrary) -> (ComplexFixed<f64>, i32) {
    let (re, p1) = value.real().to_f64_exp();
    let (im, p2) = value.imag().to_f64_exp();

    let im_new = im * 2f64.powi(p2 - p1);

    (ComplexFixed::new(re, im_new), p1)
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
    pub p_current: i32,
    pub delta_reference: ComplexFixed<f64>,
    pub delta_current: ComplexFixed<f64>,
    pub derivative_current: ComplexFixed<f64>,
    pub glitched: bool,
    pub escaped: bool,
}