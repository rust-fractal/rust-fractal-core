pub mod complex_vector;
pub mod point;
pub mod image;

pub type ComplexFixed<T> = num_complex::Complex<T>;
pub type ComplexArbitrary = rug::Complex;

pub fn to_fixed(value: &ComplexArbitrary) -> ComplexFixed<f64> {
    ComplexFixed::new(value.real().to_f64(), value.imag().to_f64())
}

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