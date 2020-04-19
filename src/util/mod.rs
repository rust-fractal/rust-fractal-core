pub mod complex_vector;
pub mod point;
pub mod image;

pub type ComplexFixed<T> = num_complex::Complex<T>;
pub type ComplexArbitrary = rug::Complex;

pub fn to_fixed(value: &ComplexArbitrary) -> ComplexFixed<f64> {
    ComplexFixed::new(value.real().to_f64(), value.imag().to_f64())
}
