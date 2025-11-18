use core::f64;

use num_complex::{Complex, ComplexFloat};
use num_traits::{Num, Float};

pub trait ToComplex<T> {
    fn to_complex(&self) -> Complex<T>;
}


impl<T: Num + Copy> ToComplex<T> for T {
    fn to_complex(&self) -> Complex<T> {
        Complex::new(*self, T::zero())
    }
}


impl<T: Num + Copy> ToComplex<T> for Complex<T> {
    fn to_complex(&self) -> Complex<T> {
        *self
    }
}