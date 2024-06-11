use std::borrow::Borrow;

use num_complex::Complex;

pub trait ApproximateEq {
    const EPSILON: f64 = 1e-2;

    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool;
}

impl ApproximateEq for f64 {
    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool {
        let other = other.borrow();
        (self - other).abs() < <f64 as ApproximateEq>::EPSILON
    }
}

impl ApproximateEq for Complex<f64> {
    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool {
        let other = other.borrow();
        self.re.approximate_eq(other.re) && self.im.approximate_eq(other.im)
    }
}
