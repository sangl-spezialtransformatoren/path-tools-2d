use std::borrow::Borrow;
use std::fmt::Display;
use std::fmt::Formatter;

use num_complex::Complex64;
use pyo3::{pyclass, pymethods, FromPyObject};

use crate::helpers::ApproximateEq;

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "python", pyclass)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Point { x, y }
    }
}

impl<T: Into<f64>, U: Into<f64>> From<(T, U)> for Point {
    fn from(value: (T, U)) -> Self {
        let x = value.0.into();
        let y = value.1.into();

        Point::new(x, y)
    }
}

impl From<Complex64> for Point {
    fn from(value: Complex64) -> Self {
        Self::new(value.re, value.im)
    }
}

impl From<Point> for Complex64 {
    fn from(value: Point) -> Self {
        Complex64::new(value.x, value.y)
    }
}

impl<'a> From<&'a Point> for Complex64 {
    fn from(value: &'a Point) -> Self {
        Complex64::new(value.x, value.y)
    }
}

impl Display for Point {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl ApproximateEq for Point {
    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool {
        let other = other.borrow();
        self.x.approximate_eq(other.x) && self.y.approximate_eq(other.y)
    }
}
