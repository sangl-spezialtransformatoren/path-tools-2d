use std::borrow::Borrow;
use std::fmt::Display;
use std::fmt::Formatter;

use num_complex::Complex64;
use num_complex::ComplexFloat;

use crate::geom::Point;
use crate::helpers::IsReal;

#[derive(Copy, Clone, Debug)]
pub struct Straight {
    pub a: Point,
    pub b: Point,
}

impl Display for Straight {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} -- {}", self.a, self.b)
    }
}

impl Straight {
    pub fn find_point<T: Borrow<Point>>(&self, point: T) -> f64 {
        let point = point.borrow();

        // Solves the equation x = a + λ (b - a) for λ.
        let a: Complex64 = self.a.into();
        let b: Complex64 = self.b.into();
        let x: Complex64 = point.into();

        // Calculate λ
        let lambda = (x - a) / (b - a);

        // Check if λ has an imaginary part
        if !lambda.is_real() {
            return std::f64::NAN;
        }

        // Return the real part of λ
        lambda.re
    }

    pub fn offset(&self, distance: f64) -> Self {
        let j = Complex64::I;
        let a: Complex64 = self.a.into();
        let b: Complex64 = self.b.into();

        let normal = j * (b - a) / (b - a).abs();
        Self {
            a: Point::from(a + distance * normal),
            b: Point::from(b + distance * normal),
        }
    }
}
