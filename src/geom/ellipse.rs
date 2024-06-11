use std::borrow::Borrow;
use std::f64::consts::PI;

use num_complex::Complex64;

use crate::errors::ValueError;
use crate::geom::point::Point;
use crate::helpers::ApproximateEq;
use crate::helpers::IsReal;

#[derive(Clone, Debug)]
pub struct Ellipse {
    pub center: Point,
    pub rx: f64,
    pub ry: f64,
    pub phi: f64,
    pub clockwise: bool,
}

impl Ellipse {
    pub fn find_point<T: Borrow<Point>>(&self, point: T) -> f64 {
        // Solves the equation
        //   x = c + u exp(iθ) + v exp(-iθ)
        // for θ.
        //
        // u and v are defined as:
        //   u := exp(iφ) (rₓ + rᵧ) / 2
        //   v := exp(iφ) (rₓ - rᵧ) / 2
        //
        // θ has two solutions, if there is an imaginary part, the point doesn't lie on the ellipse.
        let point = point.borrow();

        let exp = Complex64::exp;
        let sqrt = Complex64::sqrt;
        let j = Complex64::I;
        let ln = Complex64::ln;

        let c: Complex64 = self.center.into();
        let x: Complex64 = point.into();
        let rx = self.rx;
        let ry = self.ry;
        let phi = self.phi;

        // Convert parametrization to
        // x = c + u exp(iθ) + v exp(-iθ), θ ∈ [θ₁, θ₁ + Δθ]
        let u = exp(j * phi) * (rx + ry) / 2.0;
        let v = exp(j * phi) * (rx - ry) / 2.0;

        // Calculate resulting angles
        let a = x - c;
        let sqrt_discriminant = sqrt(a.powi(2) - 4.0 * u * v);
        let i_plus = (a + sqrt_discriminant) / (2.0 * u);
        let i_minus = (a - sqrt_discriminant) / (2.0 * u);

        let theta_plus = -j * ln(i_plus);
        let theta_minus = -j * ln(i_minus);

        // Find the real theta
        let theta = if theta_plus.is_real() {
            theta_plus.re
        } else if theta_minus.is_real() {
            theta_minus.re
        } else {
            f64::NAN
        };

        theta % (2.0 * PI)
    }

    pub fn offset(&self, distance: f64) -> Result<Self, ValueError> {
        if !self.rx.approximate_eq(self.ry) {
            return Err(ValueError(
                "Not implemented for general ellipses, just for circle segments",
            ));
        }

        let sign = if self.clockwise { -1.0 } else { 1.0 };
        Ok(Self {
            center: self.center,
            rx: self.rx + sign * distance,
            ry: self.ry + sign * distance,
            clockwise: self.clockwise,
            phi: self.phi,
        })
    }
}
