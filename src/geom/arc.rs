use std::borrow::Borrow;
use std::f64::consts::PI;
use std::fmt::Display;
use std::fmt::Formatter;

use derive_more::Constructor;
use num_complex::Complex64;
use num_complex::ComplexFloat;
use spec_math::cephes64::ellie;

use crate::algorithms::IntersectionSolution;
use crate::geom::Ellipse;
use crate::geom::Point;
use crate::geom::Segment;
use crate::geom::SegmentTrait;
use crate::helpers::ApproximateEq;

#[derive(Copy, Clone, Debug)]
pub struct Arc {
    start: Point,
    end: Point,
    rx: f64,
    ry: f64,
    phi: f64,
    large_arc: bool,
    sweep: bool,
}

#[derive(Constructor, Copy, Clone)]
pub struct CenterParameterizedArc {
    center: Point,
    rx: f64,
    ry: f64,
    phi: f64,
    theta1: f64,
    d_theta: f64,
}

impl Arc {
    pub fn new<T: Into<Point>, U: Into<Point>>(
        start: T,
        end: U,
        rx: f64,
        ry: f64,
        phi: f64,
        large_arc: bool,
        sweep: bool,
    ) -> Self {
        let start = start.into();
        let end = end.into();

        Arc {
            start,
            end,
            rx,
            ry,
            phi,
            large_arc,
            sweep,
        }
    }

    pub fn rx(&self) -> f64 {
        self.rx
    }

    pub fn ry(&self) -> f64 {
        self.ry
    }

    pub fn phi(&self) -> f64 {
        self.phi
    }

    pub fn large_arc(&self) -> bool {
        self.large_arc
    }

    pub fn sweep(&self) -> bool {
        self.sweep
    }

    pub fn center(&self) -> Point {
        let center_params = CenterParameterizedArc::from(self);
        center_params.center
    }

    fn theta_to_t(&self, theta: f64) -> f64 {
        if !self.rx.approximate_eq(self.ry) {
            todo!("Calculate inverse of elliptic integral");
        }

        let floor_mod = f64::rem_euclid;

        let params = CenterParameterizedArc::from(self);
        let d_theta = params.d_theta;
        let theta1 = params.theta1;

        let reverse = d_theta < 0.0;

        let mut t1_ = floor_mod(theta1, 2.0 * PI);
        let mut t2_ = floor_mod(theta1 + d_theta, 2.0 * PI);
        let mut theta = floor_mod(theta, 2.0 * PI);

        let (mut t1, mut t2) = if reverse { (t2_, t1_) } else { (t1_, t2_) };

        if theta < t1 {
            theta += 2.0 * PI;
        }

        let lambda = (theta - t1) / (t2 - t1);

        if reverse {
            1.0 - lambda
        } else {
            lambda
        }
    }
}

impl Display for Arc {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} -> {} (r: {})", self.start, self.end, self.rx)
    }
}

impl ApproximateEq for Arc {
    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool {
        let other = other.borrow();
        self.start.approximate_eq(other.start)
            && self.end.approximate_eq(other.end)
            && self.rx.approximate_eq(other.rx)
            && self.ry.approximate_eq(other.ry)
            && self.phi.approximate_eq(other.phi)
            && self.sweep == other.sweep
            && self.large_arc == other.large_arc
    }
}

impl<'a> From<&'a Arc> for Ellipse {
    fn from(value: &'a Arc) -> Self {
        let params = CenterParameterizedArc::from(value);

        Ellipse {
            center: params.center,
            rx: params.rx,
            ry: params.ry,
            phi: params.phi,
            clockwise: params.d_theta >= 0.0,
        }
    }
}

impl<'a> From<&'a CenterParameterizedArc> for Arc {
    fn from(value: &'a CenterParameterizedArc) -> Self {
        // Convert an arc in center parameterization to the endpoint parameterization.
        // See https://www.w3.org/TR/SVG/implnote.html#ArcConversionCenterToEndpoint

        let cos = f64::cos;
        let sin = f64::sin;

        // Extract components from the input arc
        let cx = value.center.x;
        let cy = value.center.y;
        let rx = value.rx;
        let ry = value.ry;
        let phi = value.phi;
        let theta_1 = value.theta1;
        let d_theta = value.d_theta;

        // Calculate (start_x, start_y)
        let x1 = cos(phi) * rx * cos(theta_1) - sin(phi) * ry * sin(theta_1) + cx;
        let y1 = sin(phi) * rx * cos(theta_1) + cos(phi) * ry * sin(theta_1) + cy;

        // Calculate (end_x, end_y)
        let x2 =
            cos(phi) * rx * cos(theta_1 + d_theta) - sin(phi) * ry * sin(theta_1 + d_theta) + cx;
        let y2 =
            sin(phi) * rx * cos(theta_1 + d_theta) + cos(phi) * ry * sin(theta_1 + d_theta) + cy;

        // Calculate flags large_arc_flag and sweep_flag
        let f_a = d_theta.abs() > PI;
        let f_s = d_theta > 0.0;

        // Create the output arc in endpoint parameterization
        Self {
            start: Point::new(x1, y1),
            end: Point::new(x2, y2),
            rx,
            ry,
            phi,
            large_arc: f_a,
            sweep: f_s,
        }
    }
}

impl<'a> From<&'a Arc> for CenterParameterizedArc {
    fn from(value: &'a Arc) -> Self {
        let x1 = value.start.x;
        let y1 = value.start.y;
        let x2 = value.end.x;
        let y2 = value.end.y;
        let rx = value.rx;
        let ry = value.ry;
        let phi = value.phi;
        let fa = value.large_arc;
        let fs = value.sweep;

        let a = (x1 - x2) / 2.0;
        let b = (y1 - y2) / 2.0;

        let x1_ = phi.cos() * a + phi.sin() * b;
        let y1_ = -phi.sin() * a + phi.cos() * b;

        let sign = if fa == fs { -1.0 } else { 1.0 };
        let nominator =
            (rx.powi(2) * ry.powi(2) - rx.powi(2) * y1_.powi(2) - ry.powi(2) * x1_.powi(2))
                .max(0.0);
        let denominator = rx.powi(2) * y1_.powi(2) + ry.powi(2) * x1_.powi(2);
        let factor = sign * (nominator / denominator).sqrt();

        let cx_ = factor * rx * y1_ / ry;
        let cy_ = -factor * ry * x1_ / rx;

        let cx = phi.cos() * cx_ - phi.sin() * cy_ + (x1 + x2) / 2.0;
        let cy = phi.sin() * cx_ + phi.cos() * cy_ + (y1 + y2) / 2.0;

        let ellipse = Ellipse {
            center: Point::new(cx, cy),
            rx,
            ry,
            phi,
            clockwise: true,
        };
        let theta1 = ellipse.find_point(value.start);
        let theta2 = ellipse.find_point(value.end);

        let mut delta_theta = (theta2 - theta1) % (2.0 * PI);
        if !fs && delta_theta > 0.0 {
            delta_theta -= 2.0 * PI;
        } else if fs && delta_theta < 0.0 {
            delta_theta += 2.0 * PI;
        }

        CenterParameterizedArc {
            center: Point { x: cx, y: cy },
            rx,
            ry,
            phi,
            theta1,
            d_theta: delta_theta,
        }
    }
}

impl SegmentTrait for Arc {
    fn start(&self) -> Point {
        self.start
    }

    fn end(&self) -> Point {
        self.end
    }

    fn length(&self) -> f64 {
        // Calculate the length of the ellipse segment with the incomplete elliptic integral
        // of the second kind.
        let p = CenterParameterizedArc::from(self);

        let a = p.rx;
        let b = p.ry;
        let psi_1 = p.theta1 - p.phi;
        let psi_2 = p.theta1 + p.d_theta - p.phi;
        let arg2 = (1.0 - a.powi(2) / b.powi(2)).sqrt();

        (b * (ellie(psi_2, arg2) - ellie(psi_1, arg2))).abs()
    }

    fn reversed(&self) -> Self {
        Arc {
            start: self.end,
            end: self.start,
            rx: self.rx,
            ry: self.ry,
            phi: self.phi,
            large_arc: self.large_arc,
            sweep: self.sweep,
        }
    }

    fn area_contribution(&self) -> f64 {
        let exp = Complex64::exp;
        let j = Complex64::I;
        let abs = ComplexFloat::abs;

        // Obtain ellipse in form y(t) = c + u * exp(it) + v * exp(-it)
        let c_p: CenterParameterizedArc = self.into();
        let c: Complex64 = c_p.center.into();
        let rx = c_p.rx;
        let ry = c_p.ry;
        let phi = c_p.phi;
        let u = exp(j * phi) * (rx + ry) / 2.0;
        let v = exp(j * phi) * (rx - ry) / 2.0;
        let theta1 = c_p.theta1;
        let theta2 = c_p.theta1 + c_p.d_theta;

        // Find k, l, and m for the area calculation
        let k = 2.0 * (abs(u).powi(2) - abs(v).powi(2));
        let l = c.conj() * u - c * v.conj();
        let m = c * u.conj() - c.conj() * v;

        // Calculate area contribution
        let a = 0.25 * (k * theta2 - j * l * exp(j * theta2) + j * m * exp(-j * theta2))
            - 0.25 * (k * theta1 - j * l * exp(j * theta1) + j * m * exp(-j * theta1));

        a.re
    }

    fn point_at_t(&self, t: f64) -> Point {
        let exp = Complex64::exp;
        let j = Complex64::I;

        let params = CenterParameterizedArc::from(self);

        let phi = params.phi;
        let rx = params.rx;
        let ry = params.ry;
        let theta1 = params.theta1;
        let d_theta = params.d_theta;

        let u = exp(j * phi) * (rx + ry) / 2.0;
        let v = exp(j * phi) * (rx - ry) / 2.0;

        (u * exp(j * (theta1 + t * d_theta)) + v * exp(-j * (theta1 + t * d_theta))).into()
    }

    fn t_of_point<T: Borrow<Point>>(&self, point: T) -> f64 {
        let point = point.borrow();

        if self.start.approximate_eq(point) {
            return 0.0;
        } else if self.end.approximate_eq(point) {
            return 1.0;
        }

        let ellipse = Ellipse::from(self);
        let theta = ellipse.find_point(point);

        let t = self.theta_to_t(theta);

        if (0.0..=1.0).contains(&t) {
            t
        } else {
            f64::NAN
        }
    }

    fn slice<T: Into<f64>, U: Into<f64>>(&self, t1: T, t2: U) -> Self {
        let t1 = t1.into();
        let t2 = t2.into();

        let p = CenterParameterizedArc::from(self);
        let new_params = CenterParameterizedArc {
            center: p.center,
            rx: p.rx,
            ry: p.ry,
            phi: p.phi,
            theta1: p.theta1 + t1 * p.d_theta,
            d_theta: (t2 - t1) * p.d_theta,
        };
        Self::from(&new_params)
    }

    fn intersect<T: Into<Segment>>(&self, other: T) -> Vec<IntersectionSolution> {
        let other = other.into();
        match other {
            Segment::Line(other) => vec![],
            Segment::Arc(other) => vec![],
        }
    }
}
