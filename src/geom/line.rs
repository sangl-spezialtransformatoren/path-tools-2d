use std::borrow::Borrow;
use std::f64::NAN;
use std::fmt::Display;
use std::fmt::Formatter;

use num_complex::Complex64;
use num_complex::ComplexFloat;

use crate::algorithms::IntersectionSolution;
use crate::geom::Point;
use crate::geom::Segment;
use crate::geom::SegmentTrait;
use crate::geom::Straight;
use crate::helpers::ApproximateEq;

#[derive(Copy, Clone, Debug)]
pub struct Line {
    start: Point,
    end: Point,
}

impl Line {
    pub fn new<T: Into<Point>, U: Into<Point>>(start: T, end: U) -> Self {
        let start = start.into();
        let end = end.into();
        Line { start, end }
    }
}

impl Display for Line {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} -> {}", self.start, self.end)
    }
}

impl ApproximateEq for Line {
    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool {
        let other = other.borrow();
        self.start.approximate_eq(other.start) && self.end.approximate_eq(other.end)
    }
}

impl<'a> From<&'a Line> for Straight {
    fn from(value: &'a Line) -> Self {
        Straight {
            a: value.start,
            b: value.end,
        }
    }
}

impl SegmentTrait for Line {
    fn start(&self) -> Point {
        self.start
    }

    fn end(&self) -> Point {
        self.end
    }

    fn length(&self) -> f64 {
        let start: Complex64 = self.start.into();
        let end: Complex64 = self.end.into();

        (end - start).abs()
    }

    fn reversed(&self) -> Line {
        Line {
            start: self.end,
            end: self.start,
        }
    }

    fn area_contribution(&self) -> f64 {
        let (px, py) = (self.start.x, self.start.y);
        let (qx, qy) = (self.end.x, self.end.y);

        0.5 * (px * qy - py * qx)
    }

    fn point_at_t(&self, t: f64) -> Point {
        let start: Complex64 = self.start.into();
        let end: Complex64 = self.end.into();

        Point::from(start + t * (end - start))
    }

    fn t_of_point<T: Borrow<Point>>(&self, point: T) -> f64 {
        let point = point.borrow();

        let straight = Straight::from(self);
        let lambda = straight.find_point(point);

        if !lambda.is_nan() && (0.0..=1.0).contains(&lambda) {
            lambda
        } else {
            NAN
        }
    }

    fn slice<T: Into<f64>, U: Into<f64>>(&self, t1: T, t2: U) -> Self {
        let t1 = t1.into();
        let t2 = t2.into();

        let start: Complex64 = self.start.into();
        let end: Complex64 = self.end.into();

        Self {
            start: (start + t1 * (end - start)).into(),
            end: (start + t2 * (end - start)).into(),
        }
    }

    fn intersect<T: Into<Segment>>(&self, other: T) -> Vec<IntersectionSolution> {
        let other = other.into();
        match other {
            Segment::Line(other) => vec![],
            Segment::Arc(other) => vec![],
        }
    }
}
