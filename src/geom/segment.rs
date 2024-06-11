use std::borrow::Borrow;

use derive_more::Display;
use enum_dispatch::enum_dispatch;

use crate::algorithms::IntersectionSolution;
use crate::geom::Arc;
use crate::geom::Line;
use crate::geom::Point;
use crate::helpers::ApproximateEq;

#[enum_dispatch]
pub trait SegmentTrait: ApproximateEq {
    fn start(&self) -> Point;
    fn end(&self) -> Point;
    fn length(&self) -> f64;
    fn reversed(&self) -> Self;
    fn area_contribution(&self) -> f64;
    fn point_at_t(&self, t: f64) -> Point;
    fn t_of_point<T: Borrow<Point>>(&self, point: T) -> f64;
    fn slice<T: Into<f64>, U: Into<f64>>(&self, t1: T, t2: U) -> Self;
    fn intersect<T: Into<Segment>>(&self, other: T) -> Vec<IntersectionSolution>;

    fn contains<T: Into<Point>>(&self, point: T) -> bool {
        let point = point.into();
        let t = self.t_of_point(point);
        !t.is_nan()
    }
}

#[derive(Display, Copy, Clone, Debug)]
#[enum_dispatch(SegmentTrait)]
pub enum Segment {
    Line(Line),
    Arc(Arc),
}

impl ApproximateEq for Segment {
    fn approximate_eq<T: Borrow<Self>>(&self, other: T) -> bool {
        let other = other.borrow();
        match (self, other) {
            (Segment::Line(a), Segment::Line(b)) => a.approximate_eq(b),
            (Segment::Arc(a), Segment::Arc(b)) => a.approximate_eq(b),
            _ => false,
        }
    }
}
