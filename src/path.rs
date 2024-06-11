use crate::geom::Point;
use crate::geom::Segment;
use crate::geom::SegmentTrait;
use crate::helpers::ApproximateEq;

#[derive(Clone, Debug)]
pub struct Path {
    segments: Vec<Segment>,
    closed: bool,
}

impl Path {
    pub fn new(segments: &[Segment], closed: bool) -> Self {
        Path {
            segments: Vec::from(segments),
            closed,
        }
    }

    pub fn segments(&self) -> &[Segment] {
        self.segments.as_slice()
    }

    pub fn length(&self) -> f64 {
        self.segments
            .iter()
            .fold(0.0, |acc, segment| acc + segment.length())
    }

    pub fn area(&self) -> f64 {
        self.segments
            .iter()
            .fold(0.0, |acc, segment| acc + segment.area_contribution())
    }

    pub fn is_clockwise(&self) -> bool {
        self.area() >= 0.0
    }

    pub fn slice(&self, t1: f64, t2: f64) -> Path {
        if t1.approximate_eq(t2) {
            return Path {
                segments: vec![],
                closed: false,
            };
        }

        let length = self.length();
        let mut l = 0.0;
        let mut result_segments: Vec<Segment> = Vec::new();

        for segment in &self.segments {
            let segment_length = segment.length();

            if t2 * length < l || t1 * length > l + segment_length {
                l += segment_length;
                continue;
            }

            let mut seg_t1 = 0.0;
            let mut seg_t2 = 1.0;

            if t1 * length >= l && t1 * length < l + segment_length {
                seg_t1 = (t1 * length - l) / segment_length;
            }

            if t2 * length >= l && t2 * length < l + segment_length {
                seg_t2 = (t2 * length - l) / segment_length;
            }

            result_segments.push(segment.slice(seg_t1, seg_t2));
            l += segment_length;
        }

        Path {
            segments: result_segments,
            closed: false,
        }
    }
    pub fn to_points(&self) -> Vec<Point> {
        let mut result: Vec<Point> = Vec::with_capacity(self.segments.len() + 1);
        result.push(self.segments[0].start());
        for segment in self.segments.iter() {
            result.push(segment.end())
        }
        result
    }
}

#[macro_export]
macro_rules! path {
    // Case: closed keyword argument provided
    ($($x:expr),*; closed: $closed:expr) => (
        crate::path::Path::new(&[$($x.into()),+], $closed)
    );

    // Case: no closed keyword argument provided (default to false)
    ($($x:expr),+) => (
        crate::path::Path::new(&[$($x.into()),+], false)
    );
}
