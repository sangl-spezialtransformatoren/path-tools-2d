use std::panic::resume_unwind;

use geo::{coord, line_string, Coord, LineString, Simplify};
use pyo3::types::{PyAny, PyIterator, PyList};
use pyo3::{pyclass, pyfunction, pymethods, PyObject, PyResult, Python};
use spec_math::cephes64::polevl;

use crate::algorithms::combine_segments;
use crate::geom::SegmentTrait;
use crate::geom::{Line, Segment};
use crate::interfaces::python::helpers::CastPyList;
use crate::interfaces::python::line::PyLine;
use crate::interfaces::python::PySegment;
#[pyclass]
pub struct PathBuilder {
    segments: Vec<Segment>,
}

#[pymethods]
impl PathBuilder {
    #[new]
    pub fn py_new(size: usize) -> Self {
        PathBuilder {
            segments: Vec::with_capacity(size),
        }
    }

    pub fn append(&mut self, segment: &PyAny) -> PyResult<()> {
        let py_segment: PySegment = segment.extract()?;
        self.segments.push(py_segment.segment);

        return Ok(());
    }

    pub fn segments(&self, py: Python) -> PyResult<Vec<PyObject>> {
        PyList::from_segment_vec(&self.segments, py)
    }
}

#[pyfunction(name = "combine_segments")]
pub fn py_combine_segments(segments: &PyList, py: Python) -> PyResult<Vec<Vec<PyObject>>> {
    let segments: Vec<Segment> = segments.to_segment_vec()?;
    let res = combine_segments(&segments);

    res.into_iter()
        .map(|strain| PyList::from_segment_vec(&strain, py))
        .collect()
}

#[pyfunction]
pub fn ramer_douglas_peucker(segments: &PyList, tol: f64, py: Python) -> PyResult<Vec<PyObject>> {
    let segments: Vec<Segment> = segments.to_segment_vec()?;
    let mut coords: Vec<geo::Coord> = Vec::with_capacity(segments.len() + 1);

    segments.iter().for_each(|segment| {
        coords.push(Coord {
            x: segment.start().x,
            y: segment.start().y,
        })
    });
    coords.push(Coord {
        x: segments.last().unwrap().end().x,
        y: segments.last().unwrap().end().y,
    });

    let line_string = LineString::new(coords);

    let simplified_line_string = line_string.simplify(&tol);
    let lines: Vec<Segment> = simplified_line_string
        .lines()
        .map(|line| Line::new(line.start.x_y(), line.end.x_y()).into())
        .collect();

    PyList::from_segment_vec(&lines, py)
}

#[pyfunction]
pub fn points_to_beziers(
    points: Vec<(f64, f64)>,
    tol: f64,
) -> Vec<((f64, f64), (f64, f64), (f64, f64), (f64, f64))> {
    let points: Vec<simplify_rs::Point> = points
        .into_iter()
        .map(|(x, y)| simplify_rs::Point { x, y })
        .collect();
    let simplified_points = simplify_rs::simplify(&points, tol);
    let mut result: Vec<((f64, f64), (f64, f64), (f64, f64), (f64, f64))> =
        Vec::with_capacity(simplified_points.len() / 4);

    for i in 0..simplified_points.len() / 4 {
        let p1 = simplified_points[i * 4];
        let p2 = simplified_points[i * 4 + 1];
        let p3 = simplified_points[i * 4 + 2];
        let p4 = simplified_points[i * 4 + 3];
        result.push(((p1.x, p1.y), (p2.x, p2.y), (p3.x, p3.y), (p4.x, p4.y)))
    }

    result
}
