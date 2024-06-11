use std::ops::Deref;

use pyo3::{Py, PyAny, pyclass, PyClassInitializer, pymethods, PyRef, PyResult};

use crate::geom::{Line, Point, Segment, SegmentTrait};
use crate::interfaces::python::point::PointLike;
use crate::interfaces::python::PySegment;

#[pyclass(name="Line", extends=PySegment)]
pub struct PyLine {}

impl PyLine {
    pub fn py_initializer(line: Line) -> PyClassInitializer<Self> {
        PyClassInitializer::from(PySegment::new(Segment::Line(line))).add_subclass(PyLine {})
    }
}

impl From<PyRef<'_, PyLine>> for Line {
    fn from(value: PyRef<'_, PyLine>) -> Self {
        let segment = value.into_super().deref().segment;
        match segment {
            Segment::Line(line) => line,
            _ => unreachable!(),
        }
    }
}

#[pymethods]
impl PyLine {
    #[new]
    fn py_new(a: &PyAny, b: &PyAny) -> PyClassInitializer<Self> {
        let a: Point = a
            .extract::<PointLike>()
            .expect("Can't cast to point.")
            .into();
        let b: Point = b
            .extract::<PointLike>()
            .expect("Can't cast to point.")
            .into();

        let line = Line::new(a, b);
        PyLine::py_initializer(line)
    }

    fn __repr__(self_: PyRef<'_, PyLine>) -> String {
        let line: Line = self_.into();
        let start = line.start();
        let end = line.end();

        format!("Line(({}, {}) -> ({}, {}))", start.x, start.y, end.x, end.y)
    }
}
