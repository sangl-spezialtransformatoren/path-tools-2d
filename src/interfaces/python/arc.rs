use std::ops::Deref;

use pyo3::{
    Bound, Py, PyAny, pyclass, PyClassInitializer, pymethods, PyObject, PyRef, PyResult, Python,
};
use pyo3::types::PyType;

use crate::geom::{Arc, CenterParameterizedArc, Line, Point, Segment, SegmentTrait};
use crate::interfaces::python::line::PyLine;
use crate::interfaces::python::point::PointLike;
use crate::interfaces::python::PySegment;

#[pyclass(name="Arc", extends=PySegment)]
pub struct PyArc {}

impl PyArc {
    pub fn py_initializer(arc: Arc) -> PyClassInitializer<Self> {
        PyClassInitializer::from(PySegment::new(Segment::Arc(arc))).add_subclass(PyArc {})
    }
}

impl From<PyRef<'_, PyArc>> for Arc {
    fn from(value: PyRef<'_, PyArc>) -> Self {
        let segment = (*value.into_super()).segment;
        match segment {
            Segment::Arc(arc) => arc,
            _ => unreachable!(), // This should never happen if used correctly
        }
    }
}

#[pymethods]
impl PyArc {
    #[new]
    fn py_new(
        a: &PyAny,
        b: &PyAny,
        rx: f64,
        ry: f64,
        phi: f64,
        large_arc: bool,
        sweep: bool,
    ) -> PyClassInitializer<Self> {
        let a: Point = a
            .extract::<PointLike>()
            .expect("Can't cast to point.")
            .into();

        let b: Point = b
            .extract::<PointLike>()
            .expect("Can't cast to point.")
            .into();

        let arc = Arc::new(a, b, rx, ry, phi, large_arc, sweep);
        PyArc::py_initializer(arc)
    }

    #[classmethod]
    fn from_center_params(
        cls: &Bound<'_, PyType>,
        center: &PyAny,
        rx: f64,
        ry: f64,
        phi: f64,
        theta1: f64,
        d_theta: f64,
        py: Python,
    ) -> PyResult<PyObject> {
        let center: Point = center
            .extract::<PointLike>()
            .expect("Can't cast to point.")
            .into();
        let arc: Arc = (&CenterParameterizedArc::new(center, rx, ry, phi, theta1, d_theta)).into();

        PySegment::new(Segment::Arc(arc)).to_py_object(py)
    }

    #[getter]
    fn get_rx(self_: PyRef<'_, PyArc>) -> f64 {
        let arc: Arc = self_.into();
        arc.rx()
    }

    #[getter]
    fn get_ry(self_: PyRef<'_, PyArc>) -> f64 {
        let arc: Arc = self_.into();
        arc.ry()
    }

    #[getter]
    fn get_phi(self_: PyRef<'_, PyArc>) -> f64 {
        let arc: Arc = self_.into();
        arc.phi()
    }

    #[getter]
    fn get_large_arc(self_: PyRef<'_, PyArc>) -> bool {
        let arc: Arc = self_.into();
        arc.large_arc()
    }

    #[getter]
    fn get_sweep(self_: PyRef<'_, PyArc>) -> bool {
        let arc: Arc = self_.into();
        arc.sweep()
    }
    
    fn get_center(self_: PyRef<'_, PyArc>) -> Point {
        let arc: Arc = self_.into();
        arc.center()
    }
    
    fn __repr__(self_: PyRef<'_, PyArc>) -> String {
        let arc: Arc = self_.into();
        let start = arc.start();
        let end = arc.end();
        format!(
            "Arc(({}, {}) -> ({}, {}), rx: {}, ry: {}, phi: {}, large_arc: {}, sweep: {})",
            start.x,
            start.y,
            end.x,
            end.y,
            arc.rx(),
            arc.ry(),
            arc.phi(),
            arc.large_arc(),
            arc.sweep()
        )
    }
}
