use std::borrow::Borrow;

use pyo3::IntoPy;
use pyo3::Py;
use pyo3::PyAny;
use pyo3::pyclass;
use pyo3::pyclass_init::PyObjectInit;
use pyo3::pymethods;
use pyo3::PyObject;
use pyo3::PyResult;
use pyo3::Python;
use pyo3::ToPyObject;

use crate::algorithms::IntersectionSolution;
use crate::geom::{Line, Point, Segment, SegmentTrait};
use crate::interfaces::python::arc::PyArc;
use crate::interfaces::python::line::PyLine;
use crate::interfaces::python::point::PointLike;

#[derive(Clone)]
#[pyclass(name = "Segment", subclass)]
pub struct PySegment {
    pub segment: Segment,
}

impl PySegment {
    pub fn new(segment: Segment) -> Self {
        PySegment { segment }
    }

    pub(crate) fn to_py_object(&self, py: Python) -> PyResult<PyObject> {
        match &self.segment {
            Segment::Line(line) => {
                let init = PyLine::py_initializer(line.clone());
                Ok(Py::new(py, init)?.to_object(py))
            }
            Segment::Arc(arc) => {
                let init = PyArc::py_initializer(arc.clone());
                Ok(Py::new(py, init)?.to_object(py))
            }
        }
    }
}

#[pymethods]
impl PySegment {
    #[getter]
    fn get_length(&self) -> f64 {
        self.segment.length()
    }

    #[getter]
    fn get_start(&self) -> Point {
        self.segment.start()
    }

    #[getter]
    fn get_end(&self) -> Point {
        self.segment.end()
    }

    fn reversed(&self, py: Python) -> PyResult<PyObject> {
        let reversed_segment = self.segment.reversed();
        PySegment::new(reversed_segment).to_py_object(py)
    }

    fn point_at_t(&self, t: f64) -> Point {
        self.segment.point_at_t(t)
    }

    fn py_t_of_point(&self, point: &PyAny) -> f64 {
        let point: Point = point
            .extract::<PointLike>()
            .expect("Can't cast to point")
            .into();

        self.segment.t_of_point(point)
    }

    fn slice(&self, t1: f64, t2: f64, py: Python) -> PyResult<PyObject> {
        let sliced_segment = self.segment.slice(t1, t2);
        PySegment::new(sliced_segment).to_py_object(py)
    }

    fn contains(&self, point: &PyAny) -> bool {
        let point: Point = point
            .extract::<PointLike>()
            .expect("Can't cast to point")
            .into();
        self.segment.contains(point)
    }

    fn intersect(&self, other: &PySegment) -> Vec<Point> {
        let solutions: Vec<IntersectionSolution> = self.segment.intersect(other.segment);
        solutions.iter().map(|solution| solution.point).collect()
    }
}
