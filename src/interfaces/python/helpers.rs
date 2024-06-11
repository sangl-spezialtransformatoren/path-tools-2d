use crate::geom::Segment;
use crate::interfaces::python::PySegment;
use pyo3::types::{PyIterator, PyList};
use pyo3::{PyObject, PyResult, Python};

pub trait CastPyList {
    fn to_segment_vec(&self) -> PyResult<Vec<Segment>>;
    fn from_segment_vec(segments: &[Segment], py: Python) -> PyResult<Vec<PyObject>>;
}

impl CastPyList for PyList {
    fn to_segment_vec(&self) -> PyResult<Vec<Segment>> {
        self.iter()
            .map(|element| {
                let py_segment: PySegment = element.extract()?;
                Ok(py_segment.segment.clone())
            })
            .collect()
    }

    fn from_segment_vec(segments: &[Segment], py: Python) -> PyResult<Vec<PyObject>> {
        segments
            .into_iter()
            .map(|segment| {
                let py_segment = PySegment::new(segment.clone()).to_py_object(py)?;
                Ok(py_segment)
            })
            .collect()
    }
}
