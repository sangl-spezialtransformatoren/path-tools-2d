use pyo3::impl_::pymethods::IterBaseKind;
use pyo3::prelude::*;
use pyo3::prelude::{PyModule, PyModuleMethods};
use pyo3::types::{PyCapsule, PyIterator, PyList};
use pyo3::{
    pyfunction, pymodule, wrap_pyfunction, Bound, PyAny, PyObject, PyRef, PyResult, Python,
};
use std::ffi::c_void;
use std::process::exit;

use helpers::CastPyList;
pub use segment::PySegment;

use crate::algorithms::combine_segments;
use crate::geom::{Arc, Line, Point, Segment};
use crate::interfaces::python::algorithms::*;
use crate::interfaces::python::arc::PyArc;
use crate::interfaces::python::line::PyLine;

mod algorithms;
mod arc;
mod helpers;
mod line;
mod point;
mod segment;

#[pymodule]
fn rs_paths(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Point>()?;
    m.add_class::<PyLine>()?;
    m.add_class::<PyArc>()?;
    m.add_class::<PySegment>()?;
    m.add_class::<PathBuilder>()?;
    m.add_function(wrap_pyfunction!(py_combine_segments, m)?)?;
    m.add_function(wrap_pyfunction!(ramer_douglas_peucker, m)?)?;
    m.add_function(wrap_pyfunction!(points_to_beziers, m)?)?;
    Ok(())
}
