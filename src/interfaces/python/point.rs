use pyo3::{pymethods, FromPyObject, PyAny, PyResult};

use crate::geom::Point;

#[derive(Clone)]
pub struct PointLike(pub Point);

impl From<PointLike> for Point {
    fn from(value: PointLike) -> Self {
        value.0
    }
}

impl<'a> FromPyObject<'a> for PointLike {
    fn extract(ob: &'a PyAny) -> PyResult<Self> {
        if let Ok(point) = ob.extract::<Point>() {
            Ok(PointLike(point))
        } else if let Ok((x, y)) = ob.extract::<(f64, f64)>() {
            Ok(PointLike(Point::new(x, y)))
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "Invalid input for PointLike",
            ))
        }
    }
}

#[pymethods]
impl Point {
    #[new]
    fn py_new(x: f64, y: f64) -> Self {
        Point::new(x, y)
    }

    #[getter]
    fn get_x(&self) -> f64 {
        self.x
    }

    #[getter]
    fn get_y(&self) -> f64 {
        self.y
    }

    fn to_tuple(&self) -> (f64, f64) {
        (self.x, self.y)
    }

    fn __repr__(&self) -> String {
        format!("Point({}, {})", self.x, self.y)
    }
}
