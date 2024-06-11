use num_complex::Complex;

use crate::helpers::ApproximateEq;

pub trait IsReal {
    fn is_real(&self) -> bool;
}

impl IsReal for Complex<f64> {
    fn is_real(&self) -> bool {
        self.im.approximate_eq(0.0)
    }
}
