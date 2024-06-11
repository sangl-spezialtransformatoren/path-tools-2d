extern crate nalgebra as na;
use crate::geom::Point;
use na::{DMatrix, SVD};
use nalgebra::ComplexField;
use std::f64::consts::PI;

pub fn fit_ellipse(points: &[Point]) -> [f64; 6] {
    let n = points.len();
    let mut d = DMatrix::zeros(n, 6);

    for (i, &point) in points.iter().enumerate() {
        let (x, y) = (point.x, point.y);
        d[(i, 0)] = x * x;
        d[(i, 1)] = x * y;
        d[(i, 2)] = y * y;
        d[(i, 3)] = x;
        d[(i, 4)] = y;
        d[(i, 5)] = 1.0;
    }

    let svd = SVD::new(d.clone(), true, true);
    let a = svd.v_t.unwrap().row(5).transpose();
    let a_vec = a.iter().cloned().collect::<Vec<_>>();

    [a_vec[0], a_vec[1], a_vec[2], a_vec[3], a_vec[4], a_vec[5]]
}

pub fn ellipse_parameters(a: [f64; 6]) -> (f64, f64, f64, f64, f64) {
    let [a, b, c, d, e, f] = a;

    // https://en.wikipedia.org/wiki/Ellipse#General_ellipse

    // Calculate the ellipse center
    let x0 = (2.0 * c * d - b * e) / (b * b - 4.0 * a * c);
    let y0 = (2.0 * a * e - b * d) / (b * b - 4.0 * a * c);

    // Rotation angle
    let phi = 0.5 * f64::atan2(-b, c - a);

    // Calculate the semi-major and semi-minor axes
    let dis = b * b - 4.0 * a * c;
    let u = 2.0 * (a * e * e + c * d * d - b * d * e + dis * f);
    let v = ((a - c).powi(2) + b.powi(2)).sqrt();

    let maj = -(u * ((a + c) + v)).sqrt() / dis;
    let min = -(u * ((a + c) - v)).sqrt() / dis;

    (x0, y0, maj, min, phi)
}
