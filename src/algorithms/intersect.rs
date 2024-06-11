use std::fmt::Display;
use std::fmt::Formatter;

use num_complex::Complex;
use num_complex::Complex64;
use num_complex::ComplexFloat;

use crate::geom::Ellipse;
use crate::geom::Line;
use crate::geom::Point;
use crate::geom::Segment;
use crate::geom::SegmentTrait;
use crate::geom::Straight;
use crate::helpers;
use crate::helpers::ApproximateEq;
use crate::helpers::IsReal;

#[derive(Copy, Clone, Debug)]
pub struct IntersectionSolution {
    pub(crate) point: Point,
    t1: f64,
    t2: f64,
}

impl Display for IntersectionSolution {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {}, {}]", self.point, self.t1, self.t2)
    }
}

fn swap_ts(solutions: &[IntersectionSolution]) -> Vec<IntersectionSolution> {
    solutions
        .iter()
        .map(|sol| IntersectionSolution {
            point: sol.point,
            t1: sol.t2,
            t2: sol.t1,
        })
        .collect()
}

pub fn intersect_straights<T: Into<Straight>, U: Into<Straight>>(
    v: T,
    w: U,
) -> Vec<IntersectionSolution> {
    let v = v.into();
    let w = w.into();

    let z1: Complex64 = v.a.into();
    let z2: Complex64 = v.b.into();
    let z3: Complex64 = w.a.into();
    let z4: Complex64 = w.b.into();

    let numerator = (z2.conj() * z1 - z2 * z1.conj()) * (z4 - z3)
        - (z4.conj() * z3 - z4 * z3.conj()) * (z2 - z1);

    let denominator = (z4 - z3) * (z2.conj() - z1.conj()) - (z2 - z1) * (z4.conj() - z3.conj());

    // Handle parallel lines
    if denominator.approximate_eq(Complex::from(0.0)) {
        // Check if lines are collinear
        if ((z3 - z1) / (z2 - z1)).is_real() {
            panic!("Specified straights are collinear");
        }

        // Otherwise, the straights are parallel and have no intersection
        return Vec::new();
    }

    // Lines are not parallel -> Intersection point exists
    let x = numerator / denominator;

    // Calculate parameters
    let lambda = (x - z1) / (z2 - z1);
    let kappa = (x - z3) / (z4 - z3);

    if !lambda.is_real() && kappa.is_real() {
        println!("V {}", v);
        println!("W {}", w);
        panic!("Error in straight intersection");
    }

    vec![IntersectionSolution {
        point: Point::from(x),
        t1: lambda.re,
        t2: kappa.re,
    }]
}

pub fn intersect_straight_with_ellipse(
    straight: &Straight,
    ellipse: &Ellipse,
) -> Vec<IntersectionSolution> {
    // Solves the equation
    //   a + λ (b - a) = c + u exp(iθ) + v exp(-iθ)
    // for λ and θ.
    //
    // u and v are defined as:
    //   u := exp(iφ) (rₓ + rᵧ) / 2
    //   v := exp(iφ) (rₓ - rᵧ) / 2

    // Some constants/functions
    let i = Complex64::I;
    let c_exp = Complex64::exp;
    let ln = Complex64::ln;
    let sqrt = Complex64::sqrt;

    // Our variables
    let a: Complex64 = straight.a.into();
    let b: Complex64 = straight.b.into();
    let c: Complex64 = ellipse.center.into();
    let phi = ellipse.phi;
    let rx = ellipse.rx;
    let ry = ellipse.ry;

    // Calculate alternative representation
    let u = c_exp(i * phi) * (rx + ry) / 2.0;
    let v = c_exp(i * phi) * (rx - ry) / 2.0;

    // Introduce auxiliary variables
    let a_ = (c - a) / (b - a);
    let b_ = u / (b - a);
    let c_ = v / (b - a);

    // Solve for θ
    let theta_plus = -i
        * ln(
            (-i * a_.im + sqrt(-(a_.im).powi(2) - (b_ - c_.conj()) * (c_ - b_.conj())))
                / (b_ - c_.conj()),
        );
    let theta_minus = -i
        * ln(
            (-i * a_.im - sqrt(-(a_.im).powi(2) - (b_ - c_.conj()) * (c_ - b_.conj())))
                / (b_ - c_.conj()),
        );

    // Non-real solutions mean that there is no intersection
    if !theta_plus.is_real() || !theta_minus.is_real() {
        return Vec::new();
    }

    let t_plus = theta_plus.re;
    let t_minus = theta_minus.re;

    // Calculate λᵢ and intersection points xᵢ
    let lambda_plus = a_.re
        + (b_ + c_.conj()) * c_exp(i * t_plus) / 2.0
        + (c_ + b_.conj()) * c_exp(-i * t_plus) / 2.0;
    let lambda_minus = a_.re
        + (b_ + c_.conj()) * c_exp(i * t_minus) / 2.0
        + (c_ + b_.conj()) * c_exp(-i * t_minus) / 2.0;

    let x_plus = a + lambda_plus * (b - a);
    let x_minus = a + lambda_minus * (b - a);

    if !lambda_plus.is_real() && lambda_minus.is_real() {
        println!("V: {:?}", v);
        println!("E: {:?}", ellipse);
        panic!("Error while intersecting straight and ellipse");
    }

    // V is tangential to E
    if x_plus.approximate_eq(x_minus) {
        return vec![IntersectionSolution {
            point: Point::from(x_plus),
            t1: lambda_plus.re,
            t2: theta_plus.re,
        }];
    }

    // We have two solutions
    vec![
        IntersectionSolution {
            point: x_plus.into(),
            t1: lambda_plus.re,
            t2: theta_plus.re,
        },
        IntersectionSolution {
            point: x_minus.into(),
            t1: lambda_minus.re,
            t2: theta_minus.re,
        },
    ]
}

pub fn intersect_ellipse_with_straight(
    ellipse: &Ellipse,
    straight: &Straight,
) -> Vec<IntersectionSolution> {
    swap_ts(&intersect_straight_with_ellipse(straight, ellipse))
}

pub fn intersect_ellipses(e1: &Ellipse, e2: &Ellipse) -> Vec<IntersectionSolution> {
    // To calculate the intersections of two ellipses, the following approach is used:
    //   1. Transform the space so that the first ellipse is the unit circle
    //      This is done by setting the two ellipses equal, splitting in real and imaginary part, scaling it so that both
    //      left hand sides become 1 and adding the parts back together. On the left hand side there is now a unit circle
    //      and on the right hand side a new ellipse.
    //   2. Solve |c + u exp(iθ) + v exp(-iθ)| = 1 for the new ellipse to find the solutions for θ.

    let i = Complex64::I;
    let exp = Complex64::exp;
    let ln = Complex64::ln;
    let complex_zero = Complex::new(0.0, 0.0);
    let sqrt = Complex64::sqrt;

    // Parameters of the first ellipse
    let c1: Complex64 = e1.center.into();
    let phi1 = e1.phi;
    let rx1 = e1.rx;
    let ry1 = e1.ry;

    // Parameters of the second ellipse
    let c2: Complex64 = e2.center.into();
    let phi2 = e2.phi;
    let rx2 = e2.rx;
    let ry2 = e2.ry;

    // Find alternative representation of second ellipse
    // E2(θ) = c2 + u2 exp(iθ) + v2 exp(-iθ)
    let u2 = exp(i * phi2) * (rx2 + ry2) / 2.0;
    let v2 = exp(i * phi2) * (rx2 - ry2) / 2.0;

    // Introduce auxiliary variables
    let f_t = 2.0 * exp(-i * phi1) * (c2 - c1);
    let f_u = 2.0 * exp(-i * phi1) * u2;
    let f_v = 2.0 * exp(-i * phi1) * v2;

    // Find transformed second ellipse
    let t_c = f_t.re / (2.0 * rx1) + i * f_t.im / (2.0 * ry1);
    let t_u = 1.0 / 4.0 * ((f_u + f_v.conj()) / rx1 + (f_u - f_v.conj()) / ry1);
    let t_v = 1.0 / 4.0 * ((f_u.conj() + f_v) / rx1 - (f_u.conj() - f_v) / ry1);

    if t_c.approximate_eq(complex_zero)
        && (t_u * t_v.conj()).approximate_eq(complex_zero)
        && (t_u.abs().powi(2) + t_v.abs().powi(2)).approximate_eq(1.0)
    {
        panic!("Ellipses are identical");
    }

    // Solving |c + u exp(iθ) + v exp(-iθ)| = 1 yields a quartic equation
    //   0 = Ax⁴ + Bx³ + Cx² + Dx + E
    // where x is exp(iθ).

    // The coefficients are
    let f_a = t_u * t_v.conj();
    let f_b = t_c * t_v.conj() + t_c.conj() * t_u;
    let f_c = t_c * t_c.conj() + t_u.abs() + t_v.abs() - 1.0;
    let f_d = t_c * t_u.conj() + t_c.conj() * t_v;
    let f_e = t_v * t_u.conj();

    // There are at most 4 solutions

    // If A is 0, also E is 0 as E is conjugate(A). The problem then transitions into a quadratic equation
    // Bx² + Cx + D
    let solutions: Vec<Complex64> = if f_a.approximate_eq(complex_zero) {
        if !f_e.approximate_eq(complex_zero) {
            panic!("Error while intersecting ellipses");
        }
        let x1 = (-f_c + (f_c.powi(2) - 4.0 * f_b * f_d).sqrt()) / (2.0 * f_b);
        let x2 = (-f_c - (f_c.powi(2) - 4.0 * f_b * f_d).sqrt()) / (2.0 * f_b);
        let theta1 = -i * ln(x1);
        let theta2 = -i * ln(x2);
        vec![theta1, theta2]
    } else {
        // A quartic equation has up to four solutions which can be obtained with Ferrari's formula.
        let alpha = -3.0 * f_b.powi(2) / (8.0 * f_a.powi(2)) + f_c / f_a;
        let beta = f_b.powi(3) / (8.0 * f_a.powi(3)) - f_b * f_c / (2.0 * f_a.powi(2)) + f_d / f_a;
        let gamma = -3.0 * f_b.powi(4) / (256.0 * f_a.powi(4))
            + f_b.powi(2) * f_c / (16.0 * f_a.powi(3))
            - f_b * f_d / (4.0 * f_a.powi(2))
            + f_e / f_a;

        let p = -alpha.powi(2) / 12.0 - gamma;
        let q = -alpha.powi(3) / 108.0 + alpha * gamma / 3.0 - beta.powi(2) / 8.0;
        let u = (-q / 2.0 + (q.powi(2) / 4.0 + p.powi(3) / 27.0).sqrt()).powf(1.0 / 3.0);
        let y = -5.0 / 6.0 * alpha
            + if p.approximate_eq(complex_zero) {
                -q.powf(1.0 / 3.0)
            } else {
                u - p / (3.0 * u)
            };
        let w = (alpha + 2.0 * y).sqrt();

        let x1 =
            -f_b / (4.0 * f_a) + 0.5 * (w + sqrt(-(alpha + 2.0 * y) - 2.0 * (alpha + beta / w)));
        let x2 =
            -f_b / (4.0 * f_a) + 0.5 * (w - sqrt(-(alpha + 2.0 * y) - 2.0 * (alpha + beta / w)));
        let x3 =
            -f_b / (4.0 * f_a) + 0.5 * (-w + sqrt(-(alpha + 2.0 * y) - 2.0 * (alpha - beta / w)));
        let x4 =
            -f_b / (4.0 * f_a) + 0.5 * (-w - sqrt(-(alpha + 2.0 * y) - 2.0 * (alpha - beta / w)));

        let theta1 = -i * ln(x1);
        let theta2 = -i * ln(x2);
        let theta3 = -i * ln(x3);
        let theta4 = -i * ln(x4);

        vec![theta1, theta2, theta3, theta4]
    };

    fn equals(a: &Complex64, b: &Complex64) -> bool {
        a.approximate_eq(b)
    }
    let deduplicated_solutions = helpers::deduplicate(&solutions, equals);

    let mut result: Vec<IntersectionSolution> = Vec::with_capacity(4);
    for &solution in deduplicated_solutions.iter() {
        if solution.is_real() {
            let c_point =
                Complex64::from(e2.center) + u2 * exp(i * solution) + v2 * exp(-i * solution);
            let point = Point::from(c_point);
            result.push(IntersectionSolution {
                point,
                t1: e1.find_point(point),
                t2: e2.find_point(point),
            });
        }
    }

    result
}

pub fn intersect_lines<T: Into<Line>, U: Into<Line>>(a: T, b: U) -> Vec<IntersectionSolution> {
    let a = a.into();
    let b = b.into();

    let s_a = Straight::from(&a);
    let s_b = Straight::from(&b);

    let candidates = intersect_straights(s_a, s_b);
    candidates
        .into_iter()
        .filter(|solution| a.contains(solution.point) && b.contains(solution.point))
        .collect()
}

pub fn intersect_segments<T: Into<Segment>, U: Into<Segment>>(
    a: T,
    b: U,
) -> Vec<IntersectionSolution> {
    let a = a.into();
    let b = b.into();

    match (a, b) {
        (Segment::Line(first), Segment::Line(second)) => vec![],
        (Segment::Line(first), Segment::Arc(second)) => vec![],
        (Segment::Arc(first), Segment::Line(second)) => vec![],
        (Segment::Arc(first), Segment::Arc(second)) => vec![],
    }
}
