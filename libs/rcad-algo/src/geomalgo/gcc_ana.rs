//! OCCT GccAna_Circ2d3Tan — the analytic Apollonius solver for circles
//! tangent to three circles (GccAna_Circ2d3Tan.cxx, 3-circle constructor
//! L31-1018 + accessors L1022-1210).
//!
//! 1:1 translation. The solver reduces the tangency problem to eight systems
//! of quadratics (R ± Ri tangency per circle), solves them via
//! math_DirectPolynomialRoots and filters the (X, Y, R) triplets.

use glam::DVec2;
use rcad_kernel::geom::Circle2d;
use rcad_kernel::math::direct_polynomial_roots::DirectPolynomialRoots;

use super::geom2d_gcc::GccEntPosition;
use super::geom2d_int::elclib2d;

/// OCCT GccEnt_QualifiedCirc — a circle with a position qualifier.
#[derive(Debug, Clone)]
pub struct QualifiedCirc {
    pub circle: Circle2d,
    pub qualifier: GccEntPosition,
}

impl QualifiedCirc {
    pub fn new(circle: Circle2d, qualifier: GccEntPosition) -> Self {
        QualifiedCirc { circle, qualifier }
    }

    /// OCCT Qualified().
    pub fn qualified(&self) -> Circle2d {
        self.circle
    }

    /// OCCT Qualifier().
    pub fn qualifier(&self) -> GccEntPosition {
        self.qualifier
    }

    pub fn is_unqualified(&self) -> bool {
        self.qualifier == GccEntPosition::Unqualified
    }
    pub fn is_enclosing(&self) -> bool {
        self.qualifier == GccEntPosition::Enclosing
    }
    pub fn is_enclosed(&self) -> bool {
        self.qualifier == GccEntPosition::Enclosed
    }
    pub fn is_outside(&self) -> bool {
        self.qualifier == GccEntPosition::Outside
    }
}

/// OCCT GccAna_Circ2d3Tan — circles tangent to three circles.
#[derive(Debug, Clone)]
pub struct GccAnaCirc2d3Tan {
    pub well_done: bool,
    pub nbr_sol: usize,
    pub cirsol: Vec<Circle2d>,
    qualifier1: Vec<GccEntPosition>,
    qualifier2: Vec<GccEntPosition>,
    qualifier3: Vec<GccEntPosition>,
    the_same1: Vec<bool>,
    the_same2: Vec<bool>,
    the_same3: Vec<bool>,
    pnttg1sol: Vec<DVec2>,
    pnttg2sol: Vec<DVec2>,
    pnttg3sol: Vec<DVec2>,
    par1sol: Vec<f64>,
    par2sol: Vec<f64>,
    par3sol: Vec<f64>,
    pararg1: Vec<f64>,
    pararg2: Vec<f64>,
    pararg3: Vec<f64>,
}

impl GccAnaCirc2d3Tan {
    fn new_empty() -> Self {
        GccAnaCirc2d3Tan {
            well_done: false,
            nbr_sol: 0,
            cirsol: Vec::new(),
            qualifier1: Vec::new(),
            qualifier2: Vec::new(),
            qualifier3: Vec::new(),
            the_same1: Vec::new(),
            the_same2: Vec::new(),
            the_same3: Vec::new(),
            pnttg1sol: Vec::new(),
            pnttg2sol: Vec::new(),
            pnttg3sol: Vec::new(),
            par1sol: Vec::new(),
            par2sol: Vec::new(),
            par3sol: Vec::new(),
            pararg1: Vec::new(),
            pararg2: Vec::new(),
            pararg3: Vec::new(),
        }
    }

    /// OCCT GccAna_Circ2d3Tan(Qualified1, Qualified2, Qualified3, Tolerance)
    /// (GccAna_Circ2d3Tan.cxx L31-1018) — the 3-circle Apollonius solver.
    pub fn new_3_circles(q1: &QualifiedCirc, q2: &QualifiedCirc, q3: &QualifiedCirc, tolerance: f64) -> Self {
        let mut r = GccAnaCirc2d3Tan::new_empty();
        let tol = tolerance.abs();
        r.well_done = false;
        r.nbr_sol = 0;
        if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified())
            || !(q2.is_enclosed() || q2.is_enclosing() || q2.is_outside() || q2.is_unqualified())
            || !(q3.is_enclosed() || q3.is_enclosing() || q3.is_outside() || q3.is_unqualified())
        {
            panic!("GccEnt_BadQualifier");
        }

        let cir1 = q1.qualified();
        let cir2 = q2.qualified();
        let cir3 = q3.qualified();
        let r1 = cir1.radius;
        let r2 = cir2.radius;
        let r3 = cir3.radius;
        let center1 = cir1.center;
        let center2 = cir2.center;
        let center3 = cir3.center;

        let x1 = center1.x;
        let x2 = center2.x;
        let x3 = center3.x;
        let y1 = center1.y;
        let y2 = center2.y;
        let y3 = center3.y;

        let dir2 = center1 - center2;
        let dir3 = center1 - center3;

        // OCCT L99-146: degenerate configurations (concentric, collinear
        // tangent triples) return with no solutions.
        if ((r1 - r2).abs() <= tolerance && center1.distance(center2) <= tolerance)
            || ((r1 - r3).abs() <= tolerance && center1.distance(center3) <= tolerance)
            || ((r2 - r3).abs() <= tolerance && center2.distance(center3) <= tolerance)
        {
            return r;
        } else {
            let cross = dir2.x * dir3.y - dir2.y * dir3.x;
            if cross.abs() <= tolerance {
                let dist1 = center1.distance(center2);
                let dist2 = center1.distance(center3);
                let dist3 = center2.distance(center3);
                if ((r1 - r2).abs() - dist1).abs() <= tolerance {
                    if ((r1 - r3).abs() - dist2).abs() <= tolerance {
                        if ((r2 - r3).abs() - dist3).abs() <= tolerance {
                            return r;
                        }
                    } else if (r1 + r3 - dist2).abs() <= tolerance {
                        if (r2 + r3 - dist3).abs() <= tolerance {
                            return r;
                        }
                    }
                } else if (r1 + r2 - dist1).abs() <= tolerance {
                    if ((r1 - r3).abs() - dist2).abs() <= tolerance
                        && (r2 + r3 - dist3).abs() <= tolerance
                    {
                        // No return here (OCCT L132-135 empty branch).
                    } else if ((r2 - r3).abs() - dist3).abs() <= tolerance
                        && (r1 + r3 - dist2).abs() <= tolerance
                    {
                        return r;
                    }
                }
            }
        }

        // OCCT L148-172: coefficient arrays (1..=8 systems).
        let mut a2 = [0.0; 9];
        let mut b2 = [0.0; 9];
        let mut c2 = [0.0; 9];
        let mut d2 = [0.0; 9];
        let mut e2 = [0.0; 9];
        let mut f2 = [0.0; 9];
        let mut a3 = [0.0; 9];
        let mut b3 = [0.0; 9];
        let mut c3 = [0.0; 9];
        let mut d3 = [0.0; 9];
        let mut e3 = [0.0; 9];
        let mut f3 = [0.0; 9];
        let mut beta2 = [0.0; 9];
        let mut gamma2 = [0.0; 9];
        let mut delta2 = [0.0; 9];
        let mut beta3 = [0.0; 9];
        let mut gamma3 = [0.0; 9];
        let mut delta3 = [0.0; 9];

        let is_touch = ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) - (r1 - r2) * (r1 - r2))
            .abs()
            <= tolerance
            || ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) - (r1 + r2) * (r1 + r2)).abs()
                <= tolerance
            || ((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) - (r1 - r3) * (r1 - r3)).abs()
                <= tolerance
            || ((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3) - (r1 + r3) * (r1 + r3)).abs()
                <= tolerance
            || ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) - (r2 - r3) * (r2 - r3)).abs()
                <= tolerance
            || ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) - (r2 + r3) * (r2 + r3)).abs()
                <= tolerance;

        // First step (OCCT L229-312): Beta/Gamma/Delta and conic coefficients.
        for i in 1..=8 {
            if i == 1 || i == 4 || i == 5 || i == 8 {
                if (r1 - r2).abs() > tolerance {
                    beta2[i] = (x1 - x2) / (r1 - r2);
                    gamma2[i] = (y1 - y2) / (r1 - r2);
                    delta2[i] =
                        (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + (r1 - r2) * (r1 - r2))
                            / (2.0 * (r1 - r2));
                }
            } else {
                beta2[i] = (x1 - x2) / (r1 + r2);
                gamma2[i] = (y1 - y2) / (r1 + r2);
                delta2[i] =
                    (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + (r1 + r2) * (r1 + r2))
                        / (2.0 * (r1 + r2));
            }
            if (i == 1 || i == 4 || i == 5 || i == 8) && (r1 - r2).abs() <= tolerance {
                // If R1 = R2.
                a2[i] = 0.0;
                b2[i] = 0.0;
                c2[i] = 0.0;
                d2[i] = x2 - x1;
                e2[i] = y2 - y1;
                f2[i] = x1 * x1 - x2 * x2 + y1 * y1 - y2 * y2;
            } else {
                a2[i] = beta2[i] * beta2[i] - 1.0;
                b2[i] = beta2[i] * gamma2[i];
                c2[i] = gamma2[i] * gamma2[i] - 1.0;
                d2[i] = beta2[i] * delta2[i] + x2;
                e2[i] = gamma2[i] * delta2[i] + y2;
                f2[i] = delta2[i] * delta2[i] - x2 * x2 - y2 * y2;
            }

            if i == 1 || i == 3 || i == 6 || i == 8 {
                if (r1 - r3).abs() > tolerance {
                    beta3[i] = (x1 - x3) / (r1 - r3);
                    gamma3[i] = (y1 - y3) / (r1 - r3);
                    delta3[i] =
                        (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + (r1 - r3) * (r1 - r3))
                            / (2.0 * (r1 - r3));
                }
            } else {
                beta3[i] = (x1 - x3) / (r1 + r3);
                gamma3[i] = (y1 - y3) / (r1 + r3);
                delta3[i] =
                    (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1 + (r1 + r3) * (r1 + r3))
                        / (2.0 * (r1 + r3));
            }
            if (i == 1 || i == 3 || i == 6 || i == 8) && (r1 - r3).abs() <= tolerance {
                a3[i] = 0.0;
                b3[i] = 0.0;
                c3[i] = 0.0;
                d3[i] = x3 - x1;
                e3[i] = y3 - y1;
                f3[i] = x1 * x1 - x3 * x3 + y1 * y1 - y3 * y3;
            } else {
                a3[i] = beta3[i] * beta3[i] - 1.0;
                b3[i] = beta3[i] * gamma3[i];
                c3[i] = gamma3[i] * gamma3[i] - 1.0;
                d3[i] = beta3[i] * delta3[i] + x3;
                e3[i] = gamma3[i] * delta3[i] + y3;
                f3[i] = delta3[i] * delta3[i] - x3 * x3 - y3 * y3;
            }
        }

        // Second step (OCCT L321-672): solve the 8 systems for (X, Y).
        let mut first_sol = [0usize; 10];
        let mut xs: Vec<f64> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        let mut cur_sol = 0usize;
        for i in 1..=8 {
            let (a2v, a3v, b2v, b3v, c2v, c3v, d2v, d3v, e2v, e3v, f2v, f3v) = (
                a2[i], a3[i], b2[i], b3[i], c2[i], c3[i], d2[i], d3[i], e2[i], e3[i], f2[i], f3[i],
            );
            first_sol[i] = cur_sol + 1;

            // OCCT L339-359: systems with no solutions due to qualifiers.
            if ((i == 2 || i == 5 || i == 6 || i == 8)
                && (q1.is_enclosed() || q1.is_enclosing()))
                || ((i == 1 || i == 3 || i == 4 || i == 7) && q1.is_outside())
            {
                continue;
            }
            if ((i == 3 || i == 5 || i == 7 || i == 8)
                && (q2.is_enclosed() || q2.is_enclosing()))
                || ((i == 1 || i == 2 || i == 4 || i == 6) && q2.is_outside())
            {
                continue;
            }
            if ((i == 4 || i == 6 || i == 7 || i == 8)
                && (q3.is_enclosed() || q3.is_enclosing()))
                || ((i == 1 || i == 2 || i == 3 || i == 5) && q3.is_outside())
            {
                continue;
            }

            // OCCT L363-371: Cir1 itself is a solution of this system.
            if (a2v - a3v).abs() <= tolerance
                && (b2v - b3v).abs() <= tolerance
                && (c2v - c3v).abs() <= tolerance
                && (d2v - d3v).abs() <= tolerance
                && (e2v - e3v).abs() <= tolerance
                && (f2v - f3v).abs() <= tolerance
            {
                xs.push(x1);
                ys.push(y1);
                cur_sol += 1;
                continue;
            }

            // 1) a2 = 0 (OCCT L372-488).
            if a2v.abs() <= tolerance {
                // 1.1) b2y + d2 = 0 — quadratic in y.
                let mut y_roots = DirectPolynomialRoots::new_quadratic(c2v, 2.0 * e2v, f2v);
                if y_roots.is_done() && !y_roots.infinite_roots() {
                    for k in 1..=y_roots.nb_solutions() {
                        let y = y_roots.value(k);
                        if (k != 2 || (y - y_roots.value(1)).abs() > 10.0 * tolerance)
                            && (b2v * y + d2v).abs() <= b2v * tolerance
                        {
                            let mut x_roots = DirectPolynomialRoots::new_quadratic(
                                a3v,
                                2.0 * (b3v * y + d3v),
                                c3v * (y * y) + 2.0 * e3v * y + f3v,
                            );
                            if x_roots.is_done() && !x_roots.infinite_roots() {
                                for j in 1..=x_roots.nb_solutions() {
                                    let x = x_roots.value(j);
                                    if j != 2 || (x - x_roots.value(1)).abs() > 10.0 * tolerance {
                                        xs.push(x);
                                        ys.push(y);
                                        cur_sol += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                // 1.2) b2y + d2 != 0 — quartic in y.
                let a = a3v * c2v * c2v - 4.0 * b2v * (b3v * c2v - b2v * c3v);
                let b = 4.0 * a3v * c2v * e2v - 4.0 * b3v * (c2v * d2v + 2.0 * b2v * e2v)
                    + 4.0 * b2v * (2.0 * c3v * d2v - c2v * d3v + 2.0 * b2v * e3v);
                let c = 2.0 * a3v * (c2v * f2v + 2.0 * e2v * e2v)
                    - 4.0 * b3v * (b2v * f2v + 2.0 * e2v * d2v) + 4.0 * c3v * d2v * d2v
                    - 4.0 * d3v * (c2v * d2v + 2.0 * b2v * e2v) + 16.0 * b2v * e3v * d2v
                    + 4.0 * b2v * b2v * f3v;
                let d = 4.0 * a3v * e2v * f2v - 4.0 * b3v * d2v * f2v
                    - 4.0 * d3v * (b2v * f2v + 2.0 * d2v * e2v) + 8.0 * d2v * d2v * e3v
                    + 8.0 * b2v * d2v * f3v;
                let e = a3v * f2v * f2v - 4.0 * d2v * d3v * f2v + 4.0 * d2v * d2v * f3v;

                // Special case: one circle touches another (derivative quartic).
                if is_touch {
                    let mut y_roots1 = DirectPolynomialRoots::new_cubic(4.0 * a, 3.0 * b, 2.0 * c, d);
                    if y_roots1.is_done() && !y_roots1.infinite_roots() {
                        for k in 1..=y_roots1.nb_solutions() {
                            let y = y_roots1.value(k);
                            let mut is_same = false;
                            for l in 1..k {
                                if (y - y_roots1.value(l)).abs() <= 10.0 * tolerance {
                                    is_same = true;
                                }
                            }
                            let eps = (((4.0 * a * y).abs() + (3.0 * b).abs()) * y).abs()
                                + (2.0 * c).abs();
                            let eps = (eps * y).abs() + d.abs();
                            if ((((a * y + b) * y + c) * y + d) * y + e).abs() <= eps * tolerance {
                                if !is_same && (b2v * y + d2v).abs() > b2v * tolerance {
                                    let x = -(c2v * (y * y) + 2.0 * e2v * y + f2v)
                                        / (2.0 * (b2v * y + d2v));
                                    xs.push(x);
                                    ys.push(y);
                                    cur_sol += 1;
                                }
                            }
                        }
                    }
                }

                let mut y_roots1 = DirectPolynomialRoots::new_quartic(a, b, c, d, e);
                if y_roots1.is_done() && !y_roots1.infinite_roots() {
                    for k in 1..=y_roots1.nb_solutions() {
                        let y = y_roots1.value(k);
                        let mut is_same = false;
                        let first_index = if i == 1 { 1 } else { first_sol[i] };
                        for l in (first_index - 1)..cur_sol {
                            if (y - ys[l]).abs() <= 10.0 * tolerance {
                                is_same = true;
                            }
                        }
                        if !is_same && (b2v * y + d2v).abs() > b2v * tolerance {
                            let x = -(c2v * (y * y) + 2.0 * e2v * y + f2v)
                                / (2.0 * (b2v * y + d2v));
                            xs.push(x);
                            ys.push(y);
                            cur_sol += 1;
                        }
                    }
                }
            } else {
                // 2) a2 != 0 (OCCT L490-669).
                let m = 2.0 * a3v * b2v * b2v / (a2v * a2v) - 2.0 * b2v * b3v / a2v
                    - a3v * c2v / a2v + c3v;
                let n = 4.0 * a3v * b2v * d2v / (a2v * a2v) - 2.0 * b3v * d2v / a2v
                    - 2.0 * b2v * d3v / a2v - 2.0 * a3v * e2v / a2v + 2.0 * e3v;
                let t = 2.0 * a3v * d2v * d2v / (a2v * a2v) - 2.0 * d2v * d3v / a2v
                    - a3v * f2v / a2v + f3v;
                let s = 2.0 * b3v - 2.0 * a3v * b2v / a2v;
                let v = 2.0 * d3v - 2.0 * d2v * a3v / a2v;

                // If s = v = 0 (OCCT L501-543).
                if s.abs() <= tolerance && v.abs() <= tolerance {
                    let mut y_roots = DirectPolynomialRoots::new_quadratic(m, n, t);
                    if y_roots.is_done() && !y_roots.infinite_roots() {
                        for k in 1..=y_roots.nb_solutions() {
                            let y = y_roots.value(k);
                            let p = -(b2v * y + d2v) / a2v;
                            let q = (c2v * (y * y) + 2.0 * e2v * y + f2v) / a2v;
                            let eps = 2.0
                                * (((b2v * b2v + (a2v * c2v).abs()) * y).abs() + (b2v * d2v).abs()
                                    + (a2v * e2v).abs())
                                / (a2v * a2v);
                            if (k != 2 || (y - y_roots.value(1)).abs() > 10.0 * tolerance)
                                && p * p - q >= -eps * tolerance
                            {
                                let mut x_roots = DirectPolynomialRoots::new_quadratic(
                                    a2v,
                                    2.0 * (b2v * y + d2v),
                                    c2v * y * y + 2.0 * e2v * y + f2v,
                                );
                                if x_roots.is_done() && !x_roots.infinite_roots() {
                                    for l in 1..=x_roots.nb_solutions() {
                                        let x = x_roots.value(l);
                                        if l != 2
                                            || (x - x_roots.value(1)).abs() > 10.0 * tolerance
                                        {
                                            xs.push(x);
                                            ys.push(y);
                                            cur_sol += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // (s*y + v) != 0 (OCCT L549-668).
                    let a = s * s * (b2v * b2v - a2v * c2v) - m * m * a2v * a2v;
                    let b = 2.0 * s * v * (b2v * b2v - a2v * c2v)
                        + 2.0 * s * s * (b2v * d2v - a2v * e2v) - 2.0 * m * n * a2v * a2v;
                    let c = v * v * (b2v * b2v - a2v * c2v) + 4.0 * s * v * (b2v * d2v - a2v * e2v)
                        + s * s * (d2v * d2v - a2v * f2v) - a2v * a2v * (2.0 * m * t + n * n);
                    let d = 2.0 * v * v * (b2v * d2v - a2v * e2v)
                        + 2.0 * s * v * (d2v * d2v - a2v * f2v) - 2.0 * n * t * a2v * a2v;
                    let e = v * v * (d2v * d2v - a2v * f2v) - t * t * a2v * a2v;

                    if is_touch {
                        let mut y_roots1 = DirectPolynomialRoots::new_cubic(4.0 * a, 3.0 * b, 2.0 * c, d);
                        if y_roots1.is_done() && !y_roots1.infinite_roots() {
                            for k in 1..=y_roots1.nb_solutions() {
                                let y = y_roots1.value(k);
                                let p = -(b2v * y + d2v) / a2v;
                                let q = (c2v * (y * y) + 2.0 * e2v * y + f2v) / a2v;
                                let mut is_same = false;
                                let first_index = if i == 1 { 1 } else { first_sol[i] };
                                for l in (first_index - 1)..cur_sol {
                                    if (y - ys[l]).abs() <= 10.0 * tolerance {
                                        is_same = true;
                                    }
                                }
                                let eps = (((4.0 * a * y).abs() + (3.0 * b).abs()) * y).abs()
                                    + (2.0 * c).abs();
                                let eps = (eps * y).abs() + d.abs();
                                if ((((a * y + b) * y + c) * y + d) * y + e).abs()
                                    <= eps * tolerance
                                {
                                    let eps = 2.0
                                        * (((b2v * b2v + (a2v * c2v).abs()) * y).abs()
                                            + (b2v * d2v).abs()
                                            + (a2v * e2v).abs())
                                        / (a2v * a2v);
                                    if !is_same && p * p - q >= -eps * tolerance {
                                        let mut x_roots = DirectPolynomialRoots::new_quadratic(
                                            a2v,
                                            2.0 * (b2v * y + d2v),
                                            c2v * y * y + 2.0 * e2v * y + f2v,
                                        );
                                        if x_roots.is_done() && !x_roots.infinite_roots() {
                                            for l in 1..=x_roots.nb_solutions() {
                                                let x = x_roots.value(l);
                                                if l != 2
                                                    || (x - x_roots.value(1)).abs()
                                                        > 10.0 * tolerance
                                                {
                                                    xs.push(x);
                                                    ys.push(y);
                                                    cur_sol += 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let mut y_roots = DirectPolynomialRoots::new_quartic(a, b, c, d, e);
                    if y_roots.is_done() && !y_roots.infinite_roots() {
                        for k in 1..=y_roots.nb_solutions() {
                            let y = y_roots.value(k);
                            let p = -(b2v * y + d2v) / a2v;
                            let q = (c2v * (y * y) + 2.0 * e2v * y + f2v) / a2v;
                            let mut is_same = false;
                            for l in 1..k {
                                if (y - y_roots.value(l)).abs() <= 10.0 * tolerance {
                                    is_same = true;
                                }
                            }
                            // OCCT L631-639: also check against previously added
                            // solutions from the is_touch path for this system.
                            let first_index = if i == 1 { 1 } else { first_sol[i] };
                            for l in (first_index - 1)..cur_sol {
                                if (y - ys[l]).abs() <= 10.0 * tolerance {
                                    is_same = true;
                                }
                            }
                            let eps = 2.0
                                * (((b2v * b2v + (a2v * c2v).abs()) * y).abs() + (b2v * d2v).abs()
                                    + (a2v * e2v).abs())
                                / (a2v * a2v);
                            if !is_same && p * p - q >= -eps * tolerance {
                                let mut x_roots = DirectPolynomialRoots::new_quadratic(
                                    a2v,
                                    2.0 * (b2v * y + d2v),
                                    c2v * y * y + 2.0 * e2v * y + f2v,
                                );
                                if x_roots.is_done() && !x_roots.infinite_roots() {
                                    for l in 1..=x_roots.nb_solutions() {
                                        let x = x_roots.value(l);
                                        if l != 2 || (x - x_roots.value(1)).abs() > 10.0 * tolerance {
                                            xs.push(x);
                                            ys.push(y);
                                            cur_sol += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        first_sol[9] = cur_sol + 1;

        // Third step (OCCT L676-774): compute R for each couple (X, Y).
        let mut first_sol1 = [0usize; 10];
        let mut xs1: Vec<f64> = Vec::new();
        let mut ys1: Vec<f64> = Vec::new();
        let mut rs1: Vec<f64> = Vec::new();
        cur_sol = 0;
        for i in 1..=8 {
            first_sol1[i] = cur_sol + 1;
            for j in (first_sol[i] - 1)..(first_sol[i + 1] - 1) {
                let x = xs[j];
                let y = ys[j];
                // In some cases when R1 = R2.
                if (i == 1 || i == 4 || i == 5 || i == 8) && (r1 - r2).abs() <= tolerance {
                    if i == 1 || i == 4 {
                        let mut r = r1 + ((x - x1) * (x - x1) + (y - y1) * (y - y1)).sqrt();
                        let eps = 10.0 * (2.0 * (r - r2).abs() + (x - x2).abs() + (y - y2).abs());
                        if ((r - r2) * (r - r2) - (x - x2) * (x - x2) - (y - y2) * (y - y2)).abs()
                            <= eps * tolerance
                        {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                        r = r1 - ((x - x1) * (x - x1) + (y - y1) * (y - y1)).sqrt();
                        let eps = 10.0 * (2.0 * (r - r2).abs() + (x - x2).abs() + (y - y2).abs());
                        if r > tolerance
                            && ((r - r2) * (r - r2) - (x - x2) * (x - x2) - (y - y2) * (y - y2))
                                .abs()
                                <= eps * tolerance
                        {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                    } else {
                        // i == 5 || i == 8.
                        let r = -r1 + ((x - x1) * (x - x1) + (y - y1) * (y - y1)).sqrt();
                        if r > tolerance {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                    }
                } else {
                    // Other cases.
                    if i == 1 || i == 4 {
                        let r = r2 + beta2[i] * x + gamma2[i] * y + delta2[i];
                        if r > tolerance {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                    }
                    if i == 5 || i == 8 {
                        let r = -r2 - beta2[i] * x - gamma2[i] * y - delta2[i];
                        if r > tolerance {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                    }
                    if i == 3 || i == 7 {
                        let r = -r2 + beta2[i] * x + gamma2[i] * y + delta2[i];
                        if r > tolerance {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                    }
                    if i == 2 || i == 6 {
                        let r = r2 - beta2[i] * x - gamma2[i] * y - delta2[i];
                        if r > tolerance {
                            xs1.push(x);
                            ys1.push(y);
                            rs1.push(r);
                            cur_sol += 1;
                        }
                    }
                }
            }
        }
        first_sol1[9] = cur_sol + 1;

        // Fourth step (OCCT L777-862): verify the triplets (X, Y, R).
        xs.clear();
        ys.clear();
        let mut rs: Vec<f64> = Vec::new();
        cur_sol = 0;
        for i in 1..=8 {
            first_sol[i] = cur_sol + 1;
            for j in (first_sol1[i] - 1)..(first_sol1[i + 1] - 1) {
                let x = xs1[j];
                let y = ys1[j];
                let r = rs1[j];
                // In some cases when R1 = R3.
                if (i == 1 || i == 3 || i == 6 || i == 8) && (r1 - r3).abs() <= tolerance {
                    if i == 1 || i == 3 {
                        let eps = 10.0 * (2.0 * (r - r3).abs() + (x - x3).abs() + (y - y3).abs());
                        if ((r - r3) * (r - r3) - (x - x3) * (x - x3) - (y - y3) * (y - y3)).abs()
                            <= eps * tolerance
                        {
                            xs.push(x);
                            ys.push(y);
                            rs.push(r);
                            cur_sol += 1;
                        }
                    } else {
                        // i == 6 || i == 8.
                        let eps = 10.0 * (2.0 * (r + r3) + (x - x3).abs() + (y - y3).abs());
                        if ((r + r3) * (r + r3) - (x - x3) * (x - x3) - (y - y3) * (y - y3)).abs()
                            <= eps * tolerance
                        {
                            xs.push(x);
                            ys.push(y);
                            rs.push(r);
                            cur_sol += 1;
                        }
                    }
                } else {
                    // Other cases.
                    let eps = 10.0 * (beta3[i].abs() + gamma3[i].abs() + 1.0);
                    if i == 1 || i == 3 {
                        if (r3 + beta3[i] * x + gamma3[i] * y + delta3[i] - r).abs()
                            <= eps * tolerance
                        {
                            xs.push(x);
                            ys.push(y);
                            rs.push(r);
                            cur_sol += 1;
                        }
                    }
                    if i == 6 || i == 8 {
                        if (r3 + beta3[i] * x + gamma3[i] * y + delta3[i] + r).abs()
                            <= eps * tolerance
                        {
                            xs.push(x);
                            ys.push(y);
                            rs.push(r);
                            cur_sol += 1;
                        }
                    }
                    if i == 4 || i == 7 {
                        if (beta3[i] * x + gamma3[i] * y + delta3[i] - r - r3).abs()
                            <= eps * tolerance
                        {
                            xs.push(x);
                            ys.push(y);
                            rs.push(r);
                            cur_sol += 1;
                        }
                    }
                    if i == 2 || i == 5 {
                        if (r - r3 + beta3[i] * x + gamma3[i] * y + delta3[i]).abs()
                            <= eps * tolerance
                        {
                            xs.push(x);
                            ys.push(y);
                            rs.push(r);
                            cur_sol += 1;
                        }
                    }
                }
            }
        }
        first_sol[9] = cur_sol + 1;

        // Fifth step (OCCT L866-1016): qualifier filter and solution building.
        for i in 1..=8 {
            for j in (first_sol[i] - 1)..(first_sol[i + 1] - 1) {
                if (q1.is_enclosed() && rs[j] > r1) || (q1.is_enclosing() && rs[j] < r1) {
                    continue;
                }
                if (q2.is_enclosed() && rs[j] > r2) || (q2.is_enclosing() && rs[j] < r2) {
                    continue;
                }
                if (q3.is_enclosed() && rs[j] > r3) || (q3.is_enclosing() && rs[j] < r3) {
                    continue;
                }

                r.nbr_sol += 1;
                // RLE, avoid out of range.
                if r.nbr_sol > 16 {
                    r.nbr_sol = 16;
                }

                let center = DVec2::new(xs[j], ys[j]);
                let circle = Circle2d {
                    center,
                    x_dir: DVec2::X,
                    y_dir: DVec2::Y,
                    radius: rs[j],
                };
                r.cirsol.push(circle);

                // OCCT L897-949: qualifiers (note: OCCT reuses distcc1 for the
                // second and third circles — ported as-is).
                let distcc1 = center.distance(center1);
                if !q1.is_unqualified() {
                    r.qualifier1.push(q1.qualifier());
                } else if (distcc1 + rs[j] - r1).abs() <= tol {
                    r.qualifier1.push(GccEntPosition::Enclosed);
                } else if (distcc1 - r1 - rs[j]).abs() <= tol {
                    r.qualifier1.push(GccEntPosition::Outside);
                } else {
                    r.qualifier1.push(GccEntPosition::Enclosing);
                }
                let distcc2 = center.distance(center1);
                if !q2.is_unqualified() {
                    r.qualifier2.push(q2.qualifier());
                } else if (distcc2 + rs[j] - r2).abs() <= tol {
                    r.qualifier2.push(GccEntPosition::Enclosed);
                } else if (distcc2 - r2 - rs[j]).abs() <= tol {
                    r.qualifier2.push(GccEntPosition::Outside);
                } else {
                    r.qualifier2.push(GccEntPosition::Enclosing);
                }
                let distcc3 = center.distance(center1);
                if !q3.is_unqualified() {
                    r.qualifier3.push(q3.qualifier());
                } else if (distcc3 + rs[j] - r3).abs() <= tol {
                    r.qualifier3.push(GccEntPosition::Enclosed);
                } else if (distcc3 - r3 - rs[j]).abs() <= tol {
                    r.qualifier3.push(GccEntPosition::Outside);
                } else {
                    r.qualifier3.push(GccEntPosition::Enclosing);
                }

                // OCCT L953-1014: TheSame flags + tangency points + parameters
                // (OCCT arrays are index-aligned; the tangency slot keeps its
                // index even when TheSame — placeholder here).
                if center.distance(cir1.center) <= tolerance {
                    r.the_same1.push(true);
                    r.pnttg1sol.push(DVec2::ZERO);
                    r.par1sol.push(0.0);
                    r.pararg1.push(0.0);
                } else {
                    r.the_same1.push(false);
                    let dc = if (i == 2 || i == 5 || i == 6 || i == 8) || rs[j] > r1 {
                        (cir1.center - center).normalize_or_zero()
                    } else {
                        (center - cir1.center).normalize_or_zero()
                    };
                    let pnt = center + rs[j] * dc;
                    r.pnttg1sol.push(pnt);
                    r.par1sol.push(elclib2d::circle_parameter(
                        circle.center,
                        circle.x_dir,
                        circle.y_dir,
                        pnt,
                    ));
                    r.pararg1.push(elclib2d::circle_parameter(
                        cir1.center,
                        cir1.x_dir,
                        cir1.y_dir,
                        pnt,
                    ));
                }

                if center.distance(cir2.center) <= tolerance {
                    r.the_same2.push(true);
                    r.pnttg2sol.push(DVec2::ZERO);
                    r.par2sol.push(0.0);
                    r.pararg2.push(0.0);
                } else {
                    r.the_same2.push(false);
                    let dc = if (i == 3 || i == 5 || i == 7 || i == 8) || rs[j] > r2 {
                        (cir2.center - center).normalize_or_zero()
                    } else {
                        (center - cir2.center).normalize_or_zero()
                    };
                    let pnt = center + rs[j] * dc;
                    r.pnttg2sol.push(pnt);
                    r.par2sol.push(elclib2d::circle_parameter(
                        circle.center,
                        circle.x_dir,
                        circle.y_dir,
                        pnt,
                    ));
                    r.pararg2.push(elclib2d::circle_parameter(
                        cir2.center,
                        cir2.x_dir,
                        cir2.y_dir,
                        pnt,
                    ));
                }

                if center.distance(cir3.center) <= tolerance {
                    r.the_same3.push(true);
                    r.pnttg3sol.push(DVec2::ZERO);
                    r.par3sol.push(0.0);
                    r.pararg3.push(0.0);
                } else {
                    r.the_same3.push(false);
                    let dc = if (i == 4 || i == 6 || i == 7 || i == 8) || rs[j] > r3 {
                        (cir3.center - center).normalize_or_zero()
                    } else {
                        (center - cir3.center).normalize_or_zero()
                    };
                    let pnt = center + rs[j] * dc;
                    r.pnttg3sol.push(pnt);
                    r.par3sol.push(elclib2d::circle_parameter(
                        circle.center,
                        circle.x_dir,
                        circle.y_dir,
                        pnt,
                    ));
                    r.pararg3.push(elclib2d::circle_parameter(
                        cir3.center,
                        cir3.x_dir,
                        cir3.y_dir,
                        pnt,
                    ));
                }
            }
        }
        r.well_done = true;
        r
    }

    /// OCCT IsDone() (L1022-1025).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// OCCT NbSolutions() (L1027-1030).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// OCCT ThisSolution(Index) (L1032-1045).
    pub fn this_solution(&self, index: usize) -> Circle2d {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        self.cirsol[index - 1]
    }

    /// OCCT WhichQualifier(Index, Qualif1, Qualif2, Qualif3) (L1047-1066).
    pub fn which_qualifier(&self, index: usize) -> (GccEntPosition, GccEntPosition, GccEntPosition) {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        (
            self.qualifier1[index - 1],
            self.qualifier2[index - 1],
            self.qualifier3[index - 1],
        )
    }

    /// OCCT Tangency1 (L1068-1094).
    pub fn tangency1(&self, index: usize) -> (f64, f64, DVec2) {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        if !self.the_same1[index - 1] {
            (
                self.par1sol[index - 1],
                self.pararg1[index - 1],
                self.pnttg1sol[index - 1],
            )
        } else {
            panic!("StdFail_NotDone");
        }
    }

    /// OCCT Tangency2 (L1096-1122).
    pub fn tangency2(&self, index: usize) -> (f64, f64, DVec2) {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        if !self.the_same2[index - 1] {
            (
                self.par2sol[index - 1],
                self.pararg2[index - 1],
                self.pnttg2sol[index - 1],
            )
        } else {
            panic!("StdFail_NotDone");
        }
    }

    /// OCCT Tangency3 (L1124-1150).
    pub fn tangency3(&self, index: usize) -> (f64, f64, DVec2) {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        if !self.the_same3[index - 1] {
            (
                self.par3sol[index - 1],
                self.pararg3[index - 1],
                self.pnttg3sol[index - 1],
            )
        } else {
            panic!("StdFail_NotDone");
        }
    }

    /// OCCT IsTheSame1 (L1152-1170).
    pub fn is_the_same1(&self, index: usize) -> bool {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        self.the_same1[index - 1]
    }

    /// OCCT IsTheSame2 (L1172-1190).
    pub fn is_the_same2(&self, index: usize) -> bool {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        self.the_same2[index - 1]
    }

    /// OCCT IsTheSame3 (L1192-1210).
    pub fn is_the_same3(&self, index: usize) -> bool {
        assert!(self.well_done, "StdFail_NotDone");
        assert!(index >= 1 && index <= self.nbr_sol, "Standard_OutOfRange");
        self.the_same3[index - 1]
    }
}
