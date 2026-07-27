//! ClassifyLin2d (IntTools_FaceFace.cxx L2574-2640).
//! Clips a 2D line to a UV rectangle, returning the parameter range where
//! the line passes through the rectangle.

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};

/// Line2d parameter from point: t such that P = O + t·D (D is unit).
fn line2d_param(l: &rcad_kernel::geom::Line2d, p: DVec2) -> f64 {
    (p - l.origin).dot(l.direction)
}

/// OCCT L2574-2640: ClassifyLin2d — find parameter range of a Line2d
/// intersecting the UV rectangle [xmin,xmax]×[ymin,ymax].
/// Returns Some([t_min, t_max]) if the line passes through, None if it misses.
pub fn classify_lin2d(pc: &Curve2d, uv: [f64; 4], tol: f64) -> Option<[f64; 2]> {
    let (xmin, xmax, ymin, ymax) = (uv[0], uv[1], uv[2], uv[3]);
    let (A, B, C) = match pc {
        Curve2d::Line(l) => {
            // Implicit form: A·x + B·y + C = 0
            // where A = D.y, B = -D.x, C = -(D.y·O.x - D.x·O.y)
            (
                l.direction.y,
                -l.direction.x,
                -(l.direction.y * l.origin.x - l.direction.x * l.origin.y),
            )
        }
        _ => return None,
    };

    let coinc = |a: f64, b: f64| a.abs() <= tol && b.abs() <= tol;
    // INTER: signs differ (one on each side of boundary).
    // Values within tol are treated as zero (on the boundary).
    // 0.0 signum is treated as compatible with any sign (line at corner).
    let inter = |a: f64, b: f64| {
        let sa = if a.abs() < tol { 0.0 } else { a.signum() };
        let sb = if b.abs() < tol { 0.0 } else { b.signum() };
        sa != sb
    };

    fn get_line<'a>(pc: &'a Curve2d) -> &'a rcad_kernel::geom::Line2d {
        match pc {
            Curve2d::Line(l) => l,
            _ => unreachable!(),
        }
    }
    let l = get_line(pc);

    let mut par: Vec<f64> = Vec::with_capacity(2);

    // edge x = xmin, y in [ymin, ymax]
    let d1 = A * xmin + B * ymin + C;
    let d2 = A * xmin + B * ymax + C;
    if inter(d1, d2) && B.abs() > 1e-15 {
        let y = -(C + A * xmin) / B;
        if y >= ymin - tol && y <= ymax + tol {
            par.push(line2d_param(l, DVec2::new(xmin, y)));
        }
    } else if coinc(d1, d2) {
        par.push(line2d_param(l, DVec2::new(xmin, ymin)));
        par.push(line2d_param(l, DVec2::new(xmin, ymax)));
    }
    if par.len() >= 2 {
        return Some([par[0].min(par[1]), par[0].max(par[1])]);
    }

    // edge y = ymax, x in [xmin, xmax]
    let d1 = A * xmin + B * ymax + C;
    let d2 = A * xmax + B * ymax + C;
    if inter(d1, d2) && A.abs() > 1e-15 {
        let x = -(C + B * ymax) / A;
        if x >= xmin - tol && x <= xmax + tol {
            par.push(line2d_param(l, DVec2::new(x, ymax)));
        }
    } else if coinc(d1, d2) && par.is_empty() {
        par.push(line2d_param(l, DVec2::new(xmin, ymax)));
        par.push(line2d_param(l, DVec2::new(xmax, ymax)));
    }
    if par.len() >= 2 {
        return Some([par[0].min(par[1]), par[0].max(par[1])]);
    }

    // edge x = xmax, y in [ymin, ymax]
    let d1 = A * xmax + B * ymax + C;
    let d2 = A * xmax + B * ymin + C;
    if inter(d1, d2) && B.abs() > 1e-15 {
        let y = -(C + A * xmax) / B;
        if y >= ymin - tol && y <= ymax + tol {
            par.push(line2d_param(l, DVec2::new(xmax, y)));
        }
    } else if coinc(d1, d2) && par.is_empty() {
        par.push(line2d_param(l, DVec2::new(xmax, ymin)));
        par.push(line2d_param(l, DVec2::new(xmax, ymax)));
    }
    if par.len() >= 2 {
        return Some([par[0].min(par[1]), par[0].max(par[1])]);
    }

    // edge y = ymin, x in [xmin, xmax]
    let d1 = A * xmax + B * ymin + C;
    let d2 = A * xmin + B * ymin + C;
    if inter(d1, d2) && A.abs() > 1e-15 {
        let x = -(C + B * ymin) / A;
        if x >= xmin - tol && x <= xmax + tol {
            par.push(line2d_param(l, DVec2::new(x, ymin)));
        }
    } else if coinc(d1, d2) && par.is_empty() {
        par.push(line2d_param(l, DVec2::new(xmin, ymin)));
        par.push(line2d_param(l, DVec2::new(xmax, ymin)));
    }
    if par.len() >= 2 {
        Some([par[0].min(par[1]), par[0].max(par[1])])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::classify_lin2d;
    use glam::DVec2;
    use rcad_kernel::geom::{Curve2d, Line2d};

    fn line(p: DVec2, d: DVec2) -> Curve2d {
        Curve2d::Line(Line2d::new(p, d))
    }

    #[test]
    fn horizontal_inside() {
        let r = classify_lin2d(
            &line(DVec2::new(0.0, 0.5), DVec2::new(1.0, 0.0)),
            [0.0, 1.0, 0.0, 1.0],
            1e-7,
        );
        assert!(r.is_some());
        assert!((r.unwrap()[1] - r.unwrap()[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn vertical_inside() {
        let r = classify_lin2d(
            &line(DVec2::new(0.3, 0.0), DVec2::new(0.0, 1.0)),
            [0.0, 1.0, 0.0, 1.0],
            1e-7,
        );
        assert!(r.is_some());
    }

    #[test]
    fn misses_rectangle() {
        let r = classify_lin2d(
            &line(DVec2::new(0.0, 2.0), DVec2::new(1.0, 0.0)),
            [0.0, 1.0, 0.0, 1.0],
            1e-7,
        );
        assert!(r.is_none());
    }

    #[test]
    fn diagonal_through_origin() {
        let dir = DVec2::new(1.0, 1.0).normalize();
        let r = classify_lin2d(&line(DVec2::new(0.0, 0.0), dir), [0.0, 1.0, 0.0, 1.0], 1e-7);
        assert!(r.is_some());
        let [p1, p2] = r.unwrap();
        let expected = (2.0_f64).sqrt();
        assert!(
            (p2 - p1 - expected).abs() < 1e-6,
            "span {}, expected {}",
            p2 - p1,
            expected
        );
    }

    #[test]
    fn non_unit_domain() {
        let r = classify_lin2d(
            &line(DVec2::new(2.0, 0.0), DVec2::new(1.0, 0.0)),
            [2.0, 5.0, 0.0, 5.0],
            1e-7,
        );
        assert!(r.is_some());
        assert!((r.unwrap()[1] - r.unwrap()[0] - 3.0).abs() < 1e-10);
    }
}
