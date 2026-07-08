//! OCCT GeomHash — tolerance-based geometric hashing for Curve3/Surface3.
//!
//! OCCT source: src/ModelingData/TKG3d/GeomHash/
//!
//! Provides `CurveHasher` with per-type hash + equivalence comparison,
//! matching the GeomHash_CurveHasher / GeomHash_SurfaceHasher API.

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use glam::DVec3;
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;

/// Quantize a f64 value into a u64 using the given tolerance.
fn quantize(val: f64, tol: f64) -> u64 {
    if tol <= 0.0 { return val.to_bits(); }
    let q = (val / tol).round();
    q as i64 as u64
}

/// Quantize a DVec3 into (u64, u64, u64) using the given tolerance.
fn quantize_vec3(v: DVec3, tol: f64) -> (u64, u64, u64) {
    (quantize(v.x, tol), quantize(v.y, tol), quantize(v.z, tol))
}

fn combine_hashes(h1: u64, h2: u64) -> u64 {
    h1.wrapping_mul(6364136223846793005).wrapping_add(h2).wrapping_add(1)
}

/// Compute hash for a Curve3 using tolerance-based quantization.
pub fn hash_curve(curve: &Curve3, comp_tol: f64) -> u64 {
    let mut h = 0u64;
    match curve {
        Curve3::Line(l) => {
            let (ox, oy, oz) = quantize_vec3(l.origin, comp_tol);
            let (dx, dy, dz) = quantize_vec3(l.direction, comp_tol);
            h = combine_hashes(h, ox);
            h = combine_hashes(h, oy);
            h = combine_hashes(h, oz);
            h = combine_hashes(h, dx);
            h = combine_hashes(h, dy);
            h = combine_hashes(h, dz);
        }
        Curve3::Circle(c) => {
            let (cx, cy, cz) = quantize_vec3(c.center, comp_tol);
            let (nx, ny, nz) = quantize_vec3(c.normal, comp_tol);
            let r = quantize(c.radius, comp_tol);
            h = combine_hashes(h, cx);
            h = combine_hashes(h, cy);
            h = combine_hashes(h, cz);
            h = combine_hashes(h, nx);
            h = combine_hashes(h, ny);
            h = combine_hashes(h, nz);
            h = combine_hashes(h, r);
        }
        Curve3::Ellipse(e) => {
            let (cx, cy, cz) = quantize_vec3(e.center, comp_tol);
            let (mx, my, mz) = quantize_vec3(e.major_dir, comp_tol);
            let (nx, ny, nz) = quantize_vec3(e.normal, comp_tol);
            let ma = quantize(e.major_radius, comp_tol);
            let mi = quantize(e.minor_radius, comp_tol);
            h = combine_hashes(h, cx); h = combine_hashes(h, cy); h = combine_hashes(h, cz);
            h = combine_hashes(h, mx); h = combine_hashes(h, my); h = combine_hashes(h, mz);
            h = combine_hashes(h, nx); h = combine_hashes(h, ny); h = combine_hashes(h, nz);
            h = combine_hashes(h, ma); h = combine_hashes(h, mi);
        }
        Curve3::Hyperbola(hyp) => {
            h = combine_hashes(h, quantize_vec3(hyp.center, comp_tol).0);
            h = combine_hashes(h, quantize(hyp.semi_major, comp_tol));
            h = combine_hashes(h, quantize(hyp.semi_minor, comp_tol));
        }
        Curve3::Parabola(p) => {
            h = combine_hashes(h, quantize_vec3(p.vertex, comp_tol).0);
            h = combine_hashes(h, quantize(p.focal_param, comp_tol));
        }
        Curve3::Bezier(b) => {
            for p in &b.control_points {
                let (qx, qy, qz) = quantize_vec3(*p, comp_tol);
                h = combine_hashes(h, qx);
                h = combine_hashes(h, qy);
                h = combine_hashes(h, qz);
            }
        }
        Curve3::BSpline(b) => {
            for p in &b.control_points {
                let (qx, qy, qz) = quantize_vec3(*p, comp_tol);
                h = combine_hashes(h, qx); h = combine_hashes(h, qy); h = combine_hashes(h, qz);
            }
            for w in &b.weights {
                h = combine_hashes(h, quantize(*w, comp_tol));
            }
        }
        Curve3::CircularHelix(hc) => {
            h = combine_hashes(h, quantize_vec3(hc.origin, comp_tol).0);
            h = combine_hashes(h, quantize(hc.radius, comp_tol));
            h = combine_hashes(h, quantize(hc.pitch, comp_tol));
        }
        Curve3::SineWave(sw) => {
            h = combine_hashes(h, quantize(sw.amplitude, comp_tol));
            h = combine_hashes(h, quantize(sw.frequency, comp_tol));
            h = combine_hashes(h, quantize(sw.phase, comp_tol));
        }
        Curve3::Offset(o) => {
            h = combine_hashes(h, hash_curve(&o.basis, comp_tol));
            h = combine_hashes(h, quantize(o.offset_distance, comp_tol));
        }
    }
    h
}

/// Check if two Curve3 are equivalent within tolerance.
pub fn curves_equivalent(a: &Curve3, b: &Curve3, tol: f64) -> bool {
    use Curve3::*;
    match (a, b) {
        (Line(la), Line(lb)) =>
            (la.origin - lb.origin).length() < tol && (la.direction - lb.direction).length() < tol,
        (Circle(ca), Circle(cb)) =>
            (ca.center - cb.center).length() < tol && (ca.radius - cb.radius).abs() < tol
            && (ca.normal - cb.normal).length() < tol,
        (Ellipse(ea), Ellipse(eb)) =>
            (ea.center - eb.center).length() < tol && (ea.major_radius - eb.major_radius).abs() < tol
            && (ea.minor_radius - eb.minor_radius).abs() < tol
            && (ea.major_dir - eb.major_dir).length() < tol,
        (Hyperbola(ha), Hyperbola(hb)) =>
            (ha.center - hb.center).length() < tol && (ha.semi_major - hb.semi_major).abs() < tol
            && (ha.semi_minor - hb.semi_minor).abs() < tol,
        (Parabola(pa), Parabola(pb)) =>
            (pa.vertex - pb.vertex).length() < tol && (pa.focal_param - pb.focal_param).abs() < tol,
        (Bezier(ba), Bezier(bb)) => {
            ba.control_points.len() == bb.control_points.len()
            && ba.control_points.iter().zip(&bb.control_points)
                .all(|(pa, pb)| (*pa - *pb).length() < tol)
        }
        (BSpline(ba), BSpline(bb)) => {
            ba.degree == bb.degree && ba.control_points.len() == bb.control_points.len()
            && ba.control_points.iter().zip(&bb.control_points)
                .all(|(pa, pb)| (*pa - *pb).length() < tol)
            && ba.weights.iter().zip(&bb.weights)
                .all(|(wa, wb)| (*wa - *wb).abs() < tol)
        }
        (Offset(oa), Offset(ob)) =>
            (oa.offset_distance - ob.offset_distance).abs() < tol
            && curves_equivalent(&oa.basis, &ob.basis, tol),
        _ => false,
    }
}

// =============================================================================
// Tests — GeomHash_CurveHasher_Test.cxx
// =============================================================================

#[cfg(test)]
mod curve_hasher_tests {
    use super::*;

    #[test]
    fn line_copied_same_hash() {
        let l1 = Curve3::Line(Line3 { origin: DVec3::new(1.0, 2.0, 3.0), direction: DVec3::X });
        let l2 = l1.clone();
        assert_eq!(hash_curve(&l1, TOL), hash_curve(&l2, TOL));
        assert!(curves_equivalent(&l1, &l2, TOL));
    }

    #[test]
    fn line_different_different_hash() {
        let l1 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let l2 = Curve3::Line(Line3 { origin: DVec3::new(1.0, 0.0, 0.0), direction: DVec3::X });
        assert_ne!(hash_curve(&l1, TOL), hash_curve(&l2, TOL));
        assert!(!curves_equivalent(&l1, &l2, TOL));
    }

    #[test]
    fn circle_copied_same_hash() {
        let c1 = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let c2 = c1.clone();
        assert_eq!(hash_curve(&c1, TOL), hash_curve(&c2, TOL));
        assert!(curves_equivalent(&c1, &c2, TOL));
    }

    #[test]
    fn circle_different_radius_different_hash() {
        let c1 = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let c2 = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 10.0));
        assert_ne!(hash_curve(&c1, TOL), hash_curve(&c2, TOL));
        assert!(!curves_equivalent(&c1, &c2, TOL));
    }

    #[test]
    fn ellipse_copied_same_hash() {
        let e1 = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 10.0, minor_radius: 5.0,
        });
        let e2 = e1.clone();
        assert_eq!(hash_curve(&e1, TOL), hash_curve(&e2, TOL));
        assert!(curves_equivalent(&e1, &e2, TOL));
    }

    #[test]
    fn hyperbola_copied_same_hash() {
        let h1 = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            semi_major: 5.0, semi_minor: 3.0,
        });
        let h2 = h1.clone();
        assert_eq!(hash_curve(&h1, TOL), hash_curve(&h2, TOL));
        assert!(curves_equivalent(&h1, &h2, TOL));
    }

    #[test]
    fn parabola_copied_same_hash() {
        let p1 = Curve3::Parabola(Parabola3 {
            center: DVec3::ZERO, normal: DVec3::Z, x_dir: DVec3::X,
            focal_param: 2.0,
        });
        let p2 = p1.clone();
        assert_eq!(hash_curve(&p1, TOL), hash_curve(&p2, TOL));
        assert!(curves_equivalent(&p1, &p2, TOL));
    }

    #[test]
    fn bezier_copied_same_hash() {
        let b1 = Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(3.0, 2.0, 0.0), DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        });
        let b2 = b1.clone();
        assert_eq!(hash_curve(&b1, TOL), hash_curve(&b2, TOL));
        assert!(curves_equivalent(&b1, &b2, TOL));
    }

    #[test]
    fn bspline_copied_same_hash() {
        let s1 = Curve3::BSpline(BSplineCurve3 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0),
                DVec3::new(2.0, 1.0, 0.0), DVec3::new(3.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        });
        let s2 = s1.clone();
        assert_eq!(hash_curve(&s1, TOL), hash_curve(&s2, TOL));
        assert!(curves_equivalent(&s1, &s2, TOL));
    }

    #[test]
    fn offset_curve_copied_same_hash() {
        let basis = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let o1 = Curve3::Offset(Box::new(OffsetCurve3 { basis: Box::new(basis.clone()), offset_distance: 2.0 }));
        let o2 = o1.clone();
        assert_eq!(hash_curve(&o1, TOL), hash_curve(&o2, TOL));
        assert!(curves_equivalent(&o1, &o2, TOL));
    }
}
