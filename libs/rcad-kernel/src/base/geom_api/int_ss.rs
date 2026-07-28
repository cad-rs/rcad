//! Surface-surface intersection (GeomAPI_IntSS).
//!
//! OCCT TKGeomBase GeomAPI package: GeomAPI_IntSS.
//!
//! Computes intersection curves between two surfaces.
//!
//! Algorithm (simplified OCCT GeomInt_IntSS):
//! 1. Discretize surface S2 into a UV grid
//! 2. For each grid point on S2, compute distance to S1 via ExtPS
//! 3. Find seed points where the two surfaces are close
//! 4. Cluster seeds into intersection curves
//! 5. March along each curve and fit a BSpline

#![allow(clippy::manual_clamp)]

use glam::DVec3;

use crate::geom::{Curve3, CurveEval, Surface3, SurfaceEval};

const TOL: f64 = 1e-7;

/// Surface-surface intersection algorithm.
///
/// OCCT: `GeomAPI_IntSS`.
pub struct IntSS {
    done: bool,
    lines: Vec<Curve3>,
}

impl IntSS {
    /// Default constructor.
    ///
    /// OCCT: `GeomAPI_IntSS()`.
    pub fn new() -> Self {
        IntSS {
            done: false,
            lines: Vec::new(),
        }
    }

    /// Constructor with two surfaces and tolerance.
    ///
    /// OCCT: `GeomAPI_IntSS(S1, S2, Tol)`.
    pub fn with_surfaces(s1: &Surface3, s2: &Surface3, _tol: f64) -> Self {
        let mut intss = IntSS::new();
        intss.perform(s1, s2, _tol);
        intss
    }

    /// Perform the intersection.
    ///
    /// OCCT: `Perform(S1, S2, Tol)`.
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3, _tol: f64) {
        self.lines.clear();

        let dom1 = s1.default_domain();
        let dom2 = s2.default_domain();

        let (u1_min, u1_max, v1_min, v1_max) = (dom1[0], dom1[1], dom1[2], dom1[3]);
        let (u2_min, u2_max, v2_min, v2_max) = (dom2[0], dom2[1], dom2[2], dom2[3]);

        if !u1_min.is_finite() || !u2_min.is_finite() {
            self.done = true;
            return;
        }

        // Strategy: for analytic surface pairs, use ProjLib-style projection.
        // For general pairs, use sampling-based approach.
        if let Some(curves) = intersect_analytic_pair(s1, s2) {
            self.lines = curves;
            self.done = true;
            return;
        }

        // General case: sample S2 at grid points, find proximity to S1
        const N_U: usize = 20;
        const N_V: usize = 20;
        let mut seeds: Vec<(f64, f64, f64, f64)> = Vec::new(); // (u1, v1, u2, v2)

        for i in 0..=N_U {
            let u2 = u2_min + (u2_max - u2_min) * (i as f64) / (N_U as f64);
            for j in 0..=N_V {
                let v2 = v2_min + (v2_max - v2_min) * (j as f64) / (N_V as f64);
                let p2 = s2.point_at(u2, v2);

                // Find closest point on S1
                use crate::base::extrema::ExtPS;
                let ext = ExtPS::with_domain(p2, s1, u1_min, u1_max, v1_min, v1_max, TOL, TOL);
                if ext.nb_ext() > 0 {
                    let p_on_s1 = ext.point(1);
                    let dist = (p_on_s1.point - p2).length();
                    if dist < TOL * 100.0 {
                        seeds.push((p_on_s1.u, p_on_s1.v, u2, v2));
                    }
                }
            }
        }

        if seeds.is_empty() {
            self.done = true;
            return;
        }

        // Cluster seeds into connected curves
        let clusters = cluster_seeds(&seeds);
        for cluster in &clusters {
            if cluster.len() < 3 {
                continue;
            }
            let pts_3d: Vec<DVec3> = cluster.iter().map(|&(u1, v1, _, _)| s1.point_at(u1, v1)).collect();
            let curve = fit_bspline(&pts_3d);
            self.lines.push(curve);
        }

        self.done = true;
    }

    /// Returns true if the intersection was computed.
    ///
    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Returns the number of computed intersection curves.
    ///
    /// OCCT: `NbLines()`.
    pub fn nb_lines(&self) -> usize {
        self.lines.len()
    }

    /// Returns the intersection curve at 1-based index.
    ///
    /// OCCT: `Line(Index)`.
    pub fn line(&self, index: usize) -> &Curve3 {
        assert!(index >= 1 && index <= self.lines.len(), "IntSS: index out of range");
        &self.lines[index - 1]
    }
}

impl Default for IntSS {
    fn default() -> Self {
        Self::new()
    }
}

/// Try analytic surface intersection (plane-plane, plane-cylinder, etc.)
fn intersect_analytic_pair(s1: &Surface3, s2: &Surface3) -> Option<Vec<Curve3>> {
    match (s1, s2) {
        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
            // Plane-plane intersection: line
            use crate::base::int_ana::intersect_plane_plane_intana;
            match intersect_plane_plane_intana(p1, p2) {
                crate::base::int_ana::PlnPlnResult::Line(l) => {
                    Some(vec![Curve3::Line(l)])
                }
                _ => Some(vec![]),
            }
        }
        (Surface3::Plane(pl), Surface3::Cylinder(cy)) |
        (Surface3::Cylinder(cy), Surface3::Plane(pl)) => {
            use crate::base::int_ana::intersect_plane_cylinder_intana;
            match intersect_plane_cylinder_intana(pl, cy) {
                crate::base::int_ana::PlnCylResult::Circle(c) => Some(vec![Curve3::Circle(c)]),
                crate::base::int_ana::PlnCylResult::Ellipse(e) => Some(vec![Curve3::Ellipse(e)]),
                crate::base::int_ana::PlnCylResult::TangentLine(l) |
                crate::base::int_ana::PlnCylResult::TwoLines(l, _) => {
                    Some(vec![Curve3::Line(l)])
                }
                _ => Some(vec![]),
            }
        }
        (Surface3::Plane(pl), Surface3::Sphere(sp)) |
        (Surface3::Sphere(sp), Surface3::Plane(pl)) => {
            use crate::base::int_ana::intersect_plane_sphere_intana;
            match intersect_plane_sphere_intana(pl, sp) {
                crate::base::int_ana::PlnSphResult::Circle(c) => Some(vec![Curve3::Circle(c)]),
                _ => Some(vec![]),
            }
        }
        (Surface3::Plane(pl), Surface3::Cone(co)) |
        (Surface3::Cone(co), Surface3::Plane(pl)) => {
            use crate::base::int_ana::intersect_plane_cone_intana;
            match intersect_plane_cone_intana(pl, co) {
                crate::base::int_ana::PlnConResult::Circle(c) => Some(vec![Curve3::Circle(c)]),
                crate::base::int_ana::PlnConResult::Ellipse(e) => Some(vec![Curve3::Ellipse(e)]),
                _ => Some(vec![]),
            }
        }
        _ => None, // Not an analytic pair we handle
    }
}

/// Cluster seed points into connected curves.
fn cluster_seeds(seeds: &[(f64, f64, f64, f64)]) -> Vec<Vec<(f64, f64, f64, f64)>> {
    if seeds.is_empty() {
        return vec![];
    }

    let tol_u = 0.1; // Normalized tolerance for U clustering
    let tol_v = 0.1;

    let mut clusters: Vec<Vec<(f64, f64, f64, f64)>> = Vec::new();
    let mut assigned = vec![false; seeds.len()];

    for i in 0..seeds.len() {
        if assigned[i] {
            continue;
        }

        let mut cluster = vec![seeds[i]];
        assigned[i] = true;

        // Grow cluster: add nearby seeds
        let mut changed = true;
        while changed {
            changed = false;
            for j in 0..seeds.len() {
                if assigned[j] {
                    continue;
                }
                let s = seeds[j];
                let near = cluster.iter().any(|&c| {
                    let du = (s.0 - c.0).abs();
                    let dv = (s.1 - c.1).abs();
                    du < tol_u && dv < tol_v
                });
                if near {
                    cluster.push(s);
                    assigned[j] = true;
                    changed = true;
                }
            }
        }

        // Sort cluster by parameter order
        cluster.sort_by(|a, b| {
            let a_avg = (a.0 + a.1) * 0.5;
            let b_avg = (b.0 + b.1) * 0.5;
            a_avg.partial_cmp(&b_avg).unwrap()
        });

        clusters.push(cluster);
    }

    clusters
}

/// Fit a cubic BSpline through 3D points.
fn fit_bspline(pts: &[DVec3]) -> Curve3 {
    let n = pts.len();
    if n < 2 {
        return Curve3::Line(crate::geom::Line3::new(DVec3::ZERO, DVec3::X));
    }
    if n == 2 {
        return Curve3::Line(crate::geom::Line3::new(pts[0], pts[1] - pts[0]));
    }

    // Chord-length parametrization
    let mut params = vec![0.0_f64; n];
    for i in 1..n {
        let d = (pts[i] - pts[i - 1]).length();
        params[i] = params[i - 1] + d.max(1e-15);
    }
    let total = params[n - 1];
    for p in &mut params {
        *p /= total;
    }

    let degree = 3.min(n - 1);
    let n_knots = n + degree + 1;
    let mut knots = vec![0.0_f64; n_knots];
    for k in &mut knots[..=degree] {
        *k = params[0];
    }
    for j in 1..n - degree {
        let mut sum = 0.0;
        for i in j..j + degree {
            sum += params[i];
        }
        knots[j + degree] = sum / (degree as f64);
    }
    for k in &mut knots[n_knots - degree - 1..] {
        *k = params[n - 1];
    }

    Curve3::BSpline(crate::geom::BSplineCurve3 {
        degree,
        knots,
        control_points: pts.to_vec(),
        weights: vec![],
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn test_plane_plane_intersect() {
        let p1 = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let p2 = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::X));
        let intss = IntSS::with_surfaces(&p1, &p2, 1e-7);
        assert!(intss.is_done());
        assert!(intss.nb_lines() > 0);
    }

    #[test]
    fn test_plane_cylinder_intersect() {
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let intss = IntSS::with_surfaces(&plane, &cyl, 1e-7);
        assert!(intss.is_done());
        assert!(intss.nb_lines() > 0);
    }
}
