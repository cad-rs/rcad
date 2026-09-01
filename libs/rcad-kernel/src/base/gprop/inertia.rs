//! OCCT GProp_PrincipalProps: inertia tensor computation.
//!
//! Split from the original properties.rs.

use glam::DVec3;
use serde::{Deserialize, Serialize};

use crate::topo::topods;

/// Inertia tensor of a BRep solid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InertiaTensor {
    pub ixx: f64, pub iyy: f64, pub izz: f64,
    pub ixy: f64, pub ixz: f64, pub iyz: f64,
}

/// Compute the inertia tensor of a BRep solid.
pub fn inertia_tensor(brep: &topods::BRep) -> InertiaTensor {
    let vol = super::volume::signed_volume(brep);
    let c = super::volume::centroid(brep);
    let faces = super::tri::face_flat_iter(brep);
    let mut ixx = 0.0; let mut iyy = 0.0; let mut izz = 0.0;
    let mut ixy = 0.0; let mut ixz = 0.0; let mut iyz = 0.0;

    for (fi, face) in &faces {
        let tris = super::tri::face_triangles_pub(brep, *fi);
        for [a, b, cc] in &tris {
            let v0 = *a - c; let v1 = *b - c; let v2 = *cc - c;
            let tv = super::tri::tet_signed_volume(*a, *b, *cc);
            // Approximate per-tetrahedron inertia contribution
            // For a tetrahedron with vertices relative to centroid:
            // I contribution ≈ (m/20) * Σ(2*vi*vi + vj*vj + vk*vk) / 6
            let (x0, y0, z0) = (v0.x, v0.y, v0.z);
            let (x1, y1, z1) = (v1.x, v1.y, v1.z);
            let (x2, y2, z2) = (v2.x, v2.y, v2.z);
            let m = tv.abs() * 6.0; // tetrahedron volume contribution
            // Second moment: ∫∫∫ (r²δ - r⊗r) dV over tetrahedron
            // Simplified: use vertex average × volume
            let avg_xx = (x0*x0 + x1*x1 + x2*x2) / 3.0;
            let avg_yy = (y0*y0 + y1*y1 + y2*y2) / 3.0;
            let avg_zz = (z0*z0 + z1*z1 + z2*z2) / 3.0;
            let avg_xy = (x0*y0 + x1*y1 + x2*y2) / 3.0;
            let avg_xz = (x0*z0 + x1*z1 + x2*z2) / 3.0;
            let avg_yz = (y0*z0 + y1*z1 + y2*z2) / 3.0;
            let scale = m / 20.0;
            ixx += scale * (2.0*avg_xx + avg_yy + avg_zz);
            iyy += scale * (2.0*avg_yy + avg_xx + avg_zz);
            izz += scale * (2.0*avg_zz + avg_xx + avg_yy);
            ixy += scale * avg_xy;
            ixz += scale * avg_xz;
            iyz += scale * avg_yz;
        }
    }

    InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz }
}

/// OCCT GProp_PrincipalProps — principal moments of inertia and symmetry.
#[derive(Debug, Clone)]
pub struct PrincipalProps {
    /// OCCT GProp_PrincipalProps::HasSymmetryAxis — two principal moments
    /// are (nearly) equal, i.e. the solid has an axis of symmetry.
    pub has_symmetry_axis: bool,
    /// Principal moments of inertia, sorted ascending.
    pub moments: [f64; 3],
}

/// Determinant of a 3x3 matrix (row-major).
fn det3(
    a: f64, b: f64, c: f64,
    d: f64, e: f64, f: f64,
    g: f64, h: f64, i: f64,
) -> f64 {
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

/// Eigenvalues of a 3x3 symmetric matrix (closed-form trigonometric method).
/// Matrix (row-major): [[a, d, e], [d, b, f], [e, f, c]].
fn symmetric_eigenvalues(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> [f64; 3] {
    let p1 = d * d + e * e + f * f;
    if p1 < 1e-30 {
        // Diagonal matrix.
        let mut v = [a, b, c];
        v.sort_by(|x, y| x.partial_cmp(y).unwrap());
        return v;
    }
    let q = (a + b + c) / 3.0;
    let p2 = (a - q).powi(2) + (b - q).powi(2) + (c - q).powi(2) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    // B = (A - q*I) / p — trace-free part.
    let b11 = (a - q) / p;
    let b22 = (b - q) / p;
    let b33 = (c - q) / p;
    let b12 = d / p;
    let b13 = e / p;
    let b23 = f / p;
    let r = det3(b11, b12, b13, b12, b22, b23, b13, b23, b33) / 2.0;
    let r = r.clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let pi_3 = std::f64::consts::PI / 3.0;
    let mut eig = [
        q + 2.0 * p * phi.cos(),
        q + 2.0 * p * (phi + 2.0 * pi_3).cos(),
        q + 2.0 * p * (phi + 4.0 * pi_3).cos(),
    ];
    eig.sort_by(|x, y| x.partial_cmp(y).unwrap());
    eig
}

/// OCCT GProp_GProps::PrincipalProperties — principal moments from the
/// inertia tensor.  `HasSymmetryAxis` is true when two principal moments
/// are equal within a relative tolerance.
pub fn principal_properties(brep: &topods::BRep) -> PrincipalProps {
    let t = inertia_tensor(brep);
    let moments = symmetric_eigenvalues(t.ixx, t.iyy, t.izz, t.ixy, t.ixz, t.iyz);
    let scale = moments[2].abs().max(1e-12);
    let has_symmetry_axis = (moments[0] - moments[1]).abs() <= 1e-6 * scale
        || (moments[1] - moments[2]).abs() <= 1e-6 * scale;
    PrincipalProps {
        has_symmetry_axis,
        moments,
    }
}
