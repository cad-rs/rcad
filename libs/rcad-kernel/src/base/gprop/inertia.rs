//! OCCT GProp_PrincipalProps: inertia tensor computation.
//!
//! Split from the original properties.rs.

use glam::DVec3;
use serde::{Deserialize, Serialize};

use crate::topo::topods;

use super::volume::shape_vinert;

/// Inertia tensor of a BRep solid — the OCCT GProp_GProps::MatrixOfInertia
/// (GProp_GProps.cxx L110-115): the accumulated per-face Vinert inertia about
/// the origin minus the Huygens shift
/// `HOperator(g, gp::Origin(), dim)` (GProp.cxx HOperator) from the center of
/// mass to the origin, i.e. the inertia tensor about the center of mass.
/// Matrix elements follow the gp_Mat convention of BRepGProp_Gauss::convert
/// (BRepGProp_Gauss.cxx L487-489): [[Ixx, -Ixy, -Ixz], [-Ixy, Iyy, -Iyz],
/// [-Ixz, -Iyz, Izz]] with the accumulated Ixy = -∫xy dV etc.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InertiaTensor {
    pub ixx: f64, pub iyy: f64, pub izz: f64,
    pub ixy: f64, pub ixz: f64, pub iyz: f64,
}

/// OCCT GProp::HOperator (GProp.cxx L27-40): the inertia matrix of a point
/// mass `dim` located at G about the point O.
fn h_operator(g: DVec3, dim: f64) -> [f64; 6] {
    let x = g.x;
    let y = g.y;
    let z = g.z;
    [
        dim * (y * y + z * z),
        dim * (x * x + z * z),
        dim * (x * x + y * y),
        -dim * x * y,
        -dim * x * z,
        -dim * y * z,
    ]
}

/// Compute the inertia tensor of a BRep solid.
///
/// OCCT BRepGProp::VolumeProperties (BRepGProp.cxx L555-586 → L413-441 →
/// volumePropertiesFaces L298-409): every face is integrated with
/// BRepGProp_Vinert about the shape origin (Eps = 1.0 fixed-order Gauss,
/// aCoeff = {0,0,0}, isByPoint = true) and accumulated per occurrence
/// (GProp_GProps::Add L42-64, all faces share the location so the matrices
/// add directly).  The resulting tensor is about the origin; the reported
/// matrix is shifted to the center of mass via the Huygens theorem
/// (GProp_GProps::MatrixOfInertia L110-115).
pub fn inertia_tensor(brep: &topods::BRep) -> InertiaTensor {
    let v = shape_vinert(brep);
    let dim = v.mass;
    // GProp_GProps::Add L50-58: the accumulated gravity center (relative to
    // the origin) = sum of the per-face first moments / dim.
    let g = if dim.abs() >= 1e-30 {
        DVec3::new(v.ix / dim, v.iy / dim, v.iz / dim)
    } else {
        DVec3::ZERO
    };
    // GProp_GProps::MatrixOfInertia (L110-115): inertia - HOperator(g, origin,
    // dim).  The accumulated matrix is the BRepGProp_Gauss::convert layout
    // (L487-489): [[Ixx, -Ixy, -Ixz], [-Ixy, Iyy, -Iyz], [-Ixz, -Iyz, Izz]]
    // with Ixy = -∫xy dV etc. — i.e. the matrix elements are (Ixx, Iyy, Izz,
    // -Ixy, -Ixz, -Iyz).
    let h = h_operator(g, dim);
    InertiaTensor {
        ixx: v.ixx - h[0],
        iyy: v.iyy - h[1],
        izz: v.izz - h[2],
        ixy: -v.ixy - h[3],
        ixz: -v.ixz - h[4],
        iyz: -v.iyz - h[5],
    }
}

/// OCCT GProp_PrincipalProps — principal moments of inertia and symmetry.
#[derive(Debug, Clone)]
pub struct PrincipalProps {
    /// OCCT GProp_PrincipalProps::HasSymmetryAxis (GProp_PrincipalProps.cxx
    /// L54-60): two principal moments are equal within a RELATIVE tolerance
    /// of 1.e-10, i.e. the solid has an axis of symmetry.
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

/// OCCT GProp_GProps::PrincipalProperties (GProp_GProps.cxx L154-197) +
/// GProp_PrincipalProps::HasSymmetryAxis (GProp_PrincipalProps.cxx L54-60):
/// the principal moments are the eigenvalues of the matrix of inertia about
/// the center of mass (math_Jacobi), and `HasSymmetryAxis` is true when two
/// principal moments are equal within the relative tolerance 1.e-10.
pub fn principal_properties(brep: &topods::BRep) -> PrincipalProps {
    let t = inertia_tensor(brep);
    let moments = symmetric_eigenvalues(t.ixx, t.iyy, t.izz, t.ixy, t.ixz, t.iyz);
    // GProp_PrincipalProps::HasSymmetryAxis (L54-60):
    //   Eps1 = |i1| * 1.e-10, Eps2 = |i2| * 1.e-10
    //   |i1 - i2| <= Eps1 || |i1 - i3| <= Eps1 || |i2 - i3| <= Eps2
    let a_rel_tol = 1.0e-10;
    let eps1 = moments[0].abs() * a_rel_tol;
    let eps2 = moments[1].abs() * a_rel_tol;
    let has_symmetry_axis = (moments[0] - moments[1]).abs() <= eps1
        || (moments[0] - moments[2]).abs() <= eps1
        || (moments[1] - moments[2]).abs() <= eps2;
    PrincipalProps {
        has_symmetry_axis,
        moments,
    }
}
