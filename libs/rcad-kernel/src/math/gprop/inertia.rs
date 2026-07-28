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
