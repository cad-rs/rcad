//! OCCT BRepGProp::VolumeProperties: volume computation.
//!
//! 1:1 translation of the `checkprops -v` path.  BRepGProp::VolumeProperties
//! (BRepGProp.cxx L555-588 -> volumePropertiesFaces L298-409) integrates every
//! face with BRepGProp_Vinert (BRepGProp_Vinert.cxx L197-204): aCoeff = {0,0,0}
//! and theIsByPoint = true, i.e. BRepGProp_Gauss::Compute (Vinert,
//! BRepGProp_Gauss.cxx L1215-1302) with the divergence-theorem elementary part
//! computeVInertiaOfElementaryPart isByPoint (L306-339):
//!
//!   dv = r . (N * w);   Mass += dv / 3
//!
//! over the face boundary edge pcurves (Green-theorem line integral), exactly
//! mirroring the surface-area path (base::gprop::surface, Sinert).

use glam::DVec3;

use crate::BRep;
use crate::geom::{Curve2dEval, Surface3, SurfaceEval};
use crate::topo::topods;
use crate::topo::topo_shape::Shape;
use crate::topo::topods::{self as td, BRepTool};
use crate::topo::topology::Face;
use crate::base::gprop::surface::{
    curve_integration_order, face_uv_bounds, occt_gauss, surface_u_v_integration_order, EPS_DIM,
};
use crate::base::gprop::tri::{face_flat_iter, face_triangles_pub, tet_signed_volume};

/// OCCT BRepGProp_Gauss::Compute (BRepGProp_Gauss.cxx L1215-1302, Vinert) for
/// one face with wires: line integral over the boundary edge pcurves.
///
///   Mass = sum over edges of
///     lr * sum_i { ur * sum_j ( (1/3) r(u_j, v_i) . N(u_j, v_i) * (dv/dl)_i * w_i * w_j ) }
///
/// with r = point - origin (the shape origin, L123-127 of BRepGProp.cxx),
/// N = D1U x D1V (BRepGProp_Face::Normal L201-210, mySReverse=false: rcad
/// faces are explored FORWARD, see face_flat_iter), and the v coordinate
/// clamped to the face UV bounds (OCC104, L1266-1267).
pub fn face_volume_gauss_domain(brep: &BRep, face: &Face, fi: usize) -> f64 {
    let surf_idx = match brep.tshapes.get(fi).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts {
            fd.surface.clone()
        } else {
            None
        }
    }) {
        Some(s) => s,
        None => return 0.0,
    };
    let surf = &surf_idx;
    // L1226-1228: theSurface.Bounds — face UV bounds; u1 feeds the inner
    // integration offset, u2/v are clamped per-sample (OCC104 L1266-1267).
    let [u1, u2_face, v1, v2] = face_uv_bounds(brep, face, fi, surf);
    // L1240-1243: aNbGaussgp_Pnts = min(max(edge IntegrationOrder,
    // VIntegrationOrder), GaussPointsMax).
    let (_nu, nv) = surface_u_v_integration_order(surf);

    let face_shape = Shape::from_parts(
        brep.tshapes[fi].clone(),
        fi,
        0,
        td::Orientation::Forward,
    );

    let mut an_inertia = 0.0f64;
    // Domain: every edge of the face (outer + inner wires).
    let mut edges = face.outer_wire.edges.clone();
    for w in &face.inner_wires {
        edges.extend(w.edges.iter().copied());
    }
    for we in &edges {
        // OCCT BRepGProp_Domain::Next (BRepGProp_Domain.cxx L27-38) skips
        // INTERNAL and EXTERNAL edges.
        if we.internal {
            continue;
        }
        let edge_shape = Shape::from_parts(
            brep.tshapes[we.idx].clone(),
            we.idx,
            we.location,
            if we.forward {
                td::Orientation::Forward
            } else {
                td::Orientation::Reversed
            },
        );
        // BRepGProp_Face::Load (BRepGProp_Face.cxx L164-185): the edge's
        // pcurve on the face; a REVERSED edge on a closed surface (seam) uses
        // PCurve2; a null pcurve makes Load return false and Compute return
        // immediately with the face mass left at 0.
        let pc = brep.curve_on_surface(&edge_shape, &face_shape);
        let (c, a0, b0) = if !we.forward {
            match brep.curve_on_surface_second(&edge_shape, &face_shape) {
                Some(v) => v,
                None => match pc {
                    Some(v) => v.clone(),
                    None => return 0.0,
                },
            }
        } else {
            match pc {
                Some(v) => v.clone(),
                None => return 0.0,
            }
        };
        // Load (L176-182): REVERSED edge reverses the pcurve and its range.
        let reversed = !we.forward;
        let (a, b) = if reversed {
            let x = a0;
            (c.reversed_parameter(b0), c.reversed_parameter(x))
        } else {
            (a0, b0)
        };
        // L1240-1243: aNbGaussgp_Pnts = min(max(edge IntOrder, VIntOrder), GPM)
        let nb_gauss = curve_integration_order(&c).max(nv).min(61);
        let (cp, cw) = match occt_gauss(nb_gauss) {
            Some(v) => v,
            None => return 0.0,
        };
        // L1250-1253: l1/l2 = the loaded pcurve First/LastParameter.
        let (l1, l2) = (a, b);
        let lm = 0.5 * (l2 + l1);
        let lr = 0.5 * (l2 - l1);

        let mut a_c_inertia = 0.0;
        for i in 0..nb_gauss {
            let l = lm + lr * cp[i];
            // D12d (L1260-1263): point and first derivative of the 2D curve.
            let (puv, vuv) = if reversed {
                let rp = c.reversed_parameter(l);
                (c.point_at(rp), -c.derivative_at(rp))
            } else {
                (c.point_at(l), c.derivative_at(l))
            };
            // L1265-1267 (OCC104): u2 = clamp(u2, u1, _u2); v = clamp(v, v1, v2)
            let u2 = puv.x.max(u1).min(u2_face);
            let v = puv.y.max(v1).min(v2);
            // L1269: Dul = dv/dl * w_i
            let dul = vuv.y * cw[i];
            // L1270-1271
            let um = 0.5 * (u2 + u1);
            let ur = 0.5 * (u2 - u1);
            // L1273-1291: aLocal = sum_j ( (1/3) r(u_j, v) . N(u_j, v) * Dul * wV_j )
            let mut a_local = 0.0;
            for j in 0..nb_gauss {
                let u = um + ur * cp[j];
                let a_weight = dul * cw[j];
                // BRepGProp_Face::Normal (L201-210): D1U x D1V; mySReverse is
                // false for rcad faces (explored FORWARD).
                let (p, du, dv) = surf.derivatives(u, v);
                let n = du.cross(dv);
                // computeVInertiaOfElementaryPart (L306-339, isByPoint): the
                // reference point is the shape origin — rcad faces are already
                // in world coordinates, so r = P.
                let dv = (p.x * n.x + p.y * n.y + p.z * n.z) * a_weight;
                a_local += dv / 3.0;
            }
            // L1293: aLocal *= ur
            a_c_inertia += a_local * ur;
        }
        // L1297-1298: aC *= lr; anInertia += aC
        an_inertia += a_c_inertia * lr;
    }
    // convert (L467-490): |Mass| >= EPS_DIM(1e-30) → mass else 0
    if an_inertia.abs() >= EPS_DIM {
        an_inertia
    } else {
        0.0
    }
}

/// OCCT BRepGProp_Gauss::Compute (BRepGProp_Gauss.cxx L1306-1393, Vinert) for
/// a natural-restriction face (no wires): Gauss-Legendre over the surface
/// natural parameter bounds.
///
///   Mass = vr * ur * sum_j ( wV_j * sum_i ( (1/3) r(u_i, v_j) . N(u_i, v_j) * wU_i ) )
pub fn face_volume_gauss_natural(brep: &BRep, fi: usize) -> f64 {
    let surf_idx = match brep.tshapes.get(fi).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts {
            fd.surface.clone()
        } else {
            None
        }
    }) {
        Some(s) => s,
        None => return 0.0,
    };
    let surf = &surf_idx;
    // L1314-1316: theSurface.Bounds — surface natural parameter domain.
    let [lower_u, upper_u, lower_v, upper_v] = surf.default_domain();
    // checkBounds (L418-429): an infinite bound makes the mass 0 (the convert
    // guard, |Mass| >= EPS_DIM is false for NaN).
    if !lower_u.is_finite() || !upper_u.is_finite() || !lower_v.is_finite() || !upper_v.is_finite()
    {
        return 0.0;
    }
    // L1318-1319: UIntegrationOrder / VIntegrationOrder, min GaussPointsMax(61)
    let (uo0, vo0) = surface_u_v_integration_order(surf);
    let uo = uo0.min(61);
    let vo = vo0.min(61);
    let (gpu, gwu) = match occt_gauss(uo) {
        Some(v) => v,
        None => return 0.0,
    };
    let (gpv, gwv) = match occt_gauss(vo) {
        Some(v) => v,
        None => return 0.0,
    };
    // L1332-1335
    let um = 0.5 * (upper_u + lower_u);
    let vm = 0.5 * (upper_v + lower_v);
    let ur = 0.5 * (upper_u - lower_u);
    let vr = 0.5 * (upper_v - lower_v);
    // L1341-1374: anInertia = sum_j ( wV_j * sum_i ( (1/3) r . N * wU_i ) )
    let mut an_inertia = 0.0;
    for j in 0..vo {
        let v = vm + vr * gpv[j];
        let mut an_inertia_of_elementary_part = 0.0;
        for i in 0..uo {
            let u = um + ur * gpu[i];
            let a_weight = gwu[i];
            let (p, du, dv) = surf.derivatives(u, v);
            let n = du.cross(dv);
            let dv_ = (p.x * n.x + p.y * n.y + p.z * n.z) * a_weight;
            an_inertia_of_elementary_part += dv_ / 3.0;
        }
        an_inertia += an_inertia_of_elementary_part * gwv[j];
    }
    // L1375: vr = vr * ur; L1392: mass *= vr (after the convert guard).
    let vr = vr * ur;
    let mass = if an_inertia.abs() >= EPS_DIM {
        an_inertia
    } else {
        0.0
    };
    mass * vr
}

/// Signed volume of a BRep solid — OCCT BRepGProp::VolumeProperties
/// (BRepGProp.cxx L298-409): every face integrated with BRepGProp_Vinert
/// (Eps = 1.0, fixed-order Gauss), IsNatRestr faces via the whole-surface
/// integral, the others via the boundary line integral.
pub fn signed_volume(brep: &topods::BRep) -> f64 {
    let mut vol = 0.0;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        let is_nat_restr = face.outer_wire.edges.is_empty() && face.inner_wires.is_empty();
        let v = if is_nat_restr {
            face_volume_gauss_natural(brep, *fi)
        } else {
            face_volume_gauss_domain(brep, face, *fi)
        };
        vol += v;
    }
    vol
}

/// Absolute volume of a BRep solid.
pub fn volume(brep: &topods::BRep) -> f64 {
    signed_volume(brep).abs()
}

/// Centroid of a BRep solid.
pub fn centroid(brep: &topods::BRep) -> DVec3 {
    let mut total_vol = 0.0;
    let mut center = DVec3::ZERO;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        let tris = face_triangles_pub(brep, *fi);
        for [a, b, c] in &tris {
            let tv = tet_signed_volume(*a, *b, *c);
            total_vol += tv;
            center += (*a + *b + *c) * tv * 0.25;
        }
    }
    if total_vol.abs() > 1e-15 { center / total_vol } else { DVec3::ZERO }
}
