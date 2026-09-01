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

/// OCCT BRepGProp_Gauss::Inertia (BRepGProp_Gauss.cxx L160-172) — the
/// accumulator of the Vinert integrator.  `mass` is the volume, `ix/iy/iz`
/// the first moments about the location point, and `ixx..iyz` the second
/// moments accumulated by computeVInertiaOfElementaryPart (isByPoint,
/// BRepGProp_Gauss.cxx L306-339):
///
///   dv = r . (N * w);            Mass += dv / 3
///   Ix/Iy/Iz += 0.25 * r * dv
///   x -= coeff[0]; y -= coeff[1]; z -= coeff[2];  dv *= 0.2
///   Ixy -= x*y*dv; Iyz -= y*z*dv; Ixz -= x*z*dv
///   x *= x; y *= y; z *= z
///   Ixx += (y+z)*dv; Iyy += (x+z)*dv; Izz += (x+y)*dv
///
/// with the divergence theorem: Ixx = ∫(y²+z²) dV, Ixy = −∫xy dV etc. about
/// the location point.  For the rcad VolumeProperties path the location is
/// the origin and aCoeff = {0,0,0}, so r = P and the accumulated moments are
/// about the world origin.
#[derive(Debug, Clone, Copy, Default)]
pub struct VinertFace {
    pub mass: f64,
    pub ix: f64,
    pub iy: f64,
    pub iz: f64,
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixy: f64,
    pub ixz: f64,
    pub iyz: f64,
}

impl VinertFace {
    /// Element-wise addition (GProp_GProps::Add first branch, L48-64 — all
    /// faces share the location, so the inertia matrices add directly).
    pub fn add(self, other: VinertFace) -> VinertFace {
        VinertFace {
            mass: self.mass + other.mass,
            ix: self.ix + other.ix,
            iy: self.iy + other.iy,
            iz: self.iz + other.iz,
            ixx: self.ixx + other.ixx,
            iyy: self.iyy + other.iyy,
            izz: self.izz + other.izz,
            ixy: self.ixy + other.ixy,
            ixz: self.ixz + other.ixz,
            iyz: self.iyz + other.iyz,
        }
    }

    /// Element-wise subtraction — a REVERSED occurrence flips the face
    /// normal (BRepGProp_Face::Load L196 mySReverse), so every component
    /// changes sign.
    pub fn sub(self, other: VinertFace) -> VinertFace {
        VinertFace {
            mass: self.mass - other.mass,
            ix: self.ix - other.ix,
            iy: self.iy - other.iy,
            iz: self.iz - other.iz,
            ixx: self.ixx - other.ixx,
            iyy: self.iyy - other.iyy,
            izz: self.izz - other.izz,
            ixy: self.ixy - other.ixy,
            ixz: self.ixz - other.ixz,
            iyz: self.iyz - other.iyz,
        }
    }
}

/// OCCT BRepGProp_Gauss::Compute (BRepGProp_Gauss.cxx L1215-1302, Vinert) for
/// one face with wires: line integral over the boundary edge pcurves.
///
/// Returns (mass, first moments) — the Vinert integrator accumulates
/// Mass += dv/3 and Ix/Iy/Iz += 0.25*r*dv per elementary part
/// (computeVInertiaOfElementaryPart isByPoint, BRepGProp_Gauss.cxx L306-339).
///
///   Mass = sum over edges of
///     lr * sum_i { ur * sum_j ( (1/3) r(u_j, v_i) . N(u_j, v_i) * (dv/dl)_i * w_i * w_j ) }
///
/// with r = point - origin (the shape origin, L123-127 of BRepGProp.cxx),
/// N = D1U x D1V (BRepGProp_Face::Normal L201-210, mySReverse=false: rcad
/// faces are explored FORWARD, see face_flat_iter), and the v coordinate
/// clamped to the face UV bounds (OCC104, L1266-1267).
pub fn face_volume_gauss_domain(brep: &BRep, face: &Face, fi: usize) -> (f64, DVec3) {
    let v = face_volume_gauss_domain_full(brep, face, fi);
    (v.mass, DVec3::new(v.ix, v.iy, v.iz))
}

/// Full Vinert integration of one face with wires — see
/// face_volume_gauss_domain for the OCCT mapping; this variant also
/// accumulates the second moments Ixx/Iyy/Izz/Ixy/Ixz/Iyz (the
/// computeVInertiaOfElementaryPart tail, BRepGProp_Gauss.cxx L326-338).
pub fn face_volume_gauss_domain_full(brep: &BRep, face: &Face, fi: usize) -> VinertFace {
    let surf_idx = match brep.tshapes.get(fi).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts {
            fd.surface.clone()
        } else {
            None
        }
    }) {
        Some(s) => s,
        None => return VinertFace::default(),
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

    let mut an_inertia = VinertFace::default();
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
                    None => return VinertFace::default(),
                },
            }
        } else {
            match pc {
                Some(v) => v.clone(),
                None => return VinertFace::default(),
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
            None => return VinertFace::default(),
        };
        // L1250-1253: l1/l2 = the loaded pcurve First/LastParameter.
        let (l1, l2) = (a, b);
        let lm = 0.5 * (l2 + l1);
        let lr = 0.5 * (l2 - l1);

        let mut a_c_inertia = VinertFace::default();
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
            // L1273-1291: aLocal = sum_j over the elementary parts
            // (computeVInertiaOfElementaryPart isByPoint, L306-339) at
            // (u_j, v) with weight Dul * wV_j.
            let mut a_local = VinertFace::default();
            for j in 0..nb_gauss {
                let u = um + ur * cp[j];
                let a_weight = dul * cw[j];
                // BRepGProp_Face::Normal (L201-210): D1U x D1V; mySReverse is
                // false for rcad faces (explored FORWARD).
                let (p, du, dv) = surf.derivatives(u, v);
                let n = du.cross(dv);
                // computeVInertiaOfElementaryPart (L306-339, isByPoint): the
                // reference point is the shape origin — rcad faces are already
                // in world coordinates, so r = P and aCoeff = {0,0,0}.
                let x = p.x;
                let y = p.y;
                let z = p.z;
                let xn = n.x * a_weight;
                let yn = n.y * a_weight;
                let zn = n.z * a_weight;
                let dv = x * xn + y * yn + z * zn;
                a_local.mass += dv / 3.0;
                // L314-316: first moments Ix/Iy/Iz += 0.25 * r * dv.
                a_local.ix += 0.25 * x * dv;
                a_local.iy += 0.25 * y * dv;
                a_local.iz += 0.25 * z * dv;
                // L326-338: second moments (dv *= 0.2, then the squares).
                let dv2 = dv * 0.2;
                a_local.ixy -= x * y * dv2;
                a_local.iyz -= y * z * dv2;
                a_local.ixz -= x * z * dv2;
                let x2 = x * x;
                let y2 = y * y;
                let z2 = z * z;
                a_local.ixx += (y2 + z2) * dv2;
                a_local.iyy += (x2 + z2) * dv2;
                a_local.izz += (x2 + y2) * dv2;
            }
            // L1293: aLocal *= ur
            a_c_inertia.mass += a_local.mass * ur;
            a_c_inertia.ix += a_local.ix * ur;
            a_c_inertia.iy += a_local.iy * ur;
            a_c_inertia.iz += a_local.iz * ur;
            a_c_inertia.ixx += a_local.ixx * ur;
            a_c_inertia.iyy += a_local.iyy * ur;
            a_c_inertia.izz += a_local.izz * ur;
            a_c_inertia.ixy += a_local.ixy * ur;
            a_c_inertia.ixz += a_local.ixz * ur;
            a_c_inertia.iyz += a_local.iyz * ur;
        }
        // L1297-1298: aC *= lr; anInertia += aC
        an_inertia.mass += a_c_inertia.mass * lr;
        an_inertia.ix += a_c_inertia.ix * lr;
        an_inertia.iy += a_c_inertia.iy * lr;
        an_inertia.iz += a_c_inertia.iz * lr;
        an_inertia.ixx += a_c_inertia.ixx * lr;
        an_inertia.iyy += a_c_inertia.iyy * lr;
        an_inertia.izz += a_c_inertia.izz * lr;
        an_inertia.ixy += a_c_inertia.ixy * lr;
        an_inertia.ixz += a_c_inertia.ixz * lr;
        an_inertia.iyz += a_c_inertia.iyz * lr;
    }
    // convert (L467-490): |Mass| >= EPS_DIM(1e-30) → mass else 0 (the
    // moments follow the same guard).
    if an_inertia.mass.abs() >= EPS_DIM {
        an_inertia
    } else {
        VinertFace::default()
    }
}

/// OCCT BRepGProp_Gauss::Compute (BRepGProp_Gauss.cxx L1306-1393, Vinert) for
/// a natural-restriction face (no wires): Gauss-Legendre over the surface
/// natural parameter bounds.
///
/// Returns (mass, first moments) — see face_volume_gauss_domain.
///   Mass = vr * ur * sum_j ( wV_j * sum_i ( (1/3) r(u_i, v_j) . N(u_i, v_j) * wU_i ) )
pub fn face_volume_gauss_natural(brep: &BRep, fi: usize) -> (f64, DVec3) {
    let v = face_volume_gauss_natural_full(brep, fi);
    (v.mass, DVec3::new(v.ix, v.iy, v.iz))
}

/// Full Vinert integration of a natural-restriction face (no wires) — see
/// face_volume_gauss_natural for the OCCT mapping; this variant also
/// accumulates the second moments Ixx/Iyy/Izz/Ixy/Ixz/Iyz.
pub fn face_volume_gauss_natural_full(brep: &BRep, fi: usize) -> VinertFace {
    let surf_idx = match brep.tshapes.get(fi).and_then(|ts| {
        if let topods::TShape::Face(fd) = &**ts {
            fd.surface.clone()
        } else {
            None
        }
    }) {
        Some(s) => s,
        None => return VinertFace::default(),
    };
    let surf = &surf_idx;
    // L1314-1316: theSurface.Bounds — surface natural parameter domain.
    let [lower_u, upper_u, lower_v, upper_v] = surf.default_domain();
    // checkBounds (L418-429): an infinite bound makes the mass 0 (the convert
    // guard, |Mass| >= EPS_DIM is false for NaN).
    if !lower_u.is_finite() || !upper_u.is_finite() || !lower_v.is_finite() || !upper_v.is_finite()
    {
        return VinertFace::default();
    }
    // L1318-1319: UIntegrationOrder / VIntegrationOrder, min GaussPointsMax(61)
    let (uo0, vo0) = surface_u_v_integration_order(surf);
    let uo = uo0.min(61);
    let vo = vo0.min(61);
    let (gpu, gwu) = match occt_gauss(uo) {
        Some(v) => v,
        None => return VinertFace::default(),
    };
    let (gpv, gwv) = match occt_gauss(vo) {
        Some(v) => v,
        None => return VinertFace::default(),
    };
    // L1332-1335
    let um = 0.5 * (upper_u + lower_u);
    let vm = 0.5 * (upper_v + lower_v);
    let ur = 0.5 * (upper_u - lower_u);
    let vr = 0.5 * (upper_v - lower_v);
    // L1341-1374: anInertia = sum_j ( wV_j * sum_i of the elementary parts
    // (computeVInertiaOfElementaryPart isByPoint, L306-339) at (u_i, v_j) ).
    let mut an_inertia = VinertFace::default();
    for j in 0..vo {
        let v = vm + vr * gpv[j];
        let mut an_inertia_of_elementary_part = VinertFace::default();
        for i in 0..uo {
            let u = um + ur * gpu[i];
            let a_weight = gwu[i];
            let (p, du, dv) = surf.derivatives(u, v);
            let n = du.cross(dv);
            let x = p.x;
            let y = p.y;
            let z = p.z;
            let xn = n.x * a_weight;
            let yn = n.y * a_weight;
            let zn = n.z * a_weight;
            let dv = x * xn + y * yn + z * zn;
            an_inertia_of_elementary_part.mass += dv / 3.0;
            // L314-316: first moments Ix/Iy/Iz += 0.25 * r * dv.
            an_inertia_of_elementary_part.ix += 0.25 * x * dv;
            an_inertia_of_elementary_part.iy += 0.25 * y * dv;
            an_inertia_of_elementary_part.iz += 0.25 * z * dv;
            // L326-338: second moments.
            let dv2 = dv * 0.2;
            an_inertia_of_elementary_part.ixy -= x * y * dv2;
            an_inertia_of_elementary_part.iyz -= y * z * dv2;
            an_inertia_of_elementary_part.ixz -= x * z * dv2;
            let x2 = x * x;
            let y2 = y * y;
            let z2 = z * z;
            an_inertia_of_elementary_part.ixx += (y2 + z2) * dv2;
            an_inertia_of_elementary_part.iyy += (x2 + z2) * dv2;
            an_inertia_of_elementary_part.izz += (x2 + y2) * dv2;
        }
        an_inertia.mass += an_inertia_of_elementary_part.mass * gwv[j];
        an_inertia.ix += an_inertia_of_elementary_part.ix * gwv[j];
        an_inertia.iy += an_inertia_of_elementary_part.iy * gwv[j];
        an_inertia.iz += an_inertia_of_elementary_part.iz * gwv[j];
        an_inertia.ixx += an_inertia_of_elementary_part.ixx * gwv[j];
        an_inertia.iyy += an_inertia_of_elementary_part.iyy * gwv[j];
        an_inertia.izz += an_inertia_of_elementary_part.izz * gwv[j];
        an_inertia.ixy += an_inertia_of_elementary_part.ixy * gwv[j];
        an_inertia.ixz += an_inertia_of_elementary_part.ixz * gwv[j];
        an_inertia.iyz += an_inertia_of_elementary_part.iyz * gwv[j];
    }
    // L1375: vr = vr * ur; L1392: every component *= vr (after the convert
    // guard, L467-490: |Mass| >= EPS_DIM).
    let vr = vr * ur;
    if an_inertia.mass.abs() >= EPS_DIM {
        an_inertia.mass *= vr;
        an_inertia.ix *= vr;
        an_inertia.iy *= vr;
        an_inertia.iz *= vr;
        an_inertia.ixx *= vr;
        an_inertia.iyy *= vr;
        an_inertia.izz *= vr;
        an_inertia.ixy *= vr;
        an_inertia.ixz *= vr;
        an_inertia.iyz *= vr;
        an_inertia
    } else {
        VinertFace::default()
    }
}

/// Signed volume of a BRep solid — OCCT BRepGProp::VolumeProperties
/// (BRepGProp.cxx L298-409): every face integrated with BRepGProp_Vinert
/// (Eps = 1.0, fixed-order Gauss), IsNatRestr faces via the whole-surface
/// integral, the others via the boundary line integral.
pub fn signed_volume(brep: &topods::BRep) -> f64 {
    shape_vinert(brep).mass
}

/// OCCT BRepGProp::VolumeProperties per-face accumulation (BRepGProp.cxx
/// L298-409 + GProp_GProps::Add L42-64): every face occurrence is integrated
/// with BRepGProp_Vinert about the shape origin (aCoeff = {0,0,0},
/// isByPoint = true) and accumulated with the cumulative orientation.
///
/// The faces are explored with TopExp_Explorer(S, TopAbs_FACE), i.e. PER
/// OCCURRENCE with the cumulative orientation (solid * shell * face).  A face
/// shared by two adjacent solids (the splitter's internal cell faces) appears
/// TWICE — FORWARD in one cell, REVERSED in the neighbor — and its two
/// contributions cancel.  The flat BRep pool stores each face TShape once;
/// walk the solid->shell->face structure to reproduce the occurrence
/// semantics, flipping the sign of every component for REVERSED cumulative
/// orientation (BRepGProp_Face::Load L196: mySReverse flips the normal,
/// BRepGProp_Face.cxx L206-208).
pub fn shape_vinert(brep: &topods::BRep) -> VinertFace {
    // Face data by flat pool index (explored FORWARD, BRepGProp_Domain L39-42).
    let faces = face_flat_iter(brep);
    let face_by_idx: std::collections::HashMap<usize, &Face> =
        faces.iter().map(|(i, f)| (*i, f)).collect();
    let integrate = |fi: usize, face: &Face| -> VinertFace {
        let is_nat_restr = face.outer_wire.edges.is_empty() && face.inner_wires.is_empty();
        if is_nat_restr {
            face_volume_gauss_natural_full(brep, fi)
        } else {
            face_volume_gauss_domain_full(brep, face, fi)
        }
    };
    let mut total = VinertFace::default();
    // 1. Faces inside solids: per occurrence with cumulative orientation.
    let mut referenced: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        let fi = face_sr.index;
                        referenced.insert(fi);
                        // OCCT TopExp_Explorer cumulative orientation:
                        // solid(Forward) * shell * face.
                        let ori = shell_sr.orientation.compose(face_sr.orientation);
                        if let Some(face) = face_by_idx.get(&fi) {
                            let v = integrate(fi, face);
                            if ori.is_reversed() {
                                total = total.sub(v);
                            } else {
                                total = total.add(v);
                            }
                        }
                    }
                }
            }
        }
    }
    // 2. Free shells (not inside any solid): per face occurrence, orientation
    // composed from the shell (BRepGProp.cxx L528-549 free-face compound).
    let mut shell_refs: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for sr in &sd.shells {
                shell_refs.insert(sr.index);
            }
        }
    }
    for (ti, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Shell(shd) = &**ts {
            if shell_refs.contains(&ti) {
                continue;
            }
            for face_sr in &shd.faces {
                let fi = face_sr.index;
                referenced.insert(fi);
                let ori = face_sr.orientation;
                if let Some(face) = face_by_idx.get(&fi) {
                    let v = integrate(fi, face);
                    if ori.is_reversed() {
                        total = total.sub(v);
                    } else {
                        total = total.add(v);
                    }
                }
            }
        }
    }
    // 3. Standalone faces not referenced by any shell: counted once FORWARD
    // (the flat pool has no orientation context for them).
    for (fi, face) in &faces {
        if !referenced.contains(fi) {
            total = total.add(integrate(*fi, face));
        }
    }
    total
}

/// Absolute volume of a BRep solid.
pub fn volume(brep: &topods::BRep) -> f64 {
    signed_volume(brep).abs()
}

/// Centroid of a BRep solid — OCCT BRepGProp::VolumeProperties
/// (BRepGProp.cxx L298-409): the first moments Ix/Iy/Iz from the same Vinert
/// line/domain integrals as the volume (computeVInertiaOfElementaryPart
/// isByPoint, BRepGProp_Gauss.cxx L314-316), divided by the mass.  Per-face
/// cumulative orientation (REVERSED faces flip both mass and moment signs,
/// so the ratio stays consistent).
pub fn centroid(brep: &topods::BRep) -> DVec3 {
    let v = shape_vinert(brep);
    if v.mass.abs() > 1e-15 {
        DVec3::new(v.ix / v.mass, v.iy / v.mass, v.iz / v.mass)
    } else {
        DVec3::ZERO
    }
}
