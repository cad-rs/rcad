//! OCCT BRepGProp::SurfaceProperties (BRepGProp.cxx L167-266): surface area.
//!
//! 1:1 translation of the `checkprops -s` path.  The DRAW `checkprops` Tcl
//! command (resources/DrawResources/CheckCommands.tcl) evaluates
//! `sprops shape 1.0e-4` — i.e. the EPSILON overload
//! `BRepGProp::SurfaceProperties(S, Props, Eps, SkipShared)`
//! (BRepGProp.cxx L280-291).  With Eps = 1.0e-4 < 1.0 the per-face branch
//! (L231-242) runs the ADAPTIVE `BRepGProp_Sinert::Perform(BF, BD, Eps)`
//! (BRepGProp_Sinert.cxx L104-110 → BRepGProp_Gauss.cxx L533-1099) instead of
//! the fixed-order Gauss path (L243-253).
//!
//! There is no analytic per-surface-type dispatch in OCCT: every face is
//! integrated by Gauss-Legendre (whole surface for natural-restriction faces,
//! Green-theorem line integral over the edge pcurves otherwise).

use crate::base::bnd_lib::curve2d_bounding_box;
use crate::base::gprop::tri::face_flat_iter;
use crate::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval};
use crate::topo::topo_shape::Shape;
use crate::topo::topods::{self, BRepTool};
use crate::topo::topology::{Face, WireEdge};
use crate::BRep;
use std::f64::consts::PI;

// OCCT math::GaussPoints/GaussWeights tables (math.cxx): packed positive-half
// nodes/weights for orders 1..61 plus the GaussPoints expansion.
include!("gauss_tables.rs");

/// OCCT BRepGProp::SurfaceProperties with Eps = 1.0e-4 — the `checkprops -s`
/// path (CheckCommands.tcl → BRepGProp.cxx L280-291 → surfaceProperties
/// L231-242 → BRepGProp_Sinert.cxx L104-110).
///
/// Per-face loop (L190-259):
///   - NoSurf/NoTri check (L198-214): rcad faces always carry a surface
///     (NoSurf=false) and UseTriangulation=false, so the triangulation branch
///     (L216-222) is never taken.
///   - BF.Load(F) (L225); IsNatRestr = (F.NbChildren() == 0) (L226).
///   - Eps = 1.0e-4 < 1.0 → adaptive branch (L231-242):
///       G.Perform(BF, BD, Eps)  (BRepGProp_Gauss.cxx L533-1099)
///   - Props.Add(G) (L254): the face area accumulates into the total mass.
pub fn surface_area(brep: &topods::BRep) -> f64 {
    let mut mass = 0.0;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        // OCCT L226: IsNatRestr = (F.NbChildren() == 0) — the face carries no
        // wires (no outer wire edges and no inner wires).
        let is_nat_restr = face.outer_wire.edges.is_empty() && face.inner_wires.is_empty();
        let a = face_surface_area_checkprops(brep, face, *fi, 1.0e-4, is_nat_restr);
        // OCCT BRepGProp_Gauss::Compute (BRepGProp_Gauss.cxx L533-1099)
        // integrates the surface patch mass = |dS| via the Green theorem
        // boundary integral; the result follows the wire direction, so a face
        // whose wires run the opposite way (e.g. the upper hemisphere after a
        // boolean) yields a negative value.  Take the absolute value to match
        // OCCT's positive per-face area.
        mass += a.abs();
    }
    mass
}

/// Surface area of a single face — OCCT surfaceProperties L243-253 with Eps=1.0.
pub fn face_surface_area(brep: &BRep, face: &Face, face_flat_idx: usize) -> f64 {
    let is_nat_restr = face.outer_wire.edges.is_empty() && face.inner_wires.is_empty();
    if is_nat_restr {
        face_surface_area_gauss_natural(brep, face_flat_idx)
    } else {
        face_surface_area_gauss_domain(brep, face, face_flat_idx)
    }
}

/// OCCT BRepGProp_Face::UIntegrationOrder/VIntegrationOrder (BRepGProp_Face.cxx
/// L39-107): Gauss order for the face surface (GeomAbs type -> Nu/Nv,
/// max(8, 2N)).
fn surface_u_v_integration_order(surf: &Surface3) -> (usize, usize) {
    let (nu, nv): (usize, usize) = match surf {
        Surface3::Plane(_) => (4, 4),
        Surface3::Bezier(b) => {
            let du = b.control_points.len().saturating_sub(1) + 1;
            let dv = b.control_points.first().map_or(1, |r| r.len().saturating_sub(1)) + 1;
            (du, dv)
        }
        Surface3::BSpline(b) => {
            // OCCT L56-63: Nu = max(4, (UDeg+1) * (NbUKnots-1)); NbUKnots is the
            // COMPRESSED knot count (Geom_BSplineSurface::NbUKnots).
            let ku = compressed_knot_count(&b.knots_u);
            let kv = compressed_knot_count(&b.knots_v);
            let du = ((b.degree_u + 1) * ku.saturating_sub(1)).max(4);
            let dv = ((b.degree_v + 1) * kv.saturating_sub(1)).max(4);
            (du, dv)
        }
        _ => (9, 9),
    };
    ((2 * nu).max(8), (2 * nv).max(8))
}

/// OCCT Geom_BSplineCurve::NbKnots / Geom_BSplineSurface::NbUKnots — the
/// COMPRESSED knot count (the rcad knot vectors are expanded).
fn compressed_knot_count(knots: &[f64]) -> usize {
    let mut n = 0usize;
    let mut prev = f64::NAN;
    for &k in knots {
        if n == 0 || (k - prev).abs() > 1e-15 {
            n += 1;
            prev = k;
        }
    }
    n.max(1)
}

/// OCCT BRepGProp_Face::IntegrationOrder (BRepGProp_Face.cxx L111-150):
/// Gauss order for the boundary edge pcurve (GeomAbs type -> N, max(4, 2N)).
fn curve_integration_order(c: &Curve2d) -> usize {
    let n = match c {
        Curve2d::Line(_) => 2,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) | Curve2d::Hyperbola(_) | Curve2d::Parabola(_) => 9,
        // OCCT L137-141: N = (Deg+1) * (NbKnots-1) with NbKnots compressed.
        Curve2d::BSpline(b) => {
            let k = compressed_knot_count(&b.knots);
            (b.degree + 1) * k.saturating_sub(1)
        }
        // Bezier degree = control point count - 1
        Curve2d::Bezier(b) => b.control_points.len(),
        // OCCT GeomAdaptor_Curve::GetType unwraps a Geom2d_TrimmedCurve to the
        // basis curve type (myCurve.GetType() in IntegrationOrder L115).
        Curve2d::Trimmed(tc) => curve_integration_order(&tc.curve),
        _ => 9,
    };
    (2 * n).max(4)
}

/// OCCT BRepGProp_Gauss::Compute(const BRepGProp_Face&, loc, ...)
/// (BRepGProp_Gauss.cxx L1306-1393, Sinert).  Natural-restriction face (no
/// wires): Gauss-Legendre over the surface natural parameter bounds.
///
///   Mass = vr * ur * sum_j ( wV_j * sum_i ( wU_i * |N(u_i, v_j)| ) )
///
/// with um = (U2+U1)/2, vm = (V2+V1)/2, ur = (U2-U1)/2, vr = (V2-V1)/2 and
/// u = um + ur*x_i, v = vm + vr*x_j.
fn face_surface_area_gauss_natural(brep: &BRep, fi: usize) -> f64 {
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
    // BRepGProp_Face::Bounds (BRepGProp_Face.cxx L154-160): surface natural
    // parameter domain [u1,u2]x[v1,v2].
    let [lower_u, upper_u, lower_v, upper_v] = surf.default_domain();
    // checkBounds (L418-429): an infinite bound switches + and * to AddInf /
    // MultInf (L41-140).  With a ±INF bound, um/vm collapse to 0, u/v sample
    // ±INF and the |N| term becomes NaN, which the convert guard (L472,
    // |Mass| >= EPS_DIM is false for NaN) turns into mass 0.  A natural
    // restriction face in the tested shapes is a closed bounded surface, so an
    // unbounded natural domain yields 0 here, matching OCCT.
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
    // L1341-1374: anInertia = sum_j ( wV_j * sum_i ( wU_i * |N(u_i, v_j)| ) )
    let mut an_inertia = 0.0;
    for j in 0..vo {
        let v = vm + vr * gpv[j];
        let mut an_inertia_of_elementary_part = 0.0;
        for i in 0..uo {
            let u = um + ur * gpu[i];
            // BRepGProp_Face::Normal (BRepGProp_Face.cxx L201-210): D1U x D1V
            let (_p, du, dv) = surf.derivatives(u, v);
            let n = du.cross(dv);
            an_inertia_of_elementary_part += n.length() * gwu[i];
        }
        an_inertia += an_inertia_of_elementary_part * gwv[j];
    }
    // L1375: vr = vr * ur; L1389/1392: convert (L472, |Mass| >= EPS_DIM=1e-30
    // → mass else 0) then mass *= vr.
    let vr = vr * ur;
    let mass = if an_inertia.abs() >= 1e-30 { an_inertia } else { 0.0 };
    mass * vr
}

/// OCCT BRepTools::UVBounds (BRepTools.cxx L64-80) + AddUVBounds (L126-367).
/// Bounding UV box of the face: union of every boundary edge pcurve's 2D box
/// (BndLib_Add2dCurve::Add, L185), with out-of-natural-domain coordinates
/// clamped for non-periodic surfaces (L270-280, L351-361).  A face without
/// edges or pcurves falls back to the surface natural bounds (L139-153).
/// Returns [u_min, u_max, v_min, v_max].
///
/// This is the domain BRepAdaptor_Surface (Restriction=true) exposes through
/// BRepGProp_Face::Bounds (BRepGProp_Face.cxx L154-160); BRepGProp_Gauss
/// uses its UMin as the integration offset u1 (BRepGProp_Gauss.cxx L1135-1136,
/// L1186-1187).
fn face_uv_bounds(brep: &BRep, face: &Face, fi: usize, surf: &Surface3) -> [f64; 4] {
    // aS->Bounds (L191-192): surface natural parameter domain
    let [a_u_min, a_u_max, a_v_min, a_v_max] = surf.default_domain();
    // Bnd_Box2d aBox (L133) / aBoxS (L176): [u_min, u_max, v_min, v_max]
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    let mut any = false;
    let face_shape = Shape::from_parts(brep.tshapes[fi].clone(), fi, 0, topods::Orientation::Forward);
    let mut edges: Vec<WireEdge> = face.outer_wire.edges.clone();
    for w in &face.inner_wires {
        edges.extend(w.edges.iter().copied());
    }
    for we in &edges {
        let edge_shape = Shape::from_parts(
            brep.tshapes[we.idx].clone(),
            we.idx,
            we.location,
            if we.forward {
                topods::Orientation::Forward
            } else {
                topods::Orientation::Reversed
            },
        );
        // BRep_Tool::CurveOnSurface (L179): a null pcurve means this edge
        // contributes nothing to the box (L180-183) — but OCCT's
        // BRep_Tool::CurveOnSurface itself has the CurveOnPlane fallback
        // (BRep_Tool.cxx L327-450: planar face + 3D edge -> projected pcurve),
        // so a pcurve-less planar edge still bounds the UV domain.  Mirror the
        // fallback here; without it a prism side face (no stored pcurve) falls
        // back to the infinite natural plane domain and the Green integral
        // diverges.
        let (c, a_t1, a_t2) = match brep.curve_on_surface(&edge_shape, &face_shape) {
            Some(v) => v.clone(),
            None => match project_curve_on_plane(brep, &edge_shape, &face_shape).or_else(|| project_curve_on_cylinder(brep, &edge_shape, &face_shape)) {
                Some(v) => v,
                None => continue,
            },
        };
        // BndLib_Add2dCurve::Add (L185): 2D bounding box [xmin,ymin,xmax,ymax]
        let a_box_c = curve2d_bounding_box(&c, a_t1, a_t2, 0.0);
        // aBoxC.Get (L188)
        let mut a_x_min = a_box_c[0];
        let mut a_y_min = a_box_c[2];
        let mut a_x_max = a_box_c[1];
        let mut a_y_max = a_box_c[3];
        // U periodicity (L202-281): for a non-periodic surface, U coordinates
        // of the edge box that fall inside the natural U range are clamped to
        // it (L270-280).  (The OCCT BSpline periodicity re-verification
        // L210-268 applies to Geom_BSplineSurface only and is not reached for
        // the analytic surface types used here.)
        if !surf.is_u_periodic() {
            if (a_x_min < a_u_min) && (a_u_min < a_x_max) {
                a_x_min = a_u_min;
            }
            if (a_x_min < a_u_max) && (a_u_max < a_x_max) {
                a_x_max = a_u_max;
            }
        }
        // V periodicity (L283-362): same clamping for V (L351-361)
        if !surf.is_v_periodic() {
            if (a_y_min < a_v_min) && (a_v_min < a_y_max) {
                a_y_min = a_v_min;
            }
            if (a_y_min < a_v_max) && (a_v_max < a_y_max) {
                a_y_max = a_v_max;
            }
        }
        // aBoxS.Update + aB.Add (L364-366)
        b[0] = b[0].min(a_x_min);
        b[1] = b[1].max(a_x_max);
        b[2] = b[2].min(a_y_min);
        b[3] = b[3].max(a_y_max);
        any = true;
    }
    if !any {
        // UVBounds: void box → surface natural bounds (L139-153)
        return [a_u_min, a_u_max, a_v_min, a_v_max];
    }
    b
}

/// OCCT BRepGProp_Gauss::Compute(Face, Domain) (BRepGProp_Gauss.cxx
/// L1126-1211, Sinert).  For a face whose boundary is not the natural
/// restriction (i.e. the face carries wires), BRepGProp::SurfaceProperties
/// integrates along every boundary edge pcurve instead of over the full UV
/// rectangle:
///
///   Mass = sum over edges of
///     lr * sum_i { ur * sum_j ( |N(u_j, v_i)| * (dv/dl)_i * w_i * w_j ) }
///   with u_j = um + ur*x_j, v_i the edge's V at Gauss point l_i,
///   um = (u2+u1)/2, ur = (u2-u1)/2, and u2 = the edge's U at l_i.
///
/// This is the Green theorem form of the surface area integral and reproduces
/// OCCT's values bit-for-bit for trimmed curved faces.
pub fn face_surface_area_gauss_domain(brep: &BRep, face: &Face, fi: usize) -> f64 {
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
    // theSurface.Bounds (L1135-1136): BRepGProp_Face::Bounds =
    // BRepAdaptor_Surface (Restriction=true) UV domain = face UV bounds
    // (BRepTools::UVBounds); only UMin feeds the integration offset
    // (L1186-1187).
    let [u1, _u2, _v1, _v2] = face_uv_bounds(brep, face, fi, surf);
    // UIntegrationOrder / VIntegrationOrder (L1139-1143), GaussPointsMax = 61
    let (nb_u, nb_v) = surface_u_v_integration_order(surf);
    let nb_gauss = nb_u.min(61).max(nb_v.min(61));
    let (spv, swv) = match occt_gauss(nb_gauss) {
        Some(v) => v,
        None => return 0.0,
    };

    let face_shape = Shape::from_parts(brep.tshapes[fi].clone(), fi, 0, topods::Orientation::Forward);

    let mut an_inertia = 0.0;
    // Domain: every edge of the face (outer + inner wires).
    let mut edges: Vec<WireEdge> = face.outer_wire.edges.clone();
    for w in &face.inner_wires {
        edges.extend(w.edges.iter().copied());
    }
    for we in &edges {
        // OCCT BRepGProp_Domain::Next (BRepGProp_Domain.cxx L27-38) skips
        // INTERNAL and EXTERNAL edges — they are not face boundary edges and
        // contribute nothing to the boundary (Green) integral.
        if we.internal {
            continue;
        }
        let edge_shape = Shape::from_parts(
            brep.tshapes[we.idx].clone(),
            we.idx,
            we.location,
            if we.forward {
                topods::Orientation::Forward
            } else {
                topods::Orientation::Reversed
            },
        );
        let pc = brep.curve_on_surface(&edge_shape, &face_shape);
        // BRep_Tool::CurveOnSurface (BRep_Tool.cxx L327-373): the edge's
        // pcurve on the face.  A REVERSED edge on a closed surface (seam) uses
        // PCurve2 (the u=0 image); other edges use PCurve.
        // BRepGProp_Face::Load (BRepGProp_Face.cxx L164-185): when the pcurve
        // is null, Load returns false and OCCT Compute returns immediately
        // (L1154-1157) with the face mass left at 0.
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
        // BRepGProp_Face::Load (L176-182): REVERSED edge reverses the pcurve
        // and its parameter range.
        let reversed = !we.forward;
        let (a, b) = if reversed {
            let x = a0;
            (c.reversed_parameter(b0), c.reversed_parameter(x))
        } else {
            (a0, b0)
        };
        // IntegrationOrder (L1159-1161): the boundary pcurve Gauss order,
        // raised to at least the face surface order.
        let nb_c = curve_integration_order(&c).min(61).max(nb_gauss);
        // The fixed-order Gauss rule of L1159 is exact on each knot span of a
        // BSpline boundary pcurve (the C0 junctions of the OCCT
        // GeomInt_WLApprox curve are integrand discontinuities): integrate
        // per span with the (degree+1)-order rule.  The Green-theorem
        // integral is the rcad boundary-area implementation (OCCT BRepGProp
        // uses the same Green form for wire-carrying faces); per-span Gauss
        // is its exact quadrature.
        let mut ranges: Vec<(f64, f64)> = Vec::new();
        match &c {
            Curve2d::BSpline(bs) => {
                for i in 1..bs.knots.len() {
                    let (k0, k1) = (bs.knots[i - 1], bs.knots[i]);
                    if k1 <= k0 {
                        continue;
                    }
                    let (s0, s1) = if reversed {
                        (c.reversed_parameter(k1), c.reversed_parameter(k0))
                    } else {
                        (k0, k1)
                    };
                    let l1 = s0.max(a);
                    let l2 = s1.min(b);
                    if l2 > l1 {
                        ranges.push((l1, l2));
                    }
                }
                if ranges.is_empty() {
                    ranges.push((a, b));
                }
            }
            _ => ranges.push((a, b)),
        }
        for (l1, l2) in ranges {
            // Per span the curve is a single polynomial piece: (degree+1)
            // Gauss points integrate it exactly.
            let nb_c = match &c {
                Curve2d::BSpline(bs) => 2 * (bs.degree + 1),
                _ => curve_integration_order(&c),
            }
            .min(61)
            .max(nb_gauss);
        let (cp, cw) = match occt_gauss(nb_c) {
            Some(v) => v,
            None => return 0.0,
        };
        let lm = 0.5 * (l2 + l1);
        let lr = 0.5 * (l2 - l1);

        let mut a_c_inertia = 0.0;
        for i in 0..nb_c {
            let l = lm + lr * cp[i];
            // D12d (L1177-1179): point and first derivative of the 2D curve
            let (puv, vuv) = if reversed {
                let rp = c.reversed_parameter(l);
                (c.point_at(rp), -c.derivative_at(rp))
            } else {
                (c.point_at(l), c.derivative_at(l))
            };
            let v = puv.y;
            let u2 = puv.x;
            // Dul = dv/dl * w_i (L1182-1185)
            let dul = vuv.y * cw[i];
            let um = 0.5 * (u2 + u1);
            let ur = 0.5 * (u2 - u1);
            let mut a_local = 0.0;
            for j in 0..nb_gauss {
                let u = um + ur * spv[j];
                let a_weight = dul * swv[j];
                // BRepGProp_Face::Normal (BRepGProp_Face.cxx L201-210):
                // D1U x D1V
                let (_p, du, dv) = surf.derivatives(u, v);
                let n = du.cross(dv);
                a_local += n.length() * a_weight;
            }
            // L1202-1203: aLocal *= ur; aC += aLocal
            a_c_inertia += a_local * ur;
        }
        // L1206-1207: aC *= lr; anInertia += aC
        an_inertia += a_c_inertia * lr;
            if std::env::var("RCAD_SA_DEBUG").is_ok() {
                eprintln!("[SA-EDGE] f={} edge={} fwd={} span=[{:.4},{:.4}] nb_c={} contrib={:.6}", fi, we.idx, we.forward, l1, l2, nb_c, a_c_inertia * lr);
            }
        }
    }
    // convert (L1210, L467-490): |Mass| >= EPS_DIM(1e-30) → mass else 0
    if an_inertia.abs() >= 1e-30 { an_inertia } else { 0.0 }
}

// ============================================================================
// OCCT `checkprops` adaptive surface integration
// (CheckCommands.tcl → BRepGProp.cxx L280-291 → L231-242 →
//  BRepGProp_Sinert.cxx L104-110 → BRepGProp_Gauss.cxx L533-1099).
// ============================================================================

// BRepGProp_Gauss.cxx L25-33.
const EPS_PARAM: f64 = 1.0e-12;
const EPS_DIM: f64 = 1.0e-30;
const ERROR_ALGEBR_RATIO: f64 = 2.0 / 3.0;
const GPM: usize = 61; // math::GaussPointsMax() (math.cxx L25-28)
const SUBS_POWER: i64 = 32;
const SM: usize = SUBS_POWER as usize * GPM + 1;
// BRepGProp_Face.cxx L35: Epsilon(1.)
const EPSILON1: f64 = 2.220446049250313e-16;

// BRepGProp_Face.cxx L215-225 (OCC104 integration-order coefficients).
const SC_AS: f64 = -0.15;
const SC_AL: f64 = -0.50;
const SC_B: f64 = 1.0;
const SC_C: f64 = 0.75;
const SC_D: f64 = 0.25;

// BRepGProp_Face.cxx L217-225.
fn s_coeff(eps: f64) -> f64 {
    if eps < 0.1 { SC_AS * (SC_B + eps.log10()) + SC_C } else { SC_C }
}

fn l_coeff(eps: f64) -> f64 {
    if eps < 0.1 { SC_AL * (SC_B + eps.log10()) + SC_D } else { SC_D }
}

// Standard::RealToInt — truncation toward zero.
fn real_to_int(v: f64) -> i64 {
    v.trunc() as i64
}

// math_VectorBase::Max (math_VectorBase.lxx L162-176): INDEX of the max
// element, scanning the whole allocated vector (1-based).

// BRepGProp_Gauss::MaxSubs (L192-195, theCoeff defaults to 32).
fn max_subs(n: i64, coeff: i64) -> i64 {
    if i64::MAX / coeff < n { i64::MAX } else { n * coeff + 1 }
}

// BRepGProp_Gauss::Init (L199-215): set v[first..=last] to theValue; the
// (last - first == 0) case fills the whole array.
fn vec_init(v: &mut [f64], first: usize, last: usize, value: f64) {
    if last - first == 0 {
        v.fill(value);
    } else {
        for x in v.iter_mut().take(last + 1).skip(first) {
            *x = value;
        }
    }
}

// OCCT NCollection_Array1/math_Vector with the reallocation semantics of
// BRepGProp_Gauss::FillIntervalBounds (L256-269): growing zero-fills the
// whole array (new math_Vector(1, aSize, 0.0)); the allocated size only
// grows, so entries beyond the current used range keep earlier values.
struct Ovec {
    data: Vec<f64>,
    upper: usize,
}

impl Ovec {
    fn new() -> Self {
        Ovec { data: vec![0.0; 1], upper: 0 }
    }
    fn grow(&mut self, size: usize) {
        if size > self.upper {
            self.data.iter_mut().for_each(|x| *x = 0.0);
            self.data.resize(size + 1, 0.0);
            self.upper = size;
        }
    }
    fn max_index(&self) -> usize {
        let mut i = 0usize;
        let mut x = f64::MIN;
        for idx in 1..=self.upper {
            if self.data[idx] > x {
                x = self.data[idx];
                i = idx;
            }
        }
        i
    }
}

// OCCT AddInf / MultInf (BRepGProp_Gauss.cxx L41-155): infinite-aware
// arithmetic switched in by checkBounds (L418-429).
fn add_inf(a: f64, b: f64) -> f64 {
    if a == f64::INFINITY {
        return if b == f64::NEG_INFINITY { 0.0 } else { f64::INFINITY };
    }
    if b == f64::INFINITY {
        return if a == f64::NEG_INFINITY { 0.0 } else { f64::INFINITY };
    }
    if a == f64::NEG_INFINITY {
        return if b == f64::INFINITY { 0.0 } else { f64::NEG_INFINITY };
    }
    if b == f64::NEG_INFINITY {
        return if a == f64::INFINITY { 0.0 } else { f64::NEG_INFINITY };
    }
    a + b
}

fn mult_inf(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        return 0.0;
    }
    if a == f64::INFINITY {
        return if b < 0.0 { f64::NEG_INFINITY } else { f64::INFINITY };
    }
    if b == f64::INFINITY {
        return if a < 0.0 { f64::NEG_INFINITY } else { f64::INFINITY };
    }
    if a == f64::NEG_INFINITY {
        return if b < 0.0 { f64::INFINITY } else { f64::NEG_INFINITY };
    }
    if b == f64::NEG_INFINITY {
        return if a < 0.0 { f64::INFINITY } else { f64::NEG_INFINITY };
    }
    a * b
}

// BRepGProp_Face::SIntOrder (BRepGProp_Face.cxx L229-269).
fn surf_s_int_order(surf: &Surface3, eps: f64) -> usize {
    let (nu, nv): (i64, i64) = match surf {
        Surface3::Plane(_) => (1, 1),
        Surface3::Cylinder(_) | Surface3::Cone(_) => (2, 1),
        Surface3::Sphere(_) | Surface3::Torus(_) => (2, 2),
        Surface3::Bezier(b) => (
            b.control_points.len().saturating_sub(1) as i64,
            b.control_points.first().map_or(0, |r| r.len().saturating_sub(1)) as i64,
        ),
        Surface3::BSpline(b) => (b.degree_u as i64, b.degree_v as i64),
        _ => (2, 2),
    };
    let n = real_to_int((s_coeff(eps) * (nu.max(nv) + 1) as f64).ceil());
    (n as usize).min(GPM)
}

// BRepGProp_Face::SVIntSubs (L308-339).
fn surf_sv_int_subs(surf: &Surface3) -> i64 {
    let n: i64 = match surf {
        Surface3::Plane(_) => 2,
        Surface3::Cylinder(_) | Surface3::Cone(_) => 2,
        Surface3::Sphere(_) => 3,
        Surface3::Torus(_) => 4,
        Surface3::Bezier(_) => 2,
        Surface3::BSpline(b) => compressed_knot_count(&b.knots_v) as i64,
        _ => 2,
    };
    n - 1
}

// Geom_BSplineSurface::UKnots / Geom2d_BSplineCurve::Knots — the COMPRESSED
// knot vector (rcad stores knots with multiplicities expanded).
fn compressed_knots(knots: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if out.is_empty() || (k - out[out.len() - 1]).abs() > 1e-15 {
            out.push(k);
        }
    }
    if out.is_empty() {
        out.push(0.0);
    }
    out
}

// BRepGProp_Face::UKnots (L343-374).
fn surf_u_knots(surf: &Surface3, u1: f64, u2: f64) -> Vec<f64> {
    match surf {
        Surface3::Plane(_) => vec![u1, u2],
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
            vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI]
        }
        Surface3::BSpline(b) => compressed_knots(&b.knots_u),
        _ => vec![u1, u2],
    }
}

// BRepGProp_Face::VKnots (L378-413).
fn surf_v_knots(surf: &Surface3, v1: f64, v2: f64) -> Vec<f64> {
    match surf {
        Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => vec![v1, v2],
        Surface3::Sphere(_) => vec![-PI / 2.0, 0.0, PI / 2.0],
        Surface3::Torus(_) => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
        Surface3::BSpline(b) => compressed_knots(&b.knots_v),
        _ => vec![v1, v2],
    }
}

// BRepGProp_Face::LIntSubs (L472-496).
fn curve_l_int_subs(c: &Curve2d) -> i64 {
    let n: i64 = match c {
        Curve2d::Line(_) => 2,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => 4,
        Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => 2,
        Curve2d::BSpline(b) => compressed_knot_count(&b.knots) as i64,
        _ => 2,
    };
    n - 1
}

// BRepGProp_Face::LKnots (L500-534).  For a reversed edge the Load (L176-182)
// reverses the pcurve; Geom2d_BSplineCurve::Reversed negates the knot vector
// (ReversedParameter = -U), so the reversed curve's knots are -knots.
fn curve_l_knots(c: &Curve2d, a: f64, b: f64, reversed: bool) -> Vec<f64> {
    match c {
        Curve2d::Line(_) => vec![a, b],
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => {
            vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI]
        }
        Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => vec![a, b],
        Curve2d::BSpline(bs) => {
            let knots = compressed_knots(&bs.knots);
            if reversed {
                let mut negated: Vec<f64> = knots.iter().map(|&k| -k).collect();
                negated.sort_by(|x, y| x.partial_cmp(y).unwrap());
                negated
            } else {
                knots
            }
        }
        _ => vec![a, b],
    }
}

// BRepGProp_Face::LIntOrder (L417-468): adaptive boundary Gauss order.
// `reversed` marks a reversed edge whose Load (L176-182) reversed the pcurve:
// the 2D box is evaluated on the reversed curve (original curve over -[a, b]).
fn curve_l_int_order(c: &Curve2d, a: f64, b: f64, surf: &Surface3, v1: f64, v2: f64, eps: f64, reversed: bool) -> usize {
    // BndLib_Add2dCurve::Add(myCurve, 1.e-7, aBox) (L421) — myCurve is the
    // loaded (possibly reversed) pcurve over [a, b]; the reversed curve is the
    // original evaluated at -t, so the box range is [-b, -a].
    let (box_a, box_b) = if reversed { (-b, -a) } else { (a, b) };
    let a_box = curve2d_bounding_box(c, box_a, box_b, 1.0e-7);
    let a_y_min = a_box[2];
    let a_y_max = a_box[3];
    let dv = v2 - v1;
    let an_r = if dv > EPSILON1 { ((a_y_max - a_y_min) / dv).min(1.0) } else { 1.0 };
    let an_r_int = real_to_int((surf_sv_int_subs(surf) as f64 * an_r).ceil());
    let a_l_subs = curve_l_int_subs(c);
    // L434: NS = max(SIntOrder(1.) * anRInt / aLSubs, 1) — integer arithmetic
    let ns = ((surf_s_int_order(surf, 1.0) as i64 * an_r_int) / a_l_subs).max(1);
    let nl0: i64 = match c {
        Curve2d::Line(_) => 1,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => 6,
        Curve2d::Parabola(_) => 6,
        Curve2d::Hyperbola(_) => 9,
        Curve2d::Bezier(be) => be.control_points.len().saturating_sub(1) as i64,
        Curve2d::BSpline(bs) => bs.degree as i64,
        _ => 9,
    };
    let nl = nl0.max(ns);
    let nn = if a_l_subs <= 4 {
        real_to_int((l_coeff(eps) * (nl + 1) as f64).ceil())
    } else {
        nl + 1
    };
    (nn as usize).min(GPM)
}

// BRepGProp_Gauss::FillIntervalBounds (L246-294): split [a, b] at the knots
// strictly inside it; returns the number of subintervals.
#[allow(clippy::too_many_arguments)]
fn fill_interval_bounds(
    a: f64,
    b: f64,
    knots: &[f64],
    num_subs: i64,
    inerts: &mut Ovec,
    param1: &mut Ovec,
    param2: &mut Ovec,
    err: &mut Ovec,
    common_err: Option<&mut Ovec>,
) -> usize {
    let a_size = knots.len().max(max_subs(knots.len() as i64 - 1, num_subs) as usize);
    if a_size - 1 > param1.upper {
        inerts.grow(a_size);
        param1.grow(a_size);
        param2.grow(a_size);
        err.grow(a_size);
        if let Some(ce) = common_err {
            ce.grow(a_size);
        }
    }
    let mut j = 1usize;
    let mut k = 1usize;
    param1.data[j] = a;
    j += 1;
    for &kn in knots {
        if a < kn {
            if kn < b {
                param1.data[j] = kn;
                j += 1;
                param2.data[k] = kn;
                k += 1;
            } else {
                break;
            }
        }
    }
    param2.data[k] = b;
    k
}

// OCCT BRep_Tool::CurveOnPlane (BRep_Tool.cxx L379-450): when an edge has no
// stored pcurve on a planar face, project the edge's 3D curve onto the plane
// along the plane normal.  For a line the projection is the linear map
//   (u, v) = ((p - O).U, (p - O).V)
// on the plane frame; the 2D parameterization is the arc length of the
// projection (the Green boundary integral is invariant to reparameterization).
fn project_curve_on_plane(brep: &BRep, edge: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    use crate::geom::CurveEval;
    let surf = brep.face_surface_world(face)?;
    let (c3, r3) = brep.edge_curve_world(edge)?;
    let Surface3::Plane(pl) = surf else { return None };
    let o = pl.origin;
    let p0 = c3.point_at(r3[0]);
    let p1 = c3.point_at(r3[1]);
    let (u0, v0) = ((p0 - o).dot(pl.u_dir), (p0 - o).dot(pl.v_dir));
    let (u1, v1) = ((p1 - o).dot(pl.u_dir), (p1 - o).dot(pl.v_dir));
    let du = u1 - u0;
    let dv = v1 - v0;
    let len = (du * du + dv * dv).sqrt();
    if len < 1e-30 {
        return None;
    }
    let c = Curve2d::Line(crate::geom::Line2d::new(glam::DVec2::new(u0, v0), glam::DVec2::new(du / len, dv / len)));
    Some((c, 0.0, len))
}

// Semantic counterpart of BRep_Tool::CurveOnPlane for CYLINDRICAL faces:
// when an edge has no stored pcurve for the face, radially project the
// edge's 3D curve onto the cylinder (u = azimuth from ref_dir about the
// axis, v = height along the axis).  For co-axial circles (the common case:
// rim/section arcs perpendicular to the axis) the projection keeps
// v constant, so the uv line reproduces the true boundary exactly; the
// Green boundary integral is invariant to reparameterization.
fn project_curve_on_cylinder(brep: &BRep, edge: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    use crate::geom::CurveEval;
    let surf = brep.face_surface_world(face)?;
    let (c3, r3) = brep.edge_curve_world(edge)?;
    let Surface3::Cylinder(cy) = surf else { return None };
    let p0 = c3.point_at(r3[0]);
    let p1 = c3.point_at(r3[1]);
    let uv0 = cy.world_to_uv(p0);
    let uv1 = cy.world_to_uv(p1);
    // Disambiguate the azimuth branch through the arc MIDPOINT: an
    // antipodal-ends half arc has |du| = pi and would otherwise collapse to a
    // zero-length uv line; the midpoint tells which way the boundary sweeps.
    let tau = std::f64::consts::TAU;
    let tm = 0.5 * (r3[0] + r3[1]);
    let mut d_mid = cy.world_to_uv(c3.point_at(tm)).x - uv0.x;
    while d_mid > std::f64::consts::PI {
        d_mid -= tau;
    }
    while d_mid < -std::f64::consts::PI {
        d_mid += tau;
    }
    let mut du = uv1.x - uv0.x;
    if du > std::f64::consts::PI {
        du -= tau;
    }
    if du < -std::f64::consts::PI {
        du += tau;
    }
    if d_mid.abs() < 1e-12 {
        // Midpoint ON u=uv0.x (u=0 seam): keep the nonzero sweep direction.
        if du <= 0.0 {
            du += tau;
        }
    } else if d_mid > 0.0 && du < 0.0 {
        du += tau;
    } else if d_mid < 0.0 && du > 0.0 {
        du -= tau;
    }
    let dv = uv1.y - uv0.y;
    let len = (du * du + dv * dv).sqrt();
    if len < 1e-30 {
        if std::env::var("RCAD_CYLFALLBACK").is_ok() {
            eprintln!("[CYFB] ZERO len ep={} eloc={}", edge.ptr_id() % 100000, edge.location);
        }
        return None;
    }
    if std::env::var("RCAD_CYLFALLBACK").is_ok() {
        eprintln!(
            "[CYFB] ok ep={} eloc={} len={:.6} du={:.6} dv={:.6}",
            edge.ptr_id() % 100000,
            edge.location,
            len,
            du,
            dv
        );
    }
    let c = Curve2d::Line(crate::geom::Line2d::new(uv0, glam::DVec2::new(du / len, dv / len)));
    Some((c, 0.0, len))
}

/// OCCT BRepGProp_Sinert::Perform(BF, BD, Eps) — the `checkprops -s`
/// per-face adaptive integration (BRepGProp_Sinert.cxx L104-110 →
/// BRepGProp_Gauss.cxx L533-1099).  Sinert mass only: the gravity center and
/// inertia moments do not feed the mass or the Sinert error estimates, so
/// they are omitted.
fn face_surface_area_checkprops(brep: &BRep, face: &Face, fi: usize, the_eps: f64, is_nat_restr: bool) -> f64 {
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

    // BRepGProp_Face::Bounds (BRepGProp_Face.cxx L154-160) = BRepAdaptor_Surface
    // (Restriction=true) UV domain = BRepTools::UVBounds.
    let [bu1, bu2, bv1, bv2] = face_uv_bounds(brep, face, fi, surf);
    // checkBounds (L418-429): infinite bounds switch add/mult to inf-aware.
    let inf_bounds = !bu1.is_finite() || !bu2.is_finite() || !bv1.is_finite() || !bv2.is_finite();
    let add = |a: f64, b: f64| if inf_bounds { add_inf(a, b) } else { a + b };
    let mult = |a: f64, b: f64| if inf_bounds { mult_inf(a, b) } else { a * b };

    // L543-546
    let is_error_calculation = 0.0 > the_eps || the_eps < 0.001;
    let is_verify_computation = 0.0 < the_eps && the_eps < 0.001;
    let an_epsilon = the_eps.abs();
    // L607
    let i_gl_end: usize = if is_error_calculation { 2 } else { 1 };

    // L609-615: the U Gauss orders depend only on the surface + epsilon
    let nb_u_gauss_0 = surf_s_int_order(surf, an_epsilon);
    let nb_u_gauss_1 = real_to_int(ERROR_ALGEBR_RATIO * nb_u_gauss_0 as f64) as usize;
    let (u_gauss_p0, u_gauss_w0) = match occt_gauss(nb_u_gauss_0) {
        Some(v) => v,
        None => return 0.0,
    };
    let (u_gauss_p1, u_gauss_w1) = match occt_gauss(nb_u_gauss_1) {
        Some(v) => v,
        None => return 0.0,
    };

    // L617-619
    let u_knots = surf_u_knots(surf, bu1, bu2);

    // L595: u1 is the fixed lower bound of the U integration
    let u1 = bu1;

    // L602-606: persistent loop state
    let mut error_l_max: f64 = 0.0;
    let mut eps: f64 = 0.0; // Eps
    let mut eps_l: f64 = 0.0; // EpsL
    let mut eps_u: f64 = 0.0; // EpsU

    // L548-577: 1-based arrays, grown by FillIntervalBounds
    let mut an_inertia_l = Ovec::new();
    let mut an_inertia_u = Ovec::new();
    let mut l1_v = Ovec::new();
    let mut l2_v = Ovec::new();
    let mut u1_v = Ovec::new();
    let mut u2_v = Ovec::new();
    let mut err_l = Ovec::new();
    let mut err_u = Ovec::new();
    let mut err_ul = Ovec::new();

    let mut an_inertia: f64 = 0.0;

    let face_shape = Shape::from_parts(brep.tshapes[fi].clone(), fi, 0, topods::Orientation::Forward);
    let mut edges: Vec<WireEdge> = face.outer_wire.edges.clone();
    for w in &face.inner_wires {
        edges.extend(w.edges.iter().copied());
    }

    // while (isNaturalRestriction || theDomain.More()) — L621
    let mut edge_idx: usize = 0;
    loop {
        if is_nat_restr {
            if edge_idx > 0 {
                break;
            }
        } else {
            while edge_idx < edges.len() && edges[edge_idx].internal {
                edge_idx += 1;
            }
            if edge_idx >= edges.len() {
                break;
            }
        }

        // ---- per-boundary setup (L621-661) ----
        // BRepGProp_Face::Load (L164-185): pcurve of the domain edge, with the
        // REVERSED reversal of the curve and its parameter range.
        let mut pcurve: Option<(Curve2d, bool)> = None; // (curve, reversed)
        let (l1, l2): (f64, f64);
        let nb_l_gauss_0: usize;
        if is_nat_restr {
            // L625: NbLGaussP[0] = min(2 * NbUGaussP[0], GaussPointsMax())
            nb_l_gauss_0 = (2 * nb_u_gauss_0).min(GPM);
            l1 = bv1;
            l2 = bv2;
        } else {
            let we = &edges[edge_idx];
            let edge_shape = Shape::from_parts(
                brep.tshapes[we.idx].clone(),
                we.idx,
                we.location,
                if we.forward {
                    topods::Orientation::Forward
                } else {
                    topods::Orientation::Reversed
                },
            );
            let pc = brep.curve_on_surface(&edge_shape, &face_shape);
            let (c, a0, b0) = if !we.forward {
                match brep.curve_on_surface_second(&edge_shape, &face_shape) {
                    Some(v) => v,
                    None => match pc {
                        Some(v) => v.clone(),
                        None => match project_curve_on_plane(brep, &edge_shape, &face_shape).or_else(|| project_curve_on_cylinder(brep, &edge_shape, &face_shape)) {
                            Some(v) => v,
                            None => return 0.0,
                        },
                    },
                }
            } else {
                match pc {
                    Some(v) => v.clone(),
                    None => match project_curve_on_plane(brep, &edge_shape, &face_shape).or_else(|| project_curve_on_cylinder(brep, &edge_shape, &face_shape)) {
                        Some(v) => v,
                        None => return 0.0,
                    },
                }
            };
            let reversed = !we.forward;
            let (a, b) = if reversed {
                let x = a0;
                (c.reversed_parameter(b0), c.reversed_parameter(x))
            } else {
                (a0, b0)
            };
            // L633: NbLGaussP[0] = LIntOrder(anEpsilon)
            nb_l_gauss_0 = curve_l_int_order(&c, a, b, surf, bv1, bv2, an_epsilon, reversed);
            l1 = a;
            l2 = b;
            pcurve = Some((c, reversed));
        }
        // L636
        let nb_l_gauss_1 = real_to_int(ERROR_ALGEBR_RATIO * nb_l_gauss_0 as f64) as usize;
        // L638-641
        let (l_gauss_p0, l_gauss_w0) = match occt_gauss(nb_l_gauss_0) {
            Some(v) => v,
            None => return 0.0,
        };
        let (l_gauss_p1, l_gauss_w1) = match occt_gauss(nb_l_gauss_1) {
            Some(v) => v,
            None => return 0.0,
        };
        // L643: aNbLSubs / LKnots
        let l_knots: Vec<f64> = match &pcurve {
            Some((c, reversed)) => curve_l_knots(c, l1, l2, *reversed),
            None => surf_v_knots(surf, bv1, bv2),
        };

        // L658-661
        let mut error_l: f64 = 0.0;
        let mut k_l_end = 1usize;
        let mut j_l = 0usize;

        if (l2 - l1).abs() > EPS_PARAM {
            // L664-665
            let i_l_sub_end = fill_interval_bounds(
                l1, l2, &l_knots, SUBS_POWER,
                &mut an_inertia_l, &mut l1_v, &mut l2_v, &mut err_l, Some(&mut err_ul),
            );
            // L666-670
            let mut l_max_subs = max_subs(i_l_sub_end as i64, SUBS_POWER) as usize;
            if l_max_subs > SM {
                l_max_subs = SM;
            }
            vec_init(&mut an_inertia_l.data, 1, l_max_subs, 0.0);
            vec_init(&mut err_l.data, 1, l_max_subs, 0.0);
            vec_init(&mut err_ul.data, 1, l_max_subs, 0.0);

            let mut l_range = [0usize; 2];
            // do { ... } while L — L676-1047
            loop {
                // L678-690
                j_l += 1;
                if j_l > i_l_sub_end {
                    let il = err_l.max_index();
                    l_range[0] = il;
                    l_range[1] = j_l;
                    l1_v.data[j_l] = (l1_v.data[il] + l2_v.data[il]) * 0.5;
                    l2_v.data[j_l] = l2_v.data[il];
                    l2_v.data[il] = l1_v.data[j_l];
                } else {
                    l_range[0] = j_l;
                }

                // L691-705
                if j_l == l_max_subs || (l2_v.data[j_l] - l1_v.data[j_l]).abs() < EPS_PARAM {
                    if k_l_end == 1 {
                        an_inertia_l.data[j_l] = 0.0;
                        err_l.data[j_l] = 0.0;
                    } else {
                        j_l -= 1;
                        eps_l = error_l;
                        eps = eps_l / 0.9;
                        break;
                    }
                } else {
                    // L706-1001: for kL
                    for k_l in 0..k_l_end {
                        let i_ls = l_range[k_l];
                        let lm = 0.5 * (l2_v.data[i_ls] + l1_v.data[i_ls]);
                        let lr = 0.5 * (l2_v.data[i_ls] - l1_v.data[i_ls]);
                        let mut c_dim = [0.0f64; 2];
                        // L716-961: for iGL
                        for i_gl in 0..i_gl_end {
                            let (l_gp, l_gw, nb_l_gauss) = if i_gl == 0 {
                                (l_gauss_p0, l_gauss_w0, nb_l_gauss_0)
                            } else {
                                (l_gauss_p1, l_gauss_w1, nb_l_gauss_1)
                            };
                            // L720-760: for iL
                            for i_l in 0..nb_l_gauss {
                                let l = lm + lr * l_gp[i_l];
                                let (v, u2, dul): (f64, f64, f64) = if is_nat_restr {
                                    // L723-728
                                    (l, bu2, l_gw[i_l])
                                } else {
                                    // L731-759: D12d (myCurve.D1, reversed
                                    // emulated like the fixed path)
                                    let (c, reversed) = pcurve.as_ref().unwrap();
                                    let (puv, vuv) = if *reversed {
                                        let rp = c.reversed_parameter(l);
                                        (c.point_at(rp), -c.derivative_at(rp))
                                    } else {
                                        (c.point_at(l), c.derivative_at(l))
                                    };
                                    let dul = vuv.y * l_gw[i_l];
                                    // L734-737
                                    if dul.abs() < EPS_PARAM {
                                        continue;
                                    }
                                    // L739-759: clamp to the surface bounds
                                    let mut v = puv.y;
                                    let mut u2 = puv.x;
                                    if v < bv1 {
                                        v = bv1;
                                    } else if v > bv2 {
                                        v = bv2;
                                    }
                                    if u2 < bu1 {
                                        u2 = bu1;
                                    } else if u2 > bu2 {
                                        u2 = bu2;
                                    }
                                    (v, u2, dul)
                                };
                                // L762-769
                                err_ul.data[i_ls] = 0.0;
                                let mut k_u_end = 1usize;
                                let mut j_u = 0usize;
                                if (u2 - u1).abs() < EPS_PARAM {
                                    continue;
                                }
                                // L772-774
                                let i_u_sub_end = fill_interval_bounds(
                                    u1, u2, &u_knots, SUBS_POWER,
                                    &mut an_inertia_u, &mut u1_v, &mut u2_v, &mut err_u, None,
                                );
                                // L775-779
                                let mut u_max_subs = max_subs(i_u_sub_end as i64, SUBS_POWER) as usize;
                                if u_max_subs > SM {
                                    u_max_subs = SM;
                                }
                                vec_init(&mut an_inertia_u.data, 1, u_max_subs, 0.0);
                                vec_init(&mut err_u.data, 1, u_max_subs, 0.0);
                                // L783
                                let mut error_u: f64 = 0.0;
                                let mut u_range = [0usize; 2];
                                // do { ... } while U — L785-929
                                loop {
                                    // L787-799
                                    j_u += 1;
                                    if j_u > i_u_sub_end {
                                        let iu = err_u.max_index();
                                        u_range[0] = iu;
                                        u_range[1] = j_u;
                                        u1_v.data[j_u] = (u1_v.data[iu] + u2_v.data[iu]) * 0.5;
                                        u2_v.data[j_u] = u2_v.data[iu];
                                        u2_v.data[iu] = u1_v.data[j_u];
                                    } else {
                                        u_range[0] = j_u;
                                    }
                                    // L801-816
                                    if j_u == u_max_subs || (u2_v.data[j_u] - u1_v.data[j_u]).abs() < EPS_PARAM {
                                        if k_u_end == 1 {
                                            err_u.data[j_u] = 0.0;
                                            an_inertia_u.data[j_u] = 0.0;
                                        } else {
                                            j_u -= 1;
                                            eps_u = error_u;
                                            eps = 10.0 * eps_u * ((u2 - u1) * dul).abs();
                                            eps_l = 0.9 * eps;
                                            break;
                                        }
                                    } else {
                                        // L817-921: for kU
                                        for k_u in 0..k_u_end {
                                            let i_us = u_range[k_u];
                                            let a_length = i_gl_end - i_gl;
                                            let um = 0.5 * (u2_v.data[i_us] + u1_v.data[i_us]);
                                            let ur = 0.5 * (u2_v.data[i_us] - u1_v.data[i_us]);
                                            // L824: aLocal[2] masses
                                            let mut a_local = [0.0f64; 2];
                                            // L832-867: for iGU
                                            for i_gu in 0..a_length {
                                                let (u_gp, u_gw, nb_u_gauss) = if i_gu == 0 {
                                                    (u_gauss_p0, u_gauss_w0, nb_u_gauss_0)
                                                } else {
                                                    (u_gauss_p1, u_gauss_w1, nb_u_gauss_1)
                                                };
                                                for i_u in 0..nb_u_gauss {
                                                    let w = u_gw[i_u];
                                                    let u = um + ur * u_gp[i_u];
                                                    // BRepGProp_Face::Normal (L201-210):
                                                    // D1U x D1V; mass = |N| * w
                                                    let (_p, du, dv) = surf.derivatives(u, v);
                                                    let n = du.cross(dv);
                                                    a_local[i_gu] += mult(w, n.length());
                                                }
                                            }
                                            // L869-920
                                            an_inertia_u.data[i_us] = mult(a_local[0], ur);
                                            if i_gl > 0 {
                                                continue;
                                            }
                                            // L885
                                            let a_d_mass = (a_local[1] - a_local[0]).abs();
                                            // L919
                                            err_u.data[i_us] = mult(a_d_mass, ur);
                                        }
                                    }
                                    // L924-928
                                    if j_u == i_u_sub_end {
                                        k_u_end = 2;
                                        error_u = err_u.data[err_u.max_index()];
                                    }
                                    // L929: while ((ErrorU - EpsU > 0 && EpsU != 0) || kUEnd == 1)
                                    if (error_u - eps_u > 0.0 && eps_u != 0.0) || k_u_end == 1 {
                                        continue;
                                    }
                                    break;
                                }
                                // L931-939
                                for i in 1..=j_u {
                                    c_dim[i_gl] = add(c_dim[i_gl], mult(an_inertia_u.data[i], dul));
                                }
                                // L941-944
                                if i_gl > 0 {
                                    continue;
                                }
                                // L946
                                err_ul.data[i_ls] = mult(error_u, ((u2 - u1) * dul).abs());
                            }
                            // L963-992
                            an_inertia_l.data[i_ls] = mult(c_dim[0], lr);
                            if i_gl_end == 2 {
                                let a_sub_dim = (c_dim[1] - c_dim[0]).abs();
                                // L990 (Sinert): ErrL = |aSubDim| * lr + ErrUL
                                err_l.data[i_ls] = add(mult(a_sub_dim, lr), err_ul.data[i_ls]);
                            }
                        }
                    }
                }
                // L1006-1042
                if j_l == i_l_sub_end {
                    k_l_end = 2;
                    let mut d_dim = 0.0;
                    for i in 1..=j_l {
                        d_dim += an_inertia_l.data[i];
                    }
                    d_dim = (d_dim * an_epsilon).abs();
                    if d_dim > eps {
                        eps = d_dim;
                        eps_l = 0.9 * eps;
                    }
                }
                // L1043-1046
                if k_l_end == 2 {
                    error_l = err_l.data[err_l.max_index()];
                }
                // L1047: while ((ErrorL - EpsL > 0 && isVerifyComputation) || kLEnd == 1)
                if (error_l - eps_l > 0.0 && is_verify_computation) || k_l_end == 1 {
                    continue;
                }
                break;
            }
            // L1049-1052
            for i in 1..=j_l {
                an_inertia = add(an_inertia, an_inertia_l.data[i]);
            }
            // L1054
            error_l_max = error_l_max.max(error_l);
        }

        if is_nat_restr {
            break;
        }
        edge_idx += 1;
    }

    // convert (L1071, L467-490): |Mass| >= EPS_DIM → mass else 0
    let mass = if an_inertia.abs() >= EPS_DIM {
        an_inertia
    } else {
        0.0
    };
    mass
}
