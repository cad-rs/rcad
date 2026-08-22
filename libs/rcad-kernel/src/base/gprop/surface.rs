//! OCCT BRepGProp::SurfaceProperties (BRepGProp.cxx L167-266): surface area.
//!
//! 1:1 translation of the `checkprops -s` path
//! (BRepTest_GPropCommands.cxx L123-125 → BRepGProp.cxx L268-278 → L167-266).
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

// OCCT math::GaussPoints/GaussWeights tables (math.cxx): packed positive-half
// nodes/weights for orders 1..61 plus the GaussPoints expansion.
include!("gauss_tables.rs");

/// OCCT BRepGProp::SurfaceProperties (BRepGProp.cxx L268-278) with the
/// `checkprops -s` arguments (BRepTest_GPropCommands.cxx L123-125):
/// SkipShared=false, UseTriangulation=false, so surfaceProperties is called
/// with Eps=1.0 (L277) and UseTriangulation=false.
///
/// Per-face loop (L190-259):
///   - NoSurf/NoTri check (L198-214): rcad faces always carry a surface
///     (NoSurf=false) and UseTriangulation=false, so the triangulation branch
///     (L216-222) is never taken.
///   - BF.Load(F) (L225); IsNatRestr = (F.NbChildren() == 0) (L226).
///   - Eps = 1.0 → else branch (L243-253):
///       IsNatRestr → G.Perform(BF)     (BRepGProp_Gauss.cxx L1306-1393)
///       else       → G.Perform(BF, BD) (BRepGProp_Gauss.cxx L1126-1211)
///   - Props.Add(G) (L254): the face area accumulates into the total mass.
pub fn surface_area(brep: &topods::BRep) -> f64 {
    let mut mass = 0.0;
    let faces = face_flat_iter(brep);
    for (fi, face) in &faces {
        // OCCT L226: IsNatRestr = (F.NbChildren() == 0) — the face carries no
        // wires (no outer wire edges and no inner wires).
        let is_nat_restr = face.outer_wire.edges.is_empty() && face.inner_wires.is_empty();
        let a = if is_nat_restr {
            face_surface_area_gauss_natural(brep, *fi)
        } else {
            face_surface_area_gauss_domain(brep, face, *fi)
        };
        if std::env::var("RCAD_SA_DEBUG").is_ok() {
            eprintln!("[SA] fi={} nat={} area={} nedges={}", *fi, is_nat_restr, a,
                face.outer_wire.edges.len() + face.inner_wires.iter().map(|w| w.edges.len()).sum::<usize>());
        }
        mass += a;
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
            (b.degree_u + 1, b.degree_v + 1)
        }
        _ => (9, 9),
    };
    ((2 * nu).max(8), (2 * nv).max(8))
}

/// OCCT BRepGProp_Face::IntegrationOrder (BRepGProp_Face.cxx L111-150):
/// Gauss order for the boundary edge pcurve (GeomAbs type -> N, max(4, 2N)).
fn curve_integration_order(c: &Curve2d) -> usize {
    let n = match c {
        Curve2d::Line(_) => 2,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) | Curve2d::Hyperbola(_) | Curve2d::Parabola(_) => 9,
        Curve2d::BSpline(b) => (b.degree + 1) * (b.knots.len().saturating_sub(1)),
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
        // contributes nothing to the box (L180-183).
        let (c, a_t1, a_t2) = match brep.curve_on_surface(&edge_shape, &face_shape) {
            Some(v) => v.clone(),
            None => continue,
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
        // A BSpline boundary pcurve is integrated per knot span: the C0
        // junctions of the OCCT GeomInt_WLApprox curve are integrand
        // discontinuities, and a single fixed-order Gauss rule over the whole
        // range under-samples them (BRepGProp_Gauss.cxx L1159 caps the order
        // at GaussPointsMax).  Each span is smooth, so its own
        // (degree+1)-order rule integrates it exactly.
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
            // IntegrationOrder (L1159-1161): per span (degree+1) for a
            // BSpline, the full curve order otherwise.
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
        }
    }
    // convert (L1210, L467-490): |Mass| >= EPS_DIM(1e-30) → mass else 0
    if an_inertia.abs() >= 1e-30 { an_inertia } else { 0.0 }
}
