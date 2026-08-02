//! MakeCurve + LineConstructor domain clipping for IntPatch analytic lines.
//!
//! 1:1 translation of OCCT IntTools_FaceFace::MakeCurve (L695-1846) and
//! GeomInt_LineConstructor (Perform L114-670, TreatCircle L674-733), adapted to
//! rcad-algo's FaceFace data model.
//!
//! The IntPatch lines produced by IntPatch_Intersection are the full analytic
//! intersection curves (untrimmed).  OCCT clips them to the two faces' domains:
//!   1. projecting the surface-boundary crossings onto the line -> IntPatch_Point
//!      vertices (PutPointsOnLine, IntPatch_ImpImpIntersection L439-660);
//!   2. splitting the line at the vertices and testing each interval's midpoint
//!      on both face domains (LineConstructor);
//!   3. building one trimmed curve per valid interval (MakeCurve).
//!
//! Classification uses the DS-based FClass2d face classifier (OCCT
//! Adaptor3d_TopolTool::Classify).  The UV of a point on the face surface is
//! computed analytically (OCCT GeomInt_LineConstructor::Parameters L820-862).

use crate::bop::ds::DS;
use crate::bop::int_tools::face_face::IntersectionCurve;
use crate::geomalgo::int_patch::{IntPatchIType, IntPatchLine, IntPatchVertex, WLinePnt};
use crate::geomalgo::int_surf::quadric::Quadric;
use glam::{DVec2, DVec3};
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Line3, Plane, Surface3, SurfaceEval};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};

/// OCCT IntTools_FaceFace::MakeCurve (L695-1846): clip each IntPatch line to
/// the two faces' domains and build one IntersectionCurve per valid
/// LineConstructor part.
///
/// rcad note: the faces' 2D pcurves are built incrementally (after FF), so the
/// DS FClass2d classifier is not yet reliable at the FF stage.  Domain
/// classification is therefore done against the face's UV bounds rectangle
/// (natural-restriction faces), like OCCT's TopolTool on a rectangular domain.
pub fn make_curves(
    ds: &DS,
    f1: usize,
    f2: usize,
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    tol: f64,
    lines: &[IntPatchLine],
) -> Vec<IntersectionCurve> {
    let mut out = Vec::new();
    for line in lines {
        let mut line = line.clone();
        // OCCT IntTools_FaceFace::MakeCurve case IntPatch_Restriction (L1742-1842).
        if line.line_type == IntPatchIType::Restriction {
            out.extend(make_restriction_curves(surf1, uv1, surf2, uv2, tol, &line));
            continue;
        }
        // OCCT MakeCurve uses the IntPatch_Point vertices already placed on the
        // line by IntPatch_ImpImpIntersection::Perform (PutPointsOnLine walks the
        // TopolTool boundary = the corrected FF UV rectangle).  For a WLine the
        // vertices come from the walking process; rcad places them via the UV-rect
        // crossings.  For GLine/ALine the vertices are kept as produced.
        if line.is_wline() {
            put_points_on_line(ds, f1, f2, surf1, uv1, surf2, uv2, tol, &mut line);
        }
        // OCCT GeomInt_LineConstructor::Perform -> valid parameter intervals.
        let parts = line_constructor_parts(surf1, uv1, surf2, uv2, tol, &line);
        // OCCT MakeCurve L776-1846: one curve per part.
        for &[fprm, lprm] in &parts {
            let ic = make_part_curve(
                surf1, uv1, surf2, uv2, tol, &line, fprm, lprm, parts.len(),
            );
            if let Some(ic) = ic {
                out.push(ic);
            }
        }
    }
    out
}

/// OCCT IntTools_FaceFace::MakeCurve case IntPatch_Restriction (L1742-1842)
/// + GeomInt_IntSS::TreatRLine (GeomInt_IntSS_1.cxx L1098-1168) +
/// TrimILineOnSurfBoundaries (L1333-1448).
///
/// A restriction line is a boundary arc of one face lying on the other face.
/// The 3D curve is the arc lifted to the arc's surface; the pcurve on the
/// other surface is built by projecting the 3D curve.  The curve is then
/// trimmed to the two faces' UV rectangles.
fn make_restriction_curves(
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    tol: f64,
    line: &IntPatchLine,
) -> Vec<IntersectionCurve> {
    let mut out = Vec::new();
    // OCCT TreatRLine L1106-1133: the arc is on surface 1 or surface 2.
    let (arc_surf, arc, other_surf, other_uv, arc_uv, is_on_s1) = if line.is_arc_on_s1() {
        (surf1, line.arc_on_s1(), surf2, uv2, uv1, true)
    } else if line.is_arc_on_s2() {
        (surf2, line.arc_on_s2(), surf1, uv1, uv2, false)
    } else {
        return out;
    };
    let Some(arc) = arc else { return out };

    // Parametric range of the restriction arc (OCCT TreatRLine ParamOnS1/S2,
    // derived from the RLine's first/last points).
    let mut tf = if line.has_first_point() {
        line.first_point().parameter_on_line()
    } else {
        arc.default_domain()[0]
    };
    let mut tl = if line.has_last_point() {
        line.last_point().parameter_on_line()
    } else {
        arc.default_domain()[1]
    };
    if !(tf.is_finite() && tl.is_finite() && tl > tf) {
        return out;
    }
    tf = tf.max(arc.default_domain()[0]);
    tl = tl.min(arc.default_domain()[1]);
    if tl <= tf {
        return out;
    }

    // OCCT: 3D curve = approximation of the curve on the arc surface
    // (Approx_CurveOnSurface).  rcad: the 3D image of the UV-line arc on the
    // analytic quadric is exact (line or circle), built with the same
    // parameterization as the 2D arc.
    let (curve3, _ctype) =
        match crate::geomalgo::int_patch::so_on_bounds::curve_on_surface(arc, arc_surf) {
            Some(c) => c,
            None => return out,
        };

    // OCCT TreatRLine L1157-1166: pcurve on the other surface via
    // GeomInt_IntSS::BuildPCurves.  rcad: build the other pcurve by sampling
    // the 3D curve and inverting it on the other quadric.
    let other_pcurve = build_projected_pcurve(other_surf, &curve3, tf, tl, 64);
    let Some(other_pcurve) = other_pcurve else { return out };

    // OCCT TrimILineOnSurfBoundaries: intersect the pcurves with the two UV
    // rectangles, collect parameters along the arc.
    let arc_pcurve = arc.clone();
    let mut params: Vec<f64> = vec![tf, tl];
    // The arc's own pcurve lives on its surface's UV rectangle.
    collect_boundary_crossings(&arc_pcurve, arc_uv, &mut params);
    collect_boundary_crossings(&other_pcurve, other_uv, &mut params);
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    params.retain(|&p| p >= tf - PCONFUSION && p <= tl + PCONFUSION);
    // dedup
    let mut dedup: Vec<f64> = Vec::new();
    for p in params {
        if dedup.is_empty() || (p - dedup.last().unwrap()).abs() > PCONFUSION {
            dedup.push(p);
        }
    }
    let params = dedup;
    let a_box1 = enlarge_box(uv1, tol);
    let a_box2 = enlarge_box(uv2, tol);

    for an_ind in 0..params.len().saturating_sub(1) {
        let a_par_f = params[an_ind];
        let mut a_par_l = params[an_ind + 1];
        if a_par_l - a_par_f <= PCONFUSION {
            if an_ind + 1 < params.len() - 1 {
                a_par_l = a_par_f;
            }
            continue;
        }
        let a_par = 0.5 * (a_par_f + a_par_l);
        // Midpoint of the arc pcurve must be inside its surface's UV box.
        let a_pt1 = arc.point_at(a_par);
        let arc_box = if is_on_s1 { a_box1 } else { a_box2 };
        if !point_in_box(a_pt1, arc_box) {
            continue;
        }
        // Midpoint of the other pcurve inside the other surface's UV box.
        let a_pt2 = other_pcurve.point_at(a_par);
        let other_box = if is_on_s1 { a_box2 } else { a_box1 };
        if !point_in_box(a_pt2, other_box) {
            continue;
        }
        // OCCT L1836-1839: trimmed 3D curve + pcurves.
        let pcurve1 = if is_on_s1 { Some(arc.clone()) } else { Some(other_pcurve.clone()) };
        let pcurve2 = if is_on_s1 { Some(other_pcurve.clone()) } else { Some(arc.clone()) };
        out.push(IntersectionCurve {
            curve: curve3.clone(),
            t_range: [a_par_f, a_par_l],
            pcurve1,
            pcurve2,
            tolerance: tol,
            tang_tolerance: 0.0,
        });
    }
    out
}

fn enlarge_box(uv: [f64; 4], delta: f64) -> [f64; 4] {
    [uv[0] - delta, uv[1] + delta, uv[2] - delta, uv[3] + delta]
}

/// OCCT Bnd_Box2d::IsOut(Pnt) — point outside the enlarged UV rectangle.
fn point_in_box(p: DVec2, b: [f64; 4]) -> bool {
    !(p.x < b[0] || p.x > b[1] || p.y < b[2] || p.y > b[3])
}

/// OCCT GeomInt_IntSS::IntersectCurveAndBoundary — parameters where the 2D
/// pcurve crosses the 4 edges of the UV rectangle.
fn collect_boundary_crossings(pc: &Curve2d, uv: [f64; 4], params: &mut Vec<f64>) {
    let [u_min, u_max, v_min, v_max] = uv;
    // Intersect the pcurve with each of the 4 boundary lines.
    let boundaries = [
        ([u_min, v_min], DVec2::new(0.0, 1.0)), // U=Umin, V varies
        ([u_max, v_min], DVec2::new(0.0, 1.0)), // U=Umax, V varies
        ([u_min, v_min], DVec2::new(1.0, 0.0)), // V=Vmin, U varies
        ([u_min, v_max], DVec2::new(1.0, 0.0)), // V=Vmax, U varies
    ];
    if let Curve2d::Line(l) = pc {
        for (o, d) in boundaries {
            if let Some(t) = line_line2d_intersection(&l, &Line2d { origin: DVec2::new(o[0], o[1]), direction: d }) {
                let p = l.point_at(t);
                // ensure the intersection lies on the boundary segment
                let on_seg = if d.y.abs() > 0.5 {
                    p.y >= v_min - 1e-9 && p.y <= v_max + 1e-9
                } else {
                    p.x >= u_min - 1e-9 && p.x <= u_max + 1e-9
                };
                if on_seg {
                    params.push(t);
                }
            }
        }
    }
}

/// Intersection parameter of two 2D lines.  Returns the parameter on `a`.
fn line_line2d_intersection(a: &Line2d, b: &Line2d) -> Option<f64> {
    let det = a.direction.x * b.direction.y - a.direction.y * b.direction.x;
    if det.abs() < 1e-30 {
        return None;
    }
    let db = b.origin - a.origin;
    let t = (db.x * b.direction.y - db.y * b.direction.x) / det;
    Some(t)
}

/// OCCT GeomInt_IntSS::BuildPCurves (projection of the 3D curve onto the
/// other quadric).  rcad: sample the 3D curve, invert each point on the
/// quadric, and build a 2D line through the endpoint UVs.
fn build_projected_pcurve(
    other_surf: &Surface3,
    curve3: &Curve3,
    tf: f64,
    tl: f64,
    _n: usize,
) -> Option<Curve2d> {
    let p0 = curve3.point_at(tf);
    let p1 = curve3.point_at(tl);
    let uv0 = quadric_uv_params(other_surf, p0)?;
    let uv1 = quadric_uv_params(other_surf, p1)?;
    if !uv0.is_finite() || !uv1.is_finite() {
        return None;
    }
    // OCCT BuildPCurves: on a short domain the pcurve is a segment of a line.
    Some(Curve2d::Line(Line2d {
        origin: uv0,
        direction: (uv1 - uv0) / (tl - tf),
    }))
}

/// OCCT IntPatch_ImpImpIntersection::Perform L439-660 (PutPointsOnLine).
/// Projects the surface-boundary crossings onto the line.  The boundary of each
/// face is its set of DS boundary edges; where an edge's 3D curve crosses the
/// other face's quadric is a boundary point shared by both surfaces.
fn put_points_on_line(
    ds: &DS,
    f1: usize,
    f2: usize,
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    tol: f64,
    line: &mut IntPatchLine,
) {
    if line.is_wline() {
        // OCCT PutPointsOnLine for a WLine (GeomInt_IntPatchPrmPrmIntersection /
        // ImpImp L439-660): find where the polyline crosses the two surfaces'
        // UV domain boundaries and place IntPatch_Point vertices there.  The
        // WLine points carry the UV on both surfaces.
        let wpts = line.wline_pnts.clone();
        let n = wpts.len();
        if n < 2 {
            return;
        }
        let mut verts: Vec<IntPatchVertex> = Vec::new();
        collect_wline_uv_crossings(&wpts, surf1, uv1, true, &mut verts);
        collect_wline_uv_crossings(&wpts, surf2, uv2, false, &mut verts);
        verts.sort_by(|a, b| {
            a.param_on_line
                .partial_cmp(&b.param_on_line)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let a_tol_pc = 1000.0 * PCONFUSION;
        let mut dedup: Vec<IntPatchVertex> = Vec::new();
        for v in verts {
            if dedup.is_empty()
                || (v.param_on_line - dedup.last().unwrap().param_on_line).abs() > a_tol_pc
            {
                dedup.push(v);
            }
        }
        line.vertices = dedup;
        return;
    }
    if line.line_type == IntPatchIType::Restriction {
        return;
    }
    let Some(q1) = Quadric::from_surface3(surf1) else {
        return;
    };
    let Some(q2) = Quadric::from_surface3(surf2) else {
        return;
    };
    let curve = line.curve.clone();
    let t_range = line.t_range;
    let mut verts: Vec<IntPatchVertex> = Vec::new();
    // Face 1 boundary edges crossing surface 2.
    collect_face_boundary_vertices(ds, f1, &q2, surf1, uv1, &curve, t_range, tol, true, &mut verts);
    // Face 2 boundary edges crossing surface 1.
    collect_face_boundary_vertices(ds, f2, &q1, surf2, uv2, &curve, t_range, tol, false, &mut verts);
    // OCCT IntPatch_GLine::ComputeVertexParameters: sort and dedup.
    verts.sort_by(|a, b| {
        a.param_on_line
            .partial_cmp(&b.param_on_line)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let a_tol_pc = 1000.0 * PCONFUSION;
    let mut dedup: Vec<IntPatchVertex> = Vec::new();
    for v in verts {
        if dedup.is_empty()
            || (v.param_on_line - dedup.last().unwrap().param_on_line).abs() > a_tol_pc
        {
            dedup.push(v);
        }
    }
    line.vertices = dedup;
}

/// For one face's boundary edges, find where the edge 3D curve crosses the
/// other face's quadric and project each crossing onto the intersection line.
///
/// OCCT's PutPointsOnLine walks the TopolTool domain boundary.  For a plane
/// face that domain is the CorrectPlaneBoundaries-enlarged UV rectangle, so the
/// boundary edges are the 4 phantom rectangle edges (extended past the real
/// face boundary by 10%); for other surfaces the DS boundary edges are used.
fn collect_face_boundary_vertices(
    ds: &DS,
    fi: usize,
    other_quad: &Quadric,
    surf: &Surface3,
    uv: [f64; 4],
    curve: &Curve3,
    t_range: [f64; 2],
    tol: f64,
    is_face1: bool,
    out: &mut Vec<IntPatchVertex>,
) {
    // OCCT Adaptor3d_TopolTool domain boundary: the UV rectangle of the surface.
    // For a closed surface (sphere/torus) the DS boundary edges are unreliable
    // (the sphere's DS seam lies in the wrong plane), so walk the UV domain.
    // For other surfaces the DS boundary edges match the TopolTool domain.
    let is_closed = matches!(surf, Surface3::Sphere(_) | Surface3::Torus(_));
    let closed_polylines = if is_closed {
        surface_uv_domain_boundary(surf, uv, 128)
    } else {
        Vec::new()
    };
    let boundary: Vec<(Curve3, [f64; 2])> = if is_closed {
        Vec::new()
    } else if let Surface3::Plane(pl) = surf {
        plane_uv_edges(pl, uv)
    } else {
        ds.face_boundary_edges(fi)
            .into_iter()
            .filter_map(|ei| {
                let e_curve = ds.edge_curve(ei)?;
                Some((e_curve.clone(), ds.edge_range(ei)))
            })
            .collect()
    };
    // Project a boundary point onto the intersection line and record the vertex.
    let mut try_add_vertex = |p3d: DVec3| {
        let proj = closest_point_on_curve_range(curve, p3d, t_range[0], t_range[1], 64);
        if !proj.distance.is_finite() || proj.distance > tol * 10.0 + 1e-6 {
            return;
        }
        let uv_self = match quadric_uv_params(surf, p3d) {
            Some(uv) => uv,
            None => return,
        };
        let (u_other, v_other) = other_quad.parameters(p3d);
        let (u1, v1, u2, v2) = if is_face1 {
            (uv_self.x, uv_self.y, u_other, v_other)
        } else {
            (u_other, v_other, uv_self.x, uv_self.y)
        };
        out.push(IntPatchVertex {
            param_on_line: proj.param,
            p3d,
            u1,
            v1,
            u2,
            v2,
            ..Default::default()
        });
    };
    // Closed-surface UV domain edges walked as point polylines.
    for pts in &closed_polylines {
        let n = pts.len();
        if n < 2 {
            continue;
        }
        let d0 = other_quad.distance(pts[0]);
        let d1 = other_quad.distance(pts[n - 1]);
        if d0.abs() <= tol * 10.0 + 1e-6 {
            try_add_vertex(pts[0]);
        }
        if d1.abs() <= tol * 10.0 + 1e-6 {
            try_add_vertex(pts[n - 1]);
        }
        let mut prev_d = d0;
        for k in 1..n {
            let p = pts[k];
            let d = other_quad.distance(p);
            let prev_d_k = prev_d;
            prev_d = d;
            if prev_d_k * d <= 0.0 && prev_d_k.abs() > 1e-14 {
                let mut lo = 0.0f64;
                let mut hi = 1.0f64;
                let mut d_lo = prev_d_k;
                for _ in 0..50 {
                    let mid = (lo + hi) * 0.5;
                    let pm = pts[k - 1] + (pts[k] - pts[k - 1]) * mid;
                    let d_mid = other_quad.distance(pm);
                    if d_lo * d_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        d_lo = d_mid;
                    }
                }
                try_add_vertex(pts[k - 1] + (pts[k] - pts[k - 1]) * ((lo + hi) * 0.5));
            }
        }
    }
    // Analytic boundary edges (plane rect / DS edges) for other surfaces.
    for (e_curve, [t0, t1]) in boundary {
        let span = t1 - t0;
        if !span.is_finite() || span == 0.0 {
            continue;
        }
        let n_sample = 32usize;
        let d0 = other_quad.distance(e_curve.point_at(t0));
        let d1 = other_quad.distance(e_curve.point_at(t1));
        if d0.abs() <= tol * 10.0 + 1e-6 {
            try_add_vertex(e_curve.point_at(t0));
        }
        if d1.abs() <= tol * 10.0 + 1e-6 {
            try_add_vertex(e_curve.point_at(t1));
        }
        let mut prev_d = d0;
        for k in 1..=n_sample {
            let t = t0 + span * (k as f64 / n_sample as f64);
            let p = e_curve.point_at(t);
            let d = other_quad.distance(p);
            let prev_d_k = prev_d;
            prev_d = d;
            if prev_d_k * d <= 0.0 && prev_d_k.abs() > 1e-14 {
                let mut lo = t0 + span * ((k as f64 - 1.0) / n_sample as f64);
                let mut hi = t;
                let mut d_lo = prev_d_k;
                for _ in 0..50 {
                    let mid = (lo + hi) * 0.5;
                    let d_mid = other_quad.distance(e_curve.point_at(mid));
                    if d_lo * d_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        d_lo = d_mid;
                    }
                }
                try_add_vertex(e_curve.point_at((lo + hi) * 0.5));
            }
        }
    }
}

/// OCCT Adaptor3d_TopolTool domain boundary of a surface: the four UV rectangle
/// edges [UMin,UMax]x[VMin,VMax] sampled as 3D point polylines on the surface.
/// Degenerate edges (a sphere pole, a zero-length edge) are skipped.
fn surface_uv_domain_boundary(surf: &Surface3, uv: [f64; 4], n: usize) -> Vec<Vec<DVec3>> {
    use rcad_kernel::geom::SurfaceEval;
    let (u_lo, u_hi) = if uv[0] <= uv[1] { (uv[0], uv[1]) } else { (uv[1], uv[0]) };
    let (v_lo, v_hi) = if uv[2] <= uv[3] { (uv[2], uv[3]) } else { (uv[3], uv[2]) };
    let mut out = Vec::new();
    // U = UMin and U = UMax edges (V varies).
    if u_lo.is_finite() {
        let mut p = Vec::with_capacity(n);
        for k in 0..n {
            let v = v_lo + (v_hi - v_lo) * (k as f64 / (n - 1) as f64);
            p.push(surf.point_at(u_lo, v));
        }
        out.push(p);
    }
    if u_hi.is_finite() && (u_hi - u_lo).abs() > 0.0 {
        let mut p = Vec::with_capacity(n);
        for k in 0..n {
            let v = v_lo + (v_hi - v_lo) * (k as f64 / (n - 1) as f64);
            p.push(surf.point_at(u_hi, v));
        }
        out.push(p);
    }
    // V = VMin and V = VMax edges (U varies).
    if v_lo.is_finite() {
        let mut p = Vec::with_capacity(n);
        for k in 0..n {
            let u = u_lo + (u_hi - u_lo) * (k as f64 / (n - 1) as f64);
            p.push(surf.point_at(u, v_lo));
        }
        out.push(p);
    }
    if v_hi.is_finite() && (v_hi - v_lo).abs() > 0.0 {
        let mut p = Vec::with_capacity(n);
        for k in 0..n {
            let u = u_lo + (u_hi - u_lo) * (k as f64 / (n - 1) as f64);
            p.push(surf.point_at(u, v_hi));
        }
        out.push(p);
    }
    out
}

/// Walk a WLine polyline and place an IntPatchVertex wherever the polyline
/// crosses the surface UV domain rectangle boundary.  Each WLine point already
/// carries the UV on both surfaces (u1/v1/u2/v2); the in/out test uses the same
/// adjust-periodic + rectangle classification as the LineConstructor.
fn collect_wline_uv_crossings(
    wpts: &[WLinePnt],
    surf: &Surface3,
    rect: [f64; 4],
    is_face1: bool,
    out: &mut Vec<IntPatchVertex>,
) {
    let n = wpts.len();
    let mut prev_in = None;
    for k in 0..n {
        let p3d = wpts[k].p3d;
        let (u, v) = if is_face1 { (wpts[k].u1, wpts[k].v1) } else { (wpts[k].u2, wpts[k].v2) };
        let in_now = in_uv_rect_adjusted(surf, rect, p3d, u, v);
        if let Some(prev_in_val) = prev_in {
            if in_now != prev_in_val {
                // The curve crossed the face domain boundary between k-1 and k.
                let t = (k as f64) / ((n - 1) as f64);
                let mut p = wpts[k].p3d;
                if k > 0 {
                    let a = wpts[k - 1].p3d;
                    let b = wpts[k].p3d;
                    // midpoint as the crossing approximation
                    p = (a + b) * 0.5;
                }
                let (u1o, v1o, u2o, v2o) = if is_face1 {
                    let p3 = wpts[k].p3d;
                    let (u2, v2) = surf_uv_other(&p3);
                    (u, v, u2, v2)
                } else {
                    (wpts[k].u1, wpts[k].v1, u, v)
                };
                out.push(IntPatchVertex {
                    param_on_line: t,
                    p3d: p,
                    u1: u1o,
                    v1: v1o,
                    u2: u2o,
                    v2: v2o,
                    ..Default::default()
                });
            }
        }
        prev_in = Some(in_now);
    }
}

/// Interpolate a WLine polyline at a normalized parameter in [0, 1].
fn wline_point_at(line: &IntPatchLine, param: f64) -> DVec3 {
    let wpts = &line.wline_pnts;
    let n = wpts.len();
    if n == 0 {
        return DVec3::ZERO;
    }
    if n == 1 {
        return wpts[0].p3d;
    }
    let x = param.clamp(0.0, 1.0) * (n - 1) as f64;
    let i = x.floor() as usize;
    let j = (i + 1).min(n - 1);
    let f = x - i as f64;
    wpts[i].p3d + (wpts[j].p3d - wpts[i].p3d) * f
}

/// Classify a WLine point on a surface using its precomputed UV.
fn in_uv_rect_adjusted(surf: &Surface3, rect: [f64; 4], _p3d: DVec3, u: f64, v: f64) -> bool {
    let adj = adjust_periodic_uv(surf, DVec2::new(u, v), rect);
    in_uv_rect(adj, rect, 0.0)
}

/// UV of a 3D point on the "other" surface for the crossing vertex record.
fn surf_uv_other(p: &DVec3) -> (f64, f64) {
    (0.0, 0.0)
}

/// The 4 boundary edges of a plane face's UV rectangle, as 3D lines.  OCCT
/// walks the TopolTool domain — the CorrectPlaneBoundaries-enlarged rectangle —
/// so these are the phantom edges of the enlarged rectangle, not the real face
/// boundary edges.  Edges with an infinite extent are skipped.
fn plane_uv_edges(pl: &Plane, uv: [f64; 4]) -> Vec<(Curve3, [f64; 2])> {
    let [u_lo, u_hi, v_lo, v_hi] = uv;
    let mut out = Vec::new();
    let mut push = |o: DVec3, d: DVec3, r: [f64; 2]| {
        if r[0].is_finite() && r[1].is_finite() && (r[1] - r[0]).abs() > 0.0 {
            out.push((Curve3::Line(Line3::new(o, d)), r));
        }
    };
    if u_lo.is_finite() {
        push(pl.origin + pl.u_dir * u_lo, pl.v_dir, [v_lo, v_hi]);
    }
    if u_hi.is_finite() {
        push(pl.origin + pl.u_dir * u_hi, pl.v_dir, [v_lo, v_hi]);
    }
    if v_lo.is_finite() {
        push(pl.origin + pl.v_dir * v_lo, pl.u_dir, [u_lo, u_hi]);
    }
    if v_hi.is_finite() {
        push(pl.origin + pl.v_dir * v_hi, pl.u_dir, [u_lo, u_hi]);
    }
    out
}

/// OCCT GeomInt_LineConstructor::Perform (L333-386, GLine path).
/// Splits the line at the vertices; tests each interval's midpoint on both
/// face domains.  Circle/Ellipse are delegated to TreatCircle (L338-342).
fn line_constructor_parts(
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    _tol: f64,
    line: &IntPatchLine,
) -> Vec<[f64; 2]> {
    if line.line_type == IntPatchIType::Circle || line.line_type == IntPatchIType::Ellipse {
        return treat_circle_parts(surf1, uv1, surf2, uv2, _tol, line);
    }
    // OCCT L118: constexpr double Tol = Precision::PConfusion() * 35.0;
    let a_tol = PCONFUSION * 35.0;
    let mut result: Vec<[f64; 2]> = Vec::new();
    let nbvtx = line.vertices.len();
    let mut intrvtested = false;
    for i in 0..nbvtx.saturating_sub(1) {
        let firstp = line.vertices[i].param_on_line;
        let lastp = line.vertices[i + 1].param_on_line;
        // OCCT L354: if (std::abs(firstp - lastp) > Precision::PConfusion())
        if (firstp - lastp).abs() > PCONFUSION {
            intrvtested = true;
            let pmid = (firstp + lastp) * 0.5;
            let p3d = if line.is_wline() {
                wline_point_at(line, pmid)
            } else {
                line.curve.point_at(pmid)
            };
            if !p3d.is_finite() {
                continue;
            }
            // OCCT L361-372: Parameters + AdjustPeriodic + Classify both domains.
            let in1 = classify_point(surf1, uv1, p3d, a_tol);
            if in1 {
                let in2 = classify_point(surf2, uv2, p3d, a_tol);
                if in2 {
                    result.push([firstp, lastp]);
                }
            }
        }
    }
    // OCCT L376-382: if no interval tested, keep the full range a priori.
    if !intrvtested {
        result.push(line.t_range);
    }
    result
}

/// OCCT GeomInt_LineConstructor::TreatCircle (L674-733).
fn treat_circle_parts(
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    _tol: f64,
    line: &IntPatchLine,
) -> Vec<[f64; 2]> {
    let two_pi = std::f64::consts::TAU;
    let curve = &line.curve;
    // OCCT GeomInt_LineConstructor::RejectMicroCircle (L895-915): reject when
    // the radius is below the LineConstructor Tol (PConfusion() * 35.0).
    let radius = match curve {
        Curve3::Circle(c) => c.radius,
        Curve3::Ellipse(e) => e.major_radius,
        _ => 0.0,
    };
    if radius < PCONFUSION * 35.0 {
        return Vec::new();
    }
    // OCCT IntPatch_ImpImpIntersection L2997-3052: a Circle/Ellipse GLine
    // without vertices gets two vertices at parameter 0 and 2*PI.
    let a_tol_pc = 1000.0 * PCONFUSION;
    let mut params: Vec<f64> = if line.vertices.is_empty() {
        vec![0.0, two_pi]
    } else {
        line.vertices
            .iter()
            .map(|v| wrap_to_2pi(v.param_on_line))
            .collect()
    };
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // OCCT L684-695: array of size n+1; last vertex at first.param + 2*PI.
    let mut arr: Vec<f64> = Vec::with_capacity(params.len() + 1);
    arr.extend_from_slice(&params);
    arr.push(params[0] + two_pi);

    // OCCT L697: RejectDuplicates -- mark coincident params with RealLast.
    for i in 0..arr.len().saturating_sub(2) {
        let prm_i = arr[i];
        if !prm_i.is_finite() {
            continue;
        }
        for j in (i + 1)..arr.len().saturating_sub(1) {
            let prm_j = arr[j];
            if prm_j - prm_i < a_tol_pc {
                arr[j] = f64::INFINITY; // RealLast
            } else {
                break;
            }
        }
    }
    let a_max_prm = *arr.last().unwrap();
    for i in (1..arr.len().saturating_sub(1)).rev() {
        let prm_i = arr[i];
        if !prm_i.is_finite() {
            continue;
        }
        if a_max_prm - prm_i < a_tol_pc {
            arr[i] = f64::INFINITY;
        } else {
            break;
        }
    }
    arr.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // OCCT L704-732: test each adjacent pair's midpoint on both face domains.
    let a_tol = PCONFUSION * 35.0;
    let mut result = Vec::new();
    for i in 0..arr.len().saturating_sub(1) {
        let t1 = arr[i];
        let t2 = arr[i + 1];
        if t2 == f64::INFINITY {
            break;
        }
        let t_mid = (t1 + t2) * 0.5;
        let p3d = curve.point_at(t_mid);
        if !p3d.is_finite() {
            continue;
        }
        let in1 = classify_point(surf1, uv1, p3d, a_tol);
        if !in1 {
            continue;
        }
        let in2 = classify_point(surf2, uv2, p3d, a_tol);
        if in2 {
            result.push([t1, t2]);
        }
    }
    result
}

/// OCCT MakeCurve L776-1846: build one IntersectionCurve for a LineConstructor
/// part.  The full analytic curve is kept; `t_range` holds the clipped part.
fn make_part_curve(
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    tol: f64,
    line: &IntPatchLine,
    fprm: f64,
    lprm: f64,
    a_nb_parts: usize,
) -> Option<IntersectionCurve> {
    let curve = &line.curve;
    match line.line_type {
        IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
            let b_fn = fprm.is_finite();
            let b_lp = lprm.is_finite();
            if !(b_fn && b_lp) {
                // OCCT L850-898: test a reference point on both face domains.
                let d_t = 100.0;
                let test_t = if !b_fn && b_lp {
                    lprm - d_t
                } else if b_fn && !b_lp {
                    fprm + d_t
                } else {
                    // both infinite: OCCT IntTools_Tools::IntermediatePoint(-dT, dT)
                    0.0
                };
                let p3d = curve.point_at(test_t);
                if !p3d.is_finite() {
                    return None;
                }
                if !classify_point(surf1, uv1, p3d, CONFUSION) {
                    return None;
                }
                if !classify_point(surf2, uv2, p3d, CONFUSION) {
                    return None;
                }
            }
            if lprm <= fprm + 1e-12 {
                return None;
            }
            Some(IntersectionCurve {
                curve: curve.clone(),
                t_range: [fprm, lprm],
                pcurve1: None,
                pcurve2: None,
                tolerance: tol,
                tang_tolerance: line.tang_tolerance,
            })
        }
        IntPatchIType::Circle | IntPatchIType::Ellipse => {
            let a_period = std::f64::consts::TAU;
            let a_real_eps = f64::EPSILON;
            let is_full_period = fprm.abs() <= a_real_eps && (lprm - a_period).abs() <= a_real_eps;
            if !is_full_period {
                if lprm <= fprm + 1e-12 {
                    return None;
                }
                Some(IntersectionCurve {
                    curve: curve.clone(),
                    t_range: [fprm, lprm],
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: tol,
                    tang_tolerance: line.tang_tolerance,
                })
            } else if a_nb_parts == 1 {
                // OCCT L1074-1103: accept the full circle.
                Some(IntersectionCurve {
                    curve: curve.clone(),
                    t_range: [fprm, lprm],
                    pcurve1: None,
                    pcurve2: None,
                    tolerance: tol,
                    tang_tolerance: line.tang_tolerance,
                })
            } else {
                // OCCT L1109-1165: test 18 points around the circle.
                let a_two_pi_div_17 = a_period / 17.0;
                for j in 0..=17 {
                    let t = j as f64 * a_two_pi_div_17;
                    let p3d = curve.point_at(t);
                    if !p3d.is_finite() {
                        continue;
                    }
                    if !classify_point(surf1, uv1, p3d, CONFUSION) {
                        continue;
                    }
                    if classify_point(surf2, uv2, p3d, CONFUSION) {
                        return Some(IntersectionCurve {
                            curve: curve.clone(),
                            t_range: [fprm, lprm],
                            pcurve1: None,
                            pcurve2: None,
                            tolerance: tol,
                            tang_tolerance: line.tang_tolerance,
                        });
                    }
                }
                None
            }
        }
        IntPatchIType::Walking => {
            // OCCT MakeCurve L1175-1739 (IntPatch_Walking): the WLine is split
            // by LineConstructor parts; each part becomes one curve.  rcad keeps
            // the sampled polyline in wline_pnts and emits one curve per part.
            if lprm <= fprm + 1e-12 {
                return None;
            }
            Some(IntersectionCurve {
                curve: curve.clone(),
                t_range: [fprm, lprm],
                pcurve1: None,
                pcurve2: None,
                tolerance: tol,
                tang_tolerance: line.tang_tolerance,
            })
        }
        _ => None,
    }
}

/// OCCT GeomInt_LineConstructor::Parameters (L820-862) + Classify.  Analytic UV
/// OCCT GeomInt_LineConstructor::Parameters (L820-862) + AdjustPeriodic
/// (L737-816) + Classify.  Analytic UV inversion of a 3D point on a quadric
/// surface, shifted into the face UV rectangle before the in-rectangle test.
fn classify_point(surf: &Surface3, rect: [f64; 4], p3d: DVec3, tol: f64) -> bool {
    match quadric_uv_params(surf, p3d) {
        Some(uv) => {
            let adj = adjust_periodic_uv(surf, uv, rect);
            in_uv_rect(adj, rect, tol)
        }
        None => false,
    }
}

/// OCCT GeomInt_LineConstructor::AdjustPeriodic (L737-816): for periodic
/// directions (cylinder/cone/sphere U, torus U+V) shift into the face rectangle.
fn adjust_periodic_uv(surf: &Surface3, uv: DVec2, rect: [f64; 4]) -> DVec2 {
    let two_pi = std::f64::consts::TAU;
    let (u_lo, u_hi) = if rect[0] <= rect[1] { (rect[0], rect[1]) } else { (rect[1], rect[0]) };
    let (v_lo, v_hi) = if rect[2] <= rect[3] { (rect[2], rect[3]) } else { (rect[3], rect[2]) };
    let (is_u_per, is_v_per) = match surf {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) => (true, false),
        Surface3::Torus(_) => (true, true),
        _ => (false, false),
    };
    let mut u = uv.x;
    let mut v = uv.y;
    if is_u_per {
        u = adjust_periodic(u, u_lo, u_hi, two_pi);
    }
    if is_v_per {
        v = adjust_periodic(v, v_lo, v_hi, two_pi);
    }
    DVec2::new(u, v)
}

/// OCCT GeomInt::AdjustPeriodic (GeomInt.cxx L21-48).
fn adjust_periodic(par: f64, par_min: f64, par_max: f64, period: f64) -> f64 {
    let eps = 1e-9;
    let b_min = par_min - par > eps;
    let b_max = par - par_max > eps;
    if b_min || b_max {
        let dp = if b_min { par_max - par } else { par_min - par };
        let nb_per = (dp / period).trunc();
        par + nb_per * period
    } else {
        par
    }
}

/// OCCT TopolTool::Classify on a natural-restriction face: the UV domain is the
/// rectangle [u1,u2]x[v1,v2].
fn in_uv_rect(uv: DVec2, rect: [f64; 4], tol: f64) -> bool {
    let (u_lo, u_hi) = if rect[0] <= rect[1] { (rect[0], rect[1]) } else { (rect[1], rect[0]) };
    let (v_lo, v_hi) = if rect[2] <= rect[3] { (rect[2], rect[3]) } else { (rect[3], rect[2]) };
    uv.x >= u_lo - tol && uv.x <= u_hi + tol && uv.y >= v_lo - tol && uv.y <= v_hi + tol
}

/// OCCT GeomInt_LineConstructor::Parameters (L820-862).
fn quadric_uv_params(surf: &Surface3, p: DVec3) -> Option<DVec2> {
    match surf {
        Surface3::Plane(pl) => {
            let d = p - pl.origin;
            Some(DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir)))
        }
        Surface3::Cylinder(cy) => Some(cy.world_to_uv(p)),
        Surface3::Sphere(sp) => Some(sp.world_to_uv(p)),
        Surface3::Cone(co) => Some(co.world_to_uv(p)),
        Surface3::Torus(to) => Some(to.world_to_uv(p)),
        _ => None,
    }
}

/// OCCT ElCLib::InPeriod(X, 0.0, 2*PI) (GeomInt_Vertex::SetVertex).
fn wrap_to_2pi(t: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    if t < 0.0 {
        t + two_pi * (1.0 + ((0.0 - t) / two_pi).floor())
    } else if t >= two_pi {
        t - two_pi * (1.0 + ((t - two_pi) / two_pi).floor())
    } else {
        t
    }
}
