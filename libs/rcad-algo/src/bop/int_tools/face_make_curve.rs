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
use crate::geomalgo::int_patch::{IntPatchIType, IntPatchLine, IntPatchVertex, WLinePnt, WLineType};
use crate::geomalgo::int_surf::quadric::Quadric;
use glam::{DVec2, DVec3};
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::geom::{
    BSplineCurve2, BSplineCurve3, Circle2d, Curve2d, Curve2dEval, Curve3, CurveEval, Line2d,
    Line3, Parabola3, Plane, Surface3, SurfaceEval,
};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topods::{ShapeType, Shape, TShape};

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
    approx: bool,
    approx1: bool,
    approx2: bool,
    tol_approx: f64,
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
        // OCCT IntTools_FaceFace::MakeCurve: the IntPatch line vertices (placed
        // by IntPatch_ImpImpIntersection::Perform / the walking process) are used
        // as-is; no additional UV-crossing placement is done for a WLine.
        // OCCT GeomInt_LineConstructor::Perform -> valid parameter intervals.
        let parts = line_constructor_parts(surf1, uv1, surf2, uv2, tol, &line);
        // OCCT MakeCurve L926-1000: a Circle/Ellipse interval crossing 0 is
        // divided on two intervals [fprm, 2*PI] and [0, lprm].
        let a_nul = 0.0;
        let a_period = std::f64::consts::TAU;
        let a_tol_pc = rcad_kernel::precision::APPROXIMATION; // myTolApprox
        let mut a_parts: Vec<[f64; 2]> = Vec::new();
        for &[fprm, lprm] in &parts {
            let is_circle_ellipse = matches!(
                line.line_type,
                IntPatchIType::Circle | IntPatchIType::Ellipse
            );
            if is_circle_ellipse && fprm < a_nul && lprm > a_nul {
                let mut fprm = fprm;
                let mut lprm = lprm;
                while fprm < a_nul || fprm > a_period {
                    fprm += a_period;
                }
                while lprm < a_nul || lprm > a_period {
                    lprm += a_period;
                }
                // OCCT L955-975: the [fprm, 2*PI] sub-interval.
                if (a_period - fprm) > a_tol_pc {
                    a_parts.push([fprm, a_period]);
                } else {
                    let p1 = line.curve.point_at(fprm);
                    let p2 = line.curve.point_at(a_period);
                    if p1.distance(p2) > tol {
                        a_parts.push([fprm, a_period]);
                    }
                }
                // OCCT L976-996: the [0, lprm] sub-interval.
                if (lprm - a_nul) > a_tol_pc {
                    a_parts.push([a_nul, lprm]);
                } else {
                    let p1 = line.curve.point_at(a_nul);
                    let p2 = line.curve.point_at(lprm);
                    if p1.distance(p2) > tol {
                        a_parts.push([a_nul, lprm]);
                    }
                }
            } else {
                a_parts.push([fprm, lprm]);
            }
        }
        // OCCT MakeCurve L776-1846: one curve per part.
        // For an analytic Circle line, anchor each face's pcurve to the frame
        // of the face's existing boundary pcurves: OCCT builds every pcurve of
        // a face in the face's own (BRepAdaptor) parameter frame, so the
        // section pcurve and the boundary pcurves agree.  rcad's raw analytic
        // frame (parameter origin of the circle) may differ by a phase from
        // the stored boundary pcurves, inverting UV classifications.
        let (anchor1, anchor2) = if line.line_type == IntPatchIType::Circle {
            if let Curve3::Circle(c) = &line.curve {
                (
                    circle_pcurve_frame_anchor(ds, f1, c),
                    circle_pcurve_frame_anchor(ds, f2, c),
                )
            } else {
                (None, None)
            }
        } else {
            (None, None)
        };
        for &[fprm, lprm] in &a_parts {
            // rcad note: the box-based classify_point clip cannot represent a
            // face hole (an inner wire); OCCT's TopolTool classifier rejects
            // points inside holes at this stage.  Drop parts whose midpoint
            // lies in a hole of either face (OCCT GeomInt_LineConstructor
            // midpoint classification, L164-179).
            if part_in_face_hole(ds, f1, surf1, &line, fprm, lprm)
                || part_in_face_hole(ds, f2, surf2, &line, fprm, lprm)
            {
                continue;
            }
            let ic = make_part_curve(
                surf1,
                uv1,
                surf2,
                uv2,
                tol,
                approx,
                approx1,
                approx2,
                tol_approx,
                &line,
                anchor1,
                anchor2,
                fprm,
                lprm,
                a_parts.len(),
            );
            if let Some(ic) = ic {
                out.push(ic);
            }
        }
    }
    out
}

/// Parameter frame (u0, d) of the face's stored boundary pcurve for an
/// analytic circle identical to `c` — the anchor that re-frames a section
/// curve's pcurve into the face's existing parameter frame.
///
/// OCCT keeps every pcurve of a face in the face's own BRepAdaptor parameter
/// frame: the section pcurve of a frame (IntTools_FaceFace::MakeCurve /
/// Recadre) and the face's boundary pcurves (ProjLib) are expressed in the
/// same frame, so UV classifications of the section agree with the boundary.
/// rcad's raw analytic frame may be phase-shifted from the stored boundary
/// pcurves' frame; the anchor is recovered by locating the face's boundary
/// edge whose 3D circle coincides with `c` and reading its stored pcurve.
fn circle_pcurve_frame_anchor(
    ds: &DS,
    fi: usize,
    c: &rcad_kernel::geom::Circle3,
) -> Option<(f64, f64)> {
    let face = &*ds.shapes.get(fi)?.shape.data;
    let face = match face {
        TShape::Face(fd) => fd,
        _ => return None,
    };
    let face_ptr = ds.shapes[fi].shape.ptr_id();
    let mut wires: Vec<Shape> = vec![face.outer_wire.clone()];
    wires.extend(face.inner_wires.iter().cloned());
    for w in &wires {
        let TShape::Wire(wd) = &*w.data else { continue };
        for we in &wd.edges {
            // DS index via the (ptr, loc) map — Shape.index is the ORIGINAL
            // BRep index, not a DS index.
            let Some(&ei) = ds.map_shape_index.get(&(we.ptr_id(), we.location)) else {
                continue;
            };
            let Some(si) = ds.shapes.get(ei) else { continue };
            let Some(ed) = si.shape.as_edge() else { continue };
            // The edge's 3D curve must be the same circle (center, radius,
            // axis within tolerance).
            let Some(Curve3::Circle(ec)) = &ed.curve else { continue };
            let same_circle = (ec.center - c.center).length() <= 1e-7
                && (ec.radius - c.radius).abs() <= 1e-7
                && ec.normal.dot(c.normal).abs() >= 1.0 - 1e-9;
            if !same_circle {
                continue;
            }
            // The stored pcurve on this face: keyed by (face ptr, loc id).
            for (k, (pc, _, _)) in &ed.pcurves {
                if k.0 != face_ptr {
                    continue;
                }
                let Curve2d::Line(l) = pc else { continue };
                // The pcurve maps the edge's own parameter t_e to
                // u = u0e + de * t_e.  Relate the section parameter t to
                // t_e: both circles share center/radius/axis; the phase is
                // the azimuth of the section's x_dir in the edge's frame.
                let zc_e = ec.x_dir.cross(ec.y_dir);
                let same_handed = zc_e.dot(c.normal) > 0.0;
                // delta = azimuth of c.x_dir measured in (ec.x_dir, ec.y_dir).
                let delta = c.x_dir.dot(ec.y_dir).atan2(c.x_dir.dot(ec.x_dir));
                // u(t) for the section:
                //   same handed:  t_e = t + delta -> u = u0e + de*(t + delta)
                //   opposite:     t_e = delta - t -> u = u0e + de*(delta - t)
                let (u0a, da) = if same_handed {
                    (l.origin.x + l.direction.x * delta, l.direction.x)
                } else {
                    (l.origin.x + l.direction.x * delta, -l.direction.x)
                };
                return Some((u0a, da));
            }
        }
    }
    None
}

/// True when the midpoint of the curve part `[fprm, lprm]` lies inside one of
/// the face's inner (hole) loops.
///
/// OCCT IntTools_FaceFace::MakeCurve clips intersection curves with the face
/// TopolTool, whose BRepClass_FaceClassifier rejects points inside the holes;
/// the rcad box-based `classify_point` clip cannot represent a hole, so this
/// restores the hole rejection at MakeCurve time (GeomInt_LineConstructor
/// classifies each part's midpoint, OCCT L164-179).
fn part_in_face_hole(
    ds: &DS,
    fi: usize,
    surf: &Surface3,
    line: &IntPatchLine,
    fprm: f64,
    lprm: f64,
) -> bool {
    let inner_wires = match &*ds.shapes[fi].shape.data {
        TShape::Face(fd) => fd.inner_wires.clone(),
        _ => return false,
    };
    if inner_wires.is_empty() {
        return false;
    }
    let p3d = line.curve.point_at(0.5 * (fprm + lprm));
    let (uv, _) = crate::bop::closest_point_on_surface(surf, p3d);
    for ws in &inner_wires {
        let wi = ds.index(ws);
        if wi < 0 || wi as usize >= ds.nb_shapes() {
            continue;
        }
        let wi = wi as usize;
        if ds.shapes[wi].shape_type != ShapeType::Wire {
            continue;
        }
        let wire_edges = match &*ds.shapes[wi].shape.data {
            TShape::Wire(w) => w.edges.clone(),
            _ => Vec::new(),
        };
        let mut poly: Vec<DVec2> = Vec::new();
        for eshape in &wire_edges {
            let ei = ds.index(eshape);
            if ei < 0 || ei as usize >= ds.nb_shapes() {
                continue;
            }
            let edge_data = match &*ds.shapes[ei as usize].shape.data {
                TShape::Edge(ed) => ed,
                _ => continue,
            };
            let Some(c3d) = edge_data.curve.clone() else { continue };
            let (t0, t1) = (edge_data.range[0], edge_data.range[1]);
            const NS: usize = 8;
            for k in 0..=NS {
                let t = t0 + (t1 - t0) * (k as f64 / NS as f64);
                let p = c3d.point_at(t);
                let (uvp, _) = crate::bop::closest_point_on_surface(surf, p);
                poly.push(uvp);
            }
        }
        if poly.len() >= 3 && point_in_uv_poly(uv, &poly) {
            return true;
        }
    }
    false
}

/// Even-odd ray-cast point-in-polygon test on the face UV plane.
fn point_in_uv_poly(p: DVec2, poly: &[DVec2]) -> bool {
    let mut inside = false;
    let mut j = poly.len() - 1;
    for i in 0..poly.len() {
        let (pi, pj) = (poly[i], poly[j]);
        if (pi.y > p.y) != (pj.y > p.y)
            && p.x < (pj.x - pi.x) * (p.y - pi.y) / (pj.y - pi.y) + pi.x
        {
            inside = !inside;
        }
        j = i;
    }
    inside
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
    let other_pcurve = build_projected_pcurve(other_surf, &curve3, tf, tl, 64, other_uv);
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
            pave_blocks: Vec::new(),
            bbox: None,
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

/// OCCT GeomInt_IntSS::IntersectCurveAndBoundary (GeomInt_IntSS_1.cxx L1450+) —
/// parameters where the 2D pcurve crosses the 4 edges of the UV rectangle.
fn collect_boundary_crossings(pc: &Curve2d, uv: [f64; 4], params: &mut Vec<f64>) {
    let [u_min, u_max, v_min, v_max] = uv;
    match pc {
        Curve2d::Line(l) => {
            // Intersect the pcurve with each of the 4 boundary lines.
            let boundaries = [
                ([u_min, v_min], DVec2::new(0.0, 1.0)), // U=Umin, V varies
                ([u_max, v_min], DVec2::new(0.0, 1.0)), // U=Umax, V varies
                ([u_min, v_min], DVec2::new(1.0, 0.0)), // V=Vmin, U varies
                ([u_min, v_max], DVec2::new(1.0, 0.0)), // V=Vmax, U varies
            ];
            for (o, d) in boundaries {
                if let Some(t) =
                    line_line2d_intersection(&l, &Line2d { origin: DVec2::new(o[0], o[1]), direction: d })
                {
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
        Curve2d::Circle(c) => {
            let cx = c.center.x;
            let cy = c.center.y;
            let r = c.radius;
            let a_tol = 10.0 * rcad_kernel::precision::CONFUSION;
            // Vertical edges U = u_min / u_max.
            for x0 in [u_min, u_max] {
                let rr = r * r - (x0 - cx) * (x0 - cx);
                if rr >= -a_tol {
                    let dy = rr.max(0.0).sqrt();
                    for y in [cy + dy, cy - dy] {
                        if y >= v_min - 1e-9 && y <= v_max + 1e-9 {
                            // Parameter on the circle (angle in its frame).
                            let d = DVec2::new(x0 - cx, y - cy);
                            let t = d.dot(c.y_dir).atan2(d.dot(c.x_dir));
                            params.push(t);
                        }
                    }
                }
            }
            // Horizontal edges V = v_min / v_max.
            for y0 in [v_min, v_max] {
                let rr = r * r - (y0 - cy) * (y0 - cy);
                if rr >= -a_tol {
                    let dx = rr.max(0.0).sqrt();
                    for x in [cx + dx, cx - dx] {
                        if x >= u_min - 1e-9 && x <= u_max + 1e-9 {
                            let d = DVec2::new(x - cx, y0 - cy);
                            let t = d.dot(c.y_dir).atan2(d.dot(c.x_dir));
                            params.push(t);
                        }
                    }
                }
            }
        }
        _ => {}
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

/// OCCT GeomInt_IntSS::BuildPCurves (GeomInt_IntSS_1.cxx L1172-1304) via
/// GeomProjLib::Curve2d -> ProjLib_ProjectedCurve: analytic projection of a 3D
/// curve onto the other surface.  Two cases are solved exactly:
/// - Circle on a plane (ProjLib_Plane::Project): the 2D image is a circle with
///   the same parameterization (the plane's orthonormal UV frame preserves the
///   circle frame);
/// - Circle on a sphere (ProjLib_Sphere::Project(Circle), ProjLib_Sphere.cxx
///   L97-180): a meridian (isIsoU) or latitude (isIsoV) circle.  The OCCT
///   analytic line is only same-parameter near the frame origin — for a
///   general arc it deviates and BRepLib::SameParameter (BRepLib.cxx L1237+,
///   called by BOPTools_AlgoTools::MakePCurve L1723) re-fits it as a BSpline.
///   rcad has no SameParameter, so the same-parameter line is built directly
///   from the arc (below) and the midpoint deviation check returns None so the
///   make_pcurves stage builds the sampled BSpline, matching the OCCT final
///   pcurve type.
/// Returns None when rcad has no analytic projection — the caller then falls
/// back to the sampling/projection path, as OCCT's BOPAlgo_MPC does for null
/// pcurves.
pub(crate) fn build_analytic_pcurve(
    other_surf: &Surface3,
    curve3: &Curve3,
    tf: f64,
    tl: f64,
    uv_bounds: [f64; 4],
    pcurve_anchor: Option<(f64, f64)>,
) -> Option<Curve2d> {
    if let (Curve3::Circle(c), Surface3::Plane(pl)) = (curve3, other_surf) {
        let d = c.center - pl.origin;
        let c2 = DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir));
        // Projected frame: vx/vy are the images of the circle x/y axes scaled by
        // the radius, so P2d(t) = c2 + cos(t)*vx + sin(t)*vy matches the 3D
        // circle's parameterization.
        let vx = DVec2::new(
            c.radius * c.x_dir.dot(pl.u_dir),
            c.radius * c.x_dir.dot(pl.v_dir),
        );
        let vy = DVec2::new(
            c.radius * c.y_dir.dot(pl.u_dir),
            c.radius * c.y_dir.dot(pl.v_dir),
        );
        let lx = vx.length();
        let ly = vy.length();
        if lx > 1e-12
            && ly > 1e-12
            && vx.dot(vy).abs() < 1e-12 * lx * ly
            && (lx - ly).abs() < 1e-9 * lx.max(1.0)
        {
            // Projected circle: same-parameter with the 3D circle.
            return Some(Curve2d::Circle(Circle2d {
                center: c2,
                x_dir: vx / lx,
                y_dir: vy / ly,
                radius: (lx + ly) * 0.5,
            }));
        }
    }
    if let (Curve3::Circle(c), Surface3::Cylinder(cyl)) = (curve3, other_surf) {
        // OCCT GeomInt_IntSS::BuildPCurves -> ProjLib_Cylinder::Project: a
        // circle lying on the cylinder (center on the axis, plane normal
        // parallel to the axis, radius equal) is an iso-parameter line
        // v=const; its pcurve is the line u(t) = u0 +/- t.
        let x_ax = cyl.ref_dir.normalize_or_zero();
        let y_ax = cyl.axis.cross(x_ax).normalize_or_zero();
        if x_ax.length_squared() < 0.5 || y_ax.length_squared() < 0.5 {
            return None;
        }
        let zc = c.x_dir.cross(c.y_dir);
        let r0 = c.center - cyl.origin;
        let tol = 1e-7;
        let on_axis = r0.dot(x_ax).abs() <= tol && r0.dot(y_ax).abs() <= tol;
        let parallel = zc.dot(cyl.axis).abs() > 1.0 - 1e-9;
        if on_axis && parallel && (c.radius - cyl.radius).abs() <= tol {
            let v = r0.dot(cyl.axis);
            // u0: the azimuth of the circle's t=0 point; u(t) = u0 + t for a
            // positively-oriented circle (normal parallel to the axis),
            // u(t) = u0 - t otherwise.
            let mut u0 = c.x_dir.dot(y_ax).atan2(c.x_dir.dot(x_ax));
            if u0 < 0.0 {
                u0 += std::f64::consts::TAU;
            }
            let d = if zc.dot(cyl.axis) > 0.0 { 1.0 } else { -1.0 };
            // With a frame anchor from the face's existing boundary pcurves,
            // the whole pcurve is expressed in that frame: u(t) = u0a + da*t.
            if let Some((u0a, da)) = pcurve_anchor {
                return Some(Curve2d::Line(Line2d::new(
                    DVec2::new(u0a, v),
                    DVec2::new(da, 0.0),
                )));
            }
            let _ = uv_bounds;
            return Some(Curve2d::Line(Line2d::new(
                DVec2::new(u0, v),
                DVec2::new(d, 0.0),
            )));
        }
    }
    if let (Curve3::Circle(c), Surface3::Sphere(sp)) = (curve3, other_surf) {
        // OCCT ProjLib_Sphere::Project(gp_Circ) (L97-180).
        let xs = sp.ref_dir.normalize();
        let ys = sp.axis.cross(xs).normalize();
        let zs = sp.axis.normalize();
        // Xc/Yc/Zc: the circle's local frame; O: the sphere location.
        let xc = c.x_dir;
        let yc = c.y_dir;
        let zc = xc.cross(yc);
        // Precision::Confusion() is the tolerance of IsNormal / IsEqual.
        let tol = CONFUSION;
        // isIsoU = Zc.IsNormal(Zs, Tol) && O.IsEqual(C.Location(), Tol)
        let is_iso_u = zc.dot(zs).abs() <= tol
            && (sp.center.x - c.center.x).abs() <= tol
            && (sp.center.y - c.center.y).abs() <= tol
            && (sp.center.z - c.center.z).abs() <= tol;
        let mut line: Option<Line2d> = None;
        if is_iso_u {
            // The circle is a meridian (the arc passes through both poles).
            // OCCT's analytic line (ProjLib_Sphere.cxx L124-179, u=const
            // through the frame points) is NOT same-parameter for a general
            // arc — its V slope is -1 through the frame while the true V
            // slope is +1 on the far side of a pole, and its U origin is the
            // frame's, not the arc's.  BRepLib::SameParameter (BRepLib.cxx
            // L1237+, called by BOPTools_AlgoTools::MakePCurve L1723) re-fits
            // it; rcad has no SameParameter, so the same-parameter line is
            // built directly from the arc's endpoint UVs (u = const through
            // the arc midpoint, V linear in t).
            let tm = 0.5 * (tf + tl);
            let u_const = quadric_uv_params(other_surf, curve3.point_at(tm))?.x;
            let uv0 = quadric_uv_params(other_surf, curve3.point_at(tf))?;
            let uv1 = quadric_uv_params(other_surf, curve3.point_at(tl))?;
            if !uv0.is_finite() || !uv1.is_finite() {
                return None;
            }
            let d_v = (uv1.y - uv0.y) / (tl - tf);
            if !d_v.is_finite() {
                return None;
            }
            let p2d1 = DVec2::new(u_const, uv0.y - tf * d_v);
            line = Some(Line2d::new(p2d1, DVec2::new(0.0, d_v)));
        }
        // isIsoV = Xc.IsNormal(Zs, Tol) && Yc.IsNormal(Zs, Tol)
        let is_iso_v = xc.dot(zs).abs() <= tol && yc.dot(zs).abs() <= tol;
        if is_iso_v {
            // The circle is a latitude circle: the pcurve is a v=const line.
            // U(t) = atan2(C(t).Ys, C(t).Xs) = U0 + (zc.zs)*t with
            // U0 = Xs.AngleWithRef(Xc, Xs^Ys) — the azimuth of the t=0 point
            // (ProjLib_Sphere.cxx L167).  gp_Dir::AngleWithRef
            // (gp_Dir.cxx L55-84): this=Xs, Other=Xc, Vref=Xs^Ys.  XYZ =
            // this^Other; Cosinus = this|Other; Sinus = |XYZ|; Ang =
            // acos(Cosinus) for |Cosinus| < 0.7071, else PI-asin(Sinus) /
            // asin(Sinus) by the Cosinus sign; the result is +Ang when
            // (this^Other)|Vref >= 0 else -Ang.  The acos/asin form is kept
            // (not atan2) because the last-bit difference decides whether the
            // SetInBounds tail wraps u=2*PI back to 0 (ProjLib_Sphere.cxx
            // L244-247), splitting the 4 rim arcs of a latitude circle.
            let mut u = {
                let ref_dir = xs.cross(ys);
                let xyz = xs.cross(xc);
                let cosinus = xs.dot(xc);
                let sinus = xyz.length();
                let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
                    cosinus.acos()
                } else if cosinus < 0.0 {
                    std::f64::consts::PI - sinus.asin()
                } else {
                    sinus.asin()
                };
                if xyz.dot(ref_dir) >= 0.0 { ang } else { -ang }
            };
            // OCCT ElSLib::SphereParameters normalizes U into [0, 2*PI)
            // (normalizeAngle, ElSLib.cxx L1643); the azimuth of the t=0 point
            // must land in the same domain as the seam pcurves.  Without the
            // wrap a latitude arc whose x_dir is exactly antiparallel to the
            // sphere's xs yields atan2(-0.0, -1.0) = -PI, and the WireSplitter's
            // closed-vertex 2D distance check sees u=-PI vs the seam's u=PI
            // (bopfuse_simple ZH7: sphere rotated by 180 deg).
            if u < 0.0 {
                u += std::f64::consts::TAU;
            }
            let z = (c.center - sp.center).dot(zs);
            let v = (z / sp.radius).clamp(-1.0, 1.0).asin();
            let p2d1 = DVec2::new(u, v);
            // D2d = ((Xc ^ Yc).Dot(Xs ^ Ys), 0) — +1 along U when the circle
            // plane normal is parallel to the sphere axis.
            let d2d = DVec2::new(zc.dot(zs), 0.0);
            let mut l = Line2d::new(p2d1, d2d);
            // OCCT ProjLib_Sphere::SetInBounds (ProjLib_Sphere.cxx L203-248),
            // called from ProjLib_ProjectedCurve::Perform (L419) with
            // U = myCurve->FirstParameter(): place the U of the
            // first-parameter point into [0, 2*PI].  A latitude circle's V is
            // constant inside [-PI/2, PI/2], so the Y-wrap (L207-211) and the
            // pole mirror (L213-242) never trigger; only the tail X-wrap
            // (L244-247) applies.  Without it every arc of the same 3D circle
            // projects to the same U0 line and the WireSplitter cannot tell
            // the arcs apart.
            let u_first = l.point_at(tf).x;
            let new_x = crate::geomalgo::int_patch::cycy_common::in_period(
                u_first, 0.0, std::f64::consts::TAU);
            l.origin.x += new_x - u_first;
            line = Some(l);
        }
        if let Some(l) = line {
            // SameParameter net effect: keep the analytic line only when it is
            // same-parameter with the 3D curve (BRepLib::SameParameter checks
            // 22 control points, BRepLib.cxx L1355-1367).  A meridian arc
            // crossing a pole folds V (the line's constant slope then points
            // the wrong way); a latitude arc crossing the seam wraps U.
            // Check the arc's midpoint against the line; a deviation beyond
            // the confusion tolerance means OCCT would re-approximate it
            // (BSpline) — return None so the make_pcurves stage samples.
            let tm = 0.5 * (tf + tl);
            let uv = quadric_uv_params(other_surf, curve3.point_at(tm))?;
            let p_on_line = l.point_at(tm);
            let du = (uv - p_on_line).x;
            let two_pi = std::f64::consts::TAU;
            let du = (du - two_pi * (du / two_pi).round()).abs();
            let dv = (uv - p_on_line).y.abs();
            if du > 1e-7 || dv > 1e-7 {
                return None;
            }
            return Some(Curve2d::Line(l));
        }
    }
    None
}

/// OCCT GeomInt_IntSS::BuildPCurves (GeomInt_IntSS_1.cxx L1172-1304):
/// projects the 3D curve onto the other surface.  For a restriction arc that is
/// a 3D circle lying on a plane (a cylinder/cone base circle on a box face),
/// the exact 2D image on the plane is a circle with the same parameterization
/// (the plane's orthonormal UV frame preserves the circle frame).  Other cases
/// fall back to the small-range branch (a line through the endpoint UVs).
fn build_projected_pcurve(
    other_surf: &Surface3,
    curve3: &Curve3,
    tf: f64,
    tl: f64,
    _n: usize,
    uv_bounds: [f64; 4],
) -> Option<Curve2d> {
    if let Some(c) = build_analytic_pcurve(other_surf, curve3, tf, tl, uv_bounds, None) {
        return Some(c);
    }
    // OCCT BuildPCurves small-range branch: the pcurve is a line segment.
    let p0 = curve3.point_at(tf);
    let p1 = curve3.point_at(tl);
    let uv0 = quadric_uv_params(other_surf, p0)?;
    let uv1 = quadric_uv_params(other_surf, p1)?;
    if !uv0.is_finite() || !uv1.is_finite() {
        return None;
    }
    Some(Curve2d::Line(Line2d {
        origin: uv0,
        direction: (uv1 - uv0) / (tl - tf),
    }))
}

/// Classify a WLine point on a surface using its precomputed UV.
/// OCCT GeomInt_LineConstructor WLine path: classify with the constructor
/// tolerance `Tol = Precision::PConfusion() * 35.0` (L118).
fn in_uv_rect_adjusted(surf: &Surface3, rect: [f64; 4], _p3d: DVec3, u: f64, v: f64) -> bool {
    let adj = adjust_periodic_uv(surf, DVec2::new(u, v), rect);
    in_uv_rect(adj, rect, PCONFUSION * 35.0)
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
    // OCCT L152-328: the WLine has its own path (vertices carry integer point
    // indices).
    if line.is_wline() {
        return line_constructor_wline_parts(surf1, uv1, surf2, uv2, line);
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
            let p3d = line.curve.point_at(pmid);
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

/// OCCT-aligned: GeomInt_LineConstructor::Perform WLine path (L152-328).
/// WLine vertices carry integer point indices; consecutive pairs are
/// classified on the stored UVs.
fn line_constructor_wline_parts(
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    line: &IntPatchLine,
) -> Vec<[f64; 2]> {
    let mut result: Vec<[f64; 2]> = Vec::new();
    let nbvtx = line.vertices.len();
    for i in 0..nbvtx.saturating_sub(1) {
        let firstp = line.vertices[i].param_on_line;
        let lastp = line.vertices[i + 1].param_on_line;
        // OCCT L162: if (firstp != lastp) — exact inequality (the GLine path
        // L354 uses abs > PConfusion, the WLine path does not).
        if firstp != lastp {
            // OCCT L163: if (lastp != firstp + 1) — any non-adjacent pair (incl.
            // fractional parameter gaps) uses the midpoint-classification branch.
            if (lastp - firstp) != 1.0 {
                // OCCT L164-179: non-adjacent vertices — classify the midpoint
                // polyline point.  OCCT L166-167: pmid = (int)((firstp+lastp)/2),
                // then WLine->Point(pmid) — the actual polyline point (1-based).
                let pmid = ((firstp + lastp) / 2.0) as usize;
                let p = line.wline_pnts[pmid - 1];
                let in1 = in_uv_rect_adjusted(surf1, uv1, p.p3d, p.u1, p.v1);
                if in1 {
                    let in2 = in_uv_rect_adjusted(surf2, uv2, p.p3d, p.u2, p.v2);
                    if in2 {
                        result.push([firstp, lastp]);
                    }
                }
            } else if line.wl_type == WLineType::ImpPrm {
                // OCCT L183-225: the implicit-parametric intersector does not
                // respect the quadric domain; classify the interpolated midpoint
                // of the two endpoint points (OCCT L203-204: Point((int)firstp),
                // Point((int)lastp) — the actual polyline points, 1-based).
                // OCCT L205-213: AdjustPeriodic each endpoint's params on both
                // surfaces FIRST, then average the adjusted params.
                let pf = line.wline_pnts[firstp as usize - 1];
                let pl = line.wline_pnts[lastp as usize - 1];
                let a_uvf1 = adjust_periodic_uv(surf1, DVec2::new(pf.u1, pf.v1), uv1);
                let a_uvf2 = adjust_periodic_uv(surf2, DVec2::new(pf.u2, pf.v2), uv2);
                let a_uvl1 = adjust_periodic_uv(surf1, DVec2::new(pl.u1, pl.v1), uv1);
                let a_uvl2 = adjust_periodic_uv(surf2, DVec2::new(pl.u2, pl.v2), uv2);
                let mu1 = 0.5 * (a_uvf1.x + a_uvl1.x);
                let mv1 = 0.5 * (a_uvf1.y + a_uvl1.y);
                let mu2 = 0.5 * (a_uvf2.x + a_uvl2.x);
                let mv2 = 0.5 * (a_uvf2.y + a_uvl2.y);
                let pmid = 0.5 * (pf.p3d + pl.p3d);
                let in1 = in_uv_rect(DVec2::new(mu1, mv1), uv1, PCONFUSION * 35.0);
                if in1 {
                    let in2 = in_uv_rect(DVec2::new(mu2, mv2), uv2, PCONFUSION * 35.0);
                    if in2 {
                        result.push([firstp, lastp]);
                    }
                }
            } else {
                // OCCT L226-252: both endpoint points must be inside both domains
                // (OCCT L228/237: Point((int)firstp), Point((int)lastp)).
                let pf = line.wline_pnts[firstp as usize - 1];
                let in1 = in_uv_rect_adjusted(surf1, uv1, pf.p3d, pf.u1, pf.v1);
                if in1 {
                    let in2 = in_uv_rect_adjusted(surf2, uv2, pf.p3d, pf.u2, pf.v2);
                    if in2 {
                        let pl = line.wline_pnts[lastp as usize - 1];
                        let in3 = in_uv_rect_adjusted(surf1, uv1, pl.p3d, pl.u1, pl.v1);
                        if in3 {
                            let in4 = in_uv_rect_adjusted(surf2, uv2, pl.p3d, pl.u2, pl.v2);
                            if in4 {
                                result.push([firstp, lastp]);
                            }
                        }
                    }
                }
            }
        }
    }
    // OCCT L257-326: when the pair count is > 1 and one surface is a Plane and
    // the other a SurfaceOfExtrusion/Revolution, collapse contiguous pairs that
    // share an integer point index into one interval.
    let a_nb_parts = result.len();
    if a_nb_parts > 1 {
        let mut b_cond = false;
        if matches!(surf1, Surface3::Plane(_)) {
            if matches!(surf2, Surface3::LinearExtrusion(_) | Surface3::Revolution(_)) {
                b_cond = !b_cond;
            }
        } else if matches!(surf2, Surface3::Plane(_)) {
            if matches!(surf1, Surface3::LinearExtrusion(_) | Surface3::Revolution(_)) {
                b_cond = !b_cond;
            }
        }
        if b_cond {
            // OCCT L293-312: walk the flat seqp; a param whose integer index has
            // not been seen yet is appended, otherwise the last appended param is
            // removed (dropping the shared endpoint of two contiguous pairs).
            let mut a_map: std::collections::HashSet<i64> = std::collections::HashSet::new();
            let mut a_seq_tmp: Vec<f64> = Vec::new();
            for pair in result.iter() {
                for &lastp in pair.iter() {
                    let an_index = lastp as i64;
                    if a_map.insert(an_index) {
                        a_seq_tmp.push(lastp);
                    } else {
                        a_seq_tmp.pop();
                    }
                }
            }
            // OCCT L314-324: rebuild the pairs from the merged flat sequence.
            result.clear();
            let a_nb = a_seq_tmp.len() / 2;
            for i in 0..a_nb {
                let firstp = a_seq_tmp[2 * i];
                let lastp = a_seq_tmp[2 * i + 1];
                result.push([firstp, lastp]);
            }
        }
    }
    // OCCT L342-361: a WLine whose vertex intervals are all rejected produces no
    // parts (there is no "keep the full range a priori" fallback).  The start
    // vertex required for a non-empty interval is inserted by IntPatch_WLine::
    // ComputeVertexParameters, which is ported (compute_vertex_parameters_wline).
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
    // OCCT GeomInt_LineConstructor::TreatCircle keeps the vertex parameters as
    // they are (no wrapping); the MakeCurve circle branch splits the intervals
    // crossing 0 into two.
    let mut params: Vec<f64> = if line.vertices.is_empty() {
        vec![0.0, two_pi]
    } else {
        line.vertices.iter().map(|v| v.param_on_line).collect()
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

/// OCCT IntTools_Tools::CurveTolerance (L430-464) + ParabolaTolerance
/// (L470-561): the tolerance of a trimmed parabola.  The endpoints (at `tf`/`tl`)
/// are projected onto the parabola X-axis (the line through the vertex along
/// `axis_dir`); the tolerance grows with the axis distance of the farther end.
fn parabola_curve_tolerance(p: &Parabola3, tf: f64, tl: f64, base_tol: f64) -> f64 {
    let focal = p.focal_param * 0.5;
    if focal == 0.0 {
        return base_tol;
    }
    let axis = p.axis_dir.normalize_or_zero();
    let vertex = p.vertex;
    let proj_x = |t: f64| (p.point_at(t) - vertex).dot(axis);
    let mut tol1 = base_tol;
    let x1 = proj_x(tf);
    if x1 >= 0.0 {
        tol1 = base_tol * (0.5 * x1 / focal).sqrt();
    }
    if tol1 == 0.0 {
        tol1 = base_tol;
    }
    let mut tol2 = base_tol;
    let x2 = proj_x(tl);
    if x2 >= 0.0 {
        tol2 = base_tol * (0.5 * x2 / focal).sqrt();
    }
    if tol2 == 0.0 {
        tol2 = base_tol;
    }
    tol1.max(tol2)
}

/// OCCT MakeCurve L776-1846: build one IntersectionCurve for a LineConstructor
/// part.  The full analytic curve is kept; `t_range` holds the clipped part.
fn make_part_curve(
    surf1: &Surface3,
    uv1: [f64; 4],
    surf2: &Surface3,
    uv2: [f64; 4],
    tol: f64,
    approx: bool,
    approx1: bool,
    approx2: bool,
    tol_approx: f64,
    line: &IntPatchLine,
    pcurve_anchor1: Option<(f64, f64)>,
    pcurve_anchor2: Option<(f64, f64)>,
    fprm: f64,
    lprm: f64,
    a_nb_parts: usize,
) -> Option<IntersectionCurve> {
    let curve = &line.curve;
    match line.line_type {
        IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
            // OCCT L829-833: bFNIt = Precision::IsNegativeInfinite(fprm),
            // bLPIt = Precision::IsPositiveInfinite(lprm).
            let b_fn = rcad_kernel::precision::is_negative_infinite_value(fprm);
            let b_lp = rcad_kernel::precision::is_positive_infinite_value(lprm);
            // OCCT L809: if (!bFNIt && !bLPIt) — both bounds finite: build the
            // trimmed curve without any reference-point test (L810-848).
            if !b_fn && !b_lp {
                // (the IntersectionCurve below carries the clipped t_range)
            } else {
                // OCCT L850-897: at least one bound is infinite — test a
                // reference point on both face domains before keeping the
                // curve (untrimmed).
                let d_t = 100.0;
                let test_t = if b_fn && !b_lp {
                    lprm - d_t
                } else if !b_fn && b_lp {
                    fprm + d_t
                } else {
                    // OCCT L868: IntTools_Tools::IntermediatePoint(-dT, dT).
                    intermediate_point(-d_t, d_t)
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
            // OCCT MakeCurve L816-820: a parabola's curve tolerance is not the
            // base tolerance — IntTools_Tools::CurveTolerance grows with the
            // distance of the trimmed endpoints from the vertex.
            let tolerance = if line.line_type == IntPatchIType::Parabola && !b_fn && !b_lp {
                if let Curve3::Parabola(p) = curve {
                    parabola_curve_tolerance(p, fprm, lprm, tol)
                } else {
                    tol
                }
            } else {
                tol
            };
            Some(IntersectionCurve {
                curve: curve.clone(),
                t_range: [fprm, lprm],
                pcurve1: None,
                pcurve2: None,
                tolerance,
                tang_tolerance: line.tang_tolerance,
                pave_blocks: Vec::new(),
            bbox: None,
            })
        }
        IntPatchIType::Circle | IntPatchIType::Ellipse => {
            let a_period = std::f64::consts::TAU;
            let a_real_eps = f64::EPSILON;
            // OCCT L1028-1061: the pcurves on both faces come from
            // GeomInt_IntSS::BuildPCurves (analytic for Circle on a plane / a
            // sphere meridian or latitude); rcad computes them analytically and
            // leaves them null otherwise (the make_pcurves stage projects).
            let pcurve1 = build_analytic_pcurve(surf1, curve, fprm, lprm, uv1, pcurve_anchor1);
            let pcurve2 = build_analytic_pcurve(surf2, curve, fprm, lprm, uv2, pcurve_anchor2);
            let is_full_period = fprm.abs() <= a_real_eps && (lprm - a_period).abs() <= a_real_eps;
            if !is_full_period {
                if lprm <= fprm + 1e-12 {
                    return None;
                }
                Some(IntersectionCurve {
                    curve: curve.clone(),
                    t_range: [fprm, lprm],
                    pcurve1,
                    pcurve2,
                    tolerance: tol,
                    tang_tolerance: line.tang_tolerance,
                    pave_blocks: Vec::new(),
            bbox: None,
                })
            } else if a_nb_parts == 1 {
                // OCCT L1074-1103: accept the full circle.
                Some(IntersectionCurve {
                    curve: curve.clone(),
                    t_range: [fprm, lprm],
                    pcurve1,
                    pcurve2,
                    tolerance: tol,
                    tang_tolerance: line.tang_tolerance,
                    pave_blocks: Vec::new(),
            bbox: None,
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
                        // OCCT L1117-1140: the 18-point branch also sets the
                        // analytic pcurves (BuildPCurves with the UV bounds).
                        return Some(IntersectionCurve {
                            curve: curve.clone(),
                            t_range: [fprm, lprm],
                            pcurve1: build_analytic_pcurve(surf1, curve, fprm, lprm, uv1, pcurve_anchor1),
                            pcurve2: build_analytic_pcurve(surf2, curve, fprm, lprm, uv2, pcurve_anchor2),
                            tolerance: tol,
                            tang_tolerance: line.tang_tolerance,
                            pave_blocks: Vec::new(),
            bbox: None,
                        });
                    }
                }
                None
            }
        }
        IntPatchIType::Walking => {
            // OCCT MakeCurve L1175-1391 (IntPatch_Walking): when myApprox is
            // true the part becomes the smooth GeomInt_WLApprox BSpline
            // (ApproxInt_Approx + ApproxInt_KnotTools + AppParCurves,
            // ported in geomalgo::approx_int); otherwise a degree-1 BSpline
            // through the polyline points (MakeBSpline L1911).
            if lprm <= fprm + 1e-12 {
                return None;
            }
            let wline = &line.wline_pnts;
            let n = wline.len();
            let ifprm = (fprm as usize).max(1).min(n);
            let ilprm = (lprm as usize).max(1).min(n);
            if ilprm <= ifprm {
                return None;
            }
            if std::env::var("RCAD_PCTRACE").is_ok() && ilprm - ifprm <= 30 {
                eprintln!("[PMC] WLine tail ifprm={} ilprm={} n={}", ifprm, ilprm, n);
                for (k, wp) in wline.iter().enumerate() {
                    eprintln!("  [PMC-PT] {} ({:.15},{:.15},{:.15})", k, wp.p3d.x, wp.p3d.y, wp.p3d.z);
                }
            }
            // OCCT L1316-1326: a Plane surface disables the 3D approximation
            // (the 3D curve is rebuilt from the 2D pcurve).
            let typs1 = matches!(surf1, Surface3::Plane(_));
            let typs2 = matches!(surf2, Surface3::Plane(_));
            let mut an_approx = approx; // myApprox (ToApproxC3d)
            let mut an_approx1 = approx1;
            let mut an_approx2 = approx2;
            if typs1 {
                an_approx = false;
                an_approx1 = true;
            } else if typs2 {
                an_approx = false;
                an_approx2 = true;
            }
            if an_approx {
                // OCCT L1343-1347: ApproxParameters + SetParameters + Perform.
                let (i_deg_min, i_deg_max, i_nb_iter) =
                    approx_parameters_for(surf1, surf2);
                let a_par_type = crate::geomalgo::approx_int::define_par_type(
                    &crate::geomalgo::approx_int::WLineAccess {
                        line,
                        indicemin: ifprm,
                        indicemax: ilprm,
                        nbp3d: 1,
                        nbp2d: 2,
                        approx_u1v1: an_approx1,
                        approx_u2v2: an_approx2,
                        p2d_on_first: true,
                        xo: 0.0,
                        yo: 0.0,
                        zo: 0.0,
                        u1o: 0.0,
                        v1o: 0.0,
                        u2o: 0.0,
                        v2o: 0.0,
                        s1: surf1,
                        s2: surf2,
                        uv1,
                        uv2,
                    },
                    ifprm,
                    ilprm,
                    an_approx,
                    an_approx1,
                    an_approx2,
                );
                let mut app = crate::geomalgo::approx_int::WLineApprox::new();
                app.set_parameters(
                    tol_approx,
                    tol_approx,
                    i_deg_min,
                    i_deg_max,
                    i_nb_iter,
                    30,
                    true,
                    a_par_type,
                );
                app.perform(
                    &crate::geomalgo::approx_int::WLineAccess {
                        line,
                        indicemin: ifprm,
                        indicemax: ilprm,
                        nbp3d: 1,
                        nbp2d: 2,
                        approx_u1v1: an_approx1,
                        approx_u2v2: an_approx2,
                        p2d_on_first: true,
                        xo: 0.0,
                        yo: 0.0,
                        zo: 0.0,
                        u1o: 0.0,
                        v1o: 0.0,
                        u2o: 0.0,
                        v2o: 0.0,
                        s1: surf1,
                        s2: surf2,
                        uv1,
                        uv2,
                    },
                    an_approx,
                    an_approx1,
                    an_approx2,
                    ifprm,
                    ilprm,
                );
                if app.is_done() {
                    let mbspc = app.value();
                    if std::env::var("RCAD_PCTRACE").is_ok() && ilprm - ifprm <= 30 {
                        let mut p3 = Vec::new();
                        mbspc.curve(1, &mut p3);
                        let mid = if p3.len() > 2 { p3[p3.len() / 2] } else { DVec3::ZERO };
                        eprintln!("[PMC] approx done deg={} np={} knots={:?} P0={:?} Pn={:?} midpole={:?}", mbspc.degree, p3.len(), mbspc.knots, p3.first(), p3.last(), mid);
                    }
                    // OCCT L1642-1733 (typs != Plane): curve 1 -> 3D BSpline,
                    // curves 2/3 -> the 2D pcurves (myApprox1/2 gates).  The
                    // resulting curves keep the approximation parameter domain
                    // [0, 1] (OCCT creates the BSpline from mbspc.Knots()).
                    if let Some((curve, pcurve1, pcurve2)) =
                        multibsp_to_curves(&mbspc, an_approx1, an_approx2)
                    {
                        return Some(IntersectionCurve {
                            curve,
                            t_range: [0.0, 1.0],
                            pcurve1,
                            pcurve2,
                            tolerance: tol,
                            tang_tolerance: line.tang_tolerance,
                            pave_blocks: Vec::new(),
                            bbox: None,
                        });
                    }
                }
            }
            // OCCT !myApprox / !IsDone fallback: MakeBSpline / MakeBSpline2d.
            let curve = wline_part_bspline(line, fprm, lprm).unwrap_or_else(|| curve.clone());
            let pcurve1 = wline_part_bspline2d(line, fprm, lprm, true);
            let pcurve2 = wline_part_bspline2d(line, fprm, lprm, false);
            Some(IntersectionCurve {
                curve,
                t_range: [fprm, lprm],
                pcurve1,
                pcurve2,
                tolerance: tol,
                tang_tolerance: line.tang_tolerance,
                pave_blocks: Vec::new(),
            bbox: None,
            })
        }
        _ => None,
    }
}

/// OCCT IntTools_FaceFace::ApproxParameters (IntTools_FaceFace.cxx
/// L2736-2783): degree and iteration parameters of the WLine approximation.
fn approx_parameters_for(surf1: &Surface3, surf2: &Surface3) -> (usize, usize, i32) {
    let mut i_nb_iter = 0i32;
    let mut i_deg_min = 4usize;
    let mut i_deg_max = 8usize;
    // Cylinder/Torus.
    let cyl_rad = match (surf1, surf2) {
        (Surface3::Cylinder(c), Surface3::Torus(t)) => Some((c.radius, t.minor_radius)),
        (Surface3::Torus(t), Surface3::Cylinder(c)) => Some((c.radius, t.minor_radius)),
        _ => None,
    };
    if let Some((a_rc, a_rt)) = cyl_rad {
        let d_r = (a_rc - a_rt).abs();
        if d_r < 1.0e-7 {
            i_deg_max = 6;
        }
    }
    // Cylinder + Cylinder.
    if matches!(surf1, Surface3::Cylinder(_)) && matches!(surf2, Surface3::Cylinder(_)) {
        i_nb_iter = 1;
    }
    (i_deg_min, i_deg_max, i_nb_iter)
}

/// Convert the OCCT MultiBSpCurve (shared degree/knots/multiplicities and
/// per-curve poles) into the rcad 3D curve and the 2D pcurves.  The knot
/// vector keeps the approximation domain (typically [0, 1]) — OCCT
/// IntTools_FaceFace creates the Geom_BSplineCurve directly from
/// mbspc.Knots()/Multiplicities() (IntTools_FaceFace.cxx L1646-1649), and
/// Approx_MCurvesToBSpCurve normalizes the knots to [0, 1].
fn multibsp_to_curves(
    mbspc: &crate::geomalgo::approx_int::MultiBSpCurve,
    with_pc1: bool,
    with_pc2: bool,
) -> Option<(Curve3, Option<Curve2d>, Option<Curve2d>)> {
    if mbspc.poles.is_empty() {
        return None;
    }
    let mut p3 = Vec::new();
    mbspc.curve(1, &mut p3);
    if p3.len() < 2 {
        return None;
    }
    // OCCT IntTools_FaceFace.cxx L1646-1649: the BSpline is created directly
    // from the MultiBSpCurve knots/multiplicities — the parameter domain is
    // the approximation's [0, 1] (Approx_MCurvesToBSpCurve normalizes the
    // knots to [0, 1]); it is NOT remapped to the WLine point-index range
    // [fprm, lprm] (those indices only select the points fed to the
    // approximation).
    let knots = expand_knots(&mbspc.knots, &mbspc.mults);
    let nb3d = mbspc.poles[0].nb_points();
    let nb2d = mbspc.poles[0].nb_points2d();
    let c3 = Curve3::BSpline(BSplineCurve3 {
        degree: mbspc.degree,
        knots: knots.clone(),
        control_points: p3,
        weights: vec![],
        is_periodic: false,
    });
    let pc1 = if with_pc1 && nb2d >= 1 {
        let mut p2 = Vec::new();
        mbspc.curve2d(nb3d + 1, &mut p2);
        Some(Curve2d::BSpline(BSplineCurve2 {
            degree: mbspc.degree,
            knots: knots.clone(),
            control_points: p2,
            weights: vec![],
        }))
    } else {
        None
    };
    let pc2 = if with_pc2 && nb2d >= 2 {
        let mut p2 = Vec::new();
        mbspc.curve2d(nb3d + 2, &mut p2);
        Some(Curve2d::BSpline(BSplineCurve2 {
            degree: mbspc.degree,
            knots,
            control_points: p2,
            weights: vec![],
        }))
    } else {
        None
    };
    Some((c3, pc1, pc2))
}

/// Expand the compressed (knots, multiplicities) representation into the full
/// knot vector, keeping the approximation parameter domain (typically [0, 1]).
fn expand_knots(knots: &[f64], mults: &[usize]) -> Vec<f64> {
    let mut out = Vec::new();
    for (i, &k) in knots.iter().enumerate() {
        let m = if i < mults.len() { mults[i] } else { 1 };
        for _ in 0..m {
            out.push(k);
        }
    }
    out
}

/// OCCT GeomInt_IntSS::MakeBSpline (GeomInt_IntSS.cxx L1452-1469) — build the
/// 3D curve of a WLine part from the polyline points.  The WLine point
/// parameters are 1-based point indices: the part [fprm, lprm] maps to the
/// polyline points wpts[(int)fprm-1 .. (int)lprm-1] (OCCT WL->Point(ifprm..ilprm)
/// with (int) truncation).  rcad emits a degree-1 (piecewise linear) BSpline
/// through those points, parameterized on the part range [fprm, lprm]
/// (point-index space).
fn wline_part_bspline(line: &IntPatchLine, fprm: f64, lprm: f64) -> Option<Curve3> {
    let wpts = &line.wline_pnts;
    let n = wpts.len();
    if n < 2 {
        return None;
    }
    let i0 = (fprm as usize).saturating_sub(1).min(n - 1);
    let i1 = ((lprm as usize).saturating_sub(1)).min(n - 1);
    if i1 <= i0 {
        return None;
    }
    let ctrl: Vec<DVec3> = wpts[i0..=i1].iter().map(|p| p.p3d).collect();
    let m = ctrl.len();
    let mut knots = Vec::with_capacity(m + 2);
    knots.push(fprm);
    knots.push(fprm);
    for j in 1..m - 1 {
        knots.push(fprm + j as f64);
    }
    knots.push(lprm);
    knots.push(lprm);
    Some(Curve3::BSpline(BSplineCurve3 {
        degree: 1,
        knots,
        control_points: ctrl,
        weights: vec![],
        is_periodic: false,
    }))
}

/// OCCT GeomInt_IntSS::MakeBSpline2d (GeomInt_IntSS.cxx L1473-1502) — build the
/// 2D pcurve of a WLine part on one surface from the polyline points' UV
/// parameters.  Same 1-based point-index mapping and [fprm, lprm]
/// parameterization as wline_part_bspline.
fn wline_part_bspline2d(line: &IntPatchLine, fprm: f64, lprm: f64, on_first: bool) -> Option<Curve2d> {
    let wpts = &line.wline_pnts;
    let n = wpts.len();
    if n < 2 {
        return None;
    }
    let i0 = (fprm as usize).saturating_sub(1).min(n - 1);
    let i1 = ((lprm as usize).saturating_sub(1)).min(n - 1);
    if i1 <= i0 {
        return None;
    }
    let ctrl: Vec<DVec2> = wpts[i0..=i1]
        .iter()
        .map(|p| {
            if on_first {
                DVec2::new(p.u1, p.v1)
            } else {
                DVec2::new(p.u2, p.v2)
            }
        })
        .collect();
    let m = ctrl.len();
    let mut knots = Vec::with_capacity(m + 2);
    knots.push(fprm);
    knots.push(fprm);
    for j in 1..m - 1 {
        knots.push(fprm + j as f64);
    }
    knots.push(lprm);
    knots.push(lprm);
    Some(Curve2d::BSpline(BSplineCurve2 {
        degree: 1,
        knots,
        control_points: ctrl,
        weights: vec![],
    }))
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

/// OCCT-aligned: BOPTools_AlgoTools2D::IntermediatePoint
/// (BOPTools_AlgoTools2D.cxx L404-411), identical to
/// IntTools_Tools::IntermediatePoint (IntTools_Tools.cxx L254-259):
///   the parameter dividing the range [aFirst, aLast] at the ratio
///   10 * e^(-PI) = 0.43213918.
pub(crate) fn intermediate_point(a_first: f64, a_last: f64) -> f64 {
    const PAR_T: f64 = 0.43213918;
    (1.0 - PAR_T) * a_first + PAR_T * a_last
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
