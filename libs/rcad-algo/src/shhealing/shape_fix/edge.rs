//! 1:1 translation of OCCT `ShapeFix_Edge`
//! (`TKShHealing/ShapeFix/ShapeFix_Edge.hxx` L17-249 +
//! `ShapeFix_Edge.cxx` L1-957, docket row `ShapeFix_Edge`).
//!
//! Function-count equation (OCCT ShapeFix_Edge.cxx = rcad edge.rs): 19
//! member functions — `ShapeFix_Edge()` / `Projector` / `FixRemovePCurve` x2
//! / `FixRemoveCurve3d` / `FixAddPCurve` x4 / `FixAddCurve3d` /
//! `FixVertexTolerance` x2 / `FixReversed2d` x2 / `FixSameParameter` x2 /
//! `Status` / `Context` / `SetContext` — plus the 2 file statics
//! `TranslatePCurve` (L171-325) and `TempSameRange` (L329-464):
//! 21 OCCT functions = 21 rcad functions.
//!
//! Architecture bridges (numbered, referenced by the methods below):
//! 1. `BRep` pool argument — OCCT `BRep_Tool`/`BRep_Builder` read and write
//!    TShapes through global accessors; the rcad equivalents take
//!    `brep: &BRep`/`&mut BRep` (the W1/W2 precedent).
//! 2. `TopLoc_Location` -> `u32` (the `BRep.locations` table index, 0 =
//!    identity); pcurve keys compose as `compose_pcurve_location`; the
//!    (surface, location) resolution walks the faces registered on the
//!    surface (the shape_build/edge.rs precedent).
//! 3. `Geom_Curve`/`Geom2d_Curve`/`Geom_Surface` handles -> `Curve3` /
//!    `Curve2d`/`Surface3` values; a handle mutation (`c2d->Reverse()`)
//!    writes the value back into the stored representation.
//! 4. OCCT overloads are disambiguated with `_face`/`_surface` suffixes (and
//!    `_sa` when the `ShapeAnalysis_Surface` overload is taken).
//! 5. OCCT `try { OCC_CATCH_SIGNALS } catch (Standard_Failure)` blocks — the
//!    failure arm sets FAIL2; rcad raises no exceptions in these bodies, so
//!    the arm stays annotated (unreachable) and the flow is preserved.
//! 6. `BRepLib::SameParameter(edge, tolerance)` (ShapeFix_Edge.cxx L850) —
//!    the docket section 4 gap 3 kernel completion item; the annotated
//!    carrier `shape_fix_gap_deps::brep_lib_same_parameter_edge` stands in.
//! 7. `GeomLib::SameRange` and `Geom2d_BezierCurve::Segment` (TempSameRange
//!    L395-455) — TKGeomBase leaves; GAP fns at the bottom of this file
//!    keep the call shape and the documented deviation.

use glam::DVec2;
use rcad_kernel::geom::{
    reverse_curve2d, transform_surface, translate_curve2d, Curve2d, Curve2dEval, CurveEval,
    Line2d, Surface3, TrimmedCurve2, BSplineCurve2,
};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{
    compose_pcurve_location, surface_same, BRepBuilder, BRepTool, CurveRepresentation, ShapeType,
    TShape,
};
use rcad_kernel::BRep;

use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_construct::gap_deps::ShapeAnalysisSurface as ProjectorShapeAnalysisSurface;
use crate::shhealing::shape_construct::project_curve_on_surface::ProjectCurveOnSurface;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::shape_fix::top_exp_vertices;
use crate::shhealing::shape_fix::shape_fix_gap_deps::brep_lib_same_parameter_edge;
use crate::shhealing::shape_fix::shape_tolerance::ShapeFixShapeTolerance;

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (architecture bridge #1; the shape_analysis/edge.rs
// precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Surface(face, L): the face surface and its location.
fn brep_tool_surface_loc(brep: &BRep, face: &Shape) -> (Option<Surface3>, u32) {
    match brep.tshapes[face.index].as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fd.surface_location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::Degenerated(edge).
fn brep_tool_degenerated(edge: &Shape) -> bool {
    matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.degenerated)
}

/// OCCT BRep_Tool::Range(edge, first, last) (BRep_Tool.cxx): the range of
/// the 3D curve when present, else the range of the first geometric
/// representation row.
fn brep_tool_range(brep: &BRep, edge: &Shape) -> Option<(f64, f64)> {
    match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => {
            // The 3D curve row (BRep_Curve3D) carries the TEdgeData range.
            if ed.curve.is_some() {
                return Some((ed.range[0], ed.range[1]));
            }
            // The first GCurve representation row.
            for r in &ed.representations {
                match r {
                    CurveRepresentation::CurveOnSurface { range, .. }
                    | CurveRepresentation::CurveOnClosedSurface { range, .. } => {
                        return Some((range[0], range[1]));
                    }
                    _ => {}
                }
            }
            Some((ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT BRep_Builder::Range(E, F, L, Only3d = true) (BRep_Builder.cxx
/// L1091-1117): sets the range of the 3D curve row only (the
/// `GC->IsCurve3D()` arm).
fn builder_range_only3d(brep: &mut BRep, edge: &Shape, first: f64, last: f64) {
    let ed = brep.edge_mut_inplace(edge.clone());
    ed.range = [first, last];
}

/// OCCT BRep_Builder::Range(E, F, L) (BRep_Builder.cxx L536-556, the
/// Only3d=false default): the range on all representations.
fn builder_range_all(brep: &mut BRep, edge: &Shape, f: f64, l: f64) {
    let ed = brep.edge_mut_inplace(edge.clone());
    ed.range = [f, l];
    for r in ed.representations.iter_mut() {
        match r {
            CurveRepresentation::CurveOnSurface { range, .. }
            | CurveRepresentation::CurveOnClosedSurface { range, .. } => {
                *range = [f, l];
            }
            _ => {}
        }
    }
}

/// The (face ptr, composed location) keys of the edge's pcurve rows whose
/// registered face surface matches (surf, loc) — the surface_face_keys
/// re-host of shape_build/edge.rs (kept private there; this module needs the
/// same resolution for the surface-keyed builder writes, bridge #2).
fn surface_face_keys(brep: &BRep, edge: &Shape, surf: &Surface3, loc: u32) -> Vec<(u64, u32)> {
    let _ = loc;
    let locs = brep.edge_wrapper_locations(edge);
    let mut keys = Vec::new();
    for idx in 0..brep.tshapes.len() {
        if let TShape::Face(fd) = brep.tshapes[idx].as_ref() {
            let matches_surf = fd.surface.as_ref().map_or(false, |s| surface_same(s, surf));
            if !matches_surf {
                continue;
            }
            let face = brep.shape_at(idx);
            for &el in &locs {
                keys.push((
                    face.ptr_id(),
                    compose_pcurve_location(face.location, el, &brep.locations),
                ));
            }
        }
    }
    keys
}

/// OCCT BRep_Builder::UpdateEdge(E, C, S, L, Tol) (BRep_Builder.cxx
/// L679-700): sets the pcurve representation on the (surface, location)
/// pair.
fn builder_update_edge_pcurve_surface(
    brep: &mut BRep,
    edge: &Shape,
    pcurve: &Curve2d,
    surf: &Surface3,
    loc: u32,
    tol: f64,
) {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        ed.pcurves.insert(*k, (pcurve.clone(), 0.0, 0.0));
        ed.representations.push(CurveRepresentation::CurveOnSurface {
            face: *k,
            pcurve: pcurve.clone(),
            range: [0.0, 0.0],
        });
    }
    ed.tolerance = ed.tolerance.max(tol);
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) (BRep_Builder.cxx
/// L702-757): sets the seam (closed-surface) pcurve pair on the
/// (surface, location) pair.
fn builder_update_edge_pcurves_surface(
    brep: &mut BRep,
    edge: &Shape,
    pcurve1: &Curve2d,
    pcurve2: &Curve2d,
    surf: &Surface3,
    loc: u32,
    tol: f64,
) {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        ed.pcurves.insert(*k, (pcurve1.clone(), 0.0, 0.0));
        ed.representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: *k,
                pcurve1: pcurve1.clone(),
                pcurve2: pcurve2.clone(),
                range: [0.0, 0.0],
            });
    }
    ed.tolerance = ed.tolerance.max(tol);
}

/// OCCT BRep_Builder::Range(E, S, L, F, L) (BRep_Builder.cxx L1121-1160):
/// sets the range of the pcurve representation on the (surface, location)
/// pair (OCCT throws Standard_DomainError when no pcurve is registered — the
/// rcad walk skips silently when no key matches).
fn builder_range_on_surface_loc(
    brep: &mut BRep,
    edge: &Shape,
    surf: &Surface3,
    loc: u32,
    first: f64,
    last: f64,
) {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        if let Some(entry) = ed.pcurves.get_mut(k) {
            entry.1 = first;
            entry.2 = last;
        }
    }
    for r in ed.representations.iter_mut() {
        match r {
            CurveRepresentation::CurveOnSurface { face: fk, range, .. }
            | CurveRepresentation::CurveOnClosedSurface { face: fk, range, .. } => {
                if keys.contains(fk) {
                    *range = [first, last];
                }
            }
            _ => {}
        }
    }
}

/// Write-back of a mutated pcurve value into the (surface, location)
/// representation (bridge #3: OCCT `c2d->Reverse()` mutates the shared
/// handle stored in the edge; rcad stores values).
fn set_stored_pcurve(brep: &mut BRep, edge: &Shape, surf: &Surface3, loc: u32, c2d: Curve2d) {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        if let Some(entry) = ed.pcurves.get_mut(k) {
            entry.0 = c2d.clone();
        }
    }
    for r in ed.representations.iter_mut() {
        match r {
            CurveRepresentation::CurveOnSurface { face: fk, pcurve, .. } => {
                if keys.contains(fk) {
                    *pcurve = c2d.clone();
                }
            }
            CurveRepresentation::CurveOnClosedSurface { face: fk, pcurve1, .. } => {
                if keys.contains(fk) {
                    *pcurve1 = c2d.clone();
                }
            }
            _ => {}
        }
    }
}

// ---------------------------------------------------------------------------
// File statics.
// ---------------------------------------------------------------------------

/// OCCT gp_Vec2d::IsParallel (gp_Vec2d.hxx): true when the two vectors are
/// parallel or anti-parallel within the angular tolerance (|cos| within the
/// squared tolerance of 1; zero-magnitude vectors are never parallel).
fn gp_vec2d_is_parallel(v1: DVec2, v2: DVec2, angular_tolerance: f64) -> bool {
    let m1 = v1.length();
    let m2 = v2.length();
    if m1 <= f64::EPSILON || m2 <= f64::EPSILON {
        return false;
    }
    let cos_abs = (v1.dot(v2) / (m1 * m2)).abs();
    cos_abs * cos_abs >= (1.0 - angular_tolerance * angular_tolerance).max(0.0)
}

/// Geom2d_BSplineCurve::FirstParameter — the first distinct knot.
fn bspline2d_first_parameter(b: &BSplineCurve2) -> f64 {
    let n = b.knots.len();
    if n == 0 {
        return 0.0;
    }
    let d = b.degree.min(n - 1);
    b.knots[d]
}

/// Geom2d_BSplineCurve::LastParameter — the last distinct knot.
fn bspline2d_last_parameter(b: &BSplineCurve2) -> f64 {
    let n = b.knots.len();
    if n == 0 {
        return 0.0;
    }
    let d = b.degree.min(n - 1);
    b.knots[n - 1 - d]
}

/// Geom2d_BSplineCurve::Transform(gp_Trsf2d translation) — the pole shift
/// through the translate_curve2d kernel walk.
fn translate_bspline2d(b: &BSplineCurve2, offset: DVec2) -> BSplineCurve2 {
    match translate_curve2d(&Curve2d::BSpline(b.clone()), offset) {
        Curve2d::BSpline(nb) => nb,
        _ => unreachable!("translate_curve2d keeps the BSpline kind"),
    }
}

/// OCCT static TranslatePCurve (ShapeFix_Edge.cxx L171-325): translates a
/// seam pcurve near the surface bounds to the other iso (a U- or V-parallel
/// line is re-anchored, a U- or V-parallel BSpline is shifted by the
/// period).
fn translate_pcurve(a_surf: &Surface3, a_c2d: &Curve2d, a_tol: f64) -> Curve2d {
    // L175-176: aSurf->Bounds(uf, ul, vf, vl) — the Geom_Surface::Bounds
    // range through the W2 analyzer (the raw (u1,u2,v1,v2) form).
    let mut uf = 0.0;
    let mut ul = 0.0;
    let mut vf = 0.0;
    let mut vl = 0.0;
    ShapeAnalysisSurface::new(a_surf.clone()).bounds(&mut uf, &mut ul, &mut vf, &mut vl);

    // L178-180: case of a line.
    if let Curve2d::Line(the_l2d) = a_c2d {
        // L182-183: theLoc = theL2d->Location(); theDir = theL2d->Direction().
        let the_loc = the_l2d.origin;
        let the_dir = the_l2d.direction;

        let mut new_loc;
        let mut the_new_l2d = *the_l2d;

        // L188-200: case UClosed.
        if the_dir.x.abs() <= a_tol && the_dir.y.abs() >= a_tol {
            if (the_loc.x - uf).abs() < (the_loc.x - ul).abs() {
                new_loc = DVec2::new(the_loc.x + (ul - uf), the_loc.y);
            } else {
                new_loc = DVec2::new(the_loc.x - (ul - uf), the_loc.y);
            }
            let _ = &mut new_loc;
            the_new_l2d = Line2d::new(new_loc, the_dir);
        }
        // (the commented-out L201-219 OCCT blocks stay untranslated)
        // L220-232: case VClosed.
        if the_dir.x.abs() >= a_tol && the_dir.y.abs() <= a_tol {
            if (the_loc.y - vf).abs() < (the_loc.y - vl).abs() {
                new_loc = DVec2::new(the_loc.x, the_loc.y + (vl - vf));
            } else {
                new_loc = DVec2::new(the_loc.x, the_loc.y - (vl - vf));
            }
            the_new_l2d = Line2d::new(new_loc, the_dir);
        }
        // (the commented-out L233-250 OCCT blocks stay untranslated)
        // L252-256: TODO Other case not yet implemented (the OCCT_DEBUG
        // print is compiled out).
        return Curve2d::Line(the_new_l2d);
    } else {
        // L259-267: case of BSpline curve.
        let a_bc: BSplineCurve2 = match a_c2d {
            Curve2d::BSpline(b) => b.clone(),
            _ => {
                // L264-267: Untreated curve type (the OCCT_DEBUG print is
                // compiled out).
                return a_c2d.clone();
            }
        };
        // L269: newC = down_cast<Geom2d_BSplineCurve>(aBC->Copy()).
        let mut new_c = a_bc.clone();
        // L270-272: StartPoint/EndPoint (Value(First/LastParameter())) and
        // theVector = gp_Vec2d(FirstPoint, LastPoint).
        let first_point = a_bc.point_at(bspline2d_first_parameter(&a_bc));
        let last_point = a_bc.point_at(bspline2d_last_parameter(&a_bc));
        let the_vector = last_point - first_point;
        // L273-275: p00(uf, vf), p01(uf, vl), p10(ul, vf) and the iso
        // vectors.
        let p00 = DVec2::new(uf, vf);
        let p01 = DVec2::new(uf, vl);
        let p10 = DVec2::new(ul, vf);
        let vect_iso_uf = p01 - p00;
        let vect_iso_vf = p10 - p00;

        // L277-290: gp_Trsf2d T — the U-parallel arm (a translation by the
        // period between the two U isos).
        if gp_vec2d_is_parallel(the_vector, vect_iso_uf, a_tol) {
            if (first_point.x - uf).abs() < (first_point.x - ul).abs() {
                // L282: T.SetTranslation(p00, p10).
                new_c = translate_bspline2d(&new_c, p10 - p00);
            } else {
                // L286: T.SetTranslation(p10, p00).
                new_c = translate_bspline2d(&new_c, p00 - p10);
            }
            // L288-289.
            return Curve2d::BSpline(new_c);
        }
        // (the commented-out L291-308 OCCT blocks stay untranslated)
        // L309-321: the V-parallel arm.
        else if gp_vec2d_is_parallel(the_vector, vect_iso_vf, a_tol) {
            if (first_point.y - vf).abs() < (first_point.y - vl).abs() {
                // L313: T.SetTranslation(p00, p01).
                new_c = translate_bspline2d(&new_c, p01 - p00);
            } else {
                // L317: T.SetTranslation(p01, p00).
                new_c = translate_bspline2d(&new_c, p00 - p01);
            }
            // L319-320.
            return Curve2d::BSpline(new_c);
        }
    }
    // L323-324: les courbes ne sont pas sur la couture.
    a_c2d.clone()
}

/// GAP: OCCT Geom2d_BezierCurve::Segment(U1, U2) (TKGeomBase/Geom2d_Bezier,
/// the work-around call at ShapeFix_Edge.cxx L411/L440) — segments a Bezier
/// pcurve copy to [U1, U2].  Pending the TKGeomBase Bezier batch; the GAP
/// returns the unsegmented copy (documented deviation source: the copy keeps
/// the original parametrization).
fn geom2d_bezier_segment_gap(bezier: &Curve2d, _u1: f64, _u2: f64) -> Curve2d {
    bezier.clone()
}

/// GAP: OCCT GeomLib::SameRange(Tolerance, CurvePtr, FirstOnCurve,
/// LastOnCurve, RequestedFirst, RequestedLast, NewCurvePtr)
/// (TKGeomBase/GeomLib/GeomLib.cxx L842-967) — reparametrizes/segments the
/// pcurve onto the requested range (line translation, conic rotation, or the
/// Geom2dConvert::CurveToBSplineCurve + BSplCLib::Reparametrize walk).
/// Pending the TKGeomBase GeomLib/Geom2dConvert batch (the geomalgo
/// `geom_lib_same_range.rs` GAP carries the same anchor); the GAP returns
/// the input curve unchanged (documented deviation source: no
/// reparametrization).
fn geom_lib_same_range_gap(
    _tolerance: f64,
    curve2d: &Curve2d,
    _first_on_curve: f64,
    _last_on_curve: f64,
    _requested_first: f64,
    _requested_last: f64,
) -> Curve2d {
    curve2d.clone()
}

/// OCCT static TempSameRange (ShapeFix_Edge.cxx L329-464) — a copy of
/// BRepLib::SameRange() modified to be able to fix seam edges (see the
/// :b0 comment, L329-333).
fn temp_same_range(brep: &mut BRep, an_edge: &Shape, tolerance: f64) {
    // L345-347: first_time_in and the target range accumulators.
    let mut first_time_in = true;
    let (mut current_first, mut current_last) = (0.0f64, 0.0f64);

    // L349-353: C = BRep_Tool::Curve(AnEdge, LocalLoc, current_first,
    // current_last); when present, first_time_in becomes false.
    let has_c3d = match brep.tshapes[an_edge.index].as_ref() {
        TShape::Edge(ed) => {
            current_first = ed.range[0];
            current_last = ed.range[1];
            ed.curve.is_some()
        }
        _ => false,
    };
    if has_c3d {
        first_time_in = false;
    }

    // L355-460: the representation-list walk (each GCurve row in turn).
    // The rows are snapshotted (the OCCT in-loop writes touch only pcurve
    // values, never the ranges read below — the sequential read flow is
    // preserved).
    let rows: Vec<CurveRepresentation> = match brep.tshapes[an_edge.index].as_ref() {
        TShape::Edge(ed) => ed.representations.clone(),
        _ => Vec::new(),
    };
    // The collected per-row write-backs (L425/L454 the
    // geometric_representation_ptr->PCurve()/PCurve2(NewCurve...) writes).
    let mut writes: Vec<((u64, u32), Option<Curve2d>, Option<Curve2d>)> = Vec::new();
    for row in &rows {
        // L360-362: first / last of the geometric representation row and the
        // pcurve slots (Curve2dPtr / Curve2dPtr2).
        let (first, last);
        let mut curve2d_ptr: Option<Curve2d> = None;
        let mut curve2d_ptr2: Option<Curve2d> = None;
        let (mut has_curve, mut has_closed_curve) = (false, false);
        let face_key: (u64, u32);
        match row {
            CurveRepresentation::CurveOnSurface { face, pcurve, range } => {
                // L363-366: IsCurveOnSurface -> PCurve().
                face_key = *face;
                curve2d_ptr = Some(pcurve.clone());
                (first, last) = (range[0], range[1]);
                has_curve = true;
            }
            CurveRepresentation::CurveOnClosedSurface { face, pcurve1, pcurve2, range } => {
                // L368-372: IsCurveOnClosedSurface -> PCurve()/PCurve2().
                face_key = *face;
                curve2d_ptr = Some(pcurve1.clone());
                curve2d_ptr2 = Some(pcurve2.clone());
                (first, last) = (range[0], range[1]);
                has_curve = true;
                has_closed_curve = true;
            }
            _ => {
                // The BRep_Curve3D / BRep_CurveOn2Surfaces / regularity
                // rows: no pcurve slots (has_curve = has_closed_curve =
                // false, L360).
                continue;
            }
        }

        if has_curve || has_closed_curve {
            // L375-380: the first row initializes the target range.
            if first_time_in {
                current_first = first;
                current_last = last;
                first_time_in = false;
            }

            // L382-385: the range-mismatch test (PConfusion, :b8).
            if (first - current_first).abs() > PCONFUSION
                || (last - current_last).abs() > PCONFUSION
            {
                // L386: oldFirst = 0., oldLast = 0. (skl).
                let mut old_first = 0.0f64;
                let mut old_last = 0.0f64;
                if has_curve {
                    // L389-391: pdn 20.05.99 work around — the row range.
                    old_first = first;
                    old_last = last;
                    // L392-400: 15.11.2002 PTV OCC966 — the periodic pcurve
                    // trim wrap.
                    let mut c2d = curve2d_ptr.clone().expect("has_curve row");
                    if ShapeAnalysisCurve.is_periodic_2d(&c2d) {
                        let tc = Curve2d::Trimmed(TrimmedCurve2 {
                            curve: Box::new(c2d.clone()),
                            t_min: old_first,
                            t_max: old_last,
                        });
                        // L397: shift = tc->FirstParameter() - oldFirst.
                        // Arch. diff.: the rcad TrimmedCurve2 stores the
                        // requested trim (no OCCT periodic recurl), so the
                        // shift reduces to zero — the OCCT FirstParameter
                        // wrap adjustment is a documented deviation source.
                        let shift = match &tc {
                            Curve2d::Trimmed(t) => t.t_min - old_first,
                            _ => 0.0,
                        };
                        old_first += shift;
                        old_last += shift;
                    }
                    // L401-416: pdn 30.06.2000 work around on beziers.
                    let (mut old_first_curve1, mut old_last_curve1) = (old_first, old_last);
                    if matches!(c2d, Curve2d::Bezier(_)) {
                        if old_first.abs() > PCONFUSION || (old_last - 1.0).abs() > PCONFUSION {
                            // L409-413: bezier = Copy(); Segment(oldFirst,
                            // oldLast); Curve2dPtr = bezier.
                            c2d = geom2d_bezier_segment_gap(&c2d, old_first, old_last);
                        }
                        // L414-415.
                        old_first_curve1 = 0.0;
                        old_last_curve1 = 1.0;
                    }

                    // L418-424: GeomLib::SameRange(Tolerance, Curve2dPtr,
                    // oldFirstCurve1, oldLastCurve1, current_first,
                    // current_last, NewCurve2dPtr).
                    let new_curve2d_ptr = geom_lib_same_range_gap(
                        tolerance,
                        &c2d,
                        old_first_curve1,
                        old_last_curve1,
                        current_first,
                        current_last,
                    );
                    // L425: geometric_representation_ptr->PCurve(NewCurve2dPtr).
                    curve2d_ptr = Some(new_curve2d_ptr);
                }
                if has_closed_curve {
                    // L430: oldFirstCurve2 = oldFirst, oldLastCurve2 =
                    // oldLast (the values possibly shifted above).
                    let mut c2d2 = curve2d_ptr2.clone().expect("has_closed_curve row");
                    let (mut old_first_curve2, mut old_last_curve2) = (old_first, old_last);
                    // L432-445: the Bezier work-around on the second pcurve.
                    if matches!(c2d2, Curve2d::Bezier(_)) {
                        if old_first.abs() > PCONFUSION || (old_last - 1.0).abs() > PCONFUSION {
                            // L438-442: bezier = Copy(); Segment(oldFirst,
                            // oldLast); Curve2dPtr2 = bezier.
                            c2d2 = geom2d_bezier_segment_gap(&c2d2, old_first, old_last);
                        }
                        // L443-444.
                        old_first_curve2 = 0.0;
                        old_last_curve2 = 1.0;
                    }

                    // L447-453: GeomLib::SameRange(..., NewCurve2dPtr2).
                    let new_curve2d_ptr2 = geom_lib_same_range_gap(
                        tolerance,
                        &c2d2,
                        old_first_curve2,
                        old_last_curve2,
                        current_first,
                        current_last,
                    );
                    // L454: geometric_representation_ptr->PCurve2(NewCurve2dPtr2).
                    curve2d_ptr2 = Some(new_curve2d_ptr2);
                }
            }
        }

        if has_curve || has_closed_curve {
            writes.push((face_key, curve2d_ptr, curve2d_ptr2));
        }
    }

    // The write-back pass for the collected rows.
    for (face_key, new1, new2) in writes {
        let ed = brep.edge_mut_inplace(an_edge.clone());
        for r in ed.representations.iter_mut() {
            match r {
                CurveRepresentation::CurveOnSurface { face: fk, pcurve, .. } => {
                    if *fk == face_key {
                        if let Some(nc) = &new1 {
                            *pcurve = nc.clone();
                        }
                    }
                }
                CurveRepresentation::CurveOnClosedSurface {
                    face: fk,
                    pcurve1,
                    pcurve2,
                    ..
                } => {
                    if *fk == face_key {
                        if let Some(nc) = &new1 {
                            *pcurve1 = nc.clone();
                        }
                        if let Some(nc) = &new2 {
                            *pcurve2 = nc.clone();
                        }
                    }
                }
                _ => {}
            }
        }
        if let Some(entry) = ed.pcurves.get_mut(&face_key) {
            if let Some(nc) = &new1 {
                entry.0 = nc.clone();
            }
        }
    }

    // L461-463: B.Range(TopoDS::Edge(AnEdge), current_first, current_last)
    // (the 3-arg form — all representations) and B.SameRange(AnEdge, true).
    builder_range_all(brep, an_edge, current_first, current_last);
    let mut builder = BRepBuilder::new();
    builder.set_edge_same_range(brep, an_edge.clone(), true);
}

// ---------------------------------------------------------------------------
// The class.
// ---------------------------------------------------------------------------

/// OCCT ShapeFix_Edge (ShapeFix_Edge.hxx L46-247): fixing invalid edge
/// (missing 3d curve or pcurve, mismatching orientations, incorrect
/// SameParameter flag, curves not adjacent to the vertices).
///
/// OCCT inherits Standard_Transient directly (NOT ShapeFix_Root) and owns
/// its myContext / myStatus / myProjector (hxx L244-246).
pub struct ShapeFixEdge {
    /// OCCT myContext (hxx L244).
    my_context: Option<ShapeBuildReShape>,
    /// OCCT myStatus (hxx L245).
    my_status: i32,
    /// OCCT myProjector (hxx L246).
    my_projector: ProjectCurveOnSurface,
}

impl Default for ShapeFixEdge {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixEdge {
    /// OCCT ShapeFix_Edge::ShapeFix_Edge() (cxx L63-67).
    pub fn new() -> Self {
        ShapeFixEdge {
            my_context: None,
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_projector: ProjectCurveOnSurface::new(),
        }
    }

    /// OCCT ShapeFix_Edge::Projector (cxx L71-74): the projector used for
    /// recomputing missing pcurves (adjustable by the caller).
    pub fn projector(&mut self) -> &mut ProjectCurveOnSurface {
        &mut self.my_projector
    }

    /// OCCT ShapeFix_Edge::FixRemovePCurve(edge, face) (cxx L78-83): the
    /// face form — surface and location of the face, then the (surface,
    /// location) form.
    pub fn fix_remove_pcurve_face(&mut self, brep: &mut BRep, edge: &Shape, face: &Shape) -> bool {
        // L80-82: L; S = BRep_Tool::Surface(face, L); FixRemovePCurve(edge,
        // S, L).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.fix_remove_pcurve_surface(brep, edge, &s, l),
            None => false,
        }
    }

    /// OCCT ShapeFix_Edge::FixRemovePCurve(edge, surface, location)
    /// (cxx L87-99): removes the pcurve when it does not match the vertices.
    pub fn fix_remove_pcurve_surface(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
    ) -> bool {
        // L91.
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // L92-93.
        let mut ea = ShapeAnalysisEdge::new();
        let result = ea.check_vertices_with_pcurve_surface(brep, edge, surface, location, 0.0, 0);
        if result {
            // L96: ShapeBuild_Edge().RemovePCurve(edge, surface, location).
            let sbe = ShapeBuildEdge;
            sbe.remove_pcurve_surface_loc(brep, edge, surface, location);
        }
        // L98.
        result
    }

    /// OCCT ShapeFix_Edge::FixRemoveCurve3d (cxx L103-113): removes the 3d
    /// curve when it does not match the vertices.
    pub fn fix_remove_curve3d(&mut self, brep: &mut BRep, edge: &Shape) -> bool {
        // L105.
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // L106-107.
        let mut ea = ShapeAnalysisEdge::new();
        let result = ea.check_vertices_with_curve3d(brep, edge, 0.0, 0);
        if result {
            // L110: ShapeBuild_Edge().RemoveCurve3d(edge).
            let sbe = ShapeBuildEdge;
            sbe.remove_curve_3d(brep, edge);
        }
        // L112.
        result
    }

    /// OCCT ShapeFix_Edge::FixAddPCurve(edge, face, isSeam, prec = 0.0)
    /// (cxx L117-125): the face form.
    pub fn fix_add_pcurve_face(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        face: &Shape,
        is_seam: bool,
        prec: f64,
    ) -> bool {
        // L122-124: L; S = BRep_Tool::Surface(face, L); FixAddPCurve(edge,
        // S, L, isSeam, prec).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.fix_add_pcurve_surface(brep, edge, &s, l, is_seam, prec),
            None => false,
        }
    }

    /// OCCT ShapeFix_Edge::FixAddPCurve(edge, surface, location, isSeam,
    /// prec = 0.0) (cxx L129-143): transforms the surface by the location
    /// and builds the analyzer, then the analyzer form.
    pub fn fix_add_pcurve_surface(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
        is_seam: bool,
        prec: f64,
    ) -> bool {
        // L135-140: aTransSurf = surface->Transformed(aTrsf) when the
        // location is not identity.
        let a_trans_surf = if location != 0 {
            let trsf = brep.get_location(location);
            transform_surface(surface, &trsf)
        } else {
            surface.clone()
        };
        // L141: sas = new ShapeAnalysis_Surface(aTransSurf).
        let mut sas = ShapeAnalysisSurface::new(a_trans_surf);
        // L142.
        self.fix_add_pcurve_surface_sa(brep, edge, surface, location, is_seam, &mut sas, prec)
    }

    /// OCCT ShapeFix_Edge::FixAddPCurve(edge, face, isSeam, surfana,
    /// prec = 0.0) (cxx L147-156): the face + analyzer form.
    pub fn fix_add_pcurve_face_sa(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        face: &Shape,
        is_seam: bool,
        surfana: &mut ShapeAnalysisSurface,
        prec: f64,
    ) -> bool {
        // L153-155: L; S = BRep_Tool::Surface(face, L); FixAddPCurve(edge,
        // S, L, isSeam, surfana, prec).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.fix_add_pcurve_surface_sa(brep, edge, &s, l, is_seam, surfana, prec),
            None => false,
        }
    }

    /// OCCT ShapeFix_Edge::FixAddPCurve(edge, surf, location, isSeam, sas,
    /// prec = 0.0) (cxx L470-614): adds the missing pcurve(s) by projecting
    /// the 3d curve.
    pub fn fix_add_pcurve_surface_sa(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        surf: &Surface3,
        location: u32,
        is_seam: bool,
        sas: &mut ShapeAnalysisSurface,
        prec: f64,
    ) -> bool {
        // L477-478.
        let sae = ShapeAnalysisEdge::new();
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // L479-483.
        if (!is_seam && sae.has_pcurve_surface(brep, edge, surf, location))
            || (is_seam && sae.is_seam_surface(brep, edge, surf, location))
        {
            return false;
        }

        // L485-489: PCurve on Plane not computed.
        if matches!(surf, Surface3::Plane(_)) {
            return false;
        }

        // L492-494: try { OCC_CATCH_SIGNALS — the Standard_Failure arm
        // (L601-611) sets FAIL2; rcad raises no exceptions in this body
        // (bridge #5).
        {
            // L499: preci = (prec > 0. ? prec : BRep_Tool::Tolerance(edge)).
            let preci = if prec > 0.0 { prec } else { brep.tolerance(edge) };

            // L500-506: c3d = BRep_Tool::Curve(edge, First, Last).
            let (c3d, first, last) = match edge.data.as_ref() {
                TShape::Edge(ed) => match &ed.curve {
                    Some(c) => (c.clone(), ed.range[0], ed.range[1]),
                    None => {
                        self.my_status |= encode_status(ShapeExtendStatus::Fail1);
                        return false;
                    }
                },
                _ => {
                    self.my_status |= encode_status(ShapeExtendStatus::Fail1);
                    return false;
                }
            };

            // (the commented-out L508-511 OCCT trim blocks stay untranslated)
            // L517-546: the projection or the stored pcurve.
            let mut c2d_opt: Option<Curve2d> = None;
            let a1: f64;
            let b1: f64;
            if !sae.has_pcurve_surface(brep, edge, surf, location) {
                // L521-531: TolFirst/TolLast from the extremity vertices.
                let mut tol_first = -1.0;
                let mut tol_last = -1.0;
                let (v1, v2) = top_exp_vertices(brep, edge);
                if !v1.is_null() {
                    tol_first = brep.tolerance(&v1);
                }
                if !v2.is_null() {
                    tol_last = brep.tolerance(&v2);
                }

                // L533-534: myProjector->Init(sas, preci) + Perform.  Arch.
                // diff.: the rcad projector owns its analyzer (no handle
                // sharing); the Init receives a fresh analyzer over the same
                // surface — the caller analyzer keeps serving the later
                // IsUClosed/IsVClosed/Surface queries (L564-579); both views
                // derive all state from the same surface.  (The projector is
                // still wired to the W1-4 analyzer carrier; the W2 real
                // analyzer parameter is narrowed to the same surface.)
                self.my_projector
                    .init_sa(ProjectorShapeAnalysisSurface::new(sas.surface().clone()), preci);
                let mut c2d_out: Option<Curve2d> = None;
                self.my_projector
                    .perform(&c3d, first, last, &mut c2d_out, tol_first, tol_last);
                c2d_opt = c2d_out;
                // L536-539: Status(DONE4) -> DONE2.
                if self.my_projector.status(ShapeExtendStatus::Done4) {
                    self.my_status |= encode_status(ShapeExtendStatus::Done2);
                }
                // L540-541.
                a1 = first;
                b1 = last;
            } else {
                // L545: sae.PCurve(edge, surf, location, c2d, a1, b1, false).
                let mut c2d_out: Option<Curve2d> = None;
                let mut fa1 = 0.0;
                let mut lb1 = 0.0;
                sae.pcurve_surface(
                    brep,
                    edge,
                    surf,
                    location,
                    &mut c2d_out,
                    &mut fa1,
                    &mut lb1,
                    false,
                );
                c2d_opt = c2d_out;
                a1 = fa1;
                b1 = lb1;
            }
            let c2d = match c2d_opt {
                Some(c) => c,
                None => {
                    // OCCT dereferences the handle at L554/L581; the null
                    // case cannot occur on the OCCT success path — kept as
                    // the FAIL2 guard (bridge #5).
                    self.my_status |= encode_status(ShapeExtendStatus::Fail2);
                    return false;
                }
            };

            // L550-585.
            if is_seam {
                // L553-554: c2d2 = down_cast<Geom2d_Curve>(c2d->Copy()).
                let mut c2d2 = c2d.clone();
                // L558-559: surf->Bounds(uf, ul, vf, vl).
                let mut uf = 0.0;
                let mut ul = 0.0;
                let mut vf = 0.0;
                let mut vl = 0.0;
                ShapeAnalysisSurface::new(surf.clone()).bounds(&mut uf, &mut ul, &mut vf, &mut vl);
                // L560-567: #4/#13/#78 rln — the closed-side selection (the
                // commented-out spherical-surface arm stays untranslated).
                if sas.is_u_closed(prec) && !sas.is_v_closed(prec) {
                    // L568-569: tranvec(ul - uf, 0).
                    c2d2 = translate_curve2d(&c2d2, DVec2::new(ul - uf, 0.0));
                }
                // L571-575.
                else if sas.is_v_closed(prec) && !sas.is_u_closed(prec) {
                    // L573-574: tranvec(0, vl - vf).
                    c2d2 = translate_curve2d(&c2d2, DVec2::new(0.0, vl - vf));
                }
                // L576-580: q8 abv — the doubly-closed case; the OCCT
                // IsUClosed()/IsVClosed() no-arg defaults are
                // Precision::Confusion().
                else if sas.is_u_closed(CONFUSION) && sas.is_v_closed(CONFUSION) {
                    // L579: c2d2 = TranslatePCurve(sas->Surface(), c2d2, prec).
                    c2d2 = translate_pcurve(sas.surface(), &c2d2, prec);
                }
                // L581-584: B.UpdateEdge(edge, c2d, c2d2, surf, location,
                // 0.); B.Range(edge, surf, location, a1, b1).  (The
                // commented-out L582-583 OCCT block stays untranslated.)
                builder_update_edge_pcurves_surface(brep, edge, &c2d, &c2d2, surf, location, 0.0);
                builder_range_on_surface_loc(brep, edge, surf, location, a1, b1);
            } else {
                // L588: B.UpdateEdge(edge, c2d, surf, location, 0.).
                builder_update_edge_pcurve_surface(brep, edge, &c2d, surf, location, 0.0);
            }

            // L592-599: the DONE3 conclusion — the 3d curve range enforced.
            if self.my_projector.status(ShapeExtendStatus::Done3) {
                // L595-596.
                let (g3d_c_first, g3d_c_last) = {
                    let dom = c3d.default_domain();
                    (dom[0], dom[1])
                };
                // L597: B.UpdateEdge(edge, c3d, 0.) — c3d is already the
                // stored 3d curve (BRep_Tool::Curve fetched it, L500); the
                // value write is a no-op on the rcad TShape (bridge #3).
                // L598: B.Range(edge, G3dCFirst, G3dCLast, true).
                builder_range_only3d(brep, edge, g3d_c_first, g3d_c_last);
            }
        } // end try
          // L612-613: myStatus |= DONE1; return true (the OCCT fall-through
          // sets DONE1 on the success and the catch path alike).
        self.my_status |= encode_status(ShapeExtendStatus::Done1);
        true
    }

    /// OCCT ShapeFix_Edge::FixAddCurve3d (cxx L618-638): builds the missing
    /// 3d curve.
    pub fn fix_add_curve3d(&mut self, brep: &mut BRep, edge: &Shape) -> bool {
        // L620.
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // L621-625.
        let ea = ShapeAnalysisEdge::new();
        if brep_tool_degenerated(edge) || ea.has_curve3d(brep, edge) {
            return false;
        }
        // L626-629.
        let same_range = matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.same_range);
        if !same_range {
            temp_same_range(brep, edge, PCONFUSION);
        }

        // L631-635.
        let sbe = ShapeBuildEdge;
        if !sbe.build_curve3d(brep, edge) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        // L636-637.
        self.my_status |= encode_status(ShapeExtendStatus::Done1);
        true
    }

    /// OCCT ShapeFix_Edge::FixVertexTolerance(edge, face) (cxx L642-685):
    /// raises the vertex tolerances to comprise the ends of the 3d curve and
    /// of the pcurve on the given face.
    pub fn fix_vertex_tolerance_face(&mut self, brep: &mut BRep, edge: &Shape, face: &Shape) -> bool {
        // L644.
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // L645-656: anEdgeCopy = Context()->Apply(edge) when the context is
        // set.
        let mut an_edge_copy = edge.clone();
        if let Some(ctx) = self.my_context.as_mut() {
            let a_shape = ctx.apply(brep, edge, ShapeType::Shape);
            if a_shape.is_null() || a_shape.shape_type() != ShapeType::Edge {
                return false;
            }
            an_edge_copy = a_shape;
        }

        // L658-662.
        let mut sae = ShapeAnalysisEdge::new();
        let mut toler1 = 0.0;
        let mut toler2 = 0.0;
        if !sae.check_vertex_tolerance_face(brep, &an_edge_copy, face, &mut toler1, &mut toler2) {
            return false;
        }
        // L663-666.
        if sae.status(ShapeExtendStatus::Done1) {
            self.my_status = encode_status(ShapeExtendStatus::Done1);
        }
        // L667-670.
        if sae.status(ShapeExtendStatus::Done2) {
            self.my_status = encode_status(ShapeExtendStatus::Done2);
        }
        // L671-683.
        let v1 = sae.first_vertex(brep, &an_edge_copy);
        let v2 = sae.last_vertex(brep, &an_edge_copy);
        if let Some(ctx) = self.my_context.as_mut() {
            ctx.copy_vertex(brep, &v1, toler1);
            ctx.copy_vertex(brep, &v2, toler2);
        } else {
            let mut b = BRepBuilder::new();
            b.update_vertex_tolerance(brep, v1, toler1);
            b.update_vertex_tolerance(brep, v2, toler2);
        }
        // L684.
        true
    }

    /// OCCT ShapeFix_Edge::FixVertexTolerance(edge) (cxx L689-731): the
    /// all-stored-pcurves form.
    pub fn fix_vertex_tolerance(&mut self, brep: &mut BRep, edge: &Shape) -> bool {
        // L691.
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // L692-703.
        let mut an_edge_copy = edge.clone();
        if let Some(ctx) = self.my_context.as_mut() {
            let a_shape = ctx.apply(brep, edge, ShapeType::Shape);
            if a_shape.is_null() || a_shape.shape_type() != ShapeType::Edge {
                return false;
            }
            an_edge_copy = a_shape;
        }
        // L704-708.
        let mut sae = ShapeAnalysisEdge::new();
        let mut toler1 = 0.0;
        let mut toler2 = 0.0;
        if !sae.check_vertex_tolerance_all(brep, &an_edge_copy, &mut toler1, &mut toler2) {
            return false;
        }
        // L709-712.
        if sae.status(ShapeExtendStatus::Done1) {
            self.my_status = encode_status(ShapeExtendStatus::Done1);
        }
        // L713-716.
        if sae.status(ShapeExtendStatus::Done2) {
            self.my_status = encode_status(ShapeExtendStatus::Done2);
        }
        // L717-729.
        let v1 = sae.first_vertex(brep, &an_edge_copy);
        let v2 = sae.last_vertex(brep, &an_edge_copy);
        if let Some(ctx) = self.my_context.as_mut() {
            ctx.copy_vertex(brep, &v1, toler1);
            ctx.copy_vertex(brep, &v2, toler2);
        } else {
            let mut b = BRepBuilder::new();
            b.update_vertex_tolerance(brep, v1, toler1);
            b.update_vertex_tolerance(brep, v2, toler2);
        }
        // L730.
        true
    }

    /// OCCT ShapeFix_Edge::FixReversed2d(edge, face) (cxx L735-740): the
    /// face form.
    pub fn fix_reversed2d_face(&mut self, brep: &mut BRep, edge: &Shape, face: &Shape) -> bool {
        // L738-739: L; S = BRep_Tool::Surface(face, L); FixReversed2d(edge,
        // S, L).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.fix_reversed2d_surface(brep, edge, &s, l),
            None => false,
        }
    }

    /// OCCT ShapeFix_Edge::FixReversed2d(edge, surface, location)
    /// (cxx L744-786): reverses the pcurve directed opposite to the 3d
    /// curve.
    pub fn fix_reversed2d_surface(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
    ) -> bool {
        // L748.
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        // L750-763.
        let mut ea = ShapeAnalysisEdge::new();
        ea.check_curve3d_with_pcurve_surface(brep, edge, surface, location);
        if ea.status(ShapeExtendStatus::Fail1) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        if ea.status(ShapeExtendStatus::Fail2) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail2);
        }
        if !ea.status(ShapeExtendStatus::Done) {
            return false;
        }

        // L765-767: EA.PCurve(edge, surface, location, c2d, f, l, false).
        let mut c2d_opt: Option<Curve2d> = None;
        let mut f = 0.0;
        let mut l = 0.0;
        ea.pcurve_surface(brep, edge, surface, location, &mut c2d_opt, &mut f, &mut l, false);
        let c2d = match c2d_opt {
            Some(c) => c,
            None => return false,
        };
        // L769 (the #46 rln comment stays untranslated): newf/newl =
        // c2d->ReversedParameter(l/f).
        let newf = c2d.reversed_parameter(l);
        let newl = c2d.reversed_parameter(f);
        // L770: c2d->Reverse() — the shared-handle mutation (bridge #3).
        let c2d = reverse_curve2d(&c2d);
        set_stored_pcurve(brep, edge, surface, location, c2d);
        // L771-773: B.Range(edge, surface, location, newf, newl) — will
        // break seams (the L772 comment stays untranslated).
        builder_range_on_surface_loc(brep, edge, surface, location, newf, newl);
        // L774-783: #51 rln — the SameRange/SameParameter reset when the
        // resulting edge range disagrees.
        let (first, last) = brep_tool_range(brep, edge).unwrap_or((0.0, 0.0));
        if first != newf || last != newl {
            let mut b = BRepBuilder::new();
            b.set_edge_same_range(brep, edge.clone(), false);
            b.set_edge_same_parameter(brep, edge.clone(), false);
        }
        // L784-785.
        self.my_status |= encode_status(ShapeExtendStatus::Done1);
        true
    }

    /// OCCT ShapeFix_Edge::FixSameParameter(edge, tolerance = 0.0)
    /// (cxx L790-794): the empty-face form.
    pub fn fix_same_parameter(&mut self, brep: &mut BRep, edge: &Shape, tolerance: f64) -> bool {
        // L792-793: anEmptyFace; FixSameParameter(edge, anEmptyFace,
        // tolerance).
        self.fix_same_parameter_face(brep, edge, &Shape::null(), tolerance)
    }

    /// OCCT ShapeFix_Edge::FixSameParameter(edge, face, tolerance = 0.0)
    /// (cxx L798-936): makes the edge SameParameter and sets the
    /// corresponding tolerance and flag (the header L160-230 documents the
    /// status contract).
    pub fn fix_same_parameter_face(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        face: &Shape,
        tolerance: f64,
    ) -> bool {
        // L802.
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        // L804-813: the degenerated edge.
        if brep_tool_degenerated(edge) {
            let same_range = matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.same_range);
            if !same_range {
                // L809: TempSameRange(edge, Precision::PConfusion()).
                temp_same_range(brep, edge, PCONFUSION);
            }
            // L811: B.SameParameter(edge, true).
            let mut b = BRepBuilder::new();
            b.set_edge_same_parameter(brep, edge.clone(), true);
            return false;
        }

        // L815-817.
        let sfst = ShapeFixShapeTolerance::new();

        // L819-824: the extremity vertices and the current tolerances.
        let sae0 = ShapeAnalysisEdge::new();
        let v1 = sae0.first_vertex(brep, edge);
        let v2 = sae0.last_vertex(brep, edge);
        let tol_fv = if v1.is_null() { 0.0 } else { brep.tolerance(&v1) };
        let tol_lv = if v2.is_null() { 0.0 } else { brep.tolerance(&v2) };
        let tol = brep.tolerance(edge);

        // L826: wasSP = BRep_Tool::SameParameter(edge), SP = false.
        let was_sp = matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.same_parameter);
        let mut sp = false;
        // L819: copyedge — declared here (L819) and assigned inside the try
        // block, consumed after it.
        let mut copyedge = Shape::null();
        // L827-868: the try/catch (bridge #5 — the catch arm at L858-867
        // sets FAIL2).
        {
            // L831-834.
            let same_range = matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.same_range);
            if !same_range {
                temp_same_range(brep, edge, PCONFUSION);
            }
            // L835-856: #81 rln — choose the best result for a not-SP edge
            // (the L843-845 comment stays untranslated).
            if !was_sp {
                // L841: copyedge = ShapeBuild_Edge().Copy(edge, false).
                let sbe = ShapeBuildEdge;
                copyedge = sbe.copy(brep, edge, false);
                let mut b = BRepBuilder::new();
                // L842.
                b.set_edge_same_parameter(brep, copyedge.clone(), false);
                // L847-849: enforce the original 3D range (Only3d = true).
                let (a_f, a_l) = brep_tool_range(brep, edge).unwrap_or((0.0, 0.0));
                builder_range_only3d(brep, &copyedge, a_f, a_l);
                // L850: BRepLib::SameParameter(copyedge, (tolerance >=
                // Precision::Confusion() ? tolerance : tol)) — the docket
                // section 4 gap 3 kernel GAP (bridge #6; the OCCT 8 2-arg
                // form, BRepLib.hxx L161).
                let brl_tol = if tolerance >= CONFUSION { tolerance } else { tol };
                brep_lib_same_parameter_edge(brep, &copyedge, brl_tol, false);
                // L851-855.
                sp = matches!(copyedge.data.as_ref(), TShape::Edge(ed) if ed.same_parameter);
                if !sp {
                    self.my_status |= encode_status(ShapeExtendStatus::Fail2);
                }
            }
        }

        // L871-872: compute the deviation on the original pcurves.
        let mut maxdev = 0.0f64;
        {
            let mut b = BRepBuilder::new();
            b.set_edge_same_parameter(brep, edge.clone(), true);
        }

        // L874-880: check all pcurves when the input was not SP.
        let a_face = if !was_sp { Shape::null() } else { face.clone() };

        // L882: sae.CheckSameParameter(edge, aFace, maxdev) (the NbControl
        // default is 23, ShapeAnalysis_Edge.hxx).
        let mut sae = ShapeAnalysisEdge::new();
        sae.check_same_parameter_face(brep, edge, &a_face, &mut maxdev, 23);
        // L883-886.
        if sae.status(ShapeExtendStatus::Fail2) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
        }

        // L889-912: if BRepLib was OK, compare and select the best variant.
        if sp {
            // L891: BRLTol = BRep_Tool::Tolerance(copyedge), BRLDev.
            let mut brl_tol = brep.tolerance(&copyedge);
            let mut brl_dev = 0.0f64;
            // L892: sae.CheckSameParameter(copyedge, BRLDev).
            sae.check_same_parameter(brep, &copyedge, &mut brl_dev, 23);
            // L893.
            self.my_status |= encode_status(ShapeExtendStatus::Done3);
            // L894-897.
            if brl_tol < brl_dev {
                brl_tol = brl_dev;
            }

            // L899-911: chose the best result.
            if brl_tol < maxdev {
                if sae.status(ShapeExtendStatus::Fail2) {
                    // L902-905.
                    self.my_status |= encode_status(ShapeExtendStatus::Fail1);
                }
                // L907: copy pcurves and tolerances from copyedge.
                let sbe = ShapeBuildEdge;
                sbe.copy_pcurves(brep, edge, &copyedge);
                // L908.
                maxdev = brl_tol;
                // L909: SFST.SetTolerance(edge, BRLTol, TopAbs_EDGE).
                sfst.set_tolerance(brep, edge, brl_tol, ShapeType::Edge);
                // L910.
                self.my_status |= encode_status(ShapeExtendStatus::Done5);
            }
        }

        // L914-922: restore the vertex tolerances (they could be modified by
        // BRepLib).
        if !v1.is_null() {
            sfst.set_tolerance(brep, &v1, maxdev.max(tol_fv), ShapeType::Vertex);
        }
        if !v2.is_null() {
            sfst.set_tolerance(brep, &v2, maxdev.max(tol_lv), ShapeType::Vertex);
        }

        // L924-929.
        if maxdev > tol {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
            // L927: B.UpdateEdge(edge, maxdev).
            let mut b = BRepBuilder::new();
            b.update_edge_tolerance(brep, edge.clone(), maxdev);
            // L928.
            self.fix_vertex_tolerance(brep, edge);
        }

        // L931-934.
        if !was_sp && !sp {
            self.my_status |= encode_status(ShapeExtendStatus::Done2);
        }
        // L935.
        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT ShapeFix_Edge::Status (cxx L940-943).
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    /// OCCT ShapeFix_Edge::Context (cxx L947-950).
    pub fn context(&self) -> Option<&ShapeBuildReShape> {
        self.my_context.as_ref()
    }

    /// OCCT ShapeFix_Edge::SetContext (cxx L954-957).
    pub fn set_context(&mut self, context: ShapeBuildReShape) {
        self.my_context = Some(context);
    }
}
