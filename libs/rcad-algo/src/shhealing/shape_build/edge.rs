//! OCCT ShapeBuild_Edge (TKShHealing ShapeBuild package).
//!
//! 1:1 translation of `ShapeBuild_Edge.hxx` L17-156 and
//! `ShapeBuild_Edge.cxx` L57-881: low-level operators for building an edge
//! 3d curve, copying edge with replaced vertices etc.
//!
//! Architecture bridges (numbered, referenced by the methods below):
//! 1. `BRep` pool argument - OCCT reaches the TEdge curve representations
//!    through the `TopoDS_Edge` TShape handle; rcad TShape data lives in the
//!    [`BRep`] pool, so every method takes `brep: &mut BRep` as that stand-in
//!    (shape_build/brep_tool.rs and kernel BRepTool precedents).
//! 2. `TopLoc_Location` -> `u32` (the `BRep.locations` table index,
//!    0 = identity).
//! 3. `Geom_Curve` / `Geom2d_Curve` / `Geom_Surface` handles ->
//!    `Curve3` / `Curve2d` / `Surface3` values; OCCT's `->Copy()` deep copy
//!    is rcad's value `clone`.
//! 4. `gp_Pnt` / `gp_Pnt2d` / `gp_Vec2d` -> `DVec3` / `DVec2` / `DVec2`.
//! 5. OCCT's `TEdge->ChangeCurves()` list of `BRep_CurveRepresentation` ->
//!    rcad `TEdgeData.representations: Vec<CurveRepresentation>` plus the
//!    `pcurves` IndexMap keyed `(face ptr, composed pcurve location)` — the
//!    same data under the read key (`BRepTool::curve_on_surface`). Every
//!    range/pcurve write below updates both mirrors. The 3D `BRep_Curve3D`
//!    representation maps to `TEdgeData.curve` + `TEdgeData.range` (the
//!    representation row's `curve: usize` carries only the kernel-side
//!    registration index).
//! 6. OCCT overloads are disambiguated with `_face` / `_surface` / `_loc`
//!    suffixes; output parameters use `&mut` (the AGENTS 1:1 rule).
//! 7. GAP carriers (other untranslated packages, with the OCCT failure
//!    path kept): `BRepLib::BuildCurve3d` (TKTopAlgo; reduced substitute
//!    documented on [`ShapeBuildEdge::build_curve3d`]), the per-variant
//!    `Geom2d_*::Transform` / `TransformedParameter` overrides (TKGeomBase,
//!    documented on `geom2d_curve_transform` / `geom2d_transformed_parameter`),
//!    `Geom2dConvert::CurveToBSplineCurve` + the full
//!    `Geom2dConvert_ApproxCurve` form (TKGeomBase, documented on
//!    `geom2d_convert_curve_to_bspline`), and the `BRepBuilderAPI_MakeEdge`
//!    constructors (TKTopAlgo; the docket defers the standalone Make* port
//!    decision to W4 — re-hosted on `brep_builder_api_make_edge_3d` /
//!    `brep_builder_api_make_edge_pcurve`).

use rcad_kernel::base::geom2d_convert::approx_curve_to_bspline;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, SurfaceEval, TrimmedCurve2,
    BSplineCurve2,
};
use rcad_kernel::precision::{
    is_negative_infinite_value, is_positive_infinite_value, APPROXIMATION, CONFUSION, PCONFUSION,
};
use rcad_kernel::topo::topods::{
    pcurve_location_id, surface_same, BRep, BRepTool, CurveRepresentation, Orientation, Shape,
    TShape,
};

use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, shape_is_null};

/// OCCT `gp::Resolution()` = the smallest positive double
/// (BRepLib_MakeEdge2d.rs GP_RESOLUTION precedent).
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

// ---------------------------------------------------------------------------
// File statics of ShapeBuild_Edge.cxx
// ---------------------------------------------------------------------------

/// OCCT ShapeBuild_Edge.cxx L148-162 — static AdjustByPeriod: the multiple of
/// `Period` that brings `Val` within half a period of `ToVal`.
fn adjust_by_period(val: f64, to_val: f64, period: f64) -> f64 {
    let diff = val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * ((d / p + 0.5) as i32 as f64)
}

/// OCCT ShapeBuild_Edge.cxx L164-183 — static IsPeriodic for a 3D curve:
/// ask IsPeriodic on the basis curve (unwrap Geom_OffsetCurve /
/// Geom_TrimmedCurve layers first; 15.11.2002 PTV OCC966).
fn is_periodic_curve3d(the_curve: &Curve3) -> bool {
    let mut a_tmp_curve = the_curve.clone();
    loop {
        match &a_tmp_curve {
            Curve3::Offset(off) => a_tmp_curve = off.basis.as_ref().clone(),
            Curve3::Trimmed(tc) => a_tmp_curve = tc.basis_curve().clone(),
            _ => break,
        }
    }
    a_tmp_curve.is_periodic()
}

/// OCCT ShapeBuild_Edge.cxx L185-204 — static IsPeriodic for a pcurve:
/// ask IsPeriodic on the basis curve (unwrap Geom2d_OffsetCurve /
/// Geom2d_TrimmedCurve layers first; 15.11.2002 PTV OCC966).
fn is_periodic_curve2d(the_curve: &Curve2d) -> bool {
    let mut a_tmp_curve = the_curve.clone();
    loop {
        match &a_tmp_curve {
            Curve2d::Offset(off) => a_tmp_curve = off.basis.as_ref().clone(),
            Curve2d::Trimmed(tc) => a_tmp_curve = tc.curve.as_ref().clone(),
            _ => break,
        }
    }
    a_tmp_curve.is_periodic()
}

/// OCCT Geom_Curve::Period for the CopyRanges periodic check: the conic
/// period is 2*Pi (Geom_Conic.cxx), a periodic BSpline's period is
/// last - first parameter (Geom_BSplineCurve.cxx).
fn curve3_period(the_curve: &Curve3) -> f64 {
    let mut a_tmp_curve = the_curve.clone();
    loop {
        match &a_tmp_curve {
            Curve3::Offset(off) => a_tmp_curve = off.basis.as_ref().clone(),
            Curve3::Trimmed(tc) => a_tmp_curve = tc.basis_curve().clone(),
            _ => break,
        }
    }
    match &a_tmp_curve {
        Curve3::Circle(_) | Curve3::Ellipse(_) => std::f64::consts::TAU,
        Curve3::BSpline(bs) if bs.is_periodic => {
            let [f, l] = a_tmp_curve.default_domain();
            l - f
        }
        // Not reachable through is_periodic_curve3d; a stable fallback only.
        _ => 1.0,
    }
}

/// OCCT Geom_Curve::FirstParameter / LastParameter — the natural parameter
/// domain.
fn curve3_first_last(the_curve: &Curve3) -> (f64, f64) {
    let [f, l] = the_curve.default_domain();
    (f, l)
}

/// OCCT Geom2d_Curve::FirstParameter / LastParameter — the natural parameter
/// domain of the basis curve.
fn curve2d_first_last(the_curve: &Curve2d) -> (f64, f64) {
    let mut a_tmp_curve = the_curve.clone();
    loop {
        match &a_tmp_curve {
            Curve2d::Offset(off) => a_tmp_curve = off.basis.as_ref().clone(),
            Curve2d::Trimmed(tc) => a_tmp_curve = tc.curve.as_ref().clone(),
            _ => break,
        }
    }
    let [f, l] = a_tmp_curve.default_domain();
    (f, l)
}

/// OCCT ShapeBuild_Edge.cxx L507-528 — static CountPCurves: count the exact
/// number of pcurves STORED in the edge for the face. This makes difference
/// for faces based on plane surfaces where pcurves can be not stored but
/// returned by BRep_Tools::CurveOnSurface.
fn count_pcurves(brep: &BRep, edge: &Shape, face: &Shape) -> i32 {
    // BRep_Tool::Surface(face, L): the face surface and its location.
    let (s, face_surface_loc) = match brep.tshapes[face.index].as_ref() {
        TShape::Face(fd) => match &fd.surface {
            Some(s) => (s.clone(), fd.surface_location),
            None => return 0,
        },
        _ => return 0,
    };
    // TopLoc_Location l = L.Predivided(edge.Location()) — the composed
    // pcurve-location hash stand-in (bridge #5).
    let l = compose_loc_hash(brep, face_surface_loc, edge.location);

    let reps = edge_representations(edge);
    for gc in reps {
        // GC->IsCurveOnSurface(S, l): the representation must be a
        // curve-on-surface kind whose surface matches S and location matches
        // l. rcad stores the owning face key; the surface identity resolves
        // through that face (surface_same stand-in for the handle match).
        let (face_ptr, rep_loc, is_closed) = match gc {
            CurveRepresentation::CurveOnSurface { face, .. } => (face.0, face.1, false),
            CurveRepresentation::CurveOnClosedSurface { face, .. } => (face.0, face.1, true),
            _ => continue,
        };
        let face_shape = match brep.index_by_ptr(face_ptr) {
            Some(idx) => brep.shape_at(idx),
            None => continue,
        };
        let rep_surface = match brep.tshapes[face_shape.index].as_ref() {
            TShape::Face(fd) => match &fd.surface {
                Some(s) => s.clone(),
                None => continue,
            },
            _ => continue,
        };
        if surface_same(&rep_surface, &s) && rep_loc == l {
            return if is_closed { 2 } else { 1 };
        }
    }
    0
}

/// The edge's curve representations (empty for non-edges) — bridge #5.
fn edge_representations(edg: &Shape) -> &[CurveRepresentation] {
    match edg.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => &[],
    }
}

/// OCCT `TopLoc_Location l = L.Predivided(edge.Location())` — the composed
/// pcurve-location hash stand-in (bridge #2/#5): the kernel
/// `compose_pcurve_location` is exactly the `Predivided` value hash.
fn compose_loc_hash(brep: &BRep, face_loc: u32, edge_loc: u32) -> u32 {
    rcad_kernel::topo::topods::compose_pcurve_location(face_loc, edge_loc, &brep.locations)
}

/// OCCT handle equality `c2d == pcurve0` for rcad curve values: the
/// Debug-rendering stand-in for the Geom2d_Curve handle identity match.
fn curve2d_handle_same(a: &Curve2d, b: &Curve2d) -> bool {
    format!("{:?}", a) == format!("{:?}", b)
}

// ---------------------------------------------------------------------------
// ShapeBuild_Edge class
// ---------------------------------------------------------------------------

/// OCCT ShapeBuild_Edge (`ShapeBuild_Edge.hxx` L35-154).
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeBuildEdge;

impl ShapeBuildEdge {
    /// OCCT ShapeBuild_Edge::CopyReplaceVertices (ShapeBuild_Edge.cxx L59-143):
    /// copies the edge and replaces one or both its vertices by the given
    /// one(s). Vertex V1 replaces the FORWARD vertex, and V2 - REVERSED, as
    /// they are found by TopoDS_Iterator. If V1 or V2 is null, the original
    /// vertex is taken.
    pub fn copy_replace_vertices(
        &self,
        brep: &mut BRep,
        edge: &Shape,
        v1: &Shape,
        v2: &Shape,
    ) -> Shape {
        // NCollection_Sequence<TopoDS_Shape> aNMVertices.
        let mut a_nm_vertices: Vec<Shape> = Vec::new();
        let mut new_v1 = v1.clone();
        let mut new_v2 = v2.clone();
        if shape_is_null(&new_v1) || shape_is_null(&new_v2) {
            // TopoDS_Iterator it: cumOri on for FORWARD/REVERSED edges, off
            // for INTERNAL/EXTERNAL; cumLoc always on.
            let cum_ori = edge.orientation == Orientation::Forward
                || edge.orientation == Orientation::Reversed;
            for v in iter_subshapes(brep, edge, cum_ori, true) {
                if v.orientation == Orientation::Forward {
                    if shape_is_null(&new_v1) {
                        new_v1 = v;
                    }
                } else if v.orientation == Orientation::Reversed {
                    if shape_is_null(&new_v2) {
                        new_v2 = v;
                    }
                } else if shape_is_null(v1) && shape_is_null(v2) {
                    a_nm_vertices.push(v);
                }
            }
        }
        new_v1.orientation = Orientation::Forward;
        new_v2.orientation = Orientation::Reversed;

        // szv#4:S4163:12Mar99 SGI warns — TopoDS_Shape sh = edge.EmptyCopied();
        // TopoDS_Edge E = TopoDS::Edge(sh).
        let e = brep.empty_copied(edge);

        if !shape_is_null(&new_v1) {
            builder_add(brep, &e, &new_v1);
        }
        if !shape_is_null(&new_v2) {
            builder_add(brep, &e, &new_v2);
        }

        // Addition of the internal or external vertices to edge.
        for i in 0..a_nm_vertices.len() {
            let nv = a_nm_vertices[i].clone();
            builder_add(brep, &e, &nv);
        }

        // S4054, rln 17.11.98 annie_surf.igs entity D77: 3D and pcurve have
        // different ranges, after B.Range all the ranges become as 3D.
        self.copy_ranges(brep, &e, edge, 0.0, 1.0);
        e
    }

    /// OCCT ShapeBuild_Edge::CopyRanges (ShapeBuild_Edge.cxx L206-334): copies
    /// ranges for curve3d and all common pcurves from `fromedge` into
    /// `toedge`, scaled to `[first + alpha * len, first + beta * len]`
    /// (defaults alpha = 0, beta = 1).
    pub fn copy_ranges(
        &self,
        brep: &mut BRep,
        toedge: &Shape,
        fromedge: &Shape,
        alpha: f64,
        beta: f64,
    ) {
        // The from-walk over BRep_CurveRepresentation rows (bridge #5): the
        // 3D GCurve row stands for TEdgeData.curve/range, the
        // CurveOnSurface/CurveOnClosedSurface rows carry (pcurve, range) and
        // the owning face key. fromGC null checks: skip a 3D row without a
        // curve and a pcurve row without a pcurve (rcad value rows always
        // carry the curve). CurveOn2Surfaces rows are regularity rows —
        // "only 3d curves and pcurves are treated".
        // fromGC->IsCurve3D() branch: skip when the 3D curve is null
        // (fromGC->Curve3D().IsNull()) — no 3D row on the from edge.
        if brep.edge(fromedge.clone()).curve.is_none() {
            self.copy_ranges_pcurves(brep, toedge, fromedge, alpha, beta);
            return;
        }
        let (first, last) = fromedge_first_last_3d(brep, fromedge);
        let len = last - first;
        let new_f = first + alpha * len;
        let new_l = first + beta * len;

        // PTV: 22.03.2002 fix for edge range (test-m020306-v2.step Shell #665,
        // Faces #40110, #40239): the periodic-range check against the TO
        // curve.
        let to_curve = brep.edge(toedge.clone()).curve.clone();
        let mut new_f = new_f;
        let mut new_l = new_l;
        if let Some(a_crv3d) = &to_curve {
            // 15.11.2002 PTV OCC966.
            if is_periodic_curve3d(a_crv3d) {
                let a_period = curve3_period(a_crv3d);
                let (a_crv_f, a_crv_l) = curve3_first_last(a_crv3d);
                if ((new_f - a_crv_f).abs() > PCONFUSION && new_f < a_crv_f) || new_f >= a_crv_l {
                    let a_shift = adjust_by_period(new_f, 0.5 * (a_crv_f + a_crv_l), a_period);
                    new_f += a_shift;
                    new_l += a_shift;
                    // BRep_Builder().SameRange(toedge, false);
                    // BRep_Builder().SameParameter(toedge, false);
                    let ed = brep.edge_mut_inplace(toedge.clone());
                    ed.same_range = false;
                    ed.same_parameter = false;
                }
            }
        }
        // toGC->SetRange(newF, newL) on the 3D representation.
        {
            let ed = brep.edge_mut_inplace(toedge.clone());
            if ed.curve.is_some() {
                ed.range = [new_f, new_l];
            }
        }
        self.copy_ranges_pcurves(brep, toedge, fromedge, alpha, beta);
    }

    /// The pcurve rows of the OCCT CopyRanges walk (the
    /// `!fromGC->IsCurve3D()` iterations of ShapeBuild_Edge.cxx L224-333).
    fn copy_ranges_pcurves(
        &self,
        brep: &mut BRep,
        toedge: &Shape,
        fromedge: &Shape,
        alpha: f64,
        beta: f64,
    ) {
        let from_reps: Vec<CurveRepresentation> =
            edge_representations(fromedge).to_vec();
        for from_gc in &from_reps {
            let (from_key, first, last) = match from_gc {
                CurveRepresentation::CurveOnSurface { face, range, .. } => (*face, range[0], range[1]),
                CurveRepresentation::CurveOnClosedSurface { face, range, .. } => (*face, range[0], range[1]),
                // Null-pcurve rows do not exist in rcad's value model;
                // regularity rows are skipped: only 3d curves and pcurves are
                // treated.
                _ => continue,
            };
            let len = last - first;
            let mut new_f = first + alpha * len;
            let mut new_l = first + beta * len;

            // The to-walk: match the representation with the same
            // (surface, location) — the owning face key stand-in.
            let to_reps: Vec<CurveRepresentation> = edge_representations(toedge).to_vec();
            for to_gc in &to_reps {
                let (to_key, to_pcurve) = match to_gc {
                    CurveRepresentation::CurveOnSurface { face, pcurve, .. } => (*face, Some(pcurve)),
                    CurveRepresentation::CurveOnClosedSurface { face, pcurve1, .. } => (*face, Some(pcurve1)),
                    // Not a 3d-curve/pcurve row (regularity row): continue.
                    _ => continue,
                };
                if to_key != from_key {
                    continue;
                }
                // PTV: 22.03.2002 periodic-range fix on the TO curve.
                if let Some(a_crv2d) = to_pcurve {
                    // 15.11.2002 PTV OCC966.
                    if is_periodic_curve2d(a_crv2d) {
                        let (a_crv_f, a_crv_l) = curve2d_first_last(a_crv2d);
                        let a_period = a_crv_l - a_crv_f;
                        if ((new_f - a_crv_f).abs() > PCONFUSION && new_f < a_crv_f)
                            || new_f >= a_crv_l
                        {
                            let a_shift =
                                adjust_by_period(new_f, 0.5 * (a_crv_f + a_crv_l), a_period);
                            new_f += a_shift;
                            new_l += a_shift;
                            // BRep_Builder().SameRange(toedge, false);
                            // BRep_Builder().SameParameter(toedge, false);
                            let ed = brep.edge_mut_inplace(toedge.clone());
                            ed.same_range = false;
                            ed.same_parameter = false;
                        }
                    }
                }
                // toGC->SetRange(newF, newL) on the pcurve row and its read
                // key (bridge #5: both mirrors).
                let ed = brep.edge_mut_inplace(toedge.clone());
                for r in ed.representations.iter_mut() {
                    match r {
                        CurveRepresentation::CurveOnSurface { face, range, .. } => {
                            if *face == to_key {
                                *range = [new_f, new_l];
                                break;
                            }
                        }
                        CurveRepresentation::CurveOnClosedSurface { face, range, .. } => {
                            if *face == to_key {
                                *range = [new_f, new_l];
                                break;
                            }
                        }
                        _ => {}
                    }
                }
                if let Some(entry) = ed.pcurves.get_mut(&to_key) {
                    entry.1 = new_f;
                    entry.2 = new_l;
                }
                break;
            }
        }
    }

    /// OCCT ShapeBuild_Edge::SetRange3d (ShapeBuild_Edge.cxx L338-356): sets
    /// the range on the 3d curve only.
    pub fn set_range3d(&self, brep: &mut BRep, edge: &Shape, first: f64, last: f64) {
        // The 3D representation row (bridge #5): TEdgeData.curve/range. The
        // walk breaks after the first (and only) 3D GCurve row.
        let has_c3d = brep.edge(edge.clone()).curve.is_some();
        if has_c3d {
            let ed = brep.edge_mut_inplace(edge.clone());
            ed.range = [first, last];
        }
    }

    /// OCCT ShapeBuild_Edge::CopyPCurves (ShapeBuild_Edge.cxx L360-413):
    /// makes copies of the pcurves from `fromedge` into `toedge`. Pcurves
    /// already present in `toedge` are replaced by copies, the others are
    /// copied. Ranges are also copied.
    pub fn copy_pcurves(&self, brep: &mut BRep, toedge: &Shape, fromedge: &Shape) {
        let from_loc = fromedge.location;
        let to_loc = toedge.location;
        let from_reps: Vec<CurveRepresentation> = edge_representations(fromedge).to_vec();
        for from_gc in &from_reps {
            // fromGC->IsCurveOnSurface() rows.
            let (from_key, from_pc1, from_pc2, from_range, from_closed) = match from_gc {
                CurveRepresentation::CurveOnSurface { face, pcurve, range } => {
                    (*face, Some(pcurve.clone()), None, *range, false)
                }
                CurveRepresentation::CurveOnClosedSurface { face, pcurve1, pcurve2, range } => {
                    (*face, Some(pcurve1.clone()), Some(pcurve2.clone()), *range, true)
                }
                _ => continue,
            };

            // The to-walk: find the representation on the same
            // (surface, location) — the face key stand-in.
            let to_reps: Vec<CurveRepresentation> = edge_representations(toedge).to_vec();
            let mut found: Option<usize> = None;
            for (i, to_gc) in to_reps.iter().enumerate() {
                let to_key = match to_gc {
                    CurveRepresentation::CurveOnSurface { face, .. }
                    | CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
                    _ => continue,
                };
                if to_key == from_key {
                    found = Some(i);
                }
                if found.is_some() {
                    break;
                }
            }

            // bug OCC209: invalid location of pcurve in the edge after
            // copying — newLoc = (fromLoc * L).Predivided(toLoc). Bridge #5:
            // the representation row stores its location only as the composed
            // hash, so the L component reduces to identity (the W1-1
            // wrapper-location reduction) and newLoc is composed from the
            // wrapper locations alone.
            let new_loc = pcurve_location_id(
                &(brep.get_location(from_loc) * brep.get_location(to_loc).inverse()),
            );
            let new_key = (from_key.0, new_loc);

            // OCCT: toGC = found ? the matched representation : a Copy of
            // fromGC appended to the to-list; then toGC->PCurve(pcurve copy)
            // (and PCurve2 for a closed-surface row) and toGC->Location(newLoc).
            let pc1 = from_pc1.expect("CurveOnSurface row carries a pcurve");
            let ed = brep.edge_mut_inplace(toedge.clone());
            match found {
                Some(i) => {
                    let old_key = match &ed.representations[i] {
                        CurveRepresentation::CurveOnSurface { face, .. }
                        | CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
                        _ => unreachable!("matched by key"),
                    };
                    ed.representations[i] = if from_closed {
                        CurveRepresentation::CurveOnClosedSurface {
                            face: new_key,
                            pcurve1: pc1.clone(),
                            pcurve2: from_pc2.clone().expect("closed row carries PCurve2"),
                            range: from_range,
                        }
                    } else {
                        CurveRepresentation::CurveOnSurface {
                            face: new_key,
                            pcurve: pc1.clone(),
                            range: from_range,
                        }
                    };
                    ed.pcurves.shift_remove(&old_key);
                    ed.pcurves.insert(new_key, (pc1.clone(), from_range[0], from_range[1]));
                }
                None => {
                    // toGC = fromGC->Copy(); tolist.Append(toGC).
                    ed.representations.push(if from_closed {
                        CurveRepresentation::CurveOnClosedSurface {
                            face: new_key,
                            pcurve1: pc1.clone(),
                            pcurve2: from_pc2.clone().expect("closed row carries PCurve2"),
                            range: from_range,
                        }
                    } else {
                        CurveRepresentation::CurveOnSurface {
                            face: new_key,
                            pcurve: pc1.clone(),
                            range: from_range,
                        }
                    });
                    ed.pcurves.insert(new_key, (pc1.clone(), from_range[0], from_range[1]));
                }
            }
        }
    }

    /// OCCT ShapeBuild_Edge::Copy (ShapeBuild_Edge.cxx L417-426): makes a copy
    /// of `edge` by call to CopyReplaceVertices (i.e. constructs a new TEdge
    /// with the same pcurves and vertices). If `sharepcurves` is false,
    /// pcurves are also replaced by their copies with help of CopyPCurves.
    pub fn copy(&self, brep: &mut BRep, edge: &Shape, sharepcurves: bool) -> Shape {
        let dummy1 = Shape::null();
        let dummy2 = Shape::null();
        let newedge = self.copy_replace_vertices(brep, edge, &dummy1, &dummy2);
        if !sharepcurves {
            self.copy_pcurves(brep, &newedge, edge);
        }
        newedge
    }

    /// OCCT ShapeBuild_Edge::RemovePCurve(edge, face) (ShapeBuild_Edge.cxx
    /// L430-443): removes the PCurve(s) recorded in the edge for the given
    /// face.
    pub fn remove_pcurve_face(&self, brep: &mut BRep, edge: &Shape, face: &Shape) {
        //: S4136  double tol = BRep_Tool::Tolerance ( edge );
        if brep.is_edge_closed_on_face(edge, face) {
            // B.UpdateEdge(edge, c2dNull, c2dNull, face, 0.).
            builder_update_edge_pcurves_null(brep, edge, face);
        } else {
            // B.UpdateEdge(edge, c2dNull, face, 0.).
            builder_update_edge_pcurve_null(brep, edge, face);
        }
    }

    /// OCCT ShapeBuild_Edge::RemovePCurve(edge, surf) (ShapeBuild_Edge.cxx
    /// L447-451): removes the PCurve(s) recorded in the edge for the given
    /// surface.
    pub fn remove_pcurve_surface(&self, brep: &mut BRep, edge: &Shape, surf: &Surface3) {
        self.remove_pcurve_surface_loc(brep, edge, surf, 0);
    }

    /// OCCT ShapeBuild_Edge::RemovePCurve(edge, surf, loc)
    /// (ShapeBuild_Edge.cxx L455-470): removes the PCurve(s) recorded in the
    /// edge for the given surface, with the given location.
    pub fn remove_pcurve_surface_loc(
        &self,
        brep: &mut BRep,
        edge: &Shape,
        surf: &Surface3,
        loc: u32,
    ) {
        //: S4136  double tol = BRep_Tool::Tolerance ( edge );
        if brep_tool_is_closed_on_surface(brep, edge, surf, loc) {
            builder_update_edge_pcurves_null_surface(brep, edge, surf, loc);
        } else {
            builder_update_edge_pcurve_null_surface(brep, edge, surf, loc);
        }
    }

    /// OCCT ShapeBuild_Edge::ReplacePCurve (ShapeBuild_Edge.cxx L474-503):
    /// replaces the PCurve in the edge for the given face. In case the edge
    /// is a seam (2 pcurves on that face), only the pcurve corresponding to
    /// the orientation of the edge is replaced.
    pub fn replace_pcurve(&self, brep: &mut BRep, edge: &Shape, pcurve: &Curve2d, face: &Shape) {
        // TopoDS_Shape dummy = edge.Reversed(); TopoDS_Edge edgerev =
        // TopoDS::Edge(dummy).
        let edgerev = Shape {
            orientation: reverse_orientation(edge.orientation),
            ..edge.clone()
        };
        // Reverse face to take the second pcurve for seams, like
        // SA_Edge::PCurve() does — TopoDS::Face(face.Oriented(TopAbs_FORWARD)).
        let f = Shape {
            orientation: Orientation::Forward,
            ..face.clone()
        };
        // BRep_Tool::CurveOnSurface(edge, F, f, l) and (edgerev, F, f, l).
        let (pcurve0, f0, l0) = match brep.curve_on_surface(edge, &f) {
            Some((c, a, b)) => (Some(c), a, b),
            None => (None, 0.0, 0.0),
        };
        let (c2d, _f1, _l1) = match brep.curve_on_surface(&edgerev, &f) {
            Some((c, a, b)) => (Some(c), a, b),
            None => (None, 0.0, 0.0),
        };
        let f_par = f0;
        let l_par = l0;
        // Add the pcurve to the edge (either as single, or as seam).
        match c2d {
            ref c2 if c2.is_none() || matches!(pcurve0, Some(ref p0) if curve2d_handle_same(c2.as_ref().unwrap(), p0)) => {
                // non-seam: B.UpdateEdge(edge, pcurve, face, 0)
                crate::shhealing::shape_build::edge::builder_update_edge_pcurve(
                    brep,
                    edge,
                    pcurve,
                    face,
                    0.0,
                );
            }
            c2d => {
                // seam
                if edge.orientation == Orientation::Forward {
                    // B.UpdateEdge(edge, pcurve, c2d, face, 0)
                    builder_update_edge_pcurves(
                        brep,
                        edge,
                        pcurve,
                        c2d.as_ref().unwrap(),
                        face,
                        f_par,
                        l_par,
                        0.0,
                    );
                } else {
                    // B.UpdateEdge(edge, c2d, pcurve, face, 0)
                    builder_update_edge_pcurves(
                        brep,
                        edge,
                        c2d.as_ref().unwrap(),
                        pcurve,
                        face,
                        f_par,
                        l_par,
                        0.0,
                    );
                }
            }
        }
        // B.Range(edge, face, f, l).
        builder_range_on_face(brep, edge, face, f_par, l_par);
    }

    /// OCCT ShapeBuild_Edge::ReassignPCurve (ShapeBuild_Edge.cxx L530-592):
    /// reassigns the edge pcurve lying on face `old` to another face `sub`.
    /// If the edge has two pcurves on `old`, only one is reassigned, the
    /// other is left alone. Similarly, if the edge already had a pcurve on
    /// `sub`, it will have two pcurves on it. Returns true on success, false
    /// when no pcurve lying on `old` was found.
    pub fn reassign_pcurve(&self, brep: &mut BRep, edge: &Shape, old: &Shape, sub: &Shape) -> bool {
        let mut npcurves = count_pcurves(brep, edge, old);
        // if ( npcurves <1 ) return false; //gka

        // pc = BRep_Tool::CurveOnSurface(edge, old, f, l).
        let (pc, f, l) = match brep.curve_on_surface(edge, old) {
            Some((c, a, b)) => (c, a, b),
            None => return false,
        };
        if npcurves == 0 {
            npcurves = 1; // gka
        }

        // If the pcurve was only one, remove; else leave the second one.
        if npcurves > 1 {
            // smh#8 Porting AIX — pc2 = BRep_Tool::CurveOnSurface(
            // edge.Reversed(), old, f, l).
            let erev = Shape {
                orientation: reverse_orientation(edge.orientation),
                ..edge.clone()
            };
            let (pc2, f2, l2) = match brep.curve_on_surface(&erev, old) {
                Some((c, a, b)) => (c, a, b),
                None => {
                    // OCCT would carry a null handle into UpdateEdge (the
                    // representation removal); the rcad value model has no
                    // null curve, so the removal form is used.
                    builder_update_edge_pcurve_null(brep, edge, old);
                    return false;
                }
            };
            let _ = (f2, l2);
            // B.UpdateEdge(edge, pc2, old, 0.); B.Range(edge, old, f, l).
            builder_update_edge_pcurve(brep, edge, &pc2, old, 0.0);
            builder_range_on_face(brep, edge, old, f, l);
        } else {
            self.remove_pcurve_face(brep, edge, old);
        }

        // If the edge does not have yet pcurves on sub, just add; else add as
        // first.
        let npcs = count_pcurves(brep, edge, sub);
        if npcs < 1 {
            // B.UpdateEdge(edge, pc, sub, 0.).
            builder_update_edge_pcurve(brep, edge, &pc, sub, 0.0);
        } else {
            // smh#8 Porting AIX — pcs = BRep_Tool::CurveOnSurface(
            // edge.Reversed(), sub, cf, cl).
            let erev = Shape {
                orientation: reverse_orientation(edge.orientation),
                ..edge.clone()
            };
            let (pcs, _cf, _cl) = match brep.curve_on_surface(&erev, sub) {
                Some((c, a, b)) => (c, a, b),
                None => {
                    builder_update_edge_pcurve(brep, edge, &pc, sub, 0.0);
                    builder_range_on_face(brep, edge, sub, f, l);
                    return true;
                }
            };
            if edge.orientation == Orientation::Reversed {
                // because B.UpdateEdge does not check edge orientation:
                // B.UpdateEdge(edge, pcs, pc, sub, 0.)
                builder_update_edge_pcurves(brep, edge, &pcs, &pc, sub, f, l, 0.0);
            } else {
                // B.UpdateEdge(edge, pc, pcs, sub, 0.)
                builder_update_edge_pcurves(brep, edge, &pc, &pcs, sub, f, l, 0.0);
            }
        }

        // B.Range(edge, sub, f, l).
        builder_range_on_face(brep, edge, sub, f, l);

        true
    }

    /// OCCT ShapeBuild_Edge::TransformPCurve (ShapeBuild_Edge.cxx L596-700):
    /// transforms the PCurve with the given matrix and affinity U factor.
    /// Returns the transformed curve; `a_first` / `a_last` are updated with
    /// the transformed parameters.
    pub fn transform_pcurve(
        &self,
        pcurve: &Curve2d,
        trans: &GpTrsf2d,
        u_fact: f64,
        a_first: &mut f64,
        a_last: &mut f64,
    ) -> Curve2d {
        // occ::handle<Geom2d_Curve> result = pcurve->Copy() — the value clone
        // (bridge #3).
        let mut result = pcurve.clone();
        if trans.form() != GpTrsfForm::Identity {
            result = geom2d_curve_transform(&result, trans);
            *a_first = geom2d_transformed_parameter(&result, *a_first, trans);
            *a_last = geom2d_transformed_parameter(&result, *a_last, trans);
        }
        if u_fact == 1.0 {
            return result;
        }

        // result->IsKind(Geom2d_TrimmedCurve): result = thecurve->BasisCurve().
        if let Curve2d::Trimmed(tc) = result {
            result = *tc.curve;
        }

        // gp_GTrsf2d tMatu; tMatu.SetAffinity(gp::OY2d(), uFact).
        let t_matu = GTrsf2dAffinity::oy2d(u_fact);

        if let Curve2d::Line(a_line2d) = result {
            // gp_Pnt2d Pf, Pl; aLine2d->D0(aFirst, Pf) — the line evaluation
            // Location + U * Direction.
            let mut pf = a_line2d.origin + a_line2d.direction * *a_first;
            pf = t_matu.transforms(pf);
            let mut pl = a_line2d.origin + a_line2d.direction * *a_last;
            pl = t_matu.transforms(pl);
            // gp_Lin2d line2d(Pf, gp_Dir2d(gp_Vec2d(Pf, Pl))).
            let line2d = Line2d::new(pf, pl - pf);
            // aFirst = ElCLib::Parameter(line2d, Pf); aLast = idem(Pl).
            *a_first = elclib_parameter_lin2d(&line2d, pf);
            *a_last = elclib_parameter_lin2d(&line2d, pl);
            // occ::handle<Geom2d_Line> Gline2d = new Geom2d_Line(line2d).
            return Curve2d::Line(line2d);
        }

        if let Curve2d::Bezier(mut bezier) = result {
            // Transform the poles of the Bezier curve: bezier->SetPole(i, Pt1).
            let nb_pol = bezier.control_points.len();
            for i in 0..nb_pol {
                let mut pxy = bezier.control_points[i];
                pxy = t_matu.transforms(pxy);
                bezier.control_points[i] = pxy;
            }
            return Curve2d::Bezier(bezier);
        }

        // The remaining branch: conics are approximated / converted, plain
        // BSplines are transformed in place.
        let a_bspline2d: Option<BSplineCurve2> = if result.is_conic() {
            // occ::handle<Geom2d_Curve> tcurve = new Geom2d_TrimmedCurve(
            // result, aFirst, aLast) — protection against parabolas etc.
            let tcurve = Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(result.clone()),
                t_min: *a_first,
                t_max: *a_last,
            });
            // Geom2dConvert_ApproxCurve approx(tcurve, Precision::Approximation(),
            // GeomAbs_C1, 100, 6) — the rcad approx carries the tolerance
            // argument only (GAP: the continuity/degree/segment parameters
            // are pending the TKGeomBase batch); HasResult() == Some.
            match approx_curve_to_bspline(&tcurve, APPROXIMATION) {
                Some(approx) => Some(approx),
                None => {
                    // aBSpline2d = Geom2dConvert::CurveToBSplineCurve(
                    // tcurve, Convert_QuasiAngular).
                    geom2d_convert_curve_to_bspline(&tcurve)
                }
            }
        } else {
            match &result {
                Curve2d::BSpline(bs) => Some(bs.clone()),
                other => geom2d_convert_curve_to_bspline(other),
            }
        };

        match a_bspline2d {
            Some(mut bs) => {
                if result.is_conic() {
                    // aFirst = aBSpline2d->FirstParameter(); aLast = idem.
                    let (cf, cl) = bspline2d_first_last(&bs);
                    *a_first = cf;
                    *a_last = cl;
                }
                // Transform the poles of the BSplineCurve:
                // aBSpline2d->SetPole(i, Pt1).
                let nb_pol = bs.control_points.len();
                for i in 0..nb_pol {
                    let mut pxy = bs.control_points[i];
                    pxy = t_matu.transforms(pxy);
                    bs.control_points[i] = pxy;
                }
                Curve2d::BSpline(bs)
            }
            // GAP failure path: the Geom2dConvert conversion produced no
            // result; OCCT always produces one, so this returns the
            // unconverted (affinity-free) curve.
            None => result,
        }
    }

    /// OCCT ShapeBuild_Edge::RemoveCurve3d (ShapeBuild_Edge.cxx L704-710):
    /// removes the Curve3D recorded in the edge.
    pub fn remove_curve_3d(&self, brep: &mut BRep, edge: &Shape) {
        //: S4136  double tol = BRep_Tool::Tolerance (edge);
        // B.UpdateEdge(edge, c3dNull, 0.).
        let ed = brep.edge_mut_inplace(edge.clone());
        ed.curve = None;
        ed.representations
            .retain(|r| !matches!(r, CurveRepresentation::Curve3D { .. }));
        ed.tolerance = ed.tolerance.max(0.0);
    }

    /// OCCT ShapeBuild_Edge::BuildCurve3d (ShapeBuild_Edge.cxx L714-774):
    /// calls BRepTools::BuildCurve3D.
    pub fn build_curve3d(&self, brep: &mut BRep, edge: &Shape) -> bool {
        // OCC_CATCH_SIGNALS try block; rcad has no signal machinery.
        // #48 rln 10.12.98 S4054 UKI60107-5 entity 365: use the maximum of
        // the tolerance and the default parameter 1.e-5.
        if brep_lib_build_curve3d(brep, edge, 1.0e-5_f64.max(brep.tolerance(edge))) {
            // #50 S4054 rln 14.12.98: pcurve and removed 3D curves have
            // different ranges — with SameRange, set the range explicitly for
            // all representations.
            if brep.edge_same_range(edge) {
                let [first, last] = brep.edge_range(edge);
                // BRep_Builder().Range(edge, first, last).
                builder_range_all(brep, edge, first, last);
            }
            // c3d = BRep_Tool::Curve(edge, f, l).
            let (c3d, f, l) = match brep.tshapes[edge.index].as_ref() {
                TShape::Edge(ed) => match &ed.curve {
                    Some(c) => (c.clone(), ed.range[0], ed.range[1]),
                    None => return false,
                },
                _ => return false,
            };
            // 15.11.2002 PTV OCC966.
            if !is_periodic_curve3d(&c3d) {
                let mut is_less = false;
                let mut f = f;
                let mut l = l;
                let (cf, cl) = curve3_first_last(&c3d);
                if f < cf {
                    is_less = true;
                    f = cf;
                }
                if l > cl {
                    is_less = true;
                    l = cl;
                }
                if is_less {
                    self.set_range3d(brep, edge, f, l);
                    // BRep_Builder().SameRange(edge, false).
                    let ed = brep.edge_mut_inplace(edge.clone());
                    ed.same_range = false;
                }
            }

            return true;
        }
        false
    }

    /// OCCT ShapeBuild_Edge::MakeEdge(edge, curve, L) (ShapeBuild_Edge.cxx
    /// L778-783): makes an edge with curve and location.
    pub fn make_edge_curve_loc(
        &self,
        brep: &mut BRep,
        edge: &mut Shape,
        curve: &Curve3,
        l: u32,
    ) {
        // curve->FirstParameter(), curve->LastParameter().
        let (p1, p2) = curve3_first_last(curve);
        self.make_edge_curve_loc_params(brep, edge, curve, l, p1, p2);
    }

    /// OCCT ShapeBuild_Edge::MakeEdge(edge, curve, L, p1, p2)
    /// (ShapeBuild_Edge.cxx L787-815): makes an edge with curve, location and
    /// range [p1, p2].
    pub fn make_edge_curve_loc_params(
        &self,
        brep: &mut BRep,
        edge: &mut Shape,
        curve: &Curve3,
        l: u32,
        p1: f64,
        p2: f64,
    ) {
        // BRepBuilderAPI_MakeEdge ME(curve, p1, p2).
        let me = match brep_builder_api_make_edge_3d(brep, curve, p1, p2) {
            Some(me) => me,
            None => {
                // OCCT_DEBUG warning "\nWarning: ShapeBuild_Edge::MakeEdge
                // BRepAPI_NotDone" suppressed (the OCCT_DEBUG block).
                return;
            }
        };
        let e = me.edge;
        if l != 0 {
            // L.IsIdentity() is false: B.UpdateEdge(E, curve, L, 0.);
            // B.Range(E, p1, p2).
            {
                let ed = brep.edge_mut_inplace(e.clone());
                ed.curve = Some(curve.clone());
                ed.range = [p1, p2];
            }
            // TopExp::Vertices(E, V1, V2): the FORWARD and REVERSED vertices.
            let (v1, v2) = {
                let ed = brep.edge(e.clone());
                (ed.first.clone(), ed.last.clone())
            };
            // gp_Pnt P1 = BRep_Tool::Pnt(V1), P2 = BRep_Tool::Pnt(V2);
            // B.UpdateVertex(V1, P1.Transformed(L.Transformation()), 0.) and
            // idem for V2.
            let p1 = brep.vertex_position(&v1);
            let p2 = brep.vertex_position(&v2);
            let trsf = brep.get_location(l);
            let mut builder = make_builder(brep);
            builder.update_vertex_point(brep, v1, trsf.transform_point3(p1), 0.0);
            builder.update_vertex_point(brep, v2, trsf.transform_point3(p2), 0.0);
        }
        // edge = E.
        *edge = e;
    }

    /// OCCT ShapeBuild_Edge::MakeEdge(edge, pcurve, face) (ShapeBuild_Edge.cxx
    /// L819-823): makes an edge with pcurve and face.
    pub fn make_edge_pcurve_face(
        &self,
        brep: &mut BRep,
        edge: &mut Shape,
        pcurve: &Curve2d,
        face: &Shape,
    ) {
        // pcurve->FirstParameter(), pcurve->LastParameter().
        let (p1, p2) = curve2d_first_last(pcurve);
        self.make_edge_pcurve_face_params(brep, edge, pcurve, face, p1, p2);
    }

    /// OCCT ShapeBuild_Edge::MakeEdge(edge, pcurve, face, p1, p2)
    /// (ShapeBuild_Edge.cxx L828-837): makes an edge with pcurve, face and
    /// range [p1, p2].
    pub fn make_edge_pcurve_face_params(
        &self,
        brep: &mut BRep,
        edge: &mut Shape,
        pcurve: &Curve2d,
        face: &Shape,
        p1: f64,
        p2: f64,
    ) {
        // TopLoc_Location L; const occ::handle<Geom_Surface>& S =
        // BRep_Tool::Surface(face, L).
        let (s, l) = match brep.tshapes[face.index].as_ref() {
            TShape::Face(fd) => match &fd.surface {
                Some(s) => (s.clone(), fd.surface_location),
                None => return,
            },
            _ => return,
        };
        self.make_edge_pcurve_surface_loc_params(brep, edge, pcurve, &s, l, p1, p2);
    }

    /// OCCT ShapeBuild_Edge::MakeEdge(edge, pcurve, S, L) (ShapeBuild_Edge.cxx
    /// L841-847): makes an edge with pcurve, surface and location.
    pub fn make_edge_pcurve_surface_loc(
        &self,
        brep: &mut BRep,
        edge: &mut Shape,
        pcurve: &Curve2d,
        s: &Surface3,
        l: u32,
    ) {
        let (p1, p2) = curve2d_first_last(pcurve);
        self.make_edge_pcurve_surface_loc_params(brep, edge, pcurve, s, l, p1, p2);
    }

    /// OCCT ShapeBuild_Edge::MakeEdge(edge, pcurve, S, L, p1, p2)
    /// (ShapeBuild_Edge.cxx L851-881): makes an edge with pcurve, surface,
    /// location and range [p1, p2].
    pub fn make_edge_pcurve_surface_loc_params(
        &self,
        brep: &mut BRep,
        edge: &mut Shape,
        pcurve: &Curve2d,
        s: &Surface3,
        l: u32,
        p1: f64,
        p2: f64,
    ) {
        // BRepBuilderAPI_MakeEdge ME(pcurve, S, p1, p2).
        let me = match brep_builder_api_make_edge_pcurve(brep, pcurve, s, p1, p2) {
            Some(me) => me,
            None => {
                // OCCT_DEBUG warning suppressed.
                return;
            }
        };
        let e = me.edge;
        if l != 0 {
            // RemovePCurve(E, S).
            self.remove_pcurve_surface(brep, &e, s);
            // B.UpdateEdge(E, pcurve, S, L, 0.); B.Range(E, S, L, p1, p2).
            {
                // GAP/architecture: OCCT keys the pcurve by (S, L); the rcad
                // pcurve map is face-keyed, so the surfaceless representation
                // uses the reserved (0, 0) key ("no face"). The W3
                // ShapeFix_Edge consumer revisits this key.
                let key_loc = pcurve_location_id(&brep.get_location(l));
                let ed = brep.edge_mut_inplace(e.clone());
                ed.pcurves.insert((0, key_loc), (pcurve.clone(), p1, p2));
            }
            // TopExp::Vertices(E, V1, V2); P1/P2 world points transformed by
            // L; B.UpdateVertex(V1/V2, ..., 0.).
            let (v1, v2) = {
                let ed = brep.edge(e.clone());
                (ed.first.clone(), ed.last.clone())
            };
            let p1 = brep.vertex_position(&v1);
            let p2 = brep.vertex_position(&v2);
            let trsf = brep.get_location(l);
            let mut builder = make_builder(brep);
            builder.update_vertex_point(brep, v1, trsf.transform_point3(p1), 0.0);
            builder.update_vertex_point(brep, v2, trsf.transform_point3(p2), 0.0);
        }
        // edge = E.
        *edge = e;
    }
}

/// OCCT TopAbs orientation reversal (TopoDS_Shape::Reversed).
fn reverse_orientation(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// The TEdgeData 3D range read (the fromGC->First()/Last() of the 3D row).
fn fromedge_first_last_3d(brep: &BRep, e: &Shape) -> (f64, f64) {
    match brep.tshapes[e.index].as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// Geom_BSplineCurve::FirstParameter/LastParameter — the clamped knot range.
fn bspline2d_first_last(bs: &BSplineCurve2) -> (f64, f64) {
    let n = bs.knots.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let d = bs.degree.min(n - 1);
    (bs.knots[d], bs.knots[n - 1 - d])
}

/// A fresh `BRep_Builder` stand-in threaded over the pool (bridge #1).
fn make_builder(brep: &mut BRep) -> rcad_kernel::topo::topods::BRepBuilder {
    let _ = brep;
    rcad_kernel::topo::topods::BRepBuilder::new()
}

// ---------------------------------------------------------------------------
// BRep_Builder re-hosts (architecture bridge #1; BRep_Builder.cxx anchors).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::UpdateEdge(E, C, F, Tol) with a null C
/// (BRep_Builder.cxx L660-700): the single-pcurve representation on the face
/// is removed from the edge TShape. rcad removes both mirrors (bridge #5) for
/// every wrapper-location key of the edge.
pub fn builder_update_edge_pcurve_null(brep: &mut BRep, edge: &Shape, face: &Shape) {
    let locs = brep.edge_wrapper_locations(edge);
    let keys: Vec<(u64, u32)> = locs
        .iter()
        .map(|&el| {
            (
                face.ptr_id(),
                rcad_kernel::topo::topods::compose_pcurve_location(
                    face.location,
                    el,
                    &brep.locations,
                ),
            )
        })
        .collect();
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        ed.pcurves.shift_remove(k);
    }
    ed.representations.retain(|r| match r {
        CurveRepresentation::CurveOnSurface { face: f, .. } => !keys.contains(f),
        _ => true,
    });
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, F, Tol) with null curves
/// (BRep_Builder.cxx L703-757): the closed-surface (seam) representation on
/// the face is removed from the edge TShape.
pub fn builder_update_edge_pcurves_null(brep: &mut BRep, edge: &Shape, face: &Shape) {
    let locs = brep.edge_wrapper_locations(edge);
    let keys: Vec<(u64, u32)> = locs
        .iter()
        .map(|&el| {
            (
                face.ptr_id(),
                rcad_kernel::topo::topods::compose_pcurve_location(
                    face.location,
                    el,
                    &brep.locations,
                ),
            )
        })
        .collect();
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        ed.pcurves.shift_remove(k);
    }
    ed.representations.retain(|r| match r {
        CurveRepresentation::CurveOnClosedSurface { face: f, .. } => !keys.contains(f),
        _ => true,
    });
}

/// OCCT BRep_Builder::UpdateEdge(E, C, S, L, Tol) with a null C: the pcurve
/// representation on the (surface, location) pair is removed. rcad resolves
/// the (S, L) pair through the faces registered on the surface (the
/// surface_same stand-in).
pub fn builder_update_edge_pcurve_null_surface(
    brep: &mut BRep,
    edge: &Shape,
    surf: &Surface3,
    loc: u32,
) {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        ed.pcurves.shift_remove(k);
    }
    ed.representations.retain(|r| match r {
        CurveRepresentation::CurveOnSurface { face: f, .. } => !keys.contains(f),
        _ => true,
    });
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) with null curves: the
/// closed-surface representation on the (surface, location) pair is removed.
pub fn builder_update_edge_pcurves_null_surface(
    brep: &mut BRep,
    edge: &Shape,
    surf: &Surface3,
    loc: u32,
) {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        ed.pcurves.shift_remove(k);
    }
    ed.representations.retain(|r| match r {
        CurveRepresentation::CurveOnClosedSurface { face: f, .. } => !keys.contains(f),
        _ => true,
    });
}

/// The (face ptr, composed location) keys of the edge's pcurve rows whose
/// registered face surface matches (surf, loc). OCCT keys the
/// representations by the (Geom_Surface, TopLoc_Location) pair; rcad stores
/// the owning face key — the surface identity resolves through the face
/// (the surface_same handle stand-in) and the location filter through the
/// composed pcurve-location hash (the loc argument participates through the
/// face registration, the W1-1 wrapper-location reduction).
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
                    rcad_kernel::topo::topods::compose_pcurve_location(
                        face.location,
                        el,
                        &brep.locations,
                    ),
                ));
            }
        }
    }
    keys
}

/// OCCT BRep_Builder::UpdateEdge(E, C, F, Tol) (BRep_Builder.cxx L660-700):
/// set the pcurve on the face (thin wrapper over the kernel builder method,
/// kept local so the walk above reads one consistent form).
pub fn builder_update_edge_pcurve(
    brep: &mut BRep,
    edge: &Shape,
    pcurve: &Curve2d,
    face: &Shape,
    tol: f64,
) {
    let mut b = make_builder(brep);
    b.update_edge_pcurve(brep, edge.clone(), pcurve.clone(), face.clone(), tol);
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, F, aFirst, aLast, Tol)
/// (BRep_Builder.cxx L703-757): set the seam pcurve pair on the face.
pub fn builder_update_edge_pcurves(
    brep: &mut BRep,
    edge: &Shape,
    pcurve1: &Curve2d,
    pcurve2: &Curve2d,
    face: &Shape,
    a_first: f64,
    a_last: f64,
    tol: f64,
) {
    let mut b = make_builder(brep);
    b.update_edge_pcurve_closed(
        brep,
        edge.clone(),
        pcurve1.clone(),
        pcurve2.clone(),
        face.clone(),
        a_first,
        a_last,
        tol,
    );
}

/// OCCT BRep_Builder::Range(E, F, f, l) (BRep_Builder.cxx L558-600): sets the
/// range of the pcurve representation of the edge on the face.
pub fn builder_range_on_face(brep: &mut BRep, edge: &Shape, face: &Shape, f: f64, l: f64) {
    let locs = brep.edge_wrapper_locations(edge);
    let keys: Vec<(u64, u32)> = locs
        .iter()
        .map(|&el| {
            (
                face.ptr_id(),
                rcad_kernel::topo::topods::compose_pcurve_location(
                    face.location,
                    el,
                    &brep.locations,
                ),
            )
        })
        .collect();
    let ed = brep.edge_mut_inplace(edge.clone());
    for k in &keys {
        if let Some(entry) = ed.pcurves.get_mut(k) {
            entry.1 = f;
            entry.2 = l;
        }
    }
    for r in ed.representations.iter_mut() {
        match r {
            CurveRepresentation::CurveOnSurface { face: fk, range, .. } => {
                if keys.contains(fk) {
                    *range = [f, l];
                }
            }
            CurveRepresentation::CurveOnClosedSurface { face: fk, range, .. } => {
                if keys.contains(fk) {
                    *range = [f, l];
                }
            }
            _ => {}
        }
    }
}

/// OCCT BRep_Builder::Range(E, f, l) (BRep_Builder.cxx L536-556): sets the
/// range on all curve representations of the edge ("explicit setting for all
/// reps").
pub fn builder_range_all(brep: &mut BRep, edge: &Shape, f: f64, l: f64) {
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

/// OCCT BRep_Tool::IsClosed(E, S, L) (BRep_Tool.cxx): true when the edge
/// carries two pcurves on the (surface, location) pair — the
/// CurveOnClosedSurface representation stand-in resolved through the faces
/// registered on the surface.
fn brep_tool_is_closed_on_surface(brep: &BRep, edge: &Shape, surf: &Surface3, loc: u32) -> bool {
    let keys = surface_face_keys(brep, edge, surf, loc);
    let reps = edge_representations(edge);
    for r in reps {
        if let CurveRepresentation::CurveOnClosedSurface { face, .. } = r {
            if keys.contains(face) {
                return true;
            }
        }
    }
    let _ = loc;
    false
}

// ---------------------------------------------------------------------------
// GAP carriers (other untranslated packages; dependency + anchor + OCCT
// failure path).
// ---------------------------------------------------------------------------

/// OCCT TKTopAlgo/BRepLib/BRepLib.cxx L301-456 — BuildCurve3d(AnEdge,
/// Tolerance, ...). "if the edge has a 3d curve returns true" (L319-325). The
/// pcurve reconstruction branches (L330-454) route through GeomLib::
/// BuildCurve3d — GAP: pending the TKTopAlgo/GeomLib approximation batch
/// (docket section 4, gap 3); OCCT returns false when the 3d curve cannot be
/// produced (L363-366 / L436-439 / L452).
fn brep_lib_build_curve3d(brep: &BRep, edge: &Shape, tolerance: f64) -> bool {
    let _ = (edge, tolerance);
    match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => ed.curve.is_some(),
        _ => false,
    }
}

/// OCCT gp_TrsfForm (gp_Trsf.hxx) — the forms the TransformPCurve path
/// distinguishes (gp_Identity against the rest).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpTrsfForm {
    Identity,
    Rotation,
    Translation,
    Scale,
    Mirror,
    MirrorPoint,
    MirrorLine,
    CompoundTrsf,
    Other,
}

/// GAP carrier for OCCT `gp_Trsf2d` (TKMath, `gp_Trsf.hxx` / `gp_Trsf.cxx`):
/// the 2D transformation with a separate scale factor and form tag. The
/// shape_extend GAP carrier (composite_surface.rs) keeps its fields private
/// and is owned by another batch, so the ShapeBuild package carries its own
/// instance; both collapse into the TKMath batch when it lands. Point
/// transformation follows gp_Trsf::Transformed: `P' = matrix * P * scale +
/// loc`; direction transformation follows gp_Dir2d::Transform: the matrix
/// part only (no translation).
#[derive(Debug, Clone, Copy)]
pub struct GpTrsf2d {
    /// Row-major 2x2 matrix (WITHOUT the scale folded in, like OCCT).
    pub matrix: [[f64; 2]; 2],
    /// Translation part.
    pub loc: [f64; 2],
    /// gp_Trsf scale factor.
    pub scale: f64,
    /// gp_TrsfForm tag.
    form: GpTrsfForm,
}

impl Default for GpTrsf2d {
    /// OCCT gp_Trsf default constructor: identity.
    fn default() -> Self {
        GpTrsf2d {
            matrix: [[1.0, 0.0], [0.0, 1.0]],
            loc: [0.0, 0.0],
            scale: 1.0,
            form: GpTrsfForm::Identity,
        }
    }
}

impl GpTrsf2d {
    /// OCCT gp_Trsf::Form().
    pub fn form(&self) -> GpTrsfForm {
        self.form
    }

    /// OCCT gp_Trsf::ScaleFactor().
    pub fn scale_factor(&self) -> f64 {
        self.scale
    }

    /// OCCT gp_Trsf2d::Transformed(gp_XY): `matrix * P * scale + loc`.
    pub fn transformed(&self, p: glam::DVec2) -> glam::DVec2 {
        let m = self.matrix;
        glam::DVec2::new(
            (m[0][0] * p.x + m[0][1] * p.y) * self.scale + self.loc[0],
            (m[1][0] * p.x + m[1][1] * p.y) * self.scale + self.loc[1],
        )
    }

    /// The matrix part only (gp_Dir2d::Transform — gp_Dir2d.cxx: the
    /// direction is transformed by the matrix and scale, then normalized).
    pub fn transform_dir(&self, d: glam::DVec2) -> glam::DVec2 {
        let m = self.matrix;
        glam::DVec2::new(
            (m[0][0] * d.x + m[0][1] * d.y) * self.scale,
            (m[1][0] * d.x + m[1][1] * d.y) * self.scale,
        )
        .normalize_or_zero()
    }
}

/// OCCT Geom2d_Curve::Transform dispatch (TKGeomBase/Geom2d). GAP carrier:
/// only the Geom2d_Line form is translated natively (Geom2d_Line.cxx L239-242:
/// pos.Transform(T)); the remaining Geom2d_* variants (Conics, Bezier,
/// BSpline, Offset) are pending the TKGeomBase batch and return the input
/// unchanged (the documented reduced behavior — OCCT never fails here).
fn geom2d_curve_transform(pcurve: &Curve2d, trans: &GpTrsf2d) -> Curve2d {
    match pcurve {
        Curve2d::Line(l) => Curve2d::Line(Line2d {
            origin: trans.transformed(l.origin),
            direction: trans.transform_dir(l.direction),
        }),
        // GAP: Geom2d_Conic/Bezier/BSpline/Offset ::Transform pending.
        other => other.clone(),
    }
}

/// OCCT Geom2d_Curve::TransformedParameter dispatch (TKGeomBase/Geom2d). The
/// base returns U unchanged (Geom2d_Curve.cxx L41-44); the Geom2d_Line
/// override returns U * |T.ScaleFactor()| with infinite U passed through
/// (Geom2d_Line.cxx L246-253). GAP: the remaining per-variant overrides
/// (Conics / BSpline / Bezier / Trimmed) are pending the TKGeomBase batch and
/// fall back to the base form.
fn geom2d_transformed_parameter(pcurve: &Curve2d, u: f64, trans: &GpTrsf2d) -> f64 {
    let is_line = matches!(pcurve, Curve2d::Line(_));
    if is_line {
        // OCCT Geom2d_Line.cxx L246-252: Precision::IsInfinite(U) passthrough.
        if rcad_kernel::precision::is_infinite_value(u) {
            return u;
        }
        return u * trans.scale_factor().abs();
    }
    u
}

/// OCCT Geom2dConvert::CurveToBSplineCurve(C, Convert_QuasiAngular)
/// (TKGeomBase/Geom2dConvert). GAP carrier: the conversion is pending the
/// TKGeomBase batch; the reduced form routes through the kernel 2D approx
/// (the same reduction the Geom2dConvert_ApproxCurve carrier uses) and
/// returns None when no approximation is produced — OCCT always produces a
/// curve, so callers treat None as the documented failure path.
fn geom2d_convert_curve_to_bspline(curve: &Curve2d) -> Option<BSplineCurve2> {
    approx_curve_to_bspline(curve, APPROXIMATION)
}

/// OCCT ElCLib::Parameter(L, P) for a gp_Lin2d (ElCLib.cxx L1276-1281,
/// LineParameter): the parameter of P on the line = (P - Location) . Direction
/// (the direction is unit).
fn elclib_parameter_lin2d(line: &Line2d, p: glam::DVec2) -> f64 {
    (p - line.origin).dot(line.direction)
}

/// OCCT gp_GTrsf2d::SetAffinity(gp::OY2d(), Ratio) (gp_GTrsf2d.cxx L24-38):
/// the affinity of ratio Ratio with respect to the OY axis. With the axis
/// direction (0, 1) the matrix evaluates to [[Ratio, 0], [0, 1]] and the
/// translation part vanishes, so Transforms(XY) scales the X component.
#[derive(Debug, Clone, Copy)]
struct GTrsf2dAffinity {
    ratio: f64,
}

impl GTrsf2dAffinity {
    fn oy2d(ratio: f64) -> Self {
        GTrsf2dAffinity { ratio }
    }

    /// OCCT gp_GTrsf2d::Transforms(XY) for this affinity.
    fn transforms(&self, p: glam::DVec2) -> glam::DVec2 {
        glam::DVec2::new(self.ratio * p.x, p.y)
    }
}

// ---------------------------------------------------------------------------
// BRepBuilderAPI_MakeEdge re-hosts (TKTopAlgo; the docket defers the
// standalone error-status Make* port decision to W4 — docket section 4).
// ---------------------------------------------------------------------------

/// The consumed state of BRepLib_MakeEdge after Init (BRepLib_MakeEdge.hxx
/// L254-256: myError / myVertex1 / myVertex2 over myShape).
struct BRepBuilderAPIMakeEdge {
    edge: Shape,
}

/// OCCT BRepLib_MakeEdge::Init(C, V1, V2, pp1, pp2)
/// (BRepLib_MakeEdge.cxx L603-793) reduced to the consumed constructor
/// BRepBuilderAPI_MakeEdge(curve, p1, p2) — both vertices null. Returns None
/// for every BRepLib_*Error branch (!IsDone), keeping the OCCT failure path.
fn brep_builder_api_make_edge_3d(brep: &mut BRep, cc: &Curve3, pp1: f64, pp2: f64) -> Option<BRepBuilderAPIMakeEdge> {
    // Kill trimmed curves.
    let mut c = cc.clone();
    loop {
        match &c {
            Curve3::Trimmed(tc) => c = tc.basis_curve().clone(),
            _ => break,
        }
    }

    // Check parameters.
    let mut p1 = pp1;
    let mut p2 = pp2;
    let (cf, cl) = curve3_first_last(&c);
    let epsilon = PCONFUSION;
    let periodic = is_periodic_curve3d(&c);

    let mut v1 = Shape::null();
    let mut v2 = Shape::null();
    if periodic {
        // Adjust in period.
        elclib_adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
    } else {
        // Reordonate.
        if p1 < p2 {
            // V1 = VV1; V2 = VV2 (both null here).
        } else {
            let x = p1;
            p1 = p2;
            p2 = x;
        }

        // Check range.
        if (cf - p1 > epsilon) || (p2 - cl > epsilon) {
            // BRepLib_ParameterOutOfRange.
            return None;
        }

        // Check punctuality.
        if (p2 - p1) <= GP_RESOLUTION {
            // BRepLib_LineThroughIdenticPoints.
            return None;
        }
    }

    // Compute points on the curve (GeomAdaptor_Curve::Value).
    let p1inf = is_negative_infinite_value(p1);
    let p2inf = is_positive_infinite_value(p2);
    let mut pt1 = glam::DVec3::ZERO;
    let mut pt2 = glam::DVec3::ZERO;
    if !p1inf {
        pt1 = c.point_at(p1);
    }
    if !p2inf {
        pt2 = c.point_at(p2);
    }

    let preci = CONFUSION; // BRepLib::Precision().

    // Check for closed curve.
    let closed;
    // OCCT L769-778 sets degenerated in the closed branch with provided
    // vertices (cut in this reduced re-host), so it stays false here.
    let degenerated = false;
    if !p1inf && !p2inf {
        closed = pt1.distance(pt2) <= preci;
    } else {
        closed = false;
    }

    // Check if the vertices are on the curve (both null here: the closed
    // branch makes V2 = V1 = MakeVertex(P1, preci)).
    if closed {
        v1 = brep.add_tvertex_unique(pt1);
        brep.vertex_mut(v1.clone()).tolerance = preci;
        v2 = v1.clone();
    } else {
        if !p1inf {
            // B.MakeVertex(V1, P1, preci).
            v1 = brep.add_tvertex_unique(pt1);
            brep.vertex_mut(v1.clone()).tolerance = preci;
        }
        if !p2inf {
            // B.MakeVertex(V2, P2, preci).
            v2 = brep.add_tvertex_unique(pt2);
            brep.vertex_mut(v2.clone()).tolerance = preci;
        }
    }

    v1.orientation = Orientation::Forward;
    v2.orientation = Orientation::Reversed;

    // TopoDS_Edge& E = TopoDS::Edge(myShape); B.MakeEdge(E, C, preci).
    let e = brep.add_tedge(Some(c), v1.clone(), v2.clone(), [p1, p2]);
    // B.Add(E, V1) / B.Add(E, V2): add_tedge carries the vertices in the
    // first/last slots (the OCCT Add order); B.Range(E, p1, p2) is the range
    // argument of add_tedge.
    // B.Degenerated(E, degenerated); E.Closed(closed).
    {
        let ed = brep.edge_mut_inplace(e.clone());
        ed.degenerated = degenerated;
    }
    if closed {
        crate::shhealing::shape_build::brep_tool::set_flag_inplace(
            brep,
            &e,
            rcad_kernel::topo::topods::tshape_flags::CLOSED,
            true,
        );
    }

    Some(BRepBuilderAPIMakeEdge { edge: e })
}

/// OCCT BRepLib_MakeEdge::Init(C, S, V1, V2, pp1, pp2)
/// (BRepLib_MakeEdge.cxx L904-1075) reduced to the consumed constructor
/// BRepBuilderAPI_MakeEdge(pcurve, S, p1, p2) — both vertices null. Returns
/// None for every BRepLib_*Error branch (!IsDone).
fn brep_builder_api_make_edge_pcurve(
    brep: &mut BRep,
    cc: &Curve2d,
    s: &Surface3,
    pp1: f64,
    pp2: f64,
) -> Option<BRepBuilderAPIMakeEdge> {
    // Kill trimmed curves.
    let mut c = cc.clone();
    loop {
        match &c {
            Curve2d::Trimmed(tc) => c = tc.curve.as_ref().clone(),
            _ => break,
        }
    }

    // Check parameters.
    let mut p1 = pp1;
    let mut p2 = pp2;
    let (cf, cl) = curve2d_first_last(&c);
    let epsilon = PCONFUSION;
    let periodic = is_periodic_curve2d(&c);

    let mut reverse = false;
    if periodic {
        // Adjust in period.
        elclib_adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
    } else {
        // Reordonate.
        if p1 >= p2 {
            let x = p1;
            p1 = p2;
            p2 = x;
            reverse = true;
        }

        // Check range.
        if (cf - p1 > epsilon) || (p2 - cl > epsilon) {
            // BRepLib_ParameterOutOfRange.
            return None;
        }
    }

    // Compute points on the curve: P2d = C->Value(p); P = S->Value(P2d.X(),
    // P2d.Y()).
    let p1inf = is_negative_infinite_value(p1);
    let p2inf = is_positive_infinite_value(p2);
    let mut pt1 = glam::DVec3::ZERO;
    let mut pt2 = glam::DVec3::ZERO;
    if !p1inf {
        let p2d1 = c.point_at(p1);
        pt1 = s.point_at(p2d1.x, p2d1.y);
    }
    if !p2inf {
        let p2d2 = c.point_at(p2);
        pt2 = s.point_at(p2d2.x, p2d2.y);
    }

    let preci = CONFUSION; // BRepLib::Precision().

    // Check for closed curve.
    let mut closed = false;
    if !p1inf && !p2inf {
        closed = pt1.distance(pt2) <= preci;
    }

    // Check if the vertices are on the curve (both vertices null here: the
    // closed branch makes V1 = MakeVertex(P1, preci) and V2 = V1; the open
    // branch makes V1/V2 at the end points).
    let mut v1 = Shape::null();
    let mut v2 = Shape::null();
    if closed {
        v1 = brep.add_tvertex_unique(pt1);
        brep.vertex_mut(v1.clone()).tolerance = preci;
        v2 = v1.clone();
    } else {
        if !p1inf {
            v1 = brep.add_tvertex_unique(pt1);
            brep.vertex_mut(v1.clone()).tolerance = preci;
        }
        if !p2inf {
            v2 = brep.add_tvertex_unique(pt2);
            brep.vertex_mut(v2.clone()).tolerance = preci;
        }
    }

    make_pcurve_edge(brep, &c, s, v1, v2, p1, p2, reverse, closed)
}

/// The edge construction tail of BRepLib_MakeEdge::Init(C, S, ...)
/// (BRepLib_MakeEdge.cxx L1064-1075): B.MakeEdge(E); B.UpdateEdge(E, C, S,
/// TopLoc_Location(), preci); B.Add(E, V1/V2); B.Range(E, p1, p2);
/// E.Orientation(REVERSED) when reversed.
fn make_pcurve_edge(
    brep: &mut BRep,
    c: &Curve2d,
    _s: &Surface3,
    v1: Shape,
    v2: Shape,
    p1: f64,
    p2: f64,
    reverse: bool,
    closed: bool,
) -> Option<BRepBuilderAPIMakeEdge> {
    // B.MakeEdge(E) + B.UpdateEdge(E, C, S, TopLoc_Location(), preci): GAP/
    // architecture — OCCT keys the pcurve by (S, L); the rcad pcurve map is
    // face-keyed, so the surfaceless representation uses the reserved (0, 0)
    // key ("no face"). The W3 ShapeFix_Edge consumer revisits this key.
    let e = brep.add_tedge(None, v1.clone(), v2.clone(), [p1, p2]);
    {
        let ed = brep.edge_mut_inplace(e.clone());
        ed.pcurves
            .insert((0, 0), (c.clone(), p1, p2));
    }
    if closed {
        crate::shhealing::shape_build::brep_tool::set_flag_inplace(
            brep,
            &e,
            rcad_kernel::topo::topods::tshape_flags::CLOSED,
            true,
        );
    }
    let mut e = e;
    if reverse {
        // E.Orientation(TopAbs_REVERSED).
        e.orientation = Orientation::Reversed;
    }
    Some(BRepBuilderAPIMakeEdge { edge: e })
}

/// OCCT ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)
/// (ElCLib.cxx L115-150).
fn elclib_adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    // OCCT ElCLib.cxx L121: Precision::IsInfinite(UFirst) || Precision::IsInfinite(ULast)
    // (Precision.hxx L350-353).
    if rcad_kernel::precision::is_infinite_value(u_first)
        || rcad_kernel::precision::is_infinite_value(u_last)
    {
        return;
    }

    let a_period = u_last - u_first;

    if a_period < (f64::EPSILON * u_last.abs()) {
        // In order to avoid FLT_Overflow exception (test bugs moddata_1
        // bug22757).
        return;
    }

    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}
