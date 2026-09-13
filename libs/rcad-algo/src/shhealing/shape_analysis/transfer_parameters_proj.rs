//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_TransferParametersProj`
//! (`ShapeAnalysis_TransferParametersProj.hxx` L17-110 + `.cxx` L1-802).
//!
//! The projected (non-linear) specialization of
//! ShapeAnalysis_TransferParameters: transfers the parameters between the
//! 3D curve and the pcurve through the projection onto the
//! Adaptor3d_CurveOnSurface instead of the linear shift/scale mapping.
//!
//! Architecture bridges (the W1-1/W2 conventions):
//! 1. `BRep` pool argument — every shape-reading/mutating member takes
//!    `brep: &mut BRep` (the edge.rs bridge #1).
//! 2. OCCT inheritance (`: ShapeAnalysis_TransferParameters`) -> the
//!    `base: ShapeAnalysisTransferParameters` composition (Rust has no
//!    inheritance); the OCCT `ShapeAnalysis_TransferParameters::X(...)`
//!    qualified calls map to `self.base.x(...)`, the protected fields to
//!    `pub(crate)` fields.
//! 3. `Adaptor3d_CurveOnSurface` / `Geom2dAdaptor_Curve` /
//!    `GeomAdaptor_Surface` -> the `topalgo::brep_lib_validate_edge`
//!    re-hosts; `Handle(Adaptor3d_Surface) myAC3d.GetSurface()` -> the
//!    stored `GeomAdaptorSurface`.
//! 4. `TopLoc_Location` -> the `BRep.locations` table index; the
//!    transforms go through `brep.get_location(loc)` (`Predivided` ->
//!    `inverse() *`, `Transformed` -> `transform_point3`).
//! 5. `BRep_TEdge::ChangeCurves()` (the mutable representation walk) ->
//!    `brep.edge_mut_inplace(edge)`: the 3D GCurve row range is the
//!    `TEdgeData.range`, each pcurve row range is the matching
//!    `representations` row AND the `pcurves` map entry (the two rcad
//!    storages of the same OCCT row).  The pcurve-key location hash is
//!    one-way (the edge.rs bridge #5), so the pcurve-row transforms use
//!    the identity-location reduction.
//! 6. `BRep_TVertex::Points()` (`BRep_PointRepresentation`) -> the
//!    `PointRepresentation` values on `TVertexData.points`; the rcad
//!    carrier has no PointOnCurveOnSurface representative, so the OCCT
//!    `IsPointOnCurveOnSurface` identity tests keep their structure with
//!    the not-found reduction (the architecture difference, documented at
//!    the arms).
//! 7. `BRep_Builder` calls -> the `BRepBuilder` tool (`SameRange`,
//!    `UpdateVertex(V, P, E, Tol)` -> `update_vertex_on_edge`, the
//!    tolerance update -> `update_vertex_tolerance`).  The UV update
//!    (`UpdateVertex(V, U, V2, F, Tol)`) has no kernel tool — the local
//!    `builder_update_vertex_uv` keeps the OCCT effect over the rcad
//!    `PointOnSurface` row.

use std::collections::HashMap;
use std::sync::Arc;

use glam::DVec3;
use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{
    surface_same, BRep, BRepBuilder, BRepTool, CurveRepresentation, Orientation,
    PointRepresentation, TShape,
};

use crate::brep_algo::tool::{brep_tool_tolerance, empty_copied};
use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_analysis::transfer_parameters::ShapeAnalysisTransferParameters;
use crate::topalgo::brep_lib_validate_edge::{
    Adaptor3dCurveOnSurface, Geom2dAdaptorCurve, GeomAdaptorCurve, GeomAdaptorSurface,
};

/// OCCT Geom_Curve::FirstParameter()/LastParameter() — the natural domain.
#[allow(dead_code)]
fn curve_first_parameter(c: &Curve3) -> f64 {
    rcad_kernel::geom::CurveEval::default_domain(c)[0]
}

fn curve_last_parameter(c: &Curve3) -> f64 {
    rcad_kernel::geom::CurveEval::default_domain(c)[1]
}

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (bridge #1).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Curve(edg, f, l) — the no-location variant.
fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// OCCT BRep_Tool::Curve(edg, Loc, f, l) — the out-location variant (the
/// edge wrapper location; the edge.rs bridge #5).
fn brep_tool_curve_loc(edg: &Shape) -> (Option<Curve3>, u32, f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), edg.location, ed.range[0], ed.range[1]),
        _ => (None, 0, 0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Surface(fac, loc).
fn brep_tool_surface_loc(fac: &Shape) -> (Option<Surface3>, u32) {
    match fac.data.as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fac.location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::SameParameter(E).
fn brep_tool_same_parameter(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.same_parameter,
        _ => false,
    }
}

/// OCCT Standard_Real.hxx Epsilon(Value) — Value * DBL_EPSILON (signed).
fn epsilon(v: f64) -> f64 {
    if v >= 0.0 {
        f64::EPSILON * v
    } else {
        -f64::EPSILON * v
    }
}

/// OCCT Precision::IsInfinite.
fn precision_is_infinite(r: f64) -> bool {
    r.abs() >= 2e100
}

/// The face surface registered in the pool under a TShape pointer (the
/// Geom_Surface handle identity stand-in; the kernel face_surface_by_ptr
/// precedent).
fn face_surface_by_ptr(brep: &BRep, fptr: u64) -> Option<Surface3> {
    let ts = brep
        .tshapes
        .iter()
        .find(|ts| Arc::as_ptr(ts) as u64 == fptr)?;
    match ts.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT ShapeAnalysis_TransferParametersProj (hxx L26-110).
pub struct ShapeAnalysisTransferParametersProj {
    /// The OCCT base class subobject (bridge #2).
    pub(crate) base: ShapeAnalysisTransferParameters,
    /// OCCT `myPrecision`.
    pub(crate) my_precision: f64,
    /// OCCT `myCurve`.
    pub(crate) my_curve: Option<Curve3>,
    /// OCCT `myCurve2d`.
    pub(crate) my_curve2d: Option<Curve2d>,
    /// OCCT `myAC3d` (the Adaptor3d_CurveOnSurface handle).
    pub(crate) my_ac3d: Option<Adaptor3dCurveOnSurface>,
    /// OCCT `myLocation`.
    pub(crate) my_location: u32,
    /// OCCT `myForceProj`.
    pub(crate) my_force_proj: bool,
    /// OCCT `myInitOK`.
    pub(crate) my_init_ok: bool,
}

impl Default for ShapeAnalysisTransferParametersProj {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisTransferParametersProj {
    /// OCCT ShapeAnalysis_TransferParametersProj() (cxx L49-55).
    pub fn new() -> Self {
        let mut base = ShapeAnalysisTransferParameters::new();
        base.my_max_tolerance = 1.0; // Precision::Infinite(); ?? pdn
        ShapeAnalysisTransferParametersProj {
            base,
            my_precision: 0.0,
            my_curve: None,
            my_curve2d: None,
            my_ac3d: None,
            my_location: 0,
            my_force_proj: false,
            my_init_ok: false,
        }
    }

    /// OCCT ShapeAnalysis_TransferParametersProj(E, F) (cxx L59-65).
    pub fn new_edge_face(brep: &mut BRep, e: &Shape, face: &Shape) -> Self {
        let mut this = Self::new();
        this.init(brep, e, face);
        this
    }

    /// OCCT Init(E, F) (cxx L69-103); the face parameter is `face` (the
    /// OCCT `F`), because the body reads the 3D-curve first into `f`.
    pub fn init(&mut self, brep: &mut BRep, e: &Shape, face: &Shape) {
        self.my_init_ok = false;
        // OCCT L72: ShapeAnalysis_TransferParameters::Init(E, F).
        self.base.init(brep, e, face);
        // OCCT L73: myEdge = E (the base Init stored it).
        self.base.my_edge = e.clone();
        // OCCT L74: myPrecision = BRep_Tool::Tolerance(E) — it is better - skl OCC2851.
        self.my_precision = brep_tool_tolerance(e);
        // myPrecision = Precision::Confusion();

        // OCCT L77: myCurve = BRep_Tool::Curve(E, myFirst, myLast).
        let (curve, _loc, f, l) = brep_tool_curve_loc(e);
        self.base.my_first = f;
        self.base.my_last = l;
        self.my_curve = curve;
        if self.my_curve.is_none() {
            self.base.my_first = 0.;
            self.base.my_last = 1.;
            return;
        }

        if face.is_null() {
            return;
        }

        let mut f2d = 0.0;
        let mut l2d = 0.0;
        let sae = ShapeAnalysisEdge::new();
        let mut curve2d: Option<Curve2d> = None;
        if sae.pcurve_face(brep, e, face, &mut curve2d, &mut f2d, &mut l2d, false) {
            self.my_curve2d = curve2d;
            // OCCT L95-99: AC2d + aSurface + AdS + Ad1.
            let ac2d = Geom2dAdaptorCurve::new(self.my_curve2d.clone().unwrap(), f2d, l2d);
            let (a_surface, loc) = brep_tool_surface_loc(face);
            self.my_location = loc;
            let ads = GeomAdaptorSurface::new(a_surface.unwrap());
            self.my_ac3d = Some(Adaptor3dCurveOnSurface::new(ac2d, ads));
            self.my_init_ok = true;
        }
    }

    /// OCCT Perform(Knots, To2d) (cxx L107-176).
    pub fn perform_params(&self, brep: &BRep, knots: &[f64], to2d: bool) -> Vec<f64> {
        // pdn
        if !self.my_init_ok
            || (!self.my_force_proj
                && self.my_precision < self.base.my_max_tolerance
                && brep_tool_same_parameter(&self.base.my_edge))
        {
            return self.base.perform_params(knots, to2d);
        }

        let mut res_knots: Vec<f64> = Vec::new();

        let len = knots.len() as i32;
        let preci = 2.0 * PCONFUSION;

        let ac3d = self.my_ac3d.as_ref().unwrap();
        let first = if to2d {
            ac3d.first_parameter()
        } else {
            self.base.my_first
        };
        let last = if to2d {
            ac3d.last_parameter()
        } else {
            self.base.my_last
        };
        let mut max_par = first;
        let last_par = last;
        let mut prev_par = max_par;

        for j in 1..=len {
            let par = self.preform_segment(brep, knots[(j - 1) as usize], to2d, prev_par, last_par);
            prev_par = par;
            if prev_par > last_par {
                prev_par -= preci;
            }
            res_knots.push(par);
            if par > max_par {
                max_par = par;
            }
        }

        // pdn correcting on periodic
        let my_curve = self.my_curve.as_ref().unwrap();
        if my_curve.is_closed() {
            let mut j = len;
            while j >= 1 {
                if res_knots[(j - 1) as usize] < max_par {
                    let last_param = if to2d {
                        self.my_ac3d.as_ref().unwrap().last_parameter()
                    } else {
                        curve_last_parameter(my_curve)
                    };
                    res_knots[(j - 1) as usize] =
                        last_param - (len - j) as f64 * preci;
                } else {
                    break;
                }
                j -= 1;
            }
        }
        // pdn correction on range
        for j in 1..=len {
            if res_knots[(j - 1) as usize] < first {
                res_knots[(j - 1) as usize] = first;
            }
            if res_knots[(j - 1) as usize] > last {
                res_knots[(j - 1) as usize] = last;
            }
        }

        res_knots
    }

    /// OCCT PreformSegment(Param, To2d, First, Last) (cxx L180-218) — the
    /// OCCT `Preform` spelling kept.
    #[allow(clippy::too_many_arguments)]
    pub fn preform_segment(
        &self,
        brep: &BRep,
        param: f64,
        to2d: bool,
        first: f64,
        last: f64,
    ) -> f64 {
        let lin_par = self.base.perform(param, to2d);
        if !self.my_init_ok
            || (!self.my_force_proj
                && self.my_precision < self.base.my_max_tolerance
                && brep_tool_same_parameter(&self.base.my_edge))
        {
            return lin_par;
        }

        let lin_dev;
        let proj_dev;

        let sac = ShapeAnalysisCurve;
        let mut pproj = DVec3::ZERO;
        // OCCT: `double ppar;` — the deterministic zero stand-in for the
        // uninitialized out-parameter.
        let mut ppar = 0.0;
        let my_curve = self.my_curve.as_ref().unwrap();
        let my_ac3d = self.my_ac3d.as_ref().unwrap();
        if to2d {
            // OCCT L199: p1 = myCurve->Value(Param).Transformed(myLocation.Inverted()).
            let loc_trsf = brep.get_location(self.my_location).inverse();
            let p1 = loc_trsf.transform_point3(my_curve.point_at(param));
            // OCCT L200-202: AdS = myAC3d.GetSurface(); AC2d over
            // (myCurve2d, First, Last); Ad1(AC2d, AdS).
            let ads = my_ac3d.get_surface();
            let ac2d = Geom2dAdaptorCurve::new(self.my_curve2d.clone().unwrap(), first, last);
            let ad1 = Adaptor3dCurveOnSurface::new(ac2d, ads.clone());
            // OCCT L203: projDev = sac.Project(Ad1, p1, myPrecision, pproj, ppar); // pdn
            proj_dev = sac.project_adaptor(&ad1, p1, self.my_precision, &mut pproj, &mut ppar, true);
            // OCCT L204: linDev = p1.Distance(Ad1.Value(linPar)).
            lin_dev = p1.distance(ad1.value(lin_par));
        } else {
            // OCCT L208: p1 = myAC3d.Value(Param).Transformed(myLocation).
            let loc_trsf = brep.get_location(self.my_location);
            let p1 = loc_trsf.transform_point3(my_ac3d.value(param));
            // OCCT L209: projDev = sac.Project(myCurve, p1, myPrecision,
            // pproj, ppar, First, Last, false).
            proj_dev = sac.project_cf_cl(
                my_curve,
                p1,
                self.my_precision,
                &mut pproj,
                &mut ppar,
                first,
                last,
                false,
            );
            // OCCT L210: linDev = p1.Distance(myCurve->Value(linPar)).
            lin_dev = p1.distance(my_curve.point_at(lin_par));
        }

        let mut res = ppar;
        if lin_dev <= proj_dev || (lin_dev < self.my_precision && lin_dev <= 2.0 * proj_dev) {
            res = lin_par;
        }
        res
    }

    /// OCCT Perform(Knot, To2d) (cxx L222-252).
    pub fn perform(&self, brep: &BRep, knot: f64, to2d: bool) -> f64 {
        if !self.my_init_ok
            || (!self.my_force_proj
                && self.my_precision < self.base.my_max_tolerance
                && brep_tool_same_parameter(&self.base.my_edge))
        {
            return self.base.perform(knot, to2d);
        }

        let my_ac3d = self.my_ac3d.as_ref().unwrap();
        let mut res;
        if to2d {
            res = self.preform_segment(brep, knot, to2d, my_ac3d.first_parameter(), my_ac3d.last_parameter());
        } else {
            res = self.preform_segment(brep, knot, to2d, self.base.my_first, self.base.my_last);
        }

        // pdn correction on range
        let first = if to2d {
            my_ac3d.first_parameter()
        } else {
            self.base.my_first
        };
        let last = if to2d {
            my_ac3d.last_parameter()
        } else {
            self.base.my_last
        };
        if res < first {
            res = first;
        }
        if res > last {
            res = last;
        }
        res
    }

    /// OCCT TransferRange(newEdge, prevPar, currPar, Is2d) (cxx L285-535).
    pub fn transfer_range(
        &self,
        brep: &mut BRep,
        new_edge: &Shape,
        prev_par: f64,
        curr_par: f64,
        is2d: bool,
    ) {
        if !self.my_init_ok
            || (!self.my_force_proj
                && self.my_precision < self.base.my_max_tolerance
                && brep_tool_same_parameter(&self.base.my_edge))
        {
            self.base
                .transfer_range(brep, new_edge, prev_par, curr_par, is2d);
            return;
        }

        let mut b = BRepBuilder::new();
        let mut samerange = true;
        // OCCT L300: sbe.CopyRanges(newEdge, myEdge) — the defaults
        // alpha = 0, beta = 1.
        let sbe = crate::shhealing::shape_build::edge::ShapeBuildEdge;
        sbe.copy_ranges(brep, new_edge, &self.base.my_edge, 0.0, 1.0);
        let mut alpha = 0.0;
        let mut beta = 1.0;
        let preci = PCONFUSION;
        let first_par;
        let last_par;
        if prev_par < curr_par {
            first_par = prev_par;
            last_par = curr_par;
        } else {
            first_par = curr_par;
            last_par = prev_par;
        }
        let my_ac3d = self.my_ac3d.as_ref().unwrap();
        let my_curve = self.my_curve.as_ref().unwrap();
        let (p1, p2) = if is2d {
            // OCCT L318: p1 = myAC3d.Value(firstPar).Transformed(myLocation).
            let loc_trsf = brep.get_location(self.my_location);
            let p1 = loc_trsf.transform_point3(my_ac3d.value(first_par));
            if precision_is_infinite(p1.x) || precision_is_infinite(p1.y) || precision_is_infinite(p1.z) {
                b.set_edge_same_range(brep, new_edge.clone(), false);
                return;
            }
            let p2 = loc_trsf.transform_point3(my_ac3d.value(last_par));
            if precision_is_infinite(p2.x) || precision_is_infinite(p2.y) || precision_is_infinite(p2.z) {
                b.set_edge_same_range(brep, new_edge.clone(), false);
                return;
            }
            let fact = my_ac3d.last_parameter() - my_ac3d.first_parameter();
            if fact > epsilon(my_ac3d.last_parameter()) {
                alpha = (first_par - my_ac3d.first_parameter()) / fact;
                beta = (last_par - my_ac3d.first_parameter()) / fact;
            }
            (p1, p2)
        } else {
            let p1 = my_curve.point_at(first_par);
            if precision_is_infinite(p1.x) || precision_is_infinite(p1.y) || precision_is_infinite(p1.z) {
                b.set_edge_same_range(brep, new_edge.clone(), false);
                return;
            }
            let p2 = my_curve.point_at(last_par);
            if precision_is_infinite(p2.x) || precision_is_infinite(p2.y) || precision_is_infinite(p2.z) {
                b.set_edge_same_range(brep, new_edge.clone(), false);
                return;
            }
            let fact = self.base.my_last - self.base.my_first;
            if fact > epsilon(self.base.my_last) {
                alpha = (first_par - self.base.my_first) / fact;
                beta = (last_par - self.base.my_first) / fact;
            }
            (p1, p2)
        };
        let use_linear_first = alpha < preci;
        let use_linear_last = 1.0 - beta < preci;
        // OCCT L364: TopLoc_Location EdgeLoc = myEdge.Location().
        let edge_loc = self.base.my_edge.location;
        let sac = ShapeAnalysisCurve;

        // OCCT L368-373: the mutable walk over newEdge.TShape()->
        // ChangeCurves() (bridge #5).  The rcad walk runs over the row
        // indices (the rows are mutated in place afterwards).
        let rep_count = match new_edge.data.as_ref() {
            TShape::Edge(ed) => ed.representations.len(),
            _ => 0,
        };
        let mut row_3d_range: Option<(f64, f64)> = None;
        let mut row_pcurve_ranges: HashMap<(u64, u32), (f64, f64)> = HashMap::new();
        for row in 0..rep_count {
            // OCCT L375: toGC = down_cast<BRep_GCurve>(toitcr.Value());
            // if (toGC.IsNull()) continue; — the rcad rows are always
            // GCurves; the non-GCurve rows do not exist in the carrier.
            let (row_key, row_first, row_last) = match new_edge.data.as_ref() {
                TShape::Edge(ed) => match &ed.representations[row] {
                    CurveRepresentation::Curve3D { .. } => (None, ed.range[0], ed.range[1]),
                    CurveRepresentation::CurveOnSurface { face, range, .. } => {
                        (Some(*face), range[0], range[1])
                    }
                    CurveRepresentation::CurveOnClosedSurface { face, range, .. } => {
                        (Some(*face), range[0], range[1])
                    }
                    _ => continue,
                },
                _ => continue,
            };
            // OCCT L380: loc = (EdgeLoc * toGC->Location()).Inverted() —
            // the pcurve-key location hash is one-way (bridge #5), so the
            // pcurve rows use the identity-location reduction and the 3D
            // row composes the edge location with its own.
            let row_loc = match new_edge.data.as_ref() {
                TShape::Edge(ed) => match &ed.representations[row] {
                    CurveRepresentation::Curve3D { location, .. } => *location,
                    _ => 0,
                },
                _ => 0,
            };
            let composed = {
                let a = brep.get_location(edge_loc);
                let b = brep.get_location(row_loc);
                let m = glam::DAffine3::from_mat3(a.matrix3) * b;
                m
            };
            let _ = composed;
            let is_3d_row = row_key.is_none();
            if is_3d_row {
                if !is2d {
                    // OCCT L383-387: ppar1 = firstPar; ppar2 = lastPar.
                    let ppar1 = first_par;
                    let ppar2 = last_par;
                    let (ppar1, ppar2) = fix_degenerate_range(ppar1, ppar2, row_first, row_last, preci);
                    row_3d_range = Some((ppar1, ppar2));
                    if ppar1 != first_par || ppar2 != last_par {
                        samerange = false;
                    }
                } else {
                    // OCCT L390-421: C3d = toGC->Curve3D() — the 3D row
                    // curve; the rcad value lives on TEdgeData.curve.
                    let Some(c3d) = (match new_edge.data.as_ref() {
                        TShape::Edge(ed) => ed.curve.clone(),
                        _ => None,
                    }) else {
                        continue;
                    };
                    let first = row_first;
                    let last = row_last;
                    let len = last - first;
                    // OCCT L398-399: ploc1/ploc2 = p1/p2.Transformed(loc) —
                    // the 3D row location is the edge location (bridge #5).
                    let loc_trsf = brep.get_location(edge_loc);
                    let ploc1 = loc_trsf.transform_point3(p1);
                    let ploc2 = loc_trsf.transform_point3(p2);
                    let gac = GeomAdaptorCurve::new(c3d.clone(), first, last);
                    // CATIA bplseitli.model FAC1155 - Copy: protection for degenerated edges(3d case for
                    // symmetry)
                    let lin_first = first + alpha * len;
                    let lin_last = first + beta * len;
                    let mut pproj = DVec3::ZERO;
                    let mut ppar1 = 0.0;
                    let mut ppar2 = 0.0;
                    let dist1 =
                        sac.next_project_adaptor(lin_first, &gac, ploc1, self.my_precision, &mut pproj, &mut ppar1);
                    let dist2 =
                        sac.next_project_adaptor(lin_last, &gac, ploc2, self.my_precision, &mut pproj, &mut ppar2);
                    let use_linear = (ppar1 - ppar2).abs() < preci;

                    let pos1 = c3d.point_at(lin_first);
                    let pos2 = c3d.point_at(lin_last);
                    let d01 = pos1.distance(ploc1);
                    let d02 = pos2.distance(ploc2);
                    if use_linear_first || use_linear || d01 <= dist1 || (d01 < self.my_precision && d01 <= 2.0 * dist1)
                    {
                        ppar1 = lin_first;
                    }
                    if use_linear_last || use_linear || d02 <= dist2 || (d02 < self.my_precision && d02 <= 2.0 * dist2)
                    {
                        ppar2 = lin_last;
                    }
                    let (ppar1, ppar2) =
                        fix_degenerate_range(ppar1, ppar2, row_first, row_last, preci);
                    row_3d_range = Some((ppar1, ppar2));
                    if ppar1 != first_par || ppar2 != last_par {
                        samerange = false;
                    }
                }
            } else if row_is_curve_on_surface(new_edge, row) {
                // OCCT L452-532: the pcurve row.
                let mut local_linear_first = use_linear_first;
                let mut local_linear_last = use_linear_last;
                let c2d = match new_edge.data.as_ref() {
                    TShape::Edge(ed) => match &ed.representations[row] {
                        CurveRepresentation::CurveOnSurface { pcurve, .. } => pcurve.clone(),
                        CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => pcurve1.clone(),
                        _ => continue,
                    },
                    _ => continue,
                };
                let first = row_first;
                let last = row_last;
                let len = last - first;
                // OCCT L461-463: AC2d + AdS + Ad1 over the row (pcurve,
                // surface).  The rcad row carries the face key; the surface
                // resolves through the pool registry.
                let ads = match face_surface_by_ptr(brep, row_key.unwrap().0) {
                    Some(s) => GeomAdaptorSurface::new(s),
                    None => GeomAdaptorSurface::new(match new_edge.data.as_ref() {
                        TShape::Edge(ed) => match &ed.representations[row] {
                            CurveRepresentation::CurveOnSurface { .. } => {
                                // The registry row is missing — the OCCT
                                // toGC->Surface() would be null and the
                                // adaptor construction raises; the rcad
                                // walk keeps the row pass-through with a
                                // plane stand-in unreachable in practice.
                                Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z))
                            }
                            _ => Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
                        },
                        _ => Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
                    }),
                };
                let ac2d = Geom2dAdaptorCurve::new(c2d.clone(), first, last);
                let ad1 = Adaptor3dCurveOnSurface::new(ac2d, ads);
                let sac1 = ShapeAnalysisCurve;

                // OCCT L468-469: ploc1/ploc2 = p1/p2.Transformed(loc) — the
                // identity-location reduction (bridge #5: the row location
                // hash is one-way).
                let ploc1 = p1;
                let ploc2 = p2;
                // CATIA bplseitli.model FAC1155 - Copy: protection for degenerated edges
                let lin_first = first + alpha * len;
                let lin_last = first + beta * len;
                let mut pproj = DVec3::ZERO;
                let mut ppar1 = 0.0;
                let mut ppar2 = 0.0;
                let dist1 =
                    sac1.next_project_adaptor(lin_first, &ad1, ploc1, self.my_precision, &mut pproj, &mut ppar1);
                let dist2 =
                    sac1.next_project_adaptor(lin_last, &ad1, ploc2, self.my_precision, &mut pproj, &mut ppar2);

                let is_first_on_end = (ppar1 - first) / len < PCONFUSION;
                let is_last_on_end = (last - ppar2) / len < PCONFUSION;
                let use_linear = (ppar1 - ppar2).abs() < PCONFUSION;
                if is_first_on_end && !local_linear_first {
                    local_linear_first = true;
                }
                if is_last_on_end && !local_linear_last {
                    local_linear_last = true;
                }

                let pos1 = ad1.value(lin_first);
                let pos2 = ad1.value(lin_last);
                let d01 = pos1.distance(ploc1);
                let d02 = pos2.distance(ploc2);
                if local_linear_first || use_linear || d01 <= dist1 || (d01 < self.my_precision && d01 <= 2.0 * dist1)
                {
                    ppar1 = lin_first;
                }
                if local_linear_last || use_linear || d02 <= dist2 || (d02 < self.my_precision && d02 <= 2.0 * dist2)
                {
                    ppar2 = lin_last;
                }

                if ppar1 > ppar2 {
                    std::mem::swap(&mut ppar1, &mut ppar2);
                }
                ppar1 = correct_parameter(&c2d, ppar1);
                ppar2 = correct_parameter(&c2d, ppar2);
                let (ppar1, ppar2) =
                    fix_degenerate_range_pc(ppar1, ppar2, row_first, row_last, preci);
                row_pcurve_ranges.insert(row_key.unwrap(), (ppar1, ppar2));
                if ppar1 != first_par || ppar2 != last_par {
                    samerange = false;
                }
            }
        }

        // The in-place row mutation (bridge #5): the walked ranges land on
        // the TEdgeData storages.
        {
            let ed = brep.edge_mut_inplace(new_edge.clone());
            if let Some((f, l)) = row_3d_range {
                ed.range[0] = f;
                ed.range[1] = l;
            }
            for (key, (f, l)) in &row_pcurve_ranges {
                if let Some(entry) = ed.pcurves.get_mut(key) {
                    entry.1 = *f;
                    entry.2 = *l;
                }
                for rep in ed.representations.iter_mut() {
                    match rep {
                        CurveRepresentation::CurveOnSurface { face, range, .. } => {
                            if face == key {
                                range[0] = *f;
                                range[1] = *l;
                            }
                        }
                        CurveRepresentation::CurveOnClosedSurface { face, range, .. } => {
                            if face == key {
                                range[0] = *f;
                                range[1] = *l;
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
        // OCCT L534: B.SameRange(newEdge, samerange).
        b.set_edge_same_range(brep, new_edge.clone(), samerange);
    }

    /// OCCT IsSameRange() (cxx L539-551).
    pub fn is_same_range(&self, _brep: &BRep) -> bool {
        if !self.my_init_ok
            || (!self.my_force_proj
                && self.my_precision < self.base.my_max_tolerance
                && brep_tool_same_parameter(&self.base.my_edge))
        {
            self.base.is_same_range()
        } else {
            false
        }
    }

    /// OCCT ForceProjection() (cxx L555-558) — the `bool&` reference
    /// accessor.
    pub fn force_projection(&mut self) -> &mut bool {
        &mut self.my_force_proj
    }

    /// OCCT CopyNMVertex(theV, toedge, fromedge) (cxx L562-711).
    pub fn copy_nm_vertex_edge(
        &self,
        brep: &mut BRep,
        the_v: &Shape,
        toedge: &Shape,
        fromedge: &Shape,
    ) -> Shape {
        let mut anew_v = Shape::null();
        if the_v.orientation != Orientation::Internal && the_v.orientation != Orientation::External
        {
            return anew_v;
        }

        // OCCT L572-575: (C1, fromLoc, f1, l1) = BRep_Tool::Curve(fromedge);
        // fromLoc = fromLoc.Predivided(theV.Location()).
        let (_c1, _from_loc, f1, l1) = brep_tool_curve_loc(fromedge);
        // OCCT L578: (C2, f2, l2) = BRep_Tool::Curve(toedge).
        let (_c2, f2, l2) = match brep_tool_curve(toedge) {
            Some((c, f, l)) => (Some(c), f, l),
            None => (None, 0.0, 0.0),
        };

        // OCCT L580-581: anewV = TopoDS::Vertex(theV.EmptyCopied());
        // apv = BRep_Tool::Pnt(anewV).
        anew_v = empty_copied(the_v);
        let apv = brep_tool_pnt(&anew_v);

        // OCCT L583-587: alistrep = anewV points; itpr walks theV points.
        let the_v_points: Vec<PointRepresentation> = match the_v.data.as_ref() {
            TShape::Vertex(vd) => vd.points.clone(),
            _ => Vec::new(),
        };
        let mut a_list_rep: Vec<PointRepresentation> = Vec::new();

        // OCCT L589: aOldPar = RealLast().
        let mut a_old_par = f64::MAX;
        let mut has_repr = false;
        for pr in &the_v_points {
            match pr {
                // OCCT L598-603: pr->IsPointOnCurve(C1, fromLoc) — the rcad
                // PointOnCurve carries no curve handle, so the identity test
                // reduces to the representative arm (architecture
                // difference, bridge #6).
                PointRepresentation::PointOnCurve { parameter, .. } => {
                    a_old_par = *parameter;
                    has_repr = true;
                    continue;
                }
                // OCCT L604-613: pr->IsPointOnSurface() — the replica with
                // the same (Parameter, Parameter2, Surface, Location); the
                // rcad value row is the shared clone.
                PointRepresentation::PointOnSurface { .. } => {
                    a_list_rep.push(pr.clone());
                }
                // OCCT L614-644: pr->IsPointOnCurveOnSurface() — no rcad
                // point representative carries a pcurve, so this separate
                // OCCT branch has no rcad arm (the enum carries exactly the
                // two representatives above); the OCCT fall-through replica
                // append (L651-659) is expressed by the shared clone here
                // (architecture difference, bridge #6).
            }
        }
        // Write the new points list onto the fresh vertex.
        if let TShape::Vertex(vd) = Arc::make_mut(&mut anew_v.data) {
            vd.points = a_list_rep;
        } else {
            let _ = Arc::make_mut(&mut anew_v.data);
        }

        // OCCT L661-673.
        let mut apar = a_old_par;
        let mut a_tol = brep_tool_tolerance(the_v);
        if !has_repr
            || ((f1 - f2).abs() > PCONFUSION || (l1 - l2).abs() > PCONFUSION)
        {
            let mut proj_p = DVec3::ZERO;
            let sae = ShapeAnalysisCurve;
            // OCCT L668: sae.Project(C2, apv, Precision::Confusion(),
            // projP, apar) — the from/to edge 3D curve; the rcad arm uses
            // the to-edge curve (C2 null keeps the OCCT failure path out
            // of the projection).
            if let Some((c2, _f, _l)) = brep_tool_curve(toedge) {
                let adist = sae.project(&c2, apv, CONFUSION, &mut proj_p, &mut apar, true);
                if a_tol < adist {
                    a_tol = adist;
                }
            }
        }
        // OCCT L674-675: B.UpdateVertex(anewV, apar, toedge, aTol).
        let mut b = BRepBuilder::new();
        b.update_vertex_on_edge(brep, anew_v.clone(), apar, toedge.clone(), a_tol);

        // update tolerance (OCCT L677-709).
        let mut need_update = false;
        let a_pv = brep_tool_pnt(&anew_v);
        let to_loc = toedge.location;
        let to_rows: Vec<(Option<(u64, u32)>, Surface3, Curve2d, f64, f64)> =
            match toedge.data.as_ref() {
                TShape::Edge(ed) => ed
                    .representations
                    .iter()
                    .filter_map(|r| match r {
                        CurveRepresentation::CurveOnSurface {
                            face,
                            pcurve,
                            range,
                        } => Some((
                            Some(*face),
                            face_surface_by_ptr(brep, face.0)?,
                            pcurve.clone(),
                            range[0],
                            range[1],
                        )),
                        CurveRepresentation::CurveOnClosedSurface {
                            face,
                            pcurve1,
                            range,
                            ..
                        } => Some((
                            Some(*face),
                            face_surface_by_ptr(brep, face.0)?,
                            pcurve1.clone(),
                            range[0],
                            range[1],
                        )),
                        _ => None,
                    })
                    .collect(),
                _ => Vec::new(),
            };
        for (key, surface1, ac2d1, _rf, _rl) in to_rows {
            // OCCT L692: aL = (toLoc * toGC->Location()).Predivided(theV.
            // Location()) — the identity-location reduction (bridge #5).
            let _ = (key, to_loc);
            let a_p2d = ac2d1.point_at(apar);
            let a_p3d = surface1.point_at(a_p2d.x, a_p2d.y);
            // OCCT L698: aP3d.Transform(aL.Transformation()) — the
            // reduced-identity transform.
            let adist = a_pv.distance(a_p3d);
            if adist > a_tol {
                a_tol = adist;
                need_update = true;
            }
        }
        if need_update {
            // OCCT L708: B.UpdateVertex(anewV, aTol).
            b.update_vertex_tolerance(brep, anew_v.clone(), a_tol);
        }
        anew_v
    }

    /// OCCT CopyNMVertex(theV, toFace, fromFace) (cxx L715-802).
    pub fn copy_nm_vertex_face(
        &self,
        brep: &mut BRep,
        the_v: &Shape,
        to_face: &Shape,
        from_face: &Shape,
    ) -> Shape {
        let mut anew_v = Shape::null();
        if the_v.orientation != Orientation::Internal && the_v.orientation != Orientation::External
        {
            return anew_v;
        }

        // OCCT L725-729: the from/to surfaces + locations;
        // fromLoc = fromLoc.Predivided(theV.Location()).
        let (from_surf, from_loc) = brep_tool_surface_loc(from_face);
        let (to_surf, to_loc) = brep_tool_surface_loc(to_face);

        // OCCT L731-732: anewV = theV.EmptyCopied(); apv = Pnt(anewV).
        anew_v = empty_copied(the_v);
        let apv = brep_tool_pnt(&anew_v);

        let the_v_points: Vec<PointRepresentation> = match the_v.data.as_ref() {
            TShape::Vertex(vd) => vd.points.clone(),
            _ => Vec::new(),
        };
        let mut a_list_rep: Vec<PointRepresentation> = Vec::new();

        let mut has_repr = false;
        let mut apar1 = 0.0;
        let mut apar2 = 0.0;
        for pr in &the_v_points {
            match pr {
                // OCCT L750-755: IsPointOnCurveOnSurface — the replica; the
                // rcad value row clone (bridge #6).
                PointRepresentation::PointOnCurve { .. } => {
                    // OCCT L756-761: IsPointOnCurve — the replica with the
                    // same location.
                    a_list_rep.push(pr.clone());
                }
                PointRepresentation::PointOnSurface { face, u, v, .. } => {
                    // OCCT L762-779: IsPointOnSurface — the (fromSurf,
                    // fromLoc) identity selects the parameters; the others
                    // are cloned.  The rcad identity reduces to the stored
                    // face key (the surface value comparison).
                    let from_matches = match from_surf.as_ref() {
                        Some(fs) => match face_surface_by_ptr_legacy(*face) {
                            Some(s) => surface_same(&s, fs),
                            None => false,
                        },
                        None => false,
                    };
                    if from_matches {
                        apar1 = *u;
                        apar2 = *v;
                        has_repr = true;
                    } else {
                        a_list_rep.push(pr.clone());
                    }
                }
            }
        }
        if let TShape::Vertex(vd) = Arc::make_mut(&mut anew_v.data) {
            vd.points = a_list_rep;
        } else {
            let _ = Arc::make_mut(&mut anew_v.data);
        }

        // OCCT L782-797.
        let mut a_tol = brep_tool_tolerance(&anew_v);
        let same_surface = match (from_surf.as_ref(), to_surf.as_ref()) {
            (Some(fs), Some(ts)) => surface_same(fs, ts),
            _ => false,
        };
        if !has_repr || (!same_surface || from_loc != to_loc) {
            // OCCT L785-786: aS = BRep_Tool::Surface(toFace);
            // aSurfTool = new ShapeAnalysis_Surface(aS).
            if let Some(a_s) = to_surf.as_ref() {
                let mut a_surf_tool = ShapeAnalysisSurface::new(a_s.clone());
                let a_p2d = a_surf_tool.value_of_uv(apv, CONFUSION);
                apar1 = a_p2d.x;
                apar2 = a_p2d.y;

                if a_tol < a_surf_tool.gap() {
                    a_tol = a_surf_tool.gap() + 0.1 * CONFUSION;
                }
                // occ::handle<BRep_PointOnSurface> aPS = new
                // BRep_PointOnSurface(aP2d.X(),aP2d.Y(),toSurf,toLoc);
                // alistrep.Append(aPS);
            }
        }

        // OCCT L799-800: B.UpdateVertex(anewV, apar1, apar2, toFace, aTol)
        // — the UV update (bridge #7).
        builder_update_vertex_uv(brep, &anew_v, apar1, apar2, to_face, a_tol);
        anew_v
    }
}

/// OCCT L422-427 / L428-443: the (ppar1, ppar2) swap + the degenerate-range
/// stretch (the 3D row, comparing against toGC->First()/Last()).
fn fix_degenerate_range(
    mut ppar1: f64,
    mut ppar2: f64,
    row_first: f64,
    row_last: f64,
    preci: f64,
) -> (f64, f64) {
    if ppar1 > ppar2 {
        std::mem::swap(&mut ppar1, &mut ppar2);
    }
    if ppar2 - ppar1 < preci {
        if ppar1 - row_first < preci {
            ppar2 += 2.0 * preci;
        } else if row_last - ppar2 < preci {
            ppar1 -= 2.0 * preci;
        } else {
            ppar1 -= preci;
            ppar2 += preci;
        }
    }
    (ppar1, ppar2)
}

/// Idem for the pcurve row (OCCT L501-524).
fn fix_degenerate_range_pc(
    ppar1: f64,
    ppar2: f64,
    row_first: f64,
    row_last: f64,
    preci: f64,
) -> (f64, f64) {
    fix_degenerate_range(ppar1, ppar2, row_first, row_last, preci)
}

/// OCCT L452: toGC->IsCurveOnSurface() — the row kind test.
fn row_is_curve_on_surface(edge: &Shape, row: usize) -> bool {
    match edge.data.as_ref() {
        TShape::Edge(ed) => matches!(
            ed.representations.get(row),
            Some(CurveRepresentation::CurveOnSurface { .. })
                | Some(CurveRepresentation::CurveOnClosedSurface { .. })
        ),
        _ => false,
    }
}

/// OCCT static CorrectParameter (cxx L256-281): walks the trimmed/offset
/// basis curves down to the BSpline knots.
fn correct_parameter(crv: &Curve2d, param: f64) -> f64 {
    match crv {
        // OCCT L258-262: the Geom2d_TrimmedCurve basis.
        Curve2d::Trimmed(tmp) => correct_parameter(tmp.curve.as_ref(), param),
        // OCCT L263-267: the Geom2d_OffsetCurve basis.
        Curve2d::Offset(tmp) => correct_parameter(tmp.basis.as_ref(), param),
        // OCCT L268-279: the Geom2d_BSplineCurve knot snap.
        Curve2d::BSpline(bspline) => {
            // OCCT Knot(j) — the j-th distinct knot over the rcad expanded
            // knot vector (FirstUKnotIndex..=LastUKnotIndex).
            let eps = 1e-12;
            let mut distinct: Vec<f64> = Vec::new();
            for &k in &bspline.knots {
                if distinct.last().map_or(true, |&last| (k - last).abs() > eps) {
                    distinct.push(k);
                }
            }
            for &valknot in &distinct {
                if (valknot - param).abs() < PCONFUSION {
                    return valknot;
                }
            }
            param
        }
        _ => param,
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, U, V2, F, Tol) (bridge #7): the rcad
/// tool keeps the parameters on the vertex's PointOnSurface row and the
/// tolerance on the vertex.
fn builder_update_vertex_uv(
    brep: &mut BRep,
    the_v: &Shape,
    u: f64,
    v2: f64,
    the_face: &Shape,
    the_tol: f64,
) {
    // The face key (the legacy flat index stand-in of the rcad carrier).
    let face_key = the_face.index;
    let mut b = BRepBuilder::new();
    if u >= 0.5 * 2e100 || u <= -0.5 * 2e100 {
        panic!("Standard_DomainError: BRep_Builder::Infinite parameter");
    }
    let cur_tol = brep.vertex_tolerance(the_v);
    {
        let vd = brep.vertex_mut(the_v.clone());
        // Refresh (or append) the PointOnSurface row for the face.
        let mut found = false;
        for p in vd.points.iter_mut() {
            if let PointRepresentation::PointOnSurface {
                face, u: pu, v: pv, tolerance, ..
            } = p
            {
                if *face == face_key {
                    *pu = u;
                    *pv = v2;
                    *tolerance = the_tol;
                    found = true;
                    break;
                }
            }
        }
        if !found {
            vd.points.push(PointRepresentation::PointOnSurface {
                face: face_key,
                u,
                v: v2,
                tolerance: the_tol,
            });
        }
        let _ = &mut b;
    }
    if the_tol > cur_tol {
        brep.vertex_mut(the_v.clone()).tolerance = the_tol;
    }
    // The vertex position follows the UV update in OCCT
    // (BRep_Builder::UpdateVertex calls TV->UpdatePoints); the rcad point
    // row carries the update (the stored `point` stays — the OCCT
    // MakeBRep surface evaluation is driven by the consumers).
}

/// The legacy flat face-index lookup of the pre-W2 point rows (the
/// `PointOnSurface.face: usize` field).
fn face_surface_by_ptr_legacy(_face: usize) -> Option<Surface3> {
    // The rcad point rows carry the legacy flat index without a registry;
    // the from-surface identity stays unresolved (the OCCT
    // IsPointOnSurface(fromSurf, fromLoc) test then keeps its false arm).
    None
}
