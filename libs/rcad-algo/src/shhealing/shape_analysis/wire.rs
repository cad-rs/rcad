//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_Wire`
//! (`ShapeAnalysis_Wire.hxx` + `ShapeAnalysis_Wire.cxx` L1-2599).
//!
//! Tool for analysing a wire (edge order, connectivity, small edges, seams,
//! degenerated edges, gaps, self-intersection, lacking edges, tails, loops).
//!
//! Architecture bridges (numbering follows the curve.rs / surface.rs style):
//! 1. OCCT `occ::handle<ShapeExtend_WireData>` (myWire) -> the W1-3
//!    [`WireData`] value (owned clone; the handle sharing is not modelled —
//!    consumers re-`Load`).
//! 2. OCCT `TopoDS_Face` (myFace) -> [`Shape`] (null-able; `is_null()`).
//! 3. OCCT `handle<ShapeAnalysis_Surface>` (mySurf) -> the W2
//!    [`ShapeAnalysisSurface`].
//! 4. The status fields are the OCCT `Standard_Integer` bitmask aggregates
//!    over the W1-3 `ShapeExtendStatus` encoding.
//! 5. OCCT `Geom2dInt_GInter` (the curve-curve intersection with domains)
//!    -> GAP: the kernel `g_inter` translation carries the line-curve arm
//!    only; the curve-curve `Perform` used by the self-intersection checks
//!    is not translated, so the re-host keeps the OCCT `!IsDone()` failure
//!    branch live (`return false`), which is the OCCT behavior when the
//!    intersection computation fails.
//! 6. OCCT `BndLib_Add2dCurve::Add(gac, Confusion, box)` -> the kernel
//!    `base::bnd_lib::curve2d_bounding_box` (the BRepTools::AddUVBounds
//!    vehicle).
//! 7. OCCT `BRepTools::Compare(V1, V2)` -> the local re-host below
//!    (BRepTools.cxx: IsSame or the point distance within the summed
//!    tolerances).
//! 8. OCCT `ShapeAnalysis_TransferParametersProj` (the CheckTail cut) and
//!    `BRepGProp` surface/linear properties -> GAP: the TransferParameters
//!    pair is a W2 class of its own and the property machinery belongs to
//!    the GProp batch; CheckTail returns at the GAP anchor with the OCCT
//!    failure path (the cut result `false`).
//! 9. OCCT `ShapeAnalysis::FindBounds` -> the analysis.rs statics.

use glam::{DVec2, DVec3};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo::topods::{BRep, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use super::analysis::{self};
use super::curve::ShapeAnalysisCurve;
use super::edge::ShapeAnalysisEdge;
use super::surface::ShapeAnalysisSurface;
use super::wire_order::ShapeAnalysisWireOrder;
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;

/// OCCT Bnd_Box2d::IsOut(other) — the disjoint-box test (the kernel BndBox2d
/// carries the point form only; the box-box walk below).
fn bnd2d_is_out(a: &BndBox2d, b: &BndBox2d) -> bool {
    let (Some((axmin, aymin, axmax, aymax)), Some((bxmin, bymin, bxmax, bymax))) =
        (a.get(), b.get())
    else {
        // Void boxes are out of everything (the OCCT VoidMask semantics).
        return true;
    };
    bxmin > axmax || bxmax < axmin || bymin > aymax || bymax < aymin
}

// OCCT gp.hxx L59-60: gp::Resolution() = RealSmall() = DBL_MIN.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;
// OCCT Standard_Real.hxx L182-185.
const REAL_LAST: f64 = f64::MAX;

// ---------------------------------------------------------------------------
// Local re-hosts
// ---------------------------------------------------------------------------

/// OCCT BRepTools::Compare(V1, V2) (BRepTools.cxx) — same TShape or the
/// points within the summed tolerances (bridge #7).
pub(crate) fn brep_tools_compare(brep: &BRep, v1: &Shape, v2: &Shape) -> bool {
    if std::sync::Arc::as_ptr(&v1.data) == std::sync::Arc::as_ptr(&v2.data)
        && v1.location == v2.location
    {
        return true;
    }
    let (p1, t1) = match v1.data.as_ref() {
        TShape::Vertex(vd) => (vd.point, vd.tolerance),
        _ => return false,
    };
    let (p2, t2) = match v2.data.as_ref() {
        TShape::Vertex(vd) => (vd.point, vd.tolerance),
        _ => return false,
    };
    p1.distance(p2) <= t1 + t2
}

/// OCCT BRep_Tool::Pnt(V).
pub(crate) fn brep_tool_pnt(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Degenerated(E).
pub(crate) fn brep_tool_degenerated(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::SameParameter(E).
pub(crate) fn brep_tool_same_parameter(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.same_parameter,
        _ => false,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edge, face, cf, cl) — the stored (face)
/// overload matched by the face TShape pointer.
pub(crate) fn brep_tool_curve_on_surface_face(
    the_edge: &Shape,
    the_face: &Shape,
) -> Option<(rcad_kernel::geom::Curve2d, f64, f64)> {
    let fptr = std::sync::Arc::as_ptr(&the_face.data) as u64;
    let reps = match the_edge.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return None,
    };
    for r in reps {
        let (key, pc, range) = match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        if key.0 == fptr {
            return Some((pc.clone(), range[0], range[1]));
        }
    }
    None
}

// ---------------------------------------------------------------------------
// ShapeAnalysis_Wire
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_Wire (hxx L25-140).
pub struct ShapeAnalysisWire {
    /// OCCT myWire (bridge #1).
    pub(crate) my_wire: Option<WireData>,
    /// OCCT myFace (bridge #2).
    pub(crate) my_face: Shape,
    /// OCCT mySurf (bridge #3).
    pub(crate) my_surf: Option<ShapeAnalysisSurface>,
    /// OCCT myPrecision.
    pub(crate) my_precision: f64,
    /// OCCT myStatusOrder / myStatusConnected / myStatusEdgeCurves /
    /// myStatusDegenerated / myStatusClosed / myStatusLacking /
    /// myStatusSelfIntersection / myStatusSmall / myStatusGaps3d /
    /// myStatusGaps2d / myStatusCurveGaps / myStatusLoop / myStatus.
    pub(crate) my_status_order: i32,
    pub(crate) my_status_connected: i32,
    pub(crate) my_status_edge_curves: i32,
    pub(crate) my_status_degenerated: i32,
    pub(crate) my_status_closed: i32,
    pub(crate) my_status_lacking: i32,
    pub(crate) my_status_self_intersection: i32,
    pub(crate) my_status_small: i32,
    pub(crate) my_status_gaps3d: i32,
    pub(crate) my_status_gaps2d: i32,
    pub(crate) my_status_curve_gaps: i32,
    pub(crate) my_status_loop: i32,
    pub(crate) my_status: i32,
    /// OCCT myMin3d / myMin2d / myMax3d / myMax2d.
    pub(crate) my_min3d: f64,
    pub(crate) my_min2d: f64,
    pub(crate) my_max3d: f64,
    pub(crate) my_max2d: f64,
}

impl Default for ShapeAnalysisWire {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisWire {
    /// OCCT ShapeAnalysis_Wire() (cxx L88-92).
    pub fn new() -> Self {
        let mut this = ShapeAnalysisWire {
            my_wire: None,
            my_face: Shape::null(),
            my_surf: None,
            my_precision: 0.0,
            my_status_order: 0,
            my_status_connected: 0,
            my_status_edge_curves: 0,
            my_status_degenerated: 0,
            my_status_closed: 0,
            my_status_lacking: 0,
            my_status_self_intersection: 0,
            my_status_small: 0,
            my_status_gaps3d: 0,
            my_status_gaps2d: 0,
            my_status_curve_gaps: 0,
            my_status_loop: 0,
            my_status: 0,
            my_min3d: 0.0,
            my_min2d: 0.0,
            my_max3d: 0.0,
            my_max2d: 0.0,
        };
        this.clear_statuses();
        this.my_precision = CONFUSION;
        this
    }

    /// OCCT ShapeAnalysis_Wire(wire, face, precision) (cxx L96-102).
    pub fn new_from_wire(brep: &mut BRep, wire: &Shape, face: &Shape, precision: f64) -> Self {
        let mut this = Self::new();
        this.init(brep, wire, face, precision);
        this
    }

    /// OCCT ShapeAnalysis_Wire(sbwd, face, precision) (cxx L105-111).
    pub fn new_from_wire_data(
        brep: &mut BRep,
        sbwd: WireData,
        face: &Shape,
        precision: f64,
    ) -> Self {
        let mut this = Self::new();
        this.init_wire_data(brep, sbwd, face, precision);
        this
    }

    /// OCCT Init(wire, face, precision) (cxx L114-121).
    pub fn init(&mut self, brep: &mut BRep, wire: &Shape, face: &Shape, precision: f64) {
        let sbwd = WireData::new_from_wire(brep, wire, false, true);
        self.init_wire_data(brep, sbwd, face, precision);
    }

    /// OCCT Init(sbwd, face, precision) (cxx L123-131).
    pub fn init_wire_data(
        &mut self,
        brep: &mut BRep,
        sbwd: WireData,
        face: &Shape,
        precision: f64,
    ) {
        self.load_wire_data(sbwd);
        self.set_face(brep, face);
        self.set_precision(precision);
    }

    /// OCCT Load(wire) (cxx L134-140).
    pub fn load(&mut self, brep: &mut BRep, wire: &Shape) {
        self.clear_statuses();
        self.my_wire = Some(WireData::new_from_wire(brep, wire, false, true));
    }

    /// OCCT Load(sbwd) (cxx L142-148).
    pub fn load_wire_data(&mut self, sbwd: WireData) {
        self.clear_statuses();
        self.my_wire = Some(sbwd);
    }

    /// OCCT SetFace(face) (cxx L150-165).
    pub fn set_face(&mut self, brep: &mut BRep, face: &Shape) {
        self.my_face = face.clone();
        if self.my_face.is_null() {
            return;
        }
        // BRep_Tool::Surface(myFace) — the TFace surface (bridge #2).
        let a_surface = match face.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        if !face.is_null() {
            if let Some(a_surface) = a_surface {
                self.my_surf = Some(ShapeAnalysisSurface::new(a_surface));
            }
        }
        let _ = brep;
    }

    /// OCCT SetFace(theFace, theSurfaceAnalysis) (cxx L167-174).
    pub fn set_face_with_surface(
        &mut self,
        the_face: &Shape,
        the_surface_analysis: ShapeAnalysisSurface,
    ) {
        self.my_face = the_face.clone();
        self.my_surf = Some(the_surface_analysis);
    }

    /// OCCT SetSurface(theSurfaceAnalysis) (cxx L176-181).
    pub fn set_surface(&mut self, the_surface_analysis: ShapeAnalysisSurface) {
        self.my_surf = Some(the_surface_analysis);
    }

    /// OCCT SetSurface(surface) (cxx L183-188).
    pub fn set_surface_geom(&mut self, brep: &mut BRep, surface: &rcad_kernel::geom::Surface3) {
        self.set_surface_geom_loc(brep, surface, 0);
    }

    /// OCCT SetSurface(surface, location) (cxx L190-199): builds a face from
    /// the surface at Confusion and delegates to SetFace.
    pub fn set_surface_geom_loc(
        &mut self,
        brep: &mut BRep,
        surface: &rcad_kernel::geom::Surface3,
        _location: u32,
    ) {
        let face = rcad_kernel::topo::topods::BRepBuilder::new().make_face(
            brep,
            Some(surface.clone()),
            Shape::null(),
        );
        self.set_face(brep, &face);
    }

    /// OCCT SetPrecision(precision) (cxx L201-206).
    pub fn set_precision(&mut self, precision: f64) {
        self.my_precision = precision;
    }

    /// OCCT Precision().
    pub fn precision(&self) -> f64 {
        self.my_precision
    }

    /// OCCT ClearStatuses() (cxx L208-217).
    pub fn clear_statuses(&mut self) {
        self.my_status_order = 0;
        self.my_status_connected = 0;
        self.my_status_edge_curves = 0;
        self.my_status_degenerated = 0;
        self.my_status_closed = 0;
        self.my_status_lacking = 0;
        self.my_status_self_intersection = 0;
        self.my_status_small = 0;
        self.my_status_gaps3d = 0;
        self.my_status_gaps2d = 0;
        self.my_status_curve_gaps = 0;
        self.my_status_loop = 0;
        self.my_status = 0;
        self.my_min3d = 0.0;
        self.my_min2d = 0.0;
        self.my_max3d = 0.0;
        self.my_max2d = 0.0;
    }

    /// OCCT Load(wire) state accessors.
    pub fn wire_data(&self) -> Option<&WireData> {
        self.my_wire.as_ref()
    }

    /// Mutable access to the loaded wire data.  OCCT Load(sbwd) stores the
    /// caller's ShapeExtend_WireData handle, so the local `sewd` handle and
    /// `saw->WireData()` are the same object (the free_bounds connection
    /// flow mutates it through both); the rcad value model resolves the
    /// aliasing to a single owner with this accessor.
    pub fn wire_data_mut(&mut self) -> Option<&mut WireData> {
        self.my_wire.as_mut()
    }

    pub fn face(&self) -> &Shape {
        &self.my_face
    }

    pub fn surf(&self) -> Option<&ShapeAnalysisSurface> {
        self.my_surf.as_ref()
    }

    /// OCCT NbEdges() (hxx): myWire->NbEdges().
    pub fn nb_edges(&self) -> i32 {
        match self.my_wire.as_ref() {
            Some(w) => w.nb_edges(),
            None => 0,
        }
    }

    /// OCCT IsLoaded() (hxx).
    pub fn is_loaded(&self) -> bool {
        self.nb_edges() > 0
    }

    /// OCCT IsReady() (hxx): IsLoaded && !myFace.IsNull().
    pub fn is_ready(&self) -> bool {
        self.is_loaded() && !self.my_face.is_null()
    }

    // -- status queries (hxx) ----------------------------------------------

    /// OCCT StatusOrder(flag).
    pub fn status_order(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_order, status)
    }

    /// OCCT StatusConnected(flag).
    pub fn status_connected(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_connected, status)
    }

    /// OCCT StatusEdgeCurves(flag).
    pub fn status_edge_curves(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_edge_curves, status)
    }

    /// OCCT StatusDegenerated(flag).
    pub fn status_degenerated(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_degenerated, status)
    }

    /// OCCT StatusSelfIntersection(flag).
    pub fn status_self_intersection(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_self_intersection, status)
    }

    /// OCCT StatusSmall(flag).
    pub fn status_small(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_small, status)
    }

    /// OCCT StatusLacking(flag).
    pub fn status_lacking(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_lacking, status)
    }

    /// OCCT StatusClosed(flag).
    pub fn status_closed(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_closed, status)
    }

    /// OCCT StatusGaps3d(flag).
    pub fn status_gaps3d(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_gaps3d, status)
    }

    /// OCCT StatusGaps2d(flag).
    pub fn status_gaps2d(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_gaps2d, status)
    }

    /// OCCT StatusCurveGaps(flag).
    pub fn status_curve_gaps(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_curve_gaps, status)
    }

    /// OCCT StatusLoop(flag).
    pub fn status_loop(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_loop, status)
    }

    /// OCCT LastCheckStatus(flag).
    pub fn last_check_status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    /// OCCT MinDistance3d() / MaxDistance3d() (hxx).
    pub fn min_distance_3d(&self) -> f64 {
        self.my_min3d
    }

    pub fn max_distance_3d(&self) -> f64 {
        self.my_max3d
    }

    /// OCCT MinDistance2d() / MaxDistance2d() (hxx).
    pub fn min_distance_2d(&self) -> f64 {
        self.my_min2d
    }

    pub fn max_distance_2d(&self) -> f64 {
        self.my_max2d
    }

    // -- top-level checks ---------------------------------------------------

    /// OCCT Perform() (cxx L219-233).
    pub fn perform(&mut self, brep: &mut BRep) -> bool {
        let mut result = false;
        result |= self.check_order_top(brep, false, true);
        result |= self.check_small_top(brep, 0.0);
        result |= self.check_connected_top(brep, 0.0);
        result |= self.check_edge_curves(brep);
        result |= self.check_degenerated_top(brep);
        result |= self.check_self_intersection(brep);
        result |= self.check_lacking_top(brep);
        result |= self.check_closed(brep, 0.0);
        result
    }

    /// OCCT CheckOrder(isClosed = false, mode3d = true) (cxx L235-243).
    pub fn check_order_top(&mut self, brep: &mut BRep, is_closed: bool, mode3d: bool) -> bool {
        let mut sawo = ShapeAnalysisWireOrder::new();
        self.check_order(brep, &mut sawo, is_closed, mode3d, false);
        self.my_status_order = self.my_status;
        self.status_order(ShapeExtendStatus::Done)
    }

    /// OCCT CheckSmall(precsmall = -1) (cxx L245-255).
    pub fn check_small_top(&mut self, brep: &mut BRep, precsmall: f64) -> bool {
        let precsmall = if precsmall == 0.0 { -1.0 } else { precsmall };
        let nb = self.nb_edges();
        for i in 1..=nb {
            self.check_small(brep, i, precsmall);
            self.my_status_small |= self.my_status;
        }
        self.status_small(ShapeExtendStatus::Done)
    }

    /// OCCT CheckConnected(prec = 0.0) (cxx L257-267).
    pub fn check_connected_top(&mut self, brep: &mut BRep, prec: f64) -> bool {
        let nb = self.nb_edges();
        for i in 1..=nb {
            self.check_connected(brep, i, prec);
            self.my_status_connected |= self.my_status;
        }
        self.status_connected(ShapeExtendStatus::Done)
    }

    /// OCCT CheckEdgeCurves() (cxx L269-358).
    pub fn check_edge_curves(&mut self, brep: &mut BRep) -> bool {
        self.my_status_edge_curves = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let nb = self.nb_edges();
        let mut sae = ShapeAnalysisEdge::new();

        for i in 1..=nb {
            let e = self.my_wire.as_ref().unwrap().edge(i);

            sae.check_curve3d_with_pcurve_face(brep, &e, &self.my_face);
            if sae.status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done1);
            }
            if sae.status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail1);
            }

            sae.check_vertices_with_pcurve_face(brep, &e, &self.my_face, 0.0, 0);
            if sae.status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done2);
            }
            if sae.status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail2);
            }

            sae.check_vertices_with_curve3d(brep, &e, 0.0, 0);
            if sae.status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done3);
            }
            if sae.status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail3);
            }

            self.check_seam(brep, i);
            if self.last_check_status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done4);
            }
            if self.last_check_status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail4);
            }

            self.check_gap3d(brep, i);
            if self.last_check_status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done5);
            }
            if self.last_check_status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail5);
            }

            self.check_gap2d(brep, i);
            if self.last_check_status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done6);
            }
            if self.last_check_status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail6);
            }

            let mut maxdev = 0.0f64;
            sae.check_same_parameter_face(brep, &e, &self.my_face, &mut maxdev, 23);
            if sae.status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done7);
            }
            if sae.status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail7);
            }
        }
        self.status_edge_curves(ShapeExtendStatus::Done)
    }

    /// OCCT CheckDegenerated() (cxx L360-368).
    pub fn check_degenerated_top(&mut self, brep: &mut BRep) -> bool {
        let nb = self.nb_edges();
        for i in 1..=nb {
            self.check_degenerated(brep, i);
            self.my_status_degenerated |= self.my_status;
        }
        self.status_degenerated(ShapeExtendStatus::Done)
    }

    /// OCCT CheckSelfIntersection() (cxx L372-452).
    pub fn check_self_intersection(&mut self, brep: &mut BRep) -> bool {
        self.my_status_self_intersection = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }
        let nb = self.nb_edges();
        for i in 1..=nb {
            self.check_self_intersecting_edge(brep, i);
            if self.last_check_status(ShapeExtendStatus::Done) {
                self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done1);
            }
            if self.last_check_status(ShapeExtendStatus::Fail) {
                self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Fail1);
            }

            self.check_intersecting_edges(brep, i);
            if self.last_check_status(ShapeExtendStatus::Done) {
                self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done2);
            }
            if self.last_check_status(ShapeExtendStatus::Fail) {
                self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Fail2);
            }
        }

        let mut boxes: Vec<BndBox2d> = (0..nb as usize).map(|_| BndBox2d::new()).collect();
        // BRep_Tool::Surface(Face(), L) — the face surface value.
        let s = match self.my_face.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
        let (mut cf, mut cl) = (0.0f64, 0.0f64);
        let mut sae = ShapeAnalysisEdge::new();
        for i in 1..=nb {
            let e = self.my_wire.as_ref().unwrap().edge(i);
            let pcurve_ok = match s.as_ref() {
                Some(s) => sae.pcurve_surface(brep, &e, s, 0, &mut c2d, &mut cf, &mut cl, false),
                None => false,
            };
            if pcurve_ok {
                let mut box2d = BndBox2d::new();
                if let Some(c2d) = c2d.as_ref() {
                    // BndLib_Add2dCurve::Add(gac, Confusion, box) — the
                    // kernel curve2d bounding walk (bridge #6).
                    let [xmin, ymin, xmax, ymax] = rcad_kernel::base::bnd_lib::curve2d_bounding_box(
                        c2d, cf, cl, CONFUSION,
                    );
                    box2d.update(xmin, ymin, xmax, ymax);
                }
                boxes[(i - 1) as usize] = box2d;
            }
        }

        let mut is_fail = false;
        let mut is_done = false;
        for num1 in 1..(nb - 1).max(1) {
            let mut fin = nb;
            if self.check_closed(brep, CONFUSION) && 1 == num1 {
                fin = nb - 1;
            }
            for num2 in (num1 + 2)..=fin {
                if !bnd2d_is_out(&boxes[(num1 - 1) as usize], &boxes[(num2 - 1) as usize]) {
                    self.check_intersecting_edges_pair(brep, num1, num2);
                    is_fail |= self.last_check_status(ShapeExtendStatus::Fail1);
                    is_done |= self.last_check_status(ShapeExtendStatus::Done1);
                }
            }
        }
        if is_fail {
            self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Fail3);
        }
        if is_done {
            self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done3);
        }

        self.status_self_intersection(ShapeExtendStatus::Done)
    }

    /// OCCT CheckLacking() (cxx L454-468).
    pub fn check_lacking_top(&mut self, brep: &mut BRep) -> bool {
        if !self.is_ready() || self.nb_edges() < 2 {
            return false;
        }
        let nb = self.nb_edges();
        for i in 1..=nb {
            self.check_lacking(brep, i, 0.0);
            self.my_status_lacking |= self.my_status;
        }
        self.status_lacking(ShapeExtendStatus::Done)
    }

    /// OCCT CheckClosed(prec = 0.0) (cxx L470-501).
    pub fn check_closed(&mut self, brep: &mut BRep, prec: f64) -> bool {
        self.my_status_closed = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 1 {
            return false;
        }

        self.check_connected(brep, 1, prec);
        if self.last_check_status(ShapeExtendStatus::Done) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Done1);
        }
        if self.last_check_status(ShapeExtendStatus::Fail) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Fail1);
        }

        self.check_degenerated(brep, 1);
        if self.last_check_status(ShapeExtendStatus::Done) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Done2);
        }
        if self.last_check_status(ShapeExtendStatus::Fail) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Fail2);
        }

        self.status_closed(ShapeExtendStatus::Done)
    }

    /// OCCT CheckGaps3d() (cxx L503-531).
    pub fn check_gaps3d(&mut self, brep: &mut BRep) -> bool {
        self.my_status_gaps3d = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() < 1 {
            return false; // gka IsLoaded
        }

        let mut dist;
        let mut maxdist = 0.0f64;

        for i in 1..=self.nb_edges() {
            self.check_gap3d(brep, i);
            self.my_status_gaps3d |= self.my_status;
            if !self.last_check_status(ShapeExtendStatus::Fail1) {
                dist = self.min_distance_3d();
                if maxdist < dist {
                    maxdist = dist;
                }
            }
        }
        self.my_min3d = maxdist;
        self.my_max3d = maxdist;

        self.status_gaps3d(ShapeExtendStatus::Done)
    }

    /// OCCT CheckGaps2d() (cxx L533-561).
    pub fn check_gaps2d(&mut self, brep: &mut BRep) -> bool {
        self.my_status_gaps2d = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 1 {
            return false;
        }

        let mut dist;
        let mut maxdist = 0.0f64;

        for i in 1..=self.nb_edges() {
            self.check_gap2d(brep, i);
            self.my_status_gaps2d |= self.my_status;
            if !self.last_check_status(ShapeExtendStatus::Fail1) {
                dist = self.min_distance_2d();
                if maxdist < dist {
                    maxdist = dist;
                }
            }
        }
        self.my_min2d = maxdist;
        self.my_max2d = maxdist;

        self.status_gaps2d(ShapeExtendStatus::Done)
    }

    /// OCCT CheckCurveGaps() (cxx L563-591).
    pub fn check_curve_gaps(&mut self, brep: &mut BRep) -> bool {
        self.my_status_curve_gaps = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 1 {
            return false;
        }

        let mut dist;
        let mut maxdist = 0.0f64;

        for i in 1..=self.nb_edges() {
            self.check_curve_gap(brep, i);
            self.my_status_curve_gaps |= self.my_status;
            if !self.last_check_status(ShapeExtendStatus::Fail1) {
                dist = self.min_distance_3d();
                if maxdist < dist {
                    maxdist = dist;
                }
            }
        }
        self.my_min3d = maxdist;
        self.my_max3d = maxdist;

        self.status_curve_gaps(ShapeExtendStatus::Done)
    }

    // -- per-edge checks ----------------------------------------------------

    /// OCCT CheckOrder(sawo, isClosed, theMode3D, theModeBoth)
    /// (cxx L593-690).
    pub fn check_order(
        &mut self,
        brep: &mut BRep,
        sawo: &mut ShapeAnalysisWireOrder,
        the_mode3d: bool,
        the_mode_both: bool,
        is_closed: bool,
    ) -> bool {
        let _ = brep;
        if (!the_mode3d || the_mode_both) && self.my_face.is_null() {
            self.my_status = encode_status(ShapeExtendStatus::Fail2);
            return false;
        }
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        sawo.set_mode(the_mode3d, 0.0, the_mode_both);
        let nb = self.nb_edges();
        let ea = ShapeAnalysisEdge::new();
        let mut is_all2d_edges_ok = true;
        for i in 1..=nb {
            let e = self.my_wire.as_ref().unwrap().edge(i);
            let mut a_p1xyz = DVec3::ZERO;
            let mut a_p2xyz = DVec3::ZERO;
            let mut a_p1xy = DVec2::ZERO;
            let mut a_p2xy = DVec2::ZERO;
            if the_mode3d || the_mode_both {
                let v1 = ea.first_vertex(brep, &e);
                let v2 = ea.last_vertex(brep, &e);
                if v1.is_null() || v2.is_null() {
                    self.my_status = encode_status(ShapeExtendStatus::Fail2);
                    return false;
                }
                a_p1xyz = brep_tool_pnt(&v1);
                a_p2xyz = brep_tool_pnt(&v2);
            }
            if !the_mode3d || the_mode_both {
                let mut f = 0.0f64;
                let mut l = 0.0f64;
                let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
                // myFace.Oriented(TopAbs_FORWARD) — the same face value (the
                // pcurve matching is by the face pointer).
                if !ea.pcurve_face(brep, &e, &self.my_face, &mut c2d, &mut f, &mut l, false) {
                    // if mode is 2d, then we can nothing to do, else we can
                    // switch to 3d mode
                    if !the_mode3d && !the_mode_both {
                        self.my_status = encode_status(ShapeExtendStatus::Fail2);
                        return false;
                    }
                    is_all2d_edges_ok = false;
                } else if let Some(c2d) = c2d.as_ref() {
                    use rcad_kernel::geom::Curve2dEval;
                    a_p1xy = Curve2dEval::point_at(c2d, f);
                    a_p2xy = Curve2dEval::point_at(c2d, l);
                }
            }
            if the_mode3d && !the_mode_both {
                sawo.add_3d(a_p1xyz, a_p2xyz);
            } else if !the_mode3d && !the_mode_both {
                sawo.add_2d(a_p1xy, a_p2xy);
            } else {
                sawo.add_both(a_p1xyz, a_p2xyz, a_p1xy, a_p2xy);
            }
        }
        // need to switch to 3d mode
        if the_mode_both && !is_all2d_edges_ok {
            sawo.set_mode(true, 0.0, false);
        }
        sawo.perform(is_closed);
        let stat = sawo.status();
        match stat {
            0 => self.my_status = encode_status(ShapeExtendStatus::Ok),
            1 => self.my_status = encode_status(ShapeExtendStatus::Done1),
            // this value is not returned
            2 => self.my_status = encode_status(ShapeExtendStatus::Done2),
            -1 => self.my_status = encode_status(ShapeExtendStatus::Done3),
            // this value is not returned
            -2 => self.my_status = encode_status(ShapeExtendStatus::Done4),
            // only shifted
            3 => self.my_status = encode_status(ShapeExtendStatus::Done5),
            // this value is not returned
            -10 => self.my_status = encode_status(ShapeExtendStatus::Fail1),
            _ => {}
        }
        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckConnected(num, prec) (cxx L693-762).
    pub fn check_connected(&mut self, brep: &mut BRep, num: i32, prec: f64) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() < 1 {
            return false;
        }

        let n2 = if num > 0 { num } else { self.nb_edges() };
        let n1 = if n2 > 1 { n2 - 1 } else { self.nb_edges() };

        let e1 = self.my_wire.as_ref().unwrap().edge(n1);
        let e2 = self.my_wire.as_ref().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.last_vertex(brep, &e1);
        let v2 = sae.first_vertex(brep, &e2);
        if v1.is_null() || v2.is_null() {
            self.my_status = encode_status(ShapeExtendStatus::Fail2);
            return false;
        }
        if v1.ptr_id() == v2.ptr_id() && v1.location == v2.location {
            return false;
        }

        let p1 = brep_tool_pnt(&v1);
        let mut p2 = brep_tool_pnt(&v2);
        self.my_min3d = p1.distance(p2);
        if self.my_min3d <= GP_RESOLUTION {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
        } else if self.my_min3d <= self.my_precision {
            self.my_status |= encode_status(ShapeExtendStatus::Done2);
        } else if self.my_min3d <= prec {
            self.my_status |= encode_status(ShapeExtendStatus::Done3);
        } else {
            // et en inversant la derniere edge ?
            if n1 == n2 {
                self.my_status = encode_status(ShapeExtendStatus::Fail1);
            } else {
                let v2r = sae.last_vertex(brep, &e2);
                p2 = brep_tool_pnt(&v2r);
                let dist = p1.distance(p2);
                if dist > self.my_precision {
                    self.my_status = encode_status(ShapeExtendStatus::Fail1);
                } else {
                    self.my_min3d = dist;
                    self.my_status = encode_status(ShapeExtendStatus::Fail2);
                }
            }
            return false;
        }
        true
    }

    /// OCCT CheckSmall(num, precsmall) (cxx L765-837).
    pub fn check_small(&mut self, brep: &mut BRep, num: i32, precsmall: f64) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() <= 1 {
            return false;
        }

        let e = self
            .my_wire
            .as_ref()
            .unwrap()
            .edge(if num != 0 { num } else { self.nb_edges() });
        let sae = ShapeAnalysisEdge::new();

        if brep_tool_degenerated(&e) {
            //: n2 abv 22 Jan 99: degen edge with no pcurve should be removed
            if !self.my_face.is_null() && sae.has_pcurve_face(brep, &e, &self.my_face) {
                return false;
            }
            self.my_status = encode_status(ShapeExtendStatus::Fail1);
        }

        let v1 = sae.first_vertex(brep, &e);
        let v2 = sae.last_vertex(brep, &e);
        if v1.is_null() || v2.is_null() {
            self.my_status = encode_status(ShapeExtendStatus::Fail2);
            return false;
        }
        let p1 = brep_tool_pnt(&v1);
        let p2 = brep_tool_pnt(&v2);
        let dist = p1.distance(p2);
        let prec = precsmall; // Min ( myPrecision, precsmall );
        if dist > prec {
            return false; // not small enough
        }

        // The 3D curve now: is it CLOSED or ZERO LENGTH...???
        let mut a_midpoint;
        let (mut cf, mut cl) = (0.0f64, 0.0f64);
        let mut c3d: Option<rcad_kernel::geom::Curve3> = None;
        if sae.curve3d(brep, &e, &mut c3d, &mut cf, &mut cl, false) {
            use rcad_kernel::geom::CurveEval;
            a_midpoint =
                CurveEval::point_at(c3d.as_ref().unwrap(), (cf + cl) / 2.0);
        } else {
            let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
            if !self.my_face.is_null()
                && self.my_surf.is_some()
                && sae.pcurve_face(brep, &e, &self.my_face, &mut c2d, &mut cf, &mut cl, false)
            {
                use rcad_kernel::geom::Curve2dEval;
                let p2m = Curve2dEval::point_at(c2d.as_ref().unwrap(), (cf + cl) / 2.0);
                a_midpoint = self.my_surf.as_ref().unwrap().value(p2m);
            } else {
                self.my_status = encode_status(ShapeExtendStatus::Fail1);
                a_midpoint = p1;
            }
        }
        if a_midpoint.distance(p1) > prec || a_midpoint.distance(p2) > prec {
            return false;
        }

        let same = v1.ptr_id() == v2.ptr_id() && v1.location == v2.location;
        self.my_status |= encode_status(if same {
            ShapeExtendStatus::Done1
        } else {
            ShapeExtendStatus::Done2
        });
        true
    }

    /// OCCT CheckSeam(num, C1, C2, cf, cl) (cxx L839-885).
    pub fn check_seam_full(
        &mut self,
        brep: &mut BRep,
        num: i32,
        c1: &mut Option<rcad_kernel::geom::Curve2d>,
        c2: &mut Option<rcad_kernel::geom::Curve2d>,
        cf: &mut f64,
        cl: &mut f64,
    ) -> bool {
        let _ = brep;
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }
        let n = if num == 0 { self.nb_edges() } else { num };
        let e = self.my_wire.as_ref().unwrap().edge(n);
        let sae = ShapeAnalysisEdge::new();
        if !sae.is_seam_face(brep, &e, &self.my_face) {
            return false;
        }
        // Extract the Two PCurves of the Seam (the FORWARD-oriented face; the
        // pcurve matching is by the face pointer, orientation-insensitive).
        let (pcf1, pcf2) = match e.data.as_ref() {
            TShape::Edge(ed) => {
                let mut p1 = None;
                let mut p2 = None;
                let fptr = std::sync::Arc::as_ptr(&self.my_face.data) as u64;
                for r in &ed.representations {
                    if let rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                        face,
                        pcurve1,
                        pcurve2,
                        ..
                    } = r
                    {
                        if face.0 == fptr {
                            p1 = Some(pcurve1.clone());
                            p2 = Some(pcurve2.clone());
                        }
                    }
                }
                (p1, p2)
            }
            _ => (None, None),
        };
        let (Some(pc1), Some(pc2)) = (pcf1, pcf2) else {
            return false;
        };
        *c1 = Some(pc1);
        *c2 = Some(pc2);

        //  SelectForward est destine a devenir un outil distinct
        let sac = ShapeAnalysisCurve;
        let the_curve_indice =
            sac.select_forward_seam(c1.as_ref().unwrap(), c2.as_ref().unwrap());
        if the_curve_indice != 2 {
            return false;
        }

        self.my_status = encode_status(ShapeExtendStatus::Done1);
        true
    }

    /// OCCT CheckSeam(num) (cxx L887-894).
    pub fn check_seam(&mut self, brep: &mut BRep, num: i32) -> bool {
        let mut c1: Option<rcad_kernel::geom::Curve2d> = None;
        let mut c2: Option<rcad_kernel::geom::Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        self.check_seam_full(brep, num, &mut c1, &mut c2, &mut cf, &mut cl)
    }

    // ... (the remaining checks are in wire_checks.rs)
}
