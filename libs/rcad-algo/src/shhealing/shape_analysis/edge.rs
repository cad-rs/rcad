//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_Edge`
//! (`ShapeAnalysis_Edge.hxx` L17-259 + `ShapeAnalysis_Edge.cxx` L1-1033).
//!
//! Tool for analyzing the edge: queries geometrical representations of the
//! edge (3d curve, pcurve on the given face or surface) and topological
//! sub-shapes (bounding vertices), and provides methods for analyzing
//! geometry/topology consistency.
//!
//! Architecture bridges (numbered, referenced by the methods below):
//! 1. `BRep` pool argument - OCCT `BRep_Tool` / `TopExp` read the TShape
//!    graph through global accessors over the shape document. The rcad
//!    equivalents resolve TShape data, the TopLoc_Location table and the
//!    pcurve-row face registry through `rcad_kernel::BRep`; every method
//!    takes `brep: &BRep` as that stand-in (shape_build/brep_tool.rs and
//!    kernel BRepTool precedents).
//! 2. `TopLoc_Location` -> `u32` (the `BRep.locations` table index,
//!    0 = identity).
//! 3. `Geom_Curve` / `Geom2d_Curve` / `Geom_Surface` handles ->
//!    `Curve3` / `Curve2d` / `Surface3` values (cloned).
//! 4. `gp_Pnt` / `gp_Pnt2d` / `gp_Vec2d` -> `DVec3` / `DVec2` / `DVec2`.
//! 5. Pcurve-key location component: OCCT stores a curve representation
//!    under `L.Predivided(E.Location())` (a TopLoc_Location); rcad stores
//!    the value hash `compose_pcurve_location(face_loc, edge_loc)` in the
//!    (face_ptr, loc) pcurve key. The (surface, location) representation
//!    walk compares that hash; recomposing a full location matrix from a
//!    hash is not possible, so the transformed points/curves below use the
//!    edge wrapper location only (the identity-location reduction of the
//!    loc_ope_* precedents).
//! 6. Output parameters use `&mut` (the AGENTS 1:1 rule); OCCT overloads
//!    are disambiguated with `_face` / `_surface` suffixes.

use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::base::gcpnts::abscissa_point::{abscissa_point_parameter, arc_length};
use rcad_kernel::base::geom_proj_lib::project_on_plane;
use rcad_kernel::geom::{
    transform_curve, transform_surface, Curve2d, Curve2dEval, Curve3, CurveEval, Surface3,
    SurfaceEval,
};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{
    compose_pcurve_location, surface_same, tshape_flags, CurveRepresentation, Orientation, TShape,
    TVertexData,
};
use rcad_kernel::BRep;

use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::topalgo::brep_lib_validate_edge::{
    Adaptor3dCurveOnSurface, BRepLibValidateEdge, Geom2dAdaptorCurve, GeomAdaptorCurve,
    GeomAdaptorSurface,
};

// OCCT Standard_Real.hxx L179-186: RealLast() - the biggest representable real.
const REAL_LAST: f64 = f64::MAX;

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (architecture bridge #1; the loc_ope_* precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Curve(edg, Loc, f, l) (BRep_Tool.cxx): the 3D curve and
/// parameter range; the out location is the edge wrapper location (bridge
/// #5: the 3D curve representation carries no own location in rcad).
fn brep_tool_curve_loc(_brep: &BRep, edg: &Shape) -> (Option<Curve3>, u32, f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), edg.location, ed.range[0], ed.range[1]),
        _ => (None, 0, 0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Curve(edg, f, l) - the no-location variant.
fn brep_tool_curve(_brep: &BRep, edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// The edge's curve representations (empty for non-edges).
fn edge_representations(edg: &Shape) -> &[CurveRepresentation] {
    match edg.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => &[],
    }
}

/// The face surface registered in the pool under a TShape pointer (the
/// Geom_Surface handle identity stand-in; kernel face_surface_by_ptr
/// precedent, topods.rs L2236).
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

/// OCCT BRep_Tool::Surface(fac, l) - the LOCAL surface with the out location
/// (bridge #5: the rcad TFace stores the surface already in its TFace frame,
/// so the out location is the occurrence location).
fn brep_tool_surface_loc(_brep: &BRep, fac: &Shape) -> (Option<Surface3>, u32) {
    match fac.data.as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fac.location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edge, surface, location, cf, cl)
/// (BRep_Tool.cxx L345-367, the stored-only (surface, location) overload):
/// the edge's curve representations are matched by the surface VALUE
/// (`cr->IsCurveOnSurface(S, loc)` with `loc = L.Predivided(E.Location())`)
/// and the matched BRep_GCurve supplies PCurve() (the first pcurve; a seam's
/// BRep_CurveOnClosedSurface matches with PCurve1).  Bridge #5: rcad matches
/// the pcurve-key location hash and resolves the row surface through the
/// pool registry.
fn brep_tool_curve_on_surface_stored(
    brep: &BRep,
    edg: &Shape,
    surface: &Surface3,
    location: u32,
) -> Option<(Curve2d, f64, f64)> {
    let expected_loc =
        compose_pcurve_location(location, edg.location, &brep.locations);
    // OCCT iterates the edge's curve representation list in order.
    for r in edge_representations(edg) {
        let (key, pcurve, range) = match r {
            CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        if key.1 != expected_loc {
            continue;
        }
        if let Some(s) = face_surface_by_ptr(brep, key.0) {
            if surface_same(&s, surface) {
                return Some((pcurve.clone(), range[0], range[1]));
            }
        }
    }
    // The pcurves fast index holds the same rows (some flows populate only
    // the map); the deterministic minimum key among the matching rows is
    // taken (kernel curve_on_surface precedent).
    let ed = match edg.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    let mut best: Option<((u64, u32), (Curve2d, f64, f64))> = None;
    for ((fptr, lhash), (pcurve, f, l)) in &ed.pcurves {
        if *lhash != expected_loc {
            continue;
        }
        if let Some(s) = face_surface_by_ptr(brep, *fptr) {
            if surface_same(&s, surface) {
                let k = (*fptr, *lhash);
                if best.as_ref().map_or(true, |(bk, _)| k < *bk) {
                    best = Some((k, (pcurve.clone(), *f, *l)));
                }
            }
        }
    }
    best.map(|(_, v)| v)
}

/// OCCT BRep_Tool::IsClosed(edge, face) (the seam test): true when the edge
/// carries a BRep_CurveOnClosedSurface representation for the face (two
/// pcurves on one surface).
fn brep_tool_is_closed_edge_face(brep: &BRep, edg: &Shape, fac: &Shape) -> bool {
    let expected_loc =
        compose_pcurve_location(fac.location, edg.location, &brep.locations);
    let fsurf = match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    };
    for r in edge_representations(edg) {
        if let CurveRepresentation::CurveOnClosedSurface { face: (fptr, lhash), .. } = r {
            if *lhash != expected_loc {
                continue;
            }
            // OCCT compares the surface handle; rcad compares the surface
            // value when the pointer row resolves.
            match (face_surface_by_ptr(brep, *fptr), fsurf.as_ref()) {
                (Some(s), Some(fs)) => {
                    if surface_same(&s, fs) {
                        return true;
                    }
                }
                (None, _) => return true,
                _ => {}
            }
        }
    }
    false
}

/// OCCT BRep_Tool::IsClosed(edge, surface, location) - the (surface,
/// location) seam test.
fn brep_tool_is_closed_edge_surface(
    brep: &BRep,
    edg: &Shape,
    surface: &Surface3,
    location: u32,
) -> bool {
    let expected_loc =
        compose_pcurve_location(location, edg.location, &brep.locations);
    for r in edge_representations(edg) {
        if let CurveRepresentation::CurveOnClosedSurface { face: (fptr, lhash), .. } = r {
            if *lhash != expected_loc {
                continue;
            }
            if let Some(s) = face_surface_by_ptr(brep, *fptr) {
                if surface_same(&s, surface) {
                    return true;
                }
            }
        }
    }
    false
}

/// OCCT BRep_Tool::CurveOnPlane(edge, surface, location, f, l)
/// (BRep_Tool.cxx L379-450): for a planar surface, the projection of the
/// edge's 3D curve onto the plane (never stored).  Mirrors the kernel
/// `BRep::curve_on_plane` walk with the surface and location given directly.
fn brep_tool_curve_on_plane_sl(
    brep: &BRep,
    edg: &Shape,
    surface: &Surface3,
    location: u32,
) -> (Option<Curve2d>, f64, f64) {
    // L385: First = Last = 0 (the rcad None return carries no range).
    // L388-398: one Geom_RectangularTrimmedSurface level unwrapped.
    let surf = match surface {
        Surface3::Trimmed(ts) => ts.basis.as_ref(),
        s => s,
    };
    // L400-404: not a plane -> null pcurve.
    let Surface3::Plane(pl) = surf else {
        return (None, 0.0, 0.0);
    };
    // L406-415: BRep_Tool::Curve(E, aCurveLocation, f, l).
    let ed = match edg.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return (None, 0.0, 0.0),
    };
    let Some(c3d) = ed.curve.as_ref() else {
        return (None, 0.0, 0.0);
    };
    // L417: aCurveLocation = L.Predivided(E.Location()).
    let a_curve_location = brep.get_location(location).inverse() * brep.get_location(edg.location);
    let mut f = ed.range[0];
    let mut l = ed.range[1];
    let first = f;
    let last = l;
    // L421-428: transform the curve and update the parameters by the scale
    // factor (TransformedParameter(P, T) = P / T.ScaleFactor(); the scale is
    // recovered as the image length of a unit axis - kernel precedent).
    let c3d = if a_curve_location != glam::DAffine3::IDENTITY {
        let scale = a_curve_location.transform_vector3(DVec3::X).length();
        f /= scale;
        l /= scale;
        transform_curve(c3d, &a_curve_location)
    } else {
        c3d.clone()
    };
    // L430-447: GeomProjLib::ProjectOnPlane + ProjLib_ProjectedCurve + the
    // Geom2d_TrimmedCurve basis unwrap (the landed project_on_plane
    // translation).
    (
        project_on_plane::curve_on_plane(&c3d, [f, l], pl),
        first,
        last,
    )
}

/// OCCT BRep_Tool::Pnt(vtx) - the vertex point located in the parent frame
/// (kernel BRepTool::vertex_position semantics).
fn brep_tool_pnt(brep: &BRep, vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => brep.get_location(vtx.location).transform_point3(vd.point),
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(shape) - vertex/edge/face tolerance.
fn brep_tool_tolerance(_brep: &BRep, the_s: &Shape) -> f64 {
    match the_s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(_brep: &BRep, edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Range(edg, first, last).
fn brep_tool_range(_brep: &BRep, edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

// ---------------------------------------------------------------------------
// TopExp re-hosts (CumOri = true, the OCCT default).
// ---------------------------------------------------------------------------

/// OCCT TopExp::FirstVertex(E, CumOri=true) (TopExp.cxx L182-194): the
/// composed-FORWARD vertex of the edge.
fn top_exp_first_vertex(e: &Shape) -> Shape {
    let (v_first, v_last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    let mut v = if e.orientation == Orientation::Reversed {
        v_last
    } else {
        v_first
    };
    v.orientation = e.orientation.compose(v.orientation);
    v
}

/// OCCT TopExp::LastVertex(E, CumOri=true) (TopExp.cxx L198-210).
fn top_exp_last_vertex(e: &Shape) -> Shape {
    let (v_first, v_last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    let mut v = if e.orientation == Orientation::Reversed {
        v_first
    } else {
        v_last
    };
    v.orientation = e.orientation.compose(v.orientation);
    v
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri=true) (TopExp.cxx): the
/// composed-FORWARD child is Vfirst, the composed-REVERSED child is Vlast.
fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let mut v_first = Shape::null();
    let mut v_last = Shape::null();
    let children: Vec<Shape> = match e.data.as_ref() {
        TShape::Edge(ed) => vec![ed.first.clone(), ed.last.clone()],
        _ => Vec::new(),
    };
    for a_v in &children {
        let mut v = a_v.clone();
        v.orientation = e.orientation.compose(v.orientation);
        match v.orientation {
            Orientation::Forward => v_first = v,
            Orientation::Reversed => v_last = v,
            _ => {}
        }
    }
    (v_first, v_last)
}

/// OCCT TopoDS_Shape::Reverse (TopAbs::Reverse: FORWARD<->REVERSED).
fn topods_shape_reverse(v: &mut Shape) {
    v.orientation = Orientation::Reversed.compose(v.orientation);
}

/// OCCT TopoDS_Shape::IsSame(S): same TShape and same Location.
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

// ---------------------------------------------------------------------------
// GAP carriers (untranslated other-package dependencies).
// ---------------------------------------------------------------------------

// The BRepLib_ValidateEdge GAP carrier was retired: CheckSameParameter now
// consumes the real 1:1 translation
// `crate::topalgo::brep_lib_validate_edge::BRepLibValidateEdge` (with the
// GeomAdaptor_Curve / Geom2dAdaptor_Curve / GeomAdaptor_Surface /
// Adaptor3d_CurveOnSurface re-hosts) at the OCCT L785-798 / L808-817 call
// forms.

/// OCCT BRepExtrema_SupportType (BRepExtrema_DistShapeShape.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BRepExtremaSupportType {
    IsUnknow,
    IsVertex,
    IsOnEdge,
    #[allow(dead_code)]
    IsFace,
}

/// GAP carrier for `BRepExtrema_DistShapeShape` (TKTopAlgo/BRepExtrema,
/// untranslated): the shape-to-shape minimum distance with solutions and
/// supports.  The dependency, the OCCT call anchors and the OCCT failure
/// path (`!IsDone()`) are kept; with the carrier the distance queries always
/// take the OCCT failure branch (IsOverlapPartEdges then reports overlap,
/// and the domain segment walk of CheckOverlapping is skipped) - the exact
/// OCCT behavior when the distance computation fails.  GAP: closes with the
/// BRepExtrema batch (the chfi3d ExtremaExtCC carrier precedent,
/// feat/loc_ope_wires_on_shape_b.rs architecture difference #11).
struct BRepExtremaDistShapeShape {
    #[allow(dead_code)]
    s1: Option<Shape>,
    #[allow(dead_code)]
    s2: Option<Shape>,
    #[allow(dead_code)]
    tolerance: f64,
}

impl BRepExtremaDistShapeShape {
    fn new(s1: &Shape, s2: &Shape, tolerance: f64) -> Self {
        BRepExtremaDistShapeShape {
            s1: Some(s1.clone()),
            s2: Some(s2.clone()),
            tolerance,
        }
    }

    fn new_empty() -> Self {
        BRepExtremaDistShapeShape {
            s1: None,
            s2: None,
            tolerance: 0.0,
        }
    }

    fn load_s1(&mut self, s1: &Shape) {
        self.s1 = Some(s1.clone());
    }

    fn load_s2(&mut self, s2: &Shape) {
        self.s2 = Some(s2.clone());
    }

    fn perform(&mut self) {}

    fn is_done(&self) -> bool {
        false
    }

    #[allow(dead_code)]
    fn value(&self) -> f64 {
        0.0
    }

    fn nb_solution(&self) -> usize {
        0
    }

    fn support_type_shape1(&self, _the_i: usize) -> BRepExtremaSupportType {
        BRepExtremaSupportType::IsUnknow
    }

    fn support_on_shape1(&self, _the_i: usize) -> Shape {
        Shape::null()
    }

    fn par_on_edge_s1(&self, _the_i: usize, _the_param: &mut f64) {}
}

// ---------------------------------------------------------------------------
// The class.
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_Edge (ShapeAnalysis_Edge.hxx L48-256).
pub struct ShapeAnalysisEdge {
    // OCCT `protected: int myStatus` (hxx L243).
    my_status: i32,
}

impl Default for ShapeAnalysisEdge {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisEdge {
    /// OCCT ShapeAnalysis_Edge() (cxx L54-57): initialises Status to OK.
    pub fn new() -> Self {
        ShapeAnalysisEdge {
            my_status: encode_status(ShapeExtendStatus::Ok), // ShapeExtend::EncodeStatus (ShapeExtend_OK)
        }
    }

    /// OCCT HasCurve3d (cxx L92-97): tells if the edge has a 3d curve.
    pub fn has_curve3d(&self, brep: &BRep, edge: &Shape) -> bool {
        // OCCT L94-96: c3d = BRep_Tool::Curve(edge, cf, cl); return !IsNull().
        let c3d = brep_tool_curve(brep, edge).map(|(c, _, _)| c);
        c3d.is_some()
    }

    /// OCCT Curve3d (cxx L101-125): returns the 3d curve and bounding
    /// parameters for the edge; when `orient` the reversed edge toggles
    /// cf/cl.
    pub fn curve3d(
        &self,
        brep: &BRep,
        edge: &Shape,
        c3d: &mut Option<Curve3>,
        cf: &mut f64,
        cl: &mut f64,
        orient: bool,
    ) -> bool {
        // OCCT L107-108: C3d = BRep_Tool::Curve(edge, L, cf, cl).
        let (curve, loc, f, l) = brep_tool_curve_loc(brep, edge);
        *c3d = curve;
        *cf = f;
        *cl = l;
        // OCCT L109-114: apply the location transformation.
        if c3d.is_some() && loc != 0 {
            let trsf = brep.get_location(loc);
            if let Some(c) = c3d.as_ref() {
                *c3d = Some(transform_curve(c, &trsf));
                *cf = c3d.as_ref().unwrap().transformed_parameter(*cf);
                *cl = c3d.as_ref().unwrap().transformed_parameter(*cl);
            }
        }
        if orient {
            if edge.orientation == Orientation::Reversed {
                let tmp = *cf;
                *cf = *cl;
                *cl = tmp;
            }
        }
        c3d.is_some()
    }

    /// OCCT IsClosed3d (cxx L129-142): true when the edge has a 3d curve,
    /// this curve is closed, and the edge has the same vertex at start and
    /// end.
    pub fn is_closed3d(&self, brep: &BRep, edge: &Shape) -> bool {
        let Some((c3d, _cf, _cl)) = brep_tool_curve(brep, edge) else {
            return false;
        };
        if !c3d.is_closed() {
            return false;
        }
        shape_is_same(
            &self.first_vertex(brep, edge),
            &self.last_vertex(brep, edge),
        )
    }

    /// OCCT HasPCurve(edge, face) (cxx L146-151).
    pub fn has_pcurve_face(&self, brep: &BRep, edge: &Shape, face: &Shape) -> bool {
        // OCCT L148-150: L; S = BRep_Tool::Surface(face, L); HasPCurve(edge, S, L).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.has_pcurve_surface(brep, edge, &s, l),
            None => false,
        }
    }

    /// OCCT HasPCurve(edge, surface, location) (cxx L155-167).
    pub fn has_pcurve_surface(
        &self,
        brep: &BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
    ) -> bool {
        // try { //szv#4:S4163:12Mar99 waste try
        let c2d = brep_tool_curve_on_surface_stored(brep, edge, surface, location);
        c2d.is_some()
        /* }
        catch (Standard_Failure) {
        }
        return false; */
    }

    /// OCCT PCurve(edge, face, C2d, cf, cl, orient) (cxx L171-188).
    pub fn pcurve_face(
        &self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
        c2d: &mut Option<Curve2d>,
        cf: &mut f64,
        cl: &mut f64,
        orient: bool,
    ) -> bool {
        //: abv 20.05.02: take into account face orientation
        // COMMENTED BACK - NEEDS MORE CHANGES IN ALL SHAPEHEALING
        //   C2d = BRep_Tool::CurveOnSurface (edge, face, cf, cl);
        //   if (orient && edge.Orientation() == TopAbs_REVERSED) {
        //     double tmp = cf; cf = cl; cl = tmp;
        //   }
        //   return !C2d.IsNull();
        // OCCT L185-187: L; S = BRep_Tool::Surface(face, L); PCurve(edge, S, L, ...).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.pcurve_surface(brep, edge, &s, l, c2d, cf, cl, orient),
            None => {
                *c2d = None;
                false
            }
        }
    }

    /// OCCT PCurve(edge, surface, location, C2d, cf, cl, orient)
    /// (cxx L192-208).
    pub fn pcurve_surface(
        &self,
        brep: &BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
        c2d: &mut Option<Curve2d>,
        cf: &mut f64,
        cl: &mut f64,
        orient: bool,
    ) -> bool {
        match brep_tool_curve_on_surface_stored(brep, edge, surface, location) {
            Some((pc, f, l)) => {
                *c2d = Some(pc);
                *cf = f;
                *cl = l;
            }
            None => {
                *c2d = None;
            }
        }
        if orient && edge.orientation == Orientation::Reversed {
            let tmp = *cf;
            *cf = *cl;
            *cl = tmp;
        }
        c2d.is_some()
    }

    /// OCCT BoundUV(edge, face, first, last) (cxx L61-69).
    pub fn bound_uv_face(
        &self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
        first: &mut DVec2,
        last: &mut DVec2,
    ) -> bool {
        // OCCT L66-68: L; S = BRep_Tool::Surface(face, L); BoundUV(edge, S, L, ...).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.bound_uv_surface(brep, edge, &s, l, first, last),
            None => false,
        }
    }

    /// OCCT BoundUV(edge, surface, location, first, last) (cxx L73-88):
    /// returns the ends of pcurve (calls PCurve with orient = True).
    pub fn bound_uv_surface(
        &self,
        brep: &BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
        first: &mut DVec2,
        last: &mut DVec2,
    ) -> bool {
        let mut c2d: Option<Curve2d> = None;
        let mut uf = 0.0;
        let mut ul = 0.0;
        if !self.pcurve_surface(brep, edge, surface, location, &mut c2d, &mut uf, &mut ul, true) {
            return false;
        }
        match c2d.as_ref() {
            Some(c) => {
                *first = c.point_at(uf);
                *last = c.point_at(ul);
            }
            None => return false,
        }
        true
    }

    /// OCCT IsSeam(edge, face) (cxx L212-215).
    pub fn is_seam_face(&self, brep: &BRep, edge: &Shape, face: &Shape) -> bool {
        brep_tool_is_closed_edge_face(brep, edge, face)
    }

    /// OCCT IsSeam(edge, surface, location) (cxx L219-224): true when the
    /// edge has two pcurves on one surface.
    pub fn is_seam_surface(
        &self,
        brep: &BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
    ) -> bool {
        brep_tool_is_closed_edge_surface(brep, edge, surface, location)
    }

    /// OCCT FirstVertex (cxx L228-241): start vertex of the edge (taking
    /// edge orientation into account).
    pub fn first_vertex(&self, brep: &BRep, edge: &Shape) -> Shape {
        let _ = brep;
        let mut v;
        if edge.orientation == Orientation::Reversed {
            v = top_exp_last_vertex(edge);
            topods_shape_reverse(&mut v);
        } else {
            v = top_exp_first_vertex(edge);
        }
        v
    }

    /// OCCT LastVertex (cxx L245-258): end vertex of the edge (taking edge
    /// orientation into account).
    pub fn last_vertex(&self, brep: &BRep, edge: &Shape) -> Shape {
        let _ = brep;
        let mut v;
        if edge.orientation == Orientation::Reversed {
            v = top_exp_first_vertex(edge);
            topods_shape_reverse(&mut v);
        } else {
            v = top_exp_last_vertex(edge);
        }
        v
    }

    /// OCCT Status (cxx L262-265): the status (True/False) of last Check.
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    /// OCCT GetEndTangent2d(edge, face, atEnd, pos, tang, dparam)
    /// (cxx L269-279).
    pub fn get_end_tangent2d_face(
        &self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
        atend1: bool, /* skl : change "atend" to "atend1" */
        pnt: &mut DVec2,
        v: &mut DVec2,
        dparam: f64,
    ) -> bool {
        // OCCT L276-278: L; S = BRep_Tool::Surface(face, L); GetEndTangent2d(edge, S, L, ...).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.get_end_tangent2d_surface(brep, edge, &s, l, atend1, pnt, v, dparam),
            None => {
                *v = DVec2::ZERO;
                false
            }
        }
    }

    /// OCCT GetEndTangent2d(edge, surface, location, atEnd, pos, tang,
    /// dparam) (cxx L283-366): tangent of the edge pcurve at its start
    /// (atend2 = False) or end (True), regarding the orientation of edge.
    /// If edge is REVERSED, tangent is reversed before return.  Returns True
    /// if pcurve is available and tangent is computed and is not null.
    pub fn get_end_tangent2d_surface(
        &self,
        brep: &BRep,
        edge: &Shape,
        s: &Surface3,
        l: u32,
        atend2: bool, /* skl : change "atend" to "atend2" */
        pnt: &mut DVec2,
        v: &mut DVec2,
        dparam: f64,
    ) -> bool {
        let mut cf = 0.0;
        let mut cl = 0.0;
        let mut c2d_opt: Option<Curve2d> = None;
        if !self.pcurve_surface(brep, edge, s, l, &mut c2d_opt, &mut cf, &mut cl, true) {
            *v = DVec2::ZERO;
            return false;
        }
        let c2d = match c2d_opt {
            Some(c) => c,
            None => return false,
        };
        let mut dpnew = dparam;

        if dpnew > CONFUSION {
            let ptmp;
            let par1;
            let par2;
            let delta = (cl - cf) * dpnew;
            if delta.abs() < PCONFUSION {
                dpnew = 0.0;
            } else {
                if atend2 {
                    par1 = cl;
                    par2 = cl - delta;
                    *pnt = c2d.point_at(par1);
                    ptmp = c2d.point_at(par2);
                    *v = *pnt - ptmp;
                } else {
                    par1 = cf;
                    par2 = cf + delta;
                    *pnt = c2d.point_at(par1);
                    ptmp = c2d.point_at(par2);
                    *v = ptmp - *pnt;
                }
                if v.length_squared() < PCONFUSION * PCONFUSION {
                    dpnew = 0.0;
                }
            }
        }

        if dpnew <= CONFUSION {
            // get non-null tangency searching until 3rd derivative, or as straight btw ends
            let par = if atend2 { cl } else { cf };
            // OCCT L337: c2d->D1(par, pnt, v).
            *pnt = c2d.point_at(par);
            *v = c2d.derivative_at(par);
            if v.length_squared() < PCONFUSION * PCONFUSION {
                // OCCT L340-341: gp_Vec2d d1; c2d->D2(par, pnt, d1, v).
                *pnt = c2d.point_at(par);
                let _d1 = c2d.derivative_at(par);
                *v = c2d.derivative2_at(par);
                if v.length_squared() < PCONFUSION * PCONFUSION {
                    // OCCT L344-345: gp_Vec2d d2; c2d->D3(par, pnt, d1, d2, v).
                    *pnt = c2d.point_at(par);
                    let _d1 = _d1;
                    let _d2 = c2d.derivative2_at(par);
                    *v = c2d.derivative3_at(par);
                    if v.length_squared() < PCONFUSION * PCONFUSION {
                        // OCCT L348-350: the straight-between-ends fallback.
                        *pnt = c2d.point_at(par);
                        let p2 = c2d.point_at(if atend2 { cf } else { cl });
                        *v = p2 - *pnt;
                        if v.length_squared() < PCONFUSION * PCONFUSION {
                            return false;
                        }
                    }
                }
            }
            if edge.orientation == Orientation::Reversed {
                *v = -*v;
            }
        }

        // if ( edge.Orientation() == TopAbs_REVERSED ) v.Reverse();
        true
    }

    /// OCCT CheckCurve3dWithPCurve(edge, face) (cxx L370-375).
    pub fn check_curve3d_with_pcurve_face(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
    ) -> bool {
        // OCCT L372-374: L; S = BRep_Tool::Surface(face, L); the (S, L) call.
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.check_curve3d_with_pcurve_surface(brep, edge, &s, l),
            None => false,
        }
    }

    /// OCCT CheckCurve3dWithPCurve(edge, surface, location) (cxx L379-425):
    /// checks mutual orientation of 3d curve and pcurve on the analysis of
    /// curves bounding points.
    pub fn check_curve3d_with_pcurve_surface(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        surface: &Surface3,
        location: u32,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        // OCCT L385: surface->IsKind(STANDARD_TYPE(Geom_Plane)).
        if matches!(surface, Surface3::Plane(_)) {
            return false;
        }

        let mut c2d_opt: Option<Curve2d> = None;
        let mut f2d = 0.0;
        let mut l2d = 0.0; // szv#4:S4163:12Mar99 moved down f3d, l3d
        if !self.pcurve_surface(brep, edge, surface, location, &mut c2d_opt, &mut f2d, &mut l2d, false)
        {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let c2d = match c2d_opt {
            Some(c) => c,
            None => return false,
        };

        let mut c3d_opt: Option<Curve3> = None; // szv#4:S4163:12Mar99 moved
        let mut f3d = 0.0;
        let mut l3d = 0.0; // szv#4:S4163:12Mar99 moved
        if !self.curve3d(brep, edge, &mut c3d_opt, &mut f3d, &mut l3d, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            return false;
        }
        let c3d = match c3d_opt {
            Some(c) => c,
            None => return false,
        };

        let a_first_vert = self.first_vertex(brep, edge);
        let a_last_vert = self.last_vertex(brep, edge);

        if a_first_vert.is_null() || a_last_vert.is_null() {
            return false;
        }

        let preci1 = brep_tool_tolerance(brep, &a_first_vert);
        let preci2 = brep_tool_tolerance(brep, &a_last_vert);

        let p2d1 = c2d.point_at(f2d);
        let p2d2 = c2d.point_at(l2d);

        // #39 rln 17.11.98 S4054, annie_surf.igs entity 39
        let loc_trsf = brep.get_location(location);
        self.check_points(
            c3d.point_at(f3d), /*.Transformed (location.Transformation())*/
            c3d.point_at(l3d), /*.Transformed (location.Transformation())*/
            loc_trsf.transform_point3(surface.point_at(p2d1.x, p2d1.y)),
            loc_trsf.transform_point3(surface.point_at(p2d2.x, p2d2.y)),
            preci1,
            preci2,
        )
    }

    /// OCCT CheckPoints (cxx L429-446): check points by pairs (A and A, B and
    /// B) with precisions (preci1 and preci2).  P1 are the points either from
    /// 3d curve or from vertices, P2 are the points from pcurve.
    fn check_points(
        &mut self,
        p1a: DVec3,
        p1b: DVec3,
        p2a: DVec3,
        p2b: DVec3,
        preci1: f64,
        preci2: f64,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if p1a.distance_squared(p2a) <= preci1 * preci1
            && p1b.distance_squared(p2b) <= preci2 * preci2
        {
            return false;
        } else if p1a.distance(p2b) + (p1b.distance(p2a)) < p1a.distance(p2a) + (p1b.distance(p2b))
        {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
        }
        true
    }

    /// OCCT CheckVerticesWithCurve3d (cxx L450-493): checks the start and/or
    /// end vertex of the edge for matching with 3d curve with the given
    /// precision.  vtx = 0: both (default), 1: start only, 2: end only.
    /// If preci < 0 the vertices are considered with their own tolerances.
    pub fn check_vertices_with_curve3d(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        preci: f64,
        vtx: i32,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        let v1 = self.first_vertex(brep, edge);
        let v2 = self.last_vertex(brep, edge);
        let p1v = brep_tool_pnt(brep, &v1);
        let p2v = brep_tool_pnt(brep, &v2);

        let mut cf = 0.0;
        let mut cl = 0.0;
        let mut c3d_opt: Option<Curve3> = None;
        if !self.curve3d(brep, edge, &mut c3d_opt, &mut cf, &mut cl, true) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let c3d = match c3d_opt {
            Some(c) => c,
            None => return false,
        };

        //  on va faire les checks ...
        if vtx != 2 {
            //  1er VTX
            let p13d = c3d.point_at(cf);
            // szv#4:S4163:12Mar99 optimized
            if p1v.distance(p13d)
                > (if preci < 0.0 {
                    brep_tool_tolerance(brep, &v1)
                } else {
                    preci
                })
            {
                self.my_status |= encode_status(ShapeExtendStatus::Done1);
            }
        }

        if vtx != 1 {
            //  2me VTX
            let p23d = c3d.point_at(cl);
            // szv#4:S4163:12Mar99 optimized
            if p2v.distance(p23d)
                > (if preci < 0.0 {
                    brep_tool_tolerance(brep, &v2)
                } else {
                    preci
                })
            {
                self.my_status |= encode_status(ShapeExtendStatus::Done2);
            }
        }

        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckVerticesWithPCurve(edge, face, preci, vtx) (cxx L497-507).
    pub fn check_vertices_with_pcurve_face(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
        preci: f64,
        vtx: i32,
    ) -> bool {
        // OCCT L502-506: L; S = BRep_Tool::Surface(face, L); the (S, L) call
        // (szv#4:S4163:12Mar99 `vtx,preci` wrong parameters order fixed).
        let (s, l) = brep_tool_surface_loc(brep, face);
        match s {
            Some(s) => self.check_vertices_with_pcurve_surface(brep, edge, &s, l, preci, vtx),
            None => false,
        }
    }

    /// OCCT CheckVerticesWithPCurve(edge, surface, location, preci, vtx)
    /// (cxx L511-564).
    pub fn check_vertices_with_pcurve_surface(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        surf: &Surface3,
        loc: u32,
        preci: f64,
        vtx: i32,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        let v1 = self.first_vertex(brep, edge);
        let v2 = self.last_vertex(brep, edge);
        let p1v = brep_tool_pnt(brep, &v1);
        let p2v = brep_tool_pnt(brep, &v2);

        let mut cf = 0.0;
        let mut cl = 0.0;
        let mut c2d_opt: Option<Curve2d> = None;
        if !self.pcurve_surface(brep, edge, surf, loc, &mut c2d_opt, &mut cf, &mut cl, true) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let c2d = match c2d_opt {
            Some(c) => c,
            None => return false,
        };

        // on va faire les checks ...
        if vtx != 2 {
            //  1er VTX
            let p1uv = c2d.point_at(cf);
            let mut p12d = surf.point_at(p1uv.x, p1uv.y);
            if loc != 0 {
                p12d = brep.get_location(loc).transform_point3(p12d);
            }
            // szv#4:S4163:12Mar99 optimized
            if p1v.distance(p12d)
                > (if preci < 0.0 {
                    brep_tool_tolerance(brep, &v1)
                } else {
                    preci
                })
            {
                self.my_status |= encode_status(ShapeExtendStatus::Done1);
            }
        }

        if vtx != 1 {
            //  2me VTX
            let p2uv = c2d.point_at(cl);
            let mut p22d = surf.point_at(p2uv.x, p2uv.y);
            if loc != 0 {
                p22d = brep.get_location(loc).transform_point3(p22d);
            }
            // szv#4:S4163:12Mar99 optimized
            if p2v.distance(p22d)
                > (if preci < 0.0 {
                    brep_tool_tolerance(brep, &v2)
                } else {
                    preci
                })
            {
                self.my_status |= encode_status(ShapeExtendStatus::Done2);
            }
        }

        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT static ::CheckVertexTolerance (cxx L568-668): checks if it is
    /// necessary to increase tolerances of the edge vertices to comprise the
    /// ends of 3d curve and pcurves.  `face` = None + `check_all` = true is
    /// the all-stored-pcurves variant (the OCCT null-face call).
    fn check_vertex_tolerance_impl(
        &self,
        brep: &BRep,
        edge: &Shape,
        face: Option<&Shape>,
        check_all: bool,
        toler1: &mut f64,
        toler2: &mut f64,
    ) -> i32 {
        let mut status = encode_status(ShapeExtendStatus::Ok);

        let v1 = self.first_vertex(brep, edge);
        let v2 = self.last_vertex(brep, edge);
        if v1.is_null() || v2.is_null() {
            //: p1 abv 22 Feb 99: r76sy.stp
            status |= encode_status(ShapeExtendStatus::Fail1);
            return status;
        }

        let old1 = brep_tool_tolerance(brep, &v1);
        let old2 = brep_tool_tolerance(brep, &v2);
        let pnt1 = brep_tool_pnt(brep, &v1);
        let pnt2 = brep_tool_pnt(brep, &v2);

        let mut a = 0.0;
        let mut b = 0.0;
        let mut c3d_opt: Option<Curve3> = None;
        if !self.curve3d(brep, edge, &mut c3d_opt, &mut a, &mut b, true) {
            if !brep_tool_degenerated(brep, edge) {
                status |= encode_status(ShapeExtendStatus::Fail2);
            }
            *toler1 = 0.;
            *toler2 = 0.;
            //    return false;
        } else if let Some(c3d) = c3d_opt.as_ref() {
            *toler1 = pnt1.distance_squared(c3d.point_at(a));
            *toler2 = pnt2.distance_squared(c3d.point_at(b));
        }

        if check_all {
            // OCCT L609-629: the BRep_TEdge curve-representation walk.
            for itcr in edge_representations(edge) {
                // OCCT: GC = down_cast<BRep_GCurve>(itcr.Value());
                //       if (GC.IsNull() || !GC->IsCurveOnSurface()) continue;
                let (pcurve, range, sfptr) = match itcr {
                    CurveRepresentation::CurveOnSurface {
                        face: (fptr, _),
                        pcurve,
                        range,
                    } => (pcurve, range, *fptr),
                    CurveRepresentation::CurveOnClosedSurface {
                        face: (fptr, _),
                        pcurve1,
                        range,
                        ..
                    } => (pcurve1, range, *fptr),
                    _ => continue,
                };
                // OCCT: S = GC->Surface(); L = edge.Location() * GC->Location()
                // (bridge #5: the row location hash is one-way; the transform
                // below uses the edge wrapper location).
                let Some(s) = face_surface_by_ptr(brep, sfptr) else {
                    continue;
                };
                a = range[0];
                b = range[1];
                // OCCT L622: sae.PCurve(edge, S, L, pcurve, a, b, true) - the
                // walk already holds the row the lookup would return.
                let p1 = pcurve.point_at(a);
                let p2 = pcurve.point_at(b);
                let loc_trsf = brep.get_location(edge.location);
                let p1w = loc_trsf.transform_point3(s.point_at(p1.x, p1.y));
                let p2w = loc_trsf.transform_point3(s.point_at(p2.x, p2.y));
                *toler1 = toler1.max(pnt1.distance_squared(p1w));
                *toler2 = toler2.max(pnt2.distance_squared(p2w));
            }
        }
        //: abv 10.06.02: porting C40 -> dev (CC670-12608.stp)
        // Check with given face is needed for plane surfaces (if no stored pcurves)
        else if let Some(face) = face {
            // OCCT L635-650.
            let (s, l) = brep_tool_surface_loc(brep, face);
            let Some(s) = s else {
                return status;
            };
            let mut pcurve: Option<Curve2d> = None;
            let mut fa = 0.0;
            let mut fb = 0.0;
            if self.pcurve_surface(brep, edge, &s, l, &mut pcurve, &mut fa, &mut fb, true) {
                if let Some(pc) = pcurve.as_ref() {
                    let p1 = pc.point_at(fa);
                    let p2 = pc.point_at(fb);
                    let loc_trsf = brep.get_location(l);
                    let p1w = loc_trsf.transform_point3(s.point_at(p1.x, p1.y));
                    let p2w = loc_trsf.transform_point3(s.point_at(p2.x, p2.y));
                    *toler1 = toler1.max(pnt1.distance_squared(p1w));
                    *toler2 = toler2.max(pnt2.distance_squared(p2w));
                }
            } else {
                status |= encode_status(ShapeExtendStatus::Fail3);
            }
        }

        //: o8 abv 19 Feb 99: CTS18541.stp #18559: coeff 1.0001 added
        // szv 18 Aug 99: edge tolerance is taken in consideration
        let tole = brep_tool_tolerance(brep, edge);
        *toler1 = (1.0000001 * toler1.sqrt()).max(tole);
        *toler2 = (1.0000001 * toler2.sqrt()).max(tole);
        if *toler1 > old1 {
            status |= encode_status(ShapeExtendStatus::Done1);
        }
        if *toler2 > old2 {
            status |= encode_status(ShapeExtendStatus::Done2);
        }

        status
    }

    /// OCCT CheckVertexTolerance(edge, face, toler1, toler2) (cxx L672-679).
    pub fn check_vertex_tolerance_face(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
        toler1: &mut f64,
        toler2: &mut f64,
    ) -> bool {
        self.my_status =
            self.check_vertex_tolerance_impl(brep, edge, Some(face), false, toler1, toler2);
        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckVertexTolerance(edge, toler1, toler2) (cxx L683-690): the
    /// all-stored-pcurves variant (OCCT passes a null face).
    pub fn check_vertex_tolerance_all(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        toler1: &mut f64,
        toler2: &mut f64,
    ) -> bool {
        self.my_status = self.check_vertex_tolerance_impl(brep, edge, None, true, toler1, toler2);
        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckSameParameter(edge, maxdev, NbControl) (cxx L694-700).
    pub fn check_same_parameter(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        maxdev: &mut f64,
        nb_control: i32,
    ) -> bool {
        // OCCT L698-699: TopoDS_Face anEmptyFace; the face variant.
        self.check_same_parameter_face(brep, edge, &Shape::null(), maxdev, nb_control)
    }

    /// OCCT CheckSameParameter(edge, face, maxdev, NbControl) (cxx L704-838):
    /// checks the edge to be SameParameter; calculates the maximal deviation
    /// between 3d curve and each pcurve of the edge on NbControl equidistant
    /// points (the same algorithm as in BRepCheck; default value is 23).
    /// If deviation is greater than tolerance of the edge (i.e. incorrect
    /// flag) returns False, else returns True.
    pub fn check_same_parameter_face(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        face: &Shape,
        maxdev: &mut f64,
        nb_control: i32,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if brep_tool_degenerated(brep, edge) {
            return false;
        }

        *maxdev = 0.0;

        // Get same parameter flag (OCCT L718-719: TE->SameParameter()).
        let same_parameter = match edge.data.as_ref() {
            TShape::Edge(ed) => ed.same_parameter,
            _ => false,
        };

        // Get 3D curve of the edge (OCCT L722-729).
        let (c3d0, a_curve_loc, mut a_first, mut a_last) = brep_tool_curve_loc(brep, edge);
        if c3d0.is_none() {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let mut a_c3d = c3d0.unwrap();

        if a_curve_loc != 0 {
            let trsf = brep.get_location(a_curve_loc);
            a_c3d = transform_curve(&a_c3d, &trsf);
            a_first = a_c3d.transformed_parameter(a_first);
            a_last = a_c3d.transformed_parameter(a_last);
        }

        // Create adaptor for the curve (OCCT L740: GeomAdaptor_Curve(aC3D,
        // aFirst, aLast)).
        let a_gac = GeomAdaptorCurve::new(a_c3d.clone(), a_first, a_last);

        // OCCT L742-747: the face surface + location.
        let (a_face_surf, a_face_loc) = if !face.is_null() {
            brep_tool_surface_loc(brep, face)
        } else {
            (None, 0)
        };

        let mut is_pcurve_found = false;
        let mut i = 1i32;

        // Iterate on all curve representations (OCCT L753: for(;;)).
        loop {
            // OCCT L760: BRep_Tool::CurveOnSurface(edge, aPC, aS, aLoc, f, l, i)
            // (the index overload: the i-th curve-on-surface representation;
            // aLoc is its full representation location - bridge #5 keeps the
            // pcurve-key hash).
            let (mut a_pc, mut a_s, mut a_loc, mut f, mut l) = (None, None, 0u32, 0.0, 0.0);
            let mut k = 0i32;
            for r in edge_representations(edge) {
                let (key, pc, range) = match r {
                    CurveRepresentation::CurveOnSurface {
                        face,
                        pcurve,
                        range,
                    } => (*face, pcurve, range),
                    CurveRepresentation::CurveOnClosedSurface {
                        face,
                        pcurve1,
                        range,
                        ..
                    } => (*face, pcurve1, range),
                    _ => continue,
                };
                k += 1;
                if k == i {
                    a_pc = Some(pc.clone());
                    a_s = face_surface_by_ptr(brep, key.0);
                    a_loc = key.1;
                    f = range[0];
                    l = range[1];
                    break;
                }
            }
            if k < i {
                a_pc = None;
            }

            let Some(a_pc) = a_pc else {
                // No more curves (OCCT L762-765).
                break;
            };

            i += 1;

            // If the input face is not null, check that the curve is on its
            // surface (OCCT L770-777: aFaceSurf != aS || aFaceLoc != aLoc).
            if let Some(fs) = a_face_surf.as_ref() {
                let on_face = match a_s.as_ref() {
                    Some(a_s) => {
                        surface_same(a_s, fs)
                            && a_loc
                                == compose_pcurve_location(
                                    a_face_loc,
                                    edge.location,
                                    &brep.locations,
                                )
                    }
                    None => false,
                };
                if !on_face {
                    continue;
                }
            }

            is_pcurve_found = true;

            // Apply transformations for the surface (OCCT L781-783;
            // bridge #5: the edge wrapper location).
            let a_st = match a_s.as_ref() {
                Some(a_s) => {
                    if edge.location != 0 {
                        transform_surface(a_s, &brep.get_location(edge.location))
                    } else {
                        a_s.clone()
                    }
                }
                None => continue,
            };

            // Compute deviation between curves (OCCT L785-789: the
            // Geom2dAdaptor_Curve + GeomAdaptor_Surface +
            // Adaptor3d_CurveOnSurface construction).
            let a_ghpc = Geom2dAdaptorCurve::new(a_pc, f, l);
            let a_gahs = GeomAdaptorSurface::new(a_st);
            let a_acs = Adaptor3dCurveOnSurface::new(a_ghpc, a_gahs);

            let mut a_validate_edge =
                BRepLibValidateEdge::new(a_gac.clone(), a_acs, same_parameter);
            a_validate_edge.set_control_points_number(nb_control - 1);
            a_validate_edge.process();
            a_validate_edge.update_tolerance(maxdev);
            if !a_validate_edge.is_done() {
                self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            }
        }

        // For the planar face and non-existing 2d curve
        // check the deviation for the projection of the 3d curve on plane
        // (OCCT L801-826).
        if !is_pcurve_found {
            if let Some(a_face_surf) = a_face_surf.as_ref() {
                let (a_pc_opt, _f, _l) =
                    brep_tool_curve_on_plane_sl(brep, edge, a_face_surf, a_face_loc);
                if let Some(a_pc) = a_pc_opt {
                    let a_st = if a_face_loc != 0 {
                        transform_surface(a_face_surf, &brep.get_location(a_face_loc))
                    } else {
                        a_face_surf.clone()
                    };
                    // OCCT L808-816: the Geom2dAdaptor_Curve +
                    // GeomAdaptor_Surface + Adaptor3d_CurveOnSurface
                    // construction.
                    let a_ghpc = Geom2dAdaptorCurve::new(a_pc, a_first, a_last);
                    let a_gahs = GeomAdaptorSurface::new(a_st);
                    let a_acs = Adaptor3dCurveOnSurface::new(a_ghpc, a_gahs);

                    let mut a_validate_edge_on_plane =
                        BRepLibValidateEdge::new(a_gac.clone(), a_acs, same_parameter);
                    a_validate_edge_on_plane.set_control_points_number(nb_control - 1);
                    a_validate_edge_on_plane.process();
                    a_validate_edge_on_plane.update_tolerance(maxdev);
                    if !a_validate_edge_on_plane.is_done() {
                        self.my_status |= encode_status(ShapeExtendStatus::Fail2);
                    }
                }
            }
        }

        // OCCT L828-835: the tolerance comparisons.
        let te_tolerance = match edge.data.as_ref() {
            TShape::Edge(ed) => ed.tolerance,
            _ => 0.0,
        };
        if *maxdev > te_tolerance {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
        }
        if !same_parameter {
            self.my_status |= encode_status(ShapeExtendStatus::Done2);
        }

        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckPCurveRange (cxx L999-1033): checks possibility for pcurve
    /// thePC to have range [theFirst, theLast] (edge range) having respect to
    /// real first, last parameters of thePC.
    pub fn check_pcurve_range(&self, the_first: f64, the_last: f64, the_pc: &Curve2d) -> bool {
        let eps = PCONFUSION;
        let mut is_valid = true;
        let mut is_periodic = the_pc.is_periodic();
        let mut a_period = REAL_LAST;
        if is_periodic {
            a_period = geom2d_curve_period(the_pc);
        }
        // OCCT L1011: fp = thePC->FirstParameter(); lp = thePC->LastParameter().
        let (mut fp, mut lp) = match the_pc {
            Curve2d::Trimmed(tc) => (tc.t_min, tc.t_max),
            _ => {
                let domain = the_pc.default_domain();
                (domain[0], domain[1])
            }
        };
        // OCCT L1012-1022: the Geom2d_TrimmedCurve basis unwrap.
        if let Curve2d::Trimmed(tc) = the_pc {
            let a_c = tc.curve.as_ref();
            let domain = a_c.default_domain();
            fp = domain[0];
            lp = domain[1];
            is_periodic = a_c.is_periodic();
            if is_periodic {
                a_period = geom2d_curve_period(a_c);
            }
        }
        if is_periodic && (the_last - the_first > a_period + eps) {
            is_valid = false;
        } else if !is_periodic && (the_first < fp - eps || the_last > lp + eps) {
            is_valid = false;
        }

        is_valid
    }

    /// OCCT CheckOverlapping (cxx L894-995): checks the first edge is
    /// overlapped with second edge.  If distance between two edges is less
    /// then theTolOverlap edges are overlapped.  theDomainDist - length of
    /// part of edges on which edges are overlapped.
    pub fn check_overlapping(
        &mut self,
        brep: &BRep,
        the_edge1: &Shape,
        the_edge2: &Shape,
        the_tol_overlap: &mut f64,
        the_domain_dist: f64,
    ) -> bool {
        let mut is_overlap = false;
        // OCCT L900-903: BRepAdaptor_Curve + GCPnts_AbscissaPoint::Length.
        let Some((a_ad_curve1, a_c1_first, a_c1_last)) = brep_tool_curve(brep, the_edge1) else {
            return false;
        };
        let a_length1 = arc_length(&a_ad_curve1, a_c1_first, a_c1_last).abs();
        let Some((a_ad_curve2, a_c2_first, a_c2_last)) = brep_tool_curve(brep, the_edge2) else {
            return false;
        };
        let a_length2 = arc_length(&a_ad_curve2, a_c2_first, a_c2_last).abs();
        let a_first_edge = if a_length1 >= a_length2 {
            the_edge2.clone()
        } else {
            the_edge1.clone()
        };
        let a_sec_edge = if a_length1 >= a_length2 {
            the_edge1.clone()
        } else {
            the_edge2.clone()
        };
        let a_length = a_length1.min(a_length2);

        // check overalpping between edges on whole edges
        let a_step = a_length1.min(a_length2) / 2.0;
        is_overlap = is_overlap_part_edges(
            brep,
            &a_first_edge,
            &a_sec_edge,
            *the_tol_overlap,
            a_step,
            0.,
            a_length1.min(a_length2),
        );

        if is_overlap {
            self.my_status |= encode_status(ShapeExtendStatus::Done3);
            return is_overlap;
        }
        if the_domain_dist == 0.0 {
            return is_overlap;
        }

        // check overalpping between edges on segment with length less than theDomainDist

        let a_domain_tol = if the_domain_dist > a_length1.min(a_length2) {
            a_length1.min(a_length2)
        } else {
            the_domain_dist
        };
        let a_min_dist =
            BRepExtremaDistShapeShape::new(&a_first_edge, &a_sec_edge, *the_tol_overlap);
        let mut ares_tol = *the_tol_overlap;
        if a_min_dist.is_done() {
            ares_tol = a_min_dist.value();
            if ares_tol >= *the_tol_overlap {
                return false;
            }
            let nb_sol = a_min_dist.nb_solution();
            let mut i = 1usize;
            while i <= nb_sol && !is_overlap {
                let a_type1 = a_min_dist.support_type_shape1(i);
                let a_length_p;
                if a_type1 == BRepExtremaSupportType::IsVertex {
                    let a_support_shape1 = a_min_dist.support_on_shape1(i);
                    let (a_v1, _a_v2) = top_exp_vertices(&a_first_edge);
                    if shape_is_same(&a_v1, &a_support_shape1) {
                        a_length_p = 0.0;
                    } else {
                        a_length_p = a_length;
                    }
                } else if a_type1 == BRepExtremaSupportType::IsOnEdge {
                    let mut a_param1 = 0.0;
                    a_min_dist.par_on_edge_s1(i, &mut a_param1);
                    // OCCT L963: BRep_Tool::Range(aFirstEdge, aFirst, aLast).
                    let (_a_first, _a_last) = brep_tool_range(brep, &a_first_edge);
                    // OCCT L964-965: BRepAdaptor_Curve +
                    // GCPnts_AbscissaPoint::Length(anAdaptor, aFirst, aParam1).
                    match brep_tool_curve(brep, &a_first_edge) {
                        Some((curve, first, _last)) => {
                            a_length_p = arc_length(&curve, first, a_param1).abs();
                        }
                        None => a_length_p = a_length,
                    }
                } else {
                    i += 1;
                    continue;
                }
                let mut a_start_length = a_length_p - a_domain_tol / 2.0;
                if a_start_length < 0.0 {
                    a_start_length = 0.0;
                    // OCCT L975: aEndLength = aDomainTol; (overwritten below).
                }
                let mut a_end_length = a_length_p + a_domain_tol / 2.0;
                if a_end_length > a_length {
                    a_end_length = a_length;
                    a_start_length = a_end_length - a_domain_tol;
                }
                let a_step = (a_end_length - a_start_length) / 5.0;
                is_overlap = is_overlap_part_edges(
                    brep,
                    &a_first_edge,
                    &a_sec_edge,
                    *the_tol_overlap,
                    a_step,
                    a_start_length,
                    a_end_length,
                );
                i += 1;
            }
        }
        if is_overlap {
            self.my_status |= encode_status(ShapeExtendStatus::Done4);
        }

        *the_tol_overlap = ares_tol;
        is_overlap
    }
}

/// OCCT static IsOverlapPartEdges (cxx L842-890): checks that every sampled
/// point of the first edge lies within theTolerance of the second edge.
fn is_overlap_part_edges(
    brep: &BRep,
    the_first_edge: &Shape,
    the_sec_edge: &Shape,
    the_tolerance: f64,
    the_step: f64,
    the_start_length: f64,
    the_end_length: f64,
) -> bool {
    // OCCT L849: NCollection_Sequence<int> aSeqIntervals; - declared, never
    // used in the walk (kept as this comment).
    // OCCT L850: BRepAdaptor_Curve aAdCurve1(theFirstEdge).
    let Some((a_ad_curve1, a_first_param, a_last_param)) = brep_tool_curve(brep, the_first_edge)
    else {
        // OCCT: BRepAdaptor_Curve without a 3d curve raises
        // Standard_NoSuchObject; the rcad walk treats it as "cannot disprove
        // the overlap" (the OCCT failure path of the distance loop).
        return true;
    };

    let mut a_min_dist = BRepExtremaDistShapeShape::new_empty();
    a_min_dist.load_s1(the_sec_edge);

    let mut a_s = the_start_length;
    while a_s <= the_end_length {
        let a_point: DVec3;
        if a_s <= CONFUSION {
            // OCCT L861-862: TopExp::FirstVertex(theFirstEdge, true).
            let v1 = top_exp_first_vertex(the_first_edge);
            a_point = brep_tool_pnt(brep, &v1);
        } else {
            // OCCT L866-869: GCPnts_AbscissaPoint(Precision::Confusion(),
            // aAdCurve1, aS, aAdCurve1.FirstParameter()); the rcad GCPnts
            // port always yields a parameter (IsDone() = true).
            let a_abs_point_param = abscissa_point_parameter(
                &a_ad_curve1,
                a_first_param,
                a_last_param,
                a_s,
                a_first_param + a_s,
            );
            // OCCT L872: aAdCurve1.D0(aAbsPoint.Parameter(), aPoint).
            a_point = a_ad_curve1.point_at(a_abs_point_param);
        }
        // OCCT L879-881: BRep_Builder aB; aB.MakeVertex(aV, aPoint,
        // Precision::Confusion()) on a local builder - the rcad vertex is a
        // detached TShape (never enters the shape pool).
        let a_v = Shape::new(
            Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                point: a_point,
                tolerance: CONFUSION,
                points: Vec::new(),
            })),
            0,
            Orientation::Forward,
        );
        a_min_dist.load_s2(&a_v);
        a_min_dist.perform();
        if a_min_dist.is_done() && a_min_dist.value() >= the_tolerance {
            return false;
        }
        a_s += the_step / 2.0;
    }
    true
}

/// OCCT Geom2d_Curve::Period() - the period of the periodic curve (Circle
/// and Ellipse are 2*pi; the rcad BSpline carrier stores no periodic flag,
/// matching its IsPeriodic() = false).
fn geom2d_curve_period(pc: &Curve2d) -> f64 {
    match pc {
        Curve2d::Circle(_) => std::f64::consts::TAU,
        Curve2d::Ellipse(_) => std::f64::consts::TAU,
        // Only reached when IsPeriodic() is true; the rcad periodic carriers
        // are the circle and the ellipse.
        _ => std::f64::consts::TAU,
    }
}
