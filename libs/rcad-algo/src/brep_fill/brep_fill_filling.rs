//! OCCT BRepFill_Filling (TKBool/BRepFill) — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill_Filling.cxx
//! (L64-867) + BRepFill_Filling.hxx (L67-231).  The header-only data classes
//! BRepFill_EdgeFaceAndOrder / BRepFill_FaceAndOrder (fully defined in their
//! .hxx and consumed only here) are translated below as well.
//!
//! First consumer: BRepOffsetAPI_MakeFilling (Stage 2e carrier).
//!
//! Reported gaps (plan §0.6, annotated at the call sites) — see the split
//! unit [`super::brep_fill_filling_b`]:
//! - BRepFill_CurveConstraint (TKBool/BRepFill — separate unit);
//! - Adaptor3d_CurveOnSurface / GeomAdaptor_Surface / Geom2dAdaptor_Curve;
//! - GeomPlate_PlateG0Criterion / GeomPlate_MakeApprox / the UV ctor of
//!   GeomPlate_PointConstraint;
//! - GeomPlate curve-constraint consumption (Add / Curves2d()->Value).
//!
//! Architecture notes:
//! - `NCollection_Sequence` / `NCollection_List` map to `Vec` (OCCT index i
//!   maps to `[i - 1]`); `NCollection_DataMap(TopoDS_Shape, ...)` maps to
//!   `HashMap<ShapeKey, ...>` (TopTools_ShapeMapHasher == IsSame identity).
//! - `handle<GeomPlate_BuildPlateSurface> myBuilder` maps to
//!   `Option<BuildPlateSurface>` (OCCT resets it in Build; the accessors
//!   unwrap like the OCCT null-handle dereference would).
//! - `TopLoc_Location` handling is implicit: the rcad `TEdgeData`/`TFaceData`
//!   surfaces and pcurves are stored in world space.
//! - `BRepTools_WireExplorer` maps to the wire-ordered `TWireData::edges`
//!   list; `TopExp::Vertices(E, V1, V2)` (cumOri = true) maps to
//!   `generator::top_exp_vertices`.

use std::collections::HashMap;

use glam::DVec3;

use rcad_kernel::base::geom_api::project_on_surf::ProjectPointOnSurf;
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve2d, Surface3};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, GeomAbsShape, Orientation, Shape, TShape};

use crate::brep_fill::generator::{shape_key, shape_oriented, top_exp_vertices, ShapeKey};
use crate::brep_fill::offset_wire::append_edge_to_wire;
use crate::brep_fill::offset_wire_b::{
    brep_tool_curve_on_surface, brep_tool_surface, is_closed_wire, make_vertex,
    set_wire_closed, update_vertex_tolerance, vertex_point, wire_edges,
};
use crate::geomalgo::geomplate::build_plate_surface::BuildPlateSurface;
use crate::geomalgo::geomplate::point_constraint::PointConstraint;
use crate::geomalgo::geomplate::surface::GeomPlateSurface;

use super::brep_fill_filling_b::{
    builder_add_curve_constraint, curves_on_plate_value, point_constraint_uv, Adaptor3dCurve,
    Adaptor3dCurveOnSurface, BRepAdaptorCurve, BRepFillCurveConstraint, CurveConstraintHandle,
    Geom2dAdaptorCurve, GeomAdaptorSurface, GeomPlateMakeApprox, GeomPlatePlateG0Criterion,
};

/// OCCT GeomAbs_C0 (GeomAbs_Shape.hxx).
const GEOM_ABS_C0: GeomAbsShape = GeomAbsShape::C0;

/// OCCT M_PI.
const PI: f64 = std::f64::consts::PI;

// ===========================================================================
// BRepFill_EdgeFaceAndOrder (BRepFill_EdgeFaceAndOrder.hxx L20-40 — header
// only data class)
// ===========================================================================

/// OCCT BRepFill_EdgeFaceAndOrder (hxx L20-40): an edge, its support face and
/// the constraint order.
#[derive(Clone)]
pub struct BRepFillEdgeFaceAndOrder {
    /// OCCT: myEdge.
    pub my_edge: Shape,
    /// OCCT: myFace.
    pub my_face: Shape,
    /// OCCT: myOrder.
    pub my_order: GeomAbsShape,
}

impl BRepFillEdgeFaceAndOrder {
    /// OCCT BRepFill_EdgeFaceAndOrder(AnEdge, Face, Order).
    pub fn new(an_edge: &Shape, face: &Shape, order: GeomAbsShape) -> Self {
        BRepFillEdgeFaceAndOrder {
            my_edge: an_edge.clone(),
            my_face: face.clone(),
            my_order: order,
        }
    }
}

// ===========================================================================
// BRepFill_FaceAndOrder (BRepFill_FaceAndOrder.hxx L20-38 — header only data
// class)
// ===========================================================================

/// OCCT BRepFill_FaceAndOrder (hxx L20-38): a support face and the constraint
/// order.
#[derive(Clone)]
pub struct BRepFillFaceAndOrder {
    /// OCCT: myFace.
    pub my_face: Shape,
    /// OCCT: myOrder.
    pub my_order: GeomAbsShape,
}

impl BRepFillFaceAndOrder {
    /// OCCT BRepFill_FaceAndOrder(Face, Order).
    pub fn new(face: &Shape, order: GeomAbsShape) -> Self {
        BRepFillFaceAndOrder {
            my_face: face.clone(),
            my_order: order,
        }
    }
}

// ===========================================================================
// File statics (BRepFill_Filling.cxx L64-140)
// ===========================================================================

/// OCCT static gp_Vec MakeFinVec(aWire, aVertex) (L64-83).
fn make_fin_vec(brep: &BRep, a_wire: &Shape, a_vertex: &Shape) -> DVec3 {
    // OCCT L66: Vfirst, Vlast, Origin.
    let mut origin = Shape::null();
    // OCCT L67-81: the wire-ordered walk looking for the edge whose first or
    // last vertex is aVertex; Origin is the opposite vertex.
    for e in wire_edges(brep, a_wire) {
        let (vfirst, vlast) = top_exp_vertices(brep, &e);
        if vfirst.ptr_id() == a_vertex.ptr_id() {
            origin = vlast;
            break;
        }
        if vlast.ptr_id() == a_vertex.ptr_id() {
            origin = vfirst;
            break;
        }
    }
    // OCCT L82: gp_Vec(BRep_Tool::Pnt(Origin), BRep_Tool::Pnt(aVertex)).
    vertex_point(brep, a_vertex) - vertex_point(brep, &origin)
}

/// OCCT static TopoDS_Wire WireFromList(Edges) (L85-140).
fn wire_from_list(brep: &mut BRep, edges: &mut Vec<Shape>) -> Shape {
    // OCCT L87-89: BB; aWire; BB.MakeWire(aWire).
    let a_wire = brep.add_twire(vec![]);
    // OCCT L90-92.
    let mut an_edge = edges.first().cloned().expect("WireFromList: empty edge list");
    append_edge_to_wire(brep, &a_wire, &an_edge);
    edges.remove(0);

    // OCCT L94-95: V1, V2 with orientation.
    let (mut v1, mut v2) = top_exp_vertices(brep, &an_edge);

    // OCCT L97-136.
    while !edges.is_empty() {
        let mut found_index: Option<usize> = None;
        for (idx, itl_value) in edges.iter().enumerate() {
            an_edge = itl_value.clone();
            // OCCT L103-104: V3, V4 with orientation.
            let (v3, v4) = top_exp_vertices(brep, &an_edge);
            let same = |a: &Shape, b: &Shape| a.ptr_id() == b.ptr_id();
            if same(&v1, &v3) || same(&v1, &v4) || same(&v2, &v3) || same(&v2, &v4) {
                if same(&v1, &v3) {
                    an_edge.orientation = Orientation::Reversed;
                    v1 = v4;
                } else if same(&v1, &v4) {
                    v1 = v3;
                } else if same(&v2, &v3) {
                    v2 = v4;
                } else {
                    an_edge.orientation = Orientation::Reversed;
                    v2 = v3;
                }
                found_index = Some(idx);
                break;
            }
        }
        match found_index {
            None => {
                // OCCT L128-133.
                eprintln!(
                    "Warning: WireFromList: can't find the next edge. The wire is not \
                     complete, some edges are lost."
                );
                break;
            }
            Some(idx) => {
                // OCCT L134: BB.Add(aWire, anEdge).
                append_edge_to_wire(brep, &a_wire, &an_edge);
                // OCCT L135: Edges.Remove(itl).
                edges.remove(idx);
            }
        }
    }

    // OCCT L138: aWire.Closed(true).
    set_wire_closed(brep, &a_wire, true);
    // OCCT L139.
    a_wire
}

// ===========================================================================
// BRepFill_Filling (hxx L67-231 + cxx L144-867)
// ===========================================================================

/// OCCT BRepFill_Filling — the N-side filling engine.
pub struct BRepFillFilling {
    /// OCCT: myBuilder (hxx L210) — reset in Build.
    my_builder: Option<BuildPlateSurface>,
    /// OCCT: myBoundary (hxx L211).
    my_boundary: Vec<BRepFillEdgeFaceAndOrder>,
    /// OCCT: myConstraints (hxx L212).
    my_constraints: Vec<BRepFillEdgeFaceAndOrder>,
    /// OCCT: myFreeConstraints (hxx L213).
    my_free_constraints: Vec<BRepFillFaceAndOrder>,
    /// OCCT: myPoints (hxx L214).
    my_points: Vec<PointConstraint>,
    /// OCCT: myOldNewMap (hxx L215).
    my_old_new_map: HashMap<ShapeKey, Shape>,
    /// OCCT: myGenerated (hxx L216).
    my_generated: Vec<Shape>,
    /// OCCT: myFace (hxx L217).
    my_face: Shape,
    /// OCCT: myInitFace (hxx L218).
    my_init_face: Shape,
    /// OCCT: myTol2d (hxx L219).
    my_tol2d: f64,
    /// OCCT: myTol3d (hxx L220).
    my_tol3d: f64,
    /// OCCT: myTolAng (hxx L221).
    my_tolang: f64,
    /// OCCT: myTolCurv (hxx L222).
    my_tolcurv: f64,
    /// OCCT: myMaxDeg (hxx L223).
    my_max_deg: i32,
    /// OCCT: myMaxSegments (hxx L224).
    my_max_segments: i32,
    /// OCCT: myDegree (hxx L225).
    my_degree: i32,
    /// OCCT: myNbPtsOnCur (hxx L226).
    my_nb_pts_on_cur: i32,
    /// OCCT: myNbIter (hxx L227).
    my_nb_iter: i32,
    /// OCCT: myAnisotropie (hxx L228).
    my_anisotropie: bool,
    /// OCCT: myIsInitFaceGiven (hxx L229).
    my_is_init_face_given: bool,
    /// OCCT: myIsDone (hxx L230).
    my_is_done: bool,
}

impl Default for BRepFillFilling {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFillFilling {
    /// OCCT BRepFill_Filling::BRepFill_Filling(Degree = 3, NbPtsOnCur = 15,
    /// NbIter = 2, Anisotropie = false, Tol2d = 0.00001, Tol3d = 0.0001,
    /// TolAng = 0.01, TolCurv = 0.1, MaxDeg = 8, MaxSegments = 9)
    /// (cxx L144-171) — the OCCT default values.
    pub fn new() -> Self {
        Self::new_with_args(3, 15, 2, false, 0.00001, 0.0001, 0.01, 0.1, 8, 9)
    }

    /// OCCT constructor with explicit parameters (cxx L144-171).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_args(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
        tol2d: f64,
        tol3d: f64,
        tolang: f64,
        tolcurv: f64,
        max_deg: i32,
        max_segments: i32,
    ) -> Self {
        BRepFillFilling {
            my_builder: None,
            my_boundary: Vec::new(),
            my_constraints: Vec::new(),
            my_free_constraints: Vec::new(),
            my_points: Vec::new(),
            my_old_new_map: HashMap::new(),
            my_generated: Vec::new(),
            my_face: Shape::null(),
            my_init_face: Shape::null(),
            // OCCT L155-170.
            my_tol2d: tol2d,
            my_tol3d: tol3d,
            my_tolang: tolang,
            my_tolcurv: tolcurv,
            my_max_deg: max_deg,
            my_max_segments: max_segments,
            my_degree: degree,
            my_nb_pts_on_cur: nb_pts_on_cur,
            my_nb_iter: nb_iter,
            my_anisotropie: anisotropie,
            my_is_init_face_given: false,
            my_is_done: false,
        }
    }

    /// OCCT BRepFill_Filling::SetConstrParam(Tol2d, Tol3d, TolAng, TolCurv)
    /// (cxx L175-184).
    pub fn set_constr_param(&mut self, tol2d: f64, tol3d: f64, tolang: f64, tolcurv: f64) {
        // OCCT L180-183.
        self.my_tol2d = tol2d;
        self.my_tol3d = tol3d;
        self.my_tolang = tolang;
        self.my_tolcurv = tolcurv;
    }

    /// OCCT BRepFill_Filling::SetResolParam(Degree, NbPtsOnCur, NbIter,
    /// Anisotropie) (cxx L188-197).
    pub fn set_resol_param(
        &mut self,
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
    ) {
        // OCCT L193-196.
        self.my_degree = degree;
        self.my_nb_pts_on_cur = nb_pts_on_cur;
        self.my_nb_iter = nb_iter;
        self.my_anisotropie = anisotropie;
    }

    /// OCCT BRepFill_Filling::SetApproxParam(MaxDeg, MaxSegments)
    /// (cxx L201-205).
    pub fn set_approx_param(&mut self, max_deg: i32, max_segments: i32) {
        // OCCT L203-204.
        self.my_max_deg = max_deg;
        self.my_max_segments = max_segments;
    }

    /// OCCT BRepFill_Filling::LoadInitSurface(aFace) (cxx L209-213).
    pub fn load_init_surface(&mut self, a_face: &Shape) {
        // OCCT L211-212.
        self.my_init_face = a_face.clone();
        self.my_is_init_face_given = true;
    }

    /// OCCT BRepFill_Filling::Add(anEdge, Order, IsBound = true)
    /// (cxx L217-231).
    pub fn add(&mut self, an_edge: &Shape, order: GeomAbsShape, is_bound: bool) -> i32 {
        // OCCT L219-220.
        let null_face = Shape::null();
        let edge_face_and_order = BRepFillEdgeFaceAndOrder::new(an_edge, &null_face, order);
        // OCCT L221-230.
        if is_bound {
            self.my_boundary.push(edge_face_and_order);
            self.my_boundary.len() as i32
        } else {
            self.my_constraints.push(edge_face_and_order);
            (self.my_boundary.len() + self.my_free_constraints.len() + self.my_constraints.len())
                as i32
        }
    }

    /// OCCT BRepFill_Filling::Add(anEdge, Support, Order, IsBound = true)
    /// (cxx L237-253).
    pub fn add_with_support(
        &mut self,
        an_edge: &Shape,
        support: &Shape,
        order: GeomAbsShape,
        is_bound: bool,
    ) -> i32 {
        // OCCT L242-252.
        let edge_face_and_order = BRepFillEdgeFaceAndOrder::new(an_edge, support, order);
        if is_bound {
            self.my_boundary.push(edge_face_and_order);
            self.my_boundary.len() as i32
        } else {
            self.my_constraints.push(edge_face_and_order);
            (self.my_boundary.len() + self.my_free_constraints.len() + self.my_constraints.len())
                as i32
        }
    }

    /// OCCT BRepFill_Filling::Add(Support, Order) — free constraint
    /// (cxx L259-264).
    pub fn add_free_constraint(&mut self, support: &Shape, order: GeomAbsShape) -> i32 {
        // OCCT L261-263.
        let face_and_order = BRepFillFaceAndOrder::new(support, order);
        self.my_free_constraints.push(face_and_order);
        (self.my_boundary.len() + self.my_free_constraints.len()) as i32
    }

    /// OCCT BRepFill_Filling::Add(Point) — punctual constraint (cxx L268-275).
    pub fn add_point(&mut self, point: DVec3) -> i32 {
        // OCCT L270-271.
        let a_pc = PointConstraint::new(point, GEOM_ABS_C0 as i32, self.my_tol3d);
        // OCCT L272-274.
        self.my_points.push(a_pc);
        (self.my_boundary.len()
            + self.my_free_constraints.len()
            + self.my_constraints.len()
            + self.my_points.len()) as i32
    }

    /// OCCT BRepFill_Filling::Add(U, V, Support, Order) — point constraint on
    /// a face (cxx L281-299).
    pub fn add_point_on_support(
        &mut self,
        brep: &BRep,
        u: f64,
        v: f64,
        support: &Shape,
        order: GeomAbsShape,
    ) -> i32 {
        // OCCT L286-287: HSurf = new BRepAdaptor_Surface(); Initialize(Support)
        // — mapped to the TFaceData surface access.
        let cur_surface = brep_tool_surface(brep, support)
            .expect("BRepFill_Filling::Add: the support face has no surface");
        // OCCT L288-295.
        let a_pc = point_constraint_uv(
            u,
            v,
            cur_surface,
            order,
            self.my_tol3d,
            self.my_tolang,
            self.my_tolcurv,
        );
        // OCCT L296-298.
        self.my_points.push(a_pc);
        (self.my_boundary.len()
            + self.my_free_constraints.len()
            + self.my_constraints.len()
            + self.my_points.len()) as i32
    }

    /// OCCT BRepFill_Filling::AddConstraints(SeqOfConstraints) (cxx L303-389).
    #[allow(clippy::too_many_arguments)]
    fn add_constraints(
        brep: &BRep,
        builder: &mut BuildPlateSurface,
        seq_of_constraints: &[BRepFillEdgeFaceAndOrder],
        is_init_face_given: bool,
        init_face: &Shape,
        nb_pts_on_cur: i32,
        tol3d: f64,
        tolang: f64,
        tolcurv: f64,
    ) {
        // OCCT L306-311: CurEdge, CurFace, CurOrder, Constr, i.
        for i in 1..=seq_of_constraints.len() {
            let cur_edge = seq_of_constraints[i - 1].my_edge.clone();
            let cur_face = seq_of_constraints[i - 1].my_face.clone();
            let cur_order = seq_of_constraints[i - 1].my_order;

            // OCCT L318-375.
            let mut constr: CurveConstraintHandle = if cur_face.is_null() {
                if cur_order == GEOM_ABS_C0 {
                    // OCCT L322-325.
                    let hcurve = BRepAdaptorCurve::initialize(brep, &cur_edge);
                    let a_hcurve = Adaptor3dCurve::BRep(hcurve); // to avoid ambiguity
                    CurveConstraintHandle::BRepFill(BRepFillCurveConstraint::new(
                        &a_hcurve,
                        cur_order,
                        nb_pts_on_cur,
                        tol3d,
                    ))
                } else {
                    // Pas de representation Topologique
                    // On prend une representation Geometrique : au pif !
                    // OCCT L330-334: the first CurveOnSurface of the edge —
                    // mapped to the first pcurve entry of TEdgeData.
                    let ed = brep.edge(cur_edge.clone());
                    match ed.pcurves.values().next().cloned() {
                        None => {
                            // OCCT L335-339: Surface.IsNull() -> throw.
                            panic!("Add")
                        }
                        Some((c2d, f, l)) => {
                            // OCCT L340-341: Surface->Copy() + Transform(loc)
                            // — implicit: rcad pcurves/surfaces are stored in
                            // world space (architecture note).
                            // OCCT L342-347.
                            let surf_adaptor = GeomAdaptorSurface::new(
                                brep_tool_surface_of_edge(brep, &cur_edge)
                                    .expect("Add: no surface under the edge"),
                            );
                            let curve2d_adaptor = Geom2dAdaptorCurve::new(c2d, f, l);
                            let curv_on_surf =
                                Adaptor3dCurveOnSurface::new(curve2d_adaptor, surf_adaptor);
                            // OCCT L349-354.
                            CurveConstraintHandle::GeomPlate {
                                curv_on_surf,
                                order: cur_order,
                                nb_pts: nb_pts_on_cur,
                                tol3d,
                                tolang,
                                tolcurv,
                            }
                        }
                    }
                }
            } else {
                // OCCT L359-374.
                let surf = brep_tool_surface(brep, &cur_face)
                    .expect("AddConstraints: the face has no surface");
                // OCCT L361-362: BRepAdaptor_Curve2d Initialize(CurEdge,
                // CurFace) — the pcurve of the edge on the face.
                let (c2d, f, l) = brep_tool_curve_on_surface(brep, &cur_edge, &cur_face).expect(
                    // OCCT L363-365 comment.
                    "If CurEdge has no 2d representation on CurFace, there will be \
                     exception \"Attempt to access to null object\" in this \
                     initialization (null pcurve).",
                );
                let surf_adaptor = GeomAdaptorSurface::new(surf);
                let curve2d_adaptor = Geom2dAdaptorCurve::new(c2d, f, l);
                let curv_on_surf = Adaptor3dCurveOnSurface::new(curve2d_adaptor, surf_adaptor);
                // OCCT L369-374.
                CurveConstraintHandle::BRepFill(BRepFillCurveConstraint::new_on_surface(
                    &curv_on_surf,
                    cur_order,
                    nb_pts_on_cur,
                    tol3d,
                    tolang,
                    tolcurv,
                ))
            };

            // OCCT L376-386.
            if is_init_face_given {
                // OCCT L378-380.
                if let Some((curve2d, first_par, last_par)) =
                    brep_tool_curve_on_surface(brep, &cur_edge, init_face)
                {
                    // OCCT L383: Curve2d = new Geom2d_TrimmedCurve(...).
                    let trimmed = trimmed_curve2d(curve2d, first_par, last_par);
                    // OCCT L384.
                    constr.set_curve2d_on_surf(trimmed);
                }
            }
            // OCCT L387.
            builder_add_curve_constraint(builder, constr);
        }
    }

    /// OCCT BRepFill_Filling::BuildWires(EdgeList, WireList) (cxx L393-482).
    fn build_wires(
        &mut self,
        brep: &mut BRep,
        edge_list: &mut Vec<Shape>,
        wire_list: &mut Vec<Shape>,
    ) {
        // OCCT L399-481.
        while !edge_list.is_empty() {
            // OCCT L401-404: BRepLib_MakeWire MW; FirstEdge; MW.Add; RemoveFirst.
            let mw = brep.add_twire(vec![]);
            let first_edge = edge_list.first().cloned().expect("BuildWires: empty list");
            append_edge_to_wire(brep, &mw, &first_edge);
            edge_list.remove(0);
            // OCCT L405: V_wire[2], V_edge[2].
            let mut v_wire: [Shape; 2] = [Shape::null(), Shape::null()];
            let mut v_edge: [Shape; 2] = [Shape::null(), Shape::null()];

            // OCCT L407: for (;;).
            loop {
                // OCCT L409-410.
                let (w0, w1) = top_exp_vertices(brep, &mw);
                v_wire[0] = w0;
                v_wire[1] = w1;
                let mut found = false;
                for idx in 0..edge_list.len() {
                    let cur_edge = edge_list[idx].clone();
                    let (e0, e1) = top_exp_vertices(brep, &cur_edge);
                    v_edge[0] = e0;
                    v_edge[1] = e1;
                    for i in 0..2 {
                        for j in 0..2 {
                            if v_wire[i].ptr_id() == v_edge[j].ptr_id() {
                                // OCCT L422-424.
                                append_edge_to_wire(brep, &mw, &cur_edge);
                                edge_list.remove(idx);
                                found = true;
                                break;
                            }
                        }
                        if found {
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
                if !found {
                    // try to find geometric coincidence
                    // OCCT L440-442.
                    let p_wire = [
                        vertex_point(brep, &v_wire[0]),
                        vertex_point(brep, &v_wire[1]),
                    ];
                    for idx in 0..edge_list.len() {
                        let cur_edge = edge_list[idx].clone();
                        let (e0, e1) = top_exp_vertices(brep, &cur_edge);
                        v_edge[0] = e0;
                        v_edge[1] = e1;
                        for i in 0..2 {
                            for j in 0..2 {
                                // OCCT L451-453.
                                let a_dist = p_wire[i].distance(vertex_point(brep, &v_edge[j]));
                                if a_dist < shape_tolerance(brep, &v_wire[i])
                                    && a_dist < shape_tolerance(brep, &v_edge[j])
                                {
                                    // OCCT L455-461.
                                    append_edge_to_wire(brep, &mw, &cur_edge);
                                    let new_edge = wire_single_edge(brep, &mw);
                                    self.my_old_new_map.insert(
                                        shape_key(&shape_oriented(&cur_edge, Orientation::Forward)),
                                        shape_oriented(&new_edge, Orientation::Forward),
                                    );
                                    edge_list.remove(idx);
                                    found = true;
                                    break;
                                }
                            }
                            if found {
                                break;
                            }
                        }
                        if found {
                            break;
                        }
                    }
                }
                if !found {
                    // end of current wire, begin next wire
                    // OCCT L477.
                    wire_list.push(mw.clone());
                    break;
                }
            } // end of for (;;)
        } // end of while (! EdgeList.IsEmpty())
    }

    /// OCCT BRepFill_Filling::FindExtremitiesOfHoles(WireList, VerSeq) const
    /// (cxx L486-567).
    fn find_extremities_of_holes(
        &self,
        brep: &BRep,
        wire_list: &[Shape],
        ver_seq: &mut Vec<Shape>,
    ) {
        // OCCT L489-494: WireSeq = copy of WireList.
        let mut wire_seq: Vec<Shape> = wire_list.to_vec();

        // OCCT L496-498.
        let the_wire = wire_seq
            .first()
            .cloned()
            .expect("FindExtremitiesOfHoles: empty wire list");
        wire_seq.remove(0);

        // OCCT L500-503.
        if is_closed_wire(brep, &the_wire) {
            return;
        }

        // OCCT L505-506.
        let (vfirst, vlast) = top_exp_vertices(brep, &the_wire);

        // OCCT L508-514: The Wire is closed indeed despite its
        // being not detected earlier.
        if vfirst.ptr_id() == vlast.ptr_id() {
            return;
        }

        // OCCT L516-518.
        let mut fin_vec = make_fin_vec(brep, &the_wire, &vlast);
        let mut the_vertex = vlast;
        ver_seq.push(the_vertex.clone());

        // OCCT L520-564.
        while !wire_seq.is_empty() {
            let mut min_vtx = Shape::null();
            let mut min_ind = 1usize;
            let mut is_last = false;
            let mut min_angle = PI;

            for i in 1..=wire_seq.len() {
                let cur_wire = wire_seq[i - 1].clone();
                let (vfirst_i, vlast_i) = top_exp_vertices(brep, &cur_wire);

                // OCCT L532-539.
                let mut angle = angle_between(
                    fin_vec,
                    vertex_point(brep, &the_vertex) - vertex_point(brep, &vfirst_i),
                );
                if angle < min_angle {
                    min_angle = angle;
                    min_vtx = vfirst_i;
                    min_ind = i;
                    is_last = true;
                }
                // OCCT L540-547.
                angle = angle_between(
                    fin_vec,
                    vertex_point(brep, &the_vertex) - vertex_point(brep, &vlast_i),
                );
                if angle < min_angle {
                    min_angle = angle;
                    min_vtx = vlast_i;
                    min_ind = i;
                    is_last = false;
                }
            }
            // OCCT L549-551.
            ver_seq.push(min_vtx);
            let min_wire = wire_seq[min_ind - 1].clone();
            let (vfirst_m, vlast_m) = top_exp_vertices(brep, &min_wire);
            // OCCT L552-561.
            if is_last {
                fin_vec = make_fin_vec(brep, &min_wire, &vlast_m);
                the_vertex = vlast_m;
            } else {
                fin_vec = make_fin_vec(brep, &min_wire, &vfirst_m);
                the_vertex = vfirst_m;
            }
            // OCCT L562-563.
            ver_seq.push(the_vertex.clone());
            wire_seq.remove(min_ind - 1);
        }
        // OCCT L565-566.
        let (vfirst_f, _vlast_f) = top_exp_vertices(brep, &the_wire);
        ver_seq.push(vfirst_f);
    }

    /// OCCT BRepFill_Filling::Build() (cxx L571-778).
    pub fn build(&mut self, brep: &mut BRep) {
        // OCCT L573-580: myBuilder.reset(new GeomPlate_BuildPlateSurface(...)).
        self.my_builder = Some(BuildPlateSurface::new(
            self.my_degree,
            self.my_nb_pts_on_cur,
            self.my_nb_iter,
            self.my_tol2d,
            self.my_tol3d,
            self.my_tolang,
            self.my_tolcurv,
            self.my_anisotropie,
        ));

        // OCCT L582-586.
        if self.my_boundary.is_empty() {
            self.my_is_done = false;
            return;
        }

        // OCCT L592-593: Creating array of vertices: extremities of wires.
        let mut ver_seq: Vec<Shape> = Vec::new();

        // OCCT L595-600: Building missing bounds.
        let mut edge_list: Vec<Shape> = Vec::new();
        let mut wire_list: Vec<Shape> = Vec::new();
        for i in 1..=self.my_boundary.len() {
            edge_list.push(self.my_boundary[i - 1].my_edge.clone());
        }

        self.build_wires(brep, &mut edge_list, &mut wire_list);
        self.find_extremities_of_holes(brep, &wire_list, &mut ver_seq);

        // OCCT L605-660: Searching for surfaces for missing bounds.
        for j in 1..=self.my_free_constraints.len() {
            // OCCT L608-614.
            let cur_face = self.my_free_constraints[j - 1].my_face.clone();
            let cur_surface = match brep_tool_surface(brep, &cur_face) {
                Some(s) => s,
                None => continue,
            };

            // OCCT L617-659: for (i = 1; i <= VerSeq.Length(); i += 2).
            let mut i = 1usize;
            while i + 1 <= ver_seq.len() {
                let first_vtx = ver_seq[i - 1].clone();
                let last_vtx = ver_seq[i].clone();

                // OCCT L622-627.
                let first_pnt = vertex_point(brep, &first_vtx);
                let mut projector = ProjectPointOnSurf::new_point(first_pnt, &cur_surface);
                if projector.lower_distance() > CONFUSION {
                    i += 2;
                    continue;
                }
                // OCCT L628.
                let (u1, v1) = projector.lower_distance_parameters();

                // OCCT L630-634 (commented-out classifier) — skipped as in OCCT.

                // OCCT L636-642.
                let last_pnt = vertex_point(brep, &last_vtx);
                let mut projector = ProjectPointOnSurf::new_point(last_pnt, &cur_surface);
                if projector.lower_distance() > CONFUSION {
                    i += 2;
                    continue;
                }
                // OCCT L642.
                let (u2, v2) = projector.lower_distance_parameters();

                // OCCT L644-648 (commented-out classifier) — skipped as in OCCT.

                // Making the constraint
                // OCCT L651-654: a 2d Bezier line through (U1,V1)-(U2,V2).
                let line2d = bezier2d_from_points([glam::DVec2::new(u1, v1), glam::DVec2::new(u2, v2)]);
                // OCCT L655: BRepLib_MakeEdge(Line2d, CurSurface, FirstVtx,
                // LastVtx).
                let e = make_edge_pcurve(brep, line2d, &first_vtx, &last_vtx);
                // OCCT L656.
                let order = self.my_free_constraints[j - 1].my_order;
                self.add_with_support(&e, &cur_face, order, true);
                // OCCT L657: VerSeq.Remove(i, i + 1).
                ver_seq.drain((i - 1)..(i + 1));
                break;
            } // for (i = 1; i <= VerSeq.Length(); i += 2)
        } // for (j = 1; j <= myFreeConstraints.Length(); j++)

        // OCCT L662-668: Load initial surface to myBuilder if it is given.
        if self.my_is_init_face_given {
            if let Some(surf_init) = brep_tool_surface(brep, &self.my_init_face) {
                self.my_builder
                    .as_mut()
                    .expect("myBuilder")
                    .load_init_surface(surf_init);
            }
        }

        // OCCT L670-677: Adding constraints to myBuilder.
        {
            let builder = self.my_builder.as_mut().expect("myBuilder");
            let boundary = self.my_boundary.clone();
            Self::add_constraints(
                brep,
                builder,
                &boundary,
                self.my_is_init_face_given,
                &self.my_init_face,
                self.my_nb_pts_on_cur,
                self.my_tol3d,
                self.my_tolang,
                self.my_tolcurv,
            );
            // OCCT L672.
            builder.set_nb_bounds(self.my_boundary.len() as i32);
            let constraints = self.my_constraints.clone();
            Self::add_constraints(
                brep,
                builder,
                &constraints,
                self.my_is_init_face_given,
                &self.my_init_face,
                self.my_nb_pts_on_cur,
                self.my_tol3d,
                self.my_tolang,
                self.my_tolcurv,
            );
            for i in 1..=self.my_points.len() {
                // OCCT L676: myBuilder->Add(myPoints(i)).
                let pc = self.my_points[i - 1].clone();
                builder.add_point_constraint(pc);
            }
        }

        // OCCT L679-688.
        self.my_builder.as_mut().expect("myBuilder").perform();
        if self.my_builder.as_ref().expect("myBuilder").is_done() {
            self.my_is_done = true;
        } else {
            self.my_is_done = false;
            return;
        }

        // OCCT L690-693.
        let gplate_surface: Option<&GeomPlateSurface> =
            self.my_builder.as_ref().expect("myBuilder").surface();
        let gplate = match gplate_surface {
            Some(gp) => gp,
            None => {
                self.my_is_done = false;
                return;
            }
        };
        let gplate_basis: Surface3 = gplate.basis_surface().clone();
        // Approximation
        // OCCT L693: double dmax = 1.1 * myBuilder->G0Error(); //???????????
        let dmax = 1.1 * self.my_builder.as_ref().expect("myBuilder").g0_error();
        let surface: Surface3;
        if !self.my_is_init_face_given {
            // OCCT L697-706.
            let seuil;
            let mut s2d: Vec<glam::DVec2> = Vec::new();
            let mut s3d: Vec<glam::DVec3> = Vec::new();
            // OCCT L701-702: Disc2dContour(4, S2D); Disc3dContour(4, 1, S3D).
            {
                let builder = self.my_builder.as_mut().expect("myBuilder");
                builder.disc2d_contour(&mut s2d);
                builder.disc3d_contour(1, &mut s3d);
            }
            // OCCT L703.
            seuil = self
                .my_tol3d
                .max(10.0 * self.my_builder.as_ref().expect("myBuilder").g0_error());
            let _ = &mut s2d;
            let _ = &mut s3d;
            // OCCT L704-706.
            let criterion = GeomPlatePlateG0Criterion::new(&s2d, &s3d, seuil);
            let mut approx = GeomPlateMakeApprox::new_criterion(
                &gplate_basis,
                &criterion,
                self.my_tol3d,
                self.my_max_segments,
                self.my_max_deg,
            );
            surface = approx.surface();
        } else {
            // OCCT L710-715.
            let mut approx = GeomPlateMakeApprox::new_dmax(
                &gplate_basis,
                self.my_tol3d,
                self.my_max_segments,
                self.my_max_deg,
                dmax,
                0,
            );
            surface = approx.surface();
        }

        // Build the final wire and final face
        // OCCT L719-721.
        let mut final_edges: Vec<Shape> = Vec::new();
        // OCCT L720: CurvesOnPlate = myBuilder->Curves2d() — the rcad GeomPlate
        // port keeps the point-constraint path only; the call is kept for form
        // and the per-edge access goes through the GAP helper.
        self.my_builder.as_ref().expect("myBuilder").curves2d();
        // OCCT L721.
        let mut bb = BRepBuilder::new();
        for i in 1..=self.my_boundary.len() {
            // OCCT L724-730.
            let init_edge = self.my_boundary[i - 1].my_edge.clone();
            let mut an_edge = init_edge.clone();
            an_edge.orientation = Orientation::Forward;
            if let Some(mapped) = self.my_old_new_map.get(&shape_key(&an_edge)) {
                an_edge = mapped.clone();
            }

            // OCCT L732.
            let a_curve_on_plate = curves_on_plate_value(i);

            // OCCT L734.
            let new_edge = brep.empty_copied(&an_edge);

            // OCCT L736-737.
            let (v1, v2) = top_exp_vertices(brep, &an_edge);

            // OCCT L739-749.
            let new_v1;
            if let Some(mapped) = self.my_old_new_map.get(&shape_key(&v1)) {
                new_v1 = mapped.clone();
            } else {
                let a_pnt = vertex_point(brep, &v1);
                let nv1 = make_vertex(brep, a_pnt);
                update_vertex_tolerance(brep, &nv1, dmax);
                self.my_old_new_map.insert(
                    shape_key(&shape_oriented(&v1, Orientation::Forward)),
                    nv1.clone(),
                );
                new_v1 = nv1;
            }

            // OCCT L751-761.
            let new_v2;
            if let Some(mapped) = self.my_old_new_map.get(&shape_key(&v2)) {
                new_v2 = mapped.clone();
            } else {
                let a_pnt = vertex_point(brep, &v2);
                let nv2 = make_vertex(brep, a_pnt);
                update_vertex_tolerance(brep, &nv2, dmax);
                self.my_old_new_map.insert(
                    shape_key(&shape_oriented(&v2, Orientation::Forward)),
                    nv2.clone(),
                );
                new_v2 = nv2;
            }

            // OCCT L763-766.
            let mut new_v1_fwd = new_v1.clone();
            new_v1_fwd.orientation = Orientation::Forward;
            bb.add_to_edge(brep, new_edge.clone(), new_v1_fwd);
            let mut new_v2_rev = new_v2.clone();
            new_v2_rev.orientation = Orientation::Reversed;
            bb.add_to_edge(brep, new_edge.clone(), new_v2_rev);
            // OCCT L767-768: UpdateEdge(NewEdge, aCurveOnPlate, Surface, Loc,
            // dmax) — Loc is identity in the world-space mapping; the pcurve
            // is registered under the standalone-surface key (0, 0) set at
            // make_edge_pcurve time.
            update_edge_pcurve_standalone(brep, &new_edge, a_curve_on_plate, dmax);
            // OCCT L769 (commented-out SameRange) — skipped as in OCCT.
            // OCCT L770: BRepLib::SameParameter(NewEdge, dmax, true).
            crate::topalgo::brep_lib::brep_lib::BRepLib::same_parameter(&new_edge, dmax);
            // OCCT L771.
            final_edges.push(new_edge.clone());
            // OCCT L772.
            self.my_old_new_map.insert(
                shape_key(&shape_oriented(&init_edge, Orientation::Forward)),
                shape_oriented(&new_edge, Orientation::Forward),
            );
        }

        // OCCT L775.
        let final_wire = wire_from_list(brep, &mut final_edges);

        // OCCT L777: myFace = BRepLib_MakeFace(Surface, FinalWire).
        self.my_face = bb.make_face(brep, Some(surface), final_wire);
    }

    /// OCCT BRepFill_Filling::IsDone() const (cxx L782-785).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepFill_Filling::Face() const (cxx L789-792).
    pub fn face(&self) -> Shape {
        self.my_face.clone()
    }

    /// OCCT BRepFill_Filling::Generated(S) (cxx L798-808).
    pub fn generated(&mut self, s: &Shape) -> &Vec<Shape> {
        self.my_generated.clear();

        if let Some(mapped) = self.my_old_new_map.get(&shape_key(s)) {
            self.my_generated.push(mapped.clone());
        }

        &self.my_generated
    }

    /// OCCT BRepFill_Filling::G0Error() const (cxx L814-817).
    pub fn g0_error(&self) -> f64 {
        self.my_builder.as_ref().expect("myBuilder").g0_error()
    }

    /// OCCT BRepFill_Filling::G1Error() const (cxx L824-827).
    pub fn g1_error(&self) -> f64 {
        self.my_builder.as_ref().expect("myBuilder").g1_error()
    }

    /// OCCT BRepFill_Filling::G2Error() const (cxx L834-837).
    pub fn g2_error(&self) -> f64 {
        self.my_builder.as_ref().expect("myBuilder").g2_error()
    }

    /// OCCT BRepFill_Filling::G0Error(Index) (cxx L844-847).
    pub fn g0_error_index(&mut self, index: i32) -> f64 {
        self.my_builder
            .as_mut()
            .expect("myBuilder")
            .g0_error_index(index)
    }

    /// OCCT BRepFill_Filling::G1Error(Index) (cxx L854-857).
    pub fn g1_error_index(&mut self, index: i32) -> f64 {
        self.my_builder
            .as_mut()
            .expect("myBuilder")
            .g1_error_index(index)
    }

    /// OCCT BRepFill_Filling::G2Error(Index) (cxx L864-867).
    pub fn g2_error_index(&mut self, index: i32) -> f64 {
        self.my_builder
            .as_mut()
            .expect("myBuilder")
            .g2_error_index(index)
    }
}

// ===========================================================================
// Kernel-mapping helpers
// ===========================================================================

/// OCCT BRep_Tool::Tolerance (vertex or edge).
fn shape_tolerance(brep: &BRep, s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        _ => 0.0,
    }
}

/// OCCT gp_Vec::Angle (the normalized-vector angle).
fn angle_between(a: DVec3, b: DVec3) -> f64 {
    let la = a.length();
    let lb = b.length();
    if la <= f64::EPSILON || lb <= f64::EPSILON {
        return 0.0;
    }
    let c = (a.dot(b) / (la * lb)).clamp(-1.0, 1.0);
    c.acos()
}

/// OCCT BRepLib_MakeWire::Edge() — the single edge of the built wire.
fn wire_single_edge(brep: &BRep, w: &Shape) -> Shape {
    let edges = wire_edges(brep, w);
    edges
        .first()
        .cloned()
        .expect("MakeWire::Edge: empty wire")
}

/// OCCT Geom2d_TrimmedCurve(C2d, FirstPar, LastPar).
fn trimmed_curve2d(c: Curve2d, first: f64, last: f64) -> Curve2d {
    rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
        curve: Box::new(c),
        t_min: first,
        t_max: last,
    })
}

/// OCCT Geom2d_BezierCurve(Points) through two 2d points (parameter range
/// [0, 1], as in OCCT).
fn bezier2d_from_points(pts: [glam::DVec2; 2]) -> Curve2d {
    rcad_kernel::geom::Curve2d::Bezier(rcad_kernel::geom::BezierCurve2 {
        control_points: pts.to_vec(),
        weights: vec![1.0; pts.len()],
    })
}

/// OCCT BRepLib_MakeEdge(C2d, Surface, V1, V2) — a pcurve-only edge over a
/// standalone surface (no support face yet).  The rcad port stores the pcurve
/// under the standalone-surface key (0, 0); UpdateEdge in Build reuses the
/// same key (architecture note: OCCT keys pcurves by (face, location)).
fn make_edge_pcurve(brep: &mut BRep, c2d: Curve2d, v1: &Shape, v2: &Shape) -> Shape {
    let mut v1_fwd = v1.clone();
    v1_fwd.orientation = Orientation::Forward;
    let mut v2_rev = v2.clone();
    v2_rev.orientation = Orientation::Reversed;
    // Geom2d_BezierCurve: FirstParameter = 0, LastParameter = 1.
    let e = brep.add_tedge(None, v1_fwd, v2_rev, [0.0, 1.0]);
    let ed = brep.edge_mut_inplace(e.clone());
    ed.pcurves.insert((0u64, 0u32), (c2d, 0.0, 1.0));
    e
}

/// OCCT BRep_Tool::Surface of the first pcurve's face (the L334 out-param
/// form) — recovered from the face of the edge's first pcurve representation.
fn brep_tool_surface_of_edge(brep: &BRep, e: &Shape) -> Option<Surface3> {
    let ed = brep.edge(e.clone());
    let (face_key, _loc) = ed.pcurves.keys().next().copied()?;
    for idx in 0..brep.tshapes.len() {
        let s = brep.shape_at(idx);
        if s.ptr_id() == face_key {
            return brep_tool_surface(brep, &s);
        }
    }
    None
}

/// OCCT BRep_Builder::UpdateEdge(NewEdge, PCurve, Surface, Loc, Tol) for a
/// standalone surface (no face handle in the rcad data model — the pcurve is
/// registered under the key (0, 0), see make_edge_pcurve).
fn update_edge_pcurve_standalone(brep: &mut BRep, e: &Shape, pcurve: Curve2d, tol: f64) {
    let ed = brep.edge_mut_inplace(e.clone());
    ed.pcurves.insert((0u64, 0u32), (pcurve, 0.0, 1.0));
    ed.tolerance = ed.tolerance.max(tol);
}

// ===========================================================================
// Unit tests (write-only per the batch discipline: `cargo test` is not run;
// compile-checked via `cargo check --tests`)
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT default constructor values (hxx L73-82 defaults).
    #[test]
    fn filling_defaults_and_add_indices() {
        let mut f = BRepFillFilling::new();
        assert!(!f.is_done());
        // OCCT Add(E, GeomAbs_C0, true) -> myBoundary.Length() = 1.
        let e = Shape::null();
        assert_eq!(f.add(&e, GeomAbsShape::C0, true), 1);
        // OCCT Add(E, NullFace, GeomAbs_G1, false) ->
        // Boundary + FreeConstraints + Constraints = 2.
        assert_eq!(f.add(&e, GeomAbsShape::G1, false), 2);
        // OCCT Add(Point) -> ... + myPoints.Length() = 3.
        assert_eq!(f.add_point(DVec3::ZERO), 3);
        // OCCT Add(Support, Order) -> Boundary + FreeConstraints = 2.
        let face = Shape::null();
        assert_eq!(f.add_free_constraint(&face, GeomAbsShape::C0), 2);
        // Generated of an unbound shape is empty.
        assert!(f.generated(&e).is_empty());
    }
}
