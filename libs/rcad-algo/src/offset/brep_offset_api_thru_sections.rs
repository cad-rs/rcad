//! OCCT BRepOffsetAPI_ThruSections — 1:1 translation (class + Build +
//! CreateRuled).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_ThruSections.cxx L60-702 (statics, class,
//!         Build, CreateRuled) + BRepOffsetAPI_ThruSections.hxx (L52-215).
//! CreateSmoothed / TotalSurf / EdgeToBSpline / Generated (cxx L706-1740)
//! live in brep_offset_api_thru_sections_b.rs.
//!
//! OCCT inheritance chain (hxx L52): BRepOffsetAPI_ThruSections ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base members
//! (myShape, myGenerated, the Done flag) are plain fields (the Stage 2e
//! facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_Sequence/List -> Vec; NCollection_DataMap -> HashMap
//!    keyed by the TShape identity; NCollection_List<int> -> Vec<i32>.
//! 2. The engines BRepFill_CompatibleWires / BRepFill_Generator are the
//!    D1-approved brep_fill translations; they consume the rcad BRep pool
//!    (architecture difference #4) — the class holds my_brep.
//! 3. TopExp::Vertices(E, V1, V2) -> the brep_fill::generator
//!    top_exp_vertices oriented access; TopExp::FirstVertex/LastVertex(E)
//!    (raw, no CumOri) -> the local first_vertex/last_vertex reads.
//! 4. BRepTools_WireExplorer -> the brep_offset_inter2d re-host; the
//!    CurrentVertex() access is the oriented last vertex of the current
//!    edge (the wexp traversal form).
//! 5. BRepBuilderAPI_FindPlane / BRepLib::EncodeRegularity (TKTopAlgo /
//!    TKBRep) have no rcad translation yet — the FindPlane carrier keeps
//!    the OCCT Found()=false exit (the PerformPlan flow), EncodeRegularity
//!    is the GAP leaf (reported gap).
//! 6. BRepClass3d_SolidClassifier -> the landed solid_classifier (the
//!    TopAbs state u8: 0 = IN).
//! 7. TopoDS_Shell::Closed flag -> the shape flag read/write carriers.

use std::collections::HashMap;

use rcad_kernel::geom::{Plane, Surface3};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool::{
    builder_add_face_wire, builder_add_wire_edge, builder_set_closed,
    explorer, reversed, shape_is_closed, sub_shapes,
};
use crate::brep_fill::compatible_wires::CompatibleWires;
use crate::brep_fill::generator::{
    top_exp_vertices, BRepFillGenerator, BRepFillThruSectionErrorStatus,
};
use crate::offset::brep_offset_inter2d::BRepToolsWireExplorer;
use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;

use super::brep_offset_api_thru_sections_b::GeomBSplineSurface;

/// OCCT TopAbs_State::TopAbs_IN (TopAbs_State.hxx).
pub const TOPABS_IN: u8 = 0;

// ---------------------------------------------------------------------------
// GAP carriers (architecture difference #5).
// ---------------------------------------------------------------------------

/// OCCT BRepBuilderAPI_FindPlane (TKTopAlgo/BRepBuilderAPI_FindPlane.hxx) —
/// the plane finder of PerformPlan (architecture difference #5; GAP: no rcad
/// translation yet — Found() = false carries the OCCT no-plane exit of the
/// PerformPlan flow; the constructor keeps the OCCT form).
pub struct BRepBuilderAPIFindPlane;

impl BRepBuilderAPIFindPlane {
    /// OCCT BRepBuilderAPI_FindPlane(W, Tol).
    pub fn new(_the_w: &Shape, _the_tol: f64) -> Self {
        BRepBuilderAPIFindPlane
    }

    /// OCCT BRepBuilderAPI_FindPlane::Found() — GAP (false).
    pub fn found(&self) -> bool {
        false
    }

    /// OCCT BRepBuilderAPI_FindPlane::Plane() — GAP.
    pub fn plane(&self) -> Plane {
        panic!("GAP: BRepBuilderAPI_FindPlane::Plane (TKTopAlgo not translated)")
    }
}

/// OCCT BRepLib::EncodeRegularity(S) (TKBRep/BRepLib) — GAP leaf (the
/// chfi2d_builder.rs precedent).
fn brep_lib_encode_regularity(_the_s: &Shape) {
    panic!("GAP: BRepLib::EncodeRegularity (TKBRep/BRepLib not translated)")
}

/// OCCT BRep_Tool::Degenerated(E).
pub(crate) fn brep_tool_degenerated(the_e: &Shape) -> bool {
    match the_e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT TopExp::FirstVertex(E) — the raw first vertex (no CumOri).
pub(crate) fn first_vertex(brep: &BRep, the_e: &Shape) -> Shape {
    brep.edge(the_e.clone()).first.clone()
}

/// OCCT TopExp::LastVertex(E) — the raw last vertex (no CumOri).
pub(crate) fn last_vertex(brep: &BRep, the_e: &Shape) -> Shape {
    brep.edge(the_e.clone()).last.clone()
}

// ---------------------------------------------------------------------------
// BRep_Builder carriers (the tool.rs style).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::MakeSolid(S) + Add(Solid, Shell).
fn brep_builder_make_solid_add_shell(brep: &mut BRep, the_shell: &Shape) -> Shape {
    brep.add_tsolid(vec![the_shell.clone()])
}

/// OCCT BRep_Builder::MakeFace(F, S, Tol) + Add(F, W) — the rcad face takes
/// the outer wire at creation (the generator.rs architecture note).
fn brep_builder_make_face_from_surface(
    brep: &mut BRep,
    surface: &Surface3,
    w: &Shape,
) -> Shape {
    brep.add_tface(Some(surface.clone()), w.clone(), vec![], None, None, vec![], true)
}

// ---------------------------------------------------------------------------
// OCCT statics.
// ---------------------------------------------------------------------------

/// OCCT static PreciseUpar(anUpar, aSurface) (cxx L61-79) — pins the
/// u-parameter of the surface close to a U-knot to this U-knot.
pub(crate) fn precise_upar(an_upar: f64, a_surface: &GeomBSplineSurface) -> f64 {
    // OCCT L63: constexpr double Tol = Precision::PConfusion().
    let tol = rcad_kernel::precision::PCONFUSION;
    // OCCT L66-71.
    let mut i1 = 0i32;
    let mut i2 = 0i32;
    a_surface.locate_u(an_upar, tol, &mut i1, &mut i2);
    let u1 = a_surface.u_knot(i1);
    let u2 = a_surface.u_knot(i2);

    // OCCT L73-78.
    if an_upar - u1 < u2 - an_upar {
        u1
    } else {
        u2
    }
}

/// OCCT static PerformPlan(W, presPln, theFace) (cxx L84-120) — constructs
/// a plane of filling if it exists.
pub(crate) fn perform_plan(
    brep: &mut BRep,
    w: &Shape,
    pres_pln: f64,
    the_face: &mut Shape,
) -> bool {
    // OCCT L86-97: the degenerated check.
    let mut is_degen = true;
    for value in sub_shapes(w) {
        let an_edge = value;
        if !brep_tool_degenerated(&an_edge) {
            is_degen = false;
        }
    }
    if is_degen {
        return true;
    }

    // OCCT L99-119.
    let mut ok = false;
    if !w.is_null() {
        // OCCT L102-106: BRepBuilderAPI_FindPlane Searcher(W, presPln).
        let searcher = BRepBuilderAPIFindPlane::new(w, pres_pln);
        if searcher.found() {
            // OCCT L105: theFace = BRepBuilderAPI_MakeFace(Searcher.Plane(),
            // W).
            let plane = searcher.plane();
            *the_face =
                brep_builder_make_face_from_surface(brep, &Surface3::Plane(plane), w);
            ok = true;
        } else {
            // OCCT L108-115: try to find another surface.
            let mf = super::brep_offset_api_thru_sections_b::BRepLibMakeFaceWire::from_wire(w);
            if mf.is_done() {
                *the_face = mf.face();
                ok = true;
            }
        }
    }

    ok
}

/// OCCT static IsSameOriented(aFace, aShell) (cxx L125-153) — checks whether
/// aFace is oriented to the same side as aShell.
pub(crate) fn is_same_oriented(a_face: &Shape, a_shell: &Shape) -> bool {
    // OCCT L127-131.
    let explo = explorer(a_face, ShapeType::Edge, ShapeType::Shape);
    let an_edge = explo.first().cloned().unwrap_or_else(Shape::null);
    let or1 = an_edge.orientation;

    // OCCT L133-135: EFmap — MapShapesAndAncestors(shell, EDGE, FACE).
    let mut ef_map = indexmap::IndexMap::new();
    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
        a_shell,
        ShapeType::Edge,
        ShapeType::Face,
        &mut ef_map,
    );

    // OCCT L137-145.
    let adjacent_face = ef_map
        .get(&(an_edge.ptr_id(), an_edge.location))
        .and_then(|(_, l)| l.first().cloned())
        .unwrap_or_else(Shape::null);
    let mut the_edge = Shape::null();
    for e in explorer(&adjacent_face, ShapeType::Edge, ShapeType::Shape) {
        the_edge = e;
        if the_edge.is_same(&an_edge) {
            break;
        }
    }

    // OCCT L147-152.
    let or2 = the_edge.orientation;
    or1 != or2
}

/// OCCT static MakeSolid(shell, wire1, wire2, presPln, face1, face2)
/// (cxx L156-269).
pub(crate) fn make_solid(
    brep: &mut BRep,
    shell: &mut Shape,
    wire1: &Shape,
    wire2: &Shape,
    pres_pln: f64,
    face1: &mut Shape,
    face2: &mut Shape,
) -> Shape {
    // OCCT L158-160.
    if shell.is_null() {
        panic!("Thrusections is not build");
    }
    // OCCT L161: bool B = shell.Closed().
    let mut b = shape_is_closed(shell);

    // OCCT L163-190.
    if !b {
        // It is necessary to close the extremities.
        b = perform_plan(brep, wire1, pres_pln, face1);
        if b {
            b = perform_plan(brep, wire2, pres_pln, face2);
            if b {
                if !face1.is_null() && !is_same_oriented(face1, shell) {
                    *face1 = reversed(face1);
                }
                if !face2.is_null() && !is_same_oriented(face2, shell) {
                    *face2 = reversed(face2);
                }

                if !face1.is_null() {
                    builder_add_face_wire(shell, face1);
                }
                if !face2.is_null() {
                    builder_add_face_wire(shell, face2);
                }

                builder_set_closed(shell, true);
            }
        }
    }

    // OCCT L192-198.
    let mut solid = brep_builder_make_solid_add_shell(brep, shell);

    // OCCT L200-207: verify the orientation of the solid.
    let mut clas3d = SolidClassifier::from_shape(&solid);
    clas3d.perform_infinite_point(CONFUSION);
    if clas3d.state() == TOPABS_IN {
        let local = reversed(shell);
        solid = brep_builder_make_solid_add_shell(brep, &local);
    }

    // OCCT L209-210.
    builder_set_closed(&mut solid, true);
    solid
}

/// OCCT BRepOffsetAPI_ThruSections (hxx L52-215).
pub struct BRepOffsetAPIThruSections {
    // OCCT BRepBuilderAPI base members.
    pub(crate) my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    pub(crate) my_shape: Shape,          // OCCT: myShape
    pub(crate) my_generated: Vec<Shape>, // OCCT: myGenerated
    // The rcad arena stand-in (architecture difference #2).
    pub(crate) my_brep: BRep, // rcad pool (arch. diff. #4)
    // OCCT private members (hxx L193-215 + the ctor defaults).
    pub(crate) my_input_wires: Vec<Shape>,                  // OCCT: myInputWires
    pub(crate) my_wires: Vec<Shape>,                        // OCCT: myWires
    pub(crate) my_edge_new_indices: HashMap<u64, Vec<i32>>, // OCCT: myEdgeNewIndices
    pub(crate) my_vertex_index: HashMap<u64, i32>,          // OCCT: myVertexIndex
    pub(crate) my_nb_edges_in_section: i32,                 // OCCT: myNbEdgesInSection
    pub(crate) my_is_solid: bool,                           // OCCT: myIsSolid
    pub(crate) my_is_ruled: bool,                           // OCCT: myIsRuled
    pub(crate) my_w_check: bool,                            // OCCT: myWCheck
    pub(crate) my_pres3d: f64,                              // OCCT: myPres3d
    pub(crate) my_first: Shape,                             // OCCT: myFirst
    pub(crate) my_last: Shape,                              // OCCT: myLast
    pub(crate) my_degen1: bool,                             // OCCT: myDegen1
    pub(crate) my_degen2: bool,                             // OCCT: myDegen2
    pub(crate) my_edge_face: HashMap<u64, Shape>,           // OCCT: myEdgeFace
    pub(crate) my_continuity: rcad_kernel::topods::GeomAbsShape, // OCCT: myContinuity
    pub(crate) my_crit_weights: [f64; 3],                   // OCCT: myCritWeights
    pub(crate) my_use_smoothing: bool,                      // OCCT: myUseSmoothing
    pub(crate) my_param_type: ApproxParametrizationType,    // OCCT: myParamType
    pub(crate) my_deg_max: i32,                             // OCCT: myDegMax
    pub(crate) my_status: BRepFillThruSectionErrorStatus,   // OCCT: myStatus
    pub(crate) my_bf_generator: Option<BRepFillGenerator>,  // OCCT: myBFGenerator
    pub(crate) my_mutable_input: bool,                      // OCCT: myMutableInput
}

/// OCCT Approx_ParametrizationType
/// (TKGeomBase/Approx_ParametrizationType.hxx L25-29).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxParametrizationType {
    IsoParametric,
    ChordLength,
    Centripetal,
}

impl BRepOffsetAPIThruSections {
    /// OCCT BRepOffsetAPI_ThruSections::BRepOffsetAPI_ThruSections(isSolid,
    /// ruled, pres3d) (cxx L271-292).
    pub fn new(is_solid: bool, ruled: bool, pres3d: f64) -> Self {
        BRepOffsetAPIThruSections {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_brep: BRep::new(),
            my_input_wires: Vec::new(),
            my_wires: Vec::new(),
            my_edge_new_indices: HashMap::new(),
            my_vertex_index: HashMap::new(),
            my_nb_edges_in_section: 0,
            my_is_solid: is_solid,
            my_is_ruled: ruled,
            my_w_check: true,
            my_pres3d: pres3d,
            my_first: Shape::null(),
            my_last: Shape::null(),
            my_degen1: false,
            my_degen2: false,
            my_edge_face: HashMap::new(),
            my_continuity: rcad_kernel::topods::GeomAbsShape::C2,
            my_crit_weights: [0.4, 0.2, 0.4],
            my_use_smoothing: false,
            my_param_type: ApproxParametrizationType::ChordLength,
            my_deg_max: 8,
            my_status: BRepFillThruSectionErrorStatus::NotDone,
            my_bf_generator: None,
            my_mutable_input: true,
        }
    }

    /// OCCT BRepOffsetAPI_ThruSections::Init(isSolid, ruled, pres3d) (cxx
    /// L294-314).
    pub fn init(&mut self, is_solid: bool, ruled: bool, pres3d: f64) {
        // OCCT L296-299.
        self.my_is_solid = is_solid;
        self.my_is_ruled = ruled;
        self.my_pres3d = pres3d;
        self.my_w_check = true;
        self.my_mutable_input = true;
        // OCCT L301-313.
        self.my_param_type = ApproxParametrizationType::ChordLength;
        self.my_deg_max = 6;
        self.my_continuity = rcad_kernel::topods::GeomAbsShape::C2;
        self.my_crit_weights[0] = 0.4;
        self.my_crit_weights[1] = 0.2;
        self.my_crit_weights[2] = 0.4;
        self.my_use_smoothing = false;
        self.my_status = BRepFillThruSectionErrorStatus::NotDone;
    }

    /// OCCT BRepOffsetAPI_ThruSections::AddWire(wire) (cxx L316-320).
    pub fn add_wire(&mut self, wire: &Shape) {
        self.my_wires.push(wire.clone());
        self.my_input_wires.push(wire.clone());
    }

    /// OCCT BRepOffsetAPI_ThruSections::AddVertex(aVertex) (cxx L322-341).
    pub fn add_vertex(&mut self, a_vertex: &Shape) {
        // OCCT L324-331: the degenerated edge.
        let mut deg_edge = crate::brep_algo::tool::builder_make_edge();
        crate::brep_algo::tool::builder_add_edge_vertex(
            &mut deg_edge,
            &crate::brep_algo::tool::oriented(a_vertex, Orientation::Forward),
        );
        crate::brep_algo::tool::builder_add_edge_vertex(
            &mut deg_edge,
            &crate::brep_algo::tool::oriented(a_vertex, Orientation::Reversed),
        );
        crate::brep_algo::tool::builder_set_degenerated(&mut deg_edge, true);

        // OCCT L333-336: the degenerated wire.
        let mut deg_wire = crate::brep_algo::tool::builder_make_wire();
        builder_add_wire_edge(&mut deg_wire, &deg_edge);
        builder_set_closed(&mut deg_wire, true);

        // OCCT L338-340.
        self.my_wires.push(deg_wire.clone());
        self.my_input_wires.push(deg_wire);
    }

    /// OCCT BRepOffsetAPI_ThruSections::CheckCompatibility(check) (cxx
    /// L343-346).
    pub fn check_compatibility(&mut self, check: bool) {
        self.my_w_check = check;
    }

    /// OCCT BRepOffsetAPI_ThruSections::Build(...) (cxx L348-540).
    pub fn build(&mut self) {
        // OCCT L349-350.
        self.my_status = BRepFillThruSectionErrorStatus::Done;
        self.my_bf_generator = None;
        // OCCT L351-366: the punctual-section configuration check — the
        // middle sections (i = 2 .. Length-1) must not be punctual.
        {
            let wires = self.my_wires.clone();
            let len = wires.len();
            let mut i = 2usize; // the OCCT 1-based i.
            while len >= 3 && i <= len - 1 {
                let mut wdeg = true;
                for an_edge in explorer(&wires[i - 1], ShapeType::Edge, ShapeType::Shape) {
                    wdeg = wdeg && brep_tool_degenerated(&an_edge);
                }
                if wdeg {
                    self.my_status = BRepFillThruSectionErrorStatus::WrongUsage;
                    return;
                }
                i += 1;
            }
        }
        // OCCT L367-378: with at most two sections, at least one edge of the
        // set must not be degenerated (the conjunction over all sections).
        if self.my_wires.len() <= 2 {
            let mut wdeg = true;
            for w in self.my_wires.iter() {
                for an_edge in explorer(w, ShapeType::Edge, ShapeType::Shape) {
                    wdeg = wdeg && brep_tool_degenerated(&an_edge);
                }
            }
            if wdeg {
                self.my_status = BRepFillThruSectionErrorStatus::WrongUsage;
                return;
            }
        }

        // OCCT L381.
        self.my_nb_edges_in_section = 0;

        if self.my_w_check {
            // OCCT L383-390: compute origin and orientation on wires to
            // avoid twisted results and update wires to have the same number
            // of edges — BRepFill_CompatibleWires.
            let my_wires = std::mem::take(&mut self.my_wires);
            let mut georges = CompatibleWires::new_with_sections(&self.my_brep, &my_wires);
            georges.perform(&mut self.my_brep, true);
            #[allow(unused_assignments)]
            let mut working_sections: Vec<Shape> = Vec::new();
            if georges.is_done() {
                // OCCT L391-396.
                working_sections = georges.shape().to_vec();
                let working_map = georges.generated().clone();
                let _ = working_map;
                self.my_degen1 = georges.is_degenerated_first_section();
                self.my_degen2 = georges.is_degenerated_last_section();
                // OCCT L397-403.
                let mut ind_first_sec = 1usize;
                if georges.is_degenerated_first_section() {
                    ind_first_sec = 2;
                }
                let mut a_working_section =
                    working_sections[ind_first_sec - 1].clone();
                self.my_nb_edges_in_section +=
                    sub_shapes(&a_working_section).len() as i32;
                // OCCT L404-436: for each sub-edge of each section, save its
                // splits.
                for ii in 1..=my_wires.len() {
                    let my_wire = my_wires[ii - 1].clone();
                    for an_edge in sub_shapes(&my_wire) {
                        let mut a_sign = 1i32;
                        let (vfirst, vlast) = top_exp_vertices(&self.my_brep, &an_edge);
                        let a_new_edges = georges.generated_shapes(&an_edge);
                        let mut i_list: Vec<i32> = Vec::new();
                        a_working_section = working_sections[ii - 1].clone();
                        let nb_new_edges = a_new_edges.len();
                        for (kk0, a_new_edge) in a_new_edges.iter().enumerate() {
                            let kk = kk0 + 1;
                            let mut inde = 1usize;
                            let mut wexp = BRepToolsWireExplorer::new();
                            wexp.init(&a_working_section, &Shape::null());
                            while wexp.more() {
                                let a_working_edge = wexp.current();
                                if a_working_edge.is_same(a_new_edge) {
                                    a_sign =
                                        if a_working_edge.orientation == Orientation::Forward {
                                            1
                                        } else {
                                            -1
                                        };
                                    break;
                                }
                                wexp.next();
                                inde += 1;
                            }
                            i_list.push(inde as i32);
                            if kk == 1 || kk == nb_new_edges {
                                // OCCT L426-434: for each sub-vertex, save the
                                // index of the new edge.
                                let (new_vfirst, new_vlast) =
                                    top_exp_vertices(&self.my_brep, a_new_edge);
                                if new_vfirst.is_same(&vfirst)
                                    && !self.my_vertex_index.contains_key(&vfirst.ptr_id())
                                {
                                    self.my_vertex_index
                                        .insert(vfirst.ptr_id(), a_sign * inde as i32);
                                }
                                if new_vlast.is_same(&vlast)
                                    && !self.my_vertex_index.contains_key(&vlast.ptr_id())
                                {
                                    self.my_vertex_index
                                        .insert(vlast.ptr_id(), a_sign * (-(inde as i32)));
                                }
                            }
                        }
                        self.my_edge_new_indices.insert(an_edge.ptr_id(), i_list);
                    }
                }
            } else {
                // OCCT L443-448.
                self.my_status = georges.get_status();
                self.my_done = false;
                self.my_wires = my_wires;
                return;
            }

            // OCCT L451.
            self.my_wires = working_sections;
        } else {
            // OCCT L452-481: no check.
            let mut an_edge: Shape = Shape::null();
            for ii in 1..=self.my_wires.len() {
                let mut inde = 1usize;
                for an_edge_e in
                    explorer(&self.my_wires[ii - 1], ShapeType::Edge, ShapeType::Shape)
                {
                    an_edge = an_edge_e;
                    let mut i_list: Vec<i32> = Vec::new();
                    i_list.push(inde as i32);
                    self.my_edge_new_indices.insert(an_edge.ptr_id(), i_list);
                    let (v1, v2) = top_exp_vertices(&self.my_brep, &an_edge);
                    if !self.my_vertex_index.contains_key(&v1.ptr_id()) {
                        self.my_vertex_index.insert(v1.ptr_id(), inde as i32);
                    }
                    if !self.my_vertex_index.contains_key(&v2.ptr_id()) {
                        self.my_vertex_index.insert(v2.ptr_id(), -(inde as i32));
                    }
                    inde += 1;
                }
                let inde = inde - 1;
                if inde as i32 > self.my_nb_edges_in_section {
                    self.my_nb_edges_in_section = inde as i32;
                }
                if inde == 1 && brep_tool_degenerated(&an_edge) {
                    if ii == 1 {
                        self.my_degen1 = true;
                    } else {
                        self.my_degen2 = true;
                    }
                }
            }
        }

        // OCCT L483-495: the try block (the catch Standard_Failure is the
        // Err branch — the GAP panics take that exit).
        if self.my_wires.len() == 2 || self.my_is_ruled {
            // OCCT L489: create a ruled shell.
            self.create_ruled();
        } else {
            // OCCT L491: create a smoothed shell.
            self.create_smoothed();
        }

        // OCCT L499-503.
        if self.my_status != BRepFillThruSectionErrorStatus::Done {
            self.my_done = false;
            return;
        }
        // OCCT L504-505: Encode the Regularities.
        brep_lib_encode_regularity(&self.my_shape.clone());
    }

    /// OCCT BRepBuilderAPI_Command::IsDone() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_Command.hxx: Standard_Boolean IsDone() const);
    /// the rcad field keeps the crate visibility, the accessor restores the
    /// OCCT-public surface (form restoration, no behavior change).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// Test-world input bridge for the rcad BRep-pool architecture
    /// difference #4: the engines read TShape-internal reference indices
    /// against this pool (top_exp_vertices(&my_brep, …), …), so input
    /// wire/vertex shapes must be built natively in it.  OCCT has no
    /// counterpart (a TopoDS_Shape carries no pool); no algorithm behavior
    /// is changed — the accessor only exposes the existing my_brep member.
    pub fn brep_mut(&mut self) -> &mut BRep {
        &mut self.my_brep
    }

    /// OCCT BRepBuilderAPI_MakeShape::Shape() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_MakeShape.hxx: Standard_EXPORT const TopoDS_Shape&
    /// Shape() const; raises StdFail_NotDone when not done).
    pub fn shape(&self) -> Shape {
        assert!(
            self.my_done,
            "StdFail_NotDone: BRepOffsetAPI_ThruSections::Shape()"
        );
        self.my_shape.clone()
    }

    /// Test-world extraction bridge (the FilletResult.brep pattern,
    /// algo_ext::topods_ext::extract_result_brep): flattens the
    /// (my_brep arena, my_shape root) pair into the self-contained BRep the
    /// test world consumes (StepWriter / total_surface_area).  OCCT has no
    /// equivalent (the TopoDS_Shape carries its arena implicitly) — the
    /// rcad BRep-pool architecture difference #4 glue.
    pub fn result_brep(&mut self) -> Option<BRep> {
        if !self.my_done || self.my_shape.is_null() {
            return None;
        }
        let locations = self.my_brep.locations.clone();
        Some(crate::algo_ext::topods_ext::extract_result_brep(
            &self.my_shape,
            locations,
        ))
    }

    /// OCCT BRepOffsetAPI_ThruSections::CreateRuled() (cxx L542-702).
    pub(crate) fn create_ruled(&mut self) {
        // OCCT L543-546.
        let nb_sects = self.my_wires.len();
        let mut the_generator = BRepFillGenerator::new();
        the_generator.set_mutable_input(self.my_mutable_input);
        // OCCT L548-551.
        for i in 0..nb_sects {
            the_generator.add_wire(self.my_wires[i].clone());
        }
        the_generator.perform(&mut self.my_brep);
        // OCCT L552-557.
        let a_status = the_generator.get_status();
        if a_status != BRepFillThruSectionErrorStatus::Done {
            self.my_status = a_status;
            self.my_bf_generator = Some(the_generator);
            return;
        }
        self.my_bf_generator = Some(the_generator);
        let my_bf_generator = self.my_bf_generator.as_mut().unwrap();
        // OCCT L558.
        let shell = my_bf_generator.shell().unwrap_or_else(Shape::null);
        let shell_for_history = shell.clone();

        if self.my_is_solid {
            // OCCT L562-563: check if the first wire is the same as the last.
            let v_closed =
                nb_sects > 0 && self.my_wires[0].is_same(&self.my_wires[nb_sects - 1]);

            if v_closed {
                // OCCT L565-580.
                let mut solid = brep_builder_make_solid_add_shell(&mut self.my_brep, &shell);

                // verify the orientation of the solid.
                let mut clas3d = SolidClassifier::from_shape(&solid);
                clas3d.perform_infinite_point(CONFUSION);
                if clas3d.state() == TOPABS_IN {
                    let local = reversed(&shell);
                    solid = brep_builder_make_solid_add_shell(&mut self.my_brep, &local);
                }
                self.my_shape = solid;
            } else {
                // OCCT L582-590: myBFGenerator stores the same myWires.
                let first = self.my_wires.first().cloned().unwrap_or_else(Shape::null);
                let last = self.my_wires.last().cloned().unwrap_or_else(Shape::null);
                let wire1 = my_bf_generator.result_shape(&mut self.my_brep, &first);
                let wire2 = my_bf_generator.result_shape(&mut self.my_brep, &last);

                let mut shell_mut = shell.clone();
                let mut face1 = self.my_first.clone();
                let mut face2 = self.my_last.clone();
                self.my_shape = make_solid(
                    &mut self.my_brep,
                    &mut shell_mut,
                    &wire1,
                    &wire2,
                    self.my_pres3d,
                    &mut face1,
                    &mut face2,
                );
                self.my_first = face1;
                self.my_last = face2;
            }

            // OCCT L594: Done().
            self.my_done = true;
        } else {
            // OCCT L598-602.
            self.my_shape = shell;
            self.my_done = true;
        }

        let shell = shell_for_history;
        // OCCT L604-702: the history.
        let mut an_exp1 = BRepToolsWireExplorer::new();
        let mut an_exp2 = BRepToolsWireExplorer::new();
        let mut m = indexmap::IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &shell,
            ShapeType::Edge,
            ShapeType::Face,
            &mut m,
        );

        let mut mv = indexmap::IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &shell,
            ShapeType::Vertex,
            ShapeType::Face,
            &mut mv,
        );

        for i in 1..nb_sects {
            let wire1 = self.my_wires[i - 1].clone();
            let wire2 = self.my_wires[i].clone();

            an_exp1.init(&wire1, &Shape::null());
            an_exp2.init(&wire2, &Shape::null());

            let mut tantque = an_exp1.more() && an_exp2.more();

            while tantque {
                let edge1 = an_exp1.current();
                let edge2 = an_exp2.current();
                let degen1 = brep_tool_degenerated(&edge1);
                let degen2 = brep_tool_degenerated(&edge2);

                // OCCT L621-641: MapFaces.
                let mut map_faces: HashMap<u64, Shape> = HashMap::new();
                if degen2 {
                    let vdegen = first_vertex(&self.my_brep, &edge2);
                    let vdegen = my_bf_generator.result_shape(&mut self.my_brep, &vdegen);
                    if let Some((_, l)) =
                        mv.get(&(vdegen.ptr_id(), vdegen.location))
                    {
                        for value in l.clone() {
                            map_faces.insert(value.ptr_id(), value);
                        }
                    }
                } else {
                    let res = my_bf_generator.result_shape(&mut self.my_brep, &edge2);
                    if let Some((_, l)) = m.get(&(res.ptr_id(), res.location)) {
                        for value in l.clone() {
                            map_faces.insert(value.ptr_id(), value);
                        }
                    }
                }

                if degen1 {
                    let vdegen = first_vertex(&self.my_brep, &edge1);
                    let vdegen = my_bf_generator.result_shape(&mut self.my_brep, &vdegen);
                    if let Some((_, l)) =
                        mv.get(&(vdegen.ptr_id(), vdegen.location))
                    {
                        for face in l.clone() {
                            if map_faces.contains_key(&face.ptr_id()) {
                                self.my_edge_face.insert(edge1.ptr_id(), face);
                                break;
                            }
                        }
                    }
                } else {
                    let res = my_bf_generator.result_shape(&mut self.my_brep, &edge1);
                    if let Some((_, l)) = m.get(&(res.ptr_id(), res.location)) {
                        for face in l.clone() {
                            if map_faces.contains_key(&face.ptr_id()) {
                                self.my_edge_face.insert(edge1.ptr_id(), face);
                                break;
                            }
                        }
                    }
                }

                // OCCT L667-684: the traversal.
                if !degen1 {
                    an_exp1.next();
                }
                if !degen2 {
                    an_exp2.next();
                }

                tantque = an_exp1.more() && an_exp2.more();
                if degen1 {
                    tantque = an_exp2.more();
                }
                if degen2 {
                    tantque = an_exp1.more();
                }
            }
        }
    }
}

impl Default for BRepOffsetAPIThruSections {
    fn default() -> Self {
        // The OCCT constructor requires (isSolid, ruled, pres3d); the
        // Default carries the OCCT doc defaults (false, false, 1.0e-6).
        Self::new(false, false, 1.0e-6)
    }
}
