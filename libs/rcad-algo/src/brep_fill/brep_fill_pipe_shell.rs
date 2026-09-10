//! OCCT BRepFill_PipeShell (TKBool/BRepFill) — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill_PipeShell.cxx
//! (L87-1681) + BRepFill_PipeShell.hxx (L50-260).
//!
//! First consumer: BRepOffsetAPI_MakePipeShell (Stage 2e carrier).
//!
//! Residual gaps (annotated at the call sites, carried by the split unit
//! [`super::brep_fill_pipe_shell_b`]): BRepFill_Section, BRepFill::SearchOrigin,
//! BRepAdaptor_CompCurve, BRepBuilderAPI_Transform/Copy,
//! IntCurveSurface_HInter host tools, BRepGProp::LinearProperties and
//! BRepFill_Sweep::Build (the part-B backlog of BRepFill_Sweep.cxx).
//! The GeomFill sweep machinery (CurveAndTrihedron, Fixed, ConstantBiNormal,
//! DiscreteTrihedron GAP aside, GuideTrihedronAC/Plan, LocationGuide), the
//! BRepFill law/placement family and Law_Interpol are the landed
//! translations — the carriers below are the real classes.
//!
//! Architecture notes:
//! - `NCollection_Sequence` -> `Vec`; `NCollection_List<TopoDS_Shape>` ->
//!   `Vec<Shape>`; `NCollection_DataMap(TopoDS_Shape, ...)` ->
//!   `HashMap<ShapeKey, ...>`; `NCollection_Map` -> `HashSet<ShapeKey>`.
//! - `handle<X>` members map to `Option<X>` when OCCT nullifies them
//!   (myLocation / mySection / myLaw); `handle(BRepFill_LocationLaw)` /
//!   `handle(BRepFill_SectionLaw)` map to the `Rc<RefCell<dyn ...Ops>>` slots
//!   (the trait slots of brep_fill_location_law.rs / brep_fill_section_law.rs);
//!   `handle(Law_Function)` maps to [`LawFunctionHandle`].
//! - `occ::down_cast<GeomFill_LocationGuide>(myLocation->Law(i))` maps to the
//!   [`LocationLaw::as_any`] downcast (BRepFill_PipeShell.cxx L489/L1219/L1271).
//! - `TopExp::Vertices(E, V1, V2)` maps to the traversal-order
//!   `generator::top_exp_vertices`; `BRepTools_WireExplorer` maps to the
//!   wire-ordered `TWireData::edges` list; `CurrentVertex` is the
//!   traversal-first vertex of the current edge
//!   (BRepTools_WireExplorer.cxx L381/395/739-742 mapping).
//! - OCCT quirk kept: Set(false) assigns `GeomFill_IsFrenet` in BOTH branches
//!   (cxx L262 and L267).

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::Arc;

use glam::{DAffine3, DVec3};

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::Curve3;
use rcad_kernel::math::gp::Trsf;
use rcad_kernel::topo::topods::{
    tshape_flags, BRep, GeomAbsShape, Orientation, Shape, ShapeType, TShape,
};

use crate::brep_fill::brep_fill_axe::BRepLibMakeFaceWire;
use crate::brep_fill::brep_fill_edge3d_law::{
    BRepFillACRLaw, BRepFillEdge3DLaw, BRepFillEdgeOnSurfLaw,
};
use crate::brep_fill::brep_fill_location_law::{BRepFillLocationLawOps};
use crate::brep_fill::brep_fill_nsections::BRepFillNSections;
use crate::brep_fill::brep_fill_section_law::BRepFillSectionLawOps;
use crate::brep_fill::brep_fill_section_placement::BRepFillSectionPlacement;
use crate::brep_fill::brep_fill_shape_law::BRepFillShapeLaw;
use crate::brep_fill::brep_fill_sweep::BRepFillSweep;
use crate::brep_fill::compatible_wires::CompatibleWires as BRepFillCompatibleWires;
use crate::brep_fill::generator::{
    shape_key, shape_reversed, top_exp_vertices, top_exp_wire_vertices, ShapeKey,
};
use crate::brep_fill::offset_wire::append_edge_to_wire;
use crate::brep_fill::offset_wire_b::{
    explored_children, is_closed_wire, set_wire_closed, vertex_point, wire_edges,
};
use crate::geomalgo::geomfill::constant_bi_normal::ConstantBiNormal;
use crate::geomalgo::geomfill::corrected_frenet::CorrectedFrenet;
use crate::geomalgo::geomfill::curve_and_trihedron::CurveAndTrihedron;
use crate::geomalgo::geomfill::fixed::Fixed;
use crate::geomalgo::geomfill::frenet::Frenet;
use crate::geomalgo::geomfill::guide_trihedron_ac::GuideTrihedronAC;
use crate::geomalgo::geomfill::guide_trihedron_plan::GuideTrihedronPlan;
use crate::geomalgo::geomfill::location_guide::LocationGuide;
use crate::geomalgo::geomfill::location_law::LocationLaw;
use crate::geomalgo::geomfill::trihedron_law::{PipeError, TrihedronLaw};
use crate::geomalgo::law::law_function::LawFunctionHandle;
use crate::geomalgo::law::law_interpol::LawInterpol;

use super::brep_fill_pipe_shell_b::{
    brep_fill_search_origin, brep_gprop_linear_properties_centre, BRepBuilderAPICopy,
    BRepBuilderAPITransform, BRepAdaptorCompCurve, BRepFillSection, BRepFillTransitionStyle,
    BRepFillTypeOfContact, GeomFillDiscreteTrihedron, GeomFillTrihedron, IntCurveSurfaceHInter,
    ShapeSet, ShapeToArray2Map, GEOM_FILL_LOCATION,
};

/// OCCT GeomFill_PipeOk (mapped to the rcad PipeError variant).
const PIPE_OK: PipeError = PipeError::PipeOk;
/// OCCT GeomFill_PipeNotOk.
const PIPE_NOT_OK: PipeError = PipeError::PipeNotOk;
/// OCCT RealLast().
const REAL_LAST: f64 = f64::MAX;

// ===========================================================================
// File statics (BRepFill_PipeShell.cxx L87-220)
// ===========================================================================

/// OCCT gp_Trsf assignment from the glam affine form — the reverse of
/// `rcad_kernel::math::gp::Trsf::to_daffine3`.  The real
/// BRepFill_SectionPlacement::Transformation returns the glam affine while
/// the downstream carriers (myTrsfs of BRepFill_NSections) keep the OCCT
/// gp_Trsf member split (matrix / loc / scale).
fn trsf_from_daffine3(t: &DAffine3) -> Trsf {
    let mut r = Trsf::identity();
    let scale = t.x_axis.length();
    let cols = [t.x_axis, t.y_axis, t.z_axis];
    for (j, col) in cols.iter().enumerate() {
        r.matrix[0][j] = col.x / scale;
        r.matrix[1][j] = col.y / scale;
        r.matrix[2][j] = col.z / scale;
    }
    r.loc = t.w_axis;
    r.scale = scale;
    r
}

/// OCCT static bool ComputeSection(W1, W2, p1, p2, Wres) (L97-123) —
/// constructs an intermediary section.
fn compute_section(
    brep: &mut BRep,
    w1: &Shape,
    w2: &Shape,
    p1: f64,
    p2: f64,
    wres: &mut Shape,
) -> bool {
    // OCCT L103-110.
    let mut sr: Vec<f64> = Vec::new();
    let mut s_sh: Vec<Shape> = Vec::new();
    sr.clear();
    sr.push(0.0);
    sr.push(1.0);
    s_sh.clear();
    s_sh.push(w1.clone());
    s_sh.push(w2.clone());
    // OCCT L111-117.
    let mut cw = BRepFillCompatibleWires::new_with_sections(brep, &s_sh);
    cw.set_percent(0.1);
    cw.perform(brep, true);
    if !cw.is_done() {
        panic!("Uncompatible wires");
    }
    // OCCT L118-119.
    let empty_trsfs: Vec<Trsf> = Vec::new();
    let sl = Rc::new(RefCell::new(BRepFillNSections::new_with_params(
        brep,
        cw.shape().clone(),
        empty_trsfs,
        sr,
        0.0,
        1.0,
        true,
    )));
    // OCCT L120-121.
    let us = p1 / (p1 + p2);
    let mut w = Shape::null();
    sl.borrow_mut().d0(brep, us, &mut w);
    *wres = w;
    // OCCT L122.
    true
}

/// OCCT static void PerformTransition(Mode, Loc, angmin) (L130-146) —
/// modifies a law of location depending on Transition.
fn perform_transition(
    brep: &BRep,
    mode: BRepFillTransitionStyle,
    loc: &mut Option<Rc<RefCell<dyn BRepFillLocationLawOps>>>,
    angmin: f64,
) {
    if let Some(l) = loc {
        // OCCT L136.
        l.borrow_mut().base_mut().delete_transform();
        if mode == BRepFillTransitionStyle::Modified {
            // OCCT L139.
            l.borrow_mut().base_mut().transform_in_g0_law(brep);
        } else {
            // OCCT L143.
            l.borrow_mut().base_mut().transform_in_compatible_law(angmin);
        }
    }
}

/// OCCT static bool PerformPlan(S) (L153-183) — constructs a plane of filling
/// if exists.
fn perform_plan(brep: &mut BRep, s: &mut Shape) -> bool {
    // OCCT L155-164.
    let mut is_degen = true;
    for e in explored_children(brep, s, ShapeType::Edge) {
        if !brep.edge(e).degenerated {
            is_degen = false;
        }
    }
    if is_degen {
        // OCCT L166-169.
        *s = Shape::null();
        return true;
    }

    // OCCT L171-182.
    let w = s.clone();
    let mut ok = false;
    if !w.is_null() {
        let mkplan = BRepLibMakeFaceWire::new(brep, &w, true);
        if mkplan.is_done() {
            *s = mkplan.face();
            ok = true;
        }
    }
    ok
}

/// OCCT static bool IsSameOriented(theFace, theShell) (L187-220).
fn is_same_oriented(brep: &BRep, the_face: &Shape, the_shell: &Shape) -> bool {
    // OCCT L189-191: anEdgeFaceMap — EDGE -> FACE ancestors of the shell.
    let an_edge_face_map = map_shapes_and_ancestors_ef(brep, the_shell);

    // OCCT L193-214.
    for a_current_edge in explored_children(brep, the_face, ShapeType::Edge) {
        let found = an_edge_face_map
            .iter()
            .find(|(k, _)| k.ptr_id() == a_current_edge.ptr_id());
        let Some((_, an_adj_faces)) = found else {
            // OCCT L199-200.
            continue;
        };

        let a_current_edge_orientation = a_current_edge.orientation;
        let an_adj_face = an_adj_faces
            .first()
            .cloned()
            .expect("IsSameOriented: empty ancestor list");
        for an_adj_edge in explored_children(brep, &an_adj_face, ShapeType::Edge) {
            // IsSame compares TShape and Location only; orientation may differ.
            if an_adj_edge.ptr_id() == a_current_edge.ptr_id() {
                return a_current_edge_orientation != an_adj_edge.orientation;
            }
        }
    }

    // No shared edge found. Orientation cannot be determined reliably and this
    // should not happen for valid input. Preserve the current face orientation
    // instead of reporting an opposite orientation and forcing a reversal.
    // OCCT L219.
    true
}

/// OCCT static bool BuildBoundaries(theSweep, theSection, theBottom, theTop)
/// (L1619-1681).
fn build_boundaries(
    brep: &mut BRep,
    the_sweep: &BRepFillSweep,
    the_section: &Rc<RefCell<dyn BRepFillSectionLawOps>>,
    the_bottom: &mut Shape,
    the_top: &mut Shape,
) -> bool {
    // OCCT L1625-1629.
    let a_bottom_wire = brep.add_twire(vec![]);
    let a_top_wire = brep.add_twire(vec![]);
    // OCCT L1630-1634.
    let mut bfoundbottom = false;
    let mut bfoundtop = false;
    let a_v_edges = the_sweep.sections().expect("myVEdges");
    let mut b_all_same = true;

    // OCCT L1636-1657.
    for i in 1..=the_section.borrow().base().nb_law() as usize {
        // OCCT LowerCol() / UpperCol() of the HArray2.
        let lower = 0usize;
        let upper = a_v_edges[i - 1].len() - 1;
        let a_bottom_edge = a_v_edges[i - 1][lower].clone();

        if !a_bottom_edge.is_null() && a_bottom_edge.shape_type() == ShapeType::Edge {
            append_edge_to_wire(brep, &a_bottom_wire, &a_bottom_edge);
            bfoundbottom = true;
        }
        let a_top_edge = a_v_edges[i - 1][upper].clone();

        if !a_top_edge.is_null() && a_top_edge.shape_type() == ShapeType::Edge {
            append_edge_to_wire(brep, &a_top_wire, &a_top_edge);
            bfoundtop = true;
        }

        if !a_bottom_edge.is_null()
            && !a_top_edge.is_null()
            && a_bottom_edge.ptr_id() != a_top_edge.ptr_id()
        {
            b_all_same = false;
        }
    }

    // OCCT L1659-1663.
    if the_section.borrow().base().is_uclosed() {
        set_wire_closed(brep, &a_bottom_wire, true);
        set_wire_closed(brep, &a_top_wire, true);
    }

    // OCCT L1665-1672.
    if bfoundbottom {
        *the_bottom = a_bottom_wire;
    }

    if bfoundtop {
        *the_top = a_top_wire;
    }

    // OCCT L1675-1678.
    if b_all_same && bfoundbottom && bfoundtop {
        *the_top = the_bottom.clone();
    }

    // OCCT L1680.
    bfoundbottom || bfoundtop
}

/// OCCT TopExp::MapShapesAndAncestors(S, TopAbs_EDGE, TopAbs_FACE, Map) —
/// local mapping (exploration order, Vec of pairs as IndexedDataMap).
fn map_shapes_and_ancestors_ef(brep: &BRep, s: &Shape) -> Vec<(Shape, Vec<Shape>)> {
    let mut map: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for f in explored_children(brep, s, ShapeType::Face) {
        for w in explored_children(brep, &f, ShapeType::Wire) {
            for e in explored_children(brep, &w, ShapeType::Edge) {
                if let Some(pos) = map.iter().position(|(k, _)| k.ptr_id() == e.ptr_id()) {
                    map[pos].1.push(f.clone());
                } else {
                    map.push((e.clone(), vec![f.clone()]));
                }
            }
        }
    }
    map
}

/// OCCT TopExp::MapShapesAndAncestors(S, TopAbs_VERTEX, TopAbs_EDGE, Map) —
/// local mapping (used by BuildHistory on the tape shell).
fn map_shapes_and_ancestors_ve(brep: &BRep, s: &Shape) -> Vec<(Shape, Vec<Shape>)> {
    let mut map: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for w in explored_children(brep, s, ShapeType::Wire) {
        for e in explored_children(brep, &w, ShapeType::Edge) {
            let (vf, vl) = top_exp_vertices(brep, &e);
            for v in [vf, vl] {
                if let Some(pos) = map.iter().position(|(k, _)| k.ptr_id() == v.ptr_id()) {
                    if !map[pos].1.iter().any(|x| x.ptr_id() == e.ptr_id()) {
                        map[pos].1.push(e.clone());
                    }
                } else {
                    map.push((v.clone(), vec![e.clone()]));
                }
            }
        }
    }
    map
}

/// OCCT TopoDS_Iterator over the direct children of any type.
fn direct_children(brep: &BRep, s: &Shape) -> Vec<Shape> {
    match s.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        TShape::Face(fd) => {
            let mut c = vec![fd.outer_wire.clone()];
            c.extend(fd.inner_wires.iter().cloned());
            c
        }
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Solid(sd) => sd.shells.clone(),
        _ => Vec::new(),
    }
}

// ===========================================================================
// BRepFill_PipeShell (hxx L50-260 + cxx L224-1615)
// ===========================================================================

/// OCCT BRepFill_PipeShell — the general sweeping engine.
pub struct BRepFillPipeShell {
    /// OCCT: mySpine (hxx L231).
    my_spine: Shape,
    /// OCCT: myFirst (hxx L232).
    my_first: Shape,
    /// OCCT: myLast (hxx L233).
    my_last: Shape,
    /// OCCT: myShape (hxx L234).
    my_shape: Shape,
    /// OCCT: mySeq (hxx L235).
    my_seq: Vec<BRepFillSection>,
    /// OCCT: WSeq (hxx L236).
    w_seq: Vec<Shape>,
    /// OCCT: myIndOfSec (hxx L237).
    my_ind_of_sec: Vec<i32>,
    /// OCCT: myEdgeNewEdges (hxx L238-239).
    my_edge_new_edges: HashMap<ShapeKey, Vec<Shape>>,
    /// OCCT: myGenMap (hxx L240-241).
    my_gen_map: HashMap<ShapeKey, Vec<Shape>>,
    /// OCCT: myTol3d (hxx L242).
    my_tol3d: f64,
    /// OCCT: myBoundTol (hxx L243).
    my_bound_tol: f64,
    /// OCCT: myTolAngular (hxx L244).
    my_tol_angular: f64,
    /// OCCT: angmin (hxx L245).
    angmin: f64,
    /// OCCT: angmax (hxx L246).
    angmax: f64,
    /// OCCT: myMaxDegree (hxx L247).
    my_max_degree: i32,
    /// OCCT: myMaxSegments (hxx L248).
    my_max_segments: i32,
    /// OCCT: myForceApproxC1 (hxx L249).
    my_force_approx_c1: bool,
    /// OCCT: myLaw (hxx L250) — handle(Law_Function).
    my_law: Option<LawFunctionHandle>,
    /// OCCT: myIsAutomaticLaw (hxx L251).
    my_is_automatic_law: bool,
    /// OCCT: myLocation (hxx L252) — handle(BRepFill_LocationLaw).
    my_location: Option<Rc<RefCell<dyn BRepFillLocationLawOps>>>,
    /// OCCT: mySection (hxx L253) — handle(BRepFill_SectionLaw).
    my_section: Option<Rc<RefCell<dyn BRepFillSectionLawOps>>>,
    /// OCCT: myFaces (hxx L254).
    my_faces: Option<Vec<Vec<Shape>>>,
    /// OCCT: myTrihedron (hxx L255).
    my_trihedron: GeomFillTrihedron,
    /// OCCT: myTransition (hxx L256).
    my_transition: BRepFillTransitionStyle,
    /// OCCT: myStatus (hxx L257).
    my_status: PipeError,
    /// OCCT: myErrorOnSurf (hxx L258).
    my_error_on_surf: f64,
    /// OCCT: myIsBuildHistory (hxx L259).
    my_is_build_history: bool,
}

impl BRepFillPipeShell {
    /// OCCT BRepFill_PipeShell::BRepFill_PipeShell(Spine) (cxx L224-251).
    pub fn new(brep: &mut BRep, spine: &Shape) -> Self {
        // OCCT L226-231: member initialization.
        let mut r = BRepFillPipeShell {
            my_spine: spine.clone(),
            my_first: Shape::null(),
            my_last: Shape::null(),
            my_shape: Shape::null(),
            my_seq: Vec::new(),
            w_seq: Vec::new(),
            my_ind_of_sec: Vec::new(),
            my_edge_new_edges: HashMap::new(),
            my_gen_map: HashMap::new(),
            my_tol3d: 0.0,
            my_bound_tol: 0.0,
            my_tol_angular: 0.0,
            angmin: 0.0,
            angmax: 0.0,
            my_max_degree: 0,
            my_max_segments: 0,
            my_force_approx_c1: false,
            my_law: None,
            my_is_automatic_law: false,
            my_location: None,
            my_section: None,
            my_faces: None,
            my_trihedron: GeomFillTrihedron::IsCorrectedFrenet,
            my_transition: BRepFillTransitionStyle::Modified,
            my_status: PIPE_OK,
            my_error_on_surf: 0.0,
            my_is_build_history: true,
        };
        // OCCT L233-236.
        r.my_location = None;
        r.my_section = None;
        r.my_law = None;
        r.set_tolerance(1.0e-4, 1.0e-4, 1.0e-2);

        // OCCT L238-239.
        r.my_max_degree = 11;
        r.my_max_segments = 100;

        // Attention to closed non-declared wire !
        // OCCT L242-250: TopExp::Vertices(mySpine, Vf, Vl) — the wire form.
        if !is_closed_wire(brep, &r.my_spine) {
            let (vf, vl) = top_exp_wire_vertices(brep, &r.my_spine);
            if vf.ptr_id() == vl.ptr_id() {
                set_wire_closed(brep, &r.my_spine, true);
            }
        }
        r
    }

    /// OCCT BRepFill_PipeShell::Set(const bool IsFrenet) (cxx L257-273).
    pub fn set(&mut self, brep: &BRep, is_frenet: bool) {
        // OCCT L259-269.
        let t_law: Box<dyn TrihedronLaw> = if is_frenet {
            self.my_trihedron = GeomFillTrihedron::IsFrenet;
            Box::new(Frenet::new()) // OCCT new GeomFill_Frenet()
        } else {
            // OCCT quirk kept: the else branch also sets GeomFill_IsFrenet.
            self.my_trihedron = GeomFillTrihedron::IsFrenet;
            Box::new(CorrectedFrenet::new()) // OCCT new GeomFill_CorrectedFrenet()
        };
        // OCCT L270-272.
        let loc = Rc::new(RefCell::new(CurveAndTrihedron::new(t_law)));
        self.my_location = Some(Rc::new(RefCell::new(BRepFillEdge3DLaw::new(
            brep,
            &self.my_spine,
            loc,
        ))));
        self.my_section = None; // It is required to relocalize sections.
    }

    /// OCCT BRepFill_PipeShell::SetDiscrete() (cxx L279-289).
    pub fn set_discrete(&mut self, brep: &BRep) {
        // OCCT L281-284.
        self.my_trihedron = GeomFillTrihedron::IsDiscreteTrihedron;
        let t_law: Box<dyn TrihedronLaw> = GeomFillDiscreteTrihedron::new().into_trihedron_law();

        // OCCT L286-288.
        let loc = Rc::new(RefCell::new(CurveAndTrihedron::new(t_law)));
        self.my_location = Some(Rc::new(RefCell::new(BRepFillEdge3DLaw::new(
            brep,
            &self.my_spine,
            loc,
        ))));
        self.my_section = None; // It is required to relocalize sections.
    }

    /// OCCT BRepFill_PipeShell::Set(const gp_Ax2& Axe) (cxx L293-303).
    pub fn set_with_axe(&mut self, brep: &BRep, axe: rcad_kernel::math::gp::Ax2) {
        // OCCT L295-298.
        self.my_trihedron = GeomFillTrihedron::IsFixed;
        let v1 = axe.direction;
        let v2 = axe.x_direction;
        // OCCT L299-302.
        let t_law: Box<dyn TrihedronLaw> = Box::new(Fixed::new(v1, v2));
        let loc = Rc::new(RefCell::new(CurveAndTrihedron::new(t_law)));
        self.my_location = Some(Rc::new(RefCell::new(BRepFillEdge3DLaw::new(
            brep,
            &self.my_spine,
            loc,
        ))));
        self.my_section = None; // It is required to relocalize sections.
    }

    /// OCCT BRepFill_PipeShell::Set(const gp_Dir& BiNormal) (cxx L309-317).
    pub fn set_with_binormal(&mut self, brep: &BRep, binormal: DVec3) {
        // OCCT L311-315.
        self.my_trihedron = GeomFillTrihedron::IsConstantNormal;

        let t_law: Box<dyn TrihedronLaw> = Box::new(ConstantBiNormal::new(binormal));
        let loc = Rc::new(RefCell::new(CurveAndTrihedron::new(t_law)));
        self.my_location = Some(Rc::new(RefCell::new(BRepFillEdge3DLaw::new(
            brep,
            &self.my_spine,
            loc,
        ))));
        self.my_section = None; // Il faut relocaliser les sections.
    }

    /// OCCT BRepFill_PipeShell::Set(const TopoDS_Shape& SpineSupport)
    /// (cxx L323-337).
    pub fn set_spine_support(&mut self, brep: &BRep, spine_support: &Shape) -> bool {
        // OCCT L327-328: A special law of location is required.
        let loc = BRepFillEdgeOnSurfLaw::new(brep, &self.my_spine, spine_support);
        // OCCT L329.
        let b = loc.has_result();
        // OCCT L330-335.
        if b {
            self.my_location = Some(Rc::new(RefCell::new(loc)));
            self.my_trihedron = GeomFillTrihedron::IsDarboux;
            self.my_section = None; // It is required to relocalize the sections.
        }
        // OCCT L336.
        b
    }

    /// OCCT BRepFill_PipeShell::Set(AuxiliarySpine, CurvilinearEquivalence,
    /// KeepContact) (cxx L343-432).
    pub fn set_with_auxiliary_spine(
        &mut self,
        brep: &mut BRep,
        auxiliary_spine: &Shape,
        curvilinear_equivalence: bool,
        keep_contact: BRepFillTypeOfContact,
    ) {
        // Reorganization of the guide (pb of orientation and origin)
        // OCCT L348-350.
        let mut the_guide = auxiliary_spine.clone();
        let sp_close = is_closed_wire(brep, &self.my_spine);
        let guide_close = is_closed_wire(brep, auxiliary_spine);

        // OCCT L352-355.
        if keep_contact == BRepFillTypeOfContact::ContactOnBorder {
            self.my_is_automatic_law = true;
        }

        // OCCT L357-396.
        if !sp_close && !guide_close {
            // Case open reorientation of the guide
            let sp = self.my_spine.clone();
            let mut seq: Vec<Shape> = Vec::new();
            seq.push(sp);
            seq.push(the_guide.clone());
            let mut cw = BRepFillCompatibleWires::new_with_sections(brep, &seq);
            cw.set_percent(0.1);
            cw.perform(brep, true);
            if !cw.is_done() {
                panic!("Uncompatible wires");
            }
            the_guide = cw.shape()[1].clone();
        } else if guide_close {
            // Case guide closed : Determination of the origin
            // & reorientation of the guide
            let dir: DVec3;
            let sp_or: DVec3;
            if !sp_close {
                // OCCT L381-389.
                let (vf, vl) = top_exp_vertices(brep, &self.my_spine);
                let p = vertex_point(brep, &vf);
                let p_l = vertex_point(brep, &vl);
                let v = p - p_l;
                // OCCT L387: SpOr.BaryCenter(0.5, P, 0.5) — the midpoint.
                sp_or = p * 0.5 + p_l * 0.5;
                dir = v;
            } else {
                // OCCT L392-393.
                let bc = BRepAdaptorCompCurve::new(brep, &self.my_spine);
                let mut sp_or_v = DVec3::ZERO;
                let mut dir_v = DVec3::ZERO;
                bc.d1(0.0, &mut sp_or_v, &mut dir_v);
                sp_or = sp_or_v;
                dir = dir_v;
            }
            // OCCT L395.
            brep_fill_search_origin(brep, &mut the_guide, sp_or, dir, 100.0 * self.my_tol3d);
        }

        // transform the guide in a single curve
        // OCCT L399.
        let guide = BRepAdaptorCompCurve::new(brep, &the_guide);
        let guide_curve = guide.comp_curve();

        // OCCT L401-415.
        if curvilinear_equivalence {
            // trihedron by curvilinear reduced abscissa
            if keep_contact == BRepFillTypeOfContact::Contact
                || keep_contact == BRepFillTypeOfContact::ContactOnBorder
            {
                self.my_trihedron = GeomFillTrihedron::IsGuideACWithContact; // with rotation
            } else {
                self.my_trihedron = GeomFillTrihedron::IsGuideAC; // without rotation
            }

            let t_law = GuideTrihedronAC::new(&guide_curve);
            let loc = Rc::new(RefCell::new(LocationGuide::new(Box::new(t_law))));
            self.my_location = Some(Rc::new(RefCell::new(BRepFillACRLaw::new(
                brep,
                &self.my_spine,
                &loc,
            ))));
        } else {
            // trihedron by plane
            // OCCT L416-430.
            if keep_contact == BRepFillTypeOfContact::Contact
                || keep_contact == BRepFillTypeOfContact::ContactOnBorder
            {
                self.my_trihedron = GeomFillTrihedron::IsGuidePlanWithContact; // with rotation
            } else {
                self.my_trihedron = GeomFillTrihedron::IsGuidePlan; // without rotation
            }

            let t_law = GuideTrihedronPlan::new(&guide_curve);
            let loc = Rc::new(RefCell::new(LocationGuide::new(Box::new(t_law))));
            self.my_location = Some(Rc::new(RefCell::new(BRepFillEdge3DLaw::new(
                brep,
                &self.my_spine,
                loc,
            ))));
        }
        self.my_section = None; // It is required to relocalize the sections.
    }

    /// OCCT BRepFill_PipeShell::SetMaxDegree(NewMaxDegree) (cxx L436-439).
    pub fn set_max_degree(&mut self, new_max_degree: i32) {
        self.my_max_degree = new_max_degree;
    }

    /// OCCT BRepFill_PipeShell::SetMaxSegments(NewMaxSegments) (cxx L443-446).
    pub fn set_max_segments(&mut self, new_max_segments: i32) {
        self.my_max_segments = new_max_segments;
    }

    /// OCCT BRepFill_PipeShell::SetForceApproxC1(ForceApproxC1) (cxx L454-457).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_force_approx_c1 = force_approx_c1;
    }

    /// OCCT inline SetIsBuildHistory (hxx L125-128).
    pub fn set_is_build_history(&mut self, is_build_history: bool) {
        self.my_is_build_history = is_build_history;
    }

    /// OCCT inline IsBuildHistory (hxx L132).
    pub fn is_build_history(&self) -> bool {
        self.my_is_build_history
    }

    /// OCCT BRepFill_PipeShell::Add(Profile, WithContact, WithCorrection)
    /// (cxx L461-469).
    pub fn add(
        &mut self,
        brep: &mut BRep,
        profile: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        // OCCT L465-468.
        let v = Shape::null();
        self.add_with_location(brep, profile, &v, with_contact, with_correction);
        self.reset_loc();
    }

    /// OCCT BRepFill_PipeShell::Add(Profile, Location, WithContact,
    /// WithCorrection) (cxx L473-548).
    pub fn add_with_location(
        &mut self,
        brep: &mut BRep,
        profile: &Shape,
        location: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        // OCCT L478: DeleteProfile(Profile) — No duplication.
        self.delete_profile(profile);
        // OCCT L479-540.
        if self.my_is_automatic_law {
            self.my_seq.clear();
            let mut s = BRepFillSection::new(profile, location, with_contact, with_correction);
            s.set(true);
            self.my_seq.push(s);
            self.my_section = None;
            self.reset_loc();

            // OCCT L488-489.
            let law = self
                .my_location
                .as_ref()
                .expect("myLocation")
                .borrow()
                .base()
                .law(1);
            // OCCT occ::down_cast<GeomFill_LocationGuide>(myLocation->Law(1)).
            let law_ref = law.borrow();
            let loc = law_ref
                .as_any()
                .downcast_ref::<LocationGuide>()
                .expect("occ::down_cast<GeomFill_LocationGuide>");
            // OCCT L491.
            let (_status, mut par_and_rad) = loc.compute_automatic_law();

            // Compuite initial width of section (this will be 1.)
            // OCCT L494-496.
            let bary_center = brep_gprop_linear_properties_centre(profile);

            // OCCT L498-500.
            // only plane
            let profile_face = BRepLibMakeFaceWire::new(brep, profile, true).face();
            let the_plane = crate::brep_fill::offset_wire_b::brep_tool_surface(brep, &profile_face)
                .expect("BRep_Tool::Surface(ProfileFace)");

            // OCCT L501-505.
            let mut intersector = IntCurveSurfaceHInter::default();
            let a_hcurve: [Curve3; 2] = [
                loc.get_curve().expect("null GetCurve"),
                loc.guide().expect("null Guide"),
            ];
            let mut points_on_spines: [DVec3; 2] = [DVec3::ZERO, DVec3::ZERO];

            // OCCT L508-522.
            for (i, curve) in a_hcurve.iter().enumerate() {
                intersector.perform(curve, &the_plane);
                let mut min_dist = REAL_LAST;
                for j in 1..=intersector.nb_points() {
                    let a_pint = intersector.point(j).pnt();
                    let a_dist = bary_center.distance(a_pint);
                    if a_dist < min_dist {
                        min_dist = a_dist;
                        points_on_spines[i] = a_pint;
                    }
                }
            }

            // Correct <ParAndRad> according to <InitialWidth>
            // OCCT L525-532.
            let initial_width = points_on_spines[0].distance(points_on_spines[1]);
            let nb_par_rad = par_and_rad.len();
            for par_rad in par_and_rad.iter_mut() {
                par_rad.y /= initial_width;
            }

            // OCCT L534-539.
            let interpol = Rc::new(RefCell::new(LawInterpol::new()));
            self.my_law = Some(interpol.clone());

            let is_periodic = (par_and_rad[0].y - par_and_rad[nb_par_rad - 1].y).abs() < CONFUSION;

            // OCCT L538: (occ::down_cast<Law_Interpol>(myLaw))->Set(...).
            interpol.borrow_mut().set(&par_and_rad, is_periodic);
        } else {
            // OCCT L543-546.
            let s = BRepFillSection::new(profile, location, with_contact, with_correction);
            self.my_seq.push(s);
            self.my_section = None;
            self.reset_loc();
        }
    }

    /// OCCT BRepFill_PipeShell::SetLaw(Profile, L, WithContact, WithCorrection)
    /// (cxx L554-563).
    pub fn set_law(
        &mut self,
        profile: &Shape,
        l: &LawFunctionHandle,
        with_contact: bool,
        with_correction: bool,
    ) {
        // OCCT L559-562.
        let v = Shape::null();
        self.set_law_with_location(profile, l, &v, with_contact, with_correction);
        self.reset_loc();
    }

    /// OCCT BRepFill_PipeShell::SetLaw(Profile, L, Location, WithContact,
    /// WithCorrection) (cxx L569-582).
    pub fn set_law_with_location(
        &mut self,
        profile: &Shape,
        l: &LawFunctionHandle,
        location: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        // OCCT L575-581.
        self.my_seq.clear();
        let mut s = BRepFillSection::new(profile, location, with_contact, with_correction);
        s.set(true);
        self.my_seq.push(s);
        self.my_law = Some(Rc::clone(l));
        self.my_section = None;
        self.reset_loc();
    }

    /// OCCT BRepFill_PipeShell::DeleteProfile(Profile) (cxx L586-605).
    pub fn delete_profile(&mut self, profile: &Shape) {
        // OCCT L588-598.
        let mut trouve = false;
        let mut ii = 1usize;
        while ii <= self.my_seq.len() && !trouve {
            let a_section = self.my_seq[ii - 1].original_shape();
            if profile.ptr_id() == a_section.ptr_id() {
                trouve = true;
                self.my_seq.remove(ii - 1);
            } else {
                ii += 1;
            }
        }

        // OCCT L600-604.
        if trouve {
            self.my_section = None;
        }
        self.reset_loc();
    }

    /// OCCT BRepFill_PipeShell::IsReady() const (cxx L609-612).
    pub fn is_ready(&self) -> bool {
        !self.my_seq.is_empty()
    }

    /// OCCT BRepFill_PipeShell::GetStatus() const (cxx L616-619).
    pub fn get_status(&self) -> PipeError {
        self.my_status
    }

    /// OCCT BRepFill_PipeShell::SetTolerance(Tol3d, BoundTol, TolAngular)
    /// (cxx L623-630).
    pub fn set_tolerance(&mut self, tol3d: f64, bound_tol: f64, tol_angular: f64) {
        // OCCT L627-629.
        self.my_tol3d = tol3d;
        self.my_bound_tol = bound_tol;
        self.my_tol_angular = tol_angular;
    }

    /// OCCT BRepFill_PipeShell::SetTransition(Mode, Angmin, Angmax)
    /// (cxx L636-647).
    pub fn set_transition(&mut self, mode: BRepFillTransitionStyle, angmin: f64, angmax: f64) {
        // OCCT L640-643.
        if self.my_transition != mode {
            self.my_section = None; // It is required to relocalize the sections.
        }
        // OCCT L644-646.
        self.my_transition = mode;
        self.angmin = angmin;
        self.angmax = angmax;
    }

    /// OCCT BRepFill_PipeShell::Simulate(N, List) (cxx L651-701).
    pub fn simulate(&mut self, brep: &mut BRep, n: i32, list: &mut Vec<Shape>) {
        // Preparation
        // OCCT L654-655.
        self.prepare(brep);
        list.clear();

        // OCCT L657-660.
        let mut first = 0.0f64;
        let mut last = 0.0f64;
        let mut length = 0.0f64;
        let mut delta;
        let mut u;
        let mut us = 0.0f64;
        let mut delta_s;
        let mut first_s = 0.0f64;
        let nb_l = self
            .my_location
            .as_ref()
            .expect("myLocation")
            .borrow()
            .base()
            .nb_law();
        let mut finis = false;
        let mut w = Shape::null();
        let mut ii = 1usize;

        // Calculate the parameters of digitalization
        // OCCT L663-670.
        self.my_section
            .as_ref()
            .expect("mySection")
            .borrow()
            .base()
            .law(1)
            .borrow()
            .get_domain(&mut first_s, &mut last);
        delta_s = last - first_s;
        self.my_location
            .as_ref()
            .expect("myLocation")
            .borrow_mut()
            .base_mut()
            .curvilinear_bounds(nb_l, &mut first, &mut length);
        delta = length;
        if n > 1 {
            delta /= (n - 1) as f64;
        }

        // OCCT L672-700.
        self.my_location
            .as_ref()
            .expect("myLocation")
            .borrow_mut()
            .base_mut()
            .curvilinear_bounds(1, &mut first, &mut last); // Initiation of Last
        u = 0.0;
        while !finis {
            if u >= length {
                u = length;
                finis = true;
            } else {
                if ii < nb_l as usize {
                    self.my_location
                        .as_ref()
                        .expect("myLocation")
                        .borrow_mut()
                        .base_mut()
                        .curvilinear_bounds(nb_l, &mut first, &mut last);
                }
                if u > last {
                    u = (last + first) / 2.0; // The edge is not skipped
                }
                if u > first {
                    ii += 1;
                }
            }
            us = first_s + (u / length) * delta_s;
            // Calcul d'une section
            self.my_section
                .as_ref()
                .expect("mySection")
                .borrow_mut()
                .d0(brep, us, &mut w);
            self.my_location
                .as_ref()
                .expect("myLocation")
                .borrow_mut()
                .base_mut()
                .d0(u, &mut w);
            list.push(w.clone());
            u += delta;
        }
    }

    /// OCCT BRepFill_PipeShell::Build() (cxx L707-845).
    pub fn build(&mut self, brep: &mut BRep) -> bool {
        // 1) Preparation
        // OCCT L712.
        self.prepare(brep);

        // OCCT L714-721.
        if self.my_status != PIPE_OK {
            // BRep_Builder B; TopoDS_Shell Sh; B.MakeShell(Sh); myShape = Sh;
            let sh = brep.add_tshell(vec![]);
            self.my_shape = sh; // Nullify
            return false;
        }

        // 2) Calculate myFirst and myLast
        // OCCT L723-756.
        let mut first_s = 0.0f64;
        let mut last_s = 0.0f64;
        let section = self.my_section.as_ref().expect("mySection");
        section
            .borrow()
            .base()
            .law(1)
            .borrow()
            .get_domain(&mut first_s, &mut last_s);
        let mut my_first = Shape::null();
        section.borrow_mut().d0(brep, first_s, &mut my_first);
        self.my_location
            .as_ref()
            .expect("myLocation")
            .borrow_mut()
            .base_mut()
            .d0(0.0, &mut my_first);
        self.my_first = my_first;
        let loc_handle = self.my_location.as_ref().expect("myLocation");
        if section.borrow().base().is_vclosed() && loc_handle.borrow().base().is_closed(brep) {
            // OCCT L729: myLocation->IsG1(0) — the OCCT header defaults
            // (SpatialTolerance = 1.0e-7, AngularTolerance = 1.0e-4).
            if loc_handle.borrow().base().is_g1(brep, 0, 1.0e-7, 1.0e-4) >= 0 {
                self.my_last = self.my_first.clone();
            } else {
                self.my_first = Shape::null();
                self.my_last = Shape::null();
            }
        } else {
            let mut length = 0.0f64;
            let nb_law = loc_handle.borrow().base().nb_law();
            loc_handle
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(nb_law, &mut first_s, &mut length);
            let mut my_last = Shape::null();
            self.my_section
                .as_ref()
                .expect("mySection")
                .borrow_mut()
                .d0(brep, last_s, &mut my_last);
            loc_handle.borrow_mut().base_mut().d0(length, &mut my_last);
            self.my_last = my_last;
            // eap 5 Jun 2002 occ332, myLast and myFirst must not share one TShape,
            // tolerances of shapes built on them may be quite different
            let partner = Arc::ptr_eq(&self.my_first.data, &self.my_last.data);
            if partner {
                // OCCT L749-753.
                let copy = BRepBuilderAPICopy::new(&self.my_last);
                if copy.is_done() {
                    self.my_last = copy.shape();
                }
            }
            // eap 5 Jun 2002 occ332, end modif
        }

        // 3) Construction
        // OCCT L758-786.
        let mut mk_sw = BRepFillSweep::new(
            brep,
            self.my_section.as_ref().expect("mySection").clone(),
            self.my_location.as_ref().expect("myLocation").clone(),
            true,
        );
        mk_sw.set_tolerance(brep, self.my_tol3d, self.my_bound_tol, 1.0e-5, self.my_tol_angular);
        mk_sw.set_angular_control(self.angmin, self.angmax);
        mk_sw.set_force_approx_c1(self.my_force_approx_c1);
        let first_wire = self.my_first.clone();
        let last_wire = self.my_last.clone();
        mk_sw.set_bounds(brep, &first_wire, &last_wire);

        let the_continuity = if self.my_trihedron == GeomFillTrihedron::IsDiscreteTrihedron {
            GeomAbsShape::C0
        } else {
            GeomAbsShape::C2
        };
        let mut dummy: ShapeSet = ShapeSet::new();
        let mut dummy2: ShapeToArray2Map = ShapeToArray2Map::new();
        let mut dummy3: ShapeToArray2Map = ShapeToArray2Map::new();
        mk_sw.build(
            brep,
            &mut dummy,
            &mut dummy2,
            &mut dummy3,
            self.my_transition,
            the_continuity,
            crate::geomalgo::geomfill::sweep::GeomFillApproxStyle::GeomFill_Location,
            self.my_max_degree,
            self.my_max_segments,
        );

        // OCCT L788-789.
        self.my_status = self
            .my_location
            .as_ref()
            .expect("myLocation")
            .borrow()
            .base()
            .get_status();
        let ok = mk_sw.is_done() && self.my_status == PIPE_OK;

        // OCCT L791-843.
        if ok {
            self.my_shape = mk_sw.shape();
            self.my_error_on_surf = mk_sw.error_on_surface();

            // OCCT L796-803.
            let mut a_bottom_wire = self.my_first.clone();
            let mut a_top_wire = self.my_last.clone();

            if build_boundaries(
                brep,
                &mk_sw,
                self.my_section.as_ref().expect("mySection"),
                &mut a_bottom_wire,
                &mut a_top_wire,
            ) {
                self.my_first = a_bottom_wire;
                self.my_last = a_top_wire;
            }

            // OCCT L805-826.
            let section = self.my_section.as_ref().expect("mySection");
            if section.borrow().base().is_uclosed() {
                let mut degen_first = true;
                let mut degen_last = true;

                for e in explored_children(brep, &self.my_first, ShapeType::Edge) {
                    degen_first = degen_first && brep.edge(e).degenerated;
                }

                for e in explored_children(brep, &self.my_last, ShapeType::Edge) {
                    degen_last = degen_last && brep.edge(e).degenerated;
                }

                if degen_first && degen_last {
                    set_shape_closed(brep, &self.my_shape, true);
                }
            }

            // OCCT L828-831.
            if self.my_is_build_history {
                let sweep = &mk_sw;
                self.build_history(brep, sweep);
            }
        } else {
            // OCCT L833-843.
            let sh = brep.add_tshell(vec![]);
            self.my_shape = sh; // Nullify
            if self.my_status == PIPE_OK {
                self.my_status = PIPE_NOT_OK;
            }
        }
        // OCCT L844.
        ok
    }

    /// OCCT BRepFill_PipeShell::MakeSolid() (cxx L849-914).
    pub fn make_solid(&mut self, brep: &mut BRep) -> bool {
        // OCCT L851-854.
        if self.my_shape.is_null() {
            panic!("PipeShell is not built");
        }
        // OCCT L855.
        let mut b = is_closed_wire(brep, &self.my_shape);

        // OCCT L857-895.
        if !b {
            if !self.my_first.is_null() && !self.my_last.is_null() {
                b = is_closed_wire(brep, &self.my_first) && is_closed_wire(brep, &self.my_last);
            }
            if b {
                // It is necessary to block the extremities
                // OCCT L867-870.
                b = perform_plan(brep, &mut self.my_first);
                if b {
                    b = perform_plan(brep, &mut self.my_last);
                    if b {
                        // OCCT L873-880.
                        if !self.my_first.is_null() && !is_same_oriented(brep, &self.my_first, &self.my_shape)
                        {
                            self.my_first = shape_reversed(&self.my_first);
                        }
                        if !self.my_last.is_null() && !is_same_oriented(brep, &self.my_last, &self.my_shape)
                        {
                            self.my_last = shape_reversed(&self.my_last);
                        }

                        // OCCT L882-889: BS.Add(myShape, TopoDS::Face(myFirst /
                        // myLast)) — the cap faces are appended to the sweep
                        // shell.
                        if !self.my_first.is_null() {
                            let f = self.my_first.clone();
                            add_shell_face(brep, &self.my_shape, &f);
                        }
                        if !self.my_last.is_null() {
                            let f = self.my_last.clone();
                            add_shell_face(brep, &self.my_shape, &f);
                        }

                        // OCCT L891.
                        set_shape_closed(brep, &self.my_shape, true);
                    }
                }
            }
        }

        // OCCT L897-913.
        if b {
            // OCCT L899-901: TopoDS_Solid solid; BS.MakeSolid(solid);
            // BS.Add(solid, TopoDS::Shell(myShape)).
            let solid = brep.add_tsolid(vec![self.my_shape.clone()]);
            // OCCT L902-904.
            let mut sc = crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(&solid);
            sc.perform_infinite_point(CONFUSION);
            if sc.state() == 0 {
                // TopAbs_IN
                // OCCT L906-908.
                self.my_shape = shape_reversed(&self.my_shape);
                let _solid = brep.add_tsolid(vec![self.my_shape.clone()]);
            }
            // OCCT L910-911.
            self.my_shape = solid;
            set_shape_closed(brep, &self.my_shape, true);
        }
        // OCCT L913.
        b
    }

    /// OCCT BRepFill_PipeShell::Shape() const (cxx L918-921).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepFill_PipeShell::ErrorOnSurface() const (cxx L925-928).
    pub fn error_on_surface(&self) -> f64 {
        self.my_error_on_surf
    }

    /// OCCT BRepFill_PipeShell::FirstShape() const (cxx L932-935).
    pub fn first_shape(&self) -> Shape {
        self.my_first.clone()
    }

    /// OCCT BRepFill_PipeShell::LastShape() const (cxx L939-942).
    pub fn last_shape(&self) -> Shape {
        self.my_last.clone()
    }

    /// OCCT inline Profiles(theProfiles) (hxx L204-208).
    pub fn profiles(&self) -> Vec<Shape> {
        let mut the_profiles = Vec::new();
        for i in 1..=self.my_seq.len() {
            the_profiles.push(self.my_seq[i - 1].original_shape());
        }
        the_profiles
    }

    /// OCCT inline Spine() (hxx L211).
    pub fn spine(&self) -> Shape {
        self.my_spine.clone()
    }

    /// OCCT BRepFill_PipeShell::Generated(theShape, theList) (cxx L946-955).
    pub fn generated(&self, the_shape: &Shape) -> Vec<Shape> {
        // OCCT L949-954.
        match self.my_gen_map.get(&shape_key(the_shape)) {
            Some(l) => l.clone(),
            None => Vec::new(),
        }
    }

    /// OCCT BRepFill_PipeShell::Prepare() (cxx L964-1234).
    fn prepare(&mut self, brep: &mut BRep) {
        // OCCT L966-967.
        self.w_seq.clear();
        self.my_edge_new_edges.clear();

        // OCCT L969-973.
        if !self.is_ready() {
            panic!("PipeShell");
        }
        // OCCT L974-977.
        if self.my_location.is_some() && self.my_section.is_some() {
            return; // It is ready
        }

        // Check set of section for right configuration of punctual sections
        // OCCT L979-1010.
        let mut i = 2i32;
        while i <= self.my_seq.len() as i32 - 1 {
            let mut wdeg = true;
            for an_edge in wire_edges(brep, &self.my_seq[(i - 1) as usize].wire()) {
                wdeg = wdeg && brep.edge(an_edge).degenerated;
            }
            if wdeg {
                panic!("Wrong usage of punctual sections");
            }
            i += 1;
        }
        if self.my_seq.len() <= 2 {
            let mut wdeg = true;
            for i in 1..=self.my_seq.len() {
                for an_edge in wire_edges(brep, &self.my_seq[i - 1].wire()) {
                    wdeg = wdeg && brep.edge(an_edge).degenerated;
                }
            }
            if wdeg {
                panic!("Wrong usage of punctual sections");
            }
        }

        // Construction of the law of location
        // OCCT L1012-1027.
        if self.my_location.is_none() {
            match self.my_trihedron {
                GeomFillTrihedron::IsCorrectedFrenet => {
                    let t_law: Box<dyn TrihedronLaw> = Box::new(CorrectedFrenet::new());
                    let loc = Rc::new(RefCell::new(CurveAndTrihedron::new(t_law)));
                    self.my_location = Some(Rc::new(RefCell::new(BRepFillEdge3DLaw::new(
                        brep,
                        &self.my_spine,
                        loc,
                    ))));
                }
                _ => {
                    // Not planned!
                    panic!("PipeShell");
                }
            }
        }

        // Transformation of the law (Transition Management)
        // OCCT L1030.
        perform_transition(brep, self.my_transition, &mut self.my_location, self.angmin);

        // Construction of the section law
        // OCCT L1032-1203.
        if self.my_seq.len() == 1 {
            // OCCT L1035-1037.
            let mut the_sect = Shape::null();
            let mut a_trsf = DAffine3::IDENTITY;
            let mut p1 = 0.0f64;
            let sec = self.my_seq[0].clone();
            self.place(brep, &sec, &mut the_sect, &mut a_trsf, &mut p1);
            let a_local_shape = the_sect.clone();
            // OCCT L1039-1047.
            if sec.is_law() {
                self.my_section = Some(Rc::new(RefCell::new(BRepFillShapeLaw::new_with_law(
                    brep,
                    &a_local_shape,
                    self.my_law.as_ref().expect("myLaw").clone(),
                    true,
                ))));
            } else {
                self.my_section = Some(Rc::new(RefCell::new(BRepFillShapeLaw::new(
                    brep,
                    &a_local_shape,
                    true,
                ))));
            }

            // OCCT L1050-1060.
            self.w_seq.push(the_sect.clone());
            // Simple case of single section
            self.my_ind_of_sec.push(1);
            for an_edge in wire_edges(brep, &the_sect) {
                self.my_edge_new_edges.insert(shape_key(&an_edge), vec![an_edge]);
            }
            ///////////////////////////////
        } else {
            // OCCT L1065-1202.
            let mut param: Vec<f64> = Vec::new();
            let mut ind_sec: Vec<i32> = Vec::new();
            let mut transformations: Vec<Trsf> = Vec::new();
            let nb_l = self
                .my_location
                .as_ref()
                .expect("myLocation")
                .borrow()
                .base()
                .nb_law();
            let mut v1 = 0.0f64;
            let mut v2 = 0.0f64;
            self.my_location
                .as_ref()
                .expect("myLocation")
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(nb_l, &mut v1, &mut v2);
            v1 = 0.0;
            let mut ideb = 0usize;
            let mut ifin = 0usize;
            for iseq in 1..=self.my_seq.len() {
                ind_sec.push(iseq as i32);
                let mut the_sect = Shape::null();
                let mut a_trsf = DAffine3::IDENTITY;
                let mut cur_param = 0.0f64;
                let sec = self.my_seq[iseq - 1].clone();
                self.place(brep, &sec, &mut the_sect, &mut a_trsf, &mut cur_param);
                param.push(cur_param);
                self.w_seq.push(the_sect);
                transformations.push(trsf_from_daffine3(&a_trsf));
                if cur_param == v1 {
                    ideb = iseq;
                }
                if cur_param == v2 {
                    ifin = iseq;
                }
            }

            // looping sections ?
            // OCCT L1093-1140.
            if self
                .my_location
                .as_ref()
                .expect("myLocation")
                .borrow()
                .base()
                .is_closed(brep)
            {
                if ideb > 0 {
                    // place the initial section at the final position
                    param.push(v2);
                    let w = self.w_seq[ideb - 1].clone();
                    self.w_seq.push(w);
                } else if ifin > 0 {
                    // place the final section at the initial position
                    param.push(v1);
                    let w = self.w_seq[ifin - 1].clone();
                    self.w_seq.push(w);
                } else {
                    // it is necessary to find a medium section to impose by V1 and by V2
                    let mut pmin = REAL_LAST;
                    let mut pmax = f64::MIN; // OCCT RealFirst()
                    let mut wmin = Shape::null();
                    let mut wmax = Shape::null();
                    for iseq in 1..=self.w_seq.len() {
                        if param[iseq - 1] < pmin {
                            pmin = param[iseq - 1];
                            wmin = self.w_seq[iseq - 1].clone();
                        }
                        if param[iseq - 1] > pmax {
                            pmax = param[iseq - 1];
                            wmax = self.w_seq[iseq - 1].clone();
                        }
                    }
                    // medium section between Wmin and Wmax
                    let mut wres = Shape::null();
                    let dmin = (pmin - v1).abs();
                    let dmax = (pmax - v2).abs();
                    if compute_section(brep, &wmin, &wmax, dmin, dmax, &mut wres) {
                        // impose section Wres at the beginning and the end
                        param.push(v1);
                        self.w_seq.push(wres.clone());
                        ind_sec.push(self.w_seq.len() as i32);
                        param.push(v2);
                        self.w_seq.push(wres);
                        ind_sec.push(self.w_seq.len() as i32);
                    }
                }
            }

            // parse sections by increasing parameter
            // OCCT L1142-1160.
            let mut play_again = true;
            while play_again {
                play_again = false;
                for iseq in 1..=self.w_seq.len() {
                    for jseq in (iseq + 1)..=self.w_seq.len() {
                        if param[iseq - 1] > param[jseq - 1] {
                            param.swap(iseq - 1, jseq - 1);
                            self.w_seq.swap(iseq - 1, jseq - 1);
                            ind_sec.swap(iseq - 1, jseq - 1);
                            play_again = true;
                        }
                    }
                }
            }
            // Fill the array of real indices of sections
            // OCCT L1161-1172.
            for ii in 1..=self.my_seq.len() {
                for jj in 1..=ind_sec.len() {
                    if ind_sec[jj - 1] == ii as i32 {
                        self.my_ind_of_sec.push(jj as i32);
                        break;
                    }
                }
            }

            //  Calculate work sections
            // OCCT L1174-1202.
            let mut working_sections: Vec<Shape> = Vec::new();
            working_sections.clear();
            let mut working_map: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();
            let mut georges = BRepFillCompatibleWires::new_with_sections(brep, &self.w_seq);
            georges.set_percent(0.1);
            georges.perform(brep, false);
            if georges.is_done() {
                working_sections = georges.shape().clone();
                for (k, v) in georges.generated().iter() {
                    working_map.insert(*k, v.clone());
                }
                // For each sub-edge of each section
                // we save its splits
                for ii in 1..=self.w_seq.len() {
                    for e in explored_children(brep, &self.w_seq[ii - 1], ShapeType::Edge) {
                        let a_new_edges = georges.generated_shapes(&e);
                        self.my_edge_new_edges.insert(shape_key(&e), a_new_edges);
                    }
                }
            } else {
                panic!("PipeShell : uncompatible wires");
            }
            self.my_section = Some(Rc::new(RefCell::new(BRepFillNSections::new_with_params(
                brep,
                working_sections,
                transformations,
                param,
                v1,
                v2,
                true,
            ))));
        } // else

        //  modify the law of location if contact
        // OCCT L1207-1227.
        if self.my_trihedron == GeomFillTrihedron::IsGuidePlanWithContact
            || self.my_trihedron == GeomFillTrihedron::IsGuideACWithContact
        {
            let mut fs = 0.0f64;
            let mut f = 0.0f64;
            let mut l = 0.0f64;
            let mut delta = 0.0f64;
            let mut length = 0.0f64;
            let nb_law = self
                .my_location
                .as_ref()
                .expect("myLocation")
                .borrow()
                .base()
                .nb_law();
            let sec_handle = self
                .my_section
                .as_ref()
                .expect("mySection")
                .borrow()
                .concatened_law(brep)
                .expect("null ConcatenedLaw");
            self.my_location
                .as_ref()
                .expect("myLocation")
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(nb_law, &mut f, &mut length);
            sec_handle.borrow().get_domain(&mut fs, &mut l);
            delta = (l - fs) / length;

            let mut old_angle = 0.0f64;
            for ipath in 1..=nb_law {
                self.my_location
                    .as_ref()
                    .expect("myLocation")
                    .borrow_mut()
                    .base_mut()
                    .curvilinear_bounds(ipath, &mut f, &mut l);
                let law = self
                    .my_location
                    .as_ref()
                    .expect("myLocation")
                    .borrow()
                    .base()
                    .law(ipath);
                // OCCT occ::down_cast<GeomFill_LocationGuide>(myLocation->Law(ipath)).
                let mut law_mut = law.borrow_mut();
                let loc = law_mut
                    .as_any_mut()
                    .downcast_mut::<LocationGuide>()
                    .expect("occ::down_cast<GeomFill_LocationGuide>");
                let mut angle = 0.0f64;
                // force the rotation
                loc.set(
                    Rc::clone(&sec_handle),
                    true,
                    fs + f * delta,
                    fs + l * delta,
                    old_angle,
                    &mut angle,
                );
                old_angle = angle;
            }
        }

        // OCCT L1229-1233.
        self.my_status = self
            .my_location
            .as_ref()
            .expect("myLocation")
            .borrow()
            .base()
            .get_status();
        if !self
            .my_section
            .as_ref()
            .expect("mySection")
            .borrow()
            .base()
            .is_done()
        {
            self.my_status = PIPE_NOT_OK;
        }
    }

    /// OCCT BRepFill_PipeShell::Place(Sec, W, aTrsf, param) (cxx L1241-1257).
    fn place(
        &mut self,
        brep: &mut BRep,
        sec: &BRepFillSection,
        w: &mut Shape,
        a_trsf: &mut DAffine3,
        param: &mut f64,
    ) {
        // OCCT L1246-1250.
        let mut place = BRepFillSectionPlacement::new_with_vertex(
            brep,
            Rc::clone(self.my_location.as_ref().expect("myLocation")),
            &sec.wire(),
            &sec.vertex(),
            sec.with_contact(),
            sec.with_correction(),
        );
        // OCCT L1252.
        *a_trsf = place.transformation();
        // Transform the copy
        // OCCT L1254.
        let tmp_wire = sec.wire();
        let transformed = BRepBuilderAPITransform::new(&tmp_wire, a_trsf, true);
        *w = transformed.shape();
        ////////////////////////////////////
        // OCCT L1256.
        *param = place.abscissa_on_path();
    }

    /// OCCT BRepFill_PipeShell::ResetLoc() (cxx L1263-1275).
    fn reset_loc(&mut self) {
        // OCCT L1265-1274.
        if self.my_trihedron == GeomFillTrihedron::IsGuidePlanWithContact
            || self.my_trihedron == GeomFillTrihedron::IsGuideACWithContact
        {
            let nb_law = self
                .my_location
                .as_ref()
                .expect("myLocation")
                .borrow()
                .base()
                .nb_law();
            for isec in 1..=nb_law {
                let law = self
                    .my_location
                    .as_ref()
                    .expect("myLocation")
                    .borrow()
                    .base()
                    .law(isec);
                // OCCT occ::down_cast<GeomFill_LocationGuide>(myLocation->Law(isec)).
                let mut law_mut = law.borrow_mut();
                let loc = law_mut
                    .as_any_mut()
                    .downcast_mut::<LocationGuide>()
                    .expect("occ::down_cast<GeomFill_LocationGuide>");
                loc.erase_rotation(); // remove the rotation
            }
        }
    }

    /// OCCT BRepFill_PipeShell::BuildHistory(theSweep) (cxx L1282-1615).
    fn build_history(&mut self, brep: &mut BRep, the_sweep: &BRepFillSweep) {
        // Filling of <myGenMap>
        // OCCT L1285-1288.
        let an_u_edges = the_sweep.inter_faces().expect("myUEdges");

        // OCCT L1288: IndWireMap.
        let mut ind_wire_map: HashMap<i32, Shape> = HashMap::new();

        // OCCT L1290-1313.
        let mut inde;
        for indw in 1..=self.my_seq.len() {
            let section = self.my_seq[indw - 1].original_shape();
            let mut a_section;
            let is_punctual = self.my_seq[indw - 1].is_punctual();
            if is_punctual {
                // for punctual sections (first or last)
                // we take all the wires generated along the path
                // OCCT L1302-1311.
                let e_list = self
                    .my_gen_map
                    .entry(shape_key(&section))
                    .or_default();
                for i in 1..=an_u_edges.len() {
                    for j in 1..=an_u_edges[0].len() {
                        e_list.push(an_u_edges[i - 1][j - 1].clone());
                    }
                }

                continue;
            } else {
                a_section = section.clone();
            }
            // Take the real index of section on the path
            // OCCT L1318-1334.
            let ind_of_w = self.my_ind_of_sec[indw - 1] as usize;
            let the_wire = self.w_seq[ind_of_w - 1].clone();
            let wexp_sec = wire_edges(brep, &a_section);
            inde = 1usize;
            for sec_edge in &wexp_sec {
                let an_original_edge = sec_edge.clone();
                let an_edge = self.my_seq[indw - 1].modified_shape(&an_original_edge);
                if brep.edge(an_edge.clone()).degenerated {
                    continue;
                }

                // OCCT L1331-1332.
                let a_shell = brep.add_tshell(vec![]);
                // OCCT L1333-1335: TopoDS_Vertex aVertex[2].
                let a_vertex = {
                    let (vf, vl) = top_exp_vertices(brep, &an_original_edge);
                    [vf, vl]
                };
                let sign_of_an_edge = if an_original_edge.orientation
                    == Orientation::Forward
                {
                    1
                } else {
                    -1
                };

                // For each non-degenerated inde-th edge of <aSection>
                // we find inde-th edge in <theWire>
                // OCCT L1339-1353.
                let mut the_edge = Shape::null();
                let wexp = wire_edges(brep, &the_wire);
                let mut i = 1usize;
                for we_edge in &wexp {
                    the_edge = we_edge.clone();
                    if brep.edge(an_edge.clone()).degenerated {
                        continue;
                    }
                    if i == inde {
                        break;
                    }
                    i += 1;
                }

                // Take the list of splits for <theEdge>
                // OCCT L1356-1375.
                let new_edges = self
                    .my_edge_new_edges
                    .get(&shape_key(&the_edge))
                    .cloned()
                    .unwrap_or_default();
                let mut sign_of_a_new_edge = 0i32;
                let mut sign_of_index = 0i32;
                for iter_value in &new_edges {
                    let a_new_edge = iter_value.clone();
                    sign_of_a_new_edge = if a_new_edge.orientation
                        == Orientation::Forward
                    {
                        1
                    } else {
                        -1
                    };
                    let mut an_ind_e = self
                        .my_section
                        .as_ref()
                        .expect("mySection")
                        .borrow()
                        .base()
                        .index_of_edge(&a_new_edge);
                    sign_of_index = if an_ind_e > 0 { 1 } else { -1 };
                    an_ind_e = an_ind_e.abs();
                    // For an edge generated shape is a "tape" -
                    // a shell usually containing this edge and
                    // passing from beginning of path to its end
                    // OCCT L1369-1374.
                    let a_tape = the_sweep.tape(an_ind_e);
                    for child in direct_children(brep, &a_tape) {
                        add_shell_face(brep, &a_shell, &child);
                    }
                }

                // Processing of vertices of <anEdge>
                // We should choose right index in <anUEdges>
                // for each vertex of edge
                // OCCT L1380-1406.
                let to_reverse = sign_of_an_edge * sign_of_a_new_edge * sign_of_index;
                let mut u_index: [usize; 2] = [0, 0];
                u_index[0] = (self
                    .my_section
                    .as_ref()
                    .expect("mySection")
                    .borrow()
                    .base()
                    .index_of_edge(new_edges.first().expect("NewEdges.First"))
                    .abs()) as usize;
                u_index[1] = (self
                    .my_section
                    .as_ref()
                    .expect("mySection")
                    .borrow()
                    .base()
                    .index_of_edge(new_edges.last().expect("NewEdges.Last"))
                    .abs()
                    + to_reverse) as usize;
                if to_reverse == -1 {
                    u_index[0] += 1;
                    u_index[1] += 1;
                }
                if self
                    .my_section
                    .as_ref()
                    .expect("mySection")
                    .borrow()
                    .base()
                    .is_uclosed()
                {
                    let nb_law = self
                        .my_section
                        .as_ref()
                        .expect("mySection")
                        .borrow()
                        .base()
                        .nb_law();
                    if u_index[0] > nb_law as usize {
                        u_index[0] = 1;
                    }
                    if u_index[1] > nb_law as usize {
                        u_index[1] = 1;
                    }
                }
                // if (SignOfAnEdge * SignOfANewEdge == -1)
                if sign_of_an_edge == -1 || sign_of_a_new_edge == -1 {
                    u_index.swap(0, 1);
                }

                // OCCT L1408-1412.
                let ve_map = map_shapes_and_ancestors_ve(brep, &a_shell);
                for kk in 0..2 {
                    if self.my_gen_map.contains_key(&shape_key(&a_vertex[kk])) {
                        continue;
                    }
                    if let Some(wire_bound) = ind_wire_map.get(&(u_index[kk] as i32)) {
                        let e_list = self
                            .my_gen_map
                            .entry(shape_key(&a_vertex[kk]))
                            .or_default();
                        for itw in wire_edges(brep, wire_bound) {
                            e_list.push(itw);
                        }

                        continue;
                    }

                    // Collect u-edges
                    // OCCT L1433-1438.
                    let mut seq_edges: Vec<Shape> = Vec::new();
                    for jj in 1..=an_u_edges[0].len() {
                        seq_edges.push(an_u_edges[u_index[kk] - 1][jj - 1].clone());
                    }

                    // Assemble the wire ("rail" along the path)
                    // checking for possible holes
                    //(they appear with option "Round Corner")
                    // and filling them
                    // Missed edges are taken from <aShell>
                    // OCCT L1445-1487.
                    let a_wire = brep.add_twire(vec![]);
                    let first_edge = seq_edges.first().cloned().expect("SeqEdges(1)");
                    if first_edge.is_null() {
                        continue;
                    }
                    append_edge_to_wire(brep, &a_wire, &first_edge);
                    let (first_vertex, mut cur_vertex) = top_exp_vertices(brep, &first_edge);
                    let mut cur_edge = Shape::null();
                    for jj in 2..=seq_edges.len() {
                        cur_edge = seq_edges[jj - 1].clone();
                        let (vfirst, vlast) = top_exp_vertices(brep, &cur_edge);
                        if cur_vertex.ptr_id() == vfirst.ptr_id() {
                            cur_vertex = vlast.clone();
                        } else {
                            // a hole
                            // OCCT L1467-1484.
                            let e_list = find_in_ve_map(&ve_map, &vfirst);
                            for itl in e_list {
                                let candidate = itl.clone();
                                if candidate.ptr_id() == cur_edge.ptr_id() {
                                    continue;
                                }
                                let (v1, v2) = top_exp_vertices(brep, &candidate);
                                if v1.ptr_id() == cur_vertex.ptr_id()
                                    || v2.ptr_id() == cur_vertex.ptr_id()
                                {
                                    append_edge_to_wire(brep, &a_wire, &candidate);
                                    break;
                                }
                            }
                        }
                        cur_vertex = vlast;
                        append_edge_to_wire(brep, &a_wire, &cur_edge);
                    } // for (jj = 2; jj <= SeqEdges.Length(); jj++)
                    // case of closed wire
                    // OCCT L1489-1508.
                    if self
                        .my_location
                        .as_ref()
                        .expect("myLocation")
                        .borrow()
                        .base()
                        .is_closed(brep)
                        && cur_vertex.ptr_id() != first_vertex.ptr_id()
                    {
                        let e_list = find_in_ve_map(&ve_map, &cur_vertex);
                        for itl in e_list {
                            let candidate = itl.clone();
                            if candidate.ptr_id() == cur_edge.ptr_id() {
                                continue;
                            }
                            let (v1, v2) = top_exp_vertices(brep, &candidate);
                            if v1.ptr_id() == first_vertex.ptr_id()
                                || v2.ptr_id() == first_vertex.ptr_id()
                            {
                                append_edge_to_wire(brep, &a_wire, &candidate);
                                break;
                            }
                        }
                    }

                    // OCCT L1510-1516.
                    let e_list = self
                        .my_gen_map
                        .entry(shape_key(&a_vertex[kk]))
                        .or_default();
                    for itw in wire_edges(brep, &a_wire) {
                        e_list.push(itw);
                    }

                    // Save already built wire with its index
                    // OCCT L1519.
                    ind_wire_map.insert(u_index[kk] as i32, a_wire);
                } // for (int kk = 0; kk < 2; kk++)
                ////////////////////////////////////

                // OCCT L1523-1529.
                let f_list = self
                    .my_gen_map
                    .entry(shape_key(&an_original_edge))
                    .or_default();
                for child in direct_children(brep, &a_shell) {
                    f_list.push(child);
                }
                ////////////////////////

                inde += 1;
            }
        }

        // For subshapes of spine
        // OCCT L1536-1614.
        let a_faces = the_sweep.sub_shape().expect("myFaces");
        let a_v_edges = the_sweep.sections().expect("myVEdges");

        let wexp_spine = wire_edges(brep, &self.my_spine);
        let mut inde = 0usize;
        let mut wexp_index = 0usize;
        loop {
            let mut to_exit = false;
            if wexp_index >= wexp_spine.len() {
                to_exit = true;
            }

            inde += 1;

            if !to_exit {
                let an_edge_of_spine = wexp_spine[wexp_index].clone();

                let f_list = self
                    .my_gen_map
                    .entry(shape_key(&an_edge_of_spine))
                    .or_default();

                for i in 1..=a_faces.len() {
                    let a_face = a_faces[i - 1][inde - 1].clone();
                    if a_face.shape_type() == ShapeType::Face {
                        f_list.push(a_face);
                    }
                }
            }

            // OCCT L1569: the traversal-first vertex of the current spine edge.
            let a_vertex_of_spine = if wexp_index < wexp_spine.len() {
                let (vf, _) = top_exp_vertices(brep, &wexp_spine[wexp_index]);
                vf
            } else {
                // OCCT reads CurrentVertex past the end only after ToExit; the
                // loop breaks below before using the list.
                Shape::null()
            };
            let list_vshapes = self
                .my_gen_map
                .entry(shape_key(&a_vertex_of_spine))
                .or_default();
            for i in 1..=a_v_edges.len() {
                let a_vshape = a_v_edges[i - 1][inde - 1].clone();
                if a_vshape.is_null() {
                    continue;
                }
                if a_vshape.shape_type() == ShapeType::Edge
                    || a_vshape.shape_type() == ShapeType::Face
                {
                    list_vshapes.push(a_vshape);
                } else {
                    for a_subshape in direct_children(brep, &a_vshape) {
                        if a_subshape.shape_type() == ShapeType::Edge
                            || a_subshape.shape_type() == ShapeType::Face
                        {
                            list_vshapes.push(a_subshape);
                        } else {
                            // it is wire
                            for itw in direct_children(brep, &a_subshape) {
                                list_vshapes.push(itw);
                            }
                        }
                    }
                }
            }

            if to_exit {
                break;
            }

            if wexp_index < wexp_spine.len() {
                wexp_index += 1;
            }
        }
    }
}

/// OCCT BRep_Builder::Add(shell, face) — the builder call used in
/// BuildHistory (kept as a tiny form helper over the shell child list).
fn add_shell_face(brep: &mut BRep, shell: &Shape, child: &Shape) {
    let sd = brep.shell_mut(shell.clone());
    sd.faces.push(child.clone());
    sd.my_shapes.push(child.clone());
}

/// OCCT TopoDS_Shape::Closed(true) on any shape kind (Build marks the result
/// shell closed; MakeSolid marks the solid closed).
fn set_shape_closed(brep: &mut BRep, s: &Shape, closed: bool) {
    let flag = rcad_kernel::topo::topods::tshape_flags::CLOSED;
    match s.data.as_ref() {
        TShape::Wire(_) => {
            set_wire_closed(brep, s, closed);
        }
        TShape::Shell(_) => {
            if closed {
                brep.shell_mut(s.clone()).flags |= flag;
            } else {
                brep.shell_mut(s.clone()).flags &= !flag;
            }
        }
        TShape::Solid(_) => {
            if closed {
                brep.solid_mut(s.clone()).flags |= flag;
            } else {
                brep.solid_mut(s.clone()).flags &= !flag;
            }
        }
        _ => {}
    }
}

/// OCCT NCollection_IndexedDataMap::FindFromKey for the VE ancestor map.
fn find_in_ve_map<'a>(map: &'a [(Shape, Vec<Shape>)], v: &Shape) -> &'a [Shape] {
    for (k, list) in map {
        if k.ptr_id() == v.ptr_id() {
            return list;
        }
    }
    &[]
}

// ===========================================================================
// Unit tests (write-only per the batch discipline: `cargo test` is not run;
// compile-checked via `cargo check --tests`)
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// A two-vertex spine wire (BRepLib_MakeWire form).
    fn spine_wire(brep: &mut BRep) -> Shape {
        let v1 = brep.add_tvertex_unique(DVec3::new(0.0, 0.0, 0.0));
        let v2 = brep.add_tvertex_unique(DVec3::new(1.0, 0.0, 0.0));
        let dir = DVec3::X;
        let e = brep.add_tedge(
            Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(
                DVec3::ZERO,
                dir,
            ))),
            v1,
            v2,
            [0.0, 1.0],
        );
        brep.add_twire(vec![e])
    }

    /// OCCT constructor defaults (cxx L224-251): MaxDegree = 11,
    /// MaxSegments = 100, CorrectedFrenet trihedron, Modified transition,
    /// PipeOk status, history on, no sections.
    #[test]
    fn pipe_shell_constructor_defaults() {
        let mut brep = BRep::new();
        let spine = spine_wire(&mut brep);
        let ps = BRepFillPipeShell::new(&mut brep, &spine);
        assert!(!ps.is_ready());
        assert_eq!(ps.get_status(), PipeError::PipeOk);
        assert!(ps.is_build_history());
        assert_eq!(ps.error_on_surface(), 0.0);
        assert!(ps.shape().is_null());
        assert!(ps.profiles().is_empty());
    }
}
