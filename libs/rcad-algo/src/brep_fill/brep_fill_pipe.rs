//! OCCT BRepFill_Pipe — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_Pipe.hxx (L17-150) +
//!         BRepFill_Pipe.cxx (L1-1221).
//!
//! Architecture differences (referenced from the affected functions):
//! 1. The sweep engine stack (BRepFill_Sweep, BRepFill_ShapeLaw,
//!    BRepFill_Edge3DLaw / BRepFill_LocationLaw,
//!    BRepFill_SectionPlacement, GeomFill_CurveAndTrihedron,
//!    GeomFill_DiscreteTrihedron, ShapeUpgrade_RemoveLocations,
//!    BRepBuilderAPI_Transform) is not translated yet — the OCCT-class-named
//!    GAP carriers below keep the call form (plan D3; the GeomFill trihedron
//!    laws themselves are the landed/parallel geomalgo::geomfill batch).
//! 2. NCollection_HArray2<TopoDS_Shape> -> the local [`ShapeArray2`]
//!    (1-based row/col indexing as in OCCT).
//! 3. NCollection_Map / NCollection_DataMap keyed by TopTools_ShapeMapHasher
//!    -> HashMap keyed by (TShape ptr, location) identity
//!    ([`ShapeKey`], orientation ignored).
//! 4. TopoDS_Iterator -> [`topods_iterator`] (the direct children in
//!    stored order); TopExp_Explorer over edges of wires ->
//!    [`brep_tools_wire_explorer`] (the BRepTools_WireExplorer reduction
//!    used by the offset pipeline).
//! 5. gp_Trsf -> glam::DAffine3; TopLoc_Location composition is carried by
//!    the DAffine3 values (the rcad location pool is BRep-bound; see the
//!    Perform note).
//! 6. GeomFill_Trihedron / GeomAbs_Shape / GeomFill_ApproxStyle /
//!    BRepFill_TransitionStyle are the local enums (the per-file enum
//!    precedent of brep_offset_api_make_pipe.rs).

use std::collections::HashMap;
use std::sync::Arc;

use glam::{DAffine3, DMat3, DVec3};

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve3, TrimmedCurve3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape, tshape_flags};

use crate::brep_fill::generator::ShapeKey;

/// OCCT Precision::Confusion().
const TOL_CONFUSION: f64 = CONFUSION;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored) — the generator.rs reduction.
fn shape_key(s: &Shape) -> ShapeKey {
    ShapeKey(s.ptr_id())
}

/// OCCT TopoDS_Shape::Reversed() — a copy with the reversed orientation.
pub(crate) fn shape_reversed(s: &Shape) -> Shape {
    let mut c = s.clone();
    c.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        _ => s.orientation,
    };
    c
}

/// OCCT gp_Trsf::SetValues(M(1, 1), ..., V.Z()) over the law frame — the
/// BRepFill_Pipe::Perform form (the brep_fill_location_law.rs re-host); a
/// null determinant keeps the OCCT Standard_ConstructionError failure.
fn fila_set_values(m: &GpMat, v: DVec3) -> DAffine3 {
    match crate::brep_fill::brep_fill_location_law::gp_trsf_set_values(&[
        m.mat[0][0], m.mat[0][1], m.mat[0][2], v.x, m.mat[1][0], m.mat[1][1], m.mat[1][2],
        v.y, m.mat[2][0], m.mat[2][1], m.mat[2][2], v.z,
    ]) {
        Ok(fila) => fila,
        Err(()) => panic!("gp_Trsf::SetValues, null determinant"),
    }
}

// ---------------------------------------------------------------------------
// GAP carriers — the sweep engine stack (architecture difference #1)
// ---------------------------------------------------------------------------

/// OCCT GeomFill_Trihedron (GeomFill_Trihedron.hxx L20-31).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomFillTrihedron {
    IsCorrectedFrenet,
    IsFixed,
    IsFrenet,
    IsConstantNormal,
    IsDarboux,
    IsGuideAC,
    IsGuidePlan,
    IsGuideACWithContact,
    IsGuidePlanWithContact,
    IsDiscreteTrihedron,
}

/// OCCT GeomAbs_Shape (continuity enum, GeomAbs_Shape.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsShape {
    AbsC0,
    AbsG1,
    AbsC1,
    AbsG2,
    AbsC2,
    AbsC3,
}

/// OCCT GeomFill_ApproxStyle (GeomFill_ApproxStyle.hxx L18-22).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum GeomFillApproxStyle {
    GeomFill_Section,
    GeomFill_Location,
}

/// OCCT BRepFill_TransitionStyle (BRepFill_TransitionStyle.hxx L23-27).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BRepFillTransitionStyle {
    TransitionStyle_Modified,
    TransitionStyle_Right,
    TransitionStyle_Round,
}

/// GAP: GeomFill_DiscreteTrihedron (TKGeomAlgo/GeomFill, parallel batch) —
/// the trihedron law constructed for GeomFill_IsDiscreteTrihedron (the
/// pipe_shell split-unit carrier).
use super::brep_fill_pipe_shell_b::GeomFillDiscreteTrihedron;

/// OCCT GeomFill_CurveAndTrihedron — the real batch-1 translation
/// (geomalgo/geomfill/curve_and_trihedron.rs); the ctor carries the
/// `Box<dyn TrihedronLaw>` law (the OCCT handle(GeomFill_TrihedronLaw)
/// TLaw of BRepFill_Pipe.cxx L228).
use crate::geomalgo::geomfill::curve_and_trihedron::CurveAndTrihedron as GeomFillCurveAndTrihedron;
use crate::geomalgo::geomfill::gp_mat::GpMat;
use crate::geomalgo::geomfill::location_law::LocationLaw;
use crate::geomalgo::geomfill::trihedron_law::TrihedronLaw;

/// GAP: BRepFill_LocationLaw (TKBool/BRepFill) — the spine location law
/// (BRepFill_Edge3DLaw result); not translated (plan D3).
pub struct BRepFillLocationLaw {
    /// OCCT: mySpine carried by the law.
    #[allow(dead_code)]
    my_spine: Shape,
}

impl BRepFillLocationLaw {
    /// OCCT BRepFill_LocationLaw::NbLaw().
    pub fn nb_law(&self) -> usize {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Edge(ii).
    pub fn edge(&self, _ii: usize) -> Shape {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Vertex(ii).
    pub fn vertex(&self, _ii: usize) -> Shape {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Law(ii) — the GeomFill location law.
    pub fn law(&self, _ii: usize) -> &GeomFillCurveAndTrihedron {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::IsClosed().
    pub fn is_closed(&self) -> bool {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::TransformInG0Law().
    pub fn transform_in_g0_law(&mut self) {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill not translated) — see file header")
    }
}

/// GAP: BRepFill_Edge3DLaw (TKBool/BRepFill) — the edge 3D law constructor.
pub struct BRepFillEdge3DLaw;

impl BRepFillEdge3DLaw {
    /// OCCT new BRepFill_Edge3DLaw(Spine, Law).
    pub fn new(the_spine: &Shape, _the_law: &GeomFillCurveAndTrihedron) -> BRepFillLocationLaw {
        BRepFillLocationLaw {
            my_spine: the_spine.clone(),
        }
    }
}

/// GAP: BRepFill_ShapeLaw (TKBool/BRepFill) — the section law.
pub struct BRepFillShapeLaw;

impl BRepFillShapeLaw {
    /// OCCT new BRepFill_ShapeLaw(Vertex).
    pub fn new(_the_vertex: &Shape) -> BRepFillShapeLaw {
        panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT new BRepFill_ShapeLaw(Wire, Skip) / (Wire).
    pub fn new_with_wire(_the_wire: &Shape, _skip: bool) -> BRepFillShapeLaw {
        panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_ShapeLaw::NbLaw().
    pub fn nb_law(&self) -> usize {
        panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_ShapeLaw::Edge(ii).
    pub fn edge(&self, _ii: usize) -> Shape {
        panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_ShapeLaw::Vertex(ii, Tol).
    pub fn vertex(&self, _ii: usize, _tol: f64) -> Shape {
        panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill not translated) — see file header")
    }
}

/// GAP: BRepFill_SectionLaw (TKBool/BRepFill) — the section law base class
/// (BRepFill_ShapeLaw derives from it); not translated (plan D3).
pub struct BRepFillSectionLaw;

/// OCCT handle up-cast Handle(BRepFill_ShapeLaw) -> Handle(BRepFill_SectionLaw).
impl From<BRepFillShapeLaw> for BRepFillSectionLaw {
    fn from(_: BRepFillShapeLaw) -> Self {
        BRepFillSectionLaw
    }
}

/// GAP: BRepFill_Sweep (TKBool/BRepFill) — the sweep engine; not translated
/// (plan D3).  The constructor and the Build/Shape/SubShape/Sections/
/// InterFaces/Tape surface keep the OCCT call form.
pub struct BRepFillSweep;

impl BRepFillSweep {
    /// OCCT BRepFill_Sweep(Section, Law, WithTrsfm) (BRepFill_Sweep.hxx):
    /// the section is carried through the BRepFill_SectionLaw base handle.
    pub fn new(
        _the_section: BRepFillSectionLaw,
        _the_law: &BRepFillLocationLaw,
        _with_trsfm: bool,
    ) -> Self {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::SetForceApproxC1(ForceApproxC1).
    pub fn set_force_approx_c1(&mut self, _force_approx_c1: bool) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::SetBounds(WFirst, WLast).
    pub fn set_bounds(&mut self, _w_first: &Shape, _w_last: &Shape) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::SetTolerance(Tol).
    pub fn set_tolerance(&mut self, _tol: f64) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::SetAngularControl(AngleMin, AngleMax).
    pub fn set_angular_control(&mut self, _angle_min: f64, _angle_max: f64) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::Build(...) (BRepFill_Sweep.hxx L93-104).
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        &mut self,
        _reversed_edges: &mut HashMap<ShapeKey, Shape>,
        _tapes: &mut HashMap<ShapeKey, ShapeArray2>,
        _rails: &mut HashMap<ShapeKey, ShapeArray2>,
        _transition: BRepFillTransitionStyle,
        _continuity: GeomAbsShape,
        _approx: GeomFillApproxStyle,
        _degmax: i32,
        _segmax: i32,
    ) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::ErrorOnSurface().
    pub fn error_on_surface(&self) -> f64 {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::SubShape().
    pub fn sub_shape(&self) -> Option<ShapeArray2> {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::Sections().
    pub fn sections(&self) -> Option<ShapeArray2> {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::InterFaces().
    pub fn interfaces(&self) -> Option<ShapeArray2> {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_Sweep::Tape(Ind) (BRepFill_Sweep.hxx).
    pub fn tape(&self, _ind: usize) -> Shape {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill not translated) — see file header")
    }
}

/// GAP: BRepFill_SectionPlacement (TKBool/BRepFill).
pub struct BRepFillSectionPlacement;

impl BRepFillSectionPlacement {
    /// OCCT new BRepFill_SectionPlacement(Law, Section).
    pub fn new(_the_law: &BRepFillLocationLaw, _the_section: &Shape) -> Self {
        panic!("GAP: BRepFill_SectionPlacement (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_SectionPlacement::Transformation().
    pub fn transformation(&self) -> DAffine3 {
        panic!("GAP: BRepFill_SectionPlacement (TKBool/BRepFill not translated) — see file header")
    }
}

/// GAP: ShapeUpgrade_RemoveLocations (TKShHealing/ShapeUpgrade).
pub struct ShapeUpgradeRemoveLocations {
    #[allow(dead_code)]
    my_result: Shape,
}

impl Default for ShapeUpgradeRemoveLocations {
    /// OCCT ShapeUpgrade_RemoveLocations() — the default constructor.
    fn default() -> Self {
        ShapeUpgradeRemoveLocations {
            my_result: Shape::null(),
        }
    }
}

impl ShapeUpgradeRemoveLocations {
    /// OCCT ShapeUpgrade_RemoveLocations::SetRemoveLevel(Level).
    pub fn set_remove_level(&mut self, _the_level: ShapeType) {
        // OCCT stores the level; the Remove body is the GAP.
        panic!("GAP: ShapeUpgrade_RemoveLocations (TKShHealing not translated) — see file header")
    }
    /// OCCT ShapeUpgrade_RemoveLocations::Remove(S).
    pub fn remove(&mut self, _the_s: &Shape) {
        panic!("GAP: ShapeUpgrade_RemoveLocations (TKShHealing not translated) — see file header")
    }
    /// OCCT ShapeUpgrade_RemoveLocations::GetResult().
    pub fn get_result(&self) -> Shape {
        panic!("GAP: ShapeUpgrade_RemoveLocations (TKShHealing not translated) — see file header")
    }
}

/// GAP: BRepBuilderAPI_Transform (TKBRep/BRepBuilderAPI) — the copying
/// transform (the true-copy form consumed by Perform).
pub struct BRepBuilderAPITransform;

impl BRepBuilderAPITransform {
    /// OCCT BRepBuilderAPI_Transform(S, T, Copy = true) — the copy form.
    pub fn new_copy(_the_shape: &Shape, _the_trsf: DAffine3) -> Shape {
        panic!("GAP: BRepBuilderAPI_Transform (copy form) — see file header")
    }
}

// ---------------------------------------------------------------------------
// Local carriers (architecture differences #2, #4, #5)
// ---------------------------------------------------------------------------

/// OCCT NCollection_HArray2<TopoDS_Shape> (architecture difference #2) —
/// 1-based row/col indexing as in OCCT.
#[derive(Clone)]
pub struct ShapeArray2 {
    lo_row: usize,
    hi_row: usize,
    lo_col: usize,
    hi_col: usize,
    data: Vec<Shape>,
}

impl ShapeArray2 {
    /// OCCT new NCollection_HArray2(Low1, Up1, Low2, Up2).
    pub fn new(low1: usize, up1: usize, low2: usize, up2: usize) -> Self {
        let n_rows = up1 - low1 + 1;
        let n_cols = up2 - low2 + 1;
        ShapeArray2 {
            lo_row: low1,
            hi_row: up1,
            lo_col: low2,
            hi_col: up2,
            data: vec![Shape::null(); n_rows * n_cols],
        }
    }
    /// OCCT ColLength().
    pub fn col_length(&self) -> usize {
        self.hi_row - self.lo_row + 1
    }
    /// OCCT RowLength().
    pub fn row_length(&self) -> usize {
        self.hi_col - self.lo_col + 1
    }
    /// OCCT UpperRow().
    pub fn upper_row(&self) -> usize {
        self.hi_row
    }
    /// OCCT UpperCol().
    pub fn upper_col(&self) -> usize {
        self.hi_col
    }
    /// OCCT Value(i, j).
    pub fn value(&self, i: usize, j: usize) -> Shape {
        let r = i - self.lo_row;
        let c = j - self.lo_col;
        self.data[r * self.row_length() + c].clone()
    }
    /// OCCT SetValue(i, j, v).
    pub fn set_value(&mut self, i: usize, j: usize, v: Shape) {
        let r = i - self.lo_row;
        let c = j - self.lo_col;
        let row_len = self.hi_col - self.lo_col + 1;
        self.data[r * row_len + c] = v;
    }
}

/// OCCT TopoDS_Iterator(S) — the direct children in stored order
/// (architecture difference #4).
pub(crate) fn topods_iterator(s: &Shape) -> Vec<Shape> {
    match s.data.as_ref() {
        TShape::Compound(children) => children.clone(),
        TShape::CompSolid(children) => children.clone(),
        TShape::Solid(sd) => sd.shells.clone(),
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Face(fd) => {
            let mut ws = vec![fd.outer_wire.clone()];
            ws.extend(fd.inner_wires.iter().cloned());
            ws
        }
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT BRepTools_WireExplorer(Spine) — the ordered spine edges
/// (architecture difference #4).
fn brep_tools_wire_explorer(w: &Shape) -> Vec<Shape> {
    match w.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT BRepTools_WireExplorer::CurrentVertex() — the vertex closing the
/// current edge of the walk.
fn wire_explorer_current_vertex(edges: &[Shape], index: usize) -> Shape {
    let (first, last) = match edges[index].data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    if index + 1 == edges.len() {
        // the last vertex of the wire.
        if edges[index].orientation == Orientation::Reversed {
            first
        } else {
            last
        }
    } else {
        // the first vertex of the next edge.
        match edges[index + 1].data.as_ref() {
            TShape::Edge(ed) => {
                if edges[index + 1].orientation == Orientation::Reversed {
                    ed.last.clone()
                } else {
                    ed.first.clone()
                }
            }
            _ => Shape::null(),
        }
    }
}

/// OCCT TopExp::Vertices(Edge, Vf, Vl) — the end vertices in traversal
/// order (orientation-aware).
pub(crate) fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let (first, last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => (Shape::null(), Shape::null()),
    };
    if e.orientation == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}

/// OCCT BRep_Tool::Degenerated(E).
fn brep_tool_degenerated(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Curve(E, f, l) — the 3D curve and range.
fn brep_tool_curve(e: &Shape) -> Option<(Curve3, f64, f64)> {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// OCCT BRep_Tool::IsClosed(S).
fn brep_tool_is_closed(s: &Shape) -> bool {
    match s.data.as_ref() {
        TShape::Wire(wd) => wd.flags & tshape_flags::CLOSED != 0,
        TShape::Shell(sd) => sd.flags & tshape_flags::CLOSED != 0,
        TShape::Edge(ed) => ed.flags & tshape_flags::CLOSED != 0,
        _ => false,
    }
}

/// OCCT TopAbs_IN (TopAbs_State enum value).
pub(crate) const TOPABS_IN: u8 = 0;

/// OCCT BRepFill_Pipe (hxx L49-148) — create a shape by sweeping a shape
/// (the profile) along a wire (the spine).
pub struct BRepFillPipe {
    my_spine: Shape,   // OCCT: mySpine (TopoDS_Wire)
    my_profile: Shape, // OCCT: myProfile
    my_shape: Shape,   // OCCT: myShape
    my_trsf: DAffine3, // OCCT: myTrsf (gp_Trsf)
    my_loc: Option<BRepFillLocationLaw>, // OCCT: myLoc (None = null handle)
    my_sections: Option<ShapeArray2>,    // OCCT: mySections
    my_faces: Option<ShapeArray2>,       // OCCT: myFaces
    my_edges: Option<ShapeArray2>,       // OCCT: myEdges
    my_reversed_edges: HashMap<ShapeKey, Shape>, // OCCT: myReversedEdges
    my_tapes: HashMap<ShapeKey, ShapeArray2>,    // OCCT: myTapes
    my_rails: HashMap<ShapeKey, ShapeArray2>,    // OCCT: myRails
    my_cur_index_of_section_edge: i32,   // OCCT: myCurIndexOfSectionEdge
    my_first: Shape,   // OCCT: myFirst
    my_last: Shape,    // OCCT: myLast
    my_gen_map: HashMap<ShapeKey, Vec<Shape>>, // OCCT: myGenMap
    my_degmax: i32,    // OCCT: myDegmax
    my_segmax: i32,    // OCCT: mySegmax
    my_continuity: GeomAbsShape, // OCCT: myContinuity
    my_mode: GeomFillTrihedron,  // OCCT: myMode
    my_force_approx_c1: bool,    // OCCT: myForceApproxC1
    my_error_on_surf: f64,       // OCCT: myErrorOnSurf
}

/// OCCT static UpdateMap (cxx L64-94).
fn update_map(
    the_key: &Shape,
    the_value: &Shape,
    the_map: &mut HashMap<ShapeKey, Vec<Shape>>,
) -> bool {
    // OCCT L70-74: if (!theMap.IsBound(theKey)) Bind(theKey, empty).
    let the_list = the_map.entry(shape_key(the_key)).or_default();
    // OCCT L76-86: the list scan for IsSame.
    let mut found = false;
    for v in the_list.iter() {
        if the_value.is_same(v) {
            found = true;
            break;
        }
    }

    // OCCT L88-91.
    if !found {
        the_list.push(the_value.clone());
    }

    !found
}

/// OCCT static UpdateTolFromTopOrBottomPCurve (cxx L96-144) — GAP: the
/// Adaptor3d_CurveOnSurface evaluation is not translated.
fn update_tol_from_top_or_bottom_pcurve(_a_face: &Shape, _an_edge: &Shape) {
    // OCCT L99-113: the pcurve / surface / ConS chain requires
    // BRep_Tool::CurveOnSurface + Geom2dAdaptor_Curve + GeomAdaptor_Surface
    // + Adaptor3d_CurveOnSurface (not translated — plan D3).
    panic!("GAP: UpdateTolFromTopOrBottomPCurve requires Adaptor3d_CurveOnSurface — see file header")
}

impl BRepFillPipe {
    /// OCCT BRepFill_Pipe::BRepFill_Pipe() (cxx L148-157).
    pub fn new_empty() -> Self {
        BRepFillPipe {
            my_spine: Shape::null(),
            my_profile: Shape::null(),
            my_shape: Shape::null(),
            my_trsf: DAffine3::IDENTITY,
            my_loc: None,
            my_sections: None,
            my_faces: None,
            my_edges: None,
            my_reversed_edges: HashMap::new(),
            my_tapes: HashMap::new(),
            my_rails: HashMap::new(),
            // OCCT L156: myCurIndexOfSectionEdge = 1.
            my_cur_index_of_section_edge: 1,
            my_first: Shape::null(),
            my_last: Shape::null(),
            my_gen_map: HashMap::new(),
            // OCCT L150-154.
            my_degmax: 11,
            my_segmax: 100,
            my_continuity: GeomAbsShape::AbsC2,
            my_mode: GeomFillTrihedron::IsCorrectedFrenet,
            my_force_approx_c1: false,
            my_error_on_surf: 0.0,
        }
    }

    /// OCCT BRepFill_Pipe::BRepFill_Pipe(Spine, Profile, aMode,
    /// ForceApproxC1, KPart) (cxx L161-189).
    pub fn new(
        the_spine: &Shape,
        the_profile: &Shape,
        a_mode: GeomFillTrihedron,
        force_approx_c1: bool,
        k_part: bool,
    ) -> Self {
        // OCCT L171-176.
        let my_mode = GeomFillTrihedron::IsCorrectedFrenet;
        let my_mode = if a_mode == GeomFillTrihedron::IsFrenet
            || a_mode == GeomFillTrihedron::IsCorrectedFrenet
            || a_mode == GeomFillTrihedron::IsDiscreteTrihedron
        {
            a_mode
        } else {
            my_mode
        };

        // OCCT L178-182.
        let my_continuity = if my_mode == GeomFillTrihedron::IsDiscreteTrihedron {
            GeomAbsShape::AbsC0
        } else {
            GeomAbsShape::AbsC2
        };

        // OCCT L184-186.
        let my_force_approx_c1 = force_approx_c1;
        let my_cur_index_of_section_edge = 1;

        let mut r = BRepFillPipe {
            my_spine: Shape::null(),
            my_profile: Shape::null(),
            my_shape: Shape::null(),
            my_trsf: DAffine3::IDENTITY,
            my_loc: None,
            my_sections: None,
            my_faces: None,
            my_edges: None,
            my_reversed_edges: HashMap::new(),
            my_tapes: HashMap::new(),
            my_rails: HashMap::new(),
            my_cur_index_of_section_edge,
            my_first: Shape::null(),
            my_last: Shape::null(),
            my_gen_map: HashMap::new(),
            // OCCT L168-169: myDegmax = 11, mySegmax = 100.
            my_degmax: 11,
            my_segmax: 100,
            my_continuity,
            my_mode,
            my_force_approx_c1,
            my_error_on_surf: 0.0,
        };

        // OCCT L188: Perform(Spine, Profile, KPart).
        r.perform(the_spine, the_profile, k_part);
        r
    }

    /// OCCT BRepFill_Pipe::Perform(Spine, Profile, KPart) (cxx L193-313).
    pub fn perform(&mut self, the_spine: &Shape, the_profile: &Shape, _k_part: bool) {
        // OCCT L198-200: mySections/myFaces/myEdges.Nullify().
        self.my_sections = None;
        self.my_faces = None;
        self.my_edges = None;

        // OCCT L202-203.
        self.my_spine = the_spine.clone();
        self.my_profile = the_profile.clone();

        // OCCT L205.
        self.define_real_segmax();

        // OCCT L210-224: the trihedron law of the mode.
        let t_law: Option<Box<dyn TrihedronLaw>> = match self.my_mode {
            GeomFillTrihedron::IsFrenet => Some(Box::new(
                crate::geomalgo::geomfill::frenet::Frenet::new(),
            )),
            GeomFillTrihedron::IsCorrectedFrenet => Some(Box::new(
                crate::geomalgo::geomfill::corrected_frenet::CorrectedFrenet::new(),
            )),
            GeomFillTrihedron::IsDiscreteTrihedron => {
                Some(GeomFillDiscreteTrihedron::new().into_trihedron_law())
            }
            _ => None,
        };
        // OCCT L225-231.  The OCCT ctor takes the possibly-null TLaw handle;
        // the rcad Box has no null state (the ctor guard above restricts
        // myMode to the three switch kinds).
        let loc = GeomFillCurveAndTrihedron::new(t_law.expect("null TLaw"));
        self.my_loc = Some(BRepFillEdge3DLaw::new(the_spine, &loc));
        let my_loc = self.my_loc.as_ref().expect("myLoc");
        if my_loc.nb_law() == 0 {
            return; // Degenerated case
        }

        // OCCT L233-234.
        let place = BRepFillSectionPlacement::new(
            self.my_loc.as_ref().expect("myLoc"),
            the_profile,
        );
        self.my_trsf = place.transformation();

        // OCCT L236-240: TopLoc_Location Loc2(myTrsf), Loc1 =
        // Profile.Location(); TheProf.Location(Loc2.Multiplied(Loc1)).
        // Architecture difference #5: the rcad location pool is BRep-bound,
        // so the composed location is carried by the local affine value and
        // TheProf keeps the profile shape.
        let loc1 = DAffine3::IDENTITY;
        let loc2 = self.my_trsf;
        let mut the_prof = self.my_profile.clone();
        let _loc_the_prof = loc2 * loc1;
        let _ = &mut the_prof;

        // OCCT L242-265: construct First && Last Shape — the law frame at
        // the first parameter.
        let my_loc = self.my_loc.as_ref().expect("myLoc");
        let law1 = my_loc.law(1);
        let mut m = GpMat::identity();
        let mut v = DVec3::ZERO;
        let mut first = 0.0;
        let mut last = 0.0;
        law1.get_domain(&mut first, &mut last);
        law1.d0(first, &mut m, &mut v);
        // OCCT fila.SetValues(M(1, 1), ..., V.Z()).
        let mut fila = fila_set_values(&m, v);

        // OCCT L264-271.
        fila = self.my_trsf * fila;
        let loc_first = fila;
        let mut my_first = self.my_profile.clone();
        if loc_first != DAffine3::IDENTITY {
            // OCCT L270: myFirst = BRepBuilderAPI_Transform(myProfile,
            // fila, true) // copy.
            my_first = BRepBuilderAPITransform::new_copy(&self.my_profile, fila);
        }

        // OCCT L273-276: ShapeUpgrade_RemoveLocations.
        let mut rem_loc = ShapeUpgradeRemoveLocations::default();
        rem_loc.set_remove_level(ShapeType::Compound);
        rem_loc.remove(&my_first);
        self.my_first = rem_loc.get_result();

        // OCCT L278-293: the law frame at the last parameter.
        let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
        let my_loc = self.my_loc.as_ref().expect("myLoc");
        let law_n = my_loc.law(nb_law);
        let mut m = GpMat::identity();
        let mut v = DVec3::ZERO;
        law_n.get_domain(&mut first, &mut last);
        law_n.d0(last, &mut m, &mut v);
        // OCCT fila.SetValues(M(1, 1), ..., V.Z()).
        let mut fila = fila_set_values(&m, v);

        // OCCT L293-307.
        fila = self.my_trsf * fila;
        let loc_last = fila;
        let my_loc = self.my_loc.as_ref().expect("myLoc");
        if !my_loc.is_closed() || loc_last != loc_first {
            let mut my_last = self.my_profile.clone();
            if loc_last != DAffine3::IDENTITY {
                // OCCT L301: myLast = BRepBuilderAPI_Transform(myProfile,
                // fila, true) // copy.
                my_last = BRepBuilderAPITransform::new_copy(&self.my_profile, fila);
            }
            self.my_last = my_last;
        } else {
            // OCCT L306: myLast = myFirst.
            self.my_last = self.my_first.clone();
        }

        // OCCT L309-310.
        let mut rem_loc = ShapeUpgradeRemoveLocations::default();
        rem_loc.remove(&self.my_last.clone());
        self.my_last = rem_loc.get_result();

        // OCCT L312.
        let the_prof = the_prof.clone();
        self.my_shape = self.make_shape(&the_prof, &self.my_profile.clone(), &self.my_first.clone(), &self.my_last.clone());
    }

    /// OCCT BRepFill_Pipe::Spine() (cxx L317-320).
    pub fn spine(&self) -> Shape {
        self.my_spine.clone()
    }

    /// OCCT BRepFill_Pipe::Profile() (cxx L324-327).
    pub fn profile(&self) -> Shape {
        self.my_profile.clone()
    }

    /// OCCT BRepFill_Pipe::Shape() (cxx L331-334).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepFill_Pipe::ErrorOnSurface() (cxx L338-341).
    pub fn error_on_surface(&self) -> f64 {
        self.my_error_on_surf
    }

    /// OCCT BRepFill_Pipe::FirstShape() (cxx L345-348).
    pub fn first_shape(&self) -> Shape {
        self.my_first.clone()
    }

    /// OCCT BRepFill_Pipe::LastShape() (cxx L352-355).
    pub fn last_shape(&self) -> Shape {
        self.my_last.clone()
    }

    /// OCCT BRepFill_Pipe::Generated(S, L) (cxx L359-367).
    pub fn generated(&self, the_shape: &Shape, the_list: &mut Vec<Shape>) {
        the_list.clear();

        if let Some(l) = self.my_gen_map.get(&shape_key(the_shape)) {
            *the_list = l.clone();
        }
    }

    /// OCCT BRepFill_Pipe::Face(ESpine, EProfile) (cxx L371-411).
    pub fn face(&mut self, e_spine: &Shape, e_profile: &Shape) -> Shape {
        let mut the_face = Shape::null();

        // OCCT L375-378.
        if brep_tool_degenerated(e_profile) {
            return the_face;
        }

        let mut count = 0;

        // OCCT L385-390: search if EProfile is an edge of myProfile.
        let iprof = self.find_edge(&self.my_profile.clone(), e_profile, &mut count);
        if iprof == 0 {
            panic!("Standard_DomainError: BRepFill_Pipe::Face : Edge not in the Profile");
        }

        // OCCT L396-407: search if ESpine is an edge of mySpine.
        let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
        let mut ispin = 0;
        for ii in 1..=nb_law {
            if ispin == 0 && e_spine.is_same(&self.my_loc.as_ref().expect("myLoc").edge(ii)) {
                ispin = ii;
            }
        }
        if ispin == 0 {
            panic!("Standard_DomainError: BRepFill_Pipe::Edge  : Edge not in the Spine");
        }

        // OCCT L409-410.
        the_face = self
            .my_faces
            .as_ref()
            .expect("myFaces")
            .value(iprof as usize, ispin as usize);
        the_face
    }

    /// OCCT BRepFill_Pipe::Edge(ESpine, VProfile) (cxx L415-453).
    pub fn edge(&mut self, e_spine: &Shape, v_profile: &Shape) -> Shape {
        let mut count = 0;

        // OCCT L422-426: search if VProfile is a Vertex of myProfile.
        let iprof = self.find_vertex(&self.my_profile.clone(), v_profile, &mut count);
        if iprof == 0 {
            panic!("Standard_DomainError: BRepFill_Pipe::Edge : Vertex not in the Profile");
        }

        // OCCT L433-444: search if ESpine is an edge of mySpine.
        let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
        let mut ispin = 0;
        for ii in 1..=nb_law {
            if ispin == 0 && e_spine.is_same(&self.my_loc.as_ref().expect("myLoc").edge(ii)) {
                ispin = ii;
            }
        }
        if ispin == 0 {
            panic!("Standard_DomainError: BRepFill_Pipe::Edge  : Edge not in the Spine");
        }

        // OCCT L449-452: generate the corresponding Shape.
        self.my_edges
            .as_ref()
            .expect("myEdges")
            .value(iprof as usize, ispin as usize)
    }

    /// OCCT BRepFill_Pipe::Section(VSpine) (cxx L457-493).
    pub fn section(&self, v_spine: &Shape) -> Shape {
        // OCCT L471-482: search the vertex among the spine vertices.
        let nb = self.my_loc.as_ref().expect("myLoc").nb_law() + 1;
        let mut ispin = 0;
        for ii in 1..=nb {
            if ispin == 0 && v_spine.is_same(&self.my_loc.as_ref().expect("myLoc").vertex(ii)) {
                ispin = ii;
            }
        }
        if ispin == 0 {
            panic!("Standard_DomainError: BRepFill_Pipe::Section  : Vertex not in the Spine");
        }

        // OCCT L484-492: B.MakeCompound(Comp); Add all the sections.
        let mut pool = BRep::new();
        let b = BRepBuilder::new();
        let mut comp = pool.add_tcompound(Vec::new());
        for ii in 1..=self.my_sections.as_ref().expect("mySections").col_length() {
            let s = self.my_sections.as_ref().expect("mySections").value(ii, ispin);
            add_to_compound(&mut pool, &b, &mut comp, &s);
        }

        comp
    }

    /// OCCT BRepFill_Pipe::PipeLine(Point) (cxx L500-526).
    pub fn pipe_line(&mut self, point: DVec3) -> Shape {
        // OCCT L502-505: P = Point; P.Transform(myTrsf).
        let p = self.my_trsf.transform_point3(point);

        // OCCT L507-508.
        let mut pool = BRep::new();
        let vertex_section = pool.add_tvertex_unique(p);
        let section = BRepFillShapeLaw::new(&vertex_section);
        let _ = &section;

        // OCCT L511-520: sweeping.
        let mut mk_sw = BRepFillSweep::new(
            section.into(),
            self.my_loc.as_ref().expect("myLoc"),
            true,
        );
        mk_sw.set_force_approx_c1(self.my_force_approx_c1);
        mk_sw.build(
            &mut self.my_reversed_edges,
            &mut self.my_tapes,
            &mut self.my_rails,
            BRepFillTransitionStyle::TransitionStyle_Modified,
            self.my_continuity,
            GeomFillApproxStyle::GeomFill_Location,
            self.my_degmax,
            self.my_segmax,
        );
        // OCCT L521-524.
        let a_local_shape = mk_sw.shape();
        self.my_error_on_surf = mk_sw.error_on_surface();
        self.build_history(&mk_sw, &vertex_section);
        a_local_shape
    }

    /// OCCT BRepFill_Pipe::MakeShape(S, theOriginalS, FirstShape,
    /// LastShape) (cxx L530-854).
    fn make_shape(
        &mut self,
        s: &Shape,
        the_original_s: &Shape,
        first_shape: &Shape,
        last_shape: &Shape,
    ) -> Shape {
        let mut pool = BRep::new();
        let b = BRepBuilder::new();
        let mut explode = false;
        let mut the_s = s.clone();
        let mut the_first = first_shape.clone();
        let mut the_last = last_shape.clone();
        // OCCT L543-546.
        let initial_length = match &self.my_faces {
            Some(f) => f.col_length(),
            None => 0,
        };

        // OCCT L548: TopLoc_Location BackLoc(myTrsf.Inverted()).
        let _back_loc = self.my_trsf.inverse();

        // OCCT L558-624: create the result empty.
        let mut result = Shape::null();
        match s.shape_type() {
            ShapeType::Vertex => {
                // OCCT L562: B.MakeWire(result).
                result = pool.add_twire(Vec::new());
            }
            ShapeType::Edge => {
                // OCCT L567-588.
                result = pool.add_tshell(Vec::new());
                let w = pool.add_twire(vec![s.clone()]);
                let _ = &b;
                // OCCT L571: W.Closed(BRep_Tool::IsClosed(S)).
                set_closed_flag(&mut pool, &w, brep_tool_is_closed(s));
                the_s = w;
                if !the_first.is_null() {
                    let w = pool.add_twire(vec![the_first.clone()]);
                    set_closed_flag(&mut pool, &w, brep_tool_is_closed(&the_first));
                    the_first = w;
                }
                if !the_last.is_null() {
                    let w = pool.add_twire(vec![the_last.clone()]);
                    set_closed_flag(&mut pool, &w, brep_tool_is_closed(&the_last));
                    the_last = w;
                }
                // OCCT L587.
                set_closed_flag(&mut pool, &result, brep_tool_is_closed(&result));
            }
            ShapeType::Wire => {
                // OCCT L592.
                result = pool.add_tshell(Vec::new());
            }
            ShapeType::Face => {
                // OCCT L596-604.
                result = pool.add_tshell(Vec::new());
                explode = true;
                if !brep_tool_is_closed(&self.my_spine) && !the_first.is_null() {
                    add_to_shell(&mut pool, &b, &mut result, &shape_reversed(&the_first));
                }
                set_closed_flag(&mut pool, &result, brep_tool_is_closed(&result));
            }
            ShapeType::Shell => {
                // OCCT L607-609.
                result = pool.add_tcompsolid(Vec::new());
                explode = true;
            }
            ShapeType::Solid | ShapeType::CompSolid => {
                // OCCT L614.
                panic!("Standard_DomainError: BRepFill_Pipe::profile contains solids");
            }
            ShapeType::Compound => {
                // OCCT L618-620.
                result = pool.add_tcompound(Vec::new());
                explode = true;
            }
            _ => {}
        }

        if explode {
            // OCCT L629-669: add the subshapes.
            let it_first_vals = if !the_first.is_null() {
                topods_iterator(&the_first)
            } else {
                Vec::new()
            };
            let it_last_vals = if !the_last.is_null() {
                topods_iterator(&the_last)
            } else {
                Vec::new()
            };
            let it_vals = topods_iterator(s);
            let it_orig_vals = topods_iterator(the_original_s);

            for idx in 0..it_vals.len() {
                let mut first = Shape::null();
                let mut last = Shape::null();
                if !the_first.is_null() {
                    first = it_first_vals[idx].clone();
                }
                if !the_last.is_null() {
                    last = it_last_vals[idx].clone();
                }
                if the_s.shape_type() == ShapeType::Face {
                    // OCCT L654: the face recursion drops the sub-result.
                    let _ = self.make_shape(&it_vals[idx], &it_orig_vals[idx], &first, &last);
                } else {
                    let sub = self.make_shape(&it_vals[idx], &it_orig_vals[idx], &first, &last);
                    // OCCT L658: B.Add(result, ...).
                    add_to_result(&mut pool, &b, &mut result, &sub);
                }
            }
        } else {
            // OCCT L674-702: TheS is a VERTEX — sweep the section.
            if the_s.shape_type() == ShapeType::Vertex {
                let section = BRepFillShapeLaw::new(&the_s);
                let _ = &section;
                // OCCT L677: BRepFill_Sweep MkSw(Section, myLoc, true).
                let mut mk_sw = BRepFillSweep::new(
                    section.into(),
                    self.my_loc.as_ref().expect("myLoc"),
                    true,
                );
                mk_sw.set_force_approx_c1(self.my_force_approx_c1);
                mk_sw.build(
                    &mut self.my_reversed_edges,
                    &mut self.my_tapes,
                    &mut self.my_rails,
                    BRepFillTransitionStyle::TransitionStyle_Modified,
                    self.my_continuity,
                    GeomFillApproxStyle::GeomFill_Location,
                    self.my_degmax,
                    self.my_segmax,
                );
                result = mk_sw.shape();
                update_map(the_original_s, &result, &mut self.my_gen_map);
                self.my_error_on_surf = mk_sw.error_on_surface();

                // OCCT L691-699.
                let a_sections = mk_sw.sections();
                if let Some(a_sections) = a_sections {
                    let a_v_last = a_sections.upper_col();
                    self.my_first = a_sections.value(1, 1);
                    self.my_last = a_sections.value(1, a_v_last);
                }

                self.build_history(&mk_sw, the_original_s);
            }

            // OCCT L704-790: TheS is a WIRE — sweep with bounds.
            if the_s.shape_type() == ShapeType::Wire {
                let section = BRepFillShapeLaw::new_with_wire(&the_s, true);
                let _ = &section;
                // OCCT L707-709: BRepFill_Sweep MkSw(Section, myLoc, true);
                // MkSw.SetBounds(TheFirst, TheLast).
                let mut mk_sw = BRepFillSweep::new(
                    section.into(),
                    self.my_loc.as_ref().expect("myLoc"),
                    true,
                );
                mk_sw.set_bounds(&the_first, &the_last);
                mk_sw.set_force_approx_c1(self.my_force_approx_c1);
                mk_sw.build(
                    &mut self.my_reversed_edges,
                    &mut self.my_tapes,
                    &mut self.my_rails,
                    BRepFillTransitionStyle::TransitionStyle_Modified,
                    self.my_continuity,
                    GeomFillApproxStyle::GeomFill_Location,
                    self.my_degmax,
                    self.my_segmax,
                );
                result = mk_sw.shape();
                self.my_error_on_surf = mk_sw.error_on_surface();

                // OCCT L722-727: labeling of elements (first sweep).
                if self.my_sections.is_none() && self.my_faces.is_none() && self.my_edges.is_none()
                {
                    self.my_faces = mk_sw.sub_shape();
                    self.my_sections = mk_sw.sections();
                    self.my_edges = mk_sw.interfaces();
                } else {
                    // OCCT L730-787: concatenate the arrays.
                    // myFaces = Somme.
                    let aux = mk_sw.sub_shape().expect("Aux");
                    let my_faces = self.my_faces.as_ref().expect("myFaces");
                    let length = aux.col_length() + my_faces.col_length();
                    let mut somme = ShapeArray2::new(1, length, 1, aux.row_length());
                    for jj in 1..=my_faces.row_length() {
                        for ii in 1..=my_faces.col_length() {
                            somme.set_value(ii, jj, my_faces.value(ii, jj));
                        }
                        let mut ii = my_faces.col_length() + 1;
                        for kk in 1..=aux.col_length() {
                            somme.set_value(ii, jj, aux.value(kk, jj));
                            ii += 1;
                        }
                    }
                    self.my_faces = Some(somme);

                    // mySections = Somme.
                    let aux = mk_sw.sections().expect("Aux");
                    let my_sections = self.my_sections.as_ref().expect("mySections");
                    let length = aux.col_length() + my_sections.col_length();
                    let mut somme = ShapeArray2::new(1, length, 1, aux.row_length());
                    for jj in 1..=my_sections.row_length() {
                        for ii in 1..=my_sections.col_length() {
                            somme.set_value(ii, jj, my_sections.value(ii, jj));
                        }

                        // OCCT L761.
                        self.my_cur_index_of_section_edge =
                            (my_sections.col_length() + 1) as i32;

                        let mut ii = my_sections.col_length() + 1;
                        for kk in 1..=aux.col_length() {
                            somme.set_value(ii, jj, aux.value(kk, jj));
                            ii += 1;
                        }
                    }
                    self.my_sections = Some(somme);

                    // myEdges = Somme.
                    let aux = mk_sw.interfaces().expect("Aux");
                    let my_edges = self.my_edges.as_ref().expect("myEdges");
                    let length = aux.col_length() + my_edges.col_length();
                    let mut somme = ShapeArray2::new(1, length, 1, aux.row_length());
                    for jj in 1..=my_edges.row_length() {
                        for ii in 1..=my_edges.col_length() {
                            somme.set_value(ii, jj, my_edges.value(ii, jj));
                        }
                        let mut ii = my_edges.col_length() + 1;
                        for kk in 1..=aux.col_length() {
                            somme.set_value(ii, jj, aux.value(kk, jj));
                            ii += 1;
                        }
                    }
                    self.my_edges = Some(somme);
                }

                self.build_history(&mk_sw, the_original_s);
            }
        }

        // OCCT L793-849: the FACE case builds a solid.
        if the_s.shape_type() == ShapeType::Face {
            // OCCT L797-802: rebuild the top faces.
            let result_faces = topods_iterator(&result);
            for a_face in result_faces {
                let _ = self.rebuild_top_or_bottom_face(&shape_reversed(&a_face), true);
            }

            // OCCT L804-818: add the tape faces.
            let my_faces = self.my_faces.as_ref().expect("myFaces");
            for ii in (initial_length + 1)..=my_faces.col_length() {
                for jj in 1..=my_faces.row_length() {
                    let f = my_faces.value(ii, jj);
                    if f.shape_type() == ShapeType::Face && !f.is_null() {
                        add_to_result(&mut pool, &b, &mut result, &f);
                    }
                }
            }

            // OCCT L820-829: add the last face for an open spine.
            if !brep_tool_is_closed(&self.my_spine) {
                let _ = self.rebuild_top_or_bottom_face(&the_last, false);
                add_to_result(&mut pool, &b, &mut result, &the_last);
            }

            // OCCT L831-846: make the solid.
            let mut solid = pool.add_tsolid(vec![result.clone()]);
            set_closed_flag(&mut pool, &result, true);

            let mut sc = crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::new();
            sc.load(&solid);
            sc.perform_infinite_point(TOL_CONFUSION);
            if sc.state() == TOPABS_IN {
                solid = pool.add_tsolid(Vec::new());
                let a_local_shape = shape_reversed(&result);
                add_to_solid(&mut pool, &b, &mut solid, &a_local_shape);
            }
            update_map(the_original_s, &solid, &mut self.my_gen_map);
            solid
        } else {
            result
        }
    }

    /// OCCT BRepFill_Pipe::FindEdge(S, E, InitialLength) (cxx L861-911).
    fn find_edge(&self, s: &Shape, e: &Shape, initial_length: &mut i32) -> i32 {
        let mut result = 0;

        match s.shape_type() {
            ShapeType::Edge => {
                // OCCT L869-874.
                *initial_length += 1;
                if s.is_same(e) {
                    result = *initial_length;
                }
            }
            ShapeType::Wire => {
                // OCCT L877-889.
                let section = BRepFillShapeLaw::new_with_wire(s, false);
                let nb_law = section.nb_law();
                let mut ii = 1;
                while ii <= nb_law && result == 0 {
                    if e.is_same(&section.edge(ii)) {
                        result = *initial_length + ii as i32;
                    }
                    ii += 1;
                }
                *initial_length += nb_law as i32;
            }
            ShapeType::Face | ShapeType::Shell | ShapeType::Compound => {
                // OCCT L892-900.
                for it in topods_iterator(s) {
                    if result != 0 {
                        break;
                    }
                    result = self.find_edge(&it, e, initial_length);
                }
            }
            ShapeType::Solid | ShapeType::CompSolid => {
                // OCCT L904.
                panic!("Standard_DomainError: BRepFill_Pipe::SOLID or COMPSOLID");
            }
            _ => {}
        }

        result
    }

    /// OCCT BRepFill_Pipe::FindVertex(S, V, InitialLength) (cxx L918-991).
    fn find_vertex(&self, s: &Shape, v: &Shape, initial_length: &mut i32) -> i32 {
        let mut result = 0;

        match s.shape_type() {
            ShapeType::Vertex => {
                // OCCT L926-933.
                *initial_length += 1;
                if s.is_same(v) {
                    result = *initial_length;
                }
            }
            ShapeType::Edge => {
                // OCCT L935-954.
                let (mut v_f, mut v_l) = top_exp_vertices(s);
                if s.orientation == Orientation::Reversed {
                    std::mem::swap(&mut v_f, &mut v_l);
                }
                if v_f.is_same(v) {
                    result = *initial_length + 1;
                } else if v_l.is_same(v) {
                    result = *initial_length + 2;
                }
                *initial_length += 2;
            }
            ShapeType::Wire => {
                // OCCT L957-969: ii = InitialLength + 1 is captured before
                // the InitialLength increment.
                let mut ii = *initial_length + 1;
                let section = BRepFillShapeLaw::new_with_wire(s, false);
                *initial_length += section.nb_law() as i32 + 1;

                while ii <= *initial_length && result == 0 {
                    if v.is_same(&section.vertex(ii as usize, 0.0)) {
                        result = ii;
                    }
                    ii += 1;
                }
            }
            ShapeType::Face | ShapeType::Shell | ShapeType::Compound => {
                // OCCT L972-979.
                for it in topods_iterator(s) {
                    if result != 0 {
                        break;
                    }
                    result = self.find_vertex(&it, v, initial_length);
                }
            }
            ShapeType::Solid | ShapeType::CompSolid => {
                // OCCT L984.
                panic!("Standard_DomainError: BRepFill_Pipe::SOLID or COMPSOLID");
            }
            _ => {}
        }

        result
    }

    /// OCCT BRepFill_Pipe::DefineRealSegmax() (cxx L999-1050).
    fn define_real_segmax(&mut self) {
        let mut real_segmax = 0;

        // OCCT L1003-1044.
        for it in topods_iterator(&self.my_spine.clone()) {
            let e = it;
            let Some((mut c, first, last)) = brep_tool_curve(&e) else {
                continue;
            };
            // OCCT L1013-1024: strip the trimmed / offset wrappers.
            loop {
                match &c {
                    Curve3::Trimmed(t) => {
                        let TrimmedCurve3 { curve, .. } = t;
                        c = (**curve).clone();
                    }
                    Curve3::Offset(o) => {
                        c = (*o.basis).clone();
                    }
                    _ => break,
                }
            }
            // OCCT L1025-1043: the bspline spine case.
            if let Curve3::BSpline(bc) = &c {
                let nb_knots = bc.knots_mults().0.len() as i32;
                let mut real_nb_knots = nb_knots;
                if first > bc.first_parameter() {
                    let i1 = locate_u(bc, first);
                    real_nb_knots -= i1 - 1;
                }
                if last < bc.last_parameter() {
                    let i2 = locate_u(bc, last);
                    real_nb_knots -= nb_knots - i2;
                }
                real_segmax += real_nb_knots - 1;
            }
        }

        // OCCT L1046-1049.
        if self.my_segmax < real_segmax {
            self.my_segmax = real_segmax;
        }
    }

    /// OCCT BRepFill_Pipe::RebuildTopOrBottomFace(aFace, IsTop)
    /// (cxx L1058-1093) — GAP: UpdateTolFromTopOrBottomPCurve requires
    /// Adaptor3d_CurveOnSurface (not translated); the walk keeps the OCCT
    /// form up to the GAP point.
    fn rebuild_top_or_bottom_face(&self, _a_face: &Shape, is_top: bool) {
        // OCCT L1060.
        let my_sections = self.my_sections.as_ref().expect("mySections");
        let index_of_section = if is_top {
            1
        } else {
            my_sections.row_length()
        };
        let _ = index_of_section;
        let _ = self.my_cur_index_of_section_edge;

        // OCCT L1064-1092: the wire/edge walk + BB.Remove/UpdateTol/Add.
        update_tol_from_top_or_bottom_pcurve(_a_face, _a_face);
    }

    /// OCCT BRepFill_Pipe::BuildHistory(theSweep, theSection)
    /// (cxx L1100-1221).
    fn build_history(&mut self, the_sweep: &BRepFillSweep, the_section: &Shape) {
        // OCCT L1103: the U edges (intermediate faces).
        let an_u_edges = the_sweep.interfaces().expect("InterFaces");

        // OCCT L1108-1176: the section (wire) branch.
        if the_section.shape_type() == ShapeType::Wire {
            let a_section = the_section.clone();
            let wexp_sec = brep_tools_wire_explorer(&a_section);
            let mut inde = 0;
            for an_edge in &wexp_sec {
                inde += 1;
                // OCCT L1116-1119.
                if brep_tool_degenerated(an_edge) {
                    continue;
                }
                if self.my_gen_map.contains_key(&shape_key(an_edge)) {
                    continue;
                }

                // OCCT L1125-1126.
                let a_vertex = top_exp_vertices(an_edge);
                let a_vertex = [a_vertex.0, a_vertex.1];

                // OCCT L1131.
                let a_tape = the_sweep.tape(inde);

                // OCCT L1136-1145.
                let mut u_index = [inde, inde + 1];
                if an_edge.orientation == Orientation::Reversed {
                    u_index.swap(0, 1);
                }

                // OCCT L1147-1166.
                for kk in 0..2 {
                    if self.my_gen_map.contains_key(&shape_key(&a_vertex[kk])) {
                        continue;
                    }
                    let elist = self
                        .my_gen_map
                        .entry(shape_key(&a_vertex[kk]))
                        .or_default();
                    for jj in 1..=an_u_edges.upper_col() {
                        let an_uedge = an_u_edges.value(u_index[kk], jj);
                        if !an_uedge.is_null() {
                            elist.push(an_uedge);
                        }
                    }
                }

                // OCCT L1168-1174.
                let flist = self.my_gen_map.entry(shape_key(an_edge)).or_default();
                for itsh in topods_iterator(&a_tape) {
                    flist.push(itsh);
                }
            }
        }

        // OCCT L1179-1180: the spine subshapes.
        let a_faces = the_sweep.sub_shape().expect("SubShape");
        let a_v_edges = the_sweep.sections().expect("Sections");

        // OCCT L1182-1220.
        let wexp = brep_tools_wire_explorer(&self.my_spine.clone());
        let mut inde = 0;
        loop {
            let mut to_exit = false;
            if inde >= wexp.len() {
                to_exit = true;
            }

            inde += 1;

            if !to_exit {
                let an_edge_of_spine = &wexp[inde - 1];
                for i in 1..=a_faces.upper_row() {
                    let a_face = a_faces.value(i, inde);
                    update_map(an_edge_of_spine, &a_face, &mut self.my_gen_map);
                }
            }

            let a_vertex_of_spine = wire_explorer_current_vertex(&wexp, (inde - 1).min(wexp.len().saturating_sub(1)));
            for i in 1..=a_v_edges.upper_row() {
                let a_vedge = a_v_edges.value(i, inde);
                update_map(&a_vertex_of_spine, &a_vedge, &mut self.my_gen_map);
            }

            if to_exit {
                break;
            }
        }
    }
}

/// OCCT W.Closed(flag) / result.Closed(flag) — set the CLOSED flag of the
/// shape TShape (the offset_wire_b.rs reduction).
pub(crate) fn set_closed_flag(pool: &mut BRep, s: &Shape, closed: bool) {
    let idx = s.index;
    if closed {
        pool.set_flag(s.clone(), tshape_flags::CLOSED, true);
        let _ = idx;
    } else {
        pool.set_flag(s.clone(), tshape_flags::CLOSED, false);
    }
}

/// OCCT B.Add(result, S) for a shell / compound / compsolid / solid
/// result — dispatch on the result type.
fn add_to_result(pool: &mut BRep, _b: &BRepBuilder, result: &mut Shape, s: &Shape) {
    match result.shape_type() {
        ShapeType::Shell => add_to_shell(pool, _b, result, s),
        ShapeType::Solid => add_to_solid(pool, _b, result, s),
        ShapeType::Compound => add_to_compound(pool, _b, result, s),
        ShapeType::CompSolid => add_to_compsolid(pool, _b, result, s),
        ShapeType::Wire => add_to_wire(pool, result, s),
        _ => {}
    }
}

/// OCCT B.Add(Shell, Face).
pub(crate) fn add_to_shell(pool: &mut BRep, _b: &BRepBuilder, shell: &mut Shape, f: &Shape) {
    let idx = shell.index;
    if let TShape::Shell(sd) = Arc::make_mut(&mut pool.tshapes[idx]) {
        sd.faces.push(f.clone());
    }
}

/// OCCT B.Add(Solid, Shell).
pub(crate) fn add_to_solid(pool: &mut BRep, _b: &BRepBuilder, solid: &mut Shape, sh: &Shape) {
    let idx = solid.index;
    if let TShape::Solid(sd) = Arc::make_mut(&mut pool.tshapes[idx]) {
        sd.shells.push(sh.clone());
    }
}

/// OCCT B.Add(Compound, S).
fn add_to_compound(pool: &mut BRep, _b: &BRepBuilder, comp: &mut Shape, s: &Shape) {
    let idx = comp.index;
    if let TShape::Compound(children) = Arc::make_mut(&mut pool.tshapes[idx]) {
        children.push(s.clone());
    }
}

/// OCCT B.Add(CompSolid, Solid).
fn add_to_compsolid(pool: &mut BRep, _b: &BRepBuilder, cs: &mut Shape, s: &Shape) {
    let idx = cs.index;
    if let TShape::CompSolid(children) = Arc::make_mut(&mut pool.tshapes[idx]) {
        children.push(s.clone());
    }
}

/// OCCT B.Add(Wire, E).
fn add_to_wire(pool: &mut BRep, w: &Shape, e: &Shape) {
    let idx = w.index;
    if let TShape::Wire(wd) = Arc::make_mut(&mut pool.tshapes[idx]) {
        wd.edges.push(e.clone());
        wd.my_shapes.push(e.clone());
    }
}

/// OCCT Geom_BSplineCurve::LocateU(U, Tol, I1, I2) — the knot interval
/// index of the parameter (the bspline_ops reduction).
fn locate_u(bc: &rcad_kernel::geom::BSplineCurve3, u: f64) -> i32 {
    let (knots, _mults) = bc.knots_mults();
    let mut i = 1i32;
    for k in knots.iter() {
        if u < *k {
            break;
        }
        i += 1;
    }
    i
}
