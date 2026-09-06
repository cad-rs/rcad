// OCCT HLRBRep_Data (TKHLR/HLRBRep/HLRBRep_Data.hxx L1-289 + .cxx L1-2683
// + .lxx L1-120) — the HLR data structure: per-edge/per-face records, the
// bound-sort exploration state and the classification machinery.
// The impl blocks live in the sibling modules [update] (the ctor + the
// Write/Update/exploration group) and [classify] (the interference /
// classification group); [tableau_rejection] is the file-static
// TableauRejection class of the .cxx (L77-492).

use crate::hlr::algo::edges_block::MinMaxIndices;
use crate::hlr::algo::interference::Interference;
use crate::hlr::algo::projector::Projector;
use crate::hlr::brep::cl_props::CLProps;
use crate::hlr::brep::curve::Curve;
use crate::hlr::brep::edge_data::EdgeData;
use crate::hlr::brep::face_data::FaceData;
use crate::hlr::brep::face_iterator::FaceIterator;
use crate::hlr::brep::intersector::Intersector;
use crate::hlr::brep::sl_props::SLProps;
use crate::hlr::brep::surface::Surface;
use crate::topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool;

use rcad_kernel::topo::topods::Orientation;
use crate::bop::int_tools::bean_face_intersector::GeomAbsCurveType;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;

mod classify;
mod tableau_rejection;
mod update;

pub(crate) use classify::reject1;

/// OCCT file statics of HLRBRep_Data.cxx (L40-75): the TRACE counters are
/// file-scope globals in OCCT; rcad keeps them as thread-locals (the
/// Contap_HContTool static precedent).
pub mod counters {
    use std::cell::Cell;
    thread_local! {
        pub static NB_OK_INTERSECTION: Cell<i32> = const { Cell::new(0) };
        pub static NB_PT_INTERSECTION: Cell<i32> = const { Cell::new(0) };
        pub static NB_SEG_INTERSECTION: Cell<i32> = const { Cell::new(0) };
        pub static NB_CLASSIFICATION: Cell<i32> = const { Cell::new(0) };
        pub static NB_CAL1_INTERSECTION: Cell<i32> = const { Cell::new(0) };
        pub static NB_CAL2_INTERSECTION: Cell<i32> = const { Cell::new(0) };
        pub static NB_CAL3_INTERSECTION: Cell<i32> = const { Cell::new(0) };
    }
}

/// OCCT `static const double CutLar = 2.e-1;` (cxx L46).
pub const CUT_LAR: f64 = 2.0e-1;
/// OCCT `static const double CutBig = 1.e-1;` (cxx L47).
pub const CUT_BIG: f64 = 1.0e-1;
/// OCCT `static const double DERIVEE_PREMIERE_NULLE = 0.000000000001;`
/// (cxx L52).
pub const DERIVEE_PREMIERE_NULLE: f64 = 0.000_000_000_001;
/// OCCT `static long unsigned Mask32[32]` (cxx L58-62) — the bit masks.
pub const MASK32: [u64; 32] = [
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
    65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608,
    16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824,
    2147483648,
];
/// OCCT `static const int SIZEUV = 8;` (cxx L65).
pub const SIZEUV: usize = 8;

/// OCCT HLRBRep_Data (hxx L39-287) — the field block maps the OCCT private
/// members (L134-286) in declaration order.  The OCCT element pointers
/// (`iFaceData`, `myLEData`, ...) point into the myEData/myFData arrays;
/// they keep the OCCT raw-pointer form (the HLR methodological exception,
/// HLRBRep_Surface precedent) and are dereferenced in the impl modules.
pub struct Data<'a> {
    // hxx L135-137.
    my_nb_vertices: usize,
    my_nb_edges: usize,
    my_nb_faces: usize,
    // hxx L138-139 — NCollection_IndexedMap<TopoDS_Shape, ShapeMapHasher>
    // as the insertion-ordered Vec keyed by Shape::ptr_id() (the
    // HLRTopoBRep_Data precedent).
    my_e_map: Vec<rcad_kernel::topods::Shape>,
    my_f_map: Vec<rcad_kernel::topods::Shape>,
    // hxx L140-141 — NCollection_Array1 (1-based, index-shifted access).
    my_e_data: Vec<EdgeData<'a>>,
    my_f_data: Vec<FaceData<'a>>,
    // hxx L142.
    my_edge_indices: Vec<i32>,
    // hxx L143.
    my_toler: f32,
    // hxx L144.
    my_proj: Projector,
    // hxx L145-147.
    my_l_l_props: CLProps<'a>,
    my_f_l_props: CLProps<'a>,
    my_s_l_props: SLProps<'a>,
    // hxx L148.
    my_big_size: f64,
    // hxx L149-150.
    my_face_itr1: FaceIterator<'static>,
    my_face_itr2: FaceIterator<'static>,
    // hxx L151-157 — the current hiding-face exploration state.
    i_face: usize,
    i_face_data: *mut FaceData<'a>,
    i_face_geom: *mut Surface<'a>,
    i_face_min_max: *mut MinMaxIndices,
    i_face_type: GeomAbsSurfaceType,
    i_face_back: bool,
    i_face_simp: bool,
    i_face_smpl: bool,
    i_face_test: bool,
    // hxx L158.
    my_hide_count: i32,
    // hxx L159-160.
    my_deca: [f64; 16],
    my_sur_d: [f64; 16],
    // hxx L161-166 — the bound-sort edge exploration state.
    my_cur_sort_ed: usize,
    my_nbr_sort_ed: usize,
    my_le: usize,
    my_le_out_line: bool,
    my_le_internal: bool,
    my_le_double: bool,
    my_le_iso_line: bool,
    // hxx L167-171.
    my_le_data: *mut EdgeData<'a>,
    my_le_geom: *const Curve<'a>,
    my_le_min_max: *mut MinMaxIndices,
    my_le_type: GeomAbsCurveType,
    my_le_tol: f32,
    // hxx L172-180 — the current face-edge state.
    my_fe: usize,
    my_fe_ori: Orientation,
    my_fe_out_line: bool,
    my_fe_internal: bool,
    my_fe_double: bool,
    my_fe_data: *mut EdgeData<'a>,
    my_fe_geom: *mut Curve<'a>,
    my_fe_type: GeomAbsCurveType,
    my_fe_tol: f32,
    // hxx L181-182.
    my_intersector: Intersector<'a>,
    my_classifier: BRepTopolTool,
    // hxx L183-188.
    my_same_vertex: bool,
    my_intersected: bool,
    my_nb_points: usize,
    my_nb_segments: usize,
    i_interf: usize,
    my_intf: Interference,
    my_above_intf: bool,
    // hxx L189 — the TableauRejection instance (OCCT new/delete in the
    // ctor/Destroy; the unique-ownership raw pointer keeps the OCCT form).
    my_reject: *mut tableau_rejection::TableauRejection,
    // rcad kernel-context deviation (no hxx counterpart; the DSFiller
    // insert(brep, ...) precedent): OCCT reaches the kernel through the
    // global TShape graph (myFEGeom->Curve().Edge() etc.); the rcad Data
    // carries the owning BRep explicitly and the per-edge TopoDS shapes
    // live in my_e_map (aligned with my_e_data by construction index).
    // Ruling of the K/M contract clash (session 10): uv_point keeps the
    // 7-arg landed signature; M's three call sites read the brep from this
    // field and the edge Shape from my_e_map.
    my_brep: Option<&'a rcad_kernel::topods::BRep>,
}
