//! OCCT Draft_Modification (Draft_Modification.hxx L47-216 +
//! Draft_Modification.cxx L56-510) — the BRepTools_Modification driver of
//! the draft angle.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_Modification.hxx / Draft_Modification.cxx
//!
//! OCCT inheritance chain (hxx L47): Draft_Modification -> BRepTools_Modification.
//! The rcad BRepTools_Modification interface (TKTopAlgo/BRepTools) is not a
//! trait yet — the six override methods keep the OCCT virtual surface as
//! plain methods here (the BRepOffset_SimpleOffset GAP-carrier precedent in
//! offset/brep_offset_make_simple_offset.rs).  First consumer:
//! BRepOffsetAPI_DraftAngle (Stage 2e).
//!
//! Architecture differences:
//! 1. NCollection_IndexedDataMap<TopoDS_Face/Edge/Vertex, Info,
//!    TopTools_ShapeMapHasher> -> ShapeIndexedMap<V> (below): an insertion-
//!    ordered IndexMap keyed by the shape identity (TShape pointer +
//!    location; the TopTools_ShapeMapHasher key ignores orientation), with
//!    the OCCT 1-based FindKey/ChangeFromIndex accessors.
//! 2. NCollection_List<TopoDS_Shape> conneF -> Vec<Shape>.
//! 3. The OCCT methods mutate the global TShape arena through BRep_Builder
//!    (NewCurve2d: B.Range(NewE, ...), BRepTools::EvalAndUpdateTol:
//!    B.UpdateEdge(theE, Tol)); the rcad carriers mutate the Shape handle in
//!    place (&mut Shape / Arc::make_mut — the brep_algo::tool builder_*
//!    precedent).
//! 4. TopLoc_Location -> u32 id (0 = identity); Geom_Surface::Transformed /
//!    the location re-application are identity in the rcad flows (the
//!    loc_ope_split_drafts.rs #12 precedent).

use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::base::bnd_lib::curve2d_bounding_box;
use rcad_kernel::base::geom_proj_lib;
use rcad_kernel::geom::{
    translate_curve2d, Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval,
    TrimmedCurve3,
};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{CurveRepresentation, GeomAbsShape, ShapeType, TShape};
use std::sync::Arc;

pub(crate) use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_curve_on_surface, brep_tool_parameter, brep_tool_range,
    brep_tool_surface, brep_tool_tolerance, builder_range_edge, explorer, top_exp_vertices_raw,
};
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;

use super::draft_edge_info::DraftEdgeInfo;
use super::draft_error_status::DraftErrorStatus;
use super::draft_face_info::DraftFaceInfo;
use super::draft_vertex_info::DraftVertexInfo;

// ---------------------------------------------------------------------------
// Shape-keyed indexed map (architecture difference #1).
// ---------------------------------------------------------------------------

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
pub(crate) fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT NCollection_IndexedDataMap<K=TopoDS_Shape, V, TopTools_ShapeMapHasher>
/// — an insertion-ordered (Key, Value) table with the OCCT 1-based indexed
/// accessors (architecture difference #1).
#[derive(Debug, Clone)]
pub(crate) struct ShapeIndexedMap<V> {
    entries: IndexMap<(u64, u32), (Shape, V)>,
}

impl<V> Default for ShapeIndexedMap<V> {
    fn default() -> Self {
        ShapeIndexedMap {
            entries: IndexMap::new(),
        }
    }
}

impl<V> ShapeIndexedMap<V> {
    /// OCCT NCollection_IndexedDataMap::Clear().
    pub(crate) fn clear(&mut self) {
        self.entries.clear();
    }

    /// OCCT NCollection_IndexedDataMap::Extent().
    pub(crate) fn extent(&self) -> usize {
        self.entries.len()
    }

    /// OCCT NCollection_IndexedDataMap::Contains(K).
    pub(crate) fn contains(&self, the_key: &Shape) -> bool {
        self.entries.contains_key(&shape_key(the_key))
    }

    /// OCCT NCollection_IndexedDataMap::Add(K, V) — the OCCT map raises on a
    /// duplicate key; the rcad insert asserts the same invariant.
    pub(crate) fn add(&mut self, the_key: &Shape, the_value: V) {
        assert!(
            !self.contains(the_key),
            "NCollection_IndexedDataMap::Add: duplicate key"
        );
        self.entries
            .insert(shape_key(the_key), (the_key.clone(), the_value));
    }

    /// OCCT NCollection_IndexedDataMap::FindKey(I) — the 1-based key.
    pub(crate) fn find_key(&self, the_i: usize) -> &Shape {
        &self.entries.get_index(the_i - 1).expect("FindKey: out of range").1 .0
    }

    /// OCCT NCollection_IndexedDataMap::FindFromIndex(I) — the 1-based value.
    #[allow(dead_code)]
    pub(crate) fn find_from_index(&self, the_i: usize) -> &V {
        &self
            .entries
            .get_index(the_i - 1)
            .expect("FindFromIndex: out of range")
            .1
             .1
    }

    /// OCCT NCollection_IndexedDataMap::FindFromKey(K).
    pub(crate) fn find_from_key(&self, the_key: &Shape) -> &V {
        &self
            .entries
            .get(&shape_key(the_key))
            .expect("NCollection_IndexedDataMap::FindFromKey: no such key")
            .1
    }

    /// OCCT NCollection_IndexedDataMap::ChangeFromKey(K).
    pub(crate) fn change_from_key(&mut self, the_key: &Shape) -> &mut V {
        &mut self
            .entries
            .get_mut(&shape_key(the_key))
            .expect("NCollection_IndexedDataMap::ChangeFromKey: no such key")
            .1
    }

    /// OCCT NCollection_IndexedDataMap::ChangeFromIndex(I) — the 1-based
    /// mutable value.
    pub(crate) fn change_from_index(&mut self, the_i: usize) -> &mut V {
        &mut self
            .entries
            .get_index_mut(the_i - 1)
            .expect("ChangeFromIndex: out of range")
            .1
             .1
    }

    /// OCCT NCollection_IndexedDataMap::RemoveKey(K).
    pub(crate) fn remove_key(&mut self, the_key: &Shape) {
        self.entries
            .swap_remove(&shape_key(the_key));
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool / BRepTools re-hosts (module-shared; the loc_ope precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Pnt(vtx) — the vertex point (None carries the OCCT null
/// TShape raise).
pub(crate) fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) (BRep_Tool.cxx L1180-1188 -> the
/// representation scan) — the regularity stored on the (E, F1/F2) curve
/// representations; the default is GeomAbs_C0 when nothing matches
/// (offset/brep_offset_make_simple_offset.rs re-host copy; the rcad face
/// locations are identity in this pipeline).
pub(crate) fn brep_tool_continuity(the_e: &Shape, the_f1: &Shape, the_f2: &Shape) -> GeomAbsShape {
    // OCCT: S1 = Surface(F1); S2 = Surface(F2) — the local surfaces.
    let (Some(s1), Some(s2)) = (brep_tool_surface(the_f1), brep_tool_surface(the_f2)) else {
        // The OCCT null-surface case has no matching representation.
        return GeomAbsShape::C0;
    };
    let Some(ed) = the_e.as_edge() else {
        return GeomAbsShape::C0;
    };
    for cr in &ed.representations {
        // OCCT: cr->IsRegularity(S1, S2, l1, l2) — the rcad face locations
        // are identity in this pipeline.
        if cr.is_regularity_on(&s1, &s2, 0, 0) {
            // OCCT: return cr->Continuity();
            if let CurveRepresentation::CurveOn2Surfaces { continuity, .. } = cr {
                return *continuity;
            }
        }
    }
    // OCCT: return GeomAbs_C0;
    GeomAbsShape::C0
}

/// OCCT BRep_Tool::IsClosed(E, F) — the pcurve on F is a closed-surface
/// representation (BRep_Tool.cxx; the brep_algo::tool re-host).
pub(crate) fn brep_tool_is_closed(the_e: &Shape, the_f: &Shape) -> bool {
    crate::brep_algo::tool::brep_tool_is_closed_on_surface(the_e, the_f)
}

/// OCCT BRepTools::IsReallyClosed(E, F) (BRepTools.cxx L1204-1222) — the edge
/// occurs exactly twice in the face wire structure.
pub(crate) fn brep_tools_is_really_closed(the_e: &Shape, the_f: &Shape) -> bool {
    // OCCT: if (!BRep_Tool::IsClosed(E, F)) return false;
    if !brep_tool_is_closed(the_e, the_f) {
        return false;
    }
    // OCCT: nbocc = 0; for edges of F: if IsSame(E) nbocc++;
    let mut nbocc = 0;
    for a_cur in explorer(the_f, ShapeType::Edge, ShapeType::Shape) {
        if a_cur.is_same(the_e) {
            nbocc += 1;
        }
    }
    // OCCT: return nbocc == 2;
    nbocc == 2
}

/// OCCT BRepTools::UVBounds(F, Umin, Umax, Vmin, Vmax) (BRepTools.cxx
/// L137-181 via AddUVBounds/BndLib_Add2dCurve) — the union of the 2D boxes
/// of the face's edge pcurves; returns [umin, umax, vmin, vmax] (the
/// loc_ope_wires_on_shape_b.rs #14 re-host form).
pub(crate) fn brep_tools_uv_bounds(the_f: &Shape) -> [f64; 4] {
    // OCCT: Umin = Vmin = RealLast(); Umax = Vmax = -RealLast();
    let (mut umin, mut umax, mut vmin, mut vmax) =
        (f64::MAX, -f64::MAX, f64::MAX, -f64::MAX);
    for a_edg in explorer(the_f, ShapeType::Edge, ShapeType::Shape) {
        if let Some((a_c2d, a_f, a_l)) = brep_tool_curve_on_surface(&a_edg, the_f) {
            // OCCT: BndLib_Add2dCurve::Add(C2d, f, l, Tol, B); B.Get(...).
            let b = curve2d_bounding_box(&a_c2d, a_f, a_l, 0.0);
            umin = umin.min(b[0]);
            umax = umax.max(b[1]);
            vmin = vmin.min(b[2]);
            vmax = vmax.max(b[3]);
        }
    }
    [umin, umax, vmin, vmax]
}

/// OCCT BRepTools::EvalAndUpdateTol(theE, C3d, C2d, S, f, l)
/// (BRepTools.cxx L1255-1344) — the max 3D/2D deviation updates the edge
/// tolerance (architecture difference #3: the arena UpdateEdge mutates the
/// handle in place).
pub(crate) fn brep_tools_eval_and_update_tol(
    the_e: &mut Shape,
    the_c3d: &Curve3,
    the_c2d: &Curve2d,
    the_s: &Surface3,
    f: f64,
    l: f64,
) -> f64 {
    let mut newtol = 0.0;
    let (mut first, mut last) = (f, l);
    // OCCT L1264-1274: set first/last to avoid ErrorStatus 2 in
    // CheckCurveOnSurface; the periodic curves keep the call range.
    let dom3 = the_c3d.default_domain();
    let dom2 = the_c2d.default_domain();
    let per3 = the_c3d.is_periodic();
    let per2 = the_c2d.is_periodic();
    if !per3 {
        first = first.max(dom3[0]);
        last = last.min(dom3[1]);
    }
    if !per2 {
        first = first.max(dom2[0]);
        last = last.min(dom2[1]);
    }

    // OCCT L1276-1283: GeomLib_CheckCurveOnSurface CT(...); CT.Perform(...)
    // — GAP: the TKTopAlgo/GeomLib checker is not translated; the carrier
    // stays not-done with the OCCT ErrorStatus 2 so the tail takes the OCCT
    // sampling branch for periodic curves and keeps newtol = 0 otherwise.
    let ct_is_done = false;
    let ct_error_status = 2;
    if ct_is_done {
        // OCCT L1285-1287: newtol = CT.MaxDistance(); — unreachable until the
        // GeomLib_CheckCurveOnSurface translation lands.
    } else if ct_error_status == 3 || (ct_error_status == 2 && (per3 || per2)) {
        // OCCT L1291-1332: the by-sample estimate.
        let nbint = 22;
        let mut dt = (last - first) / nbint as f64;
        dt = dt.max(CONFUSION);
        let mut dmax = 0.0f64;
        let mut cnt = 0;
        let mut t = first;
        while t <= last {
            cnt += 1;
            let a_p2d = the_c2d.point_at(t);
            let a_pc = the_c3d.point_at(t);
            let a_ps = the_s.point_at(a_p2d.x, a_p2d.y);
            let d = a_ps.distance_squared(a_pc);
            if d > dmax {
                dmax = d;
            }
            t += dt;
        }
        if cnt < nbint + 1 {
            let a_p2d = the_c2d.point_at(last);
            let a_pc = the_c3d.point_at(last);
            let a_ps = the_s.point_at(a_p2d.x, a_p2d.y);
            let d = a_ps.distance_squared(a_pc);
            if d > dmax {
                dmax = d;
            }
        }
        newtol = 1.2 * dmax.sqrt();
    }
    // OCCT L1333-1343: Tol = BRep_Tool::Tolerance(theE);
    // if (newtol > Tol) { Tol = newtol; B.UpdateEdge(theE, Tol); }
    let mut tol = brep_tool_tolerance(the_e);
    if newtol > tol {
        tol = newtol;
        if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
            ed.tolerance = tol;
        }
    }
    // OCCT L1345: return Tol;
    tol
}

// ---------------------------------------------------------------------------
// The class (Draft_Modification.hxx L47-216).
// ---------------------------------------------------------------------------

/// OCCT Draft_Modification (Draft_Modification.hxx L47-216).
pub struct DraftModification {
    my_fmap: ShapeIndexedMap<DraftFaceInfo>, // OCCT: myFMap
    my_emap: ShapeIndexedMap<DraftEdgeInfo>, // OCCT: myEMap
    my_vmap: ShapeIndexedMap<DraftVertexInfo>, // OCCT: myVMap
    my_comp: bool,                           // OCCT: myComp
    my_shape: Shape,                         // OCCT: myShape
    bad_shape: Shape,                        // OCCT: badShape
    err_stat: DraftErrorStatus,              // OCCT: errStat
    cur_face: Shape,                         // OCCT: curFace
    conne_f: Vec<Shape>,                     // OCCT: conneF
    // OCCT: myEFMap — NCollection_IndexedDataMap<TopoDS_Shape,
    // NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>.
    my_efmap: IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
}

impl DraftModification {
    /// OCCT Draft_Modification::Draft_Modification(const TopoDS_Shape& S)
    /// (cxx L56-61).
    pub fn new(the_s: &Shape) -> Self {
        let mut m = DraftModification {
            my_fmap: ShapeIndexedMap::default(),
            my_emap: ShapeIndexedMap::default(),
            my_vmap: ShapeIndexedMap::default(),
            my_comp: false,
            my_shape: the_s.clone(),
            bad_shape: Shape::null(),
            err_stat: DraftErrorStatus::NoError,
            cur_face: Shape::null(),
            conne_f: Vec::new(),
            my_efmap: IndexMap::new(),
        };
        // OCCT L60: TopExp::MapShapesAndAncestors(myShape, EDGE, FACE, myEFMap);
        map_shapes_and_ancestors(
            &m.my_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut m.my_efmap,
        );
        m
    }

    /// OCCT Draft_Modification::Clear() (cxx L65-74) — resets on the same
    /// shape.
    pub fn clear(&mut self) {
        self.my_comp = false;
        self.my_fmap.clear();
        self.my_emap.clear();
        self.my_vmap.clear();
        self.my_efmap.clear();
        self.bad_shape = Shape::null();
        self.err_stat = DraftErrorStatus::NoError;
    }

    /// OCCT Draft_Modification::Init(const TopoDS_Shape& S) (cxx L78-83) —
    /// changes the basis shape and resets.
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.clear();
        // OCCT L82: TopExp::MapShapesAndAncestors(myShape, EDGE, FACE, myEFMap);
        map_shapes_and_ancestors(
            &self.my_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut self.my_efmap,
        );
    }

    /// OCCT Draft_Modification::Add(F, Direction, Angle, NeutralPlane, Flag)
    /// (cxx L87-104).  The body (InternalAdd) lives in
    /// offset/draft_modification_1.rs.
    pub fn add(
        &mut self,
        the_f: &Shape,
        the_direction: glam::DVec3,
        the_angle: f64,
        the_neutral_plane: &rcad_kernel::geom::Plane,
        the_flag: bool,
    ) -> bool {
        // OCCT L93-96: if (!badShape.IsNull()) throw Standard_ConstructionError();
        if !self.bad_shape.is_null() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L98-101: if (myComp) Clear();
        if self.my_comp {
            self.clear();
        }
        // OCCT L102: curFace = F;
        self.cur_face = the_f.clone();
        // OCCT L103: return InternalAdd(F, Direction, Angle, NeutralPlane, Flag);
        self.internal_add(the_f, the_direction, the_angle, the_neutral_plane, the_flag)
    }

    /// OCCT Draft_Modification::Remove(const TopoDS_Face& F) (cxx L108-154).
    pub fn remove(&mut self, the_f: &Shape) {
        // OCCT L110-113: if (!myFMap.Contains(F) || myComp)
        //                  throw Standard_NoSuchObject();
        if !self.my_fmap.contains(the_f) || self.my_comp {
            panic!("Standard_NoSuchObject");
        }
        // OCCT L115: conneF.Clear();
        self.conne_f.clear();
        // OCCT L116: NCollection_List<TopoDS_Shape>::Iterator ltod; — carried
        // as an index over conne_f.
        // OCCT L118: curFace = myFMap.FindFromKey(F).RootFace();
        self.cur_face = self.my_fmap.find_from_key(the_f).root_face().clone();
        // OCCT L119-130.
        for i in 1..=self.my_fmap.extent() {
            let the_f_i = self.my_fmap.find_key(i).clone();
            if self
                .my_fmap
                .find_from_key(&the_f_i)
                .root_face()
                .is_same(&self.cur_face)
            {
                self.conne_f.push(the_f_i.clone());
                if the_f_i.is_same(&self.bad_shape) {
                    self.bad_shape = Shape::null();
                }
            }
        }
        // OCCT L132-137: remove the collected faces.
        let mut ltod = 0usize;
        while ltod < self.conne_f.len() {
            let f = self.conne_f[ltod].clone();
            self.my_fmap.remove_key(&f); // TopoDS::Face(ltod.Value())
            ltod += 1;
        }
        // OCCT L139: conneF.Clear();
        self.conne_f.clear();
        // OCCT L140-147: collect the edges of the same root face.
        for i in 1..=self.my_emap.extent() {
            let the_e = self.my_emap.find_key(i).clone();
            if self
                .my_emap
                .find_from_key(&the_e)
                .root_face()
                .is_same(&self.cur_face)
            {
                self.conne_f.push(the_e);
            }
        }
        // OCCT L148-153: remove the collected edges.
        let mut ltod = 0usize;
        while ltod < self.conne_f.len() {
            let e = self.conne_f[ltod].clone();
            self.my_emap.remove_key(&e); // TopoDS::Edge(ltod.Value())
            ltod += 1;
        }
    }

    /// OCCT Draft_Modification::IsDone() (cxx L158-161).
    pub fn is_done(&self) -> bool {
        self.my_comp && self.bad_shape.is_null()
    }

    /// OCCT Draft_Modification::Error() (cxx L165-168).
    pub fn error(&self) -> DraftErrorStatus {
        self.err_stat
    }

    /// OCCT Draft_Modification::ProblematicShape() (cxx L172-175).
    pub fn problematic_shape(&self) -> &Shape {
        &self.bad_shape
    }

    /// OCCT Draft_Modification::ConnectedFaces(F) (cxx L179-202).
    pub fn connected_faces(&mut self, the_f: &Shape) -> &Vec<Shape> {
        // OCCT L181-184.
        if !self.my_fmap.contains(the_f) {
            panic!("Standard_NoSuchObject");
        }
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        // OCCT L189: conneF.Clear();
        self.conne_f.clear();
        // OCCT L190: curFace = myFMap.FindFromKey(F).RootFace();
        self.cur_face = self.my_fmap.find_from_key(the_f).root_face().clone();
        // OCCT L192-199.
        for i in 1..=self.my_fmap.extent() {
            let the_f_i = self.my_fmap.find_key(i).clone();
            if self
                .my_fmap
                .find_from_key(&the_f_i)
                .root_face()
                .is_same(&self.cur_face)
            {
                self.conne_f.push(the_f_i);
            }
        }
        // OCCT L201: return conneF;
        &self.conne_f
    }

    /// OCCT Draft_Modification::ModifiedFaces() (cxx L206-224).
    pub fn modified_faces(&mut self) -> &Vec<Shape> {
        // OCCT L208-211.
        if !self.bad_shape.is_null() {
            panic!("StdFail_NotDone");
        }
        self.conne_f.clear();
        // OCCT L214-221.
        for i in 1..=self.my_fmap.extent() {
            let the_f = self.my_fmap.find_key(i).clone();
            if !self.my_fmap.find_from_key(&the_f).root_face().is_null() {
                self.conne_f.push(the_f);
            }
        }
        // OCCT L223: return conneF;
        &self.conne_f
    }

    /// OCCT Draft_Modification::NewSurface(F, S, L, Tol, RevWires, RevFace)
    /// (cxx L228-256) — the BRepTools_Modification override.
    pub fn new_surface(
        &mut self,
        the_f: &Shape,
        the_s: &mut Option<Surface3>,
        the_l: &mut u32,
        the_tol: &mut f64,
        the_rev_wires: &mut bool,
        the_rev_face: &mut bool,
    ) -> bool {
        // OCCT L235-238.
        if !self.is_done() {
            panic!("Standard_DomainError");
        }
        // OCCT L240-243.
        if !self.my_fmap.contains(the_f) || !self.my_fmap.find_from_key(the_f).new_geometry() {
            return false;
        }
        // OCCT L245-247.
        *the_rev_wires = false;
        *the_rev_face = false;
        *the_tol = brep_tool_tolerance(the_f);
        // OCCT L249: S = BRep_Tool::Surface(F, L) — the face surface location
        // is carried as the u32 id.
        *the_l = match the_f.data.as_ref() {
            TShape::Face(fd) => fd.surface_location,
            _ => 0,
        };
        *the_s = brep_tool_surface(the_f);
        // OCCT L251: L.Identity();
        *the_l = 0;
        // OCCT L253: S = myFMap.FindFromKey(F).Geometry();
        *the_s = self.my_fmap.find_from_key(the_f).geometry().cloned();
        // OCCT L255: return true;
        true
    }

    /// OCCT Draft_Modification::NewCurve(E, C, L, Tol) (cxx L260-287) — the
    /// BRepTools_Modification override.
    pub fn new_curve(
        &mut self,
        the_e: &Shape,
        the_c: &mut Option<Curve3>,
        the_l: &mut u32,
        the_tol: &mut f64,
    ) -> bool {
        // OCCT L265-268.
        if !self.is_done() {
            panic!("Standard_DomainError");
        }
        // OCCT L270-273.
        if !self.my_emap.contains(the_e) {
            return false;
        }
        // OCCT L275-282.
        let einf_tol = {
            let einf = self.my_emap.find_from_key(the_e);
            if !einf.new_geometry() {
                return false;
            }
            einf.tolerance()
        };
        // OCCT L281-282: Tol = Einf.Tolerance(); Tol = max(Tol, BRep_Tool::Tolerance(E));
        *the_tol = einf_tol.max(brep_tool_tolerance(the_e));
        // OCCT L283: L.Identity();
        *the_l = 0;
        // OCCT L284: C = myEMap.FindFromKey(E).Geometry();
        *the_c = self.my_emap.find_from_key(the_e).geometry().cloned();
        // OCCT L286: return true;
        true
    }

    /// OCCT Draft_Modification::NewPoint(V, P, Tol) (cxx L291-306) — the
    /// BRepTools_Modification override.
    pub fn new_point(&mut self, the_v: &Shape, the_p: &mut glam::DVec3, the_tol: &mut f64) -> bool {
        // OCCT L293-296.
        if !self.is_done() {
            panic!("Standard_DomainError");
        }
        // OCCT L298-301.
        if !self.my_vmap.contains(the_v) {
            return false;
        }
        // OCCT L303: Tol = BRep_Tool::Tolerance(V);
        *the_tol = brep_tool_tolerance(the_v);
        // OCCT L304: P = myVMap.FindFromKey(V).Geometry();
        *the_p = *self.my_vmap.find_from_key(the_v).geometry();
        // OCCT L305: return true;
        true
    }

    /// OCCT Draft_Modification::NewCurve2d(E, F, NewE, NewF, C, Tol)
    /// (cxx L310-427) — the BRepTools_Modification override.  NewE is
    /// &mut (architecture difference #3: B.Range(NewE, Fp, Lp) and the
    /// tolerance update mutate the handle in place).
    pub fn new_curve2d(
        &mut self,
        the_e: &Shape,
        the_f: &Shape,
        the_new_e: &mut Shape,
        the_new_f: &Shape,
        the_c: &mut Option<Curve2d>,
        the_tol: &mut f64,
    ) -> bool {
        let _ = the_new_f; // OCCT L313: the NewF parameter is unnamed/unused.
        // OCCT L318-321.
        if !self.is_done() {
            panic!("Standard_DomainError");
        }
        // OCCT L323-326.
        if !self.my_emap.contains(the_e) {
            return false;
        }
        // OCCT L328-329: double Fp, Lp; BRep_Tool::Range(NewE, Fp, Lp);
        let (mut fp, mut lp) = brep_tool_range(the_new_e);
        // OCCT L331: SB = myFMap.FindFromKey(F).Geometry();
        let mut sb: Option<Surface3> = self.my_fmap.find_from_key(the_f).geometry().cloned();
        // OCCT L333: Tol = BRep_Tool::Tolerance(E);
        *the_tol = brep_tool_tolerance(the_e);
        // OCCT L335-343.
        let einf = self.my_emap.find_from_key(the_e);
        let mut assigned = false;
        if einf.first_face().is_same(the_f) && einf.first_pc().is_some() {
            *the_c = einf.first_pc().cloned();
            assigned = true;
        } else if einf.second_face().is_same(the_f) && einf.second_pc().is_some() {
            *the_c = einf.second_pc().cloned();
            assigned = true;
        }
        if !assigned {
            // OCCT L347-355.
            if !self.my_emap.find_from_key(the_e).new_geometry() {
                let (fpi, lpi) = brep_tool_range(the_e);
                if fpi <= fp && fp <= lpi && fpi <= lp && lp <= lpi {
                    return false;
                }
            }
            // OCCT L357-358: BRep_Tool::Range(NewE, Fp, Lp);
            let (nfp, nlp) = brep_tool_range(the_new_e);
            fp = nfp;
            lp = nlp;
            // OCCT L359-360: TC = new Geom_TrimmedCurve(Geometry(), Fp, Lp);
            let tc = Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(
                    self.my_emap
                        .find_from_key(the_e)
                        .geometry()
                        .cloned()
                        .expect("Geom_TrimmedCurve(null)"),
                ),
                first: fp,
                last: lp,
            });
            // OCCT L361-362: Fp = TC->FirstParameter(); Lp = TC->LastParameter();
            if let Curve3::Trimmed(t) = &tc {
                fp = t.first;
                lp = t.last;
            }
            // OCCT L363-364: BRep_Builder B; B.Range(NewE, Fp, Lp);
            builder_range_edge(the_new_e, fp, lp);
            // OCCT L365: C = GeomProjLib::Curve2d(TC, Fp, Lp, SB, Tol);
            if let Curve3::Trimmed(t) = &tc {
                *the_c = match &sb {
                    Some(s) => {
                        // The natural-bounds overload (the OCCT Tol is carried
                        // by the kernel default).
                        let dom = s.default_domain();
                        geom_proj_lib::curve2d(&t.curve, fp, lp, s, dom[0], dom[1], dom[2], dom[3])
                    }
                    None => None,
                };
            }
        }

        // OCCT L368-373: unwrap one rectangular trimmed level of SB.
        if let Some(s) = &sb {
            if let Surface3::Trimmed(t) = s {
                sb = Some(t.basis.as_ref().clone());
            }
        }

        // OCCT L375-388: JeRecadre.
        let mut je_recadre = false;
        if let Some(s) = &sb {
            if let Surface3::LinearExtrusion(sle) = s {
                // OCCT L376-384: the basis curve being a circle forces the
                // re-framing.
                let a_c = sle.profile.as_ref();
                if matches!(a_c, Curve3::Circle(_)) {
                    je_recadre = true;
                }
            }
            // OCCT L386-388.
            je_recadre = je_recadre
                || matches!(s, Surface3::Cylinder(_))
                || matches!(s, Surface3::Cone(_));
            // OCCT names Geom_SphericalSurface as well; the rcad match keeps
            // the same predicate (sphere surfaces arrive as Spherical).
            je_recadre = je_recadre || matches!(s, Surface3::Sphere(_));
        }

        // OCCT L390-422.
        if je_recadre {
            let mut b_translate;
            let a_d2;
            let vectra = glam::DVec2::new(2.0 * std::f64::consts::PI, 0.0);
            let a_v2dt;
            // OCCT L398: aC2DE = BRep_Tool::CurveOnSurface(E, F, aT1, aT2);
            let a_c2de = brep_tool_curve_on_surface(the_e, the_f)
                .expect("BRep_Tool::CurveOnSurface(E, F)");
            let (a_t1, a_t2) = (a_c2de.1, a_c2de.2);
            // OCCT L400: PF = aC2DE->Value(0.5 * (aT1 + aT2));
            let pf = a_c2de.0.point_at(0.5 * (a_t1 + a_t2));
            // OCCT L402: NewPF = C->Value(0.5 * (Fp + Lp));
            let new_pf = the_c
                .as_ref()
                .expect("null Geom2d_Curve")
                .point_at(0.5 * (fp + lp));
            // OCCT L404: aD2 = NewPF.SquareDistance(PF);
            a_d2 = (new_pf - pf).length_squared();
            // OCCT L406-416.
            b_translate = false;
            if (new_pf + vectra - pf).length_squared() < a_d2 {
                a_v2dt = vectra;
                b_translate = true; // OCCT: !bTranslate
            } else if (new_pf - vectra - pf).length_squared() < a_d2 {
                a_v2dt = -vectra;
                b_translate = true; // OCCT: !bTranslate
            } else {
                a_v2dt = glam::DVec2::ZERO;
            }
            // OCCT L418-421: if (bTranslate) C->Translate(aV2DT);
            if b_translate {
                if let Some(c) = the_c {
                    *c = translate_curve2d(c, a_v2dt);
                }
            }
            let _ = a_d2;
        }
        // OCCT L424-425: aC3d = BRep_Tool::Curve(NewE, Fp, Lp);
        // Tol = BRepTools::EvalAndUpdateTol(NewE, aC3d, C, SB, Fp, Lp);
        let (a_c3d, _, _) = brep_tool_curve(the_new_e)
            .expect("BRep_Tool::Curve(NewE)");
        let a_c3d = {
            // OCCT passes the located curve; the location re-application is
            // identity (architecture difference #4).
            a_c3d
        };
        *the_tol = match the_c.as_ref() {
            Some(c) => match &sb {
                Some(s) => brep_tools_eval_and_update_tol(the_new_e, &a_c3d, c, s, fp, lp),
                None => brep_tool_tolerance(the_new_e),
            },
            None => brep_tool_tolerance(the_new_e),
        };
        // OCCT L426: return true;
        true
    }

    /// OCCT Draft_Modification::NewParameter(V, E, P, Tol) (cxx L431-498) —
    /// the BRepTools_Modification override.
    pub fn new_parameter(
        &mut self,
        the_v: &Shape,
        the_e: &Shape,
        the_p: &mut f64,
        the_tol: &mut f64,
    ) -> bool {
        // OCCT L437-440.
        if !self.is_done() {
            panic!("Standard_DomainError");
        }
        // OCCT L442-445.
        if !self.my_vmap.contains(the_v) {
            return false;
        }
        // OCCT L447: P = myVMap.ChangeFromKey(V).Parameter(E);
        *the_p = self.my_vmap.change_from_key(the_v).parameter(the_e);
        // OCCT L448-454: GC = myEMap.FindFromKey(E).Geometry(); the OCCT
        // trimmed down_cast re-assigns the same handle (the type stays
        // Geom_TrimmedCurve) — the rcad form keeps the curve unchanged.
        let gc = self
            .my_emap
            .find_from_key(the_e)
            .geometry()
            .cloned()
            .expect("null Geom_Curve");

        // OCCT L456: if (GC->IsClosed()).
        if gc.is_closed() {
            // OCCT L458: FV = TopExp::FirstVertex(E);
            let mut fv = top_exp_vertices_raw(the_e)
                .0
                .expect("TopExp::FirstVertex");
            // OCCT L459-467.
            let paramf = if self.my_vmap.contains(&fv) {
                self.my_vmap.change_from_key(&fv).parameter(the_e)
            } else {
                brep_tool_parameter(&fv, the_e)
            };

            // OCCT L469-471: Patch — FirstPar/LastPar/pconf.
            let dom = gc.default_domain();
            let first_par = dom[0];
            let last_par = dom[1];
            let pconf = PCONFUSION;
            // OCCT L472-480.
            if (paramf - last_par).abs() <= pconf {
                let paramf = first_par;
                fv.orientation = the_e.orientation;
                if the_v.is_equal(&fv) {
                    *the_p = paramf;
                }
            }

            // OCCT L482-493.
            fv.orientation = the_e.orientation;
            if !the_v.is_equal(&fv) && *the_p <= paramf {
                if gc.is_periodic() {
                    // OCCT L487: P += GC->Period(); — the period is the domain
                    // length of the periodic curve.
                    *the_p += dom[1] - dom[0];
                } else {
                    *the_p = dom[1];
                }
            }
        }

        // OCCT L496: Tol = max(Tolerance(V), Tolerance(E));
        *the_tol = brep_tool_tolerance(the_v).max(brep_tool_tolerance(the_e));
        // OCCT L497: return true;
        true
    }

    /// OCCT Draft_Modification::Continuity(E, F1, F2, NewE, NewF1, NewF2)
    /// (cxx L502-510) — the BRepTools_Modification override.
    pub fn continuity(
        &mut self,
        the_e: &Shape,
        the_f1: &Shape,
        the_f2: &Shape,
        _the_new_e: &Shape,
        _the_new_f1: &Shape,
        _the_new_f2: &Shape,
    ) -> GeomAbsShape {
        // OCCT L509: return BRep_Tool::Continuity(E, F1, F2);
        brep_tool_continuity(the_e, the_f1, the_f2)
    }

    // -----------------------------------------------------------------------
    // Accessors used by the _1 submodule (Draft_Modification.hxx L205-215).
    // -----------------------------------------------------------------------

    /// OCCT: myFMap.
    pub(crate) fn fmap(&self) -> &ShapeIndexedMap<DraftFaceInfo> {
        &self.my_fmap
    }

    /// Mutable OCCT: myFMap.
    pub(crate) fn fmap_mut(&mut self) -> &mut ShapeIndexedMap<DraftFaceInfo> {
        &mut self.my_fmap
    }

    /// OCCT: myEMap.
    pub(crate) fn emap(&self) -> &ShapeIndexedMap<DraftEdgeInfo> {
        &self.my_emap
    }

    /// Mutable OCCT: myEMap.
    pub(crate) fn emap_mut(&mut self) -> &mut ShapeIndexedMap<DraftEdgeInfo> {
        &mut self.my_emap
    }

    /// OCCT: myVMap.
    pub(crate) fn vmap(&self) -> &ShapeIndexedMap<DraftVertexInfo> {
        &self.my_vmap
    }

    /// Mutable OCCT: myVMap.
    pub(crate) fn vmap_mut(&mut self) -> &mut ShapeIndexedMap<DraftVertexInfo> {
        &mut self.my_vmap
    }

    /// OCCT: myComp.
    pub(crate) fn comp(&self) -> bool {
        self.my_comp
    }

    /// Mutable OCCT: myComp.
    pub(crate) fn set_comp(&mut self, v: bool) {
        self.my_comp = v;
    }

    /// OCCT: myShape.
    pub(crate) fn my_shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT: badShape.
    pub(crate) fn bad_shape(&self) -> &Shape {
        &self.bad_shape
    }

    /// Mutable OCCT: badShape.
    pub(crate) fn bad_shape_mut(&mut self) -> &mut Shape {
        &mut self.bad_shape
    }

    /// Mutable OCCT: errStat.
    pub(crate) fn set_err_stat(&mut self, v: DraftErrorStatus) {
        self.err_stat = v;
    }

    /// OCCT: curFace.
    pub(crate) fn cur_face(&self) -> &Shape {
        &self.cur_face
    }

    /// The disjoint (myFMap, myEMap, myVMap) borrows consumed by the Choose
    /// call of Perform (the OCCT call passes three separate objects; the rcad
    /// maps are members of one struct, so the split borrow is explicit).
    pub(crate) fn draft_maps_split(
        &mut self,
    ) -> (
        &ShapeIndexedMap<DraftFaceInfo>,
        &mut ShapeIndexedMap<DraftEdgeInfo>,
        &mut ShapeIndexedMap<DraftVertexInfo>,
    ) {
        (&self.my_fmap, &mut self.my_emap, &mut self.my_vmap)
    }

    /// OCCT: myEFMap — FindFromKey(edg) (the OCCT map raises on an absent
    /// key; the rcad get returns an empty list at the single L520 call site
    /// where the OCCT key is guaranteed present by the myEMap construction).
    pub(crate) fn efmap_find(&self, the_e: &Shape) -> Vec<Shape> {
        self.my_efmap
            .get(&shape_key(the_e))
            .map(|(_, l)| l.clone())
            .unwrap_or_default()
    }
}

// The OCCT private member functions (Draft_Modification.hxx L183-203) —
// InternalAdd / Propagate / Perform / NewSurface(S, ...) / NewCurve(C, ...) —
// live in the sibling submodules offset/draft_modification_1.rs and
// offset/draft_modification_1_b.rs (the single-file <2000-line rule).
