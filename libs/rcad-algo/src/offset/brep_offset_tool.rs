// OCCT BRepOffset_Tool.cxx L1-1437 + BRepOffset_Tool.hxx L17-211 — module a
// of the 1:1 translation (the Tool translation is split by OCCT order into
// brep_offset_tool.rs / brep_offset_tool_b.rs / brep_offset_tool_c.rs).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_Tool.cxx / .hxx
//
// Module a carries cxx L111-L1437: TheInfini, the class methods
// EdgeVertices / Gabarit / OrientSection / FindCommonShapes x2 / PipeInter,
// and the file statics up to AssembleEdge.  Module b (brep_offset_tool_b.rs)
// carries cxx L1441-L2576: Inter3D / TryProject / InterOrExtent / ExtentEdge
// / ProjectVertexOnEdge / Inter2d / SelectEdge.  Module c
// (brep_offset_tool_c.rs) carries cxx L2580-L4659: MakeFace /
// EnlargeGeometry / UpdatePCurves / CompactUVBounds / CheckBounds /
// EnLargeFace / TryParameter / MapVertexEdges / BuildNeighbour / ExtentFace
// / Deboucle3D / IsInOut / CorrectOrientation / CheckPlanesNormals /
// PerformPlanes / UpdateVertexTolerances.
//
// OCCT inheritance chain (hxx L43): none — BRepOffset_Tool is a static-only
// class.
//
// Architecture differences (numbering continues from brep_offset_offset.rs):
// 21. NCollection_List<TopoDS_Shape> -> Vec<Shape>; NCollection_Sequence ->
//     Vec<Shape>; NCollection_Map -> OcctShapeSet (HashMap keyed by
//     (TShape ptr, Location)); NCollection_IndexedMap ->
//     OcctIndexedShapeMap (insertion-ordered, OCCT 1-based operator()
//     preserved); NCollection_DataMap -> ShapeDataMap<V>; the
//     NCollection_IndexedDataMap forms -> ShapeIndexedDataMap<V>.
// 22. BRep_Tool / TopExp / BRep_Builder -> the crate::brep_algo::tool
//     re-hosts (bare-Shape Arc::make_mut edits; no global TShape arena).
// 23. The Handle forms carry no refcount in rcad — the tools are static
//     functions, so only the argument forms change (const& -> &).
// 24. Geom2dInt_GInter (2d curve/curve intersection), ProjLib_ProjectedCurve
//     and GeomProjLib::Curve2d (TKTopAlgo),
//     GeomConvert_CompCurveToBSplineCurve / Geom2dConvert_CompCurveToBSpline
//     Curve, GeomConvert_ApproxCurve / Geom2dConvert_ApproxCurve
//     (TKGeomBase), GeomAPI::To2d/To3d, GCPnts_AbscissaPoint /
//     GCPnts_QuasiUniformDeflection (TKGeomAlgo),
//     ShapeCustom_Curve2d::ConvertToLine2d (TKShHealing) have no rcad
//     translation yet — GAP carriers below close on those batches; the
//     reduced re-hosts keep the OCCT control flow and take the OCCT
//     failure/null paths (the loc_ope_wires_on_shape_b.rs #9-#11 precedent).
// 25. Geom_TrimmedCurve / Geom2d_TrimmedCurve -> TrimmedCurve3 /
//     TrimmedCurve2 wrappers; the IsInstance(STANDARD_TYPE(...)) type
//     probes map onto the Surface3 / Curve2d / Curve3 variant matches (the
//     chfi3d_builder_0.rs surface_type_of precedent).
// 26. Bnd_Box2d -> rcad_kernel::math::bnd::BndBox2d (the
//     loc_ope_wires_on_shape_b.rs #7 precedent).
// 27. GeomInt_IntSS (TKGeomAlgo/GeomInt — the PipeInter / InterOrExtent
//     engine) is not translated; both call sites take the OCCT !IsDone()
//     path (empty result lists) behind the GAP annotation.

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    translate_curve2d, Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, SurfaceEval,
    TrimmedCurve2, TrimmedCurve3,
};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::el::{
    elclib_line_parameter_2d, elslib_cone_parameters, elslib_cone_value,
    elslib_cylinder_parameters, elslib_cylinder_value, elslib_plane_parameters,
    elslib_plane_value, elslib_sphere_parameters, elslib_sphere_value, elslib_torus_parameters,
    elslib_torus_value,
};
use rcad_kernel::topo::topods::{Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_offset::GeomAbsShapeKind;
use super::brep_offset_tool_b::{
    Geom2dConvertApproxCurve, Geom2dConvertCompCurveToBSplineCurve, GeomConvertApproxCurve,
    GeomConvertCompCurveToBSplineCurve,
};

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_range, brep_tool_tolerance, ShapeKey,
};

// ---------------------------------------------------------------------------
// Shared shape-map forms (used by the Inter3d module as well — the OCCT
// dependency direction Inter3d -> Tool is kept).
// ---------------------------------------------------------------------------

/// OCCT NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>
/// (architecture difference #21): insertion-ordered identity set with the
/// OCCT 1-based operator() indexing.
pub(crate) struct OcctIndexedShapeMap {
    index: HashMap<ShapeKey, usize>, // 0-based
    items: Vec<Shape>,
}

impl OcctIndexedShapeMap {
    pub fn new() -> Self {
        OcctIndexedShapeMap { index: HashMap::new(), items: Vec::new() }
    }
    /// OCCT IndexedMap::Add — returns the index (0-based) of the added or
    /// existing key.
    pub fn add(&mut self, s: &Shape) -> usize {
        let k = bat::shape_key(s);
        if let Some(&i) = self.index.get(&k) {
            return i;
        }
        let i = self.items.len();
        self.index.insert(k, i);
        self.items.push(s.clone());
        i
    }
    /// OCCT IndexedMap::Contains.
    pub fn contains(&self, s: &Shape) -> bool {
        self.index.contains_key(&bat::shape_key(s))
    }
    /// OCCT IndexedMap::Extent.
    pub fn extent(&self) -> usize {
        self.items.len()
    }
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
    /// OCCT IndexedMap::operator(i) — 1-based.
    pub fn at_1(&self, i: usize) -> &Shape {
        &self.items[i - 1]
    }
    /// OCCT IndexedMap::FindKey(i) — 1-based.
    pub fn find_key_1(&self, i: usize) -> &Shape {
        &self.items[i - 1]
    }
    pub fn iter(&self) -> impl Iterator<Item = &Shape> {
        self.items.iter()
    }
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.index.clear();
        self.items.clear();
    }
}

/// OCCT NCollection_DataMap<K = TopoDS_Shape, V> value carrier — the key
/// Shape is stored beside the value (the cxx key reads keep working).
pub(crate) type ShapeDataMap<V> = HashMap<ShapeKey, (Shape, V)>;

/// OCCT NCollection_DataMap operations on ShapeDataMap.
pub(crate) mod shape_data_map {
    use super::{Shape, ShapeDataMap, ShapeKey};

    /// OCCT DataMap::Bind.
    pub fn bind<V>(m: &mut ShapeDataMap<V>, k: &Shape, v: V) {
        m.insert(super::bat::shape_key(k), (k.clone(), v));
    }
    /// OCCT DataMap::IsBound.
    pub fn is_bound<V>(m: &ShapeDataMap<V>, k: &Shape) -> bool {
        m.contains_key(&super::bat::shape_key(k))
    }
    /// OCCT DataMap::operator() (const) — asserts when unbound.
    pub fn value<'a, V>(m: &'a ShapeDataMap<V>, k: &Shape) -> &'a V {
        &m.get(&super::bat::shape_key(k)).expect("DataMap::operator() unbound").1
    }
    /// OCCT DataMap::Find — a clone of the value.
    pub fn find<V: Clone>(m: &ShapeDataMap<V>, k: &Shape) -> V {
        m.get(&super::bat::shape_key(k)).expect("DataMap::Find unbound").1.clone()
    }
    /// OCCT DataMap::Seek — an optional value reference.
    pub fn seek<'a, V>(m: &'a ShapeDataMap<V>, k: &Shape) -> Option<&'a V> {
        m.get(&super::bat::shape_key(k)).map(|e| &e.1)
    }
    /// OCCT DataMap::operator() (mutable) / ChangeFind.
    pub fn change_find<'a, V>(m: &'a mut ShapeDataMap<V>, k: &Shape) -> &'a mut V {
        &mut m.get_mut(&super::bat::shape_key(k)).expect("DataMap::ChangeFind unbound").1
    }
    /// OCCT DataMap::UnBind.
    pub fn un_bind<V>(m: &mut ShapeDataMap<V>, k: &Shape) -> bool {
        m.remove(&super::bat::shape_key(k)).is_some()
    }
    /// OCCT DataMap::Bound — bind a default value when absent and return the
    /// mutable value.
    pub fn bound<'a, V: Default>(m: &'a mut ShapeDataMap<V>, k: &Shape) -> &'a mut V {
        let key: ShapeKey = super::bat::shape_key(k);
        &mut m.entry(key).or_insert((k.clone(), V::default())).1
    }
}

/// OCCT NCollection_IndexedDataMap<K = TopoDS_Shape, V> — insertion-ordered
/// (indexmap) with the key Shape carried.
pub(crate) type ShapeIndexedDataMap<V> = indexmap::IndexMap<ShapeKey, (Shape, V)>;

/// OCCT NCollection_IndexedDataMap operations on ShapeIndexedDataMap.
pub(crate) mod shape_indexed_data_map {
    use super::{Shape, ShapeIndexedDataMap};

    /// OCCT IndexedDataMap::Add — binds when absent; returns the index
    /// (0-based) of the added or existing key.
    pub fn add<V>(m: &mut ShapeIndexedDataMap<V>, k: &Shape, v: V) -> usize {
        let key = super::bat::shape_key(k);
        if let Some((i, _, _)) = m.get_full(&key) {
            return i;
        }
        m.insert(key, (k.clone(), v));
        m.len() - 1
    }
    /// OCCT IndexedDataMap::Extent.
    pub fn extent<V>(m: &ShapeIndexedDataMap<V>) -> usize {
        m.len()
    }
    /// OCCT IndexedDataMap::operator(i) (const) — 1-based value.
    pub fn value_1<V>(m: &ShapeIndexedDataMap<V>, i: usize) -> &V {
        &m.get_index(i - 1).expect("IndexedDataMap::operator() out of range").1 .1
    }
    /// OCCT IndexedDataMap::FindKey(i) — 1-based key shape.
    pub fn find_key_1<V>(m: &ShapeIndexedDataMap<V>, i: usize) -> &Shape {
        &m.get_index(i - 1).expect("IndexedDataMap::FindKey out of range").1 .0
    }
    /// OCCT IndexedDataMap::operator(i) (mutable) — 1-based value.
    pub fn value_1_mut<V>(m: &mut ShapeIndexedDataMap<V>, i: usize) -> &mut V {
        &mut m.get_index_mut(i - 1).expect("IndexedDataMap::operator() out of range").1 .1
    }
    /// OCCT IndexedDataMap::Find(k) (const) — asserts when unbound.
    pub fn find<'a, V>(m: &'a ShapeIndexedDataMap<V>, k: &Shape) -> &'a V {
        &m.get(&super::bat::shape_key(k)).expect("IndexedDataMap::Find unbound").1
    }
    /// OCCT IndexedDataMap::ChangeSeek(k).
    pub fn change_seek<'a, V>(m: &'a mut ShapeIndexedDataMap<V>, k: &Shape) -> Option<&'a mut V> {
        m.get_mut(&super::bat::shape_key(k)).map(|e| &mut e.1)
    }
    /// OCCT IndexedDataMap::operator() (mutable) / ChangeFind.
    pub fn change_find<'a, V>(m: &'a mut ShapeIndexedDataMap<V>, k: &Shape) -> &'a mut V {
        &mut m.get_mut(&super::bat::shape_key(k)).expect("IndexedDataMap::ChangeFind unbound").1
    }
}

/// OCCT NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> — the plain
/// (unordered) identity set.
pub(crate) type OcctShapeSet = HashMap<ShapeKey, Shape>;

/// OCCT NCollection_Map::Add — returns true when newly added.
pub(crate) fn set_add(m: &mut OcctShapeSet, s: &Shape) -> bool {
    let k = bat::shape_key(s);
    if m.contains_key(&k) {
        return false;
    }
    m.insert(k, s.clone());
    true
}

/// OCCT NCollection_Map::Contains.
pub(crate) fn set_contains(m: &OcctShapeSet, s: &Shape) -> bool {
    m.contains_key(&bat::shape_key(s))
}

// ---------------------------------------------------------------------------
// File constants.
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool.cxx L111-121 — TheInfini (the maximal surface
/// enlargement value; see the cxx comment block for the 1e+7 rationale).
pub(crate) const THE_INFINI: f64 = 1e+7;

// ---------------------------------------------------------------------------
// TopExp / BRep_Tool / TopoDS re-hosts local to the Tool cluster.
// ---------------------------------------------------------------------------

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) — CumOri = False (the raw stored
/// extremities; null when absent).
pub(crate) fn top_exp_vertices(edg: &Shape) -> (Shape, Shape) {
    let (v1, v2) = bat::top_exp_vertices_raw(edg);
    (v1.unwrap_or_else(Shape::null), v2.unwrap_or_else(Shape::null))
}

/// OCCT TopExp::CommonVertex(E1, E2, CV) (TopExp.cxx L321-338) — the shared
/// extremity of the two edges (null when the edges share no vertex).
pub(crate) fn top_exp_common_vertex(e1: &Shape, e2: &Shape) -> Shape {
    let (e1_v1, e1_v2) = top_exp_vertices(e1);
    let (e2_v1, e2_v2) = top_exp_vertices(e2);
    if e1_v1.is_same(&e2_v1) || e1_v1.is_same(&e2_v2) {
        return e1_v1;
    }
    if e1_v2.is_same(&e2_v1) || e1_v2.is_same(&e2_v2) {
        return e1_v2;
    }
    Shape::null()
}

/// OCCT TopExp::Vertices(W, Vfirst, Vlast) — the null-Shape form of the
/// wire-vertex walk.
pub(crate) fn wire_vertices_null(w: &Shape) -> (Shape, Shape) {
    let (a, b) = bat::top_exp_vertices_wire(w);
    (a.unwrap_or_else(Shape::null), b.unwrap_or_else(Shape::null))
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
pub(crate) fn oriented(s: &Shape, the_or: Orientation) -> Shape {
    bat::oriented(s, the_or)
}

/// OCCT BRep_Tool::Curve(E, f, l) — the (curve, fpar, lpar) triple.
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    bat::brep_tool_curve(edg)
}

/// OCCT BRep_Tool::Surface(F, L) — the face surface.
pub(crate) fn face_surface_of(f: &Shape) -> Option<Surface3> {
    bat::brep_tool_surface(f)
}

/// OCCT TopoDS_Shape::EmptyCopied() — the TShape copy without sub-shapes.
pub(crate) fn empty_copied(r: &Shape) -> Shape {
    bat::empty_copied(r)
}

/// OCCT new Geom_TrimmedCurve(C, f, l) — the trimmed 3d curve form.
pub(crate) fn trimmed_curve3(c: &Curve3, f: f64, l: f64) -> Curve3 {
    Curve3::Trimmed(TrimmedCurve3 {
        curve: Box::new(c.clone()),
        first: f,
        last: l,
    })
}

/// OCCT new Geom2d_TrimmedCurve(C, f, l) — the trimmed pcurve form.
pub(crate) fn trimmed_curve2(c: &Curve2d, f: f64, l: f64) -> Curve2d {
    Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(c.clone()),
        t_min: f,
        t_max: l,
    })
}

/// OCCT Geom_Surface::IsUPeriodic()/IsVPeriodic() on the optional surface.
fn surface_is_u_or_v_periodic(s: &Option<Surface3>) -> bool {
    match s {
        Some(s) => s.is_u_periodic() || s.is_v_periodic(),
        None => false,
    }
}

/// OCCT Geom_Surface::UPeriod() — the pure-math re-host (the
/// loc_ope_wires_on_shape_b.rs surface_u_period form).
pub(crate) fn surface_u_period(s: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match s {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => TAU,
        _ => 0.0,
    }
}

/// OCCT Geom_Surface::VPeriod() — the pure-math re-host (2*PI sphere,
/// 2*minorRadius torus).
pub(crate) fn surface_v_period(s: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match s {
        Surface3::Sphere(_) => TAU,
        Surface3::Torus(t) => 2.0 * t.minor_radius,
        _ => 0.0,
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::EdgeVertices (hxx L50-52; cxx L141-151).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::EdgeVertices(E, V1, V2) (cxx L141-151) — V1 is the
/// first vertex and V2 the last vertex of E taking account of the edge
/// orientation.
pub fn edge_vertices(e: &Shape, v1: &mut Shape, v2: &mut Shape) {
    // OCCT TopExp::Vertices(E, Vfirst, Vlast) (CumOri = False).
    if e.orientation == Orientation::Reversed {
        // OCCT L145: TopExp::Vertices(E, V2, V1).
        let (a_first, a_last) = top_exp_vertices(e);
        *v2 = a_first;
        *v1 = a_last;
    } else {
        // OCCT L149: TopExp::Vertices(E, V1, V2).
        let (a_first, a_last) = top_exp_vertices(e);
        *v1 = a_first;
        *v2 = a_last;
    }
}

// ---------------------------------------------------------------------------
// OCCT static FindPeriod (cxx L155-190).
// ---------------------------------------------------------------------------

/// OCCT static FindPeriod(F, umin, umax, vmin, vmax) (cxx L155-190) — fill a
/// 2d box from a sampled pcurve of every edge of F.
fn find_period(f: &Shape, umin: &mut f64, umax: &mut f64, vmin: &mut f64, vmax: &mut f64) {
    let mut b = BndBox2d::new();
    for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L165: const Handle(Geom2d_Curve) C = BRep_Tool::CurveOnSurface(E, F, pf, pl).
        let (c, mut pf, pl) = match brep_tool_curve_on_surface(&e, f) {
            Some(t) => t,
            // OCCT L167-169: if (C.IsNull()) return.
            None => return,
        };
        // OCCT L170-175: Geom2dAdaptor_Curve PC(C, pf, pl); nbp = 20 (2 for a line).
        let mut nbp: f64 = 20.0;
        if matches!(c, Curve2d::Line(_)) {
            nbp = 2.0;
        }
        // OCCT L176-179: step = (pl - pf) / nbp; PC.D0(pf, P); B.Add(P).
        let step = (pl - pf) / nbp;
        let p = Curve2dEval::point_at(&c, pf);
        b.update(p.x, p.y, p.x, p.y);
        // OCCT L180-185: for (i = 2; i < nbp; i++) { pf += step; PC.D0(pf, P); B.Add(P); }
        let mut i: f64 = 2.0;
        while i < nbp {
            pf += step;
            let p = Curve2dEval::point_at(&c, pf);
            b.update(p.x, p.y, p.x, p.y);
            i += 1.0;
        }
        // OCCT L186-187: PC.D0(pl, P); B.Add(P).
        let p = Curve2dEval::point_at(&c, pl);
        b.update(p.x, p.y, p.x, p.y);
    }
    // OCCT L188: B.Get(umin, vmin, umax, vmax).
    if let Some((g0, g1, g2, g3)) = b.get() {
        *umin = g0;
        *vmin = g1;
        *umax = g2;
        *vmax = g3;
    }
}

// ---------------------------------------------------------------------------
// OCCT static PutInBounds (cxx L197-309).
// ---------------------------------------------------------------------------

/// OCCT static PutInBounds(F, E, C2d) (cxx L197-309) — recadre la courbe 2d
/// dans les bounds de la face.
fn put_in_bounds(f: &Shape, e: &Shape, c2d: &mut Curve2d) {
    // OCCT L201: BRep_Tool::Range(E, f, l).
    let (f_par, l_par) = brep_tool_range(e);

    // OCCT L203-204: TopLoc_Location L; S = BRep_Tool::Surface(F, L) — the
    // rcad face surface is stored location-baked (annotated).
    let mut s = match face_surface_of(f) {
        Some(v) => v,
        None => return,
    };
    // OCCT L206-209: RectangularTrimmedSurface -> BasisSurface.
    if let Surface3::Trimmed(ts) = s {
        s = (*ts.basis).clone();
    }
    //---------------
    // Recadre en U.
    //---------------
    // OCCT L213-216: if (!S->IsUPeriodic() && !S->IsVPeriodic()) return.
    if !s.is_u_periodic() && !s.is_v_periodic() {
        return;
    }

    // OCCT L218: FindPeriod(F, umin, umax, vmin, vmax).
    let (mut umin, mut umax, mut vmin, mut vmax) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    find_period(f, &mut umin, &mut umax, &mut vmin, &mut vmax);

    // OCCT L220-262.
    if s.is_u_periodic() {
        let period = surface_u_period(&s);
        let eps = period * 1e-6;
        let pf = Curve2dEval::point_at(c2d, f_par);
        let pl = Curve2dEval::point_at(c2d, l_par);
        let pm = Curve2dEval::point_at(c2d, 0.34 * f_par + 0.66 * l_par);
        let mut min_c = pf.x.min(pl.x);
        min_c = min_c.min(pm.x);
        let mut max_c = pf.x.max(pl.x);
        max_c = max_c.max(pm.x);
        let mut du = 0.0f64;
        if min_c < umin - eps {
            du = (((umin - min_c) / period) as i32 as f64 + 1.0) * period;
        }
        if min_c > umax + eps {
            du = -(((min_c - umax) / period) as i32 as f64 + 1.0) * period;
        }
        if du != 0.0 {
            // OCCT L242-244: gp_Vec2d T1(du, 0.); C2d->Translate(T1).
            *c2d = translate_curve2d(c2d, DVec2::new(du, 0.0));
            min_c += du;
            max_c += du;
        }
        // Ajuste au mieux la courbe dans le domaine.
        if max_c > umax + 100.0 * eps {
            let d1 = max_c - umax;
            let d2 = umin - min_c + period;
            if d2 < d1 {
                du = -period;
            }
            if du != 0.0 {
                // OCCT L258-260: gp_Vec2d T2(du, 0.); C2d->Translate(T2).
                *c2d = translate_curve2d(c2d, DVec2::new(du, 0.0));
            }
        }
    }
    //------------------
    // Recadre en V.
    //------------------
    // OCCT L266-308.
    if s.is_v_periodic() {
        let period = surface_v_period(&s);
        let eps = period * 1e-6;
        let pf = Curve2dEval::point_at(c2d, f_par);
        let pl = Curve2dEval::point_at(c2d, l_par);
        let pm = Curve2dEval::point_at(c2d, 0.34 * f_par + 0.66 * l_par);
        let mut min_c = pf.y.min(pl.y);
        min_c = min_c.min(pm.y);
        let mut max_c = pf.y.max(pl.y);
        max_c = max_c.max(pm.y);
        let mut dv = 0.0f64;
        if min_c < vmin - eps {
            dv = (((vmin - min_c) / period) as i32 as f64 + 1.0) * period;
        }
        if min_c > vmax + eps {
            dv = -(((min_c - vmax) / period) as i32 as f64 + 1.0) * period;
        }
        if dv != 0.0 {
            // OCCT L288-290: gp_Vec2d T1(0., dv); C2d->Translate(T1).
            *c2d = translate_curve2d(c2d, DVec2::new(0.0, dv));
            min_c += dv;
            max_c += dv;
        }
        // Ajuste au mieux la courbe dans le domaine.
        if max_c > vmax + 100.0 * eps {
            let d1 = max_c - vmax;
            let d2 = vmin - min_c + period;
            if d2 < d1 {
                dv = -period;
            }
            if dv != 0.0 {
                // OCCT L303-305: gp_Vec2d T2(0., dv); C2d->Translate(T2).
                *c2d = translate_curve2d(c2d, DVec2::new(0.0, dv));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::Gabarit (hxx L202; cxx L313-323).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::Gabarit(aCurve) (cxx L313-323) — the bounding-box
/// diameter of the curve.  GAP leaf: BndLib_Add3dCurve::Add (TKTopAlgo) is
/// not translated; the reduced re-host samples the curve (the
/// CompactUVBounds N=33 sampling precedent).
pub fn gabarit(a_curve: &Curve3) -> f64 {
    // OCCT L315-317: GeomAdaptor_Curve GC(aCurve); Bnd_Box aBox;
    // BndLib_Add3dCurve::Add(GC, Precision::Confusion(), aBox).
    // Reduced re-host: uniform 33-point sampling over the curve parameter
    // range carried by the caller's (curve, f, l) view — the OCCT adaptor
    // covers the full basis-curve domain, so the caller passes the range it
    // wants measured.
    let (f, l) = gabarit_range(a_curve);
    let n: usize = 33;
    let mut box_min = DVec3::splat(f64::INFINITY);
    let mut box_max = DVec3::splat(f64::NEG_INFINITY);
    for i in 0..n {
        let t = f + (l - f) * (i as f64) / ((n - 1) as f64);
        let p = CurveEval::point_at(a_curve, t);
        box_min = box_min.min(p);
        box_max = box_max.max(p);
    }
    let d = box_max - box_min;
    // OCCT L318-322: dist = max(aXmax - aXmin, aYmax - aYmin);
    // dist = max(dist, aZmax - aZmin).
    let dist = d.x.max(d.y).max(d.z);
    dist
}

/// The Gabarit sampling range of a curve — the full domain for a bare curve
/// (the GeomAdaptor_Curve first/last parameter form).
fn gabarit_range(c: &Curve3) -> (f64, f64) {
    match c {
        Curve3::Trimmed(tc) => (tc.first, tc.last),
        _ => (0.0, 1.0),
    }
}

// ---------------------------------------------------------------------------
// OCCT static BuildPCurves (cxx L327-482).
// ---------------------------------------------------------------------------

/// OCCT static BuildPCurves(E, F) (cxx L327-482) — compute the pcurve of E
/// on F when absent.  GAP leaves: ProjLib_ProjectedCurve and Extrema_ExtPC
/// are not translated (architecture difference #24); the reduced re-host
/// keeps the OCCT control flow, takes the OCCT "no bound edge found"
/// fall-through of the BSpline/Bezier branch, and the projection GAP keeps
/// the OCCT Standard_ConstructionError path.
pub(crate) fn build_pcurves(e: &Shape, f: &Shape) {
    // OCCT L330: C2d = BRep_Tool::CurveOnSurface(E, F, ff, ll).
    let c2d = brep_tool_curve_on_surface(e, f);
    // OCCT L331-334: if (!C2d.IsNull()) return.
    if c2d.is_some() {
        return;
    }

    // OCCT L337: constexpr double Tolerance = Precision::Confusion();
    let tolerance = rcad_kernel::precision::CONFUSION;

    // OCCT L339-340: BRepAdaptor_Surface AS(F, false); BRepAdaptor_Curve AC(E)
    // — the rcad direct surface/curve reads stand in for the adaptors.

    // OCCT L343-348: the surface type probe (OffsetSurface -> BasisSurface).
    let the_surf = face_surface_of(f);
    let mut typ_s_is_bspl_or_bezier = false;
    if let Some(s) = &the_surf {
        let mut basis = s;
        if let Surface3::Offset(os) = s {
            basis = os.basis.as_ref();
        }
        // OCCT L349: typS == STANDARD_TYPE(Geom_BezierSurface) ||
        //            typS == STANDARD_TYPE(Geom_BSplineSurface).
        typ_s_is_bspl_or_bezier = matches!(basis, Surface3::BSpline(_) | Surface3::Bezier(_));
    }

    if typ_s_is_bspl_or_bezier {
        // OCCT L350-428: the pcurve-on-bound search (Extrema_ExtPC of the
        // edge endpoints onto the face-bounds edges, then a
        // Geom2d_TrimmedCurve of the found bound pcurve + PutInBounds +
        // B.UpdateEdge + BRepLib::SameRange).  GAP leaf: Extrema_ExtPC is
        // not translated (architecture difference #24); the whole branch
        // takes the OCCT "no edge found" fall-through (theEdge stays null).
    } // OCCT L428: } // if (typS == ...

    // OCCT L430-433: ProjLib_ProjectedCurve Proj(HS, HC, Tolerance) — GAP
    // (architecture difference #24).
    //
    // OCCT L435-467: switch (Proj.GetType()) — the Line / Circle / Ellipse /
    // Parabola / Hyperbola / Bezier / BSpline pcurve constructors; with the
    // GAP the projection yields no curve, so C2d stays null and the OCCT
    // L473-481 null-C2d path is taken.
    let c2d_new: Option<Curve2d> = None;
    let _ = tolerance;

    if let Some(mut cc) = c2d_new {
        // OCCT L469-472: if (AS.IsUPeriodic() || AS.IsVPeriodic())
        // PutInBounds(F, E, C2d).
        if surface_is_u_or_v_periodic(&face_surface_of(f)) {
            put_in_bounds(f, e, &mut cc);
        }
        // OCCT L473-477: B.UpdateEdge(E, C2d, F, BRep_Tool::Tolerance(E)).
        let e_tol = brep_tool_tolerance(e);
        bat::builder_update_edge_pcurve(&mut e.clone(), &cc, f, e_tol);
        return;
    }
    // OCCT L478-481: throw Standard_ConstructionError(
    //   "BRepOffset_Tool::BuildPCurves") — the GAP keeps the OCCT failure
    // form (the plan §0.6 annotation).
    panic!("GAP: BRepOffset_Tool::BuildPCurves (ProjLib_ProjectedCurve not translated)");
}

// ---------------------------------------------------------------------------
// OCCT static ToSmall (cxx L622-638).
// ---------------------------------------------------------------------------

/// OCCT static ToSmall(C) (cxx L622-638).
fn to_small(c: &Curve3) -> bool {
    let tol = 10.0 * rcad_kernel::precision::CONFUSION;
    // OCCT L625-628: m = f * 0.668 + l * 0.332; P1 = Value(f); P2 = Value(l);
    // P3 = Value(m) — the (f, l) values are the Geom_Curve parameter range;
    // the rcad trimmed view carries it, bare curves use the [0, 1] view.
    let (f, l) = match c {
        Curve3::Trimmed(tc) => (tc.first, tc.last),
        _ => (0.0, 1.0),
    };
    let m = f * 0.668 + l * 0.332;
    let p1 = CurveEval::point_at(c, f);
    let p2 = CurveEval::point_at(c, l);
    let p3 = CurveEval::point_at(c, m);
    if p1.distance(p2) > tol {
        return false;
    }
    if p2.distance(p3) > tol {
        return false;
    }
    true
}

// ---------------------------------------------------------------------------
// OCCT static IsOnSurface (cxx L642-743).
// ---------------------------------------------------------------------------

/// OCCT static IsOnSurface(C, S, TolConf, TolReached) (cxx L642-743) — the
/// elementary-surface parameter probes through ElSLib; the rcad re-host
/// dispatches on the Surface3 variant and uses the math::el ElSLib re-hosts.
pub(crate) fn is_on_surface(c: &Curve3, s: &Surface3, tol_conf: f64, tol_reached: &mut f64) -> bool {
    // OCCT L647-651: f = C->FirstParameter(); l = C->LastParameter(); n = 5;
    // du = (f - l) / (n - 1); TolReached = 0.
    let f = 0.0f64;
    let l = 1.0f64;
    let n = 5usize;
    let du = (f - l) / (n as f64 - 1.0);
    *tol_reached = 0.0;

    // OCCT L658-740: switch (AS.GetType()).
    match s {
        Surface3::Plane(pl) => {
            // OCCT L660-673: the Plane branch (PlaneParameters / PlaneValue).
            for i in 0..n {
                let p = CurveEval::point_at(c, f + i as f64 * du);
                let (u, v) = elslib_plane_parameters(p, pl.origin, pl.u_dir, pl.v_dir);
                let sp = elslib_plane_value(u, v, pl.origin, pl.u_dir, pl.v_dir);
                *tol_reached = p.distance(sp);
                if *tol_reached > tol_conf {
                    return false;
                }
            }
        }
        Surface3::Cylinder(cy) => {
            // OCCT L674-688: the Cylinder branch.
            for i in 0..n {
                let p = CurveEval::point_at(c, f + i as f64 * du);
                let (u, v) = elslib_cylinder_parameters(
                    p, cy.origin, cy.ref_dir, cy.axis.cross(cy.ref_dir), cy.axis, cy.radius,
                );
                let sp = elslib_cylinder_value(u, v, cy.origin, cy.axis, cy.ref_dir, cy.radius);
                *tol_reached = p.distance(sp);
                if *tol_reached > tol_conf {
                    return false;
                }
            }
        }
        Surface3::Cone(co) => {
            // OCCT L689-704: the Cone branch (RefRadius / SemiAngle).
            for i in 0..n {
                let p = CurveEval::point_at(c, f + i as f64 * du);
                let (u, v) = elslib_cone_parameters(
                    p,
                    co.apex,
                    co.ref_dir,
                    co.axis.cross(co.ref_dir),
                    co.axis,
                    co.radius,
                    co.half_angle_rad,
                );
                let sp = elslib_cone_value(u, v, co.apex, co.axis, co.half_angle_rad, co.radius);
                *tol_reached = p.distance(sp);
                if *tol_reached > tol_conf {
                    return false;
                }
            }
        }
        Surface3::Sphere(sp2) => {
            // OCCT L705-719: the Sphere branch.
            for i in 0..n {
                let p = CurveEval::point_at(c, f + i as f64 * du);
                let (u, v) = elslib_sphere_parameters(
                    p,
                    sp2.center,
                    sp2.ref_dir,
                    sp2.axis.cross(sp2.ref_dir),
                    sp2.axis,
                );
                let sp = elslib_sphere_value(u, v, sp2.center, sp2.axis, sp2.ref_dir, sp2.radius);
                *tol_reached = p.distance(sp);
                if *tol_reached > tol_conf {
                    return false;
                }
            }
        }
        Surface3::Torus(to) => {
            // OCCT L720-735: the Torus branch (MajorRadius / MinorRadius).
            for i in 0..n {
                let p = CurveEval::point_at(c, f + i as f64 * du);
                let (u, v) = elslib_torus_parameters(
                    p,
                    to.center,
                    to.ref_dir,
                    to.axis.cross(to.ref_dir),
                    to.axis,
                    to.major_radius,
                    to.minor_radius,
                );
                let sp = elslib_torus_value(
                    u,
                    v,
                    to.center,
                    to.axis,
                    to.major_radius,
                    to.minor_radius,
                );
                *tol_reached = p.distance(sp);
                if *tol_reached > tol_conf {
                    return false;
                }
            }
        }
        // OCCT L737-739: default: return false.
        _ => {
            return false;
        }
    }

    // OCCT L742: return true.
    true
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::OrientSection (hxx L57-61; cxx L486-568).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::OrientSection(E, F1, F2, O1, O2) (cxx L486-568) —
/// computes O1 the orientation of E in F1 influenced by F2; idem for O2.
pub fn orient_section(
    e: &Shape,
    f1: &Shape,
    f2: &Shape,
    o1: &mut Orientation,
    o2: &mut Orientation,
) {
    // OCCT L492: TopLoc_Location L — the rcad face surface is stored
    // location-baked, so the L.Transformation() applications below are
    // identity (annotated).
    // OCCT L495-499: S1/S2 = the face surfaces; C1/C2 = the pcurves;
    // C = the edge 3d curve with (f, l).
    let s1 = match face_surface_of(f1) {
        Some(v) => v,
        None => return,
    };
    let s2 = match face_surface_of(f2) {
        Some(v) => v,
        None => return,
    };
    let c1 = match brep_tool_curve_on_surface(e, f1) {
        Some((c, _, _)) => c,
        None => return,
    };
    let c2 = match brep_tool_curve_on_surface(e, f2) {
        Some((c, _, _)) => c,
        None => return,
    };
    let c = match brep_tool_curve(e) {
        Some((c, _, _)) => c,
        None => return,
    };

    // OCCT L501-513: ParOnC = GCPnts_AbscissaPoint(BAcurve, Length/2, f) —
    // GAP leaf (architecture difference #24); the reduced re-host takes the
    // OCCT fallback branch BOPTools_AlgoTools2D::IntermediatePoint(f, l)
    // unconditionally (the mid parameter of the trimmed range).
    let (f_par, l_par) = brep_tool_range(e);
    let par_on_c = crate::bop::int_tools::face_make_curve::intermediate_point(f_par, l_par);

    // OCCT L515-519: T1 = C->DN(ParOnC, 1).Transformed(L.Transformation());
    // normalized when above gp::Resolution.
    let mut t1 = CurveEval::derivative_at(&c, par_on_c);
    if t1.length_squared() > f64::MIN_POSITIVE {
        t1 = t1.normalize();
    }

    // OCCT L521-530: DN1 = D1U ^ D1V on F1 at the pcurve point (reversed
    // for a REVERSED face).
    let p = Curve2dEval::point_at(&c1, par_on_c);
    let mut d_n1 = surface_normal(&s1, p.x, p.y);
    if f1.orientation == Orientation::Reversed {
        d_n1 = -d_n1;
    }

    // OCCT L532-538: DN2 (idem on F2).
    let p = Curve2dEval::point_at(&c2, par_on_c);
    let mut d_n2 = surface_normal(&s2, p.x, p.y);
    if f2.orientation == Orientation::Reversed {
        d_n2 = -d_n2;
    }

    // OCCT L540-549: ProVec = DN2 ^ T1; Prod = DN1 . ProVec; O1.
    let pro_vec = d_n2.cross(t1);
    let prod = d_n1.dot(pro_vec);
    if prod < 0.0 {
        *o1 = Orientation::Forward;
    } else {
        *o1 = Orientation::Reversed;
    }
    // OCCT L550-559: ProVec = DN1 ^ T1; Prod = DN2 . ProVec; O2.
    let pro_vec = d_n1.cross(t1);
    let prod = d_n2.dot(pro_vec);
    if prod < 0.0 {
        *o2 = Orientation::Forward;
    } else {
        *o2 = Orientation::Reversed;
    }
    // OCCT L560-567: reverse O1/O2 for REVERSED faces.
    if f1.orientation == Orientation::Reversed {
        *o1 = bat::top_abs_reverse(*o1);
    }
    if f2.orientation == Orientation::Reversed {
        *o2 = bat::top_abs_reverse(*o2);
    }
}

/// OCCT gp_Vec DN1(D1U ^ D1V) — the rcad stand-in evaluates the surface
/// normal through the SurfaceEval derivative vehicle (architecture
/// difference #6 of the MakeSimpleOffset module).
fn surface_normal(s: &Surface3, u: f64, v: f64) -> DVec3 {
    let (_, d1u, d1v) = s.derivatives(u, v);
    d1u.cross(d1v)
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::FindCommonShapes (hxx L63-78; cxx L572-618).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::FindCommonShapes(F1, F2, LE, LV) (cxx L572-580) —
/// looks for the common vertices and edges between the two faces.
pub fn find_common_shapes(
    the_f1: &Shape,
    the_f2: &Shape,
    the_le: &mut Vec<Shape>,
    the_lv: &mut Vec<Shape>,
) -> bool {
    let b_found_edges = find_common_shapes_of_type(the_f1, the_f2, ShapeType::Edge, the_le);
    let b_found_verts = find_common_shapes_of_type(the_f1, the_f2, ShapeType::Vertex, the_lv);
    b_found_edges || b_found_verts
}

/// OCCT BRepOffset_Tool::FindCommonShapes(S1, S2, Type, LSC) (cxx L584-618)
/// — looks for the common shapes of theType between the two shapes.
pub fn find_common_shapes_of_type(
    the_s1: &Shape,
    the_s2: &Shape,
    the_type: ShapeType,
    the_lsc: &mut Vec<Shape>,
) -> bool {
    the_lsc.clear();
    //
    // OCCT L591-596: aMS = the map of the sub-shapes of S1 of theType.
    let mut a_ms: OcctShapeSet = HashMap::new();
    for a_exs in explorer(the_s1, the_type, ShapeType::Shape) {
        set_add(&mut a_ms, &a_exs);
    }
    //
    // OCCT L598-601: if (aMS.IsEmpty()) return false.
    if a_ms.is_empty() {
        return false;
    }
    //
    // OCCT L603-615: the fence-filtered scan of S2.
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    for a_s2 in explorer(the_s2, the_type, ShapeType::Shape) {
        if set_contains(&a_ms, &a_s2) {
            if set_add(&mut a_m_fence, &a_s2) {
                the_lsc.push(a_s2);
            }
        }
    }
    //
    // OCCT L617: return !theLSC.IsEmpty().
    !the_lsc.is_empty()
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::PipeInter (hxx L105-109; cxx L747-804).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::PipeInter(F1, F2, L1, L2, Side) (cxx L747-804) —
/// the GeomInt_IntSS intersection of the two pipe faces.  GAP leaf:
/// GeomInt_IntSS is not translated (architecture difference #27); the
/// carrier keeps the OCCT structure and takes the OCCT !IsDone() path
/// (empty result lists) until the GeomInt batch closes.
pub fn pipe_inter(f1: &Shape, f2: &Shape, l1: &mut Vec<Shape>, l2: &mut Vec<Shape>, side: State) {
    // OCCT L754-757: CI; O1, O2; L1.Clear(); L2.Clear().
    l1.clear();
    l2.clear();

    // OCCT L758: BRep_Builder B.
    // OCCT L759-762: S1/S2 = the face surfaces;
    // GeomInt_IntSS Inter(S1, S2, Precision::Confusion(), true, true, true)
    // — GAP (architecture difference #27).
    let _s1 = face_surface_of(f1);
    let _s2 = face_surface_of(f2);
    let inter_done: bool = false; // GAP: GeomInt_IntSS::Perform
    let _ = inter_done;

    // OCCT L764-803: if (Inter.IsDone()) { for each line: ToSmall filter,
    // BRepLib_MakeEdge(CI), the LineOnS1/LineOnS2 pcurve attachment or
    // BuildPCurves, OrientSection + the Side reversal, the result appends }
    // — behind the OCCT !IsDone() path while the engine is a GAP.
    if inter_done {
        unreachable!("GeomInt_IntSS::IsDone unreachable while the engine is a GAP");
    }
    let _ = side;
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::InterOrExtent (hxx L117-121; cxx L2010-2068).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::InterOrExtent(F1, F2, L1, L2, Side) (cxx
/// L2010-2068) — the GeomInt_IntSS intersection after the trimmed-plane
/// basis collapse.  GAP leaf: GeomInt_IntSS (architecture difference #27);
/// the OCCT !IsDone() path is taken (empty result lists).
pub fn inter_or_extent(f1: &Shape, f2: &Shape, l1: &mut Vec<Shape>, l2: &mut Vec<Shape>, side: State) {
    // OCCT L2017-2020: CI; O1, O2; L1.Clear(); L2.Clear().
    l1.clear();
    l2.clear();
    // OCCT L2021-2041: S1/S2 with the trimmed-plane basis collapse.
    let mut s1 = face_surface_of(f1);
    let mut s2 = face_surface_of(f2);
    if let Some(Surface3::Trimmed(ts)) = &s1 {
        if matches!(ts.basis.as_ref(), Surface3::Plane(_)) {
            s1 = Some((*ts.basis).clone());
        }
    }
    if let Some(Surface3::Trimmed(ts)) = &s2 {
        if matches!(ts.basis.as_ref(), Surface3::Plane(_)) {
            s2 = Some((*ts.basis).clone());
        }
    }
    let _ = (&s1, &s2);

    // OCCT L2043: GeomInt_IntSS Inter(S1, S2, Precision::Confusion()) — GAP.
    let inter_done: bool = false; // GAP: GeomInt_IntSS::Perform

    // OCCT L2045-2067: if (Inter.IsDone()) { for each line: ToSmall filter,
    // BRepLib_MakeEdge(CI), BuildPCurves(E, F1/F2), OrientSection + the
    // Side reversal, the result appends } — behind the OCCT !IsDone() path.
    if inter_done {
        unreachable!("GeomInt_IntSS::IsDone unreachable while the engine is a GAP");
    }
    let _ = side;
}

// ---------------------------------------------------------------------------
// OCCT static SelectEdge (cxx L2536-2576) — placed in module a because
// module c (ExtentFace) consumes it; OCCT order L2536.
// ---------------------------------------------------------------------------

/// OCCT static SelectEdge(F, EF, E, LInt) (cxx L2536-2576) — detrompeur sur
/// les intersections sur les faces periodiques: keep the edge of LInt that
/// covers the initial-edge domain best.
pub(crate) fn select_edge(_f: &Shape, _ef: &Shape, e: &Shape, l_int: &mut Vec<Shape>) {
    // OCCT L2544-2549: dU = 1.0e100; GE; the (Fst, Lst) range of E and its
    // endpoint points PFirst/PLast.
    let mut d_u = 1.0e100f64;
    let mut ge: Option<Shape> = None;

    let (fst0, lst0) = brep_tool_range(e);
    let ad1 = brep_tool_curve(e);
    let p_first = ad1.as_ref().map(|(c, f, _)| CurveEval::point_at(c, *f));
    let p_last = ad1.as_ref().map(|(c, _, l)| CurveEval::point_at(c, *l));
    let _ = (fst0, lst0);

    //----------------------------------------------------------------------
    // Selection de l edge qui couvre le plus le domaine de l edge initiale.
    //----------------------------------------------------------------------
    for it_value in l_int.iter() {
        let ei = it_value;
        // OCCT L2562-2565: BRep_Tool::Range(EI, Fst, Lst); P1/P2.
        let (fst, lst) = brep_tool_range(ei);
        let (p1, p2) = match brep_tool_curve(ei) {
            Some((c, _, _)) => (CurveEval::point_at(&c, fst), CurveEval::point_at(&c, lst)),
            None => continue,
        };
        let (p_first, p_last) = match (p_first, p_last) {
            (Some(a), Some(b)) => (a, b),
            _ => return,
        };
        // OCCT L2567-2572: tmp = P1.Distance(PFirst) + P2.Distance(PLast).
        let tmp = p1.distance(p_first) + p2.distance(p_last);
        if tmp <= d_u {
            d_u = tmp;
            ge = Some(ei.clone());
        }
    }
    // OCCT L2574-2575: LInt.Clear(); LInt.Append(GE).
    l_int.clear();
    if let Some(g) = ge {
        l_int.push(g);
    }
}

// ---------------------------------------------------------------------------
// OCCT ElCLib::Parameter forms on the 2d conics (pure-math re-hosts; the
// math::el module carries only the 2d line parameter).
// ---------------------------------------------------------------------------

/// OCCT ElCLib::Parameter(gp_Circ2d, P) (ElCLib.cxx L1049-1058) —
/// U = atan2((P - Loc) . YDir, (P - Loc) . XDir) with the circle sense.
pub(crate) fn elclib_circle_parameter_2d(
    p: DVec2,
    loc: DVec2,
    x_dir: DVec2,
    y_dir: DVec2,
) -> f64 {
    let d = p - loc;
    y_dir.dot(d).atan2(x_dir.dot(d))
}

/// OCCT ElCLib::Parameter(gp_Elips2d, P) (ElCLib.cxx L1072-1086) — the
/// ellipse parametrization u from the major/minor axes.
pub(crate) fn elclib_ellipse_parameter_2d(
    p: DVec2,
    loc: DVec2,
    x_dir: DVec2,
    y_dir: DVec2,
    major: f64,
    minor: f64,
) -> f64 {
    let d = p - loc;
    let x = x_dir.dot(d);
    let y = y_dir.dot(d);
    let u0 = (y / minor).atan2(x / major);
    // OCCT: u = u0 with sense; the sign branch keeps the ellipse sense.
    if u0 < 0.0 {
        u0 + 2.0 * std::f64::consts::PI
    } else {
        u0
    }
}

/// OCCT ElCLib::Parameter(gp_Parab2d, P) (ElCLib.cxx L1110-1124) —
/// u = y² focal form: u = (P - Loc) . YDir squared over the focal length
/// form u = focal * tan²(u/2) inversion -> u = 2 * atan(sqrt(d.y *
/// focal…)); the direct ElCLib form: U = focal * ((P.Loc).Y / focal)²
/// reduced to u = d.y * d.y / (4 * focal).
pub(crate) fn elclib_parabola_parameter_2d(p: DVec2, loc: DVec2, y_dir: DVec2, focal: f64) -> f64 {
    let d = p - loc;
    let y = y_dir.dot(d);
    y * y / (4.0 * focal)
}

/// OCCT ElCLib::Parameter(gp_Hypr2d, P) (ElCLib.cxx L1160-1180) —
/// u = acosh form from the major/minor radii.
pub(crate) fn elclib_hyperbola_parameter_2d(
    p: DVec2,
    loc: DVec2,
    x_dir: DVec2,
    y_dir: DVec2,
    major: f64,
    minor: f64,
) -> f64 {
    let d = p - loc;
    let x = x_dir.dot(d) / major;
    let y = y_dir.dot(d) / minor;
    let m = (x * x - y * y).max(1.0);
    let u = m.acosh();
    if x < 0.0 {
        -u
    } else {
        u
    }
}

// ---------------------------------------------------------------------------
// OCCT static IsAutonomVertex (cxx L811-832) — the (vertex, PDS, face1,
// face2) form over the FaceInfo VerticesOn maps.
// ---------------------------------------------------------------------------

/// OCCT static IsAutonomVertex(theVertex, thePDS, theFace1, theFace2) (cxx
/// L811-832) — checks whether the vertex is "autonom" (not aVerticesOn of
/// either face).
///
/// DS-context annotation (the batch-2a form): OCCT holds the arguments by
/// TopoDS identity inside the PDS; the rcad DS deep-clones the arguments, so
/// `DS::index` resolves the caller's handle through the (ptr, location)
/// map — the caller must pass the argument handles as produced by the DS
/// (the same discipline as the bop Builder consumers).
pub(crate) fn is_autonom_vertex_on_faces(
    the_vertex: &Shape,
    the_pds: &crate::bop::ds::DS,
    the_face1: &Shape,
    the_face2: &Shape,
) -> bool {
    // OCCT L816-819: nV = thePDS->Index(theVertex); nF[2] = the face indices.
    let n_v = match the_pds.index(the_vertex) {
        i if i >= 0 => i as usize,
        _ => usize::MAX,
    };
    let n_f = [
        the_pds.index(the_face1),
        the_pds.index(the_face2),
    ];

    // OCCT L821-829: for each face, if the FaceInfo VerticesOn map contains
    // nV the vertex is not autonom.
    for i in 0..2 {
        if n_f[i] < 0 {
            continue;
        }
        let a_face_info = the_pds.face_info(n_f[i] as usize);
        if a_face_info.vertices_on.contains(&n_v) {
            return false;
        }
    }

    // OCCT L831: return true.
    true
}

// ---------------------------------------------------------------------------
// OCCT static IsAutonomVertex (cxx L839-911) — the (vertex, PDS) form over
// the VV/EE/EF interference arrays.
// ---------------------------------------------------------------------------

/// OCCT static IsAutonomVertex(aVertex, pDS) (cxx L839-911) — checks whether
/// the vertex with the given index is not created in a VV, EE or EF
/// interference.
pub(crate) fn is_autonom_vertex(a_vertex: &Shape, p_ds: &crate::bop::ds::DS) -> bool {
    // OCCT L844-859: index = pDS->Index(aVertex); the fall-back scan over
    // the new shapes.
    let mut index = p_ds.index(a_vertex);
    if index == -1 {
        let i1 = p_ds.nb_source_shapes();
        let i2 = p_ds.nb_shapes();
        for i in i1..i2 {
            let a_sx = p_ds.shape(i);
            if a_sx.is_same(a_vertex) {
                index = i as isize;
                break;
            }
        }
    }
    if index < 0 {
        // OCCT would pass -1 into IsNewShape (false) — keep the guard.
        return false;
    }
    let index = index as usize;
    //
    // OCCT L861-864: if (!pDS->IsNewShape(index)) return false.
    if !p_ds.is_new_shape(index) {
        return false;
    }
    // check if vertex with index "index" is not created in VV or EE or EF
    // interference.
    // VV
    // OCCT L867-879: aVVs = pDS->InterfVV(); HasIndexNew / IndexNew == index.
    // rcad: the InterferenceVV carries the new vertex as `merged_vertex`
    // (the fused-vertex index; usize::MAX when absent).
    for a_vv in &p_ds.interf_vv {
        if a_vv.merged_vertex == index {
            return false;
        }
    }
    // EE
    // OCCT L881-894: aEEs = pDS->InterfEE(); CommonPart Type == VERTEX and
    // IndexNew == index.  rcad: the InterferenceEE `new_vertex` field is
    // only populated for the vertex-type common parts (usize::MAX
    // otherwise).
    for a_ee in &p_ds.interf_ee {
        if a_ee.new_vertex == index {
            return false;
        }
    }
    // EF
    // OCCT L896-909: aEFs = pDS->InterfEF(); idem.
    for a_ef in &p_ds.interf_ef {
        if a_ef.new_vertex == index {
            return false;
        }
    }
    // OCCT L910: return true.
    true
}

// ---------------------------------------------------------------------------
// OCCT static AreConnex (cxx L918-925).
// ---------------------------------------------------------------------------

/// OCCT static AreConnex(W1, W2) (cxx L918-925) — whether two wires are
/// connex by a vertex.
pub(crate) fn are_connex(w1: &Shape, w2: &Shape) -> bool {
    let (v11, v12) = wire_vertices_null(w1);
    let (v21, v22) = wire_vertices_null(w2);

    v11.is_same(&v21) || v11.is_same(&v22) || v12.is_same(&v21) || v12.is_same(&v22)
}

// ---------------------------------------------------------------------------
// OCCT static AreClosed (cxx L932-939).
// ---------------------------------------------------------------------------

/// OCCT static AreClosed(E1, E2) (cxx L932-939) — whether two edges are
/// connex by two vertices.
fn are_closed(e1: &Shape, e2: &Shape) -> bool {
    let (v11, v12) = top_exp_vertices(e1);
    let (v21, v22) = top_exp_vertices(e2);

    (v11.is_same(&v21) && v12.is_same(&v22)) || (v11.is_same(&v22) && v12.is_same(&v21))
}

// ---------------------------------------------------------------------------
// OCCT static BSplineEdges (cxx L943-987).
// ---------------------------------------------------------------------------

/// OCCT static BSplineEdges(E1, E2, par1, par2, angle) (cxx L943-987) —
/// computes the angle between the two edges' BSpline basis curves at their
/// first/last parameters.
fn bspline_edges(e1: &Shape, e2: &Shape, par1: i32, par2: i32, angle: &mut f64) -> bool {
    // OCCT L951-961: C1/C2 = BRep_Tool::Curve with the trimmed basis
    // extraction.
    let (c1, first1, last1) = match brep_tool_curve(e1) {
        Some(t) => t,
        None => return false,
    };
    let c1 = basis_curve3(&c1);
    let (c2, first2, last2) = match brep_tool_curve(e2) {
        Some(t) => t,
        None => return false,
    };
    let c2 = basis_curve3(&c2);

    // OCCT L963-967: both must be BSpline.
    let (b1, b2) = match (&c1, &c2) {
        (Curve3::BSpline(b1), Curve3::BSpline(b2)) => (b1, b2),
        _ => return false,
    };
    let _ = (b1, b2);

    // OCCT L969-970: Param1/Param2 = first or last by the par flag.
    let param1 = if par1 == 0 { first1 } else { last1 };
    let param2 = if par2 == 0 { first2 } else { last2 };

    // OCCT L972-975: D1 at the parameters.
    let der1 = CurveEval::derivative_at(&c1, param1);
    let der2 = CurveEval::derivative_at(&c2, param2);

    // OCCT L977-984: the degenerate-derivative / angle form.
    if der1.length() <= f64::MIN_POSITIVE || der2.length() <= f64::MIN_POSITIVE {
        *angle = std::f64::consts::PI / 2.0;
    } else {
        *angle = der1.angle_between(der2);
    }

    // OCCT L986: return true.
    true
}

/// OCCT C->IsInstance(STANDARD_TYPE(Geom_TrimmedCurve)) ->
/// BasisCurve() — the rcad variant match.
pub(crate) fn basis_curve3(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Trimmed(tc) => (*tc.curve).clone(),
        _ => c.clone(),
    }
}

/// OCCT PCurve->IsInstance(STANDARD_TYPE(Geom2d_TrimmedCurve)) ->
/// BasisCurve() — the rcad variant match.
pub(crate) fn basis_curve2(c: &Curve2d) -> Curve2d {
    match c {
        Curve2d::Trimmed(tc) => (*tc.curve).clone(),
        _ => c.clone(),
    }
}

// ---------------------------------------------------------------------------
// OCCT static AngleWireEdge (cxx L991-1031).
// ---------------------------------------------------------------------------

/// OCCT static AngleWireEdge(aWire, anEdge) (cxx L991-1031) — the angle
/// between the wire edge adjacent to the shared vertex and anEdge.
pub(crate) fn angle_wire_edge(a_wire: &Shape, an_edge: &Shape) -> f64 {
    // OCCT L993-996: the wire/edge vertices; CV = the shared vertex.
    let (mut v11, mut v12) = bat::top_exp_vertices_wire(a_wire);
    let (v21, v22) = top_exp_vertices(an_edge);
    let w_v11 = v11.clone().unwrap_or_else(Shape::null);
    let w_v12 = v12.clone().unwrap_or_else(Shape::null);
    let mut cv = if w_v11.is_same(&v21) || w_v11.is_same(&v22) {
        w_v11.clone()
    } else {
        w_v12.clone()
    };
    // OCCT L997-1010: FirstEdge = the wire edge touching CV.
    let mut first_edge = Shape::null();
    for e in bat::sub_shapes(a_wire) {
        let cand = e;
        let (v1o, v2o) = top_exp_vertices(&cand);
        if v1o.is_same(&cv) || v2o.is_same(&cv) {
            v11 = Some(v1o.clone());
            v12 = Some(v2o.clone());
            first_edge = cand;
            break;
        }
    }
    let v11 = v11.unwrap_or_else(Shape::null);
    let v12 = v12.unwrap_or_else(Shape::null);
    // OCCT L1011-1029: the four (Vx, Vy) sharing branches.
    let mut angle = 0.0f64;
    if v11.is_same(&cv) && v21.is_same(&cv) {
        bspline_edges(&first_edge, an_edge, 0, 0, &mut angle);
        angle = std::f64::consts::PI - angle;
    } else if v11.is_same(&cv) && v22.is_same(&cv) {
        bspline_edges(&first_edge, an_edge, 0, 1, &mut angle);
    } else if v12.is_same(&cv) && v21.is_same(&cv) {
        bspline_edges(&first_edge, an_edge, 1, 0, &mut angle);
    } else {
        bspline_edges(&first_edge, an_edge, 1, 1, &mut angle);
        angle = std::f64::consts::PI - angle;
    }
    cv = cv; // the OCCT CV local stays live through the function (no-op mirror).
    angle
}

// ---------------------------------------------------------------------------
// OCCT static ReconstructPCurves (cxx L1035-1058).
// ---------------------------------------------------------------------------

/// OCCT static ReconstructPCurves(anEdge) (cxx L1035-1058) — for every
/// curve-on-surface representation of the edge, re-project the 3d curve
/// onto the representation surface.  GAP leaves (architecture difference
/// #24): GeomProjLib::Curve2d and the per-representation surface storage
/// (the rcad edge carries pcurves keyed by the face only); the reduced
/// re-host keeps the OCCT control flow as a no-op (the tolerance rework of
/// the Inter3D L1910-1917 branch keeps the OCCT SameParameter re-run
/// behind it).
pub(crate) fn reconstruct_pcurves(an_edge: &Shape) {
    // OCCT L1037-1038: C3d = BRep_Tool::Curve(anEdge, f, l).
    let _c3d = brep_tool_curve(an_edge);
    // OCCT L1040-1057: the representation loop (IsCurveOnSurface ->
    // Surface/Location -> Transformed -> GeomProjLib::Curve2d ->
    // CurveRep->PCurve) — GAP no-op (annotated above).
}

// ---------------------------------------------------------------------------
// OCCT static ConcatPCurves (cxx L1062-1161).
// ---------------------------------------------------------------------------

/// OCCT static ConcatPCurves(E1, E2, F, After, newFirst, newLast) (cxx
/// L1062-1161) — the concatenated pcurve of the two edges on F.
fn concat_pcurves(
    e1: &Shape,
    e2: &Shape,
    f: &Shape,
    after: bool,
    new_first: &mut f64,
    new_last: &mut f64,
) -> Option<Curve2d> {
    // OCCT L1069-1072: Tol = 1e-7; Continuity = GeomAbs_C1; MaxDeg = 14;
    // MaxSeg = 16.
    let tol = 1e-7;
    let continuity = GeomAbsShapeKind::C1;
    let max_deg = 14;
    let max_seg = 16;

    // OCCT L1077-1087: PCurve1/PCurve2 = CurveOnSurface with the trimmed
    // basis extraction.
    let (pcurve1_raw, first1, last1) = match brep_tool_curve_on_surface(e1, f) {
        Some(t) => t,
        None => return None,
    };
    let pcurve1 = basis_curve2(&pcurve1_raw);
    let (pcurve2_raw, mut first2, mut last2) = match brep_tool_curve_on_surface(e2, f) {
        Some(t) => t,
        None => return None,
    };
    let pcurve2 = basis_curve2(&pcurve2_raw);

    // OCCT L1089-1094: same pcurve handle -> the range union.
    // Reduced re-host: the OCCT handle identity (PCurve1 == PCurve2) maps to
    // the rcad structural identity of the variant payloads (the same TShape
    // pcurve storage produces equal clones).
    if pcurve_handle_same_2d(&pcurve1, &pcurve2) {
        *new_first = first1.min(first2);
        *new_last = last1.max(last2);
        return Some(pcurve1);
    }
    // OCCT L1095-1140: same dynamic type and (Line or Conic) -> recompute
    // the E2 parameters on the PCurve1 carrier.
    let same_kind = std::mem::discriminant(&pcurve1) == std::mem::discriminant(&pcurve2);
    let is_line = matches!(pcurve1, Curve2d::Line(_));
    let is_conic = matches!(
        pcurve1,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) | Curve2d::Parabola(_) | Curve2d::Hyperbola(_)
    );
    if same_kind && (is_line || is_conic) {
        // newPCurve = PCurve1.
        let p1 = Curve2dEval::point_at(&pcurve2, first2);
        let p2 = Curve2dEval::point_at(&pcurve2, last2);
        if let Curve2d::Line(lin1) = &pcurve1 {
            // OCCT L1103-1109: the Geom2d_Line branch.
            first2 = elclib_line_parameter_2d(p1, lin1.origin, lin1.direction);
            last2 = elclib_line_parameter_2d(p2, lin1.origin, lin1.direction);
        } else if let Curve2d::Circle(circ1) = &pcurve1 {
            // OCCT L1110-1116: the Geom2d_Circle branch.
            first2 = elclib_circle_parameter_2d(p1, circ1.center, circ1.x_dir, circ1.y_dir);
            last2 = elclib_circle_parameter_2d(p2, circ1.center, circ1.x_dir, circ1.y_dir);
        } else if let Curve2d::Ellipse(el1) = &pcurve1 {
            // OCCT L1117-1123: the Geom2d_Ellipse branch.
            first2 = elclib_ellipse_parameter_2d(
                p1,
                el1.center,
                el1.major_dir,
                el1.minor_dir,
                el1.major_radius,
                el1.minor_radius,
            );
            last2 = elclib_ellipse_parameter_2d(
                p2,
                el1.center,
                el1.major_dir,
                el1.minor_dir,
                el1.major_radius,
                el1.minor_radius,
            );
        } else if let Curve2d::Parabola(par1c) = &pcurve1 {
            // OCCT L1124-1130: the Geom2d_Parabola branch.
            let par_y_dir = par1c.axis_dir.perp();
            first2 = elclib_parabola_parameter_2d(p1, par1c.origin, par_y_dir, par1c.focal_param);
            last2 = elclib_parabola_parameter_2d(p2, par1c.origin, par_y_dir, par1c.focal_param);
        } else if let Curve2d::Hyperbola(hy1) = &pcurve1 {
            // OCCT L1131-1137: the Geom2d_Hyperbola branch.
            let hy_y_dir = hy1.major_dir.perp();
            first2 = elclib_hyperbola_parameter_2d(
                p1,
                hy1.center,
                hy1.major_dir,
                hy_y_dir,
                hy1.semi_major,
                hy1.semi_minor,
            );
            last2 = elclib_hyperbola_parameter_2d(
                p2,
                hy1.center,
                hy1.major_dir,
                hy_y_dir,
                hy1.semi_major,
                hy1.semi_minor,
            );
        }
        *new_first = first1.min(first2);
        *new_last = last1.max(last2);
        return Some(pcurve1);
    }
    // OCCT L1141-1158: the generic concatenation branch — GAP leaves
    // (architecture difference #24): Geom2dConvert_CompCurveToBSplineCurve
    // / Geom2dConvert_ApproxCurve are not translated; the OCCT carrier
    // panics carry the annotation.
    let _ = (tol, continuity, max_deg, max_seg, after);
    let tc1 = trimmed_curve2(&pcurve1, first1, last1);
    let tc2 = trimmed_curve2(&pcurve2, first2, last2);
    let mut concat2d = Geom2dConvertCompCurveToBSplineCurve::new(&tc1);
    if !concat2d.add(&tc2, rcad_kernel::precision::CONFUSION, after) {
        // OCCT keeps the (possibly partial) BSpline; the carrier panics
        // before this point — unreachable.
        return None;
    }
    let mut new_pcurve = match concat2d.bspline_curve() {
        Some(c) => c,
        None => return None,
    };
    if new_pcurve_continuity(&new_pcurve) < 1 {
        let approx2d = Geom2dConvertApproxCurve::new(&new_pcurve, tol, continuity, max_seg, max_deg);
        if approx2d.has_result() {
            new_pcurve = approx2d.curve();
        }
    }
    *new_first = 0.0;
    *new_last = 1.0;
    Some(new_pcurve)
}

/// The OCCT handle identity (PCurve1 == PCurve2) — the rcad structural
/// stand-in (see the ConcatPCurves annotation).
fn pcurve_handle_same_2d(a: &Curve2d, b: &Curve2d) -> bool {
    // The rcad geometry payloads carry no PartialEq — the Debug canonical
    // form is the structural stand-in.
    std::mem::discriminant(a) == std::mem::discriminant(b)
        && format!("{:?}", a) == format!("{:?}", b)
}

/// OCCT Geom2d_Curve::Continuity() < GeomAbs_C1 probe — the rcad variant
/// stand-in (analytic conics are C-infinite, BSpline knots decide; the
/// reduced form reports C1+ for analytics and C0 for BSpline).
fn new_pcurve_continuity(c: &Curve2d) -> i32 {
    match c {
        Curve2d::BSpline(_) => 0,
        _ => 2,
    }
}

// ---------------------------------------------------------------------------
// OCCT static Glue (cxx L1165-1260).
// ---------------------------------------------------------------------------

/// OCCT static Glue(E1, E2, Vfirst, Vlast, After, F1, addPCurve1, F2,
/// addPCurve2, theGlueTol) (cxx L1165-1260) — the glued edge of E1+E2.
#[allow(clippy::too_many_arguments)]
fn glue(
    e1: &Shape,
    e2: &Shape,
    vfirst: &Shape,
    vlast: &Shape,
    after: bool,
    f1: &Shape,
    add_pcurve1: bool,
    f2: &Shape,
    add_pcurve2: bool,
    the_glue_tol: f64,
) -> Shape {
    let new_edge = bat::builder_make_edge();

    // OCCT L1178-1181: Tol = 1e-7; Continuity = GeomAbs_C1; MaxDeg = 14;
    // MaxSeg = 16.
    let tol = 1e-7;
    let continuity = GeomAbsShapeKind::C1;
    let max_deg = 14;
    let max_seg = 16;

    // OCCT L1188-1198: C1/C2 = BRep_Tool::Curve with the trimmed basis
    // extraction.
    let (c1_raw, first1, last1) = match brep_tool_curve(e1) {
        Some(t) => t,
        None => return new_edge,
    };
    let c1 = basis_curve3(&c1_raw);
    let (c2_raw, first2, last2) = match brep_tool_curve(e2) {
        Some(t) => t,
        None => return new_edge,
    };
    let c2 = basis_curve3(&c2_raw);

    let mut is_canonic = false;
    let (new_curve, fparam, lparam) =
        // OCCT L1200-1232: the same-handle / same-canonic-type / generic
        // concatenation branches.
        if curve3_handle_same(&c1, &c2) {
            // OCCT L1200-1205: newCurve = C1; range union.
            (c1.clone(), first1.min(first2), last1.max(last2))
        } else if std::mem::discriminant(&c1) == std::mem::discriminant(&c2)
            && (matches!(c1, Curve3::Line(_)) || c1_is_conic(&c1))
        {
            // OCCT L1206-1211: IsCanonic = true; newCurve = C1.
            is_canonic = true;
            (c1.clone(), 0.0, 0.0)
        } else {
            // OCCT L1212-1232: the generic concatenation — GAP leaves
            // (architecture difference #24): GeomConvert_CompCurveToBSpline
            // Curve / GeomConvert_ApproxCurve.
            let tc1 = trimmed_curve3(&c1, first1, last1);
            let tc2 = trimmed_curve3(&c2, first2, last2);
            let mut concat = GeomConvertCompCurveToBSplineCurve::new(&tc1);
            if !concat.add(&tc2, the_glue_tol, after) {
                // OCCT L1218-1220: return newEdge (null).
                return new_edge;
            }
            let new_curve = match concat.bspline_curve() {
                Some(c) => c,
                None => return new_edge,
            };
            let mut new_curve = new_curve;
            if new_curve3_continuity(&new_curve) < 1 {
                let approx3d =
                    GeomConvertApproxCurve::new(&new_curve, tol, continuity, max_seg, max_deg);
                if approx3d.has_result() {
                    new_curve = approx3d.curve();
                }
            }
            let fparam = 0.0;
            let lparam = 1.0;
            (new_curve, fparam, lparam)
        };

    // OCCT L1234-1243: BRepLib_MakeEdge(newCurve, Vfirst, Vlast[, fparam,
    // lparam]) — the rcad BRepBuilder::add_edge form (the IsDone guard has
    // no rcad counterpart, architecture difference #5).
    //
    // The OCCT constructor is replaced by the bare-Shape assembly over the
    // new edge (the crate::brep_algo::tool layer): the vertices and the
    // range are attached directly.
    let mut new_edge = new_edge;
    if let rcad_kernel::topo::topods::TShape::Edge(ed) =
        std::sync::Arc::make_mut(&mut new_edge.data)
    {
        ed.curve = Some(new_curve.clone());
        if is_canonic {
            ed.range = [fparam.min(lparam), fparam.max(lparam)];
        } else {
            ed.range = [fparam, lparam];
        }
        ed.first = oriented(vfirst, Orientation::Forward);
        ed.last = oriented(vlast, Orientation::Reversed);
        ed.my_shapes.push(oriented(vfirst, Orientation::Forward));
        ed.my_shapes.push(oriented(vlast, Orientation::Reversed));
    }

    // OCCT L1245-1257: the pcurve re-attachment.
    let mut new_first = 0.0f64;
    let mut new_last = 0.0f64;
    if add_pcurve1 {
        if let Some(new_pcurve) = concat_pcurves(e1, e2, f1, after, &mut new_first, &mut new_last) {
            bat::builder_update_edge_pcurve(&mut new_edge, &new_pcurve, f1, 0.0);
            bat::builder_range_edge_on_face(&mut new_edge, f1, new_first, new_last);
        }
    }
    if add_pcurve2 {
        if let Some(new_pcurve) = concat_pcurves(e1, e2, f2, after, &mut new_first, &mut new_last) {
            bat::builder_update_edge_pcurve(&mut new_edge, &new_pcurve, f2, 0.0);
            bat::builder_range_edge_on_face(&mut new_edge, f2, new_first, new_last);
        }
    }

    // OCCT L1259: return newEdge.
    new_edge
}

/// OCCT C1 == C2 (the Geom_Curve handle identity) — the rcad structural
/// stand-in (see the Glue annotation).
fn curve3_handle_same(a: &Curve3, b: &Curve3) -> bool {
    // The rcad geometry payloads carry no PartialEq — the Debug canonical
    // form is the structural stand-in.
    std::mem::discriminant(a) == std::mem::discriminant(b)
        && format!("{:?}", a) == format!("{:?}", b)
}

/// OCCT C1->IsKind(STANDARD_TYPE(Geom_Conic)) — the rcad variant match.
fn c1_is_conic(c: &Curve3) -> bool {
    matches!(c, Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Parabola(_) | Curve3::Hyperbola(_))
}

/// OCCT Geom_Curve::Continuity() < GeomAbs_C1 probe — the rcad stand-in
/// (see new_pcurve_continuity).
fn new_curve3_continuity(c: &Curve3) -> i32 {
    match c {
        Curve3::BSpline(_) => 0,
        _ => 2,
    }
}

// ---------------------------------------------------------------------------
// OCCT static CheckIntersFF (cxx L1264-1341).
// ---------------------------------------------------------------------------

/// OCCT static CheckIntersFF(pDS, RefEdge, TrueEdges) (cxx L1264-1341) —
/// the FF section-edge selection around RefEdge (the reduced
/// MakeConnexityBlocks / BRepExtrema GAP forms carry the OCCT annotation).
pub(crate) fn check_inters_ff(
    p_ds: &crate::bop::ds::DS,
    ref_edge: &Shape,
    true_edges: &mut OcctIndexedShapeMap,
) {
    // OCCT L1268-1270: aFFs = pDS->InterfFF(); nbe = 0.
    let a_ffs = &p_ds.interf_ff;
    let mut nbe = 0usize;

    // OCCT L1272-1274: Edges = a compound; BB.MakeCompound(Edges).
    let mut edges = bat::builder_make_compound();

    // OCCT L1276-1299: collect the FF section edges (per InterfFF curve per
    // PaveBlock -> the DS edge).
    for a_ff in a_ffs {
        for &a_bc in &a_ff.curves {
            let a_curve = &p_ds.intersection_curves[a_bc];
            for pb in &a_curve.pave_blocks {
                let n_sect = pb.read().edge();
                let an_edge = p_ds.shape(n_sect);
                bat::builder_add_compound_shape(&mut edges, &an_edge);
                nbe += 1;
            }
        }
    }

    // OCCT L1301-1304: if (nbe == 0) return.
    if nbe == 0 {
        return;
    }

    // OCCT L1306-1307: CompList; BOPTools_AlgoTools::MakeConnexityBlocks(
    // Edges, TopAbs_VERTEX, TopAbs_EDGE, CompList).
    //
    // Reduced re-host: the (VERTEX, EDGE) instantiation of
    // MakeConnexityBlocks — BFS grouping of the edges sharing vertices (the
    // shell_splitter.rs form is the face-instantiation; annotated).
    let edge_list = bat::sub_shapes(&edges);
    let comp_list = make_connexity_blocks_edges(&edge_list);

    // OCCT L1309-1338: the nearest-compound selection.
    let nearest_compound: Shape = if comp_list.len() == 1 {
        // OCCT L1310-1313.
        comp_list[0].clone()
    } else {
        // OCCT L1315-1337: the BRepExtrema_DistShapeShape distance of the
        // RefEdge mid-vertex to each compound — GAP leaf (architecture
        // difference #24); the reduced re-host takes the OCCT
        // !IsDone() continue for every compound (NearestCompound stays
        // null) and the MapShapes of a null shape yields no entries.
        Shape::null()
    };

    // OCCT L1340: TopExp::MapShapes(NearestCompound, TopAbs_EDGE, TrueEdges).
    let _ = ref_edge;
    for e in explorer(&nearest_compound, ShapeType::Edge, ShapeType::Shape) {
        true_edges.add(&e);
    }
}

/// OCCT BOPTools_AlgoTools::MakeConnexityBlocks(S, TopAbs_VERTEX,
/// TopAbs_EDGE, LC) — the (VERTEX, EDGE) instantiation reduced re-host:
/// BFS grouping of the edges sharing extremity vertices
/// (BOPTools_AlgoTools.cxx L187-256 form).
pub(crate) fn make_connexity_blocks_edges(edges: &[Shape]) -> Vec<Shape> {
    let mut used: OcctShapeSet = HashMap::new();
    let mut comp_list: Vec<Shape> = Vec::new();
    // vertex key -> incident edge indices
    let mut vmap: HashMap<ShapeKey, Vec<usize>> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        let (v1, v2) = top_exp_vertices(e);
        if !v1.is_null() {
            vmap.entry(bat::shape_key(&v1)).or_default().push(i);
        }
        if !v2.is_null() {
            vmap.entry(bat::shape_key(&v2)).or_default().push(i);
        }
    }
    for start in 0..edges.len() {
        if set_contains(&used, &edges[start]) {
            continue;
        }
        let mut block = bat::builder_make_compound();
        let mut stack = vec![start];
        set_add(&mut used, &edges[start]);
        while let Some(i) = stack.pop() {
            bat::builder_add_compound_shape(&mut block, &edges[i]);
            let (v1, v2) = top_exp_vertices(&edges[i]);
            for v in [v1, v2] {
                if v.is_null() {
                    continue;
                }
                if let Some(nbrs) = vmap.get(&bat::shape_key(&v)) {
                    for &j in nbrs {
                        if set_add(&mut used, &edges[j]) {
                            stack.push(j);
                        }
                    }
                }
            }
        }
        comp_list.push(block);
    }
    comp_list
}

// ---------------------------------------------------------------------------
// OCCT static AssembleEdge (cxx L1345-1437).
// ---------------------------------------------------------------------------

/// OCCT static AssembleEdge(pDS, F1, F2, addPCurve1, addPCurve2,
/// EdgesForConcat) (cxx L1345-1437) — the glued assembly of the concat
/// edge sequence.
pub(crate) fn assemble_edge(
    p_ds: &crate::bop::ds::DS,
    f1: &Shape,
    f2: &Shape,
    add_pcurve1: bool,
    add_pcurve2: bool,
    edges_for_concat: &[Shape],
) -> Shape {
    let null_edge = Shape::null();
    if edges_for_concat.is_empty() {
        // OCCT reads EdgesForConcat(1) on a non-empty sequence; keep the
        // guard for the empty-argument form.
        return null_edge;
    }
    // OCCT L1353-1354: CurEdge = EdgesForConcat(1); aGlueTol = Confusion.
    let mut cur_edge = edges_for_concat[0].clone();
    let mut a_glue_tol = rcad_kernel::precision::CONFUSION;

    // OCCT L1356-1434: for (j = 2; j <= EdgesForConcat.Length(); j++).
    for j in 1..edges_for_concat.len() {
        let an_edge = edges_for_concat[j].clone();
        let mut after = false;
        let mut vfirst = Shape::null();
        let mut vlast = Shape::null();
        let are_closed_wire = are_closed(&cur_edge, &an_edge);
        if are_closed_wire {
            // OCCT L1364-1381: the closed-wire branch.
            let (v1, v2) = top_exp_vertices(&cur_edge);
            let is_autonom_v1 = is_autonom_vertex_on_faces(&v1, p_ds, f1, f2);
            let is_autonom_v2 = is_autonom_vertex_on_faces(&v2, p_ds, f1, f2);
            if is_autonom_v1 {
                after = false;
                vfirst = v2.clone();
                vlast = v2.clone();
            } else if is_autonom_v2 {
                after = true;
                vfirst = v1.clone();
                vlast = v1.clone();
            } else {
                // OCCT L1380: return NullEdge.
                return null_edge;
            }
        } else {
            // OCCT L1385-1421: the open-wire branch.
            let cv = top_exp_common_vertex(&cur_edge, &an_edge);
            let mut is_autonom_cv = false;
            if !cv.is_null() {
                is_autonom_cv = is_autonom_vertex_on_faces(&cv, p_ds, f1, f2);
            }
            if is_autonom_cv {
                a_glue_tol = brep_tool_tolerance(&cv);
                let (v11, v12) = top_exp_vertices(&cur_edge);
                let (v21, v22) = top_exp_vertices(&an_edge);
                if v11.is_same(&cv) && v21.is_same(&cv) {
                    vfirst = v22;
                    vlast = v12;
                } else if v11.is_same(&cv) && v22.is_same(&cv) {
                    vfirst = v21;
                    vlast = v12;
                } else if v12.is_same(&cv) && v21.is_same(&cv) {
                    vfirst = v11;
                    vlast = v22;
                } else {
                    vfirst = v11;
                    vlast = v21;
                }
            } else {
                // OCCT L1420: return NullEdge.
                return null_edge;
            }
        } // OCCT L1422: end of else (open wire)

        // OCCT L1424-1433: Glue + the CurEdge advance.
        let new_edge = glue(
            &cur_edge,
            &an_edge,
            &vfirst,
            &vlast,
            after,
            f1,
            add_pcurve1,
            f2,
            add_pcurve2,
            a_glue_tol,
        );
        if new_edge.is_null() {
            return null_edge;
        } else {
            cur_edge = new_edge;
        }
    } // OCCT L1434: end of for

    // OCCT L1436: return CurEdge.
    cur_edge
}
