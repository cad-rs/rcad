// OCCT BRepTools_WireExplorer.hxx L34-110 + BRepTools_WireExplorer.cxx L17-873
// — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingData/TKBRep/BRepTools/BRepTools_WireExplorer.hxx
//         $OCCT_SRC/src/ModelingData/TKBRep/BRepTools/BRepTools_WireExplorer.cxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap myMap (vertex -> list of edges) -> the keyed
//    NCollection_DataMap re-host `MyMap` (HashMap<ShapeKey, (Shape, Vec)>);
//    NCollection_Map myDoubles / emap -> the insertion-ordered ShapeKeySet
//    below (the OCCT bucket iteration order is not reproduced — the same
//    reduction as feat::loc_ope_spliter::ShapeSet).
//    NCollection_IndexedMap vmap -> the same insertion-ordered carrier, whose
//    Add returns the 1-based index of an already-bound key (the
//    NCollection_IndexedMap::Add contract read by the `currsize >= ind` test).
// 2. BRep_Tool::Surface/CurveOnSurface/Tolerance/Range/Degenerated ->
//    crate::brep_algo::tool + crate::feat::loc_ope_wires_on_shape_b re-hosts
//    (feat loc_ope_gluer.rs arch. diff. #1).
// 3. GeomAdaptor_Surface (over the face surface) UResolution / VResolution /
//    GetType / D1 -> u_resolution_for_surface / v_resolution_for_surface
//    (rcad_kernel::topo::topods) + SurfaceEval (loc_ope_gluer.rs arch.
//    diff. #2 / #4).
// 4. TopExp::Vertices(E, V1, V2, CumOri=true) / FirstVertex(E, true) /
//    LastVertex(E, true) -> top_exp_vertices / top_exp_first_vertex /
//    top_exp_last_vertex (loc_ope_wires_on_shape.rs re-hosts).
// 5. BRepTools::UVBounds(F, UMin, UMax, VMin, VMax) -> brep_tools_uv_bounds
//    (loc_ope_wires_on_shape_b.rs re-host).
// 6. The OCCT_DEBUG / OCCT_DEBUG_MESH compile-time branches are not
//    translated (not compiled in the reference build).

use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_range, brep_tool_surface, brep_tool_tolerance, oriented,
};
use crate::feat::brep_feat_builder::sub_shapes;
use crate::feat::loc_ope_wires_on_shape::{
    top_exp_first_vertex, top_exp_last_vertex, top_exp_vertices,
};
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_degenerated;
use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, Surface3, SurfaceEval};
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::precision::{CONFUSION, REAL_LAST};
use rcad_kernel::topo::topods::{u_resolution_for_surface, v_resolution_for_surface};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::HashMap;

pub(crate) type ShapeKey = (u64, u32);

/// OCCT TopTools_ShapeMapHasher identity key (TShape + Location,
/// orientation ignored; architecture difference #1).
pub(crate) fn shape_key(s: &Shape) -> ShapeKey {
    (s.ptr_id(), s.location)
}

/// OCCT NCollection_Map / NCollection_IndexedMap over TopoDS_Shape — an
/// insertion-ordered carrier keyed by the shape identity (architecture
/// difference #1).
#[derive(Default)]
pub(crate) struct ShapeKeySet {
    keys: Vec<ShapeKey>,
    items: HashMap<ShapeKey, Shape>,
}

impl ShapeKeySet {
    pub fn new() -> Self {
        ShapeKeySet {
            keys: Vec::new(),
            items: HashMap::new(),
        }
    }

    /// OCCT NCollection_Map::Add / NCollection_IndexedMap::Add — returns
    /// true when newly added (NCollection_Map) or the 1-based index of the
    /// entry (NCollection_IndexedMap; already-bound keys keep their index).
    pub fn add(&mut self, the_s: &Shape) -> usize {
        let k = shape_key(the_s);
        if let Some(pos) = self.keys.iter().position(|&x| x == k) {
            return pos + 1;
        }
        self.keys.push(k);
        self.items.insert(k, the_s.clone());
        self.keys.len()
    }

    /// OCCT NCollection_Map::Contains.
    pub fn contains(&self, the_s: &Shape) -> bool {
        self.items.contains_key(&shape_key(the_s))
    }

    /// OCCT NCollection_IndexedMap::RemoveKey.
    pub fn remove(&mut self, the_s: &Shape) {
        let k = shape_key(the_s);
        if self.items.remove(&k).is_some() {
            self.keys.retain(|&x| x != k);
        }
    }

    /// OCCT NCollection_Map::Extent / IsEmpty.
    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// OCCT NCollection_Map / IndexedMap::FindKey(i) — 1-based.
    pub fn find_key(&self, the_i: usize) -> Option<&Shape> {
        self.items.get(self.keys.get(the_i - 1)?)
    }

    /// OCCT the first key of the set (the map-iterator head).
    pub fn first(&self) -> Option<&Shape> {
        self.items.get(self.keys.first()?)
    }

    /// OCCT NCollection_Map::Clear.
    pub fn clear(&mut self) {
        self.keys.clear();
        self.items.clear();
    }
}

/// OCCT NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher> — the key Shape travels with the list (the OCCT
/// iterator reads `itm.Key()`); architecture difference #1.
type MyMap = HashMap<ShapeKey, (Shape, Vec<Shape>)>;

/// OCCT BRepTools_WireExplorer (BRepTools_WireExplorer.hxx L47-110).
pub(crate) struct WireExplorer {
    my_map: MyMap,             // OCCT: myMap
    my_edge: Shape,            // OCCT: myEdge
    my_vertex: Shape,          // OCCT: myVertex
    my_face: Shape,            // OCCT: myFace
    my_doubles: ShapeKeySet,   // OCCT: myDoubles
    my_reverse: bool,          // OCCT: myReverse
    my_tol_u: f64,             // OCCT: myTolU
    my_tol_v: f64,             // OCCT: myTolV
}

impl Default for WireExplorer {
    fn default() -> Self {
        Self::new()
    }
}

impl WireExplorer {
    /// OCCT BRepTools_WireExplorer::BRepTools_WireExplorer() (cxx L59-64).
    pub fn new() -> Self {
        WireExplorer {
            my_map: HashMap::new(),
            my_edge: Shape::null(),
            my_vertex: Shape::null(),
            my_face: Shape::null(),
            my_doubles: ShapeKeySet::new(),
            my_reverse: false,
            my_tol_u: 0.0,
            my_tol_v: 0.0,
        }
    }

    /// OCCT BRepTools_WireExplorer(W, F) (cxx L76-79).
    pub fn with_wire_face(the_w: &Shape, the_f: &Shape) -> Self {
        let mut e = WireExplorer::new();
        e.init_wire_face(the_w, the_f);
        e
    }

    /// OCCT BRepTools_WireExplorer::Init(W, F) (cxx L89-109).
    pub fn init_wire_face(&mut self, the_w: &Shape, the_f: &Shape) {
        // OCCT L91-93.
        self.my_edge = Shape::null();
        self.my_vertex = Shape::null();
        self.my_map.clear();
        self.my_doubles.clear();

        // OCCT L95-97: if (W.IsNull()) return;
        if the_w.is_null() {
            return;
        }

        // OCCT L99-108.
        let mut u_min = 0.0;
        let mut u_max = 0.0;
        let mut v_min = 0.0;
        let mut v_max = 0.0;
        if !the_f.is_null() {
            // For the faces based on Cone, BSpline and Bezier compute the
            // UV bounds to precise the UV tolerance values (OCCT L102-107).
            let a_surf_type = geom_abs_surface_type(brep_tool_surface(the_f).as_ref());
            if a_surf_type == GeomAbsSurfaceType::Cone
                || a_surf_type == GeomAbsSurfaceType::BSplineSurface
                || a_surf_type == GeomAbsSurfaceType::BezierSurface
            {
                if let Some(b) = crate::feat::loc_ope_wires_on_shape_b::brep_tools_uv_bounds(the_f)
                {
                    u_min = b[0];
                    u_max = b[1];
                    v_min = b[2];
                    v_max = b[3];
                }
            }
        }

        // OCCT L111-114: Init(W, F, UMin, UMax, VMin, VMax).
        self.init_wire_face_bounds(the_w, the_f, u_min, u_max, v_min, v_max);
    }

    /// OCCT BRepTools_WireExplorer::Init(W, F, UMin, UMax, VMin, VMax)
    /// (cxx L117-373).
    pub fn init_wire_face_bounds(
        &mut self,
        the_w: &Shape,
        the_f: &Shape,
        u_min: f64,
        u_max: f64,
        v_min: f64,
        v_max: f64,
    ) {
        // OCCT L119-127.
        self.my_edge = Shape::null();
        self.my_vertex = Shape::null();
        self.my_map.clear();
        self.my_doubles.clear();

        // OCCT L129-132.
        if the_w.is_null() {
            return;
        }

        // OCCT L134-138.
        self.my_face = the_f.clone();
        let mut df_vert_toler = 0.0;
        self.my_reverse = false;

        if !self.my_face.is_null() {
            // OCCT L141-152.
            let a_surf = brep_tool_surface(&self.my_face).expect("face surface");
            for a_v in crate::feat::brep_feat_builder::explorer(
                the_w,
                ShapeType::Vertex,
                ShapeType::Shape,
            ) {
                df_vert_toler = brep_tool_tolerance(&a_v).max(df_vert_toler);
            }
            if df_vert_toler < CONFUSION {
                // Use tolerance of edges (OCCT L156-166).
                for it in sub_shapes(the_w) {
                    df_vert_toler = brep_tool_tolerance(&it).max(df_vert_toler);
                }
                if df_vert_toler < CONFUSION {
                    // empty wire
                    return;
                }
            }
            // OCCT L170-171.
            self.my_tol_u = 2.0 * u_resolution_for_surface(&a_surf, df_vert_toler);
            self.my_tol_v = 2.0 * v_resolution_for_surface(&a_surf, df_vert_toler);

            // uresolution for cone with infinite vmin vmax is too small
            // (OCCT L173-197).
            if geom_abs_surface_type(Some(&a_surf)) == GeomAbsSurfaceType::Cone {
                let (_, d1u, _) = a_surf.derivatives(u_min, v_min);
                let mut tol1;
                let mut tol2;
                let maxtol = 0.0005 * (u_max - u_min);
                let a = d1u.length();
                if a <= CONFUSION {
                    tol1 = maxtol;
                } else {
                    tol1 = maxtol.min(df_vert_toler / a);
                }

                let (_, d1u, _) = a_surf.derivatives(u_min, v_max);
                let a = d1u.length();
                if a <= CONFUSION {
                    tol2 = maxtol;
                } else {
                    tol2 = maxtol.min(df_vert_toler / a);
                }
                self.my_tol_u = 2.0 * tol1.max(tol2);
            }

            // OCCT L199-214.
            if geom_abs_surface_type(Some(&a_surf)) == GeomAbsSurfaceType::BSplineSurface
                || geom_abs_surface_type(Some(&a_surf)) == GeomAbsSurfaceType::BezierSurface
            {
                let mut max_tol = self.my_tol_u.max(self.my_tol_v);
                let (_, a_du, a_dv) = a_surf.derivatives(0.5 * (u_max - u_min), 0.5 * (v_max - v_min));
                let md = (a_du.dot(a_du) + a_dv.dot(a_dv)).sqrt();
                if md > GP_RESOLUTION {
                    if md * max_tol / df_vert_toler < 1.5 {
                        max_tol = 1.5 * df_vert_toler / md;
                    }
                    self.my_tol_u = max_tol;
                    self.my_tol_v = max_tol;
                }
            }

            // OCCT L216.
            self.my_reverse = self.my_face.orientation == Orientation::Reversed;
        }

        // map of vertices to know if the wire is open (OCCT L219).
        let mut vmap = ShapeKeySet::new();
        // map of infinite edges (OCCT L221).
        let mut an_inf_emap = ShapeKeySet::new();

        // OCCT L225-281: list the vertices.
        for e in sub_shapes(the_w) {
            let e_ori = e.orientation;
            if e_ori == Orientation::Internal || e_ori == Orientation::External {
                continue;
            }
            // OCCT L236: TopExp::Vertices(E, V1, V2, true).
            let (v1_opt, v2_opt) = top_exp_vertices(&e);
            let v1 = v1_opt.unwrap_or_else(Shape::null);
            let v2 = v2_opt.unwrap_or_else(Shape::null);

            if !v1.is_null() {
                // OCCT L240-241.
                self.my_map
                    .entry(shape_key(&v1))
                    .or_insert_with(|| (v1.clone(), Vec::new()))
                    .1
                    .push(e.clone());
                // add or remove in the vertex map (OCCT L244-249).
                let v1f = oriented(&v1, Orientation::Forward);
                let currsize = vmap.extent();
                let ind = vmap.add(&v1f);
                if currsize >= ind {
                    vmap.remove(&v1f);
                }
            }

            if !v2.is_null() {
                // OCCT L252-257.
                let v2r = oriented(&v2, Orientation::Reversed);
                let currsize = vmap.extent();
                let ind = vmap.add(&v2r);
                if currsize >= ind {
                    vmap.remove(&v2r);
                }
            }

            if v1.is_null() || v2.is_null() {
                // OCCT L259-279.
                let (a_f, a_l) = brep_tool_range(&e);
                if e_ori == Orientation::Forward {
                    if a_f == -rcad_kernel::precision::INFINITE_VALUE {
                        an_inf_emap.add(&e);
                    }
                } else if a_l == rcad_kernel::precision::INFINITE_VALUE {
                    an_inf_emap.add(&e);
                }
            }
        }

        // Construction of the set of double edges (OCCT L283-295).
        let mut emap = ShapeKeySet::new();
        for it2 in sub_shapes(the_w) {
            if !emap.contains(&it2) {
                emap.add(&it2);
            } else {
                self.my_doubles.add(&it2);
            }
        }

        // if vmap is not empty the wire is open, let us find the first vertex
        // (OCCT L297-322).
        let mut v1 = Shape::null();
        if !vmap.is_empty() {
            for ind in 1..=vmap.extent() {
                let v = vmap.find_key(ind).expect("vmap entry");
                if v.orientation == Orientation::Forward {
                    v1 = v.clone();
                    break;
                }
            }
        } else {
            // The wire is infinite. Try to find the first vertex. It may be
            // NULL (OCCT L324-356).
            if !an_inf_emap.is_empty() {
                for ind in 1..=an_inf_emap.extent() {
                    let an_edge = an_inf_emap.find_key(ind).expect("anInfEmap entry").clone();
                    let an_ori = an_edge.orientation;
                    let (a_f, a_l) = brep_tool_range(&an_edge);
                    if (an_ori == Orientation::Forward
                        && a_f == -rcad_kernel::precision::INFINITE_VALUE)
                        || (an_ori == Orientation::Reversed
                            && a_l == rcad_kernel::precision::INFINITE_VALUE)
                    {
                        self.my_edge = an_edge;
                        self.my_vertex = Shape::null();
                        return;
                    }
                }
            }

            // use the first vertex in iterator (OCCT L350-356).
            for e in sub_shapes(the_w) {
                let e_ori = e.orientation;
                if e_ori == Orientation::Internal || e_ori == Orientation::External {
                    continue;
                }
                let (a_v1, _) = top_exp_vertices(&e);
                v1 = a_v1.unwrap_or_else(Shape::null);
                break;
            }
        }

        // OCCT L358-366.
        if v1.is_null() {
            return;
        }
        if !self.my_map.contains_key(&shape_key(&v1)) {
            return;
        }

        // OCCT L368-372.
        let l = self.my_map.get_mut(&shape_key(&v1)).expect("myMap(V1)");
        self.my_edge = l.1.first().expect("myMap(V1) non-empty").clone();
        l.1.remove(0);
        self.my_vertex = top_exp_first_vertex(&self.my_edge).expect("FirstVertex");
    }

    /// OCCT BRepTools_WireExplorer::More() (cxx L377-380).
    pub fn more(&self) -> bool {
        !self.my_edge.is_null()
    }

    /// OCCT BRepTools_WireExplorer::Next() (cxx L384-670).
    pub fn next(&mut self) {
        // OCCT L386.
        self.my_vertex = top_exp_last_vertex(&self.my_edge).unwrap_or_else(Shape::null);

        // OCCT L388-393.
        if self.my_vertex.is_null() {
            self.my_edge = Shape::null();
            return;
        }
        if !self.my_map.contains_key(&shape_key(&self.my_vertex)) {
            self.my_edge = Shape::null();
            return;
        }

        // OCCT L395.
        let l_key = shape_key(&self.my_vertex);
        let l_len = self.my_map.get(&l_key).expect("myMap(myVertex)").1.len();

        if l_len == 0 {
            // OCCT L397-400.
            self.my_edge = Shape::null();
        } else if l_len == 1 {
            // Modified by Sergey KHROMOV - Fri Jun 21 10:28:01 2002 OCC325
            // (OCCT L401-451).
            let a_next_edge = self.my_map.get(&l_key).expect("myMap").1[0].clone();
            let (a_v1_opt, a_v2_opt) = top_exp_vertices(&a_next_edge);
            let a_v1 = a_v1_opt.unwrap_or_else(Shape::null);
            let a_v2 = a_v2_opt.unwrap_or_else(Shape::null);

            if !a_v1.is_same(&self.my_vertex) {
                self.my_edge = Shape::null();
                return;
            }
            if !self.my_face.is_null() && a_v1.is_same(&a_v2) {
                let a_prev_pc = brep_tool_curve_on_surface(&self.my_edge, &self.my_face);
                let a_next_pc = brep_tool_curve_on_surface(&a_next_edge, &self.my_face);
                let (Some((a_prev_pc, a_par11, a_par12)), Some((a_next_pc, a_par21, a_par22))) =
                    (a_prev_pc, a_next_pc)
                else {
                    self.my_edge = Shape::null();
                    return;
                };

                let a_prev_par = if self.my_edge.orientation == Orientation::Forward {
                    a_par12
                } else {
                    a_par11
                };

                let (a_next_f_par, a_next_l_par) =
                    if a_next_edge.orientation == Orientation::Forward {
                        (a_par21, a_par22)
                    } else {
                        (a_par22, a_par21)
                    };

                let a_p_prev = a_prev_pc.point_at(a_prev_par);
                let a_p_next_f = a_next_pc.point_at(a_next_f_par);
                let a_p_next_l = a_next_pc.point_at(a_next_l_par);

                if (a_p_prev - a_p_next_f).length_squared()
                    > (a_p_prev - a_p_next_l).length_squared()
                {
                    self.my_edge = Shape::null();
                    return;
                }
            }
            // Modified by Sergey KHROMOV - Fri Jun 21 11:08:16 2002 End
            self.my_edge = self.my_map.get(&l_key).expect("myMap").1[0].clone();
            self.my_map.get_mut(&l_key).expect("myMap").1.clear();
        } else {
            if self.my_face.is_null() {
                // Without Face - try to return edges as logically as possible
                // (OCCT L452-489).
                let mut e = self.my_edge.clone();
                if select_degenerated(&mut self.my_map.get_mut(&l_key).expect("myMap").1, &mut e) {
                    self.my_edge = e;
                    return;
                }
                // At second double edges.
                e = self.my_edge.clone();
                if select_double(
                    &self.my_doubles,
                    &mut self.my_map.get_mut(&l_key).expect("myMap").1,
                    &mut e,
                ) {
                    self.my_edge = e;
                    return;
                }

                let l = self.my_map.get_mut(&l_key).expect("myMap");
                let mut notfound = true;
                let mut pos = None;
                for (i, v) in l.1.iter().enumerate() {
                    if !v.is_same(&self.my_edge) {
                        self.my_edge = v.clone();
                        pos = Some(i);
                        notfound = false;
                        break;
                    }
                }
                if let Some(i) = pos {
                    l.1.remove(i);
                }
                if notfound {
                    self.my_edge = Shape::null();
                }
                return;
            }

            // If we have more than one edge attached to the list probably
            // wire that we explore contains a loop or loops (OCCT L491-669).
            let Some((a_pcurvel, df_f_par, df_l_par)) =
                brep_tool_curve_on_surface(&self.my_edge, &self.my_face)
            else {
                self.my_edge = Shape::null();
                return;
            };
            let a_pcurve = a_pcurvel;

            // Get 2D point equals to <myVertex> in 2D for current edge.
            // OCCT L522-529.
            let p_ref = if self.my_edge.orientation == Orientation::Reversed {
                a_pcurve.point_at(df_f_par)
            } else {
                a_pcurve.point_at(df_l_par)
            };

            // Get next 2D point from current edge's PCurve with parameter
            // F + dP (REV) or L - dP (FOR) (OCCT L531-532).
            let isrevese = self.my_edge.orientation == Orientation::Reversed;
            let df_m_par = get_next_param_on_pc(
                &a_pcurve,
                p_ref,
                df_f_par,
                df_l_par,
                self.my_tol_u,
                self.my_tol_v,
                isrevese,
            );

            let p_refm = a_pcurve.point_at(df_m_par);
            // Get vector from PRef to PRefm (OCCT L537).
            let an_e_ref_dir = p_refm - p_ref;
            if an_e_ref_dir.length_squared() < GP_RESOLUTION {
                self.my_edge = Shape::null();
                return;
            }

            // Search the list of edges looking for the edge having nearest
            // 2D point of connected vertex to current one and smallest angle
            // (OCCT L543-544). First process all degenerated edges, then -
            // all others.
            let mut k = 1usize;
            let mut k_min = 0usize;
            let mut is_degenerated = true;
            let mut dmin = REAL_LAST;
            let mut df_min_angle = 3.0 * std::f64::consts::PI;
            let mut df_cur_angle = 3.0 * std::f64::consts::PI;

            for _i_done in 0..2 {
                let l = self.my_map.get(&l_key).expect("myMap").1.clone();
                k = 1;
                for e in l.iter() {
                    let e = e.clone();
                    if e.is_same(&self.my_edge) {
                        k += 1;
                        continue;
                    }

                    let (a_vert1_opt, a_vert2_opt) = top_exp_vertices(&e);
                    let a_vert1 = a_vert1_opt.unwrap_or_else(Shape::null);
                    let a_vert2 = a_vert2_opt.unwrap_or_else(Shape::null);
                    if a_vert1.is_null() || a_vert2.is_null() {
                        k += 1;
                        continue;
                    }

                    let Some((a_pcurve_e, df_f_par_e, df_l_par_e)) =
                        brep_tool_curve_on_surface(&e, &self.my_face)
                    else {
                        k += 1;
                        continue;
                    };

                    if a_vert1.is_same(&a_vert2) == is_degenerated {
                        // OCCT L589-594.
                        let a_p_eb = if e.orientation == Orientation::Reversed {
                            a_pcurve_e.point_at(df_l_par_e)
                        } else {
                            a_pcurve_e.point_at(df_f_par_e)
                        };

                        if (df_l_par_e - df_f_par_e).abs() > rcad_kernel::precision::PCONFUSION
                        {
                            // OCCT L596-620.
                            let mut isrev = e.orientation == Orientation::Reversed;
                            isrev = !isrev;
                            let a_e_pm = get_next_param_on_pc(
                                &a_pcurve_e,
                                a_p_eb,
                                df_f_par_e,
                                df_l_par_e,
                                self.my_tol_u,
                                self.my_tol_v,
                                isrev,
                            );

                            let mut a_p_ee = a_pcurve_e.point_at(a_e_pm);
                            if (a_p_eb - a_p_ee).length_squared() <= GP_RESOLUTION {
                                // seems to be very short curve
                                let (_, d1) = curve2d_d1(&a_pcurve_e, a_e_pm);
                                a_p_ee = if e.orientation == Orientation::Reversed {
                                    DVec2::new(a_p_eb.x - d1.x, a_p_eb.y - d1.y)
                                } else {
                                    DVec2::new(a_p_eb.x + d1.x, a_p_eb.y + d1.y)
                                };

                                if (a_p_eb - a_p_ee).length_squared() <= GP_RESOLUTION {
                                    k += 1;
                                    continue;
                                }
                            }
                            let an_e_dir = a_p_ee - a_p_eb;
                            df_cur_angle = vec2_angle(an_e_dir, an_e_ref_dir).abs();
                        }

                        if df_cur_angle <= df_min_angle {
                            let mut d = (p_ref - a_p_eb).length_squared();
                            if d <= rcad_kernel::precision::PCONFUSION {
                                d = 0.0;
                            }
                            if (a_p_eb.x - p_ref.x).abs() < self.my_tol_u
                                && (a_p_eb.y - p_ref.y).abs() < self.my_tol_v
                            {
                                if d <= dmin {
                                    df_min_angle = df_cur_angle;
                                    k_min = k;
                                    dmin = d;
                                }
                            }
                        }
                    }
                    k += 1;
                }

                // OCCT L649-657.
                if k_min == 0 {
                    is_degenerated = false;
                    k = 1;
                    dmin = REAL_LAST;
                } else {
                    break;
                }
            }

            if k_min == 0 {
                // probably unclosed in 2d space wire (OCCT L661-665).
                self.my_edge = Shape::null();
                return;
            }

            // Selection the edge (OCCT L667-...).
            let l = self.my_map.get_mut(&l_key).expect("myMap");
            if k_min >= 1 && k_min <= l.1.len() {
                self.my_edge = l.1[k_min - 1].clone();
                l.1.remove(k_min - 1);
            }
        }
    }

    /// OCCT BRepTools_WireExplorer::Current() (cxx L675-678).
    pub fn current(&self) -> &Shape {
        &self.my_edge
    }

    /// OCCT BRepTools_WireExplorer::CurrentVertex() (cxx L705-708).
    pub fn current_vertex(&self) -> &Shape {
        &self.my_vertex
    }

    /// OCCT BRepTools_WireExplorer::Clear() (cxx L712-719).
    pub fn clear(&mut self) {
        self.my_map.clear();
        self.my_doubles.clear();
        self.my_edge = Shape::null();
        self.my_face = Shape::null();
        self.my_vertex = Shape::null();
    }
}

// ---------------------------------------------------------------------------
// Static helpers (cxx L737-873).
// ---------------------------------------------------------------------------

/// OCCT static SelectDouble(Doubles, L, E) (cxx L737-755).
fn select_double(the_doubles: &ShapeKeySet, the_l: &mut Vec<Shape>, the_e: &mut Shape) -> bool {
    let mut pos = None;
    for (i, ce) in the_l.iter().enumerate() {
        if the_doubles.contains(ce) && !the_e.is_same(ce) {
            *the_e = ce.clone();
            pos = Some(i);
            break;
        }
    }
    match pos {
        Some(i) => {
            the_l.remove(i);
            true
        }
        None => false,
    }
}

/// OCCT static SelectDegenerated(L, E) (cxx L759-778).
fn select_degenerated(the_l: &mut Vec<Shape>, the_e: &mut Shape) -> bool {
    let mut pos = None;
    for (i, v) in the_l.iter().enumerate() {
        if !v.is_same(the_e) {
            *the_e = v.clone();
            if brep_tool_degenerated(the_e) {
                pos = Some(i);
                break;
            }
        }
    }
    match pos {
        Some(i) => {
            the_l.remove(i);
            true
        }
        None => false,
    }
}

/// OCCT static GetNextParamOnPC(aPC, aPRef, fP, lP, tolU, tolV, reverse)
/// (cxx L782-873).
fn get_next_param_on_pc(
    the_pc: &rcad_kernel::geom::Curve2d,
    the_p_ref: DVec2,
    the_f: f64,
    the_l: f64,
    the_tol_u: f64,
    the_tol_v: f64,
    the_reverse: bool,
) -> f64 {
    let mut result = if the_reverse { the_f } else { the_l };
    let d_par = (the_l - the_f).abs() / 1000.0;

    if the_reverse {
        let mut start_par = the_f;
        let mut next_pnt_on_edge = false;
        while !next_pnt_on_edge && start_par < the_l {
            start_par += d_par;
            let pnt = the_pc.point_at(start_par);
            if (the_p_ref.x - pnt.x).abs() < the_tol_u && (the_p_ref.y - pnt.y).abs() < the_tol_v {
                continue;
            } else {
                result = start_par;
                next_pnt_on_edge = true;
                break;
            }
        }
        if !next_pnt_on_edge {
            result = the_l;
        }
        if result > the_l {
            result = the_l;
        }
    } else {
        let mut start_par = the_l;
        let mut next_pnt_on_edge = false;
        while !next_pnt_on_edge && start_par > the_f {
            start_par -= d_par;
            let pnt = the_pc.point_at(start_par);
            if (the_p_ref.x - pnt.x).abs() < the_tol_u && (the_p_ref.y - pnt.y).abs() < the_tol_v {
                continue;
            } else {
                result = start_par;
                next_pnt_on_edge = true;
                break;
            }
        }
        if !next_pnt_on_edge {
            result = the_f;
        }
        if result < the_f {
            result = the_f;
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Local re-hosts (arch. diffs. #2/#3).
// ---------------------------------------------------------------------------

/// OCCT GeomAbs_SurfaceType of a face surface (GeomAdaptor_Surface::GetType).
#[derive(PartialEq, Eq, Clone, Copy)]
pub(crate) enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BezierSurface,
    BSplineSurface,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    OffsetSurface,
    OtherSurface,
}

/// OCCT GeomAdaptor_Surface::GetType (GeomAdaptor_Surface.cxx) over the rcad
/// Surface3 enumeration.
pub(crate) fn geom_abs_surface_type(the_surf: Option<&Surface3>) -> GeomAbsSurfaceType {
    match the_surf {
        Some(Surface3::Plane(_)) => GeomAbsSurfaceType::Plane,
        Some(Surface3::Cylinder(_)) => GeomAbsSurfaceType::Cylinder,
        Some(Surface3::Cone(_)) => GeomAbsSurfaceType::Cone,
        Some(Surface3::Sphere(_)) => GeomAbsSurfaceType::Sphere,
        Some(Surface3::Torus(_)) => GeomAbsSurfaceType::Torus,
        Some(Surface3::Bezier(_)) => GeomAbsSurfaceType::BezierSurface,
        Some(Surface3::BSpline(_)) => GeomAbsSurfaceType::BSplineSurface,
        Some(Surface3::Revolution(_)) => GeomAbsSurfaceType::SurfaceOfRevolution,
        Some(Surface3::LinearExtrusion(_)) => GeomAbsSurfaceType::SurfaceOfExtrusion,
        Some(Surface3::Offset(_)) => GeomAbsSurfaceType::OffsetSurface,
        // OCCT GeomAdaptor_Surface::GetType has no case for these derived
        // surface kinds and falls through to GeomAbs_OtherSurface.
        Some(_) => GeomAbsSurfaceType::OtherSurface,
        None => GeomAbsSurfaceType::OtherSurface,
    }
}

/// OCCT Geom2d_Curve::D1(U, P, V) — the point and first derivative of a 2D
/// curve (arch. diff. #3).
fn curve2d_d1(the_c: &rcad_kernel::geom::Curve2d, the_u: f64) -> (DVec2, DVec2) {
    (the_c.point_at(the_u), the_c.derivative_at(the_u))
}

/// OCCT gp_Vec2d::Angle(gp_Dir2d) — the angle in [-PI, PI] between the two
/// directions (gp_Vec2d.cxx Angle: atan2 of the cross/dot pair).
pub(crate) fn vec2_angle(the_a: DVec2, the_b: DVec2) -> f64 {
    let d1 = the_a.normalize_or_zero();
    let d2 = the_b.normalize_or_zero();
    let dot = (d1.x * d2.x + d1.y * d2.y).clamp(-1.0, 1.0);
    let cross = d1.x * d2.y - d1.y * d2.x;
    cross.atan2(dot)
}
