// OCCT BOPAlgo_WireSplitter (BOPAlgo_WireSplitter.cxx L91-226,
// BOPAlgo_WireSplitter_1.cxx L113-1150, BOPAlgo_WireSplitter.lxx).
//
// Builds loops (wires) from a set of edges lying on a face.
//
// Translation notes:
// - TopoDS_Edge/Vertex/Face/Wire -> rcad Shape; a wire is Vec<Shape> (edge order).
// - TopTools_ShapeMapHasher identity -> Shape::is_same (ptr_id + location).
// - BRepAdaptor_Surface -> SurfaceAdaptor (analytic UResolution/VResolution).
// - Geom2dInt_GInter (RefineAngle2D 2D curve intersection) not yet ported:
//   RefineAngle2D returns false and the iCntInt==2 fallback of RefineAngles is kept.

use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use glam::{DVec2, DVec3};
use indexmap::IndexMap;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, TShape};
use std::collections::{HashMap, HashSet};

/// OCCT BOPAlgo_EdgeInfo (BOPAlgo_WireSplitter.lxx) — per (vertex, edge) record.
#[derive(Clone)]
pub(crate) struct EdgeInfo {
    pub edge: Shape,
    pub passed: bool,
    pub in_flag: bool,
    pub is_inside: bool,
    pub angle: f64,
}

impl EdgeInfo {
    pub fn new() -> Self {
        EdgeInfo {
            edge: Shape::null(),
            passed: false,
            in_flag: false,
            is_inside: false,
            angle: -1.0,
        }
    }
    pub fn set_edge(&mut self, the_e: &Shape) {
        self.edge = the_e.clone();
    }
    pub fn edge(&self) -> &Shape {
        &self.edge
    }
    pub fn set_passed(&mut self, the_flag: bool) {
        self.passed = the_flag;
    }
    pub fn passed(&self) -> bool {
        self.passed
    }
    pub fn set_in_flag(&mut self, the_flag: bool) {
        self.in_flag = the_flag;
    }
    pub fn is_in(&self) -> bool {
        self.in_flag
    }
    pub fn set_angle(&mut self, the_angle: f64) {
        self.angle = the_angle;
    }
    pub fn angle(&self) -> f64 {
        self.angle
    }
    pub fn is_inside(&self) -> bool {
        self.is_inside
    }
    pub fn set_is_inside(&mut self, the_is_inside: bool) {
        self.is_inside = the_is_inside;
    }
}

/// OCCT BOPTools_ConnexityBlock — a connected group of edges with the loops
/// produced from it.
#[derive(Clone)]
pub(crate) struct ConnexityBlock {
    pub shapes: Vec<Shape>,
    pub loops: Vec<Vec<Shape>>,
    pub is_regular: bool,
}

impl ConnexityBlock {
    pub fn new() -> Self {
        ConnexityBlock {
            shapes: Vec::new(),
            loops: Vec::new(),
            is_regular: true,
        }
    }
}

/// OCCT BOPTools_AlgoTools::MakeConnexityBlocks (BOPTools_AlgoTools.cxx
/// L187-256) — groups edges into connected blocks by shared vertices and
/// marks a block irregular when any edge is repeated or any vertex has more
/// than two incident edges.
pub(crate) fn make_connexity_blocks(edges: &[Shape]) -> Vec<ConnexityBlock> {
    // aMFence: dedup start elements; aMNRegular: repeated (multi-connexity) edges.
    let mut a_mfence: HashSet<u64> = HashSet::new();
    let mut a_mn_regular: HashSet<u64> = HashSet::new();
    let mut a_c_start: Vec<Shape> = Vec::new();
    for a_s in edges {
        if a_mfence.insert(a_s.ptr_id()) {
            a_c_start.push(a_s.clone());
        } else {
            a_mn_regular.insert(a_s.ptr_id());
        }
    }
    // Map vertices to incident edges (MapShapesAndAncestors).
    let mut a_c_map: HashMap<u64, Vec<u64>> = HashMap::new(); // vertex ptr -> edge ptrs
    let mut a_edge_ptr_to_idx: HashMap<u64, usize> = HashMap::new(); // edge ptr -> index in a_c_start
    for (ei, e) in a_c_start.iter().enumerate() {
        a_edge_ptr_to_idx.insert(e.ptr_id(), ei);
        for v in edge_vertices(e) {
            let l = a_c_map.entry(v.ptr_id()).or_default();
            if !l.contains(&e.ptr_id()) {
                l.push(e.ptr_id());
            }
        }
    }
    // BFS blocks over edges via shared vertices.
    let n = a_c_start.len();
    let mut a_mfence2: HashSet<u64> = HashSet::new();
    let mut a_blocks: Vec<Vec<usize>> = Vec::new();
    for s in 0..n {
        if !a_mfence2.insert(a_c_start[s].ptr_id()) {
            continue;
        }
        let mut a_l_block: Vec<usize> = vec![s];
        let mut i = 0;
        while i < a_l_block.len() {
            let ei = a_l_block[i];
            for v in edge_vertices(&a_c_start[ei]) {
                if let Some(l) = a_c_map.get(&v.ptr_id()) {
                    for &ep in l {
                        if let Some(&eidx) = a_edge_ptr_to_idx.get(&ep) {
                            if a_mfence2.insert(a_c_start[eidx].ptr_id()) {
                                a_l_block.push(eidx);
                            }
                        }
                    }
                }
            }
            i += 1;
        }
        a_blocks.push(a_l_block);
    }
    // Build ConnexityBlocks.
    let mut result: Vec<ConnexityBlock> = Vec::new();
    for block in a_blocks {
        let mut a_cb = ConnexityBlock::new();
        let mut b_regular = true;
        for &bi in &block {
            let a_s = a_c_start[bi].clone();
            if a_mn_regular.contains(&a_s.ptr_id()) {
                b_regular = false;
                let mut f = a_s.clone();
                f.orientation = Orientation::Forward;
                a_cb.shapes.push(f);
                let mut r = a_s.clone();
                r.orientation = Orientation::Reversed;
                a_cb.shapes.push(r);
            } else {
                a_cb.shapes.push(a_s.clone());
                if b_regular {
                    // Check no multi-connected vertices on this edge.
                    for v in edge_vertices(&a_s) {
                        let cnt = a_c_map.get(&v.ptr_id()).map_or(0, |l| l.len());
                        if cnt != 2 {
                            b_regular = false;
                            break;
                        }
                    }
                }
            }
        }
        a_cb.is_regular = b_regular;
        result.push(a_cb);
    }
    result
}

/// OCCT BOPAlgo_WireSplitter::Perform (L91-118) + MakeWires (L164-226).
/// Returns the loops as edge sequences.
pub(crate) fn split_into_wires(
    face: &Shape,
    face_index: usize,
    edges: &[Shape],
    ds: &DS,
) -> Vec<Vec<Shape>> {
    if edges.is_empty() {
        return Vec::new();
    }
    let mut my_lcb = make_connexity_blocks(edges);
    let mut a_vcb: Vec<ConnexityBlock> = Vec::new();
    let mut result: Vec<Vec<Shape>> = Vec::new();
    for cb in my_lcb.iter_mut() {
        if cb.is_regular {
            let a_w = make_wire(&cb.shapes);
            result.push(a_w);
        } else {
            a_vcb.push(cb.clone());
        }
    }
    let a_context = IntToolsContext::new();
    for mut cb in a_vcb {
        split_block(face, face_index, &mut cb, &a_context, ds);
        for l in &cb.loops {
            result.push(l.clone());
        }
    }
    result
}

/// OCCT BOPAlgo_WireSplitter::MakeWire (BOPAlgo_WireSplitter.lxx) — the wire
/// is the list of edges in order.
fn make_wire(a_le: &[Shape]) -> Vec<Shape> {
    a_le.to_vec()
}

/// OCCT BOPAlgo_WireSplitter::SplitBlock (BOPAlgo_WireSplitter_1.cxx L113-355).
fn split_block(
    face: &Shape,
    face_index: usize,
    cb: &mut ConnexityBlock,
    the_context: &IntToolsContext,
    ds: &DS,
) {
    let my_edges = cb.shapes.clone();
    // mySmartMap: vertex ptr -> (vertex shape, list of EdgeInfo).
    let mut my_smart_map: IndexMap<u64, (Shape, Vec<EdgeInfo>)> = IndexMap::new();
    // aVertMap: vertex ptr -> closed flag.
    let mut a_vert_map: HashMap<u64, bool> = HashMap::new();
    // aMS: edge identity fence with the odd/closed semantics of OCCT.
    let mut a_ms: HashSet<u64> = HashSet::new();
    let mut a_v1 = Shape::null();

    // 1. Fill mySmartMap.
    for a_e in &my_edges {
        if !has_curve_on_surface(a_e, face_index, ds) {
            continue;
        }
        let mut b_is_closed = is_degenerated(a_e) || is_closed_on_face(a_e, face);
        // OCCT L149: if (!aMS.Add(aE) && !bIsClosed) aMS.Remove(aE);
        if !a_ms.insert(a_e.ptr_id()) && !b_is_closed {
            a_ms.remove(&a_e.ptr_id());
        }
        let verts = edge_vertices(a_e);
        for (i, a_v) in verts.iter().enumerate() {
            let vptr = a_v.ptr_id();
            let entry = my_smart_map.entry(vptr).or_insert_with(|| (a_v.clone(), Vec::new()));
            let leix = entry.1.len();
            let mut a_ei = EdgeInfo::new();
            a_ei.set_edge(a_e);
            let a_or = a_v.orientation;
            let b_is_in = a_or == Orientation::Reversed;
            a_ei.set_in_flag(b_is_in);
            entry.1.push(a_ei);
            if i == 0 {
                a_v1 = a_v.clone();
            } else {
                b_is_closed = b_is_closed || a_v1.is_same(a_v);
            }
            if a_vert_map.contains_key(&vptr) {
                if b_is_closed {
                    a_vert_map.insert(vptr, b_is_closed);
                }
            } else {
                a_vert_map.insert(vptr, b_is_closed);
            }
            let _ = leix;
        }
    }
    let a_nb = my_smart_map.len();

    // 2. bNothingToDo check: every vertex has exactly one In + one Out edge.
    let mut b_nothing_to_do = true;
    for (_, leinfo) in my_smart_map.values() {
        let mut a_cnt_in = 0;
        let mut a_cnt_out = 0;
        for a_ei in leinfo {
            if a_ei.is_in() {
                a_cnt_in += 1;
            } else {
                a_cnt_out += 1;
            }
        }
        if a_cnt_in != 1 || a_cnt_out != 1 {
            b_nothing_to_do = false;
            break;
        }
    }
    // Second part of bNothingToDo: check edges on TShape coincidence.
    if b_nothing_to_do {
        let mut a_map_ee: HashMap<u64, Vec<u64>> = HashMap::new();
        for a_e in &my_edges {
            a_map_ee.entry(a_e.ptr_id()).or_default().push(a_e.ptr_id());
        }
        let mut b_flag = true;
        for (_, a_l_ex) in &a_map_ee {
            let a_nb_e = a_l_ex.len();
            if a_nb_e == 1 {
                continue;
            } else if a_nb_e == 2 {
                // Both entries are the same edge identity -> not a single wire.
                b_flag = false;
                break;
            } else {
                b_flag = false;
                break;
            }
        }
        b_nothing_to_do = b_nothing_to_do && b_flag;
    }
    if b_nothing_to_do {
        let a_w = make_wire(&cb.shapes);
        cb.loops.push(a_w);
        return;
    }

    // 3. Angles in mySmartMap.
    let a_gas = SurfaceAdaptor::from_face(face);
    for i in 0..a_nb {
        // OCCT uses ChangeEdges() — a mutable reference into mySmartMap; the
        // angles/inside flags set here must persist for the Path building.
        let (a_v_sh, mut leinfo) = my_smart_map.get_index(i).unwrap().1.clone();
        for a_ei in leinfo.iter_mut() {
            let a_e = a_ei.edge().clone();
            a_ei.set_is_inside(!a_ms.contains(&a_e.ptr_id()));
            let mut a_vv = a_v_sh.clone();
            let b_is_in = a_ei.is_in();
            let a_or = if b_is_in { Orientation::Reversed } else { Orientation::Forward };
            a_vv.orientation = a_or;
            let a_angle = angle_2d(&a_vv, &a_e, face_index, &a_gas, b_is_in, the_context, ds);
            a_ei.set_angle(a_angle);
        }
        if let Some((_, val)) = my_smart_map.get_index_mut(i) {
            val.1 = leinfo;
        }
    }

    // OCCT L324: RefineAngles(myFace, mySmartMap, theContext).
    refine_angles(&mut my_smart_map);

    // 4. Path building.
    for i in 0..a_nb {
        let (a_va_sh, leinfo) = my_smart_map.get_index(i).unwrap().1.clone();
        for a_ei in leinfo.iter() {
            let a_ei_mut = a_ei.clone();
            let a_e_outa = a_ei_mut.edge().clone();
            let b_is_out = !a_ei_mut.is_in();
            // OCCT reads the live EdgeInfo from mySmartMap; Path() marks edges
            // passed, so re-read the current flag instead of the stale clone.
            let cur_passed = my_smart_map
                .get_index(i)
                .unwrap()
                .1
                .1
                .iter()
                .find(|e| e.edge().is_same(&a_e_outa) && e.is_in() == a_ei_mut.is_in())
                .map(|e| e.passed())
                .unwrap_or(true);
            let b_is_not_passed = !cur_passed;
            if b_is_out && b_is_not_passed {
                let mut a_ls: Vec<Shape> = Vec::new();
                let mut a_vert_va: Vec<Shape> = Vec::new();
                let mut a_coord_va: Vec<DVec2> = Vec::new();
                path(
                    &a_gas,
                    face_index,
                    &a_vert_map,
                    &a_va_sh,
                    &a_e_outa,
                    a_ei_mut,
                    &mut a_ls,
                    &mut a_vert_va,
                    &mut a_coord_va,
                    cb,
                    &mut my_smart_map,
                    ds,
                );
            }
        }
    }
}

/// OCCT Path (BOPAlgo_WireSplitter_1.cxx L359-627) — walks edges from a
/// starting (vertex, out-edge) and closes loops by the smallest clockwise angle.
#[allow(clippy::too_many_arguments)]
fn path(
    a_gas: &SurfaceAdaptor,
    face_index: usize,
    a_vert_map: &HashMap<u64, bool>,
    a_v_first: &Shape,
    a_e_first: &Shape,
    a_ei_first: EdgeInfo,
    a_ls: &mut Vec<Shape>,
    a_vert_va: &mut Vec<Shape>,
    a_coord_va: &mut Vec<DVec2>,
    cb: &mut ConnexityBlock,
    my_smart_map: &mut IndexMap<u64, (Shape, Vec<EdgeInfo>)>,
    ds: &DS,
) {
    let mut a_va = a_v_first.clone();
    let mut a_e_outa = a_e_first.clone();
    let mut an_edge_info = a_ei_first;
    let a_two_pi = std::f64::consts::PI * 2.0;
    let eps = f64::EPSILON;
    let mut an_info_seq: Vec<EdgeInfo> = Vec::new();

    loop {
        // Do not escape through the edge from which you enter.
        let a_nb = a_ls.len();
        if a_nb == 1 {
            if a_ls[a_nb - 1].is_same(&a_e_outa) {
                return;
            }
        }
        an_edge_info.set_passed(true);
        // OCCT L406: anEdgeInfo->SetPassed(true) — the edge info is a reference
        // INTO mySmartMap, so the flag must be written back to the map, not
        // only tracked locally; the outer loop reads the live flag.
        mark_edge_passed(my_smart_map, &a_va, &a_e_outa, an_edge_info.is_in());
        a_ls.push(a_e_outa.clone());
        a_vert_va.push(a_va.clone());
        an_info_seq.push(an_edge_info.clone());

        let mut p_va = a_va.clone();
        p_va.orientation = Orientation::Forward;
        let a_pa = coord_2d(&p_va, &a_e_outa, face_index, ds);
        a_coord_va.push(a_pa);

        let a_vb = get_next_vertex(&p_va, &a_e_outa);
        let a_pb = coord_2d(&a_vb, &a_e_outa, face_index, ds);

        let a_lei_opt = my_smart_map.get(&a_vb.ptr_id());
        let a_lei: Vec<EdgeInfo> = a_lei_opt.map(|(_, l)| l.clone()).unwrap_or_default();

        let a_tol2d = 2.0 * tolerance_2d(&a_vb, a_gas);
        let a_tol2d2 = a_tol2d * a_tol2d;
        let b_is_closed = a_vert_map.get(&a_vb.ptr_id()).copied().unwrap_or(false);

        // Scan back for the loop-closing vertex.
        let mut a_buf: Vec<Shape> = Vec::new();
        let mut b_has_edge = false;
        let mut closed_break = false;
        {
            let a_nb2 = a_ls.len();
            let mut i = a_nb2;
            while i > 0 {
                let a_v_prev = a_vert_va[i - 1].clone();
                let a_pa_prev = a_coord_va[i - 1];
                let a_e_prev = a_ls[i - 1].clone();
                a_buf.push(a_e_prev.clone());
                if !b_has_edge {
                    b_has_edge = !is_degenerated(&a_e_prev);
                    if !b_has_edge {
                        i -= 1;
                        continue;
                    }
                }
                let an_is_same_v = a_v_prev.is_same(&a_vb);
                let mut an_is_same_v2d = an_is_same_v;
                if an_is_same_v {
                    if b_is_closed {
                        let a_d2 = a_pa_prev.distance_squared(a_pb);
                        an_is_same_v2d = a_d2 < a_tol2d2;
                        if an_is_same_v2d {
                            let udist = (a_pa_prev.x - a_pb.x).abs();
                            let vdist = (a_pa_prev.y - a_pb.y).abs();
                            let a_tol_u = 2.0 * u_tolerance_2d(&a_vb, a_gas);
                            let a_tol_v = 2.0 * v_tolerance_2d(&a_vb, a_gas);
                            if udist > a_tol_u || vdist > a_tol_v {
                                an_is_same_v2d = false;
                            }
                        }
                    }
                }
                if an_is_same_v && an_is_same_v2d {
                    let mut i_priz = 1;
                    if a_buf.len() == 2 {
                        if a_buf[0].is_same(&a_buf[1]) {
                            i_priz = 0;
                        }
                    }
                    if i_priz != 0 {
                        let a_w = make_wire(&a_buf);
                        cb.loops.push(a_w);
                    }
                    let a_nbj = i - 1;
                    if a_nbj < 1 {
                        a_ls.clear();
                        a_vert_va.clear();
                        a_coord_va.clear();
                        return;
                    }
                    // Truncate sequences to the first aNbj entries.
                    a_ls.truncate(a_nbj);
                    a_vert_va.truncate(a_nbj);
                    a_coord_va.truncate(a_nbj);
                    an_info_seq.truncate(a_nbj);
                    let _ = a_vb.clone();
                    a_e_outa = a_ls[a_nbj - 1].clone();
                    an_edge_info = an_info_seq[a_nbj - 1].clone();
                    closed_break = true;
                    break;
                }
                i -= 1;
            }
        }
        if closed_break {
            continue;
        }

        // Select the next edge.
        let an_angle_in = angle_in(&a_e_outa, &a_lei);
        let mut a_min_angle = 100.0;
        let i_cnt = nb_ways_out(&a_lei);
        let is_boundary = !an_edge_info.is_inside();
        let mut a_nb_ways_inside = 0;
        let mut p_only_way_in: Option<EdgeInfo> = None;
        let mut p_edge_info: Option<EdgeInfo> = None;

        for an_ei in &a_lei {
            let mut an_ei_mut = an_ei.clone();
            let a_e = an_ei_mut.edge().clone();
            let an_is_out = !an_ei_mut.is_in();
            let an_is_not_passed = !an_ei_mut.passed();
            if an_is_out && an_is_not_passed {
                if i_cnt == 0 {
                    return; // no way to go
                }
                if i_cnt == 1 {
                    p_edge_info = Some(an_ei_mut.clone());
                    break;
                }
                let mut an_angle;
                if a_e.is_same(&a_e_outa) {
                    an_angle = a_two_pi;
                } else {
                    if b_is_closed {
                        let a_p2dx = coord_2d_vf(&a_e, face_index, ds);
                        let a_d2 = a_p2dx.distance_squared(a_pb);
                        if a_d2 > a_tol2d2 {
                            continue;
                        }
                    }
                    let an_angle_out = an_ei_mut.angle();
                    an_angle = clock_wise_angle(an_angle_in, an_angle_out);
                }
                if is_boundary && an_ei_mut.is_inside() {
                    a_nb_ways_inside += 1;
                    p_only_way_in = Some(an_ei_mut.clone());
                }
                if an_angle < a_min_angle - eps {
                    a_min_angle = an_angle;
                    p_edge_info = Some(an_ei_mut.clone());
                }
            }
        }
        if a_nb_ways_inside == 1 {
            if let Some(owi) = p_only_way_in {
                p_edge_info = Some(owi);
            }
        }
        let pe = match p_edge_info {
            Some(e) => e,
            None => return,
        };
        a_va = a_vb;
        a_e_outa = pe.edge().clone();
        an_edge_info = pe.clone();
    }
}

/// Mark the EdgeInfo for (vertex, edge, in-flag) as passed in mySmartMap.
/// OCCT Path (L406) sets the flag on a reference into the map.
fn mark_edge_passed(
    my_smart_map: &mut IndexMap<u64, (Shape, Vec<EdgeInfo>)>,
    a_v: &Shape,
    a_e: &Shape,
    in_flag: bool,
) {
    if let Some((_, leinfo)) = my_smart_map.get_mut(&a_v.ptr_id()) {
        if let Some(ei) = leinfo
            .iter_mut()
            .find(|ei| ei.edge().is_same(a_e) && ei.is_in() == in_flag)
        {
            ei.set_passed(true);
        }
    }
}

/// OCCT RefineAngles (BOPAlgo_WireSplitter_1.cxx L930-1054) — refines the
/// angles of the internal (section) edges at a vertex where exactly two
/// boundary edges meet, pushing the internal edges out of the boundary wedge.
/// The geometric refinement (RefineAngle2D, L1058-1149, Geom2dInt_GInter) is
/// not ported yet; the iCntInt==2 fallback of OCCT L1021-1025 is kept.
fn refine_angles(my_smart_map: &mut IndexMap<u64, (Shape, Vec<EdgeInfo>)>) {
    let a_nb = my_smart_map.len();
    for i in 0..a_nb {
        let (_, leinfo) = &mut my_smart_map[i];
        // aA1 = angle of the outgoing boundary edge; aA2 = incoming boundary.
        let mut a_a1 = 0.0;
        let mut a_a2 = 0.0;
        let mut i_cnt_bnd = 0usize;
        let mut i_cnt_int = 0usize;
        for a_ei in leinfo.iter() {
            let b_is_in = a_ei.is_in();
            let a_a = a_ei.angle();
            if !a_ei.is_inside() {
                i_cnt_bnd += 1;
                if !b_is_in {
                    a_a1 = a_a;
                } else {
                    a_a2 = a_a;
                }
            } else {
                i_cnt_int += 1;
            }
        }
        if i_cnt_bnd != 2 {
            continue;
        }
        let a_delta = clock_wise_angle(a_a2, a_a1);
        // Refine the internal OUT edges (edge ptr -> new angle, aDMSR).
        let mut a_dmsr: HashMap<u64, f64> = HashMap::new();
        for a_ei in leinfo.iter() {
            let a_e = a_ei.edge();
            let b_is_boundary = !a_ei.is_inside();
            let b_is_in = a_ei.is_in();
            if b_is_boundary || b_is_in {
                continue;
            }
            let a_a = a_ei.angle();
            let a_da = clock_wise_angle(a_a2, a_a);
            if a_da < a_delta {
                continue; // already inside
            }
            // bRefined = RefineAngle2D(aV, aE, myFace, aA1, aA2, aDelta, aA, ctx)
            // — not ported (Geom2dInt_GInter). The iCntInt==2 fallback remains.
            let b_refined = false;
            let a_a_new;
            if b_refined {
                a_a_new = a_a;
            } else if i_cnt_int == 2 {
                a_a_new = if a_a <= a_a1 {
                    a_a1 + rcad_kernel::core::precision::ANGULAR
                } else {
                    a_a2 - rcad_kernel::core::precision::ANGULAR
                };
            } else {
                continue;
            }
            a_dmsr.insert(a_e.ptr_id(), a_a_new);
        }
        if a_dmsr.is_empty() {
            continue;
        }
        // OCCT L1033-1053: update the angles.
        for a_ei in leinfo.iter_mut() {
            let a_e = a_ei.edge();
            let b_is_in = a_ei.is_in();
            let a_a = match a_dmsr.get(&a_e.ptr_id()) {
                Some(v) => *v,
                None => continue,
            };
            let mut a_a = a_a;
            if b_is_in {
                a_a = a_a + std::f64::consts::PI;
            }
            a_ei.set_angle(a_a);
        }
    }
}

/// OCCT ClockWiseAngle (L631-669).
fn clock_wise_angle(a_angle_in: f64, a_angle_out: f64) -> f64 {
    let a_two_pi = std::f64::consts::PI * 2.0;
    let mut a_in = a_angle_in;
    let mut a_out = a_angle_out;
    if a_in >= a_two_pi {
        a_in -= a_two_pi;
    }
    if a_out >= a_two_pi {
        a_out -= a_two_pi;
    }
    let a1 = a_in + std::f64::consts::PI;
    let mut a1 = if a1 >= a_two_pi { a1 - a_two_pi } else { a1 };
    let a2 = a_out;
    let mut d_a = a1 - a2;
    if d_a <= 0.0 {
        d_a = a_two_pi + d_a;
    } else if d_a <= 1e-14 {
        d_a = a_two_pi;
    }
    let _ = a1;
    d_a
}

/// OCCT Coord2d (L673-684) — 2D parameter of a vertex on an edge in a face.
fn coord_2d(a_v1: &Shape, a_e1: &Shape, face_index: usize, ds: &DS) -> DVec2 {
    let a_t = match vertex_param_on_edge(a_v1, a_e1) {
        Some(t) => t,
        None => return DVec2::ZERO,
    };
    match edge_pcurve(a_e1, face_index, ds) {
        Some((c, _, _)) => curve2d_point(&c, a_t),
        None => DVec2::ZERO,
    }
}

/// OCCT Coord2dVf (L688-707) — 2D coord of the FORWARD endpoint of an edge.
fn coord_2d_vf(a_e: &Shape, face_index: usize, ds: &DS) -> DVec2 {
    let a_coord = 99.0;
    for v in edge_vertices(a_e) {
        if v.orientation == Orientation::Forward {
            return coord_2d(&v, a_e, face_index, ds);
        }
    }
    DVec2::new(a_coord, a_coord)
}

/// OCCT NbWaysOut (L711-730).
fn nb_ways_out(a_lei_info: &[EdgeInfo]) -> usize {
    let mut i_cnt = 0;
    for an_ei in a_lei_info {
        let b_is_out = !an_ei.is_in();
        let b_is_not_passed = !an_ei.passed();
        if b_is_out && b_is_not_passed {
            i_cnt += 1;
        }
    }
    i_cnt
}

/// OCCT AngleIn (L734-755) — angle of the incoming edge at this vertex.
fn angle_in(a_e_in: &Shape, a_lei_info: &[EdgeInfo]) -> f64 {
    for an_edge_info in a_lei_info {
        let a_e = an_edge_info.edge();
        let an_is_in = an_edge_info.is_in();
        // OCCT L747: aE == aEIn — TopoDS_Shape::operator== is IsEqual
        // (orientation-sensitive); a Reversed copy of the same edge does NOT
        // match (BOPAlgo_WireSplitter_1.cxx L747).
        if an_is_in && a_e.is_equal(a_e_in) {
            return an_edge_info.angle();
        }
    }
    0.0
}

/// OCCT GetNextVertex (L759-774).
fn get_next_vertex(a_v: &Shape, a_e: &Shape) -> Shape {
    for v in edge_vertices(a_e) {
        if !v.is_equal(a_v) {
            return v;
        }
    }
    a_v.clone()
}

/// OCCT Angle2D (L778-865) — tangent angle of the edge at the vertex, in the
/// face's UV space. The p-curve direction is taken toward the vertex and
/// reversed when the vertex is the edge's REVERSED endpoint (bIsIN).
fn angle_2d(
    a_v: &Shape,
    an_edge: &Shape,
    face_index: usize,
    a_gas: &SurfaceAdaptor,
    b_is_in: bool,
    _the_context: &IntToolsContext,
    ds: &DS,
) -> f64 {
    let a_tv = match vertex_param_on_edge(a_v, an_edge) {
        Some(t) => t,
        None => return 0.0,
    };
    let (a_c2d, a_first, a_last) = match edge_pcurve(an_edge, face_index, ds) {
        Some(v) => v,
        None => return 0.0,
    };
    let tol2d = 2.0 * tolerance_2d(a_v, a_gas);
    let mut dt = curve2d_resolution(&a_c2d, tol2d).max(rcad_kernel::core::precision::PCONFUSION);

    // OCCT L820-833: for non-line curves adjust dt by curvature (GeomLProp_CLProps2d).
    // Pending: curvature-based dt adjustment for non-line p-curves.

    // OCCT L835-845: aTX = 0.05*(aLast-aFirst).
    let mut a_tx = 0.05 * (a_last - a_first);
    if a_tx < 5e-5 {
        a_tx = 5e-5f64.min((a_last - a_first) / 2.0);
    }
    if dt > a_tx {
        dt = a_tx;
    }

    let a_tv1 = if (a_tv - a_first).abs() < (a_tv - a_last).abs() {
        a_tv + dt
    } else {
        a_tv - dt
    };
    let a_pv1 = curve2d_point(&a_c2d, a_tv1);
    let a_pv = curve2d_point(&a_c2d, a_tv);
    let a_v2d = if b_is_in { a_pv1 - a_pv } else { a_pv - a_pv1 };
    let a_dir2d = a_v2d.normalize_or_zero();
    angle_from_dir(a_dir2d)
}

/// OCCT Angle (L869-880) — angle of a 2D direction against +X, in [0, 2π).
fn angle_from_dir(a_dir2d: DVec2) -> f64 {
    let mut an_angle = a_dir2d.x.atan2(a_dir2d.y).mul_add(-1.0, 0.0);
    // OCCT: gp_Dir2d::Angle(aRefDir) — angle from reference dir (+X) to aDir2d.
    an_angle = (a_dir2d.y).atan2(a_dir2d.x);
    if an_angle < 0.0 {
        an_angle += std::f64::consts::PI * 2.0;
    }
    an_angle
}

/// OCCT Tolerance2D (L884-906).
fn tolerance_2d(a_v: &Shape, a_gas: &SurfaceAdaptor) -> f64 {
    let a_tol_v3d = vertex_tolerance(a_v);
    let an_ur = a_gas.u_resolution(a_tol_v3d);
    let a_vr = a_gas.v_resolution(a_tol_v3d);
    let mut a_tol2d = if a_vr > an_ur { a_vr } else { an_ur };
    if a_tol2d < a_tol_v3d {
        a_tol2d = a_tol_v3d;
    }
    if a_gas.get_type() == SurfaceType::Bspline {
        a_tol2d *= 1.1;
    }
    a_tol2d
}

/// OCCT UTolerance2D (L910-916).
fn u_tolerance_2d(a_v: &Shape, a_gas: &SurfaceAdaptor) -> f64 {
    let a_tol_v3d = vertex_tolerance(a_v);
    a_gas.u_resolution(a_tol_v3d)
}

/// OCCT VTolerance2D (L920-926).
fn v_tolerance_2d(a_v: &Shape, a_gas: &SurfaceAdaptor) -> f64 {
    let a_tol_v3d = vertex_tolerance(a_v);
    a_gas.v_resolution(a_tol_v3d)
}

/// OCCT BRepAdaptor_Surface (simplified) — analytic U/V resolution.
struct SurfaceAdaptor {
    surf: Surface3,
    radius: f64,
    uv_domain: Option<[f64; 4]>,
}

#[derive(PartialEq, Clone, Copy)]
enum SurfaceType {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Bspline,
    Other,
}

impl SurfaceAdaptor {
    fn from_face(face: &Shape) -> Self {
        let (surf, uv_domain) = match &*face.data {
            TShape::Face(fd) => (fd.surface.clone().unwrap_or(Surface3::Plane(
                rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z),
            )), fd.uv_domain),
            _ => (
                Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
                None,
            ),
        };
        let radius = match &surf {
            Surface3::Cylinder(c) => c.radius,
            Surface3::Sphere(s) => s.radius,
            Surface3::Cone(c) => c.radius,
            _ => 1.0,
        };
        SurfaceAdaptor {
            surf,
            radius,
            uv_domain,
        }
    }
    fn get_type(&self) -> SurfaceType {
        match &self.surf {
            Surface3::Plane(_) => SurfaceType::Plane,
            Surface3::Cylinder(_) => SurfaceType::Cylinder,
            Surface3::Sphere(_) => SurfaceType::Sphere,
            Surface3::Cone(_) => SurfaceType::Cone,
            Surface3::Torus(_) => SurfaceType::Torus,
            Surface3::BSpline(_) | Surface3::Bezier(_) => SurfaceType::Bspline,
            _ => SurfaceType::Other,
        }
    }
    fn u_resolution(&self, tol: f64) -> f64 {
        match &self.surf {
            Surface3::Plane(_) => tol,
            Surface3::Cylinder(_) | Surface3::Sphere(_) => tol / self.radius.max(1e-12),
            Surface3::Cone(_) => tol / self.radius.max(1e-12),
            _ => {
                // fallback: tol / uv extent (i_walking.rs convention)
                match self.uv_domain {
                    Some(d) if (d[1] - d[0]).abs() > 1e-12 => tol / (d[1] - d[0]).abs(),
                    _ => tol,
                }
            }
        }
    }
    fn v_resolution(&self, tol: f64) -> f64 {
        match &self.surf {
            Surface3::Plane(_) => tol,
            Surface3::Cylinder(_) => tol,
            Surface3::Sphere(_) => tol / self.radius.max(1e-12),
            Surface3::Cone(_) => tol / self.radius.max(1e-12),
            _ => {
                match self.uv_domain {
                    Some(d) if (d[3] - d[2]).abs() > 1e-12 => tol / (d[3] - d[2]).abs(),
                    _ => tol,
                }
            }
        }
    }
}

/// OCCT Geom2dAdaptor_Curve::Resolution(tol) — approximate parameter step for
/// a given 2D tolerance.
fn curve2d_resolution(c: &Curve2d, tol: f64) -> f64 {
    match c {
        Curve2d::Line(_) => tol,
        Curve2d::Circle(cir) => tol / cir.radius.max(1e-12),
        Curve2d::Ellipse(el) => tol / el.minor_radius.max(1e-12),
        Curve2d::BSpline(b) => {
            let extent = match (b.knots.first(), b.knots.last()) {
                (Some(f), Some(l)) => (l - f).abs(),
                _ => 1.0,
            };
            if extent > 1e-12 {
                tol / extent.max(1e-12)
            } else {
                tol
            }
        }
        _ => tol,
    }
}

// ---- data access helpers ----

fn edge_vertices(e: &Shape) -> [Shape; 2] {
    match &*e.data {
        // OCCT TopoDS_Iterator(aE) with cumOri=true (default) composes the edge
        // orientation into the vertices (TopoDS_Iterator.cxx L35-37, L72-80):
        // the stored TEdge nodes [V1(Fwd), V2(Rev)] become [V1(Rev), V2(Fwd)]
        // for a REVERSED edge — same order, orientations composed.
        TShape::Edge(ed) => {
            if e.orientation == Orientation::Reversed {
                [Shape::new(ed.first.data.clone(), ed.first.location, Orientation::Reversed),
                 Shape::new(ed.last.data.clone(), ed.last.location, Orientation::Forward)]
            } else {
                [Shape::new(ed.first.data.clone(), ed.first.location, Orientation::Forward),
                 Shape::new(ed.last.data.clone(), ed.last.location, Orientation::Reversed)]
            }
        }
        _ => [Shape::null(), Shape::null()],
    }
}

fn edge_pcurve(e: &Shape, face_index: usize, ds: &DS) -> Option<(Curve2d, f64, f64)> {
    match &*e.data {
        // OCCT BRep_Tool::CurveOnSurface(aE, aF) — keyed by the face identity.
        // rcad keys the edge pcurve map by the DS face index.
        TShape::Edge(ed) => {
            if let Some(v) = ed.pcurves.get(&face_index) {
                return Some(v.clone());
            }
            // make_pcurves inserts pcurves through Arc::make_mut, which clones
            // shared edges: the DS shape receives the pcurve while the face-wire
            // edge (same logical edge) does not. OCCT mutates edge TShapes in
            // place, so the pcurve is visible from either reference; here fall
            // back to the DS canonical shape.
            if let Some(idx) = ds.map_shape_index.get(&(e.ptr_id(), e.location)) {
                if let Some(ed2) = ds.shape_info(*idx).shape.as_edge() {
                    if let Some(v) = ed2.pcurves.get(&face_index) {
                        return Some(v.clone());
                    }
                }
            }
            None
        }
        _ => None,
    }
}

fn has_curve_on_surface(e: &Shape, face_index: usize, ds: &DS) -> bool {
    match &*e.data {
        // OCCT BRep_Tool::CurveOnSurface(aE, aF) — keyed by the face identity.
        // rcad keys the edge pcurve map by the DS face index (make_pcurves uses
        // FaceInfo::Index()), so the DS index — not the Shape's BRep `index` —
        // is the matching key. Same DS-canonical fallback as edge_pcurve.
        TShape::Edge(ed) => {
            if ed.pcurves.contains_key(&face_index) {
                return true;
            }
            if let Some(idx) = ds.map_shape_index.get(&(e.ptr_id(), e.location)) {
                if let Some(ed2) = ds.shape_info(*idx).shape.as_edge() {
                    return ed2.pcurves.contains_key(&face_index);
                }
            }
            false
        }
        _ => false,
    }
}

fn is_degenerated(e: &Shape) -> bool {
    match &*e.data {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::IsClosed(aE, aF) — the edge has a CurveOnClosedSurface
/// representation for the face (seam edge).
fn is_closed_on_face(e: &Shape, face: &Shape) -> bool {
    let f_index = face.index;
    match &*e.data {
        TShape::Edge(ed) => ed.representations.iter().any(|r| {
            matches!(
                r,
                rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, .. }
                    if *face == f_index
            )
        }),
        _ => false,
    }
}

/// OCCT BRep_Tool::Parameter(aV, aE, aF) — vertex parameter on the edge.
fn vertex_param_on_edge(v: &Shape, e: &Shape) -> Option<f64> {
    match &*e.data {
        TShape::Edge(ed) => {
            if let Some(t) = ed.vertex_params.get(&v.index) {
                return Some(*t);
            }
            // DS-created edges (push_edge) have empty vertex_params; the param is
            // only stored on the source BRep edges. OCCT reads the stored param;
            // computing it from the geometry here is the semantic equivalent.
            let curve = ed.curve.as_ref()?;
            let p = match &*v.data {
                TShape::Vertex(vd) => vd.point,
                _ => return None,
            };
            Some(crate::bop::closest_point_on_curve(curve, p).0)
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Tolerance(aV).
fn vertex_tolerance(v: &Shape) -> f64 {
    match &*v.data {
        TShape::Vertex(vd) => vd.tolerance,
        _ => 0.0,
    }
}

fn curve2d_point(c: &Curve2d, t: f64) -> DVec2 {
    Curve2dEval::point_at(c, t)
}
