//! OCCT BRepCheck_Wire (TKTopAlgo/BRepCheck).
//!
//! Source: `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Wire.cxx`
//! (L69-2176) and `BRepCheck_Wire.hxx` (L34-122).
//!
//! GAP summary (see the module doc of `brep_check_analyzer`):
//! - `SelfIntersect` (Wire.cxx L1074-1742) needs `Geom2dInt_GInter`
//!   curve-curve 2D intersection (self-intersection of one pcurve and pairs
//!   of pcurves). The rcad `GInter` translation covers the line-vs-curve
//!   case only, so the intersection runs are neutral (no status set); the
//!   surrounding structure (edge collection, pcurve/table setup, bounding-box
//!   rejection, the NoCurveOnSurface / InvalidRange guards) is ported.

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, SurfaceEval};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, ShapeType, TShape};
use rcad_kernel::math::bnd::BndBox2d;
use std::collections::{HashMap, HashSet};

use crate::topalgo::brep_class::bnd_lib_add2d_curve::add_2d_curve;
use crate::topalgo::shape_source::FaceShapeSource;

use super::brep_check_result::{
    brep_check_add, brep_tool_tolerance_vertex, explorer,
    iterator_subshapes, oriented, ShapeKey, BRepCheckResultBase, BRepCheckStatus,
};

/// OCCT Wire.cxx L95-98: `IsOriented(S)`.
pub fn is_oriented(s: &Shape) -> bool {
    s.orientation == Orientation::Forward || s.orientation == Orientation::Reversed
}

/// OCCT `NCollection_Map<TopoDS_Shape>` (TopTools_ShapeMapHasher) with the
/// insertion order preserved for the GetOrientation key scan.
#[derive(Debug, Default)]
pub struct ShapeSet {
    items: Vec<Shape>,
    index: HashSet<ShapeKey>,
}

impl ShapeSet {
    pub fn new() -> Self {
        ShapeSet::default()
    }

    /// OCCT Map::Add — returns false when already present.
    pub fn add(&mut self, s: &Shape) -> bool {
        if self.index.contains(&ShapeKey::of(s)) {
            false
        } else {
            self.index.insert(ShapeKey::of(s));
            self.items.push(s.clone());
            true
        }
    }

    /// OCCT Map::Contains.
    pub fn contains(&self, s: &Shape) -> bool {
        self.index.contains(&ShapeKey::of(s))
    }

    pub fn clear(&mut self) {
        self.items.clear();
        self.index.clear();
    }

    /// OCCT Map::Remove.
    pub fn remove(&mut self, s: &Shape) {
        if self.index.remove(&ShapeKey::of(s)) {
            self.items.retain(|it| !it.is_equal(s));
        }
    }

    pub fn extent(&self) -> usize {
        self.items.len()
    }

    pub fn items(&self) -> &[Shape] {
        &self.items
    }
}

/// OCCT `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// (myMapVE / myMapEF) — indexed insertion order with IsEqual keys.
#[derive(Debug, Default)]
pub struct IndexedShapeMap {
    keys: Vec<Shape>,
    values: Vec<Vec<Shape>>,
    index: HashMap<ShapeKey, usize>,
}

impl IndexedShapeMap {
    pub fn new() -> Self {
        IndexedShapeMap::default()
    }

    /// OCCT FindIndex — the 1-based OCCT index mapped to Option<usize>.
    pub fn find_index(&self, s: &Shape) -> Option<usize> {
        self.index.get(&ShapeKey::of(s)).copied()
    }

    /// OCCT Add — appends with an empty list; returns the index.
    pub fn add(&mut self, s: &Shape, list: Vec<Shape>) -> usize {
        if let Some(&i) = self.index.get(&ShapeKey::of(s)) {
            self.values[i] = list;
            return i;
        }
        self.index.insert(ShapeKey::of(s), self.keys.len());
        self.keys.push(s.clone());
        self.values.push(list);
        self.keys.len() - 1
    }

    /// OCCT FindIndex == 0 → Add.
    pub fn find_index_or_add(&mut self, s: &Shape) -> usize {
        match self.find_index(s) {
            Some(i) => i,
            None => self.add(s, Vec::new()),
        }
    }

    /// OCCT operator()(i) — the list at index i.
    pub fn value(&self, i: usize) -> &Vec<Shape> {
        &self.values[i]
    }

    pub fn value_mut(&mut self, i: usize) -> &mut Vec<Shape> {
        &mut self.values[i]
    }

    /// OCCT FindKey(i).
    pub fn find_key(&self, i: usize) -> &Shape {
        &self.keys[i]
    }

    /// OCCT Seek — the list of `s` when present.
    pub fn seek(&self, s: &Shape) -> Option<&Vec<Shape>> {
        self.find_index(s).map(|i| &self.values[i])
    }

    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    pub fn clear(&mut self) {
        self.keys.clear();
        self.values.clear();
        self.index.clear();
    }
}

/// OCCT Wire.cxx L431-436: `IsDistanceIn3DTolerance`.
pub fn is_distance_in_3d_tolerance(the_pnt_f: glam::DVec3, the_pnt_l: glam::DVec3, a_tol3d: f64) -> bool {
    let dist = the_pnt_f.distance(the_pnt_l);
    dist < a_tol3d
}

/// OCCT Wire.cxx L440-519: `IsDistanceIn2DTolerance` (static). The
/// `BRepAdaptor_Surface` of the OCCT signature is carried as the face
/// surface value; the domain bounds come from
/// `surface_adaptor_basis_and_bounds` (the adaptor Load unwrap).
pub fn is_distance_in_2d_tolerance(
    a_face_surface: &rcad_kernel::geom::Surface3,
    the_pnt: DVec2,
    the_pnt_ref: DVec2,
    a_tol3d: f64,
) -> bool {
    let (_basis, bounds) =
        rcad_kernel::topods::surface_adaptor_basis_and_bounds(a_face_surface);
    // OCCT L450-451.
    let dumax = 0.01 * (bounds[1] - bounds[0]);
    let dvmax = 0.01 * (bounds[3] - bounds[2]);
    // OCCT L452-453.
    let dumin = (the_pnt.x - the_pnt_ref.x).abs();
    let dvmin = (the_pnt.y - the_pnt_ref.y).abs();

    if dumin < dumax && dvmin < dvmax {
        return true;
    }

    // OCCT L478-479: UResolution / VResolution.
    let mut dumax = a_face_surface.u_resolution(a_tol3d);
    let mut dvmax = a_face_surface.v_resolution(a_tol3d);
    // OCCT L480-484: D1 at the mid point.
    let um = (the_pnt.x + the_pnt_ref.x) / 2.;
    let vm = (the_pnt.y + the_pnt_ref.y) / 2.;
    let (_a_p, a_du, a_dv) = a_face_surface.derivatives(um, vm);
    let a_mdu = a_du.length();
    if a_mdu > rcad_kernel::CONFUSION {
        dumax = (a_tol3d / a_mdu).max(dumax);
    }
    let a_mdv = a_dv.length();
    if a_mdv > rcad_kernel::CONFUSION {
        dvmax = (a_tol3d / a_mdv).max(dvmax);
    }

    // OCCT L505.
    let a_tol2d = 2.0 * dumax.max(dvmax);

    // OCCT L516-518.
    let dist = dumin.max(dvmin);

    dist < a_tol2d
}

/// OCCT BRepCheck_Wire (Wire.hxx L34-122).
#[derive(Debug)]
pub struct BRepCheckWire {
    /// OCCT protected base (Result.hxx L82-90).
    pub base: BRepCheckResultBase,
    /// OCCT myCdone.
    pub my_cdone: bool,
    /// OCCT myCstat.
    pub my_cstat: BRepCheckStatus,
    /// OCCT myMapVE.
    pub my_map_ve: IndexedShapeMap,
    /// OCCT myGctrl.
    pub my_gctrl: bool,
}

impl BRepCheckWire {
    /// OCCT BRepCheck_Wire::BRepCheck_Wire(const TopoDS_Wire& W)
    /// (Wire.cxx L120-126).
    pub fn new(brep: &BRep, w: &Shape) -> Self {
        let mut r = BRepCheckWire {
            base: BRepCheckResultBase::new(),
            my_cdone: false,
            my_cstat: BRepCheckStatus::NoError,
            my_map_ve: IndexedShapeMap::new(),
            my_gctrl: false,
        };
        r.base.init(w);
        r.minimum(brep);
        r
    }

    /// OCCT BRepCheck_Wire::Minimum (Wire.cxx L130-188).
    pub fn minimum(&mut self, brep: &BRep) {
        // OCCT L131-132.
        self.my_cdone = false;
        self.my_gctrl = true;
        if !self.base.my_min {
            // OCCT L136-138.
            self.base.my_map.bound(&self.base.my_shape);
            let my_shape = self.base.my_shape.clone();

            // OCCT L141-160: check that the wire is "connex"; fill myMapVE.
            let edges = explorer(brep, &my_shape, ShapeType::Edge);
            let mut nbedge = 0i32;
            self.my_map_ve.clear();
            for e in &edges {
                nbedge += 1;
                for vtx in explorer(brep, e, ShapeType::Vertex) {
                    let index = self.my_map_ve.find_index_or_add(&vtx);
                    self.my_map_ve.value_mut(index).push(e.clone());
                }
            }
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            if nbedge == 0 {
                // OCCT L162-165.
                brep_check_add(lst, BRepCheckStatus::EmptyWire);
            } else if nbedge >= 2 {
                // OCCT L167-180.
                let mut map_e = ShapeSet::new();
                let first = edges[0].clone();
                propagate(&self.my_map_ve, &first, &mut map_e);
                let mut not_connected = false;
                for e in &edges {
                    if !map_e.contains(e) {
                        brep_check_add(lst, BRepCheckStatus::NotConnected);
                        not_connected = true;
                        break;
                    }
                }
                let _ = not_connected;
            }
            // OCCT L181-186.
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            if lst.is_empty() {
                lst.push(BRepCheckStatus::NoError);
            }
            self.my_map_ve.clear();
            self.base.my_min = true;
        }
    }

    /// OCCT BRepCheck_Wire::InContext (Wire.cxx L192-268).
    pub fn in_context(&mut self, brep: &BRep, s: &Shape) {
        // OCCT L194-209: bound check under the (parallel) lock.
        if self.base.my_map.is_bound(s) {
            return;
        }
        self.base.my_map.bind(s.clone(), Vec::new());

        // OCCT L213-225: check if my wire is in <S>.
        {
            let mut found = false;
            let mut more = false;
            for cur in explorer(brep, s, ShapeType::Wire) {
                more = true;
                if cur.is_same(&self.base.my_shape) {
                    found = true;
                    break;
                }
            }
            if !more || !found {
                let lst = self
                    .base
                    .my_map
                    .find_mut(s)
                    .expect("InContext: the context must be bound");
                brep_check_add(lst, BRepCheckStatus::SubshapeNotInShape);
                return;
            }
        }

        let mut st = BRepCheckStatus::NoError;
        let styp = s.shape_type();
        match styp {
            ShapeType::Face => {
                // OCCT L231-253.
                if self.my_gctrl {
                    let mut ed1 = Shape::null();
                    let mut ed2 = Shape::null();
                    st = self.self_intersect(brep, s, &mut ed1, &mut ed2, true);
                }
                if st != BRepCheckStatus::NoError {
                    // break
                } else {
                    st = self.closed(brep, false);
                    if st != BRepCheckStatus::NoError {
                        // break
                    } else {
                        st = self.orientation(brep, s, false);
                        if st != BRepCheckStatus::NoError {
                            // break
                        } else {
                            st = self.closed2d(brep, s, false);
                        }
                    }
                }
            }
            _ => {
                // OCCT L254-256: default — nothing.
            }
        }

        // OCCT L259-267.
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");
        if st != BRepCheckStatus::NoError {
            brep_check_add(lst, st);
        }
        if lst.is_empty() {
            lst.push(BRepCheckStatus::NoError);
        }
    }

    /// OCCT BRepCheck_Wire::Blind (Wire.cxx L272-279).
    pub fn blind(&mut self) {
        if !self.base.my_blind {
            self.base.my_blind = true;
        }
    }

    /// OCCT BRepCheck_Wire::Closed (Wire.cxx L283-424).
    pub fn closed(&mut self, brep: &BRep, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        if self.my_cdone {
            // OCCT L296-303.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed: myShape must be bound");
                brep_check_add(lst, self.my_cstat);
            }
            return self.my_cstat;
        }

        // OCCT L305.
        self.my_cdone = true;

        // OCCT L307-312: an already stored error is returned as-is.
        {
            let list = self
                .base
                .my_map
                .find(&my_shape)
                .expect("Closed: myShape must be bound");
            if let Some(&first) = list.first() {
                if first != BRepCheckStatus::NoError {
                    self.my_cstat = first;
                    return self.my_cstat; // already saved
                }
            }
        }

        self.my_cstat = BRepCheckStatus::NoError;

        // OCCT L316-319.
        let mut map_s = ShapeSet::new();
        let mut cradoc: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();
        let mut cradoc_order: Vec<Shape> = Vec::new();
        self.my_map_ve.clear();

        // OCCT L323-350: each oriented vertex on oriented edges found 2 times.
        let edges = explorer(brep, &my_shape, ShapeType::Edge);
        for e in &edges {
            if is_oriented(e) {
                let key = ShapeKey::of(e);
                if !cradoc.contains_key(&key) {
                    cradoc.insert(key, Vec::new());
                    cradoc_order.push(e.clone());
                }
                cradoc.get_mut(&key).expect("Cradoc entry").push(e.clone());

                map_s.add(e);
                for vtx in explorer(brep, e, ShapeType::Vertex) {
                    if is_oriented(&vtx) {
                        let index = self.my_map_ve.find_index_or_add(&vtx);
                        self.my_map_ve.value_mut(index).push(e.clone());
                    }
                }
            }
        }

        // OCCT L352-373.
        let the_nbori = map_s.extent();
        if the_nbori >= 2 {
            map_s.clear();
            let mut first: Option<Shape> = None;
            for e in &edges {
                if is_oriented(e) {
                    first = Some(e.clone());
                    break;
                }
            }
            if let Some(e) = &first {
                propagate(&self.my_map_ve, e, &mut map_s);
            }
        }
        if the_nbori != map_s.extent() {
            self.my_cstat = BRepCheckStatus::NotConnected;
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed: myShape must be bound");
                brep_check_add(lst, self.my_cstat);
            }
            return self.my_cstat;
        }

        // OCCT L375-404: the occurrence check (maximum 2, FORWARD + REVERSED).
        let mut yabug = false;
        for key in &cradoc_order {
            let occurrences = &cradoc[&ShapeKey::of(key)];
            if occurrences.len() >= 3 {
                yabug = true;
            } else if occurrences.len() == 2 {
                if occurrences.first().map(|s| s.orientation)
                    == occurrences.last().map(|s| s.orientation)
                {
                    yabug = true;
                }
            }
            if yabug {
                self.my_cstat = BRepCheckStatus::RedundantEdge;
                if update {
                    let lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Closed: myShape must be bound");
                    brep_check_add(lst, self.my_cstat);
                }
                return self.my_cstat;
            }
        }

        // OCCT L406-417.
        for i in 0..self.my_map_ve.extent() {
            if self.my_map_ve.value(i).len() % 2 != 0 {
                self.my_cstat = BRepCheckStatus::NotClosed;
                if update {
                    let lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Closed: myShape must be bound");
                    brep_check_add(lst, self.my_cstat);
                }
                return self.my_cstat;
            }
        }

        // OCCT L419-423.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Closed: myShape must be bound");
            brep_check_add(lst, self.my_cstat);
        }
        self.my_cstat
    }

    /// OCCT BRepCheck_Wire::Closed2d (Wire.cxx L523-715).
    pub fn closed2d(&mut self, brep: &BRep, the_face: &Shape, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        let mut a_closed_stat;

        // OCCT L536-545: 3d closure checked too.
        a_closed_stat = self.closed(brep, false);
        if a_closed_stat != BRepCheckStatus::NoError {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed2d: myShape must be bound");
                brep_check_add(lst, a_closed_stat);
            }
            return a_closed_stat;
        }

        // OCCT L550: BRepAdaptor_Surface aFaceSurface(theFace, false).
        let face_surface = face_surface_adaptor(brep, the_face);
        let Some(face_surface) = face_surface else {
            // No face surface — GAP: the adaptor cannot be built; neutral
            // (no status set) beyond the OCCT NoError carried so far.
            return a_closed_stat;
        };

        // OCCT L560-577: count the oriented edges.
        let mut a_nb_oriented_edges = 0i32;
        for e in explorer(brep, &my_shape, ShapeType::Edge) {
            if is_oriented(&e) {
                a_nb_oriented_edges += 1;
            }
        }

        if a_nb_oriented_edges == 0 {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed2d: myShape must be bound");
                brep_check_add(lst, a_closed_stat);
            }
            return a_closed_stat;
        }

        // OCCT L579-601: the wire explorer over the oriented edges.
        let wire_edges = explorer(brep, &my_shape, ShapeType::Edge);
        let ordered = wire_explorer(brep, the_face, &my_shape, &wire_edges);
        let a_nb_found_edges = ordered.len() as i32;
        let a_first_edge = ordered.first().cloned();
        let a_last_edge = ordered.last().cloned();

        if a_nb_found_edges != a_nb_oriented_edges {
            a_closed_stat = BRepCheckStatus::NotClosed;
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed2d: myShape must be bound");
                brep_check_add(lst, a_closed_stat);
            }
            return a_closed_stat;
        }
        let (a_first_edge, a_first_vertex) = match chain_endpoints(brep, &ordered) {
            Some((start, _head)) => {
                let e = a_first_edge.expect("aFirstEdge");
                (e, start)
            }
            None => {
                // An empty traversal with oriented edges present — GAP
                // (neutral, no status set).
                return a_closed_stat;
            }
        };

        // OCCT L605-647: infinite range checks.
        let mut is_first_infinite = false;
        let mut is_last_infinite = false;

        let mut an_ori = a_first_edge.orientation;
        let mut range = brep.edge_range(&a_first_edge);
        let (mut a_f, mut a_l) = (range[0], range[1]);
        if (an_ori == Orientation::Forward
            && rcad_kernel::precision::is_negative_infinite_value(a_f))
            || (an_ori == Orientation::Reversed
                && rcad_kernel::precision::is_positive_infinite_value(a_l))
        {
            is_first_infinite = true;
        }

        if let Some(le) = &a_last_edge {
            an_ori = le.orientation;
            range = brep.edge_range(le);
            (a_f, a_l) = (range[0], range[1]);
            if (an_ori == Orientation::Forward
                && rcad_kernel::precision::is_positive_infinite_value(a_l))
                || (an_ori == Orientation::Reversed
                    && rcad_kernel::precision::is_negative_infinite_value(a_f))
            {
                is_last_infinite = true;
            }
        }

        if is_first_infinite && is_last_infinite {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed2d: myShape must be bound");
                brep_check_add(lst, a_closed_stat);
            }
            return a_closed_stat;
        } else if a_first_vertex.is_null() {
            a_closed_stat = BRepCheckStatus::NotClosed;
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed2d: myShape must be bound");
                brep_check_add(lst, a_closed_stat);
            }
            return a_closed_stat;
        }

        // OCCT L649-656: get the last point.
        let mut a_p_first = DVec2::ZERO;
        let mut a_p_last = DVec2::ZERO;
        if let Some(le) = &a_last_edge {
            let pc = brep.curve_on_surface(le, the_face);
            if let Some((pc, f, l)) = pc {
                let (p_first_raw, p_last_raw) = (pc.point_at(f), pc.point_at(l));
                // OCCT L652-656: UVPoints(aLastEdge, theFace, aP_temp, aP_last);
                // REVERSED swaps.
                a_p_last = if le.orientation == Orientation::Reversed {
                    p_first_raw
                } else {
                    p_last_raw
                };
            }
        }

        // OCCT L666-673: get the first point.
        let pc1 = brep.curve_on_surface(&a_first_edge, the_face);
        if let Some((pc, f, l)) = pc1 {
            let (p_first_raw, p_last_raw) = (pc.point_at(f), pc.point_at(l));
            if a_first_edge.orientation == Orientation::Reversed {
                a_p_first = p_last_raw;
            } else {
                a_p_first = p_first_raw;
            }
        }

        // OCCT L677-686: the periodic-face seam check.
        if !is_closed_2d_for_periodic_face(brep, the_face, a_p_first, a_p_last, &a_first_vertex) {
            a_closed_stat = BRepCheckStatus::NotClosed;
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed2d: myShape must be bound");
                brep_check_add(lst, a_closed_stat);
            }
            return a_closed_stat;
        }

        // OCCT L694-708: the 2D/3D tolerance distance checks.
        // (BRepTools_WireExplorer::CurrentVertex() after the full traversal is
        // the traversal head — the END of the last ordered edge.)
        let current_vertex = chain_endpoints(brep, &ordered)
            .map(|(_start, head)| head)
            .unwrap_or_else(|| a_first_vertex.clone());
        let a_tol3d = brep_tool_tolerance_vertex(brep, &a_first_vertex)
            .max(brep_tool_tolerance_vertex(brep, &current_vertex));

        let a_pnt_ref = brep.vertex_position(&a_first_vertex);
        let a_pnt = brep.vertex_position(&current_vertex);

        if !is_distance_in_2d_tolerance(&face_surface, a_p_first, a_p_last, a_tol3d) {
            a_closed_stat = BRepCheckStatus::NotClosed;
        }

        if !is_distance_in_3d_tolerance(a_pnt_ref, a_pnt, a_tol3d) {
            a_closed_stat = BRepCheckStatus::NotClosed;
        }

        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Closed2d: myShape must be bound");
            brep_check_add(lst, a_closed_stat);
        }
        a_closed_stat
    }

    /// OCCT BRepCheck_Wire::Orientation (Wire.cxx L719-1070).
    pub fn orientation(&mut self, brep: &BRep, f: &Shape, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        // OCCT L721: theOstat = Closed().
        let mut the_ostat = self.closed(brep, false);

        // OCCT L733-740.
        if the_ostat != BRepCheckStatus::NotClosed && the_ostat != BRepCheckStatus::NoError {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Orientation: myShape must be bound");
                brep_check_add(lst, the_ostat);
            }
            return the_ostat;
        }

        the_ostat = BRepCheckStatus::NoError;

        // OCCT L744-749: VF, VL, ledge, ListOfPassedEdge, mapS, theEdge, theRef.
        let mut vf: Option<Shape> = None;
        let mut vl: Option<Shape> = None;
        let mut map_s = ShapeSet::new();
        let mut the_edge: Option<Shape> = None;
        let mut the_ref: Option<Shape> = None;

        // OCCT L752-783: checks the orientation of the edges.
        for e in explorer(brep, &my_shape, ShapeType::Edge) {
            let edg = e;
            let orient = edg.orientation;
            if is_oriented(&edg) {
                map_s.add(&edg);
                the_edge = Some(edg.clone());
                the_ref = Some(edg.clone());
                for vte in explorer(brep, &edg, ShapeType::Vertex) {
                    let vto = vte.orientation;
                    if vto == Orientation::Forward {
                        vf = Some(vte.clone());
                    } else if vto == Orientation::Reversed {
                        vl = Some(vte.clone());
                    }
                    if vf.is_some() && vl.is_some() {
                        break;
                    }
                }
                if vf.is_none() && vl.is_none() {
                    the_ostat = BRepCheckStatus::InvalidDegeneratedFlag;
                }
                break;
            }
            let _ = orient;
        }

        if the_ostat == BRepCheckStatus::NoError {
            // OCCT L787-1064.
            let mut index: i32 = 1;
            let nb_ori_no_degen = self.my_map_ve.extent() as i32;
            // OCCT L790-796.
            let mut is_go_fwd = true;
            if vl.is_none() {
                is_go_fwd = false;
            }

            while index < nb_ori_no_degen {
                let mut ledge: Vec<Shape> = Vec::new();
                let _list_of_passed_edge: Vec<Shape> = Vec::new();
                // OCCT L805-818: find the chain on VL (or VF).
                let ind: Option<usize> = if let Some(vl_s) = &vl {
                    self.my_map_ve.find_index(vl_s)
                } else if let Some(vf_s) = &vf {
                    self.my_map_ve.find_index(vf_s)
                } else {
                    the_ostat = BRepCheckStatus::InvalidDegeneratedFlag;
                    break;
                };

                let Some(ind) = ind else {
                    // OCCT L814-817 (the else branch): the vertex is absent
                    // from myMapVE — FindIndex returns 0 and the OCCT code
                    // would dereference an empty list; treated as the
                    // degenerated flag path.
                    the_ostat = BRepCheckStatus::InvalidDegeneratedFlag;
                    break;
                };

                let chain: Vec<Shape> = self.my_map_ve.value(ind).clone();
                let mut ortmp = Orientation::Forward;
                for itls in &chain {
                    let edg = itls;
                    let orient = edg.orientation;
                    if map_s.contains(&edg) {
                        ortmp = get_orientation(&map_s, &edg);
                    }

                    // OCCT L830-852: add already passed outcoming edges.
                    if map_s.contains(&edg) && ortmp == orient && !edg.is_same(the_edge.as_ref().expect("theEdge")) {
                        let anchor = if vl.is_some() { vl.as_ref() } else { vf.as_ref() };
                        for vte in explorer(brep, &edg, ShapeType::Vertex) {
                            let vto = vte.orientation;
                            if let Some(anchor) = anchor {
                                if vl.is_some() {
                                    if vto == Orientation::Forward && anchor.is_same(&vte) {
                                        // ListOfPassedEdge.Append(edg)
                                        break;
                                    }
                                } else if vto == Orientation::Reversed && anchor.is_same(&vte) {
                                    // ListOfPassedEdge.Append(edg)
                                    break;
                                }
                            }
                        }
                    }

                    // OCCT L855-887.
                    if !map_s.contains(&edg) || ortmp != orient {
                        for vte in explorer(brep, &edg, ShapeType::Vertex) {
                            let vto = vte.orientation;
                            if let Some(anchor) = if vl.is_some() {
                                vl.as_ref()
                            } else {
                                vf.as_ref()
                            } {
                                if vl.is_some() {
                                    if vto == Orientation::Forward && anchor.is_same(&vte) {
                                        // OCCT L866-869.
                                        if !f.is_null() || !brep.is_edge_degenerated(&edg) {
                                            ledge.push(edg.clone());
                                        }
                                        break;
                                    }
                                } else if vto == Orientation::Reversed && anchor.is_same(&vte) {
                                    // OCCT L879-883.
                                    if !f.is_null() || !brep.is_edge_degenerated(&edg) {
                                        ledge.push(edg.clone());
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                let mut nbconnex = ledge.len() as i32;
                let mut changedesens = false;
                if nbconnex == 0 {
                    // OCCT L891-919.
                    if self.my_cstat == BRepCheckStatus::NotClosed {
                        if vl.is_none() {
                            if update {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(&my_shape)
                                    .expect("Orientation: myShape must be bound");
                                brep_check_add(lst, the_ostat);
                            }
                            return the_ostat; // leave
                        } else {
                            index -= 1; // because after Index++ and no chain
                            vl = None; // chain on VF is forced
                            the_edge = the_ref.clone();
                            changedesens = true;
                        }
                    } else {
                        the_ostat = BRepCheckStatus::BadOrientationOfSubshape;
                        if update {
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("Orientation: myShape must be bound");
                            brep_check_add(lst, the_ostat);
                        }
                        return the_ostat;
                    }
                } else if !f.is_null() {
                    // OCCT L923-948: try to see in 2d.
                    let pivot = if let Some(vl_s) = &vl {
                        vl_s.clone()
                    } else {
                        vf.clone().expect("pivot")
                    };
                    choix_uv(brep, &pivot, the_edge.as_ref().expect("theEdge"), f, &mut ledge);
                    nbconnex = ledge.len() as i32;
                }

                if nbconnex >= 2 {
                    // OCCT L950-958.
                    the_ostat = BRepCheckStatus::BadOrientationOfSubshape;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Orientation: myShape must be bound");
                        brep_check_add(lst, the_ostat);
                    }
                    return the_ostat;
                } else if nbconnex == 1 {
                    // OCCT L959-995: offset the vertex.
                    let first = ledge[0].clone();
                    let mut advanced = false;
                    for vte in explorer(brep, &first, ShapeType::Vertex) {
                        let vto = vte.orientation;
                        if vl.is_some() {
                            if vto == Orientation::Reversed {
                                vl = Some(vte.clone());
                                advanced = true;
                                break;
                            }
                        } else if vto == Orientation::Forward {
                            vf = Some(vte.clone());
                            advanced = true;
                            break;
                        }
                    }
                    map_s.add(&first);
                    the_edge = Some(first);
                    if !advanced {
                        if vl.is_some() {
                            vl = None;
                        } else {
                            vf = None;
                        }
                    }
                } else if !changedesens {
                    // OCCT L996-1004: nbconnex == 0.
                    the_ostat = BRepCheckStatus::NotClosed;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Orientation: myShape must be bound");
                        brep_check_add(lst, the_ostat);
                    }
                    return the_ostat;
                }

                // OCCT L1006-1059: check the closure of the wire in 2d.
                let mut a_v_ref: Option<Shape> = None;
                let mut is_check_close = false;
                if is_go_fwd && vf.is_some() {
                    a_v_ref = vf.clone();
                    is_check_close = true;
                } else if !is_go_fwd && vl.is_some() {
                    a_v_ref = vl.clone();
                    is_check_close = true;
                }

                if index == 1
                    && self.my_cstat != BRepCheckStatus::NotClosed
                    && is_check_close
                    && !f.is_null()
                {
                    ledge.clear();
                    let a_v_ref_s = a_v_ref.clone().expect("aVRef");
                    if let Some(ind2) = self.my_map_ve.find_index(&a_v_ref_s) {
                        let chain2: Vec<Shape> = self.my_map_ve.value(ind2).clone();
                        for itlsh in &chain2 {
                            let edg = itlsh;
                            let _orient = edg.orientation;
                            if !the_ref
                                .as_ref()
                                .expect("theRef")
                                .is_same(&edg)
                            {
                                for vte in explorer(brep, &edg, ShapeType::Vertex) {
                                    let vto = vte.orientation;
                                    if vto == Orientation::Reversed && a_v_ref_s.is_same(&vte) {
                                        ledge.push(edg.clone());
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    choix_uv(
                        brep,
                        &a_v_ref_s,
                        the_ref.as_ref().expect("theRef"),
                        f,
                        &mut ledge,
                    );
                    if ledge.is_empty() {
                        the_ostat = BRepCheckStatus::NotClosed;
                        if update {
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("Orientation: myShape must be bound");
                            brep_check_add(lst, the_ostat);
                        }
                        return the_ostat;
                    }
                }

                index += 1;
            }
        }
        // OCCT L1065-1069.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Orientation: myShape must be bound");
            brep_check_add(lst, the_ostat);
        }
        the_ostat
    }

    /// OCCT BRepCheck_Wire::SelfIntersect (Wire.cxx L1074-1742).
    ///
    /// GAP: the `Geom2dInt_GInter` curve-curve intersection runs
    /// (self-intersection L1172, pairwise L1302) have no rcad translation
    /// with general curve input — the intersection result stays empty and
    /// the surrounding checks run neutrally (no status set); the structural
    /// guards (EmptyWire, NoCurveOnSurface, InvalidRange, box rejection,
    /// same-edge rejection) are ported.
    pub fn self_intersect(
        &mut self,
        brep: &BRep,
        f: &Shape,
        ret_e1: &mut Shape,
        ret_e2: &mut Shape,
        update: bool,
    ) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        let _ = ret_e2;
        // OCCT L1101-1105.
        let _tolint = 1e-10;
        // OCCT L1104-1105: HS = BRepAdaptor_Surface(F, false).
        let hs = face_surface_adaptor(brep, f);
        // OCCT L1107-1115: EMap over the direct wire children.
        let mut e_map = ShapeSet::new();
        for it1 in iterator_subshapes(brep, &my_shape) {
            if it1.shape_type() == ShapeType::Edge {
                e_map.add(&it1);
            }
        }
        let nbedges = e_map.extent();
        if nbedges == 0 {
            // OCCT L1116-1123.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("SelfIntersect: myShape must be bound");
                brep_check_add(lst, BRepCheckStatus::EmptyWire);
            }
            return BRepCheckStatus::EmptyWire;
        }

        // OCCT L1125-1127: tabDom / tabCur / boxes.
        let mut tab_cur: Vec<Option<(Curve2d, f64, f64)>> = vec![None; nbedges];
        let mut tab_dom_first: Vec<Option<DVec2>> = vec![None; nbedges];
        let mut tab_dom_last: Vec<Option<DVec2>> = vec![None; nbedges];
        let mut boxes: Vec<BndBox2d> = (0..nbedges).map(|_| BndBox2d::new()).collect();

        for i in 0..nbedges {
            let e1 = e_map.items()[i].clone();
            if i == 0 {
                // OCCT L1134-1163.
                let pcu = brep.curve_on_surface(&e1, f);
                let Some((pcu, mut first1, mut last1)) = pcu else {
                    // OCCT L1135-1144.
                    *ret_e1 = e1;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("SelfIntersect: myShape must be bound");
                        brep_check_add(lst, BRepCheckStatus::SelfIntersectingWire);
                    }
                    return BRepCheckStatus::SelfIntersectingWire;
                };
                // OCCT L1146-1158: the periodic guard on the adaptor range.
                if !Curve2dEval::is_periodic(&pcu) {
                    let dom = Curve2dEval::default_domain(&pcu);
                    if dom[0] > first1 {
                        first1 = dom[0];
                    }
                    if dom[1] < last1 {
                        last1 = dom[1];
                    }
                }
                // OCCT L1160-1161.
                let uv = uv_points(&brep, &e1, f, &pcu, first1, last1);
                tab_dom_first[0] = Some(uv.0);
                tab_dom_last[0] = Some(uv.1);
                tab_cur[0] = Some((pcu.clone(), first1, last1));
                // OCCT L1163.
                add_2d_curve(&pcu, first1, last1, rcad_kernel::precision::PCONFUSION, &mut boxes[0]);
            }

            // OCCT L1171-1172: Inter.Perform(C1, myDomain1, tolint, tolint) —
            // self-intersection of C1.
            // GAP: Geom2dInt_GInter::Perform(C, D, TolConf, Tol) — no rcad
            // translation for general curves — the result stays empty.

            // OCCT L1242-1302: the pairwise setup for j > i.
            for j in (i + 1)..nbedges {
                let e2 = e_map.items()[j].clone();
                if i == 0 {
                    // OCCT L1247-1281.
                    let pc2 = brep.curve_on_surface(&e2, f);
                    match pc2 {
                        Some((c2, mut first2, mut last2)) if last2 > first2 => {
                            if !Curve2dEval::is_periodic(&c2) {
                                let dom = Curve2dEval::default_domain(&c2);
                                if dom[0] > first2 {
                                    first2 = dom[0];
                                }
                                if dom[1] < last2 {
                                    last2 = dom[1];
                                }
                            }
                            let uv = uv_points(&brep, &e2, f, &c2, first2, last2);
                            tab_dom_first[j] = Some(uv.0);
                            tab_dom_last[j] = Some(uv.1);
                            tab_cur[j] = Some((c2.clone(), first2, last2));
                            add_2d_curve(
                                &c2,
                                first2,
                                last2,
                                rcad_kernel::precision::PCONFUSION,
                                &mut boxes[j],
                            );
                        }
                        Some((c2, _f2, _l2)) => {
                            // OCCT L1269-1281: last2 <= first2 → InvalidRange;
                            // null curve → NoCurveOnSurface.
                            let _ = c2;
                            return BRepCheckStatus::InvalidRange;
                        }
                        None => {
                            return BRepCheckStatus::NoCurveOnSurface;
                        }
                    }
                }

                // OCCT L1288-1296: box rejection and same-edge rejection.
                if boxes[i].is_out_box(&boxes[j]) {
                    continue;
                }
                if e1.is_same(&e2) {
                    continue;
                }

                // OCCT L1302: Inter.Perform(C1, myDomain1, C2, tabDom[j-1],
                // tolint, tolint).
                // GAP: Geom2dInt_GInter::Perform(C1, D1, C2, D2, ...) — no
                // rcad translation for general curves — the intersection
                // points/segments stay empty and the point/segment loops
                // (L1335-1730) do not run.
                let _ = hs;
            }
        }

        // OCCT L1736-1741.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("SelfIntersect: myShape must be bound");
            brep_check_add(lst, BRepCheckStatus::NoError);
        }
        BRepCheckStatus::NoError
    }

    /// OCCT BRepCheck_Wire::SetStatus (Wire.cxx L1746-1749).
    pub fn set_status(&mut self, the_status: BRepCheckStatus) {
        let lst = self
            .base
            .my_map
            .find_mut(&self.base.my_shape)
            .expect("SetStatus: myShape must be bound");
        brep_check_add(lst, the_status);
    }

    /// OCCT BRepCheck_Wire::GeometricControls(B) (Wire.cxx L1753-1763).
    pub fn set_geometric_controls(&mut self, b: bool) {
        if self.my_gctrl != b {
            if b {
                self.my_cdone = false;
            }
            self.my_gctrl = b;
        }
    }

    /// OCCT BRepCheck_Wire::GeometricControls() (Wire.cxx L1767-1770).
    pub fn geometric_controls(&self) -> bool {
        self.my_gctrl
    }
}

/// `BRep_Tool::UVPoints(E, F, PF, PL)` — the pcurve values at the (adjusted)
/// parameter range ends.
fn uv_points(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    pc: &Curve2d,
    first: f64,
    last: f64,
) -> (DVec2, DVec2) {
    let _ = brep;
    let _ = e;
    let _ = f;
    (pc.point_at(first), pc.point_at(last))
}

/// `BRepAdaptor_Surface(theFace, false)` — the face surface value in world
/// coordinates.
fn face_surface_adaptor(brep: &BRep, face: &Shape) -> Option<rcad_kernel::geom::Surface3> {
    brep.face_surface_world(face)
}

/// OCCT BRepTools_WireExplorer (the rcad `order_wire_edges` translation,
/// fclass2d.rs L179) — the wire's edges reordered for traversal along the
/// face pcurves. Returns the edge occurrences in traversal order.
fn wire_explorer(brep: &BRep, face: &Shape, _wire: &Shape, wire_edges: &[Shape]) -> Vec<Shape> {
    // Build the FaceShapeSource-compatible index mapping (the same
    // enumeration order as FaceShapeSource::new).
    let Some(fd) = face.as_face() else {
        return wire_edges.to_vec();
    };
    let surf = match brep.face_surface_world(face) {
        Some(s) => s,
        None => return wire_edges.to_vec(),
    };
    // DS-convention location table: slot 0 = identity.
    let mut locations: Vec<glam::DAffine3> = vec![glam::DAffine3::IDENTITY];
    locations.extend(brep.locations.iter().copied());
    let ds = FaceShapeSource::new(face, surf, &locations);

    // The exact index map replays the FaceShapeSource::new enumeration
    // (indices are 1-based over the edge list).
    let mut source_index: HashMap<ShapeKey, usize> = HashMap::new();
    let mut next = 1usize;
    for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
        if let TShape::Wire(wd) = &*w.data {
            for e in &wd.edges {
                source_index.entry(ShapeKey::of(e)).or_insert(next);
                next += 1;
            }
        }
    }

    let input: Vec<(usize, Orientation)> = wire_edges
        .iter()
        .map(|e| {
            let idx = source_index.get(&ShapeKey::of(e)).copied().unwrap_or(0);
            (idx, e.orientation)
        })
        .collect();
    let ordered = crate::topalgo::brep_top_adaptor::fclass2d::order_wire_edges(
        &ds,
        0,
        &input,
    );
    // Map back from the (source index, orientation) pairs to the wire
    // occurrences (first match).
    let mut used = vec![false; wire_edges.len()];
    let mut out = Vec::with_capacity(wire_edges.len());
    for (idx, ori) in ordered {
        let mut found = None;
        for (k, e) in wire_edges.iter().enumerate() {
            if !used[k]
                && e.orientation == ori
                && source_index.get(&ShapeKey::of(e)).copied() == Some(idx)
            {
                found = Some(k);
                break;
            }
        }
        match found {
            Some(k) => {
                used[k] = true;
                out.push(wire_edges[k].clone());
            }
            None => {
                // The explorer dropped an edge (INTERNAL/EXTERNAL or
                // degenerate traversal) — keep it out, like OCCT.
            }
        }
    }
    out
}

/// The traversal chain endpoints of the ordered wire: (start of the first
/// edge, head after the last edge), derived from the vertex-adjacency chain
/// (the OCCT BRepTools_WireExplorer walks the same adjacency through its
/// myMap; rcad primitive wires carry unreliable per-edge orientation tags, so
/// the endpoints are taken from the IsSame vertex identity chain — for a
/// closed wire both endpoints are the shared closure corner, exactly the
/// OCCT head position).
fn chain_endpoints(brep: &BRep, ordered: &[Shape]) -> Option<(Shape, Shape)> {
    if ordered.is_empty() {
        return None;
    }
    let edge_ends = |e: &Shape| -> Vec<Shape> {
        match e.as_edge() {
            Some(ed) => vec![ed.first.clone(), ed.last.clone()],
            None => Vec::new(),
        }
    };
    let same_key = |a: &Shape, b: &Shape| ShapeKey::of(a) == ShapeKey::of(b);

    let first = &ordered[0];
    let last = ordered[ordered.len() - 1].clone();
    let first_ends = edge_ends(first);
    if first_ends.is_empty() {
        return None;
    }
    let start = if ordered.len() == 1 {
        first_ends[0].clone()
    } else {
        let next_ends = edge_ends(&ordered[1]);
        first_ends
            .iter()
            .find(|v| !next_ends.iter().any(|w| same_key(w, v)))
            .cloned()
            .unwrap_or_else(|| first_ends[0].clone())
    };
    let last_ends = edge_ends(&last);
    let head = if ordered.len() == 1 {
        last_ends.last().cloned().unwrap_or_else(|| first_ends[0].clone())
    } else {
        let prev_ends = edge_ends(&ordered[ordered.len() - 2]);
        last_ends
            .iter()
            .find(|v| !prev_ends.iter().any(|w| same_key(w, v)))
            .cloned()
            .unwrap_or_else(|| last_ends.last().cloned().unwrap_or_else(|| first_ends[0].clone()))
    };
    let _ = brep;
    Some((start, head))
}

/// OCCT Wire.cxx L1777-1822: `Propagate` — fill <mapE> with edges connected
/// to <edg> through vertices contained in <mapVE>.
pub fn propagate(map_ve: &IndexedShapeMap, edg: &Shape, map_e: &mut ShapeSet) {
    let mut current_edges: Vec<Shape> = vec![edg.clone()];

    loop {
        let mut next_edges: Vec<Shape> = Vec::new();
        for edge in &current_edges {
            if !map_e.contains(edge) {
                map_e.add(edge);
            }

            for ex in explorer_of_edge_vertices(edge) {
                if let Some(indv) = map_ve.find_index(&ex) {
                    let edges = map_ve.value(indv);
                    for e in edges.clone() {
                        if !edge.is_same(&e) && !map_e.contains(&e) {
                            map_e.add(&e);
                            next_edges.push(e);
                        }
                    }
                }
            }
        }
        current_edges = next_edges;
        if current_edges.is_empty() {
            break;
        }
    }
}

/// The vertex sub-shapes of an edge (the OCCT TopExp_Explorer over
/// TopAbs_VERTEX — for an edge this is exactly its two stored vertices,
/// enumerated with the TopoDS_Iterator composition).
fn explorer_of_edge_vertices(edge: &Shape) -> Vec<Shape> {
    super::brep_check_result::child_occurrences(edge)
}

/// OCCT Wire.cxx L1826-1839: `GetOrientation` — the orientation stored in
/// <mapE> for <edg>.
pub fn get_orientation(map_e: &ShapeSet, edg: &Shape) -> Orientation {
    for item in map_e.items() {
        if item.is_same(edg) {
            return item.orientation;
        }
    }
    edg.orientation
}

/// OCCT Wire.cxx L1846-2012: `ChoixUV` — for vertex theVertex find the edge
/// along which we should go further; keeps only that edge in theLOfShape.
pub fn choix_uv(
    brep: &BRep,
    the_vertex: &Shape,
    the_edge: &Shape,
    the_face: &Shape,
    the_lof_shape: &mut Vec<Shape>,
) {
    // OCCT L1851-1862: remove theEdge occurrences.
    the_lof_shape.retain(|s| !the_edge.is_same(s));

    let a_tol3d = brep_tool_tolerance_vertex(brep, the_vertex);

    let mut an_index: i32 = 0;
    let mut an_ind_min: i32 = 0;
    let mut a_pnt_ref = DVec2::ZERO;
    let mut a_pnt = DVec2::ZERO;
    let mut a_min_angle = f64::MAX;
    let mut a_max_angle = f64::MIN;
    let a_gp_resolution = rcad_kernel::math::gp::GP_RESOLUTION;
    let mut a_param: f64 = 0.0;
    #[allow(unused_assignments)]
    let mut a_par_piv: f64 = 0.0;

    // OCCT L1874: BRepAdaptor_Surface aFaceSurface(theFace, false).
    let Some(a_face_surface) = face_surface_adaptor(brep, the_face) else {
        // GAP: the face adaptor cannot be built — the selection is neutral
        // (the candidate list is left untouched).
        return;
    };

    // OCCT L1876-1881: C2d over theEdge.
    let Some((c2d_edge, a_first_param, a_last_param)) = brep.curve_on_surface(the_edge, the_face)
    else {
        // OCCT L1878-1881: JAG 10.12.96.
        return;
    };

    let a_v_orientation = the_vertex.orientation;
    let an_edg_orientation = the_edge.orientation;

    // OCCT L1886.
    a_par_piv = if a_v_orientation == an_edg_orientation {
        a_first_param
    } else {
        a_last_param
    };

    // OCCT L1890: CurveDirForParameter(C2d, aParPiv, aPntRef, aDerRef).
    let mut a_der_ref = DVec2::ZERO;
    curve_dir_for_parameter_2d(&c2d_edge, a_par_piv, &mut a_pnt_ref, &mut a_der_ref);

    // OCCT L1892-1895.
    if a_v_orientation != an_edg_orientation {
        a_der_ref = DVec2::new(-a_der_ref.x, -a_der_ref.y);
    }

    let mut candidates: Vec<(usize, Shape)> = the_lof_shape.iter().cloned().enumerate().collect();
    for (pos, an_e) in candidates.iter_mut() {
        let _ = pos;
        an_index += 1;
        let an_e = an_e.clone();
        // OCCT L1903-1907.
        let Some((c2d_e, fp_e, lp_e)) = brep.curve_on_surface(&an_e, the_face) else {
            continue;
        };
        let _ = c2d_e.clone();

        // OCCT L1910-1911.
        a_param = if a_v_orientation != an_e.orientation {
            fp_e
        } else {
            lp_e
        };
        a_pnt = c2d_e.point_at(a_param);

        // OCCT L1913-1916.
        if !is_distance_in_2d_tolerance(&a_face_surface, a_pnt, a_pnt_ref, a_tol3d) {
            continue;
        }

        // OCCT L1918.
        let mut a_der = DVec2::ZERO;
        curve_dir_for_parameter_2d(&c2d_e, a_param, &mut a_pnt, &mut a_der);

        // OCCT L1920-1923.
        if a_v_orientation == an_e.orientation {
            a_der = DVec2::new(-a_der.x, -a_der.y);
        }

        // OCCT L1925-1929.
        if a_der_ref.length() <= a_gp_resolution || a_der.length() <= a_gp_resolution {
            continue;
        }

        // OCCT L1931-1936.
        let mut an_angle = -angle_2d(a_der_ref, a_der);
        if an_angle < 0. {
            an_angle += 2. * std::f64::consts::PI;
        }

        // OCCT L1938-1953.
        if the_face.orientation == Orientation::Forward {
            if an_angle < a_min_angle {
                an_ind_min = an_index;
                a_min_angle = an_angle;
            }
        } else if an_angle > a_max_angle {
            an_ind_min = an_index;
            a_max_angle = an_angle;
        }
    }
    let _ = &candidates;

    // OCCT L1956-2011: update the edge list.
    if an_ind_min == 0 {
        if the_lof_shape.len() == 1 {
            let is_found_initial = true;
            let an_e_found = the_lof_shape[0].clone();
            #[allow(unused_assignments)]
            let mut is_found = is_found_initial;

            if an_e_found.is_null()
                || brep.is_edge_degenerated(the_edge)
                || brep.is_edge_degenerated(&an_e_found)
            {
                is_found = false;
            } else if !is_distance_in_2d_tolerance(&a_face_surface, a_pnt, a_pnt_ref, a_tol3d) {
                is_found = false;
            } else {
                // OCCT L1975-1982: closure in 3D — the BRepAdaptor_Curve(E, F)
                // values are the surface images of the pcurves.
                let p_edg = brep
                    .curve_on_surface(the_edge, the_face)
                    .and_then(|(pc, _f, _l)| {
                        brep.face_surface_world(the_face)
                            .map(|s| s.point_at(pc.point_at(a_par_piv).x, pc.point_at(a_par_piv).y))
                    });
                let p_e_found = brep
                    .curve_on_surface(&an_e_found, the_face)
                    .and_then(|(pc, _f, _l)| {
                        brep.face_surface_world(the_face)
                            .map(|s| s.point_at(pc.point_at(a_param).x, pc.point_at(a_param).y))
                    });
                match (p_edg, p_e_found) {
                    (Some(p1), Some(p2)) => {
                        is_found = is_distance_in_3d_tolerance(p1, p2, a_tol3d);
                    }
                    _ => {
                        is_found = false;
                    }
                }
            }

            if !is_found {
                the_lof_shape.clear();
            }
        } else {
            the_lof_shape.clear();
        }
    } else {
        // OCCT L1994-2011: keep only the selected edge.
        let keep = (an_ind_min - 1) as usize;
        let selected = if keep < the_lof_shape.len() {
            Some(the_lof_shape[keep].clone())
        } else {
            None
        };
        the_lof_shape.clear();
        if let Some(s) = selected {
            the_lof_shape.push(s);
        }
    }
}

/// Signed angle of two 2D vectors (OCCT gp_Vec2d::Angle).
fn angle_2d(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// OCCT Wire.cxx L2016-2037: `CurveDirForParameter`.
pub fn curve_dir_for_parameter_2d(a_c2d: &Curve2d, a_prm: f64, pnt: &mut DVec2, a_vec2d: &mut DVec2) {
    let a_tol = rcad_kernel::math::gp::GP_RESOLUTION;

    *pnt = a_c2d.point_at(a_prm);
    *a_vec2d = a_c2d.derivative_at(a_prm);
    //
    if a_vec2d.length() <= a_tol {
        // GAP: the OCCT loop evaluates the higher derivatives DN(aPrm, i)
        // for i in 2..=100; the rcad Curve2dEval carries derivatives only up
        // to the third order — the fallback uses those orders.
        for i in 2..=3 {
            *a_vec2d = match i {
                2 => a_c2d.derivative2_at(a_prm),
                _ => a_c2d.derivative3_at(a_prm),
            };
            if a_vec2d.length() > a_tol {
                break;
            }
        }
    }
}

/// OCCT Wire.cxx L2048-2078: `GetPnt2d` — the parametric point of theVertex
/// on theFace through theEdge's pcurve.
pub fn get_pnt2d(brep: &BRep, the_vertex: &Shape, the_edge: &Shape, the_face: &Shape) -> Option<DVec2> {
    let Some(ed) = the_edge.as_edge() else {
        return None;
    };
    let a_first_vtx = &ed.first;
    let a_last_vtx = &ed.last;

    if !the_vertex.is_same(a_first_vtx) && !the_vertex.is_same(a_last_vtx) {
        return None;
    }

    let (_a_pc, _a_f_par, _a_l_par) = brep.curve_on_surface(the_edge, the_face)?;

    let a_par_on_edge = ed.vertex_params.get(&the_vertex.ptr_id()).copied()?;

    let (pc, _f, _l) = brep.curve_on_surface(the_edge, the_face)?;
    Some(pc.point_at(a_par_on_edge))
}

/// OCCT Wire.cxx L2085-2174: `IsClosed2dForPeriodicFace` — the 2D distance
/// check between the wire ends for periodic faces with a seam edge.
pub fn is_closed_2d_for_periodic_face(
    brep: &BRep,
    the_face: &Shape,
    the_p1: DVec2,
    the_p2: DVec2,
    the_vertex: &Shape,
) -> bool {
    // OCCT L2092-2117: searching for the seam edges.
    let mut a_seam_edges: Vec<Shape> = Vec::new();
    let mut not_seams = ShapeSet::new();
    let mut closed_edges = ShapeSet::new();

    for an_exp in explorer(brep, the_face, ShapeType::Edge) {
        let an_edge = an_exp;

        if not_seams.contains(&an_edge) {
            continue;
        }

        if !is_oriented(&an_edge) || !brep.is_edge_closed_on_face(&an_edge, the_face) {
            not_seams.add(&an_edge);
            continue;
        }

        if !closed_edges.add(&an_edge) {
            a_seam_edges.push(an_edge);
        }
    }

    if a_seam_edges.is_empty() {
        return true;
    }

    // OCCT L2123-2129: the vicinity tolerance.
    let Some(a_face_surface) = face_surface_adaptor(brep, the_face) else {
        return true;
    };
    let a_tol = brep_tool_tolerance_vertex(brep, the_vertex);
    let a_uresol = a_face_surface.u_resolution(a_tol);
    let a_vresol = a_face_surface.v_resolution(a_tol);
    let a_vicinity = (a_uresol * a_uresol + a_vresol * a_vresol).sqrt();
    let a_dist_p1_p2 = the_p1.distance(the_p2);

    for a_seam_edge in a_seam_edges.clone() {
        for a_vtx in explorer(brep, &a_seam_edge, ShapeType::Vertex) {
            // OCCT L2142-2169.
            if is_oriented(&a_vtx) && a_vtx.is_same(the_vertex) {
                let Some(a_pnt1) = get_pnt2d(brep, the_vertex, &a_seam_edge, the_face) else {
                    continue;
                };
                let a_seam_edge_rev = oriented(&a_seam_edge, a_seam_edge.orientation);
                let a_seam_edge_rev = {
                    // OCCT L2155: aSeamEdge = TopoDS::Edge(aSeamEdge.Reversed()).
                    let mut r = a_seam_edge_rev.clone();
                    r.orientation = match r.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        o => o,
                    };
                    r
                };
                let Some(a_pnt2) = get_pnt2d(brep, the_vertex, &a_seam_edge_rev, the_face) else {
                    continue;
                };

                let mut a2d_tol = a_pnt1.distance(a_pnt2) * 1e-2;
                a2d_tol = a2d_tol.max(a_vicinity);

                if a_dist_p1_p2 > a2d_tol {
                    return false;
                }
            }
        }
    }

    true
}
