// OCCT BRepClass_FaceExplorer (BRepClass_FaceExplorer.cxx / .hxx)
// Explores the wires/edges of a face in 2D for the ray-casting classifier.
//
// Segment/OtherSegment build a probing ray through the query point using a
// point on each boundary pcurve (BRepClass_FaceExplorer.cxx L111-280).

use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, Line2d};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, TShape};

use crate::bop::ds::DS;
use crate::topalgo::brep_class::edge::ClassEdge;

// OCCT BRepClass_FaceExplorer.cxx L30-32.
const PROBING_START: f64 = 0.123;
const PROBING_END: f64 = 0.7;
const PROBING_STEP: f64 = 0.2111;

/// OCCT BRepClass_FaceExplorer — explores a face's wires and edges in 2D.
pub struct FaceExplorer {
    face: usize,
    /// Ordered boundary edges (DS index + orientation), outer wire first.
    face_edges: Vec<(usize, Orientation)>,
    /// Per-wire edge grouping for InitWires/InitEdges iteration.
    wires: Vec<Vec<(usize, Orientation)>>,
    // Probing state (OCCT myCurEdgeInd / myCurEdgePar).
    cur_edge_ind: usize,
    cur_edge_par: f64,
    // Wire/edge iteration state.
    current_wire: usize,
    current_edge: usize,
    // Face UV bounds for CheckPoint (OCCT myUMin/myUMax/...).
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    bounds_computed: bool,
    max_tolerance: f64,
    use_bnd_box: bool,
}

impl FaceExplorer {
    /// OCCT BRepClass_FaceExplorer(F) — initialize from the face.
    pub fn new(ds: &DS, face: usize) -> Self {
        let mut fe = FaceExplorer {
            face,
            face_edges: Vec::new(),
            wires: Vec::new(),
            cur_edge_ind: 1,
            cur_edge_par: PROBING_START,
            current_wire: 0,
            current_edge: 0,
            u_min: f64::INFINITY,
            u_max: f64::NEG_INFINITY,
            v_min: f64::INFINITY,
            v_max: f64::NEG_INFINITY,
            bounds_computed: false,
            max_tolerance: 0.1,
            use_bnd_box: false,
        };
        fe.build_edges(ds);
        fe
    }

    /// Build the ordered boundary edges from the face's wires.
    fn build_edges(&mut self, ds: &DS) {
        let face_data = match &*ds.shapes[self.face].shape.data {
            TShape::Face(fd) => fd,
            _ => return,
        };
        let wires: Vec<&Shape> =
            std::iter::once(&face_data.outer_wire).chain(face_data.inner_wires.iter()).collect();
        for w in wires {
            let Some(&wi) = ds.map_shape_index.get(&(w.ptr_id(), w.location)) else {
                continue;
            };
            if wi >= ds.nb_shapes() {
                continue;
            }
            let wire_edges = match &*ds.shapes[wi].shape.data {
                TShape::Wire(wd) => wd.edges.clone(),
                _ => Vec::new(),
            };
            let mut wire_list: Vec<(usize, Orientation)> = Vec::new();
            for eshape in &wire_edges {
                if let Some(&ei) = ds.map_shape_index.get(&(eshape.ptr_id(), eshape.location)) {
                    wire_list.push((ei, eshape.orientation));
                }
            }
            self.wires.push(wire_list.clone());
            self.face_edges.extend(wire_list);
        }
    }

    pub fn set_max_tolerance(&mut self, t: f64) {
        self.max_tolerance = t;
    }

    pub fn set_use_bnd_box(&mut self, b: bool) {
        self.use_bnd_box = b;
    }

    /// OCCT BRepClass_FaceExplorer::CheckPoint (L67-100) — adjust the point if
    /// it is too far from the face's UV bounding box. Returns the (possibly
    /// adjusted) point.
    pub fn check_point(&mut self, ds: &DS, mut point: DVec2) -> DVec2 {
        if !self.bounds_computed {
            self.compute_face_bounds(ds);
        }
        if self.u_min > self.u_max
            || self.u_min.is_infinite()
            || self.u_max.is_infinite()
            || self.v_min.is_infinite()
            || self.v_max.is_infinite()
        {
            return point;
        }
        let a_center = DVec2::new((self.u_min + self.u_max) * 0.5, (self.v_min + self.v_max) * 0.5);
        let a_distance = a_center.distance(point);
        if a_distance.is_infinite() {
            return DVec2::new(
                self.u_min - (self.u_max - self.u_min),
                self.v_min - (self.v_max - self.v_min),
            );
        }
        let an_epsilon = f64::EPSILON * a_distance.abs();
        if an_epsilon > (self.u_max - self.u_min).max(self.v_max - self.v_min) {
            let a_lin_vec = point - a_center;
            let a_lin_dir = a_lin_vec.normalize_or_zero();
            point = a_center + a_lin_dir * (2.0 * an_epsilon);
        }
        point
    }

    /// Compute the face's UV bounds from the boundary pcurves.
    fn compute_face_bounds(&mut self, ds: &DS) {
        self.bounds_computed = true;
        for &(ei, _) in &self.face_edges {
            if ei >= ds.nb_shapes() {
                continue;
            }
            let edge_data = match &*ds.shapes[ei].shape.data {
                TShape::Edge(ed) => ed,
                _ => continue,
            };
            let Some((c2d, f, l)) = edge_data.pcurves.get(&self.face).cloned() else {
                continue;
            };
            const N: usize = 16;
            for i in 0..=N {
                let t = f + (l - f) * (i as f64) / (N as f64);
                let p = c2d.point_at(t);
                self.u_min = self.u_min.min(p.x);
                self.u_max = self.u_max.max(p.x);
                self.v_min = self.v_min.min(p.y);
                self.v_max = self.v_max.max(p.y);
            }
        }
    }

    /// OCCT BRepClass_FaceExplorer::Reject — always false.
    pub fn reject(&self, _p: DVec2) -> bool {
        false
    }

    /// OCCT BRepClass_FaceExplorer::Segment (L111-117) — reset probing state
    /// and find a valid ray through `p`.
    pub fn segment(&mut self, ds: &DS, p: DVec2) -> Option<(Line2d, f64)> {
        self.cur_edge_ind = 1;
        self.cur_edge_par = PROBING_START;
        self.other_segment(ds, p)
    }

    /// OCCT BRepClass_FaceExplorer::OtherSegment (L121-280) — find the next
    /// probing ray through `p` using a point on a boundary pcurve.
    pub fn other_segment(&mut self, ds: &DS, p: DVec2) -> Option<(Line2d, f64)> {
        let a_tol_par_conf2 = rcad_kernel::PCONFUSION * rcad_kernel::PCONFUSION;
        let n_edges = self.face_edges.len();
        while self.cur_edge_ind <= n_edges {
            let (ei, ori) = self.face_edges[self.cur_edge_ind - 1];
            if ori != Orientation::Forward && ori != Orientation::Reversed {
                self.cur_edge_ind += 1;
                self.cur_edge_par = PROBING_START;
                continue;
            }
            if ei >= ds.nb_shapes() {
                self.cur_edge_ind += 1;
                self.cur_edge_par = PROBING_START;
                continue;
            }
            let edge_data = match &*ds.shapes[ei].shape.data {
                TShape::Edge(ed) => ed,
                _ => {
                    self.cur_edge_ind += 1;
                    self.cur_edge_par = PROBING_START;
                    continue;
                }
            };
            let Some((a_c2d, mut a_f_par, mut a_l_par)) =
                edge_data.pcurves.get(&self.face).cloned()
            else {
                // rcad: no pcurve — this edge cannot provide a probing point.
                self.cur_edge_ind += 1;
                self.cur_edge_par = PROBING_START;
                continue;
            };

            // OCCT L150-166: infinite-range normalization.
            if a_f_par.is_infinite() && a_f_par.is_sign_negative() {
                if a_l_par.is_infinite() && a_l_par.is_sign_positive() {
                    a_f_par = -1.0;
                    a_l_par = 1.0;
                } else {
                    a_f_par = a_l_par - 1.0;
                }
            } else if a_l_par.is_infinite() && a_l_par.is_sign_positive() {
                a_l_par = a_f_par + 1.0;
            }

            let mut found = false;
            while self.cur_edge_par < PROBING_END {
                let a_param_in = self.cur_edge_par * a_f_par + (1.0 - self.cur_edge_par) * a_l_par;
                let a_p_on_c = a_c2d.point_at(a_param_in);
                let a_tan_vec = a_c2d.derivative_at(a_param_in);
                let mut par = a_p_on_c.distance_squared(p);
                if par > a_tol_par_conf2 {
                    let a_lin_vec = a_p_on_c - p;
                    let a_lin_dir = a_lin_vec.normalize_or_zero();
                    let a_tan_mod = a_tan_vec.length_squared();
                    if a_tan_mod < a_tol_par_conf2 {
                        self.cur_edge_par += PROBING_STEP;
                        continue;
                    }
                    let a_tan_norm = a_tan_vec / a_tan_mod.sqrt();
                    let a_sin_a = a_tan_norm.x * a_lin_dir.y - a_tan_norm.y * a_lin_dir.x;
                    const SMALL_ANGLE: f64 = 0.001;
                    let is_small_angle = a_sin_a.abs() < SMALL_ANGLE;
                    if is_small_angle {
                        if self.cur_edge_par + PROBING_STEP < PROBING_END {
                            self.cur_edge_par += PROBING_STEP;
                            continue;
                        }
                    }
                    let line = Line2d::new(p, a_lin_dir);

                    // Check that the curve's ends do not lie on the line
                    // (OCCT L207-264).
                    let a_f_p_on_c = a_c2d.point_at(a_f_par);
                    let a_l_p_on_c = a_c2d.point_at(a_l_par);
                    if line.distance(a_f_p_on_c) > a_tol_par_conf2.sqrt()
                        && line.distance(a_l_p_on_c) > a_tol_par_conf2.sqrt()
                    {
                        if is_small_angle {
                            // Small-angle fallback: use the closest point on the
                            // curve to P (Geom2dAPI_ProjectPointOnCurve).
                            let ext = rcad_kernel::base::extrema::ExtPC2d::new(
                                p,
                                &a_c2d,
                                rcad_kernel::PCONFUSION,
                                a_f_par,
                                a_l_par,
                            );
                            if ext.is_done() && ext.nb_ext() > 0 {
                                let a_p_on_c = ext.point(1).point;
                                let a_min_dist = ext.square_distance(1).sqrt();
                                // OCCT compares against the endpoint distances
                                // and picks the min.
                                let a_f_dist = p.distance(a_f_p_on_c);
                                let a_l_dist = p.distance(a_l_p_on_c);
                                let mut best = a_p_on_c;
                                let mut best_d = a_min_dist;
                                if a_f_dist < best_d {
                                    best_d = a_f_dist;
                                    best = a_f_p_on_c;
                                }
                                if a_l_dist < best_d {
                                    best_d = a_l_dist;
                                    best = a_l_p_on_c;
                                }
                                if best_d < par.sqrt() {
                                    par = best_d;
                                    if par < a_tol_par_conf2.sqrt() {
                                        self.cur_edge_par += PROBING_STEP;
                                        continue;
                                    }
                                    let vec = best - p;
                                    let dir = vec.normalize_or_zero();
                                    let line2 = Line2d::new(p, dir);
                                    self.cur_edge_par += PROBING_STEP;
                                    if self.cur_edge_par >= PROBING_END {
                                        self.cur_edge_ind += 1;
                                        self.cur_edge_par = PROBING_START;
                                    }
                                    return Some((line2, par));
                                }
                            }
                        }
                        self.cur_edge_par += PROBING_STEP;
                        if self.cur_edge_par >= PROBING_END {
                            self.cur_edge_ind += 1;
                            self.cur_edge_par = PROBING_START;
                        }
                        found = true;
                        return Some((line, par.sqrt()));
                    }
                }
                self.cur_edge_par += PROBING_STEP;
            }
            if found {
                return Some((Line2d::new(p, DVec2::X), f64::MAX));
            }
            // This curve is not valid for line construction: next edge.
            self.cur_edge_ind += 1;
            self.cur_edge_par = PROBING_START;
        }
        // Nothing found. OCCT sets the horizontal line and returns false —
        // the caller terminates the probing loop.
        None
    }

    /// Wire iteration (OCCT InitWires/MoreWires/NextWire/RejectWire).
    pub fn init_wires(&mut self) {
        self.current_wire = 0;
    }
    pub fn more_wires(&self) -> bool {
        self.current_wire < self.wires.len()
    }
    pub fn next_wire(&mut self) {
        self.current_wire += 1;
    }
    pub fn reject_wire(&self, _l: DVec2, _par: f64) -> bool {
        false
    }

    /// Edge iteration within the current wire.
    pub fn init_edges(&mut self) {
        self.current_edge = 0;
    }
    pub fn more_edges(&self) -> bool {
        if self.current_wire < self.wires.len() {
            self.current_edge < self.wires[self.current_wire].len()
        } else {
            false
        }
    }
    pub fn next_edge(&mut self) {
        self.current_edge += 1;
    }
    pub fn reject_edge(&self, _l: DVec2, _par: f64) -> bool {
        false
    }

    /// OCCT BRepClass_FaceExplorer::CurrentEdge(E, Or) — the current edge and
    /// its orientation, with the next-edge and tolerance flags set.
    pub fn current_edge(&self, ds: &DS) -> Option<ClassEdge> {
        if self.current_wire < self.wires.len()
            && self.current_edge < self.wires[self.current_wire].len()
        {
            let (ei, ori) = self.wires[self.current_wire][self.current_edge];
            let mut ce = ClassEdge::new(ei, self.face);
            ce.set_max_tolerance(self.max_tolerance);
            ce.set_use_bnd_box(self.use_bnd_box);
            ce.set_next_edge(ds);
            // The classifier needs the edge orientation; store it via the
            // edge's canonical orientation in the DS (the Intersector uses the
            // edge data directly). OCCT returns the orientation separately —
            // callers read it from the DS edge shape.
            let _ = ori;
            Some(ce)
        } else {
            None
        }
    }

    /// The current edge's orientation (OCCT CurrentEdge returns it).
    pub fn current_edge_orientation(&self) -> Orientation {
        if self.current_wire < self.wires.len()
            && self.current_edge < self.wires[self.current_wire].len()
        {
            self.wires[self.current_wire][self.current_edge].1
        } else {
            Orientation::Forward
        }
    }
}
