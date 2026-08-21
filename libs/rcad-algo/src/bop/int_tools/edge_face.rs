// OCCT IntTools_EdgeFace ?edge-face intersection
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Surface3};
#[derive(Debug, Clone)]
pub struct EdgeFaceHit { pub point: DVec3, pub edge_param: f64 }

pub fn intersect_line_plane(line: &rcad_kernel::geom::Line3, t_range: [f64; 2], plane: &rcad_kernel::geom::Plane) -> Option<EdgeFaceHit> {
    let denom = line.direction.dot(plane.normal);
    if denom.abs() < rcad_kernel::CONFUSION { return None; }
    let t = (plane.origin - line.origin).dot(plane.normal) / denom;
    if t < t_range[0] - rcad_kernel::CONFUSION || t > t_range[1] + rcad_kernel::CONFUSION { return None; }
    Some(EdgeFaceHit { point: line.origin + line.direction * t, edge_param: t })
}

/// Same as [`intersect_line_plane`] with explicit edge-parameter margin
/// `param_tol` (minimum CONFUSION). Parallel/near-parallel denominator
/// threshold stays strict at CONFUSION.
pub fn intersect_line_plane_with_tol(
    line: &rcad_kernel::geom::Line3,
    t_range: [f64; 2],
    plane: &rcad_kernel::geom::Plane,
    param_tol: f64,
) -> Option<EdgeFaceHit> {
    let ptol = param_tol.max(rcad_kernel::CONFUSION);
    let denom = line.direction.dot(plane.normal);
    if denom.abs() < rcad_kernel::CONFUSION { return None; }
    let t = (plane.origin - line.origin).dot(plane.normal) / denom;
    if t < t_range[0] - ptol || t > t_range[1] + ptol { return None; }
    Some(EdgeFaceHit { point: line.origin + line.direction * t, edge_param: t })
}
pub fn plane_local_basis(plane: &rcad_kernel::geom::Plane) -> (DVec3, DVec3) {
    let n = plane.normal;
    let ref_dir = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    (n.cross(ref_dir).normalize(), n.cross(n.cross(ref_dir)).normalize())
}
pub fn point_in_planar_face(point: DVec3, plane: &rcad_kernel::geom::Plane, face_verts: &[DVec3]) -> bool {
    if face_verts.len() < 3 { return false; }
    let (u_axis, v_axis) = plane_local_basis(plane);
    let project = |p: DVec3| { let d = p - plane.origin; (d.dot(u_axis), d.dot(v_axis)) };
    let (px, py) = project(point);
    let poly: Vec<(f64, f64)> = face_verts.iter().map(|v| project(*v)).collect();
    let n = poly.len(); let mut inside = false; let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i]; let (xj, yj) = poly[j];
        if (yi > 0.0) != (yj > 0.0) {
            let xint = (xj - xi) * (0.0 - yi) / (yj - yi) + xi;
            if px < xint { inside = !inside; }
        }
        j = i;
    }
    inside
}

// ============================================================================
// OCCT IntTools_EdgeFace — edge/face intersection (IntTools_EdgeFace.cxx).
// ============================================================================
use crate::bop::ds::DS;
use crate::bop::int_tools::bean_face_intersector::{
    BeanFaceIntersector, BRepAdaptorCurve, BRepAdaptorSurface, ExtremaExtCS,
    GeomAbsCurveType, GeomAbsSurfaceType, IntCurveSurfaceHInter,
};
use crate::bop::int_tools::common_prt::{CommonPrt, CommonPrtType};
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::int_tools::face_make_curve::intermediate_point;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

/// OCCT IntTools_EdgeFace — determines the intersection of an edge with a
/// face on a given parameter range, producing common parts of type VERTEX
/// or EDGE (IntTools_EdgeFace.cxx L48-875).
pub struct EdgeFace<'a> {
    ds: &'a DS,
    n_e: usize,
    n_f: usize,
    my_edge: BRepAdaptorCurve,
    my_face: BRepAdaptorSurface,
    my_range: [f64; 2],
    my_fuzzy_value: f64,
    my_criteria: f64,
    my_quick_coincidence_check: bool,
    my_is_done: bool,
    my_error_status: i32,
    my_common_parts: Vec<CommonPrt>,
    my_min_distance: f64,
    edge_tol: f64,
    face_tol: f64,
}

/// OCCT IsEqDistance (L274-315) — for Cylinder/Cone/Torus the distance to
/// the surface equals the analytic radius when the point is on/near the axis.
fn is_eq_distance(
    a_p: DVec3,
    a_bas: &BRepAdaptorSurface,
    a_tol: f64,
    a_d: &mut f64,
) -> bool {
    match a_bas.get_type() {
        GeomAbsSurfaceType::Cylinder => {
            let a_cyl = a_bas.cylinder();
            let a_lin_axis = rcad_kernel::geom::Line3::new(a_cyl.origin, a_cyl.axis);
            let a_dc = a_lin_axis.distance(a_p);
            if a_dc < a_tol {
                *a_d = a_cyl.radius;
                return true;
            }
        }
        GeomAbsSurfaceType::Cone => {
            let a_cone = a_bas.cone();
            let a_lin_axis = rcad_kernel::geom::Line3::new(a_cone.apex, a_cone.axis);
            let a_dc = a_lin_axis.distance(a_p);
            if a_dc < a_tol {
                let an_apex = a_cone.apex;
                let a_semi_angle = a_cone.half_angle_rad;
                let a_ds = a_p.distance(an_apex);
                *a_d = a_ds * a_semi_angle.tan();
                return true;
            }
        }
        GeomAbsSurfaceType::Torus => {
            let a_torus = a_bas.torus();
            let a_ploc = a_torus.center;
            let a_major_radius = a_torus.major_radius;
            let a_dc = (a_ploc.distance(a_p) - a_major_radius).abs();
            if a_dc < a_tol {
                *a_d = a_torus.minor_radius;
                return true;
            }
        }
        _ => {}
    }
    false
}

/// OCCT IsCoplanar (L780-806) — Circle curve coplanar with Plane surface.
fn is_coplanar(a_curve: &BRepAdaptorCurve, a_surface: &BRepAdaptorSurface) -> bool {
    if a_curve.get_type() == GeomAbsCurveType::Circle
        && a_surface.get_type() == GeomAbsSurfaceType::Plane
    {
        let a_circ = a_curve.circle();
        let a_dir_ax1 = a_circ.normal;
        let a_pln = a_surface.plane();
        let a_dir_pln = a_pln.normal;
        return crate::bop::int_tools::face_face::is_dirs_coinside(a_dir_ax1, a_dir_pln);
    }
    false
}

/// OCCT IsRadius (L808-838) — Circle coplanar with Plane and its center at
/// the plane distance equal to the radius (within criteria).
fn is_radius(
    a_curve: &BRepAdaptorCurve,
    a_surface: &BRepAdaptorSurface,
    a_criteria: f64,
) -> bool {
    if a_curve.get_type() == GeomAbsCurveType::Circle
        && a_surface.get_type() == GeomAbsSurfaceType::Plane
    {
        let a_circ = a_curve.circle();
        let a_center = a_circ.center;
        let a_r = a_circ.radius;
        let a_pln = a_surface.plane();
        let a_d = (a_center - a_pln.origin).dot(a_pln.normal).abs();
        if (a_d - a_r).abs() < a_criteria {
            return true;
        }
    }
    false
}

impl<'a> EdgeFace<'a> {
    /// rcad constructor — from DS edge/face indices. Returns None when the
    /// edge has no geometric curve (OCCT CheckData error status 3) — the
    /// caller skips such pairs.
    pub fn new(ds: &'a DS, n_e: usize, n_f: usize) -> Option<Self> {
        let curve = ds.edge_curve(n_e)?.clone();
        let surf = ds.face_surface(n_f)?;
        let range = ds.edge_range(n_e);
        let edge_tol = ds.edge_tolerance(n_e);
        let face_tol = ds.face_tolerance(n_f);
        Some(EdgeFace {
            ds,
            n_e,
            n_f,
            my_edge: BRepAdaptorCurve::with_range(curve, range[0], range[1]),
            my_face: BRepAdaptorSurface::new(surf),
            my_range: range,
            my_fuzzy_value: CONFUSION,
            my_criteria: 0.0,
            my_quick_coincidence_check: false,
            my_is_done: false,
            my_error_status: 1,
            my_common_parts: Vec::new(),
            my_min_distance: f64::MAX,
            edge_tol,
            face_tol,
        })
    }

    pub fn set_range(&mut self, a_first: f64, a_last: f64) {
        self.my_range = [a_first, a_last];
    }
    /// OCCT IntTools_EdgeFace::SetFuzzyValue (IntTools_EdgeFace.hxx L85-88):
    /// myFuzzyValue = max(theFuzz, Precision::Confusion()).
    pub fn set_fuzzy_value(&mut self, the_fuzz: f64) {
        self.my_fuzzy_value = the_fuzz.max(CONFUSION);
    }
    pub fn use_quick_coincidence_check(&mut self, the_flag: bool) {
        self.my_quick_coincidence_check = the_flag;
    }
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }
    pub fn error_status(&self) -> i32 {
        self.my_error_status
    }
    pub fn common_parts(&self) -> &[CommonPrt] {
        &self.my_common_parts
    }
    pub fn edge_index(&self) -> usize {
        self.n_e
    }
    pub fn face_index(&self) -> usize {
        self.n_f
    }
    pub fn min_distance(&self) -> f64 {
        self.my_min_distance
    }
    pub fn criteria(&self) -> f64 {
        self.my_criteria
    }
    /// Point on the edge at parameter t (OCCT BRepAdaptor_Curve::Value).
    pub fn edge_value(&self, t: f64) -> DVec3 {
        self.my_edge.value(t)
    }

    /// OCCT CheckData (L133-141) — degenerated edge → status 2.
    fn check_data(&mut self) {
        let si = self.ds.shape_info(self.n_e);
        let degen = match &*si.shape.data {
            rcad_kernel::topods::TShape::Edge(e) => e.degenerated,
            _ => true,
        };
        if degen {
            self.my_error_status = 2;
        }
    }

    /// OCCT IsProjectable (L143-151) — the point at t is valid for the face
    /// within myCriteria.
    fn is_projectable(&self, a_t: f64) -> bool {
        let a_pc = self.my_edge.value(a_t);
        let mut a_ctx = IntToolsContext::new();
        a_ctx.is_valid_point_for_face(a_pc, self.n_f, self.ds, self.my_criteria)
    }

    /// OCCT DistanceFunction (L153-183) — signed distance point→face minus
    /// myCriteria (analytic for Cylinder/Cone/Torus).
    fn distance_function(&self, t: f64) -> f64 {
        let p = self.my_edge.value(t);
        let mut a_d = 0.0;
        if is_eq_distance(p, &self.my_face, 1e-7, &mut a_d) {
            return a_d - self.my_criteria;
        }
        let mut a_ctx = IntToolsContext::new();
        let proj = a_ctx.proj_ps(self.ds, self.n_f);
        proj.perform(p);
        if proj.nb_points() == 0 {
            return 99.0;
        }
        proj.lower_distance() - self.my_criteria
    }

    /// OCCT IsCoincident (L55-132) — sample the edge, project to the face and
    /// classify: >50% of the samples inside the face and within criteria.
    fn is_coincident(&self) -> bool {
        let mut a_ctx = IntToolsContext::new();
        let mut a_nb_seg = 23;
        if self.my_edge.get_type() == GeomAbsCurveType::Line
            && self.my_face.get_type() == GeomAbsSurfaceType::Plane
        {
            a_nb_seg = 2; // Check only three points for Line/Plane
        }
        let a_tresh = 0.5;
        let a_tresh_idx_f = ((a_nb_seg + 1) as f64 * 0.25) as usize;
        let a_tresh_idx_l = ((a_nb_seg + 1) as f64 * 0.75) as usize;
        let mut a_t1 = self.my_range[0];
        let mut a_t2 = self.my_range[1];
        let a_bnd_shift = 0.01 * (a_t2 - a_t1);
        a_t1 += a_bnd_shift;
        a_t2 -= a_bnd_shift;
        let d_t = (a_t2 - a_t1) / a_nb_seg as f64;
        let mut is_classified = false;
        let mut i_cnt = 0;
        for i in 0..=a_nb_seg {
            let a_t = a_t1 + i as f64 * d_t;
            let a_p = self.my_edge.value(a_t);
            let proj = a_ctx.proj_ps(self.ds, self.n_f);
            proj.perform(a_p);
            if proj.nb_points() == 0 {
                continue;
            }
            let a_d = proj.lower_distance();
            if a_d > self.my_criteria {
                if a_d > 100.0 * self.my_criteria {
                    return false;
                }
                continue;
            }
            i_cnt += 1;
            if ((0 < i) && (i < a_tresh_idx_f)) || ((a_tresh_idx_l < i) && (i < a_nb_seg)) {
                continue;
            }
            if is_classified && (i != a_nb_seg) {
                continue;
            }
            let (a_u, a_v) = proj.lower_distance_parameters();
            let a_p2d = DVec2::new(a_u, a_v);
            let a_class2d = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
                self.ds,
                self.n_f,
                self.ds.face_tolerance(self.n_f),
            );
            let a_state = a_class2d.perform(self.ds, a_p2d, true);
            if a_state == crate::topalgo::brep_top_adaptor::fclass2d::State::Out {
                return false;
            }
            if i != 0 {
                is_classified = true;
            }
        }
        let a_coeff = i_cnt as f64 / (a_nb_seg as f64 + 1.0);
        a_coeff > a_tresh
    }

    /// OCCT MakeType (L317-385) — classify a common part as EDGE or VERTEX.
    fn make_type(&self, a_common_prt: &mut CommonPrt) {
        if a_common_prt.all_null_flag() {
            a_common_prt.set_type(CommonPrtType::Edge);
            return;
        }
        let (af1, al1) = (a_common_prt.range1[0], a_common_prt.range1[1]);
        let a_pf = self.my_edge.value(af1);
        let a_pl = self.my_edge.value(al1);
        let df1 = a_pf.distance(a_pl);
        let a_res = self.my_edge.resolution(self.my_criteria);
        let is_whole_range = (af1 - self.my_range[0]).abs() < a_res
            && (al1 - self.my_range[1]).abs() < a_res;
        if (df1 > self.my_criteria * 2.0) && is_whole_range {
            a_common_prt.set_type(CommonPrtType::Edge);
            return;
        }
        if is_whole_range {
            let tm = (af1 + al1) * 0.5;
            if a_pf.distance(self.my_edge.value(tm)) > self.my_criteria * 2.0 {
                a_common_prt.set_type(CommonPrtType::Edge);
                return;
            }
        }
        let mut tm = (af1 + al1) * 0.5;
        if !self.check_touch(a_common_prt, &mut tm) {
            tm = (af1 + al1) * 0.5;
        }
        a_common_prt.set_type(CommonPrtType::Vertex);
        a_common_prt.set_vertex_parameter1(tm);
        a_common_prt.set_range1(af1, al1);
    }

    /// OCCT CheckTouch (L387-485) — decide whether the range touches the face
    /// at a single point (VERTEX) rather than over a range (EDGE).
    fn check_touch(&self, a_cp: &CommonPrt, a_tx: &mut f64) -> bool {
        let (a_tf, a_tl) = (a_cp.range1[0], a_cp.range1[1]);
        let a_cr = self.my_edge.resolution(self.my_criteria);
        if (a_tf - self.my_range[0]).abs() < a_cr && (a_tl - self.my_range[1]).abs() < a_cr {
            return false; // EDGE
        }
        let tol = PCONFUSION;
        let (u1f, u1l, v1f, v1l) = (
            self.my_face.first_u_parameter(),
            self.my_face.last_u_parameter(),
            self.my_face.first_v_parameter(),
            self.my_face.last_v_parameter(),
        );
        let mut a_extrema = ExtremaExtCS::new();
        a_extrema.initialize_with_bounds(
            self.my_face.surface(),
            u1f,
            u1l,
            v1f,
            v1l,
            tol,
            tol,
        );
        a_extrema.perform(&self.my_edge, a_tf, a_tl);
        let mut a_dist2 = 1e100;
        let mut a_tx_out = *a_tx;
        if a_extrema.is_done() {
            if !a_extrema.is_parallel() {
                let a_nb_ext = a_extrema.nb_ext();
                if a_nb_ext > 0 {
                    let mut i_lower = 1;
                    for i in 1..=a_nb_ext {
                        let a_sq = a_extrema.square_distance(i);
                        if a_sq < a_dist2 {
                            a_dist2 = a_sq;
                            i_lower = i;
                        }
                    }
                    let (a_p_on_c, _a_p_on_s) = a_extrema.points(i_lower);
                    a_tx_out = a_p_on_c.parameter();
                } else {
                    // OCCT: IntCurveSurface_HInter fallback
                    let mut an_exact = IntCurveSurfaceHInter::new();
                    an_exact.perform(&self.my_edge, &self.my_face);
                    if an_exact.is_done() {
                        let a_nb = an_exact.nb_points();
                        for i in 1..=a_nb {
                            let a_point = an_exact.point(i);
                            let w = a_point.w();
                            if w >= a_tf && w <= a_tl {
                                a_dist2 = 0.0;
                                a_tx_out = w;
                            }
                        }
                    }
                }
            } else {
                return false;
            }
        }
        let mut a_boundary_dist = self.distance_function(a_tf) + self.my_criteria;
        if a_boundary_dist * a_boundary_dist < a_dist2 {
            a_dist2 = a_boundary_dist * a_boundary_dist;
            a_tx_out = a_tf;
        }
        a_boundary_dist = self.distance_function(a_tl) + self.my_criteria;
        if a_boundary_dist * a_boundary_dist < a_dist2 {
            a_dist2 = a_boundary_dist * a_boundary_dist;
            a_tx_out = a_tl;
        }
        let a_parameter = (a_tf + a_tl) * 0.5;
        a_boundary_dist = self.distance_function(a_parameter) + self.my_criteria;
        if a_boundary_dist * a_boundary_dist < a_dist2 {
            a_dist2 = a_boundary_dist * a_boundary_dist;
            a_tx_out = a_parameter;
        }
        if a_dist2 > self.my_criteria * self.my_criteria {
            return false;
        }
        if (a_tx_out - a_tf).abs() < PCONFUSION {
            return true;
        }
        if (a_tx_out - a_tl).abs() < PCONFUSION {
            return true;
        }
        if a_tx_out > a_tf && a_tx_out < a_tl {
            return true;
        }
        false
    }

    /// OCCT CheckTouchVertex (L719-777) — refine a VERTEX touch parameter.
    fn check_touch_vertex(&self, a_cp: &CommonPrt, a_tx: &mut f64) -> bool {
        let (a_tf, a_tl) = (a_cp.range1[0], a_cp.range1[1]);
        let a_type = self.my_edge.get_type();
        let mut a_eps_t = 8e-5;
        if a_type == GeomAbsCurveType::Line {
            a_eps_t = 9e-5;
        }
        let a_tm = 0.5 * (a_tf + a_tl);
        let mut a_dist2 = self.distance_function(a_tm);
        a_dist2 *= a_dist2;
        let tol = PCONFUSION;
        let (u1f, u1l, v1f, v1l) = (
            self.my_face.first_u_parameter(),
            self.my_face.last_u_parameter(),
            self.my_face.first_v_parameter(),
            self.my_face.last_v_parameter(),
        );
        let mut a_extrema = ExtremaExtCS::new();
        a_extrema.initialize_with_bounds(
            self.my_face.surface(),
            u1f,
            u1l,
            v1f,
            v1l,
            tol,
            tol,
        );
        a_extrema.perform(&self.my_edge, a_tf, a_tl);
        if !a_extrema.is_done() {
            return false;
        }
        if a_extrema.is_parallel() {
            return false;
        }
        let a_nb_ext = a_extrema.nb_ext();
        if a_nb_ext == 0 {
            return false;
        }
        let mut i_lower = 1;
        let mut a_min_dist2 = 1e100;
        for i in 1..=a_nb_ext {
            let a_sq = a_extrema.square_distance(i);
            if a_sq < a_min_dist2 {
                a_min_dist2 = a_sq;
                i_lower = i;
            }
        }
        let a_dist2_new = a_extrema.square_distance(i_lower);
        if a_dist2_new > a_dist2 {
            *a_tx = a_tm;
            return true;
        }
        if a_dist2_new > self.my_criteria * self.my_criteria {
            return false;
        }
        let (a_p_on_c, _) = a_extrema.points(i_lower);
        *a_tx = a_p_on_c.parameter();
        if (*a_tx - a_tf).abs() < a_eps_t {
            return false;
        }
        if (*a_tx - a_tl).abs() < a_eps_t {
            return false;
        }
        if *a_tx > a_tf && *a_tx < a_tl {
            return true;
        }
        false
    }

    /// OCCT Perform (L488-660) — main entry.
    pub fn perform(&mut self) {
        let mut a_common_prt = CommonPrt::new();
        self.my_error_status = 0;
        self.check_data();
        if self.my_error_status != 0 {
            return;
        }
        self.my_is_done = false;
        self.my_common_parts.clear();
        let a_curve_type = self.my_edge.get_type();
        let a_surf_type = self.my_face.get_type();
        // Prepare myCriteria (L524-537)
        let a_fuzz = self.my_fuzzy_value / 2.0;
        let a_tol_f = self.face_tol + a_fuzz;
        let a_tol_e = self.edge_tol + a_fuzz;
        if a_curve_type == GeomAbsCurveType::BSplineCurve
            || a_curve_type == GeomAbsCurveType::BezierCurve
        {
            let diff1 = a_tol_e / a_tol_f;
            let diff2 = a_tol_f / a_tol_e;
            if diff1 > 100.0 || diff2 > 100.0 {
                self.my_criteria = a_tol_e.max(a_tol_f);
            } else {
                self.my_criteria = 1.5 * a_tol_e + a_tol_f;
            }
        } else {
            self.my_criteria = a_tol_e + a_tol_f;
        }
        // Quick coincidence check (L539-548)
        if self.my_quick_coincidence_check {
            if self.is_coincident() {
                a_common_prt.set_type(CommonPrtType::Edge);
                a_common_prt.set_range1(self.my_range[0], self.my_range[1]);
                self.my_common_parts.push(a_common_prt);
                self.my_is_done = true;
                return;
            }
        }
        // BeanFaceIntersector (L550-566)
        let mut an_intersector =
            BeanFaceIntersector::with_adaptors(self.my_edge.clone(), self.my_face.clone(), a_tol_e, a_tol_f);
        an_intersector.set_bean_parameters(self.my_range[0], self.my_range[1]);
        an_intersector.perform();
        if an_intersector.minimal_square_distance() < f64::MAX {
            self.my_min_distance = an_intersector.minimal_square_distance().sqrt();
        }
        if !an_intersector.is_done() {
            return;
        }
        for a_range in an_intersector.result() {
            if self.is_projectable(intermediate_point(a_range.first(), a_range.last())) {
                a_common_prt.set_range1(a_range.first(), a_range.last());
                self.my_common_parts.push(a_common_prt.clone());
            }
        }
        let a_nb = self.my_common_parts.len();
        for i in 0..a_nb {
            let mut a_cp = self.my_common_parts[i].clone();
            let (a_tx1, a_tx2) = (a_cp.range1[0], a_cp.range1[1]);
            let a_px1 = self.my_edge.value(a_tx1);
            let a_px2 = self.my_edge.value(a_tx2);
            a_cp.set_bounding_points(a_px1, a_px2);
            self.make_type(&mut a_cp);
            self.my_common_parts[i] = a_cp;
        }
        // Line/Cylinder common parts treatment (L590-622)
        if a_curve_type == GeomAbsCurveType::Line && a_surf_type == GeomAbsSurfaceType::Cylinder {
            for i in 0..a_nb {
                let mut a_cp = self.my_common_parts[i].clone();
                let a_type = a_cp.get_type();
                let mut a_tx = 0.0;
                if a_type == CommonPrtType::Edge {
                    if self.check_touch(&a_cp, &mut a_tx) {
                        a_cp.set_type(CommonPrtType::Vertex);
                        a_cp.set_vertex_parameter1(a_tx);
                    }
                } else if a_type == CommonPrtType::Vertex {
                    if self.check_touch_vertex(&a_cp, &mut a_tx) {
                        a_cp.set_vertex_parameter1(a_tx);
                    }
                }
                self.my_common_parts[i] = a_cp;
            }
        }
        // Circle/Plane common parts treatment (L624-656)
        if a_curve_type == GeomAbsCurveType::Circle && a_surf_type == GeomAbsSurfaceType::Plane {
            let b_is_coplanar = is_coplanar(&self.my_edge, &self.my_face);
            let b_is_radius = is_radius(&self.my_edge, &self.my_face, self.my_criteria);
            if !b_is_coplanar && !b_is_radius {
                for i in 0..a_nb {
                    let mut a_cp = self.my_common_parts[i].clone();
                    let a_type = a_cp.get_type();
                    let mut a_tx = 0.0;
                    if a_type == CommonPrtType::Edge {
                        if self.check_touch(&a_cp, &mut a_tx) {
                            a_cp.set_type(CommonPrtType::Vertex);
                            a_cp.set_vertex_parameter1(a_tx);
                        }
                    } else if a_type == CommonPrtType::Vertex {
                        if self.check_touch_vertex(&a_cp, &mut a_tx) {
                            a_cp.set_vertex_parameter1(a_tx);
                        }
                    }
                    self.my_common_parts[i] = a_cp;
                }
            }
        }
        self.my_is_done = true;
    }
}

// ============================================================================
// Tests
// ============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bop::ds::DS;
    use rcad_kernel::geom::{Curve3, Surface3};
    use rcad_kernel::topods::{Orientation, ShapeType, TShape};
    use rcad_kernel::topo_shape::Shape;
    use rcad_modeling::prim::primapi::make_box_brep;

    /// Build a DS with a unit box at the origin and return (ds, top-face idx,
    /// an edge of the top face at z=1).
    fn box_ds() -> (DS, usize, usize) {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
        let mut root_ts = None;
        let mut root_idx = 0;
        for (i, ts) in brep.tshapes.iter().enumerate().rev() {
            if matches!(&**ts, TShape::Solid(_)) {
                root_ts = Some(ts.clone());
                root_idx = i;
                break;
            }
        }
        let root = Shape::from_parts(root_ts.expect("solid"), root_idx, 0, Orientation::Forward);
        let mut ds = DS::new();
        ds.set_arguments(vec![root]);
        ds.init(1e-7);
        // Top face: Plane with normal (0,0,1) at z=1.
        let n_f = (0..ds.nb_shapes())
            .find(|&i| {
                ds.shape_info(i).shape_type == ShapeType::Face
                    && ds.face_surface(i).map_or(false, |s| {
                        matches!(s, Surface3::Plane(p) if p.normal.z > 0.9 && p.origin.z > 0.9)
                    })
            })
            .expect("top face");
        // An edge of the top face: horizontal line along X at z=1.
        let n_e = (0..ds.nb_shapes())
            .find(|&i| {
                ds.shape_info(i).shape_type == ShapeType::Edge
                    && ds.edge_curve(i).map_or(false, |c| match c {
                        Curve3::Line(l) => {
                            l.origin.z > 0.9 && l.direction.x.abs() > 0.9 && l.direction.z.abs() < 1e-9
                        }
                        _ => false,
                    })
            })
            .expect("top-face edge");
        (ds, n_f, n_e)
    }

    /// A top-face edge coincides with the top face: the quick coincidence
    /// check must produce a single TopAbs_EDGE common part covering the
    /// whole range (OCCT IntTools_EdgeFace::Perform L539-548).
    #[test]
    fn edge_face_quick_coincidence_returns_edge_part() {
        let (ds, n_f, n_e) = box_ds();
        let mut ef = EdgeFace::new(&ds, n_e, n_f).expect("edge-face");
        let range = ds.edge_range(n_e);
        ef.set_range(range[0], range[1]);
        ef.use_quick_coincidence_check(true);
        ef.perform();
        assert!(ef.is_done(), "perform must complete");
        assert_eq!(ef.error_status(), 0, "no error expected");
        let parts = ef.common_parts();
        assert_eq!(parts.len(), 1, "coincident edge-face must yield one common part");
        assert_eq!(parts[0].get_type(), CommonPrtType::Edge);
        assert!((parts[0].range1[0] - range[0]).abs() < 1e-9);
        assert!((parts[0].range1[1] - range[1]).abs() < 1e-9);
    }

    /// A top-face edge vs the bottom face (z=0) is parallel and 1 unit away:
    /// no intersection, no common parts, no error.
    #[test]
    fn edge_face_disjoint_parallel_yields_no_parts() {
        let (ds, _, n_e) = box_ds();
        let n_f_bottom = (0..ds.nb_shapes())
            .find(|&i| {
                ds.shape_info(i).shape_type == ShapeType::Face
                    && ds.face_surface(i).map_or(false, |s| {
                        matches!(s, Surface3::Plane(p) if p.normal.z > 0.9 && p.origin.z < 1e-9)
                    })
            })
            .expect("bottom face");
        let mut ef = EdgeFace::new(&ds, n_e, n_f_bottom).expect("edge-face");
        let range = ds.edge_range(n_e);
        ef.set_range(range[0], range[1]);
        ef.perform();
        assert_eq!(ef.error_status(), 0);
        assert!(ef.common_parts().is_empty(), "disjoint edge/face must have no parts");
    }
}
