// OCCT IntTools_Context — intersection context for VE/EF computations.
//
// OCCT ref: IntTools_Context.hxx / IntTools_Context.cxx
//
// Provides cached projection onto curves/surfaces, point-in-face classification,
// and surface adaptors. Each tool is lazily created and cached per face index.

use std::collections::HashMap;
use std::sync::Arc;
use crate::bop::ds::DS;
use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;
use crate::topalgo::brep_class3d::solid_explorer::SolidExplorer;
use crate::topalgo::brep_top_adaptor::fclass2d::{FClass2d, State};
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval, Curve2dEval};
use rcad_kernel::topods::{ShapeType, TShape};
use glam::{DVec2, DVec3};

// ====================================================================
// GeomAPI_ProjectPointOnSurf — OCCT GeomAPI_ProjectPointOnSurf
// ====================================================================
/// OCCT GeomAPI_ProjectPointOnSurf — projects a 3D point onto a surface.
pub struct ProjectOnSurface {
    surf: Surface3,
    uv_bounds: [f64; 4],
    tolerance: f64,
    // Last projection result
    last_point: Option<DVec3>,
    last_uv: Option<DVec2>,
    last_distance: f64,
}

impl ProjectOnSurface {
    /// OCCT: Init(aS, Umin, Usup, Vmin, Vsup, Tol)
    pub fn init(&mut self, surf: Surface3, uv_bounds: [f64; 4], tolerance: f64) {
        self.surf = surf;
        self.uv_bounds = uv_bounds;
        self.tolerance = tolerance;
        self.last_point = None;
        self.last_uv = None;
        self.last_distance = f64::MAX;
    }

    /// OCCT: Perform(aP) — find closest point on surface.
    pub fn perform(&mut self, point: DVec3) {
        let (uv, proj) = crate::bop::closest_point_on_surface(&self.surf, point);
        self.last_uv = Some(uv);
        self.last_point = Some(proj);
        self.last_distance = (proj - point).length();
    }

    /// OCCT: NbPoints() — number of solutions found.
    pub fn nb_points(&self) -> usize {
        if self.last_point.is_some() { 1 } else { 0 }
    }

    /// OCCT: LowerDistance() — minimal distance from point to surface.
    pub fn lower_distance(&self) -> f64 {
        self.last_distance
    }

    /// OCCT: LowerDistanceParameters(U, V) — UV of the closest point.
    pub fn lower_distance_parameters(&self) -> (f64, f64) {
        self.last_uv.map(|uv| (uv.x, uv.y)).unwrap_or((0.0, 0.0))
    }

    /// OCCT: NearestPoint() — 3D coordinates of the closest point on surface.
    pub fn nearest_point(&self) -> DVec3 {
        self.last_point.unwrap_or(DVec3::ZERO)
    }
}

// ====================================================================
// BRepAdaptor_Surface — OCCT BRepAdaptor_Surface
// ====================================================================
/// OCCT BRepAdaptor_Surface — surface adaptor with type and derivative queries.
pub struct SurfaceAdaptor {
    surf: Surface3,
}

/// Surface type classification, analogous to OCCT GeomAbs_SurfaceType.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BSpline,
    Bezier,
    Other,
}

impl SurfaceAdaptor {
    /// Create adaptor from a surface.
    pub fn new(surf: Surface3) -> Self {
        SurfaceAdaptor { surf }
    }

    /// OCCT: GetType() — returns the surface type.
    pub fn get_type(&self) -> SurfaceType {
        match self.surf {
            Surface3::Plane(_) => SurfaceType::Plane,
            Surface3::Cylinder(_) => SurfaceType::Cylinder,
            Surface3::Cone(_) => SurfaceType::Cone,
            Surface3::Sphere(_) => SurfaceType::Sphere,
            Surface3::Torus(_) => SurfaceType::Torus,
            Surface3::BSpline(_) => SurfaceType::BSpline,
            _ => SurfaceType::Other,
        }
    }

    /// OCCT: D1(U, V) — returns (point, dU, dV).
    pub fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        self.surf.derivatives(u, v)
    }
}

// ====================================================================
// IntTools_Context
// ====================================================================
/// OCCT IntTools_Context — cached geometric tools for intersection.
///
/// Lazily creates and caches:
/// - ProjectOnSurface per face (ProjPS)
/// - SurfaceAdaptor per face (SurfaceAdaptor)
/// - Face classifier (FClass2d / IsPointInFace)
pub struct IntToolsContext {
    // OCCT: myProjPSMap — maps face → GeomAPI_ProjectPointOnSurf*
    proj_ps_cache: HashMap<usize, ProjectOnSurface>,
    // OCCT: mySurfAdaptorMap — maps face → BRepAdaptor_Surface*
    surf_adapt_cache: HashMap<usize, SurfaceAdaptor>,
    // rcad: cached UV bounds per face (OCCT: UVBounds computed on demand)
    uv_bounds_cache: HashMap<usize, [f64; 4]>,
}

impl IntToolsContext {
    /// OCCT: IntTools_Context() — default constructor.
    pub fn new() -> Self {
        IntToolsContext {
            proj_ps_cache: HashMap::new(),
            surf_adapt_cache: HashMap::new(),
            uv_bounds_cache: HashMap::new(),
        }
    }

    /// OCCT: Clear() — clear all cached data.
    pub fn clear(&mut self) {
        self.proj_ps_cache.clear();
        self.surf_adapt_cache.clear();
        self.uv_bounds_cache.clear();
    }

    // ====================================================================
    // ComputeVE — OCCT IntTools_Context::ComputeVE (IntTools_Context.cxx L499-545)
    // ====================================================================

    // OCCT IntTools_Context::ComputeVE (IntTools_Context.cxx L499-542).
    pub fn compute_ve(
        &mut self,
        n_vx: usize, n_e: usize,
        ds: &crate::bop::ds::DS,
        the_fuzz: f64,
    ) -> (i32, f64, f64) {
        // OCCT L505-508: degenerated edge check
        if ds.shapes[n_e].shape.as_edge().map_or(true, |ed| ed.degenerated) {
            return (-1, 0.0, 0.0);
        }
        // OCCT L509-512: non-geometric edge check — BRep_Tool::IsGeometric
        if ds.shapes[n_e].shape.as_edge().and_then(|ed| ed.curve.as_ref()).is_none() {
            return (-2, 0.0, 0.0);
        }
        // OCCT L517: vertex point — BRep_Tool::Pnt(theV)
        let a_p = ds.vertex_point_by_idx(n_vx);
        // OCCT L519-521: GeomAPI_ProjectPointOnCurve& aProjector = ProjPC(theE);
        let curve = match ds.edge_curve(n_e) { Some(c) => c.clone(), None => return (-3, 0.0, 0.0) };
        // OCCT ProjPC (IntTools_Context.cxx L276-281): BRep_Tool::Curve(theE, f, l)
        // then aProjector.Init(aC3D, f, l) — the projection is limited to the edge's
        // parameter range [f, l], not the full curve domain.
        let a_range = ds.shapes[n_e].shape.as_edge()
            .map(|ed| ed.range).unwrap_or([0.0, 0.0]);
        let (a_f, a_l) = (a_range[0], a_range[1]);
        let a_proj = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
            &curve, a_p, a_f, a_l, 128);
        let a_t = a_proj.param;
        let a_proj = a_proj.point;
        // OCCT L522-526: if (!aProjector.NbPoints()) return -3;
        let a_dist = (a_proj - a_p).length();
        // OCCT L530-532: tolerance sum
        let a_tol_v = ds.vertex_tolerance_by_idx(n_vx);
        let a_tol_e = ds.edge_tolerance(n_e);
        let a_tol_sum = a_tol_v + a_tol_e + the_fuzz.max(rcad_kernel::CONFUSION);
        // OCCT L534: theTol = aDist + aTolE;
        let a_tol = a_dist + a_tol_e;
        // OCCT L537-538: if (aDist > aTolSum) return -4;
        if a_dist > a_tol_sum {
            return (-4, a_t, a_tol);
        }
        // OCCT L542: return 0;
        (0, a_t, a_tol)
    }

    // ====================================================================
    // ComputeEF — OCCT IntTools_EdgeFace + IntTools_BeanFaceIntersector
    // ====================================================================

    pub fn compute_ef(
        &mut self, n_e: usize, n_f: usize, a_t1: f64, a_t2: f64,
        ds: &DS, the_fuzz: f64,
    ) -> (i32, Vec<(f64, f64, bool)>, f64) {
        // OCCT IntTools_EdgeFace::Perform (EdgeFace.cxx L529-549): myCriteria
        let a_fuzz = the_fuzz / 2.0;
        let a_tol_e = ds.edge_tolerance(n_e) + a_fuzz;
        let a_tol_f = ds.face_tolerance(n_f) + a_fuzz;
        let a_criteria = a_tol_e + a_tol_f;
        let curve = match ds.edge_curve(n_e) { Some(c) => c.clone(), None => return (-1, vec![], f64::MAX) };
        let surf = match ds.face_surface(n_f) { Some(s) => s, None => return (-1, vec![], f64::MAX) };
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            let ct = match &curve {
                Curve3::Line(_) => "Line", Curve3::Circle(_) => "Circle", Curve3::Ellipse(_) => "Ellipse",
                Curve3::BSpline(_) => "BSpline", Curve3::Bezier(_) => "Bezier", _ => "Other",
            };
            let st = match &surf {
                Surface3::Plane(_) => "Plane", Surface3::Cylinder(_) => "Cylinder", Surface3::Cone(_) => "Cone",
                Surface3::Sphere(_) => "Sphere", Surface3::Torus(_) => "Torus", _ => "Other",
            };
            eprintln!("[EF-DBG] compute_ef e={}({}) f={}({}) curve={} surf={} edgerange=[{:.4},{:.4}]", n_e, ds.rank(n_e), n_f, ds.rank(n_f), ct, st,
                ds.shapes[n_e].shape.as_edge().map(|ed| ed.range[0]).unwrap_or(0.0),
                ds.shapes[n_e].shape.as_edge().map(|ed| ed.range[1]).unwrap_or(0.0));
        }

        // OCCT IntTools_EdgeFace L551: myS = myContext->SurfaceAdaptor(myFace)
        // OCCT L565-570: IntTools_BeanFaceIntersector anIntersector(myC, myS, aTolE, aTolF);
        //   SetBeanParameters(range), SetContext, Perform.
        // myC uses the edge's real parameter range (BRepAdaptor_Curve from the
        // edge, not the bare curve's default domain).
        let edge_range = ds.shapes[n_e].shape.as_edge().map(|ed| ed.range).unwrap_or([a_t1, a_t2]);
        let mut bfi = crate::bop::int_tools::bean_face_intersector::BeanFaceIntersector::with_adaptors(
            crate::bop::int_tools::bean_face_intersector::BRepAdaptorCurve::with_range(
                curve.clone(), edge_range[0], edge_range[1]),
            crate::bop::int_tools::bean_face_intersector::BRepAdaptorSurface::new(surf.clone()),
            a_tol_e,
            a_tol_f,
        );
        let uv_bounds = ds.face_uv_boundary(n_f);
        bfi.set_surface_parameters(uv_bounds[0], uv_bounds[1], uv_bounds[2], uv_bounds[3]);
        bfi.set_bean_parameters(a_t1, a_t2);
        bfi.perform();

        // OCCT L572-575: myMinDistance = sqrt(anIntersector.MinimalSquareDistance())
        let mut min_sq_dist = bfi.minimal_square_distance();
        let min_dist = if min_sq_dist < f64::MAX { min_sq_dist.sqrt() } else { f64::MAX };

        // OCCT L577-580: if (!anIntersector.IsDone()) return;
        if !bfi.is_done() {
            return (0, vec![], min_dist);
        }

        // OCCT L582-590: for each result range — IsProjectable at the intermediate
        // point (IntTools_Tools::IntermediatePoint = midpoint) → keep the range
        // as a common part. OCCT IntTools_Context::IsValidPointForFace (L648-674)
        // projects the point onto the surface and classifies the 2D point.
        let mut common_parts: Vec<(f64, f64, bool)> = Vec::new();
        for range in bfi.result() {
            let (rf, rl) = (range.first(), range.last());
            let t_mid = (rf + rl) * 0.5;
            let p_mid = curve.point_at(t_mid);
            let (uv, _) = crate::bop::closest_point_on_surface(&surf, p_mid);
            // OCCT IntTools_EdgeFace L595-608: MakeType first, then IsProjectable.
            let (_, _, is_edge) = make_type_ef(&curve, rf, rl, a_t1, a_t2, a_criteria);
            // OCCT IsValidPointForFace uses IsPointInOnFace (ON included), but the
            // rcad polygon classifier must exclude ON-boundary points (like OCCT's
            // FClass2d does for the box×box coplanar edges) — keep IsPointInFace.
            let in_face = self.is_point_in_face(ds, n_f, uv);
            if std::env::var("RCAD_EE_DEBUG").is_ok() {
                eprintln!("[EF-DBG]   range=[{:.5},{:.5}] mid={:.5} uv=({:.4},{:.4}) inFace={} ty={}", rf, rl, t_mid, uv.x, uv.y, in_face, if is_edge { "E" } else { "V" });
            }
            if !in_face {
                continue;
            }
            common_parts.push((rf, rl, is_edge));
        }
        // Merge adjacent ranges (OCCT IntTools_BeanFaceIntersector::Perform L352-378
        // merges consecutive ranges whose gap is within PConfusion).
        if !common_parts.is_empty() {
            common_parts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let mut merged: Vec<(f64, f64, bool)> = Vec::new();
            let curve_res = curve_res_ef(&curve, a_criteria);
            for (s, e, ty) in common_parts {
                if let Some(last) = merged.last_mut() {
                    if s <= last.1 + curve_res {
                        last.1 = last.1.max(e);
                        last.2 = last.2 || ty;
                        continue;
                    }
                }
                merged.push((s, e, ty));
            }
            common_parts = merged;
        }
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            eprintln!("[EF-DBG] compute_ef e={} f={} r=[{:.4},{:.4}] kept={:?}", n_e, n_f, a_t1, a_t2,
                common_parts.iter().map(|&(s, e, ty)| format!("[{:.5},{:.5}]{}", s, e, if ty { "E" } else { "V" })).collect::<Vec<_>>());
        }
        (0, common_parts, min_dist)
    }

    /// OCCT IntTools_BeanFaceIntersector::ComputeLinePlane (L820-906).
    fn compute_line_plane_ef(
        &self, curve: &Curve3, surf: &Surface3, n_f: usize, ds: &DS,
        a_t1: f64, a_t2: f64, a_criteria: f64, bean_tol: f64, face_tol: f64,
    ) -> (Vec<(f64, f64, bool)>, f64) {
        let t_ang = 1e-9;
        let (Curve3::Line(l), Surface3::Plane(p)) = (curve, surf) else { return (vec![], f64::MAX) };
        // Plane coefficients A x + B y + C z + D = 0
        let (a, b, c) = (p.normal.x, p.normal.y, p.normal.z);
        let d = -(p.normal.dot(p.origin));
        let dir = l.direction.normalize_or_zero();
        let (al, bl, cl) = (dir.x, dir.y, dir.z);
        let direc = a * al + b * bl + c * cl;
        let dis = a * l.origin.x + b * l.origin.y + c * l.origin.z + d;
        let mut parallel = false;
        let mut inplane = false;
        if direc.abs() < t_ang {
            parallel = true;
            inplane = dis.abs() < a_criteria;
        } else {
            let p1 = l.point_at(a_t1);
            let p2 = l.point_at(a_t2);
            let mut d1 = a * p1.x + b * p1.y + c * p1.z + d;
            if d1 < 0.0 { d1 = -d1; }
            let mut d2 = a * p2.x + b * p2.y + c * p2.z + d;
            if d2 < 0.0 { d2 = -d2; }
            if d1 <= a_criteria && d2 <= a_criteria {
                inplane = true;
            }
        }
        if inplane {
            // OCCT IntTools_EdgeFace L586: IsProjectable(midpoint) — keep the
            // whole-range common part only if the edge's midpoint lies inside
            // the face (an edge coplanar with the plane but outside the face is
            // not an intersection).
            let t_mid = (a_t1 + a_t2) * 0.5;
            let p_mid = l.point_at(t_mid);
            let rel = p_mid - p.origin;
            let um = rel.dot(p.u_dir);
            let vm = rel.dot(p.v_dir);
            if self.is_point_in_face(ds, n_f, glam::DVec2::new(um, vm)) {
                return (vec![(a_t1, a_t2, true)], 0.0);
            }
            return (vec![], f64::MAX);
        }
        if parallel {
            return (vec![], f64::MAX);
        }
        let t = -dis / direc;
        if t < a_t1 || t > a_t2 {
            return (vec![], f64::MAX);
        }
        let pint = l.origin + dir * t;
        // ElSLib::Parameters(P, pint, u, v)
        let rel = pint - p.origin;
        let u = rel.dot(p.u_dir);
        let v = rel.dot(p.v_dir);
        // OCCT: myUMinParameter <= u <= myUMaxParameter etc. — the point must
        // project INSIDE the face (not just the infinite plane).
        if !self.is_point_in_face(ds, n_f, glam::DVec2::new(u, v)) {
            return (vec![], f64::MAX);
        }
        // ComputeIntRange for the correct range on the edge
        let an_angle = (std::f64::consts::FRAC_PI_2 - l.direction.normalize().angle_between(p.normal.normalize()).abs()).abs();
        let a_dt = compute_int_range_ef(bean_tol, face_tol, an_angle);
        let t1 = a_t1.max(t - a_dt);
        let t2 = a_t2.min(t + a_dt);
        (vec![(t1, t2, false)], 0.0)
    }

    /// OCCT IntTools_BeanFaceIntersector::ComputeRangeFromStartPoint (L1150-1340).
    /// Expands the intersection range from the start point `t0` in the
    /// increasing (ToIncreaseParameter=true) or decreasing direction.
    fn compute_range_from_start_point(
        &mut self,
        to_increase: bool,
        t0: f64, u0: f64, v0: f64,
        ranges: &mut Vec<(f64, f64, i32)>,
        curve: &Curve3, surf: &Surface3,
        a_criteria: f64, curve_res: f64,
    ) {
        let a_first = ranges[0].0;
        let a_last = ranges[ranges.len() - 1].1;
        let Some(mut valid_index) = ranges.iter().position(|r| r.0 <= t0 && t0 <= r.1) else { return; };
        if ranges[valid_index].2 > 0 { return; }
        let mut a_min_delta = curve_res * 0.5;
        let mut a_delta_restrictor = 0.1 * (a_last - a_first);
        if a_min_delta > a_delta_restrictor { a_min_delta = a_delta_restrictor * 0.5; }
        let ten_of_min_delta = a_min_delta * 10.0;
        let mut a_delta = curve_res;
        let mut a_cur_par = if to_increase { t0 + a_delta } else { t0 - a_delta };
        let mut a_prev_par = t0;
        let mut current_range = ranges[valid_index];
        let mut boundary_condition =
            if to_increase { a_cur_par > current_range.1 } else { a_cur_par < current_range.0 };
        if boundary_condition {
            a_cur_par = if to_increase { current_range.1 } else { current_range.0 };
            boundary_condition = false;
        }
        let mut loopcounter = 0;
        let mut u = u0;
        let mut v = v0;
        let mut another_solution_found = false;
        let mut isboundaryindex = false;
        let mut isvalidindex = true;
        while a_delta >= a_min_delta && loopcounter <= 10 {
            let mut pointfound = false;
            let a_point = curve.point_at(a_cur_par);
            let (uv, q) = crate::bop::closest_point_on_surface(surf, a_point);
            let dist = (q - a_point).length();
            if dist < a_criteria {
                u = uv.x;
                v = uv.y;
                pointfound = true;
            } else {
                // Extrema_GenLocateExtPS fallback — distance to the surface domain boundary
                pointfound = self.dist_to_surface_domain(surf, a_point) < a_criteria;
            }
            if pointfound {
                a_prev_par = a_cur_par;
                another_solution_found = true;
                if boundary_condition && (isboundaryindex || !isvalidindex) { break; }
            } else {
                a_delta_restrictor = a_delta;
            }
            a_delta = if pointfound { a_delta * 2.0 } else { a_delta * 0.5 };
            a_delta = a_delta.min(a_delta_restrictor);
            a_cur_par = if to_increase { a_prev_par + a_delta } else { a_prev_par - a_delta };
            if a_cur_par == a_prev_par { break; }
            boundary_condition =
                if to_increase { a_cur_par > current_range.1 } else { a_cur_par < current_range.0 };
            isboundaryindex = false;
            isvalidindex = true;
            if boundary_condition {
                isboundaryindex = (!to_increase && valid_index == 0)
                    || (to_increase && valid_index == ranges.len() - 1);
                if !isboundaryindex {
                    if pointfound {
                        let adj_flag = if to_increase { ranges[valid_index + 1].2 } else { ranges[valid_index - 1].2 };
                        if adj_flag == 0 {
                            valid_index = if to_increase { valid_index + 1 } else { valid_index - 1 };
                            current_range = ranges[valid_index];
                            if (to_increase && a_cur_par > current_range.1)
                                || (!to_increase && a_cur_par < current_range.0) {
                                a_cur_par = (current_range.0 + current_range.1) * 0.5;
                                a_delta *= 0.5;
                            }
                        } else {
                            isvalidindex = false;
                            a_cur_par = if to_increase { current_range.1 } else { current_range.0 };
                        }
                    }
                } else {
                    a_cur_par = if to_increase { current_range.1 } else { current_range.0 };
                }
                if a_delta < ten_of_min_delta { loopcounter += 1; } else { loopcounter = 0; }
            }
        }
        if another_solution_found {
            let (ns, ne) = if to_increase { (t0, a_prev_par) } else { (a_prev_par, t0) };
            if ne - ns > 1e-9 {
                insert_range_ef(ranges, ns, ne, 2);
            }
        }
        let _ = (u, v);
    }

    /// OCCT IntTools_BeanFaceIntersector::Distance fallback — distance from a
    /// point to the surface domain boundary corners/isoparameters.
    fn dist_to_surface_domain(&self, surf: &Surface3, p: DVec3) -> f64 {
        let [u0, u1, v0, v1] = surf.default_domain();
        let mut best = f64::MAX;
        let corners = [
            (u0, v0), (u1, v0), (u0, v1), (u1, v1),
        ];
        for (cu, cv) in corners {
            if !cu.is_finite() || !cv.is_finite() { continue; }
            let q = surf.point_at(cu, cv);
            let d = (q - p).length();
            if d < best { best = d; }
        }
        best
    }

    // ====================================================================
    // ProjPS — OCCT GeomAPI_ProjectPointOnSurf per face
    // ====================================================================

    /// OCCT IntTools_Context::ComputeVF (IntTools_Context.cxx L546-591).
    /// Projects the vertex onto the face surface; returns (flag, U, V, TolNew)
    /// where flag 0 = intersection, -1 = not projectable, -2 = distance too
    /// large, -3 = projection outside the face.
    pub fn compute_vf(&mut self, n_v: usize, n_f: usize, ds: &DS, the_fuzz: f64) -> (i32, f64, f64, f64) {
        // aP = BRep_Tool::Pnt(theVertex)
        let a_p = ds.vertex_point_by_idx(n_v);
        // 1. GeomAPI_ProjectPointOnSurf& aProjector = ProjPS(theFace); aProjector.Perform(aP)
        let (a_dist, a_u, a_v) = {
            let proj = self.proj_ps(ds, n_f);
            proj.perform(a_p);
            if proj.nb_points() == 0 {
                return (-1, 0.0, 0.0, 0.0);
            }
            let d = proj.lower_distance();
            let (u, v) = proj.lower_distance_parameters();
            (d, u, v)
        };
        // 2. aTolV/aTolF/aTolSum; theTol = aDist + aTolF
        let a_tol_v = ds.vertex_tolerance_by_idx(n_v);
        let a_tol_f = ds.face_tolerance(n_f);
        let a_tol_sum = a_tol_v + a_tol_f + the_fuzz.max(rcad_kernel::CONFUSION);
        let the_tol = a_dist + a_tol_f;
        // if (aDist > aTolSum) return -2;
        if a_dist > a_tol_sum {
            return (-2, a_u, a_v, the_tol);
        }
        // 3. IsPointInFace(theFace, aP2d)
        let pri = self.is_point_in_face(ds, n_f, glam::DVec2::new(a_u, a_v));
        if !pri {
            return (-3, a_u, a_v, the_tol);
        }
        (0, a_u, a_v, the_tol)
    }

    /// OCCT IntTools_Context::ProjPS (IntTools_Context.cxx L247-265).
    /// Returns a cached point-on-surface projector for the face `fi`.
    pub fn proj_ps(&mut self, ds: &DS, fi: usize) -> &mut ProjectOnSurface {
        if !self.proj_ps_cache.contains_key(&fi) {
            // OCCT L252-253: UVBounds + BRep_Tool::Surface
            let surf = ds.face_surface(fi)
                .expect("ProjPS: face has no surface");
            let uv_bounds = ds.face_uv_boundary(fi);
            // OCCT L257-260: new GeomAPI_ProjectPointOnSurf(); Init(aS, bounds, tol)
            let mut proj = ProjectOnSurface {
                surf: surf.clone(),
                uv_bounds,
                tolerance: 1e-12,
                last_point: None,
                last_uv: None,
                last_distance: f64::MAX,
            };
            proj.init(surf, uv_bounds, 1e-12);
            self.proj_ps_cache.insert(fi, proj);
        }
        self.proj_ps_cache.get_mut(&fi).unwrap()
    }

    // ====================================================================
    // SurfaceAdaptor — OCCT BRepAdaptor_Surface per face
    // ====================================================================

    /// OCCT IntTools_Context::SurfaceAdaptor (IntTools_Context.cxx L327-339).
    /// Returns a cached surface adaptor for the face `fi`.
    pub fn surface_adaptor(&mut self, ds: &DS, fi: usize) -> &mut SurfaceAdaptor {
        if !self.surf_adapt_cache.contains_key(&fi) {
            let surf = ds.face_surface(fi)
                .expect("SurfaceAdaptor: face has no surface");
            let adapt = SurfaceAdaptor::new(surf.clone());
            self.surf_adapt_cache.insert(fi, adapt);
        }
        self.surf_adapt_cache.get_mut(&fi).unwrap()
    }

    // ====================================================================
    // IsPointInFace — OCCT IntTools_FClass2d / IsPointInFace
    // ====================================================================

    /// OCCT IntTools_Context::IsPointInFace (IntTools_Context.hxx L155).
    /// Returns true if the 2D point (U,V) is inside the face `fi`.
    ///
    /// OCCT uses IntTools_FClass2d (IntTools_FClass2d.hxx) which builds a
    /// per-wire UV polygon of the face boundary and classifies with
    /// CSLib_Class2d. Excludes ON points (`aState != TopAbs_OUT && != ON`).
    pub fn is_point_in_face(&self, ds: &DS, fi: usize, uv: DVec2) -> bool {
        let class2d = FClass2d::new(ds, fi, Self::classifier_tol(ds, fi));
        let state = class2d.perform(ds, uv, true);
        state != State::Out && state != State::On
    }

    /// OCCT IntTools_Context::IsPointInOnFace (IntTools_Context.cxx L640-644):
    /// `aState != TopAbs_OUT` — boundary (ON) points are considered inside.
    /// Used by EF's IsProjectable (IsValidPointForFace), unlike IsPointInFace
    /// (used by VF) which excludes ON points.
    pub fn is_point_in_on_face(&self, ds: &DS, fi: usize, uv: DVec2) -> bool {
        let class2d = FClass2d::new(ds, fi, Self::classifier_tol(ds, fi));
        let state = class2d.perform(ds, uv, true);
        state != State::Out
    }

    /// OCCT BRep_Tool::Tolerance(aF) — the face tolerance used to build the
    /// FClass2d classifier (IntTools_Context::FClass2d, IntTools_Context.cxx
    /// L225-242). Floored at CONFUSION so the ON tolerance is non-zero.
    fn classifier_tol(ds: &DS, fi: usize) -> f64 {
        ds.face_tolerance(fi).max(rcad_kernel::CONFUSION)
    }

    // ====================================================================
    // UVBounds — OCCT UVBounds
    // ====================================================================

    /// OCCT IntTools_FClass2d::IsHole — checks if the face wire is a hole.
    /// Uses BRepTopAdaptor_FClass2d (brep_top_adaptor) for classification.
    pub fn fclass2d_is_hole(&self, ds: &DS, fi: usize, _surf: &rcad_kernel::geom::Surface3) -> bool {
        // OCCT: create BRepTopAdaptor_FClass2d(aF, Tol)
        FClass2d::new(ds, fi, 1e-7).is_hole()
    }

    /// OCCT IntTools_Context::SolidClassifier (IntTools_Context.cxx L312-322).
    /// Returns a point-in-solid classifier.
    /// rcad: delegates to brep_class3d::SolidClassifier.
    pub fn solid_classifier_perform(&self, ds: &DS, solid_idx: usize, point: DVec3, tol: f64) -> u8 {
        let si = ds.shape_info(solid_idx);
        let s_shape = si.shape.clone();

        // Collect face indices from solid sub-shapes for classification
        let mut explorer = SolidExplorer::new();
        for &shi in &si.sub_shapes {
            if shi >= ds.nb_shapes() { continue; }
            let sh_info = ds.shape_info(shi);
            if sh_info.shape_type != rcad_kernel::topods::ShapeType::Shell { continue; }
            for &fi in &sh_info.sub_shapes {
                if fi >= ds.nb_shapes() { continue; }
                if ds.shape_info(fi).shape_type == rcad_kernel::topods::ShapeType::Face {
                    explorer.add_face_index(fi);
                }
            }
        }

        // OCCT: create BRepClass3d_SolidClassifier with the solid
        let mut clsf = SolidClassifier::from_shape(&s_shape);
        clsf.explorer = explorer;

        // OCCT: SolidClassifier::Perform(P, Tol)
        clsf.perform(point, tol);
        clsf.my_state
    }

    /// OCCT IntTools_Context::UVBounds (IntTools_Context.cxx L220).
    /// Returns the UV boundaries of the face `fi`.
    pub fn uv_bounds(&mut self, ds: &DS, fi: usize) -> [f64; 4] {
        if !self.uv_bounds_cache.contains_key(&fi) {
            let bounds = ds.face_uv_boundary(fi);
            self.uv_bounds_cache.insert(fi, bounds);
        }
        self.uv_bounds_cache[&fi]
    }
}

/// OCCT IntTools_BeanFaceIntersector::FastComputeAnalytic (L692-816).
/// Returns (isCoincide, hasIntersection). If isCoincide, the whole range is a
/// common part; if !hasIntersection, there is no common part; otherwise the
/// general algorithm continues.
fn fast_compute_analytic(curve: &Curve3, surf: &Surface3, a_criteria: f64) -> (bool, bool) {
    let mut is_coincide = false;
    let mut has_intersection = true;

    // Plane - Circle/Ellipse/Hyperbola/Parabola
    if let Surface3::Plane(surf_plane) = surf {
        let (a_dir, a_p_loc): (glam::DVec3, glam::DVec3) = match curve {
            Curve3::Circle(c) => (c.normal, c.center),
            Curve3::Ellipse(e) => (e.normal, e.center),
            Curve3::Hyperbola(h) => (h.normal, h.center),
            Curve3::Parabola(p) => (p.normal, p.vertex),
            _ => return (false, true),
        };
        let an_angle = a_dir.angle_between(surf_plane.normal);
        if an_angle > rcad_kernel::precision::ANGULAR {
            return (false, true);
        }
        has_intersection = false;
        let a_dist = ((a_p_loc - surf_plane.origin).dot(surf_plane.normal)).abs();
        is_coincide = a_dist < a_criteria;
    }
    // Cylinder - Line/Circle
    else if let Surface3::Cylinder(cyl) = surf {
        let cyl_dir = cyl.axis;
        let cyl_radius = cyl.radius;
        if let Curve3::Line(l) = curve {
            let ldir = l.direction.normalize_or_zero();
            if ldir.cross(cyl_dir).length() > rcad_kernel::precision::ANGULAR {
                return (false, true);
            }
            has_intersection = false;
            // distance from the line to the cylinder axis minus radius
            let axis_pt = cyl.origin;
            let d_line = (l.origin - axis_pt).cross(ldir).length();
            let a_dist = (d_line - cyl_radius).abs();
            is_coincide = a_dist < a_criteria;
        } else if let Curve3::Circle(c) = curve {
            let an_angle = cyl_dir.angle_between(c.normal);
            if an_angle > rcad_kernel::precision::ANGULAR {
                return (false, true);
            }
            // distance from the cylinder axis to the circle center + radius diff
            let a_dist_loc = (c.center - cyl.origin).cross(cyl_dir).length();
            let a_dist = a_dist_loc + (c.radius - cyl_radius).abs();
            is_coincide = a_dist < a_criteria;
            if !is_coincide {
                has_intersection = (a_dist_loc - (c.radius + cyl_radius)) < a_criteria
                    && ((c.radius - cyl_radius).abs() - a_dist_loc) < a_criteria;
            }
        }
    }
    // Sphere - Line
    else if let Surface3::Sphere(sph) = surf {
        if let Curve3::Line(l) = curve {
            let d_line = (l.origin - sph.center).cross(l.direction.normalize()).length();
            let a_dist = d_line - sph.radius;
            has_intersection = a_dist < a_criteria;
        }
    }
    (is_coincide, has_intersection)
}

/// Multi-resolution refinement of the minimum-distance point near `t0`.
fn refine_min(curve: &Curve3, surf: &Surface3, t0: f64, dt: f64, a_t1: f64, a_t2: f64) -> (f64, f64) {
    let mut t_min = t0;
    let mut best_d = f64::MAX;
    let mut win = dt.max((a_t2 - a_t1) * 0.01);
    for _ in 0..8 {
        let lo = (t_min - win).max(a_t1);
        let hi = (t_min + win).min(a_t2);
        if hi - lo <= 0.0 { break; }
        let n_f = 64usize;
        let d_f = (hi - lo) / n_f as f64;
        for i in 0..=n_f {
            let m = lo + d_f * i as f64;
            let pm = curve.point_at(m);
            let (_, qm) = crate::bop::closest_point_on_surface(surf, pm);
            let dm = (qm - pm).length();
            if dm < best_d { best_d = dm; t_min = m; }
        }
        win *= 0.1;
    }
    (t_min, best_d)
}

/// OCCT IntTools_EdgeFace::MakeType (L304-359): VERTEX/EDGE type of a common part.
fn make_type_ef(curve: &Curve3, af1: f64, al1: f64, a_t1: f64, a_t2: f64, a_criteria: f64) -> (f64, f64, bool) {
    let p_f = curve.point_at(af1);
    let p_l = curve.point_at(al1);
    let df1 = (p_f - p_l).length();
    let curve_res = curve_res_ef(curve, a_criteria);
    let is_whole = (af1 - a_t1).abs() < curve_res && (al1 - a_t2).abs() < curve_res;
    if df1 > a_criteria * 2.0 && is_whole {
        return (af1, al1, true); // EDGE
    }
    if is_whole {
        // OCCT L338-347: for a whole-range part with degenerate endpoints
        // (e.g. a full circle: P(first) == P(last)), check the midpoint —
        // if it is also far from the start, the part spans the whole range → EDGE.
        let tm = (af1 + al1) * 0.5;
        if (p_f - curve.point_at(tm)).length() > a_criteria * 2.0 {
            return (af1, al1, true); // EDGE
        }
    }
    (af1, al1, false) // VERTEX
}

/// OCCT IntTools_Tools::ComputeIntRange (for line/plane range width).
fn compute_int_range_ef(tol1: f64, tol2: f64, angle: f64) -> f64 {
    if (std::f64::consts::FRAC_PI_2 - angle).abs() < rcad_kernel::precision::ANGULAR {
        tol2
    } else {
        let a = if angle > std::f64::consts::FRAC_PI_2 { std::f64::consts::PI - angle } else { angle };
        tol1 * (std::f64::consts::FRAC_PI_2 - a).tan() + tol2 / a.sin()
    }
}

/// Insert a marked range into the range manager, splitting overlapping ranges.
fn insert_range_ef(ranges: &mut Vec<(f64, f64, i32)>, s: f64, e: f64, flag: i32) {
    let mut new_ranges: Vec<(f64, f64, i32)> = Vec::new();
    for &(rf, rl, fl) in ranges.iter() {
        if e <= rf || s >= rl {
            new_ranges.push((rf, rl, fl));
        } else if s <= rf && e >= rl {
            new_ranges.push((rf, rl, flag));
        } else {
            // partial overlap: split
            if rf < s {
                new_ranges.push((rf, s, fl));
            }
            let a = s.max(rf);
            let b = e.min(rl);
            if a < b {
                new_ranges.push((a, b, flag));
            }
            if e < rl {
                new_ranges.push((e, rl, fl));
            }
        }
    }
    new_ranges.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    *ranges = new_ranges;
}

/// OCCT BRepAdaptor_Curve::Resolution — parameter step for a 3D tolerance.
fn curve_res_ef(curve: &Curve3, tol: f64) -> f64 {
    match curve {
        Curve3::Line(_) => tol,
        Curve3::Circle(c) => {
            let dt = tol / c.radius;
            if dt <= 1.0 { 2.0 * dt.asin() } else { std::f64::consts::TAU }
        }
        _ => tol,
    }
}
