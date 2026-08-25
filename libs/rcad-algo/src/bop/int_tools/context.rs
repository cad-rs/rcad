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
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::topods::{ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
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
        if std::env::var("RCAD_PROJ_DEBUG").is_ok() {
            eprintln!("[PROJ] surf={:?} point=({:.3},{:.3},{:.3}) -> uv=({:.4},{:.4}) proj=({:.3},{:.3},{:.3}) dist={:.4}",
                std::mem::discriminant(&self.surf), point.x, point.y, point.z,
                uv.x, uv.y, proj.x, proj.y, proj.z, (proj - point).length());
        }
        // OCCT ProjPS (IntTools_Context.cxx L257-260): Init(aS, Umin, Usup,
        // Vmin, Vsup) restricts the projection to the face's UV rectangle.
        // For a periodic (closed) U direction the unconstrained solution is
        // wrapped back into the principal range — u=2π+ε and u=ε are the SAME
        // point, so clamping (instead of wrapping) would move the solution to
        // a wrong 3D location. Non-periodic directions are clamped to the
        // boundary (the constrained nearest point lies there).
        use rcad_kernel::geom::SurfaceEval;
        let mut u = uv.x;
        let mut v = uv.y;
        let u0 = self.uv_bounds[0];
        let u1 = self.uv_bounds[1];
        if rcad_kernel::geom::SurfaceEval::is_u_periodic(&self.surf) {
            let period = u1 - u0;
            if period > 0.0 && period.is_finite() {
                u = u0 + (u - u0).rem_euclid(period);
            }
        } else {
            u = u.clamp(u0, u1);
        }
        let v0 = self.uv_bounds[2];
        let v1 = self.uv_bounds[3];
        if rcad_kernel::geom::SurfaceEval::is_v_periodic(&self.surf) {
            // OCCT ProjPS wraps periodic directions (same argument as U): a
            // torus v = -PI/2 and v = 3*PI/2 are the SAME point, clamping would
            // move the solution to a wrong 3D location.
            let period = v1 - v0;
            if period > 0.0 && period.is_finite() {
                v = v0 + (v - v0).rem_euclid(period);
            }
        } else {
            v = v.clamp(v0, v1);
        }
        if u == uv.x && v == uv.y {
            self.last_uv = Some(uv);
            self.last_point = Some(proj);
            self.last_distance = (proj - point).length();
        } else {
            let proj2 = rcad_kernel::geom::SurfaceEval::point_at(&self.surf, u, v);
            self.last_uv = Some(DVec2::new(u, v));
            self.last_point = Some(proj2);
            self.last_distance = (proj2 - point).length();
        }
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
        express: bool, ds: &DS, the_fuzz: f64,
    ) -> (i32, Vec<(f64, f64, bool)>, f64) {
        let curve = match ds.edge_curve(n_e) { Some(c) => c.clone(), None => return (-1, vec![], f64::MAX) };
        let surf = match ds.face_surface(n_f) { Some(s) => s, None => return (-1, vec![], f64::MAX) };
        // OCCT IntTools_EdgeFace::Perform (EdgeFace.cxx L528-549): myCriteria
        //   aFuzz = myFuzzyValue/2; aTolF = Tol(face)+aFuzz; aTolE = Tol(edge)+aFuzz
        //   BSpline/Bezier: diff1>100||diff2>100 -> max(aTolE,aTolF), else 1.5*aTolE+aTolF
        let a_fuzz = the_fuzz / 2.0;
        let a_tol_f = ds.face_tolerance(n_f) + a_fuzz;
        let a_tol_e = ds.edge_tolerance(n_e) + a_fuzz;
        let a_criteria = if matches!(curve, Curve3::BSpline(_) | Curve3::Bezier(_)) {
            let diff1 = a_tol_e / a_tol_f;
            let diff2 = a_tol_f / a_tol_e;
            if diff1 > 100.0 || diff2 > 100.0 {
                a_tol_e.max(a_tol_f)
            } else {
                1.5 * a_tol_e + a_tol_f
            }
        } else {
            a_tol_e + a_tol_f
        };
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

        // OCCT EdgeFace L553-563: Quick coincidence check (myQuickCoincidenceCheck).
        // If the edge's PB vertices are both already on the face, check whether the
        // whole edge lies on the face; if so, emit a single full-range EDGE common
        // part and skip the BeanFaceIntersector entirely.
        if express {
            if self.is_coincident_ef(ds, n_f, &curve, &surf, a_t1, a_t2, a_criteria) {
                return (0, vec![(a_t1, a_t2, true)], f64::MAX);
            }
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
            // OCCT EdgeFace L586: IsProjectable(IntermediatePoint(aRange)) =
            //   IntTools_Context::IsValidPointForFace(aPC, aF, myCriteria):
            //   project; if (Umin > aTol) return false; IsPointInOnFace(aP2D).
            //   Range merging happens inside BeanFaceIntersector (L352-378).
            let t_mid = (rf + rl) * 0.5;
            let p_mid = curve.point_at(t_mid);
            let (a_dist, a_u, a_v) = {
                let proj = self.proj_ps(ds, n_f);
                proj.perform(p_mid);
                if proj.nb_points() > 0 {
                    let (u, v) = proj.lower_distance_parameters();
                    (proj.lower_distance(), u, v)
                } else {
                    (f64::MAX, 0.0, 0.0)
                }
            };
            // OCCT L595-608: MakeType runs after the projection filter, per common part.
            let (_, _, is_edge) = make_type_ef(&curve, rf, rl, a_t1, a_t2, a_criteria);
            let in_face = a_dist <= a_criteria
                && self.is_point_in_on_face(ds, n_f, glam::DVec2::new(a_u, a_v));
            if std::env::var("RCAD_EE_DEBUG").is_ok() {
                eprintln!("[EF-DBG]   range=[{:.5},{:.5}] mid={:.5} uv=({:.4},{:.4}) dist={:.3e} inFace={} ty={}", rf, rl, t_mid, a_u, a_v, a_dist, in_face, if is_edge { "E" } else { "V" });
            }
            if !in_face {
                continue;
            }
            common_parts.push((rf, rl, is_edge));
        }
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            eprintln!("[EF-DBG] compute_ef e={} f={} r=[{:.4},{:.4}] kept={:?}", n_e, n_f, a_t1, a_t2,
                common_parts.iter().map(|&(s, e, ty)| format!("[{:.5},{:.5}]{}", s, e, if ty { "E" } else { "V" })).collect::<Vec<_>>());
        }
        (0, common_parts, min_dist)
    }

    /// OCCT IntTools_EdgeFace::IsCoincident (IntTools_EdgeFace.cxx L62-163).
    /// Quick coincidence check used when both PB vertices are already known to
    /// be on the face. Samples aNbSeg+1 points along the (boundary-shifted)
    /// edge range; if enough of them are within myCriteria of the face surface
    /// and classify as In/On (not Out), the edge is coincident with the face.
    fn is_coincident_ef(&mut self, ds: &DS, n_f: usize, curve: &Curve3, surf: &Surface3,
                        a_t1: f64, a_t2: f64, a_criteria: f64) -> bool {
        let a_nb_seg: usize = if matches!(curve, Curve3::Line(_)) && matches!(surf, Surface3::Plane(_)) {
            2 // Line/Plane: check only three points
        } else {
            23
        };
        // OCCT: aTreshIdxF = RealToInt((aNbSeg+1)*0.25), aTreshIdxL = RealToInt((aNbSeg+1)*0.75)
        let a_tresh_idx_f = ((a_nb_seg as f64 + 1.0) * 0.25).round() as usize;
        let a_tresh_idx_l = ((a_nb_seg as f64 + 1.0) * 0.75).round() as usize;
        // Shift the sample range in from the boundaries by 1% to avoid projection
        // on the surface boundary (OCCT L86-90).
        let a_bnd_shift = 0.01 * (a_t2 - a_t1);
        let a_t1 = a_t1 + a_bnd_shift;
        let a_t2 = a_t2 - a_bnd_shift;
        let d_t = (a_t2 - a_t1) / a_nb_seg as f64;
        let mut is_classified = false;
        let mut i_cnt: usize = 0;
        let class2d = FClass2d::new(ds, n_f, Self::classifier_tol(ds, n_f));
        for i in 0..=a_nb_seg {
            let a_t = a_t1 + (i as f64) * d_t;
            let a_p = curve.point_at(a_t);
            let (a_d, a_u, a_v) = {
                let proj = self.proj_ps(ds, n_f);
                proj.perform(a_p);
                if proj.nb_points() > 0 {
                    let (u, v) = proj.lower_distance_parameters();
                    (proj.lower_distance(), u, v)
                } else {
                    (f64::MAX, 0.0, 0.0)
                }
            };
            if a_d == f64::MAX {
                // OCCT L101-104: !aProjector.IsDone() -> continue
                continue;
            }
            if a_d > a_criteria {
                if a_d > 100.0 * a_criteria {
                    return false;
                }
                continue;
            }
            i_cnt += 1;
            // Only the begin, end and middle points are classified
            if ((0 < i) && (i < a_tresh_idx_f)) || ((a_tresh_idx_l < i) && (i < a_nb_seg)) {
                continue;
            }
            if is_classified && (i != a_nb_seg) {
                continue;
            }
            let state = class2d.perform(ds, glam::DVec2::new(a_u, a_v), true);
            if state == State::Out {
                return false;
            }
            if i != 0 {
                is_classified = true;
            }
        }
        let a_coeff = i_cnt as f64 / (a_nb_seg as f64 + 1.0);
        a_coeff > 0.5
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

    // ====================================================================
    // ComputePE — OCCT IntTools_Context::ComputePE (IntTools_Context.cxx L437-495)
    // ====================================================================

    /// OCCT IntTools_Context::ComputePE (IntTools_Context.cxx L437-495):
    ///   project a point onto an edge; returns 0 (ok), -2 (not geometric),
    ///   -4 (distance exceeds tolerance), -3 (no projection).
    ///   rcad: `closest_point_on_curve_range` always returns a clamped
    ///   projection, so the OCCT "point falls out of the curve" branch
    ///   (distance to the edge vertices) is covered by the clamp.
    pub fn compute_pe(
        &mut self,
        a_p1: DVec3,
        a_tol_p1: f64,
        n_e: usize,
        ds: &DS,
        a_t: &mut f64,
        a_dist: &mut f64,
    ) -> i32 {
        // OCCT L443-446: if (!BRep_Tool::IsGeometric(aE2)) return -2;
        let curve = match ds.edge_curve(n_e) { Some(c) => c.clone(), None => return -2 };
        // OCCT ProjPC (IntTools_Context.cxx L276-281): projection limited to the
        // edge's parameter range.
        let a_range = ds.edge_range(n_e);
        let a_proj = closest_point_on_curve_range(&curve, a_p1, a_range[0], a_range[1], 128);
        // OCCT L457-467: point falls on the curve.
        let a_dist_v = a_proj.distance;
        let a_tol_e = ds.edge_tolerance(n_e);
        let a_tol_sum = a_tol_p1 + a_tol_e + rcad_kernel::CONFUSION;
        if a_dist_v > a_tol_sum {
            return -4;
        }
        *a_t = a_proj.param;
        *a_dist = a_dist_v;
        0
    }

    // ====================================================================
    // IsVertexOnLine — OCCT IntTools_Context::IsVertexOnLine
    // (IntTools_Context.cxx L776-983)
    // ====================================================================

    /// OCCT IntTools_Context::IsVertexOnLine (IntTools_Context.cxx L776-983).
    /// Returns true when the vertex (point aPV, tolerance aTolV) lies on the
    /// section curve aCurve within aTolC; the parameter is returned in aT.
    pub fn is_vertex_on_line(
        &mut self,
        a_pv: DVec3,
        a_tol_v: f64,
        a_curve: &Curve3,
        a_tol_c: f64,
        a_t: &mut f64,
        a_domain: [f64; 2],
    ) -> bool {
        // OCCT L788-809: tolerance sum depending on the curve type.
        let [a_first, a_last] = a_domain;
        let mut a_tol_sum = a_tol_v + a_tol_c;
        let is_spline = matches!(a_curve, Curve3::BSpline(_) | Curve3::Bezier(_));
        a_tol_sum = 2.0 * a_tol_sum;
        if a_tol_sum < if is_spline { 1.0e-5 } else { 1.0e-6 } {
            a_tol_sum = if is_spline { 1.0e-5 } else { 1.0e-6 };
        }
        let is_inf = |v: f64| v.is_infinite();
        // OCCT L814-874: checking extremities first.
        let mut b_first_valid = false;
        let mut a_first_dist = f64::INFINITY;
        if !is_inf(a_first) {
            let a_pc_first = a_curve.point_at(a_first);
            a_first_dist = a_pv.distance(a_pc_first);
            if a_first_dist < a_tol_sum {
                b_first_valid = true;
                *a_t = a_first;
                if a_first_dist > a_tol_v {
                    // OCCT L829: Extrema_LocateExtPC(aPv, aGAC, aFirst, 1.e-10)
                    if let Some(poc) = rcad_kernel::base::extrema::extrema_locate_ext_pc(
                        a_pv, a_curve, a_first, a_first, a_last, 1e-10,
                    ) {
                        // OCCT L836-840: validate local result
                        let mid = (a_last + a_first) * 0.5;
                        if (poc.param > mid)
                            || (a_pv.distance(poc.point) > a_tol_sum)
                            || (a_pc_first.distance(poc.point) < rcad_kernel::CONFUSION)
                        {
                            *a_t = a_first;
                        } else {
                            *a_t = poc.param;
                        }
                    } else {
                        // OCCT L842-870: Extrema_LocateExtPC failed -> global fallback (Extrema_ExtPC)
                        let ext = rcad_kernel::base::extrema::ExtPC::new(
                            a_pv, a_curve, 1e-10, a_first, a_last,
                        );
                        if ext.is_done() {
                            let mut a_min_dist = f64::INFINITY;
                            let mut a_min_idx = None;
                            for i in 1..=ext.nb_ext() {
                                let sq_d = ext.square_distance(i);
                                if sq_d < a_min_dist {
                                    a_min_dist = sq_d;
                                    a_min_idx = Some(i);
                                }
                            }
                            if let Some(idx) = a_min_idx {
                                let poc = ext.point(idx);
                                let mid = (a_last + a_first) * 0.5;
                                if (poc.param > mid)
                                    || (a_pv.distance(poc.point) > a_tol_sum)
                                    || (a_pc_first.distance(poc.point) < rcad_kernel::CONFUSION)
                                {
                                    *a_t = a_first;
                                } else {
                                    *a_t = poc.param;
                                }
                            }
                        }
                    }
                }
            }
        }
        // OCCT L876-941: last extremity.
        if !is_inf(a_last) {
            let a_pc_last = a_curve.point_at(a_last);
            let a_dist = a_pv.distance(a_pc_last);
            if b_first_valid && (a_first_dist < a_dist) {
                return true;
            }
            if a_dist < a_tol_sum {
                *a_t = a_last;
                if a_dist > a_tol_v {
                    // OCCT L890: Extrema_LocateExtPC(aPv, aGAC, aLast, 1.e-10)
                    if let Some(poc) = rcad_kernel::base::extrema::extrema_locate_ext_pc(
                        a_pv, a_curve, a_last, a_first, a_last, 1e-10,
                    ) {
                        // OCCT L897-901: validate local result
                        let mid = (a_last + a_first) * 0.5;
                        if (poc.param < mid)
                            || (a_pv.distance(poc.point) > a_tol_sum)
                            || (a_pc_last.distance(poc.point) < rcad_kernel::CONFUSION)
                        {
                            *a_t = a_last;
                        } else {
                            *a_t = poc.param;
                        }
                    } else {
                        // OCCT L905-931: Extrema_LocateExtPC failed -> global fallback (Extrema_ExtPC)
                        let ext = rcad_kernel::base::extrema::ExtPC::new(
                            a_pv, a_curve, 1e-10, a_first, a_last,
                        );
                        if ext.is_done() {
                            let mut a_min_dist = f64::INFINITY;
                            let mut a_min_idx = None;
                            for i in 1..=ext.nb_ext() {
                                let sq_d = ext.square_distance(i);
                                if sq_d < a_min_dist {
                                    a_min_dist = sq_d;
                                    a_min_idx = Some(i);
                                }
                            }
                            if let Some(idx) = a_min_idx {
                                let poc = ext.point(idx);
                                let mid = (a_last + a_first) * 0.5;
                                if (poc.param < mid)
                                    || (a_pv.distance(poc.point) > a_tol_sum)
                                    || (a_pc_last.distance(poc.point) < rcad_kernel::CONFUSION)
                                {
                                    *a_t = a_last;
                                } else {
                                    *a_t = poc.param;
                                }
                            }
                        }
                    }
                }
                return true;
            }
        } else if b_first_valid {
            return true;
        }
        // OCCT L943-982: projection onto the full curve (ProjPT).
        let proj = closest_point_on_curve_range(a_curve, a_pv, a_first, a_last, 64);
        let a_dist = proj.distance;
        if a_dist > a_tol_sum {
            return false;
        }
        *a_t = proj.param;
        true
    }

    // ====================================================================
    // IsValidPointForFace(s) — OCCT IntTools_Context::IsValidPointForFace
    // (IntTools_Context.cxx L648-674, L678-692)
    // ====================================================================

    /// OCCT IntTools_Context::IsValidPointForFace (IntTools_Context.cxx L648-674).
    pub fn is_valid_point_for_face(&mut self, a_p: DVec3, n_f: usize, ds: &DS, a_tol: f64) -> bool {
        // OCCT L655-657: GeomAPI_ProjectPointOnSurf& aProjector = ProjPS(aF); Perform(aP).
        let (done, umin, uv) = {
            let proj = self.proj_ps(ds, n_f);
            proj.perform(a_p);
            (
                proj.nb_points() > 0,
                proj.lower_distance(),
                proj.lower_distance_parameters(),
            )
        };
        if !done {
            return false;
        }
        // OCCT L662-667: if (Umin > aTol) return !bFlag; (false)
        if umin > a_tol {
            return false;
        }
        // OCCT L669-671: IsPointInOnFace(aF, aP2D(U, V)).
        let r = self.is_point_in_on_face(ds, n_f, DVec2::new(uv.0, uv.1));
        r
    }

    /// OCCT IntTools_Context::IsValidPointForFaces (IntTools_Context.cxx L678-692).
    pub fn is_valid_point_for_faces(
        &mut self,
        a_p: DVec3,
        n_f1: usize,
        n_f2: usize,
        ds: &DS,
        a_tol: f64,
    ) -> bool {
        let b_flag1 = self.is_valid_point_for_face(a_p, n_f1, ds, a_tol);
        if !b_flag1 {
            return b_flag1;
        }
        self.is_valid_point_for_face(a_p, n_f2, ds, a_tol)
    }

    // ====================================================================
    // IsValidBlockForFace(s) — OCCT IntTools_Context::IsValidBlockForFaces
    // (IntTools_Context.cxx L696-756)
    // ====================================================================

    /// OCCT IntTools_Context::IsValidBlockForFace (IntTools_Context.cxx L696-714).
    pub fn is_valid_block_for_face(
        &mut self,
        a_t1: f64,
        a_t2: f64,
        a_c: &crate::bop::int_tools::face_face::IntersectionCurve,
        n_f: usize,
        ds: &DS,
        a_tol: f64,
    ) -> bool {
        // OCCT L706-712: aTInterm = IntermediatePoint; aPInterm = aC3D->D0(aTInterm).
        let a_t_interm = crate::bop::int_tools::face_make_curve::intermediate_point(a_t1, a_t2);
        let a_p_interm = a_c.curve.point_at(a_t_interm);
        self.is_valid_point_for_face(a_p_interm, n_f, ds, a_tol)
    }

    /// OCCT IntTools_Context::IsValidBlockForFaces (IntTools_Context.cxx L718-756).
    pub fn is_valid_block_for_faces(
        &mut self,
        the_t1: f64,
        the_t2: f64,
        the_c: &crate::bop::int_tools::face_face::IntersectionCurve,
        the_f1: usize,
        the_f2: usize,
        ds: &DS,
        the_tol: f64,
    ) -> bool {
        let a_mid_par = crate::bop::int_tools::face_make_curve::intermediate_point(the_t1, the_t2);
        let a_p = the_c.curve.point_at(a_mid_par);
        // OCCT IntTools_Context::IsValidBlockForFaces (IntTools_Context.cxx
        // L717-754): for each face use the curve's 2D pcurve for the mid-point
        // classification when available (IsPointInOnFace — ON counts as valid),
        // otherwise project the 3D point (IsValidPointForFace). The 2D path is
        // decisive for section blocks lying outside the face's UV domain.
        let mut b_flag = true;
        for (a_pc, a_f) in [
            (the_c.pcurve1.as_ref(), the_f1),
            (the_c.pcurve2.as_ref(), the_f2),
        ] {
            if !b_flag {
                break;
            }
            b_flag = match a_pc {
                Some(pc) => {
                    let a_pnt2d = rcad_kernel::geom::Curve2dEval::point_at(pc, a_mid_par);
                    self.is_point_in_on_face(ds, a_f, a_pnt2d)
                }
                None => self.is_valid_point_for_face(a_p, a_f, ds, the_tol),
            };
        }
        b_flag
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
    /// L225-242; aTolF = BRep_Tool::Tolerance(aFF), where aFF carries the DS
    /// face tolerance aTol = BRep_Tool::Tolerance(myFace),
    /// BOPAlgo_BuilderFace.cxx L396).
    fn classifier_tol(ds: &DS, fi: usize) -> f64 {
        ds.face_tolerance(fi)
    }

    // ====================================================================
    // UVBounds — OCCT UVBounds
    // ====================================================================

    /// OCCT IntTools_FClass2d::IsHole — checks if the wire is a hole.
    /// Uses BRepTopAdaptor_FClass2d (brep_top_adaptor) for classification.
    /// `loop_edges` is the analyzed loop wire of a temporary face
    /// (BOPAlgo_BuilderFace.cxx L437-445: aBB.MakeFace(aFace, aS, aLoc, aTol);
    /// aBB.Add(aFace, aWire); IntTools_Context::FClass2d(aFace)).
    pub fn fclass2d_is_hole(&self, ds: &DS, fi: usize, loop_edges: &[Shape]) -> bool {
        // OCCT IntTools_Context::FClass2d(aF) -- IntTools_Context.cxx L225-242,
        // aTolF = BRep_Tool::Tolerance(aFF).
        FClass2d::new_for_loop(ds, fi, Self::classifier_tol(ds, fi), loop_edges).is_hole()
    }

    /// OCCT IntTools_Context::SolidClassifier (IntTools_Context.cxx L312-322).
    /// Returns a point-in-solid classifier.
    /// rcad: delegates to brep_class3d::SolidClassifier, which explores the
    /// solid shape's faces directly (OCCT BRepClass3d explores the TopoDS
    /// shape, never BOPDS).
    pub fn solid_classifier_perform(&self, ds: &DS, solid_idx: usize, point: DVec3, tol: f64) -> u8 {
        let si = ds.shape_info(solid_idx);
        let s_shape = si.shape.clone();

        // OCCT: create BRepClass3d_SolidClassifier with the solid
        let mut clsf = SolidClassifier::from_shape(&s_shape);

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

// ============================================================================
// Tests
// ============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Plane, Surface3};

    /// Regression test for audit item #7 — ProjPS restricts the projection to
    /// the face's UV rectangle (OCCT IntTools_Context.cxx L257-260). The
    /// unconstrained nearest point (2,3) on an unbounded plane is clamped to
    /// the [0,1]x[0,1] boundary. The old code returned the out-of-bounds UV
    /// and its 3D point, so this test fails before the fix.
    #[test]
    fn project_on_surface_clamps_to_bounds() {
        let surf = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let mut proj = ProjectOnSurface {
            surf: surf.clone(),
            uv_bounds: [0.0, 1.0, 0.0, 1.0],
            tolerance: 1e-12,
            last_point: None,
            last_uv: None,
            last_distance: f64::MAX,
        };
        proj.perform(DVec3::new(2.0, 3.0, 5.0));
        let (u, v) = proj.lower_distance_parameters();
        assert_eq!(u, 1.0, "U must be clamped to the boundary");
        assert_eq!(v, 1.0, "V must be clamped to the boundary");
        let p = proj.nearest_point();
        assert!(
            (p - DVec3::new(1.0, 1.0, 0.0)).length() < 1e-12,
            "3D point must lie on the clamped UV, got {:?}",
            p
        );
    }

    /// A point whose projection already lies inside the bounds is unchanged.
    #[test]
    fn project_on_surface_inside_bounds_unchanged() {
        let surf = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let mut proj = ProjectOnSurface {
            surf,
            uv_bounds: [0.0, 1.0, 0.0, 1.0],
            tolerance: 1e-12,
            last_point: None,
            last_uv: None,
            last_distance: f64::MAX,
        };
        proj.perform(DVec3::new(0.5, 0.5, 5.0));
        let (u, v) = proj.lower_distance_parameters();
        assert!((u - 0.5).abs() < 1e-12);
        assert!((v - 0.5).abs() < 1e-12);
    }
}
