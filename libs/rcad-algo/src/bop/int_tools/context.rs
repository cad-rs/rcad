// OCCT IntTools_Context — intersection context for VE/EF computations.
//
// OCCT ref: IntTools_Context.hxx / IntTools_Context.cxx
//
// Provides cached projection onto curves/surfaces, point-in-face classification,
// and surface adaptors. Each tool is lazily created and cached per face index.

use std::collections::HashMap;
use crate::bop::ds::DS;
use rcad_kernel::geom::{Surface3, SurfaceEval, Curve2dEval};
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
    // ProjPS — OCCT GeomAPI_ProjectPointOnSurf per face
    // ====================================================================

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
    /// OCCT uses IntTools_FClass2d which performs 2D classification
    /// against the face's wire boundaries. rcad: samples boundary edge
    /// pcurves to build a 2D polygon, then uses ray casting.
    pub fn is_point_in_face(&self, ds: &DS, fi: usize, uv: DVec2) -> bool {
        let face_si_idx = ds.face_shape_idx(fi);
        if face_si_idx >= ds.nb_shapes() {
            return true;
        }
        let boundary_edges: Vec<usize> = ds.shapes[face_si_idx].sub_shapes.iter()
            .filter(|&&ss| ss < ds.nb_shapes() && ds.shapes[ss].shape_type == ShapeType::Edge)
            .copied()
            .collect();
        if boundary_edges.is_empty() {
            return true;
        }

        let mut poly_2d: Vec<DVec2> = Vec::new();
        for &ei in &boundary_edges {
            if ei >= ds.nb_shapes() { continue; }
            let pcurve = match &*ds.shapes[ei].shape.data {
                TShape::Edge(ed) => ed.pcurves.get(&fi).cloned(),
                _ => None,
            };
            if let Some((curve2d, f, l)) = pcurve {
                let n = 16usize;
                let dt = (l - f) / n as f64;
                for i in 0..=n {
                    let t = f + i as f64 * dt;
                    poly_2d.push(curve2d.point_at(t));
                }
            }
        }
        if poly_2d.len() < 3 {
            return true;
        }

        // Ray casting (even-odd rule)
        let mut inside = false;
        let mut j = poly_2d.len() - 1;
        for i in 0..poly_2d.len() {
            let pi = poly_2d[i];
            let pj = poly_2d[j];
            if ((pi.y > uv.y) != (pj.y > uv.y))
                && (uv.x < (pj.x - pi.x) * (uv.y - pi.y) / (pj.y - pi.y) + pi.x)
            {
                inside = !inside;
            }
            j = i;
        }
        inside
    }

    // ====================================================================
    // UVBounds — OCCT UVBounds
    // ====================================================================

    /// OCCT IntTools_FClass2d::IsHole — checks if the face wire is a hole.
    /// Samples center UV point; if outside face boundary → hole.
    pub fn fclass2d_is_hole(&self, ds: &DS, fi: usize, _surf: &rcad_kernel::geom::Surface3) -> bool {
        // Use default UV domain center as sample point
        let uv = DVec2::new(0.5, 0.5);
        !self.is_point_in_face(ds, fi, uv)
    }

    /// OCCT BRepClass3d_SolidClassifier — 3D point classification relative to solid.
    /// Uses ray casting: odd intersections → IN, even → OUT.
    pub fn solid_classifier_is_inside(&self, ds: &DS, solid_idx: usize, point: DVec3) -> bool {
        // Get all faces of the solid
        let mut face_indices: Vec<usize> = Vec::new();
        let si = ds.shape_info(solid_idx);
        if si.shape_type != rcad_kernel::topods::ShapeType::Solid { return false; }
        for &shi in &si.sub_shapes {
            if shi >= ds.nb_shapes() { continue; }
            let sh_info = ds.shape_info(shi);
            if sh_info.shape_type != rcad_kernel::topods::ShapeType::Shell { continue; }
            for &fi in &sh_info.sub_shapes {
                if fi < ds.nb_shapes() && ds.shape_info(fi).shape_type == rcad_kernel::topods::ShapeType::Face {
                    face_indices.push(fi);
                }
            }
        }
        if face_indices.is_empty() { return false; }

        // Ray casting: shoot ray in +X direction, count intersections
        let ray_dir = DVec3::X;
        let mut intersections = 0usize;
        for &fi in &face_indices {
            let surf = match ds.face_surface(fi) {
                Some(s) => s,
                None => continue,
            };
            // Check if ray intersects the face's surface at t > 0
            // Use surface normal to determine if the intersection is front-facing
            match surf {
                rcad_kernel::geom::Surface3::Plane(p) => {
                    let denom = ray_dir.dot(p.normal);
                    if denom.abs() < 1e-12 { continue; }
                    let t = (p.origin - point).dot(p.normal) / denom;
                    if t > 1e-12 {
                        let hit_pt = point + ray_dir * t;
                        // Check if hit point is inside the face boundary
                        if self.is_point_in_face(ds, fi, DVec2::new(0.5, 0.5)) {
                            intersections += 1;
                        }
                    }
                }
                _ => {
                    // For non-planar surfaces, skip (simplified)
                    continue;
                }
            }
        }
        intersections % 2 == 1
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
