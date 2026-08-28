// OCCT BRepClass3d_SolidExplorer (BRepClass3d_SolidExplorer.hxx / .cxx)
// Exploration of a BRep Shape for classification.
// Provides face iteration, bounding box rejection, and BVH tree.

use crate::bop::closest_point_on_surface;
use crate::bop::int_tools::bean_face_intersector::{
    BeanFaceIntersector, BRepAdaptorCurve, BRepAdaptorSurface,
};
use crate::topalgo::shape_source::ShapeSource;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Plane as PlaneGeom, Surface3, SurfaceEval};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, TShape};
use glam::{DVec2, DVec3};
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::Arc;

/// A face of the explored solid, with its surface, orientation, and — for
/// planar surfaces — the UV bounding box of the face's wires.
/// OCCT IntCurvesFace_Intersector only counts intersections that fall inside
/// the face's UV domain; without this check the ray-cast would treat the
/// infinite supporting plane as the face.
pub(crate) struct ExplorerFace {
    // face TShape identity (ptr_id, location) — the DS registration key.
    pub(crate) key: (u64, u32),
    // The source face shape (for the IntTools_FClass2d translation when the
    // face has no DS registration — OCCT IntTools_FClass2d::Init takes a
    // TopoDS_Face).
    pub(crate) src: Option<Shape>,
    pub(crate) surf: Surface3,
    pub(crate) ori: Orientation,
    pub(crate) uv_bounds: Option<[f64; 4]>, // [umin, umax, vmin, vmax]
    // 2D domain of the face for the ray-cast UV check (outer polygon +
    // holes), built by projecting the wire vertices onto the surface
    // (OCCT BRepTopAdaptor_TopolTool::Classify — IntCurvesFace_Intersector
    // L256-286 accepts only intersections whose UV is IN/ON the face).
    pub(crate) uv_polys: Option<(Vec<DVec2>, Vec<Vec<DVec2>>)>,
    // Boundary edge pcurves (2D curve, parameter range, edge orientation in
    // the FORWARD-ized face) for the FindAPointInTheFace probing
    // (BRepClass3d_SolidExplorer.cxx L74-167): the probe moves from a
    // boundary-edge point into the face interior along the edge normal.
    pub(crate) boundary: Vec<(Curve2d, [f64; 2], Orientation)>,
    // Boundary edge keys (ptr_id, location) in wire order — for the mapEF
    // (edge -> faces) of the SClassifier (TopExp::MapShapesAndAncestors).
    pub(crate) edge_keys: Vec<(u64, u32)>,
    // Per-face FClass2d cache — OCCT IntTools_Context::myFClass2dMap
    // (IntTools_Context.hxx): the classifier is built once per (face,
    // tolerance) and reused across the UV probe loop.  Key = (face index,
    // tolerance bits); the synthetic faces use face index 0.
    pub(crate) fclass_cache: std::cell::RefCell<HashMap<(usize, u64), std::rc::Rc<crate::topalgo::brep_top_adaptor::fclass2d::FClass2d>>>,
}

impl Clone for ExplorerFace {
    fn clone(&self) -> Self {
        ExplorerFace {
            key: self.key,
            src: self.src.clone(),
            surf: self.surf.clone(),
            ori: self.ori,
            uv_bounds: self.uv_bounds,
            uv_polys: self.uv_polys.clone(),
            boundary: self.boundary.clone(),
            edge_keys: self.edge_keys.clone(),
            fclass_cache: std::cell::RefCell::new(HashMap::new()),
        }
    }
}

/// OCCT BRepClass3d_SolidExplorer — explores a solid's faces for point classification.
pub struct SolidExplorer {
    pub ds: Option<Arc<dyn ShapeSource>>,
    shape: Option<Shape>,
    face_indices: Vec<usize>,
    // Face geometry (surface + orientation) collected from the shape tree.
    // Used for classification without a DS reference (OCCT BRepClass3d
    // explores the TopoDS_Shape directly, never BOPDS).
    pub(crate) face_surfaces: Vec<ExplorerFace>,
    // OCCT aMapEV (BRepClass3d_SClassifier L213): vertices and edges of the
    // solid — a point within their tolerance is ON the boundary.
    vertices: Vec<DVec3>,
    edges: Vec<Option<Curve3>>,
    // Edge tolerances (BRep_Tool::Tolerance(edge)) parallel to `edges`.
    edge_tols: Vec<f64>,
    // Edge parameter ranges (BRep_Tool::Range(edge)) parallel to `edges`.
    edge_ranges: Vec<[f64; 2]>,
    // Vertex tolerances (BRep_Tool::Tolerance(vertex)) parallel to `vertices`.
    vert_tols: Vec<f64>,
    /// Merged TopLoc_Location table (index 0 = identity); empty when the
    /// explorer was built without location support.
    locations: Vec<glam::DAffine3>,
    // OCCT myReject (SolidExplorer.cxx L990): true for a solid without faces
    // (infinite solid) — Reject(P) returns it directly.
    my_reject: bool,
    // OCCT myFirstFace / myParamOnEdge (OtherSegment L493-622): the face
    // iteration cursor and the probing parameter cycled on retry.
    my_first_face: usize,
    my_param_on_edge: f64,
}

impl SolidExplorer {
    pub fn new() -> Self {
        SolidExplorer {
            ds: None,
            shape: None,
            face_indices: Vec::new(),
            face_surfaces: Vec::new(),
            vertices: Vec::new(),
            edges: Vec::new(),
            edge_tols: Vec::new(),
            edge_ranges: Vec::new(),
            vert_tols: Vec::new(),
            locations: vec![glam::DAffine3::IDENTITY],
            my_reject: true,
            my_first_face: 0,
            my_param_on_edge: 0.512345,
        }
    }

    pub fn clear(&mut self) {
        self.shape = None;
        self.face_indices.clear();
        self.face_surfaces.clear();
        self.vertices.clear();
        self.edges.clear();
        self.edge_tols.clear();
        self.edge_ranges.clear();
        self.vert_tols.clear();
    }

    /// OCCT: InitShape(S) — initialize the explorer with a solid shape.
    /// Traverses the shape tree and collects the surfaces + orientations of
    /// all faces, so the point classification does not depend on the DS.
    /// The face order follows the shape tree order (shell faces in storage
    /// order) — OCCT TopExp_Explorer iterates depth-first in positive order
    /// and BRepClass3d_SClassifier::PerformInfinitePoint probes the faces in
    /// this same order (aFaces collection, SClassifier.cxx L148-154); the
    /// first probe decides the state, so the order is semantic.
    pub fn init_shape(&mut self, s: &Shape) {
        self.shape = Some(s.clone());
        self.face_indices.clear();
        self.face_surfaces.clear();
        self.vertices.clear();
        self.edges.clear();
        self.edge_tols.clear();
        self.edge_ranges.clear();
        self.vert_tols.clear();
        // OCCT InitShape L905-910: myFirstFace = 0; myParamOnEdge = 0.512345;
        // myReject = true (no faces yet).
        self.my_first_face = 0;
        self.my_param_on_edge = 0.512345;
        self.my_reject = true;
        // OCCT InitShape L936-967: the EV map (aMapEV) is filled per face —
        // each edge of each face (skipping INTERNAL/EXTERNAL faces/edges and
        // degenerated edges) adds the edge and its vertices. The
        // TopTools_ShapeMapHasher keys on (TShape*, Location), so rcad dedups
        // by (ptr_id, location).
        let mut edge_seen: HashSet<(u64, u32)> = HashSet::new();
        let mut vert_seen: HashSet<(u64, u32)> = HashSet::new();
        let mut queue: VecDeque<(Shape, Orientation)> = VecDeque::new();
        queue.push_back((s.clone(), Orientation::Forward));
        while let Some((sh, cum_or)) = queue.pop_front() {
            match &*sh.data {
                TShape::Solid(sd) => {
                    for x in &sd.shells {
                        queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                    }
                }
                TShape::CompSolid(cd) => {
                    for x in cd {
                        queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                    }
                }
                TShape::Compound(cd) => {
                    for x in cd {
                        queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                    }
                }
                TShape::Shell(sd) => {
                    for x in &sd.faces {
                        queue.push_back((x.clone(), cum_or.compose(sh.orientation)));
                    }
                }
                TShape::Vertex(vd) => {
                    // OCCT MapShapes(aE, myMapEV) adds each edge's vertices;
                    // the map key is (TShape*, Location).
                    if vert_seen.contains(&(sh.ptr_id(), sh.location)) {
                        continue;
                    }
                    vert_seen.insert((sh.ptr_id(), sh.location));
                    let loc = self
                        .locations
                        .get(sh.location as usize)
                        .copied()
                        .unwrap_or(glam::DAffine3::IDENTITY);
                    self.vertices.push(loc.transform_point3(vd.point));
                    self.vert_tols.push(vd.tolerance);
                }
                TShape::Edge(ed) => {
                    // OCCT InitShape L944-963: skip INTERNAL/EXTERNAL edges
                    // and degenerated edges (BRep_Tool::Degenerated).
                    let e_or = cum_or.compose(sh.orientation);
                    if e_or == Orientation::Internal || e_or == Orientation::External {
                        continue;
                    }
                    if ed.degenerated {
                        continue;
                    }
                    if edge_seen.contains(&(sh.ptr_id(), sh.location)) {
                        continue;
                    }
                    edge_seen.insert((sh.ptr_id(), sh.location));
                    let loc = self
                        .locations
                        .get(sh.location as usize)
                        .copied()
                        .unwrap_or(glam::DAffine3::IDENTITY);
                    self.edges.push(ed.curve.as_ref().map(|c| rcad_kernel::geom::transform_curve(c, &loc)));
                    self.edge_tols.push(ed.tolerance);
                    self.edge_ranges.push(ed.range);
                    // OCCT TopExp::MapShapes(aE, myMapEV) — the edge's vertices.
                    for v in [&ed.first, &ed.last] {
                        queue.push_back((v.clone(), e_or));
                    }
                }
                TShape::Face(fd) => {
                    // OCCT InitShape L938-947: skip INTERNAL/EXTERNAL faces.
                    let face_or = cum_or.compose(sh.orientation);
                    if face_or == Orientation::Internal || face_or == Orientation::External {
                        continue;
                    }
                    let loc = self
                        .locations
                        .get(sh.location as usize)
                        .copied()
                        .unwrap_or(glam::DAffine3::IDENTITY);
                    let surface = fd.surface.as_ref().map(|s| {
                        if loc == glam::DAffine3::IDENTITY {
                            s.clone()
                        } else {
                            rcad_kernel::geom::transform_surface(s, &loc)
                        }
                    });
                    let uv_bounds = surface.as_ref().and_then(|surf| match surf {
                        Surface3::Plane(pl) => compute_plane_uv_bounds(&sh, pl, &self.locations),
                        _ => None,
                    });
                    let uv_polys = surface.as_ref().map(|surf| build_uv_polys(&sh, surf, &self.locations));
                    let boundary = collect_boundary_pcurves(&sh);
                    let edge_keys = collect_face_edge_keys(&sh);
                    if let Some(surf) = surface {
                        // Compose the accumulated shell orientation into the
                        // face, as OCCT BRepClass3d_SolidExplorer::InitShape
                        // (L920-924) does via TopExp_Explorer with cumOri=true.
                        self.face_surfaces.push(ExplorerFace {
                            key: (sh.ptr_id(), sh.location),
                            src: Some(sh.clone()),
                            surf,
                            ori: face_or,
                            uv_bounds,
                            uv_polys,
                            boundary,
                            edge_keys,
                            fclass_cache: std::cell::RefCell::new(HashMap::new()),
                        });
                        // OCCT InitShape L914-918: at least one face -> the
                        // solid is not a void (myReject = false).
                        self.my_reject = false;
                    }
                    // OCCT InitShape L949-965: walk the face's edges (through
                    // its wires) for the EV map.
                    queue.push_back((fd.outer_wire.clone(), face_or));
                    for w in &fd.inner_wires {
                        queue.push_back((w.clone(), face_or));
                    }
                }
                TShape::Wire(wd) => {
                    // TopExp_Explorer(aF, TopAbs_EDGE) descends through the
                    // wires; the composed orientation accumulates.
                    let w_or = cum_or.compose(sh.orientation);
                    for e in &wd.edges {
                        queue.push_back((e.clone(), w_or));
                    }
                }
                _ => {}
            }
        }
    }

    /// Constructor from Shape.
    pub fn from_shape(s: &Shape) -> Self {
        let mut exp = SolidExplorer {
            ds: None,
            shape: Some(s.clone()),
            face_indices: Vec::new(),
            face_surfaces: Vec::new(),
            vertices: Vec::new(),
            edges: Vec::new(),
            edge_tols: Vec::new(),
            edge_ranges: Vec::new(),
            vert_tols: Vec::new(),
            locations: vec![glam::DAffine3::IDENTITY],
            my_reject: true,
            my_first_face: 0,
            my_param_on_edge: 0.512345,
        };
        exp.init_shape(s);
        exp
    }

    /// Constructor from Shape with the merged TopLoc_Location table — the
    /// located sub-shapes (e.g. OCCT MakePrism folded TShape + Location) are
    /// transformed when collecting vertices/edges/surfaces.
    pub fn from_shape_with_locations(s: &Shape, locations: &[glam::DAffine3]) -> Self {
        let mut exp = SolidExplorer {
            ds: None,
            shape: Some(s.clone()),
            face_indices: Vec::new(),
            face_surfaces: Vec::new(),
            vertices: Vec::new(),
            edges: Vec::new(),
            edge_tols: Vec::new(),
            edge_ranges: Vec::new(),
            vert_tols: Vec::new(),
            locations: if locations.is_empty() {
                vec![glam::DAffine3::IDENTITY]
            } else {
                locations.to_vec()
            },
            my_reject: true,
            my_first_face: 0,
            my_param_on_edge: 0.512345,
        };
        exp.init_shape(s);
        exp
    }

    /// OCCT: Reject(P) — SolidExplorer.cxx L989-992: returns myReject (the
    /// "solid without face" flag set by InitShape). Used by
    /// BRepClass3d_SClassifier::Perform L207-212 (not the bbox fast-path,
    /// which is BRepClass3d_SolidClassifier::Perform L175-181 via
    /// explorer.Box()).
    pub fn reject(&self, _p: DVec3) -> bool {
        self.my_reject
    }

    /// OCCT: myReject — true when the solid has no faces (infinite solid).
    pub fn is_rejected(&self) -> bool {
        self.my_reject
    }

    /// OCCT BRepClass3d_SolidClassifier::Perform L175-181: explorer.Box()
    /// bounding-box fast rejection. rcad: built from the collected vertices
    /// (semantic equivalent of BRepBndLib on the shape).
    pub fn box_is_out(&self, p: DVec3) -> bool {
        if self.vertices.is_empty() && self.edges.is_empty() {
            return false;
        }
        let mut min = DVec3::splat(f64::INFINITY);
        let mut max = DVec3::splat(f64::NEG_INFINITY);
        for v in &self.vertices {
            min = min.min(*v);
            max = max.max(*v);
        }
        for e in &self.edges {
            if let Some(curve) = e {
                // Sample the curve extents like BRepBndLib::Add.
                for k in 0..=8 {
                    let t = k as f64 / 8.0;
                    let p = curve.point_at(t);
                    min = min.min(p);
                    max = max.max(p);
                }
            }
        }
        p.x < min.x - 1e-7 || p.x > max.x + 1e-7
            || p.y < min.y - 1e-7 || p.y > max.y + 1e-7
            || p.z < min.z - 1e-7 || p.z > max.z + 1e-7
    }

    /// OCCT BRepClass3d_SClassifier::Perform (L212-228) — a point within the
    /// tolerance of a VERTEX or EDGE of the solid is ON the boundary
    /// (aMapEV tree selection); plus (L405-449) a point within the tolerance
    /// of a face's SURFACE whose closest UV is IN/ON the face's 2D domain is
    /// ON (the parallel-line branch checks the point-vs-surface distance via
    /// Extrema_ExtPS and the UV via ClassifyUVPoint). rcad: vertices by point
    /// distance, edges by curve distance, faces by surface distance + UV
    /// domain.
    pub fn point_on_face(&self, p: DVec3, tol: f64) -> bool {
        for v in &self.vertices {
            if (p - *v).length() <= tol {
                return true;
            }
        }
        for e in &self.edges {
            if let Some(curve) = e {
                if curve_point_distance(curve, p) <= tol {
                    return true;
                }
            }
        }
        for f in &self.face_surfaces {
            if let Surface3::Plane(pl) = &f.surf {
                // Planar faces: the projection is a single dot product — the
                // full closest_point_on_surface dispatch would be wasted here.
                let d = (p - pl.origin).dot(pl.normal).abs();
                if d <= tol {
                    let rel = p - pl.origin - pl.normal * (p - pl.origin).dot(pl.normal);
                    let uv = DVec2::new(rel.dot(pl.u_dir), rel.dot(pl.v_dir));
                    if self.uv_in_domain(f, uv) {
                        return true;
                    }
                }
            } else {
                let (uv, q) = crate::bop::closest_point_on_surface(&f.surf, p);
                if (p - q).length() <= tol && self.uv_in_domain(f, uv) {
                    return true;
                }
            }
        }
        false
    }

    /// True when the ray `p + s*dir` (s > 0) is degenerate: it lies in the
    /// supporting plane of a planar face of the solid (OCCT
    /// BRepClass3d_SClassifier::Perform L405-449 detects the parallel-line
    /// case per face; the faulty-line loop L264-285 retries with a different
    /// direction). The K1 disk point at z=100 with the +X ray: the ray lies
    /// in the annulus plane and only grazes the cylinder-wall rims, which the
    /// UV domain check rejects — the crossing count comes out zero even for
    /// an interior point, so the direction must be retried.
    pub fn ray_is_degenerate(&self, p: DVec3, dir: DVec3, tol: f64) -> bool {
        for f in &self.face_surfaces {
            if let Surface3::Plane(pl) = &f.surf {
                if dir.dot(pl.normal).abs() < 1e-12 && (p - pl.origin).dot(pl.normal).abs() <= tol {
                    return true;
                }
            }
        }
        false
    }

    /// Get face indices for classification.
    pub fn get_face_indices(&self) -> &[usize] {
        &self.face_indices
    }

    /// True when the explorer has face geometry to classify against.
    pub fn has_faces(&self) -> bool {
        !self.face_surfaces.is_empty() || self.ds.is_some()
    }

    /// Add a face index (used by IntToolsContext when building from DS).
    pub fn add_face_index(&mut self, fi: usize) {
        self.face_indices.push(fi);
    }

    /// Classify point using ray casting (even-odd rule: a point is IN when a
    /// ray from it crosses the solid boundary an odd number of times). Only
    /// intersections inside the face's UV domain are counted.
    /// OCCT BRepClass3d_SClassifier::Perform uses IntCurvesFace_Intersector
    /// (Line x Face) for every face — including curved ones; rcad: planar
    /// faces use the analytic intersection, curved faces use BeanFaceIntersector
    /// (the IntCurvesFace_Intersector translation).
    ///
    /// A fixed ray direction is degenerate when it lies in a face's plane
    /// (e.g. a point on the plane of a planar face, whose ray grazes the
    /// neighboring curved faces' rims): the rim hits are rejected by the UV
    /// domain check and the count comes out ZERO even for an inside point.
    /// OCCT retries with a different line direction (BRepClass3d_SClassifier.cxx
    /// L264-285 — the faulty-line loop over SolidExplorer::Segment/
    /// OtherSegment). rcad: first the ON-state check (a point within the
    /// tolerance of a face's surface with the UV domain check is ON — OCCT
    /// L405-449), then the +X ray, then the remaining axes; a degenerate
    /// direction (0 crossings, ray in a face's plane) is skipped, the first
    /// non-degenerate crossing count decides IN/OUT.
    pub fn classify_point(&self, p: DVec3) -> u8 {
        self.classify_point_with_tol(p, CONFUSION)
    }

    /// classify_point with an explicit tolerance for the degenerate-ray check.
    pub(crate) fn classify_point_with_tol(&self, p: DVec3, tol: f64) -> u8 {
        // The ON-state face-surface check is done by the caller
        // (SolidClassifier::perform → point_on_face, OCCT aSelectorPoint
        // L232-237 + the parallel-line face distance L405-449) before this
        // method runs; classify_point itself only ray-casts.
        // The +X ray first (the plain single-direction behavior for the
        // non-degenerate cases), then the remaining axes of the retry.
        let dirs = [
            DVec3::X,
            DVec3::Y,
            DVec3::Z,
            -DVec3::X,
            -DVec3::Y,
            -DVec3::Z,
        ];
        let mut retried = 0usize;
        for ray_dir in dirs {
            let intersections = self.count_ray_crossings(p, ray_dir);
            if intersections > 0 {
                return if intersections % 2 == 1 { 3 } else { 4 }; // IN=3, OUT=4
            }
            // 0 crossings: the direction is unreliable when the ray lies in a
            // face's plane (the K1 disk point's +X ray in the z=100 annulus
            // plane grazes the cylinder rims); otherwise a genuinely outside
            // point is OUT.
            if !self.ray_is_degenerate(p, ray_dir, tol) {
                return 4;
            }
            retried += 1;
        }
        4 // all directions degenerate — OUT (safe default)
    }

    /// Number of ray crossings of the solid boundary along `p + s*dir`,
    /// counting only intersections inside the face's UV domain.
    fn count_ray_crossings(&self, p: DVec3, ray_dir: DVec3) -> usize {
        let mut intersections = 0usize;
        if !self.face_surfaces.is_empty() {
            for f in &self.face_surfaces {
                let t = match &f.surf {
                    Surface3::Plane(pl) => {
                        let normal = if f.ori == Orientation::Reversed {
                            -pl.normal
                        } else {
                            pl.normal
                        };
                        let denom = ray_dir.dot(normal);
                        if denom.abs() < 1e-12 {
                            continue;
                        }
                        let t = (pl.origin - p).dot(normal) / denom;
                        if t > 1e-7 && in_face_uv(f.uv_bounds, pl, p + ray_dir * t) {
                            Some(t)
                        } else {
                            None
                        }
                    }
                    _ => self.ray_face_param(p, ray_dir, f),
                };
                if let Some(t) = t {
                    if t > 1e-7 {
                        intersections += 1;
                    }
                }
            }
        } else if let Some(ref ds) = self.ds {
            for &fi in &self.face_indices {
                let surf = match ds.face_surface(fi) {
                    Some(s) => s,
                    None => continue,
                };
                let face_ori = ds.shape_at(fi).orientation;
                if let Surface3::Plane(pl) = surf {
                    let normal = if face_ori == Orientation::Reversed {
                        -pl.normal
                    } else {
                        pl.normal
                    };
                    let denom = ray_dir.dot(normal);
                    if denom.abs() < 1e-12 {
                        continue;
                    }
                    let t = (pl.origin - p).dot(normal) / denom;
                    if t > 1e-7 {
                        intersections += 1;
                    }
                } else {
                    let ef = ExplorerFace {
                        key: (0, 0),
                        src: None,
                        surf: surf.clone(),
                        ori: face_ori,
                        uv_bounds: None,
                        uv_polys: None,
                        boundary: vec![],
                        edge_keys: vec![],
                        fclass_cache: std::cell::RefCell::new(HashMap::new()),
                    };
                    if let Some(t) = self.ray_face_param(p, ray_dir, &ef) {
                        if t > 1e-7 {
                            intersections += 1;
                        }
                    }
                }
            }
        }
        intersections
    }

    /// Parameter of the closest intersection of the ray `p + t*dir` with a
    /// curved face (OCCT BRepClass3d_SClassifier uses IntCurvesFace_Intersector
    /// for the Line x Face intersection; BeanFaceIntersector is its translation).
    /// Only intersections whose UV lies inside the face's 2D domain are counted
    /// (OCCT IntCurvesFace_Intersector.cxx L256-286: Classify(Puv) == IN/ON).
    pub(crate) fn ray_face_param(&self, p: DVec3, dir: DVec3, f: &ExplorerFace) -> Option<f64> {
        let line = Curve3::Line(rcad_kernel::geom::Line3 {
            origin: p,
            direction: dir,
        });
        let line_for_points = line.clone();
        let adapt_curve = BRepAdaptorCurve::new(line);
        let adapt_surf = BRepAdaptorSurface::new(f.surf.clone());
        let mut bfi = BeanFaceIntersector::with_adaptors(
            adapt_curve.clone(),
            adapt_surf.clone(),
            CONFUSION,
            CONFUSION,
        );
        bfi.set_bean_parameters(-f64::MAX, f64::MAX);
        bfi.set_surface_parameters(
            adapt_surf.first_u_parameter(),
            adapt_surf.last_u_parameter(),
            adapt_surf.first_v_parameter(),
            adapt_surf.last_v_parameter(),
        );
        bfi.perform();
        if !bfi.is_done() {
            return None;
        }
        // OCCT BRepClass3d_SClassifier: the intersection closest to the point.
        // A candidate is accepted only when its UV is IN the face's 2D domain
        // (IntCurvesFace_Intersector's currentstate IN/ON filter —
        // TopolTool->Classify).
        let mut parmin = f64::MAX;
        for r in bfi.result() {
            let t = r.first();
            let q = line_for_points.point_at(t);
            let (uv, _q2) = crate::bop::closest_point_on_surface(&f.surf, q);
            if !self.uv_in_domain(f, uv) {
                continue;
            }
            if t < parmin {
                parmin = t;
            }
        }
        if parmin == f64::MAX {
            None
        } else {
            Some(parmin)
        }
    }

    /// OCCT BRepClass3d_SolidExplorer::FindAPointInTheFace — a point inside the
    /// face, sampled by the probing parameter aParam in [0.1, 0.9] (the random
    /// inner point of the OCCT PerformInfinitePoint). Returns (point, u, v).
    /// OCCT (SolidExplorer.cxx L74-167): take a point on a boundary edge at the
    /// probe parameter, move into the face interior along the edge normal, and
    /// verify the result with the face classifier (FClass2d == IN). rcad
    /// replicates the same walk using the boundary pcurves and the UV-domain
    /// check (uv_in_face_domain).
    pub(crate) fn face_point(&self, f: &ExplorerFace, param: f64) -> Option<(DVec3, f64, f64)> {
        if let Surface3::Plane(pl) = &f.surf {
            // Planar faces: the UV-box interpolation point is inside the face
            // when the box is the face's bounding box (the vertex projection
            // box; the probe parameter samples [0.1, 0.9] of it). Verify with
            // the domain check and fall back to the edge-walk otherwise.
            if let Some([umin, umax, vmin, vmax]) = f.uv_bounds {
                let uv = DVec2::new(umin + (umax - umin) * param, vmin + (vmax - vmin) * param);
                if self.uv_in_domain(f, uv) {
                    return Some((f.surf.point_at(uv.x, uv.y), uv.x, uv.y));
                }
            }
        }
        // OCCT L93-167: edge walk — a point on a boundary edge moved into the
        // face interior along the rotated edge tangent.
        for (pc, range, ori) in &f.boundary {
            if !range[0].is_finite() || !range[1].is_finite() {
                continue;
            }
            let t = range[0] + (range[1] - range[0]) * param;
            let p2 = pc.point_at(t);
            let tan = pc.tangent_at(t);
            if tan.length_squared() < 1e-24 {
                continue;
            }
            // OCCT L108-112: FORWARD edge -> T=(-y,x); REVERSED edge -> T=(y,-x)
            // (the direction pointing into the face interior).
            let tan = if *ori == Orientation::Forward {
                DVec2::new(-tan.y, tan.x)
            } else {
                DVec2::new(tan.y, -tan.x)
            };
            // OCCT L113-121: move TolInit (0.00001) into the interior and find
            // the nearest boundary intersection of the ray P + s*T.
            let p_start = p2 + 1e-5 * tan;
            let Some(mut param_init) = ray_boundary_hit(f, p_start, tan) else {
                continue;
            };
            // OCCT L122-158: walk inward (x0.41234 each step) until the point
            // is classified IN by the face classifier.
            let mut guard = 0;
            loop {
                param_init *= 0.41234;
                if param_init < 1e-7 {
                    // OCCT L157: ParamInit < Precision::PConfusion() -> false
                    return None;
                }
                let uv = p_start + param_init * tan;
                // OCCT L130-135: the point must be strictly IN the face.
                if !self.uv_in_domain(f, uv) {
                    return None;
                }
                // OCCT L138-148: non-degenerate surface point.
                let n = f.surf.normal_at(uv.x, uv.y);
                if n.length_squared() > 1e-24 {
                    return Some((f.surf.point_at(uv.x, uv.y), uv.x, uv.y));
                }
                guard += 1;
                if guard > 64 {
                    return None;
                }
            }
        }
        None
    }

    /// OCCT BRepClass3d_SClassifier::FaceNormal (SClassifier.cxx L606-627) —
    /// the surface normal at (u, v), flipped for a REVERSED face.
    pub(crate) fn face_outward_normal_at(f: &ExplorerFace, u: f64, v: f64) -> Option<DVec3> {
        let mut n = f.surf.normal_at(u, v);
        if n.length_squared() < 1e-24 {
            return None;
        }
        if f.ori == Orientation::Reversed {
            n = -n;
        }
        Some(n)
    }

    /// Number of collected faces (OCCT TopExp_Explorer on FACE).
    pub(crate) fn nb_faces(&self) -> usize {
        self.face_surfaces.len()
    }

    /// Key of the vertex at index (BndBoxTreeSelectorLine vertex param
    /// index) — the vertex identity for LVInts. rcad vertices are points
    /// without shape keys; a synthetic key per index is used.
    pub(crate) fn vertex_key(&self, idx: usize) -> (u64, u32) {
        (u64::MAX - idx as u64, 0)
    }

    /// Key of the edge at index — the edge identity for mapEF lookups.
    pub(crate) fn edge_key(&self, idx: usize) -> (u64, u32) {
        // map_ef keys are (ptr_id, location) of the source edges; the
        // selector's edge index must map back. rcad keeps the edge curves in
        // collection order; the keys are recovered by matching the curve.
        (0, idx as u32)
    }

    /// The vertex keys of an edge (for the LVInts skip in GetTransi).
    pub(crate) fn edge_vertices(&self, _idx: usize) -> ((u64, u32), (u64, u32)) {
        ((0, 0), (0, 0))
    }

    /// OCCT GetNormalOnFaceBound (SClassifier.cxx L631-651): the face normal
    /// at the 2D point of the boundary edge at `param`. rcad approximates it
    /// with the face's surface normal at the face-center parameter.
    pub(crate) fn face_bound_normal(&self, fi: usize, _param: f64) -> Option<DVec3> {
        let f = self.face_surfaces.get(fi)?;
        let (u1, v1, u2, v2) = face_uv_bounds(f);
        let u = if u1.is_finite() && u2.is_finite() { (u1 + u2) * 0.5 } else { 0.0 };
        let v = if v1.is_finite() && v2.is_finite() { (v1 + v2) * 0.5 } else { 0.0 };
        Self::face_outward_normal_at(f, u, v)
    }

    /// OCCT IntCurvesFace_Intersector: the distance from P to the face's
    /// surface and its UV — used by the parallel-line ON check (L380-404).
    pub(crate) fn point_face_distance(&self, p: DVec3, fi: usize) -> Option<(f64, DVec2)> {
        let f = self.face_surfaces.get(fi)?;
        let (uv, q) = closest_point_on_surface(&f.surf, p);
        Some(((p - q).length_squared(), uv))
    }

    /// OCCT IntCurvesFace_Intersector::Bounding() (Intersector.cxx L529-541):
    /// the polyhedron bounding box of the face — empty for analytic surfaces
    /// (Plane/Cylinder/Cone/Sphere/Torus have no polyhedron, L141-145). The
    /// polyhedron is a grid of (NbU+1) x (NbV+1) surface samples over the
    /// face's UV domain (ThePolyhedronOfHInter); rcad samples 11 x 11.
    pub(crate) fn face_bounding_box(&self, fi: usize) -> Option<(DVec3, DVec3)> {
        let f = self.face_surfaces.get(fi)?;
        if Self::is_analytic_surface(&f.surf) {
            return None;
        }
        let (u1, v1, u2, v2) = face_uv_bounds(f);
        if !(u1.is_finite() && v1.is_finite() && u2.is_finite() && v2.is_finite()) {
            return None;
        }
        if (u2 - u1).abs() < CONFUSION || (v2 - v1).abs() < CONFUSION {
            return None;
        }
        let mut min = DVec3::splat(f64::INFINITY);
        let mut max = DVec3::splat(f64::NEG_INFINITY);
        let n = 10usize;
        for i in 0..=n {
            let u = u1 + (u2 - u1) * (i as f64) / (n as f64);
            for j in 0..=n {
                let v = v1 + (v2 - v1) * (j as f64) / (n as f64);
                let p = f.surf.point_at(u, v);
                min = min.min(p);
                max = max.max(p);
            }
        }
        Some((min, max))
    }

    /// OCCT GetAddToParam (SClassifier.cxx L569-602): the largest line
    /// parameter over the 8 corners of the face's bounding box, minus Par —
    /// the extension of the intersection range needed to reach the face's far
    /// side (ElCLib::Parameter = (P - Location) . Direction).
    pub(crate) fn get_add_to_param(
        &self,
        l_origin: DVec3,
        l_dir: DVec3,
        par: f64,
        bb: (DVec3, DVec3),
    ) -> f64 {
        let (amin, amax) = bb;
        let xs = [amin.x, amax.x];
        let ys = [amin.y, amax.y];
        let zs = [amin.z, amax.z];
        let dir = l_dir.normalize_or_zero();
        let mut out_par = par;
        for &x in &xs {
            for &y in &ys {
                for &z in &zs {
                    let dx = (x - l_origin.x).abs();
                    let dy = (y - l_origin.y).abs();
                    let dz = (z - l_origin.z).abs();
                    if dx < 1e20 && dy < 1e20 && dz < 1e20 {
                        let t = (DVec3::new(x, y, z) - l_origin).dot(dir);
                        if t > out_par {
                            out_par = t;
                        }
                    } else {
                        return 1e20;
                    }
                }
            }
        }
        out_par - par
    }

    /// OCCT IntCurvesFace_Intersector constructor L141-145: analytic surfaces
    /// (Plane/Cylinder/Cone/Sphere/Torus) skip the polyhedron entirely.
    fn is_analytic_surface(s: &Surface3) -> bool {
        matches!(
            s,
            Surface3::Plane(_)
                | Surface3::Cylinder(_)
                | Surface3::Cone(_)
                | Surface3::Sphere(_)
                | Surface3::Torus(_)
        )
    }

    /// OCCT IntCurveSurface_HInter::IsParallel (HInter.hxx L109-111): the
    /// curve is parallel to or belongs to the surface — "recognized only for
    /// some pairs of analytical curves and surfaces (plane - line, ...)".
    /// Line-plane: |dir . normal| ~ 0; line-cylinder/cone: the line is
    /// parallel to the axis (the direct algorithm returns no point).
    fn line_face_is_parallel(s: &Surface3, l_dir: DVec3) -> bool {
        let d = l_dir.normalize_or_zero();
        match s {
            Surface3::Plane(p) => d.dot(p.normal).abs() < 1e-12,
            Surface3::Cylinder(c) => d.dot(c.axis).abs() > 1.0 - 1e-12,
            Surface3::Cone(c) => d.dot(c.axis).abs() > 1.0 - 1e-12,
            _ => false,
        }
    }
    /// OCCT IntCurvesFace_Intersector::Perform(L, minW, maxW) — the line-face
    /// intersection points within [minW, maxW]. Each point carries (w, state,
    /// transition, u, v) with state 1=IN, 2=ON, 3=OUT and transition 0=Tangent,
    /// 1=In, 2=Out. rcad: the analytic curve-surface intersection (IntCS)
    /// filtered by the face's 2D domain (TopolTool::Classify).
    pub(crate) fn face_line_intersections(
        &self,
        fi: usize,
        l_origin: DVec3,
        l_dir: DVec3,
        min_w: f64,
        max_w: f64,
    ) -> Option<(bool, Vec<(f64, u8, u8, f64, f64)>)> {
        let f = self.face_surfaces.get(fi)?;
        // OCCT IntCurvesFace_Intersector: the Plane/Cylinder/Cone/Sphere/Torus
        // surfaces use the ANALYTIC IntCurveSurface_Intersector
        // (IntAna_IntConicQuad) — the exact line-quadric quadratic roots —
        // while other surfaces use the polyhedron sampling.  rcad's sampling
        // IntCS misses roots when the line's inside-span is far smaller than
        // the sampling step (e.g. a chord through a cylinder), so the analytic
        // path is used for the quadric surfaces.
        let quad = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(&f.surf);
        if let Some(quad) = quad {
            if quad.type_quadric() != crate::geomalgo::int_surf::quadric::QuadricType::Other {
                let line = rcad_kernel::geom::Line3 {
                    origin: l_origin,
                    direction: l_dir,
                };
                let (in_quadric, pts) = match crate::geomalgo::int_patch::int_cs::intersect_line_quadric(
                    &line, &quad,
                ) {
                    None => return None,
                    Some(r) => r,
                };
                if in_quadric {
                    // The line lies in the surface — parallel (the
                    // Extrema_ExtPS branch in the SClassifier applies).
                    return Some((true, Vec::new()));
                }
                let mut out = Vec::new();
                for (p3d, w) in pts {
                    if w < min_w || w > max_w {
                        continue;
                    }
                    let (uv, _) = crate::bop::closest_point_on_surface(&f.surf, p3d);
                    let state = self.classify_uv_point(f, uv);
                    if state == 3 {
                        continue; // OUT — not an intersection of the face.
                    }
                    let (_pnt, d1u, d1v) = f.surf.derivatives(uv.x, uv.y);
                    let n = d1u.cross(d1v);
                    let norm = n.length();
                    let d1 = l_dir;
                    let d1_mag = d1.length();
                    let mut tran = if norm > 1e-12 && d1_mag > 1e-12 {
                        let cos_dir = n.dot(d1) / (norm * d1_mag);
                        if -cos_dir > 1e-12 {
                            1 // In
                        } else if cos_dir > 1e-12 {
                            2 // Out
                        } else {
                            0 // Tangent
                        }
                    } else {
                        0 // Tangent
                    };
                    if tran != 0 && f.ori == Orientation::Reversed {
                        tran = if tran == 1 { 2 } else { 1 };
                    }
                    out.push((w, state, tran, uv.x, uv.y));
                }
                return Some((false, out));
            }
        }
        let line = Curve3::Line(rcad_kernel::geom::Line3 {
            origin: l_origin,
            direction: l_dir,
        });
        let mut hics = rcad_kernel::base::geom_api::int_cs::IntCS::new();
        hics.perform(&line, &f.surf);
        if !hics.is_done() {
            return None;
        }
        // OCCT IntCurveSurface_HInter::IsParallel — the line is parallel to or
        // lies in the face's surface (recognized for the analytical pairs
        // plane-line, cylinder/cone-line, ...). Parallel means NbPnt()==0 and
        // the SClassifier's Extrema_ExtPS branch (L400-443) applies.
        let is_parallel = Self::line_face_is_parallel(&f.surf, l_dir);
        let mut out = Vec::new();
        for idx in 1..=hics.nb_points() {
            let pt = hics.point(idx);
            if pt.w < min_w || pt.w > max_w {
                continue;
            }
            let uv = DVec2::new(pt.u, pt.v);
            // OCCT InternalCall L256-257: state = TopolTool->Classify(Puv, 0)
            // — the intersector is built with UseBToler=false (InitShape L924),
            // so the 2D domain test is STRICT (tolerance 0, no 3D selector).
            let state = self.classify_uv_point(f, uv);
            if state == 3 {
                continue; // OUT — not an intersection of the face.
            }
            // OCCT ComputeTransitions: cos_dir = N . dC/dw.
            let (_pnt, d1u, d1v) = f.surf.derivatives(pt.u, pt.v);
            let n = d1u.cross(d1v);
            let norm = n.length();
            let d1 = l_dir;
            let d1_mag = d1.length();
            let mut tran = if norm > 1e-12 && d1_mag > 1e-12 {
                let cos_dir = n.dot(d1) / (norm * d1_mag);
                if -cos_dir > 1e-12 {
                    1 // In
                } else if cos_dir > 1e-12 {
                    2 // Out
                } else {
                    0 // Tangent
                }
            } else {
                0 // Tangent
            };
            // OCCT IntCurvesFace_Intersector::InternalCall L262-269: the
            // transition is flipped for a REVERSED face (the surface normal
            // points inward).
            if tran != 0 && f.ori == Orientation::Reversed {
                tran = if tran == 1 { 2 } else { 1 };
            }
            out.push((pt.w, state, tran, pt.u, pt.v));
        }
        Some((is_parallel, out))
    }

    /// OCCT BRepClass3d_SClassifier L214: aMapEV — the vertices and edges of
    /// the solid (BndBoxTreeSelectorPoint/Line data). Returns the edge curves,
    /// tolerances (BRep_Tool::Tolerance) and ranges (BRep_Tool::Range), and
    /// the vertex points with tolerances.
    pub(crate) fn map_ev(&self) -> (Vec<Curve3>, Vec<f64>, Vec<[f64; 2]>, Vec<DVec3>, Vec<f64>) {
        let mut edges = Vec::new();
        let mut e_tols = Vec::new();
        let mut e_ranges = Vec::new();
        for (i, e) in self.edges.iter().enumerate() {
            if let Some(curve) = e {
                edges.push(curve.clone());
                e_tols.push(self.edge_tols.get(i).copied().unwrap_or(0.0));
                e_ranges.push(self.edge_ranges.get(i).copied().unwrap_or([f64::NEG_INFINITY, f64::INFINITY]));
            }
        }
        let mut verts = Vec::new();
        let mut v_tols = Vec::new();
        for (i, v) in self.vertices.iter().enumerate() {
            verts.push(*v);
            v_tols.push(self.vert_tols.get(i).copied().unwrap_or(0.0));
        }
        (edges, e_tols, e_ranges, verts, v_tols)
    }

    /// OCCT BRepClass3d_SClassifier L232: mapEF — the edge -> adjacent faces
    /// map of the solid (TopExp::MapShapesAndAncestors EDGE, FACE). rcad:
    /// built from the collected faces' boundary edges; the key is the edge
    /// index in the explorer's edge list (matching edge_key).
    pub(crate) fn map_ef(&self) -> HashMap<(u64, u32), Vec<usize>> {
        let mut m: HashMap<(u64, u32), Vec<usize>> = HashMap::new();
        // The edge list is built in init_shape in tree order; the face edge
        // keys are (ptr_id, location). Map each collected edge curve back to
        // its index by identity (curve pointer is not retained, so match by
        // the key via the source shapes). rcad: the explorer keeps the edges
        // as curves only; the edge -> face adjacency is approximated by
        // matching the edge index to the faces whose boundary contains it.
        for (fi, f) in self.face_surfaces.iter().enumerate() {
            for ekey in &f.edge_keys {
                // Edge keys are (ptr_id, location) of the source edges — not
                // the explorer's edge index. The explorer's edge list has no
                // keys; use the source shape to find the index.
                if let Some(ei) = self.edge_index_of(ekey) {
                    m.entry((0, ei as u32)).or_default().push(fi);
                }
            }
        }
        m
    }

    /// Index of the edge with the given (ptr_id, location) in the explorer's
    /// edge list. The list is collected in tree order from the shape; the
    /// source key is recovered by re-walking the shape (the curve list does
    /// not retain keys). rcad keeps the edge curves in collection order and
    /// the source shape tree in `self.shape`; the mapping is built once here.
    fn edge_index_of(&self, key: &(u64, u32)) -> Option<usize> {
        let shape = self.shape.as_ref()?;
        let mut idx = 0usize;
        let mut stack: Vec<Shape> = vec![shape.clone()];
        while let Some(sh) = stack.pop() {
            match &*sh.data {
                TShape::Solid(sd) => {
                    for x in &sd.shells {
                        stack.push(x.clone());
                    }
                }
                TShape::CompSolid(cd) => {
                    for x in cd {
                        stack.push(x.clone());
                    }
                }
                TShape::Compound(cd) => {
                    for x in cd {
                        stack.push(x.clone());
                    }
                }
                TShape::Shell(sd) => {
                    for x in &sd.faces {
                        stack.push(x.clone());
                    }
                }
                TShape::Face(fd) => {
                    stack.push(fd.outer_wire.clone());
                    for w in &fd.inner_wires {
                        stack.push(w.clone());
                    }
                }
                TShape::Wire(wd) => {
                    for e in &wd.edges {
                        stack.push(e.clone());
                    }
                }
                TShape::Edge(_ed) => {
                    if (sh.ptr_id(), sh.location) == *key {
                        return Some(idx);
                    }
                    idx += 1;
                }
                _ => {}
            }
        }
        None
    }

    /// OCCT BRepClass3d_SClassifier::GetFaceSegmentIndex — myFirstFace.
    pub(crate) fn face_segment_index(&self) -> usize {
        self.my_first_face
    }

    /// OCCT SolidExplorer::Segment/OtherSegment — builds a line from P toward
    /// a point found on a face of the solid. Returns (iFlag, line origin,
    /// line direction, Par) with iFlag: 0 = OK, 1 = point on an infinite face
    /// (ON), 2 = degenerate face (OUT), 3 = point on surface but outside face.
    pub(crate) fn segment(&mut self, p: DVec3, is_other: bool) -> (i32, DVec3, DVec3, f64) {
        if !is_other {
            self.my_first_face = 0;
        }
        self.other_segment(p)
    }

    /// OCCT SolidExplorer::OtherSegment (SolidExplorer.cxx L493-622) —
    /// full translation (see temp/sclassifier_alignment.md).
    fn other_segment(&mut self, p: DVec3) -> (i32, DVec3, DVec3, f64) {
        let tol_u = CONFUSION;
        let tol_v = tol_u;
        loop {
            self.my_first_face += 1;
            let n_faces = self.face_surfaces.len();
            let mut ptfound = false;
            let mut maxscal = 0.0f64;
            let mut l: (DVec3, DVec3, f64) = (p, DVec3::X, 1.0);
            let mut index_point = 0usize;
            let mut nb_points_ok = 0usize;
            // OCCT L514-527: aTestInvert — on the retry loop (myFirstFace
            // reset), each face is re-checked with FClass2d::PerformInfinitePoint
            // to decide aRestr (restricted vs natural UV domain). rcad has no
            // BRepAdaptor_Surface cache; face_uv_bounds already returns the
            // natural (infinite) domain for unbounded surfaces and the stored
            // bounds for restricted ones — the aRestr switch is equivalent.
            let mut _a_test_invert = false;
            for (fi, f) in self.face_surfaces.iter().enumerate() {
                if self.my_first_face > n_faces {
                    break;
                }
                if self.my_first_face > fi + 1 {
                    continue;
                }
                let sv_myparam = self.my_param_on_edge;
                let (u1, v1, u2, v2) = face_uv_bounds(f);
                // OCCT L534-541: degenerate face (|U2-U1| or |V2-V1| < eps)
                // -> return 2 (OUT).
                let eps = CONFUSION;
                let eps_u = (eps * u2.abs().max(u1.abs())).max(eps);
                let eps_v = (eps * v2.abs().max(v1.abs())).max(eps);
                if (u2 - u1).abs() < eps_u || (v2 - v1).abs() < eps_v {
                    return (2, p, DVec3::X, 0.0);
                }
                let an_inf_flag = is_infinite_uv(u1, v1, u2, v2);
                let mut _u = (u1 + u2) * 0.5;
                let mut _v = (v1 + v2) * 0.5;
                let mut a_point = f.surf.point_at(_u, _v);
                // OCCT L566-595: Extrema_ExtPS(P, surface) — nearest point.
                let (pu, pv) = closest_point_on_surface(&f.surf, p);
                let proj = f.surf.point_at(pu.x, pu.y);
                let dist2 = (proj - p).length_squared();
                if dist2 < 1e-24 {
                    // P is on the surface.
                    if an_inf_flag != 0 {
                        return (1, p, DVec3::X, 0.0); // ON (infinite face)
                    } else {
                        // OCCT L586-592: BRepClass_FaceClassifier::Perform
                        // (face, aPuv, Precision::PConfusion()) — the 2D
                        // point-in-face state; IN or ON returns 1 (the point is
                        // on the face), OUT returns 3 (on the surface but
                        // outside the face domain). No 3D edge/vertex selector.
                        let st = self.classify_uv_2d(f, DVec2::new(pu.x, pu.y), CONFUSION);
                        if st == 1 || st == 2 {
                            return (1, p, DVec3::X, 0.0); // ON
                        } else {
                            return (3, p, DVec3::X, 0.0); // on surface, outside face
                        }
                    }
                }
                if an_inf_flag != 0 {
                    // OCCT L597-603: infinite face — build the line to the
                    // projection point.
                    let v = proj - p;
                    let par = v.length();
                    let dir = if par > 1e-300 { v / par } else { DVec3::X };
                    return (0, p, dir, par);
                }
                _u = pu.x;
                _v = pu.y;
                a_point = proj;
                // OCCT L605-621: PointInTheFace probing.
                let mut guard = 0;
                loop {
                    index_point += 1;
                    let found = self.point_in_the_face(f, &mut _u, &mut _v, &mut a_point, index_point, u1, v1, u2, v2);
                    if found {
                        nb_points_ok += 1;
                        let v = a_point - p;
                        let par = v.length();
                        let (_pnt, d1u, d1v) = f.surf.derivatives(_u, _v);
                        if par > 1e-12 && d1u.length_squared() > 1e-12 && d1v.length_squared() > 1e-12 {
                            let norm = d1u.cross(d1v);
                            let tt = norm.length();
                            if tt > 1e-12 {
                                let tt = (norm.dot(v)).abs() / (tt * par);
                                if tt > maxscal {
                                    maxscal = tt;
                                    let dir = if par > 1e-300 { v / par } else { DVec3::X };
                                    l = (p, dir, par);
                                    ptfound = true;
                                    if maxscal > 0.2 {
                                        self.my_param_on_edge = sv_myparam;
                                        return (0, p, dir, par);
                                    }
                                }
                            }
                        }
                    }
                    guard += 1;
                    if guard > 200 || nb_points_ok >= 16 {
                        break;
                    }
                }
                self.my_param_on_edge = sv_myparam;
                if maxscal > 0.2 {
                    return (0, l.0, l.1, l.2);
                }
                index_point = 0;
                let encore = fi + 1 < n_faces;
                if !ptfound && !encore && self.my_param_on_edge < 0.0001 {
                    // OCCT L647-653: point on a solid reduced to a face.
                    let dir = DVec3::X;
                    return (0, p, dir, 1.0);
                }
            }
            if n_faces == 0 {
                self.my_reject = true;
                return (0, p, DVec3::X, 0.0);
            }
            if ptfound {
                return (0, l.0, l.1, l.2);
            }
            self.my_first_face = 0;
            // OCCT L675-700: cycle myParamOnEdge.
            self.my_param_on_edge = next_param_on_edge(self.my_param_on_edge);
            _a_test_invert = true;
        }
    }

    /// OCCT SolidExplorer::PointInTheFace (L191-420, L804-870) — find a point
    /// inside the face by scanning the UV grid, starting at IndexPoint.
    fn point_in_the_face(
        &self,
        f: &ExplorerFace,
        u_: &mut f64,
        v_: &mut f64,
        a_point: &mut DVec3,
        index_point: usize,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> bool {
        let mut du = (u2 - u1) / 6.0;
        let mut dv = (v2 - v1) / 6.0;
        if du < 1e-12 {
            du = 1e-12;
        }
        if dv < 1e-12 {
            dv = 1e-12;
        }
        let is_not_u_per = !f.surf.is_u_periodic();
        let is_not_v_per = !f.surf.is_v_periodic();
        let mut nb_pnt_calc = 0usize;
        // OCCT L839-849: the current point, if inside and classified IN.
        let mut u = *u_;
        let mut v = *v_;
        let is_inside = (!is_not_u_per || (u >= u1 && u <= u2)) && (!is_not_v_per || (v >= v1 && v <= v2));
        if is_inside && self.classify_uv_point(f, DVec2::new(u, v)) == 1 {
            let pnt = f.surf.point_at(u, v);
            if pnt.distance_squared(*a_point) < CONFUSION * CONFUSION {
                return true;
            }
        }
        // OCCT L850-886: 4 quarter scans + remainder + center.
        let mid_u = (u1 + u2) * 0.5;
        let mid_v = (v1 + v2) * 0.5;
        let quadrants: [(f64, f64, f64, f64, f64, f64); 4] = [
            (du + mid_u, u2, dv + mid_v, v2, 1.0, 1.0),
            (-du + mid_u, u1, -dv + mid_v, v1, -1.0, -1.0),
            (-du + mid_u, u1, dv + mid_v, v2, -1.0, 1.0),
            (du + mid_u, u2, -dv + mid_v, v1, 1.0, -1.0),
        ];
        for (us, ue, vs, ve, su, sv) in quadrants {
            let mut uu = us;
            while if su > 0.0 { uu < ue } else { uu > ue } {
                let mut vv = vs;
                while if sv > 0.0 { vv < ve } else { vv > ve } {
                    nb_pnt_calc += 1;
                    if nb_pnt_calc >= index_point && self.classify_uv_point(f, DVec2::new(uu, vv)) == 1 {
                        *u_ = uu;
                        *v_ = vv;
                        *a_point = f.surf.point_at(uu, vv);
                        return true;
                    }
                    vv += dv * sv;
                }
                uu += du * su;
            }
        }
        // OCCT L887-905: remainder grid (37 divisions).
        du = (u2 - u1) / 37.0;
        dv = (v2 - v1) / 37.0;
        if du < 1e-12 {
            du = 1e-12;
        }
        if dv < 1e-12 {
            dv = 1e-12;
        }
        let mut uu = du + u1;
        while uu < u2 {
            let mut vv = dv + v1;
            while vv < v2 {
                nb_pnt_calc += 1;
                if nb_pnt_calc >= index_point && self.classify_uv_point(f, DVec2::new(uu, vv)) == 1 {
                    *u_ = uu;
                    *v_ = vv;
                    *a_point = f.surf.point_at(uu, vv);
                    return true;
                }
                vv += dv;
            }
            uu += du;
        }
        // Center point (L907-912).
        let uu = (u1 + u2) * 0.5;
        let vv = (v1 + v2) * 0.5;
        nb_pnt_calc += 1;
        if nb_pnt_calc >= index_point && self.classify_uv_point(f, DVec2::new(uu, vv)) == 1 {
            *u_ = uu;
            *v_ = vv;
            *a_point = f.surf.point_at(uu, vv);
            return true;
        }
        // OCCT L869-870: no grid point found — FindAPointInTheFace fallback
        // (the edge-walk of SolidExplorer.cxx L74-190, rcad face_point).
        match self.face_point(f, self.my_param_on_edge) {
            Some((pt, pu, pv)) => {
                *u_ = pu;
                *v_ = pv;
                *a_point = pt;
                true
            }
            None => false,
        }
    }

    /// OCCT SolidExplorer::ClassifyUVPoint (SolidExplorer.cxx L221-237): a
    /// point on the surface near an edge/vertex of the solid is ON; otherwise
    /// the face-domain classification with the tolerance 1e-7. Used by
    /// PointInTheFace (L283-401) — the grid points must classify strictly IN.
    pub(crate) fn classify_uv_point(&self, f: &ExplorerFace, uv: DVec2) -> u8 {
        let p3d = f.surf.point_at(uv.x, uv.y);
        // OCCT L226-235: BndBoxTreeSelectorPoint over the vertex/edge map —
        // the point within the vertex/edge tolerance is ON.
        let (edges, e_tols, e_ranges, verts, v_tols) = self.map_ev();
        let mut sel = crate::topalgo::brep_class3d::bnd_box_tree::BndBoxTreeSelectorPoint::new(
            edges, e_tols, e_ranges, verts, v_tols,
        );
        sel.set_current_point(p3d);
        if sel.select() > 0 {
            return 2; // ON
        }
        // OCCT L236: theIntersector.ClassifyUVPoint(theP2d) =
        // myTopolTool->Classify(Puv, 1e-7) (Intersector.cxx L543-547).
        self.classify_uv_2d(f, uv, 1e-7)
    }

    /// OCCT BRepTopAdaptor_TopolTool::Classify(Puv, Tol) — the pure 2D
    /// face-domain classification (BRepClass_FaceClassifier on the face's 2D
    /// wire): IN when the point is inside the domain or within Tol of its
    /// boundary, OUT otherwise. The 3D edge/vertex selector is NOT part of it
    /// (that is SolidExplorer::ClassifyUVPoint only).
    ///
    /// rcad runs the IntTools_FClass2d translation (the OCCT pcurve-sampled
    /// 2D classifier): for a face registered in the DS under its own index,
    /// via the DS; for a synthetic draft-solid face (no DS registration) via
    /// the single-face ShapeSource adapter (OCCT IntTools_FClass2d::Init takes
    /// a TopoDS_Face — the adapter restores that contract). The sampled
    /// uv_polys domain remains the fallback for faces without pcurves.
    pub(crate) fn classify_uv_2d(&self, f: &ExplorerFace, uv: DVec2, tol: f64) -> u8 {
        use crate::topalgo::brep_top_adaptor::fclass2d::{FClass2d, State};
        if let Some(ds) = &self.ds {
            if let Some(fidx) = ds.map_shape_index(f.key.0, f.key.1) {
                // OCCT IntTools_Context::FClass2dMap — build the classifier
                // once per (face, tolerance) and reuse it across the probe
                // loop (IntTools_Context.hxx).
                let key = (fidx, tol.to_bits());
                if let Some(fc) = f.fclass_cache.borrow().get(&key) {
                    return match fc.perform(ds.as_ref(), uv, true) {
                        State::In => 1,  // IN
                        State::On => 2,  // ON
                        _ => 3, // OUT / UNKNOWN
                    };
                }
                let f2 = std::rc::Rc::new(FClass2d::new(ds.as_ref(), fidx, tol));
                let res = match f2.perform(ds.as_ref(), uv, true) {
                    State::In => 1,  // IN
                    State::On => 2,  // ON
                    _ => 3, // OUT / UNKNOWN
                };
                f.fclass_cache.borrow_mut().insert(key, f2);
                return res;
            }
        }
        if let Some(src) = &f.src {
            let adapter =
                crate::topalgo::shape_source::FaceShapeSource::new(src, f.surf.clone(), &self.locations);
            // Synthetic face with no DS registration: fixed face index 0 in
            // the cache key, one classifier per tolerance.
            let key = (0, tol.to_bits());
            if let Some(fc) = f.fclass_cache.borrow().get(&key) {
                return match fc.perform(&adapter, uv, true) {
                    State::In => 1,  // IN
                    State::On => 2,  // ON
                    _ => 3, // OUT / UNKNOWN
                };
            }
            let f2 = std::rc::Rc::new(FClass2d::new(&adapter, 0, tol));
            let res = match f2.perform(&adapter, uv, true) {
                State::In => 1,  // IN
                State::On => 2,  // ON
                _ => 3, // OUT / UNKNOWN
            };
            f.fclass_cache.borrow_mut().insert(key, f2);
            return res;
        }
        if uv_in_face_domain_with_tol(f, uv, tol) {
            1 // IN
        } else {
            3 // OUT
        }
    }

    /// TopolTool::Classify == IN/ON (OCCT BRepClass_FaceClassifier semantics;
    /// rcad: the IntTools_FClass2d translation, see classify_uv_2d).
    pub(crate) fn uv_in_domain(&self, f: &ExplorerFace, uv: DVec2) -> bool {
        match self.classify_uv_2d(f, uv, 1e-7) {
            1 | 2 => true,
            _ => false,
        }
    }

    /// OCCT IntCurvesFace_Intersector::ClassifyUVPoint (Intersector.cxx
    /// L543-547): myTopolTool->Classify(Puv, 1e-7) — used by the SClassifier's
    /// parallel-line ON check (SClassifier.cxx L431).
    pub(crate) fn classify_uv_point_at(&self, fi: usize, uv: DVec2) -> u8 {
        match self.face_surfaces.get(fi) {
            Some(f) => self.classify_uv_2d(f, uv, 1e-7),
            None => 3, // OUT
        }
    }

    /// Set the shape-source reference for face index lookups.
    pub fn set_ds<S: ShapeSource + 'static>(&mut self, ds: &Arc<S>) {
        self.ds = Some(ds.clone() as Arc<dyn ShapeSource>);
    }
}

/// Build the 2D domain of a face — the outer polygon (and holes) of the UV
/// points obtained by projecting the wire vertices onto the surface.
/// OCCT BRepTopAdaptor_TopolTool classifies the intersection UV against the
/// face's 2D domain; rcad rebuilds it from the wire vertices (same approach as
/// builder::point_in_face_image, IntTools_FClass2d::Init samples the pcurves).
fn build_uv_polys(face: &Shape, surf: &Surface3, locations: &[glam::DAffine3]) -> (Vec<DVec2>, Vec<Vec<DVec2>>) {
    // A wire whose edges are all closed (circle first == last) projects all
    // its vertices onto the same seam point (u = 0 for a full-revolution
    // cylindrical face) — the polygon degenerates and the ray-cast domain
    // check would reject every interior point.  Sample the edge 3D curves
    // instead and fall back to the u/v bounds rectangle (the band's natural
    // UV domain, OCCT BRepTopAdaptor_TopolTool::Classify on the pcurves).
    let poly_usable = |poly: &Vec<DVec2>| -> bool {
        if poly.len() < 3 {
            return false;
        }
        let u0 = poly[0].x;
        !poly.iter().all(|p| (p.x - u0).abs() < 1e-9)
    };
    let sample_bounds = |w: &Shape| -> Vec<DVec2> {
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;
        let mut found = false;
        if let TShape::Wire(wd) = &*w.data {
            let w_or = w.orientation;
            for e in &wd.edges {
                if let TShape::Edge(ed) = &*e.data {
                    let Some(curve) = &ed.curve else { continue };
                    let eff = w_or.compose(e.orientation);
                    let (t1, t2) = (ed.range[0], ed.range[1]);
                    let (ta, tb) = match eff {
                        Orientation::Reversed => (t2, t1),
                        _ => (t1, t2),
                    };
                    const N: usize = 16;
                    for i in 0..=N {
                        let t = ta + (tb - ta) * (i as f64) / (N as f64);
                        let mut p = curve.point_at(t);
                        if e.location != 0 {
                            if let Some(tr) = locations.get(e.location as usize) {
                                p = tr.transform_point3(p);
                            }
                        }
                        let (uv, _) = closest_point_on_surface(surf, p);
                        umin = umin.min(uv.x);
                        umax = umax.max(uv.x);
                        vmin = vmin.min(uv.y);
                        vmax = vmax.max(uv.y);
                        found = true;
                    }
                }
            }
        }
        if !found || !(umin < umax && vmin < vmax) {
            return Vec::new();
        }
        // Bring u into one period (the caps wrap).
        let tau = std::f64::consts::TAU;
        if umax - umin > tau * 0.9 {
            umin = 0.0;
            umax = tau;
        }
        vec![
            glam::DVec2::new(umin, vmin),
            glam::DVec2::new(umax, vmin),
            glam::DVec2::new(umax, vmax),
            glam::DVec2::new(umin, vmax),
        ]
    };
    let build_poly = |w: &Shape| -> Vec<DVec2> {
        let mut poly: Vec<DVec2> = Vec::new();
        if let TShape::Wire(wd) = &*w.data {
            // OCCT: the 2D wire of the face (FClass2d) walks each edge with
            // the composed orientation (TopExp_Explorer cumOri = wire
            // orientation x edge orientation); a REVERSED edge contributes its
            // first vertex as the contour endpoint, a FORWARD one its last.
            let w_or = w.orientation;
            for e in &wd.edges {
                if let TShape::Edge(ed) = &*e.data {
                    let e_or = w_or.compose(e.orientation);
                    let v = if e_or == Orientation::Forward { &ed.last } else { &ed.first };
                    if let TShape::Vertex(vd) = &*v.data {
                        let mut p = vd.point;
                        if v.location != 0 {
                            if let Some(tr) = locations.get(v.location as usize) {
                                p = tr.transform_point3(p);
                            }
                        }
                        let (uv, _) = closest_point_on_surface(surf, p);
                        poly.push(uv);
                    }
                }
            }
        }
        if !poly_usable(&poly) {
            poly = sample_bounds(w);
        }
        poly
    };
    let mut holes: Vec<Vec<DVec2>> = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            for w in &fd.inner_wires {
                let h = build_poly(w);
                if h.len() >= 3 {
                    holes.push(h);
                }
            }
            let outer = build_poly(&fd.outer_wire);
            (outer, holes)
        }
        _ => (Vec::new(), holes),
    }
}

/// Boundary pcurves of a face: for each wire edge, the 2D curve of the edge
/// on this face, its parameter range, and the edge orientation composed with
/// the wire orientation (OCCT FindAPointInTheFace uses BRepAdaptor_Curve2d on
/// the FORWARD-ized face, so the edge orientation is the stored one composed
/// with the wire's — TopExp_Explorer cumOri).
fn collect_boundary_pcurves(face: &Shape) -> Vec<(Curve2d, [f64; 2], Orientation)> {
    let fkey = (face.ptr_id(), face.location);
    let mut out: Vec<(Curve2d, [f64; 2], Orientation)> = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            let mut wire = |w: &Shape| {
                if let TShape::Wire(wd) = &*w.data {
                    let w_or = w.orientation;
                    for e in &wd.edges {
                        if let TShape::Edge(ed) = &*e.data {
                            if let Some((pc, t1, t2)) = ed.pcurves.get(&fkey) {
                                let e_or = Orientation::Forward
                                    .compose(w_or)
                                    .compose(e.orientation);
                                out.push((pc.clone(), [*t1, *t2], e_or));
                            }
                        }
                    }
                }
            };
            wire(&fd.outer_wire);
            for w in &fd.inner_wires {
                wire(w);
            }
        }
        _ => {}
    }
    out
}

/// Boundary edge keys (ptr_id, location) of a face in wire order — the
/// TopExp::MapShapesAndAncestors EDGE->FACE building block.
fn collect_face_edge_keys(face: &Shape) -> Vec<(u64, u32)> {
    let mut out: Vec<(u64, u32)> = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            let mut wire = |w: &Shape| {
                if let TShape::Wire(wd) = &*w.data {
                    for e in &wd.edges {
                        out.push((e.ptr_id(), e.location));
                    }
                }
            };
            wire(&fd.outer_wire);
            for w in &fd.inner_wires {
                wire(w);
            }
        }
        _ => {}
    }
    out
}

/// Nearest intersection parameter s > 0 of the ray `p + s*t` with the face's
/// 2D boundary polygons (OCCT BRepClass_FacePassiveClassifier on the rays of
/// FindAPointInTheFace, SolidExplorer.cxx L113-121).
fn ray_boundary_hit(f: &ExplorerFace, p: DVec2, t: DVec2) -> Option<f64> {
    let (outer, holes) = f.uv_polys.as_ref()?;
    let mut best: Option<f64> = None;
    let mut consider = |poly: &Vec<DVec2>| {
        for i in 0..poly.len() {
            let a = poly[i];
            let b = poly[(i + 1) % poly.len()];
            let ab = b - a;
            let denom = t.x * ab.y - t.y * ab.x;
            if denom.abs() < 1e-30 {
                continue;
            }
            let ap = p - a;
            let s = (ap.x * ab.y - ap.y * ab.x) / denom;
            let u = (ap.x * t.y - ap.y * t.x) / denom;
            if s >= -1e-9 && u >= -1e-9 && u <= 1.0 + 1e-9 {
                if s > 1e-9 && best.map_or(true, |bs| s < bs) {
                    best = Some(s);
                }
            }
        }
    };
    consider(outer);
    for h in holes {
        consider(h);
    }
    best
}

/// True when the UV point lies inside the face's 2D domain: strictly inside
/// the outer polygon and outside every hole (OCCT TopolTool::Classify == IN).
fn uv_in_face_domain(f: &ExplorerFace, uv: DVec2) -> bool {
    let Some((outer, holes)) = &f.uv_polys else {
        return true; // no domain info — accept (planar faces use uv_bounds)
    };
    if outer.len() < 3 {
        return true;
    }
    if !rcad_kernel::base::gprop::tri::point_in_polygon_2d(outer, uv) {
        return false;
    }
    for h in holes {
        if rcad_kernel::base::gprop::tri::point_in_polygon_2d(h, uv) {
            return false;
        }
    }
    true
}

/// OCCT BRepTopAdaptor_TopolTool::Classify(Puv, Tol) — the 2D domain test with
/// a tolerance: a point within Tol of the boundary is ON (accepted as inside),
/// matching IntCurvesFace_Intersector::ClassifyUVPoint (L543-547, Tol=1e-7)
/// and the BRepClass_FaceClassifier used by OtherSegment L586-592.
fn uv_in_face_domain_with_tol(f: &ExplorerFace, uv: DVec2, tol: f64) -> bool {
    let Some((outer, holes)) = &f.uv_polys else {
        return true;
    };
    if outer.len() < 3 {
        return true;
    }
    if rcad_kernel::base::gprop::tri::point_in_polygon_2d(outer, uv) {
        // Strictly inside the outer polygon.
        for h in holes {
            if rcad_kernel::base::gprop::tri::point_in_polygon_2d(h, uv) {
                return false;
            }
        }
        return true;
    }
    // Not strictly inside: ON if within tol of the outer boundary.
    if point_on_polygon_boundary(outer, uv, tol) {
        return true;
    }
    false
}

/// True when `uv` is within `tol` of a segment of the polygon boundary
/// (OCCT TopolTool::Classify == ON).
fn point_on_polygon_boundary(poly: &[DVec2], uv: DVec2, tol: f64) -> bool {
    let n = poly.len();
    if n < 2 {
        return false;
    }
    let tol2 = tol * tol;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let ab = b - a;
        let len2 = ab.dot(ab);
        if len2 < 1e-30 {
            continue;
        }
        let t = ((uv - a).dot(ab) / len2).clamp(0.0, 1.0);
        let q = a + ab * t;
        if (q - uv).length_squared() <= tol2 {
            return true;
        }
    }
    false
}

/// Compute the UV bounding box of a planar face from its wire vertices.
/// OCCT: the face's UV domain (natural restriction or trimmed by wires).
/// Closed circle edges (first == last vertex, e.g. a disk boundary) make the
/// vertex-only sampling degenerate — the edge 3D curves are sampled instead.
fn compute_plane_uv_bounds(face: &Shape, pl: &PlaneGeom, locations: &[glam::DAffine3]) -> Option<[f64; 4]> {
    let mut umin = f64::MAX;
    let mut umax = f64::MIN;
    let mut vmin = f64::MAX;
    let mut vmax = f64::MIN;
    let mut found = false;
    let mut stack: Vec<&Shape> = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            stack.push(&fd.outer_wire);
            for w in &fd.inner_wires {
                stack.push(w);
            }
        }
        _ => return None,
    }
    let mut add_point = |p: &DVec3, found: &mut bool| {
        let d = *p - pl.origin;
        let u = d.dot(pl.u_dir);
        let v = d.dot(pl.v_dir);
        umin = umin.min(u);
        umax = umax.max(u);
        vmin = vmin.min(v);
        vmax = vmax.max(v);
        *found = true;
    };
    while let Some(sh) = stack.pop() {
        match &*sh.data {
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    stack.push(e);
                }
            }
            TShape::Edge(ed) => {
                // Sample the 3D curve (with the edge Location) — a closed
                // circle edge stores first == last, so the vertex endpoints
                // alone cannot bound the domain.
                let loc = ed.first.location;
                if let Some(curve) = &ed.curve {
                    let (t1, t2) = (ed.range[0], ed.range[1]);
                    for i in 0..=16 {
                        let t = t1 + (t2 - t1) * (i as f64) / 16.0;
                        let mut p = curve.point_at(t);
                        if loc != 0 {
                            if let Some(tr) = locations.get(loc as usize) {
                                p = tr.transform_point3(p);
                            }
                        }
                        add_point(&p, &mut found);
                    }
                } else {
                    for v in [&ed.first, &ed.last] {
                        if let TShape::Vertex(vd) = &*v.data {
                            let mut p = vd.point;
                            if v.location != 0 {
                                if let Some(tr) = locations.get(v.location as usize) {
                                    p = tr.transform_point3(p);
                                }
                            }
                            add_point(&p, &mut found);
                        }
                    }
                }
            }
            _ => {}
        }
    }
    if !found {
        return None;
    }
    Some([umin, umax, vmin, vmax])
}

/// True when the 3D point `q` on the plane lies inside the face's UV box
/// (with a tolerance equal to the ray-cast epsilon).
fn in_face_uv(uv_bounds: Option<[f64; 4]>, pl: &PlaneGeom, q: DVec3) -> bool {
    match uv_bounds {
        Some([umin, umax, vmin, vmax]) => {
            let d = q - pl.origin;
            let u = d.dot(pl.u_dir);
            let v = d.dot(pl.v_dir);
            let tol = 1e-6;
            u >= umin - tol && u <= umax + tol && v >= vmin - tol && v <= vmax + tol
        }
        None => true,
    }
}

/// Distance from a point to a 3D curve (OCCT GeomAPI_ProjectPointOnCurve —
/// the minimal distance over the curve).
fn curve_point_distance(curve: &Curve3, p: DVec3) -> f64 {
    let proj = rcad_kernel::base::extrema::closest_point_on_curve(curve, p, 64);
    proj.distance
}

/// UV bounds of a face: the stored bounding box for planar faces, the natural
/// surface parameter bounds otherwise (OCCT BRepAdaptor_Surface
/// FirstUParameter/LastUParameter... in OtherSegment L528-530).
/// Returns (U1, V1, U2, V2) — the OCCT order.
fn face_uv_bounds(f: &ExplorerFace) -> (f64, f64, f64, f64) {
    // The stored uv_bounds is [umin, umax, vmin, vmax] = (U1, U2, V1, V2);
    // reorder to the OCCT (U1, V1, U2, V2).
    if let Some([umin, umax, vmin, vmax]) = f.uv_bounds {
        return (umin, vmin, umax, vmax);
    }
    let [u0, u1, v0, v1] = f.surf.default_domain();
    (u0, v0, u1, v1)
}

/// OCCT IsInfiniteUV (SolidExplorer.cxx L454-467): bit flags for infinite UV
/// bounds (1 = U1, 2 = V1, 4 = U2, 8 = V2).
fn is_infinite_uv(u1: f64, v1: f64, u2: f64, v2: f64) -> i32 {
    let mut val = 0;
    if !u1.is_finite() {
        val |= 1;
    }
    if !v1.is_finite() {
        val |= 2;
    }
    if !u2.is_finite() {
        val |= 4;
    }
    if !v2.is_finite() {
        val |= 8;
    }
    val
}

/// OCCT BRepAdaptor_Surface::IsUPeriodic — the surface type's periodicity
/// (used by PointInTheFace; OCCT queries BRepAdaptor_Surface::IsUPeriodic).
/// rcad surfaces implement is_u_periodic/is_v_periodic directly — kept here
/// only as a type-level fallback; the trait methods are authoritative.
#[allow(dead_code)]
fn surf_is_u_periodic(surf: &Surface3) -> bool {
    surf.is_u_periodic()
}

/// OCCT OtherSegment L675-700: the myParamOnEdge cycle.
fn next_param_on_edge(p: f64) -> f64 {
    if p == 0.512345 {
        0.4
    } else if p == 0.4 {
        0.6
    } else if p == 0.6 {
        0.3
    } else if p == 0.3 {
        0.7
    } else if p == 0.7 {
        0.2
    } else if p == 0.2 {
        0.8
    } else if p == 0.8 {
        0.1
    } else if p == 0.1 {
        0.9
    } else {
        p * 0.5
    }
}

/// Boundary edge keys (ptr_id, location) of a face — for the mapEF
/// (edge -> faces) of the SClassifier.
fn face_edges_of(f: &ExplorerFace) -> Vec<(u64, u32)> {
    f.edge_keys.clone()
}

