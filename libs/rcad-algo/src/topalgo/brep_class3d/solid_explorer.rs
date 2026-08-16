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
use std::collections::VecDeque;
use std::sync::Arc;

/// A face of the explored solid, with its surface, orientation, and — for
/// planar surfaces — the UV bounding box of the face's wires.
/// OCCT IntCurvesFace_Intersector only counts intersections that fall inside
/// the face's UV domain; without this check the ray-cast would treat the
/// infinite supporting plane as the face.
pub(crate) struct ExplorerFace {
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
}

impl Clone for ExplorerFace {
    fn clone(&self) -> Self {
        ExplorerFace {
            surf: self.surf.clone(),
            ori: self.ori,
            uv_bounds: self.uv_bounds,
            uv_polys: self.uv_polys.clone(),
            boundary: self.boundary.clone(),
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
    /// Merged TopLoc_Location table (index 0 = identity); empty when the
    /// explorer was built without location support.
    locations: Vec<glam::DAffine3>,
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
            locations: vec![glam::DAffine3::IDENTITY],
        }
    }

    pub fn clear(&mut self) {
        self.shape = None;
        self.face_indices.clear();
        self.face_surfaces.clear();
        self.vertices.clear();
        self.edges.clear();
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
        let mut queue: VecDeque<Shape> = VecDeque::new();
        queue.push_back(s.clone());
        while let Some(sh) = queue.pop_front() {
            match &*sh.data {
                TShape::Solid(sd) => {
                    for x in &sd.shells {
                        queue.push_back(x.clone());
                    }
                }
                TShape::CompSolid(cd) => {
                    for x in cd {
                        queue.push_back(x.clone());
                    }
                }
                TShape::Compound(cd) => {
                    for x in cd {
                        queue.push_back(x.clone());
                    }
                }
                TShape::Shell(sd) => {
                    for x in &sd.faces {
                        queue.push_back(x.clone());
                    }
                }
                TShape::Vertex(vd) => {
                    let loc = self
                        .locations
                        .get(sh.location as usize)
                        .copied()
                        .unwrap_or(glam::DAffine3::IDENTITY);
                    self.vertices.push(loc.transform_point3(vd.point));
                }
                TShape::Edge(ed) => {
                    let loc = self
                        .locations
                        .get(sh.location as usize)
                        .copied()
                        .unwrap_or(glam::DAffine3::IDENTITY);
                    self.edges.push(ed.curve.as_ref().map(|c| rcad_kernel::geom::transform_curve(c, &loc)));
                }
                TShape::Face(fd) => {
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
                    if let Some(surf) = surface {
                        self.face_surfaces.push(ExplorerFace {
                            surf,
                            ori: sh.orientation,
                            uv_bounds,
                            uv_polys,
                            boundary,
                        });
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
            locations: vec![glam::DAffine3::IDENTITY],
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
            locations: if locations.is_empty() {
                vec![glam::DAffine3::IDENTITY]
            } else {
                locations.to_vec()
            },
        };
        exp.init_shape(s);
        exp
    }

    /// OCCT: Reject(P) — fast bounding box rejection.
    /// Returns true if P is definitely outside the solid.
    pub fn reject(&self, _p: DVec3) -> bool {
        // OCCT uses Bnd_Box from the shape's bounding volume.
        // rcad: simplified — no bounding box check.
        false
    }

    /// OCCT BRepClass3d_SClassifier::Perform (L212-228) — a point within the
    /// tolerance of a VERTEX or EDGE of the solid is ON the boundary
    /// (aMapEV tree selection); the distance to a face interior is NOT an ON
    /// test. rcad: vertices by point distance, edges by curve distance.
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
    /// The single fixed ray direction is degenerate when it lies in a face's
    /// plane (e.g. a point on the plane of a planar face, whose ray grazes the
    /// neighboring curved faces' rims): the rim hits are rejected by the UV
    /// domain check and the count comes out ZERO even for an inside point.
    /// OCCT retries with a different line direction (BRepClass3d_SClassifier.cxx
    /// L264-285 — the faulty-line loop over SolidExplorer::Segment/
    /// OtherSegment). The multi-direction retry (K1 workstream) is not enabled
    /// yet — it needs the ON-state face-distance check first (the e1 infinite
    /// point lies ON a lateral face; see the session handover). The plain +X
    /// ray keeps the non-degenerate cases as before.
    pub fn classify_point(&self, p: DVec3) -> u8 {
        let ray_dir = DVec3::X;
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
                    _ => Self::ray_face_param(p, ray_dir, f),
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
                        surf: surf.clone(),
                        ori: face_ori,
                        uv_bounds: None,
                        uv_polys: None,
                        boundary: vec![],
                    };
                    if let Some(t) = Self::ray_face_param(p, ray_dir, &ef) {
                        if t > 1e-7 {
                            intersections += 1;
                        }
                    }
                }
            }
        }
        if intersections % 2 == 1 { 3 } else { 4 } // IN=3, OUT=4
    }

    /// Parameter of the closest intersection of the ray `p + t*dir` with a
    /// curved face (OCCT BRepClass3d_SClassifier uses IntCurvesFace_Intersector
    /// for the Line x Face intersection; BeanFaceIntersector is its translation).
    /// Only intersections whose UV lies inside the face's 2D domain are counted
    /// (OCCT IntCurvesFace_Intersector.cxx L256-286: Classify(Puv) == IN/ON).
    pub(crate) fn ray_face_param(p: DVec3, dir: DVec3, f: &ExplorerFace) -> Option<f64> {
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
        // (IntCurvesFace_Intersector's currentstate IN/ON filter).
        let mut parmin = f64::MAX;
        for r in bfi.result() {
            let t = r.first();
            let q = line_for_points.point_at(t);
            let (uv, _q2) = crate::bop::closest_point_on_surface(&f.surf, q);
            if !uv_in_face_domain(f, uv) {
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
    pub(crate) fn face_point(f: &ExplorerFace, param: f64) -> Option<(DVec3, f64, f64)> {
        if let Surface3::Plane(pl) = &f.surf {
            // Planar faces: the UV-box interpolation point is inside the face
            // when the box is the face's bounding box (the vertex projection
            // box; the probe parameter samples [0.1, 0.9] of it). Verify with
            // the domain check and fall back to the edge-walk otherwise.
            if let Some([umin, umax, vmin, vmax]) = f.uv_bounds {
                let uv = DVec2::new(umin + (umax - umin) * param, vmin + (vmax - vmin) * param);
                if uv_in_face_domain(f, uv) {
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
                if !uv_in_face_domain(f, uv) {
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
fn build_uv_polys(face: &Shape, surf: &Surface3, locations: &[glam::DAffine3]) -> (Vec<DVec2>, Vec<Vec<DVec2>>) {    let build_poly = |w: &Shape| -> Vec<DVec2> {
        let mut poly: Vec<DVec2> = Vec::new();
        if let TShape::Wire(wd) = &*w.data {
            for e in &wd.edges {
                if let TShape::Edge(ed) = &*e.data {
                    if let TShape::Vertex(vd) = &*ed.last.data {
                        let (uv, _) = closest_point_on_surface(surf, vd.point);
                        poly.push(uv);
                    }
                }
            }
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

