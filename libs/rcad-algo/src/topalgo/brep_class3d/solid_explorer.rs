// OCCT BRepClass3d_SolidExplorer (BRepClass3d_SolidExplorer.hxx / .cxx)
// Exploration of a BRep Shape for classification.
// Provides face iteration, bounding box rejection, and BVH tree.

use crate::topalgo::shape_source::ShapeSource;
use rcad_kernel::geom::{Plane as PlaneGeom, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, TShape};
use glam::DVec3;
use std::sync::Arc;

/// A face of the explored solid, with its surface, orientation, and — for
/// planar surfaces — the UV bounding box of the face's wires.
/// OCCT IntCurvesFace_Intersector only counts intersections that fall inside
/// the face's UV domain; without this check the ray-cast would treat the
/// infinite supporting plane as the face.
struct ExplorerFace {
    surf: Surface3,
    ori: Orientation,
    uv_bounds: Option<[f64; 4]>, // [umin, umax, vmin, vmax]
}

/// OCCT BRepClass3d_SolidExplorer — explores a solid's faces for point classification.
pub struct SolidExplorer {
    pub ds: Option<Arc<dyn ShapeSource>>,
    shape: Option<Shape>,
    face_indices: Vec<usize>,
    // Face geometry (surface + orientation) collected from the shape tree.
    // Used for classification without a DS reference (OCCT BRepClass3d
    // explores the TopoDS_Shape directly, never BOPDS).
    face_surfaces: Vec<ExplorerFace>,
}

impl SolidExplorer {
    pub fn new() -> Self {
        SolidExplorer {
            ds: None,
            shape: None,
            face_indices: Vec::new(),
            face_surfaces: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.shape = None;
        self.face_indices.clear();
        self.face_surfaces.clear();
    }

    /// OCCT: InitShape(S) — initialize the explorer with a solid shape.
    /// Traverses the shape tree and collects the surfaces + orientations of
    /// all faces, so the point classification does not depend on the DS.
    pub fn init_shape(&mut self, s: &Shape) {
        self.shape = Some(s.clone());
        self.face_indices.clear();
        self.face_surfaces.clear();
        let mut stack: Vec<Shape> = vec![s.clone()];
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
                    let uv_bounds = fd.surface.as_ref().and_then(|surf| match surf {
                        Surface3::Plane(pl) => compute_plane_uv_bounds(&sh, pl),
                        _ => None,
                    });
                    if let Some(surf) = fd.surface.clone() {
                        self.face_surfaces.push(ExplorerFace {
                            surf,
                            ori: sh.orientation,
                            uv_bounds,
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

    /// Classify point using ray casting (simplified).
    /// OCCT IntCurvesFace_Intersector: the face's orientation flips the
    /// effective surface normal (a reversed face bounds the solid on the
    /// opposite side of its surface). Even-odd rule: a point is IN when a ray
    /// from it crosses the solid boundary an odd number of times. Only
    /// intersections inside the face's UV domain are counted.
    pub fn classify_point(&self, p: DVec3) -> u8 {
        let ray_dir = DVec3::X;
        let mut intersections = 0usize;
        if !self.face_surfaces.is_empty() {
            for f in &self.face_surfaces {
                if let Surface3::Plane(pl) = &f.surf {
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
                }
            }
        }
        if intersections % 2 == 1 { 3 } else { 4 } // IN=3, OUT=4
    }

    /// Set the shape-source reference for face index lookups.
    pub fn set_ds<S: ShapeSource + 'static>(&mut self, ds: &Arc<S>) {
        self.ds = Some(ds.clone() as Arc<dyn ShapeSource>);
    }
}

/// Compute the UV bounding box of a planar face from its wire vertices.
/// OCCT: the face's UV domain (natural restriction or trimmed by wires).
fn compute_plane_uv_bounds(face: &Shape, pl: &PlaneGeom) -> Option<[f64; 4]> {
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
    while let Some(sh) = stack.pop() {
        match &*sh.data {
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    stack.push(e);
                }
            }
            TShape::Edge(ed) => {
                for v in [&ed.first, &ed.last] {
                    if let TShape::Vertex(vd) = &*v.data {
                        let d = vd.point - pl.origin;
                        let u = d.dot(pl.u_dir);
                        let v = d.dot(pl.v_dir);
                        umin = umin.min(u);
                        umax = umax.max(u);
                        vmin = vmin.min(v);
                        vmax = vmax.max(v);
                        found = true;
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
