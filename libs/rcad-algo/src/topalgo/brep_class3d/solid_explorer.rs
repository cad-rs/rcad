// OCCT BRepClass3d_SolidExplorer (BRepClass3d_SolidExplorer.hxx / .cxx)
// Exploration of a BRep Shape for classification.
// Provides face iteration, bounding box rejection, and BVH tree.

use crate::bop::ds::DS;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::TShape;
use glam::DVec3;
use std::sync::Arc;

/// OCCT BRepClass3d_SolidExplorer — explores a solid's faces for point classification.
pub struct SolidExplorer {
    pub ds: Option<Arc<DS>>,
    shape: Option<Shape>,
    face_indices: Vec<usize>,
    // BVH tree for fast rejection (OCCT: UBTree)
    // rcad: simplified — use direct face iteration
}

impl SolidExplorer {
    pub fn new() -> Self {
        SolidExplorer { ds: None, shape: None, face_indices: Vec::new() }
    }

    pub fn clear(&mut self) {
        self.shape = None;
        self.face_indices.clear();
    }

    /// OCCT: InitShape(S) — initialize the explorer with a solid shape.
    pub fn init_shape(&mut self, s: &Shape) {
        self.shape = Some(s.clone());
        self.face_indices.clear();
        // Collect face DS indices from the solid's sub-shapes
        // rcad: this requires the DS to look up indices
    }

    /// Constructor from Shape.
    pub fn from_shape(s: &Shape) -> Self {
        let mut exp = SolidExplorer { ds: None, shape: Some(s.clone()), face_indices: Vec::new() };
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

    /// Add a face index (used by IntToolsContext when building from DS).
    pub fn add_face_index(&mut self, fi: usize) {
        self.face_indices.push(fi);
    }

    /// Classify point using ray casting (simplified).
    /// OCCT IntCurvesFace_Intersector: the face's orientation flips the
    /// effective surface normal (a reversed face bounds the solid on the
    /// opposite side of its surface).
    pub fn classify_point(&self, p: DVec3) -> u8 {
        if let Some(ref ds) = self.ds {
            let mut intersections = 0usize;
            for &fi in &self.face_indices {
                let surf = match ds.face_surface(fi) {
                    Some(s) => s,
                    None => continue,
                };
                let face_ori = ds.shape_info(fi).shape.orientation;
                if let rcad_kernel::geom::Surface3::Plane(pl) = surf {
                    let normal = if face_ori == rcad_kernel::topods::Orientation::Reversed {
                        -pl.normal
                    } else {
                        pl.normal
                    };
                    let ray_dir = DVec3::X;
                    let denom = ray_dir.dot(normal);
                    if denom.abs() < 1e-12 { continue; }
                    let t = (pl.origin - p).dot(normal) / denom;
                    if t > 1e-7 && denom < 0.0 {
                        intersections += 1;
                    }
                }
            }
            if intersections % 2 == 1 { 3 } else { 4 } // IN=3, OUT=4
        } else {
            4 // OUT
        }
    }

    /// Set the DS reference for face index lookups.
    pub fn set_ds(&mut self, ds: &Arc<DS>) {
        self.ds = Some(ds.clone());
    }
}
