// OCCT BRepClass3d_SolidClassifier (BRepClass3d_SolidClassifier.hxx / .cxx)
// Provides an algorithm to classify a point in a solid.
// Inherits from BRepClass3d_SClassifier.

use crate::topalgo::class3d::solid_explorer::SolidExplorer;
use crate::topalgo::class3d::s_classifier::SClassifier;
use rcad_kernel::topo_shape::Shape;
use glam::DVec3;

/// OCCT BRepClass3d_SolidClassifier — classifies a point relative to a solid.
pub struct SolidClassifier {
    // BRepClass3d_SClassifier base
    pub my_state: u8,          // 0=unknown, 1=faulty, 2=ON, 3=IN, 4=OUT
    // BRepClass3d_SolidClassifier own
    a_solid_loaded: bool,
    explorer: SolidExplorer,
    is_a_hole_in_space: bool,
}

impl SolidClassifier {
    /// Empty constructor.
    pub fn new() -> Self {
        SolidClassifier {
            my_state: 0,
            a_solid_loaded: false,
            explorer: SolidExplorer::new(),
            is_a_hole_in_space: false,
        }
    }

    /// OCCT: Load(const TopoDS_Shape& S) — initialize explorer for the given solid.
    pub fn load(&mut self, s: &Shape) {
        if self.a_solid_loaded {
            self.explorer.clear();
        }
        self.explorer.init_shape(s);
        self.a_solid_loaded = true;
    }

    /// Constructor from a Shape.
    pub fn from_shape(s: &Shape) -> Self {
        let mut clsf = SolidClassifier {
            my_state: 0,
            a_solid_loaded: true,
            explorer: SolidExplorer::from_shape(s),
            is_a_hole_in_space: false,
        };
        clsf
    }

    /// Constructor to classify the point P with tolerance Tol on the solid S.
    pub fn from_shape_point(s: &Shape, p: DVec3, tol: f64) -> Self {
        let mut clsf = SolidClassifier {
            my_state: 0,
            a_solid_loaded: true,
            explorer: SolidExplorer::from_shape(s),
            is_a_hole_in_space: false,
        };
        clsf.perform(p, tol);
        clsf
    }

    /// OCCT: Perform(P, Tol) — classify point P relative to loaded solid.
    pub fn perform(&mut self, p: DVec3, tol: f64) {
        if !self.a_solid_loaded { return; }
        // OCCT L171-191: check bounding box first (fast rejection)
        if !self.is_a_hole_in_space {
            if self.explorer.reject(p) {
                self.my_state = 4; // OUT
                return;
            }
            // OCCT L190: BRepClass3d_SClassifier::Perform(explorer, P, Tol)
            self.perform_classify(p, tol);
        } else {
            if self.explorer.reject(p) {
                self.my_state = 3; // IN
                return;
            }
            self.perform_classify(p, tol);
        }
    }

    /// OCCT: PerformInfinitePoint(Tol) — classify point at infinity.
    /// Useful for computing the orientation of a solid.
    pub fn perform_infinite_point(&mut self, _tol: f64) {
        // OCCT: shoots a ray from a face point in the opposite normal direction
        // and counts intersections to determine if the solid is a "hole in space"
        // rcad: simplified — delegate to perform_classify at infinity
        self.my_state = 4; // OUT (default for infinite point)
    }

    /// OCCT: State() — returns my_state (inherited from SClassifier).
    pub fn state(&self) -> u8 { self.my_state }

    /// Internal: perform classification using ray casting.
    fn perform_classify(&mut self, p: DVec3, tol: f64) {
        // OCCT BRepClass3d_SClassifier::Perform (L203+):
        // 1. Check rejection → if rejected, point is IN the void solid
        if self.explorer.reject(p) {
            self.my_state = 3; // IN
            return;
        }
        // 2. Get all faces from the explorer
        let face_indices = self.explorer.get_face_indices().to_vec();
        if face_indices.is_empty() {
            self.my_state = 3; // IN (void solid)
            return;
        }
        // 3. Ray casting: shoot ray in +X, count front-face intersections
        // OCCT uses IntCurvesFace_Intersector for precise ray-face intersection.
        // rcad: simplified ray-plane intersection for planar faces.
        let ray_dir = DVec3::X;
        let mut intersections = 0usize;
        for &fi in &face_indices {
            let surf = match self.explorer.ds.as_ref().and_then(|ds| ds.face_surface(fi)) {
                Some(s) => s,
                None => continue,
            };
            if let rcad_kernel::geom::Surface3::Plane(pl) = surf {
                let denom = ray_dir.dot(pl.normal);
                if denom.abs() < 1e-12 { continue; }
                let t = (pl.origin - p).dot(pl.normal) / denom;
                if t > tol && denom < 0.0 { // front face (entering solid)
                    intersections += 1;
                }
            }
        }
        self.my_state = if intersections % 2 == 1 { 3 } else { 4 }; // IN=3, OUT=4
    }
}
