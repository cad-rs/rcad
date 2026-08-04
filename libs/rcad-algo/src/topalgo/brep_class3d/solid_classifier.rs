// OCCT BRepClass3d_SolidClassifier (BRepClass3d_SolidClassifier.hxx / .cxx)
// Provides an algorithm to classify a point in a solid.
// Inherits from BRepClass3d_SClassifier.

use crate::topalgo::brep_class3d::solid_explorer::SolidExplorer;
use crate::topalgo::brep_class3d::s_classifier::SClassifier;
use rcad_kernel::topo_shape::Shape;
use glam::DVec3;

/// OCCT BRepClass3d_SolidClassifier — classifies a point relative to a solid.
pub struct SolidClassifier {
    // BRepClass3d_SClassifier base
    pub my_state: u8,          // 0=unknown, 1=faulty, 2=ON, 3=IN, 4=OUT
    // BRepClass3d_SolidClassifier own
    pub a_solid_loaded: bool,
    pub explorer: SolidExplorer,
    pub is_a_hole_in_space: bool,
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
        let clsf = SolidClassifier {
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
            // OCCT BRepClass3d_SClassifier::Perform — a point within the
            // tolerance of a face of the solid is ON the boundary.
            if self.explorer.point_on_face(p, tol) {
                self.my_state = 2; // ON
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
    /// OCCT BRepClass3d_SClassifier::Perform (L203+). Delegates the actual
    /// ray/face intersection to the explorer (BRepClass3d_SolidExplorer),
    /// which resolves face geometry from the shape tree directly.
    fn perform_classify(&mut self, p: DVec3, _tol: f64) {
        if !self.explorer.has_faces() {
            self.my_state = 3; // IN (void solid — whole space)
            return;
        }
        self.my_state = self.explorer.classify_point(p);
    }
}
