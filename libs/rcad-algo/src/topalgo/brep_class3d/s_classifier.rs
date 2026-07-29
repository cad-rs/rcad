// OCCT BRepClass3d_SClassifier (BRepClass3d_SClassifier.hxx / .cxx)
// Base class for solid classification via ray casting.

use crate::topalgo::brep_class3d::solid_explorer::SolidExplorer;
use glam::DVec3;

/// OCCT BRepClass3d_SClassifier — base classification algorithm.
/// Uses ray casting to determine if a point is IN/OUT/ON relative to a solid.
pub struct SClassifier {
    pub my_state: u8,  // 0=unknown, 1=faulty, 2=ON, 3=IN, 4=OUT
}

impl SClassifier {
    pub fn new() -> Self {
        SClassifier { my_state: 0 }
    }

    /// OCCT L203-400: Perform(SolidExplorer, P, Tol)
    pub fn perform(&mut self, explorer: &SolidExplorer, p: DVec3, _tol: f64) {
        // OCCT L207-212: Reject check
        if explorer.reject(p) {
            self.my_state = 3; // IN (void solid — whole space)
            return;
        }

        // OCCT L214-230: BVH vertex/edge proximity check → ON
        // rcad: skip BVH tree, delegate to ray casting

        // OCCT L232-400: ray casting with transition tracking
        self.my_state = explorer.classify_point(p);
    }
}
