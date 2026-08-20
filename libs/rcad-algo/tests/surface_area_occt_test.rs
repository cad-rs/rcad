//! OCCT BRepGProp::SurfaceProperties (BRepGProp.cxx L167-266) alignment.
//!
//! `checkprops result -s` (BRepTest_GPropCommands.cxx L123-125) calls
//! BRepGProp::SurfaceProperties(S, G, SkipShared=false, UseTriangulation=false)
//! → surfaceProperties(S, Props, Eps=1.0, false, false).  Per face:
//!   - natural-restriction faces (no wires) → G.Perform(BF): whole-surface
//!     Gauss-Legendre over the surface natural bounds (BRepGProp_Gauss.cxx
//!     L1306-1393),
//!   - faces with wires → G.Perform(BF, BD): Green-theorem line integral over
//!     the edge pcurves with the face UV bounds as the U offset
//!     (BRepGProp_Gauss.cxx L1126-1211, BRepTools.cxx L64-367).
//!
//! The sphere test fixes the semantics of `surface_area` on a full sphere
//! (a natural-restriction face, integrated by the whole-surface Gauss path):
//! exactly 4*PI*R^2.

use glam::DVec3;
use rcad_kernel::base::gprop::surface::surface_area;

/// OCCT `psphere 1`: a full sphere is a natural-restriction face (no wires),
/// so BRepGProp::SurfaceProperties integrates the whole surface by
/// Gauss-Legendre (G.Perform(BF)); the surface area is exactly 4*PI*R^2.
#[test]
fn sphere_surface_area_matches_occt() {
    let sphere = rcad_modeling::make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let expected = 4.0 * std::f64::consts::PI;
    let got = surface_area(&sphere);
    assert!(
        (got - expected).abs() < 1e-6,
        "sphere surface_area = {got}, expected {expected}"
    );
}
