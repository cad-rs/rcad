//! TKMath GTest translations — ported from rcad-algorithms.
//! OCCT source: src/FoundationClasses/TKMath/GTests/

use glam::{DVec2, DVec3};

const TOL: f64 = 1e-12;
const TOL_EVAL: f64 = 1e-6;

// --- Minimal BoundingBox (OCCT Bnd_Box equivalent for tests) ---
#[derive(Clone)]
struct BoundingBox { min: DVec3, max: DVec3 }
impl BoundingBox {
    fn new() -> Self { Self { min: DVec3::splat(f64::INFINITY), max: DVec3::splat(f64::NEG_INFINITY) } }
    fn add_point(&mut self, p: DVec3) { self.min = self.min.min(p); self.max = self.max.max(p); }
    fn is_valid(&self) -> bool { self.min.x.is_finite() }
    fn is_empty(&self) -> bool { !self.is_valid() }
    fn center(&self) -> DVec3 { (self.min + self.max) * 0.5 }
    fn add_box(&mut self, other: &Self) { self.min = self.min.min(other.min); self.max = self.max.max(other.max); }
    fn is_out_point(&self, p: DVec3) -> bool { p.x < self.min.x || p.x > self.max.x || p.y < self.min.y || p.y > self.max.y || p.z < self.min.z || p.z > self.max.z }
    fn is_out_box(&self, other: &Self) -> bool { self.max.x < other.min.x || self.min.x > other.max.x || self.max.y < other.min.y || self.min.y > other.max.y || self.max.z < other.min.z || self.min.z > other.max.z }
    fn set_gap(&mut self, _g: f64) {}
    fn contains(&self, p: DVec3) -> bool { !self.is_out_point(p) }
}

// --- Minimal BoundingBox2d (OCCT Bnd_Box2d equivalent for tests) ---
struct BoundingBox2d { min: DVec2, max: DVec2 }
impl BoundingBox2d {
    fn new() -> Self { Self { min: DVec2::splat(f64::INFINITY), max: DVec2::splat(f64::NEG_INFINITY) } }
    fn add_point(&mut self, p: DVec2) { self.min = self.min.min(p); self.max = self.max.max(p); }
    fn is_valid(&self) -> bool { self.min.x.is_finite() }
    fn is_empty(&self) -> bool { !self.is_valid() }
    fn center(&self) -> DVec2 { (self.min + self.max) * 0.5 }
    fn area(&self) -> f64 { let d = self.max - self.min; d.x * d.y }
}

// =============================================================================
// Bnd_Box_Test.cxx
// =============================================================================
mod bnd_box_tests {
    use super::*;
    #[test] fn default_constructor_is_empty() { let bb = BoundingBox::new(); assert!(bb.is_empty()); }
    #[test] fn single_point_box() { let mut bb = BoundingBox::new(); bb.add_point(DVec3::new(1.0, 2.0, 3.0)); assert!((bb.center() - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-10); }
    #[test] fn two_point_box() { let mut bb = BoundingBox::new(); bb.add_point(DVec3::ZERO); bb.add_point(DVec3::new(10.0, 20.0, 30.0)); assert!((bb.center() - DVec3::new(5.0, 10.0, 15.0)).length() < 1e-10); }
    #[test] fn box_contains_point() { let mut bb = BoundingBox::new(); bb.add_point(DVec3::ZERO); bb.add_point(DVec3::new(10.0, 10.0, 10.0)); assert!(bb.contains(DVec3::new(5.0, 5.0, 5.0))); assert!(!bb.contains(DVec3::new(15.0, 5.0, 5.0))); }
    #[test] fn box_contains_itself() { let mut bb = BoundingBox::new(); bb.add_point(DVec3::ZERO); bb.add_point(DVec3::new(10.0, 10.0, 10.0)); assert!(!bb.is_out_box(&bb.clone())); }
    #[test] fn add_box_union() { let mut a = BoundingBox::new(); a.add_point(DVec3::ZERO); a.add_point(DVec3::new(1.0, 1.0, 1.0)); let mut b = BoundingBox::new(); b.add_point(DVec3::new(2.0, 2.0, 2.0)); b.add_point(DVec3::new(3.0, 3.0, 3.0)); a.add_box(&b); assert!(a.contains(DVec3::new(2.5, 2.5, 2.5))); }
    #[test] fn gap_expands_box() { let mut bb = BoundingBox::new(); bb.add_point(DVec3::ZERO); bb.set_gap(1.0); assert!(bb.contains(DVec3::new(0.9, 0.0, 0.0))); }
}

// =============================================================================
// Bnd_Box2d_Test.cxx
// =============================================================================
mod bnd_box2d_tests {
    use super::*;
    #[test] fn default_constructor_is_invalid() { let bb = BoundingBox2d::new(); assert!(!bb.is_valid()); assert!(bb.is_empty()); }
    #[test] fn single_point_box2d() { let mut bb = BoundingBox2d::new(); bb.add_point(DVec2::new(1.0, 2.0)); assert!(bb.is_valid()); assert!(!bb.is_empty()); }
    #[test] fn two_point_box2d() { let mut bb = BoundingBox2d::new(); bb.add_point(DVec2::ZERO); bb.add_point(DVec2::new(10.0, 20.0)); assert!((bb.center() - DVec2::new(5.0, 10.0)).length() < 1e-10); }
    #[test] fn box2d_area() { let mut bb = BoundingBox2d::new(); bb.add_point(DVec2::ZERO); bb.add_point(DVec2::new(10.0, 20.0)); assert!((bb.area() - 200.0).abs() < TOL); }
}

// =============================================================================
// BVH_Box_Test.cxx
// =============================================================================
mod bvh_box_tests {
    use super::*;
    use rcad_kernel::math::bvh::Aabb;
    #[test] fn bvh_aabb_empty() { let a = Aabb::empty(); assert!(a.surface_area() == 0.0); }
    #[test] fn bvh_aabb_single_point() { let a = Aabb::from_points(&[DVec3::new(1.0, 2.0, 3.0)]); assert!((a.min - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-10); }
    #[test] fn bvh_aabb_two_points() { let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(10.0, 20.0, 30.0)]); assert!((a.center() - DVec3::new(5.0, 10.0, 15.0)).length() < 1e-10); }
    #[test] fn bvh_aabb_intersects() { let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(10.0, 10.0, 10.0)]); let b = Aabb::from_points(&[DVec3::new(5.0, 5.0, 5.0), DVec3::new(15.0, 15.0, 15.0)]); assert!(a.intersects(&b)); }
    #[test] fn bvh_aabb_no_intersect() { let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(10.0, 10.0, 10.0)]); let b = Aabb::from_points(&[DVec3::new(20.0, 20.0, 20.0), DVec3::new(30.0, 30.0, 30.0)]); assert!(!a.intersects(&b)); }
}

// =============================================================================
// ElCLib_Test.cxx — ELCLib functions
// =============================================================================
mod elclib_tests {
    use super::*;
    use rcad_kernel::math::el::*;
    use rcad_kernel::geom::CurveEval;
    #[test] fn elclib_line_value_test() { let p = elclib_line_value(2.0, DVec3::ZERO, DVec3::X); assert!((p - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-10); }
    #[test] fn elclib_circle_value_test() { let p = elclib_circle_value(0.0, DVec3::ZERO, DVec3::X, DVec3::Y, 2.0); assert!((p - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-10); }
    #[test] fn elclib_circle_d1_test() { let (p, d1) = elclib_circle_d1(0.0, DVec3::ZERO, DVec3::X, DVec3::Y, 2.0); assert!((p - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-10); assert!((d1 - DVec3::new(0.0, 2.0, 0.0)).length() < TOL_EVAL); }
}

// =============================================================================
// GeomLib_Test.cxx — plane/cylinder/sphere to BSpline conversion
// =============================================================================
mod geom_lib_tests {
    use super::*;
    use rcad_kernel::nurbs_convert::{plane_to_bspline, cylinder_to_bspline, sphere_to_bspline};
    use rcad_kernel::geom::{Plane, CylindricalSurface, SphericalSurface};
    #[test] fn plane_to_bspline_surface() { let p = Plane::new(DVec3::ZERO, DVec3::Z); let bs = plane_to_bspline(&p); assert!(bs.control_points.len() > 0); }
    #[test] fn cylinder_to_bspline_surface() { let c = CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, radius: 2.0, ref_dir: DVec3::X }; let bs = cylinder_to_bspline(&c); assert!(bs.control_points.len() > 0); }
    #[test] fn sphere_to_bspline_surface() { let s = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 2.0, ref_dir: DVec3::X }; let bs = sphere_to_bspline(&s); assert!(bs.control_points.len() > 0); }
}
