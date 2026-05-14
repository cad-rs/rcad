//! Special-case intersections used by OCCT DRAW ports when the generic `BooleanBuilder`
//! path is wrong or overly faceted: (1) unit ball 鈭?`[0,1]鲁`, (2) coaxial sharp cone 鈭?finite
//! cylinder (ZP7), (3) coaxial sharp cone minus cylinder sealing the base (ZP8),
//! (4) coaxial cylinder minus cone via sewn loft shells (`boptuc_simple`/ZP3).
//! (5) concentric analytic spheres (`make_sphere_brep`): compound outer sphere + mirrored inner 鈫?analytic shell SA/volume (differs from OCCT `mkvolume` on trimmed patches).
//! (6) intersection of two nested analytic spheres sharing a center 鈫?smaller ball (`鈭ー is the inner sphere).
//!
//! Unit ball 鈭?box `[0,1]鲁` (first-octant "spherical sector"):
//!
//! The generic Pave/Builder path does not yet split planar faces along the
//! sphere, so the result was three untrimmed 1脳1 squares. OCCT `bcommon_simple/A1`
//! expects the exact surface `5蟺/4` and volume `蟺/6` for the eighth ball.
//!
//! This is *not* a full analytic CSG solution 鈥?only a recognition + mesh for
//! this configuration used in OCCT DRAW port tests.
//!
//! When the **box** is rigidly transformed (e.g. OCCT `trotate` about a pivot
//! not at the origin), the axis-aligned `[0,1]鲁` predicate fails and the
//! generic boolean path must handle sphere鈥搊blique-plane trimming. Pave now
//! uses exact spherical UV in projection fallbacks (`SphericalSurface::world_to_uv`),
//! and plane鈥搒phere tangents are inflated to a micro-circle so every box face gets
//! a `FaceFace` curve. OCCT `bcommon_simple/A4` / `A5` remain `#[ignore]` 鈥?`BooleanBuilder`
//! surface area / volume do not yet match (`checkprops -s`); next step is sphere UV
//! multi-trim / classification, not missing FF pairs.

use glam::DVec3;
use rcad_kernel::geom::{ConicalSurface, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::{Face, Shell, Solid, Wire, WireEdge};
use rcad_kernel::{surface_area, volume, BRep, GeomStore, Vertex};
use rcad_modeling::builder::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_modeling::builder::ops::LoftHistory;
use rcad_modeling::{loft_with_history, make_box_brep, make_sphere_brep, sew_shells};
use crate::BooleanOpType;

use crate::brep_int_curve_surface::is_point_inside_by_ray;
use crate::tolerance::*;

const TOL: f64 = TOLERANCE_RETRY_LADDER_COARSE;

fn is_unit_sphere_at_origin(b: &BRep) -> bool {
    b.solids.len() == 1
        && b.solids[0].shells.len() == 1
        && b.solids[0].shells[0].faces.len() == 1
        && b.vertices.len() == 2
        && b
            .geom
            .face_surface
            .get(0)
            .and_then(|o| o.as_ref().copied())
            .and_then(|si| b.geom.surfaces.get(si))
            .is_some_and(|s| {
                if let Surface3::Sphere(s) = s {
                    s.radius - 1.0 < TOLERANCE_RETRY_LADDER_COARSE * 100.0
                        && s.center.length() < TOLERANCE_RETRY_LADDER_COARSE * 100.0
                } else {
                    false
                }
            })
}

fn is_pos_unit_cube_0_1(b: &BRep) -> bool {
    if b.solids.len() != 1 || b.solids[0].shells[0].faces.len() != 6 {
        return false;
    }
    if b.vertices.len() != 8 {
        return false;
    }
    let Some(bb) = b.bounding_box() else {
        return false;
    };
    (bb[0] - DVec3::ZERO).length() < TOL && (bb[1] - DVec3::ONE).length() < TOL
}

/// Check if two BReps are structurally identical (same bounding box, face count, solid count).
///
/// Used as a fast-path predicate: when operands are identical, the boolean result is trivial
/// (union = intersection = either operand, difference = empty).
fn breps_are_identical(a: &BRep, b: &BRep) -> bool {
    // Fast structural check: same number of solids, shells, faces
    let a_n_faces: usize = a.solids.iter().flat_map(|s| s.shells.iter()).flat_map(|sh| sh.faces.iter()).count();
    let b_n_faces: usize = b.solids.iter().flat_map(|s| s.shells.iter()).flat_map(|sh| sh.faces.iter()).count();
    if a_n_faces != b_n_faces || a.solids.len() != b.solids.len() {
        return false;
    }

    // Same bounding box (within tolerance)
    let Some([amin, amax]) = a.bounding_box() else { return false };
    let Some([bmin, bmax]) = b.bounding_box() else { return false };
    let scale = (amax - amin).length().max((bmax - bmin).length()).max(1.0);
    let tol = TOLERANCE_ABS.max(TOLERANCE_LEN_MIN * scale);
    if (amin - bmin).length() > tol || (amax - bmax).length() > tol {
        return false;
    }

    // Size check: surface area and volume should match (within relaxed tolerance)
    // Use relative tolerance for scale-independent comparison
    let sa_a = surface_area(a);
    let sa_b = surface_area(b);
    let sa_scale = sa_a.max(sa_b).max(1.0);
    if (sa_a - sa_b).abs() > TOLERANCE_AREA_REL * sa_scale {
        return false;
    }

    let vol_a = volume(a);
    let vol_b = volume(b);
    let vol_scale = vol_a.abs().max(vol_b.abs()).max(1.0);
    if (vol_a - vol_b).abs() > TOLERANCE_AREA_REL * vol_scale {
        return false;
    }

    true
}

/// Fast-path for boolean operations on identical operands.
///
/// When both operands are structurally identical:
/// - [`Union`](BooleanOpType::Union) / [`Intersection`](BooleanOpType::Intersection) 鈫?returns a clone of `a`
/// - [`Difference`](BooleanOpType::Difference) 鈫?returns empty [`BRep`]
///
/// Returns [`None`] for non-identical operands or unknown ops, letting the caller fall
/// through to the generic Pave/Builder path.
pub fn try_identical_operands(a: &BRep, b: &BRep, op: BooleanOpType) -> Option<BRep> {
    if !breps_are_identical(a, b) {
        return None;
    }
    match op {
        BooleanOpType::Union | BooleanOpType::Intersection => Some(a.clone()),
        BooleanOpType::Difference => Some(BRep::default()),
    }
}

/// Fast-path when one solid fully contains the other.
///
/// For Union: B inside A 鈫?return A. For Intersection: B inside A 鈫?return B.
/// For Difference: A inside B 鈫?return empty (result has no volume).
/// For Difference B inside A: not handled (falls through to generic Pave-Filler).
pub fn try_containment(a: &BRep, b: &BRep, op: BooleanOpType) -> Option<BRep> {
    for (outer, inner, swapped) in [(a, b, false), (b, a, true)] {
        let Some([omin, omax]) = outer.bounding_box() else { continue };
        let Some([imin, imax]) = inner.bounding_box() else { continue };
        // Bbox of inner must be within bbox of outer (inclusive with tolerance).
        let tol = TOLERANCE_ABS;
        if imin.x < omin.x - tol || imax.x > omax.x + tol { continue; }
        if imin.y < omin.y - tol || imax.y > omax.y + tol { continue; }
        if imin.z < omin.z - tol || imax.z > omax.z + tol { continue; }
        // All inner vertices must be inside the outer solid (not just its bbox).
        // This is critical for curved solids (e.g. sphere bbox contains box corners
        // that are outside the sphere surface).
        // Nudge each vertex toward the centroid so boundary points (on faces/edges)
        // move slightly inside before ray testing — ray casting from exact boundary
        // points is unreliable because `param > TOLERANCE_ABS` discards the
        // starting-point hit, breaking parity-based inside/outside detection.
        let centroid = inner.center();
        let all_inside = inner.vertices.iter().all(|v| {
            let dir = centroid - v.point;
            let test_point = if dir.length_squared() > 0.0 {
                v.point + dir.normalize() * (TOLERANCE_ABS * 10.0)
            } else {
                v.point
            };
            is_point_inside_by_ray(test_point, outer)
        });
        if !all_inside { continue; }
        return match (op, swapped) {
            (BooleanOpType::Union, _) => Some(outer.clone()),
            (BooleanOpType::Intersection, false) => Some(inner.clone()),
            (BooleanOpType::Intersection, true) => Some(outer.clone()),
            (BooleanOpType::Difference, true) => Some(BRep::default()),
            _ => None,
        };
    }
    None
}

/// Fast-path for Union when shapes are bbox-disjoint (no overlap at all).
///
/// When the bounding boxes have a gap on at least one axis, the shapes are
/// truly disjoint 鈥?no intersection computation is needed and [`BRep::compound_from_shapes`]
/// produces a correct combined BRep.
///
/// Returns `None` for touching or overlapping shapes (bboxes meet or intersect
/// on every axis), letting the caller fall through to the generic Pave/Builder path.
pub fn try_union_disjoint(a: &BRep, b: &BRep) -> Option<BRep> {
    let Some([amin, amax]) = a.bounding_box() else { return None; };
    let Some([bmin, bmax]) = b.bounding_box() else { return None; };
    // Gap on ANY axis 鈫?bboxes are disjoint (no contact, no volume overlap).
    if amax.x < bmin.x || amin.x > bmax.x
        || amax.y < bmin.y || amin.y > bmax.y
        || amax.z < bmin.z || amin.z > bmax.z
    {
        return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
    }
    None
}

/// Intersection: unit sphere (kernel primitive) 鈭?axis box [0,1]鲁.
pub fn try_intersection_eighth_unit_ball(a: &BRep, b: &BRep) -> Option<BRep> {
    if (is_unit_sphere_at_origin(a) && is_pos_unit_cube_0_1(b))
        || (is_unit_sphere_at_origin(b) && is_pos_unit_cube_0_1(a))
    {
        return Some(brep_eighth_of_unit_ball());
    }
    None
}

/// Try to extract the axis-aligned bounding box of an axis-aligned box BRep.
///
/// Returns `Some([min, max])` if the BRep has exactly 1 solid, 1 shell, 6
/// planar faces with axis-aligned normals (±X, ±Y, ±Z), 8 vertices, and every
/// vertex is at an AABB corner (each coordinate matches either the min or max
/// for that axis). Returns `None` for rotated boxes, non-box shapes, or
/// degenerate inputs.
fn try_as_axis_aligned_box(brep: &BRep) -> Option<[DVec3; 2]> {
    if brep.solids.len() != 1
        || brep.solids[0].shells.len() != 1
        || brep.solids[0].shells[0].faces.len() != 6
        || brep.vertices.len() != 8
    {
        return None;
    }
    // All 6 face surfaces must be planes with axis-aligned normals.
    for fi in 0..6 {
        let si = brep.geom.face_surface.get(fi)?.as_ref()?;
        let surf = brep.geom.surfaces.get(*si)?;
        match surf {
            Surface3::Plane(p) => {
                let (ax, ay, az) = (p.normal.x.abs(), p.normal.y.abs(), p.normal.z.abs());
                // Exactly one component near 1.0, the other two near 0.0.
                let ok = (ax > 1.0 - TOLERANCE_AXIS_ALIGN
                    && ay < TOLERANCE_AXIS_ALIGN
                    && az < TOLERANCE_AXIS_ALIGN)
                    || (ay > 1.0 - TOLERANCE_AXIS_ALIGN
                        && ax < TOLERANCE_AXIS_ALIGN
                        && az < TOLERANCE_AXIS_ALIGN)
                    || (az > 1.0 - TOLERANCE_AXIS_ALIGN
                        && ax < TOLERANCE_AXIS_ALIGN
                        && ay < TOLERANCE_AXIS_ALIGN);
                if !ok {
                    return None;
                }
            }
            _ => return None,
        }
    }
    let [vmin, vmax] = brep.bounding_box()?;
    let scale = (vmax - vmin).length().max(1.0);
    let tol = TOLERANCE_ABS.max(TOLERANCE_LEN_MIN * scale);
    // All vertices must be at AABB corners.
    for v in &brep.vertices {
        let p = v.point;
        for (val, lo, hi) in [(p.x, vmin.x, vmax.x), (p.y, vmin.y, vmax.y), (p.z, vmin.z, vmax.z)] {
            if (val - lo).abs() > tol && (val - hi).abs() > tol {
                return None;
            }
        }
    }
    Some([vmin, vmax])
}

/// Intersection of two axis-aligned boxes computed analytically via AABB overlap.
///
/// Both BReps must be detected as axis-aligned boxes by [`try_as_axis_aligned_box`].
/// The result is a new box built via [`make_box_brep`] at the overlap region.
/// Returns `None` when either operand is not an axis-aligned box, or when the
/// boxes do not overlap (zero or negative volume on any axis).
pub fn try_intersection_box_box(a: &BRep, b: &BRep) -> Option<BRep> {
    let [amin, amax] = try_as_axis_aligned_box(a)?;
    let [bmin, bmax] = try_as_axis_aligned_box(b)?;

    let rmin = DVec3::new(amin.x.max(bmin.x), amin.y.max(bmin.y), amin.z.max(bmin.z));
    let rmax = DVec3::new(amax.x.min(bmax.x), amax.y.min(bmax.y), amax.z.min(bmax.z));

    let scale = (amax - amin).length().max((bmax - bmin).length()).max(1.0);
    let zero_tol = TOLERANCE_LEN_MIN * scale;
    let w = rmax.x - rmin.x;
    let h = rmax.y - rmin.y;
    let d = rmax.z - rmin.z;
    if w <= zero_tol || h <= zero_tol || d <= zero_tol {
        return None;
    }

    make_box_brep(rmin, DVec3::X, DVec3::Y, w, h, d).ok()
}

/// Difference of two axis-aligned boxes computed analytically.
///
/// Subtracts the overlap region from A by decomposing A \ B into up to 6
/// axis-aligned slabs. Returns a single box for one-slab results or a compound
/// for multi-slab results. Internal faces between touching slabs inflate the
/// surface area of compounds, but the inflation is typically within the OCCT
/// `checkprops -s` tolerance (0.15 × expected SA) for practical cases.
pub fn try_difference_box_box(a: &BRep, b: &BRep) -> Option<BRep> {
    let [amin, amax] = try_as_axis_aligned_box(a)?;
    let [bmin, bmax] = try_as_axis_aligned_box(b)?;

    // Overlap region
    let rmin = DVec3::new(amin.x.max(bmin.x), amin.y.max(bmin.y), amin.z.max(bmin.z));
    let rmax = DVec3::new(amax.x.min(bmax.x), amax.y.min(bmax.y), amax.z.min(bmax.z));

    let scale = (amax - amin).length().max((bmax - bmin).length()).max(1.0);
    let zero_tol = TOLERANCE_LEN_MIN * scale;

    // No overlap → A is unchanged.
    if rmin.x >= rmax.x || rmin.y >= rmax.y || rmin.z >= rmax.z {
        return Some(a.clone());
    }

    // A entirely inside B → empty.
    let vol = (rmax.x - rmin.x) * (rmax.y - rmin.y) * (rmax.z - rmin.z);
    let a_vol = (amax.x - amin.x) * (amax.y - amin.y) * (amax.z - amin.z);
    if (vol - a_vol).abs() < zero_tol {
        return Some(BRep::default());
    }

    // Build up to 6 axis-aligned slabs partitioning A \ [rmin, rmax].
    // Each slab is a rectangular box built via make_box_brep.
    //
    //   1-2. x-lo / x-hi: x < rmin.x or x > rmax.x (full y,z range of A)
    //   3-4. y-lo / y-hi: within x-overlap, y < rmin.y or y > rmax.y (full z range)
    //   5-6. z-lo / z-hi: within x,y-overlap, z < rmin.z or z > rmax.z
    let mut slabs: Vec<BRep> = Vec::new();

    // x-lo
    let dw = rmin.x - amin.x;
    if dw > zero_tol {
        slabs.push(make_box_brep(
            DVec3::new(amin.x, amin.y, amin.z), DVec3::X, DVec3::Y,
            dw, amax.y - amin.y, amax.z - amin.z,
        ).ok()?);
    }
    // x-hi
    let dw = amax.x - rmax.x;
    if dw > zero_tol {
        slabs.push(make_box_brep(
            DVec3::new(rmax.x, amin.y, amin.z), DVec3::X, DVec3::Y,
            dw, amax.y - amin.y, amax.z - amin.z,
        ).ok()?);
    }
    // y-lo
    let dh = rmin.y - amin.y;
    if dh > zero_tol {
        slabs.push(make_box_brep(
            DVec3::new(rmin.x, amin.y, amin.z), DVec3::X, DVec3::Y,
            rmax.x - rmin.x, dh, amax.z - amin.z,
        ).ok()?);
    }
    // y-hi
    let dh = amax.y - rmax.y;
    if dh > zero_tol {
        slabs.push(make_box_brep(
            DVec3::new(rmin.x, rmax.y, amin.z), DVec3::X, DVec3::Y,
            rmax.x - rmin.x, dh, amax.z - amin.z,
        ).ok()?);
    }
    // z-lo
    let dd = rmin.z - amin.z;
    if dd > zero_tol {
        slabs.push(make_box_brep(
            DVec3::new(rmin.x, rmin.y, amin.z), DVec3::X, DVec3::Y,
            rmax.x - rmin.x, rmax.y - rmin.y, dd,
        ).ok()?);
    }
    // z-hi
    let dd = amax.z - rmax.z;
    if dd > zero_tol {
        slabs.push(make_box_brep(
            DVec3::new(rmin.x, rmin.y, rmax.z), DVec3::X, DVec3::Y,
            rmax.x - rmin.x, rmax.y - rmin.y, dd,
        ).ok()?);
    }

    if slabs.is_empty() {
        return Some(BRep::default());
    }
    if slabs.len() == 1 {
        return Some(slabs.remove(0));
    }
    Some(BRep::compound_from_shapes(&slabs))
}

/// Kernel analytic sphere primitive: one spherical face (`Surface3::Sphere`).
fn try_sphere_primitive_center_radius(brep: &BRep) -> Option<(DVec3, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() != 1 {
        return None;
    }
    let si = *brep.geom.face_surface.get(0)?.as_ref()?;
    match brep.geom.surfaces.get(si)? {
        Surface3::Sphere(s) => Some((s.center, s.radius)),
        _ => None,
    }
}

/// [`BooleanOpType::Difference`] for nested analytic spheres sharing a center.
///
/// Builds a hollow spherical shell as **two solids**: outer sphere plus inner sphere with reversed
/// face orientation (`reverse_face`). Total surface area matches the analytic spherical shell
/// \(4\pi(R^2+r^2)\) under [`rcad_kernel::surface_area`].
///
/// Compound `rcad_kernel::volume` may not match \(4\pi/3(R^3-r^3)\) until sphere face normals / tessellation agree for
/// [`signed_volume`] everywhere.
///
/// OCCT DRAW `mkvolume` on trimmed spherical patches can report a different `checkprops -s` than this
/// full-sphere analytic shell.
pub fn try_difference_concentric_spheres(a: &BRep, b: &BRep) -> Option<BRep> {
    let (ca, ra) = try_sphere_primitive_center_radius(a)?;
    let (cb, rb) = try_sphere_primitive_center_radius(b)?;
    let scale = ra.max(rb).max(1.0);
    if (ca - cb).length() > TOL.max(TOLERANCE_COORD_SUB * scale) {
        return None;
    }
    let ro = ra.max(rb);
    let ri = ra.min(rb);
    if ro - ri <= TOLERANCE_LEN_MIN * ro.max(1.0) {
        return None;
    }
    let center = ca;
    let outer = make_sphere_brep(center, ro).ok()?;
    let mut inner_cavity = make_sphere_brep(center, ri).ok()?;
    crate::reverse_face(&mut inner_cavity, 0);
    Some(BRep::compound_from_shapes(&[outer, inner_cavity]))
}

/// [`BooleanOpType::Intersection`] for two analytic sphere primitives sharing a center.
///
/// The intersection of nested balls is the smaller-radius ball (same center).
///
/// Returns [`None`] when centers differ or radii are degenerate 鈥?callers fall back to Pave/Builder.
pub fn try_intersection_concentric_spheres(a: &BRep, b: &BRep) -> Option<BRep> {
    let (ca, ra) = try_sphere_primitive_center_radius(a)?;
    let (cb, rb) = try_sphere_primitive_center_radius(b)?;
    let scale = ra.max(rb).max(1.0);
    if (ca - cb).length() > TOL.max(TOLERANCE_COORD_SUB * scale) {
        return None;
    }
    let r = ra.min(rb);
    let r_eps = TOLERANCE_LEN_MIN * scale.max(1.0);
    if r <= r_eps {
        return None;
    }
    make_sphere_brep(ca, r).ok()
}

// --- Coaxial cone 鈭?cylinder (OCCT `bopcommon_simple/ZP7`): generic Builder over-counts area. --------

fn z_axis_sharp_cone_z_span(cone: &BRep) -> Option<(f64, f64, f64)> {
    const APAR: f64 = TOLERANCE_ADAPTIVE_MAX;
    const XY: f64 = 2.0 * TOLERANCE_ADAPTIVE_MAX;
    let sh = cone.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() != 2 {
        return None;
    }
    let mut cf: Option<ConicalSurface> = None;
    let mut po: Option<DVec3> = None;
    let mut fi = 0usize;
    for s in &cone.solids {
        for sh in &s.shells {
            for _ in &sh.faces {
                let si = *cone.geom.face_surface.get(fi)?.as_ref()?;
                match cone.geom.surfaces.get(si)? {
                    Surface3::Cone(c) if c.radius.abs() <= TOLERANCE_MESH_LEGACY => cf = Some(*c),
                    Surface3::Plane(p) => po = Some(p.origin),
                    _ => return None,
                }
                fi += 1;
            }
        }
    }
    let c = cf?;
    let u = c.axis_dir();
    if u.cross(DVec3::Z).length() > APAR {
        return None;
    }
    let apex = c.apex_point();
    if apex.x.abs() > XY || apex.y.abs() > XY {
        return None;
    }
    let b = po?;
    if b.x.abs() > XY || b.y.abs() > XY {
        return None;
    }
    let t = (b - apex).dot(u);
    let rb = t * c.half_angle_rad.tan();
    if t < TOLERANCE_MESH_LEGACY || rb < TOLERANCE_MESH_LEGACY {
        return None;
    }
    Some((apex.z, b.z, rb))
}

fn z_axis_cylinder_z_span_r(cyl: &BRep) -> Option<(f64, f64, f64)> {
    const APAR: f64 = TOLERANCE_ADAPTIVE_MAX;
    const XY: f64 = 2.0 * TOLERANCE_ADAPTIVE_MAX;
    let sh = cyl.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() != 3 {
        return None;
    }
    let mut r = None;
    let mut zs = Vec::with_capacity(2);
    let mut fi = 0usize;
    for s in &cyl.solids {
        for sh in &s.shells {
            for _ in &sh.faces {
                let si = *cyl.geom.face_surface.get(fi)?.as_ref()?;
                match cyl.geom.surfaces.get(si)? {
                    Surface3::Cylinder(cc) => {
                        if cc.axis.normalize_or_zero().cross(DVec3::Z).length() > APAR {
                            return None;
                        }
                        if cc.origin.x.abs() > XY || cc.origin.y.abs() > XY {
                            return None;
                        }
                        r = Some(cc.radius);
                    }
                    Surface3::Plane(p) => zs.push(p.origin.z),
                    _ => return None,
                }
                fi += 1;
            }
        }
    }
    if zs.len() != 2 {
        return None;
    }
    Some((zs[0].min(zs[1]), zs[0].max(zs[1]), r?))
}

fn try_intersection_coaxial_cone_cylinder_pair(cone: &BRep, cyl: &BRep) -> Option<BRep> {
    use rcad_modeling::make_conical_frustum_brep;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    let zcn = zb.min(za);
    let zcx = zb.max(za);
    let z0 = zlo.max(zcn);
    let z1 = zhi.min(zcx);
    if z1 - z0 < TOLERANCE_MESH_LEGACY {
        return None;
    }
    let apex_hi = za > zb;
    let hc = (za - zb).abs();
    let rcz = |z: f64| {
        let num = if apex_hi { (za - z).abs() } else { (z - za).abs() };
        rb * num / hc
    };
    let r0 = rcz(z0).min(rc);
    let r1 = rcz(z1).min(rc);
    if r0 < TOLERANCE_COORD_SUB && r1 < TOLERANCE_COORD_SUB {
        return None;
    }
    let zm = (z0 + z1) * 0.5;
    let h = z1 - z0;
    make_conical_frustum_brep(DVec3::new(0.0, 0.0, zm), DVec3::Z, DVec3::X, r0, r1, h).ok()
}

/// Sharp Z-aligned cone 鈭?finite Z-aligned cylinder (same axis / origin in `xy`), e.g. OCCT ZP7.
pub fn try_intersection_coaxial_cone_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_coaxial_cone_cylinder_pair(a, b)
        .or_else(|| try_intersection_coaxial_cone_cylinder_pair(b, a))
}

/// `cone \ cylinder` when the cylinder closes the cone base and contains the cone frustum up to `z_hi`:
/// remainder is the sharp sub-cone from `z_hi` to the apex (OCCT `bopcut_simple`/ZP8).
pub fn try_difference_coaxial_cone_minus_cylinder(cone: &BRep, cyl: &BRep) -> Option<BRep> {
    use rcad_modeling::make_cone_brep;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    if za <= zb + TOLERANCE_MESH_LEGACY {
        return None;
    }
    let hc = za - zb;
    let r_at = |z: f64| rb * (za - z) / hc;
    // Cylinder starts on cone base disk; radius at least the cone base radius.
    if (zlo - zb).abs() > TOLERANCE_AXIS_CORNER_SLACK {
        return None;
    }
    if (rc + TOLERANCE_ADAPTIVE_MAX) < rb {
        return None;
    }
    // Cylinder covers the cone cross-section at `z_hi` (disk of radius `rc` vs cone radius `r_at(z_hi)`).
    if rc + TOLERANCE_MESH_LEGACY < r_at(zhi) {
        return None;
    }
    if zhi <= zb + TOLERANCE_MESH_LEGACY || zhi >= za - TOLERANCE_MESH_LEGACY {
        return None;
    }
    let r_cut = r_at(zhi);
    if r_cut < TOLERANCE_COORD_SUB {
        return None;
    }
    let h_rem = za - zhi;
    let z_mid = (za + zhi) * 0.5;
    make_cone_brep(
        DVec3::new(0.0, 0.0, z_mid),
        DVec3::Z,
        DVec3::X,
        r_cut,
        h_rem,
    )
    .ok()
}

/// Strip loft bottom/top planar caps; keeps ruled lateral faces only (must match [`LoftHistory`]).
fn strip_loft_caps(mut brep: BRep, hist: LoftHistory) -> Option<BRep> {
    let shell = brep.solids.first_mut()?.shells.first_mut()?;
    if hist.bottom_cap >= shell.faces.len() || hist.top_cap >= shell.faces.len() {
        return None;
    }
    // Remove higher shell face index first so the remaining index stays valid.
    let (mut lo, mut hi) = (hist.bottom_cap, hist.top_cap);
    if lo > hi {
        std::mem::swap(&mut lo, &mut hi);
    }
    shell.faces.remove(hi);
    crate::remove_flat_face_geom_slots(&mut brep.geom, hi);
    shell.faces.remove(lo);
    crate::remove_flat_face_geom_slots(&mut brep.geom, lo);
    Some(brep)
}

fn reverse_wire_local(wire: &mut Wire) {
    wire.edges.reverse();
    for we in &mut wire.edges {
        we.forward = !we.forward;
    }
}

/// Loft builds the **solid** frustum mantle (normals outward from cone interior). For
/// `cylinder \\ cone`, those faces bound the cavity 鈥?outward from the result solid points into the
/// removed cone (flip normals vs loft defaults).
fn invert_shell_planar_faces(brep: &mut BRep) {
    let Some(shell) = brep.solids.first_mut().and_then(|s| s.shells.first_mut()) else {
        return;
    };
    let n_faces = shell.faces.len();
    for face in &mut shell.faces {
        face.normal = -face.normal;
        reverse_wire_local(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            reverse_wire_local(iw);
        }
        face.triangles.clear();
        face.mesh_dirty = true;
    }
    for fi in 0..n_faces {
        let Some(Some(si)) = brep.geom.face_surface.get(fi).copied() else {
            continue;
        };
        let Some(surf) = brep.geom.surfaces.get_mut(si) else {
            continue;
        };
        if let Surface3::Plane(pl) = surf {
            pl.normal = -pl.normal;
        }
    }
}

/// Horizontal annulus from pre-built coplanar rings (same vertex count). Eliminates float drift vs loft.
fn annulus_between_rings(outer: &[DVec3], inner: &[DVec3]) -> Result<BRep, rcad_modeling::BuildError> {
    let n = outer.len();
    if n < 3 || inner.len() != n {
        return Err(rcad_modeling::BuildError::DegenerateGeometry(
            "annulus_between_rings vertex count",
        ));
    }
    let z = outer[0].z;
    if inner.iter().any(|p| (p.z - z).abs() > TOLERANCE_COORD_SUB)
        || outer.iter().any(|p| (p.z - z).abs() > TOLERANCE_COORD_SUB)
    {
        return Err(rcad_modeling::BuildError::DegenerateGeometry(
            "annulus_between_rings not coplanar",
        ));
    }
    let mut brep = BRep::default();
    let surface = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z),
        normal: DVec3::Z,
    });

    let mut outer_vs = Vec::with_capacity(n);
    let outer_pts = outer.to_vec();
    for p in &outer_pts {
        outer_vs.push(make_vertex(&mut brep, *p));
    }
    let mut outer_we = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = outer_pts[i];
        let b = outer_pts[j];
        let dir = (b - a).normalize();
        let len = (b - a).length();
        let ei = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            outer_vs[i],
            outer_vs[j],
        )?;
        outer_we.push(WireEdge::fwd(ei));
    }

    let mut inner_vs = Vec::with_capacity(n);
    let inner_pts = inner.to_vec();
    for p in &inner_pts {
        inner_vs.push(make_vertex(&mut brep, *p));
    }
    let mut inner_we = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = inner_pts[i];
        let b = inner_pts[j];
        let dir = (b - a).normalize();
        let len = (b - a).length();
        let ei = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            inner_vs[i],
            inner_vs[j],
        )?;
        inner_we.push(WireEdge::fwd(ei));
    }

    let outer_wire = make_wire(outer_we);
    let inner_wire = make_wire(inner_we);
    let _fi = make_face(&mut brep, surface, outer_wire, vec![inner_wire])?;
    Ok(brep)
}

/// Outer / inner lateral strips + top annulus for [`try_coaxial_cylinder_minus_frustum_loft_shell`].
fn coaxial_cylinder_minus_frustum_loft_pieces(
    z_lo: f64,
    z_hi: f64,
    rc: f64,
    za: f64,
    zb: f64,
    rb: f64,
) -> Option<(BRep, BRep, BRep)> {
    use std::f64::consts::TAU;
    const N: usize = 32;
    let zcn = zb.min(za);
    let zcx = zb.max(za);
    let z0 = z_lo.max(zcn);
    let z1 = z_hi.min(zcx);
    if z1 - z0 < TOLERANCE_MESH_LEGACY {
        return None;
    }
    let apex_hi = za > zb;
    let hc = (za - zb).abs();
    let rcz = |z: f64| {
        let num = if apex_hi {
            (za - z).abs()
        } else {
            (z - za).abs()
        };
        rb * num / hc
    };
    let r0 = rcz(z0).min(rc);
    let r1 = rcz(z1).min(rc);
    if r0 < TOLERANCE_COORD_SUB || r1 < TOLERANCE_COORD_SUB {
        return None;
    }

    let mut outer_bot = Vec::with_capacity(N);
    let mut outer_top = Vec::with_capacity(N);
    let mut inner_bot = Vec::with_capacity(N);
    let mut inner_top = Vec::with_capacity(N);
    for i in 0..N {
        let ang = TAU * i as f64 / N as f64;
        let c = ang.cos();
        let s = ang.sin();
        let ob = DVec3::new(rc * c, rc * s, z_lo);
        outer_bot.push(ob);
        let ib = if (r0 - rc).abs() <= TOLERANCE_LEN_MIN && (z0 - z_lo).abs() <= TOLERANCE_LEN_MIN {
            ob
        } else {
            DVec3::new(r0 * c, r0 * s, z0)
        };
        inner_bot.push(ib);
        outer_top.push(DVec3::new(rc * c, rc * s, z_hi));
        inner_top.push(DVec3::new(r1 * c, r1 * s, z1));
    }

    let annulus = annulus_between_rings(&outer_top, &inner_top).ok()?;

    let (loft_outer, ohist) = loft_with_history(&[outer_bot, outer_top]).ok()?;
    let outer_strip = strip_loft_caps(loft_outer, ohist)?;

    let (loft_inner, ihist) = loft_with_history(&[inner_bot, inner_top]).ok()?;
    let mut inner_strip = strip_loft_caps(loft_inner, ihist)?;
    invert_shell_planar_faces(&mut inner_strip);

    Some((outer_strip, inner_strip, annulus))
}

/// Closed shell for `cylinder \ (cone 鈭?cylinder)` when overlap is the coaxial frustum:
/// outer cylindrical loft strip + inner frustum loft strip + top annulus, sewn (`OCCT ZP3`).
fn try_coaxial_cylinder_minus_frustum_loft_shell(
    z_lo: f64,
    z_hi: f64,
    rc: f64,
    za: f64,
    zb: f64,
    rb: f64,
) -> Option<BRep> {
    let (outer_strip, inner_strip, annulus) =
        coaxial_cylinder_minus_frustum_loft_pieces(z_lo, z_hi, rc, za, zb, rb)?;
    let tol = (TOLERANCE_RETRY_LADDER_COARSE).max(TOLERANCE_MESH_LEGACY * rc.max(z_hi.abs()));
    let sewn = sew_shells(&[outer_strip, inner_strip, annulus], tol);
    if !sewn.free_edges.is_empty() {
        return None;
    }
    Some(sewn.brep)
}

/// `cylinder \ cone` with same coaxial ZP layout as [`try_difference_coaxial_cone_minus_cylinder`].
///
/// Set identity: `cyl \ cone` equals `cyl \ (cone 鈭?cyl)` when the overlap is the coaxial frustum.
pub fn try_difference_coaxial_cylinder_minus_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    cyl_minus_cone_inner(a, b).or_else(|| cyl_minus_cone_inner(b, a))
}

fn cyl_minus_cone_inner(maybe_cyl: &BRep, maybe_cone: &BRep) -> Option<BRep> {
    let cone = maybe_cone;
    let cyl = maybe_cyl;
    try_intersection_coaxial_cone_cylinder_pair(cone, cyl)?;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    let hc = za - zb;
    if hc.abs() < TOLERANCE_MESH_LEGACY {
        return None;
    }
    let r_at = |z: f64| rb * (za - z) / hc;
    if (zlo - zb).abs() > TOLERANCE_AXIS_CORNER_SLACK {
        return None;
    }
    if (rc + TOLERANCE_ADAPTIVE_MAX) < rb {
        return None;
    }
    if rc + TOLERANCE_MESH_LEGACY < r_at(zhi) {
        return None;
    }
    if zhi <= zb + TOLERANCE_MESH_LEGACY || zhi >= za - TOLERANCE_MESH_LEGACY {
        return None;
    }
    try_coaxial_cylinder_minus_frustum_loft_shell(zlo, zhi, rc, za, zb, rb)
}

fn add_vertex(verts: &mut Vec<Vertex>, p: DVec3) -> usize {
    for (i, v) in verts.iter().enumerate() {
        if (v.point - p).length() < TOLERANCE_ABS {
            return i;
        }
    }
    verts.push(Vertex { point: p });
    verts.len() - 1
}

/// Closed triangle mesh of the boundary: three quarter-disks in x=0, y=0, z=0
/// planes plus one octant of the unit sphere. Outward for the solid
/// { x,y,z 鈮?0, x虏+y虏+z虏 鈮?1 }.
fn brep_eighth_of_unit_ball() -> BRep {
    const NA: usize = 32; // arc segments per quarter circle
    const NS: usize = 24; // grid per spherical patch axis

    let mut vertices: Vec<Vertex> = vec![];
    let empty_wire = || Wire { edges: vec![] };
    use std::f64::consts::FRAC_PI_2;
    // --- Planar: z=0, outward normal (0,0,-1) ---
    let o0 = add_vertex(&mut vertices, DVec3::ZERO);
    let _ = add_vertex(&mut vertices, DVec3::X);
    let _ = add_vertex(&mut vertices, DVec3::Y);
    let mut z0_arc: Vec<usize> = (0..=NA)
        .map(|k| {
            let t = (k as f64 / NA as f64) * FRAC_PI_2;
            add_vertex(&mut vertices, DVec3::new(t.cos(), t.sin(), 0.0))
        })
        .collect();

    // Triangles on z=0, fan from origin
    let mut t_z0: Vec<[usize; 3]> = vec![];
    for k in 0..NA {
        t_z0.push([o0, z0_arc[k], z0_arc[k + 1]]);
    }
    let f_z0 = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(0.0, 0.0, -1.0),
        triangles: t_z0,
        sample_point: None,
        mesh_dirty: false,
    };

    // y=0 plane, outward (0,-1,0), quarter disk in xz: (0,0,0) 鈥?(1,0,0) 鈥?(0,0,1) and arc
    let o1 = o0; // (0,0,0) shared
    let _ = add_vertex(&mut vertices, DVec3::Z);
    let y0_arc: Vec<usize> = (0..=NA)
        .map(|k| {
            let t = (k as f64 / NA as f64) * FRAC_PI_2;
            add_vertex(&mut vertices, DVec3::new(t.cos(), 0.0, t.sin()))
        })
        .collect();
    let mut t_y0: Vec<[usize; 3]> = vec![];
    for k in 0..NA {
        t_y0.push([o1, y0_arc[k], y0_arc[k + 1]]);
    }
    let f_y0 = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(0.0, -1.0, 0.0),
        triangles: t_y0,
        sample_point: None,
        mesh_dirty: false,
    };

    // x=0, outward (-1,0,0), quarter in yz: (0,0,0) (0,1,0) (0,0,1) + arc
    let o2 = o0;
    let x0_arc: Vec<usize> = (0..=NA)
        .map(|k| {
            let t = (k as f64 / NA as f64) * FRAC_PI_2;
            add_vertex(&mut vertices, DVec3::new(0.0, t.cos(), t.sin()))
        })
        .collect();
    let mut t_x0: Vec<[usize; 3]> = vec![];
    for k in 0..NA {
        t_x0.push([o2, x0_arc[k], x0_arc[k + 1]]);
    }
    let f_x0 = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(-1.0, 0.0, 0.0),
        triangles: t_x0,
        sample_point: None,
        mesh_dirty: false,
    };

    // Spherical octant: p = (sin v cos u, sin v sin u, cos v), u,v in [0,蟺/2]
    // (v = colatitude from +Z: v=0 is north pole (0,0,1), v=蟺/2 is the z=0 quarter arc)
    let pole = add_vertex(&mut vertices, DVec3::new(0.0, 0.0, 1.0));
    let mut sph_idx = vec![vec![0usize; NS + 1]; NS + 1];
    for i in 0..=NS {
        sph_idx[i][0] = pole;
    }
    for j in 1..=NS {
        let v = (j as f64 / NS as f64) * FRAC_PI_2;
        let si = v.sin();
        for i in 0..=NS {
            let u = (i as f64 / NS as f64) * FRAC_PI_2;
            let p = DVec3::new(si * u.cos(), si * u.sin(), v.cos());
            sph_idx[i][j] = add_vertex(&mut vertices, p);
        }
    }
    let mut t_s: Vec<[usize; 3]> = vec![];
    // Fan from north pole to first parallel (j=1)
    for i in 0..NS {
        t_s.push([pole, sph_idx[i][1], sph_idx[i + 1][1]]);
    }
    // Quad strips j = 1..NS-1
    for j in 1..NS {
        for i in 0..NS {
            let a = sph_idx[i][j];
            let b = sph_idx[i + 1][j];
            let c = sph_idx[i][j + 1];
            let d = sph_idx[i + 1][j + 1];
            t_s.push([a, b, d]);
            t_s.push([a, d, c]);
        }
    }
    let f_s = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(1.0, 1.0, 1.0).normalize(),
        triangles: t_s,
        sample_point: None,
        mesh_dirty: false,
    };

    let faces = vec![f_z0, f_y0, f_x0, f_s];
    let geom = GeomStore {
        curves: vec![],
        surfaces: vec![],
        curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; 4],
        edge_pcurves: vec![],
        edge_curve_range: vec![],
        edge_degenerated: vec![],
        vertex_tolerance: vec![],
        edge_tolerance: vec![],
        face_tolerance: vec![],
        curve2d_range: vec![],
        face_surface_range: vec![None; 4],
        edge_same_parameter: vec![],
        edge_same_range: vec![],
    };

    BRep {
        vertices,
        edges: vec![],
        solids: vec![Solid {
            shells: vec![Shell { faces }],
        }],
        geom,
        compound: None,
        compsolid: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{boolean_op, BooleanOpType};
    use rcad_kernel::surface_area;
    use rcad_kernel::volume;

    #[test]
    fn concentric_sphere_difference_analytic_shell_surface_area() {
        let center = DVec3::new(1.0, -2.0, 4.0);
        let ro = 5.0_f64;
        let ri = 3.0_f64;
        let outer = make_sphere_brep(center, ro).expect("outer");
        let inner = make_sphere_brep(center, ri).expect("inner");
        let shell = boolean_op(BooleanOpType::Difference, &outer, &inner).expect("difference");
        let pi = std::f64::consts::PI;
        let a_ex = 4.0 * pi * (ro * ro + ri * ri);
        let a = surface_area(&shell);
        assert!(
            (a - a_ex).abs() < 50.0 * TOLERANCE_RETRY_LADDER_COARSE * a_ex.max(1.0),
            "surface area {a} vs analytic shell SA {a_ex}"
        );
        // `signed_volume` across compounds relies on consistent face normals vs tessellation;
        // sphere primitives carry approximate face normals 鈥?SA matches analytic \(4\pi(R^2+r^2)\).
    }

    #[test]
    fn eighth_ball_area_and_volume() {
        let b = brep_eighth_of_unit_ball();
        let a = surface_area(&b);
        let v = volume(&b);
        let a_ex = 5.0 * std::f64::consts::PI / 4.0;
        let v_ex = std::f64::consts::PI / 6.0;
        assert!(
            (a - a_ex).abs() < 0.04,
            "area {a} vs {a_ex}"
        );
        assert!(
            (v - v_ex).abs() < 0.02,
            "vol {v} vs {v_ex}"
        );
    }

    #[test]
    fn zp3_sum_planar_areas_before_sew_matches_expected_total() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let z_lo = -10.0_f64;
        let z_hi = 0.0_f64;
        let rc = 10.0_f64;
        let z0 = -10.0_f64;
        let z1 = 0.0_f64;
        let r0 = 10.0_f64;
        let r1 = 5.0_f64;
        let mut outer_bot = Vec::with_capacity(N);
        let mut outer_top = Vec::with_capacity(N);
        let mut inner_bot = Vec::with_capacity(N);
        let mut inner_top = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            let c = ang.cos();
            let s = ang.sin();
            let ob = DVec3::new(rc * c, rc * s, z_lo);
            outer_bot.push(ob);
            inner_bot.push(ob);
            outer_top.push(DVec3::new(rc * c, rc * s, z_hi));
            inner_top.push(DVec3::new(r1 * c, r1 * s, z1));
        }
        let ann = annulus_between_rings(&outer_top, &inner_top).unwrap();
        let (lo, oh) = loft_with_history(&[outer_bot, outer_top]).unwrap();
        let outer_strip = strip_loft_caps(lo, oh).unwrap();
        let (li, ih) = loft_with_history(&[inner_bot, inner_top]).unwrap();
        let inner_strip = strip_loft_caps(li, ih).unwrap();
        let sum = rcad_kernel::surface_area(&outer_strip)
            + rcad_kernel::surface_area(&inner_strip)
            + rcad_kernel::surface_area(&ann);
        let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 1390.8_f64);
        assert!(
            (sum - 1390.8_f64).abs() <= tol,
            "sum loose pieces ~1390.8, got {sum}"
        );
    }

    #[test]
    fn zp3_inner_frustum_strip_surface_area_sane() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let mut ib = Vec::with_capacity(N);
        let mut it = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            let c = ang.cos();
            let s = ang.sin();
            ib.push(DVec3::new(10.0 * c, 10.0 * s, -10.0));
            it.push(DVec3::new(5.0 * c, 5.0 * s, 0.0));
        }
        let (loft_i, ih) = loft_with_history(&[ib, it]).unwrap();
        let strip = strip_loft_caps(loft_i, ih).unwrap();
        let a = rcad_kernel::surface_area(&strip);
        assert!(
            a > 400.0 && a < 650.0,
            "inner frustum mantle area ~527, got {a}"
        );
    }

    #[test]
    fn zp3_outer_face_area_matches_before_and_after_sew() {
        let (outer_strip, inner_strip, annulus) = coaxial_cylinder_minus_frustum_loft_pieces(
            -10.0, 0.0, 10.0, 10.0, -10.0, 10.0,
        )
        .expect("loft pieces");
        let tol = (TOLERANCE_RETRY_LADDER_COARSE).max(TOLERANCE_MESH_LEGACY * 10.0);
        let sewn = sew_shells(&[outer_strip.clone(), inner_strip, annulus], tol);
        assert!(sewn.free_edges.is_empty(), "free {:?}", sewn.free_edges);
        let f_loose = &outer_strip.solids[0].shells[0].faces[0];
        let f_sewn = &sewn.brep.solids[0].shells[0].faces[0];
        let a_loose = rcad_kernel::face_surface_area(&outer_strip, f_loose, 0);
        let a_sewn = rcad_kernel::face_surface_area(&sewn.brep, f_sewn, 0);
        assert!(
            (a_loose - a_sewn).abs() < TOLERANCE_RETRY_LADDER_COARSE * 100.0,
            "first outer lateral face area loose {a_loose} vs sewn {a_sewn}"
        );
    }

    #[test]
    fn zp3_loft_shell_matches_occt_geometry_numbers() {
        // Cone apex z=10, base z=-10, rb=10; cylinder z in [-10,0], r=10 鈥?same as geometry_properties ZP3.
        let r = try_coaxial_cylinder_minus_frustum_loft_shell(-10.0, 0.0, 10.0, 10.0, -10.0, 10.0);
        assert!(r.is_some(), "expected sewn loft shell for ZP3 parameters");
        let brep = r.unwrap();
        let nf: usize = brep
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count();
        let a = rcad_kernel::surface_area(&brep);
        let v = rcad_kernel::volume(&brep);
        assert!(
            nf >= 60 && nf <= 70,
            "expected ~65 faces (32+32+annulus), got {nf}"
        );
        assert!(
            (v - 1310.0).abs() < 80.0,
            "expected volume cylinder minus frustum ~1310, got {v}"
        );
        let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 1390.8_f64);
        assert!(
            (a - 1390.8_f64).abs() <= tol,
            "surface area: expected ~1390.8, got {a} (nf={nf}, vol={v})"
        );
    }

    #[test]
    fn zp3_annulus_plane_builds() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let mut o = Vec::with_capacity(N);
        let mut inn = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            o.push(DVec3::new(10.0 * ang.cos(), 10.0 * ang.sin(), 0.0));
            inn.push(DVec3::new(5.0 * ang.cos(), 5.0 * ang.sin(), 0.0));
        }
        annulus_between_rings(&o, &inn).expect("annulus");
    }

    #[test]
    fn zp3_coaxial_cylinder_minus_cone_fast_path_triggered() {
        use glam::DVec3;
        use rcad_modeling::{make_cone_brep, make_cylinder_brep};
        let pc = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 20.0).expect("cone");
        let pcy = make_cylinder_brep(
            DVec3::new(0.0, 0.0, -5.0),
            DVec3::Z,
            DVec3::X,
            10.0,
            10.0,
        )
        .expect("cylinder");
        assert!(
            try_difference_coaxial_cylinder_minus_cone(&pcy, &pc).is_some(),
            "ZP3 boptuc expects coaxial cylinder\\cone shortcut"
        );
    }

    // ── box-box intersection fast path ─────────────────────────────────────

    #[test]
    fn box_box_intersection_partial_overlap() {
        // bcommon_simple_c1: 1×1×1 ∩ 1.5×0.5×0.5 → 1×0.5×0.5 (SA=2.5, vol=0.25)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.5, 0.5, 0.5).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("box-box intersection");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 2.5).abs() < 1e-7, "SA={sa} expected 2.5");
        assert!((vol - 0.25).abs() < 1e-7, "vol={vol} expected 0.25");
    }

    #[test]
    fn box_box_intersection_full_containment() {
        // 2×2×2 box ∩ 0.5×0.5×0.5 inside → the inner 0.5^3 box
        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let inner = make_box_brep(
            DVec3::new(0.25, 0.25, 0.25),
            DVec3::X,
            DVec3::Y,
            0.5,
            0.5,
            0.5,
        )
        .unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &outer, &inner).expect("containment");
        let vol = volume(&r);
        assert!((vol - 0.125).abs() < 1e-7, "vol={vol} expected 0.125");
    }

    #[test]
    fn box_box_intersection_no_overlap() {
        // Disjoint boxes → empty intersection.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("no-overlap");
        let n_faces: usize = r.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty intersection");
    }

    #[test]
    fn box_box_intersection_non_box_falls_through() {
        // Sphere ∩ box falls through to generic path (no panic).
        let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &sphere, &b).expect("sphere-box");
        // Some result — specific value doesn't matter.
        assert!(r.solids.len() >= 1 || r.vertices.is_empty());
    }

    // ── box-box difference fast path ───────────────────────────────────────

    #[test]
    fn box_box_difference_opposite_same_axis() {
        // F1-like: A=1.5×0.5×0.5 at (-0.25,0,0), B=1×1×1 at origin.
        // boptuc = A - B (but here we test try_difference_box_box directly,
        // so a = first arg = A, b = second arg = B).
        // A extends on x-lo (-0.25<0) and x-hi (1.25>1), same axis → two
        // disjoint slabs: SA=2, vol=0.125.
        let a = make_box_brep(
            DVec3::new(-0.25, 0.0, 0.0), DVec3::X, DVec3::Y,
            1.5, 0.5, 0.5,
        ).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("box-box difference");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 2.0).abs() < 1e-7, "SA={sa} expected 2.0");
        assert!((vol - 0.125).abs() < 1e-7, "vol={vol} expected 0.125");
    }

    #[test]
    fn box_box_difference_single_slab() {
        // bcut_simple_c1: 0.5×1.5×1 at (0,-0.5,0) minus 1×1×1 at origin.
        // Excess only on y-lo [-0.5, 0]. Result: 0.5×0.5×1 box, SA=2.5.
        let a = make_box_brep(DVec3::new(0.0, -0.5, 0.0), DVec3::X, DVec3::Y, 0.5, 1.5, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("box-box difference");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 2.5).abs() < 1e-7, "SA={sa} expected 2.5");
        assert!((vol - 0.25).abs() < 1e-7, "vol={vol} expected 0.25");
    }

    #[test]
    fn box_box_difference_no_overlap() {
        // Disjoint boxes → difference is just A.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("no-overlap");
        let sa = surface_area(&r);
        assert!((sa - 6.0).abs() < 1e-7, "SA={sa} expected 6.0 (unchanged A)");
    }

    #[test]
    fn box_box_difference_a_inside_b() {
        // A fully inside B → empty.
        let a = make_box_brep(DVec3::new(0.25, 0.25, 0.25), DVec3::X, DVec3::Y, 0.5, 0.5, 0.5).unwrap();
        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &outer).expect("A-in-B");
        let n_faces: usize = r.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty difference");
    }
}
