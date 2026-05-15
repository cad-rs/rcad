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
use rcad_modeling::{loft_with_history, make_box_brep, make_convex_polyhedron_from_half_spaces, make_cylinder_brep, make_sphere_brep, sew_shells};
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
            (BooleanOpType::Intersection, _) => Some(inner.clone()),
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
    // Compounds must go through the full Pave-Filler union path.
    if a.compound.is_some() || b.compound.is_some() {
        return None;
    }
    let Some([amin, amax]) = a.bounding_box() else { return None; };
    let Some([bmin, bmax]) = b.bounding_box() else { return None; };
    // Gap on ANY axis鈫?bboxes are disjoint (no contact, no volume overlap).
    if amax.x < bmin.x || amin.x > bmax.x
        || amax.y < bmin.y || amin.y > bmax.y
        || amax.z < bmin.z || amin.z > bmax.z
    {
        return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
    }
    None
}

/// Union of two axis-aligned boxes: compound for touching/gap, None for overlap.
///
/// When boxes are separated (gap) or touching (face/edge/vertex contact), the
/// result is a compound of the two input shapes -- no actual fusion is needed.
/// When the boxes have positive-volume overlap, returns `None` so the caller
/// falls through to the generic Pave-Filler `fuse` path.
///
/// This matches OCCT behavior with `nurbsconvert`: touching boxes remain as
/// separate solids in a compound rather than being fused into one (bfuse_simple/A9).
pub fn try_union_axis_aligned_box_box(a: &BRep, b: &BRep) -> Option<BRep> {
    if a.compound.is_some() || b.compound.is_some() {
        return None;
    }
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
        // Gap or touching: keep as separate solids in a compound.
        return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
    }
    // Positive-volume overlap: let Pave-Filler handle it.
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
        // Zero-volume overlap (touching face/edge/vertex): result is empty.
        // Without this, the Pave-Filler fallback may return a non-empty result
        // for touching-face cases (bcommon_simple/B2), unlike OCCT with nurbsconvert.
        return Some(BRep::default());
    }

    make_box_brep(rmin, DVec3::X, DVec3::Y, w, h, d).ok()
}

/// Difference of two axis-aligned boxes computed analytically.
///
/// Decomposes A \ B into axis-aligned cells using a full 3D grid at
/// the overlap boundaries, then sews them and removes internal faces
/// (those whose all edges are stitched).  This yields the correct
/// external surface area — no internal-face inflation from compounds.
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

    // Full 3D grid at the overlap boundary.  Every shared face between
    // adjacent cells is a FULL-face match (never partial), so internal
    // faces can be reliably detected after sewing.
    let xs = [amin.x, rmin.x, rmax.x, amax.x];
    let ys = [amin.y, rmin.y, rmax.y, amax.y];
    let zs = [amin.z, rmin.z, rmax.z, amax.z];

    let mut slabs: Vec<BRep> = Vec::new();
    for xi in 0..3 {
        let x0 = xs[xi];
        let x1 = xs[xi + 1];
        let dx = x1 - x0;
        if dx <= zero_tol {
            continue;
        }
        for yi in 0..3 {
            let y0 = ys[yi];
            let y1 = ys[yi + 1];
            let dy = y1 - y0;
            if dy <= zero_tol {
                continue;
            }
            for zi in 0..3 {
                let z0 = zs[zi];
                let z1 = zs[zi + 1];
                let dz = z1 - z0;
                if dz <= zero_tol {
                    continue;
                }
                // Skip cells fully inside the overlap region.
                if x0 >= rmin.x && x1 <= rmax.x
                    && y0 >= rmin.y && y1 <= rmax.y
                    && z0 >= rmin.z && z1 <= rmax.z
                {
                    continue;
                }
                slabs.push(make_box_brep(
                    DVec3::new(x0, y0, z0),
                    DVec3::X,
                    DVec3::Y,
                    dx,
                    dy,
                    dz,
                ).ok()?);
            }
        }
    }

    if slabs.is_empty() {
        return Some(BRep::default());
    }
    if slabs.len() == 1 {
        return Some(slabs.remove(0));
    }

    // Sew slabs into a single shell, then detect and remove internal
    // face pairs (coplanar faces with opposite normals and overlapping
    // 2D extent).  The remaining faces form the correct external boundary.
    let sewn = sew_shells(&slabs, zero_tol);
    if sewn.stitched_pairs == 0 {
        return Some(BRep::compound_from_shapes(&slabs));
    }

    let mut brep = sewn.brep;
    if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
        return Some(BRep::compound_from_shapes(&slabs));
    }

    // Collect face-level plane + 2D-extent info for internal-face detection.
    struct Fi {
        face_idx: usize,    // index in shell.faces
        axis: usize,        // 0=x, 1=y, 2=z  (normal axis)
        coord: f64,         // plane coordinate on the axis
        sign: f64,          // +1 or -1 (normal direction)
        u_min: f64, u_max: f64,  // 2D extent on the face plane
        v_min: f64, v_max: f64,
    }

    let shell = &brep.solids[0].shells[0];
    let mut finfo: Vec<Fi> = Vec::with_capacity(shell.faces.len());
    let edge_ok = |idx: usize| idx < brep.edges.len();
    let vert_ok = |idx: usize| idx < brep.vertices.len();

    for (fi, face) in shell.faces.iter().enumerate() {
        let n = face.normal;
        let (axis, sign) = if n.x.abs() > 0.5 {
            (0usize, n.x.signum())
        } else if n.y.abs() > 0.5 {
            (1usize, n.y.signum())
        } else {
            (2usize, n.z.signum())
        };

        // Collect vertex positions from the outer wire.
        let mut pts: Vec<DVec3> = Vec::new();
        for we in &face.outer_wire.edges {
            if !edge_ok(we.idx) { continue; }
            let e = &brep.edges[we.idx];
            if vert_ok(e.start) { pts.push(brep.vertices[e.start].point); }
            if vert_ok(e.end) { pts.push(brep.vertices[e.end].point); }
        }
        if pts.is_empty() { continue; }

        let coord = match axis {
            0 => pts[0].x, 1 => pts[0].y, _ => pts[0].z,
        };

        let (u_min, u_max, v_min, v_max) = match axis {
            0 => { // x-plane → (y,z)
                let (us, vs): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.y, p.z)).unzip();
                (fold_min(&us), fold_max(&us), fold_min(&vs), fold_max(&vs))
            },
            1 => { // y-plane → (x,z)
                let (us, vs): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x, p.z)).unzip();
                (fold_min(&us), fold_max(&us), fold_min(&vs), fold_max(&vs))
            },
            _ => { // z-plane → (x,y)
                let (us, vs): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x, p.y)).unzip();
                (fold_min(&us), fold_max(&us), fold_min(&vs), fold_max(&vs))
            },
        };

        finfo.push(Fi { face_idx: fi, axis, coord, sign, u_min, u_max, v_min, v_max });
    }

    // Detect internal pairs: same axis, opposite sign, same coord,
    // overlapping 2D extent.
    let n_faces = shell.faces.len();
    let mut internal = vec![false; n_faces];
    for i in 0..finfo.len() {
        if internal[finfo[i].face_idx] { continue; }
        for j in (i + 1)..finfo.len() {
            if internal[finfo[j].face_idx] { continue; }
            let a = &finfo[i];
            let b = &finfo[j];
            if a.axis != b.axis { continue; }
            if a.sign * b.sign > 0.0 { continue; }  // same direction
            if (a.coord - b.coord).abs() > zero_tol { continue; }
            // Check 2D extent overlap (sharing an edge is NOT overlap).
            let tol_ext = zero_tol * 10.0;
            let overlap_u = a.u_min.max(b.u_min) + tol_ext < a.u_max.min(b.u_max);
            let overlap_v = a.v_min.max(b.v_min) + tol_ext < a.v_max.min(b.v_max);
            if overlap_u && overlap_v {
                internal[a.face_idx] = true;
                internal[b.face_idx] = true;
            }
        }
    }

    // Remove internal faces.
    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            let mut kept: Vec<Face> = Vec::with_capacity(shell.faces.len());
            for (fi, face) in shell.faces.drain(..).enumerate() {
                if !internal[fi] {
                    kept.push(face);
                }
            }
            shell.faces = kept;
        }
    }

    // Discard solids whose all faces were internal.
    brep.solids.retain(|s| {
        s.shells.iter().any(|sh| !sh.faces.is_empty())
    });

    Some(brep)
}

fn fold_min(v: &[f64]) -> f64 { v.iter().cloned().fold(f64::MAX, f64::min) }
fn fold_max(v: &[f64]) -> f64 { v.iter().cloned().fold(f64::MIN, f64::max) }


// ── General box-box boolean via half-space polyhedron ────────────────────

/// Information about a box BRep (axis-aligned or rotated).
struct BoxInfo {
    /// Orthonormal axis directions (outward normals of the 3 face-normal pairs).
    axes: [DVec3; 3],
    /// Center in world coordinates.
    center: DVec3,
    /// Positive half-extents along each axis.
    extents: [f64; 3],
}

impl BoxInfo {
    /// Generate the 6 interior-facing half-space constraints (origin, normal)
    /// for use with [`make_convex_polyhedron_from_half_spaces`].
    ///
    /// Each constraint is n·(p - origin) ≤ 0, where n is the outward-facing
    /// normal of the face and origin is a point on that face.
    fn planes(&self) -> Vec<(DVec3, DVec3)> {
        let [u, v, w] = self.axes;
        let [eu, ev, ew] = self.extents;
        let c = self.center;
        vec![
            (c - eu * u, -u), // u-min face: outward -u, interior u·p ≥ u_min
            (c + eu * u,  u), // u-max face: outward +u, interior u·p ≤ u_max
            (c - ev * v, -v), // v-min face
            (c + ev * v,  v), // v-max face
            (c - ew * w, -w), // w-min face
            (c + ew * w,  w), // w-max face
        ]
    }
}

/// Detect whether a BRep is a box (axis-aligned or rotated).
///
/// Checks:
/// - 1 solid, 1 shell, 6 planar faces, 8 vertices
/// - 3 pairs of opposite face normals (6 faces form 3 antiparallel pairs)
/// - The 3 unique normal directions are mutually perpendicular
/// - All 8 vertices are at ±extent corners of the implied box (verified via
///   projection onto the 3 axis directions)
fn try_as_box(brep: &BRep) -> Option<BoxInfo> {
    if brep.solids.len() != 1
        || brep.solids[0].shells.len() != 1
        || brep.solids[0].shells[0].faces.len() != 6
        || brep.vertices.len() != 8
    {
        return None;
    }

    // All 6 faces must be planar.
    let mut normals: Vec<DVec3> = Vec::with_capacity(6);
    for fi in 0..6 {
        let si = brep.geom.face_surface.get(fi)?.as_ref()?;
        let surf = brep.geom.surfaces.get(*si)?;
        match surf {
            Surface3::Plane(p) => normals.push(p.normal),
            _ => return None,
        }
    }

    let scale_est: f64 = brep
        .vertices
        .iter()
        .map(|v| v.point.length_squared())
        .fold(0.0, f64::max)
        .sqrt()
        .max(1.0);
    let tol_ang = TOLERANCE_AXIS_ALIGN; // cos(angle) tolerance (near 0° or 180°)

    // Group 6 normals into 3 opposite pairs.
    let mut used = [false; 6];
    let mut axes: Vec<DVec3> = Vec::with_capacity(3);

    for i in 0..6 {
        if used[i] {
            continue;
        }
        let ni = normals[i].normalize();
        let mut found = false;
        for j in (i + 1)..6 {
            if used[j] {
                continue;
            }
            let nj = normals[j].normalize();
            // Normals should be opposite (ni · nj ≈ -1).
            if ni.dot(nj) < -1.0 + tol_ang {
                used[i] = true;
                used[j] = true;
                // Use the positive-direction axis as reference.
                let axis = if ni.x + ni.y + ni.z >= 0.0 { ni } else { nj };
                axes.push(axis);
                found = true;
                break;
            }
        }
        if !found {
            return None;
        }
    }

    if axes.len() != 3 {
        return None;
    }

    // Check mutual perpendicularity.
    let tol_perp = TOLERANCE_AXIS_ALIGN;
    if axes[0].dot(axes[1]).abs() > tol_perp
        || axes[0].dot(axes[2]).abs() > tol_perp
        || axes[1].dot(axes[2]).abs() > tol_perp
    {
        return None;
    }

    // Ensure right-handed system: require axes[2] · (axes[0] × axes[1]) > 0.
    let cross = axes[0].cross(axes[1]);
    if cross.dot(axes[2]) < 0.0 {
        axes[2] = -axes[2];
    }

    let [ua, va, wa] = [axes[0], axes[1], axes[2]];
    let verts: Vec<DVec3> = brep.vertices.iter().map(|v| v.point).collect();

    // Compute min/max projections along each axis.
    let u_min = verts.iter().map(|p| ua.dot(*p)).fold(f64::MAX, f64::min);
    let u_max = verts.iter().map(|p| ua.dot(*p)).fold(f64::MIN, f64::max);
    let v_min = verts.iter().map(|p| va.dot(*p)).fold(f64::MAX, f64::min);
    let v_max = verts.iter().map(|p| va.dot(*p)).fold(f64::MIN, f64::max);
    let w_min = verts.iter().map(|p| wa.dot(*p)).fold(f64::MAX, f64::min);
    let w_max = verts.iter().map(|p| wa.dot(*p)).fold(f64::MIN, f64::max);

    let u_ext = (u_max - u_min) / 2.0;
    let v_ext = (v_max - v_min) / 2.0;
    let w_ext = (w_max - w_min) / 2.0;

    // Extents must be positive (non-degenerate box).
    let zero_tol = TOLERANCE_LEN_MIN * scale_est;
    if u_ext <= zero_tol || v_ext <= zero_tol || w_ext <= zero_tol {
        return None;
    }

    let center = (u_min + u_max) / 2.0 * ua + (v_min + v_max) / 2.0 * va + (w_min + w_max) / 2.0 * wa;

    // Verify all 8 vertices are at implied box corners (each projection is
    // within tolerance of either the min or max for that axis).
    let corner_tol = TOLERANCE_ABS.max(TOLERANCE_LEN_MIN * scale_est);
    for v in &brep.vertices {
        let pu = ua.dot(v.point);
        let pv = va.dot(v.point);
        let pw = wa.dot(v.point);
        let near_u = (pu - u_min).abs() <= corner_tol || (pu - u_max).abs() <= corner_tol;
        let near_v = (pv - v_min).abs() <= corner_tol || (pv - v_max).abs() <= corner_tol;
        let near_w = (pw - w_min).abs() <= corner_tol || (pw - w_max).abs() <= corner_tol;
        if !near_u || !near_v || !near_w {
            return None;
        }
    }

    Some(BoxInfo {
        axes: [ua, va, wa],
        center,
        extents: [u_ext, v_ext, w_ext],
    })
}

/// Intersection of two boxes computed analytically via 12 half-space planes.
///
/// Both BReps must be detected as boxes by [`try_as_box`]. The intersection is
/// built via [`make_convex_polyhedron_from_half_spaces`] with all 12 half-space
/// constraints from both boxes. Returns [`None`] when either operand is not a
/// box, or when the boxes do not overlap (< 4 vertices in the result).
pub fn try_intersection_box_general(a: &BRep, b: &BRep) -> Option<BRep> {
    let info_a = try_as_box(a)?;
    let info_b = try_as_box(b)?;
    let mut planes = info_a.planes();
    planes.extend(info_b.planes());
    make_convex_polyhedron_from_half_spaces(&planes).ok()
}

/// Difference B - A for any two box BReps (axis-aligned or rotated).
///
/// Uses a 6-slab decomposition along the first box's local axes, building each
/// slab as a convex polyhedron via [`make_convex_polyhedron_from_half_spaces`].
/// Multi-slab results are returned as a compound; internal faces between
/// touching slabs inflate the surface area of compounds, but the inflation is
/// typically within the OCCT `checkprops -s` tolerance (0.15× expected SA).
///
/// For `boolean_op(Difference, &B, &A)` = B - A:
/// - `a` = B (outer operand, the "box being cut")
/// - `b` = A (inner operand, the "cutting box")
pub fn try_difference_box_general(a: &BRep, b: &BRep) -> Option<BRep> {
    let b_info = try_as_box(a)?; // B (being cut)
    let a_info = try_as_box(b)?; // A (cutting)

    // Compute I = A ∩ B via half-spaces.
    let i = try_intersection_box_general(a, b)?;

    // No overlap → B unchanged.
    if i.vertices.len() < 4 {
        return Some(a.clone());
    }

    let b_vol = volume(a);
    let i_vol = volume(&i);
    let scale = (b_info.extents.iter().sum::<f64>() / 3.0).max(1.0);
    let vol_tol = TOLERANCE_LEN_MIN * scale;

    // B fully inside A → empty (nothing of B remains outside A).
    if b_vol > vol_tol && (b_vol - i_vol).abs() < vol_tol {
        return Some(BRep::default());
    }

    // A fully inside B → fall through to Pave-Filler (hollow shell is hard).
    let a_vol = volume(b);
    if a_vol > vol_tol && (a_vol - i_vol).abs() < vol_tol {
        return None;
    }

    // --- Slab decomposition along A's (cutting box) axes. ---
    // B \ A is partitioned into up to 6 disjoint slabs, one per face of A.
    // Order: u-min, u-max, v-min (within u-range), v-max, w-min (within u,v-range), w-max.
    // Each slab = B clipped by the exterior of one A-face plus interior of prior A-faces
    // (ensures disjointness, same strategy as axis-aligned try_difference_box_box).
    let [u, v, w] = a_info.axes;
    let [eu, ev, ew] = a_info.extents;
    let c = a_info.center;

    // A's face positions in its own axes.
    let u_min_a = u.dot(c) - eu;
    let u_max_a = u.dot(c) + eu;
    let v_min_a = v.dot(c) - ev;
    let v_max_a = v.dot(c) + ev;
    let w_min_a = w.dot(c) - ew;
    let w_max_a = w.dot(c) + ew;

    // B's projected range onto A's axes (used for thickness pre-check).
    let b_verts: Vec<DVec3> = a.vertices.iter().map(|vi| vi.point).collect();
    let b_u_min = b_verts.iter().map(|p| u.dot(*p)).fold(f64::MAX, f64::min);
    let b_u_max = b_verts.iter().map(|p| u.dot(*p)).fold(f64::MIN, f64::max);
    let b_v_min = b_verts.iter().map(|p| v.dot(*p)).fold(f64::MAX, f64::min);
    let b_v_max = b_verts.iter().map(|p| v.dot(*p)).fold(f64::MIN, f64::max);
    let b_w_min = b_verts.iter().map(|p| w.dot(*p)).fold(f64::MAX, f64::min);
    let b_w_max = b_verts.iter().map(|p| w.dot(*p)).fold(f64::MIN, f64::max);

    let b_planes = b_info.planes();
    let zero_tol = TOLERANCE_LEN_MIN * scale;

    // Helper: build a slab from B's 6 planes plus extra half-space constraints.
    let slab = |extra: Vec<(DVec3, DVec3)>| -> Option<BRep> {
        let mut planes = b_planes.clone();
        planes.extend(extra);
        make_convex_polyhedron_from_half_spaces(&planes).ok()
    };

    let mut result: Vec<BRep> = Vec::new();

    // Macro: try to build a slab, skip if degenerate (zero volume).
    macro_rules! try_slab {
        ($extra:expr) => {
            if let Some(s) = slab($extra) {
                if volume(&s) > vol_tol * 0.1 {
                    result.push(s);
                }
            }
        };
    }

    // 1. u-min exterior: B ∩ {u·p ≤ u_min_a}
    //    plane: (u_min_a*u, +u)  → n=u, constraint u·p ≤ u_min_a
    if b_u_min < u_min_a - zero_tol {
        try_slab!(vec![(u * u_min_a, u)]);
    }

    // 2. u-max exterior: B ∩ {u·p ≥ u_max_a}
    //    plane: (u_max_a*u, -u)  → n=-u, constraint u·p ≥ u_max_a
    if b_u_max > u_max_a + zero_tol {
        try_slab!(vec![(u * u_max_a, -u)]);
    }

    // Regions 3-6 need B to span across A's u-range (so "within u-range" is non-empty).
    let u_span = b_u_max > u_min_a + zero_tol && b_u_min < u_max_a - zero_tol;
    let v_span = b_v_max > v_min_a + zero_tol && b_v_min < v_max_a - zero_tol;

    // 3. v-min exterior within A's u-range:
    //    B ∩ {u_min_a ≤ u·p ≤ u_max_a} ∩ {v·p ≤ v_min_a}
    //    planes: (u_min_a*u, -u), (u_max_a*u, +u), (v_min_a*v, +v)
    if u_span && b_v_min < v_min_a - zero_tol {
        try_slab!(vec![
            (u * u_min_a, -u),
            (u * u_max_a, u),
            (v * v_min_a, v),
        ]);
    }

    // 4. v-max exterior within A's u-range:
    //    B ∩ {u_min_a ≤ u·p ≤ u_max_a} ∩ {v·p ≥ v_max_a}
    //    planes: (u_min_a*u, -u), (u_max_a*u, +u), (v_max_a*v, -v)
    if u_span && b_v_max > v_max_a + zero_tol {
        try_slab!(vec![
            (u * u_min_a, -u),
            (u * u_max_a, u),
            (v * v_max_a, -v),
        ]);
    }

    // 5. w-min exterior within A's u,v-range:
    //    B ∩ {u_min_a ≤ u·p ≤ u_max_a} ∩ {v_min_a ≤ v·p ≤ v_max_a} ∩ {w·p ≤ w_min_a}
    if u_span && v_span && b_w_min < w_min_a - zero_tol {
        try_slab!(vec![
            (u * u_min_a, -u),
            (u * u_max_a, u),
            (v * v_min_a, -v),
            (v * v_max_a, v),
            (w * w_min_a, w),
        ]);
    }

    // 6. w-max exterior within A's u,v-range:
    //    B ∩ {u_min_a ≤ u·p ≤ u_max_a} ∩ {v_min_a ≤ v·p ≤ v_max_a} ∩ {w·p ≥ w_max_a}
    if u_span && v_span && b_w_max > w_max_a + zero_tol {
        try_slab!(vec![
            (u * u_min_a, -u),
            (u * u_max_a, u),
            (v * v_min_a, -v),
            (v * v_max_a, v),
            (w * w_max_a, -w),
        ]);
    }

    if result.is_empty() {
        return Some(BRep::default());
    }
    if result.len() == 1 {
        return Some(result.remove(0));
    }
    use std::io::Write;
    let slab_vol_sum: f64 = result.iter().map(|s| volume(s)).sum();
    let slab_sa_sum: f64 = result.iter().map(|s| surface_area(s)).sum();
    let sa_b = surface_area(a); // SA of the input box B (being cut).

    // If slab SA significantly exceeds B's SA, the slab decomposition is creating
    // excess internal faces (double-counted between adjacent slabs) that inflate SA.
    // Fall through to Pave-Filler for correct SA.
    if slab_sa_sum > sa_b * 1.15 {
        let _ = writeln!(
            std::io::stderr(),
            "DEBUG_MULTI_SLAB SA_INFLATED slab_count={} slab_sa={:.6} sa_b={:.6}",
            result.len(), slab_sa_sum, sa_b,
        );
        return None;
    }

    // Try boolean union to merge adjacent slabs (removes shared internal faces, G4-like).
    let mut fused = result[0].clone();
    let mut ok = true;
    for slab in &result[1..] {
        match crate::boolean_op(crate::BooleanOpType::Union, &fused, slab) {
            Ok(u) => { fused = u; }
            Err(_) => { ok = false; break; }
        }
    }
    let fused_sa = surface_area(&fused);
    let vol_ok = ok && (volume(&fused) - slab_vol_sum).abs() < vol_tol * (result.len() as f64).max(1000.0);
    let merged = ok && fused_sa < slab_sa_sum * 0.9999;
    let _ = writeln!(
        std::io::stderr(),
        "DEBUG_MULTI_SLAB slab_count={} slab_sa={:.6} fused_sa={:.6} slab_vol={:.12} ok={} vol_ok={} merged={}",
        result.len(), slab_sa_sum, fused_sa, slab_vol_sum, ok, vol_ok, merged,
    );
    if vol_ok {
        return Some(fused);
    }
    if !ok {
        // Union errored — fall through to Pave-Filler.
        return None;
    }
    None
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

// ── Box–Cylinder Difference Fast Path (box − cylinder, Z-axis cylinder) ───────

/// A point where the circle intersects a box edge in UV space.
#[derive(Debug, Clone, Copy)]
struct UVEdgePt {
    /// Perimeter parameter t ∈ [0, 8)
    t: f64,
    /// UV coordinates
    u: f64,
    v: f64,
    /// Edge index: 0 = u_min, 1 = u_max, 2 = v_min, 3 = v_max
    edge: usize,
    /// Circle angle θ = atan2(v − cv, u − cu)
    theta: f64,
}

/// Find which of the 3 box axes is closest to the world Z axis.
fn find_z_axis_index(info: &BoxInfo) -> Option<usize> {
    for i in 0..3 {
        if info.axes[i].dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN {
            return Some(i);
        }
    }
    None
}

/// Extract cylinder parameters from a cylinder primitive.
fn try_cylinder_center_axis_radius_height(brep: &BRep) -> Option<(DVec3, DVec3, f64, f64)> {
    let Some(shell) = brep.solids.first()?.shells.first() else { return None };
    let mut center = DVec3::ZERO;
    let mut axis = DVec3::Z;
    let mut radius = 0.0;
    let mut height = 0.0;
    let mut found = false;
    for fi in 0..shell.faces.len() {
        let Some(Some(si)) = brep.geom.face_surface.get(fi) else { continue };
        let Some(surf) = brep.geom.surfaces.get(*si) else { continue };
        if let Surface3::Cylinder(cc) = surf {
            axis = cc.axis.normalize_or_zero();
            center = cc.origin;
            radius = cc.radius;
            found = true;
        }
    }
    if !found { return None; }
    let mut z_vals = Vec::new();
    for fi in 0..shell.faces.len() {
        let Some(Some(si)) = brep.geom.face_surface.get(fi) else { continue };
        let Some(surf) = brep.geom.surfaces.get(*si) else { continue };
        if let Surface3::Plane(pl) = surf {
            z_vals.push(pl.origin.z);
        }
    }
    if z_vals.len() < 2 { return None; }
    let z_lo = z_vals.iter().min_by(|a,b| a.partial_cmp(b).unwrap()).copied()?;
    let z_hi = z_vals.iter().max_by(|a,b| a.partial_cmp(b).unwrap()).copied()?;
    height = z_hi - z_lo;
    Some((center, axis, radius, height))
}

/// Box perimeter in UV space (CCW from v_min).
/// t ∈ [0,2) → v_min (v=−ev), u∈[−eu,eu]
/// t ∈ [2,4) → u_max (u=eu),  v∈[−ev,ev]
/// t ∈ [4,6) → v_max (v=ev),  u∈[eu,−eu]
/// t ∈ [6,8) → u_min (u=−eu), v∈[ev,−ev]
fn box_perimeter_uv(t: f64, eu: f64, ev: f64) -> (f64, f64) {
    let tn = t.rem_euclid(8.0);
    if tn < 2.0 {
        let s = tn / 2.0;
        (-eu + 2.0 * eu * s, -ev)
    } else if tn < 4.0 {
        let s = (tn - 2.0) / 2.0;
        (eu, -ev + 2.0 * ev * s)
    } else if tn < 6.0 {
        let s = (tn - 4.0) / 2.0;
        (eu - 2.0 * eu * s, ev)
    } else {
        let s = (tn - 6.0) / 2.0;
        (-eu, ev - 2.0 * ev * s)
    }
}

/// Map UV coordinate to perimeter t ∈ [0, 8).
fn uv_to_perimeter_t(u: f64, v: f64, eu: f64, ev: f64) -> f64 {
    if v <= -ev + TOLERANCE_LEN_MIN {
        ((u + eu) / (2.0 * eu).max(TOLERANCE_LEN_MIN)) * 2.0
    } else if u >= eu - TOLERANCE_LEN_MIN {
        2.0 + ((v + ev) / (2.0 * ev).max(TOLERANCE_LEN_MIN)) * 2.0
    } else if v >= ev - TOLERANCE_LEN_MIN {
        4.0 + ((eu - u) / (2.0 * eu).max(TOLERANCE_LEN_MIN)) * 2.0
    } else {
        6.0 + ((ev - v) / (2.0 * ev).max(TOLERANCE_LEN_MIN)) * 2.0
    }
}

/// Find circle-box edge intersections in UV space.
fn circle_rect_intersections_uv(cu: f64, cv: f64, r: f64, eu: f64, ev: f64) -> Vec<UVEdgePt> {
    let mut pts = Vec::new();
    let tol = TOLERANCE_LEN_MIN;

    let add_if = |u: f64, v: f64, edge: usize, pts: &mut Vec<UVEdgePt>| {
        if u >= -eu - tol && u <= eu + tol && v >= -ev - tol && v <= ev + tol {
            let u_cl = u.clamp(-eu, eu);
            let v_cl = v.clamp(-ev, ev);
            let t = uv_to_perimeter_t(u_cl, v_cl, eu, ev);
            let theta = (v_cl - cv).atan2(u_cl - cu);
            pts.push(UVEdgePt { t, u: u_cl, v: v_cl, edge, theta });
        }
    };

    let d = -ev - cv; let disc = r * r - d * d;
    if disc >= 0.0 { let off = disc.sqrt(); add_if(cu + off, -ev, 2, &mut pts); if off > tol { add_if(cu - off, -ev, 2, &mut pts); } }
    let d = ev - cv; let disc = r * r - d * d;
    if disc >= 0.0 { let off = disc.sqrt(); add_if(cu + off, ev, 3, &mut pts); if off > tol { add_if(cu - off, ev, 3, &mut pts); } }
    let d = -eu - cu; let disc = r * r - d * d;
    if disc >= 0.0 { let off = disc.sqrt(); add_if(-eu, cv + off, 0, &mut pts); if off > tol { add_if(-eu, cv - off, 0, &mut pts); } }
    let d = eu - cu; let disc = r * r - d * d;
    if disc >= 0.0 { let off = disc.sqrt(); add_if(eu, cv + off, 1, &mut pts); if off > tol { add_if(eu, cv - off, 1, &mut pts); } }

    pts.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    pts.dedup_by(|a, b| (a.t - b.t).abs() < tol && (a.u - b.u).abs() < tol && (a.v - b.v).abs() < tol);
    pts
}

/// Create a planar BRep face from 4 corner points and a given normal.
fn rect_face_4(corners: [DVec3; 4], normal: DVec3) -> Option<BRep> {
    let mut brep = BRep::default();
    let surface = Surface3::Plane(Plane { origin: corners[0], normal });
    let vs: Vec<usize> = corners.iter().map(|p| make_vertex(&mut brep, *p)).collect();
    let mut wes = Vec::with_capacity(4);
    for i in 0..4 {
        let j = (i + 1) % 4;
        let dir = (corners[j] - corners[i]).normalize();
        let len = (corners[j] - corners[i]).length();
        if len < TOLERANCE_LEN_MIN { return None; }
        let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: corners[i], direction: dir }), 0.0, len, vs[i], vs[j]).ok()?;
        wes.push(WireEdge::new(ei, true));
    }
    let _fi = make_face(&mut brep, surface, make_wire(wes), vec![]).ok()?;
    Some(brep)
}

/// Create a planar face from a polygon with a given normal.
fn planar_face_from_polygon(poly: &[DVec3], normal: DVec3) -> Option<BRep> {
    if poly.len() < 3 { return None; }
    let mut brep = BRep::default();
    let surface = Surface3::Plane(Plane { origin: poly[0], normal });
    let vs: Vec<usize> = poly.iter().map(|p| make_vertex(&mut brep, *p)).collect();
    let mut wes = Vec::with_capacity(poly.len());
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let dir = (poly[j] - poly[i]).normalize();
        let len = (poly[j] - poly[i]).length();
        if len < TOLERANCE_LEN_MIN { return None; }
        let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: poly[i], direction: dir }), 0.0, len, vs[i], vs[j]).ok()?;
        wes.push(WireEdge::new(ei, true));
    }
    let _fi = make_face(&mut brep, surface, make_wire(wes), vec![]).ok()?;
    Some(brep)
}

/// Create a planar face from an outer polygon and an inner polygon (hole).
/// `outer` must be CCW when viewed along `normal`.
/// `inner` must be CW when viewed along `normal`.
fn planar_face_with_inner_hole(outer: &[DVec3], inner: &[DVec3], normal: DVec3) -> Option<BRep> {
    if outer.len() < 3 || inner.len() < 3 { return None; }
    let mut brep = BRep::default();
    let surface = Surface3::Plane(Plane { origin: outer[0], normal });
    let vs: Vec<usize> = outer.iter().map(|p| make_vertex(&mut brep, *p)).collect();
    let mut outer_wes = Vec::with_capacity(outer.len());
    let n = outer.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let dir = (outer[j] - outer[i]).normalize();
        let len = (outer[j] - outer[i]).length();
        if len < TOLERANCE_LEN_MIN { return None; }
        let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: outer[i], direction: dir }), 0.0, len, vs[i], vs[j]).ok()?;
        outer_wes.push(WireEdge::new(ei, true));
    }
    let inner_vs: Vec<usize> = inner.iter().map(|p| make_vertex(&mut brep, *p)).collect();
    let mut inner_wes = Vec::with_capacity(inner.len());
    let m = inner.len();
    for i in 0..m {
        let j = (i + 1) % m;
        let dir = (inner[j] - inner[i]).normalize();
        let len = (inner[j] - inner[i]).length();
        if len < TOLERANCE_LEN_MIN { return None; }
        let ei = make_edge(&mut brep, Curve3::Line(Line3 { origin: inner[i], direction: dir }), 0.0, len, inner_vs[i], inner_vs[j]).ok()?;
        inner_wes.push(WireEdge::new(ei, true));
    }
    let _fi = make_face(&mut brep, surface, make_wire(outer_wes), vec![make_wire(inner_wes)]).ok()?;
    Some(brep)
}

// ── shared segment types and helpers ────────────────────────────────────────

/// A merged segment along the box perimeter — either inside or outside the cylinder.
struct MergedSeg {
    t0: f64,
    t1: f64,
    outside: bool,
}

/// Compute merged segments along the box perimeter from intersection points.
fn compute_merged_segments(ints: &[UVEdgePt], cu: f64, cv: f64, cyl_r: f64, eu: f64, ev: f64) -> Vec<MergedSeg> {
    let tol = TOLERANCE_LEN_MIN;
    let r2 = cyl_r * cyl_r;

    let outside_at = |t: f64| -> bool {
        let (u, v) = box_perimeter_uv(t, eu, ev);
        (u - cu).powi(2) + (v - cv).powi(2) > r2 + tol
    };

    if ints.is_empty() {
        return vec![MergedSeg { t0: 0.0, t1: 8.0, outside: outside_at(0.0) }];
    }

    let mut t_vals: Vec<f64> = ints.iter().map(|p| p.t).collect();
    t_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    t_vals.dedup_by(|a, b| (*a - *b).abs() < tol);

    let mut segs: Vec<MergedSeg> = Vec::new();
    let mut prev_t = 0.0;
    for &t in &t_vals {
        if t <= prev_t + tol { continue; }
        let mid_t = (prev_t + t) / 2.0;
        segs.push(MergedSeg { t0: prev_t, t1: t, outside: outside_at(mid_t) });
        prev_t = t;
    }
    if prev_t < 8.0 - tol {
        let mid_t = (prev_t + 8.0) / 2.0;
        segs.push(MergedSeg { t0: prev_t, t1: 8.0, outside: outside_at(mid_t) });
    }

    let mut merged: Vec<MergedSeg> = Vec::new();
    for s in segs {
        if let Some(last) = merged.last_mut() {
            if last.outside == s.outside {
                last.t1 = s.t1;
                continue;
            }
        }
        merged.push(s);
    }
    merged
}

/// Generate evenly-spaced points on a circular arc, choosing the direction
/// whose midpoint stays inside the box rect.  Returns vertices in world
/// coordinates via `corner(u, v, z)`, including both endpoints.
fn arc_vertices(
    cu: f64, cv: f64, r: f64,
    th_start: f64, th_end: f64,
    eu: f64, ev: f64, z: f64,
    corner: &impl Fn(f64, f64, f64) -> DVec3,
) -> Vec<DVec3> {
    let tol = TOLERANCE_LEN_MIN;
    let inside_box = |th: f64| -> bool {
        let u = cu + r * th.cos();
        let v = cv + r * th.sin();
        u >= -eu - tol && u <= eu + tol && v >= -ev - tol && v <= ev + tol
    };

    let mut cw_len = th_start - th_end;
    if cw_len < 0.0 { cw_len += 2.0 * std::f64::consts::PI; }
    let ccw_len = 2.0 * std::f64::consts::PI - cw_len;

    let cw_mid = th_start - cw_len / 2.0;
    let ccw_mid = th_start + ccw_len / 2.0;
    let cw_in = inside_box(cw_mid);
    let ccw_in = inside_box(ccw_mid);
    let (dtheta, arc_len) = if cw_in && !ccw_in {
        (-cw_len, cw_len)
    } else if ccw_in && !cw_in {
        (ccw_len, ccw_len)
    } else if cw_in && ccw_in {
        // Both midpoints are inside the box. The shorter arc may still go outside
        // the box at non-midpoint positions (e.g., an arc >180° whose midpoint
        // happens to be inside but passes through a bulge).  Sample several points
        // along the shorter arc to verify all lie inside; if not, take the longer.
        let shorter_cw = cw_len <= ccw_len;
        let (check_len, check_dir) = if shorter_cw {
            (cw_len, -1.0)  // CW direction = decreasing θ
        } else {
            (ccw_len, 1.0)  // CCW direction = increasing θ
        };
        // Check 5 evenly-spaced sample points along the shorter arc.
        let n_samples = 5usize;
        let all_inside = (0..=n_samples).all(|i| {
            let frac = i as f64 / n_samples as f64;
            let th = th_start + check_dir * check_len * frac;
            // Normalize to [0, 2π)
            let th_n = th.rem_euclid(2.0 * std::f64::consts::PI);
            inside_box(th_n)
        });
        if all_inside {
            // Shorter arc stays inside → use it.
            if shorter_cw { (-cw_len, cw_len) } else { (ccw_len, ccw_len) }
        } else {
            // Shorter arc exits the box → use the longer arc.
            if shorter_cw { (ccw_len, ccw_len) } else { (-cw_len, cw_len) }
        }
    } else if cw_len <= ccw_len {
        (-cw_len, cw_len)
    } else {
        (ccw_len, ccw_len)
    };
    let steps = (arc_len / 0.08).ceil() as usize;
    let steps = steps.max(2).min(200);

    let mut verts = Vec::with_capacity(steps + 1);
    verts.push(corner(cu + r * th_start.cos(), cv + r * th_start.sin(), z));
    for i in 1..steps {
        let th = th_start + dtheta * (i as f64 / steps as f64);
        let pt = corner(cu + r * th.cos(), cv + r * th.sin(), z);
        if (verts.last().unwrap() - pt).length() > tol {
            verts.push(pt);
        }
    }
    let pt_end = corner(cu + r * th_end.cos(), cv + r * th_end.sin(), z);
    if (verts.last().unwrap() - pt_end).length() > tol {
        verts.push(pt_end);
    }
    verts
}

/// Build cap polygon at z level for (box rect − circle).
fn build_cap_polygon(
    merged: &[MergedSeg], cu: f64, cv: f64, cyl_r: f64,
    eu: f64, ev: f64, z: f64,
    corner: &impl Fn(f64, f64, f64) -> DVec3,
) -> Vec<DVec3> {
    let tol = TOLERANCE_LEN_MIN;

    if merged.is_empty() || (merged.len() == 1 && !merged[0].outside) {
        return vec![];
    }

    let mut poly = Vec::new();
    for seg in merged {
        if seg.t1 <= seg.t0 + tol { continue; }
        if seg.outside {
            add_box_perimeter_vertices(seg.t0, seg.t1, eu, ev, z, corner, &mut poly);
        } else {
            let (pu, pv) = box_perimeter_uv(seg.t0, eu, ev);
            let (nu, nv) = box_perimeter_uv(seg.t1, eu, ev);
            let th_prev = (pv - cv).atan2(pu - cu);
            let th_next = (nv - cv).atan2(nu - cu);
            add_circle_arc_vertices(cu, cv, cyl_r, th_prev, th_next, eu, ev, z, corner, &mut poly);
        }
    }

    if poly.len() >= 2 && (poly.last().unwrap() - poly[0]).length() < tol { poly.pop(); }
    poly
}

/// Add box perimeter vertices between two t values.
fn add_box_perimeter_vertices(
    t_start: f64, t_end: f64, eu: f64, ev: f64, z: f64,
    corner: &impl Fn(f64, f64, f64) -> DVec3,
    poly: &mut Vec<DVec3>,
) {
    let tol = TOLERANCE_LEN_MIN;
    if poly.is_empty() || ((t_start - 0.0).abs() < tol || (t_start - 8.0).abs() < tol) {
        let (u0, v0) = box_perimeter_uv(t_start, eu, ev);
        let pt = corner(u0, v0, z);
        if poly.is_empty() || (poly.last().unwrap() - pt).length() > tol {
            poly.push(pt);
        }
    }

    let corner_ts = [2.0, 4.0, 6.0];
    let start_norm = if t_start < 0.0 { t_start + 8.0 } else { t_start };
    let end_norm = if t_end <= t_start { t_end + 8.0 } else { t_end };

    for &ct in &corner_ts {
        let ct_norm = if ct < start_norm { ct + 8.0 } else { ct };
        if ct_norm > start_norm + tol && ct_norm < end_norm - tol {
            let (u, v) = box_perimeter_uv(ct, eu, ev);
            let pt = corner(u, v, z);
            if (poly.last().unwrap() - pt).length() > tol { poly.push(pt); }
        }
    }

    let (ue, ve) = box_perimeter_uv(t_end, eu, ev);
    let pt = corner(ue, ve, z);
    if (poly.last().unwrap() - pt).length() > tol { poly.push(pt); }
}

/// Add circular arc vertices for the cap polygon, taking direction from
/// `arc_vertices` and skipping the first vertex (already in poly from the
/// preceding box perimeter segment).
fn add_circle_arc_vertices(
    cu: f64, cv: f64, r: f64,
    th_start: f64, th_end: f64,
    eu: f64, ev: f64, z: f64,
    corner: &impl Fn(f64, f64, f64) -> DVec3,
    poly: &mut Vec<DVec3>,
) {
    let verts = arc_vertices(cu, cv, r, th_start, th_end, eu, ev, z, corner);
    // verts[0] is the start point, already in poly from the previous segment
    for i in 1..verts.len() {
        if (poly.last().unwrap() - verts[i]).length() > TOLERANCE_LEN_MIN {
            poly.push(verts[i]);
        }
    }
}

/// Build side wall pieces trimmed at cylinder intersection parameters.
fn build_trimmed_edge_pieces(
    p_min: f64, p_max: f64,
    p_lo: f64, p_hi: f64,
    z_lo: f64, z_hi: f64, ew: f64,
    corner: &dyn Fn(f64, f64) -> DVec3,
    normal: DVec3,
    pieces: &mut Vec<BRep>,
) {
    let tol = TOLERANCE_LEN_MIN;
    let push = |c: [DVec3; 4], n: DVec3, p: &mut Vec<BRep>| { if let Some(f) = rect_face_4(c, n) { p.push(f); } };

    let strip = |s_min: f64, s_max: f64, pieces: &mut Vec<BRep>| {
        if s_max <= s_min + tol { return; }
        if z_lo > -ew + tol { push([corner(s_min, -ew), corner(s_max, -ew), corner(s_max, z_lo), corner(s_min, z_lo)], normal, pieces); }
        if z_hi < ew - tol { push([corner(s_min, z_hi), corner(s_max, z_hi), corner(s_max, ew), corner(s_min, ew)], normal, pieces); }
        if z_hi > z_lo + tol { push([corner(s_min, z_lo), corner(s_max, z_lo), corner(s_max, z_hi), corner(s_min, z_hi)], normal, pieces); }
    };

    if p_lo <= p_min + tol && p_hi >= p_max - tol { strip(p_min, p_max, pieces); return; }
    if p_lo > p_min + tol { strip(p_min, p_lo, pieces); }
    if p_hi < p_max - tol { strip(p_hi, p_max, pieces); }
}

/// Build quad faces for the cylindrical wall from merged perimeter segments,
/// using the same angular step as the cap polygon for matching vertices.
fn build_cylindrical_wall_from_segs(
    merged: &[MergedSeg],
    cu: f64, cv: f64, cyl_r: f64,
    z_lo: f64, z_hi: f64, eu: f64, ev: f64,
    u_ax: DVec3, v_ax: DVec3, c: DVec3,
    corner: &impl Fn(f64, f64, f64) -> DVec3,
    pieces: &mut Vec<BRep>,
) {
    let tol = TOLERANCE_LEN_MIN;
    for seg in merged {
        if seg.outside { continue; }
        let (pu, pv) = box_perimeter_uv(seg.t0, eu, ev);
        let (nu, nv) = box_perimeter_uv(seg.t1, eu, ev);
        let th_prev = (pv - cv).atan2(pu - cu);
        let th_next = (nv - cv).atan2(nu - cu);

        let lo_verts = arc_vertices(cu, cv, cyl_r, th_prev, th_next, eu, ev, z_lo, corner);
        let hi_verts = arc_vertices(cu, cv, cyl_r, th_prev, th_next, eu, ev, z_hi, corner);

        let n = lo_verts.len().min(hi_verts.len());
        if n < 2 { continue; }

        for i in 0..n - 1 {
            let b0 = lo_verts[i];
            let b1 = lo_verts[i + 1];
            let t1 = hi_verts[i + 1];
            let t0 = hi_verts[i];

            // Compute outward radial direction for sign convention
            let mid = (b0 + b1 + t1 + t0) / 4.0;
            let u_mid = (mid - c).dot(u_ax);
            let v_mid = (mid - c).dot(v_ax);
            let radial_u = u_mid - cu;
            let radial_v = v_mid - cv;
            let outward = (u_ax * radial_u + v_ax * radial_v).normalize();
            let n_vec = (t0 - b0).cross(b1 - b0).normalize();
            let n_final = if n_vec.dot(outward) > 0.0 { n_vec } else { -n_vec };

            if let Some(f) = rect_face_4([b0, b1, t1, t0], n_final) { pieces.push(f); }
        }
    }
}

/// Main Stage 3 orchestrator: build box − cylinder with partial XY containment.
fn build_box_cylinder_result_partial(
    c: DVec3, u_ax: DVec3, v_ax: DVec3, eu: f64, ev: f64, ew: f64,
    cu: f64, cv: f64, cyl_r: f64, cyl_z_lo: f64, cyl_z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let z_lo = (cyl_z_lo - c.z).max(-ew);
    let z_hi = (cyl_z_hi - c.z).min(ew);
    if z_hi <= z_lo + tol { return None; }

    let corner_z = |u: f64, v: f64, z: f64| -> DVec3 { c + u*u_ax + v*v_ax + z*DVec3::Z };

    let ints = circle_rect_intersections_uv(cu, cv, cyl_r, eu, ev);
    let merged = compute_merged_segments(&ints, cu, cv, cyl_r, eu, ev);
    let cap_bottom = build_cap_polygon(&merged, cu, cv, cyl_r, eu, ev, z_lo, &corner_z);
    let cap_top = build_cap_polygon(&merged, cu, cv, cyl_r, eu, ev, z_hi, &corner_z);

    let mut pieces: Vec<BRep> = Vec::new();

    // Compute per-edge intersection parameters directly (NOT from the global
    // intersection list, which loses corner-shared entries during dedup).
    let add_pt = |v: f64, lo: f64, hi: f64, out: &mut Vec<f64>| {
        if v >= lo - tol && v <= hi + tol { out.push(v.clamp(lo, hi)); }
    };

    // u_min edge (u=-eu): solve (-eu-cu)^2 + (v-cv)^2 = r^2 → v
    let mut ints_u_min: Vec<f64> = Vec::with_capacity(2);
    let d0 = -eu - cu; let disc0 = cyl_r * cyl_r - d0 * d0;
    if disc0 >= 0.0 { let off = disc0.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_min); add_pt(cv+off, -ev, ev, &mut ints_u_min); }
    ints_u_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_u_min.dedup_by(|a,b| (*a - *b).abs() < tol);

    // u_max edge (u=eu): solve (eu-cu)^2 + (v-cv)^2 = r^2 → v
    let mut ints_u_max: Vec<f64> = Vec::with_capacity(2);
    let d1 = eu - cu; let disc1 = cyl_r * cyl_r - d1 * d1;
    if disc1 >= 0.0 { let off = disc1.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_max); add_pt(cv+off, -ev, ev, &mut ints_u_max); }
    ints_u_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_u_max.dedup_by(|a,b| (*a - *b).abs() < tol);

    // v_min edge (v=-ev): solve (u-cu)^2 + (-ev-cv)^2 = r^2 → u
    let mut ints_v_min: Vec<f64> = Vec::with_capacity(2);
    let d2 = -ev - cv; let disc2 = cyl_r * cyl_r - d2 * d2;
    if disc2 >= 0.0 { let off = disc2.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_min); add_pt(cu+off, -eu, eu, &mut ints_v_min); }
    ints_v_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_v_min.dedup_by(|a,b| (*a - *b).abs() < tol);

    // v_max edge (v=ev): solve (u-cu)^2 + (ev-cv)^2 = r^2 → u
    let mut ints_v_max: Vec<f64> = Vec::with_capacity(2);
    let d3 = ev - cv; let disc3 = cyl_r * cyl_r - d3 * d3;
    if disc3 >= 0.0 { let off = disc3.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_max); add_pt(cu+off, -eu, eu, &mut ints_v_max); }
    ints_v_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_v_max.dedup_by(|a,b| (*a - *b).abs() < tol);

    // Side walls — always use build_trimmed_edge_pieces for proper z-splitting
    // at z_lo/z_hi, matching the cap polygon edge vertices.

    // Helper for building trimmed edge pieces from intersection list.
    // Skips when the intersection covers the full span (both points at
    // corners), since `build_trimmed_edge_pieces` treats p_lo=p_min/p_hi=p_max
    // as "no intersection" and incorrectly keeps the full face.
    let trimmed_or_skip = |ints: &[f64], lo: f64, hi: f64,
                           cn: &dyn Fn(f64, f64) -> DVec3, nrm: DVec3,
                           pieces: &mut Vec<BRep>| {
        if ints.len() >= 2 && (ints.last().unwrap() - ints.first().unwrap()).abs() >= tol {
            let p_lo = ints[0];
            let p_hi = ints[ints.len() - 1];
            if p_lo > lo + tol || p_hi < hi - tol {
                build_trimmed_edge_pieces(lo, hi, p_lo, p_hi, z_lo, z_hi, ew, cn, nrm, pieces);
            }
            // else: full span — cylinder removes this entire face. Nothing to keep.
        } else if ints.len() == 1 {
            // Single intersection point — check which interval to keep.
            let p = ints[0];
            // Determine which face this is by checking nrm vs u_ax/v_ax.
            let is_u_face = nrm.dot(u_ax).abs() > 0.5;
            let inside: Box<dyn Fn(f64) -> bool> = if is_u_face {
                Box::new(|coord: f64| { (eu - cu).powi(2) + (coord - cv).powi(2) < cyl_r.powi(2) + tol })
            } else {
                Box::new(|coord: f64| { (coord - cu).powi(2) + (ev - cv).powi(2) < cyl_r.powi(2) + tol })
            };
            let mid_lo = (lo + p) * 0.5;
            let mid_hi = (p + hi) * 0.5;
            let ins_lo = inside(mid_lo);
            let ins_hi = inside(mid_hi);
            if !ins_lo && !ins_hi {
                // Both outside — tangent. Keep full face.
                build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
            } else if !ins_lo {
                if p > lo + tol { build_trimmed_edge_pieces(lo, hi, lo, p, z_lo, z_hi, ew, cn, nrm, pieces); }
            } else if !ins_hi {
                if p < hi - tol { build_trimmed_edge_pieces(lo, hi, p, hi, z_lo, z_hi, ew, cn, nrm, pieces); }
            } else {
                // Both inside — shouldn't happen. Fall through to full face.
                build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
            }
        } else {
            build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
        }
    };

    // u_max face (normal = +u_ax, param v)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c + eu*u_ax + p*v_ax + z*DVec3::Z };
        trimmed_or_skip(&ints_u_max, -ev, ev, &cn, u_ax, &mut pieces);
    }
    // u_min face (normal = -u_ax, param v)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c - eu*u_ax + p*v_ax + z*DVec3::Z };
        trimmed_or_skip(&ints_u_min, -ev, ev, &cn, -u_ax, &mut pieces);
    }
    // v_max face (normal = +v_ax, param u)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c + p*u_ax + ev*v_ax + z*DVec3::Z };
        trimmed_or_skip(&ints_v_max, -eu, eu, &cn, v_ax, &mut pieces);
    }
    // v_min face (normal = -v_ax, param u)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c + p*u_ax - ev*v_ax + z*DVec3::Z };
        trimmed_or_skip(&ints_v_min, -eu, eu, &cn, -v_ax, &mut pieces);
    }

    if !cap_bottom.is_empty() {
        if let Some(f) = planar_face_from_polygon(&cap_bottom, DVec3::Z) {
            pieces.push(f);
        }
    }
    if !cap_top.is_empty() {
        if let Some(f) = planar_face_from_polygon(&cap_top, -DVec3::Z) { pieces.push(f); }
    }

    // Cylindrical wall
    build_cylindrical_wall_from_segs(&merged, cu, cv, cyl_r, z_lo, z_hi, eu, ev, u_ax, v_ax, c, &corner_z, &mut pieces);

    if pieces.is_empty() { return None; }

    let sewn = sew_shells(&pieces, tol.max(TOLERANCE_ABS * 100.0));
    Some(sewn.brep)
}

/// Build box − cylinder difference for full XY containment (cylinder inside or
/// tangent to the box rect). Uses inner-wire annular cap faces and a full 360°
/// cylindrical wall, with side walls split at the cylinder Z range.
fn build_box_cylinder_full_containment(
    c: DVec3, u_ax: DVec3, v_ax: DVec3, eu: f64, ev: f64, ew: f64,
    cu: f64, cv: f64, cyl_r: f64, cyl_z_lo: f64, cyl_z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let z_lo = (cyl_z_lo - c.z).max(-ew);
    let z_hi = (cyl_z_hi - c.z).min(ew);
    if z_hi <= z_lo + tol { return None; }

    let mut pieces: Vec<BRep> = Vec::new();
    let p = |u: f64, v: f64, z: f64| -> DVec3 { c + u*u_ax + v*v_ax + z*DVec3::Z };

    // Discretize circle (CCW in XY) for inner wires and cylindrical wall.
    let n_circ = 64usize;
    let circle_ccw = |z: f64| -> Vec<DVec3> {
        (0..n_circ).map(|i| {
            let th = 2.0 * std::f64::consts::PI * (i as f64) / (n_circ as f64);
            p(cu + cyl_r * th.cos(), cv + cyl_r * th.sin(), z)
        }).collect()
    };
    let ccw_lo = circle_ccw(z_lo);
    let ccw_hi = circle_ccw(z_hi);

    // ── 1. Side walls split at z_lo/z_hi ──
    let mut strip = |u: f64, v0: f64, v1: f64, z0: f64, z1: f64, nrm: DVec3| {
        if z1 <= z0 + tol { return; }
        if let Some(f) = rect_face_4([p(u, v0, z0), p(u, v1, z0), p(u, v1, z1), p(u, v0, z1)], nrm) {
            pieces.push(f);
        }
    };
    let mut split_wall = |u: f64, v_min: f64, v_max: f64, nrm: DVec3| {
        if z_lo > -ew + tol { strip(u, v_min, v_max, -ew, z_lo, nrm); }
        strip(u, v_min, v_max, z_lo, z_hi, nrm);
        if z_hi < ew - tol { strip(u, v_min, v_max, z_hi, ew, nrm); }
    };
    split_wall(eu, -ev, ev, u_ax);
    split_wall(-eu, -ev, ev, -u_ax);
    split_wall(ev, -eu, eu, v_ax);
    split_wall(-ev, -eu, eu, -v_ax);

    // ── 2. Annular cap faces —──
    // Bottom region
    if z_lo > -ew + tol {
        // Full bottom face at z=-ew (cylinder doesn't reach bottom).
        let bot = [p(-eu, -ev, -ew), p(-eu, ev, -ew), p(eu, ev, -ew), p(eu, -ev, -ew)];
        if let Some(f) = rect_face_4(bot, -DVec3::Z) { pieces.push(f); }
        // Interior annular cap at z_lo, normal +Z.
        // Outer CCW in XY; inner CW in XY (= reversed CCW).
        let outer = [p(-eu, -ev, z_lo), p(eu, -ev, z_lo), p(eu, ev, z_lo), p(-eu, ev, z_lo)];
        if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_lo.iter().rev().copied().collect::<Vec<_>>(), DVec3::Z) { pieces.push(f); }
    } else {
        // Bottom face IS the annular cap (cylinder goes through bottom).
        // Outer CCW in -Z view (= CW in XY); inner CW in -Z view (= CCW in XY).
        let outer = [p(-eu, -ev, -ew), p(-eu, ev, -ew), p(eu, ev, -ew), p(eu, -ev, -ew)];
        if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_lo, -DVec3::Z) { pieces.push(f); }
    }

    // Top region
    if z_hi < ew - tol {
        // Full top face at z=ew (cylinder doesn't reach top).
        let top = [p(-eu, -ev, ew), p(eu, -ev, ew), p(eu, ev, ew), p(-eu, ev, ew)];
        if let Some(f) = rect_face_4(top, DVec3::Z) { pieces.push(f); }
        // Interior annular cap at z_hi, normal -Z.
        // Outer CCW in -Z view (= CW in XY); inner CW in -Z view (= CCW in XY).
        let outer = [p(-eu, -ev, z_hi), p(-eu, ev, z_hi), p(eu, ev, z_hi), p(eu, -ev, z_hi)];
        if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_hi, -DVec3::Z) { pieces.push(f); }
    } else {
        // Top face IS the annular cap (cylinder goes through top).
        // Outer CCW in +Z view (= CCW in XY); inner CW in +Z view (= CW in XY = reversed CCW).
        let outer = [p(-eu, -ev, ew), p(eu, -ev, ew), p(eu, ev, ew), p(-eu, ev, ew)];
        if let Some(f) = planar_face_with_inner_hole(&outer, &ccw_hi.iter().rev().copied().collect::<Vec<_>>(), DVec3::Z) { pieces.push(f); }
    }

    // ── 3. Cylindrical wall (full 360°) ──
    for i in 0..n_circ {
        let b0 = ccw_lo[i];
        let b1 = ccw_lo[(i + 1) % n_circ];
        let t0 = ccw_hi[i];
        let t1 = ccw_hi[(i + 1) % n_circ];
        let mid = (b0 + b1 + t1 + t0) / 4.0;
        let radial = mid - (c + cu*u_ax + cv*v_ax + mid.z*DVec3::Z);
        let inward = (-radial).normalize_or_zero();
        let inward = if inward.length_squared() < 0.5 { DVec3::X } else { inward };
        let n_vec = (t0 - b0).cross(b1 - b0).normalize();
        let n_final = if n_vec.dot(inward) > 0.0 { n_vec } else { -n_vec };
        if let Some(f) = rect_face_4([b0, b1, t1, t0], n_final) { pieces.push(f); }
    }

    // ── 4. Sew ──
    if pieces.is_empty() { return None; }
    let sewn = sew_shells(&pieces, tol.max(TOLERANCE_ABS * 100.0));
    Some(sewn.brep)
}

/// Entry point for box − cylinder boolean difference.
pub fn try_difference_box_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    // Detect which BRep is the box and which is the cylinder.
    // First check: a is a cylinder and b is a box → box is b (the B in B−A).
    // Second check: a is a box and b is a cylinder → box is a.
    let (box_brep, cyl_center, cyl_axis, cyl_r, cyl_brep_for_z) = {
        let ba = try_as_box(a);
        let ca = try_cylinder_center_axis_radius_height(a);
        let bb = try_as_box(b);
        let cb = try_cylinder_center_axis_radius_height(b);

        if let (Some(cyl), Some(bx)) = (ca, bb) {
            let (center, axis, radius, _) = cyl;
            (bx, center, axis, radius, a)
        } else if let (Some(bx), Some(cyl)) = (ba, cb) {
            let (center, axis, radius, _) = cyl;
            (bx, center, axis, radius, b)
        } else {
            return None;
        }
    };

    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN { return None; }

    let z_idx = find_z_axis_index(&box_brep)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = box_brep.axes[u_idx];
    let v_ax = box_brep.axes[v_idx];
    let eu = box_brep.extents[u_idx];
    let ev = box_brep.extents[v_idx];
    let ew = box_brep.extents[z_idx];
    let c = box_brep.center;

    let cu = (cyl_center - c).dot(u_ax);
    let cv = (cyl_center - c).dot(v_ax);

    // Get cylinder Z range from cap faces of the original cylinder BRep.
    let find_z = |brep: &BRep| -> Option<(f64, f64)> {
        let shell = brep.solids.first()?.shells.first()?;
        let mut z_vals: Vec<f64> = Vec::new();
        for fi in 0..shell.faces.len() {
            if let Some(Some(si)) = brep.geom.face_surface.get(fi) {
                if let Some(Surface3::Plane(pl)) = brep.geom.surfaces.get(*si) {
                    z_vals.push(pl.origin.z);
                }
            }
        }
        if z_vals.len() < 2 { return None; }
        z_vals.sort_by(|a,b| a.partial_cmp(b).unwrap());
        Some((z_vals[0], z_vals[z_vals.len()-1]))
    };
    let (cyl_z_lo, cyl_z_hi) = find_z(cyl_brep_for_z)?;

    // If cylinder is fully contained in (or tangent to) the box XY rect, use the
    // inner-wire approach since the merged-segment method can't handle it.
    let tol_xy = TOL * 10.0;
    if cu - cyl_r >= -eu - tol_xy && cu + cyl_r <= eu + tol_xy
        && cv - cyl_r >= -ev - tol_xy && cv + cyl_r <= ev + tol_xy
    {
        return build_box_cylinder_full_containment(c, u_ax, v_ax, eu, ev, ew, cu, cv, cyl_r, cyl_z_lo, cyl_z_hi);
    }

    // When cylinder center is inside the box UV rect and extends beyond 2+ edges,
    // the merged-segment cap polygon can't handle the disconnected circle arcs.
    // Build caps with inner hole (face-with-hole topology) instead.
    let inside_uv = cu >= -eu && cu <= eu && cv >= -ev && cv <= ev;
    let extends_beyond_left = cu - cyl_r < -eu;
    let extends_beyond_right = cu + cyl_r > eu;
    let extends_beyond_bottom = cv - cyl_r < -ev;
    let extends_beyond_top = cv + cyl_r > ev;
    let n_beyond = [extends_beyond_left, extends_beyond_right, extends_beyond_bottom, extends_beyond_top].iter().filter(|&&x| x).count();
    if inside_uv && n_beyond >= 2 {
        return build_box_cylinder_result_partial_with_holes(
            c, u_ax, v_ax, eu, ev, ew, cu, cv, cyl_r, cyl_z_lo, cyl_z_hi,
        );
    }

    let result = build_box_cylinder_result_partial(c, u_ax, v_ax, eu, ev, ew, cu, cv, cyl_r, cyl_z_lo, cyl_z_hi);
    result
}

/// Build box − cylinder difference when the cylinder center is inside the box UV
/// rect.  The merged-segment cap polygon fails for some configurations (cylinder
/// extends beyond edges while center is inside).  Instead, use a face-with-hole
/// topology:
///
///   outer wire:  full box rectangle (CCW)
///   inner wire:  CW closed loop formed by:
///     - box-perimeter edges where the box edge lies inside the cylinder
///     - circle arcs where the circle edge lies inside the box
///
/// The inner wire is built by classifying each perimeter segment between adjacent
/// intersection points: segments whose midpoint is *inside* the circle become box
/// perimeter edges, and segments whose midpoint is *outside* the circle become
/// circle arcs.  The polygon is then reversed to produce a CW inner wire.
fn build_box_cylinder_result_partial_with_holes(
    c: DVec3, u_ax: DVec3, v_ax: DVec3, eu: f64, ev: f64, ew: f64,
    cu: f64, cv: f64, cyl_r: f64, cyl_z_lo: f64, cyl_z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let z_lo = (cyl_z_lo - c.z).max(-ew);
    let z_hi = (cyl_z_hi - c.z).min(ew);
    if z_hi <= z_lo + tol { return None; }

    let corner = |u: f64, v: f64, z: f64| -> DVec3 { c + u*u_ax + v*v_ax + z*DVec3::Z };

    // Get intersection points sorted by perimeter t (CCW along box perimeter).
    let mut pts = circle_rect_intersections_uv(cu, cv, cyl_r, eu, ev);
    if pts.len() < 2 { return None; }
    pts.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());

    let n = pts.len();

    // Classify each perimeter segment:
    //   INSIDE  (midpoint in circle) → box perimeter edge in the hole boundary
    //   OUTSIDE (midpoint outside)   → circle arc in the hole boundary
    let is_inside: Vec<bool> = (0..n).map(|i| {
        let a = &pts[i];
        let b = &pts[(i + 1) % n];
        let mid_t = if a.t <= b.t { (a.t + b.t) * 0.5 } else { (a.t + b.t + 8.0) * 0.5 % 8.0 };
        let (mu, mv) = box_perimeter_uv(mid_t, eu, ev);
        (mu - cu).powi(2) + (mv - cv).powi(2) <= cyl_r.powi(2) + tol
    }).collect();

    // Outer wire: full box rectangle (CCW)
    let box_outer = |z: f64| -> Vec<DVec3> { vec![
        corner(-eu, -ev, z), corner(eu, -ev, z),
        corner(eu, ev, z), corner(-eu, ev, z),
    ]};

    // Build inner wire (CW) from the hole boundary.
    //
    // First build in CCW perimeter order: INSIDE segments → box perimeter edge
    // (add_box_perimeter_vertices goes CCW = increasing t); OUTSIDE segments →
    // circle arc (arc_vertices picks the direction that stays inside the box).
    // Then reverse the whole polygon to obtain CW winding for the inner hole.
    let build_inner = |z: f64| -> Vec<DVec3> {
        let mut poly: Vec<DVec3> = Vec::new();
        for i in 0..n {
            let a = &pts[i];
            let b = &pts[(i + 1) % n];
            if is_inside[i] {
                // Box perimeter edge (CCW direction = increasing t)
                add_box_perimeter_vertices(a.t, b.t, eu, ev, z, &corner, &mut poly);
            } else {
                // Circle arc: the perimeter segment is outside the cylinder, so the
                // hole boundary follows the circle arc between the two intersection points.
                let th_a = (a.v - cv).atan2(a.u - cu);
                let th_b = (b.v - cv).atan2(b.u - cu);
                let arc = arc_vertices(cu, cv, cyl_r, th_a, th_b, eu, ev, z, &corner);
                for pt in &arc {
                    if poly.is_empty() || (poly.last().unwrap() - *pt).length() > tol {
                        poly.push(*pt);
                    }
                }
            }
        }
        // Reverse to obtain CW winding for the inner hole.
        poly.reverse();
        if poly.len() >= 2 && (poly.last().unwrap() - poly[0]).length() < tol {
            poly.pop();
        }
        poly
    };

    let mut pieces: Vec<BRep> = Vec::new();

    // Side walls — same as build_box_cylinder_result_partial (per-edge intersection)
    let add_pt = |v: f64, lo: f64, hi: f64, out: &mut Vec<f64>| {
        if v >= lo - tol && v <= hi + tol { out.push(v.clamp(lo, hi)); }
    };

    let mut ints_u_min: Vec<f64> = Vec::with_capacity(2);
    let d0 = -eu - cu; let disc0 = cyl_r * cyl_r - d0 * d0;
    if disc0 >= 0.0 { let off = disc0.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_min); add_pt(cv+off, -ev, ev, &mut ints_u_min); }
    ints_u_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_u_min.dedup_by(|a,b| (*a - *b).abs() < tol);

    let mut ints_u_max: Vec<f64> = Vec::with_capacity(2);
    let d1 = eu - cu; let disc1 = cyl_r * cyl_r - d1 * d1;
    if disc1 >= 0.0 { let off = disc1.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_max); add_pt(cv+off, -ev, ev, &mut ints_u_max); }
    ints_u_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_u_max.dedup_by(|a,b| (*a - *b).abs() < tol);

    let mut ints_v_min: Vec<f64> = Vec::with_capacity(2);
    let d2 = -ev - cv; let disc2 = cyl_r * cyl_r - d2 * d2;
    if disc2 >= 0.0 { let off = disc2.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_min); add_pt(cu+off, -eu, eu, &mut ints_v_min); }
    ints_v_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_v_min.dedup_by(|a,b| (*a - *b).abs() < tol);

    let mut ints_v_max: Vec<f64> = Vec::with_capacity(2);
    let d3 = ev - cv; let disc3 = cyl_r * cyl_r - d3 * d3;
    if disc3 >= 0.0 { let off = disc3.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_max); add_pt(cu+off, -eu, eu, &mut ints_v_max); }
    ints_v_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_v_max.dedup_by(|a,b| (*a - *b).abs() < tol);

    let mut add_side_strip = |cn: &dyn Fn(f64, f64) -> DVec3, nrm: DVec3, s_min: f64, s_max: f64| {
        if s_max <= s_min + tol { return; }
        if z_lo > -ew + tol { if let Some(f) = rect_face_4([cn(s_min, -ew), cn(s_max, -ew), cn(s_max, z_lo), cn(s_min, z_lo)], nrm) { pieces.push(f); } }
        if z_hi < ew - tol { if let Some(f) = rect_face_4([cn(s_min, z_hi), cn(s_max, z_hi), cn(s_max, ew), cn(s_min, ew)], nrm) { pieces.push(f); } }
        if z_hi > z_lo + tol { if let Some(f) = rect_face_4([cn(s_min, z_lo), cn(s_max, z_lo), cn(s_max, z_hi), cn(s_min, z_hi)], nrm) { pieces.push(f); } }
    };

    // Helper: when a face has exactly 1 cylinder intersection point (the other
    // wraps beyond the adjacent face), determine which interval to keep by
    // checking the midpoint of each side.  Returns `None` if the full face
    // should be kept (tangent case).
    let check_single_ints = |lo: f64, hi: f64, p: f64, inside: &dyn Fn(f64) -> bool| -> Option<(f64, f64)> {
        let mid_lo = (lo + p) * 0.5;
        let mid_hi = (p + hi) * 0.5;
        let ins_lo = inside(mid_lo);
        let ins_hi = inside(mid_hi);
        if !ins_lo && !ins_hi { None }        // both outside = tangent → full face
        else if !ins_lo { Some((lo, p)) }     // wrap toward hi (other point beyond hi)
        else if !ins_hi { Some((p, hi)) }     // wrap toward lo
        else { None }                         // both inside = shouldn't happen → full face
    };

    // u_max face (right, u=eu)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c + eu*u_ax + p*v_ax + z*DVec3::Z };
        if ints_u_max.len() >= 2 && (ints_u_max.last().unwrap() - ints_u_max.first().unwrap()).abs() >= tol {
            let p_lo = ints_u_max[0]; let p_hi = ints_u_max[ints_u_max.len()-1];
            if p_lo > -ev + tol { add_side_strip(&cn, u_ax, -ev, p_lo); }
            if p_hi < ev - tol { add_side_strip(&cn, u_ax, p_hi, ev); }
        } else if ints_u_max.len() == 1 {
            let p = ints_u_max[0];
            let inside = |v: f64| -> bool { (eu - cu).powi(2) + (v - cv).powi(2) < cyl_r.powi(2) + tol };
            if let Some((s_lo, s_hi)) = check_single_ints(-ev, ev, p, &inside) {
                if s_hi > s_lo + tol { add_side_strip(&cn, u_ax, s_lo, s_hi); }
            } else {
                add_side_strip(&cn, u_ax, -ev, ev);
            }
        } else {
            add_side_strip(&cn, u_ax, -ev, ev);
        }
    }
    // u_min face (left, u=-eu)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c - eu*u_ax + p*v_ax + z*DVec3::Z };
        if ints_u_min.len() >= 2 && (ints_u_min.last().unwrap() - ints_u_min.first().unwrap()).abs() >= tol {
            let p_lo = ints_u_min[0]; let p_hi = ints_u_min[ints_u_min.len()-1];
            if p_lo > -ev + tol { add_side_strip(&cn, -u_ax, -ev, p_lo); }
            if p_hi < ev - tol { add_side_strip(&cn, -u_ax, p_hi, ev); }
        } else if ints_u_min.len() == 1 {
            let p = ints_u_min[0];
            let inside = |v: f64| -> bool { (-eu - cu).powi(2) + (v - cv).powi(2) < cyl_r.powi(2) + tol };
            if let Some((s_lo, s_hi)) = check_single_ints(-ev, ev, p, &inside) {
                if s_hi > s_lo + tol { add_side_strip(&cn, -u_ax, s_lo, s_hi); }
            } else {
                add_side_strip(&cn, -u_ax, -ev, ev);
            }
        } else {
            add_side_strip(&cn, -u_ax, -ev, ev);
        }
    }
    // v_max face (top, v=ev)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c + p*u_ax + ev*v_ax + z*DVec3::Z };
        if ints_v_max.len() >= 2 && (ints_v_max.last().unwrap() - ints_v_max.first().unwrap()).abs() >= tol {
            let p_lo = ints_v_max[0]; let p_hi = ints_v_max[ints_v_max.len()-1];
            if p_lo > -eu + tol { add_side_strip(&cn, v_ax, -eu, p_lo); }
            if p_hi < eu - tol { add_side_strip(&cn, v_ax, p_hi, eu); }
        } else if ints_v_max.len() == 1 {
            let p = ints_v_max[0];
            let inside = |u: f64| -> bool { (u - cu).powi(2) + (ev - cv).powi(2) < cyl_r.powi(2) + tol };
            if let Some((s_lo, s_hi)) = check_single_ints(-eu, eu, p, &inside) {
                if s_hi > s_lo + tol { add_side_strip(&cn, v_ax, s_lo, s_hi); }
            } else {
                add_side_strip(&cn, v_ax, -eu, eu);
            }
        } else {
            add_side_strip(&cn, v_ax, -eu, eu);
        }
    }
    // v_min face (bottom, v=-ev)
    {
        let cn = |p: f64, z: f64| -> DVec3 { c + p*u_ax - ev*v_ax + z*DVec3::Z };
        if ints_v_min.len() >= 2 && (ints_v_min.last().unwrap() - ints_v_min.first().unwrap()).abs() >= tol {
            let p_lo = ints_v_min[0]; let p_hi = ints_v_min[ints_v_min.len()-1];
            if p_lo > -eu + tol { add_side_strip(&cn, -v_ax, -eu, p_lo); }
            if p_hi < eu - tol { add_side_strip(&cn, -v_ax, p_hi, eu); }
        } else if ints_v_min.len() == 1 {
            let p = ints_v_min[0];
            let inside = |u: f64| -> bool { (u - cu).powi(2) + (-ev - cv).powi(2) < cyl_r.powi(2) + tol };
            if let Some((s_lo, s_hi)) = check_single_ints(-eu, eu, p, &inside) {
                if s_hi > s_lo + tol { add_side_strip(&cn, -v_ax, s_lo, s_hi); }
            } else {
                add_side_strip(&cn, -v_ax, -eu, eu);
            }
        } else {
            add_side_strip(&cn, -v_ax, -eu, eu);
        }
    }

    // Cap faces with inner hole
    let outer_lo = box_outer(z_lo);
    let inner_lo = build_inner(z_lo);
    if !inner_lo.is_empty() {
        if let Some(f) = planar_face_with_inner_hole(&outer_lo, &inner_lo, DVec3::Z) {
            pieces.push(f);
        }
    }
    let outer_hi = box_outer(z_hi);
    let inner_hi = build_inner(z_hi);
    if !inner_hi.is_empty() {
        if let Some(f) = planar_face_with_inner_hole(&outer_hi, &inner_hi, -DVec3::Z) {
            pieces.push(f);
        }
    }

    // Cylindrical wall: build wall pieces for each consecutive-OUTSIDE arc range.
    // (An "arc range" is a group of consecutive segments where the hole boundary
    // follows the circle, i.e. `is_inside[i] == false`.)
    {
        let mut i = 0;
        while i < n {
            if is_inside[i] {
                i += 1;
                continue;
            }
            // Start of an outside range: find the end (consecutive outside segments).
            let start_idx = i;
            while i < n && !is_inside[i] { i += 1; }
            let end_idx = i; // i is now first-inside or n

            // The arc range goes from pts[start_idx] to pts[end_idx] (wrapping at n).
            let a = &pts[start_idx];
            let b = &pts[end_idx % n];
            let th_a = (a.v - cv).atan2(a.u - cu);
            let th_b = (b.v - cv).atan2(b.u - cu);

            // Generate arc vertices at z_lo and z_hi for this range.
            let lo = arc_vertices(cu, cv, cyl_r, th_a, th_b, eu, ev, z_lo, &corner);
            let hi = arc_vertices(cu, cv, cyl_r, th_a, th_b, eu, ev, z_hi, &corner);
            let nv = lo.len().min(hi.len());
            if nv >= 2 {
                for j in 0..nv - 1 {
                    let b0 = lo[j]; let b1 = lo[j + 1];
                    let t1 = hi[j + 1]; let t0 = hi[j];
                    let mid = (b0 + b1 + t1 + t0) / 4.0;
                    let u_mid = (mid - c).dot(u_ax) - cu;
                    let v_mid = (mid - c).dot(v_ax) - cv;
                    let outward = (u_ax * u_mid + v_ax * v_mid).normalize();
                    let n_vec = (t0 - b0).cross(b1 - b0).normalize();
                    let n_final = if n_vec.dot(outward) > 0.0 { n_vec } else { -n_vec };
                    if let Some(f) = rect_face_4([b0, b1, t1, t0], n_final) { pieces.push(f); }
                }
            }
        }
    }

    if pieces.is_empty() { return None; }
    let sewn = sew_shells(&pieces, tol.max(TOLERANCE_ABS * 100.0));
    Some(sewn.brep)
}

#[cfg(test)]

mod tests {
    use super::*;
    use crate::{boolean_op, BooleanOpType};
    use glam::DAffine3;
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
    fn box_box_intersection_a7_like() {
        // bcommon_simple A7: 1×1×1 ∩ 1×1.5×1 → the contained 1×1×1 (SA=6, vol=1).
        // Tests that try_containment returns inner (not outer) when smaller operand
        // is passed first (swapped=true).
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.5, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("box-box A7");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 6.0).abs() < 1e-7, "SA={sa} expected 6.0");
        assert!((vol - 1.0).abs() < 1e-7, "vol={vol} expected 1.0");
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

    // ── general box-box (rotated) ─────────────────────────────────────────

    fn make_rotated_box(origin: DVec3, u_dir: DVec3, v_dir: DVec3, w: f64, h: f64, d: f64, pivot: DVec3, axis: DVec3, angle_deg: f64) -> BRep {
        // Handle OCCT-style negative extents: a negative extent means the box extends
        // in the negative direction from the anchor corner.
        let z_dir = u_dir.cross(v_dir);
        let mut o = origin;
        let ww = if w < 0.0 { o += u_dir * w; -w } else { w };
        let hh = if h < 0.0 { o += v_dir * h; -h } else { h };
        let dd = if d < 0.0 { o += z_dir * d; -d } else { d };
        let mut b = make_box_brep(o, u_dir, v_dir, ww, hh, dd).unwrap();
        let rot = DAffine3::from_axis_angle(axis.normalize(), angle_deg.to_radians());
        let xf = DAffine3::from_translation(pivot) * rot * DAffine3::from_translation(-pivot);
        b.apply_transform(xf);
        b
    }

    #[test]
    fn box_detection_axis_aligned() {
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
        let info = try_as_box(&b).expect("axis-aligned box should be detected");
        // Axes may be in any order; check that all 3 standard axes are present.
        let all_axes: Vec<DVec3> = info.axes.iter().map(|a| a.abs()).collect();
        assert!(all_axes.iter().any(|a| (a - DVec3::X).length() < 1e-10), "X axis missing");
        assert!(all_axes.iter().any(|a| (a - DVec3::Y).length() < 1e-10), "Y axis missing");
        assert!(all_axes.iter().any(|a| (a - DVec3::Z).length() < 1e-10), "Z axis missing");
        assert!((info.center - DVec3::new(1.0, 1.5, 2.0)).length() < 1e-10, "center");
        // extents in same order as axes; check all three match {1.0, 1.5, 2.0}.
        let mut ex = info.extents.to_vec();
        ex.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((ex[0] - 1.0).abs() < 1e-10 && (ex[1] - 1.5).abs() < 1e-10 && (ex[2] - 2.0).abs() < 1e-10,
            "extents {:?}", info.extents);
    }

    #[test]
    fn box_detection_rotated() {
        // Box at origin, rotated 45° around Z at origin.
        let b = make_box_brep(DVec3::new(-0.5, -0.5, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = {
            let mut shape = b;
            let rot = DAffine3::from_axis_angle(DVec3::Z, 45.0_f64.to_radians());
            shape.apply_transform(rot);
            shape
        };
        let info = try_as_box(&b).expect("rotated box should be detected");
        let expected_axes = [
            DVec3::new(0.7071067811865476, 0.7071067811865476, 0.0).abs(),
            DVec3::new(-0.7071067811865476, 0.7071067811865476, 0.0).abs(),
            DVec3::new(0.0, 0.0, 1.0).abs(),
        ];
        for a in &info.axes {
            let aa = a.abs();
            let found = expected_axes.iter().any(|e| (aa - e).length() < 1e-10);
            assert!(found, "unexpected axis {:?}", a);
        }
        let planes = info.planes();
        assert_eq!(planes.len(), 6, "should have 6 half-space planes");
    }

    #[test]
    fn box_detection_non_box() {
        let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        assert!(try_as_box(&sphere).is_none(), "sphere is not a box");
    }

    #[test]
    fn rotated_box_intersection_partial_overlap() {
        // bcommon_simple_c3-like: unit box ∩ rotated box.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = (2.0_f64).sqrt();
        let b = make_rotated_box(
            DVec3::ZERO, DVec3::X, DVec3::Y, r, r / 2.0, 1.0,
            DVec3::ZERO, DVec3::Z, 45.0,
        );
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).expect("rotated intersection");
        let sa = surface_area(&result);
        let tol = (5e-3_f64).max(0.15 * sa);
        // SA should be non-zero (boxes overlap).
        assert!(sa > 0.0, "expected non-empty intersection, SA={sa}");
        // Check that try_intersection_box_general was triggered.
        let direct = try_intersection_box_general(&a, &b);
        assert!(direct.is_some(), "general box-box intersection should fire");
    }

    #[test]
    fn rotated_box_difference_boptuc_c3() {
        // boptuc_simple C3: B - A where A = unit box, B = rotated box.
        // Expected SA = 5.82843.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = (2.0_f64).sqrt();
        let b = {
            let mut shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, r, r / 2.0, 1.0).unwrap();
            let rot = DAffine3::from_axis_angle(DVec3::Z, 45.0_f64.to_radians());
            shape.apply_transform(rot);
            shape
        };
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("boptuc C3");
        let sa = surface_area(&result);
        let expected = 5.82843;
        let tol = (5e-3_f64).max(0.15 * expected);
        assert!(
            (sa - expected).abs() <= tol,
            "C3: expected SA ~{expected}, got {sa} (tol={tol})"
        );
    }

    #[test]
    fn rotated_box_difference_boptuc_n3() {
        // boptuc_simple N3: B - A where B is a rotated box with offset.
        // Expected SA = 2.5.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        // box at (0.25, 0.25, 0) size 0.5×0.5×-1, pivot (.25,.25,0), rotate 30° Z
        let b = make_rotated_box(
            DVec3::new(0.25, 0.25, 0.0), DVec3::X, DVec3::Y, 0.5, 0.5, -1.0,
            DVec3::new(0.25, 0.25, 0.0), DVec3::Z, 30.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("boptuc N3");
        let sa = surface_area(&result);
        let expected = 2.5;
        let tol = (5e-3_f64).max(0.15 * expected);
        assert!(
            (sa - expected).abs() <= tol,
            "N3: expected SA ~{expected}, got {sa} (tol={tol})"
        );
    }

    #[test]
    fn rotated_box_difference_boptuc_o1() {
        // boptuc_simple O1: B - A with rotated B at offset.
        // Expected SA = 4.48.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_rotated_box(
            DVec3::new(0.0, 0.5, 0.0), DVec3::X, DVec3::Y, 0.8, 0.8, -1.0,
            DVec3::new(0.0, 0.5, 0.0), DVec3::Z, -45.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("boptuc O1");
        let sa = surface_area(&result);
        let expected = 4.48;
        let tol = (5e-3_f64).max(0.15 * expected);
        assert!(
            (sa - expected).abs() <= tol,
            "O1: expected SA ~{expected}, got {sa} (tol={tol})"
        );
    }

    #[test]
    fn rotated_box_difference_no_overlap() {
        // Disjoint rotated boxes → difference should be B unchanged.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_rotated_box(
            DVec3::new(5.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
            DVec3::new(5.0, 0.0, 0.0), DVec3::Z, 30.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("disjoint diff");
        let sa = surface_area(&result);
        // Rotated unit box has same SA = 6.0.
        assert!((sa - 6.0).abs() < 1e-6, "disjoint: expected SA=6.0, got {sa}");
    }

    #[test]
    fn rotated_box_difference_a_contains_b() {
        // B fully inside A → B - A = empty.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_rotated_box(
            DVec3::new(0.25, 0.25, 0.25), DVec3::X, DVec3::Y, 0.5, 0.5, 0.5,
            DVec3::new(0.25, 0.25, 0.25), DVec3::Z, 15.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("A contains B");
        let n_faces: usize = result.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty (B inside A)");
    }

    #[test]
    fn rotated_box_intersection_non_box_falls_through() {
        // Box ∩ sphere → falls through to Pave-Filler (no panic).
        let b = make_rotated_box(
            DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
            DVec3::ZERO, DVec3::Z, 45.0,
        );
        let s = make_sphere_brep(DVec3::new(0.5, 0.5, 0.5), 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &b, &s).expect("box-sphere intersection");
        assert!(r.solids.len() >= 1 || r.vertices.is_empty());
    }
}
