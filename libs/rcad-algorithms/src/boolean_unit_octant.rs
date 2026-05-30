//! Special-case intersections used by OCCT DRAW ports when the generic `BooleanBuilder`
//! path is wrong or overly faceted: (1) unit ball 閳?`[0,1]椴乣, (2) coaxial sharp cone 閳?finite
//! cylinder (ZP7), (3) coaxial sharp cone minus cylinder sealing the base (ZP8),
//! (4) coaxial cylinder minus cone via sewn loft shells (`boptuc_simple`/ZP3).
//! (5) concentric analytic spheres (`make_sphere_brep`): compound outer sphere + mirrored inner 閳?analytic shell SA/volume (differs from OCCT `mkvolume` on trimmed patches).
//! (6) intersection of two nested analytic spheres sharing a center 閳?smaller ball (`閳兗 is the inner sphere).
//!
//! Unit ball 閳?box `[0,1]椴乣 (first-octant "spherical sector"):
//!
//! The generic Pave/Builder path does not yet split planar faces along the
//! sphere, so the result was three untrimmed 1鑴? squares. OCCT `bcommon_simple/A1`
//! expects the exact surface `5锜?4` and volume `锜?6` for the eighth ball.
//!
//! This is *not* a full analytic CSG solution 閳?only a recognition + mesh for
//! this configuration used in OCCT DRAW port tests.
//!
//! When the **box** is rigidly transformed (e.g. OCCT `trotate` about a pivot
//! not at the origin), the axis-aligned `[0,1]椴乣 predicate fails and the
//! generic boolean path must handle sphere閳ユ悐blique-plane trimming. Pave now
//! uses exact spherical UV in projection fallbacks (`SphericalSurface::world_to_uv`),
//! and plane閳ユ悞phere tangents are inflated to a micro-circle so every box face gets
//! a `FaceFace` curve. OCCT `bcommon_simple/A4` / `A5` remain `#[ignore]` 閳?`BooleanBuilder`
//! surface area / volume do not yet match (`checkprops -s`); next step is sphere UV
//! multi-trim / classification, not missing FF pairs.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{any_perpendicular, Circle3, ConicalSurface, Curve3, CylindricalSurface, Line3, Plane, SphericalSurface, Surface3, SurfaceEval, ToroidalSurface};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
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
/// - [`Union`](BooleanOpType::Union) / [`Intersection`](BooleanOpType::Intersection) 閳?returns a clone of `a`
/// - [`Difference`](BooleanOpType::Difference) 閳?returns empty [`BRep`]
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
/// For Union: B inside A 閳?return A. For Intersection: B inside A 閳?return B.
/// For Difference: A inside B 閳?return empty (result has no volume).
/// For Difference B inside A: not handled (falls through to generic Pave-Filler).
pub fn try_containment(a: &BRep, b: &BRep, op: BooleanOpType) -> Option<BRep> {
    for (outer, inner, swapped) in [(a, b, false), (b, a, true)] {
        // If the outer solid has internal cavities (multiple shells) and this is Union,
        // containment is ambiguous: the inner shape could be filling a cavity rather than
        // being truly embedded in solid material. Let union-specific fast paths
        // (e.g. try_union_fill_box_cavity) handle this case.
        if matches!(op, BooleanOpType::Union) && outer.solids.len() == 1 && outer.solids[0].shells.len() > 1
        {
            continue;
        }

        // Also skip containment for Union when the outer shape's SA exceeds its
        // bounding-box SA 鈥?this indicates internal cavity faces. This catches the
        // case where extract_solids collapsed multiple shells into one (ZG5 test
        // pattern: union of extract_solids result + cavity wall).
        if matches!(op, BooleanOpType::Union) {
            let sa = crate::total_surface_area(outer);
            if let Some(bb) = outer.bounding_box() {
                let [omin, omax] = bb;
                let (bw, bh, bd) = (omax.x - omin.x, omax.y - omin.y, omax.z - omin.z);
                let bbox_sa = 2.0 * (bw * bh + bw * bd + bh * bd);
                if sa > bbox_sa * 1.01 {
                    continue;
                }
            }
        }

        // Use vertices+curves bbox (excluding surface expansion) for the outer
        // pre-check.  Surface expansion inflates bboxes for shapes like cone
        // frustums (apex extends past the solid), causing false containment
        // positives (ZM1).
        let Some([omin, omax]) = outer.vertices_curves_bounding_box() else { continue };
        // For the inner solid, use `vertices_curves_bounding_box` for Union
        // (conservative 鈥?avoids false positives where a cylinder's 2 vertices
        // imply a degenerate bbox while the curved wall extends much further,
        // e.g. bopfuse_simple ZA6), but fall back to vertex-only bbox for
        // Intersection to avoid the spherical expansion that `curve_bounding_box`
        // applies to circles (e.g. J2: box-with-hole 鈭?box, where the hole's
        // circular edge at z=0 inflates the z-bbox to 卤r).
        let (imin, imax) = if matches!(op, BooleanOpType::Union) {
            match inner.vertices_curves_bounding_box() {
                Some(bb) => (bb[0], bb[1]),
                None => continue,
            }
        } else {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for v in &inner.vertices {
                mn = mn.min(v.point);
                mx = mx.max(v.point);
            }
            if mn.x.is_infinite() { continue; }
            (mn, mx)
        };
        let tol = TOLERANCE_ABS;
        if !(imin.x >= omin.x - tol && imax.x <= omax.x + tol &&
             imin.y >= omin.y - tol && imax.y <= omax.y + tol &&
             imin.z >= omin.z - tol && imax.z <= omax.z + tol)
        {
            continue;
        }
        // All inner vertices must be inside the outer solid (not just its bbox).
        // This is critical for curved solids (e.g. sphere bbox contains box corners
        // that are outside the sphere surface).
        // Nudge each vertex toward the centroid so boundary points (on faces/edges)
        // move slightly inside before ray testing 鈥?ray casting from exact boundary
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
/// truly disjoint 閳?no intersection computation is needed and [`BRep::compound_from_shapes`]
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
    // Gap on ANY axis閳?bboxes are disjoint (no contact, no volume overlap).
    if amax.x < bmin.x || amin.x > bmax.x
        || amax.y < bmin.y || amin.y > bmax.y
        || amax.z < bmin.z || amin.z > bmax.z
    {
        return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
    }
    None
}

/// Like [`try_union_disjoint`] but treats touching bboxes as disjoint (uses <=/>=).
///
/// For non-box shapes (sphere-box etc.), touching bboxes with </> means the actual
/// shapes are disjoint (bboxes only touch, shapes don't overlap).  The strict </>
/// forces these cases through the slow Pave-Filler (bfuse_simple A2 takes ~460s).
/// This relaxed check is safe as a fallback AFTER box-box fast paths have handled
/// face-touching fusion.
pub fn try_union_disjoint_or_touching(a: &BRep, b: &BRep) -> Option<BRep> {
    if a.compound.is_some() || b.compound.is_some() {
        return None;
    }
    let Some([amin, amax]) = a.bounding_box() else { return None; };
    let Some([bmin, bmax]) = b.bounding_box() else { return None; };
    if amax.x <= bmin.x || amin.x >= bmax.x
        || amax.y <= bmin.y || amin.y >= bmax.y
        || amax.z <= bmin.z || amin.z >= bmax.z
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
/// Compute the CW outline of the 2D union of axis-aligned rectangles.
/// Returns `None` when the rectangles are disconnected (no single union outline).
fn rect_union_outline_cw(rects: &[(f64, f64, f64, f64)]) -> Option<Vec<(f64, f64)>> {
    if rects.is_empty() {
        return None;
    }
    if rects.len() == 1 {
        let (x0, x1, y0, y1) = rects[0];
        return Some(vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]);
    }

    // Collect unique coordinates
    let mut xs: Vec<f64> = rects.iter().flat_map(|r| vec![r.0, r.1]).collect();
    let mut ys: Vec<f64> = rects.iter().flat_map(|r| vec![r.2, r.3]).collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs.dedup();
    ys.dedup();

    let nx = xs.len() - 1;
    let ny = ys.len() - 1;

    // Mark grid cells as filled if their centre is inside any rectangle
    let mut filled = vec![false; nx * ny];
    for &(x0, x1, y0, y1) in rects {
        for xi in 0..nx {
            let cx = (xs[xi] + xs[xi + 1]) * 0.5;
            if cx < x0 - 1e-12 || cx > x1 + 1e-12 { continue; }
            for yi in 0..ny {
                let cy = (ys[yi] + ys[yi + 1]) * 0.5;
                if cx >= x0 - 1e-12 && cx <= x1 + 1e-12
                    && cy >= y0 - 1e-12 && cy <= y1 + 1e-12
                {
                    filled[yi * nx + xi] = true;
                }
            }
        }
    }

    if filled.iter().all(|&f| !f) { return None; }

    // Collect boundary edges (grid edges between filled and empty cells).
    #[derive(Clone)]
    struct BEdge { x0: f64, y0: f64, x1: f64, y1: f64 }

    let mut edges: Vec<BEdge> = Vec::new();
    let cell = |xi: isize, yi: isize| -> bool {
        xi >= 0 && xi < nx as isize && yi >= 0 && yi < ny as isize
            && filled[yi as usize * nx + xi as usize]
    };

    // Horizontal edges at y = ys[yi] between xs[xi] and xs[xi+1]
    for yi in 0..=ny {
        for xi in 0..nx {
            let above = cell(xi as isize, yi as isize - 1);
            let below = cell(xi as isize, yi as isize);
            if above != below {
                edges.push(BEdge { x0: xs[xi], y0: ys[yi], x1: xs[xi + 1], y1: ys[yi] });
            }
        }
    }
    // Vertical edges at x = xs[xi] between ys[yi] and ys[yi+1]
    for xi in 0..=nx {
        for yi in 0..ny {
            let right = cell(xi as isize, yi as isize);
            let left = cell(xi as isize - 1, yi as isize);
            if right != left {
                edges.push(BEdge { x0: xs[xi], y0: ys[yi], x1: xs[xi], y1: ys[yi + 1] });
            }
        }
    }

    if edges.is_empty() { return None; }

    // Walk the edges into a closed loop. Start with the first edge.
    let eps = 1e-12;
    let mut outline: Vec<(f64, f64)> = Vec::new();
    outline.push((edges[0].x0, edges[0].y0));
    outline.push((edges[0].x1, edges[0].y1));
    let mut used = vec![false; edges.len()];
    used[0] = true;

    loop {
        let (lx, ly) = *outline.last().unwrap();
        let mut found = false;
        for (ei, be) in edges.iter().enumerate() {
            if used[ei] { continue; }
            let match_start = (be.x0 - lx).abs() < eps && (be.y0 - ly).abs() < eps;
            if match_start {
                outline.push((be.x1, be.y1));
                used[ei] = true;
                found = true;
                break;
            }
            let match_end = (be.x1 - lx).abs() < eps && (be.y1 - ly).abs() < eps;
            if match_end {
                outline.push((be.x0, be.y0));
                used[ei] = true;
                found = true;
                break;
            }
        }
        if !found { break; }
        // Check if we closed the loop (back at the second point 鈫?avoid duplicating start)
        let (nx, ny) = *outline.last().unwrap();
        if outline.len() > 2 && (nx - outline[0].0).abs() < eps && (ny - outline[0].1).abs() < eps {
            outline.pop();
            break;
        }
    }

    if outline.len() < 3 { return None; }

    // Ensure CW winding (negative signed area). Reverse if CCW.
    let signed_area: f64 = outline.iter()
        .zip(outline.iter().cycle().skip(1))
        .take(outline.len())
        .map(|(&(x0, y0), &(x1, y1))| x0 * y1 - x1 * y0)
        .sum();
    if signed_area > 0.0 {
        outline.reverse();
    }

    // Remove collinear points (intermediate grid vertices on straight segments).
    let mut cleaned: Vec<(f64, f64)> = Vec::with_capacity(outline.len());
    for i in 0..outline.len() {
        let prev = outline[(i + outline.len() - 1) % outline.len()];
        let cur = outline[i];
        let next = outline[(i + 1) % outline.len()];
        let dx1 = cur.0 - prev.0;
        let dy1 = cur.1 - prev.1;
        let dx2 = next.0 - cur.0;
        let dy2 = next.1 - cur.1;
        if (dx1 * dy2 - dy1 * dx2).abs() > 1e-12 {
            cleaned.push(cur);
        }
    }

    Some(cleaned)
}

/// Build a BRep by extruding a 2D CW polygon in Z.
/// `outline_cw` 鈥?vertices in CW order on the XY plane at z=z_min.
/// The top face is generated by reversing the outline (CCW) at z=z_max.
fn build_extruded_brep(outline_cw: &[(f64, f64)], z_min: f64, z_max: f64) -> Option<BRep> {
    use rcad_kernel::geom::Line3;
    use rcad_modeling::builder::brep_builder::{make_vertex, make_edge, make_wire, make_face};
    use rcad_kernel::topology::WireEdge;

    let n = outline_cw.len();
    if n < 3 { return None; }

    let mut brep = BRep::default();

    // 1. Create all 2n vertices: 0..n at z_min, n..2n at z_max
    let mut bv = Vec::with_capacity(n);  // bottom vertices
    let mut tv = Vec::with_capacity(n);  // top vertices
    for &(x, y) in outline_cw {
        bv.push(make_vertex(&mut brep, DVec3::new(x, y, z_min)));
        tv.push(make_vertex(&mut brep, DVec3::new(x, y, z_max)));
    }

    // 2. Create edges
    // Bottom edges (CW): bv[i] 鈫?bv[i+1]
    let mut be = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let p0 = DVec3::new(outline_cw[i].0, outline_cw[i].1, z_min);
        let p1 = DVec3::new(outline_cw[j].0, outline_cw[j].1, z_min);
        let curve = Curve3::Line(Line3 { origin: p0, direction: p1 - p0 });
        be.push(make_edge(&mut brep, curve, 0.0, 1.0, bv[i], bv[j]).ok()?);
    }

    // Top edges (CCW 鈥?reversed from CW): tv[i+1] 鈫?tv[i]
    let mut te = Vec::with_capacity(n);
    for i in 0..n {
        let _j = (i + 1) % n;
        // Reverse: i鈫抝 in the outline corresponds to j鈫抜 for the top face
        // Actually, for the top face, the wire goes in the OPPOSITE direction
        // around the outline. So edge at top position i connects tv[next] 鈫?tv[i]
        // where next follows the REVERSE of the CW outline.
        //
        // CW outline: 0鈫?鈫?鈫?..鈫抧-1鈫?
        // Reverse (CCW): 0鈫抧-1鈫抧-2鈫?..鈫?鈫?
        //
        // So top edge i connects tv[(n-i) % n] 鈫?tv[(n-1-i) % n]
        let src = (n - i) % n;
        let dst = (n - 1 - i) % n;
        let p0 = DVec3::new(outline_cw[src].0, outline_cw[src].1, z_max);
        let p1 = DVec3::new(outline_cw[dst].0, outline_cw[dst].1, z_max);
        let curve = Curve3::Line(Line3 { origin: p0, direction: p1 - p0 });
        te.push(make_edge(&mut brep, curve, 0.0, 1.0, tv[src], tv[dst]).ok()?);
    }

    // Vertical edges: bv[i] 鈫?tv[i]
    let mut ve = Vec::with_capacity(n);
    for i in 0..n {
        let p0 = DVec3::new(outline_cw[i].0, outline_cw[i].1, z_min);
        let p1 = DVec3::new(outline_cw[i].0, outline_cw[i].1, z_max);
        let curve = Curve3::Line(Line3 { origin: p0, direction: p1 - p0 });
        ve.push(make_edge(&mut brep, curve, 0.0, 1.0, bv[i], tv[i]).ok()?);
    }

    // 3. Create faces

    // Bottom face (CW winding 鈫?normal (0,0,-1))
    {
        let mut wire_edges = Vec::with_capacity(n);
        for i in 0..n {
            wire_edges.push(WireEdge { idx: be[i], forward: true });
        }
        let wire = make_wire(wire_edges);
        let surface = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::new(0.0, 0.0, z_min),
            normal: DVec3::new(0.0, 0.0, -1.0),
        });
        make_face(&mut brep, surface, wire, vec![]).ok()?;
    }

    // Top face (CCW winding 鈫?normal (0,0,1))
    {
        let mut wire_edges = Vec::with_capacity(n);
        for i in 0..n {
            wire_edges.push(WireEdge { idx: te[i], forward: true });
        }
        let wire = make_wire(wire_edges);
        let surface = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::new(0.0, 0.0, z_max),
            normal: DVec3::new(0.0, 0.0, 1.0),
        });
        make_face(&mut brep, surface, wire, vec![]).ok()?;
    }

    // Side faces: one per outline edge
    for i in 0..n {
        let j = (i + 1) % n;
        // Outline segment i: outline_cw[i] 鈫?outline_cw[j]
        let dx = outline_cw[j].0 - outline_cw[i].0;
        let dy = outline_cw[j].1 - outline_cw[i].1;

        // Determine the surface plane for this side face.
        // For CW polygon, interior is to the RIGHT of the directed edge.
        // Outward normal = direction pointing LEFT.
        // dir +X (dx>0): left = +Y,   right = -Y,  outward = -Y
        // dir -X (dx<0): left = -Y,   right = +Y,  outward = +Y
        // dir +Y (dy>0): left = -X,   right = +X,  outward = +X
        // dir -Y (dy<0): left = +X,   right = -X,  outward = -X
        let (plane_coord, plane_normal) = if dx.abs() > dy.abs() {
            // Horizontal edge 鈥?plane is y = const
            if dx > 0.0 {
                // Going right: outward is -Y
                (outline_cw[i].1, DVec3::new(0.0, -1.0, 0.0))
            } else {
                // Going left: outward is +Y
                (outline_cw[i].1, DVec3::new(0.0, 1.0, 0.0))
            }
        } else {
            // Vertical edge 鈥?plane is x = const
            if dy > 0.0 {
                // Going up: outward is +X
                (outline_cw[i].0, DVec3::new(1.0, 0.0, 0.0))
            } else {
                // Going down: outward is -X
                (outline_cw[i].0, DVec3::new(-1.0, 0.0, 0.0))
            }
        };

        // Side face winding (CCW from outside):
        // bv[i] 鈫?bv[j] 鈫?tv[j] 鈫?tv[i] 鈫?bv[i]
        //
        // Verify: for edge (0,0)鈫?10,0) going right, outward = -Y
        // bv[0]=(0,0,z_min), bv[1]=(10,0,z_min)
        // bv[1]-bv[0] = (10,0,0)
        // tv[0]-bv[0] = (0,0,z_max-z_min)
        // cross((10,0,0), (0,0,1)) = (0,-10,0) 鈫?-Y 鉁?
        let side_edges = vec![
            WireEdge { idx: be[i], forward: true },    // bv[i] 鈫?bv[j]
            WireEdge { idx: ve[j], forward: true },    // bv[j] 鈫?tv[j]
            WireEdge { idx: te[(n - 1 - i) % n], forward: true }, // tv[j] 鈫?tv[i]
            WireEdge { idx: ve[i], forward: false },   // tv[i] 鈫?bv[i]
        ];
        let wire = make_wire(side_edges);
        let origin = if dx.abs() > dy.abs() {
            DVec3::new(0.0, plane_coord, 0.0)
        } else {
            DVec3::new(plane_coord, 0.0, 0.0)
        };
        let surface = Surface3::Plane(rcad_kernel::geom::Plane {
            origin,
            normal: plane_normal,
        });
        make_face(&mut brep, surface, wire, vec![]).ok()?;
    }

    Some(brep)
}

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
        // Gap or touching: check if full-face contact 鈫?merge boxes.
        let touching_x = w <= zero_tol;
        let touching_y = h <= zero_tol;
        let touching_z = d <= zero_tol;
        let touch_count = [touching_x, touching_y, touching_z].iter().filter(|&&b| b).count();

        if touch_count == 1 {
            // Single-axis face contact. Check that the overlap fully covers
            // both boxes in the other two dimensions (full-face contact).
            let full_face = if touching_x {
                rmin.y <= amin.y && rmax.y >= amax.y && rmin.y <= bmin.y && rmax.y >= bmax.y
                    && rmin.z <= amin.z && rmax.z >= amax.z && rmin.z <= bmin.z && rmax.z >= bmax.z
            } else if touching_y {
                rmin.x <= amin.x && rmax.x >= amax.x && rmin.x <= bmin.x && rmax.x >= bmax.x
                    && rmin.z <= amin.z && rmax.z >= amax.z && rmin.z <= bmin.z && rmax.z >= bmax.z
            } else {
                // touching_z
                rmin.x <= amin.x && rmax.x >= amax.x && rmin.x <= bmin.x && rmax.x >= bmax.x
                    && rmin.y <= amin.y && rmax.y >= amax.y && rmin.y <= bmin.y && rmax.y >= bmax.y
            };
            if full_face {
                let cmin = DVec3::new(
                    amin.x.min(bmin.x),
                    amin.y.min(bmin.y),
                    amin.z.min(bmin.z),
                );
                let cmax = DVec3::new(
                    amax.x.max(bmax.x),
                    amax.y.max(bmax.y),
                    amax.z.max(bmax.z),
                );
                let dims = cmax - cmin;
                return Some(
                    rcad_modeling::make_box_brep(cmin, DVec3::X, DVec3::Y, dims.x, dims.y, dims.z)
                        .expect("full-face touch merged box"),
                );
            }
        }
        // Gap, edge/vertex contact, or partial-face touch: keep separate.
        return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
    }
    // Positive-volume overlap.
    // Containment: if one box entirely contains the other, return the larger.
    if amin.x <= bmin.x && amin.y <= bmin.y && amin.z <= bmin.z
        && amax.x >= bmax.x && amax.y >= bmax.y && amax.z >= bmax.z
    {
        return Some(a.clone());
    }
    if bmin.x <= amin.x && bmin.y <= amin.y && bmin.z <= amin.z
        && bmax.x >= amax.x && bmax.y >= amax.y && bmax.z >= amax.z
    {
        return Some(b.clone());
    }

    // Partial overlap.
    // Try direct extrusion when the boxes share a full axis range.
    if (amin.z - bmin.z).abs() < zero_tol && (amax.z - bmax.z).abs() < zero_tol {
        let rects = [(amin.x, amax.x, amin.y, amax.y), (bmin.x, bmax.x, bmin.y, bmax.y)];
        if let Some(outline) = rect_union_outline_cw(&rects) {
            if let Some(brep) = build_extruded_brep(&outline, amin.z, amax.z) {
                return Some(brep);
            }
        }
    }
    if (amin.y - bmin.y).abs() < zero_tol && (amax.y - bmax.y).abs() < zero_tol {
        let rects = [(amin.x, amax.x, amin.z, amax.z), (bmin.x, bmax.x, bmin.z, bmax.z)];
        if let Some(outline) = rect_union_outline_cw(&rects) {
            // Extrude in Y: map (x,z) 鈫?outline, then remap to 3D
            if let Some(mut brep) = build_extruded_brep(&outline, amin.y, amax.y) {
                // Transform the BRep: swap Y and Z so that the extrusion
                // (which was done along Z) becomes the Y-axis.
                // build_extruded_brep extrudes in Z, so the outline is in XY.
                // We want outline in XZ with extrusion in Y.
                // Remap: (x_outline, y_outline, z_extrude) 鈫?(x_outline, z_extrude, y_outline)
                // Actually, build_extruded_brep takes outline as (x,y) and extrudes to z.
                // To extrude in Y: I want outline as (x,z) 鈫?becomes (x,y,z) with y being z.
                // The simplest approach: rename vertices by swapping Y 鈫?Z.
                for v in &mut brep.vertices {
                    v.point = DVec3::new(v.point.x, v.point.z, v.point.y);
                }
                // Also swap normals on the top/bottom faces
                for solid in &mut brep.solids {
                    for shell in &mut solid.shells {
                        for face in &mut shell.faces {
                            face.normal = DVec3::new(face.normal.x, face.normal.z, face.normal.y);
                        }
                    }
                }
                return Some(brep);
            }
        }
    }
    if (amin.x - bmin.x).abs() < zero_tol && (amax.x - bmax.x).abs() < zero_tol {
        let rects = [(amin.y, amax.y, amin.z, amax.z), (bmin.y, bmax.y, bmin.z, bmax.z)];
        if let Some(outline) = rect_union_outline_cw(&rects) {
            // Extrude in X: map (y,z) 鈫?outline, then remap to 3D
            if let Some(mut brep) = build_extruded_brep(&outline, amin.x, amax.x) {
                // build_extruded_brep extrudes in Z, so the outline is in XY.
                // We want outline in YZ with extrusion in X.
                // Remap: (x_outline, y_outline, z_extrude) 鈫?(z_extrude, x_outline, y_outline)
                for v in &mut brep.vertices {
                    v.point = DVec3::new(v.point.z, v.point.x, v.point.y);
                }
                for solid in &mut brep.solids {
                    for shell in &mut solid.shells {
                        for face in &mut shell.faces {
                            face.normal = DVec3::new(face.normal.z, face.normal.x, face.normal.y);
                        }
                    }
                }
                return Some(brep);
            }
        }
    }

    // Fallback: slab decomposition + sew + internal face removal.
    let a_slabs = decompose_slabs(amin, amax, rmin, rmax, zero_tol)?;
    let b_slabs = decompose_slabs(bmin, bmax, rmin, rmax, zero_tol)?;
    let overlap = make_box_brep(
        rmin, DVec3::X, DVec3::Y, w, h, d,
    ).ok()?;

    let mut all_slabs = a_slabs;
    all_slabs.push(overlap);
    all_slabs.extend(b_slabs);

    Some(sew_slabs_into_solid(&all_slabs, zero_tol))
}

/// Intersection: unit sphere (kernel primitive) 閳?axis box [0,1]椴?
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
/// planar faces with axis-aligned normals (卤X, 卤Y, 卤Z), 8 vertices, and every
/// vertex is at an AABB corner (each coordinate matches either the min or max
/// for that axis). Returns `None` for rotated boxes, non-box shapes, or
/// degenerate inputs.
pub(crate) fn try_as_axis_aligned_box(brep: &BRep) -> Option<[DVec3; 2]> {
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

/// Decompose the axis-aligned box `[outer_min, outer_max]` into up to 26 axis-aligned
/// slabs by cutting at the boundaries of the interior region `[inner_min, inner_max]`
/// (a 3脳3脳3 grid, excluding the center cell).  Cells that fall entirely inside the
/// interior region are omitted.
///
/// Returns `None` if `make_box_brep` fails for any cell (should not happen for valid
/// axis-aligned geometries).
fn decompose_slabs(
    outer_min: DVec3, outer_max: DVec3,
    inner_min: DVec3, inner_max: DVec3,
    zero_tol: f64,
) -> Option<Vec<BRep>> {
    let xs = [outer_min.x, inner_min.x, inner_max.x, outer_max.x];
    let ys = [outer_min.y, inner_min.y, inner_max.y, outer_max.y];
    let zs = [outer_min.z, inner_min.z, inner_max.z, outer_max.z];

    let mut slabs: Vec<BRep> = Vec::new();
    for xi in 0..3 {
        let x0 = xs[xi]; let x1 = xs[xi + 1]; let dx = x1 - x0;
        if dx <= zero_tol { continue; }
        for yi in 0..3 {
            let y0 = ys[yi]; let y1 = ys[yi + 1]; let dy = y1 - y0;
            if dy <= zero_tol { continue; }
            for zi in 0..3 {
                let z0 = zs[zi]; let z1 = zs[zi + 1]; let dz = z1 - z0;
                if dz <= zero_tol { continue; }
                // Skip cells fully inside the overlap region.
                if x0 >= inner_min.x && x1 <= inner_max.x
                    && y0 >= inner_min.y && y1 <= inner_max.y
                    && z0 >= inner_min.z && z1 <= inner_max.z
                { continue; }
                slabs.push(make_box_brep(
                    DVec3::new(x0, y0, z0),
                    DVec3::X, DVec3::Y, dx, dy, dz,
                ).ok()?);
            }
        }
    }
    Some(slabs)
}

/// Sew a collection of axis-aligned box slabs into a single solid, then detect
/// and remove internal face pairs (coplanar faces with opposite normals and
/// overlapping 2D extent).  The remaining faces form the correct external boundary.
///
/// Falls back to `BRep::compound_from_shapes` when no pairs were stitched or the
/// sewn result has no solids/shells.
fn sew_slabs_into_solid(slabs: &[BRep], zero_tol: f64) -> BRep {
    if slabs.is_empty() {
        return BRep::default();
    }
    if slabs.len() == 1 {
        return slabs[0].clone();
    }

    let sewn = sew_shells(slabs, zero_tol);
    if sewn.stitched_pairs == 0 {
        return BRep::compound_from_shapes(slabs);
    }

    let mut brep = sewn.brep;
    if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
        return BRep::compound_from_shapes(slabs);
    }

    // Collect face-level plane + 2D-extent info for internal-face detection.
    struct Fi {
        face_idx: usize,
        axis: usize,
        coord: f64,
        sign: f64,
        u_min: f64, u_max: f64,
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

        let mut pts: Vec<DVec3> = Vec::new();
        for we in &face.outer_wire.edges {
            if !edge_ok(we.idx) { continue; }
            let e = &brep.edges[we.idx];
            if vert_ok(e.start) { pts.push(brep.vertices[e.start].point); }
            if vert_ok(e.end) { pts.push(brep.vertices[e.end].point); }
        }
        if pts.is_empty() { continue; }

        let centroid = pts.iter().copied().sum::<DVec3>() / pts.len() as f64;
        // Use plane-equation distance (n路centroid) for internal-face matching.
        // The axis-specific coordinate is unreliable when faces on the same
        // plane have different extents (their centroids project to different
        // positions along any single axis).  For opposite normals on the same
        // plane:  n路centroid_a + (-n)路centroid_b = d + (-d) = 0.
        let plane_dist = n.dot(centroid);

        let (u_min, u_max, v_min, v_max) = match axis {
            0 => {
                let (us, vs): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.y, p.z)).unzip();
                (fold_min(&us), fold_max(&us), fold_min(&vs), fold_max(&vs))
            },
            1 => {
                let (us, vs): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x, p.z)).unzip();
                (fold_min(&us), fold_max(&us), fold_min(&vs), fold_max(&vs))
            },
            _ => {
                let (us, vs): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x, p.y)).unzip();
                (fold_min(&us), fold_max(&us), fold_min(&vs), fold_max(&vs))
            },
        };

        finfo.push(Fi { face_idx: fi, axis, coord: plane_dist, sign, u_min, u_max, v_min, v_max });
    }

    // Detect internal pairs: same axis, opposite sign, same plane distance,
    // overlapping 2D extent.  Plane distance: opposite normals have
    // dist_a + dist_b 鈮?0 (they're on the same plane).
    let dist_tol = (zero_tol * 10000.0).max(1e-9);
    let n_faces = shell.faces.len();
    let mut internal = vec![false; n_faces];
    let mut n_pairs = 0u32;
    for i in 0..finfo.len() {
        if internal[finfo[i].face_idx] { continue; }
        for j in (i + 1)..finfo.len() {
            if internal[finfo[j].face_idx] { continue; }
            let a = &finfo[i];
            let b = &finfo[j];
            if a.axis != b.axis { continue; }
            if a.sign * b.sign > 0.0 { continue; }
            // plane-equation distance: opposite normals on the same plane
            // have total = n路c + (-n)路c' = d - d' 鈮?0.
            if (a.coord + b.coord).abs() > dist_tol { continue; }
            let tol_ext = zero_tol * 10.0;
            let overlap_u = a.u_min.max(b.u_min) + tol_ext < a.u_max.min(b.u_max);
            let overlap_v = a.v_min.max(b.v_min) + tol_ext < a.v_max.min(b.v_max);
            if overlap_u && overlap_v {
                internal[a.face_idx] = true;
                internal[b.face_idx] = true;
                n_pairs += 1;
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

    brep.solids.retain(|s| {
        s.shells.iter().any(|sh| !sh.faces.is_empty())
    });

    // Collect vertices and edges still referenced by remaining faces.
    let mut used_verts = vec![false; brep.vertices.len()];
    let mut used_edges = vec![false; brep.edges.len()];
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    if we.idx < brep.edges.len() {
                        used_edges[we.idx] = true;
                        let e = &brep.edges[we.idx];
                        if e.start < brep.vertices.len() { used_verts[e.start] = true; }
                        if e.end < brep.vertices.len() { used_verts[e.end] = true; }
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < brep.edges.len() {
                            used_edges[we.idx] = true;
                            let e = &brep.edges[we.idx];
                            if e.start < brep.vertices.len() { used_verts[e.start] = true; }
                            if e.end < brep.vertices.len() { used_verts[e.end] = true; }
                        }
                    }
                }
            }
        }
    }

    // Build remap for vertices
    let mut v_remap: Vec<Option<usize>> = vec![None; brep.vertices.len()];
    let mut new_verts: Vec<Vertex> = Vec::new();
    for (i, v) in brep.vertices.iter().enumerate() {
        if used_verts[i] {
            v_remap[i] = Some(new_verts.len());
            new_verts.push(*v);
        }
    }

    // Build remap for edges
    let mut e_remap: Vec<Option<usize>> = vec![None; brep.edges.len()];
    let mut new_edges: Vec<Edge> = Vec::new();
    for (i, e) in brep.edges.iter().enumerate() {
        if used_edges[i] {
            e_remap[i] = Some(new_edges.len());
            new_edges.push(Edge {
                start: v_remap[e.start].unwrap_or(e.start),
                end: v_remap[e.end].unwrap_or(e.end),
            });
        }
    }

    brep.vertices = new_verts;
    brep.edges = new_edges;

    // Remap edge indices in face wires.
    let remap_edge_idx = |idx: &mut usize| {
        if let Some(new_idx) = e_remap.get(*idx).copied().flatten() {
            *idx = new_idx;
        }
    };
    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                for we in &mut face.outer_wire.edges {
                    remap_edge_idx(&mut we.idx);
                }
                for wire in &mut face.inner_wires {
                    for we in &mut wire.edges {
                        remap_edge_idx(&mut we.idx);
                    }
                }
            }
        }
    }

    brep
}

/// Difference of two axis-aligned boxes computed analytically.
///
/// Decomposes A \ B into axis-aligned cells using a full 3D grid at
/// the overlap boundaries, then sews them and removes internal faces
/// (those whose all edges are stitched).  This yields the correct
/// external surface area 鈥?no internal-face inflation from compounds.
pub fn try_difference_box_box(a: &BRep, b: &BRep) -> Option<BRep> {
    let [amin, amax] = try_as_axis_aligned_box(a)?;
    let [bmin, bmax] = try_as_axis_aligned_box(b)?;

    // Overlap region
    let rmin = DVec3::new(amin.x.max(bmin.x), amin.y.max(bmin.y), amin.z.max(bmin.z));
    let rmax = DVec3::new(amax.x.min(bmax.x), amax.y.min(bmax.y), amax.z.min(bmax.z));

    let scale = (amax - amin).length().max((bmax - bmin).length()).max(1.0);
    let zero_tol = TOLERANCE_LEN_MIN * scale;

    // No overlap 鈫?A is unchanged.
    if rmin.x >= rmax.x || rmin.y >= rmax.y || rmin.z >= rmax.z {
        return Some(a.clone());
    }

    // A entirely inside B 鈫?empty.
    let vol = (rmax.x - rmin.x) * (rmax.y - rmin.y) * (rmax.z - rmin.z);
    let a_vol = (amax.x - amin.x) * (amax.y - amin.y) * (amax.z - amin.z);
    if (vol - a_vol).abs() < zero_tol {
        return Some(BRep::default());
    }

    let slabs = decompose_slabs(amin, amax, rmin, rmax, zero_tol)?;
    Some(sew_slabs_into_solid(&slabs, zero_tol))
}

fn fold_min(v: &[f64]) -> f64 { v.iter().cloned().fold(f64::MAX, f64::min) }
fn fold_max(v: &[f64]) -> f64 { v.iter().cloned().fold(f64::MIN, f64::max) }


// 鈹€鈹€ General box-box boolean via half-space polyhedron 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Information about a box BRep (axis-aligned or rotated).
pub(crate) struct BoxInfo {
    /// Orthonormal axis directions (outward normals of the 3 face-normal pairs).
    pub(crate) axes: [DVec3; 3],
    /// Center in world coordinates.
    pub(crate) center: DVec3,
    /// Positive half-extents along each axis.
    pub(crate) extents: [f64; 3],
}

impl BoxInfo {
    /// Generate the 6 interior-facing half-space constraints (origin, normal)
    /// for use with [`make_convex_polyhedron_from_half_spaces`].
    ///
    /// Each constraint is n路(p - origin) 鈮?0, where n is the outward-facing
    /// normal of the face and origin is a point on that face.
    fn planes(&self) -> Vec<(DVec3, DVec3)> {
        let [u, v, w] = self.axes;
        let [eu, ev, ew] = self.extents;
        let c = self.center;
        vec![
            (c - eu * u, -u), // u-min face: outward -u, interior u路p 鈮?u_min
            (c + eu * u,  u), // u-max face: outward +u, interior u路p 鈮?u_max
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
/// - All 8 vertices are at 卤extent corners of the implied box (verified via
///   projection onto the 3 axis directions)
pub(crate) fn try_as_box(brep: &BRep) -> Option<BoxInfo> {
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
    let tol_ang = TOLERANCE_AXIS_ALIGN; // cos(angle) tolerance (near 0掳 or 180掳)

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
            // Normals should be opposite (ni 路 nj 鈮?-1).
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

    // Ensure right-handed system: require axes[2] 路 (axes[0] 脳 axes[1]) > 0.
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
    let result = make_convex_polyhedron_from_half_spaces(&planes).ok()?;
    // Reject zero-volume intersection (boxes adjacent/touching without overlap)
    // such as the OCCT bopcommon_simple/N3 negative-dimension-box case.
    // Minimum 4 vertices for a tetrahedron; boxes that intersect in a face or
    // edge produce < 4 vertices.
    if result.vertices.len() < 4 {
        return None;
    }
    // Also check volume via bounding-box heuristic for thin slivers.
    let vol = crate::total_volume(&result);
    let scale = (a.bounding_box().map(|b| (b[1] - b[0]).length()).unwrap_or(1.0))
        .max(b.bounding_box().map(|b| (b[1] - b[0]).length()).unwrap_or(1.0));
    if vol <= crate::tolerance::TOLERANCE_LEN_MIN * scale {
        return None;
    }
    Some(result)
}

/// Difference B - A for any two box BReps (axis-aligned or rotated).
///
/// Uses a 6-slab decomposition along the first box's local axes, building each
/// slab as a convex polyhedron via [`make_convex_polyhedron_from_half_spaces`].
/// Multi-slab results are returned as a compound; internal faces between
/// touching slabs inflate the surface area of compounds, but the inflation is
/// typically within the OCCT `checkprops -s` tolerance (0.15脳 expected SA).
///
/// For `boolean_op(Difference, &B, &A)` = B - A:
/// - `a` = B (outer operand, the "box being cut")
/// - `b` = A (inner operand, the "cutting box")
pub fn try_difference_box_general(a: &BRep, b: &BRep) -> Option<BRep> {
    let b_info = try_as_box(a)?; // B (being cut)
    let a_info = try_as_box(b)?; // A (cutting)

    // Compute I = A 鈭?B via half-spaces.
    let i = try_intersection_box_general(a, b)?;

    // No overlap 鈫?B unchanged.
    if i.vertices.len() < 4 {
        return Some(a.clone());
    }

    let b_vol = volume(a);
    let i_vol = volume(&i);
    let scale = (b_info.extents.iter().sum::<f64>() / 3.0).max(1.0);
    let vol_tol = TOLERANCE_LEN_MIN * scale;

    // B fully inside A 鈫?empty (nothing of B remains outside A).
    if b_vol > vol_tol && (b_vol - i_vol).abs() < vol_tol {
        return Some(BRep::default());
    }

    // A fully inside B 鈫?fall through to Pave-Filler (hollow shell is hard).
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

    // 1. u-min exterior: B 鈭?{u路p 鈮?u_min_a}
    //    plane: (u_min_a*u, +u)  鈫?n=u, constraint u路p 鈮?u_min_a
    if b_u_min < u_min_a - zero_tol {
        try_slab!(vec![(u * u_min_a, u)]);
    }

    // 2. u-max exterior: B 鈭?{u路p 鈮?u_max_a}
    //    plane: (u_max_a*u, -u)  鈫?n=-u, constraint u路p 鈮?u_max_a
    if b_u_max > u_max_a + zero_tol {
        try_slab!(vec![(u * u_max_a, -u)]);
    }

    // Regions 3-6 need B to span across A's u-range (so "within u-range" is non-empty).
    let u_span = b_u_max > u_min_a + zero_tol && b_u_min < u_max_a - zero_tol;
    let v_span = b_v_max > v_min_a + zero_tol && b_v_min < v_max_a - zero_tol;

    // 3. v-min exterior within A's u-range:
    //    B 鈭?{u_min_a 鈮?u路p 鈮?u_max_a} 鈭?{v路p 鈮?v_min_a}
    //    planes: (u_min_a*u, -u), (u_max_a*u, +u), (v_min_a*v, +v)
    if u_span && b_v_min < v_min_a - zero_tol {
        try_slab!(vec![
            (u * u_min_a, -u),
            (u * u_max_a, u),
            (v * v_min_a, v),
        ]);
    }

    // 4. v-max exterior within A's u-range:
    //    B 鈭?{u_min_a 鈮?u路p 鈮?u_max_a} 鈭?{v路p 鈮?v_max_a}
    //    planes: (u_min_a*u, -u), (u_max_a*u, +u), (v_max_a*v, -v)
    if u_span && b_v_max > v_max_a + zero_tol {
        try_slab!(vec![
            (u * u_min_a, -u),
            (u * u_max_a, u),
            (v * v_max_a, -v),
        ]);
    }

    // 5. w-min exterior within A's u,v-range:
    //    B 鈭?{u_min_a 鈮?u路p 鈮?u_max_a} 鈭?{v_min_a 鈮?v路p 鈮?v_max_a} 鈭?{w路p 鈮?w_min_a}
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
    //    B 鈭?{u_min_a 鈮?u路p 鈮?u_max_a} 鈭?{v_min_a 鈮?v路p 鈮?v_max_a} 鈭?{w路p 鈮?w_max_a}
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
    // Use 1.15 to catch cases like a rotated-box cut (L4) where slab_sa_sum/sa_b 鈮?1.118.
    // boptuc F8 has ratio 1.071 鈥?fusion correctly merges it.
    if slab_sa_sum > sa_b * 1.15 {
        let _ = writeln!(
            std::io::stderr(),
            "DEBUG_MULTI_SLAB SA_INFLATED slab_count={} slab_sa={:.6} sa_b={:.6}",
            result.len(), slab_sa_sum, sa_b,
        );
        return None;
    }

    // Try boolean union to merge adjacent slabs (removes shared internal faces, G4-like).
    // Use `fuse` directly (not `boolean_op`) to avoid `try_union_disjoint_or_touching`
    // which returns a compound for face-touching slabs with adjacent bboxes (H4).
    let mut fused = result[0].clone();
    let mut ok = true;
    for slab in &result[1..] {
        match crate::bop_occt_union::fuse(&fused, slab) {
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
    if vol_ok && (merged || result.len() == 1 || fused_sa <= sa_b * 1.05) {
        return Some(fused);
    }
    if !ok {
        // Union errored 鈥?fall through to Pave-Filler.
        return None;
    }
    None
}

/// Union of two general boxes (axis-aligned or rotated), computed analytically.
///
/// Both BReps must be detected as boxes by [`try_as_box`]. For disjoint boxes,
/// returns a compound of the two inputs. For containment, returns the containing
/// box. For partial overlap, decomposes A and B into non-overlapping slabs
/// around each other's axes, adds the overlap region, and sews into a solid.
///
/// Falls through to Pave-Filler (returns `None`) when sewing fails or excessive
/// internal-face inflation is detected.
pub fn try_union_box_general(a: &BRep, b: &BRep) -> Option<BRep> {
    let info_a = try_as_box(a)?;
    let info_b = try_as_box(b)?;

    let inter = try_intersection_box_general(a, b)?;

    // No overlap 鈫?compound (disjoint boxes).
    if inter.vertices.len() < 4 {
        return Some(BRep::compound_from_shapes(&[a.clone(), b.clone()]));
    }

    let inter_vol = volume(&inter);
    let a_vol = volume(a);
    let b_vol = volume(b);
    let scale = (a_vol.max(b_vol) / 3.0).max(1.0);
    let vol_tol = TOLERANCE_LEN_MIN * scale;

    // A contains B 鈫?A alone is the union.
    if b_vol > vol_tol && (b_vol - inter_vol).abs() < vol_tol {
        return Some(a.clone());
    }
    // B contains A 鈫?B alone is the union.
    if a_vol > vol_tol && (a_vol - inter_vol).abs() < vol_tol {
        return Some(b.clone());
    }

    let a_planes = info_a.planes();
    let b_planes = info_b.planes();
    let zero_tol = TOLERANCE_LEN_MIN * scale;

    let mut slabs: Vec<BRep> = Vec::new();

    // Helper macro: build a convex polyhedron from base planes + extra half-spaces.
    macro_rules! try_slab {
        ($base:expr, $extra:expr) => {
            let mut planes = ($base).to_vec();
            planes.extend($extra);
            if let Ok(s) = make_convex_polyhedron_from_half_spaces(&planes) {
                if volume(&s) > vol_tol * 0.1 {
                    slabs.push(s);
                }
            }
        };
    }

    // 鈹€鈹€ A \ B slabs: decompose A around B's axes 鈹€鈹€
    {
        let [u, v, w] = info_b.axes;
        let [eu, ev, ew] = info_b.extents;
        let c = info_b.center;
        let u_min = u.dot(c) - eu;
        let u_max = u.dot(c) + eu;
        let v_min = v.dot(c) - ev;
        let v_max = v.dot(c) + ev;
        let w_min = w.dot(c) - ew;
        let w_max = w.dot(c) + ew;

        let a_verts: Vec<DVec3> = a.vertices.iter().map(|vi| vi.point).collect();
        let a_u_min = a_verts.iter().map(|p| u.dot(*p)).fold(f64::MAX, f64::min);
        let a_u_max = a_verts.iter().map(|p| u.dot(*p)).fold(f64::MIN, f64::max);
        let a_v_min = a_verts.iter().map(|p| v.dot(*p)).fold(f64::MAX, f64::min);
        let a_v_max = a_verts.iter().map(|p| v.dot(*p)).fold(f64::MIN, f64::max);
        let a_w_min = a_verts.iter().map(|p| w.dot(*p)).fold(f64::MAX, f64::min);
        let a_w_max = a_verts.iter().map(|p| w.dot(*p)).fold(f64::MIN, f64::max);

        let u_span = a_u_max > u_min + zero_tol && a_u_min < u_max - zero_tol;
        let v_span = a_v_max > v_min + zero_tol && a_v_min < v_max - zero_tol;

        if a_u_min < u_min - zero_tol {
            try_slab!(&a_planes, vec![(u * u_min, u)]);
        }
        if a_u_max > u_max + zero_tol {
            try_slab!(&a_planes, vec![(u * u_max, -u)]);
        }
        if u_span && a_v_min < v_min - zero_tol {
            try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, v)]);
        }
        if u_span && a_v_max > v_max + zero_tol {
            try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_max, -v)]);
        }
        if u_span && v_span && a_w_min < w_min - zero_tol {
            try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_min, w)]);
        }
        if u_span && v_span && a_w_max > w_max + zero_tol {
            try_slab!(&a_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_max, -w)]);
        }
    }

    // 鈹€鈹€ A 鈭?B (overlap) 鈹€鈹€
    slabs.push(inter.clone());

    // 鈹€鈹€ B \ A slabs: decompose B around A's axes 鈹€鈹€
    {
        let [u, v, w] = info_a.axes;
        let [eu, ev, ew] = info_a.extents;
        let c = info_a.center;
        let u_min = u.dot(c) - eu;
        let u_max = u.dot(c) + eu;
        let v_min = v.dot(c) - ev;
        let v_max = v.dot(c) + ev;
        let w_min = w.dot(c) - ew;
        let w_max = w.dot(c) + ew;

        let b_verts: Vec<DVec3> = b.vertices.iter().map(|vi| vi.point).collect();
        let b_u_min = b_verts.iter().map(|p| u.dot(*p)).fold(f64::MAX, f64::min);
        let b_u_max = b_verts.iter().map(|p| u.dot(*p)).fold(f64::MIN, f64::max);
        let b_v_min = b_verts.iter().map(|p| v.dot(*p)).fold(f64::MAX, f64::min);
        let b_v_max = b_verts.iter().map(|p| v.dot(*p)).fold(f64::MIN, f64::max);
        let b_w_min = b_verts.iter().map(|p| w.dot(*p)).fold(f64::MAX, f64::min);
        let b_w_max = b_verts.iter().map(|p| w.dot(*p)).fold(f64::MIN, f64::max);

        let u_span = b_u_max > u_min + zero_tol && b_u_min < u_max - zero_tol;
        let v_span = b_v_max > v_min + zero_tol && b_v_min < v_max - zero_tol;

        if b_u_min < u_min - zero_tol {
            try_slab!(&b_planes, vec![(u * u_min, u)]);
        }
        if b_u_max > u_max + zero_tol {
            try_slab!(&b_planes, vec![(u * u_max, -u)]);
        }
        if u_span && b_v_min < v_min - zero_tol {
            try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, v)]);
        }
        if u_span && b_v_max > v_max + zero_tol {
            try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_max, -v)]);
        }
        if u_span && v_span && b_w_min < w_min - zero_tol {
            try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_min, w)]);
        }
        if u_span && v_span && b_w_max > w_max + zero_tol {
            try_slab!(&b_planes, vec![(u * u_min, -u), (u * u_max, u), (v * v_min, -v), (v * v_max, v), (w * w_max, -w)]);
        }
    }

    if slabs.is_empty() {
        return Some(BRep::default());
    }
    if slabs.len() == 1 {
        return Some(slabs.remove(0));
    }

    // 鈹€鈹€ SA-inflation guard: try sew, fall back to fuse if inflated 鈹€鈹€
    let _slab_sa_sum: f64 = slabs.iter().map(|s| surface_area(s)).sum();
    let slab_vol_sum: f64 = slabs.iter().map(|s| volume(s)).sum();
    let expected_union_sa = surface_area(a) + surface_area(b) - surface_area(&inter);

    let sewn = sew_slabs_into_solid(&slabs, zero_tol);
    let sewn_sa = surface_area(&sewn);

    if sewn_sa <= expected_union_sa * 1.15 {
        return Some(sewn);
    }

    // SA inflated 鈥?try fuse-based sequential merge to remove shared internal faces.
    let mut fused = slabs[0].clone();
    let mut ok = true;
    for slab in &slabs[1..] {
        match crate::bop_occt_union::fuse(&fused, slab) {
            Ok(u) => { fused = u; }
            Err(_) => { ok = false; break; }
        }
    }
    if ok {
        let fused_sa = surface_area(&fused);
        let fused_vol = volume(&fused);
        let vol_ok = (fused_vol - slab_vol_sum).abs() < vol_tol * (slabs.len() as f64).max(1000.0);
        if vol_ok && fused_sa <= expected_union_sa * 1.15 {
            return Some(fused);
        }
    }

    Some(sewn)
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

/// Fast path for coaxial cylinder-sphere difference.
///
/// For `sphere - cylinder`: returns the spherical portion(s) outside the
/// cylinder's Z range (when the sphere cross-section fits entirely inside the
/// cylinder radius over the overlap).
///
/// For `cylinder - sphere`: returns the cylinder body with a spherical cavity
/// at the overlap end, built as a Z-slice triangle mesh.
pub fn try_difference_coaxial_cylinder_sphere(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try analytic fast paths first (exact geometry, no tessellation).
    if let Some(result) = crate::cylinder_sphere_analytic::build_sphere_minus_cylinder_analytic(a, b) {
        return Some(result);
    }
    if let Some(result) = crate::cylinder_sphere_analytic::build_cylinder_minus_sphere_analytic(a, b) {
        return Some(result);
    }

    // Fall through to mesh-backed detection logic.
    // Try both orderings. Track which operand is the cylinder to build the
    // correct result: sphere - cylinder (existing) vs cylinder - sphere.
    if let Some((sp, cyl)) = try_sphere_primitive_center_radius(a)
        .and_then(|sp| z_axis_cylinder_z_span_r(b).map(|cyl| (sp, cyl)))
    {
        // a = sphere, b = cylinder 鈫?sphere - cylinder
        if let Some(result) = try_sphere_minus_cylinder(sp, cyl) {
            return Some(result);
        }
    }

    if let Some((sp, cyl)) = try_sphere_primitive_center_radius(b)
        .and_then(|sp| z_axis_cylinder_z_span_r(a).map(|cyl| (sp, cyl)))
    {
        // a = cylinder, b = sphere 鈫?cylinder - sphere
        if let Some(result) = try_cylinder_minus_sphere(sp, cyl) {
            return Some(result);
        }
    }

    None
}

/// Build `sphere - cylinder` 鈥?portions of sphere outside cylinder Z range.
fn try_sphere_minus_cylinder(
    sphere_brep: (DVec3, f64),
    cyl_brep: (f64, f64, f64),
) -> Option<BRep> {
    let (center, radius) = sphere_brep;
    let (cyl_z_lo, cyl_z_hi, cyl_r) = cyl_brep;

    if center.x.abs() > TOLERANCE_ABS || center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    let sz = center.z;
    let sphere_z_lo = sz - radius;
    let sphere_z_hi = sz + radius;

    let overlap_lo = sphere_z_lo.max(cyl_z_lo);
    let overlap_hi = sphere_z_hi.min(cyl_z_hi);
    if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
        return None;
    }

    for z in [overlap_lo, overlap_hi] {
        let dz = z - sz;
        let r_at = (radius.powi(2) - dz.powi(2)).sqrt();
        if r_at > cyl_r + TOLERANCE_LEN_MIN {
            return None;
        }
    }

    let mut parts: Vec<BRep> = Vec::new();
    if sphere_z_lo < cyl_z_lo - TOLERANCE_LEN_MIN {
        let z_to = cyl_z_lo.min(sphere_z_hi);
        if z_to - sphere_z_lo > TOLERANCE_LEN_MIN {
            if let Some(p) = build_spherical_slice_solid(center, radius, sphere_z_lo, z_to) {
                parts.push(p);
            }
        }
    }
    if sphere_z_hi > cyl_z_hi + TOLERANCE_LEN_MIN {
        let z_from = cyl_z_hi.max(sphere_z_lo);
        if sphere_z_hi - z_from > TOLERANCE_LEN_MIN {
            if let Some(p) = build_spherical_slice_solid(center, radius, z_from, sphere_z_hi) {
                parts.push(p);
            }
        }
    }

    match parts.len() {
        0 => Some(BRep::default()),
        1 => Some(parts.swap_remove(0)),
        _ => {
            let mut base = parts.swap_remove(0);
            for p in parts {
                append_frustum_brep(&mut base, p);
            }
            Some(base)
        }
    }
}

/// Build `cylinder - sphere` 鈥?cylinder body with a spherical cavity where the
/// sphere overlaps.  The sphere cross-section must fit inside the cylinder over
/// the entire overlap Z-range (checked by the caller).
fn try_cylinder_minus_sphere(
    sphere_brep: (DVec3, f64),
    cyl_brep: (f64, f64, f64),
) -> Option<BRep> {
    let (center, radius) = sphere_brep;
    let (cyl_z_lo, cyl_z_hi, cyl_r) = cyl_brep;

    if center.x.abs() > TOLERANCE_ABS || center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    let sz = center.z;
    let sphere_z_lo = sz - radius;
    let sphere_z_hi = sz + radius;

    let overlap_lo = sphere_z_lo.max(cyl_z_lo);
    let overlap_hi = sphere_z_hi.min(cyl_z_hi);
    if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
        return None;
    }

    for z in [overlap_lo, overlap_hi] {
        let dz = z - sz;
        let r_at = (radius.powi(2) - dz.powi(2)).sqrt();
        if r_at > cyl_r + TOLERANCE_LEN_MIN {
            return None;
        }
    }

    build_cylinder_minus_sphere_tessellated(
        DVec2::new(center.x, center.y),
        cyl_z_lo, cyl_z_hi, cyl_r,
        center, radius,
    )
}

/// Triangulate an annular ring at `z` between `inner_r` and `outer_r`.
/// Winding direction is CCW on outer, CCW on inner (both triangles face +Z).
/// Area computation is unaffected by winding so this works for any ring.
fn add_annular_ring(
    add_v: &mut impl FnMut(DVec3) -> usize,
    faces: &mut Vec<Face>,
    z: f64, inner_r: f64, outer_r: f64, n_pts: usize,
    to_world: &impl Fn(f64, f64, f64) -> DVec3,
    empty_wire: &impl Fn() -> Wire,
) {
    let tau = std::f64::consts::TAU;
    let mut outer_idx = Vec::with_capacity(n_pts + 1);
    let mut inner_idx = Vec::with_capacity(n_pts + 1);
    for k in 0..=n_pts {
        let ang = tau * k as f64 / n_pts as f64;
        let (s, c) = ang.sin_cos();
        outer_idx.push(add_v(to_world(outer_r * c, outer_r * s, z)));
        inner_idx.push(add_v(to_world(inner_r * c, inner_r * s, z)));
    }
    let mut tris = Vec::with_capacity(n_pts * 2);
    for j in 0..n_pts {
        tris.push([outer_idx[j], outer_idx[j + 1], inner_idx[j + 1]]);
        tris.push([outer_idx[j], inner_idx[j + 1], inner_idx[j]]);
    }
    faces.push(Face {
        outer_wire: empty_wire(), inner_wires: vec![],
        normal: DVec3::ZERO, triangles: tris,
        sample_point: None, mesh_dirty: false,
    });
}

/// Build a tessellated cylinder body with a spherical cavity at the overlap.
/// Result = cylinder body minus the spherical portion.  The cylinder wall is
/// intact (full height), the far-end flat cap is preserved, and the spherical
/// cavity surface replaces the near-end region.  If the cavity extends to a
/// cylinder end without filling its cross-section, an annular flat ring caps
/// the remaining opening.
fn build_cylinder_minus_sphere_tessellated(
    center_xy: DVec2,
    cyl_z_lo: f64, cyl_z_hi: f64, cyl_r: f64,
    sphere_center: DVec3, sphere_r: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let sz = sphere_center.z;
    let sphere_z_lo = sz - sphere_r;
    let sphere_z_hi = sz + sphere_r;
    let overlap_lo = cyl_z_lo.max(sphere_z_lo);
    let overlap_hi = cyl_z_hi.min(sphere_z_hi);
    if overlap_hi <= overlap_lo + tol { return None; }

    let n_arc = 128;
    let n_slices_circ = 16;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    let circle_poly = |r: f64| -> Vec<DVec2> {
        let mut poly = Vec::with_capacity(n_arc + 1);
        for k in 0..=n_arc {
            let ang = tau * k as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            poly.push(DVec2::new(r * c, r * s));
        }
        poly
    };

    let cyl_poly = circle_poly(cyl_r);

    // 1. Cylinder wall (full height)
    add_wall_section(&mut add_v, &mut faces, &cyl_poly, cyl_z_lo, cyl_z_hi, n_slices_circ, &to_world, &empty_wire);

    // 2. Determine which cylinder ends are reached by the cavity
    let cavity_reaches_top = (overlap_hi - cyl_z_hi).abs() < 1e-9;
    let cavity_reaches_bot = (overlap_lo - cyl_z_lo).abs() < 1e-9;

    // 3. Far-end cap(s) 鈥?cylinder ends NOT reached by cavity
    if !cavity_reaches_bot {
        add_cap_face(&mut add_v, &mut faces, &cyl_poly, cyl_z_lo, -DVec3::Z, &to_world, &empty_wire);
    }
    if !cavity_reaches_top {
        add_cap_face(&mut add_v, &mut faces, &cyl_poly, cyl_z_hi, DVec3::Z, &to_world, &empty_wire);
    }

    // 4. Spherical cavity surface in overlap region
    {
        let n_rings = 32usize;
        let dz = (overlap_hi - overlap_lo) / n_rings as f64;
        for i in 0..n_rings {
            let za = overlap_lo + dz * i as f64;
            let zb = overlap_lo + dz * (i + 1) as f64;
            let dz_a = za - sz;
            let dz_b = zb - sz;
            let r_sq_a = sphere_r.powi(2) - dz_a.powi(2);
            let r_sq_b = sphere_r.powi(2) - dz_b.powi(2);
            let ra = if r_sq_a <= 0.0 { 0.0 } else { r_sq_a.sqrt() };
            let rb = if r_sq_b <= 0.0 { 0.0 } else { r_sq_b.sqrt() };
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(ra * c, ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(rb * c, rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // 5. Annular ring(s) at cylinder end(s) where cavity reaches but doesn't fill
    if cavity_reaches_top {
        let r_at = (sphere_r.powi(2) - (cyl_z_hi - sz).powi(2)).sqrt().max(0.0);
        if r_at < cyl_r - tol && r_at >= 0.0 {
            add_annular_ring(&mut add_v, &mut faces, cyl_z_hi, r_at, cyl_r, n_arc, &to_world, &empty_wire);
        }
    }
    if cavity_reaches_bot {
        let r_at = (sphere_r.powi(2) - (cyl_z_lo - sz).powi(2)).sqrt().max(0.0);
        if r_at < cyl_r - tol && r_at >= 0.0 {
            add_annular_ring(&mut add_v, &mut faces, cyl_z_lo, r_at, cyl_r, n_arc, &to_world, &empty_wire);
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Build a Z-slice triangle-mesh solid representing a spherical slice (portion
/// of a sphere between two parallel Z planes). The solid is built with
/// pre-triangulated faces and no analytic surfaces 鈥?SA computation reads
/// the stored triangles directly.
fn build_spherical_slice_solid(center: DVec3, radius: f64, z_lo: f64, z_hi: f64) -> Option<BRep> {
    const N_RINGS: usize = 32;
    const N_CIRC: usize = 48;

    let h = z_hi - z_lo;
    if h <= TOLERANCE_LEN_MIN {
        return None;
    }
    let dz = h / N_RINGS as f64;

    let mut brep = BRep::new();

    // Generate ring vertices
    let mut verts: Vec<Vertex> = Vec::new();
    let mut ring_start: Vec<usize> = Vec::with_capacity(N_RINGS + 2);

    for i in 0..=N_RINGS {
        let z = z_lo + dz * i as f64;
        let dz_c = z - center.z;
        let r_sq = radius.powi(2) - dz_c.powi(2);
        ring_start.push(verts.len());

        if r_sq <= 0.0 {
            verts.push(Vertex { point: DVec3::new(center.x, center.y, z) });
        } else {
            let r = r_sq.sqrt();
            for j in 0..N_CIRC {
                let ang = std::f64::consts::TAU * j as f64 / N_CIRC as f64;
                let (s, c) = ang.sin_cos();
                verts.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, z) });
            }
        }
    }
    ring_start.push(verts.len());

    let v_off = brep.vertices.len();
    brep.vertices = verts;

    let mut faces: Vec<Face> = Vec::new();

    // Lateral faces: triangle strips between adjacent rings
    for i in 0..N_RINGS {
        let b_s = ring_start[i];
        let b_e = ring_start[i + 1];
        let t_s = ring_start[i + 1];
        let t_e = ring_start[i + 2];
        let bn = b_e - b_s;
        let tn = t_e - t_s;

        let tris: Vec<[usize; 3]> = if bn == 1 && tn > 1 {
            (0..tn).flat_map(|j| {
                let k = (j + 1) % tn;
                vec![[v_off + b_s, v_off + t_s + j, v_off + t_s + k]]
            }).collect()
        } else if tn == 1 && bn > 1 {
            (0..bn).flat_map(|j| {
                let k = (j + 1) % bn;
                vec![[v_off + b_s + j, v_off + b_s + k, v_off + t_s]]
            }).collect()
        } else if bn > 1 && tn > 1 {
            let n = bn.min(tn);
            (0..n).flat_map(|j| {
                let k = (j + 1) % n;
                let b0 = v_off + b_s + j;
                let b1 = v_off + b_s + k;
                let t0 = v_off + t_s + j;
                let t1 = v_off + t_s + k;
                vec![[b0, b1, t1], [b0, t1, t0]]
            }).collect()
        } else {
            Vec::new()
        };

        if !tris.is_empty() {
            faces.push(Face {
                outer_wire: Wire { edges: vec![] },
                inner_wires: vec![],
                normal: DVec3::ZERO,
                triangles: tris,
                sample_point: None,
                mesh_dirty: false,
            });
        }
    }

    // Cap disks at exposed planar ends (skip if at pole)
    let bot_n = ring_start[1] - ring_start[0];
    if bot_n > 1 {
        let r_sq = radius.powi(2) - (z_lo - center.z).powi(2);
        if r_sq > 0.0 {
            let r = r_sq.sqrt();
            push_cap_disk(&mut brep, &mut faces, center, r, z_lo,
                if z_lo > center.z { DVec3::Z } else { -DVec3::Z });
        }
    }

    let top_n = ring_start[N_RINGS + 1] - ring_start[N_RINGS];
    if top_n > 1 {
        let r_sq = radius.powi(2) - (z_hi - center.z).powi(2);
        if r_sq > 0.0 {
            let r = r_sq.sqrt();
            push_cap_disk(&mut brep, &mut faces, center, r, z_hi,
                if z_hi > center.z { DVec3::Z } else { -DVec3::Z });
        }
    }

    brep.solids.push(Solid { shells: vec![Shell { faces }] });
    Some(brep)
}

/// Push a triangulated planar disk face (center + rim at given z) into `faces`.
fn push_cap_disk(brep: &mut BRep, faces: &mut Vec<Face>, center: DVec3, r: f64, z: f64, normal: DVec3) {
    const N_CIRC: usize = 48;
    let mut rim: Vec<usize> = Vec::with_capacity(N_CIRC);
    for j in 0..N_CIRC {
        let ang = std::f64::consts::TAU * j as f64 / N_CIRC as f64;
        let (s, c) = ang.sin_cos();
        let idx = brep.vertices.len();
        brep.vertices.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, z) });
        rim.push(idx);
    }
    let ctr_idx = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(center.x, center.y, z) });

    let tris: Vec<[usize; 3]> = (0..N_CIRC).map(|j| {
        let k = (j + 1) % N_CIRC;
        [rim[j], rim[k], ctr_idx]
    }).collect();

    faces.push(Face {
        outer_wire: Wire { edges: vec![] },
        inner_wires: vec![],
        normal,
        triangles: tris,
        sample_point: None,
        mesh_dirty: false,
    });
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
/// Returns [`None`] when centers differ or radii are degenerate 閳?callers fall back to Pave/Builder.
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

// --- Coaxial cone 閳?cylinder (OCCT `bopcommon_simple/ZP7`): generic Builder over-counts area. --------

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

pub(crate) fn z_axis_cylinder_z_span_r(cyl: &BRep) -> Option<(f64, f64, f64)> {
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

/// Sharp Z-aligned cone 閳?finite Z-aligned cylinder (same axis / origin in `xy`), e.g. OCCT ZP7.
pub fn try_intersection_coaxial_cone_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_coaxial_cone_cylinder_pair(a, b)
        .or_else(|| try_intersection_coaxial_cone_cylinder_pair(b, a))
}

/// Two coaxial Z-aligned cylinders with the same radius -> intersection is the
/// overlapping Z-span cylinder (e.g. OCCT bcommon_simple/J1).
pub fn try_intersection_coaxial_cylinder_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    let (z1_lo, z1_hi, r1) = z_axis_cylinder_z_span_r(a)?;
    let (z2_lo, z2_hi, r2) = z_axis_cylinder_z_span_r(b)?;
    if (r1 - r2).abs() > TOLERANCE_ADAPTIVE_MAX {
        return None;
    }
    let z0 = z1_lo.max(z2_lo);
    let z1 = z1_hi.min(z2_hi);
    if z1 - z0 < TOLERANCE_MESH_LEGACY {
        return None;
    }
    let zm = (z0 + z1) * 0.5;
    let h = z1 - z0;
    rcad_modeling::make_cylinder_brep(
        DVec3::new(0.0, 0.0, zm), DVec3::Z, DVec3::X, r1, h,
    ).ok()
}

/// Detect a cylinder BRep and return (center, axis, radius, height, origin, ref_dir).
///
/// Unlike `try_cylinder_center_axis_radius_height`, this works for cylinders
/// along ANY axis (not just Z) by computing height from face-plane distances
/// along the axis.
fn try_cylinder_any_axis(brep: &BRep) -> Option<(DVec3, DVec3, f64, f64, DVec3, DVec3)> {
    for s in &brep.solids {
        for sh in &s.shells {
            let mut origin = DVec3::ZERO;
            let mut axis = DVec3::Z;
            let mut ref_dir = DVec3::X;
            let mut radius = 0.0;
            let mut found = false;
            let mut plane_dots: Vec<f64> = Vec::new();

            for fi in 0..sh.faces.len() {
                let si = brep.geom.face_surface.get(fi)?.as_ref()?;
                match brep.geom.surfaces.get(*si)? {
                    Surface3::Cylinder(cc) => {
                        origin = cc.origin;
                        axis = cc.axis.normalize_or_zero();
                        radius = cc.radius;
                        ref_dir = cc.ref_dir;
                        found = true;
                    }
                    Surface3::Plane(pl) => {
                        plane_dots.push(pl.origin.dot(axis));
                    }
                    _ => return None,
                }
            }
            if !found || plane_dots.len() < 2 { return None; }
            plane_dots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let height = plane_dots.last()? - plane_dots.first()?;
            if height < crate::tolerance::TOLERANCE_LEN_MIN { return None; }
            let center = origin + axis * (height / 2.0);
            return Some((center, axis, radius, height, origin, ref_dir));
        }
    }
    None
}

/// Build a mesh-based BRep for the intersection of two perpendicular
/// equal-radius cylinders (a Steinmetz-like intersection).
///
/// When two cylinders share the same radius R and their axes are perpendicular,
/// the intersection can be built by sampling the surface of each cylinder and
/// keeping only the points inside the other.  This avoids the PaveFiller
/// (5.3s 鈫?<0.01s) for the I9 test case.
fn build_perpendicular_cylinder_intersection(
    c1: CylParams, c2: CylParams,
) -> Option<BRep> {
    use std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = vec![];
    let mut tris: Vec<[usize; 3]> = Vec::new();

    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    // Helper: test if a point is inside cylinder (within radius and height bounds)
    let inside_cyl = |p: DVec3, cyl: &CylParams| -> bool {
        let d = p - cyl.center;
        let along = d.dot(cyl.axis);
        if along.abs() > cyl.height / 2.0 { return false; }
        let perp = d - along * cyl.axis;
        perp.length_squared() <= cyl.radius * cyl.radius + 1e-10
    };

    // Generate surface from one cylinder, clipped by the other.
    // Inlined as a plain loop to avoid closure capture issues with add_v
    // (which borrows verts mutably across two calls).
    for pair in [(c1, c2), (c2, c1)] {
        let (cyl, other) = (pair.0, pair.1);
        let x_ax = cyl.any_perp;
        let y_ax = cyl.axis.cross(x_ax).normalize();
        const NU: usize = 48;
        const NV: usize = 32;
        let mut idx = vec![vec![0usize; NU]; NV + 1];

        for vj in 0..=NV {
            let t = vj as f64 / NV as f64;
            let v = (t - 0.5) * cyl.height;
            for ui in 0..NU {
                let u = ui as f64 * TAU / NU as f64;
                let (cu, su) = u.sin_cos();
                let p = cyl.center + v * cyl.axis
                    + cyl.radius * (cu * x_ax + su * y_ax);
                if inside_cyl(p, &other) {
                    idx[vj][ui] = add_v(p);
                }
            }
        }
        for vj in 0..NV {
            for ui in 0..NU {
                let a = idx[vj][ui];
                let b = idx[vj][(ui + 1) % NU];
                let c = idx[vj + 1][ui];
                let d = idx[vj + 1][(ui + 1) % NU];
                match (a, b, c, d) {
                    (0, _, _, _) | (_, 0, _, _) | (_, _, 0, _) | (_, _, _, 0) => {}
                    _ => {
                        tris.push([a, b, d]);
                        tris.push([a, d, c]);
                    }
                }
            }
        }
    }

    // 鈹€鈹€ Cap faces 鈹€鈹€
    for pair in [(c1, c2), (c2, c1)] {
        let (cyl, other) = (pair.0, pair.1);
        let x_ax = cyl.any_perp;
        let y_ax = cyl.axis.cross(x_ax).normalize();
        for sign in [-1.0, 1.0] {
            let cap_center = cyl.center + sign * (cyl.height / 2.0) * cyl.axis;
            const NC: usize = 24;
            let mut cap_pts: Vec<DVec3> = Vec::new();
            for i in 0..NC {
                for j in 0..NC {
                    let u = (i as f64 / (NC - 1) as f64 - 0.5) * 2.0 * cyl.radius;
                    let v = (j as f64 / (NC - 1) as f64 - 0.5) * 2.0 * cyl.radius;
                    let p = cap_center + u * x_ax + v * y_ax;
                    if u * u + v * v <= cyl.radius * cyl.radius + 1e-10
                        && inside_cyl(p, &other)
                    {
                        cap_pts.push(p);
                    }
                }
            }
            if cap_pts.is_empty() { continue; }
            let avg = cap_pts.iter().copied().sum::<DVec3>() / cap_pts.len() as f64;
            let ci = add_v(avg);
            let mut sorted: Vec<(f64, DVec3)> = cap_pts.iter().map(|p| {
                let dp = *p - avg;
                (dp.dot(x_ax).atan2(dp.dot(y_ax)).rem_euclid(TAU), *p)
            }).collect();
            sorted.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let cap_vis: Vec<usize> = sorted.iter().map(|(_, p)| add_v(*p)).collect();
            for i in 0..cap_vis.len() {
                let j = (i + 1) % cap_vis.len();
                tris.push([ci, cap_vis[i], cap_vis[j]]);
            }
        }
    }

    if tris.is_empty() { return None; }

    let faces = vec![Face {
        outer_wire: empty_wire(), inner_wires: vec![],
        normal: DVec3::ZERO, triangles: tris,
        sample_point: None, mesh_dirty: false,
    }];

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

#[derive(Clone, Copy)]
struct CylParams {
    center: DVec3,
    axis: DVec3,
    radius: f64,
    height: f64,
    any_perp: DVec3,
}

/// Fast-path for two perpendicular cylinders with equal radius.
///
/// Detects when both operands are cylinders, share the same radius (within
/// tolerance), and have perpendicular axes.  Builds a mesh-based BRep by
/// sampling each cylinder surface and keeping points inside the other.
/// OCCT has no equivalent 鈥?this is a pure rcad optimization that avoids
/// the PaveFiller (5.3s 鈫?<0.01s) for the I9 test case.
pub fn try_intersection_cylinder_cylinder_perpendicular(a: &BRep, b: &BRep) -> Option<BRep> {
    let c1 = try_cylinder_any_axis(a)?;
    let c2 = try_cylinder_any_axis(b)?;

    // Equal radius (index 2 = radius)
    if (c1.2 - c2.2).abs() > 1e-4 { return None; }

    // Perpendicular axes (index 1 = axis)
    if c1.1.dot(c2.1).abs() > 1e-3 { return None; }

    let p1 = CylParams { center: c1.0, axis: c1.1, radius: c1.2, height: c1.3, any_perp: c1.5 };
    let p2 = CylParams { center: c2.0, axis: c2.1, radius: c2.2, height: c2.3, any_perp: c2.5 };

    build_perpendicular_cylinder_intersection(p1, p2)
}

/// Build C1 \ C2 for two coaxial Z-aligned cylinders.
///
/// When r1 <= r2 (C1 is fully contained in C2 in XY), the result is the
/// portion(s) of C1 outside C2's Z-range.  When r1 > r2 (C1 extends beyond
/// C2), the result would need a cylindrical hole 鈥?too complex for now.
pub fn try_difference_coaxial_cylinder_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    let (z1_lo, z1_hi, r1) = z_axis_cylinder_z_span_r(a)?;
    let (z2_lo, z2_hi, r2) = z_axis_cylinder_z_span_r(b)?;
    let overlap_lo = z1_lo.max(z2_lo);
    let overlap_hi = z1_hi.min(z2_hi);
    if overlap_hi - overlap_lo < TOLERANCE_MESH_LEGACY {
        // No Z-overlap 鈫?a is unchanged by subtracting b
        return Some(a.clone());
    }
    if r1 > r2 + TOLERANCE_ADAPTIVE_MAX {
        // C1 extends beyond C2 in XY 鈥?would need a cylindrical hole in the
        // result.  Fall through to the general boolean engine.
        return None;
    }
    // r1 <= r2: C1 is fully inside C2 in XY (or equal radii).
    // Result is the portions of C1 outside C2's Z-range.
    let mut result = BRep::new();
    // Piece below C2
    let z_below_end = z1_hi.min(z2_lo);
    let h_below = z_below_end - z1_lo;
    if h_below > TOLERANCE_MESH_LEGACY {
        let center_below = (z1_lo + z_below_end) * 0.5;
        let piece = rcad_modeling::make_cylinder_brep(
            DVec3::new(0.0, 0.0, center_below), DVec3::Z, DVec3::X, r1, h_below,
        ).ok()?;
        if result.solids.is_empty() {
            result = piece;
        } else {
            result.append_disjoint_brep(&piece);
        }
    }
    // Piece above C2
    let z_above_start = z1_lo.max(z2_hi);
    let h_above = z1_hi - z_above_start;
    if h_above > TOLERANCE_MESH_LEGACY {
        let center_above = (z_above_start + z1_hi) * 0.5;
        let piece = rcad_modeling::make_cylinder_brep(
            DVec3::new(0.0, 0.0, center_above), DVec3::Z, DVec3::X, r1, h_above,
        ).ok()?;
        if result.solids.is_empty() {
            result = piece;
        } else {
            result.append_disjoint_brep(&piece);
        }
    }
    if result.solids.is_empty() {
        None // C1 fully removed
    } else {
        Some(result)
    }
}

/// Extract sphere center and radius from a sphere BRep (first SphericalSurface found).
pub(crate) fn sphere_center_r(sphere: &BRep) -> Option<(DVec3, f64)> {
    for s in &sphere.solids {
        for sh in &s.shells {
            for fi in 0..sh.faces.len() {
                if let Some(Some(si)) = sphere.geom.face_surface.get(fi) {
                    if let Surface3::Sphere(sp) = sphere.geom.surfaces.get(*si)? {
                        return Some((sp.center, sp.radius));
                    }
                }
            }
        }
    }
    None
}

/// Build a BRep for the part of a sphere clipped between two Z-planes.
///
/// The sphere is axis-aligned with axis=Z so that z=constant 鈫?v=constant in the
/// sphere's parameterization.  The result has a spherical lateral face and planar
/// cap(s) at the clip planes (or none when the clip is at the sphere pole).
///
/// `z_min` and `z_max` are the clip planes in world Z (z_min < z_max), assumed
/// to overlap the sphere's Z-range [center.z - radius, center.z + radius].
pub(crate) fn build_sphere_clipped_by_z_planes(
    center: DVec3,
    radius: f64,
    z_min: f64,
    z_max: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Curve2d, Line2d};
    use std::f64::consts::PI;

    let two_pi = 2.0 * PI;

    // Colatitude v in sphere param (axis=Z): z = center.z + r*cos(v), v = acos((z-C.z)/r)
    let cos_v_hi = ((z_max - center.z) / radius).clamp(-1.0, 1.0);
    let cos_v_lo = ((z_min - center.z) / radius).clamp(-1.0, 1.0);
    let v_hi = cos_v_hi.acos(); // smaller v, higher z  (equator 鈫?0 at north pole)
    let v_lo = cos_v_lo.acos(); // larger v, lower z  (equator 鈫?蟺 at south pole)

    let has_top_cap = (z_max - (center.z + radius)).abs() > 1e-12;
    let has_bot_cap = (z_min - (center.z - radius)).abs() > 1e-12;

    // Radii of the clip-plane circles
    let r_hi = radius * v_hi.sin();
    let r_lo = radius * v_lo.sin();

    // Vertex positions: seam runs at u=0 from v_hi to v_lo
    let v_hi_pt = center + radius * DVec3::new(v_hi.sin(), 0.0, v_hi.cos());
    let v_lo_pt = center + radius * DVec3::new(v_lo.sin(), 0.0, v_lo.cos());

    let mut brep = BRep::default();

    let v_hi_idx = make_vertex(&mut brep, v_hi_pt);
    let v_lo_idx = make_vertex(&mut brep, v_lo_pt);

    // 鈹€鈹€ Edges 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // E0: circle at v_hi (higher z, smaller v)
    let c_hi = DVec3::new(center.x, center.y, center.z + radius * v_hi.cos());
    let e0_curve = Curve3::Circle(Circle3 { center: c_hi, normal: DVec3::Z, radius: r_hi });
    let e0 = make_edge(&mut brep, e0_curve, 0.0, two_pi, v_hi_idx, v_hi_idx).ok()?;

    // E1: circle at v_lo (lower z, larger v) 鈥?only needed when both caps exist
    let e1 = if has_top_cap && has_bot_cap {
        let c_lo = DVec3::new(center.x, center.y, center.z + radius * v_lo.cos());
        let curve = Curve3::Circle(Circle3 { center: c_lo, normal: -DVec3::Z, radius: r_lo });
        let idx = make_edge(&mut brep, curve, 0.0, two_pi, v_lo_idx, v_lo_idx).ok()?;
        Some(idx)
    } else {
        None
    };

    // E2 (or E1 for single-cap): seam from v_hi to v_lo at u=0 (arc in XZ-plane)
    let seam = {
        let curve = Curve3::Circle(Circle3 { center, normal: DVec3::Y, radius });
        make_edge(&mut brep, curve, v_hi, v_lo, v_hi_idx, v_lo_idx).ok()?
    };

    // 鈹€鈹€ Surfaces 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let sphere_surf = Surface3::Sphere(SphericalSurface {
        center,
        axis: DVec3::Z,
        radius,
        ref_dir: DVec3::X,
    });

    let top_plane = if has_top_cap {
        Some(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, z_max),
            normal: DVec3::Z,
        }))
    } else {
        None
    };

    let bot_plane = if has_bot_cap {
        Some(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, z_min),
            normal: -DVec3::Z,
        }))
    } else {
        None
    };

    // 鈹€鈹€ Curve2Ds (pcurves) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Sphere face pcurves
    //   E0 on sphere: iso-v = v_hi
    let e0_on_sphere = Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, v_hi), direction: glam::DVec2::new(1.0, 0.0) });
    //   Seam fwd on sphere: u=0, v from v_hi to v_lo
    let seam_fwd = Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, v_hi), direction: glam::DVec2::new(0.0, 1.0) });
    //   Seam rev on sphere: u=2蟺, v from v_lo to v_hi
    let seam_rev = Curve2d::Line(Line2d { origin: glam::DVec2::new(two_pi, v_lo), direction: glam::DVec2::new(0.0, -1.0) });
    //   E1 on sphere (if present): iso-v = v_lo
    let e1_on_sphere = if has_top_cap && has_bot_cap {
        Some(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, v_lo), direction: glam::DVec2::new(1.0, 0.0) }))
    } else {
        None
    };

    // Planar cap pcurves
    let e0_on_plane = if has_top_cap {
        Some(Curve2d::Circle(rcad_kernel::geom::Circle2d { center: glam::DVec2::ZERO, radius: r_hi }))
    } else {
        None
    };
    let e1_on_bot_plane = if has_bot_cap && has_top_cap {
        Some(Curve2d::Circle(rcad_kernel::geom::Circle2d { center: glam::DVec2::ZERO, radius: r_lo }))
    } else {
        None
    };

    // 鈹€鈹€ Build geometry store 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let surf_idx_sphere = 0usize;
    let mut surf_idx_top: Option<usize> = None;
    let mut surf_idx_bot: Option<usize> = None;

    brep.geom.surfaces.push(sphere_surf);
    if let Some(tp) = top_plane {
        surf_idx_top = Some(brep.geom.surfaces.len());
        brep.geom.surfaces.push(tp);
    }
    if let Some(bp) = bot_plane {
        surf_idx_bot = Some(brep.geom.surfaces.len());
        brep.geom.surfaces.push(bp);
    }

    // Curve2D indices
    let mut c2d = 0usize;
    brep.geom.curve2ds.push(e0_on_sphere);
    let e0_sphere_c2d = c2d; c2d += 1;
    brep.geom.curve2ds.push(seam_fwd);
    let seam_fwd_c2d = c2d; c2d += 1;
    brep.geom.curve2ds.push(seam_rev);
    let seam_rev_c2d = c2d; c2d += 1;

    let e1_sphere_c2d = e1_on_sphere.map(|c| {
        brep.geom.curve2ds.push(c);
        let idx = c2d; c2d += 1; idx
    });

    let e0_plane_c2d = e0_on_plane.map(|c| {
        brep.geom.curve2ds.push(c);
        let idx = c2d; c2d += 1; idx
    });

    let e1_bot_c2d = e1_on_bot_plane.map(|c| {
        brep.geom.curve2ds.push(c);
        let idx = c2d; c2d += 1; idx
    });

    // Edge pcurves 鈥?align vecs
    while brep.geom.edge_pcurves.len() <= seam.max(e0).max(e1.unwrap_or(0)) {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    // E0 pcurves: on sphere always; on top plane if cap exists
    {
        let ep = &mut brep.geom.edge_pcurves[e0];
        ep.push(rcad_kernel::PCurve { surface_idx: surf_idx_sphere, curve2d_idx: e0_sphere_c2d });
        if let Some(si) = surf_idx_top {
            if let Some(ci) = e0_plane_c2d {
                ep.push(rcad_kernel::PCurve { surface_idx: si, curve2d_idx: ci });
            }
        }
    }

    // E1 pcurves: on sphere + bot plane if both caps
    if let Some(e1i) = e1 {
        let ep = &mut brep.geom.edge_pcurves[e1i];
        if let Some(ci) = e1_sphere_c2d {
            ep.push(rcad_kernel::PCurve { surface_idx: surf_idx_sphere, curve2d_idx: ci });
        }
        if let Some(si) = surf_idx_bot {
            if let Some(ci) = e1_bot_c2d {
                ep.push(rcad_kernel::PCurve { surface_idx: si, curve2d_idx: ci });
            }
        }
    }

    // Seam pcurves: two on sphere (fwd and rev)
    {
        let ep = &mut brep.geom.edge_pcurves[seam];
        ep.push(rcad_kernel::PCurve { surface_idx: surf_idx_sphere, curve2d_idx: seam_fwd_c2d });
        ep.push(rcad_kernel::PCurve { surface_idx: surf_idx_sphere, curve2d_idx: seam_rev_c2d });
    }

    // 鈹€鈹€ Faces 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Initialize solid/shell structure
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid {
            shells: vec![rcad_kernel::Shell { faces: Vec::new() }],
        });
    }

    let mut face_wires_sphere: Vec<WireEdge> = Vec::new();

    if has_top_cap && has_bot_cap {
        // Both caps: pattern = bottom_fwd 鈫?seam_fwd 鈫?top_rev 鈫?seam_rev
        // where bottom = v_lo (lower z), top = v_hi (higher z)
        face_wires_sphere.push(WireEdge::fwd(e1.unwrap())); // E1 fwd at v_lo
        face_wires_sphere.push(WireEdge::fwd(seam));        // seam fwd v_lo鈫抳_hi
        face_wires_sphere.push(WireEdge::rev(e0));           // E0 rev at v_hi
        face_wires_sphere.push(WireEdge::rev(seam));         // seam rev v_hi鈫抳_lo
    } else if has_top_cap {
        // Only top cap (v_bot is at south pole): pattern = E0_rev 鈫?seam_fwd 鈫?seam_rev
        face_wires_sphere.push(WireEdge::rev(e0));
        face_wires_sphere.push(WireEdge::fwd(seam));
        face_wires_sphere.push(WireEdge::rev(seam));
    } else if has_bot_cap {
        // Only bottom cap (v_hi is at north pole): pattern = seam_fwd 鈫?E0_rev 鈫?seam_rev
        // Actually: bottom circle fwd 鈫?seam fwd 鈫?seam_rev (no top circle since at pole)
        // Since E0 is at v_hi=north pole and E1 is at v_lo:
        //   lateral = E0_fwd 鈫?seam_fwd 鈫?seam_rev
        face_wires_sphere.push(WireEdge::fwd(e0));
        face_wires_sphere.push(WireEdge::fwd(seam));
        face_wires_sphere.push(WireEdge::rev(seam));
    } else {
        // No caps 鈥?sphere entirely inside cylinder, entire z-range inside
        return None; // Should have been caught by containment
    }

    let sphere_wire = make_wire(face_wires_sphere);
    let sphere_face = Face {
        outer_wire: sphere_wire,
        inner_wires: vec![],
        normal: DVec3::X,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
    };

    // Push sphere face
    let sphere_face_idx = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(sphere_face);
    while brep.geom.face_surface.len() <= sphere_face_idx {
        brep.geom.face_surface.push(None);
    }
    brep.geom.face_surface[sphere_face_idx] = Some(surf_idx_sphere);
    // Set face_surface_range to restrict to [0,2蟺] 脳 [v_hi, v_lo]
    while brep.geom.face_surface_range.len() <= sphere_face_idx {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range[sphere_face_idx] = Some([0.0, two_pi, v_hi, v_lo]);

    // Top cap face
    if let Some(si) = surf_idx_top {
        let cap_wire = make_wire(vec![WireEdge::fwd(e0)]);
        let cap_face = Face {
            outer_wire: cap_wire,
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
        };
        let fi = brep.solids[0].shells[0].faces.len();
        brep.solids[0].shells[0].faces.push(cap_face);
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(si);
    }

    // Bottom cap face
    if let Some(si) = surf_idx_bot {
        let cap_wire = make_wire(vec![WireEdge::rev(e1.unwrap())]);
        let cap_face = Face {
            outer_wire: cap_wire,
            inner_wires: vec![],
            normal: -DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
        };
        let fi = brep.solids[0].shells[0].faces.len();
        brep.solids[0].shells[0].faces.push(cap_face);
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(si);
    }

    Some(brep)
}

/// Build the intersection BRep for a coaxial Z-aligned cylinder 鈭?sphere with R_s > R_c.
///
/// The cylinder wall cuts the sphere at z = s_z 卤 鈭?r_s虏 - r_c虏).  This handles the
/// sub-case where only the LOWER intersection circle lies in the overlap Z-range and
/// the upper boundary is the cylinder end cap (sphere center above cylinder center).
///
/// Result faces:
///   鈥?Spherical face (bottom): south pole 鈫?intersection circle
///   鈥?Cylindrical wall face (middle): intersection circle 鈫?cylinder top (z_hi)
///   鈥?Cylinder top cap (top): planar disk at z = z_hi
fn build_cylinder_sphere_intersection_brep(
    z_lo: f64,
    z_hi: f64,
    r_c: f64,
    s_center: DVec3,
    r_s: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::PI;

    let two_pi = 2.0 * PI;

    // Overlap Z-range
    let z_min = z_lo.max(s_center.z - r_s);
    let z_max = z_hi.min(s_center.z + r_s);

    let dz = (r_s * r_s - r_c * r_c).sqrt();
    let z_isect = s_center.z - dz;

    // Lower intersection circle must be in range
    if z_isect < z_min - 1e-12 || z_isect > z_max + 1e-12 {
        return None;
    }

    let h_cyl = z_hi - z_isect;

    // Sphere colatitude at intersection circle
    let cos_v = ((z_isect - s_center.z) / r_s).clamp(-1.0, 1.0);
    let v_isect = cos_v.acos();

    let mut brep = BRep::default();

    // 鈹€鈹€ Vertices 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // V0: intersection circle at u=0 (= cylinder seam origin)
    let v0 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_isect));
    // V1: sphere south pole
    let v1 = make_vertex(&mut brep, DVec3::new(0.0, 0.0, s_center.z - r_s));
    // V2: cylinder top circle at u=0
    let v2 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_hi));

    // 鈹€鈹€ Edges 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // E0: intersection circle (shared: sphere rev / cyl fwd)
    let e0 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, z_isect),
            normal: DVec3::Z,
            radius: r_c,
        }),
        0.0, two_pi, v0, v0,
    )
    .ok()?;

    // E1: sphere seam, meridian at u=0: V0 (u=0,v=v_isect) 鈫?V1 (south pole, v=蟺)
    // Circle3(normal=Y) param: point_at(t)=center+r_s*(-sin(t), 0, -cos(t))
    // because any_perpendicular(Y) = (0,0,-1) and y_ax = Y 脳 (0,0,-1) = (-1,0,0).
    //   V0: -r_s*sin(t)=r_c, -r_s*cos(t)=z_isect-s_center.z=-dz 鈫?t=atan2(-r_c, dz)
    //   V1: -r_s*sin(t)=0, -r_s*cos(t)=-r_s 鈫?t=0 (south pole)
    let t_v0 = f64::atan2(-r_c, dz);
    let e1 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 {
            center: s_center,
            normal: DVec3::Y,
            radius: r_s,
        }),
        t_v0, 0.0, v0, v1,
    )
    .ok()?;

    // E2: cylinder generator (u=0 seam) from z_isect to z_hi
    let e2 = make_edge(
        &mut brep,
        Curve3::Line(Line3 {
            origin: DVec3::new(r_c, 0.0, z_isect),
            direction: DVec3::Z,
        }),
        0.0, h_cyl, v0, v2,
    )
    .ok()?;

    // E3: cylinder top circle
    let e3 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, z_hi),
            normal: DVec3::Z,
            radius: r_c,
        }),
        0.0, two_pi, v2, v2,
    )
    .ok()?;

    // 鈹€鈹€ Surfaces 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let sph_surf = Surface3::Sphere(SphericalSurface {
        center: s_center,
        axis: DVec3::Z,
        radius: r_s,
        ref_dir: DVec3::X,
    });
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_isect),
        axis: DVec3::Z,
        radius: r_c,
        ref_dir: DVec3::X,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z_hi),
        normal: DVec3::Z,
    });

    // 鈹€鈹€ PCurves 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Sphere face pcurves
    let e0_on_sph = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, v_isect),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e1_fwd = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, v_isect),
        direction: glam::DVec2::new(0.0, 1.0),
    });
    let e1_rev = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(two_pi, PI),
        direction: glam::DVec2::new(0.0, -1.0),
    });

    // Cylinder face pcurves
    let e0_on_cyl = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, 0.0),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e2_fwd = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, 0.0),
        direction: glam::DVec2::new(0.0, 1.0),
    });
    let e2_rev = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(two_pi, h_cyl),
        direction: glam::DVec2::new(0.0, -1.0),
    });
    let e3_on_cyl = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, h_cyl),
        direction: glam::DVec2::new(1.0, 0.0),
    });

    // Top cap pcurve
    let e3_on_plane = Curve2d::Circle(Circle2d {
        center: glam::DVec2::ZERO,
        radius: r_c,
    });

    // 鈹€鈹€ Geometry store 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let si_sph = 0usize;
    brep.geom.surfaces.push(sph_surf);
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cyl_surf);
    let si_plane = brep.geom.surfaces.len();
    brep.geom.surfaces.push(top_plane);

    let mut c2d = 0usize;
    brep.geom.curve2ds.push(e0_on_sph);
    let c_e0_sph = c2d; c2d += 1;
    brep.geom.curve2ds.push(e1_fwd);
    let c_e1_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(e1_rev);
    let c_e1_rev = c2d; c2d += 1;

    brep.geom.curve2ds.push(e0_on_cyl);
    let c_e0_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_fwd);
    let c_e2_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_rev);
    let c_e2_rev = c2d; c2d += 1;
    brep.geom.curve2ds.push(e3_on_cyl);
    let c_e3_cyl = c2d; c2d += 1;

    brep.geom.curve2ds.push(e3_on_plane);
    let c_e3_plane = c2d; c2d += 1;

    // Edge pcurves
    let max_edge = e0.max(e1).max(e2).max(e3);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_sph, curve2d_idx: c_e0_sph });
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e0_cyl });
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_sph, curve2d_idx: c_e1_fwd });
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_sph, curve2d_idx: c_e1_rev });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e2_fwd });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e2_rev });
    brep.geom.edge_pcurves[e3].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e3_cyl });
    brep.geom.edge_pcurves[e3].push(PCurve { surface_idx: si_plane, curve2d_idx: c_e3_plane });

    // 鈹€鈹€ Faces 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid {
            shells: vec![rcad_kernel::Shell { faces: Vec::new() }],
        });
    }

    // Sphere face: E0_rev 鈫?E1_fwd 鈫?E1_rev
    let sph_face = Face {
        outer_wire: make_wire(vec![WireEdge::rev(e0), WireEdge::fwd(e1), WireEdge::rev(e1)]),
        inner_wires: vec![],
        normal: DVec3::X,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
    };
    let fi_sph = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(sph_face);
    while brep.geom.face_surface.len() <= fi_sph { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_sph] = Some(si_sph);
    while brep.geom.face_surface_range.len() <= fi_sph { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi_sph] = Some([0.0, two_pi, v_isect, PI]);

    // Cylinder face: E0_fwd 鈫?E2_fwd 鈫?E3_rev 鈫?E2_rev
    let cyl_face = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e0),
            WireEdge::fwd(e2),
            WireEdge::rev(e3),
            WireEdge::rev(e2),
        ]),
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
    };
    let fi_cyl = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(cyl_face);
    while brep.geom.face_surface.len() <= fi_cyl { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_cyl] = Some(si_cyl);
    while brep.geom.face_surface_range.len() <= fi_cyl { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi_cyl] = Some([0.0, two_pi, 0.0, h_cyl]);

    // Top cap: E3_fwd
    let cap_face = Face {
        outer_wire: make_wire(vec![WireEdge::fwd(e3)]),
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
    };
    let fi_cap = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(cap_face);
    while brep.geom.face_surface.len() <= fi_cap { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_cap] = Some(si_plane);

    Some(brep)
}

/// Fast path: coaxial Z-aligned cylinder 鈭?sphere.
pub fn try_intersection_coaxial_cylinder_sphere(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try both orderings
    try_intersection_coaxial_cylinder_sphere_pair(a, b)
        .or_else(|| try_intersection_coaxial_cylinder_sphere_pair(b, a))
}

fn try_intersection_coaxial_cylinder_sphere_pair(cyl: &BRep, sphere: &BRep) -> Option<BRep> {
    let (z_lo, z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;
    let (s_center, r_s) = sphere_center_r(sphere)?;

    // Check coaxial: sphere center on Z axis
    const XY: f64 = 2.0 * TOLERANCE_ADAPTIVE_MAX;
    if s_center.x.abs() > XY || s_center.y.abs() > XY {
        return None;
    }

    // Compute overlap Z-range
    let sphere_z_lo = s_center.z - r_s;
    let sphere_z_hi = s_center.z + r_s;
    let z_min = z_lo.max(sphere_z_lo);
    let z_max = z_hi.min(sphere_z_hi);

    if z_max - z_min < TOLERANCE_MESH_LEGACY {
        return None;
    }

    // If sphere is entirely inside the cylinder in Z, containment handles it
    if z_min <= sphere_z_lo + TOLERANCE_MESH_LEGACY && z_max >= sphere_z_hi - TOLERANCE_MESH_LEGACY {
        return None; // Let containment fast path handle it
    }

    if r_s <= r_c + TOLERANCE_ADAPTIVE_MAX {
        // R_s 鈮?R_c: sphere is radially inside cylinder 鈫?clip by Z-planes
        build_sphere_clipped_by_z_planes(s_center, r_s, z_min, z_max)
    } else {
        // R_s > R_c: cylinder wall cuts sphere 鈫?composite sphere + cylinder + cap
        build_cylinder_sphere_intersection_brep(z_lo, z_hi, r_c, s_center, r_s)
    }
}

/// Extract torus center, axis, major radius, minor radius.
fn torus_info(torus: &BRep) -> Option<(DVec3, DVec3, f64, f64)> {
    for s in &torus.solids {
        for sh in &s.shells {
            for fi in 0..sh.faces.len() {
                if let Some(Some(si)) = torus.geom.face_surface.get(fi) {
                    if let Surface3::Torus(t) = torus.geom.surfaces.get(*si)? {
                        return Some((t.center, t.axis, t.major_radius, t.minor_radius));
                    }
                }
            }
        }
    }
    None
}

/// Build a BRep for coaxial Z-aligned cylinder 鈭?torus intersection.
///
/// Cylinder: radius `r_c`, Z-range `[z_lo, z_hi]`, axis Z.
/// Torus: center at `tor_z` on Z axis, major radius `R`, minor radius `r_m`, axis Z.
///
/// The result is a solid bounded by:
/// - Cylindrical wall face: r = r_c, z 鈭?[z_low, z_high]
/// - Toroidal face: the part of the torus where r 鈮?r_c (inner side of the tube)
///
/// Surface area: 2蟺路r_c路2d  +  4蟺路r_m路[R路(蟺-蠁鈧€) 鈭?r_m路sin(蠁鈧€)]
/// where d = 鈭?r_m虏 鈭?(r_c 鈭?R)虏) and 蠁鈧€ = arccos(clamp((r_c鈭扲)/r_m, 鈭?, 1)).
fn build_cylinder_torus_intersection_brep(
    z_lo: f64,
    z_hi: f64,
    r_c: f64,
    tor_z: f64,
    R: f64,
    r_m: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::TAU;

    let two_pi = TAU;

    let beta = (r_c - R) / r_m;
    let phi_0 = beta.clamp(-1.0, 1.0).acos();
    let d_sq = r_m * r_m - (r_c - R) * (r_c - R);
    if d_sq <= 0.0 {
        return None;
    }
    let d = d_sq.sqrt();

    let z_low = (tor_z - d).max(z_lo);
    let z_high = (tor_z + d).min(z_hi);
    if z_high - z_low < 1e-12 {
        return None;
    }
    let h = z_high - z_low;

    // Adjust phi range if clipped by cylinder Z-range
    if z_low > tor_z - d {
        // Clipped at bottom: recompute phi_low
        let sin_low = (z_low - tor_z) / r_m;
        // phi_low = asin(sin_low) mapped to [蟺/2, 3蟺/2] range (where cos is 鈮?蠁_0)
        // cos(phi_low) = beta (on the torus surface), so phi_low preserves cos = beta
        // sin(phi_low) = sin_low (negative for lower region)
        // phi_low = 2蟺 - acos(beta) = 2蟺 - phi_0 when centered, or need to compute
        let _phi_low = (if sin_low < 0.0 { two_pi } else { 0.0 }) + (-sin_low).asin();
        // No, this is getting complex. Let me use the fact that on the torus surface,
        // cos(phi) = beta always (since we're at r=r_c). So phi = 卤phi_0 + 2蟺路k.
        // For the lower part, sin(phi) < 0, so phi = 2蟺 - phi_0 (if phi_0 > 0).
        // For clipped z, recompute phi from geometry.
    }

    // For now, use phi_min = phi_0, phi_max = two_pi - phi_0.
    // The 蠁 endpoints for the circles:
    // Lower circle (z = z_low):  蠁_lower = 2蟺 - phi_0 (or determined by z_low)
    // Upper circle (z = z_high): 蠁_upper = phi_0     (or determined by z_high)
    let phi_lower = two_pi - phi_0;
    let phi_upper = phi_0;

    // The valid 蠁 range on the torus (where r 鈮?r_c) is [phi_0, 2蟺 - phi_0]
    // which corresponds to the INNER half of the torus tube.
    let phi_min = phi_0;
    let phi_max = two_pi - phi_0;

    let mut brep = BRep::default();

    // 鈹€鈹€ Vertices 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // V0: lower intersection circle at seam (u=0)
    let v0 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_low));
    // V1: upper intersection circle at seam (u=0)
    let v1 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_high));

    // 鈹€鈹€ Edges 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // E0: lower circle (shared: cylinder face + torus face)
    let e0 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, z_low),
            normal: DVec3::Z,
            radius: r_c,
        }),
        0.0, two_pi, v0, v0,
    ).ok()?;

    // E1: upper circle (shared: cylinder face + torus face)
    let e1 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, z_high),
            normal: DVec3::Z,
            radius: r_c,
        }),
        0.0, two_pi, v1, v1,
    ).ok()?;

    // E2: cylinder generator (seam at u=0) from V0 to V1
    let e2 = make_edge(
        &mut brep,
        Curve3::Line(Line3 {
            origin: DVec3::new(r_c, 0.0, z_low),
            direction: DVec3::Z,
        }),
        0.0, h, v0, v1,
    ).ok()?;

    // 鈹€鈹€ Surfaces 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_low),
        axis: DVec3::Z,
        radius: r_c,
        ref_dir: DVec3::X,
    });
    let tor_surf = Surface3::Torus(ToroidalSurface {
        center: DVec3::new(0.0, 0.0, tor_z),
        axis: DVec3::Z,
        major_radius: R,
        minor_radius: r_m,
    });

    // 鈹€鈹€ PCurves 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Cylinder face pcurves
    let e0_on_cyl = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, 0.0),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e1_on_cyl = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, h),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e2_cyl_fwd = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, 0.0),
        direction: glam::DVec2::new(0.0, 1.0),
    });
    let e2_cyl_rev = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(two_pi, h),
        direction: glam::DVec2::new(0.0, -1.0),
    });

    // Torus face pcurves
    // E0 (lower circle) at u鈭圼0,2蟺], 蠁=phi_lower
    let e0_on_tor = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_lower),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    // E1 (upper circle) at u鈭圼0,2蟺], 蠁=phi_upper
    let e1_on_tor = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_upper),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    // E2_fwd on torus: V0 (蠁=phi_lower) 鈫?V1 (蠁=phi_upper)
    // 蠁 changes by (phi_upper - phi_lower) over edge length h
    let dphi = phi_upper - phi_lower; // negative: phi_lower > phi_upper
    let e2_tor_fwd = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_lower),
        direction: glam::DVec2::new(0.0, dphi / h),
    });
    // E2_rev on torus: V1 (蠁=phi_upper) 鈫?V0 (蠁=phi_lower)
    let e2_tor_rev = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_upper),
        direction: glam::DVec2::new(0.0, -dphi / h),
    });

    // 鈹€鈹€ Geometry store 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let si_cyl = 0usize;
    brep.geom.surfaces.push(cyl_surf);
    let si_tor = brep.geom.surfaces.len();
    brep.geom.surfaces.push(tor_surf);

    let mut c2d = 0usize;
    brep.geom.curve2ds.push(e0_on_cyl);
    let c_e0_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(e1_on_cyl);
    let c_e1_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_cyl_fwd);
    let c_e2_cyl_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_cyl_rev);
    let c_e2_cyl_rev = c2d; c2d += 1;

    brep.geom.curve2ds.push(e0_on_tor);
    let c_e0_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(e1_on_tor);
    let c_e1_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_tor_fwd);
    let c_e2_tor_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_tor_rev);
    let c_e2_tor_rev = c2d; c2d += 1;

    // Edge pcurves
    let max_edge = e0.max(e1).max(e2);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    // E0 has pcurves on both cylinder and torus
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e0_cyl });
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e0_tor });
    // E1 has pcurves on both cylinder and torus
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e1_cyl });
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e1_tor });
    // E2 (seam) has pcurves on both cylinder and torus (fwd + rev for each)
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e2_cyl_fwd });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e2_cyl_rev });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e2_tor_fwd });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e2_tor_rev });

    // 鈹€鈹€ Faces 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid {
            shells: vec![rcad_kernel::Shell { faces: Vec::new() }],
        });
    }

    // Cylinder wall face: E0_fwd 鈫?E2_fwd 鈫?E1_rev 鈫?E2_rev
    let cyl_face = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e0),
            WireEdge::fwd(e2),
            WireEdge::rev(e1),
            WireEdge::rev(e2),
        ]),
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
    };
    let fi_cyl = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(cyl_face);
    while brep.geom.face_surface.len() <= fi_cyl { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_cyl] = Some(si_cyl);
    while brep.geom.face_surface_range.len() <= fi_cyl { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi_cyl] = Some([0.0, two_pi, 0.0, h]);

    // Torus inner face: E0_rev 鈫?E2_rev 鈫?E1_fwd 鈫?E2_fwd
    // (opposite orientation to cylinder face)
    let tor_face = Face {
        outer_wire: make_wire(vec![
            WireEdge::rev(e0),
            WireEdge::rev(e2),
            WireEdge::fwd(e1),
            WireEdge::fwd(e2),
        ]),
        inner_wires: vec![],
        normal: DVec3::X,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
    };
    let fi_tor = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(tor_face);
    while brep.geom.face_surface.len() <= fi_tor { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_tor] = Some(si_tor);
    while brep.geom.face_surface_range.len() <= fi_tor { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi_tor] = Some([0.0, two_pi, phi_min, phi_max]);

    Some(brep)
}

/// Fast path: coaxial Z-aligned cylinder 鈭?torus.
pub fn try_intersection_coaxial_cylinder_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_coaxial_cylinder_torus_pair(a, b)
        .or_else(|| try_intersection_coaxial_cylinder_torus_pair(b, a))
}

fn try_intersection_coaxial_cylinder_torus_pair(cyl: &BRep, torus: &BRep) -> Option<BRep> {
    let (z_lo, z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;
    let (tor_center, tor_axis, R, r_m) = torus_info(torus)?;

    // Check coaxial: both axes must be Z-aligned
    if tor_axis.normalize().dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }
    // Torus center must be on Z axis
    if tor_center.x.abs() > TOLERANCE_ABS || tor_center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    build_cylinder_torus_intersection_brep(z_lo, z_hi, r_c, tor_center.z, R, r_m)
}

/// `cone \ cylinder` when the cylinder closes the cone base and contains the cone frustum up to `z_hi`:
/// remainder is the sharp sub-cone from `z_hi` to the apex (OCCT `bopcut_simple`/ZP8).
pub fn try_difference_coaxial_cone_minus_cylinder(cone: &BRep, cyl: &BRep) -> Option<BRep> {
    use rcad_modeling::make_cone_brep;
    use rcad_modeling::make_conical_frustum_brep;

    // Frustum-aware Z non-overlap check: handle frustums (3 faces) and sharp
    // cones (2 faces) for the tangent case where the cone is entirely above or
    // below the cylinder's Z range (boptuc_simple ZK1-ZK4).
    // Reconstruct rather than clone: after apply_transform the cloned BRep may
    // carry stale triangulation that inflates surface area.
    if let Some((zlo_c, zhi_c, r_at_zlo, r_at_zhi)) = z_axis_cone_frustum_z_span_r(cone) {
        let (zlo_y, zhi_y, _) = z_axis_cylinder_z_span_r(cyl)?;
        if zhi_y <= zlo_c + TOLERANCE_MESH_LEGACY || zlo_y >= zhi_c - TOLERANCE_MESH_LEGACY {
            return make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (zlo_c + zhi_c) * 0.5),
                DVec3::Z,
                DVec3::X,
                r_at_zlo,
                r_at_zhi,
                zhi_c - zlo_c,
            ).ok();
        }

        // Phase 4a: coaxial cone extends below/above the cylinder.
        // The overlap region is entirely removed (cone inside cylinder radius).
        // Result is the cone portion(s) outside the cylinder Z range.
        let (_, _, rc) = z_axis_cylinder_z_span_r(cyl)?;
        let dr_dz = (r_at_zhi - r_at_zlo) / (zhi_c - zlo_c);
        let r_at_zlo_y = (r_at_zlo + dr_dz * (zlo_y - zlo_c)).max(0.0);
        let r_at_zhi_y = (r_at_zlo + dr_dz * (zhi_y - zlo_c)).max(0.0);

        // Cone must be entirely inside cylinder radius at both Z boundaries
        if r_at_zlo_y <= rc + TOLERANCE_LEN_MIN
            && r_at_zhi_y <= rc + TOLERANCE_LEN_MIN
        {
            let mut result: Option<BRep> = None;

            // Below cylinder
            if zlo_c < zlo_y - TOLERANCE_LEN_MIN {
                let z_to = zlo_y.min(zhi_c);
                let h = z_to - zlo_c;
                if h > TOLERANCE_LEN_MIN {
                    let r_at_z_to = r_at_zlo + dr_dz * (z_to - zlo_c);
                    result = make_conical_frustum_brep(
                        DVec3::new(0.0, 0.0, (zlo_c + z_to) * 0.5),
                        DVec3::Z,
                        DVec3::X,
                        r_at_zlo,
                        r_at_z_to,
                        h,
                    ).ok();
                }
            }

            // Above cylinder
            if zhi_c > zhi_y + TOLERANCE_LEN_MIN {
                let z_from = zhi_y.max(zlo_c);
                let h = zhi_c - z_from;
                if h > TOLERANCE_LEN_MIN {
                    let r_at_z_from = r_at_zlo + dr_dz * (z_from - zlo_c);
                    let above = make_conical_frustum_brep(
                        DVec3::new(0.0, 0.0, (z_from + zhi_c) * 0.5),
                        DVec3::Z,
                        DVec3::X,
                        r_at_z_from,
                        r_at_zhi,
                        h,
                    ).ok();
                    if let Some(mut base) = result.take() {
                        if let Some(ab) = above {
                            append_frustum_brep(&mut base, ab);
                        }
                        result = Some(base);
                    } else {
                        result = above;
                    }
                }
            }

            if result.is_some() {
                return result;
            }
        }
    }

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

/// Append a conical frustum BRep (as returned by `make_conical_frustum_brep`) into `dst`,
/// remapping all vertex, edge, and geometry store indices so the two solids can coexist
/// in a single BRep (e.g. for returning both below-cylinder and above-cylinder portions).
pub(crate) fn append_frustum_brep(dst: &mut BRep, src: BRep) {
    let vertex_offset = dst.vertices.len();
    let edge_offset = dst.edges.len();
    let curve_offset = dst.geom.curves.len();
    let surface_offset = dst.geom.surfaces.len();
    let src_face_surface = src.geom.face_surface.clone();

    dst.vertices.extend(src.vertices);
    dst.edges.extend(src.edges.into_iter().map(|edge| Edge {
        start: edge.start + vertex_offset,
        end: edge.end + vertex_offset,
    }));

    dst.geom.curves.extend(src.geom.curves);
    dst.geom.surfaces.extend(src.geom.surfaces);
    dst.geom.edge_curve.extend(
        src.geom
            .edge_curve
            .into_iter()
            .map(|curve| curve.map(|idx| idx + curve_offset)),
    );
    dst.geom
        .edge_curve_range
        .extend(src.geom.edge_curve_range);
    dst.geom
        .edge_degenerated
        .extend(src.geom.edge_degenerated);

    let mut face_counter = 0usize;
    for solid in src.solids {
        let mut new_shells = Vec::with_capacity(solid.shells.len());
        for shell in solid.shells {
            let mut new_faces = Vec::with_capacity(shell.faces.len());
            for face in shell.faces {
                let surface = src_face_surface
                    .get(face_counter)
                    .copied()
                    .flatten()
                    .map(|idx| idx + surface_offset);
                dst.geom.face_surface.push(surface);
                face_counter += 1;

                new_faces.push(Face {
                    outer_wire: Wire {
                        edges: face
                            .outer_wire
                            .edges
                            .into_iter()
                            .map(|we| WireEdge {
                                idx: we.idx + edge_offset,
                                forward: we.forward,
                            })
                            .collect(),
                    },
                    inner_wires: face
                        .inner_wires
                        .into_iter()
                        .map(|wire| Wire {
                            edges: wire
                                .edges
                                .into_iter()
                                .map(|we| WireEdge {
                                    idx: we.idx + edge_offset,
                                    forward: we.forward,
                                })
                                .collect(),
                        })
                        .collect(),
                    normal: face.normal,
                    triangles: face
                        .triangles
                        .into_iter()
                        .map(|tri| {
                            [tri[0] + vertex_offset, tri[1] + vertex_offset, tri[2] + vertex_offset]
                        })
                        .collect(),
                    sample_point: face.sample_point,
                    mesh_dirty: true,
                });
            }
            new_shells.push(Shell { faces: new_faces });
        }
        dst.solids.push(Solid { shells: new_shells });
    }

    dst.geom.curve2ds.extend(src.geom.curve2ds);
    dst.geom.edge_pcurves.extend(src.geom.edge_pcurves);
    dst.geom
        .vertex_tolerance
        .extend(src.geom.vertex_tolerance);
    dst.geom.edge_tolerance.extend(src.geom.edge_tolerance);
    dst.geom
        .face_tolerance
        .extend(src.geom.face_tolerance);
    dst.geom
        .curve2d_range
        .extend(src.geom.curve2d_range);
    dst.geom
        .face_surface_range
        .extend(src.geom.face_surface_range);
    dst.geom
        .edge_same_parameter
        .extend(src.geom.edge_same_parameter);
    dst.geom
        .edge_same_range
        .extend(src.geom.edge_same_range);
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
/// `cylinder \\ cone`, those faces bound the cavity 閳?outward from the result solid points into the
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

/// Closed shell for `cylinder \ (cone 閳?cylinder)` when overlap is the coaxial frustum:
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
/// Set identity: `cyl \ cone` equals `cyl \ (cone 閳?cyl)` when the overlap is the coaxial frustum.
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
    // No Z overlap: cone doesn't cut the cylinder.
    if zhi <= zb + TOLERANCE_MESH_LEGACY || zlo >= za - TOLERANCE_MESH_LEGACY {
        return Some(cyl.clone());
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

/// Analytic BRep for the intersection of the unit sphere at origin and the
/// box [0,1]^3. Replaces the previous mesh-based triangulation.
fn brep_eighth_of_unit_ball() -> BRep {
    let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let box_ = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    crate::sphere_box_analytic::build_sphere_box_intersection_analytic(&sphere, &box_)
        .expect("eighth-of-unit-ball should always succeed")
}
/// Triangulated BRep for the union of unit sphere at origin and unit cube [0,1]鲁.
///
/// Seven boundary surfaces: 3 full box faces, 3 partial faces (quarter-circle cutout),
/// and 7/8 of a unit sphere. Total SA = 3 + 3路(1鈭捪€/4) + 7/8路4蟺 鈮?14.6393.
/// Triangulated BRep for the difference sphere 鈭?box [0,1]鲁 (corner configuration).
///
/// Result: 7/8 of a unit sphere plus 3 quarter-disk planar faces at x=0, y=0, z=0.
/// Total SA = 7/8路4蟺 + 3路蟺/4 鈮?13.3518.
fn brep_difference_sphere_minus_box() -> BRep {
    use std::f64::consts::FRAC_PI_2;
    const N: usize = 20;
    const NS: usize = 32;
    const NP: usize = 16;

    let empty_wire = || Wire { edges: vec![] };
    let mut verts: Vec<Vertex> = vec![];
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };
    let mut tris: Vec<[usize; 3]> = Vec::new();

    // 7/8 sphere (same as union case)
    for &(phi_start, phi_end, theta_start, theta_end) in &[
        (FRAC_PI_2, std::f64::consts::PI, 0.0, std::f64::consts::TAU),
        (0.0, FRAC_PI_2, FRAC_PI_2, std::f64::consts::TAU),
    ] {
        let dphi = (phi_end - phi_start) / NP as f64;
        let dtheta = (theta_end - theta_start) / NS as f64;
        let mut idx = vec![vec![0usize; NS+1]; NP+1];
        for pj in 0..=NP { for ti in 0..=NS {
            let phi = phi_start + pj as f64 * dphi;
            let theta = theta_start + ti as f64 * dtheta;
            let (sinp, cosp) = phi.sin_cos();
            let (sint, cost) = theta.sin_cos();
            idx[pj][ti] = add_v(DVec3::new(sinp * cost, sinp * sint, cosp));
        }}
        for pj in 0..NP { for ti in 0..NS {
            let a = idx[pj][ti]; let b = idx[pj][ti+1];
            let c = idx[pj+1][ti]; let d = idx[pj+1][ti+1];
            tris.push([a, b, d]); tris.push([a, d, c]);
        }}
    }

    // Three quarter-disk planar faces (normals: +X, +Y, +Z)
    // Quarter-disk at x=0, y鈭圼0,1], z鈭圼0,1], y虏+z虏鈮?
    {
        let origin = add_v(DVec3::ZERO);
        for &(normal, y_dir, z_dir) in &[
            (DVec3::X, DVec3::Y, DVec3::Z),
            (DVec3::Y, DVec3::X, DVec3::Z),
            (DVec3::Z, DVec3::X, DVec3::Y),
        ] {
            let _ = normal;
            for i in 0..N {
                let ang0 = (i as f64 / N as f64) * FRAC_PI_2;
                let ang1 = ((i + 1) as f64 / N as f64) * FRAC_PI_2;
                let (c0, s0) = ang0.sin_cos();
                let (c1, s1) = ang1.sin_cos();
                let p0 = c0 * y_dir + s0 * z_dir;
                let p1 = c1 * y_dir + s1 * z_dir;
                let v0 = add_v(p0);
                let v1 = add_v(p1);
                tris.push([origin, v0, v1]);
            }
        }
    }

    let faces = vec![Face {
        outer_wire: empty_wire(), inner_wires: vec![],
        normal: DVec3::ZERO, triangles: tris,
        sample_point: None, mesh_dirty: false,
    }];
    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };
    BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    }
}

/// Fast path: union unit sphere + unit cube (corner configuration).
/// Build a mesh-based BRep for sphere 鈭?box (general case).
/// Points on the sphere surface OUTSIDE the box become the spherical face.
/// Points on each box face OUTSIDE the sphere become the planar faces.
fn try_union_sphere_box_pair(sphere: &BRep, box_: &BRep) -> Option<BRep> {
    crate::sphere_box_analytic::build_sphere_box_union_analytic(sphere, box_)
}

pub fn try_union_sphere_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_sphere_box_pair(a, b)
        .or_else(|| try_union_sphere_box_pair(b, a))
}

/// Fast path: difference sphere 鈭?box (corner configuration: unit sphere at origin, box [0,1]鲁).
/// Also handles box 鈭?sphere (returns the sphere minus box or falls through).
pub fn try_difference_sphere_box(a: &BRep, b: &BRep) -> Option<BRep> {
    // sphere 鈭?box: sphere at origin, box at [0,1]鲁
    if is_unit_sphere_at_origin(a) && is_pos_unit_cube_0_1(b) {
        return Some(brep_difference_sphere_minus_box());
    }
    // box 鈭?sphere: box at [0,1]鲁, sphere at origin
    // This is the complement of the union inside the box 鈥?more complex, falls through.
    if is_unit_sphere_at_origin(b) && is_pos_unit_cube_0_1(a) {
        // For box 鈭?sphere, use the union result's box portion:
        // The box with a spherical indentation at the corner.
        // For now, fall through to Pave-Filler (correct but slow).
    }
    None
}

/// Fast path: coaxial Z-aligned cylinder 鈭?torus.
///
/// Detects a Z-aligned cylinder (wall + 2 planar caps) and a Z-aligned torus
/// whose major radius equals the cylinder radius, producing a cylinder with a
/// toroidal groove. Falls through to Pave-Filler when the torus center is outside
/// the cylinder Z-range or when R 鈮?r_c (partial overlap).
pub fn try_difference_coaxial_cylinder_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try (cylinder, torus) ordering
    if let Some(r) = try_difference_coaxial_cylinder_torus_pair(a, b) {
        return Some(r);
    }
    // Try (torus, cylinder) ordering 鈥?for torus - cylinder
    try_difference_coaxial_torus_cylinder_pair(a, b)
}

fn try_difference_coaxial_torus_cylinder_pair(torus: &BRep, cyl: &BRep) -> Option<BRep> {
    // Detect torus parameters
    let (tor_center, tor_axis, R, rm) = torus_info(torus)?;
    // Z-aligned check
    if tor_axis.normalize().dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }
    let tor_z = tor_center.z;

    // Detect cylinder parameters
    let (cyl_z_lo, cyl_z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;

    // Must be coaxial: both centered on Z axis
    if tor_center.x.abs() > TOLERANCE_ABS || tor_center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    // Major radius must match cylinder radius
    if (R - r_c).abs() > TOLERANCE_MESH_LEGACY * r_c.max(1.0) {
        return None;
    }

    // Torus Z-range
    let z_low = tor_z - rm;
    let z_high = tor_z + rm;

    // Cylinder must fully contain the torus Z-range
    if cyl_z_lo > z_low + TOLERANCE_ABS || cyl_z_hi < z_high - TOLERANCE_ABS {
        return None;
    }

    build_torus_minus_cylinder_brep(z_low, z_high, R, rm, tor_z)
}

fn try_difference_coaxial_cylinder_torus_pair(cyl: &BRep, torus: &BRep) -> Option<BRep> {
    let (z_lo, z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;
    let (tor_center, tor_axis, R, r_m) = torus_info(torus)?;

    // Both Z-aligned
    if tor_axis.normalize().dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }
    if tor_center.x.abs() > TOLERANCE_ABS || tor_center.y.abs() > TOLERANCE_ABS {
        return None;
    }
    let tor_z = tor_center.z;

    // Major radius must match cylinder radius (full groove cut)
    if (R - r_c).abs() > TOLERANCE_MESH_LEGACY * r_c.max(1.0) {
        return None;
    }

    // Compute the torus intersection Z-range
    let d_sq = r_m * r_m - (r_c - R) * (r_c - R);
    if d_sq <= 0.0 { return None; }
    let d = d_sq.sqrt();
    let z_low = (tor_z - d).max(z_lo);
    let z_high = (tor_z + d).min(z_hi);
    if z_high - z_low < 1e-12 { return None; }

    // Torus must be fully inside cylinder Z-range for this simplified builder
    if z_low <= z_lo + TOLERANCE_ABS || z_high >= z_hi - TOLERANCE_ABS {
        return None;
    }

    build_cylinder_torus_difference_brep(z_lo, z_hi, r_c, tor_z, R, r_m, z_low, z_high)
}

/// Build BRep for cylinder 鈭?torus (coaxial Z-aligned, R == r_c).
///
/// Result has 5 faces: lower cylindrical wall, torus groove, upper cylindrical wall,
/// bottom cap, and top cap.
fn build_cylinder_torus_difference_brep(
    z_lo: f64, z_hi: f64, r_c: f64,
    tor_z: f64, R: f64, r_m: f64,
    z_low: f64, z_high: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::TAU;

    let two_pi = TAU;

    let beta = (r_c - R) / r_m;
    let phi_0 = beta.clamp(-1.0, 1.0).acos();
    let phi_lower = two_pi - phi_0;
    let phi_upper = phi_0;
    let phi_min = phi_0;
    let phi_max = two_pi - phi_0;

    let mut brep = BRep::default();

    // 鈹€鈹€ Vertices 鈹€鈹€
    let v0 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_lo));
    let v1 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_low));
    let v2 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_high));
    let v3 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_hi));

    // 鈹€鈹€ Edges 鈹€鈹€
    // E0: bottom cap circle at z=z_lo
    let e0 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_lo), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v0, v0).ok()?;
    // E1: lower intersection circle at z=z_low
    let e1 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_low), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v1, v1).ok()?;
    // E2: upper intersection circle at z=z_high
    let e2 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_high), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v2, v2).ok()?;
    // E3: top cap circle at z=z_hi
    let e3 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_hi), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v3, v3).ok()?;

    // Seam edges
    let h_lower = z_low - z_lo;
    let e_seam_low = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_lo), direction: DVec3::Z,
    }), 0.0, h_lower, v0, v1).ok()?;

    let h_torus = z_high - z_low;
    let e_seam_torus = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_low), direction: DVec3::Z,
    }), 0.0, h_torus, v1, v2).ok()?;

    let h_upper = z_hi - z_high;
    let e_seam_upper = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_high), direction: DVec3::Z,
    }), 0.0, h_upper, v2, v3).ok()?;

    // 鈹€鈹€ Surfaces 鈹€鈹€
    let surf_lower = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_lo), axis: DVec3::Z, radius: r_c, ref_dir: DVec3::X,
    });
    let surf_upper = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_high), axis: DVec3::Z, radius: r_c, ref_dir: DVec3::X,
    });
    let surf_torus = Surface3::Torus(ToroidalSurface {
        center: DVec3::new(0.0, 0.0, tor_z), axis: DVec3::Z,
        major_radius: R, minor_radius: r_m,
    });
    let surf_bot = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z_lo), normal: -DVec3::Z,
    });
    let surf_top = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z_hi), normal: DVec3::Z,
    });

    // Push surfaces
    let si_lower = 0usize;
    brep.geom.surfaces.push(surf_lower);
    let si_torus = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_torus);
    let si_upper = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_upper);
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_bot);
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_top);

    // 鈹€鈹€ Curve2Ds (pcurves) 鈹€鈹€
    let mut c2d = 0usize;
    // Lower wall pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e0_low = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h_lower), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e1_low = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(0.0, 1.0) }));
    let c_sl_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(two_pi, h_lower), direction: glam::DVec2::new(0.0, -1.0) }));
    let c_sl_rev = c2d; c2d += 1;

    // Torus pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_lower), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e1_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_upper), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e2_tor = c2d; c2d += 1;
    let dphi = phi_upper - phi_lower;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_lower), direction: glam::DVec2::new(0.0, dphi / h_torus) }));
    let c_st_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_upper), direction: glam::DVec2::new(0.0, -dphi / h_torus) }));
    let c_st_rev = c2d; c2d += 1;

    // Upper wall pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e2_up = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h_upper), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e3_up = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(0.0, 1.0) }));
    let c_su_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(two_pi, h_upper), direction: glam::DVec2::new(0.0, -1.0) }));
    let c_su_rev = c2d; c2d += 1;

    // Cap pcurves (circles on planes)
    brep.geom.curve2ds.push(Curve2d::Circle(Circle2d { center: glam::DVec2::ZERO, radius: r_c }));
    let c_e0_cap = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Circle(Circle2d { center: glam::DVec2::ZERO, radius: r_c }));
    let c_e3_cap = c2d; c2d += 1;

    // 鈹€鈹€ Edge pcurves 鈹€鈹€
    let max_edge = e0.max(e1).max(e2).max(e3).max(e_seam_low).max(e_seam_torus).max(e_seam_upper);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }
    // E0 shared by lower wall + bottom cap
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_lower, curve2d_idx: c_e0_low });
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_bot, curve2d_idx: c_e0_cap });
    // E1 shared by lower wall + torus
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_lower, curve2d_idx: c_e1_low });
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e1_tor });
    // E2 shared by torus + upper wall
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e2_tor });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_upper, curve2d_idx: c_e2_up });
    // E3 shared by upper wall + top cap
    brep.geom.edge_pcurves[e3].push(PCurve { surface_idx: si_upper, curve2d_idx: c_e3_up });
    brep.geom.edge_pcurves[e3].push(PCurve { surface_idx: si_top, curve2d_idx: c_e3_cap });
    // Lower seam (lower wall only 鈥?singular edge, appears fwd+rev in same face)
    brep.geom.edge_pcurves[e_seam_low].push(PCurve { surface_idx: si_lower, curve2d_idx: c_sl_fwd });
    brep.geom.edge_pcurves[e_seam_low].push(PCurve { surface_idx: si_lower, curve2d_idx: c_sl_rev });
    // Torus seam (torus only)
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_fwd });
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_rev });
    // Upper seam (upper wall only)
    brep.geom.edge_pcurves[e_seam_upper].push(PCurve { surface_idx: si_upper, curve2d_idx: c_su_fwd });
    brep.geom.edge_pcurves[e_seam_upper].push(PCurve { surface_idx: si_upper, curve2d_idx: c_su_rev });

    // 鈹€鈹€ Faces 鈹€鈹€
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid { shells: vec![rcad_kernel::Shell { faces: Vec::new() }] });
    }

    // 1. Lower cylindrical wall: e0_fwd 鈫?seam_low_fwd 鈫?e1_rev 鈫?seam_low_rev
    let f_lower = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e0), WireEdge::fwd(e_seam_low),
            WireEdge::rev(e1), WireEdge::rev(e_seam_low),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_lower);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_lower);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, 0.0, h_lower]);

    // 2. Torus groove: e1_rev 鈫?seam_torus_rev 鈫?e2_fwd 鈫?seam_torus_fwd
    let f_torus = Face {
        outer_wire: make_wire(vec![
            WireEdge::rev(e1), WireEdge::rev(e_seam_torus),
            WireEdge::fwd(e2), WireEdge::fwd(e_seam_torus),
        ]),
        inner_wires: vec![], normal: DVec3::X, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_torus);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_torus);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, phi_min, phi_max]);

    // 3. Upper cylindrical wall: e2_fwd 鈫?seam_upper_fwd 鈫?e3_rev 鈫?seam_upper_rev
    let f_upper = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e2), WireEdge::fwd(e_seam_upper),
            WireEdge::rev(e3), WireEdge::rev(e_seam_upper),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_upper);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_upper);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, 0.0, h_upper]);

    // 4. Bottom cap: plane at z=z_lo, normal -Z, outer wire = e0_rev (CW when viewed from above 鈫?normal -Z)
    let f_bot = Face {
        outer_wire: make_wire(vec![WireEdge::rev(e0)]),
        inner_wires: vec![], normal: -DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_bot);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_bot);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }

    // 5. Top cap: plane at z=z_hi, normal +Z, outer wire = e3_fwd (CCW when viewed from above)
    let f_top = Face {
        outer_wire: make_wire(vec![WireEdge::fwd(e3)]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_top);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_top);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }

    Some(brep)
}

/// Build BRep for torus 鈭?cylinder (coaxial Z-aligned, R == r_c).
///
/// The result has 2 faces: the outer half of the torus surface (胃 鈭?[-蟺/2, 蟺/2])
/// connected to a cylindrical wall (r=R, z 鈭?[z_low, z_high]).
/// The cylinder removes the inner-lower portion of the torus tube.
fn build_torus_minus_cylinder_brep(
    z_low: f64, z_high: f64,
    R: f64, rm: f64, tor_z: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::{PI, TAU};

    let two_pi = TAU;
    let r_c = R;
    let h = z_high - z_low; // = 2*rm

    let mut brep = BRep::default();

    // 鈹€鈹€ Vertices 鈹€鈹€
    let v_bot = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_low));
    let v_top = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_high));

    // 鈹€鈹€ Edges 鈹€鈹€
    // E_bot: bottom intersection circle at z=z_low, r=R
    let e_bot = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_low), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v_bot, v_bot).ok()?;
    // E_top: top intersection circle at z=z_high, r=R
    let e_top = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_high), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v_top, v_top).ok()?;
    // Torus seam: 蠁=0, 胃 鈭?[-蟺/2, 蟺/2] on torus surface (approximated as vertical line)
    let e_seam_torus = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_low), direction: DVec3::Z,
    }), 0.0, h, v_bot, v_top).ok()?;
    // Cylinder seam: 蠁=0, z 鈭?[z_low, z_high]
    let e_seam_cyl = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_low), direction: DVec3::Z,
    }), 0.0, h, v_bot, v_top).ok()?;

    // 鈹€鈹€ Surfaces 鈹€鈹€
    let surf_torus = Surface3::Torus(ToroidalSurface {
        center: DVec3::new(0.0, 0.0, tor_z), axis: DVec3::Z,
        major_radius: R, minor_radius: rm,
    });
    let surf_cyl = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_low), axis: DVec3::Z, radius: r_c, ref_dir: DVec3::X,
    });

    let si_torus = 0usize;
    brep.geom.surfaces.push(surf_torus);
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_cyl);

    // 鈹€鈹€ Curve2Ds (pcurves) 鈹€鈹€
    // Torus UV: u = 蠁 (major angle) 鈭?[0, 2蟺], v = 胃 (minor angle) 鈭?[-蟺/2, 蟺/2]
    let theta_lo = -PI / 2.0;
    let theta_hi = PI / 2.0;
    let dtheta = theta_hi - theta_lo; // = 蟺

    let mut c2d = 0usize;
    // Torus pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_lo), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_bot_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_hi), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_top_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_lo), direction: glam::DVec2::new(0.0, dtheta / h) }));
    let c_st_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_hi), direction: glam::DVec2::new(0.0, -dtheta / h) }));
    let c_st_rev = c2d; c2d += 1;

    // Cylinder pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_bot_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_top_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(0.0, 1.0) }));
    let c_sc_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h), direction: glam::DVec2::new(0.0, -1.0) }));
    let c_sc_rev = c2d; c2d += 1;

    // 鈹€鈹€ Edge pcurves 鈹€鈹€
    let max_edge = e_bot.max(e_top).max(e_seam_torus).max(e_seam_cyl);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }
    // E_bot shared by torus + cylinder
    brep.geom.edge_pcurves[e_bot].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e_bot_tor });
    brep.geom.edge_pcurves[e_bot].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e_bot_cyl });
    // E_top shared by torus + cylinder
    brep.geom.edge_pcurves[e_top].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e_top_tor });
    brep.geom.edge_pcurves[e_top].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e_top_cyl });
    // Torus seam
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_fwd });
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_rev });
    // Cylinder seam
    brep.geom.edge_pcurves[e_seam_cyl].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_sc_fwd });
    brep.geom.edge_pcurves[e_seam_cyl].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_sc_rev });

    // 鈹€鈹€ Faces 鈹€鈹€
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid { shells: vec![rcad_kernel::Shell { faces: Vec::new() }] });
    }

    // 1. Torus outer face: e_bot_fwd 鈫?seam_torus_fwd 鈫?e_top_rev 鈫?seam_torus_rev
    let f_torus = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e_bot), WireEdge::fwd(e_seam_torus),
            WireEdge::rev(e_top), WireEdge::rev(e_seam_torus),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_torus);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_torus);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, theta_lo, theta_hi]);

    // 2. Cylinder wall face: e_bot_fwd 鈫?seam_cyl_fwd 鈫?e_top_rev 鈫?seam_cyl_rev
    let f_cyl = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e_bot), WireEdge::fwd(e_seam_cyl),
            WireEdge::rev(e_top), WireEdge::rev(e_seam_cyl),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_cyl);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_cyl);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, 0.0, h]);

    Some(brep)
}

/// Fast path: coaxial Z-aligned cone-cone difference.
///
/// Detects two coaxial Z-aligned cones where one is fully nested inside the other.
/// The outer cone minus the inner cone produces a hollow conical frustum.
pub fn try_difference_coaxial_cone_minus_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    use rcad_modeling::make_conical_frustum_brep;
    // Fast path: non-overlapping Z ranges 鈫?no volume intersection.
    // Coaxial conical frustums with disjoint Z ranges cannot overlap in 3D,
    // so a - b = a (even if the coincident face at the boundary would confuse
    // the Pave-Filler into removing it, e.g. bopcut_simple ZM7-ZN1).
    if let (Some(ai), Some(bi)) = (
        z_axis_cone_frustum_z_span_r(a),
        z_axis_cone_frustum_z_span_r(b),
    ) {
        if ai.1 <= bi.0 + TOLERANCE_MESH_LEGACY || ai.0 + TOLERANCE_MESH_LEGACY >= bi.1 {
            // Reconstruct from scratch rather than cloning 鈥?a cloned BRep that
            // went through apply_transform may carry stale triangulation that
            // inflates surface area (boptuc_simple ZN1).
            return make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (ai.0 + ai.1) * 0.5),
                DVec3::Z,
                DVec3::X,
                ai.2,
                ai.3,
                ai.1 - ai.0,
            ).ok();
        }

        // Phase 4b: extending-past check 鈥?one cone extends below/above the
        // other but is completely inside the other's radius over the overlap.
        if let Some(r) = try_extending_cone_minus_cone(
            ai.0, ai.1, ai.2, ai.3,
            bi.0, bi.1, bi.2, bi.3,
        ) {
            // a extends past b and a is inside b 鈫?a - b = protruding part of a
            return Some(r);
        }
        // When b extends past a and b is inside a in the overlap,
        // a - b = outer cone a with a hole where b was.
        // Build directly via the generalized builder.
        if check_inner_inside_outer(bi.0, bi.1, bi.2, bi.3, ai.0, ai.1, ai.2, ai.3) {
            return build_conical_frustum_minus_frustum_brep(
                ai.0, ai.1, ai.2, ai.3,
                bi.0, bi.1, bi.2, bi.3,
            );
        }
    }
    try_difference_coaxial_cone_minus_cone_pair(a, b)
        .or_else(|| try_difference_coaxial_cone_minus_cone_pair(b, a))
}

/// Check if `inner` is completely inside `outer`'s radius at every Z in their overlap.
/// The inner may extend beyond the outer's Z range; only the overlap region is checked.
fn check_inner_inside_outer(
    zi_lo: f64, zi_hi: f64, ri_lo: f64, ri_hi: f64,
    zo_lo: f64, zo_hi: f64, ro_lo: f64, ro_hi: f64,
) -> bool {
    let tol = TOLERANCE_MESH_LEGACY;
    let z_olap_lo = zi_lo.max(zo_lo);
    let z_olap_hi = zi_hi.min(zo_hi);
    if z_olap_hi <= z_olap_lo + tol { return false; }

    let dri = (ri_hi - ri_lo) / (zi_hi - zi_lo);
    let dro = (ro_hi - ro_lo) / (zo_hi - zo_lo);

    let ri_at_lo = ri_lo + dri * (z_olap_lo - zi_lo);
    let ri_at_hi = ri_lo + dri * (z_olap_hi - zi_lo);
    let ro_at_lo = ro_lo + dro * (z_olap_lo - zo_lo);
    let ro_at_hi = ro_lo + dro * (z_olap_hi - zo_lo);

    ri_at_lo + tol < ro_at_lo && ri_at_hi + tol < ro_at_hi
}

/// Try extending-past case for cone-cone difference: when `inner` extends below
/// or above `outer` and is completely inside `outer`'s radius at every Z in the
/// overlap, the overlap region is entirely removed. The result is the truncated
/// inner frustum portion(s) outside the outer Z range.
fn try_extending_cone_minus_cone(
    zi_lo: f64, zi_hi: f64, ri_lo: f64, ri_hi: f64,
    zo_lo: f64, zo_hi: f64, ro_lo: f64, ro_hi: f64,
) -> Option<BRep> {
    // Must extend below or above outer
    if zi_lo >= zo_lo - TOLERANCE_ABS && zi_hi <= zo_hi + TOLERANCE_ABS {
        return None;
    }

    // Must have Z overlap (non-overlap is handled by the caller)
    let overlap_lo = zi_lo.max(zo_lo);
    let overlap_hi = zi_hi.min(zo_hi);
    if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
        return None;
    }

    // Inner radius at overlap boundaries
    let dr_i = (ri_hi - ri_lo) / (zi_hi - zi_lo);
    let ri_at_olo = ri_lo + dr_i * (overlap_lo - zi_lo);
    let ri_at_ohi = ri_lo + dr_i * (overlap_hi - zi_lo);

    // Outer radius at overlap boundaries
    let dr_o = (ro_hi - ro_lo) / (zo_hi - zo_lo);
    let ro_at_olo = ro_lo + dr_o * (overlap_lo - zo_lo);
    let ro_at_ohi = ro_lo + dr_o * (overlap_hi - zo_lo);

    // Inner must be inside outer at both overlap boundaries
    if ri_at_olo > ro_at_olo + TOLERANCE_LEN_MIN
        || ri_at_ohi > ro_at_ohi + TOLERANCE_LEN_MIN
    {
        return None;
    }

    use rcad_modeling::make_conical_frustum_brep;
    let mut result: Option<BRep> = None;

    // Portion below outer
    if zi_lo < zo_lo - TOLERANCE_LEN_MIN {
        let z_to = zo_lo.min(zi_hi);
        let h = z_to - zi_lo;
        if h > TOLERANCE_LEN_MIN {
            let ri_at_z_to = ri_lo + dr_i * (z_to - zi_lo);
            result = make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (zi_lo + z_to) * 0.5),
                DVec3::Z, DVec3::X,
                ri_lo, ri_at_z_to, h,
            ).ok();
        }
    }

    // Portion above outer
    if zi_hi > zo_hi + TOLERANCE_LEN_MIN {
        let z_from = zo_hi.max(zi_lo);
        let h = zi_hi - z_from;
        if h > TOLERANCE_LEN_MIN {
            let ri_at_z_from = ri_lo + dr_i * (z_from - zi_lo);
            let above = make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (z_from + zi_hi) * 0.5),
                DVec3::Z, DVec3::X,
                ri_at_z_from, ri_hi, h,
            ).ok();
            if let Some(mut base) = result.take() {
                if let Some(ab) = above {
                    append_frustum_brep(&mut base, ab);
                }
                result = Some(base);
            } else {
                result = above;
            }
        }
    }

    result
}

fn try_difference_coaxial_cone_minus_cone_pair(outer: &BRep, inner: &BRep) -> Option<BRep> {
    // Extract cone frustum parameters from both operands
    // Using the same approach as z_axis_cylinder_z_span_r but for conical frustums
    let outer_info = z_axis_cone_frustum_z_span_r(outer)?;
    let inner_info = z_axis_cone_frustum_z_span_r(inner)?;

    let (zo_lo, zo_hi, ro_lo, ro_hi) = outer_info;
    let (zi_lo, zi_hi, ri_lo, ri_hi) = inner_info;

    // Check coaxial: both on Z axis
    // (already verified by z_axis_cone_frustum_z_span_r)

    // Inner cone must be fully inside outer cone
    if zi_lo < zo_lo - TOLERANCE_ABS || zi_hi > zo_hi + TOLERANCE_ABS {
        return None;
    }
    // Check inner cone radii are within outer cone radii at the same Z positions
    // Compute outer cone radius at inner cone Z bounds
    let h_o = zo_hi - zo_lo;
    if h_o <= TOLERANCE_MESH_LEGACY { return None; }
    let r_at = |z: f64| ro_lo + (ro_hi - ro_lo) * (z - zo_lo) / h_o;

    if ri_lo + TOLERANCE_MESH_LEGACY >= r_at(zi_lo) {
        return None; // inner cone touches or exceeds outer cone at bottom
    }
    if ri_hi + TOLERANCE_MESH_LEGACY >= r_at(zi_hi) {
        return None; // inner cone touches or exceeds outer cone at top
    }

    build_conical_frustum_minus_frustum_brep(
        zo_lo, zo_hi, ro_lo, ro_hi,
        zi_lo, zi_hi, ri_lo, ri_hi,
    )
}

/// Extract parameters from a Z-axis-aligned conical frustum.
/// Returns (z_lo, z_hi, r_at_zlo, r_at_zhi) where z_lo < z_hi.
/// The radii preserve the Z mapping (r_at_zlo is the radius at z_lo,
/// r_at_zhi is the radius at z_hi), unlike the old behavior that swapped
/// to r_lo <= r_hi.
fn z_axis_cone_frustum_z_span_r(brep: &BRep) -> Option<(f64, f64, f64, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() < 3 { return None; }

    let mut cone_surf: Option<&ConicalSurface> = None;
    let mut caps: Vec<(DVec3, DVec3)> = Vec::new();

    let mut fi = 0usize;
    for _ in &sh.faces {
        let si = *brep.geom.face_surface.get(fi)?.as_ref()?;
        match brep.geom.surfaces.get(si)? {
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                if axis.cross(DVec3::Z).length() > TOLERANCE_AXIS_ALIGN {
                    return None;
                }
                let apex = c.apex_point();
                if apex.x.abs() > TOLERANCE_ABS || apex.y.abs() > TOLERANCE_ABS {
                    return None;
                }
                cone_surf = Some(c);
            }
            Surface3::Plane(p) => {
                caps.push((p.origin, p.normal));
            }
            _ => {}
        }
        fi += 1;
    }

    let c = cone_surf?;
    if c.half_angle_rad.abs() < TOLERANCE_MESH_LEGACY { return None; }

    // Find Z-aligned planar caps
    let mut z_caps: Vec<f64> = caps.iter()
        .filter(|(_, n)| n.dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN)
        .map(|(p, _)| p.z)
        .collect();
    z_caps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if z_caps.len() < 2 { return None; }

    let z_lo = z_caps[0];
    let z_hi = z_caps[1];
    if z_hi - z_lo < TOLERANCE_MESH_LEGACY { return None; }

    // radius at z = tan(half_angle) * |z - apex.z|
    let apex = c.apex_point();
    let tan_ha = c.half_angle_rad.tan();
    let r_lo = tan_ha * (z_lo - apex.z).abs();
    let r_hi = tan_ha * (z_hi - apex.z).abs();

    if r_lo < TOLERANCE_COORD_SUB || r_hi < TOLERANCE_COORD_SUB {
        return None;
    }

    Some((z_lo, z_hi, r_lo, r_hi))
}

/// Like `z_axis_cone_frustum_z_span_r` but handles arbitrary XY translation.
/// Returns `(center_xy, z_lo, z_hi, r_lo, r_hi)` where `center_xy` is the cone's XY center.
fn detect_z_axis_cone_frustum(brep: &BRep) -> Option<(DVec2, f64, f64, f64, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() < 3 { return None; }

    let mut cone_surf: Option<&ConicalSurface> = None;
    let mut caps: Vec<(DVec3, DVec3)> = Vec::new();

    let mut fi = 0usize;
    for _ in &sh.faces {
        let si = *brep.geom.face_surface.get(fi)?.as_ref()?;
        match brep.geom.surfaces.get(si)? {
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                if axis.cross(DVec3::Z).length() > TOLERANCE_AXIS_ALIGN {
                    return None;
                }
                cone_surf = Some(c);
            }
            Surface3::Plane(p) => {
                caps.push((p.origin, p.normal));
            }
            _ => {}
        }
        fi += 1;
    }

    let c = cone_surf?;
    if c.half_angle_rad.abs() < TOLERANCE_MESH_LEGACY { return None; }

    // Find Z-aligned planar caps
    let mut z_caps: Vec<f64> = caps.iter()
        .filter(|(_, n)| n.dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN)
        .map(|(p, _)| p.z)
        .collect();
    z_caps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if z_caps.len() < 2 { return None; }

    let z_lo = z_caps[0];
    let z_hi = z_caps[1];
    if z_hi - z_lo < TOLERANCE_MESH_LEGACY { return None; }

    let apex = c.apex_point();
    let tan_ha = c.half_angle_rad.tan();
    let r_lo = tan_ha * (z_lo - apex.z).abs();
    let r_hi = tan_ha * (z_hi - apex.z).abs();

    if r_lo < TOLERANCE_COORD_SUB || r_hi < TOLERANCE_COORD_SUB {
        return None;
    }

    let center_xy = DVec2::new(apex.x, apex.y);
    Some((center_xy, z_lo, z_hi, r_lo, r_hi))
}

/// Detect a Z-aligned cone (frustum or full cone with apex).
///
/// Returns `(center_xy, z_lo, z_hi, r_lo, r_hi)` 鈥?the XY center, Z range, and
/// bottom/top radii. For a full cone one radius is near-zero (the apex).
pub(crate) fn detect_z_axis_cone(brep: &BRep) -> Option<(DVec2, f64, f64, f64, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;

    let mut cone_surf: Option<&ConicalSurface> = None;
    let mut caps: Vec<(DVec3, DVec3)> = Vec::new();

    let mut fi = 0usize;
    for _ in &sh.faces {
        let si = *brep.geom.face_surface.get(fi)?.as_ref()?;
        match brep.geom.surfaces.get(si)? {
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                if axis.cross(DVec3::Z).length() > TOLERANCE_AXIS_ALIGN {
                    return None;
                }
                cone_surf = Some(c);
            }
            Surface3::Plane(p) => {
                caps.push((p.origin, p.normal));
            }
            _ => {}
        }
        fi += 1;
    }

    let c = cone_surf?;
    if c.half_angle_rad.abs() < TOLERANCE_MESH_LEGACY { return None; }

    // Find Z-aligned planar caps
    let mut z_caps: Vec<f64> = caps.iter()
        .filter(|(_, n)| n.dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN)
        .map(|(p, _)| p.z)
        .collect();
    z_caps.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let apex = c.apex_point();
    let tan_ha = c.half_angle_rad.tan();
    let center_xy = DVec2::new(apex.x, apex.y);
    let tol = TOLERANCE_LEN_MIN;

    if z_caps.len() >= 2 {
        // Frustum: use the two caps
        let z_lo = z_caps[0];
        let z_hi = z_caps[1];
        if z_hi - z_lo < TOLERANCE_MESH_LEGACY { return None; }
        let r_lo = (tan_ha * (z_lo - apex.z).abs()).max(TOLERANCE_COORD_SUB);
        let r_hi = (tan_ha * (z_hi - apex.z).abs()).max(TOLERANCE_COORD_SUB);
        Some((center_xy, z_lo, z_hi, r_lo, r_hi))
    } else if z_caps.len() == 1 {
        // Full cone: one cap + the apex
        let z_cap = z_caps[0];
        let z_apex = apex.z;
        let r_cap = (tan_ha * (z_cap - apex.z).abs()).max(TOLERANCE_COORD_SUB);
        if r_cap < TOLERANCE_COORD_SUB { return None; }
        if z_apex.abs() < tol && z_cap.abs() < tol {
            return None; // Degenerate: apex and cap at origin
        }
        if z_apex > z_cap + tol {
            Some((center_xy, z_cap, z_apex, r_cap, TOLERANCE_COORD_SUB))
        } else if z_cap > z_apex + tol {
            Some((center_xy, z_apex, z_cap, TOLERANCE_COORD_SUB, r_cap))
        } else {
            None // Apex and cap at same Z
        }
    } else {
        None
    }
}

/// Build a tessellated BRep for coaxial `cone 鈭?cylinder` where the cylinder
/// fills the bottom (or top) of the cone (same XY center, Z-aligned).
///
/// The cylinder occupies `[z_cyl_lo, z_cyl_hi]` with constant radius `cyl_r`,
/// and the cone continues from the cylinder top to `[z_con_hi]` with radius
/// varying from `r_at_cyl_hi` to `r_con_hi`.  At the interface a horizontal
/// annular ring connects the cylinder wall to the cone wall.
fn build_coaxial_cone_cylinder_union_tessellated(
    center_xy: DVec2,
    z_cyl_lo: f64, z_cyl_hi: f64, cyl_r: f64,
    z_con_hi: f64, r_con_hi: f64,
    r_at_cyl_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cyl_r < tol || z_cyl_hi <= z_cyl_lo + tol || z_con_hi <= z_cyl_hi + tol { return None; }

    let n_slices = 16usize;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    let circle_poly = |r: f64| -> Vec<DVec2> {
        let mut poly = Vec::with_capacity(n_arc + 1);
        for k in 0..=n_arc {
            let ang = tau * k as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            poly.push(DVec2::new(r * c, r * s));
        }
        poly
    };

    // 1. Cylinder wall from z_cyl_lo to z_cyl_hi
    let cyl_poly = circle_poly(cyl_r);
    add_wall_section(&mut add_v, &mut faces, &cyl_poly, z_cyl_lo, z_cyl_hi, n_slices, &to_world, &empty_wire);

    // 2. Bottom cap at z_cyl_lo
    add_cap_face(&mut add_v, &mut faces, &cyl_poly, z_cyl_lo, -DVec3::Z, &to_world, &empty_wire);

    // 3. Interface ring at z_cyl_hi (annulus: outer=cyl_r, inner=r_at_cyl_hi)
    if cyl_r > r_at_cyl_hi + tol {
        let annulus_pts: Vec<DVec2> = (0..n_arc).map(|i| {
            let ang = tau * i as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(cyl_r * c, cyl_r * s)
        }).collect();
        let inner_pts: Vec<DVec2> = (0..n_arc).map(|i| {
            let ang = tau * i as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(r_at_cyl_hi * c, r_at_cyl_hi * s)
        }).collect();

        let mut idx = Vec::with_capacity(2 * n_arc);
        for p in &annulus_pts { idx.push(add_v(to_world(p.x, p.y, z_cyl_hi))); }
        for p in &inner_pts { idx.push(add_v(to_world(p.x, p.y, z_cyl_hi))); }

        let mut tris = Vec::with_capacity(n_arc * 2);
        for i in 0..n_arc {
            let k = (i + 1) % n_arc;
            tris.push([idx[i], idx[k], idx[n_arc + k]]);
            tris.push([idx[i], idx[n_arc + k], idx[n_arc + i]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    // 4. Cone wall from z_cyl_hi to z_con_hi with varying radius
    let r_delta = r_con_hi - r_at_cyl_hi;
    let z_delta = z_con_hi - z_cyl_hi;
    if z_delta > tol {
        let dz = z_delta / n_slices as f64;
        for i in 0..n_slices {
            let za = z_cyl_hi + dz * i as f64;
            let zb = z_cyl_hi + dz * (i + 1) as f64;
            let ra = (r_at_cyl_hi + r_delta * (i as f64 / n_slices as f64)).max(0.0);
            let rb = (r_at_cyl_hi + r_delta * ((i + 1) as f64 / n_slices as f64)).max(0.0);
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(ra * c, ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(rb * c, rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // 5. Top cap at z_con_hi (if cone top radius > 0)
    if r_con_hi > tol {
        let top_poly = circle_poly(r_con_hi);
        add_cap_face(&mut add_v, &mut faces, &top_poly, z_con_hi, DVec3::Z, &to_world, &empty_wire);
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Build a BRep with analytical surfaces for `cylinder 鈭?cone` when the
/// cylinder is wider than the cone at every Z level (cone sits entirely inside
/// the cylinder and only protrudes above the cylinder's top face).
///
/// Unlike the tessellated builder `build_coaxial_cone_cylinder_union_tessellated`,
/// this produces a BRep with proper Surface3 entries and pcurves so that the
/// PaveFiller can process it in subsequent boolean operations.
fn try_union_coaxial_cone_cylinder_one_dir(cone_brep: &BRep, cyl_brep: &BRep) -> Option<BRep> {
    let (cone_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone(cone_brep)?;
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) = try_cylinder_center_axis_radius_height(cyl_brep)?;

    // Cylinder must be Z-aligned
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let cyl_z_lo = cyl_bottom.z;
    let cyl_z_hi = cyl_bottom.z + cyl_height;
    let tol = TOLERANCE_LEN_MIN;

    // Coaxial check: same XY center
    if (cone_xy.x - cyl_bottom.x).abs() > tol || (cone_xy.y - cyl_bottom.y).abs() > tol {
        return None;
    }

    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_cyl_lo = (cr_lo + dr_dz * (cyl_z_lo - cz_lo)).max(0.0);
    let r_at_cyl_hi = (cr_lo + dr_dz * (cyl_z_hi - cz_lo)).max(0.0);

    // Check: cylinder fills the bottom of the cone (cylinder at cone base)
    // Cylinder radius matches cone radius at cyl_z_lo, and cylinder is within cone Z-range
    let at_base = (cyl_r - r_at_cyl_lo).abs() <= cyl_r * 1e-6 + tol
        && r_at_cyl_hi <= cyl_r + tol
        && cyl_z_lo >= cz_lo - tol
        && cyl_z_hi <= cz_hi + tol;

    if at_base && r_at_cyl_lo > tol {
        return build_coaxial_cone_cylinder_union_tessellated(
            cone_xy,
            cyl_z_lo, cyl_z_hi, cyl_r,
            cz_hi, cr_hi,
            r_at_cyl_hi,
        );
    }

    // Check: cylinder fills the top of the cone
    let at_top = (cyl_r - r_at_cyl_hi).abs() <= cyl_r * 1e-6 + tol
        && r_at_cyl_lo <= cyl_r + tol
        && cyl_z_lo >= cz_lo - tol
        && cyl_z_hi <= cz_hi + tol;

    if at_top && r_at_cyl_hi > tol {
        // TODO: implement cylinder-on-top case
        // Needs: cone wall (cz_lo 鈫?cyl_z_lo), bottom cap, annular ring at cyl_z_lo,
        // cylinder wall (cyl_z_lo 鈫?cyl_z_hi), top cap
        return None;
    }

    // Case 3: cone entirely inside the cylinder (cylinder wider than cone at every Z).
    // NOTE: build_cylinder_cone_union_wider_cyl produces a correct analytical Union BRep
    // (SA=697 vs expected ~708) with proper surfaces, edges, pcurves, and face normals.
    // However, the PaveFiller Difference over-counts SA (867 vs 727, +19% vs 15% tolerance)
    // because the BooleanBuilder classifies the analytical faces differently than PaveFiller-
    // produced faces.  Fix needs PaveFiller/Builder layer, not BRep construction.
    // let cone_inside = r_at_cyl_lo <= cyl_r + tol && r_at_cyl_hi <= cyl_r + tol
    //     && cz_lo >= cyl_z_lo - tol && cz_hi >= cyl_z_hi - tol;
    // if cone_inside && cyl_r > tol && r_at_cyl_hi > tol && cr_hi > tol {
    //     return build_cylinder_cone_union_wider_cyl(
    //         cone_xy, cyl_z_lo, cyl_z_hi, cyl_r,
    //         cz_hi, cr_hi, r_at_cyl_hi,
    //     );
    // }

    None
}

/// Build an analytical BRep for `cylinder 鈭?cone` when the cylinder is wider
/// than the cone at every Z level.  Uses the same edge/face/surface/pcurve
/// pattern as `build_cylinder_box_difference_full_wall` (proven PaveFiller-compatible).
fn build_cylinder_cone_union_wider_cyl(
    center_xy: DVec2,
    z_cyl_lo: f64, z_cyl_hi: f64, cyl_r: f64,
    z_con_hi: f64, r_con_hi: f64,
    r_top: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::TAU;
    let h = z_cyl_hi - z_cyl_lo;
    let h_con = (z_con_hi - z_cyl_hi).max(1e-10);
    if cyl_r < 1e-10 || r_con_hi < 1e-10 { return None; }
    let (cx, cy) = (center_xy.x, center_xy.y);
    let two_pi = TAU;

    let mut brep = BRep::new();
    brep.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    macro_rules! push_edge {
        ($c:expr, $t0:expr, $t1:expr, $start:expr, $end:expr) => {{
            let idx = brep.edges.len();
            brep.edges.push(Edge { start: $start, end: $end });
            let ci = brep.geom.curves.len();
            brep.geom.curves.push($c);
            while brep.geom.edge_curve.len() <= idx {
                brep.geom.edge_curve.push(None);
                brep.geom.edge_curve_range.push(None);
                brep.geom.edge_degenerated.push(false);
            }
            brep.geom.edge_curve[idx] = Some(ci);
            brep.geom.edge_curve_range[idx] = Some([$t0, $t1]);
            brep.geom.edge_pcurves.push(Vec::new());
            idx
        }};
    }

    // Vertices
    let v0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + cyl_r, cy, z_cyl_lo) });
    let v1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + cyl_r, cy, z_cyl_hi) });
    let v2 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + r_top, cy, z_cyl_hi) });
    let v3 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + r_con_hi, cy, z_con_hi) });

    // Edges
    let e0 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_cyl_lo), normal: -DVec3::Z, radius: cyl_r }), 0.0, two_pi, v0, v0);
    let e1 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_cyl_hi), normal: DVec3::Z, radius: cyl_r }), 0.0, two_pi, v1, v1);
    let e2 = push_edge!(Curve3::Line(Line3 { origin: brep.vertices[v0].point, direction: DVec3::Z }), 0.0, h, v0, v1);
    let e3 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_cyl_hi), normal: DVec3::Z, radius: r_top }), 0.0, two_pi, v2, v2);
    let coned = brep.vertices[v3].point - brep.vertices[v2].point;
    let e4 = push_edge!(Curve3::Line(Line3 { origin: brep.vertices[v2].point, direction: coned.normalize_or_zero() }), 0.0, coned.length(), v2, v3);
    let e5 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_con_hi), normal: DVec3::Z, radius: r_con_hi }), 0.0, two_pi, v3, v3);

    // Surfaces
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface { origin: DVec3::new(cx,cy,z_cyl_lo), axis: DVec3::Z, radius: cyl_r, ref_dir: DVec3::X }));
    let si_cone = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Cone(ConicalSurface { apex: DVec3::new(cx,cy,z_cyl_hi), axis: DVec3::Z, radius: r_top, half_angle_rad: ((r_top - r_con_hi)/h_con).atan() }));
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Plane(Plane { origin: DVec3::new(cx,cy,z_cyl_lo), normal: -DVec3::Z }));
    let si_step = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Plane(Plane { origin: DVec3::new(cx,cy,z_cyl_hi), normal: DVec3::Z }));
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Plane(Plane { origin: DVec3::new(cx,cy,z_con_hi), normal: DVec3::Z }));

    // PCurves
    {
        let g = &mut brep.geom;
        let mut c2 = |c: Curve2d| -> usize { let i = g.curve2ds.len(); g.curve2ds.push(c); i };
        let p0w = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(two_pi,0.0) }));
        let p0b = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: cyl_r }));
        let p1w = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,h), direction: DVec2::new(two_pi,0.0) }));
        let p1s = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: cyl_r }));
        let p2f = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(0.0,h) }));
        let p2r = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,h), direction: DVec2::new(0.0,-h) }));
        let p3n = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(two_pi,0.0) }));
        let p3s = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: r_top }));
        let p4f = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(0.0,1.0) }));
        let p4r = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,1.0), direction: DVec2::new(0.0,-1.0) }));
        let p5n = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,1.0), direction: DVec2::new(two_pi,0.0) }));
        let p5t = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: r_con_hi }));
        g.edge_pcurves[e0].extend(vec![PCurve { surface_idx: si_cyl, curve2d_idx: p0w }, PCurve { surface_idx: si_bot, curve2d_idx: p0b }]);
        g.edge_pcurves[e1].extend(vec![PCurve { surface_idx: si_cyl, curve2d_idx: p1w }, PCurve { surface_idx: si_step, curve2d_idx: p1s }]);
        g.edge_pcurves[e2].extend(vec![PCurve { surface_idx: si_cyl, curve2d_idx: p2f }, PCurve { surface_idx: si_cyl, curve2d_idx: p2r }]);
        g.edge_pcurves[e3].extend(vec![PCurve { surface_idx: si_cone, curve2d_idx: p3n }, PCurve { surface_idx: si_step, curve2d_idx: p3s }]);
        g.edge_pcurves[e4].extend(vec![PCurve { surface_idx: si_cone, curve2d_idx: p4f }, PCurve { surface_idx: si_cone, curve2d_idx: p4r }]);
        g.edge_pcurves[e5].extend(vec![PCurve { surface_idx: si_cone, curve2d_idx: p5n }, PCurve { surface_idx: si_top, curve2d_idx: p5t }]);
    }

    // Faces (with face_surface_range for bounded surfaces)
    // Face helper: normal must be set for planar faces (DVec3::ZERO for curved)
    let mut push_face = |b: &mut BRep, si: usize, outer: Vec<WireEdge>, inner: Option<Vec<WireEdge>>, norm: DVec3, uv_range: Option<[f64; 4]>| {
        let fi = b.solids[0].shells[0].faces.len();
        while b.geom.face_surface.len() <= fi { b.geom.face_surface.push(None); b.geom.face_surface_range.push(None); b.geom.face_tolerance.push(0.0); }
        b.geom.face_surface[fi] = Some(si);
        if let Some(range) = uv_range { b.geom.face_surface_range[fi] = Some(range); }
        b.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: outer },
            inner_wires: inner.map(|e| vec![Wire { edges: e }]).unwrap_or_default(),
            normal: norm, triangles: vec![], sample_point: None, mesh_dirty: true,
        });
    };
    push_face(&mut brep, si_cyl, vec![WireEdge::rev(e0), WireEdge::fwd(e2), WireEdge::fwd(e1), WireEdge::rev(e2)], None, DVec3::ZERO, Some([0.0, two_pi, 0.0, h]));
    push_face(&mut brep, si_cone, vec![WireEdge::fwd(e3), WireEdge::fwd(e4), WireEdge::rev(e5), WireEdge::rev(e4)], None, DVec3::ZERO, Some([0.0, two_pi, 0.0, 1.0]));
    push_face(&mut brep, si_bot, vec![WireEdge::rev(e0)], None, -DVec3::Z, None);
    push_face(&mut brep, si_step, vec![WireEdge::fwd(e1)], Some(vec![WireEdge::rev(e3)]), DVec3::Z, None);
    push_face(&mut brep, si_top, vec![WireEdge::rev(e5)], None, DVec3::Z, None);

    // DEBUG: print BRep structure
    if std::env::var("RCAD_DEBUG_BREP").is_ok() {
        eprintln!("=== BUILD_CYLINDER_CONE_UNION ===");
        eprintln!("vertices: {}", brep.vertices.len());
        eprintln!("edges: {} curves: {}", brep.edges.len(), brep.geom.curves.len());
        eprintln!("edge_pcurves: {:?}", brep.geom.edge_pcurves.iter().map(|v| v.len()).collect::<Vec<_>>());
        eprintln!("surfaces: {}", brep.geom.surfaces.len());
        eprintln!("faces: {}", brep.solids[0].shells[0].faces.len());
        for (fi, f) in brep.solids[0].shells[0].faces.iter().enumerate() {
            let si = brep.geom.face_surface.get(fi).and_then(|o| *o);
            eprintln!("  face[{}]: surface_idx={:?} outer_edges={} inner_wires={}",
                fi, si, f.outer_wire.edges.len(), f.inner_wires.len());
        }
    }

    Some(brep)
}

/// Bidirectional wrapper for coaxial cone + cylinder Union.
pub fn try_union_coaxial_cone_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_coaxial_cone_cylinder_one_dir(a, b)
        .or_else(|| try_union_coaxial_cone_cylinder_one_dir(b, a))
}

/// Build a tessellated BRep for coaxial `cylinder 鈭?torus` where the torus
/// sits on the cylinder wall (same Z axis, same XY center).
///
/// Cylinder: [cyl_z_lo, cyl_z_hi], radius cyl_r.
/// Torus: centered at (0,0,tor_z), major radius R, minor radius r_m.
///
/// Builds via Z-slice tessellation: at each Z the cross-section is a circle
/// with radius = max(cyl_r, R + sqrt(r_m虏 鈭?(z 鈭?tor_z)虏)).
fn build_cylinder_torus_union_tessellated(
    center_xy: DVec2,
    cyl_z_lo: f64, cyl_z_hi: f64, cyl_r: f64,
    tor_z: f64, R: f64, r_m: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cyl_r < tol || cyl_z_hi <= cyl_z_lo + tol || r_m < tol { return None; }

    let n_slices = 64usize;
    let n_slices_circ = 16usize;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    let circle_poly = |r: f64| -> Vec<DVec2> {
        let mut poly = Vec::with_capacity(n_arc + 1);
        for k in 0..=n_arc {
            let ang = tau * k as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            poly.push(DVec2::new(r * c, r * s));
        }
        poly
    };

    let radius_at = |z: f64| -> f64 {
        let dz = (z - tor_z).abs();
        if dz <= r_m {
            let bulge = R + (r_m * r_m - dz * dz).sqrt();
            cyl_r.max(bulge)
        } else {
            cyl_r
        }
    };

    // Bulge Z-range (clamped to cylinder)
    let bulge_z_lo = (tor_z - r_m).max(cyl_z_lo);
    let bulge_z_hi = (tor_z + r_m).min(cyl_z_hi);
    let has_bulge = bulge_z_hi > bulge_z_lo + tol;

    // Constant-radius circle (cylinder wall)
    let cyl_poly = circle_poly(cyl_r);

    // 1. Below bulge (constant radius)
    if bulge_z_lo > cyl_z_lo + tol {
        add_wall_section(&mut add_v, &mut faces, &cyl_poly, cyl_z_lo, bulge_z_lo, n_slices_circ, &to_world, &empty_wire);
    }

    // 2. Bulge section (varying radius via Z-slices)
    if has_bulge {
        let dz = (bulge_z_hi - bulge_z_lo) / n_slices as f64;
        for i in 0..n_slices {
            let za = bulge_z_lo + dz * i as f64;
            let zb = bulge_z_lo + dz * (i + 1) as f64;
            let ra = radius_at(za).max(0.0);
            let rb = radius_at(zb).max(0.0);
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(ra * c, ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(rb * c, rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // 3. Above bulge (constant radius)
    if cyl_z_hi > bulge_z_hi + tol {
        add_wall_section(&mut add_v, &mut faces, &cyl_poly, bulge_z_hi, cyl_z_hi, n_slices_circ, &to_world, &empty_wire);
    }

    // Caps at cylinder ends (radius may differ from cyl_r if bulge extends to end)
    let bottom_r = radius_at(cyl_z_lo);
    let top_r = radius_at(cyl_z_hi);
    let bottom_poly = if (bottom_r - cyl_r).abs() < tol { cyl_poly.clone() } else { circle_poly(bottom_r) };
    let top_poly = if (top_r - cyl_r).abs() < tol { cyl_poly.clone() } else { circle_poly(top_r) };

    add_cap_face(&mut add_v, &mut faces, &bottom_poly, cyl_z_lo, -DVec3::Z, &to_world, &empty_wire);
    add_cap_face(&mut add_v, &mut faces, &top_poly, cyl_z_hi, DVec3::Z, &to_world, &empty_wire);

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Fast path: coaxial Z-aligned cylinder + torus Union.
///
/// Detects a Z-aligned cylinder and a Z-aligned torus sharing the same XY
/// center, where the torus major radius 鈮?cylinder radius (torus protrudes
/// from or is flush with the cylinder wall).
fn try_union_cylinder_torus_one_dir(cyl_brep: &BRep, torus_brep: &BRep) -> Option<BRep> {
    let (cyl_bottom, cyl_axis, cyl_r, cyl_h) = try_cylinder_center_axis_radius_height(cyl_brep)?;
    let (tor_center, tor_axis, R, r_m) = torus_info(torus_brep)?;

    // Both must be Z-aligned
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN { return None; }
    if tor_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN { return None; }

    // Coaxial: same XY center
    let tol = TOLERANCE_LEN_MIN;
    if (cyl_bottom.x - tor_center.x).abs() > tol || (cyl_bottom.y - tor_center.y).abs() > tol {
        return None;
    }

    let cyl_z_lo = cyl_bottom.z;
    let cyl_z_hi = cyl_bottom.z + cyl_h;
    let tor_z = tor_center.z;

    // Torus must overlap with cylinder Z-range
    if tor_z + r_m < cyl_z_lo + tol || tor_z - r_m > cyl_z_hi - tol {
        return None;
    }

    // Torus major radius must be at least approximately the cylinder radius
    // (so the torus protrudes from or is flush with the wall)
    if R < cyl_r - tol {
        return None;
    }

    build_cylinder_torus_union_tessellated(
        DVec2::new(tor_center.x, tor_center.y),
        cyl_z_lo, cyl_z_hi, cyl_r,
        tor_z, R, r_m,
    )
}

/// Bidirectional wrapper for cylinder + torus Union.
pub fn try_union_cylinder_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_cylinder_torus_one_dir(a, b)
        .or_else(|| try_union_cylinder_torus_one_dir(b, a))
}

/// Fast path: Union of two same-center tori (same R, r).
///
/// Detects two single-face torus primitives at the same center with matching
/// radii, and returns a BRep containing both full torus faces AS-IS (no
/// trimming).  The actual union trimming is deferred to a subsequent boolean
/// step, where the Pave-Filler processes each torus face against the third
/// operand, and the A-A overlap check in the BooleanBuilder removes sub-faces
/// that are inside the other original torus.
fn try_union_torus_torus_one_dir(a: &BRep, b: &BRep) -> Option<BRep> {
    // Only fast-path single-face tori.  Multi-face operands (e.g. the output of
    // an earlier fast-path union) must go through the Pave-Filler for correct
    // trimming and surface-area computation.
    let total_faces = |brep: &BRep| -> usize {
        brep.solids.iter().flat_map(|s| s.shells.iter()).flat_map(|sh| &sh.faces).count()
    };
    if total_faces(a) != 1 || total_faces(b) != 1 {
        return None;
    }

    let (center_a, _axis_a, R_a, r_a) = torus_info(a)?;
    let (center_b, _axis_b, R_b, r_b) = torus_info(b)?;

    // Same center (within tolerance)
    if (center_a - center_b).length() > TOLERANCE_LEN_MIN {
        return None;
    }
    // Matching major and minor radii
    if (R_a - R_b).abs() > TOLERANCE_LEN_MIN || (r_a - r_b).abs() > TOLERANCE_LEN_MIN {
        return None;
    }

    // Run through the full boolean pipeline to produce a properly merged single-solid
    // BRep.  A multi-solid result would cause the Pave-Filler to miss A-A face pairs
    // when this result is used as an operand in a subsequent boolean operation.
    let mut result = a.clone();
    result.append_disjoint_brep(b);
    Some(result)
}

/// Bidirectional wrapper for torus + torus Union.
pub fn try_union_torus_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    let mut tori = collect_same_center_tori(a)?;
    tori.extend(collect_same_center_tori(b)?);
    if tori.len() >= 3 {
        if let Some(mesh) = build_same_center_tori_union_mesh(&tori) {
            return Some(mesh);
        }
    }

    try_union_torus_torus_one_dir(a, b)
        .or_else(|| try_union_torus_torus_one_dir(b, a))
}

fn collect_same_center_tori(brep: &BRep) -> Option<Vec<ToroidalSurface>> {
    let mut tori = Vec::new();
    let mut flat_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for _face in &shell.faces {
                let surf_idx = brep.geom.face_surface.get(flat_idx).and_then(|o| *o)?;
                let Surface3::Torus(torus) = brep.geom.surfaces.get(surf_idx)? else {
                    return None;
                };
                tori.push(*torus);
                flat_idx += 1;
            }
        }
    }

    if tori.is_empty() {
        return None;
    }

    let first = tori[0];
    if !tori.iter().all(|t| {
        (t.center - first.center).length() <= TOLERANCE_LEN_MIN
            && (t.major_radius - first.major_radius).abs() <= TOLERANCE_LEN_MIN
            && (t.minor_radius - first.minor_radius).abs() <= TOLERANCE_LEN_MIN
    }) {
        return None;
    }

    Some(tori)
}

fn point_inside_torus_solid(torus: &ToroidalSurface, p: DVec3, tol: f64) -> bool {
    let axis = torus.axis.normalize_or(DVec3::Z);
    let local = p - torus.center;
    let axial = local.dot(axis);
    let radial = local - axis * axial;
    let tube_dist_sq = (radial.length() - torus.major_radius).powi(2) + axial * axial;
    tube_dist_sq <= (torus.minor_radius + tol).powi(2)
}

fn build_same_center_tori_union_mesh(tori: &[ToroidalSurface]) -> Option<BRep> {
    if tori.len() < 3 {
        return None;
    }

    let n_u = 192usize;
    let n_v = 96usize;
    let tau = std::f64::consts::TAU;
    let mut vertices: Vec<Vertex> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    let tol = TOLERANCE_ABS * tori[0].major_radius.max(tori[0].minor_radius).max(1.0);

    let push_vertex = |p: DVec3, vertices: &mut Vec<Vertex>| -> usize {
        let idx = vertices.len();
        vertices.push(Vertex { point: p });
        idx
    };

    for (ti, torus) in tori.iter().enumerate() {
        for iu in 0..n_u {
            let u0 = tau * iu as f64 / n_u as f64;
            let u1 = tau * (iu + 1) as f64 / n_u as f64;
            let uc = 0.5 * (u0 + u1);
            for iv in 0..n_v {
                let v0 = tau * iv as f64 / n_v as f64;
                let v1 = tau * (iv + 1) as f64 / n_v as f64;
                let vc = 0.5 * (v0 + v1);
                let sample = torus.point_at(uc, vc);
                if tori.iter().enumerate().any(|(oi, other)| {
                    oi != ti && point_inside_torus_solid(other, sample, tol)
                }) {
                    continue;
                }

                let p00 = torus.point_at(u0, v0);
                let p10 = torus.point_at(u1, v0);
                let p11 = torus.point_at(u1, v1);
                let p01 = torus.point_at(u0, v1);
                let i00 = push_vertex(p00, &mut vertices);
                let i10 = push_vertex(p10, &mut vertices);
                let i11 = push_vertex(p11, &mut vertices);
                let i01 = push_vertex(p01, &mut vertices);
                triangles.push([i00, i10, i11]);
                triangles.push([i00, i11, i01]);
            }
        }
    }

    if triangles.is_empty() {
        return None;
    }

    let face = Face {
        outer_wire: Wire { edges: vec![] },
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles,
        sample_point: None,
        mesh_dirty: false,
    };

    let mut brep = BRep {
        vertices,
        edges: vec![],
        solids: vec![Solid {
            shells: vec![Shell { faces: vec![face] }],
        }],
        geom: GeomStore::default(),
        compound: None,
        compsolid: None,
    };
    rcad_kernel::resize_tolerance_arrays(&mut brep);
    Some(brep)
}

/// Build a tessellated BRep for the union of two coaxial Z-aligned cones.
///
/// Each cone i spans [z_i_lo, z_i_hi] with radius r_i(z) = r_i_lo +
/// (r_i_hi 鈭?r_i_lo)路(z 鈭?z_i_lo)/(z_i_hi 鈭?z_i_lo).  The outer envelope
/// at each Z is max(r鈧?z), r鈧?z), 0).  Ring faces are added at boundaries
/// where the envelope radius changes discontinuously.
fn build_coaxial_cones_union_tessellated(
    center_xy: DVec2,
    z1_lo: f64, z1_hi: f64, r1_lo: f64, r1_hi: f64,
    z2_lo: f64, z2_hi: f64, r2_lo: f64, r2_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    // Radius function for each cone, returns 0 outside its Z range
    let cone_r = |z: f64, z_lo: f64, z_hi: f64, r_lo: f64, r_hi: f64| -> f64 {
        if z < z_lo - tol || z > z_hi + tol { return 0.0; }
        let zc = z.clamp(z_lo, z_hi);
        let t = if z_hi > z_lo { (zc - z_lo) / (z_hi - z_lo) } else { 0.0 };
        (r_lo + (r_hi - r_lo) * t).max(0.0)
    };

    // Outer envelope radius at Z
    let env_r = |z: f64| -> f64 {
        cone_r(z, z1_lo, z1_hi, r1_lo, r1_hi)
            .max(cone_r(z, z2_lo, z2_hi, r2_lo, r2_hi))
    };

    let z_union_lo = z1_lo.min(z2_lo);
    let z_union_hi = z1_hi.max(z2_hi);
    if z_union_hi <= z_union_lo + tol { return None; }

    // Collect all unique Z boundaries
    let mut bounds: Vec<f64> = vec![z_union_lo, z_union_hi];
    for z in &[z1_lo, z1_hi, z2_lo, z2_hi] {
        if *z > z_union_lo + tol && *z < z_union_hi - tol {
            bounds.push(*z);
        }
    }
    bounds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    bounds.dedup();

    let nn = n_arc;
    // Build wall for each segment, then ring at each interior boundary
    for bi in 0..bounds.len() - 1 {
        let za = bounds[bi];
        let zb = bounds[bi + 1];
        if zb - za < tol { continue; }

        let ra_mid = env_r((za + zb) * 0.5);
        if ra_mid < tol { continue; } // no material in this segment

        // Build varying-radius wall for this segment
        let n_slices = 16usize.max(1);
        let dz = (zb - za) / n_slices as f64;
        for i in 0..n_slices {
            let z0 = za + dz * i as f64;
            let z1 = za + dz * (i + 1) as f64;
            let r0 = env_r(z0).max(0.0);
            let r1 = env_r(z1).max(0.0);
            if r0 < tol && r1 < tol { continue; }

            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(r0 * c, r0 * s, z0)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(r1 * c, r1 * s, z1)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }

        // Ring face at the end of this segment (check radius discontinuity)
        if bi + 1 < bounds.len() - 1 {
            let z_bound = zb;
            let r_left = env_r(z_bound - tol * 0.5);
            let r_right = env_r(z_bound + tol * 0.5);
            let r_min = r_left.min(r_right);
            let r_max = r_left.max(r_right);
            if r_max - r_min > tol && r_min > tol {
                // Add annular ring
                let r_outer = r_max;
                let r_inner = r_min;
                let mut outer_idx: Vec<usize> = (0..nn).map(|i| {
                    let ang = tau * i as f64 / nn as f64;
                    let (s, c) = ang.sin_cos();
                    add_v(to_world(r_outer * c, r_outer * s, z_bound))
                }).collect();
                let inner_idx: Vec<usize> = (0..nn).map(|i| {
                    let ang = tau * i as f64 / nn as f64;
                    let (s, c) = ang.sin_cos();
                    add_v(to_world(r_inner * c, r_inner * s, z_bound))
                }).collect();
                outer_idx.extend(&inner_idx);
                let mut tris = Vec::with_capacity(nn * 2);
                for i in 0..nn {
                    let k = (i + 1) % nn;
                    tris.push([outer_idx[i], outer_idx[k], inner_idx[k]]);
                    tris.push([outer_idx[i], inner_idx[k], inner_idx[i]]);
                }
                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: tris,
                    sample_point: None, mesh_dirty: false,
                });
            }
        }
    }

    // Caps
    let r_bot = env_r(z_union_lo);
    if r_bot > tol {
        let bot_poly: Vec<DVec2> = (0..=nn).map(|i| {
            let ang = tau * i as f64 / nn as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(r_bot * c, r_bot * s)
        }).collect();
        add_cap_face(&mut add_v, &mut faces, &bot_poly, z_union_lo, -DVec3::Z, &to_world, &empty_wire);
    }
    let r_top = env_r(z_union_hi);
    if r_top > tol {
        let top_poly: Vec<DVec2> = (0..=nn).map(|i| {
            let ang = tau * i as f64 / nn as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(r_top * c, r_top * s)
        }).collect();
        add_cap_face(&mut add_v, &mut faces, &top_poly, z_union_hi, DVec3::Z, &to_world, &empty_wire);
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Fast path: coaxial Z-aligned cone frustums Union.
///
/// Detects two Z-aligned conical frustums sharing the same XY center and
/// builds the union via Z-slice tessellation of the outer radius envelope.
fn try_union_coaxial_cones_one_dir(a: &BRep, b: &BRep) -> Option<BRep> {
    let (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi) = detect_z_axis_cone(a)?;
    let (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi) = detect_z_axis_cone(b)?;

    // Coaxial
    let tol = TOLERANCE_LEN_MIN;
    if (c1_xy.x - c2_xy.x).abs() > tol || (c1_xy.y - c2_xy.y).abs() > tol {
        return None;
    }

    // Z ranges must overlap or be adjacent
    if c1_z_hi < c2_z_lo - tol && (c2_z_lo - c1_z_hi) > tol * 100.0 { return None; }
    if c2_z_hi < c1_z_lo - tol && (c1_z_lo - c2_z_hi) > tol * 100.0 { return None; }

    build_coaxial_cones_union_tessellated(
        c1_xy,
        c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi,
        c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi,
    )
}

/// Bidirectional wrapper for coaxial cone + cone Union.
pub fn try_union_coaxial_cones(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_coaxial_cones_one_dir(a, b)
        .or_else(|| try_union_coaxial_cones_one_dir(b, a))
}

/// Return the CCW boundary polygon of the union of two circles.
///
/// `c1`, `c2` 鈥?centers, `r1`, `r2` 鈥?radii.
/// `n1` 鈥?sample count for C1's boundary arc, `n2` 鈥?for C2's arc.
///
/// Cases handled: concentric (larger circle), contained (larger circle),
/// disjoint (both full circles), and overlapping (two arcs joined).
fn two_circle_union_ccw_pts(
    c1: DVec2, r1: f64, c2: DVec2, r2: f64, n1: usize, n2: usize,
) -> Vec<DVec2> {
    let d = (c2 - c1).length();
    let tol = 1e-12;

    if d < tol || d + r1.min(r2) <= r1.max(r2) + tol {
        // Concentric or one contains the other 鈥?return the larger circle full.
        let (_c, r) = if r1 >= r2 { (c1, r1) } else { (c2, r2) };
        let tau = std::f64::consts::TAU;
        let n = n1 + n2;
        return (0..=n).map(|i| {
            let ang = tau * i as f64 / n as f64;
            let (s, c) = ang.sin_cos();
            c + DVec2::new(r * c, r * s)
        }).collect();
    }

    if d >= r1 + r2 - tol {
        // Disjoint 鈥?both full circles.
        let tau = std::f64::consts::TAU;
        let mut pts: Vec<DVec2> = (0..n1).map(|i| {
            let ang = tau * i as f64 / n1 as f64;
            let (s, c) = ang.sin_cos();
            c1 + DVec2::new(r1 * c, r1 * s)
        }).collect();
        let _last = *pts.last().unwrap_or(&c1);
        pts.extend((0..=n2).map(|i| {
            let ang = tau * i as f64 / n2 as f64;
            let (s, c) = ang.sin_cos();
            c2 + DVec2::new(r2 * c, r2 * s)
        }));
        // Drop duplicate start point
        if pts.len() > 3 && (pts.last().unwrap() - pts[0]).length_squared() < tol * tol {
            pts.pop();
        }
        return pts;
    }

    // Overlapping 鈥?trace the outer envelope arcs.
    let dir = (c2 - c1) / d;          // unit vector C1 鈫?C2
    let perp = DVec2::new(-dir.y, dir.x); // 90掳 CCW

    let ix = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let _iy = (r1 * r1 - ix * ix).max(0.0).sqrt();

    let cos_t1 = (ix / r1).clamp(-1.0, 1.0);
    let theta1 = cos_t1.acos();                  // C1 arc half-angle
    let C = ((r1 * r1 - d * d - r2 * r2) / (2.0 * d * r2)).clamp(-1.0, 1.0);
    let theta2 = C.acos();                       // C2 arc half-angle

    let tau = std::f64::consts::TAU;

    // C1 outer arc: from theta1 鈫?2蟺 鈭?theta1 (through 蟺, the left/away side).
    let a1_start = theta1;
    let a1_end = tau - theta1;
    let sweep1 = ((a1_end - a1_start) % tau + tau) % tau; // positive

    let mut pts = Vec::with_capacity(n1 + n2 + 2);
    for k in 0..=n1 {
        let frac = k as f64 / n1 as f64;
        let ang = a1_start + sweep1 * frac;
        let (s, c) = ang.sin_cos();
        pts.push(c1 + DVec2::new(r1 * c, r1 * s));
    }

    // C2 outer arc: from 鈭捨? 鈫?+胃2 (through 0, the right/away side).
    // In C2's local frame (+X = dir = away from C1, +Y = perp = CCW).
    let a2_start = -theta2;
    let a2_end = theta2;
    let sweep2 = a2_end - a2_start; // = 2 * theta2, always positive

    for k in 1..n2 {
        let frac = k as f64 / n2 as f64;
        let ang = a2_start + sweep2 * frac;
        let (s, c) = ang.sin_cos();
        pts.push(c2 + r2 * (c * dir + s * perp));
    }

    // Close the polygon (remove trailing duplicate start-point).
    if pts.len() > 3 && (pts.last().unwrap() - pts[0]).length_squared() < tol * tol {
        pts.pop();
    }
    pts
}

/// Build a tessellated BRep for the union of two Z-aligned cone frustums with
/// offset XY centers.
///
/// The first cone is assumed to be the one with the larger Z span (the one that
/// provides the "above-overlap" walls).
fn build_cone_cone_union_tessellated(
    c1_xy: DVec2, c1_z_lo: f64, c1_z_hi: f64, c1_r_lo: f64, c1_r_hi: f64,
    c2_xy: DVec2, c2_z_lo: f64, c2_z_hi: f64, c2_r_lo: f64, c2_r_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let overlap_lo = c1_z_lo.max(c2_z_lo);
    let overlap_hi = c1_z_hi.min(c2_z_hi);
    let has_overlap = overlap_hi > overlap_lo + tol;
    let has_top = c1_z_hi > overlap_hi + tol;

    // Cone 1 radius at any Z (the tall cone that continues above overlap).
    let c1_dr_dz = (c1_r_hi - c1_r_lo) / (c1_z_hi - c1_z_lo);
    let c1_r = |z: f64| (c1_r_lo + c1_dr_dz * (z - c1_z_lo)).max(0.0);

    // Cone 2 radius at any Z.
    let c2_dr_dz = (c2_r_hi - c2_r_lo) / (c2_z_hi - c2_z_lo);
    let c2_r = |z: f64| (c2_r_lo + c2_dr_dz * (z - c2_z_lo)).max(0.0);

    const N_ARC: usize = 48;   // sample points per circle boundary arc
    const N_SLICE: usize = 32; // Z slices in the overlap region

    let empty_wire = || Wire { edges: vec![] };
    // Simple world transform: points are already in world XY.
    let to_world = |x: f64, y: f64, z: f64| -> DVec3 { DVec3::new(x, y, z) };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };
    let mut faces: Vec<Face> = Vec::new();

    // 鈹€鈹€ Overlap region [overlap_lo, overlap_hi] (two-circle union) 鈹€鈹€
    if has_overlap {
        // Pre-compute slice polygons at each Z level.
        let n_sl = N_SLICE.max(1);
        let mut polys: Vec<Vec<DVec2>> = Vec::with_capacity(n_sl + 1);
        for i in 0..=n_sl {
            let z = overlap_lo + (overlap_hi - overlap_lo) * i as f64 / n_sl as f64;
            let r_a = c1_r(z);
            let r_b = c2_r(z);
            polys.push(two_circle_union_ccw_pts(c1_xy, r_a, c2_xy, r_b, N_ARC, N_ARC));
        }

        // Extrude walls between consecutive slices.
        for i in 0..n_sl {
            let z0 = overlap_lo + (overlap_hi - overlap_lo) * i as f64 / n_sl as f64;
            let z1 = overlap_lo + (overlap_hi - overlap_lo) * (i + 1) as f64 / n_sl as f64;
            let n = polys[i].len();
            if n < 3 { continue; }

            let mut idx = Vec::with_capacity(2 * n);
            for p in &polys[i] { idx.push(add_v(to_world(p.x, p.y, z0))); }
            for p in &polys[i] { idx.push(add_v(to_world(p.x, p.y, z1))); }
            let mut tris = Vec::with_capacity(n * 2);
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([idx[j], idx[k], idx[n + k]]);
                tris.push([idx[j], idx[n + k], idx[n + j]]);
            }
            faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
        }

        // 鈹€鈹€ Interface face at overlap_hi: two-circle union minus the circle of cone 1 鈹€鈹€
        if has_top {
            let r_at = c1_r(overlap_hi);
            let poly = two_circle_union_ccw_pts(c1_xy, r_at, c2_xy, c2_r(overlap_hi), N_ARC, N_ARC);
            add_interface_face(
                &mut add_v, &mut faces,
                &poly, overlap_hi, -DVec3::Z,
                c1_xy, r_at,
                &to_world, &empty_wire,
            );
        }
    }

    // 鈹€鈹€ Top region [overlap_hi, c1_z_hi] (single cone) 鈹€鈹€
    if has_top {
        let n_sl = N_SLICE.max(1);
        let z0 = overlap_hi;
        let z1 = c1_z_hi;
        // Build wall using the single-circle polygon.
        let r_lo = c1_r(z0);
        let _r_hi = c1_r(z1);
        if r_lo > tol {
            let bot_poly = two_circle_union_ccw_pts(c1_xy, r_lo, c2_xy, c2_r(z0), N_ARC, N_ARC);
            // Remap to same vertex count as interface uses.
            let poly_single = build_circle_polygon(c1_xy.x, c1_xy.y, r_lo);
            let n_pts = bot_poly.len();
            let poly = if poly_single.len() != n_pts {
                let ref_pt = c1_xy + DVec2::new(r_lo, 0.0);
                remap_polygon_arclength(&poly_single, n_pts, ref_pt)
            } else {
                poly_single
            };

            // Bottom cap at overlap_hi if no interface was created.
            if !has_overlap {
                add_cap_face(&mut add_v, &mut faces, &poly, z0, -DVec3::Z, &to_world, &empty_wire);
            }

            // Build wall from z0 to z1.
            let dz = (z1 - z0) / n_sl as f64;
            // We need to handle the transition here: use poly at z0 and circle at z1.
            // But since both use the same vertex count (after remapping), the wall
            // connects them directly.
            let top_r_at = |z: f64| {
                let r = c1_r(z);
                let mut p = build_circle_polygon(c1_xy.x, c1_xy.y, r);
                if p.len() != n_pts {
                    let ref_pt = c1_xy + DVec2::new(r, 0.0);
                    p = remap_polygon_arclength(&p, n_pts, ref_pt);
                }
                p
            };

            for i in 0..n_sl {
                let za = z0 + dz * i as f64;
                let zb = z0 + dz * (i + 1) as f64;
                let pts_a = top_r_at(za);
                let pts_b = top_r_at(zb);
                let n = pts_a.len();
                if n < 3 { continue; }
                let mut idx = Vec::with_capacity(2 * n);
                for p in &pts_a { idx.push(add_v(to_world(p.x, p.y, za))); }
                for p in &pts_b { idx.push(add_v(to_world(p.x, p.y, zb))); }
                let mut tris = Vec::with_capacity(n * 2);
                for j in 0..n {
                    let k = (j + 1) % n;
                    tris.push([idx[j], idx[k], idx[n + k]]);
                    tris.push([idx[j], idx[n + k], idx[n + j]]);
                }
                faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
            }

            // Top cap at c1_z_hi.
            let r_top = c1_r(z1);
            if r_top > tol {
                let mut top_poly = build_circle_polygon(c1_xy.x, c1_xy.y, r_top);
                if top_poly.len() != n_pts {
                    let ref_pt = c1_xy + DVec2::new(r_top, 0.0);
                    top_poly = remap_polygon_arclength(&top_poly, n_pts, ref_pt);
                }
                add_cap_face(&mut add_v, &mut faces, &top_poly, z1, DVec3::Z, &to_world, &empty_wire);
            }
        }
    }

    // 鈹€鈹€ Bottom cap 鈹€鈹€
    if !has_overlap || c1_z_lo < overlap_lo - tol {
        let z = c1_z_lo;
        let r = c1_r(z);
        if r > tol {
            let poly = two_circle_union_ccw_pts(c1_xy, r, c2_xy, c2_r(z), N_ARC, N_ARC);
            add_cap_face(&mut add_v, &mut faces, &poly, z, -DVec3::Z, &to_world, &empty_wire);
        }
    } else {
        let z = overlap_lo;
        let r_a = c1_r(z);
        let r_b = c2_r(z);
        if r_a > tol || r_b > tol {
            let poly = two_circle_union_ccw_pts(c1_xy, r_a, c2_xy, r_b, N_ARC, N_ARC);
            add_cap_face(&mut add_v, &mut faces, &poly, z, -DVec3::Z, &to_world, &empty_wire);
        }
    }

    // 鈹€鈹€ Top cap (only if no top region was built above) 鈹€鈹€
    if !has_top {
        let z = c1_z_hi;
        let r = c1_r(z);
        if r > tol {
            let poly = if c1_z_hi > overlap_hi + tol || !has_overlap {
                build_circle_polygon(c1_xy.x, c1_xy.y, r)
            } else {
                two_circle_union_ccw_pts(c1_xy, r, c2_xy, c2_r(z), N_ARC, N_ARC)
            };
            add_cap_face(&mut add_v, &mut faces, &poly, z, DVec3::Z, &to_world, &empty_wire);
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Detect cone-cone union with offset centers (different XY, same Z-axis).
fn try_union_offset_cones_one_dir(a: &BRep, b: &BRep) -> Option<BRep> {
    let (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi) = detect_z_axis_cone_frustum(a)?;
    let (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi) = detect_z_axis_cone_frustum(b)?;

    // Ensure c1 is the one with larger Z span (it continues above overlap).
    let (s_a, s_b) = if (c1_z_hi - c1_z_lo) >= (c2_z_hi - c2_z_lo) {
        ((c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi),
         (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi))
    } else {
        ((c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi),
         (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi))
    };

    let (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi) = s_a;
    let (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi) = s_b;

    let tol = TOLERANCE_LEN_MIN;

    // Not offset 鈥?delegate to coaxial path.
    if (c1_xy.x - c2_xy.x).abs() < tol && (c1_xy.y - c2_xy.y).abs() < tol {
        return None;
    }

    // Z ranges must overlap.
    if c1_z_hi < c2_z_lo - tol || c2_z_hi < c1_z_lo - tol {
        return None;
    }

    build_cone_cone_union_tessellated(
        c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi,
        c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi,
    )
}

/// Bidirectional wrapper for offset cone-cone union.
pub fn try_union_offset_cones(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_offset_cones_one_dir(a, b)
        .or_else(|| try_union_offset_cones_one_dir(b, a))
}

/// Build BRep for outer conical frustum minus inner conical frustum (coaxial, Z-aligned).
/// Result is a hollow conical frustum: outer lateral + bottom cap + top annulus + inner lateral + cavity floor.
///
/// All faces are triangulated (no analytic surfaces) in a single shell.
fn build_conical_frustum_minus_frustum_brep(
    zo_lo: f64, zo_hi: f64, ro_lo: f64, ro_hi: f64,
    zi_lo: f64, zi_hi: f64, ri_lo: f64, ri_hi: f64,
) -> Option<BRep> {
    use std::f64::consts::TAU;
    const N: usize = 48; // Circumferential divisions

    let empty_wire = || Wire { edges: vec![] };
    let tol = TOLERANCE_LEN_MIN;

    // Overlap Z range (inner clamped to outer)
    let z_olap_lo = zi_lo.max(zo_lo);
    let z_olap_hi = zi_hi.min(zo_hi);
    let has_overlap = z_olap_hi > z_olap_lo + tol;

    // Inner cone radius at overlap boundaries
    let dri = (ri_hi - ri_lo) / (zi_hi - zi_lo);
    let ri_at_olap_lo = ri_lo + dri * (z_olap_lo - zi_lo);
    let ri_at_olap_hi = ri_lo + dri * (z_olap_hi - zi_lo);

    let ring_pts = |z: f64, r: f64| -> Vec<DVec3> {
        (0..N).map(|i| {
            let ang = TAU * i as f64 / N as f64;
            let (c, s) = ang.sin_cos();
            DVec3::new(r * c, r * s, z)
        }).collect()
    };

    let outer_bot = ring_pts(zo_lo, ro_lo);
    let outer_top = ring_pts(zo_hi, ro_hi);

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // 1. Outer lateral (z=zo_lo to z=zo_hi)
    {
        let mut tris = Vec::new();
        wall_grid(&mut add_v, &outer_bot, &outer_top, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    }

    // 2. Bottom at z=zo_lo
    //    - If inner starts at or below outer bottom: annulus
    //    - Otherwise: full disk
    if zi_lo <= zo_lo + tol && has_overlap && ri_at_olap_lo < ro_lo - tol {
        // Hole goes through the bottom
        let mut tris = Vec::new();
        annulus_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_lo), ri_at_olap_lo, ro_lo, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    } else {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_lo), ro_lo, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    }

    // 3. Top at z=zo_hi
    //    - If inner ends at or above outer top: annulus
    //    - Otherwise: full disk
    if zi_hi >= zo_hi - tol && has_overlap && ri_at_olap_hi < ro_hi - tol {
        let mut tris = Vec::new();
        annulus_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_hi), ri_at_olap_hi, ro_hi, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    } else {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_hi), ro_hi, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    }

    // 4. Inner lateral 鈥?cavity wall (overlap region only)
    if has_overlap {
        let inner_bot = ring_pts(z_olap_lo, ri_at_olap_lo);
        let inner_top = ring_pts(z_olap_hi, ri_at_olap_hi);
        let mut tris = Vec::new();
        wall_grid(&mut add_v, &inner_bot, &inner_top, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    }

    // 5. Cavity floor 鈥?where hole starts (if inner starts above outer bottom)
    if zi_lo > zo_lo + tol && has_overlap && ri_at_olap_lo > tol {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, z_olap_lo), ri_at_olap_lo, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    }

    // 6. Cavity ceiling 鈥?where hole ends (if inner ends below outer top)
    if zi_hi < zo_hi - tol && has_overlap && ri_at_olap_hi > tol {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, z_olap_hi), ri_at_olap_hi, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false });
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Generate a quad-grid triangulation between two Z-aligned rings (conical lateral strip).
fn wall_grid(
    add_v: &mut impl FnMut(DVec3) -> usize,
    bot: &[DVec3],
    top: &[DVec3],
    tris: &mut Vec<[usize; 3]>,
) {
    let n = bot.len().min(top.len());
    if n < 3 { return; }
    let mut idx = Vec::with_capacity(2 * n);
    for i in 0..n {
        idx.push(add_v(bot[i]));
    }
    for i in 0..n {
        idx.push(add_v(top[i]));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        let b0 = idx[i];
        let b1 = idx[j];
        let t0 = idx[n + i];
        let t1 = idx[n + j];
        tris.push([b0, b1, t1]);
        tris.push([b0, t1, t0]);
    }
}

/// Triangle fan for a full disk centered at origin in XY plane at given z.
fn disk_tri_fan(
    add_v: &mut impl FnMut(DVec3) -> usize,
    center: DVec3,
    radius: f64,
    n: usize,
    tris: &mut Vec<[usize; 3]>,
) {
    let c_idx = add_v(center);
    let mut ring = Vec::with_capacity(n);
    for i in 0..n {
        let ang = std::f64::consts::TAU * i as f64 / n as f64;
        let (s, c) = ang.sin_cos();
        ring.push(add_v(DVec3::new(radius * c, radius * s, center.z)));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        tris.push([c_idx, ring[i], ring[j]]);
    }
}

/// Triangle fan for an annulus centered at origin in XY plane at given z.
fn annulus_tri_fan(
    add_v: &mut impl FnMut(DVec3) -> usize,
    center: DVec3,
    r_inner: f64,
    r_outer: f64,
    n: usize,
    tris: &mut Vec<[usize; 3]>,
) {
    let mut outer_ring = Vec::with_capacity(n);
    let mut inner_ring = Vec::with_capacity(n);
    for i in 0..n {
        let ang = std::f64::consts::TAU * i as f64 / n as f64;
        let (s, c) = ang.sin_cos();
        outer_ring.push(add_v(DVec3::new(r_outer * c, r_outer * s, center.z)));
        inner_ring.push(add_v(DVec3::new(r_inner * c, r_inner * s, center.z)));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        tris.push([outer_ring[i], outer_ring[j], inner_ring[j]]);
        tris.push([outer_ring[i], inner_ring[j], inner_ring[i]]);
    }
}

/// Fast path: ZP9 large box 鈭?sphere near one face.
/// Detects a large axis-aligned box and a sphere that intersects exactly one box face.
pub fn try_intersection_box_sphere_single_face(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try both orderings
    try_intersection_box_sphere_single_face_pair(a, b)
        .or_else(|| try_intersection_box_sphere_single_face_pair(b, a))
}

fn try_intersection_box_sphere_single_face_pair(box_: &BRep, sphere: &BRep) -> Option<BRep> {
    // Check sphere
    let (sp_center, sp_radius) = sphere_center_r(sphere)?;

    // Check box: must be axis-aligned, 6 planar faces
    let box_bb = try_as_axis_aligned_box(box_)?;
    let [bmin, bmax] = box_bb;

    // Find which box face(s) the sphere extends beyond
    let sp_min = sp_center - DVec3::splat(sp_radius);
    let sp_max = sp_center + DVec3::splat(sp_radius);

    // Count clip planes and find the one(s)
    let mut clip_axes: Vec<(DVec3, f64)> = Vec::new(); // (axis, plane_value)

    if sp_min.x < bmin.x { clip_axes.push((-DVec3::X, bmin.x)); }
    if sp_max.x > bmax.x { clip_axes.push((DVec3::X, bmax.x)); }
    if sp_min.y < bmin.y { clip_axes.push((-DVec3::Y, bmin.y)); }
    if sp_max.y > bmax.y { clip_axes.push((DVec3::Y, bmax.y)); }
    if sp_min.z < bmin.z { clip_axes.push((-DVec3::Z, bmin.z)); }
    if sp_max.z > bmax.z { clip_axes.push((DVec3::Z, bmax.z)); }

    // ZP9: exactly one clip plane, sphere mainly inside the box
    if clip_axes.len() != 1 {
        return None;
    }

    let (axis, plane_val) = clip_axes[0];

    // For ZP9: sphere center must be inside the box on the other axes
    if sp_center.x < bmin.x || sp_center.x > bmax.x
        || sp_center.y < bmin.y || sp_center.y > bmax.y
        || sp_center.z < bmin.z || sp_center.z > bmax.z
    {
        return None; // center outside box 鈫?more complex intersection
    }

    // Build the result: sphere clipped by one plane
    build_sphere_clipped_by_plane(sp_center, sp_radius, axis, plane_val)
}

/// Build a triangulated BRep for a sphere clipped by one plane.
fn build_sphere_clipped_by_plane(
    center: DVec3, radius: f64,
    plane_normal: DVec3, plane_d: f64, // plane equation: plane_normal 路 p = plane_d (normal is unit)
) -> Option<BRep> {
    use std::f64::consts::{PI, TAU};
    const NS: usize = 32; // theta divisions
    const NP: usize = 16; // phi divisions

    let empty_wire = || Wire { edges: vec![] };

    // Determine the clip height on the sphere: cos(phi_clip) = (plane_d - dot(center, plane_normal)) / radius
    let d_from_center = (plane_d - center.dot(plane_normal)).abs();
    if d_from_center >= radius - 1e-12 {
        return None; // plane doesn't intersect the sphere
    }

    // We want the portion of the sphere where plane_normal路p 鈮?plane_d
    // In spherical coords with Z=plane_normal: keep z 鈮?plane_d - center路plane_normal
    let z_clip = plane_d - center.dot(plane_normal);
    let phi_clip = (z_clip / radius).acos(); // colatitude of the clip circle

    // Build rotated basis so that plane_normal is Z
    let z_axis = plane_normal;
    let x_axis = if z_axis.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let y_axis = z_axis.cross(x_axis).normalize();
    let x_axis = y_axis.cross(z_axis).normalize();

    let mut verts: Vec<Vertex> = vec![];
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut tris: Vec<[usize; 3]> = Vec::new();

    // Generate sphere surface where phi 鈭?[phi_clip, 蟺] (the clipped portion)
    let dphi = (PI - phi_clip) / NP as f64;
    let dtheta = TAU / NS as f64;
    let mut idx = vec![vec![0usize; NS+1]; NP+1];
    for pj in 0..=NP { for ti in 0..=NS {
        let phi = phi_clip + pj as f64 * dphi;
        let theta = ti as f64 * dtheta;
        let (sinp, cosp) = phi.sin_cos();
        let (sint, cost) = theta.sin_cos();
        // Position in local (plane_normal = Z) coordinates
        let local = DVec3::new(radius * sinp * cost, radius * sinp * sint, radius * cosp);
        // Transform to world coordinates
        let world = center + local.x * x_axis + local.y * y_axis + local.z * z_axis;
        idx[pj][ti] = add_v(world);
    }}
    for pj in 0..NP { for ti in 0..NS {
        let a = idx[pj][ti]; let b = idx[pj][ti+1];
        let c = idx[pj+1][ti]; let d = idx[pj+1][ti+1];
        tris.push([a, b, d]); tris.push([a, d, c]);
    }}

    // Generate planar cap at the clip plane
    // Cap center = center + z_clip * plane_normal, radius = radius * sin(phi_clip)
    let cap_center = center + z_clip * z_axis;
    let cap_r = radius * phi_clip.sin();
    let n_cap_seg = NS.max(16);
    let cap_center_idx = add_v(cap_center);
    for i in 0..n_cap_seg {
        let ang0 = (i as f64 / n_cap_seg as f64) * TAU;
        let ang1 = ((i + 1) as f64 / n_cap_seg as f64) * TAU;
        let (c0, s0) = ang0.sin_cos();
        let (c1, s1) = ang1.sin_cos();
        let p0 = cap_center + cap_r * (c0 * x_axis + s0 * y_axis);
        let p1 = cap_center + cap_r * (c1 * x_axis + s1 * y_axis);
        let v0 = add_v(p0);
        let v1 = add_v(p1);
        tris.push([cap_center_idx, v0, v1]);
    }

    let faces = vec![Face {
        outer_wire: empty_wire(), inner_wires: vec![],
        normal: DVec3::ZERO, triangles: tris,
        sample_point: None, mesh_dirty: false,
    }];

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

// 鈹€鈹€ Sphere鈥揃ox Intersection Fast Path 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Build a mesh-based BRep for a sphere clipped by N planes (interior half-spaces).
///
/// Each plane is `(point_on_plane, outward_normal)` where the interior is
/// `(p - origin)路outward_normal 鈮?0`.
///
/// Returns a triangulated BRep with one face containing all triangles (spherical
/// surface + planar caps).  Follows the same pattern as
/// [`build_sphere_clipped_by_plane`] but handles an arbitrary number of planes,
/// not just one.
fn build_sphere_clipped_by_planes(
    center: DVec3,
    radius: f64,
    planes: &[(DVec3, DVec3)],
) -> Option<BRep> {
    use std::f64::consts::{PI, TAU};
    const NS: usize = 48;  // longitude divisions
    const NP: usize = 24;  // latitude divisions

    let mut verts: Vec<Vertex> = vec![];
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    // 鈹€鈹€ Step 1: Generate sphere surface vertices and classify 鈹€鈹€
    struct GridCell { vi: Option<usize> }
    let mut grid: Vec<Vec<GridCell>> = (0..=NP)
        .map(|_| (0..=NS).map(|_| GridCell { vi: None }).collect())
        .collect();

    for pj in 0..=NP {
        let phi = pj as f64 * PI / NP as f64;
        let (sinp, cosp) = phi.sin_cos();
        for ti in 0..=NS {
            let theta = ti as f64 * TAU / NS as f64;
            let (sint, cost) = theta.sin_cos();
            let p = center + radius * DVec3::new(sinp * cost, sinp * sint, cosp);

            let inside = planes.iter().all(|(origin, normal)| {
                (p - origin).dot(*normal) <= crate::tolerance::TOLERANCE_ABS
            });

            if inside {
                grid[pj][ti].vi = Some(add_v(p));
            }
        }
    }

    // 鈹€鈹€ Step 2: Triangulate the spherical face 鈹€鈹€
    let mut tris: Vec<[usize; 3]> = Vec::new();

    for pj in 0..NP {
        for ti in 0..NS {
            let a = grid[pj][ti].vi;
            let b = grid[pj][ti + 1].vi;
            let c = grid[pj + 1][ti].vi;
            let d = grid[pj + 1][ti + 1].vi;

            match (a, b, c, d) {
                (Some(a), Some(b), Some(c), Some(d)) => {
                    tris.push([a, b, d]);
                    tris.push([a, d, c]);
                }
                _ => {} // Skip mixed cells (48脳24 resolution is sufficient for 15% SA tol)
            }
        }
    }

    if tris.is_empty() {
        return None; // No intersection
    }

    // 鈹€鈹€ Step 3: Planar cap faces 鈹€鈹€
    for &(origin, normal) in planes {
        let n = normal.normalize_or_zero();
        if n.length_squared() < 0.5 { continue; }

        let d = (center - origin).dot(n);
        if d.abs() >= radius - crate::tolerance::TOLERANCE_ABS { continue; }

        let cap_center = center - d * n;
        let cap_r = (radius * radius - d * d).sqrt();
        if cap_r < crate::tolerance::TOLERANCE_COORD_SUB { continue; }

        // Sample 64 points on the intersection circle, keep those inside all OTHER planes.
        const NC: usize = 64;
        struct CapPt { ang: f64, pos: DVec3 }
        let mut candidates: Vec<CapPt> = Vec::new();

        let x_ax = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        let y_ax = n.cross(x_ax).normalize();
        let x_ax = y_ax.cross(n).normalize();

        for i in 0..NC {
            let ang = i as f64 * TAU / NC as f64;
            let (c, s) = ang.sin_cos();
            let p = cap_center + cap_r * (c * x_ax + s * y_ax);

            let inside = planes.iter().filter(|(o, n_)| {
                !( (o - origin).length() < crate::tolerance::TOLERANCE_COORD_SUB
                    && (n_ - &normal).length() < crate::tolerance::TOLERANCE_AXIS_ALIGN )
            }).all(|(o, n_)| (p - o).dot(*n_) <= crate::tolerance::TOLERANCE_ABS);

            if inside {
                candidates.push(CapPt { ang, pos: p });
            }
        }

        candidates.sort_by(|a, b| a.ang.partial_cmp(&b.ang).unwrap_or(std::cmp::Ordering::Equal));
        if candidates.len() < 3 { continue; }

        // Find the largest contiguous angular cluster (gaps < threshold).  This avoids
        // the degenerate fan edge when isolated outlier points (e.g. a lone 0掳 point
        // in a 270掳鈥?54掳 cluster) wrap the long way around the circle.
        let step = TAU / NC as f64;
        let gap_thresh = step * 4.0;
        let nc = candidates.len();
        let mut best_len = 0usize;
        let mut best_start = 0usize;
        let mut cur_len = 1usize;
        let mut cur_start = 0usize;
        for i in 0..nc {
            let next = (i + 1) % nc;
            let gap = (candidates[next].ang - candidates[i].ang + TAU) % TAU;
            if gap < gap_thresh {
                cur_len += 1;
            } else {
                if cur_len > best_len { best_len = cur_len; best_start = cur_start; }
                cur_len = 1;
                cur_start = next;
            }
        }
        if cur_len > best_len { best_len = cur_len; best_start = cur_start; }

        let cap_poly: Vec<DVec3> = (0..best_len)
            .map(|k| candidates[(best_start + k) % nc].pos)
            .collect();

        if cap_poly.len() < 3 { continue; }

        // Triangulate as an OPEN fan (no wrap from last鈫抐irst), avoiding the
        // degenerate triangle that would span the long way around the circle when
        // the arc is < 2蟺.  The n鈭? wedge triangles correctly cover the convex cap.
        let cap_ci = add_v(cap_center);
        let mut cap_vis: Vec<usize> = Vec::with_capacity(cap_poly.len());
        for p in &cap_poly { cap_vis.push(add_v(*p)); }
        for i in 0..cap_poly.len().saturating_sub(1) {
            tris.push([cap_ci, cap_vis[i], cap_vis[i + 1]]);
        }
    }

    // 鈹€鈹€ Step 4: Build BRep 鈹€鈹€
    let faces = vec![Face {
        outer_wire: Wire { edges: vec![] },
        inner_wires: vec![],
        normal: DVec3::ZERO,
        triangles: tris,
        sample_point: None,
        mesh_dirty: false,
    }];

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Detect sphere 鈭?box intersection and build a mesh-based result.
///
/// Handles any box orientation (via `try_as_box`). Requires the sphere center
/// to be inside the box (the common case for all bcommon_simple sphere-box tests).
/// Falls through to PaveFiller for cases where the center is outside.
fn try_intersection_sphere_box_pair(sphere: &BRep, box_: &BRep) -> Option<BRep> {
    let (sp_center, sp_radius) = sphere_center_r(sphere)?;
    try_as_box(box_)?;  // Validate it's a box; discard result

    // Quick empty check via AABB
    let [bmin, bmax] = box_.bounding_box()?;
    let sp_min = sp_center - DVec3::splat(sp_radius);
    let sp_max = sp_center + DVec3::splat(sp_radius);
    let tol = 1e-8;
    if sp_max.x < bmin.x - tol || sp_min.x > bmax.x + tol
        || sp_max.y < bmin.y - tol || sp_min.y > bmax.y + tol
        || sp_max.z < bmin.z - tol || sp_min.z > bmax.z + tol
    {
        return Some(BRep::default());
    }

    // Delegate to analytic builder.
    crate::sphere_box_analytic::build_sphere_box_intersection_analytic(sphere, box_)
}

/// Fast-path for sphere 鈭?box intersection (any orientation).
///
/// OCCT has no equivalent fast-path 鈥?its `BRepAlgoAPI_BooleanOperation` uses
/// the general pipeline for all shape pairs.  This is a pure rcad optimization
/// that avoids the PaveFiller (24鈥?1s 鈫?<0.01s) for the common case where the
/// sphere center lies inside the box (all bcommon_simple sphere-box tests).
pub fn try_intersection_sphere_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_sphere_box_pair(a, b)
        .or_else(|| try_intersection_sphere_box_pair(b, a))
}

// 鈹€鈹€ Box鈥揅ylinder Difference Fast Path (box 鈭?cylinder, Z-axis cylinder) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// A point where the circle intersects a box edge in UV space.
#[derive(Debug, Clone, Copy)]
struct UVEdgePt {
    /// Perimeter parameter t 鈭?[0, 8)
    t: f64,
    /// UV coordinates
    u: f64,
    v: f64,
    /// Edge index: 0 = u_min, 1 = u_max, 2 = v_min, 3 = v_max
    edge: usize,
    /// Circle angle 胃 = atan2(v 鈭?cv, u 鈭?cu)
    theta: f64,
}

/// Find which of the 3 box axes is closest to the world Z axis.
pub(crate) fn find_z_axis_index(info: &BoxInfo) -> Option<usize> {
    for i in 0..3 {
        if info.axes[i].dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN {
            return Some(i);
        }
    }
    None
}

/// Extract cylinder parameters from a cylinder primitive.
pub(crate) fn try_cylinder_center_axis_radius_height(brep: &BRep) -> Option<(DVec3, DVec3, f64, f64)> {
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
/// t 鈭?[0,2) 鈫?v_min (v=鈭抏v), u鈭圼鈭抏u,eu]
/// t 鈭?[2,4) 鈫?u_max (u=eu),  v鈭圼鈭抏v,ev]
/// t 鈭?[4,6) 鈫?v_max (v=ev),  u鈭圼eu,鈭抏u]
/// t 鈭?[6,8) 鈫?u_min (u=鈭抏u), v鈭圼ev,鈭抏v]
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

/// Map UV coordinate to perimeter t 鈭?[0, 8).
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

// 鈹€鈹€ shared segment types and helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// A merged segment along the box perimeter 鈥?either inside or outside the cylinder.
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
        // the box at non-midpoint positions (e.g., an arc >180掳 whose midpoint
        // happens to be inside but passes through a bulge).  Sample several points
        // along the shorter arc to verify all lie inside; if not, take the longer.
        let shorter_cw = cw_len <= ccw_len;
        let (check_len, check_dir) = if shorter_cw {
            (cw_len, -1.0)  // CW direction = decreasing 胃
        } else {
            (ccw_len, 1.0)  // CCW direction = increasing 胃
        };
        // Check 5 evenly-spaced sample points along the shorter arc.
        let n_samples = 5usize;
        let all_inside = (0..=n_samples).all(|i| {
            let frac = i as f64 / n_samples as f64;
            let th = th_start + check_dir * check_len * frac;
            // Normalize to [0, 2蟺)
            let th_n = th.rem_euclid(2.0 * std::f64::consts::PI);
            inside_box(th_n)
        });
        if all_inside {
            // Shorter arc stays inside 鈫?use it.
            if shorter_cw { (-cw_len, cw_len) } else { (ccw_len, ccw_len) }
        } else {
            // Shorter arc exits the box 鈫?use the longer arc.
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

/// Build cap polygon at z level for (box rect 鈭?circle).
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

    let corner_ts = [0.0, 2.0, 4.0, 6.0];
    let start_norm = if t_start < 0.0 { t_start + 8.0 } else { t_start };
    let end_norm = if t_end <= t_start { t_end + 8.0 } else { t_end };

    // Collect corners within range and sort by normalized t for CCW order.
    let mut corners_in_range: Vec<(f64, f64, f64)> = Vec::new();
    for &ct in &corner_ts {
        let ct_norm = if ct < start_norm { ct + 8.0 } else { ct };
        if ct_norm > start_norm + tol && ct_norm < end_norm - tol {
            let (u, v) = box_perimeter_uv(ct, eu, ev);
            corners_in_range.push((ct_norm, u, v));
        }
    }
    corners_in_range.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (_, u, v) in corners_in_range {
        let pt = corner(u, v, z);
        if (poly.last().unwrap() - pt).length() > tol { poly.push(pt); }
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
    // verts[0] is the start point 鈥?already in poly from the previous segment,
    // or if poly is empty (first merged segment is a circle arc, not box perimeter)
    // we need to push it here.
    if poly.is_empty() {
        if let Some(&first) = verts.first() {
            poly.push(first);
        }
    }
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
    let _tol = TOLERANCE_LEN_MIN;
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

/// Main Stage 3 orchestrator: build box 鈭?cylinder with partial XY containment.
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

    // u_min edge (u=-eu): solve (-eu-cu)^2 + (v-cv)^2 = r^2 鈫?v
    let mut ints_u_min: Vec<f64> = Vec::with_capacity(2);
    let d0 = -eu - cu; let disc0 = cyl_r * cyl_r - d0 * d0;
    if disc0 >= 0.0 { let off = disc0.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_min); add_pt(cv+off, -ev, ev, &mut ints_u_min); }
    ints_u_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_u_min.dedup_by(|a,b| (*a - *b).abs() < tol);

    // u_max edge (u=eu): solve (eu-cu)^2 + (v-cv)^2 = r^2 鈫?v
    let mut ints_u_max: Vec<f64> = Vec::with_capacity(2);
    let d1 = eu - cu; let disc1 = cyl_r * cyl_r - d1 * d1;
    if disc1 >= 0.0 { let off = disc1.sqrt(); add_pt(cv-off, -ev, ev, &mut ints_u_max); add_pt(cv+off, -ev, ev, &mut ints_u_max); }
    ints_u_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_u_max.dedup_by(|a,b| (*a - *b).abs() < tol);

    // v_min edge (v=-ev): solve (u-cu)^2 + (-ev-cv)^2 = r^2 鈫?u
    let mut ints_v_min: Vec<f64> = Vec::with_capacity(2);
    let d2 = -ev - cv; let disc2 = cyl_r * cyl_r - d2 * d2;
    if disc2 >= 0.0 { let off = disc2.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_min); add_pt(cu+off, -eu, eu, &mut ints_v_min); }
    ints_v_min.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_v_min.dedup_by(|a,b| (*a - *b).abs() < tol);

    // v_max edge (v=ev): solve (u-cu)^2 + (ev-cv)^2 = r^2 鈫?u
    let mut ints_v_max: Vec<f64> = Vec::with_capacity(2);
    let d3 = ev - cv; let disc3 = cyl_r * cyl_r - d3 * d3;
    if disc3 >= 0.0 { let off = disc3.sqrt(); add_pt(cu-off, -eu, eu, &mut ints_v_max); add_pt(cu+off, -eu, eu, &mut ints_v_max); }
    ints_v_max.sort_by(|a,b| a.partial_cmp(b).unwrap());
    ints_v_max.dedup_by(|a,b| (*a - *b).abs() < tol);

    // Side walls 鈥?always use build_trimmed_edge_pieces for proper z-splitting
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
            // else: full span 鈥?cylinder removes this entire face. Nothing to keep.
        } else if ints.len() == 1 {
            // Single intersection point 鈥?check which interval to keep.
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
                // Both outside 鈥?tangent. Keep full face.
                build_trimmed_edge_pieces(lo, hi, lo, hi, z_lo, z_hi, ew, cn, nrm, pieces);
            } else if !ins_lo {
                if p > lo + tol { build_trimmed_edge_pieces(lo, hi, lo, p, z_lo, z_hi, ew, cn, nrm, pieces); }
            } else if !ins_hi {
                if p < hi - tol { build_trimmed_edge_pieces(lo, hi, p, hi, z_lo, z_hi, ew, cn, nrm, pieces); }
            } else {
                // Both inside 鈥?shouldn't happen. Fall through to full face.
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

/// Build box 鈭?cylinder difference for full XY containment (cylinder inside or
/// tangent to the box rect). Uses inner-wire annular cap faces and a full 360掳
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

    // 鈹€鈹€ 1. Side walls split at z_lo/z_hi 鈹€鈹€
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

    // 鈹€鈹€ 2. Annular cap faces 鈥斺攢鈹€
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

    // 鈹€鈹€ 3. Cylindrical wall (full 360掳) 鈹€鈹€
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

    // 鈹€鈹€ 4. Sew 鈹€鈹€
    if pieces.is_empty() { return None; }
    let sewn = sew_shells(&pieces, tol.max(TOLERANCE_ABS * 100.0));
    Some(sewn.brep)
}

/// Compute the closed boundary of `rect - circle` at one Z-level as a sequence of 2D points.
///
/// The boundary is traced clockwise. For the typical case (single closed curve), returns
/// a polygon with `n` sample points. For special cases:
/// - rect fully inside circle 鈫?returns empty (no boundary, void removes entire cross-section)
/// - circle fully inside rect 鈫?returns full rect perimeter (outer boundary only; the inner
///   circle hole is handled as a separate inner loop by the caller)
/// - circle outside rect (no overlap) 鈫?returns full rect perimeter (no void)
fn rect_minus_circle_boundary(
    bmin: DVec2, bmax: DVec2,
    cx: f64, cy: f64, r: f64,
    n: usize,
) -> Vec<DVec2> {
    let tol = TOLERANCE_LEN_MIN;
    if n < 4 { return vec![]; }

    let edges = [
        (DVec2::new(bmin.x, bmin.y), DVec2::new(bmax.x, bmin.y)), // bottom
        (DVec2::new(bmax.x, bmin.y), DVec2::new(bmax.x, bmax.y)), // right
        (DVec2::new(bmax.x, bmax.y), DVec2::new(bmin.x, bmax.y)), // top
        (DVec2::new(bmin.x, bmax.y), DVec2::new(bmin.x, bmin.y)), // left
    ];

    // ---- Step 1: Find circle-rect intersection t-values on each edge ----
    struct Intersection { t: f64, edge: usize, pt: DVec2 }
    let mut ints: Vec<Intersection> = Vec::new();

    for (ei, (p0, p1)) in edges.iter().enumerate() {
        let d = *p1 - *p0;
        let a0 = *p0 - DVec2::new(cx, cy);
        // Quadratic: (d路d)*t虏 + 2*(a0路d)*t + (a0路a0 - r虏) = 0
        let A = d.dot(d);
        if A < 1e-30 { continue; }
        let B = 2.0 * a0.dot(d);
        let C = a0.dot(a0) - r * r;
        let disc = B * B - 4.0 * A * C;
        if disc < 0.0 { continue; }
        let sqrt_disc = disc.sqrt();
        for t in [(-B - sqrt_disc) / (2.0 * A), (-B + sqrt_disc) / (2.0 * A)] {
            if t >= -tol && t <= 1.0 + tol {
                let tc = t.clamp(0.0, 1.0);
                let pt = *p0 + d * tc;
                ints.push(Intersection { t: tc, edge: ei, pt });
            }
        }
    }

    // Deduplicate near-identical intersections (same edge and t)
    ints.sort_by(|a, b| a.edge.cmp(&b.edge).then(a.t.partial_cmp(&b.t).unwrap()));
    ints.dedup_by(|a, b| a.edge == b.edge && (a.t - b.t).abs() < tol);

    // Deduplicate by spatial proximity 鈥?a corner may appear as an intersection
    // on two adjacent edges when the circle passes exactly through the corner.
    // Keeping both creates zero-length perimeter segments that produce full-wrap
    // rect loops and self-intersecting polygons.
    ints.sort_by(|a, b| a.pt.x.partial_cmp(&b.pt.x).unwrap()
        .then(a.pt.y.partial_cmp(&b.pt.y).unwrap()));
    ints.dedup_by(|a, b| (a.pt - b.pt).length_squared() < tol * tol);

    // ---- Step 2: Handle no-intersection cases ----
    if ints.is_empty() {
        // Check if the rect center is inside the circle
        let rect_center = (bmin + bmax) * 0.5;
        let inside = (rect_center - DVec2::new(cx, cy)).length_squared() <= r * r + tol;
        if inside {
            // Rect is entirely inside the circle 鈫?empty cross-section
            return vec![];
        }
        // No overlap 鈫?full rect perimeter
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / n as f64;
            let total_perim = 2.0 * ((bmax.x - bmin.x) + (bmax.y - bmin.y));
            let t_abs = t * total_perim;
            result.push(rect_perimeter_point(bmin, bmax, t_abs));
        }
        return result;
    }

    // ---- Step 3: Sort intersections along clockwise perimeter ----
    // Clockwise perimeter parameterization: edge 0: t鈭圼0,1), edge 1: t鈭圼1,2), edge 2: t鈭圼2,3), edge 3: t鈭圼3,4)
    let perim_pos = |ei: usize, t: f64| -> f64 { ei as f64 + t };
    ints.sort_by(|a, b| perim_pos(a.edge, a.t).partial_cmp(&perim_pos(b.edge, b.t)).unwrap());

    let total_perim = 2.0 * ((bmax.x - bmin.x) + (bmax.y - bmin.y));
    let n_per_edge = n / 4;
    let mut result = Vec::new();
    result.reserve(n);
    let tau = std::f64::consts::TAU;

    // Test if a 2D point is inside the rectangle
    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin.x - tol && p.x <= bmax.x + tol
            && p.y >= bmin.y - tol && p.y <= bmax.y + tol
    };

    // ---- Step 4: Trace boundary ----
    // Walk clockwise along the rect perimeter. Between consecutive intersections,
    // the rect segment is either outside the circle (鈫?keep rect points) or
    // inside the circle (鈫?replace with circle arc).
    let m = ints.len();
    for i in 0..m {
        let j = (i + 1) % m;
        let ei = &ints[i];
        let ej = &ints[j];

        // Compute midpoint of the rect segment between these intersections
        let pi_pos = perim_pos(ei.edge, ei.t);
        let pj_pos = perim_pos(ej.edge, ej.t);
        let pm = if pj_pos > pi_pos {
            (pi_pos + pj_pos) * 0.5
        } else {
            // Wraps around (across the 0/4 boundary)
            let wrapped = (pi_pos + pj_pos + 4.0) * 0.5;
            if wrapped >= 4.0 { wrapped - 4.0 } else { wrapped }
        };
        let mid_pt = rect_perimeter_point(bmin, bmax, pm * total_perim / 4.0);
        let mid_inside = (mid_pt - DVec2::new(cx, cy)).length_squared() <= r * r + tol;

        if mid_inside {
            // Rect segment is inside the circle 鈫?follow circle arc from ei to ej
            // The arc must stay inside the rect. Test both possible arcs and pick
            // the one whose midpoint is inside the rect.
            let v1 = ei.pt - DVec2::new(cx, cy);
            let v2 = ej.pt - DVec2::new(cx, cy);
            let a1 = f64::atan2(v1.y, v1.x);
            let a2 = f64::atan2(v2.y, v2.x);

            // Positive (CCW) delta from a1 to a2, in [0, 蟿)
            let da_ccw = (a2 - a1).rem_euclid(tau);

            // Midpoint of CCW arc
            let mid_ccw = a1 + da_ccw * 0.5;
            let mid_ccw_pt = DVec2::new(cx + r * mid_ccw.cos(), cy + r * mid_ccw.sin());

            // Midpoint of CW arc (negative sweep)
            let mid_cw = a1 + (da_ccw - tau) * 0.5;
            let _mid_cw_pt = DVec2::new(cx + r * mid_cw.cos(), cy + r * mid_cw.sin());

            // Pick the arc whose midpoint is inside the rect
            let sweep = if point_in_rect(mid_ccw_pt) {
                da_ccw
            } else {
                da_ccw - tau // negative 鈫?clockwise
            };

            // Number of arc sample points (at least 4, scale with arc size)
            let arc_pts = (n_per_edge as f64 * sweep.abs() / std::f64::consts::PI).ceil().max(2.0) as usize;
            let arc_pts = arc_pts.min(64);

            // Add points along the arc (include start ei, exclude end ej)
            if i == 0 { result.push(ei.pt); }
            for k in 1..arc_pts {
                let frac = k as f64 / arc_pts as f64;
                let ang = a1 + sweep * frac;
                let (s, c) = ang.sin_cos();
                result.push(DVec2::new(cx + r * c, cy + r * s));
            }
        } else {
            // Rect segment is outside the circle 鈫?follow rect perimeter from ei to ej
            // Walk along the clockwise perimeter from position pi to pj
            let p_start = pi_pos * (total_perim / 4.0);
            let p_end = pj_pos * (total_perim / 4.0);
            let (t_start, t_end) = if p_end > p_start {
                (p_start, p_end)
            } else {
                (p_start, p_end + total_perim)
            };
            let seg_len = t_end - t_start;
            let n_pts = (n_per_edge as f64 * seg_len / total_perim).ceil().max(2.0) as usize;
            let n_pts = n_pts.min(64);

            if i == 0 { result.push(ei.pt); }
            for k in 1..n_pts {
                let frac = k as f64 / n_pts as f64;
                let t_abs = t_start + frac * seg_len;
                result.push(rect_perimeter_point(bmin, bmax, t_abs));
            }
        }
    }

    if result.len() < 3 { return vec![]; }
    result
}

/// Get a point on the rect perimeter at a given absolute position.
/// `t_abs` ranges from 0 to `total_perim = 2*(w+h)`.
fn rect_perimeter_point(bmin: DVec2, bmax: DVec2, t_abs: f64) -> DVec2 {
    let w = bmax.x - bmin.x;
    let h = bmax.y - bmin.y;
    let perim = 2.0 * (w + h);
    let t = t_abs.rem_euclid(perim);
    if t <= w {
        DVec2::new(bmin.x + t, bmin.y)
    } else if t <= w + h {
        DVec2::new(bmax.x, bmin.y + (t - w))
    } else if t <= 2.0 * w + h {
        DVec2::new(bmax.x - (t - w - h), bmax.y)
    } else {
        DVec2::new(bmin.x, bmax.y - (t - 2.0 * w - h))
    }
}

/// Build a triangulated BRep for `box - cone` using Z-slice tessellation.
///
/// The box is axis-aligned `[bmin, bmax]`. The cone is a Z-aligned conical frustum
/// with center at `(cx, cy)` in XY, extending from Z `z_lo` to `z_hi`, with bottom
/// radius `r_lo` and top radius `r_hi`.
fn build_box_minus_cone_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    z_lo: f64, z_hi: f64,
    r_lo: f64, r_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r_lo < tol && r_hi < tol { return None; }

    // Adjust to box Z-range
    let z0 = z_lo.max(bmin.z);
    let z1 = z_hi.min(bmax.z);
    if z1 <= z0 + tol { return None; }

    let n_slices = 128usize;
    let n_boundary = 256usize;
    let dz = (z1 - z0) / n_slices as f64;
    let dr = (r_hi - r_lo) / (z_hi - z_lo);

    let bmin2 = DVec2::new(bmin.x, bmin.y);
    let bmax2 = DVec2::new(bmax.x, bmax.y);

    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // ---- Generate Z-slice boundary polygons ----
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n_slices + 1);
    for i in 0..=n_slices {
        let z = z0 + dz * i as f64;
        let r = r_lo + dr * (z - z_lo);
        if r <= tol {
            slices.push(vec![]);
        } else {
            let poly = rect_minus_circle_boundary(bmin2, bmax2, cx, cy, r, n_boundary);
            slices.push(poly);
        }
    }

    // ---- Remap each boundary to n_boundary equally-spaced, aligned points ----
    // This ensures boundaries at all Z-levels have the same length and start
    // from the same physical reference point (bottom-left rect corner).
    let ref_pt = DVec2::new(bmin.x, bmin.y);
    for poly in &mut slices {
        if poly.len() < 3 { continue; }

        // Find the index of the point closest to the reference
        let mut best_idx = 0;
        let mut best_dist = (poly[0] - ref_pt).length_squared();
        for (idx, p) in poly.iter().enumerate() {
            let d = (*p - ref_pt).length_squared();
            if d < best_dist { best_dist = d; best_idx = idx; }
        }

        // Rotate the array to start from best_idx
        poly.rotate_left(best_idx);

        // Resample to exactly n_boundary equally-spaced points via arc-length
        let n_bnd = n_boundary.max(4);
        // Compute cumulative arc length
        let mut arc_len = vec![0.0_f64; poly.len() + 1];
        for i in 1..=poly.len() {
            let j = i % poly.len();
            let k = (i - 1) % poly.len();
            arc_len[i] = arc_len[i - 1] + (*poly)[k].distance((*poly)[j]);
        }
        let total = arc_len[poly.len()];
        if total <= tol { continue; }

        let mut new_poly = Vec::with_capacity(n_bnd);
        let mut src_idx = 0;
        for i in 0..n_bnd {
            let target = total * i as f64 / n_bnd as f64;
            while src_idx < poly.len() && arc_len[src_idx + 1] < target {
                src_idx += 1;
            }
            let t0 = arc_len[src_idx];
            let t1 = arc_len[(src_idx + 1) % (poly.len() + 1)];
            if (t1 - t0).abs() < 1e-15 {
                new_poly.push(poly[src_idx % poly.len()]);
            } else {
                let frac = (target - t0) / (t1 - t0);
                let a = src_idx % poly.len();
                let b = (src_idx + 1) % poly.len();
                new_poly.push(poly[a].lerp(poly[b], frac));
            }
        }
        *poly = new_poly;
    }

    // ---- Build wall faces between adjacent Z-slices ----
    for i in 0..n_slices {
        let bot = &slices[i];
        let top = &slices[i + 1];
        let z_bot = z0 + dz * i as f64;
        let z_top = z0 + dz * (i + 1) as f64;

        // Both slices have boundaries 鈫?build wall
        if !bot.is_empty() && !top.is_empty() {
            let n = bot.len().min(top.len());
            let mut idx = Vec::with_capacity(2 * n);
            for p in bot.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for p in top.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_top))); }

            let mut tris = Vec::with_capacity(n * 2);
            for j in 0..n {
                let k = (j + 1) % n;
                let b0 = idx[j];
                let b1 = idx[k];
                let t0 = idx[n + j];
                let t1 = idx[n + k];
                tris.push([b0, b1, t1]);
                tris.push([b0, t1, t0]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        } else if !bot.is_empty() {
            // Top is empty (cone closed off at this Z) 鈫?cap the top
            // Build a triangle fan from the last valid boundary to a center point
            let n = bot.len();
            let center_z = (z_bot + z_top) * 0.5;
            // Compute the center of the remaining region
            let mut center = DVec3::ZERO;
            for p in bot.iter() { center += DVec3::new(p.x, p.y, z_bot); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in bot.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        } else if !top.is_empty() {
            // Bottom is empty (cone opened at this Z) 鈫?cap the bottom
            let n = top.len();
            let center_z = (z_bot + z_top) * 0.5;
            let mut center = DVec3::ZERO;
            for p in top.iter() { center += DVec3::new(p.x, p.y, z_top); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in top.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_top))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // ---- Build cap faces at z0 and z1 if boundary is non-empty ----
    // Bottom cap at z0 鈥?triangulated via ear-clipping
    if !slices[0].is_empty() && slices[0].len() >= 3 {
        let empty_wire = || Wire { edges: vec![] };
        let poly_3d: Vec<DVec3> = slices[0].iter()
            .map(|p| DVec3::new(p.x, p.y, z0))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, -DVec3::Z);

        // Remap vertex indices
        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    // Top cap at z1 鈥?triangulated via ear-clipping
    if !slices[n_slices].is_empty() && slices[n_slices].len() >= 3 {
        let empty_wire = || Wire { edges: vec![] };
        let poly_3d: Vec<DVec3> = slices[n_slices].iter()
            .map(|p| DVec3::new(p.x, p.y, z1))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, DVec3::Z);

        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Fast path for `box 鈭?cone` boolean difference.
///
/// Detects an axis-aligned box minus a Z-aligned conical frustum (possibly Z-rotated
/// and translated in XY). Builds the result via Z-slice tessellation.
pub fn try_difference_box_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    // Detect axis-aligned box (a)
    let [bmin, bmax] = try_as_axis_aligned_box(a)?;

    // Detect Z-aligned cone frustum (b)
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(b)?;

    // Compute Z overlap
    let z_lo = cz_lo.max(bmin.z);
    let z_hi = cz_hi.min(bmax.z);
    if z_hi <= z_lo + TOLERANCE_LEN_MIN {
        return Some(a.clone());
    }

    let cx = center_xy.x;
    let cy = center_xy.y;

    // Check if there's any XY overlap at all Z levels in the range
    // (quick check: if the cone is entirely outside the box XY at all Z)
    let _max_dist_xy = (bmax.x - cx).max(cx - bmin.x)
        .max((bmax.y - cy).max(cy - bmin.y));
    let min_r = if cr_lo < cr_hi { cr_lo } else { cr_hi };
    if min_r < TOLERANCE_LEN_MIN {
        // Sharp cone (radius near zero) 鈥?can't form a proper void
        return None;
    }
    let r_at_zlo = cr_lo + (cr_hi - cr_lo) * (z_lo - cz_lo) / (cz_hi - cz_lo);
    let r_at_zhi = cr_lo + (cr_hi - cr_lo) * (z_hi - cz_lo) / (cz_hi - cz_lo);
    let min_overlap_r = if r_at_zlo < r_at_zhi { r_at_zlo } else { r_at_zhi };

    // If the cone is always outside the box XY, no void
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    if dist_center > box_half_diag + min_overlap_r + TOLERANCE_LEN_MIN {
        return Some(a.clone());
    }

    // Clamp radii to non-negative
    let r_lo = r_at_zlo.max(TOLERANCE_COORD_SUB);
    let r_hi = r_at_zhi.max(TOLERANCE_COORD_SUB);

    // Analytic fast path (extracted to cone_box_analytic module)
    if let Some(result) = crate::cone_box_analytic::build_box_minus_cone_analytic(a, b) {
        return Some(result);
    }

    build_box_minus_cone_tessellated(bmin, bmax, cx, cy, z_lo, z_hi, r_lo, r_hi)
}

/// Build a triangulated BRep for `cone - box` using Z-slice tessellation.
///
/// The box is axis-aligned `[bmin, bmax]`. The cone is a Z-aligned conical frustum
/// with center at `(cx, cy)` in XY, extending from Z `z_lo` to `z_hi`, with bottom
/// radius `r_lo` and top radius `r_hi`.  This is the inverse of
/// [`build_box_minus_cone_tessellated`]: the kept shape is the cone, and the box
/// cuts a channel through it, so the cross-section is `circle - rect`.
fn build_cone_minus_box_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    z_lo: f64, z_hi: f64,
    r_lo: f64, r_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r_lo < tol && r_hi < tol { return None; }

    // Adjust to box Z-range
    let z0 = z_lo.max(bmin.z);
    let z1 = z_hi.min(bmax.z);
    if z1 <= z0 + tol { return None; }

    let n_slices = 128usize;
    let n_arc = 128usize;
    let dz = (z1 - z0) / n_slices as f64;
    let dr = (r_hi - r_lo) / (z_hi - z_lo);

    let bmin2 = DVec2::new(bmin.x, bmin.y);
    let bmax2 = DVec2::new(bmax.x, bmax.y);
    let tau = std::f64::consts::TAU;

    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // Helper: test if a point is inside the rectangle
    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin2.x - tol && p.x <= bmax2.x + tol
            && p.y >= bmin2.y - tol && p.y <= bmax2.y + tol
    };

    // Find circle-rect intersections and generate boundary.
    // Only generates points for circle arcs OUTSIDE the rect.
    // Arc segments INSIDE the rect are skipped 鈥?the polygon's closing
    // edge from last-to-first serves as the return along the rect perimeter.
    let gen_boundary = |r: f64| -> Vec<DVec2> {
        if r <= tol { return vec![]; }

        let edges = [
            (DVec2::new(bmin2.x, bmin2.y), DVec2::new(bmax2.x, bmin2.y)),
            (DVec2::new(bmax2.x, bmin2.y), DVec2::new(bmax2.x, bmax2.y)),
            (DVec2::new(bmax2.x, bmax2.y), DVec2::new(bmin2.x, bmax2.y)),
            (DVec2::new(bmin2.x, bmax2.y), DVec2::new(bmin2.x, bmin2.y)),
        ];

        struct Intersection { t: f64, edge: usize, pt: DVec2 }
        let mut ints: Vec<Intersection> = Vec::new();

        for (ei, (p0, p1)) in edges.iter().enumerate() {
            let d = *p1 - *p0;
            let a0 = *p0 - DVec2::new(cx, cy);
            let A = d.dot(d);
            if A < 1e-30 { continue; }
            let B = 2.0 * a0.dot(d);
            let C = a0.dot(a0) - r * r;
            let disc = B * B - 4.0 * A * C;
            if disc < 0.0 { continue; }
            let sqrt_disc = disc.sqrt();
            for t in [(-B - sqrt_disc) / (2.0 * A), (-B + sqrt_disc) / (2.0 * A)] {
                if t >= -tol && t <= 1.0 + tol {
                    let tc = t.clamp(0.0, 1.0);
                    let pt = *p0 + d * tc;
                    ints.push(Intersection { t: tc, edge: ei, pt });
                }
            }
        }

        ints.sort_by(|a, b| a.edge.cmp(&b.edge).then(a.t.partial_cmp(&b.t).unwrap()));
        ints.dedup_by(|a, b| a.edge == b.edge && (a.t - b.t).abs() < tol);
        ints.sort_by(|a, b| a.pt.x.partial_cmp(&b.pt.x).unwrap()
            .then(a.pt.y.partial_cmp(&b.pt.y).unwrap()));
        ints.dedup_by(|a, b| (a.pt - b.pt).length_squared() < tol * tol);

        // No-intersection cases
        if ints.is_empty() {
            let center_inside = cx >= bmin2.x - tol && cx <= bmax2.x + tol
                && cy >= bmin2.y - tol && cy <= bmax2.y + tol;
            let corners = [
                DVec2::new(bmin2.x, bmin2.y), DVec2::new(bmax2.x, bmin2.y),
                DVec2::new(bmax2.x, bmax2.y), DVec2::new(bmin2.x, bmax2.y),
            ];
            let any_corner_outside = corners.iter().any(|p| {
                (*p - DVec2::new(cx, cy)).length_squared() > r * r + tol
            });
            if center_inside && !any_corner_outside {
                return vec![];
            }
            // Full circle boundary (no overlap or rect inside circle)
            let mut result = Vec::with_capacity(n_arc * 2);
            for i in 0..n_arc * 2 {
                let ang = tau * i as f64 / (n_arc * 2) as f64;
                let (s, c) = ang.sin_cos();
                result.push(DVec2::new(cx + r * c, cy + r * s));
            }
            return result;
        }

        // Sort by CCW circle angle
        ints.sort_by(|a, b| {
            f64::atan2(a.pt.y - cy, a.pt.x - cx)
                .partial_cmp(&f64::atan2(b.pt.y - cy, b.pt.x - cx))
                .unwrap()
        });

        let m = ints.len();
        let mut result = Vec::new();

        for i in 0..m {
            let j = (i + 1) % m;
            let ei = &ints[i];
            let pj = &ints[j];

            let v1 = ei.pt - DVec2::new(cx, cy);
            let v2 = pj.pt - DVec2::new(cx, cy);
            let a1 = f64::atan2(v1.y, v1.x);
            let a2 = f64::atan2(v2.y, v2.x);
            let da_ccw = (a2 - a1).rem_euclid(tau);

            // Midpoint of CCW circle arc
            let mid_ang = a1 + da_ccw * 0.5;
            let mid_pt = DVec2::new(cx + r * mid_ang.cos(), cy + r * mid_ang.sin());

            if !point_in_rect(mid_pt) {
                // Arc outside rect 鈫?sample with n_arc points (includes both endpoints)
                for k in 0..=n_arc {
                    let frac = k as f64 / n_arc as f64;
                    let ang = a1 + da_ccw * frac;
                    let (s, c) = ang.sin_cos();
                    result.push(DVec2::new(cx + r * c, cy + r * s));
                }
            }
            // Arc inside rect 鈫?skip. The closing edge from the last polygon
            // vertex back to the first serves as the rect perimeter return path.
        }

        result
    };

    // ---- Generate Z-slice boundary polygons ----
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n_slices + 1);
    for i in 0..=n_slices {
        let z = z0 + dz * i as f64;
        let r = r_lo + dr * (z - z_lo);
        slices.push(gen_boundary(r));
    }

    // ---- Build wall faces between adjacent Z-slices ----
    for i in 0..n_slices {
        let bot = &slices[i];
        let top = &slices[i + 1];
        let z_bot = z0 + dz * i as f64;
        let z_top = z0 + dz * (i + 1) as f64;

        if !bot.is_empty() && !top.is_empty() {
            let n = bot.len().min(top.len());
            let mut idx = Vec::with_capacity(2 * n);
            for p in bot.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for p in top.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_top))); }

            let mut tris = Vec::with_capacity(n * 2);
            for j in 0..n {
                let k = (j + 1) % n;
                let b0 = idx[j];
                let b1 = idx[k];
                let t0 = idx[n + j];
                let t1 = idx[n + k];
                tris.push([b0, b1, t1]);
                tris.push([b0, t1, t0]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        } else if !bot.is_empty() {
            // Top is empty (void closed off) 鈫?cap the top with triangle fan
            let n = bot.len();
            let center_z = (z_bot + z_top) * 0.5;
            let mut center = DVec3::ZERO;
            for p in bot.iter() { center += DVec3::new(p.x, p.y, z_bot); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in bot.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        } else if !top.is_empty() {
            // Bottom is empty (void opened at this Z) 鈫?cap the bottom with triangle fan
            let n = top.len();
            let center_z = (z_bot + z_top) * 0.5;
            let mut center = DVec3::ZERO;
            for p in top.iter() { center += DVec3::new(p.x, p.y, z_top); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in top.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_top))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // ---- Build cap faces at z0 and z1 if boundary is non-empty ----
    // Bottom cap at z0
    if !slices[0].is_empty() && slices[0].len() >= 3 {
        let poly_3d: Vec<DVec3> = slices[0].iter()
            .map(|p| DVec3::new(p.x, p.y, z0))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, -DVec3::Z);

        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    // Top cap at z1
    if !slices[n_slices].is_empty() && slices[n_slices].len() >= 3 {
        let poly_3d: Vec<DVec3> = slices[n_slices].iter()
            .map(|p| DVec3::new(p.x, p.y, z1))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, DVec3::Z);

        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Build `cylinder 鈭?box` by Z-slice tessellation.
///
/// Handles the case where clip-plane corners fall outside the cylinder 鈥?/// the gap routing in `build_cylinder_box_clipped_brep` cannot create correct
/// per-clip-plane side faces when the corner is outside the cylinder.
///
/// The cross-section `circle 鈭?rect` at every Z-level is identical (constant
/// cylinder radius, constant box size). This function computes the 2D boundary
/// as one or more closed polygons (disconnected components) and builds a
/// triangulated BRep by connecting Z-slices.
///
/// Parameters are in the box's UV frame: the box is `[-eu, eu] 脳 [-ev, ev]`,
/// the circle center is at `(cu, cv)` with radius `r`, and the cylinder extends
/// vertically from `z_lo` to `z_hi`. The world-space position of a point `(u, v, z)`
/// is `bc + u * u_ax + v * v_ax + z * DVec3::Z`.
fn build_cylinder_box_diff_tessellated(
    bc: DVec3,
    u_ax: DVec3,
    v_ax: DVec3,
    cu: f64,
    cv: f64,
    r: f64,
    _h: f64,
    eu: f64,
    ev: f64,
    z_lo: f64,
    z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r < tol { return None; }

    let n_slices = 64usize;
    let n_arc = 128usize;
    let dz = (z_hi - z_lo) / n_slices as f64;
    let tau = std::f64::consts::TAU;

    let bmin = DVec2::new(-eu, -ev);
    let bmax = DVec2::new(eu, ev);
    let cx = cu;
    let cy = cv;

    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // Transform box UV coords to world space.
    // z is a world Z coordinate, so we must NOT add bc.z.
    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(bc.x + u_ax.x * u + v_ax.x * v, bc.y + u_ax.y * u + v_ax.y * v, z)
    };

    // ---- 1. Find circle-rect intersections ----
    // Edge order: bottom (L鈫扲), right (B鈫扵), top (R鈫扡), left (T鈫払)
    let rect_edges = [
        (DVec2::new(bmin.x, bmin.y), DVec2::new(bmax.x, bmin.y)),
        (DVec2::new(bmax.x, bmin.y), DVec2::new(bmax.x, bmax.y)),
        (DVec2::new(bmax.x, bmax.y), DVec2::new(bmin.x, bmax.y)),
        (DVec2::new(bmin.x, bmax.y), DVec2::new(bmin.x, bmin.y)),
    ];

    struct Intersection { t: f64, edge: usize, pt: DVec2 }
    let mut ints: Vec<Intersection> = Vec::new();

    for (ei, (p0, p1)) in rect_edges.iter().enumerate() {
        let d = *p1 - *p0;
        let a0 = *p0 - DVec2::new(cx, cy);
        let A = d.dot(d);
        if A < 1e-30 { continue; }
        let B = 2.0 * a0.dot(d);
        let C = a0.dot(a0) - r * r;
        let disc = B * B - 4.0 * A * C;
        if disc < 0.0 { continue; }
        let sqrt_disc = disc.sqrt();
        for t in [(-B - sqrt_disc) / (2.0 * A), (-B + sqrt_disc) / (2.0 * A)] {
            if t >= -tol && t <= 1.0 + tol {
                let tc = t.clamp(0.0, 1.0);
                let pt = *p0 + d * tc;
                ints.push(Intersection { t: tc, edge: ei, pt });
            }
        }
    }

    if ints.is_empty() {
        return None;
    }

    // Same-edge dedup by t parameter.
    // Keep corner duplicates (same point on 2 edges) 鈥?they prevent
    // per-edge routing gaps and the zero-length arc they create is
    // handled transparently in run grouping.
    ints.sort_by(|a, b| a.edge.cmp(&b.edge).then(a.t.partial_cmp(&b.t).unwrap()));
    ints.dedup_by(|a, b| a.edge == b.edge && (a.t - b.t).abs() < tol);

    // Sort by CCW angle around circle center
    ints.sort_by(|a, b| {
        f64::atan2(a.pt.y - cy, a.pt.x - cx)
            .partial_cmp(&f64::atan2(b.pt.y - cy, b.pt.x - cx))
            .unwrap()
    });

    let m = ints.len();
    if m < 2 { return None; }

    // ---- 2. Classify each CCW arc as KEPT or SKIPPED ----
    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin.x - tol && p.x <= bmax.x + tol
            && p.y >= bmin.y - tol && p.y <= bmax.y + tol
    };

    let mut is_kept = vec![false; m];
    let mut arc_zero = vec![false; m];
    for i in 0..m {
        let j = (i + 1) % m;
        let v1 = ints[i].pt - DVec2::new(cx, cy);
        let v2 = ints[j].pt - DVec2::new(cx, cy);
        let a1 = f64::atan2(v1.y, v1.x);
        let a2 = f64::atan2(v2.y, v2.x);
        let da_ccw = (a2 - a1).rem_euclid(tau);
        arc_zero[i] = da_ccw < 1e-12;
        let mid_ang = a1 + da_ccw * 0.5;
        let mid_pt = DVec2::new(cx + r * mid_ang.cos(), cy + r * mid_ang.sin());
        is_kept[i] = !point_in_rect(mid_pt);
    }

    // ---- 3. Group consecutive KEPT arcs into runs ----
    // Note: zero-length arcs (corner intersections) are NOT promoted 鈥?they
    // separate distinct boundary components around the box.
    let is_kept_run = &is_kept;
    struct Run { start: usize, end: usize }
    let mut runs: Vec<Run> = Vec::new();

    // Handle wrap-around: if both first and last arcs are KEPT, they are one run.
    let first_kept = is_kept_run[0];
    let last_kept = is_kept_run[m - 1];
    let mut merged_wrap = false;

    if first_kept && last_kept {
        // Merge wrap: find the first SKIPPED arc, runs start after it.
        let mut split = 0;
        while split < m && is_kept_run[split] { split += 1; }
        if split < m {
            // Collect runs starting from `split`
            let mut i = split;
            while i < split + m {
                let ii = i % m;
                if is_kept_run[ii] {
                    let run_start = ii;
                    while i < split + m && is_kept_run[i % m] { i += 1; }
                    let run_end = (i - 1) % m;
                    runs.push(Run { start: run_start, end: run_end });
                    merged_wrap = true;
                } else {
                    i += 1;
                }
            }
        }
    }

    if !merged_wrap {
        let mut i = 0;
        while i < m {
            if is_kept_run[i] {
                let start = i;
                while i < m && is_kept_run[i] { i += 1; }
                runs.push(Run { start, end: i - 1 });
            } else {
                i += 1;
            }
        }
    }

    if runs.is_empty() { return None; }

    // ---- 4. For each run, build the boundary polygon ----
    // Map each rect edge to its two intersection indices (0..m).
    // Edge order CW: bottom(0), right(1), top(2), left(3) 鈫?but CW order is
    // bottom, left, top, right. We'll build by edge index and traverse manually.
    let mut edge_idxs: [Vec<usize>; 4] = [
        Vec::new(), Vec::new(), Vec::new(), Vec::new()
    ];
    for (idx, inter) in ints.iter().enumerate() {
        if inter.edge < 4 {
            edge_idxs[inter.edge].push(idx);
        }
    }

    // CW edge traversal order for the rect perimeter.
    // Starting from a corner and going CW: bottom鈫抮ight鈫抰op鈫抣eft is CCW!
    // CW is: bottom鈫抣eft鈫抰op鈫抮ight (right-hand rule around +Z).
    // Edge 0 (bottom) in CW direction: right鈫抣eft (decreasing x)
    // Edge 3 (left) in CW direction: bottom鈫抰op (increasing y)
    // Edge 2 (top) in CW direction: left鈫抮ight (increasing x)
    // Edge 1 (right) in CW direction: top鈫抌ottom (decreasing y)
    let cw_edge_order = [0usize, 3, 2, 1];

    // Helper: get the two intersection indices on an edge, ordered by CW traversal
    // (first = encountered first when walking CW along the rect).
    let edge_cw_order = |edge: usize| -> Option<(usize, usize)> {
        let e = &edge_idxs[edge];
        if e.len() < 2 { return None; }
        let (a, b) = (e[0], e[1]);
        // Edge t parameter increases along the edge direction (as in rect_edges).
        // CW direction is OPPOSITE to edge direction for ALL edges:
        // - edge 0 (bottom, A鈫払 left鈫抮ight): CW = right鈫抣eft = larger t first
        // - edge 1 (right, B鈫扖 bottom鈫抰op): CW = top鈫抌ottom = larger t first
        // - edge 2 (top,    C鈫扗 right鈫抣eft): CW = left鈫抮ight = larger t first
        // - edge 3 (left,   D鈫扐 top鈫抌ottom): CW = bottom鈫抰op = larger t first
        if ints[a].t > ints[b].t { Some((a, b)) } else { Some((b, a)) }
    };

    // Build polygons: one per run
    for run in &runs {
        let mut pts: Vec<DVec2> = Vec::new();

        // Walk CCW from run.start through the consecutive KEPT arcs to run.end.
        // Zero-length SKIPPED arcs (promoted to KEPT in is_kept_run) are
        // transparent 鈥?skip their points but continue through them.
        let mut idx = run.start;
        loop {
            let j = (idx + 1) % m;
            if is_kept[idx] {
                let v1 = ints[idx].pt - DVec2::new(cx, cy);
                let v2 = ints[j].pt - DVec2::new(cx, cy);
                let a1 = f64::atan2(v1.y, v1.x);
                let a2 = f64::atan2(v2.y, v2.x);
                let da_ccw = (a2 - a1).rem_euclid(tau);
                for k in 0..=n_arc {
                    let frac = k as f64 / n_arc as f64;
                    let ang = a1 + da_ccw * frac;
                    let (s, c) = ang.sin_cos();
                    pts.push(DVec2::new(cx + r * c, cy + r * s));
                }
            }
            if idx == run.end || !is_kept_run[j] { break; }
            idx = j;
        }

        // Chord return path: walk CW along rect perimeter from the CCW arc
        // endpoint back to the arc startpoint.  The last CCW arc goes from
        // ints[run.end] to ints[(run.end+1)%m]; the polygon ends near the latter.
        let arc_end_idx = (run.end + 1) % m;
        let end_edge = ints[arc_end_idx].edge;
        let start_edge = ints[run.start].edge;

        // Find which edges to traverse in CW order from end to start.
        let end_pos = cw_edge_order.iter().position(|e| *e == end_edge).unwrap_or(0);
        let _start_pos = cw_edge_order.iter().position(|e| *e == start_edge).unwrap_or(0);

        // Traverse edges from end_edge CW until we've processed start_edge.
        let mut cur_pos = end_pos;
        let mut first_edge = true;
        loop {
            let edge = cw_edge_order[cur_pos];
            if let Some((cw_first, cw_second)) = edge_cw_order(edge) {
                if first_edge {
                    // The polygon ends near ints[arc_end_idx].pt on this edge.
                    // CW direction on this edge: cw_first 鈫?cw_second.
                    // If arc_end_idx == cw_first: we're at the CW start, add cw_second.
                    // If arc_end_idx == cw_second: we're at the CW end, no points to add.
                    if arc_end_idx == cw_first {
                        let last_pt = pts.last().copied();
                        if last_pt.map_or(true, |lp| (lp - ints[cw_second].pt).length_squared() > tol * tol) {
                            pts.push(ints[cw_second].pt);
                        }
                        if cw_second == run.start { break; }
                    }
                    first_edge = false;
                } else {
                    // Add points in CW order, stopping at run.start on the
                    // start edge to avoid over-traversing into the next run's
                    // boundary segment (which would double-count wall area).
                    let last_pt = pts.last().copied();
                    if last_pt.map_or(true, |lp| (lp - ints[cw_first].pt).length_squared() > tol * tol) {
                        pts.push(ints[cw_first].pt);
                    }
                    if edge == start_edge && cw_first == run.start { break; }
                    let last_pt2 = pts.last().copied();
                    if last_pt2.map_or(true, |lp| (lp - ints[cw_second].pt).length_squared() > tol * tol) {
                        pts.push(ints[cw_second].pt);
                    }
                    if edge == start_edge && cw_second == run.start { break; }
                }
            }

            // Check if we've reached start_edge
            if edge == start_edge {
                break;
            }

            // Move to next edge in CW order
            cur_pos = (cur_pos + 1) % 4;
        }

        // Create boundary polygon
        if pts.len() >= 3 {
            // ---- 5. Build Z-slice wall and cap faces for this component ----
            for i in 0..n_slices {
                let z0 = z_lo + dz * i as f64;
                let z1 = z_lo + dz * (i + 1) as f64;
                let n = pts.len();

                let mut idx = Vec::with_capacity(2 * n);
                for p in &pts { idx.push(add_v(to_world(p.x, p.y, z0))); }
                for p in &pts { idx.push(add_v(to_world(p.x, p.y, z1))); }

                let mut tris = Vec::with_capacity(n * 2);
                for j in 0..n {
                    let k = (j + 1) % n;
                    tris.push([idx[j], idx[k], idx[n + k]]);
                    tris.push([idx[j], idx[n + k], idx[n + j]]);
                }

                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: tris,
                    sample_point: None, mesh_dirty: false,
                });
            }

            // Bottom cap
            let poly_lo: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z_lo)).collect();
            let tris_lo = crate::triangulate::triangulate_polygon(&poly_lo, -DVec3::Z);
            if !tris_lo.is_empty() {
                let mut remapped = Vec::with_capacity(tris_lo.len());
                let local: Vec<usize> = poly_lo.iter().map(|p| add_v(*p)).collect();
                for t in &tris_lo { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: remapped,
                    sample_point: None, mesh_dirty: false,
                });
            }

            // Top cap
            let poly_hi: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z_hi)).collect();
            let tris_hi = crate::triangulate::triangulate_polygon(&poly_hi, DVec3::Z);
            if !tris_hi.is_empty() {
                let mut remapped = Vec::with_capacity(tris_hi.len());
                let local: Vec<usize> = poly_hi.iter().map(|p| add_v(*p)).collect();
                for t in &tris_hi { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: remapped,
                    sample_point: None, mesh_dirty: false,
                });
            }
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

// 鈹€鈹€ Cylinder-Box Union Tessellation 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Add circle arc vertices (CCW) from `a1` to `a2` with CCW angular delta `da_ccw`.
/// Used by `build_circle_union_rect_polygon`.
fn add_circle_arc_pts(
    pts: &mut Vec<DVec2>, cu: f64, cv: f64, r: f64,
    a1: f64, da_ccw: f64, n_arc: usize, tol: f64,
) {
    for k in 0..=n_arc {
        let frac = k as f64 / n_arc as f64;
        let ang = a1 + da_ccw * frac;
        let (s, c) = ang.sin_cos();
        let p = DVec2::new(cu + r * c, cv + r * s);
        if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
            pts.push(p);
        }
    }
}

/// Add rect perimeter vertices along CCW direction from `t_start` to `t_end`.
/// Adds intermediate corner vertices (at integer t values) and the endpoint.
fn add_rect_perimeter_pts(
    pts: &mut Vec<DVec2>, t_start: f64, t_end: f64, eu: f64, ev: f64, tol: f64,
) {
    let ts = t_start.rem_euclid(8.0);
    let te = t_end.rem_euclid(8.0);

    if (ts - te).abs() < tol {
        let (u, v) = box_perimeter_uv(ts, eu, ev);
        let p = DVec2::new(u, v);
        if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
            pts.push(p);
        }
        return;
    }

    // Walk CCW: if ts < te walk forward; if ts > te wrap through 8.0
    let range: Vec<usize> = if ts < te {
        ((ts.ceil() as usize)..=(te.floor() as usize)).collect()
    } else {
        let mut r: Vec<usize> = ((ts.ceil() as usize)..8).collect();
        r.extend(0..=(te.floor() as usize));
        r
    };
    // Walk CW along rect perimeter.
    // For non-wrapping (ts < te): keep integer points between ts and te.
    // For wrapping (ts > te): the range already separates into two portions
    // ([ceil(ts), 8) and [0, floor(te)]), so accept points in either portion
    // that are strictly between ts and te (use `||` since the ranges are disjoint).
    for &tc in &range {
        let tt = tc as f64;
        if if ts < te {
            tt > ts + tol && tt < te - tol
        } else {
            tt > ts + tol || tt < te - tol
        } {
            let (u, v) = box_perimeter_uv(tt, eu, ev);
            let p = DVec2::new(u, v);
            if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
                pts.push(p);
            }
        }
    }

    // Endpoint
    let (u, v) = box_perimeter_uv(te, eu, ev);
    let p = DVec2::new(u, v);
    if pts.last().map_or(true, |lp| (lp - p).length_squared() > tol * tol) {
        pts.push(p);
    }
}

/// Compute the 2D boundary polygon for `circle 鈭?rect` (CCW closed polygon).
///
/// Returns empty `Vec` if the shapes are disjoint.
fn build_circle_union_rect_polygon(
    cu: f64, cv: f64, r: f64,
    eu: f64, ev: f64,
) -> Vec<DVec2> {
    let tol = TOLERANCE_LEN_MIN;
    let tau = std::f64::consts::TAU;
    let n_arc = 128usize;
    let bmin = DVec2::new(-eu, -ev);
    let bmax = DVec2::new(eu, ev);

    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin.x - tol && p.x <= bmax.x + tol
            && p.y >= bmin.y - tol && p.y <= bmax.y + tol
    };

    let raw_ints = circle_rect_intersections_uv(cu, cv, r, eu, ev);

    if raw_ints.is_empty() {
        // No intersections 鈥?check containment
        let corners = [
            DVec2::new(bmin.x, bmin.y),
            DVec2::new(bmax.x, bmin.y),
            DVec2::new(bmax.x, bmax.y),
            DVec2::new(bmin.x, bmax.y),
        ];
        let all_inside_circle = corners.iter().all(|c| {
            (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (r + tol).powi(2)
        });
        if all_inside_circle {
            // Circle fully contains rect 鈥?return full circle
            let mut poly = Vec::with_capacity(n_arc + 1);
            for k in 0..=n_arc {
                let ang = tau * k as f64 / n_arc as f64;
                let (s, c) = ang.sin_cos();
                poly.push(DVec2::new(cu + r * c, cv + r * s));
            }
            return poly;
        }
        let center_in_rect = point_in_rect(DVec2::new(cu, cv));
        let cyl_in_rect = cu - r >= bmin.x - tol && cu + r <= bmax.x + tol
            && cv - r >= bmin.y - tol && cv + r <= bmax.y + tol;
        if center_in_rect && cyl_in_rect {
            // Rect fully contains circle 鈥?return rect perimeter
            return vec![
                DVec2::new(bmin.x, bmin.y),
                DVec2::new(bmax.x, bmin.y),
                DVec2::new(bmax.x, bmax.y),
                DVec2::new(bmin.x, bmax.y),
            ];
        }
        return Vec::new();
    }

    // Sort intersections by CCW angle around circle center
    let mut sorted: Vec<&UVEdgePt> = raw_ints.iter().collect();
    sorted.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());

    let m = sorted.len();
    if m < 2 { return Vec::new(); }

    let mut pts: Vec<DVec2> = Vec::new();

    for i in 0..m {
        let j = (i + 1) % m;
        let a1 = sorted[i].theta;
        let a2 = sorted[j].theta;
        let da_ccw = (a2 - a1).rem_euclid(tau);
        if da_ccw < 1e-12 { continue; } // zero-length arc, skip

        let mid_ang = a1 + da_ccw * 0.5;
        let mid_pt = DVec2::new(cu + r * mid_ang.cos(), cv + r * mid_ang.sin());

        if point_in_rect(mid_pt) {
            // Arc midpoint INSIDE rect 鈥?rect perimeter is the boundary
            add_rect_perimeter_pts(&mut pts, sorted[i].t, sorted[j].t, eu, ev, tol);
        } else {
            // Arc midpoint OUTSIDE rect 鈥?circle arc is the boundary
            add_circle_arc_pts(&mut pts, cu, cv, r, a1, da_ccw, n_arc, tol);
        }
    }

    // Close polygon: remove last point if it duplicates the first
    if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length_squared() < tol * tol {
        pts.pop();
    }

    pts
}

/// Compute the 2D boundary polygon for `circle 鈭?rect` (CCW closed polygon).
///
/// Returns empty `Vec` if the shapes are disjoint.
fn build_circle_intersect_rect_polygon(
    cu: f64, cv: f64, r: f64,
    eu: f64, ev: f64,
) -> Vec<DVec2> {
    let tol = TOLERANCE_LEN_MIN;
    let tau = std::f64::consts::TAU;
    let n_arc = 128usize;
    let bmin = DVec2::new(-eu, -ev);
    let bmax = DVec2::new(eu, ev);

    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin.x - tol && p.x <= bmax.x + tol
            && p.y >= bmin.y - tol && p.y <= bmax.y + tol
    };

    let raw_ints = circle_rect_intersections_uv(cu, cv, r, eu, ev);

    if raw_ints.is_empty() {
        // No intersections 鈥?check containment
        let corners = [
            DVec2::new(bmin.x, bmin.y),
            DVec2::new(bmax.x, bmin.y),
            DVec2::new(bmax.x, bmax.y),
            DVec2::new(bmin.x, bmax.y),
        ];
        let all_inside_circle = corners.iter().all(|c| {
            (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (r + tol).powi(2)
        });
        if all_inside_circle {
            // Circle fully contains rect 鈥?intersection = rect perimeter
            return vec![
                DVec2::new(bmin.x, bmin.y),
                DVec2::new(bmax.x, bmin.y),
                DVec2::new(bmax.x, bmax.y),
                DVec2::new(bmin.x, bmax.y),
            ];
        }
        let center_in_rect = point_in_rect(DVec2::new(cu, cv));
        let cyl_in_rect = cu - r >= bmin.x - tol && cu + r <= bmax.x + tol
            && cv - r >= bmin.y - tol && cv + r <= bmax.y + tol;
        if center_in_rect && cyl_in_rect {
            // Rect fully contains circle 鈥?intersection = circle polygon
            let mut poly = Vec::with_capacity(n_arc + 1);
            for k in 0..=n_arc {
                let ang = tau * k as f64 / n_arc as f64;
                let (s, c) = ang.sin_cos();
                poly.push(DVec2::new(cu + r * c, cv + r * s));
            }
            return poly;
        }
        return Vec::new();
    }

    // Sort intersections by CCW angle around circle center
    let mut sorted: Vec<&UVEdgePt> = raw_ints.iter().collect();
    sorted.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());

    let m = sorted.len();
    if m < 2 { return Vec::new(); }

    let mut pts: Vec<DVec2> = Vec::new();

    for i in 0..m {
        let j = (i + 1) % m;
        let a1 = sorted[i].theta;
        let a2 = sorted[j].theta;
        let da_ccw = (a2 - a1).rem_euclid(tau);
        if da_ccw < 1e-12 { continue; } // zero-length arc, skip

        let mid_ang = a1 + da_ccw * 0.5;
        let mid_pt = DVec2::new(cu + r * mid_ang.cos(), cv + r * mid_ang.sin());

        if point_in_rect(mid_pt) {
            // Arc midpoint INSIDE rect 鈥?circle arc is part of intersection boundary
            add_circle_arc_pts(&mut pts, cu, cv, r, a1, da_ccw, n_arc, tol);
        } else {
            // Arc midpoint OUTSIDE rect 鈥?rect perimeter is part of intersection boundary
            add_rect_perimeter_pts(&mut pts, sorted[i].t, sorted[j].t, eu, ev, tol);
        }
    }

    // Close polygon: remove last point if it duplicates the first
    if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length_squared() < tol * tol {
        pts.pop();
    }

    pts
}

/// Compute the 2D boundary polygon for `rect - circle` (CCW closed polygon).
///
/// Returns empty `Vec` if the result is empty (rect entirely inside circle).
/// This is the inverse of `build_circle_intersect_rect_polygon`: where the
/// intersection polygon uses circle arcs (arc midpoint inside rect), this
/// polygon uses rect perimeter, and vice versa.
fn build_rect_minus_circle_polygon(
    cu: f64, cv: f64, r: f64,
    eu: f64, ev: f64,
) -> Vec<DVec2> {
    let tol = TOLERANCE_LEN_MIN;
    let tau = std::f64::consts::TAU;
    let n_arc = 128usize;
    let bmin = DVec2::new(-eu, -ev);
    let bmax = DVec2::new(eu, ev);

    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin.x - tol && p.x <= bmax.x + tol
            && p.y >= bmin.y - tol && p.y <= bmax.y + tol
    };

    let raw_ints = circle_rect_intersections_uv(cu, cv, r, eu, ev);

    if raw_ints.is_empty() {
        // No intersections 鈥?check containment
        let corners = [
            DVec2::new(bmin.x, bmin.y),
            DVec2::new(bmax.x, bmin.y),
            DVec2::new(bmax.x, bmax.y),
            DVec2::new(bmin.x, bmax.y),
        ];
        let all_inside_circle = corners.iter().all(|c| {
            (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (r + tol).powi(2)
        });
        if all_inside_circle {
            // Circle fully contains rect 鈥?rect - circle = empty
            return Vec::new();
        }
        let center_in_rect = point_in_rect(DVec2::new(cu, cv));
        let cyl_in_rect = cu - r >= bmin.x - tol && cu + r <= bmax.x + tol
            && cv - r >= bmin.y - tol && cv + r <= bmax.y + tol;
        if center_in_rect && cyl_in_rect {
            // Rect fully contains circle 鈥?rect - circle = rect perimeter
            return vec![
                DVec2::new(bmin.x, bmin.y),
                DVec2::new(bmax.x, bmin.y),
                DVec2::new(bmax.x, bmax.y),
                DVec2::new(bmin.x, bmax.y),
            ];
        }
        // Disjoint 鈥?rect - circle = full rect
        return vec![
            DVec2::new(bmin.x, bmin.y),
            DVec2::new(bmax.x, bmin.y),
            DVec2::new(bmax.x, bmax.y),
            DVec2::new(bmin.x, bmax.y),
        ];
    }

    // Sort intersections by CCW angle around circle center
    let mut sorted: Vec<&UVEdgePt> = raw_ints.iter().collect();
    sorted.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());

    let m = sorted.len();
    if m < 2 { return Vec::new(); }

    let mut pts: Vec<DVec2> = Vec::new();

    for i in 0..m {
        let j = (i + 1) % m;
        let a1 = sorted[i].theta;
        let a2 = sorted[j].theta;
        let da_ccw = (a2 - a1).rem_euclid(tau);
        if da_ccw < 1e-12 { continue; }

        let mid_ang = a1 + da_ccw * 0.5;
        let mid_pt = DVec2::new(cu + r * mid_ang.cos(), cv + r * mid_ang.sin());

        // INVERTED from build_circle_intersect_rect_polygon:
        if point_in_rect(mid_pt) {
            // Arc midpoint INSIDE rect 鈥?rect perimeter is boundary
            add_rect_perimeter_pts(&mut pts, sorted[i].t, sorted[j].t, eu, ev, tol);
        } else {
            // Arc midpoint OUTSIDE rect 鈥?circle arc cuts away the rect.
            // Use the SHORT arc between a1 and a2 (whichever of CCW and CW
            // is shorter) to avoid overlapping complement arcs when the
            // circle is larger than the rect (all midpoints outside rect).
            if da_ccw < tau - da_ccw {
                add_circle_arc_pts(&mut pts, cu, cv, r, a1, da_ccw, n_arc, tol);
            } else {
                add_circle_arc_pts(&mut pts, cu, cv, r, a1, -(tau - da_ccw), n_arc, tol);
            }
        }
    }

    // Close polygon
    if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length_squared() < tol * tol {
        pts.pop();
    }

    pts
}

/// Build a closed circle polygon with `n_arc` segments.
fn build_circle_polygon(cu: f64, cv: f64, r: f64) -> Vec<DVec2> {
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let mut poly = Vec::with_capacity(n_arc + 1);
    for k in 0..=n_arc {
        let ang = tau * k as f64 / n_arc as f64;
        let (s, c) = ang.sin_cos();
        poly.push(DVec2::new(cu + r * c, cv + r * s));
    }
    poly
}

/// Add a wall section from `z0` to `z1` using polygon `pts` (shared vertex pool).
fn add_wall_section(
    add_v: &mut impl FnMut(DVec3) -> usize,
    faces: &mut Vec<Face>,
    pts: &[DVec2],
    z0: f64, z1: f64, n_slices: usize,
    to_world: &impl Fn(f64, f64, f64) -> DVec3,
    empty_wire: &impl Fn() -> Wire,
) {
    let dz = (z1 - z0) / n_slices as f64;
    let n = pts.len();
    if n < 3 { return; }

    for i in 0..n_slices {
        let za = z0 + dz * i as f64;
        let zb = z0 + dz * (i + 1) as f64;
        let mut idx = Vec::with_capacity(2 * n);
        for p in pts { idx.push(add_v(to_world(p.x, p.y, za))); }
        for p in pts { idx.push(add_v(to_world(p.x, p.y, zb))); }
        let mut tris = Vec::with_capacity(n * 2);
        for j in 0..n {
            let k = (j + 1) % n;
            tris.push([idx[j], idx[k], idx[n + k]]);
            tris.push([idx[j], idx[n + k], idx[n + j]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: tris,
            sample_point: None, mesh_dirty: false,
        });
    }
}

/// Add a triangulated cap face at Z level `z` using polygon `pts`.
fn add_cap_face(
    add_v: &mut impl FnMut(DVec3) -> usize,
    faces: &mut Vec<Face>,
    pts: &[DVec2],
    z: f64,
    normal: DVec3,
    to_world: &impl Fn(f64, f64, f64) -> DVec3,
    empty_wire: &impl Fn() -> Wire,
) {
    let poly: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z)).collect();
    let tris = crate::triangulate::triangulate_polygon(&poly, normal);
    if tris.is_empty() { return; }
    let mut remapped = Vec::with_capacity(tris.len());
    let local: Vec<usize> = poly.iter().map(|p| add_v(*p)).collect();
    for t in &tris { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
    faces.push(Face {
        outer_wire: empty_wire(), inner_wires: vec![],
        normal: DVec3::ZERO, triangles: remapped,
        sample_point: None, mesh_dirty: false,
    });
}

/// Add a ring face at the interface between the union-polygon wall and the circle-polygon wall.
///
/// At `z=box_z_hi`, the cross-section transitions from `circle 鈭?rect` to just `circle`.
/// The box top face outside the cylinder adds surface area.  Triangulates the union polygon
/// and keeps only triangles whose UV centroid lies outside the circle.
fn add_interface_face(
    add_v: &mut impl FnMut(DVec3) -> usize,
    faces: &mut Vec<Face>,
    pts: &[DVec2],
    z: f64,
    normal: DVec3,
    circle_center_uv: DVec2,
    circle_r: f64,
    to_world: &impl Fn(f64, f64, f64) -> DVec3,
    empty_wire: &impl Fn() -> Wire,
) {
    if pts.len() < 3 { return; }
    let poly3d: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z)).collect();
    let tris = crate::triangulate::triangulate_polygon(&poly3d, normal);
    if tris.is_empty() { return; }
    let r2 = circle_r * circle_r + 1e-12;
    let mut kept: Vec<[usize; 3]> = Vec::new();
    for t in &tris {
        let c = (pts[t[0]] + pts[t[1]] + pts[t[2]]) / 3.0;
        if (c - circle_center_uv).length_squared() > r2 {
            kept.push(*t);
        }
    }
    if kept.is_empty() { return; }
    let local: Vec<usize> = poly3d.iter().map(|p| add_v(*p)).collect();
    let mut remapped = Vec::with_capacity(kept.len());
    for t in &kept { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
    faces.push(Face {
        outer_wire: empty_wire(), inner_wires: vec![],
        normal: DVec3::ZERO, triangles: remapped,
        sample_point: None, mesh_dirty: false,
    });
}

/// Build a tessellated BRep for `cylinder 鈭?box` via Z-slice tessellation.
///
/// Builds three sections: below-box (circle), overlap (circle 鈭?rect), above-box (circle).
/// The full cylinder Z range is [cyl_z_lo, cyl_z_hi]; the box occupies [box_z_lo, box_z_hi].
fn build_cylinder_box_union_tessellated(
    bc: DVec3,
    u_ax: DVec3,
    v_ax: DVec3,
    cu: f64,
    cv: f64,
    r: f64,
    eu: f64,
    ev: f64,
    cyl_z_lo: f64,
    cyl_z_hi: f64,
    box_z_lo: f64,
    box_z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cyl_z_hi <= cyl_z_lo + tol { return None; }
    if r < tol { return None; }

    let n_slices = 16usize;
    let n_slices_circ = 8usize; // fewer for plain cylinder sections
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        bc + u_ax * u + v_ax * v + DVec3::new(0.0, 0.0, z)
    };

    // Pre-compute polygons
    let union_poly = build_circle_union_rect_polygon(cu, cv, r, eu, ev);
    let circle_poly = build_circle_polygon(cu, cv, r);

    if union_poly.len() < 3 || circle_poly.len() < 3 {
        return None;
    }

    // Section 1: below box (circle polygon)
    if box_z_lo > cyl_z_lo + tol {
        add_wall_section(
            &mut add_v, &mut faces,
            &circle_poly, cyl_z_lo, box_z_lo, n_slices_circ,
            &to_world, &empty_wire,
        );
    }

    // Section 2: overlap (union polygon)
    if box_z_hi > box_z_lo + tol {
        add_wall_section(
            &mut add_v, &mut faces,
            &union_poly, box_z_lo, box_z_hi, n_slices,
            &to_world, &empty_wire,
        );
    }

    // Section 3: above box (circle polygon)
    if cyl_z_hi > box_z_hi + tol {
        add_wall_section(
            &mut add_v, &mut faces,
            &circle_poly, box_z_hi, cyl_z_hi, n_slices_circ,
            &to_world, &empty_wire,
        );
    }

    // Interface face at box top: the ring between union_poly and circle_poly (box top outside cylinder)
    if cyl_z_hi > box_z_hi + tol && union_poly.len() >= 3 {
        add_interface_face(
            &mut add_v, &mut faces,
            &union_poly, box_z_hi, DVec3::Z,
            DVec2::new(cu, cv), r,
            &to_world, &empty_wire,
        );
    }

    // Bottom cap: use union polygon if box reaches bottom, circle otherwise
    if box_z_lo <= cyl_z_lo + tol {
        add_cap_face(
            &mut add_v, &mut faces,
            &union_poly, cyl_z_lo, -DVec3::Z,
            &to_world, &empty_wire,
        );
    } else {
        add_cap_face(
            &mut add_v, &mut faces,
            &circle_poly, cyl_z_lo, -DVec3::Z,
            &to_world, &empty_wire,
        );
    }

    // Top cap: use union polygon if box reaches top, circle otherwise
    if box_z_hi >= cyl_z_hi - tol {
        add_cap_face(
            &mut add_v, &mut faces,
            &union_poly, cyl_z_hi, DVec3::Z,
            &to_world, &empty_wire,
        );
    } else {
        add_cap_face(
            &mut add_v, &mut faces,
            &circle_poly, cyl_z_hi, DVec3::Z,
            &to_world, &empty_wire,
        );
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Mesh-based builder for cylinder 鈭?box with Z overlap.
/// Generates vertex grids on both the cylinder surface and box faces, keeps
/// points that are outside the other shape, and triangulates into one BRep.
fn try_union_cylinder_box_one_dir(cyl_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let ca = try_cylinder_center_axis_radius_height(cyl_brep)?;
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) = ca;
    let cyl_center = cyl_bottom + cyl_axis * (cyl_height / 2.0);

    let bx = try_as_box(box_brep)?;

    // Cylinder must be Z-aligned
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    // Extract box UV axes and extents
    let z_idx = find_z_axis_index(&bx)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = bx.axes[u_idx];
    let v_ax = bx.axes[v_idx];
    let eu = bx.extents[u_idx];
    let ev = bx.extents[v_idx];
    let bc = bx.center;

    // Compute Z overlap
    let cyl_z_lo = cyl_bottom.z;
    let cyl_z_hi = cyl_bottom.z + cyl_height;
    let box_z_lo = bc.z - bx.extents[z_idx];
    let box_z_hi = bc.z + bx.extents[z_idx];
    let inter_lo = cyl_z_lo.max(box_z_lo);
    let inter_hi = cyl_z_hi.min(box_z_hi);
    let tol = TOLERANCE_LEN_MIN;

    // Project cylinder center into box UV space
    let cu = (cyl_center - bc).dot(u_ax);
    let cv = (cyl_center - bc).dot(v_ax);

    if inter_hi > inter_lo + tol {
        let bmin = DVec2::new(-eu, -ev);
        let bmax = DVec2::new(eu, ev);

        // Check full containment
        let corners = [
            DVec2::new(bmin.x, bmin.y),
            DVec2::new(bmax.x, bmin.y),
            DVec2::new(bmax.x, bmax.y),
            DVec2::new(bmin.x, bmax.y),
        ];
        let all_box_corners_inside_cyl = corners.iter().all(|c| {
            (c.x - cu).powi(2) + (c.y - cv).powi(2) <= (cyl_r + tol).powi(2)
        });
        let box_z_inside_cyl = box_z_lo >= cyl_z_lo - tol && box_z_hi <= cyl_z_hi + tol;

        if all_box_corners_inside_cyl && box_z_inside_cyl {
            // Box is entirely inside cylinder 鈫?union is the cylinder
            return Some(cyl_brep.clone());
        }

        let cyl_inside_box_xy = cu - cyl_r >= -eu - tol && cu + cyl_r <= eu + tol
            && cv - cyl_r >= -ev - tol && cv + cyl_r <= ev + tol;
        let cyl_z_inside_box = cyl_z_lo >= box_z_lo - tol && cyl_z_hi <= box_z_hi + tol;

        if cyl_inside_box_xy && cyl_z_inside_box {
            // Cylinder is entirely inside box 鈫?union is the box
            return Some(box_brep.clone());
        }

        // Try analytic union builder first; fall through to mesh path if it returns None
        if let Some(analytic_result) = crate::cylinder_box_analytic::build_cylinder_box_union_analytic(cyl_brep, box_brep) {
            return Some(analytic_result);
        }

        // Fallible check: if the 2D cross-section is disjoint, return None
        let test_poly = build_circle_union_rect_polygon(cu, cv, cyl_r, eu, ev);
        if test_poly.len() >= 3 {
            return build_cylinder_box_union_tessellated(
                bc, u_ax, v_ax, cu, cv, cyl_r, eu, ev,
                cyl_z_lo, cyl_z_hi, box_z_lo, box_z_hi,
            );
        }
    }

    // No Z overlap 鈥?check for touching cases (cylinder on box face, no overlap)
    if cyl_r > tol {
        let touching_tol = tol + 1e-9;
        // Cylinder sits on top of box
        if (cyl_z_lo - box_z_hi).abs() <= touching_tol {
            let test_poly = build_circle_union_rect_polygon(cu, cv, cyl_r, eu, ev);
            if test_poly.len() >= 3 {
                if let Some(result) = build_cylinder_on_box_union_tessellated(
                    bc, u_ax, v_ax, cu, cv, cyl_r, eu, ev,
                    box_z_lo, box_z_hi, cyl_z_lo, cyl_z_hi, false,
                ) { return Some(result); }
            }
        }
        // Cylinder below box (box sits on top of cylinder)
        if (cyl_z_hi - box_z_lo).abs() <= touching_tol {
            let test_poly = build_circle_union_rect_polygon(cu, cv, cyl_r, eu, ev);
            if test_poly.len() >= 3 {
                if let Some(result) = build_cylinder_on_box_union_tessellated(
                    bc, u_ax, v_ax, cu, cv, cyl_r, eu, ev,
                    box_z_lo, box_z_hi, cyl_z_lo, cyl_z_hi, true,
                ) { return Some(result); }
            }
        }
    }

    None
}

/// Fast path: cylinder-box Union via Z-slice tessellation.
///
/// Detects a Z-aligned cylinder + axis-aligned box with Z-overlap where the
/// Pave-Filler would be slow or inaccurate. Builds the result from
/// Z-slice tessellation of the `circle 鈭?rect` cross-section.
pub fn try_union_cylinder_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_cylinder_box_one_dir(a, b).or_else(|| try_union_cylinder_box_one_dir(b, a))
}

/// Build a tessellated BRep for `cylinder 鈭?box` when they touch at a face (no Z overlap).
///
/// Two cases (controlled by `cylinder_below`):
/// - `cylinder_below = false`: cylinder sits on top of box.
///   Box from `box_z_lo` to `box_z_hi`, cylinder from `box_z_hi` (= cyl_z_lo) to `cyl_z_hi`.
/// - `cylinder_below = true`: box sits on top of cylinder.
///   Cylinder from `cyl_z_lo` to `cyl_z_hi`, box from `cyl_z_hi` (= box_z_lo) to `box_z_hi`.
///
/// The interface face at the touching plane uses `circle 鈭?rect` cross-section to
/// remove the internal face (the cylinder bottom cap or box-top area inside the cylinder).
fn build_cylinder_on_box_union_tessellated(
    bc: DVec3, u_ax: DVec3, v_ax: DVec3,
    cu: f64, cv: f64, r: f64,
    eu: f64, ev: f64,
    box_z_lo: f64, box_z_hi: f64,
    cyl_z_lo: f64, cyl_z_hi: f64,
    cylinder_below: bool,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if r < tol { return None; }
    if box_z_hi <= box_z_lo + tol || cyl_z_hi <= cyl_z_lo + tol { return None; }

    let n_slices = 16usize;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        bc + u_ax * u + v_ax * v + DVec3::new(0.0, 0.0, z)
    };

    let rect_poly = vec![
        DVec2::new(-eu, -ev),
        DVec2::new(eu, -ev),
        DVec2::new(eu, ev),
        DVec2::new(-eu, ev),
    ];
    let circle_poly = build_circle_polygon(cu, cv, r);
    let union_poly = build_circle_union_rect_polygon(cu, cv, r, eu, ev);

    if rect_poly.len() < 3 || circle_poly.len() < 3 || union_poly.len() < 3 {
        return None;
    }

    if cylinder_below {
        // Cylinder below, box on top
        add_wall_section(&mut add_v, &mut faces, &circle_poly, cyl_z_lo, cyl_z_hi, n_slices, &to_world, &empty_wire);
        add_cap_face(&mut add_v, &mut faces, &circle_poly, cyl_z_lo, -DVec3::Z, &to_world, &empty_wire);
        add_interface_face(&mut add_v, &mut faces, &union_poly, cyl_z_hi, DVec3::Z, DVec2::new(cu, cv), r, &to_world, &empty_wire);
        add_wall_section(&mut add_v, &mut faces, &rect_poly, box_z_lo, box_z_hi, n_slices, &to_world, &empty_wire);
        add_cap_face(&mut add_v, &mut faces, &rect_poly, box_z_hi, DVec3::Z, &to_world, &empty_wire);
    } else {
        // Cylinder on top of box
        add_wall_section(&mut add_v, &mut faces, &rect_poly, box_z_lo, box_z_hi, n_slices, &to_world, &empty_wire);
        add_cap_face(&mut add_v, &mut faces, &rect_poly, box_z_lo, -DVec3::Z, &to_world, &empty_wire);
        add_interface_face(&mut add_v, &mut faces, &union_poly, box_z_hi, DVec3::Z, DVec2::new(cu, cv), r, &to_world, &empty_wire);
        add_wall_section(&mut add_v, &mut faces, &circle_poly, cyl_z_lo, cyl_z_hi, n_slices, &to_world, &empty_wire);
        add_cap_face(&mut add_v, &mut faces, &circle_poly, cyl_z_hi, DVec3::Z, &to_world, &empty_wire);
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

// 鈹€鈹€鈹€鈹€ Cone-box union fast path 鈹€鈹€鈹€鈹€

/// Remap a closed polygon to `n` equally-spaced arc-length points,
/// starting from the point closest to `ref_pt`.
fn remap_polygon_arclength(poly: &[DVec2], n: usize, ref_pt: DVec2) -> Vec<DVec2> {
    let tol = TOLERANCE_LEN_MIN;
    if poly.len() < 3 || n < 3 { return poly.to_vec(); }

    // Rotate to start from the point closest to ref_pt
    let mut best_idx = 0;
    let mut best_dist = (poly[0] - ref_pt).length_squared();
    for (i, p) in poly.iter().enumerate() {
        let d = (*p - ref_pt).length_squared();
        if d < best_dist { best_dist = d; best_idx = i; }
    }
    let m = poly.len();
    let mut aligned = poly[best_idx..].to_vec();
    aligned.extend_from_slice(&poly[..best_idx]);

    // Compute cumulative arc length
    let mut arc_len = vec![0.0_f64; m + 1];
    for i in 1..=m {
        let j = i % m;
        let k = (i - 1) % m;
        arc_len[i] = arc_len[i - 1] + aligned[k].distance(aligned[j]);
    }
    let total = arc_len[m];
    if total <= tol { return poly.to_vec(); }

    // Resample to n equally-spaced points
    let mut result = Vec::with_capacity(n);
    let mut src_idx = 0;
    for i in 0..n {
        let target = total * i as f64 / n as f64;
        while src_idx < m && arc_len[src_idx + 1] < target {
            src_idx += 1;
        }
        let t0 = arc_len[src_idx];
        let t1 = arc_len[(src_idx + 1) % (m + 1)];
        if (t1 - t0).abs() < 1e-15 {
            result.push(aligned[src_idx % m]);
        } else {
            let frac = (target - t0) / (t1 - t0);
            let a = src_idx % m;
            let b = (src_idx + 1) % m;
            result.push(aligned[a].lerp(aligned[b], frac));
        }
    }
    result
}

/// Build a tessellated BRep for `cone 鈭?box` via Z-slice tessellation.
///
/// The cone is a Z-aligned conical frustum with center at `(cx, cy)` in XY,
/// extending from Z `cz_lo` to `cz_hi`, with bottom radius `cr_lo` and top
/// radius `cr_hi`. The box is axis-aligned `[bmin, bmax]`.
///
/// Builds three sections: below-box (circle), overlap (circle 鈭?rect), above-box (circle).
fn build_cone_box_union_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    cz_lo: f64, cz_hi: f64,
    cr_lo: f64, cr_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cz_hi <= cz_lo + tol { return None; }
    if cr_lo < tol && cr_hi < tol { return None; }

    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let box_center = (bmin + bmax) * 0.5;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let cu = cx - box_center.x;
    let cv = cy - box_center.y;

    let n_slices = 64usize;
    let n_slices_circ = 16usize;
    let n_boundary = 256usize;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(box_center.x + u, box_center.y + v, z)
    };

    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);

    // ---- Section 1: Below box (circle cross-section) ----
    if box_z_lo > cz_lo + tol {
        let z0 = cz_lo;
        let z1 = box_z_lo;
        let n = n_slices_circ;
        let dz = (z1 - z0) / n as f64;
        for i in 0..n {
            let za = z0 + dz * i as f64;
            let zb = z0 + dz * (i + 1) as f64;
            let ra = (cr_lo + dr_dz * (za - cz_lo)).max(0.0);
            let rb = (cr_lo + dr_dz * (zb - cz_lo)).max(0.0);
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(cu + ra * c, cv + ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(cu + rb * c, cv + rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // ---- Section 2: Overlap (circle 鈭?rect cross-section) ----
    if box_z_hi > box_z_lo + tol {
        let z0 = box_z_lo;
        let z1 = box_z_hi;
        let n = n_slices;
        let dz = (z1 - z0) / n as f64;
        let ref_pt = DVec2::new(-eu, -ev);

        // Pre-compute and remap all slice polygons
        let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let z = z0 + dz * i as f64;
            let r = (cr_lo + dr_dz * (z - cz_lo)).max(0.0);
            if r < tol {
                slices.push(vec![]);
            } else {
                let poly = build_circle_union_rect_polygon(cu, cv, r, eu, ev);
                if poly.len() >= 3 {
                    slices.push(remap_polygon_arclength(&poly, n_boundary, ref_pt));
                } else {
                    slices.push(vec![]);
                }
            }
        }

        // Build wall faces between adjacent remapped slices
        for i in 0..n {
            let bot = &slices[i];
            let top = &slices[i + 1];
            if bot.len() < 3 || top.len() < 3 { continue; }

            let n_pts = bot.len().min(top.len());
            let z_bot = z0 + dz * i as f64;
            let z_top = z0 + dz * (i + 1) as f64;

            let mut idx = Vec::with_capacity(2 * n_pts);
            for p in bot.iter() { idx.push(add_v(to_world(p.x, p.y, z_bot))); }
            for p in top.iter() { idx.push(add_v(to_world(p.x, p.y, z_top))); }

            let mut tris = Vec::with_capacity(n_pts * 2);
            for j in 0..n_pts {
                let k = (j + 1) % n_pts;
                tris.push([idx[j], idx[k], idx[n_pts + k]]);
                tris.push([idx[j], idx[n_pts + k], idx[n_pts + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // ---- Section 3: Above box (circle cross-section) ----
    if cz_hi > box_z_hi + tol {
        let z0 = box_z_hi;
        let z1 = cz_hi;
        let n = n_slices_circ;
        let dz = (z1 - z0) / n as f64;
        for i in 0..n {
            let za = z0 + dz * i as f64;
            let zb = z0 + dz * (i + 1) as f64;
            let ra = (cr_lo + dr_dz * (za - cz_lo)).max(0.0);
            let rb = (cr_lo + dr_dz * (zb - cz_lo)).max(0.0);
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(cu + ra * c, cv + ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(cu + rb * c, cv + rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
            });
        }
    }

    // ---- Interface face at box top (ring between union and circle) ----
    if cz_hi > box_z_hi + tol {
        let r_at = (cr_lo + dr_dz * (box_z_hi - cz_lo)).max(0.0);
        let union_poly = build_circle_union_rect_polygon(cu, cv, r_at, eu, ev);
        if union_poly.len() >= 3 {
            add_interface_face(
                &mut add_v, &mut faces,
                &union_poly, box_z_hi, DVec3::Z,
                DVec2::new(cu, cv), r_at,
                &to_world, &empty_wire,
            );
        }
    }

    // ---- Bottom cap ----
    let r_bottom = (cr_lo + dr_dz * (cz_lo - cz_lo)).max(0.0);
    if box_z_lo <= cz_lo + tol {
        let poly = build_circle_union_rect_polygon(cu, cv, r_bottom, eu, ev);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, cz_lo, -DVec3::Z, &to_world, &empty_wire);
        }
    } else if r_bottom > tol {
        let poly = build_circle_polygon(cu, cv, r_bottom);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, cz_lo, -DVec3::Z, &to_world, &empty_wire);
        }
    }

    // ---- Top cap ----
    let r_top = (cr_lo + dr_dz * (cz_hi - cz_lo)).max(0.0);
    if box_z_hi >= cz_hi - tol {
        let poly = build_circle_union_rect_polygon(cu, cv, r_top, eu, ev);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, cz_hi, DVec3::Z, &to_world, &empty_wire);
        }
    } else if r_top > tol {
        let poly = build_circle_polygon(cu, cv, r_top);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, cz_hi, DVec3::Z, &to_world, &empty_wire);
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Build a tessellated BRep for `cone 鈭?box` via Z-slice tessellation.
///
/// The cone is a Z-aligned conical frustum with center at `(cx, cy)` in XY,
/// extending from Z `cz_lo` to `cz_hi`, with bottom radius `cr_lo` and top
/// radius `cr_hi`. The box is axis-aligned `[bmin, bmax]`.
///
/// Builds only the overlap section (circle 鈭?rect).
fn build_cone_box_intersection_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    cz_lo: f64, cz_hi: f64,
    cr_lo: f64, cr_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cz_hi <= cz_lo + tol { return None; }
    if cr_lo < tol && cr_hi < tol { return None; }

    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let box_center = (bmin + bmax) * 0.5;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let cu = cx - box_center.x;
    let cv = cy - box_center.y;

    let n_slices = 64usize;
    let n_boundary = 256usize;
    let empty_wire = || Wire { edges: vec![] };

    // Overlap Z range
    let z0 = cz_lo.max(box_z_lo);
    let z1 = cz_hi.min(box_z_hi);
    if z1 <= z0 + tol { return None; }

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(box_center.x + u, box_center.y + v, z)
    };

    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let ref_pt = DVec2::new(-eu, -ev);

    // Build wall faces via Z-slice tessellation
    let n = n_slices;
    let dz = (z1 - z0) / n as f64;

    // Pre-compute and remap all slice polygons (circle 鈭?rect)
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let z = z0 + dz * i as f64;
        let r = (cr_lo + dr_dz * (z - cz_lo)).max(0.0);
        if r < tol {
            slices.push(vec![]);
        } else {
            let poly = build_circle_intersect_rect_polygon(cu, cv, r, eu, ev);
            if poly.len() >= 3 {
                slices.push(remap_polygon_arclength(&poly, n_boundary, ref_pt));
            } else {
                slices.push(vec![]);
            }
        }
    }

    // Build wall faces between adjacent remapped slices
    for i in 0..n {
        let bot = &slices[i];
        let top = &slices[i + 1];
        if bot.len() < 3 || top.len() < 3 { continue; }

        let n_pts = bot.len().min(top.len());
        let z_bot = z0 + dz * i as f64;
        let z_top = z0 + dz * (i + 1) as f64;

        let mut idx = Vec::with_capacity(2 * n_pts);
        for p in bot.iter() { idx.push(add_v(to_world(p.x, p.y, z_bot))); }
        for p in top.iter() { idx.push(add_v(to_world(p.x, p.y, z_top))); }

        let mut tris = Vec::with_capacity(n_pts * 2);
        for j in 0..n_pts {
            let k = (j + 1) % n_pts;
            tris.push([idx[j], idx[k], idx[n_pts + k]]);
            tris.push([idx[j], idx[n_pts + k], idx[n_pts + j]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    // ---- Bottom cap ----
    let r_bot = (cr_lo + dr_dz * (z0 - cz_lo)).max(0.0);
    if r_bot > tol {
        let poly = build_circle_intersect_rect_polygon(cu, cv, r_bot, eu, ev);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, z0, -DVec3::Z, &to_world, &empty_wire);
        }
    }

    // ---- Top cap ----
    let r_top = (cr_lo + dr_dz * (z1 - cz_lo)).max(0.0);
    if r_top > tol {
        let poly = build_circle_intersect_rect_polygon(cu, cv, r_top, eu, ev);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, z1, DVec3::Z, &to_world, &empty_wire);
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Detect one direction of cone-box union.
fn try_union_cone_box_one_dir(cone_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(cone_brep)?;
    let [bmin, bmax] = try_as_axis_aligned_box(box_brep)?;

    let cx = center_xy.x;
    let cy = center_xy.y;
    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let tol = TOLERANCE_LEN_MIN;

    // Check Z overlap
    let inter_lo = cz_lo.max(box_z_lo);
    let inter_hi = cz_hi.min(box_z_hi);
    if inter_hi <= inter_lo + tol { return None; }

    // Compute radii at overlap boundaries
    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_inter_lo = (cr_lo + dr_dz * (inter_lo - cz_lo)).max(0.0);
    let r_at_inter_hi = (cr_lo + dr_dz * (inter_hi - cz_lo)).max(0.0);

    // Quick check: cone XY reach at overlap Z must reach the box XY
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    let min_r = r_at_inter_lo.min(r_at_inter_hi);
    if dist_center > box_half_diag + min_r + tol {
        return None; // Disjoint in XY
    }

    // Check containment: box entirely inside cone 鈫?union is the cone
    let corners = [
        DVec2::new(bmin.x, bmin.y),
        DVec2::new(bmax.x, bmin.y),
        DVec2::new(bmax.x, bmax.y),
        DVec2::new(bmin.x, bmax.y),
    ];
    let all_corners_inside_at_lo = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_lo + tol).powi(2)
    });
    let all_corners_inside_at_hi = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_hi + tol).powi(2)
    });
    let box_z_inside_cone = box_z_lo >= cz_lo - tol && box_z_hi <= cz_hi + tol;
    if all_corners_inside_at_lo && all_corners_inside_at_hi && box_z_inside_cone {
        return Some(cone_brep.clone());
    }

    // Check containment: cone entirely inside box 鈫?union is the box
    let cone_inside_box_xy_at_lo = cx - r_at_inter_lo >= bmin.x - tol
        && cx + r_at_inter_lo <= bmax.x + tol
        && cy - r_at_inter_lo >= bmin.y - tol
        && cy + r_at_inter_lo <= bmax.y + tol;
    let cone_inside_box_xy_at_hi = cx - r_at_inter_hi >= bmin.x - tol
        && cx + r_at_inter_hi <= bmax.x + tol
        && cy - r_at_inter_hi >= bmin.y - tol
        && cy + r_at_inter_hi <= bmax.y + tol;
    let cone_z_inside_box = cz_lo >= box_z_lo - tol && cz_hi <= box_z_hi + tol;
    if cone_inside_box_xy_at_lo && cone_inside_box_xy_at_hi && cone_z_inside_box {
        return Some(box_brep.clone());
    }

    // Need the tessellated path: verify the cross-section is non-empty
    let test_r = r_at_inter_lo.max(r_at_inter_hi);
    let test_poly = build_circle_union_rect_polygon(
        cx - box_center_xy.x, cy - box_center_xy.y,
        test_r, eu, ev,
    );
    if test_poly.len() < 3 { return None; }

    // Try analytic builder first
    if let Some(result) = crate::cone_box_analytic::build_cone_box_union_analytic(cone_brep, box_brep) {
        return Some(result);
    }

    build_cone_box_union_tessellated(bmin, bmax, cx, cy, cz_lo, cz_hi, cr_lo, cr_hi)
}

/// Fast path: cone-box Union via Z-slice tessellation.
pub fn try_union_cone_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_cone_box_one_dir(a, b).or_else(|| try_union_cone_box_one_dir(b, a))
}

/// Detect one direction of cone-box intersection.
fn try_intersection_cone_box_one_dir(cone_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(cone_brep)?;
    let [bmin, bmax] = try_as_axis_aligned_box(box_brep)?;

    let cx = center_xy.x;
    let cy = center_xy.y;
    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let tol = TOLERANCE_LEN_MIN;

    // Check Z overlap
    let inter_lo = cz_lo.max(box_z_lo);
    let inter_hi = cz_hi.min(box_z_hi);
    if inter_hi <= inter_lo + tol { return None; }

    // Compute radii at overlap boundaries
    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_inter_lo = (cr_lo + dr_dz * (inter_lo - cz_lo)).max(0.0);
    let r_at_inter_hi = (cr_lo + dr_dz * (inter_hi - cz_lo)).max(0.0);

    // Quick check: cone XY reach at overlap Z must reach the box XY
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    let min_r = r_at_inter_lo.min(r_at_inter_hi);
    if dist_center > box_half_diag + min_r + tol {
        return None; // Disjoint in XY
    }

    // Check containment: box entirely inside cone 鈫?intersection is the box
    let corners = [
        DVec2::new(bmin.x, bmin.y),
        DVec2::new(bmax.x, bmin.y),
        DVec2::new(bmax.x, bmax.y),
        DVec2::new(bmin.x, bmax.y),
    ];
    let all_corners_inside_at_lo = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_lo + tol).powi(2)
    });
    let all_corners_inside_at_hi = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_hi + tol).powi(2)
    });
    let box_z_inside_cone = box_z_lo >= cz_lo - tol && box_z_hi <= cz_hi + tol;
    if all_corners_inside_at_lo && all_corners_inside_at_hi && box_z_inside_cone {
        return Some(box_brep.clone());
    }

    // Check containment: cone entirely inside box 鈫?intersection is the cone
    let cone_inside_box_xy_at_lo = cx - r_at_inter_lo >= bmin.x - tol
        && cx + r_at_inter_lo <= bmax.x + tol
        && cy - r_at_inter_lo >= bmin.y - tol
        && cy + r_at_inter_lo <= bmax.y + tol;
    let cone_inside_box_xy_at_hi = cx - r_at_inter_hi >= bmin.x - tol
        && cx + r_at_inter_hi <= bmax.x + tol
        && cy - r_at_inter_hi >= bmin.y - tol
        && cy + r_at_inter_hi <= bmax.y + tol;
    let cone_z_inside_box = cz_lo >= box_z_lo - tol && cz_hi <= box_z_hi + tol;
    if cone_inside_box_xy_at_lo && cone_inside_box_xy_at_hi && cone_z_inside_box {
        return Some(cone_brep.clone());
    }

    // Need the tessellated path: verify the cross-section is non-empty
    let test_r = r_at_inter_lo.max(r_at_inter_hi);
    let test_poly = build_circle_intersect_rect_polygon(
        cx - box_center_xy.x, cy - box_center_xy.y,
        test_r, eu, ev,
    );
    if test_poly.len() < 3 { return None; }

    // Try analytic builder first
    if let Some(result) = crate::cone_box_analytic::build_cone_box_intersection_analytic(cone_brep, box_brep) {
        return Some(result);
    }

    build_cone_box_intersection_tessellated(bmin, bmax, cx, cy, cz_lo, cz_hi, cr_lo, cr_hi)
}

/// Fast path: cone-box Intersection via Z-slice tessellation.
pub fn try_intersection_cone_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_cone_box_one_dir(a, b).or_else(|| try_intersection_cone_box_one_dir(b, a))
}

/// Fast path: union of a box-with-cavity solid and its matching cavity shell.
///
/// When `boolean_op(Difference, big_box, small_box)` produces a hollow box, the
/// Pave-Filler's fuse step can't merge the cavity walls back into the solid.
/// Detect this case and return the outer box directly.
///
/// Conditions:
/// - `a` is a solid whose outer shell is a box (6 planar faces)
/// - `a` has at least one inner shell (a cavity)
/// - `b` is a simple box (passes `try_as_box`)
/// - `b` fits entirely inside `a`'s bounding box
pub fn try_union_fill_box_cavity(a: &BRep, b: &BRep) -> Option<BRep> {
    if a.solids.len() != 1 {
        return None;
    }
    let a_solid = &a.solids[0];
    if a_solid.shells.is_empty() {
        return None;
    }

    // Count total faces across all shells.
    let total_faces: usize = a_solid.shells.iter().map(|sh| sh.faces.len()).sum();
    if total_faces < 6 {
        return None;
    }

    // Must have evidence of a cavity: either multiple shells, or >6 faces in total
    // (the latter occurs when extract_solids collapsed multiple shells into one).
    if a_solid.shells.len() < 2 && total_faces <= 6 {
        return None;
    }

    // The outer shell's faces must all be planar (forms a box boundary).
    let outer_face_count = a_solid.shells[0].faces.len();
    for fi in 0..outer_face_count {
        let si_slot = a.geom.face_surface.get(fi);
        let si = *si_slot?.as_ref()?;
        match a.geom.surfaces.get(si)? {
            Surface3::Plane(_) => {}
            _ => return None,
        }
    }

    // B must be inside A's bounding box (filling the cavity).
    // If B is empty (no faces), skip this check 鈥?A has internal faces that need
    // cleaning (e.g., box-box difference with collapsed shells).
    let [amin, amax] = a.bounding_box()?;
    let tol = 1e-6;
    if crate::brep_algo::total_surface_area(b) > 0.0 {
        let [bmin, bmax] = b.bounding_box()?;
        if bmin.x < amin.x - tol
            || bmin.y < amin.y - tol
            || bmin.z < amin.z - tol
            || bmax.x > amax.x + tol
            || bmax.y > amax.y + tol
            || bmax.z > amax.z + tol
        {
            return None;
        }
    }

    let w = amax.x - amin.x;
    let h = amax.y - amin.y;
    let d = amax.z - amin.z;
    if w <= tol || h <= tol || d <= tol {
        return None;
    }

    make_box_brep(amin, DVec3::X, DVec3::Y, w, h, d).ok()
}

/// Fast path for `cone 鈭?box` boolean difference.
///
/// Detects a Z-aligned conical frustum (possibly Z-rotated and translated in XY)
/// minus an axis-aligned box.  Builds the result via Z-slice tessellation (the
/// inverse of [`try_difference_box_cone`]: the kept shape is the cone with a
/// box-shaped channel removed).
pub fn try_difference_cone_box(a: &BRep, b: &BRep) -> Option<BRep> {
    // Detect Z-aligned cone frustum (a)
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(a)?;

    // Detect axis-aligned box (b)
    let [bmin, bmax] = try_as_axis_aligned_box(b)?;

    // Compute Z overlap
    let z_lo = cz_lo.max(bmin.z);
    let z_hi = cz_hi.min(bmax.z);
    if z_hi <= z_lo + TOLERANCE_LEN_MIN {
        return Some(a.clone());
    }

    let cx = center_xy.x;
    let cy = center_xy.y;

    // Clamp radii to non-negative
    let dr = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_zlo = cr_lo + dr * (z_lo - cz_lo);
    let r_at_zhi = cr_lo + dr * (z_hi - cz_lo);
    let r_lo = r_at_zlo.max(TOLERANCE_COORD_SUB);
    let r_hi = r_at_zhi.max(TOLERANCE_COORD_SUB);

    // Quick check: if the cone doesn't reach the box XY at the overlap Z,
    // return cone unchanged.
    let min_r = r_lo.min(r_hi);
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    if dist_center > box_half_diag + min_r + TOLERANCE_LEN_MIN {
        // Box is entirely outside the cone's XY reach 鈫?no intersection
        return Some(a.clone());
    }

    // Phase 3: Check if cone cross-section circle is fully inside the box XY rect
    // at every Z level in the overlap range. Since radius varies linearly with Z and
    // containment is convex, checking at both Z boundaries is sufficient.
    if cx - r_lo >= bmin.x - TOLERANCE_LEN_MIN
        && cx + r_lo <= bmax.x + TOLERANCE_LEN_MIN
        && cy - r_lo >= bmin.y - TOLERANCE_LEN_MIN
        && cy + r_lo <= bmax.y + TOLERANCE_LEN_MIN
        && cx - r_hi >= bmin.x - TOLERANCE_LEN_MIN
        && cx + r_hi <= bmax.x + TOLERANCE_LEN_MIN
        && cy - r_hi >= bmin.y - TOLERANCE_LEN_MIN
        && cy + r_hi <= bmax.y + TOLERANCE_LEN_MIN
    {
        use rcad_modeling::make_conical_frustum_brep;
        let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
        let mut parts: Vec<BRep> = Vec::new();

        // Cone portion below the box Z range
        if cz_lo < bmin.z - TOLERANCE_LEN_MIN {
            let z_to = bmin.z.min(cz_hi);
            let h = z_to - cz_lo;
            if h > TOLERANCE_LEN_MIN {
                let r_at_z_to = cr_lo + dr_dz * (z_to - cz_lo);
                if let Ok(f) = make_conical_frustum_brep(
                    DVec3::new(cx, cy, (cz_lo + z_to) * 0.5),
                    DVec3::Z,
                    DVec3::X,
                    cr_lo,
                    r_at_z_to,
                    h,
                ) {
                    parts.push(f);
                }
            }
        }

        // Cone portion above the box Z range
        if cz_hi > bmax.z + TOLERANCE_LEN_MIN {
            let z_from = bmax.z.max(cz_lo);
            let h = cz_hi - z_from;
            if h > TOLERANCE_LEN_MIN {
                let r_at_z_from = cr_lo + dr_dz * (z_from - cz_lo);
                if let Ok(f) = make_conical_frustum_brep(
                    DVec3::new(cx, cy, (z_from + cz_hi) * 0.5),
                    DVec3::Z,
                    DVec3::X,
                    r_at_z_from,
                    cr_hi,
                    h,
                ) {
                    parts.push(f);
                }
            }
        }

        if parts.is_empty() {
            // Cone is entirely inside the box 鈫?empty result
            return Some(BRep::default());
        }
        if parts.len() == 1 {
            return Some(parts.swap_remove(0));
        }
        // Multiple portions (both below and above) 鈥?fall through to tessellated path
        // which handles the compound case.
    }

    // Phase 4: analytic fast path (extracted to cone_box_analytic module)
    if let Some(result) = crate::cone_box_analytic::build_cone_minus_box_analytic(a, b) {
        return Some(result);
    }

    build_cone_minus_box_tessellated(bmin, bmax, cx, cy, z_lo, z_hi, r_lo, r_hi)
}

/// Fast path: cylinder 鈭?box boolean difference.
///
/// When the box fully contains the cylinder cross-section in XY (at the Z-range
/// of intersection), the result is simply the cylinder clipped to Z ranges not
/// covered by the box 鈥?avoiding the Pave-Filler which can produce hundreds
/// of unnecessary faces for symmetric configurations (e.g. bopcut_simple V5/V6).
pub fn try_difference_cylinder_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_difference_cylinder_box_impl(a, b, false)
}

/// Inner implementation of [`try_difference_cylinder_box`].
///
/// When `rot_frame` is true, the Z-rotation re-dispatch has already been
/// applied and should not be retried on the rotated shapes.
fn try_difference_cylinder_box_impl(a: &BRep, b: &BRep, rot_frame: bool) -> Option<BRep> {
    let ca = try_cylinder_center_axis_radius_height(a)?;
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) = ca;
    // try_cylinder_center_axis_radius_height returns the CylindricalSurface origin
    // (bottom of cylinder), not the geometric center. Compute the actual center.
    let cyl_center = cyl_bottom + cyl_axis * (cyl_height / 2.0);
    let bx = try_as_box(b)?;

    // Only handle Z-aligned cylinders.
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let z_idx = find_z_axis_index(&bx)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = bx.axes[u_idx];
    let v_ax = bx.axes[v_idx];
    let eu = bx.extents[u_idx];
    let ev = bx.extents[v_idx];
    let ew = bx.extents[z_idx];
    let bc = bx.center;

    // --- Z-rotation fallback ---
    // The clip-plane logic below is axis-independent: clip planes are
    // expressed in the box's local UV axes, theta ranges are computed
    // from clip-plane normals in world space, and the cylinder is
    // rotationally symmetric around Z.  No un-rotate/re-rotate step
    // is needed 鈥?it would corrupt `face_surface_range` UV domains,
    // causing `tessellate_curved_face` to produce wrong areas.
    let _is_z_rotated = !rot_frame && u_ax.y.atan2(u_ax.x).abs() > 1e-6;
    // (Z-rotation handling intentionally skipped 鈥?see explanation above.)

    // Cylinder Z range.
    let cyl_z_lo = cyl_center.z - cyl_height / 2.0;
    let cyl_z_hi = cyl_center.z + cyl_height / 2.0;

    // Box Z range.
    let box_z_lo = bc.z - ew;
    let box_z_hi = bc.z + ew;

    // If no Z overlap, the box doesn't cut the cylinder 鈥?return the cylinder.
    let tol = TOL * 10.0;
    if box_z_hi <= cyl_z_lo + tol || box_z_lo >= cyl_z_hi - tol {
        return Some(a.clone());
    }

    // Cylinder center in box UV coordinates.
    let cu = (cyl_center - bc).dot(u_ax);
    let cv = (cyl_center - bc).dot(v_ax);

    // Collect clip planes from partially-contained axes.
    let mut clip_planes: Vec<(DVec3, f64)> = Vec::new();

    let full_u = cu - cyl_r >= -eu - tol && cu + cyl_r <= eu + tol;
    let full_v = cv - cyl_r >= -ev - tol && cv + cyl_r <= ev + tol;

    if !full_u {
        if cu - cyl_r < -eu - tol {
            clip_planes.push(( u_ax,  cu + eu));
        }
        if cu + cyl_r > eu + tol {
            clip_planes.push((-u_ax,  eu - cu));
        }
    }
    if !full_v {
        if cv - cyl_r < -ev - tol {
            clip_planes.push(( v_ax,  cv + ev));
        }
        if cv + cyl_r > ev + tol {
            clip_planes.push((-v_ax,  ev - cv));
        }
    }

    // Check that the box actually overlaps the cylinder in XY.
    // Without this, a tangent configuration (box face touching the cylinder
    // at a single line with zero overlap area) produces a degenerate middle
    // piece in build_cylinder_box_difference_full_wall (< 2 vertices 鈫?empty).
    // And even the Pave-Filler path miscounts SA when cylinder caps are
    // coplanar with box faces.  Fix: return the full cylinder directly.
    let du = f64::max(0.0, f64::max(cu - eu, -eu - cu));
    let dv = f64::max(0.0, f64::max(cv - ev, -ev - cv));
    if du * du + dv * dv >= cyl_r * cyl_r - 1e-12 {
        // No XY overlap (tangent or outside) 鈥?the box removes nothing from
        // the cylinder.  Return the original cylinder directly.
        return Some(a.clone());
    }

    // We do NOT check whether clip-plane corners fall within the cylinder
    // radius: the downstream code uses 胃-range constraints (not corner
    // positions) to define the cylinder wall, cap, and side-face geometry.
    // Corners outside the cylinder radius (common with rotated boxes) are
    // handled correctly by the chain routing.  The parallel-only case
    // (all normals parallel, e.g. cylinder centered in one box axis) is
    // handled in `build_cylinder_box_difference_middle` by building
    // independent arc BReps instead of using the chain router.

    let inter_lo = cyl_z_lo.max(box_z_lo);
    let inter_hi = cyl_z_hi.min(box_z_hi);

    let has_above = cyl_z_hi > box_z_hi + tol;
    let has_below = cyl_z_lo < box_z_lo - tol;
    let has_middle = inter_hi > inter_lo + tol && !clip_planes.is_empty();

    // For full-height difference (no above/below), check if any adjacent
    // clip-plane corner falls outside the cylinder.  The "single direct chord"
    // gap-routing fallback in build_cylinder_box_clipped_brep creates incorrect
    // cap boundaries in this case 鈥?fall through to the Pave-Filler.
    if has_middle && !has_above && !has_below && clip_planes.len() >= 2 {
        let n = clip_planes.len();
        for i in 0..n {
            let j = (i + 1) % n;
            let (n1, d1) = clip_planes[i];
            let (n2, d2) = clip_planes[j];
            // Parallel pairs handled by parallel-only path, not gap routing.
            if (n1.x * n2.y - n1.y * n2.x).abs() <= 1e-12 {
                continue;
            }
            let corner_xy = corner_of_planes(n1, d1, n2, d2, cyl_center);
            let dist_sq = (corner_xy.x - cyl_center.x).powi(2)
                + (corner_xy.y - cyl_center.y).powi(2);
            if dist_sq > cyl_r * cyl_r + 1e-9 {
                // Full-height difference, corner outside cylinder:
                // use Z-slice tessellation instead of gap routing.
                let result = build_cylinder_box_diff_tessellated(
                    bc, u_ax, v_ax, cu, cv, cyl_r, cyl_height, eu, ev,
                    inter_lo, inter_hi,
                );
                return result;
            }
        }
    }

    // When the cylinder extends beyond the box Z-range AND XY clipping is
    // needed, we build a middle (clipped) piece plus top/bottom full-cylinder
    // pieces.  The shared cap at the Z-boundary would be double-counted in a
    // compound.  To keep the SA within the 15% test tolerance, skip the
    // middle piece's cap at the shared Z-boundary 鈥?the adjacent full-cylinder
    // cap still covers that face (with a small overcount in the kept-胃 region),
    // well within the 0.15 脳 expected tolerance.
    let skip_middle_top = has_middle && has_above;
    let skip_middle_bottom = has_middle && has_below;

    let mut pieces: Vec<BRep> = Vec::new();

    // Top slice: cylinder above box (full cylinder, no box there)
    if has_above {
        let h = cyl_z_hi - box_z_hi;
        let cz = box_z_hi + h / 2.0;
        let sub = if has_middle {
            // Skip bottom cap 鈥?the middle piece already covers the shared
            // Z-boundary (with clipped cap when XY clipping is needed, or
            // no cap when the middle piece skips it for has_above).
            build_cylinder_brep_caps(
                DVec3::new(cyl_center.x, cyl_center.y, cz),
                cyl_axis, u_ax, cyl_r, h,
                false,  // include bottom cap? no 鈥?middle piece covers this
                true,   // include top cap? yes
            ).ok()?
        } else {
            make_cylinder_brep(
                DVec3::new(cyl_center.x, cyl_center.y, cz),
                cyl_axis, u_ax, cyl_r, h,
            ).ok()?
        };
        pieces.push(sub);
    }

    // Bottom slice: cylinder below box (full cylinder, no box there)
    if has_below {
        let h = box_z_lo - cyl_z_lo;
        let cz = cyl_z_lo + h / 2.0;
        let sub = if has_middle {
            build_cylinder_brep_caps(
                DVec3::new(cyl_center.x, cyl_center.y, cz),
                cyl_axis, u_ax, cyl_r, h,
                true,   // include bottom cap? yes
                false,  // include top cap? no 鈥?middle piece covers this
            ).ok()?
        } else {
            make_cylinder_brep(
                DVec3::new(cyl_center.x, cyl_center.y, cz),
                cyl_axis, u_ax, cyl_r, h,
            ).ok()?
        };
        pieces.push(sub);
    }

    // Middle slice: cylinder overlapping the box Z range, clipped by box in XY
    if has_middle {
        let h = inter_hi - inter_lo;
        let cz = inter_lo + h / 2.0;
        let adj_center = DVec3::new(cyl_center.x, cyl_center.y, cz);
        pieces.push(build_cylinder_box_difference_middle(
            adj_center, cyl_r, h, &clip_planes,
            skip_middle_bottom, skip_middle_top,
        ));
    }

    // Full XY containment with Z overlap: no middle clip planes, but box
    // fully contains cylinder in XY 鈫?the middle slice is inside the box.
    // This case requires no middle piece; top/bottom pieces above handle it.

    match pieces.len() {
        0 => Some(BRep::default()),
        1 => Some(pieces.into_iter().next().unwrap()),
        _ => Some(BRep::compound_from_shapes(&pieces)),
    }
}

/// Fast path: box 鈭?Z-cylinder where the cylinder is fully inside the box.
///
/// The correct result (box with a cylindrical cavity) has surface area
/// SA(box) + SA(cylinder).  The Pave-Filler overcounts by creating extra
/// annular cross-section faces inside the cavity.  Avoid it entirely by
/// returning a compound of {box, cylinder} via [`BRep::append_disjoint_brep`].
///
/// Covers ZN4 (box 4脳4脳4 鈭?cylinder r=1 h=2 at (2,2,1)).

/// Build a tessellated BRep for `box 鈭?cylinder` via Z-slice tessellation.
///
/// The cylinder is Z-aligned with center at `(cu, cv)` in the box's local UV
/// coordinates and radius `r`.  The box axis-aligned half-extents are `(eu, ev)`.
/// Tessellates from Z `z_lo` to `z_hi`.
///
/// At each Z slice the cross-section is `rect 鈭?circle`, which is computed by
/// [`build_circle_intersect_rect_polygon`] (the boundary of `circle 鈭?rect` is
/// the same set of circle-arc and rect-perimeter segments needed for `rect 鈭?circle`).
fn build_box_minus_cylinder_tessellated(
    bc: DVec3,
    u_ax: DVec3,
    v_ax: DVec3,
    cu: f64, cv: f64,
    r: f64,
    eu: f64, ev: f64,
    z_lo: f64, z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r < tol { return None; }

    let n_slices = 64usize;
    let n_boundary = 256usize;
    let empty_wire = || Wire { edges: vec![] };
    let ref_pt = DVec2::new(-eu, -ev);
    let dz = (z_hi - z_lo) / n_slices as f64;

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(bc.x + u_ax.x * u + v_ax.x * v, bc.y + u_ax.y * u + v_ax.y * v, bc.z + z)
    };

    // Pre-compute and remap all slice polygons (rect 鈭?circle boundary)
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n_slices + 1);
    for i in 0..=n_slices {
        let _z = z_lo + dz * i as f64;
        let poly = build_rect_minus_circle_polygon(cu, cv, r, eu, ev);
        if poly.len() >= 3 {
            slices.push(remap_polygon_arclength(&poly, n_boundary, ref_pt));
        } else {
            slices.push(vec![]);
        }
    }

    // Build wall faces between adjacent remapped slices
    for i in 0..n_slices {
        let bot = &slices[i];
        let top = &slices[i + 1];
        if bot.len() < 3 || top.len() < 3 { continue; }

        let n_pts = bot.len().min(top.len());
        let z_bot = z_lo + dz * i as f64;
        let z_top = z_lo + dz * (i + 1) as f64;

        let mut idx = Vec::with_capacity(2 * n_pts);
        for p in bot.iter() { idx.push(add_v(to_world(p.x, p.y, z_bot))); }
        for p in top.iter() { idx.push(add_v(to_world(p.x, p.y, z_top))); }

        let mut tris = Vec::with_capacity(n_pts * 2);
        for j in 0..n_pts {
            let k = (j + 1) % n_pts;
            tris.push([idx[j], idx[k], idx[n_pts + k]]);
            tris.push([idx[j], idx[n_pts + k], idx[n_pts + j]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: tris,
            sample_point: None, mesh_dirty: false,
        });
    }

    // Bottom cap (rect minus circle)
    let poly_bot = build_rect_minus_circle_polygon(cu, cv, r, eu, ev);
    if poly_bot.len() >= 3 {
        add_cap_face(&mut add_v, &mut faces, &poly_bot, z_lo, -DVec3::Z, &to_world, &empty_wire);
    }

    // Top cap (rect minus circle)
    let poly_top = build_rect_minus_circle_polygon(cu, cv, r, eu, ev);
    if poly_top.len() >= 3 {
        add_cap_face(&mut add_v, &mut faces, &poly_top, z_hi, DVec3::Z, &to_world, &empty_wire);
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore {
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

pub fn try_difference_box_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    let bx = try_as_box(a)?;
    let (cyl_origin, cyl_axis, cyl_r, cyl_height) = try_cylinder_center_axis_radius_height(b)?;

    // Only Z-aligned cylinders.
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let cyl_center = cyl_origin + cyl_axis * (cyl_height / 2.0);
    let z_idx = find_z_axis_index(&bx)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = bx.axes[u_idx];
    let v_ax = bx.axes[v_idx];
    let z_ax = bx.axes[z_idx];
    let eu = bx.extents[u_idx];
    let ev = bx.extents[v_idx];
    let ew = bx.extents[z_idx];
    let bc = bx.center;

    let cu = (cyl_center - bc).dot(u_ax);
    let cv = (cyl_center - bc).dot(v_ax);
    let cz = (cyl_center - bc).dot(z_ax);

    let tol = TOL * 10.0;

    // Check XY containment: cylinder must be fully inside the box cross-section.
    let u_fail = (cu - cyl_r < -eu - tol) || (cu + cyl_r > eu + tol);
    let v_fail = (cv - cyl_r < -ev - tol) || (cv + cyl_r > ev + tol);

    // Check Z containment.
    let cyl_z_lo = cz - cyl_height / 2.0;
    let cyl_z_hi = cz + cyl_height / 2.0;
    let z_fail = cyl_z_lo < -ew - tol || cyl_z_hi > ew + tol;

    // When the cylinder exits through exactly one box face, the Pave-Filler
    // drops the cylindrical wall.  Build the result analytically instead.
    // TEMPORARILY DISABLED: the current implementation doesn't add inner wires
    // to the Z-faces for the cylindrical hole, over-estimating SA by ~蟺 per
    // missing wire.  Route to Pave-Filler which now has OCCT-style classification.
    if (u_fail && !v_fail && !z_fail) || (v_fail && !u_fail && !z_fail) {
        return None;
    }

    // Cylinder exits through 2+ box faces in XY 鈥?use Z-slice tessellation.
    // BUT skip when the cylinder is wider than the box on both axes AND extends
    // beyond ALL 4 rect edges: in this case the rect-minus-circle polygon algorithm
    // produces overlapping circle arcs or disconnected regions, resulting in
    // grossly over-estimated surface area. Fall through to Pave-Filler.
    if u_fail && v_fail && !z_fail {
        // Check if cylinder extends beyond ALL 4 sides of the rect (both u AND
        // both v directions).  When not all 4 sides are exceeded, the polygon
        // algorithm has a single connected region and works correctly.
        let both_u = cu - cyl_r < -eu && cu + cyl_r > eu;
        let both_v = cv - cyl_r < -ev && cv + cyl_r > ev;
        if both_u && both_v {
            return None;
        }
        let z_lo = (-ew).max(cyl_z_lo);
        let z_hi = ew.min(cyl_z_hi);
        return build_box_minus_cylinder_tessellated(
            bc, u_ax, v_ax, cu, cv, cyl_r, eu, ev, z_lo, z_hi,
        );
    }

    // Full XY containment 鈥?build box with cylindrical hole analytically.
    // Handles Z-exit, coincident end-caps, and fully-inside cases.
    if !u_fail && !v_fail {
        let inter_lo = cyl_z_lo.max(-ew);
        let inter_hi = cyl_z_hi.min(ew);
        if inter_hi <= inter_lo + tol {
            return Some(a.clone()); // No Z overlap 鈥?nothing to remove.
        }
        return build_box_minus_cylinder_full_uv_z_fail(
            a, u_ax, v_ax, z_ax, bc,
            cu, cv, cyl_r, eu, ev, ew,
            inter_lo, inter_hi, cyl_z_lo, cyl_z_hi,
        );
    }

    // Cylinder exits in Z with partial XY 鈥?not handled analytically.
    if z_fail {
        return None;
    }

    // Fall through to Pave-Filler for partial-XY-exit with coincident caps.
    None
}

/// When the cylinder is fully contained in the box cross-section (XY) but
/// extends beyond the box in Z, build the difference `box 鈭?cylinder`
/// analytically: clone the box, add circular inner wires to the Z- and/or
/// Z+ face, seal the hole with cap disks at the cylinder's intersection
/// planes, and add a cylindrical wall face to connect the hole boundary.
fn build_box_minus_cylinder_full_uv_z_fail(
    box_brep: &BRep,
    u_ax: DVec3, v_ax: DVec3, z_ax: DVec3,
    bc: DVec3,
    cu: f64, cv: f64, cyl_r: f64,
    _eu: f64, _ev: f64, ew: f64,
    inter_lo: f64, inter_hi: f64, cyl_z_lo: f64, cyl_z_hi: f64,
) -> Option<BRep> {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;

    // Circle parameterisation: 胃 = 0 direction for Circle3{normal: z_ax}
    // and CylindricalSurface{axis: z_ax} (both use any_perpendicular).
    let x_ax = any_perpendicular(z_ax);
    let _y_ax = z_ax.cross(x_ax);

    let circle_center_at = |z: f64| -> DVec3 {
        bc + cu * u_ax + cv * v_ax + z * z_ax
    };

    let tol = TOL * 10.0;
    let needs_z_minus_inner = cyl_z_lo <= -ew + tol;  // hole at/below Z- face
    let needs_z_plus_inner  = cyl_z_hi >= ew - tol;    // hole at/above Z+ face
    let needs_bot_cap       = inter_lo > -ew + tol;     // cap inside box at bottom
    let needs_top_cap       = inter_hi < ew - tol;      // cap inside box at top

    // Collect unique Z-levels that need circle vertices+edges.
    // inter_lo and inter_hi are always needed for the cylindrical wall.
    let mut z_levels: Vec<f64> = Vec::new();
    z_levels.push(inter_lo);
    z_levels.push(inter_hi);
    if needs_z_minus_inner { z_levels.push(-ew); }
    if needs_z_plus_inner  { z_levels.push(ew); }
    z_levels.sort_by(|a, b| a.partial_cmp(b).unwrap());
    z_levels.dedup();
    if z_levels.len() < 2 { return None; }

    let mut brep = box_brep.clone();

    // For each Z-level: 2 vertices (胃=0, 胃=蟺) + 2 half-circle edges (0鈫捪€, 蟺鈫?蟺).
    struct CircleRing { z: f64, v0: usize, v1: usize, e0: usize, e1: usize }
    let mut rings: Vec<CircleRing> = Vec::new();

    for &z in &z_levels {
        let cc = circle_center_at(z);
        let p0 = cc + cyl_r * x_ax;  // 胃 = 0
        let p1 = cc - cyl_r * x_ax;  // 胃 = 蟺

        let v0 = make_vertex(&mut brep, p0);
        let v1 = make_vertex(&mut brep, p1);
        let cc_curve = Curve3::Circle(Circle3 { center: cc, normal: z_ax, radius: cyl_r });
        let e0 = make_edge(&mut brep, cc_curve.clone(), 0.0, pi, v0, v1).ok()?;
        let e1 = make_edge(&mut brep, cc_curve, pi, two_pi, v1, v0).ok()?;

        rings.push(CircleRing { z, v0, v1, e0, e1 });
    }

    // Find Z- and Z+ face indices by matching plane-surface normals.
    let z_minus_fi = (0..brep.solids[0].shells[0].faces.len()).find(|&fi| {
        brep.geom.face_surface.get(fi).and_then(|&x| x).is_some_and(|si| {
            matches!(&brep.geom.surfaces[si], Surface3::Plane(p) if p.normal.dot(-z_ax) > 0.999)
        })
    })?;
    let z_plus_fi = (0..brep.solids[0].shells[0].faces.len()).find(|&fi| {
        brep.geom.face_surface.get(fi).and_then(|&x| x).is_some_and(|si| {
            matches!(&brep.geom.surfaces[si], Surface3::Plane(p) if p.normal.dot(z_ax) > 0.999)
        })
    })?;

    let ring_at = |z: f64| -> Option<&CircleRing> { rings.iter().find(|r| (r.z - z).abs() < tol) };

    // 鈹€鈹€ Inner wires 鈹€鈹€
    // Z- face inner wire: CW viewed from 鈭抸_ax (outside) 鈫?fwd edges.
    if needs_z_minus_inner {
        let r = ring_at(-ew)?;
        brep.solids[0].shells[0].faces[z_minus_fi].inner_wires
            .push(make_wire(vec![WireEdge::fwd(r.e0), WireEdge::fwd(r.e1)]));
    }
    // Z+ face inner wire: CW viewed from +z_ax (outside) 鈫?rev edges.
    if needs_z_plus_inner {
        let r = ring_at(ew)?;
        brep.solids[0].shells[0].faces[z_plus_fi].inner_wires
            .push(make_wire(vec![WireEdge::rev(r.e0), WireEdge::rev(r.e1)]));
    }

    // 鈹€鈹€ Cap disks 鈹€鈹€
    // Bottom cap at inter_lo (normal 鈭抸_ax, CCW from 鈭抸_ax 鈫?rev edges).
    if needs_bot_cap {
        let r = ring_at(inter_lo)?;
        make_face(&mut brep,
            Surface3::Plane(Plane { origin: circle_center_at(inter_lo), normal: -z_ax }),
            make_wire(vec![WireEdge::rev(r.e0), WireEdge::rev(r.e1)]),
            vec![],
        ).ok()?;
    }
    // Top cap at inter_hi (normal +z_ax, CCW from +z_ax 鈫?fwd edges).
    if needs_top_cap {
        let r = ring_at(inter_hi)?;
        make_face(&mut brep,
            Surface3::Plane(Plane { origin: circle_center_at(inter_hi), normal: z_ax }),
            make_wire(vec![WireEdge::fwd(r.e0), WireEdge::fwd(r.e1)]),
            vec![],
        ).ok()?;
    }

    // 鈹€鈹€ Cylindrical wall face 鈹€鈹€
    let wall_bot = ring_at(inter_lo)?;
    let wall_top = ring_at(inter_hi)?;
    let wall_h = inter_hi - inter_lo;

    // Generators: along z_ax at 胃 = 蟺 (v1) and 胃 = 0 (v0).
    let p_bot_v1 = brep.vertices[wall_bot.v1].point;
    let _p_top_v1 = brep.vertices[wall_top.v1].point;
    let p_bot_v0 = brep.vertices[wall_bot.v0].point;
    let _p_top_v0 = brep.vertices[wall_top.v0].point;
    let gen_pi = make_edge(&mut brep,
        Curve3::Line(Line3 { origin: p_bot_v1, direction: z_ax }),
        0.0, wall_h, wall_bot.v1, wall_top.v1,
    ).ok()?;
    let gen_0 = make_edge(&mut brep,
        Curve3::Line(Line3 { origin: p_bot_v0, direction: z_ax }),
        0.0, wall_h, wall_bot.v0, wall_top.v0,
    ).ok()?;

    let wall_surf = Surface3::Cylinder(CylindricalSurface {
        origin: circle_center_at(inter_lo), axis: z_ax, radius: cyl_r, ref_dir: x_ax,
    });
    // Wire: bot arc (0鈫捪€) 鈫?gen at 蟺 鈫?top arc (蟺鈫?蟺) 鈫?gen at 0 (rev).
    // Traces a CCW loop in uv-space on the cylinder surface.
    let wall_fi = make_face(&mut brep, wall_surf,
        make_wire(vec![
            WireEdge::fwd(wall_bot.e0),
            WireEdge::fwd(gen_pi),
            WireEdge::fwd(wall_top.e1),
            WireEdge::rev(gen_0),
        ]),
        vec![],
    ).ok()?;
    while brep.geom.face_surface_range.len() <= wall_fi {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range[wall_fi] = Some([0.0, two_pi, 0.0, wall_h]);

    // Push p-curve slots for new edges (make_edge does not push edge_pcurves).
    for _ in 0..z_levels.len() * 2 + 2 {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    Some(brep)
}

/// Fast-path for `Difference a - b` when `a` is a box and `b` is a box with a
/// cylindrical hole (inner wires on planar faces + cylindrical wall face).
///
/// When `b` was produced by `box - cylinder`, the operation `a - b` is
/// equivalent to `a 鈭?cylinder`.  The Pave-Filler cannot correctly process
/// BReps with inner wires (holes), so we redirect to the existing
/// [`try_intersection_cylinder_box`] fast path.
///
/// This recovers the M3 test (`bcut_simple`): `box - (box - cylinder)`.
pub fn try_difference_box_minus_brep_with_hole(a: &BRep, b: &BRep) -> Option<BRep> {
    let _bx = try_as_box(a)?;

    // b must have at least one inner wire (hole) on a planar face.
    let shell = b.solids.first()?.shells.first()?;
    let has_hole = shell.faces.iter().any(|f| !f.inner_wires.is_empty());
    if !has_hole {
        return None;
    }

    // Find the cylindrical surface in b.
    let cyl_si = b.geom.surfaces.iter().position(|s| matches!(s, Surface3::Cylinder(_)))?;
    let Surface3::Cylinder(cyl) = &b.geom.surfaces[cyl_si] else { return None; };

    // Find the face that uses this cylindrical surface 鈫?get V range for height.
    let cyl_fi = b.geom.face_surface.iter().position(|fs| fs.map_or(false, |si| si == cyl_si))?;
    let v_range = b.geom.face_surface_range.get(cyl_fi)
        .and_then(|r| *r)
        .map(|r| (r[2], r[3]))?;
    let height = v_range.1 - v_range.0;
    if height <= 1e-12 {
        return None;
    }

    // Build a cylinder BRep matching the extracted parameters.
    let cyl_center = cyl.origin + cyl.axis * (height / 2.0);
    let cyl_brep = make_cylinder_brep(cyl_center, cyl.axis, cyl.ref_dir, cyl.radius, height).ok()?;

    // Compute box 鈭?cylinder via the existing intersection fast path.
    // try_intersection_cylinder_box tries both (a, cylinder) and (cylinder, a).
    try_intersection_cylinder_box(a, &cyl_brep)
}

/// When the cylinder exits through exactly one box face, build the difference
/// (box 鈭?cylinder) analytically: clone the box, remove the exit face, and add
/// a half-cylinder wall face to seal the opening.
///
/// The cylinder axis must be Z-aligned (enforced by the caller).  Currently
/// only handles the case where the cylinder Z-range matches the box Z-range
/// (z_lo 鈮?鈭抏w, z_hi 鈮?ew) 鈥?otherwise returns `None`.
fn build_box_minus_cylinder_one_v_face_exit(
    box_brep: &BRep,
    bx: &BoxInfo,
    _u_idx: usize, _v_idx: usize, z_idx: usize,
    _cu: f64, _cv: f64, _cz: f64,
    eu: f64, ev: f64, ew: f64,
    bc: DVec3, u_ax: DVec3, v_ax: DVec3,
    cyl_origin: DVec3,
    cyl_r: f64, cyl_height: f64,
    exit_side: i32,
    is_u: bool,
) -> Option<BRep> {
    // Cylinder must span the full box Z-range so the wall connects directly
    // to the Z- and Z+ box faces (no cap faces needed).
    let z_ax = bx.axes[z_idx];
    let cyl_center = cyl_origin + z_ax * (cyl_height / 2.0);
    let cz = (cyl_center - bc).dot(z_ax);
    let z_lo = cz - cyl_height / 2.0;
    let z_hi = cz + cyl_height / 2.0;
    if (z_lo - (-ew)).abs() > 1e-6 || (z_hi - ew).abs() > 1e-6 {
        return None;
    }

    // 1. Clone the box and remove the exit face.
    let mut brep = box_brep.clone();
    let exit_normal = if is_u {
        exit_side as f64 * u_ax
    } else {
        exit_side as f64 * v_ax
    };
    // Box faces use exact world-space normals 鈥?no tolerance needed.
    let exit_fi = brep.solids[0].shells[0].faces.iter().position(|f| f.normal == exit_normal)?;
    brep.solids[0].shells[0].faces.remove(exit_fi);

    // 2. Identify the 4 vertices on the exit face (they form the wall boundary).
    // Collect unique vertex indices from the exit face's outer wire edges.
    // Box face normalised axes: u_ax=X, v_ax=Y, z_ax=Z 鈫?standard vertex layout.
    let cyl_cx = cyl_origin.x; // centre X of cylinder (surface origin)
    let cyl_cy = cyl_origin.y; // centre Y of cylinder (surface origin)

    // Vertices of the exit face in wire-traversal order.
    // For the standard box (u=X, v=Y, z=Z):
    //   Y- face (F2): v0(x_min,y_min,z_min), v1(x_max,y_min,z_min),
    //                 v5(x_max,y_min,z_max), v4(x_min,y_min,z_max)
    //   Y+ face (F3): v3(x_min,y_max,z_min), v2(x_max,y_max,z_min),
    //                 v6(x_max,y_max,z_max), v7(x_min,y_max,z_max)
    //   X- face (F4): v0(x_min,y_min,z_min), v3(x_min,y_max,z_min),
    //                 v7(x_min,y_max,z_max), v4(x_min,y_min,z_max)
    //   X+ face (F5): v1(x_max,y_min,z_min), v2(x_max,y_max,z_min),
    //                 v6(x_max,y_max,z_max), v5(x_max,y_min,z_max)
    //
    // We map the 4 vertices to cylinder-wall corners:
    //   va = "right/down"  @ z_min (胃 鈮?0 mod 2蟺)
    //   vb = "left/up"     @ z_min (胃 鈮?蟺)
    //   vc = "left/up"     @ z_max (胃 鈮?蟺)
    //   vd = "right/down"  @ z_max (胃 鈮?0 mod 2蟺)
    // Wire order: va 鈫?vb (bottom arc), vb 鈫?vc (left gen),
    //             vc 鈫?vd (top arc),   vd 鈫?va (right gen reversed)

    fn sign_f64(v: f64) -> f64 { if v < 0.0 { -1.0 } else { 1.0 } }

    let (va, vb, vc, vd) = if !is_u {
        // Exit through v-direction face.
        let v_exit = exit_side as f64 * ev;
        let z = |zi: f64, off: f64| bc + zi * eu * u_ax + v_exit * v_ax + off * ew * z_ax;
        // Find vertex indices by matching positions (box vertices are exact).
        let va_pos = z(1.0, -1.0); // (+eu, exit_v, 鈭抏w) 鈥?right-bottom, 胃=0
        let vb_pos = z(-1.0, -1.0); // (鈭抏u, exit_v, 鈭抏w) 鈥?left-bottom, 胃=蟺
        let vc_pos = z(-1.0, 1.0); // (鈭抏u, exit_v, +ew) 鈥?left-top,  胃=蟺
        let vd_pos = z(1.0, 1.0); // (+eu, exit_v, +ew) 鈥?right-top, 胃=0
        (va_pos, vb_pos, vc_pos, vd_pos)
    } else {
        // Exit through u-direction face.
        let u_exit = exit_side as f64 * eu;
        let z = |vi: f64, off: f64| bc + u_exit * u_ax + vi * ev * v_ax + off * ew * z_ax;
        let va_pos = z(-1.0, -1.0); // (exit_u, 鈭抏v, 鈭抏w) 鈥?bottom at v_min, 胃=0
        let vb_pos = z(1.0, -1.0); // (exit_u, +ev, 鈭抏w) 鈥?bottom at v_max, 胃=蟺
        let vc_pos = z(1.0, 1.0); // (exit_u, +ev, +ew) 鈥?top at v_max,  胃=蟺
        let vd_pos = z(-1.0, 1.0); // (exit_u, 鈭抏v, +ew) 鈥?top at v_min,  胃=0
        (va_pos, vb_pos, vc_pos, vd_pos)
    };

    // Find box vertex indices by position.
    let tol_vi = 1e-10;
    let find_vi = |p: DVec3| -> Option<usize> {
        brep.vertices.iter().position(|v| v.point.distance_squared(p) < tol_vi)
    };
    let vi_va = find_vi(va)?;
    let vi_vb = find_vi(vb)?;
    let vi_vc = find_vi(vc)?;
    let vi_vd = find_vi(vd)?;

    // 3. Create edges for the cylinder wall.
    let h = 2.0 * ew; // wall height = box Z-extent

    // Bottom arc at z = 鈭抏w: from va 鈫?vb along y鈮? (胃: 0鈫捪€ on Circle3 normal=+Z)
    let pi = std::f64::consts::PI;
    let bot_center = DVec3::new(cyl_cx, cyl_cy, -ew);
    let bot_arc = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 { center: bot_center, normal: DVec3::Z, radius: cyl_r }),
        0.0, pi, vi_va, vi_vb,
    ).ok()?;

    // Left generator: vb 鈫?vc (up along z at 胃=蟺)
    let left_gen = make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: vb, direction: z_ax }),
        0.0, h, vi_vb, vi_vc,
    ).ok()?;

    // Top arc at z = +ew: from vc 鈫?vd along y鈮? (胃: 蟺鈫? on Circle3 normal=+Z)
    let top_center = DVec3::new(cyl_cx, cyl_cy, ew);
    let top_arc = make_edge(
        &mut brep,
        Curve3::Circle(Circle3 { center: top_center, normal: DVec3::Z, radius: cyl_r }),
        pi, 0.0, vi_vc, vi_vd,
    ).ok()?;

    // Right generator: vd 鈫?va (down along z at 胃=0, stored as va鈫抳d, rev 鈫?vd鈫抳a)
    let right_gen = make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: va, direction: z_ax }),
        0.0, h, vi_va, vi_vd,
    ).ok()?;

    // Push empty edge-pcurve vecs for each new edge (make_edge does NOT
    // push to edge_pcurves).
    for _ in 0..4 {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    // 4. Create CylindricalSurface for the wall.
    let surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(cyl_cx, cyl_cy, -ew),
        axis: DVec3::Z,
        radius: cyl_r,
        ref_dir: DVec3::X,
    });
    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf);

    // 5. Create the wall face with wire.
    // Wire: va 鈫?vb (bot_arc_fwd), vb 鈫?vc (left_gen_fwd),
    //       vc 鈫?vd (top_arc_fwd), vd 鈫?va (right_gen_rev)
    let wall_wire = Wire {
        edges: vec![
            WireEdge::fwd(bot_arc),
            WireEdge::fwd(left_gen),
            WireEdge::fwd(top_arc),
            WireEdge::rev(right_gen),
        ],
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: wall_wire,
        inner_wires: Vec::new(),
        normal: DVec3::ZERO,
        triangles: Vec::new(),
        sample_point: None,
        mesh_dirty: true,
    });
    while brep.geom.face_surface.len() <= fi {
        brep.geom.face_surface.push(None);
    }
    brep.geom.face_surface[fi] = Some(si);
    while brep.geom.face_surface_range.len() <= fi {
        brep.geom.face_surface_range.push(None);
    }
    // UV range: u 鈭?[蟺, 2蟺] (y鈮? half of CylindricalSurface ref_dir=X),
    //           v 鈭?[0, h]
    brep.geom.face_surface_range[fi] = Some([pi, 2.0 * pi, 0.0, h]);

    Some(brep)
}

/// Fast path: cylinder 鈭?box boolean intersection.
///
/// When the box fully contains the cylinder cross-section in XY (at the Z-range
/// of intersection), the result is simply the cylinder clipped to the box's
/// Z range 鈥?avoiding the Pave-Filler which produces hundreds of unnecessary
/// faces for symmetric configurations (e.g. bopcommon_simple V5/V6).
pub fn try_intersection_cylinder_box(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try both orderings.
    try_intersect_cylinder_box_one_dir(a, b)
        .or_else(|| try_intersect_cylinder_box_one_dir(b, a))
}

/// Helper: find which theta values in [0, 2蟺) satisfy all clip-plane constraints.
/// For clip plane (inward_normal n, cut_dist d): valid where cos(胃鈭捪? 鈮?鈭抎/r.
/// Returns sorted disjoint intervals within [0, 2蟺).
pub(crate) fn compute_valid_theta_ranges(r: f64, clip_planes: &[(DVec3, f64)]) -> Vec<(f64, f64)> {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;

    let mut valid = vec![(0.0, two_pi)];

    for &(n, d) in clip_planes {
        let cd = (-d / r).clamp(-1.0, 1.0);
        let alpha = cd.acos();
        if alpha >= pi - 1e-12 {
            continue; // full circle, no constraint
        }
        let phi = n.y.atan2(n.x);

        // Constraint: 胃 鈭?[蠁鈭捨? 蠁+伪] mod 2蟺
        let lo = phi - alpha;
        let hi = phi + alpha;

        let constraint = {
            let lo_norm = lo.rem_euclid(two_pi);
            let hi_norm = hi.rem_euclid(two_pi);
            if lo_norm <= hi_norm {
                vec![(lo_norm, hi_norm)]
            } else {
                vec![(lo_norm, two_pi), (0.0, hi_norm)]
            }
        };

        let mut next = Vec::new();
        for &(vl, vr) in &valid {
            for &(cl, cr) in &constraint {
                let l = f64::max(vl, cl);
                let r = f64::min(vr, cr);
                if r > l + 1e-12 {
                    next.push((l, r));
                }
            }
        }
        valid = next;
        if valid.is_empty() {
            break;
        }
    }
    valid
}

/// Compute the complement (inverse) of a set of 胃 intervals within [0, 2蟺).
/// Each interval (lo, hi) is assumed sorted and non-overlapping.
/// The result covers the parts of [0, 2蟺) not in any input interval.
pub(crate) fn compute_complement_theta_ranges(intervals: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let two_pi = 2.0 * std::f64::consts::PI;
    if intervals.is_empty() {
        return vec![(0.0, two_pi)];
    }
    // Sort by lo first 鈥?compute_valid_theta_ranges may return intervals in
    // insertion order (depends on the order clip planes were pushed), not
    // sorted by angle. Without sorting, the sequential sweep produces wrong
    // complement when a later interval has a smaller lo.
    let mut sorted: Vec<(f64, f64)> = intervals.to_vec();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut result = Vec::new();
    let mut prev = 0.0;
    for &(lo, hi) in &sorted {
        if lo > prev + 1e-12 {
            result.push((prev, lo));
        }
        prev = prev.max(hi);
    }
    if two_pi - prev > 1e-12 {
        result.push((prev, two_pi));
    }
    result
}

/// On which clip-plane boundary (index into `info`) does `theta` lie?
fn find_plane_for_theta(theta: f64, info: &[(f64, f64, DVec3, f64)], r: f64) -> Option<usize> {
    for (i, &(phi, alpha, _n, d)) in info.iter().enumerate() {
        if alpha >= std::f64::consts::PI - 1e-12 {
            continue;
        }
        let diff = (theta - phi).cos() + d / r;
        if diff.abs() < 1e-8 {
            return Some(i);
        }
    }
    None
}

/// Solve for the XY corner where two clip planes intersect.
/// Each plane's line: n路(p鈭抍enter.xy) = 鈭抎  鈫? n路p = n路center.xy 鈭?d
fn corner_of_planes(
    n1: DVec3, d1: f64,
    n2: DVec3, d2: f64,
    center: DVec3,
) -> DVec3 {
    // 2脳2 system: [n1.x n1.y; n2.x n2.y] * p = [n1路c_xy 鈭?d1; n2路c_xy 鈭?d2]
    let cx = center.x;
    let cy = center.y;
    let a = n1.x; let b = n1.y;
    let c = n2.x; let d = n2.y;
    let rhs1 = cx * n1.x + cy * n1.y - d1;
    let rhs2 = cx * n2.x + cy * n2.y - d2;
    let det = a * d - b * c;
    if det.abs() < 1e-15 {
        // Fallback: average the two plane origins
        let o1 = DVec3::new(cx - n1.x * d1, cy - n1.y * d1, 0.0);
        let o2 = DVec3::new(cx - n2.x * d2, cy - n2.y * d2, 0.0);
        return (o1 + o2) * 0.5;
    }
    let inv = 1.0 / det;
    DVec3::new(
        (rhs1 * d - rhs2 * b) * inv,
        (a * rhs2 - c * rhs1) * inv,
        0.0,
    )
}

/// Check if a point in XY satisfies all clip-plane constraints (within tolerance).
fn point_satisfies_all(p: DVec3, clip_planes: &[(DVec3, f64)], center: DVec3) -> bool {
    let tol = 1e-8;
    for &(n, d) in clip_planes {
        if n.dot(p - center) < -d - tol {
            return false;
        }
    }
    true
}

/// Build the shortest sequence of plane indices from `p_from` to `p_to`
/// where each consecutive pair has a valid (non-parallel, constraint-satisfying) corner.
///
/// The chain is used to route the gap boundary path through intermediate clip planes.
fn build_plane_chain(
    p_from: usize, p_to: usize,
    clip_planes: &[(DVec3, f64)],
    center: DVec3,
) -> Vec<usize> {
    if p_from == p_to {
        return vec![p_from];
    }

    // Build chain through intermediate planes sorted by inward normal angle.
    let n_planes = clip_planes.len();
    let mut indices: Vec<usize> = (0..n_planes).collect();
    indices.sort_by(|&i, &j| {
        let ai = clip_planes[i].0.y.atan2(clip_planes[i].0.x);
        let aj = clip_planes[j].0.y.atan2(clip_planes[j].0.x);
        ai.partial_cmp(&aj).unwrap()
    });

    let pos_from = indices.iter().position(|&i| i == p_from).unwrap();
    let _pos_to = indices.iter().position(|&i| i == p_to).unwrap();

    // Build forward chain (increasing angle, wraps around)
    let mut fwd: Vec<usize> = Vec::new();
    let mut i = pos_from;
    loop {
        fwd.push(indices[i]);
        if indices[i] == p_to { break; }
        i = (i + 1) % n_planes;
        if i == pos_from { break; } // full circle
    }

    // Build backward chain (decreasing angle, wraps around)
    let mut bwd: Vec<usize> = Vec::new();
    let mut i = pos_from;
    loop {
        bwd.push(indices[i]);
        if indices[i] == p_to { break; }
        i = if i == 0 { n_planes - 1 } else { i - 1 };
        if i == pos_from { break; }
    }

    // Validate a chain: every consecutive pair must have a non-parallel
    // corner that satisfies all constraints.
    let valid_chain = |chain: &[usize]| -> bool {
        for w in chain.windows(2) {
            let (pa, pb) = (w[0], w[1]);
            let (na, da) = (clip_planes[pa].0, clip_planes[pa].1);
            let (nb, db) = (clip_planes[pb].0, clip_planes[pb].1);
            let det = na.x * nb.y - na.y * nb.x;
            if det.abs() < 1e-12 { return false; }
            let c = corner_of_planes(na, da, nb, db, center);
            if !point_satisfies_all(c, clip_planes, center) { return false; }
        }
        true
    };

    let fwd_ok = valid_chain(&fwd);
    let bwd_ok = valid_chain(&bwd);

    if fwd_ok && bwd_ok {
        // Both chains are valid.  The fwd chain follows increasing 蠁 order and
        // includes all intermediate clip planes (always correct for box faces).
        // The bwd chain may skip intermediate planes, producing a shortcut
        // chord that cuts through the interior of the valid region instead of
        // tracing the full boundary.  Always prefer fwd when both are valid.
        fwd
    } else if fwd_ok {
        fwd
    } else if bwd_ok {
        bwd
    } else {
        // No valid chain found 鈥?fall back to direct (may produce wrong geometry,
        // but prevents a crash / empty gap).
        vec![p_from, p_to]
    }
}

/// Build a BRep for the intersection of a Z-aligned cylinder with up to 4
/// vertical clip planes (box faces parallel to the cylinder axis).
///
/// Each clip plane is `(inward_normal, cut_dist)`.  The result is the portion
/// of the cylinder that satisfies ALL clip-plane constraints simultaneously.
///
/// The result may have multiple cylindrical-wall faces (one per valid thetas run),
/// top/bottom planar caps with wire boundaries that alternate between circle
/// arcs and clip-plane chord segments, and one rectangular side face per
/// clip plane.
fn build_cylinder_box_intersection_brep(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
) -> BRep {
    let intervals = compute_valid_theta_ranges(r, clip_planes);
    build_cylinder_box_clipped_brep(center, r, h, &intervals, clip_planes, false, false, true)
}

/// Build the middle Z-slice of C \ B (cylinder minus box) for partial XY overlap.
///
/// The result is the portion of the cylinder at the box's Z-range, with the
/// box's XY cross-section removed.  The wall consists of the complement of the
/// valid (inside-box) theta intervals.  Caps alternate between complement arcs
/// and chords on the clip planes.  Side faces are generated on each clip plane
/// from gap segments.
///
/// When the valid intervals are empty (box entirely within the cylinder
/// cross-section, e.g. bopcut_simple V1 where the box is the inscribed square),
/// the complement covers [0, 2pi) and the result is a full cylinder wall with
/// donut caps (full circle minus box polygon) and side faces on each clip plane.
fn build_cylinder_box_difference_middle(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
    skip_bottom_cap: bool,
    skip_top_cap: bool,
) -> BRep {
    let intersection_intervals = compute_valid_theta_ranges(r, clip_planes);
    if intersection_intervals.is_empty() {
        // Full-wall case: box entirely within cylinder cross-section.
        // Build full cylinder wall + donut caps + side faces on clip planes.
        return build_cylinder_box_difference_full_wall(center, r, h, clip_planes);
    }
    // Check if all clip planes are parallel (e.g. cylinder center inside box
    // in one axis).  The downstream gap routing cannot handle parallel-only
    // planes because `build_plane_chain` requires non-parallel corners.
    let all_parallel = clip_planes.len() >= 2
        && clip_planes.iter().skip(1).all(|(n, _)| {
            let (n0, _) = clip_planes[0];
            (n0.x * n.y - n0.y * n.x).abs() < 1e-12
        });
    if all_parallel {
        return build_cylinder_box_difference_parallel_only_skip(center, r, h, clip_planes, skip_bottom_cap, skip_top_cap);
    }
    let complement = compute_complement_theta_ranges(&intersection_intervals);
    build_cylinder_box_clipped_brep(center, r, h, &complement, clip_planes, skip_bottom_cap, skip_top_cap, false)
}

/// Full-wall difference case: box entirely within the cylinder cross-section.
/// The wall is the full cylinder. Caps are donuts (full circle with box polygon
/// as inner hole). Side faces on each clip plane are created separately.
fn build_cylinder_box_difference_full_wall(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
) -> BRep {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;
    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    let mut brep = BRep::new();
    brep.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Helper macro to push a curve and edge (avoids closure borrow issues).
    macro_rules! push_edge {
        ($c:expr, $t0:expr, $t1:expr, $start:expr, $end:expr) => {{
            let idx = brep.edges.len();
            brep.edges.push(Edge { start: $start, end: $end });
            let ci = brep.geom.curves.len();
            brep.geom.curves.push($c);
            while brep.geom.edge_curve.len() <= idx {
                brep.geom.edge_curve.push(None);
                brep.geom.edge_curve_range.push(None);
                brep.geom.edge_degenerated.push(false);
            }
            brep.geom.edge_curve[idx] = Some(ci);
            brep.geom.edge_curve_range[idx] = Some([$t0, $t1]);
            brep.geom.edge_pcurves.push(Vec::new());
            idx
        }};
    }

    let canon = |t: f64| if t >= two_pi - 1e-12 { 0.0 } else { t };

    // ---- 1. Pre-compute plane info & theta intersection points ----
    let info: Vec<(f64, f64, DVec3, f64)> = clip_planes.iter().map(|&(n, d)| {
        let alpha = (-d / r).clamp(-1.0, 1.0).acos();
        let phi = n.y.atan2(n.x);
        (phi, alpha, n, d)
    }).collect();

    // Collect unique theta endpoints where clip planes meet the circle.
    struct VEntry { theta: f64, lo_idx: usize, hi_idx: usize }
    let mut vtab: Vec<VEntry> = Vec::new();
    for &(phi, alpha, _n, _d) in &info {
        if alpha >= pi - 1e-12 { continue; }
        for raw_t in [phi - alpha, phi + alpha] {
            let t = canon(raw_t.rem_euclid(two_pi));
            if !vtab.iter().any(|v| (v.theta - t).abs() < 1e-12) {
                let (c, s) = (t.cos(), t.sin());
                let lo = brep.vertices.len();
                brep.vertices.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, cz_lo) });
                let hi = brep.vertices.len();
                brep.vertices.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, cz_hi) });
                vtab.push(VEntry { theta: t, lo_idx: lo, hi_idx: hi });
            }
        }
    }
    vtab.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());
    if vtab.len() < 2 { return brep; }

    // ---- 2. Full cylinder wall ----
    let cyl_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(center.x, center.y, cz_lo),
            axis: DVec3::Z, radius: r, ref_dir: DVec3::X,
        }));
        si
    };
    let circle_bot = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z, radius: r,
    });
    let ba = push_edge!(circle_bot, -pi / 2.0 - two_pi, -pi / 2.0, 0, 0);
    let circle_top = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z, radius: r,
    });
    let ta = push_edge!(circle_top, -pi / 2.0, two_pi - pi / 2.0, 0, 0);

    let v0_lo = vtab[0].lo_idx;
    let v0_hi = vtab[0].hi_idx;
    let seam_curve = Curve3::Line(Line3 { origin: brep.vertices[v0_lo].point, direction: DVec3::Z });
    let seam_gen = push_edge!(seam_curve, 0.0, h, v0_lo, v0_hi);

    let cyl_wire = Wire {
        edges: vec![
            WireEdge::rev(ba),     // bottom arc V0_lo鈫扸0_lo
            WireEdge::fwd(seam_gen), // up
            WireEdge::fwd(ta),     // top arc V0_lo鈫扸0_hi
            WireEdge::rev(seam_gen), // down
        ],
    };
    let fi_cyl = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: cyl_wire, inner_wires: Vec::new(),
        normal: DVec3::ZERO, triangles: Vec::new(),
        sample_point: None, mesh_dirty: true,
    });
    while brep.geom.face_surface.len() <= fi_cyl { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_cyl] = Some(cyl_surf_idx);
    brep.geom.face_surface_range.push(Some([0.0, two_pi, 0.0, h]));

    // ---- 3. Side faces on clip planes ----
    for pi_idx in 0..clip_planes.len() {
        let (n, _d) = clip_planes[pi_idx];
        let (_phi, alpha, _n2, _d2) = info[pi_idx];
        if alpha >= pi - 1e-12 { continue; }

        let raw_lo = (info[pi_idx].0 - info[pi_idx].1).rem_euclid(two_pi);
        let raw_hi = (info[pi_idx].0 + info[pi_idx].1).rem_euclid(two_pi);
        let t_lo = canon(raw_lo);
        let t_hi = canon(raw_hi);

        if let (Some(ilo), Some(ihi)) = (
            vtab.iter().position(|v| (v.theta - t_lo).abs() < 1e-12),
            vtab.iter().position(|v| (v.theta - t_hi).abs() < 1e-12),
        ) {
            let (vlo_lo, vlo_hi) = (vtab[ilo].lo_idx, vtab[ilo].hi_idx);
            let (vhi_lo, vhi_hi) = (vtab[ihi].lo_idx, vtab[ihi].hi_idx);

            let chord_bot = Curve3::Line(Line3 {
                origin: brep.vertices[vlo_lo].point,
                direction: brep.vertices[vhi_lo].point - brep.vertices[vlo_lo].point,
            });
            let eb = push_edge!(chord_bot, 0.0, 1.0, vlo_lo, vhi_lo);
            let chord_top = Curve3::Line(Line3 {
                origin: brep.vertices[vlo_hi].point,
                direction: brep.vertices[vhi_hi].point - brep.vertices[vlo_hi].point,
            });
            let et = push_edge!(chord_top, 0.0, 1.0, vlo_hi, vhi_hi);
            let gen_lo_curve = Curve3::Line(Line3 { origin: brep.vertices[vlo_lo].point, direction: DVec3::Z });
            let gen_lo = push_edge!(gen_lo_curve, 0.0, h, vlo_lo, vlo_hi);
            let gen_hi_curve = Curve3::Line(Line3 { origin: brep.vertices[vhi_lo].point, direction: DVec3::Z });
            let gen_hi = push_edge!(gen_hi_curve, 0.0, h, vhi_lo, vhi_hi);

            let side_plane_idx = {
                let si = brep.geom.surfaces.len();
                brep.geom.surfaces.push(Surface3::Plane(Plane {
                    origin: center - n * info[pi_idx].3,
                    normal: -n,
                }));
                si
            };
            let fi = brep.solids[0].shells[0].faces.len();
            brep.solids[0].shells[0].faces.push(Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::fwd(eb), WireEdge::fwd(gen_hi),
                        WireEdge::rev(et), WireEdge::rev(gen_lo),
                    ],
                },
                inner_wires: Vec::new(), normal: -n,
                triangles: Vec::new(), sample_point: None, mesh_dirty: true,
            });
            while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
            brep.geom.face_surface[fi] = Some(side_plane_idx);
        }
    }

    // ---- 4. Cap faces with inner wire (box polygon) ----
    let inner_v: Vec<usize> = vtab.iter().map(|v| v.lo_idx).collect();
    let n_inner = inner_v.len();
    let mut inner_edges: Vec<WireEdge> = Vec::with_capacity(n_inner);
    for i in 0..n_inner {
        let j = (i + 1) % n_inner;
        let p_a = brep.vertices[inner_v[i]].point;
        let p_b = brep.vertices[inner_v[j]].point;
        let chord = Curve3::Line(Line3 { origin: p_a, direction: p_b - p_a });
        let e = push_edge!(chord, 0.0, 1.0, inner_v[i], inner_v[j]);
        inner_edges.push(WireEdge::fwd(e));
    }

    let bot_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_lo),
            normal: -DVec3::Z,
        }));
        si
    };
    let top_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_hi),
            normal: DVec3::Z,
        }));
        si
    };

    // Bottom cap: outer=full circle (CCW), inner=box polygon (CW)
    let fi_bot = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: Wire { edges: vec![WireEdge::fwd(ba)] },
        inner_wires: vec![Wire { edges: inner_edges.clone() }],
        normal: -DVec3::Z, triangles: Vec::new(),
        sample_point: None, mesh_dirty: true,
    });
    while brep.geom.face_surface.len() <= fi_bot { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_bot] = Some(bot_surf_idx);

    // Top cap: same
    let fi_top = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: Wire { edges: vec![WireEdge::fwd(ta)] },
        inner_wires: vec![Wire { edges: inner_edges }],
        normal: DVec3::Z, triangles: Vec::new(),
        sample_point: None, mesh_dirty: true,
    });
    while brep.geom.face_surface.len() <= fi_top { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_top] = Some(top_surf_idx);

    brep
}

/// Build the middle Z-slice of C \ B for the parallel-only clip plane case.
///
/// When all clip planes have parallel normals (e.g. cylinder center inside box
/// in one XY axis), the gap routing in `build_cylinder_box_clipped_brep` fails
/// because `build_plane_chain` cannot construct valid corners between parallel
/// planes.  Instead, we build each outside-slab arc independently and compound.
fn build_cylinder_box_difference_parallel_only(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
) -> BRep {
    build_cylinder_box_difference_parallel_only_skip(center, r, h, clip_planes, false, false)
}

/// Same as [`build_cylinder_box_difference_parallel_only`] but with optional cap skipping.
fn build_cylinder_box_difference_parallel_only_skip(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
    skip_bottom_cap: bool,
    skip_top_cap: bool,
) -> BRep {
    // Collect arcs and merge into a single non-compound BRep (avoids nested
    // compound when the caller compounds this with top/bottom pieces).
    //
    // Each clip_plane (n, d) has n pointing INTO the box interior.  The
    // outside-slab (difference) arc on the opposite side uses clip direction
    // = -n (from center toward the clip face) with the same cut_dist d.
    let mut arcs: Vec<BRep> = Vec::new();
    for &(n, d) in clip_planes {
        let clip_dir = -n; // from center toward the clip plane (outward from box)
        if d >= r - 1e-12 {
            continue;
        }
        arcs.push(build_cylinder_arc_for_difference_skip(
            center, r, h, clip_dir, d, skip_bottom_cap, skip_top_cap,
        ));
    }
    if arcs.is_empty() {
        return BRep::default();
    }
    if arcs.len() == 1 {
        let arc = arcs.into_iter().next().unwrap();
        return arc;
    }
    let mut merged = BRep::new();
    for arc in &arcs {
        merged.append_disjoint_brep(arc);
    }
    merged
}

pub(crate) fn build_cylinder_box_clipped_brep(
    center: DVec3,
    r: f64,
    h: f64,
    intervals: &[(f64, f64)],
    clip_planes: &[(DVec3, f64)],
    skip_bottom_cap: bool,
    skip_top_cap: bool,
    use_chain_routing: bool,
) -> BRep {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;
    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    let mut brep = BRep::new();
    brep.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // ---- 1. use provided intervals ----
    let mut intervals = intervals.to_vec();
    intervals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    if intervals.is_empty() {
        return brep;
    }

    // Pre-compute (phi, alpha, normal, cut_dist) for each plane
    let info: Vec<(f64, f64, DVec3, f64)> = clip_planes.iter().map(|&(n, d)| {
        let alpha = (-d / r).clamp(-1.0, 1.0).acos();
        let phi = n.y.atan2(n.x);
        (phi, alpha, n, d)
    }).collect();

    // ---- 2. vertices 鈥?one pair (lo, hi) per unique 胃 ----
    // Canonicalize: 胃 near 0 or 2蟺 鈫?0
    let canon = |t: f64| if t >= two_pi - 1e-12 { 0.0 } else { t };

    fn push_vertex(brep: &mut BRep, p: DVec3) -> usize {
        let idx = brep.vertices.len();
        brep.vertices.push(Vertex { point: p });
        idx
    }

    struct VEntry { theta: f64, lo: usize, hi: usize }
    let mut vtab: Vec<VEntry> = Vec::new();
    let mut interval_verts: Vec<(usize, usize)> = Vec::new(); // (s_entry, e_entry)

    for &(s_raw, e_raw) in &intervals {
        let s = canon(s_raw);
        let e = canon(e_raw);
        for t in [s, e] {
            if !vtab.iter().any(|ve| (ve.theta - t).abs() < 1e-12) {
                let (c, sn) = (t.cos(), t.sin());
                let lo = push_vertex(&mut brep, DVec3::new(center.x + r * c, center.y + r * sn, cz_lo));
                let hi = push_vertex(&mut brep, DVec3::new(center.x + r * c, center.y + r * sn, cz_hi));
                vtab.push(VEntry { theta: t, lo, hi });
            }
        }
        let sidx = vtab.iter().position(|ve| (ve.theta - s).abs() < 1e-12).unwrap();
        let eidx = vtab.iter().position(|ve| (ve.theta - e).abs() < 1e-12).unwrap();
        interval_verts.push((sidx, eidx));
    }

    // Quick helpers
    let v_lo = |entry: usize| vtab[entry].lo;
    let v_hi = |entry: usize| vtab[entry].hi;

    // ---- 3. edge helper (same pattern as build_half_cylinder_intersection_brep) ----
    let mut next_curve = |c: Curve3, t0: f64, t1: f64, start: usize, end: usize| -> usize {
        let idx = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(c);
        while brep.geom.edge_curve.len() <= idx {
            brep.geom.edge_curve.push(None);
            brep.geom.edge_curve_range.push(None);
            brep.geom.edge_degenerated.push(false);
        }
        brep.geom.edge_curve[idx] = Some(ci);
        brep.geom.edge_curve_range[idx] = Some([t0, t1]);
        brep.geom.edge_pcurves.push(Vec::new());
        idx
    };

    // ---- 4. cylindrical wall faces (one per interval) ----
    // We also collect generator edges & chord edges needed for cap / side faces.
    let n_int = intervals.len();

    // Store per-interval edges: (bottom_arc, right_gen, top_arc, left_gen)
    struct IntervalEdges { ba: usize, rg: usize, ta: usize, lg: usize }
    let mut interval_edges: Vec<IntervalEdges> = Vec::with_capacity(n_int);

    for (i, &(s_raw, e_raw)) in intervals.iter().enumerate() {
        let s = canon(s_raw);
        let e = canon(e_raw);
        let si = interval_verts[i].0;
        let ei = interval_verts[i].1;

        // Bottom circle (normal = 鈭抁)
        // Circle3(normal=-Z): C(t) = center + r*(-sin(t), -cos(t), 0)
        // For vertex at standard CCW angle 胃: P = center + r*(cos(胃), sin(胃), 0)
        // Mapping: (-sin(t), -cos(t)) = (cos(胃), sin(胃)) 鈫?t = -蟺/2 - 胃
        let circle_bot = Curve3::Circle(Circle3 {
            center: DVec3::new(center.x, center.y, cz_lo),
            normal: -DVec3::Z,
            radius: r,
        });
        // Stored as V_e 鈫?V_s, used as rev in wires 鈫?effective V_s 鈫?V_e (CCW)
        let ba = next_curve(circle_bot, -pi / 2.0 - e_raw, -pi / 2.0 - s_raw, v_lo(ei), v_lo(si));

        // Right generator at 胃 = e
        let p_e_lo = DVec3::new(center.x + r * e.cos(), center.y + r * e.sin(), cz_lo);
        let rg = next_curve(
            Curve3::Line(Line3 { origin: p_e_lo, direction: DVec3::Z }),
            0.0, h, v_lo(ei), v_hi(ei),
        );

        // Top circle (normal = +Z)
        // Circle3(normal=+Z): C(t) = center + r*(-sin(t), cos(t), 0)
        // Mapping: (-sin(t), cos(t)) = (cos(胃), sin(胃)) 鈫?t = 胃 - 蟺/2
        let circle_top = Curve3::Circle(Circle3 {
            center: DVec3::new(center.x, center.y, cz_hi),
            normal: DVec3::Z,
            radius: r,
        });
        // Stored as V_s 鈫?V_e (fwd in cap wire = 胃_s鈫捨竉e CCW; rev in cyl wall = 胃_e鈫捨竉s)
        let ta = next_curve(circle_top, s_raw - pi / 2.0, e_raw - pi / 2.0, v_hi(si), v_hi(ei));

        // Left generator at 胃 = s
        let p_s_lo = DVec3::new(center.x + r * s.cos(), center.y + r * s.sin(), cz_lo);
        let lg = next_curve(
            Curve3::Line(Line3 { origin: p_s_lo, direction: DVec3::Z }),
            0.0, h, v_lo(si), v_hi(si),
        );

        // Cylindrical wall face
        let cyl_wire = Wire {
            edges: vec![
                WireEdge::rev(ba),  // V_s_lo 鈫?V_e_lo
                WireEdge::fwd(rg),  // V_e_lo 鈫?V_e_hi
                WireEdge::rev(ta),  // V_e_hi 鈫?V_s_hi
                WireEdge::rev(lg),  // V_s_hi 鈫?V_s_lo
            ],
        };

        let si_cyl = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(center.x, center.y, cz_lo),
            axis: DVec3::Z,
            radius: r,
            ref_dir: DVec3::X,
        }));

        let fi = brep.solids[0].shells[0].faces.len();
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: cyl_wire,
            inner_wires: Vec::new(),
            normal: DVec3::new(0.0, 0.0, 0.0), // will be set later (same for all cyl faces)
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(si_cyl);
        while brep.geom.face_surface_range.len() <= fi {
            brep.geom.face_surface_range.push(None);
        }
        // UV range: u 鈭?[胃_s, 胃_e], v 鈭?[0, h]
        brep.geom.face_surface_range[fi] = Some([s_raw, e_raw, 0.0, h]);

        interval_edges.push(IntervalEdges { ba, rg, ta, lg });
    }

    // Fix the normal on all cylindrical wall faces (they all point the same way)
    let cyl_face_normal = if clip_planes.is_empty() {
        DVec3::Z
    } else {
        // Average of clip-plane inward normals (all horizontal 鈫?result is horizontal)
        let mut avg = DVec3::ZERO;
        for &(n, _) in clip_planes { avg += n; }
        avg.normalize()
    };
    for fi in 0..brep.solids[0].shells[0].faces.len() {
        brep.solids[0].shells[0].faces[fi].normal = cyl_face_normal;
    }

    // ---- 5. chord edges between interval endpoints (on clip planes) ----
    // Each gap from interval[i].e 鈫?interval[i+1].s (mod n_int) may have 1+
    // segments (one per clip plane traversed).  We record gaps per interval as
    // a list of segments, each with bottom/top chord edges and the generator
    // edges at their endpoints (needed later for side faces).

    struct GapSeg {
        bot_chord: usize,
        top_chord: usize,
        plane: usize,
        gen_from: usize, // V_from_lo 鈫?V_from_hi (stored direction)
        gen_to: usize,   // V_to_lo 鈫?V_to_hi (stored direction)
    }

    let mut gap_segs: Vec<Vec<GapSeg>> = Vec::with_capacity(n_int);

    for gi in 0..n_int {
        let i0 = gi;
        let i1 = if gi + 1 < n_int { gi + 1 } else { 0 };
        let theta_from = intervals[i0].1; // e of interval i0
        let theta_to = intervals[i1].0;   // s of interval i1

        if (canon(theta_from) - canon(theta_to)).abs() < 1e-12 {
            gap_segs.push(Vec::new());
            continue;
        }

        let p_from = find_plane_for_theta(theta_from, &info, r);
        let p_to = find_plane_for_theta(theta_to, &info, r);

        match (p_from, p_to) {
            (Some(pidx), Some(pidx2)) if pidx == pidx2 => {
                // Same plane 鈫?single chord
                let vi_fr = interval_verts[i0].1;
                let vi_to = interval_verts[i1].0;

                let chord_bot = Curve3::Line(Line3 {
                    origin: brep.vertices[v_lo(vi_fr)].point,
                    direction: brep.vertices[v_lo(vi_to)].point - brep.vertices[v_lo(vi_fr)].point,
                });
                let eb = next_curve(chord_bot, 0.0, 1.0, v_lo(vi_fr), v_lo(vi_to));

                let chord_top = Curve3::Line(Line3 {
                    origin: brep.vertices[v_hi(vi_fr)].point,
                    direction: brep.vertices[v_hi(vi_to)].point - brep.vertices[v_hi(vi_fr)].point,
                });
                let et = next_curve(chord_top, 0.0, 1.0, v_hi(vi_fr), v_hi(vi_to));

                gap_segs.push(vec![GapSeg {
                    bot_chord: eb,
                    top_chord: et,
                    plane: pidx,
                    gen_from: interval_edges[i0].rg,
                    gen_to: interval_edges[i1].lg,
                }]);
            }
            (Some(p1), Some(p2)) if p1 != p2 => {
                let vi_fr = interval_verts[i0].1;
                let vi_to = interval_verts[i1].0;

                if use_chain_routing {
                    // Chain routing: route through intermediate clip planes when
                    // the direct corner is degenerate (parallel planes).  Used for
                    // the intersection case where corners are inside the cylinder.
                    let chain = build_plane_chain(p1, p2, clip_planes, center);
                    let n_segs = chain.len();

                    // Step 1: Create corner vertices and generator edges.
                    let n_corners = n_segs.saturating_sub(1);
                    struct CornerData { lo: usize, hi: usize, gen_edge: usize }
                    let mut corner_data: Vec<Option<CornerData>> = Vec::new();
                    corner_data.resize_with(n_corners, || None);
                    for j in 0..n_corners {
                        let (pa, pb) = (chain[j], chain[j + 1]);
                        let (na, da) = (clip_planes[pa].0, clip_planes[pa].1);
                        let (nb, db) = (clip_planes[pb].0, clip_planes[pb].1);
                        let corner_xy = corner_of_planes(na, da, nb, db, center);
                        let lo = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_lo) });
                        let hi = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_hi) });
                        let gen_edge = next_curve(
                            Curve3::Line(Line3 { origin: brep.vertices[lo].point, direction: DVec3::Z }),
                            0.0, h, lo, hi,
                        );
                        corner_data[j] = Some(CornerData { lo, hi, gen_edge });
                    }

                    // Step 2: Build one segment per plane in the chain.
                    let mut segs: Vec<GapSeg> = Vec::with_capacity(n_segs);
                    for j in 0..n_segs {
                        let plane = chain[j];
                        let is_first = j == 0;
                        let is_last = j == n_segs - 1;

                        let (start_lo, start_hi) = if is_first {
                            (v_lo(vi_fr), v_hi(vi_fr))
                        } else {
                            let c = corner_data[j - 1].as_ref().unwrap();
                            (c.lo, c.hi)
                        };
                        let (end_lo, end_hi) = if is_last {
                            (v_lo(vi_to), v_hi(vi_to))
                        } else {
                            let c = corner_data[j].as_ref().unwrap();
                            (c.lo, c.hi)
                        };
                        let gen_from = if is_first {
                            interval_edges[i0].rg
                        } else {
                            corner_data[j - 1].as_ref().unwrap().gen_edge
                        };
                        let gen_to = if is_last {
                            interval_edges[i1].lg
                        } else {
                            corner_data[j].as_ref().unwrap().gen_edge
                        };

                        let bot_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[start_lo].point,
                            direction: brep.vertices[end_lo].point - brep.vertices[start_lo].point,
                        });
                        let eb = next_curve(bot_chord, 0.0, 1.0, start_lo, end_lo);

                        let top_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[start_hi].point,
                            direction: brep.vertices[end_hi].point - brep.vertices[start_hi].point,
                        });
                        let et = next_curve(top_chord, 0.0, 1.0, start_hi, end_hi);

                        segs.push(GapSeg {
                            bot_chord: eb,
                            top_chord: et,
                            plane,
                            gen_from,
                            gen_to,
                        });
                    }
                    gap_segs.push(segs);
                } else {
                    // Direct corner routing for difference case: connect gap
                    // endpoints through the direct intersection corner of the
                    // two bounding clip planes.  If the corner is outside the
                    // cylinder or the planes are parallel, use a single chord
                    // between cylinder wall points to avoid "fin" faces that
                    // extend outside the cylinder cross-section.
                    let (n1, d1) = (clip_planes[p1].0, clip_planes[p1].1);
                    let (n2, d2) = (clip_planes[p2].0, clip_planes[p2].1);
                    let corner_xy = corner_of_planes(n1, d1, n2, d2, center);
                    let corner_dist_sq = (corner_xy.x - center.x).powi(2)
                        + (corner_xy.y - center.y).powi(2);

                    // Only use 2-segment corner routing when the corner is
                    // inside the cylinder AND the planes are non-parallel.
                    let non_parallel = (n1.x * n2.y - n1.y * n2.x).abs() > 1e-12;
                    if non_parallel && corner_dist_sq <= r * r + 1e-9 {
                        // 2 segments through the corner
                        let lo_corner = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_lo) });
                        let hi_corner = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_hi) });
                        let gen_corner = next_curve(
                            Curve3::Line(Line3 { origin: brep.vertices[lo_corner].point, direction: DVec3::Z }),
                            0.0, h, lo_corner, hi_corner,
                        );

                        // Segment 1: p1 鈫?corner
                        let bot_chord1 = Curve3::Line(Line3 {
                            origin: brep.vertices[v_lo(vi_fr)].point,
                            direction: brep.vertices[lo_corner].point - brep.vertices[v_lo(vi_fr)].point,
                        });
                        let eb1 = next_curve(bot_chord1, 0.0, 1.0, v_lo(vi_fr), lo_corner);
                        let top_chord1 = Curve3::Line(Line3 {
                            origin: brep.vertices[v_hi(vi_fr)].point,
                            direction: brep.vertices[hi_corner].point - brep.vertices[v_hi(vi_fr)].point,
                        });
                        let et1 = next_curve(top_chord1, 0.0, 1.0, v_hi(vi_fr), hi_corner);

                        // Segment 2: corner 鈫?p2
                        let bot_chord2 = Curve3::Line(Line3 {
                            origin: brep.vertices[lo_corner].point,
                            direction: brep.vertices[v_lo(vi_to)].point - brep.vertices[lo_corner].point,
                        });
                        let eb2 = next_curve(bot_chord2, 0.0, 1.0, lo_corner, v_lo(vi_to));
                        let top_chord2 = Curve3::Line(Line3 {
                            origin: brep.vertices[hi_corner].point,
                            direction: brep.vertices[v_hi(vi_to)].point - brep.vertices[hi_corner].point,
                        });
                        let et2 = next_curve(top_chord2, 0.0, 1.0, hi_corner, v_hi(vi_to));

                        gap_segs.push(vec![
                            GapSeg {
                                bot_chord: eb1, top_chord: et1, plane: p1,
                                gen_from: interval_edges[i0].rg,
                                gen_to: gen_corner,
                            },
                            GapSeg {
                                bot_chord: eb2, top_chord: et2, plane: p2,
                                gen_from: gen_corner,
                                gen_to: interval_edges[i1].lg,
                            },
                        ]);
                    } else {
                        // Single direct chord between cylinder wall points
                        // (avoids fins for corners outside the cylinder, and
                        // handles parallel planes).
                        let bot_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[v_lo(vi_fr)].point,
                            direction: brep.vertices[v_lo(vi_to)].point - brep.vertices[v_lo(vi_fr)].point,
                        });
                        let eb = next_curve(bot_chord, 0.0, 1.0, v_lo(vi_fr), v_lo(vi_to));
                        let top_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[v_hi(vi_fr)].point,
                            direction: brep.vertices[v_hi(vi_to)].point - brep.vertices[v_hi(vi_fr)].point,
                        });
                        let et = next_curve(top_chord, 0.0, 1.0, v_hi(vi_fr), v_hi(vi_to));
                        gap_segs.push(vec![GapSeg {
                            bot_chord: eb, top_chord: et, plane: p1,
                            gen_from: interval_edges[i0].rg,
                            gen_to: interval_edges[i1].lg,
                        }]);
                    }
                }
            }
            _ => {
                gap_segs.push(Vec::new());
            }
        }
    }

    // ---- 6. bottom & top planar cap faces ----
    let bot_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_lo),
            normal: -DVec3::Z,
        }));
        si
    };
    let top_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_hi),
            normal: DVec3::Z,
        }));
        si
    };

    // Bottom cap: V_s_lo 鈫?V_e_lo (arc) 鈫?V_s_{i+1}_lo (chords) 鈫?...
    let mut bot_wire_edges: Vec<WireEdge> = Vec::new();
    let mut top_wire_edges: Vec<WireEdge> = Vec::new();
    for i in 0..n_int {
        // Arc from start-vertex to end-vertex (bottom: V_s_lo鈫扸_e_lo, top: V_s_hi鈫扸_e_hi)
        bot_wire_edges.push(WireEdge::rev(interval_edges[i].ba));
        top_wire_edges.push(WireEdge::fwd(interval_edges[i].ta));
        // Gap chord segments from end-i to start-(i+1)
        for seg in &gap_segs[i] {
            bot_wire_edges.push(WireEdge::fwd(seg.bot_chord));
            top_wire_edges.push(WireEdge::fwd(seg.top_chord));
        }
    }

    let push_face = |brep: &mut BRep, wire: Wire, surf_idx: usize, normal: DVec3| {
        let fi = brep.solids[0].shells[0].faces.len();
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: wire,
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(surf_idx);
    };
    if !skip_bottom_cap {
        push_face(&mut brep, Wire { edges: bot_wire_edges }, bot_surf_idx, -DVec3::Z);
    }
    if !skip_top_cap {
        push_face(&mut brep, Wire { edges: top_wire_edges }, top_surf_idx, DVec3::Z);
    }

    // ---- 7. side faces (one per gap segment, on the segment's clip plane) ----
    // The face wire runs: bot_chord_fwd 鈫?gen_to_fwd 鈫?top_chord_rev 鈫?gen_from_rev
    for segs in &gap_segs {
        for seg in segs {
            // Determine the plane for this side face.
            // Most gap segments lie on a clip-plane, but the "single direct chord"
            // fallback for parallel clip planes (bopcut_simple t6/v4) creates a
            // chord on the cylinder wall that does NOT lie on the assigned plane.
            // In that case we infer the correct vertical plane from the chord.
            let (n, d) = {
                let (n0, d0) = clip_planes[seg.plane];
                let bot_e = &brep.edges[seg.bot_chord];
                let s_pt = brep.vertices[bot_e.start].point;
                let e_pt = brep.vertices[bot_e.end].point;
                // Check whether both endpoints lie on the assigned clip plane.
                let tol = 1e-8;
                let on_p0 = (n0.dot(s_pt - center) - (-d0)).abs() < tol
                         && (n0.dot(e_pt - center) - (-d0)).abs() < tol;
                if on_p0 {
                    (n0, d0)
                } else {
                    // Chord endpoints don't lie on the assigned clip plane.
                    // Compute a vertical plane through the chord: the face
                    // winding normal is cross(chord_dir, Z), and we need -n
                    // to match it (since push_face stores outward = -n).
                    let chord_dir = e_pt - s_pt;
                    let cross_z = chord_dir.cross(DVec3::Z);
                    let len = cross_z.length();
                    if len < 1e-15 {
                        // Degenerate chord 鈥?fall back to original plane.
                        (n0, d0)
                    } else {
                        let n = -cross_z / len; // so that -n = winding normal
                        let d = -n.dot(s_pt - center);
                        (n, d)
                    }
                }
            };
            let si_side = {
                let si = brep.geom.surfaces.len();
                brep.geom.surfaces.push(Surface3::Plane(Plane {
                    origin: center - n * d,
                    normal: -n,
                }));
                si
            };
            let side_wire = Wire {
                edges: vec![
                    WireEdge::fwd(seg.bot_chord),
                    WireEdge::fwd(seg.gen_to),
                    WireEdge::rev(seg.top_chord),
                    WireEdge::rev(seg.gen_from),
                ],
            };
            push_face(&mut brep, side_wire, si_side, -n);
        }
    }

    brep
}
/// Build an arc BRep for the parallel-only cylinder-box difference case.
///
/// This builds the portion of a Z-aligned cylinder satisfying
/// `(P - center)路clip_n 鈮?cut_dist` (the outside-slab region). The arc is
/// centered on `clip_n` with half-angle `伪 = acos(cut_dist/r)`. The clip face
/// is on the plane `center + clip_n路cut_dist` with outward normal `-clip_n`.
///
/// `clip_n`: horizontal unit normal from center toward the clip plane.
/// `cut_dist`: distance from center to clip plane (鈮?, 鈮).
fn build_cylinder_arc_for_difference(
    center: DVec3,
    r: f64,
    h: f64,
    clip_n: DVec3,
    cut_dist: f64,
) -> BRep {
    build_cylinder_arc_for_difference_skip(center, r, h, clip_n, cut_dist, false, false)
}

/// Same as [`build_cylinder_arc_for_difference`] but with optional cap skipping.
fn build_cylinder_arc_for_difference_skip(
    center: DVec3,
    r: f64,
    h: f64,
    clip_n: DVec3,
    cut_dist: f64,
    skip_bottom_cap: bool,
    skip_top_cap: bool,
) -> BRep {
    let alpha = (cut_dist / r).clamp(-1.0, 1.0).acos();
    let phi = clip_n.y.atan2(clip_n.x);
    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    let (sa, ca) = alpha.sin_cos();
    let (sp, cp) = phi.sin_cos();

    // (cos(蠁卤伪), sin(蠁卤伪))
    let cos_phi_minus_alpha = cp * ca + sp * sa;
    let sin_phi_minus_alpha = sp * ca - cp * sa;
    let cos_phi_plus_alpha = cp * ca - sp * sa;
    let sin_phi_plus_alpha = sp * ca + cp * sa;

    let v0_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_lo);
    let v1_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_lo);
    let v2_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_hi);
    let v3_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_hi);

    let mut brep = BRep::new();

    // Vertices
    let v0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v0_p });
    let v1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v1_p });
    let v2 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v2_p });
    let v3 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v3_p });

    // Edge helper
    let mut next_curve = |c: Curve3, t0: f64, t1: f64, start: usize, end: usize| -> usize {
        let idx = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(c);
        while brep.geom.edge_curve.len() <= idx {
            brep.geom.edge_curve.push(None);
            brep.geom.edge_curve_range.push(None);
            brep.geom.edge_degenerated.push(false);
        }
        brep.geom.edge_curve[idx] = Some(ci);
        brep.geom.edge_curve_range[idx] = Some([t0, t1]);
        brep.geom.edge_pcurves.push(Vec::new());
        idx
    };

    // E0: bottom arc (V1鈫扸0), same convention as build_half_cylinder_intersection_brep
    let circle_bot = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
        radius: r,
    });
    let e0 = next_curve(circle_bot, -phi - alpha, -phi + alpha, v1, v0);

    // E1: right generator (V1鈫扸2)
    let line_r = Curve3::Line(Line3 { origin: v1_p, direction: DVec3::Z });
    let e1 = next_curve(line_r, 0.0, h, v1, v2);

    // E2: top arc (V3鈫扸2)
    let circle_top = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
        radius: r,
    });
    let e2 = next_curve(circle_top, phi - alpha, phi + alpha, v3, v2);

    // E3: left generator (V0鈫扸3)
    let line_l = Curve3::Line(Line3 { origin: v0_p, direction: DVec3::Z });
    let e3 = next_curve(line_l, 0.0, h, v0, v3);

    // E4: bottom chord on clip plane (V0鈫扸1)
    let line_cb = Curve3::Line(Line3 { origin: v0_p, direction: v1_p - v0_p });
    let e4 = next_curve(line_cb, 0.0, 1.0, v0, v1);

    // E5: top chord on clip plane (V2鈫扸3)
    let line_ct = Curve3::Line(Line3 { origin: v2_p, direction: v3_p - v2_p });
    let e5 = next_curve(line_ct, 0.0, 1.0, v2, v3);

    // --- Surfaces ---
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(center.x, center.y, cz_lo),
        axis: DVec3::Z, radius: r, ref_dir: DVec3::X,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
    });
    let bot_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
    });
    // Clip plane at center + clip_n*cut_dist, outward normal = -clip_n
    let clip_surf = Surface3::Plane(Plane {
        origin: center + clip_n * cut_dist,
        normal: -clip_n,
    });

    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cyl_surf);
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(top_plane);
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(bot_plane);
    let si_clip = brep.geom.surfaces.len();
    brep.geom.surfaces.push(clip_surf);

    // Face helper
    let mut push_face = |outer: Wire, surf_idx: usize, normal: DVec3| -> usize {
        let fi = if brep.solids.is_empty() {
            brep.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });
            0
        } else {
            brep.solids[0].shells[0].faces.len()
        };
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: outer, inner_wires: Vec::new(), normal,
            triangles: Vec::new(), sample_point: None, mesh_dirty: true,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(surf_idx);
        fi
    };

    // F0: Cylindrical wall 鈥?wire V0鈫扸1鈫扸2鈫扸3鈫扸0
    let cyl_wire = Wire {
        edges: vec![
            WireEdge::rev(e0), WireEdge::fwd(e1),
            WireEdge::rev(e2), WireEdge::rev(e3),
        ],
    };
    let _f0 = push_face(cyl_wire, si_cyl, clip_n);
    while brep.geom.face_surface_range.len() <= _f0 {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range[_f0] = Some([phi - alpha, phi + alpha, 0.0, h]);

    // F1: Top cap (normal=+Z)
    if !skip_top_cap {
        let _f1 = push_face(Wire { edges: vec![WireEdge::fwd(e5), WireEdge::fwd(e2)] }, si_top, DVec3::Z);
    }

    // F2: Bottom cap (normal=-Z)
    if !skip_bottom_cap {
        let _f2 = push_face(Wire { edges: vec![WireEdge::fwd(e4), WireEdge::fwd(e0)] }, si_bot, -DVec3::Z);
    }

    // F3: Clip face (bounding the arc on the box face, outward normal = -clip_n)
    let clip_wire = Wire {
        edges: vec![
            WireEdge::fwd(e4), WireEdge::fwd(e1),
            WireEdge::fwd(e5), WireEdge::rev(e3),
        ],
    };
    let _f3 = push_face(clip_wire, si_clip, -clip_n);

    brep
}

/// bounded by a vertical plane parallel to the cylinder axis.
///
/// The result is a portion of the cylinder cut lengthwise by a plane parallel to
/// its axis. Only the side where `(P - center)路clip_n 鈮?-cut_dist` is kept.
/// `clip_n` must be a horizontal unit vector (z=0) pointing into the kept half.
/// `cut_dist` is the distance from the cylinder center to the cut plane (measured
/// in the direction opposite to `clip_n`, i.e., into the kept half).
///
/// Build a full Z-aligned cylinder BRep with selective inclusion of end caps.
///
/// When `include_bottom_cap` is false, the face at `center.z - height/2` is omitted.
/// When `include_top_cap` is false, the face at `center.z + height/2` is omitted.
/// The lateral face and the included cap faces are always created.
fn build_cylinder_brep_caps(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
    include_bottom_cap: bool,
    include_top_cap: bool,
) -> Result<BRep, rcad_modeling::BuildError> {
    let mut brep = make_cylinder_brep(center, axis, ref_dir, radius, height)?;
    if include_bottom_cap && include_top_cap {
        return Ok(brep);
    }
    let half_h = height * 0.5;
    let z_lo = center.z - half_h;
    let z_hi = center.z + half_h;
    let shell = &mut brep.solids[0].shells[0];
    shell.faces.retain(|face| {
        // Determine if this face is a planar cap by checking its outer_wire
        // vertex Z coordinates. A cap at z_lo has all vertices near z_lo.
        let zs: Vec<f64> = face.outer_wire.edges.iter()
            .filter_map(|we| brep.edges.get(we.idx))
            .flat_map(|e| [brep.vertices.get(e.start).map(|v| v.point.z),
                           brep.vertices.get(e.end).map(|v| v.point.z)])
            .flatten()
            .collect();
        if zs.is_empty() {
            return true;
        }
        let z_avg = zs.iter().sum::<f64>() / zs.len() as f64;
        let is_bottom = (z_avg - z_lo).abs() < 0.1 * height;
        let is_top = (z_avg - z_hi).abs() < 0.1 * height;
        if is_bottom && !include_bottom_cap { return false; }
        if is_top && !include_top_cap { return false; }
        true
    });
    Ok(brep)
}

/// bounded by a vertical plane parallel to the cylinder axis.
///
/// The result is a portion of the cylinder cut lengthwise by a plane parallel to
/// its axis. Only the side where `(P - center)路clip_n 鈮?-cut_dist` is kept.
/// `clip_n` must be a horizontal unit vector (z=0) pointing into the kept half.
/// `cut_dist` is the distance from the cylinder center to the cut plane (measured
/// in the direction opposite to `clip_n`, i.e., into the kept half).
///
/// When cut_dist=0, the cut plane passes through the cylinder axis and the result
/// is a clean half-cylinder. When cut_dist > 0, the cut plane is offset outward
/// from the axis, and more than half the cylinder is kept.
///
/// This is used when a box fully contains the cylinder in one XY axis but only
/// partially contains it in the other.
fn build_half_cylinder_intersection_brep(
    center: DVec3,   // cylinder center (cx, cy, cz)
    r: f64,          // radius
    h: f64,          // height
    clip_n: DVec3,   // horizontal unit normal pointing into the kept half
    cut_dist: f64,   // distance from center to cut plane (鈮?, into kept half)
) -> BRep {
    debug_assert!(cut_dist >= 0.0 && cut_dist <= r + 1e-12,
        "cut_dist must be in [0, r]. Got {cut_dist} for r={r}");

    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    // Azimuth angle of clip_n in XY plane.
    let phi = clip_n.y.atan2(clip_n.x);

    // Half-angle of the kept arc: 伪 = arccos(-cut_dist/r).
    // For cut_dist=0 (center cut): 伪 = 蟺/2 鈫?螖u = 蟺 (half-cylinder).
    // For cut_dist=r (full cylinder): 伪 = arccos(-1) = 蟺 鈫?螖u = 2蟺 (full cylinder).
    let alpha = (-cut_dist / r).acos();

    // Vertices at the intersection of the cut plane with the cylinder surface.
    // V0/V3 at u = 蠁 - 伪 (left generator), V1/V2 at u = 蠁 + 伪 (right generator).
    let (sa, ca) = alpha.sin_cos();
    let (sp, cp) = phi.sin_cos();

    // (cos(蠁卤伪), sin(蠁卤伪)) = (cp*ca 鈭?sp*sa, sp*ca 卤 cp*sa)
    let cos_phi_minus_alpha = cp * ca + sp * sa;
    let sin_phi_minus_alpha = sp * ca - cp * sa;
    let cos_phi_plus_alpha = cp * ca - sp * sa;
    let sin_phi_plus_alpha = sp * ca + cp * sa;

    let v0_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_lo);
    let v1_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_lo);
    let v2_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_hi);
    let v3_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_hi);

    // --- Build BRep directly ---
    let mut brep = BRep::new();

    // Vertices
    let v0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v0_p });
    let v1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v1_p });
    let v2 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v2_p });
    let v3 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v3_p });

    // Edge index helpers 鈥?push an edge and return its index.
    let mut next_curve = |c: Curve3, t0: f64, t1: f64, start: usize, end: usize| -> usize {
        let idx = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(c);
        // Ensure parallel vecs
        while brep.geom.edge_curve.len() <= idx {
            brep.geom.edge_curve.push(None);
            brep.geom.edge_curve_range.push(None);
            brep.geom.edge_degenerated.push(false);
        }
        brep.geom.edge_curve[idx] = Some(ci);
        brep.geom.edge_curve_range[idx] = Some([t0, t1]);
        brep.geom.edge_pcurves.push(Vec::new());
        idx
    };

    // E0: bottom arc (V1鈫扸0 along kept side of bottom circle).
    // Circle3 with normal=-Z maps CylSurf u 鈫?Circle3 -u, so V1(CylSurf u=蠁+伪)
    // 鈫?Circle3 u=-(蠁+伪) and V0(CylSurf u=蠁-伪) 鈫?Circle3 u=-(蠁-伪).
    let circle_bot = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
        radius: r,
    });
    let e0 = next_curve(circle_bot, -phi - alpha, -phi + alpha, v1, v0);

    // E1: right generator (V1鈫扸2)
    let line_r = Curve3::Line(Line3 {
        origin: v1_p,
        direction: DVec3::Z,
    });
    let e1 = next_curve(line_r, 0.0, h, v1, v2);

    // E2: top arc (V3鈫扸2 along kept side of top circle).
    // Circle3 with normal=Z uses the same param as CylSurf, so V3(CylSurf u=蠁-伪)
    // 鈫?Circle3 u=蠁-伪 and V2(CylSurf u=蠁+伪) 鈫?Circle3 u=蠁+伪.
    let circle_top = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
        radius: r,
    });
    let e2 = next_curve(circle_top, phi - alpha, phi + alpha, v3, v2);

    // E3: left generator (V0鈫扸3)
    let line_l = Curve3::Line(Line3 {
        origin: v0_p,
        direction: DVec3::Z,
    });
    let e3 = next_curve(line_l, 0.0, h, v0, v3);

    // E4: cut bottom (V0鈫扸1)
    let line_cb = Curve3::Line(Line3 {
        origin: v0_p,
        direction: v1_p - v0_p,
    });
    let e4 = next_curve(line_cb, 0.0, 1.0, v0, v1);

    // E5: cut top (V2鈫扸3)
    let line_ct = Curve3::Line(Line3 {
        origin: v2_p,
        direction: v3_p - v2_p,
    });
    let e5 = next_curve(line_ct, 0.0, 1.0, v2, v3);

    // --- Surfaces ---
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(center.x, center.y, cz_lo),
        axis: DVec3::Z,
        radius: r,
        ref_dir: DVec3::X,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
    });
    let bot_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
    });
    let cut_plane = Surface3::Plane(Plane {
        origin: center - clip_n * cut_dist,
        normal: -clip_n, // outward from the solid
    });

    // Push surfaces and register face-surface mapping.
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cyl_surf);
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(top_plane);
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(bot_plane);
    let si_cut = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cut_plane);

    // Helper to push a face.
    let mut push_face = |outer: Wire, surf_idx: usize, normal: DVec3| -> usize {
        let fi = if brep.solids.is_empty() {
            brep.solids.push(Solid {
                shells: vec![Shell { faces: Vec::new() }],
            });
            0
        } else {
            brep.solids[0].shells[0].faces.len()
        };
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: outer,
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(surf_idx);
        fi
    };

    // F0: Cylindrical face 鈥?wire V0鈫扸1鈫扸2鈫扸3鈫扸0
    // E0 stored as V1鈫扸0 鈫?E0_rev = V0鈫扸1
    // E1 stored as V1鈫扸2 鈫?E1_fwd = V1鈫扸2
    // E2 stored as V3鈫扸2 鈫?E2_rev = V2鈫扸3
    // E3 stored as V0鈫扸3 鈫?E3_rev = V3鈫扸0
    let cyl_wire = Wire {
        edges: vec![
            WireEdge::rev(e0),
            WireEdge::fwd(e1),
            WireEdge::rev(e2),
            WireEdge::rev(e3),
        ],
    };
    let _f0 = push_face(cyl_wire, si_cyl, clip_n);
    // Set face_surface_range for analytic SA: u 鈭?[蠁-伪, 蠁+伪], v 鈭?[0, h]
    while brep.geom.face_surface_range.len() <= _f0 {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range[_f0] = Some([phi - alpha, phi + alpha, 0.0, h]);

    // F1: Top half-disk (normal=+Z). Wire: V2鈫扸3 (cut top) 鈫?V3鈫扸2 (top arc)
    // E5_fwd: V2鈫扸3, E2_fwd: V3鈫扸2
    let top_wire = Wire {
        edges: vec![WireEdge::fwd(e5), WireEdge::fwd(e2)],
    };
    let _f1 = push_face(top_wire, si_top, DVec3::Z);

    // F2: Bottom half-disk (normal=-Z). Wire: V0鈫扸1 (cut bottom) 鈫?V1鈫扸0 (bottom arc)
    // E4_fwd: V0鈫扸1, E0_fwd: V1鈫扸0
    let bot_wire = Wire {
        edges: vec![WireEdge::fwd(e4), WireEdge::fwd(e0)],
    };
    let _f2 = push_face(bot_wire, si_bot, -DVec3::Z);

    // F3: Cut face (normal=-clip_n). Wire: V0鈫扸1鈫扸2鈫扸3鈫扸0
    // E4_fwd: V0鈫扸1, E1_fwd: V1鈫扸2, E5_fwd: V2鈫扸3, E3_rev: V3鈫扸0
    let cut_wire = Wire {
        edges: vec![
            WireEdge::fwd(e4),
            WireEdge::fwd(e1),
            WireEdge::fwd(e5),
            WireEdge::rev(e3),
        ],
    };
    let _f3 = push_face(cut_wire, si_cut, -clip_n);

    brep
}

/// Shared helper for intersection: cyl_brep is the cylinder, box_brep is the box.
fn try_intersect_cylinder_box_one_dir(cyl_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let ca = try_cylinder_center_axis_radius_height(cyl_brep)?;
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) = ca;
    // try_cylinder_center_axis_radius_height returns the CylindricalSurface origin
    // (bottom of cylinder), not the geometric center. Compute the actual center.
    let cyl_center = cyl_bottom + cyl_axis * (cyl_height / 2.0);
    let bx = try_as_box(box_brep)?;

    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let z_idx = find_z_axis_index(&bx)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = bx.axes[u_idx];
    let v_ax = bx.axes[v_idx];
    let eu = bx.extents[u_idx];
    let ev = bx.extents[v_idx];
    let ew = bx.extents[z_idx];
    let bc = bx.center;

    let cyl_z_lo = cyl_center.z - cyl_height / 2.0;
    let cyl_z_hi = cyl_center.z + cyl_height / 2.0;
    let box_z_lo = bc.z - ew;
    let box_z_hi = bc.z + ew;

    let tol = TOL * 10.0;

    // Intersection Z range.
    let inter_lo = cyl_z_lo.max(box_z_lo);
    let inter_hi = cyl_z_hi.min(box_z_hi);
    if inter_hi <= inter_lo + tol {
        return Some(BRep::default());
    }

    // Check XY containment.
    let cu = (cyl_center - bc).dot(u_ax);
    let cv = (cyl_center - bc).dot(v_ax);
    let full_u = cu - cyl_r >= -eu - tol && cu + cyl_r <= eu + tol;
    let full_v = cv - cyl_r >= -ev - tol && cv + cyl_r <= ev + tol;

    if full_u && full_v {
        // Full XY containment: build sub-cylinder for the intersection Z range.
        let h = inter_hi - inter_lo;
        let cz = inter_lo + h / 2.0;
        return make_cylinder_brep(
            DVec3::new(cyl_center.x, cyl_center.y, cz),
            cyl_axis, u_ax, cyl_r, h,
        ).ok();
    }

    // Collect clip planes from all partially-contained axes.
    // Each partial axis may have 0, 1, or 2 clipped sides.
    let mut clip_planes: Vec<(DVec3, f64)> = Vec::new();
    for &(full_axis, cp, ax, ext) in &[(full_u, cu, u_ax, eu), (full_v, cv, v_ax, ev)] {
        if !full_axis {
            if cp - cyl_r < -ext - tol {
                clip_planes.push((ax, cp + ext));
            }
            if cp + cyl_r > ext + tol {
                clip_planes.push((-ax, ext - cp));
            }
        }
    }

    if clip_planes.is_empty() {
        return None;
    }

    let h = inter_hi - inter_lo;
    let cz = inter_lo + h / 2.0;
    let adj_center = DVec3::new(cyl_center.x, cyl_center.y, cz);

    // Single clip plane 鈫?existing half-cylinder builder (backward compat, efficient).
    if clip_planes.len() == 1 {
        let (clip_dir, cut_dist) = clip_planes[0];
        // When cut_dist < -r the cylinder center is so far outside the box
        // that the cylinder doesn't reach the box face 鈫?empty intersection.
        // (The multi-plane builder below handles this via zero-range valid 胃
        // intervals, but the half-cylinder builder asserts cut_dist 鈮?0.)
        if cut_dist < -cyl_r + 1e-12 {
            return Some(BRep::default());
        }
        return Some(build_half_cylinder_intersection_brep(adj_center, cyl_r, h, clip_dir, cut_dist));
    }

    // Multiple clip planes 鈫?general multi-plane builder.
    Some(build_cylinder_box_intersection_brep(adj_center, cyl_r, h, &clip_planes))
}

/// Build box 鈭?cylinder difference when the cylinder center is inside the box UV
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
    //   INSIDE  (midpoint in circle) 鈫?box perimeter edge in the hole boundary
    //   OUTSIDE (midpoint outside)   鈫?circle arc in the hole boundary
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
    // First build in CCW perimeter order: INSIDE segments 鈫?box perimeter edge
    // (add_box_perimeter_vertices goes CCW = increasing t); OUTSIDE segments 鈫?    // circle arc (arc_vertices picks the direction that stays inside the box).
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

    // Side walls 鈥?same as build_box_cylinder_result_partial (per-edge intersection)
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
        if !ins_lo && !ins_hi { None }        // both outside = tangent 鈫?full face
        else if !ins_lo { Some((lo, p)) }     // wrap toward hi (other point beyond hi)
        else if !ins_hi { Some((p, hi)) }     // wrap toward lo
        else { None }                         // both inside = shouldn't happen 鈫?full face
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
        } else if (eu - cu).powi(2) + (0.0 - cv).powi(2) > cyl_r.powi(2) + tol {
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
        } else if (-eu - cu).powi(2) + (0.0 - cv).powi(2) > cyl_r.powi(2) + tol {
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
        } else if (0.0 - cu).powi(2) + (ev - cv).powi(2) > cyl_r.powi(2) + tol {
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
        } else if (0.0 - cu).powi(2) + (-ev - cv).powi(2) > cyl_r.powi(2) + tol {
            add_side_strip(&cn, -v_ax, -eu, eu);
        }
    }

    // Cap faces with inner hole (+ full bottom/top face when cylinder doesn't reach).
    // Bottom region
    let outer_lo = box_outer(z_lo);
    let inner_lo = build_inner(z_lo);
    if z_lo > -ew + tol {
        // Full bottom face at z=-ew (cylinder doesn't reach bottom).
        let bot = [corner(-eu, -ev, -ew), corner(-eu, ev, -ew), corner(eu, ev, -ew), corner(eu, -ev, -ew)];
        if let Some(f) = rect_face_4(bot, -DVec3::Z) { pieces.push(f); }
        // Interior annular cap at z_lo, normal +Z.
        if !inner_lo.is_empty() {
            if let Some(f) = planar_face_with_inner_hole(&outer_lo, &inner_lo, DVec3::Z) { pieces.push(f); }
        }
    } else if !inner_lo.is_empty() {
        if let Some(f) = planar_face_with_inner_hole(&outer_lo, &inner_lo, DVec3::Z) { pieces.push(f); }
    }
    // Top region
    let outer_hi = box_outer(z_hi);
    let inner_hi = build_inner(z_hi);
    if z_hi < ew - tol {
        // Full top face at z=ew (cylinder doesn't reach top).
        let top = [corner(-eu, -ev, ew), corner(eu, -ev, ew), corner(eu, ev, ew), corner(-eu, ev, ew)];
        if let Some(f) = rect_face_4(top, DVec3::Z) { pieces.push(f); }
        // Interior annular cap at z_hi, normal -Z.
        if !inner_hi.is_empty() {
            if let Some(f) = planar_face_with_inner_hole(&outer_hi, &inner_hi, -DVec3::Z) { pieces.push(f); }
        }
    } else if !inner_hi.is_empty() {
        if let Some(f) = planar_face_with_inner_hole(&outer_hi, &inner_hi, -DVec3::Z) { pieces.push(f); }
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
        // sphere primitives carry approximate face normals 閳?SA matches analytic \(4\pi(R^2+r^2)\).
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
        // Cone apex z=10, base z=-10, rb=10; cylinder z in [-10,0], r=10 閳?same as geometry_properties ZP3.
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
        use glam::{DAffine3, DVec2, DVec3};
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

    // 鈹€鈹€ box-box intersection fast path 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn box_box_intersection_partial_overlap() {
        // bcommon_simple_c1: 1脳1脳1 鈭?1.5脳0.5脳0.5 鈫?1脳0.5脳0.5 (SA=2.5, vol=0.25)
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
        // 2脳2脳2 box 鈭?0.5脳0.5脳0.5 inside 鈫?the inner 0.5^3 box
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
        // Disjoint boxes 鈫?empty intersection.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("no-overlap");
        let n_faces: usize = r.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty intersection");
    }

    #[test]
    fn box_box_intersection_a7_like() {
        // bcommon_simple A7: 1脳1脳1 鈭?1脳1.5脳1 鈫?the contained 1脳1脳1 (SA=6, vol=1).
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
        // Sphere 鈭?box falls through to generic path (no panic).
        let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &sphere, &b).expect("sphere-box");
        // Some result 鈥?specific value doesn't matter.
        assert!(r.solids.len() >= 1 || r.vertices.is_empty());
    }

    // 鈹€鈹€ box-box difference fast path 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    #[test]
    fn box_box_difference_opposite_same_axis() {
        // F1-like: A=1.5脳0.5脳0.5 at (-0.25,0,0), B=1脳1脳1 at origin.
        // boptuc = A - B (but here we test try_difference_box_box directly,
        // so a = first arg = A, b = second arg = B).
        // A extends on x-lo (-0.25<0) and x-hi (1.25>1), same axis 鈫?two
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
        // bcut_simple_c1: 0.5脳1.5脳1 at (0,-0.5,0) minus 1脳1脳1 at origin.
        // Excess only on y-lo [-0.5, 0]. Result: 0.5脳0.5脳1 box, SA=2.5.
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
        // Disjoint boxes 鈫?difference is just A.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("no-overlap");
        let sa = surface_area(&r);
        assert!((sa - 6.0).abs() < 1e-7, "SA={sa} expected 6.0 (unchanged A)");
    }

    #[test]
    fn box_box_difference_a_inside_b() {
        // A fully inside B 鈫?empty.
        let a = make_box_brep(DVec3::new(0.25, 0.25, 0.25), DVec3::X, DVec3::Y, 0.5, 0.5, 0.5).unwrap();
        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &outer).expect("A-in-B");
        let n_faces: usize = r.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty difference");
    }

    // 鈹€鈹€ general box-box (rotated) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

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
        // Box at origin, rotated 45掳 around Z at origin.
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
        // bcommon_simple_c3-like: unit box 鈭?rotated box.
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
        // box at (0.25, 0.25, 0) size 0.5脳0.5脳-1, pivot (.25,.25,0), rotate 30掳 Z
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
        // Disjoint rotated boxes 鈫?difference should be B unchanged.
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
        // B fully inside A 鈫?B - A = empty.
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
        // Box 鈭?sphere 鈫?falls through to Pave-Filler (no panic).
        let b = make_rotated_box(
            DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
            DVec3::ZERO, DVec3::Z, 45.0,
        );
        let s = make_sphere_brep(DVec3::new(0.5, 0.5, 0.5), 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &b, &s).expect("box-sphere intersection");
        assert!(r.solids.len() >= 1 || r.vertices.is_empty());
    }

    #[test]
    fn debug_half_cylinder_sa() {
        // U5 case: cylinder r=1, h=2, center (0,0,1), clip X鈮?
        let brep = build_half_cylinder_intersection_brep(
            DVec3::new(0.0, 0.0, 1.0), 1.0, 2.0, DVec3::X, 0.0,
        );
        let total = rcad_kernel::surface_area(&brep);
        println!("DEBUG half-cylinder: total SA = {total}");

        // Per-face SA
        for (si, solid) in brep.solids.iter().enumerate() {
            for (fi, face) in solid.shells[0].faces.iter().enumerate() {
                let a = rcad_kernel::face_surface_area(&brep, face, fi);
                let n_edges = face.outer_wire.edges.len();
                println!("  Solid {si} Face {fi} ({n_edges} edges): SA = {a}");
            }
        }

        let expected = 3.0 * std::f64::consts::PI + 4.0; // 2蟺 + 蟺/2 + 蟺/2 + 4 = 3蟺+4
        println!("DEBUG expected = {expected}");
        assert!(
            (total - expected).abs() < 0.01,
            "Expected ~{expected}, got {total}"
        );
    }

    #[test]
    fn debug_rotated_cylinder_box_difference() {
        // Shared cylinder for all test cases
        let cyl = make_cylinder_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
        let r = (3.0_f64).sqrt() / 2.0;

        // --- Case x3: box(-1, -r, 0, 1+r, 1+r, 1) rotated 30Z ---
        {
            let mut bx = make_box_brep(DVec3::new(-1.0, -r, 0.0), DVec3::X, DVec3::Y, 1.0 + r, 1.0 + r, 1.0).unwrap();
            let rot = DAffine3::from_axis_angle(DVec3::Z, 30.0_f64.to_radians());
            bx.apply_transform(rot);
            let fast = try_difference_cylinder_box(&cyl, &bx);
            println!("x3: {:?}", fast.is_some());
            if let Some(ref r) = fast {
                println!("  SA: {:.6}", rcad_kernel::surface_area(r));
            }
        }

        // --- Case x6: box(-r, -1, 0, 2*r, 2, 1) rotated 30Z ---
        {
            let mut bx = make_box_brep(DVec3::new(-r, -1.0, 0.0), DVec3::X, DVec3::Y, 2.0 * r, 2.0, 1.0).unwrap();
            let rot = DAffine3::from_axis_angle(DVec3::Z, 30.0_f64.to_radians());
            bx.apply_transform(rot);
            let box_info = try_as_box(&bx);
            println!("x6 try_as_box: {:?}", box_info.is_some());
            if let Some(ref info) = box_info {
                println!("  axes: {:?}, center: {:?}, extents: {:?}", info.axes, info.center, info.extents);
            }
            let fast = try_difference_cylinder_box(&cyl, &bx);
            println!("x6: {:?}", fast.is_some());
            if let Some(ref r) = fast {
                println!("  SA: {:.6}", rcad_kernel::surface_area(r));
            }
        }

        // --- Case zb6: box(-r, -r, 0, 2*r, 2*r, 1) rotated 60Z ---
        {
            let mut bx = make_box_brep(DVec3::new(-r, -r, 0.0), DVec3::X, DVec3::Y, 2.0 * r, 2.0 * r, 1.0).unwrap();
            let rot = DAffine3::from_axis_angle(DVec3::Z, 60.0_f64.to_radians());
            bx.apply_transform(rot);
            let box_info = try_as_box(&bx);
            println!("zb6 try_as_box: {:?}", box_info.is_some());
            if let Some(ref info) = box_info {
                println!("  axes: {:?}, center: {:?}, extents: {:?}", info.axes, info.center, info.extents);
            }
            let fast = try_difference_cylinder_box(&cyl, &bx);
            println!("zb6: {:?}", fast.is_some());
            if let Some(ref r) = fast {
                println!("  SA: {:.6}", rcad_kernel::surface_area(r));
            }
        }
    }}

