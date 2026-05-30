//! Analytic cylinder-box union builder.
//!
//! Builds a BRep for the union of a Z-aligned cylinder and an axis-aligned box
//! using exact analytic geometry (no tessellation / Pave-Filler).
//!
//! The result contains:
//!
//! **Part A** — Cylinder wall arcs (and caps) where the cylinder protrudes
//! outside the box in XY.  Built by computing the complement of the inside-box
//! theta intervals and delegating to `build_cylinder_box_clipped_brep`.
//!
//! **Part B** — Box faces trimmed away from the cylinder.  Each box face that
//! the cylinder passes through gets a planar face with a cylindrical cutout
//! (inner wire).  Z-faces get a circular hole; side faces get a rectangular
//! slot (the extrusion of the cylinder cross-section in Z).

use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::geom::{any_perpendicular, Circle3, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::BRep;
use rcad_modeling::builder::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use crate::tolerance::*;

// ── Public API ──────────────────────────────────────────────────────────────

/// Build the analytic union of a Z-aligned cylinder and an axis-aligned box.
///
/// Returns `None` when either operand cannot be identified or the configuration
/// is not handled by this fast path.
///
/// The result is a single BRep containing:
/// 1. Cylinder wall arcs (plus top/bottom caps) for the cylinder portion
///    outside the box XY projection.
/// 2. Six box faces — each potentially carrying a cutout where the cylinder
///    passes through the face.
pub fn build_cylinder_box_union_analytic(cyl_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    // ── 1. Detect operands ──────────────────────────────────────────────
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) =
        super::boolean_unit_octant::try_cylinder_center_axis_radius_height(cyl_brep)?;
    let cyl_center = cyl_bottom + cyl_axis * (cyl_height / 2.0);
    let bx = super::boolean_unit_octant::try_as_box(box_brep)?;

    // Only Z-aligned cylinders for now.
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let z_idx = super::boolean_unit_octant::find_z_axis_index(&bx)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = bx.axes[u_idx];
    let v_ax = bx.axes[v_idx];
    let _z_ax = bx.axes[z_idx];
    let eu = bx.extents[u_idx];
    let ev = bx.extents[v_idx];
    let ew = bx.extents[z_idx];
    let bc = bx.center;

    let cyl_z_lo = cyl_center.z - cyl_height / 2.0;
    let cyl_z_hi = cyl_center.z + cyl_height / 2.0;
    let box_z_lo = bc.z - ew;
    let box_z_hi = bc.z + ew;

    let tol = TOLERANCE_LEN_MIN;

    // ── 2. Compute clip planes (box interior constraints in XY) ────────
    let cu = (cyl_center - bc).dot(u_ax);
    let cv = (cyl_center - bc).dot(v_ax);

    let mut clip_planes: Vec<(DVec3, f64)> = Vec::new();
    for &(cp, ax, ext) in &[(cu, u_ax, eu), (cv, v_ax, ev)] {
        if cp - cyl_r < -ext - tol {
            clip_planes.push((ax, cp + ext));
        }
        if cp + cyl_r > ext + tol {
            clip_planes.push((-ax, ext - cp));
        }
    }

    if clip_planes.is_empty() {
        // Cylinder is fully inside the box XY cross-section.
        // The union is the box — nothing from the cylinder protrudes.
        return None;
    }

    // ── 3. Theta intervals where the cylinder wall is OUTSIDE the box ──
    let inside_intervals =
        super::boolean_unit_octant::compute_valid_theta_ranges(cyl_r, &clip_planes);
    let complement =
        super::boolean_unit_octant::compute_complement_theta_ranges(&inside_intervals);

    if complement.is_empty() {
        return None;
    }

    // ── 4. Build Part A: cylinder wall + caps ──────────────────────────
    //
    // We build up to three Z-segments and combine them:
    //
    //   (a) Lower segment:  cyl_z_lo  → box_z_lo   (full wall if non-empty)
    //   (b) Middle segment: box_z_lo  → box_z_hi   (complement intervals)
    //   (c) Upper segment:  box_z_hi  → cyl_z_hi   (full wall if non-empty)
    //
    // Segments (a) and (c) are full cylinders; segment (b) is the
    // cylinder wall trimmed to the complement theta intervals.

    // --- (b) Middle segment: trimmed cylinder at box Z range ---
    let mut result = BRep::new();
    let inter_lo = cyl_z_lo.max(box_z_lo);
    let inter_hi = cyl_z_hi.min(box_z_hi);
    if inter_hi > inter_lo + tol {
        let mid_h = inter_hi - inter_lo;
        let mid_cz = inter_lo + mid_h / 2.0;
        let mid_center = DVec3::new(cyl_center.x, cyl_center.y, mid_cz);

        let mid_part = super::boolean_unit_octant::build_cylinder_box_clipped_brep(
            mid_center,
            cyl_r,
            mid_h,
            &complement,
            &clip_planes,
            /*skip_bottom_cap=*/ false,
            /*skip_top_cap=*/ false,
            /*use_chain_routing=*/ false,
        );
        result.append_disjoint_brep(&mid_part);
    }

    // --- (a) Lower segment: full cylinder below box ---
    if cyl_z_lo < box_z_lo - tol {
        let lower_h = box_z_lo - cyl_z_lo;
        let lower_cz = cyl_z_lo + lower_h / 2.0;
        let lower_center = DVec3::new(cyl_center.x, cyl_center.y, lower_cz);
        if let Ok(cyl) = rcad_modeling::make_cylinder_brep(
            lower_center, cyl_axis, u_ax, cyl_r, lower_h,
        ) {
            result.append_disjoint_brep(&cyl);
        }
    }

    // --- (c) Upper segment: full cylinder above box ---
    if cyl_z_hi > box_z_hi + tol {
        let upper_h = cyl_z_hi - box_z_hi;
        let upper_cz = box_z_hi + upper_h / 2.0;
        let upper_center = DVec3::new(cyl_center.x, cyl_center.y, upper_cz);
        if let Ok(cyl) = rcad_modeling::make_cylinder_brep(
            upper_center, cyl_axis, u_ax, cyl_r, upper_h,
        ) {
            result.append_disjoint_brep(&cyl);
        }
    }

    // ── 5. Build Part B: box faces with cylindrical cutouts ────────────
    let box_part = build_box_faces_with_cylinder_cutouts(
        &bc, &bx.axes, &bx.extents,
        cyl_center, cyl_r, cyl_z_lo, cyl_z_hi, inter_lo, inter_hi,
    );
    if let Some(bp) = box_part {
        result.append_disjoint_brep(&bp);
    }

    Some(result)
}

// ── Helpers ─────────────────────────────────────────────────────────────────

/// Ensure `GeomStore` edge vectors are long enough for `edge_idx`.
fn align_edge_geom(brep: &mut BRep, edge_idx: usize) {
    while brep.geom.edge_pcurves.len() <= edge_idx {
        brep.geom.edge_pcurves.push(Vec::new());
    }
    while brep.geom.edge_same_parameter.len() <= edge_idx {
        brep.geom.edge_same_parameter.push(false);
    }
    while brep.geom.edge_same_range.len() <= edge_idx {
        brep.geom.edge_same_range.push(false);
    }
}

/// Build the 6 box faces, each trimmed by the cylinder if the cylinder
/// passes through the face.
fn build_box_faces_with_cylinder_cutouts(
    bc: &DVec3,
    axes: &[DVec3; 3],
    extents: &[f64; 3],
    cyl_center: DVec3,
    cyl_r: f64,
    cyl_z_lo: f64,
    cyl_z_hi: f64,
    inter_z_lo: f64,
    inter_z_hi: f64,
) -> Option<BRep> {
    let [ua, va, _wa] = *axes;
    let [eu, ev, ew] = *extents;
    let c = *bc;

    let z_idx = super::boolean_unit_octant::find_z_axis_index(
        &super::boolean_unit_octant::BoxInfo {
            axes: *axes,
            center: *bc,
            extents: *extents,
        }
    )?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = axes[u_idx];
    let v_ax = axes[v_idx];
    let z_ax = axes[z_idx];

    let tol = TOLERANCE_LEN_MIN;
    let two_pi = 2.0 * std::f64::consts::PI;

    // ── 1. Box corner vertices (8 unique) ──────────────────────────────
    // Indices 0..7 in the BRep vertex table
    let corners: [DVec3; 8] = [
        c - eu*u_ax - ev*v_ax - ew*z_ax,   // 0: (-, -, -)
        c + eu*u_ax - ev*v_ax - ew*z_ax,   // 1: (+, -, -)
        c + eu*u_ax + ev*v_ax - ew*z_ax,   // 2: (+, +, -)
        c - eu*u_ax + ev*v_ax - ew*z_ax,   // 3: (-, +, -)
        c - eu*u_ax - ev*v_ax + ew*z_ax,   // 4: (-, -, +)
        c + eu*u_ax - ev*v_ax + ew*z_ax,   // 5: (+, -, +)
        c + eu*u_ax + ev*v_ax + ew*z_ax,   // 6: (+, +, +)
        c - eu*u_ax + ev*v_ax + ew*z_ax,   // 7: (-, +, +)
    ];

    let mut brep = BRep::new();

    // Push all 8 vertices first (so indices 0..7 correspond to corners).
    for &p in &corners {
        make_vertex(&mut brep, p);
    }

    // ── 2. Shared box edges (12 unique) ────────────────────────────────
    // Each edge has a start/end vertex index into the 8 corner vertices.
    // We store edges keyed by (min_vertex, max_vertex) so adjacent faces
    // can share the same edge.
    let box_edge_pairs: [(usize, usize); 12] = [
        (0, 1), (0, 3), (0, 4),
        (1, 2), (1, 5),
        (2, 3), (2, 6),
        (3, 7),
        (4, 5), (4, 7),
        (5, 6),
        (6, 7),
    ];
    let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
    for &(a, b) in &box_edge_pairs {
        let p0 = corners[a];
        let p1 = corners[b];
        let curve = Curve3::Line(Line3 { origin: p0, direction: p1 - p0 });
        let ei = make_edge(&mut brep, curve, 0.0, 1.0, a, b).ok()?;
        align_edge_geom(&mut brep, ei);
        edge_map.insert((a, b), ei);
    }

    // ── 3. Helper: get a box edge, handling edge direction ─────────────
    let get_edge = |a: usize, b: usize| -> usize {
        let (ea, eb) = if a < b { (a, b) } else { (b, a) };
        edge_map[&(ea, eb)]
    };

    // ── 4. Helper: push a planar face (with optional inner wires) ──────
    let push_planar_face = |brep: &mut BRep,
                            outer_wes: Vec<WireEdge>,
                            inner_wires: Vec<Vec<WireEdge>>,
                            normal: DVec3,
                            plane_origin: DVec3| -> Option<()> {
        let outer = make_wire(outer_wes);
        let inners: Vec<rcad_kernel::topology::Wire> =
            inner_wires.into_iter().map(make_wire).collect();

        let surf = Surface3::Plane(Plane { origin: plane_origin, normal });
        make_face(brep, surf, outer, inners).ok()?;
        Some(())
    };

    // Cylinder center projected onto box axes.
    let c_cs = (cyl_center - c).dot(u_ax); // cylinder center in u
    let c_ct = (cyl_center - c).dot(v_ax); // cylinder center in v
    let c_cw = (cyl_center - c).dot(z_ax); // cylinder center in w

    let cyl_in_box_z = cyl_z_hi > c_cw - ew - tol && cyl_z_lo < c_cw + ew + tol;

    // ── 5a. -Z face (bottom) & +Z face (top) ──────────────────────────
    //   -Z (bot): corners 0→3→2→1, normal = -z_ax
    //   +Z (top): corners 4→5→6→7, normal = +z_ax
    for &(corner_inds, normal, w_val, is_top) in &[
        ([0, 3, 2, 1], -z_ax, -ew, false),
        ([4, 5, 6, 7],  z_ax,  ew, true),
    ] {
        let ci = corner_inds;
        let mut outer_wes = Vec::with_capacity(4);
        for i in 0..4 {
            let j = (i + 1) % 4;
            let a = ci[i];
            let b = ci[j];
            let ei = get_edge(a, b);
            outer_wes.push(if a < b { WireEdge::fwd(ei) } else { WireEdge::rev(ei) });
        }

        let face_z = c_cw + w_val;

        // Does the cylinder pass through this Z-face?
        let cyl_passes = cyl_in_box_z
            && ((is_top && cyl_z_hi > face_z + tol)
                || (!is_top && cyl_z_lo < face_z - tol))
            && (c_cs.abs() - cyl_r < eu + tol || c_ct.abs() - cyl_r < ev + tol);

        if !cyl_passes {
            push_planar_face(&mut brep, outer_wes, vec![], normal, DVec3::new(0.0, 0.0, face_z))?;
            continue;
        }

        // Z-face with circular hole: the cylinder cross-section circle.
        let circle_center = DVec3::new(cyl_center.x, cyl_center.y, face_z);
        let x_axis = any_perpendicular(normal);

        let circle_v0 = circle_center + cyl_r * x_axis;
        let cv = make_vertex(&mut brep, circle_v0);

        let circle_curve = Curve3::Circle(Circle3 {
            center: circle_center,
            normal,
            radius: cyl_r,
        });
        let circle_e = make_edge(&mut brep, circle_curve, 0.0, two_pi, cv, cv).ok()?;
        align_edge_geom(&mut brep, circle_e);

        // Inner wire direction is reversed from outer wire (hole winding).
        let inner_wire = vec![WireEdge::rev(circle_e)];

        push_planar_face(&mut brep, outer_wes, vec![inner_wire], normal, DVec3::new(0.0, 0.0, face_z))?;
    }

    // ── 5b. Side faces (-v, +v, -u, +u) ────────────────────────────────
    // Each side face may have a rectangular slot where the cylinder
    // protrudes through the face.  The slot is bounded by:
    //   left/right: cylinder silhouette on the face
    //   bottom/top: overlap Z range
    //
    // Side-face definitions:
    //
    //   -v (fwd): 0→1→5→4, normal = -v_ax (perp axis = u_ax)
    //   +v (bck): 3→7→6→2, normal = +v_ax (perp axis = u_ax)
    //   -u (lft): 0→4→7→3, normal = -u_ax (perp axis = v_ax)
    //   +u (rgt): 1→2→6→5, normal = +u_ax (perp axis = v_ax)
    let side_faces: [([usize; 4], DVec3, DVec3); 4] = [
        ([0, 1, 5, 4], -v_ax, u_ax),  // -v face, perp = u
        ([3, 7, 6, 2],  v_ax, u_ax),  // +v face, perp = u
        ([0, 4, 7, 3], -u_ax, v_ax),  // -u face, perp = v
        ([1, 2, 6, 5],  u_ax, v_ax),  // +u face, perp = v
    ];

    for &(ci, normal, perp_ax) in &side_faces {
        // Build outer wire from box edges.
        let mut outer_wes = Vec::with_capacity(4);
        for i in 0..4 {
            let j = (i + 1) % 4;
            let a = ci[i];
            let b = ci[j];
            let ei = get_edge(a, b);
            outer_wes.push(if a < b { WireEdge::fwd(ei) } else { WireEdge::rev(ei) });
        }

        // Determine which coordinate axis the face normal is aligned with.
        // perp_ax is the tangential axis parallel to the face (horizontal).
        let is_u_face = normal.dot(ua).abs() > 0.9;
        let is_v_face = normal.dot(va).abs() > 0.9;

        // Cylinder center projection along the face-normal direction.
        let center_on_normal = if is_u_face { c_cs } else if is_v_face { c_ct } else { 0.0 };

        // Half-extent along the normal direction.
        let half_ext = if is_u_face { eu } else if is_v_face { ev } else { 0.0 };

        // Center projection along the perpendicular (tangential) axis.
        let perp_center = match perp_ax == u_ax {
            true => c_cs,
            false => c_ct,
        };

        // Half-extent along the perpendicular axis.
        let perp_half = match perp_ax == u_ax { true => eu, false => ev };

        // Does the cylinder protrude through this face?
        // The cylinder center must be close enough that the cylinder radius
        // extends outward past the face, OR the center is already beyond the face.
        let protrudes = cyl_in_box_z
            && inter_z_hi > inter_z_lo + tol
            && (center_on_normal + cyl_r > half_ext + tol
                || center_on_normal - cyl_r < -half_ext - tol)
            && perp_center.abs() < perp_half + cyl_r + tol;

        if !protrudes {
            push_planar_face(&mut brep, outer_wes, vec![], normal, DVec3::new(0.0, 0.0, 0.0))?;
            continue;
        }

        // ── Build the rectangular slot (inner wire) ──
        // The slot is the cylinder's projection onto this face.
        // Its horizontal bounds are the cylinder radius range projected onto
        // the perpendicular axis.  Its vertical bounds are the overlap Z range.

        let perp_min = (perp_center - cyl_r).max(-perp_half);
        let perp_max = (perp_center + cyl_r).min(perp_half);

        if perp_max <= perp_min + tol || inter_z_hi <= inter_z_lo + tol {
            push_planar_face(&mut brep, outer_wes, vec![], normal, DVec3::new(0.0, 0.0, 0.0))?;
            continue;
        }

        // Slot corners in (u, v, w) coordinates as tangent × normal × Z on the face.
        // The face has constant u or v (the normal direction), and extends in
        // the perpendicular axis and Z.
        let slot_pts_uvw: [(f64, f64, f64); 4] = if is_u_face {
            // Face at constant u = ±eu.  Tangent axes: (v, w).
            // v goes from perp_min to perp_max, w goes from inter_z_lo to inter_z_hi.
            let u_const = half_ext * normal.dot(ua).signum();
            [
                (u_const, perp_min, inter_z_lo),
                (u_const, perp_max, inter_z_lo),
                (u_const, perp_max, inter_z_hi),
                (u_const, perp_min, inter_z_hi),
            ]
        } else {
            // Face at constant v = ±ev.  Tangent axes: (u, w).
            let v_const = half_ext * normal.dot(va).signum();
            [
                (perp_min, v_const, inter_z_lo),
                (perp_max, v_const, inter_z_lo),
                (perp_max, v_const, inter_z_hi),
                (perp_min, v_const, inter_z_hi),
            ]
        };

        let to_world = |(u, v, w): (f64, f64, f64)| -> DVec3 {
            c + u*u_ax + v*v_ax + w*z_ax
        };

        let slot_world: [DVec3; 4] = [
            to_world(slot_pts_uvw[0]),
            to_world(slot_pts_uvw[1]),
            to_world(slot_pts_uvw[2]),
            to_world(slot_pts_uvw[3]),
        ];

        let sv: [usize; 4] = [
            make_vertex(&mut brep, slot_world[0]),
            make_vertex(&mut brep, slot_world[1]),
            make_vertex(&mut brep, slot_world[2]),
            make_vertex(&mut brep, slot_world[3]),
        ];

        let mut slot_edges: Vec<usize> = Vec::with_capacity(4);
        for i in 0..4 {
            let j = (i + 1) % 4;
            let dir = slot_world[j] - slot_world[i];
            let len = dir.length();
            if len < tol {
                continue;
            }
            let curve = Curve3::Line(Line3 { origin: slot_world[i], direction: dir / len });
            let ei = make_edge(&mut brep, curve, 0.0, len, sv[i], sv[j]).ok()?;
            align_edge_geom(&mut brep, ei);
            slot_edges.push(ei);
        }

        if slot_edges.len() < 4 {
            // Degenerate slot → full face.
            push_planar_face(&mut brep, outer_wes, vec![], normal, DVec3::new(0.0, 0.0, 0.0))?;
            continue;
        }

        // Inner wire with reversed winding (CW instead of CCW) for hole convention.
        let inner_wes: Vec<WireEdge> = slot_edges.iter().rev().map(|&ei| WireEdge::fwd(ei)).collect();
        push_planar_face(&mut brep, outer_wes, vec![inner_wes], normal, DVec3::new(0.0, 0.0, 0.0))?;
    }

    Some(brep)
}
