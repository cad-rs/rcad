//! Analytic sphere-box union builder.
//!
//! Builds a BRep for the union of a sphere and an axis-aligned box using exact
//! analytic geometry (no tessellation).
//!
//! For each box face that intersects the sphere this produces:
//! - A planar face (outer rectangle, inner circular hole)
//! - A spherical cap face bounded by the intersection circle
//!
//! Box faces that do not intersect the sphere get a full planar face (no hole).

use glam::DVec3;
use rcad_kernel::geom::{any_perpendicular, Circle3, Curve3, Line3, Plane, SphericalSurface, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::BRep;
use rcad_modeling::builder::brep_builder::{make_edge, make_face, make_vertex, make_wire};

// ── Helpers ───────────────────────────────────────────────────────────

/// Detect sphere center and radius from a BRep by scanning all face surfaces.
fn sphere_center_r(sphere: &BRep) -> Option<(DVec3, f64)> {
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

/// Compute the axis-aligned bounding-box min/max from a BRep's vertices.
/// Returns `None` if the vertex set is degenerate.
fn compute_bbox_min_max(brep: &BRep) -> Option<(DVec3, DVec3)> {
    let mut bmin = DVec3::splat(f64::MAX);
    let mut bmax = DVec3::splat(f64::MIN);
    for v in &brep.vertices {
        bmin = bmin.min(v.point);
        bmax = bmax.max(v.point);
    }
    if bmin.x < bmax.x && bmin.y < bmax.y && bmin.z < bmax.z {
        Some((bmin, bmax))
    } else {
        None
    }
}

/// Ensure all parallel `GeomStore` edge vectors are long enough for `edge_idx`.
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

// ── Public API ────────────────────────────────────────────────────────

/// Build the analytic union of a sphere and an axis-aligned box.
///
/// The box is assumed to be axis-aligned (extents discovered from vertex
/// bounding-box).  Returns `None` when either operand cannot be identified.
pub fn build_sphere_box_union_analytic(sphere: &BRep, box_: &BRep) -> Option<BRep> {
    let (center, radius) = sphere_center_r(sphere)?;
    let (bmin, bmax) = compute_bbox_min_max(box_)?;

    let mut brep = BRep::new();
    let two_pi = 2.0 * std::f64::consts::PI;

    // ── 1. Box corner vertices (8 unique) ─────────────────────────────
    let corners: [DVec3; 8] = [
        DVec3::new(bmin.x, bmin.y, bmin.z), // 0
        DVec3::new(bmax.x, bmin.y, bmin.z), // 1
        DVec3::new(bmax.x, bmax.y, bmin.z), // 2
        DVec3::new(bmin.x, bmax.y, bmin.z), // 3
        DVec3::new(bmin.x, bmin.y, bmax.z), // 4
        DVec3::new(bmax.x, bmin.y, bmax.z), // 5
        DVec3::new(bmax.x, bmax.y, bmax.z), // 6
        DVec3::new(bmin.x, bmax.y, bmax.z), // 7
    ];
    let mut cvi = [0usize; 8];
    for (i, &p) in corners.iter().enumerate() {
        cvi[i] = make_vertex(&mut brep, p);
    }

    // ── 2. Box face definitions ──────────────────────────────────────
    // Each entry: (indices into `corners[]` in CCW order from outside,
    //              outward normal, a point on the plane).
    //
    //    Face     CCW corner order     Outward normal   Plane point
    //    ────────────────────────────────────────────────────────────
    //    -Z (bot)  0→3→2→1             ( 0, 0,-1)       (0,0,bmin.z)
    //    +Z (top)  4→5→6→7             ( 0, 0, 1)       (0,0,bmax.z)
    //    -Y (fwd)  0→1→5→4             ( 0,-1, 0)       (0,bmin.y,0)
    //    +Y (bck)  3→7→6→2             ( 0, 1, 0)       (0,bmax.y,0)
    //    -X (lft)  0→4→7→3             (-1, 0, 0)       (bmin.x,0,0)
    //    +X (rgt)  1→2→6→5             ( 1, 0, 0)       (bmax.x,0,0)

    struct FaceInfo {
        corners: [usize; 4],
        normal: DVec3,
        plane_origin: DVec3,
    }

    let faces = [
        FaceInfo { corners: [0, 3, 2, 1], normal: DVec3::NEG_Z, plane_origin: DVec3::new(0.0, 0.0, bmin.z) },
        FaceInfo { corners: [4, 5, 6, 7], normal: DVec3::Z,     plane_origin: DVec3::new(0.0, 0.0, bmax.z) },
        FaceInfo { corners: [0, 1, 5, 4], normal: DVec3::NEG_Y, plane_origin: DVec3::new(0.0, bmin.y, 0.0) },
        FaceInfo { corners: [3, 7, 6, 2], normal: DVec3::Y,     plane_origin: DVec3::new(0.0, bmax.y, 0.0) },
        FaceInfo { corners: [0, 4, 7, 3], normal: DVec3::NEG_X, plane_origin: DVec3::new(bmin.x, 0.0, 0.0) },
        FaceInfo { corners: [1, 2, 6, 5], normal: DVec3::X,     plane_origin: DVec3::new(bmax.x, 0.0, 0.0) },
    ];

    // ── 3. Shared box edges (12 unique) ──────────────────────────────
    // Pre-create all box edges once and share them between adjacent faces.
    let box_edge_pairs: [(usize, usize); 12] = [
        (0, 1), (0, 3), (0, 4),
        (1, 2), (1, 5),
        (2, 3), (2, 6),
        (3, 7),
        (4, 5), (4, 7),
        (5, 6),
        (6, 7),
    ];
    let mut edge_map: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for &(a, b) in &box_edge_pairs {
        let p0 = corners[a];
        let p1 = corners[b];
        let curve = Curve3::Line(Line3 { origin: p0, direction: p1 - p0 });
        let ei = make_edge(&mut brep, curve, 0.0, 1.0, cvi[a], cvi[b]).ok()?;
        align_edge_geom(&mut brep, ei);
        edge_map.insert((a, b), ei);
    }

    // ── 4. Build faces ───────────────────────────────────────────────
    for fi in &faces {
        let n = fi.normal;
        let pp = fi.plane_origin;

        // Signed distance from sphere centre to plane (positive = outward side).
        let d = n.dot(center - pp);

        // ── 3a. Outer wire: 4 shared line edges around the rectangle ──
        let c = fi.corners;
        let mut wire_edges = Vec::with_capacity(4);
        for i in 0..4 {
            let j = (i + 1) % 4;
            let a = c[i];
            let b = c[j];
            let (ea, eb) = if a < b { (a, b) } else { (b, a) };
            let ei = edge_map[&(ea, eb)];
            wire_edges.push(if a < b {
                WireEdge::fwd(ei)
            } else {
                WireEdge::rev(ei)
            });
        }
        let outer_wire = make_wire(wire_edges);

        if d.abs() < radius {
            // ── 3b. Sphere intersects this plane: hole + cap ──
            let circle_center = center - n * d;
            let circle_r = (radius * radius - d * d).sqrt();

            // Build a local orthonormal frame inside the plane.
            let x_axis = any_perpendicular(n);

            // Single vertex on the circle at t = 0.
            let circle_v0 = circle_center + circle_r * x_axis;
            let cv = make_vertex(&mut brep, circle_v0);

            // Circle curve and edge (closed: start == end).
            let circle_curve = Curve3::Circle(Circle3 {
                center: circle_center,
                normal: n,
                radius: circle_r,
            });
            let circle_e = make_edge(&mut brep, circle_curve, 0.0, two_pi, cv, cv).ok()?;
            align_edge_geom(&mut brep, circle_e);

            // Inner wire (hole direction is opposite to outer-wire direction).
            let inner_wire = make_wire(vec![WireEdge::rev(circle_e)]);

            // Planar face — box face with a circular hole.
            let plane_surf = Surface3::Plane(Plane { origin: pp, normal: n });
            make_face(&mut brep, plane_surf, outer_wire, vec![inner_wire]).ok()?;

            // Spherical cap face — the portion of the sphere protruding
            // outside the box through this face.
            let sphere_surf = Surface3::Sphere(SphericalSurface {
                center,
                axis: n,
                radius,
                ref_dir: x_axis,
            });
            let cap_wire = make_wire(vec![WireEdge::fwd(circle_e)]);
            make_face(&mut brep, sphere_surf, cap_wire, vec![]).ok()?;
        } else {
            // ── 3c. No intersection: full planar face (no hole) ──
            let plane_surf = Surface3::Plane(Plane { origin: pp, normal: n });
            make_face(&mut brep, plane_surf, outer_wire, vec![]).ok()?;
        }
    }

    Some(brep)
}
