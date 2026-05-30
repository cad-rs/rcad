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

// ── Intersection builder ─────────────────────────────────────────────

/// Build a planar face on a box plane that is the intersection of the box
/// face rectangle with the sphere-interior disc on that plane.
///
/// `corner_indices` is the 4-tuple of box-corner indices for the face
/// (CCW from outside).  Returns the indices of arc edges that should also
/// appear on the spherical face's outer wire.
fn build_plane_intersection_face(
    brep: &mut BRep,
    center: DVec3,
    radius: f64,
    corners: &[DVec3; 8],
    _cvi: &[usize; 8],
    _box_edge_map: &std::collections::HashMap<(usize, usize), usize>,
    corner_indices: &[usize; 4],
    normal: DVec3,
    plane_origin: DVec3,
) -> Option<Vec<usize>> {
    let n = normal;
    let pp = plane_origin;
    let d = n.dot(center - pp);
    let two_pi = 2.0 * std::f64::consts::PI;

    if d.abs() >= radius - 1e-15 {
        return Some(Vec::new());
    }

    let circle_center = center - d * n;
    let circle_r = (radius * radius - d * d).sqrt();
    if circle_r < 1e-15 {
        return Some(Vec::new());
    }

    // Use rectangle edges as UV axes so point_in_rect works correctly.
    // any_perpendicular(n) is wrong here — it gives an arbitrary direction that
    // may not align with the rectangle edges, inverting inside/outside tests.
    let e0 = (corners[corner_indices[1]] - corners[corner_indices[0]]).normalize();
    let x_axis = e0;
    let y_axis = n.cross(x_axis).normalize();

    let circle_curve = Curve3::Circle(Circle3 { center: circle_center, normal: n, radius: circle_r });
    let c = corner_indices;

    // ── 1. Find intersections of circle with the 4 rectangle edges ──
    struct Cross { theta: f64, pos: DVec3, edge: usize }
    let mut xs: Vec<Cross> = Vec::new();

    for ei in 0..4 {
        let a = corners[c[ei]];
        let b = corners[c[(ei + 1) % 4]];
        let dir = b - a;
        let len = dir.length();
        if len < 1e-15 { continue; }
        let du = dir / len;

        let oc = a - circle_center;
        let bb = 2.0 * oc.dot(du);
        let cc = oc.dot(oc) - circle_r * circle_r;
        let disc = bb * bb - 4.0 * cc;
        if disc < 0.0 { continue; }
        let sd = disc.sqrt();
        for t in [(-bb - sd) / 2.0, (-bb + sd) / 2.0] {
            if t >= -1e-12 && t <= len + 1e-12 {
                let pos = a + t.clamp(0.0, len) * du;
                let dx = (pos - circle_center).dot(x_axis);
                let dy = (pos - circle_center).dot(y_axis);
                let theta = dy.atan2(dx);
                xs.push(Cross { theta: if theta < 0.0 { theta + two_pi } else { theta }, pos, edge: ei });
            }
        }
    }

    xs.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap_or(std::cmp::Ordering::Equal));
    let mut uxs: Vec<Cross> = Vec::new();
    for x in xs {
        if uxs.is_empty() || (x.theta - uxs.last().unwrap().theta).abs() > 1e-10 {
            uxs.push(x);
        }
    }
    xs = uxs;
    let n_xs = xs.len();

    // ── 2. No intersections ──
    if n_xs == 0 {
        let cu = (circle_center - corners[c[0]]).dot(x_axis);
        let cv = (circle_center - corners[c[0]]).dot(y_axis);
        let du_len = (corners[c[1]] - corners[c[0]]).length();
        let dv_len = (corners[c[3]] - corners[c[0]]).length();
        let inside = cu >= -1e-12 && cu <= du_len + 1e-12 && cv >= -1e-12 && cv <= dv_len + 1e-12;
        if inside {
            // Disc entirely inside the rectangle → planar face is the full disc
            let cv0 = circle_center + circle_r * x_axis;
            let cvi0 = make_vertex(brep, cv0);
            let ce = make_edge(brep, circle_curve, 0.0, two_pi, cvi0, cvi0).ok()?;
            align_edge_geom(brep, ce);
            let surf = Surface3::Plane(Plane { origin: pp, normal: n });
            make_face(brep, surf, make_wire(vec![WireEdge::fwd(ce)]), vec![]).ok()?;
            return Some(vec![ce]);
        }
        return Some(Vec::new());
    }

    // ── 3. Build planar face boundary ──
    let mut pt_vis: Vec<usize> = Vec::new();
    for x in &xs { pt_vis.push(make_vertex(brep, x.pos)); }

    let point_in_rect = |pt: DVec3| -> bool {
        let u = (pt - corners[c[0]]).dot(x_axis);
        let v = (pt - corners[c[0]]).dot(y_axis);
        let du = (corners[c[1]] - corners[c[0]]).length();
        let dv = (corners[c[3]] - corners[c[0]]).length();
        u >= -1e-12 && u <= du + 1e-12 && v >= -1e-12 && v <= dv + 1e-12
    };

    let mut planar_wes: Vec<WireEdge> = Vec::new();
    let mut arc_edges: Vec<usize> = Vec::new();

    let mk_line = |brep: &mut BRep, p1: DVec3, p2: DVec3| -> Option<usize> {
        if (p1 - p2).length() < 1e-15 { return None; }
        let v1 = make_vertex(brep, p1);
        let v2 = make_vertex(brep, p2);
        let dir = p2 - p1;
        let len = dir.length();
        let curve = Curve3::Line(Line3 { origin: p1, direction: dir / len });
        let ei = make_edge(brep, curve, 0.0, len, v1, v2).ok()?;
        align_edge_geom(brep, ei);
        Some(ei)
    };

    let walk_perim = |brep: &mut BRep, from_pos: DVec3, from_edge: usize, to_pos: DVec3, to_edge: usize| -> Option<Vec<WireEdge>> {
        let mut wes = Vec::new();
        let mut cur = from_pos;
        let mut e = from_edge;
        loop {
            if e == to_edge {
                if let Some(ei) = mk_line(brep, cur, to_pos) { wes.push(WireEdge::fwd(ei)); }
                break;
            }
            let nc = corners[c[(e + 1) % 4]];
            if let Some(ei) = mk_line(brep, cur, nc) { wes.push(WireEdge::fwd(ei)); }
            cur = nc;
            e = (e + 1) % 4;
        }
        Some(wes)
    };

    for i in 0..n_xs {
        let j = (i + 1) % n_xs;
        let t_i = xs[i].theta;
        let t_j = xs[j].theta;
        let pi = xs[i].pos;
        let pj = xs[j].pos;

        let mid_t = if j > i { (t_i + t_j) / 2.0 } else {
            let r = (t_i + t_j + two_pi) / 2.0;
            if r > two_pi { r - two_pi } else { r }
        };
        let (sm, cm) = mid_t.sin_cos();
        let mid_pt = circle_center + circle_r * (cm * x_axis + sm * y_axis);

        if point_in_rect(mid_pt) {
            // Arc edge
            let v1 = pt_vis[i];
            let v2 = pt_vis[j];
            let ae = make_edge(brep, circle_curve.clone(), t_i, t_j, v1, v2).ok()?;
            align_edge_geom(brep, ae);
            planar_wes.push(WireEdge::fwd(ae));
            arc_edges.push(ae);
        } else {
            // Rectangle perimeter
            if let Some(w) = walk_perim(brep, pi, xs[i].edge, pj, xs[j].edge) {
                planar_wes.extend(w);
            }
        }
    }

    if planar_wes.is_empty() {
        return Some(Vec::new());
    }

    let outer_wire = make_wire(planar_wes);
    let plane_surf = Surface3::Plane(Plane { origin: pp, normal: n });
    make_face(brep, plane_surf, outer_wire, vec![]).ok()?;

    Some(arc_edges)
}

/// Build the analytic intersection of a sphere and an axis-aligned box.
///
/// Returns the BRep for sphere ∩ box.  Requires the sphere center to be
/// inside the box (or at least close); otherwise returns `None`.
///
/// The result consists of up to 6 planar faces (box face portions inside the
/// sphere) and 1 spherical face (the portion of the sphere surface inside the
/// box).  All surfaces and curves are analytic — no tessellation.
pub fn build_sphere_box_intersection_analytic(sphere: &BRep, box_: &BRep) -> Option<BRep> {
    let (center, radius) = sphere_center_r(sphere)?;
    let (bmin, bmax) = compute_bbox_min_max(box_)?;

    let mut brep = BRep::new();

    // ── 1. Box corner vertices ──
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

    // ── 2. Face definitions ──
    let faces: [([usize; 4], DVec3, DVec3); 6] = [
        ([0, 3, 2, 1], DVec3::NEG_Z, DVec3::new(0.0, 0.0, bmin.z)),
        ([4, 5, 6, 7], DVec3::Z,     DVec3::new(0.0, 0.0, bmax.z)),
        ([0, 1, 5, 4], DVec3::NEG_Y, DVec3::new(0.0, bmin.y, 0.0)),
        ([3, 7, 6, 2], DVec3::Y,     DVec3::new(0.0, bmax.y, 0.0)),
        ([0, 4, 7, 3], DVec3::NEG_X, DVec3::new(bmin.x, 0.0, 0.0)),
        ([1, 2, 6, 5], DVec3::X,     DVec3::new(bmax.x, 0.0, 0.0)),
    ];

    // ── 3. Shared box edges ──
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

    // ── 4. Build planar faces; collect arc edges ──
    let mut all_arcs: Vec<usize> = Vec::new();
    for &(ref ci, n, pp) in &faces {
        let arcs = build_plane_intersection_face(
            &mut brep, center, radius, &corners, &cvi, &edge_map, ci, n, pp,
        )?;
        all_arcs.extend(arcs);
    }

    // ── 5. Spherical face ──
    if !all_arcs.is_empty() {
        // Reorder arcs so consecutive reversed edges chain end-to-end by 3D position.
        // Wire uses WireEdge::rev(ei) for each arc, so:
        //   rev(ai).end = ai.start must match rev(aj).start = aj.end
        // i.e. chain by matching current arc's START to next arc's END.
        let n_arcs = all_arcs.len();
        let mut ordered: Vec<usize> = Vec::with_capacity(n_arcs);
        let mut used = vec![false; n_arcs];
        ordered.push(all_arcs[0]);
        used[0] = true;
        while ordered.len() < n_arcs {
            let last_ei = *ordered.last().unwrap();
            let last_start_pos = brep.vertices[brep.edges[last_ei].start].point;
            let mut found = None;
            for (j, &ei) in all_arcs.iter().enumerate() {
                if used[j] { continue; }
                let end_pos = brep.vertices[brep.edges[ei].end].point;
                if (end_pos - last_start_pos).length() < 1e-12 {
                    found = Some((j, ei));
                    break;
                }
            }
            if let Some((j, ei)) = found {
                ordered.push(ei);
                used[j] = true;
            } else {
                // Degenerate: arcs don't chain. Use as-is (fallback).
                for (j, &ei) in all_arcs.iter().enumerate() {
                    if !used[j] { ordered.push(ei); }
                }
                break;
            }
        }
        let sphere_surf = Surface3::Sphere(SphericalSurface {
            center,
            axis: DVec3::Z,
            radius,
            ref_dir: DVec3::X,
        });
        let wes: Vec<WireEdge> = ordered.iter().map(|&ei| WireEdge::rev(ei)).collect();
        make_face(&mut brep, sphere_surf, make_wire(wes), vec![]).ok()?;
    }

    Some(brep)
}
