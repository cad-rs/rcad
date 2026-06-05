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

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle2d, Circle3, Curve3, Line2d, Line3, Plane, SphericalSurface, Surface3, Curve2d};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::{BRep, PCurve};
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

/// Get or create a shared box edge (and its corner vertices) lazily.
/// If the edge already exists in `box_edge_map`, returns it.  Otherwise creates
/// the required corner vertices (if not yet created), builds the edge, inserts
/// it into the map, and returns the new edge index.
fn get_or_create_box_edge(
    brep: &mut BRep,
    box_edge_map: &mut std::collections::HashMap<(usize, usize), usize>,
    cvi: &mut [Option<usize>; 8],
    corners: &[DVec3; 8],
    a: usize,
    b: usize,
) -> Option<usize> {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&ei) = box_edge_map.get(&key) {
        return Some(ei);
    }
    // Create corner vertices if they don't exist yet.
    if cvi[a].is_none() {
        cvi[a] = Some(make_vertex(brep, corners[a]));
    }
    if cvi[b].is_none() {
        cvi[b] = Some(make_vertex(brep, corners[b]));
    }
    let p0 = corners[a];
    let p1 = corners[b];
    let curve = Curve3::Line(Line3 {
        origin: p0,
        direction: p1 - p0,
    });
    let ei = make_edge(brep, curve, 0.0, 1.0, cvi[a].unwrap(), cvi[b].unwrap()).ok()?;
    align_edge_geom(brep, ei);
    box_edge_map.insert(key, ei);
    Some(ei)
}

/// Compute the UV coordinate of a point on a sphere surface, propagating the
/// longitude (U) from `other_point` when `point` is at the pole (|V| ≈ π/2).
/// This ensures iso-parametric curves (meridians) get a consistent U value
/// even at the degenerate polar vertex.
fn sphere_uv_propagate(sphere: &SphericalSurface, point: DVec3, other_point: DVec3, radius: f64) -> DVec2 {
    let dp = point - sphere.center;
    let dv = other_point - sphere.center;
    let v = (dp.z / radius).asin().clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
    let u = if v.abs() >= std::f64::consts::FRAC_PI_2 - 1e-12 {
        dv.y.atan2(dv.x)
    } else {
        dp.y.atan2(dp.x)
    };
    DVec2::new(u, v)
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

    // ── 1. Box corner positions (vertices created lazily) ──
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
    let mut cvi: [Option<usize>; 8] = [None; 8]; // created lazily

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

    // ── 3. Shared box edges (created lazily) ─────────────────────────
    let mut edge_map: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();

    // ── 4. Build faces ───────────────────────────────────────────────
    // For each box face: if the sphere intersects the plane, build a trimmed
    // planar face via build_plane_intersection_face (outside=true, returning
    // arc edges).  Otherwise, emit a full rectangular planar face.
    let mut all_arcs: Vec<usize> = Vec::new();
    for fi in &faces {
        let n = fi.normal;
        let pp = fi.plane_origin;
        let d = n.dot(center - pp);

        if d.abs() < radius {
            let arcs = build_plane_intersection_face(
                &mut brep, center, radius, &corners, &mut cvi, &mut edge_map,
                &fi.corners, n, pp, true,
            )?;
            all_arcs.extend(arcs);
        } else {
            // Full planar face (sphere does not intersect this plane at all).
            let c = fi.corners;
            let mut wire_edges = Vec::with_capacity(4);
            for i in 0..4 {
                let j = (i + 1) % 4;
                let a = c[i];
                let b = c[j];
                let ei = get_or_create_box_edge(&mut brep, &mut edge_map, &mut cvi, &corners, a, b)?;
                wire_edges.push(if a < b {
                    WireEdge::fwd(ei)
                } else {
                    WireEdge::rev(ei)
                });
            }
            let plane_surf = Surface3::Plane(Plane { origin: pp, normal: n });
            make_face(&mut brep, plane_surf, make_wire(wire_edges), vec![]).ok()?;
        }
    }

    // ── 5. Single merged spherical face with seam ───────────────────
    // Build a single spherical face from all arc edges, adding a SEAM_CURVE
    // from a seam-vertex (u≈0) to the south pole so the face covers the
    // full 7/8 sphere region instead of just the 1/8 octant.
    if !all_arcs.is_empty() {
        let n_arcs = all_arcs.len();
        // Find the 3 distinct vertices of the spherical triangle.
        let mut tri_vis: Vec<usize> = Vec::with_capacity(3);
        for &ei in &all_arcs {
            for &vi in &[brep.edges[ei].start, brep.edges[ei].end] {
                if !tri_vis.contains(&vi) { tri_vis.push(vi); }
            }
        }
        // Pick the vertex whose sphere UV has u ≈ 0 (parameterization seam)
        // as the attachment point for the SEAM_CURVE — matches OCCT.
        let mut seam_vi = tri_vis[0];
        for &vi in &tri_vis {
            let pt = brep.vertices[vi].point - center;
            let u = pt.y.atan2(pt.x);
            if u.abs() < 1e-6 || (u - two_pi).abs() < 1e-6 {
                seam_vi = vi;
                break;
            }
        }
        // Build the arc loop starting & ending at seam_vi so consecutive
        // WireEdges share vertices (necessary for a valid EDGE_LOOP).
        // The traversal direction matches OCCT's reference:
        //   (u≈0, equator) → (0,0,z>0) → (0,y>0,0) → (u≈0, equator)
        // Then the SEAM goes (u≈0, equator) → south pole → back.
        let mut sphere_wes: Vec<WireEdge> = Vec::with_capacity(n_arcs + 2);
        let mut cur_vi = seam_vi;
        let mut remaining: Vec<usize> = all_arcs.clone();
        while !remaining.is_empty() {
            // From seam_vi (on the sphere's parameterization seam at u≈0,
            // z≈0 = equator), prefer the edge going NORTH (z≠0) over the
            // one staying on the equator.  This matches OCCT's loop order:
            //   (u≈0, equator) → (0,0,z>0) → (0,y>0,0) → (u≈0, equator)
            // Without this preference the loop goes equator-first, and the
            // left-hand rule selects the 1/8 octant instead of the 7/8 sphere.
            let pos = if cur_vi == seam_vi && remaining.len() > 1 {
                // Try to find an edge whose other endpoint is off-equator.
                let north_idx = remaining.iter().position(|&ei| {
                    let e = &brep.edges[ei];
                    let other_vi = if e.start == cur_vi { e.end }
                        else if e.end == cur_vi { e.start } else { return false; };
                    (brep.vertices[other_vi].point - center).z.abs() > 1e-6
                });
                // Fall back to the first matching edge if no north-going edge.
                north_idx.or_else(|| {
                    remaining.iter().position(|&ei| {
                        brep.edges[ei].start == cur_vi || brep.edges[ei].end == cur_vi
                    })
                })
            } else {
                remaining.iter().position(|&ei| {
                    brep.edges[ei].start == cur_vi || brep.edges[ei].end == cur_vi
                })
            };
            match pos {
                Some(idx) => {
                    let ei = remaining.remove(idx);
                    let e = &brep.edges[ei];
                    if e.start == cur_vi {
                        sphere_wes.push(WireEdge::fwd(ei));
                        cur_vi = e.end;
                    } else {
                        sphere_wes.push(WireEdge::rev(ei));
                        cur_vi = e.start;
                    }
                }
                None => { break; }
            }
        }
        let loop_closed = cur_vi == seam_vi;
        // Create south-pole vertex and the SEAM meridian edge.
        let seam_ei_opt: Option<usize> = if n_arcs == 3 && loop_closed {
            let south_pole_pt = center - radius * DVec3::Z;
            let south_pole_vi = make_vertex(&mut brep, south_pole_pt);
            // Meridian in the y=0 plane: P(t)=radius*(cos(t), 0, -sin(t))
            // t=0 → (radius,0,0) = seam_vi point; t=π/2 → (0,0,-radius)
            let seam_circle = Curve3::Circle(Circle3 {
                center, normal: DVec3::Y, radius,
            });
            let half_pi = std::f64::consts::FRAC_PI_2;
            let seam_ei =
                make_edge(&mut brep, seam_circle, 0.0, half_pi, seam_vi, south_pole_vi).ok()?;
            align_edge_geom(&mut brep, seam_ei);
            sphere_wes.push(WireEdge::fwd(seam_ei)); // seam_vi → south_pole
            sphere_wes.push(WireEdge::rev(seam_ei)); // south_pole → seam_vi
            Some(seam_ei)
        } else {
            None
        };
        let sphere_surf = Surface3::Sphere(SphericalSurface {
            center,
            axis: DVec3::Z,
            radius,
            ref_dir: DVec3::X,
        });
        make_face(&mut brep, sphere_surf.clone(), make_wire(sphere_wes), vec![]).ok()?;
        // Add spherical-surface pcurves to each arc edge.  The pcurve must
        // follow the 3D curve direction in the STEP file.  For CIRCLE edges,
        // the STEP writer determines same_sense based on the angle between
        // start and end in the CIRCLE's canonical frame.  When same_sense=.F.
        // (curve goes opposite to BRep edge direction), the pcurve direction
        // must be inverted.
        let sphere_surface_idx = brep.geom.surfaces.len() - 1;
        if let Surface3::Sphere(sphere) = &brep.geom.surfaces[sphere_surface_idx] {
            for &ei in &all_arcs {
                let start_pt = brep.vertices[brep.edges[ei].start].point;
                let end_pt = brep.vertices[brep.edges[ei].end].point;
                let uv_start = sphere_uv_propagate(sphere, start_pt, end_pt, radius);
                let uv_end = sphere_uv_propagate(sphere, end_pt, start_pt, radius);
                // Determine if the STEP writer will use same_sense=.F. for this
                // CIRCLE edge by checking the angle in the circle's frame.
                // This mirrors the writer's logic at writer.rs:1847-1888.
                let mut same_sense = true;
                if let Some(curve_idx) = brep.geom.edge_curve.get(ei).copied().flatten() {
                    if let Some(Curve3::Circle(circ)) = brep.geom.curves.get(curve_idx) {
                        let canon_axis = {
                            let n = circ.normal.normalize_or_zero();
                            let eps = 1e-12;
                            if n.z.abs() > eps { if n.z >= 0.0 { n } else { -n } }
                            else if n.y.abs() > eps { if n.y >= 0.0 { n } else { -n } }
                            else { if n.x >= 0.0 { n } else { -n } }
                        };
                        let ref_dir = if canon_axis.z.abs() > 0.999999 {
                            DVec3::X
                        } else {
                            let helper = if canon_axis.y.abs() < 0.9 { DVec3::Y } else { DVec3::X };
                            canon_axis.cross(helper).normalize()
                        };
                        let perp = canon_axis.cross(ref_dir);
                        let d_start = start_pt - circ.center;
                        let d_end = end_pt - circ.center;
                        let theta_start = d_start.dot(perp).atan2(d_start.dot(ref_dir));
                        let theta_end = d_end.dot(perp).atan2(d_end.dot(ref_dir));
                        let ts = theta_start.rem_euclid(std::f64::consts::TAU);
                        let te = theta_end.rem_euclid(std::f64::consts::TAU);
                        let forward = if te >= ts { te - ts } else { te + std::f64::consts::TAU - ts };
                        if forward > std::f64::consts::PI { same_sense = false; }
                    }
                }
                let (u_start, u_end) = if same_sense {
                    (uv_start, uv_end)
                } else {
                    (uv_end, uv_start)
                };
                let uv_dir = u_end - u_start;
                if uv_dir.length_squared() > 1e-24 {
                    let curve2d_idx = brep.geom.curve2ds.len();
                    brep.geom.curve2ds.push(Curve2d::Line(Line2d {
                        origin: u_start,
                        direction: uv_dir,
                    }));
                    brep.geom.edge_pcurves[ei].insert(0, PCurve {
                        surface_idx: sphere_surface_idx,
                        curve2d_idx,
                    });
                }
            }
            // Add spherical-surface pcurves for the seam edge.
            if let Some(seam_ei) = seam_ei_opt {
                let seam_start_pt = brep.vertices[brep.edges[seam_ei].start].point;
                // Pcurve at u=0 / u=2π: from equator (v=0) to south pole (v=-π/2).
                let u_at_seam = (seam_start_pt - center).y.atan2((seam_start_pt - center).x);
                let v_start = 0.0;             // equator
                let v_end = -std::f64::consts::FRAC_PI_2; // south pole
                for u in [u_at_seam, u_at_seam + two_pi] {
                    let uv_start = DVec2::new(u, v_start);
                    let uv_end = DVec2::new(u, v_end);
                    let uv_dir = uv_end - uv_start;
                    if uv_dir.length_squared() > 1e-24 {
                        let curve2d_idx = brep.geom.curve2ds.len();
                        brep.geom.curve2ds.push(Curve2d::Line(Line2d {
                            origin: uv_start,
                            direction: uv_dir,
                        }));
                        brep.geom.edge_pcurves[seam_ei].push(PCurve {
                            surface_idx: sphere_surface_idx,
                            curve2d_idx,
                        });
                    }
                }
            }
        }
    }

    // Pad auxiliary GeomStore arrays
    while brep.geom.edge_pcurves.len() < brep.edges.len() {
        brep.geom.edge_pcurves.push(vec![]);
    }
    while brep.geom.edge_curve_range.len() < brep.edges.len() {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.face_surface_range.len() < brep.solids[0].shells[0].faces.len() {
        brep.geom.face_surface_range.push(None);
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
    cvi: &mut [Option<usize>; 8],
    box_edge_map: &mut std::collections::HashMap<(usize, usize), usize>,
    corner_indices: &[usize; 4],
    normal: DVec3,
    plane_origin: DVec3,
    outside: bool,
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
                let theta = dy.atan2(dx).rem_euclid(two_pi);
                xs.push(Cross { theta, pos, edge: ei });
            }
        }
    }

    xs.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap_or(std::cmp::Ordering::Equal));
    let mut uxs: Vec<Cross> = Vec::new();
    // Use a liberal tolerance (1e-6) for floating-point duplicates in theta:
    // clamped edge endpoints at the same 3D position have theta differing by
    // ~1.5e-8, and the minimum angular separation of distinct circle-edge
    // intersections on a unit circle is ~0.93 rad → safe at 1e-6.
    let dedup_tol = 1e-6;
    for x in xs {
        let dup = uxs.last().map(|last| {
            let d = (x.theta - last.theta).abs();
            let d_circ = d.min(two_pi - d);
            d_circ <= dedup_tol
        }).unwrap_or(false);
        if !dup { uxs.push(x); }
    }
    // Also dedup the last element against the first (θ≈0 vs θ≈2π).
    if uxs.len() >= 2 {
        let d = (uxs.last().unwrap().theta - uxs[0].theta).abs();
        let d_circ = d.min(two_pi - d);
        if d_circ <= dedup_tol { uxs.pop(); }
    }
    xs = uxs;
    let n_xs = xs.len();
    if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() && n_xs > 2 {
        let label = if normal.x.abs() > 0.5 { if normal.x > 0.0 { "+X" } else { "-X" } }
            else if normal.y.abs() > 0.5 { if normal.y > 0.0 { "+Y" } else { "-Y" } }
            else { if normal.z > 0.0 { "+Z" } else { "-Z" } };
        eprintln!("[N_XS={}] {} n_xs={}", n_xs, label, n_xs);
        for (xi, x) in xs.iter().enumerate() {
            eprintln!("[N_XS={}]   xs[{}] theta={:.10e} pos=({:.2},{:.2},{:.2}) edge={}", n_xs, xi, x.theta, x.pos.x, x.pos.y, x.pos.z, x.edge);
        }
    }

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
    for x in &xs {
        // Reuse box-corner vertices when the intersection point coincides
        // with a corner — avoids duplicate VERTEX_POINT in the STEP output.
        // Creates the corner vertex lazily if it doesn't exist yet.
        let mut reused = None;
        for (ci, corner_vi) in cvi.iter_mut().enumerate() {
            if (x.pos - corners[ci]).length() < 1e-12 {
                if corner_vi.is_none() {
                    *corner_vi = Some(make_vertex(brep, corners[ci]));
                }
                reused = *corner_vi;
                break;
            }
        }
        pt_vis.push(reused.unwrap_or_else(|| make_vertex(brep, x.pos)));
    }

    let point_in_rect = |pt: DVec3| -> bool {
        let u = (pt - corners[c[0]]).dot(x_axis);
        let v = (pt - corners[c[0]]).dot(y_axis);
        let du = (corners[c[1]] - corners[c[0]]).length();
        let dv = (corners[c[3]] - corners[c[0]]).length();
        u >= -1e-12 && u <= du + 1e-12 && v >= -1e-12 && v <= dv + 1e-12
    };

    let mut planar_wes: Vec<WireEdge> = Vec::new();
    let mut arc_edges: Vec<usize> = Vec::new();

    // Find an existing vertex at `p` within tolerance, or create one.
    let find_or_make_vertex = |brep: &mut BRep, p: DVec3| -> usize {
        for (i, v) in brep.vertices.iter().enumerate() {
            if (v.point - p).length_squared() < 1e-20 {
                return i;
            }
        }
        make_vertex(brep, p)
    };
    let mk_line = |brep: &mut BRep, p1: DVec3, p2: DVec3| -> Option<usize> {
        if (p1 - p2).length() < 1e-15 { return None; }
        let v1 = find_or_make_vertex(brep, p1);
        let v2 = find_or_make_vertex(brep, p2);
        let dir = p2 - p1;
        let len = dir.length();
        let curve = Curve3::Line(Line3 { origin: p1, direction: dir / len });
        let ei = make_edge(brep, curve, 0.0, len, v1, v2).ok()?;
        align_edge_geom(brep, ei);
        Some(ei)
    };

    let mut walk_perim = |brep: &mut BRep, from_pos: DVec3, from_edge: usize, to_pos: DVec3, to_edge: usize| -> Option<Vec<WireEdge>> {
        let mut wes = Vec::new();
        let mut cur = from_pos;
        let mut e = from_edge;
        loop {
            // If already at destination, stop.
            if (cur - to_pos).length() < 1e-12 { break; }
            let ci = c[e];
            let cj = c[(e + 1) % 4];
            let nc = corners[cj];
            // If cur is already at the next corner, skip this edge and advance.
            if (cur - nc).length() < 1e-12 {
                if e == to_edge { break; }
                cur = nc;
                e = (e + 1) % 4;
                continue;
            }
            let at_start_corner = (cur - corners[ci]).length() < 1e-12;
            // Handle final segment where destination is mid-edge (not at a corner).
            if at_start_corner && e == to_edge {
                let at_end_corner = (to_pos - corners[cj]).length() < 1e-12;
                if !at_end_corner {
                    if let Some(ei) = mk_line(brep, cur, to_pos) { wes.push(WireEdge::fwd(ei)); }
                    break;
                }
            }
            let ei = if at_start_corner {
                get_or_create_box_edge(brep, box_edge_map, cvi, corners, ci, cj)
            } else if e == to_edge {
                mk_line(brep, cur, to_pos)
            } else {
                mk_line(brep, cur, nc)
            };
            if let Some(ei) = ei {
                let forward = if at_start_corner { ci < cj } else { true };
                wes.push(if forward { WireEdge::fwd(ei) } else { WireEdge::rev(ei) });
            }
            if e == to_edge { break; }
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
            // Arc edge.  For the wrapping interval (j <= i), t_j wraps
            // through 2π and needs +2π so the stored range covers the
            // short arc (e.g. [3π/2, 2π] not [3π/2, 0]).
            let v1 = pt_vis[i];
            let v2 = pt_vis[j];
            let tj_adj = if j > i { t_j } else { t_j + two_pi };
            let ae = make_edge(brep, circle_curve.clone(), t_i, tj_adj, v1, v2).ok()?;
            align_edge_geom(brep, ae);
            arc_edges.push(ae);
            if outside {
                // Union: arc is the hole boundary — walk perimeter FORWARD
                // (increasing edge index) from pi to pj so the planar face
                // excludes the sphere-interior region.
                planar_wes.push(WireEdge::fwd(ae));
                if let Some(w) = walk_perim(brep, pi, xs[i].edge, pj, xs[j].edge) {
                    planar_wes.extend(w);
                }
            } else {
                // Intersection: arc is part of the face boundary.
                planar_wes.push(WireEdge::fwd(ae));
            }
        } else if !outside {
            // Rectangle perimeter (intersection only — union skips arcs
            // whose midpoint is outside the rectangle).
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

    // Add planar pcurves for each arc edge so the STEP writer's SURFACE_CURVE
    // includes a PCURVE on the planar surface.  Without this, shared arcs
    // (reused by the spherical face) would only carry the spherical pcurve,
    // causing STEP readers to show 3/4-circle faces instead of 1/4-circle.
    if !arc_edges.is_empty() {
        let plane_surf_idx = brep.geom.surfaces.len() - 1;
        let center_u = (circle_center - pp).dot(x_axis);
        let center_v = (circle_center - pp).dot(y_axis);
        for &ae in &arc_edges {
            // Retrieve the 3D parameter range from the edge.
            let range = brep.geom.edge_curve_range.get(ae).copied().flatten().unwrap_or([0.0, std::f64::consts::TAU]);
            let curve2d_idx = brep.geom.curve2ds.len();
            brep.geom.curve2ds.push(Curve2d::Circle(Circle2d {
                center: DVec2::new(center_u, center_v),
                radius: circle_r,
            }));
            brep.geom.edge_pcurves[ae].push(PCurve {
                surface_idx: plane_surf_idx,
                curve2d_idx,
            });
        }
    }

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

    // ── 1. Box corner positions (vertices created lazily) ──
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
    let mut cvi: [Option<usize>; 8] = [None; 8]; // created lazily

    // ── 2. Face definitions ──
    let faces: [([usize; 4], DVec3, DVec3); 6] = [
        ([0, 3, 2, 1], DVec3::NEG_Z, DVec3::new(0.0, 0.0, bmin.z)),
        ([4, 5, 6, 7], DVec3::Z,     DVec3::new(0.0, 0.0, bmax.z)),
        ([0, 1, 5, 4], DVec3::NEG_Y, DVec3::new(0.0, bmin.y, 0.0)),
        ([3, 7, 6, 2], DVec3::Y,     DVec3::new(0.0, bmax.y, 0.0)),
        ([0, 4, 7, 3], DVec3::NEG_X, DVec3::new(bmin.x, 0.0, 0.0)),
        ([1, 2, 6, 5], DVec3::X,     DVec3::new(bmax.x, 0.0, 0.0)),
    ];

    // ── 3. Shared box edges (created lazily) ──
    let mut edge_map: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();

    // ── 4. Build planar faces; collect arc edges ──
    let mut all_arcs: Vec<usize> = Vec::new();
    for (fi, &(ref ci, n, pp)) in faces.iter().enumerate() {
        let arcs = build_plane_intersection_face(
            &mut brep, center, radius, &corners, &mut cvi, &mut edge_map, ci, n, pp, false,
        )?;
        all_arcs.extend(arcs);
    }
    if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
        eprintln!("[ANALYTIC] total_arcs={}", all_arcs.len());
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
        // Reuse existing arc edges for the spherical face (reversed direction)
        // so that edges are shared between planar and spherical faces.
        // This allows shell_is_closed to detect a watertight shell and emit
        // CLOSED_SHELL + MANIFOLD_SOLID_BREP instead of OPEN_SHELL.
        let mut sphere_wes: Vec<WireEdge> = Vec::with_capacity(n_arcs);
        for &ei in &ordered {
            sphere_wes.push(WireEdge::rev(ei));
        }
        // Create the spherical face first so the surface is registered in GeomStore.
        make_face(&mut brep, sphere_surf.clone(), make_wire(sphere_wes), vec![]).ok()?;
        // Add spherical-surface pcurves to each arc edge so the STEP writer can
        // emit a complete SURFACE_CURVE with PCURVE on the SPHERICAL_SURFACE.
        let sphere_surface_idx = brep.geom.surfaces.len() - 1;
        if let Surface3::Sphere(sphere) = &brep.geom.surfaces[sphere_surface_idx] {
            for &ei in &ordered {
                let start_pt = brep.vertices[brep.edges[ei].start].point;
                let end_pt = brep.vertices[brep.edges[ei].end].point;
                let uv_start = sphere_uv_propagate(sphere, start_pt, end_pt, radius);
                let uv_end = sphere_uv_propagate(sphere, end_pt, start_pt, radius);
                let uv_dir = uv_end - uv_start;
                if uv_dir.length_squared() > 1e-24 {
                    let curve2d_idx = brep.geom.curve2ds.len();
                    brep.geom.curve2ds.push(Curve2d::Line(Line2d {
                        origin: uv_start,
                        direction: uv_dir,
                    }));
                    // Insert at position 0 so the spherical pcurve is FIRST
                    // in the SURFACE_CURVE's pcurve list.  Viewers that use
                    // .PCURVE_S1. (first pcurve) for surface matching will
                    // then correctly use the spherical pcurve instead of the
                    // planar one.
                    brep.geom.edge_pcurves[ei].insert(0, PCurve {
                        surface_idx: sphere_surface_idx,
                        curve2d_idx,
                    });
                }
            }
        }
    }

    Some(brep)
}
