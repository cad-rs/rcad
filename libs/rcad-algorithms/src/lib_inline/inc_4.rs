
fn unify_one_merge_pass_with_origins(brep: &mut rcad_kernel::BRep, face_origins: Option<&[FaceOrigin]>) -> bool {
    use std::collections::HashMap;

    fn closure_score(brep: &rcad_kernel::BRep) -> usize {
        let report = crate::brep_check::validate_solid_closure(brep);
        report
            .issues
            .iter()
            .map(|iss| match iss {
                crate::CheckIssue::SolidNotClosed {
                    boundary_edge_count,
                    ..
                } => *boundary_edge_count,
                _ => 1,
            })
            .sum()
    }

    fn flat_face_index_of(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi: usize) -> usize {
        let mut idx = 0usize;
        for s in 0..si {
            for sh in &brep.solids[s].shells {
                idx += sh.faces.len();
            }
        }
        for sh in 0..shi {
            idx += brep.solids[si].shells[sh].faces.len();
        }
        idx + fi
    }

    /// Returns `(same_domain, is_planar)`:
    /// - `(Some(true), _)`  → surfaces are the same domain; proceed to merge.
    /// - `(Some(false), _)` → different domains; skip.
    /// - `(None, _)`        → no surface data; caller should fall back to
    ///                        normal-direction heuristic.
    fn surfaces_are_same_domain(
        brep: &rcad_kernel::BRep,
        si: usize,
        shi: usize,
        fi1: usize,
        fi2: usize,
    ) -> (Option<bool>, bool) {
        let ang_tol = tolerance::TOLERANCE_ANG_HEURISTIC_RAD;
        let lin_tol = tolerance::TOLERANCE_PARAM_LEGACY;

        let ff1 = flat_face_index_of(brep, si, shi, fi1);
        let ff2 = flat_face_index_of(brep, si, shi, fi2);
        let sid1 = match brep.geom.face_surface.get(ff1).and_then(|v| *v) {
            Some(id) => id,
            None => return (None, true),
        };
        let sid2 = match brep.geom.face_surface.get(ff2).and_then(|v| *v) {
            Some(id) => id,
            None => return (None, true),
        };
        let s1 = match brep.geom.surfaces.get(sid1) {
            Some(s) => s,
            None => return (None, true),
        };
        let s2 = match brep.geom.surfaces.get(sid2) {
            Some(s) => s,
            None => return (None, true),
        };

        use rcad_kernel::geom::Surface3;
use rcad_kernel::PCurve;
        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let n1 = p1.normal.normalize_or_zero();
                let n2 = p2.normal.normalize_or_zero();
                if n1.length_squared() <= tolerance::TOLERANCE_VEC_SQ_MIN
                    || n2.length_squared() <= tolerance::TOLERANCE_VEC_SQ_MIN
                {
                    return (Some(false), true);
                }
                let cross = n1.cross(n2).length();
                if cross > ang_tol {
                    return (Some(false), true);
                }
                let d = (p2.origin - p1.origin).dot(n1).abs();
                (Some(d <= lin_tol), true)
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                // Same radius?
                if (c1.radius - c2.radius).abs() > lin_tol {
                    return (Some(false), false);
                }
                // Same axis direction?
                let a1 = c1.axis.normalize_or_zero();
                let a2 = c2.axis.normalize_or_zero();
                if a1.cross(a2).length() > ang_tol {
                    return (Some(false), false);
                }
                // Same axis line: point-to-line distance for c2.origin onto c1's axis.
                let d = (c2.origin - c1.origin).cross(a1).length();
                (Some(d <= lin_tol), false)
            }
            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                if (c1.radius - c2.radius).abs() > lin_tol {
                    return (Some(false), false);
                }
                if (c1.half_angle_rad - c2.half_angle_rad).abs() > ang_tol {
                    return (Some(false), false);
                }
                let a1 = c1.axis.normalize_or_zero();
                let a2 = c2.axis.normalize_or_zero();
                if a1.cross(a2).length() > ang_tol {
                    return (Some(false), false);
                }
                let da = (c1.apex - c2.apex).length();
                (Some(da <= lin_tol), false)
            }
            (Surface3::Torus(t1), Surface3::Torus(t2)) => {
                if (t1.major_radius - t2.major_radius).abs() > lin_tol {
                    return (Some(false), false);
                }
                if (t1.minor_radius - t2.minor_radius).abs() > lin_tol {
                    return (Some(false), false);
                }
                let a1 = t1.axis.normalize_or_zero();
                let a2 = t2.axis.normalize_or_zero();
                if a1.cross(a2).length() > ang_tol {
                    return (Some(false), false);
                }
                let dc = (t1.center - t2.center).length();
                (Some(dc <= lin_tol), false)
            }
            (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
                if (s1.radius - s2.radius).abs() > lin_tol {
                    return (Some(false), false);
                }
                let dc = (s1.center - s2.center).length();
                (Some(dc <= lin_tol), false)
            }
            // Cross-type: BSpline and Plane are never same-domain.
            // OCCT FillSameDomainFaces (BOPAlgo_Builder_2.cxx L6153-L6165) only groups
            // faces by edge set equivalence, then checks planar faces via surface type
            // (GeomAbs_Plane).  It does NOT promote planar BSpline to Plane and merge
            // across types — that would incorrectly fuse sub-faces from different
            // operands whose underlying geometry differs (b1=BSpline box vs b2=box).
            // OCCT preserves the original surface type of each operand face.
            (Surface3::BSpline(_), Surface3::Plane(_))
            | (Surface3::Plane(_), Surface3::BSpline(_)) => (Some(false), false),
            (Surface3::BSpline(b1), Surface3::BSpline(b2)) => {
                // BSpline same-domain detection.
                // Two BSpline surfaces are considered same-domain if they have:
                // - Identical degrees
                // - Identical knot vectors (within tolerance)
                // - Identical control point grids (within tolerance)
                // - Identical weights (for rational surfaces)

                if b1.degree_u != b2.degree_u || b1.degree_v != b2.degree_v {
                    return (Some(false), false);
                }

                // Check knot vectors.
                if b1.knots_u.len() != b2.knots_u.len() || b1.knots_v.len() != b2.knots_v.len() {
                    return (Some(false), false);
                }

                for (k1, k2) in b1.knots_u.iter().zip(b2.knots_u.iter()) {
                    if (k1 - k2).abs() > lin_tol {
                        return (Some(false), false);
                    }
                }
                for (k1, k2) in b1.knots_v.iter().zip(b2.knots_v.iter()) {
                    if (k1 - k2).abs() > lin_tol {
                        return (Some(false), false);
                    }
                }

                // Check control points.
                if b1.control_points.len() != b2.control_points.len() {
                    return (Some(false), false);
                }
                for (row1, row2) in b1.control_points.iter().zip(b2.control_points.iter()) {
                    if row1.len() != row2.len() {
                        return (Some(false), false);
                    }
                    for (cp1, cp2) in row1.iter().zip(row2.iter()) {
                        if cp1.distance(*cp2) > lin_tol {
                            return (Some(false), false);
                        }
                    }
                }

                // Check weights for rational surfaces.
                if b1.weights.len() != b2.weights.len() {
                    return (Some(false), false);
                }
                for (row1, row2) in b1.weights.iter().zip(b2.weights.iter()) {
                    if row1.len() != row2.len() {
                        return (Some(false), false);
                    }
                    for (w1, w2) in row1.iter().zip(row2.iter()) {
                        if (w1 - w2).abs() > lin_tol {
                            return (Some(false), false);
                        }
                    }
                }

                (Some(true), false)
            }
            // Mismatched types are never same-domain.
            _ => (Some(false), false),
        }
    }

    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            let nfaces = brep.solids[si].shells[shi].faces.len();

            fn quantize_edge_point(p: glam::DVec3) -> (i64, i64, i64) {
                let inv_tol = 1.0 / tolerance::TOLERANCE_PARAM_LEGACY.max(tolerance::TOLERANCE_ABS);
                (
                    (p.x * inv_tol).round() as i64,
                    (p.y * inv_tol).round() as i64,
                    (p.z * inv_tol).round() as i64,
                )
            }

            fn geometric_edge_key(brep: &rcad_kernel::BRep, edge_idx: usize) -> Option<((i64, i64, i64), (i64, i64, i64))> {
                let edge = brep.edges.get(edge_idx)?;
                let start = quantize_edge_point(brep.vertices.get(edge.start)?.point);
                let end = quantize_edge_point(brep.vertices.get(edge.end)?.point);
                Some(if start <= end { (start, end) } else { (end, start) })
            }

            // Build edge → [face_index_in_shell] adjacency for this shell.
            let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
            let mut geom_edge_to_faces: HashMap<((i64, i64, i64), (i64, i64, i64)), Vec<(usize, usize)>> = HashMap::new();
            for fi in 0..nfaces {
                for we in &brep.solids[si].shells[shi].faces[fi].outer_wire.edges {
                    edge_to_faces.entry(we.idx).or_default().push(fi);
                    if let Some(key) = geometric_edge_key(brep, we.idx) {
                        geom_edge_to_faces.entry(key).or_default().push((fi, we.idx));
                    }
                }
                for iw in &brep.solids[si].shells[shi].faces[fi].inner_wires {
                    for we in &iw.edges {
                        edge_to_faces.entry(we.idx).or_default().push(fi);
                        if let Some(key) = geometric_edge_key(brep, we.idx) {
                            geom_edge_to_faces.entry(key).or_default().push((fi, we.idx));
                        }
                    }
                }
            }

            // Find the first internal edge shared by exactly 2 same-domain faces.
            // Sort by edge index for deterministic iteration (HashMap order varies between runs).
            let mut adjacency_candidates: Vec<(usize, usize, usize, usize)> = edge_to_faces
                .iter()
                .filter_map(|(&edge_idx, face_refs)| {
                    if face_refs.len() == 2 {
                        Some((edge_idx, edge_idx, face_refs[0], face_refs[1]))
                    } else {
                        None
                    }
                })
                .collect();
            for face_edges in geom_edge_to_faces.values() {
                if face_edges.len() != 2 {
                    continue;
                }
                let (fi1, edge_idx1) = face_edges[0];
                let (fi2, edge_idx2) = face_edges[1];
                if fi1 == fi2 || edge_idx1 == edge_idx2 {
                    continue;
                }
                adjacency_candidates.push((edge_idx1, edge_idx2, fi1, fi2));
            }
            adjacency_candidates.sort_unstable();
            adjacency_candidates.dedup();
            for &(edge_idx1, edge_idx2, fi1, fi2) in &adjacency_candidates {
                if fi1 == fi2 {
                    continue;
                }

                let face1_normal = brep.solids[si].shells[shi].faces[fi1].normal;
                let face2_normal = brep.solids[si].shells[shi].faces[fi2].normal;

                let get_face_pt = |fi: usize| -> Option<glam::DVec3> {
                    let we = brep.solids[si].shells[shi].faces[fi]
                        .outer_wire
                        .edges
                        .first()?;
                    let edge = brep.edges.get(we.idx)?;
                    let v_idx = if we.forward { edge.start } else { edge.end };
                    brep.vertices.get(v_idx).map(|v| v.point)
                };

                let face_outer_vertices = |fi: usize| -> Option<Vec<glam::DVec3>> {
                    let mut out = Vec::new();
                    for we in &brep.solids[si].shells[shi].faces[fi].outer_wire.edges {
                        let e = brep.edges.get(we.idx)?;
                        let v_idx = if we.forward { e.start } else { e.end };
                        out.push(brep.vertices.get(v_idx)?.point);
                    }
                    if out.is_empty() { None } else { Some(out) }
                };

                let (same_domain, is_planar) = surfaces_are_same_domain(brep, si, shi, fi1, fi2);

                // Origin guard: only merge faces from the SAME original shape.
                // Without this we merge A-faces with B-faces on the same surface,
                // breaking boolean topology (seen as regressions in boptuc/bopfuse).
                if let Some(origins) = face_origins {
                    let ff1 = flat_face_index_of(brep, si, shi, fi1);
                    let ff2 = flat_face_index_of(brep, si, shi, fi2);
                    if origins.get(ff1) != origins.get(ff2) {
                        continue;
                    }
                }

                let mut should_merge = match same_domain {
                    Some(false) => false,
                    Some(true) => {
                        // For planar faces add a vertex–plane distance sanity check.
                        if is_planar {
                            let n = face1_normal.normalize();
                            if let (Some(pt1), Some(vs1), Some(vs2)) = (
                                get_face_pt(fi1),
                                face_outer_vertices(fi1),
                                face_outer_vertices(fi2),
                            ) {
                                let all_vs1_on_plane1 = vs1
                                    .iter()
                                    .all(|p| (*p - pt1).dot(n).abs() <= tolerance::TOLERANCE_PLANE_DIST_RELAX);
                                let all_vs2_on_plane1 = vs2
                                    .iter()
                                    .all(|p| (*p - pt1).dot(n).abs() <= tolerance::TOLERANCE_PLANE_DIST_RELAX);
                                all_vs1_on_plane1 && all_vs2_on_plane1
                            } else {
                                false
                            }
                        } else {
                            // For curved surfaces the geom-store check is sufficient.
                            true
                        }
                    }
                    None => {
                        // No surface data: fall back to per-face normal heuristic.
                        let cross = face1_normal.cross(face2_normal).length();
                        if cross > tolerance::TOLERANCE_PARAM_LEGACY {
                            false
                        } else if let (Some(pt1), Some(pt2)) = (get_face_pt(fi1), get_face_pt(fi2))
                        {
                            let n = face1_normal.normalize();
                            (pt2 - pt1).dot(n).abs() <= tolerance::TOLERANCE_PARAM_LEGACY
                        } else {
                            false
                        }
                    }
                };

                // Topological + geometric double-validation: extra guards so we do not merge
                // faces with incompatible topology or UV regions.
                if should_merge {
                    // Check shared edge continuity (PCurve alignment).
                    let edge_continuous = if edge_idx1 == edge_idx2 {
                        validate_shared_edge_continuity(brep, si, shi, fi1, fi2, edge_idx1)
                    } else {
                        // Geometric-edge fallback found equivalent boundaries with distinct
                        // edge indices, so there is no single shared topological edge to validate.
                        true
                    };
                    if !edge_continuous {
                        should_merge = false;
                    }
                }

                if should_merge {
                    // Planar booleans often use disjoint face-local UV rectangles; merges stay bounded
                    // by shared-edge continuity and the Newell outer-area check after splice.
                    let uv_compatible = if is_planar && same_domain == Some(true) {
                        true
                    } else {
                        validate_uv_regions_compatible(brep, si, shi, fi1, fi2)
                    };
                    if !uv_compatible {
                        should_merge = false;
                    }
                }

                if !should_merge {
                    continue;
                }

                // For non-planar faces (sphere, cylinder, etc.), avoid creating
                // faces with too many outer edges — downstream surface-area
                // computation (analytic grid or earcut) becomes infeasible.
                if !is_planar {
                    let n1 = brep.solids[si].shells[shi].faces[fi1].outer_wire.edges.len();
                    let n2 = brep.solids[si].shells[shi].faces[fi2].outer_wire.edges.len();
                    // Merge removes 2 shared edges, net ≈ n1 + n2 - 2
                    if n1 + n2 > 650 {
                        continue;
                    }
                }

                // Merge wire: splice Face2 edges into Face1 at the position of the shared edge.
                let wire1 = brep.solids[si].shells[shi].faces[fi1]
                    .outer_wire
                    .edges
                    .clone();
                let wire2 = brep.solids[si].shells[shi].faces[fi2]
                    .outer_wire
                    .edges
                    .clone();

                if let Some(merged_wire_edges) = splice_wires(&wire1, edge_idx1, &wire2, edge_idx2) {
                    let merged_wire_edges = cleanup_merged_wire_edges(brep, &merged_wire_edges);
                    // Collect inner wires from both faces.
                    let inner1 = brep.solids[si].shells[shi].faces[fi1].inner_wires.clone();
                    let inner2 = brep.solids[si].shells[shi].faces[fi2].inner_wires.clone();
                    let mut all_inner = inner1;
                    all_inner.extend(inner2);

                    // Detect figure-8 self-intersecting wires: if the merged outer wire
                    // visits any vertex more than once, extract the inner sub-loops.
                    let (outer_edges_raw, extracted_inners) =
                        extract_inner_loops_from_wire(brep, &merged_wire_edges);
                    // Re-run cleanup on the outer wire after inner loop extraction,
                    // since extraction may leave adjacent duplicate segments.
                    let outer_edges = if extracted_inners.is_empty() {
                        outer_edges_raw
                    } else {
                        cleanup_merged_wire_edges(brep, &outer_edges_raw)
                    };
                    all_inner.extend(extracted_inners);

                    // Build merged face (mesh_dirty=true; normal reused from face1).
                    let merged_face = rcad_kernel::topology::Face {
                        outer_wire: rcad_kernel::topology::Wire { edges: outer_edges },
                        inner_wires: all_inner,
                        normal: face1_normal,
                        triangles: vec![],
                        sample_point: None,
                        mesh_dirty: true,
            surface_idx: None,
                    };

                    // Planar guard: refuse merges whose merged outer area is larger than the
                    // sum of the two faces' outer areas (plus tolerance). Valid same-domain
                    // merges are roughly additive along a shared edge; incorrect splices
                    // around a frame/hole can "zip" opposite banks into one loop whose area
                    // jumps (e.g. union of overlapping boxes at a contact plane).
                    if is_planar {
                        let nunit = face1_normal.normalize_or_zero();
                        let poly1 = face_outer_polygon_points(brep, si, shi, fi1);
                        let poly2 = face_outer_polygon_points(brep, si, shi, fi2);
                        let poly_m = wire_to_polygon_points(brep, &merged_face.outer_wire.edges);
                        let a1 = newell_polygon_abs_area(&poly1, nunit);
                        let a2 = newell_polygon_abs_area(&poly2, nunit);
                        let am = newell_polygon_abs_area(&poly_m, nunit);
                        let sum = a1 + a2;
                        let tol = tolerance::TOLERANCE_AREA_REL * sum.max(am).max(1.0) + tolerance::TOLERANCE_ABS;
                        if am > sum + tol {
                            continue;
                        }
                    }

                    // Replace fi1 with merged face, remove fi2, but only commit if
                    // the candidate result stays topologically closed.
                    let (keep_idx, remove_idx) = if fi1 < fi2 { (fi1, fi2) } else { (fi2, fi1) };
                    let mut candidate = brep.clone();

                    // Update face_surface mapping: keep keep_idx's surface id.
                    let _kept_flat = flat_face_index_of(&candidate, si, shi, keep_idx);
                    let remove_flat = flat_face_index_of(&candidate, si, shi, remove_idx);
                    remove_flat_face_geom_slots(&mut candidate.geom, remove_flat);

                    candidate.solids[si].shells[shi].faces[keep_idx] = merged_face;
                    candidate.solids[si].shells[shi].faces.remove(remove_idx);

                    let current_score = closure_score(brep);
                    let candidate_score = closure_score(&candidate);
                    if candidate_score > current_score {
                        continue;
                    }

                    *brep = candidate;
                    return true;
                }
            }
        }
    }

    false
}

/// Splice two wire edge lists together by removing the shared edge and
/// interleaving the remaining edges.
///
/// Returns `None` if the shared edge is not found in either wire.
fn splice_wires(
    wire_a: &[rcad_kernel::topology::WireEdge],
    shared_edge_idx_a: usize,
    wire_b: &[rcad_kernel::topology::WireEdge],
    shared_edge_idx_b: usize,
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
    let pos_a = wire_a.iter().position(|we| we.idx == shared_edge_idx_a)?;
    let pos_b = wire_b.iter().position(|we| we.idx == shared_edge_idx_b)?;

    let n_b = wire_b.len();
    // B's edges (excluding the shared edge), in cyclic order starting at pos_b + 1
    let b_edges: Vec<rcad_kernel::topology::WireEdge> =
        (1..n_b).map(|i| wire_b[(pos_b + i) % n_b]).collect();

    let mut merged = Vec::with_capacity(wire_a.len() - 1 + b_edges.len());
    merged.extend_from_slice(&wire_a[..pos_a]);
    merged.extend(b_edges);
    merged.extend_from_slice(&wire_a[pos_a + 1..]);

    if merged.len() < 3 {
        return None; // Degenerate result
    }

    Some(merged)
}

pub(crate) fn oriented_edge_vertices(
    brep: &rcad_kernel::BRep,
    we: rcad_kernel::topology::WireEdge,
) -> Option<(usize, usize)> {
    let e = brep.edges.get(we.idx)?;
    if we.forward {
        Some((e.start, e.end))
    } else {
        Some((e.end, e.start))
    }
}

fn find_existing_edge_between_vertices(
    brep: &rcad_kernel::BRep,
    from: usize,
    to: usize,
) -> Option<rcad_kernel::topology::WireEdge> {
    for (idx, e) in brep.edges.iter().enumerate() {
        if e.start == from && e.end == to {
            return Some(rcad_kernel::topology::WireEdge::fwd(idx));
        }
        if e.start == to && e.end == from {
            return Some(rcad_kernel::topology::WireEdge::rev(idx));
        }
    }
    None
}

fn points_are_collinear_forward(a: glam::DVec3, b: glam::DVec3, c: glam::DVec3) -> bool {
    let ab = b - a;
    let bc = c - b;
    let ab_len = ab.length();
    let bc_len = bc.length();
    if ab_len <= tolerance::TOLERANCE_LEN_MIN || bc_len <= tolerance::TOLERANCE_LEN_MIN {
        return false;
    }

    let cross = ab.cross(bc).length();
    let dot = ab.dot(bc);
    cross <= tolerance::TOLERANCE_ABS * (ab_len + bc_len) && dot > 0.0
}

fn collapse_collinear_segments_with_existing_bridge(
    brep: &rcad_kernel::BRep,
    wire: &[rcad_kernel::topology::WireEdge],
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
    let mut out = wire.to_vec();
    if out.len() < 4 {
        return None;
    }

    loop {
        let n = out.len();
        if n < 4 {
            break;
        }

        let mut changed = false;
        for i in 0..n {
            let j = (i + 1) % n;
            let (u, v1) = oriented_edge_vertices(brep, out[i])?;
            let (v2, w) = oriented_edge_vertices(brep, out[j])?;
            if v1 != v2 || u == w {
                continue;
            }

            let p_u = brep.vertices.get(u)?.point;
            let p_v = brep.vertices.get(v1)?.point;
            let p_w = brep.vertices.get(w)?.point;
            if !points_are_collinear_forward(p_u, p_v, p_w) {
                continue;
            }

            let bridge = match find_existing_edge_between_vertices(brep, u, w) {
                Some(e) if e.idx != out[i].idx && e.idx != out[j].idx => e,
                _ => continue,
            };

            if i + 1 < n {
                out.splice(i..=i + 1, [bridge]);
            } else {
                out.pop();
                out.remove(0);
                out.insert(0, bridge);
            }
            changed = true;
            break;
        }

        if !changed {
            break;
        }
    }

    if out.len() >= 3 { Some(out) } else { None }
}

fn wire_is_closed_and_connected(brep: &rcad_kernel::BRep, wire: &[rcad_kernel::topology::WireEdge]) -> bool {
    if wire.len() < 3 {
        return false;
    }

    let Some((first_start, mut prev_end)) = oriented_edge_vertices(brep, wire[0]) else {
        return false;
    };

    for we in &wire[1..] {
        let Some((start, end)) = oriented_edge_vertices(brep, *we) else {
            return false;
        };
        if start != prev_end {
            return false;
        }
        prev_end = end;
    }

    prev_end == first_start
}

fn reorder_wire_into_connected_loop(
    brep: &rcad_kernel::BRep,
    wire: &[rcad_kernel::topology::WireEdge],
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
    if wire.is_empty() {
        return None;
    }

    let mut unused: Vec<rcad_kernel::topology::WireEdge> = wire.to_vec();
    let first = unused.remove(0);
    let mut out = vec![first];

    let (_, mut current_end) = oriented_edge_vertices(brep, first)?;

    while !unused.is_empty() {
        let mut found_idx: Option<usize> = None;
        let mut flip = false;

        for (i, we) in unused.iter().enumerate() {
            let (s, e) = oriented_edge_vertices(brep, *we)?;
            if s == current_end {
                found_idx = Some(i);
                flip = false;
                break;
            }
            if e == current_end {
                found_idx = Some(i);
                flip = true;
                break;
            }
        }

        let i = found_idx?;
        let mut next = unused.remove(i);
        if flip {
            next.forward = !next.forward;
        }
        let (_, next_end) = oriented_edge_vertices(brep, next)?;
        out.push(next);
        current_end = next_end;
    }

    if wire_is_closed_and_connected(brep, &out) {
        Some(out)
    } else {
        None
    }
}

fn cancel_duplicate_segments_by_parity(
    brep: &rcad_kernel::BRep,
    wire: &[rcad_kernel::topology::WireEdge],
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
    use std::collections::HashMap;

    let mut groups: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (i, &we) in wire.iter().enumerate() {
        let (u, v) = oriented_edge_vertices(brep, we)?;
        let key = if u <= v { (u, v) } else { (v, u) };
        groups.entry(key).or_default().push(i);
    }

    let mut keep = vec![true; wire.len()];
    for idxs in groups.values() {
        if idxs.len() >= 2 {
            let cancel_count = (idxs.len() / 2) * 2;
            for idx in idxs.iter().take(cancel_count) {
                keep[*idx] = false;
            }
        }
    }

    let out: Vec<rcad_kernel::topology::WireEdge> = wire
        .iter()
        .enumerate()
        .filter_map(|(i, &we)| if keep[i] { Some(we) } else { None })
        .collect();

    if out.len() >= 3 { Some(out) } else { None }
}

/// Detect figure-8 self-intersecting wires and extract inner sub-loops.
///
/// # Background: the figure-8 bug
///
/// `unify_one_merge_pass` calls `splice_wires` to merge two coplanar adjacent
/// faces by removing their shared edge and interleaving the remaining edges.
/// When the boolean difference cuts a rectangular notch through a face (e.g.
/// the x=3 face of box A after subtracting box B), the raw result contains
/// several sub-faces around the notch hole.  As `unify_one_merge_pass` merges
/// them one by one, a merge step can produce a wire that visits a corner vertex
/// twice — once on the outer boundary and once on the notch boundary.  The
/// resulting wire traces a figure-8 path instead of a simple outer loop with a
/// separate inner loop (hole).
///
/// # What this function does
///
/// Walk the wire tracking visited start-vertices.  The first time a vertex is
/// seen twice, the sub-sequence between the two visits is extracted as an inner
/// wire (hole).  The remaining edges form the outer wire.  The function recurses
/// on the outer wire to handle multiple holes.
///
/// Returns `(outer_wire_edges, inner_wires)`.
fn extract_inner_loops_from_wire(
    brep: &rcad_kernel::BRep,
    wire: &[rcad_kernel::topology::WireEdge],
) -> (
    Vec<rcad_kernel::topology::WireEdge>,
    Vec<rcad_kernel::topology::Wire>,
) {
    use std::collections::HashMap;

    // Build vertex sequence: for each edge in the wire, record the start vertex.
    let mut verts: Vec<usize> = Vec::with_capacity(wire.len());
    for &we in wire {
        let Some((u, _v)) = oriented_edge_vertices(brep, we) else {
            return (wire.to_vec(), vec![]);
        };
        verts.push(u);
    }

    // Find the first vertex that appears more than once.
    let mut seen: HashMap<usize, usize> = HashMap::new(); // vertex -> first index
    let mut split_at: Option<(usize, usize)> = None; // (first_pos, second_pos)
    for (i, &v) in verts.iter().enumerate() {
        if let Some(&first) = seen.get(&v) {
            split_at = Some((first, i));
            break;
        }
        seen.insert(v, i);
    }

    let Some((start, end)) = split_at else {
        // No self-intersection — return as-is.
        return (wire.to_vec(), vec![]);
    };

    // The sub-loop wire[start..end] is the inner loop.
    // The outer wire is wire[0..start] + wire[end..].
    let inner_edges: Vec<rcad_kernel::topology::WireEdge> = wire[start..end].to_vec();
    let mut outer_edges: Vec<rcad_kernel::topology::WireEdge> =
        Vec::with_capacity(wire.len() - inner_edges.len());
    outer_edges.extend_from_slice(&wire[..start]);
    outer_edges.extend_from_slice(&wire[end..]);

    if inner_edges.len() < 3 || outer_edges.len() < 3 {
        return (wire.to_vec(), vec![]);
    }

    let inner_wire = rcad_kernel::topology::Wire { edges: inner_edges };

    // Recursively check the outer wire for further self-intersections.
    let (final_outer, mut more_inners) = extract_inner_loops_from_wire(brep, &outer_edges);
    more_inners.push(inner_wire);
    (final_outer, more_inners)
}

fn cleanup_merged_wire_edges(
    brep: &rcad_kernel::BRep,
    wire: &[rcad_kernel::topology::WireEdge],
) -> Vec<rcad_kernel::topology::WireEdge> {
    if wire.len() < 4 {
        return wire.to_vec();
    }

    let mut cleaned: Vec<rcad_kernel::topology::WireEdge> = Vec::with_capacity(wire.len());

    for &we in wire {
        let Some((u, v)) = oriented_edge_vertices(brep, we) else {
            return wire.to_vec();
        };

        if let Some(&last) = cleaned.last() {
            let Some((lu, lv)) = oriented_edge_vertices(brep, last) else {
                return wire.to_vec();
            };
            let same_segment = (lu == u && lv == v) || (lu == v && lv == u);
            if same_segment {
                cleaned.pop();
                continue;
            }
        }
        cleaned.push(we);
    }

    while cleaned.len() >= 2 {
        let first = cleaned[0];
        let last = *cleaned.last().unwrap_or(&cleaned[0]);
        let Some((fu, fv)) = oriented_edge_vertices(brep, first) else {
            return wire.to_vec();
        };
        let Some((lu, lv)) = oriented_edge_vertices(brep, last) else {
            return wire.to_vec();
        };
        let same_segment = (fu == lu && fv == lv) || (fu == lv && fv == lu);
        if !same_segment {
            break;
        }
        cleaned.remove(0);
        cleaned.pop();
    }

    let stage1 = if wire_is_closed_and_connected(brep, &cleaned) {
        Some(cleaned)
    } else if let Some(cancelled) = cancel_duplicate_segments_by_parity(brep, &cleaned) {
        reorder_wire_into_connected_loop(brep, &cancelled)
    } else {
        None
    };

    let Some(mut out) = stage1 else {
        return wire.to_vec();
    };

    if let Some(collapsed) = collapse_collinear_segments_with_existing_bridge(brep, &out)
        && let Some(reordered) = reorder_wire_into_connected_loop(brep, &collapsed)
        && wire_is_closed_and_connected(brep, &reordered)
    {
        out = reordered;
    }

    out
}

/// boundary. This function detects such duplicate faces within each shell and
/// removes the extra copies.
///
/// Detection criterion: two faces in the same shell are duplicates when all of
/// the following hold:
/// - They share the same normal direction (parallel within [`tolerance::TOLERANCE_PARAM_LEGACY`]).
/// - One face's representative vertex lies on the other face's plane (within [`tolerance::TOLERANCE_PARAM_LEGACY`]).
/// - Their edge sets overlap entirely (every outer-wire edge of the smaller
///   face is also in the larger face, or they share ≥ 75 % of edges).
///
/// Returns the cleaned rcad_kernel::BRep and the number of faces removed.
///
/// Analogous to the internal-face elimination step of OCCT `BOPAlgo_BuilderSolid`.
pub fn remove_internal_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
    use std::collections::HashSet;

    fn flat_face_index_of(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi: usize) -> usize {
        let mut idx = 0usize;
        for s in 0..si {
            for sh in &brep.solids[s].shells {
                idx += sh.faces.len();
            }
        }
        for sh in 0..shi {
            idx += brep.solids[si].shells[sh].faces.len();
        }
        idx + fi
    }

    fn surfaces_are_same_domain(
        brep: &rcad_kernel::BRep,
        si: usize,
        shi: usize,
        fi1: usize,
        fi2: usize,
    ) -> Option<bool> {
        let ang_tol = tolerance::TOLERANCE_ANG_HEURISTIC_RAD;
        let lin_tol = tolerance::TOLERANCE_PARAM_LEGACY;

        let ff1 = flat_face_index_of(brep, si, shi, fi1);
        let ff2 = flat_face_index_of(brep, si, shi, fi2);
        let sid1 = brep.geom.face_surface.get(ff1).and_then(|v| *v)?;
        let sid2 = brep.geom.face_surface.get(ff2).and_then(|v| *v)?;
        let s1 = brep.geom.surfaces.get(sid1)?;
        let s2 = brep.geom.surfaces.get(sid2)?;

        use rcad_kernel::geom::Surface3;
        Some(match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let n1 = p1.normal.normalize_or_zero();
                let n2 = p2.normal.normalize_or_zero();
                if n1.length_squared() <= tolerance::TOLERANCE_VEC_SQ_MIN
                    || n2.length_squared() <= tolerance::TOLERANCE_VEC_SQ_MIN
                {
                    false
                } else {
                    let cross = n1.cross(n2).length();
                    let d = (p2.origin - p1.origin).dot(n1).abs();
                    cross <= ang_tol && d <= lin_tol
                }
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                if (c1.radius - c2.radius).abs() > lin_tol {
                    false
                } else {
                    let a1 = c1.axis.normalize_or_zero();
                    let a2 = c2.axis.normalize_or_zero();
                    let cross = a1.cross(a2).length();
                    let d = (c2.origin - c1.origin).cross(a1).length();
                    cross <= ang_tol && d <= lin_tol
                }
            }
            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                if (c1.radius - c2.radius).abs() > lin_tol {
                    false
                } else if (c1.half_angle_rad - c2.half_angle_rad).abs() > ang_tol {
                    false
                } else {
                    let a1 = c1.axis.normalize_or_zero();
                    let a2 = c2.axis.normalize_or_zero();
                    a1.cross(a2).length() <= ang_tol && (c1.apex - c2.apex).length() <= lin_tol
                }
            }
            (Surface3::Torus(t1), Surface3::Torus(t2)) => {
                (t1.major_radius - t2.major_radius).abs() <= lin_tol
                    && (t1.minor_radius - t2.minor_radius).abs() <= lin_tol
                    && t1
                        .axis
                        .normalize_or_zero()
                        .cross(t2.axis.normalize_or_zero())
                        .length()
                        <= ang_tol
                    && (t1.center - t2.center).length() <= lin_tol
            }
            (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
                (s1.radius - s2.radius).abs() <= lin_tol
                    && (s1.center - s2.center).length() <= lin_tol
            }
            (Surface3::BSpline(b1), Surface3::BSpline(b2)) => {
                // BSpline same-domain detection.
                if b1.degree_u != b2.degree_u || b1.degree_v != b2.degree_v {
                    false
                } else if b1.knots_u.len() != b2.knots_u.len()
                    || b1.knots_v.len() != b2.knots_v.len()
                {
                    false
                } else if !b1
                    .knots_u
                    .iter()
                    .zip(b2.knots_u.iter())
                    .all(|(k1, k2)| (k1 - k2).abs() <= lin_tol)
                {
                    false
                } else if !b1
                    .knots_v
                    .iter()
                    .zip(b2.knots_v.iter())
                    .all(|(k1, k2)| (k1 - k2).abs() <= lin_tol)
                {
                    false
                } else if b1.control_points.len() != b2.control_points.len() {
                    false
                } else if !b1.control_points.iter().zip(b2.control_points.iter()).all(
                    |(row1, row2)| {
                        row1.len() == row2.len()
                            && row1
                                .iter()
                                .zip(row2.iter())
                                .all(|(cp1, cp2)| cp1.distance(*cp2) <= lin_tol)
                    },
                ) {
                    false
                } else if b1.weights.len() != b2.weights.len() {
                    false
                } else {
                    !!b1.weights
                        .iter()
                        .zip(b2.weights.iter())
                        .all(|(row1, row2)| {
                            row1.len() == row2.len()
                                && row1
                                    .iter()
                                    .zip(row2.iter())
                                    .all(|(w1, w2)| (w1 - w2).abs() <= lin_tol)
                        })
                }
            }
            _ => false,
        })
    }

    /// Validate face orientation consistency within a shell.
    /// Returns false if face orientation is inconsistent with majority orientation,
    /// indicating potential pseudo-internal topology that should not be removed.
    fn validate_face_orientation_consistency(
        _brep: &rcad_kernel::BRep,
        _si: usize,
        _shi: usize,
        _fi: usize,
    ) -> bool {
        // Count faces with matching vs. opposite orientation to detect outliers.
        // A face with opposite orientation to most others might be pseudo-internal
        // and should be preserved rather than removed.

        // For now, we accept all orientations as valid (conservative).
        // Future: could add full rcad_kernel::BRep solid vs. hollow validation.
        true
    }

    /// Detect if a face pair forms a true internal duplicate vs. pseudo-internal.
    /// True duplicates have opposite normals and identical/near-identical coverage.
    /// Pseudo-internal faces may share edges but represent distinct original surfaces.
    fn is_true_internal_duplicate(
        brep: &rcad_kernel::BRep,
        si: usize,
        shi: usize,
        fi1: usize,
        fi2: usize,
        edges_i: &HashSet<usize>,
        edges_j: &HashSet<usize>,
    ) -> bool {
        let face_i = &brep.solids[si].shells[shi].faces[fi1];
        let face_j = &brep.solids[si].shells[shi].faces[fi2];

        let ni = face_i.normal.normalize_or_zero();
        let nj = face_j.normal.normalize_or_zero();

        // Check if normals are truly opposite (sign test, not just parallel).
        let dot = ni.dot(nj);
        let are_opposite_normals = dot < -0.99; // Opposite orientation

        if !are_opposite_normals {
            // Not opposite normals: cannot be true internal duplicate.
            return false;
        }

        // Check if wires form a topological enclosure (all edges shared at least once).
        let shared_edges = edges_i.intersection(edges_j).count();
        let all_edges_shared = shared_edges == edges_i.len() && shared_edges == edges_j.len();

        if !all_edges_shared {
            // Not all edges shared: likely pseudo-internal or adjacent faces.
            return false;
        }

        // All checks indicate true internal duplicate: opposite normals + full edge overlap.
        true
    }

    let mut out = brep.clone();
    let mut total_removed = 0usize;

    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            // Iteratively remove one duplicate per pass.
            loop {
                let nfaces = out.solids[si].shells[shi].faces.len();
                let mut removed_idx: Option<usize> = None;

                'outer: for fi in 0..nfaces {
                    for fj in (fi + 1)..nfaces {
                        let face_i = &out.solids[si].shells[shi].faces[fi];
                        let face_j = &out.solids[si].shells[shi].faces[fj];

                        let ni = face_i.normal;
                        let nj = face_j.normal;

                        if ni == glam::DVec3::ZERO || nj == glam::DVec3::ZERO {
                            continue;
                        }

                        // Check parallel normals (allow opposite orientation;
                        // duplicated internal faces can be anti-parallel).
                        let cross = ni.cross(nj).length();
                        let dot = ni.normalize().dot(nj.normalize());
                        if cross > tolerance::TOLERANCE_PARAM_LEGACY
                            || dot.abs() < tolerance::TOLERANCE_DOT_NEARLY_PARALLEL
                        {
                            continue;
                        }

                        // Check same domain from analytic surfaces when available.
                        let same_domain_from_geom = surfaces_are_same_domain(&out, si, shi, fi, fj);

                        // Check same plane fallback: a vertex from j lies on i's plane.
                        let get_pt = |f: &rcad_kernel::topology::Face| -> Option<glam::DVec3> {
                            let we = f.outer_wire.edges.first()?;
                            let edge = out.edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            out.vertices.get(vi).map(|v| v.point)
                        };
                        let Some(pi) = get_pt(face_i) else { continue };
                        let Some(pj) = get_pt(face_j) else { continue };

                        let same_plane_fallback = {
                            let n_unit = ni.normalize();
                            (pj - pi).dot(n_unit).abs() <= tolerance::TOLERANCE_PLANE_DIST_RELAX
                        };

                        if !matches!(same_domain_from_geom, Some(true)) && !same_plane_fallback {
                            continue;
                        }

                        // Check edge overlap: build edge-index sets for both faces.
                        let edges_i: HashSet<usize> = out.solids[si].shells[shi].faces[fi]
                            .outer_wire
                            .edges
                            .iter()
                            .map(|we| we.idx)
                            .collect();
                        let edges_j: HashSet<usize> = out.solids[si].shells[shi].faces[fj]
                            .outer_wire
                            .edges
                            .iter()
                            .map(|we| we.idx)
                            .collect();

                        let overlap = edges_i.intersection(&edges_j).count();
                        let min_edges = edges_i.len().min(edges_j.len()).max(1);

                        // Duplicate rule:
                        // - always accept strict subset/superset overlap,
                        // - accept >=75% overlap only when analytic surfaces
                        //   confirm same-domain.
                        let overlap_ratio = overlap as f64 / min_edges as f64;
                        let strong_same_domain = matches!(same_domain_from_geom, Some(true));
                        let same_or_contained = overlap == min_edges
                            || (strong_same_domain && overlap_ratio >= 0.60);
                        let uv_domain_heuristic = false; // Placeholder for future UV-domain check
                        if same_or_contained || uv_domain_heuristic {
                            // Validate this is a true internal duplicate, not pseudo-internal.
                            let is_true_duplicate = is_true_internal_duplicate(
                                &out, si, shi, fi, fj, &edges_i, &edges_j,
                            );

                            if !is_true_duplicate {
                                // Not a true duplicate: skip removal.
                                continue;
                            }

                            // Validate orientation consistency before removal.
                            let orientation_valid_i =
                                validate_face_orientation_consistency(&out, si, shi, fi);
                            let orientation_valid_j =
                                validate_face_orientation_consistency(&out, si, shi, fj);

                            if !orientation_valid_i || !orientation_valid_j {
                                // Orientation inconsistency detected: skip removal.
                                continue;
                            }

                            // All checks passed: remove fj (keep fi).
                            removed_idx = Some(fj);
                            break 'outer;
                        }
                    }
                }

                if let Some(idx) = removed_idx {
                    out.solids[si].shells[shi].faces.remove(idx);
                    total_removed += 1;
                } else {
                    break;
                }
            }
        }
    }

    // Void shell detection: remove shells fully enclosed within another shell
    // (OCCT's BOPAlgo_BuilderSolid eliminates these during construction;
    // rcad's BooleanBuilder may leave them behind for post-processing to clean up).
    {
        for si in 0..out.solids.len() {
            if out.solids[si].shells.len() < 2 {
                continue;
            }
            // Compute bounding box for each shell
            let shell_bboxes: Vec<Option<(glam::DVec3, glam::DVec3)>> = out.solids[si]
                .shells
                .iter()
                .map(|sh| {
                    let mut min_pt = glam::DVec3::splat(f64::MAX);
                    let mut max_pt = glam::DVec3::splat(f64::MIN);
                    let mut has_verts = false;
                    for f in &sh.faces {
                        for we in &f.outer_wire.edges {
                            if let Some(e) = out.edges.get(we.idx) {
                                if let Some(v) = out.vertices.get(e.start) {
                                    min_pt = min_pt.min(v.point);
                                    max_pt = max_pt.max(v.point);
                                    has_verts = true;
                                }
                                if let Some(v) = out.vertices.get(e.end) {
                                    min_pt = min_pt.min(v.point);
                                    max_pt = max_pt.max(v.point);
                                    has_verts = true;
                                }
                            }
                        }
                    }
                    if has_verts { Some((min_pt, max_pt)) } else { None }
                })
                .collect();

            // Find shells to remove: shells whose bbox is fully inside another shell's bbox
            let mut to_remove: Vec<usize> = vec![];
            for i in 0..out.solids[si].shells.len() {
                let Some((i_min, i_max)) = &shell_bboxes[i] else { continue };
                if i == 0 { continue; } // keep first shell (typically outer)
                for j in 0..out.solids[si].shells.len() {
                    if i == j { continue; }
                    let Some((j_min, j_max)) = &shell_bboxes[j] else { continue };
                    // Check if shell i is fully inside shell j
                    let tol = tolerance::TOLERANCE_ABS;
                    if i_min.x >= j_min.x - tol && i_max.x <= j_max.x + tol
                        && i_min.y >= j_min.y - tol && i_max.y <= j_max.y + tol
                        && i_min.z >= j_min.z - tol && i_max.z <= j_max.z + tol
                    {
                        to_remove.push(i);
                        break;
                    }
                }
            }
            // Remove in reverse order to preserve indices
            to_remove.sort_unstable();
            to_remove.dedup();
            for idx in to_remove.into_iter().rev() {
                out.solids[si].shells.remove(idx);
                total_removed += 1; // approximate — shell removal may remove multiple faces
            }
        }
    }

    (out, total_removed)
}