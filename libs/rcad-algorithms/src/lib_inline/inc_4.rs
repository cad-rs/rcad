// inc_4.rs — Face merging and internal-face removal.
// Included inline from lib.rs. All crate-level imports are available.

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::topods::{BRep, ShapeRef, TShape, Orientation};
use rcad_kernel::topology::{WireEdge, Wire};

// ---------------------------------------------------------------------------
// Local TShape-access helpers for this file
// ---------------------------------------------------------------------------

fn vpoint(brep: &BRep, vi: usize) -> DVec3 {
    brep.vertex_point(vi).unwrap_or(DVec3::ZERO)
}

fn edge_start(brep: &BRep, ei: usize) -> usize {
    match &*brep.tshapes[ei] { TShape::Edge(ed) => ed.first.index, _ => 0 }
}

fn edge_end(brep: &BRep, ei: usize) -> usize {
    match &*brep.tshapes[ei] { TShape::Edge(ed) => ed.last.index, _ => 0 }
}

fn face_edge_refs(brep: &BRep, face_sr: ShapeRef) -> Vec<usize> {
    let mut out = Vec::new();
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return out };
    let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else { return out };
    for esr in &wd.edges { out.push(esr.index); }
    for iw_sr in &fd.inner_wires {
        let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { continue };
        for esr in &iwd.edges { out.push(esr.index); }
    }
    out
}

fn face_outer_edge_refs(brep: &BRep, face_sr: ShapeRef) -> Vec<(usize, bool)> {
    let mut out = Vec::new();
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return out };
    let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else { return out };
    for esr in &wd.edges {
        out.push((esr.index, esr.orientation == Orientation::Forward));
    }
    out
}

fn face_inner_edge_refs(brep: &BRep, face_sr: ShapeRef) -> Vec<Vec<(usize, bool)>> {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return vec![] };
    fd.inner_wires.iter().filter_map(|iw_sr| {
        let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { return None };
        Some(iwd.edges.iter().map(|esr| (esr.index, esr.orientation == Orientation::Forward)).collect())
    }).collect()
}

fn face_outer_wire_edges(brep: &BRep, face_sr: ShapeRef) -> Vec<WireEdge> {
    face_outer_edge_refs(brep, face_sr).into_iter()
        .map(|(idx, fwd)| WireEdge { idx, forward: fwd }).collect()
}

fn face_inner_wires(brep: &BRep, face_sr: ShapeRef) -> Vec<Wire> {
    face_inner_edge_refs(brep, face_sr).into_iter()
        .map(|edges| Wire { edges: edges.into_iter().map(|(idx, fwd)| WireEdge { idx, forward: fwd }).collect() })
        .collect()
}

fn face_normal(brep: &BRep, face_sr: ShapeRef) -> DVec3 {
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return DVec3::ZERO };
    fd.surface.as_ref().map(|s| rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0)).unwrap_or(DVec3::ZERO)
}

fn face_outer_polygon_points(brep: &BRep, si: usize, shi: usize, fi: usize) -> Vec<DVec3> {
    // Resolve shell/face from TShape hierarchy (read-only, no persistent borrow)
    let TShape::Solid(sd) = &*brep.tshapes[si] else { return vec![] };
    let shell_sr = sd.shells.get(shi).copied().unwrap_or(ShapeRef::NULL);
    let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { return vec![] };
    let face_sr = shd.faces.get(fi).copied().unwrap_or(ShapeRef::NULL);
    let edges = face_outer_edge_refs(brep, face_sr);
    edges.iter().map(|&(ei, fwd)| vpoint(brep, if fwd { edge_start(brep, ei) } else { edge_end(brep, ei) })).collect()
}

/// Surface data plus hierarchy path for a face, extracted without holding persistent refs.
struct FaceMergeInfo {
    surface: Option<rcad_kernel::geom::Surface3>,
    sample_point: Option<DVec3>,
    uv_domain: Option<[f64; 4]>,
    internal_vertices: Vec<ShapeRef>,
    natural_restriction: bool,
    normal: DVec3,
    outer_wire_edges: Vec<(usize, bool)>,
    inner_wires: Vec<Vec<(usize, bool)>>,
}

fn extract_face_info(brep: &BRep, si: usize, shi: usize, fi: usize) -> Option<FaceMergeInfo> {
    let TShape::Solid(sd) = &*brep.tshapes[si] else { return None };
    let shell_sr = sd.shells.get(shi).copied()?;
    let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { return None };
    let face_sr = shd.faces.get(fi).copied()?;
    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return None };
    Some(FaceMergeInfo {
        surface: fd.surface.clone(),
        sample_point: fd.sample_point,
        uv_domain: fd.uv_domain,
        internal_vertices: fd.internal_vertices.clone(),
        natural_restriction: fd.natural_restriction,
        normal: fd.surface.as_ref().map(|s| rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0)).unwrap_or(DVec3::ZERO),
        outer_wire_edges: {
            let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else { return None };
            wd.edges.iter().map(|esr| (esr.index, esr.orientation == Orientation::Forward)).collect()
        },
        inner_wires: fd.inner_wires.iter().filter_map(|iw_sr| {
            let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else { return None };
            Some(iwd.edges.iter().map(|esr| (esr.index, esr.orientation == Orientation::Forward)).collect())
        }).collect(),
    })
}

// ===================================================================
// Main merge pass
// ===================================================================

fn unify_one_merge_pass_with_origins(brep: &mut BRep, face_origins: Option<&[FaceOrigin]>) -> bool {
    fn closure_score(brep: &BRep) -> usize {
        let report = crate::brep_check::validate_solid_closure(brep);
        report.issues.iter().map(|iss| match iss {
            crate::CheckIssue::SolidNotClosed { boundary_edge_count, .. } => *boundary_edge_count,
            _ => 1,
        }).sum()
    }

    fn surfaces_are_same_domain(
        brep: &BRep, si: usize, shi: usize, fi1: usize, fi2: usize,
    ) -> (Option<bool>, bool) {
        let ang_tol = tolerance::TOLERANCE_ANG_HEURISTIC_RAD;
        let lin_tol = tolerance::TOLERANCE_PARAM_LEGACY;

        let fi1 = extract_face_info(brep, si, shi, fi1);
        let fi2 = extract_face_info(brep, si, shi, fi2);
        let (fi1, fi2) = match (fi1, fi2) { (Some(a), Some(b)) => (a, b), _ => return (None, true) };

        let s1 = match &fi1.surface { Some(s) => s, None => return (None, true) };
        let s2 = match &fi2.surface { Some(s) => s, None => return (None, true) };

        use rcad_kernel::geom::Surface3;
        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let n1 = p1.normal.normalize_or_zero();
                let n2 = p2.normal.normalize_or_zero();
                if n1.length_squared() <= tolerance::TOLERANCE_VEC_SQ_MIN || n2.length_squared() <= tolerance::TOLERANCE_VEC_SQ_MIN { return (Some(false), true); }
                let cross = n1.cross(n2).length();
                if cross > ang_tol { return (Some(false), true); }
                (Some((p2.origin - p1.origin).dot(n1).abs() <= lin_tol), true)
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                if (c1.radius - c2.radius).abs() > lin_tol { return (Some(false), false); }
                let (a1, a2) = (c1.axis.normalize_or_zero(), c2.axis.normalize_or_zero());
                if a1.cross(a2).length() > ang_tol { return (Some(false), false); }
                (Some((c2.origin - c1.origin).cross(a1).length() <= lin_tol), false)
            }
            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                if (c1.radius - c2.radius).abs() > lin_tol { return (Some(false), false); }
                if (c1.half_angle_rad - c2.half_angle_rad).abs() > ang_tol { return (Some(false), false); }
                let (a1, a2) = (c1.axis.normalize_or_zero(), c2.axis.normalize_or_zero());
                if a1.cross(a2).length() > ang_tol { return (Some(false), false); }
                (Some((c1.apex - c2.apex).length() <= lin_tol), false)
            }
            (Surface3::Torus(t1), Surface3::Torus(t2)) => {
                if (t1.major_radius - t2.major_radius).abs() > lin_tol { return (Some(false), false); }
                if (t1.minor_radius - t2.minor_radius).abs() > lin_tol { return (Some(false), false); }
                let (a1, a2) = (t1.axis.normalize_or_zero(), t2.axis.normalize_or_zero());
                if a1.cross(a2).length() > ang_tol { return (Some(false), false); }
                (Some((t1.center - t2.center).length() <= lin_tol), false)
            }
            (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
                if (s1.radius - s2.radius).abs() > lin_tol { return (Some(false), false); }
                (Some((s1.center - s2.center).length() <= lin_tol), false)
            }
            (Surface3::BSpline(_), Surface3::Plane(_)) | (Surface3::Plane(_), Surface3::BSpline(_)) => (Some(false), false),
            (Surface3::BSpline(b1), Surface3::BSpline(b2)) => {
                if b1.degree_u != b2.degree_u || b1.degree_v != b2.degree_v { return (Some(false), false); }
                if b1.knots_u.len() != b2.knots_u.len() || b1.knots_v.len() != b2.knots_v.len() { return (Some(false), false); }
                if !b1.knots_u.iter().zip(b2.knots_u.iter()).all(|(k1,k2)| (k1-k2).abs()<=lin_tol) { return (Some(false), false); }
                if !b1.knots_v.iter().zip(b2.knots_v.iter()).all(|(k1,k2)| (k1-k2).abs()<=lin_tol) { return (Some(false), false); }
                if b1.control_points.len() != b2.control_points.len() { return (Some(false), false); }
                if !b1.control_points.iter().zip(b2.control_points.iter()).all(|(row1,row2)| {
                    row1.len()==row2.len() && row1.iter().zip(row2.iter()).all(|(cp1,cp2)| cp1.distance(*cp2)<=lin_tol)
                }) { return (Some(false), false); }
                if b1.weights.len() != b2.weights.len() { return (Some(false), false); }
                if !b1.weights.iter().zip(b2.weights.iter()).all(|(row1,row2)| {
                    row1.len()==row2.len() && row1.iter().zip(row2.iter()).all(|(w1,w2)| (w1-w2).abs()<=lin_tol)
                }) { return (Some(false), false); }
                (Some(true), false)
            }
            _ => (Some(false), false),
        }
    }

    fn flat_face_idx(brep: &BRep, si: usize, shi: usize, fi: usize) -> usize {
        let mut idx = 0usize;
        for s in 0..brep.tshapes.len() {
            let TShape::Solid(sd) = &*brep.tshapes[s] else { continue };
            for (i, shell_sr) in sd.shells.iter().enumerate() {
                let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
                if s < si || (s == si && i < shi) { idx += shd.faces.len(); }
            }
        }
        idx + fi
    }

    fn quantize_edge_point(p: DVec3) -> (i64, i64, i64) {
        let inv_tol = 1.0 / tolerance::TOLERANCE_PARAM_LEGACY.max(tolerance::TOLERANCE_ABS);
        ((p.x * inv_tol).round() as i64, (p.y * inv_tol).round() as i64, (p.z * inv_tol).round() as i64)
    }

    // Collect shell info without holding persistent borrows
    struct ShellEntry { si: usize, shi: usize, nfaces: usize }
    let mut shells: Vec<ShellEntry> = Vec::new();
    for si in 0..brep.tshapes.len() {
        let TShape::Solid(sd) = &*brep.tshapes[si] else { continue };
        for shi in 0..sd.shells.len() {
            let TShape::Shell(shd) = &*brep.tshapes[sd.shells[shi].index] else { continue };
            if shd.faces.len() >= 2 {
                shells.push(ShellEntry { si, shi, nfaces: shd.faces.len() });
            }
        }
    }

    for se in &shells {
        let si = se.si;
        let shi = se.shi;
        let nfaces = se.nfaces;

        // Read-only phase: build adjacency
        let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut geom_edge_to_faces: HashMap<((i64,i64,i64),(i64,i64,i64)), Vec<(usize,usize)>> = HashMap::new();

        for fi in 0..nfaces {
            let fi_data = match extract_face_info(brep, si, shi, fi) { Some(d) => d, None => continue };
            for &(ei, _) in &fi_data.outer_wire_edges {
                edge_to_faces.entry(ei).or_default().push(fi);
                let qs = quantize_edge_point(vpoint(brep, edge_start(brep, ei)));
                let qe = quantize_edge_point(vpoint(brep, edge_end(brep, ei)));
                let key = if qs <= qe { (qs, qe) } else { (qe, qs) };
                geom_edge_to_faces.entry(key).or_default().push((fi, ei));
            }
            for inner in &fi_data.inner_wires {
                for &(ei, _) in inner {
                    edge_to_faces.entry(ei).or_default().push(fi);
                    let qs = quantize_edge_point(vpoint(brep, edge_start(brep, ei)));
                    let qe = quantize_edge_point(vpoint(brep, edge_end(brep, ei)));
                    let key = if qs <= qe { (qs, qe) } else { (qe, qs) };
                    geom_edge_to_faces.entry(key).or_default().push((fi, ei));
                }
            }
        }

        let mut candidates: Vec<(usize, usize, usize, usize)> = edge_to_faces.iter()
            .filter_map(|(&ei, refs)| if refs.len() == 2 { Some((ei, ei, refs[0], refs[1])) } else { None })
            .collect();
        for face_edges in geom_edge_to_faces.values() {
            if face_edges.len() != 2 { continue; }
            let (fi1, ei1) = face_edges[0];
            let (fi2, ei2) = face_edges[1];
            if fi1 != fi2 && ei1 != ei2 { candidates.push((ei1, ei2, fi1, fi2)); }
        }
        candidates.sort_unstable();
        candidates.dedup();

        for &(edge_idx1, edge_idx2, fi1, fi2) in &candidates {
            if fi1 == fi2 { continue; }
            let d1 = match extract_face_info(brep, si, shi, fi1) { Some(d) => d, None => continue };
            let d2 = match extract_face_info(brep, si, shi, fi2) { Some(d) => d, None => continue };

            let (same_domain, is_planar) = surfaces_are_same_domain(brep, si, shi, fi1, fi2);

            if let Some(origins) = face_origins {
                let ff1 = flat_face_idx(brep, si, shi, fi1);
                let ff2 = flat_face_idx(brep, si, shi, fi2);
                if origins.get(ff1) != origins.get(ff2) { continue; }
            }

            // Should-merge decision (same logic as original)
            let get_pt = |fi: usize| -> Option<DVec3> {
                let d = if fi == fi1 { &d1 } else { &d2 };
                let (ei, fwd) = *d.outer_wire_edges.first()?;
                Some(vpoint(brep, if fwd { edge_start(brep, ei) } else { edge_end(brep, ei) }))
            };
            let outer_verts = |fi: usize| -> Option<Vec<DVec3>> {
                let d = if fi == fi1 { &d1 } else { &d2 };
                let mut out = Vec::new();
                for &(ei, fwd) in &d.outer_wire_edges {
                    out.push(vpoint(brep, if fwd { edge_start(brep, ei) } else { edge_end(brep, ei) }));
                }
                if out.is_empty() { None } else { Some(out) }
            };

            let mut should_merge = match same_domain {
                Some(false) => false,
                Some(true) => {
                    if is_planar {
                        match (get_pt(fi1), outer_verts(fi1), outer_verts(fi2)) {
                            (Some(pt1), Some(vs1), Some(vs2)) => {
                                let n = d1.normal.normalize();
                                vs1.iter().all(|p| (*p-pt1).dot(n).abs() <= tolerance::TOLERANCE_PLANE_DIST_RELAX)
                                    && vs2.iter().all(|p| (*p-pt1).dot(n).abs() <= tolerance::TOLERANCE_PLANE_DIST_RELAX)
                            }
                            _ => false,
                        }
                    } else { true }
                }
                None => {
                    let cross = d1.normal.cross(d2.normal).length();
                    if cross > tolerance::TOLERANCE_PARAM_LEGACY { false }
                    else if let (Some(pt1), Some(pt2)) = (get_pt(fi1), get_pt(fi2)) {
                        (pt2 - pt1).dot(d1.normal.normalize()).abs() <= tolerance::TOLERANCE_PARAM_LEGACY
                    } else { false }
                }
            };

            if should_merge && edge_idx1 == edge_idx2 {
                if !validate_shared_edge_continuity(brep, si, shi, fi1, fi2, edge_idx1) {
                    should_merge = false;
                }
            }
            if should_merge {
                let uv_ok = if is_planar && same_domain == Some(true) { true }
                    else { validate_uv_regions_compatible(brep, si, shi, fi1, fi2) };
                if !uv_ok { should_merge = false; }
            }
            if !should_merge { continue; }
            if !is_planar {
                if d1.outer_wire_edges.len() + d2.outer_wire_edges.len() > 650 { continue; }
            }

            // Build merged wire (using old WireEdge for splice ops — still valid type)
            let wire1: Vec<WireEdge> = d1.outer_wire_edges.iter().map(|&(ei,fwd)| WireEdge { idx: ei, forward: fwd }).collect();
            let wire2: Vec<WireEdge> = d2.outer_wire_edges.iter().map(|&(ei,fwd)| WireEdge { idx: ei, forward: fwd }).collect();

            if let Some(merged_wire) = splice_wires(&wire1, edge_idx1, &wire2, edge_idx2) {
                let merged_wire = cleanup_merged_wire_edges(brep, &merged_wire);
                let mut all_inner: Vec<Wire> = d1.inner_wires.iter().map(|edges| {
                    Wire { edges: edges.iter().map(|&(ei,fwd)| WireEdge { idx: ei, forward: fwd }).collect() }
                }).collect();
                for edges in &d2.inner_wires {
                    all_inner.push(Wire { edges: edges.iter().map(|&(ei,fwd)| WireEdge { idx: ei, forward: fwd }).collect() });
                }

                let (outer_raw, extracted) = extract_inner_loops_from_wire(brep, &merged_wire);
                let outer_clean = if extracted.is_empty() { outer_raw } else { cleanup_merged_wire_edges(brep, &outer_raw) };
                all_inner.extend(extracted);

                // Planar area guard
                if is_planar {
                    let nunit = d1.normal.normalize_or_zero();
                    let poly1 = face_outer_polygon_points(brep, si, shi, fi1);
                    let poly2 = face_outer_polygon_points(brep, si, shi, fi2);
                    let a1 = newell_polygon_abs_area(&poly1, nunit);
                    let a2 = newell_polygon_abs_area(&poly2, nunit);
                    let mut poly_m: Vec<DVec3> = Vec::new();
                    for we in &outer_clean {
                        if let Some((u, _)) = oriented_edge_vertices(brep, *we) { poly_m.push(vpoint(brep, u)); }
                    }
                    let am = newell_polygon_abs_area(&poly_m, nunit);
                    let sum = a1 + a2;
                    let tol = tolerance::TOLERANCE_AREA_REL * sum.max(am).max(1.0) + tolerance::TOLERANCE_ABS;
                    if am > sum + tol { continue; }
                }

                // --- Write phase: mutate clone, then assign back ---
                let (keep_idx, remove_idx) = if fi1 < fi2 { (fi1, fi2) } else { (fi2, fi1) };
                let surface = d1.surface.clone();
                let sample_point = d1.sample_point;
                let uv_domain = d1.uv_domain;
                let internal_vertices = d1.internal_vertices.clone();
                let natural_restriction = d1.natural_restriction;
                let current_score = closure_score(brep);

                // Clone, mutate, assign back (no persistent borrows on brep)
                let mut candidate = brep.clone();
                let merged_edge_refs: Vec<ShapeRef> = outer_clean.iter().map(|we| {
                    ShapeRef::synthetic_with_orientation(we.idx, if we.forward { Orientation::Forward } else { Orientation::Reversed })
                }).collect();
                let merged_wire_sr = candidate.add_twire(merged_edge_refs);
                let merged_inner_wire_srs: Vec<ShapeRef> = all_inner.iter().map(|w| {
                    let refs: Vec<ShapeRef> = w.edges.iter().map(|we| {
                        ShapeRef::synthetic_with_orientation(we.idx, if we.forward { Orientation::Forward } else { Orientation::Reversed })
                    }).collect();
                    candidate.add_twire(refs)
                }).collect();
                let merged_face_sr = candidate.add_tface(surface, merged_wire_sr, merged_inner_wire_srs,
                    sample_point, uv_domain, internal_vertices, natural_restriction);

                // Modify candidate's shell face list
                let TShape::Solid(sd_c) = &*candidate.tshapes[si] else { continue };
                let TShape::Shell(shd_c) = &*candidate.tshapes[sd_c.shells[shi].index] else { continue };
                let mut new_faces: Vec<ShapeRef> = Vec::with_capacity(shd_c.faces.len() - 1);
                for (i, &fsr) in shd_c.faces.iter().enumerate() {
                    if i == remove_idx { continue; }
                    if i == keep_idx { new_faces.push(merged_face_sr); } else { new_faces.push(fsr); }
                }
                candidate.shell_mut(sd_c.shells[shi]).faces = new_faces;

                let candidate_score = closure_score(&candidate);
                if candidate_score > current_score { continue; }
                *brep = candidate;
                return true;
            }
        }
    }

    false
}

// ===================================================================
// Wire splicing helpers (unchanged logic, TShape-based vertex access)
// ===================================================================

fn splice_wires(wire_a: &[WireEdge], shared_idx_a: usize, wire_b: &[WireEdge], shared_idx_b: usize) -> Option<Vec<WireEdge>> {
    let pos_a = wire_a.iter().position(|we| we.idx == shared_idx_a)?;
    let pos_b = wire_b.iter().position(|we| we.idx == shared_idx_b)?;
    let b_edges: Vec<WireEdge> = (1..wire_b.len()).map(|i| wire_b[(pos_b + i) % wire_b.len()]).collect();
    let mut merged = Vec::with_capacity(wire_a.len() - 1 + b_edges.len());
    merged.extend_from_slice(&wire_a[..pos_a]);
    merged.extend(b_edges);
    merged.extend_from_slice(&wire_a[pos_a + 1..]);
    if merged.len() < 3 { None } else { Some(merged) }
}

pub(crate) fn oriented_edge_vertices(brep: &BRep, we: WireEdge) -> Option<(usize, usize)> {
    let (s, e) = (edge_start(brep, we.idx), edge_end(brep, we.idx));
    Some(if we.forward { (s, e) } else { (e, s) })
}

fn find_existing_edge_between_vertices(brep: &BRep, from: usize, to: usize) -> Option<WireEdge> {
    for (idx, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Edge(ed) = &**ts {
            if ed.first.index == from && ed.last.index == to { return Some(WireEdge::fwd(idx)); }
            if ed.first.index == to && ed.last.index == from { return Some(WireEdge::rev(idx)); }
        }
    }
    None
}

fn points_are_collinear_forward(a: DVec3, b: DVec3, c: DVec3) -> bool {
    let (ab, bc) = (b - a, c - b);
    let (ab_len, bc_len) = (ab.length(), bc.length());
    if ab_len <= tolerance::TOLERANCE_LEN_MIN || bc_len <= tolerance::TOLERANCE_LEN_MIN { return false; }
    ab.cross(bc).length() <= tolerance::TOLERANCE_ABS * (ab_len + bc_len) && ab.dot(bc) > 0.0
}

fn collapse_collinear_segments_with_existing_bridge(brep: &BRep, wire: &[WireEdge]) -> Option<Vec<WireEdge>> {
    let mut out = wire.to_vec();
    if out.len() < 4 { return None; }
    loop {
        if out.len() < 4 { break; }
        let (mut changed, n) = (false, out.len());
        for i in 0..n {
            let j = (i + 1) % n;
            let (u, v1) = oriented_edge_vertices(brep, out[i])?;
            let (v2, w) = oriented_edge_vertices(brep, out[j])?;
            if v1 != v2 || u == w { continue; }
            let (pu, pv, pw) = (vpoint(brep, u), vpoint(brep, v1), vpoint(brep, w));
            if !points_are_collinear_forward(pu, pv, pw) { continue; }
            let bridge = match find_existing_edge_between_vertices(brep, u, w) {
                Some(e) if e.idx != out[i].idx && e.idx != out[j].idx => e, _ => continue,
            };
            if i + 1 < n { out.splice(i..=i+1, [bridge]); }
            else { out.pop(); out.remove(0); out.insert(0, bridge); }
            changed = true; break;
        }
        if !changed { break; }
    }
    if out.len() >= 3 { Some(out) } else { None }
}

fn wire_is_closed_and_connected(brep: &BRep, wire: &[WireEdge]) -> bool {
    if wire.len() < 3 { return false; }
    let (first_start, mut prev) = match oriented_edge_vertices(brep, wire[0]) { Some(v) => v, None => return false };
    for we in &wire[1..] {
        let (s, e) = match oriented_edge_vertices(brep, *we) { Some(v) => v, None => return false };
        if s != prev { return false; }
        prev = e;
    }
    prev == first_start
}

fn reorder_wire_into_connected_loop(brep: &BRep, wire: &[WireEdge]) -> Option<Vec<WireEdge>> {
    if wire.is_empty() { return None; }
    let mut unused: Vec<WireEdge> = wire.to_vec();
    let first = unused.remove(0);
    let mut out = vec![first];
    let (_, mut cur_end) = oriented_edge_vertices(brep, first)?;
    while !unused.is_empty() {
        let mut found = None;
        for (i, we) in unused.iter().enumerate() {
            let (s, e) = oriented_edge_vertices(brep, *we)?;
            if s == cur_end { found = Some((i, false)); break; }
            if e == cur_end { found = Some((i, true)); break; }
        }
        let (i, flip) = found?;
        let mut next = unused.remove(i);
        if flip { next.forward = !next.forward; }
        let (_, ne) = oriented_edge_vertices(brep, next)?;
        out.push(next); cur_end = ne;
    }
    if wire_is_closed_and_connected(brep, &out) { Some(out) } else { None }
}

fn cancel_duplicate_segments_by_parity(brep: &BRep, wire: &[WireEdge]) -> Option<Vec<WireEdge>> {
    let mut groups: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (i, &we) in wire.iter().enumerate() {
        let (u, v) = oriented_edge_vertices(brep, we)?;
        let key = if u <= v { (u, v) } else { (v, u) };
        groups.entry(key).or_default().push(i);
    }
    let mut keep = vec![true; wire.len()];
    for idxs in groups.values() {
        if idxs.len() >= 2 {
            for idx in idxs.iter().take((idxs.len() / 2) * 2) { keep[*idx] = false; }
        }
    }
    let out: Vec<WireEdge> = wire.iter().enumerate().filter_map(|(i, &we)| if keep[i] { Some(we) } else { None }).collect();
    if out.len() >= 3 { Some(out) } else { None }
}

fn extract_inner_loops_from_wire(brep: &BRep, wire: &[WireEdge]) -> (Vec<WireEdge>, Vec<Wire>) {
    let verts: Vec<usize> = wire.iter().filter_map(|&we| oriented_edge_vertices(brep, we).map(|(u,_)| u)).collect();
    let mut seen: HashMap<usize, usize> = HashMap::new();
    let mut split: Option<(usize, usize)> = None;
    for (i, &v) in verts.iter().enumerate() {
        if let Some(&f) = seen.get(&v) { split = Some((f, i)); break; }
        seen.insert(v, i);
    }
    let Some((start, end)) = split else { return (wire.to_vec(), vec![]) };
    let inner: Vec<WireEdge> = wire[start..end].to_vec();
    let outer: Vec<WireEdge> = wire[..start].iter().chain(wire[end..].iter()).copied().collect();
    if inner.len() < 3 || outer.len() < 3 { return (wire.to_vec(), vec![]); }
    let inner_wire = Wire { edges: inner };
    let (final_outer, mut more) = extract_inner_loops_from_wire(brep, &outer);
    more.push(inner_wire);
    (final_outer, more)
}

fn cleanup_merged_wire_edges(brep: &BRep, wire: &[WireEdge]) -> Vec<WireEdge> {
    if wire.len() < 4 { return wire.to_vec(); }
    let mut cleaned: Vec<WireEdge> = Vec::new();
    for &we in wire {
        let (u, v) = oriented_edge_vertices(brep, we).unwrap_or((0,0));
        if let Some(&last) = cleaned.last() {
            let (lu, lv) = oriented_edge_vertices(brep, last).unwrap_or((0,0));
            if (lu==u && lv==v) || (lu==v && lv==u) { cleaned.pop(); continue; }
        }
        cleaned.push(we);
    }
    while cleaned.len() >= 2 {
        let (fu, fv) = oriented_edge_vertices(brep, cleaned[0]).unwrap_or((0,0));
        let (lu, lv) = oriented_edge_vertices(brep, *cleaned.last().unwrap()).unwrap_or((0,0));
        if !((fu==lu && fv==lv) || (fu==lv && fv==lu)) { break; }
        cleaned.remove(0); cleaned.pop();
    }
    let stage = if wire_is_closed_and_connected(brep, &cleaned) { Some(cleaned) }
        else if let Some(c) = cancel_duplicate_segments_by_parity(brep, &cleaned) { reorder_wire_into_connected_loop(brep, &c) }
        else { None };
    let Some(mut out) = stage else { return wire.to_vec() };
    if let Some(c) = collapse_collinear_segments_with_existing_bridge(brep, &out)
        && let Some(r) = reorder_wire_into_connected_loop(brep, &c)
        && wire_is_closed_and_connected(brep, &r) { out = r; }
    out
}

// ===================================================================
// remove_internal_faces — TShape-based rewrite
// ===================================================================

pub fn remove_internal_faces(brep: &BRep) -> (BRep, usize) {
    fn same_domain(brep: &BRep, si: usize, shi: usize, fi1: usize, fi2: usize) -> Option<bool> {
        let ang_tol = tolerance::TOLERANCE_ANG_HEURISTIC_RAD;
        let lin_tol = tolerance::TOLERANCE_PARAM_LEGACY;
        let fi1 = extract_face_info(brep, si, shi, fi1)?;
        let fi2 = extract_face_info(brep, si, shi, fi2)?;
        let (s1, s2) = (fi1.surface.as_ref()?, fi2.surface.as_ref()?);
        use rcad_kernel::geom::Surface3;
        Some(match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let (n1,n2)=(p1.normal.normalize_or_zero(),p2.normal.normalize_or_zero());
                if n1.length_squared()<=tolerance::TOLERANCE_VEC_SQ_MIN||n2.length_squared()<=tolerance::TOLERANCE_VEC_SQ_MIN{false}
                else { n1.cross(n2).length()<=ang_tol && (p2.origin-p1.origin).dot(n1).abs()<=lin_tol }
            }
            (Surface3::Cylinder(c1),Surface3::Cylinder(c2)) => {
                if (c1.radius-c2.radius).abs()>lin_tol{false}
                else{let(a1,a2)=(c1.axis.normalize_or_zero(),c2.axis.normalize_or_zero());
                a1.cross(a2).length()<=ang_tol&&(c2.origin-c1.origin).cross(a1).length()<=lin_tol}
            }
            (Surface3::Cone(c1),Surface3::Cone(c2)) => {
                let (a1c,a2c)=(c1.axis.normalize_or_zero(),c2.axis.normalize_or_zero());
                (c1.radius-c2.radius).abs()<=lin_tol&&(c1.half_angle_rad-c2.half_angle_rad).abs()<=ang_tol
                    &&a1c.cross(a2c).length()<=ang_tol&&(c1.apex-c2.apex).length()<=lin_tol
            }
            (Surface3::Torus(t1),Surface3::Torus(t2)) => {
                (t1.major_radius-t2.major_radius).abs()<=lin_tol&&(t1.minor_radius-t2.minor_radius).abs()<=lin_tol
                    &&t1.axis.normalize_or_zero().cross(t2.axis.normalize_or_zero()).length()<=ang_tol
                    &&(t1.center-t2.center).length()<=lin_tol
            }
            (Surface3::Sphere(s1),Surface3::Sphere(s2)) => {
                (s1.radius-s2.radius).abs()<=lin_tol&&(s1.center-s2.center).length()<=lin_tol
            }
            (Surface3::BSpline(b1),Surface3::BSpline(b2)) => {
                b1.degree_u==b2.degree_u&&b1.degree_v==b2.degree_v
                    &&b1.knots_u.len()==b2.knots_u.len()&&b1.knots_v.len()==b2.knots_v.len()
                    &&b1.knots_u.iter().zip(&b2.knots_u).all(|(k1,k2)|(k1-k2).abs()<=lin_tol)
                    &&b1.knots_v.iter().zip(&b2.knots_v).all(|(k1,k2)|(k1-k2).abs()<=lin_tol)
                    &&b1.control_points.len()==b2.control_points.len()
                    &&b1.control_points.iter().zip(&b2.control_points).all(|(r1,r2)|
                        r1.len()==r2.len()&&r1.iter().zip(r2.iter()).all(|(c1,c2)|c1.distance(*c2)<=lin_tol))
                    &&b1.weights.len()==b2.weights.len()
                    &&b1.weights.iter().zip(&b2.weights).all(|(r1,r2)|
                        r1.len()==r2.len()&&r1.iter().zip(r2.iter()).all(|(w1,w2)|(w1-w2).abs()<=lin_tol))
            }
            _=>false,
        })
    }

    let mut out = brep.clone();
    let mut total_removed = 0usize;
    struct ShellIdx { si: usize, shi: usize }
    let mut shell_list: Vec<ShellIdx> = Vec::new();
    for si in 0..out.tshapes.len() {
        let TShape::Solid(sd) = &*out.tshapes[si] else { continue };
        for shi in 0..sd.shells.len() {
            shell_list.push(ShellIdx { si, shi });
        }
    }

    for sii in &shell_list {
        let si = sii.si;
        loop {
            let TShape::Solid(sd) = &*out.tshapes[si] else { break };
            let TShape::Shell(shd) = &*out.tshapes[sd.shells[sii.shi].index] else { break };
            let nfaces = shd.faces.len();
            let mut removed = None;
            'outer: for fi in 0..nfaces {
                for fj in (fi+1)..nfaces {
                    let fi_d = match extract_face_info(&out, si, sii.shi, fi) { Some(d)=>d, None=>continue };
                    let fj_d = match extract_face_info(&out, si, sii.shi, fj) { Some(d)=>d, None=>continue };
                    if fi_d.normal==DVec3::ZERO || fj_d.normal==DVec3::ZERO { continue; }
                    let cross = fi_d.normal.cross(fj_d.normal).length();
                    let dot = fi_d.normal.normalize().dot(fj_d.normal.normalize());
                    if cross > tolerance::TOLERANCE_PARAM_LEGACY || dot.abs() < tolerance::TOLERANCE_DOT_NEARLY_PARALLEL { continue; }

                    let sdg = same_domain(&out, si, sii.shi, fi, fj);
                    let pi = vpoint(&out, edge_start(&out, fi_d.outer_wire_edges.first().map(|&(ei,_)|ei).unwrap_or(0)));
                    let pj = vpoint(&out, edge_start(&out, fj_d.outer_wire_edges.first().map(|&(ei,_)|ei).unwrap_or(0)));
                    let plane_fb = (pj-pi).dot(fi_d.normal.normalize()).abs() <= tolerance::TOLERANCE_PLANE_DIST_RELAX;
                    if !matches!(sdg, Some(true)) && !plane_fb { continue; }

                    let eis: HashSet<usize> = fi_d.outer_wire_edges.iter().chain(fi_d.inner_wires.iter().flatten()).map(|&(ei,_)|ei).collect();
                    let ejs: HashSet<usize> = fj_d.outer_wire_edges.iter().chain(fj_d.inner_wires.iter().flatten()).map(|&(ei,_)|ei).collect();
                    let overlap = eis.intersection(&ejs).count();
                    let min_edges = eis.len().min(ejs.len()).max(1);
                    let same_or_contained = overlap == min_edges || (matches!(sdg, Some(true)) && overlap as f64 / min_edges as f64 >= 0.60);
                    if !same_or_contained { continue; }

                    let nidot = fi_d.normal.normalize_or_zero().dot(fj_d.normal.normalize_or_zero());
                    let opposite = nidot < -0.99;
                    if !opposite { continue; }
                    if overlap != eis.len() || overlap != ejs.len() { continue; }

                    removed = Some(fj);
                    break 'outer;
                }
            }
            if let Some(idx) = removed {
                let TShape::Solid(sd_m) = &*out.tshapes[si] else { break };
                let shd_m = out.shell_mut(sd_m.shells[sii.shi]);
                shd_m.faces.remove(idx);
                total_removed += 1;
            } else { break; }
        }
    }

    (out, total_removed)
}
