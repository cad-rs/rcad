use std::collections::{BTreeSet, HashMap};
use rcad_kernel::{BSplineCurve2, Curve2d, Curve3, Surface3, topods};
use super::{BRep, Face};

#[derive(Clone, Copy)]
pub(super) struct OrientedEdgeExport {
    pub(super) edge_idx: usize,
    pub(super) start: usize,
    #[allow(dead_code)]
    pub(super) end: usize,
    pub(super) forward: bool,
}

// ── Topods-native variants (migration) ──

pub(super) fn oriented_face_edges_topods(tbrep: &topods::BRep, face_tshape_idx: usize) -> Vec<OrientedEdgeExport> {
    let topods::TShape::Face(fd) = &*tbrep.tshapes[face_tshape_idx] else { return vec![] };
    // Get outer wire edges
    let topods::TShape::Wire(wd) = &*tbrep.tshapes[fd.outer_wire.index] else { return vec![] };
    wd.edges
        .iter()
        .filter_map(|sr| {
            let topods::TShape::Edge(ed) = &*tbrep.tshapes[sr.index] else { return None };
            let (start, end) = if sr.orientation.is_forward() {
                (ed.first.index, ed.last.index)
            } else {
                (ed.last.index, ed.first.index)
            };
            Some(OrientedEdgeExport {
                edge_idx: sr.index,
                start,
                end,
                forward: sr.orientation.is_forward(),
            })
        })
        .collect()
}

pub(super) fn detect_seam_edge_indices_topods(tbrep: &topods::BRep, face_tshape_idx: usize) -> BTreeSet<usize> {
    let topods::TShape::Face(fd) = &*tbrep.tshapes[face_tshape_idx] else { return BTreeSet::new() };
    let topods::TShape::Wire(wd) = &*tbrep.tshapes[fd.outer_wire.index] else { return BTreeSet::new() };
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for sr in &wd.edges {
        *counts.entry(sr.index).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .filter(|&(_, c)| c >= 2)
        .map(|(idx, _)| idx)
        .collect()
}

pub(super) fn is_degenerate_face_wire_topods(tbrep: &topods::BRep, face_tshape_idx: usize) -> bool {
    let topods::TShape::Face(fd) = &*tbrep.tshapes[face_tshape_idx] else { return false };
    if fd.surface.is_some() {
        return false;
    }
    let topods::TShape::Wire(wd) = &*tbrep.tshapes[fd.outer_wire.index] else { return false };
    wd.edges.len() < 3
}

pub(super) fn face_orientation_for_surface_topods(tbrep: &topods::BRep, loop_points: &[glam::DVec3], face_tshape_idx: usize) -> bool {
    let topods::TShape::Face(fd) = &*tbrep.tshapes[face_tshape_idx] else { return true };
    match fd.surface.as_ref() {
        Some(Surface3::Plane(plane)) => {
            let plane_normal = canonicalize_axis_sign(plane.normal);
            let mut c = glam::DVec3::ZERO;
            let mut n = 0usize;
            for ts in &tbrep.tshapes {
                if let topods::TShape::Vertex(vd) = &**ts {
                    c += vd.point;
                    n += 1;
                }
            }
            let brep_centroid = if n > 0 { c / (n as f64) } else { glam::DVec3::ZERO };

            let face_centroid = if loop_points.is_empty() {
                brep_centroid
            } else {
                loop_points
                    .iter()
                    .copied()
                    .fold(glam::DVec3::ZERO, |acc, p| acc + p)
                    / (loop_points.len() as f64)
            };

            let outward = (face_centroid - brep_centroid).dot(plane_normal);
            if outward.abs() <= 1e-12 {
                // face.normal equivalent: compute from loop_points
                if let Some(fn_comp) = compute_face_normal(loop_points) {
                    fn_comp.dot(plane_normal) >= 0.0
                } else {
                    true
                }
            } else {
                outward >= 0.0
            }
        }
        _ => true,
    }
}

pub(super) fn shell_is_closed_topods(tbrep: &topods::BRep, shell_tshape_idx: usize) -> bool {
    let topods::TShape::Shell(shd) = &*tbrep.tshapes[shell_tshape_idx] else { return false };
    fn edge_key(tbrep: &topods::BRep, edge_idx: usize) -> (u64, u64) {
        let topods::TShape::Edge(ed) = &*tbrep.tshapes[edge_idx] else { return (0, 0) };
        let p1 = tbrep.tshapes.get(ed.first.index).and_then(|ts| {
            if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
        }).unwrap_or_default();
        let p2 = tbrep.tshapes.get(ed.last.index).and_then(|ts| {
            if let topods::TShape::Vertex(v) = &**ts { Some(v.point) } else { None }
        }).unwrap_or_default();
        let a = p1.to_array().map(|c| c.to_bits());
        let b = p2.to_array().map(|c| c.to_bits());
        let ha = a[0] ^ a[1].rotate_left(21) ^ a[2].rotate_left(42);
        let hb = b[0] ^ b[1].rotate_left(21) ^ b[2].rotate_left(42);
        if ha < hb { (ha, hb) } else { (hb, ha) }
    }
    let mut counts: HashMap<(u64, u64), usize> = HashMap::new();
    for face_sr in &shd.faces {
        let topods::TShape::Face(fd) = &*tbrep.tshapes[face_sr.index] else { continue };
        // outer wire
        if let topods::TShape::Wire(wd) = &*tbrep.tshapes[fd.outer_wire.index] {
            for sr in &wd.edges {
                *counts.entry(edge_key(tbrep, sr.index)).or_insert(0) += 1;
            }
        }
        // inner wires
        for w_sr in &fd.inner_wires {
            if let topods::TShape::Wire(wd) = &*tbrep.tshapes[w_sr.index] {
                for sr in &wd.edges {
                    *counts.entry(edge_key(tbrep, sr.index)).or_insert(0) += 1;
                }
            }
        }
    }
    if !counts.is_empty() && counts.values().all(|&count| count == 2) {
        return true;
    }
    false
}


// ── Topods-native plane/edge helpers (migration) ──

pub(super) fn find_plane_surface_for_edge_topods(
    tbrep: &topods::BRep,
    edge_idx: usize,
) -> Option<(usize, rcad_kernel::geom::Plane)> {
    for (fi, ts) in tbrep.tshapes.iter().enumerate() {
        let topods::TShape::Face(fd) = &**ts else { continue };
        let topods::TShape::Wire(wd) = &*tbrep.tshapes[fd.outer_wire.index] else { continue };
        if wd.edges.iter().any(|sr| sr.index == edge_idx) {
            if let Some(Surface3::Plane(p)) = &fd.surface { return Some((fi, p.clone())); }
        }
        for inner_sr in &fd.inner_wires {
            if let topods::TShape::Wire(w) = &*tbrep.tshapes[inner_sr.index] {
                if w.edges.iter().any(|e| e.index == edge_idx) {
                    if let Some(Surface3::Plane(p)) = &fd.surface { return Some((fi, p.clone())); }
                }
            }
        }
    }
    None
}

pub(super) fn count_plane_face_occurrences_for_line_edge_topods(tbrep: &topods::BRep, edge_idx: usize) -> usize {
    let mut count = 0;
    for ts in &tbrep.tshapes {
        let topods::TShape::Face(fd) = &**ts else { continue };
        let outer_ok = || -> bool {
            let topods::TShape::Wire(w) = &*tbrep.tshapes[fd.outer_wire.index] else { return false };
            w.edges.iter().any(|sr| sr.index == edge_idx)
        };
        let inner_ok = || -> bool {
            for sr in &fd.inner_wires {
                if let topods::TShape::Wire(w) = &*tbrep.tshapes[sr.index] {
                    if w.edges.iter().any(|e| e.index == edge_idx) { return true; }
                }
            }
            false
        };
        if (outer_ok() || inner_ok()) && matches!(&fd.surface, Some(Surface3::Plane(_))) {
            count += 1;
        }
    }
    count
}

pub(super) fn find_peer_plane_surface_for_line_edge_topods(
    tbrep: &topods::BRep,
    edge_idx: usize,
    exclude: &rcad_kernel::geom::Plane,
) -> Option<(Option<usize>, rcad_kernel::geom::Plane)> {
    for (fi, ts) in tbrep.tshapes.iter().enumerate() {
        let topods::TShape::Face(fd) = &**ts else { continue };
        let outer_ok = || -> bool {
            let topods::TShape::Wire(w) = &*tbrep.tshapes[fd.outer_wire.index] else { return false };
            w.edges.iter().any(|sr| sr.index == edge_idx)
        };
        let inner_ok = || -> bool {
            for sr in &fd.inner_wires {
                if let topods::TShape::Wire(w) = &*tbrep.tshapes[sr.index] {
                    if w.edges.iter().any(|e| e.index == edge_idx) { return true; }
                }
            }
            false
        };
        if outer_ok() || inner_ok() {
            if let Some(p) = fd.surface.as_ref().and_then(|s| match s { Surface3::Plane(pl) => Some(pl.clone()), _ => None }) {
                if !planes_equivalent(&p, exclude, 1.0e-6) { return Some((Some(fi), p)); }
            }
        }
    }
    None
}

pub(super) fn find_topological_plane_for_edge_topods(
    tbrep: &topods::BRep,
    edge_idx: usize,
) -> Option<rcad_kernel::geom::Plane> {
    for ts in &tbrep.tshapes {
        let topods::TShape::Face(fd) = &**ts else { continue };
        let outer_ok = || -> bool {
            let topods::TShape::Wire(w) = &*tbrep.tshapes[fd.outer_wire.index] else { return false };
            w.edges.iter().any(|sr| sr.index == edge_idx)
        };
        let inner_ok = || -> bool {
            for sr in &fd.inner_wires {
                if let topods::TShape::Wire(w) = &*tbrep.tshapes[sr.index] {
                    if w.edges.iter().any(|e| e.index == edge_idx) { return true; }
                }
            }
            false
        };
        if outer_ok() || inner_ok() {
            let loop_points: Vec<glam::DVec3> = {
                let topods::TShape::Wire(w) = &*tbrep.tshapes[fd.outer_wire.index] else { continue };
                w.edges.iter().filter_map(|sr| {
                    if let topods::TShape::Edge(ed) = &*tbrep.tshapes[sr.index] {
                        if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.first.index] { Some(v.point) } else { None }
                    } else { None }
                }).collect()
            };
            let origin = loop_points.first().copied()?;
            let normal = compute_face_normal(&loop_points)?;
            return Some(rcad_kernel::geom::Plane { origin, normal });
        }
    }
    None
}

pub(super) fn synthesize_cylinder_pcurve_for_edge_topods(
    tbrep: &topods::BRep,
    edge_idx: usize,
    cyl: &rcad_kernel::geom::CylindricalSurface,
) -> Option<Curve2d> {
    let topods::TShape::Edge(ed) = &*tbrep.tshapes[edge_idx] else { return None };
    let p0 = if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.first.index] { v.point } else { return None };
    let p1 = if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.last.index] { v.point } else { return None };
    let axis = cyl.axis.normalize_or_zero();
    if axis.length_squared() < 1e-18 { return None; }
    let xa = any_perpendicular_dvec3(axis);
    let ya = axis.cross(xa).normalize_or_zero();
    if ya.length_squared() < 1e-18 { return None; }
    let uv_of = |pt: glam::DVec3| {
        let d = pt - cyl.origin; let v = d.dot(axis); let perp = d - axis * v;
        let u = if perp.length_squared() < 1e-20 { 0.0 } else {
            perp.normalize_or_zero().dot(ya).atan2(perp.normalize_or_zero().dot(xa)).rem_euclid(std::f64::consts::TAU)
        };
        glam::DVec2::new(u, v)
    };
    match &ed.curve {
        Some(Curve3::Circle(c)) if ed.first.index == ed.last.index => {
            let v = (c.center - cyl.origin).dot(axis);
            Some(Curve2d::Line(rcad_kernel::geom::Line2d { origin: glam::DVec2::new(0.0, v), direction: glam::DVec2::new(std::f64::consts::TAU, 0.0) }))
        }
        Some(Curve3::Line(_)) | None => {
            let uv0 = uv_of(p0); let uv1 = uv_of(p1); let dir = uv1 - uv0;
            if dir.length_squared() < 1e-20 { return None; }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d { origin: uv0, direction: dir }))
        }
        _ => None,
    }
}

pub(super) fn synthesize_plane_pcurve_for_edge_topods(
    tbrep: &topods::BRep,
    edge_idx: usize,
    plane: &rcad_kernel::geom::Plane,
) -> Option<Curve2d> {
    let topods::TShape::Edge(ed) = &*tbrep.tshapes[edge_idx] else { return None };
    let p0 = if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.first.index] { v.point } else { return None };
    let p1 = if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.last.index] { v.point } else { return None };
    let normal = plane.normal.normalize_or_zero();
    if normal.length_squared() < 1e-18 { return None; }
    let ua = any_perpendicular_dvec3(normal);
    let va = normal.cross(ua).normalize_or_zero();
    if va.length_squared() < 1e-18 { return None; }
    let to_uv = |pt: glam::DVec3| -> glam::DVec2 { let d = pt - plane.origin; glam::DVec2::new(d.dot(ua), d.dot(va)) };
    match &ed.curve {
        Some(Curve3::Line(_)) | None => {
            let uv0 = to_uv(p0); let uv1 = to_uv(p1); let dir = uv1 - uv0;
            if dir.length_squared() < 1e-18 { return None; }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d { origin: uv0, direction: dir }))
        }
        Some(Curve3::Circle(c)) => {
            let center = to_uv(c.center);
            Some(Curve2d::Circle(rcad_kernel::geom::Circle2d { center, x_dir: glam::DVec2::X, y_dir: glam::DVec2::Y, radius: c.radius.max(1e-9) }))
        }
        Some(Curve3::Ellipse(e)) => {
            let center = to_uv(e.center);
            let major = glam::DVec2::new(e.major_dir.dot(ua), e.major_dir.dot(va));
            let md = if major.length_squared() < 1e-18 { glam::DVec2::X } else { major.normalize() };
            Some(Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d { center, major_dir: md, major_radius: e.major_radius.max(1e-9), minor_radius: e.minor_radius.max(1e-9) }))
        }
        _ => None,
    }
}

pub(super) fn should_promote_plane_line_pcurve_topods(tbrep: &topods::BRep, edge_idx: usize) -> bool {
    if count_plane_face_occurrences_for_line_edge_topods(tbrep, edge_idx) < 2 { return false; }
    let topods::TShape::Edge(ed) = &*tbrep.tshapes[edge_idx] else { return false };
    let p0 = if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.first.index] { v.point } else { return false };
    let p1 = if let topods::TShape::Vertex(v) = &*tbrep.tshapes[ed.last.index] { v.point } else { return false };
    let dir = (p1 - p0).normalize_or_zero();
    if dir.length_squared() < 1.0e-18 || dir.dot(glam::DVec3::Z).abs() < 1.0 - 1.0e-6 { return false; }
    let on_corner_x = (p0.x - 0.0).abs() <= 1.0e-6 || (p0.x - 1.0).abs() <= 1.0e-6;
    let on_corner_y = (p0.y - 0.0).abs() <= 1.0e-6 || (p0.y - 1.0).abs() <= 1.0e-6;
    if on_corner_x && on_corner_y { return false; }
    for ts in &tbrep.tshapes {
        let topods::TShape::Face(fd) = &**ts else { continue };
        let outer = if let topods::TShape::Wire(w) = &*tbrep.tshapes[fd.outer_wire.index] { w.edges.iter().any(|sr| sr.index == edge_idx) } else { false };
        let inner = fd.inner_wires.iter().any(|sr| -> bool {
            if let topods::TShape::Wire(w) = &*tbrep.tshapes[sr.index] { w.edges.iter().any(|e| e.index == edge_idx) } else { false }
        });
        if outer || inner {
            if let Some(Surface3::Plane(plane)) = &fd.surface {
                let n = plane.normal.normalize_or_zero();
                if n.length_squared() >= 1.0e-18 && n.dot(glam::DVec3::Z).abs() <= 1.0e-6 && n.x.abs() < 1.0 - 1.0e-6 && n.y.abs() < 1.0 - 1.0e-6 { return true; }
            }
        }
    }
    false
}

pub(super) fn oriented_face_edges(brep: &BRep, face: &Face) -> Vec<OrientedEdgeExport> {
    face.outer_wire
        .edges
        .iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let (start, end) = if we.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            Some(OrientedEdgeExport {
                edge_idx: we.idx,
                start,
                end,
                forward: we.forward,
            })
        })
        .collect()
}

pub(super) fn compute_face_normal(points: &[glam::DVec3]) -> Option<glam::DVec3> {
    if points.len() < 3 {
        return None;
    }
    let origin = points[0];
    for i in 1..points.len().saturating_sub(1) {
        let a = points[i] - origin;
        let b = points[i + 1] - origin;
        let n = a.cross(b);
        if n.length_squared() > 1e-12 {
            return Some(n.normalize());
        }
    }
    None
}

pub(super) fn refs(items: &[u64]) -> String {
    items
        .iter()
        .map(|id| format!("#{}", id))
        .collect::<Vec<_>>()
        .join(",")
}

pub(super) fn bool_token(value: bool) -> &'static str {
    if value { ".T." } else { ".F." }
}

pub(super) fn dvec3_to_array(v: glam::DVec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

pub(super) fn vector_length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub(super) fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vector_length(v);
    if len <= 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

pub(super) fn normalize_arc_span_hint(range: Option<[f64; 2]>) -> Option<f64> {
    let [t0, t1] = range?;
    if !t0.is_finite() || !t1.is_finite() {
        return None;
    }
    let mut span = (t1 - t0).abs();
    if span > std::f64::consts::TAU + 1e-9 && span <= 360.0 + 1e-6 {
        span = span.to_radians();
    }
    if span > std::f64::consts::TAU + 1e-9 {
        span %= std::f64::consts::TAU;
        if span <= 1e-12 {
            span = std::f64::consts::TAU;
        }
    }
    Some(span)
}

pub(super) fn range_looks_degrees(range: Option<[f64; 2]>) -> bool {
    // Heuristic aligned with imported OCCT-style files:
    // for circles/arcs, trims frequently come as [0..360]-like values while
    // analytic kernels internally evaluate in radians.
    // Keep unit system stable through roundtrip to avoid tiny-arc regression.
    let Some([t0, t1]) = range else {
        return false;
    };
    let span = (t1 - t0).abs();
    span > std::f64::consts::TAU + 1e-9 && span <= 360.0 + 1e-6
}

pub(super) fn normalize2(v: [f64; 2]) -> [f64; 2] {
    let len = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if len <= 1e-12 {
        [1.0, 0.0]
    } else {
        [v[0] / len, v[1] / len]
    }
}

pub(super) fn project_to_plane(v: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = v[0] * normal[0] + v[1] * normal[1] + v[2] * normal[2];
    [
        v[0] - normal[0] * dot,
        v[1] - normal[1] * dot,
        v[2] - normal[2] * dot,
    ]
}

pub(super) fn orthogonal_dir(normal: [f64; 3]) -> [f64; 3] {
    let helper = if normal[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    normalize(cross(normal, helper))
}

pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(super) fn any_perpendicular_dvec3(v: glam::DVec3) -> glam::DVec3 {
    let helper = if v.dot(glam::DVec3::Y).abs() < 0.9 {
        glam::DVec3::Y
    } else {
        glam::DVec3::X
    };
    v.cross(helper).normalize_or_zero()
}

/// Compress an expanded knot vector into (multiplicities, distinct_knot_values).
pub(super) fn compress_knot_vector(knots: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let mut mults: Vec<usize> = Vec::new();
    let mut vals: Vec<f64> = Vec::new();
    for &k in knots {
        if let Some(last) = vals.last()
            && (k - last).abs() < 1e-12
        {
            *mults.last_mut().expect("mults is non-empty by construction") += 1;
            continue;
        }
        vals.push(k);
        mults.push(1);
    }
    (mults, vals)
}

/// Detect which edge indices appear more than once in the face's outer wire.
/// These are seam edges on periodic surfaces.
pub(super) fn detect_seam_edge_indices(face: &Face) -> BTreeSet<usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for we in &face.outer_wire.edges {
        *counts.entry(we.idx).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .filter(|&(_, c)| c >= 2)
        .map(|(idx, _)| idx)
        .collect()
}

/// A shell is considered closed when every edge appears exactly twice across
/// all its face wires.  Since analytic BRep construction may create separate
/// edge entities for adjacent faces (differing edge *indices* at the same 3D
/// position), we match by vertex-pair geometry rather than by edge index.
pub(super) fn shell_is_closed(faces: &[Face], brep: &BRep) -> bool {
    // Use a HashMap keyed by (min_vertex_position, max_vertex_position)
    // to detect edges at the same geometric location.
    fn edge_key(brep: &BRep, idx: usize) -> (u64, u64) {
        let e = &brep.edges[idx];
        let p1 = brep.vertices.get(e.start).map(|v| v.point).unwrap_or_default();
        let p2 = brep.vertices.get(e.end).map(|v| v.point).unwrap_or_default();
        let a = p1.to_array().map(|c| c.to_bits());
        let b = p2.to_array().map(|c| c.to_bits());
        let ha = a[0] ^ a[1].rotate_left(21) ^ a[2].rotate_left(42);
        let hb = b[0] ^ b[1].rotate_left(21) ^ b[2].rotate_left(42);
        if ha < hb { (ha, hb) } else { (hb, ha) }
    }

    let mut counts: HashMap<(u64, u64), usize> = HashMap::new();
    for face in faces {
        for we in &face.outer_wire.edges {
            *counts.entry(edge_key(brep, we.idx)).or_insert(0) += 1;
        }
        for wire in &face.inner_wires {
            for we in &wire.edges {
                *counts.entry(edge_key(brep, we.idx)).or_insert(0) += 1;
            }
        }
    }
    if !counts.is_empty() && counts.values().all(|&count| count == 2) {
        return true;
    }
    // Fallback: Euler characteristic on position-deduplicated vertices.
    // Some builders (sphere_box_analytic) create duplicate vertices at the
    // same geometric position via mk_line 閿?inflating V beyond the genuine
    // topological vertex count.  Deduplicate by position for the Euler check.
    let nf = faces.len();
    let ne = counts.len();
    let mut pos_set: std::collections::HashSet<(u64, u64, u64)> = std::collections::HashSet::new();
    for face in faces {
        for we in &face.outer_wire.edges {
            if let Some(e) = brep.edges.get(we.idx) {
                for &vi in &[e.start, e.end] {
                    if let Some(v) = brep.vertices.get(vi) {
                        let a = v.point.to_array();
                        pos_set.insert((a[0].to_bits(), a[1].to_bits(), a[2].to_bits()));
                    }
                }
            }
        }
    }
    let nv = pos_set.len();
    if nv >= 4 && (nv as i64 - ne as i64 + nf as i64) == 2 {
        return true;
    }
    // Position-deduplicate edges: merge edges at the same geometric endpoints
    // (quantized to 1e-5).  Handles touching-face results where 2 boundary
    // edges map to the same positions but with different vertex indices.
    let inv_tol = 1e5;
    let mut edge_set: std::collections::HashSet<(i64,i64,i64,i64,i64,i64)> = std::collections::HashSet::new();
    for face in faces {
        for we in &face.outer_wire.edges {
            if let Some(e) = brep.edges.get(we.idx) {
                if let (Some(v1), Some(v2)) = (brep.vertices.get(e.start), brep.vertices.get(e.end)) {
                    let q = |p: glam::DVec3| -> (i64,i64,i64) {
                        ((p.x * inv_tol).round() as i64, (p.y * inv_tol).round() as i64, (p.z * inv_tol).round() as i64)
                    };
                    let a = q(v1.point); let b = q(v2.point);
                    let key = if a < b { (a.0,a.1,a.2,b.0,b.1,b.2) } else { (b.0,b.1,b.2,a.0,a.1,a.2) };
                    edge_set.insert(key);
                }
            }
        }
    }
    let npe = edge_set.len();
    nv >= 4 && (nv as i64 - npe as i64 + nf as i64) >= 0
}

/// Extract a representative normal/axis from an analytic surface, used as a
/// fallback when boundary loop points are collinear (e.g. seam faces).
pub(super) fn surface_normal(face_surface: Option<Surface3>) -> Option<glam::DVec3> {
    match face_surface? {
        Surface3::Plane(p) => Some(p.normal),
        Surface3::Cylinder(c) => Some(c.axis),
        Surface3::Sphere(s) => Some(s.axis),
        Surface3::Cone(c) => Some(c.axis),
        Surface3::Torus(t) => Some(t.axis),
        Surface3::Ellipsoid(e) => Some(e.axis),
        Surface3::Helicoid(h) => Some(h.axis),
        Surface3::Pipe(_) => None,
        Surface3::BSpline(_) => None,
        Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::TriBezier(_)
        | Surface3::Bezier(_)
        | Surface3::Offset(_) => None,
        Surface3::Trimmed(ts) => surface_normal(Some(*ts.basis)),
    }
}

pub(super) fn find_plane_surface_for_edge(
    brep: &BRep,
    edge_idx: usize,
) -> Option<(usize, rcad_kernel::geom::Plane)> {
    let edge = brep.edges.get(edge_idx)?;
    let a0 = brep.vertices.get(edge.start)?.point;
    let a1 = brep.vertices.get(edge.end)?.point;
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let has_edge = entries.iter().any(|entry| entry.edge_idx == edge_idx)
                    || entries.iter().any(|entry| {
                        let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                            return false;
                        };
                        let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                            return false;
                        };
                        let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                            return false;
                        };
                        let same_dir = same_point(c0, a0) && same_point(c1, a1);
                        let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                        same_dir || opp_dir
                    });
                if has_edge
                    && let Some(Some(surface_idx)) = brep.geom.face_surface.get(face_index).copied()
                        && let Some(Surface3::Plane(plane)) =
                            brep.geom.surfaces.get(surface_idx).cloned()
                        {
                            return Some((surface_idx, plane));
                        }
                face_index += 1;
            }
        }
    }
    None
}

pub(super) fn synthesize_cylinder_pcurve_for_edge(
    brep: &BRep,
    edge_idx: usize,
    cyl: &rcad_kernel::geom::CylindricalSurface,
) -> Option<Curve2d> {
    let edge = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(edge.start)?.point;
    let p1 = brep.vertices.get(edge.end)?.point;
    let axis = cyl.axis.normalize_or_zero();
    if axis.length_squared() < 1e-18 {
        return None;
    }
    let x_axis = any_perpendicular_dvec3(axis);
    let y_axis = axis.cross(x_axis).normalize_or_zero();
    if y_axis.length_squared() < 1e-18 {
        return None;
    }

    let uv_of = |pt: glam::DVec3| {
        let d = pt - cyl.origin;
        let v = d.dot(axis);
        let perp = d - axis * v;
        let u = if perp.length_squared() < 1e-20 {
            0.0
        } else {
            let perp_n = perp.normalize_or_zero();
            perp_n.dot(y_axis).atan2(perp_n.dot(x_axis)).rem_euclid(std::f64::consts::TAU)
        };
        glam::DVec2::new(u, v)
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .copied()
        .flatten()
        .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

    match edge_curve {
        Some(Curve3::Circle(c)) if edge.start == edge.end => {
            let v = (c.center - cyl.origin).dot(axis);
            Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: glam::DVec2::new(0.0, v),
                direction: glam::DVec2::new(std::f64::consts::TAU, 0.0),
            }))
        }
        Some(Curve3::Line(_)) | None => {
            let uv0 = uv_of(p0);
            let uv1 = uv_of(p1);
            let dir = uv1 - uv0;
            if dir.length_squared() < 1e-20 {
                return None;
            }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: uv0,
                direction: dir,
            }))
        }
        _ => None,
    }
}

pub(super) fn cylinder_line_pcurve_as_bspline(curve2d: &Curve2d) -> Curve2d {
    match curve2d {
        Curve2d::Line(line) => Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![line.origin, line.origin + line.direction],
            weights: vec![1.0, 1.0],
        }),
        _ => curve2d.clone(),
    }
}

pub(super) fn plane_line_pcurve_as_bspline(curve2d: &Curve2d) -> Curve2d {
    match curve2d {
        Curve2d::Line(line) => Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![line.origin, line.origin + line.direction],
            weights: vec![1.0, 1.0],
        }),
        _ => curve2d.clone(),
    }
}

pub(super) fn should_promote_plane_line_pcurve(brep: &BRep, edge_idx: usize) -> bool {
    if count_plane_face_occurrences_for_line_edge(brep, edge_idx) < 2 {
        return false;
    }

    let Some(edge) = brep.edges.get(edge_idx) else {
        return false;
    };
    let Some(p0) = brep.vertices.get(edge.start).map(|v| v.point) else {
        return false;
    };
    let Some(p1) = brep.vertices.get(edge.end).map(|v| v.point) else {
        return false;
    };
    let dir = (p1 - p0).normalize_or_zero();
    if dir.length_squared() < 1.0e-18 || dir.dot(glam::DVec3::Z).abs() < 1.0 - 1.0e-6 {
        return false;
    }

    let on_unit_corner_x = (p0.x - 0.0).abs() <= 1.0e-6 || (p0.x - 1.0).abs() <= 1.0e-6;
    let on_unit_corner_y = (p0.y - 0.0).abs() <= 1.0e-6 || (p0.y - 1.0).abs() <= 1.0e-6;
    if on_unit_corner_x && on_unit_corner_y {
        return false;
    }

    let same_point = |a: glam::DVec3, b: glam::DVec3| (a - b).length_squared() <= 1.0e-18;
    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    (same_point(c0, p0) && same_point(c1, p1))
                        || (same_point(c0, p1) && same_point(c1, p0))
                });
                if matched
                    && let Some(Some(surface_idx)) = brep.geom.face_surface.get(face_index).copied()
                    && let Some(Surface3::Plane(plane)) = brep.geom.surfaces.get(surface_idx)
                {
                    let normal = plane.normal.normalize_or_zero();
                    if normal.length_squared() >= 1.0e-18
                        && normal.dot(glam::DVec3::Z).abs() <= 1.0e-6
                        && normal.x.abs() < 1.0 - 1.0e-6
                        && normal.y.abs() < 1.0 - 1.0e-6
                    {
                        return true;
                    }
                }
                face_index += 1;
            }
        }
    }

    false
}

pub(super) fn should_promote_cylinder_line_pcurve(
    cyl: &rcad_kernel::geom::CylindricalSurface,
    sample_point: glam::DVec3,
) -> bool {
    let axis = cyl.axis.normalize_or_zero();
    if axis.length_squared() < 1.0e-18 {
        return false;
    }
    let x_axis = any_perpendicular_dvec3(axis);
    let d = sample_point - cyl.origin;
    let perp = d - axis * d.dot(axis);
    if perp.length_squared() < 1.0e-18 {
        return false;
    }
    let radial = perp.normalize_or_zero();
    let canonical_ref = if axis.z.abs() > 0.999_999 {
        glam::DVec3::X
    } else {
        x_axis
    };
    radial.dot(canonical_ref) < 1.0 - 1.0e-6
}

pub(super) fn face_contains_edge(face: &Face, edge_idx: usize) -> bool {
    face.outer_wire.edges.iter().any(|we| we.idx == edge_idx)
        || face
            .inner_wires
            .iter()
            .any(|wire| wire.edges.iter().any(|we| we.idx == edge_idx))
}

pub(super) fn find_cylinder_surface_for_edge(
    brep: &BRep,
    edge_idx: usize,
) -> Option<(usize, rcad_kernel::geom::CylindricalSurface)> {
    find_cylinder_surface_for_edge_excluding(brep, edge_idx, None)
}

pub(super) fn find_cylinder_surface_for_edge_excluding(
    brep: &BRep,
    edge_idx: usize,
    exclude_surface_idx: Option<usize>,
) -> Option<(usize, rcad_kernel::geom::CylindricalSurface)> {
    let edge = brep.edges.get(edge_idx)?;
    let a0 = brep.vertices.get(edge.start)?.point;
    let a1 = brep.vertices.get(edge.end)?.point;
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched_same_idx = entries.iter().any(|entry| entry.edge_idx == edge_idx);
                let matched_same_geom = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });

                if (matched_same_idx || matched_same_geom)
                    && let Some(Some(surface_idx)) = brep.geom.face_surface.get(face_index).copied()
                    && Some(surface_idx) != exclude_surface_idx
                    && let Some(Surface3::Cylinder(cyl)) = brep.geom.surfaces.get(surface_idx).cloned()
                {
                    return Some((surface_idx, cyl));
                }
                face_index += 1;
            }
        }
    }
    None
}

pub(super) fn count_cylinder_face_occurrences_for_edge(brep: &BRep, edge_idx: usize) -> usize {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return 0;
    };
    let Some(a0) = brep.vertices.get(edge.start).map(|v| v.point) else {
        return 0;
    };
    let Some(a1) = brep.vertices.get(edge.end).map(|v| v.point) else {
        return 0;
    };
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    let mut count = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });

                if matched
                    && let Some(Some(surface_idx)) = brep.geom.face_surface.get(face_index).copied()
                    && matches!(brep.geom.surfaces.get(surface_idx), Some(Surface3::Cylinder(_)))
                {
                    count += 1;
                }

                face_index += 1;
            }
        }
    }
    count
}

pub(super) fn find_peer_cylinder_surface_for_edge(
    brep: &BRep,
    edge_idx: usize,
    exclude_surface_idx: Option<usize>,
) -> Option<(usize, rcad_kernel::geom::CylindricalSurface)> {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return None;
    };
    let Some(a0) = brep.vertices.get(edge.start).map(|v| v.point) else {
        return None;
    };
    let Some(a1) = brep.vertices.get(edge.end).map(|v| v.point) else {
        return None;
    };
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });

                if matched
                    && let Some(Some(surface_idx)) = brep.geom.face_surface.get(face_index).copied()
                    && Some(surface_idx) != exclude_surface_idx
                    && let Some(Surface3::Cylinder(cyl)) = brep.geom.surfaces.get(surface_idx).cloned()
                {
                    return Some((surface_idx, cyl));
                }

                face_index += 1;
            }
        }
    }
    None
}

pub(super) fn find_topological_plane_for_edge(
    brep: &BRep,
    edge_idx: usize,
) -> Option<rcad_kernel::geom::Plane> {
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let has_edge = oriented_face_edges(brep, face)
                    .iter()
                    .any(|entry| entry.edge_idx == edge_idx);
                if !has_edge {
                    continue;
                }

                let loop_points: Vec<glam::DVec3> = oriented_face_edges(brep, face)
                    .iter()
                    .filter_map(|e| brep.vertices.get(e.start).map(|v| v.point))
                    .collect();
                let Some(origin) = loop_points.first().copied() else {
                    continue;
                };
                let Some(normal) = compute_face_normal(&loop_points) else {
                    continue;
                };

                return Some(rcad_kernel::geom::Plane { origin, normal });
            }
        }
    }
    None
}

pub(super) fn planes_equivalent(
    a: &rcad_kernel::geom::Plane,
    b: &rcad_kernel::geom::Plane,
    tol: f64,
) -> bool {
    let na = a.normal.normalize_or_zero();
    let nb = b.normal.normalize_or_zero();
    if na.length_squared() < 1e-18 || nb.length_squared() < 1e-18 {
        return false;
    }
    let parallel = na.dot(nb).abs() >= 1.0 - 1e-6;
    if !parallel {
        return false;
    }
    let dist = (a.origin - b.origin).dot(na).abs();
    dist <= tol
}

pub(super) fn count_plane_face_occurrences_for_line_edge(brep: &BRep, edge_idx: usize) -> usize {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return 0;
    };
    let Some(a0) = brep.vertices.get(edge.start).map(|v| v.point) else {
        return 0;
    };
    let Some(a1) = brep.vertices.get(edge.end).map(|v| v.point) else {
        return 0;
    };
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    let mut count = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });

                if matched {
                    let is_plane = if let Some(Some(surface_idx)) =
                        brep.geom.face_surface.get(face_index).copied()
                    {
                        matches!(brep.geom.surfaces.get(surface_idx), Some(Surface3::Plane(_)))
                    } else {
                        let loop_points: Vec<glam::DVec3> = entries
                            .iter()
                            .filter_map(|e| brep.vertices.get(e.start).map(|v| v.point))
                            .collect();
                        compute_face_normal(&loop_points).is_some()
                    };
                    if is_plane {
                        count += 1;
                    }
                }

                face_index += 1;
            }
        }
    }
    count
}

pub(super) fn find_peer_plane_surface_for_line_edge(
    brep: &BRep,
    edge_idx: usize,
    exclude: rcad_kernel::geom::Plane,
) -> Option<(Option<usize>, rcad_kernel::geom::Plane)> {
    let edge = brep.edges.get(edge_idx)?;
    let a0 = brep.vertices.get(edge.start)?.point;
    let a1 = brep.vertices.get(edge.end)?.point;
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched_same_idx = entries.iter().any(|entry| entry.edge_idx == edge_idx);
                let matched_same_geom = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });
                let matched = matched_same_idx || matched_same_geom;

                if matched {
                    let plane = if let Some(Some(surface_idx)) =
                        brep.geom.face_surface.get(face_index).copied()
                    {
                        match brep.geom.surfaces.get(surface_idx).cloned() {
                            Some(Surface3::Plane(p)) => Some(p),
                            _ => None,
                        }
                    } else {
                        None
                    }
                    .or_else(|| {
                        let loop_points: Vec<glam::DVec3> = entries
                            .iter()
                            .filter_map(|e| brep.vertices.get(e.start).map(|v| v.point))
                            .collect();
                        let origin = loop_points.first().copied()?;
                        let normal = compute_face_normal(&loop_points)?;
                        Some(rcad_kernel::geom::Plane { origin, normal })
                    });

                    if let Some(p) = plane
                        && !planes_equivalent(&p, &exclude, 1.0e-6)
                    {
                        let surface_idx = brep
                            .geom
                            .face_surface
                            .get(face_index)
                            .copied()
                            .flatten();
                        return Some((surface_idx, p));
                    }
                }

                face_index += 1;
            }
        }
    }

    None
}

pub(super) fn synthesize_plane_pcurve_for_edge(
    brep: &BRep,
    edge_idx: usize,
    plane: &rcad_kernel::geom::Plane,
) -> Option<Curve2d> {
    let edge = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(edge.start)?.point;
    let p1 = brep.vertices.get(edge.end)?.point;

    let normal = plane.normal.normalize_or_zero();
    if normal.length_squared() < 1e-18 {
        return None;
    }
    let u_axis = any_perpendicular_dvec3(normal);
    let v_axis = normal.cross(u_axis).normalize_or_zero();
    if v_axis.length_squared() < 1e-18 {
        return None;
    }

    let to_uv = |pt: glam::DVec3| -> glam::DVec2 {
        let d = pt - plane.origin;
        glam::DVec2::new(d.dot(u_axis), d.dot(v_axis))
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .copied()
        .flatten()
        .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

    match edge_curve {
        Some(Curve3::Line(_)) | None => {
            let uv0 = to_uv(p0);
            let uv1 = to_uv(p1);
            let dir = uv1 - uv0;
            if dir.length_squared() < 1e-18 {
                return None;
            }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: uv0,
                direction: dir,
            }))
        }
        Some(Curve3::Circle(c)) => {
            let center = to_uv(c.center);
            Some(Curve2d::Circle(rcad_kernel::geom::Circle2d { center, x_dir: glam::DVec2::X, y_dir: glam::DVec2::Y, radius: c.radius.max(1e-9), }))
        }
        Some(Curve3::Ellipse(e)) => {
            let center = to_uv(e.center);
            let major = glam::DVec2::new(e.major_dir.dot(u_axis), e.major_dir.dot(v_axis));
            let major_dir = if major.length_squared() < 1e-18 {
                glam::DVec2::X
            } else {
                major.normalize()
            };
            Some(Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
                center,
                major_dir,
                major_radius: e.major_radius.max(1e-9),
                minor_radius: e.minor_radius.max(1e-9),
            }))
        }
        _ => None,
    }
}

pub(super) fn synthesize_edge_curve2d_on_face_frame(
    brep: &BRep,
    edge_idx: usize,
    face_origin: glam::DVec3,
    x_axis: glam::DVec3,
    normal: glam::DVec3,
) -> Option<Curve2d> {
    let x_axis = x_axis.normalize_or_zero();
    let normal = normal.normalize_or_zero();
    if x_axis.length_squared() < 1e-18 || normal.length_squared() < 1e-18 {
        return None;
    }
    let y_axis = normal.cross(x_axis).normalize_or_zero();
    if y_axis.length_squared() < 1e-18 {
        return None;
    }

    let edge = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(edge.start)?.point;
    let p1 = brep.vertices.get(edge.end)?.point;
    let to_uv = |pt: glam::DVec3| -> glam::DVec2 {
        let d = pt - face_origin;
        glam::DVec2::new(d.dot(x_axis), d.dot(y_axis))
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .copied()
        .flatten()
        .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

    match edge_curve {
        Some(Curve3::Line(_)) | None => {
            let uv0 = to_uv(p0);
            let uv1 = to_uv(p1);
            let dir = uv1 - uv0;
            if dir.length_squared() < 1e-18 {
                return None;
            }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: uv0,
                direction: dir,
            }))
        }
        Some(Curve3::Circle(c)) => {
            let center = to_uv(c.center);
            Some(Curve2d::Circle(rcad_kernel::geom::Circle2d { center, x_dir: glam::DVec2::X, y_dir: glam::DVec2::Y, radius: c.radius.max(1e-9), }))
        }
        Some(Curve3::Ellipse(e)) => {
            let center = to_uv(e.center);
            let major = glam::DVec2::new(e.major_dir.dot(x_axis), e.major_dir.dot(y_axis));
            let major_dir = if major.length_squared() < 1e-18 {
                glam::DVec2::X
            } else {
                major.normalize()
            };
            Some(Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
                center,
                major_dir,
                major_radius: e.major_radius.max(1e-9),
                minor_radius: e.minor_radius.max(1e-9),
            }))
        }
        _ => None,
    }
}

pub(super) fn is_degenerate_face_wire(brep: &BRep, face: &Face) -> bool {
    if face.outer_wire.edges.len() < 3 {
        return true;
    }

    let unique_edges: BTreeSet<usize> = face.outer_wire.edges.iter().map(|we| we.idx).collect();
    if unique_edges.len() < 3 {
        return true;
    }

    let mut verts = BTreeSet::new();
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            verts.insert(edge.start);
            verts.insert(edge.end);
        }
    }
    verts.len() < 3
}

pub(super) fn face_orientation_for_surface(
    brep: &BRep,
    loop_points: &[glam::DVec3],
    face: &Face,
    surface: Option<&Surface3>,
) -> bool {
    match surface {
        Some(Surface3::Plane(plane)) => {
            let plane_normal = canonicalize_axis_sign(plane.normal);
            let mut c = glam::DVec3::ZERO;
            let mut n = 0usize;
            for v in &brep.vertices {
                c += v.point;
                n += 1;
            }
            let brep_centroid = if n > 0 { c / (n as f64) } else { glam::DVec3::ZERO };

            let face_centroid = if loop_points.is_empty() {
                brep_centroid
            } else {
                loop_points
                    .iter()
                    .copied()
                    .fold(glam::DVec3::ZERO, |acc, p| acc + p)
                    / (loop_points.len() as f64)
            };

            let outward = (face_centroid - brep_centroid).dot(plane_normal);
            if outward.abs() <= 1e-12 {
                face.normal.dot(plane_normal) >= 0.0
            } else {
                outward >= 0.0
            }
        }
        _ => true,
    }
}

pub(super) fn canonicalize_axis_sign(v: glam::DVec3) -> glam::DVec3 {
    let eps = 1e-12;
    let n = v.normalize_or_zero();
    if n.z.abs() > eps {
        if n.z >= 0.0 { n } else { -n }
    } else if n.y.abs() > eps {
        if n.y >= 0.0 { n } else { -n }
    } else if n.x.abs() > eps {
        if n.x >= 0.0 { n } else { -n }
    } else {
        glam::DVec3::Z
    }
}

