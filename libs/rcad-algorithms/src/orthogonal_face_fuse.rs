//! Fuse coplanar **axis-aligned** rectangular patches into one [`Face`] using a 2D
//! axis-aligned rectangle union on a grid, producing one outer boundary and optional
//! inner wires (holes). Complements [`unify_same_domain_faces`](crate::unify_same_domain_faces),
//! which only merges along shared **edges** and leaves corner-only adjacency split.

use glam::DVec3;
use std::collections::{HashMap, HashSet};
use rcad_kernel::geom::{Plane, Surface3};
use rcad_kernel::topology::{Edge, Face, Vertex, Wire, WireEdge};
use rcad_kernel::BRep;

use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::TOLERANCE_ABS;

type Pt = (i64, i64);

fn qpt(x: f64, y: f64, scale: f64) -> Pt {
    ((x * scale).round() as i64, (y * scale).round() as i64)
}

/// Merge groups of coplanar orthogonal faces in each shell into single faces (with holes).
pub fn fuse_orthogonal_coplanar_faces(brep: &BRep, tol: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total = 0usize;
    let t = tol.max(TOLERANCE_ABS);

    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            fuse_orthogonal_in_shell(&mut out, si, shi, t, &mut total);
        }
    }

    (out, total)
}

fn fuse_orthogonal_in_shell(brep: &mut BRep, si: usize, shi: usize, tol: f64, total: &mut usize) {
    let n = brep.solids[si].shells[shi].faces.len();
    if n < 2 {
        return;
    }

    let mut groups: std::collections::HashMap<(i64, i64, i64, i64), Vec<usize>> =
        std::collections::HashMap::new();
    for fi in 0..n {
        let face = &brep.solids[si].shells[shi].faces[fi];
        if !face.inner_wires.is_empty() {
            continue;
        }
        let Some(p0) = face_first_point(brep, face) else {
            continue;
        };
        let nrm = face.normal.normalize_or_zero();
        if nrm.length_squared() < 1e-24 {
            continue;
        }
        let d = nrm.dot(p0);
        let key = plane_key(nrm, d, tol);
        groups.entry(key).or_default().push(fi);
    }

    for (_, mut fis) in groups {
        if fis.len() < 2 {
            continue;
        }
        fis.sort_unstable();
        if try_fuse_orthogonal_group(brep, si, shi, &fis, tol) {
            *total += 1;
        }
    }
}

fn plane_key(n: DVec3, d: f64, tol: f64) -> (i64, i64, i64, i64) {
    let s = 1.0 / tol.max(1e-9);
    (
        (n.x * s).round() as i64,
        (n.y * s).round() as i64,
        (n.z * s).round() as i64,
        (d * s).round() as i64,
    )
}

fn face_first_point(brep: &BRep, face: &Face) -> Option<DVec3> {
    let we = face.outer_wire.edges.first()?;
    let e = brep.edges.get(we.idx)?;
    let vi = if we.forward { e.start } else { e.end };
    Some(brep.vertices.get(vi)?.point)
}

fn try_fuse_orthogonal_group(
    brep: &mut BRep,
    si: usize,
    shi: usize,
    fis: &[usize],
    tol: f64,
) -> bool {
    let n_faces_shell = brep.solids[si].shells[shi].faces.len();
    if fis.is_empty() || fis[0] >= n_faces_shell {
        return false;
    }
    let plane = {
        let f0 = &brep.solids[si].shells[shi].faces[fis[0]];
        let n = f0.normal.normalize_or_zero();
        let Some(p0) = face_first_point(brep, f0) else {
            return false;
        };
        Plane {
            origin: p0,
            normal: n,
        }
    };
    let (u_axis, v_axis) = plane_local_basis(&plane);
    let to_uv = |p: DVec3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(u_axis), d.dot(v_axis))
    };

    let mut rects: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(fis.len());
    for &fi in fis {
        if fi >= n_faces_shell {
            return false;
        }
        let face = &brep.solids[si].shells[shi].faces[fi];
        if !face.inner_wires.is_empty() {
            return false;
        }
        let poly = face_outer_points(brep, face);
        if poly.len() < 3 || poly.len() > 64 {
            return false;
        }
        if !polygon_is_orthogonal(&poly, &to_uv, tol) {
            return false;
        }
        let uv: Vec<(f64, f64)> = poly.iter().map(|&p| to_uv(p)).collect();
        let (umin, umax, vmin, vmax) = bbox2d(&uv);
        rects.push((umin, umax, vmin, vmax));
    }

    let Some(geo_rings) = union_rects_to_rings_grid(&rects, tol) else {
        return false;
    };
    if geo_rings.is_empty() {
        return false;
    }

    let normal = brep.solids[si].shells[shi].faces[fis[0]].normal;

    let ring_vertices = add_vertices_for_rings(brep, &geo_rings, &plane, u_axis, v_axis, tol);
    if ring_vertices.is_empty() || ring_vertices[0].len() < 3 {
        return false;
    }

    let mut edge_pairs: Vec<(usize, usize)> = Vec::new();
    for ring_v in &ring_vertices {
        let n = ring_v.len();
        if n < 3 {
            return false;
        }
        for i in 0..n {
            edge_pairs.push((ring_v[i], ring_v[(i + 1) % n]));
        }
    }

    let base_ei = brep.edges.len();
    push_new_edges(brep, edge_pairs);

    let mut ei = base_ei;
    let n0 = ring_vertices[0].len();
    let outer_wire = Wire {
        edges: (0..n0).map(|i| {
            let w = WireEdge::fwd(ei + i);
            w
        })
        .collect(),
    };
    ei += n0;

    let mut inner_wires: Vec<Wire> = Vec::new();
    for ring_v in ring_vertices.iter().skip(1) {
        let n = ring_v.len();
        if n < 3 {
            return false;
        }
        inner_wires.push(Wire {
            edges: (0..n)
                .map(|i| {
                    let w = WireEdge::fwd(ei + i);
                    w
                })
                .collect(),
        });
        ei += n;
    }

    let merged_face = Face {
        outer_wire,
        inner_wires,
        normal,
        triangles: vec![],
        mesh_dirty: true,
    };

    let surf_idx = {
        let p0 = face_first_point(brep, &merged_face).unwrap_or(plane.origin);
        let idx = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: p0,
            normal,
        }));
        idx
    };

    let flat_indices: Vec<usize> = fis
        .iter()
        .map(|&fi| flat_face_index(brep, si, shi, fi))
        .collect();

    replace_shell_faces_and_geom(brep, si, shi, fis, merged_face, surf_idx, &flat_indices);
    true
}

/// Union of axis-aligned rectangles → outer ring first, then hole rings (UV coords).
fn union_rects_to_rings_grid(rects: &[(f64, f64, f64, f64)], tol: f64) -> Option<Vec<Vec<(f64, f64)>>> {
    let t = tol.max(TOLERANCE_ABS);
    let (occ_ext, xs, ys) = build_padded_occ_grid(rects, t)?;
    let nx = occ_ext.len();
    let ny = occ_ext.first()?.len();
    if nx < 3 || ny < 3 {
        return None;
    }

    let mut outside = vec![vec![false; ny]; nx];
    let mut stack = vec![(0usize, 0usize)];
    outside[0][0] = true;
    while let Some((i, j)) = stack.pop() {
        for (di, dj) in [(0, 1isize), (0, -1), (1, 0), (-1, 0)] {
            let ni = i as isize + di;
            let nj = j as isize + dj;
            if ni < 0 || nj < 0 || ni >= nx as isize || nj >= ny as isize {
                continue;
            }
            let ni = ni as usize;
            let nj = nj as usize;
            if occ_ext[ni][nj] || outside[ni][nj] {
                continue;
            }
            outside[ni][nj] = true;
            stack.push((ni, nj));
        }
    }

    let mut outer_segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    let mut hole_segs: Vec<((f64, f64), (f64, f64))> = Vec::new();

    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            let l = occ_ext[i][j];
            let r = occ_ext[i + 1][j];
            if l == r {
                continue;
            }
            let x = xs[i + 1];
            let y0 = ys[j];
            let y1 = ys[j + 1];
            let seg = ((x, y0), (x, y1));
            if l && !r {
                if outside[i + 1][j] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            } else if !l && r {
                if outside[i][j] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            }
        }
    }

    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            let b = occ_ext[i][j];
            let t = occ_ext[i][j + 1];
            if b == t {
                continue;
            }
            let y = ys[j + 1];
            let x0 = xs[i];
            let x1 = xs[i + 1];
            let seg = ((x0, y), (x1, y));
            if b && !t {
                if outside[i][j + 1] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            } else if !b && t {
                if outside[i][j] {
                    outer_segs.push(seg);
                } else {
                    hole_segs.push(seg);
                }
            }
        }
    }

    let scale = 1.0 / t.max(1e-12);
    let mut rings = Vec::new();
    if let Some(r0) = segments_to_ring(&outer_segs, scale) {
        rings.push(r0);
    } else {
        return None;
    }
    let mut hs = segments_to_rings(&hole_segs, scale)?;
    rings.append(&mut hs);
    Some(rings)
}

fn build_padded_occ_grid(
    rects: &[(f64, f64, f64, f64)],
    tol: f64,
) -> Option<(Vec<Vec<bool>>, Vec<f64>, Vec<f64>)> {
    let mut xs: Vec<f64> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    for &(u0, u1, v0, v1) in rects {
        xs.push(u0);
        xs.push(u1);
        ys.push(v0);
        ys.push(v1);
    }
    xs.sort_by(|a, b| a.total_cmp(b));
    xs.dedup_by(|a, b| (*a - *b).abs() <= tol * 10.0);
    ys.sort_by(|a, b| a.total_cmp(b));
    ys.dedup_by(|a, b| (*a - *b).abs() <= tol * 10.0);
    if xs.len() < 2 || ys.len() < 2 {
        return None;
    }
    let dx = ((xs[xs.len() - 1] - xs[0]).abs()).max(1.0) * 1e-4;
    let dy = ((ys[ys.len() - 1] - ys[0]).abs()).max(1.0) * 1e-4;
    let mut xs_e = vec![xs[0] - dx];
    xs_e.extend(xs.iter().cloned());
    xs_e.push(xs[xs.len() - 1] + dx);
    let mut ys_e = vec![ys[0] - dy];
    ys_e.extend(ys.iter().cloned());
    ys_e.push(ys[ys.len() - 1] + dy);

    let nx = xs_e.len() - 1;
    let ny = ys_e.len() - 1;
    let mut occ = vec![vec![false; ny]; nx];
    for i in 0..nx {
        for j in 0..ny {
            let cu = (xs_e[i] + xs_e[i + 1]) * 0.5;
            let cv = (ys_e[j] + ys_e[j + 1]) * 0.5;
            occ[i][j] = rects.iter().any(|&(u0, u1, v0, v1)| {
                cu + tol >= u0 && cu - tol <= u1 && cv + tol >= v0 && cv - tol <= v1
            });
        }
    }

    let nx2 = nx + 2;
    let ny2 = ny + 2;
    let mut occ_ext = vec![vec![false; ny2]; nx2];
    for i in 0..nx {
        for j in 0..ny {
            occ_ext[i + 1][j + 1] = occ[i][j];
        }
    }

    let mut xs2 = vec![xs_e[0] - dx];
    xs2.extend(xs_e.iter().cloned());
    xs2.push(xs_e[xs_e.len() - 1] + dx);
    let mut ys2 = vec![ys_e[0] - dy];
    ys2.extend(ys_e.iter().cloned());
    ys2.push(ys_e[ys_e.len() - 1] + dy);

    Some((occ_ext, xs2, ys2))
}

fn segments_to_ring(segs: &[((f64, f64), (f64, f64))], scale: f64) -> Option<Vec<(f64, f64)>> {
    let mut rings = segments_to_rings(segs, scale)?;
    if rings.is_empty() {
        return None;
    }
    rings.sort_by(|a, b| {
        ring_area_uv(b)
            .partial_cmp(&ring_area_uv(a))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    rings.into_iter().next()
}

fn segments_to_rings(segs: &[((f64, f64), (f64, f64))], scale: f64) -> Option<Vec<Vec<(f64, f64)>>> {
    if segs.is_empty() {
        return Some(vec![]);
    }
    let mut adj: HashMap<Pt, Vec<Pt>> = HashMap::new();
    let mut undir: HashSet<(Pt, Pt)> = HashSet::new();
    for &((x0, y0), (x1, y1)) in segs {
        let a = qpt(x0, y0, scale);
        let b = qpt(x1, y1, scale);
        let k0 = if a <= b { (a, b) } else { (b, a) };
        if !undir.insert(k0) {
            continue;
        }
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }

    let mut out: Vec<Vec<(f64, f64)>> = Vec::new();
    loop {
        let Some(&v0) = adj.keys().next() else { break };
        if adj.get(&v0).map_or(true, |n| n.is_empty()) {
            adj.remove(&v0);
            continue;
        }
        let Some(ring) = trace_one_ring(&mut adj, v0, scale) else {
            return None;
        };
        if ring.len() >= 3 {
            out.push(ring);
        }
    }
    Some(out)
}

fn remove_adj(adj: &mut HashMap<Pt, Vec<Pt>>, a: Pt, b: Pt) {
    if let Some(v) = adj.get_mut(&a) {
        if let Some(p) = v.iter().position(|&x| x == b) {
            v.remove(p);
        }
        if v.is_empty() {
            adj.remove(&a);
        }
    }
    if let Some(v) = adj.get_mut(&b) {
        if let Some(p) = v.iter().position(|&x| x == a) {
            v.remove(p);
        }
        if v.is_empty() {
            adj.remove(&b);
        }
    }
}

fn trace_one_ring(adj: &mut HashMap<Pt, Vec<Pt>>, v0: Pt, scale: f64) -> Option<Vec<(f64, f64)>> {
    let v1 = adj.get_mut(&v0)?.pop()?;
    remove_adj(adj, v0, v1);
    let to_f = |p: Pt| (p.0 as f64 / scale, p.1 as f64 / scale);
    let mut pts = vec![to_f(v0), to_f(v1)];
    let mut prev = v0;
    let mut cur = v1;
    loop {
        let nbs = adj.get_mut(&cur)?;
        if nbs.is_empty() {
            return None;
        }
        let pos = nbs.iter().position(|&x| x != prev)?;
        let nxt = nbs.remove(pos);
        if nbs.is_empty() {
            adj.remove(&cur);
        }
        remove_adj(adj, cur, nxt);
        if nxt == v0 {
            return Some(pts);
        }
        pts.push(to_f(nxt));
        prev = cur;
        cur = nxt;
        if pts.len() > 100_000 {
            return None;
        }
    }
}

fn ring_area_uv(ring: &[(f64, f64)]) -> f64 {
    if ring.len() < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    for i in 0..ring.len() {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % ring.len()];
        a += x0 * y1 - x1 * y0;
    }
    0.5 * a.abs()
}

fn flat_face_index(brep: &BRep, si: usize, shi: usize, fi: usize) -> usize {
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

fn replace_shell_faces_and_geom(
    brep: &mut BRep,
    si: usize,
    shi: usize,
    remove_fis: &[usize],
    merged_face: Face,
    surf_idx: usize,
    flat_indices: &[usize],
) {
    let remove: std::collections::HashSet<usize> = remove_fis.iter().copied().collect();
    let min_keep = *remove_fis.iter().min().unwrap();
    let shell = &mut brep.solids[si].shells[shi];
    let mut new_faces: Vec<Face> = Vec::with_capacity(shell.faces.len() - remove_fis.len() + 1);
    for fi in 0..shell.faces.len() {
        if remove.contains(&fi) {
            if fi == min_keep {
                new_faces.push(merged_face.clone());
            }
        } else {
            new_faces.push(shell.faces[fi].clone());
        }
    }
    shell.faces = new_faces;

    let mut flats: Vec<usize> = flat_indices.to_vec();
    flats.sort_unstable();
    let insert_at = flats[0];
    for &f in flats.iter().rev() {
        if f != insert_at {
            crate::remove_flat_face_geom_slots(&mut brep.geom, f);
        }
    }
    crate::remove_flat_face_geom_slots(&mut brep.geom, insert_at);
    brep.geom.face_surface.insert(insert_at, Some(surf_idx));
    if brep.geom.face_surface_range.len() < brep.geom.face_surface.len() {
        brep.geom
            .face_surface_range
            .resize(brep.geom.face_surface.len(), None);
    } else {
        brep.geom.face_surface_range.insert(insert_at, None);
    }
    if brep.geom.face_tolerance.len() < brep.geom.face_surface.len() {
        brep.geom
            .face_tolerance
            .resize(brep.geom.face_surface.len(), 0.0);
    } else {
        brep.geom.face_tolerance.insert(insert_at, 0.0);
    }
}

fn push_new_edges(brep: &mut BRep, edges: Vec<(usize, usize)>) {
    for (a, b) in edges {
        brep.edges.push(Edge { start: a, end: b });
        let ei = brep.edges.len() - 1;
        if brep.geom.edge_curve.len() <= ei {
            brep.geom.edge_curve.resize(ei + 1, None);
        }
        if brep.geom.edge_curve_range.len() <= ei {
            brep.geom.edge_curve_range.resize(ei + 1, None);
        }
        if brep.geom.edge_degenerated.len() <= ei {
            brep.geom.edge_degenerated.resize(ei + 1, false);
        }
        let p0 = brep.vertices[a].point;
        let p1 = brep.vertices[b].point;
        let delta = p1 - p0;
        let len = delta.length();
        let dir = if len > 1e-12 { delta / len } else { DVec3::X };
        let curve_idx = brep.geom.curves.len();
        brep.geom.curves.push(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
            origin: p0,
            direction: dir,
        }));
        brep.geom.edge_curve[ei] = Some(curve_idx);
        brep.geom.edge_curve_range[ei] = Some([0.0, (p1 - p0).dot(dir)]);
        brep.geom.edge_degenerated[ei] = len <= 1e-12;
    }
}

fn add_vertices_for_rings(
    brep: &mut BRep,
    rings: &[Vec<(f64, f64)>],
    plane: &Plane,
    u_axis: DVec3,
    v_axis: DVec3,
    tol: f64,
) -> Vec<Vec<usize>> {
    let match_tol = tol.max(TOLERANCE_ABS) * 50.0;
    let mut out: Vec<Vec<usize>> = Vec::new();
    for ring in rings {
        let mut vi: Vec<usize> = Vec::new();
        for &(u, v) in ring {
            let p = plane.origin + u_axis * u + v_axis * v;
            let mut found = None;
            for (i, vtx) in brep.vertices.iter().enumerate() {
                if (vtx.point - p).length() <= match_tol {
                    found = Some(i);
                    break;
                }
            }
            let idx = found.unwrap_or_else(|| {
                let i = brep.vertices.len();
                brep.vertices.push(Vertex { point: p });
                i
            });
            vi.push(idx);
        }
        out.push(vi);
    }
    out
}

// Fix try_fuse: use add_vertices_for_rings first, then build wires from vertex indices only

fn face_outer_points(brep: &BRep, face: &Face) -> Vec<DVec3> {
    let mut pts = Vec::new();
    for we in &face.outer_wire.edges {
        let Some((u, _)) = crate::oriented_edge_vertices(brep, *we) else {
            continue;
        };
        if let Some(v) = brep.vertices.get(u) {
            pts.push(v.point);
        }
    }
    pts
}

fn polygon_is_orthogonal(poly: &[DVec3], to_uv: &dyn Fn(DVec3) -> (f64, f64), tol: f64) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let (u0, v0) = to_uv(poly[i]);
        let (u1, v1) = to_uv(poly[(i + 1) % n]);
        let du = (u1 - u0).abs();
        let dv = (v1 - v0).abs();
        if du > tol && dv > tol {
            return false;
        }
        if du <= tol && dv <= tol {
            return false;
        }
    }
    true
}

fn bbox2d(uv: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for &(u, v) in uv {
        umin = umin.min(u);
        umax = umax.max(u);
        vmin = vmin.min(v);
        vmax = vmax.max(v);
    }
    (umin, umax, vmin, vmax)
}

