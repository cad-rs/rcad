//! Cylindrical rcad_kernel::BRep projection (OCCT `rcad_kernel::BRepProj_Projection` with `gp_Dir`).
//!
//! Strategy (matches OCCT): translate the wire by `-mdis·D`, build a prism by extruding
//! the polyline along `2·mdis·D`, then intersect that prism with the target shape’s
//! triangle soup to obtain projected polylines as wires.


use std::collections::HashSet;

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Line3};
use rcad_kernel::topology::{Edge, Vertex};
use rcad_kernel::topods;

use crate::brep_tools;
use crate::section::{brep_triangle_soup, intersect_triangle_soups_eps};
use crate::tolerance::*;

/// Options for [`brep_proj_cylindrical`].
///
/// Default [`BrepProjOptions::tolerance`] is [`TOLERANCE_MESH_LEGACY`] for backward compatibility.
/// [`brep_proj_cylindrical`] still **raises** soup / wire eps by `max` with
/// [`max_face_tolerance_or_abs_pair`](crate::tolerance::max_face_tolerance_or_abs_pair) and
/// [`tessellation_merge_linear_from_two_breps`](crate::tolerance::tessellation_merge_linear_from_two_breps)
/// (phase C). To force a stricter user floor, set `tolerance` below those only when you accept possible
/// dropped hits on coarse models.
#[derive(Debug, Clone)]
pub struct BrepProjOptions {
    /// Geometric tolerance (segment chaining, degenerate skips, triangle–triangle intersection).
    pub tolerance: f64,
    /// Uniform samples per edge when discretizing edge curves (≥ 2).
    pub samples_per_edge: usize,
}

impl Default for BrepProjOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_MESH_LEGACY,
            samples_per_edge: 24,
        }
    }
}

/// OCCT-style `DistanceIn`-like scale: diagonal(A)+diagonal(B)+bbox separation.
fn proj_distance_scale(wire_brep: &rcad_kernel::BRep, target: &rcad_kernel::BRep) -> f64 {
    let Some(ba) = brep_tools::bounding_box(wire_brep) else {
        return 1.0;
    };
    let Some(bb) = brep_tools::bounding_box(target) else {
        return 1.0;
    };
    let da = ba[1] - ba[0];
    let db = bb[1] - bb[0];
    let la = da.length();
    let lb = db.length();
    let sep = aabb_distance(ba, bb);
    (la + lb + sep.max(0.0)).max(TOLERANCE_MESH_LEGACY)
}

fn aabb_distance(a: [DVec3; 2], b: [DVec3; 2]) -> f64 {
    let mut d2 = 0.0_f64;
    for i in 0..3 {
        let amin = a[0][i].min(a[1][i]);
        let amax = a[0][i].max(a[1][i]);
        let bmin = b[0][i].min(b[1][i]);
        let bmax = b[0][i].max(b[1][i]);
        let gap = if amax < bmin {
            bmin - amax
        } else if bmax < amin {
            amin - bmax
        } else {
            0.0
        };
        d2 += gap * gap;
    }
    d2.sqrt()
}

fn first_outer_wire_edge_chain(brep: &rcad_kernel::BRep) -> Option<Vec<usize>> {
    let face = brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .next()?;
    if face.outer_wire.edges.is_empty() {
        return None;
    }
    Some(
        face.outer_wire
            .edges
            .iter()
            .map(|e| e.idx)
            .collect(),
    )
}

fn edge_chains_greedy(brep: &rcad_kernel::BRep) -> Vec<Vec<usize>> {
    let n = brep.edges.len();
    if n == 0 {
        return Vec::new();
    }
    let adj = build_vertex_edge_adj(brep);
    let mut remaining: HashSet<usize> = (0..n).collect();
    let mut chains: Vec<Vec<usize>> = Vec::new();

    while let Some(&start_e) = remaining.iter().next() {
        let e0 = &brep.edges[start_e];
        let mut v = e0.start;
        let mut next_e = start_e;
        let mut chain: Vec<usize> = Vec::new();

        loop {
            if !remaining.remove(&next_e) {
                break;
            }
            chain.push(next_e);
            let e = &brep.edges[next_e];
            v = if e.start == v { e.end } else { e.start };
            let mut cont: Option<usize> = None;
            for &cand in &adj[v] {
                if remaining.contains(&cand) {
                    cont = Some(cand);
                    break;
                }
            }
            match cont {
                Some(ei) => next_e = ei,
                None => break,
            }
        }
        if !chain.is_empty() {
            chains.push(chain);
        }
    }
    chains
}

fn build_vertex_edge_adj(brep: &rcad_kernel::BRep) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); brep.vertices.len()];
    for (ei, e) in brep.edges.iter().enumerate() {
        if e.start < adj.len() {
            adj[e.start].push(ei);
        }
        if e.end < adj.len() {
            adj[e.end].push(ei);
        }
    }
    adj
}

fn wire_edge_chains(brep: &rcad_kernel::BRep) -> Vec<Vec<usize>> {
    if let Some(chain) = first_outer_wire_edge_chain(brep) {
        if !chain.is_empty() {
            return vec![chain];
        }
    }
    edge_chains_greedy(brep)
}

fn push_pt(pts: &mut Vec<DVec3>, p: DVec3, tol: f64) {
    if pts.last().map_or(true, |q| (*q - p).length() > tol) {
        pts.push(p);
    }
}

fn sample_edge(
    brep: &rcad_kernel::BRep,
    ei: usize,
    from_v: usize,
    to_v: usize,
    n_seg: usize,
    tol: f64,
    pts: &mut Vec<DVec3>,
) {
    let edge = match brep.edges.get(ei) {
        Some(e) => e,
        None => return,
    };
    let n_seg = n_seg.max(2);
    if let Some(ci) = brep.geom.edge_curve.get(ei).copied().flatten() {
        if let Some(curve) = brep.geom.curves.get(ci) {
            let range = brep
                .geom
                .edge_curve_range
                .get(ei)
                .copied()
                .flatten()
                .unwrap_or_else(|| curve.default_domain());
            let forward = from_v == edge.start;
            let (t0, t1) = if forward {
                (range[0], range[1])
            } else {
                (range[1], range[0])
            };
            for i in 0..n_seg {
                let t = t0 + (t1 - t0) * i as f64 / (n_seg - 1) as f64;
                let p = curve.point_at(t);
                push_pt(pts, p, tol);
            }
            return;
        }
    }
    let p0 = brep.vertices[from_v].point;
    let p1 = brep.vertices[to_v].point;
    push_pt(pts, p0, tol);
    push_pt(pts, p1, tol);
}

fn dense_polyline_from_chain(brep: &rcad_kernel::BRep, chain: &[usize], options: &BrepProjOptions) -> Vec<DVec3> {
    if chain.is_empty() {
        return Vec::new();
    }
    let mut pts: Vec<DVec3> = Vec::new();
    let mut v = brep.edges[chain[0]].start;
    for &ei in chain {
        let e = &brep.edges[ei];
        let (from_v, to_v) = if e.start == v {
            (e.start, e.end)
        } else if e.end == v {
            (e.end, e.start)
        } else {
            (e.start, e.end)
        };
        sample_edge(
            brep,
            ei,
            from_v,
            to_v,
            options.samples_per_edge,
            options.tolerance,
            &mut pts,
        );
        v = to_v;
    }
    pts
}

fn extrusion_prism_tris(polyline: &[DVec3], extrusion: DVec3) -> Vec<[DVec3; 3]> {
    let mut tris = Vec::new();
    for i in 0..polyline.len().saturating_sub(1) {
        let a = polyline[i];
        let b = polyline[i + 1];
        let c = b + extrusion;
        let d = a + extrusion;
        tris.push([a, b, c]);
        tris.push([a, c, d]);
    }
    tris
}

/// Build a rcad_kernel::BRep containing only vertices and line edges (no solids), one edge per segment.
fn wire_brep_from_polyline(poly: &[DVec3], tol: f64) -> rcad_kernel::BRep {
    let mut brep = rcad_kernel::BRep::new();
    if poly.len() < 2 {
        return brep;
    }
    for p in poly {
        brep.vertices.push(Vertex { point: *p });
    }
    for i in 0..poly.len().saturating_sub(1) {
        let vi_a = i;
        let vi_b = i + 1;
        let a = brep.vertices[vi_a].point;
        let b = brep.vertices[vi_b].point;
        let d = b - a;
        let len = d.length();
        let ei = brep.edges.len();
        brep.edges.push(Edge {
            start: vi_a,
            end: vi_b,
        });
        let dir = if len > tol {
            d / len
        } else {
            DVec3::X
        };
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 {
            origin: a,
            direction: dir,
        }));
        while brep.geom.edge_curve.len() <= ei {
            brep.geom.edge_curve.push(None);
        }
        while brep.geom.edge_curve_range.len() <= ei {
            brep.geom.edge_curve_range.push(None);
        }
        while brep.geom.edge_degenerated.len() <= ei {
            brep.geom.edge_degenerated.push(false);
        }
        brep.geom.edge_curve[ei] = Some(ci);
        brep.geom.edge_curve_range[ei] = Some([0.0, len.max(tol)]);
    }
    brep
}

/// Cylindrical projection of wire-like `shape` onto `target` along `direction`.
///
/// Returns one rcad_kernel::BRep per connected intersection chain (OCCT `rcad_kernel::BRepProj_Projection::More` /
/// `Current`), each holding line edges approximating the image on `target`.
pub fn brep_proj_cylindrical(
    shape: &rcad_kernel::BRep,
    target: &rcad_kernel::BRep,
    direction: DVec3,
    options: &BrepProjOptions,
) -> Vec<rcad_kernel::BRep> {
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
        return Vec::new();
    }
    let target_tris = brep_triangle_soup(target);
    if target_tris.is_empty() {
        return Vec::new();
    }
    let mdis = proj_distance_scale(shape, target);
    let v_sup = dir * (2.0 * mdis);
    let v_inf = dir * (-mdis);

    let chains = wire_edge_chains(shape);
    if chains.is_empty() {
        return Vec::new();
    }

    let mut out: Vec<rcad_kernel::BRep> = Vec::new();
    for chain in chains {
        let pts = dense_polyline_from_chain(shape, &chain, options);
        if pts.len() < 2 {
            continue;
        }
        let shifted: Vec<DVec3> = pts.iter().map(|p| *p + v_inf).collect();
        let prism_tris = extrusion_prism_tris(&shifted, v_sup);
        if prism_tris.is_empty() {
            continue;
        }
        let model_te = max_face_tolerance_or_abs_pair(shape, target);
        let mesh_te = tessellation_merge_linear_from_two_breps(shape, target);
        let te = options
            .tolerance
            .max(model_te)
            .max(mesh_te)
            .max(TOLERANCE_ABS);
        let loops = intersect_triangle_soups_eps(&prism_tris, &target_tris, te, te);
        for lp in loops {
            if lp.len() < 2 {
                continue;
            }
            let wb = wire_brep_from_polyline(&lp, te);
            if wb.edges.is_empty() {
                continue;
            }
            out.push(wb);
        }
    }
    out
}


