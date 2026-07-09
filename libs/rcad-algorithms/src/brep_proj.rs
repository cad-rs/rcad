//! Cylindrical rcad_kernel::BRep projection (OCCT `rcad_kernel::BRepProj_Projection` with `gp_Dir`).
//!
//! Strategy (matches OCCT): translate the wire by `-mdis·D`, build a prism by extruding
//! the polyline along `2·mdis·D`, then intersect that prism with the target shape's
//! triangle soup to obtain projected polylines as wires.


use std::collections::HashSet;

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Line3};
use rcad_kernel::topods::{self, TShape};

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
    /// Uniform samples per edge when discretizing edge curves (>= 2).
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

/// Get the first face's outer wire edge chain from TShape iteration.
fn first_outer_wire_edge_chain(brep: &rcad_kernel::BRep) -> Option<Vec<usize>> {
    for ts in &brep.tshapes {
        let TShape::Face(fd) = ts.as_ref() else { continue };
        let wire_ts = brep.tshapes.get(fd.outer_wire.index)?;
        let TShape::Wire(wd) = wire_ts.as_ref() else { continue };
        if wd.edges.is_empty() {
            return None;
        }
        return Some(wd.edges.iter().map(|er| er.index).collect());
    }
    None
}

/// Build vertex-edge adjacency from tshapes.
fn build_vertex_edge_adj(brep: &rcad_kernel::BRep) -> Vec<Vec<usize>> {
    let nv = brep.vertex_count();
    let mut adj = vec![Vec::new(); nv];
    for (ei, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(ed) = ts.as_ref() else { continue };
        let s = ed.first.index;
        let e = ed.last.index;
        if s < adj.len() {
            adj[s].push(ei);
        }
        if e < adj.len() {
            adj[e].push(ei);
        }
    }
    adj
}

/// Greedy edge chaining using vertex-edge adjacency.
fn edge_chains_greedy(brep: &rcad_kernel::BRep) -> Vec<Vec<usize>> {
    let n = brep.edge_count();
    if n == 0 {
        return Vec::new();
    }
    let adj = build_vertex_edge_adj(brep);
    let mut remaining: HashSet<usize> = (0..n).collect();
    let mut chains: Vec<Vec<usize>> = Vec::new();

    while let Some(&start_e) = remaining.iter().next() {
        let ts = &brep.tshapes[start_e];
        let TShape::Edge(ed) = ts.as_ref() else { remaining.remove(&start_e); continue };
        let mut v = ed.first.index;
        let mut next_e = start_e;
        let mut chain: Vec<usize> = Vec::new();

        loop {
            if !remaining.remove(&next_e) {
                break;
            }
            chain.push(next_e);
            let ts = &brep.tshapes[next_e];
            let TShape::Edge(ed) = ts.as_ref() else { break };
            v = if ed.first.index == v { ed.last.index } else { ed.first.index };
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
    let ed = match brep.tshapes.get(ei) {
        Some(ts) => match ts.as_ref() {
            TShape::Edge(ed) => ed,
            _ => return,
        },
        None => return,
    };

    let n_seg = n_seg.max(2);
    if let Some(curve) = &ed.curve {
        let range = ed.range;
        let forward = from_v == ed.first.index;
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
    // Fallback: vertex positions
    let p0 = brep.vertex_point(from_v).unwrap_or(DVec3::ZERO);
    let p1 = brep.vertex_point(to_v).unwrap_or(DVec3::ZERO);
    push_pt(pts, p0, tol);
    push_pt(pts, p1, tol);
}

fn dense_polyline_from_chain(brep: &rcad_kernel::BRep, chain: &[usize], options: &BrepProjOptions) -> Vec<DVec3> {
    if chain.is_empty() {
        return Vec::new();
    }
    let mut pts: Vec<DVec3> = Vec::new();
    // Get first edge's first vertex
    let first_ts = &brep.tshapes[chain[0]];
    let TShape::Edge(first_ed) = first_ts.as_ref() else { return Vec::new() };
    let mut v = first_ed.first.index;
    for &ei in chain {
        let ts = &brep.tshapes[ei];
        let TShape::Edge(ed) = ts.as_ref() else { continue };
        let (from_v, to_v) = if ed.first.index == v {
            (ed.first.index, ed.last.index)
        } else if ed.last.index == v {
            (ed.last.index, ed.first.index)
        } else {
            (ed.first.index, ed.last.index)
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

/// Build a topods::BRep containing only vertices and line edges (no solids), one edge per segment.
fn wire_brep_from_polyline(poly: &[DVec3], tol: f64) -> rcad_kernel::BRep {
    let mut brep = rcad_kernel::BRep::new();
    if poly.len() < 2 {
        return brep;
    }
    let mut vrefs = Vec::new();
    for p in poly {
        vrefs.push(brep.add_tvertex(*p));
    }
    for i in 0..poly.len().saturating_sub(1) {
        let vi_a = vrefs[i];
        let vi_b = vrefs[i + 1];
        let a = poly[i];
        let b = poly[i + 1];
        let d = b - a;
        let len = d.length();
        let dir = if len > tol {
            d / len
        } else {
            DVec3::X
        };
        let curve = Curve3::Line(Line3 {
            origin: a,
            direction: dir,
        });
        brep.add_tedge(Some(curve), vi_a, vi_b, [0.0, len.max(tol)]);
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
            if wb.edge_count() == 0 {
                continue;
            }
            out.push(wb);
        }
    }
    out
}
