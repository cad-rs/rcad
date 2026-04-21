use crate::{BRep, Edge, WireEdge};

fn oriented_edge_vertices(brep: &BRep, we: WireEdge) -> Option<(usize, usize)> {
    let e = brep.edges.get(we.idx)?;
    if we.forward {
        Some((e.start, e.end))
    } else {
        Some((e.end, e.start))
    }
}

fn points_are_collinear_forward(
    a: glam::DVec3,
    b: glam::DVec3,
    c: glam::DVec3,
    linear_tol: f64,
) -> bool {
    let ab = b - a;
    let bc = c - b;
    let ab_len = ab.length();
    let bc_len = bc.length();
    if ab_len <= 1e-12 || bc_len <= 1e-12 {
        return false;
    }
    let cross = ab.cross(bc).length();
    let dot = ab.dot(bc);
    cross <= linear_tol * (ab_len + bc_len) && dot > 0.0
}

fn find_existing_edge_between_vertices(brep: &BRep, from: usize, to: usize) -> Option<WireEdge> {
    for (idx, e) in brep.edges.iter().enumerate() {
        if e.start == from && e.end == to {
            return Some(WireEdge::fwd(idx));
        }
        if e.start == to && e.end == from {
            return Some(WireEdge::rev(idx));
        }
    }
    None
}

fn push_topology_edge(brep: &mut BRep, start: usize, end: usize) -> WireEdge {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start, end });

    if !brep.geom.edge_curve.is_empty() {
        brep.geom.edge_curve.push(None);
    }
    if !brep.geom.edge_pcurves.is_empty() {
        brep.geom.edge_pcurves.push(Vec::new());
    }
    if !brep.geom.edge_curve_range.is_empty() {
        brep.geom.edge_curve_range.push(None);
    }
    if !brep.geom.edge_degenerated.is_empty() {
        brep.geom.edge_degenerated.push(false);
    }
    if !brep.geom.edge_tolerance.is_empty() {
        brep.geom.edge_tolerance.push(0.0);
    }
    if !brep.geom.edge_same_parameter.is_empty() {
        brep.geom.edge_same_parameter.push(false);
    }
    if !brep.geom.edge_same_range.is_empty() {
        brep.geom.edge_same_range.push(false);
    }

    WireEdge::fwd(idx)
}

fn merge_collinear_consecutive_edges_in_wire(
    brep: &mut BRep,
    wire: &mut Vec<WireEdge>,
    linear_tol: f64,
) -> usize {
    if wire.len() < 4 {
        return 0;
    }

    let mut merged = 0usize;
    loop {
        let n = wire.len();
        if n < 4 {
            break;
        }

        let mut changed = false;
        for i in 0..n {
            let j = (i + 1) % n;
            let Some((u, v1)) = oriented_edge_vertices(brep, wire[i]) else {
                continue;
            };
            let Some((v2, w)) = oriented_edge_vertices(brep, wire[j]) else {
                continue;
            };
            if v1 != v2 || u == w {
                continue;
            }

            let Some(p_u) = brep.vertices.get(u).map(|v| v.point) else {
                continue;
            };
            let Some(p_v) = brep.vertices.get(v1).map(|v| v.point) else {
                continue;
            };
            let Some(p_w) = brep.vertices.get(w).map(|v| v.point) else {
                continue;
            };
            if !points_are_collinear_forward(p_u, p_v, p_w, linear_tol) {
                continue;
            }

            let bridge = find_existing_edge_between_vertices(brep, u, w)
                .unwrap_or_else(|| push_topology_edge(brep, u, w));

            if i + 1 < n {
                wire.splice(i..=i + 1, [bridge]);
            } else {
                wire.pop();
                wire.remove(0);
                wire.insert(0, bridge);
            }
            merged += 1;
            changed = true;
            break;
        }

        if !changed {
            break;
        }
    }
    merged
}

/// Merge fragmented collinear edges inside all face wires (outer + inner).
///
/// This pass is intentionally local/topological: it only collapses *consecutive*
/// collinear edge pairs in a wire into a direct bridge edge, preserving wire
/// order and loop closure.
///
/// Returns the number of edge-pair merges performed.
pub fn merge_collinear_edges_in_wires(brep: &mut BRep, linear_tol: f64) -> usize {
    let tol = linear_tol.max(1e-12);
    let mut merged_total = 0usize;

    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            for fi in 0..brep.solids[si].shells[shi].faces.len() {
                let mut outer = brep.solids[si].shells[shi].faces[fi].outer_wire.edges.clone();
                merged_total += merge_collinear_consecutive_edges_in_wire(brep, &mut outer, tol);
                brep.solids[si].shells[shi].faces[fi].outer_wire.edges = outer;

                let inner_count = brep.solids[si].shells[shi].faces[fi].inner_wires.len();
                for wi in 0..inner_count {
                    let mut inner =
                        brep.solids[si].shells[shi].faces[fi].inner_wires[wi].edges.clone();
                    merged_total += merge_collinear_consecutive_edges_in_wire(brep, &mut inner, tol);
                    brep.solids[si].shells[shi].faces[fi].inner_wires[wi].edges = inner;
                }
            }
        }
    }

    merged_total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{Face, Shell, Solid, Vertex, Wire};

    #[test]
    fn merge_collinear_edges_in_wires_collapses_fragmented_chain() {
        let mut brep = BRep::new();
        brep.vertices = vec![
            Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) }, // 0
            Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) }, // 1 split point
            Vertex { point: glam::DVec3::new(2.0, 0.0, 0.0) }, // 2
            Vertex { point: glam::DVec3::new(2.0, 1.0, 0.0) }, // 3
            Vertex { point: glam::DVec3::new(0.0, 1.0, 0.0) }, // 4
        ];
        brep.edges = vec![
            Edge { start: 0, end: 1 }, // e0
            Edge { start: 1, end: 2 }, // e1
            Edge { start: 2, end: 3 }, // e2
            Edge { start: 3, end: 4 }, // e3
            Edge { start: 4, end: 0 }, // e4
        ];
        brep.solids = vec![Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![
                            WireEdge::fwd(0),
                            WireEdge::fwd(1),
                            WireEdge::fwd(2),
                            WireEdge::fwd(3),
                            WireEdge::fwd(4),
                        ],
                    },
                    inner_wires: vec![],
                    normal: glam::DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        }];

        let merged = merge_collinear_edges_in_wires(&mut brep, 1e-6);
        assert!(merged >= 1);
        let edge_count = brep.solids[0].shells[0].faces[0].outer_wire.edges.len();
        assert_eq!(edge_count, 4);
    }
}
