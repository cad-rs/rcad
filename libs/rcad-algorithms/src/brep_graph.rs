use std::collections::HashMap;

use rcad_kernel::BRep;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NodeKind {
    Solid,
    Shell,
    Face,
    Wire,
    Edge,
    Vertex,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TopoNode {
    pub kind: NodeKind,
    pub index: usize,
}

#[derive(Debug, Clone)]
pub struct TopoGraphHistoryEvent {
    pub action: String,
}

#[derive(Debug, Clone, Default)]
pub struct TopoGraphHistory {
    pub events: Vec<TopoGraphHistoryEvent>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TopoGraphValidationIssue {
    MissingAdjacency { node: TopoNode },
    NonSymmetricAdjacency { a: TopoNode, b: TopoNode },
    InvalidEdgeVertexRef { edge_index: usize, vertex_index: usize },
}

#[derive(Debug, Clone, Default)]
pub struct TopoGraph {
    pub nodes: Vec<TopoNode>,
    pub history: TopoGraphHistory,
    adjacency: HashMap<TopoNode, Vec<TopoNode>>,
    solid_shells: Vec<Vec<usize>>,
    shell_faces: Vec<Vec<usize>>,
    face_wires: Vec<Vec<usize>>,
    wire_edges: Vec<Vec<usize>>,
    edge_vertices: Vec<[usize; 2]>,
}

impl TopoGraph {
    pub fn from_brep(brep: &BRep) -> Self {
        let mut g = Self::default();
        g.record("from_brep");

        for vi in 0..brep.vertices.len() {
            g.add_node(TopoNode {
                kind: NodeKind::Vertex,
                index: vi,
            });
        }

        for (ei, e) in brep.edges.iter().enumerate() {
            let en = TopoNode {
                kind: NodeKind::Edge,
                index: ei,
            };
            g.add_node(en);
            g.edge_vertices.push([e.start, e.end]);
            g.connect(
                en,
                TopoNode {
                    kind: NodeKind::Vertex,
                    index: e.start,
                },
            );
            g.connect(
                en,
                TopoNode {
                    kind: NodeKind::Vertex,
                    index: e.end,
                },
            );
        }

        let mut shell_idx = 0usize;
        let mut face_idx = 0usize;
        let mut wire_idx = 0usize;

        for (si, solid) in brep.solids.iter().enumerate() {
            let sn = TopoNode {
                kind: NodeKind::Solid,
                index: si,
            };
            g.add_node(sn);
            g.solid_shells.push(Vec::new());

            for shell in &solid.shells {
                let shn = TopoNode {
                    kind: NodeKind::Shell,
                    index: shell_idx,
                };
                g.add_node(shn);
                g.connect(sn, shn);
                g.solid_shells[si].push(shell_idx);
                g.shell_faces.push(Vec::new());

                for face in &shell.faces {
                    let fnn = TopoNode {
                        kind: NodeKind::Face,
                        index: face_idx,
                    };
                    g.add_node(fnn);
                    g.connect(shn, fnn);
                    g.shell_faces[shell_idx].push(face_idx);
                    g.face_wires.push(Vec::new());

                    let wires = std::iter::once(&face.outer_wire).chain(face.inner_wires.iter());
                    for wire in wires {
                        let wn = TopoNode {
                            kind: NodeKind::Wire,
                            index: wire_idx,
                        };
                        g.add_node(wn);
                        g.connect(fnn, wn);
                        g.face_wires[face_idx].push(wire_idx);
                        g.wire_edges.push(Vec::new());

                        for we in &wire.edges {
                            let en = TopoNode {
                                kind: NodeKind::Edge,
                                index: we.idx,
                            };
                            if we.idx < brep.edges.len() {
                                g.connect(wn, en);
                                if !g.wire_edges[wire_idx].contains(&we.idx) {
                                    g.wire_edges[wire_idx].push(we.idx);
                                }
                            }
                        }

                        wire_idx += 1;
                    }

                    face_idx += 1;
                }

                shell_idx += 1;
            }
        }

        g
    }

    pub fn record(&mut self, action: impl Into<String>) {
        self.history.events.push(TopoGraphHistoryEvent {
            action: action.into(),
        });
    }

    pub fn neighbors(&self, node: TopoNode) -> Vec<TopoNode> {
        self.adjacency.get(&node).cloned().unwrap_or_default()
    }

    pub fn faces_of_shell(&self, shell: TopoNode) -> Vec<TopoNode> {
        if shell.kind != NodeKind::Shell {
            return Vec::new();
        }
        self.shell_faces
            .get(shell.index)
            .map(|v| {
                v.iter()
                    .map(|&i| TopoNode {
                        kind: NodeKind::Face,
                        index: i,
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    pub fn edges_of_face(&self, face: TopoNode) -> Vec<TopoNode> {
        if face.kind != NodeKind::Face {
            return Vec::new();
        }
        let mut out: Vec<usize> = Vec::new();
        if let Some(wires) = self.face_wires.get(face.index) {
            for &wi in wires {
                if let Some(edges) = self.wire_edges.get(wi) {
                    for &ei in edges {
                        if !out.contains(&ei) {
                            out.push(ei);
                        }
                    }
                }
            }
        }
        out.into_iter()
            .map(|i| TopoNode {
                kind: NodeKind::Edge,
                index: i,
            })
            .collect()
    }

    pub fn vertices_of_edge(&self, edge: TopoNode) -> Vec<TopoNode> {
        if edge.kind != NodeKind::Edge {
            return Vec::new();
        }
        self.edge_vertices
            .get(edge.index)
            .map(|v| {
                vec![
                    TopoNode {
                        kind: NodeKind::Vertex,
                        index: v[0],
                    },
                    TopoNode {
                        kind: NodeKind::Vertex,
                        index: v[1],
                    },
                ]
            })
            .unwrap_or_default()
    }

    pub fn validate(&self) -> Vec<TopoGraphValidationIssue> {
        let mut issues = Vec::new();
        for node in &self.nodes {
            let Some(neigh) = self.adjacency.get(node) else {
                issues.push(TopoGraphValidationIssue::MissingAdjacency { node: *node });
                continue;
            };
            for n in neigh {
                if let Some(back) = self.adjacency.get(n) {
                    if !back.contains(node) {
                        issues.push(TopoGraphValidationIssue::NonSymmetricAdjacency {
                            a: *node,
                            b: *n,
                        });
                    }
                } else {
                    issues.push(TopoGraphValidationIssue::MissingAdjacency { node: *n });
                }
            }
        }

        for (ei, vv) in self.edge_vertices.iter().enumerate() {
            for &vi in vv {
                if !self.nodes.contains(&TopoNode {
                    kind: NodeKind::Vertex,
                    index: vi,
                }) {
                    issues.push(TopoGraphValidationIssue::InvalidEdgeVertexRef {
                        edge_index: ei,
                        vertex_index: vi,
                    });
                }
            }
        }

        issues
    }

    /// Compact graph storage by dropping orphan adjacency entries and
    /// deduplicating neighbor lists.
    pub fn compact(&mut self) {
        let mut node_set = std::collections::HashSet::new();
        for n in &self.nodes {
            node_set.insert(*n);
        }

        self.adjacency.retain(|node, _| node_set.contains(node));
        for neigh in self.adjacency.values_mut() {
            neigh.retain(|n| node_set.contains(n));
            neigh.sort_by_key(|n| {
                let kind_rank = match n.kind {
                    NodeKind::Solid => 0usize,
                    NodeKind::Shell => 1,
                    NodeKind::Face => 2,
                    NodeKind::Wire => 3,
                    NodeKind::Edge => 4,
                    NodeKind::Vertex => 5,
                };
                (kind_rank, n.index)
            });
            neigh.dedup();
        }
        self.record("compact");
    }

    /// Apply a mutation and run graph validation afterward.
    ///
    /// This is a lightweight baseline for mutation-guard workflows: callers
    /// can route all topology edits through this helper and reject invalid
    /// states before continuing downstream processing.
    pub fn mutate_checked<F>(
        &mut self,
        action: impl Into<String>,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph),
    {
        let action = action.into();
        mutator(self);
        let issues = self.validate();
        if issues.is_empty() {
            self.record(format!("mutate:{action}"));
            Ok(())
        } else {
            self.record(format!("mutate_invalid:{action}"));
            Err(issues)
        }
    }

    /// Apply a mutation with rollback-on-failure semantics.
    ///
    /// If validation fails after the mutation, graph state is restored to the
    /// pre-mutation snapshot and validation issues are returned.
    pub fn mutate_guarded<F>(
        &mut self,
        action: impl Into<String>,
        mutator: F,
    ) -> Result<(), Vec<TopoGraphValidationIssue>>
    where
        F: FnOnce(&mut TopoGraph),
    {
        let action = action.into();
        let before = self.clone();
        mutator(self);
        let issues = self.validate();
        if issues.is_empty() {
            self.record(format!("mutate_guarded:{action}"));
            Ok(())
        } else {
            *self = before;
            self.record(format!("mutate_guarded_rollback:{action}"));
            Err(issues)
        }
    }

    fn add_node(&mut self, node: TopoNode) {
        if !self.nodes.contains(&node) {
            self.nodes.push(node);
        }
        self.adjacency.entry(node).or_default();
    }

    fn connect(&mut self, a: TopoNode, b: TopoNode) {
        self.add_node(a);
        self.add_node(b);
        let va = self.adjacency.entry(a).or_default();
        if !va.contains(&b) {
            va.push(b);
        }
        let vb = self.adjacency.entry(b).or_default();
        if !vb.contains(&a) {
            vb.push(a);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn topo_graph_from_box_has_expected_counts() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);

        let solids = g.nodes.iter().filter(|n| n.kind == NodeKind::Solid).count();
        let shells = g.nodes.iter().filter(|n| n.kind == NodeKind::Shell).count();
        let faces = g.nodes.iter().filter(|n| n.kind == NodeKind::Face).count();
        let wires = g.nodes.iter().filter(|n| n.kind == NodeKind::Wire).count();
        let edges = g.nodes.iter().filter(|n| n.kind == NodeKind::Edge).count();
        let vertices = g.nodes.iter().filter(|n| n.kind == NodeKind::Vertex).count();

        assert_eq!(solids, 1);
        assert_eq!(shells, 1);
        assert_eq!(faces, 6);
        assert_eq!(wires, 6);
        assert_eq!(edges, 12);
        assert_eq!(vertices, 8);
    }

    #[test]
    fn topo_graph_faces_and_edges_queries_work() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);

        let shell0 = TopoNode {
            kind: NodeKind::Shell,
            index: 0,
        };
        let faces = g.faces_of_shell(shell0);
        assert_eq!(faces.len(), 6);

        let face0 = TopoNode {
            kind: NodeKind::Face,
            index: 0,
        };
        let edges = g.edges_of_face(face0);
        assert_eq!(edges.len(), 4);

        let edge0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        let verts = g.vertices_of_edge(edge0);
        assert_eq!(verts.len(), 2);
        assert_eq!(verts[0].index, brep.edges[0].start);
        assert_eq!(verts[1].index, brep.edges[0].end);
    }

    #[test]
    fn topo_graph_validate_passes_on_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let g = TopoGraph::from_brep(&brep);
        assert!(g.validate().is_empty());
        assert!(!g.history.events.is_empty());
    }

    #[test]
    fn topo_graph_compact_drops_orphans_and_dedups_neighbors() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);

        // Inject duplicate and orphan adjacency entries to simulate noisy edits.
        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        if let Some(neigh) = g.adjacency.get_mut(&v0) {
            neigh.push(e0);
            neigh.push(e0);
        }
        let orphan = TopoNode {
            kind: NodeKind::Face,
            index: 9999,
        };
        g.adjacency.insert(orphan, vec![v0]);

        g.compact();

        assert!(!g.adjacency.contains_key(&orphan));
        let neigh = g.adjacency.get(&v0).expect("vertex adjacency exists");
        let count_e0 = neigh.iter().filter(|n| **n == e0).count();
        assert_eq!(count_e0, 1);
        assert!(g.history.events.iter().any(|e| e.action == "compact"));
    }

    #[test]
    fn topo_graph_mutate_checked_reports_invalid_graph() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);

        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };

        let res = g.mutate_checked("inject_nonsymmetric", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.retain(|n| *n != v0);
            }
        });

        assert!(res.is_err());
        let issues = res.expect_err("mutation should be invalid");
        assert!(issues
            .iter()
            .any(|i| matches!(i, TopoGraphValidationIssue::NonSymmetricAdjacency { .. })));
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_invalid:inject_nonsymmetric"));
    }

    #[test]
    fn topo_graph_mutate_guarded_rolls_back_on_invalid_mutation() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);
        let before = g.clone();

        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };

        let res = g.mutate_guarded("inject_nonsymmetric", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.retain(|n| *n != v0);
            }
        });

        assert!(res.is_err());
        assert_eq!(g.nodes, before.nodes);
        assert_eq!(g.adjacency, before.adjacency);
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_guarded_rollback:inject_nonsymmetric"));
    }

    #[test]
    fn topo_graph_mutate_guarded_commits_valid_mutation() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut g = TopoGraph::from_brep(&brep);
        let e0 = TopoNode {
            kind: NodeKind::Edge,
            index: 0,
        };
        let v0 = TopoNode {
            kind: NodeKind::Vertex,
            index: 0,
        };

        let res = g.mutate_guarded("dedup_neighbors", |graph| {
            if let Some(neigh) = graph.adjacency.get_mut(&e0) {
                neigh.push(v0);
            }
            graph.compact();
        });

        assert!(res.is_ok());
        let neigh = g.adjacency.get(&e0).expect("edge adjacency exists");
        let count_v0 = neigh.iter().filter(|n| **n == v0).count();
        assert_eq!(count_v0, 1);
        assert!(g
            .history
            .events
            .iter()
            .any(|e| e.action == "mutate_guarded:dedup_neighbors"));
    }
}
