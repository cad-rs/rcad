//! Example: Non-manifold topology analysis with BRepGraph.
//!
//! Demonstrates:
//!   1. Building a BRepGraph from a manifold BRep
//!   2. O(1) adjacency queries
//!   3. Manifold detection
//!   4. Non-manifold topology analysis
//!   5. Repair hints generation
//!   6. Graph traversal (DFS/BFS)
//!   7. Mutation tracking and history
//!
//! Run:
//!   cargo run -p rcad-examples --example non_manifold_topology

use glam::DVec3;
use rcad_kernel::{
    BRep, BRepGraph, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge,
    PrimitiveSolid,
    brep_graph::{RepairHint, BRepGraphHistory},
};

fn separator(title: &str) {
    println!("\n──────────────────────────────────────────");
    println!("  {title}");
    println!("──────────────────────────────────────────");
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Building BRepGraph from manifold BRep
// ─────────────────────────────────────────────────────────────────────────────

fn demo_build_graph() {
    separator("1. Building BRepGraph from Manifold BRep");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    let graph = BRepGraph::from_brep(&brep);

    println!("  Box BRep statistics:");
    println!("    Vertices: {}", graph.vertex_count);
    println!("    Edges: {}", graph.edge_count);
    println!("    Faces: {}", graph.face_count);

    assert_eq!(graph.vertex_count, 8, "box has 8 vertices");
    assert_eq!(graph.edge_count, 12, "box has 12 edges");
    assert_eq!(graph.face_count, 6, "box has 6 faces");

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. O(1) Adjacency queries
// ─────────────────────────────────────────────────────────────────────────────

fn demo_adjacency_queries() {
    separator("2. O(1) Adjacency Queries");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    let graph = BRepGraph::from_brep(&brep);

    // Edge -> adjacent faces
    {
        let adj = graph.edge_adjacent_faces(0);
        println!("  Edge 0 adjacent faces: {:?}", adj);
        assert_eq!(adj.len(), 2, "manifold edge has 2 adjacent faces");
    }

    // Face -> edges
    {
        let edges = graph.face_edges(0);
        println!("  Face 0 edges: {:?}", edges);
        assert_eq!(edges.len(), 4, "box face has 4 edges");
    }

    // Vertex -> edges
    {
        let adj_edges = graph.vertex_adjacent_edges(0);
        println!("  Vertex 0 adjacent edges: {:?}", adj_edges);
        assert_eq!(adj_edges.len(), 3, "box vertex has 3 adjacent edges");
    }

    // Vertex -> faces
    {
        let adj_faces = graph.vertex_adjacent_faces(0);
        println!("  Vertex 0 adjacent faces: {:?}", adj_faces);
        assert_eq!(adj_faces.len(), 3, "box vertex touches 3 faces");
    }

    // Edge endpoints
    {
        let (start, end) = graph.edge_endpoints(0).unwrap();
        println!("  Edge 0 endpoints: ({}, {})", start, end);
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Manifold detection
// ─────────────────────────────────────────────────────────────────────────────

fn demo_manifold_detection() {
    separator("3. Manifold Detection");

    // Closed manifold box
    {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let graph = BRepGraph::from_brep(&brep);

        assert!(graph.is_manifold(), "closed box should be manifold");
        assert!(graph.is_closed(), "closed box should have no boundary edges");

        let nm_edges = graph.non_manifold_edges();
        let boundary = graph.boundary_edges();
        let orphan = graph.orphan_edges();

        println!("  Closed box:");
        println!("    Is manifold: {}", graph.is_manifold());
        println!("    Is closed: {}", graph.is_closed());
        println!("    Non-manifold edges: {}", nm_edges.len());
        println!("    Boundary edges: {}", boundary.len());
        println!("    Orphan edges: {}", orphan.len());
    }

    // Open shell (missing one face)
    {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        // Remove one face
        if let Some(s) = brep.solids.first_mut() {
            if let Some(sh) = s.shells.first_mut() {
                sh.faces.pop();
            }
        }

        let graph = BRepGraph::from_brep(&brep);

        println!("\n  Open box (missing one face):");
        println!("    Is manifold: {}", graph.is_manifold());
        println!("    Is closed: {}", graph.is_closed());

        let boundary = graph.boundary_edges();
        println!("    Boundary edges: {:?}", boundary);
        assert!(!boundary.is_empty(), "open shell should have boundary edges");
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Non-manifold topology analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Build a minimal non-manifold BRep where edge 0 is shared by 3 faces.
fn build_tripod() -> BRep {
    let vertices = vec![
        Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
        Vertex { point: DVec3::new(1.0, 0.0, 0.0) }, // 1
        Vertex { point: DVec3::new(0.0, 1.0, 0.0) }, // 2
        Vertex { point: DVec3::new(0.0, 0.0, 1.0) }, // 3
        Vertex { point: DVec3::new(0.0, -1.0, 0.0) }, // 4
    ];

    let edges = vec![
        Edge { start: 0, end: 1 }, // 0: shared spine
        Edge { start: 1, end: 2 }, // 1
        Edge { start: 2, end: 0 }, // 2
        Edge { start: 1, end: 3 }, // 3
        Edge { start: 3, end: 0 }, // 4
        Edge { start: 1, end: 4 }, // 5
        Edge { start: 4, end: 0 }, // 6
    ];

    let f0 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
        },
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles: vec![],
        mesh_dirty: true,
    };
    let f1 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(3), WireEdge::fwd(4)],
        },
        inner_wires: vec![],
        normal: DVec3::Y,
        triangles: vec![],
        mesh_dirty: true,
    };
    let f2 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(5), WireEdge::fwd(6)],
        },
        inner_wires: vec![],
        normal: -DVec3::Y,
        triangles: vec![],
        mesh_dirty: true,
    };

    BRep {
        vertices,
        edges,
        solids: vec![Solid {
            shells: vec![Shell {
                faces: vec![f0, f1, f2],
            }],
        }],
        geom: Default::default(),
    }
}

fn demo_non_manifold_analysis() {
    separator("4. Non-Manifold Topology Analysis");

    let brep = build_tripod();
    let graph = BRepGraph::from_brep(&brep);

    println!("  Tripod (non-manifold) statistics:");
    println!("    Vertices: {}", graph.vertex_count);
    println!("    Edges: {}", graph.edge_count);
    println!("    Faces: {}", graph.face_count);

    // Edge valence
    println!("\n  Edge valences:");
    for ei in 0..graph.edge_count {
        let valence = graph.edge_valence(ei);
        println!("    Edge {}: valence = {}", ei, valence);
    }

    // Non-manifold edges
    let multi_face = graph.multi_face_edges();
    println!("\n  Multi-face edges: {:?}", multi_face);
    assert_eq!(multi_face, vec![0], "edge 0 should be shared by 3 faces");

    // Boundary edges
    let boundary = graph.boundary_edges();
    println!("  Boundary edges: {:?}", boundary);

    // Non-manifold vertices
    let nm_vertices = graph.non_manifold_vertices();
    println!("  Non-manifold vertices: {:?}", nm_vertices);

    // Vertex degrees
    println!("\n  Vertex degrees:");
    for vi in 0..graph.vertex_count {
        println!("    Vertex {}: degree = {}", vi, graph.vertex_degree(vi));
    }

    assert!(!graph.is_manifold(), "tripod should be non-manifold");

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Repair hints
// ─────────────────────────────────────────────────────────────────────────────

fn demo_repair_hints() {
    separator("5. Repair Hints");

    let brep = build_tripod();
    let graph = BRepGraph::from_brep(&brep);
    let hints = graph.repair_hints(&brep);

    println!("  Repair hints for tripod:");
    for hint in &hints.hints {
        match hint {
            RepairHint::MultiManifoldEdge { edge_idx, face_count } => {
                println!("    MultiManifoldEdge: edge {} shared by {} faces", edge_idx, face_count);
            }
            RepairHint::NonManifoldVertex { vertex_idx, connected_multi_edges } => {
                println!("    NonManifoldVertex: vertex {} connected to multi-edges {:?}",
                    vertex_idx, connected_multi_edges);
            }
            RepairHint::UnmatchedBoundaryEdge { edge_idx, face_idx } => {
                println!("    UnmatchedBoundaryEdge: edge {} on face {}", edge_idx, face_idx);
            }
            RepairHint::StitchablePair { edge_a, edge_b, face_a, face_b } => {
                println!("    StitchablePair: edges {} and {} (faces {} and {})",
                    edge_a, edge_b, face_a, face_b);
            }
            RepairHint::OrphanEdge { edge_idx } => {
                println!("    OrphanEdge: edge {}", edge_idx);
            }
        }
    }

    // Check for expected hints
    let has_multi = hints.hints.iter().any(|h| {
        matches!(h, RepairHint::MultiManifoldEdge { edge_idx: 0, .. })
    });
    assert!(has_multi, "should have MultiManifoldEdge hint for edge 0");

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Graph traversal (DFS/BFS)
// ─────────────────────────────────────────────────────────────────────────────

fn demo_traversal() {
    separator("6. Graph Traversal (DFS/BFS)");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    let graph = BRepGraph::from_brep(&brep);

    // BFS faces
    {
        let mut visited: Vec<usize> = graph.bfs_faces(0).collect();
        visited.sort_unstable();
        println!("  BFS from face 0: {:?}", visited);
        assert_eq!(visited.len(), 6, "should visit all 6 faces");
    }

    // DFS faces
    {
        let mut visited: Vec<usize> = graph.dfs_faces(0).collect();
        visited.sort_unstable();
        println!("  DFS from face 0: {:?}", visited);
        assert_eq!(visited.len(), 6, "should visit all 6 faces");
    }

    // DFS edges from vertex
    {
        let mut visited: Vec<usize> = graph.dfs_edges_from_vertex(0).collect();
        visited.sort_unstable();
        println!("  DFS edges from vertex 0: {:?}", visited);
        assert_eq!(visited.len(), 12, "should visit all 12 edges");
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. Mutation tracking and history
// ─────────────────────────────────────────────────────────────────────────────

fn demo_mutation_tracking() {
    separator("7. Mutation Tracking and History");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    let mut graph = BRepGraph::from_brep(&brep);

    // RAII mutation guard
    {
        let mut guard = graph.begin_mutation();
        guard.graph().mark_face_modified(0);
        guard.graph().mark_edge_modified(1);
        guard.graph().mark_vertex_modified(2);

        println!("  During mutation:");
        println!("    Modified faces: {:?}", guard.graph().modified_faces());
        println!("    Modified edges: {:?}", guard.graph().modified_edges());
        println!("    Modified vertices: {:?}", guard.graph().modified_vertices());

        // Commit with history
        let mut history = BRepGraphHistory::new();
        guard.commit_with_history(&mut history, Some("test_mutation".to_string()))
            .expect("commit should succeed");

        println!("  History events: {}", history.len());
        let event = history.last().unwrap();
        println!("    Label: {:?}", event.label);
        println!("    Topology changed: {}", event.topology_changed);
    }

    // Check that changes persisted
    println!("\n  After commit:");
    println!("    Modified faces: {:?}", graph.modified_faces());
    println!("    Modified edges: {:?}", graph.modified_edges());

    // Checkpoint and rollback
    {
        let checkpoint = graph.checkpoint();
        graph.mark_face_modified(5);
        graph.mark_edge_modified(7);
        println!("\n  Before rollback:");
        println!("    Modified faces: {:?}", graph.modified_faces());

        graph.restore_from_checkpoint(&checkpoint);
        println!("  After rollback:");
        println!("    Modified faces: {:?}", graph.modified_faces());
        assert!(!graph.modified_faces().contains(&5), "face 5 should not be modified after rollback");
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    println!("=================================================");
    println!("  Non-Manifold Topology Analysis Demo");
    println!("=================================================");

    demo_build_graph();
    demo_adjacency_queries();
    demo_manifold_detection();
    demo_non_manifold_analysis();
    demo_repair_hints();
    demo_traversal();
    demo_mutation_tracking();

    println!("\n=================================================");
    println!("  Non-Manifold Topology: All demos completed successfully");
    println!("=================================================");
}
