//! Example: Phase G — Topology Query API and Curvature Analysis.
//!
//! Demonstrates:
//!   1. Topology queries: edge adjacency, vertex adjacency, shape counts on a box
//!   2. Analytic curvature: cylinder, sphere, torus principal/Gaussian/mean curvatures
//!   3. Numerical BSpline curvature: flat bilinear patch → K ≈ 0, H ≈ 0
//!
//! Run: cargo run --example phase_g_demo

use glam::DVec3;
use rcad_kernel::{
    curvature::{gaussian_curvature, mean_curvature, principal_curvatures},
    geom::{BSplineSurface, Surface3},
    topo_query::{
        edge_adjacent_faces, edge_count, face_count, face_edges, vertex_adjacent_edges,
        vertex_count,
    },
};
use rcad_modeling::{box_brep, cylinder_brep, sphere_brep, torus_brep};

// ── 1. Topology queries ───────────────────────────────────────────────────────

fn demo_topo_queries() {
    println!("\n=== 1. Topology Queries: 2×2×2 box ===");

    let brep = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

    println!(
        "  faces={}, edges={}, vertices={}",
        face_count(&brep),
        edge_count(&brep),
        vertex_count(&brep)
    );

    // Every edge of a closed convex solid must be adjacent to exactly 2 faces
    let mut bad_edges = Vec::new();
    for ei in 0..edge_count(&brep) {
        let adj = edge_adjacent_faces(&brep, ei);
        if adj.len() != 2 {
            bad_edges.push((ei, adj.len()));
        }
    }
    if bad_edges.is_empty() {
        println!("  ✓ All {} edges have exactly 2 adjacent faces", edge_count(&brep));
    } else {
        println!("  ✗ Edges with wrong adjacency count: {:?}", bad_edges);
    }

    // Every vertex of a box is shared by exactly 3 edges
    let mut bad_verts = Vec::new();
    for vi in 0..vertex_count(&brep) {
        let adj = vertex_adjacent_edges(&brep, vi);
        if adj.len() != 3 {
            bad_verts.push((vi, adj.len()));
        }
    }
    if bad_verts.is_empty() {
        println!("  ✓ All {} vertices have exactly 3 adjacent edges", vertex_count(&brep));
    } else {
        println!("  ✗ Vertices with wrong edge count: {:?}", bad_verts);
    }

    // Face 0 edge list
    let fe = face_edges(&brep, 0);
    println!("  Face 0 outer_wire edge indices: {:?}", fe);
}

// ── 2. Analytic curvature ─────────────────────────────────────────────────────

fn print_curvature(label: &str, surface: &Surface3, u: f64, v: f64) {
    let (k1, k2) = principal_curvatures(surface, u, v);
    let k = gaussian_curvature(surface, u, v);
    let h = mean_curvature(surface, u, v);
    println!("  {label:30}  k1={k1:8.4}  k2={k2:8.4}  K={k:8.4}  H={h:8.4}");
}

fn demo_analytic_curvature() {
    println!("\n=== 2. Analytic Curvature ===");

    // Cylinder r=1, h=3  (K=0, H=0.5 everywhere)
    let cyl_brep = cylinder_brep(DVec3::ZERO, DVec3::Y, DVec3::X, 1.0, 3.0).unwrap();
    if let Some(s) = cyl_brep.geom.surfaces.iter().find(|s| matches!(s, Surface3::Cylinder(_))) {
        print_curvature("Cylinder(r=1) at (0,0)", s, 0.0, 0.0);
        print_curvature("Cylinder(r=1) at (π/2,1)", s, std::f64::consts::FRAC_PI_2, 1.0);
    }

    // Sphere r=1  (K=1, H=1 everywhere)
    let sph_brep = sphere_brep(DVec3::ZERO, 1.0).unwrap();
    if let Some(s) = sph_brep.geom.surfaces.iter().find(|s| matches!(s, Surface3::Sphere(_))) {
        print_curvature("Sphere(r=1) at (0,π/2)", s, 0.0, std::f64::consts::FRAC_PI_2);
    }

    // Sphere r=2  (K=0.25, H=0.5)
    let sph2_brep = sphere_brep(DVec3::ZERO, 2.0).unwrap();
    if let Some(s) = sph2_brep.geom.surfaces.iter().find(|s| matches!(s, Surface3::Sphere(_))) {
        print_curvature("Sphere(r=2) at (0,π/2)", s, 0.0, std::f64::consts::FRAC_PI_2);
    }

    // Torus R=2, r=0.5  (outer equator v=0)
    //   k_tube = 1/r = 2.0,  k_major = 1/(R+r) = 1/2.5 = 0.4
    //   K = 0.8,  H = 1.2
    let tor_brep = torus_brep(DVec3::ZERO, DVec3::Y, DVec3::X, 2.0, 0.5).unwrap();
    if let Some(s) = tor_brep.geom.surfaces.iter().find(|s| matches!(s, Surface3::Torus(_))) {
        print_curvature("Torus(R=2,r=0.5) outer (v=0)", s, 0.0, 0.0);
        print_curvature("Torus(R=2,r=0.5) inner (v=π)", s, 0.0, std::f64::consts::PI);
        print_curvature("Torus(R=2,r=0.5) side  (v=π/2)", s, 0.0, std::f64::consts::FRAC_PI_2);
    }
}

// ── 3. BSpline numerical curvature ───────────────────────────────────────────

fn demo_bspline_curvature() {
    println!("\n=== 3. BSpline Numerical Curvature ===");

    // Flat bilinear patch in XY plane — K≈0, H≈0
    let flat_patch = Surface3::BSpline(BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
        ],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    });
    print_curvature("BSpline flat patch at (0.5,0.5)", &flat_patch, 0.5, 0.5);

    // Saddle surface z = x*y, bilinear patch lifted at one corner
    let saddle_patch = Surface3::BSpline(BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 1.0)], // corner at z=1
        ],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    });
    print_curvature("BSpline saddle-like at (0.5,0.5)", &saddle_patch, 0.5, 0.5);
    println!("  (saddle: K<0 expected for the lifted-corner patch at center)");
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔═══════════════════════════════════════════════════╗");
    println!("║              RCAD Phase G Demo                    ║");
    println!("║   Topology Queries · Curvature Analysis           ║");
    println!("╚═══════════════════════════════════════════════════╝");

    demo_topo_queries();
    demo_analytic_curvature();
    demo_bspline_curvature();

    println!("\n✓ Phase G demo complete.");
}
