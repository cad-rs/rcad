//! Example: Phase I — 2D B-Spline PCurve + Per-Entity Tolerance System.
//!
//! Demonstrates:
//!   1. BSplineCurve2 evaluation: degree-1 and degree-3 curves in UV parameter space
//!   2. Curve2d::BSpline stored in GeomStore and used as a PCurve
//!   3. Per-entity tolerance queries: default CONFUSION fallback, custom values,
//!      model_tolerance returning the maximum
//!
//! Run: cargo run --example phase_i_demo

use glam::DVec2;
use rcad_kernel::{
    geom::{BSplineCurve2, Curve2d, Line2d},
    Curve2dEval,
    BRep, PrimitiveSolid,
    CONFUSION, ANGULAR, APPROXIMATION,
    vertex_tolerance, edge_tolerance, face_tolerance, model_tolerance,
};

// ── 1. BSplineCurve2 evaluation ───────────────────────────────────────────────

fn demo_bspline_curve2() {
    println!("\n=== 1. BSplineCurve2 Evaluation ===");

    // Degree-1 (linear) BSpline: straight line (0,0) → (1,0) → (1,1)
    // knots: [0,0,0.5,1,1]  control: [(0,0), (1,0), (1,1)]
    let linear = BSplineCurve2 {
        degree: 1,
        knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
        control_points: vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(1.0, 1.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let p0 = linear.point_at(0.0);
    let p_mid = linear.point_at(0.5);
    let p1 = linear.point_at(1.0);
    println!("  Linear BSpline2 (0,0)→(1,0)→(1,1):");
    println!("    t=0.0  → ({:.4}, {:.4})  (expect 0.0, 0.0)", p0.x, p0.y);
    println!("    t=0.5  → ({:.4}, {:.4})  (expect 1.0, 0.0)", p_mid.x, p_mid.y);
    println!("    t=1.0  → ({:.4}, {:.4})  (expect 1.0, 1.0)", p1.x, p1.y);

    // Degree-3 cubic: 4 control points, uniform knots
    // Bezier-like: ctrl = [(0,0), (0,1), (1,1), (1,0)] — S-curve in UV space
    let cubic = BSplineCurve2 {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(0.0, 1.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(1.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let c0 = cubic.point_at(0.0);
    let c_half = cubic.point_at(0.5);
    let c1 = cubic.point_at(1.0);
    println!("  Cubic BSpline2 Bezier arc [(0,0),(0,1),(1,1),(1,0)]:");
    println!("    t=0.0  → ({:.4}, {:.4})  (expect 0.0, 0.0)", c0.x, c0.y);
    println!("    t=0.5  → ({:.4}, {:.4})  (midpoint by symmetry)", c_half.x, c_half.y);
    println!("    t=1.0  → ({:.4}, {:.4})  (expect 1.0, 0.0)", c1.x, c1.y);
}

// ── 2. Curve2d::BSpline in GeomStore ─────────────────────────────────────────

fn demo_curve2d_in_geomstore() {
    println!("\n=== 2. Curve2d::BSpline in GeomStore ===");

    let mut brep = BRep::new();

    // Push three 2D curves into the pool: Line, Circle (existing variants),
    // and a BSpline (new variant)
    brep.geom.curve2ds.push(Curve2d::Line(Line2d {
        origin: DVec2::ZERO,
        direction: DVec2::X,
    }));

    brep.geom.curve2ds.push(Curve2d::BSpline(BSplineCurve2 {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(0.5, 1.0),
            DVec2::new(1.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    }));

    println!("  curve2ds pool size: {}  (expect 2)", brep.geom.curve2ds.len());

    // Evaluate the BSpline PCurve at various parameters
    let bspline_pcurve = &brep.geom.curve2ds[1];
    let mid = bspline_pcurve.point_at(0.5);
    let start = bspline_pcurve.point_at(0.0);
    let end = bspline_pcurve.point_at(1.0);
    println!("  Quadratic PCurve [(0,0),(0.5,1),(1,0)]:");
    println!("    t=0.0  → ({:.4}, {:.4})  (expect 0.0, 0.0)", start.x, start.y);
    println!("    t=0.5  → ({:.4}, {:.4})  (expect 0.5, 0.5 — apex at midpoint)", mid.x, mid.y);
    println!("    t=1.0  → ({:.4}, {:.4})  (expect 1.0, 0.0)", end.x, end.y);

    println!("  ✓ Curve2d::BSpline stored and evaluated via GeomStore");
}

// ── 3. Per-entity tolerance ───────────────────────────────────────────────────

fn demo_tolerance() {
    println!("\n=== 3. Per-Entity Tolerance System ===");

    println!("  Precision constants:");
    println!("    CONFUSION    = {:.0e}  (point coincidence, OCCT default)", CONFUSION);
    println!("    ANGULAR      = {:.0e}  (angular, radians)", ANGULAR);
    println!("    APPROXIMATION= {:.0e}  (tessellation approximation)", APPROXIMATION);

    // BRep with no stored tolerances → all queries return CONFUSION
    let brep_empty = BRep::new();
    println!("\n  Empty BRep (no stored tolerances):");
    println!("    vertex_tolerance(0) = {:.2e}  (expect {:.2e})", vertex_tolerance(&brep_empty, 0), CONFUSION);
    println!("    edge_tolerance(0)   = {:.2e}  (expect {:.2e})", edge_tolerance(&brep_empty, 0), CONFUSION);
    println!("    face_tolerance(0)   = {:.2e}  (expect {:.2e})", face_tolerance(&brep_empty, 0), CONFUSION);
    println!("    model_tolerance     = {:.2e}  (expect {:.2e})", model_tolerance(&brep_empty), CONFUSION);

    // Box BRep with custom per-edge tolerance
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();
    let n_faces = brep.solids[0].shells[0].faces.len();

    // Populate tolerances: vertices at 1e-6, edges at 1e-5, faces at 1e-4
    brep.geom.vertex_tolerance = vec![1e-6; n_verts];
    brep.geom.edge_tolerance   = vec![1e-5; n_edges];
    brep.geom.face_tolerance   = vec![1e-4; n_faces];

    println!("\n  Unit box with explicit tolerances (vtx=1e-6, edge=1e-5, face=1e-4):");
    println!("    vertex_tolerance(0) = {:.2e}  (expect 1.00e-6)", vertex_tolerance(&brep, 0));
    println!("    edge_tolerance(3)   = {:.2e}  (expect 1.00e-5)", edge_tolerance(&brep, 3));
    println!("    face_tolerance(5)   = {:.2e}  (expect 1.00e-4)", face_tolerance(&brep, 5));
    println!("    model_tolerance     = {:.2e}  (expect 1.00e-4, max)", model_tolerance(&brep));

    // Zero tolerance falls back to CONFUSION
    brep.geom.vertex_tolerance[0] = 0.0;
    println!("\n  After zeroing vertex[0]:");
    println!("    vertex_tolerance(0) = {:.2e}  (expect {:.2e}, fallback)", vertex_tolerance(&brep, 0), CONFUSION);

    // Out-of-range index → CONFUSION
    println!("\n  Out-of-range edge index 999:");
    println!("    edge_tolerance(999) = {:.2e}  (expect {:.2e}, fallback)", edge_tolerance(&brep, 999), CONFUSION);
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔════════════════════════════════════════════════════╗");
    println!("║              RCAD Phase I Demo                     ║");
    println!("║   2D B-Spline PCurve · Per-Entity Tolerance        ║");
    println!("╚════════════════════════════════════════════════════╝");

    demo_bspline_curve2();
    demo_curve2d_in_geomstore();
    demo_tolerance();

    println!("\n✓ Phase I demo complete.");
}
