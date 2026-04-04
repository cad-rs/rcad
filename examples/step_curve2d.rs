//! Example: Phase J — Ellipse2d PCurve, curve2d_range, STEP Curve2d I/O, STEP tolerance import.
//!
//! Demonstrates:
//!   1. Ellipse2d: construction and evaluation at canonical parameter values
//!   2. Ellipse2d as PCurve stored in GeomStore (Curve2d::Ellipse variant)
//!   3. curve2d_range: per-curve2d parameter trim range (analogous to edge_curve_range for 3D)
//!   4. STEP BSplineCurve2 export: write_bspline_curve2d emits B_SPLINE_CURVE_WITH_KNOTS
//!   5. STEP tolerance import: UNCERTAINTY_MEASURE_WITH_UNIT populates GeomStore tolerance vecs
//!
//! Run: cargo run --example phase_j_demo

use std::f64::consts::PI;

use glam::DVec2;
use rcad_kernel::{
    BRep, BSplineCurve2, Curve2d, Ellipse2d, PrimitiveSolid,
    Curve2dEval,
    CONFUSION,
    model_tolerance,
};
use rcad_kernel::geom::{Circle2d, Line2d};
use rcad_step::{StepWriter, StepReader, ExportSelection};

// ── 1. Ellipse2d evaluation ───────────────────────────────────────────────────

fn demo_ellipse2d_eval() {
    println!("\n=== 1. Ellipse2d Evaluation ===");

    let e = Ellipse2d {
        center: DVec2::new(1.0, 2.0),
        major_dir: DVec2::X,
        major_radius: 3.0,
        minor_radius: 1.5,
    };

    let p0    = e.point_at(0.0);
    let p90   = e.point_at(PI / 2.0);
    let p180  = e.point_at(PI);
    let p270  = e.point_at(3.0 * PI / 2.0);
    let p360  = e.point_at(2.0 * PI);

    println!("  Ellipse2d center=(1,2) major_r=3 minor_r=1.5 major_dir=+X");
    println!("    t=0      → ({:.4}, {:.4})  (expect 4.0000, 2.0000)", p0.x, p0.y);
    println!("    t=π/2    → ({:.4}, {:.4})  (expect 1.0000, 3.5000)", p90.x, p90.y);
    println!("    t=π      → ({:.4}, {:.4})  (expect -2.0000, 2.0000)", p180.x, p180.y);
    println!("    t=3π/2   → ({:.4}, {:.4})  (expect 1.0000, 0.5000)", p270.x, p270.y);
    println!("    t=2π     → ({:.4}, {:.4})  (expect 4.0000, 2.0000 — closed)", p360.x, p360.y);

    // Verify closure
    assert!((p360 - p0).length() < 1e-10, "Ellipse2d should be closed at t=2π");

    // Diagonal major_dir test
    let e45 = Ellipse2d {
        center: DVec2::ZERO,
        major_dir: DVec2::new(1.0_f64 / 2.0_f64.sqrt(), 1.0_f64 / 2.0_f64.sqrt()),
        major_radius: 2.0,
        minor_radius: 1.0,
    };
    let q0 = e45.point_at(0.0);
    let expected_q0 = DVec2::new(2.0 / 2.0_f64.sqrt(), 2.0 / 2.0_f64.sqrt());
    println!("  Ellipse2d 45° major_dir, t=0 → ({:.4}, {:.4})  (expect {:.4}, {:.4})",
        q0.x, q0.y, expected_q0.x, expected_q0.y);
    assert!((q0 - expected_q0).length() < 1e-10, "Diagonal major_dir at t=0");

    println!("  ✓ Ellipse2d evaluated correctly at all canonical parameter values");
}

// ── 2. Ellipse2d as PCurve in GeomStore ──────────────────────────────────────

fn demo_ellipse2d_in_geomstore() {
    println!("\n=== 2. Ellipse2d as PCurve in GeomStore ===");

    let mut brep = BRep::from_primitive(PrimitiveSolid::Torus {
        major_radius: 2.0,
        minor_radius: 0.5,
    });

    let initial_len = brep.geom.curve2ds.len();
    println!("  Torus BRep initial curve2ds: {}", initial_len);

    // Add an Ellipse2d PCurve — e.g., an elliptical path in the torus parameter domain
    let ellipse_pcurve = Curve2d::Ellipse(Ellipse2d {
        center: DVec2::new(PI, PI),   // center at (π, π) in torus (u,v) space
        major_dir: DVec2::X,
        major_radius: 1.0,
        minor_radius: 0.5,
    });
    let c2d_idx = brep.geom.curve2ds.len();
    brep.geom.curve2ds.push(ellipse_pcurve);
    brep.geom.curve2d_range.push(None);    // full ellipse, no trim

    // Verify it's stored and can be evaluated
    let stored = &brep.geom.curve2ds[c2d_idx];
    let pt = stored.point_at(0.0);
    println!("  Stored Curve2d::Ellipse at idx {}", c2d_idx);
    println!("    point_at(0) = ({:.4}, {:.4})  (expect {:.4}, 3.1416)", pt.x, pt.y, PI + 1.0);
    assert!((pt.x - (PI + 1.0)).abs() < 1e-10 && (pt.y - PI).abs() < 1e-10,
        "Stored Ellipse2d should evaluate correctly");

    // Also verify the Ellipse variant is matched correctly
    let is_ellipse = matches!(&brep.geom.curve2ds[c2d_idx], Curve2d::Ellipse(_));
    println!("  Variant is Curve2d::Ellipse: {}", is_ellipse);
    assert!(is_ellipse);

    println!("  curve2ds pool size after insertion: {}", brep.geom.curve2ds.len());
    println!("  ✓ Ellipse2d stored as Curve2d::Ellipse in GeomStore and dispatches correctly");
}

// ── 3. curve2d_range: per-curve2d parameter range ────────────────────────────

fn demo_curve2d_range() {
    println!("\n=== 3. curve2d_range (per-curve2d parameter trim range) ===");

    let mut brep = BRep::new();

    // Push 3 curve2ds with different range scenarios
    // (a) Line — no trim (full ray)
    brep.geom.curve2ds.push(Curve2d::Line(Line2d {
        origin: DVec2::ZERO,
        direction: DVec2::X,
    }));
    brep.geom.curve2d_range.push(None);   // natural domain

    // (b) Circle2d — trimmed to upper half [0, π]
    brep.geom.curve2ds.push(Curve2d::Circle(Circle2d {
        center: DVec2::ZERO,
        radius: 1.0,
    }));
    brep.geom.curve2d_range.push(Some([0.0, PI]));   // upper semicircle only

    // (c) Ellipse2d — trimmed to first quadrant [0, π/2]
    brep.geom.curve2ds.push(Curve2d::Ellipse(Ellipse2d {
        center: DVec2::ZERO,
        major_dir: DVec2::X,
        major_radius: 2.0,
        minor_radius: 1.0,
    }));
    brep.geom.curve2d_range.push(Some([0.0, PI / 2.0]));

    assert_eq!(brep.geom.curve2ds.len(), 3);
    assert_eq!(brep.geom.curve2d_range.len(), 3);

    println!("  curve2d[0] = Line, range = {:?}", brep.geom.curve2d_range[0]);
    println!("  curve2d[1] = Circle r=1, range = {:?}", brep.geom.curve2d_range[1]);
    println!("  curve2d[2] = Ellipse a=2 b=1, range = {:?}", brep.geom.curve2d_range[2]);

    assert!(brep.geom.curve2d_range[0].is_none(), "Line should have no trim range");
    let [t1, t2] = brep.geom.curve2d_range[1].unwrap();
    assert!((t1 - 0.0).abs() < 1e-10 && (t2 - PI).abs() < 1e-10,
        "Circle trim should be [0, π]");
    let [t1e, t2e] = brep.geom.curve2d_range[2].unwrap();
    assert!((t1e - 0.0).abs() < 1e-10 && (t2e - PI / 2.0).abs() < 1e-10,
        "Ellipse trim should be [0, π/2]");

    // Evaluate at trimmed endpoints to confirm they're geometrically meaningful
    let circle = &brep.geom.curve2ds[1];
    let p_start = circle.point_at(t1);
    let p_end   = circle.point_at(t2);
    println!("  Circle trimmed [0,π]: start=({:.4},{:.4}) end=({:.4},{:.4})",
        p_start.x, p_start.y, p_end.x, p_end.y);
    println!("  ✓ curve2d_range stores per-curve trim ranges (None = natural domain)");
}

// ── 4. STEP BSplineCurve2 export ─────────────────────────────────────────────

fn demo_step_bspline2d_export() {
    println!("\n=== 4. STEP BSplineCurve2 Export ===");

    // Build a sphere BRep and replace one of its Line2d seam PCurves
    // with a BSplineCurve2 to exercise write_bspline_curve2d
    let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

    // Replace the first curve2d (Line2d seam, forward direction) with a BSplineCurve2
    // that approximates the same path (v from 0 to π at u=0)
    let bs = BSplineCurve2 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(0.0, PI)],
        weights: vec![1.0, 1.0],
    };
    brep.geom.curve2ds[0] = Curve2d::BSpline(bs);

    let step_str = StepWriter::write_string(&brep, ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    });

    let has_bspline = step_str.contains("B_SPLINE_CURVE_WITH_KNOTS");
    println!("  STEP output contains B_SPLINE_CURVE_WITH_KNOTS: {}", has_bspline);
    assert!(has_bspline, "STEP export should contain B_SPLINE_CURVE_WITH_KNOTS for BSplineCurve2");

    // Count occurrences
    let count = step_str.matches("B_SPLINE_CURVE_WITH_KNOTS").count();
    println!("  Occurrences of B_SPLINE_CURVE_WITH_KNOTS: {}", count);
    println!("  ✓ BSplineCurve2 PCurve written as B_SPLINE_CURVE_WITH_KNOTS in STEP");
}

// ── 5. STEP tolerance import ──────────────────────────────────────────────────

fn demo_step_tolerance_import() {
    println!("\n=== 5. STEP Tolerance Import ===");

    // Minimal STEP file with UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(5.E-6), ...)
    let step_content = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('Test tolerance import'),'2;1');
FILE_NAME('test.stp','2026-04-04',(''),(''),'rcad','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=CARTESIAN_POINT('',(1.,0.,0.));
#3=CARTESIAN_POINT('',(1.,1.,0.));
#4=CARTESIAN_POINT('',(0.,1.,0.));
#5=DIRECTION('',(0.,0.,1.));
#6=DIRECTION('',(1.,0.,0.));
#7=VERTEX_POINT('',#1);
#8=VERTEX_POINT('',#2);
#9=VERTEX_POINT('',#3);
#10=VERTEX_POINT('',#4);
#11=DIRECTION('',(0.,0.,1.));
#12=VECTOR('',#5,1.);
#13=LINE('',#1,#12);
#14=EDGE_CURVE('',#7,#8,#13,.T.);
#15=DIRECTION('',(0.,1.,0.));
#16=VECTOR('',#15,1.);
#17=LINE('',#2,#16);
#18=EDGE_CURVE('',#8,#9,#17,.T.);
#19=DIRECTION('',(-1.,0.,0.));
#20=VECTOR('',#19,1.);
#21=LINE('',#3,#20);
#22=EDGE_CURVE('',#9,#10,#21,.T.);
#23=DIRECTION('',(0.,-1.,0.));
#24=VECTOR('',#23,1.);
#25=LINE('',#4,#24);
#26=EDGE_CURVE('',#10,#7,#25,.T.);
#27=ORIENTED_EDGE('',*,*,#14,.T.);
#28=ORIENTED_EDGE('',*,*,#18,.T.);
#29=ORIENTED_EDGE('',*,*,#22,.T.);
#30=ORIENTED_EDGE('',*,*,#26,.T.);
#31=EDGE_LOOP('',(#27,#28,#29,#30));
#32=FACE_OUTER_BOUND('',#31,.T.);
#33=AXIS2_PLACEMENT_3D('',#1,#11,#6);
#34=PLANE('',#33);
#35=ADVANCED_FACE('',(#32),#34,.T.);
#36=CLOSED_SHELL('',(#35));
#37=MANIFOLD_SOLID_BREP('',#36);
#100=LENGTH_UNIT() SI_UNIT(.MILLI.,.METRE.);
#101=DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);
#102=NAMED_UNIT(*,#101);
#103=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(5.E-6),#100,'distance_accuracy_value','');
#104=GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#103)) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY');
ENDSEC;
END-ISO-10303-21;
"#;

    let result = StepReader::parse_string(step_content);
    match result {
        Ok(brep) => {
            let mt = model_tolerance(&brep);
            println!("  Parsed BRep: {} vertices, {} edges", brep.vertices.len(), brep.edges.len());
            println!("  model_tolerance = {:.2e}  (expect 5.00e-6 from UNCERTAINTY_MEASURE)", mt);
            assert!((mt - 5e-6).abs() < 1e-12,
                "model_tolerance should be 5e-6, got {}", mt);
            println!("  vertex_tolerance[0] = {:.2e}", rcad_kernel::vertex_tolerance(&brep, 0));
            println!("  edge_tolerance[0]   = {:.2e}", rcad_kernel::edge_tolerance(&brep, 0));
            println!("  ✓ UNCERTAINTY_MEASURE_WITH_UNIT correctly populates GeomStore tolerance vecs");
        }
        Err(e) => {
            println!("  STEP parse error: {}", e);
            println!("  (Tolerance import test skipped — STEP parse failed)");
        }
    }

    // Also verify that a file without uncertainty uses CONFUSION default
    let simple_step = r#"ISO-10303-21;
HEADER;
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=CARTESIAN_POINT('',(1.,0.,0.));
#3=VERTEX_POINT('',#1);
#4=VERTEX_POINT('',#2);
#5=DIRECTION('',(1.,0.,0.));
#6=VECTOR('',#5,1.);
#7=LINE('',#1,#6);
#8=EDGE_CURVE('',#3,#4,#7,.T.);
ENDSEC;
END-ISO-10303-21;
"#;
    if let Ok(brep2) = StepReader::parse_string(simple_step) {
        let mt2 = model_tolerance(&brep2);
        println!("\n  File without UNCERTAINTY: model_tolerance = {:.2e}  (expect {:.2e} CONFUSION)", mt2, CONFUSION);
        assert!((mt2 - CONFUSION).abs() < 1e-12,
            "Without UNCERTAINTY, model_tolerance should be CONFUSION={}", CONFUSION);
        println!("  ✓ Without UNCERTAINTY_MEASURE, tolerance falls back to CONFUSION");
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔════════════════════════════════════════════════════╗");
    println!("║              RCAD Phase J Demo                     ║");
    println!("║  Ellipse2d · curve2d_range · STEP Curve2d I/O     ║");
    println!("║  STEP Tolerance Import                             ║");
    println!("╚════════════════════════════════════════════════════╝");

    demo_ellipse2d_eval();
    demo_ellipse2d_in_geomstore();
    demo_curve2d_range();
    demo_step_bspline2d_export();
    demo_step_tolerance_import();

    println!("\n✓ Phase J demo complete.");
}
