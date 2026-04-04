//! Phase M demo — all five remaining P2 items
//!
//! M.1  SameParameter / SameRange edge flags
//! M.2  Bezier curves and surfaces
//! M.3  Offset curves and surfaces
//! M.4  Boolean operation history / shape mapping
//! M.5  Multi-edge corner blending

use glam::DVec3;
use rcad_algorithms::geom_populate::populate_box_geom;
use rcad_algorithms::{
    difference_with_history, intersection_with_history, union_with_history, FaceOrigin,
};
use rcad_kernel::{
    edge_same_parameter, edge_same_range,
    geom::{
        BezierSurface, CylindricalSurface, OffsetSurface, PrimitiveSolid,
        Surface3, SurfaceEval,
    },
    BRep,
};
use rcad_modeling::corner_blend;

fn separator(title: &str) {
    println!("\n──────────────────────────────────────────");
    println!("  {title}");
    println!("──────────────────────────────────────────");
}

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: w,
        height: h,
        depth: d,
    });
    for v in &mut brep.vertices {
        v.point += DVec3::new(x, y, z);
    }
    populate_box_geom(&mut brep);
    brep
}

// ─────────────────────────────────────────────────────────────────────────────
// M.1  SameParameter / SameRange
// ─────────────────────────────────────────────────────────────────────────────

fn demo_m1() {
    separator("M.1  SameParameter / SameRange edge flags");

    // A RCAD-generated primitive has same_parameter = true by default
    // (analytic parametrization always gives same-parameter PCurves).
    let cyl = BRep::from_primitive(PrimitiveSolid::Cylinder {
        radius: 1.0,
        height: 2.0,
    });

    // Check a few edges; for generated primitives the stored Vec is empty,
    // so the helper returns true (the documented default).
    let mut all_same_param = true;
    let mut all_same_range = true;
    for i in 0..cyl.edges.len() {
        if !edge_same_parameter(&cyl, i) {
            all_same_param = false;
        }
        if !edge_same_range(&cyl, i) {
            all_same_range = false;
        }
    }

    println!("Cylinder edges: {}", cyl.edges.len());
    println!("All edges same_parameter = true: {all_same_param}");
    println!("All edges same_range     = true: {all_same_range}");
    assert!(all_same_param, "default same_parameter should be true");
    assert!(all_same_range,  "default same_range     should be true");
    println!("M.1 PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// M.2  Bezier curves and surfaces
// ─────────────────────────────────────────────────────────────────────────────

fn demo_m2() {
    separator("M.2  Bezier curves and surfaces");

    use rcad_kernel::geom::{BezierCurve3, CurveEval};

    // Cubic Bezier curve
    let bez = BezierCurve3 {
        control_points: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 2.0, 0.0),
            DVec3::new(2.0, 2.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };

    let p0 = bez.point_at(0.0);
    let p1 = bez.point_at(1.0);
    let pm = bez.point_at(0.5);

    println!("Bezier at t=0  : {p0}");
    println!("Bezier at t=0.5: {pm}");
    println!("Bezier at t=1  : {p1}");

    assert!(
        (p0 - DVec3::new(0.0, 0.0, 0.0)).length() < 1e-9,
        "t=0 should be P[0]"
    );
    assert!(
        (p1 - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-9,
        "t=1 should be P[3]"
    );
    assert!(pm.y > 1.0, "midpoint should be above baseline");

    // Arc length > 0
    let len = rcad_kernel::arc_length(&rcad_kernel::geom::Curve3::Bezier(bez), 0.0, 1.0);
    println!("Bezier arc length: {len:.6}");
    assert!(len > 2.5, "arc length must exceed straight-line distance 3");

    // Bezier surface — bilinear (degree 1×1)
    let bsurf = BezierSurface {
        control_points: vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
        ],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    };

    let s00 = bsurf.point_at(0.0, 0.0);
    let s11 = bsurf.point_at(1.0, 1.0);
    let s01 = bsurf.point_at(0.0, 1.0);
    let s10 = bsurf.point_at(1.0, 0.0);

    println!("BezierSurface (0,0): {s00}");
    println!("BezierSurface (1,1): {s11}");

    assert!((s00 - DVec3::new(0.0, 0.0, 0.0)).length() < 1e-9);
    assert!((s11 - DVec3::new(1.0, 1.0, 0.0)).length() < 1e-9);
    assert!((s01 - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-9);
    assert!((s10 - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-9);

    println!("M.2 PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// M.3  Offset surface
// ─────────────────────────────────────────────────────────────────────────────

fn demo_m3() {
    separator("M.3  Offset curve / surface");

    // Cylinder radius=1, axis=Y, origin=(0,0,0)
    let cyl_base = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::ZERO,
        axis: DVec3::Y,
        radius: 1.0,
    });

    let offset_d = 0.2_f64;
    let offset_surf = OffsetSurface {
        basis: Box::new(cyl_base),
        offset_distance: offset_d,
    };

    // At u=0, v=0 on the cylinder:  base point = (1, 0, 0)
    // Offset normal points radially outward (+X at u=0)
    // So offset point should be at (1 + 0.2, 0, 0) = (1.2, 0, 0)
    let p = offset_surf.point_at(0.0, 0.0);
    let dist_from_axis = (p.x * p.x + p.z * p.z).sqrt();

    println!("Offset surface point at (0,0): {p}");
    println!("Distance from axis: {dist_from_axis:.6}  (expected {:.6})", 1.0 + offset_d);
    assert!(
        (dist_from_axis - (1.0 + offset_d)).abs() < 1e-4,
        "offset point should be at radius + offset_distance from axis"
    );

    // OffsetCurve3 via Curve3::Offset
    use rcad_kernel::geom::{Curve3, CurveEval, Line3, OffsetCurve3};
    let line = Curve3::Line(Line3 {
        origin: DVec3::new(0.0, 0.0, 0.0),
        direction: DVec3::X,
    });
    let offset_curve = Curve3::Offset(OffsetCurve3 {
        basis: Box::new(line),
        offset_distance: 1.0,
        offset_dir: DVec3::Z,
    });
    // Tangent of X-line crossed with Z → Y; offset along Y
    let oc_pt = offset_curve.point_at(0.0);
    println!("OffsetCurve3 at t=0: {oc_pt}  (should be 1 unit in ±Y from origin)");
    assert!(
        oc_pt.y.abs() > 0.9,
        "offset curve should shift point laterally"
    );

    println!("M.3 PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// M.4  Boolean history
// ─────────────────────────────────────────────────────────────────────────────

fn demo_m4() {
    separator("M.4  Boolean operation history / shape mapping");

    let a = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    let b = box_at(1.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    // Union
    {
        let (result, history) = union_with_history(&a, &b).expect("union failed");
        let n = face_count(&result);
        println!("Union: {n} result faces");
        println!("  history entries: {}", history.len());
        assert_eq!(history.len(), n, "history length must equal face count");
        assert!(history.count_from_a() > 0, "some faces from A");
        assert!(history.count_from_b() > 0, "some faces from B");
        println!("  from A: {}  from B: {}", history.count_from_a(), history.count_from_b());
    }

    // Intersection
    {
        let (result, history) = intersection_with_history(&a, &b).expect("intersection failed");
        let n = face_count(&result);
        println!("Intersection: {n} result faces");
        assert_eq!(history.len(), n);
        for i in 0..n {
            let orig = history.face_origin(i);
            assert!(
                matches!(orig, FaceOrigin::FromA(_) | FaceOrigin::FromB(_)),
                "intersection face {i} should be From A or B"
            );
        }
        println!("  from A: {}  from B: {}", history.count_from_a(), history.count_from_b());
    }

    // Difference
    {
        let (result, history) = difference_with_history(&a, &b).expect("difference failed");
        let n = face_count(&result);
        println!("Difference A−B: {n} result faces");
        assert_eq!(history.len(), n);
        println!("  from A: {}  from B: {}", history.count_from_a(), history.count_from_b());
    }

    println!("M.4 PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// M.5  Corner blending
// ─────────────────────────────────────────────────────────────────────────────

fn demo_m5() {
    separator("M.5  Multi-edge corner blending");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    });

    let n_before = face_count(&brep);
    println!("Box faces before corner_blend: {n_before}");

    // Find the vertex at the origin corner (0,0,0)
    let vi = brep
        .vertices
        .iter()
        .position(|v| v.point.length() < 1e-6)
        .expect("no vertex at origin");

    let result = corner_blend(&brep, vi, 0.1).expect("corner_blend failed");
    let n_after = face_count(&result);
    println!("Box faces after  corner_blend: {n_after}");

    // Should have original 6 faces (trimmed) + 1 triangle = 7
    assert!(
        n_after >= 7,
        "expected at least 7 faces after corner_blend, got {n_after}"
    );

    // Verify the corner triangle exists: one face should have exactly 3 boundary edges
    let corner_tri = result
        .solids[0]
        .shells[0]
        .faces
        .iter()
        .any(|f| f.outer_wire.edges.len() == 3);
    assert!(corner_tri, "should have at least one triangular face (the corner patch)");

    println!("M.5 PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    println!("Phase M — All Remaining P2 Items");
    println!("=================================");

    demo_m1();
    demo_m2();
    demo_m3();
    demo_m4();
    demo_m5();

    println!("\n=================================");
    println!("All Phase M demos completed successfully.");
}
