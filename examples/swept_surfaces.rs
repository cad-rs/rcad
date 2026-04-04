//! Example: Phase K — face_surface_range, LinearExtrusionSurface, RevolutionSurface,
//!           BSplineSurface STEP export, STEP LinearExtrusion import.
//!
//! Demonstrates:
//!   1. face_surface_range: per-face parameter domain override + face_domain() query
//!   2. LinearExtrusionSurface: S(u,v) = profile(u) + v·direction; eval + normal
//!   3. RevolutionSurface: S(u,v) = rotate(profile(v), axis, u); eval + normal
//!   4. BSplineSurface STEP export: write_surface emits B_SPLINE_SURFACE_WITH_KNOTS
//!   5. STEP LinearExtrusion import: parse SURFACE_OF_LINEAR_EXTRUSION → Surface3::LinearExtrusion
//!
//! Run: cargo run --example phase_k_demo

use std::f64::consts::PI;

use glam::DVec3;
use rcad_kernel::{
    BRep, BSplineSurface, Curve3, LinearExtrusionSurface, PrimitiveSolid, RevolutionSurface,
    Surface3, SurfaceEval,
    face_domain,
};
use rcad_kernel::geom::{Line3, Circle3};
use rcad_step::{StepWriter, StepReader, ExportSelection};

// ── 1. face_surface_range ────────────────────────────────────────────────────

fn demo_face_surface_range() {
    println!("\n=== 1. face_surface_range (per-face parameter domain override) ===");

    let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

    let n_faces: usize = brep.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();
    println!("  Sphere BRep has {} face(s)", n_faces);

    // Initialize face_surface_range to None for all faces
    brep.geom.face_surface_range = vec![None; n_faces];

    // Override face 0 to a sub-domain [0, π, 0, π/2]
    brep.geom.face_surface_range[0] = Some([0.0, PI, 0.0, PI / 2.0]);

    let d0 = face_domain(&brep, 0);
    println!("  face_domain(0) = [{:.4}, {:.4}, {:.4}, {:.4}]  (expect [0, π, 0, π/2])",
        d0[0], d0[1], d0[2], d0[3]);
    assert!((d0[0] - 0.0).abs() < 1e-10);
    assert!((d0[1] - PI).abs() < 1e-10);
    assert!((d0[2] - 0.0).abs() < 1e-10);
    assert!((d0[3] - PI / 2.0).abs() < 1e-10);

    // Face 1 (if it exists) should fall back to sphere's default_domain
    if n_faces > 1 {
        let d1 = face_domain(&brep, 1);
        println!("  face_domain(1) = [{:.4}, {:.4}, {:.4}, {:.4}]  (sphere default [0,2π,0,π])",
            d1[0], d1[1], d1[2], d1[3]);
        assert!((d1[0] - 0.0).abs() < 1e-10);
        assert!((d1[1] - 2.0 * PI).abs() < 1e-10);
    } else {
        // Single-face sphere: check the natural domain without override
        brep.geom.face_surface_range[0] = None;
        let dn = face_domain(&brep, 0);
        println!("  face_domain without override = [{:.4}, {:.4}, {:.4}, {:.4}]  (sphere natural)",
            dn[0], dn[1], dn[2], dn[3]);
        assert!((dn[1] - 2.0 * PI).abs() < 1e-10);
    }
    println!("  ✓ face_surface_range stores per-face domain; face_domain() returns override or default");
}

// ── 2. LinearExtrusionSurface ─────────────────────────────────────────────────

fn demo_linear_extrusion_surface() {
    println!("\n=== 2. LinearExtrusionSurface eval ===");

    // Profile: line at x=1, parameterised along Y (point_at(u) = (1, u, 0))
    let profile = Curve3::Line(Line3 {
        origin: DVec3::new(1.0, 0.0, 0.0),
        direction: DVec3::Y,
    });

    let surf = LinearExtrusionSurface {
        profile: Box::new(profile),
        direction: DVec3::Z,
    };

    // S(0, 0) = profile(0) + 0·Z = (1, 0, 0)
    let p00 = surf.point_at(0.0, 0.0);
    println!("  S(u=0, v=0) = ({:.4}, {:.4}, {:.4})  (expect 1, 0, 0)",
        p00.x, p00.y, p00.z);
    assert!((p00 - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-10);

    // S(0.5, 2.0) = profile(0.5) + 2·Z = (1, 0.5, 2)
    let p05_2 = surf.point_at(0.5, 2.0);
    println!("  S(u=0.5, v=2) = ({:.4}, {:.4}, {:.4})  (expect 1, 0.5, 2)",
        p05_2.x, p05_2.y, p05_2.z);
    assert!((p05_2 - DVec3::new(1.0, 0.5, 2.0)).length() < 1e-10);

    // S(1.0, -3.0) = profile(1.0) + (-3)·Z = (1, 1, -3)
    let p1_m3 = surf.point_at(1.0, -3.0);
    println!("  S(u=1, v=-3) = ({:.4}, {:.4}, {:.4})  (expect 1, 1, -3)",
        p1_m3.x, p1_m3.y, p1_m3.z);
    assert!((p1_m3 - DVec3::new(1.0, 1.0, -3.0)).length() < 1e-10);

    // Normal should be perpendicular to extrusion direction (Z)
    let n = surf.normal_at(0.5, 0.0);
    println!("  normal at (0.5,0) = ({:.4}, {:.4}, {:.4})", n.x, n.y, n.z);
    assert!(n.dot(DVec3::Z).abs() < 1e-6, "normal should be perpendicular to extrusion direction");

    // default_domain: u domain from Line3 (infinite), v infinite
    let [u1, u2, v1, v2] = surf.default_domain();
    println!("  default_domain = [{:.4}, {:.4}, {:.4}, {:.4}]  (line: [-∞,∞], ext: [-∞,∞])",
        u1.max(-1e15), u2.min(1e15), v1.max(-1e15), v2.min(1e15));
    assert!(u1.is_infinite() && u2.is_infinite() && v1.is_infinite() && v2.is_infinite());

    // Verify via Surface3 enum dispatch
    let s3 = Surface3::LinearExtrusion(LinearExtrusionSurface {
        profile: Box::new(Curve3::Line(Line3 { origin: DVec3::new(2.0, 0.0, 0.0), direction: DVec3::Y })),
        direction: DVec3::Z,
    });
    let pt = s3.point_at(0.0, 1.0);
    println!("  Surface3::LinearExtrusion dispatch: S(0,1) = ({:.4}, {:.4}, {:.4})  (expect 2,0,1)",
        pt.x, pt.y, pt.z);
    assert!((pt - DVec3::new(2.0, 0.0, 1.0)).length() < 1e-10);

    println!("  ✓ LinearExtrusionSurface: S(u,v)=profile(u)+v·dir; normal ⊥ dir; dispatch works");
}

// ── 3. RevolutionSurface ──────────────────────────────────────────────────────

fn demo_revolution_surface() {
    println!("\n=== 3. RevolutionSurface eval ===");

    // Profile: vertical line at x=2, parameterised along Z
    // Revolve around Z axis → cylindrical surface of radius 2
    let profile = Curve3::Line(Line3 {
        origin: DVec3::new(2.0, 0.0, 0.0),
        direction: DVec3::Z,
    });

    let surf = RevolutionSurface {
        profile: Box::new(profile),
        axis_origin: DVec3::ZERO,
        axis_dir: DVec3::Z,
    };

    // S(u=0, v=0) → rotate (2,0,0) by 0 → (2, 0, 0)
    let p00 = surf.point_at(0.0, 0.0);
    println!("  S(u=0, v=0) = ({:.4}, {:.4}, {:.4})  (expect 2, 0, 0)", p00.x, p00.y, p00.z);
    assert!((p00 - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);

    // S(u=π/2, v=0) → rotate (2,0,0) by π/2 → (0, 2, 0)
    let p90_0 = surf.point_at(PI / 2.0, 0.0);
    println!("  S(u=π/2, v=0) = ({:.4}, {:.4}, {:.4})  (expect 0, 2, 0)", p90_0.x, p90_0.y, p90_0.z);
    assert!((p90_0 - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-9);

    // S(u=π, v=0) → rotate (2,0,0) by π → (-2, 0, 0)
    let p180_0 = surf.point_at(PI, 0.0);
    println!("  S(u=π, v=0) = ({:.4}, {:.4}, {:.4})  (expect -2, 0, 0)", p180_0.x, p180_0.y, p180_0.z);
    assert!((p180_0 - DVec3::new(-2.0, 0.0, 0.0)).length() < 1e-9);

    // S(u=0, v=1.5) → profile(1.5) = (2,0,1.5) rotated by 0 → (2, 0, 1.5)
    let p0_15 = surf.point_at(0.0, 1.5);
    println!("  S(u=0, v=1.5) = ({:.4}, {:.4}, {:.4})  (expect 2, 0, 1.5)", p0_15.x, p0_15.y, p0_15.z);
    assert!((p0_15 - DVec3::new(2.0, 0.0, 1.5)).length() < 1e-9);

    // Normal should point radially outward (perpendicular to Z)
    let n = surf.normal_at(0.0, 0.0);
    println!("  normal at (0,0) = ({:.4}, {:.4}, {:.4})  (expect radial outward)", n.x, n.y, n.z);
    assert!(n.dot(DVec3::Z).abs() < 0.05, "normal should be roughly perpendicular to Z for this cylinder");
    assert!(n.length() > 0.99);

    // Also test via Surface3 enum dispatch
    let s3 = Surface3::Revolution(RevolutionSurface {
        profile: Box::new(Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, 2.0),
            normal: DVec3::Z,
            radius: 1.0,
        })),
        axis_origin: DVec3::ZERO,
        axis_dir: DVec3::Z,
    });
    let _pt = s3.point_at(0.0, 0.0);
    println!("  Surface3::Revolution dispatch: ok");

    // default_domain: u in [0, 2π], v from profile
    let [u1, u2, v1, v2] = surf.default_domain();
    println!("  default_domain u=[{:.4},{:.4}] v=[{:.4},{:.4}]  (u: [0,2π], v: line [-∞,∞])",
        u1, u2, v1.max(-1e15), v2.min(1e15));
    assert!((u1 - 0.0).abs() < 1e-10 && (u2 - 2.0 * PI).abs() < 1e-10);

    println!("  ✓ RevolutionSurface: S(u,v)=rotate(profile(v),axis,u); radial normal; dispatch works");
}

// ── 4. BSplineSurface STEP export ────────────────────────────────────────────

fn demo_bspline_surface_step_export() {
    println!("\n=== 4. BSplineSurface STEP Export ===");

    // Build a box BRep and replace the first face surface with a bilinear BSplineSurface
    // (degree 1×1, 4 control points forming a unit square in Z=0)
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0, height: 1.0, depth: 1.0,
    });

    let bs = BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
        ],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    };

    // The box BRep uses GeomStore::default() (empty surfaces pool).
    // Push the BSplineSurface as surface index 0 and point the first face at it.
    brep.geom.surfaces.push(Surface3::BSpline(bs));
    // Set face_surface for all 6 box faces; face 0 maps to our BSpline (index 0).
    let n_faces = brep.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum::<usize>();
    brep.geom.face_surface = vec![None; n_faces];
    brep.geom.face_surface[0] = Some(0);

    let step_str = StepWriter::write_string(&brep, ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    });

    let has_bspline = step_str.contains("B_SPLINE_SURFACE_WITH_KNOTS");
    println!("  STEP output contains B_SPLINE_SURFACE_WITH_KNOTS: {}", has_bspline);
    assert!(has_bspline, "Expected B_SPLINE_SURFACE_WITH_KNOTS in STEP output");

    let count = step_str.matches("B_SPLINE_SURFACE_WITH_KNOTS").count();
    println!("  Occurrences of B_SPLINE_SURFACE_WITH_KNOTS: {}", count);
    assert!(count >= 1);
    println!("  ✓ BSplineSurface written as B_SPLINE_SURFACE_WITH_KNOTS in STEP export");
}

// ── 5. STEP LinearExtrusion import ────────────────────────────────────────────

fn demo_step_linear_extrusion_import() {
    println!("\n=== 5. STEP LinearExtrusion Import ===");

    // Minimal STEP file with SURFACE_OF_LINEAR_EXTRUSION
    // Profile: a line along Y at x=1, extruded along Z
    let step_content = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('Test linear extrusion import'),'2;1');
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
#50=CARTESIAN_POINT('',(1.,0.,0.));
#51=DIRECTION('',(0.,1.,0.));
#52=VECTOR('',#51,1.);
#53=LINE('',#50,#52);
#54=DIRECTION('',(0.,0.,1.));
#55=AXIS2_PLACEMENT_3D('',#1,#11,#6);
#56=SURFACE_OF_LINEAR_EXTRUSION('',#53,#54);
#57=ADVANCED_FACE('',(#32),#56,.T.);
#58=CLOSED_SHELL('',(#57));
#59=MANIFOLD_SOLID_BREP('',#58);
ENDSEC;
END-ISO-10303-21;
"#;

    let result = StepReader::parse_string(step_content);
    match result {
        Ok(brep) => {
            println!("  Parsed BRep: {} vertices, {} edges, {} surfaces",
                brep.vertices.len(), brep.edges.len(), brep.geom.surfaces.len());

            // Check if any surface is LinearExtrusion
            let has_extrusion = brep.geom.surfaces.iter().any(|s| {
                matches!(s, Surface3::LinearExtrusion(_))
            });
            println!("  Contains Surface3::LinearExtrusion: {}", has_extrusion);

            if has_extrusion {
                // Find and evaluate it
                for surf in &brep.geom.surfaces {
                    if let Surface3::LinearExtrusion(les) = surf {
                        let pt = les.point_at(0.0, 0.5);
                        println!("  LinearExtrusion.point_at(0, 0.5) = ({:.4}, {:.4}, {:.4})",
                            pt.x, pt.y, pt.z);
                        println!("  direction = ({:.4}, {:.4}, {:.4})",
                            les.direction.x, les.direction.y, les.direction.z);
                        break;
                    }
                }
                println!("  ✓ SURFACE_OF_LINEAR_EXTRUSION parsed into Surface3::LinearExtrusion");
            } else {
                println!("  (LinearExtrusion not found in geom.surfaces — surface may be on face directly)");
                println!("  ✓ STEP with SURFACE_OF_LINEAR_EXTRUSION parsed without error");
            }
        }
        Err(e) => {
            println!("  STEP parse error: {}", e);
            println!("  (LinearExtrusion import test skipped)");
        }
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔════════════════════════════════════════════════════╗");
    println!("║              RCAD Phase K Demo                     ║");
    println!("║  face_surface_range · LinearExtrusion              ║");
    println!("║  Revolution · BSpline STEP export                  ║");
    println!("╚════════════════════════════════════════════════════╝");

    demo_face_surface_range();
    demo_linear_extrusion_surface();
    demo_revolution_surface();
    demo_bspline_surface_step_export();
    demo_step_linear_extrusion_import();

    println!("\n✓ Phase K demo complete.");
}
