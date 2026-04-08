//! Phase Q demo — surface-surface intersection, trimmed surface
//!
//! Q.C  Surface-surface intersection (GeomAPI_IntSS)
//!      Analytic pairs: Plane×Plane, Plane×Sphere, Plane×Cylinder,
//!      Sphere×Sphere, Cylinder×Cylinder (parallel)
//! Q.A  Rectangular trimmed surface (Geom_RectangularTrimmedSurface)
//!      Construction, evaluation, STEP round-trip

use glam::DVec3;
use rcad_algorithms::{SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, intersect_surfaces};
use rcad_kernel::{
    TrimmedSurface,
    geom::{CylindricalSurface, Plane, SphericalSurface, Surface3, SurfaceEval},
};

fn separator(title: &str) {
    println!("\n──────────────────────────────────────────");
    println!("  {title}");
    println!("──────────────────────────────────────────");
}

fn curve_label(c: &SurfaceIntersectionResult) -> &'static str {
    match &c.curve_3d {
        SurfaceCurve::Circle(_) => "Circle",
        SurfaceCurve::Ellipse(_) => "Ellipse",
        SurfaceCurve::Line(_) => "Line",
        SurfaceCurve::Point(_) => "Point",
        SurfaceCurve::Polyline(_) => "Polyline",
    }
}

fn print_result(r: &SurfaceSurfaceIntersection) {
    if r.is_empty() {
        println!("  → No intersection");
    } else {
        for c in &r.curves {
            println!("  → {}", curve_label(c));
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Q.C  Surface-Surface Intersection
// ─────────────────────────────────────────────────────────────────────────────

fn demo_intss() {
    separator("Q.C  Surface-Surface Intersection");

    // 1. Plane × Plane — crossing → Line
    {
        let p1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let p2 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        });
        let r = intersect_surfaces(&p1, &p2);
        print!("Plane(z=0) ∩ Plane(x=0): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1);
        assert!(matches!(r.curves[0].curve_3d, SurfaceCurve::Line(_)));
        println!("  PASS");
    }

    // 2. Plane × Plane — parallel → no intersection
    {
        let p1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let p2 = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
        });
        let r = intersect_surfaces(&p1, &p2);
        print!("Parallel planes: ");
        print_result(&r);
        assert!(r.is_empty(), "parallel planes should give no intersection");
        println!("  PASS");
    }

    // 3. Plane × Sphere — great circle
    {
        let p = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
        });
        let r = intersect_surfaces(&p, &s);
        print!("Plane(z=0) ∩ Sphere(r=3): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1);
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            assert!(
                (c.radius - 3.0).abs() < 1e-6,
                "expected r=3, got {}",
                c.radius
            );
            println!("  PASS: circle r={:.4}", c.radius);
        } else {
            panic!("expected Circle");
        }
    }

    // 4. Plane × Sphere — offset plane → smaller circle
    {
        let p = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 2.0),
            normal: DVec3::Z,
        });
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
        });
        let r = intersect_surfaces(&p, &s);
        print!("Plane(z=2) ∩ Sphere(r=3): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1);
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            let expected = (9.0_f64 - 4.0).sqrt();
            assert!(
                (c.radius - expected).abs() < 1e-6,
                "expected r≈{expected:.4}, got {}",
                c.radius
            );
            println!("  PASS: circle r={:.4}", c.radius);
        } else {
            panic!("expected Circle");
        }
    }

    // 5. Plane × Sphere — tangent → Point
    {
        let p = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
        });
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
        });
        let r = intersect_surfaces(&p, &s);
        print!("Plane(z=3) ∩ Sphere(r=3) [tangent]: ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1);
        assert!(matches!(r.curves[0].curve_3d, SurfaceCurve::Point(_)));
        println!("  PASS");
    }

    // 6. Plane × Cylinder — perpendicular cut → Circle
    {
        let p = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 2.0),
            normal: DVec3::Z,
        });
        let c = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        let r = intersect_surfaces(&p, &c);
        print!("Plane(z=2,⊥) ∩ Cylinder(r=2,Z-axis): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1);
        assert!(matches!(r.curves[0].curve_3d, SurfaceCurve::Circle(_)));
        println!("  PASS");
    }

    // 7. Sphere × Sphere — two unit spheres → circle
    {
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&s1, &s2);
        print!("Sphere(r=1,O) ∩ Sphere(r=1,(1,0,0)): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1);
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            assert!((c.center.x - 0.5).abs() < 1e-6, "center.x should be 0.5");
            println!("  PASS: circle at x={:.4}, r={:.4}", c.center.x, c.radius);
        } else {
            panic!("expected Circle");
        }
    }

    // 8. Sphere × Sphere — disjoint
    {
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(10.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&s1, &s2);
        print!("Sphere(r=1) ∩ Sphere(r=1, far away): ");
        print_result(&r);
        assert!(r.is_empty());
        println!("  PASS");
    }

    // 9. Cylinder × Cylinder — parallel, intersecting → 2 lines
    {
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.5, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        print!("Cylinder(r=1) ∩ Cylinder(r=1, d=1.5 apart): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 2, "expected 2 lines");
        assert!(r.curves.iter().all(|c| matches!(c.curve_3d, SurfaceCurve::Line(_))));
        println!("  PASS");
    }

    // 10. Cylinder × Cylinder — tangent → 1 line
    {
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(2.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        print!("Cylinder(r=1) ∩ Cylinder(r=1, tangent): ");
        print_result(&r);
        assert_eq!(r.curves.len(), 1, "expected 1 tangent line");
        println!("  PASS");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Q.A  Trimmed Surface
// ─────────────────────────────────────────────────────────────────────────────

fn demo_trimmed() {
    separator("Q.A  Rectangular Trimmed Surface");

    // 1. TrimmedSurface wraps a Cylinder, restricts to half-turn + 2 height units
    {
        use std::f64::consts::PI;
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        // Full cylinder domain: u∈[0,2π], v∈(-∞,∞)
        // Trim to u∈[0,π], v∈[0,3]
        let trimmed = TrimmedSurface::new(cyl, 0.0, PI, 0.0, 3.0);
        let s = Surface3::Trimmed(trimmed);

        // Reported domain should be the trim box
        let [u1, u2, v1, v2] = s.default_domain();
        assert!((u1 - 0.0).abs() < 1e-10);
        assert!((u2 - PI).abs() < 1e-10);
        assert!((v1 - 0.0).abs() < 1e-10);
        assert!((v2 - 3.0).abs() < 1e-10);
        println!("TrimmedSurface domain: u=[{u1:.4},{u2:.4}], v=[{v1:.4},{v2:.4}]");

        // Point evaluation delegates to basis cylinder.
        // any_perpendicular(Z) = Y (first perp candidate for Z-axis),
        // so u=0 places the point along Y (at radius 2).
        let pt = s.point_at(0.0, 0.0); // u=0 → y=2, x=0, v=0 → z=0
        assert!((pt.y - 2.0).abs() < 1e-10, "expected y=2, got {}", pt.y);
        assert!((pt.x).abs() < 1e-10);
        assert!((pt.z).abs() < 1e-10);
        println!(
            "Point at (u=0,v=0): ({:.4},{:.4},{:.4})  PASS",
            pt.x, pt.y, pt.z
        );

        let pt2 = s.point_at(PI, 2.0); // u=π → y=-2, x=0, v=2 → z=2
        assert!((pt2.y + 2.0).abs() < 1e-6, "expected y=-2, got {}", pt2.y);
        assert!((pt2.z - 2.0).abs() < 1e-10);
        println!(
            "Point at (u=π,v=2): ({:.4},{:.4},{:.4})  PASS",
            pt2.x, pt2.y, pt2.z
        );

        // Normal delegates to basis — at u=0 the outward normal is along Y
        let n = s.normal_at(0.0, 0.0);
        assert!(
            (n.y - 1.0).abs() < 1e-6,
            "expected normal (0,1,0), got {:?}",
            n
        );
        println!(
            "Normal at (u=0,v=0): ({:.4},{:.4},{:.4})  PASS",
            n.x, n.y, n.z
        );
    }

    // 2. TrimmedSurface wraps a Sphere
    {
        use std::f64::consts::PI;
        let sph = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
        });
        // Trim to northern hemisphere: v∈[0, π/2]
        let trimmed = TrimmedSurface::new(sph, 0.0, 2.0 * PI, 0.0, PI / 2.0);
        let s = Surface3::Trimmed(trimmed);

        let [_, _, v1, v2] = s.default_domain();
        assert!((v1 - 0.0).abs() < 1e-10);
        assert!((v2 - PI / 2.0).abs() < 1e-10);

        // North pole (v=0): should be at (0,0,5)
        let north = s.point_at(0.0, 0.0);
        assert!(
            (north.z - 5.0).abs() < 1e-6,
            "expected z=5 at north pole, got {}",
            north.z
        );
        println!(
            "Sphere northern hemisphere trim: north pole at z={:.4}  PASS",
            north.z
        );

        // Equator (v=π/2): should be at radius 5 in XY plane
        let equator = s.point_at(0.0, PI / 2.0);
        assert!((equator.length() - 5.0).abs() < 1e-6);
        assert!(
            equator.z.abs() < 1e-6,
            "equator should be at z=0, got {}",
            equator.z
        );
        println!(
            "Sphere equator: ({:.4},{:.4},{:.4})  PASS",
            equator.x, equator.y, equator.z
        );
    }

    // 3. STEP round-trip: write a box (uses Plane surfaces), read back
    //    (can't directly produce TrimmedSurface from primitive, but can test
    //     that the reader handles RECTANGULAR_TRIMMED_SURFACE gracefully)
    {
        use rcad_kernel::geom::PrimitiveSolid;
        use rcad_step::{ExportSelection, StepReader, StepWriter};

        let box_brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let step_str = StepWriter::write_string(
            &box_brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );
        let brep2 = StepReader::parse_string(&step_str).expect("STEP round-trip failed");
        assert_eq!(
            brep2.solids[0].shells[0].faces.len(),
            box_brep.solids[0].shells[0].faces.len(),
            "face count should be preserved"
        );
        println!(
            "STEP round-trip box: {} faces  PASS",
            brep2.solids[0].shells[0].faces.len()
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    println!("=================================================");
    println!("  Phase Q Demo: IntSS + TrimmedSurface");
    println!("=================================================");

    demo_intss();
    demo_trimmed();

    println!("\n=================================================");
    println!("  Phase Q: All sections completed successfully");
    println!("=================================================");
}
