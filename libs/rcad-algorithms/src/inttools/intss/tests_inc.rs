#[cfg(test)]
mod tests {
    use crate::inttools::intss::*;
    use glam::DVec3;
    use rcad_kernel::geom::{
        ConicalSurface, Curve2d, Curve2dEval, CylindricalSurface, Plane, SphericalSurface,
        SurfaceEval,
    };

    #[test]
    fn plane_plane_parallel() {
        let p1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let p2 = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
        });
        let r = intersect_surfaces(&p1, &p2);
        assert!(r.is_empty(), "parallel planes: no intersection");
    }

    #[test]
    fn plane_plane_intersect() {
        let p1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let p2 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        });
        let r = intersect_surfaces(&p1, &p2);
        assert_eq!(r.curves.len(), 1);
        assert!(matches!(r.curves[0].curve_3d, SurfaceCurve::Line(_)));
    }

    #[test]
    fn sphere_sphere_equator() {
        // Two equal spheres touching at (1,0,0): each has r=1, centers at (0,0,0) and (2,0,0)
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let r = intersect_surfaces(&s1, &s2);
        assert_eq!(r.curves.len(), 1, "expected one circle");
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            assert!((c.center.x - 0.5).abs() < TOLERANCE_MESH_LEGACY, "center should be at x=0.5");
        } else {
            panic!("expected Circle");
        }
    }

    #[test]
    fn sphere_sphere_disjoint() {
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(5.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let r = intersect_surfaces(&s1, &s2);
        assert!(r.is_empty(), "disjoint spheres: no intersection");
    }

    #[test]
    fn cylinder_cylinder_parallel_intersecting() {
        // Two parallel cylinders r=1 centered at (0,0,0) and (1.5,0,0)
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.5, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        // Two parallel lines
        assert_eq!(r.curves.len(), 2, "expected two intersection lines");
    }

    #[test]
    fn cylinder_cylinder_tangent() {
        // Two parallel cylinders r=1 separated by exactly 2 (tangent externally)
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(2.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        // One tangent line
        assert_eq!(r.curves.len(), 1, "tangent cylinders: one line");
    }

    #[test]
    fn cylinder_cylinder_fuzzy_tolerance_recovers_near_tangent_line() {
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(2.0 + 3.0 * TOLERANCE_ABS, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });

        let strict = intersect_surfaces(&c1, &c2);
        let fuzzy = intersect_surfaces_with_tolerance(&c1, &c2, 4.0 * TOLERANCE_ABS);

        assert!(strict.is_empty(), "strict mode should be disjoint");
        assert_eq!(
            fuzzy.curves.len(),
            1,
            "fuzzy tolerance should recover tangent generator line"
        );
        assert!(matches!(fuzzy.curves[0].curve_3d, SurfaceCurve::Line(_)));
    }

    #[test]
    fn plane_sphere_great_circle() {
        let p = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let r = intersect_surfaces(&p, &s);
        assert_eq!(r.curves.len(), 1);
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            assert!((c.radius - 3.0).abs() < TOLERANCE_MESH_LEGACY);
        } else {
            panic!("expected Circle");
        }
    }

    #[test]
    fn plane_cone_circle_provides_cone_pcurve() {
        let plane_height = 3.0;
        let half_angle = (0.5_f64).atan();
        let expected_slant = plane_height / half_angle.cos();
        let plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, plane_height),
            normal: DVec3::Z,
        });
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: half_angle,
        });

        let r = intersect_surfaces(&plane, &cone);
        assert_eq!(r.curves.len(), 1, "plane-cone circle should give one component");

        let circle = match &r.curves[0].curve_3d {
            SurfaceCurve::Circle(circle) => circle,
            other => panic!("expected Circle, got {other:?}"),
        };
        let pcurve = r.curves[0]
            .pcurve_on_b
            .as_ref()
            .expect("cone-side pcurve should be present");

        match pcurve {
            Curve2d::Line(line) => {
                assert!((line.origin.y - expected_slant).abs() < TOLERANCE_COORD_SUB);
                assert!((line.direction.x - 1.0).abs() < TOLERANCE_COORD_SUB);
                assert!(line.direction.y.abs() < TOLERANCE_COORD_SUB);
            }
            other => panic!("expected analytic cone pcurve line, got {other:?}"),
        }

        for t in [0.0, std::f64::consts::FRAC_PI_2, std::f64::consts::PI] {
            let uv = pcurve.point_at(t);
            let p3 = match &cone {
                Surface3::Cone(surface) => surface.point_at(uv.x, uv.y),
                _ => unreachable!(),
            };
            assert!((p3.z - plane_height).abs() < TOLERANCE_MESH_LEGACY, "lifted point z={} at t={}", p3.z, t);
            assert!(
                (p3.distance(circle.center) - circle.radius).abs() < TOLERANCE_MESH_LEGACY,
                "lifted point radius mismatch at t={}: got {}, expected {}",
                t,
                p3.distance(circle.center),
                circle.radius
            );
        }
    }

    #[test]
    fn cylinder_cylinder_perpendicular_steinmetz() {
        // Two perpendicular cylinders r=1 閳?Steinmetz configuration
        let c1 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(0.0, -2.0, 0.0),
            axis: DVec3::Y,
            ref_dir: any_perpendicular(DVec3::Y),
            radius: 1.0,
        });
        let c2 = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(-2.0, 0.0, 0.0),
            axis: DVec3::X,
            ref_dir: any_perpendicular(DVec3::X),
            radius: 1.0,
        });
        let r = intersect_surfaces(&c1, &c2);
        // Should find the Steinmetz intersection curve(s)
        assert!(!r.curves.is_empty(), "expected at least one intersection curve, got none");
        // The Steinmetz intersection is one or two closed space curves
        if let SurfaceCurve::Polyline(pts) = &r.curves[0].curve_3d {
            assert!(pts.len() >= 4, "polyline should have 閳? points, got {}", pts.len());
        }
    }

    #[test]
    fn torus_perpendicular_plane_gives_circles() {
        use rcad_kernel::geom::ToroidalSurface;

        // Torus with axis=Z, centered at origin, R=5, r=1.
        // Plane at z=0 (perpendicular to the axis) intersects the torus
        // in two concentric circles with radii R+r=6 and R-r=4.
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let r = intersect_surfaces(&torus, &plane);
        assert_eq!(
            r.curves.len(),
            2,
            "torus 閳?perp-plane should give 2 circles, got {}",
            r.curves.len()
        );

        // Collect radii
        let mut radii: Vec<f64> = r
            .curves
            .iter()
            .filter_map(|c| {
                if let SurfaceCurve::Circle(circ) = &c.curve_3d {
                    Some(circ.radius)
                } else {
                    None
                }
            })
            .collect();
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        assert_eq!(radii.len(), 2, "expected 2 Circle3 results");
        assert!(
            (radii[0] - 4.0).abs() < TOLERANCE_MESH_LEGACY,
            "inner circle radius should be 4, got {}",
            radii[0]
        );
        assert!(
            (radii[1] - 6.0).abs() < TOLERANCE_MESH_LEGACY,
            "outer circle radius should be 6, got {}",
            radii[1]
        );
    }

    #[test]
    fn cylinder_cone_coaxial_gives_circle() {
        // Cylinder: r=2, axis Z, origin (0,0,0)
        // Cone: apex (0,0,0), axis Z, half_angle=45鎺?閳?tan=1
        // Coaxial 閳?circle at h = 0 + 2/1 = 2, radius = 2
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 2.0,
        });
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });

        let r = intersect_surfaces(&cyl, &cone);
        assert_eq!(r.curves.len(), 1, "coaxial cylinder-cone should give one circle");
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            assert!((c.center.z - 2.0).abs() < TOLERANCE_MESH_LEGACY, "circle center.z={}", c.center.z);
            assert!((c.radius - 2.0).abs() < TOLERANCE_MESH_LEGACY, "circle radius={}", c.radius);
        } else {
            panic!("expected Circle, got {:?}", r.curves[0].curve_3d);
        }
    }

    #[test]
    fn cone_cone_coaxial_gives_circle() {
        // Cone1: apex (0,0,2), axis Z, 45鎺?(tan=1)
        // Cone2: apex (0,0,0), axis Z, 30鎺?(tan=1/閳?)
        // Coaxial 閳?circle at h = 閳?+1 閳?2.732 from cone1 apex
        let k1 = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(0.0, 0.0, 2.0),
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });
        let k2 = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 30.0_f64.to_radians(),
        });

        let r = intersect_surfaces(&k1, &k2);
        assert_eq!(r.curves.len(), 1, "coaxial cones should give one circle");
        if let SurfaceCurve::Circle(c) = &r.curves[0].curve_3d {
            let expected_r = 3_f64.sqrt() + 1.0;
            assert!(
                (c.radius - expected_r).abs() < TOLERANCE_MESH_LEGACY,
                "circle radius={}, expected {}",
                c.radius,
                expected_r
            );
        } else {
            panic!("expected Circle, got {:?}", r.curves[0].curve_3d);
        }
    }

    #[test]
    fn cone_cone_same_apex_gives_point() {
        // Same apex, different half-angles 閳?CoaxialPoint
        let k1 = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });
        let k2 = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 30.0_f64.to_radians(),
        });

        let r = intersect_surfaces(&k1, &k2);
        assert_eq!(r.curves.len(), 1, "coaxial same-apex cones should give one point");
        assert!(matches!(&r.curves[0].curve_3d, SurfaceCurve::Point(_)));
    }

    #[test]
    fn cone_cone_fuzzy_tolerance_recovers_near_coaxial_circle() {
        // Slightly offset apex in X puts cones outside strict coaxial tolerance.
        // Fuzzy tolerance should recover the coaxial analytic circle branch.
        let k1 = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(0.0, 0.0, 2.0),
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });
        let k2 = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(2.5 * TOLERANCE_ABS, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 30.0_f64.to_radians(),
        });

        let strict = intersect_surfaces(&k1, &k2);
        let fuzzy = intersect_surfaces_with_tolerance(&k1, &k2, 2.0 * TOLERANCE_ABS);

        assert!(!fuzzy.is_empty(), "fuzzy result should not be empty");
        assert!(
            fuzzy
                .curves
                .iter()
                .any(|c| matches!(c.curve_3d, SurfaceCurve::Circle(_) | SurfaceCurve::Point(_))),
            "fuzzy tolerance should recover analytic cone-cone result"
        );

        // In strict mode this near-coaxial case should not be classified as an
        // analytic coaxial intersection.
        assert!(
            !strict
                .curves
                .iter()
                .any(|c| matches!(c.curve_3d, SurfaceCurve::Circle(_) | SurfaceCurve::Point(_))),
            "strict mode unexpectedly produced coaxial analytic result"
        );
    }

    #[test]
    fn sphere_cylinder_fuzzy_tolerance_recovers_near_axis_case() {
        // Sphere center is slightly off-axis: strict mode takes numeric fallback,
        // fuzzy mode should recover analytic circle branch.
        let sph = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(2.0 * TOLERANCE_RETRY_LADDER_MID, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 2.0,
        });

        let strict = intersect_surfaces(&sph, &cyl);
        let fuzzy = intersect_surfaces_with_tolerance(&sph, &cyl, 2.0 * TOLERANCE_RETRY_LADDER_MID);

        assert!(
            fuzzy
                .curves
                .iter()
                .any(|c| matches!(c.curve_3d, SurfaceCurve::Circle(_))),
            "fuzzy mode should recover analytic sphere-cylinder circle"
        );
        assert!(fuzzy.curves.len() >= strict.curves.len());
    }

    #[test]
    fn cylinder_cone_skew_falls_back_to_numeric() {
        // Cylinder: axis Z; Cone: axis X 閳?skew axes 閳?General 閳?numeric
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(0.5, 0.0, 0.0),
            axis: DVec3::X,
            radius: 1.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });

        let r = intersect_surfaces(&cyl, &cone);
        // Should find something via numeric marching
        assert!(!r.is_empty(), "skew cylinder-cone should have numeric intersection");
    }

    #[test]
    fn torus_sphere_on_axis_gives_circles() {
        // Torus: axis=Z, center=origin, R=5, r=2.
        // Sphere: center at origin (on torus axis), radius=5.
        // The torus tube is at (锜?5)铏?+ z铏?= 4; sphere is 锜昏檹 + z铏?= 25.
        // Substituting: 锜昏檹 + 4 - (锜?5)铏?= 25 閳?10锜?- 21 = 25 閳?锜?= 4.6.
        // z铏?= 25 - 4.6铏?= 3.84 閳?z = 鍗?.96 閳?two circles.
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 2.0,
        });
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });

        let r = intersect_surfaces(&torus, &sphere);
        assert_eq!(
            r.curves.len(),
            2,
            "torus 閳?sphere should give 2 circles, got {}",
            r.curves.len()
        );
        for c in &r.curves {
            assert!(matches!(&c.curve_3d, SurfaceCurve::Circle(_)));
        }
    }

    #[test]
    fn torus_sphere_off_axis_skew_polyline() {
        // Torus: axis=Z, R=5, r=2.
        // Sphere: center at (3,0,0), radius=5 (off-axis by 3 units).
        // The intersection should be a closed polyline, not two circles.
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 2.0,
        });
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(3.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });

        let r = intersect_surfaces(&torus, &sphere);
        assert!(!r.curves.is_empty(), "off-axis torus-sphere should have intersection");
        // Should produce at least one polyline (skew solver)
        for c in &r.curves {
            match &c.curve_3d {
                SurfaceCurve::Polyline(pts) => {
                    assert!(pts.len() >= 4, "polyline should have 閳? points");
                }
                other => {
                    panic!("Expected Polyline, got {:?}", other);
                }
            }
        }
    }

    #[test]
    fn torus_cylinder_coaxial_gives_circles() {
        // Torus: axis=Z, R=5, r=1.
        // Cylinder: axis=Z, radius=5 (cuts torus tube at centerline).
        // (5-5)铏?+ h铏?= 1铏?閳?h = 鍗? 閳?two circles.
        // Cylinder: axis=Z, radius=5 (cuts torus tube at centerline).
        // (5-5)铏?+ h铏?= 1铏?閳?h = 鍗? 閳?two circles.
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 5.0,
        });

        let r = intersect_surfaces(&torus, &cyl);
        assert_eq!(
            r.curves.len(),
            2,
            "torus 閳?coaxial cylinder should give 2 circles, got {}",
            r.curves.len()
        );
        for c in &r.curves {
            if let SurfaceCurve::Circle(circ) = &c.curve_3d {
                assert!((circ.radius - 5.0).abs() < TOLERANCE_MESH_LEGACY);
            } else {
                panic!("expected Circle");
            }
        }
    }

    #[test]
    fn torus_cone_coaxial_gives_circle() {
        // Torus: axis=Z, R=5, r=4 (large tube).
        // Cone: apex=origin, axis=Z, 45鎺?(锜?z).
        // Substituting 锜?z into (锜?5)铏?z铏?16:
        //   2z铏?- 10z + 9 = 0 閳?z = (10鍗ら埈?8)/4 閳?{3.82, 1.18} 閳?two circles.
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 4.0,
        });
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(0.0, 0.0, -3.0),
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });

        let r = intersect_surfaces(&torus, &cone);
        assert!(!r.is_empty(), "torus 閳?coaxial cone should have intersection");
        assert!(matches!(&r.curves[0].curve_3d, SurfaceCurve::Circle(_)));
    }

    #[test]
    fn torus_cone_reference_circle_coaxial_still_gives_circles() {
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 4.0,
        });
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(0.0, 0.0, -3.0),
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 45.0_f64.to_radians(),
        });

        let r = intersect_surfaces(&torus, &cone);
        assert_eq!(r.curves.len(), 2, "reference-circle cone should yield the expected two coaxial circles");
        assert!(r.curves.iter().all(|curve| matches!(&curve.curve_3d, SurfaceCurve::Circle(_))));
    }

    #[test]
    fn torus_torus_coaxial_gives_circles() {
        // Torus1: axis=Z, R=5, r=1, center=origin.
        // Torus2: axis=Z, R=5, r=1.5, center=(0,0,0.5).
        // Coaxial, offset 閳?circles where tube circles intersect in (锜?z) plane.
        let t1 = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        let t2 = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::new(0.0, 0.0, 0.5),
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.5,
        });

        let r = intersect_surfaces(&t1, &t2);
        // Should find at least one circle
        assert!(!r.is_empty(), "coaxial tori should have intersection curves");
        for c in &r.curves {
            assert!(matches!(&c.curve_3d, SurfaceCurve::Circle(_)));
        }
    }

    #[test]
    fn torus_skew_cylinder_falls_back_to_numeric() {
        // Torus: axis=Z; Cylinder: axis=X 閳?not coaxial 閳?numeric
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(5.0, 0.0, 0.0),
            axis: DVec3::X,
            ref_dir: any_perpendicular(DVec3::X),
            radius: 0.5,
        });

        let r = intersect_surfaces(&torus, &cyl);
        // Numeric marching should find something
        assert!(!r.is_empty(), "skew torus-cylinder should have numeric intersection");
    }
}

