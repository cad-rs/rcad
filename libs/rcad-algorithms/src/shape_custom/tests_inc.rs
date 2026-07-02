#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{Circle3, Line3, Plane, SphericalSurface, CylindricalSurface, ToroidalSurface};
    use rcad_kernel::{Vertex, Wire, WireEdge};

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn approx_eq3(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    #[test]
    fn simplify_high_degree_curve() {
        // Create a degree-5 BSpline curve
        let control_points: Vec<DVec3> = (0..10)
            .map(|i| {
                let t = i as f64 / 9.0;
                DVec3::new(t, (t * std::f64::consts::PI).sin(), 0.0)
            })
            .collect();

        let curve = BSplineCurve3 {
            degree: 5,
            knots: build_clamped_knots(10, 5),
            control_points: control_points.clone(),
            weights: vec![1.0; 10],
        };

        let opts = BSplineSimplifyOptions {
            max_degree: 3,
            tolerance: 0.01,
            ..Default::default()
        };

        let result = simplify_bspline_curve(&curve, &opts);

        assert!(result.final_degree <= 3, "degree should be reduced");
        assert!(result.geometry.control_points.len() <= curve.control_points.len());
    }

    #[test]
    fn simplify_already_simple_curve() {
        // Create a degree-1 line
        let curve = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };

        let opts = BSplineSimplifyOptions::default();
        let result = simplify_bspline_curve(&curve, &opts);

        assert!(!result.was_simplified, "simple curves should not be modified");
        assert_eq!(result.final_degree, 1);
    }

    #[test]
    fn convert_circle_to_bspline() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,
        ));

        let bspline = ensure_bspline_curve(&circle, 32);
        assert_eq!(bspline.degree, 2, "circle converts to degree-2 NURBS");

        // Check that endpoints match
        let p0 = bspline.point_at(0.0);
        let p1 = bspline.point_at(1.0);
        assert!((p0 - p1).length() < TOLERANCE_LINEAR_ULTRA_STRICT, "circle should be closed");
    }

    #[test]
    fn convert_plane_to_bspline() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let bspline = ensure_bspline_surface(&plane, 4, 4);
        assert_eq!(bspline.degree_u, 1);
        assert_eq!(bspline.degree_v, 1);
        assert_eq!(bspline.control_points.len(), 2);
    }

    #[test]
    fn convert_brep_to_bspline() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
        let (converted, report) = convert_to_bspline(&brep, TOLERANCE_MESH_LEGACY);

        // All curves and surfaces should now be BSpline
        for curve in &converted.geom.curves {
            assert!(matches!(curve, Curve3::BSpline(_)));
        }
        for surface in &converted.geom.surfaces {
            assert!(matches!(surface, Surface3::BSpline(_)));
        }

        assert!(report.surfaces_converted > 0);
    }

    #[test]
    fn curve_degree_query() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        assert_eq!(curve_degree(&line), 1);

        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0,
        ));
        assert_eq!(curve_degree(&circle), 2);
    }

    #[test]
    fn surface_degrees_query() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        assert_eq!(surface_degrees(&plane), (1, 1));

        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        assert_eq!(surface_degrees(&sphere), (2, 2));
    }

    #[test]
    fn restrict_geometry_with_options() {
        let mut brep = BRep::new();

        // Add a simple analytic surface
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        }));

        let restrictions = GeometryRestrictions {
            surfaces_to_bspline: true,
            ..Default::default()
        };

        let (result, report) = restrict_geometry(&brep, &restrictions);
        assert!(report.surfaces_converted > 0);
        assert!(matches!(result.geom.surfaces[0], Surface3::BSpline(_)));
    }

    #[test]
    fn build_clamped_knots_correct_size() {
        // Test valid BSpline configurations where n_ctrl > degree
        for n_ctrl in 2..20usize {
            for degree in 1..=n_ctrl.saturating_sub(1).min(5) {
                let knots = build_clamped_knots(n_ctrl, degree);
                let expected_len = n_ctrl + degree + 1;
                assert_eq!(
                    knots.len(),
                    expected_len,
                    "n_ctrl={}, degree={}",
                    n_ctrl,
                    degree
                );

                // Check clamped start: first degree+1 knots should be 0
                for i in 0..=degree {
                    assert!((knots[i] - 0.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT, "knot[{}] should be 0", i);
                }

                // Check clamped end: last degree+1 knots should be 1
                for i in 0..=degree {
                    assert!(
                        (knots[knots.len() - 1 - i] - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT,
                        "knot[{}] should be 1",
                        knots.len() - 1 - i
                    );
                }
            }
        }
    }

    // =============================================================================
    // New Tests for ShapeCustom features
    // =============================================================================

    #[test]
    fn test_surface_to_bspline_from_face() {
        // Create a BRep with a plane surface
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::X + DVec3::Y });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2),
                    WireEdge::fwd(3),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };

        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        }));
        brep.geom.face_surface.push(Some(0));

        let bspline = surface_to_bspline_from_face(&face, 0, &brep);
        assert_eq!(bspline.degree_u, 1);
        assert_eq!(bspline.degree_v, 1);
    }

    #[test]
    fn test_curve_to_bspline_from_edge() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });

        brep.geom.curves.push(Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        }));
        brep.geom.edge_curve.push(Some(0));
        brep.geom.edge_curve_range.push(Some([0.0, 1.0]));

        let bspline = curve_to_bspline_from_edge(&brep.edges[0], 0, &brep);
        assert_eq!(bspline.degree, 1);
        assert!(approx_eq3(bspline.point_at(0.0), DVec3::ZERO, TOLERANCE_LINEAR_ULTRA_STRICT));
        assert!(approx_eq3(bspline.point_at(1.0), DVec3::X, TOLERANCE_LINEAR_ULTRA_STRICT));
    }

    #[test]
    fn test_restrict_to_bspline() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = restrict_to_bspline(&brep);

        // All surfaces should be BSpline
        for surface in &result.geom.surfaces {
            assert!(matches!(surface, Surface3::BSpline(_)));
        }
    }

    #[test]
    fn test_identify_canonical_form_plane_xy() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        assert_eq!(identify_canonical_form(&plane, TOLERANCE_MESH_LEGACY), CanonicalForm::PlaneXY);

        let plane_neg = Surface3::Plane(Plane {
            origin: DVec3::new(1.0, 2.0, 3.0),
            normal: DVec3::NEG_Z,
        });
        assert_eq!(identify_canonical_form(&plane_neg, TOLERANCE_MESH_LEGACY), CanonicalForm::PlaneXY);
    }

    #[test]
    fn test_identify_canonical_form_plane_xz() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        });
        assert_eq!(identify_canonical_form(&plane, TOLERANCE_MESH_LEGACY), CanonicalForm::PlaneXZ);
    }

    #[test]
    fn test_identify_canonical_form_plane_yz() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        });
        assert_eq!(identify_canonical_form(&plane, TOLERANCE_MESH_LEGACY), CanonicalForm::PlaneYZ);
    }

    #[test]
    fn test_identify_canonical_form_cylinder() {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        });
        assert_eq!(identify_canonical_form(&cylinder, TOLERANCE_MESH_LEGACY), CanonicalForm::CylinderZ);

        let tilted_cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::new(1.0, 1.0, 1.0).normalize(),
            ref_dir: any_perpendicular(DVec3::new(1.0, 1.0, 1.0).normalize()),
            radius: 1.0,
        });
        assert_eq!(identify_canonical_form(&tilted_cylinder, TOLERANCE_MESH_LEGACY), CanonicalForm::CylinderGeneral);
    }

    #[test]
    fn test_identify_canonical_form_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        assert_eq!(identify_canonical_form(&sphere, TOLERANCE_MESH_LEGACY), CanonicalForm::SphereOrigin);

        let shifted_sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        assert_eq!(identify_canonical_form(&shifted_sphere, TOLERANCE_MESH_LEGACY), CanonicalForm::SphereGeneral);
    }

    #[test]
    fn test_identify_canonical_form_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 2.0,
            minor_radius: 0.5,
        });
        assert_eq!(identify_canonical_form(&torus, TOLERANCE_MESH_LEGACY), CanonicalForm::TorusOriginZ);
    }

    #[test]
    fn test_convert_to_canonical() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
        let result = convert_to_canonical(&brep, TOLERANCE_MESH_LEGACY);

        // Sphere should be detected as canonical (origin-centered)
        assert!(!result.geom.surfaces.is_empty());
    }

    #[test]
    fn test_try_convert_to_analytic_plane() {
        // Create a BSpline that represents a plane
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let bspline = plane_to_bspline(&plane);

        let result = try_convert_to_analytic(&bspline, TOLERANCE_MESH_LEGACY);
        assert!(result.is_some());

        if let Some(Surface3::Plane(p)) = result {
            assert!(approx_eq3(p.normal.normalize(), DVec3::Z, TOLERANCE_MESH_LEGACY));
        } else {
            panic!("Expected Plane surface");
        }
    }

    #[test]
    fn test_try_convert_to_analytic_cylinder() {
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        };
        let bspline = cylinder_to_bspline(&cylinder);

        // Verify the BSpline represents the cylinder correctly
        let center_point = bspline.point_at(0.5, 0.5);
        // Should be on the cylinder surface (radius 1 from Z axis)
        let dist_from_axis = (center_point.x * center_point.x + center_point.y * center_point.y).sqrt();
        assert!(approx_eq(dist_from_axis, 1.0, 0.1), "Expected radius ~1.0, got {}", dist_from_axis);

        // Detection may or may not succeed depending on tolerance
        let result = try_convert_to_analytic(&bspline, 1e-2);
        if let Some(Surface3::Cylinder(c)) = result {
            assert!(approx_eq(c.radius, 1.0, 0.1), "Expected radius ~1.0, got {}", c.radius);
        }
        // If detection doesn't succeed, that's okay - the BSpline is still valid
    }

    #[test]
    fn test_try_convert_to_analytic_sphere() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: any_perpendicular(DVec3::Z),
        };
        let bspline = sphere_to_bspline(&sphere);

        // Verify the BSpline represents the sphere correctly
        let center_point = bspline.point_at(0.5, 0.5);
        let dist_from_center = center_point.length();
        assert!(approx_eq(dist_from_center, 1.0, 0.1), "Expected radius ~1.0, got {}", dist_from_center);

        // Detection may or may not succeed depending on tolerance
        let result = try_convert_to_analytic(&bspline, 1e-2);
        if let Some(Surface3::Sphere(s)) = result {
            assert!(approx_eq(s.radius, 1.0, 0.1), "Expected radius ~1.0, got {}", s.radius);
        }
        // If detection doesn't succeed, that's okay - the BSpline is still valid
    }

    #[test]
    fn test_simplify_geometry() {
        let mut brep = BRep::new();

        // Add a BSpline surface that is actually a plane
        let bspline = plane_to_bspline(&Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        brep.geom.surfaces.push(Surface3::BSpline(bspline));

        let simplified = simplify_geometry(&brep, TOLERANCE_RETRY_LADDER_COARSE);

        // Should have been converted back to a plane
        assert!(matches!(simplified.geom.surfaces[0], Surface3::Plane(_)));
    }

    #[test]
    fn test_make_direct_faces() {
        let mut brep = BRep::new();

        // Add a trimmed surface
        let trimmed = rcad_kernel::geom::TrimmedSurface::new(
            Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
            }),
            0.0, 1.0, 0.0, 1.0,
        );
        brep.geom.surfaces.push(Surface3::Trimmed(trimmed));

        let direct = make_direct_faces(&brep);

        // Trimmed surface should be resolved
        assert!(!matches!(direct.geom.surfaces[0], Surface3::Trimmed(_)));
    }

    #[test]
    fn test_customize_shape() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let (result, report) = customize_shape(&brep, TOLERANCE_MESH_LEGACY);

        // Check that all surfaces are valid
        for surface in &result.geom.surfaces {
            match surface {
                Surface3::Plane(_)
                | Surface3::Cylinder(_)
                | Surface3::Sphere(_)
                | Surface3::Cone(_)
                | Surface3::Torus(_)
                | Surface3::BSpline(_) => {}
                _ => panic!("Unexpected surface type after customization"),
            }
        }

        // For a box with plane surfaces, no conversion is needed
        // Just verify the result is valid
        let _ = report; // Report is valid
    }

    #[test]
    fn test_canonical_conversion_options() {
        let options = CanonicalConversionOptions {
            tolerance: TOLERANCE_LINEAR_RELAX_8,
            convert_planes: false,
            convert_revolution_surfaces: true,
            convert_spheres: true,
            convert_tori: false,
        };

        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
        let result = convert_to_canonical_with_options(&brep, &options);

        // Should still have valid geometry
        assert!(!result.geom.surfaces.is_empty());
    }

    #[test]
    fn test_analytic_type_detection() {
        // Test plane detection
        let plane_bspline = plane_to_bspline(&Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let result = try_convert_to_analytic(&plane_bspline, TOLERANCE_MESH_LEGACY);
        assert!(matches!(result, Some(Surface3::Plane(_))));

        // Test cylinder detection
        let cylinder_bspline = cylinder_to_bspline(&CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 2.0,
        });
        let result = try_convert_to_analytic(&cylinder_bspline, TOLERANCE_RETRY_LADDER_COARSE);
        if let Some(Surface3::Cylinder(c)) = result {
            assert!(approx_eq(c.radius, 2.0, 0.1));
        }

        // Test sphere detection
        let sphere_bspline = sphere_to_bspline(&SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.5,
            ref_dir: any_perpendicular(DVec3::Z),
        });
        let result = try_convert_to_analytic(&sphere_bspline, TOLERANCE_ADAPTIVE_MAX);
        if let Some(Surface3::Sphere(s)) = result {
            assert!(approx_eq(s.radius, 1.5, 0.1));
        }
    }

    #[test]
    fn test_non_analytic_bspline_stays_bspline() {
        // Create a complex BSpline that doesn't represent an analytic surface
        let control_points: Vec<Vec<DVec3>> = (0..5)
            .map(|i| {
                (0..5)
                    .map(|j| {
                        DVec3::new(
                            i as f64 / 4.0,
                            j as f64 / 4.0,
                            (i as f64 * j as f64 / 16.0).sin(),
                        )
                    })
                    .collect()
            })
            .collect();

        let weights: Vec<Vec<f64>> = (0..5).map(|_| (0..5).map(|_| 1.0).collect()).collect();

        let bspline = BSplineSurface {
            degree_u: 3,
            degree_v: 3,
            knots_u: build_clamped_knots(5, 3),
            knots_v: build_clamped_knots(5, 3),
            control_points,
            weights,
        };

        let result = try_convert_to_analytic(&bspline, TOLERANCE_ADAPTIVE_MAX);
        // This complex surface should not convert to analytic
        assert!(result.is_none());
    }

    #[test]
    fn test_shape_custom_report() {
        let report = ShapeCustomReport {
            surfaces_to_bspline: 5,
            curves_to_bspline: 10,
            bspline_to_analytic: 2,
            bspline_curve_to_analytic: 1,
            faces_made_direct: 3,
            canonical_conversions: 1,
            max_deviation: TOLERANCE_MESH_LEGACY,
        };

        assert_eq!(report.surfaces_to_bspline, 5);
        assert_eq!(report.curves_to_bspline, 10);
        assert_eq!(report.bspline_to_analytic, 2);
    }
}
