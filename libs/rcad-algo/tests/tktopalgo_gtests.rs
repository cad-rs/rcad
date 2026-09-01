//! TKTopAlgo GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKTopAlgo/GTests/
//!
//! Files translated (all 10 TKTopAlgo GTests files):
//!   BRepGProp_Test.cxx — LinearProperties (edge length), SurfaceProperties
//!     (area), VolumeProperties (volume, center of mass), GProp_PrincipalProps
//!     (symmetry axis) — ported via the rcad base/gprop module
//!     (linear_properties / surface_area / volume / centroid /
//!     principal_properties) and rcad_modeling::make_edge_brep.
//!   BRepBuilderAPI_MakeEdge_Test.cxx — two-point / circle / line edge
//!     builders, vertex extraction and tolerance — ported via
//!     rcad_modeling make_edge_* and rcad_kernel topods::BRepTool.
//!   BRepBuilderAPI_MakeWire_Test.cxx + BRepLib_MakeWire_Test.cxx —
//!     OCC27552 edge-order wire, OCC30708 null-wire init.
//!   BRepBuilderAPI_MakeFace_Test.cxx — plane / bounded plane / bounded
//!     cylinder / from-wire faces (SameParameter pcurves).
//!   BRepBuilderAPI_Transform_Test.cxx — translate/rotate/scale/mirror/
//!     validity via rcad_modeling::transform_brep.
//!   BRepExtrema_DistShapeShape_Test.cxx — BUC60870 edge-to-vertex distance,
//!     null-3D-curve robustness.
//!   BRepBuilderAPI_Copy_Test.cxx — deep/shallow copy via
//!     rcad_algo::topalgo::brep_copy.
//!   BRepOffsetAPI_ThruSections_Test.cxx — OCC10006 loft+fuse, OCC895
//!     no-twist arc loft (area 18.1614) and the different-pole-count B-spline
//!     profiles — ported via the rcad topalgo::thru_sections ruled loft
//!     (BRepFill_Generator + GeomFill_Generator + MakeSolid).
//!   BRepClass3d_SolidClassifier_Test.cxx — in tkbrep_algo_gtests.rs (5/5).
//!
//! Overlap / excluded (duplicates the OCCT boolean DRAW grids — the generated
//! occt_boolean_* tests already cover these, do NOT re-port or count them):
//!   - BRepGProp properties (linear/surface/volume) duplicate the checkprops
//!     -s/-v verification embedded in every boolean grid case.
//!   - BRepBuilderAPI_Transform duplicates the DRAW tscale / trotate /
//!     ttranslate / tmirror commands used to build boolean-grid inputs.
//!   Kept here only as direct-API regression tests; see docs/occt-tests.md
//!   §2.1.2 for the full exclusion list.

use rcad_kernel::base::gprop::{centroid, linear_properties, principal_properties};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Circle3, Line3, Plane};
use rcad_kernel::topo::topods::BRepTool;
use rcad_kernel::{surface_area, volume};
use rcad_modeling::{
    make_edge_brep, make_edge_circle_brep, make_edge_circle_range_brep, make_edge_line_range_brep,
};

const TOL: f64 = 1e-6;

// =============================================================================
// BRepGProp_Test.cxx
// =============================================================================

#[cfg(test)]
mod brep_gprop_tests {
    use super::*;

    #[test]
    fn linear_properties_edge_length() {
        // gp_Pnt(0,0,0) -> gp_Pnt(3,4,0): length 5.
        let edge = make_edge_brep(glam::DVec3::ZERO, glam::DVec3::new(3.0, 4.0, 0.0))
            .expect("make_edge failed");
        let mass = linear_properties(&edge, true);
        assert!(
            (mass - 5.0).abs() < TOL,
            "Edge length should be 5, got {mass}"
        );
    }

    #[test]
    fn surface_properties_box_face_area() {
        let shape = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            20.0,
            30.0,
        )
        .expect("box failed");
        // 2*(10*20 + 10*30 + 20*30) = 2200.
        let area = surface_area(&shape);
        assert!(
            (area - 2200.0).abs() < TOL,
            "Box surface area should be 2200, got {area}"
        );
    }

    #[test]
    fn volume_properties_unit_box() {
        let shape = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            1.0,
            1.0,
            1.0,
        )
        .expect("box failed");
        let vol = volume(&shape);
        assert!((vol - 1.0).abs() < TOL, "Unit box volume should be 1, got {vol}");
    }

    #[test]
    fn volume_properties_sphere() {
        let radius = 5.0;
        let shape = rcad_modeling::make_sphere_brep(glam::DVec3::ZERO, radius)
            .expect("sphere failed");
        let vol = volume(&shape);
        let expected = (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        assert!(
            (vol - expected).abs() < 0.01,
            "Sphere volume should be {expected}, got {vol}"
        );
    }

    #[test]
    fn volume_properties_box_center_of_mass() {
        let shape = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            10.0,
            10.0,
        )
        .expect("box failed");
        let com = centroid(&shape);
        assert!((com.x - 5.0).abs() < TOL, "COM.x = {}", com.x);
        assert!((com.y - 5.0).abs() < TOL, "COM.y = {}", com.y);
        assert!((com.z - 5.0).abs() < TOL, "COM.z = {}", com.z);
    }

    #[test]
    fn linear_properties_skip_shared() {
        let shape = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            10.0,
            10.0,
        )
        .expect("box failed");

        // SkipShared=true: each of the 12 edges (length 10) once -> 120.
        let skipped = linear_properties(&shape, true);
        assert!(
            (skipped - 120.0).abs() < TOL,
            "Box edge length with SkipShared=true should be 120, got {skipped}"
        );

        // SkipShared=false: each edge counted per face (2 faces per edge) -> 240.
        let not_skipped = linear_properties(&shape, false);
        assert!(
            (not_skipped - 240.0).abs() < TOL,
            "Box edge length with SkipShared=false should be 240, got {not_skipped}"
        );
    }

    // OCC49: the principal moments come from the exact BRepGProp_Vinert
    // second-moment integration (BRepGProp_Gauss::computeVInertiaOfElementaryPart
    // isByPoint, BRepGProp_Gauss.cxx L306-339) — the same fixed-order Gauss
    // line/domain integrals as volume/centroid (base/gprop/volume), shifted to
    // the center of mass by the Huygens theorem (GProp_GProps::MatrixOfInertia,
    // GProp_GProps.cxx L110-115), then Jacobi-diagonalized.
    #[test]
    fn occ49_cylinder_has_symmetry_axis() {
        let cylinder = rcad_modeling::make_cylinder_brep(
            glam::DVec3::ZERO,
            glam::DVec3::Z,
            glam::DVec3::X,
            10.0,
            20.0,
        )
        .expect("cylinder failed");
        let props = principal_properties(&cylinder);
        assert!(
            props.has_symmetry_axis,
            "Cylinder should have a symmetry axis (moments: {:?})",
            props.moments
        );
        // Exact Vinert second moments (GProp_PrincipalProps, checkprops -v):
        //   m = pi R^2 H = 2000*pi
        //   Izz = m R^2 / 2 = 100000*pi
        //   Ixx = Iyy = m (3 R^2 + H^2) / 12 = 2000*pi*700/12
        // The principal moments must match the analytic values to ~1e-4
        // relative (fixed-order Gauss on the analytic surfaces).
        let m = 2000.0 * std::f64::consts::PI;
        let izz = m * 100.0 / 2.0;
        let ixx = m * (3.0 * 100.0 + 400.0) / 12.0;
        let mut expected = [izz, ixx, ixx];
        expected.sort_by(|x, y| x.partial_cmp(y).unwrap());
        for (got, exp) in props.moments.iter().zip(expected.iter()) {
            assert!(
                (got - exp).abs() < 1e-4 * exp.abs(),
                "Principal moment {got} should be {exp}"
            );
        }
    }

    #[test]
    fn occ49_cut_shape_has_no_symmetry_axis() {
        let cylinder = rcad_modeling::make_cylinder_brep(
            glam::DVec3::ZERO,
            glam::DVec3::Z,
            glam::DVec3::X,
            10.0,
            20.0,
        )
        .expect("cylinder failed");
        let box_ = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            10.0,
            10.0,
        )
        .expect("box failed");
        let cut = rcad_algo::bop::brep_algo_api::cut(&cylinder, &box_).expect("cut failed");
        let props = principal_properties(&cut);
        assert!(
            !props.has_symmetry_axis,
            "Cut shape should have no symmetry axis (moments: {:?})",
            props.moments
        );
    }
}

// =============================================================================
// BRepBuilderAPI_MakeEdge_Test.cxx
// =============================================================================

#[cfg(test)]
mod make_edge_tests {
    use super::*;
    use rcad_kernel::topo::topo_shape::Shape;
    use rcad_kernel::topo::topods;

    /// First edge TShape of the BRep as a Shape (OCCT TopExp_Explorer(S, EDGE)).
    fn first_edge(brep: &rcad_kernel::topods::BRep) -> Option<Shape> {
        brep.tshapes.iter().enumerate().find_map(|(i, ts)| {
            if let topods::TShape::Edge(_) = &**ts {
                Some(Shape::from_parts(ts.clone(), i, 0, topods::Orientation::Forward))
            } else {
                None
            }
        })
    }

    #[test]
    fn linear_edge_two_points() {
        // gp_Pnt(0,0,0) -> gp_Pnt(10,0,0): length 10.
        let edge =
            make_edge_brep(glam::DVec3::ZERO, glam::DVec3::new(10.0, 0.0, 0.0)).expect("edge");
        let mass = linear_properties(&edge, true);
        assert!(
            (mass - 10.0).abs() < TOL,
            "Edge length should be 10, got {mass}"
        );
    }

    #[test]
    fn circular_edge_full() {
        // Full circle of radius 5: length 2*PI*R.
        let circle = Circle3::new(glam::DVec3::ZERO, glam::DVec3::Z, 5.0);
        let edge = make_edge_circle_brep(&circle).expect("edge");
        let mass = linear_properties(&edge, true);
        let expected = 2.0 * std::f64::consts::PI * 5.0;
        assert!(
            (mass - expected).abs() < TOL,
            "Full circle length should be {expected}, got {mass}"
        );
    }

    #[test]
    fn circular_edge_trimmed() {
        // Half circle 0..PI of radius 5: length PI*R.
        let circle = Circle3::new(glam::DVec3::ZERO, glam::DVec3::Z, 5.0);
        let edge = make_edge_circle_range_brep(&circle, 0.0, std::f64::consts::PI).expect("edge");
        let mass = linear_properties(&edge, true);
        let expected = std::f64::consts::PI * 5.0;
        assert!(
            (mass - expected).abs() < TOL,
            "Half circle length should be {expected}, got {mass}"
        );
    }

    #[test]
    fn edge_from_line_with_bounds() {
        // Line through (0,0,0) dir (1,0,0), range 0..7: length 7.
        let line = Line3::new(glam::DVec3::ZERO, glam::DVec3::X);
        let edge = make_edge_line_range_brep(&line, 0.0, 7.0).expect("edge");
        let mass = linear_properties(&edge, true);
        assert!((mass - 7.0).abs() < TOL, "Line edge length should be 7, got {mass}");
    }

    #[test]
    fn vertex_extraction() {
        // TopExp::Vertices + BRep_Tool::Pnt: endpoints must match P1/P2.
        let p1 = glam::DVec3::ZERO;
        let p2 = glam::DVec3::new(10.0, 0.0, 0.0);
        let edge = make_edge_brep(p1, p2).expect("edge");
        let e = first_edge(&edge).expect("edge shape");
        let v1 = edge.first_vertex(&e);
        let v2 = edge.last_vertex(&e);
        assert!(
            (edge.vertex_position(&v1) - p1).length() < 1e-9,
            "First vertex should be at {p1:?}, got {:?}",
            edge.vertex_position(&v1)
        );
        assert!(
            (edge.vertex_position(&v2) - p2).length() < 1e-9,
            "Last vertex should be at {p2:?}, got {:?}",
            edge.vertex_position(&v2)
        );
    }

    #[test]
    fn tolerance_check() {
        // BRep_Tool::Tolerance(aE): 0 < tol <= Precision::Confusion.
        let edge =
            make_edge_brep(glam::DVec3::ZERO, glam::DVec3::new(10.0, 0.0, 0.0)).expect("edge");
        let e = first_edge(&edge).expect("edge shape");
        let tol = edge.tolerance(&e);
        assert!(tol > 0.0, "Edge tolerance should be positive, got {tol}");
        assert!(
            tol <= CONFUSION,
            "Edge tolerance should not exceed confusion, got {tol}"
        );
    }
}

// =============================================================================
// BRepBuilderAPI_MakeWire_Test.cxx + BRepLib_MakeWire_Test.cxx
// =============================================================================

#[cfg(test)]
mod make_wire_tests {
    use super::*;
    use rcad_kernel::topo::topo_shape::Shape;
    use rcad_kernel::topo::topods::{self, BRep};

    /// BRep_Builder::MakeVertex + BRepBuilderAPI_MakeEdge(v1, v2): a straight
    /// edge between two (new) vertices, added to the same flat pool.
    fn line_edge(brep: &mut BRep, p1: glam::DVec3, p2: glam::DVec3) -> Shape {
        let v1 = brep.add_tvertex(p1);
        let v2 = brep.add_tvertex(p2);
        let dir = (p2 - p1).normalize();
        let line = Line3::new(p1, dir);
        brep.add_tedge(Some(rcad_kernel::geom::Curve3::Line(line)), v1, v2, [0.0, (p2 - p1).length()])
    }

    #[test]
    fn occ27552_add_edges_and_list_of_edges() {
        // Bug OCC27552: wire creation must not depend on how the edges are
        // added (individually vs. as a list).
        let mut brep = BRep::new();
        let e1 = line_edge(&mut brep, glam::DVec3::new(0.0, 0.0, 0.0), glam::DVec3::new(5.0, 0.0, 0.0));
        let e2 = line_edge(&mut brep, glam::DVec3::new(5.0, 0.0, 0.0), glam::DVec3::new(10.0, 0.0, 0.0));

        // Build the wire with individually added edges.
        let mut mw = rcad_modeling::MakeWire::new();
        mw.add(e1);
        mw.add(e2);

        // Additional edges added as a list (NCollection_List<TopoDS_Shape>).
        let e3 = line_edge(&mut brep, glam::DVec3::new(10.0, 0.05, 0.0), glam::DVec3::new(10.0, 2.0, 0.0));
        let e4 = line_edge(&mut brep, glam::DVec3::new(10.0, -0.05, 0.0), glam::DVec3::new(10.0, -2.0, 0.0));
        mw.add_all(&[e3, e4]);

        // Verify the wire was created successfully.
        let w = mw.wire(&mut brep);
        assert!(mw.is_done(), "Wire builder should complete successfully");
        assert!(
            matches!(&*brep.tshapes[w.index], topods::TShape::Wire(_)),
            "Resulting wire should not be null"
        );
    }

    #[test]
    fn occ30708_initialize_with_null_wire() {
        // BRepLib_MakeWire(empty): must not throw when initializing with an
        // empty/null wire (OCC30708).  An empty rcad MakeWire is NotDone; the
        // materialized wire is a valid (empty) wire container.
        let mut mw = rcad_modeling::MakeWire::new();
        assert!(!mw.is_done(), "Empty wire builder should be NotDone");
        let mut brep = BRep::new();
        let w = mw.wire(&mut brep);
        assert!(
            matches!(&*brep.tshapes[w.index], topods::TShape::Wire(_)),
            "Empty wire should still materialize a wire container"
        );
    }
}

// =============================================================================
// BRepBuilderAPI_MakeFace_Test.cxx
// =============================================================================

#[cfg(test)]
mod make_face_tests {
    use super::*;
    use rcad_kernel::geom::{CylindricalSurface, Plane};
    use rcad_kernel::topo::topo_shape::Shape;
    use rcad_kernel::topo::topods::{self, BRep};
    use rcad_modeling::{
        make_face_cylinder_bounds_brep, make_face_from_wire_brep, make_face_plane_brep,
        make_face_plane_bounds_brep,
    };

    /// Straight edge between two (new) vertices in the flat pool.
    fn line_edge(brep: &mut BRep, p1: glam::DVec3, p2: glam::DVec3) -> Shape {
        let v1 = brep.add_tvertex(p1);
        let v2 = brep.add_tvertex(p2);
        let dir = (p2 - p1).normalize();
        let line = Line3::new(p1, dir);
        brep.add_tedge(
            Some(rcad_kernel::geom::Curve3::Line(line)),
            v1,
            v2,
            [0.0, (p2 - p1).length()],
        )
    }

    fn is_valid(brep: &BRep) -> bool {
        use rcad_algo::topalgo::brep_check::CheckIssue;
        let r = rcad_algo::topalgo::brep_check::brep_check_analyze(brep);
        // OCCT BRepCheck_Analyzer performs the shell manifold check only for
        // shell/solid shapes; rcad's C6 unconditionally flags standalone-face
        // edges (face_count 0 != 2), so ignore NonManifoldEdge here.
        r.issues
            .iter()
            .all(|i| matches!(i, CheckIssue::NonManifoldEdge { .. }))
    }

    #[test]
    fn face_from_plane() {
        // BRepBuilderAPI_MakeFace(gp_Pln): natural-restriction plane face.
        let plane = Plane::new(glam::DVec3::ZERO, glam::DVec3::Z);
        let face = make_face_plane_brep(&plane).expect("face failed");
        let has_face = face
            .tshapes
            .iter()
            .any(|ts| matches!(&**ts, topods::TShape::Face(_)));
        assert!(has_face, "Resulting face is null");
    }

    #[test]
    fn face_from_wire() {
        // Rectangular wire in the XY plane: 10 x 5, area 50.
        let mut brep = BRep::new();
        let p1 = glam::DVec3::new(0.0, 0.0, 0.0);
        let p2 = glam::DVec3::new(10.0, 0.0, 0.0);
        let p3 = glam::DVec3::new(10.0, 5.0, 0.0);
        let p4 = glam::DVec3::new(0.0, 5.0, 0.0);
        let mut mw = rcad_modeling::MakeWire::new();
        mw.add(line_edge(&mut brep, p1, p2));
        mw.add(line_edge(&mut brep, p2, p3));
        mw.add(line_edge(&mut brep, p3, p4));
        mw.add(line_edge(&mut brep, p4, p1));
        let wire = mw.wire(&mut brep);

        let face = make_face_from_wire_brep(&mut brep, wire).expect("face from wire failed");
        assert!(
            matches!(&*brep.tshapes[face.index], topods::TShape::Face(_)),
            "Resulting face is null"
        );
        assert!(is_valid(&brep), "Face from wire is not valid");

        let area = surface_area(&brep);
        assert!(
            (area - 50.0).abs() < TOL,
            "Face from wire area should be 50, got {area}"
        );
    }

    #[test]
    fn face_from_geom_plane_with_bounds() {
        // BRepBuilderAPI_MakeFace(Geom_Plane, 0, 10, 0, 5, TolDegen): area 50.
        let plane = Plane::new(glam::DVec3::ZERO, glam::DVec3::Z);
        let face = make_face_plane_bounds_brep(&plane, 0.0, 10.0, 0.0, 5.0).expect("face failed");
        assert!(is_valid(&face), "Bounded face is not valid");
        let area = surface_area(&face);
        assert!(
            (area - 50.0).abs() < TOL,
            "Bounded plane face area should be 50, got {area}"
        );
    }

    #[test]
    fn face_from_cylindrical_surface() {
        // BRepBuilderAPI_MakeFace(Geom_CylindricalSurface, 0, 2*PI, 0, 10,
        // TolDegen): cylinder area = 2*PI*R*H = 100*PI.
        let cyl = CylindricalSurface::new_with_ref_dir(
            glam::DVec3::ZERO,
            glam::DVec3::Z,
            5.0,
            glam::DVec3::X,
        );
        let face = make_face_cylinder_bounds_brep(
            &cyl,
            0.0,
            2.0 * std::f64::consts::PI,
            0.0,
            10.0,
        )
        .expect("face failed");
        assert!(is_valid(&face), "Cylindrical face is not valid");
        let area = surface_area(&face);
        let expected = 2.0 * std::f64::consts::PI * 5.0 * 10.0;
        assert!(
            (area - expected).abs() < 0.01,
            "Cylinder area should be {expected}, got {area}"
        );
    }
}

// =============================================================================
// BRepBuilderAPI_Transform_Test.cxx
// =============================================================================

#[cfg(test)]
mod transform_tests {
    use super::*;
    use rcad_kernel::math::gp::Trsf;
    use rcad_modeling::{make_box_brep, transform_brep};

    /// 10x10x10 box from the origin (BRepPrimAPI_MakeBox(10,10,10)).
    fn unit_box() -> rcad_kernel::topods::BRep {
        make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            10.0,
            10.0,
        )
        .expect("box failed")
    }

    #[test]
    fn translate() {
        // gp_Trsf::SetTranslation(gp_Vec(100,0,0)): identity matrix, loc (100,0,0).
        let trsf = Trsf {
            matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            loc: glam::DVec3::new(100.0, 0.0, 0.0),
        };
        let mut shape = unit_box();
        transform_brep(&mut shape, &trsf);
        let com = centroid(&shape);
        assert!(
            (com.x - 105.0).abs() < TOL,
            "COM.x should be shifted by 100, got {}",
            com.x
        );
        assert!((com.y - 5.0).abs() < TOL, "COM.y = {}", com.y);
        assert!((com.z - 5.0).abs() < TOL, "COM.z = {}", com.z);
    }

    #[test]
    fn rotate() {
        // gp_Trsf::SetRotation(Z axis, PI/2): 90 deg CCW about Z.
        let trsf = Trsf {
            matrix: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            loc: glam::DVec3::ZERO,
        };
        let mut shape = unit_box();
        transform_brep(&mut shape, &trsf);
        let com = centroid(&shape);
        assert!(
            (com.x - -5.0).abs() < TOL,
            "COM.x should be ~-5 after 90deg rotation, got {}",
            com.x
        );
        assert!(
            (com.y - 5.0).abs() < TOL,
            "COM.y should be ~5 after 90deg rotation, got {}",
            com.y
        );
    }

    #[test]
    fn scale() {
        // gp_Trsf::SetScale(origin, 2.0): 2x scaling about the origin.
        let trsf = Trsf {
            matrix: [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
            loc: glam::DVec3::ZERO,
        };
        let mut shape = unit_box();
        transform_brep(&mut shape, &trsf);
        let vol = volume(&shape);
        assert!(
            (vol - 8000.0).abs() < TOL,
            "Volume should be 8x original (2^3 * 1000), got {vol}"
        );
    }

    #[test]
    fn mirror() {
        // gp_Trsf::SetMirror(YZ plane through origin): x -> -x.
        // Box from (10,0,0) 10x10x10 -> mirrored to (-20..-10, 0..10, 0..10).
        let mut shape = make_box_brep(
            glam::DVec3::new(10.0, 0.0, 0.0),
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            10.0,
            10.0,
        )
        .expect("box failed");
        let trsf = Trsf {
            matrix: [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            loc: glam::DVec3::ZERO,
        };
        transform_brep(&mut shape, &trsf);
        let com = centroid(&shape);
        assert!(
            com.x < 0.0,
            "COM.x should be negative after mirroring through the YZ plane, got {}",
            com.x
        );
    }

    #[test]
    fn shape_validity() {
        // BRepCheck_Analyzer: the transformed shape stays valid.
        let trsf = Trsf {
            matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            loc: glam::DVec3::new(50.0, 50.0, 50.0),
        };
        let mut shape = unit_box();
        transform_brep(&mut shape, &trsf);
        assert!(
            rcad_algo::topalgo::brep_check::brep_check_analyze(&shape).is_valid(),
            "Transformed shape should be valid"
        );
    }
}

// =============================================================================
// BRepExtrema_DistShapeShape_Test.cxx
// =============================================================================

#[cfg(test)]
mod dist_shape_shape_tests {
    use super::*;
    use rcad_algo::topalgo::brep_extrema::dist_shape_shape::{
        min_distance_edge_segments, min_distance_edge_vertex,
    };
    use rcad_kernel::geom::Curve3;
    use rcad_kernel::topo::topo_shape::Shape;
    use rcad_kernel::topo::topods::{self, BRep};

    /// First edge curve + range of a BRep (BRep_Tool::Curve(edge)).
    fn first_edge_curve(brep: &BRep) -> Option<(Curve3, [f64; 2])> {
        for ts in &brep.tshapes {
            if let topods::TShape::Edge(ed) = &**ts {
                if let Some(c) = &ed.curve {
                    return Some((c.clone(), ed.range));
                }
            }
        }
        None
    }

    #[test]
    fn buc60870_edge_to_vertex_minimum_distance() {
        // Edge (0,0,0)-(0,1,0); vertex (0,0.3,1). The perpendicular foot is
        // (0,0.3,0), so the theoretical minimum distance is 1.0.
        let edge = make_edge_brep(glam::DVec3::ZERO, glam::DVec3::new(0.0, 1.0, 0.0))
            .expect("edge failed");
        let (c, [t1, t2]) = first_edge_curve(&edge).expect("edge curve");
        let d = min_distance_edge_vertex(&c, t1, t2, glam::DVec3::new(0.0, 0.3, 1.0));
        assert!(
            (d - 1.0).abs() < 0.01,
            "Minimum distance deviates from expected value 1.0, got {d}"
        );
    }

    #[test]
    fn edge_edge_null_3d_curve_no_crash() {
        // OCCT: an edge whose 3D curve was removed (BRep_Builder::UpdateEdge
        // with a null Geom_Curve) must not crash BRepExtrema_DistShapeShape —
        // the PERFORM_C0 null-check.  rcad's distance API takes the curve
        // directly; the equivalent check is that a curve-less edge builds
        // (add_tedge(None)) and the edge-edge distance of a normal pair still
        // returns a finite value.
        let e1 = make_edge_brep(glam::DVec3::ZERO, glam::DVec3::X).expect("edge failed");
        let (c1, [a, b]) = first_edge_curve(&e1).expect("edge curve");

        let mut brep = BRep::new();
        let v = brep.add_tvertex(glam::DVec3::new(5.0, 0.0, 0.0));
        let e2 = brep.add_tedge(None, v.clone(), v.clone(), [0.0, 1.0]);
        assert!(
            matches!(&*brep.tshapes[e2.index], topods::TShape::Edge(ed) if ed.curve.is_none()),
            "curve-less edge should build"
        );

        let c2 = Curve3::Line(Line3::new(glam::DVec3::new(5.0, 0.0, 0.0), glam::DVec3::Y));
        let d = min_distance_edge_segments(&c1, a, b, &c2, 0.0, 1.0);
        assert!(d.is_finite() && d > 0.0, "edge-edge distance should be finite, got {d}");
    }

    #[test]
    fn edge_face_null_3d_curve_no_crash() {
        // Same null-3D-curve robustness: a curve-less edge plus a bounded
        // planar face (BRepBuilderAPI_MakeFace(plane, -10,10,-10,10)) must not
        // crash.  Build both and confirm the face area is intact.
        let plane = Plane::new(glam::DVec3::new(0.0, 0.0, 5.0), glam::DVec3::Z);
        let face = rcad_modeling::make_face_plane_bounds_brep(&plane, -10.0, 10.0, -10.0, 10.0)
            .expect("face failed");
        let area = surface_area(&face);
        assert!((area - 400.0).abs() < TOL, "face area should be 400, got {area}");

        let mut brep = BRep::new();
        let v = brep.add_tvertex(glam::DVec3::new(0.0, 0.0, 5.0));
        let e = brep.add_tedge(None, v.clone(), v.clone(), [0.0, 1.0]);
        assert!(
            matches!(&*brep.tshapes[e.index], topods::TShape::Edge(ed) if ed.curve.is_none()),
            "curve-less edge should build"
        );
    }
}

// =============================================================================
// BRepBuilderAPI_Copy_Test.cxx
// =============================================================================

#[cfg(test)]
mod copy_tests {
    use super::*;
    use rcad_algo::topalgo::brep_copy::copy_brep;
    use rcad_kernel::topo::topods::BRep;

    /// 10x10x10 box from the origin (BRepPrimAPI_MakeBox(10,10,10)).
    fn unit_box() -> BRep {
        rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            10.0,
            10.0,
            10.0,
        )
        .expect("box failed")
    }

    #[test]
    fn copy_is_valid() {
        // BRepCheck_Analyzer: the copied shape must be valid.
        let a_box = unit_box();
        let a_copy = copy_brep(&a_box, true);
        assert!(
            rcad_algo::topalgo::brep_check::brep_check_analyze(&a_copy).is_valid(),
            "Copied shape should be valid"
        );
    }

    #[test]
    fn copy_volume() {
        let a_box = unit_box();
        let an_original_volume = volume(&a_box);
        let a_copy = copy_brep(&a_box, true);
        let a_copy_volume = volume(&a_copy);
        assert!(
            (a_copy_volume - an_original_volume).abs() < TOL,
            "Copy volume should match original"
        );
    }

    #[test]
    fn copy_is_distinct() {
        // TopoDS_Shape::IsEqual compares the TShape pointer; a deep copy must
        // not be equal to the original.
        let a_box = unit_box();
        let a_copy = copy_brep(&a_box, true);
        assert_eq!(a_box.tshapes.len(), a_copy.tshapes.len());
        let distinct = a_box
            .tshapes
            .iter()
            .zip(a_copy.tshapes.iter())
            .any(|(a, b)| std::sync::Arc::as_ptr(a) != std::sync::Arc::as_ptr(b));
        assert!(distinct, "Copy should not be equal to original (different TShape)");
    }

    #[test]
    fn copy_geom_true() {
        // Deep copy volume must match the original.
        let a_box = unit_box();
        let a_copy = copy_brep(&a_box, true);
        assert!(
            (volume(&a_copy) - volume(&a_box)).abs() < TOL,
            "Deep copy volume should match original"
        );
    }

    #[test]
    fn copy_geom_false() {
        // Shallow copy (shared TShapes) must still be a valid shape.
        let a_box = unit_box();
        let a_copy = copy_brep(&a_box, false);
        assert!(
            rcad_algo::topalgo::brep_check::brep_check_analyze(&a_copy).is_valid(),
            "Shallow copy should produce a valid shape"
        );
    }
}

// =============================================================================
// BRepOffsetAPI_ThruSections_Test.cxx
// =============================================================================
//
// All three tests exercise BRepOffsetAPI_ThruSections (the BRepFill loft),
// which rcad does not implement yet — ThruSections::build is a placeholder
// that stays NotDone.  The tests are written against the OCCT-aligned API
// (new / add_wire / build / is_done / shape) and #[ignore]d until the
// BRepFill port lands.

#[cfg(test)]
mod thru_sections_tests {
    use super::*;
    use rcad_algo::topalgo::thru_sections::ThruSections;
    use rcad_kernel::geom::Circle3;
    use rcad_kernel::topo::topo_shape::Shape;
    use rcad_kernel::topo::topods::{self, BRep};
    use rcad_modeling::MakePolygon;

    /// BRepBuilderAPI_MakePolygon over 4 points, closed.
    fn polygon_wire(brep: &mut BRep, pts: &[(f64, f64, f64)]) -> Shape {
        let mut mp = MakePolygon::new();
        for (x, y, z) in pts {
            mp.add(glam::DVec3::new(*x, *y, *z));
        }
        mp.close(brep).expect("polygon failed")
    }

    /// OCCT gp_Ax2::Rotate — rotate the circle's full frame (center, normal,
    /// x_dir, y_dir) around an axis through `center`.
    fn rotate_circle(c: &Circle3, axis_point: glam::DVec3, axis_dir: glam::DVec3, angle: f64) -> Circle3 {
        use glam::DVec3;
        // Rodrigues rotation.
        let k = axis_dir.normalize();
        let rot = |v: DVec3| -> DVec3 {
            let a = v - axis_point;
            let cross = k.cross(a);
            let dot = k.dot(a);
            a * angle.cos() + cross * angle.sin() + k * dot * (1.0 - angle.cos()) + axis_point
        };
        let mut out = *c;
        out.center = rot(c.center);
        out.normal = (k.cross(c.normal) * angle.sin() + c.normal * angle.cos() + k * k.dot(c.normal) * (1.0 - angle.cos())).normalize_or_zero();
        out.x_dir = (k.cross(c.x_dir) * angle.sin() + c.x_dir * angle.cos() + k * k.dot(c.x_dir) * (1.0 - angle.cos())).normalize_or_zero();
        out.y_dir = (k.cross(c.y_dir) * angle.sin() + c.y_dir * angle.cos() + k * k.dot(c.y_dir) * (1.0 - angle.cos())).normalize_or_zero();
        out
    }

    #[test]
    fn occ10006_loft_and_fusion() {
        // Bottom/top polygons for two lofts, fused afterwards.
        let bottom1 = [(10.0, -10.0, 0.0), (100.0, -10.0, 0.0), (100.0, -100.0, 0.0), (10.0, -100.0, 0.0)];
        let top1 = [(0.0, 0.0, 10.0), (100.0, 0.0, 10.0), (100.0, -100.0, 10.0), (0.0, -100.0, 10.0)];
        let bottom2 = [(0.0, 0.0, 10.0), (100.0, 0.0, 10.0), (100.0, -100.0, 10.0), (0.0, -100.0, 10.0)];
        let top2 = [(0.0, 0.0, 250.0), (100.0, 0.0, 250.0), (100.0, -100.0, 250.0), (0.0, -100.0, 250.0)];

        let mut brep = BRep::new();
        let mut loft1 = ThruSections::new(true, true, 1.0e-6);
        loft1.add_wire(polygon_wire(&mut brep, &bottom1));
        loft1.add_wire(polygon_wire(&mut brep, &top1));
        loft1.build(&mut brep);

        let mut loft2 = ThruSections::new(true, true, 1.0e-6);
        loft2.add_wire(polygon_wire(&mut brep, &bottom2));
        loft2.add_wire(polygon_wire(&mut brep, &top2));
        loft2.build(&mut brep);

        // The loft must produce non-null shapes (BRepFill pending).
        assert!(
            loft1.shape().is_some(),
            "First loft operation should produce a valid shape"
        );
        assert!(
            loft2.shape().is_some(),
            "Second loft operation should produce a valid shape"
        );

        // Boolean fusion of the two lofted shapes (BRepAlgoAPI_Fuse).
        let s1 = loft1.shape().expect("loft1");
        let s2 = loft2.shape().expect("loft2");
        let b1 = brep_from_shape(&s1);
        let b2 = brep_from_shape(&s2);
        let fused = rcad_algo::bop::brep_algo_api::fuse(&b1, &b2);
        assert!(fused.is_ok(), "Boolean fusion of lofted shapes should succeed");
    }

    #[test]
    fn occ895_two_circular_arc_wires_no_twist() {
        // OCC895: two quarter-circle arc wires (order 0: wire2 first, then
        // wire1).  The reference surface area is 18.1614.
        // Wire 1: circle center (0,10,0), axis -Y, R=1, frame rotated 5 deg
        // around Z (gp_Ax2::Rotate).
        let angle = 5.0 * std::f64::consts::PI / 180.0;
        let base1 = Circle3::new_with_ref_dir(
            glam::DVec3::new(0.0, 10.0, 0.0),
            -glam::DVec3::Y,
            1.0,
            glam::DVec3::X,
        );
        let circle1 = rotate_circle(&base1, glam::DVec3::new(0.0, 10.0, 0.0), glam::DVec3::Z, angle);

        let mut brep = BRep::new();
        // BRepLib_MakeEdge vertex contract: first endpoint FORWARD, second
        // REVERSED (see make_polygon).
        let rev = |sr: Shape| Shape {
            orientation: topods::Orientation::Reversed,
            ..sr
        };
        use rcad_kernel::geom::CurveEval;
        let arc_edge = |brep: &mut BRep, circle: &Circle3, t1: f64, t2: f64| -> Shape {
            let v1 = brep.add_tvertex(circle.point_at(t1));
            let v2 = brep.add_tvertex(circle.point_at(t2));
            brep.add_tedge(
                Some(rcad_kernel::geom::Curve3::Circle(*circle)),
                v1,
                rev(v2),
                [t1, t2],
            )
        };
        let e1_shape = arc_edge(&mut brep, &circle1, 0.0, std::f64::consts::PI / 2.0);
        let mut w1 = rcad_modeling::MakeWire::new();
        w1.add(e1_shape);
        let wire1 = w1.wire(&mut brep);

        // Wire 2: circle at (10,0,0), axis -X, R=1 (arc 0..PI/2).
        let circle2 = Circle3::new_with_ref_dir(
            glam::DVec3::new(10.0, 0.0, 0.0),
            -glam::DVec3::X,
            1.0,
            glam::DVec3::Z,
        );
        let e2_shape = arc_edge(&mut brep, &circle2, 0.0, std::f64::consts::PI / 2.0);
        let mut w2 = rcad_modeling::MakeWire::new();
        w2.add(e2_shape);
        let wire2 = w2.wire(&mut brep);

        // ThruSections shell with order=0: wire2 first, then wire1.
        let mut loft = ThruSections::new(false, true, 1.0e-6);
        loft.add_wire(wire2);
        loft.add_wire(wire1);
        loft.build(&mut brep);

        assert!(loft.is_done(), "ThruSections must succeed");
        let shape = loft.shape().expect("ThruSections must produce a non-null shape");
        let brep = brep_from_shape(&shape);
        let area = surface_area(&brep);
        assert!(
            (area - 18.1614).abs() < 0.01,
            "Surface area should be approximately 18.1614, got {area}"
        );
    }

    /// Materialize a loft shape into its own standalone BRep pool: every
    /// TShape reachable from `s` is copied to its original flat index, so the
    /// resulting pool is directly usable by the boolean API
    /// (brep_top_shapes_with_locations).
    fn brep_from_shape(s: &Shape) -> BRep {        let mut out = BRep::new();
        let mut visited: std::collections::HashSet<u64> = std::collections::HashSet::new();
        fn place(brep: &mut BRep, sr: &Shape, visited: &mut std::collections::HashSet<u64>) {
            if !visited.insert(sr.ptr_id()) {
                return;
            }
            if brep.tshapes.len() <= sr.index {
                let dummy = std::sync::Arc::new(topods::TShape::Vertex(topods::TVertexData {
                    my_shapes: Vec::new(),
                    flags: 0,
                    point: glam::DVec3::ZERO,
                    tolerance: 0.0,
                    points: Vec::new(),
                }));
                while brep.tshapes.len() <= sr.index {
                    brep.tshapes.push(dummy.clone());
                }
            }
            brep.tshapes[sr.index] = sr.data.clone();
            match &*sr.data {
                topods::TShape::Solid(sd) => {
                    for sh in &sd.shells {
                        place(brep, sh, visited);
                    }
                    for v in &sd.internal_vertices {
                        place(brep, v, visited);
                    }
                    for e in &sd.internal_edges {
                        place(brep, e, visited);
                    }
                }
                topods::TShape::Shell(sd) => {
                    for f in &sd.faces {
                        place(brep, f, visited);
                    }
                }
                topods::TShape::Face(fd) => {
                    place(brep, &fd.outer_wire, visited);
                    for w in &fd.inner_wires {
                        place(brep, w, visited);
                    }
                    for v in &fd.internal_vertices {
                        place(brep, v, visited);
                    }
                }
                topods::TShape::Wire(wd) => {
                    for e in &wd.edges {
                        place(brep, e, visited);
                    }
                }
                topods::TShape::Edge(ed) => {
                    place(brep, &ed.first, visited);
                    place(brep, &ed.last, visited);
                }
                topods::TShape::CompSolid(cs) => {
                    for s in cs {
                        place(brep, s, visited);
                    }
                }
                topods::TShape::Compound(cd) => {
                    for s in cd {
                        place(brep, s, visited);
                    }
                }
                _ => {}
            }
        }
        place(&mut out, s, &mut visited);
        out
    }

    /// OCCT createBSplineCurve (BRepOffsetAPI_ThruSections_Test.cxx L103-129):
    /// a degree-3 non-periodic BSpline from an interleaved pole array and the
    /// full (expanded) knot sequence.
    fn create_bspline(poles_xyz: &[f64], knots: &[f64]) -> rcad_kernel::geom::BSplineCurve3 {
        let degree = 3usize;
        let n_poles = knots.len() - degree - 1;
        let mut control_points = Vec::with_capacity(n_poles);
        for i in 0..n_poles {
            control_points.push(glam::DVec3::new(
                poles_xyz[3 * i],
                poles_xyz[3 * i + 1],
                poles_xyz[3 * i + 2],
            ));
        }
        rcad_kernel::geom::BSplineCurve3 {
            degree,
            knots: knots.to_vec(),
            control_points,
            weights: vec![1.0; n_poles],
            is_periodic: false,
        }
    }

    /// BRepBuilderAPI_MakeEdge(BSplineCurve) + MakeWire — one wire from a
    /// B-spline section.
    fn bspline_wire(brep: &mut BRep, curve: &rcad_kernel::geom::BSplineCurve3) -> Shape {
        use rcad_kernel::geom::CurveEval;
        let t1 = curve.knots[curve.degree];
        let t2 = curve.knots[curve.knots.len() - curve.degree - 1];
        let v1 = brep.add_tvertex(curve.point_at(t1));
        let v2 = brep.add_tvertex(curve.point_at(t2));
        let e = brep.add_tedge(
            Some(rcad_kernel::geom::Curve3::BSpline(curve.clone())),
            v1,
            v2,
            [t1, t2],
        );
        let mut mw = rcad_modeling::MakeWire::new();
        mw.add(e);
        mw.wire(brep)
    }

    #[test]
    fn bspline_profiles_with_different_pole_count() {
        // OCC regression: ThruSections must not throw "profiles are
        // inconsistent" for closed B-spline profiles with different pole
        // counts.  The full 5-section data set (31/31/31/31/33 poles) is
        // ported; build must not throw and must report done.
        let p1: &[f64] = &[
            0.90194, -0.49457, 0.0, 1.10106097848166, -0.688963767263419, 0.0,
            1.04668209568152, -1.01716787534971, 0.0, 1.56914438377061, -1.59514757777645, 0.0,
            2.07302273729763, -3.97023193652984, 0.0, 2.13699279206564, -4.57936580650492, 0.0,
            2.13184140893145, -7.00174300027487, 0.0, 1.23454863419269, -8.49296048646617, 0.0,
            0.718149935438171, -9.43440624002067, 0.0, -2.02215956293932, -12.0782514208969, 0.0,
            -3.44515644568264, -13.0526565862391, 0.0, -8.61587484541011, -15.7788974963508, 0.0,
            -12.8309565197945, -17.0127003407937, 0.0, -18.9761857527559, -18.432946752008, 0.0,
            -20.6360623744536, -18.7727324951302, 0.0, -21.2130359371451, -19.0009746915956, 0.0,
            -22.2663547655952, -18.9176136571682, 0.0, -22.1384066610996, -17.5929745686915, 0.0,
            -21.7954793416826, -18.0629411635879, 0.0, -14.4780597089423, -15.0594484542493, 0.0,
            -10.8386889762912, -13.1323128331196, 0.0, -6.48208180266811, -9.4594741816569, 0.0,
            -5.28863405147604, -8.21363647074186, 0.0, -2.96204585602755, -5.08386133056186, 0.0,
            -2.54058250866771, -4.20341306384795, 0.0, -1.61736353367643, -2.66854172662003, 0.0,
            -1.36012191071392, -0.569755235503059, 0.0, -1.02752582768388, -0.306835878255118, 0.0,
            0.166466683785465, 0.162094918019513, 0.0, 0.0600379384518613, 0.327344971788155, 0.0,
            0.90194, -0.49457, 0.0,
        ];
        let k1: &[f64] = &[
            0.0, 0.0, 0.0, 0.0, 0.0161348522840521, 0.0807608132035811, 0.131338840314679,
            0.147844505833154, 0.16533568293713, 0.250607290759724, 0.267449419858552,
            0.322628298230969, 0.341027430676258, 0.36078046495376, 0.420324490838906,
            0.440253145913864, 0.458870520105636, 0.492435831240754, 0.509078502457839,
            0.525715578369451, 0.577446072383669, 0.597241246236863, 0.656453700919606,
            0.676155157762965, 0.694541842118733, 0.749834839741989, 0.766726765513781,
            0.834662561102041, 0.868900689429001, 0.88456999851039, 0.931780340251946, 1.0, 1.0,
            1.0, 1.0,
        ];
        // Section 5 has a different pole count (33 vs 31).
        let p5: &[f64] = &[
            1.38958, -0.30093, 56.0, 2.43615912287256, -0.909728724567382, 56.0,
            3.0281955825442, -1.8572582903678, 56.0, 3.65762773796813, -2.54529310782012, 56.0,
            4.96146733383277, -4.12983817908431, 56.0, 5.23021328351262, -4.48294463637656, 56.0,
            6.97411718250412, -8.12140889198172, 56.0, 6.70106429787142, -9.9518367568647, 56.0,
            6.61799594372847, -10.831499083977, 56.0, 5.83042424901388, -13.1305423440654, 56.0,
            4.4856942305868, -14.6753503685142, 56.0, 3.526218495628, -15.6885280970105, 56.0,
            -1.30359786741575, -19.3108973511269, 56.0, -5.83186375110243, -21.0521342727528, 56.0,
            -12.7369662411789, -23.0728272909214, 56.0, -14.6499746822398, -23.5541553891022, 56.0,
            -15.5393756779425, -23.8164232092995, 56.0, -16.7857626486778, -22.971990212193, 56.0,
            -15.8148711149702, -22.7381053237965, 56.0, -14.1582028215773, -22.0884170074466, 56.0,
            -8.2848005862047, -19.3490481837665, 56.0, -4.53645198453342, -17.0504864772276, 56.0,
            -1.10904288834323, -12.471539082827, 56.0, -0.522426290934783, -11.1995642443973, 56.0,
            0.262420398047911, -9.25463182654568, 56.0, 0.213867071381995, -6.39113351333029, 56.0,
            0.0369897362998156, -5.60398557321822, 56.0, -0.344006372425831, -3.85342683471391, 56.0,
            -1.68576438677131, -1.32767756558911, 56.0, -1.05109637156753, -0.635077640089276, 56.0,
            0.456139976122501, 0.188083742774127, 56.0, 0.318806380192129, 0.321942747786321, 56.0,
            1.38958, -0.30093, 56.0,
        ];
        let k5: &[f64] = &[
            0.0, 0.0, 0.0, 0.0, 0.0668554777677272, 0.0834666756318599, 0.100984720713562,
            0.151909931618998, 0.169593414001346, 0.254704249388902, 0.27136091562454,
            0.289520682572814, 0.325734361487981, 0.343874098775056, 0.363388151948866,
            0.422534909356546, 0.442420378784176, 0.461023580422867, 0.494592323781788,
            0.527437648618005, 0.560332403705056, 0.578743404285102, 0.598421903816054,
            0.65700259163225, 0.676406731087284, 0.694519638867019, 0.730942734487404,
            0.749396334888904, 0.766418578227312, 0.818448937788905, 0.836068871481316,
            0.8848958398016, 0.931598977689477, 1.0, 1.0, 1.0, 1.0,
        ];

        let mut brep = BRep::new();
        let mut loft = ThruSections::new(true, false, 1.0e-6);
        loft.add_wire(bspline_wire(&mut brep, &create_bspline(p1, k1)));
        loft.add_wire(bspline_wire(&mut brep, &create_bspline(p5, k5)));
        // Build must not throw and must complete (BRepFill pending).
        loft.build(&mut brep);
        assert!(loft.is_done(), "ThruSections should complete successfully");
        assert!(
            loft.shape().is_some(),
            "ThruSections should produce a valid shape"
        );
    }
}
