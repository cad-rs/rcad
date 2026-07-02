#[cfg(test)]

mod tests {
    use super::*;
    use crate::{boolean_op, BooleanOpType};
    use glam::DAffine3;
    use rcad_kernel::surface_area;
    use rcad_kernel::volume;

    #[test]
    fn concentric_sphere_difference_analytic_shell_surface_area() {
        let center = DVec3::new(1.0, -2.0, 4.0);
        let ro = 5.0_f64;
        let ri = 3.0_f64;
        let outer = make_sphere_brep(center, ro).expect("outer");
        let inner = make_sphere_brep(center, ri).expect("inner");
        let shell = boolean_op(BooleanOpType::Difference, &outer, &inner).expect("difference");
        let pi = std::f64::consts::PI;
        let a_ex = 4.0 * pi * (ro * ro + ri * ri);
        let a = surface_area(&shell);
        assert!(
            (a - a_ex).abs() < 50.0 * TOLERANCE_RETRY_LADDER_COARSE * a_ex.max(1.0),
            "surface area {a} vs analytic shell SA {a_ex}"
        );
        // `signed_volume` across compounds relies on consistent face normals vs tessellation;
        // sphere primitives carry approximate face normals 闁?SA matches analytic \(4\pi(R^2+r^2)\).
    }

    #[test]
    fn zp3_sum_planar_areas_before_sew_matches_expected_total() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let z_lo = -10.0_f64;
        let z_hi = 0.0_f64;
        let rc = 10.0_f64;
        let z0 = -10.0_f64;
        let z1 = 0.0_f64;
        let r0 = 10.0_f64;
        let r1 = 5.0_f64;
        let mut outer_bot = Vec::with_capacity(N);
        let mut outer_top = Vec::with_capacity(N);
        let mut inner_bot = Vec::with_capacity(N);
        let mut inner_top = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            let c = ang.cos();
            let s = ang.sin();
            let ob = DVec3::new(rc * c, rc * s, z_lo);
            outer_bot.push(ob);
            inner_bot.push(ob);
            outer_top.push(DVec3::new(rc * c, rc * s, z_hi));
            inner_top.push(DVec3::new(r1 * c, r1 * s, z1));
        }
        let ann = annulus_between_rings(&outer_top, &inner_top).unwrap();
        let (lo, oh) = loft_with_history(&[outer_bot, outer_top]).unwrap();
        let outer_strip = strip_loft_caps(lo, oh).unwrap();
        let (li, ih) = loft_with_history(&[inner_bot, inner_top]).unwrap();
        let inner_strip = strip_loft_caps(li, ih).unwrap();
        let sum = rcad_kernel::surface_area(&outer_strip)
            + rcad_kernel::surface_area(&inner_strip)
            + rcad_kernel::surface_area(&ann);
        let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 1390.8_f64);
        assert!(
            (sum - 1390.8_f64).abs() <= tol,
            "sum loose pieces ~1390.8, got {sum}"
        );
    }

    #[test]
    fn zp3_inner_frustum_strip_surface_area_sane() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let mut ib = Vec::with_capacity(N);
        let mut it = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            let c = ang.cos();
            let s = ang.sin();
            ib.push(DVec3::new(10.0 * c, 10.0 * s, -10.0));
            it.push(DVec3::new(5.0 * c, 5.0 * s, 0.0));
        }
        let (loft_i, ih) = loft_with_history(&[ib, it]).unwrap();
        let strip = strip_loft_caps(loft_i, ih).unwrap();
        let a = rcad_kernel::surface_area(&strip);
        assert!(
            a > 400.0 && a < 650.0,
            "inner frustum mantle area ~527, got {a}"
        );
    }

    #[test]
    fn zp3_outer_face_area_matches_before_and_after_sew() {
        let (outer_strip, inner_strip, annulus) = coaxial_cylinder_minus_frustum_loft_pieces(
            -10.0, 0.0, 10.0, 10.0, -10.0, 10.0,
        )
        .expect("loft pieces");
        let tol = (TOLERANCE_RETRY_LADDER_COARSE).max(TOLERANCE_MESH_LEGACY * 10.0);
        let sewn = sew_shells(&[outer_strip.clone(), inner_strip, annulus], tol);
        assert!(sewn.free_edges.is_empty(), "free {:?}", sewn.free_edges);
        let f_loose = &outer_strip.solids[0].shells[0].faces[0];
        let f_sewn = &sewn.brep.solids[0].shells[0].faces[0];
        let a_loose = rcad_kernel::face_surface_area(&outer_strip, f_loose, 0);
        let a_sewn = rcad_kernel::face_surface_area(&sewn.brep, f_sewn, 0);
        assert!(
            (a_loose - a_sewn).abs() < TOLERANCE_RETRY_LADDER_COARSE * 100.0,
            "first outer lateral face area loose {a_loose} vs sewn {a_sewn}"
        );
    }

    #[test]
    fn zp3_loft_shell_matches_occt_geometry_numbers() {
        // Cone apex z=10, base z=-10, rb=10; cylinder z in [-10,0], r=10 闁?same as geometry_properties ZP3.
        let r = try_coaxial_cylinder_minus_frustum_loft_shell(-10.0, 0.0, 10.0, 10.0, -10.0, 10.0);
        assert!(r.is_some(), "expected sewn loft shell for ZP3 parameters");
        let brep = r.unwrap();
        let nf: usize = brep
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count();
        let a = rcad_kernel::surface_area(&brep);
        let v = rcad_kernel::volume(&brep);
        assert!(
            nf >= 60 && nf <= 70,
            "expected ~65 faces (32+32+annulus), got {nf}"
        );
        assert!(
            (v - 1310.0).abs() < 80.0,
            "expected volume cylinder minus frustum ~1310, got {v}"
        );
        let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.0625 * 1390.8_f64);
        assert!(
            (a - 1390.8_f64).abs() <= tol,
            "surface area: expected ~1390.8, got {a} (nf={nf}, vol={v})"
        );
    }

    #[test]
    fn zp3_annulus_plane_builds() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let mut o = Vec::with_capacity(N);
        let mut inn = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            o.push(DVec3::new(10.0 * ang.cos(), 10.0 * ang.sin(), 0.0));
            inn.push(DVec3::new(5.0 * ang.cos(), 5.0 * ang.sin(), 0.0));
        }
        annulus_between_rings(&o, &inn).expect("annulus");
    }

    #[test]
    fn zp3_coaxial_cylinder_minus_cone_fast_path_triggered() {
        use glam::{DAffine3, DVec2, DVec3};
        use rcad_modeling::{make_cone_brep, make_cylinder_brep};
        let pc = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 20.0).expect("cone");
        let pcy = make_cylinder_brep(
            DVec3::new(0.0, 0.0, -5.0),
            DVec3::Z,
            DVec3::X,
            10.0,
            10.0,
        )
        .expect("cylinder");
        assert!(
            try_difference_coaxial_cylinder_minus_cone(&pcy, &pc).is_some(),
            "ZP3 boptuc expects coaxial cylinder\\cone shortcut"
        );
    }

    // 閳光偓閳光偓 box-box intersection fast path 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

    #[test]
    fn box_box_intersection_partial_overlap() {
        // bcommon_simple_c1: 1鑴?鑴? 閳?1.5鑴?.5鑴?.5 閳?1鑴?.5鑴?.5 (SA=2.5, vol=0.25)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.5, 0.5, 0.5).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("box-box intersection");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 2.5).abs() < 1e-7, "SA={sa} expected 2.5");
        assert!((vol - 0.25).abs() < 1e-7, "vol={vol} expected 0.25");
    }

    #[test]
    fn box_box_intersection_full_containment() {
        // 2鑴?鑴? box 閳?0.5鑴?.5鑴?.5 inside 閳?the inner 0.5^3 box
        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let inner = make_box_brep(
            DVec3::new(0.25, 0.25, 0.25),
            DVec3::X,
            DVec3::Y,
            0.5,
            0.5,
            0.5,
        )
        .unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &outer, &inner).expect("containment");
        let vol = volume(&r);
        assert!((vol - 0.125).abs() < 1e-7, "vol={vol} expected 0.125");
    }

    #[test]
    fn box_box_intersection_no_overlap() {
        // Disjoint boxes 閳?empty intersection.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("no-overlap");
        let n_faces: usize = r.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty intersection");
    }

    #[test]
    fn box_box_intersection_a7_like() {
        // bcommon_simple A7: 1鑴?鑴? 閳?1鑴?.5鑴? 閳?the contained 1鑴?鑴? (SA=6, vol=1).
        // Tests that try_containment returns inner (not outer) when smaller operand
        // is passed first (swapped=true).
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.5, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &a, &b).expect("box-box A7");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 6.0).abs() < 1e-7, "SA={sa} expected 6.0");
        assert!((vol - 1.0).abs() < 1e-7, "vol={vol} expected 1.0");
    }

    #[test]
    fn box_box_intersection_non_box_falls_through() {
        // Sphere 閳?box falls through to generic path (no panic).
        let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &sphere, &b).expect("sphere-box");
        // Some result 閳?specific value doesn't matter.
        assert!(r.solids.len() >= 1 || r.vertices.is_empty());
    }

    // 閳光偓閳光偓 box-box difference fast path 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

    #[test]
    fn box_box_difference_opposite_same_axis() {
        // F1-like: A=1.5鑴?.5鑴?.5 at (-0.25,0,0), B=1鑴?鑴? at origin.
        // boptuc = A - B (but here we test try_difference_box_box directly,
        // so a = first arg = A, b = second arg = B).
        // A extends on x-lo (-0.25<0) and x-hi (1.25>1), same axis 閳?two
        // disjoint slabs: SA=2, vol=0.125.
        let a = make_box_brep(
            DVec3::new(-0.25, 0.0, 0.0), DVec3::X, DVec3::Y,
            1.5, 0.5, 0.5,
        ).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("box-box difference");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 2.0).abs() < 1e-7, "SA={sa} expected 2.0");
        assert!((vol - 0.125).abs() < 1e-7, "vol={vol} expected 0.125");
    }

    #[test]
    fn box_box_difference_single_slab() {
        // bcut_simple_c1: 0.5鑴?.5鑴? at (0,-0.5,0) minus 1鑴?鑴? at origin.
        // Excess only on y-lo [-0.5, 0]. Result: 0.5鑴?.5鑴? box, SA=2.5.
        let a = make_box_brep(DVec3::new(0.0, -0.5, 0.0), DVec3::X, DVec3::Y, 0.5, 1.5, 1.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("box-box difference");
        let sa = surface_area(&r);
        let vol = volume(&r);
        assert!((sa - 2.5).abs() < 1e-7, "SA={sa} expected 2.5");
        assert!((vol - 0.25).abs() < 1e-7, "vol={vol} expected 0.25");
    }

    #[test]
    fn box_box_difference_no_overlap() {
        // Disjoint boxes 閳?difference is just A.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &b).expect("no-overlap");
        let sa = surface_area(&r);
        assert!((sa - 6.0).abs() < 1e-7, "SA={sa} expected 6.0 (unchanged A)");
    }

    #[test]
    fn box_box_difference_a_inside_b() {
        // A fully inside B 閳?empty.
        let a = make_box_brep(DVec3::new(0.25, 0.25, 0.25), DVec3::X, DVec3::Y, 0.5, 0.5, 0.5).unwrap();
        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let r = boolean_op(BooleanOpType::Difference, &a, &outer).expect("A-in-B");
        let n_faces: usize = r.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty difference");
    }

    // 閳光偓閳光偓 general box-box (rotated) 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

    fn make_rotated_box(origin: DVec3, u_dir: DVec3, v_dir: DVec3, w: f64, h: f64, d: f64, pivot: DVec3, axis: DVec3, angle_deg: f64) -> BRep {
        // Handle OCCT-style negative extents: a negative extent means the box extends
        // in the negative direction from the anchor corner.
        let z_dir = u_dir.cross(v_dir);
        let mut o = origin;
        let ww = if w < 0.0 { o += u_dir * w; -w } else { w };
        let hh = if h < 0.0 { o += v_dir * h; -h } else { h };
        let dd = if d < 0.0 { o += z_dir * d; -d } else { d };
        let mut b = make_box_brep(o, u_dir, v_dir, ww, hh, dd).unwrap();
        let rot = DAffine3::from_axis_angle(axis.normalize(), angle_deg.to_radians());
        let xf = DAffine3::from_translation(pivot) * rot * DAffine3::from_translation(-pivot);
        b.apply_transform(xf);
        b
    }

    #[test]
    fn box_detection_axis_aligned() {
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
        let info = try_as_box(&b).expect("axis-aligned box should be detected");
        // Axes may be in any order; check that all 3 standard axes are present.
        let all_axes: Vec<DVec3> = info.axes.iter().map(|a| a.abs()).collect();
        assert!(all_axes.iter().any(|a| (a - DVec3::X).length() < 1e-10), "X axis missing");
        assert!(all_axes.iter().any(|a| (a - DVec3::Y).length() < 1e-10), "Y axis missing");
        assert!(all_axes.iter().any(|a| (a - DVec3::Z).length() < 1e-10), "Z axis missing");
        assert!((info.center - DVec3::new(1.0, 1.5, 2.0)).length() < 1e-10, "center");
        // extents in same order as axes; check all three match {1.0, 1.5, 2.0}.
        let mut ex = info.extents.to_vec();
        ex.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((ex[0] - 1.0).abs() < 1e-10 && (ex[1] - 1.5).abs() < 1e-10 && (ex[2] - 2.0).abs() < 1e-10,
            "extents {:?}", info.extents);
    }

    #[test]
    fn box_detection_rotated() {
        // Box at origin, rotated 45鎺?around Z at origin.
        let b = make_box_brep(DVec3::new(-0.5, -0.5, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = {
            let mut shape = b;
            let rot = DAffine3::from_axis_angle(DVec3::Z, 45.0_f64.to_radians());
            shape.apply_transform(rot);
            shape
        };
        let info = try_as_box(&b).expect("rotated box should be detected");
        let expected_axes = [
            DVec3::new(0.7071067811865476, 0.7071067811865476, 0.0).abs(),
            DVec3::new(-0.7071067811865476, 0.7071067811865476, 0.0).abs(),
            DVec3::new(0.0, 0.0, 1.0).abs(),
        ];
        for a in &info.axes {
            let aa = a.abs();
            let found = expected_axes.iter().any(|e| (aa - e).length() < 1e-10);
            assert!(found, "unexpected axis {:?}", a);
        }
        let planes = info.planes();
        assert_eq!(planes.len(), 6, "should have 6 half-space planes");
    }

    #[test]
    fn box_detection_non_box() {
        let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        assert!(try_as_box(&sphere).is_none(), "sphere is not a box");
    }

    #[test]
    fn rotated_box_intersection_partial_overlap() {
        // bcommon_simple_c3-like: unit box 閳?rotated box.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = (2.0_f64).sqrt();
        let b = make_rotated_box(
            DVec3::ZERO, DVec3::X, DVec3::Y, r, r / 2.0, 1.0,
            DVec3::ZERO, DVec3::Z, 45.0,
        );
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).expect("rotated intersection");
        let sa = surface_area(&result);
        let tol = (5e-3_f64).max(0.15 * sa);
        // SA should be non-zero (boxes overlap).
        assert!(sa > 0.0, "expected non-empty intersection, SA={sa}");
        // Check that try_intersection_box_general was triggered.
        let direct = try_intersection_box_general(&a, &b);
        assert!(direct.is_some(), "general box-box intersection should fire");
    }

    #[test]
    fn rotated_box_difference_boptuc_c3() {
        // boptuc_simple C3: B - A where A = unit box, B = rotated box.
        // Expected SA = 5.82843.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let r = (2.0_f64).sqrt();
        let b = {
            let mut shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, r, r / 2.0, 1.0).unwrap();
            let rot = DAffine3::from_axis_angle(DVec3::Z, 45.0_f64.to_radians());
            shape.apply_transform(rot);
            shape
        };
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("boptuc C3");
        let sa = surface_area(&result);
        let expected = 5.82843;
        let tol = (5e-3_f64).max(0.15 * expected);
        assert!(
            (sa - expected).abs() <= tol,
            "C3: expected SA ~{expected}, got {sa} (tol={tol})"
        );
    }

    #[test]
    fn rotated_box_difference_boptuc_n3() {
        // boptuc_simple N3: B - A where B is a rotated box with offset.
        // Expected SA = 2.5.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        // box at (0.25, 0.25, 0) size 0.5鑴?.5鑴?1, pivot (.25,.25,0), rotate 30鎺?Z
        let b = make_rotated_box(
            DVec3::new(0.25, 0.25, 0.0), DVec3::X, DVec3::Y, 0.5, 0.5, -1.0,
            DVec3::new(0.25, 0.25, 0.0), DVec3::Z, 30.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("boptuc N3");
        let sa = surface_area(&result);
        let expected = 2.5;
        let tol = (5e-3_f64).max(0.15 * expected);
        assert!(
            (sa - expected).abs() <= tol,
            "N3: expected SA ~{expected}, got {sa} (tol={tol})"
        );
    }

    #[test]
    fn rotated_box_difference_boptuc_o1() {
        // boptuc_simple O1: B - A with rotated B at offset.
        // Expected SA = 4.48.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_rotated_box(
            DVec3::new(0.0, 0.5, 0.0), DVec3::X, DVec3::Y, 0.8, 0.8, -1.0,
            DVec3::new(0.0, 0.5, 0.0), DVec3::Z, -45.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("boptuc O1");
        let sa = surface_area(&result);
        let expected = 4.48;
        let tol = (5e-3_f64).max(0.15 * expected);
        assert!(
            (sa - expected).abs() <= tol,
            "O1: expected SA ~{expected}, got {sa} (tol={tol})"
        );
    }

    #[test]
    fn rotated_box_difference_no_overlap() {
        // Disjoint rotated boxes 閳?difference should be B unchanged.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_rotated_box(
            DVec3::new(5.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
            DVec3::new(5.0, 0.0, 0.0), DVec3::Z, 30.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("disjoint diff");
        let sa = surface_area(&result);
        // Rotated unit box has same SA = 6.0.
        assert!((sa - 6.0).abs() < 1e-6, "disjoint: expected SA=6.0, got {sa}");
    }

    #[test]
    fn rotated_box_difference_a_contains_b() {
        // B fully inside A 閳?B - A = empty.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_rotated_box(
            DVec3::new(0.25, 0.25, 0.25), DVec3::X, DVec3::Y, 0.5, 0.5, 0.5,
            DVec3::new(0.25, 0.25, 0.25), DVec3::Z, 15.0,
        );
        let result = boolean_op(BooleanOpType::Difference, &b, &a).expect("A contains B");
        let n_faces: usize = result.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        assert_eq!(n_faces, 0, "expected empty (B inside A)");
    }

    #[test]
    fn rotated_box_intersection_non_box_falls_through() {
        // Box 閳?sphere 閳?falls through to Pave-Filler (no panic).
        let b = make_rotated_box(
            DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0,
            DVec3::ZERO, DVec3::Z, 45.0,
        );
        let s = make_sphere_brep(DVec3::new(0.5, 0.5, 0.5), 1.0).unwrap();
        let r = boolean_op(BooleanOpType::Intersection, &b, &s).expect("box-sphere intersection");
        assert!(r.solids.len() >= 1 || r.vertices.is_empty());
    }

    #[test]
    fn debug_half_cylinder_sa() {
        // U5 case: cylinder r=1, h=2, center (0,0,1), clip X閳?
        let brep = build_half_cylinder_intersection_brep(
            DVec3::new(0.0, 0.0, 1.0), 1.0, 2.0, DVec3::X, 0.0,
        );
        let total = rcad_kernel::surface_area(&brep);
        println!("DEBUG half-cylinder: total SA = {total}");

        // Per-face SA
        for (si, solid) in brep.solids.iter().enumerate() {
            for (fi, face) in solid.shells[0].faces.iter().enumerate() {
                let a = rcad_kernel::face_surface_area(&brep, face, fi);
                let n_edges = face.outer_wire.edges.len();
                println!("  Solid {si} Face {fi} ({n_edges} edges): SA = {a}");
            }
        }

        let expected = 3.0 * std::f64::consts::PI + 4.0; // 2锜?+ 锜?2 + 锜?2 + 4 = 3锜?4
        println!("DEBUG expected = {expected}");
        assert!(
            (total - expected).abs() < 0.01,
            "Expected ~{expected}, got {total}"
        );
    }


    fn count_face_surface_kinds(brep: &BRep) -> (usize, usize) {
        let mut n_cyl = 0;
        let mut n_plan = 0;
        let mut global_fi = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for _face in &shell.faces {
                    let Some(si) = brep.geom.face_surface.get(global_fi).and_then(|s| *s) else {
                        global_fi += 1;
                        continue;
                    };
                    match &brep.geom.surfaces[si] {
                        Surface3::Cylinder(_) => n_cyl += 1,
                        Surface3::Plane(_) => n_plan += 1,
                        _ => {}
                    }
                    global_fi += 1;
                }
            }
        }
        (n_cyl, n_plan)
    }

    #[test]
    fn cylinder_box_intersection_splits_half_cylinder_at_internal_tangent() {
        let cyl = make_cylinder_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
        let bx = make_box_brep(DVec3::new(-1.0, -1.0, 0.0), DVec3::X, DVec3::Y, 3.0, 1.0, 1.0).unwrap();
        let brep = try_intersect_cylinder_box_one_dir(&cyl, &bx).expect("half-cylinder tangent split");
        let n_faces: usize = brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        let (n_cyl, n_plan) = count_face_surface_kinds(&brep);
        assert_eq!((n_faces, n_cyl, n_plan), (5, 2, 3));
    }

    #[test]
    fn cylinder_box_intersection_splits_full_cylinder_at_tangent_generators() {
        let cyl = make_cylinder_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
        let bx = make_box_brep(DVec3::new(-1.0, -1.0, 0.0), DVec3::X, DVec3::Y, 2.5, 3.0, 1.0).unwrap();
        let brep = try_intersect_cylinder_box_one_dir(&cyl, &bx).expect("full-cylinder tangent split");
        let n_faces: usize = brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count();
        let (n_cyl, n_plan) = count_face_surface_kinds(&brep);
        assert_eq!((n_faces, n_cyl, n_plan), (5, 3, 2));
    }

    #[test]
    fn cylinder_box_intersection_external_tangent_is_empty() {
        let cyl = make_cylinder_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
        let bx = make_box_brep(DVec3::new(1.0, -1.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let brep = try_intersection_cylinder_box(&cyl, &bx).expect("external tangent handled by cylinder-box fast path");
        assert!(brep.solids.is_empty(), "expected empty tangent intersection, got {} solids", brep.solids.len());
    }

}

/// Build a cylinder split into 4 lateral quadrant faces + 2 planar caps.
fn build_cylinder_quadrant_brep(
    center: DVec3, axis: DVec3, ref_dir: DVec3, r: f64, h: f64,
) -> BRep {
    use std::f64::consts::PI;
    let half_h = h * 0.5;
    let v_dir = axis.cross(ref_dir).normalize();
    let lo = center - half_h * axis;
    let hi = center + half_h * axis;
    let mut brep = BRep::new();
    macro_rules! v { ($p:expr) => {{ let i = brep.vertices.len(); brep.vertices.push(Vertex { point: $p }); i }}; }
    let angles = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0];
    let mut b = [0usize; 4]; let mut t = [0usize; 4];
    for (i, &ang) in angles.iter().enumerate() {
        let p_on_cyl = |z: DVec3| -> DVec3 { z + r * (ang.cos() * ref_dir + ang.sin() * v_dir) };
        b[i] = v!(p_on_cyl(lo)); t[i] = v!(p_on_cyl(hi));
    }
    macro_rules! e {
        ($s:expr, $e:expr, $curve:expr, $t0:expr, $t1:expr) => {{
            let idx = brep.edges.len();
            brep.edges.push(Edge { start: $s, end: $e });
            let ci = brep.geom.curves.len(); brep.geom.curves.push($curve);
            while brep.geom.edge_curve.len() <= idx {
                brep.geom.edge_curve.push(None); brep.geom.edge_curve_range.push(None);
                brep.geom.edge_degenerated.push(false);
            }
            brep.geom.edge_curve[idx] = Some(ci);
            brep.geom.edge_curve_range[idx] = Some([$t0, $t1]);
            brep.geom.edge_pcurves.push(Vec::new());
            idx
        }};
    }
    let mut gens = [0usize; 4];
    for i in 0..4 {
        let dir = (brep.vertices[t[i]].point - brep.vertices[b[i]].point).normalize();
        gens[i] = e!(b[i], t[i],
            Curve3::Line(Line3 { origin: brep.vertices[b[i]].point, direction: dir }), 0.0, h);
    }
    let cyl_surf_3 = Surface3::Cylinder(CylindricalSurface { origin: lo, axis, ref_dir, radius: r });
    let bot_c = lo; let top_c = hi;
    let mut bot_arcs = [0usize; 4]; let mut top_arcs = [0usize; 4];
    for i in 0..4 {
        let j = (i + 1) % 4;
        let a0 = angles[i]; let a1 = angles[j];
        bot_arcs[i] = e!(b[i], b[j],
            Curve3::Circle(Circle3::new(bot_c, axis, r )), a0, a1);
        top_arcs[i] = e!(t[i], t[j],
            Curve3::Circle(Circle3::new(top_c, axis, r )), a0, a1);
    }
    let cyl_id = brep.geom.surfaces.len(); brep.geom.surfaces.push(cyl_surf_3.clone());
    let bp_id = brep.geom.surfaces.len(); brep.geom.surfaces.push(Surface3::Plane(Plane { origin: bot_c, normal: -axis }));
    let tp_id = brep.geom.surfaces.len(); brep.geom.surfaces.push(Surface3::Plane(Plane { origin: top_c, normal: axis }));
    let mut faces = Vec::new();
    for i in 0..4 {
        let j = (i + 1) % 4;
        let wire_edges = vec![
            WireEdge::fwd(gens[i]), WireEdge::rev(top_arcs[i]),
            WireEdge::rev(gens[j]), WireEdge::fwd(bot_arcs[i]),
        ];
        let fs = brep.geom.surfaces.len(); brep.geom.surfaces.push(cyl_surf_3.clone());
        brep.geom.face_surface.push(Some(fs)); brep.geom.face_surface_range.push(None);
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: DVec3::ZERO, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None });
    }
    // Bottom cap (normal = -axis, outward)
    {
        let wire_edges: Vec<WireEdge> = (0..4).map(|i| WireEdge::fwd(bot_arcs[i])).collect();
        brep.geom.face_surface.push(Some(bp_id)); brep.geom.face_surface_range.push(None);
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: -axis, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None });
    }
    // Top cap (normal = axis, outward)
    {
        let wire_edges: Vec<WireEdge> = (0..4).map(|i| WireEdge::fwd(top_arcs[i])).collect();
        brep.geom.face_surface.push(Some(tp_id)); brep.geom.face_surface_range.push(None);
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: axis, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None });
    }
    brep.solids.push(Solid { shells: vec![Shell { faces }] });
    brep
}
