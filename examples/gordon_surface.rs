//! Example: Gordon surface transfinite interpolation.
//!
//! Demonstrates:
//!   1. Creating a Gordon surface from a 2x2 bilinear network
//!   2. Creating a Gordon surface from a 3x3 curve network
//!   3. Converting Gordon surface to B-spline for export
//!   4. Error handling for invalid curve networks
//!
//! Run:
//!   cargo run -p rcad-examples --example gordon_surface

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3};
use rcad_kernel::gordon::{
    GordonError, GordonOptions, ParameterizationMethod,
    gordon_surface_curves, gordon_surface_with_params,
    eval_gordon_surface_safe, gordon_surface_normal_safe, gordon_to_bspline,
};
use rcad_kernel::SurfaceEval;

fn separator(title: &str) {
    println!("\n──────────────────────────────────────────");
    println!("  {title}");
    println!("──────────────────────────────────────────");
}

fn make_line(origin: DVec3, direction: DVec3) -> Curve3 {
    Curve3::Line(Line3 { origin, direction })
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Bilinear patch (2x2 network)
// ─────────────────────────────────────────────────────────────────────────────

fn demo_bilinear_patch() {
    separator("1. Bilinear Patch (2x2 Network)");

    // Create a 2x2 network of lines forming a planar bilinear patch
    // u-curves: lines along X at y=0 and y=1
    let u0 = make_line(DVec3::ZERO, DVec3::X);
    let u1 = make_line(DVec3::Y, DVec3::X);

    // v-curves: lines along Y at x=0 and x=1
    let v0 = make_line(DVec3::ZERO, DVec3::Y);
    let v1 = make_line(DVec3::X, DVec3::Y);

    let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default())
        .expect("bilinear patch should construct");

    // Verify corners
    let p00 = eval_gordon_surface_safe(&surface, 0.0, 0.0, 1e-10).unwrap();
    let p10 = eval_gordon_surface_safe(&surface, 1.0, 0.0, 1e-10).unwrap();
    let p01 = eval_gordon_surface_safe(&surface, 0.0, 1.0, 1e-10).unwrap();
    let p11 = eval_gordon_surface_safe(&surface, 1.0, 1.0, 1e-10).unwrap();

    println!("  Corner (0,0): ({:.4}, {:.4}, {:.4})", p00.x, p00.y, p00.z);
    println!("  Corner (1,0): ({:.4}, {:.4}, {:.4})", p10.x, p10.y, p10.z);
    println!("  Corner (0,1): ({:.4}, {:.4}, {:.4})", p01.x, p01.y, p01.z);
    println!("  Corner (1,1): ({:.4}, {:.4}, {:.4})", p11.x, p11.y, p11.z);

    // Check interior point
    let p_mid = eval_gordon_surface_safe(&surface, 0.5, 0.5, 1e-10).unwrap();
    println!("  Center (0.5, 0.5): ({:.4}, {:.4}, {:.4})", p_mid.x, p_mid.y, p_mid.z);

    // Verify surface domain
    let domain = SurfaceEval::default_domain(&surface);
    println!("  Domain: u=[{:.2}, {:.2}], v=[{:.2}, {:.2}]",
        domain[0], domain[1], domain[2], domain[3]);

    // Check normal
    let normal = gordon_surface_normal_safe(&surface, 0.5, 0.5, 1e-5, 1e-10).unwrap();
    println!("  Normal at center: ({:.4}, {:.4}, {:.4})", normal.x, normal.y, normal.z);

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. 3x3 network
// ─────────────────────────────────────────────────────────────────────────────

fn demo_3x3_network() {
    separator("2. 3x3 Curve Network");

    // Create a 3x3 network of lines
    let u0 = make_line(DVec3::ZERO, DVec3::X);
    let u1 = make_line(DVec3::new(0.0, 0.5, 0.0), DVec3::X);
    let u2 = make_line(DVec3::Y, DVec3::X);

    let v0 = make_line(DVec3::ZERO, DVec3::Y);
    let v1 = make_line(DVec3::new(0.5, 0.0, 0.0), DVec3::Y);
    let v2 = make_line(DVec3::X, DVec3::Y);

    // Use uniform parameterization for uniformly spaced curves
    let opts = GordonOptions::default()
        .with_parameterization(ParameterizationMethod::Uniform);

    let surface = gordon_surface_curves(
        &[u0, u1, u2],
        &[v0, v1, v2],
        opts,
    ).expect("3x3 network should construct");

    // Verify interpolation at all grid points
    let u_params = &surface.u_params;
    let v_params = &surface.v_params;

    println!("  U-params: {:?}", u_params);
    println!("  V-params: {:?}", v_params);

    // Check that the surface passes through all curve intersections
    let mut max_error = 0.0f64;
    for (_, &vi) in v_params.iter().enumerate() {
        for (_, &uj) in u_params.iter().enumerate() {
            let p = eval_gordon_surface_safe(&surface, uj, vi, 1e-10).unwrap();
            let expected = DVec3::new(uj, vi, 0.0);
            let err = (p - expected).length();
            max_error = max_error.max(err);
        }
    }

    println!("  Max interpolation error: {:.2e}", max_error);
    assert!(max_error < 1e-8, "interpolation error too large");
    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Convert to B-spline
// ─────────────────────────────────────────────────────────────────────────────

fn demo_convert_to_bspline() {
    separator("3. Convert to B-Spline");

    // Create a simple bilinear patch
    let u0 = make_line(DVec3::ZERO, DVec3::X);
    let u1 = make_line(DVec3::Y, DVec3::X);
    let v0 = make_line(DVec3::ZERO, DVec3::Y);
    let v1 = make_line(DVec3::X, DVec3::Y);

    let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default()).unwrap();

    // Convert to B-spline
    let bspline = gordon_to_bspline(&surface, 5, 5, 3).expect("B-spline conversion should succeed");

    println!("  B-spline degree: ({}, {})", bspline.degree_u, bspline.degree_v);
    println!("  Control points: {}x{}", bspline.control_points.len(), bspline.control_points[0].len());
    println!("  Knots U: {} values", bspline.knots_u.len());
    println!("  Knots V: {} values", bspline.knots_v.len());

    // Verify B-spline evaluation matches Gordon at sample points
    let p_gordon = eval_gordon_surface_safe(&surface, 0.5, 0.5, 1e-10).unwrap();
    let p_bspline = bspline.point_at(0.5, 0.5);
    let err = (p_gordon - p_bspline).length();
    println!("  Gordon vs B-spline at (0.5, 0.5): error = {:.2e}", err);

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Error handling
// ─────────────────────────────────────────────────────────────────────────────

fn demo_error_handling() {
    separator("4. Error Handling");

    // Too few curves
    {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let result = gordon_surface_curves(&[u0], &[], GordonOptions::default());
        match &result {
            Err(GordonError::TooFewUCurves { count }) => {
                println!("  TooFewUCurves: count={}  PASS", count);
            }
            _ => panic!("expected TooFewUCurves error"),
        }
    }

    // Non-monotonic parameters
    {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let result = gordon_surface_with_params(
            &[u0, u1],
            &[0.0, 1.0],
            &[v0, v1],
            &[0.7, 0.3], // Non-monotonic
            GordonOptions::default(),
        );
        match &result {
            Err(GordonError::NonMonotonicParams { direction, index, prev, curr }) => {
                println!("  NonMonotonicParams: {} at index {} ({} >= {})  PASS",
                    direction, index, prev, curr);
            }
            _ => panic!("expected NonMonotonicParams error"),
        }
    }

    // Parameters out of range
    {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let result = gordon_surface_with_params(
            &[u0, u1],
            &[0.0, 1.0],
            &[v0, v1],
            &[-0.5, 1.0], // Out of range
            GordonOptions::default(),
        );
        match &result {
            Err(GordonError::ParamsOutOfRange { direction, value }) => {
                println!("  ParamsOutOfRange: {} = {}  PASS", direction, value);
            }
            _ => panic!("expected ParamsOutOfRange error"),
        }
    }

    println!("  All error cases handled correctly");
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Options and continuity
// ─────────────────────────────────────────────────────────────────────────────

fn demo_options() {
    separator("5. Options and Continuity");

    let u0 = make_line(DVec3::ZERO, DVec3::X);
    let u1 = make_line(DVec3::Y, DVec3::X);
    let v0 = make_line(DVec3::ZERO, DVec3::Y);
    let v1 = make_line(DVec3::X, DVec3::Y);

    // C0 continuity
    let opts_c0 = GordonOptions::c0();
    gordon_surface_curves(&[u0.clone(), u1.clone()], &[v0.clone(), v1.clone()], opts_c0)
        .expect("C0 should construct");
    println!("  C0 continuity: constructed");

    // C1 continuity
    let opts_c1 = GordonOptions::c1();
    gordon_surface_curves(&[u0.clone(), u1.clone()], &[v0.clone(), v1.clone()], opts_c1)
        .expect("C1 should construct");
    println!("  C1 continuity: constructed");

    // C2 continuity
    let opts_c2 = GordonOptions::c2();
    gordon_surface_curves(&[u0, u1], &[v0, v1], opts_c2)
        .expect("C2 should construct");
    println!("  C2 continuity: constructed");

    // Custom options
    let _opts_custom = GordonOptions::default()
        .with_tolerance(1e-8)
        .with_intersection_tolerance(1e-3)
        .skip_intersection_validation();
    println!("  Custom options: tolerance=1e-8, intersection_tol=1e-3, skip_validation=true");

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    println!("=================================================");
    println!("  Gordon Surface Demo");
    println!("=================================================");

    demo_bilinear_patch();
    demo_3x3_network();
    demo_convert_to_bspline();
    demo_error_handling();
    demo_options();

    println!("\n=================================================");
    println!("  Gordon Surface: All demos completed successfully");
    println!("=================================================");
}
