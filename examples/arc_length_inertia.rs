//! Example: Arc length and moment of inertia tensor.
//!
//! Demonstrates:
//!   1. Curve arc length: Line3 (analytic), Circle3 (analytic), Ellipse3 (GL16),
//!      BSplineCurve3 (GL16)
//!   2. Moment of inertia tensor: box and cylinder BReps
//!   3. Combination: sum of all 12 box-edge arc lengths to verify the box perimeter
//!
//! Run: cargo run --example phase_h_demo

use std::f64::consts::PI;

use glam::DVec3;
use rcad_kernel::{
    arc_length,
    geom::{BSplineCurve3, Circle3, Curve3, Ellipse3, Line3},
    inertia_tensor,
};
use rcad_modeling::box_brep;

// ── 1. Arc length ─────────────────────────────────────────────────────────────

fn demo_arc_length() {
    println!("\n=== 1. Curve Arc Length ===");

    // Line3: unit X direction, length = 5 for t ∈ [0, 5]
    let line = Curve3::Line(Line3 {
        origin: DVec3::ZERO,
        direction: DVec3::X,
    });
    let l_line = arc_length(&line, 0.0, 5.0);
    println!(
        "  Line3(dir=X, t=0→5):          {:.6}  (expected 5.000000)",
        l_line
    );

    // Line3 with 3-4-5 direction vector (must be unit)
    let dir_345 = DVec3::new(3.0, 4.0, 0.0).normalize();
    let line345 = Curve3::Line(Line3 {
        origin: DVec3::ZERO,
        direction: dir_345,
    });
    let l_345 = arc_length(&line345, 0.0, 7.0);
    println!(
        "  Line3(3-4 dir, t=0→7):         {:.6}  (expected 7.000000)",
        l_345
    );

    // Circle3 r=1: half circumference [0, π] → π
    let circle1 = Curve3::Circle(Circle3 {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        radius: 1.0,
    });
    let l_half = arc_length(&circle1, 0.0, PI);
    println!(
        "  Circle3(r=1, 0→π):             {:.6}  (expected {:.6})",
        l_half, PI
    );

    // Circle3 r=2: full circumference → 4π
    let circle2 = Curve3::Circle(Circle3 {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        radius: 2.0,
    });
    let l_full = arc_length(&circle2, 0.0, 2.0 * PI);
    println!(
        "  Circle3(r=2, 0→2π):            {:.6}  (expected {:.6})",
        l_full,
        4.0 * PI
    );

    // Ellipse3 a=3, b=1: quarter arc [0, π/2] — numerical
    // Known quarter-perimeter ≈ 2.4221 (from complete elliptic integral, scaled)
    // Full perimeter by Ramanujan: π(3(a+b) - sqrt((3a+b)(a+3b))) ≈ 9.6884
    let ellipse = Curve3::Ellipse(Ellipse3 {
        center: DVec3::ZERO,
        normal: DVec3::Z,
        major_dir: DVec3::X,
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let l_quarter = arc_length(&ellipse, 0.0, PI / 2.0);
    let l_ellipse_full = arc_length(&ellipse, 0.0, 2.0 * PI);
    let a = 3.0_f64;
    let b = 1.0_f64;
    let ramanujan = PI * (3.0 * (a + b) - ((3.0 * a + b) * (a + 3.0 * b)).sqrt());
    println!(
        "  Ellipse3(a=3,b=1, 0→π/2):      {:.6}  (quarter arc, GL16)",
        l_quarter
    );
    println!(
        "  Ellipse3(a=3,b=1, 0→2π):        {:.6}  (Ramanujan≈{:.6})",
        l_ellipse_full, ramanujan
    );

    // BSplineCurve3: degree-1 segment (0,0,0)→(3,4,0), length = 5
    let bspline_seg = Curve3::BSpline(BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(3.0, 4.0, 0.0)],
        weights: vec![1.0, 1.0],
    });
    let l_bspline = arc_length(&bspline_seg, 0.0, 1.0).abs();
    println!(
        "  BSpline seg (0,0,0)→(3,4,0):   {:.6}  (expected 5.000000)",
        l_bspline
    );

    // Signed arc length: reversed direction gives negative
    let l_rev = arc_length(&line, 5.0, 0.0);
    println!(
        "  Line3 reversed (t=5→0):         {:.6}  (signed, expected -5.000000)",
        l_rev
    );
}

// ── 2. Inertia tensor ─────────────────────────────────────────────────────────

fn print_tensor(label: &str, ixx: f64, iyy: f64, izz: f64) {
    println!("  {label:35}  Ixx={ixx:.4}  Iyy={iyy:.4}  Izz={izz:.4}");
}

fn demo_inertia_tensor() {
    println!("\n=== 2. Moment of Inertia Tensor ===");

    // Unit box [0,1]^3 about world origin:
    //   Ixx = ∫(y²+z²)dV = 1*(1/3+1/3) = 2/3 ≈ 0.6667
    let box_111 = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let it1 = inertia_tensor(&box_111);
    print_tensor(
        "Unit box [0,1]³  (expect 0.6667)",
        it1.ixx,
        it1.iyy,
        it1.izz,
    );
    println!(
        "    Off-diagonal: Ixy={:.4}  Ixz={:.4}  Iyz={:.4}",
        it1.ixy, it1.ixz, it1.iyz
    );
    println!("    (Off-diagonal about origin for [0,1]³: Ixy = -∫xy dV = -0.25)");

    // 2×1×1 box [0,2]×[0,1]×[0,1]:
    //   Ixx = 2*(1/3+1/3) = 4/3 ≈ 1.3333
    //   Iyy = Izz = 8/3 + 2/3 = 10/3 ≈ 3.3333
    let box_211 = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 1.0, 1.0).unwrap();
    let it2 = inertia_tensor(&box_211);
    print_tensor("Box 2×1×1 (Ixx≈1.333,Iyy≈3.333)", it2.ixx, it2.iyy, it2.izz);

    // Symmetric 2×2×2 box about its own centroid (shifted to center):
    // For a cube of side 2 centered at (1,1,1): Ixx=Iyy=Izz = m*(a²+b²)/6 = 8*(4+4)/6 = 32/3
    // Note: this box is at [0,2]^3 so about world origin Ixx = 8*(1/3+1/3) + 8*(1²+1²+1²)*m/V...
    // Just print it as a cross-check:
    let box_222 = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let it3 = inertia_tensor(&box_222);
    // Ixx = ∫₀²∫₀²∫₀² (y²+z²) dV = 8*(4/3+4/3) = 8*(8/3) = 64/3 ≈ 21.333
    print_tensor("Box 2×2×2 (expect Ixx≈21.333)", it3.ixx, it3.iyy, it3.izz);
}

// ── 3. Box edge arc lengths ───────────────────────────────────────────────────

fn demo_box_edge_arc_lengths() {
    println!("\n=== 3. Box Edge Arc Lengths (combination demo) ===");

    let w = 3.0_f64;
    let h = 2.0_f64;
    let d = 1.5_f64;
    let brep = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, w, h, d).unwrap();

    // For edges with analytic curves in GeomStore, use arc_length.
    // For box edges (Line3 segments stored as vertex pairs without explicit curves),
    // compute length from vertex positions directly.
    let mut total = 0.0_f64;
    for (ei, edge) in brep.edges.iter().enumerate() {
        let v_start = &brep.vertices[edge.start];
        let v_end = &brep.vertices[edge.end];

        let l = if let (Some(Some(ci)), Some(Some([t1, t2]))) = (
            brep.geom.edge_curve.get(ei),
            brep.geom.edge_curve_range.get(ei),
        ) {
            // Analytic curve available — use arc_length
            arc_length(&brep.geom.curves[*ci], *t1, *t2).abs()
        } else {
            // Fall back to Euclidean vertex distance (correct for straight edges)
            (v_end.point - v_start.point).length()
        };

        total += l;
        println!(
            "  Edge {ei}: v{}→v{}  length = {l:.4}",
            edge.start, edge.end
        );
    }

    // Box with w=3, h=2, d=1.5: total = 4*(3+2+1.5) = 4*6.5 = 26
    let expected = 4.0 * (w + h + d);
    println!("  Total edge length: {total:.4}  (expected {expected:.4})");
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔════════════════════════════════════════════════════╗");
    println!("║        RCAD Arc Length / Inertia Demo               ║");
    println!("║   Arc Length · Moment of Inertia Tensor            ║");
    println!("╚════════════════════════════════════════════════════╝");

    demo_arc_length();
    demo_inertia_tensor();
    demo_box_edge_arc_lengths();

    println!("\n✓ Arc length / inertia demo complete.");
}
