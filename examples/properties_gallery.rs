//! Example: Properties gallery — surface area, volume, centroid for many shapes.
//!
//! Prints a comparison table for boxes, prisms, beams, and other extruded shapes.
//!
//! Run: cargo run --example properties_gallery

use glam::DVec3;
use rcad_algorithms::check;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::{BRep, PrimitiveSolid, centroid, surface_area, volume};
use rcad_modeling::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_modeling::extrude;

fn main() {
    let mut rows: Vec<Row> = Vec::new();

    // ── Primitive boxes ────────────────────────────────────────────────
    // Analytic formulas: SA = 2(wh+wd+hd), V = whd, C = (w/2, h/2, d/2)
    let b1 = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    add(&mut rows, "Box 1×1×1", "SA=6, V=1", &b1);

    let b2 = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 3.0,
        depth: 4.0,
    });
    add(&mut rows, "Box 2×3×4", "SA=52, V=24", &b2);

    let b3 = BRep::from_primitive(PrimitiveSolid::Box {
        width: 5.0,
        height: 0.5,
        depth: 0.5,
    });
    add(&mut rows, "Flat bar 5×½×½", "SA=12.5, V=1.25", &b3);

    // ── Square / rectangular prisms ───────────────────────────────────
    let sq = polygon_profile(&[
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
    ]);
    add(
        &mut rows,
        "Square pillar 1×1×3",
        "SA=14, V=3",
        &extrude(&sq, 0, DVec3::Z, 3.0).unwrap(),
    );

    let wide = polygon_profile(&[
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(4.0, 0.0, 0.0),
        DVec3::new(4.0, 0.2, 0.0),
        DVec3::new(0.0, 0.2, 0.0),
    ]);
    add(
        &mut rows,
        "Flat plank 4×0.2×1",
        "SA=11.2, V=0.8",
        &extrude(&wide, 0, DVec3::Z, 1.0).unwrap(),
    );

    // ── Equilateral triangle prism ─────────────────────────────────────
    // Equilateral triangle with side=2: area=√3, perimeter=6
    let tri = polygon_profile(&[
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(1.0, 1.732, 0.0),
    ]);
    add(
        &mut rows,
        "Eq. prism side=2 h=4",
        "SA≈27.46, V≈6.93",
        &extrude(&tri, 0, DVec3::Z, 4.0).unwrap(),
    );

    // ── Hexagon bolt head ──────────────────────────────────────────────
    // Regular hexagon r=1: area=3√3/2≈2.598, perimeter=6
    let hex_pts: Vec<DVec3> = (0..6)
        .map(|i| {
            let a = std::f64::consts::FRAC_PI_3 * i as f64;
            DVec3::new(a.cos(), a.sin(), 0.0)
        })
        .collect();
    let hex = polygon_profile(&hex_pts);
    add(
        &mut rows,
        "Hex bolt head r=1 h=0.6",
        "SA≈8.80, V≈1.56",
        &extrude(&hex, 0, DVec3::Z, 0.6).unwrap(),
    );

    // ── L-beam ────────────────────────────────────────────────────────
    // L: 3×3 outer, t=0.5 legs; cross-section = 3*0.5 + 2.5*0.5 = 2.75
    let l_pts = [
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(3.0, 0.0, 0.0),
        DVec3::new(3.0, 0.5, 0.0),
        DVec3::new(0.5, 0.5, 0.0),
        DVec3::new(0.5, 3.0, 0.0),
        DVec3::new(0.0, 3.0, 0.0),
    ];
    add(
        &mut rows,
        "L-beam 3×3×5 t=0.5",
        "SA=65.5, V=13.75",
        &extrude(&polygon_profile(&l_pts), 0, DVec3::Z, 5.0).unwrap(),
    );

    // ── T-profile beam ────────────────────────────────────────────────
    // T: flange 4 wide × 0.5 thick; web 0.5 wide × 3.5 tall
    let t_pts = [
        DVec3::new(-2.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(2.0, 0.5, 0.0),
        DVec3::new(0.25, 0.5, 0.0),
        DVec3::new(0.25, 4.0, 0.0),
        DVec3::new(-0.25, 4.0, 0.0),
        DVec3::new(-0.25, 0.5, 0.0),
        DVec3::new(-2.0, 0.5, 0.0),
    ];
    add(
        &mut rows,
        "T-beam 4×4×6 t=0.5",
        "cross≈3.5",
        &extrude(&polygon_profile(&t_pts), 0, DVec3::Z, 6.0).unwrap(),
    );

    // ── I-beam profile ────────────────────────────────────────────────
    // I: two flanges 3 wide × 0.4 thick, web 0.3 wide × 2.2 tall, total h=3
    let hw = 1.5_f64; // half-flange width
    let ft = 0.4_f64; // flange thickness
    let wh = 0.15_f64; // half-web width
    let tot = 3.0_f64; // total height
    let i_pts = [
        DVec3::new(-hw, 0.0, 0.0),
        DVec3::new(hw, 0.0, 0.0),
        DVec3::new(hw, ft, 0.0),
        DVec3::new(wh, ft, 0.0),
        DVec3::new(wh, tot - ft, 0.0),
        DVec3::new(hw, tot - ft, 0.0),
        DVec3::new(hw, tot, 0.0),
        DVec3::new(-hw, tot, 0.0),
        DVec3::new(-hw, tot - ft, 0.0),
        DVec3::new(-wh, tot - ft, 0.0),
        DVec3::new(-wh, ft, 0.0),
        DVec3::new(-hw, ft, 0.0),
    ];
    add(
        &mut rows,
        "I-beam h=3 len=8",
        "sym. section",
        &extrude(&polygon_profile(&i_pts), 0, DVec3::Z, 8.0).unwrap(),
    );

    // ── Regular polygon prisms ─────────────────────────────────────────
    for (n, label) in [
        (5, "Pentagon prism r=1 h=2"),
        (8, "Octagon prism r=1 h=2"),
        (12, "Dodecagon prism r=1 h=2"),
    ] {
        let pts: Vec<DVec3> = (0..n)
            .map(|i| {
                let a = std::f64::consts::TAU * i as f64 / n as f64;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let prof = polygon_profile(&pts);
        add(
            &mut rows,
            label,
            "",
            &extrude(&prof, 0, DVec3::Z, 2.0).unwrap(),
        );
    }

    // ── Star / hexagram profile ────────────────────────────────────────
    let r_out = 1.5_f64;
    let r_in = 0.75_f64;
    let star_pts: Vec<DVec3> = (0..12)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / 12.0;
            let r = if i % 2 == 0 { r_out } else { r_in };
            DVec3::new(r * a.cos(), r * a.sin(), 0.0)
        })
        .collect();
    add(
        &mut rows,
        "Star (12-pt) h=1",
        "",
        &extrude(&polygon_profile(&star_pts), 0, DVec3::Z, 1.0).unwrap(),
    );

    // ── Print table ────────────────────────────────────────────────────
    let w = 92;
    println!(
        "{:<28}  {:>10}  {:>10}  {:>26}  {:<16}  {}",
        "Shape", "Area", "Volume", "Centroid (x,y,z)", "Expected", "Valid?"
    );
    println!("{}", "─".repeat(w));
    for r in &rows {
        println!(
            "{:<28}  {:>10.4}  {:>10.4}  ({:>6.3},{:>6.3},{:>6.3})  {:<16}  {}",
            r.name,
            r.area,
            r.vol,
            r.c.x,
            r.c.y,
            r.c.z,
            r.expected,
            if r.valid { "✓" } else { "✗" }
        );
    }
    println!("{}", "─".repeat(w));
    println!(
        "{} shapes. All geometry computed from triangulated faces.",
        rows.len()
    );
}

struct Row {
    name: &'static str,
    area: f64,
    vol: f64,
    c: DVec3,
    expected: &'static str,
    valid: bool,
}

fn add(rows: &mut Vec<Row>, name: &'static str, expected: &'static str, brep: &BRep) {
    rows.push(Row {
        name,
        area: surface_area(brep),
        vol: volume(brep),
        c: centroid(brep),
        expected,
        valid: check(brep).is_valid(),
    });
}

fn polygon_profile(pts: &[DVec3]) -> BRep {
    let n = pts.len();
    let mut brep = BRep::default();
    let vis: Vec<usize> = pts.iter().map(|&p| make_vertex(&mut brep, p)).collect();
    let mut wires = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = pts[i];
        let b = pts[j];
        let dir = (b - a).normalize_or_zero();
        let len = (b - a).length();
        let eidx = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            vis[i],
            vis[j],
        )
        .unwrap();
        wires.push(WireEdge::fwd(eidx));
    }
    let surface = Surface3::Plane(Plane {
        origin: pts[0],
        normal: DVec3::Z,
    });
    make_face(&mut brep, surface, make_wire(wires), vec![]).unwrap();
    brep
}
