// Quick topology checker for bfuse_simple B1-B9
// Run with: cargo test --test bfuse_topo_check -- --nocapture  (from rcad/ dir)
//
// OCCT Draw script pattern:
//   box b1 <x> <y> <z> <dx> <dy> <dz>
//   nurbsconvert b1 b1        ← b1 = BSpline box
//   box b2 <x> <y> <z> <dx> <dy> <dz>   ← b2 = Plane box
//   bfuse result b2 b1        ← b2=ShapeA, b1=ShapeB

use rcad_algorithms::{boolean_op, BooleanOpType};
use rcad_algorithms::geom_convert::{surface_to_bspline, ConvertParams};
use rcad_kernel::{topods, Surface3, BRep};
use rcad_modeling::make_box_brep;

/// Convert all surfaces in a BRep to BSpline (matching OCCT `nurbsconvert`).
fn nurbsconvert(mut brep: BRep) -> BRep {
    let params = ConvertParams::default();
    brep.geom.surfaces = brep.geom.surfaces.into_iter()
        .map(|s| Surface3::BSpline(surface_to_bspline(&s, &params))).collect();
    brep
}

fn box_brep(x: f64, y: f64, z: f64, dx: f64, dy: f64, dz: f64) -> BRep {
    make_box_brep(glam::DVec3::new(x, y, z), glam::DVec3::X, glam::DVec3::Y, dx, dy, dz).expect("box")
}

fn count_topology(brep: &BRep) -> (usize, usize, usize) {
    use std::collections::BTreeSet;
    // OCCT checknbshapes: count only vertices referenced by edges (not triangulation vertices)
    let v = brep.edges.iter().flat_map(|e| [e.start, e.end]).collect::<BTreeSet<_>>().len();
    let e = brep.edges.len();
    let f: usize = brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    (v, e, f)
}

fn count_surface_kinds(brep: &BRep) -> (usize, usize) {
    let mut np = 0; let mut nb = 0;
    for (fi, _) in brep.solids.iter().flat_map(|s|&s.shells).flat_map(|sh|&sh.faces).enumerate() {
        match brep.geom.face_surface.get(fi).copied().flatten().and_then(|si| brep.geom.surfaces.get(si)) {
            Some(Surface3::Plane(_)) => np += 1,
            Some(Surface3::BSpline(_)) => nb += 1,
            _ => {}
        }
    }
    (np, nb)
}

// OCCT reference topology for bfuse_simple B1-B9
// From tests/occt/step_reference/occt_boolean_bfuse_simple_b*.json
const REF: &[(&str, usize, usize, usize)] = &[
    ("B1", 14, 22, 10),
    ("B2", 14, 23, 11),
    ("B3", 16, 28, 14),
    ("B4", 16, 28, 14),
    ("B5", 16, 28, 14),
    ("B6", 14, 21,  9),
    ("B7", 15, 24, 11),
    ("B8", 15, 24, 12),
    ("B9", 15, 24, 12),
];

/// Run union matching OCCT's `bfuse result <plane_box> <bspline_box>`.
/// plane = ShapeA (b2 in OCCT draw), bspline = ShapeB (b1 in OCCT, nurbsconverted).
fn run_case(name: &str, ref_v: usize, ref_e: usize, ref_f: usize,
            plane: BRep, bspline: BRep, label: &str) {
    let r = boolean_op(BooleanOpType::Union, &plane, &bspline).expect(label);
    let (v, e, f) = count_topology(&r);
    let (np, nb) = count_surface_kinds(&r);
    let v_ok = if v == ref_v { "✅" } else { "❌" };
    let e_ok = if e == ref_e { "✅" } else { "❌" };
    let f_ok = if f == ref_f { "✅" } else { "❌" };
    println!("{name}: V={v} E={e} F={f}   {np}P+{nb}BS   vs  V={ref_v} E={ref_e} F={ref_f}  {v_ok}V {e_ok}E {f_ok}F");
}

// OCCT script for each case:
//   box b1 <bx1> <by1> <bz1> <bdx1> <bdy1> <bdz1>
//   nurbsconvert b1 b1
//   box b2 <bx2> <by2> <bz2> <bdx2> <bdy2> <bdz2>
//   bfuse result b2 b1

#[test]
fn bfuse_b1_box_on_corner() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 0 0 0.5 1 0.5
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,0.,0., 0.5,1.,0.5);
    run_case("B1", 14,22,10, plane, bspline, "fuse B1");
}

#[test]
fn bfuse_b2_box_side_overhang() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 -0.5 0 0.5 0.5 1
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,-0.5,0., 0.5,0.5,1.);
    run_case("B2", 14,23,11, plane, bspline, "fuse B2");
}

#[test]
fn bfuse_b3_box_side_overhang_taller() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 -0.5 0 0.5 1.5 1
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,-0.5,0., 0.5,1.5,1.);
    run_case("B3", 16,28,14, plane, bspline, "fuse B3");
}

#[test]
fn bfuse_b4_box_side_overlap() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 0.5 0 1 1 1
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,0.5,0., 1.,1.,1.);
    run_case("B4", 16,28,14, plane, bspline, "fuse B4");
}

#[test]
fn bfuse_b5_box_insert() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 0.25 0 1 0.5 1
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,0.25,0., 1.,0.5,1.);
    run_case("B5", 16,28,14, plane, bspline, "fuse B5");
}

#[test]
fn bfuse_b6_box_on_bottom() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 0 0 1 1 0.5
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,0.,0., 1.,1.,0.5);
    run_case("B6", 14,21,9, plane, bspline, "fuse B6");
}

#[test]
fn bfuse_b7_box_insert_thin() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 0 0 1 0.3 1
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,0.,0., 1.,0.3,1.);
    run_case("B7", 15,24,11, plane, bspline, "fuse B7");
}

#[test]
fn bfuse_b8_box_insert_thin_half() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): 0 0 0 1 0.3 0.5
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(0.,0.,0., 1.,0.3,0.5);
    run_case("B8", 15,24,12, plane, bspline, "fuse B8");
}

#[test]
fn bfuse_b9_box_wrap() {
    // b1(BSpline): 0 0 0 1 1 1    b2(Plane): -0.5 0 -0.5 2 1 2
    let bspline = nurbsconvert(box_brep(0.,0.,0., 1.,1.,1.));
    let plane = box_brep(-0.5,0.,-0.5, 2.,1.,2.);
    run_case("B9", 15,24,12, plane, bspline, "fuse B9");
}
