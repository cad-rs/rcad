// Sanity test for the OCCT-aligned translation of
// BOPTools_AlgoTools::AreFacesSameDomain (BOPTools_AlgoTools.cxx L1139-1205)
// and its dependency BOPTools_AlgoTools3D::PointInFace (L906-988).
//
// Verifies the full chain on box faces: PointInFace finds a point strictly
// inside the first face, and IsValidPointForFace classifies it on the second
// face.

use glam::DVec3;
use rcad_algo::bop::algo::pave_filler::PaveFiller;
use rcad_algo::bop::int_tools::context::IntToolsContext;
use rcad_algo::bop::tools::algo_tools::are_faces_same_domain;
use rcad_algo::topalgo::shape_source::ShapeSource;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topods::{self, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_modeling::prim::primapi::make_box_brep;

fn root_shape(brep: &topods::BRep, location: u32) -> Shape {
    for (i, ts) in brep.tshapes.iter().enumerate().rev() {
        match &**ts {
            TShape::Solid(_) | TShape::Shell(_) => {
                return Shape::from_parts(ts.clone(), i, location, Orientation::Forward);
            }
            _ => {}
        }
    }
    panic!("no root Solid/Shell in BRep");
}

#[test]
fn are_faces_same_domain_box() {
    // Two identical unit boxes: A's faces are coplanar with B's.
    let a = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("box a");
    let b = make_box_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("box b");

    let mut filler = PaveFiller::new();
    filler.set_arguments(vec![root_shape(&a, 0), root_shape(&b, 1)]);
    filler.set_fuzzy_value(0.0);
    let prog = NoopProgress;
    let ps = ProgressScope::new(&prog, "test", 100);
    filler.perform(&ps);
    let ds = filler.ds();

    let mut faces: Vec<usize> = Vec::new();
    for i in 0..ds.nb_shapes() {
        if ds.shape_type(i) == ShapeType::Face {
            faces.push(i);
        }
    }
    // Box A (6) + box B (6) source faces.
    assert!(faces.len() >= 12, "expected >=12 faces, got {}", faces.len());

    let mut ctx = IntToolsContext::new();
    // A face is same-domain with itself.
    for &f in faces.iter().take(6) {
        let ok = are_faces_same_domain(f, f, &mut ctx, ds, 0.0);
        assert!(ok, "box face {} not same-domain with itself", f);
    }
    // Two perpendicular faces of box A are not same-domain.
    let ok = are_faces_same_domain(faces[0], faces[1], &mut ctx, ds, 0.0);
    assert!(!ok, "perpendicular box faces {} and {} reported SD", faces[0], faces[1]);
    // Coplanar faces of the two identical boxes are same-domain.
    for i in 0..6 {
        let ok = are_faces_same_domain(faces[i], faces[i + 6], &mut ctx, ds, 0.0);
        assert!(ok, "coplanar box faces {} and {} not same-domain", faces[i], faces[i + 6]);
    }
}
