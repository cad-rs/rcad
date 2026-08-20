//! OCCT BRepGProp_Gauss::Compute(Face, Domain) alignment.
//!
//! OCCT BRepGProp::SurfaceProperties (BRepGProp.cxx L367-395) integrates every
//! non-natural-restriction face (NbChildren() != 0, i.e. a face with wires)
//! via BRepGProp_Gauss::Compute(Face, Domain) (BRepGProp_Gauss.cxx L1126-1211):
//! the Green-theorem line integral of |Su x Sv| over the face domain bounded
//! by its edge pcurves.  rcad's `face_surface_area_gauss_domain` is the 1:1
//! translation.
//!
//! Reference: OCCT `pcone 1 0.5 2` lateral face (r1=1, r2=0.5, h=2):
//!   slant l = sqrt(h^2 + (r1-r2)^2) = sqrt(4.25) = 2.0615528128...
//!   area = PI * (r1+r2) * l = PI * 1.5 * 2.0615528128 = 9.7139956...

use rcad_kernel::base::gprop::surface::face_surface_area_gauss_domain;
use rcad_kernel::base::gprop::tri::face_flat_iter;
use rcad_kernel::topods;
use glam::DVec3;

#[test]
fn cone_lateral_gauss_domain_matches_occt() {
    // make_cone_brep(center, axis, ref_dir, base_radius, top_radius, height):
    // base circle at z=0 (make_cone.rs "base at z=0"), matching OCCT
    // BRepPrimAPI_MakeCone(R1, R2, H) with gp::XOY().
    let cone = rcad_modeling::make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 0.5, 2.0).unwrap();
    let faces = face_flat_iter(&cone);
    let mut lateral: Option<f64> = None;
    let mut lateral_idx = None;
    for (ti, f) in &faces {
        let is_lateral = matches!(&*cone.tshapes[*ti], topods::TShape::Face(fd) if matches!(&fd.surface, Some(rcad_kernel::geom::Surface3::Cone(_))));
        if is_lateral {
            lateral = Some(face_surface_area_gauss_domain(&cone, f, *ti));
            lateral_idx = Some(*ti);
        }
    }
    let expected = std::f64::consts::PI * 1.5 * (4.0f64 + 0.25).sqrt();
    let got = lateral.expect("gauss_domain on cone lateral face returned None");
    assert!(
        (got - expected).abs() < 1e-6,
        "cone lateral gauss_domain = {}, expected {} (face idx {:?})",
        got, expected, lateral_idx
    );
}
