// P013 pattern: two same-axis cylinders (r=1, h=2), the second rotated 45 deg
// about the shared Z axis.  OCCT reference: V=4 E=6 F=4 S=1 for the common.
use glam::DVec3;
use rcad_algo::bop::brep_algo_api::{common, fuse};
use rcad_kernel::topods::{self, TShape};
use rcad_modeling::make_cylinder_brep;

fn make_pair() -> (topods::BRep, topods::BRep) {
    let a = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).expect("cyl a");
    let rad = 45.0_f64.to_radians();
    let b = make_cylinder_brep(
        DVec3::ZERO,
        DVec3::Z,
        DVec3::X * rad.cos() + DVec3::Y * rad.sin(),
        1.0,
        2.0,
    )
    .expect("cyl b");
    (a, b)
}

fn count_topology(brep: &topods::BRep) -> (usize, usize, usize, usize) {
    let mut n_v = 0usize;
    let mut n_e = 0usize;
    let mut n_f = 0usize;
    let mut n_s = 0usize;
    for ts in &brep.tshapes {
        match &**ts {
            TShape::Vertex(_) => n_v += 1,
            TShape::Edge(_) => n_e += 1,
            TShape::Face(_) => n_f += 1,
            TShape::Solid(_) => n_s += 1,
            _ => {}
        }
    }
    (n_v, n_e, n_f, n_s)
}

#[test]
fn p013_cyl_cyl_rotated_common() {
    let (a, b) = make_pair();
    let result = common(&a, &b).expect("common");
    let (v, e, f, s) = count_topology(&result);
    eprintln!("P013 common: V={v} E={e} F={f} S={s}");
    // OCCT reference (same-axis cylinders rotated 45 deg): V=4 E=6 F=4 S=1.
    assert_eq!((v, e, f, s), (4, 6, 4, 1), "P013 common topology");
}

#[test]
fn p013_cyl_cyl_rotated_fuse() {
    let (a, b) = make_pair();
    let result = fuse(&a, &b).expect("fuse");
    let (v, e, f, s) = count_topology(&result);
    eprintln!("P013 fuse: V={v} E={e} F={f} S={s}");
    assert_eq!((v, e, f, s), (4, 6, 4, 1), "P013 fuse topology");
}
