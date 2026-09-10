//! TEMP probe: per-pair FF curve counts vs OCCT ref (delete after use).

use glam::DVec3;
use rcad_algo::bop::algo::pave_filler::PaveFiller;
use rcad_algo::bop::ds::DS;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_modeling::prim::primapi::make_cylinder_brep;

fn root_shape(brep: &rcad_kernel::topods::BRep, location: u32) -> Shape {
    for (i, ts) in brep.tshapes.iter().enumerate().rev() {
        match &**ts {
            TShape::Solid(_) | TShape::Shell(_) => {
                return Shape::from_parts(ts.clone(), i, location, Orientation::Forward);
            }
            _ => {}
        }
    }
    panic!("no root");
}

fn ff_pairs(ds: &DS) -> Vec<(usize, usize, usize)> {
    let mut out = Vec::new();
    for ff in &ds.interf_ff {
        let (a, b) = (ff.f1, ff.f2);
        let n = ff.curves.len();
        out.push((a, b, n));
    }
    out.sort();
    out
}

#[test]
fn probe_ff_pair_curves_rotated_cyl() {
    let a = make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 2.0)
        .expect("cylinder a");
    let rad_b = 45.0_f64.to_radians();
    let b = make_cylinder_brep(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::Z,
        DVec3::X * rad_b.cos() + DVec3::Y * rad_b.sin(),
        1.0,
        2.0,
    )
    .expect("cylinder b");

    let mut filler = PaveFiller::new();
    filler.set_arguments(vec![root_shape(&a, 0), root_shape(&b, 1)]);
    filler.stop_after = Some("after_PerformFF");
    let prog = NoopProgress;
    let ps = ProgressScope::new(&prog, "probe", 100);
    filler.perform(&ps);

    let ds = filler.ds();
    let n_ic = ds.intersection_curves.len();
    let pairs = ff_pairs(ds);
    eprintln!("PROBE nIC={n_ic}");
    for (a, b, n) in &pairs {
        eprintln!("PROBE pair ({a},{b}) curves={n}");
    }
    for ff in &ds.interf_ff {
        for (k, &ci) in ff.curves.iter().enumerate() {
            let c = &ds.intersection_curves[ci];
            match &c.curve {
                rcad_kernel::geom::Curve3::Circle(circ) => {
                    eprintln!(
                        "PROBE curve ff=({},{}) #{k} range={:?} CIRCLE center=({:.9},{:.9},{:.9}) normal=({:.9},{:.9},{:.9}) x_dir=({:.9},{:.9},{:.9}) r={:.9}",
                        ff.f1, ff.f2, c.t_range,
                        circ.center.x, circ.center.y, circ.center.z,
                        circ.normal.x, circ.normal.y, circ.normal.z,
                        circ.x_dir.x, circ.x_dir.y, circ.x_dir.z,
                        circ.radius
                    );
                    // params of key points: rim crossing (0.7071,0.7071,2),
                    // seam point (1,0,2), (0,1,2)
                    for (name, p) in [
                        ("rim45", glam::DVec3::new(0.7071067811865476, 0.7071067811865476, 2.0)),
                        ("seamX", glam::DVec3::new(1.0, 0.0, 2.0)),
                        ("seamY", glam::DVec3::new(0.0, 1.0, 2.0)),
                    ] {
                        let rel = p - circ.center;
                        let u = rel.dot(circ.x_dir);
                        let v = rel.dot(circ.y_dir);
                        let par = v.atan2(u);
                        eprintln!("PROBE   param({name}) = {par:.9}  (dist to circle plane = {:.3e})", rel.dot(circ.normal));
                    }
                }
                _ => {
                    eprintln!("PROBE curve ff=({},{}) #{k} range={:?} type={:?}", ff.f1, ff.f2, c.t_range, std::mem::discriminant(&c.curve));
                }
            }
        }
    }
    // OCCT ref: (2,22)=1 (2,24)=1 (9,15)=1 (11,15)=1; others 0; nIC=4.
    assert_eq!(n_ic, 4, "PROBE nIC: rcad={n_ic} ref=4");
}
