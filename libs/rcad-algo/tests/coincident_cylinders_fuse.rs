// Regression test for bfuse_simple L2: two coaxial same-radius cylinders
// (pcylinder(20,100) | pcylinder(20,100) translated +50 Z) fused.
// OCCT reference topology: V=4 E=7 F=5 S=1. Guards the CurveOnClosedSurface
// rep on the shared seam edge (a coincident-edge CS rep for the second face).
use glam::DVec3;
use rcad_algo::bop::brep_algo_api::fuse;
use rcad_kernel::topods::{self, TShape};
use rcad_modeling::make_cylinder_brep;

fn make_cyl(z0: f64) -> topods::BRep {
    make_cylinder_brep(DVec3::new(0.0, 0.0, z0), DVec3::Z, DVec3::X, 20.0, 100.0).unwrap()
}

#[test]
fn coincident_cylinders_fuse() {
    let c1 = make_cyl(0.0);
    let c2 = make_cyl(50.0);
    let result = fuse(&c1, &c2).expect("fuse");
    let solids = result.solids();
    let geom = result.geom();
    let mut ekeys: Vec<(usize, u32)> = Vec::new();
    let mut n_f = 0usize;
    let mut vpts: Vec<DVec3> = Vec::new();
    for s in &solids {
        for sh in &s.shells {
            for f in &sh.faces {
                n_f += 1;
                let _st = f.surface_idx.and_then(|i| geom.surfaces.get(i)).map(|sf| match sf {
                    rcad_kernel::geom::Surface3::Plane(p) => format!("PL(z={:.1})", p.origin.z),
                    rcad_kernel::geom::Surface3::Cylinder(cy) => format!("CY(r={:.1})", cy.radius),
                    _ => "??".into(),
                }).unwrap_or_else(|| "?".into());
                for w in std::iter::once(&f.outer_wire).chain(f.inner_wires.iter()) {
                    for e in &w.edges {
                        if !ekeys.iter().any(|&(i, l)| i == e.idx && l == e.location) {
                            ekeys.push((e.idx, e.location));
                        }
                        if let Some(ts) = result.tshapes.get(e.idx) {
                            if let TShape::Edge(ed) = &**ts {
                                for p in [&ed.first, &ed.last] {
                                    if let Some(v) = p.as_vertex() {
                                        if !vpts.iter().any(|q| q.distance(v.point) < 1e-7) {
                                            vpts.push(v.point);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    let n_v = vpts.len();
    let n_e = ekeys.len();
    // OCCT reference for L2: V=4 E=7 F=5 S=1.
    assert_eq!(n_v, 4, "V must match OCCT reference");
    assert_eq!(n_e, 7, "E must match OCCT reference");
    assert_eq!(n_f, 5, "F must match OCCT reference");
    assert_eq!(solids.len(), 1, "S must match OCCT reference");
}
