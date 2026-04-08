use criterion::{criterion_group, criterion_main, Criterion};
use glam::DVec3;
use rcad_algorithms::{boolean_op, BooleanOpType};
use rcad_modeling::{make_box_brep, make_sphere_brep};

fn boolean_union_boxes(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let b = make_box_brep(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    c.bench_function("boolean_union_boxes", |bench| {
        bench.iter(|| boolean_op(BooleanOpType::Union, &a, &b).unwrap())
    });
}

fn boolean_diff_box_sphere(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 0.8).unwrap();
    c.bench_function("boolean_diff_box_sphere", |bench| {
        bench.iter(|| boolean_op(BooleanOpType::Difference, &a, &b).unwrap())
    });
}

fn intss_plane_sphere(c: &mut Criterion) {
    use rcad_algorithms::inttools::intss::intersect_surfaces;
    use rcad_kernel::geom::{Plane, SphericalSurface, Surface3};
    let s1 = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
    let s2 = Surface3::Sphere(SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0 });
    c.bench_function("intss_plane_sphere", |bench| {
        bench.iter(|| intersect_surfaces(&s1, &s2))
    });
}

criterion_group!(benches, boolean_union_boxes, boolean_diff_box_sphere, intss_plane_sphere);
criterion_main!(benches);
