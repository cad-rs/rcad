//! Debug: sphere-box union — raw builder vs post-processed
use glam::DVec3;
use rcad_algorithms::geom_populate::populate_box_geom;
use rcad_algorithms::bopds;
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_algorithms::builder::BooleanBuilder;
use rcad_algorithms::brep_repair::merge_close_vertices;
use rcad_algorithms::{
    boolean_op, boolean_op_with_history, BooleanOpType,
    total_surface_area, total_volume,
};
use rcad_algorithms::tolerance::TOLERANCE_ABS;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_modeling::{make_box_brep, make_sphere_brep};

fn face_count(brep: &BRep) -> usize {
    brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count()
}

fn describe(name: &str, brep: &BRep) {
    let sa = total_surface_area(brep);
    let vol = total_volume(brep);
    let nf = face_count(brep);
    let ns = brep.solids.iter().map(|s| s.shells.len()).sum::<usize>();
    let nso = brep.solids.len();
    eprintln!("{name}: SA={sa:.6} vol={vol:.6} faces={nf} shells={ns} solids={nso}");
}

#[test]
fn debug_box_sphere_raw_builder_vs_postproc() {
    let s = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere");
    let mut b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    populate_box_geom(&mut b);

    describe("sphere", &s);
    describe("box", &b);
    eprintln!();

    // === PATH 1: Raw builder (PaveFiller + BooleanBuilder, no post-processing) ===
    let mut ds = bopds::ds::DS::new(&s, &b);
    let bvh_s = rcad_algorithms::bvh::Bvh::build(&s);
    let bvh_b = rcad_algorithms::bvh::Bvh::build(&b);
    let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_s, &bvh_b);
    filler.perform();
    let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
    let raw = builder.build().expect("raw union build");
    describe("1  raw builder (no postproc)", &raw);

    // DEBUG: analyze shared boundary vertices
    {
        eprintln!("\n=== Face analysis for PATH 1 ===");
        let brep = &raw;
        for (fi, face) in brep.solids[0].shells[0].faces.iter().enumerate() {
            let n_edges = face.outer_wire.edges.len();
            let edge_verts: Vec<(usize, usize)> = face.outer_wire.edges.iter().map(|we| {
                let e = &brep.edges[we.idx];
                (e.start, e.end)
            }).collect();
            let vert_str: Vec<String> = edge_verts.iter().flat_map(|(s,e)| vec![*s, *e]).take(12).map(|vi| {
                let p = brep.vertices[vi].point;
                format!("({:.4},{:.4},{:.4})", p.x, p.y, p.z)
            }).collect();
            let n_tris = face.triangles.len();
            eprintln!("  face[{fi}]: edges={} tris={} first_verts: {:?}", n_edges, n_tris, vert_str);
        }
    }

    // === PATH 2: raw builder + merge_close_vertices ===
    let (sewn, _n) = merge_close_vertices(&raw, TOLERANCE_ABS * 64.0);
    describe("2  +merge_close_vertices", &sewn);

    // === PATH 3: boolean_op_with_history (fuse_with_history — no postproc) ===
    let (hist_result, _history) = boolean_op_with_history(BooleanOpType::Union, &s, &b)
        .expect("boolean_op_with_history");
    describe("3  boolean_op_with_history", &hist_result);

    // === PATH 4: boolean_op (fuse = full postproc) ===
    let fused = boolean_op(BooleanOpType::Union, &s, &b).expect("fuse");
    describe("4  boolean_op (fuse, full postproc)", &fused);

    // === Box-box union (known working baseline) ===
    {
        let mut b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
        let mut b2 = make_box_brep(DVec3::splat(0.5), DVec3::X, DVec3::Y, 0.5, 0.5, 0.5).expect("box2");
        populate_box_geom(&mut b1);
        populate_box_geom(&mut b2);
        let mut ds_bb = bopds::ds::DS::new(&b1, &b2);
        let bvh_bb1 = rcad_algorithms::bvh::Bvh::build(&b1);
        let bvh_bb2 = rcad_algorithms::bvh::Bvh::build(&b2);
        let mut filler_bb = PaveFiller::with_bvh(&mut ds_bb, &bvh_bb1, &bvh_bb2);
        filler_bb.perform();
        let builder_bb = BooleanBuilder::new(&ds_bb, BooleanOpType::Union);
        let res_bb = builder_bb.build().expect("box-box union");
        describe("5  box-box union (baseline)", &res_bb);
    }

    // === Sphere fully inside box (no intersection expected) ===
    {
        let s_in = make_sphere_brep(DVec3::splat(0.5), 0.4).expect("sphere inside");
        let mut b_out = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
        populate_box_geom(&mut b_out);
        // Union should be just the box
        let mut ds = bopds::ds::DS::new(&s_in, &b_out);
        let bvh_s = rcad_algorithms::bvh::Bvh::build(&s_in);
        let bvh_b = rcad_algorithms::bvh::Bvh::build(&b_out);
        let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_s, &bvh_b);
        filler.perform();
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let res = builder.build().expect("contained sphere union");
        describe("6  sphere-in-box union (contained)", &res);
        // Intersection should be just the sphere
        let mut ds2 = bopds::ds::DS::new(&s_in, &b_out);
        let mut filler2 = PaveFiller::with_bvh(&mut ds2, &bvh_s, &bvh_b);
        filler2.perform();
        let builder2 = BooleanBuilder::new(&ds2, BooleanOpType::Intersection);
        let res2 = builder2.build().expect("contained sphere intersection");
        describe("7  sphere-in-box intersection (contained)", &res2);
    }

    // === Sphere-box with NO face intersection (sphere corner intersection only) ===
    {
        // Sphere touches box corner, no face intersection
        let s_corner = make_sphere_brep(DVec3::splat(1.0), 0.5).expect("sphere corner");
        let mut b_corner = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
        populate_box_geom(&mut b_corner);
        let mut ds = bopds::ds::DS::new(&s_corner, &b_corner);
        let bvh_s = rcad_algorithms::bvh::Bvh::build(&s_corner);
        let bvh_b = rcad_algorithms::bvh::Bvh::build(&b_corner);
        let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_s, &bvh_b);
        filler.perform();
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let res = builder.build().expect("corner sphere union");
        describe("8  sphere-touching-corner union", &res);
    }

    eprintln!();
    eprintln!("Expected (A1) SA ~14.6394, vol ~4.665");
}
