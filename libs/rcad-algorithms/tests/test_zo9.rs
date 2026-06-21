use glam::{DAffine3, DVec3};
use rcad_kernel::Curve2dEval;
use rcad_modeling::{make_box_brep, make_cylinder_brep};
use rcad_algorithms::bopds::ds::{DS, ShapeOrigin};
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_algorithms::builder::BooleanBuilder;
use rcad_algorithms::{boolean_op, BooleanOpType, total_surface_area};

fn surf_type_str(surface: &rcad_kernel::geom::Surface3) -> &'static str {
    match surface {
        rcad_kernel::geom::Surface3::Plane(_) => "Plane",
        rcad_kernel::geom::Surface3::Cylinder(_) => "Cylinder",
        rcad_kernel::geom::Surface3::Sphere(_) => "Sphere",
        rcad_kernel::geom::Surface3::Cone(_) => "Cone",
        _ => "Other",
    }
}

fn curve_type_str(curve: &rcad_kernel::geom::Curve3) -> &'static str {
    match curve {
        rcad_kernel::geom::Curve3::Line(_) => "Line",
        rcad_kernel::geom::Curve3::Circle(_) => "Circle",
        rcad_kernel::geom::Curve3::Ellipse(_) => "Ellipse",
        _ => "Other",
    }
}

fn analyze_ds(ds: &DS) {
    println!("  DS vertices: {}", ds.vertices.len());
    println!("  DS edges: {}", ds.edges.len());
    println!("  DS faces: {}", ds.faces.len());
    println!("  Intersection curves: {}", ds.intersection_curves.len());
    println!("  Interferences: {}", ds.interferences.len());
    for (i, ic) in ds.intersection_curves.iter().enumerate() {
        let ct = curve_type_str(&ic.curve);
        println!("    IC[{i}]: {ct} t=[{:.3},{:.3}] polyline_len={}", ic.t_range[0], ic.t_range[1], ic.polyline.len());
    }
    for (fi, face) in ds.faces.iter().enumerate() {
        let origin = match face.origin {
            ShapeOrigin::ShapeA => "A",
            ShapeOrigin::ShapeB => "B",
        };
        let surf = surf_type_str(&face.surface);
        println!("  Face[{fi}] ({origin}) {surf}: curves_sc={}, verts_in={}, bnd_verts={}",
            face.face_info.curves_sc.len(),
            face.face_info.vertices_in.len(),
            face.boundary_verts.len());
        for &ci in &face.face_info.curves_sc {
            let ic = &ds.intersection_curves[ci];
            println!("      curves_sc: IC[{ci}] {} polyline_len={}", curve_type_str(&ic.curve), ic.polyline.len());
        }
    }
}

/// Reproduce OCCT bopfuse_simple ZO9:
/// box b1 4 4 4
/// pcylinder b2 1 4
/// ttranslate b2 3 1 4
/// bop b1 b2; bopfuse result
/// checkprops result -s 121.133
#[test]
fn debug_zo9() {
    let b1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();

    // pcylinder b2 1 4 — cylinder centered at (0,0,2), radius 1, height 4, axis Z
    let b2_base = make_cylinder_brep(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, DVec3::X, 1.0, 4.0).unwrap();
    let mut b2 = b2_base;
    // ttranslate b2 3 1 4
    b2.apply_transform(DAffine3::from_translation(DVec3::new(3.0, 1.0, 4.0)));

    println!("=== Input box ===");
    println!("  Vertices: {}  Edges: {}", b1.vertices.len(), b1.edges.len());

    println!("=== Input cylinder ===");
    println!("  Vertices: {}  Edges: {}", b2.vertices.len(), b2.edges.len());
    for (i, v) in b2.vertices.iter().enumerate() {
        println!("  Vertex[{}]: ({:.3}, {:.3}, {:.3})", i, v.point.x, v.point.y, v.point.z);
    }

    // DS + pave_fill diagnostics
    let mut ds = DS::new(&b1, &b2);
    let bvh_a = rcad_algorithms::bvh::Bvh::build(&b1);
    let bvh_b = rcad_algorithms::bvh::Bvh::build(&b2);
    let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_a, &bvh_b);
    filler.perform();

    println!("\n=== DS after pave_fill ===");
    analyze_ds(&ds);

    // Diagnose cylinder wall face UV boundary
    for (fi, face) in ds.faces.iter().enumerate() {
        if let rcad_kernel::geom::Surface3::Cylinder(_) = &face.surface {
            println!("  Cylinder Face[{fi}] UV boundary:");
            if let Some(ref uv_bnd) = face.uv_boundary {
                let u_min = uv_bnd.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                let u_max = uv_bnd.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                let v_min = uv_bnd.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                let v_max = uv_bnd.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                println!("    U range: [{:.6}, {:.6}] span={:.6}", u_min, u_max, u_max - u_min);
                println!("    V range: [{:.6}, {:.6}] span={:.6}", v_min, v_max, v_max - v_min);
                println!("    Boundary points ({})", uv_bnd.len());
            }
            // Print pcurve samples for each intersection curve
            for &ci in &face.face_info.curves_sc {
                let ic = &ds.intersection_curves[ci];
                let [t0, t1] = ic.t_range;
                // Cylinder is ShapeB, use pcurve_on_b
                if let Some(ref pcurve) = ic.pcurve_on_b {
                    println!("    IC[{ci}] (pcurve_on_b): t=[{:.6},{:.6}]", t0, t1);
                    for frac in [0.0, 0.25, 0.5, 0.75, 1.0] {
                        let t = t0 + (t1 - t0) * frac;
                        let uv = pcurve.point_at(t);
                        println!("      t={:.2}: u={:.6} v={:.6}", frac, uv.x, uv.y);
                    }
                }
                if let Some(ref pcurve) = ic.pcurve_on_a {
                    let uv = pcurve.point_at(t0);
                    println!("    IC[{ci}] pcurve_on_a t0: u={:.6} v={:.6}", uv.x, uv.y);
                }
            }
        }
    }
    // Diagnose A-face surfaces to figure out which is box top vs bottom
    for (fi, face) in ds.faces.iter().enumerate() {
        if face.origin == ShapeOrigin::ShapeA {
            if let rcad_kernel::geom::Surface3::Plane(p) = &face.surface {
                let d = p.normal.dot(p.origin);
                println!("  Face[{fi}] (A) Plane: d={:.6}, normal=({:.2},{:.2},{:.2}) origin=({:.2},{:.2},{:.2}) curves_sc={}",
                    d, p.normal.x, p.normal.y, p.normal.z,
                    p.origin.x, p.origin.y, p.origin.z,
                    face.face_info.curves_sc.len());
            }
        }
    }

    // Check which split path the builder uses
    let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
    let result = builder.build().unwrap();
    let sa_raw = total_surface_area(&result);
    let n_vert = result.vertices.len();
    let n_edges = result.edges.len();
    let n_faces: usize = result.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();
    println!("\n=== Builder result (raw, before fuse post-processing) ===");
    println!("Vertices: {}  Edges: {}  Faces: {}", n_vert, n_edges, n_faces);
    println!("Surface area: {:.3} (expected 121.133)", sa_raw);

    // Now try the full fuse path
    let r = boolean_op(BooleanOpType::Union, &b1, &b2);
    match r {
        Ok(ref result) => {
            let sa = total_surface_area(result);
            let n_vert = result.vertices.len();
            let n_edges = result.edges.len();
            let n_faces: usize = result.solids.iter()
                .flat_map(|s| &s.shells)
                .flat_map(|sh| &sh.faces)
                .count();
            println!("\n=== Full fuse result ===");
            println!("Vertices: {}  Edges: {}  Faces: {}", n_vert, n_edges, n_faces);
            println!("Surface area: {:.3} (expected 121.133)", sa);
            println!("Error: {:.3} ({:.1}%)", sa - 121.133, (sa - 121.133) / 121.133 * 100.0);

            if (sa - 121.133).abs() > 0.15 * 121.133 {
                println!("FAIL: SA mismatch exceeds 15% tolerance");
            } else {
                println!("PASS: SA within tolerance");
            }
        }
        Err(e) => println!("Full boolean_op failed: {:?}", e),
    }
}
