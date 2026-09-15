//! A0 discriminating tests: the Poly data layer, the BRep mount semantics and
//! the triangulation persistence record set.
//!
//! Coverage:
//! 1. `plane_face_triangulation_round_trip_is_lossless` — a programmatically
//!    triangulated plane face is attached to the BRep, serialized through the
//!    OCCT v2/v3-aligned record set and read back; every node / triangle /
//!    UV / normal value compares EXACTLY (bit-equal f64 / f32), including the
//!    deflection, the mesh purpose and the triangulation parameters.
//! 2. `cylinder_seam_face_carries_double_polygon_on_triangulation` — the seam
//!    edge of a cylindrical face carries the BRep_PolygonOnClosedTriangulation
//!    double-string structure; BRep_Tool::PolygonOnTriangulation returns the
//!    second string for a REVERSED edge wrapper (BRep_Tool.cxx L700-703), and
//!    BRepTools::Clean removes the mount exactly like OCCT.
//! 3. `poly_connect_triangle_fan_adjacency` — Poly_Connect adjacency on a
//!    known fan triangulation returns the OCCT-documented neighbor/order
//!    semantics (implementation order: the adjacent triangle slots run across
//!    the edges (n1,n2), (n2,n3), (n3,n1), Poly_Connect.cxx L162-190).

use glam::{DVec2, DVec3};
use rcad_brep::poly::{
    PolyConnect, PolyTriangle, PolyTriangulation, PolyTriangulationParameters,
    POLY_MESH_PURPOSE_CALCULATION, POLY_MESH_PURPOSE_NONE, POLY_MESH_PURPOSE_PRESENTATION,
};
use rcad_brep::tools::{
    brep_tools_clean, read_brep_with_triangles_from_string,
    write_brep_with_triangles_to_string, BrepFormatVersion,
};
use rcad_kernel::geom::{CylindricalSurface, Line3, Plane};
use rcad_kernel::topods::{BRep, Orientation, Shape};
use rcad_kernel::{Curve3, Surface3};

/// Build a unit quad face on the z=0 plane: 4 vertices, 4 edges, 1 wire.
fn build_plane_quad_face(brep: &mut BRep) -> Shape {
    let corners = [
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
    ];
    let v: Vec<Shape> = corners.iter().map(|&p| brep.add_tvertex(p)).collect();
    let mut edge_refs = Vec::new();
    for j in 0..4 {
        let a = j;
        let b = (j + 1) % 4;
        let edge_sr = brep.add_tedge(
            Some(Curve3::Line(Line3::new(corners[a], corners[b]))),
            v[a].clone(),
            v[b].clone(),
            [0.0, 1.0],
        );
        let orient = if a < b { Orientation::Forward } else { Orientation::Reversed };
        edge_refs.push(Shape::synthetic(edge_sr.index, orient));
    }
    let wire = brep.add_twire(edge_refs);
    brep.add_tface(
        Some(Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))),
        wire,
        vec![],
        None,
        None,
        vec![],
        true,
    )
}

/// Build the OCCT-style triangulation of the unit quad: 4 nodes, 2 triangles,
/// UV nodes, single-precision normals.
fn quad_triangulation() -> PolyTriangulation {
    let mut tri = PolyTriangulation::with_sizes(4, 2, true, true);
    tri.set_node(1, DVec3::new(0.0, 0.0, 0.0));
    tri.set_node(2, DVec3::new(1.0, 0.0, 0.0));
    tri.set_node(3, DVec3::new(1.0, 1.0, 0.0));
    tri.set_node(4, DVec3::new(0.0, 1.0, 0.0));
    tri.set_triangle(1, PolyTriangle::new(1, 2, 3));
    tri.set_triangle(2, PolyTriangle::new(1, 3, 4));
    tri.set_uv_node(1, DVec2::new(0.0, 0.0));
    tri.set_uv_node(2, DVec2::new(1.0, 0.0));
    tri.set_uv_node(3, DVec2::new(1.0, 1.0));
    tri.set_uv_node(4, DVec2::new(0.0, 1.0));
    tri.set_normal(1, [0.0, 0.0, 1.0]);
    tri.set_normal(2, [0.0, 0.0, 1.0]);
    tri.set_normal(3, [0.0, 0.0, 1.0]);
    tri.set_normal(4, [0.25f32, 0.0, 0.968245836551854326f32]); // non-trivial f32
    tri.set_deflection(0.001);
    tri.set_mesh_purpose(POLY_MESH_PURPOSE_CALCULATION);
    tri.set_parameters(Some(PolyTriangulationParameters::new(0.001, 0.5, 1.0e-3)));
    tri
}

fn find_face(brep: &BRep) -> Shape {
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(&**ts, rcad_kernel::topods::TShape::Face(_)) {
            return Shape::from_parts(ts.clone(), i, 0, Orientation::Forward);
        }
    }
    panic!("no face in brep");
}

fn find_first_edge(brep: &BRep) -> Shape {
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(&**ts, rcad_kernel::topods::TShape::Edge(_)) {
            return Shape::from_parts(ts.clone(), i, 0, Orientation::Forward);
        }
    }
    panic!("no edge in brep");
}

/// Test (a): plane face triangulated programmatically, attached, serialized,
/// deserialized — node/triangle/UV/normal values compare EXACTLY.
#[test]
fn plane_face_triangulation_round_trip_is_lossless() {
    let mut brep = BRep::new();
    let face = build_plane_quad_face(&mut brep);

    // ---- mount the triangulation (BRep_Builder::UpdateFace) ----
    let tri_idx = brep.add_triangulation(quad_triangulation());
    brep.update_face_triangulation(&face, Some(tri_idx), true);

    // The active-triangulation bit (BRep_TFace.cxx L95) is set on mount.
    assert_eq!(
        brep.triangulations[tri_idx].mesh_purpose() & POLY_MESH_PURPOSE_CALCULATION,
        POLY_MESH_PURPOSE_CALCULATION
    );

    // ---- attach one PolygonOnTriangulation on edge 0 ----
    let edge = find_first_edge(&brep);
    let mut poly = PolyTriangulationPoly::from_nodes(&[1, 2]);
    poly.set_parameters(Some(vec![0.0, 1.0]));
    
    poly.set_deflection(0.001);
    brep.update_edge_polygon_on_triangulation(
        &edge,
        Some(poly),
        tri_idx,
        face.location,
    );

    // ---- serialize (v3, with triangles and normals) ----
    let json = write_brep_with_triangles_to_string(
        &brep,
        true,
        true,
        BrepFormatVersion::Version3,
    )
    .expect("write failed");

    // ---- deserialize and compare EXACTLY ----
    let doc = read_brep_with_triangles_from_string(&json).expect("read failed");
    let restored = &doc.shape_set;

    // The pool triangulation survives bit-for-bit (lossless round trip).
    assert_eq!(restored.triangulations.len(), 1);
    let before = &brep.triangulations[tri_idx];
    let after = &restored.triangulations[0];
    assert_eq!(before.nb_nodes(), after.nb_nodes());
    assert_eq!(before.nb_triangles(), after.nb_triangles());
    for j in 1..=before.nb_nodes() {
        assert_eq!(before.node(j), after.node(j), "node {} bit-exact", j);
        assert_eq!(before.uv_node(j), after.uv_node(j), "uv node {} bit-exact", j);
        assert_eq!(before.normal(j), after.normal(j), "normal {} f32-exact", j);
    }
    for j in 1..=before.nb_triangles() {
        assert_eq!(before.triangle(j), after.triangle(j), "triangle {} exact", j);
    }
    assert_eq!(before.deflection(), after.deflection());
    assert_eq!(before.mesh_purpose(), after.mesh_purpose());
    // rcad superset: the OCCT brep text format never persists these two, the
    // rcad document does (inside shape_set).
    assert_eq!(before.parameters(), after.parameters());

    // The face mount survives (TFaceData list + active entry).
    let rface = find_face(restored);
    let (tris, _loc) = restored.face_triangulations(&rface);
    assert_eq!(tris, &[0usize], "face -> section triangulation binding");
    let mounted = restored.triangulation(&rface, POLY_MESH_PURPOSE_NONE);
    assert!(mounted.is_some(), "active triangulation restored");

    // The edge polygon-on-triangulation string survives exactly.
    redge_string_exact(restored);

    // ---- OCCT record-set field checks (v3) ----
    assert_eq!(doc.sections.triangulations.len(), 1);
    let rec = &doc.sections.triangulations[0];
    assert_eq!(rec.nb_nodes, 4);
    assert_eq!(rec.nb_triangles, 2);
    assert_eq!(rec.has_uv_nodes, 1);
    assert_eq!(rec.has_normals, 1);
    assert_eq!(rec.deflection, 0.001);
    assert_eq!(rec.nodes[2], [1.0, 1.0, 0.0]);
    assert_eq!(rec.uv_nodes.as_ref().unwrap()[3], [0.0, 1.0]);
    assert_eq!(rec.triangles, &[[1, 2, 3], [1, 3, 4]]);
    assert_eq!(rec.normals.as_ref().unwrap()[3][0], 0.25f32);

    // The polygon-on-triangulation record: nodes + deflection + parameters.
    assert_eq!(doc.sections.polygon_on_triangulations.len(), 1);
    let prec = &doc.sections.polygon_on_triangulations[0];
    assert_eq!(prec.nb_nodes, 2);
    assert_eq!(prec.nodes, &[1, 2]);
    assert_eq!(prec.deflection, 0.001);
    assert_eq!(prec.has_parameters, 1);
    assert_eq!(prec.parameters.as_ref().unwrap(), &[0.0, 1.0]);

    // The edge record references kind 6 (single string) and triangulation 1.
    assert_eq!(doc.sections.references.edges.len(), 1);
    assert_eq!(doc.sections.references.edges[0].kind, 6);
    assert_eq!(doc.sections.references.edges[0].triangulation, 1);
    assert!(doc.sections.references.edges[0].polygon2.is_none());

    // The face record carries the `2 <tri-index>` tail.
    assert_eq!(doc.sections.references.faces.len(), 1);
    assert_eq!(doc.sections.references.faces[0].triangulation, 1);

    // ---- v2 discriminator: the normals flag does not exist before v3
    // (WriteTriangulation L1605-1608 gates hasNormals on FormatNb >= V3) ----
    let json_v2 = write_brep_with_triangles_to_string(
        &brep,
        true,
        true,
        BrepFormatVersion::Version2,
    )
    .expect("write v2 failed");
    let doc_v2 = read_brep_with_triangles_from_string(&json_v2).expect("read v2 failed");
    assert_eq!(doc_v2.sections.triangulations[0].has_normals, 0);
    assert!(doc_v2.sections.triangulations[0].normals.is_none());
}

/// Re-find the edge in a restored brep and assert the polygon string is
/// bit-exact.
fn redge_string_exact(brep: &BRep) {
    let edge = find_first_edge(brep);
    let p = brep.polygon_on_triangulation(&edge, 0, 0).expect("polygon restored");
    assert_eq!(p.node(1), 1);
    assert_eq!(p.node(2), 2);
    assert_eq!(p.parameter(1), 0.0);
    assert_eq!(p.parameter(2), 1.0);
    assert_eq!(p.deflection(), 0.001);
}

/// Helper alias so the test file reads like the OCCT class name.
use rcad_kernel::poly::PolyPolygonOnTriangulation as PolyTriangulationPoly;

/// Build a cylinder seam face: one lateral face whose wire holds the single
/// seam edge (u = 0 line of the cylinder), like the OCCT cylindrical face
/// after BRepPrimAPI_MakeCylinder (the seam edge appears twice on the face).
fn build_cylinder_seam_face(brep: &mut BRep) -> (Shape, Shape) {
    let radius = 2.0;
    let height = 3.0;
    let cyl = Surface3::Cylinder(CylindricalSurface::new(DVec3::ZERO, DVec3::Z, radius));
    let v_bottom = brep.add_tvertex(DVec3::new(radius, 0.0, 0.0));
    let v_top = brep.add_tvertex(DVec3::new(radius, 0.0, height));
    let seam = brep.add_tedge(
        Some(Curve3::Line(Line3::new(
            DVec3::new(radius, 0.0, 0.0),
            DVec3::new(radius, 0.0, height),
        ))),
        v_bottom.clone(),
        v_top.clone(),
        [0.0, height],
    );
    // A closed seam wire walks the single edge twice (FORWARD up, REVERSED
    // down) — the topological shape of a seam loop.
    let wire = brep.add_twire(vec![
        Shape::synthetic(seam.index, Orientation::Forward),
        Shape::synthetic(seam.index, Orientation::Reversed),
    ]);
    let face = brep.add_tface(Some(cyl), wire, vec![], None, None, vec![], true);
    (face, seam)
}

/// Test (b): the cylinder seam edge carries the DOUBLE
/// PolygonOnTriangulation (BRep_PolygonOnClosedTriangulation); the accessor
/// returns the second string for REVERSED wrappers; persistence keeps both
/// strings; Clean removes them like OCCT.
#[test]
fn cylinder_seam_face_carries_double_polygon_on_triangulation() {
    let mut brep = BRep::new();
    let (face, seam) = build_cylinder_seam_face(&mut brep);

    // Face triangulation (the same table both seam strings index into).
    let mut tri = PolyTriangulation::with_sizes(4, 2, true, false);
    tri.set_node(1, DVec3::new(2.0, 0.0, 0.0));
    tri.set_node(2, DVec3::new(-2.0, 0.0, 0.0));
    tri.set_node(3, DVec3::new(2.0, 0.0, 3.0));
    tri.set_node(4, DVec3::new(-2.0, 0.0, 3.0));
    tri.set_triangle(1, PolyTriangle::new(1, 2, 3));
    tri.set_triangle(2, PolyTriangle::new(2, 4, 3));
    tri.set_uv_node(1, DVec2::new(0.0, 0.0));
    tri.set_uv_node(2, DVec2::new(std::f64::consts::PI, 0.0));
    tri.set_uv_node(3, DVec2::new(0.0, 3.0));
    tri.set_uv_node(4, DVec2::new(std::f64::consts::PI, 3.0));
    tri.set_deflection(0.1);
    let tri_idx = brep.add_triangulation(tri.clone());
    brep.update_face_triangulation(&face, Some(tri_idx), true);

    // The seam double-string: P1 (FORWARD side, u = 0) and P2 (REVERSED
    // side, u = 2pi). Same node indices, different parameter arrays — the
    // two strings are distinct values.
    let mut p1 = PolyTriangulationPoly::from_nodes(&[1, 3]);
    p1.set_parameters(Some(vec![0.0, 3.0]));
    
    p1.set_deflection(0.1);
    let mut p2 = PolyTriangulationPoly::from_nodes(&[1, 3]);
    p2.set_parameters(Some(vec![0.0, 3.0]));
    
    p2.set_deflection(0.2);
    brep.update_edge_polygon_on_closed_triangulation(
        &seam,
        Some(p1),
        Some(p2),
        tri_idx,
        face.location,
    );

    // ---- accessor semantics (BRep_Tool::PolygonOnTriangulation, L683-714):
    // FORWARD wrapper -> P1, REVERSED wrapper -> P2. ----
    let got_p1 = brep
        .polygon_on_triangulation(&seam, tri_idx, face.location)
        .expect("P1 present");
    assert_eq!(got_p1.deflection(), 0.1, "FORWARD orientation reads string 1");

    let seam_rev = Shape::from_parts(seam.data.clone(), seam.index, seam.location, Orientation::Reversed);
    let got_p2 = brep
        .polygon_on_triangulation(&seam_rev, tri_idx, face.location)
        .expect("P2 present");
    assert_eq!(got_p2.deflection(), 0.2, "REVERSED orientation reads string 2");

    // ---- persistence: kind 7 record with BOTH string indices ----
    let json = write_brep_with_triangles_to_string(
        &brep,
        true,
        false,
        BrepFormatVersion::Version3,
    )
    .expect("write failed");
    let doc = read_brep_with_triangles_from_string(&json).expect("read failed");
    let seam_binding = doc
        .sections
        .references
        .edges
        .iter()
        .find(|b| b.kind == 7)
        .expect("kind 7 (closed) record");
    assert_ne!(
        seam_binding.polygon,
        seam_binding.polygon2.unwrap(),
        "the two strings are distinct section records"
    );
    let rec1 = &doc.sections.polygon_on_triangulations[seam_binding.polygon - 1];
    let rec2 = &doc.sections.polygon_on_triangulations[seam_binding.polygon2.unwrap() - 1];
    assert_eq!(rec1.deflection, 0.1);
    assert_eq!(rec2.deflection, 0.2);

    // Both strings restored bit-exact.
    let restored = &doc.shape_set;
    let rseam = find_first_edge(restored);
    let rp1 = restored
        .polygon_on_triangulation(&rseam, 0, 0)
        .expect("P1 restored");
    let rseam_rev =
        Shape::from_parts(rseam.data.clone(), rseam.index, rseam.location, Orientation::Reversed);
    let rp2 = restored
        .polygon_on_triangulation(&rseam_rev, 0, 0)
        .expect("P2 restored");
    assert_eq!(rp1, got_p1);
    assert_eq!(rp2, got_p2);

    // ---- BRepTools::Clean (theForce = false): the face is geometric and
    // carries the triangulation, so the seam polygons AND the face
    // triangulation are removed (BRepTools.cxx L754-767). ----
    let mut brep2 = BRep::new();
    let (face2, seam2) = build_cylinder_seam_face(&mut brep2);
    let tri_idx2 = brep2.add_triangulation(tri.clone());
    brep2.update_face_triangulation(&face2, Some(tri_idx2), true);
    brep2.update_edge_polygon_on_closed_triangulation(
        &seam2,
        Some(PolyTriangulationPoly::from_nodes(&[1, 3])),
        Some(PolyTriangulationPoly::from_nodes(&[1, 3])),
        tri_idx2,
        face2.location,
    );
    brep_tools_clean(&mut brep2, false);
    let seam2_after = find_first_edge(&brep2);
    let reps = match &*brep2.tshapes[seam2_after.index] {
        rcad_kernel::topods::TShape::Edge(ed) => &ed.representations,
        _ => panic!("not an edge"),
    };
    assert!(
        !reps.iter().any(|cr| cr.is_polygon_on_triangulation()),
        "Clean nullifies the edge polygons"
    );
    let rface2 = find_face(&brep2);
    let (tris2, _) = brep2.face_triangulations(&rface2);
    assert!(tris2.is_empty(), "Clean resets the face triangulation list");
    assert!(brep2.triangulation(&rface2, POLY_MESH_PURPOSE_NONE).is_none());

    // ---- BRepTools::Clean force semantics (L793-812): without a mounted
    // face triangulation the polygons survive theForce=false and are removed
    // by theForce=true. ----
    let mut brep3 = BRep::new();
    let (_face3, seam3) = build_cylinder_seam_face(&mut brep3);
    let tri_idx3 = brep3.add_triangulation(tri.clone());
    brep3.update_edge_polygon_on_triangulation(
        &seam3,
        Some(PolyTriangulationPoly::from_nodes(&[1, 2])),
        tri_idx3,
        0,
    );
    brep_tools_clean(&mut brep3, false);
    let seam3_kept = find_first_edge(&brep3);
    let kept = match &*brep3.tshapes[seam3_kept.index] {
        rcad_kernel::topods::TShape::Edge(ed) => &ed.representations,
        _ => panic!("not an edge"),
    };
    assert!(
        kept.iter().any(|cr| cr.is_polygon_on_triangulation()),
        "without force the polygons are kept when no face triangulation exists"
    );
    brep_tools_clean(&mut brep3, true);
    let seam3_forced = find_first_edge(&brep3);
    let forced = match &*brep3.tshapes[seam3_forced.index] {
        rcad_kernel::topods::TShape::Edge(ed) => &ed.representations,
        _ => panic!("not an edge"),
    };
    assert!(
        !forced.iter().any(|cr| cr.is_polygon_on_triangulation()),
        "theForce removes all polygon-on-triangulation representations"
    );
}

/// Test (c): Poly_Connect adjacency on a triangle fan returns the
/// OCCT-documented neighbor/order semantics.
///
/// Fan: hub node 1, rim nodes 2..6,
/// T1(1,2,3) T2(1,3,4) T3(1,4,5) T4(1,5,6).
/// Interior shared edges: (1,3) T1-T2, (1,4) T2-T3, (1,5) T3-T4.
#[test]
fn poly_connect_triangle_fan_adjacency() {
    let mut tri = PolyTriangulation::with_sizes(6, 4, false, false);
    tri.set_node(1, DVec3::ZERO);
    for (i, ang) in (2..=6).zip([0.4f64, 0.9, 1.4, 1.9, 2.4]) {
        tri.set_node(i, DVec3::new(ang.cos(), ang.sin(), 0.0));
    }
    tri.set_triangle(1, PolyTriangle::new(1, 2, 3));
    tri.set_triangle(2, PolyTriangle::new(1, 3, 4));
    tri.set_triangle(3, PolyTriangle::new(1, 4, 5));
    tri.set_triangle(4, PolyTriangle::new(1, 5, 6));

    let conn = PolyConnect::new(&tri);

    // Node -> a containing triangle (the LAST triangle written wins,
    // Poly_Connect.cxx L104-107).
    assert_eq!(conn.triangle(1), 4);
    assert_eq!(conn.triangle(2), 1);
    assert_eq!(conn.triangle(6), 4);

    // T2 = (1,3,4): slot order runs across edges (1,3), (3,4), (4,1)
    // (implementation order, Poly_Connect.cxx L162-190).
    // across (1,3): T1, third node of T1 = 2.
    // across (3,4): boundary -> 0/0.
    // across (4,1): T3, third node of T3 = 5.
    assert_eq!(conn.triangles(2), (1, 0, 3));
    assert_eq!(conn.nodes(2), (2, 0, 5));

    // T1 = (1,2,3): two boundary edges.
    assert_eq!(conn.triangles(1), (0, 0, 2));
    assert_eq!(conn.nodes(1), (0, 0, 4));

    // T3 = (1,4,5): across (1,4): T2 (third node 3); across (4,5): boundary;
    // across (5,1): T4 (third node 6).
    assert_eq!(conn.triangles(3), (2, 0, 4));
    assert_eq!(conn.nodes(3), (3, 0, 6));

    // ---- iterator walk around the hub node 1 (Initialize/More/Next/Value,
    // Poly_Connect.cxx L208-296): the OCCT walk goes right around the fan
    // then left, visiting each containing triangle exactly once. ----
    let mut conn = PolyConnect::new(&tri);
    conn.initialize(1);
    assert!(conn.more());
    assert_eq!(conn.value(), 4, "starts at the last triangle containing node 1");
    conn.next();
    assert!(conn.more());
    assert_eq!(conn.value(), 3);
    conn.next();
    assert!(conn.more());
    assert_eq!(conn.value(), 2);
    conn.next();
    assert!(conn.more());
    assert_eq!(conn.value(), 1);
    conn.next();
    assert!(!conn.more(), "the fan walk terminates after visiting T1..T4");
}
