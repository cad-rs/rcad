/// Convert a DS (DataSource) into a topods BRep.
///
/// This is the transitional bridge between the index-based DS data model
/// and the TopoDS-style BRep shape hierarchy.  PaveFiller writes to DS;
/// `ds_to_brep` converts the result into a BRep so downstream pipeline
/// stages can use BRepTool queries instead of direct DS access.
use rcad_kernel::geom::*;
use rcad_kernel::topods::*;

/// Convert the entire DS into a BRep shape hierarchy.
pub fn ds_to_brep(ds: &crate::bopds::ds::DS) -> BRep {
    let mut br = BRep::new();

    // Step 1: collect surfaces from faces (dedup by identity)
    let mut surf_to_idx: Vec<Option<usize>> = vec![None; ds.faces.len()];
    for (fi, face) in ds.faces.iter().enumerate() {
        let si = find_or_add_surface(&mut br, &face.surface);
        surf_to_idx[fi] = Some(si);
    }

    // Step 2: collect 3D curves from edges (dedup by identity)
    let mut edge_curve_idx: Vec<Option<usize>> = vec![None; ds.edges.len()];
    for (ei, edge) in ds.edges.iter().enumerate() {
        let ci = find_or_add_curve3(&mut br, &edge.curve);
        edge_curve_idx[ei] = Some(ci);
    }

    // Step 3: vertices
    for v in &ds.vertices {
        let r = br.add_tvertex(v.point);
        br.vertex_mut(r).tolerance = v.geom_tol.max(crate::tolerance::TOLERANCE_ABS);
    }

    // Step 4: edges with pcurves, vertex_params, degenerated
    for (ei, edge) in ds.edges.iter().enumerate() {
        let curve_idx = edge_curve_idx[ei];
        let sv = ShapeRef::new(edge.start_vertex);
        let ev = ShapeRef::with_orientation(edge.end_vertex, Orientation::Reversed);
        let e = br.add_tedge(curve_idx, sv, ev, edge.t_range);
        let ed = br.edge_mut(e);
        ed.degenerated = ds.is_edge_degenerated(ei);
        // Copy pcurves from face_reps
        for rep in &edge.face_reps {
            ed.pcurves.insert(rep.face_idx,
                (rep.pcurve.clone(), rep.start_param, rep.end_param));
            // Also add second pcurve if present
            if let Some(ref pc2) = rep.pcurve2 {
                // For periodic surfaces, store second pcurve with shifted face index
                let shifted_fi = rep.face_idx + ds.faces.len();
                ed.pcurves.insert(shifted_fi, (pc2.clone(), rep.start_param, rep.end_param));
            }
        }
        // Copy vertex_params
        ed.vertex_params = edge.vertex_params.clone();
    }

    // Step 5: intersection curves → section edges (stored as standalone edges, not in wires)
    for (ci, ic) in ds.intersection_curves.iter().enumerate() {
        let curve_idx = find_or_add_curve3(&mut br, &ic.curve);
        let sv = ShapeRef::new(ic.start_vertex);
        let ev = ShapeRef::with_orientation(ic.end_vertex, Orientation::Reversed);
        let e = br.add_tedge(Some(curve_idx), sv, ev, ic.t_range);
        let ed = br.edge_mut(e);
        if let Some(ref pc) = ic.pcurve_on_a {
            ed.pcurves.insert(0, (pc.clone(), ic.t_range[0], ic.t_range[1]));
        }
        if let Some(ref pc) = ic.pcurve_on_b {
            ed.pcurves.insert(1, (pc.clone(), ic.t_range[0], ic.t_range[1]));
        }
        ed.degenerated = false;
    }

    // Step 6: build wires and faces from DS topology
    for (fi, face) in ds.faces.iter().enumerate() {
        // Build outer wire from boundary_edges
        let mut wire_edges: Vec<ShapeRef> = Vec::new();
        for &ei in &face.boundary_edges {
            // Determine orientation: forward if edge's start_vertex matches the wire chain
            wire_edges.push(ShapeRef::new(ei));
        }
        let outer_wire = br.add_twire(wire_edges);

        // Build inner wires
        let mut inner_wires: Vec<ShapeRef> = Vec::new();
        for iw in &face.inner_boundary_edges {
            let iw_edges: Vec<ShapeRef> = iw.iter().map(|(ei, _)| ShapeRef::new(*ei)).collect();
            inner_wires.push(br.add_twire(iw_edges));
        }

        let surface = surf_to_idx[fi];
        let internal_vertices: Vec<ShapeRef> = face.face_info.vertices_in.iter()
            .map(|&vi| ShapeRef::new(vi)).collect();

        br.add_tface(surface, outer_wire, inner_wires, None, None, internal_vertices);
    }

    // Step 7: build shells from DS.shells
    for (_shi, shell) in ds.shells.iter().enumerate() {
        let mut faces: Vec<ShapeRef> = Vec::new();
        for &fi in &shell.faces {
            faces.push(ShapeRef::new(fi));
        }
        br.add_tshell(faces);
    }

    br
}

/// Find or add a Surface3 to the BRep's surface pool.
fn find_or_add_surface(br: &mut BRep, surf: &Surface3) -> usize {
    // OCCT handles surface identity via TShapes; rcad uses pointer identity.
    // Use surface type + geometric comparison for dedup.
    for (i, s) in br.surfaces.iter().enumerate() {
        if surfaces_equal(s, surf) {
            return i;
        }
    }
    let idx = br.surfaces.len();
    br.surfaces.push(surf.clone());
    idx
}

/// Find or add a Curve3 to the BRep's curve pool.
fn find_or_add_curve3(br: &mut BRep, curve: &Curve3) -> usize {
    for (i, c) in br.curves.iter().enumerate() {
        if curves_equal(c, curve) {
            return i;
        }
    }
    let idx = br.curves.len();
    br.curves.push(curve.clone());
    idx
}

fn surfaces_equal(a: &Surface3, b: &Surface3) -> bool {
    match (a, b) {
        (Surface3::Plane(pa), Surface3::Plane(pb)) =>
            (pa.origin - pb.origin).length_squared() < 1e-20
            && pa.normal.dot(pb.normal) > 0.999999,
        (Surface3::Sphere(sa), Surface3::Sphere(sb)) =>
            (sa.center - sb.center).length_squared() < 1e-20
            && (sa.radius - sb.radius).abs() < 1e-20,
        (Surface3::Cylinder(ca), Surface3::Cylinder(cb)) =>
            (ca.origin - cb.origin).length_squared() < 1e-20
            && ca.axis.dot(cb.axis) > 0.999999
            && (ca.radius - cb.radius).abs() < 1e-20,
        _ => std::ptr::eq(a as *const _, b as *const _),
    }
}

fn curves_equal(a: &Curve3, b: &Curve3) -> bool {
    match (a, b) {
        (Curve3::Line(la), Curve3::Line(lb)) =>
            (la.origin - lb.origin).length_squared() < 1e-20
            && la.direction.dot(lb.direction) > 0.999999,
        (Curve3::Circle(ca), Curve3::Circle(cb)) =>
            (ca.center - cb.center).length_squared() < 1e-20
            && (ca.radius - cb.radius).abs() < 1e-20
            && ca.normal.dot(cb.normal) > 0.999999,
        (Curve3::Ellipse(ea), Curve3::Ellipse(eb)) =>
            ea.major_radius == eb.major_radius && ea.minor_radius == eb.minor_radius,
        _ => std::ptr::eq(a as *const _, b as *const _),
    }
}
