/// Convert a DS (DataSource) into a topods BRep.
///
/// This is the transitional bridge between the index-based DS data model
/// and the TopoDS-style BRep shape hierarchy.  PaveFiller writes to DS;
/// `ds_to_brep` converts the result into a BRep so downstream pipeline
/// stages can use BRepTool queries instead of direct DS access.
use rcad_kernel::geom::*;
use rcad_kernel::topods::*;

/// Convert the entire DS into a BRep shape hierarchy.
/// Returns (BRep, Vec<ShapeRef>) where the Vec maps DS face_idx -> BRep face ShapeRef.
pub fn ds_to_brep(ds: &crate::bopds::ds::DS) -> (BRep, Vec<ShapeRef>) {
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

    let v_count = ds.vertices.len();
    let e_base = v_count;
    let ic_base = e_base + ds.edges.len();

    // Step 4: edges with pcurves (keyed by DS face_idx — will be rekeyed in Step 7)
    for (ei, edge) in ds.edges.iter().enumerate() {
        let curve_idx = edge_curve_idx[ei];
        let sv = ShapeRef::new(edge.start_vertex);
        let ev = ShapeRef::with_orientation(edge.end_vertex, Orientation::Reversed);
        let e = br.add_tedge(curve_idx, sv, ev, edge.t_range);
        let ed = br.edge_mut(e);
        ed.degenerated = ds.is_edge_degenerated(ei);
        for rep in &edge.face_reps {
            ed.pcurves.insert(rep.face_idx,
                (rep.pcurve.clone(), rep.start_param, rep.end_param));
            if let Some(ref pc2) = rep.pcurve2 {
                let shifted_fi = rep.face_idx + ds.faces.len();
                ed.pcurves.insert(shifted_fi, (pc2.clone(), rep.start_param, rep.end_param));
            }
        }
        ed.vertex_params = edge.vertex_params.clone();
    }

    // Step 5: intersection curves -> section edges (keyed by DS face_idx, rekeyed in Step 7)
    for (ci, ic) in ds.intersection_curves.iter().enumerate() {
        let curve_idx = find_or_add_curve3(&mut br, &ic.curve);
        let sv = ShapeRef::new(ic.start_vertex);
        let ev = ShapeRef::with_orientation(ic.end_vertex, Orientation::Reversed);
        let e = br.add_tedge(Some(curve_idx), sv, ev, ic.t_range);
        let ed = br.edge_mut(e);
        // Find the two faces for this IC to use correct DS face_idx keys
        let face_idxs = find_face_idxs_for_curve(ds, ci);
        if let Some(ref pc) = ic.pcurve_on_a {
            if face_idxs[0] != usize::MAX {
                ed.pcurves.insert(face_idxs[0], (pc.clone(), ic.t_range[0], ic.t_range[1]));
            }
        }
        if let Some(ref pc) = ic.pcurve_on_b {
            if face_idxs[1] != usize::MAX {
                ed.pcurves.insert(face_idxs[1], (pc.clone(), ic.t_range[0], ic.t_range[1]));
            }
        }
        ed.degenerated = false;
    }

    // Step 6: build wires and faces from DS topology, collect face_refs mapping
    let mut face_refs: Vec<ShapeRef> = Vec::with_capacity(ds.faces.len());
    for (fi, face) in ds.faces.iter().enumerate() {
        // Build outer wire from boundary_edges
        let outer_edges: Vec<ShapeRef> = face.boundary_edges.iter()
            .map(|ei| ShapeRef::new(e_base + *ei)).collect();
        let outer_wire = br.add_twire(outer_edges);

        // Build inner wires
        let inner_wires: Vec<ShapeRef> = face.inner_boundary_edges.iter().map(|iw| {
            let iw_edges: Vec<ShapeRef> = iw.iter()
                .map(|(ei, _)| ShapeRef::new(e_base + *ei)).collect();
            br.add_twire(iw_edges)
        }).collect();

        let surface = surf_to_idx[fi].unwrap();
        let internal_vertices: Vec<ShapeRef> = face.face_info.vertices_in.iter()
            .map(|&vi| ShapeRef::new(vi)).collect();

        let face_ref = br.add_tface(Some(surface), outer_wire, inner_wires,
            None, None, internal_vertices);
        face_refs.push(face_ref);
    }

    // Step 7: rekey pcurves — replace DS face_idx keys with BRep face TShape indices
    // DS edges
    for (ei, _edge) in ds.edges.iter().enumerate() {
        let e_ref = ShapeRef::new(e_base + ei);
        let ed = br.edge_mut(e_ref);
        let old_pcurves = std::mem::take(&mut ed.pcurves);
        for (ds_fi, pc) in old_pcurves {
            if ds_fi < face_refs.len() {
                ed.pcurves.insert(face_refs[ds_fi].index, pc);
            } else {
                // Shifted key (pcurve2 for periodic surfaces): subtract ds.faces.len()
                let original_fi = ds_fi - ds.faces.len();
                if original_fi < face_refs.len() {
                    ed.pcurves.insert(face_refs[original_fi].index, pc);
                }
            }
        }
    }
    // IC edges
    for (ci, _ic) in ds.intersection_curves.iter().enumerate() {
        let e_ref = ShapeRef::new(ic_base + ci);
        let ed = br.edge_mut(e_ref);
        let old_pcurves = std::mem::take(&mut ed.pcurves);
        for (ds_fi, pc) in old_pcurves {
            if ds_fi < face_refs.len() {
                ed.pcurves.insert(face_refs[ds_fi].index, pc);
            }
        }
    }

    // Step 8: build shells from DS.shells
    for shell in &ds.shells {
        let faces: Vec<ShapeRef> = shell.faces.iter()
            .map(|fi| face_refs[*fi]).collect();
        br.add_tshell(faces);
    }

    (br, face_refs)
}

/// Find which two DS faces are connected by an intersection curve.
/// Returns [faceA, faceB], where the indices correspond to pcurve_on_a / pcurve_on_b.
fn find_face_idxs_for_curve(ds: &crate::bopds::ds::DS, ci: usize) -> [usize; 2] {
    let mut result = [usize::MAX; 2];
    let mut idx = 0;
    for (fi, face) in ds.faces.iter().enumerate() {
        if face.face_info.curves_sc.contains(&ci) {
            if idx < 2 {
                result[idx] = fi;
                idx += 1;
            }
        }
    }
    result
}

/// Find or add a Surface3 to the BRep's surface pool.
fn find_or_add_surface(br: &mut BRep, surf: &Surface3) -> usize {
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
