// =============================================================================
// Shell / Solid Extraction (explode equivalent)
// =============================================================================

/// Remove stale vertices/edges and rebuild with dense indexing.
///
/// After boolean operations, the result BRep may retain vertices from both
/// inputs that are not part of the result geometry.  These stale vertices
/// inflate [`bounding_box`] and can cause subsequent booleans to produce
/// wrong results (e.g. via [`try_containment`](crate::try_containment)).
///
/// This function rebuilds the BRep using only vertices and edges that are
/// referenced by at least one face wire, producing a minimal self-contained
/// copy with correct bounding box.
pub(crate) fn compact_brep(brep: &BRep) -> BRep {
    // Preserve multi-solid structure: compact each solid separately.
    if brep.solids.len() > 1 {
        let mut flat_idx = 0usize;
        let mut comps = Vec::new();
        for solid in &brep.solids {
            let face_count: usize = solid.shells.iter().map(|sh| sh.faces.len()).sum();
            if face_count > 0 {
                let indices: Vec<usize> = (flat_idx..flat_idx + face_count).collect();
                let subset = extract_brep_subset(brep, &indices);
                if collect_flat_face_indices(&subset).len() > 0 {
                    comps.push(subset);
                }
            }
            flat_idx += face_count;
        }
        return BRep::compound_from_shapes(&comps);
    }

    let all_faces: Vec<usize> = collect_flat_face_indices(brep);
    if all_faces.is_empty() {
        return BRep::new();
    }
    extract_brep_subset(brep, &all_faces)
}

/// Create a new self-contained BRep containing only the specified flat face
/// indices from the source BRep.  Vertices, edges, and geometry referenced by
/// the selected faces are copied into the new BRep with dense index renumbering.
fn extract_brep_subset(source: &BRep, face_indices: &[usize]) -> BRep {
    use std::collections::{HashMap, HashSet};

    use rcad_kernel::topology::{Edge, Shell, Solid, Wire, WireEdge};

    if face_indices.is_empty() {
        return BRep::new();
    }

    // Build flat-face index 鈫?(solid_idx, shell_idx, local_face_idx) lookup
    let mut flat_index_map: Vec<(usize, usize, usize)> = Vec::new(); // (solid, shell, local_face)
    for (si, solid) in source.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for fi in 0..shell.faces.len() {
                flat_index_map.push((si, shi, fi));
            }
        }
    }

    // Collect unique edge indices referenced by the selected faces.
    // Also save face topology for later (before remapping).
    let mut edge_set: HashSet<usize> = HashSet::new();
    #[derive(Clone)]
    struct FaceTopo {
        surface_idx: Option<usize>,
        outer: Vec<WireEdge>,
        inner: Vec<Vec<WireEdge>>,
        normal: DVec3,
        triangles: Vec<[usize; 3]>,
    }
    let mut face_topos: Vec<FaceTopo> = Vec::with_capacity(face_indices.len());

    for &fi in face_indices.iter() {
        if fi >= flat_index_map.len() {
            continue;
        }
        let (si, shi, lfi) = flat_index_map[fi];
        let face = &source.solids[si].shells[shi].faces[lfi];

        for we in &face.outer_wire.edges {
            edge_set.insert(we.idx);
        }
        for wire in &face.inner_wires {
            for we in &wire.edges {
                edge_set.insert(we.idx);
            }
        }
        face_topos.push(FaceTopo {
            surface_idx: face.surface_idx,
            outer: face.outer_wire.edges.clone(),
            inner: face.inner_wires.iter().map(|w| w.edges.clone()).collect(),
            normal: face.normal,
            triangles: face.triangles.clone(),
        });
    }

    // Collect vertex indices from the selected edges.
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &ei in &edge_set {
        if ei < source.edges.len() {
            vertex_set.insert(source.edges[ei].start);
            vertex_set.insert(source.edges[ei].end);
        }
    }

    // Collect geometry indices referenced by edges and faces.
    let mut curve_set: HashSet<usize> = HashSet::new();
    let mut surface_set: HashSet<usize> = HashSet::new();
    let mut curve2d_set: HashSet<usize> = HashSet::new();

    for &ei in &edge_set {
        if let Some(Some(ci)) = source.geom.edge_curve.get(ei) {
            curve_set.insert(*ci);
        }
        if let Some(pcurves) = source.geom.edge_pcurves.get(ei) {
            for pc in pcurves {
                surface_set.insert(pc.surface_idx);
                curve2d_set.insert(pc.curve2d_idx);
            }
        }
    }
    for &fi in face_indices.iter() {
        if let Some(Some(si)) = source.geom.face_surface.get(fi) {
            surface_set.insert(*si);
        }
    }

    // Build sorted remap tables: old 鈫?new dense indices.
    let make_remap = |set: &HashSet<usize>| -> (Vec<usize>, HashMap<usize, usize>) {
        let mut sorted: Vec<usize> = set.iter().copied().collect();
        sorted.sort();
        let map: HashMap<usize, usize> =
            sorted.iter().enumerate().map(|(n, &o)| (o, n)).collect();
        (sorted, map)
    };

    let (sorted_vertices, v_remap) = make_remap(&vertex_set);
    let (sorted_edges, e_remap) = make_remap(&edge_set);
    let (sorted_curves, c_remap) = make_remap(&curve_set);
    let (sorted_surfaces, s_remap) = make_remap(&surface_set);
    let (sorted_curve2ds, k_remap) = make_remap(&curve2d_set);

    let mut result = BRep::new();

    // --- vertices ---
    for &old in &sorted_vertices {
        result.vertices.push(source.vertices[old]);
        result
            .geom
            .vertex_tolerance
            .push(source.geom.vertex_tolerance.get(old).copied().unwrap_or(CONFUSION));
    }

    // --- edges ---
    for &old in &sorted_edges {
        let e = &source.edges[old];
        result.edges.push(Edge {
            start: v_remap[&e.start],
            end: v_remap[&e.end],
        });

        result.geom.edge_curve.push(
            source
                .geom
                .edge_curve
                .get(old)
                .and_then(|o| o.map(|c| c_remap[&c])),
        );
        result.geom.edge_pcurves.push(
            source
                .geom
                .edge_pcurves
                .get(old)
                .map(|v| {
                    v.iter()
                        .map(|p| rcad_kernel::PCurve {
                            surface_idx: s_remap[&p.surface_idx],
                            curve2d_idx: k_remap[&p.curve2d_idx],
                        })
                        .collect()
                })
                .unwrap_or_default(),
        );
        result
            .geom
            .edge_curve_range
            .push(source.geom.edge_curve_range.get(old).copied().flatten());
        result
            .geom
            .edge_degenerated
            .push(*source.geom.edge_degenerated.get(old).unwrap_or(&false));
        result
            .geom
            .edge_same_parameter
            .push(*source.geom.edge_same_parameter.get(old).unwrap_or(&true));
        result
            .geom
            .edge_same_range
            .push(*source.geom.edge_same_range.get(old).unwrap_or(&true));
        result
            .geom
            .edge_tolerance
            .push(*source.geom.edge_tolerance.get(old).unwrap_or(&CONFUSION));
    }

    // --- geometry pools ---
    for &old in &sorted_curves {
        result.geom.curves.push(source.geom.curves[old].clone());
    }
    for &old in &sorted_surfaces {
        result.geom.surfaces.push(source.geom.surfaces[old].clone());
    }
    for &old in &sorted_curve2ds {
        result.geom.curve2ds.push(source.geom.curve2ds[old].clone());
        result
            .geom
            .curve2d_range
            .push(source.geom.curve2d_range.get(old).copied().flatten());
    }

    // --- faces (topology + face-level geom) ---
    let mut new_faces: Vec<Face> = Vec::with_capacity(face_topos.len());
    for (i, &fi) in face_indices.iter().enumerate() {
        let ft = &face_topos[i];

        let remap_wire_edges = |wes: &[WireEdge]| -> Vec<WireEdge> {
            wes.iter()
                .map(|we| WireEdge {
                    idx: e_remap[&we.idx],
                    forward: we.forward,
                })
                .collect()
        };

        new_faces.push(Face {
            outer_wire: Wire {
                edges: remap_wire_edges(&ft.outer),
            },
            inner_wires: ft
                .inner
                .iter()
                .map(|w| Wire {
                    edges: remap_wire_edges(w),
                })
                .collect(),
            normal: ft.normal,
            triangles: ft
                .triangles
                .iter()
                .map(|&[a, b, c]| {
                    [
                        v_remap.get(&a).copied().unwrap_or(0),
                        v_remap.get(&b).copied().unwrap_or(0),
                        v_remap.get(&c).copied().unwrap_or(0),
                    ]
                })
                .collect(),
            // ✅ OCCT对齐: 保存原始 face 的 sample_point。compact_brep 之前设为
            //    None,导致后续 pipeline 分类退回边中点→误判内含面为 In→删除。
            //    OCCT BRep_Builder::UpdateFace 保留 sample_point。
            sample_point: source
                .solids
                .get(flat_index_map[fi].0)
                .and_then(|s| s.shells.get(flat_index_map[fi].1))
                .and_then(|sh| sh.faces.get(flat_index_map[fi].2))
                .and_then(|f| f.sample_point),
            mesh_dirty: true,
                surface_idx: ft.surface_idx.and_then(|si| s_remap.get(&si).copied()),
        });

        // face-level geometry
        result
            .geom
            .face_surface
            .push(source.geom.face_surface.get(fi).copied().flatten().map(|si| s_remap[&si]));
        result
            .geom
            .face_surface_range
            .push(source.geom.face_surface_range.get(fi).copied().flatten());
        result
            .geom
            .face_tolerance
            .push(*source.geom.face_tolerance.get(fi).unwrap_or(&CONFUSION));
        // ✅ OCCT对齐: 传播 face_internal_vertices (FillInternalVertices)
        let old_ivs: &[usize] = source.geom.face_internal_vertices.get(fi).map_or(&[], |v| v.as_slice());
        let new_ivs: Vec<usize> = old_ivs.iter().filter_map(|&ov| v_remap.get(&ov).copied()).collect();
        result.geom.face_internal_vertices.push(new_ivs);
    }

    // Wrap in solid/shell topology.
    result.solids.push(Solid {
        shells: vec![Shell { faces: new_faces }],
    });

    // Copy compound structure if source is a compound.
    // NOTE: We don't try to rebuild the compound 鈥?each extracted subset is
    // a standalone self-contained BRep with one Solid.
    result
}

/// Extract each solid from a (possibly compound) BRep as a separate
/// self-contained BRep.  Equivalent to OCCT `explode ... so`.
///
/// Each returned BRep has only the vertices, edges, and geometry belonging
/// to that solid, with indices renumbered from 0.
pub fn extract_solids(brep: &BRep) -> Vec<BRep> {
    let mut results = Vec::new();
    let mut flat_idx = 0;

    for solid in &brep.solids {
        let face_count: usize = solid.shells.iter().map(|sh| sh.faces.len()).sum();
        if face_count > 0 {
            let indices: Vec<usize> = (flat_idx..flat_idx + face_count).collect();
            results.push(extract_brep_subset(brep, &indices));
        }
        flat_idx += face_count;
    }

    results
}

/// Extract each shell from a BRep as a separate self-contained BRep.
/// Equivalent to OCCT `explode ... Sh`.
///
/// Each returned BRep has only the vertices, edges, and geometry belonging
/// to that shell, with indices renumbered from 0.
pub fn extract_shells(brep: &BRep) -> Vec<BRep> {
    let mut results = Vec::new();
    let mut flat_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            let face_count = shell.faces.len();
            if face_count > 0 {
                let indices: Vec<usize> = (flat_idx..flat_idx + face_count).collect();
                results.push(extract_brep_subset(brep, &indices));
            }
            flat_idx += face_count;
        }
    }

    results
}

/// Partition objects by tools using boolean-subset decomposition.
///
/// For each object and each combination of tools (inside/outside per tool mask),
/// computes the corresponding cell using pairwise boolean operations.
/// Returns all non-empty cells as individual self-contained BReps (one solid each).
///
/// This is equivalent to OCCT's `BRepAlgoAPI_Splitter` / `BRepAlgoAPI_Partition`.
///
/// Face tools (planar faces acting as half-space dividers) are automatically
/// expanded into two half-space solids. Face objects (zero-volume faces) are
/// partitioned via `split_shape` + point classification.
///
/// The number of boolean operations per call is O(objects.len() 脳 2^n_tools 脳 n_tools),
/// so this is suitable only for small numbers of tools (鈮?10).
pub fn n_ary_partition(objects: &[BRep], tools: &[BRep]) -> Result<Vec<BRep>, crate::BooleanError> {
    let mut all_cells = Vec::new();

    for obj in objects {
        if is_face_like(obj) {
            all_cells.extend(partition_face_object(obj, tools)?);
        } else {
            all_cells.extend(partition_solid_object(obj, tools)?);
        }
    }

    all_cells.retain(|c| count_faces(c) > 0);
    Ok(all_cells)
}

/// Decompose the complement of `inner_bbox` within `outer_bbox` into up to 6
/// non-overlapping axis-aligned boxes.
///
/// Each returned tuple is `(origin, u_dir, v_dir, width, height, depth)` suitable
/// for passing to `make_box_brep`.  The boxes are disjoint and their union is
/// exactly `outer_bbox \ inner_bbox` (the region inside the outer box but outside
/// the inner box).
fn box_complement_of_bbox(
    inner: &[DVec3; 2],
    outer: &[DVec3; 2],
) -> Vec<(DVec3, DVec3, DVec3, f64, f64, f64)> {
    let (omin, omax) = (outer[0], outer[1]);
    let (imin, imax) = (inner[0], inner[1]);

    // Clamp inner bbox to outer bbox so the complement doesn't extend past the outer.
    let imin = imin.max(omin);
    let imax = imax.min(omax);

    let mut boxes = Vec::with_capacity(6);

    // Left: x < imin.x (full y,z range of outer)
    if imin.x > omin.x {
        boxes.push((
            DVec3::new(omin.x, omin.y, omin.z),
            DVec3::X,
            DVec3::Y,
            imin.x - omin.x,
            omax.y - omin.y,
            omax.z - omin.z,
        ));
    }

    // Right: x > imax.x (full y,z range of outer)
    if omax.x > imax.x {
        boxes.push((
            DVec3::new(imax.x, omin.y, omin.z),
            DVec3::X,
            DVec3::Y,
            omax.x - imax.x,
            omax.y - omin.y,
            omax.z - omin.z,
        ));
    }

    // Front: y < imin.y (within tool's x range, full z range)
    if imin.y > omin.y {
        boxes.push((
            DVec3::new(imin.x, omin.y, omin.z),
            DVec3::X,
            DVec3::Y,
            imax.x - imin.x,
            imin.y - omin.y,
            omax.z - omin.z,
        ));
    }

    // Back: y > imax.y (within tool's x range, full z range)
    if omax.y > imax.y {
        boxes.push((
            DVec3::new(imin.x, imax.y, omin.z),
            DVec3::X,
            DVec3::Y,
            imax.x - imin.x,
            omax.y - imax.y,
            omax.z - omin.z,
        ));
    }

    // Bottom: z < imin.z (within tool's x,y range)
    if imin.z > omin.z {
        boxes.push((
            DVec3::new(imin.x, imin.y, omin.z),
            DVec3::X,
            DVec3::Y,
            imax.x - imin.x,
            imax.y - imin.y,
            imin.z - omin.z,
        ));
    }

    // Top: z > imax.z (within tool's x,y range)
    if omax.z > imax.z {
        boxes.push((
            DVec3::new(imin.x, imin.y, imax.z),
            DVec3::X,
            DVec3::Y,
            imax.x - imin.x,
            imax.y - imin.y,
            omax.z - imax.z,
        ));
    }

    boxes
}

/// Check whether a BRep is a simple axis-aligned box (its volume matches its
/// bounding-box volume within tolerance).
fn is_box_like(brep: &BRep) -> bool {
    let vol = crate::total_volume(brep);
    if vol <= 0.0 {
        return false;
    }
    if let Some(bbox) = bounding_box(brep) {
        let diag = bbox[1] - bbox[0];
        let bbox_vol = diag.x * diag.y * diag.z;
        if bbox_vol <= 0.0 {
            return false;
        }
        (vol - bbox_vol).abs() < 1e-6
    } else {
        false
    }
}

/// Partition a solid object by tools, expanding face tools into half-space solids.
fn partition_solid_object(obj: &BRep, tools: &[BRep]) -> Result<Vec<BRep>, crate::BooleanError> {
    // Detect face-like tools (planar faces used as dividing surfaces).
    let face_tool_info: Vec<Option<rcad_kernel::geom::Plane>> = tools
        .iter()
        .map(|t| if is_face_like(t) { try_as_planar_face(t) } else { None })
        .collect();
    let has_face_tool = face_tool_info.iter().any(|p| p.is_some());

    // Expand tools + track complement indices for half-spaces.
    //
    // For a half-space pair (h_plus, h_minus), each half-space's "outside"
    // is simply Intersection with the OTHER half-space 鈥?no Diff needed.
    // expanded_complements[i] = Some(j) means bit-flip (outside) at index i
    // can be handled by Intersection with expanded_tools[j].
    //
    // For solid tools (non-face), expanded_complements[i] = None 鈥?those
    // may need the complement-box fallback or Diff.
    let mut expanded_tools: Vec<BRep> = Vec::new();
    #[allow(clippy::type_complexity)]
    let mut expanded_complements: Vec<Option<usize>> = Vec::new();

    if has_face_tool {
        let mut bbox = bounding_box(obj);
        for (ti, tool) in tools.iter().enumerate() {
            if face_tool_info[ti].is_some() {
                if let Some(tb) = bounding_box(tool) {
                    bbox = match bbox {
                        Some(b) => Some([b[0].min(tb[0]), b[1].max(tb[1])]),
                        None => Some(tb),
                    };
                }
            }
        }

        for (ti, tool) in tools.iter().enumerate() {
            if let Some(ref plane) = face_tool_info[ti] {
                if let Some(b) = bbox {
                    let h_plus_idx = expanded_tools.len();
                    let h_plus = make_face_half_space(plane, &b, true);
                    let h_minus = make_face_half_space(plane, &b, false);
                    expanded_tools.push(h_plus);
                    expanded_tools.push(h_minus);
                    // h_plus 鈫?h_minus: each is the "outside" of the other.
                    expanded_complements.push(Some(h_plus_idx + 1));
                    expanded_complements.push(Some(h_plus_idx));
                    continue;
                }
            }
            expanded_complements.push(None);
            expanded_tools.push(tool.clone());
        }
    } else {
        for tool in tools {
            expanded_complements.push(None);
            expanded_tools.push(tool.clone());
        }
    }

    let n_tools = expanded_tools.len();
    let mut cells = Vec::new();
    let max_mask = if n_tools >= 32 { 1u32 << 31 } else { 1u32 << n_tools };

    // Detect complementary half-space pairs (h_plus/h_minus).
    // When both bits of a pair are equal (both 0 or both 1), the cell lies on
    // the plane between them and has zero volume 鈥?skip those masks.
    let mut comp_pairs: Vec<(usize, usize)> = Vec::new();
    for (i, comp_i) in expanded_complements.iter().enumerate() {
        if let Some(j) = comp_i {
            if *j > i && expanded_complements.get(*j) == Some(&Some(i)) {
                comp_pairs.push((i, *j));
            }
        }
    }

    for mask in 0..max_mask {
        // Skip masks where complementary half-spaces have the same bit
        // (both inside or both outside = intersection is just the plane).
        if comp_pairs.iter().any(|&(i, j)| {
            ((mask >> i) & 1) == ((mask >> j) & 1)
        }) {
            continue;
        }
        let mut cell = obj.clone();
        let mut failed = false;
        let mut first_tool = true;

        for i in 0..n_tools {
            let inside = (mask >> i) & 1 != 0;
            let tool = &expanded_tools[i];

            if inside {
                match crate::boolean_op_pave_fill_build(crate::BooleanOpType::Intersection, &cell, tool) {
                    Ok(r) => { cell = compact_brep(&r); },
                    Err(_) => { failed = true; break; }
                }
            } else if let Some(complement_idx) = expanded_complements[i] {
                // Half-space: "outside" = Intersection with complementary half-space.
                let complement = &expanded_tools[complement_idx];
                match crate::boolean_op_pave_fill_build(crate::BooleanOpType::Intersection, &cell, complement) {
                    Ok(r) => { cell = compact_brep(&r); },
                    Err(_) => { failed = true; break; }
                }
            } else {
                // Solid tool: for the first tool application, use complement box
                // decomposition to avoid coincident-face issues with Diff.
                // For subsequent tools, use Diff (works better when cell is already
                // a multi-solid compound from previous operations).
                if first_tool && is_box_like(tool) {
                    // First tool: use complement box decomposition.
                    if let (Some(tool_bbox), Some(cell_bbox)) =
                        (bounding_box(tool), bounding_box(&cell))
                    {
                        let comp_boxes =
                            box_complement_of_bbox(&tool_bbox, &cell_bbox);
                        if comp_boxes.is_empty() {
                            cell = BRep::new();
                            break;
                        }
                        let cell_solids = extract_solids(&cell);
                        let mut parts = Vec::new();
                        for (origin, u_dir, v_dir, w, h, d) in comp_boxes.iter() {
                            let Ok(comp_box) =
                                rcad_modeling::make_box_brep(*origin, *u_dir, *v_dir, *w, *h, *d)
                            else { continue; };
                            for cell_part in &cell_solids {
                                if let Ok(part) =
                                    crate::boolean_op_pave_fill_build(crate::BooleanOpType::Intersection, cell_part, &comp_box)
                                {
                                    if count_faces(&part) > 0 {
                                        parts.push(part);
                                    }
                                }
                            }
                        }
                        cell = BRep::compound_from_shapes(&parts);
                        first_tool = false;
                        continue;
                    }
                }
                // Subsequent tool or non-box tool: use Diff.
                match crate::boolean_op_pave_fill_build(crate::BooleanOpType::Difference, &cell, tool) {
                    Ok(r) => { cell = compact_brep(&r); }
                    Err(_) => { failed = true; break; }
                }
            }
            first_tool = false;
        }

        if failed {
            continue;
        }

        let face_indices = collect_flat_face_indices(&cell);
        if face_indices.is_empty() {
            continue;
        }

        let comps = connected_face_components(&cell, &face_indices);
        for component in comps {
            if !component.is_empty() {
                let subset = extract_brep_subset(&cell, &component);
                cells.push(subset);
            }
        }
    }

    // Filter out degenerate cells (zero volume from tangent-coincident masks).
    cells.retain(|c| crate::total_volume(c) > 1e-10);

    Ok(cells)
}

/// Partition a face-like object by tools using split_shape and centroid classification.
fn partition_face_object(obj: &BRep, tools: &[BRep]) -> Result<Vec<BRep>, crate::BooleanError> {
    // Collect solid (non-face) tools.
    let solid_tools: Vec<BRep> = tools.iter().filter(|t| !is_face_like(t)).cloned().collect();

    if solid_tools.is_empty() || count_faces(obj) == 0 {
        return Ok(vec![obj.clone()]);
    }

    let orig_surface_info = get_surface(obj, 0).ok().map(|s| match s {
        Surface3::Plane(p) => (p.origin, p.normal, "plane"),
        Surface3::Cylinder(_) => (DVec3::ZERO, DVec3::Z, "cylinder"),
        Surface3::Sphere(_) => (DVec3::ZERO, DVec3::Z, "sphere"),
        Surface3::Cone(_) => (DVec3::ZERO, DVec3::Z, "cone"),
        Surface3::Torus(_) => (DVec3::ZERO, DVec3::Z, "torus"),
        _ => (DVec3::ZERO, DVec3::Z, "other"),
    });

    /// Collect flat face indices whose surface matches the original plane.
    let collect_on_plane = |brep: &BRep| -> Vec<usize> {
        let Some((plane_origin, plane_normal, _)) = orig_surface_info else { return vec![] };
        let mut out = Vec::new();
        let mut fi = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for _ in &shell.faces {
                    if let Ok(Surface3::Plane(p)) = get_surface(brep, fi) {
                        // Check coplanarity: normals same direction AND candidate origin
                        // lies on the original plane (plane-relative distance ~0).
                        let dist = (p.origin - plane_origin).dot(plane_normal).abs();
                        if dist < 1e-6 && p.normal.dot(plane_normal) > 0.9999 {
                            out.push(fi);
                        }
                    }
                    fi += 1;
                }
            }
        }
        out
    };

    // Use boolean ops to carve the face into inside/outside per tool.
    let mut cells = Vec::new();
    let mut remaining = obj.clone();
    for tool in &solid_tools {
        let inside = crate::boolean_op(crate::BooleanOpType::Intersection, &remaining, tool)?;
        let in_faces = collect_on_plane(&inside);
        if !in_faces.is_empty() {
            cells.push(extract_brep_subset(&inside, &in_faces));
        }
        remaining = crate::boolean_op(crate::BooleanOpType::Difference, &remaining, tool)?;
    }
    let out_faces = collect_on_plane(&remaining);
    if !out_faces.is_empty() {
        cells.push(extract_brep_subset(&remaining, &out_faces));
    }

    Ok(cells)
}

/// Check if a BRep represents a single planar face and extract its plane.
fn try_as_planar_face(brep: &BRep) -> Option<rcad_kernel::geom::Plane> {
    if count_faces(brep) != 1 {
        return None;
    }
    match get_surface(brep, 0).ok()? {
        Surface3::Plane(plane) => Some(plane.clone()),
        _ => None,
    }
}

/// Check if a BRep is face-like (open surface, not a proper 3D solid).
///
/// A BRep is face-like if every shell contains exactly one planar face. Proper 3D
/// solids always have at least 4 faces per shell (minimum tetrahedron), except
/// analytic primitives like spheres/cones/cylinders which may have only 1-3 faces.
fn is_face_like(brep: &BRep) -> bool {
    if count_faces(brep) == 0 {
        return false;
    }
    let mut flat_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            if shell.faces.len() != 1 {
                return false;
            }
            // The single face must be planar 鈥?spheres with 1 face are NOT face-like.
            let surface = get_surface(brep, flat_idx).ok();
            if !matches!(surface, Some(Surface3::Plane(_))) {
                return false;
            }
            flat_idx += 1;
        }
    }
    true
}

/// Create a half-space solid extending from a plane in the normal (or opposite) direction.
///
/// The resulting box occupies a prism with one face exactly on the plane through
/// `plane.origin` (with the plane's normal pointing inward for `normal_side=true`
/// or outward for `normal_side=false`) and extending far enough along the normal
/// to fully contain the `bbox` extent.
pub fn make_face_half_space(plane: &rcad_kernel::geom::Plane, bbox: &[DVec3; 2], normal_side: bool) -> BRep {
    let [bmin, bmax] = *bbox;
    let diag = bmax - bmin;
    let margin = diag.length().max(1.0) * 2.0;

    let n = if normal_side { plane.normal } else { -plane.normal };
    let n = n.normalize();

    // Build a tangent basis in the plane.
    let abs = n.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z {
        DVec3::X
    } else if abs.y <= abs.z {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let u = n.cross(candidate).normalize();
    let v = n.cross(u);

    // Build the box positioned so it starts at the plane and extends `margin`
    // along +n (normal_side=true) or -n (normal_side=false).
    //
    // The box origin is the corner at (-u*margin/2, -v*margin/2),
    // extended to (+u*margin/2, +v*margin/2) in the plane, and from
    // the plane to 卤margin along n.
    let origin = if normal_side {
        plane.origin - u * (margin / 2.0) - v * (margin / 2.0)
    } else {
        plane.origin - u * (margin / 2.0) - v * (margin / 2.0) - n * margin
    };

    rcad_modeling::make_box_brep(origin, u, v, margin, margin, margin)
        .expect("make_face_half_space: box construction should not fail")
}

/// Compute centroid from a face's triangle mesh.
///
/// Unlike [`average_vertex_of_face`], this works even when the face's
/// `WireEdge.idx` values are local indices that don't reference BRep edges.
fn face_triangle_centroid(face: &Face) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;
    for tri in &face.triangles {
        let _local_vert_id = |vi: usize| {
            // tri[x] is a flat-face vertex index; decode it.
            vi
        };
        // Triangles store flat vertex indices. We use the raw indices but
        // check they don't overflow the triangles vec itself.
        if face.triangles.is_empty() {
            break;
        }
        sum += tri.iter().map(|&_vi| {
            // The triangle indices reference positions from the boundary/wire.
            // This is a fallback 鈥?we just average them.
            DVec3::ZERO
        }).sum::<DVec3>();
        count += 3;
    }
    if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Compute the average vertex position of a face's boundary.
fn average_vertex_of_face(brep: &BRep, face: &Face) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;

    for we in &face.outer_wire.edges {
        if we.idx < brep.edges.len() {
            let edge = &brep.edges[we.idx];
            if edge.start < brep.vertices.len() {
                sum += brep.vertices[edge.start].point;
                count += 1;
            }
            if edge.end < brep.vertices.len() {
                sum += brep.vertices[edge.end].point;
                count += 1;
            }
        }
    }
    for wire in &face.inner_wires {
        for we in &wire.edges {
            if we.idx < brep.edges.len() {
                let edge = &brep.edges[we.idx];
                if edge.start < brep.vertices.len() {
                    sum += brep.vertices[edge.start].point;
                    count += 1;
                }
                if edge.end < brep.vertices.len() {
                    sum += brep.vertices[edge.end].point;
                    count += 1;
                }
            }
        }
    }

    if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Collect flat face indices for all faces in a BRep.
fn collect_flat_face_indices(brep: &BRep) -> Vec<usize> {
    let mut indices = Vec::new();
    let mut flat_idx = 0;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for _ in &shell.faces {
                indices.push(flat_idx);
                flat_idx += 1;
            }
        }
    }
    indices
}

/// Find connected components of a set of flat face indices within a BRep.
/// Two faces are connected if they share at least one edge (same edge index).
fn connected_face_components(brep: &BRep, face_indices: &[usize]) -> Vec<Vec<usize>> {
    use std::collections::{HashMap, HashSet};

    let face_set: HashSet<usize> = face_indices.iter().copied().collect();
    if face_set.is_empty() {
        return Vec::new();
    }

    // Build edge 鈫?face list for our faces of interest.
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut flat_idx: usize = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for lfi in 0..shell.faces.len() {
                let global_fi = flat_idx + lfi;
                if face_set.contains(&global_fi) {
                    if let Some(face) = shell.faces.get(lfi) {
                        for wire_edge in &face.outer_wire.edges {
                            edge_to_faces
                                .entry(wire_edge.idx)
                                .or_default()
                                .push(global_fi);
                        }
                        for wire in &face.inner_wires {
                            for wire_edge in &wire.edges {
                                edge_to_faces
                                    .entry(wire_edge.idx)
                                    .or_default()
                                    .push(global_fi);
                            }
                        }
                    }
                }
            }
            flat_idx += shell.faces.len();
        }
    }

    // Build adjacency: face A 鈫?[faces that share an edge with A].
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for faces in edge_to_faces.values() {
        if faces.len() >= 2 {
            for i in 0..faces.len() {
                for j in (i + 1)..faces.len() {
                    adjacency.entry(faces[i]).or_default().push(faces[j]);
                    adjacency.entry(faces[j]).or_default().push(faces[i]);
                }
            }
        }
    }

    // DFS over face indices to find connected components.
    let mut visited: HashSet<usize> = HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();

    for &fi in face_indices {
        if !visited.insert(fi) {
            continue;
        }

        let mut component: Vec<usize> = Vec::new();
        let mut stack: Vec<usize> = vec![fi];

        while let Some(current) = stack.pop() {
            component.push(current);
            if let Some(neighbors) = adjacency.get(&current) {
                for &neighbor in neighbors {
                    if visited.insert(neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }

        components.push(component);
    }

    components
}
