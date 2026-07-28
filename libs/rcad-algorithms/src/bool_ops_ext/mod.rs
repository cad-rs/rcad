// Boolean-ops extension: solid/shell extraction, n-ary partition, half-space.
// These were formerly include!ed from brep_tools; kept here as rcad-algorithms
// internal utilities (not part of TKBRep BRepTools).

use rcad_brep::tools::{bounding_box, count_faces, get_surface};
use rcad_kernel::APPROXIMATION;
use rcad_kernel::topods;
use glam::DVec3;

// Topods-native extraction helpers (same module scope via include):
// compact_brep_topods, extract_solids_topods, extract_shells_topods
include!("topods_ext.rs");

/// Extract each solid from a BRep as a separate self-contained BRep.
pub fn extract_solids(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
    extract_solids_topods(brep)
}

/// Extract each shell from a BRep as a separate self-contained BRep.
pub fn extract_shells(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
    extract_shells_topods(brep)
}

/// Partition objects by tools using boolean-subset decomposition.
pub fn n_ary_partition(
    objects: &[rcad_kernel::BRep],
    tools: &[rcad_kernel::BRep],
) -> Result<Vec<rcad_kernel::BRep>, crate::BooleanError> {
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

pub fn make_face_half_space(
    plane: &rcad_kernel::geom::Plane,
    bbox: &[DVec3; 2],
    normal_side: bool,
) -> topods::BRep {
    let [bmin, bmax] = *bbox;
    let diag = bmax - bmin;
    let margin = diag.length().max(1.0) * 2.0;
    let n = if normal_side { plane.normal } else { -plane.normal }.normalize();
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
    let origin = if normal_side {
        plane.origin - u * (margin / 2.0) - v * (margin / 2.0)
    } else {
        plane.origin - u * (margin / 2.0) - v * (margin / 2.0) - n * margin
    };
    rcad_modeling::make_box_brep(origin, u, v, margin, margin, margin)
        .expect("make_face_half_space: box construction failed")
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn compact_brep(brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
    compact_brep_topods(brep)
}

fn extract_brep_subset(source: &rcad_kernel::BRep, _face_indices: &[usize]) -> rcad_kernel::BRep {
    source.clone()
}

fn is_face_like(brep: &rcad_kernel::BRep) -> bool {
    if count_faces(brep) == 0 {
        return false;
    }
    let mut flat_idx = 0usize;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    if shd.faces.len() != 1 {
                        return false;
                    }
                    let surface = get_surface(brep, flat_idx).ok();
                    if !matches!(surface, Some(rcad_kernel::geom::Surface3::Plane(_))) {
                        return false;
                    }
                    flat_idx += 1;
                }
            }
        }
    }
    true
}

fn try_as_planar_face(brep: &rcad_kernel::BRep) -> Option<rcad_kernel::geom::Plane> {
    if count_faces(brep) != 1 {
        return None;
    }
    match get_surface(brep, 0).ok()? {
        rcad_kernel::geom::Surface3::Plane(plane) => Some(plane.clone()),
        _ => None,
    }
}

fn is_box_like(brep: &rcad_kernel::BRep) -> bool {
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

fn box_complement_of_bbox(
    inner: &[DVec3; 2],
    outer: &[DVec3; 2],
) -> Vec<(DVec3, DVec3, DVec3, f64, f64, f64)> {
    let (omin, omax) = (outer[0], outer[1]);
    let (imin, imax) = (inner[0], inner[1]);
    let imin = imin.max(omin);
    let imax = imax.min(omax);
    let mut boxes = Vec::with_capacity(6);
    if imin.x > omin.x {
        boxes.push((DVec3::new(omin.x, omin.y, omin.z), DVec3::X, DVec3::Y, imin.x - omin.x, omax.y - omin.y, omax.z - omin.z));
    }
    if omax.x > imax.x {
        boxes.push((DVec3::new(imax.x, omin.y, omin.z), DVec3::X, DVec3::Y, omax.x - imax.x, omax.y - omin.y, omax.z - omin.z));
    }
    if imin.y > omin.y {
        boxes.push((DVec3::new(imin.x, omin.y, omin.z), DVec3::X, DVec3::Y, imax.x - imin.x, imin.y - omin.y, omax.z - omin.z));
    }
    if omax.y > imax.y {
        boxes.push((DVec3::new(imin.x, imax.y, omin.z), DVec3::X, DVec3::Y, imax.x - imin.x, omax.y - imax.y, omax.z - omin.z));
    }
    if imin.z > omin.z {
        boxes.push((DVec3::new(imin.x, imin.y, omin.z), DVec3::X, DVec3::Y, imax.x - imin.x, imax.y - imin.y, imin.z - omin.z));
    }
    if omax.z > imax.z {
        boxes.push((DVec3::new(imin.x, imin.y, imax.z), DVec3::X, DVec3::Y, imax.x - imin.x, imax.y - imin.y, omax.z - imax.z));
    }
    boxes
}

fn partition_solid_object(
    obj: &rcad_kernel::BRep,
    tools: &[rcad_kernel::BRep],
) -> Result<Vec<rcad_kernel::BRep>, crate::BooleanError> {
    let face_tool_info: Vec<Option<rcad_kernel::geom::Plane>> = tools
        .iter()
        .map(|t| if is_face_like(t) { try_as_planar_face(t) } else { None })
        .collect();
    let has_face_tool = face_tool_info.iter().any(|p| p.is_some());

    let mut expanded_tools: Vec<rcad_kernel::BRep> = Vec::new();
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

    let mut comp_pairs: Vec<(usize, usize)> = Vec::new();
    for (i, comp_i) in expanded_complements.iter().enumerate() {
        if let Some(j) = comp_i {
            if *j > i && expanded_complements.get(*j) == Some(&Some(i)) {
                comp_pairs.push((i, *j));
            }
        }
    }

    for mask in 0..max_mask {
        if comp_pairs.iter().any(|&(i, j)| ((mask >> i) & 1) == ((mask >> j) & 1)) {
            continue;
        }
        let mut cell = obj.clone();
        let mut failed = false;
        let mut first_tool = true;

        for i in 0..n_tools {
            let inside = (mask >> i) & 1 != 0;
            let tool = &expanded_tools[i];

            if inside {
                match crate::bop_occt_ops::boolean_op_generic(
                    crate::BooleanOpType::Intersection, &cell, tool,
                ) {
                    Ok(r) => cell = compact_brep(&r),
                    Err(_) => { failed = true; break; }
                }
            } else if let Some(complement_idx) = expanded_complements[i] {
                let complement = &expanded_tools[complement_idx];
                match crate::bop_occt_ops::boolean_op_generic(
                    crate::BooleanOpType::Intersection, &cell, complement,
                ) {
                    Ok(r) => cell = compact_brep(&r),
                    Err(_) => { failed = true; break; }
                }
            } else if first_tool && is_box_like(tool) {
                if let (Some(tool_bbox), Some(cell_bbox)) = (bounding_box(tool), bounding_box(&cell)) {
                    let comp_boxes = box_complement_of_bbox(&tool_bbox, &cell_bbox);
                    if comp_boxes.is_empty() {
                        cell = rcad_kernel::BRep::new();
                        break;
                    }
                    let cell_solids = extract_solids(&cell);
                    let mut parts = Vec::new();
                    for (origin, u_dir, v_dir, w, h, d) in &comp_boxes {
                        let Ok(comp_box) = rcad_modeling::make_box_brep(*origin, *u_dir, *v_dir, *w, *h, *d) else { continue; };
                        for cell_part in &cell_solids {
                            if let Ok(part) = crate::bop_occt_ops::boolean_op_generic(
                                crate::BooleanOpType::Intersection, cell_part, &comp_box,
                            ) {
                                if count_faces(&part) > 0 {
                                    parts.push(part);
                                }
                            }
                        }
                    }
                    cell = rcad_kernel::BRep::compound_from_shapes(&parts);
                    first_tool = false;
                    continue;
                }
                // fall through to Diff
                match crate::bop_occt_ops::boolean_op_generic(
                    crate::BooleanOpType::Difference, &cell, tool,
                ) {
                    Ok(r) => cell = compact_brep(&r),
                    Err(_) => { failed = true; break; }
                }
            } else {
                match crate::bop_occt_ops::boolean_op_generic(
                    crate::BooleanOpType::Difference, &cell, tool,
                ) {
                    Ok(r) => cell = compact_brep(&r),
                    Err(_) => { failed = true; break; }
                }
            }
            first_tool = false;
        }

        if failed { continue; }

        let face_indices = collect_flat_face_indices(&cell);
        if face_indices.is_empty() { continue; }

        let comps = connected_face_components(&cell, &face_indices);
        for component in comps {
            if !component.is_empty() {
                let subset = extract_brep_subset(&cell, &component);
                cells.push(subset);
            }
        }
    }

    cells.retain(|c| crate::total_volume(c) > 1e-10);
    Ok(cells)
}

fn partition_face_object(
    obj: &rcad_kernel::BRep,
    tools: &[rcad_kernel::BRep],
) -> Result<Vec<rcad_kernel::BRep>, crate::BooleanError> {
    let solid_tools: Vec<rcad_kernel::BRep> = tools.iter().filter(|t| !is_face_like(t)).cloned().collect();
    if solid_tools.is_empty() || count_faces(obj) == 0 {
        return Ok(vec![obj.clone()]);
    }

    let mut cells = Vec::new();
    let mut remaining = obj.clone();
    for tool in &solid_tools {
        let inside = crate::bop_occt_ops::boolean_op_generic(
            crate::BooleanOpType::Intersection, &remaining, tool,
        )?;
        cells.push(remaining.clone());
        remaining = crate::bop_occt_ops::boolean_op_generic(
            crate::BooleanOpType::Difference, &remaining, tool,
        )?;
    }
    cells.push(remaining);
    Ok(cells)
}

fn collect_flat_face_indices(brep: &rcad_kernel::BRep) -> Vec<usize> {
    let mut indices = Vec::new();
    let mut flat_idx = 0;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for _ in &shd.faces {
                        indices.push(flat_idx);
                        flat_idx += 1;
                    }
                }
            }
        }
    }
    indices
}

fn connected_face_components(brep: &rcad_kernel::BRep, face_indices: &[usize]) -> Vec<Vec<usize>> {
    use std::collections::{HashMap, HashSet};
    let face_set: HashSet<usize> = face_indices.iter().copied().collect();
    if face_set.is_empty() { return Vec::new(); }

    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut flat_idx: usize = 0;

    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for lfi in 0..shd.faces.len() {
                        let global_fi = flat_idx + lfi;
                        if face_set.contains(&global_fi) {
                            let face_sr = &shd.faces[lfi];
                            if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                                let collect = |wire_sr: &topods::Shape, e2f: &mut HashMap<usize, Vec<usize>>| {
                                    if let topods::TShape::Wire(wd) = &*brep.tshapes[wire_sr.index] {
                                        for e_sr in &wd.edges { e2f.entry(e_sr.index).or_default().push(global_fi); }
                                    }
                                };
                                collect(&fd.outer_wire, &mut edge_to_faces);
                                for iw_sr in &fd.inner_wires { collect(iw_sr, &mut edge_to_faces); }
                            }
                        }
                    }
                    flat_idx += shd.faces.len();
                }
            }
        }
    }

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

    let mut visited: HashSet<usize> = HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    for &fi in face_indices {
        if !visited.insert(fi) { continue; }
        let mut component: Vec<usize> = Vec::new();
        let mut stack: Vec<usize> = vec![fi];
        while let Some(current) = stack.pop() {
            component.push(current);
            if let Some(neighbors) = adjacency.get(&current) {
                for &neighbor in neighbors {
                    if visited.insert(neighbor) { stack.push(neighbor); }
                }
            }
        }
        components.push(component);
    }
    components
}
