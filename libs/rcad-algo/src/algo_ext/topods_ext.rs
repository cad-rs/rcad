// =============================================================================
// Topods-native Shell / Solid Extraction + Partition (migration)
// =============================================================================

use rcad_kernel::topods;

/// Remove stale vertices/edges from a topods::BRep, rebuild with only
/// topology-referenced data.
pub fn compact_brep_topods(brep: &topods::BRep) -> topods::BRep {
    // Collect all face tshape indices
    let face_indices: Vec<usize> = brep.tshapes.iter().enumerate()
        .filter(|(_, ts)| matches!(&***ts, topods::TShape::Face(_)))
        .map(|(i, _)| i)
        .collect();
    compact_brep_face_subset(brep, &face_indices)
}

/// Extract each solid from a topods::BRep as a separate self-contained BRep.
pub fn extract_solids_topods(brep: &topods::BRep) -> Vec<topods::BRep> {
    let mut groups = Vec::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            let mut faces = Vec::new();
            for sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                    for fsr in &shd.faces { faces.push(fsr.index); }
                }
            }
            if !faces.is_empty() { groups.push(faces); }
        }
    }
    groups.into_iter().map(|f| compact_brep_face_subset(brep, &f)).collect()
}

/// Extract each shell from a topods::BRep as a separate self-contained BRep.
pub fn extract_shells_topods(brep: &topods::BRep) -> Vec<topods::BRep> {
    let mut groups = Vec::new();
    for ts in &brep.tshapes {
        let faces = match &**ts {
            topods::TShape::Solid(sd) => {
                let mut all = Vec::new();
                for sr in &sd.shells {
                    if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                        let f: Vec<usize> = shd.faces.iter().map(|fsr| fsr.index).collect();
                        if !f.is_empty() { all.push(f); }
                    }
                }
                all
            }
            topods::TShape::Shell(shd) => {
                let f: Vec<usize> = shd.faces.iter().map(|fsr| fsr.index).collect();
                if !f.is_empty() { vec![f] } else { vec![] }
            }
            _ => vec![],
        };
        groups.extend(faces);
    }
    groups.into_iter().map(|f| compact_brep_face_subset(brep, &f)).collect()
}

/// Build a self-contained topods::BRep containing only the specified face
/// tshape indices. Copies only referenced edges/vertices/geometry.
fn compact_brep_face_subset(brep: &topods::BRep, face_indices: &[usize]) -> topods::BRep {
    use std::collections::{HashMap, HashSet};
    if face_indices.is_empty() { return topods::BRep::new(); }

    // Collect unique edge tshape indices from face wires
    let mut edge_set: HashSet<usize> = HashSet::new();
    for &fi in face_indices {
        if let topods::TShape::Face(fd) = &*brep.tshapes[fi] {
            if let topods::TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                for sr in &wd.edges { edge_set.insert(sr.index); }
            }
            for isr in &fd.inner_wires {
                if let topods::TShape::Wire(wd) = &*brep.tshapes[isr.index] {
                    for sr in &wd.edges { edge_set.insert(sr.index); }
                }
            }
        }
    }

    // Collect vertex tshape indices from edges
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &ei in &edge_set {
        if let topods::TShape::Edge(ed) = &*brep.tshapes[ei] {
            vertex_set.insert(ed.first.index);
            vertex_set.insert(ed.last.index);
        }
    }

    // Sort for deterministic output
    let mut sorted_v: Vec<usize> = vertex_set.iter().copied().collect(); sorted_v.sort();
    let mut sorted_e: Vec<usize> = edge_set.iter().copied().collect(); sorted_e.sort();

    // Old tshape index -> new Shape
    let mut v_map: HashMap<usize, topods::Shape> = HashMap::new();
    let mut e_map: HashMap<usize, topods::Shape> = HashMap::new();

    let mut r = topods::BRep::new();

    // Add vertices
    for &old in &sorted_v {
        if let topods::TShape::Vertex(vd) = &*brep.tshapes[old] {
            let sr = r.add_tvertex(vd.point);
            r.vertex_mut(sr.clone()).tolerance = vd.tolerance;
            v_map.insert(old, sr);
        }
    }

    // Add edges with remapped vertex refs
    for &old in &sorted_e {
        if let topods::TShape::Edge(ed) = &*brep.tshapes[old] {
            let first = v_map[&ed.first.index].clone();
            let last = v_map[&ed.last.index].clone();
            let sr = r.add_tedge(ed.curve.clone(), first, last, ed.range);
            let em = r.edge_mut(sr.clone());
            em.tolerance = ed.tolerance;
            em.representations = ed.representations.clone();
            em.pcurves = ed.pcurves.clone();
            em.degenerated = ed.degenerated;
            em.same_parameter = ed.same_parameter;
            em.same_range = ed.same_range;
            em.vertex_params = ed.vertex_params.clone();
            e_map.insert(old, sr);
        }
    }

    // Build wire refs helper
    let wire_ref = |brep: &topods::BRep, r: &mut topods::BRep, old_wire_sr: &topods::Shape| -> topods::Shape {
        if let topods::TShape::Wire(wd) = &*brep.tshapes[old_wire_sr.index] {
            let edges: Vec<topods::Shape> = wd.edges.iter().map(|sr| e_map[&sr.index].clone()).collect();
            r.add_twire(edges)
        } else {
            old_wire_sr.clone()
        }
    };

    // Add faces
    let mut face_srs = Vec::new();
    for &fi in face_indices {
        if let topods::TShape::Face(fd) = &*brep.tshapes[fi] {
            let ow = wire_ref(brep, &mut r, &fd.outer_wire);
            let iw: Vec<topods::Shape> = fd.inner_wires.iter().map(|sr| wire_ref(brep, &mut r, sr)).collect();

            // Build internal_vertices: map old vertex tshape indices to new ShapeRefs
            let iv: Vec<topods::Shape> = fd.internal_vertices.iter()
                .filter_map(|iv_sr| v_map.get(&iv_sr.index).cloned())
                .collect();

            let sr = r.add_tface(fd.surface.clone(), ow, iw, fd.sample_point, fd.uv_domain, iv, fd.natural_restriction);
            r.face_mut(sr.clone()).tolerance = fd.tolerance;
            face_srs.push(sr);
        }
    }

    // Wrap in Shell -> Solid
    if !face_srs.is_empty() {
        let shell = r.add_tshell(face_srs);
        r.add_tsolid(vec![shell]);
    }

    r
}
