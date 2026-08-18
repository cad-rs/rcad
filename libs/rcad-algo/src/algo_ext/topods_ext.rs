// =============================================================================
// Topods-native Shell / Solid Extraction + Partition (migration)
// =============================================================================

use rcad_kernel::topods;

/// Remove stale vertices/edges from a topods::BRep, rebuild with only
/// topology-referenced data. Face/shell orientations of the source solids are
/// preserved (OCCT TopExp explode semantics: sub-shapes keep their
/// orientation); a BRep without solids falls back to one FORWARD shell over
/// all face TShapes.
pub fn compact_brep_topods(brep: &topods::BRep) -> topods::BRep {
    let mut shells = Vec::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                    if !shd.faces.is_empty() {
                        shells.push((sr.orientation, shd.faces.clone()));
                    }
                }
            }
        }
    }
    if shells.is_empty() {
        let faces: Vec<topods::Shape> = brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(&***ts, topods::TShape::Face(_)))
            .map(|(i, ts)| topods::Shape {
                data: ts.clone(),
                index: i,
                orientation: topods::Orientation::Forward,
                location: 0,
            })
            .collect();
        if faces.is_empty() {
            return topods::BRep::new();
        }
        shells.push((topods::Orientation::Forward, faces));
    }
    compact_brep_face_subset(brep, &shells)
}

/// Extract each solid from a topods::BRep as a separate self-contained BRep.
/// OCCT `explode ... SOLID`: the returned solid keeps its shells and the face
/// orientations referenced by them.
pub fn extract_solids_topods(brep: &topods::BRep) -> Vec<topods::BRep> {
    let mut groups: Vec<Vec<(topods::Orientation, Vec<topods::Shape>)>> = Vec::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            let mut shells = Vec::new();
            for sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                    if !shd.faces.is_empty() {
                        shells.push((sr.orientation, shd.faces.clone()));
                    }
                }
            }
            if !shells.is_empty() { groups.push(shells); }
        }
    }
    groups.into_iter().map(|sh| compact_brep_face_subset(brep, &sh)).collect()
}

/// Extract each shell from a topods::BRep as a separate self-contained BRep.
pub fn extract_shells_topods(brep: &topods::BRep) -> Vec<topods::BRep> {
    let mut groups: Vec<Vec<(topods::Orientation, Vec<topods::Shape>)>> = Vec::new();
    for ts in &brep.tshapes {
        match &**ts {
            topods::TShape::Solid(sd) => {
                for sr in &sd.shells {
                    if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                        if !shd.faces.is_empty() {
                            groups.push(vec![(sr.orientation, shd.faces.clone())]);
                        }
                    }
                }
            }
            topods::TShape::Shell(shd) => {
                if !shd.faces.is_empty() {
                    groups.push(vec![(topods::Orientation::Forward, shd.faces.clone())]);
                }
            }
            _ => {}
        }
    }
    groups.into_iter().map(|sh| compact_brep_face_subset(brep, &sh)).collect()
}

/// Build a self-contained topods::BRep containing only the specified shells.
/// Each shell carries its face references (with orientations) and the shell
/// orientation; the copies preserve the orientations of the referenced
/// vertices, edges and wires too.
fn compact_brep_face_subset(
    brep: &topods::BRep,
    shells: &[(topods::Orientation, Vec<topods::Shape>)],
) -> topods::BRep {
    use std::collections::{HashMap, HashSet};
    if shells.is_empty() { return topods::BRep::new(); }

    // Collect unique edge tshape indices from face wires
    let mut edge_set: HashSet<usize> = HashSet::new();
    for (_, faces) in shells {
        for fsr in faces {
            if let topods::TShape::Face(fd) = &*brep.tshapes[fsr.index] {
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
            let mut first = v_map[&ed.first.index].clone();
            first.orientation = ed.first.orientation;
            first.location = ed.first.location;
            let mut last = v_map[&ed.last.index].clone();
            last.orientation = ed.last.orientation;
            last.location = ed.last.location;
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
            let edges: Vec<topods::Shape> = wd.edges.iter().map(|sr| {
                let mut e = e_map[&sr.index].clone();
                e.orientation = sr.orientation;
                e.location = sr.location;
                e
            }).collect();
            let mut w = r.add_twire(edges);
            w.orientation = old_wire_sr.orientation;
            w.location = old_wire_sr.location;
            w
        } else {
            old_wire_sr.clone()
        }
    };

    // Add shells with remapped face refs
    let mut shell_srs = Vec::new();
    for (sh_or, faces) in shells {
        let mut face_srs = Vec::new();
        for fsr in faces {
            if let topods::TShape::Face(fd) = &*brep.tshapes[fsr.index] {
                let ow = wire_ref(brep, &mut r, &fd.outer_wire);
                let iw: Vec<topods::Shape> = fd.inner_wires.iter().map(|sr| wire_ref(brep, &mut r, sr)).collect();

                // Build internal_vertices: map old vertex tshape indices to new ShapeRefs
                let iv: Vec<topods::Shape> = fd.internal_vertices.iter()
                    .filter_map(|iv_sr| {
                        let mut s = v_map.get(&iv_sr.index)?.clone();
                        s.orientation = iv_sr.orientation;
                        s.location = iv_sr.location;
                        Some(s)
                    })
                    .collect();

                let mut sr = r.add_tface(fd.surface.clone(), ow, iw, fd.sample_point, fd.uv_domain, iv, fd.natural_restriction);
                sr.orientation = fsr.orientation;
                sr.location = fsr.location;
                r.face_mut(sr.clone()).tolerance = fd.tolerance;
                face_srs.push(sr);
            }
        }
        if !face_srs.is_empty() {
            let mut sh = r.add_tshell(face_srs);
            sh.orientation = *sh_or;
            shell_srs.push(sh);
        }
    }

    // Wrap in Shell -> Solid
    if !shell_srs.is_empty() {
        r.add_tsolid(shell_srs);
    }

    r
}
