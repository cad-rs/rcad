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
///
/// The TShape Arcs are SHARED with the source BRep (OCCT TopoDS reference
/// semantics 閳?extracting sub-shapes never clones the TShapes), only the
/// Shape.index fields are re-pointed from the source array positions to the
/// result array positions, in place on the shared TShapes.  Sharing keeps the
/// edge pcurves (keyed by the owning face TShape pointer) addressable from
/// the extracted BRep: a rebuild (new TShapes) would orphan every pcurve row
/// whose owner face stayed in the source BRep, breaking BRepGProp face
/// integration and later boolean steps on the extracted solid.
fn compact_brep_face_subset(
    brep: &topods::BRep,
    shells: &[(topods::Orientation, Vec<topods::Shape>)],
) -> topods::BRep {
    if shells.is_empty() { return topods::BRep::new(); }

    use std::collections::HashMap;
    let mut r = topods::BRep::new();
    r.locations = brep.locations.clone();
    // Source tshape index -> result tshapes index (Arc shared, index re-pointed).
    let mut remap: HashMap<usize, usize> = HashMap::new();

    // Push a source Shape (by its tshapes index) into the result, sharing the
    // TShape Arc; sub-shapes are pushed first so their result indices exist.
    // The returned Shape carries the caller's orientation/location.
    fn push(
        brep: &topods::BRep,
        r: &mut topods::BRep,
        remap: &mut HashMap<usize, usize>,
        src: usize,
    ) -> Option<topods::Shape> {
        use topods::Orientation;
        if src >= brep.tshapes.len() { return None; }
        if let Some(&i) = remap.get(&src) {
            return Some(topods::Shape::from_parts(
                r.tshapes[i].clone(),
                i,
                0,
                Orientation::Forward,
            ));
        }
        let idx = r.tshapes.len();
        r.tshapes.push(brep.tshapes[src].clone());
        remap.insert(src, idx);
        let base = topods::Shape::from_parts(
            r.tshapes[idx].clone(),
            idx,
            0,
            Orientation::Forward,
        );
        // Recurse into the sub-shapes (same reachability as
        // builder::push_shape_recursive).
        match &*brep.tshapes[src] {
            topods::TShape::Edge(ed) => {
                push(brep, r, remap, ed.first.index);
                push(brep, r, remap, ed.last.index);
            }
            topods::TShape::Wire(wd) => {
                for e in &wd.edges { push(brep, r, remap, e.index); }
            }
            topods::TShape::Face(fd) => {
                push(brep, r, remap, fd.outer_wire.index);
                for w in &fd.inner_wires { push(brep, r, remap, w.index); }
                for v in &fd.internal_vertices { push(brep, r, remap, v.index); }
            }
            topods::TShape::Shell(sd) => {
                for f in &sd.faces { push(brep, r, remap, f.index); }
            }
            topods::TShape::Solid(sd) => {
                for s in &sd.shells { push(brep, r, remap, s.index); }
                for v in &sd.internal_vertices { push(brep, r, remap, v.index); }
                for e in &sd.internal_edges { push(brep, r, remap, e.index); }
            }
            _ => {}
        }
        // Re-point the shared TShape's internal reference indices from the
        // source array positions to the result array positions (single
        // ownership: the extracted BReps are not read while another extraction
        // mutates the same Arc).
        let raw = std::sync::Arc::as_ptr(&brep.tshapes[src]) as *mut topods::TShape;
        unsafe {
            match &mut *raw {
                topods::TShape::Vertex(vd) => {
                    for s in vd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                }
                topods::TShape::Edge(ed) => {
                    for s in ed.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                    if let Some(&i) = remap.get(&ed.first.index) { ed.first.index = i; }
                    if let Some(&i) = remap.get(&ed.last.index) { ed.last.index = i; }
                }
                topods::TShape::Wire(wd) => {
                    for s in wd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                    for e in wd.edges.iter_mut() {
                        if let Some(&i) = remap.get(&e.index) { e.index = i; }
                    }
                }
                topods::TShape::Face(fd) => {
                    for s in fd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                    if let Some(&i) = remap.get(&fd.outer_wire.index) { fd.outer_wire.index = i; }
                    for w in fd.inner_wires.iter_mut() {
                        if let Some(&i) = remap.get(&w.index) { w.index = i; }
                    }
                    for v in fd.internal_vertices.iter_mut() {
                        if let Some(&i) = remap.get(&v.index) { v.index = i; }
                    }
                }
                topods::TShape::Shell(sd) => {
                    for s in sd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                    for f in sd.faces.iter_mut() {
                        if let Some(&i) = remap.get(&f.index) { f.index = i; }
                    }
                }
                topods::TShape::Solid(sd) => {
                    for s in sd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                    for s in sd.shells.iter_mut() {
                        if let Some(&i) = remap.get(&s.index) { s.index = i; }
                    }
                }
                _ => {}
            }
        }
        Some(base)
    }

    // Rebuild each shell over the shared faces; the shell/solid containers are
    // new (the source shell/solid references were decomposed into the
    // (orientation, faces) groups by the caller).
    let mut shell_srs = Vec::new();
    for (sh_or, faces) in shells {
        let mut face_srs = Vec::new();
        for fsr in faces {
            let Some(base) = push(brep, &mut r, &mut remap, fsr.index) else { continue };
            let mut f = base;
            f.orientation = fsr.orientation;
            f.location = fsr.location;
            face_srs.push(f);
        }
        if !face_srs.is_empty() {
            let mut sh = r.add_tshell(face_srs);
            sh.orientation = *sh_or;
            shell_srs.push(sh);
        }
    }

    // Wrap in Shell -> Solid.
    if !shell_srs.is_empty() {
        r.add_tsolid(shell_srs);
    }

    r
}
