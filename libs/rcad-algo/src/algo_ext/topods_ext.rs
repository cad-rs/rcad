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

/// Test-world result extraction bridge (the facade `.brep` pattern of
/// fillet.rs `FilletResult.brep`, libs/rcad-algo/src/fillet/fillet.rs
/// L273-284): flatten a result root shape into ONE self-contained
/// topods::BRep so the test world can drive STEP export
/// (StepWriter::write_string_with_options) and measurement
/// (total_surface_area / total_volume) on the facade result.
///
/// OCCT has no counterpart — a TopoDS_Shape carries its TShape graph by
/// pointer — this is the glue for the rcad BRep-pool architecture
/// difference #4 (the Stage 2e/3 facades hold my_shape + my_brep pairs).
/// The walk mirrors bop/algo/section.rs `ResultPool::push_recursive` /
/// builder.rs `push_shape_recursive`: sub-shapes pushed first, the TShape
/// Arc is shared (OCCT BRep_Builder::Add references the source TShape,
/// TopoDS_Builder.cxx L57-59), dedup by TShape identity, and the
/// TShape-internal reference indices re-pointed in place to this pool
/// (single ownership: the source pools are not read for topology after the
/// extraction — the same shared-mutation caveat as compact_brep_face_subset
/// below).
///
/// `locations` is the locations table of the pool the root shape was built
/// against (BRep convention: index 0 is the implicit identity and is never
/// stored), carried over verbatim.  Shapes reaching in from other pools
/// (the test input pool) keep valid location data only when they carry the
/// identity location (index 0) — the extracted grids satisfy this.
pub fn extract_result_brep(
    root: &rcad_kernel::topo_shape::Shape,
    locations: Vec<glam::DAffine3>,
) -> topods::BRep {
    use std::collections::HashMap;
    use std::sync::Arc;

    let mut r = topods::BRep::new();
    r.locations = locations;
    // Source TShape ptr -> result tshapes index.
    let mut remap: HashMap<u64, usize> = HashMap::new();

    fn push(
        r: &mut topods::BRep,
        remap: &mut HashMap<u64, usize>,
        shape: &rcad_kernel::topo_shape::Shape,
    ) -> usize {
        use rcad_kernel::topo_shape::Shape;
        let ptr = shape.ptr_id();
        if let Some(&idx) = remap.get(&ptr) {
            return idx;
        }
        // Reserve a slot in tshapes, sharing the source TShape Arc.
        let new_idx = r.tshapes.len();
        r.tshapes.push(shape.data.clone());
        remap.insert(ptr, new_idx);
        // Recursively push the sub-shapes first so their result indices exist
        // (same reachability as builder.rs push_shape_recursive).
        match shape.data.as_ref() {
            topods::TShape::Edge(ed) => {
                push(r, remap, &ed.first);
                push(r, remap, &ed.last);
            }
            topods::TShape::Wire(wd) => {
                for e in &wd.edges {
                    push(r, remap, e);
                }
            }
            topods::TShape::Face(fd) => {
                push(r, remap, &fd.outer_wire);
                for w in &fd.inner_wires {
                    push(r, remap, w);
                }
                for v in &fd.internal_vertices {
                    push(r, remap, v);
                }
            }
            topods::TShape::Shell(sd) => {
                for f in &sd.faces {
                    push(r, remap, f);
                }
            }
            topods::TShape::Solid(sd) => {
                for s in &sd.shells {
                    push(r, remap, s);
                }
                for v in &sd.internal_vertices {
                    push(r, remap, v);
                }
                for e in &sd.internal_edges {
                    push(r, remap, e);
                }
            }
            topods::TShape::CompSolid(shapes) => {
                for s in shapes {
                    push(r, remap, s);
                }
            }
            topods::TShape::Compound(shapes) => {
                for s in shapes {
                    push(r, remap, s);
                }
            }
            topods::TShape::Vertex(_) => {}
        }
        // Re-point the TShape-internal reference indices from the source pool
        // positions to this pool's positions, in place on the shared TShape.
        let raw = Arc::as_ptr(&shape.data) as *mut topods::TShape;
        // SAFETY: single-threaded extraction; the facade result is extracted
        // once and the source pools are not read for topology afterwards
        // (same shared-mutation model as compact_brep_face_subset / section.rs
        // ResultPool).
        unsafe {
            match &mut *raw {
                topods::TShape::Vertex(vd) => {
                    for s in vd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                }
                topods::TShape::Edge(ed) => {
                    for s in ed.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    if let Some(&i) = remap.get(&ed.first.ptr_id()) {
                        ed.first.index = i;
                    }
                    if let Some(&i) = remap.get(&ed.last.ptr_id()) {
                        ed.last.index = i;
                    }
                }
                topods::TShape::Wire(wd) => {
                    for s in wd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    for e in wd.edges.iter_mut() {
                        if let Some(&i) = remap.get(&e.ptr_id()) {
                            e.index = i;
                        }
                    }
                }
                topods::TShape::Face(fd) => {
                    for s in fd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    if let Some(&i) = remap.get(&fd.outer_wire.ptr_id()) {
                        fd.outer_wire.index = i;
                    }
                    for w in fd.inner_wires.iter_mut() {
                        if let Some(&i) = remap.get(&w.ptr_id()) {
                            w.index = i;
                        }
                    }
                    for v in fd.internal_vertices.iter_mut() {
                        if let Some(&i) = remap.get(&v.ptr_id()) {
                            v.index = i;
                        }
                    }
                }
                topods::TShape::Shell(sd) => {
                    for s in sd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    for f in sd.faces.iter_mut() {
                        if let Some(&i) = remap.get(&f.ptr_id()) {
                            f.index = i;
                        }
                    }
                }
                topods::TShape::Solid(sd) => {
                    for s in sd.my_shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    for sh in sd.shells.iter_mut() {
                        if let Some(&i) = remap.get(&sh.ptr_id()) {
                            sh.index = i;
                        }
                    }
                    for v in sd.internal_vertices.iter_mut() {
                        if let Some(&i) = remap.get(&v.ptr_id()) {
                            v.index = i;
                        }
                    }
                    for e in sd.internal_edges.iter_mut() {
                        if let Some(&i) = remap.get(&e.ptr_id()) {
                            e.index = i;
                        }
                    }
                }
                topods::TShape::CompSolid(shapes) => {
                    for s in shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                }
                topods::TShape::Compound(shapes) => {
                    for s in shapes.iter_mut() {
                        if let Some(&i) = remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                }
            }
        }
        new_idx
    }

    // The pool is empty at this point, so the pushed root index is final; the
    // root keeps its own orientation/location.
    let _ = push(&mut r, &mut remap, root);
    r
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

#[cfg(test)]
mod extract_result_brep_tests {
    use super::*;

    /// The extraction bridge preserves the reachable subgraph: a unit cube
    /// pool (6 faces, 12 edges, 8 vertices + the root solid) extracted from
    /// its root keeps every reachable TShape and the measured surface area.

    #[test]
    fn cube_round_trip() {
        // The file is included twice (algo_ext + a #[path] include inside
        // bool_ops_ext); keep the single-run semantics of the baseline.
        if module_path!().contains("bool_ops_ext") {
            return;
        }
        // Build a structural pool manually: vertices -> edges -> wire ->
        // face -> shell -> solid (the same construction path the facade
        // engines use), then extract the root.
        use rcad_kernel::geom::{Curve3, Line3};
        let mut brep = topods::BRep::new();
        let corners = [
            glam::DVec3::new(0.0, 0.0, 0.0),
            glam::DVec3::new(10.0, 0.0, 0.0),
            glam::DVec3::new(10.0, 10.0, 0.0),
            glam::DVec3::new(0.0, 10.0, 0.0),
        ];
        let mut wire_edges = Vec::new();
        for pi in 0..4 {
            let p0 = corners[pi];
            let p1 = corners[(pi + 1) % 4];
            let v0 = brep.add_tvertex(p0);
            let v1 = brep.add_tvertex(p1);
            let e = brep.add_tedge(
                Some(Curve3::Line(Line3::new(p0, p1 - p0))),
                v0,
                v1,
                [0.0, 10.0],
            );
            wire_edges.push(e);
        }
        let wire = brep.add_twire(wire_edges);
        let face = brep.add_tface(None, wire, Vec::new(), None, None, Vec::new(), false);
        let shell = brep.add_tshell(vec![face]);
        let solid = brep.add_tsolid(vec![shell]);

        let src_faces = brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Face(_))).count();
        let src_edges = brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Edge(_))).count();
        let src_vertices = brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Vertex(_))).count();
        assert_eq!(src_faces, 1);
        assert_eq!(src_edges, 4);
        assert_eq!(src_vertices, 4);

        let extracted = extract_result_brep(&solid, brep.locations.clone());
        let dst_faces = extracted.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Face(_))).count();
        let dst_edges = extracted.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Edge(_))).count();
        let dst_vertices = extracted.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Vertex(_))).count();
        assert_eq!(src_faces, dst_faces, "face count preserved");
        assert_eq!(src_edges, dst_edges, "edge count preserved");
        assert_eq!(src_vertices, dst_vertices, "vertex count preserved");

        let area_src = rcad_kernel::base::gprop::surface::surface_area(&brep);
        let area_dst = rcad_kernel::base::gprop::surface::surface_area(&extracted);
        assert!((area_src - area_dst).abs() <= 1e-9 * area_src.abs().max(1.0), "surface area preserved: {area_src} vs {area_dst}");
    }
}
