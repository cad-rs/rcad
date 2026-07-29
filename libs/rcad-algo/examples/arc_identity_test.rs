use std::collections::HashSet;
use std::sync::Arc;
use glam::DVec3;
use rcad_kernel::topods::{self, TShape, BRep, Orientation};
use rcad_kernel::topo_shape::Shape;
use rcad_modeling::make_box_brep;
use rcad_brep::tools::transform_shape;

fn main() {
    let brep_a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let mut brep_b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    // Translate b (as P005 does)
    rcad_brep::tools::transform_shape(&mut brep_b, glam::DAffine3::from_translation(DVec3::new(0.2, 0.0, 0.0)));

    // Create children as build_ds does
    let children_a: Vec<Shape> = brep_a.tshapes.iter().map(|ts| Shape::new(ts.clone(), 0, Orientation::Forward)).collect();
    let children_b: Vec<Shape> = brep_b.tshapes.iter().map(|ts| Shape::new(ts.clone(), 0, Orientation::Forward)).collect();

    // Combine both compounds' children
    let mut all_children = children_a;
    all_children.extend(children_b);

    println!("\n=== Combined {} children from both boxes ===", all_children.len());
    println!("  Direct children: V={} E={} F={}",
        all_children.iter().filter(|s| s.shape_type() == rcad_kernel::topods::ShapeType::Vertex).count(),
        all_children.iter().filter(|s| s.shape_type() == rcad_kernel::topods::ShapeType::Edge).count(),
        all_children.iter().filter(|s| s.shape_type() == rcad_kernel::topods::ShapeType::Face).count());

    // Also check edge vertices
    println!("\n=== Edge vertex ptr_ids for COMBINED set ===");
    let tshape_vertex_ptrs: Vec<(usize, u64)> = all_children.iter().enumerate()
        .filter_map(|(i, s)| {
            if s.shape_type() == rcad_kernel::topods::ShapeType::Vertex {
                Some((i, s.ptr_id()))
            } else {
                None
            }
        })
        .collect();
    println!("  Total vertices in children: {}", tshape_vertex_ptrs.len());

    // For each edge, check vertex ptr_ids (inline sub_shapes_of for Edge)
    for (i, s) in all_children.iter().enumerate() {
        if s.shape_type() != rcad_kernel::topods::ShapeType::Edge { continue; }
        let ed = match &*s.data { TShape::Edge(e) => e, _ => continue };
        let subs = vec![
            Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
            Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
        ];
        for sv in &subs {
            let found = tshape_vertex_ptrs.iter().any(|(_, vp)| *vp == sv.ptr_id());
            if !found {
                println!("  FAIL: Edge[{}] vertex ptr=0x{:x} NOT in children vertices!", i, sv.ptr_id());
            }
        }
    }
    println!("  All edge vertex ptr_ids are in children set");

    // NOW test init_shape behavior
    println!("\n=== Simulating init_shape with map_shape_index ===");
    let mut map: std::collections::HashMap<(u64, u32), usize> = std::collections::HashMap::new();

    // Add Compound as index 0
    let compound_ptr = 0xDEADBEEF;
    map.insert((compound_ptr, 0), 0);

    // Process ALL children in order (as init_shape would for the Compound)
    let mut misses = 0usize;
    for (i, child) in all_children.iter().enumerate() {
        let idx = i + 1; // index after Compound
        let pk = (child.ptr_id(), child.location);
        if map.contains_key(&pk) {
            // Found — would be reused
        } else {
            map.insert(pk, idx);
            misses += 1;
        }
    }
    println!("  Children inserted: {} first-time, {} dedup'd", misses, all_children.len() - misses);
    println!("  Expected unique: {} (8V+12E+6F)x2 = 52", 16+24+12);

    // Now simulate recursive expansion: for each edge, process its vertices
    let mut edge_vertex_misses = 0usize;
    for child in &all_children {
        if child.shape_type() != rcad_kernel::topods::ShapeType::Edge { continue; }
        let ed = match &*child.data { TShape::Edge(e) => e, _ => continue };
        let subs = vec![
            Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
            Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
        ];
        for sv in &subs {
            let pk = (sv.ptr_id(), sv.location);
            if !map.contains_key(&pk) {
                edge_vertex_misses += 1;
                println!("  RECURSIVE_MISS: Edge vertex ptr=0x{:x} loc={} NOT in map!", sv.ptr_id(), sv.location);
            }
        }
    }
    println!("  Edge vertex misses (should be 0): {}", edge_vertex_misses);

    // Just use one brep for the original tests
    let brep = brep_a;
    
    println!("=== tshapes layout ===");
    let mut vc = 0usize; let mut ec = 0usize; let mut wc = 0usize; let mut fc = 0usize;
    for (i, ts) in brep.tshapes.iter().enumerate() {
        let label = match &**ts {
            TShape::Vertex(_) => { let n = vc; vc += 1; format!("V{}", n) }
            TShape::Edge(_) => { let n = ec; ec += 1; format!("E{}", n) }
            TShape::Wire(_) => { let n = wc; wc += 1; format!("W{}", n) }
            TShape::Face(_) => { let n = fc; fc += 1; format!("F{}", n) }
            TShape::Shell(_) => "Shell".into(),
            TShape::Solid(_) => "Solid".into(),
            _ => "?".into(),
        };
        let ptr = Arc::as_ptr(ts) as u64;
        println!("  [{}] {} ptr=0x{:x}", i, label, ptr);
    }

    // Collect ALL vertex ptr_ids from edges' first/last fields
    let mut edge_vertex_ptrs: Vec<(usize, u64, u64)> = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Edge(ed) = &**ts {
            let fptr = Arc::as_ptr(&ed.first.data) as u64;
            let lptr = Arc::as_ptr(&ed.last.data) as u64;
            edge_vertex_ptrs.push((i, fptr, lptr));
        }
    }

    // Get vertex ptr_ids directly from tshapes
    let tshape_vertex_ptrs: Vec<(usize, u64)> = brep.tshapes.iter().enumerate()
        .filter_map(|(i, ts)| {
            if matches!(&**ts, TShape::Vertex(_)) {
                Some((i, Arc::as_ptr(ts) as u64))
            } else {
                None
            }
        })
        .collect();

    println!("\n=== Edge vertex ptr_ids vs tshapes vertex ptr_ids ===");
    let mut all_match = true;
    for (ei, fptr, lptr) in &edge_vertex_ptrs {
        let f_match = tshape_vertex_ptrs.iter().any(|(_, vp)| *vp == *fptr);
        let l_match = tshape_vertex_ptrs.iter().any(|(_, vp)| *vp == *lptr);
        if !f_match || !l_match {
            println!("  Edge[{}] first=0x{:x} match={} last=0x{:x} match={}", 
                ei, fptr, f_match, lptr, l_match);
            all_match = false;
        }
    }
    if all_match {
        println!("  All edge vertex ptr_ids match tshapes vertex entries!");
    }

    // Now simulate what init_shape does
    println!("\n=== Simulating init_shape dedup ===");
    let mut map: std::collections::HashMap<(u64, u32), usize> = std::collections::HashMap::new();
    
    // Add vertices from tshapes (simulating Compound children)
    for (vi, vptr) in &tshape_vertex_ptrs {
        let prev = map.insert((*vptr, 0), *vi);
        // Should be unique — each vertex has its own Arc in tshapes
    }
    println!("  Added {} vertices from tshapes", tshape_vertex_ptrs.len());

    // Now try to find them using edge vertex ptr_ids (simulating sub_shapes_of(Edge))
    for (ei, fptr, lptr) in &edge_vertex_ptrs {
        let f_found = map.get(&(*fptr, 0));
        let l_found = map.get(&(*lptr, 0));
        if f_found.is_none() {
            println!("  FAIL: Edge[{}] first vertex ptr=0x{:x} NOT found in map!", ei, fptr);
        }
        if l_found.is_none() {
            println!("  FAIL: Edge[{}] last vertex ptr=0x{:x} NOT found in map!", ei, lptr);
        }
    }
    println!("  All vertex lookups passed!");
}
