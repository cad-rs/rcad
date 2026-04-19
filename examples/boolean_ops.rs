//! Example: Boolean operations (union, intersection, difference) between boxes.
//!
//! Run: cargo run --example boolean_ops

use glam::DVec3;
use rcad_algorithms::{
    BooleanOpType,
    SimplifyOptions,
    boolean_op_simplified,
    geom_populate::populate_box_geom,
};
use rcad_kernel::BRep;
use rcad_modeling::*;
use rcad_step::writer::{ExportSelection, StepWriteOptions, StepWriter};
use std::collections::HashMap;

fn main() {
    // Keep default simplification for stable boolean geometry.
    let simplify_opts = SimplifyOptions::default();
    let union_simplify_opts = simplify_opts;

    // ── 1. Union of two overlapping boxes ──────────────────────────────
    println!("1. Union of two overlapping boxes");
    let mut a = make_box_at(0.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    let mut b = make_box_at(1.5, 0.5, 0.5, 3.0, 1.0, 1.0);
    populate_box_geom(&mut a);
    populate_box_geom(&mut b);

    let (union_raw, report) =
        boolean_op_simplified(BooleanOpType::Union, &a, &b, union_simplify_opts).expect("union");
    println!(
        "   Simplified: merges={}, internal_removed={}, small_edges_removed={}, wires_fixed={}, vertices_merged={}",
        report.same_domain_face_merges,
        report.internal_faces_removed,
        report.small_edges_removed,
        report.wires_fixed,
        report.vertices_merged,
    );
    inspect_union_intersection_face_edges(&union_raw);
    write_step_with_mode(&union_raw, "output_bool_union.step", false);

    // ── 2. Intersection of two overlapping boxes ───────────────────────
    println!("2. Intersection of two overlapping boxes");
    let (intersection, _) = boolean_op_simplified(BooleanOpType::Intersection, &a, &b, simplify_opts).expect("intersection");
    write_step(&intersection, "output_bool_intersection.step");

    // ── 3. Difference A - B ────────────────────────────────────────────
    println!("3. Difference A - B");
    let (difference, report) = boolean_op_simplified(BooleanOpType::Difference, &a, &b, simplify_opts).expect("difference");
    println!(
        "   Simplified: merges={}, internal_removed={}, small_edges_removed={}, wires_fixed={}, vertices_merged={}",
        report.same_domain_face_merges,
        report.internal_faces_removed,
        report.small_edges_removed,
        report.wires_fixed,
        report.vertices_merged,
    );
    write_step(&difference, "output_bool_difference.step");

    // ── 4. Difference B - A (asymmetric) ──────────────────────────────
    println!("4. Difference B - A");
    let (diff_ba, _) = boolean_op_simplified(BooleanOpType::Difference, &b, &a, simplify_opts).expect("difference B-A");
    write_step(&diff_ba, "output_bool_difference_ba.step");

    // ── 5. Box with a rectangular hole (contained subtraction) ────────
    println!("5. Box with rectangular slot");
    let mut outer = make_box_at(0.0, 0.0, 0.0, 6.0, 4.0, 4.0);
    let mut slot = make_box_at(2.0, 1.0, -0.5, 2.0, 2.0, 5.0);
    populate_box_geom(&mut outer);
    populate_box_geom(&mut slot);

    let (slotted, _) = boolean_op_simplified(BooleanOpType::Difference, &outer, &slot, simplify_opts).expect("slot");
    write_step(&slotted, "output_bool_slot.step");

    // ── 6. Three-box union (chained) ──────────────────────────────────
    println!("6. Three-box cross (chained union)");
    let mut bx = make_box_at(-0.5, -2.0, -0.5, 1.0, 4.0, 1.0);
    let mut by = make_box_at(-2.0, -0.5, -0.5, 4.0, 1.0, 1.0);
    let mut bz = make_box_at(-0.5, -0.5, -2.0, 1.0, 1.0, 4.0);
    populate_box_geom(&mut bx);
    populate_box_geom(&mut by);
    populate_box_geom(&mut bz);

    let (mut cross, _) = boolean_op_simplified(BooleanOpType::Union, &bx, &by, simplify_opts).expect("cross xy");
    // The result of union doesn't have GeomStore populated for further booleans,
    // so we export the two-arm cross and the third arm separately.
    // For full chaining, populate_box_geom would need to be generalized.
    // Instead, combine visually via append_brep:
    rcad_scene::append_brep(&mut cross, bz);
    write_step(&cross, "output_bool_cross.step");

    println!("Exported 6 boolean operation STEP files.");
}

/// Helper: create an axis-aligned box at (x, y, z) with given dimensions.
fn make_box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
    let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, w, h, d).expect("make_box_brep");
    for v in &mut brep.vertices {
        v.point += DVec3::new(x, y, z);
    }
    brep
}

fn write_step(brep: &BRep, path: &str) {
    write_step_with_mode(brep, path, true);
}

fn write_step_with_mode(brep: &BRep, path: &str, gmsh_strict: bool) {
    let options = StepWriteOptions {
        gmsh_strict,
        ..Default::default()
    };
    let step = StepWriter::write_string_with_options(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
        &options,
    );
    std::fs::write(path, step).expect("write STEP file");
    println!("  -> {path}");
}

fn inspect_union_intersection_face_edges(brep: &BRep) {
    let mut patch_faces: Vec<(usize, usize, usize)> = Vec::new();
    let mut nearest: Option<(usize, usize, usize, f64)> = None;

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                let n = face.normal.normalize_or_zero();
                if n.x.abs() < 0.9 {
                    continue;
                }
                let mut xs: Vec<f64> = Vec::with_capacity(face.outer_wire.edges.len());
                for we in &face.outer_wire.edges {
                    if let Some(e) = brep.edges.get(we.idx) {
                        let vi = if we.forward { e.start } else { e.end };
                        if let Some(v) = brep.vertices.get(vi) {
                            xs.push(v.point.x);
                        }
                    }
                }
                if xs.len() < 3 {
                    continue;
                }
                let x_min = xs.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                let x_max = xs.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                let x_avg = xs.iter().sum::<f64>() / xs.len() as f64;

                let score = (x_avg - 3.0).abs();
                match nearest {
                    Some((_, _, _, best)) if score >= best => {}
                    _ => nearest = Some((si, shi, fi, score)),
                }

                if (x_max - x_min) <= 1e-5 && (x_avg - 3.0).abs() <= 1e-5 {
                    patch_faces.push((si, shi, fi));
                }
            }
        }
    }

    if patch_faces.is_empty() {
        if let Some((si, shi, fi, score)) = nearest {
            let face = &brep.solids[si].shells[shi].faces[fi];
            println!(
                "   [diag] 未找到严格x=3补丁；改为分析最近候选 face={}/{}/{} score={:.6} outer_edges={}",
                si,
                shi,
                fi,
                score,
                face.outer_wire.edges.len()
            );
            patch_faces.push((si, shi, fi));
        } else {
            println!("   [diag] 未找到目标相交面补丁（x≈3, normal≈+X）");
            return;
        }
    }

    let mut all_seg_counts: HashMap<String, usize> = HashMap::new();
    let mut raw_outer_edges = 0usize;
    let mut y_line_intervals: HashMap<i64, Vec<(f64, f64)>> = HashMap::new();
    let mut z_line_intervals: HashMap<i64, Vec<(f64, f64)>> = HashMap::new();

    for (si, shi, fi) in &patch_faces {
        let face = &brep.solids[*si].shells[*shi].faces[*fi];
        raw_outer_edges += face.outer_wire.edges.len();
        for we in &face.outer_wire.edges {
            let Some(edge) = brep.edges.get(we.idx) else { continue };
            let (sv, ev) = if we.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            let Some(sp) = brep.vertices.get(sv).map(|v| v.point) else { continue };
            let Some(ep) = brep.vertices.get(ev).map(|v| v.point) else { continue };
            let k = segment_key(sp, ep);
            *all_seg_counts.entry(k).or_insert(0) += 1;

            let dy = (ep.y - sp.y).abs();
            let dz = (ep.z - sp.z).abs();
            if dy >= dz {
                let line_z = ((sp.z + ep.z) * 0.5 * 1e9).round() as i64;
                let a = sp.y.min(ep.y);
                let b = sp.y.max(ep.y);
                y_line_intervals.entry(line_z).or_default().push((a, b));
            } else {
                let line_y = ((sp.y + ep.y) * 0.5 * 1e9).round() as i64;
                let a = sp.z.min(ep.z);
                let b = sp.z.max(ep.z);
                z_line_intervals.entry(line_y).or_default().push((a, b));
            }
        }
    }

    let unique_segments = all_seg_counts.len();
    let repeated_groups = all_seg_counts.values().filter(|&&c| c > 1).count();
    let repeated_total: usize = all_seg_counts
        .values()
        .filter(|&&c| c > 1)
        .map(|&c| c - 1)
        .sum();
    let boundary_segments = all_seg_counts.values().filter(|&&c| c == 1).count();

    let merged_y_edges = y_line_intervals
        .values()
        .map(|itvs| merge_intervals_count(itvs))
        .sum::<usize>();
    let merged_z_edges = z_line_intervals
        .values()
        .map(|itvs| merge_intervals_count(itvs))
        .sum::<usize>();
    let merged_collinear_total = merged_y_edges + merged_z_edges;

    println!("   [diag] Union相交面补丁(x=3): face片数量 = {}", patch_faces.len());
    println!("   [diag] 补丁原始外环边总数 = {raw_outer_edges}");
    println!("   [diag] 几何唯一边段数 = {unique_segments}");
    println!("   [diag] 重复边段组数 = {repeated_groups}, 重复边段总数 = {repeated_total}");
    println!("   [diag] 补丁边界边段数(去内部共享) = {boundary_segments}");
    println!("   [diag] 共线分段合并后边数 = {merged_collinear_total} (期望 8)");
}

fn segment_key(a: DVec3, b: DVec3) -> String {
    let qa = (
        (a.x * 1e9).round() as i64,
        (a.y * 1e9).round() as i64,
        (a.z * 1e9).round() as i64,
    );
    let qb = (
        (b.x * 1e9).round() as i64,
        (b.y * 1e9).round() as i64,
        (b.z * 1e9).round() as i64,
    );
    let (p, q) = if qa <= qb { (qa, qb) } else { (qb, qa) };
    format!(
        "{}:{}:{}|{}:{}:{}",
        p.0, p.1, p.2, q.0, q.1, q.2
    )
}

fn merge_intervals_count(intervals: &[(f64, f64)]) -> usize {
    if intervals.is_empty() {
        return 0;
    }
    let mut v: Vec<(f64, f64)> = intervals
        .iter()
        .map(|&(a, b)| if a <= b { (a, b) } else { (b, a) })
        .collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut merged = 0usize;
    let mut cur_b = v[0].1;
    for (a, b) in v.into_iter().skip(1) {
        if a <= cur_b + 1e-9 {
            if b > cur_b {
                cur_b = b;
            }
        } else {
            merged += 1;
            cur_b = b;
        }
    }
    merged + 1
}
