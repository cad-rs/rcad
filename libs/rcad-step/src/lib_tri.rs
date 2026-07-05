fn collect_shell_faces(parsed: &ParsedStep) -> Vec<Vec<u64>> {
    let mut shells = Vec::new();

    if !parsed.manifold_solids.is_empty() {
        for shell_id in &parsed.manifold_solids {
            if let Some(face_ids) = parsed.closed_shells.get(shell_id) {
                shells.push(face_ids.clone());
            }
        }
    }

    // BREP_WITH_VOIDS: insert all shells (outer + voids) in order so that
    // the caller can detect multi-shell groupings from parsed.brep_with_voids.
    if !parsed.brep_with_voids.is_empty() {
        for (_, (outer, voids)) in &parsed.brep_with_voids {
            if let Some(face_ids) = parsed.closed_shells.get(outer) {
                shells.push(face_ids.clone());
            }
            for void_ref in voids {
                if let Some(face_ids) = parsed.closed_shells.get(void_ref) {
                    shells.push(face_ids.clone());
                }
            }
        }
    }

    if !parsed.shell_based_surface_models.is_empty() {
        for shell_refs in &parsed.shell_based_surface_models {
            for shell_id in shell_refs {
                if let Some(face_ids) = parsed.open_shells.get(shell_id) {
                    shells.push(face_ids.clone());
                } else if let Some(face_ids) = parsed.closed_shells.get(shell_id) {
                    shells.push(face_ids.clone());
                }
            }
        }
    }

    if shells.is_empty() && !parsed.closed_shells.is_empty() {
        for face_ids in parsed.closed_shells.values() {
            shells.push(face_ids.clone());
        }
    }

    if shells.is_empty() && !parsed.open_shells.is_empty() {
        for face_ids in parsed.open_shells.values() {
            shells.push(face_ids.clone());
        }
    }

    if shells.is_empty() && !parsed.advanced_faces.is_empty() {
        let mut all_faces: Vec<u64> = parsed.advanced_faces.keys().copied().collect();
        all_faces.sort_unstable();
        shells.push(all_faces);
    }

    shells
}

fn collect_used_vertices(
    parsed: &ParsedStep,
    shell_face_sets: &[Vec<u64>],
) -> Result<BTreeSet<u64>, StepError> {
    let mut used = BTreeSet::new();

    for shell in shell_face_sets {
        for face_id in shell {
            let bound_ids = parsed
                .advanced_faces
                .get(face_id)
                .ok_or(StepError::MissingEntity {
                    entity_type: "ADVANCED_FACE",
                    id: Some(*face_id),
                })?;
            for bound_id in &bound_ids.bounds {
                let (loop_id, _) = parsed
                    .face_bounds
                    .get(bound_id)
                    .ok_or(StepError::MissingEntity {
                        entity_type: "FACE_BOUND",
                        id: Some(*bound_id),
                    })?;
                if let Some(oriented_ids) = parsed.edge_loops.get(loop_id) {
                    for oriented_id in oriented_ids {
                        let (edge_curve_id, _) =
                            parsed
                                .oriented_edges
                                .get(oriented_id)
                                .ok_or(StepError::MissingEntity {
                                    entity_type: "ORIENTED_EDGE",
                                    id: Some(*oriented_id),
                                })?;
                        let (start, end, _, _) =
                            parsed
                                .edge_curves
                                .get(edge_curve_id)
                                .ok_or(StepError::MissingEntity {
                                    entity_type: "EDGE_CURVE",
                                    id: Some(*edge_curve_id),
                                })?;
                        used.insert(*start);
                        used.insert(*end);
                    }
                } else if let Some(vp_id) = parsed.vertex_loops.get(loop_id) {
                    used.insert(*vp_id);
                } else {
                    return Err(StepError::MissingEntity {
                        entity_type: "EDGE_LOOP",
                        id: Some(*loop_id),
                    });
                }
            }
        }
    }

    Ok(used)
}

