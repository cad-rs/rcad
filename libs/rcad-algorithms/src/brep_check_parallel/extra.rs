// extra.rs — included inline via `include!("extra.rs")` in mod.rs.
// All module-level imports from mod.rs are available here.
// Helper functions (vpoint, edge_start, edge_end, etc.) defined in mod.rs.

pub fn validate_solids_parallel(brep: &BRep) -> Vec<SolidValidationResult> {
    let solid_indices: Vec<usize> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(i, ts)| {
            if matches!(&**ts as &TShape, TShape::Solid(_)) {
                Some(i)
            } else {
                None
            }
        })
        .collect();
    solid_indices
        .par_iter()
        .map(|&si| validate_single_solid(brep, si))
        .collect()
}

/// Validate a single solid.
fn validate_single_solid(brep: &BRep, si: usize) -> SolidValidationResult {
    use std::collections::HashSet;

    let sd = match &*brep.tshapes[si] {
        TShape::Solid(s) => s,
        _ => panic!("not a solid at index {}", si),
    };

    let n_edges = brep.edge_count();

    let mut errors = Vec::new();
    let mut warnings = Vec::new();

    // Get shell results
    let shell_results: Vec<ShellValidationResult> = sd
        .shells
        .iter()
        .enumerate()
        .map(|(shi, _)| validate_single_shell(brep, si, shi))
        .collect();

    // Aggregate counts
    let face_count: usize = sd
        .shells
        .iter()
        .map(|shell_sr| {
            if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                shd.faces.len()
            } else {
                0
            }
        })
        .sum();

    let edge_count: usize;
    let vertex_count: usize;
    {
        let mut edges: HashSet<usize> = HashSet::new();
        let mut verts: HashSet<usize> = HashSet::new();

        for shell_sr in &sd.shells {
            if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                for face_sr in &shd.faces {
                    for ei in face_edge_refs(brep, *face_sr) {
                        if ei < n_edges {
                            edges.insert(ei);
                            verts.insert(edge_start(brep, ei));
                            verts.insert(edge_end(brep, ei));
                        }
                    }
                }
            }
        }

        edge_count = edges.len();
        vertex_count = verts.len();
    }

    // Compute Euler characteristic
    let euler_characteristic = vertex_count as i64 - edge_count as i64 + face_count as i64;

    // Check if all shells are closed and manifold
    let is_closed = shell_results.iter().all(|s| s.is_closed);
    let is_manifold = shell_results.iter().all(|s| s.is_manifold);
    let orientation_valid = shell_results.iter().all(|s| s.orientation_consistent);

    // Compute volume (approximate using shell volumes)
    let volume: f64 = sd
        .shells
        .iter()
        .map(|shell_sr| compute_shell_volume(brep, *shell_sr))
        .sum();

    let has_positive_volume = volume > 0.0;

    // Compute genus
    let genus = if is_closed && is_manifold {
        let g = (2 - euler_characteristic) / 2;
        if (2 - euler_characteristic) % 2 == 0 && g >= 0 {
            Some(g)
        } else {
            None
        }
    } else {
        None
    };

    // Generate errors
    if !is_closed {
        errors.push("Solid has unclosed shells".to_string());
    }
    if !is_manifold {
        errors.push("Solid has non-manifold topology".to_string());
    }
    if !has_positive_volume {
        warnings.push("Solid has zero or negative volume".to_string());
    }

    let is_valid = errors.is_empty() && shell_results.iter().all(|s| s.is_valid);

    SolidValidationResult {
        solid_idx: si,
        is_valid,
        shell_count: sd.shells.len(),
        face_count,
        edge_count,
        vertex_count,
        euler_characteristic,
        is_closed,
        is_manifold,
        orientation_valid,
        has_positive_volume,
        volume,
        genus,
        shell_results,
        errors,
        warnings,
    }
}

/// Compute the volume of a shell using signed volume method.
fn compute_shell_volume(brep: &BRep, shell_sr: ShapeRef) -> f64 {
    let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
        return 0.0;
    };

    let mut volume = 0.0_f64;

    for face_sr in &shd.faces {
        let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };

        // Compute face normal from surface
        let normal = fd
            .surface
            .as_ref()
            .map(|s| rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0))
            .unwrap_or(DVec3::ZERO);

        // Get vertices of the outer wire
        let mut verts: Vec<DVec3> = Vec::new();
        let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else { continue };
        for esr in &wd.edges {
            let fi = esr.index;
            let vi = if esr.orientation == Orientation::Forward {
                edge_start(brep, fi)
            } else {
                edge_end(brep, fi)
            };
            if vi < brep.vertex_count() {
                verts.push(vpoint(brep, vi));
            }
        }

        // Compute signed volume contribution using triangulation
        if verts.len() >= 3 {
            let origin = verts[0];
            for i in 1..verts.len() - 1 {
                let v1 = verts[i] - origin;
                let v2 = verts[i + 1] - origin;
                let signed_vol = v1.cross(v2).dot(normal) / 6.0;
                volume += signed_vol;
            }
        }
    }

    volume.abs()
}

// ===========================================================?
// Comprehensive Parallel Check
// ===========================================================?

/// Perform a comprehensive parallel check of a brep.
pub fn check_brep_parallel(brep: &BRep, config: &ParallelCheckConfig) -> ParallelCheckReport {
    use std::time::Instant;

    let start_time = Instant::now();
    let mut phase_timings: Vec<CheckPhaseTiming> = Vec::new();
    let mut structural_issues = Vec::new();
    let mut parallel_issues = Vec::new();

    let threads_used = if config.num_threads > 0 {
        config.num_threads
    } else {
        rayon::current_num_threads()
    };

    // Count totals via TShape iteration
    let total_solids = brep.solid_count();
    let total_shells: usize = brep
        .tshapes
        .iter()
        .filter(|ts| matches!(&**ts as &TShape, TShape::Shell(_)))
        .count();
    let total_faces = brep.face_count();
    let total_edges = brep.edge_count();
    let total_vertices = brep.vertex_count();

    let use_parallel = total_faces >= config.parallel_threshold
        || total_edges >= config.parallel_threshold
        || total_vertices >= config.parallel_threshold;

    // Face checking
    let mut face_results = Vec::new();
    if config.check_faces {
        let phase_start = Instant::now();
        face_results = if use_parallel {
            check_faces_parallel(brep, threads_used)
        } else {
            check_faces_sequential(brep)
        };
        phase_timings.push(CheckPhaseTiming {
            phase: "faces".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_faces,
        });

        for fr in &face_results {
            if !fr.is_valid {
                structural_issues.push(CheckIssue::DegenerateFace {
                    solid: fr.solid_idx,
                    shell: fr.shell_idx,
                    face: fr.face_idx,
                });
            }
        }
    }

    // Edge checking
    let mut edge_results = Vec::new();
    if config.check_edges {
        let phase_start = Instant::now();
        edge_results = if use_parallel {
            check_edges_parallel(brep, threads_used)
        } else {
            check_edges_sequential(brep)
        };
        phase_timings.push(CheckPhaseTiming {
            phase: "edges".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_edges,
        });

        for er in &edge_results {
            for issue in &er.issues {
                match issue {
                    EdgeCheckIssue::InvalidVertexIndex { vertex_idx } => {
                        structural_issues.push(CheckIssue::InvalidVertexIndex {
                            edge: er.edge_idx,
                            vertex_idx: *vertex_idx,
                        });
                    }
                    EdgeCheckIssue::NonManifold { face_count } => {
                        structural_issues.push(CheckIssue::NonManifoldEdge {
                            edge_idx: er.edge_idx,
                            face_count: *face_count,
                        });
                    }
                    _ => {}
                }
            }
        }
    }

    // Vertex checking
    if config.check_vertices {
        let phase_start = Instant::now();

        if config.check_finite_vertices {
            for (vidx, ts) in brep.tshapes.iter().enumerate() {
                if let TShape::Vertex(vd) = &**ts {
                    if !vd.point.is_finite() {
                        parallel_issues.push(ParallelCheckIssue::NonFiniteVertex {
                            vertex_idx: vidx,
                        });
                    }
                }
            }
        }

        if config.check_isolated_vertices {
            let n_verts = brep.vertex_count();
            let mut referenced = vec![false; n_verts];
            for ts in &brep.tshapes {
                if let TShape::Edge(ed) = &**ts {
                    if ed.first.index < n_verts {
                        referenced[ed.first.index] = true;
                    }
                    if ed.last.index < n_verts {
                        referenced[ed.last.index] = true;
                    }
                }
            }
            for (vidx, &is_ref) in referenced.iter().enumerate() {
                if !is_ref {
                    parallel_issues.push(ParallelCheckIssue::IsolatedVertex {
                        vertex_idx: vidx,
                    });
                }
            }
        }

        if config.check_duplicate_vertices {
            let vert_data: Vec<DVec3> = brep
                .tshapes
                .iter()
                .filter_map(|ts| {
                    if let TShape::Vertex(vd) = &**ts {
                        Some(vd.point)
                    } else {
                        None
                    }
                })
                .collect();
            let duplicates = find_duplicate_vertices_parallel(&vert_data, config.tolerance);
            parallel_issues.extend(duplicates);
        }

        phase_timings.push(CheckPhaseTiming {
            phase: "vertices".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_vertices,
        });
    }

    // Shell validation
    let mut shell_results = Vec::new();
    if config.check_shells {
        let phase_start = Instant::now();
        shell_results = validate_shells_parallel(brep);
        phase_timings.push(CheckPhaseTiming {
            phase: "shells".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_shells,
        });
    }

    // Solid validation
    let mut solid_results = Vec::new();
    if config.check_solids {
        let phase_start = Instant::now();
        solid_results = validate_solids_parallel(brep);
        phase_timings.push(CheckPhaseTiming {
            phase: "solids".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_solids,
        });
    }

    let total_duration_ms = start_time.elapsed().as_millis() as u64;

    let is_valid = structural_issues.is_empty()
        && parallel_issues.is_empty()
        && shell_results.iter().all(|s| s.is_valid)
        && solid_results.iter().all(|s| s.is_valid);

    let stats = ParallelCheckStats {
        face_count: total_faces,
        edge_count: total_edges,
        vertex_count: total_vertices,
        issue_count: structural_issues.len() + parallel_issues.len(),
        is_valid,
        was_parallel: use_parallel,
        thread_count: threads_used,
    };

    ParallelCheckReport {
        is_valid,
        total_faces,
        total_edges,
        total_vertices,
        total_solids,
        total_shells,
        threads_used,
        was_parallel: use_parallel,
        total_duration_ms,
        phase_timings,
        face_results,
        edge_results,
        shell_results,
        solid_results,
        structural_issues,
        parallel_issues,
        stats,
    }
}

/// Sequential face checking fallback.
fn check_faces_sequential(brep: &BRep) -> Vec<FaceCheckResult> {
    let n_edges = brep.edge_count();
    let tolerance = TOLERANCE_MESH_LEGACY;

    let mut results = Vec::new();
    for (si, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Solid(sd) = &**ts else { continue };
        for (shi, shell_sr) in sd.shells.iter().enumerate() {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
            for fi in 0..shd.faces.len() {
                results.push(check_single_face_detailed(brep, si, shi, fi, n_edges, tolerance));
            }
        }
    }
    results
}

/// Sequential edge checking fallback.
fn check_edges_sequential(brep: &BRep) -> Vec<EdgeCheckResult> {
    let n_verts = brep.vertex_count();
    let tolerance = TOLERANCE_MESH_LEGACY;

    let n_edges = brep.edge_count();
    let mut edge_face_counts = vec![0usize; n_edges];
    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { continue };
            for face_sr in &shd.faces {
                for ei in face_edge_refs(brep, *face_sr) {
                    if ei < n_edges {
                        edge_face_counts[ei] += 1;
                    }
                }
            }
        }
    }

    // Build edges data for iteration
    let edge_data: Vec<(usize, usize, usize)> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(ei, ts)| {
            if let TShape::Edge(ed) = &**ts {
                Some((ei, ed.first.index, ed.last.index))
            } else {
                None
            }
        })
        .collect();

    edge_data
        .iter()
        .map(|&(ei, start, end)| {
            let edge = rcad_kernel::topology::Edge { start, end };
            check_single_edge(brep, ei, &edge, n_verts, edge_face_counts[ei], tolerance)
        })
        .collect()
}

/// Perform parallel check and return detailed statistics.
pub fn check_parallel_with_stats(brep: &BRep) -> (CheckResult, ParallelCheckStats) {
    let face_count = brep.face_count();
    let edge_count = brep.edge_count();
    let vertex_count = brep.vertex_count();

    let options = ParallelCheckOptions::default();
    let result = check_parallel_with_options(brep, &options);

    let stats = ParallelCheckStats {
        face_count,
        edge_count,
        vertex_count,
        issue_count: result.issues.len() + result.parallel_issues.len(),
        is_valid: result.is_valid(),
        was_parallel: result.was_parallel,
        thread_count: result.thread_count,
    };

    (result.to_check_result(), stats)
}
