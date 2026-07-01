
/// (fan-triangulation from outer-wire vertices) is <= `max_area`.
///
/// Returns a sorted, deduplicated list of local face indices.
///
/// Note: the area estimate is a polygon fan-triangulation; it is exact for
/// planar convex faces and an approximation for curved faces.
pub fn identify_small_faces(brep: &BRep, max_area: f64) -> Vec<usize> {
    if max_area <= 0.0 {
        return Vec::new();
    }
    let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
        return Vec::new();
    };

    let mut result = Vec::new();

    for (fi, face) in shell.faces.iter().enumerate() {
        // Collect outer-wire vertex positions (in order).
        let mut pts: Vec<DVec3> = Vec::new();
        for we in &face.outer_wire.edges {
            let Some(edge) = brep.edges.get(we.idx) else {
                continue;
            };
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                pts.push(v.point);
            }
        }

        if pts.len() < 3 {
            // Degenerate -> counts as small.
            result.push(fi);
            continue;
        }

        // Fan-triangulation area from pts[0].
        let mut area = 0.0f64;
        let p0 = pts[0];
        for i in 1..pts.len() - 1 {
            area += (pts[i] - p0).cross(pts[i + 1] - p0).length() * 0.5;
        }

        if area <= max_area {
            result.push(fi);
        }
    }

    result
}

// -- Fill helpers ------------------------------------------------------------

/// Build a fill cylinder BRep that covers a cylindrical hole, extended by
/// `margin` on each side.
fn make_fill_cylinder(
    feature: &CylindricalFeature,
    margin: f64,
) -> Result<BRep, rcad_modeling::BuildError> {
    let ax = feature.axis;
    let height = feature.height() + 2.0 * margin;
    // Base center of the fill cylinder (slightly below t_min).
    let base_pt = feature.origin + ax * (feature.t_min - margin);
    // A reference direction perpendicular to the axis (needed for seam placement).
    let ref_dir = any_perpendicular(ax);
    // Expand radius slightly (10x TOLERANCE_ABS) so the boolean unambiguously
    // fills the hole even at analytic floating-point surfaces.
    let expanded_r = feature.radius + TOLERANCE_ABS * 10.0;
    make_cylinder_brep(base_pt, ax, ref_dir, expanded_r, height)
}

/// Build a boss cylinder BRep to subtract from the host for boss removal.
fn make_boss_cylinder(
    feature: &CylindricalFeature,
    margin: f64,
) -> Result<BRep, rcad_modeling::BuildError> {
    // Same geometry as hole fill -> boolean Difference is used instead of Union.
    make_fill_cylinder(feature, margin)
}

// -- Main API ----------------------------------------------------------------

/// Perform a defeaturing pass on `brep`, suppressing small cylindrical holes
/// and bosses according to `options`.
///
/// Returns the modified BRep and a [`DefeaturingReport`] describing the
/// changes.  The input BRep is not modified.
///
/// # Errors
///
/// Returns [`DefeaturingError::EmptyInput`] if `brep` has no solids/shells.
///
/// # Notes
///
/// * Only `solids[0].shells[0]` is inspected for features.  Multi-solid BReps
///   are processed as a whole through boolean operations.
/// * A feature that causes a boolean failure is counted in
///   [`DefeaturingReport::failed_features`]; the pass continues with
///   remaining features.
/// * When `enable_retry` is enabled, failed boolean operations are retried
///   with increased fuzzy tolerance according to `retry_fuzzy_multiplier`.
/// * When `run_post_healing` is enabled, `make_connected_enhanced` is called
///   after all features are processed to repair connectivity.
pub fn defeature_brep(
    brep: &BRep,
    options: &DefeaturingOptions,
) -> Result<(BRep, DefeaturingReport), DefeaturingError> {
    if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
        return Err(DefeaturingError::EmptyInput);
    }

    let mut report = DefeaturingReport::default();
    let mut current = brep.clone();

    // -- Small-face identification ------------------------------------------
    if options.max_small_face_area > 0.0 {
        report.small_faces_identified =
            identify_small_faces(&current, options.max_small_face_area).len();
    }

    // -- Cylindrical holes and bosses ---------------------------------------
    let needs_cyl = options.max_hole_radius > 0.0 || options.max_boss_radius > 0.0;
    if needs_cyl {
        let features = detect_cylindrical_features(
            &current,
            options.max_hole_radius,
            options.max_boss_radius,
        );

        let margin = if options.fill_margin > 0.0 {
            options.fill_margin
        } else {
            DEFAULT_FILL_MARGIN
        };

        for feature in &features {
            // Guard each operation by the applicable threshold; a feature may
            // be in the detection pool (<= effective_max) yet outside the
            // specific threshold for its operation type.
            if feature.is_hole {
                if options.max_hole_radius <= 0.0 || feature.radius > options.max_hole_radius {
                    continue;
                }
                match make_fill_cylinder(feature, margin) {
                    Ok(fill) => {
                        let result = if options.enable_retry {
                            try_boolean_with_retry(
                                BooleanOpType::Union,
                                &current,
                                &fill,
                                options.retry_fuzzy_multiplier,
                                options.max_retries,
                                &mut report,
                            )
                        } else {
                            boolean_op(BooleanOpType::Union, &current, &fill)
                                .map(|b| (b, false))
                        };
                        match result {
                            Ok((new_brep, retried)) => {
                                current = new_brep;
                                report.holes_removed += 1;
                                if retried {
                                    report.succeeded_after_retry += 1;
                                }
                            }
                            Err(_) => {
                                report.failed_features += 1;
                            }
                        }
                    }
                    Err(_) => {
                        report.failed_features += 1;
                    }
                }
            } else {
                if options.max_boss_radius <= 0.0 || feature.radius > options.max_boss_radius {
                    continue;
                }
                match make_boss_cylinder(feature, margin) {
                    Ok(boss) => {
                        let result = if options.enable_retry {
                            try_boolean_with_retry(
                                BooleanOpType::Difference,
                                &current,
                                &boss,
                                options.retry_fuzzy_multiplier,
                                options.max_retries,
                                &mut report,
                            )
                        } else {
                            boolean_op(BooleanOpType::Difference, &current, &boss)
                                .map(|b| (b, false))
                        };
                        match result {
                            Ok((new_brep, retried)) => {
                                current = new_brep;
                                report.bosses_removed += 1;
                                if retried {
                                    report.succeeded_after_retry += 1;
                                }
                            }
                            Err(_) => {
                                report.failed_features += 1;
                            }
                        }
                    }
                    Err(_) => {
                        report.failed_features += 1;
                    }
                }
            }
        }
    }

    // -- Conical features ---------------------------------------------------
    if options.enable_conical_features && options.max_conical_hole_radius > 0.0 {
        let features = detect_conical_features(&current, options.max_conical_hole_radius);

        for feature in &features {
            if !feature.is_hole {
                // Boss removal for cones not yet implemented.
                continue;
            }

            // Build a fill cone using a cylinder approximation for now.
            // A proper implementation would construct a conical solid.
            match make_fill_cone(feature, options.fill_margin) {
                Ok(fill) => {
                    let result = if options.enable_retry {
                        try_boolean_with_retry(
                            BooleanOpType::Union,
                            &current,
                            &fill,
                            options.retry_fuzzy_multiplier,
                            options.max_retries,
                            &mut report,
                        )
                    } else {
                        boolean_op(BooleanOpType::Union, &current, &fill)
                            .map(|b| (b, false))
                    };
                    match result {
                        Ok((new_brep, retried)) => {
                            current = new_brep;
                            report.conical_features_removed += 1;
                            if retried {
                                report.succeeded_after_retry += 1;
                            }
                        }
                        Err(_) => {
                            report.failed_features += 1;
                        }
                    }
                }
                Err(_) => {
                    report.failed_features += 1;
                }
            }
        }
    }

    // -- Post-defeature healing ---------------------------------------------
    if options.run_post_healing {
        let (healed_brep, heal_report) =
            make_connected_enhanced(&current, options.healing_tolerance, 3);
        current = healed_brep;
        report.healing_performed = true;
        report.healing_vertices_merged = heal_report.vertices_merged;
        report.healing_small_edges_removed = heal_report.small_edges_removed;
    }

    Ok((current, report))
}

/// Try a boolean operation with retry using increased fuzzy tolerance.
///
/// Returns `Ok((brep, true))` if succeeded after retry, `Ok((brep, false))` if
/// succeeded on first try, or `Err` if all attempts failed.
fn try_boolean_with_retry(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
    fuzzy_multiplier: f64,
    max_retries: usize,
    report: &mut DefeaturingReport,
) -> Result<(BRep, bool), crate::BooleanError> {
    // First attempt with default fuzzy tolerance.
    match boolean_op(op, a, b) {
        Ok(result) => Ok((result, false)),
        Err(first_err) => {
            // Build retry ladder based on multiplier.
            let ladder: Vec<f64> = (1..=max_retries)
                .map(|i| TOLERANCE_ABS * fuzzy_multiplier * (i as f64))
                .collect();

            let robust_opts = BooleanRobustOptions {
                base: BooleanOptions::default(),
                fuzzy_retry_ladder: ladder,
                retry_policy: BooleanRetryPolicy::Aggressive,
                extreme_geometry: crate::ExtremeGeometryRetryConfig::default(),
            };

            report.retry_attempts += 1;

            match boolean_op_robust(op, a, b, robust_opts) {
                Ok((result, _exec_report)) => {
                    report.retry_attempts += _exec_report.retry_count;
                    Ok((result, true))
                }
                Err(_) => {
                    Err(first_err)
                }
            }
        }
    }
}

/// Build a fill solid for a conical feature.
///
/// Currently uses a cylinder approximation. A proper implementation would
/// construct a conical solid matching the feature geometry.
fn make_fill_cone(
    feature: &ConicalFeature,
    margin: f64,
) -> Result<BRep, rcad_modeling::BuildError> {
    // Use the reference radius and height to build an approximating cylinder.
    let ax = feature.axis;
    let height = (feature.t_max - feature.t_min).max(0.0) + 2.0 * margin;
    let base_pt = feature.apex + ax * (feature.t_min - margin);
    let ref_dir = any_perpendicular(ax);
    // Expand radius slightly to ensure coverage.
    let expanded_r = feature.reference_radius + TOLERANCE_ABS * 10.0;
    make_cylinder_brep(base_pt, ax, ref_dir, expanded_r, height)
}

// -- Enhanced Defeaturing with Feature Groups and Robust Healing ------------

/// Enhanced defeaturing options with deeper control.
#[derive(Debug, Clone)]
pub struct DefeaturingOptionsEnhanced {
    /// Base defeaturing options.
    pub base: DefeaturingOptions,
    /// Process features in connected groups.
    pub process_feature_groups: bool,
    /// Run enhanced post-healing with strategy.
    pub enhanced_healing: bool,
    /// MakeConnectedStrategy for post-healing.
    pub healing_strategy: crate::brep_repair::MakeConnectedStrategy,
    /// Maximum number of features to process in a single boolean operation.
    pub max_features_per_operation: usize,
    /// Enable adaptive tolerance for difficult features.
    pub adaptive_tolerance: bool,
    /// Tolerance growth factor for adaptive mode.
    pub tolerance_growth_factor: f64,
}

impl Default for DefeaturingOptionsEnhanced {
    fn default() -> Self {
        Self {
            base: DefeaturingOptions::default(),
            process_feature_groups: true,
            enhanced_healing: true,
            healing_strategy: crate::brep_repair::MakeConnectedStrategy::standard(),
            max_features_per_operation: 5,
            adaptive_tolerance: true,
            tolerance_growth_factor: 2.0,
        }
    }
}

impl DefeaturingOptionsEnhanced {
    /// Create conservative options (slower but safer).
    pub fn conservative() -> Self {
        Self {
            base: DefeaturingOptions {
                enable_retry: true,
                max_retries: 5,
                retry_fuzzy_multiplier: 5.0,
                ..Default::default()
            },
            process_feature_groups: false,
            enhanced_healing: true,
            healing_strategy: crate::brep_repair::MakeConnectedStrategy::conservative(),
            max_features_per_operation: 1,
            adaptive_tolerance: true,
            tolerance_growth_factor: 1.5,
        }
    }

    /// Create aggressive options (faster but may miss edge cases).
    pub fn aggressive() -> Self {
        Self {
            base: DefeaturingOptions {
                enable_retry: true,
                max_retries: 3,
                retry_fuzzy_multiplier: 20.0,
                run_post_healing: false, // We'll use enhanced healing
                ..Default::default()
            },
            process_feature_groups: true,
            enhanced_healing: true,
            healing_strategy: crate::brep_repair::MakeConnectedStrategy::aggressive(),
            max_features_per_operation: 10,
            adaptive_tolerance: true,
            tolerance_growth_factor: 3.0,
        }
    }

    /// Create options optimized for injection molding.
    pub fn for_injection_molding() -> Self {
        Self {
            base: DefeaturingOptions {
                max_hole_radius: 3.0,
                max_boss_radius: 2.0,
                enable_conical_features: true,
                max_conical_hole_radius: 3.0,
                enable_slot_features: true,
                max_slot_width: 5.0,
                max_slot_depth: 10.0,
                enable_blend_features: true,
                max_blend_radius: 2.0,
                enable_retry: true,
                max_retries: 3,
                run_post_healing: false,
                ..Default::default()
            },
            process_feature_groups: true,
            enhanced_healing: true,
            healing_strategy: crate::brep_repair::MakeConnectedStrategy::for_injection_molding(),
            max_features_per_operation: 5,
            adaptive_tolerance: true,
            tolerance_growth_factor: 2.0,
        }
    }
}

/// Enhanced report with additional details.
#[derive(Debug, Clone, Default)]
pub struct DefeaturingReportEnhanced {
    /// Base report.
    pub base: DefeaturingReport,
    /// Number of feature groups processed.
    pub groups_processed: usize,
    /// Number of features processed in groups.
    pub features_in_groups: usize,
    /// Post-healing report.
    pub healing_report: Option<crate::brep_repair::MakeConnectedReport>,
    /// Adaptive tolerance escalations.
    pub tolerance_escalations: usize,
    /// Features that required multiple attempts.
    pub multi_attempt_features: usize,
}

/// Enhanced defeaturing with feature group processing and robust healing.
///
/// This function extends `defeature_brep` with:
/// - Feature group detection and batch processing
/// - Integration with `MakeConnectedStrategy` for post-healing
/// - Adaptive tolerance escalation for difficult features
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `options` - Enhanced options
///
/// # Returns
/// Defeatured B-Rep and detailed report.
pub fn defeature_brep_enhanced(
    brep: &BRep,
    options: &DefeaturingOptionsEnhanced,
) -> Result<(BRep, DefeaturingReportEnhanced), DefeaturingError> {
    if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
        return Err(DefeaturingError::EmptyInput);
    }

    let mut report = DefeaturingReportEnhanced::default();
    let mut current = brep.clone();

    // Detect all feature types
    let cylindrical_features = if options.base.max_hole_radius > 0.0 || options.base.max_boss_radius > 0.0 {
        detect_cylindrical_features(&current, options.base.max_hole_radius, options.base.max_boss_radius)
    } else {
        Vec::new()
    };

    let conical_features = if options.base.enable_conical_features && options.base.max_conical_hole_radius > 0.0 {
        detect_conical_features(&current, options.base.max_conical_hole_radius)
    } else {
        Vec::new()
    };

    let slot_features = if options.base.enable_slot_features && options.base.max_slot_width > 0.0 {
        detect_slot_features(&current, options.base.max_slot_width, options.base.max_slot_depth)
    } else {
        Vec::new()
    };

    let pocket_features = if options.base.enable_pocket_features && options.base.max_pocket_diameter > 0.0 {
        detect_pocket_features(&current, options.base.max_pocket_diameter, options.base.max_pocket_depth)
    } else {
        Vec::new()
    };

    let blend_features = if options.base.enable_blend_features {
        detect_blend_features(&current, options.base.max_blend_radius, options.base.max_chamfer_distance)
    } else {
        Vec::new()
    };

    // Small face identification
    if options.base.max_small_face_area > 0.0 {
        report.base.small_faces_identified =
            identify_small_faces(&current, options.base.max_small_face_area).len();
    }

    // Process feature groups if enabled
    if options.process_feature_groups {
        let (groups, _face_to_group) = detect_connected_feature_groups(
            &current,
            &cylindrical_features,
            &conical_features,
            &slot_features,
            &pocket_features,
            &blend_features,
        );

        report.groups_processed = groups.len();

        // Process each group
        for group in &groups {
            let group_result = process_feature_group(
                &current,
                group,
                &cylindrical_features,
                &conical_features,
                options,
                &mut report,
            );

            if let Ok(new_brep) = group_result {
                current = new_brep;
                report.features_in_groups += group.total_faces;
            }
        }
    } else {
        // Process features individually (use base function)
        let (new_brep, base_report) = defeature_brep(&current, &options.base)?;
        current = new_brep;
        report.base = base_report;
    }

    // Enhanced post-healing with strategy
    if options.enhanced_healing {
        let (healed, healing_report) = options.healing_strategy.apply(&current);
        current = healed;
        report.healing_report = Some(healing_report);
    }

    Ok((current, report))
}

/// Process a feature group as a batch.
fn process_feature_group(
    brep: &BRep,
    group: &FeatureGroup,
    cylindrical_features: &[CylindricalFeature],
    conical_features: &[ConicalFeature],
    options: &DefeaturingOptionsEnhanced,
    report: &mut DefeaturingReportEnhanced,
) -> Result<BRep, DefeaturingError> {
    let mut current = brep.clone();
    let margin = if options.base.fill_margin > 0.0 {
        options.base.fill_margin
    } else {
        DEFAULT_FILL_MARGIN
    };

    // Process cylindrical features in this group
    for &idx in &group.cylindrical_indices {
        if let Some(feature) = cylindrical_features.get(idx) {
            let fill_result = if feature.is_hole {
                make_fill_cylinder(feature, margin)
            } else {
                make_boss_cylinder(feature, margin)
            };

            if let Ok(fill) = fill_result {
                let op = if feature.is_hole {
                    BooleanOpType::Union
                } else {
                    BooleanOpType::Difference
                };

                let result = if options.base.enable_retry {
                    let ladder: Vec<f64> = (1..=options.base.max_retries)
                        .map(|i| TOLERANCE_ABS * options.base.retry_fuzzy_multiplier * (i as f64))
                        .collect();

                    let robust_opts = BooleanRobustOptions {
                        base: BooleanOptions::default(),
                        fuzzy_retry_ladder: ladder,
                        retry_policy: BooleanRetryPolicy::AdaptiveByFailureClass,
                        extreme_geometry: crate::ExtremeGeometryRetryConfig::default(),
                    };

                    report.base.retry_attempts += 1;
                    boolean_op_robust(op, &current, &fill, robust_opts)
                        .map(|(b, _)| b)
                } else {
                    boolean_op(op, &current, &fill)
                };

                match result {
                    Ok(new_brep) => {
                        current = new_brep;
                        if feature.is_hole {
                            report.base.holes_removed += 1;
                        } else {
                            report.base.bosses_removed += 1;
                        }
                    }
                    Err(_) => {
                        report.base.failed_features += 1;

                        // Try adaptive tolerance escalation
                        if options.adaptive_tolerance {
                            let tol = TOLERANCE_ABS * options.tolerance_growth_factor;
                            let (retried, mc_report) = crate::brep_repair::MakeConnectedStrategy {
                                merge_tolerance: tol,
                                ..crate::brep_repair::MakeConnectedStrategy::default()
                            }.apply(&current);

                            current = retried;
                            report.tolerance_escalations += 1;
                            let _ = mc_report;
                        }
                    }
                }
            }
        }
    }

    // Process conical features
    for &idx in &group.conical_indices {
        if let Some(feature) = conical_features.get(idx)
            && feature.is_hole
                && let Ok(fill) = make_fill_cone(feature, margin)
                    && boolean_op(BooleanOpType::Union, &current, &fill).is_ok() {
                        report.base.conical_features_removed += 1;
                    }
    }

    Ok(current)
}

// -- Tests -------------------------------------------------------------------
