/// Counts of key STEP entities for validation.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct StepEntityCounts {
    pub manifold_solid_brep: usize,
    pub plane: usize,
    pub cylindrical_surface: usize,
    pub conical_surface: usize,
    pub spherical_surface: usize,
    pub toroidal_surface: usize,
    pub surface_of_revolution: usize,
    pub surface_of_linear_extrusion: usize,
    pub offset_surface: usize,
    pub b_spline_surface_with_knots: usize,
}

/// Count STEP entity types in a STEP string.
pub fn count_step_entities_from_str(content: &str) -> StepEntityCounts {
    let mut counts = StepEntityCounts::default();

    for line in content.lines() {
        // Skip header and empty lines
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with("ISO-10303-21")
            || trimmed.starts_with("HEADER")
            || trimmed.starts_with("ENDSEC")
            || trimmed.starts_with("DATA")
            || trimmed.starts_with("END-ISO")
        {
            continue;
        }
        // Skip entity ID prefix (#N=)
        let after_id = if let Some(eq_pos) = trimmed.find('=') {
            &trimmed[eq_pos + 1..]
        } else {
            continue;
        };

        if after_id.starts_with("MANIFOLD_SOLID_BREP") {
            counts.manifold_solid_brep += 1;
        } else if after_id.starts_with("PLANE(") || after_id.starts_with("PLANE (") {
            counts.plane += 1;
        } else if after_id.starts_with("CYLINDRICAL_SURFACE") {
            counts.cylindrical_surface += 1;
        } else if after_id.starts_with("CONICAL_SURFACE") {
            counts.conical_surface += 1;
        } else if after_id.starts_with("SPHERICAL_SURFACE") {
            counts.spherical_surface += 1;
        } else if after_id.starts_with("TOROIDAL_SURFACE") {
            counts.toroidal_surface += 1;
        } else if after_id.starts_with("SURFACE_OF_REVOLUTION") {
            counts.surface_of_revolution += 1;
        } else if after_id.starts_with("SURFACE_OF_LINEAR_EXTRUSION") {
            counts.surface_of_linear_extrusion += 1;
        } else if after_id.starts_with("OFFSET_SURFACE") {
            counts.offset_surface += 1;
        } else if after_id.starts_with("B_SPLINE_SURFACE_WITH_KNOTS") {
            counts.b_spline_surface_with_knots += 1;
        }
    }

    counts
}
