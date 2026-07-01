// =============================================================================
// Shape Modification Utilities
// =============================================================================

/// Apply an affine transformation to a shape in-place.
///
/// This transforms all vertices, curves, surfaces, and face normals.
///
/// # Example
///
/// ```
/// # use rcad_algorithms::tolerance::*;
/// use rcad_algorithms::brep_tools::transform_shape;
/// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
/// use rcad_kernel::BRep;
/// use glam::{DAffine3, DVec3};
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let translation = DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0));
/// transform_shape(&mut brep, translation);
/// // The box is now centered at (5.5, 0.5, 0.5)
/// assert!((brep.vertices[0].point.x - 5.0).abs() < TOLERANCE_COORD_SUB);
/// ```
pub fn transform_shape(brep: &mut BRep, transform: DAffine3) {
    brep.apply_transform(transform);
}

/// Mirror a shape across a plane.
///
/// # Arguments
///
/// * `brep` - The BRep to mirror (modified in place)
/// * `plane_origin` - A point on the mirror plane
/// * `plane_normal` - Normal vector of the mirror plane (will be normalized)
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::mirror_shape;
/// use rcad_kernel::BRep;
/// use glam::DVec3;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Mirror across the YZ plane (x = 0)
/// mirror_shape(&mut brep, DVec3::ZERO, DVec3::X);
/// // The box is now in the negative X half-space
/// ```
pub fn mirror_shape(brep: &mut BRep, plane_origin: DVec3, plane_normal: DVec3) {
    let normal = plane_normal.normalize_or(DVec3::X);

    // Reflection matrix: R = I - 2 * n * n^T
    // Where n is the normalized plane normal
    let mat = DMat4::from_cols(
        DVec4::new(1.0 - 2.0 * normal.x * normal.x, -2.0 * normal.x * normal.y, -2.0 * normal.x * normal.z, 0.0),
        DVec4::new(-2.0 * normal.y * normal.x, 1.0 - 2.0 * normal.y * normal.y, -2.0 * normal.y * normal.z, 0.0),
        DVec4::new(-2.0 * normal.z * normal.x, -2.0 * normal.z * normal.y, 1.0 - 2.0 * normal.z * normal.z, 0.0),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    );

    // Combine: translate to origin, reflect, translate back
    let to_origin = DMat4::from_translation(-plane_origin);
    let from_origin = DMat4::from_translation(plane_origin);
    let transform_mat = from_origin * mat * to_origin;

    // Convert DMat4 to DAffine3
    let transform = DAffine3::from_cols(
        transform_mat.x_axis.truncate(),
        transform_mat.y_axis.truncate(),
        transform_mat.z_axis.truncate(),
        transform_mat.w_axis.truncate(),
    );

    brep.apply_transform(transform);
}

/// Scale a shape about a center point.
///
/// # Arguments
///
/// * `brep` - The BRep to scale (modified in place)
/// * `factor` - Uniform scale factor
/// * `center` - Center point for scaling
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::scale_shape;
/// use rcad_kernel::BRep;
/// use glam::DVec3;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Scale by 2x about the origin
/// scale_shape(&mut brep, 2.0, DVec3::ZERO);
/// // The box is now 2x2x2
/// ```
pub fn scale_shape(brep: &mut BRep, factor: f64, center: DVec3) {
    let _transform = DAffine3::from_scale(glam::DVec3::splat(factor))
        * DAffine3::from_translation(-center)
        * DAffine3::from_translation(center);

    // Actually we need: translate to origin, scale, translate back
    let to_origin = DAffine3::from_translation(-center);
    let scale = DAffine3::from_scale(glam::DVec3::splat(factor));
    let from_origin = DAffine3::from_translation(center);
    let final_transform = from_origin * scale * to_origin;

    brep.apply_transform(final_transform);
}

/// Rotate a shape about an axis.
///
/// # Arguments
///
/// * `brep` - The BRep to rotate (modified in place)
/// * `axis_origin` - A point on the rotation axis
/// * `axis_direction` - Direction of the rotation axis
/// * `angle` - Rotation angle in radians
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::rotate_shape;
/// use rcad_kernel::BRep;
/// use glam::DVec3;
/// use std::f64::consts::PI;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Rotate 90 degrees about the Z axis
/// rotate_shape(&mut brep, DVec3::ZERO, DVec3::Z, PI / 2.0);
/// ```
pub fn rotate_shape(brep: &mut BRep, axis_origin: DVec3, axis_direction: DVec3, angle: f64) {
    let axis = axis_direction.normalize_or(DVec3::Z);

    // Rotation about an arbitrary axis through a point:
    // Translate to origin, rotate, translate back
    let to_origin = DAffine3::from_translation(-axis_origin);
    let rotation = DAffine3::from_axis_angle(axis, angle);
    let from_origin = DAffine3::from_translation(axis_origin);
    let transform = from_origin * rotation * to_origin;

    brep.apply_transform(transform);
}
