// =============================================================================
// Topods-native Shape Modification Utilities (migration)
// =============================================================================

/// Apply an affine transformation to a topods::BRep in-place.
pub fn transform_shape_topods(brep: &mut topods::BRep, transform: DAffine3) {
    brep.apply_transform(transform);
}

/// Mirror a topods::BRep across a plane.
pub fn mirror_shape_topods(brep: &mut topods::BRep, plane_origin: DVec3, plane_normal: DVec3) {
    let normal = plane_normal.normalize_or(DVec3::X);
    let mat = DMat4::from_cols(
        DVec4::new(1.0 - 2.0 * normal.x * normal.x, -2.0 * normal.x * normal.y, -2.0 * normal.x * normal.z, 0.0),
        DVec4::new(-2.0 * normal.y * normal.x, 1.0 - 2.0 * normal.y * normal.y, -2.0 * normal.y * normal.z, 0.0),
        DVec4::new(-2.0 * normal.z * normal.x, -2.0 * normal.z * normal.y, 1.0 - 2.0 * normal.z * normal.z, 0.0),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    );
    let to_origin = DMat4::from_translation(-plane_origin);
    let from_origin = DMat4::from_translation(plane_origin);
    let transform_mat = from_origin * mat * to_origin;
    let transform = DAffine3::from_cols(
        transform_mat.x_axis.truncate(),
        transform_mat.y_axis.truncate(),
        transform_mat.z_axis.truncate(),
        transform_mat.w_axis.truncate(),
    );
    brep.apply_transform(transform);
}

/// Scale a topods::BRep about a center point.
pub fn scale_shape_topods(brep: &mut topods::BRep, factor: f64, center: DVec3) {
    let to_origin = DAffine3::from_translation(-center);
    let scale = DAffine3::from_scale(glam::DVec3::splat(factor));
    let from_origin = DAffine3::from_translation(center);
    brep.apply_transform(from_origin * scale * to_origin);
}

/// Rotate a topods::BRep about an axis.
pub fn rotate_shape_topods(brep: &mut topods::BRep, axis_origin: DVec3, axis_direction: DVec3, angle: f64) {
    let axis = axis_direction.normalize_or(DVec3::Z);
    let to_origin = DAffine3::from_translation(-axis_origin);
    let rotation = DAffine3::from_axis_angle(axis, angle);
    let from_origin = DAffine3::from_translation(axis_origin);
    brep.apply_transform(from_origin * rotation * to_origin);
}
