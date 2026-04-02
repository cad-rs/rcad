use super::{
    basis_from_axis_ref, basis_from_x_y, transform_brep, translate_brep, validate_point,
    validate_positive, validate_segments, BuildError,
};
use glam::DVec3;
use rcad_kernel::{BRep, PrimitiveSolid};

pub fn box_primitive(width: f64, height: f64, depth: f64) -> Result<PrimitiveSolid, BuildError> {
    let width = validate_positive("width", width)?;
    let height = validate_positive("height", height)?;
    let depth = validate_positive("depth", depth)?;
    Ok(PrimitiveSolid::Box {
        width,
        height,
        depth,
    })
}

pub fn make_box_primitive(
    width: f64,
    height: f64,
    depth: f64,
) -> Result<PrimitiveSolid, BuildError> {
    box_primitive(width, height, depth)
}

pub fn box_brep(
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    width: f64,
    height: f64,
    depth: f64,
) -> Result<BRep, BuildError> {
    let origin = validate_point("origin", origin)?;
    let primitive = box_primitive(width, height, depth)?;
    let (x_axis, y_axis, z_axis) = basis_from_x_y(x_dir, y_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, origin, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_box_brep(
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    width: f64,
    height: f64,
    depth: f64,
) -> Result<BRep, BuildError> {
    box_brep(origin, x_dir, y_dir, width, height, depth)
}

pub fn sphere_primitive(
    radius: f64,
    u_segments: usize,
    v_segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    validate_segments("u_segments", u_segments, 3)?;
    validate_segments("v_segments", v_segments, 2)?;
    Ok(PrimitiveSolid::Sphere {
        radius,
        u_segments,
        v_segments,
    })
}

pub fn make_sphere_primitive(
    radius: f64,
    u_segments: usize,
    v_segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    sphere_primitive(radius, u_segments, v_segments)
}

pub fn sphere_brep(
    center: DVec3,
    radius: f64,
    u_segments: usize,
    v_segments: usize,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = sphere_primitive(radius, u_segments, v_segments)?;
    let mut brep = BRep::from_primitive(primitive);
    translate_brep(&mut brep, center);
    Ok(brep)
}

pub fn make_sphere_brep(
    center: DVec3,
    radius: f64,
    u_segments: usize,
    v_segments: usize,
) -> Result<BRep, BuildError> {
    sphere_brep(center, radius, u_segments, v_segments)
}

pub fn cylinder_primitive(
    radius: f64,
    height: f64,
    segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    let height = validate_positive("height", height)?;
    validate_segments("segments", segments, 3)?;
    Ok(PrimitiveSolid::Cylinder {
        radius,
        height,
        segments,
    })
}

pub fn make_cylinder_primitive(
    radius: f64,
    height: f64,
    segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    cylinder_primitive(radius, height, segments)
}

pub fn cylinder_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
    segments: usize,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = cylinder_primitive(radius, height, segments)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_cylinder_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
    segments: usize,
) -> Result<BRep, BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height, segments)
}

pub fn cone_primitive(
    base_radius: f64,
    height: f64,
    segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    let base_radius = validate_positive("base_radius", base_radius)?;
    let height = validate_positive("height", height)?;
    validate_segments("segments", segments, 3)?;
    Ok(PrimitiveSolid::Cone {
        base_radius,
        height,
        segments,
    })
}

pub fn make_cone_primitive(
    base_radius: f64,
    height: f64,
    segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    cone_primitive(base_radius, height, segments)
}

pub fn cone_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    base_radius: f64,
    height: f64,
    segments: usize,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = cone_primitive(base_radius, height, segments)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_cone_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    base_radius: f64,
    height: f64,
    segments: usize,
) -> Result<BRep, BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, height, segments)
}

pub fn torus_primitive(
    major_radius: f64,
    minor_radius: f64,
    major_segments: usize,
    minor_segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    let major_radius = validate_positive("major_radius", major_radius)?;
    let minor_radius = validate_positive("minor_radius", minor_radius)?;
    validate_segments("major_segments", major_segments, 3)?;
    validate_segments("minor_segments", minor_segments, 3)?;
    Ok(PrimitiveSolid::Torus {
        major_radius,
        minor_radius,
        major_segments,
        minor_segments,
    })
}

pub fn make_torus_primitive(
    major_radius: f64,
    minor_radius: f64,
    major_segments: usize,
    minor_segments: usize,
) -> Result<PrimitiveSolid, BuildError> {
    torus_primitive(major_radius, minor_radius, major_segments, minor_segments)
}

pub fn torus_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
    major_segments: usize,
    minor_segments: usize,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = torus_primitive(major_radius, minor_radius, major_segments, minor_segments)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_torus_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
    major_segments: usize,
    minor_segments: usize,
) -> Result<BRep, BuildError> {
    torus_brep(
        center,
        axis,
        ref_dir,
        major_radius,
        minor_radius,
        major_segments,
        minor_segments,
    )
}