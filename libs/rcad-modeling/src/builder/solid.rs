use super::{
    BuildError, basis_from_axis_ref, basis_from_x_y, do_mirror_brep, normalize_vector,
    transform_brep, translate_brep, validate_point, validate_positive,
};
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
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

pub fn sphere_primitive(radius: f64) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    Ok(PrimitiveSolid::Sphere { radius })
}

pub fn make_sphere_primitive(radius: f64) -> Result<PrimitiveSolid, BuildError> {
    sphere_primitive(radius)
}

pub fn sphere_brep(center: DVec3, radius: f64) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = sphere_primitive(radius)?;
    let mut brep = BRep::from_primitive(primitive);
    translate_brep(&mut brep, center);
    Ok(brep)
}

pub fn make_sphere_brep(center: DVec3, radius: f64) -> Result<BRep, BuildError> {
    sphere_brep(center, radius)
}

pub fn cylinder_primitive(radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    let height = validate_positive("height", height)?;
    Ok(PrimitiveSolid::Cylinder { radius, height })
}

pub fn make_cylinder_primitive(radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    cylinder_primitive(radius, height)
}

pub fn cylinder_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = cylinder_primitive(radius, height)?;
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
) -> Result<BRep, BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height)
}

pub fn cone_primitive(base_radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    let base_radius = validate_positive("base_radius", base_radius)?;
    let height = validate_positive("height", height)?;
    Ok(PrimitiveSolid::Cone {
        base_radius,
        height,
    })
}

pub fn make_cone_primitive(base_radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    cone_primitive(base_radius, height)
}

pub fn cone_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    base_radius: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = cone_primitive(base_radius, height)?;
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
) -> Result<BRep, BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, height)
}

/// Right conical frustum (truncated cone), matching OCCT `pcone` when both bottom and top radii are positive.
///
/// `center` is the **midpoint** between the circular face centers; `axis` points from the bottom face
/// toward the top face; `r_bottom` / `r_top` are radii in those end planes; `height` is the distance
/// between the planes. Built with [`super::ops::loft`] between two regular polygon approximations.
pub fn make_conical_frustum_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    r_bottom: f64,
    r_top: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let height = validate_positive("height", height)?;
    let rb = validate_positive("r_bottom", r_bottom)?;
    let rt = validate_positive("r_top", r_top)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let half_h = height * 0.5;
    let bottom_c = center - y_axis * half_h;
    let top_c = center + y_axis * half_h;
    const N: usize = 32;
    let mut lo = Vec::with_capacity(N);
    let mut hi = Vec::with_capacity(N);
    use std::f64::consts::TAU;
    for i in 0..N {
        let ang = TAU * (i as f64) / (N as f64);
        let c = ang.cos();
        let s = ang.sin();
        lo.push(bottom_c + x_axis * (c * rb) + z_axis * (s * rb));
        hi.push(top_c + x_axis * (c * rt) + z_axis * (s * rt));
    }
    super::ops::loft(&[lo, hi])
}

pub fn torus_primitive(major_radius: f64, minor_radius: f64) -> Result<PrimitiveSolid, BuildError> {
    let major_radius = validate_positive("major_radius", major_radius)?;
    let minor_radius = validate_positive("minor_radius", minor_radius)?;
    Ok(PrimitiveSolid::Torus {
        major_radius,
        minor_radius,
    })
}

pub fn make_torus_primitive(
    major_radius: f64,
    minor_radius: f64,
) -> Result<PrimitiveSolid, BuildError> {
    torus_primitive(major_radius, minor_radius)
}

pub fn torus_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = torus_primitive(major_radius, minor_radius)?;
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
) -> Result<BRep, BuildError> {
    torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}

/// Create a BRep with a single rectangular planar face.
///
/// `origin` is the plane origin; `u_axis`/`v_axis` are the parametric axes
/// (the normal is `u_axis × v_axis`). The corners are `origin + u_axis*u + v_axis*v`.
pub fn make_planar_rect_brep(
    origin: DVec3,
    u_axis: DVec3,
    v_axis: DVec3,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> Result<BRep, BuildError> {
    let normal = u_axis.cross(v_axis).normalize();
    let surface = Surface3::Plane(Plane { origin, normal });

    let c0 = origin + u_axis * umin + v_axis * vmin;
    let c1 = origin + u_axis * umax + v_axis * vmin;
    let c2 = origin + u_axis * umax + v_axis * vmax;
    let c3 = origin + u_axis * umin + v_axis * vmax;

    let mut brep = BRep::default();
    let v0 = super::brep_builder::make_vertex(&mut brep, c0);
    let v1 = super::brep_builder::make_vertex(&mut brep, c1);
    let v2 = super::brep_builder::make_vertex(&mut brep, c2);
    let v3 = super::brep_builder::make_vertex(&mut brep, c3);

    let len01 = (c1 - c0).length();
    let len12 = (c2 - c1).length();
    let len23 = (c3 - c2).length();
    let len30 = (c0 - c3).length();

    let dir01 = (c1 - c0).normalize();
    let dir12 = (c2 - c1).normalize();
    let dir23 = (c3 - c2).normalize();
    let dir30 = (c0 - c3).normalize();

    let e0 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c0, direction: dir01 }),
        0.0, len01, v0, v1,
    )?;
    let e1 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c1, direction: dir12 }),
        0.0, len12, v1, v2,
    )?;
    let e2 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c2, direction: dir23 }),
        0.0, len23, v2, v3,
    )?;
    let e3 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c3, direction: dir30 }),
        0.0, len30, v3, v0,
    )?;

    let w = super::brep_builder::make_wire(vec![
        WireEdge::new(e0, true),
        WireEdge::new(e1, true),
        WireEdge::new(e2, true),
        WireEdge::new(e3, true),
    ]);

    super::brep_builder::make_face(&mut brep, surface, w, vec![])?;
    Ok(brep)
}

/// Create a BRep with a single planar face bounded by a polygon, with optional inner holes.
///
/// `origin` and `normal` define the plane. `outer_polygon` gives the outer boundary
/// vertices (≥3 points in CCW order). `inner_polygons` are hole boundaries (each ≥3 points
/// in CW order).
pub fn make_planar_polygon_brep(
    origin: DVec3,
    normal: DVec3,
    outer_polygon: &[DVec3],
    inner_polygons: &[Vec<DVec3>],
) -> Result<BRep, BuildError> {
    let surface = Surface3::Plane(Plane { origin, normal });
    let mut brep = BRep::default();
    let outer_wire = make_polygon_wire(&mut brep, outer_polygon)?;
    let mut inner_wires = Vec::new();
    for poly in inner_polygons {
        inner_wires.push(make_polygon_wire(&mut brep, poly)?);
    }
    super::brep_builder::make_face(&mut brep, surface, outer_wire, inner_wires)?;
    Ok(brep)
}

fn make_polygon_wire(brep: &mut BRep, points: &[DVec3]) -> Result<rcad_kernel::topology::Wire, BuildError> {
    let n = points.len();
    if n < 3 {
        return Err(BuildError::DegenerateGeometry("polygon needs at least 3 points"));
    }
    let mut vertex_indices = Vec::with_capacity(n);
    for p in points {
        vertex_indices.push(super::brep_builder::make_vertex(brep, *p));
    }
    let mut wire_edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let vi = vertex_indices[i];
        let vj = vertex_indices[j];
        let dir = (points[j] - points[i]).normalize();
        let len = (points[j] - points[i]).length();
        let ei = super::brep_builder::make_edge(
            brep,
            Curve3::Line(Line3 { origin: points[i], direction: dir }),
            0.0, len, vi, vj,
        )?;
        wire_edges.push(WireEdge::new(ei, true));
    }
    Ok(super::brep_builder::make_wire(wire_edges))
}

/// Mirror a BRep across a plane defined by `origin` and `normal`.
///
/// The mirrored BRep has inverted face normals and reversed wire orientations
/// to maintain consistent outward-facing normals.
pub fn mirror_brep(brep: &BRep, plane_origin: DVec3, plane_normal: DVec3) -> Result<BRep, BuildError> {
    let _ = validate_point("plane_origin", plane_origin)?;
    let n = normalize_vector("plane_normal", plane_normal)?;
    Ok(do_mirror_brep(brep, plane_origin, n))
}
