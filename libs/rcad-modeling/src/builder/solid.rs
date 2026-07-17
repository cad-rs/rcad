use super::{
    BuildError, basis_from_axis_ref, basis_from_x_y, do_mirror_brep, normalize_vector,
    transform_brep, translate_brep, validate_point, validate_positive,
};
use glam::{DMat3, DVec3};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::{topods, BRep, PrimitiveSolid};

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
) -> Result<topods::BRep, BuildError> {
    make_box_brep(origin, x_dir, y_dir, width, height, depth)
}

pub fn make_box_brep(
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    width: f64,
    height: f64,
    depth: f64,
) -> Result<topods::BRep, BuildError> {
    let w = validate_positive("width", width)?;
    let h = validate_positive("height", height)?;
    let d = validate_positive("depth", depth)?;
    let (x_axis, y_axis, z_axis) = basis_from_x_y(x_dir, y_dir)?;
    let o = validate_point("origin", origin)?;
    use rcad_kernel::topods::Orientation;
    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };

    let p = |dx: f64, dy: f64, dz: f64| o + x_axis * dx + y_axis * dy + z_axis * dz;

    let mut t = topods::BRep::new();
    let v = [
        t.add_tvertex(p(0.0,0.0,0.0)), t.add_tvertex(p(w,0.0,0.0)),
        t.add_tvertex(p(w,h,0.0)), t.add_tvertex(p(0.0,h,0.0)),
        t.add_tvertex(p(0.0,0.0,d)), t.add_tvertex(p(w,0.0,d)),
        t.add_tvertex(p(w,h,d)), t.add_tvertex(p(0.0,h,d)),
    ];

    let e01 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(0.0,0.0,0.0), direction: x_axis })), v[0], v[1], [0.0, w]);
    let e12 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(w,0.0,0.0), direction: y_axis })), v[1], v[2], [0.0, h]);
    let e23 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(w,h,0.0), direction: -x_axis })), v[2], v[3], [0.0, w]);
    let e30 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(0.0,h,0.0), direction: -y_axis })), v[3], v[0], [0.0, h]);
    let e04 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(0.0,0.0,0.0), direction: z_axis })), v[0], v[4], [0.0, d]);
    let e15 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(w,0.0,0.0), direction: z_axis })), v[1], v[5], [0.0, d]);
    let e26 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(w,h,0.0), direction: z_axis })), v[2], v[6], [0.0, d]);
    let e37 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(0.0,h,0.0), direction: z_axis })), v[3], v[7], [0.0, d]);
    let e45 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(0.0,0.0,d), direction: x_axis })), v[4], v[5], [0.0, w]);
    let e56 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(w,0.0,d), direction: y_axis })), v[5], v[6], [0.0, h]);
    let e67 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(w,h,d), direction: -x_axis })), v[6], v[7], [0.0, w]);
    let e74 = t.add_tedge(Some(Curve3::Line(Line3 { origin: p(0.0,h,d), direction: -y_axis })), v[7], v[4], [0.0, h]);

    let wires = [
        t.add_twire(vec![e01, e12, e23, e30]),
        t.add_twire(vec![e45, e56, e67, e74]),
        t.add_twire(vec![e01, e15, rev(e45), rev(e04)]),
        t.add_twire(vec![rev(e23), e26, e67, rev(e37)]),
        t.add_twire(vec![rev(e30), e37, rev(e74), e04]),
        t.add_twire(vec![e12, e26, rev(e56), rev(e15)]),
    ];
    let normals = [-z_axis, z_axis, -y_axis, y_axis, -x_axis, x_axis];
    let centers = [
        p(w/2.0, h/2.0, 0.0), p(w/2.0, h/2.0, d),
        p(w/2.0, 0.0, d/2.0), p(w/2.0, h, d/2.0),
        p(0.0, h/2.0, d/2.0), p(w, h/2.0, d/2.0),
    ];

    let mut face_refs = Vec::new();
    for i in 0..6 {
        let surface = Surface3::Plane(Plane::new(centers[i], normals[i]));
        face_refs.push(t.add_tface(Some(surface), wires[i], vec![], Some(centers[i]), None, vec![], true));
    }
    let shell = t.add_tshell(face_refs);
    t.add_tsolid(vec![shell]);
    t.same_parameter();
    Ok(t)
}

pub fn sphere_primitive(radius: f64) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    Ok(PrimitiveSolid::Sphere { radius })
}

pub fn make_sphere_primitive(radius: f64) -> Result<PrimitiveSolid, BuildError> {
    sphere_primitive(radius)
}

pub fn sphere_brep(center: DVec3, radius: f64) -> Result<topods::BRep, BuildError> {
    let c = validate_point("center", center)?;
    let r = validate_positive("radius", radius)?;
    use rcad_kernel::topods::Orientation;
    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };

    let mut t = topods::BRep::new();
    let north = t.add_tvertex(c + DVec3::Z * r);
    let south = t.add_tvertex(c - DVec3::Z * r);

    // OCCT: sphere seam = circle on sphere surface (meridian at X=0 half-plane).
    // rcad was using Line through sphere center which causes false EE intersections.
    let seam_curve = Curve3::Circle(rcad_kernel::geom::Circle3 {
        center: c, radius: r,
        normal: DVec3::X,
        x_dir: DVec3::Z,
        y_dir: DVec3::Y,
    });
    // Degenerate edge at north pole (start == end)
    let e_top = t.add_tedge(None, north, north, [0.0, std::f64::consts::PI * r]);
    // Seam edge north→south (parameter range: half-circle, 0 = north, π = south)
    let e_seam = t.add_tedge(Some(seam_curve.clone()), north, south, [0.0, std::f64::consts::PI]);
    // Degenerate edge at south pole
    let e_bot = t.add_tedge(None, south, south, [0.0, std::f64::consts::PI * r]);

    let wire = t.add_twire(vec![e_top, e_seam, e_bot, rev(e_seam)]);
    let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface::new(c, DVec3::Y, r));
    let face = t.add_tface(Some(surface), wire, vec![], Some(c + DVec3::Z * r), None, vec![], true);
    let shell = t.add_tshell(vec![face]);
    t.add_tsolid(vec![shell]);
    t.same_parameter();
    Ok(t)
}

pub fn make_sphere_brep(center: DVec3, radius: f64) -> Result<topods::BRep, BuildError> {
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
) -> Result<topods::BRep, BuildError> {
    let c = validate_point("center", center)?;
    let r = validate_positive("radius", radius)?;
    let h = validate_positive("height", height)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    use rcad_kernel::topods::Orientation;
    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };
    let p = |dx: f64, dy: f64, dz: f64| c + x_axis * dx + y_axis * dy + z_axis * dz;

    let half = h * 0.5;
    let mut t = topods::BRep::new();
    let top_v = t.add_tvertex(p(0.0, half, 0.0));
    let bot_v = t.add_tvertex(p(0.0, -half, 0.0));

    // Degenerate edges at centers
    let e_top = t.add_tedge(None, top_v, top_v, [0.0, r * std::f64::consts::TAU]);
    let e_bot = t.add_tedge(None, bot_v, bot_v, [0.0, r * std::f64::consts::TAU]);
    // Seam edge
    let seam_curve = Curve3::Line(Line3 { origin: p(0.0, half, 0.0), direction: z_axis });
    let e_seam = t.add_tedge(Some(seam_curve), top_v, bot_v, [-half, half]);

    // Faces
    let cyl_surface = Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
        origin: c, axis: z_axis, ref_dir: x_axis, radius: r,
    });
    let lateral_wire = t.add_twire(vec![e_seam, rev(e_bot), rev(e_seam), e_top]);
    let top_wire = t.add_twire(vec![e_top]);
    let bot_wire = t.add_twire(vec![rev(e_bot)]);

    let cyl_surface = Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
        origin: c, axis: z_axis, ref_dir: x_axis, radius: r,
    });
    let lateral_face = t.add_tface(Some(cyl_surface), lateral_wire, vec![], Some(p(0.0, 0.0, r)), None, vec![], true);

    let top_plane = Surface3::Plane(Plane::new(p(0.0, half, 0.0), z_axis));
    let top_face = t.add_tface(Some(top_plane), top_wire, vec![], Some(p(0.0, half, 0.0)), None, vec![], false);

    let bot_plane = Surface3::Plane(Plane::new(p(0.0, -half, 0.0), -z_axis));
    let bot_face = t.add_tface(Some(bot_plane), bot_wire, vec![], Some(p(0.0, -half, 0.0)), None, vec![], false);

    let shell = t.add_tshell(vec![lateral_face, top_face, bot_face]);
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn make_cylinder_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
) -> Result<topods::BRep, BuildError> {
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
) -> Result<topods::BRep, BuildError> {
    let c = validate_point("center", center)?;
    let r = validate_positive("base_radius", base_radius)?;
    let h = validate_positive("height", height)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    use rcad_kernel::topods::Orientation;
    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };
    let p = |dx: f64, dy: f64, dz: f64| c + x_axis * dx + y_axis * dy + z_axis * dz;

    let half = h * 0.5;
    let mut t = topods::BRep::new();
    let apex = t.add_tvertex(p(0.0, half, 0.0));
    let base_c = t.add_tvertex(p(0.0, -half, 0.0));

    let e_apex = t.add_tedge(None, apex, apex, [0.0, r * std::f64::consts::TAU]);
    let e_base = t.add_tedge(None, base_c, base_c, [0.0, r * std::f64::consts::TAU]);
    let seam_curve = Curve3::Line(Line3 { origin: p(0.0, half, 0.0), direction: z_axis });
    let e_seam = t.add_tedge(Some(seam_curve), apex, base_c, [-half, half]);

    let top_wire = t.add_twire(vec![e_seam, rev(e_base), rev(e_seam), e_apex]);
    let bot_wire = t.add_twire(vec![e_base]);

    let cone_surf = Surface3::Cone(rcad_kernel::geom::ConicalSurface {
        apex: p(0.0, half, 0.0), axis: z_axis, radius: 0.0, half_angle_rad: (r / h).atan(),
    });
    let lateral = t.add_tface(Some(cone_surf), top_wire, vec![], Some(p(0.0, 0.0, r)), None, vec![], true);

    let bot_surf = Surface3::Plane(Plane::new(p(0.0, -half, 0.0), -z_axis));
    let bottom = t.add_tface(Some(bot_surf), bot_wire, vec![], Some(p(0.0, -half, 0.0)), None, vec![], false);

    let shell = t.add_tshell(vec![lateral, bottom]);
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn make_cone_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    base_radius: f64,
    height: f64,
) -> Result<topods::BRep, BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, height)
}

/// Right conical frustum (truncated cone), matching OCCT `pcone` when both bottom and top radii are positive.
///
/// `center` is the **midpoint** between the circular face centers; `axis` points from the bottom face
/// toward the top face; `r_bottom` / `r_top` are radii in those end planes; `height` is the distance
/// between the planes.
///
/// Built as an analytic BRep using `Surface3::Cone` for the lateral face and `Surface3::Plane`
/// for each cap, matching the topology of [`rcad_kernel::BRep::create_cone`].
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

    // Near-cylinder case: delegate to cylinder_brep
    if (rb - rt).abs() < 1e-12 {
        return cylinder_brep(center, axis, ref_dir, rb, height);
    }

    // Reuse the stable "narrower bottom 闁?wider top" cone parameterization for
    // the inverted frustum case by reversing the construction axis. This keeps
    // the cone V direction and apex-side semantics consistent instead of using a
    // separate mirrored branch that can drift in area/topology behavior.
    if rb > rt {
        return make_conical_frustum_brep(center, -axis, ref_dir, rt, rb, height);
    }

    use std::f64::consts::PI;
    use glam::DVec2;
    use rcad_kernel::geom::{self, *};
    use rcad_kernel::topods::Orientation;

    let half_angle = ((rb - rt).abs() / height).atan();
    let cos_ha = half_angle.cos();
    let tan_ha = half_angle.tan();
    let seam_len = height / cos_ha;

    let d_bottom = rb / tan_ha;
    let v_bottom = d_bottom / cos_ha;
    let v_top = (d_bottom + height) / cos_ha;

    let bottom_pt = DVec3::new(rb, -half_h, 0.0);
    let top_pt = DVec3::new(rt, half_h, 0.0);

    let mut t = topods::BRep::new();
    let v0 = t.add_tvertex(bottom_pt);
    let v1 = t.add_tvertex(top_pt);

    let bottom_circle = Curve3::Circle(Circle3::new(DVec3::new(0.0, -half_h, 0.0), -DVec3::Y, rb));
    let top_circle = Curve3::Circle(Circle3::new(DVec3::new(0.0, half_h, 0.0), DVec3::Y, rt));
    let seam_dir = (top_pt - bottom_pt).normalize();
    let seam_curve = Curve3::Line(Line3 { origin: bottom_pt, direction: seam_dir });

    let e0 = t.add_tedge(Some(bottom_circle), v0, v0, [0.0, 2.0 * PI]);
    let e1 = t.add_tedge(Some(top_circle), v1, v1, [0.0, 2.0 * PI]);
    let e2 = t.add_tedge(Some(seam_curve), v0, v1, [0.0, seam_len]);

        // Surfaces
    let apex = DVec3::new(0.0, -half_h - d_bottom, 0.0);
    let cone_surf = Surface3::Cone(geom::ConicalSurface {
        apex,
        axis: DVec3::Y,
        radius: 0.0,
        half_angle_rad: half_angle,
    });
    let bottom_plane = Surface3::Plane(Plane::new(DVec3::new(0.0, -half_h, 0.0), -DVec3::Y));
    let top_plane = Surface3::Plane(Plane::new(DVec3::new(0.0, half_h, 0.0), DVec3::Y));

    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };

    // Wires and faces
    let w0 = t.add_twire(vec![e0, e2, rev(e1), rev(e2)]);
    let f0 = t.add_tface(Some(cone_surf), w0, vec![], Some(DVec3::new(0.0, 0.0, rb)), None, vec![], true);

    let w1 = t.add_twire(vec![rev(e0)]);
    let f1 = t.add_tface(Some(bottom_plane), w1, vec![], Some(DVec3::new(0.0, -half_h, 0.0)), None, vec![], false);

    let w2 = t.add_twire(vec![e1]);
    let f2 = t.add_tface(Some(top_plane), w2, vec![], Some(DVec3::new(0.0, half_h, 0.0)), None, vec![], false);

    // PCurves on edges (keyed by face.index)
    let e0_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_bottom),
        direction: DVec2::new(1.0, 0.0),
    });
    let e0_on_plane = Curve2d::Circle(Circle2d::new(DVec2::ZERO, rb));
    let e1_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_top),
        direction: DVec2::new(1.0, 0.0),
    });
    let e1_on_plane = Curve2d::Circle(Circle2d::new(DVec2::ZERO, rt));
    let e2_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_bottom),
        direction: DVec2::new(0.0, (v_top - v_bottom) / seam_len),
    });

    t.edge_mut(e0).pcurves.insert(f0.index, (e0_on_cone, 0.0, 2.0 * PI));
    t.edge_mut(e0).pcurves.insert(f1.index, (e0_on_plane, 0.0, 2.0 * PI));
    t.edge_mut(e1).pcurves.insert(f0.index, (e1_on_cone, 0.0, 2.0 * PI));
    t.edge_mut(e1).pcurves.insert(f2.index, (e1_on_plane, 0.0, 2.0 * PI));
    t.edge_mut(e2).pcurves.insert(f0.index, (e2_on_cone, 0.0, seam_len));

    // Shell and solid
    let shell = t.add_tshell(vec![f0, f1, f2]);
    t.add_tsolid(vec![shell]);

    // Transform from local Y-up coordinates to target frame
    let mat = glam::DAffine3::from_cols(x_axis, y_axis, z_axis, center);
    t.apply_transform(mat);
    Ok(t)
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
) -> Result<topods::BRep, BuildError> {
    let c = validate_point("center", center)?;
    let maj = validate_positive("major_radius", major_radius)?;
    let min = validate_positive("minor_radius", minor_radius)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    use rcad_kernel::topods::Orientation;
    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };
    let p = |dx: f64, dy: f64, dz: f64| c + x_axis * dx + y_axis * dy + z_axis * dz;

    let mut t = topods::BRep::new();
    let seam_v = t.add_tvertex(p(maj + min, 0.0, 0.0));

    let major_circle = Curve3::Circle(rcad_kernel::geom::Circle3::new(c, z_axis, maj));
    let minor_circle = Curve3::Circle(rcad_kernel::geom::Circle3::new(p(maj, 0.0, 0.0), x_axis, min));

    let e_major = t.add_tedge(Some(major_circle), seam_v, seam_v, [0.0, std::f64::consts::TAU]);
    let e_minor = t.add_tedge(Some(minor_circle), seam_v, seam_v, [0.0, std::f64::consts::TAU]);

    let wire = t.add_twire(vec![e_major, e_minor, rev(e_major), rev(e_minor)]);
    let surface = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
        center: c, axis: z_axis, major_radius: maj, minor_radius: min,
    });
    let face = t.add_tface(Some(surface), wire, vec![], Some(p(maj + min, 0.0, 0.0)), None, vec![], true);
    let shell = t.add_tshell(vec![face]);
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn make_torus_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<topods::BRep, BuildError> {
    torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}

#[cfg(test)]
mod tests {
    use super::make_conical_frustum_brep;
    use glam::DVec3;
    use rcad_kernel::properties::{signed_volume, surface_area};

    #[test]
    fn conical_frustum_swapped_radii_preserve_area_and_volume() {
        let wide_bottom = make_conical_frustum_brep(
            DVec3::ZERO,
            DVec3::Z,
            DVec3::X,
            6.0,
            1.0,
            10.0,
        )
        .expect("wide-bottom frustum");
        let wide_top = make_conical_frustum_brep(
            DVec3::ZERO,
            DVec3::Z,
            DVec3::X,
            1.0,
            6.0,
            10.0,
        )
        .expect("wide-top frustum");

        let area_bottom = surface_area(&wide_bottom);
        let area_top = surface_area(&wide_top);
        assert!(
            (area_bottom - area_top).abs() <= 1e-6,
            "surface area mismatch: bottom={} top={}",
            area_bottom,
            area_top
        );

        let volume_bottom = signed_volume(&wide_bottom).abs();
        let volume_top = signed_volume(&wide_top).abs();
        assert!(
            (volume_bottom - volume_top).abs() <= 1e-6,
            "volume mismatch: bottom={} top={}",
            volume_bottom,
            volume_top
        );
    }
}

/// Create a BRep with a single rectangular planar face.
///
/// `origin` is the plane origin; `u_axis`/`v_axis` are the parametric axes
/// (the normal is `u_axis 閼?v_axis`). The corners are `origin + u_axis*u + v_axis*v`.
pub fn make_planar_rect_brep(
    origin: DVec3,
    u_axis: DVec3,
    v_axis: DVec3,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> Result<topods::BRep, BuildError> {
    let normal = u_axis.cross(v_axis).normalize();
    let surface = Surface3::Plane(Plane::new(origin, normal));

    let c0 = origin + u_axis * umin + v_axis * vmin;
    let c1 = origin + u_axis * umax + v_axis * vmin;
    let c2 = origin + u_axis * umax + v_axis * vmax;
    let c3 = origin + u_axis * umin + v_axis * vmax;

    let mut t = topods::BRep::new();
    let v0 = t.add_tvertex(c0);
    let v1 = t.add_tvertex(c1);
    let v2 = t.add_tvertex(c2);
    let v3 = t.add_tvertex(c3);

    let dir01 = (c1 - c0).normalize();
    let dir12 = (c2 - c1).normalize();
    let dir23 = (c3 - c2).normalize();
    let dir30 = (c0 - c3).normalize();

    let len01 = (c1 - c0).length();
    let len12 = (c2 - c1).length();
    let len23 = (c3 - c2).length();
    let len30 = (c0 - c3).length();

    let e0 = t.add_tedge(Some(Curve3::Line(Line3 { origin: c0, direction: dir01 })), v0, v1, [0.0, len01]);
    let e1 = t.add_tedge(Some(Curve3::Line(Line3 { origin: c1, direction: dir12 })), v1, v2, [0.0, len12]);
    let e2 = t.add_tedge(Some(Curve3::Line(Line3 { origin: c2, direction: dir23 })), v2, v3, [0.0, len23]);
    let e3 = t.add_tedge(Some(Curve3::Line(Line3 { origin: c3, direction: dir30 })), v3, v0, [0.0, len30]);

    let wire = t.add_twire(vec![e0, e1, e2, e3]);
    t.add_tface(Some(surface), wire, vec![], None, None, vec![], true);
    Ok(t)
}

/// Create a BRep with a single planar face bounded by a polygon, with optional inner holes.
///
/// `origin` and `normal` define the plane. `outer_polygon` gives the outer boundary
/// vertices (闁? points in CCW order). `inner_polygons` are hole boundaries (each 闁? points
/// in CW order).
pub fn make_planar_polygon_brep(
    origin: DVec3,
    normal: DVec3,
    outer_polygon: &[DVec3],
    inner_polygons: &[Vec<DVec3>],
) -> Result<topods::BRep, BuildError> {
    let surface = Surface3::Plane(Plane::new(origin, normal));
    let mut t = topods::BRep::new();
    let outer_wire = make_polygon_wire_topods(&mut t, outer_polygon)?;
    let mut inner_wires = Vec::new();
    for poly in inner_polygons {
        inner_wires.push(make_polygon_wire_topods(&mut t, poly)?);
    }
    t.add_tface(Some(surface), outer_wire, inner_wires, None, None, vec![], true);
    Ok(t)
}

fn make_polygon_wire_topods(t: &mut topods::BRep, points: &[DVec3]) -> Result<topods::ShapeRef, BuildError> {
    let n = points.len();
    if n < 3 { return Err(BuildError::DegenerateGeometry("polygon needs at least 3 points")); }
    let verts: Vec<_> = points.iter().map(|&p| t.add_tvertex(p)).collect();
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let dir = (points[j] - points[i]).normalize();
        let len = (points[j] - points[i]).length();
        edges.push(t.add_tedge(Some(Curve3::Line(Line3 { origin: points[i], direction: dir })), verts[i], verts[j], [0.0, len]));
    }
    Ok(t.add_twire(edges))
}

/// Create a half-space solid bounded by a plane: a large box extending from `origin` in the
/// `-normal` direction (the "interior" side). Used by OCCT `mkvolume` translation to intersect
/// planar half-spaces.
///
/// `size` controls the extents: 閸楊槞ize in the plane directions, size in the -normal direction.
pub fn make_half_space_brep(
    origin: DVec3,
    normal: DVec3,
    size: f64,
) -> Result<topods::BRep, BuildError> {
    let normal = normalize_vector("normal", normal)?;
    let origin = validate_point("origin", origin)?;
    let size = validate_positive("size", size)?;

    // Choose u_axis perpendicular to normal
    let abs = normal.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z {
        DVec3::X
    } else if abs.y <= abs.z {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let u_axis = normal.cross(candidate).normalize();
    // v_axis = u_axis 閼?normal gives u_axis 閼?v_axis = -normal (the interior direction).
    let v_axis = u_axis.cross(normal).normalize();

    // box_brep creates a corner-box at `origin - u_axis*size - v_axis*size` that extends
    // +u_axis*2*size, +v_axis*2*size, +z_axis*size, giving 閸楊槞ize UV coverage around `origin`
    // and `size` units in the -normal (interior) direction.
    let corner = origin - u_axis * size - v_axis * size;
    box_brep(corner, u_axis, v_axis, 2.0 * size, 2.0 * size, size)
}

/// Build a convex polyhedron BRep directly from half-space plane equations (no boolean ops).
///
/// Each half-space is `(origin, normal)` where the interior is the side opposite the normal,
/// i.e. `normal 鐠?(p - origin) 闁?0`. Input planes may contain coplanar duplicates 闁?they
/// are merged into a single face by keeping the tightest constraint.
///
/// The algorithm computes all vertices by intersecting 3 planes at a time, keeps those
// satisfying all half-spaces, then builds faces from coplanar vertex sets.
pub fn make_convex_polyhedron_from_half_spaces(
    planes: &[(DVec3, DVec3)],
) -> Result<BRep, BuildError> {
    if planes.len() < 4 {
        return Err(BuildError::DegenerateGeometry(
            "need at least 4 planes for a convex polyhedron",
        ));
    }

    // Geometry-adaptive tolerance
    let extent: f64 = planes
        .iter()
        .map(|(o, _)| o.length_squared())
        .fold(0.0, f64::max)
        .sqrt()
        .max(1.0);
    let tol = extent * 1e-12 + 1e-10;

    // --- 1. Deduplicate coplanar planes, keep tightest constraint ---
    struct PlaneEq {
        n: DVec3,
        d: f64,
        origin: DVec3,
    }
    let mut eqs: Vec<PlaneEq> = Vec::new();

    for &(origin, normal) in planes {
        let n = normal.normalize();
        let d = n.dot(origin);
        let mut found = false;
        for eq in &mut eqs {
            let cos = n.dot(eq.n);
            if cos > 1.0 - tol {
                // Same direction 闁?keep the tighter constraint (smaller d for n鐠虹棷 闁?d)
                if d < eq.d {
                    eq.d = d;
                    eq.origin = origin;
                }
                found = true;
                break;
            }
        }
        if !found {
            eqs.push(PlaneEq { n, d, origin });
        }
    }

    if eqs.len() < 4 {
        return Err(BuildError::DegenerateGeometry(
            "fewer than 4 unique planes after deduplication",
        ));
    }

    // --- 2. Compute all intersection points of 3 planes ---
    let mut verts: Vec<DVec3> = Vec::new();
    let np = eqs.len();

    for i in 0..np {
        for j in (i + 1)..np {
            for k in (j + 1)..np {
                let (ni, nj, nk) = (eqs[i].n, eqs[j].n, eqs[k].n);
                // Rows are plane normals (ni鐠虹棷=di, 闁?; columns store x/y/z components per plane.
                let m = DMat3::from_cols(
                    DVec3::new(ni.x, nj.x, nk.x),
                    DVec3::new(ni.y, nj.y, nk.y),
                    DVec3::new(ni.z, nj.z, nk.z),
                );
                if m.determinant().abs() < tol {
                    continue; // singular 闁?planes meet in a line or are parallel
                }
                let p = m.inverse() * DVec3::new(eqs[i].d, eqs[j].d, eqs[k].d);

                // Check against all half-space constraints
                let valid = eqs.iter().all(|eq| eq.n.dot(p) - eq.d <= tol * 100.0);
                if !valid {
                    continue;
                }

                // Deduplicate
                let dup = verts.iter().any(|v| (v - p).length_squared() < tol * tol);
                if !dup {
                    verts.push(p);
                }
            }
        }
    }

    if verts.len() < 4 {
        return Err(BuildError::DegenerateGeometry(
            "fewer than 4 valid vertices 闁?half-spaces may not form a bounded convex polyhedron",
        ));
    }

    // --- 3. Group vertices by plane ---
    let mut face_verts: Vec<Vec<usize>> = vec![Vec::new(); np];
    for (vi, &v) in verts.iter().enumerate() {
        for (ei, eq) in eqs.iter().enumerate() {
            if (eq.n.dot(v) - eq.d).abs() < tol * 100.0 {
                face_verts[ei].push(vi);
            }
        }
    }

    // --- 4. Build BRep ---
    use super::brep_builder::*;
    let mut brep = BRep::default();

    // Create all vertices upfront so edges on adjacent faces share indices
    let bv: Vec<usize> = verts.iter().map(|&v| make_vertex(&mut brep, v)).collect();

    for fi in 0..np {
        if face_verts[fi].len() < 3 {
            continue;
        }
        let n = eqs[fi].n;
        let origin = eqs[fi].origin;

        // Centroid of face vertices
        let centroid: DVec3 = face_verts[fi]
            .iter()
            .map(|&vi| verts[vi])
            .sum::<DVec3>()
            / face_verts[fi].len() as f64;

        // UV basis in the plane (u 閼?v = n for CCW winding from +n)
        let abs = n.abs();
        let candidate = if abs.x <= abs.y && abs.x <= abs.z {
            DVec3::X
        } else if abs.y <= abs.z {
            DVec3::Y
        } else {
            DVec3::Z
        };
        let u_axis = n.cross(candidate).normalize();
        let v_axis = n.cross(u_axis).normalize(); // u 閼?v = n

        // Sort vertices by angle around centroid in the UV plane
        let mut sorted: Vec<usize> = face_verts[fi].clone();
        sorted.sort_by(|&a, &b| {
            let da = verts[a] - centroid;
            let db = verts[b] - centroid;
            let aa = da.dot(u_axis).atan2(da.dot(v_axis));
            let ab = db.dot(u_axis).atan2(db.dot(v_axis));
            aa.partial_cmp(&ab).unwrap()
        });

        // Build wire
        let mut wire_edges = Vec::with_capacity(sorted.len());
        for w in 0..sorted.len() {
            let vi0 = bv[sorted[w]];
            let vi1 = bv[sorted[(w + 1) % sorted.len()]];
            let p0 = verts[sorted[w]];
            let p1 = verts[sorted[(w + 1) % sorted.len()]];
            let dir = (p1 - p0).normalize();
            let len = (p1 - p0).length();
            let e = make_edge(
                &mut brep,
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: dir,
                }),
                0.0,
                len,
                vi0,
                vi1,
            )?;
            wire_edges.push(WireEdge::new(e, true));
        }

        let wire = make_wire(wire_edges);
        let surface = Surface3::Plane(Plane::new(origin, n));
        make_face(&mut brep, surface, wire, vec![])?;
    }

    Ok(brep)
}

/// Mirror a BRep across a plane defined by `origin` and `normal`.
///
/// The mirrored BRep has inverted face normals and reversed wire orientations
/// to maintain consistent outward-facing normals.
pub fn mirror_brep(brep: &topods::BRep, plane_origin: DVec3, plane_normal: DVec3) -> Result<topods::BRep, BuildError> {
    let _ = validate_point("plane_origin", plane_origin)?;
    let n = normalize_vector("plane_normal", plane_normal)?;
    let mirror_p = |p: DVec3| { let d = (p - plane_origin).dot(n); p - n * (2.0 * d) };
    let mirror_v = |v: DVec3| v - n * (2.0 * v.dot(n));

    let mut out = topods::BRep::new();
    use rcad_kernel::topods::{ShapeRef, Orientation};
    use rcad_kernel::topods::TShape;
    let rev = |sr: ShapeRef| ShapeRef { orientation: Orientation::Reversed, ..sr };

    // Phase 1: mirror all vertices, build ptr_id → new ShapeRef map
    let mut vmap: std::collections::HashMap<u64, ShapeRef> = std::collections::HashMap::new();
    for ts in &brep.tshapes {
        if let TShape::Vertex(vd) = &**ts {
            let sr = out.add_tvertex(mirror_p(vd.point));
            vmap.insert(tshape_ptr_id(ts), sr);
        }
    }

    // Phase 2: mirror curves (transform each curve component)
    let mirror_curve = |c: &Curve3| -> Curve3 {
        match c {
            Curve3::Line(l) => Curve3::Line(Line3 { origin: mirror_p(l.origin), direction: mirror_v(l.direction) }),
            Curve3::Circle(c3) => Curve3::Circle(rcad_kernel::geom::Circle3::new(mirror_p(c3.center), mirror_v(c3.normal), c3.radius)),
            _ => c.clone(), // fallback for complex curves — clone as-is
        }
    };

    // Phase 3: create edges, map old edge ptr_id → new ShapeRef
    let mut emap: std::collections::HashMap<u64, ShapeRef> = std::collections::HashMap::new();
    for ts in &brep.tshapes {
        if let TShape::Edge(ed) = &**ts {
            let first = *vmap.get(&ed.first.ptr_id).unwrap_or(&ed.first);
            let last = *vmap.get(&ed.last.ptr_id).unwrap_or(&ed.last);
            let curve = ed.curve.as_ref().map(mirror_curve);
            let sr = out.add_tedge(curve, first, last, ed.range);
            emap.insert(tshape_ptr_id(ts), sr);
        }
    }

    // Phase 4: create wires, faces, shells, solids — traverse in order
    // topods::BRep stores TShapes in flat order. We iterate and build mirror
    // structures, storing new ShapeRef for wires and faces as we go.
    let mut wmap: std::collections::HashMap<u64, ShapeRef> = std::collections::HashMap::new();
    let mut fmap: std::collections::HashMap<u64, ShapeRef> = std::collections::HashMap::new();
    let mut shmap: std::collections::HashMap<u64, ShapeRef> = std::collections::HashMap::new();

    for ts in &brep.tshapes {
        match &**ts {
            TShape::Wire(wd) => {
                let mirrored: Vec<ShapeRef> = wd.edges.iter().map(|&sr| rev(sr)).collect();
                let sr = out.add_twire(mirrored);
                wmap.insert(tshape_ptr_id(ts), sr);
            }
            TShape::Face(fd) => {
                let outer = *wmap.get(&fd.outer_wire.ptr_id).unwrap_or(&fd.outer_wire);
                let inner: Vec<ShapeRef> = fd.inner_wires.iter().map(|sr| *wmap.get(&sr.ptr_id).unwrap_or(sr)).collect();
                let surface = fd.surface.as_ref().map(|s| mirror_surface(s, &mirror_p, &mirror_v));
                let sr = out.add_tface(surface, outer, inner, fd.sample_point, fd.uv_domain, fd.internal_vertices.clone(), fd.natural_restriction);
                fmap.insert(tshape_ptr_id(ts), sr);
            }
            TShape::Shell(sd) => {
                let faces: Vec<ShapeRef> = sd.faces.iter().map(|sr| *fmap.get(&sr.ptr_id).unwrap_or(sr)).collect();
                let sr = out.add_tshell(faces);
                shmap.insert(tshape_ptr_id(ts), sr);
            }
            TShape::Solid(sd) => {
                let shells: Vec<ShapeRef> = sd.shells.iter().map(|sr| *shmap.get(&sr.ptr_id).unwrap_or(sr)).collect();
                out.add_tsolid(shells);
            }
            _ => {}
        }
    }
    Ok(out)
}

fn tshape_ptr_id(ts: &std::sync::Arc<rcad_kernel::topods::TShape>) -> u64 {
    std::sync::Arc::as_ptr(ts) as u64
}

fn mirror_surface(s: &Surface3, mirror_p: impl Fn(DVec3) -> DVec3, mirror_v: impl Fn(DVec3) -> DVec3) -> Surface3 {
    match s {
        Surface3::Plane(p) => Surface3::Plane(Plane::new(mirror_p(p.origin), mirror_v(p.normal))),
        Surface3::Cylinder(c) => Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: mirror_p(c.origin), axis: mirror_v(c.axis), radius: c.radius, ref_dir: mirror_v(c.ref_dir),
        }),
        Surface3::Sphere(s) => Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: mirror_p(s.center), axis: mirror_v(s.axis), radius: s.radius, ref_dir: mirror_v(s.ref_dir),
        }),
        Surface3::Cone(c) => Surface3::Cone(rcad_kernel::geom::ConicalSurface {
            apex: mirror_p(c.apex), axis: mirror_v(c.axis), radius: c.radius, half_angle_rad: c.half_angle_rad,
        }),
        Surface3::Torus(t) => Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: mirror_p(t.center), axis: mirror_v(t.axis), major_radius: t.major_radius, minor_radius: t.minor_radius,
        }),
        _ => s.clone(),
    }
}

pub fn make_conical_frustum_brep_topods(center: DVec3, axis: DVec3, ref_dir: DVec3, r_bottom: f64, r_top: f64, height: f64) -> Result<topods::BRep, BuildError> {
    make_conical_frustum_brep(center, axis, ref_dir, r_bottom, r_top, height)
}

pub fn make_convex_polyhedron_from_half_spaces_topods(planes: &[(DVec3, DVec3)]) -> Result<topods::BRep, BuildError> {
    make_convex_polyhedron_from_half_spaces(planes)
}

// 闁冲厜鍋撻柍鍏夊亾 Topods-native wrappers 闁冲厜鍋撻柍鍏夊亾

