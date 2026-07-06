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

pub fn sphere_brep(center: DVec3, radius: f64) -> Result<topods::BRep, BuildError> {
    let c = validate_point("center", center)?;
    let r = validate_positive("radius", radius)?;
    use rcad_kernel::topods::Orientation;
    let rev = |sr: rcad_kernel::topods::ShapeRef| rcad_kernel::topods::ShapeRef { orientation: Orientation::Reversed, ..sr };

    let mut t = topods::BRep::new();
    let north = t.add_tvertex(c + DVec3::Z * r);
    let south = t.add_tvertex(c - DVec3::Z * r);

    let seam_curve = Curve3::Line(Line3 { origin: c, direction: DVec3::Z });
    // Degenerate edge at north pole (start == end)
    let e_top = t.add_tedge(None, north, north, [0.0, std::f64::consts::PI * r]);
    // Seam edge north→south
    let e_seam = t.add_tedge(Some(seam_curve.clone()), north, south, [-r, r]);
    // Degenerate edge at south pole
    let e_bot = t.add_tedge(None, south, south, [0.0, std::f64::consts::PI * r]);

    let wire = t.add_twire(vec![e_top, e_seam, e_bot, rev(e_seam)]);
    let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface::new(c, DVec3::Y, r));
    let face = t.add_tface(Some(surface), wire, vec![], Some(c + DVec3::Z * r), None, vec![], true);
    let shell = t.add_tshell(vec![face]);
    t.add_tsolid(vec![shell]);
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

    // Near-cylinder case: fall back to kernel cylinder primitive
    if (rb - rt).abs() < 1e-12 {
        let primitive = PrimitiveSolid::Cylinder {
            radius: rb,
            height,
        };
        let mut brep = BRep::from_primitive(primitive);
        transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
        return Ok(brep);
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
    use rcad_kernel::{
        geom::{self, *},
        Edge, Face, GeomStore, PCurve, Shell, Solid, Vertex, Wire, WireEdge,
    };

    let half_angle = ((rb - rt).abs() / height).atan();
    let cos_ha = half_angle.cos();
    let tan_ha = half_angle.tan();
    let seam_len = height / cos_ha;

    // Compute apex position and cone orientation.
    // The apex is on the central axis, on the side of the narrower end.
    let d_bottom = rb / tan_ha;
    let apex_y = -half_h - d_bottom;
    let axis_dir = DVec3::Y;
    let v_bottom = d_bottom / cos_ha;
    let v_top = (d_bottom + height) / cos_ha;

    // Vertices (in local coordinates: axis = Y, bottom at -half_h, top at +half_h)
    let bottom_pt = DVec3::new(rb, -half_h, 0.0);
    let top_pt = DVec3::new(rt, half_h, 0.0);
    let vertices = vec![Vertex { point: bottom_pt }, Vertex { point: top_pt }];

    // Edges: E0 = bottom circle, E1 = top circle, E2 = seam
    let edges = vec![
        Edge { start: 0, end: 0 },
        Edge { start: 1, end: 1 },
        Edge { start: 0, end: 1 },
    ];

    // 闁冲厜鍋撻柍鍏夊亾 3D curves 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋?
    let bottom_circle = Curve3::Circle(Circle3::new(DVec3::new(0.0, -half_h, 0.0), -DVec3::Y, rb));
    let top_circle = Curve3::Circle(Circle3::new(DVec3::new(0.0, half_h, 0.0), DVec3::Y, rt));
    let seam_curve = Curve3::Line(Line3 {
        origin: bottom_pt,
        direction: (top_pt - bottom_pt).normalize(),
    });

    // 闁冲厜鍋撻柍鍏夊亾 Surfaces 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾
    let apex = DVec3::new(0.0, apex_y, 0.0);
    let cone_surf = Surface3::Cone(geom::ConicalSurface {
        apex,
        axis: axis_dir,
        radius: 0.0,
        half_angle_rad: half_angle,
    });
    let bottom_plane = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, -half_h, 0.0),
        normal: -DVec3::Y,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, half_h, 0.0),
        normal: DVec3::Y,
    });

    // 闁冲厜鍋撻柍鍏夊亾 PCurves 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋?
    // E0 bottom circle on cone face: iso-V at V=v_bottom
    let e0_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_bottom),
        direction: DVec2::new(1.0, 0.0),
    });
    // E0 on bottom plane: full circle in UV
    let e0_on_plane = Curve2d::Circle(Circle2d::new(DVec2::ZERO, rb));
    // E1 top circle on cone face: iso-V at V=v_top
    let e1_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_top),
        direction: DVec2::new(1.0, 0.0),
    });
    // E1 on top plane: full circle in UV
    let e1_on_plane = Curve2d::Circle(Circle2d::new(DVec2::ZERO, rt));
    // E2 seam on cone face: iso-U at U=0, V ranges v_bottom 闁?v_top
    let e2_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_bottom),
        direction: DVec2::new(0.0, (v_top - v_bottom) / seam_len),
    });

    // 闁冲厜鍋撻柍鍏夊亾 Faces 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋?
    // F0 lateral (cone): bottom circle fwd 闁?seam fwd 闁?top circle rev 闁?seam rev
    let f0 = Face {
        outer_wire: Wire {
            edges: vec![
                WireEdge::fwd(0),
                WireEdge::fwd(2),
                WireEdge::rev(1),
                WireEdge::rev(2),
            ],
        },
        inner_wires: vec![],
        normal: DVec3::X,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };
    // F1 bottom cap (plane): E0 rev (CCW when viewed from -Y)
    let f1 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::rev(0)],
        },
        inner_wires: vec![],
        normal: -DVec3::Y,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };
    // F2 top cap (plane): E1 fwd (CCW when viewed from +Y)
    let f2 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(1)],
        },
        inner_wires: vec![],
        normal: DVec3::Y,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };

    let solid = Solid {
        shells: vec![Shell {
            faces: vec![f0, f1, f2],
        }],
    };

    // 闁冲厜鍋撻柍鍏夊亾 Geometry store 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾
    let geom = GeomStore { face_internal_vertices: vec![], edge_vertex_params: vec![],
        curves: vec![bottom_circle, top_circle, seam_curve],
        surfaces: vec![cone_surf, bottom_plane, top_plane],
        curve2ds: vec![e0_on_cone, e0_on_plane, e1_on_cone, e1_on_plane, e2_on_cone],
        edge_curve: vec![Some(0), Some(1), Some(2)],
        face_surface: vec![Some(0), Some(1), Some(2)],
        edge_pcurves: vec![
            // E0: on cone (surface 0) + bottom plane (surface 1)
            vec![
                PCurve { surface_idx: 0, curve2d_idx: 0 },
                PCurve { surface_idx: 1, curve2d_idx: 1 },
            ],
            // E1: on cone (surface 0) + top plane (surface 2)
            vec![
                PCurve { surface_idx: 0, curve2d_idx: 2 },
                PCurve { surface_idx: 2, curve2d_idx: 3 },
            ],
            // E2 seam: on cone only
            vec![PCurve { surface_idx: 0, curve2d_idx: 4 }],
        ],
        edge_curve_range: vec![
            Some([0.0, 2.0 * PI]),
            Some([0.0, 2.0 * PI]),
            Some([0.0, seam_len]),
        ],
        edge_degenerated: vec![false, false, false],
        vertex_tolerance: Vec::new(),
        edge_tolerance: Vec::new(),
        face_tolerance: Vec::new(),
        curve2d_range: Vec::new(),
        face_surface_range: Vec::new(),
        edge_same_parameter: Vec::new(),
        edge_same_range: Vec::new(),
    };

    let mut brep = BRep {
        vertices,
        edges,
        solids: vec![solid],
        geom,
        compound: None,
        compsolid: None,
    };

    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
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
/// vertices (闁? points in CCW order). `inner_polygons` are hole boundaries (each 闁? points
/// in CW order).
pub fn make_planar_polygon_brep(
    origin: DVec3,
    normal: DVec3,
    outer_polygon: &[DVec3],
    inner_polygons: &[Vec<DVec3>],
) -> Result<topods::BRep, BuildError> {
    let surface = Surface3::Plane(Plane { origin, normal });
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
) -> Result<BRep, BuildError> {
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
        let surface = Surface3::Plane(Plane { origin, normal: n });
        make_face(&mut brep, surface, wire, vec![])?;
    }

    Ok(brep)
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

// 闁冲厜鍋撻柍鍏夊亾 Topods-native wrappers 闁冲厜鍋撻柍鍏夊亾

