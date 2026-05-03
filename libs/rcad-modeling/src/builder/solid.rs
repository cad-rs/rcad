use super::{
    BuildError, basis_from_axis_ref, basis_from_x_y, do_mirror_brep, normalize_vector,
    transform_brep, translate_brep, validate_point, validate_positive,
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
    let (apex_y, axis_dir, v_bottom, v_top) = if rt > rb {
        // Wider at top → apex below bottom, axis = +Y
        let d_bottom = rb / tan_ha;
        let ay = -half_h - d_bottom;
        (ay, DVec3::Y, d_bottom / cos_ha, (d_bottom + height) / cos_ha)
    } else {
        // Wider at bottom → apex above top, axis = -Y
        let d_top = rt / tan_ha;
        let ay = half_h + d_top;
        (ay, -DVec3::Y, (d_top + height) / cos_ha, d_top / cos_ha)
    };

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

    // ── 3D curves ───────────────────────────────────────────────
    let bottom_circle = Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, -half_h, 0.0),
        normal: -DVec3::Y,
        radius: rb,
    });
    let top_circle = Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, half_h, 0.0),
        normal: DVec3::Y,
        radius: rt,
    });
    let seam_curve = Curve3::Line(Line3 {
        origin: bottom_pt,
        direction: (top_pt - bottom_pt).normalize(),
    });

    // ── Surfaces ────────────────────────────────────────────────
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

    // ── PCurves ─────────────────────────────────────────────────
    // E0 bottom circle on cone face: iso-V at V=v_bottom
    let e0_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_bottom),
        direction: DVec2::new(1.0, 0.0),
    });
    // E0 on bottom plane: full circle in UV
    let e0_on_plane = Curve2d::Circle(Circle2d {
        center: DVec2::ZERO,
        radius: rb,
    });
    // E1 top circle on cone face: iso-V at V=v_top
    let e1_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_top),
        direction: DVec2::new(1.0, 0.0),
    });
    // E1 on top plane: full circle in UV
    let e1_on_plane = Curve2d::Circle(Circle2d {
        center: DVec2::ZERO,
        radius: rt,
    });
    // E2 seam on cone face: iso-U at U=0, V ranges v_bottom → v_top
    let e2_on_cone = Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, v_bottom),
        direction: DVec2::new(0.0, (v_top - v_bottom) / seam_len),
    });

    // ── Faces ───────────────────────────────────────────────────
    // F0 lateral (cone): bottom circle fwd → seam fwd → top circle rev → seam rev
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
        mesh_dirty: true,
    };
    // F1 bottom cap (plane): E0 rev (CCW when viewed from -Y)
    let f1 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::rev(0)],
        },
        inner_wires: vec![],
        normal: -DVec3::Y,
        triangles: vec![],
        mesh_dirty: true,
    };
    // F2 top cap (plane): E1 fwd (CCW when viewed from +Y)
    let f2 = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(1)],
        },
        inner_wires: vec![],
        normal: DVec3::Y,
        triangles: vec![],
        mesh_dirty: true,
    };

    let solid = Solid {
        shells: vec![Shell {
            faces: vec![f0, f1, f2],
        }],
    };

    // ── Geometry store ──────────────────────────────────────────
    let geom = GeomStore {
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

/// Mirror a BRep across a plane defined by `origin` and `normal`.
///
/// The mirrored BRep has inverted face normals and reversed wire orientations
/// to maintain consistent outward-facing normals.
pub fn mirror_brep(brep: &BRep, plane_origin: DVec3, plane_normal: DVec3) -> Result<BRep, BuildError> {
    let _ = validate_point("plane_origin", plane_origin)?;
    let n = normalize_vector("plane_normal", plane_normal)?;
    Ok(do_mirror_brep(brep, plane_origin, n))
}
