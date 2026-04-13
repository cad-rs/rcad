//! Array (pattern) operations — linear and circular repetition of BRep solids.
//!
//! Analogous to OCCT `BRepOffsetAPI_MakeThickSolid`-style patterns and
//! `BRepFeat_MakeLinearForm` / `BRepFeat_MakeRevol` for feature repetition.
//!
//! # Operations
//!
//! - **Linear pattern**: repeat along a direction with uniform spacing
//! - **Circular pattern**: rotate around an axis with uniform angular spacing

use glam::{DMat4, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{
    Circle3, ConicalSurface, Curve3, CylindricalSurface, Ellipse3, Hyperbola3, Line3,
    LinearExtrusionSurface, OffsetSurface, Plane, RevolutionSurface, SphericalSurface, Surface3,
    ToroidalSurface, TrimmedSurface,
};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

/// Parameters for a linear pattern.
#[derive(Debug, Clone)]
pub struct LinearPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// Number of copies (including the original). Must be >= 1.
    pub count: usize,
    /// Spacing between consecutive copies.
    pub spacing: f64,
}

/// Parameters for a circular pattern.
#[derive(Debug, Clone)]
pub struct CircularPatternParams {
    /// A point on the rotation axis.
    pub axis_origin: DVec3,
    /// Normalized rotation axis direction.
    pub axis_direction: DVec3,
    /// Number of copies (including the original). Must be >= 1.
    pub count: usize,
    /// Total angle in radians for the full pattern (copies are evenly spaced).
    pub total_angle: f64,
}

/// Error type for pattern operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PatternError {
    /// Count must be at least 1.
    InvalidCount,
    /// Spacing must be positive.
    InvalidSpacing,
    /// Direction vector must be non-zero.
    ZeroDirection,
    /// Axis direction must be non-zero.
    ZeroAxis,
    /// Total angle must be non-zero and <= 2*pi.
    InvalidAngle,
    /// Input BRep has no solids.
    NoSolids,
}

impl std::fmt::Display for PatternError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCount => write!(f, "pattern count must be >= 1"),
            Self::InvalidSpacing => write!(f, "pattern spacing must be > 0"),
            Self::ZeroDirection => write!(f, "pattern direction must be non-zero"),
            Self::ZeroAxis => write!(f, "pattern axis must be non-zero"),
            Self::InvalidAngle => write!(f, "pattern angle must be > 0 and <= 2*pi"),
            Self::NoSolids => write!(f, "input BRep has no solids"),
        }
    }
}

/// Apply a linear pattern to a BRep — repeat copies along a direction.
///
/// Returns a new BRep containing all copies merged into a single solid.
/// The original is included as the first copy (offset 0).
pub fn linear_pattern(
    brep: &BRep,
    params: &LinearPatternParams,
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }
    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    for i in 0..params.count {
        let offset = dir * (i as f64 * params.spacing);
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

/// Apply a circular pattern to a BRep — rotate copies around an axis.
///
/// Returns a new BRep containing all copies merged into a single solid.
/// The original is included as the first copy (angle 0).
pub fn circular_pattern(
    brep: &BRep,
    params: &CircularPatternParams,
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_angle <= 0.0 || params.total_angle > std::f64::consts::TAU {
        return Err(PatternError::InvalidAngle);
    }
    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();
    let angle_step = params.total_angle / params.count as f64;

    for i in 0..params.count {
        let angle = i as f64 * angle_step;
        let mat = rotation_matrix(params.axis_origin, axis, angle);
        append_transformed_brep(&mut out, brep, &mat)?;
    }

    Ok(out)
}

// ── Internal helpers ─────────────────────────────────────────────────────────

fn translation_matrix(offset: DVec3) -> DMat4 {
    DMat4::from_translation(offset)
}

fn rotation_matrix(origin: DVec3, axis: DVec3, angle: f64) -> DMat4 {
    DMat4::from_translation(origin)
        * DMat4::from_axis_angle(axis, angle)
        * DMat4::from_translation(-origin)
}

fn append_transformed_brep(
    target: &mut BRep,
    source: &BRep,
    mat: &DMat4,
) -> Result<(), PatternError> {
    let v_offset = target.vertices.len();
    let e_offset = target.edges.len();
    let curve_offset = target.geom.curves.len();
    let surface_offset = target.geom.surfaces.len();

    // Transform and copy vertices
    for v in &source.vertices {
        let p = mat.transform_point3(v.point.into());
        target.vertices.push(Vertex {
            point: DVec3::new(p.x, p.y, p.z),
        });
    }

    // Transform and copy curves
    for curve in &source.geom.curves {
        target.geom.curves.push(transform_curve(curve, mat));
    }

    // Transform and copy surfaces
    for surface in &source.geom.surfaces {
        target.geom.surfaces.push(transform_surface(surface, mat));
    }

    // Copy edges with remapped vertex indices
    for e in &source.edges {
        target.edges.push(Edge {
            start: e.start + v_offset,
            end: e.end + v_offset,
        });
    }

    // Remap edge geometry references
    for &ec in &source.geom.edge_curve {
        target
            .geom
            .edge_curve
            .push(ec.map(|c| c + curve_offset));
    }
    for &ecr in &source.geom.edge_curve_range {
        target.geom.edge_curve_range.push(ecr);
    }
    for &ed in &source.geom.edge_degenerated {
        target.geom.edge_degenerated.push(ed);
    }

    // Copy faces with remapped edge indices
    for solid in &source.solids {
        let mut shell = Shell { faces: Vec::new() };
        for face in &solid.shells[0].faces {
            let wire_edges: Vec<WireEdge> = face
                .outer_wire
                .edges
                .iter()
                .map(|we| WireEdge {
                    idx: we.idx + e_offset,
                    forward: we.forward,
                })
                .collect();

            // Remap inner wire edge indices too
            let inner_wires: Vec<Wire> = face
                .inner_wires
                .iter()
                .map(|wire| Wire {
                    edges: wire
                        .edges
                        .iter()
                        .map(|we| WireEdge {
                            idx: we.idx + e_offset,
                            forward: we.forward,
                        })
                        .collect(),
                })
                .collect();

            shell.faces.push(Face {
                outer_wire: Wire { edges: wire_edges },
                inner_wires,
                normal: {
                    let rotated = mat.transform_vector3(face.normal.into());
                    DVec3::new(rotated.x, rotated.y, rotated.z).normalize_or(face.normal)
                },
                triangles: face.triangles.iter().map(|[i, j, k]| [i + v_offset, j + v_offset, k + v_offset]).collect(),
                mesh_dirty: face.mesh_dirty,
            });
        }
        target.solids.push(Solid {
            shells: vec![shell],
        });
    }

    // Remap face_surface references
    for &fs in &source.geom.face_surface {
        target
            .geom
            .face_surface
            .push(fs.map(|s| s + surface_offset));
    }

    Ok(())
}

fn transform_curve(curve: &Curve3, mat: &DMat4) -> Curve3 {
    let transform_point = |p: DVec3| {
        let r = mat.transform_point3(p.into());
        DVec3::new(r.x, r.y, r.z)
    };
    let transform_direction = |v: DVec3| {
        let r = mat.transform_vector3(v.into());
        DVec3::new(r.x, r.y, r.z).normalize_or(v)
    };

    match curve {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: transform_point(l.origin),
            direction: transform_direction(l.direction),
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3 {
            center: transform_point(c.center),
            normal: transform_direction(c.normal),
            radius: c.radius,
        }),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: transform_point(e.center),
            normal: transform_direction(e.normal),
            major_dir: transform_direction(e.major_dir),
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }),
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: transform_point(h.center),
            normal: transform_direction(h.normal),
            major_dir: transform_direction(h.major_dir),
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }),
        Curve3::BSpline(b) => {
            let mut nb = b.clone();
            for cp in &mut nb.control_points {
                *cp = transform_point(*cp);
            }
            Curve3::BSpline(nb)
        }
        Curve3::Bezier(b) => {
            let mut nb = b.clone();
            for cp in &mut nb.control_points {
                *cp = transform_point(*cp);
            }
            Curve3::Bezier(nb)
        }
        _ => curve.clone(),
    }
}

fn transform_surface(surface: &Surface3, mat: &DMat4) -> Surface3 {
    let transform_point = |p: DVec3| {
        let r = mat.transform_point3(p.into());
        DVec3::new(r.x, r.y, r.z)
    };
    let transform_direction = |v: DVec3| {
        let r = mat.transform_vector3(v.into());
        DVec3::new(r.x, r.y, r.z).normalize_or(v)
    };

    match surface {
        Surface3::Plane(p) => Surface3::Plane(Plane {
            origin: transform_point(p.origin),
            normal: transform_direction(p.normal),
        }),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: transform_point(c.origin),
            axis: transform_direction(c.axis),
            radius: c.radius,
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: transform_point(s.center),
            axis: transform_direction(s.axis),
            radius: s.radius,
        }),
        Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
            apex: transform_point(c.apex),
            axis: transform_direction(c.axis),
            radius: c.radius,
            half_angle_rad: c.half_angle_rad,
        }),
        Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
            center: transform_point(t.center),
            axis: transform_direction(t.axis),
            major_radius: t.major_radius,
            minor_radius: t.minor_radius,
        }),
        Surface3::BSpline(b) => {
            let mut nb = b.clone();
            for row in &mut nb.control_points {
                for cp in row {
                    *cp = transform_point(*cp);
                }
            }
            Surface3::BSpline(nb)
        }
        Surface3::LinearExtrusion(le) => Surface3::LinearExtrusion(LinearExtrusionSurface {
            profile: Box::new(transform_curve(&le.profile, mat)),
            direction: le.direction,
        }),
        Surface3::Revolution(r) => Surface3::Revolution(RevolutionSurface {
            profile: Box::new(transform_curve(&r.profile, mat)),
            axis_origin: transform_point(r.axis_origin),
            axis_dir: transform_direction(r.axis_dir),
        }),
        Surface3::Offset(o) => Surface3::Offset(OffsetSurface {
            basis: Box::new(transform_surface(&o.basis, mat)),
            offset_distance: o.offset_distance,
        }),
        Surface3::Trimmed(t) => Surface3::Trimmed(TrimmedSurface {
            basis: Box::new(transform_surface(&t.basis, mat)),
            trim: t.trim,
        }),
        _ => surface.clone(),
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_modeling::make_box_brep;

    fn make_box() -> BRep {
        let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn linear_pattern_count_1_returns_original() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 1,
            spacing: 2.0,
        };
        let result = linear_pattern(&brep, &params).unwrap();
        let v_result = rcad_kernel::properties::volume(&result);

        assert!((v_result - v_orig).abs() < 1e-9);
    }

    #[test]
    fn linear_pattern_count_3_produces_3x_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: 2.0,
        };
        let result = linear_pattern(&brep, &params).unwrap();
        let v_result = rcad_kernel::properties::volume(&result);

        assert!(
            (v_result - 3.0 * v_orig).abs() < 0.01,
            "expected 3x volume, got {v_result} vs expected {}",
            3.0 * v_orig
        );
    }

    #[test]
    fn linear_pattern_invalid_spacing_returns_error() {
        let brep = make_box();
        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: -1.0,
        };
        assert!(matches!(
            linear_pattern(&brep, &params),
            Err(PatternError::InvalidSpacing)
        ));
    }

    #[test]
    fn linear_pattern_zero_direction_returns_error() {
        let brep = make_box();
        let params = LinearPatternParams {
            direction: DVec3::ZERO,
            count: 3,
            spacing: 1.0,
        };
        assert!(matches!(
            linear_pattern(&brep, &params),
            Err(PatternError::ZeroDirection)
        ));
    }

    #[test]
    fn circular_pattern_count_4_produces_4x_volume() {
        // Use a box offset from the rotation axis so copies don't overlap the origin
        let mut brep = make_box_brep(DVec3::new(3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: std::f64::consts::TAU,
        };
        let result = circular_pattern(&brep, &params).unwrap();

        let total_solids = result.solids.len();
        assert_eq!(total_solids, 4, "expected 4 solids, got {total_solids}");

        let v_result = rcad_kernel::properties::volume(&result);
        assert!(
            (v_result - 4.0 * v_orig).abs() < 0.01,
            "expected 4x volume, got {v_result} vs expected {}",
            4.0 * v_orig
        );
    }

    #[test]
    fn circular_pattern_half_turn_produces_2x_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 2,
            total_angle: std::f64::consts::PI,
        };
        let result = circular_pattern(&brep, &params).unwrap();
        let v_result = rcad_kernel::properties::volume(&result);

        assert!(
            (v_result - 2.0 * v_orig).abs() < 0.01,
            "expected 2x volume, got {v_result} vs expected {}",
            2.0 * v_orig
        );
    }

    #[test]
    fn circular_pattern_invalid_angle_returns_error() {
        let brep = make_box();
        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: 0.0,
        };
        assert!(matches!(
            circular_pattern(&brep, &params),
            Err(PatternError::InvalidAngle)
        ));
    }

    #[test]
    fn circular_pattern_angle_too_large_returns_error() {
        let brep = make_box();
        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: 7.0, // > 2*pi
        };
        assert!(matches!(
            circular_pattern(&brep, &params),
            Err(PatternError::InvalidAngle)
        ));
    }

    #[test]
    fn linear_pattern_zero_count_returns_error() {
        let brep = make_box();
        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 0,
            spacing: 1.0,
        };
        assert!(matches!(
            linear_pattern(&brep, &params),
            Err(PatternError::InvalidCount)
        ));
    }
}
