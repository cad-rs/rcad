//! User-facing geometry construction helpers.
//!
//! This module is the public modeling entry layer for RCAD.
//! The API intentionally prefers OCCT-style direct constructor functions
//! over fluent builder structs.

pub mod brep_builder;
mod curve;
pub mod fillet;
pub mod ops;
mod solid;
mod surface;
pub mod wire_ops;

pub use curve::*;
pub use fillet::{CornerBlendHistory, FilletHistory, MultiFilletHistory, SafeFilletResult};
pub use fillet::{
    chamfer_edge, chamfer_edge_angle, chamfer_edge_safe, corner_blend, fillet_edge,
    fillet_edge_safe, fillet_edge_variable_radius, fillet_edges,
};
pub use fillet::{
    chamfer_edge_angle_topods, chamfer_edge_topods, fillet_edge_topods,
    fillet_edge_variable_radius_topods,
};
pub use fillet::{
    chamfer_edge_angle_with_history, chamfer_edge_with_history, corner_blend_with_history,
    fillet_edge_variable_radius_with_history, fillet_edge_with_history, fillet_edges_with_history,
};
pub use ops::*;
pub use solid::*;
pub use surface::*;
pub use wire_ops::{chamfer_wire_2d, fillet_wire_2d, project_wire_onto_surface};

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Surface3};
use rcad_kernel::topods::{Orientation, Shape, TShape};
use std::error::Error;
use std::fmt;
use std::sync::Arc;

const EPS: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    NonFiniteValue(&'static str),
    NonPositiveValue(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    DegenerateGeometry(&'static str),
    InvalidIndex(usize),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteValue(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveValue(name) => write!(f, "{name} must be > 0"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::DegenerateGeometry(msg) => write!(f, "degenerate geometry: {msg}"),
            Self::InvalidIndex(idx) => write!(f, "invalid index: {idx}"),
        }
    }
}

impl Error for BuildError {}

fn validate_point(name: &'static str, point: DVec3) -> Result<DVec3, BuildError> {
    if point.is_finite() {
        Ok(point)
    } else {
        Err(BuildError::NonFiniteValue(name))
    }
}

fn validate_positive(name: &'static str, value: f64) -> Result<f64, BuildError> {
    if !value.is_finite() {
        Err(BuildError::NonFiniteValue(name))
    } else if value <= 0.0 {
        Err(BuildError::NonPositiveValue(name))
    } else {
        Ok(value)
    }
}

fn normalize_vector(name: &'static str, vector: DVec3) -> Result<DVec3, BuildError> {
    validate_point(name, vector)?;
    if vector.length_squared() <= EPS {
        Err(BuildError::ZeroVector(name))
    } else {
        Ok(vector.normalize())
    }
}

fn normalize_rejection(
    name: &'static str,
    vector: DVec3,
    reference_name: &'static str,
    reference: DVec3,
) -> Result<DVec3, BuildError> {
    let vector = normalize_vector(name, vector)?;
    let rejected = vector - reference * vector.dot(reference);
    if rejected.length_squared() <= EPS {
        Err(BuildError::ParallelVectors(name, reference_name))
    } else {
        Ok(rejected.normalize())
    }
}

fn basis_from_x_y(x_dir: DVec3, y_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BuildError> {
    let x_axis = normalize_vector("x_dir", x_dir)?;
    let y_axis = normalize_rejection("y_dir", y_dir, "x_dir", x_axis)?;
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

fn basis_from_axis_ref(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BuildError> {
    let y_axis = normalize_vector("axis", axis)?;
    let x_axis = normalize_rejection("ref_dir", ref_dir, "axis", y_axis)?;
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

fn translate_brep(brep: &mut BRep, offset: DVec3) {
    for ts in &mut brep.tshapes {
        match &mut *Arc::make_mut(ts) {
            TShape::Vertex(v) => v.point += offset,
            TShape::Edge(e) => {
                if let Some(ref mut curve) = e.curve {
                    match curve {
                        Curve3::Line(l) => l.origin += offset,
                        Curve3::Circle(c) => c.center += offset,
                        Curve3::Ellipse(el) => el.center += offset,
                        Curve3::Hyperbola(h) => h.center += offset,
                        Curve3::BSpline(b) => {
                            for cp in &mut b.control_points {
                                *cp += offset;
                            }
                        }
                        Curve3::Bezier(b) => {
                            for cp in &mut b.control_points {
                                *cp += offset;
                            }
                        }
                        _ => {}
                    }
                }
            }
            TShape::Face(f) => {
                if let Some(ref mut surface) = f.surface {
                    match surface {
                        Surface3::Plane(p) => p.origin += offset,
                        Surface3::Cylinder(c) => c.origin += offset,
                        Surface3::Sphere(s) => s.center += offset,
                        Surface3::Cone(c) => c.apex += offset,
                        Surface3::Torus(t) => t.center += offset,
                        Surface3::BSpline(b) => {
                            for row in &mut b.control_points {
                                for cp in row {
                                    *cp += offset;
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }
}

fn do_mirror_brep(brep: &BRep, plane_origin: DVec3, plane_normal: DVec3) -> BRep {
    let n = plane_normal.normalize();
    let mirror_point = |p: DVec3| -> DVec3 {
        let d = (p - plane_origin).dot(n);
        p - n * (2.0 * d)
    };
    let mirror_vec = |v: DVec3| -> DVec3 { v - n * (2.0 * v.dot(n)) };

    let mut out = BRep::new();

    // Build old-index → new-Shape mapping
    let n_shapes = brep.tshapes.len();
    let mut old_to_new: Vec<Option<Shape>> = vec![None; n_shapes];

    // Pass 1: mirror vertices (no dependencies)
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Vertex(v) = &**ts {
            old_to_new[i] = Some(out.add_tvertex(mirror_point(v.point)));
        }
    }

    // Pass 2: mirror edges (depends on vertices)
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Edge(e) = &**ts {
            let first = old_to_new[e.first.index].clone().unwrap_or(e.first.clone());
            let last = old_to_new[e.last.index].clone().unwrap_or(e.last.clone());
            let curve = e.curve.as_ref().map(|c| {
                use rcad_kernel::geom::{Circle3, Ellipse3, Hyperbola3, Line3};
                match c {
                    Curve3::Line(l) => Curve3::Line(Line3 {
                        origin: mirror_point(l.origin),
                        direction: mirror_vec(l.direction),
                    }),
                    Curve3::Circle(c) => Curve3::Circle(Circle3::new(
                        mirror_point(c.center),
                        mirror_vec(c.normal),
                        c.radius,
                    )),
                    Curve3::Ellipse(ec) => Curve3::Ellipse(Ellipse3 {
                        center: mirror_point(ec.center),
                        normal: mirror_vec(ec.normal),
                        major_dir: mirror_vec(ec.major_dir),
                        major_radius: ec.major_radius,
                        minor_radius: ec.minor_radius,
                    }),
                    Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
                        center: mirror_point(h.center),
                        normal: mirror_vec(h.normal),
                        major_dir: mirror_vec(h.major_dir),
                        semi_major: h.semi_major,
                        semi_minor: h.semi_minor,
                    }),
                    Curve3::BSpline(b) => {
                        let mut nb = b.clone();
                        for cp in &mut nb.control_points {
                            *cp = mirror_point(*cp);
                        }
                        Curve3::BSpline(nb)
                    }
                    Curve3::Bezier(b) => {
                        let mut nb = b.clone();
                        for cp in &mut nb.control_points {
                            *cp = mirror_point(*cp);
                        }
                        Curve3::Bezier(nb)
                    }
                    _ => c.clone(),
                }
            });
            old_to_new[i] = Some(out.add_tedge(curve, first, last, e.range));
        }
    }

    // Pass 3: mirror wires (depends on edges, flip edge orientation)
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Wire(w) = &**ts {
            let edges: Vec<Shape> = w
                .edges
                .iter()
                .map(|e_sr| Shape {
                    data: e_sr.data.clone(),
                    index: e_sr.index,
                    orientation: match e_sr.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        other => other,
                    },
                    location: e_sr.location,
                })
                .collect();
            old_to_new[i] = Some(out.add_twire(edges));
        }
    }

    // Pass 4: mirror faces (depends on wires, mirror surface)
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Face(f) = &**ts {
            let surface = f.surface.as_ref().map(|s| {
                use rcad_kernel::geom::{
                    ConicalSurface, CylindricalSurface, Plane, SphericalSurface, ToroidalSurface,
                };
                match s {
                    Surface3::Plane(p) => {
                        Surface3::Plane(Plane::new(mirror_point(p.origin), mirror_vec(p.normal)))
                    }
                    Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
                        origin: mirror_point(c.origin),
                        axis: mirror_vec(c.axis),
                        radius: c.radius,
                        ref_dir: mirror_vec(c.ref_dir),
                    }),
                    Surface3::Sphere(sp) => Surface3::Sphere(SphericalSurface {
                        center: mirror_point(sp.center),
                        axis: mirror_vec(sp.axis),
                        radius: sp.radius,
                        ref_dir: mirror_vec(sp.ref_dir),
                    }),
                    Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
                        apex: mirror_point(c.apex),
                        axis: mirror_vec(c.axis),
                        radius: c.radius,
                        half_angle_rad: c.half_angle_rad,
                    }),
                    Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
                        center: mirror_point(t.center),
                        axis: mirror_vec(t.axis),
                        major_radius: t.major_radius,
                        minor_radius: t.minor_radius,
                    }),
                    Surface3::BSpline(b) => {
                        let mut nb = b.clone();
                        for row in &mut nb.control_points {
                            for cp in row {
                                *cp = mirror_point(*cp);
                            }
                        }
                        Surface3::BSpline(nb)
                    }
                    _ => s.clone(),
                }
            });
            let outer_wire = old_to_new[f.outer_wire.index].clone().unwrap_or(f.outer_wire.clone());
            let inner_wires: Vec<Shape> = f
                .inner_wires
                .iter()
                .map(|w_sr| old_to_new[w_sr.index].clone().unwrap_or_else(|| w_sr.clone()))
                .collect();
            let sample_point = f.sample_point.map(mirror_point);
            old_to_new[i] = Some(out.add_tface(
                surface,
                outer_wire,
                inner_wires,
                sample_point,
                f.uv_domain,
                Vec::new(),
                f.natural_restriction,
            ));
        }
    }

    // Pass 5: mirror shells (depends on faces)
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Shell(s) = &**ts {
            let faces: Vec<Shape> = s
                .faces
                .iter()
                .map(|f_sr| old_to_new[f_sr.index].clone().unwrap_or_else(|| f_sr.clone()))
                .collect();
            old_to_new[i] = Some(out.add_tshell(faces));
        }
    }

    // Pass 6: mirror solids (depends on shells)
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Solid(s) = &**ts {
            let shells: Vec<Shape> = s
                .shells
                .iter()
                .map(|sh_sr| old_to_new[sh_sr.index].clone().unwrap_or_else(|| sh_sr.clone()))
                .collect();
            old_to_new[i] = Some(out.add_tsolid(shells));
        }
    }

    out
}

fn transform_brep(brep: &mut BRep, origin: DVec3, x_axis: DVec3, y_axis: DVec3, z_axis: DVec3) {
    let xform_point = |p: DVec3| -> DVec3 { origin + x_axis * p.x + y_axis * p.y + z_axis * p.z };
    let xform_vec =
        |v: DVec3| -> DVec3 { (x_axis * v.x + y_axis * v.y + z_axis * v.z).normalize_or_zero() };

    for ts in &mut brep.tshapes {
        match &mut *Arc::make_mut(ts) {
            TShape::Vertex(v) => v.point = xform_point(v.point),
            TShape::Edge(e) => {
                if let Some(ref mut curve) = e.curve {
                    match curve {
                        Curve3::Line(l) => {
                            l.origin = xform_point(l.origin);
                            l.direction = xform_vec(l.direction);
                        }
                        Curve3::Circle(c) => {
                            c.center = xform_point(c.center);
                            c.normal = xform_vec(c.normal);
                        }
                        Curve3::Ellipse(el) => {
                            el.center = xform_point(el.center);
                            el.normal = xform_vec(el.normal);
                            el.major_dir = xform_vec(el.major_dir);
                        }
                        Curve3::Hyperbola(h) => {
                            h.center = xform_point(h.center);
                            h.normal = xform_vec(h.normal);
                            h.major_dir = xform_vec(h.major_dir);
                        }
                        Curve3::BSpline(b) => {
                            for cp in &mut b.control_points {
                                *cp = xform_point(*cp);
                            }
                        }
                        Curve3::Bezier(b) => {
                            for cp in &mut b.control_points {
                                *cp = xform_point(*cp);
                            }
                        }
                        _ => {}
                    }
                }
            }
            TShape::Face(f) => {
                if let Some(ref mut surface) = f.surface {
                    match surface {
                        Surface3::Plane(p) => {
                            p.origin = xform_point(p.origin);
                            p.normal = xform_vec(p.normal);
                        }
                        Surface3::Cylinder(c) => {
                            c.origin = xform_point(c.origin);
                            c.axis = xform_vec(c.axis);
                            c.ref_dir = xform_vec(c.ref_dir);
                        }
                        Surface3::Sphere(s) => {
                            s.center = xform_point(s.center);
                            s.axis = xform_vec(s.axis);
                        }
                        Surface3::Cone(c) => {
                            c.apex = xform_point(c.apex);
                            c.axis = xform_vec(c.axis);
                        }
                        Surface3::Torus(t) => {
                            t.center = xform_point(t.center);
                            t.axis = xform_vec(t.axis);
                        }
                        Surface3::BSpline(b) => {
                            for row in &mut b.control_points {
                                for cp in row {
                                    *cp = xform_point(*cp);
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{PrimitiveSolid, Surface3, topods::TShape};

    #[test]
    fn line_rejects_zero_direction() {
        let err = line(DVec3::ZERO, DVec3::ZERO).unwrap_err();
        assert_eq!(err, BuildError::ZeroVector("direction"));
    }

    #[test]
    fn ellipse_rejects_parallel_major_direction() {
        let err = ellipse(DVec3::ZERO, DVec3::Z, DVec3::Z, 2.0, 1.0).unwrap_err();
        assert_eq!(err, BuildError::ParallelVectors("major_dir", "normal"));
    }

    #[test]
    fn box_brep_builds_transformed_vertices() {
        let brep = box_brep(DVec3::new(1.0, 2.0, 3.0), DVec3::Y, DVec3::Z, 2.0, 3.0, 4.0).unwrap();

        let pts: Vec<DVec3> = brep
            .tshapes
            .iter()
            .filter_map(|ts| match ts.as_ref() {
                TShape::Vertex(v) => Some(v.point),
                _ => None,
            })
            .collect();
        assert_eq!(pts.len(), 8);
        assert!(pts.contains(&DVec3::new(1.0, 2.0, 3.0)));
        assert!(pts.contains(&DVec3::new(5.0, 4.0, 6.0)));
    }

    #[test]
    fn sphere_brep_translates_bounds() {
        let brep = sphere_brep(DVec3::new(10.0, -2.0, 4.0), 2.0).unwrap();

        let y_vals: Vec<f64> = brep
            .tshapes
            .iter()
            .filter_map(|ts| match ts.as_ref() {
                TShape::Vertex(v) => Some(v.point.y),
                _ => None,
            })
            .collect();
        let min_y = y_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_y = y_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        assert!((min_y - (-4.0)).abs() < 1e-6);
        assert!((max_y - 0.0).abs() < 1e-6);
    }

    #[test]
    fn cylinder_primitive_returns_expected_shape() {
        let primitive = cylinder_primitive(3.0, 5.0).unwrap();

        match primitive {
            PrimitiveSolid::Cylinder { radius, height } => {
                assert_eq!(radius, 3.0);
                assert_eq!(height, 5.0);
            }
            other => panic!("expected cylinder primitive, got {other:?}"),
        }
    }

    #[test]
    fn make_plane_alias_matches_plane_constructor() {
        let surface = make_plane(DVec3::new(1.0, 2.0, 3.0), DVec3::Z).unwrap();

        match surface {
            Surface3::Plane(plane) => {
                assert_eq!(plane.origin, DVec3::new(1.0, 2.0, 3.0));
                assert_eq!(plane.normal, DVec3::Z);
            }
            other => panic!("expected plane surface, got {other:?}"),
        }
    }

    #[test]
    fn mirror_box_across_xy_plane() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let mirrored = do_mirror_brep(&brep, DVec3::ZERO, DVec3::Z);
        let v_mirrored = rcad_kernel::properties::volume(&mirrored);

        assert!(
            (v_mirrored - v_orig).abs() < 0.01,
            "mirror should preserve volume: {v_orig} vs {v_mirrored}"
        );

        let orig_pts: Vec<DVec3> = brep
            .tshapes
            .iter()
            .filter_map(|ts| match ts.as_ref() {
                TShape::Vertex(v) => Some(v.point),
                _ => None,
            })
            .collect();
        let mir_pts: Vec<DVec3> = mirrored
            .tshapes
            .iter()
            .filter_map(|ts| match ts.as_ref() {
                TShape::Vertex(v) => Some(v.point),
                _ => None,
            })
            .collect();
        for (i, pt) in orig_pts.iter().enumerate() {
            let mp = mir_pts[i];
            assert!(
                (mp.z - (-pt.z)).abs() < 1e-9,
                "vertex {i} z should be negated"
            );
        }
    }
}
