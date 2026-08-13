//! Revolve a closed planar polygon around an axis into a solid (topods).
//!
//! Replaces the legacy `rcad_modeling::revolve` (removed with the old builder
//! API). Lateral faces are exact `RevolutionSurface`s; the swept boundary arcs
//! are `Circle3` edges; the start/end profile caps are planar faces. The
//! polygon must lie in a plane containing the axis (the OCCT `revol` usage).

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Line3, Plane, RevolutionSurface, Surface3};
use rcad_kernel::topods::{self, Orientation};

/// Rotate `p` around the axis `(origin, dir)` by `angle` radians.
fn rotate_point(p: DVec3, origin: DVec3, dir: DVec3, angle: f64) -> DVec3 {
    let v = p - origin;
    let v_para = dir * v.dot(dir);
    let v_perp = v - v_para;
    let (s, c) = angle.sin_cos();
    v_para + v_perp * c + dir.cross(v_perp) * s + origin
}

/// Revolve a closed planar polygon (coplanar with the axis) around an axis.
///
/// `angle_rad` in `(0, 2π]`. For `2π` the start/end profile coincide and the
/// result is a closed solid of revolution without caps.
pub fn revolve_polygon(
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
) -> Result<topods::BRep, String> {
    let n = profile_verts.len();
    if n < 3 {
        return Err("profile_verts needs >= 3 vertices".into());
    }
    if !axis_origin.is_finite() || !axis_dir.is_finite() {
        return Err("axis must be finite".into());
    }
    let dir = axis_dir.normalize_or_zero();
    if dir.length_squared() < 1e-24 {
        return Err("axis_dir must be non-zero".into());
    }
    if !(angle_rad > 0.0 && angle_rad <= std::f64::consts::TAU + 1e-9) {
        return Err("angle_rad must be in (0, 2*PI]".into());
    }
    let angle = angle_rad.min(std::f64::consts::TAU);

    let mut brep = topods::BRep::new();

    // Project each vertex onto the axis to get the arc center / radius.
    let centers: Vec<DVec3> = profile_verts
        .iter()
        .map(|&p| axis_origin + dir * (p - axis_origin).dot(dir))
        .collect();
    let radii: Vec<f64> = profile_verts
        .iter()
        .zip(&centers)
        .map(|(&p, &o)| (p - o).length())
        .collect();

    // Swept positions of each vertex (end angle).
    let end_verts: Vec<DVec3> = profile_verts
        .iter()
        .map(|&p| rotate_point(p, axis_origin, dir, angle))
        .collect();

    // Vertices: start profile (0..n) then end profile (n..2n). When angle ==
    // 2π the end positions coincide with the start positions; share vertices.
    let full_turn = (angle - std::f64::consts::TAU).abs() < 1e-9;
    let start_srs: Vec<topods::Shape> = profile_verts
        .iter()
        .map(|&p| brep.add_tvertex(p))
        .collect();
    let end_srs: Vec<topods::Shape> = if full_turn {
        start_srs.clone()
    } else {
        end_verts.iter().map(|&p| brep.add_tvertex(p)).collect()
    };

    // Arc edges: vertex i sweeps from start to end.
    let mut arc_edges: Vec<topods::Shape> = Vec::with_capacity(n);
    for i in 0..n {
        let r = radii[i];
        if r < 1e-12 {
            // Degenerate (on-axis) vertex: no arc edge; reuse a degenerate
            // edge via the start vertex so wires stay closed.
            let sr = brep.add_tedge(None, start_srs[i].clone(), end_srs[i].clone(), [0.0, 0.0]);
            arc_edges.push(sr);
            continue;
        }
        let center = centers[i];
        let x_dir = (profile_verts[i] - center) / r;
        let y_dir = dir.cross(x_dir).normalize();
        let curve = Curve3::Circle(Circle3 {
            center,
            normal: dir,
            x_dir,
            y_dir,
            radius: r,
        });
        let sr = brep.add_tedge(Some(curve), start_srs[i].clone(), end_srs[i].clone(), [0.0, angle]);
        arc_edges.push(sr);
    }

    // Profile edge helper: line from p0 to p1 at the given vertex refs.
    let add_profile_edge = |brep: &mut topods::BRep,
                            p0: DVec3,
                            p1: DVec3,
                            v0: &topods::Shape,
                            v1: &topods::Shape|
     -> topods::Shape {
        let d = p1 - p0;
        let len = d.length();
        let curve = if len > 1e-12 {
            Some(Curve3::Line(Line3 {
                origin: p0,
                direction: d / len,
            }))
        } else {
            None
        };
        brep.add_tedge(curve, v0.clone(), v1.clone(), [0.0, len])
    };

    // Profile edges at start and end positions.
    let mut start_edges: Vec<topods::Shape> = Vec::with_capacity(n);
    let mut end_edges: Vec<topods::Shape> = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        start_edges.push(add_profile_edge(
            &mut brep,
            profile_verts[i],
            profile_verts[j],
            &start_srs[i],
            &start_srs[j],
        ));
        if !full_turn {
            end_edges.push(add_profile_edge(
                &mut brep,
                end_verts[i],
                end_verts[j],
                &end_srs[i],
                &end_srs[j],
            ));
        }
    }

    // Lateral faces: one per profile edge.
    let mut face_srs = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let profile_line = Curve3::Line(Line3 {
            origin: profile_verts[i],
            direction: (profile_verts[j] - profile_verts[i]).normalize_or_zero(),
        });
        let surface = Surface3::Revolution(RevolutionSurface {
            profile: Box::new(profile_line),
            axis_origin,
            axis_dir: dir,
        });
        // Wire: start_edge_i (fwd), arc_j (fwd), end_edge_i (rev), arc_i (rev).
        let rev = |sr: &topods::Shape| topods::Shape {
            data: sr.data.clone(),
            index: sr.index,
            orientation: Orientation::Reversed,
            location: sr.location,
        };
        let wire_edges = if full_turn {
            vec![
                start_edges[i].clone(),
                arc_edges[j].clone(),
                rev(&start_edges[i]),
                rev(&arc_edges[i]),
            ]
        } else {
            vec![
                start_edges[i].clone(),
                arc_edges[j].clone(),
                rev(&end_edges[i]),
                rev(&arc_edges[i]),
            ]
        };
        let wire = brep.add_twire(wire_edges);
        let face = brep.add_tface(Some(surface), wire, vec![], None, None, vec![], false);
        face_srs.push(face);
    }

    // Caps (only when not a full turn).
    if !full_turn {
        // Bottom cap: start profile polygon.
        let normal = {
            let a = profile_verts[1] - profile_verts[0];
            let b = profile_verts[2] - profile_verts[0];
            let n = a.cross(b);
            if n.length_squared() > 1e-24 {
                n.normalize()
            } else {
                dir
            }
        };
        let bottom_wire = brep.add_twire(start_edges.clone());
        let bottom_face = brep.add_tface(
            Some(Surface3::Plane(Plane::new(profile_verts[0], -normal))),
            bottom_wire,
            vec![],
            None,
            None,
            vec![],
            false,
        );
        face_srs.push(bottom_face);

        let end_normal = rotate_point(normal, axis_origin, dir, angle);
        let top_wire = brep.add_twire(end_edges.clone());
        let top_face = brep.add_tface(
            Some(Surface3::Plane(Plane::new(end_verts[0], end_normal))),
            top_wire,
            vec![],
            None,
            None,
            vec![],
            false,
        );
        face_srs.push(top_face);
    }

    let shell = brep.add_tshell(face_srs);
    brep.add_tsolid(vec![shell]);
    Ok(brep)
}
