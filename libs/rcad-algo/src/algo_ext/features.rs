//! Polygon feature helpers (extrude / revolve) for the OCCT grid tests.
//!
//! Migrated from the legacy `rcad-algorithms::features` module; only the
//! functions the generated DRAW tests use. `revolve` uses the local
//! [`super::revolve::revolve_polygon`] (the old `rcad_modeling::revolve` was
//! removed with the legacy builder API).

use crate::algo_ext::revolve::revolve_polygon;
use crate::algo_ext::tolerance::{TOLERANCE_LEN_MIN, TOLERANCE_VEC_SQ_MIN};
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::topods;

/// Errors returned by feature operations.
#[derive(Debug)]
pub enum FeatureError {
    NonFiniteInput(&'static str),
    NonPositiveInput(&'static str),
    InvalidInput(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    Modeling(String),
}

impl std::fmt::Display for FeatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteInput(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveInput(name) => write!(f, "{name} must be > 0"),
            Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} and {b} must not be parallel"),
            Self::Modeling(msg) => write!(f, "modeling error: {msg}"),
        }
    }
}

impl std::error::Error for FeatureError {}

impl From<String> for FeatureError {
    fn from(value: String) -> Self {
        Self::Modeling(value)
    }
}

const EPS: f64 = TOLERANCE_LEN_MIN;

fn validate_finite(name: &'static str, v: f64) -> Result<f64, FeatureError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(FeatureError::NonFiniteInput(name))
    }
}

fn validate_positive(name: &'static str, v: f64) -> Result<f64, FeatureError> {
    let v = validate_finite(name, v)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(FeatureError::NonPositiveInput(name))
    }
}

fn normalize(name: &'static str, v: DVec3) -> Result<DVec3, FeatureError> {
    if !v.is_finite() {
        return Err(FeatureError::NonFiniteInput(name));
    }
    if v.length_squared() <= EPS {
        return Err(FeatureError::ZeroVector(name));
    }
    Ok(v.normalize())
}

/// Extrude a closed planar polygon into a solid prism.
pub fn extrude_polygon_solid(
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
) -> Result<topods::BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput(
            "profile_verts needs >= 3 vertices",
        ));
    }
    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;
    build_polygon_prism(profile_verts, dir, depth)
}

/// Revolve a closed planar polygon around an axis into a solid of revolution.
pub fn revolve_polygon_solid(
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
) -> Result<topods::BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput(
            "profile_verts needs >= 3 vertices",
        ));
    }
    if !axis_origin.is_finite() {
        return Err(FeatureError::NonFiniteInput("axis_origin"));
    }
    let axis_dir = normalize("axis_dir", axis_dir)?;
    let angle_rad = validate_positive("angle_rad", angle_rad)?;

    revolve_polygon(profile_verts, axis_origin, axis_dir, angle_rad).map_err(Into::into)
}

fn build_polygon_face_brep(profile_verts: &[DVec3]) -> Result<(topods::BRep, usize), FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput(
            "profile_verts needs >= 3 vertices",
        ));
    }

    let n = profile_verts.len();
    let mut brep = topods::BRep::new();

    // Add vertices
    let verts: Vec<topods::Shape> = profile_verts.iter().map(|&p| brep.add_tvertex(p)).collect();

    // Add edges: closed polygon loop
    let mut edge_refs = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let p0 = profile_verts[i];
        let p1 = profile_verts[j];
        let d = p1 - p0;
        let len = d.length();
        let curve = if len > EPS {
            Some(Curve3::Line(Line3 {
                origin: p0,
                direction: d / len,
            }))
        } else {
            None
        };
        edge_refs.push(brep.add_tedge(curve, verts[i].clone(), verts[j].clone(), [0.0, len]));
    }

    // Wire from all edges
    let wire = brep.add_twire(edge_refs);

    // Face with the wire
    let face = brep.add_tface(None, wire, vec![], None, None, vec![], false);
    let face_idx = face.index;

    // Shell and solid
    let shell = brep.add_tshell(vec![face]);
    brep.add_tsolid(vec![shell]);

    Ok((brep, face_idx))
}

/// Build a solid BRep prism from a polygon profile (n vertices, coplanar) extruded
fn build_polygon_prism(
    profile_verts: &[DVec3],
    dir: DVec3,
    depth: f64,
) -> Result<topods::BRep, FeatureError> {
    let bot: Vec<DVec3> = profile_verts.to_vec();
    let top: Vec<DVec3> = bot.iter().map(|&p| p + dir * depth).collect();
    build_prism_from_sections(&bot, &top, dir)
}

fn build_prism_from_sections(
    bot: &[DVec3],
    top: &[DVec3],
    dir: DVec3,
) -> Result<topods::BRep, FeatureError> {
    let n = bot.len();
    if n < 3 || top.len() != n {
        return Err(FeatureError::InvalidInput("section vertex count mismatch"));
    }

    let mut brep = topods::BRep::new();

    // Add vertices: bot[0..n] then top[0..n]
    let bot_sr: Vec<topods::Shape> = bot.iter().map(|&p| brep.add_tvertex(p)).collect();
    let top_sr: Vec<topods::Shape> = top.iter().map(|&p| brep.add_tvertex(p)).collect();

    /// Add a line edge from start to end and return its Shape.
    fn add_line_edge(
        brep: &mut topods::BRep,
        p0: DVec3,
        p1: DVec3,
        start_sr: topods::Shape,
        end_sr: topods::Shape,
    ) -> topods::Shape {
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > EPS { d / len } else { DVec3::X };
        brep.add_tedge(
            Some(Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            })),
            start_sr,
            end_sr,
            [0.0, len],
        )
    }

    // Bottom-cap edges: bot[i] -> bot[(i+1)%n]
    let bot_edges: Vec<topods::Shape> = (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            add_line_edge(&mut brep, bot[i].clone(), bot[j].clone(), bot_sr[i].clone(), bot_sr[j].clone())
        })
        .collect();
    // Top-cap edges: top[i] -> top[(i+1)%n]
    let top_edges: Vec<topods::Shape> = (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            add_line_edge(&mut brep, top[i].clone(), top[j].clone(), top_sr[i].clone(), top_sr[j].clone())
        })
        .collect();
    // Vertical edges: bot[i] -> top[i]
    let vert_edges: Vec<topods::Shape> = (0..n)
        .map(|i| add_line_edge(&mut brep, bot[i].clone(), top[i].clone(), bot_sr[i].clone(), top_sr[i].clone()))
        .collect();

    // Bottom cap (outward normal = -dir): reverse traversal of bot edges
    let bot_wire: Vec<topods::Shape> = (0..n)
        .map(|i| {
            let ei = bot_edges[n - 1 - i].clone();
            // Reversed orientation
            topods::Shape {
                data: ei.data.clone(),
                index: ei.index,
                orientation: topods::Orientation::Reversed,
                location: 0,
            }
        })
        .collect();
    let bot_wire_sr = brep.add_twire(bot_wire);
    brep.add_tface(
        Some(Surface3::Plane(rcad_kernel::geom::Plane::new(bot[0], -dir))),
        bot_wire_sr,
        vec![],
        None,
        None,
        vec![],
        false,
    );

    // Top cap (outward normal = +dir): forward traversal of top edges
    let top_wire_edges: Vec<topods::Shape> = (0..n).map(|i| top_edges[i].clone()).collect();
    let top_wire_sr = brep.add_twire(top_wire_edges);
    brep.add_tface(
        Some(Surface3::Plane(rcad_kernel::geom::Plane::new(top[0], dir))),
        top_wire_sr,
        vec![],
        None,
        None,
        vec![],
        false,
    );

    // Lateral quad faces
    for i in 0..n {
        let j = (i + 1) % n;
        let a = bot[i];
        let b = bot[j];
        let c = top[j];
        let face_normal = {
            let ab = b - a;
            let ac = c - a;
            let nv = ac.cross(ab);
            if nv.length_squared() > TOLERANCE_VEC_SQ_MIN {
                nv.normalize()
            } else {
                dir.cross(ab).normalize()
            }
        };
        // wire: bot[i]->bot[j] (fwd), bot[j]->top[j] (fwd), top[j]->top[i] (rev), top[i]->bot[i] (rev)
        let wire_edges = vec![
            bot_edges[i].clone(),
            vert_edges[j].clone(),
            topods::Shape {
                data: top_edges[i].data.clone(),
                index: top_edges[i].index,
                orientation: topods::Orientation::Reversed,
                location: 0,
            },
            topods::Shape {
                data: vert_edges[i].data.clone(),
                index: vert_edges[i].index,
                orientation: topods::Orientation::Reversed,
                location: 0,
            },
        ];
        let wire_sr = brep.add_twire(wire_edges);
        brep.add_tface(
            Some(Surface3::Plane(rcad_kernel::geom::Plane::new(
                a,
                face_normal,
            ))),
            wire_sr,
            vec![],
            None,
            None,
            vec![],
            false,
        );
    }

    // Build shell and solid
    let face_refs: Vec<topods::Shape> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter(|(_, ts)| matches!(ts.as_ref(), topods::TShape::Face(_)))
        .map(|(i, ts)| topods::Shape {
            data: ts.clone(),
            index: i,
            orientation: topods::Orientation::Forward,
            location: 0,
        })
        .collect();
    let shell = brep.add_tshell(face_refs);
    brep.add_tsolid(vec![shell]);

    Ok(brep)
}
