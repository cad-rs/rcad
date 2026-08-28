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
        // OCCT BRepLib_MakeEdge forces V1=FWD, V2=REV (BRepLib_MakeEdge.cxx
        // L772-774) — the WireSplitter smart map keys the in/out flags by the
        // composed vertex orientation, so the last vertex must carry REVERSED.
        let end_sr = topods::Shape {
            orientation: topods::Orientation::Reversed,
            ..verts[j].clone()
        };
        edge_refs.push(brep.add_tedge(curve, verts[i].clone(), end_sr, [0.0, len]));
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

    // OCCT BRepPrimAPI_MakePrism: the swept top is a located copy of the
    // source face — same TShapes + a Location (BRepSweep_Trsf::Process). The
    // prism therefore contributes n vertex TShapes, not 2n. Separate top
    // TShapes inflate the nbshapes VERTEX/EDGE counts of OCCT boolean tests
    // (bfuse_simple/D9: 18V/30E vs reference 12V/28E).
    let trans = top[0] - bot[0];
    let loc = brep.add_location(glam::DAffine3::from_translation(trans));
    // Add vertices: bot[0..n]; the top reuses the same TShapes + Location.
    let bot_sr: Vec<topods::Shape> = bot.iter().map(|&p| brep.add_tvertex(p)).collect();
    let top_sr: Vec<topods::Shape> = bot_sr
        .iter()
        .map(|s| topods::Shape {
            location: loc,
            ..s.clone()
        })
        .collect();

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
        // OCCT BRepLib_MakeEdge forces V1=FWD, V2=REV (BRepLib_MakeEdge.cxx
        // L772-774) — the WireSplitter smart map keys the in/out flags by the
        // composed vertex orientation, so the last vertex must carry REVERSED.
        let end_sr = topods::Shape {
            orientation: topods::Orientation::Reversed,
            ..end_sr
        };
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
    // Top-cap edges: located copies of the bottom edges.
    let top_edges: Vec<topods::Shape> = bot_edges
        .iter()
        .map(|e| topods::Shape {
            location: loc,
            ..e.clone()
        })
        .collect();
    // Vertical edges: bot[i] -> top[i] (the top vertex is the located copy).
    let vert_edges: Vec<topods::Shape> = (0..n)
        .map(|i| add_line_edge(&mut brep, bot[i].clone(), top[i].clone(), bot_sr[i].clone(), top_sr[i].clone()))
        .collect();

    // Caps. OCCT BRepSweep_Translation: the bottom cap IS the source face —
    // natural plane = the profile winding normal n_src, wire = profile edges
    // FORWARD, shell reference FORWARD. The top cap is the located copy with
    // the SAME natural plane and a REVERSED shell reference. The whole shell
    // carries DirectSolid's orientation: REVERSED when sweep . n_src > 0
    // (BRepSweep_Translation.cxx DirectSolid: myVec.DotCross(du, dv) > 0).
    let mut n_src = DVec3::ZERO;
    for i in 0..n {
        let d0 = bot[(i + 1) % n] - bot[i];
        let d1 = bot[(i + 2) % n] - bot[(i + 1) % n];
        let c = d0.cross(d1);
        if c.length_squared() > TOLERANCE_VEC_SQ_MIN {
            n_src = c.normalize();
            break;
        }
    }
    let bot_wire_sr = brep.add_twire(bot_edges.clone());
    let bot_face = brep.add_tface(
        Some(Surface3::Plane(rcad_kernel::geom::Plane::new(bot[0], n_src))),
        bot_wire_sr,
        vec![],
        None,
        None,
        vec![],
        false,
    );

    // Top cap: the located copy of the source face (same TShape + Location,
    // OCCT BRepSweep_Trsf::Process, BRepSweep_Trsf.cxx L83-90: newShape.Move(loc)
    // — the wire is NOT reversed).
    let top_wire_edges: Vec<topods::Shape> = (0..n).map(|i| top_edges[i].clone()).collect();
    let top_wire_sr = brep.add_twire(top_wire_edges);
    let top_face = brep.add_tface(
        Some(Surface3::Plane(rcad_kernel::geom::Plane::new(top[0], n_src))),
        top_wire_sr,
        vec![],
        None,
        None,
        vec![],
        false,
    );

    // Lateral quad faces
    let mut lateral_faces: Vec<topods::Shape> = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = bot[i];
        let b = bot[j];
        // OCCT BRepSweep_Translation::MakeEmptyFace plane (via
        // GeomAdaptor_SurfaceOfLinearExtrusion::Plane): u = the profile edge
        // direction, v = -sweep direction (the base edge at v=0, the extruded
        // copy at v=-depth), origin at the profile edge start. Runner
        // [DBG-WIRE] bfuse_simple D9: p2 f1 o=(-100,0,100) u=(0,1,0)
        // v=(0,0,-1).
        let u_dir = (b - a).normalize_or_zero();
        let v_dir = -dir;
        let n_nat = u_dir.cross(v_dir).normalize_or_zero();
        // OCCT BRepSweep_NumLinearRegularSweep lateral wire (the intermediate
        // newWire of BuildShell): [sweep@V1 F, sweep@V2 R, profile R, top F].
        // The traversal is CCW in the natural frame (runner [DBG-WIRE] p2 f1 /
        // p1 f0). The previous [bot F, vertV2 F, top R, vertV1 R] winding was
        // the bitwise inverse, so the FClass2d loop verdict of every split
        // lateral face inverted to is_hole and BuilderFace dropped the face
        // image (bfuse_simple D9: p2's x=-100 / y=100 sides vanished from the
        // draft solid, a_sfs 14 -> 12, empty result).
        let wire_edges = vec![
            vert_edges[i].clone(),
            topods::Shape {
                orientation: topods::Orientation::Reversed,
                ..vert_edges[j].clone()
            },
            topods::Shape {
                orientation: topods::Orientation::Reversed,
                ..bot_edges[i].clone()
            },
            top_edges[i].clone(),
        ];
        let wire_sr = brep.add_twire(wire_edges);
        let lat_face = brep.add_tface(
            Some(Surface3::Plane(rcad_kernel::geom::Plane {
                origin: a,
                normal: n_nat,
                u_dir,
                v_dir,
            })),
            wire_sr,
            vec![],
            None,
            None,
            vec![],
            false,
        );
        lateral_faces.push(lat_face);
    }

    // Build shell and solid. OCCT BRepSweep_NumLinearRegularSweep::BuildShell:
    // the source (bottom) face and the lateral faces are referenced FORWARD,
    // the translated (top) face REVERSED; the whole shell carries DirectSolid's
    // orientation (REVERSED when sweep . n_src > 0).
    let mut face_refs: Vec<topods::Shape> = Vec::with_capacity(n + 2);
    face_refs.push(bot_face.clone());
    face_refs.extend(lateral_faces.iter().cloned());
    face_refs.push(topods::Shape {
        orientation: topods::Orientation::Reversed,
        ..top_face.clone()
    });
    let mut shell = brep.add_tshell(face_refs);
    if dir.dot(n_src) > 0.0 {
        shell.orientation = topods::Orientation::Reversed;
    }
    brep.add_tsolid(vec![shell]);

    Ok(brep)
}
