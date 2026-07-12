//! First-stage feature operations (TKFeat-like APIs).
//!
//! This module builds practical feature workflows on top of the existing
//! boolean kernel. The first shipped feature is a cylindrical hole.

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::topods;
use rcad_kernel::geom::{Curve3, Line3, Surface3};

use crate::{BooleanError, BooleanOpType};
use crate::bop_occt_ops::boolean_op_generic as boolean_op;

/// Errors returned by feature operations.
#[derive(Debug)]
pub enum FeatureError {
    NonFiniteInput(&'static str),
    NonPositiveInput(&'static str),
    InvalidInput(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    Modeling(String),
    Boolean(BooleanError),
}

impl std::fmt::Display for FeatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteInput(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveInput(name) => write!(f, "{name} must be > 0"),
            Self::InvalidInput(name) => write!(f, "{name} is invalid"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::Modeling(msg) => write!(f, "modeling operation failed: {msg}"),
            Self::Boolean(err) => write!(f, "boolean operation failed: {err}"),
        }
    }
}

impl std::error::Error for FeatureError {}

impl From<BooleanError> for FeatureError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

impl From<rcad_modeling::BuildError> for FeatureError {
    fn from(value: rcad_modeling::BuildError) -> Self {
        Self::Modeling(value.to_string())
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

fn axis_ref_basis(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), FeatureError> {
    let y_axis = normalize("axis", axis)?;
    let ref_dir = normalize("ref_dir", ref_dir)?;
    let x_reject = ref_dir - y_axis * ref_dir.dot(y_axis);
    if x_reject.length_squared() <= EPS {
        return Err(FeatureError::ParallelVectors("ref_dir", "axis"));
    }
    let x_axis = x_reject.normalize();
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

/// Create a cylindrical through/blind hole by subtracting an oriented cylinder
/// from `target`.
///
/// - `center`: center of the tool cylinder.
/// - `axis`: cylinder axis direction.
/// - `ref_dir`: reference direction used to build local orientation.
/// - `radius`: hole radius.
/// - `depth`: tool cylinder height.
///
/// For through holes, pass a `depth` larger than the part thickness along
/// `axis`.
pub fn make_cylindrical_hole(
    target: &topods::BRep,
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    depth: f64,
) -> Result<topods::BRep, FeatureError> {
    if !center.is_finite() {
        return Err(FeatureError::NonFiniteInput("center"));
    }
    let radius = validate_positive("radius", radius)?;
    let depth = validate_positive("depth", depth)?;

    let (x_axis, y_axis, z_axis) = axis_ref_basis(axis, ref_dir)?;

    let tool_center = center - axis * (depth * 0.5);
    let tool = rcad_modeling::make_cylinder_brep(tool_center, z_axis, x_axis, radius, depth)
        .map_err(|_| FeatureError::InvalidInput("failed to build tool cylinder"))?;
    Ok(boolean_op(BooleanOpType::Difference, target, &tool)?)
}

/// Extrude a closed planar polygon into a solid prism (no boolean against a target).
///
/// Uses the same section loft as [`make_prism`], but returns the tool solid directly.
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

/// Revolve a closed planar polygon (wires only; no boolean against a target).
///
/// Analogous to OCCT DRAW `revol` of an `mkplane` face for solids of revolution.
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

    let (profile, face_idx) = build_polygon_face_brep(profile_verts)?;
    rcad_modeling::revolve(&profile, face_idx, axis_origin, axis_dir, angle_rad).map_err(Into::into)
}

/// Create a prismatic boss or pocket by extruding a polygon profile and
/// performing a boolean union (boss) or difference (pocket) with `target`.
///
/// - `profile_verts`: 3D coplanar polygon vertices in CCW order when viewed
///   along the extrusion direction.  Minimum 3 vertices required.
/// - `direction`: extrusion direction (unit vector is computed internally).
/// - `depth`: extrusion length (must be > 0).
/// - `op`: [`BooleanOpType::Union`] = boss; [`BooleanOpType::Difference`] = pocket.
///
/// Analogous to OCCT `BRepFeat_MakePrism` for linear boss/pocket features.
pub fn make_prism(
    target: &topods::BRep,
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
    op: BooleanOpType,
) -> Result<topods::BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::NonPositiveInput("profile_verts needs >= 3 vertices"));
    }
    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;

    let tool = build_polygon_prism(profile_verts, dir, depth)?;
    Ok(boolean_op(op, target, &tool)?)
}

/// Create a drafted prismatic boss or pocket by extruding a polygon profile
/// with radial taper and applying a boolean operation.
///
/// Positive `draft_angle_rad` expands the top profile outward; negative values
/// shrink it inward.
///
/// Analogous to OCCT `BRepFeat_MakeDPrism` (linear draft prism).
pub fn make_draft_prism(
    target: &topods::BRep,
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
    draft_angle_rad: f64,
    op: BooleanOpType,
) -> Result<topods::BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }
    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;
    let angle = validate_finite("draft_angle_rad", draft_angle_rad)?;
    if angle.abs() >= std::f64::consts::FRAC_PI_2 - TOLERANCE_MESH_LEGACY {
        return Err(FeatureError::InvalidInput("draft_angle_rad must be in (-pi/2, pi/2)"));
    }

    let bot: Vec<DVec3> = profile_verts.to_vec();
    let centroid = bot.iter().copied().fold(DVec3::ZERO, |acc, p| acc + p) / bot.len() as f64;
    let axial = dir * depth;
    let taper = depth * angle.tan();

    let top: Vec<DVec3> = bot
        .iter()
        .map(|&p| {
            let v = p - centroid;
            let radial = v - dir * v.dot(dir);
            let radial_dir = if radial.length_squared() > EPS {
                radial.normalize()
            } else {
                DVec3::ZERO
            };
            p + axial + radial_dir * taper
        })
        .collect();

    let tool = build_prism_from_sections(&bot, &top, dir)?;
    Ok(boolean_op(op, target, &tool)?)
}

/// Create a revolution boss/pocket feature from a planar profile.
///
/// The profile polygon is revolved around `axis_origin + t * axis_dir` by
/// `angle_rad`, then combined with `target` by boolean `op`.
///
/// Analogous to OCCT `BRepFeat_MakeRevol` for linear profile faces.
pub fn make_revolution(
    target: &topods::BRep,
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
    op: BooleanOpType,
) -> Result<topods::BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }
    if !axis_origin.is_finite() {
        return Err(FeatureError::NonFiniteInput("axis_origin"));
    }
    let axis_dir = normalize("axis_dir", axis_dir)?;
    let angle_rad = validate_positive("angle_rad", angle_rad)?;

    let (profile, face_idx) = build_polygon_face_brep(profile_verts)?;
    let tool = rcad_modeling::revolve(&profile, face_idx, axis_origin, axis_dir, angle_rad)?;
    Ok(boolean_op(op, target, &tool)?)
}

fn build_polygon_face_brep(profile_verts: &[DVec3]) -> Result<(topods::BRep, usize), FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }

    let n = profile_verts.len();
    let mut brep = topods::BRep::new();

    // Add vertices
    let verts: Vec<topods::ShapeRef> = profile_verts.iter().map(|&p| brep.add_tvertex(p)).collect();

    // Add edges: closed polygon loop
    let mut edge_refs = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let p0 = profile_verts[i];
        let p1 = profile_verts[j];
        let d = p1 - p0;
        let len = d.length();
        let curve = if len > EPS {
            Some(Curve3::Line(Line3 { origin: p0, direction: d / len }))
        } else {
            None
        };
        edge_refs.push(brep.add_tedge(curve, verts[i], verts[j], [0.0, len]));
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
fn build_polygon_prism(profile_verts: &[DVec3], dir: DVec3, depth: f64) -> Result<topods::BRep, FeatureError> {
    let bot: Vec<DVec3> = profile_verts.to_vec();
    let top: Vec<DVec3> = bot.iter().map(|&p| p + dir * depth).collect();
    build_prism_from_sections(&bot, &top, dir)
}

fn build_prism_from_sections(bot: &[DVec3], top: &[DVec3], dir: DVec3) -> Result<topods::BRep, FeatureError> {
    let n = bot.len();
    if n < 3 || top.len() != n {
        return Err(FeatureError::InvalidInput("section vertex count mismatch"));
    }

    let mut brep = topods::BRep::new();

    // Add vertices: bot[0..n] then top[0..n]
    let bot_sr: Vec<topods::ShapeRef> = bot.iter().map(|&p| brep.add_tvertex(p)).collect();
    let top_sr: Vec<topods::ShapeRef> = top.iter().map(|&p| brep.add_tvertex(p)).collect();

    /// Add a line edge from start to end and return its ShapeRef.
    fn add_line_edge(brep: &mut topods::BRep, p0: DVec3, p1: DVec3, start_sr: topods::ShapeRef, end_sr: topods::ShapeRef) -> topods::ShapeRef {
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > EPS { d / len } else { DVec3::X };
        brep.add_tedge(Some(Curve3::Line(Line3 { origin: p0, direction: dir })), start_sr, end_sr, [0.0, len])
    }

    // Bottom-cap edges: bot[i] -> bot[(i+1)%n]
    let bot_edges: Vec<topods::ShapeRef> = (0..n).map(|i| {
        let j = (i + 1) % n;
        add_line_edge(&mut brep, bot[i], bot[j], bot_sr[i], bot_sr[j])
    }).collect();
    // Top-cap edges: top[i] -> top[(i+1)%n]
    let top_edges: Vec<topods::ShapeRef> = (0..n).map(|i| {
        let j = (i + 1) % n;
        add_line_edge(&mut brep, top[i], top[j], top_sr[i], top_sr[j])
    }).collect();
    // Vertical edges: bot[i] -> top[i]
    let vert_edges: Vec<topods::ShapeRef> = (0..n).map(|i| {
        add_line_edge(&mut brep, bot[i], top[i], bot_sr[i], top_sr[i])
    }).collect();

    // Bottom cap (outward normal = -dir): reverse traversal of bot edges
    let bot_wire: Vec<topods::ShapeRef> = (0..n).map(|i| {
        let ei = bot_edges[n - 1 - i];
        // Reversed orientation
        topods::ShapeRef { ptr_id: 0, index: ei.index, orientation: topods::Orientation::Reversed, location: 0 }
    }).collect();
    let bot_wire_sr = brep.add_twire(bot_wire);
    brep.add_tface(
        Some(Surface3::Plane(rcad_kernel::geom::Plane::new(bot[0], -dir))),
        bot_wire_sr, vec![], None, None, vec![], false,
    );

    // Top cap (outward normal = +dir): forward traversal of top edges
    let top_wire_edges: Vec<topods::ShapeRef> = (0..n).map(|i| top_edges[i]).collect();
    let top_wire_sr = brep.add_twire(top_wire_edges);
    brep.add_tface(
        Some(Surface3::Plane(rcad_kernel::geom::Plane::new(top[0], dir))),
        top_wire_sr, vec![], None, None, vec![], false,
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
            if nv.length_squared() > TOLERANCE_VEC_SQ_MIN { nv.normalize() } else { dir.cross(ab).normalize() }
        };
        // wire: bot[i]->bot[j] (fwd), bot[j]->top[j] (fwd), top[j]->top[i] (rev), top[i]->bot[i] (rev)
        let wire_edges = vec![
            bot_edges[i],
            vert_edges[j],
            topods::ShapeRef { ptr_id: 0, index: top_edges[i].index, orientation: topods::Orientation::Reversed, location: 0 },
            topods::ShapeRef { ptr_id: 0, index: vert_edges[i].index, orientation: topods::Orientation::Reversed, location: 0 },
        ];
        let wire_sr = brep.add_twire(wire_edges);
        brep.add_tface(
            Some(Surface3::Plane(rcad_kernel::geom::Plane::new(a, face_normal))),
            wire_sr, vec![], None, None, vec![], false,
        );
    }

    // Build shell and solid
    let face_refs: Vec<topods::ShapeRef> = brep.tshapes.iter().enumerate()
        .filter(|(_, ts)| matches!(ts.as_ref(), topods::TShape::Face(_)))
        .map(|(i, _)| topods::ShapeRef { ptr_id: 0, index: i, orientation: topods::Orientation::Forward, location: 0 })
        .collect();
    let shell = brep.add_tshell(face_refs);
    brep.add_tsolid(vec![shell]);

    Ok(brep)
}



// ─── SplitShape: split a face by a cutting wire ──────────────────────────────

/// Error returned by [`split_face_by_wire`].
#[derive(Debug)]
pub enum SplitShapeError {
    FaceNotFound,
    CutPathTooShort,
    CutVertexNotOnWire { vertex_idx: usize },
    CutPathClosedLoop,
    DegenerateResult,
}

impl std::fmt::Display for SplitShapeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::FaceNotFound => write!(f, "face index out of range"),
            Self::CutPathTooShort => write!(f, "cut path needs at least 2 vertices"),
            Self::CutVertexNotOnWire { vertex_idx } => {
                write!(f, "cut vertex {vertex_idx} is not on the face outer wire")
            }
            Self::CutPathClosedLoop => write!(f, "cut path start and end are the same wire vertex"),
            Self::DegenerateResult => write!(f, "split produced a degenerate wire"),
        }
    }
}

impl std::error::Error for SplitShapeError {}

/// Split a face by a cutting wire (path of vertex indices already in `brep.vertices`).
///
/// - `cut_path`: at least 2 vertex indices; first and last must appear as the
///   *start* vertex of some edge in the face's outer wire.
/// - New line edges are inserted for each segment of `cut_path`.
/// - The face is replaced by two sub-faces.
///
/// Analogous to OCCT `BRepFeat_SplitShape`.
pub fn split_face_by_wire(
    brep: &mut topods::BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    cut_path: &[usize],
) -> Result<usize, SplitShapeError> {
    use std::sync::Arc;
    if cut_path.len() < 2 {
        return Err(SplitShapeError::CutPathTooShort);
    }

    // Collect all solid TShapes.
    let solid_ts_indices: Vec<usize> = brep.tshapes.iter().enumerate()
        .filter(|(_, ts)| matches!(ts.as_ref(), topods::TShape::Solid(_)))
        .map(|(i, _)| i)
        .collect();
    let solid_tsi = *solid_ts_indices.get(solid_idx).ok_or(SplitShapeError::FaceNotFound)?;
    let solid_sd = match &*brep.tshapes[solid_tsi] {
        topods::TShape::Solid(s) => s.clone(),
        _ => return Err(SplitShapeError::FaceNotFound),
    };

    let shell_sr = solid_sd.shells.get(shell_idx).ok_or(SplitShapeError::FaceNotFound)?;
    let shell_sd = match &*brep.tshapes[shell_sr.index] {
        topods::TShape::Shell(s) => s.clone(),
        _ => return Err(SplitShapeError::FaceNotFound),
    };

    let face_sr = shell_sd.faces.get(face_idx).ok_or(SplitShapeError::FaceNotFound)?;
    let face_sd = match &*brep.tshapes[face_sr.index] {
        topods::TShape::Face(f) => f.clone(),
        _ => return Err(SplitShapeError::FaceNotFound),
    };

    let wire_sr = face_sd.outer_wire;
    let wire_sd = match &*brep.tshapes[wire_sr.index] {
        topods::TShape::Wire(w) => w.clone(),
        _ => return Err(SplitShapeError::FaceNotFound),
    };

    let outer_edges = wire_sd.edges.clone();

    let start_v = cut_path[0];
    let end_v = *cut_path.last().unwrap();

    // Build ordered vertex sequence of the outer wire.
    let wire_verts: Vec<usize> = outer_edges.iter().map(|edge_sr| {
        let ed = match &*brep.tshapes[edge_sr.index] {
            topods::TShape::Edge(e) => e,
            _ => unreachable!(),
        };
        if edge_sr.orientation.is_forward() { ed.first.index } else { ed.last.index }
    }).collect();

    let pos_start = wire_verts.iter().position(|&v| v == start_v)
        .ok_or(SplitShapeError::CutVertexNotOnWire { vertex_idx: start_v })?;
    let pos_end = wire_verts.iter().position(|&v| v == end_v)
        .ok_or(SplitShapeError::CutVertexNotOnWire { vertex_idx: end_v })?;
    if pos_start == pos_end {
        return Err(SplitShapeError::CutPathClosedLoop);
    }

    let n = outer_edges.len();

    /// Build a ShapeRef for vertex by index.
    fn vert_ref(brep: &topods::BRep, idx: usize) -> topods::ShapeRef {
        let ts = &brep.tshapes[idx];
        topods::ShapeRef {
            ptr_id: Arc::as_ptr(ts) as u64,
            index: idx,
            orientation: topods::Orientation::Forward,
            location: 0,
        }
    }

    // Add line edges for the cut path segments.
    let cut_edge_refs: Vec<topods::ShapeRef> = cut_path.windows(2).map(|w| {
        let (sv, ev) = (w[0], w[1]);
        let p0 = brep.vertex_point(sv).unwrap();
        let p1 = brep.vertex_point(ev).unwrap();
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > TOLERANCE_LEN_SQ_DIV_SAFE { d / len } else { DVec3::X };
        brep.add_tedge(
            Some(Curve3::Line(Line3 { origin: p0, direction: dir })),
            vert_ref(brep, sv), vert_ref(brep, ev), [0.0, len],
        )
    }).collect();

    // Half A: outer[pos_start..pos_end] + cut forward.
    let mut half_a: Vec<topods::ShapeRef> = (0..(pos_end - pos_start))
        .map(|i| {
            let ei = &outer_edges[(pos_start + i) % n];
            *ei
        })
        .collect();
    half_a.extend(cut_edge_refs.iter().copied());
    let wire_a = brep.add_twire(half_a);

    // Half B: outer[pos_end..] + outer[..pos_start] + cut reversed.
    let half_b_len = n - (pos_end - pos_start);
    let mut half_b: Vec<topods::ShapeRef> = (0..half_b_len)
        .map(|i| outer_edges[(pos_end + i) % n])
        .collect();
    half_b.extend(cut_edge_refs.iter().rev().map(|e| topods::ShapeRef {
        ptr_id: e.ptr_id,
        index: e.index,
        orientation: topods::Orientation::Reversed,
        location: e.location,
    }));
    let wire_b = brep.add_twire(half_b);

    if brep.tshapes.get(wire_a.index).and_then(|ts| {
        if let topods::TShape::Wire(w) = &**ts { Some(w.edges.len()) } else { None }
    }).unwrap_or(0) < 3
        || brep.tshapes.get(wire_b.index).and_then(|ts| {
            if let topods::TShape::Wire(w) = &**ts { Some(w.edges.len()) } else { None }
        }).unwrap_or(0) < 3
    {
        return Err(SplitShapeError::DegenerateResult);
    }

    // Create the two sub-faces
    let orig_surface = face_sd.surface.clone();
    let orig_inner = face_sd.inner_wires.clone();

    let face_a_ref = brep.add_tface(orig_surface.clone(), wire_a, orig_inner.clone(), face_sd.sample_point, face_sd.uv_domain, face_sd.internal_vertices.clone(), face_sd.natural_restriction);
    let face_b_ref = brep.add_tface(orig_surface, wire_b, orig_inner, face_sd.sample_point, face_sd.uv_domain, face_sd.internal_vertices.clone(), face_sd.natural_restriction);

    // Add faces to the shell
    let sd = brep.shell_mut(*shell_sr);
    sd.faces[face_idx] = face_a_ref;
    sd.faces.insert(face_idx + 1, face_b_ref);
    sd.my_shapes.push(face_b_ref);

    Ok(1)
}

// ─── Linear rib / slot ───────────────────────────────────────────────────────

/// Create a linear rib (or slot) feature via prism boolean.
///
/// Analogous to OCCT `BRepFeat_MakeLinearForm`.
pub fn make_linear_rib(
    target: &topods::BRep,
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
    op: BooleanOpType,
) -> Result<topods::BRep, FeatureError> {
    make_prism(target, profile_verts, direction, depth, op)
}

/// Create a revolution rib/slot feature via revolve boolean.
///
/// Analogous to OCCT `BRepFeat_MakeRevolutionForm`.
pub fn make_revolution_rib(
    target: &topods::BRep,
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
    op: BooleanOpType,
) -> Result<topods::BRep, FeatureError> {
    make_revolution(target, profile_verts, axis_origin, axis_dir, angle_rad, op)
}


