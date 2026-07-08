//! First-stage feature operations (TKFeat-like APIs).
//!
//! This module builds practical feature workflows on top of the existing
//! boolean kernel. The first shipped feature is a cylindrical hole.

use crate::tolerance::*;
use glam::{DAffine3, DMat3, DVec3};
use rcad_kernel::{topods, BRep, GeomStore, PrimitiveSolid};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::{BooleanError, BooleanOpType, boolean_op};

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
    let target_old = rcad_kernel::BRep::from_topods(target);
    let tool_old = rcad_kernel::BRep::from_topods(&tool);
    Ok(boolean_op(BooleanOpType::Difference, &target_old, &tool_old)?)
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

    let profile = build_polygon_face_brep(profile_verts)?;
    rcad_modeling::revolve(&profile, 0, axis_origin, axis_dir, angle_rad).map_err(Into::into)
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
    let target_old = rcad_kernel::BRep::from_topods(target);
    let tool_old = rcad_kernel::BRep::from_topods(&tool);
    Ok(boolean_op(op, &target_old, &tool_old)?)
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
    let target_old = rcad_kernel::BRep::from_topods(target);
    let tool_old = rcad_kernel::BRep::from_topods(&tool);
    Ok(boolean_op(op, &target_old, &tool_old)?)
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

    let profile = build_polygon_face_brep(profile_verts)?;
    let tool = rcad_modeling::revolve(&profile, 0, axis_origin, axis_dir, angle_rad)?;
    let target_old = rcad_kernel::BRep::from_topods(target);
    let tool_old = rcad_kernel::BRep::from_topods(&tool);
    Ok(boolean_op(op, &target_old, &tool_old)?)
}

fn build_polygon_face_brep(profile_verts: &[DVec3]) -> Result<topods::BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }

    let n = profile_verts.len();
    let mut brep = BRep {
        vertices: Vec::with_capacity(n),
        edges: Vec::with_capacity(n),
        solids: Vec::new(),
        geom: GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    for &p in profile_verts {
        brep.vertices.push(Vertex { point: p });
    }

    for i in 0..n {
        let j = (i + 1) % n;
        brep.edges.push(Edge { start: i, end: j });
    }

    // Single shell with one face containing all edges as a single outer wire.
    let wire = Wire {
        edges: (0..n).map(|i| WireEdge { idx: i, forward: true }).collect(),
    };
    let face = Face {
        outer_wire: wire,
        inner_wires: Vec::new(),
        triangles: Vec::new(),
        normal: DVec3::Z,
        sample_point: None,
        mesh_dirty: true,
        surface_idx: None,
    };
    brep.solids.push(Solid {
        shells: vec![Shell { faces: vec![face] }],
    });

    Ok(brep.to_topods())
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

    let mut brep = BRep {
        vertices: Vec::with_capacity(2 * n),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    // Add vertices: bot[0..n] then top[0..n]
    // bot vertex index: i; top vertex index: n + i
    for &p in bot { brep.vertices.push(Vertex { point: p }); }
    for &p in top { brep.vertices.push(Vertex { point: p }); }

    /// Add a line edge from start to end and return its index.
    fn add_line_edge(brep: &mut BRep, start: usize, end: usize) -> usize {
        let p0 = brep.vertices[start].point;
        let p1 = brep.vertices[end].point;
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > 0.0 { d / len } else { DVec3::X };
        let ei = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, len]));
        brep.geom.edge_degenerated.push(false);
        ei
    }

    // Bottom-cap edges: bot[i] -> bot[(i+1)%n]
    let bot_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, i, (i + 1) % n)).collect();
    // Top-cap edges: top[i] -> top[(i+1)%n]
    let top_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, n + i, n + (i + 1) % n)).collect();
    // Vertical edges: bot[i] -> top[i]
    let vert_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, i, n + i)).collect();

    let mut faces = Vec::with_capacity(n + 2);

    // Bottom cap (outward normal = -dir): reverse traversal of bot edges
    {
        let wire_edges: Vec<WireEdge> = (0..n)
            .map(|i| WireEdge { idx: bot_edges[n - 1 - i], forward: false })
            .collect();
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: -dir, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: bot[0], normal: -dir }));
        brep.geom.face_surface.push(Some(si));
    }

    // Top cap (outward normal = +dir): forward traversal of top edges
    {
        let wire_edges: Vec<WireEdge> = (0..n)
            .map(|i| WireEdge { idx: top_edges[i], forward: true })
            .collect();
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: dir, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: top[0], normal: dir }));
        brep.geom.face_surface.push(Some(si));
    }

    // Lateral quad faces: quad bot[i] -> bot[j] -> top[j] -> top[i] for each edge i
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
            WireEdge { idx: bot_edges[i],  forward: true },
            WireEdge { idx: vert_edges[j], forward: true },
            WireEdge { idx: top_edges[i],  forward: false },
            WireEdge { idx: vert_edges[i], forward: false },
        ];
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: face_normal, triangles: vec![], sample_point: None, mesh_dirty: true, surface_idx: None });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: a, normal: face_normal }));
        brep.geom.face_surface.push(Some(si));
    }

    brep.solids.push(Solid { shells: vec![Shell { faces }] });
    Ok(brep.to_topods())
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
    brep: &mut BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    cut_path: &[usize],
) -> Result<usize, SplitShapeError> {
    if cut_path.len() < 2 {
        return Err(SplitShapeError::CutPathTooShort);
    }
    if solid_idx >= brep.solids.len()
        || shell_idx >= brep.solids[solid_idx].shells.len()
        || face_idx >= brep.solids[solid_idx].shells[shell_idx].faces.len()
    {
        return Err(SplitShapeError::FaceNotFound);
    }

    let start_v = cut_path[0];
    let end_v = *cut_path.last().unwrap();

    let outer_edges = brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
        .outer_wire.edges.clone();

    // Build ordered vertex sequence of the outer wire.
    let wire_verts: Vec<usize> = outer_edges.iter().map(|we| {
        let e = &brep.edges[we.idx];
        if we.forward { e.start } else { e.end }
    }).collect();

    let pos_start = wire_verts.iter().position(|&v| v == start_v)
        .ok_or(SplitShapeError::CutVertexNotOnWire { vertex_idx: start_v })?;
    let pos_end = wire_verts.iter().position(|&v| v == end_v)
        .ok_or(SplitShapeError::CutVertexNotOnWire { vertex_idx: end_v })?;
    if pos_start == pos_end {
        return Err(SplitShapeError::CutPathClosedLoop);
    }

    let n = outer_edges.len();

    // Add line edges for the cut path segments.
    let cut_edge_indices: Vec<usize> = cut_path.windows(2).map(|w| {
        let (sv, ev) = (w[0], w[1]);
        let ei = brep.edges.len();
        let p0 = brep.vertices[sv].point;
        let p1 = brep.vertices[ev].point;
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > TOLERANCE_LEN_SQ_DIV_SAFE { d / len } else { DVec3::X };
        brep.edges.push(Edge { start: sv, end: ev });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, len]));
        brep.geom.edge_degenerated.push(false);
        ei
    }).collect();

    // Half A: outer[pos_start..pos_end] + cut forward.
    let half_a: Vec<WireEdge> = (0..(pos_end - pos_start))
        .map(|i| outer_edges[(pos_start + i) % n]).collect();
    let cut_fwd: Vec<WireEdge> = cut_edge_indices.iter().map(|&ei| WireEdge::fwd(ei)).collect();
    let mut wire_a = half_a;
    wire_a.extend_from_slice(&cut_fwd);

    // Half B: outer[pos_end..] + outer[..pos_start] + cut reversed.
    let half_b_len = n - (pos_end - pos_start);
    let half_b: Vec<WireEdge> = (0..half_b_len)
        .map(|i| outer_edges[(pos_end + i) % n]).collect();
    let cut_rev: Vec<WireEdge> = cut_edge_indices.iter().rev().map(|&ei| WireEdge::rev(ei)).collect();
    let mut wire_b = half_b;
    wire_b.extend_from_slice(&cut_rev);

    if wire_a.len() < 3 || wire_b.len() < 3 {
        return Err(SplitShapeError::DegenerateResult);
    }

    let orig_normal = brep.solids[solid_idx].shells[shell_idx].faces[face_idx].normal;
    let orig_inner = brep.solids[solid_idx].shells[shell_idx].faces[face_idx].inner_wires.clone();

    let face_a = Face {
        outer_wire: Wire { edges: wire_a },
        inner_wires: orig_inner.clone(),
        normal: orig_normal,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,

    };
    let face_b = Face {
        outer_wire: Wire { edges: wire_b },
        inner_wires: orig_inner,
        normal: orig_normal,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,

    };

    // Update GeomStore face_surface flat index.
    let flat_idx: usize = brep.solids[..solid_idx].iter()
        .flat_map(|s| s.shells.iter()).map(|sh| sh.faces.len()).sum::<usize>()
        + brep.solids[solid_idx].shells[..shell_idx].iter()
            .map(|sh| sh.faces.len()).sum::<usize>()
        + face_idx;
    let orig_surf = brep.geom.face_surface.get(flat_idx).copied().flatten();
    if flat_idx < brep.geom.face_surface.len() {
        brep.geom.face_surface.insert(flat_idx + 1, orig_surf);
    }
    if flat_idx < brep.geom.face_tolerance.len() {
        let ft = brep.geom.face_tolerance[flat_idx];
        brep.geom.face_tolerance.insert(flat_idx + 1, ft);
    }

    brep.solids[solid_idx].shells[shell_idx].faces[face_idx] = face_a;
    brep.solids[solid_idx].shells[shell_idx].faces.insert(face_idx + 1, face_b);
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


