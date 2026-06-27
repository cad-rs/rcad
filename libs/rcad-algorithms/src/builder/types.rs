use std::collections::HashMap;
use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods::{Orientation, ShapeRef};
use crate::tolerance::*;
use crate::history::{FaceOrigin, ShellOrigin, SolidOrigin};
use crate::bopds::ds::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    Union,
    Intersection,
    Difference,
}

/// OCCT-aligned: TopAbs_ShapeEnum subset used by the Builder pipeline.
/// Matches OCCT's dimension-by-dimension ordering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
}

#[derive(Debug)]
pub enum BooleanError {
    /// Operation type is not valid for the boolean operation.
    InvalidOperation,
    EmptyInput,
    MissingGeometry(&'static str),
    DegenerateResult,
    /// A numeric operation produced a non-finite or NaN value.
    NumericalFailure(&'static str),
    /// An expected non-empty collection was empty (e.g. polyline with no points).
    EmptyCollection(&'static str),
    /// Result fails validity checks (non-manifold, open shells, invalid orientation).
    InvalidResult(&'static str),
    /// Missing intersection curves between surfaces that should intersect.
    IncompleteIntersection(&'static str),
    /// Result contains self-intersecting geometry.
    SelfIntersection(&'static str),
    /// Result shell has edges with incorrect face reference counts (orphan or over-shared).
    OpenShell {
        orphan_edges: Vec<usize>,
        over_shared_edges: Vec<usize>,
    },
}

impl std::fmt::Display for BooleanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidOperation => write!(f, "invalid operation"),
            Self::EmptyInput => write!(f, "empty input"),
            Self::MissingGeometry(msg) => write!(f, "missing geometry: {msg}"),
            Self::DegenerateResult => write!(f, "degenerate result"),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {msg}"),
            Self::EmptyCollection(msg) => write!(f, "unexpected empty collection: {msg}"),
            Self::InvalidResult(msg) => write!(f, "invalid result: {msg}"),
            Self::IncompleteIntersection(msg) => write!(f, "incomplete intersection: {msg}"),
            Self::SelfIntersection(msg) => write!(f, "self-intersection: {msg}"),
            Self::OpenShell { orphan_edges, over_shared_edges } => {
                write!(f, "open shell: {} orphan edges, {} over-shared edges",
                    orphan_edges.len(), over_shared_edges.len())
            }
        }
    }
}

impl std::error::Error for BooleanError {}

/// ✅ OCCT-aligned: classify 闃舵闇€瑕佺殑鏁版嵁,鏇夸唬 FaceSampleData銆?
///    浠?WireFace + WireSegments + DS + face_idx 鎻愬彇銆?
///    sample_point() / surface / normal / boundary 绛?classify 渚濊禆鐨勫瓧娈点€?
#[derive(Debug, Clone)]
pub struct FaceSampleData {
    pub boundary: Vec<DVec3>,
    pub surface: Surface3,
    pub normal: DVec3,
    pub inner_wires: Vec<Vec<DVec3>>,
    pub uv_domain: Option<[f64; 4]>,
    pub uv_centroid: Option<DVec2>,
    pub sample_override: Option<DVec3>,
    pub outer_circle_edges: Vec<(usize, Curve3)>,
    pub seam_edge: Option<(usize, Curve3)>,
    pub inner_wire_circle: Option<(usize, Curve3)>,
}

impl FaceSampleData {
    /// ✅ OCCT-aligned: construct from sub-face data (transitional shim).
    fn from_sub_face(sub: &FaceSampleData) -> Self {
        sub.clone()
    }

    /// Returns a point slightly INSIDE the surface (toward the interior of the solid).
    /// 浠?FaceSampleData::sample_point 绉绘,浣跨敤 WireFace 鐨勬暟鎹簮銆?
    pub(crate) fn sample_point(&self) -> DVec3 {
        if let Some(pt) = self.sample_override {
            return pt;
        }
        match &self.surface {
            Surface3::Sphere(s) => {
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    let sp = s.point_at(uv.x, uv.y);
                    eprintln!("[SAMPLE_PT] sphere uv_centroid=({:.4},{:.4}) 鈫?3D=({:.4},{:.4},{:.4})",
                        uv.x, uv.y, sp.x, sp.y, sp.z);
                    sp
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else {
                    s.center + s.radius * DVec3::X
                };
                let to_center = (s.center - surface_pt).normalize_or_zero();
                let inward = if to_center.length_squared() > 0.5 { to_center } else { -self.normal };
                surface_pt + inward * (TOLERANCE_ABS * 10.0)
            }
            Surface3::Cylinder(c) => {
                use rcad_kernel::geom::SurfaceEval;
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    c.point_at(uv.x, uv.y)
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else {
                    c.origin + c.axis.normalize() * 0.5
                };
                let axis = c.axis.normalize();
                let to_axis = c.origin + axis * (surface_pt - c.origin).dot(axis) - surface_pt;
                let inward = to_axis.normalize_or_zero();
                surface_pt + inward * (TOLERANCE_ABS * 5000.0)
            }
            Surface3::Torus(t) => {
                use rcad_kernel::geom::SurfaceEval;
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    t.point_at(uv.x, uv.y)
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else {
                    t.center + (t.major_radius + t.minor_radius) * DVec3::X
                };
                let axis = t.axis.normalize_or_zero();
                let local = surface_pt - t.center;
                let axial = local.dot(axis);
                let radial = local - axial * axis;
                let inward = if radial.length_squared() > TOLERANCE_FLOAT_ULTRA {
                    let tube_center = t.center + axial * axis + radial.normalize() * t.major_radius;
                    (tube_center - surface_pt).normalize_or_zero()
                } else { -self.normal };
                surface_pt + inward * (TOLERANCE_ABS * 10.0)
            }
            Surface3::Cone(c) => {
                use rcad_kernel::geom::SurfaceEval;
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    c.point_at(uv.x, uv.y)
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else { c.point_at(0.0, 1.0) };
                let axis = c.axis_dir();
                let local = surface_pt - c.apex;
                let axial = local.dot(axis);
                let axis_pt = c.apex + axis * axial;
                let inward = (axis_pt - surface_pt).normalize_or_zero();
                let inward = if inward.length_squared() > 0.5 { inward } else { -self.normal };
                surface_pt + inward * (TOLERANCE_ABS * 5000.0)
            }
            _ => {
                let centroid = if self.boundary.len() >= 3 {
                    planar_polygon_centroid(&self.boundary, self.normal)
                } else if self.boundary.is_empty() { DVec3::ZERO } else {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                };
                centroid + self.normal * TOLERANCE_ABS * 10.0
            }
        }
    }
}

/// DEPRECATED: 鍐呴儴閬楃暀绫诲瀷銆備笉褰卞搷 OCCT 瀵归綈 鈥?浠呭湪 split_face 鍐呴儴 + emit 鍥為€€浣跨敤銆?
/// 澶栭儴鎺ュ彛缁熶竴浣跨敤 FaceSampleData (classify) 鍜?WireFace (emit)銆?
/// OCCT-aligned: wire grouping result — ordered segment chains forming a face boundary.
#[derive(Clone)]
pub struct WireFace {
    pub outer_wire: Vec<usize>,
    pub inner_wires: Vec<Vec<usize>>,
    /// OCCT-aligned: Internal wires from PerformShapesToAvoid (BOPAlgo_BuilderFace.cxx L327-382).
    pub internal_wires: Vec<Vec<usize>>,
}

/// ✅ OCCT-aligned: collected sub-face result before classification.
/// Holds either a wire-pipeline result (to emit via emit_wire_face) or
/// a legacy split_face result (to emit via emit_face_with_origin).
/// Used to defer classification until after all faces are split.
#[derive(Clone)]
pub(crate) enum CollectedFaceResult {
    Wire {
        wf: WireFace,
        segments: Vec<WireSegment>,
        vertex_positions: std::collections::HashMap<usize, DVec3>,
        fi: usize,
        flip: bool,
        origin: FaceOrigin,
    },
    Legacy(FaceSampleData, bool, FaceOrigin),
}

/// OCCT-aligned: Source of a virtual edge segment in the edge-to-wire pipeline.
#[derive(Debug, Clone)]
pub(crate) enum WireEdgeSource {
    DsEdge(usize),
    IntersectionCurve(usize),
    SeamEdge,
}

/// OCCT-aligned: TopAbs_Orientation for WireSegment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum WireOrientation {
    Forward,
    Reversed,
    Internal,
    External,
}

/// OCCT-aligned: Virtual edge used in the edge-to-wire pipeline.
#[derive(Debug, Clone)]
pub(crate) struct WireSegment {
    pub(crate) start_vertex: usize,
    pub(crate) end_vertex: usize,
    pub(crate) source: WireEdgeSource,
    pub(crate) orientation: WireOrientation,
    pub(crate) is_seam: bool,
    pub(crate) tangent_start: Option<f64>,
    pub(crate) tangent_end: Option<f64>,
    /// OCCT DoSplitSEAMOnFace: second pcurve with U shifted by the surface
    /// period (e.g. 2*PI for sphere). Used by refine_angle_2d to project IC
    /// edges onto the other side of the parametric seam, preventing figure-8
    /// wires. Set for seam segments on split-seam periodic surfaces.
    pub(crate) second_pcurve: Option<Curve2d>,
    pub(crate) first_pcurve: Option<Curve2d>,
    /// ✅ OCCT-aligned: vertex parameters on the pcurve (BRep_Tool::Parameter,
    ///   WireSplitter_1.cxx L669). t_range[0] = start_vertex param,
    ///   t_range[1] = end_vertex param.  vertex_uv evaluates pc.point_at(t).
    pub(crate) t_range: [f64; 2],
}

impl WireSegment {
    fn reversed(&self) -> Self {
        WireSegment {
            start_vertex: self.end_vertex,
            end_vertex: self.start_vertex,
            source: match &self.source {
                WireEdgeSource::DsEdge(i) => WireEdgeSource::DsEdge(*i),
                WireEdgeSource::IntersectionCurve(i) => WireEdgeSource::IntersectionCurve(*i),
                WireEdgeSource::SeamEdge => WireEdgeSource::SeamEdge,
            },
            orientation: match self.orientation {
                WireOrientation::Forward => WireOrientation::Reversed,
                WireOrientation::Reversed => WireOrientation::Forward,
                o => o,
            },
            is_seam: self.is_seam,
            second_pcurve: None, first_pcurve: None,
            t_range: [self.t_range[1], self.t_range[0]],
            tangent_start: self.tangent_end
                .map(|a| (a + std::f64::consts::PI) % std::f64::consts::TAU),
            tangent_end: self.tangent_start
                .map(|a| (a + std::f64::consts::PI) % std::f64::consts::TAU),
        }
    }
}

/// Compute the true area centroid of a planar polygon in 3D by projecting onto
/// the plane's 2D orthonormal basis and using the shoelace formula.
/// Guaranteed to lie inside a convex polygon and close to the interior of a
/// concave polygon, unlike the boundary-vertex centroid which can be arbitrarily
/// biased by uneven vertex distribution along the boundary.
fn planar_polygon_centroid(boundary: &[DVec3], normal: DVec3) -> DVec3 {
    if boundary.len() < 3 {
        return if boundary.is_empty() {
            DVec3::ZERO
        } else {
            boundary.iter().copied().sum::<DVec3>() / boundary.len() as f64
        };
    }

    let n = normal.normalize();
    let ref_vec = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let u = n.cross(ref_vec).normalize();
    let v = n.cross(u).normalize();

    let origin = boundary[0];
    let count = boundary.len();

    let mut area2 = 0.0_f64;
    let mut cx6 = 0.0_f64;
    let mut cy6 = 0.0_f64;

    for i in 0..count {
        let j = (i + 1) % count;
        let xi = (boundary[i] - origin).dot(u);
        let yi = (boundary[i] - origin).dot(v);
        let xj = (boundary[j] - origin).dot(u);
        let yj = (boundary[j] - origin).dot(v);
        let cross = xi * yj - xj * yi;
        area2 += cross;
        cx6 += (xi + xj) * cross;
        cy6 += (yi + yj) * cross;
    }

    if area2.abs() < 1e-30 {
        return boundary.iter().copied().sum::<DVec3>() / count as f64;
    }

    let inv = 1.0 / (3.0 * area2);
    origin + u * (cx6 * inv) + v * (cy6 * inv)
}

pub(crate) type FaceEntry = (
    Vec<(usize, bool)>,        // outer wire: (edge_idx, forward)
    Vec<Vec<(usize, bool)>>,   // inner wires: each is Vec<(edge_idx, forward)>
    Vec<[usize; 3]>,
    DVec3,
    Surface3,
    Option<[f64; 4]>,
    DVec3,
    f64,
    DVec3,
    Vec<Vec<(usize, bool)>>,   // internal wire edges (TopAbs_INTERNAL)
);

/// ✅ OCCT-aligned: intermediate result of the LOW-D phase (V+E+W creation)
/// in the dimension-by-dimension pipeline.  Carries the data needed for
/// HIGH-D face assembly from build_face_edges_and_wires to
/// build_face_from_wire_edges, matching OCCT's separation of edge/wire
/// construction from face triangulation/assembly.
pub(crate) struct FaceWireEdges {
    outer_edges: Vec<(usize, bool)>,
    inner_wires_edges: Vec<Vec<(usize, bool)>>,
    internal_wire_edges: Vec<Vec<(usize, bool)>>,
    normal: DVec3,
    surface: Surface3,
    sphere_uv: Option<[f64; 4]>,
    centroid: DVec3,
    area: f64,
    sample_pt: DVec3,
    outer_boundary: Vec<DVec3>,
    iw_boundaries: Vec<Vec<DVec3>>,
    all_vert_indices: Vec<usize>,
    outer_sig: Vec<usize>,
}

/// OCCT-aligned: Source of a virtual edge segment, TopoDS variant.
#[derive(Debug, Clone)]
pub(crate) enum WireEdgeSourceTopoDS {
    DsEdge(ShapeRef),
    IntersectionCurve(ShapeRef),
    SeamEdge,
}

/// OCCT-aligned: Virtual edge using ShapeRef handles instead of usize indices.
/// Designed to carry the same information as WireSegment but with TopoDS handles
/// readable through BRepTool queries.
#[derive(Debug, Clone)]
pub(crate) struct WireSegmentTopoDS {
    pub(crate) edge: ShapeRef,
    pub(crate) face: ShapeRef,
    pub(crate) start_vertex: ShapeRef,
    pub(crate) end_vertex: ShapeRef,
    pub(crate) source: WireEdgeSourceTopoDS,
    pub(crate) orientation: Orientation,
    pub(crate) is_seam: bool,
    pub(crate) tangent_start: Option<f64>,
    pub(crate) tangent_end: Option<f64>,
    pub(crate) first_pcurve: Option<Curve2d>,
    pub(crate) second_pcurve: Option<Curve2d>,
    /// ✅ OCCT-aligned: vertex parameters on the pcurve (BRep_Tool::Parameter).
    ///   t_range[0] = start_vertex param, t_range[1] = end_vertex param.
    pub(crate) t_range: [f64; 2],
}