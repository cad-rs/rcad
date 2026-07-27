use crate::bopds::ds::*;
use crate::history::{FaceOrigin, ShellOrigin, SolidOrigin};
use crate::tolerance::*;
use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods::{Orientation, ShapeRef};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    Union,
    Intersection,
    Difference,
}

/// TopAbs_ShapeEnum subset used by the Builder pipeline.
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
    /// AlertTooFewArguments 闁?fewer than 2 arguments.
    TooFewArguments,
    /// AlertNoFiller 闁?PaveFiller not initialized.
    NoFiller,
    /// AlertBOPNotAllowed 闁?non-licit operation for the arguments.
    BOPNotAllowed,
    /// AlertBOPNotSet 闁?operation type not set.
    BOPNotSet,
    /// AlertEmptyShape 闁?one argument is empty.
    EmptyShape,
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
            Self::TooFewArguments => write!(f, "too few arguments"),
            Self::NoFiller => write!(f, "no pave filler"),
            Self::BOPNotAllowed => write!(f, "operation not allowed for arguments"),
            Self::BOPNotSet => write!(f, "operation type not set"),
            Self::EmptyShape => write!(f, "empty shape argument"),
            Self::EmptyInput => write!(f, "empty input"),
            Self::MissingGeometry(msg) => write!(f, "missing geometry: {msg}"),
            Self::DegenerateResult => write!(f, "degenerate result"),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {msg}"),
            Self::EmptyCollection(msg) => write!(f, "unexpected empty collection: {msg}"),
            Self::InvalidResult(msg) => write!(f, "invalid result: {msg}"),
            Self::IncompleteIntersection(msg) => write!(f, "incomplete intersection: {msg}"),
            Self::SelfIntersection(msg) => write!(f, "self-intersection: {msg}"),
            Self::OpenShell {
                orphan_edges,
                over_shared_edges,
            } => {
                write!(
                    f,
                    "open shell: {} orphan edges, {} over-shared edges",
                    orphan_edges.len(),
                    over_shared_edges.len()
                )
            }
        }
    }
}

impl std::error::Error for BooleanError {}

/// (FaceSampleData removed 閳?OCCT FClass2d + WireFace used instead)

/// DEPRECATED: OCCT = ?  split_face + emit  闁逞屽厴閸? ?
///  闁?=  FaceSampleData (classify)  ?WireFace (emit) ?
/// wire grouping result 闁?ordered segment chains forming a face boundary.
#[derive(Clone)]
pub struct WireFace {
    pub outer_wire: Vec<usize>,
    pub inner_wires: Vec<Vec<usize>>,
    /// Internal wires from PerformShapesToAvoid (BOPAlgo_BuilderFace.cxx L327-382).
    pub internal_wires: Vec<Vec<usize>>,
}

/// 闁collected sub-face result before classification.
/// Holds either a wire-pipeline result (to emit via emit_wire_face) or
/// a legacy split_face result (to emit via emit_face_with_origin).

/// Source of a virtual edge segment in the edge-to-wire pipeline.
#[derive(Debug, Clone)]
pub(crate) enum WireEdgeSource {
    DsEdge(usize),
    IntersectionCurve(usize),
    SeamEdge,
}

/// TopAbs_Orientation for WireSegment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum WireOrientation {
    Forward,
    Reversed,
    Internal,
    External,
}

/// Virtual edge used in the edge-to-wire pipeline.
#[derive(Debug, Clone)]
pub(crate) struct WireSegment {
    pub(crate) start_vertex: usize,
    pub(crate) end_vertex: usize,
    pub(crate) source: WireEdgeSource,
    pub(crate) orientation: WireOrientation,
    pub(crate) is_closed_on_face: bool,
    /// OCCT DoSplitSEAMOnFace: second pcurve with U shifted by the surface
    /// period (e.g. 2*PI for sphere). Used by refine_angle_2d to project IC
    /// edges onto the other side of the parametric seam, preventing figure-8
    /// wires. Set for seam segments on split-seam periodic surfaces.
    pub(crate) second_pcurve: Option<Curve2d>,
    pub(crate) first_pcurve: Option<Curve2d>,
    /// 闁vertex parameters on the pcurve (BRep_Tool::Parameter,
    /// WireSplitter_1.cxx L669). t_range[0] = start_vertex param,
    /// t_range[1] = end_vertex param.  vertex_uv evaluates pc.point_at(t).
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
            is_closed_on_face: self.is_closed_on_face,
            second_pcurve: None,
            first_pcurve: None,
            t_range: [self.t_range[1], self.t_range[0]],
        }
    }
}

/// Compute the true area centroid of a planar polygon in 3D by projecting onto
/// the plane's 2D orthonormal basis and using the shoelace formula.
/// Guaranteed to lie inside a convex polygon and close to the interior of a
/// concave polygon, unlike the boundary-vertex centroid which can be arbitrarily
/// biased by uneven vertex distribution along the boundary.

pub(crate) type FaceEntry = (
    Vec<(usize, bool)>,      // outer wire: (edge_idx, forward)
    Vec<Vec<(usize, bool)>>, // inner wires: each is Vec<(edge_idx, forward)>
    Vec<[usize; 3]>,
    DVec3,
    Surface3,
    Option<[f64; 4]>,
    DVec3,
    f64,
    DVec3,
    Vec<Vec<(usize, bool)>>, // internal wire edges (TopAbs_INTERNAL)
);

/// 闁intermediate result of the LOW-D phase (V+E+W creation)
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

/// Source of a virtual edge segment, TopoDS variant.
#[derive(Debug, Clone)]
pub(crate) enum WireEdgeSourceTopoDS {
    DsEdge(ShapeRef),
    IntersectionCurve(ShapeRef),
    SeamEdge,
}

/// Virtual edge using ShapeRef handles instead of usize indices.
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
    pub(crate) is_closed_on_face: bool,
    pub(crate) first_pcurve: Option<Curve2d>,
    pub(crate) second_pcurve: Option<Curve2d>,
    /// 闁vertex parameters on the pcurve (BRep_Tool::Parameter).
    /// t_range[0] = start_vertex param, t_range[1] = end_vertex param.
    pub(crate) t_range: [f64; 2],
}
