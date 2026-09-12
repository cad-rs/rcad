//! OCCT BRepCheck_Analyzer (TKTopAlgo/BRepCheck) — the strict 1:1 port.
//!
//! Sources:
//! - `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Analyzer.cxx` (L44-539)
//! - `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Analyzer.hxx` (L37-166)
//! - the BRepCheck_Result / BRepCheck_{Vertex,Edge,Wire,Face,Shell,Solid} ports
//!   in the sibling files of this directory.
//!
//! The sibling modules are declared here with `#[path]` so the existing
//! `brep_check/mod.rs` keeps its single `pub mod brep_check_analyzer;` line.
//!
//! # GAP list (sub-checks whose geometric machinery is not yet translated;
//! they are neutral — no status set, matching the OCCT pass-success path —
//! and never fabricate results)
//!
//! - `Geom2dInt_GInter` curve-curve 2D intersection (general curve input):
//!   needed by `BRepCheck_Wire::SelfIntersect` (Wire.cxx L1074-1742) and the
//!   static `BRepCheck_Face::Intersect` (Face.cxx L626-798). The rcad
//!   `GInter` covers the line-vs-curve case only; the intersection runs are
//!   neutral, the surrounding structure is ported.
//! - `BRepCheck_Edge::InContext` geometric control when the reference curve
//!   is a pcurve (Edge.cxx L411-453 with myHCurve = Adaptor3d_CurveOnSurface):
//!   the rcad `BRepLibValidateEdge` translation accepts only GeomAdaptorCurve
//!   references — neutral.
//! - `BRepCheck_Vertex` point-representation identity matching
//!   (Vertex.cxx L166-240, L260-274): the rcad point representations carry
//!   pool indices without a handle-identity mapping — the pr loops are
//!   neutral; the edge First/Last parameter control (L181-202) is ported.
//!
//! # Architecture deviations (documented at their sites)
//!
//! - OCCT `TopoDS_Shape` carries the Location in the handle; rcad stores a
//!   location index resolved through the `BRep` table, so the Result methods
//!   take `&BRep`.
//! - OCCT `BRep_TEdge::Curves()` is flattened from `TEdgeData`
//!   (`curve` + `representations` + the `pcurves` fast index) by
//!   `edge_curve_reps` (see brep_check_result.rs).
//! - The multi-threaded `BRepCheck_ParallelAnalyzer` functor and the
//!   `OSD_ThreadPool` task split of `Perform` collapse into the sequential
//!   per-shape pass (the OCCT default `theIsParallel = false` runs the same
//!   functor sequentially through `OSD_Parallel::For(..., !myIsParallel)`).
//! - OCCT mutex locking is a no-op for the single-threaded analyzer and is
//!   not modelled.
//! - OCCT `try/catch` failure paths (`SetFailStatus`) only trigger on C++
//!   exceptions; the rcad None-data branches follow the in-try null-handle
//!   paths of the OCCT code instead.

#[path = "brep_check_result.rs"]
pub mod brep_check_result;
#[path = "brep_check_vertex.rs"]
pub mod brep_check_vertex;
#[path = "brep_check_edge.rs"]
pub mod brep_check_edge;
#[path = "brep_check_wire.rs"]
pub mod brep_check_wire;
#[path = "brep_check_wire_self_intersect.rs"]
pub mod brep_check_wire_self_intersect;
#[path = "brep_check_face.rs"]
pub mod brep_check_face;
#[path = "brep_check_shell.rs"]
pub mod brep_check_shell;
#[path = "brep_check_solid.rs"]
pub mod brep_check_solid;

pub use brep_check_edge::{BRepCheckEdge, HCurveAdaptor};
pub use brep_check_face::BRepCheckFace;
pub use brep_check_result::{BRepCheckResultBase, BRepCheckStatus, StatusMap};
pub use brep_check_shell::BRepCheckShell;
pub use brep_check_solid::BRepCheckSolid;
pub use brep_check_vertex::BRepCheckVertex;
pub use brep_check_wire::BRepCheckWire;

use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, ShapeType, TShape};

use brep_check_edge::BRepCheckEdge as EdgeRes;
use brep_check_face::BRepCheckFace as FaceRes;
use brep_check_shell::BRepCheckShell as ShellRes;
use brep_check_solid::BRepCheckSolid as SolidRes;
use brep_check_vertex::BRepCheckVertex as VertexRes;
use brep_check_wire::BRepCheckWire as WireRes;

/// OCCT `occ::handle<BRepCheck_Result>` — the typed variant of the per-shape
/// result (`occ::down_cast` maps to the enum match).
#[derive(Debug)]
pub enum AnyResult {
    /// The null handle (COMPSOLID / COMPOUND — OCCT myMap stores HR = NULL).
    Null,
    Vertex(Box<VertexRes>),
    Edge(Box<EdgeRes>),
    Wire(Box<WireRes>),
    Face(Box<FaceRes>),
    Shell(Box<ShellRes>),
    Solid(Box<SolidRes>),
}

impl AnyResult {
    fn base(&self) -> Option<&BRepCheckResultBase> {
        match self {
            AnyResult::Null => None,
            AnyResult::Vertex(r) => Some(&r.base),
            AnyResult::Edge(r) => Some(&r.base),
            AnyResult::Wire(r) => Some(&r.base),
            AnyResult::Face(r) => Some(&r.base),
            AnyResult::Shell(r) => Some(&r.base),
            AnyResult::Solid(r) => Some(&r.base),
        }
    }

    fn base_mut(&mut self) -> Option<&mut BRepCheckResultBase> {
        match self {
            AnyResult::Null => None,
            AnyResult::Vertex(r) => Some(&mut r.base),
            AnyResult::Edge(r) => Some(&mut r.base),
            AnyResult::Wire(r) => Some(&mut r.base),
            AnyResult::Face(r) => Some(&mut r.base),
            AnyResult::Shell(r) => Some(&mut r.base),
            AnyResult::Solid(r) => Some(&mut r.base),
        }
    }

    /// OCCT BRepCheck_Result::InContext — the virtual dispatch.
    fn in_context(&mut self, brep: &BRep, s: &Shape) {
        match self {
            AnyResult::Null => {}
            AnyResult::Vertex(r) => r.in_context(brep, s),
            AnyResult::Edge(r) => r.in_context(brep, s),
            AnyResult::Wire(r) => r.in_context(brep, s),
            AnyResult::Face(r) => r.in_context(brep, s),
            AnyResult::Shell(r) => r.in_context(brep, s),
            AnyResult::Solid(r) => r.in_context(brep, s),
        }
    }

    /// OCCT BRepCheck_Result::SetFailStatus.
    #[allow(dead_code)]
    fn set_fail_status(&mut self, s: &Shape) {
        if let Some(base) = self.base_mut() {
            base.set_fail_status(s);
        }
    }
}

/// OCCT BRepCheck_Analyzer (Analyzer.hxx L37-166) — a framework to check the
/// overall validity of a shape.
pub struct BRepCheckAnalyzer {
    /// OCCT myShape.
    my_shape: Shape,
    /// OCCT myMap (NCollection_IndexedDataMap<TopoDS_Shape,
    /// Handle(BRepCheck_Result)>) — insertion order preserved.
    my_map: Vec<(Shape, AnyResult)>,
    my_map_index: std::collections::HashMap<brep_check_result::ShapeKey, usize>,
    /// OCCT myIsParallel (the analyzer runs the sequential default).
    my_is_parallel: bool,
    /// OCCT myIsExact.
    my_is_exact: bool,
}

impl BRepCheckAnalyzer {
    /// OCCT BRepCheck_Analyzer::BRepCheck_Analyzer(S, GeomControls = true,
    /// theIsParallel = false, theIsExact = false) (Analyzer.hxx L59-67).
    pub fn new(brep: &BRep, s: &Shape, geom_controls: bool) -> Self {
        let mut a = BRepCheckAnalyzer {
            my_shape: Shape::null(),
            my_map: Vec::new(),
            my_map_index: std::collections::HashMap::new(),
            my_is_parallel: false,
            my_is_exact: false,
        };
        a.init(brep, s, geom_controls);
        a
    }

    /// OCCT BRepCheck_Analyzer::SetExactMethod (Analyzer.hxx L92).
    pub fn set_exact_method(&mut self, the_is_exact: bool) {
        self.my_is_exact = the_is_exact;
    }

    /// OCCT BRepCheck_Analyzer::IsExactMethod (Analyzer.hxx L95).
    pub fn is_exact_method(&self) -> bool {
        self.my_is_exact
    }

    /// OCCT BRepCheck_Analyzer::SetParallel (Analyzer.hxx L98).
    pub fn set_parallel(&mut self, the_is_parallel: bool) {
        self.my_is_parallel = the_is_parallel;
    }

    /// OCCT BRepCheck_Analyzer::IsParallel (Analyzer.hxx L101).
    pub fn is_parallel(&self) -> bool {
        self.my_is_parallel
    }

    /// OCCT BRepCheck_Analyzer::Result(theSubS) (Analyzer.hxx L149-152).
    pub fn result(&self, the_sub_s: &Shape) -> Option<&AnyResult> {
        self.my_map_index
            .get(&brep_check_result::ShapeKey::of(the_sub_s))
            .map(|&i| &self.my_map[i].1)
    }

    /// OCCT BRepCheck_Analyzer::Init(S, B) (Analyzer.cxx L352-363).
    pub fn init(&mut self, brep: &BRep, the_shape: &Shape, b: bool) {
        if the_shape.is_null() {
            panic!("BRepCheck_Analyzer::Init() - NULL shape");
        }

        self.my_shape = the_shape.clone();
        self.my_map.clear();
        self.my_map_index.clear();
        self.put(brep, the_shape, b);
        self.perform(brep);
    }

    /// OCCT BRepCheck_Analyzer::Put(S, B) (Analyzer.cxx L367-416).
    fn put(&mut self, brep: &BRep, the_shape: &Shape, b: bool) {
        if self.my_map_index.contains_key(&brep_check_result::ShapeKey::of(the_shape)) {
            return;
        }

        // OCCT L374-404: HR by shape type.
        let hr: Option<AnyResult> = match the_shape.shape_type() {
            ShapeType::Vertex => Some(AnyResult::Vertex(Box::new(VertexRes::new(brep, the_shape)))),
            ShapeType::Edge => {
                let mut r = EdgeRes::new(brep, the_shape);
                // OCCT L382-383.
                r.set_geometric_controls(b);
                r.set_exact_method(self.my_is_exact);
                Some(AnyResult::Edge(Box::new(r)))
            }
            ShapeType::Wire => {
                let mut r = WireRes::new(brep, the_shape);
                // OCCT L387.
                r.set_geometric_controls(b);
                Some(AnyResult::Wire(Box::new(r)))
            }
            ShapeType::Face => {
                let mut r = FaceRes::new(brep, the_shape);
                // OCCT L391.
                r.set_geometric_controls(b);
                Some(AnyResult::Face(Box::new(r)))
            }
            ShapeType::Shell => Some(AnyResult::Shell(Box::new(ShellRes::new(brep, the_shape)))),
            ShapeType::Solid => Some(AnyResult::Solid(Box::new(SolidRes::new(brep, the_shape)))),
            ShapeType::CompSolid | ShapeType::Compound | ShapeType::Shape => None,
        };

        // OCCT L406-410.
        if let Some(mut hr) = hr {
            if let Some(base) = hr.base_mut() {
                base.my_is_parallel = self.my_is_parallel;
            }
            let key = brep_check_result::ShapeKey::of(the_shape);
            self.my_map_index.insert(key, self.my_map.len());
            self.my_map.push((the_shape.clone(), hr));
        } else {
            // OCCT L410: myMap.Add(theShape, HR) — OCCT also stores the null
            // handle (COMPSOLID/COMPOUND); the rcad map registers a null
            // marker so Put never revisits the subtree.
            let key = brep_check_result::ShapeKey::of(the_shape);
            self.my_map_index.insert(key, self.my_map.len());
            self.my_map.push((the_shape.clone(), AnyResult::Null));
        }

        // OCCT L412-415: the recursion over the children.
        for it in brep_check_result::iterator_subshapes(brep, the_shape) {
            self.put(brep, &it, b); // performs minimum on each shape
        }
    }

    /// OCCT BRepCheck_Analyzer::Perform (Analyzer.cxx L420-454).
    ///
    /// The OCCT body splits the shapes over thread tasks and runs the
    /// `BRepCheck_ParallelAnalyzer` functor through `OSD_Parallel::For` —
    /// with the default `!myIsParallel` the functor runs sequentially in the
    /// map order, which is what this port does.
    fn perform(&mut self, brep: &BRep) {
        // OCCT L422-453: the task-size arithmetic is thread bookkeeping and
        // does not change the sequential result.
        for i in 0..self.my_map.len() {
            let (a_shape, _) = {
                let (s, _) = &self.my_map[i];
                (s.clone(), ())
            };
            self.perform_shape(brep, &a_shape);
        }
    }

    /// OCCT BRepCheck_ParallelAnalyzer::operator() body for one shape
    /// (Analyzer.cxx L56-339, the per-shape switch).
    fn perform_shape(&mut self, brep: &BRep, a_shape: &Shape) {
        let a_type = a_shape.shape_type();
        let a_result_key = brep_check_result::ShapeKey::of(a_shape);
        match a_type {
            ShapeType::Vertex => {
                // OCCT L68-76: no check needed.
            }
            ShapeType::Edge => {
                // OCCT L77-95: CheckPolygonOnTriangulation.
                let ste = {
                    let res = self
                        .my_map_index
                        .get(&a_result_key)
                        .map(|&i| &mut self.my_map[i].1);
                    match res {
                        Some(AnyResult::Edge(r)) => {
                            r.check_polygon_on_triangulation(brep, a_shape)
                        }
                        _ => BRepCheckStatus::NoError,
                    }
                };
                if ste != BRepCheckStatus::NoError {
                    if let Some(&i) = self.my_map_index.get(&a_result_key) {
                        if let AnyResult::Edge(r) = &mut self.my_map[i].1 {
                            r.set_status(ste);
                        }
                    }
                }

                // OCCT L97-124: the vertices in the edge context.
                let mut map_s: std::collections::HashSet<
                    brep_check_result::ShapeKey,
                > = std::collections::HashSet::new();
                for a_vertex in brep_check_result::explorer(brep, a_shape, ShapeType::Vertex)
                {
                    let a_res_of_vertex_key =
                        brep_check_result::ShapeKey::of(&a_vertex);
                    if map_s.insert(a_res_of_vertex_key) {
                        let res = self
                            .my_map_index
                            .get(&a_res_of_vertex_key)
                            .map(|&i| &mut self.my_map[i].1);
                        if let Some(a_res_of_vertex) = res {
                            a_res_of_vertex.in_context(brep, a_shape);
                        }
                    }
                }
            }
            ShapeType::Wire => {
                // OCCT L127-129: nothing.
            }
            ShapeType::Face => {
                // OCCT L130-156: the face vertices in the face context.
                let mut map_s: std::collections::HashSet<
                    brep_check_result::ShapeKey,
                > = std::collections::HashSet::new();
                for a_face_vertex in
                    brep_check_result::explorer(brep, a_shape, ShapeType::Vertex)
                {
                    let key = brep_check_result::ShapeKey::of(&a_face_vertex);
                    if map_s.insert(key) {
                        let res = self
                            .my_map_index
                            .get(&key)
                            .map(|&i| &mut self.my_map[i].1);
                        if let Some(a_face_vertex_res) = res {
                            a_face_vertex_res.in_context(brep, a_shape);
                        }
                    }
                }

                // OCCT L158-210: the face edges in the face context with the
                // `performwire` downgrade.
                let mut performwire = true;
                let is_invalid_tolerance = false; // OCCT L159: never set true.
                map_s.clear();
                for a_face_edge in
                    brep_check_result::explorer(brep, a_shape, ShapeType::Edge)
                {
                    let key = brep_check_result::ShapeKey::of(&a_face_edge);
                    if map_s.insert(key) {
                        let res = self
                            .my_map_index
                            .get(&key)
                            .map(|&i| &mut self.my_map[i].1);
                        if let Some(a_face_edge_res) = res {
                            a_face_edge_res.in_context(brep, a_shape);

                            // OCCT L171-194.
                            if performwire {
                                let base_has = a_face_edge_res
                                    .base()
                                    .map(|b| b.is_status_on_shape(a_shape))
                                    .unwrap_or(false);
                                if base_has {
                                    let statuses = a_face_edge_res
                                        .base()
                                        .map(|b| b.status_on_shape(a_shape).clone())
                                        .unwrap_or_default();
                                    for ste in statuses {
                                        if ste == BRepCheckStatus::NoCurveOnSurface
                                            || ste
                                                == BRepCheckStatus::InvalidCurveOnSurface
                                            || ste == BRepCheckStatus::InvalidRange
                                            || ste
                                                == BRepCheckStatus::InvalidCurveOnClosedSurface
                                        {
                                            performwire = false;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L212-258: the face wires in the face context with the
                // `orientofwires` downgrade.
                let mut orientofwires = performwire;
                for a_face_wire in
                    brep_check_result::explorer(brep, a_shape, ShapeType::Wire)
                {
                    let key = brep_check_result::ShapeKey::of(&a_face_wire);
                    let res = self.my_map_index.get(&key).map(|&i| &mut self.my_map[i].1);
                    if let Some(a_face_wire_res) = res {
                        a_face_wire_res.in_context(brep, a_shape);

                        // OCCT L221-243.
                        if orientofwires {
                            let base_has = a_face_wire_res
                                .base()
                                .map(|b| b.is_status_on_shape(a_shape))
                                .unwrap_or(false);
                            if base_has {
                                let statuses = a_face_wire_res
                                    .base()
                                    .map(|b| b.status_on_shape(a_shape).clone())
                                    .unwrap_or_default();
                                for ste in statuses {
                                    if ste != BRepCheckStatus::NoError {
                                        orientofwires = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L260-302: the face result update.
                let key = brep_check_result::ShapeKey::of(a_shape);
                let res = self.my_map_index.get(&key).map(|&i| &mut self.my_map[i].1);
                if let Some(AnyResult::Face(a_face_res)) = res {
                    if is_invalid_tolerance {
                        a_face_res.set_status(BRepCheckStatus::InvalidToleranceValue);
                    } else if performwire {
                        if orientofwires {
                            a_face_res.orientation_of_wires(brep, true); // on enregistre
                        } else {
                            a_face_res.set_unorientable();
                        }
                    } else {
                        a_face_res.set_unorientable();
                    }
                }
            }
            ShapeType::Shell => {
                // OCCT L305-307: nothing.
            }
            ShapeType::Solid => {
                // OCCT L308-333: the shells in the solid context.
                for a_shell in brep_check_result::explorer(brep, a_shape, ShapeType::Shell)
                {
                    let key = brep_check_result::ShapeKey::of(&a_shell);
                    let res = self.my_map_index.get(&key).map(|&i| &mut self.my_map[i].1);
                    if let Some(a_solid_res) = res {
                        a_solid_res.in_context(brep, a_shape);
                    }
                }
            }
            _ => {
                // OCCT L335-337: default — nothing.
            }
        }
    }

    /// OCCT BRepCheck_Analyzer::IsValid(S) (Analyzer.cxx L458-505).
    pub fn is_valid(&self, brep: &BRep, s: &Shape) -> bool {
        if let Some(res) = self.result(s) {
            if let Some(base) = res.base() {
                let statuses = base.status();
                // OCCT L462-467: itl.Value() — the first status of the list.
                if let Some(&first) = statuses.first() {
                    if first != BRepCheckStatus::NoError {
                        // a voir
                        return false;
                    }
                }
            }
        }

        // OCCT L470-476: the direct children.
        for it in brep_check_result::iterator_subshapes(brep, s) {
            if !self.is_valid(brep, &it) {
                return false;
            }
        }

        // OCCT L478-502.
        match s.shape_type() {
            ShapeType::Edge => return self.valid_sub(brep, s, ShapeType::Vertex),
            ShapeType::Face => {
                let mut valid = self.valid_sub(brep, s, ShapeType::Wire);
                valid = valid && self.valid_sub(brep, s, ShapeType::Edge);
                valid = valid && self.valid_sub(brep, s, ShapeType::Vertex);
                return valid;
            }
            ShapeType::Shell => {
                // OCCT L492-494: commented out in OCCT.
            }
            ShapeType::Solid => {
                // OCCT L495-499.
                return self.valid_sub(brep, s, ShapeType::Shell);
            }
            _ => {}
        }

        true
    }

    /// OCCT BRepCheck_Analyzer::IsValid() (Analyzer.hxx L147).
    pub fn is_valid_all(&self, brep: &BRep) -> bool {
        self.is_valid(brep, &self.my_shape.clone())
    }

    /// OCCT BRepCheck_Analyzer::ValidSub(S, SubType) (Analyzer.cxx L509-539).
    fn valid_sub(&self, brep: &BRep, s: &Shape, sub_type: ShapeType) -> bool {
        for exp in brep_check_result::explorer(brep, s, sub_type) {
            let Some(rv) = self.result(&exp) else {
                continue;
            };
            // OCCT L517-523: the context iterator scan for S (OCCT mutates
            // the iterator through the const handle; rcad keeps the iterator
            // position in a Cell).
            let has_context = {
                let base = match rv.base() {
                    Some(b) => b,
                    None => continue,
                };
                base.init_context_iterator();
                let mut found_ctx = false;
                while base.more_shape_in_context() {
                    if base.contextual_shape().is_same(s) {
                        found_ctx = true;
                        break;
                    }
                    base.next_shape_in_context();
                }
                found_ctx
            };

            // OCCT L525-528: if (!RV->MoreShapeInContext()) break.
            if !has_context {
                break;
            }

            // OCCT L530-536.
            if let Some(base) = rv.base() {
                for itl in base.status_on_shape_in_context() {
                    if *itl != BRepCheckStatus::NoError {
                        return false;
                    }
                }
            }
        }
        true
    }
}

impl BRepCheckAnalyzer {
    /// Read access to the result map (the diagnostics helper).
    pub fn map_ref(&self) -> &[(Shape, AnyResult)] {
        &self.my_map
    }
}

/// Convenience: run the analyzer on the whole BRep pool for the root solid
/// and report the validity (the common rcad entry form).
pub fn analyze_brep_valid(brep: &BRep, root: &Shape) -> bool {
    let analyzer = BRepCheckAnalyzer::new(brep, root, true);
    analyzer.is_valid_all(brep)
}

/// The first TShape::Solid of the pool as the root shape (the rcad BRep
/// always stores its result solid; OCCT callers pass the TopoDS_Shape
/// directly).
pub fn brep_root_shape(brep: &BRep) -> Option<Shape> {
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(ts.as_ref(), TShape::Solid(_)) {
            return Some(Shape {
                data: ts.clone(),
                index: i,
                location: 0,
                orientation: rcad_kernel::topods::Orientation::Forward,
            });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_box() -> (rcad_kernel::BRep, Shape) {
        let brep = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            1.0,
            1.0,
            1.0,
        )
        .expect("make_box_brep");
        let root = brep_root_shape(&brep).expect("root solid");
        (brep, root)
    }

    /// Test 1: a unit box must be VALID through the OCCT Perform/IsValid
    /// structure.
    #[test]
    fn unit_box_is_valid() {
        let (brep, root) = unit_box();
        let analyzer = BRepCheckAnalyzer::new(&brep, &root, true);
        assert!(
            analyzer.is_valid(&brep, &root),
            "the unit box must be valid; statuses: {:?}",
            collect_status_summary(&brep, &analyzer)
        );
    }

    /// Test 2: a translated box stays VALID.
    #[test]
    fn transformed_box_is_valid() {
        let (mut brep, _root) = unit_box();
        let mut t = rcad_kernel::math::gp::Trsf::identity();
        t.set_translation_part(glam::DVec3::new(10.0, -5.0, 2.5));
        crate::topalgo::brep_builderapi_transform::transform_brep(&mut brep, &t);
        let root = brep_root_shape(&brep).expect("root solid");
        let analyzer = BRepCheckAnalyzer::new(&brep, &root, true);
        assert!(
            analyzer.is_valid(&brep, &root),
            "the transformed box must be valid; statuses: {:?}",
            collect_status_summary(&brep, &analyzer)
        );
    }

    /// Test 3: the BRepCheck_Wire OCCT-equivalent API shape on a real face
    /// wire: Closed / Orientation / Closed2d return NoError, and the result
    /// objects are reachable through the analyzer map.
    #[test]
    fn wire_subcheck_api_on_box_face() {
        let (brep, root) = unit_box();
        // Pick the first face and its outer wire (the traversal yields the
        // same composed forms the analyzer uses).
        let shell = root
            .as_solid()
            .and_then(|sd| sd.shells.first().cloned())
            .expect("solid shell");
        let face = shell
            .as_shell()
            .and_then(|shd| shd.faces.first().cloned())
            .expect("shell face");
        let wire = face
            .as_face()
            .map(|fd| fd.outer_wire.clone())
            .expect("outer wire");

        let mut chk = BRepCheckWire::new(&brep, &wire);
        assert!(
            chk.base
                .status()
                .iter()
                .all(|st| *st == BRepCheckStatus::NoError),
            "the wire minimum must be NoError"
        );
        assert_eq!(chk.closed(&brep, false), BRepCheckStatus::NoError);
        assert_eq!(
            chk.orientation(&brep, &face, false),
            BRepCheckStatus::NoError
        );
        assert_eq!(chk.closed2d(&brep, &face, false), BRepCheckStatus::NoError);

        // The analyzer structure exposes the wire result with the shape key.
        let analyzer = BRepCheckAnalyzer::new(&brep, &root, true);
        assert!(analyzer.result(&wire).is_some(), "wire result registered");
    }

    /// Diagnostics helper for the failure messages.
    fn collect_status_summary(brep: &rcad_kernel::BRep, analyzer: &BRepCheckAnalyzer) -> Vec<String> {
        let mut out = Vec::new();
        for (s, res) in analyzer.map_ref() {
            let Some(base) = res.base() else {
                continue;
            };
            let list = base.status();
            let bad: Vec<&BRepCheckStatus> =
                list.iter().filter(|st| **st != BRepCheckStatus::NoError).collect();
            if !bad.is_empty() {
                out.push(format!(
                    "{:?} (type {:?}): {:?}",
                    s.index,
                    s.shape_type(),
                    bad
                ));
            }
        }
        let _ = brep;
        out
    }
}
