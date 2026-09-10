//! OCCT HLRBRep_EdgeInterferenceTool (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_EdgeInterferenceTool.hxx` (L33-88) +
//! `.cxx` (L26-100) + `.lxx` (L21-77).  Implements the methods required to
//! instantiate the EdgeInterferenceList from HLRAlgo.
//!
//! Documented deferral: OCCT passes `occ::handle<HLRBRep_Data>`; the
//! HLRBRep_Data structure lands with Stage 3f.  The tool is therefore
//! generic over the [`Data`] trait, which mirrors exactly the four
//! HLRBRep_Data members consumed here (EDataArray, Edge,
//! LocalLEGeometry2D, LocalFEGeometry2D); HLRBRep_Data will implement it
//! when it lands.

use glam::{DVec2, DVec3};
use rcad_kernel::topods::Orientation;

use crate::hlr::algo::intersection::Intersection;
use crate::hlr::algo::interference::Interference;
use super::edge_data::EdgeData;

/// OCCT `HLRBRep_Data` (the slice consumed by this tool).  The Data array
/// access `ED(myDS->Edge())` is 1-based in OCCT; the rcad Data convention
/// indexes `[edge() - 1]`.
pub trait Data {
    /// OCCT HLRBRep_Data::EDataArray() — the per-edge records.
    fn e_data_array(&self) -> &[EdgeData<'static>];

    /// OCCT HLRBRep_Data::Edge() — the current edge index (1-based).
    fn edge(&self) -> i32;

    /// OCCT HLRBRep_Data::LocalLEGeometry2D(Param, TgLE, NmLE, CrLE).
    fn local_le_geometry_2d(
        &self,
        param: f64,
        tg_le: &mut DVec2,
        nm_le: &mut DVec2,
        cr_le: &mut f64,
    );

    /// OCCT HLRBRep_Data::LocalFEGeometry2D(FE, Param, TgFE, NmFE, CrFE).
    fn local_fe_geometry_2d(
        &self,
        fe: i32,
        param: f64,
        tg_fe: &mut DVec2,
        nm_fe: &mut DVec2,
        cr_fe: &mut f64,
    );
}

/// OCCT HLRBRep_EdgeInterferenceTool.
pub struct EdgeInterferenceTool<'a, D: Data + ?Sized> {
    my_ds: &'a D,
    inter: [Intersection; 2],
    cur: i32,
}

impl<'a, D: Data + ?Sized> Clone for EdgeInterferenceTool<'a, D> {
    fn clone(&self) -> Self {
        EdgeInterferenceTool {
            my_ds: self.my_ds,
            inter: self.inter,
            cur: self.cur,
        }
    }
}

impl<'a, D: Data + ?Sized> Copy for EdgeInterferenceTool<'a, D> {}

impl<'a, D: Data + ?Sized> EdgeInterferenceTool<'a, D> {
    /// OCCT HLRBRep_EdgeInterferenceTool(DS) — cxx L26-29.  `cur` is left
    /// uninitialized in OCCT (InitVertices sets it before any read); the
    /// neutral default is 0.
    pub fn new(ds: &'a D) -> Self {
        EdgeInterferenceTool {
            my_ds: ds,
            inter: [Intersection::new(), Intersection::new()],
            cur: 0,
        }
    }

    /// OCCT LoadEdge — cxx L33-46.
    pub fn load_edge(&mut self) {
        let ds: &'a D = self.my_ds;
        let ed = &ds.e_data_array()[(ds.edge() - 1) as usize];
        let (p1, t1, p2, t2) = ed.status_ref().bounds();
        let v_sta = ed.v_sta();
        let v_end = ed.v_end();
        self.inter[0].set_parameter(p1);
        self.inter[0].set_tolerance(t1);
        self.inter[0].set_index(v_sta);
        self.inter[1].set_parameter(p2);
        self.inter[1].set_tolerance(t2);
        self.inter[1].set_index(v_end);
    }

    /// OCCT EdgeGeometry(Param, Tgt, Nrm, Curv) — cxx L50-59: the local
    /// geometric description of the Edge at parameter `Param`.  See method
    /// Reset of class EdgeFaceTransition from TopCnx for other arguments.
    pub fn edge_geometry(&self, param: f64, tgt: &mut DVec3, nrm: &mut DVec3, cr_le: &mut f64) {
        // OCCT gp_Dir2d locals (uninitialized); the neutral default is 0.
        let mut tg_le = DVec2::ZERO;
        let mut nm_le = DVec2::ZERO;
        self.my_ds
            .local_le_geometry_2d(param, &mut tg_le, &mut nm_le, cr_le);
        *tgt = DVec3::new(tg_le.x, tg_le.y, 0.0);
        *nrm = DVec3::new(nm_le.x, nm_le.y, 0.0);
    }

    /// OCCT SameInterferences(I1, I2) — cxx L63-73: true if the two
    /// interferences are on the same geometric locus.
    pub fn same_interferences(&self, i1: &Interference, i2: &Interference) -> bool {
        let ind1 = i1.intersection().index();
        let ind2 = i2.intersection().index();
        if ind1 != 0 && ind2 != 0 {
            return ind1 == ind2;
        }
        false
    }

    /// OCCT SameVertexAndInterference(I) — cxx L77-84: true if the
    /// Interference and the current Vertex are on the same geometric locus.
    pub fn same_vertex_and_interference(&self, i: &Interference) -> bool {
        if i.intersection().index() == self.inter[self.cur as usize].index() {
            return true;
        }
        i.intersection().orientation() == if self.cur == 0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }

    /// OCCT InterferenceBoundaryGeometry(I, Tang, Norm, Curv) — cxx L88-100:
    /// the geometry of the boundary at the interference `I`.  See the
    /// AddInterference method of the class EdgeFaceTransition from TopCnx
    /// for the other arguments.
    pub fn interference_boundary_geometry(
        &self,
        i: &Interference,
        tang: &mut DVec3,
        norm: &mut DVec3,
        cr_fe: &mut f64,
    ) {
        // OCCT gp_Dir2d locals (uninitialized); the neutral default is 0.
        let mut tg_fe = DVec2::ZERO;
        let mut nm_fe = DVec2::ZERO;
        let (fe, param) = i.boundary().value_2d();
        self.my_ds
            .local_fe_geometry_2d(fe, param, &mut tg_fe, &mut nm_fe, cr_fe);
        *tang = DVec3::new(tg_fe.x, tg_fe.y, 0.0);
        *norm = DVec3::new(nm_fe.x, nm_fe.y, 0.0);
    }

    /// OCCT InitVertices — lxx L21-24.
    pub fn init_vertices(&mut self) {
        self.cur = 0;
    }

    /// OCCT MoreVertices — lxx L28-31.
    pub fn more_vertices(&self) -> bool {
        self.cur < 2
    }

    /// OCCT NextVertex — lxx L35-38.
    pub fn next_vertex(&mut self) {
        self.cur += 1;
    }

    /// OCCT CurrentVertex — lxx L42-45.
    pub fn current_vertex(&self) -> &Intersection {
        &self.inter[self.cur as usize]
    }

    /// OCCT CurrentOrientation — lxx L49-55.
    pub fn current_orientation(&self) -> Orientation {
        if self.cur == 0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }

    /// OCCT CurrentParameter — lxx L59-62.
    pub fn current_parameter(&self) -> f64 {
        self.inter[self.cur as usize].parameter()
    }

    /// OCCT IsPeriodic — lxx L66-69.
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT ParameterOfInterference — lxx L73-77.
    pub fn parameter_of_interference(&self, i: &Interference) -> f64 {
        i.intersection().parameter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::algo::coincidence::Coincidence;
    use crate::hlr::brep::edge_data::EdgeData;

    /// Minimal HLRBRep_Data stand-in: one edge with bounds [0, 1] and the
    /// fixed local geometries used by the anchors.
    struct TestData {
        edges: Vec<EdgeData<'static>>,
        edge: i32,
    }

    impl TestData {
        fn new() -> Self {
            let mut ed = EdgeData::new();
            ed.set_v_sta(11);
            ed.set_v_end(22);
            ed.status().initialize(0.0, 1.0e-4, 1.0, 1.0e-4);
            TestData {
                edges: vec![ed],
                edge: 1,
            }
        }
    }

    impl Data for TestData {
        fn e_data_array(&self) -> &[EdgeData<'static>] {
            &self.edges
        }
        fn edge(&self) -> i32 {
            self.edge
        }
        fn local_le_geometry_2d(
            &self,
            _param: f64,
            tg_le: &mut DVec2,
            nm_le: &mut DVec2,
            cr_le: &mut f64,
        ) {
            *tg_le = DVec2::new(1.0, 0.0);
            *nm_le = DVec2::new(0.0, 1.0);
            *cr_le = 0.5;
        }
        fn local_fe_geometry_2d(
            &self,
            _fe: i32,
            _param: f64,
            tg_fe: &mut DVec2,
            nm_fe: &mut DVec2,
            cr_fe: &mut f64,
        ) {
            *tg_fe = DVec2::new(0.0, -1.0);
            *nm_fe = DVec2::new(1.0, 0.0);
            *cr_fe = 2.0;
        }
    }

    fn interference_at(param: f64, index: i32, orient: Orientation) -> Interference {
        let mut inters = Intersection::new();
        inters.set_parameter(param);
        inters.set_index(index);
        inters.set_orientation(orient);
        let mut bound = Coincidence::new();
        bound.set_2d(3, param * 2.0);
        Interference::from_parts(
            inters,
            bound,
            Orientation::Forward,
            Orientation::Forward,
            Orientation::Forward,
        )
    }

    /// OCCT LoadEdge (cxx L33-46): the tool vertices get the edge bounds
    /// and the VSta/VEnd indices.
    #[test]
    fn edge_interference_tool_load_edge_bounds() {
        let data = TestData::new();
        let mut tool = EdgeInterferenceTool::new(&data);
        tool.load_edge();
        tool.init_vertices();
        assert!(tool.more_vertices());
        assert_eq!(tool.current_parameter(), 0.0);
        assert_eq!(tool.current_vertex().index(), 11);
        assert_eq!(tool.current_orientation(), Orientation::Forward);
        tool.next_vertex();
        assert!(tool.more_vertices());
        assert_eq!(tool.current_parameter(), 1.0);
        assert_eq!(tool.current_vertex().index(), 22);
        assert_eq!(tool.current_orientation(), Orientation::Reversed);
        tool.next_vertex();
        assert!(!tool.more_vertices());
        assert!(!tool.is_periodic());
    }

    /// OCCT SameInterferences / SameVertexAndInterference / ParameterOfInterference
    /// (cxx L63-84, lxx L73-77).
    #[test]
    fn edge_interference_tool_same_semantics() {
        let data = TestData::new();
        let mut tool = EdgeInterferenceTool::new(&data);
        tool.load_edge();
        tool.init_vertices();
        let i_a = interference_at(0.4, 7, Orientation::Forward);
        let i_b = interference_at(0.6, 7, Orientation::Reversed);
        let i_c = interference_at(0.5, 0, Orientation::Forward);
        let i_d = interference_at(0.5, 8, Orientation::Forward);
        // Same non-zero index -> same locus; zero index -> never same.
        assert!(tool.same_interferences(&i_a, &i_b));
        assert!(!tool.same_interferences(&i_a, &i_c));
        assert!(!tool.same_interferences(&i_a, &i_d));
        assert_eq!(tool.parameter_of_interference(&i_a), 0.4);
        // cur = 0: index 11 matches inter[0].index -> same locus.
        assert!(tool.same_vertex_and_interference(&interference_at(0.1, 11, Orientation::Reversed)));
        // cur = 0: index mismatch, orientation must be FORWARD to match.
        assert!(tool.same_vertex_and_interference(&i_a));
        tool.next_vertex();
        // cur = 1: index mismatch, orientation must be REVERSED to match.
        assert!(tool.same_vertex_and_interference(&i_b));
        assert!(!tool.same_vertex_and_interference(&i_a));
        assert!(tool.same_vertex_and_interference(&interference_at(0.9, 22, Orientation::Forward)));
    }

    /// OCCT EdgeGeometry / InterferenceBoundaryGeometry (cxx L50-59,
    /// L88-100): the 2D local geometry is lifted to z = 0.
    #[test]
    fn edge_interference_tool_geometry_lift_to_z0() {
        let data = TestData::new();
        let tool = EdgeInterferenceTool::new(&data);
        let mut tgt = DVec3::ZERO;
        let mut nrm = DVec3::ZERO;
        let mut curv: f64 = 0.0;
        tool.edge_geometry(0.25, &mut tgt, &mut nrm, &mut curv);
        assert_eq!(tgt, DVec3::new(1.0, 0.0, 0.0));
        assert_eq!(nrm, DVec3::new(0.0, 1.0, 0.0));
        assert_eq!(curv, 0.5);
        let i = interference_at(0.4, 7, Orientation::Forward);
        let mut tang = DVec3::ZERO;
        let mut norm = DVec3::ZERO;
        let mut crv: f64 = 0.0;
        tool.interference_boundary_geometry(&i, &mut tang, &mut norm, &mut crv);
        assert_eq!(tang, DVec3::new(0.0, -1.0, 0.0));
        assert_eq!(norm, DVec3::new(1.0, 0.0, 0.0));
        assert_eq!(crv, 2.0);
    }
}
