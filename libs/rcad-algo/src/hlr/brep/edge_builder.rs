//! OCCT HLRBRep_EdgeBuilder (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_EdgeBuilder.hxx` (L31-113) + `.cxx`
//! (L32-521).  At creation the EdgeBuilder explores the VertexList and
//! uses it to build a list of "AreaLimit" on the edge; an area is a part
//! of the curve between two consecutive vertices.
//!
//! `occ::handle<HLRBRep_AreaLimit>` maps to [`SharedAreaLimit`]
//! (`Arc<AreaLimit>`, see the area_limit module notes on interior
//! mutability); the OCCT handle comparison `cur == myLimits` maps to
//! pointer equality.

use std::sync::Arc;

use rcad_kernel::topods::{Orientation, State};

use crate::hlr::algo::intersection::Intersection;

use super::area_limit::{AreaLimit, SharedAreaLimit};
use super::edge_interference_tool::Data;
use super::vertex_list::VertexList;

/// OCCT handle equality (`h1 == h2`): pointer identity, with
/// null == null true (the OCCT operator compares the pointers).
fn same_handle(a: &Option<SharedAreaLimit>, b: &Option<SharedAreaLimit>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => Arc::ptr_eq(x, y),
        (None, None) => true,
        _ => false,
    }
}

/// OCCT HLRBRep_EdgeBuilder.
pub struct EdgeBuilder {
    to_build: State,
    my_limits: Option<SharedAreaLimit>,
    left: Option<SharedAreaLimit>,
    right: Option<SharedAreaLimit>,
    current: i32,
}

impl EdgeBuilder {
    /// OCCT HLRBRep_EdgeBuilder(VList) — cxx L32-230: creates an
    /// EdgeBuilder algorithm.  `VList` describes the edge and the
    /// interferences.  AreaLimits are created from the vertices.
    /// Builds(IN) is automatically called.
    pub fn new<'a, D: Data + ?Sized>(vlist: &mut VertexList<'a, D>) -> Self {
        // at creation the EdgeBuilder explore the VertexList
        // and use it to build a list of "AreaLimit" on the edge.
        // An area is a part of the curve between
        // two consecutive vertices
        if !vlist.more() {
            panic!("Standard_DomainError: EdgeBuilder  : Empty vertex list");
        }

        let mut eb = EdgeBuilder {
            to_build: State::Unknown,
            my_limits: None,
            left: None,
            right: None,
            // OCCT leaves `current` uninitialized; Builds() sets it before
            // any read.  The neutral default is 0.
            current: 0,
        };
        let mut last: Option<SharedAreaLimit> = None;
        let mut cur: Option<SharedAreaLimit>;

        // loop on the Vertices
        while vlist.more() {
            let mut before = State::Unknown;
            let mut after = State::Unknown;
            let mut ebefore = State::Unknown;
            let mut eafter = State::Unknown;
            // compute the states
            if vlist.is_boundary() {
                match vlist.orientation() {
                    Orientation::Forward => {
                        ebefore = State::Out;
                        eafter = State::In;
                    }
                    Orientation::Reversed => {
                        ebefore = State::In;
                        eafter = State::Out;
                    }
                    Orientation::Internal => {
                        ebefore = State::In;
                        eafter = State::In;
                    }
                    Orientation::External => {
                        ebefore = State::Out;
                        eafter = State::Out;
                    }
                }
            }

            if vlist.is_interference() {
                match vlist.transition() {
                    Orientation::Forward => {
                        before = State::Out;
                        after = State::In;
                    }
                    Orientation::Reversed => {
                        before = State::In;
                        after = State::Out;
                    }
                    Orientation::Internal => {
                        before = State::In;
                        after = State::In;
                    }
                    Orientation::External => {
                        before = State::Out;
                        after = State::Out;
                    }
                }

                match vlist.boundary_transition() {
                    Orientation::Forward => {
                        after = State::On;
                    }
                    Orientation::Reversed => {
                        before = State::On;
                    }
                    Orientation::Internal => {
                        before = State::On;
                        after = State::On;
                    }
                    Orientation::External => {}
                }
            }

            // create the Limit and connect to list
            let v: Intersection = *vlist.current();
            cur = Some(Arc::new(AreaLimit::new(
                v,
                vlist.is_boundary(),
                vlist.is_interference(),
                before,
                after,
                ebefore,
                eafter,
            )));
            if eb.my_limits.is_none() {
                eb.my_limits = cur.clone();
                last = cur;
            } else {
                last.as_ref().expect("last limit").set_next(cur.clone());
                cur
                    .as_ref()
                    .expect("current limit")
                    .set_previous(last.clone());
                last = cur;
            }
            vlist.next();
        }

        // periodicity, make a circular list
        if vlist.is_periodic() {
            last
                .as_ref()
                .expect("last limit")
                .set_next(eb.my_limits.clone());
            eb.my_limits
                .as_ref()
                .expect("first limit")
                .set_previous(last);
        }

        // process UNKNOWN areas
        let mut stat = State::Unknown;
        let mut estat = State::Unknown;

        cur = eb.my_limits.clone();
        while let Some(c) = cur {
            if stat == State::Unknown {
                stat = c.state_before();
                if stat == State::Unknown {
                    stat = c.state_after();
                }
            }
            if estat == State::Unknown {
                estat = c.edge_before();
                if estat == State::Unknown {
                    estat = c.edge_after();
                }
            }
            cur = c.next();
            // test for periodicicity
            if same_handle(&cur, &eb.my_limits) {
                break;
            }
        }

        // error if no interferences
        if stat == State::Unknown {
            panic!("Standard_DomainError: EdgeBuilder : No interferences");
        }
        // if no boundary the edge covers the whole curve
        if estat == State::Unknown {
            estat = State::In;
        }

        // propagate states
        cur = eb.my_limits.clone();
        while let Some(c) = cur {
            if c.state_before() == State::Unknown {
                c.set_state_before(stat);
            } else {
                stat = c.state_after();
            }
            if c.state_after() == State::Unknown {
                c.set_state_after(stat);
            }
            if c.edge_before() == State::Unknown {
                c.set_edge_before(estat);
            } else {
                estat = c.edge_after();
            }
            if c.edge_after() == State::Unknown {
                c.set_edge_after(estat);
            }

            cur = c.next();
            if same_handle(&cur, &eb.my_limits) {
                break;
            }
        }

        // initialise with IN parts
        eb.builds(State::In);
        eb
    }

    /// OCCT InitAreas — cxx L234-238: initialize an iteration on the areas.
    pub fn init_areas(&mut self) {
        self.left = self
            .my_limits
            .as_ref()
            .and_then(|m| m.previous());
        self.right = self.my_limits.clone();
    }

    /// OCCT NextArea — cxx L242-249: set the current area to the next area.
    pub fn next_area(&mut self) {
        self.left = self.right.clone();
        if let Some(r) = self.right.clone() {
            self.right = r.next();
        }
    }

    /// OCCT PreviousArea — cxx L253-260: set the current area to the
    /// previous area.
    pub fn previous_area(&mut self) {
        self.right = self.left.clone();
        if let Some(l) = self.left.clone() {
            self.left = l.previous();
        }
    }

    /// OCCT HasArea — cxx L264-278: returns true if there is a current
    /// area.
    pub fn has_area(&self) -> bool {
        if self.left.is_none() {
            if self.right.is_none() {
                return false;
            }
        }
        if same_handle(&self.right, &self.my_limits) {
            return false;
        }
        true
    }

    /// OCCT AreaState — cxx L282-294: the state of the current area.
    pub fn area_state(&self) -> State {
        let mut stat = State::Unknown;
        if let Some(l) = &self.left {
            stat = l.state_after();
        }
        if let Some(r) = &self.right {
            stat = r.state_before();
        }
        stat
    }

    /// OCCT AreaEdgeState — cxx L298-310: the edge state of the current
    /// area.
    pub fn area_edge_state(&self) -> State {
        let mut stat = State::Unknown;
        if let Some(l) = &self.left {
            stat = l.edge_after();
        }
        if let Some(r) = &self.right {
            stat = r.edge_before();
        }
        stat
    }

    /// OCCT LeftLimit — cxx L314-317: the AreaLimit beginning the current
    /// area.  This is a NULL handle when the area is infinite on the left.
    pub fn left_limit(&self) -> Option<SharedAreaLimit> {
        self.left.clone()
    }

    /// OCCT RightLimit — cxx L321-324: the AreaLimit ending the current
    /// area.  This is a NULL handle when the area is infinite on the right.
    pub fn right_limit(&self) -> Option<SharedAreaLimit> {
        self.right.clone()
    }

    /// OCCT Builds(ToBuild) — cxx L328-349: reinitialize the results
    /// iteration to the parts with State `ToBuild`.  If this method is not
    /// called after construction the default is `ToBuild` = IN.
    pub fn builds(&mut self, to_build: State) {
        self.to_build = to_build;
        self.init_areas();
        loop {
            if (self.area_state() == self.to_build) && (self.area_edge_state() == State::In) {
                if self.left.is_none() {
                    self.current = 2;
                } else {
                    self.current = 1;
                }
                return;
            }
            self.next_area();
            if !self.has_area() {
                break;
            }
        }
        self.current = 3;
    }

    /// OCCT MoreEdges — cxx L353-356: true if there are more new edges to
    /// build.
    pub fn more_edges(&self) -> bool {
        self.has_area()
    }

    /// OCCT NextEdge — cxx L360-384: proceeds to the next edge to build.
    /// Skip all remaining vertices on the current edge.
    pub fn next_edge(&mut self) {
        // clean the current edge
        while self.area_state() == self.to_build {
            self.next_area();
        }
        // go to the next edge
        while self.has_area() {
            if (self.area_state() == self.to_build)
                && (self.area_edge_state() == State::In)
            {
                if self.left.is_none() {
                    self.current = 2;
                } else {
                    self.current = 1;
                }
                return;
            }
            self.next_area();
        }
    }

    /// OCCT MoreVertices — cxx L388-391: true if there are more vertices in
    /// the current new edge.
    pub fn more_vertices(&self) -> bool {
        self.current < 3
    }

    /// OCCT NextVertex — cxx L395-421: proceeds to the next vertex of the
    /// current edge.
    pub fn next_vertex(&mut self) {
        if self.current == 1 {
            self.current = 2;
            if self.right.is_none() {
                self.current = 3;
            }
        } else if self.current == 2 {
            self.next_area();
            if (self.area_state() == self.to_build)
                && (self.area_edge_state() == State::In)
            {
                self.current = 2;
            } else {
                self.current = 3;
            }
        } else {
            panic!("Standard_NoSuchObject: EdgeBuilder::NextVertex : No current edge");
        }
    }

    /// OCCT Current — cxx L425-439: the current vertex of the current edge.
    pub fn current(&self) -> &Intersection {
        if self.current == 1 {
            self.left.as_ref().expect("left limit").vertex()
        } else if self.current == 2 {
            self.right.as_ref().expect("right limit").vertex()
        } else {
            panic!("Standard_NoSuchObject: EdgeBuilder::Current : No current vertex");
        }
    }

    /// OCCT IsBoundary — cxx L443-457: true if the current vertex comes
    /// from the boundary of the edge.
    pub fn is_boundary(&self) -> bool {
        if self.current == 1 {
            self.left.as_ref().expect("left limit").is_boundary()
        } else if self.current == 2 {
            self.right.as_ref().expect("right limit").is_boundary()
        } else {
            panic!("Standard_NoSuchObject: EdgeBuilder::IsBoundary : No current vertex");
        }
    }

    /// OCCT IsInterference — cxx L461-475: true if the current vertex was
    /// an interference.
    pub fn is_interference(&self) -> bool {
        if self.current == 1 {
            self.left.as_ref().expect("left limit").is_interference()
        } else if self.current == 2 {
            self.right.as_ref().expect("right limit").is_interference()
        } else {
            panic!("Standard_NoSuchObject: EdgeBuilder::IsInterference : No current vertex");
        }
    }

    /// OCCT Orientation — cxx L479-505: the new orientation of the current
    /// vertex.
    pub fn orientation(&self) -> Orientation {
        if self.current == 1 {
            let l = self.left.as_ref().expect("left limit");
            if (l.state_before() == l.state_after()) && (l.edge_before() == l.edge_after()) {
                Orientation::Internal
            } else {
                Orientation::Forward
            }
        } else if self.current == 2 {
            let r = self.right.as_ref().expect("right limit");
            if (r.state_before() == r.state_after()) && (r.edge_before() == r.edge_after()) {
                Orientation::Internal
            } else {
                Orientation::Reversed
            }
        } else {
            Orientation::External // only for WNT.
        }
    }

    /// OCCT Destroy — cxx L509-521.
    pub fn destroy(&mut self) {
        let mut cur = self.my_limits.clone();
        while let Some(c) = cur {
            let n = c.next();
            c.clear();
            cur = n;
        }
        self.left = None;
        self.right = None;
        self.my_limits = None;
    }
}

impl Drop for EdgeBuilder {
    /// OCCT ~HLRBRep_EdgeBuilder() { Destroy(); } — hxx L105.
    fn drop(&mut self) {
        self.destroy();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::algo::coincidence::Coincidence;
    use crate::hlr::algo::interference::Interference;
    use crate::hlr::brep::edge_data::EdgeData;
    use crate::hlr::brep::edge_interference_tool::EdgeInterferenceTool;
    use crate::hlr::brep::edge_ilist::EdgeIList;
    use crate::hlr::brep::vertex_list::VertexList;
    use glam::DVec2;

    /// Minimal HLRBRep_Data stand-in: one edge with bounds [0, 1].
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
            _tg: &mut DVec2,
            _nm: &mut DVec2,
            _cr: &mut f64,
        ) {
        }
        fn local_fe_geometry_2d(
            &self,
            _fe: i32,
            _param: f64,
            _tg: &mut DVec2,
            _nm: &mut DVec2,
            _cr: &mut f64,
        ) {
        }
    }

    fn interference_at(
        param: f64,
        index: i32,
        trans: Orientation,
        b_trans: Orientation,
    ) -> Interference {
        let mut inters = Intersection::new();
        inters.set_parameter(param);
        inters.set_index(index);
        // The intersection orientation is INTERNAL so that OCCT's
        // SameVertexAndInterference (cxx L77-84) never merges the
        // interference with an edge boundary vertex (the merge is exercised
        // by the OCCT semantics themselves, not by these anchors).
        inters.set_orientation(Orientation::Internal);
        let mut bound = Coincidence::new();
        bound.set_2d(1, param);
        Interference::from_parts(inters, bound, Orientation::Forward, trans, b_trans)
    }

    /// The OCCT Hider::Hide sequence on a loaded edge: LoadEdge, sort the
    /// interferences, merge the vertices, build.
    fn build_for(data: &TestData, il: &mut Vec<Interference>) -> EdgeBuilder {
        let mut tool = EdgeInterferenceTool::new(data);
        tool.load_edge();
        let sorted: Vec<Interference> = {
            let mut s: Vec<Interference> = Vec::new();
            for i in il.drain(..) {
                EdgeIList::add_interference(&mut s, &i, &tool);
            }
            s
        };
        let mut vl = VertexList::new(&tool, &sorted);
        EdgeBuilder::new(&mut vl)
    }

    /// One interference at 0.4 entering the face (FORWARD transition, no
    /// boundary patch) on an edge with boundary vertices at 0 (FORWARD) and
    /// 1 (REVERSED): the IN parts of the curve are [0.4, 1.0] — one built
    /// edge with the vertices 0.4 (interference) and 1.0 (boundary,
    /// REVERSED).
    #[test]
    fn edge_builder_splits_in_parts() {
        let data = TestData::new();
        let mut il = vec![interference_at(
            0.4,
            7,
            Orientation::Forward,
            Orientation::External,
        )];
        let mut b = build_for(&data, &mut il);

        // Builds(IN) ran in the constructor: first area with state IN and
        // edge state IN starts at 0.4.
        assert!(b.more_edges());
        assert!(b.more_vertices());
        // current == 1 (left = the 0.4 limit).
        assert!(b.left_limit().is_some());
        assert!(b.is_interference());
        assert!(!b.is_boundary());
        assert_eq!(b.current().parameter(), 0.4);
        assert_eq!(b.current().index(), 7);
        // The interference entered the face: FORWARD orientation.
        assert_eq!(b.orientation(), Orientation::Forward);
        b.next_vertex();
        // current == 2 (right = the 1.0 boundary vertex).
        assert!(b.more_vertices());
        assert!(!b.is_interference());
        assert!(b.is_boundary());
        assert_eq!(b.current().parameter(), 1.0);
        assert_eq!(b.current().index(), 22);
        // State before == after (IN/IN) but edge before != after -> REVERSED.
        assert_eq!(b.orientation(), Orientation::Reversed);
        b.next_vertex();
        assert!(!b.more_vertices());
        b.next_edge();
        assert!(!b.more_edges());
    }

    /// The same area list rebuilt for OUT: the OUT part with edge state IN
    /// is [0.0, 0.4] — the vertices 0.0 (boundary, FORWARD) and 0.4.
    #[test]
    fn edge_builder_builds_out_part() {
        let data = TestData::new();
        let mut il = vec![interference_at(
            0.4,
            7,
            Orientation::Forward,
            Orientation::External,
        )];
        let mut b = build_for(&data, &mut il);
        b.builds(State::Out);
        assert!(b.more_edges());
        assert!(b.more_vertices());
        assert!(b.is_boundary());
        assert_eq!(b.current().parameter(), 0.0);
        assert_eq!(b.current().index(), 11);
        assert_eq!(b.orientation(), Orientation::Forward);
        b.next_vertex();
        assert!(b.is_interference());
        assert_eq!(b.current().parameter(), 0.4);
        b.next_vertex();
        assert!(!b.more_vertices());
        b.next_edge();
        assert!(!b.more_edges());
    }

    /// Area walk API: InitAreas / NextArea / PreviousArea / HasArea /
    /// AreaState / AreaEdgeState / LeftLimit / RightLimit (cxx L234-324).
    #[test]
    fn edge_builder_area_walk() {
        let data = TestData::new();
        let mut il = vec![
            interference_at(0.3, 7, Orientation::Forward, Orientation::Forward),
            interference_at(0.6, 8, Orientation::Reversed, Orientation::Forward),
        ];
        let mut b = build_for(&data, &mut il);
        // Propagated areas: [0, 0.3] OUT, [0.3, 0.6] IN (the FORWARD
        // boundary transitions turn the states after 0.3 and 0.6 into ON
        // patches, overridden by the next limit's state-before).
        b.init_areas();
        // Area before 0.0: left null, right = the 0.0 boundary limit.
        assert!(b.left_limit().is_none());
        assert!(b.right_limit().is_some());
        // OCCT HasArea (cxx L264-278) is false while right == myLimits —
        // the first area is only visited by the Builds/NextEdge loops.
        assert!(!b.has_area());
        assert_eq!(b.area_edge_state(), State::Out);
        b.next_area();
        // [0.0, 0.3]: OUT, edge IN.
        assert_eq!(b.area_state(), State::Out);
        assert_eq!(b.area_edge_state(), State::In);
        b.next_area();
        // [0.3, 0.6]: IN, edge IN.
        assert_eq!(b.area_state(), State::In);
        assert_eq!(b.area_edge_state(), State::In);
        // PreviousArea walks back.
        b.previous_area();
        assert_eq!(b.area_state(), State::Out);
    }

    /// OCCT raises when no interference gives a state
    /// (cxx L185, Standard_DomainError_Raise_if).
    #[test]
    #[should_panic(expected = "No interferences")]
    fn edge_builder_no_interferences_raises() {
        let data = TestData::new();
        let mut tool = EdgeInterferenceTool::new(&data);
        tool.load_edge();
        let il: Vec<Interference> = Vec::new();
        // Boundary vertices only: after propagation stat stays UNKNOWN.
        let mut vl = VertexList::new(&tool, &il);
        let _ = EdgeBuilder::new(&mut vl);
    }
}
