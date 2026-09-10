//! OCCT HLRBRep_VertexList (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_VertexList.hxx` (L32-76) + `.cxx` (L28-155).
//! An iterator merging the two edge-boundary vertices (held by the tool)
//! with the sorted interference list.  The OCCT
//! `NCollection_List<HLRAlgo_Interference>::Iterator myIterator` maps to a
//! slice reference plus the index position (the landed HLRAlgo precedent
//! for NCollection_List).

use rcad_kernel::topods::Orientation;

use crate::hlr::algo::intersection::Intersection;
use crate::hlr::algo::interference::Interference;

use super::edge_interference_tool::{Data, EdgeInterferenceTool};

/// OCCT HLRBRep_VertexList.
pub struct VertexList<'a, D: Data + ?Sized> {
    /// OCCT myIterator (NCollection_List<HLRAlgo_Interference>::Iterator):
    /// the referenced list plus the current index position.
    my_iterator: &'a [Interference],
    my_index: usize,
    my_tool: EdgeInterferenceTool<'a, D>,
    from_edge: bool,
    from_interf: bool,
}

impl<'a, D: Data + ?Sized> VertexList<'a, D> {
    /// OCCT HLRBRep_VertexList(T, I) — cxx L28-37.  The OCCT iterator `I`
    /// is positioned on the first element of the interference list.
    pub fn new(t: &EdgeInterferenceTool<'a, D>, i: &'a [Interference]) -> Self {
        let mut vl = VertexList {
            my_iterator: i,
            my_index: 0,
            my_tool: *t,
            from_edge: false,
            from_interf: false,
        };
        vl.my_tool.init_vertices();
        vl.next();
        vl
    }

    /// OCCT IsPeriodic — cxx L41-44: returns true when the curve is
    /// periodic.
    pub fn is_periodic(&self) -> bool {
        self.my_tool.is_periodic()
    }

    /// OCCT More — cxx L48-51: returns true when there are more vertices.
    pub fn more(&self) -> bool {
        self.from_edge || self.from_interf
    }

    /// OCCT Next — cxx L55-81: proceeds to the next vertex.
    pub fn next(&mut self) {
        if self.from_interf {
            // myIterator.Next()
            self.my_index += 1;
        }
        if self.from_edge {
            self.my_tool.next_vertex();
        }
        // myIterator.More()
        self.from_interf = self.my_index < self.my_iterator.len();
        self.from_edge = self.my_tool.more_vertices();
        if self.from_edge && self.from_interf {
            if !self
                .my_tool
                .same_vertex_and_interference(&self.my_iterator[self.my_index])
            {
                if self.my_tool.current_parameter()
                    < self
                        .my_tool
                        .parameter_of_interference(&self.my_iterator[self.my_index])
                {
                    self.from_interf = false;
                } else {
                    self.from_edge = false;
                }
            }
        }
    }

    /// OCCT Current — cxx L85-99: returns the current vertex.
    pub fn current(&self) -> &Intersection {
        if self.from_edge {
            self.my_tool.current_vertex()
        } else if self.from_interf {
            self.my_iterator[self.my_index].intersection()
        } else {
            panic!("Standard_NoSuchObject: HLRBRep_VertexList::Current");
        }
    }

    /// OCCT IsBoundary — cxx L103-106: true if the current vertex is on the
    /// boundary of the edge.
    pub fn is_boundary(&self) -> bool {
        self.from_edge
    }

    /// OCCT IsInterference — cxx L110-113: true if the current vertex is an
    /// interference.
    pub fn is_interference(&self) -> bool {
        self.from_interf
    }

    /// OCCT Orientation — cxx L117-127: the orientation of the current
    /// vertex if it is on the boundary of the edge.
    pub fn orientation(&self) -> Orientation {
        if self.from_edge {
            self.my_tool.current_orientation()
        } else {
            panic!("Standard_DomainError: HLRBRep_VertexList::Orientation");
        }
    }

    /// OCCT Transition — cxx L131-141: the transition of the current vertex
    /// if it is an interference.
    pub fn transition(&self) -> Orientation {
        if self.from_interf {
            self.my_iterator[self.my_index].transition()
        } else {
            panic!("Standard_DomainError: HLRBRep_VertexList::Transition");
        }
    }

    /// OCCT BoundaryTransition — cxx L145-155: the transition of the
    /// current vertex relative to the boundary if it is an interference.
    pub fn boundary_transition(&self) -> Orientation {
        if self.from_interf {
            self.my_iterator[self.my_index].boundary_transition()
        } else {
            panic!("Standard_DomainError: HLRBRep_VertexList::BoundaryTransition");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::algo::coincidence::Coincidence;
    use crate::hlr::brep::edge_data::EdgeData;
    use crate::hlr::brep::edge_interference_tool::EdgeInterferenceTool;
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

    fn interference_at(param: f64, index: i32) -> Interference {
        let mut inters = Intersection::new();
        inters.set_parameter(param);
        inters.set_index(index);
        // INTERNAL intersection orientation: SameVertexAndInterference never
        // merges the interference with a boundary vertex, so the traversal
        // exercises the parameter ordering branch.
        inters.set_orientation(Orientation::Internal);
        let mut bound = Coincidence::new();
        bound.set_2d(1, param);
        Interference::from_parts(
            inters,
            bound,
            Orientation::Forward,
            Orientation::Forward,
            Orientation::Forward,
        )
    }

    /// OCCT ctor + Next traversal (cxx L28-81): the boundary vertices and
    /// the interferences merge in parameter order: 0.0 (boundary),
    /// 0.4, 0.7 (interferences), 1.0 (boundary).
    #[test]
    fn vertex_list_traversal_order_merges_edge_and_interferences() {
        let data = TestData::new();
        let mut tool = EdgeInterferenceTool::new(&data);
        tool.load_edge();
        let il = vec![interference_at(0.4, 7), interference_at(0.7, 8)];
        let mut vl = VertexList::new(&tool, &il);
        // Vertex 1: edge start at 0.0.
        assert!(vl.more());
        assert!(vl.is_boundary());
        assert!(!vl.is_interference());
        assert_eq!(vl.current().parameter(), 0.0);
        assert_eq!(vl.current().index(), 11);
        assert_eq!(vl.orientation(), Orientation::Forward);
        vl.next();
        // Vertex 2: interference at 0.4.
        assert!(vl.more());
        assert!(!vl.is_boundary());
        assert!(vl.is_interference());
        assert_eq!(vl.current().parameter(), 0.4);
        assert_eq!(vl.transition(), Orientation::Forward);
        assert_eq!(vl.boundary_transition(), Orientation::Forward);
        vl.next();
        // Vertex 3: interference at 0.7.
        assert!(vl.more());
        assert!(vl.is_interference());
        assert_eq!(vl.current().parameter(), 0.7);
        vl.next();
        // Vertex 4: edge end at 1.0.
        assert!(vl.more());
        assert!(vl.is_boundary());
        assert_eq!(vl.current().parameter(), 1.0);
        assert_eq!(vl.current().index(), 22);
        assert_eq!(vl.orientation(), Orientation::Reversed);
        vl.next();
        assert!(!vl.more());
        assert!(!vl.is_periodic());
    }
}
