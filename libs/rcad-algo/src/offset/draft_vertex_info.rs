//! OCCT Draft_VertexInfo (Draft_VertexInfo.hxx L30-60 + Draft_VertexInfo.cxx
//! L25-117) — the per-vertex record of Draft_Modification.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_VertexInfo.hxx / Draft_VertexInfo.cxx
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> myEdges + NCollection_List<double>
//!    myParams -> parallel Vec<Shape> / Vec<f64> (insertion order preserved).
//! 2. NCollection_List<TopoDS_Shape>::Iterator myItEd -> a 1-past-the-end
//!    style index into my_edges (my_edges.len() == no more). The OCCT
//!    Add/Parameter loops leave the member iterator at the found element (or
//!    at the end) — the index reproduces the same observable state.

use glam::DVec3;
use rcad_kernel::topo_shape::Shape;

/// OCCT Draft_VertexInfo (Draft_VertexInfo.hxx L30-60).
#[derive(Debug, Clone)]
pub struct DraftVertexInfo {
    my_geom: DVec3,        // OCCT: myGeom
    my_edges: Vec<Shape>,  // OCCT: myEdges (NCollection_List<TopoDS_Shape>)
    my_params: Vec<f64>,   // OCCT: myParams (NCollection_List<double>)
    my_it_ed: usize,       // OCCT: myItEd (list iterator; index into my_edges)
}

impl Default for DraftVertexInfo {
    fn default() -> Self {
        // OCCT Draft_VertexInfo::Draft_VertexInfo() = default (cxx L25).
        DraftVertexInfo {
            my_geom: DVec3::ZERO,
            my_edges: Vec::new(),
            my_params: Vec::new(),
            my_it_ed: 0,
        }
    }
}

impl DraftVertexInfo {
    /// OCCT Draft_VertexInfo::Draft_VertexInfo() = default (cxx L25).
    pub fn new() -> Self {
        DraftVertexInfo {
            my_geom: DVec3::ZERO,
            my_edges: Vec::new(),
            my_params: Vec::new(),
            my_it_ed: 0,
        }
    }

    /// OCCT Draft_VertexInfo::Add(const TopoDS_Edge& E) (cxx L29-43).
    pub fn add(&mut self, the_e: &Shape) {
        // OCCT: for (myItEd.Initialize(myEdges); myItEd.More(); myItEd.Next())
        //         if (E.IsSame(myItEd.Value())) break;
        self.my_it_ed = 0;
        while self.my_it_ed < self.my_edges.len() {
            if the_e.is_same(&self.my_edges[self.my_it_ed]) {
                break;
            }
            self.my_it_ed += 1;
        }
        if self.my_it_ed >= self.my_edges.len() {
            // OCCT: myEdges.Append(E); myParams.Append(RealLast());
            self.my_edges.push(the_e.clone());
            self.my_params.push(f64::MAX);
        }
    }

    /// OCCT Draft_VertexInfo::Geometry() (cxx L47-50).
    pub fn geometry(&self) -> &DVec3 {
        &self.my_geom
    }

    /// OCCT Draft_VertexInfo::ChangeGeometry() (cxx L54-57).
    pub fn change_geometry(&mut self) -> &mut DVec3 {
        &mut self.my_geom
    }

    /// OCCT Draft_VertexInfo::Parameter(const TopoDS_Edge& E) (cxx L61-73) —
    /// the member iterator walks the edges, the local iterator walks the
    /// parameters in parallel.
    pub fn parameter(&mut self, the_e: &Shape) -> f64 {
        let mut itp = 0usize; // OCCT: NCollection_List<double>::Iterator itp(myParams);
        self.my_it_ed = 0;
        while self.my_it_ed < self.my_edges.len() {
            if self.my_edges[self.my_it_ed].is_same(the_e) {
                return self.my_params[itp];
            }
            self.my_it_ed += 1;
            itp += 1;
        }
        // OCCT: throw Standard_DomainError();
        panic!("Standard_DomainError");
    }

    /// OCCT Draft_VertexInfo::ChangeParameter(const TopoDS_Edge& E)
    /// (cxx L77-89) — the mutable parameter reference is the &mut carrier of
    /// the OCCT `double&` return.
    pub fn change_parameter(&mut self, the_e: &Shape) -> &mut f64 {
        let mut itp = 0usize;
        self.my_it_ed = 0;
        while self.my_it_ed < self.my_edges.len() {
            if self.my_edges[self.my_it_ed].is_same(the_e) {
                return &mut self.my_params[itp];
            }
            self.my_it_ed += 1;
            itp += 1;
        }
        // OCCT: throw Standard_DomainError();
        panic!("Standard_DomainError");
    }

    /// OCCT Draft_VertexInfo::InitEdgeIterator() (cxx L93-96).
    pub fn init_edge_iterator(&mut self) {
        self.my_it_ed = 0;
    }

    /// OCCT Draft_VertexInfo::Edge() (cxx L100-103) — TopoDS::Edge cast is
    /// the identity on rcad Shape handles.
    pub fn edge(&self) -> Shape {
        self.my_edges[self.my_it_ed].clone()
    }

    /// OCCT Draft_VertexInfo::MoreEdge() (cxx L107-110).
    pub fn more_edge(&self) -> bool {
        self.my_it_ed < self.my_edges.len()
    }

    /// OCCT Draft_VertexInfo::NextEdge() (cxx L114-117).
    pub fn next_edge(&mut self) {
        self.my_it_ed += 1;
    }
}
