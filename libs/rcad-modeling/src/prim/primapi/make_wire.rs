//! OCCT BRepBuilderAPI_MakeWire (TKTopAlgo) — wire builder.
//!
//! OCCT BRepBuilderAPI_MakeWire collects edges (individually or as a list)
//! and materializes a TopoDS_Wire; IsDone() reports whether the wire was
//! built successfully.  The rcad port keeps the OCCT method surface
//! (new / add / add_all / is_done / wire); the full OCCT merge-and-connect
//! logic (BRepBuilderAPI_MakeWire::Add edge-connection handling, OCC27552)
//! is a later work item — the current implementation appends edges in the
//! given order (the flat pool wire container is order-independent).

use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::BRep;

/// OCCT BRepBuilderAPI_MakeWire — edge collection with IsDone/Wire accessors.
#[derive(Debug, Default)]
pub struct MakeWire {
    edges: Vec<Shape>,
    done: bool,
}

impl MakeWire {
    /// OCCT BRepBuilderAPI_MakeWire() — empty wire builder, NotDone().
    pub fn new() -> Self {
        MakeWire {
            edges: Vec::new(),
            done: false,
        }
    }

    /// OCCT BRepBuilderAPI_MakeWire::Add(const TopoDS_Edge&).
    pub fn add(&mut self, edge: Shape) {
        self.edges.push(edge);
    }

    /// OCCT BRepBuilderAPI_MakeWire::Add(const NCollection_List<TopoDS_Shape>&).
    pub fn add_all(&mut self, edges: &[Shape]) {
        self.edges.extend_from_slice(edges);
    }

    /// OCCT BRepBuilderAPI_MakeWire::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BRepBuilderAPI_MakeWire::Wire() — materialize the wire into
    /// `brep` (the flat BRep pool owns the shapes).
    pub fn wire(&mut self, brep: &mut BRep) -> Shape {
        let w = brep.add_twire(std::mem::take(&mut self.edges));
        self.done = true;
        w
    }
}
