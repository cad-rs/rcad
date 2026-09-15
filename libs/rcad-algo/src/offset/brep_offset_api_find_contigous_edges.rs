//! OCCT BRepOffsetAPI_FindContigousEdges — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_FindContigousEdges.cxx (L29-316) +
//!         BRepOffsetAPI_FindContigousEdges.hxx (L359-452).
//!
//! OCCT inheritance chain (hxx L359): none — BRepOffsetAPI_FindContigousEdges
//! is a standalone class holding a BRepBuilderAPI_Sewing handle.
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member BRepBuilderAPI_Sewing — the real 1:1 body lives in
//!    crate::topalgo::brep_builder_api_sewing; the facade holds it by value
//!    plus the owning `BRep` pool the sewing methods consume (the
//!    brep_tools_quilt.rs pool convention).  The local GAP carrier is
//!    deleted.
//! 3. Standard_OutOfRange_Raise_if / Standard_NoSuchObject_Raise_if map to
//!    the panic! guards at the OCCT call sites.

use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;

use crate::topalgo::brep_builder_api_sewing::BRepBuilderAPISewing;

/// OCCT BRepOffsetAPI_FindContigousEdges (hxx L359-452).
pub struct BRepOffsetAPIFindContigousEdges {
    // OCCT private member (hxx L448): Handle(BRepBuilderAPI_Sewing) mySewing.
    my_sewing: BRepBuilderAPISewing, // OCCT: mySewing
    // rcad arena stand-in (the sewing pool convention): the OCCT facade
    // rides the global shape space, the rcad one needs the owning pool.
    my_brep: BRep,
}

impl BRepOffsetAPIFindContigousEdges {
    /// OCCT BRepOffsetAPI_FindContigousEdges::BRepOffsetAPI_FindContigousEdges
    /// (tolerance, option) (cxx L29-34; defaults 1.0e-06 / true).
    pub fn new(tolerance: f64, option: bool) -> Self {
        // OCCT L32: mySewing = new BRepBuilderAPI_Sewing;
        let mut r = BRepOffsetAPIFindContigousEdges {
            my_sewing: BRepBuilderAPISewing::new(1.0e-06),
            my_brep: BRep::new(),
        };
        r.init(tolerance, option);
        r
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Init(tolerance, option)
    /// (cxx L38-41).
    pub fn init(&mut self, tolerance: f64, option: bool) {
        // OCCT L40: mySewing->Init(tolerance, option, false, true).
        self.my_sewing
            .init_full(&mut self.my_brep, tolerance, option, false, true, false);
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Add(aShape) (cxx L45-48).
    pub fn add(&mut self, a_shape: &Shape) {
        self.my_sewing.add(&mut self.my_brep, a_shape);
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Perform() (cxx L52-55).
    pub fn perform(&mut self) {
        self.my_sewing.perform(&mut self.my_brep);
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::NbContigousEdges()
    /// (cxx L59-62).
    pub fn nb_contigous_edges(&self) -> i32 {
        self.my_sewing.nb_contigous_edges()
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::ContigousEdge(index)
    /// (cxx L66-71).
    pub fn contigous_edge(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_contigous_edges()),
            "BRepOffsetAPI_FindContigousEdges::ContigousEdge"
        );
        self.my_sewing.contigous_edge(index)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::ContigousEdgeCouple(index)
    /// (cxx L75-81).
    pub fn contigous_edge_couple(&self, index: i32) -> Vec<Shape> {
        assert!(
            !(index < 0 || index > self.nb_contigous_edges()),
            "BRepOffsetAPI_FindContigousEdges::ContigousEdgeCouple"
        );
        self.my_sewing.contigous_edge_couple(index)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::SectionToBoundary(section)
    /// (cxx L85-91).
    pub fn section_to_boundary(&self, section: &Shape) -> Shape {
        assert!(
            self.my_sewing.is_section_bound(section),
            "BRepOffsetAPI_FindContigousEdges::SectionToBoundary"
        );
        self.my_sewing.section_to_boundary(section)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::NbDegeneratedShapes()
    /// (cxx L95-98).
    pub fn nb_degenerated_shapes(&self) -> i32 {
        self.my_sewing.nb_degenerated_shapes()
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::DegeneratedShape(index)
    /// (cxx L102-107).
    pub fn degenerated_shape(&self, index: i32) -> Shape {
        assert!(
            !(index < 0 || index > self.nb_degenerated_shapes()),
            "BRepOffsetAPI_FindContigousEdges::DegereratedShape"
        );
        self.my_sewing.degenerated_shape(index)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::IsDegenerated(aShape)
    /// (cxx L111-114).
    pub fn is_degenerated(&mut self, a_shape: &Shape) -> bool {
        self.my_sewing.is_degenerated(&mut self.my_brep, a_shape)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::IsModified(aShape)
    /// (cxx L118-121).
    pub fn is_modified(&self, a_shape: &Shape) -> bool {
        self.my_sewing.is_modified(a_shape)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Modified(aShape)
    /// (cxx L125-129).
    pub fn modified(&self, a_shape: &Shape) -> Shape {
        assert!(
            self.is_modified(a_shape),
            "BRepOffsetAPI_FindContigousEdges::Modified"
        );
        self.my_sewing.modified(a_shape)
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Dump() (cxx L133-136).
    pub fn dump(&mut self) {
        self.my_sewing.dump(&mut self.my_brep);
    }
}
