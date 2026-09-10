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
//! 2. The engine member BRepBuilderAPI_Sewing (TKShHealing/BRepBuilderAPI)
//!    belongs to the shhealing domain, excluded from this package batch —
//!    the BRepBuilderAPISewing carrier below keeps the OCCT method surface
//!    with GAP panics (port plan section 0.6); every facade body around the
//!    engine calls is translated 1:1.
//! 3. Standard_OutOfRange_Raise_if / Standard_NoSuchObject_Raise_if map to
//!    the panic! guards at the OCCT call sites.

use rcad_kernel::topo_shape::Shape;

// ---------------------------------------------------------------------------
// GAP carrier (architecture difference #2).
// ---------------------------------------------------------------------------

/// OCCT BRepBuilderAPI_Sewing (TKShHealing/BRepBuilderAPI_Sewing.hxx) — the
/// sewing engine of FindContigousEdges (architecture difference #2; GAP: the
/// shhealing-domain sewing is excluded from this batch — the GAP panics are
/// the section 0.6 annotation; the method surface keeps the OCCT form).
pub struct BRepBuilderAPISewing;

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::Init(tolerance, option, sewing, analysis)
    /// — GAP.
    pub fn init(&mut self, _tolerance: f64, _option: bool, _sewing: bool, _analysis: bool) {
        panic!("GAP: BRepBuilderAPI_Sewing::Init (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::Add(shape) — GAP.
    pub fn add(&mut self, _shape: &Shape) {
        panic!("GAP: BRepBuilderAPI_Sewing::Add (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::Perform() — GAP.
    pub fn perform(&mut self) {
        panic!("GAP: BRepBuilderAPI_Sewing::Perform (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::NbContigousEdges() — GAP.
    pub fn nb_contigous_edges(&self) -> i32 {
        panic!("GAP: BRepBuilderAPI_Sewing::NbContigousEdges (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::ContigousEdge(index) — GAP.
    pub fn contigous_edge(&self, _index: i32) -> Shape {
        panic!("GAP: BRepBuilderAPI_Sewing::ContigousEdge (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::ContigousEdgeCouple(index) — GAP.
    pub fn contigous_edge_couple(&self, _index: i32) -> Vec<Shape> {
        panic!(
            "GAP: BRepBuilderAPI_Sewing::ContigousEdgeCouple (TKShHealing excluded from this batch)"
        );
    }

    /// OCCT BRepBuilderAPI_Sewing::IsSectionBound(section) — GAP.
    pub fn is_section_bound(&self, _section: &Shape) -> bool {
        panic!("GAP: BRepBuilderAPI_Sewing::IsSectionBound (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::SectionToBoundary(section) — GAP.
    pub fn section_to_boundary(&self, _section: &Shape) -> Shape {
        panic!(
            "GAP: BRepBuilderAPI_Sewing::SectionToBoundary (TKShHealing excluded from this batch)"
        );
    }

    /// OCCT BRepBuilderAPI_Sewing::NbDegeneratedShapes() — GAP.
    pub fn nb_degenerated_shapes(&self) -> i32 {
        panic!(
            "GAP: BRepBuilderAPI_Sewing::NbDegeneratedShapes (TKShHealing excluded from this batch)"
        );
    }

    /// OCCT BRepBuilderAPI_Sewing::DegeneratedShape(index) — GAP.
    pub fn degenerated_shape(&self, _index: i32) -> Shape {
        panic!(
            "GAP: BRepBuilderAPI_Sewing::DegeneratedShape (TKShHealing excluded from this batch)"
        );
    }

    /// OCCT BRepBuilderAPI_Sewing::IsDegenerated(shape) — GAP.
    pub fn is_degenerated(&self, _shape: &Shape) -> bool {
        panic!("GAP: BRepBuilderAPI_Sewing::IsDegenerated (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::IsModified(shape) — GAP.
    pub fn is_modified(&self, _shape: &Shape) -> bool {
        panic!("GAP: BRepBuilderAPI_Sewing::IsModified (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::Modified(shape) — GAP.
    pub fn modified(&self, _shape: &Shape) -> Shape {
        panic!("GAP: BRepBuilderAPI_Sewing::Modified (TKShHealing excluded from this batch)");
    }

    /// OCCT BRepBuilderAPI_Sewing::Dump() — GAP.
    pub fn dump(&self) {
        panic!("GAP: BRepBuilderAPI_Sewing::Dump (TKShHealing excluded from this batch)");
    }
}

/// OCCT BRepOffsetAPI_FindContigousEdges (hxx L359-452).
pub struct BRepOffsetAPIFindContigousEdges {
    // OCCT private member (hxx L448): Handle(BRepBuilderAPI_Sewing) mySewing.
    my_sewing: BRepBuilderAPISewing, // OCCT: mySewing
}

impl BRepOffsetAPIFindContigousEdges {
    /// OCCT BRepOffsetAPI_FindContigousEdges::BRepOffsetAPI_FindContigousEdges
    /// (tolerance, option) (cxx L29-34; defaults 1.0e-06 / true).
    pub fn new(tolerance: f64, option: bool) -> Self {
        // OCCT L32: mySewing = new BRepBuilderAPI_Sewing;
        let mut r = BRepOffsetAPIFindContigousEdges {
            my_sewing: BRepBuilderAPISewing,
        };
        r.init(tolerance, option);
        r
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Init(tolerance, option)
    /// (cxx L38-41).
    pub fn init(&mut self, tolerance: f64, option: bool) {
        // OCCT L40: mySewing->Init(tolerance, option, false, true);
        self.my_sewing.init(tolerance, option, false, true);
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Add(aShape) (cxx L45-48).
    pub fn add(&mut self, a_shape: &Shape) {
        self.my_sewing.add(a_shape);
    }

    /// OCCT BRepOffsetAPI_FindContigousEdges::Perform() (cxx L52-55).
    pub fn perform(&mut self) {
        self.my_sewing.perform();
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
    pub fn is_degenerated(&self, a_shape: &Shape) -> bool {
        self.my_sewing.is_degenerated(a_shape)
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
    pub fn dump(&self) {
        self.my_sewing.dump();
    }
}
