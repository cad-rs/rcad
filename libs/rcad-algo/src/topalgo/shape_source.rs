// Shape source abstraction for the face/solid classifiers.
//
// OCCT: BRepClass / BRepClass3d (TKTopAlgo) classify points against
// TopoDS_Shape using the BRepAdaptor_* adaptors (TKBRep). They never touch
// BOPDS (TKBO). rcad's boolean pipeline stores its shapes in the DS (BOPDS)
// at flat indices, so the classifier code needs a way to read shape data
// without depending on the DS type.
//
// This trait exposes exactly the BOPDS-style shape access the classifiers use
// (NbShapes, ShapeInfo::Shape/ShapeType/SubShapes, ShapeIndex, vertex-edge
// map, BRep_Tool-style Surface/Tolerance). The bop::ds::DS implements it, so
// the dependency direction matches OCCT: topalgo (TKTopAlgo) depends only on
// rcad-kernel (TKBRep), and bop (TKBO) depends on topalgo.

use rcad_kernel::geom::Surface3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;

/// BOPDS-style shape data source for the classifiers.
pub trait ShapeSource: Send + Sync {
    /// BOPDS_DS::NbShapes — total number of shapes in the DS.
    fn nb_shapes(&self) -> usize;
    /// BOPDS_ShapeInfo::Shape at a flat DS index.
    fn shape_at(&self, i: usize) -> Shape;
    /// BOPDS_ShapeInfo::ShapeType at a flat DS index.
    fn shape_type(&self, i: usize) -> ShapeType;
    /// BOPDS_ShapeInfo::SubShapes at a flat DS index (child shape indices).
    fn sub_shapes(&self, i: usize) -> &[usize];
    /// BOPDS_DS::ShapeIndex — map a (ptr_id, location) to its flat index.
    fn map_shape_index(&self, ptr_id: u64, location: u32) -> Option<usize>;
    /// BOPDS_DS vertex → incident edges map.
    fn map_ve(&self, vertex: usize) -> Option<&Vec<usize>>;
    /// BRep_Tool::Surface(face) — face surface at a flat DS index.
    fn face_surface(&self, i: usize) -> Option<Surface3>;
    /// BRep_Tool::Tolerance(vertex) — vertex tolerance at a flat DS index.
    fn vertex_tolerance(&self, i: usize) -> f64;
    /// BRep_Tool::Degenerated(edge) — true when the edge at a flat DS index is
    /// degenerated (no 3D curve, coincident vertices).
    fn is_edge_degenerated(&self, i: usize) -> bool;
}
