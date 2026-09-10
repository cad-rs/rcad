// OCCT LocOpe_GeneratedShape.hxx L29-52 + LocOpe_GeneratedShape.cxx L17-20 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_GeneratedShape.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_GeneratedShape.cxx
//
// OCCT inheritance chain (LocOpe_GeneratedShape.hxx L29):
//   LocOpe_GeneratedShape : Standard_Transient
//     (abstract base of LocOpe_GluedShape, LocOpe_Generator's operand,
//      LocOpe_Prism / LocOpe_Pipe / LocOpe_Revol / LocOpe_DPrism)
// Rust has no inheritance: the abstract base maps to a trait. The OCCT
// protected data members (myGEdges/myList) cannot live on a Rust trait —
// every implementor carries them as its own fields (see
// loc_ope_glued_shape::LocOpeGluedShape) and returns them by reference from
// the accessors, which keeps the OCCT data flow (the base methods return
// references to the protected lists).
//
// The .cxx only carries the RTTI boilerplate (IMPLEMENT_STANDARD_RTTI) —
// no method bodies to translate.
//
// first consumer: BRepFeat_Form family (3b) — BRepFeat_Form::Perform uses
// the GeneratedShape interface through LocOpe_Generator.

use rcad_kernel::topo_shape::Shape;

/// OCCT LocOpe_GeneratedShape (LocOpe_GeneratedShape.hxx L29-52).
pub trait LocOpeGeneratedShape {
    /// OCCT LocOpe_GeneratedShape::GeneratingEdges() (hxx L33) — returns the
    /// list of generating edges (the base myGEdges protected member).
    fn generating_edges(&mut self) -> &Vec<Shape>;

    /// OCCT LocOpe_GeneratedShape::Generated(const TopoDS_Vertex&) (hxx
    /// L37) — returns the edge created by the vertex <v>. If none, must
    /// return a null shape.
    fn generated_vertex(&mut self, v: &Shape) -> Shape;

    /// OCCT LocOpe_GeneratedShape::Generated(const TopoDS_Edge&) (hxx L41) —
    /// returns the face created by the edge <e>. If none, must return a null
    /// shape. Rust has no overloading: the `_edge` suffix carries the OCCT
    /// parameter-type distinction.
    fn generated_edge(&mut self, e: &Shape) -> Shape;

    /// OCCT LocOpe_GeneratedShape::OrientedFaces() (hxx L45) — returns the
    /// list of correctly oriented generated faces (the base myList protected
    /// member).
    fn oriented_faces(&mut self) -> &Vec<Shape>;
}
