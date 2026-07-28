//! Strongly-typed entity and reference IDs.
//!
//! OCCT BRepGraphInc: BRepGraphInc_RepId.hxx

use std::fmt;

macro_rules! typed_id {
    ($name:ident, $doc:literal) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[doc = concat!("Strongly-typed ID for ", $doc, ".")]
        pub struct $name(pub u32);

        impl $name {
            pub const INVALID: Self = $name(u32::MAX);
            pub fn is_valid(self) -> bool { self.0 != u32::MAX }
            pub fn index(self) -> usize { self.0 as usize }
        }
        impl From<u32> for $name { fn from(v: u32) -> Self { $name(v) } }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}({})", stringify!($name), self.0)
            }
        }
    };
}

typed_id!(VertexId, "vertex entities");
typed_id!(EdgeId, "edge entities");
typed_id!(CoEdgeId, "co-edge entities");
typed_id!(WireId, "wire entities");
typed_id!(FaceId, "face entities");
typed_id!(ShellId, "shell entities");
typed_id!(SolidId, "solid entities");
typed_id!(CompoundId, "compound entities");
typed_id!(CompSolidId, "comp-solid entities");
typed_id!(ProductId, "product (assembly) entities");
typed_id!(OccurrenceId, "occurrence (assembly placement) entities");

typed_id!(ShellRefId, "shell reference entries");
typed_id!(FaceRefId, "face reference entries");
typed_id!(WireRefId, "wire reference entries");
typed_id!(VertexRefId, "vertex reference entries");
typed_id!(SolidRefId, "solid reference entries");
typed_id!(ChildRefId, "child (compound) reference entries");
typed_id!(OccurrenceRefId, "occurrence reference entries");

typed_id!(SurfaceRepId, "face surface representation");
typed_id!(Curve3DRepId, "edge 3D curve representation");
typed_id!(Curve2DRepId, "co-edge 2D curve (pcurve) representation");
typed_id!(TriangulationRepId, "face triangulation representation");
typed_id!(Polygon3DRepId, "edge 3D polygon representation");
typed_id!(Polygon2DRepId, "co-edge 2D polygon representation");
typed_id!(PolygonOnTriRepId, "co-edge polygon-on-triangulation representation");
