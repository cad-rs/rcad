use std::sync::Arc;
use rcad_kernel::topods::{self, BRep, ShapeRef, TShape, TVertexData, TEdgeData};

// ---- Helper views providing old-style BRep field access ----

pub(crate) struct SolidsView<'a> { brep: &'a BRep }
pub(crate) struct EdgesView<'a> { brep: &'a BRep }
pub(crate) struct VerticesView<'a> { brep: &'a BRep }

impl<'a> SolidsView<'a> {
    fn indices(&self) -> Vec<ShapeRef> {
        self.brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(&**ts, TShape::Solid(_)))
            .map(|(i, _)| ShapeRef::synthetic(i)).collect()
    }
    pub fn is_empty(&self) -> bool { self.indices().is_empty() }
    pub fn len(&self) -> usize { self.indices().len() }
    pub fn iter(&self) -> Vec<ShapeRef> { self.indices() }
}

impl<'a> EdgesView<'a> {
    fn indices(&self) -> Vec<ShapeRef> {
        self.brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(&**ts, TShape::Edge(_)))
            .map(|(i, _)| ShapeRef::synthetic(i)).collect()
    }
    pub fn is_empty(&self) -> bool { self.indices().is_empty() }
    pub fn len(&self) -> usize { self.indices().len() }
    pub fn iter(&self) -> Vec<ShapeRef> { self.indices() }
    pub fn get(&self, idx: usize) -> Option<ShapeRef> { self.indices().get(idx).copied() }
}

impl<'a> VerticesView<'a> {
    fn indices(&self) -> Vec<ShapeRef> {
        self.brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(&**ts, TShape::Vertex(_)))
            .map(|(i, _)| ShapeRef::synthetic(i)).collect()
    }
    pub fn is_empty(&self) -> bool { self.indices().is_empty() }
    pub fn len(&self) -> usize { self.indices().len() }
    pub fn iter(&self) -> Vec<ShapeRef> { self.indices() }
    pub fn get(&self, idx: usize) -> Option<ShapeRef> { self.indices().get(idx).copied() }
}

/// Extension trait: old-style BRep field access for `topods::BRep`.
pub(crate) trait BRepExt {
    fn solids(&self) -> SolidsView<'_>;
    fn edges(&self) -> EdgesView<'_>;
    fn vertices(&self) -> VerticesView<'_>;
    fn vertex_data(&self, sr: ShapeRef) -> &TVertexData;
    fn edge_data(&self, sr: ShapeRef) -> &TEdgeData;
    fn edge_start(&self, sr: ShapeRef) -> ShapeRef;
    fn edge_end(&self, sr: ShapeRef) -> ShapeRef;
    fn edge_curve(&self, sr: ShapeRef) -> Option<rcad_kernel::geom::Curve3>;
    fn edge_t_range(&self, sr: ShapeRef) -> [f64; 2];
}

impl BRepExt for BRep {
    fn solids(&self) -> SolidsView<'_> { SolidsView { brep: self } }
    fn edges(&self) -> EdgesView<'_> { EdgesView { brep: self } }
    fn vertices(&self) -> VerticesView<'_> { VerticesView { brep: self } }
    fn vertex_data(&self, sr: ShapeRef) -> &TVertexData { self.vertex(sr) }
    fn edge_data(&self, sr: ShapeRef) -> &TEdgeData { self.edge(sr) }
    fn edge_start(&self, sr: ShapeRef) -> ShapeRef { self.edge(sr).first }
    fn edge_end(&self, sr: ShapeRef) -> ShapeRef { self.edge(sr).last }
    fn edge_curve(&self, sr: ShapeRef) -> Option<rcad_kernel::geom::Curve3> { self.edge(sr).curve.clone() }
    fn edge_t_range(&self, sr: ShapeRef) -> [f64; 2] { self.edge(sr).t_range }
}
