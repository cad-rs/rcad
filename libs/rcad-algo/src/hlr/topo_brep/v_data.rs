// OCCT HLRTopoBRep_VData (TKHLR/HLRTopoBRep/HLRTopoBRep_VData.hxx L26-41
// + .cxx L26-31 + .lxx L19-38) — the (parameter, vertex) pair stored in the
// per-edge vertex lists of HLRTopoBRep_Data.

use rcad_kernel::topods::Shape;

/// OCCT HLRTopoBRep_VData.
#[derive(Debug, Clone)]
pub struct VData {
    /// OCCT myParameter (hxx L40).
    my_parameter: f64,
    /// OCCT myVertex (hxx L41).
    my_vertex: Shape,
}

impl VData {
    /// OCCT HLRTopoBRep_VData() (lxx L19-23) — myParameter(0.0); the vertex
    /// handle stays default (null).
    pub fn new() -> Self {
        VData {
            my_parameter: 0.0,
            my_vertex: Shape::null(),
        }
    }

    /// OCCT HLRTopoBRep_VData(const double P, const TopoDS_Shape& V) (cxx
    /// L26-31) — myParameter(P), myVertex(V).
    pub fn new_with(p: f64, v: Shape) -> Self {
        VData {
            my_parameter: p,
            my_vertex: v,
        }
    }

    /// OCCT Parameter() (lxx L26-30).
    pub fn parameter(&self) -> f64 {
        self.my_parameter
    }

    /// OCCT Vertex() (lxx L33-37).
    pub fn vertex(&self) -> &Shape {
        &self.my_vertex
    }
}

impl Default for VData {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use rcad_kernel::topods::TShape;

    /// The i-th wire edge of the test face (the same extraction the
    /// topol_tool_brep tests use).
    pub(crate) fn face_outer_edge(face: &Shape, i: usize) -> Shape {
        match &*face.data {
            TShape::Face(fd) => match &*fd.outer_wire.data {
                TShape::Wire(wd) => wd.edges[i].clone(),
                _ => panic!("not a wire"),
            },
            _ => panic!("not a face"),
        }
    }

    /// OCCT anchor: the (parameter, vertex) round-trip — the default ctor
    /// zeroes myParameter (lxx L19-23) and the (P, V) ctor stores both
    /// members (cxx L26-31).
    #[test]
    fn v_data_round_trip() {
        let empty = VData::new();
        assert_eq!(empty.parameter(), 0.0);
        assert!(empty.vertex().is_null());

        let filled = VData::new_with(0.75, Shape::null());
        assert_eq!(filled.parameter(), 0.75);
        assert!(filled.vertex().is_null());

        let (_brep, face) =
            crate::topalgo::brep_top_adaptor::topol_tool_brep::tests::square_face();
        let edge = face_outer_edge(&face, 0);
        let filled = VData::new_with(-1.25, edge.clone());
        assert_eq!(filled.parameter(), -1.25);
        assert!(filled.vertex().is_same(&edge));
    }
}
