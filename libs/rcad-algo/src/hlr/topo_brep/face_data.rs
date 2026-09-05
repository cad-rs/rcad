// OCCT HLRTopoBRep_FaceData (TKHLR/HLRTopoBRep/HLRTopoBRep_FaceData.hxx
// L29-51 + .cxx L21-23 + .lxx L19-58) — the three shape lists of a face
// (Internal OutLines, OutLines on restriction and IsoLines).
//
// OCCT NCollection_List<TopoDS_Shape> maps to Vec<Shape> (the established
// list mapping, e.g. geomalgo/hatch).

use rcad_kernel::topods::Shape;

/// OCCT HLRTopoBRep_FaceData — contains the 3 ListOfShape of a Face.
#[derive(Debug, Clone)]
pub struct FaceData {
    /// OCCT myIntL (hxx L49) — the internal OutLines.
    my_int_l: Vec<Shape>,
    /// OCCT myOutL (hxx L50) — the OutLines on restriction.
    my_out_l: Vec<Shape>,
    /// OCCT myIsoL (hxx L51) — the IsoLines.
    my_iso_l: Vec<Shape>,
}

impl FaceData {
    /// OCCT HLRTopoBRep_FaceData() (cxx L21-23) — the default ctor; the
    /// three lists start empty.
    pub fn new() -> Self {
        FaceData {
            my_int_l: Vec::new(),
            my_out_l: Vec::new(),
            my_iso_l: Vec::new(),
        }
    }

    /// OCCT FaceIntL() (lxx L19-23).
    pub fn face_int_l(&self) -> &Vec<Shape> {
        &self.my_int_l
    }

    /// OCCT FaceOutL() (lxx L26-30).
    pub fn face_out_l(&self) -> &Vec<Shape> {
        &self.my_out_l
    }

    /// OCCT FaceIsoL() (lxx L33-37).
    pub fn face_iso_l(&self) -> &Vec<Shape> {
        &self.my_iso_l
    }

    /// OCCT AddIntL() (lxx L40-44).
    pub fn add_int_l(&mut self) -> &mut Vec<Shape> {
        &mut self.my_int_l
    }

    /// OCCT AddOutL() (lxx L47-51).
    pub fn add_out_l(&mut self) -> &mut Vec<Shape> {
        &mut self.my_out_l
    }

    /// OCCT AddIsoL() (lxx L54-58).
    pub fn add_iso_l(&mut self) -> &mut Vec<Shape> {
        &mut self.my_iso_l
    }
}

impl Default for FaceData {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::topo_brep::v_data::tests::face_outer_edge;

    /// OCCT anchor: the per-face three-list round-trip — the default ctor
    /// leaves the lists empty (cxx L21-23) and the Add/Face accessor pairs
    /// alias the same list (lxx L19-58).
    #[test]
    fn face_data_round_trip() {
        let empty = FaceData::new();
        assert!(empty.face_int_l().is_empty());
        assert!(empty.face_out_l().is_empty());
        assert!(empty.face_iso_l().is_empty());

        let (_brep, face) =
            crate::topalgo::brep_top_adaptor::topol_tool_brep::tests::square_face();
        let e0 = face_outer_edge(&face, 0);
        let e1 = face_outer_edge(&face, 1);

        let mut fd = FaceData::new();
        fd.add_int_l().push(e0.clone());
        fd.add_out_l().push(e1.clone());
        fd.add_out_l().push(e0.clone());
        fd.add_iso_l().push(e1.clone());

        assert_eq!(fd.face_int_l().len(), 1);
        assert_eq!(fd.face_out_l().len(), 2);
        assert_eq!(fd.face_iso_l().len(), 1);
        assert!(fd.face_int_l()[0].is_same(&e0));
        assert!(fd.face_out_l()[0].is_same(&e1));
        assert!(fd.face_out_l()[1].is_same(&e0));
        assert!(fd.face_iso_l()[0].is_same(&e1));
    }
}
