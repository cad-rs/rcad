// OCCT BRepFeat_Gluer.hxx L17-108 + BRepFeat_Gluer.cxx L17-63 +
// BRepFeat_Gluer.lxx L17-90 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Gluer.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Gluer.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Gluer.lxx
//
// OCCT inheritance chain (BRepFeat_Gluer.hxx L49):
//   BRepFeat_Gluer : BRepBuilderAPI_MakeShape
//     BRepBuilderAPI_MakeShape : BRepBuilderAPI_Command
// Rust has no inheritance -> composition + delegation. The
// BRepBuilderAPI_MakeShape base sub-object is carried by the BRepFeatGluer
// fields listed below; the LocOpe_Gluer member is rcad's LocOpeGluer
// (crate::feat::loc_ope_gluer::LocOpeGluer).
//
// Field mapping:
//   BRepBuilderAPI_Command::myDone         -> my_done
//   BRepBuilderAPI_MakeShape::myShape      -> my_shape (None = null)
//   BRepBuilderAPI_MakeShape::myGenerated  -> my_generated
//   BRepFeat_Gluer::myGluer                -> my_gluer (LocOpeGluer)
//
// Architecture differences (referenced from the affected functions):
// 1. LocOpe_Gluer::Perform (LocOpe_Gluer.cxx L156-334) is a deferred
//    translation (needs LocOpe_WiresOnShape / LocOpe_Spliter /
//    LocOpe_Generator / LocOpe::TgtFaces — see the loc_ope_gluer.rs header).
//    The myGluer.Perform() call of Build is carried at its spot as the
//    marked gap; with the deferred Perform the myDone stays false and the
//    OCCT NotDone() branch is the one taken (same guarded-skip model as the
//    BOPAlgo_BOP::Perform step of BRepFeat_MakeCylindricalHole).
// 2. The function static `NCollection_List<TopoDS_Shape> LIM` of Modified
//    (cxx L60) is carried as the my_lim field (the OCCT static is process
//    wide and always empty; the Rust carrier is a per-object empty list,
//    same vehicle as LocOpeGluer::my_null_list).

use crate::feat::loc_ope_gluer::LocOpeGluer;
use rcad_kernel::topo_shape::Shape;

/// OCCT BRepFeat_Gluer — glueing of a shape on a basis shape through local
/// operations (BRepFeat_Gluer.hxx L33-49).
pub struct BRepFeatGluer {
    my_done: bool,           // BRepBuilderAPI_Command::myDone
    my_shape: Option<Shape>, // BRepBuilderAPI_MakeShape::myShape (None = null)
    // myGenerated is the MakeShape base member; BRepFeat_Gluer never fills it
    // (the OCCT Modified override reads myGluer.DescendantFaces instead).
    #[allow(dead_code)]
    my_generated: Vec<Shape>,
    my_gluer: LocOpeGluer,   // BRepFeat_Gluer::myGluer
    my_lim: Vec<Shape>,      // OCCT static LIM of Modified (arch. diff. #2)
}

impl Default for BRepFeatGluer {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFeatGluer {
    /// OCCT BRepBuilderAPI_Command::Done().
    fn done(&mut self) {
        self.my_done = true;
    }

    /// OCCT BRepBuilderAPI_Command::NotDone().
    fn not_done(&mut self) {
        self.my_done = false;
    }

    /// OCCT BRepBuilderAPI_Command::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepBuilderAPI_MakeShape::Shape() — the result shape; None
    /// carries the OCCT null shape.
    pub fn shape(&self) -> Option<&Shape> {
        self.my_shape.as_ref()
    }

    /// OCCT BRepFeat_Gluer::BRepFeat_Gluer() (lxx L19) — = default.
    pub fn new() -> Self {
        BRepFeatGluer {
            my_done: false,
            my_shape: None,
            my_generated: Vec::new(),
            my_gluer: LocOpeGluer::new(),
            my_lim: Vec::new(),
        }
    }

    /// OCCT BRepFeat_Gluer::BRepFeat_Gluer(Snew, Sbase) (lxx L23-26) —
    /// myGluer(Sbase, Snew) ("Attention a l'inversion").
    pub fn new_with_shapes(the_snew: &Shape, the_sbase: &Shape) -> Self {
        BRepFeatGluer {
            my_done: false,
            my_shape: None,
            my_generated: Vec::new(),
            my_gluer: LocOpeGluer::with_shapes(the_sbase, the_snew),
            my_lim: Vec::new(),
        }
    }

    /// OCCT BRepFeat_Gluer::Init(Snew, Sbase) (lxx L30-33).
    pub fn init(&mut self, the_snew: &Shape, the_sbase: &Shape) {
        self.my_gluer.init(the_sbase, the_snew);
    }

    /// OCCT BRepFeat_Gluer::Bind(Fnew, Fbase) (lxx L37-40) — the face-face
    /// overload.
    pub fn bind_face(&mut self, the_fnew: &Shape, the_fbase: &Shape) {
        self.my_gluer.bind_face(the_fnew, the_fbase);
    }

    /// OCCT BRepFeat_Gluer::Bind(Enew, Ebase) (lxx L44-47) — the edge-edge
    /// overload.
    pub fn bind_edge(&mut self, the_enew: &Shape, the_ebase: &Shape) {
        self.my_gluer.bind_edge(the_enew, the_ebase);
    }

    /// OCCT BRepFeat_Gluer::OpeType() (lxx L51-54).
    pub fn ope_type(&self) -> crate::feat::loc_ope_operation::LocOpeOperation {
        self.my_gluer.ope_type()
    }

    /// OCCT BRepFeat_Gluer::BasisShape() (lxx L58-61).
    pub fn basis_shape(&self) -> &Shape {
        self.my_gluer.basis_shape()
    }

    /// OCCT BRepFeat_Gluer::GluedShape() (lxx L65-68).
    pub fn glued_shape(&self) -> &Shape {
        self.my_gluer.glued_shape()
    }

    /// OCCT BRepFeat_Gluer::Build (cxx L24-36).
    pub fn build(&mut self) {
        // OCCT cxx L26: myGluer.Perform(); — deferred dependency (architecture
        // difference #1); with it myDone stays false below.
        // OCCT cxx L27: if (myGluer.IsDone()).
        if self.my_gluer.is_done() {
            self.done();
            // OCCT cxx L30: myShape = myGluer.ResultingShape().
            self.my_shape = self.my_gluer.resulting_shape().cloned();
        } else {
            // OCCT cxx L33: NotDone().
            self.not_done();
        }
    }

    /// OCCT BRepFeat_Gluer::IsDeleted(F) (cxx L40-43).
    pub fn is_deleted(&self, the_f: &Shape) -> bool {
        self.my_gluer.descendant_faces(the_f).is_empty()
    }

    /// OCCT BRepFeat_Gluer::Modified(F) (cxx L47-62).
    pub fn modified(&mut self, the_f: &Shape) -> &Vec<Shape> {
        if the_f.shape_type() == rcad_kernel::topods::ShapeType::Face {
            // OCCT cxx L51: LS = myGluer.DescendantFaces(Face(F)).
            if !self.my_gluer.descendant_faces(the_f).is_empty() {
                // OCCT cxx L54: if (!LS.First().IsSame(F)).
                let ls = self.my_gluer.descendant_faces(the_f);
                if !shape_is_same(&ls[0], the_f) {
                    return self.my_gluer.descendant_faces(the_f);
                }
            }
        }
        // OCCT cxx L60: static NCollection_List LIM.
        &self.my_lim
    }
}

/// OCCT TopoDS_Shape::IsSame(S) — same TShape + Location (orientation
/// ignored; the TopTools_ShapeMapHasher identity).
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}
