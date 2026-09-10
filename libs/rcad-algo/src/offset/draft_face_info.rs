//! OCCT Draft_FaceInfo (Draft_FaceInfo.hxx L29-65 + Draft_FaceInfo.cxx
//! L24-117) — the per-face record of Draft_Modification.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_FaceInfo.hxx / Draft_FaceInfo.cxx

use rcad_kernel::geom::{Curve3, Surface3};
use rcad_kernel::topo_shape::Shape;

/// OCCT Draft_FaceInfo (Draft_FaceInfo.hxx L29-65).
#[derive(Debug, Clone)]
pub struct DraftFaceInfo {
    my_new_geom: bool,       // OCCT: myNewGeom
    my_geom: Option<Surface3>, // OCCT: myGeom (Handle(Geom_Surface))
    my_root_face: Shape,     // OCCT: myRootFace
    my_f1: Shape,            // OCCT: myF1
    my_f2: Shape,            // OCCT: myF2
    my_curv: Option<Curve3>, // OCCT: myCurv (Handle(Geom_Curve))
}

impl Default for DraftFaceInfo {
    fn default() -> Self {
        // OCCT Draft_FaceInfo::Draft_FaceInfo() = default (cxx L24).
        DraftFaceInfo {
            my_new_geom: false,
            my_geom: None,
            my_root_face: Shape::null(),
            my_f1: Shape::null(),
            my_f2: Shape::null(),
            my_curv: None,
        }
    }
}

impl DraftFaceInfo {
    /// OCCT Draft_FaceInfo::Draft_FaceInfo() = default (cxx L24).
    pub fn new() -> Self {
        DraftFaceInfo {
            my_new_geom: false,
            my_geom: None,
            my_root_face: Shape::null(),
            my_f1: Shape::null(),
            my_f2: Shape::null(),
            my_curv: None,
        }
    }

    /// OCCT Draft_FaceInfo::Draft_FaceInfo(S, HasNewGeometry) (cxx L28-40) —
    /// the rectangular trimmed surface is replaced by its basis; the OCCT
    /// null-handle constructor (S = null) is the None carrier.
    pub fn new_with_geometry(the_s: Option<&Surface3>, has_new_geometry: bool) -> Self {
        let mut info = DraftFaceInfo {
            my_new_geom: has_new_geometry,
            my_geom: None,
            my_root_face: Shape::null(),
            my_f1: Shape::null(),
            my_f2: Shape::null(),
            my_curv: None,
        };
        // OCCT cxx L31-39: T = down_cast<Geom_RectangularTrimmedSurface>(S);
        // if (!T.IsNull()) myGeom = T->BasisSurface(); else myGeom = S;
        match the_s {
            Some(Surface3::Trimmed(t)) => {
                info.my_geom = Some(t.basis.as_ref().clone());
            }
            Some(s) => {
                info.my_geom = Some(s.clone());
            }
            None => {}
        }
        info
    }

    /// OCCT Draft_FaceInfo::RootFace(const TopoDS_Face& F) (cxx L44-47).
    pub fn set_root_face(&mut self, the_f: &Shape) {
        self.my_root_face = the_f.clone();
    }

    /// OCCT Draft_FaceInfo::Add(const TopoDS_Face& F) (cxx L51-61).
    pub fn add(&mut self, the_f: &Shape) {
        if self.my_f1.is_null() {
            self.my_f1 = the_f.clone();
        } else if self.my_f2.is_null() {
            self.my_f2 = the_f.clone();
        }
    }

    /// OCCT Draft_FaceInfo::FirstFace() (cxx L65-68).
    pub fn first_face(&self) -> &Shape {
        &self.my_f1
    }

    /// OCCT Draft_FaceInfo::SecondFace() (cxx L72-75).
    pub fn second_face(&self) -> &Shape {
        &self.my_f2
    }

    /// OCCT Draft_FaceInfo::NewGeometry() (cxx L79-82).
    pub fn new_geometry(&self) -> bool {
        self.my_new_geom
    }

    /// OCCT Draft_FaceInfo::Geometry() (cxx L86-89).
    pub fn geometry(&self) -> Option<&Surface3> {
        self.my_geom.as_ref()
    }

    /// OCCT Draft_FaceInfo::ChangeGeometry() (cxx L93-96).
    pub fn change_geometry(&mut self) -> &mut Option<Surface3> {
        &mut self.my_geom
    }

    /// OCCT Draft_FaceInfo::Curve() (cxx L100-103).
    pub fn curve(&self) -> Option<&Curve3> {
        self.my_curv.as_ref()
    }

    /// OCCT Draft_FaceInfo::ChangeCurve() (cxx L107-110).
    pub fn change_curve(&mut self) -> &mut Option<Curve3> {
        &mut self.my_curv
    }

    /// OCCT Draft_FaceInfo::RootFace() (cxx L114-117).
    pub fn root_face(&self) -> &Shape {
        &self.my_root_face
    }
}
