//! OCCT Draft_EdgeInfo (Draft_EdgeInfo.hxx L29-83 + Draft_EdgeInfo.cxx
//! L26-173) — the per-edge record of Draft_Modification.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_EdgeInfo.hxx / Draft_EdgeInfo.cxx

use glam::DVec3;
use rcad_kernel::geom::{Curve2d, Curve3};
use rcad_kernel::topo_shape::Shape;

use super::draft_modification::brep_tool_tolerance;

/// OCCT Draft_EdgeInfo (Draft_EdgeInfo.hxx L29-83).
#[derive(Debug, Clone)]
pub struct DraftEdgeInfo {
    my_new_geom: bool,          // OCCT: myNewGeom
    my_geom: Option<Curve3>,    // OCCT: myGeom (Handle(Geom_Curve))
    my_first_f: Shape,          // OCCT: myFirstF
    my_secon_f: Shape,          // OCCT: mySeconF
    my_first_pc: Option<Curve2d>, // OCCT: myFirstPC (Handle(Geom2d_Curve))
    my_secon_pc: Option<Curve2d>, // OCCT: mySeconPC (Handle(Geom2d_Curve))
    my_root_face: Shape,        // OCCT: myRootFace
    my_tgt: bool,               // OCCT: myTgt
    my_pt: DVec3,               // OCCT: myPt
    my_tol: f64,                // OCCT: myTol
}

impl Default for DraftEdgeInfo {
    fn default() -> Self {
        // OCCT Draft_EdgeInfo::Draft_EdgeInfo() (cxx L26-31).
        DraftEdgeInfo {
            my_new_geom: false,
            my_geom: None,
            my_first_f: Shape::null(),
            my_secon_f: Shape::null(),
            my_first_pc: None,
            my_secon_pc: None,
            my_root_face: Shape::null(),
            my_tgt: false,
            my_pt: DVec3::ZERO,
            my_tol: 0.0,
        }
    }
}

impl DraftEdgeInfo {
    /// OCCT Draft_EdgeInfo::Draft_EdgeInfo() (cxx L26-31).
    pub fn new() -> Self {
        DraftEdgeInfo {
            my_new_geom: false,
            my_geom: None,
            my_first_f: Shape::null(),
            my_secon_f: Shape::null(),
            my_first_pc: None,
            my_secon_pc: None,
            my_root_face: Shape::null(),
            my_tgt: false,
            my_pt: DVec3::ZERO,
            my_tol: 0.0,
        }
    }

    /// OCCT Draft_EdgeInfo::Draft_EdgeInfo(const bool HasNewGeometry)
    /// (cxx L35-40).
    pub fn new_with_geometry(has_new_geometry: bool) -> Self {
        DraftEdgeInfo {
            my_new_geom: has_new_geometry,
            my_geom: None,
            my_first_f: Shape::null(),
            my_secon_f: Shape::null(),
            my_first_pc: None,
            my_secon_pc: None,
            my_root_face: Shape::null(),
            my_tgt: false,
            my_pt: DVec3::ZERO,
            my_tol: 0.0,
        }
    }

    /// OCCT Draft_EdgeInfo::Add(const TopoDS_Face& F) (cxx L44-55).
    pub fn add(&mut self, the_f: &Shape) {
        if self.my_first_f.is_null() {
            self.my_first_f = the_f.clone();
        } else if !self.my_first_f.is_same(the_f) && self.my_secon_f.is_null() {
            self.my_secon_f = the_f.clone();
        }
        // OCCT: myTol = std::max(myTol, BRep_Tool::Tolerance(F));
        self.my_tol = self.my_tol.max(brep_tool_tolerance(the_f));
    }

    /// OCCT Draft_EdgeInfo::RootFace(const TopoDS_Face& F) (cxx L59-62).
    pub fn set_root_face(&mut self, the_f: &Shape) {
        self.my_root_face = the_f.clone();
    }

    /// OCCT Draft_EdgeInfo::Tangent(const gp_Pnt& P) (cxx L66-70).
    pub fn tangent(&mut self, the_p: DVec3) {
        self.my_tgt = true;
        self.my_pt = the_p;
    }

    /// OCCT Draft_EdgeInfo::IsTangent(gp_Pnt& P) (cxx L74-78) — the output
    /// parameter is carried as &mut.
    pub fn is_tangent(&self, the_p: &mut DVec3) -> bool {
        *the_p = self.my_pt;
        self.my_tgt
    }

    /// OCCT Draft_EdgeInfo::NewGeometry() (cxx L82-85).
    pub fn new_geometry(&self) -> bool {
        self.my_new_geom
    }

    /// OCCT Draft_EdgeInfo::SetNewGeometry(const bool NewGeom) (cxx L89-92).
    pub fn set_new_geometry(&mut self, new_geom: bool) {
        self.my_new_geom = new_geom;
    }

    /// OCCT Draft_EdgeInfo::Geometry() (cxx L96-99).
    pub fn geometry(&self) -> Option<&Curve3> {
        self.my_geom.as_ref()
    }

    /// OCCT Draft_EdgeInfo::FirstFace() (cxx L103-106).
    pub fn first_face(&self) -> &Shape {
        &self.my_first_f
    }

    /// OCCT Draft_EdgeInfo::SecondFace() (cxx L110-113).
    pub fn second_face(&self) -> &Shape {
        &self.my_secon_f
    }

    /// OCCT Draft_EdgeInfo::FirstPC() (cxx L127-130).
    pub fn first_pc(&self) -> Option<&Curve2d> {
        self.my_first_pc.as_ref()
    }

    /// OCCT Draft_EdgeInfo::SecondPC() (cxx L137-140).
    pub fn second_pc(&self) -> Option<&Curve2d> {
        self.my_secon_pc.as_ref()
    }

    /// OCCT Draft_EdgeInfo::ChangeGeometry() (cxx L117-120).
    pub fn change_geometry(&mut self) -> &mut Option<Curve3> {
        &mut self.my_geom
    }

    /// OCCT Draft_EdgeInfo::ChangeFirstPC() (cxx L144-147).
    pub fn change_first_pc(&mut self) -> &mut Option<Curve2d> {
        &mut self.my_first_pc
    }

    /// OCCT Draft_EdgeInfo::ChangeSecondPC() (cxx L151-154).
    pub fn change_second_pc(&mut self) -> &mut Option<Curve2d> {
        &mut self.my_secon_pc
    }

    /// OCCT Draft_EdgeInfo::RootFace() (cxx L158-161).
    pub fn root_face(&self) -> &Shape {
        &self.my_root_face
    }

    /// OCCT Draft_EdgeInfo::Tolerance(const double tol) (cxx L165-168) — the
    /// setter is named set_tolerance (no Rust overloading).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.my_tol = tol;
    }

    /// OCCT Draft_EdgeInfo::Tolerance() (cxx L170-173).
    pub fn tolerance(&self) -> f64 {
        self.my_tol
    }
}
