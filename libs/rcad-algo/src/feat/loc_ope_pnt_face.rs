// OCCT LocOpe_PntFace.hxx L56-105 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_PntFace.hxx
//
// LocOpe_PntFace is a header-only value class (no .cxx). It carries a point
// on a face together with the parameter along the intersecting line and the
// surface U,V parameters at the point.
//
// Null-ability mapping: the empty constructor leaves myFace a NULL
// TopoDS_Face (gp_Pnt/gp scalars default to zero); rcad has no null Shape,
// so Option<Shape> carries the OCCT null face (the same mapping the 3a
// skeleton used in brep_feat_make_cylindrical_hole.rs).
//
// first consumer: BRepFeat_MakeCylindricalHole (3a) — GetOffset/Perform read
// Pnt/Face/Parameter/UParameter/VParameter off the intersector results;
// BRepFeat_Form family (3b) through the LocOpe_CSIntersector.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::Orientation;

/// OCCT LocOpe_PntFace (LocOpe_PntFace.hxx L56-105).
#[derive(Debug, Clone)]
pub struct LocOpePntFace {
    my_pnt: Option<glam::DVec3>, // OCCT: myPnt (gp_Pnt)
    my_face: Option<Shape>,      // OCCT: myFace (TopoDS_Face)
    my_ori: Orientation,         // OCCT: myOri (TopAbs_Orientation)
    my_par: f64,                 // OCCT: myPar
    my_u_par: f64,               // OCCT: myUPar
    my_v_par: f64,               // OCCT: myVPar
}

impl Default for LocOpePntFace {
    fn default() -> Self {
        LocOpePntFace {
            my_pnt: None,
            my_face: None,
            my_ori: Orientation::Forward,
            my_par: 0.0,
            my_u_par: 0.0,
            my_v_par: 0.0,
        }
    }
}

impl LocOpePntFace {
    /// OCCT LocOpe_PntFace::LocOpe_PntFace() (hxx L62-67) — empty constructor.
    /// Useful only for the list.
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT LocOpe_PntFace::LocOpe_PntFace(P, F, Or, Param, UPar, VPar)
    /// (hxx L69-82) — Rust has no overloading: the value constructor carries
    /// the `_full` suffix.
    pub fn new_full(
        the_p: glam::DVec3,
        the_f: &Shape,
        the_or: Orientation,
        the_param: f64,
        the_upar: f64,
        the_vpar: f64,
    ) -> Self {
        LocOpePntFace {
            my_pnt: Some(the_p),
            my_face: Some(the_f.clone()),
            my_ori: the_or,
            my_par: the_param,
            my_u_par: the_upar,
            my_v_par: the_vpar,
        }
    }

    /// OCCT LocOpe_PntFace::Pnt() (hxx L84).
    pub fn pnt(&self) -> glam::DVec3 {
        self.my_pnt.unwrap_or(glam::DVec3::ZERO)
    }

    /// OCCT LocOpe_PntFace::Face() (hxx L86) — None carries the OCCT null face.
    pub fn face(&self) -> Option<&Shape> {
        self.my_face.as_ref()
    }

    /// OCCT LocOpe_PntFace::Orientation() (hxx L88).
    pub fn orientation(&self) -> Orientation {
        self.my_ori
    }

    /// OCCT LocOpe_PntFace::ChangeOrientation() (hxx L90) — the mutable
    /// reference form keeps the OCCT "in-place change" shape.
    pub fn change_orientation(&mut self) -> &mut Orientation {
        &mut self.my_ori
    }

    /// OCCT LocOpe_PntFace::Parameter() (hxx L92).
    pub fn parameter(&self) -> f64 {
        self.my_par
    }

    /// OCCT LocOpe_PntFace::UParameter() (hxx L94).
    pub fn u_parameter(&self) -> f64 {
        self.my_u_par
    }

    /// OCCT LocOpe_PntFace::VParameter() (hxx L96).
    pub fn v_parameter(&self) -> f64 {
        self.my_v_par
    }
}
