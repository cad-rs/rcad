// OCCT IntPatch_Line — the abstract intersection line of two surfaces.
//
// OCCT IntPatch_Line.hxx / .cxx / .lxx. A line is either geometric (line,
// circle, ellipse, parabola, hyperbola: GLine), analytic (ALine), or a set of
// points from a walking algorithm (WLine), or a restriction arc (RLine).
//
// rcad data-model notes:
// - IntSurf_TypeTrans -> crate::topalgo::int_res2d::TypeTrans.
// - IntSurf_Situation  -> crate::topalgo::int_res2d::Situation.

use crate::topalgo::int_patch::i_type::IntPatchIType;
use crate::topalgo::int_res2d::{Situation, TypeTrans};

/// OCCT IntPatch_Line — base class.
#[derive(Debug, Clone)]
pub struct IntPatchLine {
    typ: IntPatchIType,
    tg: bool,
    t_s1: TypeTrans,
    t_s2: TypeTrans,
    sit1: Situation,
    sit2: Situation,
    u_s1: bool,
    v_s1: bool,
    u_s2: bool,
    v_s2: bool,
}

impl IntPatchLine {
    /// OCCT IntPatch_Line(Tang, Trans1, Trans2) — transitions In or Out.
    pub fn new_with_trans(tang: bool, trans1: TypeTrans, trans2: TypeTrans) -> Self {
        IntPatchLine {
            typ: IntPatchIType::Walking,
            tg: tang,
            t_s1: trans1,
            t_s2: trans2,
            sit1: Situation::Unknown,
            sit2: Situation::Unknown,
            u_s1: false,
            v_s1: false,
            u_s2: false,
            v_s2: false,
        }
    }

    /// OCCT IntPatch_Line(Tang, Situ1, Situ2) — transitions Touch.
    pub fn new_with_situ(tang: bool, situ1: Situation, situ2: Situation) -> Self {
        IntPatchLine {
            typ: IntPatchIType::Walking,
            tg: tang,
            t_s1: TypeTrans::Touch,
            t_s2: TypeTrans::Touch,
            sit1: situ1,
            sit2: situ2,
            u_s1: false,
            v_s1: false,
            u_s2: false,
            v_s2: false,
        }
    }

    /// OCCT IntPatch_Line(Tang) — transitions Undecided.
    pub fn new(tang: bool) -> Self {
        IntPatchLine {
            typ: IntPatchIType::Walking,
            tg: tang,
            t_s1: TypeTrans::Undecided,
            t_s2: TypeTrans::Undecided,
            sit1: Situation::Unknown,
            sit2: Situation::Unknown,
            u_s1: false,
            v_s1: false,
            u_s2: false,
            v_s2: false,
        }
    }

    /// OCCT SetValue(Uiso1, Viso1, Uiso2, Viso2) — default False.
    pub fn set_value(&mut self, uiso1: bool, viso1: bool, uiso2: bool, viso2: bool) {
        self.u_s1 = uiso1;
        self.v_s1 = viso1;
        self.u_s2 = uiso2;
        self.v_s2 = viso2;
    }

    /// OCCT ArcType().
    pub fn arc_type(&self) -> IntPatchIType {
        self.typ
    }

    /// OCCT IsTangent().
    pub fn is_tangent(&self) -> bool {
        self.tg
    }

    /// OCCT TransitionOnS1().
    pub fn transition_on_s1(&self) -> TypeTrans {
        self.t_s1
    }

    /// OCCT TransitionOnS2().
    pub fn transition_on_s2(&self) -> TypeTrans {
        self.t_s2
    }

    /// OCCT SituationS1() — raises if TransitionOnS1 is not Touch.
    pub fn situation_s1(&self) -> Situation {
        if self.t_s1 != TypeTrans::Touch {
            panic!("IntPatch_Line::SituationS1(): TransitionOnS1 is not Touch");
        }
        self.sit1
    }

    /// OCCT SituationS2() — raises if TransitionOnS2 is not Touch.
    pub fn situation_s2(&self) -> Situation {
        if self.t_s2 != TypeTrans::Touch {
            panic!("IntPatch_Line::SituationS2(): TransitionOnS2 is not Touch");
        }
        self.sit2
    }

    /// OCCT IsUIsoOnS1().
    pub fn is_u_iso_on_s1(&self) -> bool {
        self.u_s1
    }

    /// OCCT IsVIsoOnS1().
    pub fn is_v_iso_on_s1(&self) -> bool {
        self.v_s1
    }

    /// OCCT IsUIsoOnS2().
    pub fn is_u_iso_on_s2(&self) -> bool {
        self.u_s2
    }

    /// OCCT IsVIsoOnS2().
    pub fn is_v_iso_on_s2(&self) -> bool {
        self.v_s2
    }
}
