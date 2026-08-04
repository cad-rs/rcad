// OCCT IntSurf_Transition / IntSurf / IntSurf::MakeTransition / Recadre
// IntSurf_Transition.hxx/.lxx/.cxx, IntSurf.hxx/.cxx,
// IntPatch_ImpImpIntersection.cxx Recadre (L316-435).
//
// 1:1 Rust translation. rcad data-model notes:
//   - gp_Vec / gp_Dir -> DVec3.
//   - The rcad surfaces are Surface3; Recadre uses periodic flags + first/last
//     parameters which rcad exposes via SurfaceEval (is_u_periodic, etc.).

use glam::DVec3;
use rcad_kernel::geom::{Surface3, SurfaceEval};

/// OCCT IntSurf_TypeTrans.hxx.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeTrans {
    In,
    Out,
    Touch,
    Undecided,
}

/// OCCT IntSurf_Situation.hxx.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Situation {
    Inside,
    Outside,
    Unknown,
}

/// OCCT IntSurf_Transition.hxx — definition of the transition at the
/// intersection between an intersection line and a restriction curve.
#[derive(Debug, Clone, Copy)]
pub struct Transition {
    tangent: bool,
    typetra: TypeTrans,
    situat: Situation,
    oppos: bool,
}

impl Transition {
    /// OCCT empty constructor — creates an UNDECIDED transition.
    pub fn new() -> Self {
        Transition {
            tangent: false,
            typetra: TypeTrans::Undecided,
            situat: Situation::Unknown,
            oppos: false,
        }
    }

    /// OCCT IntSurf_Transition(Tangent, Type).
    pub fn new_in_out(tangent: bool, typ: TypeTrans) -> Self {
        Transition {
            tangent,
            typetra: typ,
            situat: Situation::Unknown,
            oppos: false,
        }
    }

    /// OCCT IntSurf_Transition(Tangent, Situ, Oppos).
    pub fn new_touch(tangent: bool, situ: Situation, oppos: bool) -> Self {
        Transition {
            tangent,
            typetra: TypeTrans::Touch,
            situat: situ,
            oppos,
        }
    }

    /// OCCT SetValue(Tangent, Type).
    pub fn set_value_in_out(&mut self, tangent: bool, typ: TypeTrans) {
        self.tangent = tangent;
        self.typetra = typ;
    }

    /// rcad helper — build a transition from an IntSurf_TypeTrans alone.
    pub fn from_type(typ: TypeTrans) -> Self {
        Transition {
            tangent: false,
            typetra: typ,
            situat: Situation::Unknown,
            oppos: false,
        }
    }

    /// OCCT SetValue(Tangent, Situ, Oppos).
    pub fn set_value_touch(&mut self, tangent: bool, situ: Situation, oppos: bool) {
        self.tangent = tangent;
        self.typetra = TypeTrans::Touch;
        self.situat = situ;
        self.oppos = oppos;
    }

    /// OCCT SetValue().
    pub fn set_value_undecided(&mut self) {
        self.typetra = TypeTrans::Undecided;
    }

    /// OCCT TransitionType().
    pub fn transition_type(&self) -> TypeTrans {
        self.typetra
    }

    /// OCCT IsTangent().
    pub fn is_tangent(&self) -> bool {
        if self.typetra == TypeTrans::Undecided {
            panic!("Transition::IsTangent on UNDECIDED");
        }
        self.tangent
    }

    /// OCCT IsTangent() — non-throwing variant used when the transition type
    /// may be Undecided (the ALine-level transition in MakeWLine).
    pub fn tangent(&self) -> bool {
        self.tangent
    }

    /// OCCT Situation().
    pub fn situation(&self) -> Situation {
        if self.typetra != TypeTrans::Touch {
            panic!("Transition::Situation on non-TOUCH");
        }
        self.situat
    }

    /// OCCT IsOpposite().
    pub fn is_opposite(&self) -> bool {
        if self.typetra != TypeTrans::Touch {
            panic!("Transition::IsOpposite on non-TOUCH");
        }
        self.oppos
    }
}

/// OCCT IntSurf::MakeTransition (IntSurf.cxx L28-75).
///
/// Computes the transition of the intersection point between two lines.
/// TgFirst is the tangent of the intersection line, TgSecond the tangent of
/// the restriction, Normale the normal used to orient the cross product.
pub fn make_transition(
    tg_first: DVec3,
    tg_second: DVec3,
    normale: DVec3,
    t_first: &mut Transition,
    t_second: &mut Transition,
) {
    // Compute the mixed product of normal, tangent 1, tangent 2.
    let pvect = tg_second.cross(tg_first);

    let n_tg_second = tg_second.length();
    let n_tg_first = tg_first.length();
    let n_tg_second_n_tg_first_angular = n_tg_second * n_tg_first * rcad_kernel::precision::ANGULAR;

    if n_tg_first <= rcad_kernel::precision::CONFUSION {
        t_first.set_value_in_out(true, TypeTrans::Undecided);
        t_second.set_value_in_out(true, TypeTrans::Undecided);
    } else if (n_tg_second <= rcad_kernel::precision::CONFUSION)
        || (pvect.length() <= n_tg_second_n_tg_first_angular)
    {
        t_first.set_value_touch(true, Situation::Unknown, tg_first.dot(tg_second) < 0.0);
        t_second.set_value_touch(true, Situation::Unknown, tg_first.dot(tg_second) < 0.0);
    } else {
        let mut yu = pvect.dot(normale);
        yu /= n_tg_second * n_tg_first;
        if yu > 0.0001 {
            t_first.set_value_in_out(false, TypeTrans::In);
            t_second.set_value_in_out(false, TypeTrans::Out);
        } else if yu < -0.0001 {
            t_first.set_value_in_out(false, TypeTrans::Out);
            t_second.set_value_in_out(false, TypeTrans::In);
        } else {
            t_first.set_value_in_out(true, TypeTrans::Undecided);
            t_second.set_value_in_out(true, TypeTrans::Undecided);
        }
    }
}

/// OCCT Recadre (IntPatch_ImpImpIntersection.cxx L316-435).
///
/// Shifts the parameters (u1,v1) and (u2,v2) into the natural domains of the
/// periodic parametric directions of the two surfaces.
pub fn recadre(
    s1: &Surface3,
    s2: &Surface3,
    u1: &mut f64,
    v1: &mut f64,
    u2: &mut f64,
    v2: &mut f64,
) {
    let lmf = std::f64::consts::TAU;

    // Surface 1.
    let s1_u_periodic = s1.is_u_periodic();
    let s1_v_periodic = s1.is_v_periodic();
    let d1 = s1.default_domain();
    if s1_u_periodic {
        let f = d1[0];
        let l = d1[1];
        let fpls2 = 0.5 * (f + l);
        while (*u1 < f) && ((fpls2 - *u1) > (*u1 + lmf - fpls2)) {
            *u1 += lmf;
        }
        while (*u1 > l) && ((*u1 - fpls2) > (fpls2 - (*u1 - lmf))) {
            *u1 -= lmf;
        }
    }
    if s1_v_periodic {
        let f = d1[2];
        let l = d1[3];
        let fpls2 = 0.5 * (f + l);
        while (*v1 < f) && ((fpls2 - *v1) > (*v1 + lmf - fpls2)) {
            *v1 += lmf;
        }
        while (*v1 > l) && ((*v1 - fpls2) > (fpls2 - (*v1 - lmf))) {
            *v1 -= lmf;
        }
    }

    // Surface 2.
    let s2_u_periodic = s2.is_u_periodic();
    let s2_v_periodic = s2.is_v_periodic();
    let d2 = s2.default_domain();
    if s2_u_periodic {
        let f = d2[0];
        let l = d2[1];
        let fpls2 = 0.5 * (f + l);
        while (*u2 < f) && ((fpls2 - *u2) > (*u2 + lmf - fpls2)) {
            *u2 += lmf;
        }
        while (*u2 > l) && ((*u2 - fpls2) > (fpls2 - (*u2 - lmf))) {
            *u2 -= lmf;
        }
    }
    if s2_v_periodic {
        let f = d2[2];
        let l = d2[3];
        let fpls2 = 0.5 * (f + l);
        while (*v2 < f) && ((fpls2 - *v2) > (*v2 + lmf - fpls2)) {
            *v2 += lmf;
        }
        while (*v2 > l) && ((*v2 - fpls2) > (fpls2 - (*v2 - lmf))) {
            *v2 -= lmf;
        }
    }
}
