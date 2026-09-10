//! OCCT TopTrans_SurfaceTransition (TopTrans_SurfaceTransition.cxx / .hxx) —
//! complex transition of a surface relative to a boundary near an
//! interference.
//!
//! 1:1 translation of TopTrans_SurfaceTransition.cxx L25-672 (occ102/occ227
//! touch case included, eap Mar 25 2002).  The rcad CurveTransition (this
//! module) is TopTrans_CurveTransition; this file is the SURFACE variant used
//! by TopClass/TopOpeBRep to combine the transitions of the faces meeting at
//! an intersection line.

use glam::DVec3;
use rcad_kernel::topods::{Orientation, State};

// OCCT L36-38: M_REVERSED / M_INTERNAL / M_UNKNOWN.
fn m_reversed(st: Orientation) -> bool {
    st == Orientation::Reversed
}
fn m_internal(st: Orientation) -> bool {
    st == Orientation::Internal
}

// OCCT L111-119.
const M_UNKNOWN: i32 = -100;
const M_NOUPDATE: i32 = 0;
const M_UPDATE_REF: i32 = 1;
const M_OINTERNAL: i32 = 10;

// OCCT L551-552: BEFORE/AFTER are 1-based array indices (2 / 1); rcad arrays
// are 0-based, so AFTER -> 0, BEFORE -> 1.
const AFTER: usize = 0;
const BEFORE: usize = 1;

/// OCCT gp_Dir::AngleWithRef (gp_Dir.cxx L55-84) — signed angle in [-PI, PI]
/// with the sign from the cross-product against Vref.
fn angle_with_ref(d1: DVec3, d2: DVec3, vref: DVec3) -> f64 {
    let xyz = d1.cross(d2);
    let cosinus = d1.dot(d2);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT FUN_nCinsideS (L25-34): normal to C, tangent to S, oriented INSIDE S.
fn fun_n_cinside_s(tg_c: DVec3, ng_s: DVec3) -> DVec3 {
    ng_s.cross(tg_c)
}

/// OCCT FUN_OO (L40-51): 1 <-> 2 (the OCCT 1-based cos/sin sign values; 0 =
/// null stays 0).
fn fun_oo(i: usize) -> usize {
    if i == 1 {
        2
    } else if i == 2 {
        1
    } else {
        0
    }
}

/// OCCT FUN_Ang (L54-68).  The first parameter (Normref) is unused in OCCT.
fn fun_ang(beafter: DVec3, tg_c: DVec3, norm: DVec3, o: Orientation) -> f64 {
    let mut diron_f = fun_n_cinside_s(tg_c, norm);
    if m_reversed(o) {
        diron_f = -diron_f;
    }
    angle_with_ref(beafter, diron_f, tg_c)
}

/// OCCT FUN_getSTA (L70-92): i = cos sign (0=null), j = sin sign (0=null).
fn fun_get_sta(ang: f64, tola: f64) -> (usize, usize) {
    let cos = ang.cos();
    let sin = ang.sin();
    let nullcos = cos.abs() < tola;
    let nullsin = sin.abs() < tola;
    let i = if nullcos {
        0
    } else if cos > 0.0 {
        1
    } else {
        2
    };
    let j = if nullsin {
        0
    } else if sin > 0.0 {
        1
    } else {
        2
    };
    (i, j)
}

/// OCCT FUN_refnearest (L120-151) — first overload (no curvature).
fn fun_refnearest(
    angref: f64,
    oriref: Orientation,
    ang: f64,
    ori: Orientation,
    tola: f64,
) -> i32 {
    let undef = angref == 100.0;
    if undef {
        return M_UPDATE_REF;
    }
    let cosref = angref.cos();
    let cos = ang.cos();
    let dcos = cosref.abs() - cos.abs();
    if dcos.abs() < tola {
        // Analysis for tangent cases : if two boundary faces are same sided
        // and have tangent normals, if they have opposite orientations
        // we choose INTERNAL as resulting complex transition (case EXTERNAL
        // referring to no logical case).
        if complement(ori) == oriref {
            return M_OINTERNAL;
        } else {
            return M_UNKNOWN; // nyi FUN_RAISE
        }
    }
    if dcos > 0.0 {
        M_NOUPDATE
    } else {
        M_UPDATE_REF
    }
}

/// OCCT FUN_refnearest (L155-271) — second overload with curvature.
#[allow(clippy::too_many_arguments)]
fn fun_refnearest_full(
    i: usize,
    j: usize,
    curv_sref: f64,
    angref: f64,
    oriref: Orientation,
    curvref: f64,
    ang: f64,
    ori: Orientation,
    curv: f64,
    tola: f64,
    touch_flag: &mut bool,
) -> i32 {
    let iisj = i == j;
    let abscos = ang.cos().abs();
    let i0 = (1.0 - abscos).abs() < tola;
    let j0 = abscos < tola;
    let nullcurv = curv == 0.0;
    let curvpos = curv > tola;
    let curvneg = curv < -tola;
    let nullcsref = curv_sref == 0.0;

    let undef = angref == 100.0;
    if undef {
        if i0 {
            if iisj && curvneg {
                return M_NOUPDATE;
            }
            if !iisj && curvpos {
                return M_NOUPDATE;
            }
        }
        if j0 {
            if !nullcsref && (j == 1) && iisj && (curvpos || nullcurv) {
                return M_UPDATE_REF;
            }
            if !nullcsref && (j == 1) && !iisj && (curvneg || nullcurv) {
                return M_UPDATE_REF;
            }
            if iisj && curvpos {
                return M_NOUPDATE;
            }
            if !iisj && curvneg {
                return M_NOUPDATE;
            }
        }
        return M_UPDATE_REF;
    } // undef

    let cosref = angref.cos();
    let cos = ang.cos();
    let dcos = cosref.abs() - cos.abs();
    let samecos = dcos.abs() < tola;
    if samecos {
        // Analysis for tangent cases : if two boundary faces are same sided
        // and have sma dironF.
        if (curvref - curv).abs() < 1.0e-4 {
            if complement(ori) == oriref {
                return M_OINTERNAL;
            } else {
                return M_UNKNOWN; // nyi FUN_RAISE
            }
        }
        let mut noupdate = false;
        if iisj && (curvref > curv) {
            noupdate = true;
        }
        if !iisj && (curvref < curv) {
            noupdate = true;
        }
        let mut updateref = if noupdate { M_NOUPDATE } else { M_UPDATE_REF };
        if !j0 {
            return updateref;
        }
        if !noupdate && !nullcsref {
            // check for (j==1) the face is ABOVE Sref
            // check for (j==2) the face is BELOW Sref
            if (j == 2) && (curv.abs() < curv_sref) {
                updateref = M_NOUPDATE;
            }
            if (j == 1) && (curv.abs() > curv_sref) {
                updateref = M_NOUPDATE;
            }
        }
        return updateref;
    } // samecos

    let updateref = if dcos > 0.0 { M_NOUPDATE } else { M_UPDATE_REF };
    if oriref != ori {
        *touch_flag = true; // eap Mar 25 2002
    }
    updateref
}

/// OCCT TopAbs::Complement — FORWARD<->REVERSED, INTERNAL<->EXTERNAL.
fn complement(or: Orientation) -> Orientation {
    match or {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// OCCT TopTrans_SurfaceTransition (TopTrans_SurfaceTransition.hxx) — complex
/// transition of a surface relative to a boundary near an interference.
pub struct SurfaceTransition {
    my_curv_ref: f64,
    // OCCT NCollection_Array2 (1..2, 1..2) — rcad 0-based [2][2].
    my_ang: [[f64; 2]; 2],
    my_curv: [[f64; 2]; 2],
    my_ori: [[Orientation; 2]; 2],
    my_is_defined: bool,
    my_touch_flag: bool,
    my_norm: DVec3,
    my_tgt: DVec3,
    beafter: DVec3,
}

impl SurfaceTransition {
    /// OCCT TopTrans_SurfaceTransition() (L277-285).
    pub fn new() -> Self {
        SurfaceTransition {
            my_curv_ref: 0.0,
            my_ang: [[100.0; 2]; 2],
            my_curv: [[0.0; 2]; 2],
            my_ori: [[Orientation::Forward; 2]; 2],
            my_is_defined: false,
            my_touch_flag: false,
            my_norm: DVec3::Z,
            my_tgt: DVec3::X,
            beafter: DVec3::Y,
        }
    }

    /// OCCT Reset(Tgt, Norm, MaxD, MinD, MaxCurv, MinCurv) (L287-349).
    pub fn reset_full(
        &mut self,
        tgt: DVec3,
        norm: DVec3,
        max_d: DVec3,
        min_d: DVec3,
        max_curv: f64,
        min_curv: f64,
    ) {
        self.my_is_defined = false;
        self.my_norm = norm;
        self.my_tgt = tgt;
        self.beafter = norm.cross(tgt);

        let tola = rcad_kernel::core::precision::ANGULAR;
        let curismax = max_d.dot(tgt).abs() < tola;
        let curismin = min_d.dot(tgt).abs() < tola;

        if max_curv.abs() < tola && min_curv.abs() < tola {
            self.reset(tgt, norm);
            return;
        }

        if !curismax && !curismin {
            // In the plane normal to <myTgt>, we see the boundary face as
            // a boundary curve.
            // NYIxpu : compute the curvature of the curve if not MaxCurv
            //          nor MinCurv.
            return;
        }

        if curismax {
            self.my_curv_ref = max_curv.abs();
        }
        if curismin {
            self.my_curv_ref = min_curv.abs();
        }
        if self.my_curv_ref < tola {
            self.my_curv_ref = 0.0;
        }

        // ============================================================
        // recall : <Norm> is oriented OUTSIDE the "geometric matter" described
        //          by the surface
        //          -  if (myCurvRef != 0.) Sref is UNDER axis (sin = 0)
        //             referential (beafter,myNorm,myTgt)  -
        // ============================================================

        for i in 0..2 {
            for j in 0..2 {
                self.my_ang[i][j] = 100.0;
            }
        }

        self.my_touch_flag = false; // eap Mar 25 2002
        self.my_is_defined = true;
    }

    /// OCCT Reset(Tgt, Norm) (L351-370).
    pub fn reset(&mut self, tgt: DVec3, norm: DVec3) {
        self.my_is_defined = false;
        // beafter oriented (before, after) the intersection on the reference
        // surface.
        self.my_norm = norm;
        self.my_tgt = tgt;
        self.beafter = norm.cross(tgt);
        for i in 0..2 {
            for j in 0..2 {
                self.my_ang[i][j] = 100.0;
            }
        }
        self.my_curv_ref = 0.0;
        self.my_touch_flag = false; // eap Mar 25 2002
        self.my_is_defined = true;
    }

    /// OCCT Compare(Tole, Norm, MaxD, MinD, MaxCurv, MinCurv, S, O)
    /// (L372-490).
    #[allow(clippy::too_many_arguments)]
    pub fn compare_full(
        &mut self,
        tole: f64,
        norm: DVec3,
        max_d: DVec3,
        min_d: DVec3,
        max_curv: f64,
        min_curv: f64,
        s: Orientation,
        o: Orientation,
    ) {
        if !self.my_is_defined {
            return;
        }
        let mut curv = 0.0;
        // ------
        let tola = if tole > 0.0 { tole } else { rcad_kernel::core::precision::ANGULAR };
        let curismax = max_d.dot(self.my_tgt).abs() < tola;
        let curismin = min_d.dot(self.my_tgt).abs() < tola;
        if !curismax && !curismin {
            // In the plane normal to <myTgt>, we see the boundary face as
            // a boundary curve.
            // NYIxpu : compute the curvature of the curve if not MaxCurv
            //          nor MinCurv.
            self.my_is_defined = false;
            return;
        }
        if curismax {
            curv = max_curv.abs();
        }
        if curismin {
            curv = min_curv.abs();
        }
        if self.my_curv_ref < tola {
            curv = 0.0;
        }
        let diron_f = fun_n_cinside_s(self.my_tgt, norm);
        let prod = diron_f.cross(norm).dot(self.my_tgt);
        if prod < 0.0 {
            curv = -curv;
        }

        let ang = fun_ang(self.beafter, self.my_tgt, norm, o);

        // i = 0,1,2 : cos = 0,>0,<0
        // j = 0,1,2 : sin = 0,>0,<0
        let (mut i, mut j) = fun_get_sta(ang, tola);

        // update nearest :
        // ---------------
        let kmax = if m_internal(o) { 2 } else { 1 };
        for k in 0..kmax {
            if k == 1 {
                // get the opposite Ang
                i = fun_oo(i);
                j = fun_oo(j);
            }
            let i0 = i == 0;
            let j0 = j == 0;
            let nmax = if i0 || j0 { 2 } else { 1 };
            for _n in 0..nmax {
                let n = _n + 1;
                if i0 {
                    i = n;
                }
                if j0 {
                    j = n;
                }
                // i/j are the OCCT 1-based cos/sin signs {1,2} here; the
                // arrays are 0-based.
                let (ia, ja) = (i - 1, j - 1);

                let refn = fun_refnearest_full(
                    i,
                    j,
                    self.my_curv_ref,
                    self.my_ang[ia][ja],
                    self.my_ori[ia][ja],
                    self.my_curv[ia][ja],
                    ang,
                    s,
                    curv,
                    tola,
                    &mut self.my_touch_flag,
                ); // eap Mar 25 2002
                if refn == M_UNKNOWN {
                    self.my_is_defined = false;
                    return;
                }
                if refn > 0 {
                    self.my_ang[ia][ja] = ang;
                    self.my_ori[ia][ja] = if refn == M_OINTERNAL {
                        Orientation::Internal
                    } else {
                        s
                    };
                    self.my_curv[ia][ja] = curv;
                }
            } // n=1..nmax
        } // k=1..kmax
    }

    /// OCCT Compare(Tole, Norm, S, O) (L492-549).
    pub fn compare(&mut self, tole: f64, norm: DVec3, s: Orientation, o: Orientation) {
        if !self.my_is_defined {
            return;
        }
        // oriented Ang(beafter,dironF),
        // dironF normal to the curve, oriented INSIDE F, the added oriented
        // support.
        let ang = fun_ang(self.beafter, self.my_tgt, norm, o);
        let tola = if tole > 0.0 { tole } else { rcad_kernel::core::precision::ANGULAR };

        // i = 0,1,2 : cos = 0,>0,<0
        // j = 0,1,2 : sin = 0,>0,<0
        let (mut i, mut j) = fun_get_sta(ang, tola);

        let kmax = if m_internal(o) { 2 } else { 1 };
        for k in 0..kmax {
            if k == 1 {
                // get the opposite Ang
                i = fun_oo(i);
                j = fun_oo(j);
            }
            let i0 = i == 0;
            let j0 = j == 0;
            let nmax = if i0 || j0 { 2 } else { 1 };
            for _n in 0..nmax {
                let n = _n + 1;
                if i0 {
                    i = n;
                }
                if j0 {
                    j = n;
                }
                // i/j are the OCCT 1-based cos/sin signs {1,2} here; the
                // arrays are 0-based.
                let (ia, ja) = (i - 1, j - 1);

                let refn = fun_refnearest(self.my_ang[ia][ja], self.my_ori[ia][ja], ang, s, tola);
                if refn == M_UNKNOWN {
                    self.my_is_defined = false;
                    return;
                }
                if refn > 0 {
                    self.my_ang[ia][ja] = ang;
                    self.my_ori[ia][ja] = if refn == M_OINTERNAL {
                        Orientation::Internal
                    } else {
                        s
                    };
                }
            } // n=1..nmax
        } // k=1..kmax
    }

    /// OCCT FUN_getstate (L554-586): state from the angle/orientation arrays
    /// at row iSTA with the before/after index iINDEX.
    fn get_state(&self, i_sta: usize, i_index: usize) -> State {
        let a1 = self.my_ang[i_sta][0];
        let a2 = self.my_ang[i_sta][1];
        let undef1 = a1 == 100.0;
        let undef2 = a2 == 100.0;
        let undef = undef1 && undef2;
        if undef {
            return State::Unknown;
        }
        if undef1 || undef2 {
            let jok = if undef1 { 1 } else { 0 };
            let o = self.my_ori[i_sta][jok];
            return if i_index == BEFORE {
                Self::get_before(o)
            } else {
                Self::get_after(o)
            };
        }
        let o1 = self.my_ori[i_sta][0];
        let o2 = self.my_ori[i_sta][1];
        let st1 = if i_index == BEFORE { Self::get_before(o1) } else { Self::get_after(o1) };
        let st2 = if i_index == BEFORE { Self::get_before(o2) } else { Self::get_after(o2) };
        if st1 != st2 {
            return State::Unknown; // Incoherent data
        }
        st1
    }

    /// OCCT StateBefore() (L588-616).
    pub fn state_before(&self) -> State {
        if !self.my_is_defined {
            return State::Unknown;
        }
        // we take the state before of before orientations
        let mut before = self.get_state(BEFORE, BEFORE);
        if before == State::Unknown {
            // looking back in before for defined states
            // we take the state before of after orientations
            before = self.get_state(AFTER, BEFORE);
            // eap Mar 25 2002
            if self.my_touch_flag {
                if before == State::Out {
                    before = State::In;
                } else if before == State::In {
                    before = State::Out;
                }
            }
        }
        before
    }

    /// OCCT StateAfter() (L618-644).
    pub fn state_after(&self) -> State {
        if !self.my_is_defined {
            return State::Unknown;
        }
        let mut after = self.get_state(AFTER, AFTER);
        if after == State::Unknown {
            // looking back in before for defined states
            after = self.get_state(BEFORE, AFTER);
            // eap Mar 25 2002
            if self.my_touch_flag {
                if after == State::Out {
                    after = State::In;
                } else if after == State::In {
                    after = State::Out;
                }
            }
        }
        after
    }

    /// OCCT GetBefore(Tran) (L646-658).
    pub fn get_before(tran: Orientation) -> State {
        match tran {
            Orientation::Forward | Orientation::External => State::Out,
            Orientation::Reversed | Orientation::Internal => State::In,
            _ => State::Out,
        }
    }

    /// OCCT GetAfter(Tran) (L660-672).
    pub fn get_after(tran: Orientation) -> State {
        match tran {
            Orientation::Forward | Orientation::Internal => State::In,
            Orientation::Reversed | Orientation::External => State::Out,
            _ => State::Out,
        }
    }
}

impl Default for SurfaceTransition {
    fn default() -> Self {
        Self::new()
    }
}
