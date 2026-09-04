// OCCT AppDef (TKGeomBase) — 1:1 Rust translation of
//   AppDef_MultiPointConstraint.hxx/.cxx (L26-317)
//   AppDef_MultiLine.hxx/.cxx (L23-103).
//
// AppDef_MultiPointConstraint extends AppParCurves_MultiPoint (rcad
// `MultiPoint` in approx_int.rs) with tangency and curvature constraints.
// Rust has no inheritance, so the base is held as a `base` member
// (composition), matching the trait-object mapping used elsewhere.
//
// Note: AppDef_MultiLine.hxx L95 declares SetParameter(Index, U) but the
// OCCT library contains no definition of it (dead declaration, link error if
// called) — it is intentionally not ported.

use glam::{DVec2, DVec3};

use super::approx_int::{ApproxStatus, MultiPoint};

/// OCCT AppDef_MyLineTool (AppDef_MyLineTool.hxx/.cxx) — "Example of
/// MultiLine tool corresponding to the tools of the packages AppParCurves
/// and Approx. For Approx, the tool will not add points if the algorithms
/// want some." This is the ToolLine binding of the AppDef instantiations of
/// AppParCurves_LeastSquare.gxx / AppParCurves_ComputeLine.gxx; the OCCT
/// statics are mapped to free functions of this module. Output arrays are
/// passed as slices sized for the multipoint (the OCCT NCollection_Array1
/// Lower() offset is carried by the slice index).
pub mod my_line_tool {
    use super::{ApproxStatus, DVec2, DVec3, MultiLine};

    /// OCCT AppDef_MyLineTool::FirstPoint (AppDef_MyLineTool.cxx L16-19).
    pub fn first_point(_ml: &MultiLine) -> usize {
        1
    }

    /// OCCT AppDef_MyLineTool::LastPoint (AppDef_MyLineTool.cxx L21-24).
    pub fn last_point(ml: &MultiLine) -> usize {
        ml.nb_multi_points()
    }

    /// OCCT AppDef_MyLineTool::NbP2d (AppDef_MyLineTool.cxx L26-29).
    pub fn nb_p2d(ml: &MultiLine) -> usize {
        ml.value(1).base.nb_points2d()
    }

    /// OCCT AppDef_MyLineTool::NbP3d (AppDef_MyLineTool.cxx L31-34).
    pub fn nb_p3d(ml: &MultiLine) -> usize {
        ml.value(1).base.nb_points()
    }

    /// OCCT AppDef_MyLineTool::Value(ML, MPointIndex, tabPt) — the 3d points
    /// of the multipoint MPointIndex when only 3d points exist
    /// (AppDef_MyLineTool.cxx L36-45).
    pub fn value_3d(ml: &MultiLine, mpoint_index: usize, tab_pt: &mut [DVec3]) {
        let mpc = ml.value(mpoint_index);
        let nbp3d = mpc.base.nb_points();
        for i in 1..=nbp3d {
            tab_pt[i - 1] = mpc.base.point(i);
        }
    }

    /// OCCT AppDef_MyLineTool::Value(ML, MPointIndex, tabPt2d) — the 2d
    /// points of the multipoint MPointIndex when only 2d points exist
    /// (AppDef_MyLineTool.cxx L47-58).
    pub fn value_2d(ml: &MultiLine, mpoint_index: usize, tab_pt2d: &mut [DVec2]) {
        let mpc = ml.value(mpoint_index);
        let nbp3d = mpc.base.nb_points();
        let nbp2d = mpc.base.nb_points2d();
        for i in 1..=nbp2d {
            tab_pt2d[i - 1] = mpc.base.point2d(nbp3d + i);
        }
    }

    /// OCCT AppDef_MyLineTool::Value(ML, MPointIndex, tabPt, tabPt2d) — the
    /// 3d and 2d points of the multipoint MPointIndex
    /// (AppDef_MyLineTool.cxx L60-76).
    pub fn value_3d_2d(
        ml: &MultiLine,
        mpoint_index: usize,
        tab_pt: &mut [DVec3],
        tab_pt2d: &mut [DVec2],
    ) {
        let mpc = ml.value(mpoint_index);
        let nbp3d = mpc.base.nb_points();
        let nbp2d = mpc.base.nb_points2d();
        for i in 1..=nbp3d {
            tab_pt[i - 1] = mpc.base.point(i);
        }
        for i in 1..=nbp2d {
            tab_pt2d[i - 1] = mpc.base.point2d(nbp3d + i);
        }
    }

    /// OCCT AppDef_MyLineTool::Tangency(ML, MPointIndex, tabV) — the 3d
    /// tangencies of the multipoint MPointIndex when only 3d points exist
    /// (AppDef_MyLineTool.cxx L78-93).
    pub fn tangency_3d(ml: &MultiLine, mpoint_index: usize, tab_v: &mut [DVec3]) -> bool {
        let mpc = ml.value(mpoint_index);
        if mpc.is_tangency_point() {
            let nbp3d = mpc.base.nb_points();
            for i in 1..=nbp3d {
                tab_v[i - 1] = mpc.tang(i);
            }
            true
        } else {
            false
        }
    }

    /// OCCT AppDef_MyLineTool::Tangency(ML, MPointIndex, tabV2d) — the 2d
    /// tangencies of the multipoint MPointIndex when only 2d points exist
    /// (AppDef_MyLineTool.cxx L95-110).
    pub fn tangency_2d(ml: &MultiLine, mpoint_index: usize, tab_v2d: &mut [DVec2]) -> bool {
        let mpc = ml.value(mpoint_index);
        if mpc.is_tangency_point() {
            let nbp3d = mpc.base.nb_points();
            let nbp2d = mpc.base.nb_points2d();
            for i in 1..=nbp2d {
                tab_v2d[i - 1] = mpc.tang2d(nbp3d + i);
            }
            true
        } else {
            false
        }
    }

    /// OCCT AppDef_MyLineTool::Tangency(ML, MPointIndex, tabV, tabV2d) — the
    /// 3d and 2d tangencies of the multipoint MPointIndex
    /// (AppDef_MyLineTool.cxx L112-131).
    pub fn tangency_3d_2d(
        ml: &MultiLine,
        mpoint_index: usize,
        tab_v: &mut [DVec3],
        tab_v2d: &mut [DVec2],
    ) -> bool {
        let mpc = ml.value(mpoint_index);
        if mpc.is_tangency_point() {
            let nbp3d = mpc.base.nb_points();
            let nbp2d = mpc.base.nb_points2d();
            for i in 1..=nbp3d {
                tab_v[i - 1] = mpc.tang(i);
            }
            for i in 1..=nbp2d {
                tab_v2d[i - 1] = mpc.tang2d(nbp3d + i);
            }
            true
        } else {
            false
        }
    }

    /// OCCT AppDef_MyLineTool::MakeMLBetween (AppDef_MyLineTool.cxx L152-158)
    /// — "Is never called in the algorithms. Nothing is done." (stub
    /// returning theML).
    pub fn make_ml_between(ml: &MultiLine, _i1: usize, _i2: usize, _nb_p_min: usize) -> MultiLine {
        ml.clone() // stub
    }

    /// OCCT AppDef_MyLineTool::MakeMLOneMorePoint
    /// (AppDef_MyLineTool.cxx L160-164) — "Is never called in the
    /// algorithms. Nothing is done."
    pub fn make_ml_one_more_point(
        _ml: &MultiLine,
        _i1: usize,
        _i2: usize,
        _indbad: usize,
        _other_line: &mut MultiLine,
    ) -> bool {
        false
    }

    /// OCCT AppDef_MyLineTool::WhatStatus (AppDef_MyLineTool.cxx L166-170) —
    /// returns Approx_NoPointsAdded.
    pub fn what_status(_ml: &MultiLine, _i1: usize, _i2: usize) -> ApproxStatus {
        ApproxStatus::NoPointsAdded
    }

    /// OCCT AppDef_MyLineTool::Curvature(ML, MPointIndex, tabV) — the 3d
    /// curvatures of the multipoint MPointIndex when only 3d points exist
    /// (AppDef_MyLineTool.cxx L172-187).
    pub fn curvature_3d(ml: &MultiLine, mpoint_index: usize, tab_v: &mut [DVec3]) -> bool {
        let mpc = ml.value(mpoint_index);
        if mpc.is_curvature_point() {
            let nbp3d = mpc.base.nb_points();
            for i in 1..=nbp3d {
                tab_v[i - 1] = mpc.curv(i);
            }
            true
        } else {
            false
        }
    }

    /// OCCT AppDef_MyLineTool::Curvature(ML, MPointIndex, tabV2d) — the 2d
    /// curvatures of the multipoint MPointIndex when only 2d points exist
    /// (AppDef_MyLineTool.cxx L189-204).
    pub fn curvature_2d(ml: &MultiLine, mpoint_index: usize, tab_v2d: &mut [DVec2]) -> bool {
        let mpc = ml.value(mpoint_index);
        if mpc.is_curvature_point() {
            let nbp3d = mpc.base.nb_points();
            let nbp2d = mpc.base.nb_points2d();
            for i in 1..=nbp2d {
                tab_v2d[i - 1] = mpc.curv2d(nbp3d + i);
            }
            true
        } else {
            false
        }
    }

    /// OCCT AppDef_MyLineTool::Curvature(ML, MPointIndex, tabV, tabV2d) — the
    /// 3d and 2d curvatures of the multipoint MPointIndex
    /// (AppDef_MyLineTool.cxx L206-225).
    pub fn curvature_3d_2d(
        ml: &MultiLine,
        mpoint_index: usize,
        tab_v: &mut [DVec3],
        tab_v2d: &mut [DVec2],
    ) -> bool {
        let mpc = ml.value(mpoint_index);
        if mpc.is_curvature_point() {
            let nbp3d = mpc.base.nb_points();
            let nbp2d = mpc.base.nb_points2d();
            for i in 1..=nbp3d {
                tab_v[i - 1] = mpc.curv(i);
            }
            for i in 1..=nbp2d {
                tab_v2d[i - 1] = mpc.curv2d(nbp3d + i);
            }
            true
        } else {
            false
        }
    }
}

// ---------------------------------------------------------------------------
// AppDef_MultiPointConstraint
// ---------------------------------------------------------------------------

/// OCCT AppDef_MultiPointConstraint — describes a MultiPointConstraint used
/// in a Multiline. MultiPointConstraints are composed of several two or
/// three-dimensional points. The purpose is to define the corresponding
/// points that share a common constraint in order to compute the
/// approximation of several lines in parallel.
///
/// Notes (AppDef_MultiPointConstraint.hxx L41-48):
/// - The order of points of a MultiPointConstraint is very important. Users
///   must give 3D points first, and then 2D points.
/// - The constraints for the points included in a MultiPointConstraint are
///   always identical for all points, including the parameter.
/// - If a MultiPointConstraint is a "tangency" point, the point is also a
///   "passing" point.
#[derive(Debug, Clone, Default)]
pub struct MultiPointConstraint {
    /// OCCT AppParCurves_MultiPoint base part.
    pub base: MultiPoint,
    /// OCCT tabTang (null handle == empty vec).
    tab_tang: Vec<DVec3>,
    /// OCCT tabCurv.
    tab_curv: Vec<DVec3>,
    /// OCCT tabTang2d.
    tab_tang2d: Vec<DVec2>,
    /// OCCT tabCurv2d.
    tab_curv2d: Vec<DVec2>,
}

impl MultiPointConstraint {
    /// OCCT AppDef_MultiPointConstraint() (cxx L26) — creates an undefined
    /// MultiPointConstraint.
    pub fn new() -> Self {
        MultiPointConstraint {
            base: MultiPoint::default(),
            tab_tang: Vec::new(),
            tab_curv: Vec::new(),
            tab_tang2d: Vec::new(),
            tab_curv2d: Vec::new(),
        }
    }

    /// OCCT AppDef_MultiPointConstraint(NbPoles, NbPoles2d) (cxx L28-31).
    pub fn new_nb(nb_poles: usize, nb_poles2d: usize) -> Self {
        MultiPointConstraint {
            base: MultiPoint::new(nb_poles, nb_poles2d),
            tab_tang: Vec::new(),
            tab_curv: Vec::new(),
            tab_tang2d: Vec::new(),
            tab_curv2d: Vec::new(),
        }
    }

    /// OCCT AppDef_MultiPointConstraint(tabP) (cxx L33-36) — a MultiPoint
    /// only composed of 3D points.
    pub fn new_tab_p(tab_p: &[DVec3]) -> Self {
        MultiPointConstraint {
            base: MultiPoint::new_tab_p3d(tab_p),
            ..MultiPointConstraint::new()
        }
    }

    /// OCCT AppDef_MultiPointConstraint(tabP2d) (cxx L38-41) — a MultiPoint
    /// only composed of 2D points.
    pub fn new_tab_p2d(tab_p2d: &[DVec2]) -> Self {
        MultiPointConstraint {
            base: MultiPoint::new_tab_p2d(tab_p2d),
            ..MultiPointConstraint::new()
        }
    }

    /// OCCT AppDef_MultiPointConstraint(tabP, tabP2d) (cxx L43-47).
    pub fn new_tab_p_p2d(tab_p: &[DVec3], tab_p2d: &[DVec2]) -> Self {
        MultiPointConstraint {
            base: MultiPoint::new_tab_p3d_p2d(tab_p, tab_p2d),
            ..MultiPointConstraint::new()
        }
    }

    /// OCCT AppDef_MultiPointConstraint(tabP, tabP2d, tabVec, tabVec2d,
    /// tabCur, tabCur2d) (cxx L49-92) — creates a MultiPointConstraint with
    /// a constraint of Curvature.
    #[allow(clippy::too_many_arguments)]
    pub fn new_curvature(
        tab_p: &[DVec3],
        tab_p2d: &[DVec2],
        tab_vec: &[DVec3],
        tab_vec2d: &[DVec2],
        tab_cur: &[DVec3],
        tab_cur2d: &[DVec2],
    ) -> Self {
        let mut r = MultiPointConstraint {
            base: MultiPoint::new_tab_p3d_p2d(tab_p, tab_p2d),
            ..MultiPointConstraint::new()
        };

        if (tab_p.len() != tab_vec.len())
            || (tab_p2d.len() != tab_vec2d.len())
            || (tab_cur.len() != tab_p.len())
            || (tab_cur2d.len() != tab_p2d.len())
        {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }

        r.tab_tang = tab_vec.to_vec();
        r.tab_tang2d = tab_vec2d.to_vec();
        // OCCT loops the tabCur fill with `i <= tabVec.Length()` (cxx L83) —
        // the tabCur length equals tabP length == tabVec length by the check
        // above, so the copies agree; rcad copies tab_cur directly.
        r.tab_curv = tab_cur.to_vec();
        r.tab_curv2d = tab_cur2d.to_vec();
        r
    }

    /// OCCT AppDef_MultiPointConstraint(tabP, tabP2d, tabVec, tabVec2d)
    /// (cxx L94-120) — creates a MultiPointConstraint with a constraint of
    /// Tangency.
    pub fn new_tangency(
        tab_p: &[DVec3],
        tab_p2d: &[DVec2],
        tab_vec: &[DVec3],
        tab_vec2d: &[DVec2],
    ) -> Self {
        let mut r = MultiPointConstraint {
            base: MultiPoint::new_tab_p3d_p2d(tab_p, tab_p2d),
            ..MultiPointConstraint::new()
        };

        if (tab_p.len() != tab_vec.len()) || (tab_p2d.len() != tab_vec2d.len()) {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }

        r.tab_tang = tab_vec.to_vec();
        r.tab_tang2d = tab_vec2d.to_vec();
        r
    }

    /// OCCT AppDef_MultiPointConstraint(tabP, tabVec) (cxx L122-139) — 3d
    /// points with constraints of tangency.
    pub fn new_p3d_tangency(tab_p: &[DVec3], tab_vec: &[DVec3]) -> Self {
        let mut r = MultiPointConstraint {
            base: MultiPoint::new_tab_p3d(tab_p),
            ..MultiPointConstraint::new()
        };

        if tab_p.len() != tab_vec.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }

        r.tab_tang = tab_vec.to_vec();
        r
    }

    /// OCCT AppDef_MultiPointConstraint(tabP, tabVec, tabCur) (cxx L141-165)
    /// — 3d points with constraints of curvature.
    pub fn new_p3d_curvature(tab_p: &[DVec3], tab_vec: &[DVec3], tab_cur: &[DVec3]) -> Self {
        let mut r = MultiPointConstraint {
            base: MultiPoint::new_tab_p3d(tab_p),
            ..MultiPointConstraint::new()
        };

        if (tab_p.len() != tab_vec.len()) || (tab_p.len() != tab_cur.len()) {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }

        r.tab_tang = tab_vec.to_vec();
        r.tab_curv = tab_cur.to_vec();
        r
    }

    /// OCCT AppDef_MultiPointConstraint(tabP2d, tabVec2d) (cxx L167-186) —
    /// 2d points with constraints of tangency.
    pub fn new_p2d_tangency(tab_p2d: &[DVec2], tab_vec2d: &[DVec2]) -> Self {
        let mut r = MultiPointConstraint {
            base: MultiPoint::new_tab_p2d(tab_p2d),
            ..MultiPointConstraint::new()
        };

        if tab_p2d.len() != tab_vec2d.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }

        r.tab_tang2d = tab_vec2d.to_vec();
        r
    }

    /// OCCT AppDef_MultiPointConstraint(tabP2d, tabVec2d, tabCur2d)
    /// (cxx L188-213) — 2d points with constraints of curvature.
    pub fn new_p2d_curvature(tab_p2d: &[DVec2], tab_vec2d: &[DVec2], tab_cur2d: &[DVec2]) -> Self {
        let mut r = MultiPointConstraint {
            base: MultiPoint::new_tab_p2d(tab_p2d),
            ..MultiPointConstraint::new()
        };

        if (tab_p2d.len() != tab_vec2d.len()) || (tab_cur2d.len() != tab_p2d.len()) {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }

        r.tab_tang2d = tab_vec2d.to_vec();
        r.tab_curv2d = tab_cur2d.to_vec();
        r
    }

    /// OCCT SetTang(Index, Tang) (cxx L215-226) — sets the value of the
    /// tangency of the point of range Index (1-based 3D index).
    pub fn set_tang(&mut self, index: usize, tang: DVec3) {
        if self.tab_tang.is_empty() {
            self.tab_tang = vec![DVec3::ZERO; self.base.nb_p3d];
        }
        assert!(
            index >= 1 && index <= self.base.nb_p3d,
            "Standard_OutOfRange: AppDef_MultiPointConstraint::SetTang"
        );
        self.tab_tang[index - 1] = tang;
    }

    /// OCCT Tang(Index) (cxx L228-235) — the tangency value of the point of
    /// range Index.
    pub fn tang(&self, index: usize) -> DVec3 {
        assert!(
            index >= 1 && index <= self.base.nb_p3d,
            "Standard_OutOfRange: AppDef_MultiPointConstraint::Tang"
        );
        self.tab_tang[index - 1]
    }

    /// OCCT SetTang2d(Index, Tang2d) (cxx L237-249) — Index is the GLOBAL
    /// point index (nbP+1 .. nbP+nbP2d).
    pub fn set_tang2d(&mut self, index: usize, tang2d: DVec2) {
        if self.tab_tang2d.is_empty() {
            self.tab_tang2d = vec![DVec2::ZERO; self.base.p2d.len()];
        }
        let nbp = self.base.nb_p3d;
        assert!(
            index > nbp && index <= nbp + self.base.p2d.len(),
            "Standard_OutOfRange: AppDef_MultiPointConstraint::SetTang2d"
        );
        self.tab_tang2d[index - nbp - 1] = tang2d;
    }

    /// OCCT Tang2d(Index) (cxx L251-258) — Index is the GLOBAL point index.
    pub fn tang2d(&self, index: usize) -> DVec2 {
        let nbp = self.base.nb_p3d;
        assert!(
            index > nbp && index <= nbp + self.base.p2d.len(),
            "Standard_OutOfRange: AppDef_MultiPointConstraint::Tang2d"
        );
        self.tab_tang2d[index - nbp - 1]
    }

    /// OCCT SetCurv(Index, Curv) (cxx L260-271) — sets the value of the
    /// normal vector at the point of index Index; the norm of the normal
    /// vector at the point of position Index is set to the normal curvature.
    pub fn set_curv(&mut self, index: usize, curv: DVec3) {
        if self.tab_curv.is_empty() {
            self.tab_curv = vec![DVec3::ZERO; self.base.nb_p3d];
        }
        assert!(
            index >= 1 && index <= self.base.nb_p3d,
            "Standard_OutOfRange: AppDef_MultiPointConstraint::SetCurv"
        );
        self.tab_curv[index - 1] = curv;
    }

    /// OCCT Curv(Index) (cxx L273-280) — the normal vector at the point of
    /// range Index.
    pub fn curv(&self, index: usize) -> DVec3 {
        assert!(
            index >= 1 && index <= self.base.nb_p3d,
            "Standard_OutOfRange: AppDef_MultiPointConstraint::Curv"
        );
        self.tab_curv[index - 1]
    }

    /// OCCT SetCurv2d(Index, Curv2d) (cxx L282-293) — Index is the GLOBAL
    /// point index.
    pub fn set_curv2d(&mut self, index: usize, curv2d: DVec2) {
        if self.tab_curv2d.is_empty() {
            self.tab_curv2d = vec![DVec2::ZERO; self.base.p2d.len()];
        }
        let nbp = self.base.nb_p3d;
        assert!(
            index > nbp && index <= nbp + self.base.p2d.len(),
            "Standard_OutOfRange: AppDef_MultiPointConstraint::SetCurv2d"
        );
        self.tab_curv2d[index - nbp - 1] = curv2d;
    }

    /// OCCT Curv2d(Index) (cxx L295-302) — Index is the GLOBAL point index.
    pub fn curv2d(&self, index: usize) -> DVec2 {
        let nbp = self.base.nb_p3d;
        assert!(
            index > nbp && index <= nbp + self.base.p2d.len(),
            "Standard_OutOfRange: AppDef_MultiPointConstraint::Curv2d"
        );
        self.tab_curv2d[index - nbp - 1]
    }

    /// OCCT IsTangencyPoint() (cxx L304-307) — returns True if the MultiPoint
    /// has a tangency value.
    pub fn is_tangency_point(&self) -> bool {
        !(self.tab_tang.is_empty() && self.tab_tang2d.is_empty())
    }

    /// OCCT IsCurvaturePoint() (cxx L309-312) — returns True if the
    /// MultiPoint has a curvature value.
    pub fn is_curvature_point(&self) -> bool {
        !(self.tab_curv.is_empty() && self.tab_curv2d.is_empty())
    }

    /// OCCT Dump(o) (cxx L314-317).
    pub fn dump(&self) {
        println!("AppDef_MultiPointConstraint dump:");
    }
}

// ---------------------------------------------------------------------------
// AppDef_MultiLine
// ---------------------------------------------------------------------------

/// OCCT AppDef_MultiLine — the organized set of points used in the
/// approximations. A MultiLine is composed of n MultiPointConstraints
/// (AppDef_MultiLine.hxx L33-53).
#[derive(Debug, Clone, Default)]
pub struct MultiLine {
    /// OCCT tabMult (HArray1 of AppDef_MultiPointConstraint).
    tab_mult: Vec<MultiPointConstraint>,
}

impl MultiLine {
    /// OCCT AppDef_MultiLine() (cxx L23) — creates an undefined MultiLine.
    pub fn new() -> Self {
        MultiLine {
            tab_mult: Vec::new(),
        }
    }

    /// OCCT AppDef_MultiLine(NbMult) (cxx L25-33) — given the number NbMult
    /// of MultiPointConstraints, initializes all the fields. SetValue must be
    /// called for the values of the multipoint constraint to be taken into
    /// account. Standard_ConstructionError if NbMult < 0 (unrepresentable
    /// for usize — the panic arm is kept for the signed overflow case).
    pub fn new_nb_mult(nb_mult: usize) -> Self {
        MultiLine {
            tab_mult: vec![MultiPointConstraint::new(); nb_mult],
        }
    }

    /// OCCT AppDef_MultiLine(tabMultiP) (cxx L35-43).
    pub fn new_tab_multi_p(tab_multi_p: &[MultiPointConstraint]) -> Self {
        MultiLine {
            tab_mult: tab_multi_p.to_vec(),
        }
    }

    /// OCCT AppDef_MultiLine(tabP3d) (cxx L45-55) — a MultiLine with one
    /// line of 3d points without their tangencies.
    pub fn new_tab_p3d(tab_p3d: &[DVec3]) -> Self {
        let mut tab_mult = Vec::with_capacity(tab_p3d.len());
        for p in tab_p3d {
            let mut mp = MultiPointConstraint::new_nb(1, 0);
            mp.base.set_point(1, *p);
            tab_mult.push(mp);
        }
        MultiLine { tab_mult }
    }

    /// OCCT AppDef_MultiLine(tabP2d) (cxx L57-67) — a MultiLine with one
    /// line of 2d points without their tangencies.
    pub fn new_tab_p2d(tab_p2d: &[DVec2]) -> Self {
        let mut tab_mult = Vec::with_capacity(tab_p2d.len());
        for p in tab_p2d {
            let mut mp = MultiPointConstraint::new_nb(0, 1);
            mp.base.set_point2d(1, *p);
            tab_mult.push(mp);
        }
        MultiLine { tab_mult }
    }

    /// OCCT NbMultiPoints() (cxx L69-72) — the number of
    /// MultiPointConstraints of the MultiLine.
    pub fn nb_multi_points(&self) -> usize {
        self.tab_mult.len()
    }

    /// OCCT NbPoints() (cxx L74-77) — the number of Points from MultiPoints
    /// composing the MultiLine.
    pub fn nb_points(&self) -> usize {
        self.tab_mult[0].base.nb_points() + self.tab_mult[0].base.nb_points2d()
    }

    /// OCCT SetValue(Index, MPoint) (cxx L79-86) — sets the
    /// MultiPointConstraint of range Index (1-based). Standard_OutOfRange
    /// outside [1, NbMultiPoints].
    pub fn set_value(&mut self, index: usize, mpoint: &MultiPointConstraint) {
        assert!(
            index >= 1 && index <= self.tab_mult.len(),
            "Standard_OutOfRange: AppDef_MultiLine::SetValue"
        );
        self.tab_mult[index - 1] = mpoint.clone();
    }

    /// OCCT Value(Index) (cxx L88-95) — the MultiPointConstraint of range
    /// Index (1-based).
    pub fn value(&self, index: usize) -> &MultiPointConstraint {
        assert!(
            index >= 1 && index <= self.tab_mult.len(),
            "Standard_OutOfRange: AppDef_MultiLine::Value"
        );
        &self.tab_mult[index - 1]
    }

    /// OCCT Dump(o) (cxx L97-103).
    pub fn dump(&self) {
        println!("AppDef_MultiLine dump:");
        println!(
            "It contains {} MultiPointConstraint",
            self.tab_mult.len()
        );
    }
}
