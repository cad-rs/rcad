// OCCT LocOpe_CurveShapeIntersector.hxx L36-131 +
// LocOpe_CurveShapeIntersector.cxx L29-429 +
// LocOpe_CurveShapeIntersector.lxx L19-72 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_CurveShapeIntersector.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_CurveShapeIntersector.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_CurveShapeIntersector.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. BRepIntCurveSurface_Inter (TKTopAlgo) — the shape/curve iterator with
//    per-face IntCurvesFace_Intersector is not translated yet; the re-host
//    below (BRepIntCurveSurfaceInter) drives the same face loop over rcad's
//    reduced vehicle (topalgo::brep_int_curve_surface::Inter through the
//    loc_ope_cs_intersector::IntCurvesFaceIntersector re-host). The reduced
//    intersection semantics (UV-domain containment, In transition only) are
//    a GAP to close when the full IntCurvesFace translation lands.
// 2. OCCT GeomAdaptor_Curve over Geom_Circle maps to the rcad
//    geom::Circle3 carrier (same convention as fillet/chfi_ds.rs); the
//    [0, 2*PI] parameter bound becomes the Perform interval.
// 3. NCollection_Sequence<LocOpe_PntFace> is Vec<LocOpePntFace>; the OCCT
//    1-based sequence indexing becomes index-1 at every access site.
// 4. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate).
//
// first consumer: BRepFeat_MakeCylindricalHole (3a) — Perform(BRepFeat) L67,
// L495, L549, L693, L789, L897 construct LocOpe_CurveShapeIntersector(Axis,
// aObject) and read IsDone/NbPoints/LocalizeAfter/LocalizeBefore/Point;
// BRepFeat_Form family (3b).

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_cs_intersector::IntCurvesFaceIntersector;
use crate::feat::loc_ope_pnt_face::LocOpePntFace;
use crate::topalgo::brep_int_curve_surface::inter::TransitionOnCurve;
use rcad_kernel::geom::{ Circle3, Curve3 };
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType };

/// One iteration record of the re-hosted BRepIntCurveSurface_Inter
/// (architecture difference #1).
struct InterRecord {
    pnt: glam::DVec3,
    face: Shape,
    w: f64,
    u: f64,
    v: f64,
    transition: TransitionOnCurve,
}

/// OCCT BRepIntCurveSurface_Inter re-host (architecture difference #1) —
/// the OCCT surface used by LocOpe: Init(S, curve, Tol), More/Next, Pnt,
/// Face, W, U, V, Transition.
pub(crate) struct BRepIntCurveSurfaceInter {
    my_records: Vec<InterRecord>,
    my_index: usize,
}

impl BRepIntCurveSurfaceInter {
    pub(crate) fn new() -> Self {
        BRepIntCurveSurfaceInter {
            my_records: Vec::new(),
            my_index: 0,
        }
    }

    /// OCCT BRepIntCurveSurface_Inter::Init(S, curve, Tol) — explores the
    /// faces of S and intersects each one (the OCCT internal face loop).
    pub(crate) fn init(&mut self, the_s: &Shape, the_curve: &Curve3, the_tol: f64) {
        self.my_records.clear();
        self.my_index = 0;
        // OCCT BRepIntCurveSurface_Inter performs over the full parameter
        // interval (the interval overloads are not used here).
        let binf = f64::MIN; // RealFirst
        let bsup = f64::MAX; // RealLast
        for theface in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
            let mut the_int = IntCurvesFaceIntersector::new(&theface, the_tol);
            the_int.perform_curve(the_curve, binf, bsup);
            if !the_int.is_done() {
                continue;
            }
            for j in 1..=the_int.nb_pnt() {
                self.my_records.push(InterRecord {
                    pnt: the_int.pnt(j),
                    face: theface.clone(),
                    w: the_int.w_parameter(j),
                    u: the_int.u_parameter(j),
                    v: the_int.v_parameter(j),
                    transition: the_int.transition(j),
                });
            }
        }
    }

    /// OCCT BRepIntCurveSurface_Inter::More().
    pub(crate) fn more(&self) -> bool {
        self.my_index < self.my_records.len()
    }

    /// OCCT BRepIntCurveSurface_Inter::Next().
    pub(crate) fn next(&mut self) {
        self.my_index += 1;
    }

    /// OCCT BRepIntCurveSurface_Inter::Pnt().
    pub(crate) fn pnt(&self) -> glam::DVec3 {
        self.my_records[self.my_index].pnt
    }

    /// OCCT BRepIntCurveSurface_Inter::Face().
    pub(crate) fn face(&self) -> Shape {
        self.my_records[self.my_index].face.clone()
    }

    /// OCCT BRepIntCurveSurface_Inter::W().
    pub(crate) fn w(&self) -> f64 {
        self.my_records[self.my_index].w
    }

    /// OCCT BRepIntCurveSurface_Inter::U().
    pub(crate) fn u(&self) -> f64 {
        self.my_records[self.my_index].u
    }

    /// OCCT BRepIntCurveSurface_Inter::V().
    pub(crate) fn v(&self) -> f64 {
        self.my_records[self.my_index].v
    }

    /// OCCT BRepIntCurveSurface_Inter::Transition().
    pub(crate) fn transition(&self) -> TransitionOnCurve {
        self.my_records[self.my_index].transition
    }
}

/// OCCT LocOpe_CurveShapeIntersector (LocOpe_CurveShapeIntersector.hxx
/// L36-131).
pub struct LocOpeCurveShapeIntersector {
    my_done: bool,                    // OCCT: myDone
    my_points: Vec<LocOpePntFace>,    // OCCT: myPoints (arch. diff. #3)
}

impl Default for LocOpeCurveShapeIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeCurveShapeIntersector {
    /// OCCT LocOpe_CurveShapeIntersector::LocOpe_CurveShapeIntersector()
    /// (lxx L21-25) — empty constructor.
    pub fn new() -> Self {
        LocOpeCurveShapeIntersector {
            my_done: false,
            my_points: Vec::new(),
        }
    }

    /// OCCT LocOpe_CurveShapeIntersector::LocOpe_CurveShapeIntersector(Axis,
    /// S) (lxx L27-31) — Rust has no overloading: the `_axis` suffix.
    pub fn new_axis(the_axis: &Ax1, the_s: &Shape) -> Self {
        let mut res = LocOpeCurveShapeIntersector::new();
        res.init(the_axis, the_s);
        res
    }

    /// OCCT LocOpe_CurveShapeIntersector::LocOpe_CurveShapeIntersector(C, S)
    /// (lxx L33-37) — Rust has no overloading: the `_circ` suffix.
    pub fn new_circ(the_c: &Circle3, the_s: &Shape) -> Self {
        let mut res = LocOpeCurveShapeIntersector::new();
        res.init_circ(the_c, the_s);
        res
    }

    /// OCCT LocOpe_CurveShapeIntersector::Init(Axis, S) (cxx L33-47).
    pub fn init(&mut self, the_axis: &Ax1, the_s: &Shape) {
        self.my_done = false;
        self.my_points.clear();
        // OCCT cxx L37-40: if (S.IsNull()) return.
        if the_s.is_null() {
            return;
        }
        // OCCT cxx L41: Tol = Precision::Confusion().
        let tol = rcad_kernel::precision::CONFUSION;

        // OCCT cxx L43-45: theInt.Init(S, gp_Lin(Axis), Tol); Perform.
        let mut the_int = BRepIntCurveSurfaceInter::new();
        let the_lin = rcad_kernel::math::gp::Lin::from_ax1(the_axis);
        let curve = Curve3::Line(rcad_kernel::geom::Line3 {
            origin: the_lin.pos,
            direction: the_lin.dir,
        });
        the_int.init(the_s, &curve, tol);
        perform(&mut the_int, &mut self.my_points);
        // OCCT cxx L46.
        self.my_done = true;
    }

    /// OCCT LocOpe_CurveShapeIntersector::Init(C, S) (cxx L51-69).
    pub fn init_circ(&mut self, the_c: &Circle3, the_s: &Shape) {
        self.my_done = false;
        self.my_points.clear();
        // OCCT cxx L55-58: if (S.IsNull()) return.
        if the_s.is_null() {
            return;
        }
        // OCCT cxx L59: Tol = Precision::Confusion().
        let tol = rcad_kernel::precision::CONFUSION;

        // OCCT cxx L61-62: GC = Geom_Circle(C); AC = GeomAdaptor_Curve(GC,
        // 0., 2.*M_PI) (architecture difference #2).
        let curve = Curve3::Circle(*the_c);
        let binf = 0.0;
        let bsup = 2.0 * std::f64::consts::PI;

        // OCCT cxx L64-67: theInt.Init(S, AC, Tol); Perform.
        let mut the_int = BRepIntCurveSurfaceInter::new();
        the_int.init(the_s, &curve, tol);
        let _ = (binf, bsup); // the interval bounds live in the re-hosted
                              // perform loop (IntCurvesFaceIntersector
                              // performs over RealFirst..RealLast; the OCCT
                              // GeomAdaptor_Curve bounds the circle to
                              // [0, 2PI] which is the rcad circle domain).
        perform(&mut the_int, &mut self.my_points);
        // OCCT cxx L68.
        self.my_done = true;
    }

    /// OCCT LocOpe_CurveShapeIntersector::IsDone() (lxx L39-43).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_CurveShapeIntersector::NbPoints() (lxx L45-55).
    pub fn nb_points(&self) -> i32 {
        if !self.my_done {
            // OCCT lxx L47-50: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_points.len() as i32
    }

    /// OCCT LocOpe_CurveShapeIntersector::Point(Index) (lxx L57-72).
    pub fn point(&self, the_i: i32) -> LocOpePntFace {
        if !self.my_done {
            // OCCT lxx L59-62: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_points[(the_i - 1) as usize].clone()
    }

    /// OCCT LocOpe_CurveShapeIntersector::LocalizeAfter(From, Or, IndFrom,
    /// IndTo) (cxx L73-133) — the double From overload.
    pub fn localize_after(
        &self,
        the_from: f64,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L78-81: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        // OCCT cxx L82-84.
        let eps = rcad_kernel::precision::CONFUSION;
        let mut param;
        let fmeps = the_from - eps;
        let nbpoints = self.my_points.len() as i32;
        // OCCT cxx L85-91.
        let mut ifirst = 1;
        while ifirst <= nbpoints {
            if self.my_points[(ifirst - 1) as usize].parameter() >= fmeps {
                break;
            }
            ifirst += 1;
        }
        // OCCT cxx L92-130.
        let mut ret_val = false;
        if ifirst <= nbpoints {
            let mut i = ifirst;
            *the_ind_from = ifirst;
            let mut found = false;
            while !found {
                // OCCT cxx L100-102.
                *the_or = self.my_points[(i - 1) as usize].orientation();
                param = self.my_points[(i - 1) as usize].parameter();
                i += 1;
                // OCCT cxx L103-117.
                while i <= nbpoints {
                    if self.my_points[(i - 1) as usize].parameter() - param <= eps {
                        if *the_or != Orientation::External
                            && *the_or != self.my_points[(i - 1) as usize].orientation()
                        {
                            *the_or = Orientation::External;
                        }
                        i += 1;
                    } else {
                        break;
                    }
                }
                // OCCT cxx L118-128.
                if *the_or == Orientation::External {
                    found = i > nbpoints;
                    *the_ind_from = i;
                } else {
                    // on a une intersection franche
                    *the_ind_to = i - 1;
                    found = true;
                    ret_val = true;
                }
            }
        }

        ret_val
    }

    /// OCCT LocOpe_CurveShapeIntersector::LocalizeBefore(From, Or, IndFrom,
    /// IndTo) (cxx L137-197) — the double From overload.
    pub fn localize_before(
        &self,
        the_from: f64,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L142-145: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        // OCCT cxx L146-148.
        let eps = rcad_kernel::precision::CONFUSION;
        let mut param;
        let fpeps = the_from + eps;
        let nbpoints = self.my_points.len() as i32;
        // OCCT cxx L149-155.
        let mut ifirst = nbpoints;
        while ifirst >= 1 {
            if self.my_points[(ifirst - 1) as usize].parameter() <= fpeps {
                break;
            }
            ifirst -= 1;
        }
        // OCCT cxx L156-194.
        let mut ret_val = false;
        if ifirst >= 1 {
            let mut i = ifirst;
            *the_ind_to = ifirst;
            let mut found = false;
            while !found {
                // OCCT cxx L164-166.
                *the_or = self.my_points[(i - 1) as usize].orientation();
                param = self.my_points[(i - 1) as usize].parameter();
                i -= 1;
                // OCCT cxx L167-181.
                while i >= 1 {
                    if param - self.my_points[(i - 1) as usize].parameter() <= eps {
                        if *the_or != Orientation::External
                            && *the_or != self.my_points[(i - 1) as usize].orientation()
                        {
                            *the_or = Orientation::External;
                        }
                        i -= 1;
                    } else {
                        break;
                    }
                }
                // OCCT cxx L182-193.
                if *the_or == Orientation::External {
                    found = i < 1;
                    *the_ind_to = i;
                } else {
                    // on a une intersection franche
                    *the_ind_from = i + 1;
                    found = true;
                    ret_val = true;
                }
            }
        }

        ret_val
    }

    /// OCCT LocOpe_CurveShapeIntersector::LocalizeAfter(FromInd, Or, IndFrom,
    /// IndTo) (cxx L201-275) — the int FromInd overload.
    pub fn localize_after_index(
        &self,
        the_from_ind: i32,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L206-209: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        let nbpoints = self.my_points.len() as i32;
        // OCCT cxx L211-214.
        if the_from_ind >= nbpoints {
            return false;
        }

        // OCCT cxx L216-233.
        let eps = rcad_kernel::precision::CONFUSION;
        let mut param;
        let ifirst;
        if the_from_ind >= 1 {
            let fmeps = self.my_points[(the_from_ind - 1) as usize].parameter() - eps;
            let mut it = the_from_ind + 1;
            loop {
                if it > nbpoints {
                    break;
                }
                if self.my_points[(it - 1) as usize].parameter() >= fmeps {
                    break;
                }
                it += 1;
            }
            ifirst = it;
        } else {
            ifirst = 1;
        }

        // OCCT cxx L235-273.
        let mut ret_val = false;
        if ifirst <= nbpoints {
            let mut i = ifirst;
            *the_ind_from = ifirst;
            let mut found = false;
            while !found {
                // OCCT cxx L243-245.
                *the_or = self.my_points[(i - 1) as usize].orientation();
                param = self.my_points[(i - 1) as usize].parameter();
                i += 1;
                // OCCT cxx L246-260.
                while i <= nbpoints {
                    if self.my_points[(i - 1) as usize].parameter() - param <= eps {
                        if *the_or != Orientation::External
                            && *the_or != self.my_points[(i - 1) as usize].orientation()
                        {
                            *the_or = Orientation::External;
                        }
                        i += 1;
                    } else {
                        break;
                    }
                }
                // OCCT cxx L261-272.
                if *the_or == Orientation::External {
                    found = i > nbpoints;
                    *the_ind_from = i;
                } else {
                    // on a une intersection franche
                    *the_ind_to = i - 1;
                    found = true;
                    ret_val = true;
                }
            }
        }
        ret_val
    }

    /// OCCT LocOpe_CurveShapeIntersector::LocalizeBefore(FromInd, Or,
    /// IndFrom, IndTo) (cxx L279-353) — the int FromInd overload.
    pub fn localize_before_index(
        &self,
        the_from_ind: i32,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L284-287: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        let nbpoints = self.my_points.len() as i32;
        // OCCT cxx L289-292.
        if the_from_ind <= 1 {
            return false;
        }

        // OCCT cxx L294-311.
        let eps = rcad_kernel::precision::CONFUSION;
        let mut param;
        let ifirst;
        if the_from_ind <= nbpoints {
            let fpeps = self.my_points[(the_from_ind - 1) as usize].parameter() + eps;
            let mut it = the_from_ind - 1;
            loop {
                if it < 1 {
                    break;
                }
                if self.my_points[(it - 1) as usize].parameter() <= fpeps {
                    break;
                }
                it -= 1;
            }
            ifirst = it;
        } else {
            ifirst = nbpoints;
        }

        // OCCT cxx L313-351.
        let mut ret_val = false;
        if ifirst >= 1 {
            let mut i = ifirst;
            *the_ind_to = ifirst;
            let mut found = false;
            while !found {
                // OCCT cxx L321-323.
                *the_or = self.my_points[(i - 1) as usize].orientation();
                param = self.my_points[(i - 1) as usize].parameter();
                i -= 1;
                // OCCT cxx L324-338.
                while i >= 1 {
                    if param - self.my_points[(i - 1) as usize].parameter() <= eps {
                        if *the_or != Orientation::External
                            && *the_or != self.my_points[(i - 1) as usize].orientation()
                        {
                            *the_or = Orientation::External;
                        }
                        i -= 1;
                    } else {
                        break;
                    }
                }
                // OCCT cxx L339-350.
                if *the_or == Orientation::External {
                    found = i < 1;
                    *the_ind_to = i;
                } else {
                    // on a une intersection franche
                    *the_ind_from = i + 1;
                    found = true;
                    ret_val = true;
                }
            }
        }
        ret_val
    }
}

/// OCCT static Perform(BRepIntCurveSurface_Inter&, NCollection_Sequence<
/// LocOpe_PntFace>&) (cxx L357-429).
fn perform(the_int: &mut BRepIntCurveSurfaceInter, the_points: &mut Vec<LocOpePntFace>) {
    // OCCT cxx L360-363.
    let mut param;
    let mut paramu;
    let mut paramv;
    let mut nbpoints: i32 = 0;

    let mut theor;
    // OCCT cxx L365-428.
    while the_int.more() {
        // OCCT cxx L367-372.
        let thept = the_int.pnt();
        let theface = the_int.face();
        let orface = theface.orientation;
        param = the_int.w();
        paramu = the_int.u();
        paramv = the_int.v();

        // OCCT cxx L374-407: transition switch.
        theor = match the_int.transition() {
            TransitionOnCurve::In => {
                if orface == Orientation::Forward {
                    Orientation::Forward
                } else if orface == Orientation::Reversed {
                    Orientation::Reversed
                } else {
                    Orientation::External
                }
            }
            TransitionOnCurve::Out => {
                if orface == Orientation::Forward {
                    Orientation::Reversed
                } else if orface == Orientation::Reversed {
                    Orientation::Forward
                } else {
                    Orientation::External
                }
            }
            TransitionOnCurve::Tangent => Orientation::External,
        };

        // OCCT cxx L409.
        let newpt = LocOpePntFace::new_full(thept, &theface, theor, param, paramu, paramv);

        // OCCT cxx L411-417: insertion point search.
        let mut i = 1;
        while i <= nbpoints {
            if the_points[(i - 1) as usize].parameter() - param > 0.0 {
                break;
            }
            i += 1;
        }
        // OCCT cxx L418-425.
        if i <= nbpoints {
            the_points.insert((i - 1) as usize, newpt);
        } else {
            the_points.push(newpt);
        }
        nbpoints += 1;
        // OCCT cxx L427: theInt.Next().
        the_int.next();
    }
}
