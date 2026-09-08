//! OCCT GeomFill_LocationGuide (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_LocationGuide.hxx (members) + GeomFill_LocationGuide.cxx (whole
//! file L58-1471; the OCCT_DEBUG TraceRevol helper and the D1/D2 commented
//! rotation bodies are not compiled there either).
//!
//! Architecture differences:
//! - `handle(GeomFill_TrihedronWithGuide) myLaw` maps to
//!   `Box<dyn TrihedronWithGuide>` (the ctor consumes the law; Copy()
//!   re-creates one through `TrihedronWithGuide::copy_with_guide` — the
//!   `down_cast<GeomFill_TrihedronWithGuide>` view of myLaw->Copy()).
//! - `handle(GeomFill_SectionLaw) mySec` maps to `Rc<dyn SectionLaw>` (the
//!   OCCT handle shares the section with the callers of Set()).
//! - `myStatus` is written from the const D0/D1/D2 evaluators (the OCCT
//!   methods are non-const; the rcad trait takes &self) — a `Cell`.
//! - `Extrema_ExtCS(*myGuide, GArevol, Confusion, Confusion)` maps to the
//!   [`ExtremaGenExtCS`] search with the OCCT Extrema_ExtCS generic
//!   dispatch parameters (NbT = 12, NbU = NbV = 10 / 13 on the periodic
//!   surface axes).  The Line-specific BndLib bounding-box range clamp is
//!   not re-hosted — the guide is searched over its full range.
//! - `GeomAdaptor_Curve::Resolution` maps to [`curve_resolution`] (the
//!   Line/Circle/Ellipse arms; the BSpline arm falls back to the OCCT
//!   adaptor default — GAP: BSplCLib::Resolution, TKMath/BSplCLib, is not
//!   translated).
//! - `math_Vector TolRes/Inf/Sup/X/R` map to the local `[f64; 3]` members
//!   (3 variables / 3 equations, 1-based in OCCT).

use std::cell::{Cell, RefCell};
use std::f64::consts::PI;
use std::rc::Rc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::geom_lib::fuse_intervals;
use rcad_kernel::geom::{
    BSplineCurve3, Curve3, CurveEval, RevolutionSurface, Surface3, SurfaceEval, TrimmedCurve3,
};
use rcad_kernel::math::el::in_period;
use rcad_kernel::math::function_set_root::FunctionSetRoot;
use rcad_kernel::math::gp::{Ax1, Ax3};
use rcad_kernel::math::GeomAbsShape;

use crate::bop::int_tools::bean_face_intersector::BRepAdaptorCurve;
use crate::bop::int_tools::extrema_gen_ext_cs::ExtremaGenExtCS;

use super::function_guide::{
    curve_transformed_by_trsf, gp_trsf_set_transformation_between, FunctionGuide,
};
use super::frenet::{curve_intervals, curve_nb_intervals};
use super::gp_mat::GpMat;
use super::location_law::LocationLaw;
use super::section_law::SectionLaw;
use super::trihedron_law::{curve_first_parameter, curve_last_parameter, PipeError, TrihedronLaw};
use super::trihedron_with_guide::TrihedronWithGuide;

/// OCCT Precision::Confusion().
const CONFUSION: f64 = 1.0e-7;
/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-9;

/// OCCT static InGoodPeriod (GeomFill_LocationGuide.cxx L128-142).
fn in_good_period(prec: f64, period: f64, current: &mut f64) {
    let mut diff = *current - prec;
    let nb = (diff / period).trunc() as i32;
    *current -= nb as f64 * period;
    diff = *current - prec;
    if diff > period / 2.0 {
        *current -= period;
    } else if diff < -period / 2.0 {
        *current += period;
    }
}

/// OCCT gp_Mat::Column(theCol) (gp_Mat.hxx) — the column read back as an XYZ
/// (pure math helper).
fn gp_mat_column(m: &GpMat, the_col: usize) -> DVec3 {
    DVec3::new(
        m.mat[0][the_col - 1],
        m.mat[1][the_col - 1],
        m.mat[2][the_col - 1],
    )
}

/// OCCT Adaptor3d_Curve::Period — LastParameter - FirstParameter (pure math
/// helper; the guide is periodic at the call sites).
fn curve_period_of(c: &Curve3) -> f64 {
    curve_last_parameter(c) - curve_first_parameter(c)
}

/// OCCT GeomAdaptor_Curve::Resolution (GeomAdaptor_Curve.cxx L1116-1148) —
/// the Line/Circle/Ellipse arms; the BSpline/Bezier arm falls back to the
/// OCCT adaptor default (GAP: BSplCLib::Resolution, TKMath/BSplCLib, is not
/// translated).
fn curve_resolution(c: &Curve3, r3d: f64) -> f64 {
    match c {
        Curve3::Line(_) => r3d,
        Curve3::Circle(circle) => {
            let r = circle.radius;
            if r > r3d / 2.0 {
                2.0 * (r3d / (2.0 * r)).asin()
            } else {
                2.0 * PI
            }
        }
        Curve3::Ellipse(ellipse) => r3d / ellipse.major_radius,
        _ => r3d * 0.01, // Precision::Parametric(R3D).
    }
}

/// OCCT GeomFill_LocationGuide (GeomFill_LocationGuide.hxx L190-218).
pub struct LocationGuide {
    /// OCCT handle(GeomFill_TrihedronWithGuide) myLaw.
    my_law: Box<dyn TrihedronWithGuide>,
    /// OCCT handle(GeomFill_SectionLaw) mySec — shared and later mutated
    /// in place through the same handle (the RefCell carries SetInterval /
    /// SetTolerance from the outer consumers).
    my_sec: Option<Rc<RefCell<dyn SectionLaw>>>,
    /// OCCT handle(Adaptor3d_Curve) myCurve.
    my_curve: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myGuide.
    my_guide: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myTrimmed.
    my_trimmed: Option<Curve3>,
    /// OCCT int myNbPts.
    my_nb_pts: usize,
    /// OCCT bool rotation.
    rotation: bool,
    /// OCCT double OrigParam1.
    orig_param1: f64,
    /// OCCT double OrigParam2.
    orig_param2: f64,
    /// OCCT double Uf.
    uf: f64,
    /// OCCT double Ul.
    ul: f64,
    /// OCCT double myFirstS.
    my_first_s: f64,
    /// OCCT double myLastS.
    my_last_s: f64,
    /// OCCT double ratio.
    ratio: f64,
    /// OCCT bool WithTrans.
    with_trans: bool,
    /// OCCT gp_Mat Trans.
    trans: GpMat,
    /// OCCT math_Vector TolRes(1, 3).
    tol_res: [f64; 3],
    /// OCCT math_Vector Inf(1, 3).
    inf: [f64; 3],
    /// OCCT math_Vector Sup(1, 3).
    sup: [f64; 3],
    /// OCCT math_Vector X(1, 3) — the FunctionSetRoot seed, written from
    /// InitX inside the const D0 evaluators (Cell).
    x: Cell<[f64; 3]>,
    /// OCCT math_Vector R(1, 3) — the FunctionSetRoot root, written from the
    /// const D0 evaluators (Cell).
    r: Cell<[f64; 3]>,
    /// OCCT GeomFill_PipeError myStatus (written from the const evaluators).
    my_status: Cell<PipeError>,
    /// OCCT handle(NCollection_HArray2<gp_Pnt2d>) myPoles2d (1..2, 1..myNbPts).
    my_poles2d: Vec<DVec2>,
}

impl LocationGuide {
    /// OCCT GeomFill_LocationGuide::GeomFill_LocationGuide (L146-180).
    pub fn new(triedre: Box<dyn TrihedronWithGuide>) -> Self {
        // OCCT: TolRes(1, 3) with TolRes.Init(1.e-6).
        let tol_res = [1.0e-6; 3];
        // OCCT: myLaw = Triedre (loi de triedre).
        let my_law = triedre;
        // OCCT: mySec.Nullify() (loi de section); myCurve.Nullify();
        // myFirstS = myLastS = -505e77.
        let my_sec: Option<Rc<RefCell<dyn SectionLaw>>> = None;
        let my_curve = None;
        let my_first_s = -505e77;
        let my_last_s = -505e77;

        let my_nb_pts = 21; // nb points pour les calculs
        // OCCT: myGuide = myLaw->Guide() (courbe guide).
        let mut my_guide = my_law.guide().expect("null guide");
        if !my_guide.is_periodic() {
            let mut f = curve_first_parameter(&my_guide);
            let mut l = curve_last_parameter(&my_guide);
            let delta = (l - f) / 100.0;
            f -= delta;
            l += delta;
            // OCCT: myGuide = myGuide->Trim(f, l, delta * 1.e-7).
            my_guide = Curve3::Trimmed(TrimmedCurve3::new(my_guide, f, l));
        } // if

        let my_poles2d = vec![DVec2::ZERO; 2 * my_nb_pts];
        let rotation = false; // contact ou non
        let orig_param1 = 0.0; // param pour ACR quand trajectoire
        let orig_param2 = 1.0; // et guide pas meme sens de parcourt
        let trans = GpMat::identity();
        let with_trans = false;

        LocationGuide {
            my_law,
            my_sec,
            my_curve,
            my_guide: Some(my_guide),
            my_trimmed: None,
            my_nb_pts,
            rotation,
            orig_param1,
            orig_param2,
            uf: 0.0,
            ul: 0.0,
            my_first_s,
            my_last_s,
            ratio: 0.0,
            with_trans,
            trans,
            tol_res,
            inf: [0.0; 3],
            sup: [0.0; 3],
            x: Cell::new([0.0; 3]),
            r: Cell::new([0.0; 3]),
            my_status: Cell::new(PipeError::PipeOk),
            my_poles2d,
        }
    }

    /// OCCT myPoles2d->Value(Row, Col) — the (1..2, 1..myNbPts) array
    /// (1-based access helper).
    #[inline]
    fn poles2d_value(&self, row: usize, col: usize) -> DVec2 {
        self.my_poles2d[(row - 1) * self.my_nb_pts + (col - 1)]
    }

    /// OCCT SetRotation (L184-466).
    fn set_rotation(&mut self, prec_angle: f64, last_angle: &mut f64) {
        if self.my_curve.is_none() {
            panic!("GeomFill_LocationGuide::The path is not set !!");
        }

        // repere fixe: gp_Ax3 Rep(gp::Origin(), gp::DZ(), gp::DX()).
        let rep = Ax3::new();

        let mut t = DVec3::ZERO;
        let mut n = DVec3::ZERO;
        let mut b = DVec3::ZERO;
        let mut isconst = false;
        let mut israt = false;
        let mut old_angle = 0.0;
        let mut cur_angle = prec_angle;

        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        let my_sec = self.my_sec.as_ref().expect("null mySec").borrow();
        let f = curve_first_parameter(my_curve);
        let l = curve_last_parameter(my_curve);
        let uperiodic = my_sec.is_u_periodic();

        let my_guide = self.my_guide.as_ref().expect("null myGuide");
        let delta_g =
            (curve_last_parameter(my_guide) - curve_first_parameter(my_guide)) / 5.0;
        let mut my_section: Option<Curve3> = None;
        let tol = 1.0e-9;

        let mut nb_poles = 0usize;
        let mut nb_knots = 0usize;
        let mut deg = 0usize;
        my_sec.section_shape(&mut nb_poles, &mut nb_knots, &mut deg);

        let mut mults: Vec<i32> = Vec::new();
        let mut knots: Vec<f64> = Vec::new();
        let mut poles: Vec<DVec3> = Vec::new();
        let mut weights: Vec<f64> = Vec::new();

        // OCCT: if (mySec->IsConstant(Tol)) { mySection = ConstantSection();
        // Uf/Ul = mySection->FirstParameter()/LastParameter(); isconst =
        // true; } else { israt = mySec->IsRational(); Mults/Knots/Poles/
        // Weights; Uf = Knots(1); Ul = Knots(NbKnots); }.
        let mut tol_const = tol;
        if my_sec.is_constant(&mut tol_const) {
            let section = my_sec.constant_section();
            self.uf = curve_first_parameter(&section);
            self.ul = curve_last_parameter(&section);
            my_section = Some(section);
            isconst = true;
        } else {
            isconst = false;
            israt = my_sec.is_rational();
            mults = vec![0i32; nb_knots];
            my_sec.mults(&mut mults);
            knots = vec![0.0; nb_knots];
            my_sec.knots(&mut knots);
            poles = vec![DVec3::ZERO; nb_poles];
            weights = vec![1.0; nb_poles];
            self.uf = knots[0];
            self.ul = knots[nb_knots - 1];
        }

        // Bornes de calculs.
        let mut delta = curve_last_parameter(my_guide) - curve_first_parameter(my_guide);
        self.inf[0] = curve_first_parameter(my_guide) - delta / 10.0;
        self.sup[0] = curve_last_parameter(my_guide) + delta / 10.0;

        self.inf[1] = -PI;
        self.sup[1] = 3.0 * PI;

        delta = self.ul - self.uf;
        self.inf[2] = self.uf - delta / 10.0;
        self.sup[2] = self.ul + delta / 10.0;

        // JALONNEMENT.
        let mut u_period = 0.0;
        if uperiodic {
            u_period = self.ul - self.uf;
        }

        let mut u = 0.0;
        for ii in 1..=self.my_nb_pts {
            let mut tt = (self.my_nb_pts - ii) as f64 * f + (ii - 1) as f64 * l;
            tt /= (self.my_nb_pts - 1) as f64;
            // myCurve->D0(tt, P).
            let p = my_curve.point_at(tt);
            // myLaw->D0(tt, T, N, B).
            let ok = TrihedronLaw::d0(self.my_law.as_ref(), tt, &mut t, &mut n, &mut b);
            if !ok {
                self.my_status.set(self.my_law.error_status());
                return; // Y a rien a faire.
            }
            // gp_Dir D = T; if (WithTrans) { gp_Mat M(N, B, T); M *= Trans;
            // D = M.Column(3); }.
            let mut d = t.normalize_or_zero();
            if self.with_trans {
                let mut m = GpMat::identity();
                m.set_cols(n, b, t);
                m = m.multiplied_mat(&self.trans);
                d = gp_mat_column(&m, 3).normalize_or_zero();
            }
            // gp_Ax1 Ax(P, D) — axe pour la surface de revolution.
            let ax = Ax1::new(p, d);

            // calculer transfo entre triedre et Oxyz:
            // gp_Dir N2 = N; gp_Ax3 N3(P, D, N2); gp_Trsf Transfo;
            // Transfo.SetTransformation(N3, Rep).
            let n3 = Ax3::from_pnt_n_vx(p, d, n.normalize_or_zero());
            let transfo = gp_trsf_set_transformation_between(&n3, &rep);

            // transformer la section.
            let s_curve: Curve3;
            if !isconst {
                u = self.my_first_s + (tt - curve_first_parameter(my_curve)) * self.ratio;
                my_sec.d0(u, &mut poles, &mut weights);
                let mut bs =
                    BSplineCurve3::from_knots_mults(deg, knots.clone(), mults.clone(), poles.clone());
                if israt {
                    bs.weights = weights.clone();
                }
                bs.is_periodic = my_sec.is_u_periodic();
                // S = new Geom_TrimmedCurve(mySection, Uf, Ul).
                s_curve = Curve3::Trimmed(TrimmedCurve3::new(
                    Curve3::BSpline(bs),
                    self.uf,
                    self.ul,
                ));
            } else {
                // S = new Geom_TrimmedCurve(mySection->Copy(), Uf, Ul).
                let section = my_section.as_ref().expect("null mySection");
                s_curve =
                    Curve3::Trimmed(TrimmedCurve3::new(section.clone(), self.uf, self.ul));
            }
            // S->Transform(Transfo).
            let s_transformed = curve_transformed_by_trsf(&s_curve, &transfo);

            // Surface de revolution: Revol = new
            // Geom_SurfaceOfRevolution(S, Ax); Extrema_ExtCS DistMini(
            // *myGuide, GArevol, Precision::Confusion(),
            // Precision::Confusion()) — the revolution surface dispatches to
            // Extrema_GenExtCS (Extrema_ExtCS.cxx L117-250): NbT = 12,
            // NbU = NbV = 10, bumped to 13 on the periodic surface axes.
            let the_u: f64;
            let the_v: f64;
            let mut pc_parameter = 0.0; // OCCT Extrema_POnCurv default state.
            {
                let revol = Surface3::Revolution(RevolutionSurface {
                    profile: Box::new(s_transformed),
                    axis_origin: ax.location,
                    axis_dir: ax.direction,
                });
                let revol_domain = revol.default_domain();
                let nb_u = if revol.is_u_periodic() { 13 } else { 10 };
                let nb_v = if revol.is_v_periodic() { 13 } else { 10 };
                let mut dist_mini = ExtremaGenExtCS::new();
                dist_mini.initialize(
                    &revol,
                    nb_u,
                    nb_v,
                    revol_domain[0],
                    revol_domain[1],
                    revol_domain[2],
                    revol_domain[3],
                    CONFUSION,
                );
                let mut adp_guide = BRepAdaptorCurve::new(my_guide.clone());
                dist_mini.perform(
                    &adp_guide,
                    12,
                    curve_first_parameter(my_guide),
                    curve_last_parameter(my_guide),
                    CONFUSION,
                );

                if !dist_mini.is_done() || dist_mini.nb_ext() == 0 {
                    let mut sos = false;
                    if ii > 1 {
                        // Intersection de secour entre surf revol et guide
                        // equation.
                        let mut x_seed = self.x.get();
                        x_seed[0] = self.poles2d_value(1, ii - 1).y;
                        x_seed[1] = self.poles2d_value(2, ii - 1).x;
                        x_seed[2] = self.poles2d_value(2, ii - 1).y;
                        self.x.set(x_seed);
                        let mut e = FunctionGuide::new(&*my_sec, my_guide, u);
                        e.set_param(u, p, t, n);
                        // resolution   =>  angle.
                        let mut result = FunctionSetRoot::new(&e, &self.tol_res, 100);
                        result.perform(&mut e, &x_seed, &self.inf, &self.sup, false);

                        if result.is_done() {
                            // OCCT: Result.FunctionSetErrors().Norm() <
                            // TolRes(1) * TolRes(1).
                            let errors = result.function_set_errors();
                            let norm = errors.iter().map(|v| v * v).sum::<f64>().sqrt();
                            if norm < self.tol_res[0] * self.tol_res[0] {
                                sos = true;
                                let rr = result.root();
                                // PInt.SetValues(P, RR(2), RR(3), RR(1),
                                // Out); theU = PInt.U(); theV = PInt.V().
                                the_u = rr[1];
                                the_v = rr[2];
                                self.r.set([rr[0], rr[1], rr[2]]);
                            } else {
                                the_u = 0.0;
                                the_v = 0.0;
                            }
                        } else {
                            the_u = 0.0;
                            the_v = 0.0;
                        }
                    } else {
                        the_u = 0.0;
                        the_v = 0.0;
                    }
                    if !sos {
                        self.my_status.set(PipeError::ImpossibleContact);
                        return;
                    }
                } else {
                    // on prend le point d'intersection d'angle le plus
                    // proche de P.
                    let mut min_dist = f64::MAX;
                    let mut jref = 0usize;
                    for j in 1..=dist_mini.nb_ext() {
                        let a_dist = dist_mini.square_distance(j);
                        if a_dist < min_dist {
                            min_dist = a_dist;
                            jref = j;
                        }
                    }
                    min_dist = min_dist.sqrt();
                    let _ = min_dist;
                    let (t_ext, u_ext, v_ext) = dist_mini.point(jref);
                    pc_parameter = t_ext;
                    the_u = u_ext;
                    the_v = v_ext;
                    let mut a1 = the_u;
                    in_good_period(cur_angle, 2.0 * PI, &mut a1);
                    // OCCT quirk: a1 is adjusted but never read afterwards.
                    let _ = a1;
                } // else
            }

            // Controle de w.
            let mut w = pc_parameter;
            if ii > 1 {
                let mut diff = w - self.poles2d_value(1, ii - 1).y;
                if diff.abs() > delta_g {
                    if my_guide.is_periodic() {
                        in_good_period(
                            self.poles2d_value(1, ii - 1).y,
                            curve_period_of(my_guide),
                            &mut w,
                        );
                        diff = w - self.poles2d_value(1, ii - 1).y;
                    }
                }
                let _ = diff;
            }

            // Recadrage de l'angle.
            let mut angle = the_u;
            if ii > 1 {
                let mut diff = angle - old_angle;
                if diff.abs() > PI {
                    in_good_period(old_angle, 2.0 * PI, &mut angle);
                    diff = angle - old_angle;
                }
                let _ = diff;
            }

            // Recadrage du V.
            let mut v = the_v;
            if ii > 1 {
                if uperiodic {
                    in_good_period(self.poles2d_value(2, ii - 1).y, u_period, &mut v);
                }
                let diff = v - self.poles2d_value(2, ii - 1).y;
                let _ = diff;
            }

            // on stocke les parametres.
            let p1 = DVec2::new(tt, w);
            let p2 = DVec2::new(angle, v);
            cur_angle = angle;
            self.my_poles2d[ii - 1] = p1;
            self.my_poles2d[self.my_nb_pts + ii - 1] = p2;
            old_angle = angle;
        }

        *last_angle = cur_angle;
        self.rotation = true; // C'est pret !
    }

    /// OCCT Set (L472-501) — init loi de section et force la Rotation.
    pub fn set(
        &mut self,
        section: Rc<RefCell<dyn SectionLaw>>,
        rotat: bool,
        s_first: f64,
        s_last: f64,
        prec_angle: f64,
        last_angle: &mut f64,
    ) {
        self.my_status.set(PipeError::PipeOk);
        self.my_first_s = s_first;
        self.my_last_s = s_last;
        *last_angle = prec_angle;
        match self.my_curve.as_ref() {
            None => self.ratio = 0.0,
            Some(my_curve) => {
                self.ratio = (s_last - s_first)
                    / (curve_last_parameter(my_curve) - curve_first_parameter(my_curve));
            }
        }
        self.my_sec = Some(section);

        if rotat {
            self.set_rotation(prec_angle, last_angle);
        } else {
            self.rotation = false;
        }
    }

    /// OCCT EraseRotation (L505-512).
    pub fn erase_rotation(&mut self) {
        self.rotation = false;
        if self.my_status.get() == PipeError::ImpossibleContact {
            self.my_status.set(PipeError::PipeOk);
        }
    }

    /// OCCT SetTrsf (L563-580).
    pub fn set_trsf(&mut self, transfo: GpMat) {
        self.trans = transfo;
        // OCCT: gp_Mat Aux; Aux.SetIdentity(); Aux -= Trans.
        let mut aux = GpMat::identity();
        aux = aux.added(&self.trans.multiplied_scalar(-1.0));
        self.with_trans = false; // Au cas ou Trans = I
        let mut ii = 1;
        while ii <= 3 && !self.with_trans {
            let mut jj = 1;
            while jj <= 3 && !self.with_trans {
                if aux.mat[ii - 1][jj - 1].abs() > 1.0e-14 {
                    self.with_trans = true;
                }
                jj += 1;
            }
            ii += 1;
        }
    }

    /// OCCT SetOrigine (L1430-1434) — utilise pour ACR dans le cas ou la
    /// trajectoire est multi-edges.
    pub fn set_origine(&mut self, param1: f64, param2: f64) {
        self.orig_param1 = param1;
        self.orig_param2 = param2;
    }

    /// OCCT Section (L1311-1314).
    pub fn section(&self) -> Curve3 {
        self.my_sec
            .as_ref()
            .expect("null mySec")
            .borrow()
            .constant_section()
    }

    /// OCCT Guide (L1318-1321).
    pub fn guide(&self) -> Option<Curve3> {
        self.my_guide.clone()
    }

    /// OCCT InitX (L1351-1424) — recherche par interpolation d'une valeur
    /// initiale (the OCCT method is non-const; the rcad form writes the
    /// X seed through the Cell).
    fn init_x(&self, param: f64) {
        let mut ideb = 1usize;
        let mut ifin = self.my_nb_pts; // myPoles2d->RowLength()
        let mut valeur;

        valeur = self.poles2d_value(1, ideb).x;
        if param == valeur {
            ifin = ideb + 1;
        }

        valeur = self.poles2d_value(1, ifin).x;
        if param == valeur {
            ideb = ifin - 1;
        }

        while ideb + 1 != ifin {
            let idemi = (ideb + ifin) / 2;
            valeur = self.poles2d_value(1, idemi).x;
            if valeur < param {
                ideb = idemi;
            } else {
                if valeur > param {
                    ifin = idemi;
                } else {
                    ideb = idemi;
                    ifin = ideb + 1;
                }
            }
        }

        let t1 = self.poles2d_value(1, ideb).x;
        let t2 = self.poles2d_value(1, ifin).x;
        let diff = t2 - t1;

        let w1 = self.poles2d_value(1, ideb).y;
        let w2 = self.poles2d_value(1, ifin).y;
        let p1 = self.poles2d_value(2, ideb);
        let p2 = self.poles2d_value(2, ifin);

        let mut x;
        if diff > 1.0e-7 {
            let b = (param - t1) / diff;
            let a = (t2 - param) / diff;
            x = [
                a * w1 + b * w2,
                a * p1.x + b * p2.x, // angle
                a * p1.y + b * p2.y, // param isov
            ];
        } else {
            x = [
                (w1 + w2) / 2.0,
                (p1.x + p2.x) / 2.0,
                (p1.y + p2.y) / 2.0,
            ];
        }

        let my_guide = self.my_guide.as_ref().expect("null myGuide");
        let my_sec = self.my_sec.as_ref().expect("null mySec").borrow();
        if my_guide.is_periodic() {
            x[0] = in_period(
                x[0],
                curve_first_parameter(my_guide),
                curve_last_parameter(my_guide),
            );
        }
        x[1] = in_period(x[1], 0.0, 2.0 * PI);
        if my_sec.is_u_periodic() {
            x[2] = in_period(x[2], self.uf, self.ul);
        }
        self.x.set(x);
    }

    /// OCCT ComputeAutomaticLaw (L1438-1471).
    pub fn compute_automatic_law(&self) -> (PipeError, Vec<DVec2>) {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        let f = curve_first_parameter(my_curve);
        let l = curve_last_parameter(my_curve);

        let mut par_and_rad = vec![DVec2::ZERO; self.my_nb_pts];
        for ii in 1..=self.my_nb_pts {
            let mut t = (self.my_nb_pts - ii) as f64 * f + (ii - 1) as f64 * l;
            t /= (self.my_nb_pts - 1) as f64;
            // myCurve->D0(t, P).
            let p = my_curve.point_at(t);
            let mut tv = DVec3::ZERO;
            let mut nv = DVec3::ZERO;
            let mut bv = DVec3::ZERO;
            let ok = TrihedronLaw::d0(self.my_law.as_ref(), t, &mut tv, &mut nv, &mut bv);
            if !ok {
                let the_status = self.my_law.error_status();
                return (the_status, par_and_rad);
            }
            let point_on_guide = self.my_law.current_point_on_guide();
            let cur_width = p.distance(point_on_guide);

            par_and_rad[ii - 1] = DVec2::new(t, cur_width);
        }

        (PipeError::PipeOk, par_and_rad)
    }

    /// OCCT D0 core (L584-650 / L656-728) — the shared rotation resolution;
    /// the two OCCT D0 overloads inline the same seed/resolution sequence
    /// (only the debug print differs).  The seed X and the R root are read
    /// and written through the Cells.
    fn d0_rotation(&self, param: f64, p: DVec3, m: &mut GpMat) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        let u = self.my_first_s + (param - curve_first_parameter(my_curve)) * self.ratio;
        // initialisations germe.
        self.init_x(param);
        let x_seed = self.x.get();

        let iter = 100;
        let t = gp_mat_column(m, 3);
        let b = gp_mat_column(m, 2);
        let n = gp_mat_column(m, 1);

        // Intersection entre surf revol et guide equation.
        let sec_guard = self.my_sec.as_ref().expect("null mySec").borrow();
        let mut e = FunctionGuide::new(
            &*sec_guard,
            self.my_guide.as_ref().expect("null myGuide"),
            u,
        );
        e.set_param(param, p, t, n);
        // resolution   =>  angle.
        let mut result = FunctionSetRoot::new(&e, &self.tol_res, iter);
        result.perform(&mut e, &x_seed, &self.inf, &self.sup, false);

        if result.is_done() {
            // solution.
            let r = result.root();
            self.r.set([r[0], r[1], r[2]]);

            // rotation: gp_Mat Rot; Rot.SetRotation(t, R(2)); b *= Rot;
            // n *= Rot.
            let mut rot = GpMat::identity();
            rot.set_rotation(t, r[1]);
            let b_rotated = GpMat::multiply_xyz_row(b, &rot);
            let n_rotated = GpMat::multiply_xyz_row(n, &rot);

            m.set_cols(n_rotated, b_rotated, t);
        } else {
            self.my_status.set(PipeError::ImpossibleContact);
            return false;
        }
        true
    }
}

impl LocationLaw for LocationGuide {
    /// OCCT Copy (L516-527).
    fn copy_law(&self) -> Box<dyn LocationLaw> {
        let mut la = 0.0f64;
        // OCCT: L = down_cast<GeomFill_TrihedronWithGuide>(myLaw->Copy()).
        let l = self.my_law.copy_with_guide();
        let mut copy = Box::new(LocationGuide::new(l));
        copy.set_origine(self.orig_param1, self.orig_param2);
        // OCCT: copy->Set(mySec, rotation, myFirstS, myLastS,
        // myPoles2d->Value(1, 1).X(), la).
        copy.set(
            self.my_sec.clone().expect("null mySec"),
            self.rotation,
            self.my_first_s,
            self.my_last_s,
            self.poles2d_value(1, 1).x,
            &mut la,
        );
        copy.set_trsf(self.trans);

        copy
    }

    /// OCCT SetCurve (L534-552).
    fn set_curve(&mut self, c: Curve3) -> bool {
        let mut last_angle = 0.0;
        self.my_curve = Some(c.clone());
        self.my_trimmed = Some(c.clone());

        if !self.my_curve.is_none() {
            self.my_law.set_curve(c);
            self.my_law.origine(self.orig_param1, self.orig_param2);
            self.my_status.set(self.my_law.error_status());

            if self.rotation {
                let seed = self.poles2d_value(1, 1).x;
                self.set_rotation(seed, &mut last_angle);
            }
        }
        self.my_status.get() == PipeError::PipeOk
    }

    /// OCCT GetCurve (L556-559).
    fn get_curve(&self) -> Option<Curve3> {
        self.my_curve.clone()
    }

    /// OCCT SetTrsf — the trait entry; the body lives in the inherent
    /// method above.
    fn set_trsf(&mut self, transfo: GpMat) {
        LocationGuide::set_trsf(self, transfo);
    }

    /// OCCT D0(Param, M, V) (L584-650).
    fn d0(&self, param: f64, m: &mut GpMat, v: &mut DVec3) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        // myCurve->D0(Param, P); V.SetXYZ(P.XYZ()).
        let p = my_curve.point_at(param);
        *v = p;
        // myLaw->D0(Param, T, N, B).
        let mut t = DVec3::ZERO;
        let mut n = DVec3::ZERO;
        let mut b = DVec3::ZERO;
        let ok = TrihedronLaw::d0(self.my_law.as_ref(), param, &mut t, &mut n, &mut b);
        if !ok {
            self.my_status.set(self.my_law.error_status());
            return ok;
        }
        m.set_cols(n, b, t);

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
        }

        if self.rotation {
            return self.d0_rotation(param, p, m);
        }

        true
    }

    /// OCCT D0(Param, M, V, Poles2d) (L656-728).
    fn d0_2d(&self, param: f64, m: &mut GpMat, v: &mut DVec3, _pnts2d: &mut [DVec2]) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        // myCurve->D0(Param, P); V.SetXYZ(P.XYZ()).
        let p = my_curve.point_at(param);
        *v = p;
        // myLaw->D0(Param, T, N, B).
        let mut t = DVec3::ZERO;
        let mut n = DVec3::ZERO;
        let mut b = DVec3::ZERO;
        let ok = TrihedronLaw::d0(self.my_law.as_ref(), param, &mut t, &mut n, &mut b);
        if !ok {
            self.my_status.set(self.my_law.error_status());
            return ok;
        }
        m.set_cols(n, b, t);

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
        }

        if self.rotation {
            // initialisation du germe + resolution (the shared rotation
            // core; the OCCT overload inlines U).
            return self.d0_rotation(param, p, m);
        }

        true
    }

    /// OCCT D1 (L734-890) — the OCCT rotation body is commented out and the
    /// literal early return stands.
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        m: &mut GpMat,
        v: &mut DVec3,
        dm: &mut GpMat,
        dv: &mut DVec3,
        _pnts2d: &mut [DVec2],
        _vecs2d: &mut [DVec2],
    ) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        // myCurve->D1(Param, P, DV).
        let p = my_curve.point_at(param);
        let d_p = my_curve.derivative_at(param);
        *v = p;
        *dv = d_p;
        // myLaw->D1(Param, T, DT, N, DN, B, DB).
        let (mut t, mut dt, mut n, mut dn, mut b, mut db) =
            (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let ok = self.my_law.d1(
            param,
            &mut t,
            &mut dt,
            &mut n,
            &mut dn,
            &mut b,
            &mut db,
        );
        if !ok {
            self.my_status.set(self.my_law.error_status());
            return ok;
        }
        m.set_cols(n, b, t);
        dm.set_cols(dn, db, dt);

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
            *dm = dm.multiplied_mat(&self.trans);
        }

        if self.rotation {
            // OCCT keeps the rotation body commented out — the literal
            // early return.
            return false;
        } // if_rotation

        true
    }

    /// OCCT D2 (L896-1145) — the OCCT rotation body is commented out and the
    /// literal early return stands; the WithTrans branch multiplies M/DM/D2M
    /// before they are filled (the C++ reads uninitialized caller storage —
    /// the rcad form keeps the literal sequence over the incoming values).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        m: &mut GpMat,
        v: &mut DVec3,
        dm: &mut GpMat,
        dv: &mut DVec3,
        d2m: &mut GpMat,
        d2v: &mut DVec3,
        _pnts2d: &mut [DVec2],
        _d1vecs2d: &mut [DVec2],
        _d2vecs2d: &mut [DVec2],
    ) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        // myCurve->D2(Param, P, DV, D2V).
        let p = my_curve.point_at(param);
        let d_p = my_curve.derivative_at(param);
        let d2_p = my_curve.derivative2_at(param);
        *v = p;
        *dv = d_p;
        *d2v = d2_p;
        // myLaw->D2(Param, T, DT, D2T, N, DN, D2N, B, DB, D2B).
        let (mut t, mut dt, mut d2t) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let (mut n, mut dn, mut d2n) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let (mut b, mut db, mut d2b) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let ok = TrihedronLaw::d2(
            self.my_law.as_ref(),
            param,
            &mut t,
            &mut dt,
            &mut d2t,
            &mut n,
            &mut dn,
            &mut d2n,
            &mut b,
            &mut db,
            &mut d2b,
        );
        if !ok {
            self.my_status.set(self.my_law.error_status());
            return ok;
        }

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
            *dm = dm.multiplied_mat(&self.trans);
            *d2m = d2m.multiplied_mat(&self.trans);
        }

        if self.rotation {
            // OCCT keeps the rotation body commented out — the literal
            // early return.
            return false;
        } else {
            m.set_cols(n, b, t);
            dm.set_cols(dn, db, dt);
            d2m.set_cols(d2n, d2b, d2t);
        }

        true
    }

    /// OCCT HasFirstRestriction (L1149-1152).
    fn has_first_restriction(&self) -> bool {
        false
    }

    /// OCCT HasLastRestriction (L1156-1159).
    fn has_last_restriction(&self) -> bool {
        false
    }

    /// OCCT TraceNumber (L1163-1166).
    fn trace_number(&self) -> usize {
        0
    }

    /// OCCT ErrorStatus (L1170-1173).
    fn error_status(&self) -> PipeError {
        self.my_status.get()
    }

    /// OCCT NbIntervals (L1177-1200).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        let nb_sec = curve_nb_intervals(my_trimmed, s);
        let nb_law = TrihedronLaw::nb_intervals(self.my_law.as_ref(), s);

        if nb_sec == 1 {
            return nb_law;
        } else if nb_law == 1 {
            return nb_sec;
        }

        let int_c = curve_intervals(my_trimmed, s);
        let mut int_l = vec![0.0; nb_law + 1];
        TrihedronLaw::intervals(self.my_law.as_ref(), &mut int_l, s);

        let mut inter: Vec<f64> = Vec::new();
        fuse_intervals(&int_c, &int_l, &mut inter, P_CONFUSION * 0.99, true);
        inter.len() - 1
    }

    /// OCCT Intervals (L1204-1232).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        let nb_sec = curve_nb_intervals(my_trimmed, s);
        let nb_law = TrihedronLaw::nb_intervals(self.my_law.as_ref(), s);

        if nb_sec == 1 {
            TrihedronLaw::intervals(self.my_law.as_ref(), t, s);
            return;
        } else if nb_law == 1 {
            let disc = curve_intervals(my_trimmed, s);
            t[..disc.len()].copy_from_slice(&disc);
            return;
        }

        let int_c = curve_intervals(my_trimmed, s);
        let mut int_l = vec![0.0; nb_law + 1];
        TrihedronLaw::intervals(self.my_law.as_ref(), &mut int_l, s);

        let mut inter: Vec<f64> = Vec::new();
        fuse_intervals(&int_c, &int_l, &mut inter, P_CONFUSION * 0.99, true);
        // OCCT: for (ii = 1; ii <= Inter.Length(); ii++) T(ii) = Inter(ii).
        for (ii, value) in inter.iter().enumerate() {
            t[ii] = *value;
        }
    }

    /// OCCT SetInterval (L1236-1240).
    fn set_interval(&mut self, first: f64, last: f64) {
        self.my_law.set_interval(first, last);
        // OCCT: myTrimmed = myCurve->Trim(First, Last, 0).
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.my_trimmed = Some(Curve3::Trimmed(TrimmedCurve3::new(
            my_curve.clone(),
            first,
            last,
        )));
    }

    /// OCCT GetInterval (L1244-1248).
    fn get_interval(&self, first: &mut f64, last: &mut f64) {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        *first = curve_first_parameter(my_trimmed);
        *last = curve_last_parameter(my_trimmed);
    }

    /// OCCT GetDomain (L1252-1256).
    fn get_domain(&self, first: &mut f64, last: &mut f64) {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        *first = curve_first_parameter(my_curve);
        *last = curve_last_parameter(my_curve);
    }

    /// OCCT SetTolerance (L1260-1264).
    fn set_tolerance(&mut self, tol3d: f64, _tol2d: f64) {
        // OCCT: TolRes(1) = myGuide->Resolution(Tol3d).
        self.tol_res[0] = curve_resolution(self.my_guide.as_ref().expect("null myGuide"), tol3d);
        // OCCT: Resolution(1, Tol3d, TolRes(2), TolRes(3)).
        let mut tolu = 0.0;
        let mut tolv = 0.0;
        LocationLaw::resolution(self, 1, tol3d, &mut tolu, &mut tolv);
        self.tol_res[1] = tolu;
        self.tol_res[2] = tolv;
    }

    /// OCCT Resolution (L1269-1276).
    fn resolution(&self, _index: usize, tol: f64, tolu: &mut f64, tolv: &mut f64) {
        *tolu = tol / 100.0;
        *tolv = tol / 100.0;
    }

    /// OCCT GetMaximalNorm (L1282-1285) — On suppose les triedres normes
    /// => return 1.
    fn get_maximal_norm(&self) -> f64 {
        1.0
    }

    /// OCCT GetAverageLaw (L1289-1307).
    fn get_average_law(&self, am: &mut GpMat, av: &mut DVec3) {
        let mut v1 = DVec3::ZERO;
        let mut v2 = DVec3::ZERO;
        let mut v3 = DVec3::ZERO;

        TrihedronLaw::get_average_law(self.my_law.as_ref(), &mut v1, &mut v2, &mut v3);
        am.set_cols(v1, v2, v3);

        *av = DVec3::ZERO;
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        let delta = (curve_last_parameter(my_trimmed) - curve_first_parameter(my_trimmed)) / 10.0;
        let mut u = curve_first_parameter(my_trimmed);
        for _ii in 0..=self.my_nb_pts {
            // OCCT: V.SetXYZ(myTrimmed->Value(U).XYZ()); AV += V.
            *av += my_trimmed.point_at(u);
            u += delta;
        }
        *av /= (self.my_nb_pts + 1) as f64;
    }

    /// OCCT IsTranslation (L1342-1345).
    fn is_translation(&self, _error: &mut f64) -> bool {
        false
    }

    /// OCCT IsRotation (L1326-1329).
    fn is_rotation(&self, _error: &mut f64) -> bool {
        false
    }

    /// OCCT Rotation (L1334-1337).
    fn rotation(&self, _centre: &mut DVec3) {
        panic!("Standard_NotImplemented: GeomFill_LocationGuide::Rotation");
    }
}
