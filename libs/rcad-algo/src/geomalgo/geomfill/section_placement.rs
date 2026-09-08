//! OCCT GeomFill_SectionPlacement (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SectionPlacement.hxx (members) + GeomFill_SectionPlacement.cxx
//! (whole file L54-983; the commented-out legacy Perform body L526-594 is
//! not compiled in the OCCT either).
//!
//! Architecture differences:
//! - `GeomAdaptor_Curve myAdpSection` maps to the rcad `Curve3` view plus
//!   the first/last parameter helpers.
//! - OCCT Perform overloads map to `perform_confusion` (Perform(Tol)),
//!   `perform_with_path` (Perform(Path, Tol)) and `perform_at`
//!   (Perform(ParamOnPath, Tol)).
//! - `Extrema_ExtPC myExt` maps to the kernel [`ExtPC`] (the OCCT
//!   Initialize+Perform pair becomes one construction at each Perform site;
//!   TrimmedSquareDistances is recomputed from the query point — the same
//!   two endpoint distances).
//! - GAP carriers: `ExtCCGap` (OCCT Extrema_ExtCC over two bounded curves —
//!   the general curve-curve extremum is not translated) and the shared
//!   [`IntCurveSurfaceHInter`] plane-curve intersection carrier (staged
//!   batch-1 GAP, see its header).
//! - `gp_Trsf` results map to `DAffine3` (rigid transforms here).

use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

use glam::{DAffine3, DVec3};

use rcad_kernel::base::bnd_lib::curve_bounding_box_range;
use rcad_kernel::base::extrema::ExtPC;
use rcad_kernel::base::geom_lib::axe_of_inertia;
use rcad_kernel::base::geom_lprop::CLProps;
use rcad_kernel::core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use rcad_kernel::geom::{transform_curve, Curve3, CurveEval, Surface3};
use rcad_kernel::math::gp::{Ax1, Ax2, Ax3};

use super::gp_mat::GpMat;
use super::int_curve_surface_h_inter::IntCurveSurfaceHInter;
use super::location_law::LocationLaw;

/// OCCT Precision::Infinite().
const INFINITE: f64 = 2.0e100;
/// OCCT Precision::Angular().
const ANGULAR: f64 = 1.0e-12;
/// OCCT gp::Resolution() (Precision::RealSmall()).
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// The OCCT knots/multiplicities arrays recovered from the rcad flat knot
/// vector (pure tool helper).
fn flat_knots_to_mults(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i == 0 || *k != knots[knots.len() - 1] {
            knots.push(*k);
            mults.push(1);
        } else {
            *mults.last_mut().expect("nonempty") += 1;
        }
    }
    (knots, mults)
}

/// OCCT GeomAdaptor_Curve::Resolution (GeomAdaptor_Curve.cxx L1116-1148) —
/// the Line/Circle/Ellipse arms; the BSpline/Bezier arm falls back to the
/// OCCT adaptor default (GAP: BSplCLib::Resolution is not translated).
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

/// OCCT static Tangente (L54-69) — the normalized D1 (or the first
/// non-vanishing higher derivative).  Architecture: the rcad CurveEval
/// carries derivatives up to the 3rd order; DN beyond the 3rd evaluates the
/// 3rd.
fn tangente(path: &Curve3, param: f64, p: &mut DVec3, tang: &mut DVec3) {
    *p = path.point_at(param);
    *tang = path.derivative_at(param);
    let mut norm = tang.length();

    let mut ii = 2;
    while ii < 12 && norm < CONFUSION {
        *tang = match ii {
            2 => path.derivative2_at(param),
            _ => path.derivative3_at(param),
        };
        norm = tang.length();
        ii += 1;
    }

    if norm > 100.0 * GP_RESOLUTION {
        *tang /= norm;
    }
}

/// OCCT static Penalite (L71-98).
fn penalite(angle: f64, dist: f64) -> f64 {
    let mut penal;

    if dist < 1.0 {
        penal = dist.sqrt();
    } else if dist < 2.0 {
        penal = dist.powi(2);
    } else {
        penal = dist + 2.0;
    }

    if angle > 1.0e-3 {
        penal += 1.0 / angle - 2.0 / PI;
    } else {
        penal += 1.0e3;
    }

    penal
}

/// OCCT gp_Vec::AngleWithRef (gp_Dir.cxx L55-84) — signed angle from V to
/// Other measured in the plane reference Vref (pure math helper).
fn gp_vec_angle_with_ref(v: DVec3, other: DVec3, vref: DVec3) -> f64 {
    let xyz = v.cross(other);
    let cosinus = v.dot(other);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT static EvalAngle (L100-109).
fn eval_angle(v1: DVec3, v2: DVec3) -> f64 {
    let mut angle = gp_vec_angle(v1, v2);
    if angle > PI / 2.0 {
        angle = PI - angle;
    }
    angle
}

/// OCCT gp_Vec::Angle (gp_XYZ::Angle) — pure math helper.
fn gp_vec_angle(v1: DVec3, v2: DVec3) -> f64 {
    let an_norm = v1.length();
    let a_no_norm = v2.length();
    let mut value = v1.dot(v2) / (an_norm * a_no_norm);
    if value > 1.0 {
        value = 1.0;
    } else if value < -1.0 {
        value = -1.0;
    }
    value.acos()
}

/// OCCT static DistMini (L115-152) — examine an extrema to update <Dist> &
/// <Param>.  Architecture: the OCCT TrimmedSquareDistances values are
/// recomputed from the query point (the same two endpoint square distances).
fn dist_mini(ext: &ExtPC, c: &Curve3, query: DVec3, dist: &mut f64, param: &mut f64) {
    let mut dist2_var = f64::MAX;

    // OCCT: Ext.TrimmedSquareDistances(dist1, dist2, P1, P2).
    let dist1 = query.distance_squared(c.point_at(curve_first_parameter(c)));
    let dist2 = query.distance_squared(c.point_at(curve_last_parameter(c)));
    if (dist1 < dist2_var) || (dist2 < dist2_var) {
        if dist1 < dist2 {
            dist2_var = dist1;
            *param = curve_first_parameter(c);
        } else {
            dist2_var = dist2;
            *param = curve_last_parameter(c);
        }
    }

    if ext.is_done() {
        for ii in 1..=ext.nb_ext() {
            if ext.square_distance(ii) < dist2_var {
                dist2_var = ext.square_distance(ii);
                *param = ext.point(ii).param;
            }
        }
    }
    *dist = dist2_var.sqrt();
}

/// OCCT Geom_BSplineCurve::LocateU(U, Tolerance, I1, I2) over the rcad knot
/// arrays (1-based knot indexes; I1 == I2 when U is on a knot within Tol).
fn locate_u(knots: &[f64], u: f64, tol: f64) -> (usize, usize) {
    let nb = knots.len();
    for (i, k) in knots.iter().enumerate() {
        if (u - k).abs() <= tol {
            return (i + 1, i + 1);
        }
    }
    if u <= knots[0] {
        return (1, 1);
    }
    if u >= knots[nb - 1] {
        return (nb, nb);
    }
    for i in 1..nb {
        if knots[i - 1] < u && u < knots[i] {
            return (i, i + 1);
        }
    }
    (1, nb)
}

/// OCCT curve_first/last parameter (Adaptor3d_Curve over the rcad Curve3).
fn curve_first_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.first,
        other => other.default_domain()[0],
    }
}

fn curve_last_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.last,
        other => other.default_domain()[1],
    }
}

/// GAP carrier: OCCT Extrema_ExtCC over two bounded curves (TKGeomBase/
/// Extrema) — the general curve-curve extremum is not translated; the
/// construction keeps the OCCT failure path.
struct ExtCCGap;

impl ExtCCGap {
    /// OCCT Extrema_ExtCC(C1, C2, U1, U2, V1, V2, Tol1, Tol2).
    #[allow(clippy::too_many_arguments)]
    fn new(
        _c1: &Curve3,
        _c2: &Curve3,
        _u1: f64,
        _u2: f64,
        _v1: f64,
        _v2: f64,
        _tol1: f64,
        _tol2: f64,
    ) -> Self {
        panic!("GAP: Extrema_ExtCC (TKGeomBase/Extrema) is not translated — see file header")
    }

    /// OCCT IsDone().
    fn is_done(&self) -> bool {
        panic!("GAP: Extrema_ExtCC (TKGeomBase/Extrema) is not translated — see file header")
    }

    /// OCCT IsParallel().
    fn is_parallel(&self) -> bool {
        panic!("GAP: Extrema_ExtCC (TKGeomBase/Extrema) is not translated — see file header")
    }

    /// OCCT NbExt().
    fn nb_ext(&self) -> usize {
        panic!("GAP: Extrema_ExtCC (TKGeomBase/Extrema) is not translated — see file header")
    }

    /// OCCT SquareDistance(ii).
    fn square_distance(&self, _ii: usize) -> f64 {
        panic!("GAP: Extrema_ExtCC (TKGeomBase/Extrema) is not translated — see file header")
    }

    /// OCCT Points(ii, P1, P2) — (param1, point1, param2, point2).
    #[allow(clippy::type_complexity)]
    fn points(&self, _ii: usize) -> (f64, DVec3, f64, DVec3) {
        panic!("GAP: Extrema_ExtCC (TKGeomBase/Extrema) is not translated — see file header")
    }
}

/// OCCT GeomFill_SectionPlacement (hxx private members L84-98).
pub struct SectionPlacement {
    /// OCCT bool done.
    done: bool,
    /// OCCT bool isplan.
    isplan: bool,
    /// OCCT gp_Ax1 TheAxe.
    the_axe: Ax1,
    /// OCCT double Gabarit.
    gabarit: f64,
    /// OCCT handle(GeomFill_LocationLaw) myLaw.
    my_law: Rc<RefCell<dyn LocationLaw>>,
    /// OCCT GeomAdaptor_Curve myAdpSection (the rcad Curve3 view).
    my_adp_section: Option<Curve3>,
    /// OCCT handle(Geom_Curve) mySection.
    my_section: Option<Curve3>,
    /// OCCT double SecParam.
    sec_param: f64,
    /// OCCT double PathParam.
    path_param: f64,
    /// OCCT double Dist.
    dist: f64,
    /// OCCT double AngleMax.
    angle_max: f64,
    /// OCCT Extrema_ExtPC myExt (the last Perform result).
    my_ext: Option<ExtPC>,
    /// OCCT bool myIsPoint.
    my_is_point: bool,
    /// OCCT gp_Pnt myPoint.
    my_point: DVec3,
}

impl SectionPlacement {
    /// OCCT GeomFill_SectionPlacement(L, Section) (L156-380) — the section
    /// is a curve (the Geom_CartesianPoint form is carried by
    /// [`SectionPlacement::new_point`]).
    pub fn new_loc(l: Rc<RefCell<dyn LocationLaw>>, section: &Curve3) -> Self {
        let mut place = SectionPlacement {
            done: false,
            isplan: false,
            the_axe: Ax1::new(DVec3::ZERO, DVec3::Z),
            gabarit: 0.0,
            my_law: l,
            my_adp_section: Some(section.clone()),
            my_section: Some(section.clone()),
            sec_param: 0.0,
            path_param: 0.0,
            dist: f64::MAX, // RealLast()
            angle_max: 0.0,
            my_ext: None,
            my_is_point: false,
            my_point: DVec3::ZERO,
        };

        // Boite d'encombrement de la section pour en deduire le gabarit.
        // OCCT: BndLib_Add3dCurve::Add(myAdpSection, 1.e-4, box).
        let adp = place.my_adp_section.as_ref().expect("null myAdpSection");
        let (mut a_xmin, mut a_ymin, mut a_zmin, mut a_xmax, mut a_ymax, mut a_zmax) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        match curve_bounding_box_range(
            adp,
            curve_first_parameter(adp),
            curve_last_parameter(adp),
            1.0e-4,
        ) {
            Some([pmin, pmax]) => {
                a_xmin = pmin.x;
                a_ymin = pmin.y;
                a_zmin = pmin.z;
                a_xmax = pmax.x;
                a_ymax = pmax.y;
                a_zmax = pmax.z;
            }
            None => {
                // Degenerate bounding box — the whole box collapses on the
                // first point.
                let p = adp.point_at(curve_first_parameter(adp));
                a_xmin = p.x;
                a_ymin = p.y;
                a_zmin = p.z;
                a_xmax = p.x;
                a_ymax = p.y;
                a_zmax = p.z;
            }
        }

        let dx = a_xmax - a_xmin;
        let dy = a_ymax - a_ymin;
        let dz = a_zmax - a_zmin;
        place.gabarit = (dx * dx + dy * dy + dz * dz).sqrt() / 2.0;

        place.gabarit += CONFUSION; // Cas des toute petite

        // Initialisation de TheAxe pour les cas singulier.
        let mut nb_poles = 21usize;
        {
            let mut p = DVec3::ZERO;
            let mut v = DVec3::ZERO;
            tangente(
                adp,
                (curve_first_parameter(adp) + curve_last_parameter(adp)) / 2.0,
                &mut p,
                &mut v,
            );
            place.the_axe = Ax1::new(p, v);

            // y a t'il un Plan moyen ?
            match adp {
                Curve3::Circle(c) => {
                    place.isplan = true;
                    place.the_axe = Ax1::new(c.center, c.normal);
                }
                Curve3::Ellipse(e) => {
                    place.isplan = true;
                    place.the_axe = Ax1::new(e.center, e.normal);
                }
                Curve3::Hyperbola(h) => {
                    place.isplan = true;
                    place.the_axe = Ax1::new(h.center, h.normal);
                }
                Curve3::Parabola(pa) => {
                    place.isplan = true;
                    place.the_axe = Ax1::new(pa.vertex, pa.normal);
                }
                Curve3::Line(_) => {
                    nb_poles = 0; // Pas de Plan !!
                }
                Curve3::BSpline(bs) => {
                    nb_poles = bs.control_points.len();
                }
                _ => {
                    nb_poles = 21;
                }
            }
        }

        if !place.isplan && nb_poles > 2 {
            // Calcul d'un plan moyen.
            let mut pnts: Vec<DVec3> = Vec::new();
            let mut first = curve_first_parameter(adp);
            let mut last = curve_last_parameter(adp);
            if adp.is_periodic() {
                // Correct boundaries to avoid mistake of LocateU.
                // OCCT: the trimmed basis curve recovery — the rcad form
                // reads the underlying curve domain through the same
                // helpers.
                let ufirst = curve_first_parameter(adp);
                let a_period = curve_last_parameter(adp) - curve_first_parameter(adp);
                let u1 = ufirst + ((first - ufirst) / a_period).floor() * a_period;
                let u2 = u1 + a_period;
                if (first - u1).abs() <= PCONFUSION {
                    first = u1;
                }
                if (last - u2).abs() <= PCONFUSION {
                    last = u2;
                }
            }
            let mut t;
            if matches!(adp, Curve3::BSpline(_)) {
                let bc = match adp {
                    Curve3::BSpline(bs) => bs.clone(),
                    _ => unreachable!(),
                };
                let (knots, _mults) = flat_knots_to_mults(&bc.knots);
                let tol = CONFUSION;
                let (i1, i2) = locate_u(&knots, first, tol);
                let (i3, i4) = locate_u(&knots, last, tol);
                let nb_knots = i3 as i64 - i2 as i64 + 1;

                let nb_local_pnts = 10usize;
                let mut nb_pnts = ((nb_knots - 1) as usize) * nb_local_pnts;
                if i1 != i2 {
                    nb_pnts += nb_local_pnts;
                }
                if i3 != i4 && first < knots[i3 - 1] {
                    nb_pnts += nb_local_pnts;
                }
                if !adp.is_closed() {
                    nb_pnts += 1;
                }
                pnts = vec![DVec3::ZERO; nb_pnts];
                let mut nb = 0usize;
                if i1 != i2 {
                    let locallast = if knots[i2 - 1] < last {
                        knots[i2 - 1]
                    } else {
                        last
                    };
                    let delta = (locallast - first) / nb_local_pnts as f64;
                    for j in 0..nb_local_pnts {
                        t = first + j as f64 * delta;
                        pnts[nb] = adp.point_at(t);
                        nb += 1;
                    }
                }
                for i in i2..i3 {
                    let mut t = knots[i - 1];
                    let delta = (knots[i] - t) / nb_local_pnts as f64;
                    for _j in 0..nb_local_pnts {
                        pnts[nb] = adp.point_at(t);
                        nb += 1;
                        t += delta;
                    }
                }
                if i3 != i4 && first < knots[i3 - 1] {
                    let mut t = knots[i3 - 1];
                    let delta = (last - t) / nb_local_pnts as f64;
                    for _j in 0..nb_local_pnts {
                        pnts[nb] = adp.point_at(t);
                        nb += 1;
                        t += delta;
                    }
                }
                if !adp.is_closed() {
                    pnts[nb] = adp.point_at(last);
                }
            } else {
                // other type
                let mut nb_pnts = nb_poles - 1;
                if !adp.is_closed() {
                    nb_pnts += 1;
                }
                pnts = vec![DVec3::ZERO; nb_pnts];
                let delta = (last - first) / (nb_poles - 1) as f64;
                for i in 0..nb_poles - 1 {
                    let t = first + i as f64 * delta;
                    pnts[i] = adp.point_at(t);
                }
                if !adp.is_closed() {
                    pnts[nb_pnts - 1] = adp.point_at(last);
                }
                let _ = t;
            }

            let mut issing = false;
            let mut axe = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
            axe_of_inertia(&pnts, &mut axe, &mut issing, CONFUSION);
            if !issing {
                place.isplan = true;
                place.the_axe = Ax1::new(axe.location, axe.direction);
            }
        }

        // OCCT: myExt.Initialize(myAdpSection, First, Last, Confusion) — the
        // rcad ExtPC is constructed at each Perform site (the same range).
        place.my_ext = None;

        place
    }

    /// OCCT GeomFill_SectionPlacement(L, Point) — the Geom_CartesianPoint
    /// section form (myIsPoint).
    pub fn new_point(l: Rc<RefCell<dyn LocationLaw>>, point: DVec3) -> Self {
        let mut place = SectionPlacement {
            done: false,
            isplan: true,
            the_axe: Ax1::new(DVec3::ZERO, DVec3::Z),
            gabarit: 0.0,
            my_law: l,
            my_adp_section: None,
            my_section: None,
            sec_param: 0.0,
            path_param: 0.0,
            dist: f64::MAX, // RealLast()
            angle_max: 0.0,
            my_ext: None,
            my_is_point: true,
            my_point: point,
        };
        // OCCT: box.Add(myPoint); Gabarit = diag / 2 + Confusion.
        place.gabarit = 0.0 + CONFUSION;
        place
    }

    /// OCCT SetLocation (L384-387).
    pub fn set_location(&mut self, l: Rc<RefCell<dyn LocationLaw>>) {
        self.my_law = l;
    }

    /// OCCT myExt.Perform(PonPath) over the stored section range.
    fn ext_perform(&mut self, p: DVec3) -> Option<ExtPC> {
        let section = self.my_adp_section.clone()?;
        Some(ExtPC::new(
            p,
            &section,
            CONFUSION,
            curve_first_parameter(&section),
            curve_last_parameter(&section),
        ))
    }

    /// OCCT Perform(Tol) (L391-396).
    pub fn perform_confusion(&mut self, tol: f64) {
        let path = self.my_law.borrow().get_curve();
        self.perform_with_path_impl(path, tol);
    }

    /// OCCT Perform(Path, Tol) (L400-705).
    pub fn perform_with_path(&mut self, path: Option<Curve3>, tol: f64) {
        self.perform_with_path_impl(path, tol);
    }

    /// OCCT Perform(ParamOnPath, Tol) (L711-754).
    pub fn perform_at(&mut self, param: f64, tol: f64) {
        self.done = true;
        let path = self.my_law.borrow().get_curve();
        let path = match path {
            Some(p) => p,
            None => return,
        };

        self.path_param = param;
        if self.my_is_point {
            let pon_path = path.point_at(self.path_param);
            self.dist = pon_path.distance(self.my_point);
            self.angle_max = PI / 2.0;
        } else {
            self.sec_param = curve_first_parameter(self.my_adp_section.as_ref().expect("null section"));

            let mut pon_path;
            let mut pon_sec;
            let v_ref;
            let mut dp1 = DVec3::ZERO;
            v_ref = self.the_axe.direction;

            let mut p = DVec3::ZERO;
            tangente(&path, self.path_param, &mut p, &mut dp1);
            pon_path = p;
            let adp = self.my_adp_section.as_ref().expect("null section").clone();
            let adp = &adp;
            pon_sec = adp.point_at(self.sec_param);
            self.dist = pon_path.distance(pon_sec);
            if self.dist > tol {
                // On Cherche un meilleur point sur la section.
                if let Some(ext) = self.ext_perform(pon_path) {
                    if ext.is_done() {
                        let mut d = self.dist;
                        let mut sp = self.sec_param;
                        dist_mini(&ext, adp, pon_path, &mut d, &mut sp);
                        self.dist = d;
                        self.sec_param = sp;
                        pon_sec = adp.point_at(self.sec_param);
                    }
                }
            }
            self.angle_max = eval_angle(v_ref, dp1);
            if self.isplan {
                self.angle_max = PI / 2.0 - self.angle_max;
            }
            let _ = (pon_path, pon_sec);
        }

        self.done = true;
    }

    /// The shared Perform(Path, Tol) body (L400-705).
    fn perform_with_path_impl(&mut self, path: Option<Curve3>, tol: f64) {
        let path = match path {
            Some(p) => p,
            None => return,
        };
        let int_tol = 1.0e-5;
        let mut dist_center = INFINITE;

        if self.my_is_point {
            let section = path.clone();
            let projector = ExtPC::new(
                self.my_point,
                &section,
                CONFUSION,
                curve_first_parameter(&section),
                curve_last_parameter(&section),
            );
            let mut d = self.dist;
            let mut pp = self.path_param;
            dist_mini(&projector, &path, self.my_point, &mut d, &mut pp);
            self.dist = d;
            self.path_param = pp;
            self.angle_max = PI / 2.0;
        } else {
            self.path_param = curve_first_parameter(&path);
            self.sec_param =
                curve_first_parameter(self.my_adp_section.as_ref().expect("null section"));

            let mut distaux = 0.0;
            let mut taux = 0.0;
            let mut pon_path = DVec3::ZERO;
            let mut pon_sec;
            let v_ref = self.the_axe.direction;
            let mut dp1 = DVec3::ZERO;

            tangente(&path, self.path_param, &mut pon_path, &mut dp1);
            let adp = self.my_adp_section.as_ref().expect("null section").clone();
            let adp = &adp;
            pon_sec = adp.point_at(self.sec_param);
            self.dist = pon_path.distance(pon_sec);
            if self.dist > tol {
                // On Cherche un meilleur point sur la section.
                if let Some(ext) = self.ext_perform(pon_path) {
                    if ext.is_done() {
                        let mut d = self.dist;
                        let mut sp = self.sec_param;
                        dist_mini(&ext, adp, pon_path, &mut d, &mut sp);
                        self.dist = d;
                        self.sec_param = sp;
                        pon_sec = adp.point_at(self.sec_param);
                    }
                }
            }
            self.angle_max = eval_angle(v_ref, dp1);
            if self.isplan {
                self.angle_max = PI / 2.0 - self.angle_max;
            }

            let mut trouve = false;

            if self.isplan {
                // (1.1) Distances Point-Plan.
                let v1 = self.the_axe.location - pon_path;
                let dist_plan = v1.dot(v_ref).abs();
                if dist_plan <= int_tol {
                    dist_center = v1.length();
                }

                let plast = path.point_at(curve_last_parameter(&path));
                let v1 = self.the_axe.location - plast;
                let dist_plan = v1.dot(v_ref).abs();
                if dist_plan <= int_tol {
                    let a_dist = v1.length();
                    if a_dist < dist_center {
                        dist_center = a_dist;
                        pon_path = plast;
                        self.path_param = curve_last_parameter(&path);
                    }
                }

                // (1.2) Intersection Plan-courbe.
                // OCCT: gp_Ax3 axe(TheAxe.Location(), TheAxe.Direction());
                // plan = new Geom_Plane(axe).
                let axe = Ax3::from_pnt_n_vx(
                    self.the_axe.location,
                    self.the_axe.direction.normalize_or_zero(),
                    DVec3::X,
                );
                let plan = Surface3::Plane(rcad_kernel::geom::Plane {
                    origin: axe.axis.location,
                    normal: axe.direction(),
                    u_dir: axe.x_direction,
                    v_dir: axe.y_direction,
                });
                let mut intersector = IntCurveSurfaceHInter::default();
                intersector.perform(&path, &plan);
                let intersector_done = intersector.is_done();
                if intersector_done {
                    for ii in 1..=intersector.nb_points() {
                        let w = intersector.point(ii).w();
                        let p = path.point_at(w);
                        let a_dist = p.distance(self.the_axe.location);
                        if a_dist < dist_center {
                            dist_center = a_dist;
                            pon_path = p;
                            self.path_param = w;
                        }
                    }
                }
                if !intersector_done || intersector.nb_points() == 0 {
                    // Comparing the distances from the path's endpoints to
                    // the best matching plane of the profile.
                    let first_point = path.point_at(curve_first_parameter(&path));
                    let last_point = path.point_at(curve_last_parameter(&path));
                    let plane_origin = axe.axis.location;
                    let plane_normal = axe.direction();
                    let first_distance = (first_point - plane_origin).dot(plane_normal).powi(2);
                    let last_distance = (last_point - plane_origin).dot(plane_normal).powi(2);

                    if (first_distance.abs() < SQUARE_CONFUSION
                        && last_distance.abs() < SQUARE_CONFUSION)
                        || first_distance < last_distance
                    {
                        self.path_param = curve_first_parameter(&path);
                    } else {
                        self.path_param = curve_last_parameter(&path);
                        tangente(&path, self.path_param, &mut pon_path, &mut dp1);
                        pon_sec = adp.point_at(self.sec_param);
                        self.dist = pon_path.distance(pon_sec);
                        if self.dist > tol {
                            // On Cherche un meilleur point sur la section.
                            if let Some(ext) = self.ext_perform(pon_path) {
                                if ext.is_done() {
                                    let mut d = self.dist;
                                    let mut sp = self.sec_param;
                                    dist_mini(&ext, adp, pon_path, &mut d, &mut sp);
                                    self.dist = d;
                                    self.sec_param = sp;
                                    pon_sec = adp.point_at(self.sec_param);
                                }
                            }
                        }
                        self.angle_max = eval_angle(v_ref, dp1);
                        self.angle_max = PI / 2.0 - self.angle_max;
                    }
                }
            }

            // Cas General.
            if !self.isplan {
                // (2.1) Distance avec les extremites ...
                if let Some(ext) = self.ext_perform(pon_path) {
                    if ext.is_done() {
                        let mut d = self.dist;
                        let mut sp = self.sec_param;
                        dist_mini(&ext, adp, pon_path, &mut d, &mut sp);
                        distaux = d;
                        taux = sp;
                        if distaux < self.dist {
                            self.dist = distaux;
                            self.sec_param = taux;
                        }
                    }
                }
                trouve = self.dist <= tol;
                if !trouve {
                    let mut plast = DVec3::ZERO;
                    tangente(&path, curve_last_parameter(&path), &mut plast, &mut dp1);
                    let alpha = eval_angle(v_ref, dp1);
                    if let Some(ext) = self.ext_perform(plast) {
                        if ext.is_done() {
                            if ext.is_done() {
                                let mut d = self.dist;
                                let mut sp = self.sec_param;
                                dist_mini(&ext, adp, plast, &mut d, &mut sp);
                                distaux = d;
                                taux = sp;
                                if self.choix(distaux, alpha) {
                                    self.dist = distaux;
                                    self.sec_param = taux;
                                    self.angle_max = alpha;
                                    pon_path = plast;
                                    self.path_param = curve_last_parameter(&path);
                                }
                            }
                        }
                    }
                    trouve = self.dist <= tol;
                }

                // (2.2) Distance courbe-courbe.
                if !trouve {
                    let ext = ExtCCGap::new(
                        &path,
                        adp,
                        curve_first_parameter(&path),
                        curve_last_parameter(&path),
                        curve_first_parameter(adp),
                        curve_last_parameter(adp),
                        curve_resolution(&path, tol / 100.0),
                        curve_resolution(adp, tol / 100.0),
                    );
                    if ext.is_done() && !ext.is_parallel() {
                        for ii in 1..=ext.nb_ext() {
                            distaux = ext.square_distance(ii).sqrt();
                            let (p1_param, _p1_val, p2_param, p2_val) = ext.points(ii);
                            let mut pp = DVec3::ZERO;
                            tangente(&path, p1_param, &mut pp, &mut dp1);
                            let alpha = eval_angle(v_ref, dp1);
                            if self.choix(distaux, alpha) {
                                trouve = true;
                                self.dist = distaux;
                                self.path_param = p1_param;
                                self.sec_param = p2_param;
                                pon_sec = p2_val;
                                pon_path = pp;
                                self.angle_max = alpha;
                            }
                        }
                    }
                    if !trouve {
                        // Si l'on a toujours rien, on essai une distance
                        // point/path c'est la derniere chance.
                        // OCCT: PExt.Initialize(*Path, First, Last, Confusion)
                        // + Perform(PonSec).
                        if let Some(pext) = self.ext_perform(pon_sec) {
                            if pext.is_done() {
                                // modified for OCC13595: DistMini(PExt, *Path, ...).
                                let mut d = self.dist;
                                let mut sp = self.sec_param;
                                dist_mini(&pext, &path, pon_sec, &mut d, &mut sp);
                                distaux = d;
                                taux = sp;
                                let mut pp = DVec3::ZERO;
                                tangente(&path, taux, &mut pp, &mut dp1);
                                let alpha = eval_angle(v_ref, dp1);
                                if self.choix(distaux, alpha) {
                                    self.dist = distaux;
                                    pon_path = pp;
                                    self.angle_max = alpha;
                                    self.path_param = taux;
                                }
                            }
                        }
                    }
                }
            }
        }

        self.done = true;
    }

    /// OCCT IsDone (L758-761).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT ParameterOnPath (L765-768).
    pub fn parameter_on_path(&self) -> f64 {
        self.path_param
    }

    /// OCCT ParameterOnSection (L772-775).
    pub fn parameter_on_section(&self) -> f64 {
        self.sec_param
    }

    /// OCCT Distance (L779-782).
    pub fn distance(&self) -> f64 {
        self.dist
    }

    /// OCCT Angle (L786-789).
    pub fn angle(&self) -> f64 {
        self.angle_max
    }

    /// OCCT Transformation (L793-886).
    pub fn transformation(&self, with_translation: bool, with_correction: bool) -> DAffine3 {
        let mut v = DVec3::ZERO;
        let mut m = GpMat::identity();
        let mut p = DVec3::ZERO;
        let mut p_section = DVec3::ZERO;

        // Calcul des reperes: myLaw->D0(PathParam, M, V).
        self.my_law.borrow().d0(self.path_param, &mut m, &mut v);

        p = v;
        let d = gp_mat_column(&m, 3);
        let dn = gp_mat_column(&m, 1);
        let paxe = Ax3::from_pnt_n_vx(p, d.normalize_or_zero(), dn.normalize_or_zero());

        if with_translation || with_correction {
            if self.my_is_point {
                p_section = self.my_point;
            } else {
                p_section = self
                    .my_section
                    .as_ref()
                    .expect("null mySection")
                    .point_at(self.sec_param);
            }
        }

        // OCCT: gp_Trsf Rot.
        let mut rot = DAffine3::IDENTITY;

        if with_correction && !self.my_is_point {
            if !self.isplan {
                panic!("Illegal usage: can't rotate non-planar profile");
            }

            let profile_normal = self.the_axe.direction.normalize_or_zero();
            let spine_start_dir = paxe.direction();
            if !(profile_normal.cross(spine_start_dir).length()
                <= ANGULAR * profile_normal.length() * spine_start_dir.length())
            {
                let dir_axe_of_rotation = profile_normal.cross(spine_start_dir);
                let angle = gp_vec_angle_with_ref(
                    profile_normal,
                    spine_start_dir,
                    dir_axe_of_rotation,
                );
                // OCCT: Rot.SetRotation(AxeOfRotation, angle) — the rotation
                // about the axis through TheAxe.Location.
                rot = trsf_rotation_ax1(
                    self.the_axe.location,
                    dir_axe_of_rotation,
                    angle,
                );
            }
            p_section = rot.transform_point3(p_section);
        }

        if with_translation {
            // OCCT: P.ChangeCoord().SetLinearForm(-1, PSection.XYZ(), V.XYZ()).
            p = v - p_section;
        } else {
            p = DVec3::ZERO;
        }

        let saxe = Ax3::from_pnt_n_vx(p, DVec3::Z, DVec3::X);

        // Transfo: Tf.SetTransformation(Saxe, Paxe).
        let mut tf = trsf_set_transformation_between(&saxe, &paxe);

        if with_correction {
            // OCCT: Tf *= Rot.
            tf = tf * rot;
        }

        tf
    }

    /// OCCT Section (L890-895).
    pub fn section(&self, with_translation: bool) -> Curve3 {
        let the_section = self.my_section.as_ref().expect("null mySection").clone();
        transform_curve(&the_section, &self.transformation(with_translation, false))
    }

    /// OCCT ModifiedSection (L899-904).
    pub fn modified_section(&self, with_translation: bool) -> Curve3 {
        let the_section = self.my_section.as_ref().expect("null mySection").clone();
        transform_curve(&the_section, &self.transformation(with_translation, true))
    }

    /// OCCT SectionAxis (L908-948).
    fn section_axis(&self, m: &GpMat, t: &mut DVec3, n: &mut DVec3, bn: &mut DVec3) {
        let eps = 1.0e-10;
        let mut path_normal;
        let my_section = self.my_section.as_ref().expect("null mySection");
        let mut cp = CLProps::with_param(my_section, self.sec_param, 2, eps);
        if cp.is_tangent_defined() {
            let d = cp.tangent().expect("tangent");
            *t = d.normalize_or_zero();
            if cp.curvature() > eps {
                if let Some(dn) = cp.normal() {
                    *n = dn;
                }
            } else {
                // Cas ambigu, on essai de recuperer la normal a la trajectoire.
                path_normal = gp_mat_column(m, 1);
                path_normal = path_normal.normalize_or_zero();
                *bn = t.cross(path_normal);
                if bn.length() > eps {
                    *bn = bn.normalize_or_zero();
                }
                *n = bn.cross(*t);
            }
        } else {
            // Cas indefinie, on prend le triedre complet sur la trajectoire.
            *t = gp_mat_column(m, 3);
            *n = gp_mat_column(m, 2);
        }
        *bn = t.cross(*n);
    }

    /// OCCT Choix (L955-983) — decide si le couple (dist, angle) est
    /// "meilleur" que le couple courant.
    fn choix(&self, dist: f64, angle: f64) -> bool {
        let evoldist = dist - self.dist;
        let evolangle = angle - self.angle_max;
        // (1) Si la gain en distance est > que le gabarit, on prend.
        if evoldist < -self.gabarit {
            return true;
        }

        //  (2) si l'ecart en distance est de l'ordre du gabarit.
        if evoldist.abs() < self.gabarit {
            //  (2.1) si le gain en angle est important on garde.
            if evolangle > 0.5 {
                return true;
            }
            //  (2.2) si la variation d'angle est moderee on evalue une
            //  fonction de penalite.
            if penalite(angle, dist / self.gabarit) < penalite(self.angle_max, self.dist / self.gabarit) {
                return true;
            }
        }

        false
    }
}

/// OCCT gp_Mat::Column(theCol) — pure math helper.
fn gp_mat_column(m: &GpMat, the_col: usize) -> DVec3 {
    DVec3::new(
        m.mat[0][the_col - 1],
        m.mat[1][the_col - 1],
        m.mat[2][the_col - 1],
    )
}

/// OCCT gp_Trsf::SetTransformation(FromT1, ToT2) — pure math re-host (the
/// function_guide.rs precedent).
fn trsf_set_transformation_between(from_t1: &Ax3, to_t2: &Ax3) -> DAffine3 {
    let frame = |a: &Ax3| -> DAffine3 {
        let mm = glam::DMat3::from_cols(a.x_direction, a.y_direction, a.axis.direction);
        let mut f = DAffine3::from_mat3(mm);
        f.translation = a.axis.location;
        f
    };
    frame(to_t2).inverse() * frame(from_t1)
}

/// OCCT gp_Trsf::SetRotation(Ax1, Angle) — pure math re-host (Rodrigues).
fn trsf_rotation_ax1(loc: DVec3, axis_dir: DVec3, angle: f64) -> DAffine3 {
    let n = axis_dir.normalize_or_zero();
    let c = angle.cos();
    let s = angle.sin();
    let rotate = |v: DVec3| -> DVec3 {
        v * c + n.cross(v) * s + n * (n.dot(v)) * (1.0 - c)
    };
    let mut f = DAffine3::from_mat3(glam::DMat3::from_cols(
        rotate(DVec3::X),
        rotate(DVec3::Y),
        rotate(DVec3::Z),
    ));
    f.translation = loc - rotate(loc);
    f
}
