//! OCCT Bisector_BisecAna — the analytic (line/circle/point) bisector,
//! 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_BisecAna.hxx (L40-183) / .cxx (L47-1826).
//!
//! The five GccAna bisector constructors and the GccInt results are GAP
//! carriers (see super::deps_gap).

use std::sync::{Arc, RwLock};

use glam::DVec2;
use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE};
use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
use rcad_kernel::math::GeomAbsShape;

use super::bisector::{cross2, dir2d_parallel, pnt2d_equal, rotate_half_pi, GeomAbsJoinType};
use super::bisector_curve::{
    BisectorCurve, CurveKind, Geom2dCurveAdaptor, Geom2dCurveHandle, ResD1, ResD2, ResD3,
    TrimmedCurve,
};
use crate::geomalgo::geom2d_int::elclib2d;
use super::deps_gap::{
    GccAnaCirc2dBisec, GccAnaCircLin2dBisec, GccAnaCircPnt2dBisec, GccAnaLinPnt2dBisec,
    GccAnaPnt2dBisec, GccIntBisec, GccIntIType,
};

/// OCCT M_PI.
const PI: f64 = std::f64::consts::PI;

/// OCCT Bisector_BisecAna (Bisector_BisecAna.hxx L40-183).
pub struct BisectorBisecAna {
    /// OCCT thebisector (Handle(Geom2d_TrimmedCurve)); null handle -> None.
    /// Interior mutability mirrors the OCCT handle mutation semantics
    /// (SetTrim/Reverse mutate the shared object).
    pub(crate) thebisector: RwLock<Option<TrimmedCurve>>,
}

impl Default for BisectorBisecAna {
    /// OCCT Bisector_BisecAna() (L54).
    fn default() -> Self {
        BisectorBisecAna { thebisector: RwLock::new(None) }
    }
}

// ---------------------------------------------------------------------------
// Kernel-curve helpers (the DynamicType / down_cast mirrors).  A
// `handle<Geom2d_Curve>` maps to `Arc<dyn BisectorCurve>`; a plain kernel
// curve enters as a Geom2dCurveHandle.
// ---------------------------------------------------------------------------

/// OCCT DynamicType comparisons on a kernel curve.
pub(crate) fn kind_of(curve: &Curve2d) -> CurveKind {
    match curve {
        Curve2d::Line(_) => CurveKind::Line,
        Curve2d::Circle(_) => CurveKind::Circle,
        Curve2d::Ellipse(_) => CurveKind::Ellipse,
        Curve2d::Parabola(_) => CurveKind::Parabola,
        Curve2d::Hyperbola(_) => CurveKind::Hyperbola,
        _ => CurveKind::Other,
    }
}

/// OCCT `CurveF = BasisCurve()` mirror —
/// `if (Type == STANDARD_TYPE(Geom2d_TrimmedCurve)) CurveF = BasisCurve();`.
pub(crate) fn unwrap_trimmed(curve: &Arc<dyn BisectorCurve>) -> Arc<dyn BisectorCurve> {
    match curve.as_any().downcast_ref::<TrimmedCurve>() {
        Some(t) => t.basis().clone(),
        None => curve.clone(),
    }
}

/// Extract the kernel curve payload of a Geom2dCurveHandle.
pub(crate) fn kernel_curve_of(curve: &Arc<dyn BisectorCurve>) -> &Curve2d {
    &downcast_handle(curve).curve
}

pub(crate) fn circle_of(curve: &Arc<dyn BisectorCurve>) -> Circle2d {
    match kernel_curve_of(curve) {
        Curve2d::Circle(c) => *c,
        _ => unreachable!("circle_of on non-circle"),
    }
}

pub(crate) fn line_of(curve: &Arc<dyn BisectorCurve>) -> Line2d {
    match kernel_curve_of(curve) {
        Curve2d::Line(l) => *l,
        _ => unreachable!("line_of on non-line"),
    }
}

pub(crate) fn geom_value(curve: &Curve2d, u: f64) -> DVec2 {
    crate::geomalgo::geom2d_int::geom2d_curve_tool::value(curve, u)
}

/// Kernel-level line extraction (a Geom2dCurveHandle payload).
pub(crate) fn line_of_kernel(curve: &Curve2d) -> Line2d {
    match curve {
        Curve2d::Line(l) => *l,
        _ => unreachable!("line_of on non-line"),
    }
}

/// Kernel-level circle extraction (a Geom2dCurveHandle payload).
pub(crate) fn circle_of_kernel(curve: &Curve2d) -> Circle2d {
    match curve {
        Curve2d::Circle(c) => *c,
        _ => unreachable!("circle_of on non-circle"),
    }
}

pub(crate) fn hyperbola_of(curve: &Curve2d) -> rcad_kernel::geom::Hyperbola2d {
    match curve {
        Curve2d::Hyperbola(h) => *h,
        _ => unreachable!("hyperbola_of on non-hyperbola"),
    }
}

pub(crate) fn parabola_of(curve: &Curve2d) -> rcad_kernel::geom::Parabola2d {
    match curve {
        Curve2d::Parabola(pp) => *pp,
        _ => unreachable!("parabola_of on non-parabola"),
    }
}

pub(crate) fn ellipse_of(curve: &Curve2d) -> rcad_kernel::geom::Ellipse2d {
    match curve {
        Curve2d::Ellipse(e) => *e,
        _ => unreachable!("ellipse_of on non-ellipse"),
    }
}

/// Down-cast a `Geom2dCurveHandle` out of the trait object.
pub(crate) fn downcast_handle(basis: &Arc<dyn BisectorCurve>) -> &Geom2dCurveHandle {
    basis
        .as_any()
        .downcast_ref::<Geom2dCurveHandle>()
        .expect("Geom2dCurveHandle downcast")
}

/// OCCT `down_cast<Bisector_BisecAna>(handle)` mirror.
pub(crate) fn downcast_to_ana(h: &Arc<dyn BisectorCurve>) -> &BisectorBisecAna {
    h.as_any()
        .downcast_ref::<BisectorBisecAna>()
        .expect("Bisector_BisecAna down-cast")
}

/// OCCT Geom2dAdaptor_Curve constructor mirror.
pub(crate) fn adaptor_of(curve: Arc<dyn BisectorCurve>) -> Geom2dCurveAdaptor {
    Geom2dCurveAdaptor::new(curve)
}

/// `down_cast<Bisector_BisecAna>(h)->Geom2dCurve()` — the kernel curve
/// payload of the analytic bisector, if any (reconciled mat2d surface).
pub fn geom2d_curve_of(h: &Arc<dyn BisectorCurve>) -> Option<Curve2d> {
    let bis = h.as_any().downcast_ref::<BisectorBisecAna>()?;
    let basis = bis.geom2d_curve();
    basis
        .as_any()
        .downcast_ref::<Geom2dCurveHandle>()
        .map(|handle| handle.curve.clone())
}

/// `down_cast<Bisector_BisecAna>(h)->SetTrim(uf, ul)` (reconciled mat2d
/// surface; mutation through the OCCT handle).
pub fn set_trim_of(h: &Arc<dyn BisectorCurve>, uf: f64, ul: f64) {
    if let Some(bis) = h.as_any().downcast_ref::<BisectorBisecAna>() {
        bis.set_trim(uf, ul);
    }
}

/// OCCT Degenerate(aBisector, Tolerance) (L1737-1803) — replaces a
/// degenerate hyperbola/parabola/ellipse by its axis line.
pub(crate) fn degenerate(a_bisector: &mut GccIntBisec, tolerance: f64) -> bool {
    let mut degeneree = false;

    let typ = a_bisector.arc_type();

    if typ == GccIntIType::Hpr {
        let gphyperbola = a_bisector.hyperbola();

        // If the Hyperbola is degenerated, it is replaced by the straight
        // line with direction to the axis if symmetry.
        if gphyperbola.semi_major < tolerance {
            // OCCT gp_Lin2d gpline(gphyperbola.YAxis()).
            let y_axis = DVec2::new(-gphyperbola.major_dir.y, gphyperbola.major_dir.x);
            let gpline = Line2d::new(gphyperbola.center, y_axis);
            *a_bisector = GccIntBisec::Line(gpline);
            degeneree = true;
        }
        if gphyperbola.semi_minor < tolerance {
            // OCCT gp_Lin2d gpline(gphyperbola.XAxis()).
            let gpline = Line2d::new(gphyperbola.center, gphyperbola.major_dir);
            *a_bisector = GccIntBisec::Line(gpline);
            degeneree = true;
        }
    } else if typ == GccIntIType::Par {
        let gpparabola = a_bisector.parabola();

        // If the parabole is degenerated, it is replaces by the straight
        // line starting at the Top and with direction of the axis of
        // symmetry.
        if gpparabola.focal_param < tolerance {
            // OCCT gp_Lin2d gpline(gpparabola.MirrorAxis()).
            let gpline = Line2d::new(gpparabola.origin, gpparabola.axis_dir);
            *a_bisector = GccIntBisec::Line(gpline);
            degeneree = true;
        }
    } else if typ == GccIntIType::Ell {
        let gpellipse = a_bisector.ellipse();

        // If the ellipse is degenerated, it is replaced by the straight line
        // defined by the great axis.
        if gpellipse.minor_radius < tolerance {
            // OCCT gp_Lin2d gpline(gpellipse.XAxis()).
            let gpline = Line2d::new(gpellipse.center, gpellipse.major_dir);
            *a_bisector = GccIntBisec::Line(gpline);
            degeneree = true;
        }
    }
    degeneree
}

impl BisectorBisecAna {
    /// OCCT Bisector_BisecAna() (L54).
    pub fn new() -> Self {
        BisectorBisecAna::default()
    }

    /// OCCT Distance (L68-220) — distance between the point and the
    /// bissectrice, and orientation of the bissectrice.
    #[allow(clippy::too_many_arguments)]
    fn distance(
        &self,
        apoint: DVec2,
        abisector: &GccIntBisec,
        afirstvector: DVec2,
        asecondvector: DVec2,
        vec_ref: DVec2,
        adirection: f64,
        aparameter: &mut f64,
        asense: &mut bool,
        astatus: &mut bool,
        is_bisec_of_two_lines: bool,
    ) -> f64 {
        *astatus = true;

        let mut point;
        let mut tangent;

        let typ = abisector.arc_type();

        if typ == GccIntIType::Lin {
            let gpline = abisector.line();
            *aparameter = elclib2d::line_parameter(gpline.origin, gpline.direction, apoint);
            let (p, t) = elclib2d::line_d1(gpline.origin, gpline.direction, *aparameter);
            point = p;
            tangent = t;
        } else if typ == GccIntIType::Cir {
            let gpcircle = abisector.circle();
            *aparameter = elclib2d::circle_parameter(
                gpcircle.center,
                gpcircle.x_dir,
                gpcircle.y_dir,
                apoint,
            );
            let (p, t) = elclib2d::circle_d1(
                gpcircle.center,
                gpcircle.x_dir,
                gpcircle.y_dir,
                gpcircle.radius,
                *aparameter,
            );
            point = p;
            tangent = t;
        } else if typ == GccIntIType::Hpr {
            let gphyperbola = abisector.hyperbola();
            let ydir = DVec2::new(-gphyperbola.major_dir.y, gphyperbola.major_dir.x);
            *aparameter = elclib2d::hyperbola_parameter(
                gphyperbola.center,
                ydir,
                gphyperbola.semi_minor,
                apoint,
            );
            let (p, t) = elclib2d::hyperbola_d1(
                gphyperbola.center,
                gphyperbola.major_dir,
                ydir,
                gphyperbola.semi_major,
                gphyperbola.semi_minor,
                *aparameter,
            );
            point = p;
            tangent = t;
        } else if typ == GccIntIType::Par {
            let gpparabola = abisector.parabola();
            let ydir = DVec2::new(-gpparabola.axis_dir.y, gpparabola.axis_dir.x);
            *aparameter = elclib2d::parabola_parameter(gpparabola.origin, ydir, apoint);
            let (p, t) = elclib2d::parabola_d1(
                gpparabola.origin,
                gpparabola.axis_dir,
                ydir,
                gpparabola.focal_param,
                *aparameter,
            );
            point = p;
            tangent = t;
        } else if typ == GccIntIType::Ell {
            let gpellipse = abisector.ellipse();
            *aparameter = elclib2d::ellipse_parameter(
                gpellipse.center,
                gpellipse.major_dir,
                gpellipse.minor_dir,
                gpellipse.major_radius,
                gpellipse.minor_radius,
                apoint,
            );
            let (p, t) = elclib2d::ellipse_d1(
                gpellipse.center,
                gpellipse.major_dir,
                gpellipse.minor_dir,
                gpellipse.major_radius,
                gpellipse.minor_radius,
                *aparameter,
            );
            point = p;
            tangent = t;
        } else {
            point = DVec2::ZERO;
            tangent = DVec2::ZERO;
        }

        let distance = apoint.distance(point);

        // OCCT gp_Dir2d constructions (normalized).
        let afirstdir = afirstvector.normalize_or_zero();
        let aseconddir = asecondvector.normalize_or_zero();
        let tangdir = tangent.normalize_or_zero();
        let secdirrev = -aseconddir;

        // 1st passage to learn if the curve is in the proper sector.
        if *asense {
            // the status is determined only in case on curve ie:
            // tangent to the bissectrice is bisectrice of two vectors.
            let sin_plat = 1.0e-3;
            if cross2(afirstdir, aseconddir).abs() < sin_plat {
                // flat
                if afirstdir.dot(aseconddir) >= 0.0 {
                    // tangent mixed
                    // correct if the scalar product is close to 1.
                    if tangdir.dot(afirstdir).abs() < 0.5 {
                        *astatus = false;
                    }
                } else {
                    // opposed tangents.
                    // correct if the scalar product is close to 0.
                    if tangdir.dot(afirstdir).abs() > 0.5 {
                        *astatus = false;
                    }
                }
            } else if cross2(afirstdir, tangdir) * cross2(tangdir, aseconddir) < -1.0e-8 {
                *astatus = false;
            }
        } else {
            *asense = true;

            if !is_bisec_of_two_lines {
                // Modified by Sergey KHROMOV — replacement of -1.E-8 for a
                // tolerance 1.e-4.
                let a_tol = 1.0e-4;

                if cross2(afirstdir, secdirrev) * adirection < -0.1 {
                    // input
                    if cross2(afirstdir, tangdir) * adirection < a_tol
                        && cross2(secdirrev, tangdir) * adirection < a_tol
                    {
                        *asense = false;
                    }
                } else if cross2(afirstdir, secdirrev) * adirection > 0.1 {
                    // output
                    if cross2(afirstdir, tangdir) * adirection < a_tol
                        || cross2(secdirrev, tangdir) * adirection < a_tol
                    {
                        *asense = false;
                    }
                } else {
                    // flat
                    if afirstdir.dot(secdirrev) > 0.0 {
                        // tangent
                        if cross2(afirstdir, tangdir) * adirection < 0.0 {
                            *asense = false;
                        }
                    } else {
                        // turn back
                        if afirstdir.dot(tangdir) < 0.0 {
                            *asense = false;
                        }
                    }
                }
                // jgv: for OCC26185
                if vec_ref.length_squared() != 0.0 {
                    let dir_ref = vec_ref.normalize_or_zero();
                    if tangdir.dot(dir_ref) < 0.0 {
                        *asense = false;
                    }
                }
            }
        }
        distance
    }

    /// OCCT Perform(Cu1, Cu2, P, V1, V2, Sense, jointype, Tolerance, oncurve)
    /// (L233-963) — the bissectrice between two curves coming from a point.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        afirstcurve: &Arc<dyn BisectorCurve>,
        asecondcurve: &Arc<dyn BisectorCurve>,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        ajointype: GeomAbsJoinType,
        tolerance: f64,
        oncurve: bool,
    ) {
        let mut ok = false;
        let mut distanceptsol;
        let mut parameter: f64 = 0.0;
        let mut firstparameter: f64 = 0.0;
        let mut thesense = false;
        let mut sense: bool;
        let mut distancemini;
        let pre_conf = CONFUSION;

        let mut thesol: Option<GccIntBisec> = None;

        // jgv: for OCC26296
        let mut line_bis_vec;
        let (pnt1, mut tan1) = afirstcurve.d1(afirstcurve.last_parameter());
        let (pnt2, mut tan2) = asecondcurve.d1(asecondcurve.first_parameter());
        if !oncurve {
            line_bis_vec = pnt2 - pnt1;
            line_bis_vec = rotate_half_pi(line_bis_vec);
        } else {
            line_bis_vec = DVec2::ZERO;
        }

        tan1 = -tan1;

        // OCCT: a trimmed input contributes its BasisCurve type.
        let curve_f = unwrap_trimmed(afirstcurve);
        let curve_e = unwrap_trimmed(asecondcurve);

        let mut circle1 = Circle2d::new(DVec2::ZERO, 0.0);
        let mut circle2 = Circle2d::new(DVec2::ZERO, 0.0);
        let mut line1 = Line2d::new(DVec2::ZERO, DVec2::X);
        let mut line2 = Line2d::new(DVec2::ZERO, DVec2::X);

        //==================================================================
        // Determination of the nature of arguments (L296-349).
        //==================================================================
        let mut cas = 0;
        match curve_f.kind() {
            CurveKind::Circle => match curve_e.kind() {
                CurveKind::Circle => {
                    cas = 1;
                    circle1 = circle_of(&curve_f);
                    circle2 = circle_of(&curve_e);
                }
                CurveKind::Line => {
                    cas = 2;
                    circle1 = circle_of(&curve_f);
                    line2 = line_of(&curve_e);
                }
                _ => {
                    println!("Not yet implemented");
                }
            },
            CurveKind::Line => match curve_e.kind() {
                CurveKind::Circle => {
                    // OCCT loads the circle into circle1 and the line into
                    // line2 (swapped operands).
                    cas = 2;
                    circle1 = circle_of(&curve_e);
                    line2 = line_of(&curve_f);
                }
                CurveKind::Line => {
                    cas = 3;
                    line1 = line_of(&curve_f);
                    line2 = line_of(&curve_e);
                }
                _ => {
                    println!("Not yet implemented");
                }
            },
            _ => {
                println!("Not yet implemented");
            }
        }

        match cas {
            //================================================================
            // Bissectrice circle - circle (OCCT case 1, L358-679).
            //================================================================
            1 => {
                let mut radius1 = circle1.radius;
                let mut radius2 = circle2.radius;

                //--------------------------------------------------
                // Particular case when two circles are mixed.
                //--------------------------------------------------
                if pnt2d_equal(circle1.center, circle2.center, pre_conf)
                    && (radius1 - radius2).abs() <= pre_conf
                {
                    let p1 = afirstcurve.value(afirstcurve.last_parameter());
                    let p2 = asecondcurve.value(asecondcurve.first_parameter());
                    let p_mil = DVec2::new((p1.x + p2.x) / 2., (p1.y + p2.y) / 2.);
                    let line;
                    if !pnt2d_equal(circle1.center, p_mil, pre_conf) {
                        // PMil doesn't coincide with the circle location.
                        line = Line2d::new(
                            p_mil,
                            DVec2::new(circle1.center.x - p_mil.x, circle1.center.y - p_mil.y),
                        );
                    } else if radius1 >= pre_conf {
                        // PMil coincides with the circle location and radius
                        // is greater then 0.
                        line = Line2d::new(
                            circle1.center,
                            DVec2::new(p1.y - circle1.center.y, circle1.center.x - p1.x),
                        );
                    } else {
                        // radius is equal to 0. No matter what direction to
                        // chose.
                        line = Line2d::new(circle1.center, DVec2::X);
                    }
                    let solution = GccIntBisec::Line(line);
                    sense = false;
                    if oncurve {
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            tan2,
                            tan1,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut ok,
                            false,
                        );
                    } else {
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            afirstvector,
                            asecondvector,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut ok,
                            false,
                        );
                    }
                    let _ = distanceptsol;
                    let bisectorcurve: Arc<dyn BisectorCurve> =
                        Arc::new(Geom2dCurveHandle::line(&line));
                    if !sense {
                        *self.thebisector.write().unwrap() =
                            Some(TrimmedCurve::new(bisectorcurve, parameter, -INFINITE_VALUE));
                    } else {
                        let mut parameter2 =
                            elclib2d::line_parameter(line.origin, line.direction, circle1.center);
                        parameter2 += 1.0e-8;
                        *self.thebisector.write().unwrap() =
                            Some(TrimmedCurve::new(bisectorcurve, parameter, parameter2));
                    }
                    // OCCT: break (end of the mixed-circles case).
                    return;
                }

                if radius1 < radius2 {
                    let circle = circle1;
                    circle1 = circle2;
                    circle2 = circle;

                    let radius = radius1;
                    radius1 = radius2;
                    radius2 = radius;
                }

                // small reframing of circles. in the case when the circles
                // are OnCurve , if they are almost tangent they become
                // tangent.
                let entre_axe = circle1.center.distance(circle2.center);
                let mut d1 = 0.5 * (radius1 - entre_axe - radius2);
                let mut circles_tangent = false;

                if oncurve && d1.abs() < pre_conf && dir2d_parallel(tan1, tan2, 1.0e-8) {
                    // C2 included in C1 and tangent.
                    circle1.radius = radius1 - d1;
                    circle2.radius = radius2 + d1;
                    circles_tangent = true;
                } else {
                    d1 = 0.5 * (radius1 - entre_axe + radius2);
                    if oncurve && d1.abs() < pre_conf && dir2d_parallel(tan1, tan2, 1.0e-8) {
                        // C2 and C1 tangent and disconnected.
                        circle1.radius = radius1 - d1;
                        circle2.radius = radius2 - d1;
                        circles_tangent = true;
                    }
                } // end of reframing.

                let bisector = GccAnaCirc2dBisec::new(circle1, circle2, tolerance);

                distancemini = INFINITE_VALUE;

                if bisector.is_done() {
                    let nbsolution = bisector.nb_solutions();
                    for i in 1..=nbsolution {
                        let mut solution = bisector.this_solution(i);
                        degenerate(&mut solution, tolerance);
                        sense = true;
                        if oncurve {
                            distanceptsol = self.distance(
                                apoint,
                                &solution,
                                tan1,
                                tan2,
                                line_bis_vec,
                                adirection,
                                &mut parameter,
                                &mut sense,
                                &mut ok,
                                false,
                            );
                        } else {
                            ok = true;
                        }
                        let _ = distanceptsol;

                        if ok {
                            sense = false;
                            if oncurve {
                                distanceptsol = self.distance(
                                    apoint,
                                    &solution,
                                    tan2,
                                    tan1,
                                    line_bis_vec,
                                    adirection,
                                    &mut parameter,
                                    &mut sense,
                                    &mut ok,
                                    false,
                                );
                            } else {
                                distanceptsol = self.distance(
                                    apoint,
                                    &solution,
                                    afirstvector,
                                    asecondvector,
                                    line_bis_vec,
                                    adirection,
                                    &mut parameter,
                                    &mut sense,
                                    &mut ok,
                                    false,
                                );
                            }
                            if distanceptsol <= distancemini {
                                thesol = Some(solution);
                                firstparameter = parameter;
                                thesense = sense;
                                distancemini = distanceptsol;
                            }
                        }
                    }
                    if let Some(ref mut the_sol) = thesol {
                        let typ = the_sol.arc_type();
                        if typ == GccIntIType::Lin {
                            let gpline = the_sol.line();
                            let mut bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::line(&gpline));

                            let mut secondparameter = INFINITE_VALUE;
                            if !thesense {
                                secondparameter = -INFINITE_VALUE;
                            }

                            if oncurve {
                                // bisectrice right and oncurve
                                // is cut between two circle of the same
                                // radius if circles are tangent.

                                // if tangent flat and the bissectrice at the
                                // side of the concavity of one of the
                                // circles. the bissectrice is a segment of
                                // the point common to first of 2 centers of
                                // circle that it meets. in this case it is
                                // important to set a segmnent for
                                // intersection in Tool2d.
                                if circles_tangent {
                                    // Trying to correct the line if the
                                    // distance between it and the reference
                                    // point is too big.
                                    if distancemini > tolerance {
                                        let a_ploc = gpline.origin;
                                        let a_new_dir = apoint - a_ploc;
                                        let a_new_lin = Line2d::new(a_ploc, a_new_dir);
                                        let a_cc2 = circle2.center;
                                        let a_new_d_min = a_new_lin.distance(apoint);
                                        let a_tol_conf = 1.0e-3;
                                        // Hope, aNewDMin is equal to 0...

                                        if a_new_lin.distance(a_cc2) <= a_tol_conf {
                                            distancemini = a_new_d_min;
                                            firstparameter = elclib2d::line_parameter(
                                                a_new_lin.origin,
                                                a_new_lin.direction,
                                                apoint,
                                            );
                                            bisectorcurve =
                                                Arc::new(Geom2dCurveHandle::line(&a_new_lin));
                                        }
                                    }
                                    if tan1.dot(tan2) < 0.0 {
                                        // flat and not turn back.
                                        let par1 = elclib2d::line_parameter(
                                            gpline.origin,
                                            gpline.direction,
                                            circle1.center,
                                        );
                                        let par2 = elclib2d::line_parameter(
                                            gpline.origin,
                                            gpline.direction,
                                            circle2.center,
                                        );
                                        let min_par = par1.min(par2);
                                        let max_par = par1.max(par2);

                                        if !thesense {
                                            if max_par < firstparameter {
                                                secondparameter = max_par - 1.0e-8;
                                            } else if min_par < firstparameter {
                                                secondparameter = min_par - 1.0e-8;
                                            }
                                        } else {
                                            if min_par > firstparameter {
                                                secondparameter = min_par + 1.0e-8;
                                            } else if max_par > firstparameter {
                                                secondparameter = max_par + 1.0e-8;
                                            }
                                        }
                                    }
                                }
                            }

                            *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                bisectorcurve,
                                firstparameter,
                                secondparameter,
                            ));
                        } else if typ == GccIntIType::Cir {
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::circle(&the_sol.circle()));
                            if !thesense {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter - 2.0 * PI,
                                    firstparameter,
                                    thesense,
                                ));
                            } else {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter,
                                    firstparameter + 2.0 * PI,
                                    thesense,
                                ));
                            }
                        } else if typ == GccIntIType::Hpr {
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::hyperbola(&the_sol.hyperbola()));
                            if !thesense {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    -INFINITE_VALUE,
                                ));
                            } else {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    INFINITE_VALUE,
                                ));
                            }
                        } else if typ == GccIntIType::Ell {
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::ellipse(&the_sol.ellipse()));
                            if !thesense {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter - 2.0 * PI,
                                    firstparameter,
                                    thesense,
                                ));
                            } else {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter,
                                    firstparameter + 2.0 * PI,
                                    thesense,
                                ));
                            }
                        }
                    }
                }
            }

            //================================================================
            // Bissectrice circle - straight (OCCT case 2, L685-834).
            //================================================================
            2 => {
                // small reframing of circles. in case OnCurve.
                // If the circle and the straight line are almost tangent
                // they become tangent.
                if oncurve {
                    let radius1 = circle1.radius;
                    let d1 = line2.distance(circle1.center) - radius1;
                    if d1.abs() < pre_conf && dir2d_parallel(tan1, tan2, 1.0e-8) {
                        circle1.radius = radius1 + d1;
                    }
                }

                let bisector = GccAnaCircLin2dBisec::new(circle1, line2);

                distancemini = INFINITE_VALUE;

                if bisector.is_done() {
                    let nbsolution = bisector.nb_solutions();
                    for i in 1..=nbsolution {
                        let mut solution = bisector.this_solution(i);
                        degenerate(&mut solution, tolerance);
                        sense = true;
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            tan1,
                            tan2,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut ok,
                            false,
                        );
                        let _ = distanceptsol;
                        if ok || !oncurve {
                            sense = false;
                            if oncurve {
                                distanceptsol = self.distance(
                                    apoint,
                                    &solution,
                                    tan2,
                                    tan1,
                                    line_bis_vec,
                                    adirection,
                                    &mut parameter,
                                    &mut sense,
                                    &mut ok,
                                    false,
                                );
                            } else {
                                distanceptsol = self.distance(
                                    apoint,
                                    &solution,
                                    afirstvector,
                                    asecondvector,
                                    line_bis_vec,
                                    adirection,
                                    &mut parameter,
                                    &mut sense,
                                    &mut ok,
                                    false,
                                );
                            }
                            if distanceptsol <= distancemini {
                                thesol = Some(solution);
                                firstparameter = parameter;
                                thesense = sense;
                                distancemini = distanceptsol + 1.0e-8;
                            }
                        }
                    }
                    if let Some(ref mut the_sol) = thesol {
                        let typ = the_sol.arc_type();
                        if typ == GccIntIType::Lin {
                            // If the bisectrice is a line
                            //       => the straight line is tangent to the
                            //          circle.
                            //       It the part of bisectrice concerned is at
                            //       the side of the center.
                            //       => the bisectrice is limited by the point
                            //          and the center of the circle.
                            // Note : In the latter case the bisectrice is a
                            // degenerated parabole.
                            let circlecenter = circle1.center;
                            let gpline = the_sol.line();
                            let mut secondparameter = elclib2d::line_parameter(
                                gpline.origin,
                                gpline.direction,
                                circlecenter,
                            );
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::line(&gpline));

                            if !thesense {
                                if secondparameter > firstparameter {
                                    secondparameter = -INFINITE_VALUE;
                                } else {
                                    secondparameter = secondparameter - 1.0e-8;
                                }
                            } else {
                                if secondparameter < firstparameter {
                                    secondparameter = INFINITE_VALUE;
                                } else {
                                    secondparameter = secondparameter + 1.0e-8;
                                }
                            }

                            *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                bisectorcurve,
                                firstparameter,
                                secondparameter,
                            ));
                        } else if typ == GccIntIType::Par {
                            let the_parabola = the_sol.parabola();
                            let parabola_curve = Curve2d::Parabola(the_parabola);
                            let apex = geom_value(&parabola_curve, 0.0);
                            let firstpnt = geom_value(&parabola_curve, firstparameter);
                            let chord_len = apex.distance(firstpnt);
                            let tol_par = 1.0e-5;
                            let mut secondparameter = INFINITE_VALUE;
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::parabola(&the_parabola));
                            if !thesense {
                                if ajointype == GeomAbsJoinType::Intersection
                                    && tol_par < firstparameter
                                    && chord_len >= circle1.radius
                                {
                                    // first parameter is too far from peak of
                                    // parabola
                                    secondparameter = 0.0;
                                }
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    -secondparameter,
                                ));
                            } else {
                                if ajointype == GeomAbsJoinType::Intersection
                                    && firstparameter < -tol_par
                                    && chord_len >= circle1.radius
                                {
                                    // first parameter is too far from peak of
                                    // parabola
                                    secondparameter = 0.0;
                                }
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    secondparameter,
                                ));
                            }
                        }
                    }
                }
            }

            //================================================================
            // Bissectrice straight - straight (OCCT case 3, L839-957).
            //================================================================
            3 => {
                let direc1 = line1.direction;
                let direc2 = line2.direction;
                let line;
                distancemini = INFINITE_VALUE;

                // Change to the same criterion as in MAT2d_Circuit.cxx:
                // method MAT2d_Circuit::InitOpen(..)
                if dir2d_parallel(direc1, direc2, 1.0e-8) {
                    if line1.distance(line2.origin) / 2. <= CONFUSION {
                        line = Line2d::new(apoint, rotate_half_pi(line1.direction));
                    } else {
                        line = Line2d::new(apoint, line2.direction);
                    }

                    let solution = GccIntBisec::Line(line);
                    sense = false;
                    if oncurve {
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            tan2,
                            tan1,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut ok,
                            false,
                        );
                    } else {
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            afirstvector,
                            asecondvector,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut ok,
                            false,
                        );
                    }
                    let _ = distanceptsol;
                    firstparameter = parameter;
                    let bisectorcurve: Arc<dyn BisectorCurve> =
                        Arc::new(Geom2dCurveHandle::line(&line));
                    if !sense {
                        *self.thebisector.write().unwrap() =
                            Some(TrimmedCurve::new(bisectorcurve, firstparameter, -INFINITE_VALUE));
                    } else {
                        *self.thebisector.write().unwrap() =
                            Some(TrimmedCurve::new(bisectorcurve, firstparameter, INFINITE_VALUE));
                    }
                } else {
                    let l = Line2d::new(apoint, direc2 - direc1);
                    let solution = GccIntBisec::Line(l);
                    let mut is_ok = false;
                    sense = false;
                    if oncurve {
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            tan2,
                            tan1,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut is_ok,
                            false,
                        );
                    } else {
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            afirstvector,
                            asecondvector,
                            line_bis_vec,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut is_ok,
                            true,
                        );
                    }
                    let _ = distanceptsol;
                    if is_ok || !oncurve {
                        thesense = sense;
                        distancemini = distanceptsol;
                    }
                    let _ = distancemini;
                    let the_sol_line = l;
                    let bisectorcurve: Arc<dyn BisectorCurve> =
                        Arc::new(Geom2dCurveHandle::line(&the_sol_line));
                    if !thesense {
                        *self.thebisector.write().unwrap() =
                            Some(TrimmedCurve::new(bisectorcurve, 0.0, -INFINITE_VALUE));
                    } else {
                        *self.thebisector.write().unwrap() =
                            Some(TrimmedCurve::new(bisectorcurve, 0.0, INFINITE_VALUE));
                    }
                }
            }

            _ => {
                // OCCT default: throw StdFail_NotDone().
                panic!("StdFail_NotDone: Bisector_BisecAna::Perform");
            }
        }
    }

    /// OCCT Perform(Cu, Pnt, P, V1, V2, Sense, Tolerance, oncurve) (L977-1223)
    /// — the bissectrice between a curve and a point.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_point(
        &mut self,
        afirstcurve: &Arc<dyn BisectorCurve>,
        asecondpoint: DVec2,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) {
        let mut ok = false;
        let mut distanceptsol;
        let mut parameter: f64 = 0.0;
        let mut firstparameter: f64 = 0.0;
        let mut secondparameter;
        let mut thesense = false;
        let mut sense: bool;
        let vec_ref = DVec2::ZERO;
        let mut thesol: Option<GccIntBisec> = None;

        let mut circle = Circle2d::new(DVec2::ZERO, 0.0);
        let mut line = Line2d::new(DVec2::ZERO, DVec2::X);

        let mut cas = 0;

        let curve = unwrap_trimmed(afirstcurve);

        match curve.kind() {
            CurveKind::Circle => {
                cas = 1;
                circle = circle_of(&curve);
            }
            CurveKind::Line => {
                cas = 2;
                line = line_of(&curve);
            }
            _ => {
                println!("Not yet implemented");
            }
        }

        match cas {
            //============================================================
            // Bissectrice point - circle (OCCT case 1, L1036-1163).
            //============================================================
            1 => {
                let bisector = GccAnaCircPnt2dBisec::new(circle, asecondpoint, tolerance);
                let mut distancemini = INFINITE_VALUE;
                if bisector.is_done() {
                    let nbsolution = bisector.nb_solutions();
                    for i in 1..=nbsolution {
                        let mut solution = bisector.this_solution(i);
                        degenerate(&mut solution, tolerance);
                        sense = false;
                        distanceptsol = self.distance(
                            apoint,
                            &solution,
                            afirstvector,
                            asecondvector,
                            vec_ref,
                            adirection,
                            &mut parameter,
                            &mut sense,
                            &mut ok,
                            false,
                        );

                        if distanceptsol <= distancemini {
                            thesol = Some(solution);
                            firstparameter = parameter;
                            thesense = sense;
                            distancemini = distanceptsol;
                        }
                    }
                    if let Some(ref mut the_sol) = thesol {
                        let sol_type = the_sol.arc_type();
                        if sol_type == GccIntIType::Lin {
                            // If the bisectrice is a line
                            //       => the point is on the circle.
                            //       If the part of bisectrice concerned is
                            //       at the side of the center.
                            //       => the bisectrice is limited by the
                            //          point and the center of the circle.
                            // Note : In this latter case the bisectrice is
                            // actually an ellipse of small null axis.
                            let circlecenter = circle.center;
                            line = the_sol.line();
                            secondparameter = elclib2d::line_parameter(
                                line.origin,
                                line.direction,
                                circlecenter,
                            );
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::line(&line));

                            if !thesense {
                                if secondparameter > firstparameter {
                                    secondparameter = -INFINITE_VALUE;
                                } else {
                                    secondparameter = secondparameter - 1.0e-8;
                                }
                            } else {
                                if secondparameter < firstparameter {
                                    secondparameter = INFINITE_VALUE;
                                } else {
                                    secondparameter = secondparameter + 1.0e-8;
                                }
                            }

                            *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                bisectorcurve,
                                firstparameter,
                                secondparameter,
                            ));
                        } else if sol_type == GccIntIType::Cir {
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::circle(&the_sol.circle()));
                            if !thesense {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter - 2.0 * PI,
                                    firstparameter,
                                    thesense,
                                ));
                            } else {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter,
                                    firstparameter + 2.0 * PI,
                                    thesense,
                                ));
                            }
                        } else if sol_type == GccIntIType::Hpr {
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::hyperbola(&the_sol.hyperbola()));
                            if !thesense {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    -INFINITE_VALUE,
                                ));
                            } else {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new(
                                    bisectorcurve,
                                    firstparameter,
                                    INFINITE_VALUE,
                                ));
                            }
                        } else if sol_type == GccIntIType::Ell {
                            let bisectorcurve: Arc<dyn BisectorCurve> =
                                Arc::new(Geom2dCurveHandle::ellipse(&the_sol.ellipse()));
                            if !thesense {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter - 2.0 * PI,
                                    firstparameter,
                                    thesense,
                                ));
                            } else {
                                *self.thebisector.write().unwrap() = Some(TrimmedCurve::new_with_sense(
                                    bisectorcurve,
                                    firstparameter,
                                    firstparameter + 2.0 * PI,
                                    thesense,
                                ));
                            }
                        }
                    }
                }
            }

            //============================================================
            // Bissectrice point - straight (OCCT case 2, L1168-1216).
            //============================================================
            2 => {
                let bisector = GccAnaLinPnt2dBisec::new(line, asecondpoint);

                let mut solution = bisector.this_solution();
                degenerate(&mut solution, tolerance);
                let typ = solution.arc_type();
                let bisectorcurve: Arc<dyn BisectorCurve> = if typ == GccIntIType::Lin {
                    Arc::new(Geom2dCurveHandle::line(&solution.line()))
                } else if typ == GccIntIType::Par {
                    Arc::new(Geom2dCurveHandle::parabola(&solution.parabola()))
                } else {
                    // OCCT: bisectorcurve stays null for other types.
                    *self.thebisector.write().unwrap() = None;
                    return;
                };
                sense = false;
                distanceptsol = self.distance(
                    apoint,
                    &solution,
                    afirstvector,
                    asecondvector,
                    vec_ref,
                    adirection,
                    &mut parameter,
                    &mut sense,
                    &mut ok,
                    false,
                );
                let _ = distanceptsol;

                if ok || !oncurve {
                    firstparameter = parameter;
                    thesense = sense;
                }

                if !thesense {
                    *self.thebisector.write().unwrap() =
                        Some(TrimmedCurve::new(bisectorcurve, firstparameter, -INFINITE_VALUE));
                } else {
                    *self.thebisector.write().unwrap() =
                        Some(TrimmedCurve::new(bisectorcurve, firstparameter, INFINITE_VALUE));
                }
            }

            _ => {
                println!("Not yet implemented");
            }
        }
    }

    /// OCCT Perform(Pnt, Cu, P, V1, V2, Sense, Tolerance, oncurve) (L1237-1259)
    /// — delegates with swapped arguments and reversed direction.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_curve(
        &mut self,
        afirstpoint: DVec2,
        asecondcurve: &Arc<dyn BisectorCurve>,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        _tolerance: f64,
        oncurve: bool,
    ) {
        let adirectionreverse = -adirection;
        self.perform_curve_point(
            asecondcurve,
            afirstpoint,
            apoint,
            asecondvector,
            afirstvector,
            adirectionreverse,
            0.0,
            oncurve,
        );
    }

    /// OCCT Perform(Pnt1, Pnt2, P, V1, V2, Sense, Tolerance, oncurve)
    /// (L1272-1304) — the bissectrice between two points.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_point(
        &mut self,
        afirstpoint: DVec2,
        asecondpoint: DVec2,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        _tolerance: f64,
        oncurve: bool,
    ) {
        let mut sense;
        let mut ok = false;
        let mut parameter: f64 = 0.0;
        let vec_ref = DVec2::ZERO;

        let bisector = GccAnaPnt2dBisec::new(afirstpoint, asecondpoint);
        let line = bisector.this_solution();
        let solution = GccIntBisec::Line(line);

        sense = false;
        self.distance(
            apoint,
            &solution,
            afirstvector,
            asecondvector,
            vec_ref,
            adirection,
            &mut parameter,
            &mut sense,
            &mut ok,
            false,
        );
        if ok || !oncurve {
            let bisectorcurve: Arc<dyn BisectorCurve> = Arc::new(Geom2dCurveHandle::line(&line));
            if !sense {
                *self.thebisector.write().unwrap() =
                    Some(TrimmedCurve::new(bisectorcurve, parameter, -INFINITE_VALUE));
            } else {
                *self.thebisector.write().unwrap() =
                    Some(TrimmedCurve::new(bisectorcurve, parameter, INFINITE_VALUE));
            }
        }
    }

    /// OCCT IsExtendAtStart() (L1308-1311).
    pub fn is_extend_at_start(&self) -> bool {
        false
    }

    /// OCCT IsExtendAtEnd() (L1315-1318).
    pub fn is_extend_at_end(&self) -> bool {
        false
    }

    /// OCCT SetTrim(const occ::handle<Geom2d_Curve>&) (L1329-1516) — the
    /// OCCT body is fully commented out (the function is void).
    pub fn set_trim_by_curve(&mut self, _cu: &Arc<dyn BisectorCurve>) {
        // OCCT body commented out — kept empty.
    }

    /// OCCT SetTrim(const double uf, const double ul) (L1518-1521).
    pub fn set_trim(&self, uf: f64, ul: f64) {
        if let Some(ref mut thebisector) = self.thebisector.write().unwrap().as_mut() {
            thebisector.set_trim(uf, ul);
        }
    }

    /// OCCT Reverse() (L1525-1528).
    pub fn reverse(&self) {
        if let Some(ref mut thebisector) = self.thebisector.write().unwrap().as_mut() {
            thebisector.reverse();
        }
    }

    /// OCCT ReversedParameter(U) (L1532-1535).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .reversed_parameter(u)
    }

    /// OCCT IsCN(N) (L1539-1542).
    pub fn is_cn(&self, n: i32) -> bool {
        self.thebisector.read().unwrap().as_ref().expect("thebisector").is_cn(n)
    }

    /// OCCT Copy() (L1546-1551).
    pub fn copy_bisec_ana(&self) -> BisectorBisecAna {
        let mut c = BisectorBisecAna::new();
        if let Some(ref thebisector) = self.thebisector.read().unwrap().as_ref() {
            // OCCT: C->Init(down_cast<Geom2d_TrimmedCurve>(thebisector->Copy())).
            let copied = thebisector.copy_curve();
            let trimmed = copied
                .as_any_arc()
                .downcast::<TrimmedCurve>()
                .expect("TrimmedCurve downcast");
            c.init_arc(trimmed);
        }
        c
    }

    /// OCCT Transform(T) (L1555-1558).
    pub fn transform(&self, t: &glam::DAffine2) {
        if let Some(ref mut thebisector) = self.thebisector.write().unwrap().as_mut() {
            thebisector.transform(t);
        }
    }

    /// OCCT FirstParameter() (L1562-1568).
    pub fn first_parameter(&self) -> f64 {
        // OCCT: return thebisector->FirstParameter().
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .first_parameter()
    }

    /// OCCT LastParameter() (L1572-1575).
    pub fn last_parameter(&self) -> f64 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .last_parameter()
    }

    /// OCCT IsClosed() (L1579-1582).
    pub fn is_closed(&self) -> bool {
        self.thebisector.read().unwrap().as_ref().expect("thebisector").is_closed()
    }

    /// OCCT IsPeriodic() (L1586-1589).
    pub fn is_periodic(&self) -> bool {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .is_periodic()
    }

    /// OCCT Continuity() (L1593-1596).
    pub fn continuity(&self) -> GeomAbsShape {
        self.thebisector.read().unwrap().as_ref().expect("thebisector").continuity()
    }

    /// OCCT EvalD0(U) (L1600-1603).
    pub fn eval_d0(&self, u: f64) -> DVec2 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .basis()
            .eval_d0(u)
    }

    /// OCCT EvalD1(U) (L1607-1610).
    pub fn eval_d1(&self, u: f64) -> ResD1 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .basis()
            .eval_d1(u)
    }

    /// OCCT EvalD2(U) (L1614-1617).
    pub fn eval_d2(&self, u: f64) -> ResD2 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .basis()
            .eval_d2(u)
    }

    /// OCCT EvalD3(U) (L1621-1624).
    pub fn eval_d3(&self, u: f64) -> ResD3 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .basis()
            .eval_d3(u)
    }

    /// OCCT EvalDN(U, N) (L1628-1631).
    pub fn eval_dn(&self, u: f64, n: i32) -> DVec2 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .basis()
            .eval_dn(u, n)
    }

    /// OCCT Geom2dCurve() (L1635-1638) — the basis curve handle.
    pub fn geom2d_curve(&self) -> Arc<dyn BisectorCurve> {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .basis()
            .clone()
    }

    /// OCCT ParameterOfStartPoint() (L1642-1645).
    pub fn parameter_of_start_point(&self) -> f64 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .first_parameter()
    }

    /// OCCT ParameterOfEndPoint() (L1649-1652).
    pub fn parameter_of_end_point(&self) -> f64 {
        self.thebisector.read().unwrap().as_ref().expect("thebisector")
            .last_parameter()
    }

    /// OCCT Parameter(P) (L1656-1693).
    pub fn parameter(&self, p: DVec2) -> f64 {
        let guard = self.thebisector.read().unwrap();
        let basis = guard.as_ref().expect("thebisector").basis().clone();
        drop(guard);

        match basis.kind() {
            CurveKind::Line => {
                let handle = downcast_handle(&basis);
                let l = line_of_kernel(&handle.curve);
                elclib2d::line_parameter(l.origin, l.direction, p)
            }
            CurveKind::Circle => {
                let handle = downcast_handle(&basis);
                let c = circle_of_kernel(&handle.curve);
                elclib2d::circle_parameter(c.center, c.x_dir, c.y_dir, p)
            }
            CurveKind::Hyperbola => {
                let handle = downcast_handle(&basis);
                let h = hyperbola_of(&handle.curve);
                elclib2d::hyperbola_parameter(
                    h.center,
                    DVec2::new(-h.major_dir.y, h.major_dir.x),
                    h.semi_minor,
                    p,
                )
            }
            CurveKind::Parabola => {
                let handle = downcast_handle(&basis);
                let pp = parabola_of(&handle.curve);
                elclib2d::parabola_parameter(
                    pp.origin,
                    DVec2::new(-pp.axis_dir.y, pp.axis_dir.x),
                    p,
                )
            }
            CurveKind::Ellipse => {
                let handle = downcast_handle(&basis);
                let e = ellipse_of(&handle.curve);
                elclib2d::ellipse_parameter(
                    e.center,
                    e.major_dir,
                    e.minor_dir,
                    e.major_radius,
                    e.minor_radius,
                    p,
                )
            }
            _ => 0.0,
        }
    }

    /// OCCT NbIntervals() (L1697-1700).
    pub fn nb_intervals(&self) -> i32 {
        1
    }

    /// OCCT IntervalFirst(I) (L1704-1711).
    pub fn interval_first(&self, i: i32) -> f64 {
        if i != 1 {
            panic!("Standard_OutOfRange: IntervalFirst");
        }
        self.first_parameter()
    }

    /// OCCT IntervalLast(I) (L1715-1722).
    pub fn interval_last(&self, i: i32) -> f64 {
        if i != 1 {
            panic!("Standard_OutOfRange: IntervalLast");
        }
        self.last_parameter()
    }

    /// OCCT Init(Bis) (L1727-1730).
    pub fn init(&mut self, bis: TrimmedCurve) {
        *self.thebisector.write().unwrap() = Some(bis);
    }

    /// Init over an Arc handle (used by Copy — the OCCT handle aliasing).
    pub(crate) fn init_arc(&mut self, bis: Arc<TrimmedCurve>) {
        *self.thebisector.write().unwrap() = Some((*bis).clone());
    }

    /// OCCT Dump(Deep, Offset) (L1819-1825).
    pub fn dump(&self, _deep: i32, _offset: i32) {
        println!("Bisector_BisecAna");
    }
}

// ---------------------------------------------------------------------------
// Bisector_Curve trait implementation for BisectorBisecAna
// (OCCT: class Bisector_BisecAna : public Bisector_Curve).
// ---------------------------------------------------------------------------

impl BisectorCurve for BisectorBisecAna {
    fn value(&self, u: f64) -> DVec2 {
        self.eval_d0(u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        let r = self.eval_d1(u);
        (r.point, r.d1)
    }

    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        let r = self.eval_d2(u);
        (r.point, r.d1, r.d2)
    }

    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        let r = self.eval_d3(u);
        (r.point, r.d1, r.d2, r.d3)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        self.eval_dn(u, n)
    }

    fn first_parameter(&self) -> f64 {
        BisectorBisecAna::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        BisectorBisecAna::last_parameter(self)
    }

    fn is_closed(&self) -> bool {
        BisectorBisecAna::is_closed(self)
    }

    fn is_periodic(&self) -> bool {
        BisectorBisecAna::is_periodic(self)
    }

    fn continuity(&self) -> GeomAbsShape {
        BisectorBisecAna::continuity(self)
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        BisectorBisecAna::reversed_parameter(self, u)
    }

    fn reverse(&mut self) {
        BisectorBisecAna::reverse(self)
    }

    fn is_cn(&self, n: i32) -> bool {
        BisectorBisecAna::is_cn(self, n)
    }

    fn transform(&mut self, t: &glam::DAffine2) {
        BisectorBisecAna::transform(self, t)
    }

    fn copy_curve(&self) -> Arc<dyn BisectorCurve> {
        Arc::new(self.copy_bisec_ana())
    }

    fn parameter(&self, p: DVec2) -> f64 {
        BisectorBisecAna::parameter(self, p)
    }

    fn is_extend_at_start(&self) -> bool {
        BisectorBisecAna::is_extend_at_start(self)
    }

    fn is_extend_at_end(&self) -> bool {
        BisectorBisecAna::is_extend_at_end(self)
    }

    fn nb_intervals(&self) -> i32 {
        BisectorBisecAna::nb_intervals(self)
    }

    fn interval_first(&self, index: i32) -> f64 {
        BisectorBisecAna::interval_first(self, index)
    }

    fn interval_last(&self, index: i32) -> f64 {
        BisectorBisecAna::interval_last(self, index)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::BisecAna
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_arc(self: Arc<Self>) -> Arc<dyn std::any::Any + Send + Sync> {
        self
    }
}
