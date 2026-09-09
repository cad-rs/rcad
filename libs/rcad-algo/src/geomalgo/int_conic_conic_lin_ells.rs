// OCCT IntCurve_IntConicConic::Perform(gp_Lin2d, gp_Elips2d) (TKGeomAlgo
// IntCurve) — the line/ellipse closed-form intersection and its two
// file-static helpers of IntCurve_IntConicConic_1.cxx.
//
// 1:1 translation of:
//   - LineEllipseGeometricIntersection          (_1.cxx L2657-2778)
//   - ProjectOnLAndIntersectWithLDomain         (_1.cxx L2780-2859)
//   - Perform(gp_Lin2d, gp_Elips2d)             (_1.cxx L2861-3215)
// plus the gp_Trsf2d::SetTransformation(XAxis) canonical-frame mapping
// (gp_Trsf2d.cxx L72-84, applied inline — rcad has no gp_Trsf2d) and the
// Extrema_ExtElC2d(gp_Lin2d, gp_Elips2d) constructor consumed by the
// degenerate branch (Extrema_ExtElC2d.cxx L170-219).

use glam::DVec2;
use rcad_kernel::geom::{Ellipse2d, Line2d};
use rcad_kernel::precision::INFINITE_VALUE;

use super::geom2d_int::elclib2d;
use super::int_conic_conic::IntConicConic;
use super::int_conic_conic_tool::{
    determine_transition_lc, normalize_on_circle_domain, Interval, PeriodicInterval, PI_PPI,
};
use super::int_imp_par_gen::determine_position;
use super::int_res2d::{
    Domain as Res2dDomain, IntersectionPoint, IntersectionSegment, Position, Transition,
};

/// OCCT Epsilon(Value) (Standard_Real.hxx L241-246) — the absolute value of
/// the difference between Value and its nearest representable neighbor.
fn epsilon_of(value: f64) -> f64 {
    if value >= 0.0 {
        value.next_up() - value
    } else {
        value - value.next_down()
    }
}

/// OCCT gp_Elips2d::YAxis().Direction() — the major direction rotated +90
/// degrees (rcad's Ellipse2d derives the minor direction from major_dir).
fn ellipse_ydir(ellipse: &Ellipse2d) -> DVec2 {
    DVec2::new(-ellipse.major_dir.y, ellipse.major_dir.x)
}

/// OCCT Extrema_ExtElC2d(const gp_Lin2d& C1, const gp_Elips2d& C2)
/// (Extrema_ExtElC2d.cxx L170-219) — the two stationary-distance solutions
/// of a line and an ellipse; consumed by LineEllipseGeometricIntersection
/// when the closed-form discriminant is negative. The consumer only reads
/// NbExt()/SquareDistance(i)/Points(i,.,.) through the ellipse extremum
/// parameter, so only those fields are kept (myDone ends true,
/// myIsPar false, as in OCCT).
struct ExtremaExtElC2dLineEllipse {
    /// OCCT mySqDist (RealLast-initialized, only the first myNbExt entries
    /// are valid).
    sq_dist: [f64; 2],
    /// OCCT myPoint[.][1].Parameter() — the ellipse parameter of each
    /// extremum.
    param2: [f64; 2],
    /// OCCT myNbExt.
    nb_ext: usize,
}

impl ExtremaExtElC2dLineEllipse {
    fn new(
        line_loc: DVec2,
        line_dir: DVec2,
        ellipse_loc: DVec2,
        ellipse_xdir: DVec2,
        ellipse_ydir: DVec2,
        major: f64,
        minor: f64,
    ) -> Self {
        let mut r = ExtremaExtElC2dLineEllipse {
            sq_dist: [f64::MAX; 2],
            param2: [0.0; 2],
            nb_ext: 0,
        };

        // OCCT: gp_Dir2d D = C1.Direction(); x2 = C2.XAxis().Direction();
        // y2 = C2.YAxis().Direction(); Dx = D.Dot(x2); Dy = D.Dot(y2);
        // r1 = C2.MajorRadius(); r2 = C2.MinorRadius(); O1 = C1.Location().
        let d_x = line_dir.dot(ellipse_xdir);
        let d_y = line_dir.dot(ellipse_ydir);

        let mut teta0 = if d_y.abs() <= f64::EPSILON {
            std::f64::consts::FRAC_PI_2
        } else {
            (-d_x * minor / (d_y * major)).atan()
        };

        let teta1 = teta0 + std::f64::consts::PI;
        if teta0 < 0.0 {
            teta0 += std::f64::consts::TAU;
        }

        // First extremum.
        let p2 = elclib2d::ellipse_value(ellipse_loc, ellipse_xdir, ellipse_ydir, major, minor, teta0);
        let u1 = (p2 - line_loc).dot(line_dir);
        let p1 = elclib2d::line_value(line_loc, line_dir, u1);
        r.sq_dist[r.nb_ext] = p1.distance_squared(p2);
        r.param2[r.nb_ext] = teta0;
        r.nb_ext += 1;

        // Second extremum.
        let p2 = elclib2d::ellipse_value(ellipse_loc, ellipse_xdir, ellipse_ydir, major, minor, teta1);
        let u1 = (p2 - line_loc).dot(line_dir);
        let p1 = elclib2d::line_value(line_loc, line_dir, u1);
        r.sq_dist[r.nb_ext] = p1.distance_squared(p2);
        r.param2[r.nb_ext] = teta1;
        r.nb_ext += 1;

        r
    }
}

/// OCCT LineEllipseGeometricIntersection (_1.cxx L2657-2778).
// The `= 0.` initializers of x1/y1/x2/y2 are part of the OCCT source
// (`double x1 = 0., y1 = 0., x2 = 0., y2 = 0.;`) and are kept 1:1, hence
// the unused_assignments allowance.
#[allow(clippy::too_many_arguments, unused_assignments)]
fn line_ellipse_geometric_intersection(
    line: &Line2d,
    ellipse: &Ellipse2d,
    _tol_conf: f64,
    tol_tang: f64,
    e_int1: &mut PeriodicInterval,
    e_int2: &mut PeriodicInterval,
    nbsol: &mut i32,
) {
    // OCCT: gp_Ax22d anElAxis = Ellipse.Axis(); gp_Trsf2d aTr;
    // aTr.SetTransformation(anElAxis.XAxis()); gp_Elips2d aTEllipse =
    // Ellipse.Transformed(aTr); gp_Lin2d aTLine = Line.Transformed(aTr).
    //
    // SetTransformation(A) (gp_Trsf2d.cxx L72-84) builds V1 = A.Direction(),
    // V2 = V1 rotated +90 degrees, then loc = -(M^T * A.Location()) with
    // M = [V1|V2]; a transformed point is M^T * (P - O), i.e. the canonical
    // coordinates in the direct frame of the ellipse major axis. The
    // transformed ellipse is the canonical one (its axis location maps to
    // exactly (0,0), its directions to (1,0)/(0,1)); the radii are unchanged.
    let v1 = ellipse.major_dir;
    let v2 = ellipse_ydir(ellipse);
    let a_t_ellipse_loc = DVec2::ZERO;
    // gp_Dir2d::Transform (gp_Dir2d.cxx L80-107) renormalizes directions.
    let a_t_ellipse_xdir = DVec2::new(v1.dot(v1), v2.dot(v1)).normalize_or_zero();
    let a_t_ellipse_ydir = DVec2::new(v1.dot(v2), v2.dot(v2)).normalize_or_zero();
    let a_t_line_loc = DVec2::new(
        v1.dot(line.origin - ellipse.center),
        v2.dot(line.origin - ellipse.center),
    );
    let a_t_line_dir_raw = DVec2::new(v1.dot(line.direction), v2.dot(line.direction));
    let a_t_line_dir = a_t_line_dir_raw / a_t_line_dir_raw.length();

    let a_dy = a_t_line_dir.y;
    let is_vert = a_dy.abs() > 1.0 - 2.0 * epsilon_of(1.0);
    //
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;
    let a2 = a * a;
    let b2 = b * b;
    //
    let mut eps0 = 1.0e-12;
    if b / a < 1.0e-5 {
        eps0 = 1.0e-6;
    }
    //
    // OCCT aTLine.Coefficients(anA, aB, aC) (gp_Lin2d.hxx L91-96).
    let an_a = a_t_line_dir.y;
    let mut a_b = -a_t_line_dir.x;
    let mut a_c = -(an_a * a_t_line_loc.x + a_b * a_t_line_loc.y);
    if is_vert {
        a_c += a_b * a_t_line_loc.y;
        a_b = 0.0;
    }
    //
    let mut x1 = 0.0;
    let mut y1 = 0.0;
    let mut x2 = 0.0;
    let mut y2 = 0.0;
    if a_b.abs() > eps0 {
        let m = -an_a / a_b;
        let m2 = m * m;
        let c = -a_c / a_b;
        let c2 = c * c;
        let mut d = a2 * m2 + b2 - c2;
        if d < 0.0 {
            // OCCT: Extrema_ExtElC2d anExt(aTLine, aTEllipse)
            // (Extrema_ExtElC2d.cxx L170-219); the closest stationary
            // distance is kept when within TolTang.
            let an_ext = ExtremaExtElC2dLineEllipse::new(
                a_t_line_loc,
                a_t_line_dir,
                a_t_ellipse_loc,
                a_t_ellipse_xdir,
                a_t_ellipse_ydir,
                a,
                b,
            );
            let mut imin = 0usize;
            let mut dmin = f64::MAX;
            for i in 1..=an_ext.nb_ext {
                if an_ext.sq_dist[i - 1] < dmin {
                    dmin = an_ext.sq_dist[i - 1];
                    imin = i;
                }
            }
            if imin > 0 && dmin <= tol_tang * tol_tang {
                *nbsol = 1;
                let pe1 = an_ext.param2[imin - 1];
                e_int1.set_values(pe1, pe1);
            } else {
                *nbsol = 0;
            }
            return;
        }
        d = d.sqrt();
        let n = a2 * m2 + b2;
        let k = a * b * d / n;
        let l = -a2 * m * c / n;
        x1 = l + k;
        y1 = m * x1 + c;
        x2 = l - k;
        y2 = m * x2 + c;
        *nbsol = 2;
    } else {
        x1 = -a_c / an_a;
        if x1.abs() > a + tol_tang {
            *nbsol = 0;
            return;
        } else if x1.abs() >= a - epsilon_of(1.0 + a) {
            *nbsol = 1;
            y1 = 0.0;
        } else {
            y1 = b * (1.0 - x1 * x1 / a2).sqrt();
            x2 = x1;
            y2 = -y1;
            *nbsol = 2;
        }
    }

    let a_p1 = DVec2::new(x1, y1);
    let a_p2 = DVec2::new(x2, y2);
    let mut pe1 = 0.0;
    let mut pe2 = 0.0;
    // OCCT ElCLib::Parameter(aTEllipse, aP1).
    pe1 = elclib2d::ellipse_parameter(
        a_t_ellipse_loc,
        a_t_ellipse_xdir,
        a_t_ellipse_ydir,
        a,
        b,
        a_p1,
    );
    if *nbsol > 1 {
        pe2 = elclib2d::ellipse_parameter(
            a_t_ellipse_loc,
            a_t_ellipse_xdir,
            a_t_ellipse_ydir,
            a,
            b,
            a_p2,
        );
        if pe2 < pe1 {
            let t = pe1;
            pe1 = pe2;
            pe2 = t;
        }
        e_int2.set_values(pe2, pe2);
    }
    e_int1.set_values(pe1, pe1);
}

/// OCCT ProjectOnLAndIntersectWithLDomain (_1.cxx L2780-2859).
#[allow(clippy::too_many_arguments)]
fn project_on_l_and_intersect_with_l_domain(
    ellipse: &Ellipse2d,
    line: &Line2d,
    e_domain_and_res: &PeriodicInterval,
    l_domain: &Interval,
    ellipse_solution: &mut [PeriodicInterval; 4],
    line_solution: &mut [Interval; 4],
    nb_sol_total: &mut usize,
    ref_line_domain: &Res2dDomain,
    _ref_ell_domain: &Res2dDomain,
) {
    if e_domain_and_res.is_null() {
        return;
    }
    //-------------------------------------------------------------------------
    //--  On cherche l intervalle correspondant sur C2
    //--  Puis on intersecte l intervalle avec le domaine de C2
    //--  Enfin, on cherche l intervalle correspondant sur C1
    //--

    let linf = elclib2d::line_parameter(
        line.origin,
        line.direction,
        elclib2d::ellipse_value(
            ellipse.center,
            ellipse.major_dir,
            ellipse_ydir(ellipse),
            ellipse.major_radius,
            ellipse.minor_radius,
            e_domain_and_res.binf,
        ),
    );
    let lsup = elclib2d::line_parameter(
        line.origin,
        line.direction,
        elclib2d::ellipse_value(
            ellipse.center,
            ellipse.major_dir,
            ellipse_ydir(ellipse),
            ellipse.major_radius,
            ellipse.minor_radius,
            e_domain_and_res.bsup,
        ),
    );

    let l_inter = Interval::new_bounded(linf, lsup); //-- Necessairement Borne

    let mut l_inter_and_domain = l_domain.intersection_with_bounded(&l_inter);

    if !l_inter_and_domain.is_null {
        let dom_linf = if ref_line_domain.has_first_point() {
            ref_line_domain.first_parameter()
        } else {
            -INFINITE_VALUE
        };
        let dom_lsup = if ref_line_domain.has_last_point() {
            ref_line_domain.last_parameter()
        } else {
            INFINITE_VALUE
        };

        let mut linf = l_inter_and_domain.binf;
        let mut lsup = l_inter_and_domain.bsup;

        if linf < dom_linf {
            linf = dom_linf;
        }
        if lsup < dom_linf {
            lsup = dom_linf;
        }

        if linf > dom_lsup {
            linf = dom_lsup;
        }
        if lsup > dom_lsup {
            lsup = dom_lsup;
        }

        l_inter_and_domain.binf = linf;
        l_inter_and_domain.bsup = lsup;

        let mut einf = e_domain_and_res.binf;
        let mut esup = e_domain_and_res.bsup;

        if einf >= esup {
            einf = e_domain_and_res.binf;
            esup = e_domain_and_res.bsup;
        }
        ellipse_solution[*nb_sol_total] = PeriodicInterval::new_ab(einf, esup);
        if ellipse_solution[*nb_sol_total].length() > std::f64::consts::PI {
            ellipse_solution[*nb_sol_total].complement();
        }

        line_solution[*nb_sol_total] = l_inter_and_domain;
        *nb_sol_total += 1;
    }
}

impl IntConicConic {
    /// OCCT Perform(const gp_Lin2d& L, const IntRes2d_Domain& DL,
    /// const gp_Elips2d& E, const IntRes2d_Domain& DE,
    /// const double TolConf, const double Tol)
    /// (IntCurve_IntConicConic_1.cxx L2861-3215).
    // The working variables (P1a/P2a/P1b/P2b, Pos1a..) mirror the OCCT
    // declarations and are always (re)assigned before read, hence the
    // unused_assignments allowance.
    #[allow(unused_assignments)]
    pub fn perform_line_ellipse(
        &mut self,
        l: &Line2d,
        dl: &Res2dDomain,
        e: &Ellipse2d,
        de: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let the_reversed_parameters = self.base.reversed_parameters();
        self.base.reset_fields();
        self.base.set_reversed_parameters(the_reversed_parameters);

        let mut nbsol = 0i32;
        let mut e_int1 = PeriodicInterval::new();
        let mut e_int2 = PeriodicInterval::new();

        line_ellipse_geometric_intersection(
            l,
            e,
            tol_conf,
            tol,
            &mut e_int1,
            &mut e_int2,
            &mut nbsol,
        );
        self.base.done = true;
        if nbsol == 0 {
            return;
        }
        //
        if nbsol == 2 && e_int2.bsup == e_int1.binf + PI_PPI {
            let first_bound = de.first_parameter();
            let last_bound = de.last_parameter();
            let first_tol = de.first_tolerance();
            let last_tol = de.last_tolerance();
            if e_int1.binf == 0.0 && first_bound - first_tol > e_int1.bsup {
                nbsol = 1;
                e_int1.set_values(e_int2.binf, e_int2.bsup);
            } else if e_int2.bsup == PI_PPI && last_bound + last_tol < e_int2.binf {
                nbsol = 1;
            }
        }
        //
        let mut e_domain = PeriodicInterval::from_domain(de);
        let mut deltat = e_domain.bsup - e_domain.binf;
        while e_domain.binf >= PI_PPI {
            e_domain.binf -= PI_PPI;
        }
        while e_domain.binf < 0.0 {
            e_domain.binf += PI_PPI;
        }
        e_domain.bsup = e_domain.binf + deltat;
        //
        let mut binf_modif = e_domain.binf;
        let mut bsup_modif = e_domain.bsup;
        binf_modif -= de.first_tolerance() / e.minor_radius;
        bsup_modif += de.last_tolerance() / e.minor_radius;
        deltat = bsup_modif - binf_modif;
        if deltat <= PI_PPI {
            e_domain.binf = binf_modif;
            e_domain.bsup = bsup_modif;
        } else {
            let mut t = PI_PPI - deltat;
            t *= 0.5;
            e_domain.binf = binf_modif + t;
            e_domain.bsup = bsup_modif - t;
        }
        deltat = e_domain.bsup - e_domain.binf;
        while e_domain.binf >= PI_PPI {
            e_domain.binf -= PI_PPI;
        }
        while e_domain.binf < 0.0 {
            e_domain.binf += PI_PPI;
        }
        e_domain.bsup = e_domain.binf + deltat;
        //
        let l_domain = Interval::from_domain(dl);

        let mut nb_sol_total = 0usize;

        let mut solution_ellipse = [PeriodicInterval::new(); 4];
        let mut solution_line = [Interval::new(); 4];

        //----------------------------------------------------------------------
        //----------- Treatment of first geometric interval EInt1           ----
        //----------------------------------------------------------------------
        // OCCT passes EInt1 by non-const reference: FirstIntersection slides
        // it in place, and the shifted value is what SecondIntersection reads.
        let mut e_domain_and_res = e_domain.first_intersection(&mut e_int1);

        project_on_l_and_intersect_with_l_domain(
            e,
            l,
            &e_domain_and_res,
            &l_domain,
            &mut solution_ellipse,
            &mut solution_line,
            &mut nb_sol_total,
            dl,
            de,
        );

        e_domain_and_res = e_domain.second_intersection(&e_int1);

        project_on_l_and_intersect_with_l_domain(
            e,
            l,
            &e_domain_and_res,
            &l_domain,
            &mut solution_ellipse,
            &mut solution_line,
            &mut nb_sol_total,
            dl,
            de,
        );

        //----------------------------------------------------------------------
        //----------- Treatment of second geometric interval EInt2          ----
        //----------------------------------------------------------------------
        if nbsol == 2 {
            e_domain_and_res = e_domain.first_intersection(&mut e_int2);

            project_on_l_and_intersect_with_l_domain(
                e,
                l,
                &e_domain_and_res,
                &l_domain,
                &mut solution_ellipse,
                &mut solution_line,
                &mut nb_sol_total,
                dl,
                de,
            );

            e_domain_and_res = e_domain.second_intersection(&e_int2);

            project_on_l_and_intersect_with_l_domain(
                e,
                l,
                &e_domain_and_res,
                &l_domain,
                &mut solution_ellipse,
                &mut solution_line,
                &mut nb_sol_total,
                dl,
                de,
            );
        }

        //----------------------------------------------------------------------
        //-- Calculation of Transitions at Positions.
        //----------------------------------------------------------------------
        let r = e.minor_radius;
        let mut max_tol = tol_conf;
        if max_tol < tol {
            max_tol = tol;
        }
        if max_tol < 1.0e-10 {
            max_tol = 1.0e-10;
        }

        for i in 0..nb_sol_total {
            if (r * solution_ellipse[i].length()) < max_tol && (solution_line[i].length()) < max_tol
            {
                let mut t = (solution_ellipse[i].binf + solution_ellipse[i].bsup) * 0.5;
                solution_ellipse[i].binf = t;
                solution_ellipse[i].bsup = t;

                t = (solution_line[i].binf + solution_line[i].bsup) * 0.5;
                solution_line[i].binf = t;
                solution_line[i].bsup = t;
            }
        }
        //
        if nb_sol_total > 0 {
            // OCCT: gp_Ax22d EllipseAxis = E.Axis(); gp_Ax2d LineAxis =
            // L.Position() — rcad keeps the frame as two directions.
            let ellipse_xdir = e.major_dir;
            let ellipse_ydir = ellipse_ydir(e);
            let mut p1a;
            let mut p2a;
            let mut p1b;
            let mut p2b;
            let mut tan1;
            let norm2 = DVec2::ZERO;
            let mut t1a = Transition::empty();
            let mut t2a = Transition::empty();
            let mut t1b = Transition::empty();
            let mut t2b = Transition::empty();
            let mut pos1a;
            let mut pos1b;
            let mut pos2a;
            let mut pos2b;

            let (pa, ta) = elclib2d::ellipse_d1(
                e.center,
                ellipse_xdir,
                ellipse_ydir,
                e.major_radius,
                e.minor_radius,
                solution_ellipse[0].binf,
            );
            p1a = pa;
            tan1 = ta;
            let (pb, tb) = elclib2d::line_d1(l.origin, l.direction, solution_line[0].binf);
            p2a = pb;
            let mut tan2 = tb;
            p1b = DVec2::ZERO;
            p2b = DVec2::ZERO;
            let mut norm1 = DVec2::ZERO;
            pos1a = Position::Middle;
            pos1b = Position::Middle;
            pos2a = Position::Middle;
            pos2b = Position::Middle;

            let is_opposite = tan1.dot(tan2) < 0.0;
            for i in 0..nb_sol_total {
                //-- 7 aout 97
                //-- On recentre Bin et Bsup de facon a avoir une portion commune avec DE
                let mut p1 = solution_ellipse[i].binf;
                let mut p2 = solution_ellipse[i].bsup;
                let q1 = de.first_parameter();
                let q2 = de.last_parameter();

                if p1 > q2 {
                    while p1 > q2 {
                        p1 -= PI_PPI;
                        p2 -= PI_PPI;
                    }
                } else if p2 < q1 {
                    while p2 < q1 {
                        p1 += PI_PPI;
                        p2 += PI_PPI;
                    }
                }
                if p1 < q1 && p2 > q1 {
                    p1 = q1;
                }
                if p1 < q2 && p2 > q2 {
                    p2 = q2;
                }

                solution_ellipse[i].binf = p1;
                solution_ellipse[i].bsup = p2;

                //-- Fin 7 aout 97

                let mut linf = if is_opposite {
                    solution_line[i].bsup
                } else {
                    solution_line[i].binf
                };
                let mut lsup = if is_opposite {
                    solution_line[i].binf
                } else {
                    solution_line[i].bsup
                };

                //---------------------------------------------------------------
                //-- Si les parametres sur l ellipse sont en premier
                //-- On doit retourner ces parametres dans l ordre croissant
                //---------------------------------------------------------------
                if linf > lsup {
                    let t = solution_ellipse[i].binf;
                    solution_ellipse[i].binf = solution_ellipse[i].bsup;
                    solution_ellipse[i].bsup = t;

                    let t = linf;
                    linf = lsup;
                    lsup = t;
                }

                let (pa, ta, na) = elclib2d::ellipse_d2(
                    e.center,
                    ellipse_xdir,
                    ellipse_ydir,
                    e.major_radius,
                    e.minor_radius,
                    solution_ellipse[i].binf,
                );
                p1a = pa;
                tan1 = ta;
                norm1 = na;
                let (pb, tb) = elclib2d::line_d1(l.origin, l.direction, linf);
                p2a = pb;
                tan2 = tb;

                determine_position(&mut pos1a, de, p1a, solution_ellipse[i].binf);
                determine_position(&mut pos2a, dl, p2a, linf);
                determine_transition_lc(
                    pos1a,
                    &mut tan1,
                    norm1,
                    &mut t1a,
                    pos2a,
                    &mut tan2,
                    norm2,
                    &mut t2a,
                    tol,
                );
                let einf;
                if pos1a == Position::End {
                    einf = de.last_parameter();
                    p1a = de.last_point();
                    linf = elclib2d::line_parameter(l.origin, l.direction, p1a);

                    let (pa, ta, na) = elclib2d::ellipse_d2(
                        e.center,
                        ellipse_xdir,
                        ellipse_ydir,
                        e.major_radius,
                        e.minor_radius,
                        einf,
                    );
                    p1a = pa;
                    tan1 = ta;
                    norm1 = na;
                    let (pb, tb) = elclib2d::line_d1(l.origin, l.direction, linf);
                    p2a = pb;
                    tan2 = tb;
                    determine_position(&mut pos1a, de, p1a, einf);
                    determine_position(&mut pos2a, dl, p2a, linf);
                    determine_transition_lc(
                        pos1a,
                        &mut tan1,
                        norm1,
                        &mut t1a,
                        pos2a,
                        &mut tan2,
                        norm2,
                        &mut t2a,
                        tol,
                    );
                } else if pos1a == Position::Head {
                    einf = de.first_parameter();
                    p1a = de.first_point();
                    linf = elclib2d::line_parameter(l.origin, l.direction, p1a);

                    let (pa, ta, na) = elclib2d::ellipse_d2(
                        e.center,
                        ellipse_xdir,
                        ellipse_ydir,
                        e.major_radius,
                        e.minor_radius,
                        einf,
                    );
                    p1a = pa;
                    tan1 = ta;
                    norm1 = na;
                    let (pb, tb) = elclib2d::line_d1(l.origin, l.direction, linf);
                    p2a = pb;
                    tan2 = tb;
                    determine_position(&mut pos1a, de, p1a, einf);
                    determine_position(&mut pos2a, dl, p2a, linf);
                    determine_transition_lc(
                        pos1a,
                        &mut tan1,
                        norm1,
                        &mut t1a,
                        pos2a,
                        &mut tan2,
                        norm2,
                        &mut t2a,
                        tol,
                    );
                } else {
                    einf = normalize_on_circle_domain(solution_ellipse[i].binf, de);
                }

                let mut new_point1 = IntersectionPoint::empty();
                new_point1.set_values(p1a, linf, einf, t2a, t1a, the_reversed_parameters);

                if (solution_line[i].length() + solution_ellipse[i].length()) > 0.0 {
                    let (pb1, tb1, nb1) = elclib2d::ellipse_d2(
                        e.center,
                        ellipse_xdir,
                        ellipse_ydir,
                        e.major_radius,
                        e.minor_radius,
                        solution_ellipse[i].bsup,
                    );
                    p1b = pb1;
                    tan1 = tb1;
                    norm1 = nb1;
                    let (pb2, tb2) = elclib2d::line_d1(l.origin, l.direction, lsup);
                    p2b = pb2;
                    tan2 = tb2;

                    determine_position(&mut pos1b, de, p1b, solution_ellipse[i].bsup);
                    determine_position(&mut pos2b, dl, p2b, lsup);
                    determine_transition_lc(
                        pos1b,
                        &mut tan1,
                        norm1,
                        &mut t1b,
                        pos2b,
                        &mut tan2,
                        norm2,
                        &mut t2b,
                        tol,
                    );
                    let esup;
                    if pos1b == Position::End {
                        // OCCT as written: `Esup = DL.LastParameter();` — the
                        // source reads the line domain here (the circ overload
                        // reads CIRC_Domain); kept 1:1.
                        esup = dl.last_parameter();
                        p1b = de.last_point();
                        lsup = elclib2d::line_parameter(l.origin, l.direction, p1b);

                        let (pa, ta, na) = elclib2d::ellipse_d2(
                            e.center,
                            ellipse_xdir,
                            ellipse_ydir,
                            e.major_radius,
                            e.minor_radius,
                            esup,
                        );
                        p1b = pa;
                        tan1 = ta;
                        norm1 = na;
                        let (pb, tb) = elclib2d::line_d1(l.origin, l.direction, lsup);
                        p2b = pb;
                        tan2 = tb;

                        determine_position(&mut pos1b, de, p1b, esup);
                        determine_position(&mut pos2b, dl, p2b, lsup);
                        determine_transition_lc(
                            pos1b,
                            &mut tan1,
                            norm1,
                            &mut t1b,
                            pos2b,
                            &mut tan2,
                            norm2,
                            &mut t2b,
                            tol,
                        );
                    } else if pos1b == Position::Head {
                        esup = de.first_parameter();
                        p1b = de.first_point();
                        lsup = elclib2d::line_parameter(l.origin, l.direction, p1b);

                        let (pa, ta, na) = elclib2d::ellipse_d2(
                            e.center,
                            ellipse_xdir,
                            ellipse_ydir,
                            e.major_radius,
                            e.minor_radius,
                            esup,
                        );
                        p1b = pa;
                        tan1 = ta;
                        norm1 = na;
                        let (pb, tb) = elclib2d::line_d1(l.origin, l.direction, lsup);
                        p2b = pb;
                        tan2 = tb;

                        determine_position(&mut pos1b, de, p1b, esup);
                        determine_position(&mut pos2b, dl, p2b, lsup);
                        determine_transition_lc(
                            pos1b,
                            &mut tan1,
                            norm1,
                            &mut t1b,
                            pos2b,
                            &mut tan2,
                            norm2,
                            &mut t2b,
                            tol,
                        );
                    } else {
                        esup = normalize_on_circle_domain(solution_ellipse[i].bsup, de);
                    }

                    let mut new_point2 = IntersectionPoint::empty();
                    new_point2.set_values(p1b, lsup, esup, t2b, t1b, the_reversed_parameters);

                    if (((esup - einf).abs() * r > max_tol) && ((lsup - linf).abs() > max_tol))
                        || (t1a.transition_type() != t2a.transition_type())
                    {
                        //-- Verifier egalement les transitions

                        let new_seg = IntersectionSegment::with_points(
                            &new_point1,
                            &new_point2,
                            is_opposite,
                            the_reversed_parameters,
                        );
                        self.base.append_segment(&new_seg);
                    } else {
                        // OCCT _1.cxx L3201/L3205: Insert(NewPoint1); Insert(NewPoint2);
                        if pos1a != Position::Middle || pos2a != Position::Middle {
                            self.base.insert(&new_point1);
                        }
                        if pos1b != Position::Middle || pos2b != Position::Middle {
                            self.base.insert(&new_point2);
                        }
                    }
                } else {
                    // OCCT _1.cxx L3211: Insert(NewPoint1);
                    self.base.insert(&new_point1);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::int_conic_conic::IntConicConic;
    use crate::geomalgo::int_res2d::Domain as Res2dDomain;

    /// The X-axis segment [-3, 3] (parameter u in [0, 6]).
    fn x_axis_line() -> (Line2d, Res2dDomain) {
        let line = Line2d {
            origin: DVec2::new(-3.0, 0.0),
            direction: DVec2::new(1.0, 0.0),
        };
        let mut d = Res2dDomain::infinite();
        d.set_values_bounded(
            line.origin,
            0.0,
            1.0e-9,
            DVec2::new(3.0, 0.0),
            6.0,
            1.0e-9,
        );
        (line, d)
    }

    /// The standard ellipse: center origin, major radius 2 along X, minor
    /// radius 1 (parameter u in [0, 2*PI], closed — the
    /// IntCurveCurveGen::ComputeDomain convention).
    fn std_ellipse() -> (Ellipse2d, Res2dDomain) {
        let ellipse = Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::new(1.0, 0.0),
            minor_dir: DVec2::new(0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 1.0,
        };
        let mut d = Res2dDomain::infinite();
        d.set_values_bounded(
            DVec2::new(2.0, 0.0),
            0.0,
            1.0e-9,
            DVec2::new(2.0, 0.0),
            std::f64::consts::TAU,
            1.0e-9,
        );
        d.set_equivalent_parameters(0.0, std::f64::consts::TAU);
        (ellipse, d)
    }

    /// OCCT anchor: the X-axis crossing the standard ellipse (a=2, b=1) at
    /// (-2, 0) and (2, 0) — the closed form gives x = +/-a, the ellipse
    /// parameters 0 and PI, and the Perform loop reports one point per
    /// geometric interval (Perform(Lin, Elps) _1.cxx L2861-3215 through
    /// LineEllipseGeometricIntersection).
    #[test]
    fn line_ellipse_two_crossings() {
        let (line, dl) = x_axis_line();
        let (ellipse, de) = std_ellipse();

        let mut ic = IntConicConic::new();
        ic.perform_line_ellipse(&line, &dl, &ellipse, &de, 1.0e-9, 1.0e-9);

        assert!(ic.base.done);
        assert_eq!(ic.base.nb_points(), 2, "nb_pt={}", ic.base.nb_points());
        assert_eq!(ic.base.nb_segments(), 0, "nb_seg={}", ic.base.nb_segments());

        // The two crossings (-2, 0) and (2, 0), sorted by x to be
        // order-independent.
        let mut pts = [
            (
                ic.base.point(1).value().x,
                ic.base.point(1).value().y,
                ic.base.point(1).param_on_first(),
                ic.base.point(1).param_on_second(),
            ),
            (
                ic.base.point(2).value().x,
                ic.base.point(2).value().y,
                ic.base.point(2).param_on_first(),
                ic.base.point(2).param_on_second(),
            ),
        ];
        pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // (-2, 0): line parameter 1 (origin -3, unit direction), ellipse
        // parameter PI.
        assert!((pts[0].0 + 2.0).abs() < 1.0e-6, "x={}", pts[0].0);
        assert!(pts[0].1.abs() < 1.0e-6, "y={}", pts[0].1);
        assert!((pts[0].2 - 1.0).abs() < 1.0e-6, "l={}", pts[0].2);
        assert!(
            (pts[0].3 - std::f64::consts::PI).abs() < 1.0e-6,
            "e={}",
            pts[0].3
        );
        // (2, 0): line parameter 5, ellipse parameter 0.
        assert!((pts[1].0 - 2.0).abs() < 1.0e-6, "x={}", pts[1].0);
        assert!(pts[1].1.abs() < 1.0e-6, "y={}", pts[1].1);
        assert!((pts[1].2 - 5.0).abs() < 1.0e-6, "l={}", pts[1].2);
        assert!(
            pts[1].3.abs() < 1.0e-6,
            "e={}",
            pts[1].3
        );
    }

    /// OCCT anchor: the tangent line y = 1 touches the standard ellipse at
    /// (0, 1) — the closed form discriminant vanishes (D == 0), so nbsol = 2
    /// with both solutions at the same parameter PI/2: Perform reports that
    /// point once per geometric interval (duplicated, the OCCT as-written
    /// behavior for the degenerate tangent), each with the line parameter 3
    /// and the TOUCH transitions of Determine_Transition_LC.
    #[test]
    fn line_ellipse_tangent_point() {
        let (ellipse, de) = std_ellipse();

        let line = Line2d {
            origin: DVec2::new(-3.0, 1.0),
            direction: DVec2::new(1.0, 0.0),
        };
        let mut dl = Res2dDomain::infinite();
        dl.set_values_bounded(
            line.origin,
            0.0,
            1.0e-9,
            DVec2::new(3.0, 1.0),
            6.0,
            1.0e-9,
        );

        let mut ic = IntConicConic::new();
        ic.perform_line_ellipse(&line, &dl, &ellipse, &de, 1.0e-9, 1.0e-9);

        assert!(ic.base.done);
        assert_eq!(ic.base.nb_segments(), 0, "nb_seg={}", ic.base.nb_segments());
        assert!(ic.base.nb_points() >= 1, "nb_pt={}", ic.base.nb_points());

        // Every reported point is the tangency (0, 1).
        for i in 1..=ic.base.nb_points() {
            let p = ic.base.point(i);
            let v = p.value();
            assert!(v.x.abs() < 1.0e-6, "pt{} x={}", i, v.x);
            assert!((v.y - 1.0).abs() < 1.0e-6, "pt{} y={}", i, v.y);
            // The line parameter of the tangency (origin x = -3 -> u = 3).
            assert!(
                (p.param_on_first() - 3.0).abs() < 1.0e-6,
                "pt{} l={}",
                i,
                p.param_on_first()
            );
            // The ellipse parameter of the tangency (PI/2), recadre into the
            // domain period.
            let c = p.param_on_second();
            let pi_2 = std::f64::consts::FRAC_PI_2;
            assert!(
                (c - pi_2).abs() < 1.0e-6
                    || (c - pi_2 - std::f64::consts::TAU).abs() < 1.0e-6
                    || (c - pi_2 + std::f64::consts::TAU).abs() < 1.0e-6,
                "pt{} c={}",
                i,
                c
            );
            // TOUCH: both transitions share the tangent type
            // (Determine_Transition_LC).
            assert_eq!(
                p.transition_of_first().transition_type(),
                p.transition_of_second().transition_type(),
                "pt{} transitions",
                i
            );
        }
    }
}
