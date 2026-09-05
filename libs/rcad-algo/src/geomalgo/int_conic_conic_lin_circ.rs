// OCCT IntCurve_IntConicConic::Perform(gp_Lin2d, gp_Circ2d) (TKGeomAlgo
// IntCurve) — the line/circle closed-form intersection and its two
// file-static helpers of IntCurve_IntConicConic_1.cxx.
//
// 1:1 translation of:
//   - LineCircleGeometricIntersection          (_1.cxx L451-640)
//   - ProjectOnLAndIntersectWithLDomain        (_1.cxx L360-443)
//   - Perform(gp_Lin2d, gp_Circ2d)             (_1.cxx L2236-2652)
// plus the gp helpers consumed by the body (gp_Lin2d::Coefficients
// gp_Lin2d.hxx L91-96, gp_Dir2d::Angle gp_Dir2d.cxx L26-59).

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Line2d};
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

/// OCCT gp_Dir2d::Angle(Other) (gp_Dir2d.cxx L26-59) — the angle from self
/// to other in [-PI, PI], with the acos/asin precision switching.
fn dir2d_angle(d1: DVec2, d2: DVec2) -> f64 {
    let cosinus = d1.dot(d2);
    let sinus = d1.x * d2.y - d1.y * d2.x;
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        if sinus > 0.0 {
            cosinus.acos()
        } else {
            -cosinus.acos()
        }
    } else if cosinus > 0.0 {
        sinus.asin()
    } else if sinus > 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        -std::f64::consts::PI - sinus.asin()
    }
}

/// OCCT LineCircleGeometricIntersection (_1.cxx L451-640).
fn line_circle_geometric_intersection(
    line: &Line2d,
    circle: &Circle2d,
    tol: f64,
    tol_tang: f64,
    c_int1: &mut PeriodicInterval,
    c_int2: &mut PeriodicInterval,
    nbsol: &mut i32,
) {
    let d_o1_o2 = line.distance(circle.center);
    let r = circle.radius;
    let r_m_tol = r - tol;
    let mut binf1;
    let mut binf2;
    let mut bsup1;
    let mut bsup2;

    //----------------------------------------------------------------
    if d_o1_o2 > (r + tol) {
        //-- pas d intersection avec le 'tuyau'
        if d_o1_o2 > (r + tol_tang) {
            *nbsol = 0;
            return;
        } else {
            binf1 = 0.0;
            bsup1 = 0.0;
            binf2 = 0.0;
            bsup2 = 0.0;
            *nbsol = 1;
        }
    } else {
        //----------------------------------------------------------------
        let mut b2_sol;
        let mut d_alpha1;
        //---------------------------------------------------------------
        //-- Line coupe le cercle Circle+ (=C(x1,y1,R1+Tol))
        b2_sol = false;
        if r > d_o1_o2 + tol_tang {
            let a_tol2 = tol * tol;
            let a_x2 = 4.0 * (r * r - d_o1_o2 * d_o1_o2);
            if a_x2 > a_tol2 {
                b2_sol = !b2_sol;
            }
        }
        if d_o1_o2 > r_m_tol && !b2_sol {
            // if(dO1O2 > RmTol) {
            let dx = d_o1_o2;
            let mut dy = 0.0f64; //(RpTol*RpTol-dx*dx); //Patch !!!
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            binf1 = -d_alpha1;
            bsup1 = d_alpha1;
            binf2 = 0.0;
            bsup2 = 0.0;
            *nbsol = 1;
        }
        //--------------------------------------------------------------------
        //--           2 segments donnes par Inter Line avec Circle-  Circle+
        //--
        else {
            //------------------- Intersection Line Circle+ --------------------------
            let dx = d_o1_o2;
            let mut dy = r * r - dx * dx; //(RpTol*RpTol-dx*dx); //Patch !!!
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };

            d_alpha1 = dy.atan2(dx);
            binf1 = -d_alpha1;
            bsup2 = d_alpha1; //--  |...?     ?...|   Sur C1

            //------------------ Intersection Line Circle-  -------------------------
            let mut dy2 = r * r - dx * dx; //(RmTol*RmTol-dx*dx); //Patch !!!
            dy2 = if dy2 >= 0.0 { dy2.sqrt() } else { 0.0 };
            d_alpha1 = dy2.atan2(dx);

            binf2 = d_alpha1;
            bsup1 = -d_alpha1; //--  |...x     x...|   Sur C1

            if (d_alpha1 * r) < (tol.max(tol_tang)) {
                bsup1 = bsup2;
                *nbsol = 1;
            } else {
                *nbsol = 2;
            }
        }
    }
    //--------------------------------------------------------------
    //-- Mise en forme des resultats :
    //--    Les calculs ont ete fait dans le repere x1,y1, (O1,O2)
    //--    On se ramene au repere propre a C1

    let mut d_angle1 = dir2d_angle(circle.x_dir, line.direction);

    // OCCT gp_Lin2d::Coefficients(a, b, c) (gp_Lin2d.hxx L91-96).
    let a = line.direction.y;
    let b = -line.direction.x;
    let c = -(a * line.origin.x + b * line.origin.y);

    let d = a * circle.center.x + b * circle.center.y + c;

    if d > 0.0 {
        d_angle1 += std::f64::consts::FRAC_PI_2;
    } else {
        d_angle1 -= std::f64::consts::FRAC_PI_2;
    }

    if d_angle1 < 0.0 {
        d_angle1 += PI_PPI;
    } else if d_angle1 > PI_PPI {
        d_angle1 -= PI_PPI;
    }

    binf1 += d_angle1;
    bsup1 += d_angle1;

    //-- par construction aucun des segments ne peut exceder PI
    //-- (permet de ne pas gerer trop de cas differents)

    // OCCT Circle.IsDirect(): the sign of the (XDirection ^ YDirection) z.
    let is_direct = circle.x_dir.x * circle.y_dir.y - circle.x_dir.y * circle.y_dir.x > 0.0;
    if !is_direct {
        let t = binf1;
        binf1 = bsup1;
        bsup1 = t;
        binf1 = -binf1;
        bsup1 = -bsup1;
    }

    c_int1.set_values(binf1, bsup1);
    if c_int1.length() > std::f64::consts::PI {
        c_int1.complement();
    }

    if *nbsol == 2 {
        binf2 += d_angle1;
        bsup2 += d_angle1;

        if !is_direct {
            let t = binf2;
            binf2 = bsup2;
            bsup2 = t;
            binf2 = -binf2;
            bsup2 = -bsup2;
        }

        c_int2.set_values(binf2, bsup2);
        if c_int2.length() > std::f64::consts::PI {
            c_int2.complement();
        }
    }
    //  Modified by Sergey KHROMOV - Thu Oct 26 17:51:05 2000 Begin
    else if c_int1.bsup > PI_PPI && c_int1.binf < PI_PPI {
        *nbsol = 2;
        binf2 = c_int1.binf;
        bsup2 = PI_PPI;
        binf1 = 0.0;
        c_int1.set_values(binf1, c_int1.bsup - PI_PPI);
        if c_int1.length() > std::f64::consts::PI {
            c_int1.complement();
        }
        c_int2.set_values(binf2, bsup2);
        if c_int2.length() > std::f64::consts::PI {
            c_int2.complement();
        }
    }
    //  Modified by Sergey KHROMOV - Thu Oct 26 17:51:13 2000 End
}

/// OCCT ProjectOnLAndIntersectWithLDomain (_1.cxx L360-443).
#[allow(clippy::too_many_arguments)]
fn project_on_l_and_intersect_with_l_domain(
    circle: &Circle2d,
    line: &Line2d,
    c_domain_and_res: &PeriodicInterval,
    l_domain: &Interval,
    circle_solution: &mut [PeriodicInterval; 4],
    line_solution: &mut [Interval; 4],
    nb_sol_total: &mut usize,
    ref_line_domain: &Res2dDomain,
    _ref_circ_domain: &Res2dDomain,
) {
    if c_domain_and_res.is_null() {
        return;
    }
    //-------------------------------------------------------------------------
    //--  On cherche l intervalle correspondant sur C2
    //--  Puis on intersecte l intervalle avec le domaine de C2
    //--  Enfin on cherche l intervalle correspondant sur C1
    //--

    let linf = elclib2d::line_parameter(
        line.origin,
        line.direction,
        elclib2d::circle_value(
            circle.center,
            circle.x_dir,
            circle.y_dir,
            circle.radius,
            c_domain_and_res.binf,
        ),
    );
    let lsup = elclib2d::line_parameter(
        line.origin,
        line.direction,
        elclib2d::circle_value(
            circle.center,
            circle.x_dir,
            circle.y_dir,
            circle.radius,
            c_domain_and_res.bsup,
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

        let mut cinf = c_domain_and_res.binf;
        let mut csup = c_domain_and_res.bsup;
        if cinf >= csup {
            cinf = c_domain_and_res.binf;
            csup = c_domain_and_res.bsup;
        }
        circle_solution[*nb_sol_total] = PeriodicInterval::new_ab(cinf, csup);
        if circle_solution[*nb_sol_total].length() > std::f64::consts::PI {
            circle_solution[*nb_sol_total].complement();
        }

        line_solution[*nb_sol_total] = l_inter_and_domain;
        *nb_sol_total += 1;
    }
}

impl IntConicConic {
    /// OCCT Perform(const gp_Lin2d& Line, const IntRes2d_Domain& LIG_Domain,
    /// const gp_Circ2d& Circle, const IntRes2d_Domain& CIRC_Domain,
    /// const double TolConf, const double Tol)
    /// (IntCurve_IntConicConic_1.cxx L2236-2652).
    pub fn perform_line_circle(
        &mut self,
        line: &Line2d,
        lig_domain: &Res2dDomain,
        circle: &Circle2d,
        circ_domain: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let the_reversed_parameters = self.base.reversed_parameters();
        self.base.reset_fields();
        self.base.set_reversed_parameters(the_reversed_parameters);

        let mut nbsol = 0i32;
        let mut c_int1 = PeriodicInterval::new();
        let mut c_int2 = PeriodicInterval::new();

        line_circle_geometric_intersection(
            line,
            circle,
            tol_conf,
            tol,
            &mut c_int1,
            &mut c_int2,
            &mut nbsol,
        );

        self.base.done = true;

        if nbsol == 0 {
            //-- Pas de solutions
            return;
        }

        //  Modified by Sergey KHROMOV - Mon Dec 18 11:13:18 2000 Begin
        if nbsol == 2 && c_int2.bsup == c_int1.binf + PI_PPI {
            let first_bound = circ_domain.first_parameter();
            let last_bound = circ_domain.last_parameter();
            let first_tol = circ_domain.first_tolerance();
            let last_tol = circ_domain.last_tolerance();
            if c_int1.binf == 0.0 && first_bound - first_tol > c_int1.bsup {
                nbsol = 1;
                c_int1.set_values(c_int2.binf, c_int2.bsup);
            } else if c_int2.bsup == PI_PPI && last_bound + last_tol < c_int2.binf {
                nbsol = 1;
            }
        }
        //  Modified by Sergey KHROMOV - Mon Dec 18 11:13:20 2000 End

        let mut c_domain = PeriodicInterval::from_domain(circ_domain);
        let mut deltat = c_domain.bsup - c_domain.binf;
        while c_domain.binf >= PI_PPI {
            c_domain.binf -= PI_PPI;
        }
        while c_domain.binf < 0.0 {
            c_domain.binf += PI_PPI;
        }
        c_domain.bsup = c_domain.binf + deltat;

        //------------------------------------------------------------
        //-- Ajout : Jeudi 28 mars 96
        //-- On agrandit artificiellement les domaines
        let mut binf_modif = c_domain.binf;
        let mut bsup_modif = c_domain.bsup;
        binf_modif -= circ_domain.first_tolerance() / circle.radius;
        bsup_modif += circ_domain.last_tolerance() / circle.radius;
        deltat = bsup_modif - binf_modif;
        if deltat <= PI_PPI {
            c_domain.binf = binf_modif;
            c_domain.bsup = bsup_modif;
        } else {
            let mut t = PI_PPI - deltat;
            t *= 0.5;
            c_domain.binf = binf_modif + t;
            c_domain.bsup = bsup_modif - t;
        }
        deltat = c_domain.bsup - c_domain.binf;
        while c_domain.binf >= PI_PPI {
            c_domain.binf -= PI_PPI;
        }
        while c_domain.binf < 0.0 {
            c_domain.binf += PI_PPI;
        }
        c_domain.bsup = c_domain.binf + deltat;
        //-- ------------------------------------------------------------

        let l_domain = Interval::from_domain(lig_domain);

        let mut nb_sol_total = 0usize;

        let mut solution_circle = [PeriodicInterval::new(); 4];
        let mut solution_line = [Interval::new(); 4];

        //----------------------------------------------------------------------
        //----------- Traitement du premier intervalle Geometrique  CInt1   ----
        //----------------------------------------------------------------------
        //-- NbSolTotal est incremente a chaque Intervalle solution.
        //-- On stocke les intervalles dans les tableaux : SolutionCircle[4]
        //--                                            et SolutionLine[4]
        //-- des Exemples faciles donnent 3 Intersections
        //-- des Problemes numeriques peuvent peut etre en donner 4 ??????
        //--
        let mut c_domain_and_res = c_domain.first_intersection(&mut c_int1);

        project_on_l_and_intersect_with_l_domain(
            circle,
            line,
            &c_domain_and_res,
            &l_domain,
            &mut solution_circle,
            &mut solution_line,
            &mut nb_sol_total,
            lig_domain,
            circ_domain,
        );

        c_domain_and_res = c_domain.second_intersection(&c_int1);

        project_on_l_and_intersect_with_l_domain(
            circle,
            line,
            &c_domain_and_res,
            &l_domain,
            &mut solution_circle,
            &mut solution_line,
            &mut nb_sol_total,
            lig_domain,
            circ_domain,
        );

        //----------------------------------------------------------------------
        //----------- Traitement du second intervalle Geometrique   C1_Int2 ----
        //----------------------------------------------------------------------
        if nbsol == 2 {
            // OCCT passes CInt2 by non-const reference: FirstIntersection
            // slides it in place, and the shifted value is what
            // SecondIntersection reads below.
            let c_domain_and_res = c_domain.first_intersection(&mut c_int2);

            project_on_l_and_intersect_with_l_domain(
                circle,
                line,
                &c_domain_and_res,
                &l_domain,
                &mut solution_circle,
                &mut solution_line,
                &mut nb_sol_total,
                lig_domain,
                circ_domain,
            );

            //--------------------------------------------------------------------
            let c_domain_and_res = c_domain.second_intersection(&c_int1);

            project_on_l_and_intersect_with_l_domain(
                circle,
                line,
                &c_domain_and_res,
                &l_domain,
                &mut solution_circle,
                &mut solution_line,
                &mut nb_sol_total,
                lig_domain,
                circ_domain,
            );
        }

        //----------------------------------------------------------------------
        //-- Calcul de toutes les transitions et Positions.
        //--
        //-- On determine si des intervalles sont reduit a des points
        //--      ( Rayon * Intervalle.Length()    <    TolConf   )   ### Modif 19 Nov Tol-->TolConf
        //--
        let r = circle.radius;
        let mut max_tol = tol_conf;
        if max_tol < tol {
            max_tol = tol;
        }
        if max_tol < 1.0e-10 {
            max_tol = 1.0e-10;
        }

        for item in solution_circle.iter_mut().enumerate().take(nb_sol_total) {
            let (i, sol) = item;
            if (r * sol.length()) < max_tol && (solution_line[i].length()) < max_tol {
                let mut t = (sol.binf + sol.bsup) * 0.5;
                sol.binf = t;
                sol.bsup = t;

                t = (solution_line[i].binf + solution_line[i].bsup) * 0.5;
                solution_line[i].binf = t;
                solution_line[i].bsup = t;
            }
        }
        //----------------------------------------------------------------------
        //-- Traitement des intervalles (ou des points obtenus)
        //--
        if nb_sol_total > 0 {
            let mut p1a;
            let mut p2a;
            let mut p1b;
            let mut p2b;
            let mut tan1;
            let mut norm1;
            let norm2 = DVec2::ZERO;
            let mut t1a = Transition::empty();
            let mut t2a = Transition::empty();
            let mut t1b = Transition::empty();
            let mut t2b = Transition::empty();
            let mut pos1a;
            let mut pos1b;
            let mut pos2a;
            let mut pos2b;

            let (pa, ta) =
                elclib2d::circle_d1(circle.center, circle.x_dir, circle.y_dir, r, solution_circle[0].binf);
            p1a = pa;
            tan1 = ta;
            let (pb, tb) = elclib2d::line_d1(line.origin, line.direction, solution_line[0].binf);
            p2a = pb;
            let mut tan2 = tb;
            p1b = DVec2::ZERO;
            p2b = DVec2::ZERO;
            norm1 = DVec2::ZERO;
            pos1a = Position::Middle;
            pos1b = Position::Middle;
            pos2a = Position::Middle;
            pos2b = Position::Middle;

            let is_opposite = tan1.dot(tan2) < 0.0;

            for i in 0..nb_sol_total {
                //-- 7 aout 97
                //-- On recentre Bin et Bsup de facon a avoir une portion commune avec CIRC_Domain
                let mut p1 = solution_circle[i].binf;
                let mut p2 = solution_circle[i].bsup;
                let q1 = circ_domain.first_parameter();
                let q2 = circ_domain.last_parameter();
                //--          |------ CircDomain ------|   [-- Sol --]
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

                solution_circle[i].binf = p1;
                solution_circle[i].bsup = p2;

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
                //-- Si les parametres sur le cercle sont en premier
                //-- On doit retourner ces parametres dans l ordre croissant
                //---------------------------------------------------------------
                if linf > lsup {
                    let t = solution_circle[i].binf;
                    solution_circle[i].binf = solution_circle[i].bsup;
                    solution_circle[i].bsup = t;

                    let t = linf;
                    linf = lsup;
                    lsup = t;
                }

                let (pa, ta, na) = elclib2d::circle_d2(
                    circle.center,
                    circle.x_dir,
                    circle.y_dir,
                    r,
                    solution_circle[i].binf,
                );
                p1a = pa;
                tan1 = ta;
                norm1 = na;
                let (pb, tb) = elclib2d::line_d1(line.origin, line.direction, linf);
                p2a = pb;
                tan2 = tb;

                determine_position(&mut pos1a, circ_domain, p1a, solution_circle[i].binf);
                determine_position(&mut pos2a, lig_domain, p2a, linf);
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
                let cinf;
                if pos1a == Position::End {
                    cinf = circ_domain.last_parameter();
                    p1a = circ_domain.last_point();
                    linf = elclib2d::line_parameter(line.origin, line.direction, p1a);

                    let (pa, ta, na) =
                        elclib2d::circle_d2(circle.center, circle.x_dir, circle.y_dir, r, cinf);
                    p1a = pa;
                    tan1 = ta;
                    norm1 = na;
                    let (pb, tb) = elclib2d::line_d1(line.origin, line.direction, linf);
                    p2a = pb;
                    tan2 = tb;
                    determine_position(&mut pos1a, circ_domain, p1a, cinf);
                    determine_position(&mut pos2a, lig_domain, p2a, linf);
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
                    cinf = circ_domain.first_parameter();
                    p1a = circ_domain.first_point();
                    linf = elclib2d::line_parameter(line.origin, line.direction, p1a);

                    let (pa, ta, na) =
                        elclib2d::circle_d2(circle.center, circle.x_dir, circle.y_dir, r, cinf);
                    p1a = pa;
                    tan1 = ta;
                    norm1 = na;
                    let (pb, tb) = elclib2d::line_d1(line.origin, line.direction, linf);
                    p2a = pb;
                    tan2 = tb;
                    determine_position(&mut pos1a, circ_domain, p1a, cinf);
                    determine_position(&mut pos2a, lig_domain, p2a, linf);
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
                    cinf = normalize_on_circle_domain(solution_circle[i].binf, circ_domain);
                }

                let mut new_point1 = IntersectionPoint::empty();
                new_point1.set_values(p1a, linf, cinf, t2a, t1a, the_reversed_parameters);

                if (solution_line[i].length() + solution_circle[i].length()) > 0.0 {
                    let (pb1, tb1, nb1) = elclib2d::circle_d2(
                        circle.center,
                        circle.x_dir,
                        circle.y_dir,
                        r,
                        solution_circle[i].bsup,
                    );
                    p1b = pb1;
                    tan1 = tb1;
                    norm1 = nb1;
                    let (pb2, tb2) = elclib2d::line_d1(line.origin, line.direction, lsup);
                    p2b = pb2;
                    tan2 = tb2;

                    determine_position(&mut pos1b, circ_domain, p1b, solution_circle[i].bsup);
                    determine_position(&mut pos2b, lig_domain, p2b, lsup);
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
                    let csup;
                    if pos1b == Position::End {
                        csup = circ_domain.last_parameter();
                        p1b = circ_domain.last_point();
                        lsup = elclib2d::line_parameter(line.origin, line.direction, p1b);
                        let (pa, ta, na) =
                            elclib2d::circle_d2(circle.center, circle.x_dir, circle.y_dir, r, csup);
                        p1b = pa;
                        tan1 = ta;
                        norm1 = na;
                        let (pb, tb) = elclib2d::line_d1(line.origin, line.direction, lsup);
                        p2b = pb;
                        tan2 = tb;

                        determine_position(&mut pos1b, circ_domain, p1b, csup);
                        determine_position(&mut pos2b, lig_domain, p2b, lsup);
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
                        csup = circ_domain.first_parameter();
                        p1b = circ_domain.first_point();
                        lsup = elclib2d::line_parameter(line.origin, line.direction, p1b);
                        let (pa, ta, na) =
                            elclib2d::circle_d2(circle.center, circle.x_dir, circle.y_dir, r, csup);
                        p1b = pa;
                        tan1 = ta;
                        norm1 = na;
                        let (pb, tb) = elclib2d::line_d1(line.origin, line.direction, lsup);
                        p2b = pb;
                        tan2 = tb;

                        determine_position(&mut pos1b, circ_domain, p1b, csup);
                        determine_position(&mut pos2b, lig_domain, p2b, lsup);
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
                        csup = normalize_on_circle_domain(solution_circle[i].bsup, circ_domain);
                    }

                    let mut new_point2 = IntersectionPoint::empty();
                    new_point2.set_values(p1b, lsup, csup, t2b, t1b, the_reversed_parameters);

                    if ((csup - cinf).abs() * r > max_tol && (lsup - linf).abs() > max_tol)
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
                        if pos1a != Position::Middle || pos2a != Position::Middle {
                            self.base.append_point(&new_point1);
                        }
                        if pos1b != Position::Middle || pos2b != Position::Middle {
                            self.base.append_point(&new_point2);
                        }
                    }
                } else {
                    self.base.append_point(&new_point1);
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

    /// The X-axis segment [-2, 2] (parameter u in [0, 4]).
    fn x_axis_line() -> (Line2d, Res2dDomain) {
        let line = Line2d {
            origin: DVec2::new(-2.0, 0.0),
            direction: DVec2::new(1.0, 0.0),
        };
        let mut d = Res2dDomain::infinite();
        d.set_values_bounded(
            line.origin,
            0.0,
            1.0e-9,
            DVec2::new(2.0, 0.0),
            4.0,
            1.0e-9,
        );
        (line, d)
    }

    /// The unit circle centered at the origin (parameter u in [0, 2*PI],
    /// closed — the IntCurveCurveGen::ComputeDomain convention).
    fn unit_circle() -> (Circle2d, Res2dDomain) {
        let circle = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 1.0,
        };
        let mut d = Res2dDomain::infinite();
        d.set_values_bounded(
            DVec2::new(1.0, 0.0),
            0.0,
            1.0e-9,
            DVec2::new(1.0, 0.0),
            std::f64::consts::TAU,
            1.0e-9,
        );
        d.set_equivalent_parameters(0.0, std::f64::consts::TAU);
        (circle, d)
    }

    /// OCCT anchor: the X-axis crossing the unit circle at (-1, 0) and
    /// (1, 0) — with the center ON the line the Circle+/Circle- bands are
    /// both degenerate, so the two solutions collapse to two points
    /// (Perform(Lin, Circ) _1.cxx L2236-2652 through
    /// LineCircleGeometricIntersection).
    #[test]
    fn line_circle_two_crossings() {
        let (line, dl) = x_axis_line();
        let (circle, dc) = unit_circle();

        let mut ic = IntConicConic::new();
        ic.perform_line_circle(&line, &dl, &circle, &dc, 1.0e-9, 1.0e-9);

        assert!(ic.base.done);
        assert_eq!(ic.base.nb_points(), 2, "nb_pt={}", ic.base.nb_points());
        assert_eq!(ic.base.nb_segments(), 0, "nb_seg={}", ic.base.nb_segments());

        // The two crossing points (-1, 0) and (1, 0).
        let p1 = ic.base.point(1).value();
        let p2 = ic.base.point(2).value();
        assert!(p1.y.abs() < 1.0e-6 && p2.y.abs() < 1.0e-6, "{:?} {:?}", p1, p2);
        let xs = [p1.x, p2.x];
        assert!(
            (xs[0] + 1.0).abs() < 1.0e-6 && (xs[1] - 1.0).abs() < 1.0e-6,
            "xs={:?}",
            xs
        );
        // Line parameters 1 and 3 (origin -2, unit direction).
        let ls = [
            ic.base.point(1).param_on_first(),
            ic.base.point(2).param_on_first(),
        ];
        assert!(
            (ls[0] - 1.0).abs() < 1.0e-6 && (ls[1] - 3.0).abs() < 1.0e-6,
            "ls={:?}",
            ls
        );
    }

    /// OCCT anchor: the tangent line y = 1 touches the unit circle at
    /// (0, 1) — the degenerate solution collapses to a single point with
    /// the line parameter 2 and the circle parameter PI/2, and the TOUCH
    /// transitions of Determine_Transition_LC.
    #[test]
    fn line_circle_tangent_point() {
        let circle = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 1.0,
        };
        let mut dc = Res2dDomain::infinite();
        dc.set_values_bounded(
            DVec2::new(1.0, 0.0),
            0.0,
            1.0e-9,
            DVec2::new(1.0, 0.0),
            std::f64::consts::TAU,
            1.0e-9,
        );
        dc.set_equivalent_parameters(0.0, std::f64::consts::TAU);

        let line = Line2d {
            origin: DVec2::new(-2.0, 1.0),
            direction: DVec2::new(1.0, 0.0),
        };
        let mut dl = Res2dDomain::infinite();
        dl.set_values_bounded(line.origin, 0.0, 1.0e-9, DVec2::new(2.0, 1.0), 4.0, 1.0e-9);

        let mut ic = IntConicConic::new();
        ic.perform_line_circle(&line, &dl, &circle, &dc, 1.0e-9, 1.0e-9);

        assert!(ic.base.done);
        assert_eq!(ic.base.nb_points(), 1, "nb_pt={}", ic.base.nb_points());
        assert_eq!(ic.base.nb_segments(), 0, "nb_seg={}", ic.base.nb_segments());

        let p = ic.base.point(1);
        let v = p.value();
        assert!(v.x.abs() < 1.0e-6, "x={}", v.x);
        assert!((v.y - 1.0).abs() < 1.0e-6, "y={}", v.y);
        // The line parameter of the tangency (origin x = -2 -> u = 2).
        assert!((p.param_on_first() - 2.0).abs() < 1.0e-6, "l={}", p.param_on_first());
        // The circle parameter of the tangency (PI/2), recadre into the domain.
        let c = p.param_on_second();
        let pi_2 = std::f64::consts::FRAC_PI_2;
        assert!(
            (c - pi_2).abs() < 1.0e-6 || (c - pi_2 - std::f64::consts::TAU).abs() < 1.0e-6,
            "c={}",
            c
        );
        // TOUCH: both transitions share the tangent type (Determine_Transition_LC).
        assert_eq!(
            p.transition_of_first().transition_type(),
            p.transition_of_second().transition_type()
        );
    }
}
