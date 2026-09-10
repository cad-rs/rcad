// OCCT IntCurve_IntConicConic::Perform(gp_Circ2d, gp_Circ2d) (TKGeomAlgo
// IntCurve) — the circle/circle closed-form intersection and the two
// file-static helpers it consumes from IntCurve_IntConicConic_1.cxx.
//
// 1:1 translation of:
//   - ProjectOnC2AndIntersectWithC2Domain       (_1.cxx L58-150)
//   - CircleCircleGeometricIntersection         (_1.cxx L154-356)
//   - Perform(gp_Circ2d, gp_Circ2d)             (_1.cxx L807-1206)
// plus the gp helpers consumed by the body (gp_Circ2d::IsDirect
// gp_Circ2d.hxx L181-184, gp_Circ2d::Reversed gp_Circ2d.hxx L288-295,
// gp_Vec2d::Angle gp_Vec2d.cxx L47-92).

use glam::DVec2;
use rcad_kernel::geom::Circle2d;

use super::geom2d_int::elclib2d;
use super::int_conic_conic::IntConicConic;
use super::int_conic_conic_tool::{
    determine_transition_lc, normalize_on_circle_domain, GP_RESOLUTION, PeriodicInterval, PI_PPI,
};
use super::int_imp_par_gen::determine_position;
use super::int_res2d::{
    Domain as Res2dDomain, IntersectionPoint, IntersectionSegment, Position, Transition,
};

/// `std::nextafter(x, to)` (not yet stable in Rust std) — the representable
/// value adjacent to `x` in the direction of `to`.
fn next_after(x: f64, to: f64) -> f64 {
    if x.is_nan() || to.is_nan() {
        return f64::NAN;
    }
    if x == to {
        return x;
    }
    if x == 0.0 {
        return if to > 0.0 {
            f64::MIN_POSITIVE
        } else {
            -f64::MIN_POSITIVE
        };
    }
    let bits = x.to_bits();
    let next = if (to > x) == (x > 0.0) {
        bits + 1
    } else {
        bits - 1
    };
    f64::from_bits(next)
}

/// OCCT gp_Circ2d::IsDirect() (gp_Circ2d.hxx L181-184): the sign of the
/// crossed XDirection ^ YDirection.
fn circle_is_direct(circle: &Circle2d) -> bool {
    circle.x_dir.x * circle.y_dir.y - circle.x_dir.y * circle.y_dir.x >= 0.0
}

/// OCCT gp_Circ2d::Reversed() (gp_Circ2d.hxx L288-295): the same circle with
/// the reversed YDirection.
fn circle_reversed(circle: &Circle2d) -> Circle2d {
    let mut r = *circle;
    r.y_dir = -r.y_dir;
    r
}

/// OCCT gp_Vec2d::Angle(Other) (gp_Vec2d.cxx L47-92) — the angle from self
/// to other in [-PI, PI], with the magnitude normalization and the
/// acos/asin precision switching. OCCT raises gp_VectorWithNullMagnitude
/// for null vectors; both call sites pass unit directions or the
/// gp::Resolution() checked AxeO1O2 vector, so the throw is unreachable.
fn vec2d_angle(v1: DVec2, v2: DVec2) -> f64 {
    let a_norm = v1.length();
    let an_other_norm = v2.length();
    let a_d = a_norm * an_other_norm;
    let a_cosinus = v1.dot(v2) / a_d;
    let a_sinus = (v1.x * v2.y - v1.y * v2.x) / a_d;

    // Use M_SQRT1_2 (1/sqrt(2) approximately 0.7071067811865476) for better readability and precision
    const A_COS_45_DEG: f64 = std::f64::consts::FRAC_1_SQRT_2;

    if a_cosinus > -A_COS_45_DEG && a_cosinus < A_COS_45_DEG {
        // For angles near +/-90 degrees, use acos for better precision
        if a_sinus > 0.0 {
            a_cosinus.acos()
        } else {
            -a_cosinus.acos()
        }
    } else {
        // For angles near 0 degrees or +/-180 degrees, use asin for better precision
        if a_cosinus > 0.0 {
            a_sinus.asin()
        } else if a_sinus > 0.0 {
            std::f64::consts::PI - a_sinus.asin()
        } else {
            -std::f64::consts::PI - a_sinus.asin()
        }
    }
}

/// OCCT ProjectOnC2AndIntersectWithC2Domain (_1.cxx L58-150).
#[allow(clippy::too_many_arguments)]
fn project_on_c2_and_intersect_with_c2_domain(
    circle1: &Circle2d,
    circle2: &Circle2d,
    c1_domain_and_res: &PeriodicInterval,
    domain_c2: &PeriodicInterval,
    solution_c1: &mut [PeriodicInterval; 4],
    solution_c2: &mut [PeriodicInterval; 4],
    nb_sol_total: &mut usize,
    ident_circles: bool,
) {
    if c1_domain_and_res.is_null() {
        return;
    }
    //-------------------------------------------------------------------------
    //--  On cherche l intervalle correspondant sur C2
    //--  Puis on intersecte l intervalle avec le domaine de C2
    //--  Enfin, on cherche l intervalle correspondant sur C1
    //--
    let c2inf = elclib2d::circle_parameter(
        circle2.center,
        circle2.x_dir,
        circle2.y_dir,
        elclib2d::circle_value(
            circle1.center,
            circle1.x_dir,
            circle1.y_dir,
            circle1.radius,
            c1_domain_and_res.binf,
        ),
    );
    let c2sup = elclib2d::circle_parameter(
        circle2.center,
        circle2.x_dir,
        circle2.y_dir,
        elclib2d::circle_value(
            circle1.center,
            circle1.x_dir,
            circle1.y_dir,
            circle1.radius,
            c1_domain_and_res.bsup,
        ),
    );

    let mut c2_inter = PeriodicInterval::new_ab(c2inf, c2sup);

    if !ident_circles {
        if c2_inter.length() > std::f64::consts::PI {
            c2_inter.complement();
        }
    } else {
        let mut c2inf = c2inf;
        let mut c2sup = c2sup;
        if c2sup <= c2inf {
            c2sup += PI_PPI;
        }
        if c2inf >= PI_PPI {
            c2sup -= PI_PPI;
            c2inf -= PI_PPI;
        }
        c2_inter.binf = c2inf;
        c2_inter.bsup = c2sup; //-- Verifier la longueur de l'intervalle sur C2
        c2_inter.bsup = c2inf + c1_domain_and_res.bsup - c1_domain_and_res.binf;
    }

    for i in 0..2 {
        // OCCT FirstIntersection takes C2Inter by non-const reference: the
        // slide is in place, and the shifted value is what
        // SecondIntersection reads for i == 1.
        let c2_inter_and_domain = if i == 0 {
            domain_c2.first_intersection(&mut c2_inter)
        } else {
            domain_c2.second_intersection(&c2_inter)
        };

        if !c2_inter_and_domain.is_null() {
            let c1inf = elclib2d::circle_parameter(
                circle1.center,
                circle1.x_dir,
                circle1.y_dir,
                elclib2d::circle_value(
                    circle2.center,
                    circle2.x_dir,
                    circle2.y_dir,
                    circle2.radius,
                    c2_inter_and_domain.binf,
                ),
            );
            let c1sup = elclib2d::circle_parameter(
                circle1.center,
                circle1.x_dir,
                circle1.y_dir,
                elclib2d::circle_value(
                    circle2.center,
                    circle2.x_dir,
                    circle2.y_dir,
                    circle2.radius,
                    c2_inter_and_domain.bsup,
                ),
            );

            solution_c1[*nb_sol_total] = PeriodicInterval::new_ab(c1inf, c1sup);
            if !ident_circles {
                if solution_c1[*nb_sol_total].length() > std::f64::consts::PI {
                    solution_c1[*nb_sol_total].complement();
                }
            } else {
                if solution_c1[*nb_sol_total].bsup <= solution_c1[*nb_sol_total].binf {
                    solution_c1[*nb_sol_total].bsup += PI_PPI;
                }
                if solution_c1[*nb_sol_total].binf >= PI_PPI {
                    solution_c1[*nb_sol_total].binf -= PI_PPI;
                    solution_c1[*nb_sol_total].bsup -= PI_PPI;
                }
            }
            solution_c2[*nb_sol_total] = c2_inter_and_domain;
            *nb_sol_total += 1;
        }
    }
}

/// OCCT CircleCircleGeometricIntersection (_1.cxx L154-356).
fn circle_circle_geometric_intersection(
    c1: &Circle2d,
    c2: &Circle2d,
    tol: f64,
    tol_tang: f64,
    c1_res1: &mut PeriodicInterval,
    c1_res2: &mut PeriodicInterval,
    nbsol: &mut i32,
) {
    // OCCT: double C1_binf1, C1_binf2 = 0, C1_bsup1, C1_bsup2 = 0;
    let mut c1_binf1;
    let mut c1_binf2 = 0.0f64;
    let mut c1_bsup1;
    let mut c1_bsup2 = 0.0f64;
    let d_o1_o2 = c1.center.distance(c2.center);
    let r1 = c1.radius;
    let r2 = c2.radius;
    let abs_r1m_r2 = (r1 - r2).abs();
    //----------------------------------------------------------------
    if d_o1_o2 > (r1 + r2 + tol) {
        if d_o1_o2 > (r1 + r2 + tol_tang) {
            *nbsol = 0;
            return;
        } else {
            c1_binf1 = 0.0;
            c1_bsup1 = 0.0;
            *nbsol = 1;
        }
    }
    //----------------------------------------------------------------
    else if d_o1_o2 <= tol && abs_r1m_r2 <= tol {
        *nbsol = 3;
        return;
    } else {
        //----------------------------------------------------------------
        let r1p_r2 = r1 + r2;
        let r1p_tol = r1 + tol;
        let r1m_tol = r1 - tol;
        //    double R1R1=R1*R1;
        let r2r2 = r2 * r2;
        let r1p_tol_r1p_tol = r1p_tol * r1p_tol;
        let r1m_tol_r1m_tol = r1m_tol * r1m_tol;
        let d_o1_o2d_o1_o2 = d_o1_o2 * d_o1_o2;
        let mut d_alpha1;
        //--------------------------------------------------------------- Cas
        //-- C2 coupe le cercle C1+ (=C(x1,y1,R1+Tol))
        //--            1 seul segment donne par Inter C2 C1+
        //--
        if d_o1_o2 > r1p_r2 - tol {
            let dx = (r1p_tol_r1p_tol + d_o1_o2d_o1_o2 - r2r2) / (d_o1_o2 + d_o1_o2);
            let mut dy = r1p_tol_r1p_tol - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            c1_binf1 = -d_alpha1;
            c1_bsup1 = d_alpha1;
            *nbsol = 1;
        }
        //--------------------------------------------------------------------
        //--           2 segments donnes par Inter C2 avec C1- C1 C1+
        //-- Seul le signe de dx change si dO1O2 < std::max(R1,R2)
        //--
        else if d_o1_o2 > abs_r1m_r2 - tol {
            // -- +
            //------------------- Intersection C2 C1+ --------------------------
            let mut dx = (r1p_tol_r1p_tol + d_o1_o2d_o1_o2 - r2r2) / (d_o1_o2 + d_o1_o2);
            let mut dy = r1p_tol_r1p_tol - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };

            d_alpha1 = dy.atan2(dx);
            c1_binf1 = -d_alpha1;
            c1_bsup2 = d_alpha1; //--  |...?     ?...|   Sur C1

            //------------------ Intersection C2 C1- -------------------------
            dx = (r1m_tol_r1m_tol + d_o1_o2d_o1_o2 - r2r2) / (d_o1_o2 + d_o1_o2);
            dy = r1m_tol_r1m_tol - dx * dx;
            dy = if dy >= 0.0 { dy.sqrt() } else { 0.0 };
            d_alpha1 = dy.atan2(dx);

            c1_binf2 = d_alpha1;
            c1_bsup1 = -d_alpha1; //--  |...x     x...|   Sur C1
            *nbsol = 2;
            //------------------------------
            //-- Les 2 intervalles sont ils
            //-- en fait un seul inter ?
            //--
            if dy == 0.0 {
                //-- Les 2 bornes internes sont identiques
                c1_bsup1 = c1_bsup2;
                *nbsol = 1;
            } else {
                if c1_binf1 > c1_bsup1 {
                    d_alpha1 = c1_binf1;
                    c1_binf1 = c1_bsup1;
                    c1_bsup1 = d_alpha1;
                }
                if c1_binf2 > c1_bsup2 {
                    d_alpha1 = c1_binf2;
                    c1_binf2 = c1_bsup2;
                    c1_bsup2 = d_alpha1;
                }
                if ((c1_binf1 <= c1_bsup2) && (c1_binf1 >= c1_binf2))
                    || ((c1_bsup1 <= c1_bsup2) && (c1_bsup1 >= c1_binf2))
                {
                    if c1_binf1 > c1_binf2 {
                        c1_binf1 = c1_binf2;
                    }
                    if c1_binf1 > c1_bsup2 {
                        c1_binf1 = c1_bsup2;
                    }
                    if c1_bsup1 < c1_binf2 {
                        c1_bsup1 = c1_binf2;
                    }
                    if c1_bsup1 < c1_bsup2 {
                        c1_bsup1 = c1_bsup2;
                    }
                    *nbsol = 1;
                }
            }
        }
        //--------------------------------------------------------------
        else {
            if (d_o1_o2 > abs_r1m_r2 - tol_tang) && (abs_r1m_r2 - tol_tang) > 0.0 {
                c1_binf1 = 0.0;
                c1_bsup1 = 0.0;
                *nbsol = 1;
            } else {
                *nbsol = 0;
                return;
            }
        }
    }

    //----------------------------------------------------------------
    //-- Mise en forme des resultats :
    //--    Les calculs ont ete fait dans le repere x1,y1, (O1,O2)
    //--    On se ramene au repere propre a C1

    let axe1 = c1.x_dir; // gp_Vec2d Axe1 = C1.XAxis().Direction();
    let axe_o1_o2 = c2.center - c1.center; // gp_Vec2d(C1.Location(), C2.Location())

    let mut d_angle1;
    if axe_o1_o2.length() <= GP_RESOLUTION {
        d_angle1 = vec2d_angle(axe1, c2.x_dir);
    } else {
        d_angle1 = vec2d_angle(axe1, axe_o1_o2);
    }

    if !circle_is_direct(c1) {
        d_angle1 = -d_angle1;
    }

    c1_binf1 += d_angle1;
    c1_bsup1 += d_angle1;

    //-- par construction aucun des segments ne peut exceder PI
    //-- (permet de ne pas gerer trop de cas differents)

    c1_res1.set_values(c1_binf1, c1_bsup1);
    if c1_res1.length() > std::f64::consts::PI {
        c1_res1.complement();
    }

    if *nbsol == 2 {
        c1_binf2 += d_angle1;
        c1_bsup2 += d_angle1;
        c1_res2.set_values(c1_binf2, c1_bsup2);
        if c1_res2.length() > std::f64::consts::PI {
            c1_res2.complement();
        }
    } else {
        c1_res2.set_null();
    }
}

impl IntConicConic {
    /// OCCT Perform(const gp_Circ2d& Circle1, const IntRes2d_Domain& DomainCirc1,
    /// const gp_Circ2d& _Circle2, const IntRes2d_Domain& _DomainCirc2,
    /// const double TolConf, const double Tol)
    /// (IntCurve_IntConicConic_1.cxx L807-1206).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle_circle(
        &mut self,
        circle1: &Circle2d,
        domain_circ1: &Res2dDomain,
        _circle2: &Circle2d,
        _domain_circ2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        //-- TRES TRES MAL FAIT    A REPRENDRE UN JOUR ....   (lbr Octobre 98)
        let mut circle2 = *_circle2;
        let mut domain_circ2 = _domain_circ2.clone();
        let mut indirect_circles = false;
        if circle_is_direct(circle1) != circle_is_direct(_circle2) {
            indirect_circles = true;
            circle2 = circle_reversed(_circle2);
            // OCCT SetValues(Pnt1, Par1, Tol1, Pnt2, Par2, Tol2) — the Rust
            // Domain names the 6-argument overload set_values_bounded.
            domain_circ2.set_values_bounded(
                _domain_circ2.last_point(),
                PI_PPI - _domain_circ2.last_parameter(),
                _domain_circ2.last_tolerance(),
                _domain_circ2.first_point(),
                PI_PPI - _domain_circ2.first_parameter(),
                _domain_circ2.first_tolerance(),
            );
            domain_circ2.set_equivalent_parameters(0.0, PI_PPI);
        }

        self.base.reset_fields();
        let mut nbsol = 0i32;
        let mut c1_int1 = PeriodicInterval::new();
        let mut c1_int2 = PeriodicInterval::new();

        //------- Intersection sans tenir compte du domaine  ----> nbsol=0,1,2,3
        circle_circle_geometric_intersection(
            circle1,
            &circle2,
            tol_conf,
            tol,
            &mut c1_int1,
            &mut c1_int2,
            &mut nbsol,
        );
        self.base.done = true;

        if nbsol == 0 {
            //-- Pas de solutions
            return;
        }

        let mut c1_domain = PeriodicInterval::from_domain(domain_circ1);
        //-- On se ramene entre 0 et 2PI
        let mut deltat = c1_domain.bsup - c1_domain.binf;
        if deltat >= PI_PPI {
            // make deltat not including the upper limit
            deltat = next_after(PI_PPI, 0.0);
        }

        while c1_domain.binf >= PI_PPI {
            c1_domain.binf -= PI_PPI;
        }
        while c1_domain.binf < 0.0 {
            c1_domain.binf += PI_PPI;
        }

        c1_domain.bsup = c1_domain.binf + deltat;

        let mut c2_domain = PeriodicInterval::from_domain(&domain_circ2);
        deltat = c2_domain.bsup - c2_domain.binf;
        if deltat >= PI_PPI {
            deltat = next_after(PI_PPI, 0.0);
        }

        while c2_domain.binf >= PI_PPI {
            c2_domain.binf -= PI_PPI;
        }
        while c2_domain.binf < 0.0 {
            c2_domain.binf += PI_PPI;
        }

        c2_domain.bsup = c2_domain.binf + deltat;

        let mut ident_circles = false;

        if nbsol > 2 {
            //-- Les 2 cercles sont confondus a Tol pres
            c1_int1.set_values(0.0, PI_PPI);
            c1_int2.set_null();
            //---------------------------------------------------------------
            //-- Flag utilise pour specifier que les intervalles manipules
            //--   peuvent etre de longueur superieure a pi.
            //-- Pour des cercles non identiques, on a necessairement cette
            //--   condition sur les resultats de l intersection geometrique
            //--   ce qui permet de normaliser rapidement les intervalles.
            //--   ex: -1 4 -> longueur > PI
            //--        donc -1 4 devient  4 , 2*pi-1
            //---------------------------------------------------------------
            ident_circles = true;
        }

        let mut nb_sol_total = 0usize;
        let mut solution_c1 = [PeriodicInterval::new(); 4];
        let mut solution_c2 = [PeriodicInterval::new(); 4];

        //----------------------------------------------------------------------
        //----------- Traitement du premier intervalle Geometrique  C1_Int1 ----
        //----------------------------------------------------------------------
        //-- NbSolTotal est incremente a chaque Intervalle solution.
        //-- On stocke les intervalles dans les tableaux : SolutionC1(C2)
        //-- Dimensionnes a 4 elements.
        //-- des Exemples faciles donnent 3 Intersections
        //-- des Problemes numeriques peuvent en donner 4 ??????
        //--
        let mut c1_domain_and_res = c1_domain.first_intersection(&mut c1_int1);

        project_on_c2_and_intersect_with_c2_domain(
            circle1,
            &circle2,
            &c1_domain_and_res,
            &c2_domain,
            &mut solution_c1,
            &mut solution_c2,
            &mut nb_sol_total,
            ident_circles,
        );
        //----------------------------------------------------------------------
        //-- Seconde Intersection :  Par exemple :     2*PI-1  2*PI+1
        //--                         Intersecte avec     0.5   2*PI-0.5
        //--     Donne les intervalles : 0.5,1    et  2*PI-1,2*PI-0.5
        //--
        c1_domain_and_res = c1_domain.second_intersection(&c1_int1);

        project_on_c2_and_intersect_with_c2_domain(
            circle1,
            &circle2,
            &c1_domain_and_res,
            &c2_domain,
            &mut solution_c1,
            &mut solution_c2,
            &mut nb_sol_total,
            ident_circles,
        );

        //----------------------------------------------------------------------
        //----------- Traitement du second intervalle Geometrique   C1_Int2 ----
        //----------------------------------------------------------------------
        if nbsol == 2 {
            c1_domain_and_res = c1_domain.first_intersection(&mut c1_int2);

            project_on_c2_and_intersect_with_c2_domain(
                circle1,
                &circle2,
                &c1_domain_and_res,
                &c2_domain,
                &mut solution_c1,
                &mut solution_c2,
                &mut nb_sol_total,
                ident_circles,
            );
            //--------------------------------------------------------------------
            c1_domain_and_res = c1_domain.second_intersection(&c1_int2);

            project_on_c2_and_intersect_with_c2_domain(
                circle1,
                &circle2,
                &c1_domain_and_res,
                &c2_domain,
                &mut solution_c1,
                &mut solution_c2,
                &mut nb_sol_total,
                ident_circles,
            );
        }
        //----------------------------------------------------------------------
        //-- Calcul de toutes les transitions et Positions.
        //--
        //----------------------------------------------------------------------
        //-- On determine si des intervalles sont reduit a des points
        //--      ( Rayon * Intervalle.Length()    <    Tol   )
        //--
        let r1 = circle1.radius;
        let r2 = circle2.radius;
        let mut tol2 = tol + tol; //---- Pour eviter de toujours retourner
        // des segments
        if tol < 1e-10 {
            tol2 = 1e-10;
        }

        for i in 0..nb_sol_total {
            if (r1 * solution_c1[i].length() <= tol2) && (r2 * solution_c2[i].length() <= tol2) {
                let t = (solution_c1[i].binf + solution_c1[i].bsup) * 0.5;
                solution_c1[i].binf = t;
                solution_c1[i].bsup = t;

                let t = (solution_c2[i].binf + solution_c2[i].bsup) * 0.5;
                solution_c2[i].binf = t;
                solution_c2[i].bsup = t;
            }
        }

        //----------------------------------------------------------------------
        //-- Traitement des intervalles (ou des points obtenus)
        //--
        let mut p1a;
        let mut p1b;
        let mut p2a;
        let mut p2b;
        let mut tan1;
        let mut tan2;
        let mut norm1;
        let mut norm2;
        let mut t1a = Transition::empty();
        let mut t1b = Transition::empty();
        let mut t2a = Transition::empty();
        let mut t2b = Transition::empty();
        let mut pos1a = Position::Middle;
        let mut pos1b = Position::Middle;
        let mut pos2a = Position::Middle;
        let mut pos2b = Position::Middle;

        let is_opposite =
            (circle1.center - circle2.center).length_squared() > (r1 * r1 + r2 * r2);

        for i in 0..nb_sol_total {
            let mut c2inf = if is_opposite {
                solution_c2[i].bsup
            } else {
                solution_c2[i].binf
            };
            let mut c2sup = if is_opposite {
                solution_c2[i].binf
            } else {
                solution_c2[i].bsup
            };
            let c1tinf = solution_c1[i].binf;
            let c2tinf = c2inf;
            let mut c1inf = normalize_on_circle_domain(c1tinf, domain_circ1);
            c2inf = normalize_on_circle_domain(c2tinf, &domain_circ2);

            let mut is_out_of_range = false;
            if c1inf < domain_circ1.first_parameter() {
                if c1tinf < domain_circ1.first_parameter() {
                    c1inf = domain_circ1.first_parameter();
                    is_out_of_range = true;
                } else {
                    c1inf = c1tinf;
                }
            }

            if c1inf > domain_circ1.last_parameter() {
                if c1tinf > domain_circ1.last_parameter() {
                    c1inf = domain_circ1.last_parameter();
                    is_out_of_range = true;
                } else {
                    c1inf = c1tinf;
                }
            }

            if c2inf < domain_circ2.first_parameter() {
                if c2tinf < domain_circ2.first_parameter() {
                    c2inf = domain_circ2.first_parameter();
                    is_out_of_range = true;
                } else {
                    c2inf = c2tinf;
                }
            }

            if c2inf > domain_circ2.last_parameter() {
                if c2tinf > domain_circ2.last_parameter() {
                    c2inf = domain_circ2.last_parameter();
                    is_out_of_range = true;
                } else {
                    c2inf = c2tinf;
                }
            }

            if is_out_of_range {
                let (a_p1, _a_v11, _a_v12) =
                    elclib2d::circle_d2(circle1.center, circle1.x_dir, circle1.y_dir, r1, c1inf);
                let (a_p2, _a_v21, _a_v22) =
                    elclib2d::circle_d2(circle2.center, circle2.x_dir, circle2.y_dir, r2, c2inf);

                if a_p1.distance_squared(a_p2) > tol2 * tol2 {
                    // there are not any solutions in given parametric range.
                    continue;
                }
            }

            if indirect_circles {
                let (pa, ta, na) =
                    elclib2d::circle_d2(circle1.center, circle1.x_dir, circle1.y_dir, r1, c1inf);
                p1a = pa;
                tan1 = ta;
                norm1 = na;
                let (pa, ta, na) =
                    elclib2d::circle_d2(circle2.center, circle2.x_dir, circle2.y_dir, r2, c2inf);
                p2a = pa;
                tan2 = ta;
                norm2 = na;
                tan2 = -tan2; // Tan2.Reverse();

                determine_position(&mut pos1a, domain_circ1, p1a, c1inf);
                determine_position(&mut pos2a, _domain_circ2, p2a, PI_PPI - c2inf);
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

                let new_point1 =
                    IntersectionPoint::new(p1a, c1inf, PI_PPI - c2inf, t1a, t2a, false);

                if (solution_c1[i].length() > 0.0) || (solution_c2[i].length() > 0.0) {
                    //-- On traite un intervalle non reduit a un point
                    let mut c1sup = normalize_on_circle_domain(solution_c1[i].bsup, domain_circ1);
                    if c1sup < c1inf {
                        c1sup += PI_PPI;
                    }
                    c2sup = normalize_on_circle_domain(c2sup, &domain_circ2);

                    let (pa, ta, na) =
                        elclib2d::circle_d2(circle1.center, circle1.x_dir, circle1.y_dir, r1, c1sup);
                    p1b = pa;
                    tan1 = ta;
                    norm1 = na;
                    let (pa, ta, na) =
                        elclib2d::circle_d2(circle2.center, circle2.x_dir, circle2.y_dir, r2, c2sup);
                    p2b = pa;
                    tan2 = ta;
                    norm2 = na;
                    tan2 = -tan2; // Tan2.Reverse();

                    determine_position(&mut pos1b, domain_circ1, p1b, c1sup);
                    determine_position(&mut pos2b, _domain_circ2, p2b, PI_PPI - c2sup);
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

                    //--------------------------------------------------

                    if is_opposite {
                        if nbsol != 3 {
                            if c2inf < c2sup {
                                c2inf += PI_PPI;
                            }
                        }
                    } else {
                        if nbsol != 3 {
                            if c2sup < c2inf {
                                c2sup += PI_PPI;
                            }
                        }
                    }

                    let new_point2 =
                        IntersectionPoint::new(p1b, c1sup, PI_PPI - c2sup, t1b, t2b, false);
                    let new_seg =
                        IntersectionSegment::with_points(&new_point1, &new_point2, !is_opposite, false);
                    self.base.append_segment(&new_seg);
                } else {
                    self.base.append_point(&new_point1);
                }
            } else {
                let (pa, ta, na) =
                    elclib2d::circle_d2(circle1.center, circle1.x_dir, circle1.y_dir, r1, c1inf);
                p1a = pa;
                tan1 = ta;
                norm1 = na;
                let (pa, ta, na) =
                    elclib2d::circle_d2(circle2.center, circle2.x_dir, circle2.y_dir, r2, c2inf);
                p2a = pa;
                tan2 = ta;
                norm2 = na;

                determine_position(&mut pos1a, domain_circ1, p1a, c1inf);
                determine_position(&mut pos2a, &domain_circ2, p2a, c2inf);
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

                let new_point1 = IntersectionPoint::new(p1a, c1inf, c2inf, t1a, t2a, false);

                if (solution_c1[i].length() > 0.0) || (solution_c2[i].length() > 0.0) {
                    //-- On traite un intervalle non reduit a un point
                    let mut c1sup = normalize_on_circle_domain(solution_c1[i].bsup, domain_circ1);
                    if c1sup < c1inf {
                        c1sup += PI_PPI;
                    }
                    c2sup = normalize_on_circle_domain(c2sup, &domain_circ2);

                    let (pa, ta, na) =
                        elclib2d::circle_d2(circle1.center, circle1.x_dir, circle1.y_dir, r1, c1sup);
                    p1b = pa;
                    tan1 = ta;
                    norm1 = na;
                    let (pa, ta, na) =
                        elclib2d::circle_d2(circle2.center, circle2.x_dir, circle2.y_dir, r2, c2sup);
                    p2b = pa;
                    tan2 = ta;
                    norm2 = na;

                    determine_position(&mut pos1b, domain_circ1, p1b, c1sup);
                    determine_position(&mut pos2b, &domain_circ2, p2b, c2sup);
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

                    //--------------------------------------------------

                    if is_opposite {
                        if c2inf < c2sup {
                            c2inf += PI_PPI;
                        }
                    } else {
                        if c2sup < c2inf {
                            c2sup += PI_PPI;
                        }
                    }

                    let new_point2 = IntersectionPoint::new(p1b, c1sup, c2sup, t1b, t2b, false);
                    let new_seg =
                        IntersectionSegment::with_points(&new_point1, &new_point2, is_opposite, false);
                    self.base.append_segment(&new_seg);
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

    /// The unit circle centered at `center` with the full closed domain
    /// [0, 2*PI] (the IntCurveCurveGen::ComputeDomain convention).
    fn unit_circle_full_domain(center: DVec2) -> (Circle2d, Res2dDomain) {
        let r = 1.0;
        let circle = Circle2d {
            center,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: r,
        };
        let mut d = Res2dDomain::infinite();
        d.set_values_bounded(
            center + DVec2::new(r, 0.0),
            0.0,
            1.0e-9,
            center + DVec2::new(r, 0.0),
            std::f64::consts::TAU,
            1.0e-9,
        );
        d.set_equivalent_parameters(0.0, std::f64::consts::TAU);
        (circle, d)
    }

    /// OCCT anchor: the unit circles centered (0,0) and (1,0) cross at
    /// (0.5, +/-sqrt(3)/2) — Perform(Circ, Circ) _1.cxx L807-1206 through
    /// CircleCircleGeometricIntersection with TolConf = Tol = 1e-9 collapses
    /// each crossing (R * Interval.Length() <= 2*Tol) to a single point.
    /// Both argument orientations must produce the same geometry.
    #[test]
    fn circle_circle_two_crossings() {
        let (c1, d1) = unit_circle_full_domain(DVec2::new(0.0, 0.0));
        let (c2, d2) = unit_circle_full_domain(DVec2::new(1.0, 0.0));

        let sq3_2 = (3.0f64).sqrt() / 2.0;
        let expected = [
            DVec2::new(0.5, sq3_2),
            DVec2::new(0.5, -sq3_2),
        ];

        // Both orientations: (C1, C2) and the swapped (C2, C1).
        for (first, df, second, ds) in [
            (&c1, &d1, &c2, &d2),
            (&c2, &d2, &c1, &d1),
        ] {
            let mut ic = IntConicConic::new();
            ic.perform_circle_circle(first, df, second, ds, 1.0e-9, 1.0e-9);

            assert!(ic.base.done);
            assert!(
                ic.base.nb_points() == 2 && ic.base.nb_segments() == 0,
                "nb_pt={} nb_seg={}",
                ic.base.nb_points(),
                ic.base.nb_segments()
            );

            let pts = [
                ic.base.point(1).value(),
                ic.base.point(2).value(),
            ];
            for e in expected.iter() {
                let ok = pts.iter().any(|p| p.distance(*e) < 1.0e-6);
                assert!(ok, "missing {:?} in {:?}", e, pts);
            }

            // Parameters on the first circle of the orientation with C1 first:
            // the crossings at PI/3 and 5*PI/3.
            if std::ptr::eq(first, &c1) {
                let mut ps = [
                    ic.base.point(1).param_on_first(),
                    ic.base.point(2).param_on_first(),
                ];
                ps.sort_by(|a, b| a.partial_cmp(b).unwrap());
                assert!(
                    (ps[0] - std::f64::consts::FRAC_PI_3).abs() < 1.0e-6,
                    "ps={:?}",
                    ps
                );
                assert!(
                    (ps[1] - 5.0 * std::f64::consts::FRAC_PI_3).abs() < 1.0e-6,
                    "ps={:?}",
                    ps
                );
            }
        }
    }

    /// OCCT anchor: two identical circles are confounded to Tol pres — the
    /// geometric intersection answers nbsol = 3, IdentCircles is set and the
    /// result is a segment running along both circles.
    #[test]
    fn circle_circle_confounded() {
        let (c1, d1) = unit_circle_full_domain(DVec2::ZERO);
        let (c2, d2) = unit_circle_full_domain(DVec2::ZERO);

        let mut ic = IntConicConic::new();
        ic.perform_circle_circle(&c1, &d1, &c2, &d2, 1.0e-9, 1.0e-9);

        assert!(ic.base.done);
        assert!(
            ic.base.nb_segments() >= 1,
            "nb_seg={}",
            ic.base.nb_segments()
        );

        // The confounded segment spans the full circle on both curves.
        let seg = ic.base.segment(1);
        assert!(seg.has_first_point() && seg.has_last_point());
        let pf = seg.first_point().param_on_first();
        let pl = seg.last_point().param_on_first();
        assert!(
            (pl - pf).abs() > std::f64::consts::PI,
            "pf={} pl={}",
            pf,
            pl
        );
    }
}
