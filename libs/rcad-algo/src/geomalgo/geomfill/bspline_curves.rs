//! OCCT GeomFill_BSplineCurves (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_BSplineCurves.cxx (whole file L28-597).
//!
//! Mapping: `Handle(Geom_BSplineCurve)` -> rcad `BSplineCurve3` (flat knot
//! vector; (knots, mults) are recovered via `knots_mults()`), out-parameter
//! handles -> `&mut`, `Geom_BSplineSurface` -> rcad `BSplineSurface` (flat
//! knots; rcad has no periodic flag — the fillings produce non-periodic
//! surfaces, matching the OCCT ctor defaults).

use glam::DVec3;

use rcad_kernel::geom::{BSplineCurve3, BSplineSurface};
use rcad_kernel::math::bspl_lib::{
    insert_knots as bspl_insert_knots, prepare_insert_knots as bspl_prepare_insert_knots,
    reparametrize as bspl_reparametrize,
};

// OCCT Precision::Confusion() / PConfusion().
const CONFUSION: f64 = 1e-7;
const PCONFUSION: f64 = 1e-12;

/// OCCT GeomFill_FillingStyle (GeomFill_FillingStyle.hxx L27-29).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillingStyle {
    StretchStyle,
    CoonsStyle,
    CurvedStyle,
}

/// OCCT static Arrange (GeomFill_BSplineCurves.cxx L28-130) — reorders and
/// orients four boundary curves into a closed contour (CC1, CC2, CC3
/// reversed, CC4 reversed).  Returns false when the contour cannot be
/// chained.
#[allow(clippy::many_single_char_names)]
fn arrange(
    c1: &BSplineCurve3,
    c2: &BSplineCurve3,
    c3: &BSplineCurve3,
    c4: &BSplineCurve3,
    cc1: &mut BSplineCurve3,
    cc2: &mut BSplineCurve3,
    cc3: &mut BSplineCurve3,
    cc4: &mut BSplineCurve3,
    tol: f64,
) -> bool {
    let mut gc: [BSplineCurve3; 4] = [c1.clone(), c2.clone(), c3.clone(), c4.clone()];

    for i in 1..=3usize {
        let mut trouve = false;

        // OCCT: GC[] is a 0-based C array indexed by the 1-based i/j loop
        // variables — candidate GC[j] / previous GC[i-1] map to gc[j] /
        // gc[i-1] unchanged.
        //
        // search for a degenerated curve = point, which would match first
        let mut j = i;
        while j <= 3 && !trouve {
            let start = ops_start(&gc[j]);
            let end = ops_end(&gc[j]);
            if start.distance(end) < tol {
                // this is a degenerated line, does it match the last endpoint?
                if start.distance(ops_end(&gc[i - 1])) < tol {
                    gc.swap(i, j);
                    trouve = true;
                }
            }
            j += 1;
        }

        // if no degenerated curve matched, try an ordinary one as next curve
        if !trouve {
            let mut j = i;
            while j <= 3 && !trouve {
                if ops_start(&gc[j]).distance(ops_end(&gc[i - 1])) < tol {
                    gc.swap(i, j);
                    trouve = true;
                } else if ops_end(&gc[j]).distance(ops_end(&gc[i - 1])) < tol {
                    gc[j] = gc[j].reversed();
                    gc.swap(i, j);
                    trouve = true;
                }
                j += 1;
            }
        }

        // if still non matched -> error, the algorithm cannot finish
        if !trouve {
            return false;
        }
    }

    *cc1 = gc[0].clone();
    *cc2 = gc[1].clone();
    *cc3 = gc[2].reversed();
    *cc4 = gc[3].reversed();
    true
}

/// OCCT Geom_BSplineCurve::StartPoint.
fn ops_start(c: &BSplineCurve3) -> DVec3 {
    use rcad_kernel::geom::CurveEval;
    c.point_at(c.first_parameter())
}

/// OCCT Geom_BSplineCurve::EndPoint.
fn ops_end(c: &BSplineCurve3) -> DVec3 {
    use rcad_kernel::geom::CurveEval;
    c.point_at(c.last_parameter())
}

/// OCCT static SetSameDistribution (GeomFill_BSplineCurves.cxx L134-268) —
/// harmonizes the knot distributions of two boundary curves and returns the
/// common pole count.
fn set_same_distribution(c1: &mut BSplineCurve3, c2: &mut BSplineCurve3) -> usize {
    let (mut k1, m1) = c1.knots_mults();
    let p1 = c1.control_points.clone();
    let w1 = c1.weights.clone();
    let (mut k2, m2) = c2.knots_mults();
    let p2 = c2.control_points.clone();
    let w2 = c2.weights.clone();

    let k11 = k1[0];
    let k12 = k1[k1.len() - 1];
    let k21 = k2[0];
    let k22 = k2[k2.len() - 1];

    if (k12 - k11) > (k22 - k21) {
        bspl_reparametrize(k11, k12, &mut k2);
        c2.set_knots(&k2, &m2);
    } else if (k12 - k11) < (k22 - k21) {
        bspl_reparametrize(k21, k22, &mut k1);
        c1.set_knots(&k1, &m1);
    } else if (k12 - k11).abs() > PCONFUSION {
        bspl_reparametrize(k11, k12, &mut k2);
        c2.set_knots(&k2, &m2);
    }

    let mut np = 0i32;
    let mut nk = 0i32;
    if bspl_prepare_insert_knots(
        c1.degree,
        false,
        &k1,
        &m1,
        &k2,
        Some(&m2),
        &mut np,
        &mut nk,
        PCONFUSION,
        false,
    ) {
        // Homogeneous flattening (OCCT passes Poles + WeightsArray).
        let (poles_flat, dim) = flatten_poles(&p1, &w1);
        let mut new_p = vec![0.0f64; np as usize * dim];
        let mut new_k = vec![0.0f64; nk as usize];
        let mut new_m = vec![0i32; nk as usize];
        // AddKnots: OCCT inserts the knots of the other curve; the flat
        // insertion reads only the knot values from the homogeneous array,
        // so the flat knot slice of C2 is passed directly.
        bspl_insert_knots(
            c1.degree,
            false,
            dim,
            &poles_flat,
            &k1,
            &m1,
            &k2,
            Some(&m2),
            &mut new_p,
            &mut new_k,
            &mut new_m,
            PCONFUSION,
            false,
        );
        let (new_poles, new_weights) = unflatten_poles(&new_p, dim);
        *c1 = BSplineCurve3 {
            degree: c1.degree,
            knots: expand_knots(&new_k, &new_m),
            control_points: new_poles,
            weights: new_weights,
            is_periodic: false,
        };

        let (poles_flat, dim) = flatten_poles(&p2, &w2);
        let mut new_p = vec![0.0f64; np as usize * dim];
        let mut new_k = vec![0.0f64; nk as usize];
        let mut new_m = vec![0i32; nk as usize];
        bspl_insert_knots(
            c2.degree,
            false,
            dim,
            &poles_flat,
            &k2,
            &m2,
            &k1,
            Some(&m1),
            &mut new_p,
            &mut new_k,
            &mut new_m,
            PCONFUSION,
            false,
        );
        let (new_poles, new_weights) = unflatten_poles(&new_p, dim);
        *c2 = BSplineCurve3 {
            degree: c2.degree,
            knots: expand_knots(&new_k, &new_m),
            control_points: new_poles,
            weights: new_weights,
            is_periodic: false,
        };
    } else {
        panic!("Standard_ConstructionError: GeomFill_BSplineCurves");
    }

    c1.control_points.len()
}

/// Flatten poles + weights the BSplCLib homogeneous way (dim = 4) — OCCT
/// passes Poles and WeightsArray as separate arrays; the rcad flat kernel
/// interleaves them.
fn flatten_poles(poles: &[DVec3], weights: &[f64]) -> (Vec<f64>, usize) {
    let mut flat = Vec::with_capacity(poles.len() * 4);
    for (p, w) in poles.iter().zip(weights.iter()) {
        flat.extend([p.x, p.y, p.z, *w]);
    }
    (flat, 4)
}

/// Inverse of [`flatten_poles`].
fn unflatten_poles(flat: &[f64], dim: usize) -> (Vec<DVec3>, Vec<f64>) {
    let mut poles = Vec::with_capacity(flat.len() / dim);
    let mut weights = Vec::with_capacity(flat.len() / dim);
    for chunk in flat.chunks(dim) {
        poles.push(DVec3::new(chunk[0], chunk[1], chunk[2]));
        weights.push(if dim == 4 { chunk[3] } else { 1.0 });
    }
    (poles, weights)
}

/// Expand a (knots, mults) pair into the rcad flat knot vector.
fn expand_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut flat = Vec::with_capacity(mults.iter().map(|&m| m as usize).sum());
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            flat.push(*k);
        }
    }
    flat
}

/// OCCT GeomFill_BSplineCurves.
#[derive(Debug, Clone, Default)]
pub struct BSplineCurves {
    my_surface: Option<BSplineSurface>,
}

impl BSplineCurves {
    /// OCCT ctor (C1, C2, C3, C4, Type) (L138-146).
    pub fn new(
        c1: &BSplineCurve3,
        c2: &BSplineCurve3,
        c3: &BSplineCurve3,
        c4: &BSplineCurve3,
        filling_type: FillingStyle,
    ) -> Self {
        let mut curves = BSplineCurves { my_surface: None };
        curves.init(c1, c2, c3, c4, filling_type);
        curves
    }

    /// OCCT ctor (C1, C2, C3, Type) (L148-155).
    pub fn new3(
        c1: &BSplineCurve3,
        c2: &BSplineCurve3,
        c3: &BSplineCurve3,
        filling_type: FillingStyle,
    ) -> Self {
        let mut curves = BSplineCurves { my_surface: None };
        curves.init3(c1, c2, c3, filling_type);
        curves
    }

    /// OCCT ctor (C1, C2, Type) (L157-163).
    pub fn new2(c1: &BSplineCurve3, c2: &BSplineCurve3, filling_type: FillingStyle) -> Self {
        let mut curves = BSplineCurves { my_surface: None };
        curves.init2(c1, c2, filling_type);
        curves
    }

    /// OCCT Init(C1, C2, C3, C4, Type) (L165-390).
    pub fn init(
        &mut self,
        c1: &BSplineCurve3,
        c2: &BSplineCurve3,
        c3: &BSplineCurve3,
        c4: &BSplineCurve3,
        filling_type: FillingStyle,
    ) {
        // On ordonne les courbes
        let mut cc1 = c1.clone();
        let mut cc2 = c2.clone();
        let mut cc3 = c3.clone();
        let mut cc4 = c4.clone();
        let tol = CONFUSION;
        let is_ok = arrange(c1, c2, c3, c4, &mut cc1, &mut cc2, &mut cc3, &mut cc4, tol);
        assert!(
            is_ok,
            "GeomFill_BSplineCurves: Courbes non jointives"
        );

        // Mise en conformite des degres
        let deg1 = cc1.degree;
        let deg2 = cc2.degree;
        let deg3 = cc3.degree;
        let deg4 = cc4.degree;
        let degu = deg1.max(deg3);
        let degv = deg2.max(deg4);
        if deg1 < degu {
            cc1.increase_degree(degu);
        }
        if deg2 < degv {
            cc2.increase_degree(degv);
        }
        if deg3 < degu {
            cc3.increase_degree(degu);
        }
        if deg4 < degv {
            cc4.increase_degree(degv);
        }

        // Mise en conformite des distributions de noeuds
        let nbu_poles = set_same_distribution(&mut cc1, &mut cc3);
        let nbv_poles = set_same_distribution(&mut cc2, &mut cc4);
        if filling_type == FillingStyle::CoonsStyle && (nbu_poles < 4 || nbv_poles < 4) {
            panic!("GeomFill_BSplineCurves: invalid filling style");
        }

        let p1 = cc1.control_points.clone();
        let p2 = cc2.control_points.clone();
        let p3 = cc3.control_points.clone();
        let p4 = cc4.control_points.clone();

        // Traitement des courbes rationelles
        let is_rat = cc1.is_rational() || cc2.is_rational() || cc3.is_rational() || cc4.is_rational();
        let w1 = cc1.weights.clone();
        let w2 = cc2.weights.clone();
        let w3 = cc3.weights.clone();
        let w4 = cc4.weights.clone();

        // GeomFill_Filling Caro — the OCCT base-class handle; note the Coons
        // argument permutation (P1, P4, P3, P2 / W1, W4, W3, W2).
        let caro: FillingVariant = match filling_type {
            FillingStyle::StretchStyle => {
                if is_rat {
                    FillingVariant::Stretch(super::stretch::Stretch::new_rational(
                        &p1, &p2, &p3, &p4, &w1, &w2, &w3, &w4,
                    ))
                } else {
                    FillingVariant::Stretch(super::stretch::Stretch::new(&p1, &p2, &p3, &p4))
                }
            }
            FillingStyle::CoonsStyle => {
                if is_rat {
                    FillingVariant::Coons(super::coons::Coons::new_rational(
                        &p1, &p4, &p3, &p2, &w1, &w4, &w3, &w2,
                    ))
                } else {
                    FillingVariant::Coons(super::coons::Coons::new(&p1, &p4, &p3, &p2))
                }
            }
            FillingStyle::CurvedStyle => {
                if is_rat {
                    FillingVariant::Curved(super::curved::Curved::new_rational(
                        &p1, &p2, &p3, &p4, &w1, &w2, &w3, &w4,
                    ))
                } else {
                    FillingVariant::Curved(super::curved::Curved::new(&p1, &p2, &p3, &p4))
                }
            }
        };

        let nbu_poles = caro.nb_u_poles();
        let nbv_poles = caro.nb_v_poles();

        // Creation de la surface
        let (u_knots, u_mults) = cc1.knots_mults();
        let (v_knots, v_mults) = cc2.knots_mults();

        let weights = if caro.is_rational() {
            caro.weights().clone()
        } else {
            vec![vec![1.0f64; nbv_poles]; nbu_poles]
        };

        self.my_surface = Some(BSplineSurface {
            degree_u: cc1.degree,
            degree_v: cc2.degree,
            knots_u: expand_knots(&u_knots, &u_mults),
            knots_v: expand_knots(&v_knots, &v_mults),
            control_points: caro.poles().clone(),
            weights,
        });
    }

    /// OCCT Init(C1, C2, C3, Type) (L392-437) — builds a degenerate fourth
    /// curve and delegates to the 4-curve Init.
    pub fn init3(
        &mut self,
        c1: &BSplineCurve3,
        c2: &BSplineCurve3,
        c3: &BSplineCurve3,
        filling_type: FillingStyle,
    ) {
        let mut poles = [DVec3::ZERO; 2];
        // Knots(1, 2); Mults(1, 2);
        let mut tol = CONFUSION;
        tol *= tol;
        if ops_start(c1).distance(ops_start(c2)) > tol
            && ops_start(c1).distance(ops_end(c2)) > tol
        {
            poles[0] = ops_start(c1);
        } else {
            poles[0] = ops_end(c1);
        }
        if ops_start(c3).distance(ops_start(c2)) > tol
            && ops_start(c3).distance(ops_end(c2)) > tol
        {
            poles[1] = ops_start(c3);
        } else {
            poles[1] = ops_end(c3);
        }
        // Knots(1) = C2->Knot(C2->FirstUKnotIndex());
        // Knots(2) = C2->Knot(C2->LastUKnotIndex());
        let knots = vec![c2.knots[0], c2.knots[c2.knots.len() - 1]];
        let mults = vec![2, 2];
        let c4 = BSplineCurve3::from_knots_mults(1, knots, mults, poles.to_vec());
        self.init(c1, c2, c3, &c4, filling_type);
    }

    /// OCCT Init(C1, C2, Type) (L439-597).
    pub fn init2(
        &mut self,
        c1: &BSplineCurve3,
        c2: &BSplineCurve3,
        filling_type: FillingStyle,
    ) {
        let mut cc1 = c1.clone();
        let mut cc2 = c2.clone();
        // OCCT keeps a mutable handle to the INPUT C1 (see below, L474).
        let mut c1_mut = c1.clone();

        let deg1 = cc1.degree;
        let deg2 = cc2.degree;
        let is_rat = cc1.is_rational() || cc2.is_rational();

        if filling_type != FillingStyle::CurvedStyle {
            let degu = deg1.max(deg2);
            if cc1.degree < degu {
                cc1.increase_degree(degu);
            }
            if cc2.degree < degu {
                cc2.increase_degree(degu);
            }

            // Mise en conformite des distributions de noeuds
            let nb_poles = set_same_distribution(&mut cc1, &mut cc2);
            let p1 = cc1.control_points.clone();
            let p2 = cc2.control_points.clone();
            let mut poles = vec![vec![DVec3::ZERO; 2]; nb_poles];
            for i in 1..=nb_poles {
                poles[i - 1][0] = p1[i - 1];
                poles[i - 1][1] = p2[i - 1];
            }
            let (u_knots, u_mults) = cc1.knots_mults();
            let v_knots = [0.0f64, 1.0];
            let v_mults = [2, 2];

            let weights = if is_rat {
                let w1 = cc1.weights.clone();
                let w2 = cc2.weights.clone();
                let mut weights = vec![vec![0.0f64; 2]; nb_poles];
                for i in 1..=nb_poles {
                    weights[i - 1][0] = w1[i - 1];
                    weights[i - 1][1] = w2[i - 1];
                }
                weights
            } else {
                vec![vec![1.0f64; 2]; nb_poles]
            };

            self.my_surface = Some(BSplineSurface {
                degree_u: cc1.degree,
                degree_v: 1,
                knots_u: expand_knots(&u_knots, &u_mults),
                knots_v: expand_knots(&v_knots, &v_mults),
                control_points: poles,
                weights,
            });
        } else {
            let eps = CONFUSION;
            let mut is_ok = false;
            if ops_start(&cc1).distance(ops_start(&cc2)) <= eps {
                is_ok = true;
            } else if ops_start(&cc1).distance(ops_end(&cc2)) <= eps {
                cc2 = cc2.reversed();
                is_ok = true;
            } else if ops_end(&cc1).distance(ops_start(&cc2)) <= eps {
                // OCCT quirk (L474): the INPUT curve C1 is reversed here,
                // not the local copy CC1 — the fill below keeps using the
                // unreversed CC1 poles.
                c1_mut = c1_mut.reversed();
                is_ok = true;
            } else if ops_end(&cc1).distance(ops_end(&cc2)) <= eps {
                cc1 = cc1.reversed();
                cc2 = cc2.reversed();
                is_ok = true;
            }
            let _ = c1_mut; // reversal is caller-visible only in OCCT
            if !is_ok {
                panic!("GeomFill_BSplineCurves: Courbes non jointives");
            }

            let p1 = cc1.control_points.clone();
            let p2 = cc2.control_points.clone();
            let w1 = cc1.weights.clone();
            let w2 = cc2.weights.clone();

            let (u_knots, u_mults) = cc1.knots_mults();
            let (v_knots, v_mults) = cc2.knots_mults();

            let caro: FillingVariant = if is_rat {
                FillingVariant::Curved(super::curved::Curved::new_two_rational(&p1, &p2, &w1, &w2))
            } else {
                FillingVariant::Curved(super::curved::Curved::new_two(&p1, &p2))
            };

            let nbu_poles = caro.nb_u_poles();
            let nbv_poles = caro.nb_v_poles();
            let weights = if caro.is_rational() {
                caro.weights().clone()
            } else {
                vec![vec![1.0f64; nbv_poles]; nbu_poles]
            };

            self.my_surface = Some(BSplineSurface {
                degree_u: cc1.degree,
                degree_v: cc2.degree,
                knots_u: expand_knots(&u_knots, &u_mults),
                knots_v: expand_knots(&v_knots, &v_mults),
                control_points: caro.poles().clone(),
                weights,
            });
        }
    }

    /// OCCT Surface() — the filling result (null before a successful Init).
    pub fn surface(&self) -> Option<&BSplineSurface> {
        self.my_surface.as_ref()
    }
}

/// OCCT `GeomFill_Filling Caro;` — the base-class handle holding one of the
/// concrete fillings (GeomFill_Filling.hxx L33-53).
enum FillingVariant {
    Stretch(super::stretch::Stretch),
    Coons(super::coons::Coons),
    Curved(super::curved::Curved),
}

impl FillingVariant {
    fn nb_u_poles(&self) -> usize {
        match self {
            FillingVariant::Stretch(f) => f.base.nb_u_poles(),
            FillingVariant::Coons(f) => f.base.nb_u_poles(),
            FillingVariant::Curved(f) => f.base.nb_u_poles(),
        }
    }

    fn nb_v_poles(&self) -> usize {
        match self {
            FillingVariant::Stretch(f) => f.base.nb_v_poles(),
            FillingVariant::Coons(f) => f.base.nb_v_poles(),
            FillingVariant::Curved(f) => f.base.nb_v_poles(),
        }
    }

    fn is_rational(&self) -> bool {
        match self {
            FillingVariant::Stretch(f) => f.base.is_rational(),
            FillingVariant::Coons(f) => f.base.is_rational(),
            FillingVariant::Curved(f) => f.base.is_rational(),
        }
    }

    fn poles(&self) -> &Vec<Vec<DVec3>> {
        match self {
            FillingVariant::Stretch(f) => f.base.poles(),
            FillingVariant::Coons(f) => f.base.poles(),
            FillingVariant::Curved(f) => f.base.poles(),
        }
    }

    fn weights(&self) -> &Vec<Vec<f64>> {
        match self {
            FillingVariant::Stretch(f) => f.base.weights(),
            FillingVariant::Coons(f) => f.base.weights(),
            FillingVariant::Curved(f) => f.base.weights(),
        }
    }
}
