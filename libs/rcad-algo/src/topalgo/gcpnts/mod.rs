// OCCT GCPnts_QuasiUniformDeflection (GCPnts_QuasiUniformDeflection.cxx)
// Discretize a 2D curve with bounded chordal deflection.
//
// Used by IntTools_FClass2d::Init's re-discretization loop (the "while"
// loop that tightens the polygon when its deflection exceeds the expected
// thickness for its area/perimeter ratio).

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::CONFUSION;

// OCCT GCPnts_QuasiUniformDeflection.cxx L23: MyMaxQuasiFleshe
const MY_MAX_QUASI_FLESHE: i32 = 2000;

// OCCT GCPnts_DeflectionType (GCPnts_DeflectionType.hxx)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DeflectionType {
    Linear,
    Circular,
    Curved,
    DefComposite,
}

// OCCT GCPnts_QuasiUniformDeflection.cxx L322-348: GetDefType
fn get_def_type(the_c: &Curve2d) -> DeflectionType {
    match the_c {
        Curve2d::Line(_) => DeflectionType::Linear,
        Curve2d::Circle(_) => DeflectionType::Circular,
        Curve2d::BSpline(a_bs) => {
            if a_bs.control_points.len() == 2 {
                DeflectionType::Linear
            } else {
                DeflectionType::Curved
            }
        }
        Curve2d::Bezier(a_bz) => {
            if a_bz.control_points.len() == 2 {
                DeflectionType::Linear
            } else {
                DeflectionType::Curved
            }
        }
        Curve2d::Trimmed(a_t) => get_def_type(&a_t.curve),
        // DefComposite (multiple C1 intervals) not distinguished — the
        // recursive QuasiFleche path handles all remaining types.
        _ => DeflectionType::Curved,
    }
}

fn circle_radius(the_c: &Curve2d) -> f64 {
    match the_c {
        Curve2d::Circle(a_c) => a_c.radius,
        Curve2d::Trimmed(a_t) => circle_radius(&a_t.curve),
        _ => 1.0,
    }
}

// OCCT GCPnts_QuasiUniformDeflection.cxx L277-291: PerformLinear
fn perform_linear(
    the_c: &Curve2d,
    the_parameters: &mut Vec<f64>,
    the_points: &mut Vec<DVec2>,
    the_u1: f64,
    the_u2: f64,
) -> bool {
    the_parameters.push(the_u1);
    the_points.push(the_c.point_at(the_u1));
    the_parameters.push(the_u2);
    the_points.push(the_c.point_at(the_u2));
    true
}

// OCCT GCPnts_QuasiUniformDeflection.cxx L295-317: PerformCircular
fn perform_circular(
    the_c: &Curve2d,
    the_parameters: &mut Vec<f64>,
    the_points: &mut Vec<DVec2>,
    the_deflection: f64,
    the_u1: f64,
    the_u2: f64,
) -> bool {
    let mut an_angle = (1.0 - (the_deflection / circle_radius(the_c))).max(0.0);
    an_angle = 2.0 * an_angle.acos();
    let mut a_nb_points = ((the_u2 - the_u1) / an_angle) as i32;
    a_nb_points += 2;
    an_angle = (the_u2 - the_u1) / (a_nb_points - 1) as f64;
    let mut u = the_u1;
    for _i in 1..=a_nb_points {
        the_parameters.push(u);
        the_points.push(the_c.point_at(u));
        u += an_angle;
    }
    true
}

// OCCT GCPnts_QuasiUniformDeflection.cxx L56-186: QuasiFleche (D1 overload)
struct FlecheCtx<'a> {
    the_c: &'a Curve2d,
    the_deflection2: f64,
    the_eps: f64,
    the_parameters: &'a mut Vec<f64>,
    the_points: &'a mut Vec<DVec2>,
    the_nb_calls: i32,
}

impl FlecheCtx<'_> {
    fn value(&self, u: f64) -> DVec2 {
        self.the_c.point_at(u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        (self.the_c.point_at(u), self.the_c.derivative_at(u))
    }

    fn quasi_fleche(
        &mut self,
        the_u_deb: f64,
        the_p_deb: DVec2,
        the_v_deb: DVec2,
        the_u_fin: f64,
        the_p_fin: DVec2,
        the_v_fin: DVec2,
        the_nbmin: i32,
    ) {
        self.the_nb_calls += 1;
        if self.the_nb_calls >= MY_MAX_QUASI_FLESHE {
            return;
        }
        let a_ptslength = self.the_points.len();
        if self.the_nb_calls > 100 && a_ptslength < 2 {
            return;
        }
        let mut a_udelta = the_u_fin - the_u_deb;
        let mut a_pdelta: DVec2;
        let mut a_vdelta: DVec2;
        if the_nbmin > 2 {
            a_udelta /= (the_nbmin - 1) as f64;
            let (p, v) = self.d1(the_u_deb + a_udelta);
            a_pdelta = p;
            a_vdelta = v;
        } else {
            a_pdelta = the_p_fin;
            a_vdelta = the_v_fin;
        }
        // Square length of chord.
        let a_norme = (the_p_deb - a_pdelta).length_squared();
        let mut a_fleche = 0.0;
        let mut is_fleche_ok = false;
        if a_norme > self.the_eps && a_norme > 16.0 * self.the_deflection2 {
            // Evaluation of the deflection by interpolation.
            // See IntWalk_IWalking::TestDeflection.
            let n1 = the_v_deb.length_squared();
            let n2 = a_vdelta.length_squared();
            if n1 > self.the_eps && n2 > self.the_eps {
                // Square distance between ends of two normalized vectors [0; 4].
                let a_norme_diff = (the_v_deb.normalize() - a_vdelta.normalize()).length_squared();
                if a_norme_diff > self.the_eps {
                    a_fleche = a_norme_diff * a_norme / 64.0;
                    is_fleche_ok = true;
                }
            }
        }
        let a_p_mid = (the_p_deb + a_pdelta) * 0.5;
        let a_p_verif = self.value(the_u_deb + a_udelta * 0.5);
        let a_fleche_mid_mid = (a_p_mid - a_p_verif).length_squared();
        if is_fleche_ok {
            // The interpolation-based algorithm can give a false-positive
            // result, so check also the Pmid-Pverif distance. But
            // aFlecheMidMid gives a worse result for non-uniform
            // parameterisation.
            if a_fleche_mid_mid > a_norme / 4.0 + self.the_deflection2 {
                a_fleche = a_fleche_mid_mid;
            }
        } else {
            a_fleche = a_fleche_mid_mid;
        }
        if a_fleche < self.the_deflection2 {
            self.the_parameters.push(the_u_deb + a_udelta);
            self.the_points.push(a_pdelta);
        } else {
            self.quasi_fleche(
                the_u_deb,
                the_p_deb,
                the_v_deb,
                the_u_deb + a_udelta,
                a_pdelta,
                a_vdelta,
                3,
            );
        }
        if the_nbmin > 2 {
            self.quasi_fleche(
                the_u_deb + a_udelta,
                a_pdelta,
                a_vdelta,
                the_u_fin,
                the_p_fin,
                the_v_fin,
                the_nbmin - (self.the_points.len() - a_ptslength) as i32,
            );
        }
        self.the_nb_calls -= 1;
    }
}

// OCCT GCPnts_QuasiUniformDeflection.cxx L352-410: PerformCurve.
// The FClass2d usage constructs with default continuity GeomAbs_C1, so
// myCont = C1 and the D1-based branch (derivatives) is taken.
fn perform_curve(
    the_c: &Curve2d,
    the_parameters: &mut Vec<f64>,
    the_points: &mut Vec<DVec2>,
    the_deflection: f64,
    the_u1: f64,
    the_u2: f64,
    the_epsilon: f64,
) -> bool {
    let a_nbmin = 2i32;
    let (a_p_deb, a_d_deb) = (the_c.point_at(the_u1), the_c.derivative_at(the_u1));
    the_parameters.push(the_u1);
    the_points.push(a_p_deb);
    // OCCT: const double aDecreasedU2 = theU2 - Epsilon(theU2) * 10.0;
    let a_decreased_u2 = the_u2 - f64::EPSILON * the_u2.abs() * 10.0;
    let (a_p_fin, a_d_fin) = (the_c.point_at(a_decreased_u2), the_c.derivative_at(a_decreased_u2));
    let mut a_ctx = FlecheCtx {
        the_c,
        the_deflection2: the_deflection * the_deflection,
        the_eps: the_epsilon * the_epsilon,
        the_parameters,
        the_points,
        the_nb_calls: 0,
    };
    a_ctx.quasi_fleche(the_u1, a_p_deb, a_d_deb, the_u2, a_p_fin, a_d_fin, a_nbmin);
    true
}

/// OCCT GCPnts_QuasiUniformDeflection (2D) — discretize `the_c` on
/// `[the_u1, the_u2]` so the chordal deflection stays below `the_deflection`.
pub struct QuasiUniformDeflection {
    params: Vec<f64>,
    points: Vec<DVec2>,
    done: bool,
}

impl QuasiUniformDeflection {
    /// OCCT constructor with explicit parameter range.
    pub fn new(the_c: &Curve2d, the_deflection: f64, the_u1: f64, the_u2: f64) -> Self {
        let mut q = QuasiUniformDeflection {
            params: Vec::new(),
            points: Vec::new(),
            done: false,
        };
        q.initialize(the_c, the_deflection, the_u1, the_u2);
        q
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// OCCT aDiscr.Parameter(i) — 1-based index.
    pub fn parameter(&self, the_index: usize) -> f64 {
        self.params[the_index - 1]
    }

    // OCCT GCPnts_QuasiUniformDeflection.cxx L573-623: initialize
    fn initialize(&mut self, the_c: &Curve2d, the_deflection: f64, the_u1: f64, the_u2: f64) {
        // myCont = C1 (the FClass2d caller uses the default GeomAbs_C1).
        self.done = false;
        self.params.clear();
        self.points.clear();
        // OCCT: anEPSILON = min(theC.Resolution(Confusion()), 1e50).
        // Resolution(Confusion()) is the parameter interval for a 3D move of
        // CONFUSION; approximated here (a guard for the recursion).
        let an_epsilon = CONFUSION;
        let a_type = get_def_type(the_c);
        let a_u1 = the_u1.min(the_u2);
        let a_u2 = the_u1.max(the_u2);
        if a_type == DeflectionType::Curved || a_type == DeflectionType::DefComposite {
            if matches!(the_c, Curve2d::BSpline(_) | Curve2d::Bezier(_)) {
                let a_max_par = a_u1.abs().max(a_u2.abs());
                if an_epsilon < f64::EPSILON * a_max_par.abs() {
                    return;
                }
            }
        }
        self.done = match a_type {
            DeflectionType::Linear => {
                perform_linear(the_c, &mut self.params, &mut self.points, a_u1, a_u2)
            }
            DeflectionType::Circular => perform_circular(
                the_c,
                &mut self.params,
                &mut self.points,
                the_deflection,
                a_u1,
                a_u2,
            ),
            DeflectionType::Curved | DeflectionType::DefComposite => perform_curve(
                the_c,
                &mut self.params,
                &mut self.points,
                the_deflection,
                a_u1,
                a_u2,
                an_epsilon,
            ),
        };
    }
}
