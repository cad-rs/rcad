//! OCCT Extrema_GenExtCS (TKGeomBase/Extrema/Extrema_GenExtCS.cxx) with its
//! math support classes math_PSO / math_PSOParticlesPool /
//! math_BullardGenerator (FoundationClasses/TKMath/math), the objective
//! functions Extrema_GlobOptFuncCS / Extrema_GlobOptFuncConicS and the
//! solution function set Extrema_FuncExtCS (TKGeomBase/Extrema).
//!
//! The final refinement of Perform is the N-variable math_FunctionSetRoot
//! (translated 1:1 in `rcad-kernel/src/math/function_set_root.rs`) applied to
//! the 3-variable Extrema_FuncExtCS function set.

use crate::bop::int_tools::bean_face_intersector::BRepAdaptorCurve;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};

const MAX_PARAM_VAL: f64 = 1.0e+10;
const BORDER_DIVISOR: f64 = 1.0e+4;
// OCCT HyperbolaLimit = 23. (ln(MaxParamVal))
const HYPERBOLA_LIMIT: f64 = 23.0;

// OCCT math_BullardGenerator (math_BullardGenerator.hxx) — Ian C. Bullard's
// fast RNG, two interleaved 32-bit state words.
pub(crate) struct BullardGenerator {
    state_hi: u32,
    state_lo: u32,
}

impl BullardGenerator {
    pub fn new() -> Self {
        BullardGenerator {
            state_hi: 0x9E3779B9,
            state_lo: 0x12345678,
        }
    }

    fn next_int(&mut self) -> u32 {
        self.state_hi = self.state_hi.wrapping_shr(2).wrapping_add(self.state_hi.wrapping_shl(2));
        self.state_hi = self.state_hi.wrapping_add(self.state_lo);
        self.state_lo = self.state_lo.wrapping_add(self.state_hi);
        self.state_hi
    }

    pub fn next_real(&mut self) -> f64 {
        self.next_int() as f64 / 0xFFFF_FFFFu32 as f64
    }
}

// OCCT PSO_Particle / math_PSOParticlesPool (math_PSOParticlesPool.hxx).
pub(crate) struct PsoParticle {
    pub position: Vec<f64>,
    pub best_position: Vec<f64>,
    pub velocity: Vec<f64>,
    pub distance: f64,
    pub best_distance: f64,
}

pub(crate) struct PsoParticlesPool {
    particles: Vec<PsoParticle>,
}

impl PsoParticlesPool {
    pub fn new(nb_particles: usize, nb_vars: usize) -> Self {
        let mut particles = Vec::with_capacity(nb_particles);
        for _ in 0..nb_particles {
            particles.push(PsoParticle {
                position: vec![0.0; nb_vars],
                best_position: vec![0.0; nb_vars],
                velocity: vec![0.0; nb_vars],
                distance: f64::MAX,
                best_distance: f64::MAX,
            });
        }
        PsoParticlesPool { particles }
    }

    /// OCCT GetWorstParticle — the particle with the greatest Distance.
    pub fn worst(&mut self) -> &mut PsoParticle {
        let mut idx = 0;
        let mut worst = f64::MIN;
        for (i, p) in self.particles.iter().enumerate() {
            if p.distance > worst {
                worst = p.distance;
                idx = i;
            }
        }
        &mut self.particles[idx]
    }

    pub fn best(&self) -> &PsoParticle {
        let mut idx = 0;
        let mut best = f64::MAX;
        for (i, p) in self.particles.iter().enumerate() {
            if p.distance < best {
                best = p.distance;
                idx = i;
            }
        }
        &self.particles[idx]
    }
}

/// OCCT math_PSO::performPSOWithGivenParticles (math_PSO.cxx L125-268) for a
/// func returning Option<f64> (invalid points yield Precision::Infinite).
pub(crate) fn pso_with_given_particles(
    func: &mut dyn FnMut(&[f64]) -> f64,
    low_border: &[f64],
    upp_border: &[f64],
    steps: &[f64],
    particles: &mut PsoParticlesPool,
    nb_iter: i32,
) -> (f64, Vec<f64>) {
    let n = low_border.len();
    let min_uv: Vec<f64> = (0..n)
        .map(|i| low_border[i] + (upp_border[i] - low_border[i]) / BORDER_DIVISOR)
        .collect();
    let max_uv: Vec<f64> = (0..n)
        .map(|i| upp_border[i] - (upp_border[i] - low_border[i]) / BORDER_DIVISOR)
        .collect();

    let mut rng = BullardGenerator::new();
    for p in particles.particles.iter_mut() {
        for d in 0..n {
            let ksi = rng.next_real();
            p.velocity[d] = steps[d] * (ksi - 0.5) * 2.0;
        }
    }

    let best = particles.best();
    let mut best_global_position = best.position.clone();
    let mut best_global_distance = best.distance;

    let termination_velocity: Vec<f64> = steps.iter().map(|s| s / 2048.0).collect();

    let mut step = 1i32;
    while step < nb_iter {
        let mut minimal_velocity = vec![f64::MAX; n];

        for p in particles.particles.iter_mut() {
            let ksi1 = rng.next_real();
            let ksi2 = rng.next_real();

            const RETENT_WEIGHT: f64 = 0.72900;
            const PERSON_WEIGHT: f64 = 1.49445;
            const SOCIAL_WEIGHT: f64 = 1.49445;

            for d in 0..n {
                p.velocity[d] = p.velocity[d] * RETENT_WEIGHT
                    + (p.best_position[d] - p.position[d]) * (PERSON_WEIGHT * ksi1)
                    + (best_global_position[d] - p.position[d]) * (SOCIAL_WEIGHT * ksi2);

                p.position[d] += p.velocity[d];
                p.position[d] = p.position[d].clamp(min_uv[d], max_uv[d]);
                minimal_velocity[d] = minimal_velocity[d].min(p.velocity[d].abs());
            }

            p.distance = func(&p.position);
            if p.distance < p.best_distance {
                p.best_distance = p.distance;
                p.best_position.copy_from_slice(&p.position);
                if p.distance < best_global_distance {
                    best_global_distance = p.distance;
                    best_global_position.copy_from_slice(&p.position);
                }
            }
        }

        let is_terminal_velocity_reached = minimal_velocity
            .iter()
            .zip(termination_velocity.iter())
            .all(|(v, t)| *v <= *t);
        if is_terminal_velocity_reached {
            // OCCT: minimum number of steps
            const MIN_STEPS: i32 = 16;
            if step > MIN_STEPS {
                break;
            }
            for p in particles.particles.iter_mut() {
                let ksi = rng.next_real();
                for d in 0..n {
                    if p.position[d] == min_uv[d] || p.position[d] == max_uv[d] {
                        p.velocity[d] = if p.position[d] == min_uv[d] {
                            steps[d] * ksi
                        } else {
                            -steps[d] * ksi
                        };
                    } else {
                        p.velocity[d] = steps[d] * (ksi - 0.5) * 2.0;
                    }
                }
            }
        }
        step += 1;
    }
    (best_global_distance, best_global_position)
}

/// OCCT Extrema_GlobOptFuncCS — F(cu,su,sv) = |C(cu) - S(su,sv)|^2 with the
/// analytic gradient (GlobOptFuncCS.cxx L30-56).
struct GlobOptFuncCS<'a> {
    curve: &'a Curve3,
    surface: &'a Surface3,
    tf: f64,
    tl: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
}

impl<'a> GlobOptFuncCS<'a> {
    fn check_input_data(&self, cu: f64, su: f64, sv: f64) -> bool {
        cu >= self.tf && cu <= self.tl && su >= self.umin && su <= self.umax && sv >= self.vmin && sv <= self.vmax
    }

    fn value(&self, cu: f64, su: f64, sv: f64) -> f64 {
        self.curve.point_at(cu).distance_squared(self.surface.point_at(su, sv))
    }

    fn gradient(&self, cu: f64, su: f64, sv: f64) -> [f64; 3] {
        let cd0 = self.curve.point_at(cu);
        let cd1 = self.curve.derivative_at(cu);
        let sd0 = self.surface.point_at(su, sv);
        let (sd1u, sd1v) = surface_d1(self.surface, su, sv);
        let diff = cd0 - sd0;
        [diff.dot(cd1), -diff.dot(sd1u), -diff.dot(sd1v)]
    }
}

/// Surface D1 (SurfaceEval has no combined D1 in rcad; derivatives via the
/// point/derivative primitives — the Revolution basis uses finite differences
/// in normal_at, so mirror Geom_SurfaceOfRevolution D1 by central differences
/// with the OCCT adaptors' step scaled by the parameter ranges).
fn surface_d1(s: &Surface3, u: f64, v: f64) -> (DVec3, DVec3) {
    let domain = s.default_domain();
    let eps_u = ((domain[1] - domain[0]).max(1e-3)) * 1e-7;
    let eps_v = ((domain[3] - domain[2]).max(1e-3)) * 1e-7;
    let du = (s.point_at(u + eps_u, v) - s.point_at(u - eps_u, v)) / (2.0 * eps_u);
    let dv = (s.point_at(u, v + eps_v) - s.point_at(u, v - eps_v)) / (2.0 * eps_v);
    (du, dv)
}

/// OCCT Extrema_GlobOptFuncConicS — F(su,sv) = |S(su,sv) - C(ct)|^2 with
/// ct = ElCLib::Parameter(conic, S) (GlobOptFuncConicS.cxx L23-62).
struct GlobOptFuncConicS<'a> {
    surface: &'a Surface3,
    ctype: ConicType,
    line: Option<rcad_kernel::geom::Line3>,
    circle: Option<rcad_kernel::geom::Circle3>,
    ellipse: Option<rcad_kernel::geom::Ellipse3>,
    hyperbola: Option<rcad_kernel::geom::Hyperbola3>,
    parabola: Option<rcad_kernel::geom::Parabola3>,
    tf: f64,
    tl: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    _marker: std::marker::PhantomData<&'a ()>,
}

#[derive(Clone, Copy, PartialEq)]
enum ConicType {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    None,
}

impl<'a> GlobOptFuncConicS<'a> {
    fn conic_parameter(&self, su: f64, sv: f64) -> f64 {
        let ps = self.surface.point_at(su, sv);
        match self.ctype {
            ConicType::Line => {
                let l = self.line.as_ref().unwrap();
                crate::geomalgo::int_patch::elclib::line_parameter(l, ps)
            }
            ConicType::Circle => {
                let c = self.circle.as_ref().unwrap();
                crate::geomalgo::int_patch::elclib::circle_parameter(c, ps)
            }
            ConicType::Ellipse => {
                let e = self.ellipse.as_ref().unwrap();
                crate::geomalgo::int_patch::elclib::ellipse_parameter(e, ps)
            }
            ConicType::Hyperbola => {
                let h = self.hyperbola.as_ref().unwrap();
                crate::geomalgo::int_patch::elclib::hyperbola_parameter(h, ps)
            }
            ConicType::Parabola => {
                let p = self.parabola.as_ref().unwrap();
                crate::geomalgo::int_patch::elclib::parabola_parameter(p, ps)
            }
            ConicType::None => self.tf,
        }
    }

    /// OCCT value(su, sv, F) with the circle/ellipse period adjustment
    /// (GlobOptFuncConicS.cxx L57-90): a parameter outside [tf, tl] on a
    /// periodic conic is folded back before evaluating the distance.
    fn value(&self, su: f64, sv: f64) -> f64 {
        if su < self.umin || su > self.umax || sv < self.vmin || sv > self.vmax {
            return f64::INFINITY;
        }
        let ps = self.surface.point_at(su, sv);
        let mut ct = self.conic_parameter(su, sv);
        let pc = match self.ctype {
            ConicType::Line => {
                let l = self.line.as_ref().unwrap();
                l.point_at(ct)
            }
            ConicType::Circle => {
                let c = self.circle.as_ref().unwrap();
                if self.tl > 2.0 * std::f64::consts::PI + rcad_kernel::core::precision::PCONFUSION {
                    ct += 2.0 * std::f64::consts::PI;
                }
                c.point_at(ct)
            }
            ConicType::Ellipse => {
                let e = self.ellipse.as_ref().unwrap();
                if self.tl > 2.0 * std::f64::consts::PI + rcad_kernel::core::precision::PCONFUSION {
                    ct += 2.0 * std::f64::consts::PI;
                }
                e.point_at(ct)
            }
            ConicType::Hyperbola => self.hyperbola.as_ref().unwrap().point_at(ct),
            ConicType::Parabola => self.parabola.as_ref().unwrap().point_at(ct),
            ConicType::None => return f64::INFINITY,
        };
        ps.distance_squared(pc)
    }
}

/// OCCT Extrema_GenExtCS — the global curve-surface extremum search
/// (Extrema_GenExtCS.cxx): PSO seed over the (u,v) (or (t,u,v)) grid followed
/// by the math_FunctionSetRoot refinement on the GlobOptFuncCS system.
pub struct ExtremaGenExtCS {
    usample: usize,
    vsample: usize,
    tsample: usize,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    tol2: f64,
    surface: Option<Surface3>,
    // Curve type of the last Perform (conic dispatch of GlobMin).
    is_done: bool,
    // (t, u, v, sq_dist)
    extrema: Vec<(f64, f64, f64, f64)>,
}

impl Default for ExtremaGenExtCS {
    fn default() -> Self {
        Self::new()
    }
}

impl ExtremaGenExtCS {
    pub fn new() -> Self {
        ExtremaGenExtCS {
            usample: 0,
            vsample: 0,
            tsample: 0,
            umin: 0.0,
            umax: 0.0,
            vmin: 0.0,
            vmax: 0.0,
            tol2: 0.0,
            surface: None,
            is_done: false,
            extrema: Vec::new(),
        }
    }

    /// OCCT Initialize(S, NbU, NbV, Umin, Usup, Vmin, Vsup, Tol2)
    /// (Extrema_GenExtCS.cxx L208-286).
    pub fn initialize(
        &mut self,
        surface: &Surface3,
        nb_u: usize,
        nb_v: usize,
        u_min: f64,
        u_max: f64,
        v_min: f64,
        v_max: f64,
        tol2: f64,
    ) {
        self.surface = Some(surface.clone());
        self.usample = nb_u;
        self.vsample = nb_v;
        self.umin = u_min;
        self.umax = u_max;
        self.vmin = v_min;
        self.vmax = v_max;
        self.tol2 = tol2;

        // OCCT L226-244: infinite bounds are restricted to MaxParamVal (and to
        // HyperbolaLimit for hyperbola-based surfaces); the tests use
        // Precision::IsInfinite (Precision.hxx L350-353).
        let vmaxpar = match surface {
            Surface3::Revolution(_) => HYPERBOLA_LIMIT,
            _ => MAX_PARAM_VAL,
        };
        if rcad_kernel::precision::is_infinite_value(self.umax) {
            self.umax = MAX_PARAM_VAL;
        }
        if rcad_kernel::precision::is_infinite_value(self.umin) {
            self.umin = -MAX_PARAM_VAL;
        }
        if rcad_kernel::precision::is_infinite_value(self.vmax) {
            self.vmax = vmaxpar;
        }
        if rcad_kernel::precision::is_infinite_value(self.vmin) {
            self.vmin = -vmaxpar;
        }
    }

    /// OCCT Perform(C, NbT, tmin, tsup, Tol1) (L299-423).
    pub fn perform(
        &mut self,
        curve: &BRepAdaptorCurve,
        nb_t: usize,
        tmin: f64,
        tsup: f64,
        tol1: f64,
    ) {
        self.is_done = false;
        self.extrema.clear();
        self.tsample = nb_t;
        let surf = match self.surface.clone() {
            Some(s) => s,
            None => return,
        };

        // OCCT L305: myF.Initialize(C, *myS).
        let mut my_f = FuncExtCS::new(curve.curve(), &surf);

        // OCCT L312-320: infinite t bounds are trimmed to GetCurvMaxParamVal(C).
        let mut my_tmin = tmin;
        let mut my_tsup = tsup;
        let a_c_max_val = get_curv_max_param_val(curve);
        if rcad_kernel::precision::is_infinite_value(my_tsup) {
            my_tsup = a_c_max_val;
        }
        if rcad_kernel::precision::is_infinite_value(my_tmin) {
            my_tmin = -a_c_max_val;
        }

        // OCCT L324-337: aNbVar = 1 for quadric surfaces, 2 for conic curves,
        // 3 otherwise.
        let ctype = curve_type_of(curve);
        let stype = surface_type_of(&surf);
        let is_quadric = matches!(
            stype,
            SType::Plane | SType::Cylinder | SType::Cone | SType::Sphere | SType::Torus
        );
        let is_conic = !matches!(ctype, ConicType::None);
        let nb_var = if is_quadric { 1 } else if is_conic { 2 } else { 3 };

        let tol = [tol1, self.tol2, self.tol2];

        // OCCT L348: aNbParticles = 48.
        const NB_PARTICLES: usize = 48;

        // OCCT L347-356: closed/periodic curves spanning > 2/3 period split in
        // two intervals (the "false extrema" of a single sweep).
        let mut a_nb_int_c = 1usize;
        if curve.is_closed() || curve.is_periodic() {
            let a_period = curve.period();
            if curve.last_parameter() - curve.first_parameter() > 2. * a_period / 3. {
                a_nb_int_c = 2;
            }
        }

        let d_t = (my_tsup - my_tmin) / a_nb_int_c as f64;
        for an_int in 1..=a_nb_int_c {
            let tuv_inf = [my_tmin + (an_int as f64 - 1.) * d_t, self.umin, self.vmin];
            let tuv_sup = [tuv_inf[0] + d_t, self.umax, self.vmax];
            let mut tuv = [tuv_inf[0] + 0.5 * d_t, self.umin, self.vmin];

            match nb_var {
                1 => {
                    Self::glob_min_c_quadric(curve, &surf, self.tsample, self.usample,
                        self.vsample, NB_PARTICLES, tuv_inf, tuv_sup, &mut tuv);
                }
                2 => {
                    Self::glob_min_conic_s(curve, &surf, self.tsample, self.usample,
                        self.vsample, NB_PARTICLES, tuv_inf, tuv_sup, &mut tuv);
                }
                _ => {
                    Self::glob_min_gen_cs(curve, &surf, self.tsample, self.usample,
                        self.vsample, NB_PARTICLES, tuv_inf, tuv_sup, &mut tuv);
                }
            }

            // OCCT L384-385: the math_FunctionSetRoot refinement of myF seeded
            // by the PSO result (default NbIterations = 100).
            let mut an_a = FunctionSetRoot::new(&my_f, &tol, 100);
            an_a.perform(&mut my_f, &tuv, &tuv_inf, &tuv_sup, false);
        }

        // OCCT L389-418: remove the "false" extrema caused by dividing the
        // curve interval (keep the solutions within SquareConfusion of the
        // minimum).
        if a_nb_int_c > 1 && my_f.nb_ext() > 1 {
            my_f.filter_false_extrema();
        }

        for i in 1..=my_f.nb_ext() {
            let (t, _p1) = my_f.point_on_curve(i);
            let (u, v, _p2) = my_f.point_on_surface(i);
            self.extrema.push((t, u, v, my_f.square_distance(i)));
        }
        self.is_done = true;
    }

    /// OCCT GlobMinConicS (Extrema_GenExtCS.cxx L522-708) — 2-variable PSO
    /// over the surface UV for a conic curve.
    #[allow(clippy::too_many_arguments)]
    fn glob_min_conic_s(
        curve: &BRepAdaptorCurve,
        surf: &Surface3,
        t_sample: usize,
        u_sample: usize,
        v_sample: usize,
        nb_particles: usize,
        tuvinf: [f64; 3],
        tuvsup: [f64; 3],
        tuv: &mut [f64; 3],
    ) {
        let uvinf = [tuvinf[1], tuvinf[2]];
        let uvsup = [tuvsup[1], tuvsup[2]];
        let mut particles = PsoParticlesPool::new(nb_particles, 2);

        let min_uv: Vec<f64> = (0..2)
            .map(|i| uvinf[i] + (uvsup[i] - uvinf[i]) / BORDER_DIVISOR)
            .collect();
        let max_uv: Vec<f64> = (0..2)
            .map(|i| uvsup[i] - (uvsup[i] - uvinf[i]) / BORDER_DIVISOR)
            .collect();

        // OCCT L546-549: extra UV samples improve the global-minimum search.
        let an_addsample = (t_sample / 2).max(3);
        let an_usample = u_sample + an_addsample;
        let a_vsample = v_sample + an_addsample;
        let step_su = (max_uv[0] - min_uv[0]) / an_usample as f64;
        let step_sv = (max_uv[1] - min_uv[1]) / a_vsample as f64;

        let c3 = curve.curve().clone();
        let mut func = GlobOptFuncConicS {
            surface: surf,
            ctype: curve_type_of(curve).clone_conic(),
            line: match &c3 { Curve3::Line(l) => Some(*l), _ => None },
            circle: match &c3 { Curve3::Circle(c) => Some(*c), _ => None },
            ellipse: match &c3 { Curve3::Ellipse(e) => Some(*e), _ => None },
            hyperbola: match &c3 { Curve3::Hyperbola(h) => Some(*h), _ => None },
            parabola: match &c3 { Curve3::Parabola(p) => Some(*p), _ => None },
            tf: tuvinf[0],
            tl: tuvsup[0],
            umin: uvinf[0],
            umax: uvsup[0],
            vmin: uvinf[1],
            vmax: uvsup[1],
            _marker: std::marker::PhantomData,
        };

        let mut psu = min_uv[0];
        for _sui in 0..=an_usample {
            let mut psv = min_uv[1];
            for _svi in 0..=a_vsample {
                let sq_dist = func.value(psu, psv);
                let a_particle = particles.worst();
                if sq_dist < a_particle.distance {
                    a_particle.position[0] = psu;
                    a_particle.position[1] = psv;
                    a_particle.best_position[0] = psu;
                    a_particle.best_position[1] = psv;
                    a_particle.distance = sq_dist;
                    a_particle.best_distance = sq_dist;
                }
                psv += step_sv;
            }
            psu += step_su;
        }

        let steps = [step_su, step_sv];
        let (_value, uv) = pso_with_given_particles(
            &mut |x: &[f64]| func.value(x[0], x[1]),
            &uvinf,
            &uvsup,
            &steps,
            &mut particles,
            100,
        );

        // OCCT L599: conic parameter of the PSO (u,v) point. The periodic fold
        // (L600-607) is not applicable to a line bean.
        let ct = func.conic_parameter(uv[0], uv[1]);

        tuv[0] = ct;
        tuv[1] = uv[0];
        tuv[2] = uv[1];

        // OCCT L613-707: bad-solution check — the curve point may sit on the
        // surface while far from the PSO surface point (perpendicular
        // PcPs vs normal); refine via a 3x3 neighborhood and the point
        // projection Newton (Extrema_GenLocateExtPS -> rcad
        // closest_point_on_surface_near on the restricted window).
        let cd0 = curve.curve().point_at(tuv[0]);
        let cd1 = curve.curve().derivative_at(tuv[0]);
        let sd0 = surf.point_at(tuv[1], tuv[2]);
        let (sd1u, sd1v) = surface_d1(surf, tuv[1], tuv[2]);
        let sq_dist = cd0.distance_squared(sd0);
        if sq_dist <= rcad_kernel::core::precision::SQUARE_CONFUSION {
            return;
        }
        let n = sd1u.cross(sd1v);
        if n.length_squared() < rcad_kernel::core::precision::SQUARE_CONFUSION {
            return;
        }
        let pcps = cd0 - sd0;
        let ang_min = std::f64::consts::FRAC_PI_2 - std::f64::consts::FRAC_PI_2 / 10.0;
        let ang_max = std::f64::consts::FRAC_PI_2 + std::f64::consts::FRAC_PI_2 / 10.0;
        let ang_n = angle_between(pcps, n);
        if ang_n >= ang_min && ang_n <= ang_max {
            let mut u = tuv[1];
            let mut v = tuv[2];
            let near = crate::extrema::closest_point_on_surface_near(surf, cd0, u, v);
            if near.distance * near.distance < sq_dist
                || angle_between(pcps, n) < ang_min
                || angle_between(pcps, n) > ang_max
            {
                u = near.params.0;
                v = near.params.1;
                tuv[1] = u;
                tuv[2] = v;
            }
        }
    }

    /// OCCT GlobMinGenCS (Extrema_GenExtCS.cxx L427-518) — 3-variable PSO over
    /// (t,u,v) with the pre-computed surface grid.
    #[allow(clippy::too_many_arguments)]
    fn glob_min_gen_cs(
        curve: &BRepAdaptorCurve,
        surf: &Surface3,
        t_sample: usize,
        u_sample: usize,
        v_sample: usize,
        nb_particles: usize,
        tuvinf: [f64; 3],
        tuvsup: [f64; 3],
        tuv: &mut [f64; 3],
    ) {
        let mut particles = PsoParticlesPool::new(nb_particles, 3);

        let min_tuv: Vec<f64> = (0..3)
            .map(|i| tuvinf[i] + (tuvsup[i] - tuvinf[i]) / BORDER_DIVISOR)
            .collect();
        let max_tuv: Vec<f64> = (0..3)
            .map(|i| tuvsup[i] - (tuvsup[i] - tuvinf[i]) / BORDER_DIVISOR)
            .collect();

        let step_cu = (max_tuv[0] - min_tuv[0]) / t_sample as f64;
        let step_su = (max_tuv[1] - min_tuv[1]) / u_sample as f64;
        let step_sv = (max_tuv[2] - min_tuv[2]) / v_sample as f64;

        // OCCT L467-474: pre-compute the curve sample points.
        let mut curv_pnts: Vec<DVec3> = Vec::with_capacity(t_sample + 1);
        let mut cu = min_tuv[0];
        for _ in 0..=t_sample {
            curv_pnts.push(curve.curve().point_at(cu));
            cu += step_cu;
        }

        let mut psu = min_tuv[1];
        for _sui in 0..=u_sample {
            let mut psv = min_tuv[2];
            for _svi in 0..=v_sample {
                let sp = surf.point_at(psu, psv);
                let mut cu2 = min_tuv[0];
                for cp in curv_pnts.iter() {
                    let sq_dist = sp.distance_squared(*cp);
                    let a_particle = particles.worst();
                    if sq_dist < a_particle.distance {
                        a_particle.position[0] = cu2;
                        a_particle.position[1] = psu;
                        a_particle.position[2] = psv;
                        a_particle.best_position.copy_from_slice(&a_particle.position);
                        a_particle.distance = sq_dist;
                        a_particle.best_distance = sq_dist;
                    }
                    cu2 += step_cu;
                }
                psv += step_sv;
            }
            psu += step_su;
        }

        let steps = [step_cu, step_su, step_sv];
        let mut func = |x: &[f64]| -> f64 {
            curve
                .curve()
                .point_at(x[0])
                .distance_squared(surf.point_at(x[1], x[2]))
        };
        let (_value, point) = pso_with_given_particles(
            &mut func,
            &tuvinf,
            &tuvsup,
            &steps,
            &mut particles,
            100,
        );
        tuv.copy_from_slice(&point);
    }

    /// OCCT GlobMinCQuadric (Extrema_GenExtCS.cxx L712-823) — 1-variable PSO
    /// over the curve parameter for quadric surfaces (the quadric point is
    /// obtained in closed form — rcad's analytic closest_point_on_surface).
    #[allow(clippy::too_many_arguments)]
    fn glob_min_c_quadric(
        curve: &BRepAdaptorCurve,
        surf: &Surface3,
        t_sample: usize,
        u_sample: usize,
        _v_sample: usize,
        nb_particles: usize,
        tuvinf: [f64; 3],
        tuvsup: [f64; 3],
        tuv: &mut [f64; 3],
    ) {
        let mut particles = PsoParticlesPool::new(nb_particles, 1);

        let min_t = tuvinf[0] + (tuvsup[0] - tuvinf[0]) / BORDER_DIVISOR;
        let max_t = tuvsup[0] - (tuvsup[0] - tuvinf[0]) / BORDER_DIVISOR;

        // OCCT L732-738: extra curve samples (dimensionality reduced).
        const MAX_NB_NODES: usize = 50;
        let an_addsample = (u_sample / 2).max(3);
        let new_csample = (t_sample + an_addsample).min(MAX_NB_NODES);
        let step_ct = (max_t - min_t) / new_csample as f64;

        let func = |t: f64| -> f64 {
            let pc = curve.curve().point_at(t);
            let proj = crate::extrema::closest_point_on_surface_near(surf, pc, tuvinf[1], tuvinf[2]);
            pc.distance_squared(proj.point)
        };
        let mut pct = min_t;
        for _cui in 0..=new_csample {
            let sq_dist = func(pct);
            let a_particle = particles.worst();
            if sq_dist < a_particle.distance {
                a_particle.position[0] = pct;
                a_particle.best_position[0] = pct;
                a_particle.distance = sq_dist;
                a_particle.best_distance = sq_dist;
            }
            pct += step_ct;
        }
        let steps = [step_ct];
        let (_value, point) = pso_with_given_particles(
            &mut |x: &[f64]| func(x[0]),
            &[tuvinf[0]],
            &[tuvsup[0]],
            &steps,
            &mut particles,
            100,
        );
        // OCCT L800-818: quadric parameters of the PSO t (periodic folds).
        let pc = curve.curve().point_at(point[0]);
        let proj = crate::extrema::closest_point_on_surface_near(surf, pc, tuvinf[1], tuvinf[2]);
        tuv[0] = point[0];
        tuv[1] = proj.params.0;
        tuv[2] = proj.params.1;
    }

    pub fn is_done(&self) -> bool {
        self.is_done
    }

    pub fn nb_ext(&self) -> usize {
        self.extrema.len()
    }

    pub fn square_distance(&self, idx_1based: usize) -> f64 {
        self.extrema[idx_1based - 1].3
    }

    pub fn point(&self, idx_1based: usize) -> (f64, f64, f64) {
        let (t, u, v, _f) = self.extrema[idx_1based - 1];
        (t, u, v)
    }
}

/// OCCT Extrema_GenExtCS::GetCurvMaxParamVal (Extrema_GenExtCS.cxx L94-114)
/// — HyperbolaLimit for a hyperbola-based curve, MaxParamVal otherwise.
fn get_curv_max_param_val(curve: &BRepAdaptorCurve) -> f64 {
    if curve_type_of(curve) == ConicType::Hyperbola {
        return HYPERBOLA_LIMIT;
    }
    MAX_PARAM_VAL
}

/// OCCT Extrema_FuncExtCS (TKGeomBase/Extrema/Extrema_FuncExtCS.cxx L1-238)
/// — the function set
/// { F1(t,u,v) = (C(t)-S(u,v)).Dtc(t) }
/// { F2(t,u,v) = (C(t)-S(u,v)).Dus(u,v) }
/// { F3(t,u,v) = (C(t)-S(u,v)).Dvs(u,v) }
/// whose roots are the curve/surface extremalities.  GetStateNumber records
/// each accepted solution (deduplicated by curve parameter).
pub(crate) struct FuncExtCS<'a> {
    curve: &'a Curve3,
    surface: &'a Surface3,
    my_t: f64,
    my_u: f64,
    my_v: f64,
    my_p1: DVec3,
    my_p2: DVec3,
    my_sq_dist: Vec<f64>,
    // Extrema_POnCurv (parameter, point) / Extrema_POnSurf (u, v, point).
    my_point1: Vec<(f64, DVec3)>,
    my_point2: Vec<(f64, f64, DVec3)>,
}

impl<'a> FuncExtCS<'a> {
    /// OCCT Extrema_FuncExtCS(C, S) -> Initialize(C, S) (cxx L59-77).
    pub fn new(curve: &'a Curve3, surface: &'a Surface3) -> Self {
        let mut f = FuncExtCS {
            curve,
            surface,
            my_t: 0.0,
            my_u: 0.0,
            my_v: 0.0,
            my_p1: DVec3::ZERO,
            my_p2: DVec3::ZERO,
            my_sq_dist: Vec::new(),
            my_point1: Vec::new(),
            my_point2: Vec::new(),
        };
        f.initialize();
        f
    }

    /// OCCT Initialize(C, S) (cxx L66-77).
    pub fn initialize(&mut self) {
        self.my_sq_dist.clear();
        self.my_point1.clear();
        self.my_point2.clear();
    }

    /// OCCT NbExt() (cxx L189-193).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// OCCT SquareDistance(N) (cxx L196-206).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist[n - 1]
    }

    /// OCCT PointOnCurve(N) — (parameter, point) (cxx L209-216).
    pub fn point_on_curve(&self, n: usize) -> (f64, DVec3) {
        self.my_point1[n - 1]
    }

    /// OCCT PointOnSurface(N) — (u, v, point) (cxx L219-226).
    pub fn point_on_surface(&self, n: usize) -> (f64, f64, DVec3) {
        self.my_point2[n - 1]
    }

    /// OCCT GetStateNumber() (cxx L156-186): record the current solution,
    /// deduplicated by curve parameter with Precision::SquarePConfusion.
    fn record_solution(&mut self) {
        let tol2d = rcad_kernel::precision::PCONFUSION * rcad_kernel::precision::PCONFUSION;
        let nb_sol = self.my_sq_dist.len();
        for i in 1..=nb_sol {
            let (a_t, _) = self.my_point1[i - 1];
            let mut d = a_t - self.my_t;
            d *= d;
            if d <= tol2d {
                return;
            }
        }
        self.my_sq_dist.push(self.my_p1.distance_squared(self.my_p2));
        self.my_point1.push((self.my_t, self.my_p1));
        self.my_point2.push((self.my_u, self.my_v, self.my_p2));
    }

    /// OCCT Extrema_GenExtCS::Perform L391-418 — remove the "false" extrema
    /// caused by dividing the curve interval (keep the solutions within
    /// Precision::SquareConfusion of the minimum).
    pub fn filter_false_extrema(&mut self) {
        let sq_dists1 = self.my_sq_dist.clone();
        let pnts_on_crv1 = self.my_point1.clone();
        let pnts_on_surf1 = self.my_point2.clone();

        let mut a_min_dist = f64::MAX;
        for &a_dist in &sq_dists1 {
            if a_dist < a_min_dist {
                a_min_dist = a_dist;
            }
        }
        self.my_sq_dist.clear();
        self.my_point1.clear();
        self.my_point2.clear();
        const A_TOL: f64 = rcad_kernel::precision::SQUARE_CONFUSION;
        for (i, &a_dist) in sq_dists1.iter().enumerate() {
            if (a_dist - a_min_dist).abs() <= A_TOL {
                self.my_sq_dist.push(a_dist);
                self.my_point1.push(pnts_on_crv1[i]);
                self.my_point2.push(pnts_on_surf1[i]);
            }
        }
    }
}

impl<'a> FunctionSetWithDerivatives for FuncExtCS<'a> {
    /// OCCT NbVariables() (cxx L80-83).
    fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() (cxx L86-89).
    fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(UV, F) (cxx L92-122): D1-based evaluations; the vector
    /// P1P2 = myP1 - myP2 = C(t) - S(u,v) (gp_Vec(myP2, myP1)).
    fn value(&mut self, uv: &[f64], f: &mut [f64]) -> bool {
        self.my_t = uv[0];
        self.my_u = uv[1];
        self.my_v = uv[2];

        let dtc = CurveEval::derivative_at(self.curve, self.my_t);
        let (_, dus, dvs) = SurfaceEval::derivatives(self.surface, self.my_u, self.my_v);
        self.my_p1 = CurveEval::point_at(self.curve, self.my_t);
        self.my_p2 = SurfaceEval::point_at(self.surface, self.my_u, self.my_v);

        let p1p2 = self.my_p1 - self.my_p2;

        f[0] = p1p2.dot(dtc);
        f[1] = p1p2.dot(dus);
        f[2] = p1p2.dot(dvs);

        true
    }

    /// OCCT Derivatives(UV, DF) (cxx L127-132) — delegates to Values.
    fn derivatives(&mut self, uv: &[f64], df: &mut [Vec<f64>]) -> bool {
        let mut f = [0.0f64; 3];
        self.values(uv, &mut f, df);
        true
    }

    /// OCCT Values(UV, F, Df) (cxx L137-186): D2-based function set and
    /// Jacobian.
    fn values(&mut self, uv: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        self.my_t = uv[0];
        self.my_u = uv[1];
        self.my_v = uv[2];

        let (dtc, dttc) = (
            CurveEval::derivative_at(self.curve, self.my_t),
            CurveEval::derivative2_at(self.curve, self.my_t),
        );
        // kernel derivatives2 = (P, D1U, D1V, D2UU, D2UV, D2VV); the OCCT D2
        // evaluation order is Dus, Dvs, Duus, Dvvs, Duvs.
        let (_, dus, dvs, duus, duvs, dvvs) =
            SurfaceEval::derivatives2(self.surface, self.my_u, self.my_v);

        self.my_p1 = CurveEval::point_at(self.curve, self.my_t);
        self.my_p2 = SurfaceEval::point_at(self.surface, self.my_u, self.my_v);

        let p1p2 = self.my_p1 - self.my_p2;

        f[0] = p1p2.dot(dtc);
        f[1] = p1p2.dot(dus);
        f[2] = p1p2.dot(dvs);

        df[0][0] = dtc.length_squared() + p1p2.dot(dttc);
        df[0][1] = -dus.dot(dtc);
        df[0][2] = -dvs.dot(dtc);

        df[1][0] = -df[0][1]; // Dtc.Dot(Dus);
        df[1][1] = -dus.length_squared() + p1p2.dot(duus);
        df[1][2] = -dvs.dot(dus) + p1p2.dot(duvs);

        df[2][0] = -df[0][2]; // -Df(1,3)
        df[2][1] = df[1][2]; // -Dus.Dot(Dvs)+P1P2.Dot(Duvs);
        df[2][2] = -dvs.length_squared() + p1p2.dot(dvvs);

        true
    }

    /// OCCT GetStateNumber() (cxx L156-186) — the solver calls this at every
    /// solution-acceptance return; the function records the current state.
    fn get_state_number(&mut self) -> i32 {
        self.record_solution();
        0
    }
}

#[derive(Clone, Copy, PartialEq)]
enum SType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Revolution,
    Other,
}

fn surface_type_of(s: &Surface3) -> SType {
    match s {
        Surface3::Plane(_) => SType::Plane,
        Surface3::Cylinder(_) => SType::Cylinder,
        Surface3::Cone(_) => SType::Cone,
        Surface3::Sphere(_) => SType::Sphere,
        Surface3::Torus(_) => SType::Torus,
        Surface3::Revolution(_) => SType::Revolution,
        _ => SType::Other,
    }
}

fn curve_type_of(c: &BRepAdaptorCurve) -> ConicType {
    match c.curve() {
        Curve3::Line(_) => ConicType::Line,
        Curve3::Circle(_) => ConicType::Circle,
        Curve3::Ellipse(_) => ConicType::Ellipse,
        Curve3::Hyperbola(_) => ConicType::Hyperbola,
        Curve3::Parabola(_) => ConicType::Parabola,
        _ => ConicType::None,
    }
}

impl ConicType {
    fn clone_conic(self) -> ConicType {
        self
    }
}

/// OCCT gp_Vec::Angle (gp_Vec.cxx) — the angle between two vectors.
fn angle_between(a: DVec3, b: DVec3) -> f64 {
    let denom = a.length() * b.length();
    if denom < 1e-30 {
        return 0.0;
    }
    let c = a.dot(b) / denom;
    c.clamp(-1.0, 1.0).acos()
}
