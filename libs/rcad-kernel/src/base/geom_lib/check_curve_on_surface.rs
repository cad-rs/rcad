//! OCCT GeomLib_CheckCurveOnSurface (TKGeomBase/GeomLib —
//! `GeomLib_CheckCurveOnSurface.hxx` L26-82 + `GeomLib_CheckCurveOnSurface.cxx`
//! L33-774): computes the maximal distance between a 3D curve and its
//! representation on a surface (the 2D curve lifted through
//! `Adaptor3d_CurveOnSurface`).
//!
//! Placement (`docs/module-map.md` L22): `GeomLib` is a TKGeomBase package and
//! TKGeomBase packages map to `rcad-kernel/src/base/`; every dependency this
//! class consumes is kernel-native — `base/proj_lib/adaptor.rs`
//! (`Adaptor3d_Curve` / `Adaptor3d_CurveOnSurface` / `Adaptor2d_Curve2d`),
//! `base/extrema_math_opt.rs` (`math_NewtonMinimum`) and `math/bspl_lib.rs`
//! (`BSplCLib::FirstUKnotIndex` / `LastUKnotIndex`).
//!
//! OCCT structure translated 1:1 (the anchor at each item is the OCCT range):
//! - `GeomLib_CheckCurveOnSurface_TargetFunc` (cxx L52-206)
//! - `GeomLib_CheckCurveOnSurface_Local` (cxx L210-283)
//! - the two ctors and the `Init` overloads (cxx L287-331)
//! - `Perform` (cxx L335-426)
//! - `FillSubIntervals` (cxx L430-640), `PSO_Perform` (cxx L644-694) and
//!   `MinComputing` (cxx L698-774)
//! - the TKMath classes the minimisation runs through (`math_PSO`,
//!   `math_PSOParticlesPool`, `math_BullardGenerator`) are re-hosted in the
//!   section marked below.
//!
//! Encodings:
//! - `occ::handle(Adaptor3d_Curve)` -> `GeomCurveHandle`; `occ::handle(Adaptor2d_Curve2d)`
//!   -> `Curve2dHandle`; `occ::handle(Adaptor3d_CurveOnSurface)` ->
//!   `Arc<CurveOnSurface>` (the kernel `Adaptor3d_CurveOnSurface` encoding).
//! - `NCollection_Array1<double>` -> `Vector` (the kernel `math_Vector`
//!   encoding: `Value`/`ChangeValue` -> `get`/`set`, `Lower`/`Upper`/`Length`
//!   identical, so every index expression stays verbatim);
//!   `handle(NCollection_HArray1<double>)` -> `Vector` (OCCT never leaves that
//!   handle null here).
//! - `RealLast()` -> `REAL_LAST`, `Precision::PConfusion()` -> `p_confusion()`,
//!   `Precision::Confusion()` -> `CONFUSION`.
//! - The OCCT `try { } catch (Standard_Failure const&) { }` arms are not
//!   modeled (rcad has no C++-exception channel); the guarded statements are
//!   carried and each catch body is documented where it stood.

use std::sync::Arc;

use glam::DVec3;

use crate::base::extrema_math_opt::{
    vec_copy_assign, MathNewtonMinimum, MathStatus, MultipleVarFunctionWithHessian,
};
use crate::base::proj_lib::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, CurveHandle, CurveOnSurface,
};
use crate::base::proj_lib::geom_adaptor_curve::Adaptor3dCurveGeom;
use crate::base::proj_lib::proj_lib_projected_curve_b::GeomCurveHandle;
use crate::base::proj_lib::CurveType;
use crate::core::precision::{p_confusion, CONFUSION, REAL_LAST};
use crate::geom::{BSplineCurve2, BSplineCurve3};
use crate::math::bspl_lib::{first_uknot_index_mults, last_uknot_index_mults};
use crate::math::math_bfgs::{MultipleVarFunction, MultipleVarFunctionWithGradient};
use crate::math::math_matrix::{Matrix, Vector};

// ===========================================================================
// Re-hosts of the OCCT accessors this translation needs
// ===========================================================================

/// The `NCollection_HArray1<double>` copy of an OCCT knot array
/// (`new NCollection_HArray1<double>(aBS3DCurv->Knots())`, cxx L462 and
/// L487): the rcad 1-based `Vector` carrying the same bounds.
fn array_from_slice(the_values: &[f64]) -> Vector {
    let mut a_result = Vector::new(1, the_values.len() as i32);
    for (an_idx, a_value) in the_values.iter().enumerate() {
        a_result.set(an_idx as i32 + 1, *a_value);
    }
    a_result
}

/// OCCT `Geom2d_BSplineCurve` flat knot sequence -> `(myKnots, myMults)`.
///
/// `Geom2d_BSplineCurve` stores the distinct knots plus their multiplicities
/// and maintains `myFlatKnots` (`KnotSequence()`); the rcad `BSplineCurve2`
/// stores the flat sequence, so the split is the 2D twin of
/// `BSplineCurve3::knots_mults` (kernel `geom/bspline_ops.rs`).  Equal
/// consecutive knots group into one entry (`Geom_BSplineCurve`'s
/// `updateKnots`/`KnotAnalysis` invariant).
fn bspline2_knots_mults(the_curve: &BSplineCurve2) -> (Vec<f64>, Vec<i32>) {
    let mut a_knots = Vec::new();
    let mut a_mults = Vec::new();
    for (an_idx, a_knot) in the_curve.knots.iter().enumerate() {
        if an_idx > 0 && *a_knot == a_knots.last().copied().unwrap_or(f64::NAN) {
            *a_mults.last_mut().unwrap() += 1;
        } else {
            a_knots.push(*a_knot);
            a_mults.push(1);
        }
    }
    (a_knots, a_mults)
}

/// OCCT `Geom_BSplineCurve::FirstUKnotIndex()` (Geom_BSplineCurve_1.cxx
/// L334-344) / the `Geom2d_BSplineCurve` twin (Geom2d_BSplineCurve_1.cxx
/// L335-345).
fn bspline_first_u_knot_index(my_periodic: bool, my_deg: usize, my_mults: &[i32]) -> i32 {
    if my_periodic {
        1
    } else {
        first_uknot_index_mults(my_deg, my_mults)
    }
}

/// OCCT `Geom_BSplineCurve::LastUKnotIndex()` (Geom_BSplineCurve_1.cxx
/// L404-412) / the `Geom2d_BSplineCurve` twin (Geom2d_BSplineCurve_1.cxx
/// L405-415); the periodic arm answers `myKnots.Length() == NbKnots()`.
fn bspline_last_u_knot_index(
    my_periodic: bool,
    my_nb_knots: i32,
    my_deg: usize,
    my_mults: &[i32],
) -> i32 {
    if my_periodic {
        my_nb_knots
    } else {
        last_uknot_index_mults(my_deg, my_mults)
    }
}

/// OCCT `Geom2d_BSplineCurve::IsPeriodic()` derived from the knot structure:
/// a periodic (unclamped) B-spline has the first and the last multiplicity
/// equal to the degree (not Degree + 1).
fn bspline2_is_periodic(the_degree: usize, the_mults: &[i32]) -> bool {
    if the_mults.len() < 3 {
        return false;
    }
    the_mults[0] == the_degree as i32 && *the_mults.last().unwrap() == the_degree as i32
}

/// OCCT `OSD_ThreadPool::DefaultPool()->NbDefaultThreadsToLaunch()` — the
/// default pool launches one task per hardware thread.  rcad carries no
/// thread pool (see `Perform`); the value only sizes the per-thread
/// `ShallowCopy` arrays, so it is re-hosted from the standard library.
fn nb_default_threads_to_launch() -> usize {
    std::thread::available_parallelism()
        .map(|a_count| a_count.get())
        .unwrap_or(1)
}

/// OCCT `Adaptor3d_CurveOnSurface::ShallowCopy()` (Adaptor3d_CurveOnSurface.cxx
/// L901-925) — the missing override in the kernel `CurveOnSurface`.
///
/// OCCT copies the adaptor and then shallow-copies each member handle; the
/// rcad members are already shared `Arc` handles (the `Adaptor3d_CurveOnSurface`
/// handle copies), so the shallow copy is a new adaptor over the same 2D curve
/// and the same surface — the evaluation and the `EvalKPart` cache are
/// recomputed identically.
///
/// Note: the kernel `CurveOnSurface` does not override
/// `Adaptor3d_Curve::ShallowCopy`, whose default raises
/// `Standard_NotImplemented` (Adaptor3d_Curve.cxx L38-42).  The override
/// belongs to `base/proj_lib/adaptor.rs`; until it lands there, the launcher
/// arm of `Perform` reaches it through this function so the parallel and the
/// sequential arms stay behaviourally identical.
fn curve_on_surface_shallow_copy(the_curve_on_surface: &Arc<CurveOnSurface>) -> CurveHandle {
    let a_copy = CurveOnSurface::new(
        Arc::clone(&the_curve_on_surface.my2d_curve),
        Arc::clone(&the_curve_on_surface.my_surface),
    );
    Arc::new(a_copy)
}

/// OCCT `gp_Pnt::SquareDistance` (gp_Pnt.cxx) — the squared modulus of the
/// difference.
fn square_distance(the_p1: DVec3, the_p2: DVec3) -> f64 {
    (the_p1 - the_p2).length_squared()
}

// ===========================================================================
// OCCT math_BullardGenerator / math_PSOParticlesPool / math_PSO (TKMath,
// FoundationClasses/TKMath/math) — the minimiser `MinComputing` runs through.
//
// Interface request: these three classes belong to the TKMath home
// (`rcad-kernel/src/math/`); the partial copy of the same algorithm in
// `rcad-algo/src/bop/int_tools/extrema_gen_ext_cs.rs` (BullardGenerator /
// PsoParticlesPool / pso_with_given_particles) should redirect here once that
// domain opens.
// ===========================================================================

/// OCCT math_BullardGenerator (math_BullardGenerator.hxx L22-54).
struct BullardGenerator {
    /// OCCT unsigned int myStateHi.
    my_state_hi: u32,
    /// OCCT unsigned int myStateLo.
    my_state_lo: u32,
}

impl BullardGenerator {
    /// OCCT math_BullardGenerator(theSeed = 1) (hxx L26-30) followed by
    /// SetSeed(theSeed) (hxx L32-37).
    fn new(the_seed: u32) -> Self {
        let mut a_this = BullardGenerator {
            my_state_hi: the_seed,
            my_state_lo: 0,
        };
        a_this.set_seed(the_seed);
        a_this
    }

    /// OCCT SetSeed(theSeed = 1) (hxx L32-37).
    fn set_seed(&mut self, the_seed: u32) {
        self.my_state_hi = the_seed;
        self.my_state_lo = the_seed ^ 0x49616E42;
    }

    /// OCCT NextInt() (hxx L40-48) — unsigned 32-bit wrapping arithmetic.
    fn next_int(&mut self) -> u32 {
        self.my_state_hi = self
            .my_state_hi
            .wrapping_shr(2)
            .wrapping_add(self.my_state_hi.wrapping_shl(2));
        self.my_state_hi = self.my_state_hi.wrapping_add(self.my_state_lo);
        self.my_state_lo = self.my_state_lo.wrapping_add(self.my_state_hi);
        self.my_state_hi
    }

    /// OCCT NextReal() (hxx L51).
    fn next_real(&mut self) -> f64 {
        self.next_int() as f64 / 0xFFFF_FFFFu32 as f64
    }
}

/// OCCT PSO_Particle (math_PSOParticlesPool.hxx L24-45).
struct PsoParticle {
    /// OCCT double* Position.
    position: Vec<f64>,
    /// OCCT double* Velocity.
    velocity: Vec<f64>,
    /// OCCT double* BestPosition.
    best_position: Vec<f64>,
    /// OCCT double Distance.
    distance: f64,
    /// OCCT double BestDistance.
    best_distance: f64,
}

impl PsoParticle {
    /// OCCT PSO_Particle() (hxx L32-39): both distances start at RealLast().
    fn new(the_dimension_count: usize) -> Self {
        PsoParticle {
            position: vec![0.0; the_dimension_count],
            velocity: vec![0.0; the_dimension_count],
            best_position: vec![0.0; the_dimension_count],
            distance: REAL_LAST,
            best_distance: REAL_LAST,
        }
    }
}

/// OCCT math_PSOParticlesPool (math_PSOParticlesPool.hxx L49-71).
struct PsoParticlesPool {
    /// OCCT NCollection_Array1<PSO_Particle> myParticlesPool — the pool is
    /// 1-based in OCCT (`myParticlesPool(1, theParticlesCount)`); the
    /// accessors below keep the OCCT indices.
    particles: Vec<PsoParticle>,
}

impl PsoParticlesPool {
    /// OCCT math_PSOParticlesPool(theParticlesCount, theDimensionCount)
    /// (math_PSOParticlesPool.cxx L19-39): myMemory is initialized to 0.
    fn new(the_particles_count: usize, the_dimension_count: usize) -> Self {
        let mut a_particles = Vec::with_capacity(the_particles_count);
        for _ in 0..the_particles_count {
            a_particles.push(PsoParticle::new(the_dimension_count));
        }
        PsoParticlesPool {
            particles: a_particles,
        }
    }

    /// OCCT GetParticle(theIdx) (cxx L45-48); 1 <= theIdx <= myParticlesCount.
    fn get_particle(&mut self, the_idx: i32) -> &mut PsoParticle {
        &mut self.particles[(the_idx - 1) as usize]
    }

    /// OCCT GetBestParticle() (cxx L52-55) — `std::min_element`, i.e. the
    /// first particle of minimal Distance (the pool orders by Distance,
    /// math_PSOParticlesPool.hxx L42).
    fn get_best_particle(&self) -> &PsoParticle {
        let mut a_best = 0usize;
        for an_idx in 1..self.particles.len() {
            if self.particles[an_idx].distance < self.particles[a_best].distance {
                a_best = an_idx;
            }
        }
        &self.particles[a_best]
    }

    /// OCCT GetWorstParticle() (cxx L58-61) — `std::max_element`.
    fn get_worst_particle(&mut self) -> &mut PsoParticle {
        let mut a_worst = 0usize;
        for an_idx in 1..self.particles.len() {
            if self.particles[an_idx].distance > self.particles[a_worst].distance {
                a_worst = an_idx;
            }
        }
        &mut self.particles[a_worst]
    }
}

/// OCCT math_PSO.cxx L21: `const double aBorderDivisor = 1.0e+4`.
const A_BORDER_DIVISOR: f64 = 1.0e+4;

/// OCCT math_PSO (math_PSO.hxx L55-104) — a variation of Particle Swarm
/// Optimization.  Only the pool-based `Perform` overload (cxx L56-63) is
/// carried: it is the one `GeomLib_CheckCurveOnSurface::PSO_Perform` calls.
/// The regular-grid `Perform(theSteps, theValue, theOutPnt, theNbIter)`
/// overload (cxx L68-123) is not on this code path.
struct MathPso<'a> {
    /// OCCT math_MultipleVarFunction* myFunc.
    my_func: &'a mut dyn MultipleVarFunction,
    /// OCCT math_Vector myLowBorder.
    my_low_border: Vector,
    /// OCCT math_Vector myUppBorder.
    my_upp_border: Vector,
    /// OCCT math_Vector mySteps.
    my_steps: Vector,
    /// OCCT int myN.
    my_n: usize,
    /// OCCT int myNbParticles.
    my_nb_particles: usize,
    /// OCCT int myNbIter.
    my_nb_iter: i32,
}

impl<'a> MathPso<'a> {
    /// OCCT math_PSO(theFunc, theLowBorder, theUppBorder, theSteps,
    /// theNbParticles = 32, theNbIter = 100) (math_PSO.cxx L31-52).
    fn new(
        the_func: &'a mut dyn MultipleVarFunction,
        the_low_border: &Vector,
        the_upp_border: &Vector,
        the_steps: &Vector,
        the_nb_particles: usize,
        the_nb_iter: i32,
    ) -> Self {
        let my_n = the_func.nb_variables() as usize;
        let mut my_low_border = Vector::new(1, my_n as i32);
        let mut my_upp_border = Vector::new(1, my_n as i32);
        let mut my_steps = Vector::new(1, my_n as i32);
        vec_copy_assign(&mut my_low_border, the_low_border);
        vec_copy_assign(&mut my_upp_border, the_upp_border);
        vec_copy_assign(&mut my_steps, the_steps);
        MathPso {
            my_func: the_func,
            my_low_border,
            my_upp_border,
            my_steps,
            my_n,
            my_nb_particles: the_nb_particles,
            my_nb_iter: the_nb_iter,
        }
    }

    /// OCCT math_PSO::Perform(theParticles, theNbParticles, theValue,
    /// theOutPnt, theNbIter = 100) (math_PSO.cxx L56-63).
    fn perform_given_particles(
        &mut self,
        the_particles: &mut PsoParticlesPool,
        the_nb_particles: usize,
        the_value: &mut f64,
        the_out_pnt: &mut Vector,
        the_nb_iter: i32,
    ) {
        self.perform_pso_with_given_particles(
            the_particles,
            the_nb_particles,
            the_value,
            the_out_pnt,
            the_nb_iter,
        );
    }

    /// OCCT math_PSO::performPSOWithGivenParticles (math_PSO.cxx L125-268).
    fn perform_pso_with_given_particles(
        &mut self,
        the_particles: &mut PsoParticlesPool,
        the_nb_particles: usize,
        the_value: &mut f64,
        the_out_pnt: &mut Vector,
        the_nb_iter: i32,
    ) {
        // cxx L131-135.
        let mut a_min_uv = Vector::new(1, self.my_n as i32);
        let mut a_max_uv = Vector::new(1, self.my_n as i32);
        for an_idx in 1..=self.my_n as i32 {
            let a_low = self.my_low_border.get(an_idx);
            let a_upp = self.my_upp_border.get(an_idx);
            a_min_uv.set(an_idx, a_low + (a_upp - a_low) / A_BORDER_DIVISOR);
            a_max_uv.set(an_idx, a_upp - (a_upp - a_low) / A_BORDER_DIVISOR);
        }
        self.my_nb_iter = the_nb_iter;
        self.my_nb_particles = the_nb_particles;

        // cxx L137-139.
        let mut a_curr_point = Vector::new(1, self.my_n as i32);
        let mut a_best_global_position = Vector::new(1, self.my_n as i32);

        // cxx L148-158: math_BullardGenerator aRandom; — the default seed 1.
        let mut a_random = BullardGenerator::new(1);
        for a_part_idx in 1..=self.my_nb_particles {
            let a_particle = the_particles.get_particle(a_part_idx as i32);
            for a_dim_idx in 0..self.my_n {
                let a_ksi = a_random.next_real();
                a_particle.velocity[a_dim_idx] =
                    self.my_steps.get(a_dim_idx as i32 + 1) * (a_ksi - 0.5) * 2.0;
            }
        }

        // cxx L160-166.
        let a_particle = the_particles.get_best_particle();
        for a_dim_idx in 0..self.my_n {
            a_best_global_position.set(a_dim_idx as i32 + 1, a_particle.position[a_dim_idx]);
        }
        let mut a_best_global_distance = a_particle.distance;

        // cxx L168-174.
        let mut a_termination_velocity = Vector::new(1, self.my_n as i32);
        for a_dim_idx in 1..=self.my_n as i32 {
            a_termination_velocity.set(a_dim_idx, self.my_steps.get(a_dim_idx) / 2048.0);
        }
        let mut a_minimal_velocity = Vector::new(1, self.my_n as i32);

        // cxx L177-179: for (int aStep = 1; aStep < myNbIter; ++aStep).
        let mut a_step = 1i32;
        while a_step < self.my_nb_iter {
            // cxx L181: aMinimalVelocity.Init(RealLast()).
            for a_dim_idx in 1..=self.my_n as i32 {
                a_minimal_velocity.set(a_dim_idx, REAL_LAST);
            }

            for a_part_idx in 1..=self.my_nb_particles {
                // cxx L185-186.
                let a_ksi1 = a_random.next_real();
                let a_ksi2 = a_random.next_real();

                let a_particle = the_particles.get_particle(a_part_idx as i32);

                // cxx L190-192.
                const A_RETENT_WEIGHT: f64 = 0.72900;
                const A_PERSON_WEIGHT: f64 = 1.49445;
                const A_SOCIAL_WEIGHT: f64 = 1.49445;

                for a_dim_idx in 0..self.my_n {
                    a_particle.velocity[a_dim_idx] = a_particle.velocity[a_dim_idx]
                        * A_RETENT_WEIGHT
                        + (a_particle.best_position[a_dim_idx] - a_particle.position[a_dim_idx])
                            * (A_PERSON_WEIGHT * a_ksi1)
                        + (a_best_global_position.get(a_dim_idx as i32 + 1)
                            - a_particle.position[a_dim_idx])
                            * (A_SOCIAL_WEIGHT * a_ksi2);

                    a_particle.position[a_dim_idx] += a_particle.velocity[a_dim_idx];
                    a_particle.position[a_dim_idx] = a_particle.position[a_dim_idx]
                        .min(a_max_uv.get(a_dim_idx as i32 + 1))
                        .max(a_min_uv.get(a_dim_idx as i32 + 1));
                    a_curr_point.set(a_dim_idx as i32 + 1, a_particle.position[a_dim_idx]);

                    let a_abs_velocity = a_particle.velocity[a_dim_idx].abs();
                    let a_minimum = a_abs_velocity.min(a_minimal_velocity.get(a_dim_idx as i32 + 1));
                    a_minimal_velocity.set(a_dim_idx as i32 + 1, a_minimum);
                }

                // cxx L211: myFunc->Value(aCurrPoint, aParticle->Distance); —
                // OCCT discards the returned flag; the out value is only
                // written when the objective succeeds.
                let _ = self.my_func.value(&a_curr_point, &mut a_particle.distance);
                if a_particle.distance < a_particle.best_distance {
                    a_particle.best_distance = a_particle.distance;
                    for a_dim_idx in 0..self.my_n {
                        a_particle.best_position[a_dim_idx] = a_particle.position[a_dim_idx];
                    }

                    if a_particle.distance < a_best_global_distance {
                        a_best_global_distance = a_particle.distance;
                        for a_dim_idx in 0..self.my_n {
                            a_best_global_position
                                .set(a_dim_idx as i32 + 1, a_particle.position[a_dim_idx]);
                        }
                    }
                }
            }

            // cxx L231-241.
            let mut is_terminal_velocity_reached = true;
            for a_dim_idx in 1..=self.my_n as i32 {
                if a_minimal_velocity.get(a_dim_idx) > a_termination_velocity.get(a_dim_idx) {
                    is_terminal_velocity_reached = false;
                    break;
                }
            }

            if is_terminal_velocity_reached {
                // cxx L244-252: the minimum number of steps.
                const A_MIN_STEPS: i32 = 16;

                if a_step > A_MIN_STEPS {
                    break;
                }

                for a_part_idx in 1..=self.my_nb_particles {
                    let a_ksi = a_random.next_real();

                    let a_particle = the_particles.get_particle(a_part_idx as i32);

                    for a_dim_idx in 0..self.my_n {
                        if a_particle.position[a_dim_idx] == a_min_uv.get(a_dim_idx as i32 + 1)
                            || a_particle.position[a_dim_idx]
                                == a_max_uv.get(a_dim_idx as i32 + 1)
                        {
                            a_particle.velocity[a_dim_idx] = if a_particle.position[a_dim_idx]
                                == a_min_uv.get(a_dim_idx as i32 + 1)
                            {
                                self.my_steps.get(a_dim_idx as i32 + 1) * a_ksi
                            } else {
                                -self.my_steps.get(a_dim_idx as i32 + 1) * a_ksi
                            };
                        } else {
                            a_particle.velocity[a_dim_idx] =
                                self.my_steps.get(a_dim_idx as i32 + 1) * (a_ksi - 0.5) * 2.0;
                        }
                    }
                }
            }

            a_step += 1;
        }

        // cxx L266-267.
        *the_value = a_best_global_distance;
        vec_copy_assign(the_out_pnt, &a_best_global_position);
    }
}

// ===========================================================================
// GeomLib_CheckCurveOnSurface_TargetFunc (cxx L52-206)
// ===========================================================================

/// OCCT GeomLib_CheckCurveOnSurface_TargetFunc — the objective of the
/// minimisation: `F(X) = -|C3D(X) - C2D(X)|^2` (minimising the negated square
/// distance maximises the distance).
struct GeomLibCheckCurveOnSurfaceTargetFunc<'a> {
    /// OCCT const Adaptor3d_Curve& myCurve1.
    my_curve1: &'a dyn Adaptor3dCurve,
    /// OCCT const Adaptor3d_Curve& myCurve2.
    my_curve2: &'a dyn Adaptor3dCurve,
    /// OCCT const double myFirst.
    my_first: f64,
    /// OCCT const double myLast.
    my_last: f64,
}

impl<'a> GeomLibCheckCurveOnSurfaceTargetFunc<'a> {
    /// OCCT TargetFunc(theC3D, theCurveOnSurface, theFirst, theLast)
    /// (cxx L55-64).
    fn new(
        the_c3d: &'a dyn Adaptor3dCurve,
        the_curve_on_surface: &'a dyn Adaptor3dCurve,
        the_first: f64,
        the_last: f64,
    ) -> Self {
        GeomLibCheckCurveOnSurfaceTargetFunc {
            my_curve1: the_c3d,
            my_curve2: the_curve_on_surface,
            my_first: the_first,
            my_last: the_last,
        }
    }

    /// OCCT TargetFunc::Value(theX, theFVal) const (cxx L75-94).  The OCCT
    /// `try { } catch (Standard_Failure const&) { return false; }` arm maps to
    /// the `false` returns of the adaptor calls.
    fn value_scalar(&self, the_x: f64, the_f_val: &mut f64) -> bool {
        if !self.check_parameter(the_x) {
            return false;
        }

        let a_p1 = self.my_curve1.value(the_x);
        let a_p2 = self.my_curve2.value(the_x);

        *the_f_val = -1.0 * square_distance(a_p1, a_p2);
        true
    }

    /// OCCT TargetFunc::Derive(theX, theDeriv1, theDeriv2 = nullptr) const
    /// (cxx L108-148).
    fn derive(&self, the_x: f64, the_deriv1: &mut f64, the_deriv2: Option<&mut f64>) -> bool {
        if !self.check_parameter(the_x) {
            return false;
        }

        // OCCT: gp_Pnt aP1, aP2; gp_Vec aDC1, aDC2, aDCC1, aDCC2; — the unused
        // pair stays zero for the D1 branch.
        let (a_p1, a_dc1, a_dcc1);
        let (a_p2, a_dc2, a_dcc2);
        if the_deriv2.is_none() {
            let (a_data_p1, a_data_dc1) = self.my_curve1.d1(the_x);
            let (a_data_p2, a_data_dc2) = self.my_curve2.d1(the_x);
            a_p1 = a_data_p1;
            a_dc1 = a_data_dc1;
            a_dcc1 = DVec3::ZERO;
            a_p2 = a_data_p2;
            a_dc2 = a_data_dc2;
            a_dcc2 = DVec3::ZERO;
        } else {
            let (a_data_p1, a_data_dc1, a_data_dcc1) = self.my_curve1.d2(the_x);
            let (a_data_p2, a_data_dc2, a_data_dcc2) = self.my_curve2.d2(the_x);
            a_p1 = a_data_p1;
            a_dc1 = a_data_dc1;
            a_dcc1 = a_data_dcc1;
            a_p2 = a_data_p2;
            a_dc2 = a_data_dc2;
            a_dcc2 = a_data_dcc2;
        }

        let a_vec1 = a_p2 - a_p1;
        let a_vec2 = a_dc2 - a_dc1;

        *the_deriv1 = -2.0 * a_vec1.dot(a_vec2);

        if let Some(a_deriv2) = the_deriv2 {
            let a_vec3 = a_dcc2 - a_dcc1;
            *a_deriv2 = -2.0 * (a_vec2.length_squared() + a_vec1.dot(a_vec3));
        }

        true
    }

    /// OCCT TargetFunc::FirstParameter() (cxx L186).
    fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT TargetFunc::LastParameter() (cxx L189).
    fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT TargetFunc::CheckParameter(theParam) const (cxx L197-200).
    fn check_parameter(&self, the_param: f64) -> bool {
        (self.my_first <= the_param) && (the_param <= self.my_last)
    }
}

impl MultipleVarFunction for GeomLibCheckCurveOnSurfaceTargetFunc<'_> {
    /// OCCT TargetFunc::NbVariables() (cxx L68) — the function is
    /// one-dimensional.
    fn nb_variables(&self) -> i32 {
        1
    }

    /// OCCT TargetFunc::Value(theX, theFVal) (cxx L71) — the math_Vector
    /// overload delegates to the one-dimensional form with theX(1).
    fn value(&mut self, the_x: &Vector, the_f_val: &mut f64) -> bool {
        self.value_scalar(the_x.get(1), the_f_val)
    }
}

impl MultipleVarFunctionWithGradient for GeomLibCheckCurveOnSurfaceTargetFunc<'_> {
    /// OCCT TargetFunc::Gradient(theX, theGrad) (cxx L101-104).
    fn gradient(&mut self, the_x: &Vector, the_grad: &mut Vector) -> bool {
        let mut a_deriv1 = 0.0;
        let is_done = self.derive(the_x.get(1), &mut a_deriv1, None);
        if is_done {
            the_grad.set(1, a_deriv1);
        }
        is_done
    }

    /// OCCT TargetFunc::Values(theX, theVal, theGrad) (cxx L151-164).
    fn values(&mut self, the_x: &Vector, the_val: &mut f64, the_grad: &mut Vector) -> bool {
        if !self.value(the_x, the_val) {
            return false;
        }
        if !self.gradient(the_x, the_grad) {
            return false;
        }
        true
    }
}

impl MultipleVarFunctionWithHessian for GeomLibCheckCurveOnSurfaceTargetFunc<'_> {
    /// OCCT TargetFunc::Values(theX, theVal, theGrad, theHessian)
    /// (cxx L167-183).
    fn values_hessian(
        &mut self,
        the_x: &Vector,
        the_val: &mut f64,
        the_grad: &mut Vector,
        the_hessian: &mut Matrix,
    ) -> bool {
        if !self.value(the_x, the_val) {
            return false;
        }

        let mut a_deriv1 = 0.0;
        let mut a_deriv2 = 0.0;
        if !self.derive(the_x.get(1), &mut a_deriv1, Some(&mut a_deriv2)) {
            return false;
        }
        the_grad.set(1, a_deriv1);
        the_hessian.set(1, 1, a_deriv2);
        true
    }
}

// ===========================================================================
// GeomLib_CheckCurveOnSurface_Local (cxx L210-283)
// ===========================================================================

/// OCCT GeomLib_CheckCurveOnSurface_Local — the per-sub-interval work item.
///
/// The OCCT `mutable NCollection_Array1<double> myArrOfDist / myArrOfParam`
/// (written through the const functor) map to `&mut self` accessors here.
struct GeomLibCheckCurveOnSurfaceLocal<'a> {
    /// OCCT const Array1OfHCurve& myCurveArray.
    my_curve_array: &'a [CurveHandle],
    /// OCCT const Array1OfHCurve& myCurveOnSurfaceArray.
    my_curve_on_surface_array: &'a [CurveHandle],
    /// OCCT const NCollection_Array1<double>& mySubIntervals.
    my_sub_intervals: &'a Vector,
    /// OCCT const double myEpsilonRange.
    my_epsilon_range: f64,
    /// OCCT const int myNbParticles.
    my_nb_particles: i32,
    /// OCCT mutable NCollection_Array1<double> myArrOfDist.
    my_arr_of_dist: Vector,
    /// OCCT mutable NCollection_Array1<double> myArrOfParam.
    my_arr_of_param: Vector,
}

impl<'a> GeomLibCheckCurveOnSurfaceLocal<'a> {
    /// OCCT Local(theCurveArray, theCurveOnSurfaceArray, theIntervalsArr,
    /// theEpsilonRange, theNbParticles) (cxx L213-226).
    fn new(
        the_curve_array: &'a [CurveHandle],
        the_curve_on_surface_array: &'a [CurveHandle],
        the_intervals_arr: &'a Vector,
        the_epsilon_range: f64,
        the_nb_particles: i32,
    ) -> Self {
        GeomLibCheckCurveOnSurfaceLocal {
            my_curve_array: the_curve_array,
            my_curve_on_surface_array: the_curve_on_surface_array,
            my_sub_intervals: the_intervals_arr,
            my_epsilon_range: the_epsilon_range,
            my_nb_particles: the_nb_particles,
            my_arr_of_dist: Vector::new(
                the_intervals_arr.lower(),
                the_intervals_arr.upper() - 1,
            ),
            my_arr_of_param: Vector::new(
                the_intervals_arr.lower(),
                the_intervals_arr.upper() - 1,
            ),
        }
    }

    /// OCCT Local::operator()(theThreadIndex, theElemIndex) const
    /// (cxx L228-251).
    fn perform(&mut self, the_thread_index: i32, the_elem_index: i32) {
        let a_curve_array: &'a [CurveHandle] = self.my_curve_array;
        let a_curve_on_surface_array: &'a [CurveHandle] = self.my_curve_on_surface_array;
        let a_sub_intervals: &'a Vector = self.my_sub_intervals;

        let mut a_func = GeomLibCheckCurveOnSurfaceTargetFunc::new(
            &*a_curve_array[the_thread_index as usize],
            &*a_curve_on_surface_array[the_thread_index as usize],
            a_sub_intervals.get(the_elem_index),
            a_sub_intervals.get(the_elem_index + 1),
        );

        let mut a_min_dist = REAL_LAST;
        let mut a_par = 0.0;
        if !min_computing(
            &mut a_func,
            self.my_epsilon_range,
            self.my_nb_particles,
            &mut a_min_dist,
            &mut a_par,
        ) {
            self.my_arr_of_dist.set(the_elem_index, REAL_LAST);
            self.my_arr_of_param
                .set(the_elem_index, a_func.first_parameter());
            return;
        }

        self.my_arr_of_dist.set(the_elem_index, a_min_dist);
        self.my_arr_of_param.set(the_elem_index, a_par);
    }

    /// OCCT Local::OptimalValues(theMinimalValue, theParameter) const
    /// (cxx L254-269).
    fn optimal_values(&self, the_minimal_value: &mut f64, the_parameter: &mut f64) {
        let a_start_ind = self.my_arr_of_dist.lower();
        *the_minimal_value = self.my_arr_of_dist.get(a_start_ind);
        *the_parameter = self.my_arr_of_param.get(a_start_ind);
        let mut an_idx = a_start_ind + 1;
        while an_idx <= self.my_arr_of_dist.upper() {
            if self.my_arr_of_dist.get(an_idx) < *the_minimal_value {
                *the_minimal_value = self.my_arr_of_dist.get(an_idx);
                *the_parameter = self.my_arr_of_param.get(an_idx);
            }
            an_idx += 1;
        }
    }
}

// ===========================================================================
// GeomLib_CheckCurveOnSurface (hxx L26-82)
// ===========================================================================

/// OCCT GeomLib_CheckCurveOnSurface (hxx L26-82) — computes the max distance
/// between the 3D curve `myCurve` and the curve-on-surface passed to
/// `Perform`.
pub struct GeomLibCheckCurveOnSurface {
    /// OCCT occ::handle(Adaptor3d_Curve) myCurve — the null handle is `None`.
    my_curve: Option<GeomCurveHandle>,
    /// OCCT int myErrorStatus.
    my_error_status: i32,
    /// OCCT double myMaxDistance.
    my_max_distance: f64,
    /// OCCT double myMaxParameter.
    my_max_parameter: f64,
    /// OCCT double myTolRange.
    my_tol_range: f64,
    /// OCCT bool myIsParallel.
    my_is_parallel: bool,
}

impl GeomLibCheckCurveOnSurface {
    /// OCCT GeomLib_CheckCurveOnSurface() (cxx L287-294).
    pub fn new() -> Self {
        GeomLibCheckCurveOnSurface {
            my_curve: None,
            my_error_status: 0,
            my_max_distance: REAL_LAST,
            my_max_parameter: 0.0,
            my_tol_range: p_confusion(),
            my_is_parallel: false,
        }
    }

    /// OCCT GeomLib_CheckCurveOnSurface(theCurve, theTolRange =
    /// Precision::PConfusion()) (cxx L298-308).  Rust has no default
    /// arguments: the OCCT default is `p_confusion()`.
    pub fn with_curve(the_curve: &GeomCurveHandle, the_tol_range: f64) -> Self {
        GeomLibCheckCurveOnSurface {
            my_curve: Some(Arc::clone(the_curve)),
            my_error_status: 0,
            my_max_distance: REAL_LAST,
            my_max_parameter: 0.0,
            my_tol_range: the_tol_range,
            my_is_parallel: false,
        }
    }

    /// OCCT GeomLib_CheckCurveOnSurface::Init() (cxx L312-319) — initializes
    /// all members by default values.  `myIsParallel` is left untouched, as in
    /// OCCT.
    pub fn init(&mut self) {
        self.my_curve = None;
        self.my_error_status = 0;
        self.my_max_distance = REAL_LAST;
        self.my_max_parameter = 0.0;
        self.my_tol_range = p_confusion();
    }

    /// OCCT GeomLib_CheckCurveOnSurface::Init(theCurve, theTolRange =
    /// Precision::PConfusion()) (cxx L323-331) — sets the data for the
    /// algorithm.
    pub fn init_curve(&mut self, the_curve: &GeomCurveHandle, the_tol_range: f64) {
        self.my_curve = Some(Arc::clone(the_curve));
        self.my_error_status = 0;
        self.my_max_distance = REAL_LAST;
        self.my_max_parameter = 0.0;
        self.my_tol_range = the_tol_range;
    }

    /// OCCT GeomLib_CheckCurveOnSurface::SetParallel(theIsParallel)
    /// (hxx L51).
    pub fn set_parallel(&mut self, the_is_parallel: bool) {
        self.my_is_parallel = the_is_parallel;
    }

    /// OCCT GeomLib_CheckCurveOnSurface::IsParallel() (hxx L54) — the OCCT
    /// declaration is non-const.
    pub fn is_parallel(&mut self) -> bool {
        self.my_is_parallel
    }

    /// OCCT GeomLib_CheckCurveOnSurface::Perform(theCurveOnSurface)
    /// (cxx L335-426).
    pub fn perform(&mut self, the_curve_on_surface: &Arc<CurveOnSurface>) {
        // cxx L338-342: the rcad `Arc<CurveOnSurface>` handle is never null,
        // so only the 3D curve slot can be empty.
        if self.my_curve.is_none() {
            self.my_error_status = 1;
            return;
        }

        // cxx L344-349.
        let my_curve: GeomCurveHandle = self.my_curve.as_ref().unwrap().clone();
        if (my_curve.first_parameter() - the_curve_on_surface.first_parameter()
            > self.my_tol_range)
            || (my_curve.last_parameter() - the_curve_on_surface.last_parameter()
                < -self.my_tol_range)
        {
            self.my_error_status = 2;
            return;
        }

        // cxx L351-359.
        let an_epsilon_range = 1.0e-3;

        let mut a_nb_particles: i32 = 3;

        // cxx L362-366.  The polynomial-function remark above the call.
        let a_nb_sub_intervals = fill_sub_intervals(
            &*my_curve,
            &**the_curve_on_surface.get_curve(),
            my_curve.first_parameter(),
            my_curve.last_parameter(),
            &mut a_nb_particles,
            None,
        );

        // cxx L368-372.
        if a_nb_sub_intervals == 0 {
            self.my_error_status = 3;
            return;
        }

        // cxx L374-421 — the OCCT try block; its
        // `catch (Standard_Failure const&) { myErrorStatus = 3; }` arm has no
        // rcad channel (see the module header).
        //
        // cxx L378: NCollection_Array1<double> anIntervals(1, aNbSubIntervals + 1);
        let mut an_intervals = Vector::new(1, a_nb_sub_intervals + 1);
        fill_sub_intervals(
            &*my_curve,
            &**the_curve_on_surface.get_curve(),
            my_curve.first_parameter(),
            my_curve.last_parameter(),
            &mut a_nb_particles,
            Some(&mut an_intervals),
        );

        // cxx L386-389.
        let a_nb_threads: usize = if self.my_is_parallel {
            std::cmp::min(
                an_intervals.length() as usize,
                nb_default_threads_to_launch(),
            )
        } else {
            1
        };

        // cxx L390-399.  The OCCT `aNbThreads > 1` arm replaces both curves by
        // their ShallowCopy; the single-thread arm hands over the handles
        // themselves (a handle copy to the same object).
        let mut a_curve_array: Vec<CurveHandle> = Vec::with_capacity(a_nb_threads);
        let mut a_curve_on_surface_array: Vec<CurveHandle> = Vec::with_capacity(a_nb_threads);
        for _an_i in 0..a_nb_threads {
            if a_nb_threads > 1 {
                a_curve_array.push(my_curve.shallow_copy());
                a_curve_on_surface_array.push(curve_on_surface_shallow_copy(the_curve_on_surface));
            } else {
                // OCCT: static_cast<const occ::handle(Adaptor3d_Curve)&>(myCurve)
                // and the handle copy of theCurveOnSurface — the same objects,
                // seen through the base handle type (the rcad trait-object
                // upcast keeps the handle pointing at the same adaptor).
                let a_curve_handle: CurveHandle = my_curve.clone();
                let a_curve_on_surface_copy = Arc::clone(the_curve_on_surface);
                let a_curve_on_surface_handle: CurveHandle = a_curve_on_surface_copy;
                a_curve_array.push(a_curve_handle);
                a_curve_on_surface_array.push(a_curve_on_surface_handle);
            }
        }

        // cxx L400-404.
        let mut a_comp = GeomLibCheckCurveOnSurfaceLocal::new(
            &a_curve_array,
            &a_curve_on_surface_array,
            &an_intervals,
            an_epsilon_range,
            a_nb_particles,
        );

        // cxx L405-417.  OCCT splits the interval range between aNbThreads
        // OSD_ThreadPool tasks (Launcher::Perform) or walks it sequentially;
        // rcad carries no thread pool, so the launcher's contiguous task split
        // collapses into the sequential pass over the same range with launcher
        // index 0 (the per-interval results are thread-local and
        // OptimalValues scans every entry, so both arms agree).
        let mut an_i = an_intervals.lower();
        while an_i < an_intervals.upper() {
            a_comp.perform(0, an_i);
            an_i += 1;
        }

        // cxx L418-420.
        a_comp.optimal_values(&mut self.my_max_distance, &mut self.my_max_parameter);

        self.my_max_distance = self.my_max_distance.abs().sqrt();
    }

    /// OCCT GeomLib_CheckCurveOnSurface::IsDone() (hxx L57).
    pub fn is_done(&self) -> bool {
        self.my_error_status == 0
    }

    /// OCCT GeomLib_CheckCurveOnSurface::ErrorStatus() (hxx L65).
    /// 0 - OK; 1 - null curve or surface or 2d curve; 2 - invalid parametric
    /// range; 3 - error in calculations.
    pub fn error_status(&self) -> i32 {
        self.my_error_status
    }

    /// OCCT GeomLib_CheckCurveOnSurface::MaxDistance() (hxx L68).
    pub fn max_distance(&self) -> f64 {
        self.my_max_distance
    }

    /// OCCT GeomLib_CheckCurveOnSurface::MaxParameter() (hxx L71).
    pub fn max_parameter(&self) -> f64 {
        self.my_max_parameter
    }
}

impl Default for GeomLibCheckCurveOnSurface {
    fn default() -> Self {
        Self::new()
    }
}

// ===========================================================================
// static FillSubIntervals (cxx L430-640)
// ===========================================================================

/// OCCT static FillSubIntervals(theCurve3d, theCurve2d, theFirst, theLast,
/// theNbParticles, theSubIntervals = nullptr) (cxx L430-640) — builds the
/// sorted sub-interval boundaries of the two curves' knot tables.
fn fill_sub_intervals(
    the_curve3d: &dyn Adaptor3dCurveGeom,
    the_curve2d: &dyn Adaptor2dCurve2d,
    the_first: f64,
    the_last: f64,
    the_nb_particles: &mut i32,
    mut the_sub_intervals: Option<&mut Vector>,
) -> i32 {
    // cxx L437-441.
    const A_MAX_KNOTS: i32 = 101;
    let mut an_arr_temp = Vector::new(1, 2);
    an_arr_temp.set(1, the_first);
    an_arr_temp.set(2, the_last);

    *the_nb_particles = 3;
    let mut a_bs2d_curv: Option<BSplineCurve2> = None;
    let mut a_bs3d_curv: Option<BSplineCurve3> = None;
    let is_trimmed3d = false;
    let is_trimmed2d = false;

    // cxx L447-454.
    if the_curve3d.get_type() == CurveType::BSpline {
        a_bs3d_curv = Some(the_curve3d.bspline());
    }
    if the_curve2d.get_type() == CurveType::BSpline {
        a_bs2d_curv = the_curve2d.bspline();
    }

    // cxx L456-507: the 3D knot array.
    let mut an_arr_knots3d: Vector;
    if let Some(a_curve) = &a_bs3d_curv {
        let (a_knots, a_mults) = a_curve.knots_mults();
        let a_nb_knots = a_knots.len() as i32;
        if a_nb_knots <= A_MAX_KNOTS {
            // cxx L462: anArrKnots3D = new NCollection_HArray1<double>(aBS3DCurv->Knots());
            an_arr_knots3d = array_from_slice(&a_knots);
        } else {
            // cxx L466-484.
            let a_knot_count: i32;
            if is_trimmed3d {
                // cxx L468-480.
                a_knot_count = {
                    let mut a_count = 0i32;
                    let mut an_i = bspline_first_u_knot_index(
                        a_curve.is_periodic,
                        a_curve.degree,
                        &a_mults,
                    );
                    let a_last_u_knot = bspline_last_u_knot_index(
                        a_curve.is_periodic,
                        a_nb_knots,
                        a_curve.degree,
                        &a_mults,
                    );
                    while an_i <= a_last_u_knot {
                        let a_knot = a_knots[(an_i - 1) as usize];
                        if a_knot > the_first && a_knot < the_last {
                            a_count += 1;
                        }
                        an_i += 1;
                    }
                    a_count + 2
                };
            } else {
                // cxx L483: KnotCount = LastUKnotIndex() - FirstUKnotIndex() + 1.
                a_knot_count = bspline_last_u_knot_index(
                    a_curve.is_periodic,
                    a_nb_knots,
                    a_curve.degree,
                    &a_mults,
                ) - bspline_first_u_knot_index(
                    a_curve.is_periodic,
                    a_curve.degree,
                    &a_mults,
                ) + 1;
            }
            if a_knot_count <= A_MAX_KNOTS {
                // cxx L487.
                an_arr_knots3d = array_from_slice(&a_knots);
            } else {
                // cxx L491-500.
                an_arr_knots3d = Vector::new(1, A_MAX_KNOTS);
                an_arr_knots3d.set(1, the_first);
                an_arr_knots3d.set(A_MAX_KNOTS, the_last);
                let a_dt = (the_last - the_first) / (A_MAX_KNOTS - 1) as f64;
                let mut t = the_first + a_dt;
                let mut an_i = 2;
                while an_i < A_MAX_KNOTS {
                    an_arr_knots3d.set(an_i, t);
                    an_i += 1;
                    t += a_dt;
                }
            }
        }
    } else {
        // cxx L506: anArrKnots3D = new NCollection_HArray1<double>(anArrTemp);
        an_arr_knots3d = an_arr_temp.clone();
    }

    // cxx L508-557: the 2D knot array.
    let mut an_arr_knots2d: Vector;
    if let Some(a_curve) = &a_bs2d_curv {
        let (a_knots, a_mults) = bspline2_knots_mults(a_curve);
        let a_nb_knots = a_knots.len() as i32;
        if a_nb_knots <= A_MAX_KNOTS {
            // cxx L512.
            an_arr_knots2d = array_from_slice(&a_knots);
        } else {
            // cxx L516-534.
            let a_knot_count: i32;
            if is_trimmed2d {
                // cxx L518-530.
                a_knot_count = {
                    let mut a_count = 0i32;
                    let mut an_i = bspline_first_u_knot_index(
                        bspline2_is_periodic(a_curve.degree, &a_mults),
                        a_curve.degree,
                        &a_mults,
                    );
                    let a_last_u_knot = bspline_last_u_knot_index(
                        bspline2_is_periodic(a_curve.degree, &a_mults),
                        a_nb_knots,
                        a_curve.degree,
                        &a_mults,
                    );
                    while an_i <= a_last_u_knot {
                        let a_knot = a_knots[(an_i - 1) as usize];
                        if a_knot > the_first && a_knot < the_last {
                            a_count += 1;
                        }
                        an_i += 1;
                    }
                    a_count + 2
                };
            } else {
                // cxx L533.
                a_knot_count = bspline_last_u_knot_index(
                    bspline2_is_periodic(a_curve.degree, &a_mults),
                    a_nb_knots,
                    a_curve.degree,
                    &a_mults,
                ) - bspline_first_u_knot_index(
                    bspline2_is_periodic(a_curve.degree, &a_mults),
                    a_curve.degree,
                    &a_mults,
                ) + 1;
            }
            if a_knot_count <= A_MAX_KNOTS {
                // cxx L537.
                an_arr_knots2d = array_from_slice(&a_knots);
            } else {
                // cxx L541-550.
                an_arr_knots2d = Vector::new(1, A_MAX_KNOTS);
                an_arr_knots2d.set(1, the_first);
                an_arr_knots2d.set(A_MAX_KNOTS, the_last);
                let a_dt = (the_last - the_first) / (A_MAX_KNOTS - 1) as f64;
                let mut t = the_first + a_dt;
                let mut an_i = 2;
                while an_i < A_MAX_KNOTS {
                    an_arr_knots2d.set(an_i, t);
                    an_i += 1;
                    t += a_dt;
                }
            }
        }
    } else {
        // cxx L556: anArrKnots2D = new NCollection_HArray1<double>(anArrTemp);
        an_arr_knots2d = an_arr_temp.clone();
    }

    // cxx L559.
    let mut a_nb_sub_intervals = 1;

    // cxx L561-627 — the OCCT try block; its
    // `catch (Standard_Failure const&) { aNbSubIntervals = 0; }` arm has no
    // rcad channel (see the module header).
    let an_ind_max3d = an_arr_knots3d.upper();
    let an_ind_max2d = an_arr_knots2d.upper();

    let mut an_index3d = an_arr_knots3d.lower();
    let mut an_index2d = an_arr_knots2d.lower();

    if let Some(a_sub_intervals) = the_sub_intervals.as_mut() {
        a_sub_intervals.set(a_nb_sub_intervals, the_first);
    }

    while (an_index3d <= an_ind_max3d) && (an_index2d <= an_ind_max2d) {
        let a_val3d = an_arr_knots3d.get(an_index3d);
        let a_val2d = an_arr_knots2d.get(an_index2d);
        let a_delta = a_val3d - a_val2d;

        if a_delta < p_confusion() {
            // aVal3D <= aVal2D
            if (a_val3d > the_first) && (a_val3d < the_last) {
                a_nb_sub_intervals += 1;

                if let Some(a_sub_intervals) = the_sub_intervals.as_mut() {
                    a_sub_intervals.set(a_nb_sub_intervals, a_val3d);
                }
            }

            an_index3d += 1;

            if -a_delta < p_confusion() {
                // aVal3D == aVal2D
                an_index2d += 1;
            }
        } else {
            // aVal2D < aVal3D
            if (a_val2d > the_first) && (a_val2d < the_last) {
                a_nb_sub_intervals += 1;

                if let Some(a_sub_intervals) = the_sub_intervals.as_mut() {
                    a_sub_intervals.set(a_nb_sub_intervals, a_val2d);
                }
            }

            an_index2d += 1;
        }
    }

    if let Some(a_sub_intervals) = the_sub_intervals.as_mut() {
        a_sub_intervals.set(a_nb_sub_intervals + 1, the_last);
    }

    // cxx L618-626.
    if let Some(a_curve) = &a_bs3d_curv {
        *the_nb_particles = (*the_nb_particles).max(a_curve.degree as i32);
    }

    if let Some(a_curve) = &a_bs2d_curv {
        *the_nb_particles = (*the_nb_particles).max(a_curve.degree as i32);
    }

    // cxx L639.
    a_nb_sub_intervals
}

// ===========================================================================
// PSO_Perform (cxx L644-694)
// ===========================================================================

/// OCCT PSO_Perform(theFunction, theParInf, theParSup, theEpsilon,
/// theNbParticles, theBestValue, theOutputParam) (cxx L644-694) — seeds the
/// particle pool from a regular grid of `3 * theNbParticles` control points
/// and runs `math_PSO` on it.
fn pso_perform(
    the_function: &mut GeomLibCheckCurveOnSurfaceTargetFunc<'_>,
    the_par_inf: &Vector,
    the_par_sup: &Vector,
    the_epsilon: f64,
    the_nb_particles: i32,
    the_best_value: &mut f64,
    the_output_param: &mut Vector,
) -> bool {
    // cxx L652-656.
    let a_delta_param = the_par_sup.get(1) - the_par_inf.get(1);
    if a_delta_param < p_confusion() {
        return false;
    }

    // cxx L658-659.
    let mut a_step_par = Vector::new(1, 1);
    a_step_par.set(1, the_epsilon * a_delta_param);

    let mut a_particles = PsoParticlesPool::new(the_nb_particles as usize, 1);

    // cxx L662-664.
    let a_nb_control_points = 3 * the_nb_particles;

    // cxx L666-688.
    let a_step = a_delta_param / (a_nb_control_points - 1) as f64;
    let mut a_count = 1i32;
    let mut a_prm = the_par_inf.get(1);
    while a_count <= a_nb_control_points {
        let mut a_val = REAL_LAST;
        if the_function.value_scalar(a_prm, &mut a_val) {
            let a_particle = a_particles.get_worst_particle();

            // OCCT: if (aVal > aParticle->BestDistance) continue; — the loop
            // increment below must still run, hence the negated form.
            if a_val <= a_particle.best_distance {
                a_particle.position[0] = a_prm;
                a_particle.best_position[0] = a_prm;
                a_particle.distance = a_val;
                a_particle.best_distance = a_val;
            }
        }

        a_count += 1;
        a_prm = if a_count == a_nb_control_points {
            the_par_sup.get(1)
        } else {
            a_prm + a_step
        };
    }

    // cxx L690-691.  The OCCT defaults are theNbParticles = 32 and
    // theNbIter = 100; Perform receives the NbIter default explicitly.
    let mut a_pso = MathPso::new(the_function, the_par_inf, the_par_sup, &a_step_par, 32, 100);
    a_pso.perform_given_particles(
        &mut a_particles,
        the_nb_particles as usize,
        the_best_value,
        the_output_param,
        100,
    );

    true
}

// ===========================================================================
// MinComputing (cxx L698-774)
// ===========================================================================

/// OCCT static MinComputing(theFunction, theEpsilon, theNbParticles,
/// theBestValue, theBestParameter) (cxx L698-774) — locates the minimum of the
/// target function with `math_PSO` and refines it with `math_NewtonMinimum`.
fn min_computing(
    the_function: &mut GeomLibCheckCurveOnSurfaceTargetFunc<'_>,
    the_epsilon: f64,
    the_nb_particles: i32,
    the_best_value: &mut f64,
    the_best_parameter: &mut f64,
) -> bool {
    // cxx L704-771 — the OCCT try block; its
    // `catch (Standard_Failure const&) { return false; }` arm has no rcad
    // channel (see the module header).

    // cxx L708-713.
    let mut a_par_inf = Vector::new(1, 1);
    let mut a_par_sup = Vector::new(1, 1);
    let mut an_output_param = Vector::new(1, 1);
    a_par_inf.set(1, the_function.first_parameter());
    a_par_sup.set(1, the_function.last_parameter());
    *the_best_parameter = a_par_inf.get(1);
    *the_best_value = REAL_LAST;

    // cxx L715-727.
    if !pso_perform(
        the_function,
        &a_par_inf,
        &a_par_sup,
        the_epsilon,
        the_nb_particles,
        the_best_value,
        &mut an_output_param,
    ) {
        // "BRepLib_CheckCurveOnSurface::Compute(): math_PSO is failed!"
        return false;
    }

    *the_best_parameter = an_output_param.get(1);

    // cxx L731-734.  math_NewtonMinimum aMinSol(theFunction) — the OCCT
    // defaults are theTolerance = Precision::Confusion(), theNbIterations =
    // 40, theConvexity = 1.0e-6, theWithSingularity = true.
    let mut a_min_sol = MathNewtonMinimum::new(&*the_function, CONFUSION, 40, 1.0e-6, true);
    a_min_sol.perform(the_function, &an_output_param);

    if a_min_sol.is_done() && (a_min_sol.get_status() == MathStatus::Ok) {
        // cxx L736-741: math_NewtonMinimum has precised the value.
        let a_location = a_min_sol.location();
        vec_copy_assign(&mut an_output_param, a_location);
        *the_best_parameter = an_output_param.get(1);
        *the_best_value = a_min_sol.minimum();
    } else {
        // cxx L743-762: use math_PSO again but on a smaller range.
        let a_step = the_epsilon * (a_par_sup.get(1) - a_par_inf.get(1));
        a_par_inf.set(1, *the_best_parameter - 0.5 * a_step);
        a_par_sup.set(1, *the_best_parameter + 0.5 * a_step);

        let mut a_value = REAL_LAST;
        if pso_perform(
            the_function,
            &a_par_inf,
            &a_par_sup,
            the_epsilon,
            the_nb_particles,
            &mut a_value,
            &mut an_output_param,
        ) {
            if a_value < *the_best_value {
                *the_best_value = a_value;
                *the_best_parameter = an_output_param.get(1);
            }
        }
    }

    // cxx L773.
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::proj_lib::geom_adaptor_surface::GeomSurfaceAdaptor;
    use crate::base::proj_lib::{Geom2dCurveAdaptor, GeomCurveAdaptor};
    use crate::geom::{Curve2d, Curve3, Line2d, Line3, Plane, Surface3};

    /// The rcad encoding of `handle(Adaptor3d_Curve)`: a
    /// `GeomAdaptor_Curve(C, f, l)` handle.
    fn curve_handle(the_curve: Curve3, the_first: f64, the_last: f64) -> GeomCurveHandle {
        Arc::new(GeomCurveAdaptor::with_range(
            the_curve,
            the_first,
            the_last,
        ))
    }

    /// The rcad encoding of `handle(Adaptor3d_CurveOnSurface)`: a
    /// `Geom2dAdaptor_Curve(C2d, f, l)` over a `GeomAdaptor_Surface(S)`.
    fn curve_on_surface_handle(
        the_curve2d: Curve2d,
        the_first: f64,
        the_last: f64,
        the_surface: Surface3,
    ) -> Arc<CurveOnSurface> {
        Arc::new(CurveOnSurface::new(
            Arc::new(Geom2dCurveAdaptor::with_range(
                the_curve2d,
                the_first,
                the_last,
            )),
            Arc::new(GeomSurfaceAdaptor::new(the_surface)),
        ))
    }

    /// An exact pcurve on the surface: the maximal distance is zero.
    #[test]
    fn check_curve_on_surface_exact_pcurve() {
        let a_curve = curve_handle(
            Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X)),
            0.0,
            1.0,
        );
        let a_curve_on_surface = curve_on_surface_handle(
            Curve2d::Line(Line2d::new(glam::DVec2::ZERO, glam::DVec2::X)),
            0.0,
            1.0,
            Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)),
        );

        let mut a_check = GeomLibCheckCurveOnSurface::with_curve(&a_curve, p_confusion());
        a_check.perform(&a_curve_on_surface);

        assert!(a_check.is_done());
        assert_eq!(a_check.error_status(), 0);
        assert!(a_check.max_distance() < CONFUSION);
    }

    /// A pcurve whose lift is offset by 0.1 along the surface normal: the
    /// maximal distance is that offset.
    #[test]
    fn check_curve_on_surface_known_offset() {
        let a_curve = curve_handle(
            Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 0.1), DVec3::X)),
            0.0,
            1.0,
        );
        let a_curve_on_surface = curve_on_surface_handle(
            Curve2d::Line(Line2d::new(glam::DVec2::ZERO, glam::DVec2::X)),
            0.0,
            1.0,
            Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)),
        );

        let mut a_check = GeomLibCheckCurveOnSurface::with_curve(&a_curve, p_confusion());
        a_check.perform(&a_curve_on_surface);

        assert!(a_check.is_done());
        assert!((a_check.max_distance() - 0.1).abs() < 1.0e-6);
        assert!(a_check.max_parameter() >= 0.0 && a_check.max_parameter() <= 1.0);
    }

    /// The OCCT failure paths: a default-constructed checker has no curve
    /// (myErrorStatus = 1) and a parametric-range mismatch answers 2.
    #[test]
    fn check_curve_on_surface_error_paths() {
        let a_curve_on_surface = curve_on_surface_handle(
            Curve2d::Line(Line2d::new(glam::DVec2::ZERO, glam::DVec2::X)),
            0.0,
            1.0,
            Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)),
        );

        let mut a_check = GeomLibCheckCurveOnSurface::new();
        a_check.perform(&a_curve_on_surface);
        assert!(!a_check.is_done());
        assert_eq!(a_check.error_status(), 1);

        // OCCT L344-349: the parametric range is invalid when the
        // curve-on-surface sticks out of the 3D curve range.
        let a_curve = curve_handle(
            Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X)),
            0.0,
            1.0,
        );
        let a_wider_curve_on_surface = curve_on_surface_handle(
            Curve2d::Line(Line2d::new(glam::DVec2::ZERO, glam::DVec2::X)),
            0.0,
            2.0,
            Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)),
        );
        let mut a_check = GeomLibCheckCurveOnSurface::with_curve(&a_curve, p_confusion());
        a_check.perform(&a_wider_curve_on_surface);
        assert!(!a_check.is_done());
        assert_eq!(a_check.error_status(), 2);

        // Init() clears the curve slot again.
        a_check.init();
        assert_eq!(a_check.error_status(), 0);
        assert_eq!(a_check.max_distance(), REAL_LAST);
    }
}
