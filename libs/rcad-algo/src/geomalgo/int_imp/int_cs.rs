//! OCCT IntImp_IntCS (TKGeomAlgo IntImp package).
//!
//! 1:1 translation of `IntImp_IntCS.gxx` (L26-198) — the curve/surface
//! root solver over [`super::zer_cs_par_func::ZerCSParFunc`], including the
//! 3-attempt w-restart loop (w, then w0, then w1).

use glam::DVec3;
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};

use super::{
    is_infinite, CurveTool3d, PSurfaceTool, PRECISION_CONFUSION, PRECISION_SQUARE_CONFUSION,
};

/// OCCT IntImp_IntCS — root of the function set F(u,v,w) = 0, i.e. a
/// curve/surface intersection point.  `PS`/`CT` are the OCCT
/// ThePSurfaceTool/TheCurveTool template parameters (kept as phantom type
/// parameters, exactly like the OCCT template instantiation).
pub struct IntCS<'a, S: ?Sized, C: ?Sized, PS, CT, F>
where
    PS: PSurfaceTool<Surface = S>,
    CT: CurveTool3d<Curve = C>,
    F: FunctionSetWithDerivatives + ZerCSAccessors<S, C>,
{
    done: bool,
    empty: bool,
    my_function: &'a mut F,
    w: f64,
    u: f64,
    v: f64,
    tol: f64,
    _types: std::marker::PhantomData<fn(&S, &C) -> (PS, CT)>,
}

/// Access to the ZerCS function's cached state and its auxillary geometry
/// (OCCT: myFunction.Root()/.Point()/.AuxillarSurface()/.AuxillarCurve()).
pub trait ZerCSAccessors<S: ?Sized, C: ?Sized> {
    fn root(&self) -> f64;
    fn point(&self) -> DVec3;
    fn auxillar_surface(&self) -> &S;
    fn auxillar_curve(&self) -> &C;
}

impl<'a, S: ?Sized, C: ?Sized, PS, CT, F> IntCS<'a, S, C, PS, CT, F>
where
    PS: PSurfaceTool<Surface = S>,
    CT: CurveTool3d<Curve = C>,
    F: FunctionSetWithDerivatives + ZerCSAccessors<S, C>,
{
    /// OCCT IntImp_IntCS(U, V, W, F, TolTangency, MarginCoef) — gxx L26-78:
    /// solves immediately with the surface bounds widened by MarginCoef.
    pub fn with_margin(
        u: f64,
        v: f64,
        w: f64,
        my_function: &'a mut F,
        tol_tangency: f64,
        margin_coef: f64,
    ) -> Self {
        let mut cs = IntCS {
            done: true,
            empty: true,
            my_function,
            w: 0.0,
            u: 0.0,
            v: 0.0,
            tol: tol_tangency * tol_tangency,
            _types: std::marker::PhantomData,
        };
        if cs.tol < PRECISION_SQUARE_CONFUSION {
            cs.tol = PRECISION_SQUARE_CONFUSION;
        }
        let mut rsnld = FunctionSetRoot::new(cs.my_function, &[0.0; 3], 100);
        let s = cs.my_function.auxillar_surface();
        let c = cs.my_function.auxillar_curve();

        let w0 = CT::first_parameter(c);
        let w1 = CT::last_parameter(c);

        let mut u0 = PS::first_u_parameter(s);
        let mut v0 = PS::first_v_parameter(s);
        let mut u1 = PS::last_u_parameter(s);
        let mut v1 = PS::last_v_parameter(s);

        if margin_coef > 0.0 {
            if !is_infinite(u0) && !is_infinite(u1) {
                let mut marg = (u1 - u0) * margin_coef;
                if u0 > u1 {
                    marg = -marg;
                }
                u0 -= marg;
                u1 += marg;
            }
            if !is_infinite(v0) && !is_infinite(v1) {
                let mut marg = (v1 - v0) * margin_coef;
                if v0 > v1 {
                    marg = -marg;
                }
                v0 -= marg;
                v1 += marg;
            }
        }

        cs.perform(u, v, w, &mut rsnld, u0, u1, v0, v1, w0, w1);
        cs
    }

    /// OCCT IntImp_IntCS(F, TolTangency) — gxx L80-90: deferred Perform.
    pub fn new(my_function: &'a mut F, tol_tangency: f64) -> Self {
        let mut cs = IntCS {
            done: true,
            empty: true,
            my_function,
            w: 0.0,
            u: 0.0,
            v: 0.0,
            tol: tol_tangency * tol_tangency,
            _types: std::marker::PhantomData,
        };
        if cs.tol < PRECISION_SQUARE_CONFUSION {
            cs.tol = PRECISION_SQUARE_CONFUSION;
        }
        cs
    }

    /// OCCT Perform(U, V, W, Rsnld, u0, u1, v0, v1, w0, w1) — gxx L92-153.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        u: f64,
        v: f64,
        w: f64,
        rsnld: &mut FunctionSetRoot,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        w0: f64,
        w1: f64,
    ) {
        self.done = true;
        let mut uvap = [u, v, w];
        let s = self.my_function.auxillar_surface();
        let c = self.my_function.auxillar_curve();

        let born_inf = [u0, v0, w0];
        let born_sup = [u1, v1, w1];

        let tolerance = [
            PS::u_resolution(&s, PRECISION_CONFUSION),
            PS::v_resolution(&s, PRECISION_CONFUSION),
            CT::resolution(&c, PRECISION_CONFUSION),
        ];
        rsnld.set_tolerance(&tolerance);
        let mut autretentative = 0;
        self.done = false;
        loop {
            if autretentative == 1 {
                uvap[2] = w0;
            } else if autretentative == 2 {
                uvap[2] = w1;
            }
            autretentative += 1;
            rsnld.perform(
                self.my_function,
                &uvap,
                &born_inf,
                &born_sup,
                false,
            );
            if rsnld.is_done() {
                let absmy_functionroot = self.my_function.root().abs();
                if absmy_functionroot <= self.tol {
                    self.root_internal(&mut uvap, rsnld);
                    self.u = uvap[0];
                    self.v = uvap[1];
                    self.w = uvap[2];
                    self.empty = false;
                    self.done = true;
                }
            }
            if !(self.done == false && autretentative < 3) {
                break;
            }
        }
    }

    /// OCCT `Rsnld.Root(UVap)` inside the loop — reads the solved root into
    /// UVap.  The solver is re-borrowed from the caller.
    fn root_internal(&mut self, uvap: &mut [f64; 3], rsnld: &mut FunctionSetRoot) {
        if rsnld.is_done() {
            let r = rsnld.root();
            uvap.copy_from_slice(&r);
        }
    }

    /// OCCT IsDone() — gxx L155-158.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsEmpty() — gxx L160-165.
    pub fn is_empty(&self) -> bool {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        self.empty
    }

    /// OCCT Point() — gxx L167-174.
    pub fn point(&self) -> DVec3 {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        if self.empty {
            panic!("Standard_DomainError");
        }
        self.my_function.point()
    }

    /// OCCT ParameterOnSurface(U, V) — gxx L176-184.
    pub fn parameter_on_surface(&self) -> (f64, f64) {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        if self.empty {
            panic!("Standard_DomainError");
        }
        (self.u, self.v)
    }

    /// OCCT ParameterOnCurve() — gxx L186-193.
    pub fn parameter_on_curve(&self) -> f64 {
        if !self.done {
            panic!("StdFail_NotDone");
        }
        if self.empty {
            panic!("Standard_DomainError");
        }
        self.w
    }

    /// OCCT Function() — gxx L195-198.
    pub fn function(&mut self) -> &mut F {
        self.my_function
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::int_imp::zer_cs_par_func::ZerCSParFunc;

    /// Plane z = 0 parametrised as S(u, v) = (u, v, 0).
    struct TestPlane;
    struct PlaneTool;
    impl PSurfaceTool for PlaneTool {
        type Surface = TestPlane;
        fn value(_s: &TestPlane, u: f64, v: f64) -> DVec3 {
            DVec3::new(u, v, 0.0)
        }
        // OCCT contract: D1's P output IS the surface point S(U, V).
        fn d1(s: &TestPlane, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
            (PlaneTool::value(s, u, v), DVec3::X, DVec3::Y)
        }
        fn first_u_parameter(_s: &TestPlane) -> f64 {
            -1.0e3
        }
        fn last_u_parameter(_s: &TestPlane) -> f64 {
            1.0e3
        }
        fn first_v_parameter(_s: &TestPlane) -> f64 {
            -1.0e3
        }
        fn last_v_parameter(_s: &TestPlane) -> f64 {
            1.0e3
        }
        fn u_resolution(_s: &TestPlane, r3d: f64) -> f64 {
            r3d
        }
        fn v_resolution(_s: &TestPlane, r3d: f64) -> f64 {
            r3d
        }
    }

    /// 3D line C(w) = (w, 0, w).
    struct TestLine;
    struct LineTool;
    impl CurveTool3d for LineTool {
        type Curve = TestLine;
        fn value(_c: &TestLine, w: f64) -> DVec3 {
            DVec3::new(w, 0.0, w)
        }
        // OCCT contract: D1's P output IS the curve point C(W).
        fn d1(c: &TestLine, w: f64) -> (DVec3, DVec3) {
            (LineTool::value(c, w), DVec3::new(1.0, 0.0, 1.0))
        }
        fn first_parameter(_c: &TestLine) -> f64 {
            -1.0e3
        }
        fn last_parameter(_c: &TestLine) -> f64 {
            1.0e3
        }
        fn resolution(_c: &TestLine, r3d: f64) -> f64 {
            r3d
        }
    }

    /// OCCT IntImp_IntCS root: the line (w,0,w) pierces z=0 at w=0, so the
    /// solution is (u,v,w) = (0,0,0) and the point is the origin.
    #[test]
    fn int_cs_solves_plane_line_hit() {
        let plane = TestPlane;
        let line = TestLine;
        let mut fct = ZerCSParFunc::<TestPlane, TestLine, PlaneTool, LineTool>::new(&plane, &line);
        // OCCT keeps Rsnld and myFunction as separate handles; the solver's
        // ctor borrow ends immediately (it only reads the dimensions).
        let mut rsnld = FunctionSetRoot::new(&fct, &[0.0; 3], 100);
        let mut cs =
            IntCS::<TestPlane, TestLine, PlaneTool, LineTool, _>::new(&mut fct, 1e-7);
        cs.perform(0.3, 0.2, 0.4, &mut rsnld, -1e3, 1e3, -1e3, 1e3, -1e3, 1e3);
        assert!(cs.is_done());
        assert!(!cs.is_empty());
        let (u, v) = cs.parameter_on_surface();
        let w = cs.parameter_on_curve();
        assert!((u - w).abs() < 1e-7, "u={u} v={v} w={w}");
        assert!(v.abs() < 1e-7);
        assert!(w.abs() < 1e-5, "w={w}");
        assert!(cs.point().z.abs() < 1e-5);
    }

    /// Direct solver test (bypasses the IntCS wrapper) to localise a
    /// failure.
    #[test]
    fn dbg_direct_root() {
        let plane = TestPlane;
        let line = TestLine;
        let mut fct =
            ZerCSParFunc::<TestPlane, TestLine, PlaneTool, LineTool>::new(&plane, &line);
        let mut rsnld = FunctionSetRoot::new(&fct, &[1e-7, 1e-7, 1e-7], 100);
        rsnld.perform(
            &mut fct,
            &[0.3, 0.2, 0.4],
            &[-1e3, -1e3, -1e3],
            &[1e3, 1e3, 1e3],
            false,
        );
        assert!(rsnld.is_done(), "not done");
        eprintln!("kount={}", rsnld.nb_iterations());
        let r = rsnld.root();
        assert!((r[0] - r[2]).abs() < 1e-7 && r[2].abs() < 1e-5, "r={r:?}");
    }
}
