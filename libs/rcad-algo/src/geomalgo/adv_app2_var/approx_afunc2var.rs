//! OCCT AdvApp2Var_ApproxAFunc2Var (AdvApp2Var_ApproxAFunc2Var.hxx + .cxx +
//! .lxx) - performs the approximation of F(U,V).
//!
//! Encoding notes (architecture):
//! - `occ::handle<NCollection_HArray1<double>>` members ->
//!   `Option<Vec<f64>>`; `handle<NCollection_HArray1<handle<Geom_Surface>>>`
//!   -> `Option<Vec<Option<BSplineSurface>>>` (the inner None models a null
//!   Geom_Surface handle / a failed down_cast to BSpline).
//! - `AdvApprox_Cutting::Value(A, B, D)` (bool out + cutting value) maps to
//!   the rcad kernel `Cutting::value(A, B) -> Option<f64>`; `(more, dec)` =
//!   `(is_some(), unwrapped)`.
//! - `throw Standard_ConstructionError / Standard_OutOfRange` -> `panic!`.
//! - `Dump(Standard_OStream&)` -> `dump(&self) -> String` (rcad dump
//!   convention, cf. bop/ds/common_block.rs).

use std::fmt::Write;

use rcad_kernel::base::proj_lib::proj_lib_projected_curve::IsoType;
use rcad_kernel::geom::BSplineSurface;
use rcad_kernel::math::adv_approx::Cutting;
use rcad_kernel::math::GeomAbsShape;

use super::approxf2var_c::EvaluatorFunc2Var;
use super::context::Context;
use super::criterion::{Criterion, CriterionRepartition, CriterionType};
use super::framework::Framework;
use super::nc_array::Array2;
use super::network::Network;
use super::node::Node;
use super::patch::Patch;
use glam::DVec2;

/// OCCT AdvApp2Var_ApproxAFunc2Var (AdvApp2Var_ApproxAFunc2Var.hxx L85-283).
pub struct ApproxAFunc2Var {
    /// hxx L244: int myNumSubSpaces[3].
    my_num_sub_spaces: [i32; 3],
    /// hxx L245: my1DTolerances.
    my1d_tolerances: Vec<f64>,
    /// hxx L246: my2DTolerances.
    my2d_tolerances: Vec<f64>,
    /// hxx L247: my3DTolerances.
    my3d_tolerances: Vec<f64>,
    /// hxx L248: my1DTolOnFront.
    my1d_tol_on_front: Array2<f64>,
    /// hxx L249: my2DTolOnFront.
    my2d_tol_on_front: Array2<f64>,
    /// hxx L250: my3DTolOnFront.
    my3d_tol_on_front: Array2<f64>,
    /// hxx L251: myFirstParInU.
    my_first_par_in_u: f64,
    /// hxx L252: myLastParInU.
    my_last_par_in_u: f64,
    /// hxx L253: myFirstParInV.
    my_first_par_in_v: f64,
    /// hxx L254: myLastParInV.
    my_last_par_in_v: f64,
    /// hxx L255: myFavoriteIso.
    my_favorite_iso: IsoType,
    /// hxx L256: myContInU.
    my_cont_in_u: GeomAbsShape,
    /// hxx L257: myContInV.
    my_cont_in_v: GeomAbsShape,
    /// hxx L258: myPrecisionCode.
    my_precision_code: i32,
    /// hxx L259: myMaxDegInU.
    my_max_deg_in_u: i32,
    /// hxx L260: myMaxDegInV.
    my_max_deg_in_v: i32,
    /// hxx L261: myMaxPatches.
    my_max_patches: i32,
    /// hxx L262: myConditions.
    my_conditions: Context,
    /// hxx L263: myResult.
    my_result: Network,
    /// hxx L264: myConstraints.
    my_constraints: Framework,
    /// hxx L265: myDone.
    my_done: bool,
    /// hxx L266: myHasResult.
    my_has_result: bool,
    /// hxx L267: mySurfaces (HArray1 of Geom_Surface handles).
    my_surfaces: Option<Vec<Option<BSplineSurface>>>,
    /// hxx L268: myDegreeInU.
    my_degree_in_u: i32,
    /// hxx L269: myDegreeInV.
    my_degree_in_v: i32,
    /// hxx L270: my1DMaxError.
    my1d_max_error: Option<Vec<f64>>,
    /// hxx L271: my1DAverageError.
    my1d_average_error: Option<Vec<f64>>,
    /// hxx L272: my1DUFrontError.
    my1d_u_front_error: Option<Vec<f64>>,
    /// hxx L273: my1DVFrontError.
    my1d_v_front_error: Option<Vec<f64>>,
    /// hxx L274: my2DMaxError.
    my2d_max_error: Option<Vec<f64>>,
    /// hxx L275: my2DAverageError.
    my2d_average_error: Option<Vec<f64>>,
    /// hxx L276: my2DUFrontError.
    my2d_u_front_error: Option<Vec<f64>>,
    /// hxx L277: my2DVFrontError.
    my2d_v_front_error: Option<Vec<f64>>,
    /// hxx L278: my3DMaxError.
    my3d_max_error: Option<Vec<f64>>,
    /// hxx L279: my3DAverageError.
    my3d_average_error: Option<Vec<f64>>,
    /// hxx L280: my3DUFrontError.
    my3d_u_front_error: Option<Vec<f64>>,
    /// hxx L281: my3DVFrontError.
    my3d_v_front_error: Option<Vec<f64>>,
    /// hxx L282: myCriterionError.
    my_criterion_error: f64,
}

impl ApproxAFunc2Var {
    /// OCCT ctor without criterion (AdvApp2Var_ApproxAFunc2Var.cxx L47-101).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        num1dss: i32,
        num2dss: i32,
        num3dss: i32,
        one_dtol: &[f64],
        two_dtol: &[f64],
        three_dtol: &[f64],
        one_dtolf_r: &Array2<f64>,
        two_dtolf_r: &Array2<f64>,
        three_dtolf_r: &Array2<f64>,
        first_in_u: f64,
        last_in_u: f64,
        first_in_v: f64,
        last_in_v: f64,
        favor_iso: IsoType,
        cont_in_u: GeomAbsShape,
        cont_in_v: GeomAbsShape,
        precis_code: i32,
        max_deg_in_u: i32,
        max_deg_in_v: i32,
        max_patch: i32,
        func: &dyn EvaluatorFunc2Var,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
    ) -> Self {
        // Member init-list (L71-92); the handle members copy the inputs.
        let mut r = Self::base_members(
            num1dss,
            num2dss,
            num3dss,
            one_dtol,
            two_dtol,
            three_dtol,
            one_dtolf_r,
            two_dtolf_r,
            three_dtolf_r,
            first_in_u,
            last_in_u,
            first_in_v,
            last_in_v,
            favor_iso,
            cont_in_u,
            cont_in_v,
            precis_code,
            max_deg_in_u,
            max_deg_in_v,
            max_patch,
        );

        r.init();
        r.perform(u_choice, v_choice, func);
        r.convert_bs();
        r
    }

    /// OCCT ctor with criterion (AdvApp2Var_ApproxAFunc2Var.cxx L105-160).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_criterion(
        num1dss: i32,
        num2dss: i32,
        num3dss: i32,
        one_dtol: &[f64],
        two_dtol: &[f64],
        three_dtol: &[f64],
        one_dtolf_r: &Array2<f64>,
        two_dtolf_r: &Array2<f64>,
        three_dtolf_r: &Array2<f64>,
        first_in_u: f64,
        last_in_u: f64,
        first_in_v: f64,
        last_in_v: f64,
        favor_iso: IsoType,
        cont_in_u: GeomAbsShape,
        cont_in_v: GeomAbsShape,
        precis_code: i32,
        max_deg_in_u: i32,
        max_deg_in_v: i32,
        max_patch: i32,
        func: &dyn EvaluatorFunc2Var,
        crit: &dyn Criterion,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
    ) -> Self {
        let mut r = Self::base_members(
            num1dss,
            num2dss,
            num3dss,
            one_dtol,
            two_dtol,
            three_dtol,
            one_dtolf_r,
            two_dtolf_r,
            three_dtolf_r,
            first_in_u,
            last_in_u,
            first_in_v,
            last_in_v,
            favor_iso,
            cont_in_u,
            cont_in_v,
            precis_code,
            max_deg_in_u,
            max_deg_in_v,
            max_patch,
        );

        r.init();
        r.perform_with_criterion(u_choice, v_choice, func, crit);
        r.convert_bs();
        r
    }

    /// The shared member init-list of both OCCT ctors (cxx L71-92 /
    /// L130-151 plus myNumSubSpaces L94-96 / L153-155).
    #[allow(clippy::too_many_arguments)]
    fn base_members(
        num1dss: i32,
        num2dss: i32,
        num3dss: i32,
        one_dtol: &[f64],
        two_dtol: &[f64],
        three_dtol: &[f64],
        one_dtolf_r: &Array2<f64>,
        two_dtolf_r: &Array2<f64>,
        three_dtolf_r: &Array2<f64>,
        first_in_u: f64,
        last_in_u: f64,
        first_in_v: f64,
        last_in_v: f64,
        favor_iso: IsoType,
        cont_in_u: GeomAbsShape,
        cont_in_v: GeomAbsShape,
        precis_code: i32,
        max_deg_in_u: i32,
        max_deg_in_v: i32,
        max_patch: i32,
    ) -> Self {
        ApproxAFunc2Var {
            my_num_sub_spaces: [num1dss, num2dss, num3dss],
            my1d_tolerances: one_dtol.to_vec(),
            my2d_tolerances: two_dtol.to_vec(),
            my3d_tolerances: three_dtol.to_vec(),
            my1d_tol_on_front: one_dtolf_r.clone(),
            my2d_tol_on_front: two_dtolf_r.clone(),
            my3d_tol_on_front: three_dtolf_r.clone(),
            my_first_par_in_u: first_in_u,
            my_last_par_in_u: last_in_u,
            my_first_par_in_v: first_in_v,
            my_last_par_in_v: last_in_v,
            my_favorite_iso: favor_iso,
            my_cont_in_u: cont_in_u,
            my_cont_in_v: cont_in_v,
            my_precision_code: precis_code,
            my_max_deg_in_u: max_deg_in_u,
            my_max_deg_in_v: max_deg_in_v,
            my_max_patches: max_patch,
            my_conditions: Context::new(),
            my_result: Network::new(),
            my_constraints: Framework::new(),
            my_done: false,
            my_has_result: false,
            my_surfaces: None,
            my_degree_in_u: 0,
            my_degree_in_v: 0,
            my1d_max_error: None,
            my1d_average_error: None,
            my1d_u_front_error: None,
            my1d_v_front_error: None,
            my2d_max_error: None,
            my2d_average_error: None,
            my2d_u_front_error: None,
            my2d_v_front_error: None,
            my3d_max_error: None,
            my3d_average_error: None,
            my3d_u_front_error: None,
            my3d_v_front_error: None,
            my_criterion_error: 0.0,
        }
    }

    /// OCCT Init() (AdvApp2Var_ApproxAFunc2Var.cxx L164-235) -
    /// initialisation of the approximation; used by Create.  The dead
    /// initializers of iu/iv (= 0, overwritten by every switch arm) are
    /// OCCT-faithful (cxx L166).
    #[allow(unused_assignments)]
    fn init(&mut self) {
        let ifav: i32;
        let mut iu: i32 = 0;
        let mut iv: i32 = 0;
        let ndu: i32;
        let ndv: i32;
        // switch (myFavoriteIso): the OCCT default arm (ifav = 2) collapses
        // with IsoV - the rcad kernel IsoType only has IsoU/IsoV.
        ifav = match self.my_favorite_iso {
            IsoType::IsoU => 1,
            IsoType::IsoV => 2,
        };
        match self.my_cont_in_u {
            GeomAbsShape::C0 => iu = 0,
            GeomAbsShape::C1 => iu = 1,
            GeomAbsShape::C2 => iu = 2,
            _ => panic!("AdvApp2Var_ApproxAFunc2Var : UContinuity Error"),
        }
        match self.my_cont_in_v {
            GeomAbsShape::C0 => iv = 0,
            GeomAbsShape::C1 => iv = 1,
            GeomAbsShape::C2 => iv = 2,
            _ => panic!("AdvApp2Var_ApproxAFunc2Var : VContinuity Error"),
        }
        ndu = (self.my_max_deg_in_u + 1).max(2 * iu + 2);
        ndv = (self.my_max_deg_in_v + 1).max(2 * iv + 2);
        if ndu < 2 * iu + 2 {
            panic!("AdvApp2Var_ApproxAFunc2Var : UMaxDegree Error");
        }
        if ndv < 2 * iv + 2 {
            panic!("AdvApp2Var_ApproxAFunc2Var : VMaxDegree Error");
        }
        self.my_precision_code = self.my_precision_code.clamp(0, 3);
        let conditions = Context::new_with_params(
            ifav,
            iu,
            iv,
            ndu,
            ndv,
            self.my_precision_code,
            self.my_num_sub_spaces[0],
            self.my_num_sub_spaces[1],
            self.my_num_sub_spaces[2],
            &self.my1d_tolerances,
            &self.my2d_tolerances,
            &self.my3d_tolerances,
            &self.my1d_tol_on_front,
            &self.my2d_tol_on_front,
            &self.my3d_tol_on_front,
        );
        self.my_conditions = conditions;
        self.init_grid(1);
    }

    /// OCCT InitGrid(NbInt) (AdvApp2Var_ApproxAFunc2Var.cxx L239-332) -
    /// initialisation of the approximation with a grid of regular cuttings;
    /// used by Init and Perform.
    fn init_grid(&mut self, nb_int: i32) {
        let iu = self.my_conditions.u_order();
        let iv = self.my_conditions.v_order();

        // occ::handle<AdvApp2Var_Patch> M0 = new AdvApp2Var_Patch(...).
        let m0 = Patch::new_with_domain(
            self.my_first_par_in_u,
            self.my_last_par_in_u,
            self.my_first_par_in_v,
            self.my_last_par_in_v,
            iu,
            iv,
        );

        let mut net: Vec<Patch> = Vec::new();
        net.push(m0);

        let mut the_u: Vec<f64> = Vec::new();
        let mut the_v: Vec<f64> = Vec::new();
        the_u.push(self.my_first_par_in_u);
        the_v.push(self.my_first_par_in_v);
        the_u.push(self.my_last_par_in_u);
        the_v.push(self.my_last_par_in_v);

        let mut result = Network::new_from(net, the_u, the_v);

        let uv1 = DVec2::new(self.my_first_par_in_u, self.my_first_par_in_v);
        let c1 = Node::new_with_coord(uv1, iu, iv);
        let uv2 = DVec2::new(self.my_last_par_in_u, self.my_first_par_in_v);
        let c2 = Node::new_with_coord(uv2, iu, iv);
        let uv4 = DVec2::new(self.my_last_par_in_u, self.my_last_par_in_v);
        let c4 = Node::new_with_coord(uv4, iu, iv);
        let uv3 = DVec2::new(self.my_first_par_in_u, self.my_last_par_in_v);
        let c3 = Node::new_with_coord(uv3, iu, iv);
        let mut bag: Vec<Node> = Vec::new();
        bag.push(c1);
        bag.push(c2);
        bag.push(c3);
        bag.push(c4);

        // OCCT: new AdvApp2Var_Iso(GeomAbs_IsoV, cte, Ufirst, Ulast, Vfirst,
        // Vlast, pos, iu, iv) (AdvApp2Var_ApproxAFunc2Var.cxx L271-288).
        let v0 = super::iso::Iso::new_full(
            IsoType::IsoV,
            self.my_first_par_in_v,
            self.my_first_par_in_u,
            self.my_last_par_in_u,
            self.my_first_par_in_v,
            self.my_last_par_in_v,
            1,
            iu,
            iv,
        );
        let v1 = super::iso::Iso::new_full(
            IsoType::IsoV,
            self.my_last_par_in_v,
            self.my_first_par_in_u,
            self.my_last_par_in_u,
            self.my_first_par_in_v,
            self.my_last_par_in_v,
            2,
            iu,
            iv,
        );
        let u0 = super::iso::Iso::new_full(
            IsoType::IsoU,
            self.my_first_par_in_u,
            self.my_first_par_in_u,
            self.my_last_par_in_u,
            self.my_first_par_in_v,
            self.my_last_par_in_v,
            3,
            iu,
            iv,
        );
        let u1 = super::iso::Iso::new_full(
            IsoType::IsoU,
            self.my_last_par_in_u,
            self.my_first_par_in_u,
            self.my_last_par_in_u,
            self.my_first_par_in_v,
            self.my_last_par_in_v,
            4,
            iu,
            iv,
        );

        let bu0: Vec<super::iso::Iso> = vec![v0, v1];
        let bv0: Vec<super::iso::Iso> = vec![u0, u1];

        let u_strip: Vec<Vec<super::iso::Iso>> = vec![bu0];
        let v_strip: Vec<Vec<super::iso::Iso>> = vec![bv0];

        let mut constraints = Framework::new_from(bag, u_strip, v_strip);

        // regular cutting if NbInt>1
        let deltu = (self.my_last_par_in_u - self.my_first_par_in_u) / nb_int as f64;
        let deltv = (self.my_last_par_in_v - self.my_first_par_in_v) / nb_int as f64;
        let mut iint: i32 = 1;
        while iint <= nb_int - 1 {
            result.update_in_u(self.my_first_par_in_u + iint as f64 * deltu);
            constraints.update_in_u(self.my_first_par_in_u + iint as f64 * deltu);
            result.update_in_v(self.my_first_par_in_v + iint as f64 * deltv);
            constraints.update_in_v(self.my_first_par_in_v + iint as f64 * deltv);
            iint += 1;
        }
        self.my_result = result;
        self.my_constraints = constraints;
    }

    /// OCCT Perform(UChoice, VChoice, Func) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L336-343) - computation of the approximation result; used by Create.
    fn perform(
        &mut self,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
        func: &dyn EvaluatorFunc2Var,
    ) {
        self.compute_patches(u_choice, v_choice, func);
        self.my_has_result = true;
        self.my_done = true;
        self.compute3d_errors();
    }

    /// OCCT Perform(UChoice, VChoice, Func, Crit)
    /// (AdvApp2Var_ApproxAFunc2Var.cxx L347-356).
    fn perform_with_criterion(
        &mut self,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
        func: &dyn EvaluatorFunc2Var,
        crit: &dyn Criterion,
    ) {
        self.compute_patches_with_criterion(u_choice, v_choice, func, crit);
        self.my_has_result = true;
        self.my_done = true;
        self.compute3d_errors();
        self.compute_crit_error();
    }

    /// OCCT ComputePatches(UChoice, VChoice, Func)
    /// (AdvApp2Var_ApproxAFunc2Var.cxx L360-475).
    fn compute_patches(
        &mut self,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
        func: &dyn EvaluatorFunc2Var,
    ) {
        // OCCT declares these at function scope (uninitialized) and
        // re-assigns them per loop iteration; the rcad bindings are mut for
        // the same reason.  first_na = 0 is the rcad requirement for the
        // by-address out parameter.
        let mut udec: f64;
        let mut vdec: f64;
        let mut umore: bool;
        let mut vmore: bool;
        let mut nb_patch: i32;
        let mut nb_u: i32;
        let mut nb_v: i32;
        let mut num_dec: i32;
        let mut first_na: i32 = 0;

        // OCCT: while (myResult.FirstNotApprox(FirstNA)).
        while self.my_result.first_not_approx(&mut first_na) {
            // complete the set of constraints
            self.compute_constraints(u_choice, v_choice, func);

            // discretization of constraints relative to the square
            self.my_result.change_patch(first_na).discretise(
                &self.my_conditions,
                &self.my_constraints,
                func,
            );
            if !self.my_result.change_patch(first_na).is_discretised() {
                self.my_has_result = false;
                self.my_done = false;
                panic!("AdvApp2Var_ApproxAFunc2Var : Surface Discretisation Error");
            }

            // calculate the number and the type of authorized cuts
            // depending on the max number of squares and the validity of
            // next cuts.
            nb_u = self.my_result.nb_patch_in_u();
            nb_v = self.my_result.nb_patch_in_v();
            nb_patch = nb_u * nb_v;
            // OCCT: Umore = UChoice.Value(myResult(FirstNA).U0(),
            // myResult(FirstNA).U1(), Udec) - the rcad Cutting collapses
            // (more, dec) into Option<f64>; the patch bounds are read
            // before the call (OCCT evaluates them as arguments).
            let (pa_u0, pa_u1, pa_v0, pa_v1) = {
                let a_patch = self.my_result.change_patch(first_na);
                (a_patch.u0(), a_patch.u1(), a_patch.v0(), a_patch.v1())
            };
            let u_choice_val = u_choice.value(pa_u0, pa_u1);
            umore = u_choice_val.is_some();
            udec = u_choice_val.unwrap_or(0.0);
            let v_choice_val = v_choice.value(pa_v0, pa_v1);
            vmore = v_choice_val.is_some();
            vdec = v_choice_val.unwrap_or(0.0);

            num_dec = 0;
            if ((nb_patch + nb_v) <= self.my_max_patches)
                && ((nb_patch + nb_u) > self.my_max_patches)
                && umore
            {
                num_dec = 1;
            }
            if ((nb_patch + nb_v) > self.my_max_patches)
                && ((nb_patch + nb_u) <= self.my_max_patches)
                && vmore
            {
                num_dec = 2;
            }
            if ((nb_patch + nb_v) <= self.my_max_patches)
                && ((nb_patch + nb_u) <= self.my_max_patches)
            {
                if umore {
                    num_dec = 3;
                }
                if (nb_v > nb_u) && vmore {
                    num_dec = 4;
                }
            }
            if (nb_u + 1) * (nb_v + 1) <= self.my_max_patches {
                if !umore && !vmore {
                    num_dec = 0;
                }
                if umore && !vmore {
                    num_dec = 3;
                }
                if !umore && vmore {
                    num_dec = 4;
                }
                if umore && vmore {
                    num_dec = 5;
                }
            }

            // approximation of the square
            self.my_result
                .change_patch(first_na)
                .make_approx(&self.my_conditions, &self.my_constraints, num_dec);

            if !self.my_result.change_patch(first_na).is_approximated() {
                match self.my_result.change_patch(first_na).cut_sense() {
                    0 => {
                        //	It is not possible to cut : the result is preserved
                        if self.my_result.change_patch(first_na).has_result() {
                            self.my_result.change_patch(first_na).overwrite_approx();
                        } else {
                            self.my_has_result = false;
                            self.my_done = false;
                            panic!("AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error");
                        }
                    }
                    1 => {
                        //      It is necessary to cut in U
                        self.my_result.update_in_u(udec);
                        self.my_constraints.update_in_u(udec);
                    }
                    2 => {
                        //      It is necessary to cut in V
                        self.my_result.update_in_v(vdec);
                        self.my_constraints.update_in_v(vdec);
                    }
                    3 => {
                        //      It is necessary to cut in U and V
                        self.my_result.update_in_u(udec);
                        self.my_constraints.update_in_u(udec);
                        self.my_result.update_in_v(vdec);
                        self.my_constraints.update_in_v(vdec);
                    }
                    _ => {
                        self.my_has_result = false;
                        self.my_done = false;
                        panic!("AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error");
                    }
                }
            }
        }
    }

    /// OCCT ComputePatches(UChoice, VChoice, Func, Crit)
    /// (AdvApp2Var_ApproxAFunc2Var.cxx L479-630).
    fn compute_patches_with_criterion(
        &mut self,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
        func: &dyn EvaluatorFunc2Var,
        crit: &dyn Criterion,
    ) {
        // OCCT declares these at function scope (uninitialized) and
        // re-assigns them per loop iteration; the rcad bindings are mut for
        // the same reason.  m1 = 0. / decision = 0 are OCCT-faithful
        // initializers (cxx L484, L487); first_na = 0 is the rcad
        // requirement for the by-address out parameter.
        let mut udec: f64;
        let mut vdec: f64;
        let mut crit_value: f64;
        let mut m1: f64 = 0.;
        let mut umore: bool;
        let mut vmore: bool;
        let crit_abs = crit.crit_type() == CriterionType::Absolute;
        let mut nb_patch: i32;
        let mut nb_u: i32;
        let mut nb_v: i32;
        let mut nb_int: i32;
        let mut num_dec: i32;
        let mut first_na: i32 = 0;
        let mut decision: i32 = 0;

        while self.my_result.first_not_approx(&mut first_na) {
            // complete the set of constraints
            self.compute_constraints_with_criterion(u_choice, v_choice, func, crit);
            if decision > 0 {
                m1 = 0.;
            }

            // discretize the constraints relative to the square
            self.my_result.change_patch(first_na).discretise(
                &self.my_conditions,
                &self.my_constraints,
                func,
            );
            if !self.my_result.change_patch(first_na).is_discretised() {
                self.my_has_result = false;
                self.my_done = false;
                panic!("AdvApp2Var_ApproxAFunc2Var : Surface Discretisation Error");
            }

            // calculate the number and type of authorized cuts
            // depending on the max number of squares and the validity of
            // next cuts
            nb_u = self.my_result.nb_patch_in_u();
            nb_v = self.my_result.nb_patch_in_v();
            nb_patch = nb_u * nb_v;
            nb_int = nb_u;
            // OCCT: Umore = UChoice.Value(..., Udec); Vmore = VChoice.Value(
            //   ..., Vdec) - the rcad Cutting collapses (more, dec) into
            // Option<f64>.
            let (pa_u0, pa_u1, pa_v0, pa_v1) = {
                let a_patch = self.my_result.change_patch(first_na);
                (a_patch.u0(), a_patch.u1(), a_patch.v0(), a_patch.v1())
            };
            let u_choice_val = u_choice.value(pa_u0, pa_u1);
            umore = u_choice_val.is_some();
            udec = u_choice_val.unwrap_or(0.0);
            let v_choice_val = v_choice.value(pa_v0, pa_v1);
            vmore = v_choice_val.is_some();
            vdec = v_choice_val.unwrap_or(0.0);

            num_dec = 0;
            if ((nb_patch + nb_v) <= self.my_max_patches)
                && ((nb_patch + nb_u) > self.my_max_patches)
                && umore
            {
                num_dec = 1;
            }
            if ((nb_patch + nb_v) > self.my_max_patches)
                && ((nb_patch + nb_u) <= self.my_max_patches)
                && vmore
            {
                num_dec = 2;
            }
            if ((nb_patch + nb_v) <= self.my_max_patches)
                && ((nb_patch + nb_u) <= self.my_max_patches)
            {
                if umore {
                    num_dec = 3;
                }
                if (nb_v > nb_u) && vmore {
                    num_dec = 4;
                }
            }
            if (nb_u + 1) * (nb_v + 1) <= self.my_max_patches {
                if !umore && !vmore {
                    num_dec = 0;
                }
                if umore && !vmore {
                    num_dec = 1;
                }
                if !umore && vmore {
                    num_dec = 2;
                }
                if umore && vmore {
                    num_dec = 5;
                }
            }

            // approximation of the square
            if crit_abs {
                self.my_result
                    .change_patch(first_na)
                    .make_approx(&self.my_conditions, &self.my_constraints, 0);
            } else {
                self.my_result
                    .change_patch(first_na)
                    .make_approx(&self.my_conditions, &self.my_constraints, num_dec);
            }
            if num_dec >= 3 {
                num_dec -= 2;
            }

            // evaluation of the criterion on the square
            if self.my_result.change_patch(first_na).has_result() {
                // OCCT: Crit.Value(myResult(FirstNA), myConditions) - the
                // criterion mutates the patch through the handle; the
                // rcad trait mutates the patch slot (disjoint field
                // borrows of self).
                crit.value(
                    self.my_result.change_patch(first_na),
                    &self.my_conditions,
                );
                crit_value = self.my_result.change_patch(first_na).crit_value();
                if m1 < crit_value {
                    m1 = crit_value;
                }
            }
            // is it necessary to cut ?
            decision = self
                .my_result
                .change_patch(first_na)
                .cut_sense_criterion(crit, num_dec);
            let regular = crit.repartition() == CriterionRepartition::Regular;
            //    bool Regular = true;
            if regular && decision > 0 {
                nb_int += 1;
                self.init_grid(nb_int);
            } else {
                match decision {
                    0 => {
                        //	Impossible to cut : the result is preserved
                        if self.my_result.change_patch(first_na).has_result() {
                            self.my_result.change_patch(first_na).overwrite_approx();
                        } else {
                            self.my_has_result = false;
                            self.my_done = false;
                            panic!("AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error");
                        }
                    }
                    1 => {
                        //      It is necessary to cut in U
                        self.my_result.update_in_u(udec);
                        self.my_constraints.update_in_u(udec);
                    }
                    2 => {
                        //      It is necessary to cut in V
                        self.my_result.update_in_v(vdec);
                        self.my_constraints.update_in_v(vdec);
                    }
                    3 => {
                        //      It is necessary to cut in U and V
                        self.my_result.update_in_u(udec);
                        self.my_constraints.update_in_u(udec);
                        self.my_result.update_in_v(vdec);
                        self.my_constraints.update_in_v(vdec);
                    }
                    _ => {
                        self.my_has_result = false;
                        self.my_done = false;
                        panic!("AdvApp2Var_ApproxAFunc2Var : Surface Approximation Error");
                    }
                }
            }
        }
    }

    /// OCCT ComputeConstraints(UChoice, VChoice, Func)
    /// (AdvApp2Var_ApproxAFunc2Var.cxx L634-719).  The N1/N2 ctor values
    /// are overwritten before use in the loop (OCCT-faithful).
    #[allow(unused_assignments)]
    fn compute_constraints(
        &mut self,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
        func: &dyn EvaluatorFunc2Var,
    ) {
        // OCCT declares these at function scope (uninitialized) and
        // re-assigns them per loop iteration; the rcad bindings are mut for
        // the same reason.
        let mut dec: f64;
        let mut more: bool;
        let mut ind1: i32 = 0;
        let mut ind2: i32 = 0;
        let mut nb_patch: i32;
        let mut nb_u: i32;
        let mut nb_v: i32;
        let iu = self.my_conditions.u_order();
        let iv = self.my_conditions.v_order();
        let mut n1 = Node::new_orders(iu, iv);
        let mut n2 = Node::new_orders(iu, iv);

        // OCCT: for (anIso = myConstraints.FirstNotApprox(ind1, ind2);
        //   !anIso.IsNull(); anIso = myConstraints.FirstNotApprox(ind1, ind2))
        loop {
            let an_iso = match self.my_constraints.first_not_approx(&mut ind1, &mut ind2) {
                Some(an_iso) => an_iso,
                None => break, // anIso.IsNull()
            };

            // approximation of iso and calculation of constraints at
            // extremities
            let ind_n1 = self
                .my_constraints
                .first_node(an_iso.type_(), ind1, ind2);
            n1 = self.my_constraints.node_index(ind_n1).clone();
            let ind_n2 = self
                .my_constraints
                .last_node(an_iso.type_(), ind1, ind2);
            n2 = self.my_constraints.node_index(ind_n2).clone();

            // note that old code attempted to make copy of anIso here (but
            // copy was incomplete)
            let mut an_iso = an_iso;
            an_iso.make_approx(
                &self.my_conditions,
                self.my_first_par_in_u,
                self.my_last_par_in_u,
                self.my_first_par_in_v,
                self.my_last_par_in_v,
                func,
                &mut n1,
                &mut n2,
            );
            if an_iso.is_approximated() {
                // iso is approached at the required tolerance
                self.my_constraints.change_iso(ind1, ind2, an_iso);
                *self.my_constraints.node_index_mut(ind_n1) = n1.clone();
                *self.my_constraints.node_index_mut(ind_n2) = n2.clone();
            } else {
                // Approximation is not satisfactory
                nb_u = self.my_result.nb_patch_in_u();
                nb_v = self.my_result.nb_patch_in_v();
                // OCCT: more = UChoice.Value(anIso->T0(), anIso->T1(), dec)
                // (or VChoice) - the rcad Cutting collapses (more, dec)
                // into Option<f64>.
                let a_t0 = an_iso.t0();
                let a_t1 = an_iso.t1();
                if an_iso.type_() == IsoType::IsoV {
                    nb_patch = (nb_u + 1) * nb_v;
                    let v = u_choice.value(a_t0, a_t1);
                    more = v.is_some();
                    dec = v.unwrap_or(0.0);
                } else {
                    nb_patch = (nb_v + 1) * nb_u;
                    let v = v_choice.value(a_t0, a_t1);
                    more = v.is_some();
                    dec = v.unwrap_or(0.0);
                }

                if nb_patch <= self.my_max_patches && more {
                    // It is possible to cut iso
                    if an_iso.type_() == IsoType::IsoV {
                        self.my_result.update_in_u(dec);
                        self.my_constraints.update_in_u(dec);
                    } else {
                        self.my_result.update_in_v(dec);
                        self.my_constraints.update_in_v(dec);
                    }
                } else {
                    // It is not possible to cut : the result is preserved
                    if an_iso.has_result() {
                        let mut an_iso = an_iso;
                        an_iso.overwrite_approx();
                        self.my_constraints.change_iso(ind1, ind2, an_iso);
                        *self.my_constraints.node_index_mut(ind_n1) = n1.clone();
                        *self.my_constraints.node_index_mut(ind_n2) = n2.clone();
                    } else {
                        self.my_has_result = false;
                        self.my_done = false;
                        panic!("AdvApp2Var_ApproxAFunc2Var : Curve Approximation Error");
                    }
                }
            }
        }
    }

    /// OCCT ComputeConstraints(UChoice, VChoice, Func, Crit)
    /// (AdvApp2Var_ApproxAFunc2Var.cxx L723-814).  The N1/N2 ctor values
    /// are overwritten before use in the loop (OCCT-faithful).
    #[allow(unused_assignments)]
    fn compute_constraints_with_criterion(
        &mut self,
        u_choice: &dyn Cutting,
        v_choice: &dyn Cutting,
        func: &dyn EvaluatorFunc2Var,
        crit: &dyn Criterion,
    ) {
        // OCCT declares these at function scope (uninitialized) and
        // re-assigns them per loop iteration; the rcad bindings are mut for
        // the same reason.
        let mut dec: f64;
        let mut more: bool;
        let crit_rel = crit.crit_type() == CriterionType::Relative;
        let mut ind1: i32 = 0;
        let mut ind2: i32 = 0;
        let mut nb_patch: i32;
        let mut nb_u: i32;
        let mut nb_v: i32;
        let mut ind_n1: i32;
        let mut ind_n2: i32;
        let iu = self.my_conditions.u_order();
        let iv = self.my_conditions.v_order();
        let mut n1 = Node::new_orders(iu, iv);
        let mut n2 = Node::new_orders(iu, iv);

        loop {
            let an_iso = match self.my_constraints.first_not_approx(&mut ind1, &mut ind2) {
                Some(an_iso) => an_iso,
                None => break, // anIso.IsNull()
            };

            // approximation of the iso and calculation of constraints at the
            // extremities
            ind_n1 = self
                .my_constraints
                .first_node(an_iso.type_(), ind1, ind2);
            n1 = self.my_constraints.node_index(ind_n1).clone();
            ind_n2 = self
                .my_constraints
                .last_node(an_iso.type_(), ind1, ind2);
            n2 = self.my_constraints.node_index(ind_n2).clone();

            // note that old code attempted to make copy of anIso here (but
            // copy was incomplete)
            let mut an_iso = an_iso;
            an_iso.make_approx(
                &self.my_conditions,
                self.my_first_par_in_u,
                self.my_last_par_in_u,
                self.my_first_par_in_v,
                self.my_last_par_in_v,
                func,
                &mut n1,
                &mut n2,
            );

            if an_iso.is_approximated() {
                // iso is approached at the required tolerance
                self.my_constraints.change_iso(ind1, ind2, an_iso);
                *self.my_constraints.node_index_mut(ind_n1) = n1.clone();
                *self.my_constraints.node_index_mut(ind_n2) = n2.clone();
            } else {
                // Approximation is not satisfactory
                nb_u = self.my_result.nb_patch_in_u();
                nb_v = self.my_result.nb_patch_in_v();
                // OCCT: more = UChoice.Value(anIso->T0(), anIso->T1(), dec)
                // (or VChoice) - the rcad Cutting collapses (more, dec)
                // into Option<f64>.
                // OCCT: more = UChoice.Value(anIso->T0(), anIso->T1(), dec)
                // (or VChoice) - the rcad Cutting collapses (more, dec)
                // into Option<f64>.
                let a_t0 = an_iso.t0();
                let a_t1 = an_iso.t1();
                if an_iso.type_() == IsoType::IsoV {
                    nb_patch = (nb_u + 1) * nb_v;
                    let v = u_choice.value(a_t0, a_t1);
                    more = v.is_some();
                    dec = v.unwrap_or(0.0);
                } else {
                    nb_patch = (nb_v + 1) * nb_u;
                    let v = v_choice.value(a_t0, a_t1);
                    more = v.is_some();
                    dec = v.unwrap_or(0.0);
                }

                // To force Overwrite if the criterion is Absolute
                more = more && crit_rel;

                if nb_patch <= self.my_max_patches && more {
                    // It is possible to cut iso
                    if an_iso.type_() == IsoType::IsoV {
                        self.my_result.update_in_u(dec);
                        self.my_constraints.update_in_u(dec);
                    } else {
                        self.my_result.update_in_v(dec);
                        self.my_constraints.update_in_v(dec);
                    }
                } else {
                    // It is not possible to cut: the result is preserved
                    if an_iso.has_result() {
                        let mut an_iso = an_iso;
                        an_iso.overwrite_approx();
                        self.my_constraints.change_iso(ind1, ind2, an_iso);
                        *self.my_constraints.node_index_mut(ind_n1) = n1.clone();
                        *self.my_constraints.node_index_mut(ind_n2) = n2.clone();
                    } else {
                        self.my_has_result = false;
                        self.my_done = false;
                        panic!("AdvApp2Var_ApproxAFunc2Var : Curve Approximation Error");
                    }
                }
            }
        }
    }

    /// OCCT Compute3DErrors() (AdvApp2Var_ApproxAFunc2Var.cxx L818-864).
    fn compute3d_errors(&mut self) {
        let mut iesp: i32;
        let mut ipat: i32;
        let mut error_max: f64;
        let mut error_moy: f64;
        let mut error_u0: f64;
        let mut error_v0: f64;
        let mut error_u1: f64;
        let mut error_v1: f64;
        let mut tol: f64;
        let mut f1tol: f64;
        let mut f2tol: f64;
        let mut f3tol: f64;
        let mut f4tol: f64;
        if self.my_num_sub_spaces[2] > 0 {
            // new (NCollection_HArray1<double>)(1, myNumSubSpaces[2]) x4.
            self.my3d_max_error = Some(vec![0.0; self.my_num_sub_spaces[2] as usize]);
            self.my3d_average_error = Some(vec![0.0; self.my_num_sub_spaces[2] as usize]);
            self.my3d_u_front_error = Some(vec![0.0; self.my_num_sub_spaces[2] as usize]);
            self.my3d_v_front_error = Some(vec![0.0; self.my_num_sub_spaces[2] as usize]);
            iesp = 1;
            while iesp <= self.my_num_sub_spaces[2] {
                error_max = 0.0;
                error_moy = 0.;
                error_u0 = 0.;
                error_v0 = 0.;
                error_u1 = 0.;
                error_v1 = 0.;
                tol = self.my3d_tolerances[(iesp - 1) as usize];
                f1tol = self.my3d_tol_on_front.value(iesp, 1);
                f2tol = self.my3d_tol_on_front.value(iesp, 2);
                f3tol = self.my3d_tol_on_front.value(iesp, 3);
                f4tol = self.my3d_tol_on_front.value(iesp, 4);
                ipat = 1;
                while ipat <= self.my_result.nb_patch() {
                    // OCCT: (myResult(ipat).MaxErrors())->Value(iesp) - the
                    // HArray1 flat entry (iesp).
                    let a_patch = self.my_result.patch(ipat, 1);
                    error_max = a_patch
                        .max_errors()
                        .expect("Compute3DErrors : MaxErrors")[(iesp - 1) as usize]
                        .max(error_max);
                    error_u0 = a_patch
                        .iso_errors()
                        .expect("Compute3DErrors : IsoErrors")
                        .value(iesp, 3)
                        .max(error_u0);
                    error_u1 = a_patch
                        .iso_errors()
                        .expect("Compute3DErrors : IsoErrors")
                        .value(iesp, 4)
                        .max(error_u1);
                    error_v0 = a_patch
                        .iso_errors()
                        .expect("Compute3DErrors : IsoErrors")
                        .value(iesp, 1)
                        .max(error_v0);
                    error_v1 = a_patch
                        .iso_errors()
                        .expect("Compute3DErrors : IsoErrors")
                        .value(iesp, 2)
                        .max(error_v1);
                    error_moy += a_patch
                        .average_errors()
                        .expect("Compute3DErrors : AverageErrors")[(iesp - 1) as usize];
                    ipat += 1;
                }
                let my3d_max_error = self.my3d_max_error.as_mut().expect("Compute3DErrors");
                my3d_max_error[(iesp - 1) as usize] = error_max;
                let my3d_u_front_error = self.my3d_u_front_error.as_mut().expect("Compute3DErrors");
                my3d_u_front_error[(iesp - 1) as usize] = error_u0.max(error_u1);
                let my3d_v_front_error = self.my3d_v_front_error.as_mut().expect("Compute3DErrors");
                my3d_v_front_error[(iesp - 1) as usize] = error_v0.max(error_v1);
                error_moy /= self.my_result.nb_patch() as f64;
                let my3d_average_error = self.my3d_average_error.as_mut().expect("Compute3DErrors");
                my3d_average_error[(iesp - 1) as usize] = error_moy;
                if error_max > tol || error_u0 > f3tol || error_u1 > f4tol || error_v0 > f1tol
                    || error_v1 > f2tol
                {
                    self.my_done = false;
                }
                iesp += 1;
            }
        }
    }

    /// OCCT ComputeCritError() (AdvApp2Var_ApproxAFunc2Var.cxx L868-885).
    fn compute_crit_error(&mut self) {
        let mut iesp: i32;
        let mut ipat: i32;
        let mut crit_max: f64;
        if self.my_num_sub_spaces[2] > 0 {
            iesp = 1;
            while iesp <= self.my_num_sub_spaces[2] {
                crit_max = 0.;
                ipat = 1;
                while ipat <= self.my_result.nb_patch() {
                    crit_max = self.my_result.patch(ipat, 1).crit_value().max(crit_max);
                    ipat += 1;
                }
                self.my_criterion_error = crit_max;
                iesp += 1;
            }
        }
    }

    /// OCCT ConvertBS() (AdvApp2Var_ApproxAFunc2Var.cxx L889-990) -
    /// conversion of the approximation result in BSpline; used by Create.
    fn convert_bs(&mut self) {
        // Homogeneization of degrees
        let iu = self.my_conditions.u_order();
        let iv = self.my_conditions.v_order();
        let mut ncfu = self.my_conditions.u_limit();
        let mut ncfv = self.my_conditions.v_limit();
        self.my_result.same_degree(iu, iv, &mut ncfu, &mut ncfv);
        self.my_degree_in_u = ncfu - 1;
        self.my_degree_in_v = ncfv - 1;

        // Calculate resulting surfaces
        // mySurfaces = new (HArray1<handle<Geom_Surface>>)(1,
        //   myNumSubSpaces[2]) - entries start as null handles.
        self.my_surfaces = Some(vec![None; self.my_num_sub_spaces[2] as usize]);

        let mut j: i32;
        // NCollection_Array1<double> UKnots(1, NbPatchInU() + 1).
        let mut u_knots = vec![0.0f64; (self.my_result.nb_patch_in_u() + 1) as usize];
        j = 1;
        while j <= u_knots.len() as i32 {
            u_knots[(j - 1) as usize] = self.my_result.u_parameter(j);
            j += 1;
        }

        let mut v_knots = vec![0.0f64; (self.my_result.nb_patch_in_v() + 1) as usize];
        j = 1;
        while j <= v_knots.len() as i32 {
            v_knots[(j - 1) as usize] = self.my_result.v_parameter(j);
            j += 1;
        }

        // Prepare data for conversion grid of polynoms --> poles
        // Uint1 = ...(1, 2) = (-1, 1); Vint1 = ...(1, 2) = (-1, 1).
        let u_int1 = vec![-1.0f64, 1.0];
        let v_int1 = vec![-1.0f64, 1.0];

        let mut u_int2 = vec![0.0f64; (self.my_result.nb_patch_in_u() + 1) as usize];
        j = 1;
        while j <= u_int2.len() as i32 {
            u_int2[(j - 1) as usize] = self.my_result.u_parameter(j);
            j += 1;
        }
        let mut v_int2 = vec![0.0f64; (self.my_result.nb_patch_in_v() + 1) as usize];
        j = 1;
        while j <= v_int2.len() as i32 {
            v_int2[(j - 1) as usize] = self.my_result.v_parameter(j);
            j += 1;
        }

        let nmax = (self.my_result.nb_patch_in_u() * self.my_result.nb_patch_in_v()) as usize;
        let size_eq =
            (self.my_conditions.u_limit() * self.my_conditions.v_limit() * 3) as usize;

        // NbCoeff = new (HArray2<int>)(1, nmax, 1, 2); Poly = new
        //   (HArray1<double>)(1, nmax * Size_eq).
        let mut nb_coeff = Array2::new(1, nmax as i32, 1, 2);
        let mut poly = vec![0.0f64; nmax * size_eq];

        let mut ssp: i32;
        let mut i: i32;
        ssp = 1;
        while ssp <= self.my_num_sub_spaces[2] {
            // Creation of the grid of polynoms
            let mut n: i32 = 0;
            let mut icf: i32 = 1;
            let mut ieq: i32;
            j = 1;
            while j <= self.my_result.nb_patch_in_v() {
                i = 1;
                while i <= self.my_result.nb_patch_in_u() {
                    n += 1;
                    nb_coeff.set_value(n, 1, self.my_result.patch(i, j).nb_coeff_in_u());
                    nb_coeff.set_value(n, 2, self.my_result.patch(i, j).nb_coeff_in_v());
                    ieq = 1;
                    while ieq <= size_eq as i32 {
                        poly[(icf - 1) as usize] = self
                            .my_result
                            .patch(i, j)
                            .coefficients(ssp, &self.my_conditions)[(ieq - 1) as usize];
                        icf += 1;
                        ieq += 1;
                    }
                    i += 1;
                }
                j += 1;
            }

            // Conversion into poles.
            // GAP: Convert_GridPolynomialToPoles (TKMath/Convert, the
            // 12-argument grid constructor, AdvApp2Var_ApproxAFunc2Var.cxx
            // L964-975: Convert_GridPolynomialToPoles
            // CvP(NbPatchInU(), NbPatchInV(), iu, iv, myMaxDegInU,
            // myMaxDegInV, NbCoeff->Array2(), Poly->Array1(), Uint1, Vint1,
            // Uint2, Vint2)) is not yet translated - it is owned by the
            // convert_comp_polynomial_to_poles.rs batch.  The OCCT failure
            // branch `if (!CvP.IsDone()) myDone = false;` is preserved: the
            // un-translated conversion yields IsDone() == false, myDone is
            // cleared and the surface entry keeps its null handle (the
            // OCCT `new Geom_BSplineSurface(CvP...)` anchor L982-988 never
            // executes with a valid conversion).  The data prepared above
            // (NbCoeff, Poly, Uint1, Vint1, Uint2, Vint2) is exactly the
            // ctor argument set of the un-translated converter.
            let _ = (
                &nb_coeff,
                &poly,
                &u_int1,
                &v_int1,
                &u_int2,
                &v_int2,
            );
            self.my_done = false;

            ssp += 1;
        }
    }

    /// OCCT IsDone() (AdvApp2Var_ApproxAFunc2Var.lxx L20-23) - true if the
    /// approximation succeeded within the imposed tolerances and the wished
    /// continuities.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT HasResult() (AdvApp2Var_ApproxAFunc2Var.lxx L25-28).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT Surface(SSPIndex) (AdvApp2Var_ApproxAFunc2Var.lxx L30-34) - the
    /// occ::down_cast<Geom_BSplineSurface> becomes the inner Option; returns
    /// the BSplineSurface of range Index.
    pub fn surface(&self, ssp_index: i32) -> Option<&BSplineSurface> {
        self.my_surfaces
            .as_ref()
            .and_then(|a| a[(ssp_index - 1) as usize].as_ref())
    }

    /// OCCT UDegree() (AdvApp2Var_ApproxAFunc2Var.lxx L36-39).
    pub fn u_degree(&self) -> i32 {
        self.my_degree_in_u
    }

    /// OCCT VDegree() (AdvApp2Var_ApproxAFunc2Var.lxx L41-44).
    pub fn v_degree(&self) -> i32 {
        self.my_degree_in_v
    }

    /// OCCT NumSubSpaces(Dimension) (AdvApp2Var_ApproxAFunc2Var.lxx L46-49).
    pub fn num_sub_spaces(&self, dimension: i32) -> i32 {
        self.my_num_sub_spaces[(dimension - 1) as usize]
    }

    /// OCCT MaxError(Dimension) (AdvApp2Var_ApproxAFunc2Var.cxx L994-1016) -
    /// returns the errors max.
    pub fn max_error(&self, dimension: i32) -> Option<Vec<f64>> {
        if !(1..=3).contains(&dimension) {
            panic!("AdvApp2Var_ApproxAFunc2Var::MaxError : Dimension must be equal to 1,2 or 3 !");
        }
        match dimension {
            1 => self.my1d_max_error.clone(),
            2 => self.my2d_max_error.clone(),
            3 => self.my3d_max_error.clone(),
            _ => unreachable!(),
        }
    }

    /// OCCT AverageError(Dimension) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1020-1042) - returns the average errors.
    pub fn average_error(&self, dimension: i32) -> Option<Vec<f64>> {
        if !(1..=3).contains(&dimension) {
            panic!("AdvApp2Var_ApproxAFunc2Var::AverageError : Dimension must be equal to 1,2 or 3 !");
        }
        match dimension {
            1 => self.my1d_average_error.clone(),
            2 => self.my2d_average_error.clone(),
            3 => self.my3d_average_error.clone(),
            _ => unreachable!(),
        }
    }

    /// OCCT UFrontError(Dimension) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1046-1068) - returns the errors max on UFrontiers.
    pub fn u_front_error(&self, dimension: i32) -> Option<Vec<f64>> {
        if !(1..=3).contains(&dimension) {
            panic!("AdvApp2Var_ApproxAFunc2Var::UFrontError : Dimension must be equal to 1,2 or 3 !");
        }
        match dimension {
            1 => self.my1d_u_front_error.clone(),
            2 => self.my2d_u_front_error.clone(),
            3 => self.my3d_u_front_error.clone(),
            _ => unreachable!(),
        }
    }

    /// OCCT VFrontError(Dimension) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1072-1094) - returns the errors max on VFrontiers.
    pub fn v_front_error(&self, dimension: i32) -> Option<Vec<f64>> {
        if dimension <= 0 || dimension > 3 {
            panic!("AdvApp2Var_ApproxAFunc2Var::VFrontError : Dimension must be equal to 1,2 or 3 !");
        }
        match dimension {
            1 => self.my1d_v_front_error.clone(),
            2 => self.my2d_v_front_error.clone(),
            3 => self.my3d_v_front_error.clone(),
            _ => unreachable!(),
        }
    }

    /// OCCT MaxError(Dimension, SSPIndex) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1098-1106) - returns the error max of the BSplineSurface of range
    /// SSPIndex.
    pub fn max_error_at(&self, dimension: i32, ssp_index: i32) -> f64 {
        if dimension != 3 || ssp_index != 1 {
            panic!("AdvApp2Var_ApproxAFunc2Var::MaxError: ONE Surface 3D only !");
        }
        let e_ptr = self.max_error(dimension).expect("MaxError");
        e_ptr[(ssp_index - 1) as usize]
    }

    /// OCCT AverageError(Dimension, SSPIndex) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1110-1118).
    pub fn average_error_at(&self, dimension: i32, ssp_index: i32) -> f64 {
        if dimension != 3 || ssp_index != 1 {
            panic!("AdvApp2Var_ApproxAFunc2Var::AverageError : ONE Surface 3D only !");
        }
        let e_ptr = self.average_error(dimension).expect("AverageError");
        e_ptr[(ssp_index - 1) as usize]
    }

    /// OCCT UFrontError(Dimension, SSPIndex) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1122-1130).
    pub fn u_front_error_at(&self, dimension: i32, ssp_index: i32) -> f64 {
        if dimension != 3 || ssp_index != 1 {
            panic!("AdvApp2Var_ApproxAFunc2Var::UFrontError : ONE Surface 3D only !");
        }
        let e_ptr = self.u_front_error(dimension).expect("UFrontError");
        e_ptr[(ssp_index - 1) as usize]
    }

    /// OCCT VFrontError(Dimension, SSPIndex) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1134-1142).
    pub fn v_front_error_at(&self, dimension: i32, ssp_index: i32) -> f64 {
        if dimension != 3 || ssp_index != 1 {
            panic!("AdvApp2Var_ApproxAFunc2Var::VFrontError : ONE Surface 3D only !");
        }
        let e_ptr = self.v_front_error(dimension).expect("VFrontError");
        e_ptr[(ssp_index - 1) as usize]
    }

    /// OCCT CritError(Dimension, SSPIndex) (AdvApp2Var_ApproxAFunc2Var.cxx
    /// L1146-1153).
    pub fn crit_error(&self, dimension: i32, ssp_index: i32) -> f64 {
        if dimension != 3 || ssp_index != 1 {
            panic!("AdvApp2Var_ApproxAFunc2Var::CritError: ONE Surface 3D only !");
        }
        self.my_criterion_error
    }

    /// OCCT Dump(o) (AdvApp2Var_ApproxAFunc2Var.cxx L1157-1207) - prints on
    /// the stream 'o' information on the current state of the object; the
    /// rcad dump convention returns the formatted String.
    pub fn dump(&self) -> String {
        let iesp: i32 = 1;
        let mut o = String::new();
        let _ = o.write_char('\n');
        if !self.my_has_result {
            let _ = o.write_str("No result\n");
        } else {
            o.write_str("There is a result").ok();
            if self.my_done {
                let _ = write!(
                    o,
                    " within the requested tolerance {}\n",
                    self.my3d_tolerances[(iesp - 1) as usize]
                );
            } else if self.my3d_max_error.as_ref().expect("Dump")[(iesp - 1) as usize]
                > self.my3d_tolerances[(iesp - 1) as usize]
            {
                let _ = write!(
                    o,
                    " WITHOUT the requested tolerance {}\n",
                    self.my3d_tolerances[(iesp - 1) as usize]
                );
            } else {
                let _ = o.write_str(" WITHOUT the requested continuities \n");
            }
            let _ = o.write_char('\n');
            let _ = writeln!(
                o,
                "Result max error :{}",
                self.my3d_max_error.as_ref().expect("Dump")[(iesp - 1) as usize]
            );
            let _ = writeln!(
                o,
                "Result average error :{}",
                self.my3d_average_error.as_ref().expect("Dump")[(iesp - 1) as usize]
            );
            let _ = writeln!(
                o,
                "Result max error on U frontiers :{}",
                self.my3d_u_front_error.as_ref().expect("Dump")[(iesp - 1) as usize]
            );
            let _ = writeln!(
                o,
                "Result max error on V frontiers :{}",
                self.my3d_v_front_error.as_ref().expect("Dump")[(iesp - 1) as usize]
            );
            let _ = o.write_char('\n');
            let _ = writeln!(
                o,
                "Degree of Bezier patches in U : {}  in V : {}",
                self.my_degree_in_u, self.my_degree_in_v
            );
            let _ = o.write_char('\n');
            let s = self.surface(iesp).expect("Dump : surface");
            // Architecture note: the rcad kernel BSplineSurface stores the
            // EXPANDED knot vectors (no separate multiplicities), so the
            // OCCT NbUPoles/NbUKnots/UKnot(ik)/UMultiplicity(ik) printout
            // degrades to pole counts and the flat knot values.
            let _ = writeln!(
                o,
                "Number of poles in U : {}  in V : {}",
                s.control_points.len(),
                s.control_points.first().map(|r| r.len()).unwrap_or(0)
            );
            let _ = o.write_char('\n');
            let nb_ku = s.knots_u.len();
            let nb_kv = s.knots_v.len();
            let _ = writeln!(o, "Number of knots in U : {}", nb_ku);
            for ik in 0..nb_ku {
                let _ = writeln!(o, "   {} : {}", ik + 1, s.knots_u[ik]);
            }
            let _ = o.write_char('\n');
            let _ = writeln!(o, "Number of knots in V : {}", nb_kv);
            for ik in 0..nb_kv {
                let _ = writeln!(o, "   {} : {}", ik + 1, s.knots_v[ik]);
            }
            let _ = o.write_char('\n');
        }
        o
    }
}
