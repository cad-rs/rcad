// OCCT BRepApprox PrmPrm chain — the parametric-surface x parametric-surface
// intersection engine feeding the intersection-line approximation.
//
// 1:1 translation of the template family instantiated by the BRepApprox
// alias tables (TKTopAlgo/BRepApprox/*_0.cxx):
//   IntImp_ZerParFunc.gxx (L15-382) + .lxx (L15-61)
//     -> [`ZerParFunc`]  (= BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox)
//   IntImp_Int2S.gxx (L15-268) + .lxx (L15-94)
//     -> [`Int2S`]       (= BRepApprox_TheInt2SOfThePrmPrmSvSurfacesOfApprox)
//   ApproxInt_PrmPrmSvSurfaces.gxx (L17-373)
//     -> [`PrmPrmSvSurfaces`] (= BRepApprox_ThePrmPrmSvSurfacesOfApprox)
//   IntImp_ComputeTangence.cxx (L19-168) + IntImp_ConstIsoparametric.hxx
//     -> [`choix_ref`] / [`compute_tangence`] / [`ConstIsoparametric`]
//
// Alias-table bindings (the three _0.cxx files):
//   ThePSurface     = BRepAdaptor_Surface    -> BRepAdaptorSurface<'a>
//   ThePSurfaceTool = BRepApprox_SurfaceTool -> brep_approx::surface_tool
//   TheLine         = BRepApprox_ApproxLine  -> brep_approx::ApproxLine
//                     (declared by the alias table; none of the three gxx
//                      bodies reference TheLine, so nothing binds it)
//   IntImp_ZerParFunc / IntImp_Int2S / ApproxInt_PrmPrmSvSurfaces ->
//   the corresponding BRepApprox_The...OfApprox names.
//
// The OCCT `void* surf1/surf2` raw-pointer members (the borrowed surfaces)
// map to shared references with an explicit lifetime; the non-const
// ApproxInt_SvSurfaces interface (called through `&self` receivers in the
// rcad trait encoding, like the OCCT non-const methods through the void*)
// maps to Cell/RefCell interior mutability.

use std::cell::{Cell, RefCell};

use glam::{DVec2, DVec3};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use rcad_kernel::{ANGULAR, CONFUSION};

use crate::geomalgo::brep_approx::{surface_tool, SvSurfaces};
use crate::geomalgo::int_surf::PntOn2S;
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

// ---------------------------------------------------------------------------
// IntImp_ConstIsoparametric + IntImp_ComputeTangence (package-level free
// members of OCCT IntImp; shared with the future ImpPrm chain)
// ---------------------------------------------------------------------------

/// OCCT IntImp_ConstIsoparametric (IntImp_ConstIsoparametric.hxx L18-23) —
/// the constant isoparametric kept fixed while the other three parameters
/// are solved (IntImp_UIsoparametricOnCaro1 / V...OnCaro1 / U...OnCaro2 /
/// V...OnCaro2 in OCCT order).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstIsoparametric {
    UIsoparametricOnCaro1 = 0,
    VIsoparametricOnCaro1 = 1,
    UIsoparametricOnCaro2 = 2,
    VIsoparametricOnCaro2 = 3,
}

// OCCT IntImp_ComputeTangence.cxx L19-24: staticChoixRef.
static STATIC_CHOIX_REF: [ConstIsoparametric; 4] = [
    ConstIsoparametric::UIsoparametricOnCaro1,
    ConstIsoparametric::VIsoparametricOnCaro1,
    ConstIsoparametric::UIsoparametricOnCaro2,
    ConstIsoparametric::VIsoparametricOnCaro2,
];

/// OCCT ChoixRef (IntImp_ComputeTangence.cxx L26-30) — the reference choice
/// by index; Standard_OutOfRange below 0 / above 3 (the negative index
/// cannot occur at the call sites, which pass ints >= 0).
pub fn choix_ref(the_index: usize) -> ConstIsoparametric {
    // OCCT Standard_OutOfRange_Raise_if(theIndex < 0 || theIndex > 3, ...).
    if the_index > 3 {
        panic!("Standard_OutOfRange: ChoixRef()");
    }
    STATIC_CHOIX_REF[the_index]
}

/// OCCT IntImp_ComputeTangence (IntImp_ComputeTangence.cxx L34-168) —
/// computes the tangent of the intersection in the two tangent planes and
/// reports whether the surfaces are tangent there.
/// Inputs (cxx L38-46): DPuv[0..1] = dP/du, dP/dv on caro 1, DPuv[2..3] on
/// caro 2; EpsUV[0..3] = u/v tolerances caro 1 then caro 2.
/// Outputs (cxx L47-53): Tgduv[0..1] = the tangent components in the caro-1
/// tangent plane, Tgduv[2..3] = caro 2; TabIso[0..3] = best iso choices.
// The gxx files using this routine define No_Standard_RangeError /
// No_Standard_OutOfRange; no bounds checks are added.
pub fn compute_tangence(
    dp_uv: &[DVec3; 4],
    eps_uv: &[f64; 4],
    tgduv: &mut [f64; 4],
    tab_iso: &mut [ConstIsoparametric; 4],
) -> bool {
    // OCCT L68-71.
    let mut norm_duv = [0.0f64; 4];
    let mut a_m2;
    let a_tol2 = 1.0e-32;

    // OCCT L73-80.
    for i in 0..4 {
        norm_duv[i] = dp_uv[i].length_squared();
        if norm_duv[i] <= a_tol2 {
            return true;
        }
    }

    // OCCT L83-84.
    let mut n1 = dp_uv[0];
    n1 = n1.cross(dp_uv[1]);

    // OCCT L87-92 (NIZNHY-PKV 2011).
    a_m2 = n1.length_squared();
    if a_m2 < a_tol2 {
        return true;
    }
    // OCCT L93.
    n1 = n1.normalize();

    // OCCT L95-96.
    let mut n2 = dp_uv[2];
    n2 = n2.cross(dp_uv[3]);

    // OCCT L98-102 (NIZNHY-PKV 2011).
    a_m2 = n2.length_squared();
    if a_m2 < a_tol2 {
        return true;
    }
    // OCCT L104.
    n2 = n2.normalize();

    // OCCT L107-110 (NIZNHY-PKV 2011).
    for i in 0..4 {
        norm_duv[i] = norm_duv[i].sqrt();
    }

    // OCCT L112-115.
    tgduv[0] = -dp_uv[1].dot(n2);
    tgduv[1] = dp_uv[0].dot(n2);
    tgduv[2] = dp_uv[3].dot(n1);
    tgduv[3] = -dp_uv[2].dot(n1);

    // OCCT L117-119.
    let mut tangent = tgduv[0].abs() <= eps_uv[0] * norm_duv[1]
        && tgduv[1].abs() <= eps_uv[1] * norm_duv[0]
        && tgduv[2].abs() <= eps_uv[2] * norm_duv[3]
        && tgduv[3].abs() <= eps_uv[3] * norm_duv[2];

    // OCCT L120-131.
    if !tangent {
        let mut t = n1.dot(n2);
        if t < 0.0 {
            t = -t;
        }
        if t > 0.999999999 {
            tangent = true;
        }
    }

    // OCCT L133-166: rank the isoparametric candidates (bubble sort of the
    // normalized tangent components, in parallel with ChoixRef).
    if !tangent {
        norm_duv[0] = tgduv[1].abs() / norm_duv[0]; // iso u sur caro1
        norm_duv[1] = tgduv[0].abs() / norm_duv[1]; // iso v sur caro1
        norm_duv[2] = tgduv[3].abs() / norm_duv[2]; // iso u sur caro2
        norm_duv[3] = tgduv[2].abs() / norm_duv[3]; // iso v sur caro2

        let mut tri_ok;
        for i in 0..=3 {
            tab_iso[i] = STATIC_CHOIX_REF[i];
        }
        loop {
            tri_ok = true;
            for i in 1..=3 {
                if norm_duv[i - 1] > norm_duv[i] {
                    tri_ok = false;
                    let t = norm_duv[i];
                    norm_duv[i] = norm_duv[i - 1];
                    norm_duv[i - 1] = t;

                    let ti = tab_iso[i];
                    tab_iso[i] = tab_iso[i - 1];
                    tab_iso[i - 1] = ti;
                }
            }
            if tri_ok {
                break;
            }
        }
    }
    tangent
}

// OCCT ApproxInt_PrmPrmSvSurfaces.gxx L17: #define TOLTANGENCY 0.0000000001.
const TOLTANGENCY: f64 = 0.0000000001;

// ---------------------------------------------------------------------------
// ZerParFunc — BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox
// ---------------------------------------------------------------------------

/// OCCT IntImp_ZerParFunc (IntImp_ZerParFunc.gxx) instantiated as
/// BRepApprox_TheFunctionOfTheInt2SOfThePrmPrmSvSurfacesOfApprox
/// (_0.cxx L25-31) — the function set F = S1(..) - S2(..) whose zero is a
/// surface/surface intersection point, with one parameter constrained to
/// an isoparametric of one of the surfaces.  The OCCT `void* surf1/surf2`
/// members map to shared references (the borrowed BRepAdaptor_Surface
/// objects must outlive this function).
pub struct ZerParFunc<'a> {
    // OCCT hxx private section (the generated GeomInt_TheFunctionOfTheInt2S...hxx
    // has the same members): void* surf1; void* surf2;
    surf1: &'a BRepAdaptorSurface<'a>,
    surf2: &'a BRepAdaptorSurface<'a>,
    // OCCT gp_Pnt pntsol1; gp_Pnt pntsol2;
    pntsol1: DVec3,
    pntsol2: DVec3,
    // OCCT double f[3];
    f: [f64; 3],
    // OCCT bool compute; — set by the constructor, never read (vestigial
    // in OCCT as well); kept 1:1.
    #[allow(dead_code)]
    compute: bool,
    // OCCT bool tangent;
    tangent: bool,
    // OCCT double tgduv[4]; — left uninitialized by OCCT; zeroed here.
    tgduv: [f64; 4],
    // OCCT gp_Vec dpuv[4]; — default-constructed (zero) vectors in OCCT.
    dpuv: [DVec3; 4],
    // OCCT IntImp_ConstIsoparametric chxIso; — left uninitialized by the
    // OCCT constructor (set by ComputeParameters); the first enum value
    // stands in for it.
    chx_iso: ConstIsoparametric,
    // OCCT double paramConst;
    param_const: f64,
    // OCCT double ua0, va0, ua1, va1, ub0, vb0, ub1, vb1, ures1, ures2,
    // vres1, vres2;
    ua0: f64,
    va0: f64,
    ua1: f64,
    va1: f64,
    ub0: f64,
    vb0: f64,
    ub1: f64,
    vb1: f64,
    ures1: f64,
    ures2: f64,
    vres1: f64,
    vres2: f64,
}

impl<'a> ZerParFunc<'a> {
    /// OCCT IntImp_ZerParFunc(S1, S2) — gxx L27-53.
    pub fn new(s1: &'a BRepAdaptorSurface<'a>, s2: &'a BRepAdaptorSurface<'a>) -> Self {
        // OCCT L32-33: surf1 = (void*)(&S1); surf2 = (void*)(&S2);
        let surf1 = s1;
        let surf2 = s2;

        // OCCT L35-43.
        let ua0 = surface_tool::first_u_parameter(surf1);
        let va0 = surface_tool::first_v_parameter(surf1);
        let ua1 = surface_tool::last_u_parameter(surf1);
        let va1 = surface_tool::last_v_parameter(surf1);

        let ub0 = surface_tool::first_u_parameter(surf2);
        let vb0 = surface_tool::first_v_parameter(surf2);
        let ub1 = surface_tool::last_u_parameter(surf2);
        let vb1 = surface_tool::last_v_parameter(surf2);

        // OCCT L45-49.
        let ures1 = surface_tool::u_resolution(surf1, CONFUSION);
        let vres1 = surface_tool::v_resolution(surf1, CONFUSION);

        let ures2 = surface_tool::u_resolution(surf2, CONFUSION);
        let vres2 = surface_tool::v_resolution(surf2, CONFUSION);

        ZerParFunc {
            surf1,
            surf2,
            pntsol1: DVec3::ZERO, // OCCT gp_Pnt default (0, 0, 0)
            pntsol2: DVec3::ZERO, // OCCT gp_Pnt default (0, 0, 0)
            f: [0.0; 3],          // OCCT L50: memset(f, 0, sizeof(f));
            compute: false,       // OCCT L51: compute = false;
            tangent: false,       // OCCT L52: tangent = false;
            tgduv: [0.0; 4],      // OCCT leaves tgduv[4] uninitialized
            dpuv: [DVec3::ZERO; 4], // OCCT gp_Vec default (0, 0, 0)
            chx_iso: ConstIsoparametric::UIsoparametricOnCaro1,
            param_const: 0.0, // OCCT L30: paramConst(0.0)
            ua0,
            va0,
            ua1,
            va1,
            ub0,
            vb0,
            ub1,
            vb1,
            ures1,
            ures2,
            vres1,
            vres2,
        }
    }

    /// OCCT ComputeParameters — gxx L231-328: stores ChoixIso, distributes
    /// the 4 input parameters over the fixed isoparametric and the 3
    /// unknowns, sets the bounds and the per-variable tolerances.
    /// (Rust mapping: `param` = OCCT const NCollection_Array1 Param (1..4);
    /// `uvap`/`born_inf`/`born_sup`/`tolerance` = the math_Vector (1..3).)
    pub fn compute_parameters(
        &mut self,
        choix_iso: ConstIsoparametric,
        param: &[f64; 4],
        uvap: &mut [f64; 3],
        born_inf: &mut [f64; 3],
        born_sup: &mut [f64; 3],
        tolerance: &mut [f64; 3],
    ) {
        // OCCT L239.
        self.chx_iso = choix_iso;
        // OCCT L240-317: switch (chxIso).  The statement order inside each
        // case is preserved verbatim (case 1 fills BornInf(2)/BornInf(3)
        // before BornSup(2)/BornSup(3)).
        match self.chx_iso {
            ConstIsoparametric::UIsoparametricOnCaro1 => {
                self.param_const = param[0];
                uvap[0] = param[1];
                uvap[1] = param[2];
                uvap[2] = param[3];

                born_inf[0] = self.va0;
                born_sup[0] = self.va1;

                born_inf[1] = self.ub0;
                born_inf[2] = self.vb0;
                born_sup[1] = self.ub1;
                born_sup[2] = self.vb1;

                tolerance[0] = self.vres1;
                tolerance[1] = self.ures2;
                tolerance[2] = self.vres2;
            }
            ConstIsoparametric::VIsoparametricOnCaro1 => {
                self.param_const = param[1];
                uvap[0] = param[0];
                uvap[1] = param[2];
                uvap[2] = param[3];
                born_inf[0] = self.ua0;
                born_sup[0] = self.ua1;

                born_inf[1] = self.ub0;
                born_sup[1] = self.ub1;
                born_inf[2] = self.vb0;
                born_sup[2] = self.vb1;

                tolerance[0] = self.ures1;
                tolerance[1] = self.ures2;
                tolerance[2] = self.vres2;
            }
            ConstIsoparametric::UIsoparametricOnCaro2 => {
                self.param_const = param[2];
                uvap[0] = param[0];
                uvap[1] = param[1];
                uvap[2] = param[3];

                born_inf[0] = self.ua0;
                born_sup[0] = self.ua1;
                born_inf[1] = self.va0;
                born_sup[1] = self.va1;

                born_inf[2] = self.vb0;
                born_sup[2] = self.vb1;

                tolerance[0] = self.ures1;
                tolerance[1] = self.vres1;
                tolerance[2] = self.vres2;
            }
            ConstIsoparametric::VIsoparametricOnCaro2 => {
                self.param_const = param[3];
                uvap[0] = param[0];
                uvap[1] = param[1];
                uvap[2] = param[2];

                born_inf[0] = self.ua0;
                born_sup[0] = self.ua1;
                born_inf[1] = self.va0;
                born_sup[1] = self.va1;

                born_inf[2] = self.ub0;
                born_sup[2] = self.ub1;

                tolerance[0] = self.ures1;
                tolerance[1] = self.vres1;
                tolerance[2] = self.ures2;
            }
        }

        // OCCT L319-327.
        let incr1 = (born_sup[0] - born_inf[0]) * 0.01;
        let incr2 = (born_sup[1] - born_inf[1]) * 0.01;
        let incr3 = (born_sup[2] - born_inf[2]) * 0.01;
        born_inf[0] -= incr1;
        born_sup[0] += incr1;
        born_inf[1] -= incr2;
        born_sup[1] += incr2;
        born_inf[2] -= incr3;
        born_sup[2] += incr3;
    }

    /// OCCT IsTangent — gxx L330-379: re-distributes the solved unknowns
    /// over Param, computes the tangent with IntImp_ComputeTangence and, if
    /// not tangent, updates chxIso to the best next choice.
    pub fn is_tangent(
        &mut self,
        uvap: &[f64; 3],
        param: &mut [f64; 4],
        best_choix: &mut ConstIsoparametric,
    ) -> bool {
        // OCCT L334-364: switch (chxIso).
        match self.chx_iso {
            ConstIsoparametric::UIsoparametricOnCaro1 => {
                param[0] = self.param_const;
                param[1] = uvap[0];
                param[2] = uvap[1];
                param[3] = uvap[2];
            }
            ConstIsoparametric::VIsoparametricOnCaro1 => {
                param[1] = self.param_const;
                param[0] = uvap[0];
                param[2] = uvap[1];
                param[3] = uvap[2];
            }
            ConstIsoparametric::UIsoparametricOnCaro2 => {
                param[2] = self.param_const;
                param[0] = uvap[0];
                param[1] = uvap[1];
                param[3] = uvap[2];
            }
            ConstIsoparametric::VIsoparametricOnCaro2 => {
                param[3] = self.param_const;
                param[0] = uvap[0];
                param[1] = uvap[1];
                param[2] = uvap[2];
            }
        }

        // OCCT L366: IntImp_ConstIsoparametric TabIso[4]; — uninitialized in
        // OCCT; the first enum value stands in for it.
        let mut tab_iso = [ConstIsoparametric::UIsoparametricOnCaro1; 4];
        // OCCT L367-372.
        let mut eps_uv = [0.0f64; 4];
        eps_uv[0] = self.ures1;
        eps_uv[1] = self.vres1;

        eps_uv[2] = self.ures2;
        eps_uv[3] = self.vres2;

        // OCCT L374-377.
        self.tangent = compute_tangence(&self.dpuv, &eps_uv, &mut self.tgduv, &mut tab_iso);
        if !self.tangent {
            self.chx_iso = tab_iso[0];
        }
        *best_choix = self.chx_iso;
        self.tangent
    }

    /// OCCT Root() — lxx L21-25: the sum of the squared function values.
    pub fn root(&self) -> f64 {
        self.f[0] * self.f[0] + self.f[1] * self.f[1] + self.f[2] * self.f[2]
    }

    /// OCCT Point() — lxx L27-30: the mid-point of the two surface points.
    pub fn point(&self) -> DVec3 {
        (self.pntsol1 + self.pntsol2) / 2.0
    }

    /// OCCT Direction() — lxx L32-37: the 3d tangent built from the caro-1
    /// components (raises StdFail_UndefinedDerivative when tangent).
    pub fn direction(&self) -> DVec3 {
        if self.tangent {
            panic!("StdFail_UndefinedDerivative");
        }
        (self.dpuv[0] * self.tgduv[0] + self.dpuv[1] * self.tgduv[1]).normalize()
    }

    /// OCCT DirectionOnS1() — lxx L39-44.
    pub fn direction_on_s1(&self) -> DVec2 {
        if self.tangent {
            panic!("StdFail_UndefinedDerivative");
        }
        DVec2::new(self.tgduv[0], self.tgduv[1]).normalize()
    }

    /// OCCT DirectionOnS2() — lxx L46-51.
    pub fn direction_on_s2(&self) -> DVec2 {
        if self.tangent {
            panic!("StdFail_UndefinedDerivative");
        }
        DVec2::new(self.tgduv[2], self.tgduv[3]).normalize()
    }

    /// OCCT AuxillarSurface1() — lxx L53-56.
    pub fn auxillar_surface1(&self) -> &'a BRepAdaptorSurface<'a> {
        self.surf1
    }

    /// OCCT AuxillarSurface2() — lxx L58-61.
    pub fn auxillar_surface2(&self) -> &'a BRepAdaptorSurface<'a> {
        self.surf2
    }
}

// The math_FunctionSetWithDerivatives base class of IntImp_ZerParFunc.
impl FunctionSetWithDerivatives for ZerParFunc<'_> {
    /// OCCT NbVariables() — gxx L55-58.
    fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() — gxx L60-63.
    fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) — gxx L65-96 (X(1..3) -> x[0..2], F(1..3) -> f[0..2]).
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L68-90: switch (chxIso) — the two surface evaluations.
        match self.chx_iso {
            ConstIsoparametric::UIsoparametricOnCaro1 => {
                self.pntsol1 = surface_tool::value(self.surf1, self.param_const, x[0]);
                self.pntsol2 = surface_tool::value(self.surf2, x[1], x[2]);
            }
            ConstIsoparametric::VIsoparametricOnCaro1 => {
                self.pntsol1 = surface_tool::value(self.surf1, x[0], self.param_const);
                self.pntsol2 = surface_tool::value(self.surf2, x[1], x[2]);
            }
            ConstIsoparametric::UIsoparametricOnCaro2 => {
                self.pntsol1 = surface_tool::value(self.surf1, x[0], x[1]);
                self.pntsol2 = surface_tool::value(self.surf2, self.param_const, x[2]);
            }
            ConstIsoparametric::VIsoparametricOnCaro2 => {
                self.pntsol1 = surface_tool::value(self.surf1, x[0], x[1]);
                self.pntsol2 = surface_tool::value(self.surf2, x[2], self.param_const);
            }
        }

        // OCCT L92-95: f[0] = F(1) = pntsol1.X() - pntsol2.X(); ...
        let f1 = self.pntsol1.x - self.pntsol2.x;
        let f2 = self.pntsol1.y - self.pntsol2.y;
        let f3 = self.pntsol1.z - self.pntsol2.z;
        self.f[0] = f1;
        f[0] = f1;
        self.f[1] = f2;
        f[1] = f2;
        self.f[2] = f3;
        f[2] = f3;
        true
    }

    /// OCCT Derivatives(X, D) — gxx L98-161 (D(i, j) -> df[i-1][j-1]).
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        // OCCT L101-159: switch (chxIso).
        match self.chx_iso {
            ConstIsoparametric::UIsoparametricOnCaro1 => {
                // OCCT L104-105.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, self.param_const, x[0]);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, x[1], x[2]);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L106-114.
                df[0][0] = self.dpuv[1].x;
                df[0][1] = -self.dpuv[2].x;
                df[0][2] = -self.dpuv[3].x;
                df[1][0] = self.dpuv[1].y;
                df[1][1] = -self.dpuv[2].y;
                df[1][2] = -self.dpuv[3].y;
                df[2][0] = self.dpuv[1].z;
                df[2][1] = -self.dpuv[2].z;
                df[2][2] = -self.dpuv[3].z;
            }
            ConstIsoparametric::VIsoparametricOnCaro1 => {
                // OCCT L118-119.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, x[0], self.param_const);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, x[1], x[2]);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L120-128.
                df[0][0] = self.dpuv[0].x;
                df[0][1] = -self.dpuv[2].x;
                df[0][2] = -self.dpuv[3].x;
                df[1][0] = self.dpuv[0].y;
                df[1][1] = -self.dpuv[2].y;
                df[1][2] = -self.dpuv[3].y;
                df[2][0] = self.dpuv[0].z;
                df[2][1] = -self.dpuv[2].z;
                df[2][2] = -self.dpuv[3].z;
            }
            ConstIsoparametric::UIsoparametricOnCaro2 => {
                // OCCT L132-133.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, x[0], x[1]);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, self.param_const, x[2]);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L134-142.
                df[0][0] = self.dpuv[0].x;
                df[0][1] = self.dpuv[1].x;
                df[0][2] = -self.dpuv[3].x;
                df[1][0] = self.dpuv[0].y;
                df[1][1] = self.dpuv[1].y;
                df[1][2] = -self.dpuv[3].y;
                df[2][0] = self.dpuv[0].z;
                df[2][1] = self.dpuv[1].z;
                df[2][2] = -self.dpuv[3].z;
            }
            ConstIsoparametric::VIsoparametricOnCaro2 => {
                // OCCT L146-147.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, x[0], x[1]);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, x[2], self.param_const);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L148-156.
                df[0][0] = self.dpuv[0].x;
                df[0][1] = self.dpuv[1].x;
                df[0][2] = -self.dpuv[2].x;
                df[1][0] = self.dpuv[0].y;
                df[1][1] = self.dpuv[1].y;
                df[1][2] = -self.dpuv[2].y;
                df[2][0] = self.dpuv[0].z;
                df[2][1] = self.dpuv[1].z;
                df[2][2] = -self.dpuv[2].z;
            }
        }
        true
    }

    /// OCCT Values(X, F, D) — gxx L163-229: the derivatives first (the same
    /// switch as Derivatives), then the function values.
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        // OCCT L166-224: switch (chxIso).
        match self.chx_iso {
            ConstIsoparametric::UIsoparametricOnCaro1 => {
                // OCCT L169-170.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, self.param_const, x[0]);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, x[1], x[2]);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L171-179.
                df[0][0] = self.dpuv[1].x;
                df[0][1] = -self.dpuv[2].x;
                df[0][2] = -self.dpuv[3].x;
                df[1][0] = self.dpuv[1].y;
                df[1][1] = -self.dpuv[2].y;
                df[1][2] = -self.dpuv[3].y;
                df[2][0] = self.dpuv[1].z;
                df[2][1] = -self.dpuv[2].z;
                df[2][2] = -self.dpuv[3].z;
            }
            ConstIsoparametric::VIsoparametricOnCaro1 => {
                // OCCT L183-184.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, x[0], self.param_const);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, x[1], x[2]);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L185-193.
                df[0][0] = self.dpuv[0].x;
                df[0][1] = -self.dpuv[2].x;
                df[0][2] = -self.dpuv[3].x;
                df[1][0] = self.dpuv[0].y;
                df[1][1] = -self.dpuv[2].y;
                df[1][2] = -self.dpuv[3].y;
                df[2][0] = self.dpuv[0].z;
                df[2][1] = -self.dpuv[2].z;
                df[2][2] = -self.dpuv[3].z;
            }
            ConstIsoparametric::UIsoparametricOnCaro2 => {
                // OCCT L197-198.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, x[0], x[1]);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, self.param_const, x[2]);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L199-207.
                df[0][0] = self.dpuv[0].x;
                df[0][1] = self.dpuv[1].x;
                df[0][2] = -self.dpuv[3].x;
                df[1][0] = self.dpuv[0].y;
                df[1][1] = self.dpuv[1].y;
                df[1][2] = -self.dpuv[3].y;
                df[2][0] = self.dpuv[0].z;
                df[2][1] = self.dpuv[1].z;
                df[2][2] = -self.dpuv[3].z;
            }
            ConstIsoparametric::VIsoparametricOnCaro2 => {
                // OCCT L211-212.
                let (p1, d1u, d1v) = surface_tool::d1(self.surf1, x[0], x[1]);
                self.pntsol1 = p1;
                self.dpuv[0] = d1u;
                self.dpuv[1] = d1v;
                let (p2, d1u, d1v) = surface_tool::d1(self.surf2, x[2], self.param_const);
                self.pntsol2 = p2;
                self.dpuv[2] = d1u;
                self.dpuv[3] = d1v;
                // OCCT L213-221.
                df[0][0] = self.dpuv[0].x;
                df[0][1] = self.dpuv[1].x;
                df[0][2] = -self.dpuv[2].x;
                df[1][0] = self.dpuv[0].y;
                df[1][1] = self.dpuv[1].y;
                df[1][2] = -self.dpuv[2].y;
                df[2][0] = self.dpuv[0].z;
                df[2][1] = self.dpuv[1].z;
                df[2][2] = -self.dpuv[2].z;
            }
        }

        // OCCT L225-228.
        let f1 = self.pntsol1.x - self.pntsol2.x;
        let f2 = self.pntsol1.y - self.pntsol2.y;
        let f3 = self.pntsol1.z - self.pntsol2.z;
        self.f[0] = f1;
        f[0] = f1;
        self.f[1] = f2;
        f[1] = f2;
        self.f[2] = f3;
        f[2] = f3;
        true
    }
}

// ---------------------------------------------------------------------------
// Int2S — BRepApprox_TheInt2SOfThePrmPrmSvSurfacesOfApprox
// ---------------------------------------------------------------------------

/// OCCT IntImp_Int2S (IntImp_Int2S.gxx) instantiated as
/// BRepApprox_TheInt2SOfThePrmPrmSvSurfacesOfApprox (_0.cxx L31-39) — the
/// driver solving ZerParFunc for the intersection point of two surfaces.
pub struct Int2S<'a> {
    // OCCT bool done;
    done: bool,
    // OCCT bool empty;
    empty: bool,
    // OCCT IntSurf_PntOn2S pint;
    pint: PntOn2S,
    // OCCT bool tangent;
    tangent: bool,
    // OCCT gp_Dir d3d; gp_Dir2d d2d1; gp_Dir2d d2d2;
    d3d: DVec3,
    d2d1: DVec2,
    d2d2: DVec2,
    // OCCT IntImp_TheFunction myZerParFunc;
    my_zer_par_func: ZerParFunc<'a>,
    // OCCT double tol;
    tol: f64,
    // OCCT double ua0, va0, ua1, va1, ub0, vb0, ub1, vb1, ures1, ures2,
    // vres1, vres2;  (the four resolution members are write-only in the
    // OCCT class as well — Perform recomputes Epsuv through the tool.)
    ua0: f64,
    va0: f64,
    ua1: f64,
    va1: f64,
    ub0: f64,
    vb0: f64,
    ub1: f64,
    vb1: f64,
    #[allow(dead_code)]
    ures1: f64,
    #[allow(dead_code)]
    ures2: f64,
    #[allow(dead_code)]
    vres1: f64,
    #[allow(dead_code)]
    vres2: f64,
}

impl<'a> Int2S<'a> {
    /// OCCT IntImp_Int2S(surf1, surf2, TolTangency) — gxx L27-51: initialize
    /// the parameters to compute the solution point.
    pub fn new(
        surf1: &'a BRepAdaptorSurface<'a>,
        surf2: &'a BRepAdaptorSurface<'a>,
        tol_tangency: f64,
    ) -> Self {
        // OCCT L36-44.
        let ua0 = surface_tool::first_u_parameter(surf1);
        let va0 = surface_tool::first_v_parameter(surf1);
        let ua1 = surface_tool::last_u_parameter(surf1);
        let va1 = surface_tool::last_v_parameter(surf1);

        let ub0 = surface_tool::first_u_parameter(surf2);
        let vb0 = surface_tool::first_v_parameter(surf2);
        let ub1 = surface_tool::last_u_parameter(surf2);
        let vb1 = surface_tool::last_v_parameter(surf2);

        // OCCT L46-50.
        let ures1 = surface_tool::u_resolution(surf1, CONFUSION);
        let vres1 = surface_tool::v_resolution(surf1, CONFUSION);

        let ures2 = surface_tool::u_resolution(surf2, CONFUSION);
        let vres2 = surface_tool::v_resolution(surf2, CONFUSION);

        Int2S {
            done: true,           // OCCT L30
            empty: true,          // OCCT L31
            pint: PntOn2S::new(), // OCCT IntSurf_PntOn2S default
            tangent: false,       // OCCT L32
            d3d: DVec3::ZERO,     // OCCT gp_Dir default
            d2d1: DVec2::ZERO,    // OCCT gp_Dir2d default
            d2d2: DVec2::ZERO,    // OCCT gp_Dir2d default
            my_zer_par_func: ZerParFunc::new(surf1, surf2), // OCCT L33
            tol: tol_tangency * tol_tangency, // OCCT L34
            ua0,
            va0,
            ua1,
            va1,
            ub0,
            vb0,
            ub1,
            vb1,
            ures1,
            ures2,
            vres1,
            vres2,
        }
    }

    /// OCCT IntImp_Int2S(Param, surf1, surf2, TolTangency) — gxx L53-79:
    /// compute the solution point with the close point.
    pub fn new_with_param(
        param: &[f64; 4],
        surf1: &'a BRepAdaptorSurface<'a>,
        surf2: &'a BRepAdaptorSurface<'a>,
        tol_tangency: f64,
    ) -> Self {
        // OCCT L57-60: done(true), empty(true), myZerParFunc(surf1, surf2),
        // tol(TolTangency * TolTangency).  The OCCT init list leaves
        // `tangent` uninitialized; the false default stands in for it.
        let my_zer_par_func = ZerParFunc::new(surf1, surf2);
        // OCCT L62: math_FunctionSetRoot Rsnld(myZerParFunc, 15);
        //          //-- Modif lbr 18 MAI ?????????????
        // (the 2-argument OCCT constructor leaves Tol zero-initialized).
        let mut rsnld = FunctionSetRoot::new(&my_zer_par_func, &[0.0, 0.0, 0.0], 15);

        // OCCT L63-71.
        let ua0 = surface_tool::first_u_parameter(surf1);
        let va0 = surface_tool::first_v_parameter(surf1);
        let ua1 = surface_tool::last_u_parameter(surf1);
        let va1 = surface_tool::last_v_parameter(surf1);

        let ub0 = surface_tool::first_u_parameter(surf2);
        let vb0 = surface_tool::first_v_parameter(surf2);
        let ub1 = surface_tool::last_u_parameter(surf2);
        let vb1 = surface_tool::last_v_parameter(surf2);

        // OCCT L73-77.
        let ures1 = surface_tool::u_resolution(surf1, CONFUSION);
        let vres1 = surface_tool::v_resolution(surf1, CONFUSION);

        let ures2 = surface_tool::u_resolution(surf2, CONFUSION);
        let vres2 = surface_tool::v_resolution(surf2, CONFUSION);

        let mut int2s = Int2S {
            done: true,
            empty: true,
            pint: PntOn2S::new(),
            tangent: false,
            d3d: DVec3::ZERO,
            d2d1: DVec2::ZERO,
            d2d2: DVec2::ZERO,
            my_zer_par_func,
            tol: tol_tangency * tol_tangency,
            ua0,
            va0,
            ua1,
            va1,
            ub0,
            vb0,
            ub1,
            vb1,
            ures1,
            ures2,
            vres1,
            vres2,
        };
        // OCCT L78: Perform(Param, Rsnld);
        int2s.perform(param, &mut rsnld);
        int2s
    }

    /// OCCT Perform(Param, Rsnld, ChoixIso) — gxx L81-123: solve with the
    /// isoparametric choice given; stores the solution point.
    #[allow(unused_assignments)]
    pub fn perform_with_choice(
        &mut self,
        param: &[f64; 4],
        rsnld: &mut FunctionSetRoot,
        choix_iso: ConstIsoparametric,
    ) -> ConstIsoparametric {
        // OCCT L85-89: the stack buffers BornInfBuf[3] ... UVapBuf[3] and
        // UvresBuf[4] with their 1-based views (math_Vector / NCollection).
        let mut born_inf = [0.0f64; 3];
        let mut born_sup = [0.0f64; 3];
        let mut tolerance = [0.0f64; 3];
        let mut uvap = [0.0f64; 3];
        let mut uvres = [0.0f64; 4];

        // OCCT L91: IntImp_ConstIsoparametric BestChoix; — uninitialized in
        // OCCT until L96 assigns it.
        let mut best_choix = choix_iso;

        // OCCT L93.
        self.my_zer_par_func.compute_parameters(
            choix_iso,
            param,
            &mut uvap,
            &mut born_inf,
            &mut born_sup,
            &mut tolerance,
        );
        // OCCT L94-95.
        rsnld.set_tolerance(&tolerance);
        rsnld.perform(&mut self.my_zer_par_func, &uvap, &born_inf, &born_sup, false);
        // OCCT L96.
        best_choix = choix_iso;
        if rsnld.is_done() {
            // OCCT L99: if (std::abs(myZerParFunc.Root()) <= tol) — distance
            // des 2 points dans la tolerance.
            if self.my_zer_par_func.root().abs() <= self.tol {
                // OCCT L102: Rsnld.Root(UVap);
                uvap.copy_from_slice(&rsnld.root());
                // OCCT L103.
                self.empty = false;
                // OCCT L104.
                self.tangent = self
                    .my_zer_par_func
                    .is_tangent(&uvap, &mut uvres, &mut best_choix);
                // OCCT L105: pint.SetValue(myZerParFunc.Point(), Uvres(1),
                // Uvres(2), Uvres(3), Uvres(4));
                let p = self.my_zer_par_func.point();
                self.pint
                    .set_value_all(p, uvres[0], uvres[1], uvres[2], uvres[3]);
                if !self.tangent {
                    // OCCT L108-110.
                    self.d3d = self.my_zer_par_func.direction();
                    self.d2d1 = self.my_zer_par_func.direction_on_s1();
                    self.d2d2 = self.my_zer_par_func.direction_on_s2();
                }
            } else {
                // OCCT L115.
                self.empty = true;
            }
        } else {
            // OCCT L120.
            self.empty = true;
        }
        // OCCT L122.
        best_choix
    }

    /// OCCT Perform(Param, Rsnld) — gxx L125-268: returns the best constant
    /// isoparametric to find the next intersection point and stores the
    /// solution point (found from the close point; the choice of the
    /// isoparametric is calculated).
    // The OCCT statement order is preserved verbatim (CurrentChoix is
    // initialized from BestChoix at L158 and read from the first loop
    // iteration on).
    #[allow(unused_assignments)]
    pub fn perform(
        &mut self,
        param: &[f64; 4],
        rsnld: &mut FunctionSetRoot,
    ) -> ConstIsoparametric {
        // OCCT L128-134: gp_Vec DPUV[4]; gp_Pnt P1, P2; double Epsuv[4];
        // NCollection_Array1 Duv(1, 4); double UVd[4], UVf[4]; ChoixIso[4]
        // (uninitialized in OCCT; zero/first-enum defaults stand in).
        let mut dpuv = [DVec3::ZERO; 4];
        let mut epsuv = [0.0f64; 4];
        let mut duv = [0.0f64; 4];
        let mut uvd = [0.0f64; 4];
        let mut uvf = [0.0f64; 4];
        let mut choix_iso_tab = [ConstIsoparametric::UIsoparametricOnCaro1; 4];

        // OCCT L135.
        let mut best_choix = choix_ref(0);
        // OCCT L136-137: the auxillary surfaces.
        let caro1 = self.my_zer_par_func.auxillar_surface1();
        let caro2 = self.my_zer_par_func.auxillar_surface2();

        // OCCT L139-140.
        let (_p1, d1u1, d1v1) = surface_tool::d1(caro1, param[0], param[1]);
        dpuv[0] = d1u1;
        dpuv[1] = d1v1;
        let (_p2, d1u2, d1v2) = surface_tool::d1(caro2, param[2], param[3]);
        dpuv[2] = d1u2;
        dpuv[3] = d1v2;

        // OCCT L142-146.
        epsuv[0] = surface_tool::u_resolution(caro1, CONFUSION);
        epsuv[1] = surface_tool::v_resolution(caro1, CONFUSION);

        epsuv[2] = surface_tool::u_resolution(caro2, CONFUSION);
        epsuv[3] = surface_tool::v_resolution(caro2, CONFUSION);

        // OCCT L148-149.
        for j in 0..=3 {
            uvd[j] = param[j];
        }

        // OCCT L151.
        self.empty = true;

        // OCCT L153-155: bool Tangent = IntImp_ComputeTangence(DPUV, Epsuv,
        // UVd, ChoixIso); (UVd doubles as the Tgduv output buffer).
        let tangent = compute_tangence(&dpuv, &epsuv, &mut uvd, &mut choix_iso_tab);
        if tangent {
            return best_choix;
        }

        // OCCT L157-158.
        let mut i = 0usize;
        let mut current_choix = best_choix; //-- Modif 17 Mai 93

        // OCCT L160-168.
        while self.empty && i <= 3 {
            current_choix = self.perform_with_choice(param, rsnld, choix_iso_tab[i]);
            if !self.empty {
                best_choix = current_choix;
            }
            i += 1;
        }
        if !self.empty {
            // verifier que l on ne deborde pas les frontieres
            // OCCT L171: pint.Parameters(Duv(1), Duv(2), Duv(3), Duv(4));
            let (d1, d2, d3, d4) = self.pint.parameters();
            duv[0] = d1;
            duv[1] = d2;
            duv[2] = d3;
            duv[3] = d4;

            // OCCT L173-181.
            uvd[0] = self.ua0;
            uvd[1] = self.va0;
            uvf[0] = self.ua1;
            uvf[1] = self.va1;

            uvd[2] = self.ub0;
            uvd[3] = self.vb0;
            uvf[2] = self.ub1;
            uvf[3] = self.vb1;

            // OCCT L183-233: the boundary-overflow chain (Nc = 0 for the
            // first surface, 2 for the second; Iiso = the iso index 0..3).
            let nc: usize;
            let iiso: usize;
            if duv[0] <= uvd[0] - epsuv[0] {
                duv[0] = uvd[0];
                nc = 0;
                iiso = 0;
            } else if duv[0] >= uvf[0] + epsuv[0] {
                duv[0] = uvf[0];
                nc = 0;
                iiso = 0;
            } else if duv[1] <= uvd[1] - epsuv[1] {
                duv[1] = uvd[1];
                nc = 0;
                iiso = 1;
            } else if duv[1] >= uvf[1] + epsuv[1] {
                duv[1] = uvf[1];
                nc = 0;
                iiso = 1;
            } else if duv[2] <= uvd[2] - epsuv[2] {
                duv[2] = uvd[2];
                nc = 2;
                iiso = 2;
            } else if duv[2] >= uvf[2] + epsuv[2] {
                duv[2] = uvf[2];
                nc = 2;
                iiso = 2;
            } else if duv[3] <= uvd[3] - epsuv[3] {
                duv[3] = uvd[3];
                nc = 2;
                iiso = 3;
            } else if duv[3] >= uvf[3] + epsuv[3] {
                duv[3] = uvf[3];
                nc = 2;
                iiso = 3;
            } else {
                return best_choix; // on a gagne
            }

            // OCCT L234-236.
            self.empty = true;
            best_choix = choix_ref(iiso); // en attendant
            best_choix = self.perform_with_choice(&duv, rsnld, best_choix);
            if !self.empty {
                // OCCT L238-256: verification sur le carreau reciproque.
                // Nc flips to the reciprocal surface (0<->3, 2<->1); Duv(...)
                // is 1-based while UVd/UVf/Epsuv are 0-based C arrays.
                let mut nc = 3 - nc;
                if duv[nc - 1] <= uvd[nc - 1] - epsuv[nc - 1] {
                    duv[nc - 1] = uvd[nc - 1];
                } else if duv[nc - 1] >= uvf[nc - 1] + epsuv[nc - 1] {
                    duv[nc - 1] = uvf[nc - 1];
                } else if duv[nc] <= uvd[nc] {
                    nc = nc + 1;
                    duv[nc - 1] = uvd[nc - 1];
                } else if duv[nc] >= uvf[nc] {
                    nc = nc + 1;
                    duv[nc - 1] = uvf[nc - 1];
                } else {
                    return best_choix;
                }

                // OCCT L258.
                self.empty = true;

                // OCCT L260-261: if (Nc == 4) Nc = 0;
                let mut nc = nc;
                if nc == 4 {
                    nc = 0;
                }

                // OCCT L263-264.
                best_choix = choix_ref(nc); // en attendant
                best_choix = self.perform_with_choice(&duv, rsnld, best_choix);
            }
        }
        // OCCT L267.
        best_choix
    }

    /// OCCT IsDone() — lxx L19-22.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsEmpty() — lxx L24-29.
    pub fn is_empty(&self) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: IntImp_Int2S::IsEmpty()");
        }
        self.empty
    }

    /// OCCT Point() — lxx L31-38 (the OCCT const-reference return maps to a
    /// clone of the stored point).
    pub fn point(&self) -> PntOn2S {
        if !self.done {
            panic!("StdFail_NotDone: IntImp_Int2S::Point()");
        }
        if self.empty {
            panic!("Standard_DomainError: IntImp_Int2S::Point()");
        }
        self.pint.clone()
    }

    /// OCCT IsTangent() — lxx L40-48.
    pub fn is_tangent(&self) -> bool {
        if !self.done {
            panic!("StdFail_NotDone: IntImp_Int2S::IsTangent()");
        }
        if self.empty {
            panic!("Standard_DomainError: IntImp_Int2S::IsTangent()");
        }
        self.tangent
    }

    /// OCCT Direction() — lxx L50-60.
    pub fn direction(&self) -> DVec3 {
        if !self.done {
            panic!("StdFail_NotDone: IntImp_Int2S::Direction()");
        }
        if self.empty {
            panic!("Standard_DomainError: IntImp_Int2S::Direction()");
        }
        if self.tangent {
            panic!("StdFail_UndefinedDerivative: IntImp_Int2S::Direction()");
        }
        self.d3d
    }

    /// OCCT DirectionOnS1() — lxx L62-72.
    pub fn direction_on_s1(&self) -> DVec2 {
        if !self.done {
            panic!("StdFail_NotDone: IntImp_Int2S::DirectionOnS1()");
        }
        if self.empty {
            panic!("Standard_DomainError: IntImp_Int2S::DirectionOnS1()");
        }
        if self.tangent {
            panic!("StdFail_UndefinedDerivative: IntImp_Int2S::DirectionOnS1()");
        }
        self.d2d1
    }

    /// OCCT DirectionOnS2() — lxx L74-84.
    pub fn direction_on_s2(&self) -> DVec2 {
        if !self.done {
            panic!("StdFail_NotDone: IntImp_Int2S::DirectionOnS2()");
        }
        if self.empty {
            panic!("Standard_DomainError: IntImp_Int2S::DirectionOnS2()");
        }
        if self.tangent {
            panic!("StdFail_UndefinedDerivative: IntImp_Int2S::DirectionOnS2()");
        }
        self.d2d2
    }

    /// OCCT Function() — lxx L86-89.
    pub fn function(&mut self) -> &mut ZerParFunc<'a> {
        &mut self.my_zer_par_func
    }

    /// OCCT ChangePoint() — lxx L91-94.
    pub fn change_point(&mut self) -> &mut PntOn2S {
        &mut self.pint
    }
}

// ---------------------------------------------------------------------------
// PrmPrmSvSurfaces — BRepApprox_ThePrmPrmSvSurfacesOfApprox
// ---------------------------------------------------------------------------

/// OCCT ApproxInt_PrmPrmSvSurfaces (ApproxInt_PrmPrmSvSurfaces.gxx)
/// instantiated as BRepApprox_ThePrmPrmSvSurfacesOfApprox (_0.cxx L28-47)
/// — the ApproxInt_SvSurfaces for two parametric surfaces.  The OCCT
/// non-const Compute/Pnt/SeekPoint/Tangency methods mutate the cached
/// state and the embedded Int2S; through the rcad [`SvSurfaces`] trait
/// (whose receivers take `&self`, like the OCCT calls through the `void*
/// PtrOnmySvSurfaces`) that state is Cell/RefCell.
pub struct PrmPrmSvSurfaces<'a> {
    // OCCT ApproxInt_SvSurfaces base: bool myUseSolver (hxx L91; the base
    // constructor sets it to false).
    my_use_solver: Cell<bool>,
    // OCCT hxx private members, in declaration order.
    my_par_on_s1: RefCell<DVec2>,        // gp_Pnt2d MyParOnS1
    my_par_on_s2: RefCell<DVec2>,        // gp_Pnt2d MyParOnS2
    my_pnt: RefCell<DVec3>,              // gp_Pnt MyPnt
    my_tguv1: RefCell<DVec2>,            // gp_Vec2d MyTguv1
    my_tguv2: RefCell<DVec2>,            // gp_Vec2d MyTguv2
    my_tg: RefCell<DVec3>,               // gp_Vec MyTg
    my_is_tangent: Cell<bool>,           // bool MyIsTangent
    my_has_been_computed: Cell<bool>,    // bool MyHasBeenComputed
    my_par_on_s1bis: RefCell<DVec2>,     // gp_Pnt2d MyParOnS1bis
    my_par_on_s2bis: RefCell<DVec2>,     // gp_Pnt2d MyParOnS2bis
    my_pntbis: RefCell<DVec3>,           // gp_Pnt MyPntbis
    my_tguv1bis: RefCell<DVec2>,         // gp_Vec2d MyTguv1bis
    my_tguv2bis: RefCell<DVec2>,         // gp_Vec2d MyTguv2bis
    my_tgbis: RefCell<DVec3>,            // gp_Vec MyTgbis
    my_is_tangentbis: Cell<bool>,        // bool MyIsTangentbis
    my_has_been_computedbis: Cell<bool>, // bool MyHasBeenComputedbis
    // OCCT TheInt2S MyIntersectionOn2S;
    my_intersection_on_2s: RefCell<Int2S<'a>>,
}

impl<'a> PrmPrmSvSurfaces<'a> {
    /// OCCT ApproxInt_PrmPrmSvSurfaces(Surf1, Surf2) — gxx L29-37 (with the
    /// base-class myUseSolver(false) of ApproxInt_SvSurfaces.hxx).
    pub fn new(surf1: &'a BRepAdaptorSurface<'a>, surf2: &'a BRepAdaptorSurface<'a>) -> Self {
        PrmPrmSvSurfaces {
            my_use_solver: Cell::new(false),
            my_par_on_s1: RefCell::new(DVec2::ZERO), // OCCT gp_Pnt2d default
            my_par_on_s2: RefCell::new(DVec2::ZERO), // OCCT gp_Pnt2d default
            my_pnt: RefCell::new(DVec3::ZERO),       // OCCT gp_Pnt default
            my_tguv1: RefCell::new(DVec2::ZERO),     // OCCT gp_Vec2d default
            my_tguv2: RefCell::new(DVec2::ZERO),     // OCCT gp_Vec2d default
            my_tg: RefCell::new(DVec3::ZERO),        // OCCT gp_Vec default
            // OCCT L31-34.
            my_is_tangent: Cell::new(false),
            my_has_been_computed: Cell::new(false),
            my_is_tangentbis: Cell::new(false),
            my_has_been_computedbis: Cell::new(false),
            my_par_on_s1bis: RefCell::new(DVec2::ZERO), // OCCT gp_Pnt2d default
            my_par_on_s2bis: RefCell::new(DVec2::ZERO), // OCCT gp_Pnt2d default
            my_pntbis: RefCell::new(DVec3::ZERO),       // OCCT gp_Pnt default
            my_tguv1bis: RefCell::new(DVec2::ZERO),     // OCCT gp_Vec2d default
            my_tguv2bis: RefCell::new(DVec2::ZERO),     // OCCT gp_Vec2d default
            my_tgbis: RefCell::new(DVec3::ZERO),        // OCCT gp_Vec default
            // OCCT L35: MyIntersectionOn2S(Surf1, Surf2, TOLTANGENCY).
            my_intersection_on_2s: RefCell::new(Int2S::new(surf1, surf2, TOLTANGENCY)),
        }
    }

    /// OCCT Compute(u1, v1, u2, v2, P, Tg, Tguv1, Tguv2) — gxx L44-228:
    /// computes the point on the curve, the 3d and 2d tangents of the
    /// intersection curve and the parameters on the surfaces; caches the
    /// last result (with the OCCT two-slot current/bis swap logic).
    // The OCCT work locals (DeltaU/DeltaV, Pbid, TU/TV, TUTU..DIS) are uninitialized then assigned.
    #[allow(clippy::too_many_arguments, unused_assignments)]
    pub fn compute(
        &self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        p: &mut DVec3,
        tg: &mut DVec3,
        tguv1: &mut DVec2,
        tguv2: &mut DVec2,
    ) -> bool {
        // OCCT L54-57.
        let tu1 = *u1;
        let tu2 = *u2;
        let tv1 = *v1;
        let tv2 = *v2;

        // OCCT L59-77.
        if self.my_has_been_computed.get() {
            let my_par_on_s1 = *self.my_par_on_s1.borrow();
            let my_par_on_s2 = *self.my_par_on_s2.borrow();
            if (my_par_on_s1.x == *u1)
                && (my_par_on_s1.y == *v1)
                && (my_par_on_s2.x == *u2)
                && (my_par_on_s2.y == *v2)
            {
                return self.my_is_tangent.get();
            } else if self.my_has_been_computedbis.get() == false {
                *self.my_tgbis.borrow_mut() = *self.my_tg.borrow();
                *self.my_tguv1bis.borrow_mut() = *self.my_tguv1.borrow();
                *self.my_tguv2bis.borrow_mut() = *self.my_tguv2.borrow();
                *self.my_pntbis.borrow_mut() = *self.my_pnt.borrow();
                *self.my_par_on_s1bis.borrow_mut() = *self.my_par_on_s1.borrow();
                *self.my_par_on_s2bis.borrow_mut() = *self.my_par_on_s2.borrow();
                self.my_is_tangentbis.set(self.my_is_tangent.get());
                self.my_has_been_computedbis
                    .set(self.my_has_been_computed.get());
            }
        }

        // OCCT L78-110.
        if self.my_has_been_computedbis.get() {
            let my_par_on_s1bis = *self.my_par_on_s1bis.borrow();
            let my_par_on_s2bis = *self.my_par_on_s2bis.borrow();
            if (my_par_on_s1bis.x == *u1)
                && (my_par_on_s1bis.y == *v1)
                && (my_par_on_s2bis.x == *u2)
                && (my_par_on_s2bis.y == *v2)
            {
                // OCCT L84-90: TV, TV1, TV2, TP, TP1, TP2, TB — the current
                // values.
                let tv = *self.my_tg.borrow();
                let tv1 = *self.my_tguv1.borrow();
                let tv2 = *self.my_tguv2.borrow();
                let tp = *self.my_pnt.borrow();
                let tp1 = *self.my_par_on_s1.borrow();
                let tp2 = *self.my_par_on_s2.borrow();
                let tb = self.my_is_tangent.get();

                // OCCT L92-98: current <- bis.
                *self.my_tg.borrow_mut() = *self.my_tgbis.borrow();
                *self.my_tguv1.borrow_mut() = *self.my_tguv1bis.borrow();
                *self.my_tguv2.borrow_mut() = *self.my_tguv2bis.borrow();
                *self.my_pnt.borrow_mut() = *self.my_pntbis.borrow();
                *self.my_par_on_s1.borrow_mut() = *self.my_par_on_s1bis.borrow();
                *self.my_par_on_s2.borrow_mut() = *self.my_par_on_s2bis.borrow();
                self.my_is_tangent.set(self.my_is_tangentbis.get());

                // OCCT L100-106: bis <- saved.
                *self.my_tgbis.borrow_mut() = tv;
                *self.my_tguv1bis.borrow_mut() = tv1;
                *self.my_tguv2bis.borrow_mut() = tv2;
                *self.my_pntbis.borrow_mut() = tp;
                *self.my_par_on_s1bis.borrow_mut() = tp1;
                *self.my_par_on_s2bis.borrow_mut() = tp2;
                self.my_is_tangentbis.set(tb);

                return self.my_is_tangent.get();
            }
        }

        // OCCT L112.
        self.my_is_tangent.set(true);

        // OCCT L114-119: NCollection_Array1<double> Param(aParam[0], 1, 4).
        let mut param = [0.0f64; 4];
        param[0] = *u1;
        param[1] = *v1;
        param[2] = *u2;
        param[3] = *v2;

        // OCCT L120: math_FunctionSetRoot Rsnld(MyIntersectionOn2S.Function());
        // (the 1-argument OCCT ctor: NbIterations = 100, Tol zero-initialized).
        let mut rsnld = {
            let mut my_int = self.my_intersection_on_2s.borrow_mut();
            FunctionSetRoot::new(&*my_int.function(), &[0.0, 0.0, 0.0], 100)
        };

        // OCCT L121.
        self.my_intersection_on_2s
            .borrow_mut()
            .perform(&param, &mut rsnld);

        // OCCT L122-126.
        if !self.my_intersection_on_2s.borrow().is_done() {
            self.my_has_been_computed.set(false);
            self.my_has_been_computedbis.set(false);
            return false;
        }

        // OCCT L127-134.
        if self.my_intersection_on_2s.borrow().is_empty() {
            self.my_is_tangent.set(false);
            self.my_has_been_computed.set(false);
            self.my_has_been_computedbis.set(false);
            return false;
        }

        // OCCT L135-136: MyPnt = P = MyIntersectionOn2S.Point().Value();
        self.my_has_been_computed.set(true);
        let my_pnt = self.my_intersection_on_2s.borrow().point().value();
        *self.my_pnt.borrow_mut() = my_pnt;
        *p = my_pnt;

        // OCCT L138-140.  The OCCT statement order is preserved verbatim:
        // the outputs u1..v2 receive the solution parameters, while
        // MyParOnS1/MyParOnS2 store the SAVED INPUT values tu1..tv2.
        let (nu1, nv1, nu2, nv2) = self.my_intersection_on_2s.borrow().point().parameters();
        *u1 = nu1;
        *v1 = nv1;
        *u2 = nu2;
        *v2 = nv2;
        *self.my_par_on_s1.borrow_mut() = DVec2::new(tu1, tv1);
        *self.my_par_on_s2.borrow_mut() = DVec2::new(tu2, tv2);

        // OCCT L142-147.
        if self.my_intersection_on_2s.borrow().is_tangent() {
            self.my_is_tangent.set(false);
            self.my_has_been_computed.set(false);
            self.my_has_been_computedbis.set(false);
            return false;
        }

        // OCCT L148-150.
        let d3d = self.my_intersection_on_2s.borrow().direction();
        *tg = d3d;
        *self.my_tg.borrow_mut() = d3d;
        let d2d1 = self.my_intersection_on_2s.borrow().direction_on_s1();
        *tguv1 = d2d1;
        *self.my_tguv1.borrow_mut() = d2d1;
        let d2d2 = self.my_intersection_on_2s.borrow().direction_on_s2();
        *tguv2 = d2d2;
        *self.my_tguv2.borrow_mut() = d2d2;

        // OCCT L172-173: Tg.Normalize(); MyTg = Tg;
        let tg_norm = tg.normalize();
        *tg = tg_norm;
        *self.my_tg.borrow_mut() = tg_norm;

        // OCCT L175-178 (the work declarations; OCCT leaves them
        // uninitialized, zero defaults stand in).
        let (mut delta_u, mut delta_v) = (0.0, 0.0);
        let (mut _p_bid, mut t_u, mut t_v) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO); // Pbid, TU, TV
        let (mut t_u_tu, mut t_v_tv, mut t_u_tv, mut t_gtu, mut t_gt_v, mut dis) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        // OCCT L179-182: calcul de Tguv1.
        let my_surf1 = self
            .my_intersection_on_2s
            .borrow_mut()
            .function()
            .auxillar_surface1();
        let (pb, du, dv) = surface_tool::d1(my_surf1, *u1, *v1);
        _p_bid = pb;
        t_u = du;
        t_v = dv;

        // OCCT L184-189.
        t_u_tu = t_u.dot(t_u);
        t_v_tv = t_v.dot(t_v);
        t_u_tv = t_u.dot(t_v);
        t_gtu = tg.dot(t_u);
        t_gt_v = tg.dot(t_v);
        dis = t_u_tu * t_v_tv - t_u_tv * t_u_tv;

        // OCCT L190-195.
        if dis.abs() < ANGULAR {
            self.my_is_tangent.set(false);
            self.my_has_been_computed.set(false);
            self.my_has_been_computedbis.set(false);
            return false;
        }

        // OCCT L197-198.
        delta_u = (t_gtu * t_v_tv - t_gt_v * t_u_tv) / dis;
        delta_v = (t_gt_v * t_u_tu - t_gtu * t_u_tv) / dis;

        // OCCT L200-201.
        *tguv1 = DVec2::new(delta_u, delta_v);
        *self.my_tguv1.borrow_mut() = *tguv1;

        // OCCT L203-206: calcul de Tguv2.
        let my_surf2 = self
            .my_intersection_on_2s
            .borrow_mut()
            .function()
            .auxillar_surface2();
        let (pb, du, dv) = surface_tool::d1(my_surf2, *u2, *v2);
        _p_bid = pb;
        t_u = du;
        t_v = dv;

        // OCCT L208-213.
        t_u_tu = t_u.dot(t_u);
        t_v_tv = t_v.dot(t_v);
        t_u_tv = t_u.dot(t_v);
        t_gtu = tg.dot(t_u);
        t_gt_v = tg.dot(t_v);
        dis = t_u_tu * t_v_tv - t_u_tv * t_u_tv;

        // OCCT L214-219.
        if dis.abs() < ANGULAR {
            self.my_is_tangent.set(false);
            self.my_has_been_computed.set(false);
            self.my_has_been_computedbis.set(false);
            return false;
        }

        // OCCT L221-222.
        delta_u = (t_gtu * t_v_tv - t_gt_v * t_u_tv) / dis;
        delta_v = (t_gt_v * t_u_tu - t_gtu * t_u_tv) / dis;

        // OCCT L224-225.
        *tguv2 = DVec2::new(delta_u, delta_v);
        *self.my_tguv2.borrow_mut() = *tguv2;

        // OCCT L227.
        true
    }
}

// The ApproxInt_SvSurfaces virtuals (the Compute above plus the forwarding
// wrappers, gxx L230-332).
impl SvSurfaces for PrmPrmSvSurfaces<'_> {
    /// OCCT ApproxInt_PrmPrmSvSurfaces::Compute — gxx L44-228.
    fn compute(
        &self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        pt: &mut DVec3,
        tg: &mut DVec3,
        tguv1: &mut DVec2,
        tguv2: &mut DVec2,
    ) -> bool {
        PrmPrmSvSurfaces::compute(self, u1, v1, u2, v2, pt, tg, tguv1, tguv2)
    }

    /// OCCT Pnt(u1, v1, u2, v2, P) — gxx L232-247.
    fn pnt(&self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut DVec3) {
        // OCCT L238-244: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        // double tu1..tv2.
        let (mut a_p, mut a_t) = (DVec3::ZERO, DVec3::ZERO);
        let (mut a_ts1, mut a_ts2) = (DVec2::ZERO, DVec2::ZERO);
        let (mut tu1, mut tu2, mut tv1, mut tv2) = (u1, u2, v1, v2);
        self.compute(
            &mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2,
        );
        // OCCT L246: P = MyPnt;
        *p = *self.my_pnt.borrow();
    }

    /// OCCT SeekPoint(u1, v1, u2, v2, Point) — gxx L254-272: computes the
    /// point on the curve and the parameters on the surfaces.
    fn seek_point(&self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool {
        // OCCT L260-266: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        // double tu1..tv2.
        let (mut a_p, mut a_t) = (DVec3::ZERO, DVec3::ZERO);
        let (mut a_ts1, mut a_ts2) = (DVec2::ZERO, DVec2::ZERO);
        let (mut tu1, mut tu2, mut tv1, mut tv2) = (u1, u2, v1, v2);
        if !self.compute(
            &mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2,
        ) {
            return false;
        }

        // OCCT L270: Point.SetValue(aP, tu1, tv1, tu2, tv2); — the SAVED
        // INPUT parameters (the OCCT literals, not the refined outputs).
        point.set_value_all(a_p, tu1, tv1, tu2, tv2);
        true
    }

    /// OCCT Tangency(u1, v1, u2, v2, T) — gxx L276-292.
    fn tangency(&self, u1: f64, v1: f64, u2: f64, v2: f64, t: &mut DVec3) -> bool {
        // OCCT L282-289: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        // double tu1..tv2.
        let (mut a_p, mut a_t) = (DVec3::ZERO, DVec3::ZERO);
        let (mut a_ts1, mut a_ts2) = (DVec2::ZERO, DVec2::ZERO);
        let (mut tu1, mut tu2, mut tv1, mut tv2) = (u1, u2, v1, v2);
        let t_flag = self.compute(
            &mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2,
        );
        // OCCT L290: T = MyTg;
        *t = *self.my_tg.borrow();
        t_flag
    }

    /// OCCT TangencyOnSurf1(u1, v1, u2, v2, T) — gxx L296-312.
    fn tangency_on_surf_1(&self, u1: f64, v1: f64, u2: f64, v2: f64, t: &mut DVec2) -> bool {
        // OCCT L302-309: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        // double tu1..tv2.
        let (mut a_p, mut a_t) = (DVec3::ZERO, DVec3::ZERO);
        let (mut a_ts1, mut a_ts2) = (DVec2::ZERO, DVec2::ZERO);
        let (mut tu1, mut tu2, mut tv1, mut tv2) = (u1, u2, v1, v2);
        let t_flag = self.compute(
            &mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2,
        );
        // OCCT L310: T = MyTguv1;
        *t = *self.my_tguv1.borrow();
        t_flag
    }

    /// OCCT TangencyOnSurf2(u1, v1, u2, v2, T) — gxx L316-332.
    fn tangency_on_surf_2(&self, u1: f64, v1: f64, u2: f64, v2: f64, t: &mut DVec2) -> bool {
        // OCCT L321-329: gp_Pnt aP; gp_Vec aT; gp_Vec2d aTS1, aTS2;
        // double tu1..tv2.
        let (mut a_p, mut a_t) = (DVec3::ZERO, DVec3::ZERO);
        let (mut a_ts1, mut a_ts2) = (DVec2::ZERO, DVec2::ZERO);
        let (mut tu1, mut tu2, mut tv1, mut tv2) = (u1, u2, v1, v2);
        let t_flag = self.compute(
            &mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2,
        );
        // OCCT L330: T = MyTguv2;
        *t = *self.my_tguv2.borrow();
        t_flag
    }

    /// OCCT ApproxInt_SvSurfaces::SetUseSolver (hxx L93).
    fn set_use_solver(&self, the_use_sol: bool) {
        self.my_use_solver.set(the_use_sol);
    }

    /// OCCT ApproxInt_SvSurfaces::GetUseSolver (hxx L95).
    fn get_use_solver(&self) -> bool {
        self.my_use_solver.get()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Line3, Plane, Point3, Surface3, Vec3};
    use rcad_kernel::topo::topods::{BRepBuilder, Shape};

    /// A plane face over the given UV window; the plane is
    /// P(u, v) = origin + u*u_dir + v*v_dir.
    fn add_plane_face(
        brep: &mut rcad_kernel::BRep,
        origin: Point3,
        normal: Vec3,
        u_dir: Vec3,
        v_dir: Vec3,
        uv: [f64; 4],
    ) -> Shape {
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(brep, origin, 1e-7);
        let v2 = b.add_vertex(brep, origin + u_dir, 1e-7);
        let e = b.add_edge(
            brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin,
                direction: u_dir,
            })),
            v1,
            v2,
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e]);
        brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin,
                normal,
                u_dir,
                v_dir,
            })),
            wire,
            Vec::new(),
            None,
            Some(uv),
            Vec::new(),
            true,
        )
    }

    /// Two orthogonal planes over [0, 1] x [0, 1]:
    /// S1(u1, v1) = (u1, v1, 0)  (normal +z),
    /// S2(u2, v2) = (0, u2, v2) (+x normal); line (0, t, 0): u1=0, v2=0, v1=u2.
    fn crossing_planes_fixture() -> (rcad_kernel::BRep, Shape, Shape) {
        let mut brep = rcad_kernel::BRep::new();
        let f1 = add_plane_face(&mut brep, Point3::ZERO, Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), [0.0, 1.0, 0.0, 1.0]);
        let f2 = add_plane_face(&mut brep, Point3::ZERO, Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0), [0.0, 1.0, 0.0, 1.0]);
        (brep, f1, f2)
    }

    /// OCCT anchor: ComputeParameters distributes Param over the fixed
    /// isoparametric and the unknowns, sets the widened bounds and the
    /// resolution tolerances (IntImp_ZerParFunc.gxx L231-328), and Value /
    /// Derivatives produce the analytic residual and Jacobian (L65-161).
    #[test]
    fn zer_par_func_value_derivatives_planes() {
        let (brep, f1, f2) = crossing_planes_fixture();
        let s1 = BRepAdaptorSurface::initialize_face(&brep, &f1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep, &f2, true);
        let mut fct = ZerParFunc::new(&s1, &s2);

        // The intervals read from the surfaces (gxx L35-49).
        assert_eq!(fct.auxillar_surface1().first_u_parameter(), 0.0);
        assert_eq!(fct.auxillar_surface2().first_u_parameter(), 0.0);

        // UIsoparametricOnCaro1: paramConst = Param(1); unknowns (v1, u2, v2).
        let param = [0.25, 0.5, 0.3, 0.2];
        let mut uvap = [0.0f64; 3];
        let mut born_inf = [0.0f64; 3];
        let mut born_sup = [0.0f64; 3];
        let mut tolerance = [0.0f64; 3];
        fct.compute_parameters(
            ConstIsoparametric::UIsoparametricOnCaro1,
            &param,
            &mut uvap,
            &mut born_inf,
            &mut born_sup,
            &mut tolerance,
        );
        assert_eq!(uvap, [0.5, 0.3, 0.2]);
        // Bounds: v of S1 then u/v of S2, widened by 1% (gxx L248-254, L319-327).
        assert!((born_inf[0] - (-0.01)).abs() < 1e-15);
        assert!((born_sup[0] - 1.01).abs() < 1e-15);
        assert!((born_inf[1] - (-0.01)).abs() < 1e-15);
        assert!((born_sup[2] - 1.01).abs() < 1e-15);
        // Tolerances = vres1, ures2, vres2 (gxx L256-258); a plane resolves
        // r3d identically.
        assert_eq!(tolerance, [CONFUSION, CONFUSION, CONFUSION]);

        // Value at X = (0.5, 0.3, 0.2): F = S1(0.25, 0.5) - S2(0.3, 0.2).
        let x = [0.5, 0.3, 0.2];
        let mut fv = [0.0f64; 3];
        assert!(FunctionSetWithDerivatives::value(&mut fct, &x, &mut fv));
        assert!((fv[0] - 0.25).abs() < 1e-15);
        assert!((fv[1] - 0.2).abs() < 1e-15);
        assert!((fv[2] + 0.2).abs() < 1e-15);
        // Root() = f0^2 + f1^2 + f2^2 (lxx L21-25).
        assert!((fct.root() - (0.0625 + 0.04 + 0.04)).abs() < 1e-15);
        // Point() = the mid-point of the two surface points (lxx L27-30):
        // ((0.25, 0.5, 0) + (0, 0.3, 0.2)) / 2 = (0.125, 0.4, 0.1).
        let mid = fct.point();
        assert!((mid.x - 0.125).abs() < 1e-15);
        assert!((mid.y - 0.4).abs() < 1e-15);
        assert!((mid.z - 0.1).abs() < 1e-15);

        // Derivatives for UIsoparametricOnCaro1 (gxx L104-114):
        // D = [[0, 0, 0], [1, -1, 0], [0, 0, -1]].
        let mut df = vec![vec![0.0f64; 3]; 3];
        assert!(FunctionSetWithDerivatives::derivatives(&mut fct, &x, &mut df));
        assert_eq!(
            df,
            vec![vec![0.0, 0.0, 0.0], vec![1.0, -1.0, 0.0], vec![0.0, 0.0, -1.0]]
        );

        // Values(X, F, D) — same Jacobian plus the residual (gxx L163-229).
        let mut fv2 = [0.0f64; 3];
        let mut df2 = vec![vec![0.0f64; 3]; 3];
        assert!(FunctionSetWithDerivatives::values(&mut fct, &x, &mut fv2, &mut df2));
        assert!((fv2[0] - 0.25).abs() < 1e-15);
        assert_eq!(df2[1][1], -1.0);
        assert_eq!(df2[2][2], -1.0);
    }

    /// OCCT anchor: IntImp_ComputeTangence on the two crossing planes —
    /// tgduv = (0, 1, 1, 0), not tangent, TabIso[0] = VIsoparametricOnCaro1
    /// (IntImp_ComputeTangence.cxx L112-166).
    #[test]
    fn compute_tangence_crossing_planes_ranking() {
        let dpuv = [
            DVec3::new(1.0, 0.0, 0.0), // dP/du caro 1
            DVec3::new(0.0, 1.0, 0.0), // dP/dv caro 1
            DVec3::new(0.0, 1.0, 0.0), // dP/du caro 2
            DVec3::new(0.0, 0.0, 1.0), // dP/dv caro 2
        ];
        let eps_uv = [CONFUSION; 4];
        let mut tgduv = [0.0f64; 4];
        let mut tab_iso = [ConstIsoparametric::UIsoparametricOnCaro1; 4];
        assert!(!compute_tangence(&dpuv, &eps_uv, &mut tgduv, &mut tab_iso));
        assert!(tgduv[0].abs() < 1e-15);
        assert!((tgduv[1] - 1.0).abs() < 1e-15);
        assert!((tgduv[2] - 1.0).abs() < 1e-15);
        assert!(tgduv[3].abs() < 1e-15);
        // The best next choice (TabIso[0]) ranks V1 first for this geometry.
        assert_eq!(tab_iso[0], ConstIsoparametric::VIsoparametricOnCaro1);
        assert_eq!(choix_ref(0), ConstIsoparametric::UIsoparametricOnCaro1);

        // A degenerate derivative (zero dP/dv on caro 1) is tangent (cxx L73-80).
        let degenerate = [
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::ZERO,
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
        ];
        assert!(compute_tangence(
            &degenerate,
            &eps_uv,
            &mut tgduv,
            &mut tab_iso
        ));
    }

    /// OCCT anchor: Int2S::Perform from a close point converges to the
    /// exact intersection point of the two planes (u1 = 0, v2 = 0, v1 = u2
    /// = 0.4 here), with the 2d/3d tangents (IntImp_Int2S.gxx L81-268).
    #[test]
    fn int2s_perform_converges_to_plane_line() {
        let (brep, f1, f2) = crossing_planes_fixture();
        let s1 = BRepAdaptorSurface::initialize_face(&brep, &f1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep, &f2, true);

        let mut int2s = Int2S::new(&s1, &s2, TOLTANGENCY);
        assert!(int2s.is_done());
        assert!(int2s.is_empty());

        // OCCT usage (IntImp_Int2S.hxx L57-64): rsnld(inter.Function()).
        let mut rsnld = FunctionSetRoot::new(int2s.function(), &[0.0, 0.0, 0.0], 15);
        let param = [0.3, 0.4, 0.2, 0.1];
        let best = int2s.perform(&param, &mut rsnld);

        assert!(!int2s.is_empty());
        let p = int2s.point();
        let v = p.value();
        // The solution point on the line: (0, 0.4, 0).
        assert!(v.x.abs() < 1e-12, "x={}", v.x);
        assert!((v.y - 0.4).abs() < 1e-12, "y={}", v.y);
        assert!(v.z.abs() < 1e-12, "z={}", v.z);
        // pint parameters come from IsTangent's Uvres in the canonical
        // (u1, v1, u2, v2) order: (0, 0.4, 0.4, 0).
        let (pu1, pv1, pu2, pv2) = p.parameters();
        assert!(pu1.abs() < 1e-9, "u1={}", pu1);
        assert!((pv1 - 0.4).abs() < 1e-9, "v1={}", pv1);
        assert!((pu2 - 0.4).abs() < 1e-9, "u2={}", pu2);
        assert!(pv2.abs() < 1e-9, "v2={}", pv2);

        // Transverse intersection: the 3d tangent is (0, 1, 0); the 2d
        // tangents are (0, 1) on S1 and (1, 0) on S2.
        assert!(!int2s.is_tangent());
        let d3 = int2s.direction();
        assert!(d3.x.abs() < 1e-12 && (d3.y - 1.0).abs() < 1e-12 && d3.z.abs() < 1e-12);
        let d1 = int2s.direction_on_s1();
        assert!(d1.x.abs() < 1e-12 && (d1.y - 1.0).abs() < 1e-12);
        let d2 = int2s.direction_on_s2();
        assert!((d2.x - 1.0).abs() < 1e-12 && d2.y.abs() < 1e-12);

        // The ranked VIsoparametricOnCaro1 (fixed iso of the last Perform).
        assert_eq!(best, ConstIsoparametric::VIsoparametricOnCaro1);

        // ChangePoint gives write access to the stored point (lxx L91-94).
        int2s.change_point().set_value_pt(DVec3::new(0.0, 0.5, 0.0));
        assert_eq!(int2s.point().value(), DVec3::new(0.0, 0.5, 0.0));
    }

    /// OCCT anchor: IntImp_Int2S(Param, S1, S2, TolTangency) solves in the
    /// constructor (gxx L53-79).
    #[test]
    fn int2s_ctor_with_param_solves() {
        let (brep, f1, f2) = crossing_planes_fixture();
        let s1 = BRepAdaptorSurface::initialize_face(&brep, &f1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep, &f2, true);

        let param = [0.05, 0.75, 0.6, 0.3];
        let int2s = Int2S::new_with_param(&param, &s1, &s2, TOLTANGENCY);
        assert!(int2s.is_done());
        assert!(!int2s.is_empty());
        let v = int2s.point().value();
        // The converged point: (0, 0.75, 0).
        assert!(v.x.abs() < 1e-12, "x={}", v.x);
        assert!((v.y - 0.75).abs() < 1e-12, "y={}", v.y);
        assert!(v.z.abs() < 1e-12, "z={}", v.z);
    }

    /// OCCT anchor: PrmPrmSvSurfaces::Compute refines the parameters to the
    /// intersection and fills the 3d point/tangent and the 2d tangents
    /// (ApproxInt_PrmPrmSvSurfaces.gxx L44-228); Pnt/SeekPoint/Tangency/
    /// TangencyOnSurf1/2 forward through the cache (L232-332) and
    /// Set/GetUseSolver store the base-class flag (hxx L93-95).
    #[test]
    fn prm_prm_sv_surfaces_compute_analytic() {
        let (brep, f1, f2) = crossing_planes_fixture();
        let s1 = BRepAdaptorSurface::initialize_face(&brep, &f1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep, &f2, true);
        let sv = PrmPrmSvSurfaces::new(&s1, &s2);

        assert!(!sv.get_use_solver());
        sv.set_use_solver(true);
        assert!(sv.get_use_solver());
        sv.set_use_solver(false);

        // First Compute (cache miss) from the close point (0.3, 0.4, 0.2,
        // 0.1): the refined parameters land on the line u1=0, v2=0, v1=u2.
        let mut u1 = 0.3;
        let mut v1 = 0.4;
        let mut u2 = 0.2;
        let mut v2 = 0.1;
        let mut p = DVec3::ZERO;
        let mut tg = DVec3::ZERO;
        let mut tguv1 = DVec2::ZERO;
        let mut tguv2 = DVec2::ZERO;
        assert!(sv.compute(
            &mut u1, &mut v1, &mut u2, &mut v2, &mut p, &mut tg, &mut tguv1, &mut tguv2
        ));
        assert!(u1.abs() < 1e-9, "u1={}", u1);
        assert!((v1 - 0.4).abs() < 1e-12, "v1={}", v1);
        assert!((u2 - 0.4).abs() < 1e-9, "u2={}", u2);
        assert!(v2.abs() < 1e-9, "v2={}", v2);
        assert!(p.x.abs() < 1e-12 && (p.y - 0.4).abs() < 1e-12 && p.z.abs() < 1e-12);
        // Tg = (0, 1, 0); Tguv1 = (0, 1); Tguv2 = (1, 0) (gxx L148-225).
        assert!(tg.x.abs() < 1e-12 && (tg.y - 1.0).abs() < 1e-12 && tg.z.abs() < 1e-12);
        assert!(tguv1.x.abs() < 1e-12 && (tguv1.y - 1.0).abs() < 1e-12);
        assert!((tguv2.x - 1.0).abs() < 1e-12 && tguv2.y.abs() < 1e-12);

        // Pnt returns the cached MyPnt (gxx L246).
        let mut pnt = DVec3::ZERO;
        sv.pnt(0.3, 0.4, 0.2, 0.1, &mut pnt);
        assert_eq!(pnt, p);

        // SeekPoint from a fresh close point (0.15, 0.6, 0.5, 0.2): the
        // local tu1..tv2 are passed BY REFERENCE into Compute (OCCT
        // double&), so Compute overwrites them with the refined parameters
        // and Point.SetValue receives (aP, refined u1, v1, u2, v2) =
        // (0, 0.6, 0, 0.6, 0) — the OCCT literal behaviour of gxx L267-271.
        let mut point = PntOn2S::new();
        assert!(sv.seek_point(0.15, 0.6, 0.5, 0.2, &mut point));
        let pv = point.value();
        assert!(pv.x.abs() < 1e-12 && (pv.y - 0.6).abs() < 1e-12 && pv.z.abs() < 1e-12);
        assert_eq!(point.parameters(), (0.0, 0.6, 0.6, 0.0));

        // Tangency/OnSurf1/2 forward MyTg/MyTguv1/MyTguv2 (gxx L290/310/330).
        let mut t = DVec3::ZERO;
        assert!(sv.tangency(0.3, 0.4, 0.2, 0.1, &mut t));
        assert_eq!(t, tg);
        let mut t2d = DVec2::ZERO;
        assert!(sv.tangency_on_surf_1(0.3, 0.4, 0.2, 0.1, &mut t2d));
        assert_eq!(t2d, tguv1);
        assert!(sv.tangency_on_surf_2(0.3, 0.4, 0.2, 0.1, &mut t2d));
        assert_eq!(t2d, tguv2);
    }

    /// OCCT anchor: Compute on two parallel offset planes finds no zero of
    /// the residual (Root stays above the squared tangency tolerance), so
    /// IsEmpty holds and Compute returns False with MyIsTangent = False
    /// (ApproxInt_PrmPrmSvSurfaces.gxx L127-134).
    #[test]
    fn prm_prm_sv_surfaces_parallel_planes_rejected() {
        let mut brep = rcad_kernel::BRep::new();
        let f1 = add_plane_face(&mut brep, Point3::ZERO, Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), [0.0, 1.0, 0.0, 1.0]);
        let f2 = add_plane_face(&mut brep, Point3::new(0.0, 0.0, 0.5),
            Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0),
            [0.0, 1.0, 0.0, 1.0]);
        let s1 = BRepAdaptorSurface::initialize_face(&brep, &f1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep, &f2, true);
        let sv = PrmPrmSvSurfaces::new(&s1, &s2);

        let mut u1 = 0.3;
        let mut v1 = 0.4;
        let mut u2 = 0.2;
        let mut v2 = 0.1;
        let mut p = DVec3::ZERO;
        let mut tg = DVec3::ZERO;
        let mut tguv1 = DVec2::ZERO;
        let mut tguv2 = DVec2::ZERO;
        assert!(!sv.compute(
            &mut u1, &mut v1, &mut u2, &mut v2, &mut p, &mut tg, &mut tguv1, &mut tguv2
        ));
        // The failure paths reset the flags (gxx L124, L132): the retry fails.
        let mut point = PntOn2S::new();
        assert!(!sv.seek_point(0.3, 0.4, 0.2, 0.1, &mut point));
    }
}
