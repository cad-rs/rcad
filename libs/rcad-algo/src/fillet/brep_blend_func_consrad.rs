//! OCCT BlendFunc_ConstRad (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_ConstRad.hxx (L27-279) + BlendFunc_ConstRad.cxx (whole file
//! L45-1956).  [`BlendFuncConstRadInv`]
//! (BlendFunc_ConstRadInv.hxx L26-78 + BlendFunc_ConstRadInv.cxx L25-646)
//! lives in [`super::brep_blend_func_consrad_c`]; ComputeValues and the four
//! Section overloads live in [`super::brep_blend_func_consrad_b`].
//!
//! Architecture mappings: `class BlendFunc_ConstRad : public Blend_Function`
//! is expressed by implementing the [`BlendFunction`] trait and the
//! `math_FunctionSetWithDerivatives` base; `tcurv` (OCCT handle aliasing
//! `curv` until Set(First, Last)) maps to an owned trimmed [`Curve3`] copy
//! plus an accessor; `math_Vector`/`math_Matrix` map to `[f64; 4]` /
//! `Vec<Vec<f64>>` (OCCT D(i, j) -> d[i - 1][j - 1]); `BlendFunc_Tensor`
//! maps to [`BlendFuncTensor`]; `gp_Circ` maps to the kernel [`Circle3`].
//! Pending kernel dependencies (marked GAP, plan 0.6): math_SVD and
//! GeomFill::GetCircle.

use glam::{DVec2, DVec3};

use rcad_kernel::base::convert::{geom_convert_curve_to_bspline_curve, ConvertParameterisation};
use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Circle3, Curve2d, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{GeomAbsShape, MatD, VecD};

use super::brep_blend_func::{
    blend_func_get_minimal_weights, blend_func_get_shape, blend_func_next_shape,
    BlendFuncSectionShape, ConvertParameterisationType,
};
use super::brep_blend_function::{BlendAppFunction, BlendFunction};
use super::brep_blend_point::BlendPoint;

/// OCCT Eps constant (BlendFunc_ConstRad.cxx L41).
pub(crate) const EPS: f64 = 1.0e-15;

/// OCCT gp::Resolution() (gp.hxx) — used by ElCLib::CircleParameter.
const GP_RESOLUTION: f64 = 1.0e-12;

/// OCCT BlendFunc_Tensor — used to store the "gradient of gradient"
/// (BlendFunc_Tensor.hxx L25-70, BlendFunc_Tensor.cxx L20-52 and
/// BlendFunc_Tensor.lxx L19-27).
pub struct BlendFuncTensor {
    tab: Vec<f64>,
    nbrow: i32,
    nbcol: i32,
    nbmat: i32,
    nbmtcl: i32,
}

impl BlendFuncTensor {
    /// OCCT BlendFunc_Tensor(NbRow, NbCol, NbMat) (BlendFunc_Tensor.cxx
    /// L20-27).
    pub fn new(nb_row: i32, nb_col: i32, nb_mat: i32) -> Self {
        BlendFuncTensor {
            tab: vec![0.0; (nb_row * nb_mat * nb_col) as usize],
            nbrow: nb_row,
            nbcol: nb_col,
            nbmat: nb_mat,
            nbmtcl: nb_mat * nb_col,
        }
    }

    /// OCCT Init(InitialValue) (BlendFunc_Tensor.cxx L29-34) — initialize
    /// all the elements of a Tensor to InitialValue.
    pub fn init(&mut self, initial_value: f64) {
        for v in self.tab.iter_mut() {
            *v = initial_value;
        }
    }

    /// OCCT Value(Row, Col, Mat) (BlendFunc_Tensor.lxx L19-22) — 1-based
    /// indexing, read access.
    #[inline]
    pub fn value(&self, row: i32, col: i32, mat: i32) -> f64 {
        self.tab[(mat + self.nbmat * (col - 1) + self.nbmtcl * (row - 1)) as usize]
    }

    /// OCCT ChangeValue(Row, Col, Mat) (BlendFunc_Tensor.lxx L24-27) —
    /// 1-based indexing, write access.
    #[inline]
    pub fn change_value(&mut self, row: i32, col: i32, mat: i32) -> &mut f64 {
        &mut self.tab[(mat + self.nbmat * (col - 1) + self.nbmtcl * (row - 1)) as usize]
    }

    /// OCCT Multiply(Right, M) (BlendFunc_Tensor.cxx L36-52).
    pub fn multiply(&self, right: &[f64], product: &mut [Vec<f64>]) {
        for i in 1..=self.nbrow {
            for j in 1..=self.nbcol {
                let mut somme = 0.0;
                for k in 1..=self.nbmat {
                    somme += self.value(i, j, k) * right[(k - 1) as usize];
                }
                product[(i - 1) as usize][(j - 1) as usize] = somme;
            }
        }
    }
}

/// OCCT GeomFill::Knots(TConv, TKnots) (GeomFill.cxx L268-286) — translated
/// here as a private helper (the GeomFill package is not ported yet; only the
/// statics consumed by BlendFunc_ConstRad are translated, 1:1).
pub(crate) fn geomfill_knots(t_conv: ConvertParameterisationType, tknots: &mut [f64]) {
    if (t_conv != ConvertParameterisationType::QuasiAngular)
        && (t_conv != ConvertParameterisationType::Polynomial)
    {
        let mut val = 0.0;
        for k in tknots.iter_mut() {
            *k = val;
            val += 1.0;
        }
    } else {
        tknots[0] = 0.0;
        tknots[1] = 1.0;
    }
}

/// OCCT GeomFill::Mults(TConv, TMults) (GeomFill.cxx L289-317).
pub(crate) fn geomfill_mults(t_conv: ConvertParameterisationType, tmults: &mut [i32]) {
    match t_conv {
        ConvertParameterisationType::QuasiAngular => {
            tmults[0] = 7;
            tmults[1] = 7;
        }
        ConvertParameterisationType::Polynomial => {
            tmults[0] = 8;
            tmults[1] = 8;
        }
        _ => {
            // Cas rational classique (GeomFill.cxx L306-315).
            let low = 0usize; // OCCT: TMults.Lower()
            let upp = tmults.len() - 1; // OCCT: TMults.Upper()
            tmults[low] = 3;
            for i in (low + 1)..=upp.saturating_sub(1) {
                tmults[i] = 2;
            }
            tmults[upp] = 3;
        }
    }
}

/// OCCT GeomFill::GetTolerance(TConv, AngleMin, Radius, AngularTol,
/// SpatialTol) (GeomFill.cxx L323-345).
pub(crate) fn geomfill_get_tolerance(
    t_conv: ConvertParameterisationType,
    angle_min: f64,
    radius: f64,
    angular_tol: f64,
    spatial_tol: f64,
) -> f64 {
    let tconv = match t_conv {
        ConvertParameterisationType::TgtThetaOver2 => ConvertParameterisation::TgtThetaOver2,
        ConvertParameterisationType::TgtThetaOver2_1 => ConvertParameterisation::TgtThetaOver2_1,
        ConvertParameterisationType::TgtThetaOver2_2 => ConvertParameterisation::TgtThetaOver2_2,
        ConvertParameterisationType::TgtThetaOver2_3 => ConvertParameterisation::TgtThetaOver2_3,
        ConvertParameterisationType::TgtThetaOver2_4 => ConvertParameterisation::TgtThetaOver2_4,
        ConvertParameterisationType::QuasiAngular => ConvertParameterisation::QuasiAngular,
        ConvertParameterisationType::RationalC1 => ConvertParameterisation::RationalC1,
        ConvertParameterisationType::Polynomial => ConvertParameterisation::Polynomial,
    };
    // OCCT: gp_Circ C(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir::Z), Radius).
    let circle = Circle3::new(DVec3::ZERO, DVec3::Z, radius);
    // OCCT: Sect = Geom_TrimmedCurve(popCircle, 0., max(AngleMin, 0.02)) —
    // 0.02 est proche d'1 degree, en desous on ne se preocupe pas de la
    // tngence afin d'eviter des tolerances d'approximation tendant vers 0 !
    let sect = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
        curve: Box::new(Curve3::Circle(circle)),
        first: 0.0,
        last: angle_min.max(0.02),
    });
    // OCCT: CtoBspl = GeomConvert::CurveToBSplineCurve(Sect, TConv).
    let ctobspl = geom_convert_curve_to_bspline_curve(&sect, tconv);
    // OCCT: Dist = CtoBspl->Pole(1).Distance(CtoBspl->Pole(2)) + SpatialTol.
    let dist = ctobspl.control_points[0].distance(ctobspl.control_points[1]) + spatial_tol;
    dist * angular_tol / 2.0
}

/// GAP (plan 0.6): GeomFill::GetCircle (TKGeomAlgo/GeomFill/GeomFill.cxx
/// L350-1900, five overloads) is a pending dependency of the ConstRad
/// sections; it needs its own GeomFill alignment batch.
pub(crate) fn geomfill_get_circle_pending() -> ! {
    unimplemented!("GeomFill::GetCircle pending in rcad (plan 0.6 kernel gap)")
}

/// OCCT normalizeAngle (ElCLib.cxx L43-72) — normalize angle to [0, 2*PI]
/// range, with special handling for values very close to zero to avoid
/// discontinuity.
fn elclib_normalize_angle(angle: &mut f64) {
    const PIPI: f64 = std::f64::consts::PI + std::f64::consts::PI;
    // OCCT: NEGATIVE_RESOLUTION = -Precision::Computational().
    const NEGATIVE_RESOLUTION: f64 = -f64::EPSILON;
    while *angle < NEGATIVE_RESOLUTION {
        *angle += PIPI;
    }
    // Only normalize angles strictly greater than 2*PI (with small tolerance)
    // to preserve the closing seam value of exactly 2*PI.
    while *angle > PIPI * (1.0 + GP_RESOLUTION) {
        *angle -= PIPI;
    }
    if *angle < 0.0 {
        *angle = 0.0;
    }
}

/// OCCT gp_Dir::AngleWithRef(Other, Vref) (gp_Dir.cxx L55-80) — the angle
/// between two directions, signed by Vref.
fn gp_dir_angle_with_ref(dir: DVec3, other: DVec3, vref: DVec3) -> f64 {
    let xyz = dir.cross(other);
    let cosinus = dir.dot(other);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT ElCLib::CircleParameter(Pos, P) (ElCLib.cxx L1199-1222) — the
/// parameter of P on the circle frame (1-based gp_Ax2 = the Circle3 fields).
pub(crate) fn elclib_circle_parameter(c: &Circle3, p: DVec3) -> f64 {
    let a_vec = p - c.center;
    // coinciding points -> infinite number of parameters
    if a_vec.length_squared() < GP_RESOLUTION {
        return 0.0;
    }

    let dir = c.normal;
    // Project vector on circle's plane.
    // OCCT: aVProj = dir.XYZ().CrossCrossed(aVec.XYZ(), dir.XYZ())
    //       (gp_XYZ::CrossCrossed(A, B) = this.Crossed(A.Crossed(B))).
    let a_vproj = dir.cross(a_vec.cross(dir));

    if a_vproj.length_squared() < GP_RESOLUTION {
        return 0.0;
    }

    // Angle between X direction and projected vector.
    let mut teta = gp_dir_angle_with_ref(c.x_dir, a_vproj, dir);
    elclib_normalize_angle(&mut teta);
    teta
}

/// OCCT BlendFunc_ConstRad — function for a constant-radius rolling-ball
/// blending surface (BlendFunc_ConstRad.hxx L27).
pub struct BlendFuncConstRad<'a> {
    // OCCT BlendFunc_ConstRad.hxx fields (L185-278).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tcurv — a handle aliasing curv until
    // Set(First, Last) trims it.  The rcad port owns the trimmed copy
    // (architecture mapping: Handle copy -> owned enum + accessor).
    pub(crate) tcurv: Option<Curve3>,
    pub(crate) pts1: DVec3,
    pub(crate) pts2: DVec3,
    pub(crate) istangent: bool,
    pub(crate) tg1: DVec3,
    pub(crate) tg12d: DVec2,
    pub(crate) tg2: DVec3,
    pub(crate) tg22d: DVec2,
    pub(crate) param: f64,
    pub(crate) ray1: f64,
    pub(crate) ray2: f64,
    pub(crate) choix: i32,
    pub(crate) my_x_order: i32,
    pub(crate) my_t_order: i32,
    pub(crate) xval: [f64; 4],
    pub(crate) tval: f64,
    // OCCT ComputeValues `static` locals (ConstRad.cxx L143-146) promoted to
    // fields: the guide point / derivatives and 1/normtg are a cache keyed on
    // tval / myTOrder (refreshed only when !t_OK).
    pub(crate) ptgui: DVec3,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) invnormtg: f64,
    pub(crate) d1u1: DVec3,
    pub(crate) d1u2: DVec3,
    pub(crate) d1v1: DVec3,
    pub(crate) d1v2: DVec3,
    pub(crate) d2u1: DVec3,
    pub(crate) d2v1: DVec3,
    pub(crate) d2uv1: DVec3,
    pub(crate) d2u2: DVec3,
    pub(crate) d2v2: DVec3,
    pub(crate) d2uv2: DVec3,
    pub(crate) dn1w: DVec3,
    pub(crate) dn2w: DVec3,
    pub(crate) d2n1w: DVec3,
    pub(crate) d2n2w: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) nsurf1: DVec3,
    pub(crate) dns1u1: DVec3,
    pub(crate) dns1u2: DVec3,
    pub(crate) dns1v1: DVec3,
    pub(crate) dns1v2: DVec3,
    pub(crate) nsurf2: DVec3,
    pub(crate) dnplan: DVec3,
    pub(crate) d2nplan: DVec3,
    // OCCT BlendFunc_ConstRad.hxx L240-241 — declared but never accessed by
    // the ConstRad bodies (vestigial fields, kept for structure parity).
    #[allow(dead_code)]
    pub(crate) dnsurf1: DVec3,
    #[allow(dead_code)]
    pub(crate) dnsurf2: DVec3,
    pub(crate) dndu1: DVec3,
    pub(crate) dndu2: DVec3,
    pub(crate) dndv1: DVec3,
    pub(crate) dndv2: DVec3,
    pub(crate) d2ndu1: DVec3,
    pub(crate) d2ndu2: DVec3,
    pub(crate) d2ndv1: DVec3,
    pub(crate) d2ndv2: DVec3,
    pub(crate) d2nduv1: DVec3,
    pub(crate) d2nduv2: DVec3,
    pub(crate) d2ndtu1: DVec3,
    pub(crate) d2ndtu2: DVec3,
    pub(crate) d2ndtv1: DVec3,
    pub(crate) d2ndtv2: DVec3,
    pub(crate) e: [f64; 4],
    pub(crate) dedx: Vec<Vec<f64>>,
    pub(crate) dedt: [f64; 4],
    pub(crate) d2edx2: BlendFuncTensor,
    pub(crate) d2edxdt: Vec<Vec<f64>>,
    pub(crate) d2edt2: [f64; 4],
    pub(crate) maxang: f64,
    pub(crate) minang: f64,
    pub(crate) distmin: f64,
    pub(crate) my_s_shape: BlendFuncSectionShape,
    pub(crate) my_t_conv: ConvertParameterisationType,
}

/// OCCT BlendFunc_ConstRadInv — inverse of the constant-radius function used
/// to find a solution on a restriction of one of the surfaces
/// (BlendFunc_ConstRadInv.hxx L26).  The vector X is t, w, U, V.  Methods
/// live in [`super::brep_blend_func_consrad_c`].
pub struct BlendFuncConstRadInv<'a> {
    // OCCT BlendFunc_ConstRadInv.hxx fields (L68-76).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) csurf: Option<&'a Curve2d>,
    pub(crate) ray1: f64,
    pub(crate) ray2: f64,
    pub(crate) choix: i32,
    pub(crate) first: bool,
}

impl<'a> BlendFuncConstRad<'a> {
    /// Architecture mapping: OCCT `tcurv` is a handle aliasing `curv` until
    /// Set(First, Last) replaces it with a trimmed copy (ConstRad.cxx L793).
    #[inline]
    pub(crate) fn tcurv(&self) -> &Curve3 {
        self.tcurv.as_ref().unwrap_or(self.curv)
    }

    /// OCCT BlendFunc_ConstRad(S1, S2, C) (BlendFunc_ConstRad.cxx L45-75).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncConstRad {
            surf1: s1,
            surf2: s2,
            curv: c,
            tcurv: None, // OCCT: tcurv(C) — alias, see tcurv().
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            istangent: true,
            tg1: DVec3::ZERO,
            tg12d: DVec2::ZERO,
            tg2: DVec3::ZERO,
            tg22d: DVec2::ZERO,
            param: 0.0,
            ray1: 0.0,
            ray2: 0.0,
            choix: 0,
            my_x_order: -1,
            my_t_order: -1,
            xval: [-9.876e100; 4], // OCCT: xval.Init(-9.876e100)
            tval: -9.876e100,      // OCCT: tval = -9.876e100
            ptgui: DVec3::ZERO,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            invnormtg: 0.0,
            d1u1: DVec3::ZERO,
            d1u2: DVec3::ZERO,
            d1v1: DVec3::ZERO,
            d1v2: DVec3::ZERO,
            d2u1: DVec3::ZERO,
            d2v1: DVec3::ZERO,
            d2uv1: DVec3::ZERO,
            d2u2: DVec3::ZERO,
            d2v2: DVec3::ZERO,
            d2uv2: DVec3::ZERO,
            dn1w: DVec3::ZERO,
            dn2w: DVec3::ZERO,
            d2n1w: DVec3::ZERO,
            d2n2w: DVec3::ZERO,
            nplan: DVec3::ZERO,
            nsurf1: DVec3::ZERO,
            dns1u1: DVec3::ZERO,
            dns1u2: DVec3::ZERO,
            dns1v1: DVec3::ZERO,
            dns1v2: DVec3::ZERO,
            nsurf2: DVec3::ZERO,
            dnplan: DVec3::ZERO,
            d2nplan: DVec3::ZERO,
            dnsurf1: DVec3::ZERO,
            dnsurf2: DVec3::ZERO,
            dndu1: DVec3::ZERO,
            dndu2: DVec3::ZERO,
            dndv1: DVec3::ZERO,
            dndv2: DVec3::ZERO,
            d2ndu1: DVec3::ZERO,
            d2ndu2: DVec3::ZERO,
            d2ndv1: DVec3::ZERO,
            d2ndv2: DVec3::ZERO,
            d2nduv1: DVec3::ZERO,
            d2nduv2: DVec3::ZERO,
            d2ndtu1: DVec3::ZERO,
            d2ndtu2: DVec3::ZERO,
            d2ndtv1: DVec3::ZERO,
            d2ndtv2: DVec3::ZERO,
            e: [0.0; 4],
            dedx: vec![vec![0.0; 4]; 4],
            dedt: [0.0; 4],
            d2edx2: BlendFuncTensor::new(4, 4, 4),
            d2edxdt: vec![vec![0.0; 4]; 4],
            d2edt2: [0.0; 4],
            maxang: -f64::MAX, // OCCT: RealFirst()
            minang: f64::MAX,  // OCCT: RealLast()
            distmin: f64::MAX, // OCCT: RealLast()
            my_s_shape: BlendFuncSectionShape::Rational,
            my_t_conv: ConvertParameterisationType::TgtThetaOver2,
        }
    }

    /// OCCT NbEquations() (BlendFunc_ConstRad.cxx L79-82) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT Set(Radius, Choix) (BlendFunc_ConstRad.cxx L86-118) — inits the
    /// value of radius, and the "quadrant".
    pub fn set(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.ray1 = -radius;
                self.ray2 = -radius;
            }
            3 | 4 => {
                self.ray1 = radius;
                self.ray2 = -radius;
            }
            5 | 6 => {
                self.ray1 = radius;
                self.ray2 = radius;
            }
            7 | 8 => {
                self.ray1 = -radius;
                self.ray2 = radius;
            }
            _ => {
                self.ray1 = -radius;
                self.ray2 = -radius;
            }
        }
    }

    /// OCCT Set(TypeSection) (BlendFunc_ConstRad.cxx L122-125) — sets the
    /// type of section generation for the approximations.
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_s_shape = type_section;
    }

    /// OCCT Set(Param) (BlendFunc_ConstRad.cxx L780-783).
    pub fn set_param(&mut self, param: f64) {
        self.param = param;
    }

    /// OCCT Set(First, Last) (BlendFunc_ConstRad.cxx L791-794) — segmente la
    /// courbe a sa partie utile.
    pub fn set_interval(&mut self, first: f64, last: f64) {
        // OCCT: tcurv = curv->Trim(First, Last, 1.e-12);
        self.tcurv = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(self.curv.clone()),
            first,
            last,
        }));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_ConstRad.cxx L798-804).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        tolerance[0] = self.surf1.u_resolution(tol);
        tolerance[1] = self.surf1.v_resolution(tol);
        tolerance[2] = self.surf2.u_resolution(tol);
        tolerance[3] = self.surf2.v_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_ConstRad.cxx L808-828).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.surf1.default_domain()[0]; // FirstUParameter
        inf_bound[1] = self.surf1.default_domain()[2]; // FirstVParameter
        inf_bound[2] = self.surf2.default_domain()[0];
        inf_bound[3] = self.surf2.default_domain()[2];
        sup_bound[0] = self.surf1.default_domain()[1]; // LastUParameter
        sup_bound[1] = self.surf1.default_domain()[3]; // LastVParameter
        sup_bound[2] = self.surf2.default_domain()[1];
        sup_bound[3] = self.surf2.default_domain()[3];

        for i in 0..4 {
            if !is_infinite_value(inf_bound[i]) && !is_infinite_value(sup_bound[i]) {
                let range = sup_bound[i] - inf_bound[i];
                inf_bound[i] -= range;
                sup_bound[i] += range;
            }
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ConstRad.cxx L832-963).
    #[allow(unreachable_code)] // the math_SVD fallback is a pending kernel gap
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let ok;

        ok = self.compute_values(sol, 1, true, self.param);

        if self.e[0].abs() <= tol
            && self.e[1] * self.e[1] + self.e[2] * self.e[2] + self.e[3] * self.e[3] <= tol * tol
        {
            // ns1, ns2 and  np are copied locally to avoid crushing the fields !
            let mut ns1 = self.nsurf1;
            let mut ns2 = self.nsurf2;
            let np = self.nplan;

            let mut norm = self.nplan.cross(ns1).length();
            if norm < EPS {
                norm = 1.0; // Unsatisfactory, but it is not necessary to stop
            }
            // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) / norm, nplan, -1. / norm, ns1);
            ns1 = (self.nplan.dot(ns1) / norm) * self.nplan + (-1.0 / norm) * ns1;

            let mut norm = self.nplan.cross(ns2).length();
            if norm < EPS {
                norm = 1.0; // Unsatisfactory, but it is not necessary to stop
            }
            ns2 = (self.nplan.dot(ns2) / norm) * self.nplan + (-1.0 / norm) * ns2;

            // OCCT: math_Vector controle(1, 4), solution(1, 4), tolerances(1, 4);
            //       GetTolerance(tolerances, Tol);
            let mut tolerances = [0.0f64; 4];
            self.get_tolerance(&mut tolerances, tol);

            self.istangent = true;
            // OCCT: math_Gauss Resol(DEDX, maxpiv).
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, self.dedx[r - 1][c - 1]);
                }
            }
            // OCCT: math_Vector solution(1, 4) — the solver output shared by
            // the Gauss and SVD branches.
            let mut solution = [0.0f64; 4];
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                // OCCT: Resol.Solve(-DEDT, solution);
                let mut x = VecD::new(4);
                for i in 1..=4 {
                    x.set(i, -self.dedt[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=4 {
                    solution[i - 1] = x.get(i);
                }
                self.istangent = false;
                // OCCT: controle = DEDT.Added(DEDX.Multiplied(solution));
                let mut controle = [0.0f64; 4];
                for i in 1..=4 {
                    let mut somme = 0.0;
                    for j in 1..=4 {
                        somme += self.dedx[i - 1][j - 1] * solution[j - 1];
                    }
                    controle[i - 1] = self.dedt[i - 1] + somme;
                }
                if controle[0].abs() > tolerances[0]
                    || controle[1].abs() > tolerances[1]
                    || controle[2].abs() > tolerances[2]
                    || controle[3].abs() > tolerances[3]
                {
                    self.istangent = true;
                }
            }

            if self.istangent {
                // OCCT L880-897: math_SVD SingRS(DEDX); if (SingRS.IsDone()) {
                // SingRS.Solve(-DEDT, solution, 1.e-6); controle = ...; }.
                // GAP (plan 0.6): rcad-kernel exposes no public math_SVD yet
                // (only private svd helpers in function_set_root.rs); the SVD
                // fallback is pending kernel support and keeps the tangency
                // state like the OCCT !IsDone() path.
            }

            if !self.istangent {
                // OCCT: tg1.SetLinearForm(solution(1), d1u1, solution(2), d1v1);
                self.tg1 = solution[0] * self.d1u1 + solution[1] * self.d1v1;
                self.tg2 = solution[2] * self.d1u2 + solution[3] * self.d1v2;
                self.tg12d = DVec2::new(solution[0], solution[1]);
                self.tg22d = DVec2::new(solution[2], solution[3]);
            }

            // update of maxang

            if self.ray1 > 0.0 {
                ns1 = -ns1;
            }
            if self.ray2 > 0.0 {
                ns2 = -ns2;
            }
            let mut cosa = ns1.dot(ns2);
            let mut sina = np.dot(ns1.cross(ns2));
            if self.choix % 2 != 0 {
                sina = -sina; // nplan is changed in -nplan
            }

            if cosa > 1.0 {
                cosa = 1.0;
                sina = 0.0;
            }
            let mut angle = cosa.acos();

            // Reframing on  ]-pi/2, 3pi/2]
            if sina < 0.0 {
                if cosa > 0.0 {
                    angle = -angle;
                } else {
                    angle = 2.0 * std::f64::consts::PI - angle;
                }
            }

            if angle.abs() > self.maxang {
                self.maxang = angle.abs();
            }
            if angle.abs() < self.minang {
                self.minang = angle.abs();
            }
            self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

            return ok;
        }
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BlendFunc_ConstRad.cxx L967-970).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT Value(X, F) (BlendFunc_ConstRad.cxx L974-979).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let ok = self.compute_values(x, 0, false, 0.0);
        f.copy_from_slice(&self.e);
        ok
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ConstRad.cxx L983-988).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let ok = self.compute_values(x, 1, false, 0.0);
        // OCCT: D = DEDX;
        for (drow, xrow) in d.iter_mut().zip(self.dedx.iter()) {
            drow.clone_from(xrow);
        }
        ok
    }

    /// OCCT Values(X, F, D) (BlendFunc_ConstRad.cxx L992-998).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        let ok = self.compute_values(x, 1, false, 0.0);
        f.copy_from_slice(&self.e);
        // OCCT: D = DEDX;
        for (drow, xrow) in d.iter_mut().zip(self.dedx.iter()) {
            drow.clone_from(xrow);
        }
        ok
    }

    /// OCCT PointOnS1() (BlendFunc_ConstRad.cxx L1002-1005).
    pub fn point_on_s1(&self) -> DVec3 {
        self.pts1
    }

    /// OCCT PointOnS2() (BlendFunc_ConstRad.cxx L1009-1012).
    pub fn point_on_s2(&self) -> DVec3 {
        self.pts2
    }

    /// OCCT IsTangencyPoint() (BlendFunc_ConstRad.cxx L1016-1019).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS1() (BlendFunc_ConstRad.cxx L1023-1030).
    pub fn tangent_on_s1(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ConstRad::TangentOnS1");
        }
        self.tg1
    }

    /// OCCT TangentOnS2() (BlendFunc_ConstRad.cxx L1034-1041).
    pub fn tangent_on_s2(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ConstRad::TangentOnS2");
        }
        self.tg2
    }

    /// OCCT Tangent2dOnS1() (BlendFunc_ConstRad.cxx L1045-1052).
    pub fn tangent_2d_on_s1(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ConstRad::Tangent2dOnS1");
        }
        self.tg12d
    }

    /// OCCT Tangent2dOnS2() (BlendFunc_ConstRad.cxx L1057-1063).
    pub fn tangent_2d_on_s2(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ConstRad::Tangent2dOnS2");
        }
        self.tg22d
    }

    /// OCCT Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)
    /// (BlendFunc_ConstRad.cxx L1067-1115).
    #[allow(clippy::too_many_arguments)]
    pub fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_f: &mut DVec3,
        tg_l: &mut DVec3,
        nm_f: &mut DVec3,
        nm_l: &mut DVec3,
    ) {
        let ns1: DVec3;

        if (u1 != self.xval[0]) || (v1 != self.xval[1]) || (u2 != self.xval[2]) || (v2 != self.xval[3])
        {
            // OCCT: surf1->D1(U1, V1, bid, d1u, d1v); NmF = ns1 = d1u.Crossed(d1v);
            let (_, d1u_1, d1v_1) = self.surf1.derivatives(u1, v1);
            let s = d1u_1.cross(d1v_1);
            *nm_f = s;
            ns1 = s;
            // OCCT: surf2->D1(U2, V2, bid, d1u, d1v); NmL = d1u.Crossed(d1v);
            let (_, d1u_2, d1v_2) = self.surf2.derivatives(u2, v2);
            *nm_l = d1u_2.cross(d1v_2);
        } else {
            *nm_f = self.nsurf1;
            ns1 = self.nsurf1;
            *nm_l = self.nsurf2;
        }

        // OCCT: invnorm1 = nplan.Crossed(ns1).Magnitude();
        //       if (invnorm1 < Eps) invnorm1 = 1; else invnorm1 = 1. / invnorm1;
        let invnorm1 = {
            let n = self.nplan.cross(ns1).length();
            if n < EPS {
                1.0
            } else {
                1.0 / n
            }
        };

        // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) * invnorm1, nplan, -invnorm1, ns1);
        let ns1 = (self.nplan.dot(ns1) * invnorm1) * self.nplan + (-invnorm1) * ns1;
        // OCCT: Center.SetXYZ(pts1.XYZ() + ray1 * ns1.XYZ());
        let center = self.pts1 + self.ray1 * ns1;

        // OCCT: TgF = nplan.Crossed(gp_Vec(Center, pts1));
        *tg_f = self.nplan.cross(self.pts1 - center);
        *tg_l = self.nplan.cross(self.pts2 - center);
        if self.choix % 2 == 1 {
            *tg_f = -*tg_f;
            *tg_l = -*tg_l;
        }
    }

    /// OCCT TwistOnS1() (BlendFunc_ConstRad.cxx L1119-1126).
    pub fn twist_on_s1(&self) -> bool {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ConstRad::TwistOnS1");
        }
        self.tg1.dot(self.nplan) < 0.0
    }

    /// OCCT TwistOnS2() (BlendFunc_ConstRad.cxx L1130-1137).
    pub fn twist_on_s2(&self) -> bool {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ConstRad::TwistOnS2");
        }
        self.tg2.dot(self.nplan) < 0.0
    }

    /// OCCT IsRational() (BlendFunc_ConstRad.cxx L1203-1206).
    pub fn is_rational(&self) -> bool {
        self.my_s_shape == BlendFuncSectionShape::Rational
            || self.my_s_shape == BlendFuncSectionShape::QuasiAngular
    }

    /// OCCT GetSectionSize() (BlendFunc_ConstRad.cxx L1210-1213).
    pub fn get_section_size(&self) -> f64 {
        self.maxang * self.ray1.abs()
    }

    /// OCCT GetMinimalWeight(Weights) (BlendFunc_ConstRad.cxx L1217-1221) —
    /// it is supposed that it does not depend on the Radius!
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        blend_func_get_minimal_weights(self.my_s_shape, self.my_t_conv, self.minang, self.maxang, weigths);
    }

    /// OCCT NbIntervals(S) (BlendFunc_ConstRad.cxx L1225-1228).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.curv.nb_intervals(blend_func_next_shape(s))
    }

    /// OCCT Intervals(T, S) (BlendFunc_ConstRad.cxx L1232-1235).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let mut intervals = Vec::new();
        self.curv.intervals(&mut intervals, blend_func_next_shape(s));
        for (dst, src) in t.iter_mut().zip(intervals) {
            *dst = src;
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BlendFunc_ConstRad.cxx L1239-1243).
    pub fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        *nb_poles_2d = 2;
        blend_func_get_shape(
            self.my_s_shape,
            self.maxang,
            nb_poles,
            nb_knots,
            degree,
            &mut self.my_t_conv,
        );
    }

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)
    /// (BlendFunc_ConstRad.cxx L1249-1262) — tolerances used for
    /// approximations.
    pub fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        let low = 0usize; // OCCT: Tol3d.Lower()
        let up = tol3d.len() - 1; // OCCT: Tol3d.Upper()
        let tol = geomfill_get_tolerance(self.my_t_conv, self.minang, self.ray1.abs(), angle_tol, surf_tol);
        for v in tol1d.iter_mut() {
            *v = surf_tol;
        }
        for v in tol3d.iter_mut() {
            *v = surf_tol;
        }
        tol3d[low + 1] = tol.min(surf_tol);
        tol3d[up - 1] = tol.min(surf_tol);
        tol3d[low] = tol.min(bound_tol);
        tol3d[up] = tol.min(bound_tol);
    }

    /// OCCT Knots(TKnots) (BlendFunc_ConstRad.cxx L1266-1269).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        geomfill_knots(self.my_t_conv, tknots);
    }

    /// OCCT Mults(TMults) (BlendFunc_ConstRad.cxx L1273-1276).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        geomfill_mults(self.my_t_conv, tmults);
    }

    /// OCCT AxeRot(Prm) (BlendFunc_ConstRad.cxx L1907-1939).
    pub fn axe_rot(&mut self, prm: f64) -> Ax1 {
        let mut axrot = Ax1::new(DVec3::ZERO, DVec3::ZERO);

        // OCCT: curv->D2(Prm, ptgui, d1gui, d2gui);
        let ptgui = self.curv.point_at(prm);
        let d1gui = self.curv.derivative_at(prm);
        let d2gui = self.curv.derivative2_at(prm);

        let normtg = d1gui.length();
        let np = d1gui.normalize_or_zero();
        // OCCT: dnp.SetLinearForm(1. / normtg, d2gui, -1. / normtg * (np.Dot(d2gui)), np);
        let dnp = (1.0 / normtg) * d2gui + (-(1.0 / normtg) * np.dot(d2gui)) * np;

        let dirax = np.cross(dnp);
        if dirax.length() >= GP_RESOLUTION {
            axrot.set_direction(dirax);
        } else {
            axrot.set_direction(np); // To avoid stop
        }
        let oriax = if dnp.length() >= GP_RESOLUTION {
            ptgui + (normtg / dnp.length()) * dnp.normalize_or_zero()
        } else {
            ptgui
        };
        axrot.set_location(oriax);
        axrot
    }

    /// OCCT Resolution(IC2d, Tol, TolU, TolV) (BlendFunc_ConstRad.cxx
    /// L1941-1956).
    pub fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic_2d == 1 {
            *tol_u = self.surf1.u_resolution(tol);
            *tol_v = self.surf1.v_resolution(tol);
        } else {
            *tol_u = self.surf2.u_resolution(tol);
            *tol_v = self.surf2.v_resolution(tol);
        }
    }
}

impl<'a> FunctionSetWithDerivatives for BlendFuncConstRad<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_Function::NbVariables (returns 4).
        BlendFunction::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncConstRad::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncConstRad::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncConstRad::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncConstRad::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendFuncConstRad<'a> {
    fn set_param(&mut self, param: f64) {
        BlendFuncConstRad::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendFuncConstRad::set_interval(self, first, last)
    }

    // OCCT Blend_Function.cxx L24-32 — Pnt1/Pnt2 delegate to
    // PointOnS1/PointOnS2.
    fn pnt1(&self) -> DVec3 {
        BlendFunction::pnt1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendFunction::pnt2(self)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendFuncConstRad::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncConstRad::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncConstRad::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendFuncConstRad::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendFuncConstRad::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendFuncConstRad::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendFuncConstRad::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendFuncConstRad::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendFuncConstRad::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendFuncConstRad::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendFuncConstRad::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendFuncConstRad::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendFuncConstRad::mults(self, tmults)
    }

    fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        BlendFuncConstRad::section_d1(
            self, p, poles, d_poles, poles_2d, d_poles_2d, weigths, d_weigths,
        )
    }

    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        BlendFuncConstRad::section_simple(self, p, poles, poles_2d, weigths)
    }

    fn section_d2(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        d2_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        d2_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
        d2_weigths: &mut [f64],
    ) -> bool {
        BlendFuncConstRad::section_d2(
            self,
            p,
            poles,
            d_poles,
            d2_poles,
            poles_2d,
            d_poles_2d,
            d2_poles_2d,
            weigths,
            d_weigths,
            d2_weigths,
        )
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendFuncConstRad::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendFunction for BlendFuncConstRad<'a> {
    fn point_on_s1(&self) -> DVec3 {
        BlendFuncConstRad::point_on_s1(self)
    }

    fn point_on_s2(&self) -> DVec3 {
        BlendFuncConstRad::point_on_s2(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendFuncConstRad::is_tangency_point(self)
    }

    fn tangent_on_s1(&self) -> DVec3 {
        BlendFuncConstRad::tangent_on_s1(self)
    }

    fn tangent_2d_on_s1(&self) -> DVec2 {
        BlendFuncConstRad::tangent_2d_on_s1(self)
    }

    fn tangent_on_s2(&self) -> DVec3 {
        BlendFuncConstRad::tangent_on_s2(self)
    }

    fn tangent_2d_on_s2(&self) -> DVec2 {
        BlendFuncConstRad::tangent_2d_on_s2(self)
    }

    fn twist_on_s1(&self) -> bool {
        BlendFuncConstRad::twist_on_s1(self)
    }

    fn twist_on_s2(&self) -> bool {
        BlendFuncConstRad::twist_on_s2(self)
    }

    fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut DVec3,
        tg_last: &mut DVec3,
        norm_first: &mut DVec3,
        norm_last: &mut DVec3,
    ) {
        BlendFuncConstRad::tangent(self, u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last)
    }
}
