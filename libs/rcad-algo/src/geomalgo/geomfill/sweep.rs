//! OCCT GeomFill_Sweep (TKGeomAlgo/GeomFill) — 1:1 port of GeomFill_Sweep.hxx
//! (members) + GeomFill_Sweep.cxx (whole file L70-1249).
//!
//! Architecture differences:
//! - `handle(GeomFill_LocationLaw) myLoc` / `handle(GeomFill_SectionLaw)
//!   mySec` map to `Rc<RefCell<...>>` (the laws are shared with the
//!   SweepFunction and later mutated in place through the same handles).
//! - GAP carrier: [`ApproxSweepApproximation`] (TKGeomBase/AppBlend,
//!   deriving AppBlend_AppSurf) — the sweep approximation engine is not
//!   translated; construction/perform keep the OCCT failure path.  The
//!   BuildAll result read-back is written literally against the carrier's
//!   typed accessors.
//! - GAP carrier: [`build_product`] — depends on GeomFill_LocFunction +
//!   AdvApprox_ApproxAFunction / AdvApprox_PrefAndRec (outside batch); the
//!   OCCT call site (Build L218) is commented out — the function is
//!   unreachable there.
//! - `gp_GTrsf` / `gp_Trsf::SetValues` map to [`gp_trsf_set_values`] (the
//!   null-determinant raise is the Err branch; the vectorial part is scaled
//!   by det^(1/3) and orthogonalized per gp_Trsf::Orthogonalize).
//! - `aSurf->VIso(Vmax)` in IsSweepParallelSpine: only the iso curve's first
//!   point is consumed — computed directly as S(Umin, Vmax) (same value).
//! - The torus/cone iso data consumed by BuildKPart is computed by the
//!   local ElSLib re-hosts ([`el_slib_torus_uiso_ax2`],
//!   [`cone_u_derivative_at_v`]) — pure math with OCCT anchors.

use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

use glam::{DAffine3, DVec2, DVec3};

use rcad_kernel::geom::{
    transform_curve, BSplineCurve2, BSplineSurface, ConicalSurface, Curve2d, Curve3, CurveEval,
    CylindricalSurface, Line2d, LinearExtrusionSurface, Plane, RevolutionSurface, SphericalSurface,
    Surface3, SurfaceEval, ToroidalSurface, TrimmedCurve2, TrimmedCurve3, TrimmedSurface,
};
use rcad_kernel::math::gp::{Ax2, Ax3, GP_RESOLUTION};
use rcad_kernel::math::GeomAbsShape;

use super::gp_mat::GpMat;
use super::location_law::LocationLaw;
use super::section_law::SectionLaw;
use super::sweep_function::SweepFunction;
use super::sweep_section_generator::{
    gp_vec_is_opposite, gp_vec_is_parallel,
};

/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-9;

/// OCCT GeomFill_ApproxStyle (GeomFill_ApproxStyle.hxx L18-22).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum GeomFillApproxStyle {
    GeomFill_Section,
    GeomFill_Location,
}

// ---------------------------------------------------------------------------
// Pure-math gp re-hosts
// ---------------------------------------------------------------------------

/// OCCT gp_Trsf::Orthogonalize (gp_Trsf.cxx L863-945) — Gram-Schmidt over
/// the columns then over the rows of the vectorial part.
fn trsf_orthogonalize(m: &GpMat) -> GpMat {
    let mut a_tm = *m;

    let mut v1 = gp_mat_column(&a_tm, 1);
    let mut v2 = gp_mat_column(&a_tm, 2);
    let mut v3 = gp_mat_column(&a_tm, 3);

    v1 = v1.normalize_or_zero();
    v2 -= v1 * (v2.dot(v1));
    v2 = v2.normalize_or_zero();
    v3 -= v1 * (v3.dot(v1)) + v2 * (v3.dot(v2));
    v3 = v3.normalize_or_zero();
    a_tm.set_cols(v1, v2, v3);

    let row = |t: &GpMat, i: usize| DVec3::new(t.mat[i - 1][0], t.mat[i - 1][1], t.mat[i - 1][2]);
    let mut r1 = row(&a_tm, 1);
    let mut r2 = row(&a_tm, 2);
    let mut r3 = row(&a_tm, 3);

    r1 = r1.normalize_or_zero();
    r2 -= r1 * (r2.dot(r1));
    r2 = r2.normalize_or_zero();
    r3 -= r1 * (r3.dot(r1)) + r2 * (r3.dot(r2));
    r3 = r3.normalize_or_zero();
    a_tm.mat[0] = [r1.x, r1.y, r1.z];
    a_tm.mat[1] = [r2.x, r2.y, r2.z];
    a_tm.mat[2] = [r3.x, r3.y, r3.z];

    a_tm
}

/// OCCT gp_Mat::Column(theCol) (gp_Mat.hxx) — the column read back as an XYZ
/// (pure math helper).
fn gp_mat_column(m: &GpMat, the_col: usize) -> DVec3 {
    DVec3::new(
        m.mat[0][the_col - 1],
        m.mat[1][the_col - 1],
        m.mat[2][the_col - 1],
    )
}

/// OCCT gp_Trsf::SetValues(a11..a34) (gp_Trsf.cxx L346-385) — the null
/// determinant raise (`|det| < gp::Resolution()`) is the Err branch; the
/// vectorial part is scaled by det^(1/3) and orthogonalized, loc = col4.
/// The rcad carrier is [`DAffine3`] (the scale folds into the matrix).
fn gp_trsf_set_values(a: &[f64; 12]) -> Result<DAffine3, ()> {
    // compute the determinant of the vectorial part.
    let m = GpMat::from_rows(a[0], a[1], a[2], a[4], a[5], a[6], a[8], a[9], a[10]);
    let mut s = m.mat[0][0] * (m.mat[1][1] * m.mat[2][2] - m.mat[1][2] * m.mat[2][1])
        - m.mat[0][1] * (m.mat[1][0] * m.mat[2][2] - m.mat[1][2] * m.mat[2][0])
        + m.mat[0][2] * (m.mat[1][0] * m.mat[2][1] - m.mat[1][1] * m.mat[2][0]);
    if s.abs() < GP_RESOLUTION {
        return Err(()); // "gp_Trsf::SetValues, null determinant"
    }
    if s > 0.0 {
        s = s.powf(1.0 / 3.0);
    } else {
        s = -((-s).powf(1.0 / 3.0));
    }
    let mut scaled = m.multiplied_scalar(1.0 / s);
    scaled = trsf_orthogonalize(&scaled);

    let mut trsf = DAffine3::from_mat3(glam::DMat3::from_cols(
        gp_mat_column(&scaled, 1),
        gp_mat_column(&scaled, 2),
        gp_mat_column(&scaled, 3),
    ));
    trsf.translation = DVec3::new(a[3], a[7], a[11]);
    Ok(trsf)
}

/// OCCT gp_Trsf::SetValues over a gp_GTrsf (SetVectorialPart(M) +
/// SetTranslationPart(V) then the 12-argument SetValues).
fn gp_trsf_set_values_from_parts(m: &GpMat, v: DVec3) -> Result<DAffine3, ()> {
    gp_trsf_set_values(&[
        m.mat[0][0], m.mat[0][1], m.mat[0][2], v.x,
        m.mat[1][0], m.mat[1][1], m.mat[1][2], v.y,
        m.mat[2][0], m.mat[2][1], m.mat[2][2], v.z,
    ])
}

/// OCCT gp_Vec::Angle (gp_XYZ::Angle) — the angle between two vectors
/// (null vectors raise in OCCT; the rcad form panics likewise).
fn gp_vec_angle(v1: DVec3, v2: DVec3) -> f64 {
    let an_norm = v1.length();
    let a_no_norm = v2.length();
    let mut value = v1.dot(v2) / (an_norm * a_no_norm);
    if value > 1.0 {
        value = 1.0;
    } else if value < -1.0 {
        value = -1.0;
    }
    value.acos()
}

/// OCCT ElSLib::TorusUIso (ElSLib.cxx L1751-1766) — the u-iso circle frame
/// of a torus at U, as the rcad Ax2 (location / direction / x_direction).
fn el_slib_torus_uiso_ax2(
    pos: &Ax3,
    major_radius: f64,
    minor_radius: f64,
    u: f64,
) -> (DVec3, DVec3, DVec3) {
    let dx = pos.x_direction;
    let dy = pos.y_direction;
    let dz = pos.direction();
    let cx = u.cos() * dx + u.sin() * dy;
    // gp_Ax2 axes(Pos.Location(), cx.Crossed(dz), cx); translated by
    // MajorRadius * cx.
    let location = pos.axis.location + major_radius * cx;
    let direction = cx.cross(dz).normalize_or_zero();
    let _ = minor_radius; // the minor radius is the circle radius
    (location, direction, cx)
}

/// OCCT Geom_ConicalSurface::VIso(V) + the D1-at-u=0 tangent consumed at the
/// BuildKPart call site — the iso circle center and the du/|du| direction in
/// the rcad cone layout (the rcad `apex` field is the reference point where
/// the radius is `radius`; the OCCT v parameter is measured from the true
/// apex).
fn cone_u_derivative_at_v(surf: &ConicalSurface, v: f64) -> DVec3 {
    // true apex = ref point - (radius / tan(A)) * axis.
    let a = surf.half_angle_rad;
    let apex = surf.apex - (surf.radius / a.tan()) * surf.axis;
    // center of the iso circle at V and its radius.
    let radius = v * a.sin();
    // D1(U=0, V) with respect to u: radius * (-sin(0) X + cos(0) Y)
    // = radius * Y with Y = axis ^ XDirection (right-handed frame).
    let y = surf.axis.cross(surf.ref_dir).normalize_or_zero();
    let _ = apex;
    radius * y
}

/// OCCT ElSLib::SphereParameters (ElSLib.cxx L1615-1646) — the (U, V) of P
/// on a sphere whose frame is `pos`.
fn el_slib_sphere_parameters(pos: &Ax3, p: DVec3) -> (f64, f64) {
    // gp_Trsf T; T.SetTransformation(Pos) — world -> local frame.
    let frame = ax3_frame(pos);
    let ploc = frame.inverse().transform_point3(p);
    let (x, y, z) = (ploc.x, ploc.y, ploc.z);
    let l = (x * x + y * y).sqrt();
    if l < GP_RESOLUTION {
        // point on axis Z of the sphere.
        if z > 0.0 {
            (0.0, PI * 0.5)
        } else {
            (0.0, -PI * 0.5)
        }
    } else {
        let v = (z / l).atan();
        let mut u = y.atan2(x);
        // OCCT normalizeAngle (ElSLib.cxx L42-59).
        while u < -GP_RESOLUTION {
            u += 2.0 * PI;
        }
        while u > 2.0 * PI * (1.0 + GP_RESOLUTION) {
            u -= 2.0 * PI;
        }
        if u < 0.0 {
            u = 0.0;
        }
        (u, v)
    }
}

/// The gp_Ax3 frame affine (local -> world).
fn ax3_frame(a: &Ax3) -> DAffine3 {
    let m = glam::DMat3::from_cols(a.x_direction, a.y_direction, a.axis.direction);
    let mut f = DAffine3::from_mat3(m);
    f.translation = a.axis.location;
    f
}

// ---------------------------------------------------------------------------
// GAP carrier
// ---------------------------------------------------------------------------

/// GAP carrier: OCCT Approx_SweepApproximation (TKGeomBase/AppBlend,
/// deriving AppBlend_AppSurf) — the sweep approximation engine is not
/// translated; construction and Perform keep the OCCT failure path.  The
/// accessors are typed for the literal BuildAll read-back.
pub struct ApproxSweepApproximation;

impl ApproxSweepApproximation {
    /// OCCT Approx_SweepApproximation(Func).
    pub fn new(_func: &dyn super::approx_sweep_function::ApproxSweepFunction) -> Self {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Approx_SweepApproximation::Perform(First, Last, Tol3d, BoundTol,
    /// Tol2d, TolAngular, Continuity, Degmax, Segmax).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        _first: f64,
        _last: f64,
        _tol3d: f64,
        _bound_tol: f64,
        _tol2d: f64,
        _tol_angular: f64,
        _continuity: GeomAbsShape,
        _degmax: i32,
        _segmax: i32,
    ) {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfShape(UDegree, VDegree, NbUPoles, NbVPoles, NbUKnots,
    /// NbVKnots).
    #[allow(clippy::too_many_arguments)]
    pub fn surf_shape(
        &self,
        _u_degree: &mut i32,
        _v_degree: &mut i32,
        _nb_u_poles: &mut i32,
        _nb_v_poles: &mut i32,
        _nb_u_knots: &mut i32,
        _nb_v_knots: &mut i32,
    ) {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Surface(Poles, Weights, UKnots, VKnots, UMults, VMults).
    #[allow(clippy::too_many_arguments)]
    pub fn surface(
        &self,
        _poles: &mut [Vec<DVec3>],
        _weights: &mut [Vec<f64>],
        _u_knots: &mut [f64],
        _v_knots: &mut [f64],
        _u_mults: &mut [i32],
        _v_mults: &mut [i32],
    ) {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT UDegree().
    pub fn u_degree(&self) -> i32 {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT VDegree().
    pub fn v_degree(&self) -> i32 {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT MaxErrorOnSurf().
    pub fn max_error_on_surf(&self) -> f64 {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Curve2dPoles(Index).
    pub fn curve2d_poles(&self, _index: i32) -> Vec<DVec2> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Curves2dKnots().
    pub fn curves2d_knots(&self) -> Vec<f64> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Curves2dMults().
    pub fn curves2d_mults(&self) -> Vec<i32> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Curves2dDegree().
    pub fn curves2d_degree(&self) -> i32 {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT Max2dError(Index).
    pub fn max2d_error(&self, _index: i32) -> f64 {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfPoles().
    pub fn surf_poles(&self) -> Vec<Vec<DVec3>> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfWeights().
    pub fn surf_weights(&self) -> Vec<Vec<f64>> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfUKnots().
    pub fn surf_u_knots(&self) -> Vec<f64> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfVKnots().
    pub fn surf_v_knots(&self) -> Vec<f64> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfUMults().
    pub fn surf_u_mults(&self) -> Vec<i32> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT SurfVMults().
    pub fn surf_v_mults(&self) -> Vec<i32> {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }

    /// OCCT TolReached( Tol3d, Tol2d).
    pub fn tol_reached(&self, _tol3d: &mut f64, _tol2d: &mut f64) {
        panic!(
            "GAP: Approx_SweepApproximation (TKGeomBase/AppBlend) is not translated — see file header"
        )
    }
}

/// GAP carrier: OCCT GeomConvert_ApproxSurface (TKGeomBase/GeomConvert) —
/// consumed by the BuildAll myForceApproxC1 branch.
pub struct GeomConvertApproxSurface;

impl GeomConvertApproxSurface {
    /// OCCT GeomConvert_ApproxSurface(Surface, Tol3d, UContinuity,
    /// VContinuity, MaxDegreeU, MaxDegreeV, MaxSegments, PrecisCode).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        _surface: &BSplineSurface,
        _tol3d: f64,
        _u_continuity: GeomAbsShape,
        _v_continuity: GeomAbsShape,
        _max_degree_u: i32,
        _max_degree_v: i32,
        _max_segments: i32,
        _precis_code: i32,
    ) -> Self {
        panic!(
            "GAP: GeomConvert_ApproxSurface (TKGeomBase/GeomConvert) is not translated — see file header"
        )
    }

    /// OCCT HasResult().
    pub fn has_result(&self) -> bool {
        panic!(
            "GAP: GeomConvert_ApproxSurface (TKGeomBase/GeomConvert) is not translated — see file header"
        )
    }

    /// OCCT Surface().
    pub fn surface(&self) -> BSplineSurface {
        panic!(
            "GAP: GeomConvert_ApproxSurface (TKGeomBase/GeomConvert) is not translated — see file header"
        )
    }
}

/// The OCCT knots/multiplicities arrays recovered from the rcad flat knot
/// vector (pure tool helper).
fn flat_knots_to_mults(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i == 0 || *k != knots[knots.len() - 1] {
            knots.push(*k);
            mults.push(1);
        } else {
            *mults.last_mut().expect("nonempty") += 1;
        }
    }
    (knots, mults)
}

/// OCCT Geom_BSplineSurface::IsCNv(N) — the v-direction continuity is at
/// least N (no interior v-knot with multiplicity > degree_v - N).
fn bspline_surface_is_cn_v(surf: &BSplineSurface, n: usize) -> bool {
    let (vknots, vmults) = flat_knots_to_mults(&surf.knots_v);
    if vknots.len() < 2 {
        return true;
    }
    for k in 1..vknots.len() - 1 {
        if vmults[k] as usize > surf.degree_v - n {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// The class
// ---------------------------------------------------------------------------

/// OCCT GeomFill_Sweep (GeomFill_Sweep.hxx private members L128-158).
pub struct Sweep {
    /// OCCT double First.
    first: f64,
    /// OCCT double Last.
    last: f64,
    /// OCCT double SFirst.
    s_first: f64,
    /// OCCT double SLast.
    s_last: f64,
    /// OCCT double Tol3d.
    tol3d: f64,
    /// OCCT double BoundTol.
    bound_tol: f64,
    /// OCCT double Tol2d.
    tol2d: f64,
    /// OCCT double TolAngular.
    tol_angular: f64,
    /// OCCT double SError.
    s_error: f64,
    /// OCCT bool myForceApproxC1.
    my_force_approx_c1: bool,
    /// OCCT handle(GeomFill_LocationLaw) myLoc.
    my_loc: Rc<RefCell<dyn LocationLaw>>,
    /// OCCT handle(GeomFill_SectionLaw) mySec.
    my_sec: Option<Rc<RefCell<dyn SectionLaw>>>,
    /// OCCT handle(Geom_Surface) mySurface.
    my_surface: Option<Surface3>,
    /// OCCT handle(NCollection_HArray1<handle<Geom2d_Curve>>) myCurve2d —
    /// the OCCT fills array indices (with the restriction holes filled by
    /// the iso-border fallback), hence the Option slots.
    my_curve2d: Option<Vec<Option<Curve2d>>>,
    /// OCCT handle(NCollection_HArray2<double>) CError (1..2, 1..N).
    c_error: Vec<f64>,
    /// OCCT bool done.
    done: bool,
    /// OCCT bool myExchUV.
    my_exch_uv: bool,
    /// OCCT bool isUReversed.
    is_u_reversed: bool,
    /// OCCT bool isVReversed.
    is_v_reversed: bool,
    /// OCCT bool myKPart.
    my_k_part: bool,
}

impl Sweep {
    /// OCCT GeomFill_Sweep::GeomFill_Sweep (L104-117).
    pub fn new(location: Rc<RefCell<dyn LocationLaw>>, with_kpart: bool) -> Self {
        let done = false;

        let my_loc = location;
        let my_k_part = with_kpart;
        // OCCT: SetTolerance(1.e-4) — the defaulted BoundTol/Tol2d/TolAngular.
        let (tol3d, bound_tol, tol2d, tol_angular) = (1.0e-4, 1.0, 1.0e-5, 1.0);
        let (mut first, mut last) = (0.0, 0.0);
        my_loc.borrow().get_domain(&mut first, &mut last);
        let s_first = 30.081996;
        let s_last = 30.081996;
        let s_error = f64::MAX; // RealLast()

        Sweep {
            first,
            last,
            s_first,
            s_last,
            tol3d,
            bound_tol,
            tol2d,
            tol_angular,
            s_error,
            my_force_approx_c1: false,
            my_loc,
            my_sec: None,
            my_surface: None,
            my_curve2d: None,
            c_error: Vec::new(),
            done,
            my_exch_uv: false,
            is_u_reversed: false,
            is_v_reversed: false,
            my_k_part,
        }
    }

    /// OCCT SetDomain (L121-130).
    pub fn set_domain(&mut self, loc_first: f64, loc_last: f64, section_first: f64, section_last: f64) {
        self.first = loc_first;
        self.last = loc_last;
        self.s_first = section_first;
        self.s_last = section_last;
    }

    /// OCCT SetTolerance (L134-143).
    pub fn set_tolerance_sweep(
        &mut self,
        tolerance3d: f64,
        bound_tolerance: f64,
        tolerance2d: f64,
        tolerance_angular: f64,
    ) {
        self.tol3d = tolerance3d;
        self.bound_tol = bound_tolerance;
        self.tol2d = tolerance2d;
        self.tol_angular = tolerance_angular;
    }

    /// OCCT SetForceApproxC1 (L151-154).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_force_approx_c1 = force_approx_c1;
    }

    /// OCCT ExchangeUV (L158-161).
    pub fn exchange_uv(&self) -> bool {
        self.my_exch_uv
    }

    /// OCCT UReversed (L165-168).
    pub fn u_reversed(&self) -> bool {
        self.is_u_reversed
    }

    /// OCCT VReversed (L172-175).
    pub fn v_reversed(&self) -> bool {
        self.is_v_reversed
    }

    /// OCCT Build (L179-232).
    pub fn build(
        &mut self,
        section: Rc<RefCell<dyn SectionLaw>>,
        methode: GeomFillApproxStyle,
        continuity: GeomAbsShape,
        degmax: i32,
        segmax: i32,
    ) {
        // Inits.
        self.done = false;
        self.my_exch_uv = false;
        self.is_u_reversed = false;
        self.is_v_reversed = false;
        self.my_sec = Some(section);

        if (self.s_first == self.s_last) && (self.s_last == 30.081996) {
            let mut f = self.s_first;
            let mut l = self.s_last;
            self.my_sec
                .as_ref()
                .expect("null mySec")
                .borrow()
                .get_domain(&mut f, &mut l);
            self.s_first = f;
            self.s_last = l;
        }

        let mut is_k_part = false;
        let mut is_product = false;

        // Traitement des KPart.
        if self.my_k_part {
            is_k_part = self.build_k_part();
        }

        if !is_k_part {
            self.my_exch_uv = false;
            self.is_u_reversed = false;
            self.is_v_reversed = false;
        }

        // Traitement des produits Formelles.
        if !is_k_part && methode == GeomFillApproxStyle::GeomFill_Location {
            let bs = self
                .my_sec
                .as_ref()
                .expect("null mySec")
                .borrow()
                .bspline_surface()
                .is_some();
            if bs {
                // Approx de la loi
                //    isProduct = BuildProduct(Continuity, Degmax, Segmax);
            }
        }

        if is_k_part || is_product {
            // Approx du 2d.
            self.done = self.build_2d(continuity, degmax, segmax);
        } else {
            // Approx globale.
            self.done = self.build_all(continuity, degmax, segmax);
        }
    }

    /// OCCT Build2d (L236-244).
    fn build_2d(&self, _continuity: GeomAbsShape, _degmax: i32, _segmax: i32) -> bool {
        let ok;
        if self.my_loc.borrow().nb_2d_curves() == 0 {
            ok = true;
        } else {
            ok = false;
        }
        ok
    }

    /// OCCT BuildAll (L248-388).
    fn build_all(&mut self, continuity: GeomAbsShape, degmax: i32, segmax: i32) -> bool {
        let func = SweepFunction::new(
            self.my_sec.clone().expect("null mySec"),
            self.my_loc.clone(),
            self.first,
            self.s_first,
            (self.s_last - self.s_first) / (self.last - self.first),
        );
        let mut approx = ApproxSweepApproximation::new(&func);

        approx.perform(
            self.first,
            self.last,
            self.tol3d,
            self.bound_tol,
            self.tol2d,
            self.tol_angular,
            continuity,
            degmax,
            segmax,
        );

        let mut ok = false;
        if approx.is_done() {
            ok = true;

            // La surface.
            let (mut u_degree, mut v_degree) = (0i32, 0i32);
            let (mut nb_u_poles, mut nb_v_poles) = (0i32, 0i32);
            let (mut nb_u_knots, mut nb_v_knots) = (0i32, 0i32);
            approx.surf_shape(
                &mut u_degree,
                &mut v_degree,
                &mut nb_u_poles,
                &mut nb_v_poles,
                &mut nb_u_knots,
                &mut nb_v_knots,
            );

            let mut poles = vec![Vec::new(); nb_u_poles as usize];
            let mut weights = vec![Vec::new(); nb_u_poles as usize];
            let mut u_knots = vec![0.0; nb_u_knots as usize];
            let mut v_knots = vec![0.0; nb_v_knots as usize];
            let mut u_mults = vec![0i32; nb_u_knots as usize];
            let mut v_mults = vec![0i32; nb_v_knots as usize];

            approx.surface(&mut poles, &mut weights, &mut u_knots, &mut v_knots, &mut u_mults, &mut v_mults);

            // OCCT: new Geom_BSplineSurface(Poles, Weights, UKnots, VKnots,
            // UMults, VMults, UDegree, VDegree, mySec->IsUPeriodic()) — the
            // rcad BSplineSurface carries the periodic flag by the knot
            // layout.
            let mut flat_u = Vec::new();
            for (k, m) in u_knots.iter().zip(u_mults.iter()) {
                for _ in 0..*m {
                    flat_u.push(*k);
                }
            }
            let mut flat_v = Vec::new();
            for (k, m) in v_knots.iter().zip(v_mults.iter()) {
                for _ in 0..*m {
                    flat_v.push(*k);
                }
            }
            let is_u_periodic = self
                .my_sec
                .as_ref()
                .expect("null mySec")
                .borrow()
                .is_u_periodic();
            let bs = BSplineSurface {
                degree_u: u_degree as usize,
                degree_v: v_degree as usize,
                knots_u: flat_u,
                knots_v: flat_v,
                control_points: poles,
                weights,
            };
            self.my_surface = Some(Surface3::BSpline(bs.clone()));
            self.s_error = approx.max_error_on_surf();

            if self.my_force_approx_c1 && !bspline_surface_is_cn_v(&bs, 1) {
                let the_tol = 1.0e-4;
                let the_u_cont = GeomAbsShape::C1;
                let the_v_cont = GeomAbsShape::C1;
                let deg_u = 14;
                let deg_v = 14;
                let nmax = 16;
                let the_prec = 1;

                let convert_approx = GeomConvertApproxSurface::new(
                    &bs, the_tol, the_u_cont, the_v_cont, deg_u, deg_v, nmax, the_prec,
                );
                if convert_approx.has_result() {
                    let new_surface = convert_approx.surface();
                    self.my_surface = Some(Surface3::BSpline(new_surface.clone()));
                    self.my_curve2d = Some(Vec::new()); // (1..2) — filled below
                    self.c_error = vec![0.0; 4]; // (1..2, 1..2)

                    // OCCT: the two boundary iso 2d lines at the U knots.
                    let u_first_knot = flat_knots_to_mults(&new_surface.knots_u).0[0];
                    let v_last_knot =
                        flat_knots_to_mults(&new_surface.knots_v).0.last().copied().unwrap_or(0.0);
                    let d = DVec2::new(0.0, 1.0); // gp_Dir2d::D::Y
                    let p = DVec2::new(u_first_knot, 0.0);
                    let tc1 = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(Curve2d::Line(Line2d::new(p, d))),
                        t_min: 0.0,
                        t_max: v_last_knot,
                    });
                    self.my_curve2d.as_mut().expect("myCurve2d").push(Some(tc1));
                    self.c_error_set(1, 1, 0.0);
                    self.c_error_set(2, 1, 0.0);

                    let u_last_knot =
                        flat_knots_to_mults(&new_surface.knots_u).0.last().copied().unwrap_or(0.0);
                    let p2 = DVec2::new(u_last_knot, 0.0);
                    let tc2 = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(Curve2d::Line(Line2d::new(p2, d))),
                        t_min: 0.0,
                        t_max: v_last_knot,
                    });
                    self.my_curve2d.as_mut().expect("myCurve2d").push(Some(tc2));
                    let length = self.my_curve2d.as_ref().expect("myCurve2d").len();
                    self.c_error_set(1, length, 0.0);
                    self.c_error_set(2, length, 0.0);

                    self.s_error = the_tol;
                }
            } // if (!mySurface->IsCNv(1))

            // Les Courbes 2d.
            if self.my_curve2d.is_none() {
                let trace_number = self.my_loc.borrow().trace_number();
                self.my_curve2d = Some(vec![None; 2 + trace_number]);
                self.c_error = vec![0.0; 2 * (2 + trace_number)];
                let mut ifin = 1usize;
                let ideb;

                if self.my_loc.borrow().has_first_restriction() {
                    ideb = 1;
                } else {
                    ideb = 2;
                }
                ifin += self.my_loc.borrow().trace_number();
                if self.my_loc.borrow().has_last_restriction() {
                    ifin += 1;
                }

                let mut kk = 1i32;
                let mut ii = ideb;
                while ii <= ifin {
                    // OCCT: C = new Geom2d_BSplineCurve(Approx.Curve2dPoles(kk),
                    // Curves2dKnots, Curves2dMults, Curves2dDegree).
                    let knots = approx.curves2d_knots();
                    let mults = approx.curves2d_mults();
                    let mut flat = Vec::new();
                    for (kv, m) in knots.iter().zip(mults.iter()) {
                        for _ in 0..*m {
                            flat.push(*kv);
                        }
                    }
                    let poles = approx.curve2d_poles(kk);
                    let c = BSplineCurve2 {
                        degree: approx.curves2d_degree() as usize,
                        knots: flat,
                        control_points: poles.clone(),
                        weights: vec![1.0; poles.len()],
                    };
                    self.my_curve2d.as_mut().expect("myCurve2d")[ii - 1] =
                        Some(Curve2d::BSpline(c));
                    self.c_error_set(1, ii, approx.max2d_error(kk));
                    self.c_error_set(2, ii, approx.max2d_error(kk));
                    kk += 1;
                    ii += 1;
                }

                // Si les courbes de restriction, ne sont pas calcules, on
                // prend les iso Bords.
                if !self.my_loc.borrow().has_first_restriction() {
                    let d = DVec2::new(0.0, 1.0); // gp_Dir2d::D::Y
                    let p = DVec2::new(u_knots[0], 0.0);
                    let tc = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(Curve2d::Line(Line2d::new(p, d))),
                        t_min: self.first,
                        t_max: self.last,
                    });
                    self.my_curve2d.as_mut().expect("myCurve2d")[0] = Some(tc);
                    self.c_error_set(1, 1, 0.0);
                    self.c_error_set(2, 1, 0.0);
                }

                if !self.my_loc.borrow().has_last_restriction() {
                    let d = DVec2::new(0.0, 1.0); // gp_Dir2d::D::Y
                    let p = DVec2::new(u_knots[u_knots.len() - 1], 0.0);
                    let tc = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(Curve2d::Line(Line2d::new(p, d))),
                        t_min: self.first,
                        t_max: self.last,
                    });
                    let length = self.my_curve2d.as_ref().expect("myCurve2d").len();
                    self.my_curve2d.as_mut().expect("myCurve2d")[length - 1] = Some(tc);
                    self.c_error_set(1, length, 0.0);
                    self.c_error_set(2, length, 0.0);
                }
            } // if (myCurve2d.IsNull())
        }
        ok
    }

    /// OCCT CError->SetValue(Row, Col, Value) — the (1..2, 1..N) array
    /// (1-based access helper).
    fn c_error_set(&mut self, row: usize, col: usize, value: f64) {
        let nb_col = self.c_error.len() / 2;
        self.c_error[(row - 1) * nb_col + (col - 1)] = value;
    }

    /// OCCT CError->Value(Row, Col).
    fn c_error_value(&self, row: usize, col: usize) -> f64 {
        let nb_col = self.c_error.len() / 2;
        self.c_error[(row - 1) * nb_col + (col - 1)]
    }

    /// OCCT BuildProduct (L392-477) — GAP carrier: the body depends on
    /// GeomFill_LocFunction + AdvApprox_ApproxAFunction /
    /// AdvApprox_PrefAndRec (outside the translated batch); the OCCT call
    /// site (Build L218) is commented out — the function is unreachable
    /// there and keeps the OCCT failure path.
    fn build_product(&mut self, _continuity: GeomAbsShape, _degmax: i32, _segmax: i32) -> bool {
        panic!(
            "GAP: GeomFill_Sweep::BuildProduct (GeomFill_LocFunction + AdvApprox_ApproxAFunction / AdvApprox_PrefAndRec) is not translated — see file header"
        )
    }

    /// OCCT BuildKPart (L577-1165).
    fn build_k_part(&mut self) -> bool {
        let mut ok = false;
        let mut is_u_periodic = false;
        let mut is_v_periodic = false;
        let mut is_trsf = true;

        is_u_periodic = self
            .my_sec
            .as_ref()
            .expect("null mySec")
            .borrow()
            .is_u_periodic();
        let mut s: Option<Surface3> = None;
        let mut v = DVec3::ZERO;
        let mut m = GpMat::identity();
        let mut levier;
        let mut error = 0.0;
        let mut u_first = 0.0;
        let mut v_first = self.first;
        let mut u_last = 0.0;
        let mut v_last = self.last;
        let tol = self.tol3d.min(self.bound_tol);

        // (1) Trajectoire Rectilignes -------------------------
        if self.my_loc.borrow().is_translation(&mut error) {
            // Donne de la translation.
            let mut ds = DVec3::ZERO;
            let mut dp;
            self.my_loc.borrow().d0(1.0, &mut m, &mut ds);
            self.my_loc.borrow().d0(0.0, &mut m, &mut v);
            dp = ds - v;
            dp = dp.normalize_or_zero();
            // OCCT: gp_GTrsf Tf; Tf.SetVectorialPart(M);
            // Tf.SetTranslationPart(V.XYZ()); Tf2.SetValues(...) with the
            // Standard_ConstructionError catch.
            let tf2 = match gp_trsf_set_values_from_parts(&m, v) {
                Ok(t) => t,
                Err(_) => {
                    is_trsf = false;
                    DAffine3::IDENTITY
                }
            };
            if !is_trsf {
                return false;
            }

            // (1.1) Cas Extrusion.
            let mut is_constant = false;
            let mut conical_error = 0.0;
            if self
                .my_sec
                .as_ref()
                .expect("null mySec")
                .borrow()
                .is_constant(&mut error)
            {
                is_constant = true;
                let section = self
                    .my_sec
                    .as_ref()
                    .expect("null mySec")
                    .borrow()
                    .constant_section();
                // GeomAdaptor_Curve AC(Section).
                u_first = curve_first_parameter(&section);
                u_last = curve_last_parameter(&section);

                // (1.1.a) Cas Plan.
                if matches!(section, Curve3::Line(_)) && is_trsf {
                    // OCCT OCC5073: the parallel-spine guard.
                    if !is_sweep_parallel_spine(
                        &self.my_loc,
                        self.my_sec.as_ref().expect("null mySec"),
                        tol,
                    ) {
                        return false;
                    }
                    let (l_origin, l_dir) = match &section {
                        Curve3::Line(l) => (l.origin, l.direction),
                        _ => unreachable!(),
                    };
                    // L.Transform(Tf2).
                    let l_dir = tf2.transform_vector3(l_dir).normalize_or_zero();
                    let l_origin = tf2.transform_point3(l_origin);
                    ds = l_dir;
                    levier = ds.dot(dp).abs();
                    self.s_error = error + levier * (self.last - self.first).abs();
                    if self.s_error <= tol {
                        ok = true;
                        // gp_Ax2 AxisOfPlane(L.Location(), DS ^ DP, DS);
                        // S = new Geom_Plane(AxisOfPlane).
                        let axis_of_plane = Ax2::new(l_origin, ds, ds.cross(dp));
                        s = Some(Surface3::Plane(Plane {
                            origin: axis_of_plane.location,
                            normal: axis_of_plane.direction,
                            u_dir: axis_of_plane.x_direction,
                            v_dir: axis_of_plane.y_direction,
                        }));
                    } else {
                        self.s_error = 0.0;
                    }
                }

                // (1.1.b) Cas Cylindrique.
                if matches!(section, Curve3::Circle(_)) && is_trsf {
                    let tol_prod = 1.0e-6;

                    // C.Transform(Tf2).
                    let c = match transform_curve(&section, &tf2) {
                        Curve3::Circle(cc) => cc,
                        _ => unreachable!(),
                    };

                    ds = c.normal;
                    levier = ds.cross(dp).length() * c.radius;
                    self.s_error = levier * (self.last - self.first).abs();
                    if self.s_error <= tol_prod {
                        ok = true;
                        // gp_Ax3 axe(C.Location(), DP, C.XDirection());
                        // S = new Geom_CylindricalSurface(axe, C.Radius()).
                        let axe = Ax3::from_pnt_n_vx(c.center, dp, c.x_dir);
                        s = Some(Surface3::Cylinder(CylindricalSurface {
                            origin: axe.axis.location,
                            axis: axe.direction(),
                            radius: c.radius,
                            ref_dir: axe.x_direction,
                            y_dir: None,
                        }));
                        if gp_vec_is_opposite(c.normal, axe.direction(), 0.1) {
                            // L'orientation parametrique est inversee.
                            let l = 2.0 * PI - u_first;
                            let f = 2.0 * PI - u_last;
                            u_first = f;
                            u_last = l;
                            self.is_u_reversed = true;
                        }
                    } else {
                        self.s_error = 0.0;
                    }
                }

                // (1.1.c) C'est bien une extrusion.
                if !ok {
                    if is_trsf {
                        // Section->Transform(Tf2); S = new
                        // Geom_SurfaceOfLinearExtrusion(Section, DP).
                        let transformed = transform_curve(&section, &tf2);
                        s = Some(Surface3::LinearExtrusion(LinearExtrusionSurface {
                            profile: Box::new(transformed),
                            direction: dp,
                        }));
                        self.s_error = 0.0;
                        ok = true;
                    } else {
                        // extrusion sur BSpline
                    }
                }
            }

            // (1.2) Cas conique.
            let is_conical = {
                let sec = self.my_sec.as_ref().expect("null mySec");
                sec.borrow().is_conical_law(&mut conical_error)
            };
            error = conical_error;
            if !is_constant && is_conical {
                let sec = self.my_sec.as_ref().expect("null mySec");

                let mut section = sec.borrow().circl_section(self.s_last);
                section = transform_curve(&section, &tf2);
                // Section->Translate(Last * DP).
                let last_dp = DAffine3::from_translation(self.last * dp);
                section = transform_curve(&section, &last_dp);
                let c2 = match &section {
                    Curve3::Circle(c) => c.clone(),
                    _ => panic!("GeomAdaptor_Curve::Circle"),
                };
                let centre2 = c2.center;
                let p2 = section.point_at(0.0);
                let dsection = section.derivative_at(0.0);
                let r2 = c2.radius;

                let mut section = sec.borrow().circl_section(self.s_first);
                section = transform_curve(&section, &tf2);
                let first_dp = DAffine3::from_translation(self.first * dp);
                section = transform_curve(&section, &first_dp);
                let c1 = match &section {
                    Curve3::Circle(c) => c.clone(),
                    _ => panic!("GeomAdaptor_Curve::Circle"),
                };
                let centre1 = c1.center;
                let p1 = section.point_at(0.0);
                let r1 = c1.radius;

                let mut section = sec
                    .borrow()
                    .circl_section(self.s_first - self.first * (self.s_last - self.s_first) / (self.last - self.first));
                section = transform_curve(&section, &tf2);
                let c0 = match &section {
                    Curve3::Circle(c) => c.clone(),
                    _ => panic!("GeomAdaptor_Curve::Circle"),
                };
                let centre0 = c0.center;

                let mut angle;
                let mut n = p1 - centre1;
                if n.length() < 1.0e-9 {
                    let bis = p2 - centre2;
                    n = bis;
                }
                let lv = p2 - p1;
                let dir = centre2 - centre1;

                angle = gp_vec_angle(lv, dir);
                if angle > 0.01 && angle < PI / 2.0 - 0.01 {
                    if r2 < r1 {
                        angle = -angle;
                    }
                    self.s_error = error;
                    // gp_Ax3 Axis(Centre0, Dir, N); S = new
                    // Geom_ConicalSurface(Axis, Angle, C.Radius()).
                    let axis = Ax3::from_pnt_n_vx(centre0, dir.normalize_or_zero(), n.normalize_or_zero());
                    s = Some(Surface3::Cone(ConicalSurface {
                        apex: axis.axis.location,
                        axis: axis.direction(),
                        radius: c0.radius,
                        half_angle_rad: angle,
                        ref_dir: axis.x_direction,
                    }));
                    // Calcul du glissement parametrique.
                    v_first = self.first / angle.cos();
                    v_last = self.last / angle.cos();

                    // Bornes en U.
                    u_first = curve_first_parameter(&section);
                    u_last = curve_last_parameter(&section);
                    // S->VIso(VLast)->D1(0, pbis, diso).
                    let conical_surface = match &s {
                        Some(Surface3::Cone(c)) => c.clone(),
                        _ => unreachable!(),
                    };
                    let diso = cone_u_derivative_at_v(&conical_surface, v_last);
                    if diso.length() > 1.0e-9 && dsection.length() > 1.0e-9 {
                        self.is_u_reversed = gp_vec_is_opposite(diso, dsection, 0.1);
                    }
                    if self.is_u_reversed {
                        // L'orientation parametrique est inversee.
                        let l = 2.0 * PI - u_first;
                        let f = 2.0 * PI - u_last;
                        u_first = f;
                        u_last = l;
                    }

                    // C'est un cone.
                    ok = true;
                }
            }
        }

        // (2) Trajectoire Circulaire.
        let mut rot_error = 0.0;
        if self.my_loc.borrow().is_rotation(&mut rot_error) {
            error = rot_error;
            let mut is_constant = false;
            {
                let sec = self.my_sec.as_ref().expect("null mySec");
                is_constant = sec.borrow().is_constant(&mut error);
            }
            if is_constant {
                // La trajectoire.
                let mut centre = DVec3::ZERO;
                is_v_periodic = (self.last - self.first - 2.0 * PI).abs() < 1.0e-15;
                let mut ds = DVec3::ZERO;
                let mut dp;
                let mut dn;
                self.my_loc.borrow().d0(0.1, &mut m, &mut ds);
                self.my_loc.borrow().d0(0.0, &mut m, &mut v);
                self.my_loc.borrow().rotation(&mut centre);

                dp = ds - v;
                ds = v - centre;
                let rot_radius = ds.length();
                if rot_radius > 1.0e-15 {
                    ds = ds.normalize_or_zero();
                } else {
                    return false; // Pas de KPart, rotation degeneree
                }
                dn = ds.cross(dp).normalize_or_zero();
                dp = dn.cross(ds).normalize_or_zero();

                // OCCT: gp_GTrsf Tf; Tf.SetVectorialPart(M);
                // Tf.SetTranslationPart(V.XYZ()); Tf2.SetValues(...) — the
                // exception handling is commented out unlike the translation
                // branch.  TODO: finding #9 - restoring the try/catch caused
                // blend regressions in the OCCT source itself.
                let tf2 =
                    gp_trsf_set_values_from_parts(&m, v).expect("Standard_ConstructionError: gp_Trsf::SetValues");

                // La section.
                let section = self
                    .my_sec
                    .as_ref()
                    .expect("null mySec")
                    .borrow()
                    .constant_section();
                u_first = curve_first_parameter(&section);
                u_last = curve_last_parameter(&section);

                // (2.1) Tore/Sphere ?
                if matches!(section, Curve3::Circle(_)) && is_trsf {
                    let is_good_side;
                    // C.Transform(Tf2).
                    let c = match transform_curve(&section, &tf2) {
                        Curve3::Circle(cc) => cc,
                        _ => unreachable!(),
                    };
                    // On calcul le centre eventuel.
                    let mut dc = c.center - centre;
                    centre += dc.dot(dn) * dn;
                    dc = c.center - centre;
                    let radius = dc.length(); // grand Rayon du tore
                    if radius > tol && dc.dot(ds) < 0.0 {
                        is_good_side = false;
                    } else {
                        is_good_side = true;
                    }
                    let mut dc = dc;
                    if radius < tol / 100.0 {
                        dc = ds; // Pour definir le tore
                    }

                    // On verifie d'abord que le plan de la section est // a
                    // l'axe de rotation.
                    let nc = c.normal.normalize_or_zero();
                    error = nc.dot(dn).abs();
                    // Puis on evalue l'erreur commise sur la section, en
                    // pivotant son plan ( pour contenir l'axe de rotation).
                    error += nc.dot(ds).abs();
                    error *= c.radius;
                    if error <= tol {
                        self.s_error = error;
                        error += radius;
                        if radius <= tol {
                            // (2.1.a) Sphere.
                            let mut f = u_first;
                            let mut l = u_last;
                            self.s_error = error;
                            // Centre.BaryCenter(1.0, C.Location(), 1.0).
                            centre = (centre + c.center) / 2.0;
                            // gp_Ax3 AxisOfSphere(Centre, DN, DS).
                            let axis_of_sphere =
                                Ax3::from_pnt_n_vx(centre, dn, ds);
                            let a_radius = c.radius;
                            let the_sphere = SphericalSurface {
                                center: axis_of_sphere.axis.location,
                                axis: axis_of_sphere.direction(),
                                radius: a_radius,
                                ref_dir: axis_of_sphere.x_direction,
                            };
                            // ElSLib::Parameters over the trimmed section
                            // end points (the OCCT transforms a trimmed copy
                            // then samples First/Last).
                            let fpar = curve_first_parameter(&section);
                            let lpar = curve_last_parameter(&section);
                            let the_section = transform_curve(
                                &Curve3::Trimmed(TrimmedCurve3::new(section.clone(), fpar, lpar)),
                                &tf2,
                            );
                            let first_point = the_section.point_at(curve_first_parameter(&the_section));
                            let last_point = the_section.point_at(curve_last_parameter(&the_section));
                            let sphere_frame = Ax3::from_pnt_n_vx(
                                the_sphere.center,
                                the_sphere.axis,
                                the_sphere.ref_dir,
                            );
                            let (ufirst_on_sec, vfirst_on_sec) =
                                el_slib_sphere_parameters(&sphere_frame, first_point);
                            let (ulast_on_sec, vlast_on_sec) =
                                el_slib_sphere_parameters(&sphere_frame, last_point);
                            if vfirst_on_sec < vlast_on_sec {
                                f = vfirst_on_sec;
                                l = vlast_on_sec;
                            } else {
                                // L'orientation parametrique est inversee.
                                f = vlast_on_sec;
                                l = vfirst_on_sec;
                                self.is_u_reversed = true;
                            }

                            if (l - f).abs() <= P_CONFUSION
                                || (ulast_on_sec - ufirst_on_sec).abs() > PI / 2.0
                            {
                                // l == f - "degenerated" surface;
                                // UlastOnSec - UfirstOnSec > M_PI_2 -
                                // "twisted" surface, it is impossible to
                                // represent with help of trimmed sphere.
                                self.is_u_reversed = false;
                                return ok;
                            }

                            if f >= -PI / 2.0 && l <= PI / 2.0 {
                                ok = true;
                                self.my_exch_uv = true;
                                u_first = f;
                                u_last = l;
                            } else {
                                // On restaure ce qu'il faut.
                                self.is_u_reversed = false;
                            }
                        } else if is_good_side {
                            // (2.1.b) Tore.
                            // gp_Ax3 AxisOfTore(Centre, DN, DC).
                            let axis_of_tore = Ax3::from_pnt_n_vx(centre, dn, dc);
                            let torus = ToroidalSurface {
                                center: axis_of_tore.axis.location,
                                axis: axis_of_tore.direction(),
                                ref_dir: axis_of_tore.x_direction,
                                major_radius: radius,
                                minor_radius: c.radius,
                            };
                            s = Some(Surface3::Torus(torus));

                            // OCCT: Iso = down_cast<Geom_Circle>(S->UIso(0.));
                            // axeiso = Iso->Circ().Position().
                            let torus_frame = Ax3::from_pnt_n_vx(
                                torus.center,
                                torus.axis,
                                torus.ref_dir,
                            );
                            let (iso_loc, iso_dir, iso_xdir) =
                                el_slib_torus_uiso_ax2(&torus_frame, torus.major_radius, torus.minor_radius, 0.0);

                            if gp_vec_is_opposite(c.normal, iso_dir, 0.1) {
                                // L'orientation parametrique est inversee.
                                let l = 2.0 * PI - u_first;
                                let f = 2.0 * PI - u_last;
                                u_first = f;
                                u_last = l;
                                self.is_u_reversed = true;
                            }
                            // On calcul le "glissement" parametrique.
                            let rot = super::sweep_section_generator::gp_vec_angle_with_ref(
                                c.x_dir,
                                iso_xdir,
                                iso_dir,
                            );
                            u_first -= rot;
                            u_last -= rot;

                            self.my_exch_uv = true;
                            // Attention l'arete de couture dans le cas
                            // periodique n'est peut etre pas a la bonne
                            // place...
                            if is_u_periodic && u_first.abs() > P_CONFUSION {
                                is_u_periodic = false; // Pour trimmer la surface...
                            }
                            ok = true;
                        }
                    } else {
                        self.s_error = 0.0;
                    }
                }
                // (2.2) Cone / Cylindre.
                if matches!(section, Curve3::Line(_)) && is_trsf {
                    // L.Transform(Tf2).
                    let (l_origin, l_dir) = match &section {
                        Curve3::Line(l) => (l.origin, l.direction),
                        _ => unreachable!(),
                    };
                    let l_origin = tf2.transform_point3(l_origin);
                    let l_dir = tf2.transform_vector3(l_dir);
                    let dl = l_dir;
                    levier = curve_first_parameter(&section).abs().max(curve_last_parameter(&section));
                    // si la line est ortogonale au cercle de rotation.
                    self.s_error = error + levier * dl.dot(dp).abs();
                    if self.s_error <= tol {
                        let dir_line = (centre, dn); // gp_Lin Dir(Centre, DN)
                        let mut aux = dl.dot(dn);
                        let reverse = aux < 0.0; // On choisit ici le sens de parametrisation

                        // Calcul du centre du vecteur supportant la
                        // "XDirection".
                        let o1o2 = centre - l_origin;
                        let trans = dn * dn.dot(o1o2);
                        let centre_of_surf = centre + trans;
                        let ds = l_origin - centre_of_surf;

                        error = self.s_error;
                        error += dl.cross(dn).length() * levier;
                        if error <= tol {
                            // (2.2.a) Cylindre — si la line est orthogonale
                            // au plan de rotation.
                            self.s_error = error;
                            //
                            let mut axis = Ax3::from_pnt_n_vx(centre_of_surf, (dir_line.1).normalize_or_zero(), ds);
                            if ds.length_squared() > 1.0e-24 {
                                // OCCT: gp::Resolution() squared comparison —
                                // SquareMagnitude() > gp::Resolution().
                                axis.set_x_direction(ds);
                            }
                            // L.Distance(CentreOfSurf).
                            let distance = (l_origin - centre_of_surf)
                                .cross(l_dir)
                                .length();
                            s = Some(Surface3::Cylinder(CylindricalSurface {
                                origin: axis.axis.location,
                                axis: axis.direction(),
                                radius: distance,
                                ref_dir: axis.x_direction,
                                y_dir: None,
                            }));
                            ok = true;
                            self.my_exch_uv = true;
                        } else {
                            // On evalue l'angle du cone.
                            let mut angle = gp_vec_angle(dir_line.1, l_dir);
                            if angle > PI / 2.0 {
                                angle = PI - angle;
                            }
                            if reverse {
                                angle = -angle;
                            }
                            aux = ds.dot(dl);
                            if aux < 0.0 {
                                angle = -angle;
                            }
                            if ((angle.abs()) - PI / 2.0).abs() > 0.01 {
                                // (2.2.b) Cone — si les 2 droites ne sont
                                // pas orthogonales.
                                let radius = (l_origin - centre_of_surf).length();
                                let axis = Ax3::from_pnt_n_vx(
                                    centre_of_surf,
                                    (dir_line.1).normalize_or_zero(),
                                    ds,
                                );
                                s = Some(Surface3::Cone(ConicalSurface {
                                    apex: axis.axis.location,
                                    axis: axis.direction(),
                                    radius,
                                    half_angle_rad: angle,
                                    ref_dir: axis.x_direction,
                                }));
                                self.my_exch_uv = true;
                                ok = true;
                            } else {
                                // On n'as pas conclue, on remet l'erreur a 0.
                                self.s_error = 0.0;
                            }
                        }
                        if ok && reverse {
                            // On reverse le parametre.
                            // OCCT: uf = CL->ReversedParameter(ULast) with
                            // CL = new Geom_Line(L) — the line reverses the
                            // parameter about its origin.
                            let uf = reversed_line_parameter(l_origin, l_dir, u_last);
                            let ul = reversed_line_parameter(l_origin, l_dir, u_first);
                            u_first = ul;
                            u_last = uf;
                        }
                    } else {
                        self.s_error = 0.0;
                    }
                }

                // (2.3) Revolution.
                if !ok {
                    if is_trsf {
                        let transformed = transform_curve(&section, &tf2);
                        // gp_Ax1 Axis(Centre, DN); S = new
                        // Geom_SurfaceOfRevolution(Section, Axis).
                        s = Some(Surface3::Revolution(RevolutionSurface {
                            profile: Box::new(transformed),
                            axis_origin: centre,
                            axis_dir: dn.normalize_or_zero(),
                        }));
                        self.my_exch_uv = true;
                        self.s_error = 0.0;
                        ok = true;
                    }
                }
            }
        }

        if ok {
            // On trimme la surface.
            if self.my_exch_uv {
                std::mem::swap(&mut is_u_periodic, &mut is_v_periodic);
                std::mem::swap(&mut u_first, &mut v_first);
                std::mem::swap(&mut u_last, &mut v_last);
            }

            if !is_u_periodic && !is_v_periodic {
                let basis = s.expect("null S");
                self.my_surface = Some(Surface3::Trimmed(TrimmedSurface {
                    basis: Box::new(basis),
                    trim: [u_first, u_last, v_first, v_last],
                }));
            } else if is_u_periodic {
                if is_v_periodic {
                    self.my_surface = s;
                } else {
                    // OCCT: RectangularTrimmedSurface(S, VFirst, VLast, False)
                    // — the V-bounds only (the U range keeps its natural
                    // domain).
                    let basis = s.expect("null S");
                    let d = basis.default_domain();
                    self.my_surface = Some(Surface3::Trimmed(TrimmedSurface {
                        basis: Box::new(basis),
                        trim: [d[0], d[1], v_first, v_last],
                    }));
                }
            } else {
                // OCCT: RectangularTrimmedSurface(S, UFirst, ULast, true) —
                // the U-bounds only.
                let basis = s.expect("null S");
                let d = basis.default_domain();
                self.my_surface = Some(Surface3::Trimmed(TrimmedSurface {
                    basis: Box::new(basis),
                    trim: [u_first, u_last, d[2], d[3]],
                }));
            }
        }

        ok
    }

    /// OCCT IsDone (L1169-1172).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT ErrorOnSurface (L1176-1179).
    pub fn error_on_surface(&self) -> f64 {
        self.s_error
    }

    /// OCCT ErrorOnRestriction (L1183-1197).
    pub fn error_on_restriction(&self, is_first: bool, u_error: &mut f64, v_error: &mut f64) {
        let ind;
        if is_first {
            ind = 1;
        } else {
            ind = self.my_curve2d.as_ref().map(|c| c.len()).unwrap_or(0);
        }

        *u_error = self.c_error_value(1, ind);
        *v_error = self.c_error_value(2, ind);
    }

    /// OCCT ErrorOnTrace (L1201-1211).
    pub fn error_on_trace(&self, index_of_trace: usize, u_error: &mut f64, v_error: &mut f64) {
        let ind = index_of_trace + 1;
        if index_of_trace > self.my_loc.borrow().trace_number() {
            panic!(" GeomFill_Sweep::ErrorOnTrace");
        }

        *u_error = self.c_error_value(1, ind);
        *v_error = self.c_error_value(2, ind);
    }

    /// OCCT Surface (L1215-1218).
    pub fn surface(&self) -> Option<&Surface3> {
        self.my_surface.as_ref()
    }

    /// OCCT Restriction (L1222-1229).
    pub fn restriction(&self, is_first: bool) -> Option<Curve2d> {
        let curves = self.my_curve2d.as_ref()?;
        if is_first {
            return curves[0].clone();
        }
        curves[curves.len() - 1].clone()
    }

    /// OCCT NumberOfTrace (L1233-1236).
    pub fn number_of_trace(&self) -> usize {
        self.my_loc.borrow().trace_number()
    }

    /// OCCT Trace (L1240-1248).
    pub fn trace(&self, index_of_trace: usize) -> Option<Curve2d> {
        let ind = index_of_trace + 1;
        if index_of_trace > self.my_loc.borrow().trace_number() {
            panic!(" GeomFill_Sweep::Trace");
        }
        self.my_curve2d
            .as_ref()
            .and_then(|c| c[ind - 1].clone())
    }
}

/// OCCT section_type scratch consumer (the GetType value is carried by the
/// Curve3 match arms; this keeps the unused-write literal form).
fn curve_first_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.first,
        other => other.default_domain()[0],
    }
}

fn curve_last_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.last,
        other => other.default_domain()[1],
    }
}

/// OCCT Geom_Line::ReversedParameter(U) — the parameter of the mirrored
/// point (2 * origin-projection - U): the Geom_Line form is
/// -U + 2 * (P - Location).Dot(Direction) evaluated as the reversal about
/// the line origin, i.e. ReversedParameter(U) = A - U with A = the origin
/// projection of U...  The literal Geom_Line definition is
/// ReversedParameter(U) = -U (Geom_Line.hxx) with the line evaluated from
/// its own origin — the reversal maps U -> -U about Location.
fn reversed_line_parameter(origin: DVec3, direction: DVec3, u: f64) -> f64 {
    // Geom_Line::ReversedParameter(U) = -U — but the reversed line is the
    // same geometry reparametrised about the origin; the consumed value is
    // the parametric mirror: 2 * 0 - U (the line origin is the reference).
    let _ = (origin, direction);
    -u
}

/// OCCT static IsSweepParallelSpine (L485-571) — the OCC5073 guard: the
/// constant Line section swept along a translation stays parallel to the
/// spine.
fn is_sweep_parallel_spine(
    the_loc: &Rc<RefCell<dyn LocationLaw>>,
    the_sec: &Rc<RefCell<dyn SectionLaw>>,
    the_tol: f64,
) -> bool {
    // Get the first and last transformations of the location.
    let (mut a_first, mut a_last) = (0.0, 0.0);
    let mut v_begin = DVec3::ZERO;
    let mut v_end = DVec3::ZERO;
    let mut m = GpMat::identity();

    the_loc.borrow().get_domain(&mut a_first, &mut a_last);

    // Get the first transformation.
    the_loc.borrow().d0(a_first, &mut m, &mut v_begin);
    // OCCT: GTfBegin.SetVectorialPart(M); SetTranslationPart(VBegin);
    // TfBegin.SetValues(...) — a raise propagates (no catch).
    let tf_begin = gp_trsf_set_values_from_parts(&m, v_begin)
        .expect("Standard_ConstructionError: gp_Trsf::SetValues");

    // Get the last transformation.
    the_loc.borrow().d0(a_last, &mut m, &mut v_end);
    let tf_end = gp_trsf_set_values_from_parts(&m, v_end)
        .expect("Standard_ConstructionError: gp_Trsf::SetValues");

    // The rcad bspline_surface() borrows the guard; the surface is cloned
    // out (the OCCT handle copies the shared pointer).
    let a_surf = the_sec.borrow().bspline_surface().cloned();
    let (umin, _umax, _vmin, vmax) = a_surf
        .as_ref()
        .map(|s| {
            let d = s.default_domain();
            (d[0], d[1], d[2], d[3])
        })
        .expect("null BSplineSurface");

    // Get and transform the first section.
    let first_section = the_sec.borrow().constant_section();
    let u_first = curve_first_parameter(&first_section);
    let (l_origin, l_dir) = match &first_section {
        Curve3::Line(l) => (l.origin, l.direction),
        _ => panic!("GeomAdaptor_Curve::Line"),
    };
    // L.Transform(TfBegin).
    let l_origin = tf_begin.transform_point3(l_origin);
    let l_dir = tf_begin.transform_vector3(l_dir);

    // Get and transform the last section — the VIso(Vmax) first point is
    // S(Umin, Vmax) (architecture note: only that point is consumed).
    let a_pnt_last_sec = a_surf.as_ref().expect("null surface").point_at(umin, vmax);
    let a_pnt_last_sec = tf_end.transform_point3(a_pnt_last_sec);

    // ElCLib::Value(UFirst, L).
    let a_pnt_first_sec = l_origin + u_first * l_dir;
    let a_vec_sec = a_pnt_last_sec - a_pnt_first_sec;
    let a_vec_spine = v_end - v_begin;

    gp_vec_is_parallel(a_vec_sec, a_vec_spine, the_tol)
}
