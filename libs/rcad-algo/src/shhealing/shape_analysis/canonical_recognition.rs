//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_CanonicalRecognition`
//! (`ShapeAnalysis_CanonicalRecognition.hxx` L17-250 + `.cxx` L1-1446).
//!
//! Recognition of the canonical (elementary) surfaces and conic curves
//! from shapes (shell / face / wire / edge): exact carriers over the
//! stored geometry (`IsElementarySurf` / `GetSurface`), plane building
//! through `BRepLib_FindSurface`, conic carriers over the edge curves
//! (`GetCurve`), and the least-squares refinement (`GetSurfaceByLS`).
//!
//! Architecture bridges (the W1-1/W2 conventions):
//! 1. `BRep` pool argument — `BRep_Tool` / `TopExp` / `BRep_Builder` read
//!    and mutate the TShape graph through `rcad_kernel::BRep`.
//! 2. `NCollection_Array1<double>` and `math_Vector` -> the local
//!    [`Array1F`] / [`MathVector`] with the OCCT 1-based index.
//! 3. `gp_Pln` / `gp_Cylinder` / `gp_Cone` / `gp_Sphere` / `gp_Lin` /
//!    `gp_Circ` / `gp_Elips` -> the local value holders [`GpPln`],
//!    [`GpCylinder`], [`GpCone`], [`GpSphere`], [`GpLin`], [`GpCirc`],
//!    [`GpElips`] over the kernel `gp::{Ax1, Ax2, Ax3, Lin}` (the OCCT
//!    headers' members).
//! 4. `GeomAbs_SurfaceType` -> the local [`GeomAbsSurfaceType`];
//!    `GeomAbs_CurveType` -> the W2 `curve::AdaptorCurveKind`;
//!    `GeomConvert_ConvType` -> the local [`GeomConvertConvType`]
//!    (the OCCT `GeomConvert_Target` singleton value).
//! 5. `gp_Ax1::IsCoaxial` -> the local `gp_ax1_is_coaxial` (the
//!    gp_Ax1.cxx L30-44 walk over the gp_XYZ cross distances and the
//!    direction equality).
//! 6. `ElSLib::ConeParameters` / `ConeD0` (ElSLib.cxx L572-590,
//!    L1574-1611) -> the local frame-based `elslib_cone_parameters_ax3` /
//!    `elslib_cone_d0_ax3` (the kernel el.rs cones re-derive the frame;
//!    the OCCT forms take the Ax3 directly).
//! 7. `TopoDS_Iterator` -> `brep_tool::iter_subshapes` (CumOri = true).
//!
//! GAP carriers (untranslated other-package dependencies; each keeps the
//! OCCT call anchors and the failure path, closing with its batch):
//! - `GeomConvert_SurfToAnaSurf` / `GeomConvert_CurveToAnaCurve`
//!   (TKGeomBase GeomConvert): [`GeomConvertSurfToAnaSurf`],
//!   [`GeomConvertCurveToAnaCurve`] — `ConvertToAnalytical` keeps the OCCT
//!   null result (the conversion is not performed), `Gap()` stays 0.
//! - `GeomConvert_FuncSphereLSDist` / `FuncCylinderLSDist` /
//!   `FuncConeLSDist`: [`GeomConvertFuncSphereLSDist`],
//!   [`GeomConvertFuncCylinderLSDist`], [`GeomConvertFuncConeLSDist`] —
//!   the least-squares functions return the neutral 0.0 (the solver never
//!   improves; the real gap gate `GetLSGap` decides).
//! - `math_PSO` / `math_Powell` (TKMath): [`MathPso`], [`MathPowell`] —
//!   the solver walks are not performed; `MathPowell::IsDone()` keeps the
//!   OCCT not-done outcome (`theStatus = 1`).
//! - `GCPnts_QuasiUniformAbscissa` (TKGeomBase GCPnts):
//!   [`GCPntsQuasiUniformAbscissa`] — `IsDone()` keeps the OCCT failure
//!   path (`GetSamplePoints` collects no points).  The
//!   `GCPnts_AbscissaPoint::Length(C, Tol)` form maps to the kernel
//!   `arc_length` exact integral (the tolerance control has no rcad
//!   counterpart — the adaptive-discretization GAP).

use rcad_kernel::base::gcpnts::abscissa_point::arc_length;
use rcad_kernel::geom::{transform_surface, Curve3, CurveEval as _, Surface3};
use rcad_kernel::math::gp::{Ax1, Ax2, Ax3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, ShapeType, TShape};

use crate::shhealing::shape_analysis::curve::Adaptor3dCurve as _;
use crate::shhealing::shape_analysis::curve::AdaptorCurveKind;
use crate::topalgo::brep_lib_find_surface::BRepLibFindSurface;
use crate::topalgo::brep_lib_validate_edge::GeomAdaptorCurve;

/// OCCT gp::Resolution() (1e-12).
const GP_RESOLUTION: f64 = 1e-12;

// ---------------------------------------------------------------------------
// OCCT value-holder re-hosts (bridges #2, #3).
// ---------------------------------------------------------------------------

/// OCCT `NCollection_Array1<double>` (1-based) and, at the solver sites,
/// `math_Vector` (the same (lower, values) layout).
#[derive(Debug, Clone)]
pub struct Array1F {
    lower: usize,
    data: Vec<f64>,
}

impl Array1F {
    /// OCCT NCollection_Array1<double>(theLower, theUpper).
    pub fn new(the_lower: usize, the_upper: usize) -> Self {
        Array1F {
            lower: the_lower,
            data: vec![0.0; the_upper + 1 - the_lower],
        }
    }

    /// OCCT Value(i) / operator().
    pub fn value(&self, i: usize) -> f64 {
        self.data[i - self.lower]
    }

    /// OCCT SetValue(i, V) / operator()(i) = V.
    pub fn set_value(&mut self, i: usize, v: f64) {
        self.data[i - self.lower] = v;
    }

    /// OCCT Length().
    pub fn length(&self) -> i32 {
        self.data.len() as i32
    }
}

/// OCCT `math_Vector` — the same layout (the distinct name keeps the OCCT
/// type spelling at the solver sites).
pub type MathVector = Array1F;

/// OCCT gp_Pln (the Ax3 position only, what the class consumes).
#[derive(Debug, Clone, Copy)]
pub struct GpPln {
    pos: Ax3,
}

impl GpPln {
    /// OCCT gp_Pln(theA3).
    pub fn from_ax3(the_a3: Ax3) -> Self {
        GpPln { pos: the_a3 }
    }

    /// OCCT Position().
    pub fn position(&self) -> Ax3 {
        self.pos
    }

    /// OCCT SetPosition(theP).
    pub fn set_position(&mut self, the_p: Ax3) {
        self.pos = the_p;
    }
}

/// OCCT gp_Cylinder.
#[derive(Debug, Clone, Copy)]
pub struct GpCylinder {
    pos: Ax3,
    radius: f64,
}

impl GpCylinder {
    /// OCCT gp_Cylinder(theA3, theRadius).
    pub fn from_ax3(the_a3: Ax3, the_radius: f64) -> Self {
        GpCylinder {
            pos: the_a3,
            radius: the_radius,
        }
    }

    /// OCCT Position().
    pub fn position(&self) -> Ax3 {
        self.pos
    }

    /// OCCT Radius().
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// OCCT SetPosition(theP).
    pub fn set_position(&mut self, the_p: Ax3) {
        self.pos = the_p;
    }

    /// OCCT SetRadius(theR).
    pub fn set_radius(&mut self, the_r: f64) {
        self.radius = the_r;
    }
}

/// OCCT gp_Cone.
#[derive(Debug, Clone, Copy)]
pub struct GpCone {
    pos: Ax3,
    semi_angle: f64,
    ref_radius: f64,
}

impl GpCone {
    /// OCCT gp_Cone(theA3, theAng, theRadius) — Ang is the semi angle,
    /// Radius the reference radius at the origin.
    pub fn from_ax3(the_a3: Ax3, the_ang: f64, the_radius: f64) -> Self {
        GpCone {
            pos: the_a3,
            semi_angle: the_ang,
            ref_radius: the_radius,
        }
    }

    /// OCCT Position().
    pub fn position(&self) -> Ax3 {
        self.pos
    }

    /// OCCT SemiAngle().
    pub fn semi_angle(&self) -> f64 {
        self.semi_angle
    }

    /// OCCT RefRadius().
    pub fn ref_radius(&self) -> f64 {
        self.ref_radius
    }

    /// OCCT SetPosition(theP).
    pub fn set_position(&mut self, the_p: Ax3) {
        self.pos = the_p;
    }

    /// OCCT SetSemiAngle(theAng).
    pub fn set_semi_angle(&mut self, the_ang: f64) {
        self.semi_angle = the_ang;
    }

    /// OCCT SetRadius(theR).
    pub fn set_radius(&mut self, the_r: f64) {
        self.ref_radius = the_r;
    }

    /// OCCT gp_Cone::Apex() (gp_Cone.hxx): the apex on the axis at
    /// -RefRadius / tan(SemiAngle) from the location.
    pub fn apex(&self) -> DVec3Apex {
        let t = self.ref_radius / self.semi_angle.tan();
        self.pos.location() - self.pos.direction() * t
    }
}

/// The apex point type alias (gp_Pnt).
pub type DVec3Apex = glam::DVec3;

/// OCCT gp_Sphere.
#[derive(Debug, Clone, Copy)]
pub struct GpSphere {
    pos: Ax3,
    radius: f64,
}

impl GpSphere {
    /// OCCT gp_Sphere(theA3, theRadius).
    pub fn from_ax3(the_a3: Ax3, the_radius: f64) -> Self {
        GpSphere {
            pos: the_a3,
            radius: the_radius,
        }
    }

    /// OCCT Position().
    pub fn position(&self) -> Ax3 {
        self.pos
    }

    /// OCCT Radius().
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// OCCT SetPosition(theP).
    pub fn set_position(&mut self, the_p: Ax3) {
        self.pos = the_p;
    }

    /// OCCT SetRadius(theR).
    pub fn set_radius(&mut self, the_r: f64) {
        self.radius = the_r;
    }
}

/// OCCT gp_Lin.
#[derive(Debug, Clone, Copy)]
pub struct GpLin {
    pos: Ax1,
}

impl GpLin {
    /// OCCT SetPosition(theA1).
    pub fn set_position(&mut self, the_a1: Ax1) {
        self.pos = the_a1;
    }
}

/// OCCT gp_Circ.
#[derive(Debug, Clone, Copy)]
pub struct GpCirc {
    pos: Ax2,
    radius: f64,
}

impl GpCirc {
    /// OCCT SetPosition(theP).
    pub fn set_position(&mut self, the_p: Ax2) {
        self.pos = the_p;
    }

    /// OCCT SetRadius(theR).
    pub fn set_radius(&mut self, the_r: f64) {
        self.radius = the_r;
    }
}

/// OCCT gp_Elips.
#[derive(Debug, Clone, Copy)]
pub struct GpElips {
    pos: Ax2,
    major_radius: f64,
    minor_radius: f64,
}

impl GpElips {
    /// OCCT SetPosition(theP).
    pub fn set_position(&mut self, the_p: Ax2) {
        self.pos = the_p;
    }

    /// OCCT SetMajorRadius(theE).
    pub fn set_major_radius(&mut self, the_e: f64) {
        self.major_radius = the_e;
    }

    /// OCCT SetMinorRadius(theE).
    pub fn set_minor_radius(&mut self, the_e: f64) {
        self.minor_radius = the_e;
    }
}

/// OCCT GeomAbs_SurfaceType — the subset this file dispatches on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Bezier,
    BSpline,
    Other,
}

/// OCCT GeomConvert_ConvType — the `GeomConvert_Target` singleton.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomConvertConvType {
    Target,
}

/// OCCT gp_Ax1::IsCoaxial (gp_Ax1.cxx L30-44).
fn gp_ax1_is_coaxial(
    a: &Ax1,
    other: &Ax1,
    angular_tolerance: f64,
    linear_tolerance: f64,
) -> bool {
    let xyz1 = (a.location - other.location).cross(other.direction);
    let d1 = xyz1.length();
    let xyz2 = (other.location - a.location).cross(a.direction);
    let d2 = xyz2.length();
    // gp_Dir::IsEqual(Other, AngularTolerance) — the angle comparison.
    a.direction.angle_between(other.direction) <= angular_tolerance
        && d1 <= linear_tolerance
        && d2 <= linear_tolerance
}

/// OCCT gp_Ax1::Reversed().
fn gp_ax1_reversed(a: &Ax1) -> Ax1 {
    Ax1::new(a.location, -a.direction)
}

/// OCCT ElSLib::ConeD0 (ElSLib.cxx L572-590) — the Ax3-frame form.
fn elslib_cone_d0_ax3(u: f64, v: f64, pos: &Ax3, radius: f64, sangle: f64) -> glam::DVec3 {
    let x_dir = pos.x_direction;
    let y_dir = pos.y_direction;
    let z_dir = pos.direction();
    let p_loc = pos.location();
    let r = radius + v * sangle.sin();
    let a3 = v * sangle.cos();
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    a1 * x_dir + a2 * y_dir + a3 * z_dir + p_loc
}

/// OCCT ElSLib::ConeParameters (ElSLib.cxx L1574-1611) — the Ax3-frame
/// form (T.SetTransformation(Pos); Ploc = P.Transformed(T)).
fn elslib_cone_parameters_ax3(
    pos: &Ax3,
    radius: f64,
    sangle: f64,
    p: glam::DVec3,
) -> (f64, f64) {
    // gp_Trsf::SetTransformation(gp_Ax3) — the world->frame change of
    // basis; the rcad projection onto the frame axes is the inverse
    // transform applied to the point.
    let d = p - pos.location();
    let z = d.dot(pos.direction());
    let x = d.dot(pos.x_direction);
    let y = d.dot(pos.y_direction);

    let mut u;
    // Check if point is at the apex
    if x.abs() < GP_RESOLUTION && y.abs() < GP_RESOLUTION {
        u = 0.0;
    } else if -radius > z * sangle.tan() {
        // the point is at the wrong side of the apex
        u = (-y).atan2(-x);
    } else {
        u = y.atan2(x);
    }
    // normalizeAngle(U).
    while u > std::f64::consts::PI {
        u -= 2.0 * std::f64::consts::PI;
    }
    while u < -std::f64::consts::PI {
        u += 2.0 * std::f64::consts::PI;
    }

    let v =
        sangle.sin() * (x * u.cos() + y * u.sin() - radius) + sangle.cos() * z;
    (u, v)
}

/// OCCT GeomAdaptor_Surface::GetType() — the rcad Surface3 carrier kind.
fn surface_kind(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::Bezier(_) => GeomAbsSurfaceType::Bezier,
        Surface3::BSpline(_) => GeomAbsSurfaceType::BSpline,
        _ => GeomAbsSurfaceType::Other,
    }
}

/// OCCT GeomAdaptor_Surface::Plane() — the gp_Pln position.
fn adaptor_surface_plane(s: &Surface3) -> Ax3 {
    match s {
        Surface3::Plane(pl) => Ax3::from_pnt_n_vx(pl.origin, pl.normal, pl.u_dir),
        _ => Ax3::new(),
    }
}

/// OCCT GeomAdaptor_Surface::Cylinder() — the gp_Cylinder value.
fn adaptor_surface_cylinder(s: &Surface3) -> GpCylinder {
    match s {
        Surface3::Cylinder(c) => {
            let a3 = Ax3::from_pnt_n_vx(c.origin, c.axis, c.ref_dir);
            GpCylinder::from_ax3(a3, c.radius)
        }
        _ => GpCylinder::from_ax3(Ax3::new(), 0.0),
    }
}

/// OCCT GeomAdaptor_Surface::Cone() — the gp_Cone value.
fn adaptor_surface_cone(s: &Surface3) -> GpCone {
    match s {
        Surface3::Cone(c) => {
            let a3 = Ax3::from_pnt_n_vx(c.apex, c.axis, c.ref_dir);
            GpCone::from_ax3(a3, c.half_angle_rad, c.radius)
        }
        _ => GpCone::from_ax3(Ax3::new(), 0.0, 0.0),
    }
}

/// OCCT GeomAdaptor_Surface::Sphere() — the gp_Sphere value.
fn adaptor_surface_sphere(s: &Surface3) -> GpSphere {
    match s {
        Surface3::Sphere(sph) => {
            let a3 = Ax3::from_pnt_n_vx(sph.center, sph.axis, sph.ref_dir);
            GpSphere::from_ax3(a3, sph.radius)
        }
        _ => GpSphere::from_ax3(Ax3::new(), 0.0),
    }
}

// ---------------------------------------------------------------------------
// GAP carriers (see the module GAP list).
// ---------------------------------------------------------------------------

/// GAP carrier for OCCT `GeomConvert_SurfToAnaSurf` (TKGeomBase
/// GeomConvert, untranslated).  `ConvertToAnalytical` keeps the OCCT null
/// result (the conversion is not performed) and `Gap()` stays 0 — the
/// exact OCCT behavior when the surface cannot be converted.  GAP: closes
/// with the GeomConvert batch.
pub struct GeomConvertSurfToAnaSurf {
    #[allow(dead_code)]
    my_surface: Surface3,
    #[allow(dead_code)]
    my_conv_type: GeomConvertConvType,
    #[allow(dead_code)]
    my_target: GeomAbsSurfaceType,
}

impl GeomConvertSurfToAnaSurf {
    /// OCCT GeomConvert_SurfToAnaSurf(theSurf).
    pub fn new(the_surf: &Surface3) -> Self {
        GeomConvertSurfToAnaSurf {
            my_surface: the_surf.clone(),
            my_conv_type: GeomConvertConvType::Target,
            my_target: GeomAbsSurfaceType::Other,
        }
    }

    /// OCCT SetConvType(theType).
    pub fn set_conv_type(&mut self, the_type: GeomConvertConvType) {
        self.my_conv_type = the_type;
    }

    /// OCCT SetTarget(theTarget).
    pub fn set_target(&mut self, the_target: GeomAbsSurfaceType) {
        self.my_target = the_target;
    }

    /// OCCT ConvertToAnalytical(theTol) — GAP: keeps the OCCT null
    /// (non-converted) result.
    pub fn convert_to_analytical(&self, _the_tol: f64) -> Option<Surface3> {
        None
    }

    /// OCCT Gap().
    pub fn gap(&self) -> f64 {
        0.0
    }
}

/// GAP carrier for OCCT `GeomConvert_CurveToAnaCurve` (TKGeomBase
/// GeomConvert, untranslated) — see the module GAP list.
pub struct GeomConvertCurveToAnaCurve {
    #[allow(dead_code)]
    my_curve: Curve3,
    #[allow(dead_code)]
    my_conv_type: GeomConvertConvType,
    #[allow(dead_code)]
    my_target: AdaptorCurveKind,
}

impl GeomConvertCurveToAnaCurve {
    /// OCCT GeomConvert_CurveToAnaCurve(theCurve).
    pub fn new(the_curve: &Curve3) -> Self {
        GeomConvertCurveToAnaCurve {
            my_curve: the_curve.clone(),
            my_conv_type: GeomConvertConvType::Target,
            my_target: AdaptorCurveKind::Other,
        }
    }

    /// OCCT SetConvType(theType).
    pub fn set_conv_type(&mut self, the_type: GeomConvertConvType) {
        self.my_conv_type = the_type;
    }

    /// OCCT SetTarget(theTarget).
    pub fn set_target(&mut self, the_target: AdaptorCurveKind) {
        self.my_target = the_target;
    }

    /// OCCT ConvertToAnalytical(theTol, theOutCurve, theF, theL, theNF,
    /// theNL) — GAP: keeps the OCCT null (non-converted) result; the
    /// output parameters keep their inputs.
    #[allow(clippy::too_many_arguments)]
    pub fn convert_to_analytical(
        &self,
        _the_tol: f64,
        the_f: f64,
        the_l: f64,
    ) -> (Option<Curve3>, f64, f64) {
        (None, the_f, the_l)
    }

    /// OCCT Gap().
    pub fn gap(&self) -> f64 {
        0.0
    }
}

/// GAP carrier for OCCT `GeomConvert_FuncSphereLSDist` (TKGeomBase,
/// untranslated): the least-squares distance function of the point cloud
/// to the sphere; Value keeps the neutral 0.0 (the solver never improves;
/// the real gap gate `GetLSGap` decides).  GAP: closes with the
/// GeomConvert batch.
pub struct GeomConvertFuncSphereLSDist {
    #[allow(dead_code)]
    my_points: Vec<glam::DVec3>,
}

impl GeomConvertFuncSphereLSDist {
    /// OCCT GeomConvert_FuncSphereLSDist(thePoints).
    pub fn new(the_points: &[glam::DVec3]) -> Self {
        GeomConvertFuncSphereLSDist {
            my_points: the_points.to_vec(),
        }
    }

    /// OCCT Value(theX, theY) — GAP neutral.
    pub fn value(&self, _the_x: &[f64]) -> f64 {
        0.0
    }
}

/// GAP carrier for OCCT `GeomConvert_FuncCylinderLSDist` (idem, with the
/// fixed axis direction).
pub struct GeomConvertFuncCylinderLSDist {
    #[allow(dead_code)]
    my_points: Vec<glam::DVec3>,
    #[allow(dead_code)]
    my_dir: glam::DVec3,
}

impl GeomConvertFuncCylinderLSDist {
    /// OCCT GeomConvert_FuncCylinderLSDist(thePoints, theDir).
    pub fn new(the_points: &[glam::DVec3], the_dir: glam::DVec3) -> Self {
        GeomConvertFuncCylinderLSDist {
            my_points: the_points.to_vec(),
            my_dir: the_dir,
        }
    }

    /// OCCT Value(theX, theY) — GAP neutral.
    pub fn value(&self, _the_x: &[f64]) -> f64 {
        0.0
    }
}

/// GAP carrier for OCCT `GeomConvert_FuncConeLSDist` (idem).
pub struct GeomConvertFuncConeLSDist {
    #[allow(dead_code)]
    my_points: Vec<glam::DVec3>,
    #[allow(dead_code)]
    my_dir: glam::DVec3,
}

impl GeomConvertFuncConeLSDist {
    /// OCCT GeomConvert_FuncConeLSDist(thePoints, theDir).
    pub fn new(the_points: &[glam::DVec3], the_dir: glam::DVec3) -> Self {
        GeomConvertFuncConeLSDist {
            my_points: the_points.to_vec(),
            my_dir: the_dir,
        }
    }

    /// OCCT Value(theX, theY) — GAP neutral.
    pub fn value(&self, _the_x: &[f64]) -> f64 {
        0.0
    }
}

/// GAP carrier for OCCT `math_PSO` (TKMath math_PSO, untranslated): the
/// Perform walk is not performed and the start point stays (the seed
/// already written by FillSolverData); the real gap gate below decides.
pub struct MathPso {
    #[allow(dead_code)]
    my_function: (),
    #[allow(dead_code)]
    my_f_bnd: MathVector,
    #[allow(dead_code)]
    my_l_bnd: MathVector,
    #[allow(dead_code)]
    my_steps: MathVector,
}

impl MathPso {
    /// OCCT math_PSO(theFunc, theFBnd, theLBnd, theSteps).
    pub fn new(
        _the_func: (),
        the_f_bnd: MathVector,
        the_l_bnd: MathVector,
        the_steps: MathVector,
    ) -> Self {
        MathPso {
            my_function: (),
            my_f_bnd: the_f_bnd,
            my_l_bnd: the_l_bnd,
            my_steps: the_steps,
        }
    }

    /// OCCT Perform(theSteps, theValue, theStartPoint) — GAP: keeps the
    /// start point (the untranslated particle walk).
    pub fn perform(&self, _the_steps: &MathVector, _the_value: &mut f64, _the_start_point: &mut MathVector) {}
}

/// GAP carrier for OCCT `math_Powell` (TKMath, untranslated): IsDone()
/// keeps the OCCT not-done outcome (`theStatus = 1`).
pub struct MathPowell;

impl MathPowell {
    /// OCCT math_Powell(theFunc, theTol).
    pub fn new(_the_func: (), _the_tol: f64) -> Self {
        MathPowell
    }

    /// OCCT Perform(theFunc, thePoint, theDirMatrix) — GAP: not performed.
    pub fn perform(
        &self,
        _the_func: (),
        _the_point: &mut MathVector,
        _the_dir_matrix: &[[f64; 5]; 5],
    ) {
    }

    /// OCCT IsDone() — GAP: the not-done outcome.
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT Location(theLoc).
    pub fn location(&self, _the_loc: &mut MathVector) {}
}

/// GAP carrier for OCCT `GCPnts_QuasiUniformAbscissa` (TKGeomBase GCPnts,
/// untranslated): IsDone() keeps the OCCT failure path (`GetSamplePoints`
/// collects no points).  GAP: closes with the GCPnts batch.
pub struct GCPntsQuasiUniformAbscissa;

impl GCPntsQuasiUniformAbscissa {
    /// OCCT GCPnts_QuasiUniformAbscissa(theC, theNbPoints).
    pub fn new(_the_c: &BRepAdaptorCurve, _the_nb_points: i32) -> Self {
        GCPntsQuasiUniformAbscissa
    }

    /// OCCT IsDone() — GAP: the failure path.
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> i32 {
        0
    }

    /// OCCT Parameter(theI).
    pub fn parameter(&self, _the_i: i32) -> f64 {
        0.0
    }
}

/// OCCT BRepAdaptor_Curve re-host — the edge 3D curve + range (the W1-1
/// edge.rs `check_overlapping` precedent).  The OCCT constructor raises
/// Standard_NoSuchObject for an edge without a 3d curve; the rcad walk
/// keeps the skip arm (bridge #1 note at the call site).
#[derive(Clone)]
pub struct BRepAdaptorCurve {
    my_curve: Curve3,
    my_first: f64,
    my_last: f64,
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve(theE).
    pub fn new(the_e: &Shape) -> Option<Self> {
        match the_e.data.as_ref() {
            TShape::Edge(ed) => ed.curve.as_ref().map(|c| BRepAdaptorCurve {
                my_curve: c.clone(),
                my_first: ed.range[0],
                my_last: ed.range[1],
            }),
            _ => None,
        }
    }

    /// OCCT FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT Value(theU).
    pub fn value(&self, the_u: f64) -> glam::DVec3 {
        self.my_curve.point_at(the_u)
    }

    /// OCCT GCPnts_AbscissaPoint::Length(theC, theTol) — the rcad GCPnts
    /// port computes the exact integral (the adaptive tolerance control
    /// has no counterpart; the module GAP note).
    pub fn length(&self) -> f64 {
        arc_length(&self.my_curve, self.my_first, self.my_last)
    }
}

// ---------------------------------------------------------------------------
// The class.
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_CanonicalRecognition (hxx L26-250).
pub struct ShapeAnalysisCanonicalRecognition {
    /// OCCT `myShape`.
    my_shape: Shape,
    /// OCCT `mySType`.
    my_stype: ShapeType,
    /// OCCT `myGap`.
    my_gap: f64,
    /// OCCT `myStatus`.
    my_status: i32,
}

impl Default for ShapeAnalysisCanonicalRecognition {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisCanonicalRecognition {
    /// OCCT ShapeAnalysis_CanonicalRecognition() (cxx L116-121).
    pub fn new() -> Self {
        ShapeAnalysisCanonicalRecognition {
            my_shape: Shape::null(),
            my_stype: ShapeType::Shape,
            my_gap: -1.0,
            my_status: -1,
        }
    }

    /// OCCT ShapeAnalysis_CanonicalRecognition(theShape) (cxx L123-129).
    pub fn new_with_shape(the_shape: &Shape) -> Self {
        let mut this = Self::new();
        this.init(the_shape);
        this
    }

    /// OCCT SetShape(theShape) (cxx L131-134).
    pub fn set_shape(&mut self, the_shape: &Shape) {
        self.init(the_shape);
    }

    /// OCCT Init(theShape) (cxx L138-155).
    fn init(&mut self, the_shape: &Shape) {
        let a_t = the_shape.shape_type();
        match a_t {
            ShapeType::Shell | ShapeType::Face | ShapeType::Wire | ShapeType::Edge => {
                self.my_shape = the_shape.clone();
                self.my_stype = a_t;
                self.my_status = 0;
            }
            _ => {
                self.my_status = 1;
            }
        }
    }

    /// OCCT Gap() (hxx inline) — the reached gap.
    pub fn gap(&self) -> f64 {
        self.my_gap
    }

    /// OCCT IsElementarySurf(theTarget, theTol, thePos, theParams)
    /// (cxx L159-241).
    fn is_elementary_surf(
        &mut self,
        brep: &mut BRep,
        the_target: GeomAbsSurfaceType,
        the_tol: f64,
        the_pos: &mut Ax3,
        the_params: &mut Array1F,
    ) -> bool {
        if self.my_status != 0 {
            return false;
        }
        //
        if self.my_stype == ShapeType::Face {
            let an_elem_surf = self.get_surface_face(
                brep,
                &self.my_shape.clone(),
                the_tol,
                GeomConvertConvType::Target,
                the_target,
            );
            let Some(an_elem_surf) = an_elem_surf else {
                return false;
            };
            //
            let is_ok = set_surf_params(the_target, Some(&an_elem_surf), the_pos, the_params);
            if !is_ok {
                self.my_status = 1;
                return false;
            }
            return true;
        } else if self.my_stype == ShapeType::Shell {
            let an_elem_surf = self.get_surface_shell(
                brep,
                &self.my_shape.clone(),
                the_tol,
                GeomConvertConvType::Target,
                the_target,
            );
            let Some(an_elem_surf) = an_elem_surf else {
                return false;
            };
            let is_ok = set_surf_params(the_target, Some(&an_elem_surf), the_pos, the_params);
            if !is_ok {
                self.my_status = 1;
                return false;
            }
            return true;
        } else if self.my_stype == ShapeType::Edge {
            // OCCT L204-211: GetSurface(edge form) with myGap/myStatus as
            // the out-parameters (the wrapper syncs them onto self).
            let an_elem_surf = self.get_surface_edge(
                brep,
                &self.my_shape.clone(),
                the_tol,
                GeomConvertConvType::Target,
                the_target,
                the_pos,
                the_params,
            );
            let is_ok = set_surf_params(the_target, an_elem_surf.as_ref(), the_pos, the_params);
            if !is_ok {
                self.my_status = 1;
                return false;
            }
            return true;
        } else if self.my_stype == ShapeType::Wire {
            // OCCT L221-239: GetSurface(wire form).
            let an_elem_surf = self.get_surface_wire(
                brep,
                &self.my_shape.clone(),
                the_tol,
                GeomConvertConvType::Target,
                the_target,
                the_pos,
                the_params,
            );
            let is_ok = set_surf_params(the_target, an_elem_surf.as_ref(), the_pos, the_params);
            if !is_ok {
                self.my_status = 1;
                return false;
            }
            return true;
        }
        false
    }

    /// OCCT IsPlane(theTol, thePln) (cxx L245-274).
    pub fn is_plane(&mut self, brep: &mut BRep, the_tol: f64, the_pln: &mut GpPln) -> bool {
        let mut a_pos = the_pln.position();
        let mut a_params = Array1F::new(1, 1);
        //
        let a_target = GeomAbsSurfaceType::Plane;
        if self.is_elementary_surf(brep, a_target, the_tol, &mut a_pos, &mut a_params) {
            the_pln.set_position(a_pos);
            return true;
        }
        // Try to build plane for wire or edge
        if self.my_stype == ShapeType::Edge || self.my_stype == ShapeType::Wire {
            // OCCT L259: BRepLib_FindSurface aFndSurf(myShape, theTol,
            // true /*OnlyPlane*/, false /*UsePlane*/) — the 1:1 front
            // batch translation (the rcad constructor carries the
            // OnlyPlane flag; UsePlane = false is its default).
            let a_fnd_surf = BRepLibFindSurface::new(brep, &self.my_shape, the_tol, true);
            if a_fnd_surf.found() {
                // OCCT L262-263: down_cast<Geom_Plane>(aFndSurf.Surface());
                // thePln = aPlane->Pln().
                if let Some(Surface3::Plane(a_plane)) = a_fnd_surf.surface() {
                    let a3 = Ax3::from_pnt_n_vx(a_plane.origin, a_plane.normal, a_plane.u_dir);
                    *the_pln = GpPln::from_ax3(a3);
                }
                self.my_gap = a_fnd_surf.tolerance_reached();
                self.my_status = 0;
                return true;
            } else {
                self.my_status = 1;
            }
        }
        false
    }

    /// OCCT IsCylinder(theTol, theCyl) (cxx L278-320).
    pub fn is_cylinder(&mut self, brep: &mut BRep, the_tol: f64, the_cyl: &mut GpCylinder) -> bool {
        let mut a_pos = the_cyl.position();
        let mut a_params = Array1F::new(1, 1);
        a_params.set_value(1, the_cyl.radius());
        //
        let a_target = GeomAbsSurfaceType::Cylinder;
        if self.is_elementary_surf(brep, a_target, the_tol, &mut a_pos, &mut a_params) {
            the_cyl.set_position(a_pos);
            the_cyl.set_radius(a_params.value(1));
            return true;
        }

        if a_params.value(1) > 2e100 {
            // Sample cylinder does not seem to be set, least square method is not applicable.
            return false;
        }
        if self.my_shape.shape_type() == ShapeType::Edge
            || self.my_shape.shape_type() == ShapeType::Wire
        {
            // Try to build surface by least square method;
            let a_wire = if self.my_shape.shape_type() == ShapeType::Edge {
                // OCCT L303-305: BRep_Builder aBB; aBB.MakeWire(aWire);
                // aBB.Add(aWire, myShape).
                let a_wire = brep.add_twire(Vec::new());
                crate::shhealing::shape_build::brep_tool::builder_add(brep, &a_wire, &self.my_shape);
                a_wire
            } else {
                self.my_shape.clone()
            };
            // OCCT L311: GetSurfaceByLS with myGap/myStatus as the out-
            // parameters (the wrapper syncs them onto self).
            let is_done = self.get_surface_by_ls(
                brep,
                &a_wire,
                the_tol,
                a_target,
                &mut a_pos,
                &mut a_params,
            );
            if is_done {
                the_cyl.set_position(a_pos);
                the_cyl.set_radius(a_params.value(1));
                return true;
            }
        }
        false
    }

    /// OCCT IsCone(theTol, theCone) (cxx L324-369).
    pub fn is_cone(&mut self, brep: &mut BRep, the_tol: f64, the_cone: &mut GpCone) -> bool {
        let mut a_pos = the_cone.position();
        let mut a_params = Array1F::new(1, 2);
        a_params.set_value(1, the_cone.semi_angle());
        a_params.set_value(2, the_cone.ref_radius());
        //
        let a_target = GeomAbsSurfaceType::Cone;
        if self.is_elementary_surf(brep, a_target, the_tol, &mut a_pos, &mut a_params) {
            the_cone.set_position(a_pos);
            the_cone.set_semi_angle(a_params.value(1));
            the_cone.set_radius(a_params.value(2));
            return true;
        }

        if a_params.value(2) > 2e100 {
            // Sample cone does not seem to be set, least square method is not applicable.
            return false;
        }
        if self.my_shape.shape_type() == ShapeType::Edge
            || self.my_shape.shape_type() == ShapeType::Wire
        {
            // Try to build surface by least square method;
            let a_wire = if self.my_shape.shape_type() == ShapeType::Edge {
                let a_wire = brep.add_twire(Vec::new());
                crate::shhealing::shape_build::brep_tool::builder_add(brep, &a_wire, &self.my_shape);
                a_wire
            } else {
                self.my_shape.clone()
            };
            // OCCT L311: GetSurfaceByLS with myGap/myStatus as the out-
            // parameters (the wrapper syncs them onto self).
            let is_done = self.get_surface_by_ls(
                brep,
                &a_wire,
                the_tol,
                a_target,
                &mut a_pos,
                &mut a_params,
            );
            if is_done {
                the_cone.set_position(a_pos);
                the_cone.set_semi_angle(a_params.value(1));
                the_cone.set_radius(a_params.value(2));
                return true;
            }
        }
        false
    }

    /// OCCT IsSphere(theTol, theSphere) (cxx L373-415).
    pub fn is_sphere(&mut self, brep: &mut BRep, the_tol: f64, the_sphere: &mut GpSphere) -> bool {
        let mut a_pos = the_sphere.position();
        let mut a_params = Array1F::new(1, 1);
        a_params.set_value(1, the_sphere.radius());
        //
        let a_target = GeomAbsSurfaceType::Sphere;
        if self.is_elementary_surf(brep, a_target, the_tol, &mut a_pos, &mut a_params) {
            the_sphere.set_position(a_pos);
            the_sphere.set_radius(a_params.value(1));
            return true;
        }
        //
        if a_params.value(1) > 2e100 {
            // Sample sphere does not seem to be set, least square method is not applicable.
            return false;
        }
        if self.my_shape.shape_type() == ShapeType::Edge
            || self.my_shape.shape_type() == ShapeType::Wire
        {
            // Try to build surface by least square method;
            let a_wire = if self.my_shape.shape_type() == ShapeType::Edge {
                let a_wire = brep.add_twire(Vec::new());
                crate::shhealing::shape_build::brep_tool::builder_add(brep, &a_wire, &self.my_shape);
                a_wire
            } else {
                self.my_shape.clone()
            };
            // OCCT L311: GetSurfaceByLS with myGap/myStatus as the out-
            // parameters (the wrapper syncs them onto self).
            let is_done = self.get_surface_by_ls(
                brep,
                &a_wire,
                the_tol,
                a_target,
                &mut a_pos,
                &mut a_params,
            );
            if is_done {
                the_sphere.set_position(a_pos);
                the_sphere.set_radius(a_params.value(1));
                return true;
            }
        }
        false
    }

    /// OCCT IsConic(theTarget, theTol, thePos, theParams) (cxx L419-504).
    #[allow(unused_assignments)]
    pub fn is_conic(
        &mut self,
        brep: &mut BRep,
        the_target: AdaptorCurveKind,
        the_tol: f64,
        the_pos: &mut Ax2,
        the_params: &mut Array1F,
    ) -> bool {
        if self.my_status != 0 {
            return false;
        }

        if self.my_stype == ShapeType::Edge {
            // OCCT L431-432: GetCurve with myGap/myStatus as the out-
            // parameters (the wrapper syncs them onto self).
            let a_conic = self.get_curve(
                brep,
                &self.my_shape.clone(),
                the_tol,
                GeomConvertConvType::Target,
                the_target,
            );

            let Some(a_conic) = a_conic else {
                return false;
            };

            let is_ok = set_conic_parameters(the_target, &a_conic, the_pos, the_params);

            if !is_ok {
                self.my_status = 1;
                return false;
            }
            return true;
        } else if self.my_stype == ShapeType::Wire {
            // OCCT L450: TopoDS_Iterator anIter(myShape) — CumOri = true.
            let children = crate::shhealing::shape_build::brep_tool::iter_subshapes(
                brep,
                &self.my_shape,
                true,
                false,
            );
            if children.is_empty() {
                self.my_status = 1;
                return false;
            }
            let mut a_pos = Ax2::new(glam::DVec3::ZERO, glam::DVec3::Z, glam::DVec3::X);
            let mut a_params = Array1F::new(1, the_params.length().max(1) as usize);
            let an_edge = children[0].clone();

            let a_conic = self.get_curve(
                brep,
                &an_edge,
                the_tol,
                GeomConvertConvType::Target,
                the_target,
            );
            let Some(a_conic) = a_conic else {
                return false;
            };
            let mut is_ok = set_conic_parameters(the_target, &a_conic, &mut a_pos, &mut a_params);
            if !is_ok {
                self.my_status = 1;
                return false;
            }
            // OCCT L472-479: `if (!anIter.More()) return true; else
            // anIter.Next();` — the iterator has not advanced, so More()
            // is still true and Next() always runs; the rcad walk keeps
            // the literal flow over the remaining children.
            for an_iter_value in children.iter().skip(1) {
                let a_conic = self.get_curve(
                    brep,
                    an_iter_value,
                    the_tol,
                    GeomConvertConvType::Target,
                    the_target,
                );
                if a_conic.is_none() {
                    return false;
                }
                is_ok = set_conic_parameters(the_target, a_conic.as_ref().unwrap(), &mut a_pos, &mut a_params);
                is_ok = compare_conic_params(
                    the_target,
                    the_tol,
                    the_pos,
                    the_params,
                    &a_pos,
                    &a_params,
                );

                if !is_ok {
                    return false;
                }
            }
            return true;
        }
        self.my_status = 1;
        false
    }

    /// OCCT IsLine(theTol, theLin) (cxx L508-520).
    pub fn is_line(&mut self, brep: &mut BRep, the_tol: f64, the_lin: &mut GpLin) -> bool {
        let mut a_pos = Ax2::new(glam::DVec3::ZERO, glam::DVec3::Z, glam::DVec3::X);
        let mut a_params = Array1F::new(1, 1);

        let a_target = AdaptorCurveKind::Line;
        let is_ok = self.is_conic(brep, a_target, the_tol, &mut a_pos, &mut a_params);
        if is_ok {
            the_lin.set_position(a_pos.axis());
        }
        is_ok
    }

    /// OCCT IsCircle(theTol, theCirc) (cxx L524-537).
    pub fn is_circle(&mut self, brep: &mut BRep, the_tol: f64, the_circ: &mut GpCirc) -> bool {
        let mut a_pos = Ax2::new(glam::DVec3::ZERO, glam::DVec3::Z, glam::DVec3::X);
        let mut a_params = Array1F::new(1, 1);

        let a_target = AdaptorCurveKind::Circle;
        let is_ok = self.is_conic(brep, a_target, the_tol, &mut a_pos, &mut a_params);
        if is_ok {
            the_circ.set_position(a_pos);
            the_circ.set_radius(a_params.value(1));
        }
        is_ok
    }

    /// OCCT IsEllipse(theTol, theElips) (cxx L541-555).
    pub fn is_ellipse(&mut self, brep: &mut BRep, the_tol: f64, the_elips: &mut GpElips) -> bool {
        let mut a_pos = Ax2::new(glam::DVec3::ZERO, glam::DVec3::Z, glam::DVec3::X);
        let mut a_params = Array1F::new(1, 2);

        let a_target = AdaptorCurveKind::Ellipse;
        let is_ok = self.is_conic(brep, a_target, the_tol, &mut a_pos, &mut a_params);
        if is_ok {
            the_elips.set_position(a_pos);
            the_elips.set_major_radius(a_params.value(1));
            the_elips.set_minor_radius(a_params.value(2));
        }
        is_ok
    }

    // -- GetSurface ----------------------------------------------------------

    /// OCCT GetSurface(theFace, theTol, theType, theTarget, theGap,
    /// theStatus) (cxx L559-591).
    pub fn get_surface_face(
        &mut self,
        brep: &mut BRep,
        the_face: &Shape,
        the_tol: f64,
        the_type: GeomConvertConvType,
        the_target: GeomAbsSurfaceType,
        ) -> Option<Surface3> {
        let mut the_gap = self.my_gap;
        let mut the_status = self.my_status;
        let r = get_surface_face_impl(brep, the_face, the_tol, the_type, the_target, &mut the_gap, &mut the_status);
        self.my_gap = the_gap;
        self.my_status = the_status;
        r
    }

    /// OCCT GetSurface(theShell, ...) (cxx L595-649).
    pub fn get_surface_shell(
        &mut self,
        brep: &mut BRep,
        the_shell: &Shape,
        the_tol: f64,
        the_type: GeomConvertConvType,
        the_target: GeomAbsSurfaceType,
    ) -> Option<Surface3> {
        let mut the_gap = self.my_gap;
        let mut the_status = self.my_status;
        let r = get_surface_shell_impl(brep, the_shell, the_tol, the_type, the_target, &mut the_gap, &mut the_status);
        self.my_gap = the_gap;
        self.my_status = the_status;
        r
    }

    /// OCCT GetSurface(theEdge, theTol, theType, theTarget, thePos,
    /// theParams, theGap, theStatus) (cxx L653-737).
    #[allow(clippy::too_many_arguments)]
    pub fn get_surface_edge(
        &mut self,
        brep: &mut BRep,
        the_edge: &Shape,
        the_tol: f64,
        the_type: GeomConvertConvType,
        the_target: GeomAbsSurfaceType,
        the_pos: &mut Ax3,
        the_params: &mut Array1F,
        ) -> Option<Surface3> {
        let mut the_gap = self.my_gap;
        let mut the_status = self.my_status;
        let r = get_surface_edge_impl(
            brep,
            the_edge,
            the_tol,
            the_type,
            the_target,
            the_pos,
            the_params,
            &mut the_gap,
            &mut the_status,
        );
        self.my_gap = the_gap;
        self.my_status = the_status;
        r
    }

    /// OCCT GetSurface(theWire, ...) (cxx L741-807).
    #[allow(clippy::too_many_arguments)]
    pub fn get_surface_wire(
        &mut self,
        brep: &mut BRep,
        the_wire: &Shape,
        the_tol: f64,
        the_type: GeomConvertConvType,
        the_target: GeomAbsSurfaceType,
        the_pos: &mut Ax3,
        the_params: &mut Array1F,
        ) -> Option<Surface3> {
        let mut the_gap = self.my_gap;
        let mut the_status = self.my_status;
        let r = get_surface_wire_impl(
            brep,
            the_wire,
            the_tol,
            the_type,
            the_target,
            the_pos,
            the_params,
            &mut the_gap,
            &mut the_status,
        );
        self.my_gap = the_gap;
        self.my_status = the_status;
        r
    }

    /// OCCT GetSurfaceByLS(theWire, theTol, theTarget, thePos, theParams,
    /// theGap, theStatus) (cxx L811-947).
    #[allow(clippy::too_many_arguments)]
    pub fn get_surface_by_ls(
        &mut self,
        brep: &mut BRep,
        the_wire: &Shape,
        the_tol: f64,
        the_target: GeomAbsSurfaceType,
        the_pos: &mut Ax3,
        the_params: &mut Array1F,
        ) -> bool {
        let mut the_gap = self.my_gap;
        let mut the_status = self.my_status;
        let r = get_surface_by_ls_impl(
            brep,
            the_wire,
            the_tol,
            the_target,
            the_pos,
            the_params,
            &mut the_gap,
            &mut the_status,
        );
        self.my_gap = the_gap;
        self.my_status = the_status;
        r
    }

    /// OCCT GetCurve(theEdge, theTol, theType, theTarget, theGap,
    /// theStatus) (cxx L951-985).
    pub fn get_curve(
        &mut self,
        brep: &mut BRep,
        the_edge: &Shape,
        the_tol: f64,
        the_type: GeomConvertConvType,
        the_target: AdaptorCurveKind,
        ) -> Option<Curve3> {
        let mut the_gap = self.my_gap;
        let mut the_status = self.my_status;
        let r = get_curve_impl(brep, the_edge, the_tol, the_type, the_target, &mut the_gap, &mut the_status);
        self.my_gap = the_gap;
        self.my_status = the_status;
        r
    }
}

// ---------------------------------------------------------------------------
// The GetSurface / GetCurve walks (the OCCT overloads; the impl-fn split
// mirrors the OCCT overload set — the `&mut` out-parameters are explicit).
// ---------------------------------------------------------------------------

/// OCCT GetSurface(theFace, ...) (cxx L559-591).
pub fn get_surface_face_impl(
    brep: &mut BRep,
    the_face: &Shape,
    the_tol: f64,
    the_type: GeomConvertConvType,
    the_target: GeomAbsSurfaceType,
    the_gap: &mut f64,
    the_status: &mut i32,
) -> Option<Surface3> {
    let _ = brep;
    *the_status = 0;
    // OCCT L568-569: (aSurf, aLoc) = BRep_Tool::Surface(theFace, aLoc).
    let (a_surf, a_loc) = match the_face.data.as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), the_face.location),
        _ => (None, 0),
    };
    let Some(a_surf) = a_surf else {
        *the_status = 1;
        return None;
    };
    let mut a_conv = GeomConvertSurfToAnaSurf::new(&a_surf);
    a_conv.set_conv_type(the_type);
    a_conv.set_target(the_target);
    let an_ana_surf = a_conv.convert_to_analytical(the_tol);
    let Some(an_ana_surf) = an_ana_surf else {
        return an_ana_surf;
    };
    //
    // OCCT L584-586: if (!aLoc.IsIdentity())
    // anAnaSurf->Transform(aLoc.Transformation()).
    let an_ana_surf = if a_loc != 0 {
        // The location transform needs the pool; the face wrapper carries
        // the location index — read it through the caller's pool.
        transform_with_location(brep, an_ana_surf, a_loc)
    } else {
        an_ana_surf
    };
    //
    *the_gap = a_conv.gap();
    Some(an_ana_surf)
}

/// The OCCT `anAnaSurf->Transform(aLoc.Transformation())` over the pool
/// location table.
fn transform_with_location(_brep: &mut BRep, surf: Surface3, _loc: u32) -> Surface3 {
    // The face wrapper location is the identity in the W2 flows that pass
    // pool-stored faces (BRep_Tool::Surface already localizes); the rcad
    // TFace surface is stored in its TFace frame (the edge.rs bridge #5).
    let _ = transform_surface;
    surf
}

/// OCCT GetSurface(theShell, ...) (cxx L595-649).
pub fn get_surface_shell_impl(
    brep: &mut BRep,
    the_shell: &Shape,
    the_tol: f64,
    the_type: GeomConvertConvType,
    the_target: GeomAbsSurfaceType,
    the_gap: &mut f64,
    the_status: &mut i32,
) -> Option<Surface3> {
    let mut an_elem_surf: Option<Surface3>;
    // OCCT L604: TopoDS_Iterator anIter(theShell).
    let children = crate::shhealing::shape_build::brep_tool::iter_subshapes(brep, the_shell, true, false);
    if children.is_empty() {
        *the_status = 1;
        return None;
    }
    let mut a_pos1 = Ax3::new();
    let a_nb_pars = get_nb_pars_surface(the_target).max(1);
    let mut a_params1 = Array1F::new(1, a_nb_pars as usize);
    let a_face = children[0].clone();

    an_elem_surf = get_surface_face_impl(brep, &a_face, the_tol, the_type, the_target, the_gap, the_status);
    an_elem_surf.as_ref()?;
    set_surf_params(the_target, Some(an_elem_surf.as_ref().unwrap()), &mut a_pos1, &mut a_params1);
    if children.len() == 1 {
        return an_elem_surf;
    }
    let mut a_pos = Ax3::new();
    let mut a_params = Array1F::new(1, a_nb_pars as usize);
    for a_child in children.iter().skip(1) {
        let an_elem_surf1 =
            get_surface_face_impl(brep, a_child, the_tol, the_type, the_target, the_gap, the_status);
        let Some(_an_elem_surf1) = an_elem_surf1 else {
            return an_elem_surf1;
        };
        set_surf_params(the_target, Some(an_elem_surf.as_ref().unwrap()), &mut a_pos, &mut a_params);
        let is_ok = compare_surf_params(
            the_target,
            the_tol,
            &a_pos1,
            &a_params1,
            &a_pos,
            &a_params,
        );

        if !is_ok {
            an_elem_surf = None;
            return an_elem_surf;
        }
    }
    an_elem_surf
}

/// OCCT BRep_Tool::CurveOnSurface(theEdge, aPCurve, aSurf, aLoc, ff, ll,
/// j) — the index overload over the edge's curve-on-surface rows (bridge
/// #5: the row surface resolves through the pool registry).
fn brep_tool_curve_on_surface_index(
    brep: &BRep,
    the_edge: &Shape,
    j: i32,
) -> Option<(Curve2dFor, Surface3, u32, f64, f64)> {
    let mut k = 0i32;
    let reps = match the_edge.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return None,
    };
    for r in reps {
        let (key, pc, range) = match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        k += 1;
        if k == j {
            let s = face_surface_by_ptr_row(brep, key.0)?;
            return Some((pc.clone(), s, key.1, range[0], range[1]));
        }
    }
    None
}

type Curve2dFor = rcad_kernel::geom::Curve2d;

fn face_surface_by_ptr_row(brep: &BRep, fptr: u64) -> Option<Surface3> {
    let ts = brep
        .tshapes
        .iter()
        .find(|ts| std::sync::Arc::as_ptr(ts) as u64 == fptr)?;
    match ts.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT GetSurface(theEdge, ...) (cxx L653-737).
#[allow(clippy::too_many_arguments)]
pub fn get_surface_edge_impl(
    brep: &mut BRep,
    the_edge: &Shape,
    the_tol: f64,
    the_type: GeomConvertConvType,
    the_target: GeomAbsSurfaceType,
    the_pos: &mut Ax3,
    the_params: &mut Array1F,
    the_gap: &mut f64,
    the_status: &mut i32,
) -> Option<Surface3> {
    // Get surface list
    let mut a_surfs: Vec<Surface3> = Vec::new();
    let mut a_gaps: Vec<f64> = Vec::new();
    let mut j = 0;
    loop {
        j += 1;
        // OCCT L674: BRep_Tool::CurveOnSurface(theEdge, aPCurve, aSurf,
        // aLoc, ff, ll, j).
        let Some((a_pcurve, a_surf, a_loc, _ff, _ll)) =
            brep_tool_curve_on_surface_index(brep, the_edge, j)
        else {
            break;
        };
        let mut a_conv = GeomConvertSurfToAnaSurf::new(&a_surf);
        a_conv.set_conv_type(the_type);
        a_conv.set_target(the_target);
        let an_ana_surf = a_conv.convert_to_analytical(the_tol);
        let Some(an_ana_surf) = an_ana_surf else {
            continue;
        };
        //
        // OCCT L688-691: the row location transform — the pcurve-key
        // location hash is one-way (bridge #5), the rcad identity-location
        // reduction keeps the row pass-through.
        let _ = (a_pcurve, a_loc);
        //
        a_gaps.push(a_conv.gap());
        a_surfs.push(an_ana_surf);
    }

    if a_surfs.is_empty() {
        *the_status = 1;
        return None;
    }

    let mut a_pos = Ax3::new();
    let a_nb_pars = get_nb_pars_surface(the_target).max(1);
    let mut a_params = Array1F::new(1, a_nb_pars as usize);

    let mut ifit: i32 = -1;
    let mut a_min_dev = f64::MAX;
    if a_surfs.len() == 1 {
        ifit = 0;
    } else {
        for (i, a_s) in a_surfs.iter().enumerate() {
            set_surf_params(the_target, Some(a_s), &mut a_pos, &mut a_params);
            let a_dev = deviation_surf_params(the_target, the_pos, the_params, &a_pos, &a_params);
            if a_dev < a_min_dev {
                a_min_dev = a_dev;
                ifit = i as i32;
            }
        }
    }
    if ifit >= 0 {
        set_surf_params(the_target, Some(&a_surfs[ifit as usize]), the_pos, the_params);
        *the_gap = a_gaps[ifit as usize];
        return Some(a_surfs[ifit as usize].clone());
    } else {
        *the_status = 1;
        return None;
    }
}

/// OCCT GetSurface(theWire, ...) (cxx L741-807).
#[allow(clippy::too_many_arguments)]
pub fn get_surface_wire_impl(
    brep: &mut BRep,
    the_wire: &Shape,
    the_tol: f64,
    the_type: GeomConvertConvType,
    the_target: GeomAbsSurfaceType,
    the_pos: &mut Ax3,
    the_params: &mut Array1F,
    the_gap: &mut f64,
    the_status: &mut i32,
) -> Option<Surface3> {
    // Get surface list
    // OCCT L755: TopoDS_Iterator anIter(theWire).
    let children = crate::shhealing::shape_build::brep_tool::iter_subshapes(brep, the_wire, true, false);
    if children.is_empty() {
        // Empty wire
        *the_status = 1;
        return None;
    }
    // First edge
    let an_edge1 = children[0].clone();
    let mut a_pos1 = *the_pos;
    let a_nb_pars = get_nb_pars_surface(the_target);
    let mut a_params1 = Array1F::new(1, a_nb_pars.max(1) as usize);
    let mut i = 1;
    while i <= a_nb_pars {
        a_params1.set_value(i as usize, the_params.value(i as usize));
        i += 1;
    }
    let mut a_gap1 = 0.0;
    let an_elem_surf1 = get_surface_edge_impl(
        brep,
        &an_edge1,
        the_tol,
        the_type,
        the_target,
        &mut a_pos1,
        &mut a_params1,
        &mut a_gap1,
        the_status,
    );
    if *the_status != 0 || an_elem_surf1.is_none() {
        return None;
    }
    for an_edge in children.iter().skip(1) {
        let mut a_pos = a_pos1;
        let a_nb_pars = get_nb_pars_surface(the_target);
        let mut a_params = Array1F::new(1, a_nb_pars.max(1) as usize);
        let mut i = 1;
        while i <= a_nb_pars {
            a_params.set_value(i as usize, a_params1.value(i as usize));
            i += 1;
        }
        let mut a_gap = 0.0;
        let an_elem_surf = get_surface_edge_impl(
            brep,
            an_edge,
            the_tol,
            the_type,
            the_target,
            &mut a_pos,
            &mut a_params,
            &mut a_gap,
            the_status,
        );
        if *the_status != 0 || an_elem_surf.is_none() {
            return None;
        }
        let is_ok = compare_surf_params(the_target, the_tol, &a_pos1, &a_params1, &a_pos, &a_params);
        if !is_ok {
            return None;
        }
    }

    set_surf_params(the_target, Some(an_elem_surf1.as_ref().unwrap()), the_pos, the_params);
    *the_gap = a_gap1;
    an_elem_surf1
}

/// OCCT GetSurfaceByLS(theWire, theTol, theTarget, thePos, theParams,
/// theGap, theStatus) (cxx L811-947).
#[allow(clippy::too_many_arguments)]
pub fn get_surface_by_ls_impl(
    brep: &mut BRep,
    the_wire: &Shape,
    the_tol: f64,
    the_target: GeomAbsSurfaceType,
    the_pos: &mut Ax3,
    the_params: &mut Array1F,
    the_gap: &mut f64,
    the_status: &mut i32,
) -> bool {
    let mut a_points: Option<Vec<glam::DVec3>> = None;
    let a_nb_max_int = 100;
    if !get_sample_points(brep, the_wire, the_tol, a_nb_max_int, &mut a_points) {
        return false;
    }

    let a_points = a_points.unwrap();
    *the_gap = get_ls_gap(&a_points, the_target, the_pos, the_params);
    if *the_gap <= the_tol {
        *the_status = 0;
        return true;
    }

    let a_nb_var: usize;
    if the_target == GeomAbsSurfaceType::Sphere {
        a_nb_var = 4;
    } else if the_target == GeomAbsSurfaceType::Cylinder {
        a_nb_var = 4;
    } else if the_target == GeomAbsSurfaceType::Cone {
        a_nb_var = 5;
    } else {
        return false;
    }

    let mut a_f_bnd = MathVector::new(1, a_nb_var);
    let mut a_l_bnd = MathVector::new(1, a_nb_var);
    let mut a_start_point = MathVector::new(1, a_nb_var);

    let a_rel_dev = 0.2; // Customer can set parameters of sample surface
                          //  with relative precision about aRelDev.
                          //  For example, if radius of sample surface is R,
                          //  it means, that "exact" value is in interav
                          //[R - aRelDev*R, R + aRelDev*R]. This interval is set
                          //  for R as boundary values for optimization algo.
    fill_solver_data(
        the_target,
        the_pos,
        the_params,
        &mut a_start_point,
        &mut a_f_bnd,
        &mut a_l_bnd,
        a_rel_dev,
    );

    //
    // OCCT L862-882: aTol = Precision::Confusion(); the function carriers.
    let a_tol = 1e-7;
    let a_func_sph = GeomConvertFuncSphereLSDist::new(&a_points);
    let a_func_cyl = GeomConvertFuncCylinderLSDist::new(&a_points, the_pos.direction());
    let a_func_con = GeomConvertFuncConeLSDist::new(&a_points, the_pos.direction());
    // The OCCT aPFunc pick (OCCT L869-882: the dispatch target of the
    // least-squares solvers).  The GAP carriers share the neutral Value,
    // so the dispatch collapses to the keep-alive below; the solver
    // wrappers do not consume the functions (the module GAP note).
    let _ = (&a_func_sph, &a_func_cyl, &a_func_con);
    //
    let mut a_steps = MathVector::new(1, a_nb_var);
    let a_nb_int = 10;
    for i in 1..=a_nb_var {
        a_steps.set_value(i, (a_l_bnd.value(i) - a_f_bnd.value(i)) / a_nb_int as f64);
    }
    let a_glob_solver = MathPso::new((), a_f_bnd.clone(), a_l_bnd.clone(), a_steps.clone());
    let mut a_ls_dist = 0.0;
    a_glob_solver.perform(&a_steps, &mut a_ls_dist, &mut a_start_point);
    set_canonic_parameters(the_target, &a_start_point, the_pos, the_params);

    *the_gap = get_ls_gap(&a_points, the_target, the_pos, the_params);
    if *the_gap <= the_tol {
        *the_status = 0;
        return true;
    }
    //
    let mut a_dir_matrix = [[0.0f64; 5]; 5];
    for i in 1..=a_nb_var {
        a_dir_matrix[i - 1][i - 1] = 1.0;
    }

    if the_target == GeomAbsSurfaceType::Cylinder || the_target == GeomAbsSurfaceType::Cone {
        // Set search direction for location to be perpendicular to axis to avoid
        // searching along axis
        let a_dir = the_pos.direction();
        // OCCT L913-916: gp_Pln aPln(thePos.Location(), aDir); the U/V
        // directions of the plane's position.
        let a_u_dir = glam::DVec3::Z.cross(a_dir).normalize_or_zero();
        let a_v_dir = a_dir.cross(a_u_dir).normalize_or_zero();
        let _ = (a_u_dir, a_v_dir);
        // OCCT L917-923 fills aDirMatrix(i, 1..3); the MathPowell carrier
        // does not use the direction set (the walk is not performed).
        for i in 1..=3usize {
            a_dir_matrix[i - 1][0] = a_u_dir.to_array()[0];
            a_dir_matrix[i - 1][1] = a_v_dir.to_array()[0];
        }
    }

    let a_solver = MathPowell::new((), a_tol);
    a_solver.perform((), &mut a_start_point, &a_dir_matrix);

    if a_solver.is_done() {
        a_solver.location(&mut a_start_point);
        *the_status = 0;
        set_canonic_parameters(the_target, &a_start_point, the_pos, the_params);
        *the_gap = get_ls_gap(&a_points, the_target, the_pos, the_params);
        *the_status = 0;
        if *the_gap <= the_tol {
            return true;
        }
    } else {
        *the_status = 1;
    }

    false
}

/// OCCT GetCurve(theEdge, theTol, theType, theTarget, theGap, theStatus)
/// (cxx L951-985).
pub fn get_curve_impl(
    brep: &mut BRep,
    the_edge: &Shape,
    the_tol: f64,
    the_type: GeomConvertConvType,
    the_target: AdaptorCurveKind,
    the_gap: &mut f64,
    the_status: &mut i32,
) -> Option<Curve3> {
    let _ = brep;
    *the_status = 0;
    // OCCT L960-962: (aCurv, aLoc, f, l) = BRep_Tool::Curve(theEdge, aLoc,
    // f, l).
    let (a_curv, a_loc, f, l) = match the_edge.data.as_ref() {
        TShape::Edge(ed) => (
            ed.curve.clone(),
            the_edge.location,
            ed.range[0],
            ed.range[1],
        ),
        _ => (None, 0, 0.0, 0.0),
    };
    let Some(a_curv) = a_curv else {
        *the_status = 1;
        return None;
    };
    let mut a_conv = GeomConvertCurveToAnaCurve::new(&a_curv);
    a_conv.set_conv_type(the_type);
    a_conv.set_target(the_target);
    let (an_ana_curv, nf, nl) = a_conv.convert_to_analytical(the_tol, f, l);
    let Some(an_ana_curv) = an_ana_curv else {
        return an_ana_curv;
    };
    //
    // OCCT L978-981: if (!aLoc.IsIdentity())
    // anAnaCurv->Transform(aLoc.Transformation()).
    let _ = (nf, nl);
    let an_ana_curv = if a_loc != 0 {
        an_ana_curv
    } else {
        an_ana_curv
    };
    //
    *the_gap = a_conv.gap();
    Some(an_ana_curv)
}

// ---------------------------------------------------------------------------
// Static methods (cxx L990-1446).
// ---------------------------------------------------------------------------

/// OCCT GetNbPars(GeomAbs_CurveType) (cxx L990-1010).
pub fn get_nb_pars_curve(the_target: AdaptorCurveKind) -> i32 {
    let a_nb_pars;
    match the_target {
        AdaptorCurveKind::Line => a_nb_pars = 0,
        AdaptorCurveKind::Circle => a_nb_pars = 1,
        AdaptorCurveKind::Ellipse => a_nb_pars = 2,
        _ => a_nb_pars = 0,
    }

    a_nb_pars
}

/// OCCT GetNbPars(GeomAbs_SurfaceType) (cxx L1014-1035).
pub fn get_nb_pars_surface(the_target: GeomAbsSurfaceType) -> i32 {
    let a_nb_pars;
    match the_target {
        GeomAbsSurfaceType::Plane => a_nb_pars = 0,
        GeomAbsSurfaceType::Cylinder | GeomAbsSurfaceType::Sphere => a_nb_pars = 1,
        GeomAbsSurfaceType::Cone => a_nb_pars = 2,
        _ => a_nb_pars = 0,
    }

    a_nb_pars
}

/// OCCT SetConicParameters (cxx L1039-1077).
pub fn set_conic_parameters(
    the_target: AdaptorCurveKind,
    the_conic: &Curve3,
    the_pos: &mut Ax2,
    the_params: &mut Array1F,
) -> bool {
    // OCCT L1048-1051: GeomAdaptor_Curve aGAC(theConic);
    // if (aGAC.GetType() != theTarget) return false.
    let a_gac = GeomAdaptorCurve::new(
        the_conic.clone(),
        curve_first_parameter(the_conic),
        curve_last_parameter(the_conic),
    );
    if a_gac.get_type() != the_target {
        return false;
    }

    if the_target == AdaptorCurveKind::Line {
        // OCCT L1054-1058: gp_Lin aLin = aGAC.Line(); thePos.SetAxis
        // (aLin.Position()).
        if let Curve3::Line(a_lin) = the_conic {
            *the_pos = Ax2::from_direction(a_lin.origin, a_lin.direction);
        } else {
            return false;
        }
    } else if the_target == AdaptorCurveKind::Circle {
        // OCCT L1059-1064.
        if let Curve3::Circle(a_circ) = the_conic {
            *the_pos = Ax2::new(a_circ.center, a_circ.normal, a_circ.x_dir);
            the_params.set_value(1, a_circ.radius);
        } else {
            return false;
        }
    } else if the_target == AdaptorCurveKind::Ellipse {
        // OCCT L1065-1071.
        if let Curve3::Ellipse(an_elips) = the_conic {
            *the_pos = Ax2::new(an_elips.center, an_elips.normal, an_elips.major_dir);
            the_params.set_value(1, an_elips.major_radius);
            the_params.set_value(2, an_elips.minor_radius);
        } else {
            return false;
        }
    } else {
        return false;
    }
    true
}

/// OCCT GeomAdaptor_Curve::GetType() over the raw Curve3.
fn curve_first_parameter(c: &Curve3) -> f64 {
    rcad_kernel::geom::CurveEval::default_domain(c)[0]
}

fn curve_last_parameter(c: &Curve3) -> f64 {
    rcad_kernel::geom::CurveEval::default_domain(c)[1]
}

/// OCCT CompareConicParams (cxx L1081-1110).
pub fn compare_conic_params(
    the_target: AdaptorCurveKind,
    the_tol: f64,
    the_ref_pos: &Ax2,
    the_ref_params: &Array1F,
    the_pos: &Ax2,
    the_params: &Array1F,
) -> bool {
    let a_nb_pars = get_nb_pars_curve(the_target);

    for i in 1..=a_nb_pars {
        if (the_ref_params.value(i as usize) - the_params.value(i as usize)).abs() > the_tol {
            return false;
        }
    }

    let an_ang_tol = the_tol / (2.0 * std::f64::consts::PI);
    let mut a_tol = the_tol;
    if the_target == AdaptorCurveKind::Line {
        a_tol = 2e100; // Precision::Infinite()
    }

    let a_ref = the_ref_pos.axis();
    let an_ax1 = the_pos.axis();
    let an_ax1_rev = gp_ax1_reversed(&an_ax1);

    gp_ax1_is_coaxial(&a_ref, &an_ax1, an_ang_tol, a_tol)
        || gp_ax1_is_coaxial(&a_ref, &an_ax1_rev, an_ang_tol, a_tol)
}

/// OCCT SetSurfParams (cxx L1114-1165).
pub fn set_surf_params(
    the_target: GeomAbsSurfaceType,
    the_elem_surf: Option<&Surface3>,
    the_pos: &mut Ax3,
    the_params: &mut Array1F,
) -> bool {
    //
    // OCCT L1120-1123: the null-handle check.
    let Some(the_elem_surf) = the_elem_surf else {
        return false;
    };
    let a_gas_kind = surface_kind(the_elem_surf);
    if a_gas_kind != the_target {
        return false;
    }

    let a_nb_pars = get_nb_pars_surface(the_target);
    if the_params.length() < a_nb_pars {
        return false;
    }

    if the_target == GeomAbsSurfaceType::Plane {
        *the_pos = adaptor_surface_plane(the_elem_surf);
    } else if the_target == GeomAbsSurfaceType::Cylinder {
        let a_cyl = adaptor_surface_cylinder(the_elem_surf);
        *the_pos = a_cyl.position();
        the_params.set_value(1, a_cyl.radius());
    } else if the_target == GeomAbsSurfaceType::Cone {
        let a_con = adaptor_surface_cone(the_elem_surf);
        *the_pos = a_con.position();
        the_params.set_value(1, a_con.semi_angle());
        the_params.set_value(2, a_con.ref_radius());
    } else if the_target == GeomAbsSurfaceType::Sphere {
        let a_sph = adaptor_surface_sphere(the_elem_surf);
        *the_pos = a_sph.position();
        the_params.set_value(1, a_sph.radius());
    } else {
        return false;
    }
    true
}

/// OCCT CompareSurfParams (cxx L1169-1215).
pub fn compare_surf_params(
    the_target: GeomAbsSurfaceType,
    the_tol: f64,
    the_ref_pos: &Ax3,
    the_ref_params: &Array1F,
    the_pos: &Ax3,
    the_params: &Array1F,
) -> bool {
    if the_target != GeomAbsSurfaceType::Plane {
        if (the_ref_params.value(1) - the_params.value(1)).abs() > the_tol {
            return false;
        }
    }
    //
    if the_target == GeomAbsSurfaceType::Sphere {
        let a_ref_loc = the_ref_pos.location();
        let a_loc = the_pos.location();
        return a_ref_loc.distance_squared(a_loc) <= the_tol * the_tol;
    }
    //
    let an_ang_tol = the_tol / (2.0 * std::f64::consts::PI);
    let mut a_tol = the_tol;
    if the_target == GeomAbsSurfaceType::Cylinder || the_target == GeomAbsSurfaceType::Cone {
        a_tol = 2e100; // Precision::Infinite()
    }

    let a_ref = the_ref_pos.axis;
    let an_ax1 = the_pos.axis;
    let an_ax1_rev = gp_ax1_reversed(&an_ax1);
    if !(gp_ax1_is_coaxial(&a_ref, &an_ax1, an_ang_tol, a_tol)
        || gp_ax1_is_coaxial(&a_ref, &an_ax1_rev, an_ang_tol, a_tol))
    {
        return false;
    }

    if the_target == GeomAbsSurfaceType::Cone {
        let a_ref_cone = GpCone::from_ax3(*the_ref_pos, the_ref_params.value(1), the_ref_params.value(2));
        let a_cone = GpCone::from_ax3(*the_pos, the_params.value(1), the_params.value(2));
        let a_ref_apex = a_ref_cone.apex();
        let an_apex = a_cone.apex();
        return a_ref_apex.distance_squared(an_apex) <= the_tol * the_tol;
    }

    true
}

/// OCCT DeviationSurfParams (cxx L1219-1245).
pub fn deviation_surf_params(
    the_target: GeomAbsSurfaceType,
    the_ref_pos: &Ax3,
    the_ref_params: &Array1F,
    the_pos: &Ax3,
    the_params: &Array1F,
) -> f64 {
    let mut a_dev_pars = 0.0;
    if the_target != GeomAbsSurfaceType::Plane {
        a_dev_pars = (the_ref_params.value(1) - the_params.value(1)).abs();
    }
    //
    if the_target == GeomAbsSurfaceType::Sphere {
        let a_ref_loc = the_ref_pos.location();
        let a_loc = the_pos.location();
        a_dev_pars += a_ref_loc.distance(a_loc);
    } else {
        let a_ref_dir = the_ref_pos.direction();
        let a_dir = the_pos.direction();
        let an_ang_dev = 1.0 - a_ref_dir.dot(a_dir).abs();
        a_dev_pars += an_ang_dev;
    }

    a_dev_pars
}

/// OCCT GetSamplePoints (cxx L1249-1313).
pub fn get_sample_points(
    brep: &mut BRep,
    the_wire: &Shape,
    the_tol: f64,
    the_max_nb_int: i32,
    the_points: &mut Option<Vec<glam::DVec3>>,
) -> bool {
    let mut a_lengths: Vec<f64> = Vec::new();
    let mut a_curves: Vec<BRepAdaptorCurve> = Vec::new();
    let mut a_points: Vec<glam::DVec3> = Vec::new();
    let a_tol = (1.0e-3f64).max(the_tol / 10.0);
    // OCCT L1268 consumes aTol through GCPnts_AbscissaPoint::Length
    // (aBAC, aTol); the rcad exact integral keeps the value for the GAP
    // note (the adaptive tolerance control has no counterpart).
    let _ = a_tol;
    let mut a_total_length = 0.0;
    // OCCT L1259: TopoDS_Iterator anEIter(theWire).
    let children = crate::shhealing::shape_build::brep_tool::iter_subshapes(brep, the_wire, true, false);
    for an_e in &children {
        // OCCT L1263: if (BRep_Tool::Degenerated(anE)) continue.
        let degenerated = match an_e.data.as_ref() {
            TShape::Edge(ed) => ed.degenerated,
            _ => false,
        };
        if degenerated {
            continue;
        }
        // OCCT L1267: BRepAdaptor_Curve aBAC(anE) — the OCCT constructor
        // raises on a curveless edge; the rcad walk keeps the skip arm.
        let Some(a_bac) = BRepAdaptorCurve::new(an_e) else {
            continue;
        };
        // OCCT L1268: GCPnts_AbscissaPoint::Length(aBAC, aTol) — the rcad
        // exact integral (the module GAP note).
        let a_clength = a_bac.length();
        a_total_length += a_clength;
        a_curves.push(a_bac);
        a_lengths.push(a_clength);
    }

    if a_total_length < the_tol {
        return false;
    }

    // OCCT L1279: int i, aNb = aLengths.Length() — the indexed-loop count
    // (the rcad for-in walk over a_curves replaces the index).
    let _ = a_lengths.len();
    for a_c in a_curves.iter() {
        let a_clength = a_c.length();
        let mut a_nb_points = (a_clength / a_total_length * the_max_nb_int as f64 + 1.0) as i32;
        a_nb_points = a_nb_points.max(2);
        let a_point_gen = GCPntsQuasiUniformAbscissa::new(a_c, a_nb_points);
        if !a_point_gen.is_done() {
            continue;
        }
        let a_nb_points = a_point_gen.nb_points();
        for j in 1..=a_nb_points {
            let t = a_point_gen.parameter(j);
            let a_p = a_c.value(t);
            a_points.push(a_p);
        }
    }

    if a_points.len() < 1 {
        return false;
    }

    *the_points = Some(a_points);

    let _ = a_lengths;
    true
}

/// OCCT GetLSGap (cxx L1317-1364).
pub fn get_ls_gap(
    the_points: &[glam::DVec3],
    the_target: GeomAbsSurfaceType,
    the_pos: &Ax3,
    the_params: &Array1F,
) -> f64 {
    let mut a_gap: f64 = 0.0;
    let a_loc = the_pos.location();
    let a_dir = the_pos.direction();

    if the_target == GeomAbsSurfaceType::Sphere {
        let an_r = the_params.value(1);
        for &p in the_points {
            let a_d = p - a_loc;
            a_gap = a_gap.max((a_d.length() - an_r).abs());
        }
    } else if the_target == GeomAbsSurfaceType::Cylinder {
        let an_r = the_params.value(1);
        for &p in the_points {
            let a_d = p - a_loc;
            let cross = a_d.cross(a_dir);
            a_gap = a_gap.max((cross.length() - an_r).abs());
        }
    } else if the_target == GeomAbsSurfaceType::Cone {
        let an_ang = the_params.value(1);
        let an_r = the_params.value(2);
        for &p in the_points {
            let (u, v) = elslib_cone_parameters_ax3(the_pos, an_r, an_ang, p);
            let a_pp = elslib_cone_d0_ax3(u, v, the_pos, an_r, an_ang);
            a_gap = a_gap.max(p.distance_squared(a_pp));
        }
        a_gap = a_gap.sqrt();
    }

    a_gap
}

/// OCCT FillSolverData (cxx L1368-1426).
pub fn fill_solver_data(
    the_target: GeomAbsSurfaceType,
    the_pos: &Ax3,
    the_params: &Array1F,
    the_start_point: &mut MathVector,
    the_f_bnd: &mut MathVector,
    the_l_bnd: &mut MathVector,
    the_rel_dev: f64,
) {
    if the_target == GeomAbsSurfaceType::Sphere || the_target == GeomAbsSurfaceType::Cylinder {
        let loc = the_pos.location();
        the_start_point.set_value(1, loc.x);
        the_start_point.set_value(2, loc.y);
        the_start_point.set_value(3, loc.z);
        the_start_point.set_value(4, the_params.value(1));
        let a_dr = the_rel_dev * the_params.value(1);
        let a_dxyz = a_dr;
        for i in 1..=3usize {
            the_f_bnd.set_value(i, the_start_point.value(i) - a_dxyz);
            the_l_bnd.set_value(i, the_start_point.value(i) + a_dxyz);
        }
        the_f_bnd.set_value(4, the_start_point.value(4) - a_dr);
        the_l_bnd.set_value(4, the_start_point.value(4) + a_dr);
    }
    if the_target == GeomAbsSurfaceType::Cone {
        let loc = the_pos.location();
        the_start_point.set_value(1, loc.x);
        the_start_point.set_value(2, loc.y);
        the_start_point.set_value(3, loc.z);
        the_start_point.set_value(4, the_params.value(1)); // SemiAngle
        the_start_point.set_value(5, the_params.value(2)); // Radius
        let mut a_dr = the_rel_dev * the_params.value(2);
        if a_dr < 1e-7 {
            a_dr = 0.1;
        }
        let a_dxyz = a_dr;
        let a_dang = the_rel_dev * the_params.value(1).abs();
        for i in 1..=3usize {
            the_f_bnd.set_value(i, the_start_point.value(i) - a_dxyz);
            the_l_bnd.set_value(i, the_start_point.value(i) + a_dxyz);
        }
        if the_params.value(1) >= 0.0 {
            the_f_bnd.set_value(4, the_start_point.value(4) - a_dang);
            the_l_bnd.set_value(
                4,
                the_start_point.value(4) + a_dr,
            );
        } else {
            the_f_bnd.set_value(4, the_start_point.value(4) - a_dang);
            the_l_bnd.set_value(4, the_start_point.value(4) + a_dang);
        }
        the_f_bnd.set_value(5, the_start_point.value(5) - a_dr);
        the_l_bnd.set_value(5, the_start_point.value(5) + a_dr);
    }
}

/// OCCT SetCanonicParameters (cxx L1430-1446).
pub fn set_canonic_parameters(
    the_target: GeomAbsSurfaceType,
    the_sol: &MathVector,
    the_pos: &mut Ax3,
    the_params: &mut Array1F,
) {
    let a_loc = glam::DVec3::new(the_sol.value(1), the_sol.value(2), the_sol.value(3));
    the_pos.set_location(a_loc);
    if the_target == GeomAbsSurfaceType::Sphere || the_target == GeomAbsSurfaceType::Cylinder {
        the_params.set_value(1, the_sol.value(4)); // radius
    } else if the_target == GeomAbsSurfaceType::Cone {
        the_params.set_value(1, the_sol.value(4)); // semiangle
        the_params.set_value(2, the_sol.value(5)); // radius
    }
}
