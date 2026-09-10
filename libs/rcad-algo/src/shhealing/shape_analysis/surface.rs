//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_Surface`
//! (`ShapeAnalysis_Surface.hxx` + `ShapeAnalysis_Surface.cxx` L1-1946).
//!
//! Tool for analysing surfaces: singularities (cone apex, sphere poles, ...),
//! U/V closedness, point projection (Newton + iso walk) and UV bound boxes.
//!
//! Architecture bridges (numbering follows the W1 edge.rs / curve.rs style):
//! 1. OCCT `occ::handle<Geom_Surface>` -> `Surface3` (clone-on-handle);
//!    `mySurf->Bounds()` -> `SurfaceEval::default_domain`; `IsUClosed()` /
//!    `IsVClosed()` -> the kernel `SurfaceEval` re-hosts.
//! 2. OCCT `GeomAdaptor_Surface` (myAdSur) -> the kernel `Surface3` methods
//!    (`Value` -> `point_at`, `D1`/`D2` -> `derivatives`/`derivatives2`,
//!    `UResolution`/`VResolution`).
//! 3. OCCT `Geom_RectangularTrimmedSurface` -> `Surface3::Trimmed` (the
//!    surftype -> OtherSurface demotion reads the trimmed wrapper).
//! 4. OCCT `Geom_BSplineSurface` accessors -> the kernel `BSplineSurface`
//!    flat-knot vectors: `UKnot(i)` -> the distinct-knot table,
//!    `UMultiplicity(i)` -> the repeat count, `Pole(i,j)` ->
//!    `control_points[i-1][j-1]`; `IsUPeriodic()` -> the SurfaceEval flag.
//! 5. OCCT `Extrema_ExtPS` (myExtPS, the deferred Initialize/Perform pair) ->
//!    the kernel `ExtPS::with_domain` construction at the myExtOK guard and
//!    at each Perform call; `SetFlag(Extrema_ExtFlag_MIN)` is the kernel's
//!    minima-only contract.
//! 6. GAP: `Geom_Surface::UIso/VIso` (the per-surface iso curve
//!    construction, TKG3d) and `Adaptor3d_IsoCurve` (the offset-surface iso
//!    adaptor) have no kernel port; `ComputeIso` yields the null handle and
//!    the OCCT `if (... myIsoUF.IsNull() ...) return theMin;` ("no
//!    isolines") path of UVFromIso stays the live path.
//! 7. OCCT `Bnd_Box` / `BndLib_Add3dCurve::Add` -> the kernel `BndBox` fed
//!    by the base::bnd_lib curve walk; `Bnd_Box::Distance(other)` -> the
//!    local lower-distance re-host.
//! 8. OCCT `try/catch` arms (ComputeIso, ValueOfUV, UVFromIso) -> the rcad
//!    evaluation is panic-free; the OCCT catch results (the null iso, the
//!    `theMin = RealLast()`, the half-bounds S/T fallback) are kept at the
//!    annotated sites.

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema::{ExtPS, POnSurface};
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::precision::{is_infinite_value, CONFUSION, INFINITE_VALUE, PCONFUSION};

use super::analysis::adjust_by_period;
use super::curve::ShapeAnalysisCurve;

// OCCT Standard_Real.hxx L182-185.
const REAL_LAST: f64 = f64::MAX;

/// The GeomAdaptor_Surface::GetType() dispatch over the rcad Surface3
/// (bridge #2) — the surface kinds the switch statements distinguish.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Bezier,
    BSpline,
    Revolution,
    Extrusion,
    Offset,
    Other,
}

fn surf_kind(s: &Surface3) -> SurfKind {
    match s {
        Surface3::Plane(_) => SurfKind::Plane,
        Surface3::Cylinder(_) => SurfKind::Cylinder,
        Surface3::Cone(_) => SurfKind::Cone,
        Surface3::Sphere(_) => SurfKind::Sphere,
        Surface3::Torus(_) => SurfKind::Torus,
        Surface3::Bezier(_) => SurfKind::Bezier,
        Surface3::BSpline(_) => SurfKind::BSpline,
        Surface3::Revolution(_) => SurfKind::Revolution,
        Surface3::LinearExtrusion(_) => SurfKind::Extrusion,
        Surface3::Offset(_) => SurfKind::Offset,
        _ => SurfKind::Other,
    }
}

/// OCCT Precision::IsNegativeInfinite / IsPositiveInfinite (Precision.hxx
/// L350-367).
fn is_neg_inf(v: f64) -> bool {
    v <= -0.5 * INFINITE_VALUE
}

fn is_pos_inf(v: f64) -> bool {
    v >= 0.5 * INFINITE_VALUE
}

/// OCCT RestrictBounds(theFirst, theLast) (cxx L60-78) — replaces the
/// infinite bounds with the +-1000 / +-2000 spans.
fn restrict_bounds_1(the_first: &mut f64, the_last: &mut f64) {
    let is_f_inf = is_neg_inf(*the_first);
    let is_l_inf = is_pos_inf(*the_last);
    if is_f_inf || is_l_inf {
        if is_f_inf && is_l_inf {
            *the_first = -1000.0;
            *the_last = 1000.0;
        } else if is_f_inf {
            *the_first = *the_last - 2000.0;
        } else {
            *the_last = *the_first + 2000.0;
        }
    }
}

/// OCCT RestrictBounds(theUf, theUl, theVf, theVl) (cxx L81-86).
fn restrict_bounds(the_uf: &mut f64, the_ul: &mut f64, the_vf: &mut f64, the_vl: &mut f64) {
    restrict_bounds_1(the_uf, the_ul);
    restrict_bounds_1(the_vf, the_vl);
}

/// The distinct-knot table + multiplicities of a flat knot vector (bridge
/// #4): the indices are 1-based like the OCCT Geom_BSplineSurface knot
/// accessors.
struct KnotTable {
    knots: Vec<f64>,
    mults: Vec<i32>,
}

impl KnotTable {
    fn new(flat: &[f64]) -> Self {
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for k in flat {
            if let Some(last) = knots.last() {
                if *k == *last {
                    *mults.last_mut().unwrap() += 1;
                    continue;
                }
            }
            knots.push(*k);
            mults.push(1);
        }
        KnotTable { knots, mults }
    }
    #[allow(dead_code)]
    fn len(&self) -> i32 {
        self.knots.len() as i32
    }
    fn knot(&self, i: i32) -> f64 {
        self.knots[(i - 1) as usize]
    }
    fn mult(&self, i: i32) -> i32 {
        self.mults[(i - 1) as usize]
    }
}

// ---------------------------------------------------------------------------
// ShapeAnalysis_Surface
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_Surface (hxx L25-180).
pub struct ShapeAnalysisSurface {
    /// OCCT mySurf.
    my_surf: Surface3,
    /// OCCT myUF / myUL / myVF / myVL.
    my_uf: f64,
    my_ul: f64,
    my_vf: f64,
    my_vl: f64,
    /// OCCT myExtOK.
    my_ext_ok: bool,
    /// OCCT myExtPS (the kernel ExtPS re-host; bridge #5).
    my_ext_ps: Option<ExtPS>,
    /// The domain of the myExtPS Initialize call (the kernel ExtPS carries
    /// its domain at construction).
    my_ext_domain: [f64; 4],
    /// OCCT myPreci[4].
    my_preci: [f64; 4],
    /// OCCT myP3d[4].
    my_p3d: [DVec3; 4],
    /// OCCT myFirstP2d[4] / myLastP2d[4].
    my_first_p2d: [DVec2; 4],
    my_last_p2d: [DVec2; 4],
    /// OCCT myFirstPar[4] / myLastPar[4].
    my_first_par: [f64; 4],
    my_last_par: [f64; 4],
    /// OCCT myUIsoDeg[4].
    my_uiso_deg: [bool; 4],
    /// OCCT myNbDeg.
    my_nb_deg: i32,
    /// OCCT myIsoUF / myIsoUL / myIsoVF / myIsoVL.
    my_iso_uf: Option<Curve3>,
    my_iso_ul: Option<Curve3>,
    my_iso_vf: Option<Curve3>,
    my_iso_vl: Option<Curve3>,
    /// OCCT myIsos / myIsoBoxes.
    my_isos: bool,
    my_iso_boxes: bool,
    /// OCCT myBndUF / myBndUL / myBndVF / myBndVL.
    my_bnd_uf: BndBox,
    my_bnd_ul: BndBox,
    my_bnd_vf: BndBox,
    my_bnd_vl: BndBox,
    /// OCCT myGap.
    my_gap: f64,
    /// OCCT myUDelt / myVDelt.
    my_udelt: f64,
    my_vdelt: f64,
    /// OCCT myUCloseVal / myVCloseVal.
    my_uclose_val: f64,
    my_vclose_val: f64,
}

impl ShapeAnalysisSurface {
    /// OCCT ShapeAnalysis_Surface(S) (cxx L99-119).
    pub fn new(s: Surface3) -> Self {
        let mut this = ShapeAnalysisSurface {
            my_surf: s,
            my_uf: 0.0,
            my_ul: 0.0,
            my_vf: 0.0,
            my_vl: 0.0,
            my_ext_ok: false, //: 30
            my_ext_ps: None,
            my_ext_domain: [0.0; 4],
            my_preci: [0.0; 4],
            my_p3d: [DVec3::ZERO; 4],
            my_first_p2d: [DVec2::ZERO; 4],
            my_last_p2d: [DVec2::ZERO; 4],
            my_first_par: [0.0; 4],
            my_last_par: [0.0; 4],
            my_uiso_deg: [false; 4],
            my_nb_deg: -1,
            my_iso_uf: None,
            my_iso_ul: None,
            my_iso_vf: None,
            my_iso_vl: None,
            my_isos: false,
            my_iso_boxes: false,
            my_bnd_uf: BndBox::new(),
            my_bnd_ul: BndBox::new(),
            my_bnd_vf: BndBox::new(),
            my_bnd_vl: BndBox::new(),
            my_gap: 0.0,
            my_udelt: 0.01,
            my_vdelt: 0.01,
            my_uclose_val: -1.0,
            my_vclose_val: -1.0,
        };
        // Bug 33895: the null-surface guard has no rcad counterpart (the
        // value is non-null by construction).
        let dom = SurfaceEval::default_domain(&this.my_surf);
        // mySurf->Bounds(myUF, myUL, myVF, myVL).
        this.my_uf = dom[0];
        this.my_ul = dom[1];
        this.my_vf = dom[2];
        this.my_vl = dom[3];
        // myAdSur = new GeomAdaptor_Surface(mySurf) — bridge #2.
        this
    }

    /// OCCT Surface() — the surface value.
    pub fn surface(&self) -> &Surface3 {
        &self.my_surf
    }

    /// OCCT Gap().
    pub fn gap(&self) -> f64 {
        self.my_gap
    }

    /// OCCT Bounds(U1, U2, V1, V2) — the stored bounds.
    pub fn bounds(&self, u1: &mut f64, u2: &mut f64, v1: &mut f64, v2: &mut f64) {
        *u1 = self.my_uf;
        *u2 = self.my_ul;
        *v1 = self.my_vf;
        *v2 = self.my_vl;
    }

    /// OCCT Value(P2d) — the 3D point at (u, v).
    pub fn value(&self, p2d: DVec2) -> DVec3 {
        self.my_surf.point_at(p2d.x, p2d.y)
    }

    /// OCCT Init(theSurface) (cxx L124-151).
    pub fn init(&mut self, the_surface: Surface3) {
        // OCCT: if (mySurf == theSurface) return; — the handle equality.
        self.my_ext_ok = false; //: 30
        self.my_ext_ps = None;
        self.my_surf = the_surface;
        self.my_nb_deg = -1;
        self.my_uclose_val = -1.0;
        self.my_vclose_val = -1.0;
        self.my_gap = 0.0;
        let dom = SurfaceEval::default_domain(&self.my_surf);
        // mySurf->Bounds(myUF, myUL, myVF, myVL).
        self.my_uf = dom[0];
        self.my_ul = dom[1];
        self.my_vf = dom[2];
        self.my_vl = dom[3];
        self.my_isos = false;
        self.my_iso_boxes = false;
        self.my_iso_uf = None;
        self.my_iso_ul = None;
        self.my_iso_vf = None;
        self.my_iso_vl = None;
    }

    /// OCCT Init(other) (cxx L154-170) — direct transmission of the
    /// singularity data (rln S4135).
    pub fn init_from_other(&mut self, other: &ShapeAnalysisSurface) {
        self.init(other.my_surf.clone());
        self.my_nb_deg = other.my_nb_deg;
        for i in 0..self.my_nb_deg.max(0) as usize {
            self.my_preci[i] = other.my_preci[i];
            self.my_p3d[i] = other.my_p3d[i];
            self.my_first_p2d[i] = other.my_first_p2d[i];
            self.my_last_p2d[i] = other.my_last_p2d[i];
            self.my_first_par[i] = other.my_first_par[i];
            self.my_last_par[i] = other.my_last_par[i];
            self.my_uiso_deg[i] = other.my_uiso_deg[i];
        }
    }

    /// OCCT ComputeSingularities() (cxx L181-296).
    pub fn compute_singularities(&mut self) {
        // rln S4135
        if self.my_nb_deg >= 0 {
            return;
        }
        let su1 = self.my_uf;
        let sv1 = self.my_vf;
        let su2 = self.my_ul;
        let sv2 = self.my_vl;

        self.my_nb_deg = 0; //: r3

        match &self.my_surf {
            Surface3::Cone(conic_s) => {
                // vApex = -RefRadius / sin(SemiAngle).
                let v_apex = -conic_s.radius / conic_s.half_angle_rad.sin();
                self.my_preci[0] = 0.0;
                // Geom_ConicalSurface::Apex() — the reference point sits on
                // the axis at radius RefRadius; the apex is at
                // radius / tan(SemiAngle) along the axis.
                let apex =
                    conic_s.apex + conic_s.axis * (conic_s.radius / conic_s.half_angle_rad.tan());
                self.my_p3d[0] = apex;
                self.my_first_p2d[0] = DVec2::new(su1, v_apex);
                self.my_last_p2d[0] = DVec2::new(su2, v_apex);
                self.my_first_par[0] = su1;
                self.my_last_par[0] = su2;
                self.my_uiso_deg[0] = false;
                self.my_nb_deg = 1;
            }
            Surface3::Torus(toroid_s) => {
                let minor_r = toroid_s.minor_radius;
                let major_r = toroid_s.major_radius;
                // szv#4:S4163:12Mar99 warning - possible div by zero
                let ang = (1.0f64).min(major_r / minor_r).acos();
                self.my_preci[0] = (major_r - minor_r).max(0.0);
                self.my_preci[1] = self.my_preci[0];
                self.my_p3d[0] = self.my_surf.point_at(0.0, std::f64::consts::PI - ang);
                self.my_first_p2d[0] = DVec2::new(su1, std::f64::consts::PI - ang);
                self.my_last_p2d[0] = DVec2::new(su2, std::f64::consts::PI - ang);
                self.my_p3d[1] = self.my_surf.point_at(0.0, std::f64::consts::PI + ang);
                self.my_first_p2d[1] = DVec2::new(su2, std::f64::consts::PI + ang);
                self.my_last_p2d[1] = DVec2::new(su1, std::f64::consts::PI + ang);
                self.my_first_par[0] = su1;
                self.my_first_par[1] = su1;
                self.my_last_par[0] = su2;
                self.my_last_par[1] = su2;
                self.my_uiso_deg[0] = false;
                self.my_uiso_deg[1] = false;
                self.my_nb_deg = if major_r > minor_r { 1 } else { 2 };
            }
            Surface3::Sphere(_) => {
                self.my_preci[0] = 0.0;
                self.my_preci[1] = 0.0;
                self.my_p3d[0] = self.my_surf.point_at(su1, sv2); // Northern pole is first
                self.my_p3d[1] = self.my_surf.point_at(su1, sv1);
                self.my_first_p2d[0] = DVec2::new(su2, sv2);
                self.my_last_p2d[0] = DVec2::new(su1, sv2);
                self.my_first_p2d[1] = DVec2::new(su1, sv1);
                self.my_last_p2d[1] = DVec2::new(su2, sv1);
                self.my_first_par[0] = su1;
                self.my_first_par[1] = su1;
                self.my_last_par[0] = su2;
                self.my_last_par[1] = su2;
                self.my_uiso_deg[0] = false;
                self.my_uiso_deg[1] = false;
                self.my_nb_deg = 2;
            }
            _ => {
                // Geom_BoundedSurface || Geom_SurfaceOfRevolution ||
                // Geom_OffsetSurface (rln S4135): the rcad Surface3 kinds
                // with bounded parameter spaces.
                if matches!(
                    self.my_surf,
                    Surface3::BSpline(_)
                        | Surface3::Bezier(_)
                        | Surface3::Trimmed(_)
                        | Surface3::Coons(_)
                        | Surface3::Ruled(_)
                        | Surface3::Revolution(_)
                        | Surface3::Offset(_)
                ) {
                    self.my_p3d[0] = self.my_surf.point_at(su1, 0.5 * (sv1 + sv2));
                    self.my_first_p2d[0] = DVec2::new(su1, sv2);
                    self.my_last_p2d[0] = DVec2::new(su1, sv1);

                    self.my_p3d[1] = self.my_surf.point_at(su2, 0.5 * (sv1 + sv2));
                    self.my_first_p2d[1] = DVec2::new(su2, sv1);
                    self.my_last_p2d[1] = DVec2::new(su2, sv2);

                    self.my_p3d[2] = self.my_surf.point_at(0.5 * (su1 + su2), sv1);
                    self.my_first_p2d[2] = DVec2::new(su1, sv1);
                    self.my_last_p2d[2] = DVec2::new(su2, sv1);

                    self.my_p3d[3] = self.my_surf.point_at(0.5 * (su1 + su2), sv2);
                    self.my_first_p2d[3] = DVec2::new(su2, sv2);
                    self.my_last_p2d[3] = DVec2::new(su1, sv2);

                    self.my_first_par[0] = sv1;
                    self.my_first_par[1] = sv1;
                    self.my_last_par[0] = sv2;
                    self.my_last_par[1] = sv2;
                    self.my_uiso_deg[0] = true;
                    self.my_uiso_deg[1] = true;

                    self.my_first_par[2] = su1;
                    self.my_first_par[3] = su1;
                    self.my_last_par[2] = su2;
                    self.my_last_par[3] = su2;
                    self.my_uiso_deg[2] = false;
                    self.my_uiso_deg[3] = false;

                    let corner1 = self.my_surf.point_at(su1, sv1);
                    let corner2 = self.my_surf.point_at(su1, sv2);
                    let corner3 = self.my_surf.point_at(su2, sv1);
                    let corner4 = self.my_surf.point_at(su2, sv2);

                    self.my_preci[0] = corner1.distance(corner2).max(
                        self.my_p3d[0]
                            .distance(corner1)
                            .max(self.my_p3d[0].distance(corner2)),
                    );
                    self.my_preci[1] = corner3.distance(corner4).max(
                        self.my_p3d[1]
                            .distance(corner3)
                            .max(self.my_p3d[1].distance(corner4)),
                    );
                    self.my_preci[2] = corner1.distance(corner3).max(
                        self.my_p3d[2]
                            .distance(corner1)
                            .max(self.my_p3d[2].distance(corner3)),
                    );
                    self.my_preci[3] = corner2.distance(corner4).max(
                        self.my_p3d[3]
                            .distance(corner2)
                            .max(self.my_p3d[3].distance(corner4)),
                    );

                    self.my_nb_deg = 4;
                }
            }
        }
        self.sort_singularities();
    }

    /// OCCT HasSingularities(preci) (cxx L298-300).
    pub fn has_singularities(&mut self, preci: f64) -> bool {
        self.nb_singularities(preci) > 0
    }

    /// OCCT NbSingularities(preci) (cxx L302-313).
    pub fn nb_singularities(&mut self, preci: f64) -> i32 {
        if self.my_nb_deg < 0 {
            self.compute_singularities();
        }
        let mut nb = 0;
        for i in 1..=self.my_nb_deg {
            if self.my_preci[(i - 1) as usize] <= preci {
                nb += 1;
            }
        }
        nb
    }

    /// OCCT Singularity(num, preci, P3d, firstP2d, lastP2d, firstpar,
    /// lastpar, uisodeg) (cxx L324-345).
    #[allow(clippy::too_many_arguments)]
    pub fn singularity(
        &mut self,
        num: i32,
        preci: &mut f64,
        p3d: &mut DVec3,
        first_p2d: &mut DVec2,
        last_p2d: &mut DVec2,
        firstpar: &mut f64,
        lastpar: &mut f64,
        uisodeg: &mut bool,
    ) -> bool {
        //  ATTENTION, les champs sont des tableaux C, n0s partent de 0. num
        //  part de 1
        if self.my_nb_deg < 0 {
            self.compute_singularities();
        }
        if num < 1 || num > self.my_nb_deg {
            return false;
        }
        let i = (num - 1) as usize;
        *p3d = self.my_p3d[i];
        *preci = self.my_preci[i];
        *first_p2d = self.my_first_p2d[i];
        *last_p2d = self.my_last_p2d[i];
        *firstpar = self.my_first_par[i];
        *lastpar = self.my_last_par[i];
        *uisodeg = self.my_uiso_deg[i];
        true
    }

    /// OCCT IsDegenerated(P3d, preci) (cxx L354-367).
    pub fn is_degenerated(&mut self, p3d: DVec3, preci: f64) -> bool {
        if self.my_nb_deg < 0 {
            self.compute_singularities();
        }
        let mut i = 0usize;
        while i < self.my_nb_deg.max(0) as usize && self.my_preci[i] <= preci {
            self.my_gap = self.my_p3d[i].distance(p3d);
            // rln S4135
            if self.my_gap <= preci {
                return true;
            }
            i += 1;
        }
        false
    }

    /// OCCT DegeneratedValues(P3d, preci, firstP2d, lastP2d, firstPar,
    /// lastPar, forward) (cxx L374-412).
    #[allow(clippy::too_many_arguments)]
    pub fn degenerated_values(
        &mut self,
        p3d: DVec3,
        preci: f64,
        first_p2d: &mut DVec2,
        last_p2d: &mut DVec2,
        first_par: &mut f64,
        last_par: &mut f64,
        _forward: bool,
    ) -> bool {
        if self.my_nb_deg < 0 {
            self.compute_singularities();
        }
        // #77 rln S4135: returning singularity which has minimum gap between
        // singular point and input 3D point
        let mut ind_min: i32 = -1;
        let mut gap_min = REAL_LAST;
        let mut i = 0usize;
        while i < self.my_nb_deg.max(0) as usize && self.my_preci[i] <= preci {
            self.my_gap = self.my_p3d[i].distance(p3d);
            // rln S4135
            if self.my_gap <= preci && gap_min > self.my_gap {
                gap_min = self.my_gap;
                ind_min = i as i32;
            }
            i += 1;
        }
        if ind_min >= 0 {
            let ind = ind_min as usize;
            self.my_gap = gap_min;
            *first_p2d = self.my_first_p2d[ind];
            *last_p2d = self.my_last_p2d[ind];
            *first_par = self.my_first_par[ind];
            *last_par = self.my_last_par[ind];
            return true;
        }
        false
    }

    /// OCCT ProjectDegenerated(P3d, preci, neighbour, result)
    /// (cxx L417-452).
    pub fn project_degenerated(
        &mut self,
        p3d: DVec3,
        preci: f64,
        neighbour: DVec2,
        result: &mut DVec2,
    ) -> bool {
        if self.my_nb_deg < 0 {
            self.compute_singularities();
        }
        // added by rln on 03/12/97
        //: c1 abv 23 Feb 98: preci (3d) -> Resolution (2d)
        // #77 rln S4135
        let mut ind_min: i32 = -1;
        let mut gap_min = REAL_LAST;
        let mut i = 0usize;
        while i < self.my_nb_deg.max(0) as usize && self.my_preci[i] <= preci {
            let mut gap2 = self.my_p3d[i].distance_squared(p3d);
            if gap2 > preci * preci {
                gap2 = gap2.min(self.my_p3d[i].distance_squared(self.value(*result)));
            }
            // rln S4135
            if gap2 <= preci * preci && gap_min > gap2 {
                gap_min = gap2;
                ind_min = i as i32;
            }
            i += 1;
        }
        if ind_min < 0 {
            return false;
        }
        let ind = ind_min as usize;
        self.my_gap = gap_min.sqrt();
        if !self.my_uiso_deg[ind] {
            result.x = neighbour.x;
        } else {
            result.y = neighbour.y;
        }
        true
    }

    /// OCCT ProjectDegenerated(nbrPnt, points, pnt2d, preci, direct)
    /// (cxx L464-545) — pdn %12 11.02.99 PRO9234 entity 15402.
    pub fn project_degenerated_seq(
        &mut self,
        nbr_pnt: i32,
        points: &[DVec3],
        pnt2d: &mut [DVec2],
        preci: f64,
        direct: bool,
    ) -> bool {
        if self.my_nb_deg < 0 {
            self.compute_singularities();
        }

        let step: i32 = if direct { 1 } else { -1 };
        // #77 rln S4135
        let mut ind_min: i32 = -1;
        let mut gap_min = REAL_LAST;
        let prec2 = preci * preci;
        let mut j: i32 = if direct { 1 } else { nbr_pnt };
        let mut i = 0usize;
        while i < self.my_nb_deg.max(0) as usize && self.my_preci[i] <= preci {
            let mut gap2 = self.my_p3d[i].distance_squared(points[(j - 1) as usize]);
            if gap2 > prec2 {
                gap2 = gap2
                    .min(self.my_p3d[i].distance_squared(self.value(pnt2d[(j - 1) as usize])));
            }
            if gap2 <= prec2 && gap_min > gap2 {
                gap_min = gap2;
                ind_min = i as i32;
            }
            i += 1;
        }
        if ind_min < 0 {
            return false;
        }
        let ind_min = ind_min as usize;

        self.my_gap = gap_min.sqrt();

        let mut k: i32 = j + step;
        while k <= nbr_pnt && k >= 1 {
            let pk = pnt2d[(k - 1) as usize];
            let p1 = points[(k - 1) as usize];
            if self.my_p3d[ind_min].distance_squared(p1) > prec2
                && self.my_p3d[ind_min].distance_squared(self.value(pk)) > prec2
            {
                break;
            }
            k += step;
        }

        //: p8 abv 11 Mar 99: PRO7226 #489490: if whole pcurve is degenerate,
        // distribute evenly
        if k < 1 || k > nbr_pnt {
            let x1 = if self.my_uiso_deg[ind_min] {
                pnt2d[0].y
            } else {
                pnt2d[0].x
            };
            let x2 = if self.my_uiso_deg[ind_min] {
                pnt2d[(nbr_pnt - 1) as usize].y
            } else {
                pnt2d[(nbr_pnt - 1) as usize].x
            };
            for jj in 1..=nbr_pnt {
                // szv#4:S4163:12Mar99 warning - possible div by zero
                let x = (x1 * ((nbr_pnt - jj) as f64) + x2 * ((jj - 1) as f64))
                    / ((nbr_pnt - 1) as f64);
                if !self.my_uiso_deg[ind_min] {
                    pnt2d[(jj - 1) as usize].x = x;
                } else {
                    pnt2d[(jj - 1) as usize].y = x;
                }
            }
            return true;
        }

        let pk = pnt2d[(k - step - 1) as usize];
        j = k - step;
        while j <= nbr_pnt && j >= 1 {
            if !self.my_uiso_deg[ind_min] {
                pnt2d[(j - 1) as usize].x = pk.x;
            } else {
                pnt2d[(j - 1) as usize].y = pk.y;
            }
            j -= step;
        }
        true
    }

    /// OCCT IsDegenerated(p2d1, p2d2, tol, ratio) (cxx L549-575).
    pub fn is_degenerated_uv(&self, p2d1: DVec2, p2d2: DVec2, tol: f64, ratio: f64) -> bool {
        let p1 = self.value(p2d1);
        let p2 = self.value(p2d2);
        let pm = self.value(0.5 * (p2d1 + p2d2));
        let mut max3d = p1.distance(p2).max(pm.distance(p1).max(pm.distance(p2)));
        if max3d > tol {
            return false;
        }

        // GeomAdaptor_Surface& SA = *Adaptor3d().
        let ru = self.my_surf.u_resolution(1.0);
        let rv = self.my_surf.v_resolution(1.0);

        if ru < PCONFUSION || rv < PCONFUSION {
            return false;
        }
        let du = (p2d1.x - p2d2.x).abs() / ru;
        let dv = (p2d1.y - p2d2.y).abs() / rv;
        max3d *= ratio;
        du * du + dv * dv > max3d * max3d
    }

    // -- Isos ---------------------------------------------------------------

    /// OCCT ComputeBoundIsos() (cxx L612-622).
    pub fn compute_bound_isos(&mut self) {
        if self.my_isos {
            return;
        }
        self.my_isos = true;
        let (uf, ul, vf, vl) = (self.my_uf, self.my_ul, self.my_vf, self.my_vl);
        self.my_iso_uf = compute_iso(&self.my_surf, true, uf);
        self.my_iso_ul = compute_iso(&self.my_surf, true, ul);
        self.my_iso_vf = compute_iso(&self.my_surf, false, vf);
        self.my_iso_vl = compute_iso(&self.my_surf, false, vl);
    }

    /// OCCT UIso(U) (cxx L627-637).
    pub fn uiso(&mut self, u: f64) -> Option<Curve3> {
        if u == self.my_uf {
            self.compute_bound_isos();
            return self.my_iso_uf.clone();
        }
        if u == self.my_ul {
            self.compute_bound_isos();
            return self.my_iso_ul.clone();
        }
        compute_iso(&self.my_surf, true, u)
    }

    /// OCCT VIso(V) (cxx L644-654).
    pub fn viso(&mut self, v: f64) -> Option<Curve3> {
        if v == self.my_vf {
            self.compute_bound_isos();
            return self.my_iso_vf.clone();
        }
        if v == self.my_vl {
            self.compute_bound_isos();
            return self.my_iso_vl.clone();
        }
        compute_iso(&self.my_surf, false, v)
    }

    // -- Closedness ---------------------------------------------------------

    /// OCCT IsUClosed(preci) (cxx L661-866).
    pub fn is_u_closed(&mut self, preci: f64) -> bool {
        let prec = preci.max(CONFUSION);
        let mut an_umid_val = -1.0f64;
        if self.my_uclose_val < 0.0 {
            //    Faut calculer : calculs minimaux
            let (mut uf, mut ul, mut vf, mut vl) = (0.0, 0.0, 0.0, 0.0);
            self.bounds(&mut uf, &mut ul, &mut vf, &mut vl);
            restrict_bounds(&mut uf, &mut ul, &mut vf, &mut vl);
            self.my_udelt = (ul - uf).abs() / 20.0;
            if SurfaceEval::is_u_closed(&self.my_surf) {
                self.my_uclose_val = 0.0;
                self.my_udelt = 0.0;
                self.my_gap = 0.0;
                return true;
            }

            // Calculs adaptes
            // #67 rln S4135
            let mut surftype = surf_kind(&self.my_surf);
            if matches!(self.my_surf, Surface3::Trimmed(_)) {
                surftype = SurfKind::Other; // bridge #3
            }

            match surftype {
                SurfKind::Plane => {
                    self.my_uclose_val = REAL_LAST;
                }
                SurfKind::Extrusion => {
                    //: c8 abv 03 Mar 98: UKI60094 #753
                    let Surface3::LinearExtrusion(extr) = &self.my_surf else {
                        unreachable!()
                    };
                    let dom = CurveEval::default_domain(&*extr.profile);
                    let (f, l) = (dom[0], dom[1]);
                    //: r3 abv (smh) 30 Mar 99: protect against unexpected signals
                    if !is_infinite_value(f) && !is_infinite_value(l) {
                        let p1 = CurveEval::point_at(&*extr.profile, f);
                        let p2 = CurveEval::point_at(&*extr.profile, l);
                        self.my_uclose_val = p1.distance_squared(p2);
                        let pm = CurveEval::point_at(&*extr.profile, (f + l) / 2.0);
                        an_umid_val = p1.distance_squared(pm);
                    } else {
                        self.my_uclose_val = REAL_LAST;
                    }
                }
                SurfKind::BSpline => {
                    let Surface3::BSpline(bs) = &self.my_surf else {
                        unreachable!()
                    };
                    let u_tab = KnotTable::new(&bs.knots_u);
                    let v_tab = KnotTable::new(&bs.knots_v);
                    let nbup = bs.control_points.len() as i32;
                    let mut distmin = REAL_LAST;
                    if SurfaceEval::is_u_periodic(&self.my_surf) {
                        self.my_uclose_val = 0.0;
                        self.my_udelt = 0.0;
                    } else if nbup < 3 {
                        // modified by rln on 12/11/97
                        self.my_uclose_val = REAL_LAST;
                    } else if is_u_rational(bs)
                        || u_tab.mult(1) != bs.degree_u as i32 + 1
                        || u_tab.mult(u_tab.len()) != bs.degree_u as i32 + 1
                    {
                        // #6 //:h4
                        let nbvk = v_tab.len();
                        let mut v = v_tab.knot(1);
                        let mut p1 = self.my_surf.point_at(uf, v);
                        let mut p2 = self.my_surf.point_at(ul, v);
                        self.my_uclose_val = p1.distance_squared(p2);
                        let mut pm = self.my_surf.point_at((uf + ul) / 2.0, v);
                        an_umid_val = p1.distance_squared(pm);
                        distmin = self.my_uclose_val;
                        for i in 2..=nbvk {
                            v = 0.5 * (v_tab.knot(i - 1) + v_tab.knot(i));
                            p1 = self.my_surf.point_at(uf, v);
                            p2 = self.my_surf.point_at(ul, v);
                            let a_dist = p1.distance_squared(p2);
                            if a_dist > self.my_uclose_val {
                                self.my_uclose_val = a_dist;
                                pm = self.my_surf.point_at((uf + ul) / 2.0, v);
                                an_umid_val = p1.distance_squared(pm);
                            } else {
                                distmin = distmin.min(a_dist);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.my_udelt =
                            self.my_udelt.min(0.5 * self.my_surf.u_resolution(distmin)); // #4 smh
                    } else {
                        let nbvp = v_tab.len();
                        let pole =
                            |i: i32, j: i32| bs.control_points[(i - 1) as usize][(j - 1) as usize];
                        self.my_uclose_val = pole(1, 1).distance_squared(pole(nbup, 1));
                        an_umid_val = pole(1, 1).distance_squared(pole(nbup / 2 + 1, 1));
                        distmin = self.my_uclose_val;
                        for i in 2..=nbvp {
                            let a_dist = pole(1, i).distance_squared(pole(nbup, i));
                            if a_dist > self.my_uclose_val {
                                self.my_uclose_val = a_dist;
                                an_umid_val = pole(1, i).distance_squared(pole(nbup / 2 + 1, i));
                            } else {
                                distmin = distmin.min(a_dist);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.my_udelt =
                            self.my_udelt.min(0.5 * self.my_surf.u_resolution(distmin)); // #4 smh
                    }
                }
                SurfKind::Bezier => {
                    let Surface3::Bezier(bz) = &self.my_surf else {
                        unreachable!()
                    };
                    let nbup = bz.control_points.len() as i32;
                    let mut distmin = REAL_LAST;
                    if nbup < 3 {
                        self.my_uclose_val = REAL_LAST;
                    } else {
                        let nbvp = bz.control_points[0].len() as i32;
                        let pole =
                            |i: i32, j: i32| bz.control_points[(i - 1) as usize][(j - 1) as usize];
                        self.my_uclose_val = pole(1, 1).distance_squared(pole(nbup, 1));
                        an_umid_val = pole(1, 1).distance_squared(pole(nbup / 2 + 1, 1));
                        distmin = self.my_uclose_val;
                        for i in 1..=nbvp {
                            let a_dist = pole(1, i).distance_squared(pole(nbup, i));
                            if a_dist > self.my_uclose_val {
                                self.my_uclose_val = a_dist;
                                an_umid_val = pole(1, i).distance_squared(pole(nbup / 2 + 1, i));
                            } else {
                                distmin = distmin.min(a_dist);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.my_udelt =
                            self.my_udelt.min(0.5 * self.my_surf.u_resolution(distmin)); // #4 smh
                    }
                }
                _ => {
                    // Geom_RectangularTrimmedSurface and Geom_OffsetSurface
                    let mut distmin = REAL_LAST;
                    let nbpoints = 101; // can be revised
                    let mut p1 = self.my_surf.point_at(uf, vf);
                    let mut p2 = self.my_surf.point_at(ul, vf);
                    self.my_uclose_val = p1.distance_squared(p2);
                    let mut pm = self.my_surf.point_at((uf + ul) / 2.0, vf);
                    an_umid_val = p1.distance_squared(pm);
                    distmin = self.my_uclose_val;
                    for i in 1..nbpoints {
                        let vparam = vf + (vl - vf) * i as f64 / (nbpoints - 1) as f64;
                        p1 = self.my_surf.point_at(uf, vparam);
                        p2 = self.my_surf.point_at(ul, vparam);
                        let a_dist = p1.distance_squared(p2);
                        if a_dist > self.my_uclose_val {
                            self.my_uclose_val = a_dist;
                            pm = self.my_surf.point_at((uf + ul) / 2.0, vparam);
                            an_umid_val = p1.distance_squared(pm);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    distmin = distmin.sqrt();
                    self.my_udelt = self.my_udelt.min(0.5 * self.my_surf.u_resolution(distmin)); // #4 smh
                }
            } // switch
            self.my_gap = self.my_uclose_val.sqrt();
            self.my_uclose_val = self.my_gap;
        }

        if an_umid_val > 0.0 && self.my_uclose_val > an_umid_val.sqrt() {
            self.my_uclose_val = REAL_LAST;
            return false;
        }

        self.my_uclose_val <= prec
    }

    /// OCCT IsVClosed(preci) (cxx L868-1075) — the V mirror of IsUClosed.
    pub fn is_v_closed(&mut self, preci: f64) -> bool {
        let prec = preci.max(CONFUSION);
        let mut a_vmid_val = -1.0f64;
        if self.my_vclose_val < 0.0 {
            let (mut uf, mut ul, mut vf, mut vl) = (0.0, 0.0, 0.0, 0.0);
            self.bounds(&mut uf, &mut ul, &mut vf, &mut vl);
            restrict_bounds(&mut uf, &mut ul, &mut vf, &mut vl);
            self.my_vdelt = (vl - vf).abs() / 20.0;
            if SurfaceEval::is_v_closed(&self.my_surf) {
                self.my_vclose_val = 0.0;
                self.my_vdelt = 0.0;
                self.my_gap = 0.0;
                return true;
            }

            //    Calculs adaptes
            // #67 rln S4135
            let mut surftype = surf_kind(&self.my_surf);
            if matches!(self.my_surf, Surface3::Trimmed(_)) {
                surftype = SurfKind::Other; // bridge #3
            }

            match surftype {
                SurfKind::Plane
                | SurfKind::Cone
                | SurfKind::Cylinder
                | SurfKind::Sphere
                | SurfKind::Extrusion => {
                    self.my_vclose_val = REAL_LAST;
                }
                SurfKind::Revolution => {
                    let Surface3::Revolution(revol) = &self.my_surf else {
                        unreachable!()
                    };
                    let crv = &revol.profile;
                    let dom = CurveEval::default_domain(&**crv);
                    let p1 = CurveEval::point_at(&**crv, dom[0]);
                    let p2 = CurveEval::point_at(&**crv, dom[1]);
                    self.my_vclose_val = p1.distance_squared(p2);
                }
                SurfKind::BSpline => {
                    let Surface3::BSpline(bs) = &self.my_surf else {
                        unreachable!()
                    };
                    let u_tab = KnotTable::new(&bs.knots_u);
                    let v_tab = KnotTable::new(&bs.knots_v);
                    let nbvp = bs.control_points[0].len() as i32;
                    let mut distmin = REAL_LAST;
                    if SurfaceEval::is_v_periodic(&self.my_surf) {
                        self.my_vclose_val = 0.0;
                        self.my_vdelt = 0.0;
                    } else if nbvp < 3 {
                        // modified by rln on 12/11/97
                        self.my_vclose_val = REAL_LAST;
                    } else if is_v_rational(bs)
                        || v_tab.mult(1) != bs.degree_v as i32 + 1
                        || v_tab.mult(v_tab.len()) != bs.degree_v as i32 + 1
                    {
                        // #6 //:h4
                        let nbuk = u_tab.len();
                        let mut u = u_tab.knot(1);
                        let mut p1 = self.my_surf.point_at(u, vf);
                        let mut p2 = self.my_surf.point_at(u, vl);
                        self.my_vclose_val = p1.distance_squared(p2);
                        let mut pm = self.my_surf.point_at(u, (vf + vl) / 2.0);
                        a_vmid_val = p1.distance_squared(pm);
                        distmin = self.my_vclose_val;
                        for i in 2..=nbuk {
                            u = 0.5 * (u_tab.knot(i - 1) + u_tab.knot(i));
                            p1 = self.my_surf.point_at(u, vf);
                            p2 = self.my_surf.point_at(u, vl);
                            let a_dist = p1.distance_squared(p2);
                            if a_dist > self.my_vclose_val {
                                self.my_vclose_val = a_dist;
                                pm = self.my_surf.point_at(u, (vf + vl) / 2.0);
                                a_vmid_val = p1.distance_squared(pm);
                            } else {
                                distmin = distmin.min(a_dist);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.my_vdelt =
                            self.my_vdelt.min(0.5 * self.my_surf.v_resolution(distmin)); // #4 smh
                    } else {
                        let nbup = bs.control_points.len() as i32;
                        let pole =
                            |i: i32, j: i32| bs.control_points[(i - 1) as usize][(j - 1) as usize];
                        self.my_vclose_val = pole(1, 1).distance_squared(pole(1, nbvp));
                        a_vmid_val = pole(1, 1).distance_squared(pole(1, nbvp / 2 + 1));
                        distmin = self.my_vclose_val;
                        for i in 2..=nbup {
                            let a_dist = pole(i, 1).distance_squared(pole(i, nbvp));
                            if a_dist > self.my_vclose_val {
                                self.my_vclose_val = a_dist;
                                a_vmid_val = pole(i, 1).distance_squared(pole(i, nbvp / 2 + 1));
                            } else {
                                distmin = distmin.min(a_dist);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.my_vdelt =
                            self.my_vdelt.min(0.5 * self.my_surf.v_resolution(distmin)); // #4 smh
                    }
                }
                SurfKind::Bezier => {
                    let Surface3::Bezier(bz) = &self.my_surf else {
                        unreachable!()
                    };
                    let nbvp = bz.control_points[0].len() as i32;
                    let mut distmin = REAL_LAST;
                    if nbvp < 3 {
                        self.my_vclose_val = REAL_LAST;
                    } else {
                        let nbup = bz.control_points.len() as i32;
                        let pole =
                            |i: i32, j: i32| bz.control_points[(i - 1) as usize][(j - 1) as usize];
                        self.my_vclose_val = pole(1, 1).distance_squared(pole(1, nbvp));
                        a_vmid_val = pole(1, 1).distance_squared(pole(1, nbvp / 2 + 1));
                        distmin = self.my_vclose_val;
                        for i in 2..=nbup {
                            let a_dist = pole(i, 1).distance_squared(pole(i, nbvp));
                            if a_dist > self.my_vclose_val {
                                self.my_vclose_val = a_dist;
                                a_vmid_val = pole(i, 1).distance_squared(pole(i, nbvp / 2 + 1));
                            } else {
                                distmin = distmin.min(a_dist);
                            }
                        }
                        distmin = distmin.sqrt();
                        self.my_vdelt =
                            self.my_vdelt.min(0.5 * self.my_surf.v_resolution(distmin)); // #4 smh
                    }
                }
                _ => {
                    // Geom_RectangularTrimmedSurface and Geom_OffsetSurface
                    let mut distmin = REAL_LAST;
                    let nbpoints = 101; // can be revised
                    let mut p1 = self.my_surf.point_at(uf, vf);
                    let mut p2 = self.my_surf.point_at(uf, vl);
                    self.my_vclose_val = p1.distance_squared(p2);
                    let mut pm = self.my_surf.point_at(uf, (vf + vl) / 2.0);
                    a_vmid_val = p1.distance_squared(pm);
                    distmin = self.my_vclose_val;
                    for i in 1..nbpoints {
                        let uparam = uf + (ul - uf) * i as f64 / (nbpoints - 1) as f64;
                        p1 = self.my_surf.point_at(uparam, vf);
                        p2 = self.my_surf.point_at(uparam, vl);
                        let a_dist = p1.distance_squared(p2);
                        if a_dist > self.my_vclose_val {
                            self.my_vclose_val = a_dist;
                            pm = self.my_surf.point_at(uparam, (vf + vl) / 2.0);
                            a_vmid_val = p1.distance_squared(pm);
                        } else {
                            distmin = distmin.min(a_dist);
                        }
                    }
                    distmin = distmin.sqrt();
                    self.my_vdelt = self.my_vdelt.min(0.5 * self.my_surf.v_resolution(distmin)); // #4 smh
                }
            } // switch
            self.my_gap = self.my_vclose_val.sqrt();
            self.my_vclose_val = self.my_gap;
        }

        if a_vmid_val > 0.0 && self.my_vclose_val > a_vmid_val.sqrt() {
            self.my_vclose_val = REAL_LAST;
            return false;
        }

        self.my_vclose_val <= prec
    }

    // -- Projection ---------------------------------------------------------

    /// OCCT SurfaceNewton(p2dPrev, P3D, preci, sol) (cxx L1065-1143) —
    /// Newton algo (S4030).
    pub fn surface_newton(&self, p2d_prev: DVec2, p3d: DVec3, preci: f64, sol: &mut DVec2) -> i32 {
        let (mut uf, mut ul, mut vf, mut vl) = (0.0, 0.0, 0.0, 0.0);
        self.bounds(&mut uf, &mut ul, &mut vf, &mut vl);
        let du = self.my_surf.u_resolution(preci);
        let dv = self.my_surf.v_resolution(preci);
        let uf = uf - du;
        let ul = ul + du;
        let vf = vf - dv;
        let vl = vl + dv;

        let tol = CONFUSION;
        let tol2 = tol * tol;
        let mut u = p2d_prev.x;
        let mut v = p2d_prev.y;
        let mut du;
        let mut dv;
        let rsfirst = p3d - self.value(DVec2::new(u, v)); // pdn
        for _ in 0..25 {
            // SurfAdapt.D2(U, V, pnt, ru, rv, ruu, rvv, ruv).
            let (pnt, ru, rv, ruu, rvv, ruv) = self.my_surf.derivatives2(u, v);

            // normal
            let ru2 = ru.dot(ru);
            let rv2 = rv.dot(rv);
            let n = ru.cross(rv);
            let nrm2 = n.length_squared();
            if nrm2 < 1e-10 || is_pos_inf(nrm2) {
                break; // n == 0, use standard
            }

            // descriminant
            let rs = p3d - pnt;
            let r_suu = rs.dot(ruu);
            let r_svv = rs.dot(rvv);
            let r_suv = rs.dot(ruv);
            let d =
                -nrm2 + rv2 * r_suu + ru2 * r_svv - 2.0 * r_suv * ru.dot(rv) + r_suv * r_suv
                    - r_suu * r_svv;
            if d.abs() < 1e-10 {
                break; // bad case; use standard
            }

            // compute step
            let fract = 1.0 / d;
            du = rs.dot(n.cross(rv) + ru * r_svv - rv * r_suv) * fract;
            dv = rs.dot(ru.cross(n) + rv * r_suu - ru * r_suv) * fract;
            u += du;
            v += dv;
            if u < uf || u > ul || v < vf || v > vl {
                break;
            }
            // test the step by uv and deviation from the solution
            let a_resolution = 1e-12f64.max((u + v) * 10e-16);
            if du.abs() + dv.abs() > a_resolution {
                continue; // Precision::PConfusion()  continue;
            }

            // pdn PRO10109 4517: protect against wrong result
            let rs2 = rs.length_squared();
            if rs2 > rsfirst.length_squared() {
                break;
            }

            let rsn = rs.dot(n);
            if rs2 - rsn * rsn / nrm2 > tol2 {
                break;
            }

            // OK, return the result
            *sol = DVec2::new(u, v);

            //: q6
            return if nrm2 < 0.01 * ru2 * rv2 { 2 } else { 1 };
        }
        0
    }

    /// OCCT NextValueOfUV(p2dPrev, P3D, preci, maxpreci) (cxx L1164-1243) —
    /// optimizing projection by Newton algo (S4030).
    pub fn next_value_of_uv(
        &mut self,
        p2d_prev: DVec2,
        p3d: DVec3,
        preci: f64,
        maxpreci: f64,
    ) -> DVec2 {
        let surftype = surf_kind(&self.my_surf);

        match surftype {
            SurfKind::Bezier
            | SurfKind::BSpline
            | SurfKind::Extrusion
            | SurfKind::Revolution
            | SurfKind::Offset => {
                if surftype == SurfKind::BSpline {
                    // Check near to knot position ~ near to C0 points on U
                    // isoline.
                    if geom_adaptor_u_continuity_c0(&self.my_surf)
                        && bspline_has_uknot_near(&self.my_surf, p2d_prev.x)
                    {
                        return self.value_of_uv(p3d, preci);
                    }
                    // Check near to knot position ~ near to C0 points on V
                    // isoline.
                    if geom_adaptor_v_continuity_c0(&self.my_surf)
                        && bspline_has_vknot_near(&self.my_surf, p2d_prev.y)
                    {
                        return self.value_of_uv(p3d, preci);
                    }
                }

                let mut sol = DVec2::ZERO;
                let res = self.surface_newton(p2d_prev, p3d, preci, &mut sol);
                if res != 0 {
                    let gap = p3d.distance(self.value(sol));
                    if res == 2 || (maxpreci > 0.0 && gap - maxpreci > CONFUSION) {
                        //: q6 abv 19 Mar 99 / :q1: check with maxpreci
                        let mut u = sol.x;
                        let mut v = sol.y;
                        self.uv_from_iso(p3d, preci, &mut u, &mut v);
                        // OCCT: if (gap >= myGap) return gp_Pnt2d(U, V).
                        if gap >= self.my_gap {
                            return DVec2::new(u, v);
                        }
                    }
                    self.my_gap = gap;
                    return sol;
                }
            }
            _ => {}
        }
        self.value_of_uv(p3d, preci)
    }

    /// OCCT ValueOfUV(P3D, preci) (cxx L1245-1520).
    pub fn value_of_uv(&mut self, p3d: DVec3, preci: f64) -> DVec2 {
        let mut s = 0.0f64;
        let mut t = 0.0f64;
        self.my_gap = -1.0; // devra etre calcule
        let mut computed = true; // a priori

        let (mut uf, mut ul, mut vf, mut vl) = (0.0, 0.0, 0.0, 0.0);
        self.bounds(&mut uf, &mut ul, &mut vf, &mut vl);

        { //: c9 abv 3 Mar 98: UKI60107-1 #350: to prevent 'catch' from
            // catching exception raising below it
            let surftype = surf_kind(&self.my_surf);
            match surftype {
                SurfKind::Plane
                | SurfKind::Cylinder
                | SurfKind::Cone
                | SurfKind::Sphere
                | SurfKind::Torus => {
                    // ElSLib::Parameters(..., P3D, S, T) — the kernel
                    // elslib_*_parameters re-hosts.
                    match &self.my_surf {
                        Surface3::Plane(pl) => {
                            let (ss, tt) = rcad_kernel::math::el::elslib_plane_parameters(
                                p3d, pl.origin, pl.u_dir, pl.v_dir,
                            );
                            s = ss;
                            t = tt;
                        }
                        Surface3::Cylinder(cyl) => {
                            let x_ax = cyl.ref_dir.normalize_or_zero();
                            let y_ax = cyl.y_axis();
                            let (ss, tt) = rcad_kernel::math::el::elslib_cylinder_parameters(
                                p3d, cyl.origin, x_ax, y_ax, cyl.axis, cyl.radius,
                            );
                            s = ss;
                            t = tt;
                            s += adjust_by_period(s, 0.5 * (uf + ul), 2.0 * std::f64::consts::PI);
                        }
                        Surface3::Cone(cone) => {
                            let x_ax = cone.ref_dir.normalize_or_zero();
                            let y_ax = x_ax.cross(cone.axis);
                            let (ss, tt) = rcad_kernel::math::el::elslib_cone_parameters(
                                p3d,
                                cone.apex,
                                x_ax,
                                y_ax,
                                cone.axis,
                                cone.radius,
                                cone.half_angle_rad,
                            );
                            s = ss;
                            t = tt;
                            s += adjust_by_period(s, 0.5 * (uf + ul), 2.0 * std::f64::consts::PI);
                        }
                        Surface3::Sphere(sph) => {
                            let x_ax = sph.ref_dir.normalize_or_zero();
                            let y_ax = x_ax.cross(sph.axis);
                            let (ss, tt) = rcad_kernel::math::el::elslib_sphere_parameters(
                                p3d, sph.center, x_ax, y_ax, sph.axis,
                            );
                            s = ss;
                            t = tt;
                            s += adjust_by_period(s, 0.5 * (uf + ul), 2.0 * std::f64::consts::PI);
                        }
                        Surface3::Torus(tor) => {
                            let x_ax = tor.ref_dir.normalize_or_zero();
                            let y_ax = x_ax.cross(tor.axis);
                            let (ss, tt) = rcad_kernel::math::el::elslib_torus_parameters(
                                p3d,
                                tor.center,
                                x_ax,
                                y_ax,
                                tor.axis,
                                tor.major_radius,
                                tor.minor_radius,
                            );
                            s = ss;
                            t = tt;
                            s += adjust_by_period(s, 0.5 * (uf + ul), 2.0 * std::f64::consts::PI);
                            t += adjust_by_period(t, 0.5 * (vf + vl), 2.0 * std::f64::consts::PI);
                        }
                        _ => unreachable!(),
                    }
                }
                SurfKind::Bezier
                | SurfKind::BSpline
                | SurfKind::Extrusion
                | SurfKind::Revolution
                | SurfKind::Offset => {
                    //: d0 abv 3 Mar 98: UKI60107-1 #350
                    s = (uf + ul) / 2.0;
                    t = (vf + vl) / 2.0; // yaura aumoins qqchose
                    // pdn to fix hangs PRO17015
                    if surftype == SurfKind::Extrusion
                        && is_infinite_value(uf)
                        && is_infinite_value(ul)
                    {
                        // conic case
                        let prev = DVec2::new(s, t);
                        let mut solution = DVec2::ZERO;
                        if self.surface_newton(prev, p3d, preci, &mut solution) != 0 {
                            return solution;
                        }
                        uf = -500.0;
                        ul = 500.0;
                    }

                    restrict_bounds(&mut uf, &mut ul, &mut vf, &mut vl);

                    //: 30 by abv 2.12.97: speed optimization
                    // code is taken from GeomAPI_ProjectPointOnSurf
                    if !self.my_ext_ok {
                        //  Forcer appel a IsU-VClosed
                        if self.my_uclose_val < 0.0 {
                            self.is_u_closed(preci);
                        }
                        if self.my_vclose_val < 0.0 {
                            self.is_v_closed(preci);
                        }
                        let mut du = 0.0f64;
                        let mut dv = 0.0f64;
                        // extension of the surface range is limited to
                        // non-offset surfaces (see id23943)
                        if !matches!(self.my_surf, Surface3::Offset(_)) {
                            // modified by rln during fixing CSR # BUC60035
                            du = self.my_udelt.min(self.my_surf.u_resolution(preci));
                            dv = self.my_vdelt.min(self.my_surf.v_resolution(preci));
                        }
                        let tol = PCONFUSION;
                        // myExtPS.SetFlag(Extrema_ExtFlag_MIN) — the kernel
                        // minima-only contract (bridge #5).
                        // myExtPS.Initialize(SurfAdapt, uf-du, ul+du,
                        //                    vf-dv, vl+dv, Tol, Tol).
                        self.my_ext_domain = [uf - du, ul + du, vf - dv, vl + dv];
                        self.my_ext_ps = Some(ExtPS::with_domain(
                            p3d,
                            &self.my_surf.clone(),
                            uf - du,
                            ul + du,
                            vf - dv,
                            vl + dv,
                            tol,
                            tol,
                        ));
                        self.my_ext_ok = true;
                    }
                    // myExtPS.Perform(P3D) — the kernel re-construction at
                    // the stored domain (bridge #5).
                    let dom = self.my_ext_domain;
                    let mut ext = ExtPS::with_domain(
                        p3d,
                        &self.my_surf.clone(),
                        dom[0],
                        dom[1],
                        dom[2],
                        dom[3],
                        PCONFUSION,
                        PCONFUSION,
                    );
                    let n_p_surf = if ext.is_done() { ext.nb_ext() } else { 0 };
                    self.my_ext_ps = Some(ext);

                    if n_p_surf > 0 {
                        let ext_ref = self.my_ext_ps.as_ref().unwrap();
                        let mut dist2_min = ext_ref.square_distance(1);
                        let mut ind_min = 1usize;
                        for sol in 2..=n_p_surf {
                            let dist2 = ext_ref.square_distance(sol);
                            if dist2_min > dist2 {
                                dist2_min = dist2;
                                ind_min = sol;
                            }
                        }
                        let p_on_s: POnSurface = ext_ref.point(ind_min).clone();
                        // myExtPS.Point(indMin).Parameter(S, T).
                        s = p_on_s.u;
                        t = p_on_s.v;
                        // PTV 26.06.2002 WORKAROUND protect OCC486.
                        let a_check_pnt = self.my_surf.point_at(s, t);
                        dist2_min = p3d.distance_squared(a_check_pnt);
                        // end of WORKAROUND
                        let mut dis_surf = dist2_min.sqrt();

                        // Test de projection merdeuse sur les bords :
                        let mut uu = s;
                        let mut vv = t;
                        let mut dist_min_on_iso = REAL_LAST;
                        let mut poss_lockal = false; //: study S4030 (optimizing)
                        if dis_surf > preci {
                            let mut pp = DVec2::new(uu, vv);
                            if self.surface_newton(pp, p3d, preci, &mut pp) != 0 {
                                //: q2 abv 16 Mar 99: PRO7226 #412920
                                let dist = p3d.distance(self.value(pp));
                                if dist < dis_surf {
                                    dis_surf = dist;
                                    s = pp.x;
                                    uu = pp.x;
                                    t = pp.y;
                                    vv = pp.y;
                                }
                            }
                            if dis_surf < 10.0 * preci {
                                if !surface_continuity_is_c0(&self.my_surf) {
                                    let tol = CONFUSION;
                                    let (pnt, d1u, d1v) = self.my_surf.derivatives(uu, vv);
                                    let b = d1u.cross(d1v);
                                    let a = p3d - pnt;
                                    let ab = a.dot(b);
                                    let nrm2 = b.length_squared();
                                    if nrm2 > 1e-10 {
                                        let dist = a.length_squared() - (ab * ab) / nrm2;
                                        poss_lockal = dist < tol * tol;
                                    }
                                }
                            }
                            if !poss_lockal {
                                dist_min_on_iso =
                                    self.uv_from_iso(p3d, preci, &mut uu, &mut vv);
                            }
                        }

                        if dis_surf > dist_min_on_iso {
                            // On prend les parametres UU et VV;
                            s = uu;
                            t = vv;
                            self.my_gap = dist_min_on_iso;
                        } else {
                            self.my_gap = dis_surf;
                        }
                    } else {
                        // on essai sur les bords
                        let mut uu = s;
                        let mut vv = t;
                        self.uv_from_iso(p3d, preci, &mut uu, &mut vv);
                        s = uu;
                        t = vv;
                    }
                }
                _ => {
                    computed = false;
                }
            }
        } //: c9 end Try ValueOfUV (CKY 30-DEC-1997)
        // The OCCT catch arm (bridge #8): S/T fall back to the half bounds.
        if computed {
            if self.my_gap <= 0.0 {
                self.my_gap = p3d.distance(self.my_surf.point_at(s, t));
            }
        } else {
            self.my_gap = -1.0;
            s = 0.0;
            t = 0.0;
        }
        DVec2::new(s, t)
    }

    /// OCCT UVFromIso(P3d, preci, U, V) (cxx L1522-1843) — the iso-walk
    /// projection.
    pub fn uv_from_iso(&mut self, p3d: DVec3, preci: f64, u: &mut f64, v: &mut f64) -> f64 {
        //  Projection qui considere les isos ... comme suit :
        //  Les 4 bords, plus les isos en U et en V
        //  En effet, souvent, un des deux est bon ...
        let mut the_min;

        //  Initialisation des recherches : point deja trouve (?)
        let mut uu = *u;
        let mut vv = *v;
        let depart = self.my_surf.point_at(*u, *v);
        the_min = depart.distance(p3d);

        if the_min < preci / 10.0 {
            return the_min; // c etait deja OK
        }
        self.compute_boxes();
        if self.my_iso_uf.is_none()
            || self.my_iso_ul.is_none()
            || self.my_iso_vf.is_none()
            || self.my_iso_vl.is_none()
        {
            // no isolines
            // no more precise computation
            // (bridge #6: the ComputeIso GAP keeps this the live path)
            return the_min;
        }
        // OCC_CATCH_SIGNALS (bridge #8) — the OCCT catch sets
        // theMin = RealLast().
        // pdn Create BndBox containing point;
        let mut a_pbox = BndBox::new();
        a_pbox.add_point(p3d);

        // modified by rln on 04/12/97 in order to use these variables later
        let mut uv = true;
        let mut par;
        let mut other = 0.0f64;
        let mut dist;
        let sac = ShapeAnalysisCurve;
        let is_offset = matches!(self.my_surf, Surface3::Offset(_));
        for num in 0..6i32 {
            uv = num < 3; // 0-1-2 : iso-U  3-4-5 : iso-V
            if !is_offset {
                let (an_iso_box, iso): (Option<&BndBox>, Option<Curve3>) = match num {
                    0 => {
                        par = self.my_uf;
                        (Some(&self.my_bnd_uf), self.my_iso_uf.clone())
                    }
                    1 => {
                        par = self.my_ul;
                        (Some(&self.my_bnd_ul), self.my_iso_ul.clone())
                    }
                    2 => {
                        par = *u;
                        (None, self.uiso(*u))
                    }
                    3 => {
                        par = self.my_vf;
                        (Some(&self.my_bnd_vf), self.my_iso_vf.clone())
                    }
                    4 => {
                        par = self.my_vl;
                        (Some(&self.my_bnd_vl), self.my_iso_vl.clone())
                    }
                    _ => {
                        par = *v;
                        (None, self.viso(*v))
                    }
                };

                //    On y va la-dessus
                if !is_infinite_value(par) && iso.is_some() {
                    if let Some(an_iso_box) = an_iso_box {
                        if bnd_distance_lower(an_iso_box, &a_pbox) > the_min {
                            continue;
                        }
                    }

                    let iso_v = iso.clone().unwrap();
                    let mut cf = CurveEval::default_domain(&iso_v)[0];
                    let mut cl = CurveEval::default_domain(&iso_v)[1];

                    restrict_bounds_1(&mut cf, &mut cl);
                    let mut pntres = DVec3::ZERO;
                    dist = sac.project_cf_cl(
                        &iso_v,
                        p3d,
                        preci,
                        &mut pntres,
                        &mut other,
                        cf,
                        cl,
                        false,
                    );
                    if dist < the_min {
                        the_min = dist;
                        //: q6  Selon une isoU, on calcule le meilleur V;
                        //  et lycee de Versailles
                        uu = if uv { par } else { other };
                        vv = if uv { other } else { par };
                    }
                }
            } else {
                // The OffsetSurface arm (Adaptor3d_IsoCurve, bridge #6 GAP):
                // the iso adaptor load demotes to the ComputeIso handles.
                let (an_iso_box, iso): (Option<&BndBox>, Option<Curve3>) = match num {
                    0 => {
                        par = self.my_uf;
                        (Some(&self.my_bnd_uf), self.my_iso_uf.clone())
                    }
                    1 => {
                        par = self.my_ul;
                        (Some(&self.my_bnd_ul), self.my_iso_ul.clone())
                    }
                    2 => {
                        par = *u;
                        (None, None) // anIsoCurve.Load(GeomAbs_IsoU, U) — GAP
                    }
                    3 => {
                        par = self.my_vf;
                        (Some(&self.my_bnd_vf), self.my_iso_vf.clone())
                    }
                    4 => {
                        par = self.my_vl;
                        (Some(&self.my_bnd_vl), self.my_iso_vl.clone())
                    }
                    _ => {
                        par = *v;
                        (None, None) // anIsoCurve.Load(GeomAbs_IsoV, V) — GAP
                    }
                };
                if let Some(an_iso_box) = an_iso_box {
                    if bnd_distance_lower(an_iso_box, &a_pbox) > the_min {
                        continue;
                    }
                }
                if let Some(iso_v) = iso {
                    let adaptor = GeomIsoAdaptor::new(iso_v);
                    let mut pntres = DVec3::ZERO;
                    let mut param = 0.0f64;
                    dist = sac.project_adaptor(
                        &adaptor,
                        p3d,
                        preci,
                        &mut pntres,
                        &mut param,
                        false,
                    );
                    other = param;
                    if dist < the_min {
                        the_min = dist;
                        uu = if uv { par } else { other };
                        vv = if uv { other } else { par };
                    }
                }
            }
        }

        // added by rln on 04/12/97 iterational process
        let mut prev_u = *u;
        let mut prev_v = *v;
        let max_iters = 5;
        let mut iters = 0;
        if !is_offset {
            while (prev_u != uu || prev_v != vv) && iters < max_iters && the_min > preci {
                prev_u = uu;
                prev_v = vv;
                let iso = if uv { self.uiso(uu) } else { self.viso(vv) };
                if let Some(iso_v) = iso {
                    let mut cf = CurveEval::default_domain(&iso_v)[0];
                    let mut cl = CurveEval::default_domain(&iso_v)[1];
                    restrict_bounds_1(&mut cf, &mut cl);
                    let mut pntres = DVec3::ZERO;
                    dist = sac.project_cf_cl(
                        &iso_v,
                        p3d,
                        preci,
                        &mut pntres,
                        &mut other,
                        cf,
                        cl,
                        false,
                    );
                    if dist < the_min {
                        the_min = dist;
                        if uv {
                            vv = other;
                        } else {
                            uu = other;
                        }
                    }
                }
                uv = !uv;
                let iso = if uv { self.uiso(uu) } else { self.viso(vv) };
                if let Some(iso_v) = iso {
                    let mut cf = CurveEval::default_domain(&iso_v)[0];
                    let mut cl = CurveEval::default_domain(&iso_v)[1];
                    restrict_bounds_1(&mut cf, &mut cl);
                    let mut pntres = DVec3::ZERO;
                    dist = sac.project_cf_cl(
                        &iso_v,
                        p3d,
                        preci,
                        &mut pntres,
                        &mut other,
                        cf,
                        cl,
                        false,
                    );
                    if dist < the_min {
                        the_min = dist;
                        if uv {
                            vv = other;
                        } else {
                            uu = other;
                        }
                    }
                }
                uv = !uv;
                iters += 1;
            }
        } else {
            // The offset arm (bridge #6 GAP): the Adaptor3d_IsoCurve loads
            // demote to the ComputeIso handles (null with the GAP above);
            // the loop keeps the OCCT structure over the available isos.
            while (prev_u != uu || prev_v != vv) && iters < max_iters && the_min > preci {
                prev_u = uu;
                prev_v = vv;
                let mut iso = if uv { self.uiso(uu) } else { self.viso(vv) };
                let mut cf;
                let mut cl;
                if let Some(iso_v) = iso.take() {
                    let adaptor = GeomIsoAdaptor::new(iso_v);
                    cf = super::curve::Adaptor3dCurve::first_parameter(&adaptor);
                    cl = super::curve::Adaptor3dCurve::last_parameter(&adaptor);
                    restrict_bounds_1(&mut cf, &mut cl);
                    let mut pntres = DVec3::ZERO;
                    let mut param = 0.0f64;
                    dist = sac.project_adaptor(
                        &adaptor,
                        p3d,
                        preci,
                        &mut pntres,
                        &mut param,
                        false,
                    );
                    other = param;
                    if dist < the_min {
                        the_min = dist;
                        if uv {
                            vv = other;
                        } else {
                            uu = other;
                        }
                    }
                }
                uv = !uv;
                let iso = if uv { self.uiso(uu) } else { self.viso(vv) };
                if let Some(iso_v) = iso {
                    cf = CurveEval::default_domain(&iso_v)[0];
                    cl = CurveEval::default_domain(&iso_v)[1];
                    restrict_bounds_1(&mut cf, &mut cl);
                    let mut pntres = DVec3::ZERO;
                    let mut param = 0.0f64;
                    // OCCT second pass: ShapeAnalysis_Curve::ProjectAct.
                    dist = sac.project_act(
                        &GeomIsoAdaptor::new(iso_v),
                        p3d,
                        preci,
                        &mut pntres,
                        &mut param,
                    );
                    other = param;
                    if dist < the_min {
                        the_min = dist;
                        if uv {
                            vv = other;
                        } else {
                            uu = other;
                        }
                    }
                }
                uv = !uv;
                iters += 1;
            }
        }

        *u = uu;
        *v = vv;

        the_min
    }

    /// OCCT SortSingularities() (cxx L1845-1885).
    pub fn sort_singularities(&mut self) {
        for i in 0..(self.my_nb_deg - 1).max(0) as usize {
            let mut min_preci = self.my_preci[i];
            let mut min_index = i;
            for j in (i + 1)..self.my_nb_deg.max(0) as usize {
                if min_preci > self.my_preci[j] {
                    min_preci = self.my_preci[j];
                    min_index = j;
                }
            }
            if min_index != i {
                self.my_preci[min_index] = self.my_preci[i];
                self.my_preci[i] = min_preci;
                self.my_p3d.swap(min_index, i);
                self.my_first_p2d.swap(min_index, i);
                self.my_last_p2d.swap(min_index, i);
                self.my_first_par.swap(min_index, i);
                self.my_last_par.swap(min_index, i);
                self.my_uiso_deg.swap(min_index, i);
            }
        }
    }

    /// OCCT SetDomain(U1, U2, V1, V2) (cxx L1887-1896).
    pub fn set_domain(&mut self, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.my_uf = u1;
        self.my_ul = u2;
        self.my_vf = v1;
        self.my_vl = v2;
    }

    /// OCCT ComputeBoxes() (cxx L1898-1917).
    pub fn compute_boxes(&mut self) {
        if self.my_iso_boxes {
            return;
        }
        self.my_iso_boxes = true;
        self.compute_bound_isos();
        // BndLib_Add3dCurve::Add(GeomAdaptor_Curve(iso), Confusion(), bnd)
        // — the base::bnd_lib curve walk (bridge #7).
        if let Some(iso) = self.my_iso_uf.as_ref() {
            self.my_bnd_uf = iso_bnd_box(iso);
        }
        if let Some(iso) = self.my_iso_ul.as_ref() {
            self.my_bnd_ul = iso_bnd_box(iso);
        }
        if let Some(iso) = self.my_iso_vf.as_ref() {
            self.my_bnd_vf = iso_bnd_box(iso);
        }
        if let Some(iso) = self.my_iso_vl.as_ref() {
            self.my_bnd_vl = iso_bnd_box(iso);
        }
    }

    /// OCCT GetBoxUF() (cxx L1921-1924).
    pub fn get_box_uf(&mut self) -> &BndBox {
        self.compute_boxes();
        &self.my_bnd_uf
    }

    /// OCCT GetBoxUL() (cxx L1926-1929).
    pub fn get_box_ul(&mut self) -> &BndBox {
        self.compute_boxes();
        &self.my_bnd_ul
    }

    /// OCCT GetBoxVF() (cxx L1931-1934).
    pub fn get_box_vf(&mut self) -> &BndBox {
        self.compute_boxes();
        &self.my_bnd_vf
    }

    /// OCCT GetBoxVL() (cxx L1936-1939).
    pub fn get_box_vl(&mut self) -> &BndBox {
        self.compute_boxes();
        &self.my_bnd_vl
    }
}

// ---------------------------------------------------------------------------
// Free helpers
// ---------------------------------------------------------------------------

/// OCCT Geom_BSplineSurface::IsURational() — any weight differs from 1.
fn is_u_rational(bs: &rcad_kernel::geom::BSplineSurface) -> bool {
    bs.weights.iter().flatten().any(|w| *w != 1.0)
}

/// OCCT Geom_BSplineSurface::IsVRational() — the rcad weight grid is a
/// single table; the U/V split follows the same values.
fn is_v_rational(bs: &rcad_kernel::geom::BSplineSurface) -> bool {
    is_u_rational(bs)
}

/// OCCT GeomAdaptor_Surface::UContinuity() (the C0 test): true when an
/// interior U knot carries the full degree multiplicity (the flat-knot
/// re-host; GeomAdaptor_Surface.cxx Load BSpline arm).
fn geom_adaptor_u_continuity_c0(s: &Surface3) -> bool {
    let Surface3::BSpline(bs) = s else {
        return false;
    };
    let tab = KnotTable::new(&bs.knots_u);
    for i in 2..tab.len() {
        if tab.mult(i) as usize >= bs.degree_u {
            return true;
        }
    }
    false
}

/// OCCT GeomAdaptor_Surface::VContinuity() (the C0 test).
fn geom_adaptor_v_continuity_c0(s: &Surface3) -> bool {
    let Surface3::BSpline(bs) = s else {
        return false;
    };
    let tab = KnotTable::new(&bs.knots_v);
    for i in 2..tab.len() {
        if tab.mult(i) as usize >= bs.degree_v {
            return true;
        }
    }
    false
}

/// The UKnot within Confusion of the value (the C0-knot check of
/// NextValueOfUV, the FirstUKnotIndex..LastUKnotIndex walk).
fn bspline_has_uknot_near(s: &Surface3, x: f64) -> bool {
    let Surface3::BSpline(bs) = s else {
        return false;
    };
    let tab = KnotTable::new(&bs.knots_u);
    for i in 1..=tab.len() {
        if (tab.knot(i) - x).abs() < CONFUSION {
            return true;
        }
    }
    false
}

/// The VKnot within Confusion of the value.
fn bspline_has_vknot_near(s: &Surface3, y: f64) -> bool {
    let Surface3::BSpline(bs) = s else {
        return false;
    };
    let tab = KnotTable::new(&bs.knots_v);
    for i in 1..=tab.len() {
        if (tab.knot(i) - y).abs() < CONFUSION {
            return true;
        }
    }
    false
}

/// OCCT mySurf->Continuity() != GeomAbs_C0 (the possLockal test) — the
/// canonical surfaces are CN; a BSpline is C0 exactly when the U/V C0
/// condition holds.
fn surface_continuity_is_c0(s: &Surface3) -> bool {
    match surf_kind(s) {
        SurfKind::BSpline => {
            geom_adaptor_u_continuity_c0(s) || geom_adaptor_v_continuity_c0(s)
        }
        SurfKind::Bezier
        | SurfKind::Plane
        | SurfKind::Cylinder
        | SurfKind::Cone
        | SurfKind::Sphere
        | SurfKind::Torus => false,
        _ => false,
    }
}

/// OCCT ComputeIso(surf, utype, par) (cxx L579-603) — the try/catch around
/// surf->UIso/VIso (bridge #8: the catch keeps the null handle).
/// GAP (bridge #6): the per-surface iso construction (TKG3d
/// Geom_Surface::UIso/VIso) has no kernel port — the null handle path is
/// kept.
fn compute_iso(_surf: &Surface3, _utype: bool, _par: f64) -> Option<Curve3> {
    None
}

/// OCCT Bnd_Box::Distance(other) — the lower distance between two boxes (0
/// when they touch or overlap).
fn bnd_distance_lower(a: &BndBox, b: &BndBox) -> f64 {
    let (Some((axmin, aymin, azmin, axmax, aymax, azmax)), Some((bxmin, bymin, bzmin, bxmax, bymax, bzmax))) =
        (a.get(), b.get())
    else {
        return 0.0;
    };
    let dx = (axmin - bxmax).max(bxmin - axmax).max(0.0);
    let dy = (aymin - bymax).max(bymin - aymax).max(0.0);
    let dz = (azmin - bzmax).max(bzmin - azmax).max(0.0);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// OCCT BndLib_Add3dCurve::Add(GeomAdaptor_Curve(iso), Confusion(), bnd) —
/// the base::bnd_lib curve walk feeding the kernel BndBox (bridge #7).
fn iso_bnd_box(c: &Curve3) -> BndBox {
    let mut bnd = BndBox::new();
    if let Some([p1, p2]) = rcad_kernel::base::bnd_lib::curve_bounding_box(c) {
        bnd.add_point(p1);
        bnd.add_point(p2);
    }
    bnd
}

/// The Adaptor3d_Curve over a concrete iso curve (bridge #6: the
/// GeomAdaptor_Curve::Load form of the offset arm).
struct GeomIsoAdaptor {
    c: Curve3,
}

impl GeomIsoAdaptor {
    fn new(c: Curve3) -> Self {
        GeomIsoAdaptor { c }
    }
}

impl CurveEval for GeomIsoAdaptor {
    fn point_at(&self, t: f64) -> DVec3 {
        CurveEval::point_at(&self.c, t)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        CurveEval::tangent_at(&self.c, t)
    }
    fn default_domain(&self) -> [f64; 2] {
        CurveEval::default_domain(&self.c)
    }
}

impl super::curve::Adaptor3dCurve for GeomIsoAdaptor {
    fn first_parameter(&self) -> f64 {
        CurveEval::default_domain(&self.c)[0]
    }
    fn last_parameter(&self) -> f64 {
        CurveEval::default_domain(&self.c)[1]
    }
    fn value(&self, the_u: f64) -> DVec3 {
        CurveEval::point_at(&self.c, the_u)
    }
    fn resolution(&self, the_r3d: f64) -> f64 {
        self.c.resolution(the_r3d)
    }
    fn is_closed(&self) -> bool {
        CurveEval::is_closed(&self.c)
    }
    fn is_kind_bounded(&self) -> bool {
        matches!(
            self.c,
            Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
        )
    }
    fn get_type(&self) -> super::curve::AdaptorCurveKind {
        match self.c {
            Curve3::Line(_) => super::curve::AdaptorCurveKind::Line,
            Curve3::Circle(_) => super::curve::AdaptorCurveKind::Circle,
            Curve3::Ellipse(_) => super::curve::AdaptorCurveKind::Ellipse,
            Curve3::Hyperbola(_) => super::curve::AdaptorCurveKind::Hyperbola,
            Curve3::Parabola(_) => super::curve::AdaptorCurveKind::Parabola,
            Curve3::Bezier(_) => super::curve::AdaptorCurveKind::Bezier,
            Curve3::BSpline(_) => super::curve::AdaptorCurveKind::BSpline,
            Curve3::Offset(_) => super::curve::AdaptorCurveKind::Offset,
            _ => super::curve::AdaptorCurveKind::Other,
        }
    }
    fn curve3d(&self) -> Option<Curve3> {
        Some(self.c.clone())
    }
}
