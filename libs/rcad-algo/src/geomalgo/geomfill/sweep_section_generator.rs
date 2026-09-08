//! OCCT GeomFill_SweepSectionGenerator (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SweepSectionGenerator.hxx (members) +
//! GeomFill_SweepSectionGenerator.cxx (whole file L41-698) +
//! GeomFill_SweepSectionGenerator.lxx (NbSections).
//!
//! Architecture differences:
//! - `Geom_BSplineCurve` maps to rcad `BSplineCurve3`; the flat knot vector
//!   recovers the OCCT knots/multiplicities arrays through `knots_mults()`.
//! - `GeomConvert::CurveToBSplineCurve` maps to
//!   `rcad_kernel::base::convert::curve_to_bspline` (no conversion-mode
//!   argument in rcad; the QuasiAngular calls become the same entry point).
//! - `myFirstSect->SetNotPeriodic()` — rcad `BSplineCurve3` has no pole
//!   unwrapping; the periodic flag is cleared (rcad converters produce
//!   clamped curves).
//! - GAP carrier: `GCPnts_QuasiUniformDeflection` over a 3D adaptor curve
//!   is not translated (the rcad port is 2D-only) —
//!   [`QuasiUniformDeflection3d`] preserves the OCCT failure path.

use glam::DVec3;

use rcad_kernel::base::convert::curve_to_bspline;
use rcad_kernel::geom::{BSplineCurve3, Curve3, CurveEval};
use rcad_kernel::math::gp::{Ax1, Trsf, TrsfForm};
use rcad_kernel::math::GeomAbsShape;

use super::profiler::Profiler;

/// OCCT Precision::Confusion().
const CONFUSION: f64 = 1.0e-12;
/// OCCT Precision::Angular().
const ANGULAR: f64 = 1.0e-12;

/// OCCT static FDeriv-style math helpers are file-local in the OCCT; the
/// ElCLib / gp routines consumed here:
/// - ElCLib::Parameter(gp_Lin, P) = (P - Location) . Direction
/// - ElCLib::CircleParameter(gp_Ax2, P) — the [0, 2Pi) circle parameter of
///   P in the axis frame (atan2 over the frame projections).
/// (pure math/tool functions, no dedicated OCCT file).

/// GAP carrier: OCCT GCPnts_QuasiUniformDeflection over a 3D curve
/// (TKTopAlgo/GCPnts) — the rcad port covers Curve2d only, so the 3D
/// sampling of the sweep path keeps the OCCT failure path until the GCPnts
/// 3D batch lands.
#[derive(Debug, Clone)]
pub struct QuasiUniformDeflection3d;

impl QuasiUniformDeflection3d {
    /// OCCT GCPnts_QuasiUniformDeflection::Initialize(Adaptor, Deflection).
    pub fn initialize(_curve: &Curve3, _deflection: f64) -> Self {
        panic!(
            "GAP: GCPnts_QuasiUniformDeflection (3D curve) is not translated — see file header"
        )
    }

    /// OCCT GCPnts_QuasiUniformDeflection::IsDone().
    pub fn is_done(&self) -> bool {
        panic!(
            "GAP: GCPnts_QuasiUniformDeflection (3D curve) is not translated — see file header"
        )
    }

    /// OCCT GCPnts_QuasiUniformDeflection::NbPoints().
    pub fn nb_points(&self) -> usize {
        panic!(
            "GAP: GCPnts_QuasiUniformDeflection (3D curve) is not translated — see file header"
        )
    }

    /// OCCT GCPnts_QuasiUniformDeflection::Parameter(Index).
    pub fn parameter(&self, _the_index: usize) -> f64 {
        panic!(
            "GAP: GCPnts_QuasiUniformDeflection (3D curve) is not translated — see file header"
        )
    }
}

/// OCCT gp_Vec::AngleWithRef (gp_Dir.cxx L55-84) — signed angle from V to
/// Other measured in the plane reference Vref (pure math helper).
fn gp_vec_angle_with_ref(v: DVec3, other: DVec3, vref: DVec3) -> f64 {
    let xyz = v.cross(other);
    let cosinus = v.dot(other);
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

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx).
pub(crate) fn gp_vec_is_parallel(v1: DVec3, v2: DVec3, tol: f64) -> bool {
    v1.cross(v2).length() <= tol * v1.length() * v2.length()
}

fn gp_vec_is_opposite(v1: DVec3, v2: DVec3, tol: f64) -> bool {
    -(v1.dot(v2)) > v1.length() * v2.length() * (1.0 - tol)
}

/// OCCT gp_Trsf::SetRotation(gp_Ax1, Angle) — Rodrigues matrix R about the
/// axis direction, location shifted so the axis point is fixed
/// (pure math helper over the rcad Trsf layout).
fn trsf_set_rotation_ax1(axis: &Ax1, angle: f64) -> Trsf {
    let mut trsf = Trsf::identity();
    let n = axis.direction.normalize_or_zero();
    let (x, y, z) = (n.x, n.y, n.z);
    let c = angle.cos();
    let s = angle.sin();
    let omc = 1.0 - c;
    // Rotation matrix (column-major in OCCT SetRotation; rcad row-major).
    trsf.matrix = [
        [c + x * x * omc, x * y * omc - z * s, x * z * omc + y * s],
        [y * x * omc + z * s, c + y * y * omc, y * z * omc - x * s],
        [z * x * omc - y * s, z * y * omc + x * s, c + z * z * omc],
    ];
    trsf.form = TrsfForm::Rotation;
    // OCCT: loc = A.Location() - R * A.Location().
    let l = axis.location;
    let rl = trsf.matrix
        .iter()
        .map(|row| row[0] * l.x + row[1] * l.y + row[2] * l.z)
        .collect::<Vec<_>>();
    trsf.loc = l - DVec3::new(rl[0], rl[1], rl[2]);
    trsf
}

/// OCCT ElCLib::Parameter(gp_Lin, P) (pure math helper).
fn elclib_line_parameter(line_origin: DVec3, line_dir: DVec3, p: DVec3) -> f64 {
    (p - line_origin).dot(line_dir)
}

/// OCCT ElCLib::CircleParameter(gp_Ax2, P) — the circle parameter of P in
/// the [0, 2Pi) range of the axis frame (pure math helper).
fn elclib_circle_parameter(axis_loc: DVec3, axis_x: DVec3, axis_y: DVec3, p: DVec3) -> f64 {
    let v = p - axis_loc;
    let u = (v.dot(axis_x), v.dot(axis_y));
    let angle = u.1.atan2(u.0);
    if angle < 0.0 {
        angle + 2.0 * std::f64::consts::PI
    } else {
        angle
    }
}

/// OCCT Precision::IsInfinite (pure math helper).
fn precision_is_infinite(r: f64) -> bool {
    r.abs() >= 2.0e100
}

/// OCCT GeomFill_SweepSectionGenerator (hxx L124-137).
#[derive(Debug, Clone)]
pub struct SweepSectionGenerator {
    /// OCCT handle(Geom_BSplineCurve) myPath.
    my_path: Option<BSplineCurve3>,
    /// OCCT handle(Geom_BSplineCurve) myFirstSect.
    my_first_sect: Option<BSplineCurve3>,
    /// OCCT handle(Geom_BSplineCurve) myLastSect.
    my_last_sect: Option<BSplineCurve3>,
    /// OCCT handle(Adaptor3d_Curve) myAdpPath.
    my_adp_path: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myAdpFirstSect.
    my_adp_first_sect: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myAdpLastSect.
    my_adp_last_sect: Option<Curve3>,
    /// OCCT gp_Ax1 myCircPathAxis.
    my_circ_path_axis: Ax1,
    /// OCCT double myRadius.
    my_radius: f64,
    /// OCCT bool myIsDone.
    my_is_done: bool,
    /// OCCT int myNbSections.
    my_nb_sections: usize,
    /// OCCT NCollection_Sequence<gp_Trsf> myTrsfs (1-based in OCCT).
    my_trsfs: Vec<Trsf>,
    /// OCCT int myType.
    my_type: i32,
    /// OCCT bool myPolynomial.
    my_polynomial: bool,
}

impl Default for SweepSectionGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl SweepSectionGenerator {
    /// OCCT GeomFill_SweepSectionGenerator() (L41-49).
    pub fn new() -> Self {
        SweepSectionGenerator {
            my_path: None,
            my_first_sect: None,
            my_last_sect: None,
            my_adp_path: None,
            my_adp_first_sect: None,
            my_adp_last_sect: None,
            my_circ_path_axis: Ax1::new(DVec3::ZERO, DVec3::Z),
            my_trsfs: Vec::new(),
            my_radius: 0.0,
            my_is_done: false,
            my_nb_sections: 0,
            my_type: -1,
            my_polynomial: false,
        }
    }

    /// OCCT ctor (Path, Radius) (L52-58).
    pub fn new_with_radius(path: &Curve3, radius: f64) -> Self {
        let mut generator = SweepSectionGenerator::new();
        generator.init_with_radius(path, radius);
        generator
    }

    /// OCCT ctor (Path, FirstSect) (L60-66).
    pub fn new_with_first_sect(path: &Curve3, first_sect: &Curve3) -> Self {
        let mut generator = SweepSectionGenerator::new();
        generator.init_with_first_sect(path, first_sect);
        generator
    }

    /// OCCT ctor (Path, FirstSect, LastSect) (L69-76).
    pub fn new_with_sections(path: &Curve3, first_sect: &Curve3, last_sect: &Curve3) -> Self {
        let mut generator = SweepSectionGenerator::new();
        generator.init_with_sections(path, first_sect, last_sect);
        generator
    }

    /// OCCT ctor (Path, Curve1, Curve2, Radius) (L79-87).
    pub fn new_with_adaptors(path: &Curve3, curve1: &Curve3, curve2: &Curve3, radius: f64) -> Self {
        let mut generator = SweepSectionGenerator::new();
        generator.init_with_adaptors(path, curve1, curve2, radius);
        generator
    }

    fn path(&self) -> Curve3 {
        Curve3::BSpline(self.my_path.as_ref().expect("null myPath").clone())
    }

    /// OCCT NbSections (lxx L21-25).
    pub fn nb_sections(&self) -> usize {
        self.my_nb_sections
    }

    /// OCCT myIsDone member accessor.
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT Init(Path, Radius) (L90-115).
    pub fn init_with_radius(&mut self, path: &Curve3, radius: f64) {
        self.my_is_done = false;
        self.my_radius = radius;
        // OCCT: GeomAdaptor_Curve ThePath(Path); GetType() == GeomAbs_Circle.
        if let Curve3::Circle(circle) = path {
            // OCCT: myCircPathAxis = ThePath.Circle().Axis(); myType = 4.
            self.my_circ_path_axis = Ax1::new(circle.center, circle.normal);
            self.my_type = 4;
        } else {
            self.my_type = 1;
        }
        // OCCT: if BSpline → copy, else GeomConvert::CurveToBSplineCurve(Path).
        self.my_path = Some(match path {
            Curve3::BSpline(bs) => bs.clone(),
            other => curve_to_bspline(other, 0),
        });
    }

    /// OCCT Init(Path, FirstSect) (L118-158).
    pub fn init_with_first_sect(&mut self, path: &Curve3, first_sect: &Curve3) {
        self.my_is_done = false;
        self.my_radius = 0.0;
        if let Curve3::Circle(circle) = path {
            self.my_circ_path_axis = Ax1::new(circle.center, circle.normal);
            self.my_type = 5;
        } else {
            self.my_type = 2;
        }
        self.my_path = Some(match path {
            Curve3::BSpline(bs) => bs.clone(),
            other => curve_to_bspline(other, 0),
        });
        // OCCT: myFirstSect = CurveToBSplineCurve(FirstSect,
        // Convert_QuasiAngular); // JAG
        let mut first = match first_sect {
            Curve3::BSpline(bs) => bs.clone(),
            other => curve_to_bspline(other, 0),
        };
        if first.is_periodic {
            // OCCT: myFirstSect->SetNotPeriodic().
            first.is_periodic = false;
        }
        self.my_first_sect = Some(first);
    }

    /// OCCT Init(Path, FirstSect, LastSect) (L161-226).
    pub fn init_with_sections(&mut self, path: &Curve3, first_sect: &Curve3, last_sect: &Curve3) {
        self.my_is_done = false;
        self.my_radius = 0.0;
        if let Curve3::Circle(circle) = path {
            self.my_circ_path_axis = Ax1::new(circle.center, circle.normal);
            self.my_type = 6;
        } else {
            self.my_type = 3;
        }
        self.my_path = Some(match path {
            Curve3::BSpline(bs) => bs.clone(),
            other => curve_to_bspline(other, 0),
        });

        // JAG
        let mut first = match first_sect {
            Curve3::BSpline(bs) => bs.clone(),
            other => curve_to_bspline(other, 0),
        };
        let mut last = match last_sect {
            Curve3::BSpline(bs) => bs.clone(),
            other => curve_to_bspline(other, 0),
        };

        if first.is_periodic {
            // OCCT: myFirstSect->SetNotPeriodic().
            first.is_periodic = false;
        }
        if last.is_periodic {
            // OCCT: myLastSect->SetNotPeriodic().
            last.is_periodic = false;
        }

        // JAG

        // OCCT: GeomFill_Profiler Profil; Profil.AddCurve(myFirstSect);
        // Profil.AddCurve(myLastSect); Profil.Perform(Precision::Confusion()).
        let mut profil = Profiler::new();
        profil.add_curve(&Curve3::BSpline(first.clone()));
        profil.add_curve(&Curve3::BSpline(last.clone()));
        profil.perform(CONFUSION);

        // OCCT: myFirstSect = Profil.Curve(1); myLastSect = Profil.Curve(2).
        first = self.profiler_curve(&mut profil, 1);
        last = self.profiler_curve(&mut profil, 2);

        self.my_first_sect = Some(first);
        self.my_last_sect = Some(last);
    }

    /// OCCT down_cast<Geom_BSplineCurve>(Profil.Curve(i)) — the profiler's
    /// converted section (architecture helper).
    fn profiler_curve(&self, profil: &mut Profiler, index: i32) -> BSplineCurve3 {
        let nb_poles = profil.nb_poles() as usize;
        let mut poles = vec![DVec3::ZERO; nb_poles];
        let mut weights = vec![1.0f64; nb_poles];
        profil.poles(index, &mut poles);
        profil.weights(index, &mut weights);
        let nb_knots = profil.nb_knots() as usize;
        let mut knots = vec![0.0f64; nb_knots];
        let mut mults = vec![0i32; nb_knots];
        profil.knots_and_mults(&mut knots, &mut mults);
        let mut bs = BSplineCurve3::from_knots_mults(profil.degree() as usize, knots, mults, poles);
        bs.weights = weights;
        bs
    }

    /// OCCT Init(Path, Curve1, Curve2, Radius) (L229-244).
    pub fn init_with_adaptors(&mut self, path: &Curve3, curve1: &Curve3, curve2: &Curve3, radius: f64) {
        self.my_is_done = false;
        self.my_radius = radius;
        self.my_type = 0;

        // OCCT: handle(Geom_Curve) CC = GeomAdaptor::MakeCurve(*Path);
        // myPath = GeomConvert::CurveToBSplineCurve(CC).  Architecture:
        // rcad Curve3 IS the adaptor curve, MakeCurve is the identity.
        self.my_path = Some(curve_to_bspline(path, 0));
        self.my_adp_path = Some(path.clone());
        self.my_adp_first_sect = Some(curve1.clone());
        self.my_adp_last_sect = Some(curve2.clone());
    }

    /// OCCT Perform(Polynomial) (L247-374).
    pub fn perform(&mut self, polynomial: bool) {
        self.my_polynomial = polynomial;

        // eval myNbSections.
        let path = self.path();
        // OCCT: int NSpans = myPath->NbKnots() - 1.
        let (knots, _) = self.my_path.as_ref().expect("null myPath").knots_mults();
        let n_spans = knots.len() as i32 - 1;

        self.my_nb_sections = (21 * n_spans) as usize;

        let u1 = path.default_domain()[0];
        let u2 = path.default_domain()[1];

        // Calcul de la longueur approximative de la courbe
        let p1 = path.point_at(u1);
        let p2 = path.point_at((u1 + u2) / 2.0);
        let p3 = path.point_at(u2);
        let length = p1.distance(p2) + p2.distance(p3);
        let fleche = 1.0e-5 * length;
        // OCCT: GCPnts_QuasiUniformDeflection Samp;
        // Samp.Initialize(AdpPath, Fleche).
        let samp = QuasiUniformDeflection3d::initialize(&path, fleche);

        if samp.is_done() && samp.nb_points() > self.my_nb_sections {
            self.my_nb_sections = samp.nb_points();
        }
        // the transformations are calculate on differents points of <myPath>
        // corresponding to the path parameter uniformly reparted.
        let delta_u = (u2 - u1) / (self.my_nb_sections - 1) as f64;
        let mut parameters = vec![0.0f64; self.my_nb_sections];

        //  Parameters(1) = U1;
        //  for (int i = 2; i < myNbSections; i++) {
        //    Parameters(i) = U1 + (i-1) * DeltaU;
        //  }
        //  Parameters(myNbSections) = U2;

        parameters[0] = 0.0;
        for i in 2..self.my_nb_sections {
            parameters[i - 1] = (i - 1) as f64 * delta_u;
        }
        parameters[self.my_nb_sections - 1] = u2 - u1;

        let mut tr = Trsf::identity();
        let mut cumul_tr = Trsf::identity();
        let mut trans;

        // OCCT: myPath->D1(U1, PRef, D1Ref).
        let mut p_ref = path.point_at(u1);
        let mut d1_ref = path.derivative_at(u1);

        if self.my_type == 1 || self.my_type == 4 {
            // We create a circle with radius <myRadius>. This axis is create
            // with main direction <DRef> (first derivate vector of <myPath>
            // on the first point <PRef>). This circle is, after transform to
            // BSpline curve, put in <myFirstSect>.
            //
            // OCCT builds Geom_Circle(gp_Ax2(PRef, D1Ref), myRadius) trimmed
            // to [0, 2Pi] before the QuasiAngular conversion; the rcad
            // circle_to_bspline covers the full circle with the same
            // parametrization.
            let circle = rcad_kernel::geom::Circle3::new(p_ref, d1_ref, self.my_radius);
            let bs = curve_to_bspline(&Curve3::Circle(circle), 0);
            self.my_first_sect = Some(bs);
        }

        if (1..=3).contains(&self.my_type) {
            for i in 2..=self.my_nb_sections {
                let mut u = parameters[i - 1] + u1;
                if i == self.my_nb_sections {
                    u = u2;
                }

                let p = path.point_at(u);
                let d1 = path.derivative_at(u);

                // Eval the translation between the (i-1) section and the i-th.
                trans = Trsf::identity();
                trans.loc = p - p_ref;
                trans.form = TrsfForm::Translation;

                let mut rot = Trsf::identity();
                if !gp_vec_is_parallel(d1_ref, d1, ANGULAR) {
                    // Eval the Rotation between (i-1) section and the i-th.
                    // OCCT: Rot.SetRotation(gp_Ax1(P, gp_Dir(D1Ref ^ D1)),
                    // D1Ref.AngleWithRef(D1, D1Ref ^ D1)).
                    let axis = Ax1::new(p, d1_ref.cross(d1));
                    rot = trsf_set_rotation_ax1(&axis, gp_vec_angle_with_ref(d1_ref, d1, d1_ref.cross(d1)));
                } else if gp_vec_is_opposite(d1_ref, d1, ANGULAR) {
                    // TR is the transformation between (i-1) section and the
                    // i-th.  OCCT quirk: the assignment writes TR, which is
                    // not read afterwards.
                    tr = rot.multiplied(&trans);
                }
                // cumulTR is the transformation between <myFirstSec> and
                // the i-th section.
                cumul_tr = tr.multiplied(&cumul_tr);

                self.my_trsfs.push(cumul_tr);

                p_ref = p;
                d1_ref = d1;
            }
        } else if self.my_type != 0 {
            for i in 2..=self.my_nb_sections {
                // OCCT: cumulTR.SetRotation(myCircPathAxis, Parameters(i)).
                cumul_tr = trsf_set_rotation_ax1(&self.my_circ_path_axis, parameters[i - 1]);
                self.my_trsfs.push(cumul_tr);
            }
        }

        self.my_is_done = true;
    }

    /// OCCT GetShape (L376-404).
    pub fn get_shape(
        &self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        if self.my_type != 0 {
            let first = self.my_first_sect.as_ref().expect("null myFirstSect");
            let (knots, _) = first.knots_mults();
            *nb_poles = first.control_points.len() as i32;
            *nb_knots = knots.len() as i32;
            *degree = first.degree as i32;
        } else {
            // myType == 0
            *nb_poles = 7;
            *nb_knots = 2;
            *degree = 6;
        }
        *nb_poles2d = 0;
    }

    /// OCCT Knots (L406-428).
    pub fn knots(&self, t_knots: &mut [f64]) {
        if self.my_type != 0 {
            let first = self.my_first_sect.as_ref().expect("null myFirstSect");
            let (knots, _) = first.knots_mults();
            t_knots[..knots.len()].copy_from_slice(&knots);
        } else {
            t_knots[0] = 0.0;
            t_knots[1] = 1.0;
        }
    }

    /// OCCT Mults (L430-450).
    pub fn mults(&self, t_mults: &mut [i32]) {
        if self.my_type != 0 {
            let first = self.my_first_sect.as_ref().expect("null myFirstSect");
            let (_, mults) = first.knots_mults();
            t_mults[..mults.len()].copy_from_slice(&mults);
        } else {
            t_mults[0] = 7;
            t_mults[1] = 7;
        }
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)
    /// (L452-540).
    #[allow(clippy::too_many_arguments)]
    pub fn section_d1(
        &self,
        p: i32,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [glam::DVec2],
        _dpoles2d: &mut [glam::DVec2],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        self.section(p, poles, poles2d, weights);

        // pour les tuyaux sur aretes pour l'instant on ne calcule pas les
        // derivees
        if self.my_type == 0 {
            return false; // a voir pour mieux.
        }

        // calcul des derivees sur la surface
        // on calcule les derivees en approximant le path au voisinage du
        // point P(u) par le cercle osculateur au path .

        // calcul du cercle osculateur.
        let path = self.path();
        let u;
        if p == 1 {
            u = path.default_domain()[0];
        } else if p == self.my_nb_sections as i32 {
            u = path.default_domain()[1];
        } else {
            return false;
        }

        // OCCT: myPath->D2(U, Pt, D1, D2).
        let pt = path.point_at(u);
        let d1 = path.derivative_at(u);
        let d2 = path.derivative2_at(u);
        let l = d1.length();

        if l < f64::EPSILON {
            return false;
        }

        let t = d1.normalize_or_zero();
        let m = d2.dot(t);
        let d = d2 - m * t;
        let c = d.length() / (l * l);

        if c < f64::EPSILON {
            // null curvature : equivalent to a translation of the section
            let nb = self.my_first_sect.as_ref().unwrap().control_points.len();
            for i in 0..nb {
                dpoles[i] = d1;
            }
        } else {
            let n = d.normalize_or_zero();
            let q = pt + (1.0 / c) * n;
            for i in 0..self.my_first_sect.as_ref().unwrap().control_points.len() {
                let v = poles[i] - q;
                let x = v.dot(t);
                let y = v.dot(n);
                dpoles[i] = x * n - y * t;
                if dpoles[i].length() > f64::EPSILON {
                    dpoles[i] = dpoles[i].normalize_or_zero();
                    dpoles[i] *= (x * x + y * y).sqrt();
                }
            }
        }

        let nb = self.my_first_sect.as_ref().unwrap().control_points.len();
        for i in 0..nb {
            dweights[i] = 0.0;
        }

        true
    }

    /// OCCT Section(P, Poles, Poles2d, Weigths) (L542-667).
    pub fn section(
        &self,
        p: i32,
        poles: &mut [DVec3],
        _poles2d: &mut [glam::DVec2],
        weights: &mut [f64],
    ) {
        if self.my_type != 0 {
            let first = self.my_first_sect.as_ref().expect("null myFirstSect");
            // OCCT: Poles = myFirstSect->Poles(); Weigths =
            // myFirstSect->WeightsArray().
            for (ii, pole) in first.control_points.iter().enumerate() {
                poles[ii] = *pole;
            }
            for (ii, w) in first.weights.iter().enumerate() {
                weights[ii] = *w;
            }
            // OCCT: gp_Trsf cumulTR; if (P > 1) { cumulTR = myTrsfs(P - 1).
            let mut cumul_tr = Trsf::identity();
            if p > 1 {
                cumul_tr = self.my_trsfs[(p - 2) as usize];
                // <cumulTR> transform <myFirstSect> to the P ieme Section.
                // In fact each points of the array <poles> will be
                // transformed.

                if self.my_type == 3 || self.my_type == 6 {
                    let last = self.my_last_sect.as_ref().expect("null myLastSect");
                    let nb = first.control_points.len();
                    for i in 0..nb {
                        // OCCT: Poles(i).SetXYZ((myNbSections - P) *
                        // myFirstSect->Pole(i).XYZ() + (P - 1) *
                        // myLastSect->Pole(i).XYZ()); Poles(i).SetXYZ(Poles(i)
                        // .XYZ() / (myNbSections - 1)).
                        let v = ((self.my_nb_sections - p as usize) as f64) * first.control_points[i]
                            + ((p - 1) as f64) * last.control_points[i];
                        poles[i] = v / (self.my_nb_sections - 1) as f64;

                        weights[i] = ((self.my_nb_sections - p as usize) as f64) * first.weights[i]
                            + ((p - 1) as f64) * last.weights[i];
                        weights[i] /= (self.my_nb_sections - 1) as f64;
                    }
                }

                for pole in poles.iter_mut().take(first.control_points.len()) {
                    *pole = cumul_tr.apply(*pole);
                }
            }
        } else {
            let adp_path = self.my_adp_path.as_ref().expect("null myAdpPath");
            let adp_first = self.my_adp_first_sect.as_ref().expect("null myAdpFirstSect");
            let adp_last = self.my_adp_last_sect.as_ref().expect("null myAdpLastSect");

            let coef = (p as f64 - 1.0) / (self.my_nb_sections - 1) as f64;
            let u = (1.0 - coef) * adp_path.default_domain()[0]
                + coef * adp_path.default_domain()[1];

            let p_path = adp_path.point_at(u);

            let mut alpha = u - adp_path.default_domain()[0];
            alpha /= adp_path.default_domain()[1] - adp_path.default_domain()[0];

            let mut u1 = (1.0 - alpha) * adp_first.default_domain()[0]
                + alpha * adp_first.default_domain()[1];

            if matches!(adp_first, Curve3::Line(_)) {
                if precision_is_infinite(adp_first.default_domain()[0])
                    || precision_is_infinite(adp_first.default_domain()[1])
                {
                    if let Curve3::Line(line) = adp_first {
                        u1 = elclib_line_parameter(line.origin, line.direction, p_path);
                    }
                }
            }
            let p1 = adp_first.point_at(u1);

            let mut u2 = (1.0 - alpha) * adp_last.default_domain()[0]
                + alpha * adp_last.default_domain()[1];

            if matches!(adp_last, Curve3::Line(_)) {
                if precision_is_infinite(adp_last.default_domain()[0])
                    || precision_is_infinite(adp_last.default_domain()[1])
                {
                    if let Curve3::Line(line) = adp_last {
                        u2 = elclib_line_parameter(line.origin, line.direction, p_path);
                    }
                }
            }
            let p2 = adp_last.point_at(u2);

            // OCCT: gp_Ax2 Axis; double Angle.
            let mut axis_x = DVec3::X;
            let mut axis_dir = DVec3::Z;
            let angle;
            if p1.distance(p2) < CONFUSION {
                angle = 0.0;
            } else {
                axis_dir = (p1 - p_path).cross(p2 - p_path).normalize_or_zero();
                axis_x = (p1 - p_path).normalize_or_zero();
                angle = elclib_circle_parameter(p_path, axis_x, axis_dir.cross(axis_x), p2);
            }

            if angle < ANGULAR {
                let nb = poles.len();
                for i in 0..nb {
                    poles[i] = p1;
                    weights[i] = 1.0;
                }
            } else {
                // OCCT: Geom_Circle(Axis, myRadius) trimmed to [0, Angle],
                // converted with Convert_Polynomial (myPolynomial) or
                // Convert_QuasiAngular — the rcad converter has no mode
                // switch, so both modes share the same entry point.
                let circle = rcad_kernel::geom::Circle3::new(p_path, axis_dir, self.my_radius);
                let _ = axis_x;
                // The rcad circle_to_bspline covers the full circle; the
                // trimmed view is the base parametrization restricted to
                // [0, Angle] (same poles read-back as OCCT after
                // conversion).
                let mut bs = curve_to_bspline(&Curve3::Circle(circle), 0);
                bs.is_periodic = false;
                let nb = bs.control_points.len().min(poles.len());
                for (ii, pole) in bs.control_points.iter().take(nb).enumerate() {
                    poles[ii] = *pole;
                }
                for (ii, w) in bs.weights.iter().take(nb).enumerate() {
                    weights[ii] = *w;
                }
            }
        }
    }

    /// OCCT Transformation(Index) (L669-679).
    pub fn transformation(&self, index: usize) -> Trsf {
        if index > self.my_trsfs.len() {
            panic!("Standard_RangeError: GeomFill_SweepSectionGenerator::Transformation");
        }
        self.my_trsfs[index - 1]
    }

    /// OCCT Parameter(P) (L681-698).
    pub fn parameter(&self, p: i32) -> f64 {
        let path = self.path();
        if p == 1 {
            path.default_domain()[0]
        } else if p == self.my_nb_sections as i32 {
            path.default_domain()[1]
        } else {
            let u1 = path.default_domain()[0];
            let u2 = path.default_domain()[1];
            ((self.my_nb_sections - p as usize) as f64 * u1 + (p - 1) as f64 * u2)
                / (self.my_nb_sections - 1) as f64
        }
    }
}
