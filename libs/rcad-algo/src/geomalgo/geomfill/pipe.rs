//! OCCT GeomFill_Pipe (TKGeomAlgo/GeomFill) — 1:1 port of GeomFill_Pipe.hxx
//! (members) + GeomFill_Pipe.cxx (whole file L79-1102) + GeomFill_Pipe.lxx
//! (the inline accessors).
//!
//! Architecture differences:
//! - `handle(Adaptor3d_Curve) myAdpPath` maps to `Option<Curve3>`; the
//!   Darboux Init overload (Path2d + Support) has no Curve3 view of the
//!   Adaptor3d_CurveOnSurface in the rcad data model — the overload keeps
//!   the OCCT failure path (GAP note at the body, darboux.rs precedent).
//! - `handle(GeomFill_LocationLaw/SectionLaw)` map to
//!   `Rc<RefCell<...>>` (shared mutable handles, the sweep_function.rs
//!   convention).
//! - GAP carriers: [`GeomFillAppSweep`] (TKGeomAlgo/GeomFill) for the
//!   ApproxSurf approximation and the shared [`ApproxSweepApproximation`]
//!   (sweep.rs) for the Perform(Tol) path; [`GeomFillLine`] is the trivial
//!   GeomFill_Line data re-host.
//! - `GeomFill_Trihedron` is the existing [`Trihedron`] enum
//!   (corrected_frenet.rs).

use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

use glam::{DAffine3, DVec3};

use rcad_kernel::base::geom_lib::axe_of_inertia;
use rcad_kernel::base::geom_lprop::CLProps;
use rcad_kernel::core::precision::{APPROXIMATION, CONFUSION, PCONFUSION};
use rcad_kernel::geom::{
    transform_curve, Circle3, Curve3, CurveEval, CylindricalSurface, Line3, Surface3, SurfaceEval,
    ToroidalSurface, TrimmedCurve3, TrimmedSurface,
};
use rcad_kernel::math::gp::{Ax1, Ax2, Ax3};

use super::constant_bi_normal::ConstantBiNormal;
use super::corrected_frenet::{CorrectedFrenet, Trihedron};
use super::curve_and_trihedron::CurveAndTrihedron;
use super::circular_blend_func::CircularBlendFunc;
use super::frenet::Frenet;
use super::fixed::Fixed;
use super::location_law::LocationLaw;
use super::nsections::NSections;
use super::section_law::SectionLaw;
use super::sweep::Sweep;
use super::sweep::ApproxSweepApproximation;
use super::sweep::GeomFillApproxStyle;
use super::sweep_section_generator::SweepSectionGenerator;
use super::trihedron_law::{curve_first_parameter, curve_last_parameter, PipeError, TrihedronLaw};
use super::trihedron_with_guide::TrihedronWithGuide;
use super::uniform_section::UniformSection;

// ---------------------------------------------------------------------------
// File statics and pure-math gp re-hosts
// ---------------------------------------------------------------------------

/// OCCT GeomFill_Line (TKGeomAlgo/GeomFill) — the trivial section-count
/// carrier consumed by GeomFill_AppSweep::Perform.
pub struct GeomFillLine {
    /// OCCT int myNbSections.
    my_nb_sections: usize,
}

impl GeomFillLine {
    /// OCCT GeomFill_Line::GeomFill_Line(NbSections).
    pub fn new(nb_sections: usize) -> Self {
        GeomFillLine {
            my_nb_sections: nb_sections,
        }
    }

    /// OCCT GeomFill_Line::NbSections().
    pub fn nb_sections(&self) -> usize {
        self.my_nb_sections
    }
}

/// GAP carrier: OCCT GeomFill_AppSweep (TKGeomAlgo/GeomFill) — the
/// sweeping-line approximation is not translated; construction/Perform keep
/// the OCCT failure path.
pub struct GeomFillAppSweep;

impl GeomFillAppSweep {
    /// OCCT GeomFill_AppSweep(DegMin, DegMax, T3d, T2d, NbIt, WithParameters).
    pub fn new(
        _deg_min: i32,
        _deg_max: i32,
        _t3d: f64,
        _t2d: f64,
        _nb_it: i32,
        _with_parameters: bool,
    ) -> Self {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::Perform(Line, Section, NbIterations).
    pub fn perform(&mut self, _line: &GeomFillLine, _section: &mut SweepSectionGenerator, _nb_iterations: i32) {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfShape(...).
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
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfPoles().
    pub fn surf_poles(&self) -> Vec<Vec<DVec3>> {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfWeights().
    pub fn surf_weights(&self) -> Vec<Vec<f64>> {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfUKnots().
    pub fn surf_u_knots(&self) -> Vec<f64> {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfVKnots().
    pub fn surf_v_knots(&self) -> Vec<f64> {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfUMults().
    pub fn surf_u_mults(&self) -> Vec<i32> {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::SurfVMults().
    pub fn surf_v_mults(&self) -> Vec<i32> {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::UDegree().
    pub fn u_degree(&self) -> i32 {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::VDegree().
    pub fn v_degree(&self) -> i32 {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }

    /// OCCT GeomFill_AppSweep::TolReached(Tol3d, Tol2d).
    pub fn tol_reached(&self, _tol3d: &mut f64, _tol2d: &mut f64) {
        panic!("GAP: GeomFill_AppSweep (TKGeomAlgo/GeomFill) is not translated — see file header")
    }
}

/// OCCT gp_Vec::AngleWithRef (gp_Dir.cxx L55-84) — the same signed-angle
/// re-host as sweep_section_generator (local copy over foreign values).
fn gp_vec_angle_with_ref(v: DVec3, other: DVec3, vref: DVec3) -> f64 {
    let xyz = v.cross(other);
    let cosinus = v.dot(other);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT gp_Dir::IsParallel + same sense (gp_Dir.hxx IsEqual = IsParallel &&
/// !IsOpposite).
fn gp_dir_is_equal(d1: DVec3, d2: DVec3, tol: f64) -> bool {
    let parallel = d1.cross(d2).length() <= tol * d1.length() * d2.length();
    let same_sense = d1.dot(d2) >= 0.0;
    parallel && same_sense
}

/// OCCT gp_Lin::Contains(P, Tol) — the distance of P to the line is within
/// Tol.
fn gp_lin_contains(origin: DVec3, direction: DVec3, p: DVec3, tol: f64) -> bool {
    (p - origin).cross(direction).length() <= tol
}

/// OCCT gp_Ax3::Rotate(Ax1, Angle) — the frame rotated about the axis
/// (Rodrigues over the location and the direction vectors).
fn ax3_rotated(a: &Ax3, axis: &Ax1, angle: f64) -> Ax3 {
    let rotate = |v: DVec3| -> DVec3 {
        let n = axis.direction.normalize_or_zero();
        let c = angle.cos();
        let s = angle.sin();
        v * c + n.cross(v) * s + n * (n.dot(v)) * (1.0 - c)
    };
    let loc = axis.location + rotate(a.axis.location - axis.location);
    Ax3::from_pnt_n_vx(loc, rotate(a.direction()), rotate(a.x_direction))
}

/// OCCT ElCLib::CircleDN(U, Pos, Radius, 1) — the first derivative of the
/// circle: Radius * (-sin(U) X + cos(U) Y) over the Ax2 frame (pure math).
fn el_clib_circle_dn(pos: &Ax2, radius: f64, u: f64) -> DVec3 {
    let x = pos.x_direction;
    let y = pos.y_direction;
    radius * (u.cos() * y - u.sin() * x)
}

/// OCCT Geom_Curve::Reverse() over the rcad Curve3 forms consumed by
/// CheckSense (Line flips the direction; Circle negates YDirection; BSpline
/// reverses poles/knots; Trimmed reverses the basis and swaps the bounds).
fn reverse_curve3(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin,
            direction: -l.direction,
        }),
        Curve3::Circle(circ) => Curve3::Circle(Circle3 {
            y_dir: -circ.y_dir,
            ..*circ
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(bs.reversed()),
        Curve3::Trimmed(tc) => Curve3::Trimmed(TrimmedCurve3::new(
            reverse_curve3(&tc.curve),
            tc.last,
            tc.first,
        )),
        other => other.clone(),
    }
}

/// OCCT static CheckSense (GeomFill_Pipe.cxx L79-221) — the section
/// orientation check/reversal over a sequence of curves.
fn check_sense(seq1: &[Curve3], seq2: &mut Vec<Curve3>) -> bool {
    // initialisation
    let mut no_sing = true;
    seq2.clear();

    let c1 = &seq1[0];
    let mut f = curve_first_parameter(c1);
    let mut l = curve_last_parameter(c1);
    let np = 21usize;
    let mut tab = vec![DVec3::ZERO; np];
    let mut u = f;
    let h = (f - l).abs() / 20.0;
    for cell in tab.iter_mut() {
        *cell = c1.point_at(u);
        u += h;
        if (u - f) * (u - l) > 0.0 {
            u = l;
        }
    }
    let mut axe_ref = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
    let mut axe = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
    let mut pos;
    let mut sing = false;
    axe_of_inertia(&tab, &mut axe_ref, &mut sing, CONFUSION);

    // si la section est une droite, ca ne marche pas
    if sing {
        no_sing = false;
    }

    pos = axe_ref.location;
    let (mut alpha1, mut alpha2, mut alpha3);
    let mut p1;
    let mut p2;
    u = (f + l - h) / 2.0 - h;
    p1 = c1.point_at(u);
    u += h;
    p2 = c1.point_at(u);
    alpha1 = gp_vec_angle_with_ref(p1 - pos, p2 - pos, axe_ref.direction);
    p1 = p2;
    u += h;
    p2 = c1.point_at(u);
    alpha2 = gp_vec_angle_with_ref(p1 - pos, p2 - pos, axe_ref.direction);
    p1 = p2;
    u += h;
    p2 = c1.point_at(u);
    alpha3 = gp_vec_angle_with_ref(p1 - pos, p2 - pos, axe_ref.direction);
    seq2.push(c1.clone());

    for iseq in 1..seq1.len() {
        // discretisation de C2
        let c2 = &seq1[iseq];
        f = curve_first_parameter(c2);
        l = curve_last_parameter(c2);
        u = f;
        for cell in tab.iter_mut() {
            *cell = c2.point_at(u);
            u += h;
            if (u - f) * (u - l) > 0.0 {
                u = l;
            }
        }
        axe_of_inertia(&tab, &mut axe, &mut sing, CONFUSION);

        // si la section est une droite, ca ne marche pas
        if sing {
            no_sing = false;
        }

        pos = axe.location;
        u = (f + l - h) / 2.0 - h;
        p1 = c2.point_at(u);
        u += h;
        p2 = c2.point_at(u);
        let beta1 = gp_vec_angle_with_ref(p1 - pos, p2 - pos, axe_ref.direction);
        p1 = p2;
        u += h;
        p2 = c2.point_at(u);
        let beta2 = gp_vec_angle_with_ref(p1 - pos, p2 - pos, axe_ref.direction);
        p1 = p2;
        u += h;
        p2 = c2.point_at(u);
        let beta3 = gp_vec_angle_with_ref(p1 - pos, p2 - pos, axe_ref.direction);

        // meme sens ?
        let mut ok = true;
        let pasnul1 = alpha1.abs() > CONFUSION && beta1.abs() > CONFUSION;
        let pasnul2 = alpha2.abs() > CONFUSION && beta2.abs() > CONFUSION;
        let pasnul3 = alpha3.abs() > CONFUSION && beta3.abs() > CONFUSION;
        if pasnul1 && pasnul2 && pasnul3 {
            if alpha1 * beta1 > 0.0 {
                ok = alpha2 * beta2 > 0.0 || alpha3 * beta3 > 0.0;
            } else {
                ok = alpha2 * beta2 > 0.0 && alpha3 * beta3 > 0.0;
            }
        } else if pasnul1 && pasnul2 && !pasnul3 {
            ok = alpha1 * beta1 > 0.0 || alpha2 * beta2 > 0.0;
        } else if pasnul1 && !pasnul2 && pasnul3 {
            ok = alpha1 * beta1 > 0.0 || alpha3 * beta3 > 0.0;
        } else if !pasnul1 && pasnul2 && pasnul3 {
            ok = alpha2 * beta2 > 0.0 || alpha3 * beta3 > 0.0;
        } else if pasnul1 {
            ok = alpha1 * beta1 > 0.0;
        } else if pasnul2 {
            ok = alpha2 * beta2 > 0.0;
        } else if pasnul3 {
            ok = alpha3 * beta3 > 0.0;
        }

        let mut c2 = c2.clone();
        if no_sing && !ok {
            // OCCT: C2->Reverse().
            c2 = reverse_curve3(&c2);
        }
        seq2.push(c2);
    }

    no_sing
}

/// OCCT GeomFill_Pipe (GeomFill_Pipe.hxx members L304-316).
pub struct Pipe {
    /// OCCT GeomFill_PipeError myStatus.
    my_status: PipeError,
    /// OCCT double myRadius.
    my_radius: f64,
    /// OCCT double myError.
    my_error: f64,
    /// OCCT handle(Adaptor3d_Curve) myAdpPath.
    my_adp_path: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myAdpFirstSect.
    my_adp_first_sect: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myAdpLastSect.
    my_adp_last_sect: Option<Curve3>,
    /// OCCT handle(Geom_Surface) mySurface.
    my_surface: Option<Surface3>,
    /// OCCT handle(GeomFill_LocationLaw) myLoc.
    my_loc: Option<Rc<RefCell<dyn LocationLaw>>>,
    /// OCCT handle(GeomFill_SectionLaw) mySec.
    my_sec: Option<Rc<RefCell<dyn SectionLaw>>>,
    /// OCCT int myType.
    my_type: i32,
    /// OCCT bool myExchUV.
    my_exch_uv: bool,
    /// OCCT bool myKPart.
    my_k_part: bool,
    /// OCCT bool myPolynomial.
    my_polynomial: bool,
}

impl Pipe {
    /// OCCT GeomFill_Pipe::GeomFill_Pipe() (L228-234) + Init() (L414-426).
    pub fn new() -> Self {
        let mut pipe = Pipe {
            my_status: PipeError::PipeNotOk,
            my_radius: 0.0,
            my_error: 0.0,
            my_adp_path: None,
            my_adp_first_sect: None,
            my_adp_last_sect: None,
            my_surface: None,
            my_loc: None,
            my_sec: None,
            my_type: 0,
            my_exch_uv: false,
            my_k_part: false,
            my_polynomial: false,
        };
        pipe.init_defaults();
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, Radius) (L238-245).
    pub fn new_with_radius(path: &Curve3, radius: f64) -> Self {
        let mut pipe = Pipe::new();
        pipe.init_with_radius(path, radius);
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, FirstSect, Option) (L249-258).
    pub fn new_with_trihedron(path: &Curve3, first_sect: &Curve3, option: Trihedron) -> Self {
        let mut pipe = Pipe::new();
        pipe.init_with_trihedron(path, first_sect, option);
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, FirstSect, LastSect) (L275-284).
    pub fn new_with_sections(path: &Curve3, first_sect: &Curve3, last_sect: &Curve3) -> Self {
        let mut pipe = Pipe::new();
        pipe.init_with_sections(path, first_sect, last_sect);
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, NSections) (L288-296).
    pub fn new_with_nsections(path: &Curve3, nsections: &[Curve3]) -> Self {
        let mut pipe = Pipe::new();
        pipe.init_with_nsections(path, nsections);
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, Curve1, Direction) (L300-308).
    pub fn new_with_direction(path: &Curve3, curve1: &Curve3, direction: DVec3) -> Self {
        let mut pipe = Pipe {
            my_status: PipeError::PipeNotOk,
            my_radius: 0.0,
            my_error: 0.0,
            my_adp_path: None,
            my_adp_first_sect: None,
            my_adp_last_sect: None,
            my_surface: None,
            my_loc: None,
            my_sec: None,
            my_type: 0,
            my_exch_uv: false,
            my_k_part: false,
            my_polynomial: false,
        };
        pipe.init_with_direction(path, curve1, direction);
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, Curve1, Curve2, Radius) (L312-326).
    pub fn new_with_adaptors(
        path: &Curve3,
        curve1: &Curve3,
        curve2: &Curve3,
        radius: f64,
    ) -> Self {
        let mut pipe = Pipe::new();
        pipe.init_with_adaptors(path, curve1, curve2, radius);
        pipe
    }

    /// OCCT GeomFill_Pipe(Path, Guide, FirstSect, byACR, rotat) (L346-361).
    pub fn new_with_guide(
        path: &Curve3,
        guide: &Curve3,
        first_sect: &Curve3,
        by_acr: bool,
        rotat: bool,
    ) -> Self {
        let mut pipe = Pipe::new();
        pipe.init_with_guide(path, guide, first_sect, by_acr, rotat);
        pipe
    }

    /// OCCT Init(Path, Guide, FirstSect, byACR, rotat) (L365-410) — Path:
    /// trajectoire, Guide: courbe guide, FirstSect: section, rotat: vrai si
    /// on veu la rotation.
    pub fn init_with_guide(
        &mut self,
        path: &Curve3,
        guide: &Curve3,
        first_sect: &Curve3,
        by_acr: bool,
        rotat: bool,
    ) {
        let mut angle = 0.0;
        self.my_adp_path = Some(path.clone());

        // loi de triedre.
        let t_law: Box<dyn TrihedronWithGuide> = if by_acr {
            let mut law = super::guide_trihedron_ac::GuideTrihedronAC::new(guide);
            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            TrihedronLaw::set_curve(&mut law, adp_path);
            Box::new(law)
        } else {
            let mut law = super::guide_trihedron_plan::GuideTrihedronPlan::new(guide);
            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            TrihedronLaw::set_curve(&mut law, adp_path);
            Box::new(law)
        };

        // loi de positionnement.
        let the_loc = Rc::new(RefCell::new(super::location_guide::LocationGuide::new(t_law)));
        {
            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            let mut loc = the_loc.borrow_mut();
            LocationLaw::set_curve(&mut *loc, adp_path);
        }

        let mut place = super::section_placement::SectionPlacement::new_loc(
            the_loc.clone(),
            first_sect,
        );
        place.perform_confusion(CONFUSION);

        // loi de section.
        self.my_sec = Some(Rc::new(RefCell::new(UniformSection::new(
            &place.section(false),
            curve_first_parameter(self.my_adp_path.as_ref().expect("null myAdpPath")),
            curve_last_parameter(self.my_adp_path.as_ref().expect("null myAdpPath")),
        ))));

        if rotat {
            let f = curve_first_parameter(self.my_adp_path.as_ref().expect("null myAdpPath"));
            let l = curve_last_parameter(self.my_adp_path.as_ref().expect("null myAdpPath"));
            the_loc.borrow_mut().set(
                self.my_sec.as_ref().expect("null mySec").clone(),
                rotat,
                f,
                l,
                0.0,
                &mut angle,
            );
        }
        self.my_loc = Some(the_loc);
    }

    /// OCCT Init() (L414-426).
    pub fn init_defaults(&mut self) {
        self.my_type = 0;
        self.my_error = 0.0;
        self.my_radius = 0.0;
        self.my_k_part = true;
        self.my_polynomial = false;
        self.my_adp_path = None;
        self.my_adp_first_sect = None;
        self.my_adp_last_sect = None;
        self.my_loc = None;
        self.my_sec = None;
    }

    /// OCCT Init(Path, Radius) (L430-446).
    pub fn init_with_radius(&mut self, path: &Curve3, radius: f64) {
        // Ancienne methode
        self.my_type = 1;
        self.my_error = 0.0;
        self.my_radius = radius;

        // Nouvelle methode
        self.my_adp_path = Some(path.clone());
        // OCCT: C = new Geom_Circle(gp::XOY(), Radius); C->Rotate(gp::OZ(),
        // M_PI / 2.) — the X direction rotated onto Y.
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::Y,
            y_dir: -DVec3::X,
            radius,
        };
        self.my_sec = Some(Rc::new(RefCell::new(UniformSection::new(
            &Curve3::Circle(c),
            curve_first_parameter(path),
            curve_last_parameter(path),
        ))));
        let t_law = CorrectedFrenet::new();
        self.my_loc = Some(Rc::new(RefCell::new(CurveAndTrihedron::new(Box::new(
            t_law,
        )))));
        if let Some(adp_path) = self.my_adp_path.clone() {
            let mut loc = self.my_loc.as_mut().expect("null myLoc").borrow_mut();
            LocationLaw::set_curve(&mut *loc, adp_path);
        }
    }

    /// OCCT Init(Path, FirstSect, Option) (L450-557).
    pub fn init_with_trihedron(&mut self, path: &Curve3, first_sect: &Curve3, option: Trihedron) {
        let mut sect: Option<Curve3> = None;
        self.my_adp_path = Some(path.clone());
        let mut param = curve_first_parameter(path);

        // Construction de la loi de triedre.
        let t_law: Option<Box<dyn TrihedronLaw>> = match option {
            Trihedron::IsCorrectedFrenet => Some(Box::new(CorrectedFrenet::new())),

            Trihedron::IsDarboux | Trihedron::IsFrenet => Some(Box::new(Frenet::new())),

            Trihedron::IsFixed => {
                let eps = 1.0e-9;
                let mut v1 = DVec3::new(0.0, 0.0, 1.0);
                let mut v2 = DVec3::new(0.0, 1.0, 0.0);
                let mut cp = CLProps::with_param(path, param, 2, eps);
                if cp.is_tangent_defined() {
                    let d = cp.tangent().expect("tangent");
                    v1 = d;
                    v1 = v1.normalize_or_zero();
                    if cp.curvature() > eps {
                        if let Some(n) = cp.normal() {
                            v2 = n.normalize_or_zero();
                        }
                    } else {
                        let p0 = DVec3::ZERO;
                        let axe = Ax2::new(p0, d, DVec3::X);
                        let d = axe.x_direction;
                        v2 = d.normalize_or_zero();
                    }
                }
                Some(Box::new(Fixed::new(v1, v2)))
            }

            Trihedron::IsConstantNormal => {
                let frenet = Frenet::new();
                let mut loc = CurveAndTrihedron::new(Box::new(frenet));
                let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
                LocationLaw::set_curve(&mut loc, adp_path);
                let mut place = super::section_placement::SectionPlacement::new_loc(
                    Rc::new(RefCell::new(loc)),
                    first_sect,
                );
                place.perform_confusion(CONFUSION);
                let ponsec = place.parameter_on_section();

                let eps = 1.0e-9;
                let mut v = DVec3::new(0.0, 1.0, 0.0);
                let mut cp = CLProps::with_param(first_sect, ponsec, 2, eps);
                if cp.is_tangent_defined() {
                    let d = cp.tangent().expect("tangent");
                    if cp.curvature() > eps {
                        if let Some(n) = cp.normal() {
                            v = n.normalize_or_zero();
                        }
                    } else {
                        let p0 = DVec3::ZERO;
                        let axe = Ax2::new(p0, d, DVec3::X);
                        let d = axe.x_direction;
                        v = d.normalize_or_zero();
                    }
                }
                Some(Box::new(ConstantBiNormal::new(v)))
            }

            _ => {
                panic!("Standard_ConstructionError: GeomFill::Init : Unknown Option");
            }
        };

        if let Some(t_law) = t_law {
            self.my_loc = Some(Rc::new(RefCell::new(CurveAndTrihedron::new(t_law))));
            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            {
                let mut loc = self.my_loc.as_mut().expect("null myLoc").borrow_mut();
                LocationLaw::set_curve(&mut *loc, adp_path);
            }
            let mut place = super::section_placement::SectionPlacement::new_loc(
                self.my_loc.as_ref().expect("null myLoc").clone(),
                first_sect,
            );
            place.perform_confusion(CONFUSION);
            param = place.parameter_on_path();
            sect = Some(place.section(false));

            self.my_sec = Some(Rc::new(RefCell::new(UniformSection::new(
                &sect.expect("null Sect"),
                curve_first_parameter(path),
                curve_last_parameter(path),
            ))));
        }
        let _ = param;
    }

    /// OCCT Init(Path2d, Support, FirstSect) (L564-581) — the Darboux sweep
    /// along a pcurve.  GAP note: the Adaptor3d_CurveOnSurface has no Curve3
    /// view in the rcad data model (darboux.rs carries the pair separately
    /// on the law), so the myAdpPath installation and the downstream
    /// Perform(myAdpPath, ...) overload cannot be represented — the OCCT
    /// failure path is kept pending the kernel Curve3 carrier.
    pub fn init_with_pcurve(&mut self, _path2d: &rcad_kernel::geom::Curve2d, _support: &Surface3, _first_sect: &Curve3) {
        panic!(
            "GAP: GeomFill_Pipe::Init(Geom2d_Curve, Geom_Surface, Geom_Curve) — Adaptor3d_CurveOnSurface has no rcad Curve3 view (darboux.rs precedent)"
        )
    }

    /// OCCT Init(Path, FirstSect, Direction) (L585-604).
    pub fn init_with_direction(&mut self, path: &Curve3, first_sect: &Curve3, direction: DVec3) {
        self.init_defaults();

        let sect;
        self.my_adp_path = Some(path.clone());
        let v = direction;
        let t_law = ConstantBiNormal::new(v);

        let loc = Rc::new(RefCell::new(CurveAndTrihedron::new(Box::new(t_law))));
        {
            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            let mut loc_mut = loc.borrow_mut();
            LocationLaw::set_curve(&mut *loc_mut, adp_path);
        }
        self.my_loc = Some(loc.clone());
        let mut place =
            super::section_placement::SectionPlacement::new_loc(loc, first_sect);
        place.perform_confusion(CONFUSION);
        sect = place.section(false);

        self.my_sec = Some(Rc::new(RefCell::new(UniformSection::new(
            &sect,
            curve_first_parameter(path),
            curve_last_parameter(path),
        ))));
    }

    /// OCCT Init(Path, NSections) (L608-675).
    pub fn init_with_nsections(&mut self, path: &Curve3, nsections: &[Curve3]) {
        self.my_type = 3;
        self.my_error = 0.0;
        self.my_radius = 0.0;

        let t_law: Box<dyn TrihedronLaw> = Box::new(CorrectedFrenet::new());
        self.my_adp_path = Some(path.clone());
        {
            let mut loc = CurveAndTrihedron::new(t_law);
            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            LocationLaw::set_curve(&mut loc, adp_path);
            let loc = Rc::new(RefCell::new(loc));
            self.my_loc = Some(loc.clone());
            let mut seq_c: Vec<Curve3> = Vec::new();
            let mut seq_p: Vec<f64> = Vec::new();
            for nsec in nsections {
                let mut place = super::section_placement::SectionPlacement::new_loc(
                    loc.clone(),
                    nsec,
                );
                place.perform_confusion(CONFUSION);
                seq_p.push(place.parameter_on_path());
                seq_c.push(place.section(false));
            }

            // verification des orientations
            let mut new_seq: Vec<Curve3> = Vec::new();
            if check_sense(&seq_c, &mut new_seq) {
                seq_c = new_seq;
            }

            // verification des parametres
            let mut play_again = true;
            while play_again {
                play_again = false;
                for i in 0..nsections.len() {
                    for j in i..nsections.len() {
                        if seq_p[i] > seq_p[j] {
                            seq_p.swap(i, j);
                            seq_c.swap(i, j);
                            play_again = true;
                        }
                    }
                }
            }
            for i in 0..nsections.len() - 1 {
                if (seq_p[i + 1] - seq_p[i]).abs() < PCONFUSION {
                    panic!(
                        "Standard_ConstructionError: GeomFill_Pipe::Init with NSections : invalid parameters"
                    );
                }
            }

            // creation de la NSections
            let first = curve_first_parameter(path);
            let last = curve_last_parameter(path);
            let deb = curve_first_parameter(&seq_c[0]);
            let fin = curve_last_parameter(&seq_c[0]);
            self.my_sec = Some(Rc::new(RefCell::new(NSections::new_with_bounds(
                seq_c, seq_p, deb, fin, first, last,
            ))));
        }
    }

    /// OCCT Init(Path, FirstSect, LastSect) (L679-729).
    pub fn init_with_sections(&mut self, path: &Curve3, first_sect: &Curve3, last_sect: &Curve3) {
        self.my_type = 3;
        self.my_error = 0.0;
        self.my_radius = 0.0;
        let first = curve_first_parameter(path);
        let last = curve_last_parameter(path);
        let t_law: Box<dyn TrihedronLaw> = Box::new(CorrectedFrenet::new());
        self.my_adp_path = Some(path.clone());

        {
            let mut loc = CurveAndTrihedron::new(t_law);

            let adp_path = self.my_adp_path.clone().expect("null myAdpPath");
            if !LocationLaw::set_curve(&mut loc, adp_path) {
                self.my_status = PipeError::ImpossibleContact;
                return;
            }
            let loc = Rc::new(RefCell::new(loc));
            self.my_loc = Some(loc.clone());

            let mut seq_c: Vec<Curve3> = Vec::new();
            let mut seq_p: Vec<f64> = Vec::new();
            // sequence of sections
            let mut pl1 = super::section_placement::SectionPlacement::new_loc(loc.clone(), first_sect);
            pl1.perform_at(first, CONFUSION);
            seq_c.push(pl1.section(false));
            let mut pl2 = super::section_placement::SectionPlacement::new_loc(loc.clone(), last_sect);
            pl2.perform_at(first, CONFUSION);
            seq_c.push(pl2.section(false));
            // sequence of associated parameters
            seq_p.push(first);
            seq_p.push(last);

            // orientation verification
            let mut new_seq: Vec<Curve3> = Vec::new();
            if check_sense(&seq_c, &mut new_seq) {
                seq_c = new_seq;
            }

            // creation of the NSections
            let deb = curve_first_parameter(&seq_c[0]);
            let fin = curve_last_parameter(&seq_c[0]);
            self.my_sec = Some(Rc::new(RefCell::new(NSections::new_with_bounds(
                seq_c, seq_p, deb, fin, first, last,
            ))));
        }
    }

    /// OCCT Init(Path, Curve1, Curve2, Radius) (L733-744).
    pub fn init_with_adaptors(&mut self, path: &Curve3, curve1: &Curve3, curve2: &Curve3, radius: f64) {
        self.my_type = 4;
        self.my_error = 0.0;
        self.my_radius = radius;
        self.my_adp_path = Some(path.clone());
        self.my_adp_first_sect = Some(curve1.clone());
        self.my_adp_last_sect = Some(curve2.clone());
    }

    /// OCCT Perform(WithParameters, Polynomial) (L748-768).
    pub fn perform(&mut self, with_parameters: bool, polynomial: bool) {
        if self.my_loc.is_some() && self.my_sec.is_some() {
            // OCCT: Perform(1.e-4, Polynomial) — the defaulted
            // Conti = GeomAbs_C1, MaxDegree = 11, NbMaxSegment = 30.
            self.perform_with_tolerance(
                1.0e-4,
                polynomial,
                rcad_kernel::math::GeomAbsShape::C1,
                11,
                30,
            );
            return;
        }

        self.my_polynomial = polynomial;
        // on traite la cas tuyau sur arete Type = 4
        if self.my_polynomial {
            self.approx_surf(with_parameters);
            return;
        }
        if !self.k_part_t4() {
            self.approx_surf(with_parameters);
        }
    }

    /// OCCT Perform(Tol, Polynomial, Conti, DegMax, NbMaxSegment) (L772-865).
    pub fn perform_with_tolerance(
        &mut self,
        tol: f64,
        polynomial: bool,
        conti: rcad_kernel::math::GeomAbsShape,
        degmax: i32,
        nb_max_segment: i32,
    ) {
        use rcad_kernel::math::GeomAbsShape;
        if self.my_status == PipeError::ImpossibleContact {
            return;
        }

        // OCCT: the G1/G1 and G2/G2 pairs fold onto C1/C2 (the rcad
        // GeomAbsShape carries the C-forms only).
        let the_conti = match conti {
            GeomAbsShape::C0 => GeomAbsShape::C0,
            GeomAbsShape::C1 => GeomAbsShape::C1,
            GeomAbsShape::C2 => GeomAbsShape::C2,
            _ => GeomAbsShape::C2, // On ne sait pas faire mieux !
        };

        if self.my_type == 4 {
            if !self.k_part_t4() {
                let func = CircularBlendFunc::new(
                    self.my_adp_path.as_ref().expect("null myAdpPath"),
                    self.my_adp_first_sect.as_ref().expect("null myAdpFirstSect"),
                    self.my_adp_last_sect.as_ref().expect("null myAdpLastSect"),
                    self.my_radius,
                    polynomial,
                );

                let mut app = ApproxSweepApproximation::new(&func);
                let adp_path = self.my_adp_path.as_ref().expect("null myAdpPath");
                app.perform(
                    curve_first_parameter(adp_path),
                    curve_last_parameter(adp_path),
                    tol,
                    tol,
                    0.0,
                    0.01,
                    the_conti,
                    degmax,
                    nb_max_segment,
                );
                if app.is_done() {
                    // OCCT: new Geom_BSplineSurface(App.SurfPoles(), ...).
                    let u_knots = app.surf_u_knots();
                    let v_knots = app.surf_v_knots();
                    let mut flat_u = Vec::new();
                    for (k, m) in u_knots.iter().zip(app.surf_u_mults().iter()) {
                        for _ in 0..*m {
                            flat_u.push(*k);
                        }
                    }
                    let mut flat_v = Vec::new();
                    for (k, m) in v_knots.iter().zip(app.surf_v_mults().iter()) {
                        for _ in 0..*m {
                            flat_v.push(*k);
                        }
                    }
                    self.my_surface = Some(Surface3::BSpline(rcad_kernel::geom::BSplineSurface {
                        degree_u: app.u_degree() as usize,
                        degree_v: app.v_degree() as usize,
                        knots_u: flat_u,
                        knots_v: flat_v,
                        control_points: app.surf_poles(),
                        weights: app.surf_weights(),
                        is_periodic_u: false,
                        is_periodic_v: false,
                    }));
                    self.my_error = app.max_error_on_surf();
                    self.my_status = PipeError::PipeOk;
                }
            }
        } else if self.my_loc.is_some() && self.my_sec.is_some() {
            let mut sweep = Sweep::new(self.my_loc.as_ref().expect("null myLoc").clone(), self.my_k_part);
            sweep.set_tolerance_sweep(tol, 1.0, 1.0e-5, 1.0);
            sweep.build(
                self.my_sec.as_ref().expect("null mySec").clone(),
                // (the Rc handle itself is cloned)
                GeomFillApproxStyle::GeomFill_Location,
                the_conti,
                degmax,
                nb_max_segment,
            );
            if sweep.is_done() {
                self.my_surface = sweep.surface().cloned();
                self.my_error = sweep.error_on_surface();
                self.my_status = PipeError::PipeOk;
            }
        } else {
            self.perform(true, polynomial);
        }
    }

    /// OCCT KPartT4 (L871-1010).
    pub fn k_part_t4(&mut self) -> bool {
        let mut ok = false;
        let adp_path = self.my_adp_path.as_ref().expect("null myAdpPath");
        let adp_first = self.my_adp_first_sect.as_ref().expect("null myAdpFirstSect");
        let adp_last = self.my_adp_last_sect.as_ref().expect("null myAdpLastSect");

        // -------    Cas du Cylindre  --------------------------
        let is_line3 = |c: &Curve3| matches!(c, Curve3::Line(_));
        let is_circle3 = |c: &Curve3| matches!(c, Curve3::Circle(_));
        if is_line3(adp_path) && is_line3(adp_first) && is_line3(adp_last) {
            // try to generate a cylinder.
            let (o0, d0v) = match adp_path {
                Curve3::Line(l) => (l.origin, l.direction),
                _ => unreachable!(),
            };
            let (o1, d1v) = match adp_first {
                Curve3::Line(l) => (l.origin, l.direction),
                _ => unreachable!(),
            };
            let (o2, d2v) = match adp_last {
                Curve3::Line(l) => (l.origin, l.direction),
                _ => unreachable!(),
            };
            let d0 = d0v.normalize_or_zero();
            let d1 = d1v.normalize_or_zero();
            let d2 = d2v.normalize_or_zero();
            // direction must be the same.
            if !gp_dir_is_equal(d0, d1, 1.0e-12) || !gp_dir_is_equal(d1, d2, 1.0e-12) {
                return ok;
            }

            // the length of the line must be te same
            let l0 = curve_last_parameter(adp_path) - curve_first_parameter(adp_path);
            let l1 = curve_last_parameter(adp_first) - curve_first_parameter(adp_first);
            let l2 = curve_last_parameter(adp_last) - curve_first_parameter(adp_last);
            if (l1 - l0).abs() > CONFUSION || (l2 - l0).abs() > CONFUSION {
                return ok;
            }

            // the first points must be normal to the path.
            let p0 = adp_path.point_at(curve_first_parameter(adp_path));
            let p1 = adp_first.point_at(curve_first_parameter(adp_first));
            let p2 = adp_last.point_at(curve_first_parameter(adp_last));
            let v1 = (p1 - p0).normalize_or_zero();
            let v2 = (p2 - p0).normalize_or_zero();
            if v1.dot(d0).abs() > CONFUSION || v2.dot(d0).abs() > CONFUSION {
                return ok;
            }

            // the result is a cylindrical surface.
            let x = v1;
            let y = v2;
            let z_ref = x.cross(y);

            let mut axis = Ax3::from_pnt_n_vx(o0, d0, x);
            if z_ref.dot(d0) < 0.0 {
                // OCCT: Axis.YReverse().
                axis.y_direction = -axis.y_direction;
            }

            // rotate the surface to set the iso U = 0 not in the result.
            axis = ax3_rotated(&axis, &Ax1::new(p0, z_ref), -PI / 2.0);

            let surface = Surface3::Cylinder(CylindricalSurface {
                origin: axis.axis.location,
                axis: axis.direction(),
                radius: self.my_radius,
                ref_dir: axis.x_direction,
                y_dir: None,
            });
            let alpha = gp_vec_angle_with_ref(v1, v2, z_ref);
            self.my_surface = Some(Surface3::Trimmed(TrimmedSurface {
                basis: Box::new(surface),
                trim: [
                    PI / 2.0,
                    PI / 2.0 + alpha,
                    curve_first_parameter(adp_path),
                    curve_last_parameter(adp_path),
                ],
            }));
            ok = true; // C'est bien un cylindre
            self.my_status = PipeError::PipeOk;
        }
        // -----------    Cas du tore  ----------------------------------
        else if is_circle3(adp_path) && is_circle3(adp_first) && is_circle3(adp_last) {
            // try to generate a toroidal surface.
            // les 3 cercles doivent avoir meme angle d'ouverture
            let alp0 = curve_first_parameter(adp_path) - curve_last_parameter(adp_path);
            let alp1 = curve_first_parameter(adp_first) - curve_last_parameter(adp_first);
            let alp2 = curve_first_parameter(adp_last) - curve_last_parameter(adp_last);

            if (alp0 - alp1).abs() > 1.0e-12 || (alp0 - alp2).abs() > 1.0e-12 {
                return ok;
            }

            let c_path = match adp_path {
                Curve3::Circle(c) => c.clone(),
                _ => unreachable!(),
            };
            let a0 = Ax2::new(c_path.center, c_path.normal, c_path.x_dir);
            let c1 = match adp_first {
                Curve3::Circle(c) => c.clone(),
                _ => unreachable!(),
            };
            let a1 = Ax2::new(c1.center, c1.normal, c1.x_dir);
            let c2 = match adp_last {
                Curve3::Circle(c) => c.clone(),
                _ => unreachable!(),
            };
            let a2 = Ax2::new(c2.center, c2.normal, c2.x_dir);
            let d0 = a0.direction;
            let d1 = a1.direction;
            let d2 = a2.direction;
            let p0 = adp_path.point_at(curve_first_parameter(adp_path));
            let p1 = adp_first.point_at(curve_first_parameter(adp_first));
            let p2 = adp_last.point_at(curve_first_parameter(adp_last));

            // les 3 directions doivent etre egales.
            if !gp_dir_is_equal(d0, d1, 1.0e-12) || !gp_dir_is_equal(d1, d2, 1.0e-12) {
                return ok;
            }

            // les 3 ax1 doivent etre confondus.
            if !gp_lin_contains(a0.location, a0.direction, a1.location, CONFUSION)
                || !gp_lin_contains(a0.location, a0.direction, a2.location, CONFUSION)
            {
                return ok;
            }

            // les 3 premiers points doivent etre dans la meme section.
            let v1 = (p1 - p0).normalize_or_zero();
            let v2 = (p2 - p0).normalize_or_zero();
            let y_ref =
                el_clib_circle_dn(&a0, c_path.radius, curve_first_parameter(adp_path));
            if v1.dot(y_ref).abs() > CONFUSION || v2.dot(y_ref).abs() > CONFUSION {
                return ok;
            }

            // OK it`s a Toroidal Surface !!  OUF !!
            let mut torus = ToroidalSurface {
                center: a0.location,
                axis: a0.direction,
                ref_dir: a0.x_direction,
                major_radius: c_path.radius,
                minor_radius: self.my_radius,
            };
            let x_ref = p0 - a0.location;
            // au maximum on fait un tore d`ouverture en V = PI
            let mut vv1 = gp_vec_angle_with_ref(v1, x_ref, y_ref);
            let mut vv2 = gp_vec_angle_with_ref(v2, x_ref, y_ref);
            let delta_v = gp_vec_angle_with_ref(v2, v1, y_ref);
            if delta_v < 0.0 {
                // OCCT: T.VReverse() — the minor-direction sense flips.
                torus.ref_dir = -torus.ref_dir;
                vv1 = -vv1;
                vv2 = 2.0 * PI + vv1 - delta_v;
            }
            self.my_surface = Some(Surface3::Trimmed(TrimmedSurface {
                basis: Box::new(Surface3::Torus(torus)),
                trim: [
                    curve_first_parameter(adp_path),
                    curve_last_parameter(adp_path),
                    vv1,
                    vv2,
                ],
            }));
            self.my_exch_uv = true;
            ok = true;
            self.my_status = PipeError::PipeOk;
        }

        ok
    }

    /// OCCT ApproxSurf (L1014-1102).
    fn approx_surf(&mut self, with_parameters: bool) {
        // Traitment of general case.  generate a sequence of the section by
        // <SweepSectionGenerator> and approximate this sequence.

        if self.my_type != 4 {
            panic!("Standard_ConstructionError: GeomFill_Pipe");
        }
        let mut section = SweepSectionGenerator::new_with_adaptors(
            self.my_adp_path.as_ref().expect("null myAdpPath"),
            self.my_adp_first_sect.as_ref().expect("null myAdpFirstSect"),
            self.my_adp_last_sect.as_ref().expect("null myAdpLastSect"),
            self.my_radius,
        );

        section.perform(self.my_polynomial);

        let line = GeomFillLine::new(section.nb_sections());
        let nb_it = 0i32;
        let t3d = APPROXIMATION;
        let t2d = APPROXIMATION * 0.01; // Precision::PApproximation()
        let mut app = GeomFillAppSweep::new(4, 8, t3d, t2d, nb_it, with_parameters);

        app.perform(&line, &mut section, 30);

        if !app.is_done() {
            // OCCT keeps the throw commented — the literal empty body.
        } else {
            let (mut u_degree, mut v_degree) = (0i32, 0i32);
            let (mut nb_u_poles, mut nb_v_poles) = (0i32, 0i32);
            let (mut nb_u_knots, mut nb_v_knots) = (0i32, 0i32);
            app.surf_shape(
                &mut u_degree,
                &mut v_degree,
                &mut nb_u_poles,
                &mut nb_v_poles,
                &mut nb_u_knots,
                &mut nb_v_knots,
            );

            let u_knots = app.surf_u_knots();
            let v_knots = app.surf_v_knots();
            let mut flat_u = Vec::new();
            for (k, m) in u_knots.iter().zip(app.surf_u_mults().iter()) {
                for _ in 0..*m {
                    flat_u.push(*k);
                }
            }
            let mut flat_v = Vec::new();
            for (k, m) in v_knots.iter().zip(app.surf_v_mults().iter()) {
                for _ in 0..*m {
                    flat_v.push(*k);
                }
            }
            self.my_surface = Some(Surface3::BSpline(rcad_kernel::geom::BSplineSurface {
                degree_u: app.u_degree() as usize,
                degree_v: app.v_degree() as usize,
                knots_u: flat_u,
                knots_v: flat_v,
                control_points: app.surf_poles(),
                weights: app.surf_weights(),
                is_periodic_u: false,
                is_periodic_v: false,
            }));
            let mut t2d = 0.0;
            let mut my_error = 0.0;
            app.tol_reached(&mut my_error, &mut t2d);
            self.my_error = my_error;
            self.my_status = PipeError::PipeOk;
        }
    }

    // -----------------------------------------------------------------
    // Accessors (GeomFill_Pipe.lxx / hxx)
    // -----------------------------------------------------------------

    /// OCCT Surface() (lxx L18-22).
    pub fn surface(&self) -> Option<&Surface3> {
        self.my_surface.as_ref()
    }

    /// OCCT ExchangeUV() (lxx L25-29).
    pub fn exchange_uv(&self) -> bool {
        self.my_exch_uv
    }

    /// OCCT GenerateParticularCase(B) (lxx L32-35).
    pub fn set_generate_particular_case(&mut self, b: bool) {
        self.my_k_part = b;
    }

    /// OCCT GenerateParticularCase() (lxx L38-41).
    pub fn generate_particular_case(&self) -> bool {
        self.my_k_part
    }

    /// OCCT ErrorOnSurf() (lxx L44-48).
    pub fn error_on_surf(&self) -> f64 {
        self.my_error
    }

    /// OCCT IsDone() (lxx L51-55).
    pub fn is_done(&self) -> bool {
        self.my_status == PipeError::PipeOk
    }

    /// OCCT GetStatus() (hxx L290).
    pub fn get_status(&self) -> PipeError {
        self.my_status
    }
}

impl Default for Pipe {
    fn default() -> Self {
        Self::new()
    }
}
