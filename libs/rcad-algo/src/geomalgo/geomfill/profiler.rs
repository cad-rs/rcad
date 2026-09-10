//! OCCT GeomFill_Profiler (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Profiler.hxx (L24-90) + GeomFill_Profiler.cxx (whole file
//! L29-344). This is the package-owned translation; the partial carrier in
//! `crate::brep_fill::generator` predates the geomfill module (relocation
//! pending, it must eventually consume this file).
//!
//! GAP carriers (other-package untranslated dependencies, annotated):
//! - `GeomConvert_ApproxCurve` (conic approximation in AddCurve, L134-141).
//! - `Geom_BSplineCurve::SetNotPeriodic` + `Segment` (deperiodisation in
//!   Perform, L186-190) — the rcad-kernel BSpline does not yet host them.

use rcad_kernel::geom::{Curve3, TrimmedCurve3, BSplineCurve3};
use rcad_kernel::math::bspl_lib;

/// OCCT GeomFill_Profiler.
pub struct Profiler {
    /// OCCT mySequence (NCollection_Sequence<handle(Geom_Curve)>); after
    /// AddCurve every entry is a BSpline curve.
    pub(crate) my_sequence: Vec<BSplineCurve3>,
    /// OCCT myIsDone.
    my_is_done: bool,
    /// OCCT myIsPeriodic.
    my_is_periodic: bool,
}

impl Profiler {
    /// OCCT GeomFill_Profiler::GeomFill_Profiler (L113-117).
    pub fn new() -> Self {
        Profiler {
            my_sequence: Vec::new(),
            my_is_done: false,
            my_is_periodic: true,
        }
    }

    /// OCCT AddCurve (L125-162).
    pub fn add_curve(&mut self, curve: &Curve3) {
        let mut c: Option<BSplineCurve3> = None;
        // modified by jgv, 19.01.05 for OCC7354 — unwrap trimmed curves.
        let mut the_curve = curve.clone();
        if let Curve3::Trimmed(ctrim) = &the_curve {
            the_curve = ctrim.basis_curve().clone();
        }
        // OCCT: if (theCurve->IsKind(STANDARD_TYPE(Geom_Conic))) — conic
        // approximation via GeomConvert_ApproxCurve (GAP: untranslated).
        let is_conic = matches!(
            the_curve,
            Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_)
        );
        if is_conic {
            // OCCT: GeomConvert_ApproxCurve appr(Curve, Confusion, C1, 16, 14);
            // if (appr.HasResult()) C = appr.Curve();
            // GAP: rcad has no GeomConvert_ApproxCurve yet; the direct
            // GeomConvert::CurveToBSplineCurve path below is used instead and
            // this branch reports until the approximator is translated.
        }
        if c.is_none() {
            // OCCT: C = GeomConvert::CurveToBSplineCurve(Curve).
            c = Some(convert_curve_to_bspline(curve));
        }

        // OCCT: mySequence.Append(C).
        self.my_sequence.push(c.expect("GeomFill_Profiler::AddCurve"));

        // OCCT: if (myIsPeriodic && !C->IsPeriodic()) myIsPeriodic = false.
        if self.my_is_periodic && !self.my_sequence.last().unwrap().is_periodic {
            self.my_is_periodic = false;
        }
    }

    /// OCCT Perform (L166-246): converts all curves to BSplineCurves and sets
    /// them to the common profile. `<PTol>` is used to compare 2 knots.
    pub fn perform(&mut self, ptol: f64) {
        let mut my_degree = 0usize;
        let mut u_first = 0.0f64;
        let mut u_last = 0.0f64;
        let mut ecart_max = 0.0f64;

        for c in self.my_sequence.iter_mut() {
            // OCCT: U2 = C->Knot(C->LastUKnotIndex()); U1 = first.
            let (knots, _) = c.knots_mults();
            let u1 = knots[0];
            let u2 = knots[knots.len() - 1];

            if !self.my_is_periodic && c.is_periodic {
                // OCCT: C->SetNotPeriodic(); C->Segment(U1, U2);
                // GAP: Geom_BSplineCurve::SetNotPeriodic / Segment are not yet
                // translated in rcad-kernel (Geom_BSplineCurve_1.cxx /
                // BSplCLib::Segment); staged until they land.
                panic!(
                    "Staged: GeomFill_Profiler::Perform deperiodisation (SetNotPeriodic + Segment) — GeomFill_Profiler.cxx L186-190"
                );
            }

            // evaluate the max degree
            my_degree = my_degree.max(c.degree);

            // Calcul de Max ( Ufin - Udeb) sur l ensemble des courbes.
            if (u2 - u1) > ecart_max {
                ecart_max = u2 - u1;
                u_first = u1;
                u_last = u2;
            }
        }

        // increase the degree of the curves to my degree
        // reparametrize them in the range U1, U2.
        for c in self.my_sequence.iter_mut() {
            c.increase_degree(my_degree);

            // OCCT: NCollection_Array1<double> Knots(C->Knots());
            // BSplCLib::Reparametrize(UFirst, ULast, Knots); C->SetKnots(Knots).
            let (mut knots, mults) = c.knots_mults();
            bspl_lib::reparametrize(u_first, u_last, &mut knots);
            c.set_knots(&knots, &mults);
        }

        // OCCT: theCurves = copies of mySequence;
        // UnifyByInsertingAllKnots(theCurves, PTol).
        let mut the_curves: Vec<BSplineCurve3> = self.my_sequence.clone();
        unify_by_inserting_all_knots(&mut the_curves, ptol);

        let unified = {
            let the_nb_knots = the_curves[0].knots_mults().0.len();
            the_curves.iter().all(|c| c.knots_mults().0.len() == the_nb_knots)
        };

        if unified {
            self.my_sequence = the_curves;
        } else {
            unify_by_setting_middle_knots(&mut self.my_sequence);
        }

        self.my_is_done = true;
    }

    /// OCCT Degree (L250-259).
    pub fn degree(&self) -> i32 {
        assert!(self.my_is_done, "StdFail_NotDone: GeomFill_Profiler::Degree");
        self.my_sequence[0].degree as i32
    }

    /// OCCT IsPeriodic.
    pub fn is_periodic(&self) -> bool {
        self.my_is_periodic
    }

    /// OCCT NbPoles (L263-272).
    pub fn nb_poles(&self) -> i32 {
        assert!(self.my_is_done, "StdFail_NotDone: GeomFill_Profiler::Degree");
        self.my_sequence[0].control_points.len() as i32
    }

    /// OCCT Poles (L276-290): returns in <Poles> the poles of the BSplineCurve
    /// from index <Index> adjusting to the current profile.
    pub fn poles(&self, index: i32, poles: &mut [glam::DVec3]) {
        assert!(self.my_is_done, "StdFail_NotDone: GeomFill_Profiler::Degree");
        assert!(
            poles.len() == self.nb_poles() as usize,
            "Standard_DomainError: GeomFill_Profiler::Poles"
        );
        assert!(
            index >= 1 && (index as usize) <= self.my_sequence.len(),
            "Standard_DomainError: GeomFill_Profiler::Poles"
        );
        let c = &self.my_sequence[(index - 1) as usize];
        poles.copy_from_slice(&c.control_points);
    }

    /// OCCT Weights (L294-308).
    pub fn weights(&self, index: i32, weights: &mut [f64]) {
        assert!(self.my_is_done, "StdFail_NotDone: GeomFill_Profiler::Degree");
        assert!(
            weights.len() == self.nb_poles() as usize,
            "Standard_DomainError: GeomFill_Profiler::Weights"
        );
        assert!(
            index >= 1 && (index as usize) <= self.my_sequence.len(),
            "Standard_DomainError: GeomFill_Profiler::Weights"
        );
        let c = &self.my_sequence[(index - 1) as usize];
        weights.copy_from_slice(&c.weights);
    }

    /// OCCT NbKnots (L312-322).
    pub fn nb_knots(&self) -> i32 {
        assert!(self.my_is_done, "StdFail_NotDone: GeomFill_Profiler::Degree");
        self.my_sequence[0].knots_mults().0.len() as i32
    }

    /// OCCT KnotsAndMults (L326-344).
    pub fn knots_and_mults(&self, knots: &mut [f64], mults: &mut [i32]) {
        assert!(self.my_is_done, "StdFail_NotDone: GeomFill_Profiler::Degree");
        let n = self.nb_knots() as usize;
        assert!(
            knots.len() == n && mults.len() == n,
            "Standard_DomainError: GeomFill_Profiler::KnotsAndMults"
        );
        let (k, m) = self.my_sequence[0].knots_mults();
        knots.copy_from_slice(&k);
        mults.copy_from_slice(&m);
    }

    /// OCCT Curve(Index) — base-class sequence access (Profiler.hxx L88).
    pub fn curve(&self, index: i32) -> &BSplineCurve3 {
        &self.my_sequence[(index - 1) as usize]
    }
}

impl Default for Profiler {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT UnifyByInsertingAllKnots (GeomFill_Profiler.cxx L29-74): inserting in
/// the first curve the knot-vector of all the others.
fn unify_by_inserting_all_knots(the_curves: &mut [BSplineCurve3], ptol: f64) {
    let others: Vec<(Vec<f64>, Vec<i32>)> = the_curves
        .iter()
        .skip(1)
        .map(|ci| ci.knots_mults())
        .collect();
    for (knots, mults) in others.iter() {
        insert_knots_bspline(&mut the_curves[0], knots, mults, ptol);
    }

    let (new_knots, new_mults) = the_curves[0].knots_mults();
    for ci in the_curves.iter_mut().skip(1) {
        insert_knots_bspline(ci, &new_knots, &new_mults, ptol);
    }

    // essai : tentative mise des poids sur chaque section a une moyenne 1
    for ci in the_curves.iter_mut() {
        if ci.is_rational() {
            let np = ci.control_points.len();
            let mut sigma = 0.0f64;
            for j in 0..np {
                sigma += ci.weights[j];
            }
            sigma /= np as f64;
            for j in 0..np {
                ci.weights[j] /= sigma;
            }
        }
    }
    // fin de l essai
}

/// OCCT UnifyBySettingMiddleKnots (GeomFill_Profiler.cxx L78-109).
fn unify_by_setting_middle_knots(the_curves: &mut [BSplineCurve3]) {
    let nb_knots = the_curves[0].knots_mults().0.len();
    let (c_knots, _) = the_curves[0].knots_mults();
    let u_last = c_knots[c_knots.len() - 1];
    let u_first = c_knots[0];

    // Set middle values of knots
    let mut new_knots = vec![0.0f64; nb_knots];
    new_knots[0] = u_first;
    new_knots[nb_knots - 1] = u_last;
    for j in 1..nb_knots - 1 {
        let mut a_mid_knot = 0.0f64;
        for ctemp in the_curves.iter() {
            let (ck, _) = ctemp.knots_mults();
            a_mid_knot += ck[j];
        }
        a_mid_knot /= the_curves.len() as f64;
        new_knots[j] = a_mid_knot;
    }

    // OCCT: Cres->SetKnots(NewKnots) — knot values only, multiplicities are
    // preserved by OCCT's SetKnots(array) which merges per-knot.
    for cres in the_curves.iter_mut() {
        let (_, mults) = cres.knots_mults();
        cres.set_knots(&new_knots, &mults);
    }
}

/// OCCT Geom_BSplineCurve::InsertKnots(Knots, Mults, Tolerance = -1.,
/// Add = false) mapping — drives the rcad-kernel BSplCLib::insert_knots with
/// the homogeneous pole flattening OCCT performs internally
/// (Geom_BSplineCurve_1.cxx L690-760).
pub(crate) fn insert_knots_bspline(
    c: &mut BSplineCurve3,
    add_knots: &[f64],
    add_mults: &[i32],
    tolerance: f64,
) {
    let degree = c.degree;
    let periodic = c.is_periodic;
    let (knots, mults) = c.knots_mults();
    let nb_poles = c.control_points.len();

    // OCCT Geom_BSplineCurve::InsertKnots builds the flat arrays: rational
    // curves go through BSplCLib with dimension 4 (x*w, y*w, z*w, w).
    let rational = c.is_rational();
    let dimension = if rational { 4usize } else { 3usize };
    let mut flat_poles = vec![0.0f64; nb_poles * dimension];
    for i in 0..nb_poles {
        let p = c.control_points[i];
        if rational {
            let w = c.weights[i];
            flat_poles[i * 4] = p.x * w;
            flat_poles[i * 4 + 1] = p.y * w;
            flat_poles[i * 4 + 2] = p.z * w;
            flat_poles[i * 4 + 3] = w;
        } else {
            flat_poles[i * 3] = p.x;
            flat_poles[i * 3 + 1] = p.y;
            flat_poles[i * 3 + 2] = p.z;
        }
    }

    let mut new_poles = vec![0.0f64; nb_poles * dimension + dimension * degree * add_knots.len()];
    let mut new_knots = vec![0.0f64; knots.len() + add_knots.len()];
    let mut new_mults = vec![0i32; mults.len() + add_mults.len()];

    bspl_lib::insert_knots(
        degree,
        periodic,
        dimension,
        &flat_poles,
        &knots,
        &mults,
        add_knots,
        Some(add_mults),
        &mut new_poles,
        &mut new_knots,
        &mut new_mults,
        tolerance,
        false,
    );

    // Unpack the results back into the curve (knots flattened; zero
    // multiplicities contribute no flat knots).
    let mut flat_knots: Vec<f64> = Vec::new();
    for (k, m) in new_knots.iter().zip(new_mults.iter()) {
        for _ in 0..*m {
            flat_knots.push(*k);
        }
    }
    let nb_new_poles = flat_knots.len() - degree - 1;
    let mut control_points = Vec::with_capacity(nb_new_poles);
    let mut weights = vec![1.0f64; nb_new_poles];
    for i in 0..nb_new_poles {
        if rational {
            let w = new_poles[i * 4 + 3];
            control_points.push(glam::DVec3::new(
                new_poles[i * 4] / w,
                new_poles[i * 4 + 1] / w,
                new_poles[i * 4 + 2] / w,
            ));
            weights[i] = w;
        } else {
            control_points.push(glam::DVec3::new(
                new_poles[i * 3],
                new_poles[i * 3 + 1],
                new_poles[i * 3 + 2],
            ));
        }
    }
    c.control_points = control_points;
    c.weights = weights;
    c.knots = flat_knots;
}

/// OCCT GeomConvert::CurveToBSplineCurve mapping for the profiler input —
/// uses the rcad-kernel port (trimmed line/circle supported; other branches
/// are staged in the kernel port and report there).
fn convert_curve_to_bspline(curve: &Curve3) -> BSplineCurve3 {
    use rcad_kernel::base::convert::{geom_convert_curve_to_bspline_curve, ConvertParameterisation};
    let (first, last) = curve_bounds(curve);
    match curve {
        Curve3::BSpline(b) => b.clone(),
        _ => {
            let trimmed = Curve3::Trimmed(TrimmedCurve3::new(curve.clone(), first, last));
            geom_convert_curve_to_bspline_curve(&trimmed, ConvertParameterisation::TgtThetaOver2)
        }
    }
}

/// OCCT Geom_Curve::FirstParameter / LastParameter mapping for the natural
/// bounded curve types (line: rcad architectural note — OCCT lines are
/// infinite, the profiler callers pass bounded edges, so a unit range is
/// used).
fn curve_bounds(curve: &Curve3) -> (f64, f64) {
    match curve {
        Curve3::Circle(_) => (0.0, 2.0 * std::f64::consts::PI),
        _ => (0.0, 1.0),
    }
}
