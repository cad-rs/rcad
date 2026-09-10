//! OCCT Geom2dHatch_Hatching (TKGeomAlgo/Geom2dHatch).
//!
//! Geom2dHatch_Hatching.hxx L38-140 + .cxx L27-332 — one hatching curve with
//! its trimmed domains and the intersection points collected on it.

use rcad_kernel::geom::{Curve2d, Line2d};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::hatch::hatch_gen::{Domain, ErrorStatus, PointOnHatching};

/// OCCT Geom2dHatch_Hatching — a hatching: its curve, trimming state, the
/// intersection points on the hatching and the resulting domains.
#[derive(Debug, Clone)]
pub struct Hatching {
    /// OCCT Geom2dAdaptor_Curve myCurve.
    my_curve: Curve2d,
    /// OCCT bool myTrimDone.
    my_trim_done: bool,
    /// OCCT bool myTrimFailed.
    my_trim_failed: bool,
    /// OCCT NCollection_Sequence<HatchGen_PointOnHatching> myPoints.
    my_points: Vec<PointOnHatching>,
    /// OCCT bool myIsDone.
    my_is_done: bool,
    /// OCCT HatchGen_ErrorStatus myStatus.
    my_status: ErrorStatus,
    /// OCCT NCollection_Sequence<HatchGen_Domain> myDomains.
    my_domains: Vec<Domain>,
}

impl Hatching {
    /// OCCT Geom2dHatch_Hatching() (cxx L27-33) — the default hatching.  OCCT
    /// leaves the curve adaptor with a null curve handle (using it is
    /// undefined); rcad initializes a neutral degenerate line so the value
    /// stays deterministic.
    pub fn empty() -> Self {
        Hatching {
            // OCCT leaves the adaptor with a null curve handle; rcad uses a
            // neutral degenerate line (same default as HatchElement::empty).
            my_curve: Curve2d::Line(Line2d::new(glam::DVec2::ZERO, glam::DVec2::X)),
            my_trim_done: false,
            my_trim_failed: false,
            my_points: Vec::new(),
            my_is_done: false,
            my_status: ErrorStatus::NoProblem,
            my_domains: Vec::new(),
        }
    }

    /// OCCT Geom2dHatch_Hatching(Curve) (cxx L37-44) — creates a hatching.
    pub fn new(curve: Curve2d) -> Self {
        Hatching {
            my_curve: curve,
            my_trim_done: false,
            my_trim_failed: false,
            my_points: Vec::new(),
            my_is_done: false,
            my_status: ErrorStatus::NoProblem,
            my_domains: Vec::new(),
        }
    }

    /// OCCT Curve (cxx L51-54) — the curve associated to the hatching.
    pub fn curve(&self) -> &Curve2d {
        &self.my_curve
    }

    /// OCCT ChangeCurve (cxx L61-64).
    pub fn change_curve(&mut self) -> &mut Curve2d {
        &mut self.my_curve
    }

    /// OCCT TrimDone(Flag) setter (cxx L72-75).
    pub fn set_trim_done(&mut self, flag: bool) {
        self.my_trim_done = flag;
    }

    /// OCCT TrimDone() getter (cxx L82-85).
    pub fn trim_done(&self) -> bool {
        self.my_trim_done
    }

    /// OCCT TrimFailed(Flag) setter (cxx L92-100) — a trimming failure also
    /// sets the status.
    pub fn set_trim_failed(&mut self, flag: bool) {
        self.my_trim_failed = flag;
        if self.my_trim_failed {
            self.my_status = ErrorStatus::TrimFailure;
        }
    }

    /// OCCT TrimFailed() getter (cxx L107-110).
    pub fn trim_failed(&self) -> bool {
        self.my_trim_failed
    }

    /// OCCT IsDone(Flag) setter (cxx L118-121).
    pub fn set_is_done(&mut self, flag: bool) {
        self.my_is_done = flag;
    }

    /// OCCT IsDone() getter (cxx L128-131).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT Status(theStatus) setter (cxx L135-138).
    pub fn set_status(&mut self, status: ErrorStatus) {
        self.my_status = status;
    }

    /// OCCT Status() getter (cxx L142-145).
    pub fn status(&self) -> ErrorStatus {
        self.my_status
    }

    /// OCCT AddPoint (cxx L152-189) — adds an intersection point to the
    /// hatching (kept sorted by parameter; equal points are merged by adding
    /// their points on element).
    pub fn add_point(&mut self, point: &PointOnHatching, confusion: f64) {
        let nb_points = self.my_points.len();
        // for (IPntH = 1; IPntH <= NbPoints; IPntH++) { if
        // (!PntH.IsLower(Point, Confusion)) break; }
        let mut ipnth = 1usize;
        while ipnth <= nb_points {
            let pnth = &self.my_points[ipnth - 1];
            if !pnth.is_lower(point, confusion) {
                break;
            }
            ipnth += 1;
        }
        if ipnth > nb_points {
            // myPoints.Append(Point);
            self.my_points.push(point.clone());
        } else {
            // PntH.IsGreater(Point, Confusion) is evaluated before the insert
            // to split the OCCT ChangeValue borrow (same logic).
            let a_is_greater = self.my_points[ipnth - 1].is_greater(point, confusion);
            if a_is_greater {
                // myPoints.InsertBefore(IPntH, Point);
                self.my_points.insert(ipnth - 1, point.clone());
            } else {
                for ipnte in 1..=point.nb_points() {
                    let pnte = point.point(ipnte);
                    let pnth = &mut self.my_points[ipnth - 1];
                    pnth.add_point(pnte, confusion);
                }
            }
        }
        if self.my_is_done {
            self.clr_domains();
        }
    }

    /// OCCT NbPoints (cxx L196-199) — the number of intersection points on
    /// the hatching.
    pub fn nb_points(&self) -> usize {
        self.my_points.len()
    }

    /// OCCT Point (cxx L206-209) — the Index-th intersection point on the
    /// hatching (1-based; OCCT raises OutOfRange when out of bounds).
    pub fn point(&self, index: usize) -> &PointOnHatching {
        &self.my_points[index - 1]
    }

    /// OCCT ChangePoint (cxx L216-219) — the Index-th intersection point
    /// (1-based).
    pub fn change_point(&mut self, index: usize) -> &mut PointOnHatching {
        &mut self.my_points[index - 1]
    }

    /// OCCT RemPoint (cxx L226-233) — removes the Index-th intersection
    /// point of the hatching (1-based).
    pub fn rem_point(&mut self, index: usize) {
        if self.my_is_done {
            self.clr_domains();
        }
        self.my_points.remove(index - 1);
    }

    /// OCCT ClrPoints (cxx L240-254) — removes all the intersection points
    /// of the hatching.
    pub fn clr_points(&mut self) {
        if self.my_is_done {
            self.clr_domains();
        }
        for ipnth in 1..=self.my_points.len() {
            let point = &mut self.my_points[ipnth - 1];
            point.clr_points();
        }
        self.my_points.clear();
        self.my_trim_done = false;
        self.my_trim_failed = false;
    }

    /// OCCT AddDomain (cxx L261-264) — adds a domain to the hatching.
    pub fn add_domain(&mut self, domain: &Domain) {
        self.my_domains.push(domain.clone());
    }

    /// OCCT NbDomains (cxx L271-274) — the number of domains on the hatching.
    pub fn nb_domains(&self) -> usize {
        self.my_domains.len()
    }

    /// OCCT Domain (cxx L281-284) — the Index-th domain on the hatching
    /// (1-based).
    pub fn domain(&self, index: usize) -> &Domain {
        &self.my_domains[index - 1]
    }

    /// The mutable counterpart of Domain (1-based).  OCCT reaches the
    /// domains of a hatching through a const_cast when a caller has to write
    /// one back — ChFi3d_Builder_SpKP.cxx Tri L608-614:
    /// `HatchGen_Domain* Dom =
    /// ((HatchGen_Domain*)(void*)&H.Domain(iH, Ind(iSansFirst)));`.
    /// HatchGen has no ChangeDomain at all, so the accessor lives here.
    pub fn change_domain(&mut self, index: usize) -> &mut Domain {
        &mut self.my_domains[index - 1]
    }

    /// OCCT RemDomain (cxx L291-294) — removes the Index-th domain of the
    /// hatching (1-based).
    pub fn rem_domain(&mut self, index: usize) {
        self.my_domains.remove(index - 1);
    }

    /// OCCT ClrDomains (cxx L301-305) — removes all the domains of the
    /// hatching.
    pub fn clr_domains(&mut self) {
        self.my_domains.clear();
        self.my_is_done = false;
    }

    /// OCCT ClassificationPoint (cxx L311-332) — returns a point on the
    /// curve; this point will be used for the classification.
    pub fn classification_point(&self) -> glam::DVec2 {
        // OCCT double t, a, b; a = myCurve.FirstParameter(); b =
        // myCurve.LastParameter();
        let a = self.my_curve.first_parameter();
        let b = self.my_curve.last_parameter();
        let t;
        if b >= rcad_kernel::INFINITE_VALUE {
            if a <= -rcad_kernel::INFINITE_VALUE {
                t = 0.0;
            } else {
                t = a;
            }
        } else {
            t = b;
        }
        self.my_curve.value(t)
    }
}

impl Default for Hatching {
    fn default() -> Self {
        Self::empty()
    }
}
