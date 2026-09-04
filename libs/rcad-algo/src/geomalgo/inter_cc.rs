// OCCT Geom2dAPI_InterCurveCurve (TKGeomAlgo/Geom2dAPI) — 1:1 Rust
// translation of Geom2dAPI_InterCurveCurve.hxx (class + members),
// .cxx L31-210 (all methods) and .lxx (Intersector()).
//
// This class implements methods for computing
// - the intersections between two 2D curves,
// - the self-intersections of a 2D curve.

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, TrimmedCurve2};

use super::geom2d_int::Curve2dAdaptor;
use super::int_curve_curve_gen::IntCurveCurveGen;
use super::int_res2d::IntersectionPoint;

/// OCCT Geom2dAPI_InterCurveCurve.
#[derive(Debug, Clone)]
pub struct InterCurveCurve {
    /// OCCT myIsDone.
    my_is_done: bool,
    /// OCCT myCurve1.
    my_curve1: Option<Curve2d>,
    /// OCCT myCurve2.
    my_curve2: Option<Curve2d>,
    /// OCCT myIntersector.
    my_intersector: IntCurveCurveGen,
}

impl InterCurveCurve {
    /// OCCT Geom2dAPI_InterCurveCurve() (cxx L31-34) — create an empty
    /// intersector. Use the function Init for further initialization of the
    /// intersection algorithm by curves or curve.
    pub fn new() -> Self {
        InterCurveCurve {
            my_is_done: false,
            my_curve1: None,
            my_curve2: None,
            my_intersector: IntCurveCurveGen::new(),
        }
    }

    /// OCCT Geom2dAPI_InterCurveCurve(C1, C2, Tol = 1.0e-6) (cxx L38-43) —
    /// creates an object and computes the intersections between the curves
    /// C1 and C2. Standard_NullObject is raised when a curve is None.
    pub fn new_curves(c1: Option<&Curve2d>, c2: Option<&Curve2d>, tol: f64) -> Self {
        let mut r = InterCurveCurve::new();
        r.init(c1, c2, tol);
        r
    }

    /// OCCT Geom2dAPI_InterCurveCurve(C1, Tol = 1.0e-6) (cxx L47-51) —
    /// creates an object and computes the self-intersections of the curve C1.
    pub fn new_curve(c1: Option<&Curve2d>, tol: f64) -> Self {
        let mut r = InterCurveCurve::new();
        r.init_curve(c1, tol);
        r
    }

    /// OCCT Init(C1, C2, Tol = 1.0e-6) (cxx L55-68).
    ///
    /// Standard_NullObject_Raise_if(C1.IsNull()) /
    /// Standard_NullObject_Raise_if(C2.IsNull()) map to a panic on None.
    pub fn init(&mut self, c1: Option<&Curve2d>, c2: Option<&Curve2d>, tol: f64) {
        let c1 = c1.expect("Geom2dAPI_InterCurveCurve::Init - C1 is null");
        let c2 = c2.expect("Geom2dAPI_InterCurveCurve::Init - C2 is null");
        // myCurve1 = occ::down_cast<Geom2d_Curve>(C1->Copy());
        self.my_curve1 = Some(c1.clone());
        // myCurve2 = occ::down_cast<Geom2d_Curve>(C2->Copy());
        self.my_curve2 = Some(c2.clone());

        // Geom2dAdaptor_Curve AC1(C1); Geom2dAdaptor_Curve AC2(C2);
        // myIntersector = Geom2dInt_GInter(AC1, AC2, Tol, Tol);
        self.my_intersector = IntCurveCurveGen::new_cc(c1, c2, tol, tol);
        self.my_is_done = self.my_intersector.is_done();
    }

    /// OCCT Init(C1, Tol = 1.0e-6) (cxx L72-81).
    ///
    /// Standard_NullObject_Raise_if(C1.IsNull()) maps to a panic on None.
    pub fn init_curve(&mut self, c1: Option<&Curve2d>, tol: f64) {
        let c1 = c1.expect("Geom2dAPI_InterCurveCurve::Init - C1 is null");
        self.my_curve1 = Some(c1.clone());
        self.my_curve2 = None;

        // Geom2dAdaptor_Curve AC1(C1);
        // myIntersector = Geom2dInt_GInter(AC1, Tol, Tol);
        self.my_intersector = IntCurveCurveGen::new_c(c1, tol, tol);
        self.my_is_done = self.my_intersector.is_done();
    }

    /// OCCT NbPoints() (cxx L85-95) — the number of intersection-points in
    /// case of cross intersections; returns 0 if no intersections were found.
    pub fn nb_points(&self) -> usize {
        if self.my_is_done {
            self.my_intersector.nb_points()
        } else {
            0
        }
    }

    /// OCCT Point(Index) (cxx L99-104) — the intersection point of index
    /// Index (1-based). Standard_OutOfRange if index is not in the range
    /// [1, NbPoints].
    pub fn point(&self, index: usize) -> DVec2 {
        assert!(
            index >= 1 && index <= self.nb_points(),
            "Geom2dAPI_InterCurveCurve::Point"
        );
        self.my_intersector.point(index).value()
    }

    /// OCCT NbSegments() (cxx L108-118) — the number of tangential
    /// intersections; returns 0 if no intersections were found.
    pub fn nb_segments(&self) -> usize {
        if self.my_is_done {
            self.my_intersector.nb_segments()
        } else {
            0
        }
    }

    /// OCCT Segment(theIndex, theCurve1, theCurve2) (cxx L145-210) — use this
    /// syntax only to get solutions of tangential intersection between two
    /// curves. Output values theCurve1 and theCurve2 are the intersection
    /// segments on the first curve and on the second curve accordingly.
    ///
    /// Standard_OutOfRange if theIndex is not in [1, NbSegments];
    /// Standard_NullObject if the algorithm is initialized for the computing
    /// of self-intersections on a curve.
    pub fn segment(
        &self,
        the_index: usize,
        the_curve1: &mut Option<Curve2d>,
        the_curve2: &mut Option<Curve2d>,
    ) {
        assert!(
            the_index >= 1 && the_index <= self.nb_segments(),
            "Geom2dAPI_InterCurveCurve::Segment"
        );

        // Standard_NullObject_Raise_if(myCurve1.IsNull(), "...").
        let my_curve1 = self
            .my_curve1
            .as_ref()
            .expect("Geom2dAPI_InterCurveCurve::Segment");

        let mut a_u1;
        let mut a_u2;
        let mut a_v1;
        let mut a_v2;
        a_u1 = my_curve1.first_parameter();
        a_u2 = my_curve1.last_parameter();
        match self.my_curve2.as_ref() {
            None => {
                a_v1 = a_u1;
                a_v2 = a_u2;
            }
            Some(my_curve2) => {
                a_v1 = my_curve2.first_parameter();
                a_v2 = my_curve2.last_parameter();
            }
        }

        let a_seg = self.my_intersector.segment(the_index);
        let is_opposite = a_seg.is_opposite();

        if a_seg.has_first_point() {
            let an_ipf: &IntersectionPoint = a_seg.first_point();
            a_u1 = an_ipf.param_on_first();

            if is_opposite {
                a_v2 = an_ipf.param_on_second();
            } else {
                a_v1 = an_ipf.param_on_second();
            }
        }

        if a_seg.has_last_point() {
            let an_ipl: &IntersectionPoint = a_seg.last_point();
            a_u2 = an_ipl.param_on_first();

            if is_opposite {
                a_v1 = an_ipl.param_on_second();
            } else {
                a_v2 = an_ipl.param_on_second();
            }
        }

        // theCurve1 = new Geom2d_TrimmedCurve(myCurve1, aU1, aU2);
        *the_curve1 = Some(Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(my_curve1.clone()),
            t_min: a_u1,
            t_max: a_u2,
        }));
        match self.my_curve2.as_ref() {
            None => {
                // theCurve2 = new Geom2d_TrimmedCurve(myCurve1, aV1, aV2);
                *the_curve2 = Some(Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(my_curve1.clone()),
                    t_min: a_v1,
                    t_max: a_v2,
                }));
            }
            Some(my_curve2) => {
                // theCurve2 = new Geom2d_TrimmedCurve(myCurve2, aV1, aV2);
                *the_curve2 = Some(Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(my_curve2.clone()),
                    t_min: a_v1,
                    t_max: a_v2,
                }));
            }
        }
    }

    /// OCCT Intersector() (lxx L19-22) — return the algorithmic object from
    /// Intersection.
    pub fn intersector(&self) -> &IntCurveCurveGen {
        &self.my_intersector
    }
}

impl Default for InterCurveCurve {
    fn default() -> Self {
        InterCurveCurve::new()
    }
}
