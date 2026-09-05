//! OCCT Geom2dHatch_Intersector (TKGeomAlgo/Geom2dHatch).
//!
//! Geom2dHatch_Intersector.hxx L28-82 + .cxx L27-101 + .lxx L19-77 — the
//! hatch-line / element-curve intersector.  OCCT derives it from
//! Geom2dInt_GInter (the general 2D curve/curve intersection engine); rcad
//! embeds the landed Geom2dInt GInter instantiation as the `base` field and
//! delegates the IntRes2d_Intersection result accessors to it.

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Line2d};
use rcad_kernel::precision::p_intersection;
use rcad_kernel::PCONFUSION;

use crate::geomalgo::geom2d_int::{elclib2d, Curve2dAdaptor, GInter};
use crate::geomalgo::int_res2d::{Domain, IntersectionPoint, IntersectionSegment};

/// OCCT Geom2dHatch_Intersector — intersects a hatching line segment with an
/// element curve.
pub struct HatchIntersector {
    /// OCCT Geom2dInt_GInter base class subobject.
    base: GInter,
    /// OCCT double myConfusionTolerance.
    my_confusion_tolerance: f64,
    /// OCCT double myTangencyTolerance.
    my_tangency_tolerance: f64,
}

impl HatchIntersector {
    /// OCCT Geom2dHatch_Intersector() (cxx L27-31).
    pub fn new() -> Self {
        HatchIntersector {
            base: GInter::new(),
            my_confusion_tolerance: 0.0,
            my_tangency_tolerance: 0.0,
        }
    }

    /// OCCT Geom2dHatch_Intersector(Confusion, Tangency) (lxx L19-24).
    pub fn with_tolerances(confusion: f64, tangency: f64) -> Self {
        HatchIntersector {
            base: GInter::new(),
            my_confusion_tolerance: confusion,
            my_tangency_tolerance: tangency,
        }
    }

    /// OCCT ConfusionTolerance (lxx L31-34).
    pub fn confusion_tolerance(&self) -> f64 {
        self.my_confusion_tolerance
    }

    /// OCCT SetConfusionTolerance (lxx L41-44).
    pub fn set_confusion_tolerance(&mut self, confusion: f64) {
        self.my_confusion_tolerance = confusion;
    }

    /// OCCT TangencyTolerance (lxx L51-54).
    pub fn tangency_tolerance(&self) -> f64 {
        self.my_tangency_tolerance
    }

    /// OCCT SetTangencyTolerance (lxx L61-64).
    pub fn set_tangency_tolerance(&mut self, tangency: f64) {
        self.my_tangency_tolerance = tangency;
    }

    /// OCCT Perform (cxx L35-70) — intersects the 2d line segment (<L>, <P>)
    /// with the curve <C>.  The line segment is the part of <L> of parameter
    /// range [0, <P>] (<P> is positive and can be RealLast()).  <Tol> is the
    /// tolerance on the segment.  The order is relevant: the first argument
    /// is the segment, the second the edge.
    pub fn perform(&mut self, l: &Line2d, p: f64, tol: f64, c: &Curve2d) {
        // IntRes2d_Domain DL;
        let mut dl = Domain::infinite();
        if p != f64::MAX {
            // DL.SetValues(L.Location(), 0., Tol, ElCLib::Value(P, L), P, Tol);
            let a_p_end = elclib2d::line_value(l.origin, l.direction, p);
            dl.set_values_bounded(l.origin, 0.0, tol, a_p_end, p, tol);
        } else {
            // DL.SetValues(L.Location(), 0., Tol, true);
            dl.set_values_semi(l.origin, 0.0, tol, true);
        }

        // IntRes2d_Domain DE(C.Value(C.FirstParameter()), C.FirstParameter(),
        //                    Precision::PIntersection(), C.Value(C.LastParameter()),
        //                    C.LastParameter(), Precision::PIntersection());
        let de = Domain::bounded(
            c.value(c.first_parameter()),
            c.first_parameter(),
            p_intersection(),
            c.value(c.last_parameter()),
            c.last_parameter(),
            p_intersection(),
        );

        // occ::handle<Geom2d_Line> GL = new Geom2d_Line(L);
        // Geom2dAdaptor_Curve CGA(GL);
        let cga = Curve2d::Line(*l);
        // OCCT: void* ptrpoureviterlesproblemesdeconst = (void*)(&C); — a C++
        // const-cast workaround; the same `c` reference is used below.

        // Geom2dInt_GInter Inter(CGA, DL, C, DE, PConfusion(), PIntersection());
        let inter = GInter::new_cd_cd(&cga, &dl, c, &de, PCONFUSION, p_intersection());
        // this->SetValues(Inter);
        self.base = inter;
    }

    /// OCCT Intersect (lxx L73-77) — intersects the curves C1 and C2 with the
    /// intersector tolerances.  The results are retrieved by the usual
    /// methods of IntRes2d_Intersection.
    pub fn intersect(&mut self, c1: &Curve2d, c2: &Curve2d) {
        self.base
            .perform_cc(c1, c2, self.my_confusion_tolerance, self.my_tangency_tolerance);
    }

    // -- IntRes2d_Intersection (Geom2dInt_GInter base) result accessors ------

    /// OCCT IsDone() — from the Geom2dInt_GInter base.
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// OCCT NbPoints() — from the Geom2dInt_GInter base (1-based accessors).
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// OCCT Point(N) — from the Geom2dInt_GInter base.
    pub fn point(&self, n: usize) -> &IntersectionPoint {
        self.base.point(n)
    }

    /// OCCT NbSegments() — from the Geom2dInt_GInter base.
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// OCCT Segment(N) — from the Geom2dInt_GInter base.
    pub fn segment(&self, n: usize) -> &IntersectionSegment {
        self.base.segment(n)
    }

    /// OCCT LocalGeometry (cxx L74-101) — returns in <Tang>, <Norm> and <C>
    /// the tangent, normal and curvature of the edge <E> at parameter value
    /// <U> (OCCT writes the out-parameters Tang/Norm/C; rcad returns them as
    /// a tuple).
    pub fn local_geometry(&self, e: &Curve2d, u: f64) -> (DVec2, DVec2, f64) {
        // GeomLProp_CLProps2d Prop(E.Curve(), U, 2, Precision::PConfusion());
        let mut prop = rcad_kernel::base::geom_lprop::ClProps2d::with_param(e, u, 2, PCONFUSION);

        let mut c = 0.0;
        let mut tang = DVec2::new(1.0, 0.0);
        if prop.is_tangent_defined() {
            // OCCT Prop.Tangent(Tang).
            tang = prop.tangent().unwrap_or(tang);
            c = prop.curvature();
        } else {
            tang = DVec2::new(1.0, 0.0);
        }

        let mut norm = DVec2::ZERO;
        if c > PCONFUSION && c < f64::MAX {
            // OCCT Prop.Normal(Norm).
            norm = prop.normal().unwrap_or(norm);
        } else {
            // OCCT Norm.SetCoord(Tang.Y(), -Tang.X()).
            norm = DVec2::new(tang.y, -tang.x);
        }
        (tang, norm, c)
    }
}

impl Default for HatchIntersector {
    fn default() -> Self {
        Self::new()
    }
}
