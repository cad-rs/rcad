// OCCT HLRBRep_Intersector (TKHLR) — the edge-intersection layer of the
// HLR chain: 2D intersections of the projections of 3D curves over the
// HLRBRep_CInter engine, and the 3D line/surface intersection over the
// HLRBRep_InterCSurf engine.
//
// HLRBRep_Intersector.hxx L26-98 + .cxx L86-646.  The trailing commented
// block of the OCCT file (the "sauvegarde de l etat du 23 janvier 98") is
// an OCCT comment and is not translated; the PERF counters are an OCCT
// #ifdef block and are not translated.  The OCCT `HLRBRep_Surface*` /
// `HLRBRep_Curve*` raw-pointer members keep the same raw-pointer shape
// (the HLRBRep_Data object graph owns the pointees).

use glam::DVec2;
use rcad_kernel::geom::Line3;
use rcad_kernel::math::bnd::BoundSortBox;

use crate::geomalgo::int_curve_surface::IntersectionPoint as CSIntersectionPoint;
use crate::geomalgo::int_curve_surface::IntersectionSegment as CSIntersectionSegment;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};
use crate::geomalgo::int_patch::elclib;
use crate::geomalgo::int_res2d::{
    Domain as Res2dDomain, IntersectionPoint, IntersectionSegment, Position, Transition,
};

use super::c_inter::CInter;
use super::curve::Curve;
use super::edge_data::EdgeData;
use super::inter_csurf::InterCSurf;
use super::surface::Surface;
use super::surface_tool::{LineTool, SurfaceTool};

/// OCCT Precision::Confusion() (Precision.hxx).
const PRECISION_CONFUSION: f64 = 1.0e-7;
/// OCCT Standard_Real RealLast() (Standard_Real.hxx).
const REAL_LAST: f64 = f64::MAX;

/// OCCT `myC1->D1(u, pa, va)` — the 2D D1 of the projected curve.
fn curve_d1_2d(c: &Curve<'_>, u: f64) -> (DVec2, DVec2) {
    let (mut p, mut v) = (DVec2::ZERO, DVec2::ZERO);
    c.d1_2d(u, &mut p, &mut v);
    (p, v)
}

/// OCCT HLRBRep_Intersector.
pub struct Intersector<'a> {
    my_single_point: IntersectionPoint,
    my_type_perform: i32,
    my_intersector: CInter<'a>,
    my_cs_intersector: InterCSurf,
    my_surface: Option<*const Surface<'a>>,
    my_polyhedron: Option<Box<ThePolyhedronOfHInter>>,
}

impl<'a> Intersector<'a> {
    /// OCCT HLRBRep_Intersector() (cxx L86-121): the minimal HLR polygon
    /// sample count (aMinNbHLRSamples = 4).
    pub fn new() -> Self {
        const A_MIN_NB_HLR_SAMPLES: i32 = 4;
        let mut r = Intersector {
            my_single_point: IntersectionPoint::empty(),
            my_type_perform: 0,
            my_intersector: CInter::new(),
            my_cs_intersector: InterCSurf::new(),
            my_surface: None,
            my_polyhedron: None,
        };
        r.my_intersector.set_min_nb_samples(A_MIN_NB_HLR_SAMPLES);
        r
    }

    /// OCCT Perform(theEdge1, theDa1, theDb1) (cxx L124-168) — the auto
    /// intersection of an edge.  The edge domain is cut at start with
    /// da1*(b-a) and at end with db1*(b-a).
    pub fn perform_auto(&mut self, the_edge1: &mut EdgeData<'a>, the_da1: f64, the_db1: f64) {
        self.my_type_perform = 1;

        let (mut a, ta, mut b, tb) = the_edge1.status_ref().bounds();
        let my_c1 = the_edge1.geometry();
        let d = b - a;
        if the_da1 != 0.0 {
            a += d * the_da1;
        }
        if the_db1 != 0.0 {
            b -= d * the_db1;
        }
        let pa = my_c1.value(a);
        let pb = my_c1.value(b);
        let a = my_c1.parameter_2d(a);
        let b = my_c1.parameter_2d(b);
        let mut d1 = Res2dDomain::infinite();
        d1.set_values_bounded(pa, a, ta as f64, pb, b, tb as f64);

        // modified by jgv, 18.04.2016 for OCC27341: the edge tolerance is
        // not used — Precision::Confusion().
        let tol = PRECISION_CONFUSION;

        self.my_intersector.perform_cd(my_c1, &d1, tol, tol);
    }

    /// OCCT Perform(theNA, theEdge1, theDa1, theDb1, theNB, theEdge2,
    /// theDa2, theDb2, theNoBound) (cxx L171-497) — the intersection between
    /// two edges, with the decalage walk-away loop on the common bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        _the_na: i32,
        the_edge1: &EdgeData<'a>,
        the_da1: f64,
        the_db1: f64,
        _the_nb: i32,
        the_edge2: &EdgeData<'a>,
        the_da2: f64,
        the_db2: f64,
        the_en_bout: bool,
    ) {
        let my_c1 = the_edge1.geometry();
        let my_c2 = the_edge2.geometry();

        self.my_type_perform = 1;

        // modified by jgv, 18.04.2016 for OCC27341: Precision::Confusion().
        let tol1 = PRECISION_CONFUSION;
        let tol2 = PRECISION_CONFUSION;
        let tol = if tol1 > tol2 { tol1 } else { tol2 };

        let mut a_decalagea1 = 100.0f64;
        let mut a_decalagea2 = 100.0f64;
        let mut a_decalageb1 = 100.0f64;
        let mut a_decalageb2 = 100.0f64;

        loop {
            let mut a_pas_bon = false;

            let (mut a1, ta1, mut b1, tb1) = the_edge1.status_ref().bounds(); // -- Parametres 3d
            let mut ta = ta1;
            let mut tb = tb1;
            let mut mtol = tol;
            if mtol < ta as f64 {
                mtol = ta as f64;
            }
            if mtol < tb as f64 {
                mtol = tb as f64;
            }
            let _ = mtol;
            let d = b1 - a1;

            let mut pdist = tol;
            if pdist < 0.0000001 {
                pdist = 0.0000001;
            }

            if the_da1 != 0.0 {
                // -- a = a + d * theDa1;
                let (_, va1) = curve_d1_2d(my_c1, a1);
                let qwe = va1.length();
                if qwe > 1e-12 {
                    let dd = pdist * a_decalagea1 / qwe;
                    if dd < d * 0.4 {
                        a1 += dd;
                    } else {
                        a1 += d * the_da1;
                        a_decalagea1 = -1.0;
                    }
                } else {
                    a1 += d * the_da1;
                    a_decalagea1 = -1.0;
                }
            }

            if the_db1 != 0.0 {
                // -- b = b - d * theDb1;
                let (_, vb1) = curve_d1_2d(my_c1, b1);
                let qwe = vb1.length();
                if qwe > 1e-12 {
                    let dd = pdist * a_decalageb1 / qwe;
                    if dd < d * 0.4 {
                        b1 -= dd;
                    } else {
                        b1 -= d * the_db1;
                        a_decalageb1 = -1.0;
                    }
                } else {
                    b1 -= d * the_db1;
                    a_decalageb1 = -1.0;
                }
            }

            //    if(theEnBout) {  //-- ******** (OCCT comment)
            //      double d=b1-a1;
            //      a1+=d*0.45;
            //      b1-=d*0.45;
            //    }

            let pa1 = my_c1.value(a1);
            let pb1 = my_c1.value(b1);

            let a1 = my_c1.parameter_2d(a1);
            let b1 = my_c1.parameter_2d(b1);

            if the_en_bout {
                ta = -1.0;
                tb = -1.0;
            }

            if ta as f64 > tol {
                ta = tol as f32;
            }
            if tb as f64 > tol {
                tb = tol as f32;
            }

            let mut d1 = Res2dDomain::infinite();
            d1.set_values_bounded(pa1, a1, ta as f64, pb1, b1, tb as f64);

            let (mut a2, ta2, mut b2, tb2) = the_edge2.status_ref().bounds();
            ta = ta2;
            tb = tb2;
            let mut mtol = tol;
            if mtol < ta as f64 {
                mtol = ta as f64;
            }
            if mtol < tb as f64 {
                mtol = tb as f64;
            }
            let _ = mtol;

            let d = b2 - a2;

            if the_da2 != 0.0 {
                // -- a = a + d * theDa2;
                let (_, va2) = curve_d1_2d(my_c2, a2);
                let qwe = va2.length();
                if qwe > 1e-12 {
                    let dd = pdist * a_decalagea2 / qwe;
                    if dd < d * 0.4 {
                        a2 += dd;
                    } else {
                        a2 += d * the_da2;
                        a_decalagea2 = -1.0;
                    }
                } else {
                    a2 += d * the_da2;
                    a_decalagea2 = -1.0;
                }
            }

            if the_db2 != 0.0 {
                // -- b = b - d * theDb2;
                let (_, vb2) = curve_d1_2d(my_c2, b2);
                let qwe = vb2.length();
                if qwe > 1e-12 {
                    let dd = pdist * a_decalageb2 / qwe;
                    if dd < d * 0.4 {
                        b2 -= dd;
                    } else {
                        b2 -= d * the_db2;
                        a_decalageb2 = -1.0;
                    }
                } else {
                    b2 -= d * the_db2;
                    a_decalageb2 = -1.0;
                }
            }

            //    if(theEnBout) { //-- ******** (OCCT comment)
            //      double d=b2-a2;
            //      a2+=d*0.45;
            //      b2-=d*0.45;
            //    }

            let pa2 = my_c2.value(a2);
            let pb2 = my_c2.value(b2);

            let a2 = my_c2.parameter_2d(a2);
            let b2 = my_c2.parameter_2d(b2);

            if the_en_bout {
                ta = -1.0;
                tb = -1.0;
            }

            if ta as f64 > tol {
                ta = tol as f32;
            }
            if tb as f64 > tol {
                tb = tol as f32;
            }

            let mut d2 = Res2dDomain::infinite();
            d2.set_values_bounded(pa2, a2, ta as f64, pb2, b2, tb as f64);

            if the_en_bout {
                let a1a2 = if the_da1 != 0.0 || the_da2 != 0.0 {
                    pa1.distance(pa2)
                } else {
                    REAL_LAST
                };
                let a1b2 = if the_da1 != 0.0 || the_db2 != 0.0 {
                    pa1.distance(pb2)
                } else {
                    REAL_LAST
                };
                let b1a2 = if the_db1 != 0.0 || the_da2 != 0.0 {
                    pb1.distance(pa2)
                } else {
                    REAL_LAST
                };
                let b1b2 = if the_db1 != 0.0 || the_db2 != 0.0 {
                    pb1.distance(pb2)
                } else {
                    REAL_LAST
                };

                let mut cote = 1usize;
                let mut mindist = a1a2; // -- cas 1
                if mindist > a1b2 {
                    mindist = a1b2;
                    cote = 2;
                }
                if mindist > b1a2 {
                    mindist = b1a2;
                    cote = 3;
                }
                if mindist > b1b2 {
                    mindist = b1b2;
                    cote = 4;
                }

                if mindist < tol * 1000.0 {
                    a_pas_bon = true;
                    match cote {
                        1 => {
                            a_decalagea1 *= 2.0;
                            a_decalagea2 *= 2.0;
                        }
                        2 => {
                            a_decalagea1 *= 2.0;
                            a_decalageb2 *= 2.0;
                        }
                        3 => {
                            a_decalageb1 *= 2.0;
                            a_decalagea2 *= 2.0;
                        }
                        _ => {
                            a_decalageb1 *= 2.0;
                            a_decalageb2 *= 2.0;
                        }
                    }
                    if a_decalagea1 < 0.0
                        || a_decalagea2 < 0.0
                        || a_decalageb1 < 0.0
                        || a_decalageb2 <= 0.0
                    {
                        a_pas_bon = false;
                    }
                }
            }
            if !a_pas_bon {
                self.my_intersector
                    .perform_cd_cd(my_c1, &d1, my_c2, &d2, tol, tol);
                break;
            }
        }
    }

    /// OCCT SimulateOnePoint(theEdge1, theU, theEdge2, theV) (cxx L500-537)
    /// — a single IntersectionPoint (U on theEdge1) (V on theEdge2); the
    /// point is the middle on both curves.
    pub fn simulate_one_point(
        &mut self,
        the_edge1: &EdgeData<'a>,
        the_u: f64,
        the_edge2: &EdgeData<'a>,
        the_v: f64,
    ) {
        let my_c1 = the_edge1.geometry();
        let my_c2 = the_edge2.geometry();

        let u3 = my_c1.parameter_3d(the_u);
        let v3 = my_c2.parameter_3d(the_v);
        let (p13, mut t13) = curve_d1_2d(my_c1, u3);
        let (p23, mut t23) = curve_d1_2d(my_c2, v3);

        let pos1 = Position::Middle;
        let pos2 = Position::Middle;
        let mut tr1 = Transition::empty();
        let mut tr2 = Transition::empty();
        crate::geomalgo::int_imp_par_gen::determine_transition_in_out(
            pos1, &mut t13, &mut tr1, pos2, &mut t23, &mut tr2, 0.0,
        );
        self.my_type_perform = 0;
        self.my_single_point
            .set_values(p13, the_u, the_v, tr1, tr2, false);
    }

    /// OCCT Load(theSurface) (cxx L540-546).
    pub fn load(&mut self, the_surface: *const Surface<'a>) {
        self.my_surface = Some(the_surface);
        self.my_polyhedron = None;
    }

    /// OCCT Perform(theL, theP) (cxx L549-727) — the 3D line/surface
    /// intersection.
    pub fn perform_line(&mut self, l: &Line3, p: f64) {
        self.my_type_perform = 2;
        let my_surface = unsafe {
            &*self
                .my_surface
                .expect("HLRBRep_Intersector::Load not called")
        };
        use crate::geomalgo::int_curve_surface::HSurfaceTool;
        use crate::geomalgo::int_patch::GeomAbsSurfaceType;
        let typ = <SurfaceTool<'_> as HSurfaceTool>::get_type(my_surface);
        match typ {
            GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus => {
                self.my_cs_intersector.perform(l, my_surface);
            }
            _ => {
                if self.my_polyhedron.is_none() {
                    let u1 = <SurfaceTool<'_> as HSurfaceTool>::first_u_parameter(my_surface);
                    let v1 = <SurfaceTool<'_> as HSurfaceTool>::first_v_parameter(my_surface);
                    let u2 = <SurfaceTool<'_> as HSurfaceTool>::last_u_parameter(my_surface);
                    let v2 = <SurfaceTool<'_> as HSurfaceTool>::last_v_parameter(my_surface);
                    let nbsu =
                        <SurfaceTool<'_> as HSurfaceTool>::nb_samples_u(my_surface, u1, u2);
                    let nbsv =
                        <SurfaceTool<'_> as HSurfaceTool>::nb_samples_v(my_surface, v1, v2);
                    self.my_polyhedron = Some(Box::new(ThePolyhedronOfHInter::new_tool::<
                        Surface<'a>,
                        SurfaceTool<'a>,
                    >(
                        my_surface, nbsu, nbsv, u1, v1, u2, v2
                    )));
                }
                let my_polyhedron = self.my_polyhedron.as_ref().unwrap();
                let (x0, y0, z0, x1, y1, z1) = my_polyhedron
                    .bounding()
                    .get()
                    .expect("ThePolyhedronOfInterCSurf::Bounding");
                //-- On va rejeter tous les points de parametres > P
                // (OCCT uses lowercase `p` as the corner temp and uppercase
                // `P` for the parameter — renamed `pp` here to avoid
                // shadowing the wLim parameter.)
                let mut pp;
                pp = elclib::line_parameter(l, glam::DVec3::new(x0, y0, z0));
                let mut pmin = pp;
                let mut pmax = pp;
                pp = elclib::line_parameter(l, glam::DVec3::new(x0, y0, z1));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pp = elclib::line_parameter(l, glam::DVec3::new(x1, y0, z0));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pp = elclib::line_parameter(l, glam::DVec3::new(x1, y0, z1));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pp = elclib::line_parameter(l, glam::DVec3::new(x0, y1, z0));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pp = elclib::line_parameter(l, glam::DVec3::new(x0, y1, z1));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pp = elclib::line_parameter(l, glam::DVec3::new(x1, y1, z0));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pp = elclib::line_parameter(l, glam::DVec3::new(x1, y1, z1));
                if pmin > pp {
                    pmin = pp;
                }
                if pmax < pp {
                    pmax = pp;
                }
                pmin -= 0.000001;
                pmax += 0.000001;

                if pmin > p {
                    //-- on va rejeter avec les boites
                    pmin = pmax + 1.0;
                    pmax += 2.0;
                } else if pmax > p {
                    pmax = p + 0.0000001;
                }
                let polygon =
                    ThePolygonOfHInter::new_tool_range::<Line3, LineTool>(l, pmin, pmax, 3);
                self.my_cs_intersector.perform_polygon_polyhedron(
                    l,
                    &polygon,
                    my_surface,
                    my_polyhedron,
                );
            }
        }
    }

    /// OCCT IsDone (cxx L730-745).
    pub fn is_done(&self) -> bool {
        if self.my_type_perform == 1 {
            self.my_intersector.is_done()
        } else if self.my_type_perform == 2 {
            self.my_cs_intersector.is_done()
        } else {
            true
        }
    }

    /// OCCT NbPoints (cxx L748-767).
    pub fn nb_points(&self) -> usize {
        if self.my_type_perform == 43 {
            return 0;
        }
        if self.my_type_perform == 1 {
            self.my_intersector.nb_points()
        } else if self.my_type_perform == 2 {
            self.my_cs_intersector.nb_points()
        } else {
            1
        }
    }

    /// OCCT Point(N) (cxx L770-780) — 1-based.
    pub fn point(&self, n: usize) -> &IntersectionPoint {
        if self.my_type_perform == 0 {
            &self.my_single_point
        } else {
            self.my_intersector.point(n)
        }
    }

    /// OCCT CSPoint(N) (cxx L783-786).
    pub fn cs_point(&self, n: usize) -> &CSIntersectionPoint {
        self.my_cs_intersector.point(n)
    }

    /// OCCT NbSegments (cxx L789-804).
    pub fn nb_segments(&self) -> usize {
        if self.my_type_perform == 1 {
            self.my_intersector.nb_segments()
        } else if self.my_type_perform == 2 {
            self.my_cs_intersector.nb_segments()
        } else {
            0
        }
    }

    /// OCCT Segment(N) (cxx L807-810).
    pub fn segment(&self, n: usize) -> &IntersectionSegment {
        self.my_intersector.segment(n)
    }

    /// OCCT CSSegment(N) (cxx L813-816).
    pub fn cs_segment(&self, n: usize) -> &CSIntersectionSegment {
        self.my_cs_intersector.segment(n)
    }

    /// OCCT Destroy (cxx L819-823).
    pub fn destroy(&mut self) {
        self.my_polyhedron = None;
    }
}

impl Default for Intersector<'_> {
    fn default() -> Self {
        Intersector::new()
    }
}

// The BoundSortBox parity reference for the polyhedron-rejection path (the
// OCCT box culling stays inside the Intf engines).
#[allow(unused)]
fn _bsb_parity() -> BoundSortBox {
    BoundSortBox::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::base::proj_lib::CurveType;
    use rcad_kernel::geom::{Circle3, Ellipse3, Point3, Surface3, Vec3};
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::topods::{BRepBuilder, Shape};

    use crate::hlr::algo::projector::Projector;
    use crate::hlr::brep::b_curve_tool::CurveView;

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double).
    pub(crate) struct SegEdge {
        pub origin: Point3,
        pub dir: Vec3,
        pub len: f64,
    }

    impl CurveView for SegEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.len
        }
        fn d0(&self, u: f64) -> Point3 {
            self.origin + self.dir * u
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (self.origin + self.dir * u, self.dir)
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (self.origin + self.dir * u, self.dir, Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> rcad_kernel::geom::Line3 {
            rcad_kernel::geom::Line3::new(self.origin, self.dir)
        }
        fn circle(&self) -> Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn degree(&self) -> i32 {
            0
        }
        fn nb_poles(&self) -> i32 {
            0
        }
        fn nb_knots(&self) -> i32 {
            0
        }
        fn is_closed(&self) -> bool {
            false
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            0.0
        }
        fn resolution(&self, r3d: f64) -> f64 {
            r3d
        }
        fn parameter_3d(&self, p2d: f64) -> f64 {
            p2d
        }
        fn poles(&self) -> Vec<Point3> {
            Vec::new()
        }
    }

    /// The OCCT edge-data setup: Set(...) with the loaded projected curve
    /// and the full-parameter domain (the Data's Set is the Stage 3f
    /// boundary; the geometry/tolerance split is documented on the rcad
    /// EdgeData::set).
    fn line_edge_data(proj: &'static Projector, edge: &'static SegEdge) -> EdgeData<'static> {
        let mut c = Curve::new();
        c.projector(proj);
        c.load(edge);
        c.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);
        let mut e = EdgeData::new();
        e.set(true, false, c, 0.0, 0, 0, false, false, false, false, 0.0, 0.0, 1.0, 0.0);
        e
    }

    /// OCCT anchor: Perform(edge1, edge2) over two straight edges whose top
    /// view images cross at (0.5, 0) — the intersector builds the 2D
    /// domains from the edge status bounds (Parameter2d conversions, cxx
    /// L171-497) and reports the single crossing point.
    #[test]
    fn intersector_two_line_edges_cross() {
        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(|| {
            Projector::from_ax2(&Ax2::new(
                glam::DVec3::ZERO,
                glam::DVec3::new(0.0, 0.0, 1.0),
                glam::DVec3::new(1.0, 0.0, 0.0),
            ))
        });
        static E1: std::sync::OnceLock<SegEdge> = std::sync::OnceLock::new();
        let e1: &'static SegEdge = E1.get_or_init(|| SegEdge {
            origin: Point3::ZERO,
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 1.0,
        });
        static E2: std::sync::OnceLock<SegEdge> = std::sync::OnceLock::new();
        let e2: &'static SegEdge = E2.get_or_init(|| SegEdge {
            origin: Point3::new(0.5, -0.5, 0.0),
            dir: Vec3::new(0.0, 1.0, 0.0),
            len: 1.0,
        });

        let ed1 = line_edge_data(proj, e1);
        let ed2 = line_edge_data(proj, e2);

        let mut inter = Intersector::new();
        inter.perform(0, &ed1, 0.0, 0.0, 0, &ed2, 0.0, 0.0, false);

        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 1, "nb={}", inter.nb_points());
        let p = inter.point(1);
        assert!((p.value().x - 0.5).abs() < 1e-9, "x={}", p.value().x);
        assert!(p.value().y.abs() < 1e-9, "y={}", p.value().y);
        assert!((p.param_on_first() - 0.5).abs() < 1e-9);
        assert!((p.param_on_second() - 0.5).abs() < 1e-9);
        assert_eq!(inter.nb_segments(), 0);
    }

    /// OCCT anchor: SimulateOnePoint (cxx L500-537) — the single point at
    /// the middle of both curves with the IN/OUT transitions from the
    /// tangents.
    #[test]
    fn intersector_simulate_one_point() {
        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(|| {
            Projector::from_ax2(&Ax2::new(
                glam::DVec3::ZERO,
                glam::DVec3::new(0.0, 0.0, 1.0),
                glam::DVec3::new(1.0, 0.0, 0.0),
            ))
        });
        static E1: std::sync::OnceLock<SegEdge> = std::sync::OnceLock::new();
        let e1: &'static SegEdge = E1.get_or_init(|| SegEdge {
            origin: Point3::ZERO,
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 1.0,
        });
        static E2: std::sync::OnceLock<SegEdge> = std::sync::OnceLock::new();
        let e2: &'static SegEdge = E2.get_or_init(|| SegEdge {
            origin: Point3::new(0.5, -0.5, 0.0),
            dir: Vec3::new(0.0, 1.0, 0.0),
            len: 1.0,
        });

        let ed1 = line_edge_data(proj, e1);
        let ed2 = line_edge_data(proj, e2);

        let mut inter = Intersector::new();
        inter.simulate_one_point(&ed1, 0.5, &ed2, 0.5);

        assert_eq!(inter.nb_points(), 1);
        let p = inter.point(1);
        assert!((p.param_on_first() - 0.5).abs() < 1e-12);
        assert!((p.param_on_second() - 0.5).abs() < 1e-12);
        // The tangents (+X, +Y) cross positively: Out on the first, In on
        // the second (IntImpParGen::DetermineTransition).
        assert_eq!(
            p.transition_of_first().transition_type(),
            crate::geomalgo::int_res2d::TypeTrans::Out
        );
        assert_eq!(
            p.transition_of_second().transition_type(),
            crate::geomalgo::int_res2d::TypeTrans::In
        );
    }

    /// OCCT anchor: Load + Perform(gp_Lin, P) on a plane face — the CS
    /// branch (cxx L552-557) routes through InterCSurf; the ray
    /// (0, 0.5, 2) + (-Z) hits z = 0.5 at w = 1.5.
    #[test]
    fn intersector_line_vs_plane_cs() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, glam::DVec3::new(0.0, 0.0, 0.5), 1e-7);
        let v2 = b.add_vertex(&mut brep, glam::DVec3::new(2.0, 0.0, 0.5), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(glam::DVec3::new(0.0, 0.0, 0.5), glam::DVec3::new(1.0, 0.0, 0.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face: Shape = brep.add_tface(
            Some(Surface3::Plane(rcad_kernel::geom::Plane {
                origin: glam::DVec3::new(0.0, 0.0, 0.5),
                normal: glam::DVec3::new(0.0, 0.0, 1.0),
                u_dir: glam::DVec3::new(1.0, 0.0, 0.0),
                v_dir: glam::DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 2.0]),
            Vec::new(),
            true,
        );

        let mut hsurf = Surface::new();
        hsurf.load(&brep, &face);

        let mut inter = Intersector::new();
        inter.load(&hsurf);
        let line = Line3::new(glam::DVec3::new(0.0, 0.5, 2.0), glam::DVec3::new(0.0, 0.0, -1.0));
        inter.perform_line(&line, 100.0);

        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 1, "nb={}", inter.nb_points());
        let pt = inter.cs_point(1);
        assert!((pt.pnt().z - 0.5).abs() < 1e-9);
        assert!((pt.w() - 1.5).abs() < 1e-9);
    }
}
