// OCCT HLRBRep_EdgeFaceTool (TKHLR/HLRBRep/HLRBRep_EdgeFaceTool.hxx L1-49
// + .cxx L1-105) — the static UV-point / curvature-value helpers over a
// projected curve and surface.
//
// The OCCT `HLRBRep_CurvePtr` / `HLRBRep_SurfacePtr` parameters keep the
// raw-pointer form (the HLR methodological exception, HLRBRep_Surface::myProj
// precedent).  The two kernel contexts the OCCT body reaches through the
// pointers — `E->Curve().Edge()` (the TopoDS_Edge) and the `BRep_Tool` /
// `BRepExtrema_ExtPF` global BRep — are explicit parameters (`edge`, `brep`):
// the rcad `CurveView` trait does not carry the TopoDS_Edge, and the kernel
// boundary reference follows the HLRBRep_EdgeData::Set precedent.

use rcad_kernel::base::extrema::{ExtPS, POnSurface};
use rcad_kernel::geom::{Point3, Vec3};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo::topods::{BRepTool, Shape, State, TShape};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;
use crate::topalgo::brep_class::face_classifier::FClassifier;
use crate::topalgo::shape_source::FaceShapeSource;

use super::curve::Curve;
use super::surface::Surface;

/// OCCT BRepExtrema_ExtPF (TKTopAlgo/BRepExtrema/BRepExtrema_ExtPF.hxx
/// L28-67 + .cxx L33-107) — the point-to-face extrema with the
/// face-domain classification of the candidates.  Local translation: this
/// class is the only TKTopAlgo dependency of HLRBRep_EdgeFaceTool::UVPoint.
struct ExtPF {
    /// OCCT myExtPS — None is the OtherSurface early-out of Initialize.
    my_ext_ps: Option<ExtPS>,
    /// OCCT mySqDist — the classified square distances (the OCCT 1-based
    /// NCollection_Sequence; rcad Vec + 1-based indexing).
    my_sq_dist: Vec<f64>,
    /// OCCT myPoints — the classified extremum points.
    my_points: Vec<POnSurface>,
}

impl ExtPF {
    /// OCCT BRepExtrema_ExtPF(TheVertex, TheFace) (cxx L33-40) — the
    /// constructor runs Initialize(TheFace) + Perform(TheVertex, TheFace).
    /// `the_pnt` is OCCT `BRep_Tool::Pnt(TheVertex)`: the
    /// BRepLib_MakeVertex round-trip of UVPoint is the identity on the
    /// point, so the vertex shape is not built (the kernel-boundary
    /// precedent of HLRBRep_EdgeData::Set).
    fn new(brep: &rcad_kernel::BRep, the_pnt: Point3, the_face: &Shape) -> Self {
        // ---- Initialize(TheFace) (cxx L43-68) ----
        // mySurf.Initialize(TheFace, false);
        let my_surf = BRepAdaptorSurface::initialize_face(brep, the_face, false);
        let mut proj = ExtPF {
            my_ext_ps: None,
            my_sq_dist: Vec::new(),
            my_points: Vec::new(),
        };
        if my_surf.get_type() == GeomAbsSurfaceType::OtherSurface {
            // protect against non-geometric type (e.g. triangulation)
            return proj;
        }
        // double Tol = std::min(BRep_Tool::Tolerance(TheFace), Precision::Confusion());
        let tol = brep.tolerance(the_face).min(CONFUSION);
        // aTolU = std::max(mySurf.UResolution(Tol), Precision::PConfusion());
        let a_tol_u = my_surf.u_resolution(tol).max(PCONFUSION);
        let a_tol_v = my_surf.v_resolution(tol).max(PCONFUSION);
        // BRepTools::UVBounds(TheFace, U1, U2, V1, V2) — the face UV window
        // (TFaceData.uv_domain); the RealLast/RealFirst box is the OCCT
        // initial state of BRepTools::UVBounds for a face without edges.
        let [u1, u2, v1, v2] = match &*the_face.data {
            TShape::Face(fd) => fd.uv_domain,
            _ => None,
        }
        .unwrap_or([f64::MAX, f64::MIN, f64::MAX, f64::MIN]);
        // myExtPS.Initialize(mySurf, U1, U2, V1, V2, aTolU, aTolV) — the
        // rcad ExtPS fuses the Initialize with the Perform below.

        // ---- Perform(TheVertex, TheFace) (cxx L70-107) ----
        // const gp_Pnt P = BRep_Tool::Pnt(TheVertex); (the_pnt)
        // myExtPS.Perform(P);
        let world = my_surf.adaptor_surface().surface3().clone();
        let my_ext_ps = ExtPS::with_domain(the_pnt, &world, u1, u2, v1, v2, a_tol_u, a_tol_v);
        // Exploration of points and classification
        if my_ext_ps.is_done() {
            // BRepClass_FaceClassifier classifier; — the FaceShapeSource
            // view of the face (index 0 = the face) with the kernel ->
            // DS location table convention (slot 0 = identity).
            let locations_ds: Vec<glam::DAffine3> =
                std::iter::once(glam::DAffine3::IDENTITY)
                    .chain(brep.locations.iter().copied())
                    .collect();
            let fss = FaceShapeSource::new(the_face, world, &locations_ds);
            // const double Tol = BRep_Tool::Tolerance(TheFace);
            let tol_cls = brep.tolerance(the_face);
            for i in 1..=my_ext_ps.nb_ext() {
                // myExtPS.Point(i).Parameter(U1, U2);
                let puv = glam::DVec2::new(my_ext_ps.point(i).u, my_ext_ps.point(i).v);
                // classifier.Perform(TheFace, Puv, Tol);
                let mut classifier = FClassifier::new();
                classifier.perform(&fss, 0, puv, tol_cls);
                // const TopAbs_State state = classifier.State();
                let state = classifier.state();
                if state == State::On || state == State::In {
                    proj.my_sq_dist.push(my_ext_ps.square_distance(i));
                    proj.my_points.push(my_ext_ps.point(i).clone());
                }
            }
        }
        proj.my_ext_ps = Some(my_ext_ps);
        proj
    }

    /// OCCT NbExt() (hxx L59) — myPoints.Length() (the classified points).
    fn nb_ext(&self) -> usize {
        self.my_points.len()
    }

    /// OCCT SquareDistance(N) (hxx L61) — the value of the <N>th extremum
    /// square distance.
    fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist[n - 1]
    }

    /// OCCT Parameter(N, U, V) (hxx L63) — the parameters on the face of
    /// the <N>th extremum distance.
    fn parameter(&self, n: usize) -> (f64, f64) {
        (self.my_points[n - 1].u, self.my_points[n - 1].v)
    }
}

/// OCCT CurvatureValue(F, U, V, Tg) (cxx L29-57) — the signed curvature
/// value in the direction <Tg> at the <U,V> point on the surface <F>.
pub fn curvature_value(f: *mut Surface<'_>, u: f64, v: f64, tg: Vec3) -> f64 {
    // ((HLRBRep_Surface*)F)->D2(U, V, P, D1U, D1V, D2U, D2V, D2UV);
    let f_ref: &Surface = unsafe { &*f };
    let (_p, d1u, d1v, d2u, d2v, d2uv) = f_ref.d2(u, v);
    let _ = _p;
    let d1ut = d1u.dot(tg);
    let d1vt = d1v.dot(tg);
    let d1ud1v = d1u.dot(d1v);
    let nmu2 = d1u.dot(d1u);
    let nmv2 = d1v.dot(d1v);
    let det = nmu2 * nmv2 - d1ud1v * d1ud1v;
    if det > f64::MIN_POSITIVE {
        // gp::Resolution() (gp.hxx L60: RealSmall() = DBL_MIN).
        let alfa = (d1ut * nmv2 - d1vt * d1ud1v) / det;
        let beta = (d1vt * nmu2 - d1ut * d1ud1v) / det;
        let alfa2 = alfa * alfa;
        let beta2 = beta * beta;
        let alfabeta = alfa * beta;
        let mut nm = d1u.cross(d1v);
        nm = nm.normalize();
        let n = (nm.dot(d2u)) * alfa2 + 2.0 * (nm.dot(d2uv)) * alfabeta + (nm.dot(d2v)) * beta2;
        let d = nmu2 * alfa2 + 2.0 * d1ud1v * alfabeta + nmv2 * beta2;
        return n / d;
    }
    0.
}

/// OCCT UVPoint(Par, E, F, U, V) (cxx L61-105) — return true if U and V are
/// found: the UV coordinates at the parameter <Par> of the curve <E> on the
/// surface <F>.  `edge` is the OCCT `E->Curve().Edge()` and `brep` the
/// BRep_Tool context (see the module header note).
pub fn uv_point(
    brep: &rcad_kernel::BRep,
    par: f64,
    e: *mut Curve<'_>,
    edge: &Shape,
    f: *mut Surface<'_>,
    u: &mut f64,
    v: &mut f64,
) -> bool {
    // double pfbid, plbid; — the unused range outs of BRep_Tool::
    // CurveOnSurface (left unassigned in OCCT; the 0 defaults are the
    // neutral stand-in).
    let mut pfbid = 0.0;
    let mut plbid = 0.0;
    let e_ref: &Curve = unsafe { &*e };
    let f_ref: &Surface = unsafe { &*f };
    // if (BRep_Tool::CurveOnSurface(E->Curve().Edge(), F->Surface().Face(),
    //                               pfbid, plbid).IsNull())
    let c_on_s = brep.curve_on_surface(edge, f_ref.my_surface().face());
    if c_on_s.is_none() {
        // BRepExtrema_ExtPF proj(BRepLib_MakeVertex(E->Value3D(Par)),
        //                        F->Surface().Face());
        let proj = ExtPF::new(brep, e_ref.value_3d(par), f_ref.my_surface().face());
        // int i, index = 0;
        let mut index = 0usize;
        // double dist2 = RealLast();
        let mut dist2 = f64::MAX;
        // const int n = proj.NbExt();
        let n = proj.nb_ext();
        // for (i = 1; i <= n; i++)
        for i in 1..=n {
            // const double newdist2 = proj.SquareDistance(i);
            let newdist2 = proj.square_distance(i);
            if newdist2 < dist2 {
                dist2 = newdist2;
                index = i;
            }
        }
        let _ = dist2;
        if index == 0 {
            return false;
        }
        // proj.Parameter(index, U, V);
        let (pu, pv) = proj.parameter(index);
        *u = pu;
        *v = pv;
    } else {
        // BRepAdaptor_Curve2d PC(E->Curve().Edge(), F->Surface().Face());
        // — the adaptor loads the same BRep_Tool::CurveOnSurface
        // representation, and PC.D0(Par, P2d) is the pcurve D0.
        let (pc, _t1, _t2) = c_on_s.unwrap();
        pfbid = _t1;
        plbid = _t2;
        // gp_Pnt2d P2d; PC.D0(Par, P2d);
        let p2d = Curve2dAdaptor::value(&pc, par);
        // U = P2d.X(); V = P2d.Y();
        *u = p2d.x;
        *v = p2d.y;
    }
    let _ = (pfbid, plbid);
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::base::proj_lib::CurveType;
    use rcad_kernel::geom::{
        Circle3, Curve3, CylindricalSurface, Ellipse3, Line2d, Line3, Plane, Surface3,
    };
    use rcad_kernel::topo::topods::BRepBuilder;

    use crate::hlr::brep::b_curve_tool::CurveView;

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double).
    struct LineEdge {
        origin: Point3,
        len: f64,
    }

    impl CurveView for LineEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.len
        }
        fn d0(&self, u: f64) -> Point3 {
            self.origin + Vec3::new(u, 0.0, 0.0)
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (self.d0(u), Vec3::new(1.0, 0.0, 0.0))
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> Line3 {
            Line3::new(self.origin, Vec3::new(1.0, 0.0, 0.0))
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

    /// A 2x2 square planar face at z = `z` with per-edge pcurves (the
    /// BRep_Tool::CurveOnSurface representations).
    fn square_face(brep: &mut rcad_kernel::BRep, z: f64, with_pcurves: bool) -> Shape {
        let mut b = BRepBuilder::new();
        let p0 = DVec3::new(0.0, 0.0, z);
        let p1 = DVec3::new(2.0, 0.0, z);
        let p2 = DVec3::new(2.0, 2.0, z);
        let p3 = DVec3::new(0.0, 2.0, z);
        let vs = [
            b.add_vertex(brep, p0, 1e-7),
            b.add_vertex(brep, p1, 1e-7),
            b.add_vertex(brep, p2, 1e-7),
            b.add_vertex(brep, p3, 1e-7),
        ];
        let seg = |a: Point3, bb: Point3| Curve3::Line(Line3::new(a, (bb - a).normalize()));
        let e0 = b.add_edge(brep, Some(seg(p0, p1)), vs[0].clone(), vs[1].clone(), [0.0, 2.0]);
        let e1 = b.add_edge(brep, Some(seg(p1, p2)), vs[1].clone(), vs[2].clone(), [0.0, 2.0]);
        let e2 = b.add_edge(brep, Some(seg(p2, p3)), vs[2].clone(), vs[3].clone(), [0.0, 2.0]);
        let e3 = b.add_edge(brep, Some(seg(p3, p0)), vs[3].clone(), vs[0].clone(), [0.0, 2.0]);
        let wire = brep.add_twire(vec![e0.clone(), e1.clone(), e2.clone(), e3.clone()]);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: p0,
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 2.0]),
            Vec::new(),
            true,
        );
        if with_pcurves {
            // OCCT BRep_Builder::UpdateEdge(E, C, F, Tol) — the pcurve of
            // each wire edge on the face.
            let pcs = [
                Line2d {
                    origin: glam::DVec2::new(0.0, 0.0),
                    direction: glam::DVec2::new(1.0, 0.0),
                },
                Line2d {
                    origin: glam::DVec2::new(2.0, 0.0),
                    direction: glam::DVec2::new(0.0, 1.0),
                },
                Line2d {
                    origin: glam::DVec2::new(2.0, 2.0),
                    direction: glam::DVec2::new(-1.0, 0.0),
                },
                Line2d {
                    origin: glam::DVec2::new(0.0, 2.0),
                    direction: glam::DVec2::new(0.0, -1.0),
                },
            ];
            for (e, pc) in [&e0, &e1, &e2, &e3].into_iter().zip(pcs) {
                b.add_pcurve(brep, (*e).clone(), face.clone(), rcad_kernel::geom::Curve2d::Line(pc), 0.0, 2.0);
            }
        }
        face
    }

    /// OCCT anchor: the pcurve branch (cxx L97-103) — an edge carrying a
    /// BRep_Tool::CurveOnSurface representation evaluates PC.D0(Par, P2d)
    /// exactly: Par 1.5 of the x-edge of the z=5 plane is (1.5, 0).
    #[test]
    fn uv_point_plane_pcurve() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let face = square_face(&mut brep, 5.0, true);
        // An edge of the face: (0,0,5) -> (2,0,5) with pcurve (0,0)+X.
        let v1 = b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 5.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 5.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 5.0), DVec3::new(1.0, 0.0, 0.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        b.add_pcurve(
            &mut brep,
            e.clone(),
            face.clone(),
            rcad_kernel::geom::Curve2d::Line(Line2d {
                origin: glam::DVec2::new(0.0, 0.0),
                direction: glam::DVec2::new(1.0, 0.0),
            }),
            0.0,
            2.0,
        );

        let view = LineEdge {
            origin: Point3::new(0.0, 0.0, 5.0),
            len: 2.0,
        };
        let mut c = Curve::new();
        c.load(&view);
        let mut s = Surface::new();
        s.load(&brep, &face);

        let mut u = 0.0;
        let mut v = 0.0;
        let found = uv_point(&brep, 1.5, &mut c, &e, &mut s, &mut u, &mut v);
        assert!(found);
        assert!((u - 1.5).abs() < 1e-12, "u={u}");
        assert!(v.abs() < 1e-12, "v={v}");
    }

    /// OCCT anchor: the pcurve branch on a cylindrical face — the
    /// generatrix edge at angle 0 carries the pcurve u=0, v=Par; Par 1.0
    /// reads (0, 1).
    #[test]
    fn uv_point_cylinder_pcurve() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::new(1.5, 0.0, 0.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(1.5, 0.0, 2.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::new(1.5, 0.0, 0.0), DVec3::new(0.0, 0.0, 1.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e.clone()]);
        let face = brep.add_tface(
            Some(Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                radius: 1.5,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
                y_dir: None,
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, std::f64::consts::TAU, 0.0, 2.0]),
            Vec::new(),
            true,
        );
        b.add_pcurve(
            &mut brep,
            e.clone(),
            face.clone(),
            rcad_kernel::geom::Curve2d::Line(Line2d {
                origin: glam::DVec2::new(0.0, 0.0),
                direction: glam::DVec2::new(0.0, 1.0),
            }),
            0.0,
            2.0,
        );

        let view = LineEdge {
            origin: Point3::new(1.5, 0.0, 0.0),
            len: 0.0,
        };
        let mut c = Curve::new();
        c.load(&view);
        let mut s = Surface::new();
        s.load(&brep, &face);

        let mut u = -1.0;
        let mut v = -1.0;
        let found = uv_point(&brep, 1.0, &mut c, &e, &mut s, &mut u, &mut v);
        assert!(found);
        assert!(u.abs() < 1e-12, "u={u}");
        assert!((v - 1.0).abs() < 1e-12, "v={v}");
    }

    /// OCCT anchor: the extrema fallback (cxx L74-94) — an edge without a
    /// pcurve on the face projects E->Value3D(Par) through
    /// BRepExtrema_ExtPF and takes the closest classified candidate:
    /// (1.5, 1, 6) lands on the z=5 plane at (1.5, 1), inside the face.
    #[test]
    fn uv_point_extrema_fallback() {
        let mut brep = rcad_kernel::BRep::new();
        let face = square_face(&mut brep, 5.0, true);
        // The query edge floats above the plane and carries no pcurve on it;
        // Value3D(0.5) = (1.5, 0, 6) projects onto the boundary edge
        // (the classified-ON candidate is the one BRepExtrema_ExtPF keeps).
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::new(1.0, 0.0, 6.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 6.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::new(1.0, 0.0, 6.0), DVec3::new(1.0, 0.0, 0.0)))),
            v1,
            v2,
            [0.0, 1.0],
        );

        let view = LineEdge {
            origin: Point3::new(1.0, 0.0, 6.0),
            len: 1.0,
        };
        let mut c = Curve::new();
        c.load(&view);
        let mut s = Surface::new();
        s.load(&brep, &face);

        let mut u = 0.0;
        let mut v = 0.0;
        let found = uv_point(&brep, 0.5, &mut c, &e, &mut s, &mut u, &mut v);
        assert!(found);
        assert!((u - 1.5).abs() < 1e-9, "u={u}");
        assert!(v.abs() < 1e-9, "v={v}");
    }

    /// OCCT anchor: the planar CurveOnPlane fallback routes the fixture to
    /// the pcurve branch of UVPoint (HLRBRep_EdgeFaceTool.cxx L68-81).  The
    /// edge carries no STORED pcurve, but BRep_Tool::CurveOnSurface computes
    /// one on the fly (BRep_Tool.cxx L367-372 -> CurveOnPlane L379-450:
    /// planar face + 3D curve -> ProjLib projection), so the IsNull branch
    /// is NOT taken: BRepAdaptor_Curve2d D0(1.0) = (1, 0) and UVPoint
    /// returns true.
    #[test]
    fn uv_point_plane_curve_on_plane_fallback_pcurve() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 5.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 5.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 5.0), DVec3::new(1.0, 0.0, 0.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e.clone()]);
        // No uv_domain: the face carries no UV window.
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::new(0.0, 0.0, 5.0),
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );

        let view = LineEdge {
            origin: Point3::new(0.0, 0.0, 5.0),
            len: 2.0,
        };
        let mut c = Curve::new();
        c.load(&view);
        let mut s = Surface::new();
        s.load(&brep, &face);

        let mut u = 0.0;
        let mut v = 0.0;
        let found = uv_point(&brep, 1.0, &mut c, &e, &mut s, &mut u, &mut v);
        assert!(found);
        assert!((u - 1.0).abs() < 1e-12, "u={u}");
        assert!(v.abs() < 1e-12, "v={v}");
    }

    /// OCCT anchor: no candidate classifies IN/ON -> index stays 0 and
    /// UVPoint returns false (cxx L88-91).  Reachable only for a NON-planar
    /// face: CurveOnPlane returns null there (BRep_Tool.cxx L400-404), so
    /// the IsNull extrema branch runs; with no UV window the extrema list
    /// stays empty.
    #[test]
    fn uv_point_not_found_returns_false() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 5.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 7.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::new(2.0, 0.0, 5.0), DVec3::new(0.0, 0.0, 1.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e.clone()]);
        // No uv_domain: the face carries no UV window.
        let face = brep.add_tface(
            Some(Surface3::Cylinder(CylindricalSurface {
                origin: DVec3::new(0.0, 0.0, 0.0),
                axis: DVec3::new(0.0, 0.0, 1.0),
                radius: 2.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
                y_dir: None,
            })),
            wire,
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );

        let view = LineEdge {
            origin: Point3::new(2.0, 0.0, 5.0),
            len: 2.0,
        };
        let mut c = Curve::new();
        c.load(&view);
        let mut s = Surface::new();
        s.load(&brep, &face);

        let mut u = 0.0;
        let mut v = 0.0;
        let found = uv_point(&brep, 1.0, &mut c, &e, &mut s, &mut u, &mut v);
        assert!(!found);
    }

    /// OCCT anchor: CurvatureValue on a plane (cxx L29-57) — the second
    /// derivatives vanish, so N = 0 and the returned value is exactly 0.
    #[test]
    fn curvature_value_plane_is_zero() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 0.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e.clone()]);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 2.0]),
            Vec::new(),
            true,
        );
        let mut s = Surface::new();
        s.load(&brep, &face);

        // Any tangent direction: the normal curvature of a plane is 0.
        let curv = curvature_value(&mut s, 0.5, 0.5, Vec3::new(1.0, 0.0, 0.0));
        assert!(curv.abs() < 1e-12, "curv={curv}");
    }

    /// OCCT anchor: CurvatureValue on a unit cylinder (the Stage 3b
    /// sl_props fixture) — the circular direction carries the principal
    /// curvature -1/R = -1 and the axial direction 0.
    #[test]
    fn curvature_value_cylinder_principal_directions() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 2.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(Curve3::Line(Line3::new(DVec3::ZERO, DVec3::new(0.0, 0.0, 1.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e.clone()]);
        let face = brep.add_tface(
            Some(Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
                y_dir: None,
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, std::f64::consts::TAU, 0.0, 2.0]),
            Vec::new(),
            true,
        );
        let mut s = Surface::new();
        s.load(&brep, &face);

        // At (u,v) = (0,1): D1U = (0,1,0) (circular), D1V = (0,0,1)
        // (axial); the Tg = D1U direction reads -1/R = -1.
        let curv_circ = curvature_value(&mut s, 0.0, 1.0, Vec3::new(0.0, 1.0, 0.0));
        assert!(
            (curv_circ + 1.0).abs() < 1e-5,
            "circular curvature = {curv_circ}"
        );
        // The Tg = D1V direction reads 0 (the axial curvature).
        let curv_axial = curvature_value(&mut s, 0.0, 1.0, Vec3::new(0.0, 0.0, 1.0));
        assert!(curv_axial.abs() < 1e-5, "axial curvature = {curv_axial}");
    }
}
