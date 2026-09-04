// OCCT HLRBRep_Surface (TKHLR) — the face surface to be projected.
//
// HLRBRep_Surface.hxx L41-193 + .cxx L28-335 + .lxx L22-271 (delegates to
// the BRepAdaptor_Surface through HLRBRep_BSurfaceTool).
//
// `myProj` mirrors the OCCT raw pointer `const HLRAlgo_Projector*` — set
// externally by the owning Data (the projector outlives the surface, same
// object graph as OCCT HLRBRep_Data).
//
// The `SurfaceAdapter` impl lets the Contap engine consume an
// HLRBRep_Surface directly (the Stage 3b wiring goal).

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::base::gprop::pequation::PEquation;
use rcad_kernel::geom::{Plane, Point3, Surface3, Vec3};

use crate::hlr::algo::projector::Projector;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

use super::b_curve_tool::CurveView;

/// OCCT HLRBRep_Surface.
#[derive(Clone)]
pub struct Surface<'a> {
    /// OCCT mySurf (BRepAdaptor_Surface, restricted).
    my_surf: BRepAdaptorSurface<'a>,
    /// OCCT myType — the classification after Surface(F): canonic types are
    /// kept, a degree-1 Bezier becomes Plane, everything else OtherSurface.
    my_type: GeomAbsSurfaceType,
    /// OCCT `const HLRAlgo_Projector* myProj`.
    my_proj: *const Projector,
}

// The raw projector pointer is owned by the caller Data (never mutated
// through the surface); the OCCT code only const-casts it for cache-free
// reads.
unsafe impl Send for Surface<'_> {}

impl<'a> Surface<'a> {
    /// OCCT HLRBRep_Surface() (cxx L28-32) — an undefined surface with no
    /// face loaded.
    pub fn new() -> Self {
        Surface {
            my_surf: BRepAdaptorSurface::new(),
            my_type: GeomAbsSurfaceType::OtherSurface,
            my_proj: std::ptr::null(),
        }
    }

    /// OCCT Projector(const HLRAlgo_Projector* Proj) (hxx L49).
    pub fn projector(&mut self, proj: *const Projector) {
        self.my_proj = proj;
    }

    /// OCCT Surface() (lxx L22-25) — the face adaptor.
    pub fn my_surface(&self) -> &BRepAdaptorSurface<'a> {
        &self.my_surf
    }

    /// OCCT Surface(const TopoDS_Face& F) (cxx L36-68) — loads the face and
    /// classifies: canonic types stay, a degree-1 Bezier becomes Plane,
    /// everything else is OtherSurface.
    pub fn load(&mut self, brep: &'a rcad_kernel::BRep, f: &rcad_kernel::topo::topods::Shape) {
        // mySurf.Initialize(F, true);
        self.my_surf = BRepAdaptorSurface::initialize_face(brep, f, true);
        let typ = self.my_surf.get_type();
        match typ {
            GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus => {
                // unchanged type
                self.my_type = typ;
            }

            GeomAbsSurfaceType::BezierSurface => {
                if self.my_surf.u_degree() == 1 && self.my_surf.v_degree() == 1 {
                    self.my_type = GeomAbsSurfaceType::Plane;
                } else {
                    self.my_type = typ;
                }
            }

            _ => {
                self.my_type = GeomAbsSurfaceType::OtherSurface;
            }
        }
    }

    fn proj(&self) -> &Projector {
        assert!(
            !self.my_proj.is_null(),
            "HLRBRep_Surface: no projector set"
        );
        unsafe { &*self.my_proj }
    }

    /// OCCT SideRowsOfPoles (cxx L72-145) — the projected poles build side
    /// iso-rows or sit in a Z-parallel plane. `pnt` is the pole grid
    /// (rows = U index; the OCCT Array2).
    fn side_rows_of_poles(&self, tol: f64, pnt: &[Vec<Point3>]) -> bool {
        let nbu_poles = pnt.len();
        let nbv_poles = pnt.first().map_or(0, |r| r.len());
        let t = self.proj().transformation().clone();
        let mut transformed = vec![vec![Point3::ZERO; nbv_poles]; nbu_poles];
        for (iu, row) in pnt.iter().enumerate() {
            for (iv, p) in row.iter().enumerate() {
                transformed[iu][iv] = t.apply(*p);
            }
        }
        let mut result = true;

        for iu in 0..nbu_poles {
            if !result {
                break;
            }
            // Side iso u ?
            let x0 = transformed[iu][0].x;
            let y0 = transformed[iu][0].y;
            for iv in 1..nbv_poles {
                if !result {
                    break;
                }
                let (x, y, _) = (transformed[iu][iv].x, transformed[iu][iv].y, 0.0);
                result = (x - x0).abs() < tol && (y - y0).abs() < tol;
            }
        }
        if result {
            return result;
        }
        result = true;

        for iv in 0..nbv_poles {
            if !result {
                break;
            }
            // Side iso v ?
            let x0 = transformed[0][iv].x;
            let y0 = transformed[0][iv].y;
            for iu in 1..nbu_poles {
                if !result {
                    break;
                }
                let (x, y, _) = (transformed[iu][iv].x, transformed[iu][iv].y, 0.0);
                result = (x - x0).abs() < tol && (y - y0).abs() < tol;
            }
        }
        if result {
            return result;
        }

        // Are the Poles in a Side Plane ? (cxx L125-144)
        let mut poles: Vec<Point3> = Vec::with_capacity(nbu_poles * nbv_poles);
        for row in &transformed {
            poles.extend_from_slice(row);
        }
        let pl = PEquation::new(&poles, tol);
        if pl.is_planar() {
            result = pl.plane().normal.z.abs() < 0.0001;
        }

        result
    }

    /// OCCT IsSide (cxx L149-240).
    pub fn is_side(&self, tol_f: f64, toler: f64) -> bool {
        let proj = self.proj();
        if self.my_type == GeomAbsSurfaceType::Plane {
            let pl = self.plane();
            // gp_Ax1 A = Pl.Axis().
            let pt = pl.origin;
            let mut d = pl.normal;
            let trsf = proj.transformation().clone();
            let pt = trsf.apply(pt);
            d = trsf.transform_vec(d);
            let r = if proj.perspective() {
                d.z * proj.focus() - (d.x * pt.x + d.y * pt.y + d.z * pt.z)
            } else {
                d.z
            };
            r.abs() < toler
        } else if self.my_type == GeomAbsSurfaceType::Cylinder {
            if proj.perspective() {
                return false;
            }
            let cyl = self.my_surf.cylinder();
            let mut d = cyl.axis;
            let trsf = proj.transformation().clone();
            d = trsf.transform_vec(d);
            let r = (d.x * d.x + d.y * d.y).sqrt();
            r < toler
        } else if self.my_type == GeomAbsSurfaceType::Cone {
            if !proj.perspective() {
                return false;
            }
            let con = self.my_surf.cone();
            let trsf = proj.transformation().clone();
            let pt = trsf.apply(con.apex);
            let tol = 0.001;
            // Pt.IsEqual(gp_Pnt(0, 0, Focus()), tol)
            pt.distance(Point3::new(0.0, 0.0, proj.focus())) <= tol
        } else if self.my_type == GeomAbsSurfaceType::BezierSurface
            || self.my_type == GeomAbsSurfaceType::BSplineSurface
        {
            if proj.perspective() {
                return false;
            }
            let pnt = self.my_surf.poles_grid();
            if pnt.is_empty() {
                return false;
            }
            // int nu = NbUPoles(mySurf); int nv = NbVPoles(mySurf);
            // ... NCollection_Array2 Pnt(1, nu, 1, nv) ... (cxx L204-215 /
            // L223-234)
            let _ = (pnt.len(), pnt[0].len());
            self.side_rows_of_poles(tol_f, &pnt)
        } else {
            false
        }
    }

    /// OCCT IsAbove (cxx L244-306) — the planar-face vs curve side test.
    pub fn is_above(&self, back: bool, a: &dyn CurveView, tol: f64) -> bool {
        let planar = self.my_type == GeomAbsSurfaceType::Plane;
        if planar {
            let pl = self.plane();
            // gp_Pln::Coefficients(a, b, c, d): n . X + d = 0.
            let (a_c, b_c, c_c, d_c) = (
                pl.normal.x,
                pl.normal.y,
                pl.normal.z,
                -pl.normal.dot(pl.origin),
            );
            let mut u;
            let u1 = a.parameter_3d(a.first_parameter());
            let u2 = a.parameter_3d(a.last_parameter());
            u = u1;
            let p = a.d0(u);
            let mut dd = a_c * p.x + b_c * p.y + c_c * p.z + d_c;
            if back {
                dd = -dd;
            }
            if dd < -tol {
                return false;
            }
            if a.get_type() != CurveType::Line {
                let nb_pnt = 30;
                let step = (u2 - u1) / (nb_pnt + 1) as f64;
                for _ in 1..=nb_pnt {
                    u += step;
                    let p = a.d0(u);
                    let mut dd = a_c * p.x + b_c * p.y + c_c * p.z + d_c;
                    if back {
                        dd = -dd;
                    }
                    if dd < -tol {
                        return false;
                    }
                }
            }
            u = u2;
            let p = a.d0(u);
            let mut dd = a_c * p.x + b_c * p.y + c_c * p.z + d_c;
            if back {
                dd = -dd;
            }
            if dd < -tol {
                return false;
            }
            true
        } else {
            false
        }
    }

    /// OCCT Value(U, V) (cxx L310-315).
    pub fn value(&self, u: f64, v: f64) -> Point3 {
        self.my_surf.value(u, v)
    }

    /// OCCT Plane (cxx L319-335) — the Bezier branch fits the plane at
    /// (0.5, 0.5): gp_Pln(P, gp_Dir(D1U.Crossed(D1V))).
    pub fn plane(&self) -> Plane {
        let typ = self.my_surf.get_type();
        match typ {
            GeomAbsSurfaceType::BezierSurface => {
                let (p, d1u, d1v) = self.my_surf.d1(0.5, 0.5);
                Plane::new(p, d1u.cross(d1v))
            }
            _ => self.my_surf.plane(),
        }
    }

    // ---- lxx one-liners (through BSurfaceTool) ----

    /// OCCT FirstUParameter (lxx L36-41).
    pub fn first_u_parameter(&self) -> f64 {
        self.my_surf.first_u_parameter()
    }
    /// OCCT LastUParameter (lxx L43-48).
    pub fn last_u_parameter(&self) -> f64 {
        self.my_surf.last_u_parameter()
    }
    /// OCCT FirstVParameter (lxx L50-55).
    pub fn first_v_parameter(&self) -> f64 {
        self.my_surf.first_v_parameter()
    }
    /// OCCT LastVParameter (lxx L57-62).
    pub fn last_v_parameter(&self) -> f64 {
        self.my_surf.last_v_parameter()
    }
    /// OCCT IsUClosed (lxx L92-97).
    pub fn is_u_closed(&self) -> bool {
        self.my_surf.is_u_closed()
    }
    /// OCCT IsVClosed (lxx L99-104).
    pub fn is_v_closed(&self) -> bool {
        self.my_surf.is_v_closed()
    }
    /// OCCT IsUPeriodic (lxx L106-111).
    pub fn is_u_periodic(&self) -> bool {
        self.my_surf.is_u_periodic()
    }
    /// OCCT UPeriod (lxx L113-118).
    pub fn u_period(&self) -> f64 {
        self.my_surf.u_period()
    }
    /// OCCT IsVPeriodic (lxx L120-125).
    pub fn is_v_periodic(&self) -> bool {
        self.my_surf.is_v_periodic()
    }
    /// OCCT VPeriod (lxx L127-132).
    pub fn v_period(&self) -> f64 {
        self.my_surf.v_period()
    }
    /// OCCT D0 (lxx L134-139).
    pub fn d0(&self, u: f64, v: f64) -> Point3 {
        self.my_surf.value(u, v)
    }
    /// OCCT D1 (lxx L141-150).
    pub fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        self.my_surf.d1(u, v)
    }
    /// OCCT D2 (lxx L152-163).
    pub fn d2(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        self.my_surf.d2(u, v)
    }
    /// OCCT GetType (lxx L191-194).
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        self.my_type
    }
    /// OCCT Cylinder (lxx L198-203).
    pub fn cylinder(&self) -> rcad_kernel::geom::CylindricalSurface {
        self.my_surf.cylinder()
    }
    /// OCCT Cone (lxx L205-210).
    pub fn cone(&self) -> rcad_kernel::geom::ConicalSurface {
        self.my_surf.cone()
    }
    /// OCCT Sphere (lxx L212-217).
    pub fn sphere(&self) -> rcad_kernel::geom::SphericalSurface {
        self.my_surf.sphere()
    }
    /// OCCT Torus (lxx L219-224).
    pub fn torus(&self) -> rcad_kernel::geom::ToroidalSurface {
        self.my_surf.torus()
    }
    /// OCCT Axis (lxx L268-273) — AxeOfRevolution.
    pub fn axis(&self) -> (Point3, Vec3) {
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::Cylinder => {
                let c = self.my_surf.cylinder();
                (c.origin, c.axis)
            }
            GeomAbsSurfaceType::Cone => {
                let c = self.my_surf.cone();
                (c.apex, c.axis)
            }
            GeomAbsSurfaceType::Sphere => {
                let c = self.my_surf.sphere();
                (c.center, c.axis)
            }
            GeomAbsSurfaceType::Torus => {
                let c = self.my_surf.torus();
                (c.center, c.axis)
            }
            _ => panic!("Standard_NoSuchObject: HLRBRep_Surface::Axis"),
        }
    }
}

impl Default for Surface<'_> {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT Adaptor3d_Surface interface — HLRBRep_Surface feeds the Contap
/// engine directly (Stage 3b wiring).
impl crate::hlr::contap::surface_adaptor::SurfaceAdapter for Surface<'_> {
    fn first_u_parameter(&self) -> f64 {
        Surface::first_u_parameter(self)
    }
    fn last_u_parameter(&self) -> f64 {
        Surface::last_u_parameter(self)
    }
    fn first_v_parameter(&self) -> f64 {
        Surface::first_v_parameter(self)
    }
    fn last_v_parameter(&self) -> f64 {
        Surface::last_v_parameter(self)
    }
    fn value(&self, u: f64, v: f64) -> Point3 {
        Surface::value(self, u, v)
    }
    fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        Surface::d1(self, u, v)
    }
    fn d2(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        Surface::d2(self, u, v)
    }
    fn u_resolution(&self, r3d: f64) -> f64 {
        self.my_surf.u_resolution(r3d)
    }
    fn v_resolution(&self, r3d: f64) -> f64 {
        self.my_surf.v_resolution(r3d)
    }
    fn get_type(&self) -> GeomAbsSurfaceType {
        Surface::get_type(self)
    }
    fn plane(&self) -> Plane {
        Surface::plane(self)
    }
    fn cylinder(&self) -> rcad_kernel::geom::CylindricalSurface {
        Surface::cylinder(self)
    }
    fn cone(&self) -> rcad_kernel::geom::ConicalSurface {
        Surface::cone(self)
    }
    fn sphere(&self) -> rcad_kernel::geom::SphericalSurface {
        Surface::sphere(self)
    }
    fn torus(&self) -> rcad_kernel::geom::ToroidalSurface {
        Surface::torus(self)
    }
    fn nb_u_poles(&self) -> usize {
        self.my_surf.nb_u_poles()
    }
    fn nb_v_poles(&self) -> usize {
        self.my_surf.nb_v_poles()
    }
    fn nb_u_knots(&self) -> usize {
        self.my_surf.nb_u_knots()
    }
    fn nb_v_knots(&self) -> usize {
        self.my_surf.nb_v_knots()
    }
    fn u_degree(&self) -> usize {
        self.my_surf.u_degree()
    }
    fn v_degree(&self) -> usize {
        self.my_surf.v_degree()
    }
    fn is_u_periodic(&self) -> bool {
        Surface::is_u_periodic(self)
    }
    fn u_period(&self) -> f64 {
        Surface::u_period(self)
    }
    fn is_v_periodic(&self) -> bool {
        Surface::is_v_periodic(self)
    }
    fn v_period(&self) -> f64 {
        Surface::v_period(self)
    }
    fn kernel_surface(&self) -> Option<&Surface3> {
        // The adaptor owns the located world surface by value.
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{CylindricalSurface, Line3, Plane as GPlane, Surface3};
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::topo::topods::{BRepBuilder, Shape};

    use crate::hlr::brep::b_curve_tool::CurveView;

    /// A straight 3D segment as the CurveView test double.
    struct TestLine(f64);

    impl CurveView for TestLine {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.0
        }
        fn d0(&self, u: f64) -> Point3 {
            Point3::new(u, 0.0, 6.0)
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (Point3::new(u, 0.0, 6.0), Vec3::new(1.0, 0.0, 0.0))
        }
        fn get_type(&self) -> rcad_kernel::base::proj_lib::CurveType {
            rcad_kernel::base::proj_lib::CurveType::Line
        }
        fn line(&self) -> rcad_kernel::geom::Line3 {
            Line3 {
                origin: Point3::ZERO,
                direction: Vec3::new(1.0, 0.0, 0.0),
            }
        }
        fn circle(&self) -> rcad_kernel::geom::Circle3 {
            panic!("no circle");
        }
        fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
            panic!("no ellipse");
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

    fn plane_face(brep: &mut rcad_kernel::BRep, z: f64) -> rcad_kernel::topo::topods::Shape {
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(brep, DVec3::new(0.0, 0.0, z), 1e-7);
        let v2 = b.add_vertex(brep, DVec3::new(2.0, 0.0, z), 1e-7);
        let e = b.add_edge(
            brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, z),
                direction: DVec3::new(1.0, 0.0, 0.0),
            })),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep.add_tface(
            Some(Surface3::Plane(GPlane {
                origin: DVec3::new(0.0, 0.0, z),
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
        face
    }

    /// OCCT anchor: a plane viewed along its own normal is a side face
    /// (IsSide plane branch: the transformed normal has D.Z() = 1, cxx
    /// L155-172); a curve above the plane passes IsAbove, below fails
    /// (cxx L244-306).
    #[test]
    fn hlr_surface_is_side_and_is_above() {
        let mut brep = rcad_kernel::BRep::new();
        let face = plane_face(&mut brep, 5.0);

        let mut s = Surface::new();
        s.load(&brep, &face);
        assert_eq!(s.get_type(), GeomAbsSurfaceType::Plane);

        let proj = Projector::from_ax2(&Ax2::new(DVec3::ZERO, DVec3::new(0.0, 0.0, 1.0), DVec3::new(1.0, 0.0, 0.0)));
        s.projector(&proj);

        // IsSide: the transformed plane normal is +Z (r = D.Z() = 1 — not
        // within toler of 0 for the direct view; the SIDE planes face the
        // viewer edge-on, n.z = 0 — a flipped projector makes r = -1).
        assert!(!s.is_side(1e-3, 1e-3));

        // IsAbove: the curve at z = 6 is above the z = 5 plane.
        let curve = TestLine(2.0);
        assert!(s.is_above(false, &curve, 1e-6));
        // Back view: the same curve is then below (dd flips sign).
        assert!(!s.is_above(true, &curve, 1e-6));

        // A curve below the plane (z = 4): at the start point dd = -1.
        struct Below(f64);
        impl CurveView for Below {
            fn first_parameter(&self) -> f64 {
                0.0
            }
            fn last_parameter(&self) -> f64 {
                self.0
            }
            fn d0(&self, u: f64) -> Point3 {
                Point3::new(u, 0.0, 4.0)
            }
            fn d1(&self, _u: f64) -> (Point3, Vec3) {
                (Point3::ZERO, Vec3::X)
            }
            fn get_type(&self) -> rcad_kernel::base::proj_lib::CurveType {
                rcad_kernel::base::proj_lib::CurveType::Line
            }
            fn line(&self) -> rcad_kernel::geom::Line3 {
                Line3 {
                    origin: Point3::ZERO,
                    direction: Vec3::new(1.0, 0.0, 0.0),
                }
            }
            fn circle(&self) -> rcad_kernel::geom::Circle3 {
                panic!("no circle");
            }
            fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
                panic!("no ellipse");
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
        let below = Below(2.0);
        assert!(!s.is_above(false, &below, 1e-6));
    }

    /// OCCT anchor: a cylinder viewed along its own axis is a side face —
    /// the transformed axis is Z, r = 0 < toler (cxx L173-185).
    #[test]
    fn hlr_surface_cylinder_is_side() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 2.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::new(0.0, 0.0, 1.0),
            })),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep
            .add_tface(
                Some(Surface3::Cylinder(CylindricalSurface {
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

        let mut s = Surface::new();
        s.load(&brep, &face);
        assert_eq!(s.get_type(), GeomAbsSurfaceType::Cylinder);

        let proj = Projector::from_ax2(&Ax2::new(DVec3::ZERO, DVec3::new(0.0, 0.0, 1.0), DVec3::new(1.0, 0.0, 0.0)));
        s.projector(&proj);
        // The axis transforms to (0, 0, 1): r = sqrt(0) = 0 < toler.
        assert!(s.is_side(1e-3, 1e-3));
    }
}
