// OCCT HLRBRep_Curve (TKHLR) — a 3D curve seen as its projection on the
// viewing plane, with an optional perspective.
//
// HLRBRep_Curve.hxx L47-211 + .cxx L39-632 + .lxx (the 3D/2D delegates).
//
// `myProj` mirrors the OCCT raw pointer; `myCurve` is the edge curve
// adaptor through the [`CurveView`] (HLRBRep_BCurveTool) interface.
//
// Documented deferral: `tangent` (cxx L301-313) depends on
// HLRBRep_CLProps (the GeomLProp_CLProps instantiation over the projected
// curve) — Stage 3a-3.

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Circle2d, Circle3, Curve2d, Ellipse2d, Line2d, Point3, Vec3};
use rcad_kernel::precision::ANGULAR;

use crate::hlr::algo::projector::Projector;
use crate::hlr::brep::b_curve_tool::CurveView;

/// OCCT HLRBRep_Curve.
pub struct Curve<'a> {
    /// OCCT myCurve (BRepAdaptor_Curve).
    my_curve: Option<&'a dyn CurveView>,
    /// OCCT myType — the projected-curve classification from Update().
    my_type: CurveType,
    /// OCCT `const HLRAlgo_Projector* myProj`.
    my_proj: *const Projector,
    my_ox: f64,
    my_oz: f64,
    my_vx: f64,
    my_vz: f64,
    my_of: f64,
}

unsafe impl Send for Curve<'_> {}

impl<'a> Curve<'a> {
    /// OCCT HLRBRep_Curve() (cxx L39).
    pub fn new() -> Self {
        Curve {
            my_curve: None,
            my_type: CurveType::Other,
            my_proj: std::ptr::null(),
            my_ox: 0.0,
            my_oz: 0.0,
            my_vx: 0.0,
            my_vz: 0.0,
            my_of: 0.0,
        }
    }

    /// OCCT Projector(const HLRAlgo_Projector* Proj) (hxx L55).
    pub fn projector(&mut self, proj: *const Projector) {
        self.my_proj = proj;
    }

    /// OCCT Curve(const TopoDS_Edge& E) (cxx L43-46) — loads the edge curve
    /// adaptor.
    pub fn load(&mut self, c: &'a dyn CurveView) {
        self.my_curve = Some(c);
    }

    /// OCCT GetCurve (hxx L64) — the loaded curve adaptor.
    pub fn view(&self) -> &'a dyn CurveView {
        self.my_curve.expect("HLRBRep_Curve::Curve not loaded")
    }

    fn proj(&self) -> &Projector {
        assert!(!self.my_proj.is_null(), "HLRBRep_Curve: no projector set");
        unsafe { &*self.my_proj }
    }

    fn c(&self) -> &'a dyn CurveView {
        self.view()
    }

    /// OCCT Parameter2d (cxx L50-75).
    pub fn parameter_2d(&self, p3d: f64) -> f64 {
        match self.my_type {
            CurveType::Line => {
                if self.proj().perspective() {
                    let fm_oz = self.my_of - self.my_oz;
                    return self.my_of * p3d * (self.my_vx * fm_oz + self.my_ox * self.my_vz)
                        / (fm_oz * (fm_oz - p3d * self.my_vz));
                }
                p3d * self.my_vx
            }
            CurveType::Ellipse => p3d + self.my_ox,
            _ => p3d,
        }
    }

    /// OCCT Parameter3d (cxx L79-104).
    pub fn parameter_3d(&self, p2d: f64) -> f64 {
        if self.my_type == CurveType::Line {
            if self.proj().perspective() {
                let fm_oz = self.my_of - self.my_oz;
                return p2d * fm_oz * fm_oz
                    / (fm_oz * (self.my_of * self.my_vx + p2d * self.my_vz)
                        + self.my_of * self.my_ox * self.my_vz);
            }
            // gp::Resolution()
            return if self.my_vx <= 1e-15 { p2d } else { p2d / self.my_vx };
        } else if self.my_type == CurveType::Ellipse {
            return p2d - self.my_ox;
        }
        p2d
    }

    /// OCCT Update (cxx L108-224) — the projected-curve classification and
    /// the line-view coefficients.
    pub fn update(&mut self, tot_min: &mut [f64; 16], tot_max: &mut [f64; 16]) -> f64 {
        let typ = self.c().get_type();
        self.my_type = CurveType::Other;

        match typ {
            CurveType::Line => {
                self.my_type = typ;
            }

            CurveType::Circle => {
                if !self.proj().perspective() {
                    let circ = self.c().circle();
                    // gp_Dir D1 = Circle.Axis().Direction(), transformed.
                    let trsf = self.proj().transformation().clone();
                    let mut d1 = trsf.transform_vec(circ.normal);
                    d1 = d1.normalize();
                    // D1.IsParallel(gp::DZ(), Precision::Angular())
                    if (d1.z.abs() - 1.0).abs() < ANGULAR {
                        self.my_type = CurveType::Circle;
                    } else if d1.z.abs() < ANGULAR * 10.0 {
                        //*10: The minor radius of ellipse should not be too small.
                        self.my_type = CurveType::Other;
                    } else {
                        self.my_type = CurveType::Ellipse;
                        // compute the angle offset
                        // gp_Dir D3 = D1.Crossed(gp::DZ());
                        let d3 = d1.cross(Vec3::new(0.0, 0.0, 1.0)).normalize();
                        // gp_Dir D2 = Circle.XAxis().Direction(), transformed.
                        let mut d2 = trsf.transform_vec(circ.x_dir);
                        d2 = d2.normalize();
                        // myOX = D3.AngleWithRef(D2, D1)
                        self.my_ox = angle_with_ref(d2, d3, d1);
                    }
                }
            }

            CurveType::Ellipse => {
                if !self.proj().perspective() {
                    let ell = self.c().ellipse();
                    let trsf = self.proj().transformation().clone();
                    let mut d1 = trsf.transform_vec(ell.normal);
                    d1 = d1.normalize();
                    if (d1.z.abs() - 1.0).abs() < ANGULAR {
                        self.my_ox = 0.0; // no offset on the angle
                        self.my_type = CurveType::Ellipse;
                    }
                }
            }

            CurveType::Bezier => {
                if self.c().degree() == 1 {
                    self.my_type = CurveType::Line;
                } else if !self.proj().perspective() {
                    self.my_type = typ;
                }
            }

            CurveType::BSpline => {
                if !self.proj().perspective() {
                    self.my_type = typ;
                }
            }

            _ => {}
        }

        if self.my_type == CurveType::Line {
            // compute the values for a line
            let l3d; // length of the 3d bezier curve
            let (mut p, mut v);
            if self.c().get_type() == CurveType::Line {
                let l = self.c().line();
                p = l.origin;
                v = l.direction;
                l3d = 1.0;
            } else {
                // bezier degree 1
                let (pl, vl) = self.c().d1(0.0);
                p = pl;
                v = vl;
                l3d = p.distance(self.c().d0(1.0));
            }
            let trsf = self.proj().transformation().clone();
            p = trsf.apply(p);
            v = trsf.transform_vec(v);
            if self.proj().perspective() {
                let mut f = glam::DVec2::ZERO;
                let mut vfx = glam::DVec2::ZERO;
                self.d1_2d(0.0, &mut f, &mut vfx);
                vfx = vfx.normalize();
                self.my_vx = (vfx.x * v.x + vfx.y * v.y) * l3d;
                let l = -(vfx.x * f.x + vfx.y * f.y);
                f = glam::DVec2::new(f.x + vfx.x * l, f.y + vfx.y * l);
                self.my_ox = vfx.x * (p.x - f.x) + vfx.y * (p.y - f.y);
                let mut vfz = Vec3::new(-f.x, -f.y, self.proj().focus());
                self.my_of = vfz.length();
                vfz /= self.my_of;
                self.my_vz = vfz.dot(v);
                self.my_vz *= l3d;
                self.my_oz = vfz.dot(Vec3::new(p.x - f.x, p.y - f.y, p.z));
            } else {
                self.my_vx = (v.x * v.x + v.y * v.y).sqrt() * l3d;
            }
        }
        self.update_min_max(tot_min, tot_max)
    }

    /// OCCT UpdateMinMax (cxx L228-287).
    pub fn update_min_max(&self, tot_min: &mut [f64; 16], tot_max: &mut [f64; 16]) -> f64 {
        let mut a = self.c().first_parameter();
        let b = self.c().last_parameter();
        let mut tol_min_max = 0.0f64;
        let (mut x, mut y, mut z) = self.proj_xyz(self.value_3d(a));
        crate::hlr::algo::hlr_algo::HLRAlgo::update_min_max(x, y, z, tot_min, tot_max);

        if self.my_type != CurveType::Line {
            let nb_pnt = 30;
            let step = (b - a) / (nb_pnt + 1) as f64;
            let (mut xa, mut ya, mut za);
            let (mut xb, mut yb, mut zb) = (0.0, 0.0, 0.0);

            for i in 1..=nb_pnt {
                a += step;
                xa = xb;
                ya = yb;
                za = zb;
                xb = x;
                yb = y;
                zb = z;
                let (x2, y2, z2) = self.proj_xyz(self.value_3d(a));
                x = x2;
                y = y2;
                z = z2;
                crate::hlr::algo::hlr_algo::HLRAlgo::update_min_max(x, y, z, tot_min, tot_max);
                if i >= 2 {
                    let (dx1, dy1, dz1) = (x - xa, y - ya, z - za);
                    let dd1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt();
                    if dd1 > 0.0 {
                        let (dx2, dy2, dz2) = (xb - xa, yb - ya, zb - za);
                        let dd2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt();
                        if dd2 > 0.0 {
                            let p = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (dd1 * dd2);
                            let (ex1, ey1, ez1) = (xa + p * dx1 - xb, ya + p * dy1 - yb, za + p * dz1 - zb);
                            let dd1b = (ex1 * ex1 + ey1 * ey1 + ez1 * ez1).sqrt();
                            if dd1b > tol_min_max {
                                tol_min_max = dd1b;
                            }
                        }
                    }
                }
            }
        }
        let (x3, y3, z3) = self.proj_xyz(self.value_3d(b));
        crate::hlr::algo::hlr_algo::HLRAlgo::update_min_max(x3, y3, z3, tot_min, tot_max);
        let _ = (x, y, z);
        tol_min_max
    }

    /// OCCT Z (cxx L291-297) — the viewing-system Z of the parameter point.
    pub fn z(&self, u: f64) -> f64 {
        let p3d = self.c().d0(u);
        let trsf = self.proj().transformation().clone();
        trsf.apply(p3d).z
    }

    /// OCCT Value3D (lxx) / D0 3D (lxx L~180) — the 3D point of parameter U.
    pub fn value_3d(&self, u: f64) -> Point3 {
        self.c().d0(u)
    }

    /// OCCT Tangent (cxx L301-313) — the 2D point and tangent at the start
    /// (or end); a null first derivative falls to the higher-order
    /// significant derivative through the CLProps.
    pub fn tangent(&self, at_start: bool) -> (glam::DVec2, glam::DVec2) {
        let u = if at_start {
            self.c().first_parameter()
        } else {
            self.c().last_parameter()
        };

        let mut p = glam::DVec2::ZERO;
        self.d0_2d(u, &mut p);
        // HLRBRep_CLProps CLP(2, Epsilon(1.)); SetCurve; SetParameter(U).
        let mut clp = super::cl_props::CLProps::new(2, f64::EPSILON);
        clp.set_curve(self);
        clp.set_parameter(u);
        // StdFail_UndefinedDerivative_Raise_if(!CLP.IsTangentDefined()).
        let d = clp.tangent().expect("StdFail_UndefinedDerivative");
        (p, d)
    }

    /// OCCT FirstParameter (lxx).
    pub fn first_parameter(&self) -> f64 {
        self.c().first_parameter()
    }

    /// OCCT GetType (lxx L138-141) — the projected-curve classification.
    pub fn get_type(&self) -> CurveType {
        self.my_type
    }

    /// OCCT LastParameter (lxx).
    pub fn last_parameter(&self) -> f64 {
        self.c().last_parameter()
    }

    /// OCCT IsClosed (lxx).
    pub fn is_closed(&self) -> bool {
        self.c().is_closed()
    }

    /// OCCT IsPeriodic (lxx).
    pub fn is_periodic(&self) -> bool {
        self.c().is_periodic()
    }

    /// OCCT Period (lxx).
    pub fn period(&self) -> f64 {
        self.c().period()
    }

    /// OCCT Value(U) — the 2D projected point (lxx: D0(U, P)).
    pub fn value(&self, u: f64) -> glam::DVec2 {
        let mut p = glam::DVec2::ZERO;
        self.d0_2d(u, &mut p);
        p
    }

    /// OCCT D0(U, P) (cxx L317-330) — Project(P3d, P).
    pub fn d0_2d(&self, u: f64, p: &mut glam::DVec2) {
        let p3d = self.c().d0(u);
        self.proj().project_pnt(p3d, p);
    }

    /// OCCT D1(U, P, V) (cxx L334-363).
    pub fn d1_2d(&self, u: f64, p: &mut glam::DVec2, v: &mut glam::DVec2) {
        let (mut p3d, mut v13d) = self.c().d1(u);
        if self.proj().perspective() {
            let trsf = self.proj().transformation().clone();
            p3d = trsf.apply(p3d);
            v13d = trsf.transform_vec(v13d);

            let f = self.proj().focus();
            let r = 1.0 - p3d.z / f;
            let e = v13d.z / (f * r * r);
            *p = glam::DVec2::new(p3d.x / r, p3d.y / r);
            *v = glam::DVec2::new(v13d.x / r + p3d.x * e, v13d.y / r + p3d.y * e);
        } else {
            // OCC155 — Project(P3D, V13D, P, V).
            self.proj().project_pnt_dir(p3d, v13d, p, v);
        }
    }

    /// OCCT D2(U, P, V1, V2) (cxx L367-402) — the 3D D2 needs the edge
    /// adaptor second derivative; the CurveView carries D2 through
    /// [`CurveView::d2`] (default Other-surface panic is the OCCT
    /// Standard_NoSuchObject behavior on non-C2 curves).
    pub fn d2_2d(&self, u: f64, p: &mut glam::DVec2, v1: &mut glam::DVec2, v2: &mut glam::DVec2) {
        let (mut p3d, mut v13d, mut v23d) = self.c().d2(u);
        let trsf = self.proj().transformation().clone();
        p3d = trsf.apply(p3d);
        v13d = trsf.transform_vec(v13d);
        v23d = trsf.transform_vec(v23d);
        if self.proj().perspective() {
            let f = self.proj().focus();
            let r = 1.0 - p3d.z / f;
            let q = f * r * r;
            let e = v13d.z / q;
            let c = e * v13d.z / (f * r);
            *p = glam::DVec2::new(p3d.x / r, p3d.y / r);
            *v1 = glam::DVec2::new(v13d.x / r + p3d.x * e, v13d.y / r + p3d.y * e);
            *v2 = glam::DVec2::new(
                v23d.x / r + 2.0 * v13d.x * e + p3d.x * v23d.z / q + 2.0 * p3d.x * c,
                v23d.y / r + 2.0 * v13d.y * e + p3d.y * v23d.z / q + 2.0 * p3d.y * c,
            );
        } else {
            *p = glam::DVec2::new(p3d.x, p3d.y);
            *v1 = glam::DVec2::new(v13d.x, v13d.y);
            *v2 = glam::DVec2::new(v23d.x, v23d.y);
        }
    }

    /// OCCT D3 (cxx L406) — empty body verbatim.
    pub fn d3_2d(&self, _u: f64, _p: &mut glam::DVec2, _v1: &mut glam::DVec2, _v2: &mut glam::DVec2, _v3: &mut glam::DVec2) {}

    /// OCCT DN (cxx L410-413) — zero vector verbatim.
    pub fn dn(&self, _u: f64, _n: i32) -> glam::DVec2 {
        glam::DVec2::ZERO
    }

    /// OCCT Line (cxx L417-423).
    pub fn line(&self) -> Line2d {
        let mut p = glam::DVec2::ZERO;
        let mut v = glam::DVec2::ZERO;
        self.d1_2d(0.0, &mut p, &mut v);
        Line2d {
            origin: p,
            direction: v.normalize_or_zero(),
        }
    }

    /// OCCT Circle (cxx L427-432) — ProjLib::Project(XOY, transformed
    /// circle).
    pub fn circle(&self) -> Circle2d {
        let mut c = self.c().circle();
        let trsf = self.proj().transformation().clone();
        c.center = trsf.apply(c.center);
        c.normal = trsf.transform_vec(c.normal);
        c.x_dir = trsf.transform_vec(c.x_dir);
        c.y_dir = trsf.transform_vec(c.y_dir);
        // ProjLib::Project(gp_Pln(gp::XOY()), C): the image on z = 0.
        Circle2d {
            center: glam::DVec2::new(c.center.x, c.center.y),
            x_dir: glam::DVec2::new(c.x_dir.x, c.x_dir.y).normalize_or_zero(),
            y_dir: glam::DVec2::new(c.y_dir.x, c.y_dir.y).normalize_or_zero(),
            radius: c.radius,
        }
    }

    /// OCCT Ellipse (cxx L436-459).
    pub fn ellipse(&self) -> Ellipse2d {
        if self.c().get_type() == CurveType::Ellipse {
            let mut e = self.c().ellipse();
            let trsf = self.proj().transformation().clone();
            e.center = trsf.apply(e.center);
            e.normal = trsf.transform_vec(e.normal);
            e.major_dir = trsf.transform_vec(e.major_dir);
            // ProjLib::Project(gp_Pln(gp::XOY()), E).
            let major_dir = glam::DVec2::new(e.major_dir.x, e.major_dir.y).normalize_or_zero();
            return Ellipse2d {
                center: glam::DVec2::new(e.center.x, e.center.y),
                major_dir,
                major_radius: e.major_radius,
                minor_radius: e.minor_radius,
            };
        }
        // this is a circle
        let mut c = self.c().circle();
        let trsf = self.proj().transformation().clone();
        c.center = trsf.apply(c.center);
        c.normal = trsf.transform_vec(c.normal).normalize();
        let d1 = c.normal;
        let d3 = d1.cross(Vec3::new(0.0, 0.0, 1.0)).normalize();
        let d2 = d1.cross(d3).normalize();
        let rap = (d2.x * d2.x + d2.y * d2.y).sqrt();
        let d = glam::DVec2::new(d1.y, -d1.x);
        let p = glam::DVec2::new(c.center.x, c.center.y);
        let mut el = Ellipse2d {
            center: p,
            major_dir: d.normalize_or_zero(),
            major_radius: c.radius,
            minor_radius: c.radius * rap,
        };
        if d1.z < 0.0 {
            // El.Reverse(): swap the radii around the minor axis — the rcad
            // Ellipse2d keeps the frame; the OCCT reversal flips the major
            // axis direction.
            el.major_dir = -el.major_dir;
        }
        el
    }

    /// OCCT Hyperbola (cxx L463-466) — default-constructed verbatim.
    pub fn hyperbola(&self) -> () {
        // return gp_Hypr2d();
    }

    /// OCCT Parabola (cxx L470-473) — default-constructed verbatim.
    pub fn parabola(&self) -> () {
        // return gp_Parab2d();
    }

    /// OCCT IsRational (lxx) — through the poles' weights; the rcad edge
    /// view carries rationality with the curve data (Stage 3f wiring).
    pub fn is_rational(&self) -> bool {
        false
    }

    /// OCCT Degree (lxx).
    pub fn degree(&self) -> i32 {
        self.c().degree()
    }

    /// OCCT NbPoles (lxx).
    pub fn nb_poles(&self) -> i32 {
        self.c().nb_poles()
    }

    /// OCCT NbKnots (lxx).
    pub fn nb_knots(&self) -> i32 {
        self.c().nb_knots()
    }

    /// OCCT Resolution (lxx).
    pub fn resolution(&self, r3d: f64) -> f64 {
        self.c().resolution(r3d)
    }

    fn proj_xyz(&self, p: Point3) -> (f64, f64, f64) {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        self.proj().project_xyz(p, &mut x, &mut y, &mut z);
        (x, y, z)
    }
}

/// gp_Dir::AngleWithRef(A, V) — the angle from `a` to `other` signed by the
// reference direction `v` (the triple product sign).
fn angle_with_ref(other: Vec3, a: Vec3, v: Vec3) -> f64 {
    let dot = (a.dot(other) / (a.length() * other.length())).clamp(-1.0, 1.0);
    let angle = dot.acos();
    if a.cross(other).dot(v) < 0.0 {
        -angle
    } else {
        angle
    }
}

// The 2D circle shape reference for the projected image.
#[allow(unused)]
fn _shapes(c: Curve2d) -> Curve2d {
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::math::gp::Ax2;

    use crate::hlr::brep::b_curve_tool::CurveView;

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double):
    /// the edge from (0,0,0) to (len,0,0).
    struct LineEdge(f64);

    impl CurveView for LineEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.0
        }
        fn d0(&self, u: f64) -> Point3 {
            Point3::new(u, 0.0, 0.0)
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (Point3::new(u, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0))
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (Point3::new(u, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> rcad_kernel::geom::Line3 {
            rcad_kernel::geom::Line3 {
                origin: Point3::ZERO,
                direction: Vec3::new(1.0, 0.0, 0.0),
            }
        }
        fn circle(&self) -> Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
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

    /// OCCT anchor: the top view of an X-aligned segment — Update classifies
    /// the projection as a line with myVX = 1 (cxx L180-222), the 2D image
    /// of U is (U, 0) (Project), Z = 0, and the 2d/3d parameter maps are
    /// the identity (Parameter2d/3d through myVX, cxx L50-104).
    #[test]
    fn hlr_curve_line_projection() {
        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(|| {
            Projector::from_ax2(&Ax2::new(
                DVec3::ZERO,
                DVec3::new(0.0, 0.0, 1.0),
                DVec3::new(1.0, 0.0, 0.0),
            ))
        });
        static EDGE: std::sync::OnceLock<LineEdge> = std::sync::OnceLock::new();
        let edge: &'static LineEdge = EDGE.get_or_init(|| LineEdge(3.0));

        let mut c = Curve::new();
        c.projector(proj);
        c.load(edge);

        let mut tot_min = [0.0f64; 16];
        let mut tot_max = [0.0f64; 16];
        let tol = c.update(&mut tot_min, &mut tot_max);
        assert_eq!(c.my_type_for_test(), CurveType::Line);
        assert!(tol.abs() < 1e-12);

        // The 2D image: (u, 0).
        let p = c.value(2.0);
        assert!((p.x - 2.0).abs() < 1e-12, "p={:?}", p);
        assert!(p.y.abs() < 1e-12);
        // Z of the parameter point: 0 (in-plane).
        assert!(c.z(2.0).abs() < 1e-12);

        // The parameter maps: identity through myVX = 1.
        assert!((c.parameter_2d(2.0) - 2.0).abs() < 1e-12);
        assert!((c.parameter_3d(2.0) - 2.0).abs() < 1e-12);

        // D1 of the projected line: point (2, 0), tangent (1, 0).
        let (mut p2, mut v2) = (glam::DVec2::ZERO, glam::DVec2::ZERO);
        c.d1_2d(2.0, &mut p2, &mut v2);
        assert!((p2.x - 2.0).abs() < 1e-12 && p2.y.abs() < 1e-12);
        assert!((v2.x.abs() - 1.0).abs() < 1e-12 && v2.y.abs() < 1e-12);

        // Line() rebuilds the projected support.
        let l = c.line();
        assert!((l.direction.x.abs() - 1.0).abs() < 1e-12);
    }

    impl Curve<'_> {
        fn my_type_for_test(&self) -> CurveType {
            self.my_type
        }
    }
}

#[cfg(test)]
mod tangent_tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::math::gp::Ax2;

    use crate::hlr::brep::b_curve_tool::CurveView;

    // Re-export the line fixture shape from the primary test module via a
    // local copy (the CurveView for a segment on X).
    struct LineEdge(f64);

    impl CurveView for LineEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.0
        }
        fn d0(&self, u: f64) -> Point3 {
            Point3::new(u, 0.0, 0.0)
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (Point3::new(u, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0))
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (Point3::new(u, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> rcad_kernel::geom::Line3 {
            rcad_kernel::geom::Line3 {
                origin: Point3::ZERO,
                direction: Vec3::new(1.0, 0.0, 0.0),
            }
        }
        fn circle(&self) -> Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
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

    /// OCCT anchor: Tangent(AtStart) on the projected segment — the point
    /// (0, 0) with tangent (1, 0); Tangent at end — (3, 0), (1, 0)
    /// (cxx L301-313 through HLRBRep_CLProps).
    #[test]
    fn hlr_curve_tangent_at_ends() {
        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(|| {
            Projector::from_ax2(&Ax2::new(
                DVec3::ZERO,
                DVec3::new(0.0, 0.0, 1.0),
                DVec3::new(1.0, 0.0, 0.0),
            ))
        });
        static EDGE: std::sync::OnceLock<LineEdge> = std::sync::OnceLock::new();
        let edge: &'static LineEdge = EDGE.get_or_init(|| LineEdge(3.0));

        let mut c = Curve::new();
        c.projector(proj);
        c.load(edge);

        let (p0, d0) = c.tangent(true);
        assert!((p0.x - 0.0).abs() < 1e-12 && p0.y.abs() < 1e-12);
        assert!((d0.x - 1.0).abs() < 1e-12 && d0.y.abs() < 1e-12);

        let (p1, d1) = c.tangent(false);
        assert!((p1.x - 3.0).abs() < 1e-12 && p1.y.abs() < 1e-12);
        assert!((d1.x - 1.0).abs() < 1e-12 && d1.y.abs() < 1e-12);
    }
}
