//! Adaptor3d_TopolTool (TKG3d) — the default topological tool over a
//! surface domain: restriction lines from the UV bounds, the classification
//! state machine, and the (adaptive) sampling of the surface.
//!
//! 1:1 translation of OCCT `Adaptor3d_TopolTool.hxx` (L17-179) + `.cxx`
//! (L33-1697).  The BRepTopAdaptor_TopolTool subclass lives in
//! `topalgo/brep_top_adaptor/topol_tool.rs`; the identical file-static
//! `Analyse` helper (Adaptor3d_TopolTool.cxx L639-748 =
//! BRepTopAdaptor_TopolTool.cxx L232-341) is shared here as
//! [`analyse`].

use glam::DVec2;
use rcad_kernel::geom::{ConicalSurface, Point3};
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::precision::{
    is_negative_infinite_value, is_positive_infinite_value, INFINITE_VALUE,
};
use rcad_kernel::topods::{Orientation, State};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::int_curve_surface::{HSurfaceTool, SurfaceType};

use super::hvertex::HVertex;
use super::super::adaptor2d::line2d::Line2dAdaptor;

/// OCCT Adaptor3d_TopolTool.cxx L31: `#define myInfinite Precision::Infinite()`.
const MY_INFINITE: f64 = INFINITE_VALUE;
/// OCCT Adaptor3d_TopolTool.cxx L635-637: `#define myMinPnts 4` — absolute
/// possible minimum of sample points (restriction of IntPolyh).
pub(crate) const MY_MIN_PNTS: usize = 4;

/// OCCT Adaptor3d_TopolTool.cxx L639-748 `static void Analyse` (duplicated
/// file-statically in BRepTopAdaptor_TopolTool.cxx L232-341) — the
/// pole-grid curvature-sign-change analysis used to size the sample grid.
/// Returns (myNbSamplesU, myNbSamplesV).
pub(crate) fn analyse(
    array2: &[Vec<Point3>],
    nbup: usize,
    nbvp: usize,
) -> (usize, usize) {
    let mut vip1;
    let mut sh = 1;
    let mut nbch = 0;
    if nbvp > 2 {
        for i in 2..nbup {
            let a = &array2[i][0];
            let b = &array2[i][1];
            let c = &array2[i][2];
            let mut v3 = c - b - b + a;
            let mut locnbch = 0;
            for j in 3..nbvp {
                let a1 = &array2[i][j - 1];
                let b1 = &array2[i][j];
                let c1 = &array2[i][j + 1];
                vip1 = c1 - b1 - b1 + a1;
                let pd = v3.dot(vip1);
                v3 = vip1;
                if pd > 1.0e-7 || pd < -1.0e-7 {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            locnbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        locnbch += 1;
                    }
                }
            }
            if locnbch > nbch {
                nbch = locnbch;
            }
        }
    }
    let my_nb_samples_v = nbch + 5;

    nbch = 0;
    if nbup > 2 {
        for j in 2..nbvp {
            let a = &array2[0][j];
            let b = &array2[1][j];
            let c = &array2[2][j];
            let mut v3 = c - b - b + a;
            let mut locnbch = 0;
            for i in 3..nbup {
                let a1 = &array2[i - 1][j];
                let b1 = &array2[i][j];
                let c1 = &array2[i + 1][j];
                vip1 = c1 - b1 - b1 + a1;
                let pd = v3.dot(vip1);
                v3 = vip1;
                if pd > 1.0e-7 || pd < -1.0e-7 {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            locnbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        locnbch += 1;
                    }
                }
            }
            if locnbch > nbch {
                nbch = locnbch;
            }
        }
    }
    (nbch + 5, my_nb_samples_v)
}

/// OCCT Adaptor3d_TopolTool — "a default topological tool, based on the
/// Umin,Vmin,Umax,Vmax of an HSurface from Adaptor3d".
pub struct TopolTool<'a, S, ST: HSurfaceTool<Surface = S>> {
    /// OCCT myS.
    pub(crate) my_s: &'a S,
    /// OCCT myNbSamplesU / myNbSamplesV.
    pub(crate) my_nb_samples_u: i32,
    pub(crate) my_nb_samples_v: i32,
    /// OCCT myUPars / myVPars (HArray1 contents).
    pub(crate) my_u_pars: Option<Vec<f64>>,
    pub(crate) my_v_pars: Option<Vec<f64>>,
    // private:
    nb_restrict: usize,
    id_restrict: usize,
    u_inf: f64,
    u_sup: f64,
    v_inf: f64,
    v_sup: f64,
    /// OCCT myRestr[4].
    my_restrict: Vec<Line2dAdaptor>,
    nb_vtx: usize,
    id_vtx: usize,
    /// OCCT myVtx[2].
    my_vtx: Vec<HVertex>,
    _tool: std::marker::PhantomData<fn(&ST)>,
}

impl<'a, S, ST: HSurfaceTool<Surface = S>> TopolTool<'a, S, ST> {
    /// OCCT Adaptor3d_TopolTool() (cxx L33-45).
    pub fn new(my_s: &'a S) -> Self {
        let mut t = TopolTool {
            my_s,
            my_nb_samples_u: -1,
            my_nb_samples_v: -1,
            my_u_pars: None,
            my_v_pars: None,
            nb_restrict: 0,
            id_restrict: 0,
            u_inf: 0.0,
            u_sup: 0.0,
            v_inf: 0.0,
            v_sup: 0.0,
            my_restrict: Vec::new(),
            nb_vtx: 0,
            id_vtx: 0,
            my_vtx: Vec::new(),
            _tool: std::marker::PhantomData,
        };
        t.initialize_surface(my_s);
        t
    }

    /// OCCT Initialize() (cxx L52-55) — the base class throws
    /// Standard_NotImplemented.
    pub fn initialize() -> ! {
        panic!("Standard_NotImplemented: Adaptor3d_TopolTool::Initialize ()");
    }

    /// OCCT Initialize(S) (cxx L57-209) — build the restriction lines from
    /// the finite domain bounds (plus the cone apex line).
    pub fn initialize_surface(&mut self, s: &'a S) {
        let mut pinf;
        let mut psup;
        let mut deltap;

        self.my_nb_samples_u = -1;
        self.u_inf = <ST as HSurfaceTool>::first_u_parameter(s);
        self.v_inf = <ST as HSurfaceTool>::first_v_parameter(s);
        self.u_sup = <ST as HSurfaceTool>::last_u_parameter(s);
        self.v_sup = <ST as HSurfaceTool>::last_v_parameter(s);
        self.nb_restrict = 0;
        self.id_restrict = 0;
        self.my_restrict.clear();

        let u_inf_infinite = is_negative_infinite_value(self.u_inf);
        let u_sup_infinite = is_positive_infinite_value(self.u_sup);
        let v_inf_infinite = is_negative_infinite_value(self.v_inf);
        let v_sup_infinite = is_positive_infinite_value(self.v_sup);

        if !v_inf_infinite {
            deltap = (self.u_sup - self.u_inf).min(2.0 * MY_INFINITE);
            if self.u_inf >= -MY_INFINITE {
                pinf = self.u_inf;
                psup = pinf + deltap;
            } else if self.u_sup <= MY_INFINITE {
                psup = self.u_sup;
                pinf = psup - deltap;
            } else {
                pinf = -MY_INFINITE;
                psup = MY_INFINITE;
            }

            self.my_restrict.push(Line2dAdaptor::new_pnt_dir(
                DVec2::new(0.0, self.v_inf),
                DVec2::X,
                pinf,
                psup,
            ));
            self.nb_restrict += 1;
        }

        if !u_sup_infinite {
            deltap = (self.v_sup - self.v_inf).min(2.0 * MY_INFINITE);
            if self.v_inf >= -MY_INFINITE {
                pinf = self.v_inf;
                psup = pinf + deltap;
            } else if self.v_sup <= MY_INFINITE {
                psup = self.v_sup;
                pinf = psup - deltap;
            } else {
                pinf = -MY_INFINITE;
                psup = MY_INFINITE;
            }

            self.my_restrict.push(Line2dAdaptor::new_pnt_dir(
                DVec2::new(self.u_sup, 0.0),
                DVec2::Y,
                pinf,
                psup,
            ));
            self.nb_restrict += 1;
        }

        if !v_sup_infinite {
            deltap = (self.u_sup - self.u_inf).min(2.0 * MY_INFINITE);
            if -self.u_sup >= -MY_INFINITE {
                pinf = -self.u_sup;
                psup = pinf + deltap;
            } else if -self.u_inf <= MY_INFINITE {
                psup = -self.u_inf;
                pinf = psup - deltap;
            } else {
                pinf = -MY_INFINITE;
                psup = MY_INFINITE;
            }

            self.my_restrict.push(Line2dAdaptor::new_pnt_dir(
                DVec2::new(0.0, self.v_sup),
                -DVec2::X,
                pinf,
                psup,
            ));
            self.nb_restrict += 1;
        }

        if !u_inf_infinite {
            deltap = (self.v_sup - self.v_inf).min(2.0 * MY_INFINITE);
            if -self.v_sup >= -MY_INFINITE {
                pinf = -self.v_sup;
                psup = pinf + deltap;
            } else if -self.v_inf <= MY_INFINITE {
                psup = -self.v_inf;
                pinf = psup - deltap;
            } else {
                pinf = -MY_INFINITE;
                psup = MY_INFINITE;
            }

            self.my_restrict.push(Line2dAdaptor::new_pnt_dir(
                DVec2::new(self.u_inf, 0.0),
                -DVec2::Y,
                pinf,
                psup,
            ));
            self.nb_restrict += 1;
        }

        self.my_s = s;

        if self.nb_restrict == 2 && <ST as HSurfaceTool>::get_type(s) == SurfaceType::Cone {
            let (mut u, mut v) = (0.0f64, 0.0f64);
            get_cone_apex_param(&<ST as HSurfaceTool>::cone(s), &mut u, &mut v);

            deltap = (self.u_sup - self.u_inf).min(2.0 * MY_INFINITE);
            if self.u_inf >= -MY_INFINITE {
                pinf = self.u_inf;
                psup = pinf + deltap;
            } else if self.u_sup <= MY_INFINITE {
                psup = self.u_sup;
                pinf = psup - deltap;
            } else {
                pinf = -MY_INFINITE;
                psup = MY_INFINITE;
            }

            self.my_restrict.push(Line2dAdaptor::new_pnt_dir(
                DVec2::new(u, v),
                DVec2::X,
                pinf,
                psup,
            ));
            self.nb_restrict += 1;
        }
    }

    /// OCCT Initialize(C) (cxx L235-254) — the arc's end vertices.
    pub fn initialize_curve(&mut self, c: &dyn Curve2dAdaptor) {
        self.nb_vtx = 0;
        self.id_vtx = 0;
        self.my_vtx.clear();
        let the_uinf = c.first_parameter();
        let the_usup = c.last_parameter();
        if the_uinf > -MY_INFINITE {
            self.my_vtx
                .push(HVertex::new_with(c.value(the_uinf), Orientation::Forward, 1.0e-8));
            self.nb_vtx += 1;
        }
        if the_usup < MY_INFINITE {
            self.my_vtx.push(HVertex::new_with(
                c.value(the_usup),
                Orientation::Reversed,
                1.0e-8,
            ));
            self.nb_vtx += 1;
        }
    }

    /// OCCT Init() (cxx L211-214).
    pub fn init(&mut self) {
        self.id_restrict = 0;
    }

    /// OCCT More() (cxx L216-219).
    pub fn more(&self) -> bool {
        self.id_restrict < self.nb_restrict
    }

    /// OCCT Value() (cxx L221-228) — the current restriction line.
    pub fn value(&self) -> &Line2dAdaptor {
        if self.id_restrict >= self.nb_restrict {
            panic!("Standard_DomainError");
        }
        &self.my_restrict[self.id_restrict]
    }

    /// OCCT Next() (cxx L230-233).
    pub fn next(&mut self) {
        self.id_restrict += 1;
    }

    /// OCCT InitVertexIterator() (cxx L256-259).
    pub fn init_vertex_iterator(&mut self) {
        self.id_vtx = 0;
    }

    /// OCCT MoreVertex() (cxx L261-264).
    pub fn more_vertex(&self) -> bool {
        self.id_vtx < self.nb_vtx
    }

    /// OCCT Vertex() (cxx L266-273).
    pub fn vertex(&self) -> &HVertex {
        if self.id_vtx >= self.nb_vtx {
            panic!("Standard_DomainError");
        }
        &self.my_vtx[self.id_vtx]
    }

    /// OCCT NextVertex() (cxx L275-278).
    pub fn next_vertex(&mut self) {
        self.id_vtx += 1;
    }

    /// OCCT Classify(P, Tol, RecadreOnPeriodic) (cxx L280-441) — the domain
    /// state machine (the RecadreOnPeriodic argument is unused in the base
    /// class).
    pub fn classify(&self, p: DVec2, tol: f64, _recadre_on_periodic: bool) -> State {
        let u = p.x;
        let v = p.y;

        if self.nb_restrict == 4 {
            if (u < self.u_inf - tol)
                || (u > self.u_sup + tol)
                || (v < self.v_inf - tol)
                || (v > self.v_sup + tol)
            {
                return State::Out;
            }
            if (u - self.u_inf).abs() <= tol
                || (u - self.u_sup).abs() <= tol
                || (v - self.v_inf).abs() <= tol
                || (v - self.v_sup).abs() <= tol
            {
                return State::On;
            }
            State::In
        } else if self.nb_restrict == 0 {
            State::In
        } else {
            let (dansu, surumin, surumax);
            if is_negative_infinite_value(self.u_inf) && is_positive_infinite_value(self.u_sup) {
                dansu = true;
                surumin = false;
                surumax = false;
            } else if is_negative_infinite_value(self.u_inf) {
                surumin = false;
                if u >= self.u_sup + tol {
                    dansu = false;
                    surumax = false;
                } else {
                    dansu = true;
                    surumax = (u - self.u_sup).abs() <= tol;
                }
            } else if is_positive_infinite_value(self.u_sup) {
                surumax = false;
                if u < self.u_inf - tol {
                    dansu = false;
                    surumin = false;
                } else {
                    dansu = true;
                    surumin = (u - self.u_inf).abs() <= tol;
                }
            } else if (u < self.u_inf - tol) || (u > self.u_sup + tol) {
                dansu = false;
                surumin = false;
                surumax = false;
            } else {
                dansu = true;
                surumin = (u - self.u_inf).abs() <= tol;
                surumax = if (u - self.u_inf).abs() <= tol {
                    false
                } else {
                    (u - self.u_sup).abs() <= tol
                };
            }

            let (dansv, survmin, survmax);
            if is_negative_infinite_value(self.v_inf) && is_positive_infinite_value(self.v_sup) {
                dansv = true;
                survmin = false;
                survmax = false;
            } else if is_negative_infinite_value(self.v_inf) {
                survmin = false;
                if v >= self.v_sup + tol {
                    dansv = false;
                    survmax = false;
                } else {
                    dansv = true;
                    survmax = (v - self.v_sup).abs() <= tol;
                }
            } else if is_positive_infinite_value(self.v_sup) {
                survmax = false;
                if v < self.v_inf - tol {
                    dansv = false;
                    survmin = false;
                } else {
                    dansv = true;
                    survmin = (v - self.v_inf).abs() <= tol;
                }
            } else if (v < self.v_inf - tol) || (v > self.v_sup + tol) {
                dansv = false;
                survmin = false;
                survmax = false;
            } else {
                dansv = true;
                survmin = (v - self.v_inf).abs() <= tol;
                survmax = if (v - self.v_inf).abs() <= tol {
                    false
                } else {
                    (v - self.v_sup).abs() <= tol
                };
            }

            if !dansu || !dansv {
                return State::Out;
            }
            if surumin || survmin || surumax || survmax {
                return State::On;
            }
            State::In
        }
    }

    /// OCCT IsThePointOn(P, Tol, RecadreOnPeriodic) (cxx L443-604) — true
    /// when the point is on the domain restriction.
    pub fn is_the_point_on(&self, p: DVec2, tol: f64, _recadre_on_periodic: bool) -> bool {
        let u = p.x;
        let v = p.y;

        if self.nb_restrict == 4 {
            if (u < self.u_inf - tol)
                || (u > self.u_sup + tol)
                || (v < self.v_inf - tol)
                || (v > self.v_sup + tol)
            {
                return false;
            }
            if (u - self.u_inf).abs() <= tol
                || (u - self.u_sup).abs() <= tol
                || (v - self.v_inf).abs() <= tol
                || (v - self.v_sup).abs() <= tol
            {
                return true;
            }
            return false;
        } else if self.nb_restrict == 0 {
            return false;
        } else {
            let (dansu, surumin, surumax);
            if is_negative_infinite_value(self.u_inf) && is_positive_infinite_value(self.u_sup) {
                dansu = true;
                surumin = false;
                surumax = false;
            } else if is_negative_infinite_value(self.u_inf) {
                surumin = false;
                if u >= self.u_sup + tol {
                    dansu = false;
                    surumax = false;
                } else {
                    dansu = true;
                    surumax = (u - self.u_sup).abs() <= tol;
                }
            } else if is_positive_infinite_value(self.u_sup) {
                surumax = false;
                if u < self.u_inf - tol {
                    dansu = false;
                    surumin = false;
                } else {
                    dansu = true;
                    surumin = (u - self.u_inf).abs() <= tol;
                }
            } else if (u < self.u_inf - tol) || (u > self.u_sup + tol) {
                dansu = false;
                surumin = false;
                surumax = false;
            } else {
                dansu = true;
                surumin = (u - self.u_inf).abs() <= tol;
                surumax = if (u - self.u_inf).abs() <= tol {
                    false
                } else {
                    (u - self.u_sup).abs() <= tol
                };
            }

            let (dansv, survmin, survmax);
            if is_negative_infinite_value(self.v_inf) && is_positive_infinite_value(self.v_sup) {
                dansv = true;
                survmin = false;
                survmax = false;
            } else if is_negative_infinite_value(self.v_inf) {
                survmin = false;
                if v >= self.v_sup + tol {
                    dansv = false;
                    survmax = false;
                } else {
                    dansv = true;
                    survmax = (v - self.v_sup).abs() <= tol;
                }
            } else if is_positive_infinite_value(self.v_sup) {
                survmax = false;
                if v < self.v_inf - tol {
                    dansv = false;
                    survmin = false;
                } else {
                    dansv = true;
                    survmin = (v - self.v_inf).abs() <= tol;
                }
            } else if (v < self.v_inf - tol) || (v > self.v_sup + tol) {
                dansv = false;
                survmin = false;
                survmax = false;
            } else {
                dansv = true;
                survmin = (v - self.v_inf).abs() <= tol;
                survmax = if (v - self.v_inf).abs() <= tol {
                    false
                } else {
                    (v - self.v_sup).abs() <= tol
                };
            }

            if !dansu || !dansv {
                return false;
            }
            if surumin || survmin || surumax || survmax {
                return true;
            }
            false
        }
    }

    /// OCCT Orientation(C) (cxx L606-609) — the base restriction arcs are
    /// FORWARD.
    pub fn orientation_curve(&self, _c: &dyn Curve2dAdaptor) -> Orientation {
        Orientation::Forward
    }

    /// OCCT Orientation(V) (cxx L611-614).
    pub fn orientation_vertex(&self, v: &HVertex) -> Orientation {
        v.orientation()
    }

    /// OCCT Identical(V1, V2) (cxx L616-620).
    pub fn identical(&self, v1: &HVertex, v2: &HVertex) -> bool {
        v1.is_same(v2)
    }

    /// OCCT Has3d() (cxx L991-994).
    pub fn has_3d(&self) -> bool {
        false
    }

    /// OCCT Tol3d(C) (cxx L998-1001).
    pub fn tol3d_curve(&self, _c: &dyn Curve2dAdaptor) -> f64 {
        panic!("Standard_DomainError: Adaptor3d_TopolTool: has no 3d representation");
    }

    /// OCCT Tol3d(V) (cxx L1005-1008).
    pub fn tol3d_vertex(&self, _v: &HVertex) -> f64 {
        panic!("Standard_DomainError: Adaptor3d_TopolTool: has no 3d representation");
    }

    /// OCCT Pnt(V) (cxx L1012-1015).
    pub fn pnt(&self, _v: &HVertex) -> Point3 {
        panic!("Standard_DomainError: Adaptor3d_TopolTool: has no 3d representation");
    }

    /// OCCT ComputeSamplePoints() (cxx L750-897).
    pub fn compute_sample_points(&mut self) {
        const A_MAX_NB_SAMPLE: i32 = 50;

        let mut uinf = <ST as HSurfaceTool>::first_u_parameter(self.my_s);
        let mut usup = <ST as HSurfaceTool>::last_u_parameter(self.my_s);
        let mut vinf = <ST as HSurfaceTool>::first_v_parameter(self.my_s);
        let mut vsup = <ST as HSurfaceTool>::last_v_parameter(self.my_s);
        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }
        if uinf == f64::MIN && usup == f64::MAX {
            uinf = -1.0e5;
            usup = 1.0e5;
        } else if uinf == f64::MIN {
            uinf = usup - 2.0e5;
        } else if usup == f64::MAX {
            usup = uinf + 2.0e5;
        }

        if vinf == f64::MIN && vsup == f64::MAX {
            vinf = -1.0e5;
            vsup = 1.0e5;
        } else if vinf == f64::MIN {
            vinf = vsup - 2.0e5;
        } else if vsup == f64::MAX {
            vsup = vinf + 2.0e5;
        }

        let typ_s = <ST as HSurfaceTool>::get_type(self.my_s);
        let (mut nbsu, mut nbsv) = match typ_s {
            SurfaceType::Plane => (2usize, 2usize),
            SurfaceType::BezierSurface => (
                3 + <ST as HSurfaceTool>::nb_u_poles(self.my_s),
                3 + <ST as HSurfaceTool>::nb_v_poles(self.my_s),
            ),
            SurfaceType::BSplineSurface => {
                let mut nbsv2 = <ST as HSurfaceTool>::nb_v_knots(self.my_s)
                    * <ST as HSurfaceTool>::v_degree(self.my_s);
                if nbsv2 < 4 {
                    nbsv2 = 4;
                }
                let mut nbsu2 = <ST as HSurfaceTool>::nb_u_knots(self.my_s)
                    * <ST as HSurfaceTool>::u_degree(self.my_s);
                if nbsu2 < 4 {
                    nbsu2 = 4;
                }
                (nbsu2, nbsv2)
            }
            SurfaceType::Cylinder
            | SurfaceType::Cone
            | SurfaceType::Sphere
            | SurfaceType::Torus
            | SurfaceType::SurfaceOfRevolution
            | SurfaceType::SurfaceOfExtrusion => (15usize, 15usize),
            _ => (10usize, 10usize),
        };

        //-- If the number of points is too great... analyze
        if nbsu < 6 {
            nbsu = 6;
        }
        if nbsv < 6 {
            nbsv = 6;
        }

        if typ_s == SurfaceType::BSplineSurface {
            if nbsu > 8 || nbsv > 8 {
                let bspl = <ST as HSurfaceTool>::bspline(self.my_s);
                let array2 = &bspl.control_points;
                let nbup = array2.len();
                let nbvp = array2.first().map_or(0, |r| r.len());
                let (u2, v2) = analyse(array2, nbup, nbvp);
                nbsu = u2;
                nbsv = v2;
            }
            // Check anisotropy
            let an_u_len = (usup - uinf) / <ST as HSurfaceTool>::u_resolution(self.my_s, 1.0);
            let an_v_len = (vsup - vinf) / <ST as HSurfaceTool>::v_resolution(self.my_s, 1.0);
            let a_ratio = an_u_len / an_v_len;
            if a_ratio >= 10.0 {
                nbsu = (nbsu as i32 * 2).min(A_MAX_NB_SAMPLE) as usize;
            } else if a_ratio <= 0.1 {
                nbsv = (nbsv as i32 * 2).min(A_MAX_NB_SAMPLE) as usize;
            }
        } else if typ_s == SurfaceType::BezierSurface && (nbsu > 8 || nbsv > 8) {
            let bez = <ST as HSurfaceTool>::bezier(self.my_s);
            let array2 = &bez.control_points;
            let nbup = array2.len();
            let nbvp = array2.first().map_or(0, |r| r.len());
            let (u2, v2) = analyse(array2, nbup, nbvp);
            nbsu = u2;
            nbsv = v2;
        }

        self.my_nb_samples_u = nbsu as i32;
        self.my_nb_samples_v = nbsv as i32;
    }

    /// OCCT NbSamplesU() (cxx L899-906).
    pub fn nb_samples_u(&mut self) -> i32 {
        if self.my_nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.my_nb_samples_u
    }

    /// OCCT NbSamplesV() (cxx L908-915).
    pub fn nb_samples_v(&mut self) -> i32 {
        if self.my_nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.my_nb_samples_v
    }

    /// OCCT NbSamples() (cxx L917-924).
    pub fn nb_samples(&mut self) -> i32 {
        if self.my_nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.my_nb_samples_u * self.my_nb_samples_v
    }

    /// OCCT UParameters(theArray) (cxx L926-929) — the sample U parameters
    /// set by SamplePnts.
    pub fn u_parameters(&self, the_array: &mut [f64]) {
        let pars = self.my_u_pars.as_ref().expect("myUPars is null");
        the_array.copy_from_slice(&pars[..the_array.len()]);
    }

    /// OCCT VParameters(theArray) (cxx L931-934).
    pub fn v_parameters(&self, the_array: &mut [f64]) {
        let pars = self.my_v_pars.as_ref().expect("myVPars is null");
        the_array.copy_from_slice(&pars[..the_array.len()]);
    }

    /// OCCT SamplePoint(i, P2d, P3d) (cxx L936-959).
    pub fn sample_point(&self, i: usize) -> (DVec2, Point3) {
        let (u, v);
        if self.my_u_pars.is_none() {
            let my_du = (self.u_sup - self.u_inf) / (self.my_nb_samples_u + 1) as f64;
            let my_dv = (self.v_sup - self.v_inf) / (self.my_nb_samples_v + 1) as f64;
            let iv = 1 + i as i32 / self.my_nb_samples_u;
            let iu = 1 + i as i32 - (iv - 1) * self.my_nb_samples_u;
            u = self.u_inf + iu as f64 * my_du;
            v = self.v_inf + iv as f64 * my_dv;
        } else {
            let u_pars = self.my_u_pars.as_ref().unwrap();
            let v_pars = self.my_v_pars.as_ref().unwrap();
            let iv = (i as i32 - 1) / self.my_nb_samples_u + 1;
            let iu = (i as i32 - 1) % self.my_nb_samples_u + 1;
            u = u_pars[(iu - 1) as usize];
            v = v_pars[(iv - 1) as usize];
        }

        (DVec2::new(u, v), <ST as HSurfaceTool>::d0(self.my_s, u, v))
    }

    /// OCCT DomainIsInfinite() (cxx L961-980).
    pub fn domain_is_infinite(&self) -> bool {
        if is_negative_infinite_value(self.u_inf) {
            return true;
        }
        if is_positive_infinite_value(self.u_sup) {
            return true;
        }
        if is_negative_infinite_value(self.v_inf) {
            return true;
        }
        is_positive_infinite_value(self.v_sup)
    }

    /// OCCT Edge() (cxx L984-987) — the base tool has no BRep edge.
    pub fn edge(&self) -> Option<()> {
        None
    }

    /// OCCT SamplePnts(theDefl, theNUmin, theNVmin) (cxx L1019-1130).
    pub fn sample_pnts(&mut self, the_defl: f64, the_numin: i32, the_nvmin: i32) {
        let mut uinf = <ST as HSurfaceTool>::first_u_parameter(self.my_s);
        let mut usup = <ST as HSurfaceTool>::last_u_parameter(self.my_s);
        let mut vinf = <ST as HSurfaceTool>::first_v_parameter(self.my_s);
        let mut vsup = <ST as HSurfaceTool>::last_v_parameter(self.my_s);
        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }
        if uinf == f64::MIN && usup == f64::MAX {
            uinf = -1.0e5;
            usup = 1.0e5;
        } else if uinf == f64::MIN {
            uinf = usup - 2.0e5;
        } else if usup == f64::MAX {
            usup = uinf + 2.0e5;
        }

        if vinf == f64::MIN && vsup == f64::MAX {
            vinf = -1.0e5;
            vsup = 1.0e5;
        } else if vinf == f64::MIN {
            vinf = vsup - 2.0e5;
        } else if vsup == f64::MAX {
            vsup = vinf + 2.0e5;
        }

        let typ_s = <ST as HSurfaceTool>::get_type(self.my_s);
        if typ_s == SurfaceType::BSplineSurface {
            // Processing BSpline surface
            self.bspl_sample_pnts(the_defl, the_numin, the_nvmin);
            return;
        }
        self.compute_sample_points();

        let n = self.my_nb_samples_u as usize;
        let mut my_u_pars = vec![0.0f64; n];
        let mut dt = (usup - uinf) / (self.my_nb_samples_u - 1) as f64;
        my_u_pars[0] = uinf;
        my_u_pars[n - 1] = usup;
        let mut t = uinf + dt;
        for p in my_u_pars.iter_mut().skip(1).take(n - 2) {
            *p = t;
            t += dt;
        }

        let nv = self.my_nb_samples_v as usize;
        let mut my_v_pars = vec![0.0f64; nv];
        dt = (vsup - vinf) / (self.my_nb_samples_v - 1) as f64;
        my_v_pars[0] = vinf;
        my_v_pars[nv - 1] = vsup;
        let mut t = vinf + dt;
        for p in my_v_pars.iter_mut().skip(1).take(nv - 2) {
            *p = t;
            t += dt;
        }

        self.my_u_pars = Some(my_u_pars);
        self.my_v_pars = Some(my_v_pars);
    }

    /// OCCT BSplSamplePnts(theDefl, theNUmin, theNVmin) (cxx L1134-1645).
    pub fn bspl_sample_pnts(&mut self, the_defl: f64, the_numin: i32, the_nvmin: i32) {
        const A_MAX_PNTS: i32 = 1001;
        let a_bs = <ST as HSurfaceTool>::bspline(self.my_s);
        let uinf_t = <ST as HSurfaceTool>::first_u_parameter(self.my_s);
        let usup_t = <ST as HSurfaceTool>::last_u_parameter(self.my_s);
        let vinf_t = <ST as HSurfaceTool>::first_v_parameter(self.my_s);
        let vsup_t = <ST as HSurfaceTool>::last_v_parameter(self.my_s);

        let u_knots = crate::geomalgo::int_curve_surface::distinct_knots(&a_bs.knots_u);
        let v_knots = crate::geomalgo::int_curve_surface::distinct_knots(&a_bs.knots_v);

        let mut j = 1usize;
        let mut ui1 = 1usize;
        let mut ui2 = u_knots.len();
        let mut vi1 = 1usize;
        let mut vi2 = v_knots.len();

        for i in ui1..ui2 {
            if uinf_t >= u_knots[i - 1] && uinf_t < u_knots[i] {
                ui1 = i;
                break;
            }
        }

        for i in (ui1 + 1..=ui2).rev() {
            if usup_t <= u_knots[i - 1] && usup_t > u_knots[i - 2] {
                ui2 = i;
                break;
            }
        }

        for i in vi1..vi2 {
            if vinf_t >= v_knots[i - 1] && vinf_t < v_knots[i] {
                vi1 = i;
                break;
            }
        }

        for i in (vi1 + 1..=vi2).rev() {
            if vsup_t <= v_knots[i - 1] && vsup_t > v_knots[i - 2] {
                vi2 = i;
                break;
            }
        }

        let mut nbsu = ui2 - ui1 + 1;
        nbsu += (nbsu - 1) * (a_bs.degree_u.saturating_sub(1));
        let mut nbsv = vi2 - vi1 + 1;
        nbsv += (nbsv - 1) * (a_bs.degree_v.saturating_sub(1));
        let mut b_uuniform = false;
        let mut b_vuniform = false;

        // modified by NIZHNY-EMV Mon Jun 10 14:19:04 2013
        if (nbsu as i32) < the_numin || (nbsv as i32) < the_nvmin {
            if nbsu < nbsv {
                let mut a_nb = ((nbsv as f64) * (the_numin as f64) / (nbsu as f64)) as i32;
                a_nb = a_nb.min(30);
                b_vuniform = (a_nb as usize > nbsv) || b_vuniform;
                if b_vuniform {
                    nbsv = a_nb as usize;
                }
            } else {
                let mut a_nb = ((nbsu as f64) * (the_nvmin as f64) / (nbsv as f64)) as i32;
                a_nb = a_nb.min(30);
                b_uuniform = (a_nb as usize > nbsu) || b_uuniform;
                if b_uuniform {
                    nbsu = a_nb as usize;
                }
            }
        }
        // modified by NIZHNY-EMV Mon Jun 10 14:19:05 2013

        if (nbsu as i32) < the_numin {
            nbsu = the_numin as usize;
            b_uuniform = true;
        } else if nbsu as i32 > A_MAX_PNTS {
            nbsu = A_MAX_PNTS as usize;
            b_uuniform = true;
        }
        if (nbsv as i32) < the_nvmin {
            nbsv = the_nvmin as usize;
            b_vuniform = true;
        } else if nbsv as i32 > A_MAX_PNTS {
            nbsv = A_MAX_PNTS as usize;
            b_vuniform = true;
        }

        let mut an_u_pars = vec![0.0f64; nbsu];
        let mut an_u_flg = vec![false; nbsu];
        let mut a_v_pars = vec![0.0f64; nbsv];
        let mut a_v_flg = vec![false; nbsv];

        // Filling of sample parameters
        if b_uuniform {
            let mut t1 = uinf_t;
            let t2 = usup_t;
            let dt = (t2 - t1) / (nbsu - 1) as f64;
            an_u_pars[0] = t1;
            an_u_flg[0] = false;
            an_u_pars[nbsu - 1] = t2;
            an_u_flg[nbsu - 1] = false;
            for p in an_u_pars.iter_mut().skip(1).take(nbsu - 2) {
                t1 += dt;
                *p = t1;
            }
        } else {
            let nbi = a_bs.degree_u;
            let mut k = 0usize;
            let mut t1 = uinf_t;
            for i in ui1 + 1..=ui2 {
                let t2 = if i == ui2 { usup_t } else { u_knots[i - 1] };
                let dt = (t2 - t1) / nbi as f64;
                j = 1;
                loop {
                    k += 1;
                    an_u_pars[k - 1] = t1;
                    an_u_flg[k - 1] = false;
                    t1 += dt;
                    j += 1;
                    if j > nbi {
                        break;
                    }
                }
                t1 = t2;
            }
            k += 1;
            an_u_pars[k - 1] = t1;
        }

        if b_vuniform {
            let mut t1 = vinf_t;
            let t2 = vsup_t;
            let dt = (t2 - t1) / (nbsv - 1) as f64;
            a_v_pars[0] = t1;
            a_v_flg[0] = false;
            a_v_pars[nbsv - 1] = t2;
            a_v_flg[nbsv - 1] = false;
            for p in a_v_pars.iter_mut().skip(1).take(nbsv - 2) {
                t1 += dt;
                *p = t1;
            }
        } else {
            let nbi = a_bs.degree_v;
            let mut k = 0usize;
            let mut t1 = vinf_t;
            for i in vi1 + 1..=vi2 {
                let t2 = if i == vi2 { vsup_t } else { v_knots[i - 1] };
                let dt = (t2 - t1) / nbi as f64;
                j = 1;
                loop {
                    k += 1;
                    a_v_pars[k - 1] = t1;
                    a_v_flg[k - 1] = false;
                    t1 += dt;
                    j += 1;
                    if j > nbi {
                        break;
                    }
                }
                t1 = t2;
            }
            k += 1;
            a_v_pars[k - 1] = t1;
        }

        // Analysis of deflection
        let a_defl2 = (the_defl * the_defl).max(1.0e-9);
        let tol = (0.01 * a_defl2).max(1.0e-9);

        an_u_flg[0] = true;
        an_u_flg[nbsu - 1] = true;
        for i in 1..=nbsv {
            let t1 = a_v_pars[i - 1];
            let mut jj = 1usize;
            let mut b_cont = true;
            while jj < nbsu - 1 && b_cont {
                if an_u_flg[jj] {
                    jj += 1;
                    continue;
                }

                let mut t2 = an_u_pars[jj - 1];
                let p1 = <ST as HSurfaceTool>::d0(self.my_s, t2, t1);
                let mut broke = false;
                for k in jj + 2..=nbsu {
                    t2 = an_u_pars[k - 1];
                    let p2 = <ST as HSurfaceTool>::d0(self.my_s, t2, t1);

                    if p1.distance_squared(p2) <= tol {
                        continue;
                    }

                    // gp_Lin lin(p1, gp_Dir(gp_Vec(p1, p2)));
                    let d = (p2 - p1).normalize_or_zero();
                    let mut ok = true;
                    for l in jj + 1..k {
                        if an_u_flg[l - 1] {
                            ok = false;
                            break;
                        }

                        let pp = <ST as HSurfaceTool>::d0(self.my_s, an_u_pars[l - 1], t1);
                        let dist = d.cross(pp - p1).length_squared();

                        if dist <= a_defl2 {
                            continue;
                        }

                        ok = false;
                        break;
                    }

                    if !ok {
                        jj = k - 1;
                        an_u_flg[jj - 1] = true;
                        broke = true;
                        break;
                    }

                    if an_u_flg[k - 1] {
                        jj = k;
                        broke = true;
                        break;
                    }
                }

                if !broke {
                    b_cont = false;
                }
            }
        }

        let mut my_nb_samples_u = an_u_flg.iter().filter(|&&f| f).count();

        if my_nb_samples_u < MY_MIN_PNTS {
            if my_nb_samples_u == 2 {
                // "uniform" distribution
                let nn = nbsu / MY_MIN_PNTS;
                an_u_flg[nn] = true;
                an_u_flg[nbsu - 1 - nn] = true;
            } else {
                // myNbSamplesU == 3: insert in bigger segment
                let mut i = 2usize;
                while !an_u_flg[i - 1] {
                    i += 1;
                }
                let j = if i < nbsu / 2 {
                    (i + (nbsu - i) / 2).min(nbsu - 1)
                } else {
                    (i / 2).max(2)
                };
                an_u_flg[j - 1] = true;
            }
            my_nb_samples_u = MY_MIN_PNTS;
        }

        a_v_flg[0] = true;
        a_v_flg[nbsv - 1] = true;
        for i in 1..=nbsu {
            let t1 = an_u_pars[i - 1];
            let mut jj = 1usize;
            let mut b_cont = true;
            while jj < nbsv - 1 && b_cont {
                if a_v_flg[jj] {
                    jj += 1;
                    continue;
                }

                let mut t2 = a_v_pars[jj - 1];
                let p1 = <ST as HSurfaceTool>::d0(self.my_s, t1, t2);
                let mut broke = false;
                for k in jj + 2..=nbsv {
                    t2 = a_v_pars[k - 1];
                    let p2 = <ST as HSurfaceTool>::d0(self.my_s, t1, t2);

                    if p1.distance_squared(p2) <= tol {
                        continue;
                    }

                    let d = (p2 - p1).normalize_or_zero();
                    let mut ok = true;
                    for l in jj + 1..k {
                        if a_v_flg[l - 1] {
                            ok = false;
                            break;
                        }

                        let pp = <ST as HSurfaceTool>::d0(self.my_s, t1, a_v_pars[l - 1]);
                        let dist = d.cross(pp - p1).length_squared();

                        if dist <= a_defl2 {
                            continue;
                        }

                        ok = false;
                        break;
                    }

                    if !ok {
                        jj = k - 1;
                        a_v_flg[jj - 1] = true;
                        broke = true;
                        break;
                    }

                    if a_v_flg[k - 1] {
                        jj = k;
                        broke = true;
                        break;
                    }
                }

                if !broke {
                    b_cont = false;
                }
            }
        }

        let mut my_nb_samples_v = a_v_flg.iter().filter(|&&f| f).count();

        if my_nb_samples_v < MY_MIN_PNTS {
            if my_nb_samples_v == 2 {
                // "uniform" distribution
                let nn = nbsv / MY_MIN_PNTS;
                a_v_flg[nn] = true;
                a_v_flg[nbsv - 1 - nn] = true;
                my_nb_samples_v = MY_MIN_PNTS;
            } else {
                // myNbSamplesU == 3: insert in bigger segment
                let mut i = 2usize;
                while !a_v_flg[i - 1] {
                    i += 1;
                }
                let j = if i < nbsv / 2 {
                    (i + (nbsv - i) / 2).min(nbsv - 1)
                } else {
                    (i / 2).max(2)
                };
                a_v_flg[j - 1] = true;
            }
            my_nb_samples_v = MY_MIN_PNTS;
        }

        // modified by NIZNHY-PKV Fri Dec 16 10:05:01 2011f
        // U
        let b_flag = (my_nb_samples_u as i32) < the_numin;
        if b_flag {
            my_nb_samples_u = nbsu;
        }
        let mut my_u_pars = Vec::with_capacity(my_nb_samples_u);
        if b_flag {
            my_u_pars = an_u_pars.clone();
        } else {
            for (i, &f) in an_u_flg.iter().enumerate() {
                if f {
                    my_u_pars.push(an_u_pars[i]);
                }
            }
        }
        // V
        let b_flag = (my_nb_samples_v as i32) < the_nvmin;
        if b_flag {
            my_nb_samples_v = nbsv;
        }
        let mut my_v_pars = Vec::with_capacity(my_nb_samples_v);
        if b_flag {
            my_v_pars = a_v_pars.clone();
        } else {
            for (i, &f) in a_v_flg.iter().enumerate() {
                if f {
                    my_v_pars.push(a_v_pars[i]);
                }
            }
        }
        // modified by NIZNHY-PKV Mon Dec 26 12:25:35 2011t

        self.my_nb_samples_u = my_nb_samples_u as i32;
        self.my_nb_samples_v = my_nb_samples_v as i32;
        self.my_u_pars = Some(my_u_pars);
        self.my_v_pars = Some(my_v_pars);
    }

    /// OCCT IsUniformSampling() (cxx L1649-1654).
    pub fn is_uniform_sampling(&self) -> bool {
        <ST as HSurfaceTool>::get_type(self.my_s) != SurfaceType::BSplineSurface
    }
}

/// OCCT GetConeApexParam(theC, theU, theV) (cxx L1660-1696) — computes the
/// cone's apex parameters.
pub fn get_cone_apex_param(the_c: &ConicalSurface, the_u: &mut f64, the_v: &mut f64) {
    // OCCT P = theC.Apex() — the true apex (R = 0 point on the axis), the
    // reference point shifted by -R/tan(SemiAngle) along the axis.
    let apex = the_c.apex - the_c.axis * (the_c.radius / the_c.half_angle_rad.tan());
    // OCCT T.SetTransformation(Pos) — the cone frame (X = ref_dir, Y =
    // axis ^ X, Z = axis), located at the reference point.
    let x_dir = the_c.ref_dir;
    let y_dir = the_c.axis.cross(x_dir);
    let rel = apex - the_c.apex;
    let ploc_x = rel.dot(x_dir);
    let ploc_y = rel.dot(y_dir);
    let ploc_z = rel.dot(the_c.axis);

    let radius = the_c.radius;
    let s_angle = the_c.half_angle_rad;

    if ploc_x == 0.0 && ploc_y == 0.0 {
        *the_u = 0.0;
    } else if -radius > ploc_z * s_angle.tan() {
        // The point is at the `wrong` side of the apex.
        *the_u = (-ploc_y).atan2(-ploc_x);
    } else {
        *the_u = ploc_y.atan2(ploc_x);
    }

    if *the_u < -1.0e-16 {
        *the_u += std::f64::consts::PI + std::f64::consts::PI;
    } else if *the_u < 0.0 {
        *the_u = 0.0;
    }

    *the_v = s_angle.sin() * (ploc_x * the_u.cos() + ploc_y * the_u.sin() - radius)
        + s_angle.cos() * ploc_z;
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{
        Circle3, ConicalSurface, CylindricalSurface, Ellipse3, Hyperbola3, Line3, Parabola3,
        Plane, SphericalSurface, ToroidalSurface,
    };
    use rcad_kernel::base::proj_lib::CurveType;

    /// The plane z = 0 over a stored window (S(u, v) = (u, v, 0)) — the
    /// OCCT adaptor model with UTrim/VTrim narrowing.
    #[derive(Clone)]
    struct TestPlane {
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    }

    struct TestPlaneTool;

    impl HSurfaceTool for TestPlaneTool {
        type Surface = TestPlane;
        type BasisCurve = TestBasisCurve;
        type BasisSurface = TestBasisSurface;

        fn first_u_parameter(s: &TestPlane) -> f64 {
            s.u0
        }
        fn first_v_parameter(s: &TestPlane) -> f64 {
            s.v0
        }
        fn last_u_parameter(s: &TestPlane) -> f64 {
            s.u1
        }
        fn last_v_parameter(s: &TestPlane) -> f64 {
            s.v1
        }
        fn nb_u_intervals(_s: &TestPlane, _sh: GeomAbsShape) -> usize {
            1
        }
        fn nb_v_intervals(_s: &TestPlane, _sh: GeomAbsShape) -> usize {
            1
        }
        fn u_intervals(s: &TestPlane, tab: &mut [f64], _sh: GeomAbsShape) {
            tab[0] = s.u0;
            tab[1] = s.u1;
        }
        fn v_intervals(s: &TestPlane, tab: &mut [f64], _sh: GeomAbsShape) {
            tab[0] = s.v0;
            tab[1] = s.v1;
        }
        fn is_u_closed(_s: &TestPlane) -> bool {
            false
        }
        fn is_v_closed(_s: &TestPlane) -> bool {
            false
        }
        fn is_u_periodic(_s: &TestPlane) -> bool {
            false
        }
        fn u_period(_s: &TestPlane) -> f64 {
            panic!("Standard_NoSuchObject");
        }
        fn is_v_periodic(_s: &TestPlane) -> bool {
            false
        }
        fn v_period(_s: &TestPlane) -> f64 {
            panic!("Standard_NoSuchObject");
        }
        fn value(s: &TestPlane, u: f64, v: f64) -> Point3 {
            s.point_at(u, v)
        }
        fn d0(s: &TestPlane, u: f64, v: f64) -> Point3 {
            s.point_at(u, v)
        }
        fn d1(s: &TestPlane, u: f64, v: f64) -> (Point3, rcad_kernel::geom::Vec3, rcad_kernel::geom::Vec3) {
            (s.point_at(u, v), DVec3::X, DVec3::Y)
        }
        fn d2(
            s: &TestPlane,
            u: f64,
            v: f64,
        ) -> (
            Point3,
            rcad_kernel::geom::Vec3,
            rcad_kernel::geom::Vec3,
            rcad_kernel::geom::Vec3,
            rcad_kernel::geom::Vec3,
            rcad_kernel::geom::Vec3,
        ) {
            let (p, du, dv) = (s.point_at(u, v), DVec3::X, DVec3::Y);
            (p, du, dv, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO)
        }
        fn dn(_s: &TestPlane, _u: f64, _v: f64, _nu: usize, _nv: usize) -> rcad_kernel::geom::Vec3 {
            panic!("Standard_NoSuchObject");
        }
        fn u_resolution(_s: &TestPlane, r3d: f64) -> f64 {
            r3d
        }
        fn v_resolution(_s: &TestPlane, r3d: f64) -> f64 {
            r3d
        }
        fn get_type(_s: &TestPlane) -> SurfaceType {
            SurfaceType::Plane
        }
        fn plane(_s: &TestPlane) -> Plane {
            panic!("Standard_NoSuchObject");
        }
        fn cylinder(_s: &TestPlane) -> CylindricalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn cone(_s: &TestPlane) -> ConicalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn torus(_s: &TestPlane) -> ToroidalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn sphere(_s: &TestPlane) -> SphericalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn axe_of_revolution(_s: &TestPlane) -> (Point3, rcad_kernel::geom::Vec3) {
            panic!("Standard_NoSuchObject");
        }
        fn direction(_s: &TestPlane) -> rcad_kernel::geom::Vec3 {
            panic!("Standard_NoSuchObject");
        }
        fn basis_curve(_s: &TestPlane) -> TestBasisCurve {
            panic!("Standard_NoSuchObject");
        }
        fn basis_surface(_s: &TestPlane) -> TestBasisSurface {
            panic!("Standard_NoSuchObject");
        }
        fn offset_value(_s: &TestPlane) -> f64 {
            panic!("Standard_NoSuchObject");
        }
        fn u_trim(s: &TestPlane, first: f64, last: f64, _tol: f64) -> TestPlane {
            let mut t = s.clone();
            t.u0 = first;
            t.u1 = last;
            t
        }
        fn v_trim(s: &TestPlane, first: f64, last: f64, _tol: f64) -> TestPlane {
            let mut t = s.clone();
            t.v0 = first;
            t.v1 = last;
            t
        }
        fn bezier(_s: &TestPlane) -> &rcad_kernel::geom::BezierSurface {
            panic!("Standard_NoSuchObject");
        }
        fn bspline(_s: &TestPlane) -> &rcad_kernel::geom::BSplineSurface {
            panic!("Standard_NoSuchObject");
        }
        fn nb_u_poles(_s: &TestPlane) -> usize {
            panic!("Standard_NoSuchObject");
        }
        fn nb_v_poles(_s: &TestPlane) -> usize {
            panic!("Standard_NoSuchObject");
        }
        fn nb_u_knots(_s: &TestPlane) -> usize {
            panic!("Standard_NoSuchObject");
        }
        fn nb_v_knots(_s: &TestPlane) -> usize {
            panic!("Standard_NoSuchObject");
        }
        fn u_degree(_s: &TestPlane) -> usize {
            panic!("Standard_NoSuchObject");
        }
        fn v_degree(_s: &TestPlane) -> usize {
            panic!("Standard_NoSuchObject");
        }
    }

    impl TestPlane {
        fn point_at(&self, u: f64, v: f64) -> Point3 {
            DVec3::new(u, v, 0.0)
        }
    }

    struct TestBasisCurve;
    impl crate::geomalgo::int_curve_surface::Adaptor3dCurveBasis for TestBasisCurve {
        fn get_type(&self) -> CurveType {
            panic!("Standard_NoSuchObject");
        }
        fn value(&self, _u: f64) -> Point3 {
            panic!("Standard_NoSuchObject");
        }
        fn line(&self) -> Line3 {
            panic!("Standard_NoSuchObject");
        }
        fn parabola(&self) -> Parabola3 {
            panic!("Standard_NoSuchObject");
        }
        fn hyperbola(&self) -> Hyperbola3 {
            panic!("Standard_NoSuchObject");
        }
    }

    struct TestBasisSurface;
    impl crate::geomalgo::int_curve_surface::Adaptor3dSurfaceBasis<TestBasisCurve> for TestBasisSurface {
        fn get_type(&self) -> SurfaceType {
            panic!("Standard_NoSuchObject");
        }
        fn plane(&self) -> Plane {
            panic!("Standard_NoSuchObject");
        }
        fn cylinder(&self) -> CylindricalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn cone(&self) -> ConicalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn direction(&self) -> rcad_kernel::geom::Vec3 {
            panic!("Standard_NoSuchObject");
        }
        fn basis_curve(&self) -> TestBasisCurve {
            panic!("Standard_NoSuchObject");
        }
    }

    /// OCCT anchor: the restriction iteration over a finite domain — four
    /// arcs in the OCCT order (Vinf along +X, Usup along +Y, Vsup along
    /// -X, Uinf along -Y) with their trim windows (cxx L57-209).
    #[test]
    fn topol_tool_restrictions_anchor() {
        let surface = TestPlane {
            u0: -10.0,
            u1: 10.0,
            v0: -5.0,
            v1: 5.0,
        };
        let mut tool = TopolTool::<TestPlane, TestPlaneTool>::new(&surface);

        assert_eq!(tool.more(), true);
        // Arc 1: the V = Vinf iso along +X, trimmed to the U range.
        {
            let r = tool.value();
            assert_eq!(r.line().origin, DVec2::new(0.0, -5.0));
            assert_eq!(r.line().direction, DVec2::X);
            assert_eq!(r.first_parameter(), -10.0);
            assert_eq!(r.last_parameter(), 10.0);
        }
        tool.next();
        // Arc 2: the U = Usup iso along +Y, trimmed to the V range.
        {
            let r = tool.value();
            assert_eq!(r.line().origin, DVec2::new(10.0, 0.0));
            assert_eq!(r.line().direction, DVec2::Y);
            assert_eq!(r.first_parameter(), -5.0);
            assert_eq!(r.last_parameter(), 5.0);
        }
        tool.next();
        // Arc 3: the V = Vsup iso along -X (OCCT gp_Dir2d::NX).
        {
            let r = tool.value();
            assert_eq!(r.line().origin, DVec2::new(0.0, 5.0));
            assert_eq!(r.line().direction, -DVec2::X);
        }
        tool.next();
        // Arc 4: the U = Uinf iso along -Y (OCCT gp_Dir2d::NY).
        {
            let r = tool.value();
            assert_eq!(r.line().origin, DVec2::new(-10.0, 0.0));
            assert_eq!(r.line().direction, -DVec2::Y);
        }
        tool.next();
        assert!(!tool.more());
    }

    /// OCCT anchor: Classify over the four-restriction domain (cxx L287-299).
    #[test]
    fn topol_tool_classify_anchor() {
        let surface = TestPlane {
            u0: -10.0,
            u1: 10.0,
            v0: -5.0,
            v1: 5.0,
        };
        let tool = TopolTool::<TestPlane, TestPlaneTool>::new(&surface);

        assert_eq!(tool.classify(DVec2::new(1.0, -2.0), 1e-9, true), State::In);
        assert_eq!(tool.classify(DVec2::new(20.0, 0.0), 1e-9, true), State::Out);
        assert_eq!(tool.classify(DVec2::new(10.0, 3.0), 1e-9, true), State::On);
        assert_eq!(tool.classify(DVec2::new(10.0 + 1e-4, 3.0), 1e-2, true), State::On);
        assert_eq!(tool.is_the_point_on(DVec2::new(-10.0, 0.0), 1e-9, true), true);
        assert_eq!(tool.is_the_point_on(DVec2::new(1.0, 1.0), 1e-9, true), false);
        assert!(!tool.domain_is_infinite());
        assert!(tool.has_3d() == false);
    }

    /// OCCT anchor: Initialize(C) end vertices (cxx L235-254) — FORWARD at
    /// the first parameter, REVERSED at the last.
    #[test]
    fn topol_tool_curve_vertices_anchor() {
        let surface = TestPlane {
            u0: -10.0,
            u1: 10.0,
            v0: -5.0,
            v1: 5.0,
        };
        let mut tool = TopolTool::<TestPlane, TestPlaneTool>::new(&surface);

        let arc = Line2dAdaptor::new_pnt_dir(DVec2::new(0.0, -5.0), DVec2::X, 2.0, 8.0);
        tool.initialize_curve(&arc);
        tool.init_vertex_iterator();

        assert!(tool.more_vertex());
        let v1 = HVertex::new_with(tool.vertex().value(), tool.vertex().orientation(), tool.vertex().resolution(&arc));
        assert_eq!(v1.value(), DVec2::new(2.0, -5.0));
        assert_eq!(v1.orientation(), Orientation::Forward);
        assert!((v1.parameter(&arc) - 2.0).abs() < 1e-12);
        tool.next_vertex();

        assert!(tool.more_vertex());
        let v2 = HVertex::new_with(tool.vertex().value(), tool.vertex().orientation(), tool.vertex().resolution(&arc));
        assert_eq!(v2.value(), DVec2::new(8.0, -5.0));
        assert_eq!(v2.orientation(), Orientation::Reversed);
        assert!(v1.is_same(&v1));
        assert!(!v1.is_same(&v2));
        tool.next_vertex();
        assert!(!tool.more_vertex());
    }

    /// OCCT anchor: SamplePnts on a plane — the Plane branch of
    /// ComputeSamplePoints gives 2, then the minimum-6 clamp (cxx L849-856)
    /// raises it to 6; SamplePnts fills uniform arrays with the domain
    /// endpoints at both ends (cxx L1019-1130).
    #[test]
    fn topol_tool_sample_pnts_plane_anchor() {
        let surface = TestPlane {
            u0: -10.0,
            u1: 10.0,
            v0: -5.0,
            v1: 5.0,
        };
        let mut tool = TopolTool::<TestPlane, TestPlaneTool>::new(&surface);
        assert_eq!(tool.nb_samples_u(), 6);
        assert_eq!(tool.nb_samples_v(), 6);
        assert_eq!(tool.nb_samples(), 36);
        assert!(tool.is_uniform_sampling());

        tool.sample_pnts(0.1, 10, 10);
        let mut upars = [0.0f64; 6];
        let mut vpars = [0.0f64; 6];
        tool.u_parameters(&mut upars);
        tool.v_parameters(&mut vpars);
        assert_eq!(upars[0], -10.0);
        assert_eq!(upars[5], 10.0);
        assert!((upars[1] - -6.0).abs() < 1e-12, "u1={}", upars[1]);
        assert!((upars[2] - -2.0).abs() < 1e-12, "u2={}", upars[2]);
        assert_eq!(vpars[0], -5.0);
        assert_eq!(vpars[5], 5.0);
        assert!((vpars[1] - -3.0).abs() < 1e-12, "v1={}", vpars[1]);
    }

    /// OCCT anchor: GetConeApexParam (cxx L1660-1696) — the apex of a cone
    /// (R = 1, alpha = 30 deg, axis +Z at the reference origin) lies on the
    /// axis so U = 0 and V = -R/sin(alpha) = -2.
    #[test]
    fn topol_tool_cone_apex_param_anchor() {
        let cone = ConicalSurface::new_with_ref_dir(
            DVec3::ZERO,
            DVec3::Z,
            1.0,
            std::f64::consts::PI / 6.0,
            DVec3::X,
        );
        let mut u = 0.0f64;
        let mut v = 0.0f64;
        get_cone_apex_param(&cone, &mut u, &mut v);
        assert!(u.abs() < 1e-12, "u={}", u);
        assert!((v - (-2.0)).abs() < 1e-12, "v={}", v);
    }
}
