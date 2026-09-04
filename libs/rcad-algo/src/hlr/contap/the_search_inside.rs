// OCCT Contap_TheSearchInside = IntStart_SearchInside.gxx instantiated for
// Contap (ThePSurface = handle<Adaptor3d_Surface>,
// ThePSurfaceTool = Adaptor3d_HSurfaceTool -> the SurfaceAdapter calls,
// TheTopolTool = Adaptor3d_TopolTool, TheSITool = Contap_HContTool,
// TheFunction = Contap_SurfFunction).
//
// IntStart_SearchInside.gxx L26-313.  The solver is math_FunctionSetRoot
// (rcad imp_prm::FunctionSetRoot, the 2-variable specialization).

use rcad_kernel::geom::Point3;
use rcad_kernel::precision::CONFUSION;

use crate::geomalgo::int_patch::imp_prm::function_set_root::FunctionSetRoot;
use crate::geomalgo::int_patch::imp_prm::path_point::InteriorPoint;

use crate::hlr::contap::domain::ContapDomain;
use crate::hlr::contap::surf_function::SurfFunction;
use crate::hlr::contap::surface_adaptor::SurfaceHandle;

/// OCCT IntStart_SearchInside.
#[derive(Debug, Clone)]
pub struct TheSearchInside {
    done: bool,
    list: Vec<InteriorPoint>,
}

impl TheSearchInside {
    /// OCCT IntStart_SearchInside() (gxx L26-29).
    pub fn new() -> Self {
        TheSearchInside {
            done: false,
            list: Vec::new(),
        }
    }

    /// OCCT IntStart_SearchInside(Func, PS, T, Epsilon) (gxx L31-38).
    pub fn with_perform(
        func: &mut SurfFunction,
        surf: &SurfaceHandle,
        t: &mut dyn ContapDomain,
        epsilon: f64,
    ) -> Self {
        let mut s = TheSearchInside::new();
        s.perform(func, surf, t, epsilon);
        s
    }

    /// OCCT Perform(Func, PS, T, Epsilon) (gxx L42-266).
    pub fn perform(
        &mut self,
        func: &mut SurfFunction,
        ps: &SurfaceHandle,
        t: &mut dyn ContapDomain,
        epsilon: f64,
    ) {
        self.done = false;
        self.list.clear();

        let umin = ps.first_u_parameter();
        let vmin = ps.first_v_parameter();
        let umax = ps.last_u_parameter();
        let vmax = ps.last_v_parameter();
        let mut binf = [umin, vmin];
        let mut bsup = [umax, vmax];

        let nbsample_u = t.nb_samples_u();
        let nbsample_v = t.nb_samples_v();
        let nbsample = t.nb_samples();

        let mut du = bsup[0] - binf[0];
        let mut dv = bsup[1] - binf[1];
        du /= nbsample_u as f64 * 0.5;
        dv /= nbsample_v as f64 * 0.5;

        let toler1 = ps.u_resolution(CONFUSION);
        let toler2 = ps.v_resolution(CONFUSION);
        let mut maxtoler1toler2 = toler1;
        if toler2 > maxtoler1toler2 {
            maxtoler1toler2 = toler2;
        }

        //-- lbr le 15 mai 97
        //-- on interdit aux points d'etre trop prets des restrictions
        maxtoler1toler2 *= 1000.0;
        if maxtoler1toler2 > du * 0.001 {
            maxtoler1toler2 = du * 0.001;
        }
        if maxtoler1toler2 > dv * 0.001 {
            maxtoler1toler2 = dv * 0.001;
        }

        func.set(ps.clone());
        let tol = func.tolerance();

        let toler = [toler1, toler2];
        let mut rs_nld = FunctionSetRoot::new(func, toler);

        let mut umin = umin + du * 0.01;
        let mut vmin = vmin + dv * 0.01;
        let mut umax = umax - du * 0.01;
        let mut vmax = vmax - dv * 0.01;

        //-- lbr le 30 octobre 97 (coin comments) --
        let mut i = 1;
        while i <= nbsample + 12 {
            let mut s2d = glam::DVec2::ZERO;
            let mut nepastester = false;
            let mut uvap = [0.0f64; 2];
            if i <= nbsample {
                let (s2d_p, _s3d) = t.sample_point(i as usize);
                s2d = s2d_p;
                uvap[0] = s2d.x;
                uvap[1] = s2d.y;

                binf[0] = umin.max(uvap[0] - du);
                binf[1] = vmin.max(uvap[1] - dv);
                bsup[0] = umax.min(uvap[0] + du);
                bsup[1] = vmax.min(uvap[1] + dv);
                let u1 = binf[0];
                let v1 = binf[1];
                let u2 = bsup[0];
                let v2 = bsup[1];

                //-- gp_Pnt Pmilieu = ThePSurfaceTool::Value(PS,0.5*(u1+u2),0.5*(v1+v2));
                let pextrm1 = ps.value(u1, v1);
                let pextrm2 = ps.value(u2, v2);
                let Some(rvalf) = func.value(&uvap) else {
                    i += 1;
                    continue;
                };
                let dist_pp = pextrm1.distance_squared(pextrm2);
                if rvalf * rvalf > 3.0 * dist_pp {
                    nepastester = true;
                }
            } else {
                if i == nbsample + 1 {
                    s2d = glam::DVec2::new(umin + du * 0.02, vmin + dv * 0.02);
                } else if i == nbsample + 2 {
                    s2d = glam::DVec2::new(umax - du * 0.02, vmin + dv * 0.02);
                } else if i == nbsample + 3 {
                    s2d = glam::DVec2::new(umin + du * 0.02, vmax - dv * 0.02);
                } else if i == nbsample + 4 {
                    s2d = glam::DVec2::new(umax - du * 0.02, vmax - dv * 0.02);
                } else if i == nbsample + 5 {
                    s2d = glam::DVec2::new(umin + du * 0.02, vmin + dv * 0.02);
                } else if i == nbsample + 6 {
                    s2d = glam::DVec2::new(umax - du * 0.02, vmin + dv * 0.02);
                } else if i == nbsample + 7 {
                    s2d = glam::DVec2::new(umin + du * 0.02, vmax - dv * 0.02);
                } else if i == nbsample + 8 {
                    s2d = glam::DVec2::new(umax - du * 0.02, vmax - dv * 0.02);
                } else if i == nbsample + 9 {
                    s2d = glam::DVec2::new(umin + du * 0.005, vmin + dv * 0.005);
                } else if i == nbsample + 10 {
                    s2d = glam::DVec2::new(umax - du * 0.005, vmin + dv * 0.005);
                } else if i == nbsample + 11 {
                    s2d = glam::DVec2::new(umin + du * 0.005, vmax - dv * 0.005);
                } else {
                    s2d = glam::DVec2::new(umax - du * 0.005, vmax - dv * 0.005);
                }

                uvap[0] = s2d.x;
                uvap[1] = s2d.y;

                binf[0] = umin.max(uvap[0] - du);
                binf[1] = vmin.max(uvap[1] - dv);
                bsup[0] = umax.min(uvap[0] + du);
                bsup[1] = vmax.min(uvap[1] + dv);
            }

            if !nepastester {
                rs_nld.perform(func, uvap, binf, bsup);
                if rs_nld.is_done() {
                    if func.root().abs() <= tol {
                        if !func.is_tangent() {
                            let psol = func.point();
                            uvap = rs_nld.root();
                            // On regarde si le point trouve est bien un
                            // nouveau point.
                            let mut j = 1usize;
                            let nbpt = self.list.len();
                            let mut testpnt = j <= nbpt;

                            while testpnt {
                                let ipj = &self.list[j - 1];
                                let pj = ipj.value();
                                if (pj.x - psol.x).abs() <= epsilon
                                    && (pj.y - psol.y).abs() <= epsilon
                                    && (pj.z - psol.z).abs() <= epsilon
                                    && (uvap[0] - ipj.u_parameter()).abs() <= toler1
                                    && (uvap[1] - ipj.v_parameter()).abs() <= toler2
                                {
                                    testpnt = false;
                                } else {
                                    j += 1;
                                    testpnt = j <= nbpt;
                                }
                            }
                            if j > nbpt {
                                //    situ = TheSITool::Classify(PS,UVap(1),UVap(2));
                                let situ = t.classify(
                                    glam::DVec2::new(uvap[0], uvap[1]),
                                    maxtoler1toler2,
                                    false,
                                ); //-- ,false pour ne pas recadrer on Periodic
                                if situ == rcad_kernel::topo::topods::State::In {
                                    let d2d = func.direction_2d();
                                    self.list.push(InteriorPoint::new_full(
                                        psol,
                                        uvap[0],
                                        uvap[1],
                                        func.direction_3d(),
                                        d2d,
                                    ));
                                }
                            }
                        }
                    }
                }
            }
            i += 1;
        }
        done_marker(&mut self.done);
    }

    /// OCCT Perform(Func, PS, UStart, VStart) (gxx L270-313).
    pub fn perform_start(
        &mut self,
        func: &mut SurfFunction,
        ps: &SurfaceHandle,
        ustart: f64,
        vstart: f64,
    ) {
        self.done = false;
        self.list.clear();

        let binf = [ps.first_u_parameter(), ps.first_v_parameter()];
        let bsup = [ps.last_u_parameter(), ps.last_v_parameter()];

        let toler = [
            ps.u_resolution(CONFUSION),
            ps.v_resolution(CONFUSION),
        ];

        if ustart - binf[0] > -toler[0]
            && ustart - bsup[0] < toler[0]
            && vstart - binf[1] > -toler[1]
            && vstart - bsup[1] < toler[1]
        {
            func.set(ps.clone());
            let mut uvap = [ustart, vstart];

            let mut rs_nld = FunctionSetRoot::new(func, toler);
            rs_nld.perform(func, uvap, binf, bsup);
            if rs_nld.is_done() {
                let tol = func.tolerance();
                let valf = func.root();
                if valf.abs() <= tol && !func.is_tangent() {
                    let psol = func.point();
                    uvap = rs_nld.root();
                    let d2d = func.direction_2d();
                    let intp = InteriorPoint::new_full(
                        psol,
                        uvap[0],
                        uvap[1],
                        func.direction_3d(),
                        d2d,
                    );
                    self.list.push(intp);
                }
            }
        }

        self.done = true;
    }

    /// OCCT IsDone (IntStart_SearchInside.lxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT NbPoints.
    pub fn nb_points(&self) -> usize {
        self.list.len()
    }

    /// OCCT Value(Index) — 1-based.
    pub fn value(&self, index: usize) -> &InteriorPoint {
        &self.list[index - 1]
    }
}

impl Default for TheSearchInside {
    fn default() -> Self {
        Self::new()
    }
}

// gxx L265: `done = true;`
fn done_marker(done: &mut bool) {
    *done = true;
}

// The Point3 shape reference for the s3d out-value of SamplePoint.
#[allow(unused)]
fn _pnt3_shape() -> Point3 {
    Point3::ZERO
}
