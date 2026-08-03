// OCCT IntStart_SearchInside.gxx (IntPatch_TheSearchInside) 1:1 Rust
// translation — search for interior starting points (points where
// F(u,v) = 0 strictly inside the parametric surface domain).
//
// rcad adaptation: the OCCT TopolTool (T) is replaced by the surface itself:
// NbSamplesU/V/NbSamples come from a uniform grid, SamplePoint from the
// domain, and Classify by the rectangle membership test.

use rcad_kernel::geom::{Surface3, SurfaceEval};

use super::function_set_root::FunctionSetRoot;
use super::path_point::InteriorPoint;
use super::surf_function::SurfFunction;

/// OCCT IntPatch_TheSearchInside.
pub struct SearchInside {
    done: bool,
    list: Vec<InteriorPoint>,
}

impl SearchInside {
    /// OCCT default constructor (IntStart_SearchInside.gxx L26-29).
    pub fn new() -> Self {
        SearchInside {
            done: false,
            list: Vec::new(),
        }
    }

    /// OCCT Perform(F, PS, T, Epsilon) (L42-266).  `samples_u` / `samples_v`
    /// stand in for the TopolTool NbSamplesU/V.
    pub fn perform(
        &mut self,
        func: &mut SurfFunction,
        surf: &Surface3,
        epsilon: f64,
    ) {
        self.done = false;
        self.list.clear();

        let d = surf.default_domain();
        let umin0 = d[0];
        let umax0 = d[1];
        let vmin0 = d[2];
        let vmax0 = d[3];

        let nbsample_u = 10i32;
        let nbsample_v = 10i32;
        let nbsample = (nbsample_u * nbsample_v) as i32;

        let mut binf = [umin0, vmin0];
        let mut bsup = [umax0, vmax0];
        let umin = umin0;
        let umax = umax0;
        let vmin = vmin0;
        let vmax = vmax0;

        let mut du = bsup[0] - binf[0];
        let mut dv = bsup[1] - binf[1];
        du /= (nbsample_u as f64) * 0.5;
        dv /= (nbsample_v as f64) * 0.5;

        // Parametric resolutions.
        let resol_u = u_resolution(surf, rcad_kernel::precision::CONFUSION);
        let resol_v = v_resolution(surf, rcad_kernel::precision::CONFUSION);
        let toler = [resol_u, resol_v];
        let mut maxtoler1toler2 = if toler[0] > toler[1] { toler[0] } else { toler[1] };

        // On interdit aux points d'etre trop proches des restrictions.
        maxtoler1toler2 *= 1000.0;
        if maxtoler1toler2 > du * 0.001 {
            maxtoler1toler2 = du * 0.001;
        }
        if maxtoler1toler2 > dv * 0.001 {
            maxtoler1toler2 = dv * 0.001;
        }

        func.set_surface(surf.clone());
        let tol = func.tolerance();

        let mut rs_nld = FunctionSetRoot::new(func, toler);

        let mut umin_s = umin + du * 0.01;
        let mut vmin_s = vmin + dv * 0.01;
        let mut umax_s = umax - du * 0.01;
        let mut vmax_s = vmax - dv * 0.01;

        let mut i: i32 = 1;
        while i <= nbsample + 12 {
            let mut s2d = glam::DVec2::ZERO;
            let mut nepastester = false;
            let mut u_vap = [0.0f64; 2];

            if i <= nbsample {
                // T->SamplePoint(i, s2d, s3d) — uniform grid.
                let idx = (i - 1) as usize;
                let iu = idx % nbsample_u as usize;
                let iv = idx / nbsample_u as usize;
                let fu = if nbsample_u == 1 {
                    0.0
                } else {
                    iu as f64 / (nbsample_u - 1) as f64
                };
                let fv = if nbsample_v == 1 {
                    0.0
                } else {
                    iv as f64 / (nbsample_v - 1) as f64
                };
                s2d = glam::DVec2::new(umin + fu * (umax - umin), vmin + fv * (vmax - vmin));
                u_vap[0] = s2d.x;
                u_vap[1] = s2d.y;

                binf[0] = umin_s.max(u_vap[0] - du);
                binf[1] = vmin_s.max(u_vap[1] - dv);
                bsup[0] = umax_s.min(u_vap[0] + du);
                bsup[1] = vmax_s.min(u_vap[1] + dv);
                let u1 = binf[0];
                let v1 = binf[1];
                let u2 = bsup[0];
                let v2 = bsup[1];

                let p_extrm1 = surf.point_at(u1, v1);
                let p_extrm2 = surf.point_at(u2, v2);
                let Some(rvalf) = func.value(&u_vap) else {
                    i += 1;
                    continue;
                };
                let dist_pp = p_extrm1.distance_squared(p_extrm2);
                if rvalf * rvalf > 3.0 * dist_pp {
                    nepastester = true;
                }
            } else {
                // The 12 corner / edge points with small offsets.
                if i == nbsample + 1 {
                    s2d = glam::DVec2::new(umin_s + du * 0.02, vmin_s + dv * 0.02);
                } else if i == nbsample + 2 {
                    s2d = glam::DVec2::new(umax_s - du * 0.02, vmin_s + dv * 0.02);
                } else if i == nbsample + 3 {
                    s2d = glam::DVec2::new(umin_s + du * 0.02, vmax_s - dv * 0.02);
                } else if i == nbsample + 4 {
                    s2d = glam::DVec2::new(umax_s - du * 0.02, vmax_s - dv * 0.02);
                } else if i == nbsample + 5 {
                    s2d = glam::DVec2::new(umin_s + du * 0.02, vmin_s + dv * 0.02);
                } else if i == nbsample + 6 {
                    s2d = glam::DVec2::new(umax_s - du * 0.02, vmin_s + dv * 0.02);
                } else if i == nbsample + 7 {
                    s2d = glam::DVec2::new(umin_s + du * 0.02, vmax_s - dv * 0.02);
                } else if i == nbsample + 8 {
                    s2d = glam::DVec2::new(umax_s - du * 0.02, vmax_s - dv * 0.02);
                } else if i == nbsample + 9 {
                    s2d = glam::DVec2::new(umin_s + du * 0.005, vmin_s + dv * 0.005);
                } else if i == nbsample + 10 {
                    s2d = glam::DVec2::new(umax_s - du * 0.005, vmin_s + dv * 0.005);
                } else if i == nbsample + 11 {
                    s2d = glam::DVec2::new(umin_s + du * 0.005, vmax_s - dv * 0.005);
                } else {
                    s2d = glam::DVec2::new(umax_s - du * 0.005, vmax_s - dv * 0.005);
                }
                u_vap[0] = s2d.x;
                u_vap[1] = s2d.y;

                binf[0] = umin_s.max(u_vap[0] - du);
                binf[1] = vmin_s.max(u_vap[1] - dv);
                bsup[0] = umax_s.min(u_vap[0] + du);
                bsup[1] = vmax_s.min(u_vap[1] + dv);
            }

            if !nepastester {
                rs_nld.perform(func, u_vap, binf, bsup);
                if rs_nld.is_done() {
                    if func.root().abs() <= tol {
                        if !func.is_tangent() {
                            let psol = func.point();
                            let uvap = rs_nld.root();
                            // On regarde si le point trouve est bien un nouveau point.
                            let mut j = 0usize;
                            let nbpt = self.list.len();
                            let mut testpnt = j < nbpt;
                            while testpnt {
                                let ipj = &self.list[j];
                                let pj = ipj.value();
                                if (pj.x - psol.x).abs() <= epsilon
                                    && (pj.y - psol.y).abs() <= epsilon
                                    && (pj.z - psol.z).abs() <= epsilon
                                    && (uvap[0] - ipj.u_parameter()).abs() <= toler[0]
                                    && (uvap[1] - ipj.v_parameter()).abs() <= toler[1]
                                {
                                    testpnt = false;
                                } else {
                                    j += 1;
                                    testpnt = j < nbpt;
                                }
                            }
                            if j >= nbpt {
                                // situ = Classify(UV, Maxtoler1toler2, false) —
                                // inside the (margin-inset) rectangle.
                                if is_inside_domain(uvap[0], uvap[1], maxtoler1toler2, &d) {
                                    self.list.push(InteriorPoint::new_full(
                                        psol,
                                        uvap[0],
                                        uvap[1],
                                        func.direction_3d(),
                                        func.direction_2d(),
                                    ));
                                }
                            }
                        }
                    }
                }
            }
            i += 1;
        }
        self.done = true;
    }

    /// OCCT Perform(F, PS, UStart, VStart) (L270-313) — from a given start point.
    pub fn perform_from_point(
        &mut self,
        func: &mut SurfFunction,
        surf: &Surface3,
        u_start: f64,
        v_start: f64,
    ) {
        self.done = false;
        self.list.clear();

        let d = surf.default_domain();
        let binf = [d[0], d[2]];
        let bsup = [d[1], d[3]];

        let resol_u = u_resolution(surf, rcad_kernel::precision::CONFUSION);
        let resol_v = v_resolution(surf, rcad_kernel::precision::CONFUSION);
        let toler = [resol_u, resol_v];

        if u_start - binf[0] > -toler[0]
            && u_start - bsup[0] < toler[0]
            && v_start - binf[1] > -toler[1]
            && v_start - bsup[1] < toler[1]
        {
            func.set_surface(surf.clone());
            let u_vap = [u_start, v_start];

            let mut rs_nld = FunctionSetRoot::new(func, toler);
            rs_nld.perform(func, u_vap, binf, bsup);
            if rs_nld.is_done() {
                let tol = func.tolerance();
                let valf = func.root();
                if valf.abs() <= tol && !func.is_tangent() {
                    let psol = func.point();
                    let uvap = rs_nld.root();
                    let intp = InteriorPoint::new_full(
                        psol,
                        uvap[0],
                        uvap[1],
                        func.direction_3d(),
                        func.direction_2d(),
                    );
                    self.list.push(intp);
                }
            }
        }
        self.done = true;
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.list.len()
    }
    /// OCCT Value(Index) — 1-based.
    pub fn value(&self, index: usize) -> &InteriorPoint {
        &self.list[index - 1]
    }
}

impl Default for SearchInside {
    fn default() -> Self {
        Self::new()
    }
}

/// rcad adaptation of Adaptor3d_HSurfaceTool::UResolution — parametric
/// resolution in the U direction for a 3D tolerance (default: the domain
/// extent scaled by the tolerance ratio).
fn u_resolution(surf: &Surface3, tol3d: f64) -> f64 {
    let d = surf.default_domain();
    let u_extent = (d[1] - d[0]).abs();
    if u_extent > 1e-12 {
        tol3d.max(1e-9) / u_extent
    } else {
        rcad_kernel::precision::PCONFUSION
    }
}

fn v_resolution(surf: &Surface3, tol3d: f64) -> f64 {
    let d = surf.default_domain();
    let v_extent = (d[3] - d[2]).abs();
    if v_extent > 1e-12 {
        tol3d.max(1e-9) / v_extent
    } else {
        rcad_kernel::precision::PCONFUSION
    }
}

/// rcad adaptation of IntPatch_HInterTool::Classify — the point is IN when it
/// lies inside the domain, inset by the classification tolerance.
fn is_inside_domain(u: f64, v: f64, margin: f64, d: &[f64; 4]) -> bool {
    u >= d[0] + margin && u <= d[1] - margin && v >= d[2] + margin && v <= d[3] - margin
}
