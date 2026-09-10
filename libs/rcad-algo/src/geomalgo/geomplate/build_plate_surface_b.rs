//! Continuation of build_plate_surface.rs — the curve-path member bodies of
//! OCCT GeomPlate_BuildPlateSurface (the #[path]-declared child module `b`;
//! the split keeps both files under the 2000-line limit).
//!
//! Anchors (GeomPlate_BuildPlateSurface.cxx):
//! ProjectCurve (L254-303), ProjectedCurve (L307-349), CourbeJointive
//! (L1349-1445), Intersect (L1913-2145), Discretise (L2158-2358),
//! CalculNbPtsInit (L2366-2400), LoadCurve (L2407-2516), VerifSurface
//! (L2588-2732), EcartContraintesMil (L751-866), Disc2dContour (L871-1019),
//! Disc3dContour (L1024-1163), G0/G1/G2Error(Index) (L1260-1322).

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor3dCurve, Curve2dHandle, CurveHandle, CurveOnSurface, SurfaceHandle,
};
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomSurfaceAdaptor;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::core::precision::p_confusion;
use rcad_kernel::geom::{Curve2d, Surface3, SurfaceEval};
use rcad_kernel::Curve2dEval;
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::gcpnts_abscissa_point::gcpnts_length_2d_range;
use crate::geomalgo::gcpnts_curve::{GCPntsCurve, GCPntsCurve2d};
use crate::geomalgo::geom2d_int::{Curve2dAdaptor, GInter};
use crate::geomalgo::law::law_function::LawFunction;
use crate::geomalgo::law::law_interpol::LawInterpol;
use crate::geomalgo::plate::{FreeGtoCConstraint, GtoCConstraint, PinpointConstraint, PlateD1, PlateD2};

use super::{vec2d_angle, BuildPlateSurface};

impl BuildPlateSurface {
    /// OCCT ProjectCurve (L254-303) — projects a curve on the initial
    /// surface (a plane in the OCCT usage).
    pub(super) fn project_curve(&self, curv: &CurveHandle) -> Option<Curve2d> {
        // hsur = new GeomAdaptor_Surface(mySurfInit).
        let hsur: SurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(
            self.my_surf_init.as_ref().unwrap().clone(),
        ));

        // occ::handle<ProjLib_HCompProjectedCurve> HProjector =
        //     new ProjLib_HCompProjectedCurve(hsur, Curv, myTol3d / 10, myTol3d / 10);
        //
        // GAP leaf: ProjLib_CompProjectedCurve lives in the
        // geomalgo/proj_lib_h_comp_projected_curve pair, which is not
        // registered in the build yet (mid-integration); the anchor
        // preserves the dependency failure path.  The rest of the OCCT body
        // (single-pnt Bezier fallback, the Approx_CurveOnSurface call with
        // MaxSeg = 20 + NbIntervals(GeomAbs_C3), MaxDegree = 10,
        // Continuity = GeomAbs_C1) lands together with it.
        let _ = (&hsur, curv);
        unimplemented!(
            "ProjLib_HCompProjectedCurve (GeomPlate ProjectCurve) is not registered in the build yet"
        );
    }

    /// OCCT ProjectedCurve (L307-349) — the projection of a curve on the
    /// initial surface returned as the projected-curve adaptor.
    pub fn projected_curve(&self, curv: &CurveHandle) -> Option<Curve2dHandle> {
        // hsur = new GeomAdaptor_Surface(mySurfInit).
        let hsur: SurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(
            self.my_surf_init.as_ref().unwrap().clone(),
        ));

        // occ::handle<ProjLib_HCompProjectedCurve> HProjector =
        //     new ProjLib_HCompProjectedCurve(hsur, Curv, myTolU / 10, myTolV / 10);
        //
        // GAP leaf: same ProjLib_CompProjectedCurve dependency as
        // ProjectCurve (unregistered mid-integration); the OCCT
        // NbCurves/Bounds/Trim flow lands together with it.
        let _ = (&hsur, curv);
        unimplemented!(
            "ProjLib_HCompProjectedCurve (GeomPlate ProjectedCurve) is not registered in the build yet"
        );
    }

    /// OCCT CourbeJointive (L1349-1445) — creates a chain of curves to
    /// calculate the initial surface with the method of max flow; returns
    /// true if it is a closed contour.
    pub(super) fn courbe_jointive(&mut self, tolerance: f64) -> bool {
        let nbf = self.my_lin_cont.len();
        // mySense = new NCollection_HArray1<int>(1, nbf, 0).
        self.my_sense = Some(vec![0i32; nbf]);
        let mut result = true;
        let mut j = 1usize;

        while j <= (self.my_nb_bounds as usize - 1) {
            let mut a = 0;
            let mut i = j + 1;
            if i > self.my_nb_bounds as usize {
                result = false;
                a = 2;
            }
            while a < 1 {
                if i > self.my_nb_bounds as usize {
                    result = false;
                    a = 2;
                } else {
                    let mut uinit1 = self.my_lin_cont[j - 1].first_parameter();
                    let mut ufinal1 = self.my_lin_cont[j - 1].last_parameter();
                    let uinit2 = self.my_lin_cont[i - 1].first_parameter();
                    let ufinal2 = self.my_lin_cont[i - 1].last_parameter();
                    if self.my_sense.as_ref().unwrap()[j - 1] == 1 {
                        ufinal1 = uinit1;
                    }
                    let p1 = self.my_lin_cont[j - 1].d0(ufinal1);
                    let mut p2 = self.my_lin_cont[i - 1].d0(uinit2);
                    if p1.distance(p2) < tolerance {
                        if i != j + 1 {
                            // Swap the constraints j+1 <-> i and the same
                            // slots of myInitOrder (see TrierTab).
                            let tampon = self.my_lin_cont[j].clone();
                            self.my_lin_cont[j] = self.my_lin_cont[i - 1].clone();
                            self.my_lin_cont[i - 1] = tampon;
                            let tmp = self.my_init_order.as_ref().unwrap()[j];
                            self.my_init_order.as_mut().unwrap()[j] =
                                self.my_init_order.as_ref().unwrap()[i - 1];
                            self.my_init_order.as_mut().unwrap()[i - 1] = tmp;
                        }
                        a = 2;
                        self.my_sense.as_mut().unwrap()[j] = 0;
                    } else {
                        p2 = self.my_lin_cont[i - 1].d0(ufinal2);

                        if p1.distance(p2) < tolerance {
                            if i != j + 1 {
                                let tampon = self.my_lin_cont[j].clone();
                                self.my_lin_cont[j] = self.my_lin_cont[i - 1].clone();
                                self.my_lin_cont[i - 1] = tampon;
                                let tmp = self.my_init_order.as_ref().unwrap()[j];
                                self.my_init_order.as_mut().unwrap()[j] =
                                    self.my_init_order.as_ref().unwrap()[i - 1];
                                self.my_init_order.as_mut().unwrap()[i - 1] = tmp;
                            }
                            a = 2;
                            self.my_sense.as_mut().unwrap()[j] = 1;
                        }
                    }
                }
                i += 1;
            }
            j += 1;
        }
        let mut uinit1 = self.my_lin_cont[self.my_nb_bounds as usize - 1].first_parameter();
        let ufinal1 = self.my_lin_cont[self.my_nb_bounds as usize - 1].last_parameter();
        let uinit2 = self.my_lin_cont[0].first_parameter();
        let _ufinal2 = self.my_lin_cont[0].last_parameter();
        let p1 = self.my_lin_cont[self.my_nb_bounds as usize - 1].d0(ufinal1);
        let p2 = self.my_lin_cont[0].d0(uinit2);
        let _ = &mut uinit1;
        if (self.my_sense.as_ref().unwrap()[self.my_nb_bounds as usize - 1] == 0)
            && (p1.distance(p2) < tolerance)
        {
            return result;
        }
        let p1 = self.my_lin_cont[self.my_nb_bounds as usize - 1].d0(uinit1);
        if (self.my_sense.as_ref().unwrap()[self.my_nb_bounds as usize - 1] == 1)
            && (p1.distance(p2) < tolerance)
        {
            return result;
        }
        false
    }

    /// OCCT Intersect (L1913-2145) — finds the intersections between the 2d
    /// curves; a compatible intersection (G1-G1 / G0-G1) removes the zone on
    /// one of the two curves.
    pub(super) fn intersect(
        &mut self,
        pnt_inter: &mut [Vec<f64>],
        pnt_g1g1: &mut [Vec<f64>],
    ) {
        let ntlincont = self.my_lin_cont.len();
        // Geom2dInt_GInter Intersection; Geom2dAdaptor_Curve Ci, Cj.
        let mut intersection: GInter = GInter::new();
        let mut ci: Option<Curve2d> = None;
        let mut cj: Option<Curve2d> = None;

        for i in 1..=ntlincont {
            // Find the intersection with each curve including the curve itself.
            ci = Some(self.my_lin_cont[i - 1]
                .curve2d_on_surf()
                .expect("Intersect: Curve2dOnSurf is null"));
            for j in i..=ntlincont {
                cj = Some(self.my_lin_cont[j - 1]
                    .curve2d_on_surf()
                    .expect("Intersect: Curve2dOnSurf is null"));
                let ci_ref = ci.as_ref().unwrap();
                let cj_ref = cj.as_ref().unwrap();
                if i == j {
                    // Intersection.Perform(Ci, myTol2d * 10, myTol2d * 10).
                    intersection.perform_c(ci_ref, self.my_tol2d * 10.0, self.my_tol2d * 10.0);
                } else {
                    // Intersection.Perform(Ci, Cj, myTol2d * 10, myTol2d * 10).
                    intersection.perform_cc(
                        ci_ref,
                        cj_ref,
                        self.my_tol2d * 10.0,
                        self.my_tol2d * 10.0,
                    );
                }

                // if (!Intersection.IsEmpty()) — there is one intersection.
                if !intersection.base.is_empty() {
                    let nbpt = intersection.nb_points();
                    // number of points of intersection
                    for k in 1..=nbpt {
                        let int2d = intersection.point(k);
                        let param_on_first = int2d.param_on_first();
                        let param_on_second = int2d.param_on_second();
                        let p1 = self.my_lin_cont[i - 1].d0(param_on_first);
                        let p2 = self.my_lin_cont[j - 1].d0(param_on_second);
                        if p1.distance(p2) < self.my_tol3d {
                            // 2D intersection corresponds to close 3D points.
                            // The point on curve i is removed; the point on
                            // curve j is preserved; the length of interval is
                            // a 2d length corresponding in 3d to myTol3d.
                            let mut tolint = ci_ref.resolution(self.my_tol3d);
                            let (_p2d, v2d) = ci_ref.d1(param_on_first);
                            let mut aux = v2d.length();
                            if aux > 1.0e-7 {
                                aux = self.my_tol3d / aux;
                                if aux > 100.0 * tolint {
                                    tolint *= 100.0;
                                } else {
                                    tolint = aux;
                                }
                            } else {
                                tolint *= 100.0;
                            }

                            pnt_inter[i - 1].push(param_on_first - tolint);
                            pnt_inter[i - 1].push(param_on_first + tolint);
                            // If G1-G1
                            if (self.my_lin_cont[i - 1].order() == 1)
                                && (self.my_lin_cont[j - 1].order() == 1)
                            {
                                let (_pa, v11, v12, _v13, _v14, v15) =
                                    self.my_lin_cont[i - 1].d2(param_on_first);
                                let (_pb, v21, v22, _v23, _v24, v25) =
                                    self.my_lin_cont[j - 1].d2(param_on_second);
                                let v16 = v11.cross(v12);
                                let v26 = v21.cross(v22);
                                let mut ant = vec3_angle(v16, v26);
                                if ant > (std::f64::consts::PI / 2.0) {
                                    ant = std::f64::consts::PI - ant;
                                }
                                if ((v16.dot(v15) - v16.dot(v25)).abs() > (self.my_tol3d / 1000.0))
                                    || (ant.abs() > self.my_tol3d / 1000.0)
                                {
                                    // Non-compatible ==> remove zone in
                                    // constraint G1 corresponding to 3D
                                    // tolerance of 0.01.
                                    let mut tol = 100.0 * self.my_tol3d;
                                    let (_p1temp, v1) = ci_ref.d1(param_on_first);
                                    let (_p2temp, v2) = cj_ref.d1(param_on_second);
                                    let mut a1 = vec2d_angle(v1, v2);
                                    if a1 > (std::f64::consts::PI / 2.0) {
                                        a1 = std::f64::consts::PI - a1;
                                    }
                                    if ((a1.abs() - std::f64::consts::PI).abs()) < self.my_tolang
                                    {
                                        tol = 100000.0 * self.my_tol3d;
                                    }

                                    let coin = ci_ref.resolution(tol);
                                    let par1 = param_on_first - coin;
                                    let par2 = param_on_first + coin;
                                    // Storage of the interval for curve i.
                                    pnt_g1g1[i - 1].push(par1);
                                    pnt_g1g1[i - 1].push(par2);
                                    let coin = cj_ref.resolution(tol);
                                    let par1 = param_on_second - coin;
                                    let par2 = param_on_second + coin;
                                    // Storage of the interval for curve j.
                                    pnt_g1g1[j - 1].push(par1);
                                    pnt_g1g1[j - 1].push(par2);
                                }
                            }
                            // If G0-G1
                            if (self.my_lin_cont[i - 1].order() == 0
                                && self.my_lin_cont[j - 1].order() == 1)
                                || (self.my_lin_cont[i - 1].order() == 1
                                    && self.my_lin_cont[j - 1].order() == 0)
                            {
                                let vec: DVec3;
                                let (vec_u, vec_v);
                                if self.my_lin_cont[i - 1].order() == 0 {
                                    let the_curve = self.my_lin_cont[i - 1].curve3d();
                                    let (p1d, v) = the_curve.unwrap().d1(param_on_first);
                                    vec = v;
                                    let _ = p1d;
                                    let (_p2d, u, w) = self.my_lin_cont[j - 1].d1(param_on_second);
                                    vec_u = u;
                                    vec_v = w;
                                } else {
                                    let the_curve = self.my_lin_cont[j - 1].curve3d();
                                    let (p2d, v) = the_curve.unwrap().d1(param_on_second);
                                    vec = v;
                                    let _ = p2d;
                                    let (_p1d, u, w) = self.my_lin_cont[i - 1].d1(param_on_first);
                                    vec_u = u;
                                    vec_v = w;
                                }
                                let n = vec_u.cross(vec_v);
                                let mut angle = vec3_angle(vec, n);
                                angle = (std::f64::consts::PI / 2.0 - angle).abs();
                                if angle > self.my_tolang / 10.0 {
                                    // Non-compatible ==> one removes zone in
                                    // constraint G0 and G1 corresponding to
                                    // 3D tolerance of 0.01.
                                    let mut tol = 100.0 * self.my_tol3d;
                                    let (_p1temp, v1) = ci_ref.d1(param_on_first);
                                    let (_p2temp, v2) = cj_ref.d1(param_on_second);
                                    let mut a1 = vec2d_angle(v1, v2);
                                    if a1 > std::f64::consts::PI / 2.0 {
                                        a1 = std::f64::consts::PI - a1;
                                    }
                                    if ((a1.abs() - std::f64::consts::PI).abs()) < self.my_tolang
                                    {
                                        tol = 100000.0 * self.my_tol3d;
                                    }
                                    if self.my_lin_cont[i - 1].order() == 1 {
                                        let mut coin = ci_ref.resolution(tol);
                                        coin *= angle / self.my_tolang * 10.0;
                                        let par1 = param_on_first - coin;
                                        let par2 = param_on_first + coin;
                                        // Storage of the interval for curve i.
                                        pnt_g1g1[i - 1].push(par1);
                                        pnt_g1g1[i - 1].push(par2);
                                    } else {
                                        let mut coin = cj_ref.resolution(tol);
                                        coin *= angle / self.my_tolang * 10.0;
                                        let par1 = param_on_second - coin;
                                        let par2 = param_on_second + coin;
                                        // Storage of the interval for curve j.
                                        pnt_g1g1[j - 1].push(par1);
                                        pnt_g1g1[j - 1].push(par2);
                                    }
                                }
                            }
                        } else {
                            // 2D intersection corresponds to extended 3D
                            // points.  The point on curve i is removed, the
                            // point on curve j is preserved.
                            let dist = p1.distance(p2);
                            let tolint = ci_ref.resolution(dist);
                            pnt_inter[i - 1].push(param_on_first - tolint);
                            pnt_inter[i - 1].push(param_on_first + tolint);
                            if j != i {
                                let tolint = cj_ref.resolution(dist);
                                pnt_inter[j - 1].push(param_on_second - tolint);
                                pnt_inter[j - 1].push(param_on_second + tolint);
                            }
                        }
                    }
                }
            }
        }
    }

    /// OCCT Discretise (L2158-2358) — discretizes the curves according to
    /// the parameters; myPlateCont excludes duplicate points and incompatible
    /// zones.
    pub(super) fn discretise(&mut self, pnt_inter: &[Vec<f64>], pnt_g1g1: &[Vec<f64>]) {
        let ntlincont = self.my_lin_cont.len();
        // Law_Interpol acrlaw = new Law_Interpol().
        let mut acrlaw = LawInterpol::new();
        self.my_plate_cont = Some(vec![Vec::new(); ntlincont]);
        self.my_par_cont = Some(vec![Vec::new(); ntlincont]);

        // Construction of the table containing parameters of constraint points.
        let mut length2d = 0.0f64;

        for i in 1..=ntlincont {
            let lin_cont = self.my_lin_cont[i - 1].clone();
            let uinit = lin_cont.first_parameter();
            let ufinal = lin_cont.last_parameter();
            // C2d = LinCont->Curve2dOnSurf(); ACR = (!C2d.IsNull()).
            let c2d = lin_cont.curve2d_on_surf();
            let acr = c2d.is_some();
            if acr {
                let c2d = c2d.as_ref().unwrap();
                // Construct a law close to curvilinear abscissa.
                let nbint = 20usize;
                let mut tab_p2d: Vec<DVec2> = vec![DVec2::ZERO; nbint + 1];
                // tabP2d(1).SetY(Uinit); tabP2d(1).SetX(0.).
                tab_p2d[0].y = uinit;
                tab_p2d[0].x = 0.0;
                tab_p2d[nbint].y = ufinal;
                // Length2d = GCPnts_AbscissaPoint::Length(AC2d, Uinit, Ufinal).
                let ac2d = GCPntsCurve2d::new(c2d);
                length2d = gcpnts_length_2d_range(&ac2d, uinit, ufinal);
                tab_p2d[nbint].x = length2d;
                for ii in 2..nbint {
                    let u = uinit + (ufinal - uinit)
                        * ((1.0 - ((ii as f64 - 1.0) * std::f64::consts::PI / (nbint as f64)).cos())
                            / 2.0);
                    tab_p2d[ii - 1].y = u;
                    // tabP2d(ii).SetX(GCPnts_AbscissaPoint::Length(AC2d, Uinit, U)).
                    tab_p2d[ii - 1].x = gcpnts_length_2d_range(&ac2d, uinit, u);
                }
                acrlaw.set(&tab_p2d, false);
            }

            let nb_pnt_i = lin_cont.nb_points();
            let nb_pt_inter = pnt_inter[i - 1].len();
            let nb_pt_g1g1 = pnt_g1g1[i - 1].len();

            for j in 1..=nb_pnt_i {
                // Distribution of points in cosine following ACR 2D.
                let inter;
                if j == nb_pnt_i {
                    inter = ufinal; // to avoid bug on Sun
                } else if acr {
                    let c2d = c2d.as_ref().unwrap();
                    let cur_length = length2d
                        * (1.0 - (((j as f64 - 1.0) * std::f64::consts::PI
                            / (nb_pnt_i as f64 - 1.0))
                            .cos()))
                            / 2.0;
                    inter = acrlaw.value(cur_length);
                    let _ = c2d;
                } else {
                    inter = uinit
                        + (ufinal - uinit)
                            * ((1.0
                                - ((j as f64 - 1.0) * std::f64::consts::PI
                                    / (nb_pnt_i as f64 - 1.0))
                                    .cos())
                                / 2.0);
                }
                self.my_par_cont.as_mut().unwrap()[i - 1].push(inter); // add a point
                if nb_pt_inter != 0 {
                    let mut l = 1usize;
                    while l <= nb_pt_inter {
                        // check if the point Inter is in the interval
                        // PntInter[i] PntInter[i+1] (duplicates problem).
                        if (inter > pnt_inter[i - 1][l - 1])
                            && (inter < pnt_inter[i - 1][l])
                        {
                            // leave the loop without storing the point.
                            break;
                        } else {
                            if l + 1 >= nb_pt_inter {
                                // one has parsed the entire table: the point
                                // does not belong to a common point interval.
                                if nb_pt_g1g1 != 0 {
                                    // if there exists an incompatible interval
                                    let mut k = 1usize;
                                    while k <= nb_pt_g1g1 {
                                        if (inter > pnt_g1g1[i - 1][k - 1])
                                            && (inter < pnt_g1g1[i - 1][k])
                                        {
                                            // to leave the loop — add points
                                            // of constraint G0.
                                            let c2d = c2d.as_ref().unwrap();
                                            let p2d = c2d.point_at(inter);
                                            let p3d = lin_cont.d0(inter);
                                            let pp = self
                                                .my_surf_init
                                                .as_ref()
                                                .unwrap()
                                                .point_at(p2d.x, p2d.y);
                                            let pdif = DVec3::new(
                                                -pp.x + p3d.x,
                                                -pp.y + p3d.y,
                                                -pp.z + p3d.z,
                                            );
                                            let pc = PinpointConstraint::new(p2d, pdif, 0, 0);
                                            self.my_plate.load_pinpoint(pc);
                                            break;
                                        } else {
                                            // the point does not belong to
                                            // interval G1.
                                            if k + 1 >= nb_pt_g1g1 {
                                                self.my_plate_cont.as_mut().unwrap()[i - 1]
                                                    .push(inter);
                                                // add the point
                                            }
                                        }
                                        k += 2;
                                    }
                                } else {
                                    self.my_plate_cont.as_mut().unwrap()[i - 1].push(inter);
                                    // add the point
                                }
                            }
                        }
                        l += 2;
                    }
                } else {
                    if nb_pt_g1g1 != 0 {
                        // there exist an incompatible interval
                        let mut k = 1usize;
                        while k <= nb_pt_g1g1 {
                            if (inter > pnt_g1g1[i - 1][k - 1]) && (inter < pnt_g1g1[i - 1][k]) {
                                // to leave the loop — add points of
                                // constraint G0.
                                let c2d = c2d.as_ref().unwrap();
                                let p2d = c2d.point_at(inter);
                                let p3d = lin_cont.d0(inter);
                                let pp = self
                                    .my_surf_init
                                    .as_ref()
                                    .unwrap()
                                    .point_at(p2d.x, p2d.y);
                                let pdif =
                                    DVec3::new(-pp.x + p3d.x, -pp.y + p3d.y, -pp.z + p3d.z);
                                let pc = PinpointConstraint::new(p2d, pdif, 0, 0);
                                self.my_plate.load_pinpoint(pc);
                                break;
                            } else {
                                // the point does not belong to interval G1.
                                if k + 1 >= nb_pt_g1g1 {
                                    self.my_plate_cont.as_mut().unwrap()[i - 1].push(inter);
                                    // add the point
                                }
                            }
                            k += 2;
                        }
                    } else {
                        // Geom2dAdaptor_Curve(LinCont->Curve2dOnSurf()).GetType()
                        // != GeomAbs_Circle — a null 2d curve reports
                        // OtherCurve, so "!= Circle" holds.
                        let type_not_circle = lin_cont
                            .curve2d_on_surf()
                            .map(|c| GCPntsCurve2d::new(&c).get_type() != CurveType::Circle)
                            .unwrap_or(true);
                        let exclude_extremities =
                            (!self.my_surf_init_is_give && type_not_circle)
                                || ((j > 1) && (j < nb_pnt_i));
                        if exclude_extremities {
                            self.my_plate_cont.as_mut().unwrap()[i - 1].push(inter);
                            // add the point
                        }
                    }
                }
            }
        }
    }

    /// OCCT CalculNbPtsInit (L2366-2400) — the number of points per curve,
    /// depending on the length, for the first iteration.
    pub(super) fn calcul_nb_pts_init(&mut self) {
        let mut len_t = 0.0f64;
        let ntlincont = self.my_lin_cont.len();
        let nt_point = self.my_nb_pts_on_cur * ntlincont as i32;

        for i in 1..=ntlincont {
            len_t += self.my_lin_cont[i - 1].length();
        }
        for i in 1..=ntlincont {
            let cont = self.my_lin_cont[i - 1].order();
            match cont {
                0 => {
                    // Case G0 *1.2
                    let nb = (1.2
                        * nt_point as f64
                        * (self.my_lin_cont[i - 1].length())
                        / len_t) as i32;
                    self.my_lin_cont[i - 1].set_nb_points(nb);
                }
                1 => {
                    // Case G1 *1
                    let nb = (nt_point as f64 * (self.my_lin_cont[i - 1].length()) / len_t) as i32;
                    self.my_lin_cont[i - 1].set_nb_points(nb);
                }
                2 => {
                    // Case G2 *0.7
                    let nb = (0.7
                        * nt_point as f64
                        * (self.my_lin_cont[i - 1].length())
                        / len_t) as i32;
                    self.my_lin_cont[i - 1].set_nb_points(nb);
                }
                _ => {}
            }
            if self.my_lin_cont[i - 1].nb_points() < 3 {
                self.my_lin_cont[i - 1].set_nb_points(3);
            }
        }
    }

    /// OCCT LoadCurve (L2407-2516) — starting from table myPlateCont load
    /// all the points noted in plate.
    pub(super) fn load_curve(&mut self, nb_boucle: i32, order_max: i32) {
        let ntlincont = self.my_lin_cont.len();

        for i in 1..=ntlincont {
            let cc = self.my_lin_cont[i - 1].clone();
            if cc.order() != -1 {
                let tang = cc.order().min(order_max);
                let nt = self.my_plate_cont.as_ref().unwrap()[i - 1].len();
                if tang != -1 {
                    for j in 1..=nt {
                        let u = self.my_plate_cont.as_ref().unwrap()[i - 1][j - 1];
                        // Loading of points G0 on boundaries.
                        let p3d = cc.d0(u);
                        let p2d;
                        if cc.projected_curve().is_some() {
                            let pc = cc.projected_curve().unwrap();
                            p2d = pc.value(u);
                        } else {
                            let c2d = cc.curve2d_on_surf();
                            if c2d.is_some() {
                                p2d = c2d.unwrap().point_at(u);
                            } else {
                                p2d = self.project_point(p3d);
                            }
                        }
                        let pp = self
                            .my_surf_init
                            .as_ref()
                            .unwrap()
                            .point_at(p2d.x, p2d.y);
                        let pdif = DVec3::new(-pp.x + p3d.x, -pp.y + p3d.y, -pp.z + p3d.z);
                        let pc = PinpointConstraint::new(p2d, pdif, 0, 0);
                        self.my_plate.load_pinpoint(pc);

                        // Loading of points G1.
                        if tang == 1 {
                            // ==1
                            let (_pp, v1, v2) = cc.d1(u);
                            let (_pp2, v3, v4) = self
                                .my_surf_init
                                .as_ref()
                                .unwrap()
                                .derivatives(p2d.x, p2d.y);

                            let d1final = PlateD1::new(v1, v2);
                            let d1init = PlateD1::new(v3, v4);
                            if !self.my_free {
                                let gcc = GtoCConstraint::new(p2d, &d1init, &d1final);
                                self.my_plate.load_gto_c(&gcc);
                            } else if nb_boucle == 1 {
                                let free_gcc =
                                    FreeGtoCConstraint::new(p2d, &d1init, &d1final, 1.0, 0);
                                self.my_plate.load_free_gto_c(&free_gcc);
                            } else {
                                // Normal = V1 ^ V2; (not normalized)
                                let mut normal = v1.cross(v2);
                                let norm = normal.length();
                                if norm > 1.0e-12 {
                                    normal /= norm;
                                }
                                let der_plate_u = self
                                    .my_prev_plate
                                    .evaluate_derivative(p2d, 1, 0);
                                let der_plate_v = self
                                    .my_prev_plate
                                    .evaluate_derivative(p2d, 0, 1);

                                let du = der_plate_u
                                    - normal * (-(v3 + der_plate_u).dot(normal));
                                let dv = der_plate_v
                                    - normal * (-(v4 + der_plate_v).dot(normal));
                                let pin_u = PinpointConstraint::new(p2d, du, 1, 0);
                                let pin_v = PinpointConstraint::new(p2d, dv, 0, 1);
                                self.my_plate.load_pinpoint(pin_u);
                                self.my_plate.load_pinpoint(pin_v);
                            }
                        }
                        // Loading of points G2.
                        if tang == 2 {
                            // ==2
                            let (_pp, v1, v2, v5, v6, v7) = cc.d2(u);
                            let (_pp2, v3, v4, v8, v9, v10) = self
                                .my_surf_init
                                .as_ref()
                                .unwrap()
                                .derivatives2(p2d.x, p2d.y);

                            let d1final = PlateD1::new(v1, v2);
                            let d1init = PlateD1::new(v3, v4);
                            let d2final = PlateD2::new(v5, v6, v7);
                            let d2init = PlateD2::new(v8, v9, v10);
                            let gcc =
                                GtoCConstraint::new_g2(p2d, &d1init, &d1final, &d2init, &d2final);
                            self.my_plate.load_gto_c(&gcc);
                        }
                    }
                }
            }
        }
    }

    /// OCCT VerifSurface (L2588-2732) — the error evaluation and the point
    /// count increase loop; returns false when another iteration is needed.
    pub(super) fn verif_surface(&mut self, nb_boucle: i32) -> bool {
        // Calculate errors.
        let ntlincont = self.my_lin_cont.len();
        let mut result = true;

        // variable for error calculation
        self.my_g0_error = 0.0;
        self.my_g1_error = 0.0;
        self.my_g2_error = 0.0;

        for i in 1..=ntlincont {
            let lin_cont = self.my_lin_cont[i - 1].clone();
            if lin_cont.order() != -1 {
                let mut nb_pts_i = self.my_par_cont.as_ref().unwrap()[i - 1].len();
                if nb_pts_i < 3 {
                    nb_pts_i = 4;
                }
                let mut tdist = vec![0.0f64; nb_pts_i - 1];
                let mut tang = vec![0.0f64; nb_pts_i - 1];
                let mut tcourb = vec![0.0f64; nb_pts_i - 1];

                self.ecart_contraintes_mil(i, &mut tdist, &mut tang, &mut tcourb);

                let mut diff_dist_max = 0.0f64;
                let mut diff_ang_max = 0.0f64;
                let mut n_diff_dist = 0;
                let mut n_diff_ang = 0;

                for j in 1..nb_pts_i {
                    if tdist[j - 1] > self.my_g0_error {
                        self.my_g0_error = tdist[j - 1];
                    }
                    if tang[j - 1] > self.my_g1_error {
                        self.my_g1_error = tang[j - 1];
                    }
                    if tcourb[j - 1] > self.my_g2_error {
                        self.my_g2_error = tcourb[j - 1];
                    }
                    let u;
                    if self.my_par_cont.as_ref().unwrap()[i - 1].len() > 3 {
                        u = (self.my_par_cont.as_ref().unwrap()[i - 1][j - 1]
                            + self.my_par_cont.as_ref().unwrap()[i - 1][j])
                            / 2.0;
                    } else {
                        u = lin_cont.first_parameter()
                            + (lin_cont.last_parameter() - lin_cont.first_parameter())
                                * (j as f64 - 1.0)
                                / (nb_pts_i as f64 - 2.0);
                    }
                    let mut diff_dist = tdist[j - 1] - lin_cont.g0_criterion(u);
                    let diff_ang;
                    if lin_cont.order() > 0 {
                        diff_ang = tang[j - 1] - lin_cont.g1_criterion(u);
                    } else {
                        diff_ang = 0.0;
                    }
                    // find the maximum variation of error and calculate the
                    // average.
                    if diff_dist > 0.0 {
                        diff_dist /= lin_cont.g0_criterion(u);
                        if diff_dist > diff_dist_max {
                            diff_dist_max = diff_dist;
                        }
                        n_diff_dist += 1;
                    } else if (diff_ang > 0.0) && (lin_cont.order() == 1) {
                        let diff_ang = diff_ang / lin_cont.g1_criterion(u);
                        if diff_ang > diff_ang_max {
                            diff_ang_max = diff_ang;
                        }
                        n_diff_ang += 1;
                    }
                    let _ = diff_ang;
                }

                if n_diff_dist > 0 {
                    // at least one point is not acceptable in G0
                    let mut coef;
                    if lin_cont.order() == 0 {
                        coef = 0.6 * (diff_dist_max + 7.4).ln();
                        // 7.4 corresponds to the calculation of min.
                        // coefficient = 1.2 is e^1.2/0.6
                    } else {
                        coef = (diff_dist_max + 3.3).ln();
                    }
                    // 3.3 corresponds to calculation of min. coefficient =
                    // 1.2 donc e^1.2
                    if coef > 3.0 {
                        coef = 3.0;
                    }
                    // experimentally after the coefficient becomes bad for
                    // L cases
                    if (nb_boucle > 1) && (diff_dist_max > 2.0) {
                        coef = 1.6;
                    }

                    if lin_cont.nb_points() as f64 >= (lin_cont.nb_points() as f64 * coef).floor()
                    {
                        coef = 2.0; // to provide increase of the number of points
                    }

                    self.my_lin_cont[i - 1].set_nb_points((lin_cont.nb_points() as f64 * coef) as i32);
                    result = false;
                } else if n_diff_ang > 0 {
                    // at least 1 point is not acceptable in G1
                    let mut coef = 1.5f64;
                    if ((lin_cont.nb_points() + 1) as f64)
                        >= (lin_cont.nb_points() as f64 * coef).floor()
                    {
                        coef = 2.0;
                    }

                    self.my_lin_cont[i - 1].set_nb_points((lin_cont.nb_points() as f64 * coef) as i32);
                    result = false;
                }
            }
        }
        if !result {
            if self.my_free && nb_boucle == 1 {
                self.my_prev_plate = self.my_plate.clone();
            }
            self.my_plate.init();
        }
        result
    }

    /// OCCT EcartContraintesMil (L751-866) — the G0/G1/G2 gaps at the middle
    /// of each constraint interval of curve c.
    pub(super) fn ecart_contraintes_mil(
        &self,
        c: usize,
        d: &mut [f64],
        an: &mut [f64],
        courb: &mut [f64],
    ) {
        let mut nb_pt = self.my_par_cont.as_ref().unwrap()[c - 1].len();
        if nb_pt < 3 {
            nb_pt = 4;
        } else {
            nb_pt = self.my_par_cont.as_ref().unwrap()[c - 1].len();
        }
        let lin_cont = self.my_lin_cont[c - 1].clone();
        match lin_cont.order() {
            0 => {
                for i in 1..nb_pt {
                    let u = (self.my_par_cont.as_ref().unwrap()[c - 1][i - 1]
                        + self.my_par_cont.as_ref().unwrap()[c - 1][i])
                        / 2.0;
                    let pi = lin_cont.d0(u);
                    let p2d;
                    if lin_cont.projected_curve().is_some() {
                        p2d = lin_cont.projected_curve().unwrap().value(u);
                    } else {
                        let c2d = lin_cont.curve2d_on_surf();
                        if c2d.is_some() {
                            p2d = c2d.unwrap().point_at(u);
                        } else {
                            p2d = self.project_point(pi);
                        }
                    }
                    let pf = self
                        .my_geom_plate_surface
                        .as_ref()
                        .unwrap()
                        .eval_d0(p2d.x, p2d.y);
                    an.fill(0.0);
                    courb.fill(0.0);
                    d[i - 1] = pf.distance(pi);
                }
            }
            1 => {
                for i in 1..nb_pt {
                    let u = (self.my_par_cont.as_ref().unwrap()[c - 1][i - 1]
                        + self.my_par_cont.as_ref().unwrap()[c - 1][i])
                        / 2.0;
                    let (pi, v1i, v2i) = lin_cont.d1(u);
                    let p2d;
                    if lin_cont.projected_curve().is_some() {
                        p2d = lin_cont.projected_curve().unwrap().value(u);
                    } else {
                        let c2d = lin_cont.curve2d_on_surf();
                        if c2d.is_some() {
                            p2d = c2d.unwrap().point_at(u);
                        } else {
                            p2d = self.project_point(pi);
                        }
                    }
                    let (pf, v1f, v2f) = self
                        .my_geom_plate_surface
                        .as_ref()
                        .unwrap()
                        .eval_d1(p2d.x, p2d.y);
                    d[i - 1] = pf.distance(pi);
                    let v3i = v1i.cross(v2i);
                    let v3f = v1f.cross(v2f);
                    let angle = vec3_angle(v3f, v3i);
                    if angle > (std::f64::consts::PI / 2.0) {
                        an[i - 1] = std::f64::consts::PI - angle;
                    } else {
                        an[i - 1] = angle;
                    }
                    courb.fill(0.0);
                }
            }
            2 => {
                // occ::handle<Geom_Surface> Splate(myGeomPlateSurface);
                // LocalAnalysis_SurfaceContinuity CG2;
                // ... CG2.ComputeAnalysis(Prop, myLinCont->Value(c)->LPropSurf(U),
                //                         GeomAbs_G2);
                //
                // GAP leaf: LocalAnalysis_SurfaceContinuity +
                // GeomLProp_SLProps (untranslated packages); the anchor
                // preserves the dependency failure path.
                unimplemented!(
                    "EcartContraintesMil case 2 needs LocalAnalysis_SurfaceContinuity + GeomLProp_SLProps (untranslated)"
                );
            }
            _ => {}
        }
    }

    /// OCCT Disc2dContour (L871-1019) — the sampling in "cosine" plus 3
    /// points on each interval (2d contour).
    pub fn disc2d_contour(&self, seq2d: &mut Vec<DVec2>) {
        // initialization
        seq2d.clear();

        // sampling in "cosine" + 3 points on each interval
        let ntcurve = self.my_lin_cont.len();
        let ntpntcont = self.my_pnt_cont.len();

        for i in 1..=ntpntcont {
            if self.my_pnt_cont[i - 1].order() != -1 {
                let p2d = self.my_pnt_cont[i - 1].pnt2d_on_surf();
                seq2d.push(p2d);
            }
        }
        for i in 1..=ntcurve {
            let lin_cont = self.my_lin_cont[i - 1].clone();
            if lin_cont.order() != -1 {
                let nb_pt = self.my_par_cont.as_ref().unwrap()[i - 1].len();
                // first point of constraint (j=0)
                let p2d;
                if lin_cont.projected_curve().is_some() {
                    p2d = lin_cont
                        .projected_curve()
                        .unwrap()
                        .value(self.my_par_cont.as_ref().unwrap()[i - 1][0]);
                } else {
                    let c2d = lin_cont.curve2d_on_surf();
                    if c2d.is_some() {
                        p2d = c2d
                            .unwrap()
                            .point_at(self.my_par_cont.as_ref().unwrap()[i - 1][0]);
                    } else {
                        let pp = lin_cont.d0(self.my_par_cont.as_ref().unwrap()[i - 1][0]);
                        p2d = self.project_point(pp);
                    }
                }

                seq2d.push(p2d);
                for j in 2..nb_pt {
                    let uj = self.my_par_cont.as_ref().unwrap()[i - 1][j - 1];
                    let ujp1 = self.my_par_cont.as_ref().unwrap()[i - 1][j];
                    // point 1/4 previous
                    let p2d;
                    if lin_cont.projected_curve().is_some() {
                        p2d = lin_cont.projected_curve().unwrap().value((ujp1 + 3.0 * uj) / 4.0);
                    } else {
                        let c2d = lin_cont.curve2d_on_surf();
                        if c2d.is_some() {
                            p2d = c2d.unwrap().point_at((ujp1 + 3.0 * uj) / 4.0);
                        } else {
                            let pp = lin_cont.d0((ujp1 + 3.0 * uj) / 4.0);
                            p2d = self.project_point(pp);
                        }
                    }
                    seq2d.push(p2d);
                    // point 1/2 previous
                    let p2d;
                    if lin_cont.projected_curve().is_some() {
                        p2d = lin_cont.projected_curve().unwrap().value((ujp1 + uj) / 2.0);
                    } else {
                        let c2d = lin_cont.curve2d_on_surf();
                        if c2d.is_some() {
                            p2d = c2d.unwrap().point_at((ujp1 + uj) / 2.0);
                        } else {
                            let pp = lin_cont.d0((ujp1 + uj) / 2.0);
                            p2d = self.project_point(pp);
                        }
                    }
                    seq2d.push(p2d);
                    // point 3/4 previous
                    let p2d;
                    if lin_cont.projected_curve().is_some() {
                        p2d = lin_cont.projected_curve().unwrap().value((3.0 * ujp1 + uj) / 4.0);
                    } else {
                        let c2d = lin_cont.curve2d_on_surf();
                        if c2d.is_some() {
                            p2d = c2d.unwrap().point_at((3.0 * ujp1 + uj) / 4.0);
                        } else {
                            let pp = lin_cont.d0((3.0 * ujp1 + uj) / 4.0);
                            p2d = self.project_point(pp);
                        }
                    }
                    seq2d.push(p2d);
                    // current constraint point
                    let p2d;
                    if lin_cont.projected_curve().is_some() {
                        p2d = lin_cont.projected_curve().unwrap().value(ujp1);
                    } else {
                        let c2d = lin_cont.curve2d_on_surf();
                        if c2d.is_some() {
                            p2d = c2d.unwrap().point_at(ujp1);
                        } else {
                            let pp = lin_cont.d0(ujp1);
                            p2d = self.project_point(pp);
                        }
                    }
                    seq2d.push(p2d);
                }
            }
        }
    }

    /// OCCT Disc3dContour (L1024-1163) — the sampling in "cosine" plus 3
    /// points on each interval (3d contour; iordre 0 = positions, 1 =
    /// cross-product normals).
    pub fn disc3d_contour(&self, iordre: i32, seq3d: &mut Vec<DVec3>) {
        // initialization
        seq3d.clear();
        // sampling in "cosine" + 3 points on each interval
        let ntcurve = self.my_lin_cont.len();
        let ntpntcont = self.my_pnt_cont.len();

        for i in 1..=ntpntcont {
            if self.my_pnt_cont[i - 1].order() != -1 {
                if iordre == 0 {
                    let p3d = self.my_pnt_cont[i - 1].d0();
                    seq3d.push(p3d);
                } else {
                    let (_p3d, v1h, v2h) = self.my_pnt_cont[i - 1].d1();
                    let v3h = v1h.cross(v2h);
                    seq3d.push(v3h);
                }
            }
        }

        for i in 1..=ntcurve {
            if self.my_lin_cont[i - 1].order() != -1 {
                let nb_pt = self.my_par_cont.as_ref().unwrap()[i - 1].len();
                // first constraint point (j=0)
                if iordre == 0 {
                    let p3d =
                        self.my_lin_cont[i - 1].d0(self.my_par_cont.as_ref().unwrap()[i - 1][0]);
                    seq3d.push(p3d);
                } else {
                    let (_p3d, v1h, v2h) = self.my_lin_cont[i - 1]
                        .d1(self.my_par_cont.as_ref().unwrap()[i - 1][0]);
                    let v3h = v1h.cross(v2h);
                    seq3d.push(v3h);
                }

                for j in 2..nb_pt {
                    let uj = self.my_par_cont.as_ref().unwrap()[i - 1][j - 1];
                    let ujp1 = self.my_par_cont.as_ref().unwrap()[i - 1][j];
                    if iordre == 0 {
                        // point 1/4 previous
                        let p3d = self.my_lin_cont[i - 1].d0((ujp1 + 3.0 * uj) / 4.0);
                        seq3d.push(p3d);
                        // point 1/2 previous
                        let p3d = self.my_lin_cont[i - 1].d0((ujp1 + uj) / 2.0);
                        seq3d.push(p3d);
                        // point 3/4 previous
                        let p3d = self.my_lin_cont[i - 1].d0((3.0 * ujp1 + uj) / 4.0);
                        seq3d.push(p3d);
                        // current constraint point
                        let p3d = self.my_lin_cont[i - 1].d0(ujp1);
                        seq3d.push(p3d);
                    } else {
                        // point 1/4 previous
                        let (_p3d, v1h, v2h) = self.my_lin_cont[i - 1].d1((ujp1 + 3.0 * uj) / 4.0);
                        let v3h = v1h.cross(v2h);
                        seq3d.push(v3h);
                        // point 1/2 previous
                        let (_p3d, v1h, v2h) = self.my_lin_cont[i - 1].d1((ujp1 + uj) / 2.0);
                        let v3h = v1h.cross(v2h);
                        seq3d.push(v3h);
                        // point 3/4 previous
                        let (_p3d, v1h, v2h) = self.my_lin_cont[i - 1].d1((3.0 * ujp1 + uj) / 4.0);
                        let v3h = v1h.cross(v2h);
                        seq3d.push(v3h);
                        // current constraint point
                        let (_p3d, v1h, v2h) = self.my_lin_cont[i - 1].d1(ujp1);
                        let v3h = v1h.cross(v2h);
                        seq3d.push(v3h);
                    }
                }
            }
        }
    }

    /// OCCT G0Error(Index) (L1260-1278).
    pub fn g0_error_index(&mut self, index: i32) -> f64 {
        let n = self.my_nb_pts_on_cur as usize;
        let mut tdistance = vec![0.0f64; n];
        let mut tangle = vec![0.0f64; n];
        let mut tcurvature = vec![0.0f64; n];
        self.ecart_contraintes_mil(index as usize, &mut tdistance, &mut tangle, &mut tcurvature);
        let mut max_distance = 0.0f64;
        for i in 1..=n {
            if tdistance[i - 1] > max_distance {
                max_distance = tdistance[i - 1];
            }
        }
        max_distance
    }

    /// OCCT G1Error(Index) (L1282-1300).
    pub fn g1_error_index(&mut self, index: i32) -> f64 {
        let n = self.my_nb_pts_on_cur as usize;
        let mut tdistance = vec![0.0f64; n];
        let mut tangle = vec![0.0f64; n];
        let mut tcurvature = vec![0.0f64; n];
        self.ecart_contraintes_mil(index as usize, &mut tdistance, &mut tangle, &mut tcurvature);
        let mut max_angle = 0.0f64;
        for i in 1..=n {
            if tangle[i - 1] > max_angle {
                max_angle = tangle[i - 1];
            }
        }
        max_angle
    }

    /// OCCT G2Error(Index) (L1304-1322).
    pub fn g2_error_index(&mut self, index: i32) -> f64 {
        let n = self.my_nb_pts_on_cur as usize;
        let mut tdistance = vec![0.0f64; n];
        let mut tangle = vec![0.0f64; n];
        let mut tcurvature = vec![0.0f64; n];
        self.ecart_contraintes_mil(index as usize, &mut tdistance, &mut tangle, &mut tcurvature);
        let mut max_curvature = 0.0f64;
        for i in 1..=n {
            if tcurvature[i - 1] > max_curvature {
                max_curvature = tcurvature[i - 1];
            }
        }
        max_curvature
    }
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx) — the angle in [0, PI]; raises
/// VectorWithNullMagnitude on a null vector.
pub(crate) fn vec3_angle(v: DVec3, the_ref: DVec3) -> f64 {
    if v.length() <= f64::MIN_POSITIVE || the_ref.length() <= f64::MIN_POSITIVE {
        panic!("gp_VectorWithNullMagnitude");
    }
    (v.dot(the_ref) / (v.length() * the_ref.length())).acos()
}
