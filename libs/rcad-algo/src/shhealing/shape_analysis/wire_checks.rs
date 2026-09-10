//! OCCT ShapeAnalysis_Wire — the per-edge checks and the tail/loop checks
//! (`ShapeAnalysis_Wire.cxx` L896-2599).  The impl submodule of
//! [`super::wire::ShapeAnalysisWire`] (the OCCT continued-file convention).
//!
//! GAP (wire.rs header bridge #5): `Geom2dInt_GInter` curve-curve
//! intersection — the re-host keeps the OCCT `!IsDone()` branch live.
//! GAP (bridge #8): `ShapeAnalysis_TransferParametersProj` and the GProp
//! property machinery of CheckTail / CheckSmallArea.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, CurveEval};
use rcad_kernel::precision::PCONFUSION;
use rcad_kernel::topo::topods::{BRep, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use super::analysis;
use super::curve::ShapeAnalysisCurve;
use super::edge::ShapeAnalysisEdge;
use super::wire::{
    brep_tool_degenerated, brep_tool_pnt, brep_tool_same_parameter, brep_tool_tolerance,
    brep_tools_compare, ShapeAnalysisWire,
};
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};

// OCCT gp.hxx L59-60: gp::Resolution() = RealSmall() = DBL_MIN.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT Geom2dInt_GInter (the curve-curve intersection) re-host (wire.rs
/// header bridge #5): the kernel carries the line-curve arm only; the
/// general curve-curve `Perform` is not translated, so the result is the
/// OCCT `!IsDone()` state — every consumer takes the OCCT failure branch.
struct GInterCurveCurve;

impl GInterCurveCurve {
    fn new() -> Self {
        GInterCurveCurve
    }
    /// OCCT IsDone() — false from the GAP.
    fn is_done(&self) -> bool {
        false
    }
}

impl ShapeAnalysisWire {
    /// OCCT CheckDegenerated(num, p2d1, p2d2) (cxx L896-1114).
    pub fn check_degenerated_full(
        &mut self,
        brep: &mut BRep,
        num: i32,
        p2d1: &mut DVec2,
        p2d2: &mut DVec2,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 1 {
            return false;
        }

        let nb_edges = self.nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let n3 = if n2 < nb_edges { n2 + 1 } else { 1 };
        let e1 = self.my_wire.as_ref().unwrap().edge(n1);
        let e2 = self.my_wire.as_ref().unwrap().edge(n2);
        let e3 = self.my_wire.as_ref().unwrap().edge(n3);

        let sae = ShapeAnalysisEdge::new();

        // skip if edge is already marked as degenerated and has pcurve
        if brep_tool_degenerated(&e2) && sae.has_pcurve_face(brep, &e2, &self.my_face) {
            // skl 30.12.2004 for OCC7630 - we have to check pcurve
            if sae.has_pcurve_face(brep, &e1, &self.my_face)
                && sae.has_pcurve_face(brep, &e3, &self.my_face)
            {
                let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
                let (mut fp, mut lp) = (0.0f64, 0.0f64);
                sae.pcurve_face(brep, &e2, &self.my_face, &mut c2d, &mut fp, &mut lp, true);
                let (p21, p22) = {
                    let c = c2d.as_ref().unwrap();
                    (
                        Curve2dEval::point_at(c, fp),
                        Curve2dEval::point_at(c, lp),
                    )
                };
                sae.pcurve_face(brep, &e1, &self.my_face, &mut c2d, &mut fp, &mut lp, true);
                let p12 = Curve2dEval::point_at(c2d.as_ref().unwrap(), lp);
                sae.pcurve_face(brep, &e3, &self.my_face, &mut c2d, &mut fp, &mut lp, true);
                let p31 = Curve2dEval::point_at(c2d.as_ref().unwrap(), fp);
                if (p12.distance(p31) - p21.distance(p22)).abs() > 2.0 * PCONFUSION {
                    // pcurve is bad => we can remove this edge in ShapeFix
                    // if set needed status
                    self.my_status = encode_status(ShapeExtendStatus::Fail2);
                }
            }
            return false;
        }

        // pdn allows to insert two sequences of degenerated edges (on
        // separate bounds of surfaces)
        if n1 != n2 && brep_tool_degenerated(&e1) && !sae.has_pcurve_face(brep, &e1, &self.my_face)
        {
            //: abv 13.05.02: OCC320 - fail (to remove edge) if two consecutive
            //: degenerated edges w/o pcurves
            if brep_tool_degenerated(&e2) {
                self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            }
            return false;
        }

        let vp = sae.first_vertex(brep, &e1); //: i9
        let v0 = sae.last_vertex(brep, &e1);
        let v1 = sae.first_vertex(brep, &e2);
        let v2 = sae.last_vertex(brep, &e2);

        if vp.is_null() || v0.is_null() || v1.is_null() || v2.is_null() {
            return false;
        }

        let pp = brep_tool_pnt(&vp); //: i9
        let p0 = brep_tool_pnt(&v0);
        let p1 = brep_tool_pnt(&v1);
        let p2 = brep_tool_pnt(&v2);
        let mut par1 = 0.0f64;
        let mut par2 = 0.0f64;
        let mut lack = false;
        let mut dgnr = false;
        // pdn 12.03.99 minimal value processing first
        let v1_tol = brep_tool_tolerance(&v1);
        let prec_first = self.my_precision.min(v1_tol);
        let prec_fin = self.my_precision.max(v1_tol);
        let prec_vtx = if self.my_precision < v1_tol {
            2.0 * prec_fin
        } else {
            prec_fin
        };
        //  forward : si Edge <num> FWD/REV. Si LACK, toujours True
        let mut forward = e2.orientation == Orientation::Forward;
        //  FIX FEV 1998 : recompute singularity according precision

        if p1.distance(p2) <= prec_first {
            // edge DGNR
            let surf = self.my_surf.as_mut().unwrap();
            dgnr = surf.degenerated_values(p1, prec_vtx, p2d1, p2d2, &mut par1, &mut par2, forward); // smh#9
            if dgnr {
                // abv 24 Feb 00: trj3_as1-ac-214.stp #6065: avoid making
                // closed edge degenerated
                if let Some(c3d) = (match e2.data.as_ref() {
                    TShape::Edge(ed) => ed.curve.clone(),
                    _ => None,
                }) {
                    let p = CurveEval::point_at(&c3d, 0.5 * (CurveEval::default_domain(&c3d)[0] + CurveEval::default_domain(&c3d)[1]));
                    if p.distance_squared(p1) > prec_vtx * prec_vtx {
                        dgnr = false;
                    }
                }
            }
        }
        if !dgnr {
            //: i9 abv 23 Sep 98: CTS20315-2 #63231: check that previous edge
            //: is not degenerated
            if n1 != n2
                && p1.distance(pp) <= prec_first
                && self
                    .my_surf
                    .as_mut()
                    .unwrap()
                    .is_degenerated(pp, prec_first)
                && !brep_tool_degenerated(&e1)
            {
                return false;
            }
            // rln S4135 / :45 / :51 — the singularity walk
            if p0.distance(p1) <= prec_fin {
                // ou DGNR manquante ?
                // #77 rln S4135: using singularity which has minimum gap
                let mut ind_min: i32 = -1;
                let mut gap_min2 = f64::MAX;
                let surf = self.my_surf.as_mut().unwrap();
                let nb_sing = surf.nb_singularities(prec_vtx);
                let mut tmp_preci = 0.0f64;
                let mut tmp_p3d = DVec3::ZERO;
                let mut tmp_uiso_deg = false;
                for i in 1..=nb_sing {
                    surf.singularity(
                        i,
                        &mut tmp_preci,
                        &mut tmp_p3d,
                        p2d1,
                        p2d2,
                        &mut par1,
                        &mut par2,
                        &mut tmp_uiso_deg,
                    );
                    let gap2 = p1.distance_squared(tmp_p3d);
                    if gap2 <= prec_vtx * prec_vtx && gap_min2 > gap2 {
                        gap_min2 = gap2;
                        ind_min = i;
                    }
                }
                if ind_min >= 1 {
                    let surf = self.my_surf.as_mut().unwrap();
                    surf.singularity(
                        ind_min,
                        &mut tmp_preci,
                        &mut tmp_p3d,
                        p2d1,
                        p2d2,
                        &mut par1,
                        &mut par2,
                        &mut tmp_uiso_deg,
                    );
                    lack = true;
                }
            }
        }

        // there you go, we either have degenerate or lack
        if !lack && !dgnr {
            //: abv 29.08.01: if singularity not detected but edge is marked
            // as degenerated, report fail
            if brep_tool_degenerated(&e2) && !sae.has_pcurve_face(brep, &e2, &self.my_face) {
                self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            }
            return false;
        }

        // OK, degenerated case detected; we will find its start and end in 2d

        if lack {
            forward = true;
        }
        let _ = forward;

        //: 24 by abv 28 Nov 97: make degenerative pcurve parametrized exactly
        // from end of pcurve of the previous edge to the start of the next one
        if lack || n1 != n2 {
            //: i8 abv 18 Sep 98
            let (mut a, mut b) = (0.0f64, 0.0f64);
            let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
            if sae.pcurve_face(brep, &e1, &self.my_face, &mut c2d, &mut a, &mut b, true) {
                *p2d1 = Curve2dEval::point_at(c2d.as_ref().unwrap(), b);
            } else {
                self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            }
            // pdn pcurves (fixing regression in f0 in degenerated case)
            let mut c2d_b: Option<rcad_kernel::geom::Curve2d> = None;
            let (mut a2, mut b2) = (0.0f64, 0.0f64);
            if sae.pcurve_face(
                brep,
                if dgnr { &e3 } else { &e2 },
                &self.my_face,
                &mut c2d_b,
                &mut a2,
                &mut b2,
                true,
            ) {
                *p2d2 = Curve2dEval::point_at(c2d_b.as_ref().unwrap(), a2);
            } else {
                self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            }
        }

        // #84 rln 18.03.99 if pcurve is not degenerate anymore, the fix is
        // postponned to ShapeFix_Wire::FixLacking
        if !self
            .my_surf
            .as_ref()
            .unwrap()
            .is_degenerated_uv(*p2d1, *p2d2, prec_vtx, 10.0)
        {
            //: s1 abv 22 Apr 99: PRO7226 #489490 //smh#9
            //: abv 24.05.02: OCC320
            if brep_tool_degenerated(&e2) {
                self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            }
            return false;
        }

        // added by rln 18/12/97 CSR# CTS18544 entity 13638
        // the situation when degenerated edge already exists but flag is not
        // set (i.e. the parametric space is closed)
        let ads_u = self.my_surf.as_ref().unwrap().surface().u_resolution(self.my_precision);
        let ads_v = self.my_surf.as_ref().unwrap().surface().v_resolution(self.my_precision);
        let max = ads_u.max(ads_v);
        if p2d1.distance(*p2d2) <= max + GP_RESOLUTION {
            return false;
        }

        self.my_status = encode_status(if dgnr {
            ShapeExtendStatus::Done2
        } else {
            ShapeExtendStatus::Done1
        });
        true
    }

    /// OCCT CheckDegenerated(num) (cxx L1116-1122).
    pub fn check_degenerated(&mut self, brep: &mut BRep, num: i32) -> bool {
        let mut p2d1 = DVec2::ZERO;
        let mut p2d2 = DVec2::ZERO;
        self.check_degenerated_full(brep, num, &mut p2d1, &mut p2d2)
    }

    /// OCCT CheckGap3d(num) (cxx L1124-1154).
    pub fn check_gap3d(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // szv#4:S4163:12Mar99 optimized
        if !self.is_loaded() || self.nb_edges() < 1 {
            return false; // szvsh was nbedges < 2
        }
        let nb_edges = self.nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let e1 = self.my_wire.as_ref().unwrap().edge(n1);
        let e2 = self.my_wire.as_ref().unwrap().edge(n2);
        let (mut uf1, mut ul1, mut uf2, mut ul2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut c1, mut c2): (
            Option<rcad_kernel::geom::Curve3>,
            Option<rcad_kernel::geom::Curve3>,
        ) = (None, None);
        let sae = ShapeAnalysisEdge::new();
        if !sae.curve3d(brep, &e1, &mut c1, &mut uf1, &mut ul1, false)
            || !sae.curve3d(brep, &e2, &mut c2, &mut uf2, &mut ul2, false)
        {
            self.my_status = encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let p1 = CurveEval::point_at(c1.as_ref().unwrap(), ul1);
        let p2 = CurveEval::point_at(c2.as_ref().unwrap(), uf2);
        self.my_min3d = p1.distance(p2);
        self.my_max3d = self.my_min3d;
        if self.my_min3d > self.my_precision {
            self.my_status = encode_status(ShapeExtendStatus::Done1);
        }
        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckGap2d(num) (cxx L1156-1194).
    pub fn check_gap2d(&mut self, brep: &mut BRep, num: i32) -> bool {
        if self.my_face.is_null() || self.my_surf.is_none() {
            self.my_status = encode_status(ShapeExtendStatus::Fail1);
            return false;
        }

        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 1 {
            return false; // szvsh was nbedges < 2
        }
        let nb_edges = self.nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let e1 = self.my_wire.as_ref().unwrap().edge(n1);
        let e2 = self.my_wire.as_ref().unwrap().edge(n2);
        let (mut uf1, mut ul1, mut uf2, mut ul2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut c1, mut c2): (
            Option<rcad_kernel::geom::Curve2d>,
            Option<rcad_kernel::geom::Curve2d>,
        ) = (None, None);
        let sae = ShapeAnalysisEdge::new();
        if !sae.pcurve_face(brep, &e1, &self.my_face, &mut c1, &mut uf1, &mut ul1, false)
            || !sae.pcurve_face(brep, &e2, &self.my_face, &mut c2, &mut uf2, &mut ul2, false)
        {
            self.my_status = encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let p1 = Curve2dEval::point_at(c1.as_ref().unwrap(), ul1);
        let p2 = Curve2dEval::point_at(c2.as_ref().unwrap(), uf2);
        self.my_min2d = p1.distance(p2);
        self.my_max2d = self.my_min2d;
        let surf = self.my_surf.as_ref().unwrap().surface();
        let sa_u = surf.u_resolution(self.my_precision);
        let sa_v = surf.v_resolution(self.my_precision);
        if self.my_min2d > (sa_u.max(sa_v) + PCONFUSION) {
            self.my_status = encode_status(ShapeExtendStatus::Done1);
        }
        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckCurveGap(num) (cxx L1196-1246).
    pub fn check_curve_gap(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() < 1 {
            return false;
        }
        let n = if num > 0 { num } else { self.nb_edges() };
        let e = self.my_wire.as_ref().unwrap().edge(n);
        let (mut cuf, mut cul, mut pcuf, mut pcul) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut c: Option<rcad_kernel::geom::Curve3> = None;
        let sae = ShapeAnalysisEdge::new();
        if !sae.curve3d(brep, &e, &mut c, &mut cuf, &mut cul, false) {
            self.my_status = encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        let mut pc: Option<rcad_kernel::geom::Curve2d> = None;
        if !sae.pcurve_face(brep, &e, &self.my_face, &mut pc, &mut pcuf, &mut pcul, false) {
            self.my_status = encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        // Adaptor3d_CurveOnSurface ACS(AC, AS): the pcurve lifted on the
        // surface (the SurfaceEval value pair; bridge #2 of wire.rs).
        let surf = self.my_surf.as_ref().unwrap().surface().clone();
        let pc = pc.unwrap();
        let nbp = 45;
        let mut dist;
        let mut maxdist = 0.0f64;
        for i in 0..nbp {
            let cpnt = CurveEval::point_at(c.as_ref().unwrap(), cuf + i as f64 * (cul - cuf) / (nbp - 1) as f64);
            let uv = Curve2dEval::point_at(&pc, pcuf + i as f64 * (pcul - pcuf) / (nbp - 1) as f64);
            let pcpnt = rcad_kernel::geom::SurfaceEval::point_at(&surf, uv.x, uv.y);
            dist = cpnt.distance_squared(pcpnt);
            if maxdist < dist {
                maxdist = dist;
            }
        }
        self.my_min3d = maxdist.sqrt();
        self.my_max3d = self.my_min3d;
        if self.my_min3d > self.my_precision {
            self.my_status = encode_status(ShapeExtendStatus::Done1);
        }
        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckSelfIntersectingEdge(num, points2d, points3d)
    /// (cxx L1269-1344) — the intersection points are collected as
    /// (param-on-first, point) pairs (the IntRes2d_IntersectionPoint
    /// projection).
    pub fn check_self_intersecting_edge_full(
        &mut self,
        brep: &mut BRep,
        num: i32,
        _points2d: &mut Vec<DVec2>,
        _points3d: &mut Vec<DVec3>,
    ) -> bool {
        _points2d.clear();
        _points3d.clear();
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let edge = self.my_wire.as_ref().unwrap().edge(if num > 0 {
            num
        } else {
            self.nb_edges()
        });
        let sae = ShapeAnalysisEdge::new();

        let (mut a, mut b) = (0.0f64, 0.0f64);
        let mut crv: Option<rcad_kernel::geom::Curve2d> = None;
        if !sae.pcurve_face(brep, &edge, &self.my_face, &mut crv, &mut a, &mut b, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        if (a - b).abs() <= PCONFUSION {
            return false;
        }

        let tolint = 1.0e-10f64;
        // szv#4:S4163:12Mar99 warning
        // IntRes2d_Domain domain(Crv->Value(a), a, tolint, Crv->Value(b), b, tolint);
        // Geom2dInt_GInter Inter(AC, domain, tolint, tolint) — the curve
        // self-intersection (bridge #5 GAP: the !IsDone branch stays live).
        let inter = GInterCurveCurve::new();

        if !inter.is_done() {
            return false;
        }

        let v1 = sae.first_vertex(brep, &edge);
        let v2 = sae.last_vertex(brep, &edge);
        if v1.is_null() || v2.is_null() {
            self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            return false;
        }
        let tol1 = brep_tool_tolerance(&v1);
        let tol2 = brep_tool_tolerance(&v2);

        let _pnt1 = brep_tool_pnt(&v1);
        let _pnt2 = brep_tool_pnt(&v2);

        // The OCCT intersection-point walk is unreachable while the GAP keeps
        // IsDone false; the tol1/tol2 comparisons anchor the loop body.
        let _ = (tol1, tol2);

        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckSelfIntersectingEdge(num) (cxx L1346-1352).
    pub fn check_self_intersecting_edge(&mut self, brep: &mut BRep, num: i32) -> bool {
        let mut points2d = Vec::new();
        let mut points3d = Vec::new();
        self.check_self_intersecting_edge_full(brep, num, &mut points2d, &mut points3d)
    }

    /// OCCT CheckIntersectingEdges(num, points2d, points3d, errors)
    /// (cxx L1360-1545).
    pub fn check_intersecting_edges_full(
        &mut self,
        brep: &mut BRep,
        num: i32,
        _points2d: &mut Vec<DVec2>,
        _points3d: &mut Vec<DVec3>,
        _errors: &mut Vec<f64>,
    ) -> bool {
        _points2d.clear();
        _points3d.clear();
        _errors.clear();
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 2 {
            return false;
        }

        // szv#4:S4163:12Mar99 optimized
        let nb_edges = self.nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let edge1 = self.my_wire.as_ref().unwrap().edge(n1);
        let edge2 = self.my_wire.as_ref().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.last_vertex(brep, &edge1);
        let v2 = sae.first_vertex(brep, &edge2);
        if v1.is_null() || v2.is_null() {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        if !brep_tools_compare(brep, &v1, &v2) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            return false;
        }

        let vp = sae.first_vertex(brep, &edge1);
        let vn = sae.last_vertex(brep, &edge2);
        let _ = (vp, vn);

        let (mut a1, mut b1, mut a2, mut b2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut crv1: Option<rcad_kernel::geom::Curve2d> = None;
        let mut crv2: Option<rcad_kernel::geom::Curve2d> = None;
        if !sae.pcurve_face(brep, &edge1, &self.my_face, &mut crv1, &mut a1, &mut b1, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }
        if !sae.pcurve_face(brep, &edge2, &self.my_face, &mut crv2, &mut a2, &mut b2, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }
        if (a1 - b1).abs() <= PCONFUSION || (a2 - b2).abs() <= PCONFUSION {
            return false; //:f7 abv 6 May 98: BUC50070 on #42276
        }

        let is_forward1 = edge1.orientation == Orientation::Forward;
        let is_forward2 = edge2.orientation == Orientation::Forward;
        let _ = (is_forward1, is_forward2);

        let tol0 = brep_tool_tolerance(&v1).max(brep_tool_tolerance(&v2));
        let tol = tol0;

        let _pnt = brep_tool_pnt(&v1);

        let tolint = 1.0e-10f64;

        // Geom2dInt_GInter Inter; Inter.Perform(C1, d1, C2, d2, tolint,
        // tolint) — the curve-curve intersection (bridge #5 GAP: the
        // !IsDone branch stays live).
        let inter = GInterCurveCurve::new();
        let _ = tolint;
        if !inter.is_done() {
            return false;
        }

        //: 86 abv 22 Jan 98
        let tole = (
            if brep_tool_same_parameter(&edge1) {
                brep_tool_tolerance(&edge1)
            } else {
                tol0
            },
            if brep_tool_same_parameter(&edge2) {
                brep_tool_tolerance(&edge2)
            } else {
                tol0
            },
        );
        let tolt = tol.min(tole.0.max(tole.1));
        let _ = tolt;
        //: l0 abv: CATIA01 #1727
        let _is_lacking: i32 = -1;


        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckIntersectingEdges(num) (cxx L1547-1555).
    pub fn check_intersecting_edges(&mut self, brep: &mut BRep, num: i32) -> bool {
        let mut points2d = Vec::new();
        let mut points3d = Vec::new();
        let mut errors = Vec::new();
        self.check_intersecting_edges_full(brep, num, &mut points2d, &mut points3d, &mut errors)
    }

    /// OCCT CheckIntersectingEdges(num1, num2, points2d, points3d, errors)
    /// (cxx L1557-1695).
    pub fn check_intersecting_edges_full_pair(
        &mut self,
        brep: &mut BRep,
        num1: i32,
        num2: i32,
        _points2d: &mut Vec<DVec2>,
        _points3d: &mut Vec<DVec3>,
        _errors: &mut Vec<f64>,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }
        let nb = self.nb_edges();
        let n2 = if num2 > 0 { num2 } else { nb };
        let n1 = if num1 > 0 { num1 } else { nb };

        let edge1 = self.my_wire.as_ref().unwrap().edge(n1);
        let edge2 = self.my_wire.as_ref().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        let (mut a1, mut b1, mut a2, mut b2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut crv1: Option<rcad_kernel::geom::Curve2d> = None;
        let mut crv2: Option<rcad_kernel::geom::Curve2d> = None;
        if !sae.pcurve_face(brep, &edge1, &self.my_face, &mut crv1, &mut a1, &mut b1, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }
        if !sae.pcurve_face(brep, &edge2, &self.my_face, &mut crv2, &mut a2, &mut b2, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }

        if (a1 - b1).abs() <= PCONFUSION || (a2 - b2).abs() <= PCONFUSION {
            return false;
        }

        _points2d.clear();
        _points3d.clear();
        _errors.clear();

        // Geom2dInt_GInter Inter; Inter.Perform(C1, d1, C2, d2, tolint,
        // tolint) — bridge #5 GAP: the !IsDone branch stays live.
        let inter = GInterCurveCurve::new();
        if !inter.is_done() {
            return false;
        }

        // #83 rln 19.03.99 — the point/segment walk is unreachable while the
        // GAP keeps IsDone false.
        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckIntersectingEdges(num1, num2) (cxx L1697-1705).
    pub fn check_intersecting_edges_pair(&mut self, brep: &mut BRep, num1: i32, num2: i32) -> bool {
        let mut points2d = Vec::new();
        let mut points3d = Vec::new();
        let mut errors = Vec::new();
        self.check_intersecting_edges_full_pair(
            brep,
            num1,
            num2,
            &mut points2d,
            &mut points3d,
            &mut errors,
        )
    }

    /// OCCT CheckLacking(num, Tolerance, p2d1, p2d2) (cxx L1711-1795) —
    /// tests if two edges are disconnected in 2d according to the
    /// Adaptor_Surface::Resolution.
    pub fn check_lacking_full(
        &mut self,
        brep: &mut BRep,
        num: i32,
        tolerance: f64,
        p2d1: &mut DVec2,
        p2d2: &mut DVec2,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        // szv#4:S4163:12Mar99 optimized
        let nb_edges = self.nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let e1 = self.my_wire.as_ref().unwrap().edge(n1);
        let e2 = self.my_wire.as_ref().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.last_vertex(brep, &e1);
        let v2 = sae.first_vertex(brep, &e2);
        // CKY 4 MAR 1998 : protection against null vertex
        if v1.is_null() || v2.is_null() {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        if !brep_tools_compare(brep, &v1, &v2) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            return false;
        }

        let (mut a, mut b) = (0.0f64, 0.0f64);
        let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
        if !sae.pcurve_face(brep, &e1, &self.my_face, &mut c2d, &mut a, &mut b, true) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }
        // anAdapt.D1(b, p2d1, v1).
        let (p_1, v_1) = Curve2dAdaptorD1::d1(c2d.as_ref().unwrap(), b);
        *p2d1 = p_1;
        let mut v1v = v_1;
        if e1.orientation == Orientation::Reversed {
            v1v = -v1v;
        }
        if !sae.pcurve_face(brep, &e2, &self.my_face, &mut c2d, &mut a, &mut b, true) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }
        let (p_2, v_2) = Curve2dAdaptorD1::d1(c2d.as_ref().unwrap(), a);
        *p2d2 = p_2;
        let mut v2v = v_2;
        if e2.orientation == Orientation::Reversed {
            v2v = -v2v;
        }
        let v12 = *p2d2 - *p2d1;
        self.my_max2d = v12.length_squared();

        // test like in BRepCheck
        let mut tol = brep_tool_tolerance(&v1).max(brep_tool_tolerance(&v2));
        tol = if tolerance > GP_RESOLUTION && tolerance < tol {
            tolerance
        } else {
            tol
        };
        let surf = self.my_surf.as_ref().unwrap().surface();
        let tol2d = 2.0 * surf.u_resolution(tol).max(surf.v_resolution(tol));
        if self.my_max2d < tol2d * tol2d {
            return false;
        }

        self.my_max2d = self.my_max2d.sqrt();
        self.my_max3d = tol * self.my_max2d / tol2d.max(GP_RESOLUTION);
        self.my_status |= encode_status(ShapeExtendStatus::Done1);

        if self.my_max2d < PCONFUSION
            || (v1v.length_squared() > GP_RESOLUTION
                && signed_angle(v12, v1v).abs() > 0.9 * std::f64::consts::PI)
            || (v2v.length_squared() > GP_RESOLUTION
                && signed_angle(v12, v2v).abs() > 0.9 * std::f64::consts::PI)
        {
            self.my_status |= encode_status(ShapeExtendStatus::Done2);
        }
        true
    }

    /// OCCT CheckLacking(num, Tolerance) (cxx L1797-1803).
    pub fn check_lacking(&mut self, brep: &mut BRep, num: i32, tolerance: f64) -> bool {
        let mut p1 = DVec2::ZERO;
        let mut p2 = DVec2::ZERO;
        self.check_lacking_full(brep, num, tolerance, &mut p1, &mut p2)
    }

    /// OCCT CheckOuterBound(APIMake) (cxx L1805-1830).
    pub fn check_outer_bound(&mut self, brep: &mut BRep, apimake: bool) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        // TopoDS_Wire wire = myWire->Wire()/WireAPIMake() — the rcad WireData
        // wire assembly (the wire as the edge sequence).
        let wire = WireDataWire::wire(self.my_wire.as_ref().unwrap(), brep, apimake);

        // myFace.EmptyCopied() — the face with no wires (the pool face shell).
        let face = rcad_kernel::topo::topods::BRepBuilder::new().make_face(
            brep,
            match self.my_face.data.as_ref() {
                TShape::Face(fd) => fd.surface.clone(),
                _ => None,
            },
            wire.clone(),
        );
        if analysis::is_outer_bound(brep, &face) {
            return false;
        }
        self.my_status = encode_status(ShapeExtendStatus::Done1);
        true
    }

    /// OCCT CheckNotchedEdges(num, shortNum, param, Tolerance)
    /// (cxx L1858-2000).
    pub fn check_notched_edges(
        &mut self,
        brep: &mut BRep,
        num: i32,
        short_num: &mut i32,
        param: &mut f64,
        tolerance: f64,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let nb_edges = self.nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let e1 = self.my_wire.as_ref().unwrap().edge(n1);
        let e2 = self.my_wire.as_ref().unwrap().edge(n2);

        if brep_tool_degenerated(&e1) || brep_tool_degenerated(&e2) {
            return false;
        }

        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.last_vertex(brep, &e1);
        let v2 = sae.first_vertex(brep, &e2);

        if v1.is_null() || v2.is_null() {
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }
        if !brep_tools_compare(brep, &v1, &v2) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail2);
            return false;
        }

        let (mut a1, mut b1, mut a2, mut b2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut c2d1: Option<rcad_kernel::geom::Curve2d> = None;
        let mut c2d2: Option<rcad_kernel::geom::Curve2d> = None;
        if !sae.pcurve_face(brep, &e1, &self.my_face, &mut c2d1, &mut a1, &mut b1, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }

        let (mut p2d1, mut v1v) = (DVec2::ZERO, DVec2::ZERO);
        if e1.orientation == Orientation::Reversed {
            let (p, v) = Curve2dAdaptorD1::d1(c2d1.as_ref().unwrap(), a1);
            p2d1 = p;
            v1v = v;
        } else {
            let (p, v) = Curve2dAdaptorD1::d1(c2d1.as_ref().unwrap(), b1);
            p2d1 = p;
            v1v = -v;
        }

        if !sae.pcurve_face(brep, &e2, &self.my_face, &mut c2d2, &mut a2, &mut b2, false) {
            self.my_status |= encode_status(ShapeExtendStatus::Fail3);
            return false;
        }
        let (mut p2d2, mut v2v) = (DVec2::ZERO, DVec2::ZERO);
        if e2.orientation == Orientation::Reversed {
            let (p, v) = Curve2dAdaptorD1::d1(c2d2.as_ref().unwrap(), b2);
            p2d2 = p;
            v2v = -v;
        } else {
            let (p, v) = Curve2dAdaptorD1::d1(c2d2.as_ref().unwrap(), a2);
            p2d2 = p;
            v2v = v;
        }

        if v2v.length() < GP_RESOLUTION || v1v.length() < GP_RESOLUTION {
            return false;
        }

        // |v2.Angle(v1)| > 0.1 — the oriented 2D angle.
        let angle = signed_angle(v1v, v2v);
        if angle.abs() > 0.1 || p2d1.distance(p2d2) > tolerance {
            return false;
        }

        // The Adaptor3d_CurveOnSurface over the default XY plane (the OCCT
        // Geom_Plane(gp_Pln()) planes) — the lifted Value is the identity
        // embedding of the pcurve point.
        let plane = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
            u_dir: DVec3::X,
            v_dir: DVec3::Y,
        });
        let ad1 = CurveOnPlane {
            c2d: c2d1.clone().unwrap(),
            surface: plane.clone(),
        };
        let ad2 = CurveOnPlane {
            c2d: c2d2.clone().unwrap(),
            surface: plane,
        };

        let sac = ShapeAnalysisCurve;

        let mut proj1 = DVec3::ZERO;
        let mut proj2 = DVec3::ZERO;
        let mut param1 = 0.0f64;
        let mut param2 = 0.0f64;
        let end2 = Curve2dEval::point_at(
            c2d2.as_ref().unwrap(),
            if e2.orientation == Orientation::Forward {
                b2
            } else {
                a2
            },
        );
        let start1 = Curve2dEval::point_at(
            c2d1.as_ref().unwrap(),
            if e1.orientation == Orientation::Forward {
                a1
            } else {
                b1
            },
        );
        let dist1 = project_inside(
            &ad1,
            DVec3::new(end2.x, end2.y, 0.0),
            tolerance,
            &mut proj1,
            &mut param1,
            false,
        );
        let dist2 = project_inside(
            &ad2,
            DVec3::new(start1.x, start1.y, 0.0),
            tolerance,
            &mut proj2,
            &mut param2,
            false,
        );

        if dist1 > tolerance && dist2 > tolerance {
            return false;
        }

        let (short_ad, long_ad, first_p) = if dist1 < dist2 {
            *short_num = n2;
            *param = param1;
            (ad2, ad1, a2)
        } else {
            *short_num = n1;
            *param = param2;
            (ad1, ad2, a1)
        };
        let len_p = if dist1 < dist2 { b2 - a2 } else { b1 - a1 };

        let step = len_p / 23.0;
        for i in 1..23 {
            // OCCT: for (i = 1; i < 23; i++, firstP += step).
            let fp = first_p + step * i as f64;
            let val = CurveOnPlane::value(&short_ad, fp);
            let mut prj = DVec3::ZERO;
            let mut prm = 0.0f64;
            let d1 = project_inside(&long_ad, val, tolerance, &mut prj, &mut prm, true);
            if d1 > tolerance {
                return false;
            }
        }

        true
    }

    /// OCCT CheckSmallArea(theWire) (cxx L2006-2115).
    pub fn check_small_area(&mut self, brep: &mut BRep, the_wire: &Shape) -> bool {
        let _ = the_wire;
        self.my_status = encode_status(ShapeExtendStatus::Fail1);
        let a_nb_control = 23;
        let nb_edges = self.nb_edges();
        if !self.is_ready() || nb_edges < 1 {
            return false;
        }
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        // The 2D center accumulation walk (the sampling loop).
        let an_inv = 1.0 / (a_nb_control - 1) as f64;
        let mut a_center2d = DVec2::ZERO;
        let sae = ShapeAnalysisEdge::new();
        for j in 1..=nb_edges {
            let e = self.my_wire.as_ref().unwrap().edge(j);
            let mut a_curve2d: Option<rcad_kernel::geom::Curve2d> = None;
            let (mut af, mut al) = (0.0f64, 0.0f64);
            if !sae.pcurve_face(brep, &e, &self.my_face, &mut a_curve2d, &mut af, &mut al, false) {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
                return false;
            }
            let c2d = a_curve2d.as_ref().unwrap();
            for i in 1..a_nb_control {
                let a_v = an_inv * (((a_nb_control - 1 - i) as f64 * af + i as f64 * al));
                a_center2d += Curve2dEval::point_at(c2d, a_v);
            }
        }
        a_center2d /= (nb_edges * (a_nb_control - 1)) as f64;

        // check approximated area in 3D
        let a_center3d = self.my_surf.as_ref().unwrap().value(a_center2d);
        let mut a_prev3d = DVec3::ZERO;
        let mut a_cross = DVec3::ZERO;
        let mut a_length = 0.0f64;
        for j in 1..=nb_edges {
            let e = self.my_wire.as_ref().unwrap().edge(j);
            let mut a_curve3d: Option<rcad_kernel::geom::Curve3> = None;
            let (mut af, mut al) = (0.0f64, 0.0f64);
            if !sae.curve3d(brep, &e, &mut a_curve3d, &mut af, &mut al, false) {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
                return false;
            }
            let c3d = a_curve3d.as_ref().unwrap();
            if rcad_kernel::precision::is_infinite_value(af)
                || rcad_kernel::precision::is_infinite_value(al)
            {
                continue;
            }

            let mut a_begin = 0usize;
            let mut a_pnt3d = DVec3::ZERO;
            if j == 1 {
                a_begin = 1;
                a_pnt3d = CurveEval::point_at(c3d, af);
                a_prev3d = a_pnt3d - a_center3d;
            }
            for i in a_begin..a_nb_control as usize {
                let an_u = an_inv * (((a_nb_control - 1) as usize - i) as f64 * af + i as f64 * al);
                let a_pnt = CurveEval::point_at(c3d, an_u);
                let a_vec = a_pnt - a_center3d;

                a_cross += a_prev3d.cross(a_vec);
                a_length += a_pnt3d.distance(a_pnt);

                a_pnt3d = a_pnt;
                a_prev3d = a_vec;
            }
        }

        let a_tolerance = a_length * self.my_precision;
        if a_cross.length() < a_tolerance {
            // check real area in 3D — GAP (bridge #8): BRepGProp
            // SurfaceProperties/LinearProperties are untranslated; the OCCT
            // failure path (the check reports not-small) is preserved.
            return false;
        }

        false
    }

    /// OCCT CheckShapeConnect(tailhead, tailtail, headtail, headhead, shape,
    /// prec) (cxx L2130-2266).
    pub fn check_shape_connect_full(
        &mut self,
        brep: &mut BRep,
        tailhead: &mut f64,
        tailtail: &mut f64,
        headtail: &mut f64,
        headhead: &mut f64,
        shape: &Shape,
        prec: f64,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Fail1);
        if !self.is_loaded() || shape.is_null() {
            return false;
        }
        let sae = ShapeAnalysisEdge::new();
        let mut v1: Option<Shape> = None;
        let mut v2: Option<Shape> = None;
        if shape.shape_type() == ShapeType::Edge {
            v1 = Some(sae.first_vertex(brep, shape));
            v2 = Some(sae.last_vertex(brep, shape));
        } else if shape.shape_type() == ShapeType::Wire {
            analysis::find_bounds(brep, shape, &mut v1, &mut v2);
        } else {
            return false;
        }
        if v1.is_none() || v2.is_none() {
            return false;
        }
        let v1 = v1.unwrap();
        let v2 = v2.unwrap();
        self.my_status = encode_status(ShapeExtendStatus::Done1);
        //  on va comparer les points avec ceux de thevfirst et thevlast
        let p1 = brep_tool_pnt(&v1);
        let p2 = brep_tool_pnt(&v2);

        let vfirst = sae.first_vertex(brep, &self.my_wire.as_ref().unwrap().edge(1));
        let vlast = sae.last_vertex(
            brep,
            &self
                .my_wire
                .as_ref()
                .unwrap()
                .edge(self.nb_edges()),
        );
        let pf = brep_tool_pnt(&vfirst);
        let pl = brep_tool_pnt(&vlast);

        *tailhead = p1.distance(pl);
        *tailtail = p2.distance(pl);
        *headtail = p1.distance(pf);
        *headhead = p2.distance(pf);
        let mut dm1 = *tailhead;
        let mut dm2 = *headtail;
        let mut res1 = 0;
        let mut res2 = 0;

        if *tailhead > *tailtail {
            res1 = 1;
            dm1 = *tailtail;
        }
        if *headtail > *headhead {
            res2 = 1;
            dm2 = *headhead;
        }
        let mut result = res1;
        self.my_min3d = dm1.min(dm2);
        self.my_max3d = dm1.max(dm2);
        if dm1 > dm2 {
            dm1 = dm2;
            result = res2 + 2;
        }
        match result {
            1 => self.my_status = encode_status(ShapeExtendStatus::Done2),
            2 => self.my_status = encode_status(ShapeExtendStatus::Done3),
            3 => self.my_status = encode_status(ShapeExtendStatus::Done4),
            _ => {}
        }
        if res1 == 0 {
            self.my_status |= encode_status(ShapeExtendStatus::Done5);
        }
        if res2 == 0 {
            self.my_status |= encode_status(ShapeExtendStatus::Done6);
        }

        if self.my_min3d > self.my_precision.max(prec) {
            self.my_status = encode_status(ShapeExtendStatus::Fail2);
        }
        self.last_check_status(ShapeExtendStatus::Done)
    }

    /// OCCT CheckShapeConnect(shape, prec) (cxx L2117-2125).
    pub fn check_shape_connect(&mut self, brep: &mut BRep, shape: &Shape, prec: f64) -> bool {
        let mut tailhead = 0.0f64;
        let mut tailtail = 0.0f64;
        let mut headtail = 0.0f64;
        let mut headhead = 0.0f64;
        self.check_shape_connect_full(
            brep,
            &mut tailhead,
            &mut tailtail,
            &mut headtail,
            &mut headhead,
            shape,
            prec,
        )
    }

    /// OCCT CheckLoop(aMapLoopVertices, aMapVertexEdges, aMapSmallEdges,
    /// aMapSeemEdges) (cxx L2276-2388) — the vertex-edge maps keyed by the
    /// vertex identity (ptr, location).  The out maps carry the OCCT NCollection
    /// containers: the loop vertices the IndexedMap (ordered, unique), the
    /// vertex-edge lists the DataMap, the small/seam edges the Maps (unique).
    pub fn check_loop(
        &mut self,
        brep: &mut BRep,
        a_map_loop_vertices: &mut Vec<(u64, u32)>,
        a_map_vertex_edges: &mut std::collections::HashMap<(u64, u32), Vec<Shape>>,
        a_map_small_edges: &mut Vec<(u64, u32)>,
        a_map_seem_edges: &mut Vec<(u64, u32)>,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() < 2 {
            return false;
        }
        let a_sav_preci = self.precision();
        self.set_precision(rcad_kernel::precision::INFINITE_VALUE);

        for i in 1..=self.nb_edges() {
            let aedge = self.my_wire.as_ref().unwrap().edge(i);
            // TopExp::Vertices(aedge, aV1, aV2).
            let (a_v1, a_v2) = edge_vertices(&aedge);
            if a_v1.is_null() || a_v2.is_null() {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
                return false;
            }
            let is_same = a_v1.ptr_id() == a_v2.ptr_id() && a_v1.location == a_v2.location;
            if self.my_wire.as_mut().unwrap().is_seam(i) {
                let k = (aedge.ptr_id(), aedge.location);
                if !a_map_seem_edges.contains(&k) {
                    a_map_seem_edges.push(k);
                }
            } else if brep_tool_degenerated(&aedge) {
                let k = (aedge.ptr_id(), aedge.location);
                if !a_map_small_edges.contains(&k) {
                    a_map_small_edges.push(k);
                }
            } else if is_same && self.check_small(brep, i, brep_tool_tolerance(&a_v1)) {
                let k = (aedge.ptr_id(), aedge.location);
                if !a_map_small_edges.contains(&k) {
                    a_map_small_edges.push(k);
                }
            }

            let k1 = (a_v1.ptr_id(), a_v1.location);
            let k2 = (a_v2.ptr_id(), a_v2.location);
            a_map_vertex_edges.entry(k1).or_default();
            a_map_vertex_edges.entry(k2).or_default();
            if is_same {
                let alshape = a_map_vertex_edges.get_mut(&k1).unwrap();
                alshape.push(aedge.clone());
                alshape.push(aedge.clone());
                if alshape.len() > 2
                    && is_multi_vertex(alshape, a_map_small_edges, a_map_seem_edges)
                {
                    if !a_map_loop_vertices.contains(&k1) {
                        a_map_loop_vertices.push(k1);
                    }
                }
            } else {
                let alshape = a_map_vertex_edges.get_mut(&k1).unwrap();
                alshape.push(aedge.clone());
                if alshape.len() > 2
                    && is_multi_vertex(alshape, &a_map_small_edges, &a_map_seem_edges)
                {
                    if !a_map_loop_vertices.contains(&k1) {
                        a_map_loop_vertices.push(k1);
                    }
                }
                let alshape2 = a_map_vertex_edges.get_mut(&k2).unwrap();
                alshape2.push(aedge.clone());
                if alshape2.len() > 2
                    && is_multi_vertex(alshape2, a_map_small_edges, a_map_seem_edges)
                {
                    if !a_map_loop_vertices.contains(&k2) {
                        a_map_loop_vertices.push(k2);
                    }
                }
            }
        }
        self.set_precision(a_sav_preci);
        if !a_map_loop_vertices.is_empty() {
            self.my_status = encode_status(ShapeExtendStatus::Done1);
            self.my_status_loop |= self.my_status;
            return true;
        }
        false
    }

    /// OCCT CheckTail(theEdge1, theEdge2, theMaxSine, theMaxWidth,
    /// theMaxTolerance, theEdge11, theEdge12, theEdge21, theEdge22)
    /// (cxx L2346-2599).
    #[allow(clippy::too_many_arguments)]
    pub fn check_tail(
        &mut self,
        brep: &mut BRep,
        the_edge1: &Shape,
        the_edge2: &Shape,
        the_max_sine: f64,
        the_max_width: f64,
        the_max_tolerance: f64,
        the_edge11: &mut Option<Shape>,
        the_edge12: &mut Option<Shape>,
        the_edge21: &mut Option<Shape>,
        the_edge22: &mut Option<Shape>,
    ) -> bool {
        let a_es = [the_edge1.clone(), the_edge2.clone()];
        if !self.is_ready() || brep_tool_degenerated(&a_es[0]) || brep_tool_degenerated(&a_es[1]) {
            return false;
        }

        // Check the distance between the edge common ends.
        let a_tol2 = the_max_width + 0.5 * CONFUSION_C;
        let a_tol3 = the_max_width + CONFUSION_C;
        let _a_tol4 = the_max_width + 1.5 * CONFUSION_C;
        let a_sq_tol2 = a_tol2 * a_tol2;
        let _a_sq_tol3 = a_tol3 * a_tol3;
        let mut a_cs: [Option<rcad_kernel::geom::Curve3>; 2] = [None, None];
        let mut a_ls = [[0.0f64; 2]; 2];
        let mut a_vis = [0usize; 2];
        let mut a_vps = [DVec3::ZERO; 2];
        {
            let sae = ShapeAnalysisEdge::new();
            for (a_ei, a_e) in a_es.iter().enumerate() {
                let (mut f0, mut l0) = (0.0f64, 0.0f64);
                let mut c0: Option<rcad_kernel::geom::Curve3> = None;
                if !sae.curve3d(brep, a_e, &mut c0, &mut f0, &mut l0, false) {
                    return false;
                }
                a_cs[a_ei] = c0;
                a_ls[a_ei][0] = f0;
                a_ls[a_ei][1] = l0;

                a_vis[a_ei] = if a_e.orientation == Orientation::Reversed {
                    a_ei
                } else {
                    1 - a_ei
                };
                a_vps[a_ei] =
                    CurveEval::point_at(a_cs[a_ei].as_ref().unwrap(), a_ls[a_ei][a_vis[a_ei]]);
            }
            if a_vps[0].distance_squared(a_vps[1]) > a_sq_tol2 {
                return false;
            }
        }

        // Check the angle between the edges.
        if the_max_sine >= 0.0 {
            let a_sq_max_sine = the_max_sine * the_max_sine;
            let mut a_ds = [DVec3::ZERO; 2];
            let mut a_reverse = false;
            for (a_ei, a_e) in a_es.iter().enumerate() {
                let c = a_cs[a_ei].as_ref().unwrap();
                // GCPnts_AbscissaPoint::Length(aCA, f, l, 0.25*Confusion)
                // — the kernel arc-length walk.
                let len = rcad_kernel::base::gcpnts::abscissa_point::arc_length(c, a_ls[a_ei][0], a_ls[a_ei][1]);
                if len < 0.5 * CONFUSION_C {
                    return false;
                }

                // GCPnts_AbscissaPoint aAP(0.25*Confusion, aCA,
                // 0.5*Confusion*(1 - 2*aVIs[aEI]), aLs[aEI][aVIs[aEI]]).
                let a_param =
                    rcad_kernel::base::gcpnts::abscissa_point::abscissa_point_parameter(
                        c,
                        a_ls[a_ei][0],
                        a_ls[a_ei][1],
                        0.5 * CONFUSION_C * (1.0 - 2.0 * a_vis[a_ei] as f64),
                        a_ls[a_ei][a_vis[a_ei]],
                    );
                // OCCT: if (!aAP.IsDone()) return false — the kernel form
                // signals failure by a NaN parameter.
                if a_param.is_nan() {
                    return false;
                }

                let mut a_ps = [DVec3::ZERO; 2];
                a_ps[a_vis[a_ei]] = a_vps[a_ei];
                a_ps[1 - a_vis[a_ei]] = CurveEval::point_at(c, a_param);
                a_ds[a_ei] = a_ps[1] - a_ps[0];
                let a_dn = a_ds[a_ei].length();
                if a_dn < 0.1 * CONFUSION_C {
                    return false;
                }

                a_ds[a_ei] *= 1.0 / a_dn;
                a_reverse ^= a_vis[a_ei] == 1;
            }
            if a_reverse {
                a_ds[0] = -a_ds[0];
            }
            if a_ds[0].dot(a_ds[1]) < 0.0
                || a_ds[0].cross(a_ds[1]).length_squared() > a_sq_max_sine
            {
                return false;
            }
        }

        // The tail-bounds walk and the final cut consume the
        // ShapeAnalysis_TransferParametersProj / GProp machinery (bridge #8
        // GAP): the OCCT failure path (no tail reported) is preserved.
        false
    }
}

// ---------------------------------------------------------------------------
// Statics and helpers (the cxx anonymous-namespace helpers)
// ---------------------------------------------------------------------------

// OCCT Standard_Real.hxx: the Confusion constant of the cxx literals.
const CONFUSION_C: f64 = 1.0e-7;

/// OCCT GetPointOnEdge(edge, surf, Crv2d, param) (cxx L1248-1266) — the
/// point is taken from the 3d curve when the edge is SameParameter.
fn get_point_on_edge(
    edge: &Shape,
    surf: &super::surface::ShapeAnalysisSurface,
    crv2d: &rcad_kernel::geom::Curve2d,
    param: f64,
) -> DVec3 {
    if brep_tool_same_parameter(edge) {
        // BRep_Tool::Curve(edge, L, f, l) — the world curve.
        if let TShape::Edge(ed) = edge.data.as_ref() {
            if let Some(con_s) = ed.curve.as_ref() {
                return CurveEval::point_at(con_s, param);
            }
        }
    }
    let a_p2d = Curve2dEval::point_at(crv2d, param);
    surf.value(a_p2d)
}

/// OCCT ProjectInside(AD, pnt, preci, proj, param, adjustToEnds = true)
/// (cxx L1837-1856).
fn project_inside(
    ad: &CurveOnPlane,
    pnt: DVec3,
    preci: f64,
    proj: &mut DVec3,
    param: &mut f64,
    adjust_to_ends: bool,
) -> f64 {
    let sac = ShapeAnalysisCurve;
    let adaptor = super::curve::Adaptor3dCurve::value(ad, 0.0); // placeholder, unused
    let _ = adaptor;
    let mut dist = {
        // sac.Project(AD, pnt, preci, proj, param, adjustToEnds) — the
        // Adaptor3d_Curve overload through the CurveOnPlane adaptor.
        let mut p = DVec3::ZERO;
        let mut prm = 0.0f64;
        let d = sac.project_adaptor(ad, pnt, preci, &mut p, &mut prm, adjust_to_ends);
        *proj = p;
        *param = prm;
        d
    };
    let u_first = super::curve::Adaptor3dCurve::first_parameter(ad);
    let u_last = super::curve::Adaptor3dCurve::last_parameter(ad);
    if *param < u_first {
        *param = u_first;
        *proj = super::curve::Adaptor3dCurve::value(ad, u_first);
        return proj.distance(pnt);
    }
    if *param > u_last {
        *param = u_last;
        *proj = super::curve::Adaptor3dCurve::value(ad, u_last);
        return proj.distance(pnt);
    }
    let _ = &mut dist;
    dist
}

/// OCCT Project(theCurve, theFirstParameter, theLastParameter, thePoint,
/// thePrecision, theParameter, theProjection) (cxx L2316-2344).
#[allow(clippy::too_many_arguments)]
fn project(
    the_curve: &rcad_kernel::geom::Curve3,
    the_first_parameter: f64,
    the_last_parameter: f64,
    the_point: DVec3,
    the_precision: f64,
    the_parameter: &mut f64,
    the_projection: &mut DVec3,
) -> f64 {
    let sac = ShapeAnalysisCurve;
    let a_dist = sac.project_cf_cl(
        the_curve,
        the_point,
        the_precision,
        the_projection,
        the_parameter,
        the_first_parameter,
        the_last_parameter,
        true,
    );
    if *the_parameter >= the_first_parameter && *the_parameter <= the_last_parameter {
        return a_dist;
    }

    let a_params = [the_first_parameter, the_last_parameter];
    let a_prjs = [
        CurveEval::point_at(the_curve, a_params[0]),
        CurveEval::point_at(the_curve, a_params[1]),
    ];
    let a_dists = [
        the_point.distance(a_prjs[0]),
        the_point.distance(a_prjs[1]),
    ];
    let a_pi = if a_dists[0] <= a_dists[1] { 0usize } else { 1usize };
    *the_parameter = a_params[a_pi];
    *the_projection = a_prjs[a_pi];
    a_dists[a_pi]
}

/// OCCT isMultiVertex(alshape, aMapSmallEdges, aMapSeemEdges)
/// (cxx L2268-2288).
fn is_multi_vertex(
    alshape: &[Shape],
    a_map_small_edges: &[(u64, u32)],
    a_map_seem_edges: &[(u64, u32)],
) -> bool {
    let mut nb_not_account = 0;
    for s in alshape {
        let key = (s.ptr_id(), s.location);
        if a_map_small_edges.contains(&key) {
            nb_not_account += 1;
        } else if a_map_seem_edges.contains(&key) {
            nb_not_account += 1;
        }
    }
    alshape.len() - nb_not_account > 2
}

/// OCCT TopExp::Vertices(edge, V1, V2) — the orientation-composed end
/// vertices.
fn edge_vertices(edge: &Shape) -> (Shape, Shape) {
    let ed = match edge.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return (Shape::null(), Shape::null()),
    };
    let first = if edge.orientation.compose(ed.first.orientation) == Orientation::Reversed {
        ed.last.clone()
    } else {
        ed.first.clone()
    };
    let last = if edge.orientation.compose(ed.last.orientation) == Orientation::Reversed {
        ed.first.clone()
    } else {
        ed.last.clone()
    };
    (first, last)
}

/// The oriented 2D angle (gp_Vec2d::Angle).
fn signed_angle(a: DVec2, b: DVec2) -> f64 {
    a.y.atan2(a.x) - b.y.atan2(b.x)
}

/// The pcurve D1 (the Geom2dAdaptor_Curve::D1 re-host over the value).
struct Curve2dAdaptorD1;

impl Curve2dAdaptorD1 {
    fn d1(c: &rcad_kernel::geom::Curve2d, u: f64) -> (DVec2, DVec2) {
        use rcad_kernel::geom::Curve2dEval;
        let p = Curve2dEval::point_at(c, u);
        let d = Curve2dEval::derivative_at(c, u);
        (p, d)
    }
}

/// OCCT Adaptor3d_CurveOnSurface over a pcurve and the default plane (the
/// CheckNotchedEdges Geom_Plane(gp_Pln()) construction).
struct CurveOnPlane {
    c2d: rcad_kernel::geom::Curve2d,
    surface: rcad_kernel::geom::Surface3,
}

impl CurveOnPlane {
    fn value(&self, u: f64) -> DVec3 {
        let p = Curve2dEval::point_at(&self.c2d, u);
        rcad_kernel::geom::SurfaceEval::point_at(&self.surface, p.x, p.y)
    }
}

impl CurveEval for CurveOnPlane {
    fn point_at(&self, t: f64) -> DVec3 {
        let p = Curve2dEval::point_at(&self.c2d, t);
        rcad_kernel::geom::SurfaceEval::point_at(&self.surface, p.x, p.y)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        // The plane embedding derivative: d.x * Udir + d.y * Vdir (the
        // OCCT plane iso embedding of the pcurve).
        let d = Curve2dEval::derivative_at(&self.c2d, t);
        if let rcad_kernel::geom::Surface3::Plane(pl) = &self.surface {
            (d.x * pl.u_dir + d.y * pl.v_dir).normalize_or_zero()
        } else {
            DVec3::ZERO
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        Curve2dEval::default_domain(&self.c2d)
    }
}

impl super::curve::Adaptor3dCurve for CurveOnPlane {
    fn first_parameter(&self) -> f64 {
        Curve2dEval::default_domain(&self.c2d)[0]
    }
    fn last_parameter(&self) -> f64 {
        Curve2dEval::default_domain(&self.c2d)[1]
    }
    fn value(&self, the_u: f64) -> DVec3 {
        CurveOnPlane::value(self, the_u)
    }
    fn resolution(&self, the_r3d: f64) -> f64 {
        let ru = self.surface.u_resolution(the_r3d);
        let rv = self.surface.v_resolution(the_r3d);
        // The default arm of Geom2dAdaptor_Curve::Resolution
        // (Precision::Parametric(Ruv)).
        ru.min(rv) * 0.01
    }
    fn is_closed(&self) -> bool {
        Curve2dEval::is_closed(&self.c2d)
    }
    fn is_kind_bounded(&self) -> bool {
        true
    }
    fn get_type(&self) -> super::curve::AdaptorCurveKind {
        super::curve::AdaptorCurveKind::Other
    }
    fn curve3d(&self) -> Option<rcad_kernel::geom::Curve3> {
        None
    }
}

// The State import feeds the CheckOuterBound path (TopAbs_OUT comparison
// lives in analysis::is_outer_bound; kept here for the module contract).

// The WireData wire assembly used by CheckOuterBound (the OCCT
// Wire()/WireAPIMake()).
impl WireDataWire for crate::shhealing::shape_extend::wire_data::WireData {
    fn wire(&self, brep: &mut BRep, apimake: bool) -> Shape {
        let _ = apimake;
        // BRepBuilder::make_wire over the edge sequence.
        let edges: Vec<Shape> = (1..=self.nb_edges()).map(|i| self.edge(i)).collect();
        brep.add_twire(edges)
    }
}

trait WireDataWire {
    fn wire(&self, brep: &mut BRep, apimake: bool) -> Shape;
}
