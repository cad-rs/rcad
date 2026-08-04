// OCCT IntPatch_ImpPrmIntersection.cxx 1:1 Rust translation — intersection
// between a natural quadric patch (Plane/Cone/Cylinder/Sphere) and a
// bi-parametrised surface.
//
// The algorithm combines the boundary-search (IntPatch_TheSOnBounds), the
// interior-point search (IntPatch_TheSearchInside) and the implicit-surface
// walking (IntPatch_TheIWalking), then converts the walking polylines into
// WLine objects and the SOnBounds segments into RLine objects.
//
// OCCT IntPatch_ImpPrmIntersection.cxx L181-1964 (Perform + ComputeTangency /
// Recadre / GetLocalStep), L1968-3892 (seam/pole decomposition helpers).

use std::f64::consts::TAU;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Line2d, Surface3, SurfaceEval};

use crate::geomalgo::int_surf::quadric::{Quadric, QuadricType};
use crate::geomalgo::int_surf::PntOn2S;
use crate::geomalgo::int_patch::so_on_bounds::{ArcFunction, Domain, SOnBounds};
use crate::geomalgo::int_patch::transitions::{make_transition, Transition, TypeTrans};
use crate::geomalgo::int_patch::{IntPatchIType, IntPatchLine, IntPatchVertex};
use crate::geomalgo::top_trans::CurveTransition;
use rcad_kernel::topods::State;

use super::i_walking::IWalking;
use super::path_point::{InteriorPoint, PathPoint};
use super::search_inside::SearchInside;
use super::surf_function::SurfFunction;

/// OCCT IntSurf_QuadricTool::Tolerance (IntSurf_QuadricTool.cxx L17-28).
fn quadric_tool_tolerance(q: &Quadric) -> f64 {
    match q.type_quadric() {
        QuadricType::Sphere => 2e-6 * q.radius(),
        QuadricType::Cylinder => 2e-6 * q.radius(),
        _ => 1e-6,
    }
}

/// OCCT IsSeamOrPole (IntPatch_ImpPrmIntersection.cxx L85-177).
#[allow(clippy::too_many_arguments)]
fn is_seam_or_pole(
    q_surf: &Surface3,
    line: &[PntOn2S],
    is_reversed: bool,
    ref_index: usize,
    tol_3d: f64,
    delta_max: f64,
) -> i32 {
    // OCCT IntPatch_SpecPntType: 0 = None, 1 = Pole, 2 = SeamU, 3 = SeamV,
    // 4 = SeamUV, 5 = PoleSeamU.
    const SPNT_NONE: i32 = 0;
    const SPNT_POLE: i32 = 1;
    const SPNT_SEAM_U: i32 = 2;
    const SPNT_SEAM_V: i32 = 3;
    const SPNT_SEAM_UV: i32 = 4;
    const SPNT_POLE_SEAM_U: i32 = 5;

    // OCCT: (theRefIndex < 1) || (theRefIndex >= NbPoints) — with 1-based
    // NbPoints, theRefIndex+1 stays valid, so in 0-based: ref_index+1 < len.
    if ref_index == 0 || ref_index + 1 >= line.len() {
        return SPNT_NONE;
    }

    let a_uq_ref;
    let a_vq_ref;
    let a_up_ref;
    let a_vp_ref;
    let a_uq_next;
    let a_vq_next;
    let a_up_next;
    let a_vp_next;

    let a_p3d = line[ref_index].value();

    if is_reversed {
        let (u1, v1, u2, v2) = line[ref_index].parameters();
        (a_up_ref, a_vp_ref, a_uq_ref, a_vq_ref) = (u1, v1, u2, v2);
        let (u1, v1, u2, v2) = line[ref_index + 1].parameters();
        (a_up_next, a_vp_next, a_uq_next, a_vq_next) = (u1, v1, u2, v2);
    } else {
        let (u1, v1, u2, v2) = line[ref_index].parameters();
        (a_uq_ref, a_vq_ref, a_up_ref, a_vp_ref) = (u1, v1, u2, v2);
        let (u1, v1, u2, v2) = line[ref_index + 1].parameters();
        (a_uq_next, a_vq_next, a_up_next, a_vp_next) = (u1, v1, u2, v2);
    }

    let a_type = q_surf;

    if let Surface3::Cone(c) = a_type {
        let apex = c.apex;
        if apex.distance_squared(a_p3d) < tol_3d * tol_3d {
            return SPNT_POLE_SEAM_U;
        }
    } else if let Surface3::Sphere(s) = a_type {
        let sq_tol = tol_3d * tol_3d;
        let a_p = s.point_at(0.0, std::f64::consts::FRAC_PI_2);
        if a_p.distance_squared(a_p3d) < sq_tol {
            return SPNT_POLE_SEAM_U;
        }
        let a_p = s.point_at(0.0, -std::f64::consts::FRAC_PI_2);
        if a_p.distance_squared(a_p3d) < sq_tol {
            return SPNT_POLE_SEAM_U;
        }
    }

    let a_delta_u = (a_uq_ref - a_uq_next).abs();

    if !matches!(a_type, Surface3::Torus(_)) && a_delta_u < delta_max {
        return SPNT_NONE;
    }

    match a_type {
        Surface3::Cylinder(_) => SPNT_SEAM_U,
        Surface3::Torus(t) => {
            let a_delta_v = (a_vq_ref - a_vq_next).abs();
            if a_delta_u >= delta_max && a_delta_v >= delta_max {
                SPNT_SEAM_UV
            } else if a_delta_u >= delta_max {
                SPNT_SEAM_U
            } else if a_delta_v >= delta_max {
                SPNT_SEAM_V
            } else {
                SPNT_NONE
            }
        }
        Surface3::Sphere(_) | Surface3::Cone(_) => SPNT_POLE_SEAM_U,
        _ => SPNT_NONE,
    }
}

/// OCCT Recadre (IntPatch_ImpPrmIntersection.cxx L473-541).
#[allow(clippy::too_many_arguments)]
fn recadre_impprm(
    type_s1: GeomAbsSurfaceTypeAlias,
    type_s2: GeomAbsSurfaceTypeAlias,
    pt: &mut IntPatchVertex,
    iwline: &[PntOn2S],
    param: usize,
    mut u1: f64,
    mut v1: f64,
    mut u2: f64,
    mut v2: f64,
) {
    // `iwline` is the 0-based point list of the walking line; OCCT reads
    // Value(Param) with 1-based indexing, so the callers pass param 0 for the
    // first vertex and Nbpts-2 (0-based) for the last (OCCT Param = Nbpts-1).
    let (u1p, v1p, u2p, v2p) = iwline[param].parameters();
    let half3 = 1.5 * std::f64::consts::PI;
    match type_s1 {
        GeomAbsSurfaceTypeAlias::Torus => {
            while v1 < v1p - half3 {
                v1 += TAU;
            }
            while v1 > v1p + half3 {
                v1 -= TAU;
            }
            // fallthrough
            while u1 < u1p - half3 {
                u1 += TAU;
            }
            while u1 > u1p + half3 {
                u1 -= TAU;
            }
        }
        GeomAbsSurfaceTypeAlias::Cylinder
        | GeomAbsSurfaceTypeAlias::Cone
        | GeomAbsSurfaceTypeAlias::Sphere => {
            while u1 < u1p - half3 {
                u1 += TAU;
            }
            while u1 > u1p + half3 {
                u1 -= TAU;
            }
        }
        _ => {}
    }
    match type_s2 {
        GeomAbsSurfaceTypeAlias::Torus => {
            while v2 < v2p - half3 {
                v2 += TAU;
            }
            while v2 > v2p + half3 {
                v2 -= TAU;
            }
            // fallthrough
            while u2 < u2p - half3 {
                u2 += TAU;
            }
            while u2 > u2p + half3 {
                u2 -= TAU;
            }
        }
        GeomAbsSurfaceTypeAlias::Cylinder
        | GeomAbsSurfaceTypeAlias::Cone
        | GeomAbsSurfaceTypeAlias::Sphere => {
            while u2 < u2p - half3 {
                u2 += TAU;
            }
            while u2 > u2p + half3 {
                u2 -= TAU;
            }
        }
        _ => {}
    }
    pt.set_parameters(u1, v1, u2, v2);
}

/// OCCT GeomAbs_SurfaceType — the surface-kind used by the ImpPrm Recadre.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GeomAbsSurfaceTypeAlias {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Other,
}

fn classify_alias(s: &Surface3) -> GeomAbsSurfaceTypeAlias {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceTypeAlias::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceTypeAlias::Cylinder,
        Surface3::Cone(_) => GeomAbsSurfaceTypeAlias::Cone,
        Surface3::Sphere(_) => GeomAbsSurfaceTypeAlias::Sphere,
        Surface3::Torus(_) => GeomAbsSurfaceTypeAlias::Torus,
        _ => GeomAbsSurfaceTypeAlias::Other,
    }
}

/// OCCT GetLocalStep (IntPatch_ImpPrmIntersection.cxx L545-613).
///
/// rcad data-model note: the OCCT UContinuity/VContinuity and C1-interval
/// branches (L567-609) need Adaptor3d_Surface continuity tables which the rcad
/// Surface3 model does not expose.  The analytic and LinearExtrusion/Revolution
/// surfaces exercised by the boolean pipeline are all C-infinity, so those
/// branches never fire for them.  The Bezier/BSpline branch (L552-565, fired
/// when both continuities > C0) is translated with the available degree and
/// UResolution/VResolution API.
fn get_local_step(surf: &Surface3, step: f64) -> f64 {
    let mut a_local_step = step;
    if let Surface3::BSpline(b) = surf {
        let a_min_res = u_res(surf, rcad_kernel::precision::CONFUSION)
            .min(v_res(surf, rcad_kernel::precision::CONFUSION));
        let a_max_deg = b.degree_u.max(b.degree_v);
        if a_min_res < 1e-10 && a_max_deg > 3 {
            a_local_step = 0.0001;
        }
    } else if let Surface3::Bezier(b) = surf {
        let a_min_res = u_res(surf, rcad_kernel::precision::CONFUSION)
            .min(v_res(surf, rcad_kernel::precision::CONFUSION));
        let a_max_deg = b.control_points.len().saturating_sub(1).max(
            b.control_points.first().map_or(0, |row| row.len().saturating_sub(1)),
        ) as f64;
        if a_min_res < 1e-10 && a_max_deg > 3.0 {
            a_local_step = 0.0001;
        }
    }
    a_local_step.min(step)
}

/// OCCT IntTools_TopolTool::ComputeSamplePoints (IntTools_TopolTool.cxx L95-394):
/// the number of U and V grid sample points of the face's TopolTool.
///
/// rcad note: the analytic quadric branches (Cylinder/Cone/Sphere/Torus) are
/// unreachable in the ImpPrm SearchInside sampling — the sampled surface is
/// always the parametric one (LinearExtrusion/Revolution/BSpline/Bezier/Other)
/// — so they are not translated here.
fn topol_tool_nb_samples(surf: &Surface3, dom: [f64; 4]) -> (usize, usize) {
    const A_MAX_NB_SAMPLE: i32 = 50;
    let (_uinf, _usup, vinf, vsup) = (dom[0], dom[1], dom[2], dom[3]);
    let mut nbsu = 0i32;
    let mut nbsv = 0i32;
    match surf {
        Surface3::Plane(_) => {
            // OCCT GeomAbs_Plane (L140-144).
            nbsu = 10;
            nbsv = 10;
        }
        Surface3::LinearExtrusion(_) => {
            // OCCT GeomAbs_SurfaceOfExtrusion (L361-373): U (profile) is fixed
            // at 15; V (extrusion) grows with the parameter span.
            nbsu = 15;
            nbsv = ((vsup - vinf) as i32) / 10;
            if nbsv < 15 {
                nbsv = 15;
            }
            if nbsv > A_MAX_NB_SAMPLE {
                nbsv = A_MAX_NB_SAMPLE;
            }
        }
        Surface3::Revolution(_) => {
            // OCCT GeomAbs_SurfaceOfRevolution (L375-378).
            nbsu = 15;
            nbsv = 15;
        }
        Surface3::BSpline(b) => {
            // OCCT GeomAbs_BSplineSurface (L318-359): nbsv = NbVKnots * VDegree.
            let n_u_knots = distinct_knot_count(&b.knots_u);
            let n_v_knots = distinct_knot_count(&b.knots_v);
            nbsv = (n_v_knots as i32).saturating_mul(b.degree_v as i32);
            if nbsv < 4 {
                nbsv = 4;
            }
            nbsu = (n_u_knots as i32).saturating_mul(b.degree_u as i32);
            if nbsu < 4 {
                nbsu = 4;
            }
            if nbsu < 10 {
                nbsu = 10;
            }
            if nbsv < 10 {
                nbsv = 10;
            }
            if nbsu > A_MAX_NB_SAMPLE {
                nbsu = A_MAX_NB_SAMPLE;
            }
            if nbsv > A_MAX_NB_SAMPLE {
                nbsv = A_MAX_NB_SAMPLE;
            }
        }
        Surface3::Bezier(b) => {
            // OCCT GeomAbs_BezierSurface (L298-316): nbsv = 3 + NbVPoles.
            let n_u_poles = b.control_points.first().map_or(0, |row| row.len());
            let n_v_poles = b.control_points.len();
            nbsv = (3 + n_v_poles) as i32;
            nbsu = (3 + n_u_poles) as i32;
            if nbsu < 10 {
                nbsu = 10;
            }
            if nbsv < 10 {
                nbsv = 10;
            }
            if nbsu > A_MAX_NB_SAMPLE {
                nbsu = A_MAX_NB_SAMPLE;
            }
            if nbsv > A_MAX_NB_SAMPLE {
                nbsv = A_MAX_NB_SAMPLE;
            }
        }
        _ => {
            // OCCT default (L380-384).
            nbsu = 10;
            nbsv = 10;
        }
    }
    (nbsu.max(1) as usize, nbsv.max(1) as usize)
}

/// Number of distinct knot values in a (multiplicity-expanded) knot vector.
fn distinct_knot_count(knots: &[f64]) -> usize {
    if knots.is_empty() {
        return 0;
    }
    let mut count = 1usize;
    for i in 1..knots.len() {
        if (knots[i] - knots[i - 1]).abs() > 1e-12 {
            count += 1;
        }
    }
    count
}

/// OCCT IntPatch_ImpPrmIntersection (IntPatch_ImpPrmIntersection.hxx L32-90).
pub struct ImpPrmIntersection {
    done: bool,
    empt: bool,
    spnt: Vec<IntPatchVertex>,
    slin: Vec<IntPatchLine>,
    solrst: SOnBounds,
    solins: SearchInside,
    my_is_start_pnt: bool,
    my_u_start: f64,
    my_v_start: f64,
}

impl ImpPrmIntersection {
    /// OCCT default constructor (L181-188).
    pub fn new() -> Self {
        ImpPrmIntersection {
            done: false,
            empt: false,
            spnt: Vec::new(),
            slin: Vec::new(),
            solrst: SOnBounds::new(),
            solins: SearchInside::new(),
            my_is_start_pnt: false,
            my_u_start: 0.0,
            my_v_start: 0.0,
        }
    }

    /// OCCT SetStartPoint (L212-217).
    pub fn set_start_point(&mut self, u: f64, v: f64) {
        self.my_is_start_pnt = true;
        self.my_u_start = u;
        self.my_v_start = v;
    }

    /// OCCT IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT IsEmpty.
    pub fn is_empty(&self) -> bool {
        self.empt
    }
    /// OCCT NbPnts.
    pub fn nb_points(&self) -> usize {
        self.spnt.len()
    }
    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> &IntPatchVertex {
        &self.spnt[index - 1]
    }
    /// OCCT NbLines.
    pub fn nb_lines(&self) -> usize {
        self.slin.len()
    }
    /// OCCT Line(Index) — 1-based.
    pub fn line(&self, index: usize) -> &IntPatchLine {
        &self.slin[index - 1]
    }

    /// OCCT Perform (L617-1964).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        uv1: [f64; 4],
        uv2: [f64; 4],
        tol_arc: f64,
        tol_tang: f64,
        fleche: f64,
        pas: f64,
    ) {
        let mut reversed = false;
        let mut procf = false;
        let mut procl = false;
        let mut dofirst = false;
        let mut dolast = false;
        let mut indfirst = 0usize;
        let mut indlast = 0usize;
        let mut ind2 = 0usize;
        let mut paramf = 0.0;
        let mut paraml = 0.0;
        let mut currentparam = 0.0;
        let mut u1 = 0.0;
        let mut v1 = 0.0;
        let mut u2 = 0.0;
        let mut v2 = 0.0;

        let type_s1 = classify_alias(s1);
        let type_s2 = classify_alias(s2);

        let mut trans1 = TypeTrans::Undecided;
        let mut trans2 = TypeTrans::Undecided;

        self.done = false;
        self.empt = true;
        self.slin.clear();
        self.spnt.clear();

        // Identify the quadric surface.
        let mut quad = Quadric::new();
        match type_s1 {
            GeomAbsSurfaceTypeAlias::Plane
            | GeomAbsSurfaceTypeAlias::Cylinder
            | GeomAbsSurfaceTypeAlias::Sphere
            | GeomAbsSurfaceTypeAlias::Cone => {
                quad = quad_of(s1);
            }
            _ => {
                reversed = true;
                match type_s2 {
                    GeomAbsSurfaceTypeAlias::Plane
                    | GeomAbsSurfaceTypeAlias::Cylinder
                    | GeomAbsSurfaceTypeAlias::Sphere
                    | GeomAbsSurfaceTypeAlias::Cone => {
                        quad = quad_of(s2);
                    }
                    _ => {
                        // OCCT throws Standard_ConstructionError.
                        self.done = true;
                        self.empt = true;
                        return;
                    }
                }
            }
        }

        let a_local_pas = if reversed {
            get_local_step(s1, pas)
        } else {
            get_local_step(s2, pas)
        };

        let mut func = SurfFunction::with_quadric(quad.clone());
        func.set_implicit_surface(quad.clone());
        func.set_tolerance(quadric_tool_tolerance(&quad));
        let mut a_func = ArcFunction::new();
        a_func.set_quadric(quad.clone());

        if !reversed {
            func.set_surface(s2.clone());
            a_func.set_surface(s2.clone());
        } else {
            func.set_surface(s1.clone());
            a_func.set_surface(s1.clone());
        }

        // SOnBounds.
        let p_domain = if !reversed { uv2 } else { uv1 };
        let mut dom = Domain::new(p_domain[0], p_domain[1], p_domain[2], p_domain[3]);
        self.solrst.perform(&mut a_func, &mut dom, tol_arc, tol_tang, false);
        if !self.solrst.is_done() {
            return;
        }

        // ComputeTangency → seqpdep.
        let mut seqpdep: Vec<PathPoint> = Vec::new();
        let nb_point_rst = self.solrst.nb_points();
        let mut destination = vec![0i32; nb_point_rst + 1];
        if nb_point_rst > 0 {
            let p_surf = if !reversed { s2 } else { s1 };
            compute_tangency(&self.solrst, &mut seqpdep, p_surf, &mut func, &mut destination, &mut dom);
        }

        // Decide whether SearchInside is needed.
        let mut search_ins = true;
        if quad.type_quadric() == QuadricType::Plane && self.solrst.nb_segments() > 0 {
            // OCCT L772-811: for a plane quadric with boundary segments the
            // parametric surface may lie entirely on one side of the plane
            // (only touching it), so no inner points exist.  The TopolTool grid
            // (NbSamples()/SamplePoint) samples the parametric face UV
            // rectangle to detect a sign change of F(u,v) — the grid covers the
            // whole face, not just the 4 corners.
            search_ins = false;
            // OCCT: T = reversed ? D1 : D2, the TopolTool of the parametric
            // face, whose UV domain is the corrected FF rectangle (uv1/uv2).
            let p_surf = if !reversed { s2 } else { s1 };
            let p_domain = if !reversed { uv2 } else { uv1 };
            let (nbsu, nbsv) = topol_tool_nb_samples(p_surf, p_domain);
            let du = (p_domain[1] - p_domain[0]) / (nbsu as f64 + 1.0);
            let dv = (p_domain[3] - p_domain[2]) / (nbsv as f64 + 1.0);
            let mut rvalf = 0.0f64;
            let mut first = true;
            // OCCT SamplePoint(Index): iv = 1 + Index/NbSamplesU, iu = 1 +
            // Index - (iv-1)*NbSamplesU; u = U0 + iu*DU, v = V0 + iv*DV.
            'outer: for iv in 1..=nbsv {
                for iu in 1..=nbsu {
                    let u = p_domain[0] + iu as f64 * du;
                    let v = p_domain[2] + iv as f64 * dv;
                    let x = [u, v];
                    if let Some(val) = func.value(&x) {
                        if first {
                            rvalf = val.copysign(1.0);
                            first = false;
                        } else if rvalf * val < 0.0 {
                            search_ins = true;
                            break 'outer;
                        }
                    }
                }
            }
        }

        // Interior points.
        let mut seqpins: Vec<InteriorPoint> = Vec::new();
        let mut nb_point_ins = 0usize;
        if search_ins {
            let p_surf = if !reversed { s2 } else { s1 };
            let p_domain = if !reversed { uv2 } else { uv1 };
            if self.my_is_start_pnt {
                self.solins.perform_from_point(
                    &mut func,
                    p_surf,
                    p_domain,
                    self.my_u_start,
                    self.my_v_start,
                );
            } else {
                self.solins.perform(&mut func, p_surf, p_domain, tol_tang);
            }
            nb_point_ins = self.solins.nb_points();
            for i in 1..=nb_point_ins {
                seqpins.push(self.solins.value(i).clone());
            }
        }

        let nb_point_dep = seqpdep.len();

        if nb_point_dep > 0 || nb_point_ins > 0 {
            let param_surf = if reversed { s1 } else { s2 };
            let p_domain = if !reversed { uv2 } else { uv1 };
            let mut iwalk = IWalking::new(tol_tang, fleche, a_local_pas, false);
            iwalk.perform(&seqpdep, &seqpins, &mut func, param_surf, p_domain, reversed);

            if !iwalk.is_done() {
                return;
            }

            let (vmin, vmax) = if !reversed {
                (uv1[2], uv1[3])
            } else {
                (uv2[2], uv2[3])
            };
            let tol_v = 1e-14;

            let nblines = iwalk.nb_lines();
            for j in 1..=nblines {
                let iwline = iwalk.value(j);
                let thelin = iwline.line().clone();

                let nbpts = thelin.nb_points();
                if nbpts >= 2 {
                    // OCCT L878-886: TangentVector(k) sets k = indextg (1-based);
                    // the point at the 1-based index k is used for the normal /
                    // transition computation.  rcad's thelin is 0-based, so read
                    // thelin.value(k - 1).
                    let mut k = 0usize;
                    let (mut tgline, k_tg) = iwline.tangent_vector();
                    k = k_tg.max(0) as usize;
                    if k >= 1 && k <= nbpts {
                    } else {
                        k = nbpts >> 1;
                    }
                    let valpt = thelin.value(k - 1).value();

                    let (u2v, v2v);
                    let norm1;
                    let norm2;
                    if !reversed {
                        let (u, v) = thelin.value(k - 1).parameters_on_surface(false);
                        (u2v, v2v) = (u, v);
                        norm1 = quad.normale(valpt);
                        let (_, d1u, d1v) = s2.derivatives(u2v, v2v);
                        norm2 = d1u.cross(d1v);
                    } else {
                        let (u, v) = thelin.value(k - 1).parameters_on_surface(true);
                        (u2v, v2v) = (u, v);
                        norm2 = quad.normale(valpt);
                        let (_, d1u, d1v) = s1.derivatives(u2v, v2v);
                        norm1 = d1u.cross(d1v);
                    }
                    if tgline.dot(norm2.cross(norm1)) > 0.0 {
                        trans1 = TypeTrans::Out;
                        trans2 = TypeTrans::In;
                    } else {
                        trans1 = TypeTrans::In;
                        trans2 = TypeTrans::Out;
                    }

                    let mut an_u1;
                    let mut an_u2;
                    let mut an_v2;
                    let mut v1q;

                    let typ_quad = quad.type_quadric();
                    let mut arecadr = false;
                    let valpt = thelin.value(0).value();
                    let (an_u1o, v1o) = quad.parameters(valpt);
                    an_u1 = an_u1o;
                    v1q = v1o;

                    if (v1q < vmin) && (vmin - v1q < tol_v) {
                        v1q = vmin;
                    }
                    if (v1q > vmax) && (v1q - vmax < tol_v) {
                        v1q = vmax;
                    }

                    let mut thelin = thelin;
                    if reversed {
                        thelin.set_uv(0, false, an_u1, v1q);
                        let (u, v) = thelin.value(0).parameters_on_surface(true);
                        (an_u2, an_v2) = (u, v);
                    } else {
                        thelin.set_uv(0, true, an_u1, v1q);
                        let (u, v) = thelin.value(0).parameters_on_surface(false);
                        (an_u2, an_v2) = (u, v);
                    }

                    if typ_quad == QuadricType::Cylinder
                        || typ_quad == QuadricType::Cone
                        || typ_quad == QuadricType::Sphere
                    {
                        arecadr = true;
                    }

                    for k in 1..nbpts {
                        let valpt = thelin.value(k).value();
                        let (u, v) = quad.parameters(valpt);
                        u1 = u;
                        v1q = v;

                        if (v1q < vmin) && (vmin - v1q < tol_v) {
                            v1q = vmin;
                        }
                        if (v1q > vmax) && (v1q - vmax < tol_v) {
                            v1q = vmax;
                        }

                        if arecadr {
                            let a_cf = 0.0f64;
                            let a_two_pi = TAU;
                            if (u1 - an_u1) > 1.5 * std::f64::consts::PI {
                                let mut cf = a_cf;
                                while (u1 - an_u1) > (1.5 * std::f64::consts::PI + cf * a_two_pi) {
                                    cf += 1.0;
                                }
                                u1 -= cf * a_two_pi;
                            } else {
                                let mut cf = a_cf;
                                while (u1 - an_u1) < (-1.5 * std::f64::consts::PI - cf * a_two_pi) {
                                    cf += 1.0;
                                }
                                u1 += cf * a_two_pi;
                            }
                        }

                        if reversed {
                            thelin.set_uv(k, false, u1, v1q);
                            let (u, v) = thelin.value(k).parameters_on_surface(true);
                            u2 = u;
                            v2 = v;
                            match type_s1 {
                                GeomAbsSurfaceTypeAlias::Cylinder
                                | GeomAbsSurfaceTypeAlias::Cone
                                | GeomAbsSurfaceTypeAlias::Sphere
                                | GeomAbsSurfaceTypeAlias::Torus => {
                                    while u2 < an_u2 - 1.5 * std::f64::consts::PI {
                                        u2 += TAU;
                                    }
                                    while u2 > an_u2 + 1.5 * std::f64::consts::PI {
                                        u2 -= TAU;
                                    }
                                }
                                _ => {}
                            }
                            if type_s2 == GeomAbsSurfaceTypeAlias::Torus {
                                while v2 < an_v2 - 1.5 * std::f64::consts::PI {
                                    v2 += TAU;
                                }
                                while v2 > an_v2 + 1.5 * std::f64::consts::PI {
                                    v2 -= TAU;
                                }
                            }
                            thelin.set_uv(k, true, u2, v2);
                        } else {
                            thelin.set_uv(k, true, u1, v1q);
                            let (u, v) = thelin.value(k).parameters_on_surface(false);
                            u2 = u;
                            v2 = v;
                            match type_s2 {
                                GeomAbsSurfaceTypeAlias::Cylinder
                                | GeomAbsSurfaceTypeAlias::Cone
                                | GeomAbsSurfaceTypeAlias::Sphere
                                | GeomAbsSurfaceTypeAlias::Torus => {
                                    while u2 < an_u2 - 1.5 * std::f64::consts::PI {
                                        u2 += TAU;
                                    }
                                    while u2 > an_u2 + 1.5 * std::f64::consts::PI {
                                        u2 -= TAU;
                                    }
                                }
                                _ => {}
                            }
                            if type_s2 == GeomAbsSurfaceTypeAlias::Torus {
                                while v2 < an_v2 - 1.5 * std::f64::consts::PI {
                                    v2 += TAU;
                                }
                                while v2 > an_v2 + 1.5 * std::f64::consts::PI {
                                    v2 -= TAU;
                                }
                            }
                            thelin.set_uv(k, false, u2, v2);
                        }

                        an_u1 = u1;
                        an_u2 = u2;
                        an_v2 = v2;
                    }

                    let wline_pnts: Vec<crate::geomalgo::int_patch::WLinePnt> = (0..thelin.nb_points())
                        .map(|i| {
                            let p = thelin.value(i);
                            let (u1, v1, u2, v2) = p.parameters();
                            crate::geomalgo::int_patch::WLinePnt {
                                p3d: p.value(),
                                u1,
                                v1,
                                u2,
                                v2,
                            }
                        })
                        .collect();
                    let mut wline = IntPatchLine::walking(wline_pnts, crate::geomalgo::int_patch::WLineType::ImpPrm);
                    wline.trans1 = Some(Transition::from_type(trans1));
                    wline.trans2 = Some(Transition::from_type(trans2));

                    // Vertex at the first point.
                    if iwline.has_first_point() && !iwline.is_tangent_at_begining() {
                        indfirst = iwline.first_point_index() as usize;
                        let ppoint = seqpdep[indfirst - 1].clone();
                        tgline = ppoint.direction_3d();
                        let mut themult = ppoint.multiplicity();
                        let mut i = nb_point_rst;
                        while i >= 1 {
                            if destination[i - 1] as usize == indfirst {
                                if !reversed {
                                    let (u, v) = quad.parameters(ppoint.value());
                                    u1 = u;
                                    v1q = v;
                                    if (v1q < vmin) && (vmin - v1q < tol_v) {
                                        v1q = vmin;
                                    }
                                    if (v1q > vmax) && (v1q - vmax < tol_v) {
                                        v1q = vmax;
                                    }
                                    ppoint.parameters(themult, &mut u2, &mut v2);
                                    let (_, d1u, d1v) = s2.derivatives(u2, v2);
                                    let vec_normale = d1u.cross(d1v);
                                    let mut ptdeb = IntPatchVertex::default();
                                    ptdeb.set_value(ppoint.value(), tol_arc, false);
                                    ptdeb.set_parameters(u1, v1q, u2, v2);
                                    ptdeb.set_parameter(1.0);
                                    recadre_impprm(
                                        type_s1,
                                        type_s2,
                                        &mut ptdeb,
                                        &thelin_points(&thelin),
                                        0,
                                        u1,
                                        v1q,
                                        u2,
                                        v2,
                                    );
                                    let current_arc = self.solrst.point(i).arc().clone();
                                    currentparam = self.solrst.point(i).parameter();
                                    let p2d = current_arc.point_at(currentparam);
                                    let d2d = current_arc.derivative_at(currentparam);
                                    let (_, d1u2, d1v2) = if !reversed { s2.derivatives(u2, v2) } else { s1.derivatives(u1, v1q) };
                                    let tgrst = d2d.x * d1u2 + d2d.y * d1v2;
                                    let mut t_line = Transition::new();
                                    let mut t_arc = Transition::new();
                                    if vec_normale.length_squared() > 1e-13 {
                                        make_transition(tgline, tgrst, vec_normale, &mut t_line, &mut t_arc);
                                    } else {
                                        t_line.set_value_in_out(true, TypeTrans::Undecided);
                                        t_arc.set_value_in_out(true, TypeTrans::Undecided);
                                    }
                                    ptdeb.set_arc(reversed, current_arc, currentparam, t_line, t_arc);
                                    if !self.solrst.point(i).is_new() {
                                        ptdeb.set_vertex(reversed);
                                    }
                                    wline.add_vertex(ptdeb);
                                    if themult == 0 {
                                        wline.set_first_point(wline.nb_vertex());
                                    }
                                    themult -= 1;
                                } else {
                                    let (u, v) = quad.parameters(ppoint.value());
                                    u2 = u;
                                    v2 = v;
                                    if (v2 < vmin) && (vmin - v2 < tol_v) {
                                        v2 = vmin;
                                    }
                                    if (v2 > vmax) && (v2 - vmax < tol_v) {
                                        v2 = vmax;
                                    }
                                    ppoint.parameters(themult, &mut u1, &mut v1);
                                    let (_, d1u, d1v) = s1.derivatives(u1, v1);
                                    let vec_normale = d1u.cross(d1v);
                                    let mut ptdeb = IntPatchVertex::default();
                                    ptdeb.set_value(ppoint.value(), tol_arc, false);
                                    ptdeb.set_parameters(u1, v1, u2, v2);
                                    ptdeb.set_parameter(1.0);
                                    recadre_impprm(
                                        type_s1,
                                        type_s2,
                                        &mut ptdeb,
                                        &thelin_points(&thelin),
                                        0,
                                        u1,
                                        v1,
                                        u2,
                                        v2,
                                    );
                                    let current_arc = self.solrst.point(i).arc().clone();
                                    currentparam = self.solrst.point(i).parameter();
                                    let p2d = current_arc.point_at(currentparam);
                                    let d2d = current_arc.derivative_at(currentparam);
                                    let (_, d1u2, d1v2) = s1.derivatives(u1, v1);
                                    let tgrst = d2d.x * d1u2 + d2d.y * d1v2;
                                    let mut t_line = Transition::new();
                                    let mut t_arc = Transition::new();
                                    if vec_normale.length_squared() > 1e-13 {
                                        make_transition(tgline, tgrst, vec_normale, &mut t_line, &mut t_arc);
                                    } else {
                                        t_line.set_value_in_out(true, TypeTrans::Undecided);
                                        t_arc.set_value_in_out(true, TypeTrans::Undecided);
                                    }
                                    ptdeb.set_arc(reversed, current_arc, currentparam, t_line, t_arc);
                                    if !self.solrst.point(i).is_new() {
                                        ptdeb.set_vertex(reversed);
                                    }
                                    wline.add_vertex(ptdeb);
                                    if themult == 0 {
                                        wline.set_first_point(wline.nb_vertex());
                                    }
                                    themult -= 1;
                                }
                            }
                            i = i.saturating_sub(1);
                        }
                    } else if iwline.is_tangent_at_begining() {
                        let psol = thelin.value(0).value();
                        let (u, v) = thelin.value(0).parameters_on_surface(true);
                        u1 = u;
                        v1 = v;
                        let (u, v) = thelin.value(0).parameters_on_surface(false);
                        u2 = u;
                        v2 = v;
                        let mut ptdeb = IntPatchVertex::default();
                        ptdeb.set_value(psol, tol_arc, true);
                        ptdeb.set_parameters(u1, v1, u2, v2);
                        ptdeb.set_parameter(1.0);
                        wline.add_vertex(ptdeb);
                        wline.set_first_point(wline.nb_vertex());
                    } else {
                        let psol = thelin.value(0).value();
                        let (u, v) = thelin.value(0).parameters_on_surface(true);
                        u1 = u;
                        v1 = v;
                        let (u, v) = thelin.value(0).parameters_on_surface(false);
                        u2 = u;
                        v2 = v;
                        let mut ptdeb = IntPatchVertex::default();
                        ptdeb.set_value(psol, tol_arc, false);
                        ptdeb.set_parameters(u1, v1, u2, v2);
                        ptdeb.set_parameter(1.0);
                        wline.add_vertex(ptdeb);
                        wline.set_first_point(wline.nb_vertex());
                    }

                    // Vertex at the last point.
                    if iwline.has_last_point() && !iwline.is_tangent_at_end() {
                        indlast = iwline.last_point_index() as usize;
                        let ppoint = seqpdep[indlast - 1].clone();
                        tgline = -ppoint.direction_3d();
                        let mut themult = ppoint.multiplicity();
                        let mut i = nb_point_rst;
                        while i >= 1 {
                            if destination[i - 1] as usize == indlast {
                                if !reversed {
                                    let (u, v) = quad.parameters(ppoint.value());
                                    u1 = u;
                                    v1q = v;
                                    if (v1q < vmin) && (vmin - v1q < tol_v) {
                                        v1q = vmin;
                                    }
                                    if (v1q > vmax) && (v1q - vmax < tol_v) {
                                        v1q = vmax;
                                    }
                                    ppoint.parameters(themult, &mut u2, &mut v2);
                                    let (_, d1u, d1v) = s2.derivatives(u2, v2);
                                    let vec_normale = d1u.cross(d1v);
                                    let mut ptfin = IntPatchVertex::default();
                                    ptfin.set_value(ppoint.value(), tol_arc, false);
                                    ptfin.set_parameters(u1, v1q, u2, v2);
                                    ptfin.set_parameter(nbpts as f64);
                                    recadre_impprm(
                                        type_s1,
                                        type_s2,
                                        &mut ptfin,
                                        &thelin_points(&thelin),
                                        nbpts - 2,
                                        u1,
                                        v1q,
                                        u2,
                                        v2,
                                    );
                                    let current_arc = self.solrst.point(i).arc().clone();
                                    currentparam = self.solrst.point(i).parameter();
                                    let p2d = current_arc.point_at(currentparam);
                                    let d2d = current_arc.derivative_at(currentparam);
                                    let (_, d1u2, d1v2) = s2.derivatives(u2, v2);
                                    let tgrst = d2d.x * d1u2 + d2d.y * d1v2;
                                    let mut t_line = Transition::new();
                                    let mut t_arc = Transition::new();
                                    if vec_normale.length_squared() > 1e-13 {
                                        make_transition(tgline, tgrst, vec_normale, &mut t_line, &mut t_arc);
                                    } else {
                                        t_line.set_value_in_out(true, TypeTrans::Undecided);
                                        t_arc.set_value_in_out(true, TypeTrans::Undecided);
                                    }
                                    ptfin.set_arc(reversed, current_arc, currentparam, t_line, t_arc);
                                    if !self.solrst.point(i).is_new() {
                                        ptfin.set_vertex(reversed);
                                    }
                                    wline.add_vertex(ptfin);
                                    if themult == 0 {
                                        wline.set_last_point(wline.nb_vertex());
                                    }
                                    themult -= 1;
                                } else {
                                    let (u, v) = quad.parameters(ppoint.value());
                                    u2 = u;
                                    v2 = v;
                                    if (v2 < vmin) && (vmin - v2 < tol_v) {
                                        v2 = vmin;
                                    }
                                    if (v2 > vmax) && (v2 - vmax < tol_v) {
                                        v2 = vmax;
                                    }
                                    ppoint.parameters(themult, &mut u1, &mut v1);
                                    let (_, d1u, d1v) = s1.derivatives(u1, v1);
                                    let vec_normale = d1u.cross(d1v);
                                    let mut ptfin = IntPatchVertex::default();
                                    ptfin.set_value(ppoint.value(), tol_arc, false);
                                    ptfin.set_parameters(u1, v1, u2, v2);
                                    ptfin.set_parameter(nbpts as f64);
                                    recadre_impprm(
                                        type_s1,
                                        type_s2,
                                        &mut ptfin,
                                        &thelin_points(&thelin),
                                        nbpts - 2,
                                        u1,
                                        v1,
                                        u2,
                                        v2,
                                    );
                                    let current_arc = self.solrst.point(i).arc().clone();
                                    currentparam = self.solrst.point(i).parameter();
                                    let p2d = current_arc.point_at(currentparam);
                                    let d2d = current_arc.derivative_at(currentparam);
                                    let (_, d1u2, d1v2) = s1.derivatives(u1, v1);
                                    let tgrst = d2d.x * d1u2 + d2d.y * d1v2;
                                    let mut t_line = Transition::new();
                                    let mut t_arc = Transition::new();
                                    if vec_normale.length_squared() > 1e-13 {
                                        make_transition(tgline, tgrst, vec_normale, &mut t_line, &mut t_arc);
                                    } else {
                                        t_line.set_value_in_out(true, TypeTrans::Undecided);
                                        t_arc.set_value_in_out(true, TypeTrans::Undecided);
                                    }
                                    ptfin.set_arc(reversed, current_arc, currentparam, t_line, t_arc);
                                    if !self.solrst.point(i).is_new() {
                                        ptfin.set_vertex(reversed);
                                    }
                                    wline.add_vertex(ptfin);
                                    if themult == 0 {
                                        wline.set_last_point(wline.nb_vertex());
                                    }
                                    themult -= 1;
                                }
                            }
                            i = i.saturating_sub(1);
                        }
                    } else if iwline.is_tangent_at_end() {
                        let psol = thelin.value(nbpts - 1).value();
                        let (u, v) = thelin.value(nbpts - 1).parameters_on_surface(true);
                        u1 = u;
                        v1 = v;
                        let (u, v) = thelin.value(nbpts - 1).parameters_on_surface(false);
                        u2 = u;
                        v2 = v;
                        let mut ptfin = IntPatchVertex::default();
                        ptfin.set_value(psol, tol_arc, true);
                        ptfin.set_parameters(u1, v1, u2, v2);
                        ptfin.set_parameter(nbpts as f64);
                        wline.add_vertex(ptfin);
                        wline.set_last_point(wline.nb_vertex());
                    } else {
                        let psol = thelin.value(nbpts - 1).value();
                        let (u, v) = thelin.value(nbpts - 1).parameters_on_surface(true);
                        u1 = u;
                        v1 = v;
                        let (u, v) = thelin.value(nbpts - 1).parameters_on_surface(false);
                        u2 = u;
                        v2 = v;
                        let mut ptfin = IntPatchVertex::default();
                        ptfin.set_value(psol, tol_arc, false);
                        ptfin.set_parameters(u1, v1, u2, v2);
                        ptfin.set_parameter(nbpts as f64);
                        wline.add_vertex(ptfin);
                        wline.set_last_point(wline.nb_vertex());
                    }

                    // Il faut traiter les points de passage.
                    self.slin.push(wline);
                }
            }

            // Connect tangent points between lines.
            let nblines = self.slin.len();
            let mut j = 1usize;
            while j <= nblines.saturating_sub(1) {
                let mut dofirst = false;
                let mut dolast = false;
                let mut ptdeb = IntPatchVertex::default();
                let mut ptfin = IntPatchVertex::default();
                {
                    let slinj = &self.slin[j - 1];
                    if slinj.has_first_point() {
                        let fp = slinj.first_point().clone();
                        if fp.is_tangency_point() {
                            dofirst = true;
                            ptdeb = fp;
                        }
                    }
                    if slinj.has_last_point() {
                        let lp = slinj.last_point().clone();
                        if lp.is_tangency_point() {
                            dolast = true;
                            ptfin = lp;
                        }
                    }
                }
                if dofirst || dolast {
                    let mut k = j + 1;
                    while k <= nblines {
                        let (fp_idx, fp2, lp_idx, lp2) = {
                            let slink = &self.slin[k - 1];
                            let fp2 = slink.has_first_point().then(|| slink.first_point().clone());
                            let lp2 = slink.has_last_point().then(|| slink.last_point().clone());
                            (slink.first_point, fp2, slink.last_point, lp2)
                        };
                        if let Some(ptbis) = fp2 {
                            if ptbis.is_tangency_point() {
                                if dofirst {
                                    if ptdeb.p3d.distance(ptbis.p3d) <= tol_arc {
                                        ptdeb.set_multiple(true);
                                        if !ptbis.is_multiple() {
                                            let mut ptbis_m = ptbis.clone();
                                            ptbis_m.set_multiple(true);
                                            self.slin[k - 1].replace_vertex(fp_idx.unwrap(), ptbis_m);
                                        }
                                    }
                                }
                                if dolast {
                                    if ptfin.p3d.distance(ptbis.p3d) <= tol_arc {
                                        ptfin.set_multiple(true);
                                        if !ptbis.is_multiple() {
                                            let mut ptbis_m = ptbis.clone();
                                            ptbis_m.set_multiple(true);
                                            self.slin[k - 1].replace_vertex(fp_idx.unwrap(), ptbis_m);
                                        }
                                    }
                                }
                            }
                        }
                        if let Some(ptbis) = lp2 {
                            if ptbis.is_tangency_point() {
                                if dofirst {
                                    if ptdeb.p3d.distance(ptbis.p3d) <= tol_arc {
                                        ptdeb.set_multiple(true);
                                        if !ptbis.is_multiple() {
                                            let mut ptbis_m = ptbis.clone();
                                            ptbis_m.set_multiple(true);
                                            self.slin[k - 1].replace_vertex(lp_idx.unwrap(), ptbis_m);
                                        }
                                    }
                                }
                                if dolast {
                                    if ptfin.p3d.distance(ptbis.p3d) <= tol_arc {
                                        ptfin.set_multiple(true);
                                        if !ptbis.is_multiple() {
                                            let mut ptbis_m = ptbis.clone();
                                            ptbis_m.set_multiple(true);
                                            self.slin[k - 1].replace_vertex(lp_idx.unwrap(), ptbis_m);
                                        }
                                    }
                                }
                            }
                        }
                        k += 1;
                    }
                    if dofirst {
                        if let Some(fp_idx) = self.slin[j - 1].first_point {
                            self.slin[j - 1].replace_vertex(fp_idx, ptdeb);
                        }
                    }
                    if dolast {
                        if let Some(lp_idx) = self.slin[j - 1].last_point {
                            self.slin[j - 1].replace_vertex(lp_idx, ptfin);
                        }
                    }
                }
                j += 1;
            }
        }

        // Treatment of the segments.
        let nb_segm = self.solrst.nb_segments();
        if nb_segm > 0 {
            for i in 1..=nb_segm {
                let thesegm = self.solrst.segment(i).clone();
                // Check if segment is degenerated.
                if thesegm.has_first_point() && thesegm.has_last_point() {
                    let tol2 = rcad_kernel::precision::CONFUSION;
                    let tol2 = tol2 * tol2;
                    let a_pf = thesegm.first_point().value();
                    let a_pl = thesegm.last_point().value();
                    if a_pf.distance_squared(a_pl) <= tol2 {
                        // Segment can be degenerated — check inner point.
                        paramf = thesegm.first_point().parameter();
                        paraml = thesegm.last_point().parameter();
                        let p2d = thesegm.curve().point_at(0.57735 * paramf + 0.42265 * paraml);
                        let a_pm = if reversed {
                            s1.point_at(p2d.x, p2d.y)
                        } else {
                            s2.point_at(p2d.x, p2d.y)
                        };
                        if a_pm.distance_squared(a_pf) <= tol2 {
                            // Degenerated.
                            continue;
                        }
                    }
                }

                dofirst = false;
                dolast = false;
                procf = false;
                procl = false;

                let mut the_point_at_beg = IntPatchVertex::default();
                let mut the_point_at_end = IntPatchVertex::default();
                let mut transition_ok = false;
                let mut trans1 = TypeTrans::Undecided;
                let mut trans2 = TypeTrans::Undecided;

                if thesegm.has_first_point() {
                    dofirst = true;
                    let p_startf = thesegm.first_point().clone();
                    paramf = p_startf.parameter();
                    let p2d = thesegm.curve().point_at(paramf);
                    let pp = p_startf.value();
                    the_point_at_beg.set_value(pp, p_startf.tolerance(), false);
                    if !reversed {
                        let (u, v) = quad.parameters(pp);
                        u1 = u;
                        v1 = v;
                        u2 = p2d.x;
                        v2 = p2d.y;
                    } else {
                        let (u, v) = quad.parameters(pp);
                        u2 = u;
                        v2 = v;
                        u1 = p2d.x;
                        v1 = p2d.y;
                    }
                    the_point_at_beg.set_parameters(u1, v1, u2, v2);
                    the_point_at_beg.set_parameter(paramf);
                    if !p_startf.is_new() {
                        the_point_at_beg.set_vertex(reversed);
                    }
                    the_point_at_beg.set_arc(
                        reversed,
                        thesegm.curve().clone(),
                        paramf,
                        Transition::new(),
                        Transition::new(),
                    );

                    let (_, d1u1, d1v1) = s1.derivatives(u1, v1);
                    let norm1 = d1u1.cross(d1v1);
                    let (_, d1u2, d1v2) = s2.derivatives(u2, v2);
                    let norm2 = d1u2.cross(d1v2);
                    let d2d = thesegm.curve().derivative_at(paramf);
                    let tgline = if reversed {
                        d2d.x * d1u1 + d2d.y * d1v1
                    } else {
                        d2d.x * d1u2 + d2d.y * d1v2
                    };
                    let u1v = tgline.dot(norm2.cross(norm1));
                    transition_ok = true;
                    if u1v > 0.00000001 {
                        trans1 = TypeTrans::Out;
                        trans2 = TypeTrans::In;
                    } else if u1v < -0.00000001 {
                        trans1 = TypeTrans::In;
                        trans2 = TypeTrans::Out;
                    } else {
                        transition_ok = false;
                    }
                }
                if thesegm.has_last_point() {
                    dolast = true;
                    let p_startl = thesegm.last_point().clone();
                    paraml = p_startl.parameter();
                    let p2d = thesegm.curve().point_at(paraml);
                    let pp = p_startl.value();
                    the_point_at_end.set_value(pp, p_startl.tolerance(), false);
                    if !reversed {
                        let (u, v) = quad.parameters(pp);
                        u1 = u;
                        v1 = v;
                        u2 = p2d.x;
                        v2 = p2d.y;
                    } else {
                        let (u, v) = quad.parameters(pp);
                        u2 = u;
                        v2 = v;
                        u1 = p2d.x;
                        v1 = p2d.y;
                    }
                    the_point_at_end.set_parameters(u1, v1, u2, v2);
                    the_point_at_end.set_parameter(paraml);
                    if !p_startl.is_new() {
                        the_point_at_end.set_vertex(reversed);
                    }
                    the_point_at_end.set_arc(
                        reversed,
                        thesegm.curve().clone(),
                        paraml,
                        Transition::new(),
                        Transition::new(),
                    );

                    let (_, d1u1, d1v1) = s1.derivatives(u1, v1);
                    let norm1 = d1u1.cross(d1v1);
                    let (_, d1u2, d1v2) = s2.derivatives(u2, v2);
                    let norm2 = d1u2.cross(d1v2);
                    let d2d = thesegm.curve().derivative_at(paraml);
                    let tgline = if reversed {
                        d2d.x * d1u1 + d2d.y * d1v1
                    } else {
                        d2d.x * d1u2 + d2d.y * d1v2
                    };
                    let u1v = tgline.dot(norm2.cross(norm1));
                    transition_ok = true;
                    if u1v > 0.00000001 {
                        trans1 = TypeTrans::Out;
                        trans2 = TypeTrans::In;
                    } else if u1v < -0.00000001 {
                        trans1 = TypeTrans::In;
                        trans2 = TypeTrans::Out;
                    } else {
                        transition_ok = false;
                    }
                }

                // Create the RLine.
                let mut rline = IntPatchLine::analytic(IntPatchIType::Restriction, rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }), [0.0, 1.0]);
                rline.line_type = IntPatchIType::Restriction;
                // OCCT L1599-1623: the RLine carries the transitions when
                // TransitionOK (IntPatch_RLine(false, trans1, trans2)); otherwise
                // the bare RLine constructor is used and no transitions are set.
                if transition_ok {
                    rline.trans1 = Some(Transition::from_type(trans1));
                    rline.trans2 = Some(Transition::from_type(trans2));
                }
                if reversed {
                    rline.set_arc_on_s1(thesegm.curve().clone());
                } else {
                    rline.set_arc_on_s2(thesegm.curve().clone());
                }

                if thesegm.has_first_point() {
                    rline.add_vertex(the_point_at_beg);
                    rline.set_first_point(rline.nb_vertex());
                }
                if thesegm.has_last_point() {
                    rline.add_vertex(the_point_at_end);
                    rline.set_last_point(rline.nb_vertex());
                }

                // Polygone sur restriction solution.
                if dofirst && dolast {
                    let nbsample = 100usize;
                    for j in 1..=nbsample {
                        let prm = paramf + (j - 1) as f64 * (paraml - paramf) / (nbsample - 1) as f64;
                        let p2d = thesegm.curve().point_at(prm);
                        let ptpoly = if reversed { s1.point_at(p2d.x, p2d.y) } else { s2.point_at(p2d.x, p2d.y) };
                        let (u, v) = quad.parameters(ptpoly);
                        rline.wline_pnts.push(crate::geomalgo::int_patch::WLinePnt {
                            p3d: ptpoly,
                            u1: if reversed { p2d.x } else { u },
                            v1: if reversed { p2d.y } else { v },
                            u2: if reversed { u } else { p2d.x },
                            v2: if reversed { v } else { p2d.y },
                        });
                    }
                }

                // Attach nearby vertices from existing lines.
                if dofirst || dolast {
                    let nblines = self.slin.len();
                    let mut j = 1usize;
                    while j <= nblines {
                        let typ = self.slin[j - 1].line_type;
                        let nbpts = if typ == IntPatchIType::Walking {
                            self.slin[j - 1].nb_vertex()
                        } else {
                            self.slin[j - 1].nb_vertex()
                        };
                        let mut k = 1usize;
                        while k <= nbpts {
                            let ptdeb = self.slin[j - 1].vertex(k).clone();
                            if dofirst {
                                let p_startf = thesegm.first_point().clone();
                                if ptdeb.p3d.distance(p_startf.value()) <= tol_arc {
                                    let mut ptdeb_m = ptdeb.clone();
                                    ptdeb_m.set_multiple(true);
                                    self.slin[j - 1].replace_vertex(k, ptdeb_m.clone());
                                    ptdeb_m.set_parameter(paramf);
                                    rline.add_vertex(ptdeb_m);
                                    if !procf {
                                        procf = true;
                                        rline.set_first_point(rline.nb_vertex());
                                    }
                                }
                            }
                            if dolast {
                                let p_startl = thesegm.last_point().clone();
                                if ptdeb.p3d.distance(p_startl.value()) <= tol_arc {
                                    let mut ptdeb_m = ptdeb.clone();
                                    ptdeb_m.set_multiple(true);
                                    self.slin[j - 1].replace_vertex(k, ptdeb_m.clone());
                                    ptdeb_m.set_parameter(paraml);
                                    rline.add_vertex(ptdeb_m);
                                    if !procl {
                                        procl = true;
                                        rline.set_last_point(rline.nb_vertex());
                                    }
                                }
                            }
                            k += 1;
                        }
                        j += 1;
                    }
                }
                self.slin.push(rline);
            }
        }

        // On traite les restrictions de la surface implicite: remove short
        // lines (<= 2 coincident points) and move Walking lines to the end.
        // OCCT IntPatch_ImpPrmIntersection.cxx L1770-1808.
        let mut a_nb_lin = self.slin.len();
        let mut i: isize = 0;
        while (i as usize) < a_nb_lin {
            let is_wline = self.slin[i as usize].is_wline();
            let a_cond = if is_wline {
                let n = self.slin[i as usize].wline_pnts.len();
                if n < 2 {
                    true
                } else {
                    let p1 = self.slin[i as usize].wline_pnts[0].p3d;
                    let p2 = self.slin[i as usize].wline_pnts[1].p3d;
                    p1.distance_squared(p2) <= rcad_kernel::precision::SQUARE_CONFUSION
                }
            } else {
                self.slin[i as usize].vertices.len() < 2
            };
            if a_cond {
                self.slin.remove(i as usize);
                a_nb_lin -= 1; // OCCT: aNbLin--
                continue;
            }

            if self.slin[i as usize].line_type == IntPatchIType::Walking {
                let wl = self.slin.remove(i as usize);
                self.slin.push(wl);
                i -= 1; // OCCT: i--
                a_nb_lin -= 1; // OCCT: aNbLin-- (remove+push keeps len() unchanged)
            }
            i += 1; // OCCT: i++
        }

        // OCCT L1817-1910: IsCoincide loop.
        let a_tol_3d = func.tolerance().max(tol_tang);
        let an_other_surf = func.p_surface().clone();
        let mut i = 0usize;
        while i < self.slin.len() {
            // aL1 must be a Restriction line; Walking-Walking cases are not
            // supported (OCCT L1822-1826).
            if self.slin[i].line_type != IntPatchIType::Restriction {
                break;
            }
            let an_arc = self.slin[i]
                .arc_on_s1
                .as_ref()
                .or(self.slin[i].arc_on_s2.as_ref())
                .cloned();
            let Some(an_arc) = an_arc else { break };
            if !matches!(an_arc, Curve2d::Line(_)) {
                break; // Restriction line must be isoline.
            }
            let is_arc_on_s1 = self.slin[i].arc_on_s1.is_some();
            let a_arc_range = line_arc_range(&self.slin[i]);

            let mut is_first_deleted = false;
            let mut j = i + 1;
            while j < self.slin.len() {
                let is_rline2 = self.slin[j].line_type == IntPatchIType::Restriction;
                if is_rline2 {
                    let an_arc2 = self.slin[j]
                        .arc_on_s1
                        .as_ref()
                        .or(self.slin[j].arc_on_s2.as_ref());
                    let Some(an_arc2) = an_arc2 else { j += 1; continue };
                    if !matches!(an_arc2, Curve2d::Line(_)) {
                        j += 1;
                        continue;
                    }
                }

                // aDir can be one of following four values only (restriction
                // line is boundary of rectangular surface).
                let a_dir = match &an_arc {
                    Curve2d::Line(l) => l.direction,
                    _ => unreachable!(),
                };
                let mut a_tol_2d = u_res(&an_other_surf, a_tol_3d);
                let mut a_period = if an_other_surf.is_v_periodic() { TAU } else { 0.0 };
                if a_dir.x.abs() < 0.5 {
                    // Restriction directs along V-direction.
                    a_tol_2d = v_res(&an_other_surf, a_tol_3d);
                    a_period = if an_other_surf.is_u_periodic() { TAU } else { 0.0 };
                }

                let is_coincide = is_coincide(
                    &mut func,
                    &self.slin[j],
                    &an_arc,
                    a_arc_range,
                    is_arc_on_s1,
                    a_tol_3d,
                    a_tol_2d,
                    a_period,
                );
                if is_coincide {
                    if !is_rline2 {
                        // Delete Walking-line.
                        self.slin.remove(j);
                        continue;
                    }
                    // Restriction-Restriction: keep the longer.
                    let a_range2 = {
                        let r2 = line_arc_range(&self.slin[j]);
                        r2[1] - r2[0]
                    };
                    let a_range1 = a_arc_range[1] - a_arc_range[0];
                    if a_range2 > a_range1 {
                        is_first_deleted = true;
                        break;
                    }
                    self.slin.remove(j);
                    continue;
                }
                j += 1;
            }
            if is_first_deleted {
                // OCCT L1909-1911: slin.Remove(i--) then the for-loop i++ —
                // net i unchanged, so the element shifted into position i is
                // processed next.
                self.slin.remove(i);
            } else {
                i += 1;
            }
        }

        self.empt = self.slin.is_empty() && self.spnt.is_empty();
        self.done = true;

        if self.slin.is_empty() {
            return;
        }

        // Post processing for cones and spheres.
        let is_decompose_required = matches!(
            quad.type_quadric(),
            QuadricType::Cone | QuadricType::Sphere | QuadricType::Cylinder | QuadricType::Torus
        );
        if !is_decompose_required {
            return;
        }

        let q_surf = if reversed { s2 } else { s1 };
        let p_surf = if reversed { s1 } else { s2 };
        let mut dslin: Vec<IntPatchLine> = Vec::new();
        let mut is_decompose = false;
        for i in 1..=self.slin.len() {
            let line = self.slin[i - 1].clone();
            if decompose_result(
                &line,
                reversed,
                &quad,
                q_surf,
                p_surf,
                tol_arc,
                a_tol_3d,
                &mut dslin,
            ) {
                is_decompose = true;
            }
        }
        if is_decompose {
            self.slin = dslin;
        }
        self.empt = self.slin.is_empty() && self.spnt.is_empty();
    }
}

impl Default for ImpPrmIntersection {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT theLine->Line() — collect the line points as a Vec for Recadre.
fn thelin_points(thelin: &crate::geomalgo::int_surf::LineOn2S) -> Vec<PntOn2S> {
    (0..thelin.nb_points())
        .map(|i| thelin.value(i).clone())
        .collect()
}

/// Extract the quadric from a Surface3 (rcad: only the analytic kinds).
fn quad_of(s: &Surface3) -> Quadric {
    Quadric::from_surface3(s).unwrap_or_default()
}

/// OCCT ComputeTangency (IntPatch_ImpPrmIntersection.cxx L221-469).
fn compute_tangency(
    solrst: &SOnBounds,
    seqpdep: &mut Vec<PathPoint>,
    p_surf: &Surface3,
    func: &mut SurfFunction,
    destination: &mut [i32],
    _dom: &mut Domain,
) {
    let nb_points = solrst.nb_points();
    let mut seqlength = 0usize;
    for i in 1..=nb_points {
        if destination[i - 1] == 0 {
            let p_start = solrst.point(i).clone();
            let thearc = p_start.arc().clone();
            let theparam = p_start.parameter();
            // OCCT L256-259: arcorien + ispassing computed before the tangent check.
            let mut arcorien = _dom.orientation_arc(&thearc);
            let mut ispassing = arcorien == rcad_kernel::topods::Orientation::Internal
                || arcorien == rcad_kernel::topods::Orientation::External;

            let p2d = thearc.point_at(theparam);
            let x = [p2d.x, p2d.y];
            let mut p_point = PathPoint::new_uv(p_start.value(), x[0], x[1]);

            let _ = func.values(&x);
            if func.is_tangent() {
                p_point.set_tangency(true);
                destination[i - 1] = (seqlength + 1) as i32;
                if !p_start.is_new() {
                    let vtx = p_start.vertex();
                    let mut k = i + 1;
                    while k <= nb_points {
                        if destination[k - 1] == 0 {
                            let p_start2 = solrst.point(k).clone();
                            if !p_start2.is_new() {
                                let vtxbis = p_start2.vertex();
                                if _dom.identical(vtx, vtxbis) {
                                    let thearc2 = p_start2.arc().clone();
                                    let theparam2 = p_start2.parameter();
                                    let arcorien2 = _dom.orientation_arc(&thearc2);
                                    ispassing = ispassing
                                        && (arcorien2 == rcad_kernel::topods::Orientation::Internal
                                            || arcorien2
                                                == rcad_kernel::topods::Orientation::External);
                                    let p2d2 = thearc2.point_at(theparam2);
                                    p_point.add_uv(p2d2.x, p2d2.y);
                                    destination[k - 1] = (seqlength + 1) as i32;
                                }
                            }
                        }
                        k += 1;
                    }
                }
                p_point.set_passing(ispassing);
                seqpdep.push(p_point);
                seqlength += 1;
            } else {
                // On a un point de depart potentiel.
                let mut vectg = func.direction_3d();
                let mut dirtg = func.direction_2d();

                let (_, d1u, d1v) = p_surf.derivatives(x[0], x[1]);
                let d2d = thearc.derivative_at(theparam);
                let mut v2 = d2d.x * d1u + d2d.y * d1v;
                let v1 = d1u.cross(d1v);

                let mut test = vectg.dot(v1.cross(v2));
                if p_start.is_new() {
                    // OCCT L316-321: reverse when (test<0 AND FORWARD) OR
                    // (test>0 AND REVERSED).
                    if (test < 0.0 && arcorien == rcad_kernel::topods::Orientation::Forward)
                        || (test > 0.0 && arcorien == rcad_kernel::topods::Orientation::Reversed)
                    {
                        p_point.set_directions(-vectg, -dirtg);
                    } else {
                        p_point.set_directions(vectg, dirtg);
                    }
                    p_point.set_passing(ispassing);
                    destination[i - 1] = (seqlength + 1) as i32;
                    seqpdep.push(p_point);
                    seqlength += 1;
                } else {
                    // Traiter la transition complexe (OCCT L328-465).
                    let bidnorm = DVec3::new(1.0, 1.0, 1.0).normalize(); // gp_Dir(1., 1., 1.)
                    let tole = 1.0e-8;
                    let mut comptrans = CurveTransition::<DVec3>::new();
                    comptrans.reset_3d(vectg, bidnorm, 0.0);
                    let vtx = p_start.vertex();
                    let mut vtxorien = _dom.orientation_vertex(vtx);
                    let mut loc_trans; // OCCT L332: declared uninitialized.
                    if arcorien == rcad_kernel::topods::Orientation::Forward
                        || arcorien == rcad_kernel::topods::Orientation::Reversed
                    {
                        // Pour essai.
                        if test.abs() <= tole {
                            loc_trans = rcad_kernel::topods::Orientation::External;
                        } else {
                            if (test > 0.0
                                && arcorien == rcad_kernel::topods::Orientation::Forward)
                                || (test < 0.0
                                    && arcorien == rcad_kernel::topods::Orientation::Reversed)
                            {
                                loc_trans = rcad_kernel::topods::Orientation::Forward;
                            } else {
                                loc_trans = rcad_kernel::topods::Orientation::Reversed;
                            }
                            if arcorien == rcad_kernel::topods::Orientation::Reversed {
                                v2 = -v2; // v2.Reverse()
                            }
                        }
                        comptrans.compare(
                            tole,
                            v2.normalize(),
                            bidnorm,
                            0.0,
                            loc_trans,
                            vtxorien,
                        );
                    }
                    destination[i - 1] = (seqlength + 1) as i32;
                    let mut k = i + 1;
                    while k <= nb_points {
                        if destination[k - 1] == 0 {
                            let p_start2 = solrst.point(k).clone();
                            if !p_start2.is_new() {
                                let vtxbis = p_start2.vertex();
                                if _dom.identical(vtx, vtxbis) {
                                    let thearc2 = p_start2.arc().clone();
                                    let theparam2 = p_start2.parameter();
                                    arcorien = _dom.orientation_arc(&thearc2);
                                    p_point.add_uv(x[0], x[1]);
                                    let p2d2 = thearc2.point_at(theparam2);
                                    p_point.add_uv(p2d2.x, p2d2.y);
                                    if arcorien == rcad_kernel::topods::Orientation::Forward
                                        || arcorien
                                            == rcad_kernel::topods::Orientation::Reversed
                                    {
                                        ispassing = false;
                                        let d2d2 = thearc2.derivative_at(theparam2);
                                        v2 = d2d2.x * d1u + d2d2.y * d1v; // v2.SetLinearForm
                                        test = vectg.dot(v1.cross(v2));
                                        vtxorien = _dom.orientation_vertex(p_start2.vertex());
                                        if test.abs() <= tole {
                                            loc_trans = rcad_kernel::topods::Orientation::External;
                                        } else {
                                            if (test > 0.0
                                                && arcorien
                                                    == rcad_kernel::topods::Orientation::Forward)
                                                || (test < 0.0
                                                    && arcorien
                                                        == rcad_kernel::topods::Orientation::Reversed)
                                            {
                                                loc_trans =
                                                    rcad_kernel::topods::Orientation::Forward;
                                            } else {
                                                loc_trans =
                                                    rcad_kernel::topods::Orientation::Reversed;
                                            }
                                            if arcorien
                                                == rcad_kernel::topods::Orientation::Reversed
                                            {
                                                v2 = -v2;
                                            }
                                        }
                                        comptrans.compare(
                                            tole,
                                            v2.normalize(),
                                            bidnorm,
                                            0.0,
                                            loc_trans,
                                            vtxorien,
                                        );
                                    }
                                    destination[k - 1] = (seqlength + 1) as i32;
                                }
                            }
                        }
                        k += 1;
                    }
                    let mut fairpt = true;
                    if !ispassing {
                        let before = comptrans.state_before();
                        let after = comptrans.state_after();
                        if before == State::Unknown || after == State::Unknown {
                            fairpt = false;
                        } else if before == State::In {
                            if after == State::In {
                                ispassing = true;
                            } else {
                                vectg = -vectg;
                                dirtg = -dirtg;
                            }
                        } else if after != State::In {
                            fairpt = false;
                        }
                    }
                    if fairpt {
                        p_point.set_directions(vectg, dirtg);
                        p_point.set_passing(ispassing);
                        seqpdep.push(p_point);
                        seqlength += 1;
                    } else {
                        for k in i..=nb_points {
                            if destination[k - 1] as usize == seqlength + 1 {
                                destination[k - 1] = -destination[k - 1];
                            }
                        }
                    }
                }
            }
        }
    }
}

/// OCCT CheckSegmSegm (IntPatch_ImpPrmIntersection.cxx L3733-3753) — true when
/// the segment [theParF, theParL] is included in [theRefParF, theRefParL].
fn check_segm_segm(the_ref_par_f: f64, the_ref_par_l: f64, the_par_f: f64, the_par_l: f64) -> bool {
    if (the_par_f < the_ref_par_f) || (the_par_f > the_ref_par_l) {
        return false;
    }
    if (the_par_l < the_ref_par_f) || (the_par_l > the_ref_par_l) {
        return false;
    }
    true
}

/// OCCT ElCLib::Parameter(gp_Lin2d, gp_Pnt2d) — parameter of the orthogonal
/// projection of a point onto a unit-direction 2D line.
fn line2d_parameter(lin: &Line2d, p: DVec2) -> f64 {
    (p - lin.origin).dot(lin.direction)
}

/// OCCT gp_Lin2d::Distance(gp_Lin2d) (gp_Lin2d.hxx L228-236) — perpendicular
/// distance between two parallel lines (0 when not parallel).
fn line2d_distance(l1: &Line2d, l2: &Line2d) -> f64 {
    let a_d = l1.origin - l2.origin;
    (a_d.x * l2.direction.y - a_d.y * l2.direction.x).abs()
}

/// OCCT gp_Dir2d::IsParallel (gp_Dir2d.hxx L422-430).
fn dirs_parallel_2d(d1: DVec2, d2: DVec2, tol: f64) -> bool {
    let an_ang = (d1.x * d2.y - d1.y * d2.x).atan2(d1.x * d2.x + d1.y * d2.y).abs();
    an_ang <= tol || std::f64::consts::PI - an_ang <= tol
}

/// OCCT anArc->FirstParameter()/LastParameter() — the finite restriction-arc
/// parameter range, carried by the RLine's first/last points.
/// OCCT IsCoincide uses theArc->FirstParameter()/LastParameter() (the trimmed
/// 2D arc of the restriction line).  rcad's Curve2d::Line is untrimmed
/// (default_domain = [-inf, inf]), so the bounded range is recovered from the
/// line's first/last vertex parameters (which store the arc parameter).
fn line_arc_range(line: &IntPatchLine) -> [f64; 2] {
    if line.has_first_point() && line.has_last_point() {
        [
            line.first_point().parameter_on_line(),
            line.last_point().parameter_on_line(),
        ]
    } else if let Some(a) = line.arc_on_s1.as_ref().or(line.arc_on_s2.as_ref()) {
        a.default_domain()
    } else {
        [0.0, 0.0]
    }
}

/// OCCT IsCoincide (IntPatch_ImpPrmIntersection.cxx L3756-3892) — check if
/// `line` coincides with the 2D restriction arc `arc` (of finite range
/// `arc_range`) in 2D-space.
fn is_coincide(
    func: &mut SurfFunction,
    line: &IntPatchLine,
    arc: &Curve2d,
    arc_range: [f64; 2],
    is_surface1_using: bool,
    tol_3d: f64,
    tol_2d: f64,
    period: f64,
) -> bool {
    const A_COEFFS: [f64; 7] = [
        0.02447174185,
        0.09549150281,
        0.20610737385,
        0.34549150281,
        0.5,
        0.65450849719,
        0.79389262615,
    ];

    if line.line_type == IntPatchIType::Restriction {
        // Restriction-restriction processing.
        let Some(arc2) = line.arc_on_s1.as_ref().or(line.arc_on_s2.as_ref()) else {
            return false;
        };
        let (Curve2d::Line(lin1), Curve2d::Line(lin2)) = (arc, arc2) else {
            return false;
        };
        if !dirs_parallel_2d(lin1.direction, lin2.direction, rcad_kernel::precision::ANGULAR) {
            return false;
        }
        let a_dist = line2d_distance(lin1, lin2);
        if (a_dist < tol_2d) || ((a_dist - period).abs() < tol_2d) {
            let (a_rf, a_rl) = (arc_range[0], arc_range[1]);
            let (a_parf, a_parl) = {
                let r2 = line_arc_range(line);
                (r2[0], r2[1])
            };
            let a_p1 = lin2.point_at(a_parf);
            let a_p2 = lin2.point_at(a_parl);
            let a_param1 = line2d_parameter(lin1, a_p1);
            let a_param2 = line2d_parameter(lin1, a_p2);
            if check_segm_segm(a_rf, a_rl, a_param1, a_param2) {
                return true;
            }
            return check_segm_segm(a_param1, a_param2, a_rf, a_rl);
        }
        false
    } else {
        // Walking line vs restriction arc.
        let Curve2d::Line(an_arc_lin) = arc else {
            return false;
        };
        let (a_u_af, a_u_al) = (arc_range[0], arc_range[1]);
        for wpt in &line.wline_pnts {
            let (a_u_f, a_v_f) = if is_surface1_using {
                (wpt.u1, wpt.v1)
            } else {
                (wpt.u2, wpt.v2)
            };
            let a_ploc = DVec2::new(a_u_f, a_v_f);
            let a_r_param = line2d_parameter(an_arc_lin, a_ploc);
            if (a_r_param < a_u_af) || (a_r_param > a_u_al) {
                return false;
            }
            let a_pmin = an_arc_lin.point_at(a_r_param);
            let a_dist = a_ploc.distance(a_pmin);
            if (a_dist < tol_2d) || ((a_dist - period).abs() < tol_2d) {
                // Considered point is in Restriction line.
                continue;
            }
            // Check if intermediate points between aPloc and theArc are
            // intersection point (i.e. if aPloc is in tangent zone between
            // two intersected surfaces).
            let (a_u_l, a_v_l) = (a_pmin.x, a_pmin.y);
            let d_u = a_u_l - a_u_f;
            let d_v = a_v_l - a_v_f;
            let mut is_on_line = true;
            for i in 0..7 {
                let a_u = a_u_f + A_COEFFS[i] * d_u;
                let a_v = a_v_f + A_COEFFS[i] * d_v;
                match func.value(&[a_u, a_v]) {
                    Some(a_val) => {
                        if a_val.abs() > tol_3d {
                            is_on_line = false;
                            break;
                        }
                    }
                    None => {
                        is_on_line = false;
                        break;
                    }
                }
            }
            if !is_on_line {
                return false;
            }
        }
        true
    }
}

// ===========================================================================
// DecomposeResult (L3146-3730) — seam/pole splitting.  Ported structurally;
// the SpecialPoints refinements reuse the existing rcad-algo special_points.
// ===========================================================================

#[allow(clippy::too_many_arguments)]
fn decompose_result(
    the_line: &IntPatchLine,
    is_reversed: bool,
    the_quad: &Quadric,
    the_q_surf: &Surface3,
    the_p_surf: &Surface3,
    the_arc_tol: f64,
    the_tol_tang: f64,
    the_lines: &mut Vec<IntPatchLine>,
) -> bool {
    let a_delta_umax = std::f64::consts::FRAC_PI_2;

    let a_s_line = &the_line.wline_pnts;
    if a_s_line.len() <= 2 {
        return false;
    }

    // Deletes repeated vertices.
    let mut a_v_line = get_vertices(the_line);

    // The walking/restriction polygon, adjusted for the quadric periodicity.
    let mut a_ss_line: Vec<PntOn2S> = a_s_line
        .iter()
        .map(|w| {
            let mut p = PntOn2S::new();
            p.set_value(w.p3d, true, w.u1, w.v1);
            p.set_value_uv(false, w.u2, w.v2);
            p
        })
        .collect();
    if a_ss_line.len() <= 1 {
        return false;
    }
    adjust_line(&mut a_ss_line, is_reversed, the_q_surf);

    let mut a_l_index = a_ss_line.len();
    let mut a_f_index = 1usize;
    let mut a_b_index = 0usize;

    let mut fl_next_line = true;
    let mut has_been_decomposed = false;
    let mut a_pre_point_exist = 0i32;

    let mut pre_point = PntOn2S::new();
    while fl_next_line {
        fl_next_line = false;
        let mut is_decomposited = false;

        let mut sline: Vec<PntOn2S> = Vec::new();

        if (a_l_index <= a_f_index) && a_pre_point_exist == 0 {
            break;
        }

        if a_pre_point_exist != 0 {
            let a_ref_pt = a_ss_line[a_f_index - 1].clone();
            let (a_u_res, a_v_res) = (u_res(the_q_surf, the_arc_tol), v_res(the_q_surf, the_arc_tol));
            let a_tol_2d = if a_pre_point_exist == 1 {
                -1.0
            } else if a_pre_point_exist == 3 {
                a_v_res
            } else if a_pre_point_exist == 4 {
                a_u_res.max(a_v_res)
            } else {
                a_u_res
            };

            if continue_after_special_point(
                the_q_surf,
                the_p_surf,
                &a_ref_pt,
                a_pre_point_exist,
                a_tol_2d,
                &mut pre_point,
                is_reversed,
            ) {
                sline.push(pre_point.clone());
                while a_f_index <= a_l_index {
                    if !pre_point.is_same(&a_ss_line[a_f_index - 1], the_tol_tang, -1.0) {
                        break;
                    }
                    a_f_index += 1;
                }
            } else {
                break;
            }
        }

        a_pre_point_exist = 0;

        let mut k = a_f_index;
        while k <= a_l_index {
            if k == a_f_index {
                pre_point = a_ss_line[k - 1].clone();
                sline.push(pre_point.clone());
                k += 1;
                continue;
            }

            let mut is_on_boundary = false;
            detect_of_boundary_achievement(
                the_q_surf,
                is_reversed,
                &a_ss_line,
                k,
                &mut sline,
                &mut is_on_boundary,
            );

            a_pre_point_exist = is_seam_or_pole(the_q_surf, &a_ss_line, is_reversed, k - 1, the_tol_tang, a_delta_umax);

            if is_on_boundary && a_pre_point_exist != 5 {
                a_pre_point_exist = 0;
            }

            if a_pre_point_exist != 0 {
                a_b_index = k;
                is_decomposited = true;

                let a_ref_pt = a_ss_line[a_b_index - 2].clone();
                let mut a_new_point = a_ref_pt.clone();
                let mut a_last_type = 0i32;

                if a_pre_point_exist == 4 {
                    a_pre_point_exist = 0;
                    a_last_type = 4;
                    let mut sp_new = to_sp_pnt(&a_ref_pt);
                    crate::geomalgo::int_patch::special_points::add_cross_uv_iso_point(
                        the_q_surf,
                        the_p_surf,
                        &to_sp_pnt(&a_ref_pt),
                        the_tol_tang,
                        &mut sp_new,
                        is_reversed,
                    );
                    a_new_point = from_sp_pnt(&sp_new);
                } else if a_pre_point_exist == 3 {
                    // WLine goes through V-seam (OCCT L3332-3388).
                    a_pre_point_exist = 0;
                    a_last_type = 3;
                    let a_ref_params = a_ref_pt.parameters();
                    let (a_u0, a_v0, a_u_quad_ref, a_v_quad_ref) = if is_reversed {
                        (a_ref_params.0, a_ref_params.1, a_ref_params.2, a_ref_params.3)
                    } else {
                        (a_ref_params.2, a_ref_params.3, a_ref_params.0, a_ref_params.1)
                    };
                    let a_next_params = a_ss_line[a_b_index - 1].parameters();
                    let (a_up, a_vp, a_uq, a_vq) = if is_reversed {
                        (a_next_params.0, a_next_params.1, a_next_params.2, a_next_params.3)
                    } else {
                        (a_next_params.2, a_next_params.3, a_next_params.0, a_next_params.1)
                    };
                    let a_tol = [
                        u_res(the_p_surf, the_arc_tol),
                        v_res(the_p_surf, the_arc_tol),
                        u_res(the_q_surf, the_arc_tol),
                    ];
                    let a_start_point = [
                        0.5 * (a_u0 + a_up),
                        0.5 * (a_v0 + a_vp),
                        0.5 * (a_u_quad_ref + a_uq),
                    ];
                    let p_dom = the_p_surf.default_domain();
                    let q_dom = the_q_surf.default_domain();
                    let an_inf_bound = [p_dom[0], p_dom[2], q_dom[0]];
                    let a_sup_bound = [p_dom[1], p_dom[3], q_dom[1]];
                    let mut sp_new = to_sp_pnt(&a_ref_pt);
                    if crate::geomalgo::int_patch::special_points::add_point_on_u_or_v_iso(
                        the_q_surf,
                        the_p_surf,
                        &to_sp_pnt(&a_ref_pt),
                        false,
                        0.0,
                        &a_tol,
                        &a_start_point,
                        &an_inf_bound,
                        &a_sup_bound,
                        &mut sp_new,
                        is_reversed,
                    ) {
                        a_new_point = from_sp_pnt(&sp_new);
                    }
                } else if a_pre_point_exist == 5 {
                    a_pre_point_exist = 0;
                    let mut a_vert = crate::geomalgo::int_patch::special_points::PatchPoint::new();
                    a_vert.pnt = to_sp_pnt(&a_ref_pt);
                    a_vert.tolerance = the_tol_tang;
                    let mut sp_new = to_sp_pnt(&a_ref_pt);
                    if crate::geomalgo::int_patch::special_points::add_singular_pole(
                        the_q_surf,
                        the_p_surf,
                        &to_sp_pnt(&a_ref_pt),
                        &mut a_vert,
                        &mut sp_new,
                        is_reversed,
                    ) {
                        a_new_point = from_sp_pnt(&sp_new);
                        a_pre_point_exist = 1;
                        a_last_type = 1;
                        if is_on_boundary {
                            is_on_boundary = false;
                            sline.pop();
                        }
                    } else {
                        a_pre_point_exist = 2;
                    }
                }

                if a_pre_point_exist == 2 {
                    // WLine goes through U-seam (OCCT L3423-3479).
                    a_pre_point_exist = 0;
                    a_last_type = 2;
                    let a_ref_params = a_ref_pt.parameters();
                    let (a_u0, a_v0, a_u_quad_ref, a_v_quad_ref) = if is_reversed {
                        (a_ref_params.0, a_ref_params.1, a_ref_params.2, a_ref_params.3)
                    } else {
                        (a_ref_params.2, a_ref_params.3, a_ref_params.0, a_ref_params.1)
                    };
                    let a_next_params = a_ss_line[a_b_index - 1].parameters();
                    let (a_up, a_vp, a_uq, a_vq) = if is_reversed {
                        (a_next_params.0, a_next_params.1, a_next_params.2, a_next_params.3)
                    } else {
                        (a_next_params.2, a_next_params.3, a_next_params.0, a_next_params.1)
                    };
                    let a_tol = [
                        u_res(the_p_surf, the_arc_tol),
                        v_res(the_p_surf, the_arc_tol),
                        v_res(the_q_surf, the_arc_tol),
                    ];
                    let a_start_point = [
                        0.5 * (a_u0 + a_up),
                        0.5 * (a_v0 + a_vp),
                        0.5 * (a_v_quad_ref + a_vq),
                    ];
                    let p_dom = the_p_surf.default_domain();
                    let q_dom = the_q_surf.default_domain();
                    let an_inf_bound = [p_dom[0], p_dom[2], q_dom[2]];
                    let a_sup_bound = [p_dom[1], p_dom[3], q_dom[3]];
                    let mut sp_new = to_sp_pnt(&a_ref_pt);
                    if crate::geomalgo::int_patch::special_points::add_point_on_u_or_v_iso(
                        the_q_surf,
                        the_p_surf,
                        &to_sp_pnt(&a_ref_pt),
                        true,
                        0.0,
                        &a_tol,
                        &a_start_point,
                        &an_inf_bound,
                        &a_sup_bound,
                        &mut sp_new,
                        is_reversed,
                    ) {
                        a_new_point = from_sp_pnt(&sp_new);
                    }
                }

                if !a_new_point.is_same(
                    &a_ref_pt,
                    rcad_kernel::precision::CONFUSION,
                    rcad_kernel::precision::PCONFUSION,
                ) {
                    if is_on_boundary {
                        break;
                    }
                    sline.push(a_new_point.clone());
                    a_pre_point_exist = a_last_type;
                    pre_point = a_new_point.clone();
                } else {
                    if is_on_boundary || sline.len() == 1 {
                        pre_point = a_ref_pt.clone();
                        a_pre_point_exist = if is_on_boundary { 0 } else { a_last_type };
                    }
                }
                break;
            }

            pre_point = a_ss_line[k - 1].clone();
            if is_on_boundary {
                a_b_index = k;
                is_decomposited = true;
                a_pre_point_exist = 0;
                break;
            } else {
                sline.push(a_ss_line[k - 1].clone());
            }
            k += 1;
        }

        if sline.len() == 1 {
            fl_next_line = true;
            if a_f_index < a_b_index {
                a_f_index = a_b_index;
            }
            continue;
        }

        let mut a_v_f = PntOn2S::new();
        let mut a_v_l = PntOn2S::new();
        let mut add_v_f = false;
        let mut add_v_l = false;
        let p_dom = the_p_surf.default_domain();
        verify_vertices(
            &sline,
            is_reversed,
            &mut a_v_line,
            rcad_kernel::precision::PCONFUSION,
            the_arc_tol,
            &p_dom,
            &mut a_v_f,
            &mut add_v_f,
            &mut a_v_l,
            &mut add_v_l,
        );

        let has_internals = has_internals(&sline, &a_v_line);

        let mut d3_f = 0.0;
        let mut d3_l = 0.0;
        to_smooth(&mut sline, is_reversed, the_quad, true, &mut d3_f);
        to_smooth(&mut sline, is_reversed, the_quad, false, &mut d3_l);

        if add_v_f || add_v_l {
            let is_added = add_vertices(&mut sline, &a_v_f, add_v_f, &a_v_l, add_v_l, d3_f, d3_l);
            if is_added {
                to_smooth(&mut sline, is_reversed, the_quad, true, &mut d3_f);
                to_smooth(&mut sline, is_reversed, the_quad, false, &mut d3_l);
            }
        }

        let mut wline_pnts: Vec<crate::geomalgo::int_patch::WLinePnt> = sline
            .iter()
            .map(|p| {
                let (u1, v1, u2, v2) = p.parameters();
                crate::geomalgo::int_patch::WLinePnt {
                    p3d: p.value(),
                    u1,
                    v1,
                    u2,
                    v2,
                }
            })
            .collect();

        if the_line.line_type == IntPatchIType::Walking {
            let mut wline = IntPatchLine::walking(wline_pnts, crate::geomalgo::int_patch::WLineType::ImpPrm);
            let mut a_tpnt_f = IntPatchVertex::default();
            let a_s_pnt = wline.point(0).p3d;
            a_tpnt_f.set_value(a_s_pnt, the_arc_tol, false);
            let (u1, v1, u2, v2) = (wline.point(0).u1, wline.point(0).v1, wline.point(0).u2, wline.point(0).v2);
            a_tpnt_f.set_parameters(u1, v1, u2, v2);
            a_tpnt_f.set_parameter(1.0);
            wline.add_vertex(a_tpnt_f);
            wline.set_first_point(1);

            if has_internals {
                put_int_vertices(&mut wline, &sline, is_reversed, &a_v_line, the_arc_tol);
            }

            let n = wline.nb_points();
            let mut a_tpnt_l = IntPatchVertex::default();
            let a_s_pnt = wline.point(n - 1).p3d;
            a_tpnt_l.set_value(a_s_pnt, the_arc_tol, false);
            let (u1, v1, u2, v2) = (wline.point(n - 1).u1, wline.point(n - 1).u2, wline.point(n - 1).v2, wline.point(n - 1).u2);
            a_tpnt_l.set_parameters(u1, v1, u2, v2);
            a_tpnt_l.set_parameter(n as f64);
            wline.add_vertex(a_tpnt_l);
            wline.set_last_point(wline.nb_vertex());

            the_lines.push(wline);
        } else {
            // Restriction line.
            if !is_decomposited && !has_been_decomposed {
                the_lines.push(the_line.clone());
                return has_been_decomposed;
            }
            let mut a_r_line = the_line.clone();
            a_r_line.vertices.clear();
            a_r_line.wline_pnts = wline_pnts;
            if has_internals {
                put_int_vertices(&mut a_r_line, &sline, is_reversed, &a_v_line, the_arc_tol);
            }
            the_lines.push(a_r_line);
        }

        if is_decomposited {
            a_f_index = a_b_index;
            fl_next_line = true;
            has_been_decomposed = true;
        }
    }

    has_been_decomposed
}

/// OCCT GeomAdaptor_Surface::UResolution (GeomAdaptor_Surface.cxx L1818-1896).
/// Parametric resolution of the surface in the U direction for a 3D tolerance.
fn u_res(s: &Surface3, tol3d: f64) -> f64 {
    match s {
        Surface3::Plane(_) => tol3d,
        Surface3::Cylinder(c) => resolution_quadric_angular(tol3d, c.radius),
        Surface3::Sphere(sp) => resolution_quadric_angular(tol3d, sp.radius),
        Surface3::Torus(t) => resolution_quadric_angular(tol3d, t.major_radius + t.minor_radius),
        Surface3::Cone(co) => {
            let d = s.default_domain();
            if d[3] - d[2] > 1.0e10 {
                // Not truly bounded => unknown resolution.
                rcad_kernel::precision::parametric_default(tol3d)
            } else {
                // OCCT: R = max radius of the VIso circles at the V bounds.
                let r1 = co.radius_at_slant(d[2]);
                let r2 = co.radius_at_slant(d[3]);
                let r = if r1 > r2 { r1 } else { r2 };
                if r > rcad_kernel::precision::CONFUSION {
                    tol3d / r
                } else {
                    0.0
                }
            }
        }
        Surface3::LinearExtrusion(le) => curve_resolution(&le.profile, tol3d),
        _ => rcad_kernel::precision::parametric_default(tol3d),
    }
}

/// OCCT GeomAdaptor_Surface::VResolution (GeomAdaptor_Surface.cxx L1900-1958).
fn v_res(s: &Surface3, tol3d: f64) -> f64 {
    match s {
        Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::LinearExtrusion(_) => {
            tol3d
        }
        Surface3::Sphere(sp) => resolution_quadric_angular(tol3d, sp.radius),
        Surface3::Torus(t) => resolution_quadric_angular(tol3d, t.minor_radius),
        Surface3::Revolution(rv) => curve_resolution(&rv.profile, tol3d),
        _ => rcad_kernel::precision::parametric_default(tol3d),
    }
}

/// OCCT L1890-1895 / L1952-1957: Res = R3d/(2R) -> 2*asin(Res) when Res <= 1,
/// else the period 2*PI; R <= Confusion -> 0.
fn resolution_quadric_angular(tol3d: f64, radius: f64) -> f64 {
    if radius > rcad_kernel::precision::CONFUSION {
        let res = tol3d / (2.0 * radius);
        if res <= 1.0 {
            2.0 * res.asin()
        } else {
            std::f64::consts::TAU
        }
    } else {
        0.0
    }
}

/// OCCT GeomAdaptor_Curve::Resolution (GeomAdaptor_Curve.cxx L1116-1149) —
/// parametric resolution of a curve for a 3D tolerance.
fn curve_resolution(c: &Curve3, tol3d: f64) -> f64 {
    match c {
        Curve3::Line(_) => tol3d,
        Curve3::Circle(cr) => {
            let r = cr.radius;
            if r > tol3d / 2.0 {
                2.0 * (tol3d / (2.0 * r)).asin()
            } else {
                std::f64::consts::TAU
            }
        }
        Curve3::Ellipse(e) => tol3d / e.major_radius,
        _ => rcad_kernel::precision::parametric_default(tol3d),
    }
}

/// OCCT AdjustU (L2172-2193).
fn adjust_u(mut u: f64) -> f64 {
    let dblpi = TAU;
    if u < 0.0 || u > dblpi {
        if u < 0.0 {
            while u < 0.0 {
                u += dblpi;
            }
        } else {
            while u > dblpi {
                u -= dblpi;
            }
        }
    }
    u
}

/// OCCT GetVertices (L2057-2142) — collect line vertices, rejecting equals.
fn get_vertices(the_p_line: &IntPatchLine) -> Vec<PntOn2S> {
    let tol_3d = 1e-10;
    let tol_2d = rcad_kernel::precision::PCONFUSION;
    let mut vertices: Vec<PntOn2S> = Vec::new();
    let nb_vrt = the_p_line.nb_vertex();
    let mut an_vrts = vec![0i32; nb_vrt];
    for i in 1..=nb_vrt {
        if an_vrts[i - 1] == -1 {
            continue;
        }
        let pi = the_p_line.vertex(i).clone();
        let mut k = i + 1;
        while k <= nb_vrt {
            if an_vrts[k - 1] == -1 {
                k += 1;
                continue;
            }
            let pk = the_p_line.vertex(k).clone();
            if pi.p3d.distance(pk.p3d) <= tol_3d {
                let (u1, v1) = pi.parameters_on_s1();
                let (u2, v2) = pk.parameters_on_s1();
                let same_u1 = (u1 - u2).abs() <= tol_2d;
                let same_v1 = (v1 - v2).abs() <= tol_2d;
                let (u1, v1) = pi.parameters_on_s2();
                let (u2, v2) = pk.parameters_on_s2();
                let same_u2 = (u1 - u2).abs() <= tol_2d;
                let same_v2 = (v1 - v2).abs() <= tol_2d;
                if (same_u1 && same_v1) && (same_u2 && same_v2) {
                    an_vrts[k - 1] = -1;
                }
            }
            k += 1;
        }
    }
    for i in 1..=nb_vrt {
        if an_vrts[i - 1] == -1 {
            continue;
        }
        let v = the_p_line.vertex(i);
        let mut p = PntOn2S::new();
        p.set_value(v.p3d, true, v.u1, v.v1);
        p.set_value_uv(false, v.u2, v.v2);
        vertices.push(p);
    }
    vertices
}

/// OCCT AdjustLine (L2226-2255).
fn adjust_line(line: &mut Vec<PntOn2S>, is_reversed: bool, q_surf: &Surface3) {
    let d = q_surf.default_domain();
    let (uf, ul, vf, vl) = (d[0], d[1], d[2], d[3]);
    let nbp = line.len();
    for ip in 0..nbp {
        if is_reversed {
            let (u, v) = line[ip].parameters_on_surface(false);
            let u = adjust_u(u);
            let (u, v) = correct_2d_bounds(uf, ul, vf, vl, rcad_kernel::precision::PCONFUSION, u, v);
            line[ip].set_value_uv(false, u, v);
        } else {
            let (u, v) = line[ip].parameters_on_surface(true);
            let u = adjust_u(u);
            let (u, v) = correct_2d_bounds(uf, ul, vf, vl, rcad_kernel::precision::PCONFUSION, u, v);
            line[ip].set_value_uv(true, u, v);
        }
    }
}

/// OCCT Correct2DBounds (L2195-2224).
fn correct_2d_bounds(uf: f64, ul: f64, vf: f64, vl: f64, tol_2d: f64, mut u: f64, mut v: f64) -> (f64, f64) {
    let eps = 1e-16;
    let d_uf = (u - uf).abs();
    let d_ul = (u - ul).abs();
    let d_vf = (v - vf).abs();
    let d_vl = (v - vl).abs();
    if d_uf <= tol_2d && d_uf > eps {
        u = uf;
    }
    if d_ul <= tol_2d && d_ul > eps {
        u = ul;
    }
    if d_vf <= tol_2d && d_vf > eps {
        v = vf;
    }
    if d_vl <= tol_2d && d_vl > eps {
        v = vl;
    }
    (u, v)
}

/// OCCT IsSeamParameter (L2167-2170).
#[allow(dead_code)]
fn is_seam_parameter(u: f64, tol_2d: f64) -> bool {
    u.abs() <= tol_2d || (TAU - u).abs() <= tol_2d
}

/// OCCT IsPointOnBoundary (IntPatch_ImpPrmIntersection.cxx L3043-3060) — TRUE
/// when theParam matches theBoundary +/- thePeriod within theToler2D.
fn is_point_on_boundary(tol_2d: f64, boundary: f64, period: f64, param: f64) -> bool {
    let mut a_delta = (param - boundary).abs();
    if period != 0.0 {
        a_delta = a_delta % period;
        return a_delta < tol_2d || (period - a_delta) < tol_2d;
    }
    a_delta < tol_2d
}

/// OCCT DetectOfBoundaryAchievement (L3067-3137): the WLine reaches the quadric
/// domain boundary when the current point is on it and the previous is not;
/// the boundary point is then adjusted to avoid U "jumping" across the seam.
fn detect_of_boundary_achievement(
    q_surf: &Surface3,
    is_reversed: bool,
    line: &[PntOn2S],
    k: usize,
    sline: &mut Vec<PntOn2S>,
    is_on_boundary: &mut bool,
) {
    let a_u_period = if q_surf.is_u_periodic() { TAU } else { 0.0 };
    let a_v_period = if q_surf.is_v_periodic() { TAU } else { 0.0 };
    let d = q_surf.default_domain();
    let (a_uf, a_ul, a_vf, a_vl) = (d[0], d[1], d[2], d[3]);
    let a_p_prev = &line[k - 2];
    let a_p_curr = &line[k - 1];
    let (mut a_u_prev, a_v_prev, mut a_u_curr, mut a_v_curr) = if is_reversed {
        let (u1, v1) = a_p_prev.parameters_on_surface(false);
        let (u2, v2) = a_p_curr.parameters_on_surface(false);
        (u1, v1, u2, v2)
    } else {
        let (u1, v1) = a_p_prev.parameters_on_surface(true);
        let (u2, v2) = a_p_curr.parameters_on_surface(true);
        (u1, v1, u2, v2)
    };
    let tol = rcad_kernel::precision::PCONFUSION;
    if is_point_on_boundary(tol, a_uf, a_u_period, a_u_curr)
        && !is_point_on_boundary(tol, a_uf, a_u_period, a_u_prev)
    {
        *is_on_boundary = true;
    } else if is_point_on_boundary(tol, a_ul, a_u_period, a_u_curr)
        && !is_point_on_boundary(tol, a_ul, a_u_period, a_u_prev)
    {
        *is_on_boundary = true;
    } else if is_point_on_boundary(tol, a_vf, a_v_period, a_v_curr)
        && !is_point_on_boundary(tol, a_vf, a_v_period, a_v_prev)
    {
        *is_on_boundary = true;
    } else if is_point_on_boundary(tol, a_vl, a_v_period, a_v_curr)
        && !is_point_on_boundary(tol, a_vl, a_v_period, a_v_prev)
    {
        *is_on_boundary = true;
    }
    if *is_on_boundary {
        // Adjust, to avoid bad jumping of the WLine.
        let a_du = a_u_prev - a_u_curr;
        let a_dv = a_v_prev - a_v_curr;
        if a_u_period > 0.0 && (2.0 * a_du.abs() > a_u_period) {
            a_u_curr += a_u_period.copysign(a_du);
        }
        if a_v_period > 0.0 && (2.0 * a_dv.abs() > a_v_period) {
            a_v_curr += a_v_period.copysign(a_dv);
        }
        let mut a_point = a_p_curr.clone();
        a_point.set_value_uv(!is_reversed, a_u_curr, a_v_curr);
        sline.push(a_point);
    }
}

/// OCCT TestMiddleOnPrm (IntPatch_ImpPrmIntersection.cxx L2527-2554): the
/// midpoint between aP and aV (on the parametric surface) must be IN or ON the
/// parametric domain.
fn test_middle_on_prm(a_p: &PntOn2S, a_v: &PntOn2S, is_reversed: bool, arc_tol: f64, p_domain: &[f64; 4]) -> bool {
    let (up, vp) = if is_reversed {
        a_p.parameters_on_surface(true)
    } else {
        a_p.parameters_on_surface(false)
    };
    let (uv, vv) = if is_reversed {
        a_v.parameters_on_surface(true)
    } else {
        a_v.parameters_on_surface(false)
    };
    let um = 0.5 * (up + uv);
    let vm = 0.5 * (vp + vv);
    // rcad: the parametric domain is a rectangle; Classify returns ON when the
    // point is within ArcTol of a bound, IN when strictly inside.
    classify_rect(um, vm, arc_tol, p_domain) != 3
}

/// rcad Adaptor3d_TopolTool::Classify for a rectangular domain:
/// 1 = IN, 2 = ON, 3 = OUT.
fn classify_rect(u: f64, v: f64, tol: f64, d: &[f64; 4]) -> i32 {
    if u < d[0] - tol || u > d[1] + tol || v < d[2] - tol || v > d[3] + tol {
        return 3;
    }
    if (u - d[0]).abs() <= tol || (u - d[1]).abs() <= tol || (v - d[2]).abs() <= tol || (v - d[3]).abs() <= tol {
        return 2;
    }
    1
}

/// OCCT ToSmooth (IntPatch_ImpPrmIntersection.cxx L2399-2525): move the
/// first/last point of the line off the seam by a small step so that the
/// line's U-parameters do not "jump" across the seam.
fn to_smooth(line: &mut [PntOn2S], is_reversed: bool, the_quad: &Quadric, is_first: bool, d3d: &mut f64) {
    let nbp = line.len();
    if nbp <= 10 {
        return;
    }
    *d3d = 0.0;
    let mut nb_test_pnts = nbp / 5;
    if nb_test_pnts < 5 {
        nb_test_pnts = 5;
    }
    let startp = if is_first { 2 } else { nbp - nb_test_pnts - 2 };
    let mut ddu = 0.0;
    for ip in startp..=nb_test_pnts {
        let (uc, vc) = if is_reversed {
            line[ip - 1].parameters_on_surface(false)
        } else {
            line[ip - 1].parameters_on_surface(true)
        };
        let (un, vn) = if is_reversed {
            line[ip].parameters_on_surface(false)
        } else {
            line[ip].parameters_on_surface(true)
        };
        ddu += (uc.abs() - un.abs()).abs();
        if ip > startp {
            *d3d += line[ip - 1].value().distance(line[ip - 2].value());
        }
    }
    ddu /= (nb_test_pnts + 1) as f64;
    *d3d /= (nb_test_pnts + 1) as f64;

    let (index1, index2, index3) = if is_first {
        (0usize, 1usize, 2usize)
    } else {
        (nbp - 1, nbp - 2, nbp - 3)
    };
    let (u1, v1, u2, v2, u3, v3) = if is_reversed {
        let (a, b) = line[index1].parameters_on_surface(false);
        let (c, d) = line[index2].parameters_on_surface(false);
        let (e, f) = line[index3].parameters_on_surface(false);
        (a, b, c, d, e, f)
    } else {
        let (a, b) = line[index1].parameters_on_surface(true);
        let (c, d) = line[index2].parameters_on_surface(true);
        let (e, f) = line[index3].parameters_on_surface(true);
        (a, b, c, d, e, f)
    };
    let mut do_u = false;
    let _ = v1;
    let _ = v2;
    let _ = v3;
    if the_quad.type_quadric() == QuadricType::Sphere {
        if (u1.abs() - u2.abs()).abs() > (std::f64::consts::PI / 16.0) {
            do_u = true;
        }
        if do_u && (u1.abs() <= 1e-9 || (u1 - TAU).abs() <= 1e-9) {
            if (v1 - std::f64::consts::FRAC_PI_2).abs() <= 1e-9 || (v1 + std::f64::consts::FRAC_PI_2).abs() <= 1e-9 {
            } else {
                do_u = false;
            }
        }
    }
    if the_quad.type_quadric() == QuadricType::Cone {
        let (u_apx, v_apx) = the_quad.parameters(the_quad.cone().apex_point());
        if (u1.abs() - u2.abs()).abs() > (std::f64::consts::PI / 32.0) {
            do_u = true;
        }
        if do_u && (u1.abs() <= 1e-9 || (u1 - TAU).abs() <= 1e-9) {
            if (v1 - v_apx).abs() <= 1e-9 {
            } else {
                do_u = false;
            }
        }
        let _ = (u_apx, v_apx);
    }
    if do_u {
        let d_u = (ddu / 10.0).min(5e-8);
        let u = if u2 > u3 { u2 + d_u } else { u2 - d_u };
        if is_reversed {
            line[index1].set_value_uv(false, u, v1);
        } else {
            line[index1].set_value_uv(true, u, v1);
        }
    }
}

/// OCCT AddVertices (L2826-2854): insert VrtF/VrtL before/after the line when
/// the 3D distance condition is satisfied.
fn add_vertices(
    line: &mut Vec<PntOn2S>,
    vrt_f: &PntOn2S,
    add_first: bool,
    vrt_l: &PntOn2S,
    add_last: bool,
    d3d_f: f64,
    d3d_l: f64,
) -> bool {
    let mut result = false;
    if add_first {
        let df = line[0].value().distance(vrt_f.value());
        if (d3d_f * 2.0) > df && df > 1.5e-7 {
            line.insert(0, vrt_f.clone());
            result = true;
        }
    }
    if add_last {
        let dl = line[line.len() - 1].value().distance(vrt_l.value());
        if (d3d_l * 2.0) > dl && dl > 1.5e-7 {
            line.push(vrt_l.clone());
            result = true;
        }
    }
    result
}

/// OCCT VerifyVertices (L2556-2824).
#[allow(clippy::too_many_arguments)]
fn verify_vertices(
    line: &[PntOn2S],
    is_reversed: bool,
    vertices: &mut [PntOn2S],
    tol_2d: f64,
    arc_tol: f64,
    p_domain: &[f64; 4],
    vrt_f: &mut PntOn2S,
    add_first: &mut bool,
    vrt_l: &mut PntOn2S,
    add_last: &mut bool,
) {
    let nbp = line.len();
    let nbv = vertices.len();
    let mut f_index_same = 0usize;
    let mut f_index_near = 0usize;
    let mut l_index_same = 0usize;
    let mut l_index_near = 0usize;
    let a_pf = &line[0];
    let a_pl = &line[nbp - 1];
    let (uf, vf) = if is_reversed {
        a_pf.parameters_on_surface(false)
    } else {
        a_pf.parameters_on_surface(true)
    };
    let (ul, vl) = if is_reversed {
        a_pl.parameters_on_surface(false)
    } else {
        a_pl.parameters_on_surface(true)
    };
    let a2d_pf = DVec2::new(uf, vf);
    let a2d_pl = DVec2::new(ul, vl);
    let mut dist_min_f = 1e+100;
    let mut dist_min_l = 1e+100;
    let mut f_conjugated = 0usize;
    let mut l_conjugated = 0usize;

    // AdjustU each vertex U and store it back.
    for iv in 0..nbv {
        let (uv, vv) = if is_reversed {
            vertices[iv].parameters_on_surface(false)
        } else {
            vertices[iv].parameters_on_surface(true)
        };
        let uv = adjust_u(uv);
        if is_reversed {
            vertices[iv].set_value_uv(false, uv, vv);
        } else {
            vertices[iv].set_value_uv(true, uv, vv);
        }
    }

    // Find the vertex matching the first point.
    for iv in 0..nbv {
        let a_v = &vertices[iv];
        if a_pf.is_same(a_v, rcad_kernel::precision::CONFUSION, rcad_kernel::precision::PCONFUSION) {
            f_index_same = iv + 1;
            break;
        } else {
            let (uv, vv) = if is_reversed {
                a_v.parameters_on_surface(false)
            } else {
                a_v.parameters_on_surface(true)
            };
            let a2d_v = DVec2::new(uv, vv);
            let dist = a2d_v.distance(a2d_pf);
            if dist < dist_min_f {
                dist_min_f = dist;
                f_index_near = iv + 1;
                if f_conjugated != 0 {
                    f_conjugated = 0;
                }
            }
            if is_seam_parameter(uv, tol_2d) {
                let ucv = if uv.abs() < (TAU - uv).abs() { TAU } else { 0.0 };
                let a2d_cv = DVec2::new(ucv, vv);
                let c_dist = a2d_cv.distance(a2d_pf);
                if c_dist < dist_min_f {
                    dist_min_f = c_dist;
                    f_conjugated = iv + 1;
                    f_index_near = iv + 1;
                }
            }
        }
    }

    // Find the vertex matching the last point.
    for iv in 0..nbv {
        let a_v = &vertices[iv];
        if a_pl.is_same(a_v, rcad_kernel::precision::CONFUSION, rcad_kernel::precision::PCONFUSION) {
            l_index_same = iv + 1;
            break;
        } else {
            let (uv, vv) = if is_reversed {
                a_v.parameters_on_surface(false)
            } else {
                a_v.parameters_on_surface(true)
            };
            let a2d_v = DVec2::new(uv, vv);
            let dist = a2d_v.distance(a2d_pl);
            if dist < dist_min_l {
                dist_min_l = dist;
                l_index_near = iv + 1;
                if l_conjugated != 0 {
                    l_conjugated = 0;
                }
            }
            if is_seam_parameter(uv, tol_2d) {
                let ucv = if uv.abs() < (TAU - uv).abs() { TAU } else { 0.0 };
                let a2d_cv = DVec2::new(ucv, vv);
                let c_dist = a2d_cv.distance(a2d_pl);
                if c_dist < dist_min_l {
                    dist_min_l = c_dist;
                    l_conjugated = iv + 1;
                    l_index_near = iv + 1;
                }
            }
        }
    }

    *add_first = false;
    *add_last = false;

    if f_index_same == 0 {
        if f_index_near != 0 {
            let a_v = &vertices[f_index_near - 1];
            let (uv, vv) = if is_reversed {
                a_v.parameters_on_surface(false)
            } else {
                a_v.parameters_on_surface(true)
            };
            if is_seam_parameter(uv, tol_2d) {
                let ucv = if uv.abs() < (TAU - uv).abs() { TAU } else { 0.0 };
                let test = test_middle_on_prm(a_pf, a_v, is_reversed, arc_tol, p_domain);
                if test {
                    vrt_f.set_value_pt(a_v.value());
                    if is_reversed {
                        let (u2, v2) = a_v.parameters_on_surface(true);
                        vrt_f.set_value_uv(true, u2, v2);
                        if f_conjugated == 0 {
                            vrt_f.set_value_uv(false, uv, vv);
                        } else {
                            vrt_f.set_value_uv(false, ucv, vv);
                        }
                    } else {
                        let (u2, v2) = a_v.parameters_on_surface(false);
                        vrt_f.set_value_uv(false, u2, v2);
                        if f_conjugated == 0 {
                            vrt_f.set_value_uv(true, uv, vv);
                        } else {
                            vrt_f.set_value_uv(true, ucv, vv);
                        }
                    }
                    let dist_3d = vrt_f.value().distance(a_pf.value());
                    if dist_3d > 1.5e-7 && dist_min_f > tol_2d {
                        *add_first = true;
                    }
                }
            }
        }
    }

    if l_index_same == 0 {
        if l_index_near != 0 {
            let a_v = &vertices[l_index_near - 1];
            let (uv, vv) = if is_reversed {
                a_v.parameters_on_surface(false)
            } else {
                a_v.parameters_on_surface(true)
            };
            if is_seam_parameter(uv, tol_2d) {
                let ucv = if uv.abs() < (TAU - uv).abs() { TAU } else { 0.0 };
                let test = test_middle_on_prm(a_pl, a_v, is_reversed, arc_tol, p_domain);
                if test {
                    vrt_l.set_value_pt(a_v.value());
                    if is_reversed {
                        let (u2, v2) = a_v.parameters_on_surface(true);
                        vrt_l.set_value_uv(true, u2, v2);
                        if l_conjugated == 0 {
                            vrt_l.set_value_uv(false, uv, vv);
                        } else {
                            vrt_l.set_value_uv(false, ucv, vv);
                        }
                    } else {
                        let (u2, v2) = a_v.parameters_on_surface(false);
                        vrt_l.set_value_uv(false, u2, v2);
                        if l_conjugated == 0 {
                            vrt_l.set_value_uv(true, uv, vv);
                        } else {
                            vrt_l.set_value_uv(true, ucv, vv);
                        }
                    }
                    let dist_3d = vrt_l.value().distance(a_pl.value());
                    if dist_3d > 1.5e-7 && dist_min_l > tol_2d {
                        *add_last = true;
                    }
                }
            }
        }
    }
}

/// OCCT HasInternals — true when any line point coincides with a vertex.
fn has_internals(line: &[PntOn2S], vertices: &[PntOn2S]) -> bool {
    for a_p in line.iter().skip(1).take(line.len().saturating_sub(2)) {
        for a_v in vertices {
            if a_p.is_same(a_v, rcad_kernel::precision::CONFUSION, rcad_kernel::precision::PCONFUSION) {
                return true;
            }
        }
    }
    false
}

/// OCCT PutIntVertices (IntPatch_ImpPrmIntersection.cxx L2856-2920): add the
/// vertices that coincide with the internal line points (indices 2..nbp-1) to
/// the line, with the parameter from the line index (or the RLine arc
/// parameter when the line is a Restriction line).
fn put_int_vertices(
    line: &mut IntPatchLine,
    result: &[PntOn2S],
    is_reversed: bool,
    vertices: &[PntOn2S],
    the_arc_tol: f64,
) {
    let nbp = result.len();
    if nbp < 3 {
        return;
    }
    let nbv = vertices.len();
    let is_rline = line.line_type == IntPatchIType::Restriction;
    for ip in 1..nbp - 1 {
        let a_p = &result[ip];
        for iv in 0..nbv {
            let a_v = &vertices[iv];
            if a_p.is_same(a_v, rcad_kernel::precision::CONFUSION, rcad_kernel::precision::PCONFUSION) {
                let a_pnt = result[ip].value();
                let (u1, v1) = result[ip].parameters_on_surface(true);
                let (u2, v2) = result[ip].parameters_on_surface(false);
                let mut the_pnt = IntPatchVertex::default();
                the_pnt.set_value(a_pnt, the_arc_tol, false);
                the_pnt.set_parameters(u1, v1, u2, v2);
                let mut a_param = (ip + 1) as f64;
                if is_rline {
                    let an_arc = line.arc_on_s1.as_ref().or(line.arc_on_s2.as_ref());
                    if let Some(Curve2d::Line(a_lin)) = an_arc {
                        let a_p_surf = if is_reversed {
                            DVec2::new(u1, v1)
                        } else {
                            DVec2::new(u2, v2)
                        };
                        a_param = line2d_parameter(a_lin, a_p_surf);
                    }
                }
                the_pnt.set_parameter(a_param);
                line.add_vertex(the_pnt);
            }
        }
    }
}

/// Convert an int_surf::PntOn2S to the special_points::PntOn2S representation.
fn to_sp_pnt(p: &PntOn2S) -> crate::geomalgo::int_patch::special_points::PntOn2S {
    let (u1, v1, u2, v2) = p.parameters();
    crate::geomalgo::int_patch::special_points::PntOn2S {
        p: p.value(),
        u1,
        v1,
        u2,
        v2,
    }
}

/// Convert a special_points::PntOn2S back to the int_surf::PntOn2S.
fn from_sp_pnt(p: &crate::geomalgo::int_patch::special_points::PntOn2S) -> PntOn2S {
    let mut r = PntOn2S::new();
    r.set_value(p.p, true, p.u1, p.v1);
    r.set_value_uv(false, p.u2, p.v2);
    r
}

/// OCCT IntPatch_SpecialPoints::ContinueAfterSpecialPoint (IntPatch_SpecialPoints.cxx
/// L990-1078): if the last point of the line is the pole/seam of the quadric
/// then the line was broken there; the new line must start from this point with
/// recomputed 2D-coordinates (U of the quadric changes by +/-PI at a pole).
#[allow(clippy::too_many_arguments)]
fn continue_after_special_point(
    q_surf: &Surface3,
    p_surf: &Surface3,
    ref_pt: &PntOn2S,
    pre_point_type: i32,
    tol_2d: f64,
    pre_point: &mut PntOn2S,
    is_reversed: bool,
) -> bool {
    if pre_point_type == 0 {
        return false;
    }
    if pre_point.is_same(ref_pt, rcad_kernel::precision::CONFUSION, tol_2d) {
        return false;
    }
    if (pre_point_type == 1) && matches!(q_surf, Surface3::Cone(_)) {
        // Check if the condition b) is satisfied.  Repeat the same steps as in
        // IntPatch_SpecialPoints::AddSingularPole(...) (L1009-1045).
        let (a_u0, a_v0, a_u_quad, a_v_quad) = if is_reversed {
            let (u1, v1, u2, v2) = pre_point.parameters();
            (u1, v1, u2, v2)
        } else {
            let (u1, v1, u2, v2) = pre_point.parameters();
            (u2, v2, u1, v1)
        };
        let (_, du, dv) = p_surf.derivatives(a_u0, a_v0);
        // Transforms the parametric surface in the coordinate-system of the cone.
        let cone = match q_surf {
            Surface3::Cone(c) => *c,
            _ => unreachable!(),
        };
        let frame = crate::geomalgo::int_patch::special_points::quadric_frame(q_surf);
        let du_t = crate::geomalgo::int_patch::special_points::transform_vec(du, frame);
        let dv_t = crate::geomalgo::int_patch::special_points::transform_vec(dv, frame);
        let mut u_quad = a_u_quad;
        let mut is_iso_chosen = false;
        crate::geomalgo::int_patch::special_points::process_cone(
            &to_sp_pnt(pre_point),
            du_t,
            dv_t,
            &cone,
            is_reversed,
            &mut u_quad,
            &mut is_iso_chosen,
        );
        // OCCT L1044: theNewPoint.SetValue(!theIsReversed, aUquad, aVquad).
        pre_point.set_value_uv(!is_reversed, u_quad, a_v_quad);
    }
    // Periods: PI/2 for a pole (to avoid "jumping" between neighbor points),
    // 2PI for seams.
    let a_period = if pre_point_type == 1 {
        std::f64::consts::FRAC_PI_2
    } else {
        std::f64::consts::TAU
    };
    let a_up_period = if p_surf.is_u_periodic() { std::f64::consts::TAU } else { 0.0 };
    let a_uq_period = if q_surf.is_u_periodic() { a_period } else { 0.0 };
    let a_vp_period = if p_surf.is_v_periodic() { std::f64::consts::TAU } else { 0.0 };
    let a_vq_period = if q_surf.is_v_periodic() { a_period } else { 0.0 };
    let an_arr_of_period = [
        if is_reversed { a_up_period } else { a_uq_period },
        if is_reversed { a_vp_period } else { a_vq_period },
        if is_reversed { a_uq_period } else { a_up_period },
        if is_reversed { a_vq_period } else { a_vp_period },
    ];
    let mut sp_new = to_sp_pnt(pre_point);
    crate::geomalgo::int_patch::special_points::adjust_point_and_vertex(
        &to_sp_pnt(ref_pt),
        &an_arr_of_period,
        &mut sp_new,
        None,
    );
    *pre_point = from_sp_pnt(&sp_new);
    true
}
