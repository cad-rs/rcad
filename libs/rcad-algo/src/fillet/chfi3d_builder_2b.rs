//! OCCT ChFi3d_Builder_2.cxx — 1:1 translation, second half (Stage 1f).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_Builder_2.cxx, L2205-3900.
//!
//! Coverage: PerformSetOfSurfOnElSpine (L2207-3000), ChFi3d_BoxDiag
//! (L3282-3294), PerformSetOfKGen (L3298-3878) and the GeomFill form
//! carriers the KGen reprocessing references.
//!
//! Pre-existing translations in `chfi3d.rs` (name collisions on
//! `impl ChFi3dBuilder`, not redefined here):
//!   - PerformSetOfKPart (OCCT L3004-3280) -> perform_set_of_k_part
//!   - PerformSetOfSurf   (OCCT L3882-3900) -> perform_set_of_surf

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema::extrema_locate_ext_pc;
use rcad_kernel::geom::{Curve2dEval as _, Surface3};
use rcad_kernel::topo::topods::Shape;

use super::chfi3d::next_side;
use super::chfi3d_builder_0::{chfi3d_enlarge_box_dstr, chfi3d_enlarge_box_edge_faces, chfi3d_reparam_pcurv, BndBox, BRepAdaptorSurface};
use super::chfi3d_builder_2::{chfi3d_build_plane, BRepAdaptorCurve2d, BRepTopAdaptorTopolTool};
use super::chfi3d_builder_2c::{chfi3d_purge, insert_after, insert_before, remove_sd};
use super::chfi_ds::{
    elclib_in_period, ChFiDS_ChamfMode, ChFiDS_ErrorStatus, ChFiDS_State, ChFiDSSurfData,
    SharedStripe, SharedSurfData,
};

// =========================================================================
// OCCT GeomFill_Boundary (TKGeomAlgo/GeomFill) — pending translation; the
// empty carrier keeps the mkbound call sites of PerformSetOfKGen typed.
// =========================================================================
#[derive(Debug, Clone, Default)]
pub struct GeomFillBoundary;

// =========================================================================
// OCCT GeomFill_ConstrainedFilling (TKGeomAlgo/GeomFill) — pending
// translation.  Init records the boundary configuration; Surface() reports
// None until the filling machinery lands (the OCCT flow always produces a
// surface — see the PerformSetOfKGen skip comment).
// =========================================================================
#[derive(Debug, Clone)]
pub struct GeomFillConstrainedFilling {
    pub deg: i32,
    pub maxseg: i32,
}

impl GeomFillConstrainedFilling {
    /// OCCT GeomFill_ConstrainedFilling(Deg, MaxSeg).
    pub fn new(deg: i32, maxseg: i32) -> Self {
        GeomFillConstrainedFilling { deg, maxseg }
    }

    /// OCCT GeomFill_ConstrainedFilling::Init(B1, B2, B3, Guess).
    pub fn init3(&mut self, _b1: &GeomFillBoundary, _b2: &GeomFillBoundary, _b3: &GeomFillBoundary, _guess: bool) {
    }

    /// OCCT GeomFill_ConstrainedFilling::Init(B1, B2, B3, B4, Guess).
    pub fn init4(
        &mut self,
        _b1: &GeomFillBoundary,
        _b2: &GeomFillBoundary,
        _b3: &GeomFillBoundary,
        _b4: &GeomFillBoundary,
        _guess: bool,
    ) {
    }

    /// OCCT GeomFill_ConstrainedFilling::Surface() — pending.
    pub fn surface(&self) -> Option<Surface3> {
        None
    }
}

// =========================================================================
// Pending-leaf stand-ins of ChFi3d_Builder::SimulSurf / PerformSurf (the
// curve-curve and curve-surface overloads used by
// PerformSetOfSurfOnElSpine).  The owning translations live outside
// Builder_2.cxx; the OCCT flow continues regardless (return values are
// ignored at these call sites).
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn simul_surf_cc_pending(
    _sd: &SharedSurfData,
    _hguide: &super::chfi_ds::ChFiDSElSpine,
    _spine: &super::chfi_ds::ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _hc1: &Option<BRepAdaptorCurve2d>,
    _hsref1: &BRepAdaptorSurface,
    _hcref1: &Option<BRepAdaptorCurve2d>,
    _decroch1: bool,
    _or1: rcad_kernel::topo::topods::Orientation,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _hc2: &Option<BRepAdaptorCurve2d>,
    _hsref2: &BRepAdaptorSurface,
    _hcref2: &Option<BRepAdaptorCurve2d>,
    _decroch2: bool,
    _or2: rcad_kernel::topo::topods::Orientation,
    _locfleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_p1: bool,
    _rec_rst1: bool,
    _rec_p2: bool,
    _rec_rst2: bool,
    _soldepcc: &[f64; 2],
) {
}

/// Pending PerformSurf — curve/curve overload (see simul_surf_cc_pending).
#[allow(clippy::too_many_arguments)]
fn perform_surf_cc_pending(
    _seqsd: &mut Vec<SharedSurfData>,
    _hguide: &super::chfi_ds::ChFiDSElSpine,
    _spine: &super::chfi_ds::ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _hc1: &Option<BRepAdaptorCurve2d>,
    _hsref1: &BRepAdaptorSurface,
    _hcref1: &Option<BRepAdaptorCurve2d>,
    _decroch1: bool,
    _or1: rcad_kernel::topo::topods::Orientation,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _hc2: &Option<BRepAdaptorCurve2d>,
    _hsref2: &BRepAdaptorSurface,
    _hcref2: &Option<BRepAdaptorCurve2d>,
    _decroch2: bool,
    _or2: rcad_kernel::topo::topods::Orientation,
    _max_step: f64,
    _locfleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_p1: bool,
    _rec_rst1: bool,
    _rec_p2: bool,
    _rec_rst2: bool,
    _soldepcc: &[f64; 2],
) {
}

/// Pending SimulSurf — curve-on-surface (obstacle on S1) overload.
#[allow(clippy::too_many_arguments)]
fn simul_surf_cs1_pending(
    _sd: &SharedSurfData,
    _hguide: &super::chfi_ds::ChFiDSElSpine,
    _spine: &super::chfi_ds::ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _hc1: &Option<BRepAdaptorCurve2d>,
    _hsref1: &BRepAdaptorSurface,
    _hcref1: &Option<BRepAdaptorCurve2d>,
    _decroch1: bool,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _or2: rcad_kernel::topo::topods::Orientation,
    _locfleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_p1: bool,
    _rec_s2: bool,
    _rec_rst1: bool,
    _soldepcs: &[f64; 3],
) {
}

/// Pending PerformSurf — curve-on-surface (obstacle on S1) overload.
#[allow(clippy::too_many_arguments)]
fn perform_surf_cs1_pending(
    _seqsd: &mut Vec<SharedSurfData>,
    _hguide: &super::chfi_ds::ChFiDSElSpine,
    _spine: &super::chfi_ds::ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _hc1: &Option<BRepAdaptorCurve2d>,
    _hsref1: &BRepAdaptorSurface,
    _hcref1: &Option<BRepAdaptorCurve2d>,
    _decroch1: bool,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _or2: rcad_kernel::topo::topods::Orientation,
    _max_step: f64,
    _locfleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_p1: bool,
    _rec_s2: bool,
    _rec_rst1: bool,
    _soldepcs: &[f64; 3],
) {
}

/// Pending SimulSurf — curve-on-surface (obstacle on S2) overload.
#[allow(clippy::too_many_arguments)]
fn simul_surf_cs2_pending(
    _sd: &SharedSurfData,
    _hguide: &super::chfi_ds::ChFiDSElSpine,
    _spine: &super::chfi_ds::ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _or1: rcad_kernel::topo::topods::Orientation,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _hc2: &Option<BRepAdaptorCurve2d>,
    _hsref2: &BRepAdaptorSurface,
    _hcref2: &Option<BRepAdaptorCurve2d>,
    _decroch2: bool,
    _locfleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_p2: bool,
    _rec_s1: bool,
    _rec_rst2: bool,
    _soldepcs: &[f64; 3],
) {
}

/// Pending PerformSurf — curve-on-surface (obstacle on S2) overload.
#[allow(clippy::too_many_arguments)]
fn perform_surf_cs2_pending(
    _seqsd: &mut Vec<SharedSurfData>,
    _hguide: &super::chfi_ds::ChFiDSElSpine,
    _spine: &super::chfi_ds::ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _or1: rcad_kernel::topo::topods::Orientation,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _hc2: &Option<BRepAdaptorCurve2d>,
    _hsref2: &BRepAdaptorSurface,
    _hcref2: &Option<BRepAdaptorCurve2d>,
    _decroch2: bool,
    _max_step: f64,
    _locfleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_p2: bool,
    _rec_s1: bool,
    _rec_rst2: bool,
    _soldepcs: &[f64; 3],
) {
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx ChFi3d_MKBound overloads (used by
// PerformSetOfKGen L3468/L3485/L3490/L3583) — pending TKGeomAlgo
// translation; the GeomFill_Boundary carriers are returned empty.
// =========================================================================

/// OCCT ChFi3d_mkbound(S, PC, Pref, P1, V1, Pref2, P2, V2, tapp, tg2d) —
/// the 2d-tangent variant.
#[allow(clippy::too_many_arguments)]
fn chfi3d_mkbound_c2d_tangents(
    _s: &BRepAdaptorSurface,
    _pc: &mut Option<rcad_kernel::geom::Curve2d>,
    _pref: i32,
    _p1: DVec2,
    _v1: DVec2,
    _pref2: i32,
    _p2: DVec2,
    _v2: DVec2,
    _tapp: f64,
    _tg2d: f64,
) -> GeomFillBoundary {
    GeomFillBoundary
}

/// OCCT ChFi3d_mkbound(S, PC, Pref, P1, V1, Pref2, P2, V2, tapp, tg3d) —
/// the 3d-tangent variant.
#[allow(clippy::too_many_arguments)]
fn chfi3d_mkbound_3d_tangents(
    _s: &BRepAdaptorSurface,
    _pc: &mut Option<rcad_kernel::geom::Curve2d>,
    _pref: i32,
    _p1: DVec2,
    _v1: DVec3,
    _pref2: i32,
    _p2: DVec2,
    _v2: DVec3,
    _tapp: f64,
    _tg3d: f64,
) -> GeomFillBoundary {
    GeomFillBoundary
}

/// OCCT ChFi3d_mkbound(S, PC, tapp, tg2d) — the simple iso-boundary variant.
fn chfi3d_mkbound_simple(
    _s: &BRepAdaptorSurface,
    _pc: &mut Option<rcad_kernel::geom::Curve2d>,
    _tapp: f64,
    _tg2d: f64,
) -> GeomFillBoundary {
    GeomFillBoundary
}

/// OCCT ChFi3d_mkbound(S, P1, P2, tapp, tg2d) — the two-point surface
/// boundary variant.
fn chfi3d_mkbound_two_points(
    _s: &Surface3,
    _p1: DVec2,
    _p2: DVec2,
    _tapp: f64,
    _tg2d: f64,
) -> GeomFillBoundary {
    GeomFillBoundary
}

/// OCCT ChFi3d_Builder_0.cxx ChFi3d_IntTraces — pending translation; the
/// OCCT failure path (traces do not intersect) is reported.  See the
/// PerformSetOfKGen reprocessing block for the call-site shape.
#[allow(clippy::too_many_arguments)]
fn chfi3d_int_traces_pending(
    _prevsd: &ChFiDSSurfData,
    _prevpar: f64,
    _nprevpar: &mut f64,
    _pref1: i32,
    _precs: i32,
    _nextsd: &ChFiDSSurfData,
    _nextpar: f64,
    _nnextpar: &mut f64,
    _nref1: i32,
    _nrefs: i32,
    _p2d: &mut DVec2,
    _isref: bool,
    _precaution: bool,
) -> bool {
    false
}

/// OCCT ChFi3d_Builder::CompleteData (declared in ChFi3d_Builder.hxx,
/// defined outside Builder_2.cxx) — pending owning translation; the OCCT
/// failure path (Standard_False) is reported.
#[allow(clippy::too_many_arguments)]
fn complete_data_pending(
    _cursd: &ChFiDSSurfData,
    _newsurf: &Surface3,
    _s1: &BRepAdaptorSurface,
    _pc1: &Option<rcad_kernel::geom::Curve2d>,
    _s2: &BRepAdaptorSurface,
    _pc2: &Option<rcad_kernel::geom::Curve2d>,
    _face_orientation: rcad_kernel::topo::topods::Orientation,
    _b1: bool,
    _b2: bool,
    _b3: bool,
    _b4: bool,
    _b5: bool,
) -> bool {
    false
}

impl super::chfi3d::ChFi3dBuilder {
    /// OCCT ChFi3d_Builder_2.cxx L2207-3000 — PerformSetOfSurfOnElSpine.
    pub fn perform_set_of_surf_on_el_spine(
        &mut self,
        hguide: &mut super::chfi_ds::ChFiDSElSpine,
        stripe: &SharedStripe,
        it1: &mut BRepTopAdaptorTopolTool,
        it2: &mut BRepTopAdaptorTopolTool,
        simul: bool,
    ) {
        // OCCT L2213: ChFiDS_ElSpine& Guide = *HGuide.
        // OCCT L2216: occ::handle<ChFiDS_Spine>& Spine = Stripe->ChangeSpine()
        // (the rcad spine is taken from the stripe and stored back on exit).
        let mut spine = {
            let mut st = stripe.write().expect("stripe lock");
            st.my_spine.take()
        }
        .expect("null spine");
        // OCCT L2236: the DS stays on self until the working copy below.
        let mut dstr = self.my_ds.take().expect("DS");

        // OCCT L2215/L2219-2229: OffsetHGuide — the parallel offset ElSpine
        // of HGuide located by handle identity in ChangeElSpines /
        // ChangeOffsetElSpines.  rcad stores the ElSpine lists by value
        // (no handle identity) — the pairing is a reported architecture gap
        // (it only matters for the ConstThroatWithPenetrationChamfer mode).
        let mut offset_hguide: Option<super::chfi_ds::ChFiDSElSpine> = None;
        if spine.base().mode() == ChFiDS_ChamfMode::ConstThroatWithPenetrationChamfer {
            // pending: handle-identity pairing.
        }

        let guide = hguide;
        let mut wf = guide.firstparam;
        let mut wl = guide.lastparam;
        let locfleche = (wl - wf) * self.fleche;
        let (wfsav, wlsav) = (wf, wl);
        if !guide.is_periodic() {
            // Now the ElSpine is artificially extended to help rsnld.
            let prab = 0.01;
            guide.set_first_parameter(wf - prab * (wl - wf));
            guide.set_last_parameter(wl + prab * (wl - wf));
            if let Some(og) = offset_hguide.as_mut() {
                og.set_first_parameter(wf - prab * (wl - wf));
                og.set_last_parameter(wl + prab * (wl - wf));
            }
        }
        // occ::handle<ChFiDS_Spine>&  Spine = Stripe->ChangeSpine();
        let mut nbed = spine.base().nb_edges() as i32;
        let lastedlastp = spine.base().last_parameter_of(spine.base().nb_edges());

        let mut ref_ = guide.previous.clone(); // OCCT L2254: Guide.Previous()
        let mut refbis: Option<SharedSurfData> = None;
        let raf = guide.next.clone(); // OCCT L2256: Guide.Next()
        {
            let mut st = stripe.write().expect("stripe lock");
            remove_sd(&mut st, &ref_, &raf);
        }

        let mut hs1 = BRepAdaptorSurface::empty();
        let mut hs2 = BRepAdaptorSurface::empty();
        let mut hs3: Option<BRepAdaptorSurface> = None;
        let mut hs4: Option<BRepAdaptorSurface> = None;
        let mut hsref1 = BRepAdaptorSurface::empty();
        let mut hsref2 = BRepAdaptorSurface::empty();
        let mut hc1: Option<BRepAdaptorCurve2d> = None;
        let mut hc2: Option<BRepAdaptorCurve2d> = None;
        let mut hcref1 = Some(BRepAdaptorCurve2d::new());
        let mut hcref2 = Some(BRepAdaptorCurve2d::new());
        let (mut decroch1, mut decroch2) = (false, false);
        let (mut rec_p1, mut rec_s1, mut rec_rst1, mut obstacleon1) =
            (false, false, false, false);
        let (mut rec_p2, mut rec_s2, mut rec_rst2, mut obstacleon2) =
            (false, false, false, false);
        let (mut pp1, mut pp2, mut pp3, mut pp4) =
            (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        let (mut w1, mut w2) = (0.0f64, 0.0f64);
        let mut soldep = [0.0f64; 4];
        let mut soldepcs = [0.0f64; 3];
        let mut soldepcc = [0.0f64; 2];

        // Restore a neighboring KPart.
        // If no neighbor calculation start point.
        let mut forward = true;
        let mut inside = false;
        let mut first = wf;
        let mut last = wl;
        let (mut ok1, mut ok2) = (true, true);
        // Restore the next KPart if it exists
        let mut vref = Shape::null();
        if ref_.is_none() && raf.is_none() {
            // sinon solution approchee.
            inside = true;

            self.start_sol_on_stripe(
                stripe, guide, &mut hs1, &mut hs2, it1, it2, &mut pp1, &mut pp2, &mut first,
            );

            last = wf;
            if guide.is_periodic() {
                last = first - guide.period();
                guide.save_first_parameter();
                guide.set_first_parameter(last);
                guide.save_last_parameter();
                guide.set_last_parameter(first * 1.1); // Extension to help rsnld.
                if let Some(og) = offset_hguide.as_mut() {
                    og.save_first_parameter();
                    og.set_first_parameter(last);
                    og.save_last_parameter();
                    og.set_last_parameter(first * 1.1); // Extension to help rsnld.
                }
            }
        } else {
            if !spine.base().is_periodic() && (wl - lastedlastp > -self.tolesp) {
                vref = spine.base().last_vertex();
            }
            if ref_.is_none() {
                if !spine.base().is_periodic() && (wf < self.tolesp) {
                    vref = spine.base().first_vertex();
                }
                ref_ = raf.clone();
                forward = false;
                first = wl;
                last = guide.firstparam;
            }

            let refs = ref_.as_ref().unwrap().read().expect("surfdata lock");
            ok1 = self.start_sol(
                &mut spine, &mut hs1, &mut pp1, &mut hc1, &mut w1, &refs, !forward, 1,
                &mut hsref1, &mut hcref1, &mut rec_p1, &mut rec_s1, &mut rec_rst1,
                &mut obstacleon1, &mut hs3, &mut pp3, decroch1, &vref,
            );
            ok2 = self.start_sol(
                &mut spine, &mut hs2, &mut pp2, &mut hc2, &mut w2, &refs, !forward, 2,
                &mut hsref2, &mut hcref2, &mut rec_p2, &mut rec_s2, &mut rec_rst2,
                &mut obstacleon2, &mut hs4, &mut pp4, decroch2, &vref,
            );
            drop(refs);
            hc1 = None;
            hc2 = None;

            if ok1 && ok2 {
                // OCCT L2366: Ok1 == 1 && Ok2 == 1.
                if forward {
                    guide.set_first_parameter(wf);
                    if let Some(og) = offset_hguide.as_mut() {
                        og.set_first_parameter(wf);
                    }
                } else {
                    guide.set_last_parameter(wl);
                    if let Some(og) = offset_hguide.as_mut() {
                        og.set_last_parameter(wl);
                    }
                }
            }
        }
        let mut fini = false;
        let mut complete = inside;
        if !guide.is_periodic() {
            let indf = spine.base().index_of_param(wf, true) as i32;
            let mut indl = spine.base().index_of_param(wl, false) as i32;
            if spine.base().is_periodic() && (indl < indf) {
                indl += nbed;
            }
            nbed = indl - indf + 1;
        }
        // No Max at the touch : 20 points by edge at average without
        // counting the extensions.

        let mut bidf = wf;
        let mut bidl = wl;
        if !spine.base().is_periodic() {
            bidf = 0.0f64.max(wf);
            bidl = wl.min(spine.base().last_parameter_of(spine.base().nb_edges()));
            // PMN 20/07/98 : Attention in case if there is only extension
            if (bidl - bidf) < 0.01 * spine.base().last_parameter_of(spine.base().nb_edges()) {
                bidf = wf;
                bidl = wl;
            }
        }
        let max_step = (bidl - bidf) * 0.05 / nbed as f64;
        let mut firstsov = 0.0f64;
        let (mut intf, mut intl) = (0i32, 0i32);
        while !fini {
            // are these the ends (no extension on periodic).
            ok1 = true;
            ok2 = true;
            if !spine.base().is_periodic() {
                if wf < self.tolesp && (complete == inside) {
                    if spine.base().first_status() == ChFiDS_State::OnSame {
                        intf = 2;
                    } else {
                        intf = 1;
                    }
                }
                if spine.base().is_tangency_extremity(true) {
                    intf = 4;
                    guide.set_first_parameter(wfsav);
                    if let Some(og) = offset_hguide.as_mut() {
                        og.set_first_parameter(wfsav);
                    }
                }
                if wl - lastedlastp > -self.tolesp {
                    if spine.base().last_status() == ChFiDS_State::OnSame {
                        intl = 2;
                    } else {
                        intl = 1;
                    }
                }
                if spine.base().is_tangency_extremity(false) {
                    intl = 4;
                    guide.set_last_parameter(wlsav);
                    if let Some(og) = offset_hguide.as_mut() {
                        og.set_last_parameter(wlsav);
                    }
                }
            }
            if intf != 0 && !forward {
                vref = spine.base().first_vertex();
            }
            if intl != 0 && forward {
                vref = spine.base().last_vertex();
            }
            if ref_.is_some() {
                let refs = ref_.as_ref().unwrap().read().expect("surfdata lock");
                ok1 = self.start_sol(
                    &mut spine, &mut hs1, &mut pp1, &mut hc1, &mut w1, &refs, !forward, 1,
                    &mut hsref1, &mut hcref1, &mut rec_p1, &mut rec_s1, &mut rec_rst1,
                    &mut obstacleon1, &mut hs3, &mut pp3, decroch1, &vref,
                );
                ok2 = self.start_sol(
                    &mut spine, &mut hs2, &mut pp2, &mut hc2, &mut w2, &refs, !forward, 2,
                    &mut hsref2, &mut hcref2, &mut rec_p2, &mut rec_s2, &mut rec_rst2,
                    &mut obstacleon2, &mut hs4, &mut pp4, decroch2, &vref,
                );
            }

            // No more connected faces. Construction of the tangent plane to
            // continue the path till the output on the other face.
            if (!ok1 && hc1.is_none()) || (!ok2 && hc2.is_none()) {
                if (intf != 0 && !forward) || (intl != 0 && forward) {
                    let refs = ref_.as_ref().unwrap().read().expect("surfdata lock");
                    if !ok1 {
                        chfi3d_build_plane(&mut dstr, &mut hs1, &mut pp1, &refs, !forward, 1);
                    }
                    if !ok2 {
                        chfi3d_build_plane(&mut dstr, &mut hs2, &mut pp2, &refs, !forward, 2);
                    }
                    if intf != 0 {
                        intf = 5;
                    } else if intl != 0 {
                        intl = 5;
                    }
                    drop(refs);
                    if forward {
                        guide.set_first_parameter(wf);
                        if let Some(og) = offset_hguide.as_mut() {
                            og.set_first_parameter(wf);
                        }
                    } else {
                        guide.set_last_parameter(wl);
                        if let Some(og) = offset_hguide.as_mut() {
                            og.set_last_parameter(wl);
                        }
                    }
                } else {
                    panic!(
                        "Standard_Failure: PerformSetOfSurfOnElSpine : Chaining is impossible."
                    );
                }
            }

            // Definition of the domain of patch It1, It2
            it1.initialize_brep(&hs1);
            it2.initialize_brep(&hs2);

            // Calculate one (several if singularity) SurfaData
            let sd = std::sync::Arc::new(std::sync::RwLock::new(ChFiDSSurfData::default()));
            let mut seqsd: Vec<SharedSurfData> = Vec::new();
            seqsd.push(sd.clone());

            if obstacleon1 && obstacleon2 {
                let mut or1 = hsref1.face.orientation;
                let mut or2 = hsref2.face.orientation;
                let (so1, so2, sc) = {
                    let st = stripe.read().expect("stripe lock");
                    (
                        st.orientation_on_face1(),
                        st.orientation_on_face2(),
                        st.choix(),
                    )
                };
                let mut choix = next_side(&mut or1, &mut or2, so1, so2, sc);

                // Calculate the criterion of Choice edge / edge
                if choix % 2 == 0 {
                    choix = 4;
                } else {
                    choix = 1;
                }

                soldepcc[0] = w1;
                soldepcc[1] = w2;
                if simul {
                    simul_surf_cc_pending(
                        &sd, guide, &spine, choix, &hs1, it1, &hc1, &hsref1, &hcref1, decroch1,
                        or1, &hs2, it2, &hc2, &hsref2, &hcref2, decroch2, or2, locfleche,
                        self.tolesp, &mut first, &mut last, inside, inside, forward, rec_p1,
                        rec_rst1, rec_p2, rec_rst2, &soldepcc,
                    );
                } else {
                    perform_surf_cc_pending(
                        &mut seqsd, guide, &spine, choix, &hs1, it1, &hc1, &hsref1, &hcref1,
                        decroch1, or1, &hs2, it2, &hc2, &hsref2, &hcref2, decroch2, or2, max_step,
                        locfleche, self.tolesp, &mut first, &mut last, inside, inside, forward,
                        rec_p1, rec_rst1, rec_p2, rec_rst2, &soldepcc,
                    );
                }
                {
                    let mut sdw = sd.write().expect("surfdata lock");
                    sdw.change_index_of_s1(dstr.add_shape(&hs1.face));
                    sdw.change_index_of_s2(dstr.add_shape(&hs2.face));
                }
            } else if obstacleon1 {
                let mut or1 = hsref1.face.orientation;
                let mut or2 = hs2.face.orientation;
                let (so1, so2, sc) = {
                    let st = stripe.read().expect("stripe lock");
                    (
                        st.orientation_on_face1(),
                        st.orientation_on_face2(),
                        -st.choix(),
                    )
                };
                let mut choix = next_side(&mut or1, &mut or2, so1, so2, sc);
                if choix % 2 == 1 {
                    choix += 1;
                } else {
                    choix -= 1;
                }
                soldepcs[2] = w1;
                soldepcs[0] = pp2.x;
                soldepcs[1] = pp2.y;
                if simul {
                    simul_surf_cs1_pending(
                        &sd, guide, &spine, choix, &hs1, it1, &hc1, &hsref1, &hcref1, decroch1,
                        &hs2, it2, or2, locfleche, self.tolesp, &mut first, &mut last, inside,
                        inside, forward, rec_p1, rec_s2, rec_rst1, &soldepcs,
                    );
                } else {
                    perform_surf_cs1_pending(
                        &mut seqsd, guide, &spine, choix, &hs1, it1, &hc1, &hsref1, &hcref1,
                        decroch1, &hs2, it2, or2, max_step, locfleche, self.tolesp, &mut first,
                        &mut last, inside, inside, forward, rec_p1, rec_s2, rec_rst1, &soldepcs,
                    );
                }
                {
                    let mut sdw = sd.write().expect("surfdata lock");
                    sdw.change_index_of_s1(dstr.add_shape(&hs1.face));
                    sdw.change_index_of_s2(dstr.add_shape(&hs2.face));
                }
                decroch2 = false;
            } else if obstacleon2 {
                let mut or1 = hs1.face.orientation;
                let mut or2 = hsref2.face.orientation;
                let (so1, so2, sc) = {
                    let st = stripe.read().expect("stripe lock");
                    (
                        st.orientation_on_face1(),
                        st.orientation_on_face2(),
                        st.choix(),
                    )
                };
                let choix = next_side(&mut or1, &mut or2, so1, so2, sc);
                soldepcs[2] = w2;
                soldepcs[0] = pp1.x;
                soldepcs[1] = pp1.y;
                if simul {
                    simul_surf_cs2_pending(
                        &sd, guide, &spine, choix, &hs1, it1, or1, &hs2, it2, &hc2, &hsref2,
                        &hcref2, decroch2, locfleche, self.tolesp, &mut first, &mut last, inside,
                        inside, forward, rec_p2, rec_s1, rec_rst2, &soldepcs,
                    );
                } else {
                    perform_surf_cs2_pending(
                        &mut seqsd, guide, &spine, choix, &hs1, it1, or1, &hs2, it2, &hc2,
                        &hsref2, &hcref2, decroch2, max_step, locfleche, self.tolesp, &mut first,
                        &mut last, inside, inside, forward, rec_p2, rec_s1, rec_rst2, &soldepcs,
                    );
                }
                {
                    let mut sdw = sd.write().expect("surfdata lock");
                    sdw.change_index_of_s1(dstr.add_shape(&hs1.face));
                    sdw.change_index_of_s2(dstr.add_shape(&hs2.face));
                }
                decroch1 = false;
            } else {
                // OCCT L2842-2843: the Surf1/Surf2 out handles ARE HS1/HS2
                // (the same handles); rcad reassigns through local slots.
                let mut surf1 = Some(hs1.clone());
                let mut surf2 = Some(hs2.clone());
                self.call_perform_surf(
                    stripe, simul, &mut seqsd, &sd, guide, &spine, &mut hs1, &mut hs3, pp1, pp3,
                    it1, &mut hs2, &mut hs4, pp2, pp4, it2, max_step, locfleche, &mut first,
                    &mut last, inside, forward, rec_s1, rec_s2, &mut soldep, &mut intf,
                    &mut intl, &mut surf1, &mut surf2,
                );
                if let Some(s) = surf1.take() {
                    hs1 = s;
                }
                if let Some(s) = surf2.take() {
                    hs2 = s;
                }
                decroch1 = false;
                decroch2 = false;
            }

            if !self.done {
                // Case of fail
                if (!ok1 && !obstacleon1) || (!ok2 && !obstacleon2) {
                    // Fail in a part of extension is not serious
                    // Here one stops.
                    self.done = true;
                    inside = false;
                    if forward {
                        intl = 1;
                    } else {
                        intf = 1;
                    }
                } else {
                    // Otherwise invalidation of the stripe.
                    spine
                        .base_mut()
                        .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
                    panic!("Standard_Failure: CallPerformSurf : Path failed!");
                }
            } else {
                refbis = ref_.clone();
                let sd_used = if forward {
                    seqsd.last().unwrap().clone()
                } else {
                    seqsd.first().unwrap().clone()
                };
                if forward {
                    let mut ii = 1usize; // OCCT 1-based
                    while ii <= seqsd.len() {
                        let sdi = seqsd[ii - 1].clone();
                        {
                            let mut w = sdi.write().expect("surfdata lock");
                            w.change_index_of_s1(dstr.add_shape(&hs1.face));
                            if obstacleon1 {
                                w.set_index_of_c1(
                                    dstr.add_shape(hc1.as_ref().expect("HC1").edge()),
                                );
                            }
                            w.change_index_of_s2(dstr.add_shape(&hs2.face));
                            if obstacleon2 {
                                w.set_index_of_c2(
                                    dstr.add_shape(hc2.as_ref().expect("HC2").edge()),
                                );
                            }
                        }
                        {
                            let mut st = stripe.write().expect("stripe lock");
                            insert_after(&mut st, &refbis, &sdi);
                        }
                        refbis = Some(sdi);
                        ii += 1;
                    }
                } else {
                    let mut ii = seqsd.len(); // OCCT 1-based: SeqSD.Length()
                    while ii >= 1 {
                        let sdi = seqsd[ii - 1].clone();
                        {
                            let mut w = sdi.write().expect("surfdata lock");
                            w.change_index_of_s1(dstr.add_shape(&hs1.face));
                            if obstacleon1 {
                                w.set_index_of_c1(
                                    dstr.add_shape(hc1.as_ref().expect("HC1").edge()),
                                );
                            }
                            w.change_index_of_s2(dstr.add_shape(&hs2.face));
                            if obstacleon2 {
                                w.set_index_of_c2(
                                    dstr.add_shape(hc2.as_ref().expect("HC2").edge()),
                                );
                            }
                        }
                        {
                            let mut st = stripe.write().expect("stripe lock");
                            insert_before(&mut st, &refbis, &sdi);
                        }
                        refbis = Some(sdi);
                        ii -= 1;
                    }
                }

                // OCCT L2916/L2922: the SD argument is the loop variable of
                // the insertion loop above (its last processed value).
                if !ok1 && !obstacleon1 {
                    // clean infos on the plane of extension.
                    let vref_cp = ref_
                        .as_ref()
                        .unwrap()
                        .read()
                        .expect("surfdata lock")
                        .vertex(!forward, 1)
                        .clone();
                    let mut st = stripe.write().expect("stripe lock");
                    chfi3d_purge(&mut st, &sd_used, &vref_cp, !forward, 1, &mut intf, &mut intl);
                }

                if !ok2 && !obstacleon2 {
                    // clean infos on the plane of extension.
                    let vref_cp = ref_
                        .as_ref()
                        .unwrap()
                        .read()
                        .expect("surfdata lock")
                        .vertex(!forward, 2)
                        .clone();
                    let mut st = stripe.write().expect("stripe lock");
                    chfi3d_purge(&mut st, &sd_used, &vref_cp, !forward, 2, &mut intf, &mut intl);
                }

                // The end. The reference is changed.
                ref_ = refbis;
            }

            if inside {
                // There are starting solutions for the next.
                inside = false;
                firstsov = first;
                if guide.is_periodic() {
                    complete = false;
                    wf = guide.firstparam;
                    wl = guide.lastparam;
                }
            }
            if forward {
                fini = (wl - last) <= 10.0 * self.tolesp
                    || (intl != 0 && !(obstacleon1 || obstacleon2)); // General case

                if !fini && guide.is_periodic() && ((wl - last) < guide.period() * 1.0e-3) {
                    // It is tested if reframing of extremes is done at the
                    // same edge — Loop Condition
                    let (thefirst, thelast) = {
                        let st = stripe.read().expect("stripe lock");
                        (
                            st.set_of_surf_data().first().expect("empty").clone(),
                            st.set_of_surf_data().last().expect("empty").clone(),
                        )
                    };

                    let cond1 = {
                        let (tf, tl) = (
                            thefirst.read().expect("surfdata lock"),
                            thelast.read().expect("surfdata lock"),
                        );
                        tf.vertex_first_on_s1().is_on_arc()
                            && tl.vertex_last_on_s1().is_on_arc()
                    };
                    if cond1 {
                        fini = {
                            let (tf, tl) = (
                                thefirst.read().expect("surfdata lock"),
                                thelast.read().expect("surfdata lock"),
                            );
                            tf.vertex_first_on_s1()
                                .arc()
                                .is_same(tl.vertex_last_on_s1().arc())
                        };
                    }
                    if !fini {
                        let cond2 = {
                            let (tf, tl) = (
                                thefirst.read().expect("surfdata lock"),
                                thelast.read().expect("surfdata lock"),
                            );
                            tf.vertex_first_on_s2().is_on_arc()
                                && tl.vertex_last_on_s2().is_on_arc()
                        };
                        if cond2 {
                            fini = {
                                let (tf, tl) = (
                                    thefirst.read().expect("surfdata lock"),
                                    thelast.read().expect("surfdata lock"),
                                );
                                tf.vertex_first_on_s2()
                                    .arc()
                                    .is_same(tl.vertex_last_on_s2().arc())
                            };
                        }
                    }

                    if fini {
                        // It is ended!
                        {
                            let mut st = stripe.write().expect("stripe lock");
                            st.my_spine = Some(spine);
                        }
                        self.my_ds = Some(dstr);
                        return;
                    }
                }

                if fini && complete {
                    // restart in the opposite direction.
                    ref_ = {
                        let st = stripe.read().expect("stripe lock");
                        Some(st.set_of_surf_data().first().expect("empty").clone())
                    };
                    forward = false;
                    fini = false;
                    first = firstsov;
                } else {
                    first = last;
                    last = wl;
                }
            }
            if !forward {
                fini = (first - wf) <= 10.0 * self.tolesp
                    || (intf != 0 && !(obstacleon1 || obstacleon2));
                complete = false;
                last = wf;
            }
        }
        // The initial state is restored
        if !guide.is_periodic() {
            guide.set_first_parameter(wfsav);
            guide.set_last_parameter(wlsav);
            if let Some(og) = offset_hguide.as_mut() {
                og.set_first_parameter(wfsav);
                og.set_last_parameter(wlsav);
            }
        }
        {
            let mut st = stripe.write().expect("stripe lock");
            st.my_spine = Some(spine);
        }
        self.my_ds = Some(dstr);
    }

    /// OCCT ChFi3d_Builder_2.cxx L3298-3878 — PerformSetOfKGen.
    pub fn perform_set_of_k_gen(&mut self, stripe: &SharedStripe, simul: bool) {
        let mut it1 = BRepTopAdaptorTopolTool::default();
        let mut it2 = BRepTopAdaptorTopolTool::default();
        let mut spine = {
            let mut st = stripe.write().expect("stripe lock");
            st.my_spine.take()
        }
        .expect("null spine");
        // OCCT L3303: NCollection_List<handle<ChFiDS_ElSpine>>& ll =
        //   Spine->ChangeElSpines();  (rcad: taken by value, restored below)
        let mut ll = std::mem::take(&mut spine.base_mut().elspines);
        {
            let mut iles = 0usize;
            while iles < ll.len() {
                self.perform_set_of_surf_on_el_spine(&mut ll[iles], stripe, &mut it1, &mut it2, simul);
                iles += 1;
            }
        }
        if !simul {
            let mut dstr = self.my_ds.take().expect("DS");
            let brep = self.my_brep.clone();
            let seqsurf: Vec<SharedSurfData> = {
                let st = stripe.read().expect("stripe lock");
                st.set_of_surf_data().clone()
            };
            let len = seqsurf.len();
            let mut lastn = len;
            let periodic = spine.base().is_periodic();
            if periodic {
                lastn += 1;
            }
            // It is attempted to reprocess the squares that bore.
            let mut i = 1usize; // OCCT 1-based
            while i <= len {
                let cursd = seqsurf[i - 1].clone();
                let (tw1, tw2) = {
                    let g = cursd.read().expect("surfdata lock");
                    (g.twist_on_s1(), g.twist_on_s2())
                };
                let mut prevsd: Option<SharedSurfData> = None;
                let mut nextsd: Option<SharedSurfData> = None;
                let mut iprev = i as i64 - 1;
                if iprev == 0 && periodic {
                    iprev = len as i64;
                }
                let mut inext = i as i64 + 1;
                if inext > len as i64 {
                    if periodic {
                        inext = 1;
                    } else {
                        inext = 0;
                    }
                }

                // For the moment only the surfaces where the twist is
                // detected at the path are corrected, it is necessary to
                // control more subtly the ugly traces (size, curvature,
                // inflexion... )
                if !tw1 && !tw2 {
                    i += 1;
                    continue;
                }

                // It is decided (fairly at random) if the extended surface is
                // ready for the filling.
                let (cursurf1, cursurf2, ddeb, dfin, don1, don2) = {
                    let g = cursd.read().expect("surfdata lock");
                    let cpd1 = g.vertex_first_on_s1().point();
                    let cpd2 = g.vertex_first_on_s2().point();
                    let cpf1 = g.vertex_last_on_s1().point();
                    let cpf2 = g.vertex_last_on_s2().point();
                    (
                        g.index_of_s1(),
                        g.index_of_s2(),
                        cpd1.distance(cpd2),
                        cpf1.distance(cpf2),
                        cpd1.distance(cpf1),
                        cpd2.distance(cpf2),
                    )
                };
                let possibleon1 = don1 < 2.0 * (ddeb + dfin);
                let possibleon2 = don2 < 2.0 * (ddeb + dfin);
                if (tw1 && !possibleon1) || (tw2 && !possibleon2) {
                    spine
                        .base_mut()
                        .set_error_status(ChFiDS_ErrorStatus::TwistedSurface);
                    panic!(
                        "Standard_Failure: adjustment by reprocessing the non-written points"
                    );
                }

                // It is checked if there are presentable neighbors
                let mut yaprevon1 = false;
                let mut yaprevon2 = false;
                let mut samesurfon1 = false;
                let mut samesurfon2 = false;
                if iprev != 0 {
                    prevsd = Some(seqsurf[(iprev - 1) as usize].clone());
                    let g = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                    yaprevon1 = !g.twist_on_s1();
                    samesurfon1 = g.index_of_s1() == cursurf1;
                    yaprevon2 = !g.twist_on_s2();
                    samesurfon2 = g.index_of_s2() == cursurf2;
                }
                let mut yanexton1 = false;
                let mut yanexton2 = false;
                if inext != 0 {
                    nextsd = Some(seqsurf[(inext - 1) as usize].clone());
                    let g = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                    yanexton1 = !g.twist_on_s1();
                    if samesurfon1 {
                        samesurfon1 = g.index_of_s1() == cursurf1;
                    }
                    yanexton2 = !g.twist_on_s2();
                    if samesurfon2 {
                        samesurfon2 = g.index_of_s2() == cursurf2;
                    }
                }
                // A contour of filling is constructed
                let (mut pc1, mut pc2, f1, f2) = {
                    let g = cursd.read().expect("surfdata lock");
                    (
                        g.interference_on_s1().pcurve_on_face().cloned(),
                        g.interference_on_s2().pcurve_on_face().cloned(),
                        dstr.shape(cursurf1).clone(),
                        dstr.shape(cursurf2).clone(),
                    )
                };
                let s1 = BRepAdaptorSurface::initialize(&brep, &f1);
                let s2 = BRepAdaptorSurface::initialize(&brep, &f2);
                let mut bon1 = GeomFillBoundary;
                let mut bon2 = GeomFillBoundary;
                let mut bdeb = GeomFillBoundary;
                let mut bfin = GeomFillBoundary;
                let mut pointuon1 = false;
                let mut pointuon2 = false;
                if tw1 {
                    if !yaprevon1 || !yanexton1 {
                        spine
                            .base_mut()
                            .set_error_status(ChFiDS_ErrorStatus::TwistedSurface);
                        panic!(
                            "Standard_Failure: adjustment by reprocessing the non-written points: no neighbor"
                        );
                    }
                    let (prevpar1, nextpar1) = {
                        let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                        let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                        (
                            pg.interference_on_s1().parameter_last(),
                            ng.interference_on_s1().parameter_first(),
                        )
                    };
                    if samesurfon1 {
                        // It is checked if it is possible to intersect traces
                        // of neighbors to create a sharp end.
                        let (pcprev1, pcnext1) = {
                            let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                            let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                            (
                                pg.interference_on_s1().pcurve_on_face().cloned(),
                                ng.interference_on_s1().pcurve_on_face().cloned(),
                            )
                        };
                        let mut nprevpar1 = 0.0f64;
                        let mut nnextpar1 = 0.0f64;
                        let mut p2d = DVec2::ZERO;
                        if chfi3d_int_traces_pending(
                            &prevsd.as_ref().unwrap().read().expect("surfdata lock"),
                            prevpar1,
                            &mut nprevpar1,
                            1,
                            1,
                            &nextsd.as_ref().unwrap().read().expect("surfdata lock"),
                            nextpar1,
                            &mut nnextpar1,
                            1,
                            -1,
                            &mut p2d,
                            false,
                            true,
                        ) {
                            prevsd
                                .as_ref()
                                .unwrap()
                                .write()
                                .expect("surfdata lock")
                                .change_interference_on_s1()
                                .set_last_parameter(nprevpar1);
                            nextsd
                                .as_ref()
                                .unwrap()
                                .write()
                                .expect("surfdata lock")
                                .change_interference_on_s1()
                                .set_first_parameter(nnextpar1);
                            pointuon1 = true;
                            pc1 = None;
                        } else {
                            let (pdeb1, vdeb1) = {
                                let pc = pcprev1.as_ref().expect("pcurve");
                                (pc.point_at(prevpar1), pc.derivative_at(prevpar1))
                            };
                            let (pfin1, vfin1) = {
                                let pc = pcnext1.as_ref().expect("pcurve");
                                (pc.point_at(nextpar1), pc.derivative_at(nextpar1))
                            };
                            bon1 = chfi3d_mkbound_c2d_tangents(
                                &s1,
                                &mut pc1,
                                -1,
                                pdeb1,
                                vdeb1,
                                1,
                                pfin1,
                                vfin1,
                                self.tolapp3d,
                                2.0e-4,
                            );
                        }
                    } else {
                        // here the base is on 3D tangents of neighbors.
                        let (c3dprev1, c3dnext1) = {
                            let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                            let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                            (
                                dstr.curve(pg.interference_on_s1().lineindex).curve.clone(),
                                dstr.curve(ng.interference_on_s1().lineindex).curve.clone(),
                            )
                        };
                        let _ = (c3dprev1, c3dnext1);
                        // OCCT L3478-3479: c3dprev1->D1(prevpar1, Pdeb1,
                        // Vdeb1) / c3dnext1->D1(nextpar1, Pfin1, Vfin1) —
                        // the DS curve handles are nullable until the
                        // walking stage fills them; the pending tangent
                        // queries stand in as zero vectors.
                        let (vdeb1_3d, vfin1_3d) = (DVec3::ZERO, DVec3::ZERO); // pending c3d
                        let (pardeb1, parfin1) = {
                            let g = cursd.read().expect("surfdata lock");
                            (
                                g.interference_on_s1().parameter_first(),
                                g.interference_on_s1().parameter_last(),
                            )
                        };
                        let (pdeb1, pfin1) = {
                            let pc = pc1.as_ref().expect("pcurve");
                            (pc.point_at(pardeb1), pc.point_at(parfin1))
                        };
                        bon1 = chfi3d_mkbound_3d_tangents(
                            &s1,
                            &mut pc1,
                            -1,
                            pdeb1,
                            vdeb1_3d,
                            1,
                            pfin1,
                            vfin1_3d,
                            self.tolapp3d,
                            2.0e-4,
                        );
                    }
                } else {
                    bon1 = chfi3d_mkbound_simple(&s1, &mut pc1, self.tolapp3d, 2.0e-4);
                }
                if tw2 {
                    if !yaprevon2 || !yanexton2 {
                        panic!(
                            "Standard_Failure: adjustment by reprocessing the non-written points: no neighbor"
                        );
                    }
                    let (prevpar2, nextpar2) = {
                        let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                        let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                        (
                            pg.interference_on_s2().parameter_last(),
                            ng.interference_on_s2().parameter_first(),
                        )
                    };
                    if samesurfon2 {
                        // It is checked if it is possible to intersect traces
                        // of neighbors to create a sharp end.
                        let (pcprev2, pcnext2) = {
                            let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                            let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                            (
                                pg.interference_on_s2().pcurve_on_face().cloned(),
                                ng.interference_on_s2().pcurve_on_face().cloned(),
                            )
                        };
                        let mut nprevpar2 = 0.0f64;
                        let mut nnextpar2 = 0.0f64;
                        let mut p2d = DVec2::ZERO;
                        if chfi3d_int_traces_pending(
                            &prevsd.as_ref().unwrap().read().expect("surfdata lock"),
                            prevpar2,
                            &mut nprevpar2,
                            2,
                            1,
                            &nextsd.as_ref().unwrap().read().expect("surfdata lock"),
                            nextpar2,
                            &mut nnextpar2,
                            2,
                            -1,
                            &mut p2d,
                            false,
                            true,
                        ) {
                            prevsd
                                .as_ref()
                                .unwrap()
                                .write()
                                .expect("surfdata lock")
                                .change_interference_on_s2()
                                .set_last_parameter(nprevpar2);
                            nextsd
                                .as_ref()
                                .unwrap()
                                .write()
                                .expect("surfdata lock")
                                .change_interference_on_s2()
                                .set_first_parameter(nnextpar2);
                            pointuon2 = true;
                            pc2 = None;
                        } else {
                            let (pdeb2, vdeb2) = {
                                let pc = pcprev2.as_ref().expect("pcurve");
                                (pc.point_at(prevpar2), pc.derivative_at(prevpar2))
                            };
                            let (pfin2, vfin2) = {
                                let pc = pcnext2.as_ref().expect("pcurve");
                                (pc.point_at(nextpar2), pc.derivative_at(nextpar2))
                            };
                            bon2 = chfi3d_mkbound_c2d_tangents(
                                &s2,
                                &mut pc2,
                                -1,
                                pdeb2,
                                vdeb2,
                                1,
                                pfin2,
                                vfin2,
                                self.tolapp3d,
                                2.0e-4,
                            );
                        }
                    } else {
                        // here the base is on 3D tangents of neighbors.
                        let (c3dprev2, c3dnext2) = {
                            let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                            let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                            (
                                dstr.curve(pg.interference_on_s2().lineindex).curve.clone(),
                                dstr.curve(ng.interference_on_s2().lineindex).curve.clone(),
                            )
                        };
                        let _ = (c3dprev2, c3dnext2);
                        // OCCT L3549-3550: c3dprev2->D1 / c3dnext2->D1 —
                        // pending DS curve handles (see the tw1 twin).
                        let (vdeb2_3d, vfin2_3d) = (DVec3::ZERO, DVec3::ZERO); // pending c3d
                        let (pardeb2, parfin2) = {
                            let g = cursd.read().expect("surfdata lock");
                            (
                                g.interference_on_s2().parameter_first(),
                                g.interference_on_s2().parameter_last(),
                            )
                        };
                        let (pdeb2, pfin2) = {
                            let pc = pc2.as_ref().expect("pcurve");
                            (pc.point_at(pardeb2), pc.point_at(parfin2))
                        };
                        bon2 = chfi3d_mkbound_3d_tangents(
                            &s2,
                            &mut pc2,
                            -1,
                            pdeb2,
                            vdeb2_3d,
                            1,
                            pfin2,
                            vfin2_3d,
                            self.tolapp3d,
                            2.0e-4,
                        );
                    }
                } else {
                    bon2 = chfi3d_mkbound_simple(&s2, &mut pc2, self.tolapp3d, 2.0e-4);
                }
                // The parameters of neighbor traces are updated, so
                // straight lines uv are pulled.
                let (sprev, snext, pdebs1, pdebs2, pfins1, pfins2) = {
                    let pg = prevsd.as_ref().unwrap().read().expect("surfdata lock");
                    let ng = nextsd.as_ref().unwrap().read().expect("surfdata lock");
                    let pcsprev1 = pg.interference_on_s1().pcurve_on_surf().cloned();
                    let pcsnext1 = ng.interference_on_s1().pcurve_on_surf().cloned();
                    let prevpar1 = pg.interference_on_s1().parameter_last();
                    let nextpar1 = ng.interference_on_s1().parameter_first();
                    let pcsprev2 = pg.interference_on_s2().pcurve_on_surf().cloned();
                    let pcsnext2 = ng.interference_on_s2().pcurve_on_surf().cloned();
                    let prevpar2 = pg.interference_on_s2().parameter_last();
                    let nextpar2 = ng.interference_on_s2().parameter_first();
                    (
                        dstr.surface(pg.surf()).surface.clone(),
                        dstr.surface(ng.surf()).surface.clone(),
                        pcsprev1
                            .map(|pc| pc.point_at(prevpar1))
                            .unwrap_or(DVec2::ZERO),
                        pcsprev2
                            .map(|pc| pc.point_at(prevpar2))
                            .unwrap_or(DVec2::ZERO),
                        pcsnext1
                            .map(|pc| pc.point_at(nextpar1))
                            .unwrap_or(DVec2::ZERO),
                        pcsnext2
                            .map(|pc| pc.point_at(nextpar2))
                            .unwrap_or(DVec2::ZERO),
                    )
                };
                bdeb = chfi3d_mkbound_two_points(&sprev, pdebs1, pdebs2, self.tolapp3d, 2.0e-4);
                bfin = chfi3d_mkbound_two_points(&snext, pfins1, pfins2, self.tolapp3d, 2.0e-4);

                let mut fil = GeomFillConstrainedFilling::new(11, 20);
                if pointuon1 {
                    fil.init3(&bon2, &bfin, &bdeb, true);
                } else if pointuon2 {
                    fil.init3(&bon1, &bfin, &bdeb, true);
                } else {
                    fil.init4(&bon1, &bfin, &bon2, &bdeb, true);
                }

                pc1 = pc1.map(|pc| chfi3d_reparam_pcurv(0.0, 1.0, pc));
                pc2 = pc2.map(|pc| chfi3d_reparam_pcurv(0.0, 1.0, pc));
                // OCCT L3602-3648: newsurf = fil.Surface(); the
                // GeomFill_ConstrainedFilling machinery is pending
                // (TKGeomAlgo) — until it lands Surface() reports None and
                // the reprocessing of the square stops here (see the stage
                // report).  The OCCT flow always produces a surface.
                if let Some(newsurf) = fil.surface() {
                    if pointuon1 {
                        let ns = newsurf.clone();
                        // OCCT: newsurf->VReverse() — the V parametrization
                        // flip is a pending Surface3 operation; this branch
                        // only becomes reachable once the GeomFill machinery
                        // lands (fil.surface() is None until then).
                        let _ = &ns;
                        self.done = {
                            let mut data = cursd.write().expect("surfdata lock");
                            self.complete_data_surfcoin(
                                &mut data,
                                &ns,
                                &s1,
                                pc1.as_ref(),
                                &s2,
                                pc2.as_ref(),
                                f2.orientation,
                                false,
                                false,
                                false,
                                false,
                                false,
                            )
                        };
                        cursd
                            .write()
                            .expect("surfdata lock")
                            .change_index_of_s1(0);
                    } else {
                        self.done = {
                            let mut data = cursd.write().expect("surfdata lock");
                            self.complete_data_surfcoin(
                                &mut data,
                                &newsurf,
                                &s1,
                                pc1.as_ref(),
                                &s2,
                                pc2.as_ref(),
                                f1.orientation,
                                true,
                                false,
                                false,
                                false,
                                false,
                            )
                        };
                        if pointuon2 {
                            cursd
                                .write()
                                .expect("surfdata lock")
                                .change_index_of_s2(0);
                        }
                    }
    let (cpd1, cpd2, cpf1, cpf2) = {
                        let g = cursd.read().expect("surfdata lock");
                        (
                            g.vertex_first_on_s1().point(),
                            g.vertex_first_on_s2().point(),
                            g.vertex_last_on_s1().point(),
                            g.vertex_last_on_s2().point(),
                        )
                    };
                    if tw1 {
                        prevsd
                            .as_ref()
                            .unwrap()
                            .write()
                            .expect("surfdata lock")
                            .change_vertex_last_on_s1()
                            .set_point(cpd1);
                        nextsd
                            .as_ref()
                            .unwrap()
                            .write()
                            .expect("surfdata lock")
                            .change_vertex_first_on_s1()
                            .set_point(cpf1);
                    }
                    if tw2 {
                        prevsd
                            .as_ref()
                            .unwrap()
                            .write()
                            .expect("surfdata lock")
                            .change_vertex_last_on_s2()
                            .set_point(cpd2);
                        nextsd
                            .as_ref()
                            .unwrap()
                            .write()
                            .expect("surfdata lock")
                            .change_vertex_first_on_s2()
                            .set_point(cpf2);
                    }
                }
                i += 1;
            }
            // The tolerance of points is updated.
            let mut i = 1usize; // OCCT 1-based: for (i = 1; i < last; i++)
            while i < lastn {
                let j = i % len + 1;
                let cursd = seqsurf[i - 1].clone();
                let nextsd = seqsurf[j - 1].clone();
                let (curs1, curs2, mut curp1, mut curp2) = {
                    let g = cursd.read().expect("surfdata lock");
                    let curs1 = if g.is_on_curve1() {
                        g.index_of_c1()
                    } else {
                        g.index_of_s1()
                    };
                    let curs2 = if g.is_on_curve2() {
                        g.index_of_c2()
                    } else {
                        g.index_of_s2()
                    };
                    (
                        curs1,
                        curs2,
                        g.vertex_last_on_s1().clone(),
                        g.vertex_last_on_s2().clone(),
                    )
                };
                let (nexts1, nexts2, mut nextp1, mut nextp2, next_on_c1, next_on_c2) = {
                    let g = nextsd.read().expect("surfdata lock");
                    let nexts1 = if g.is_on_curve1() {
                        g.index_of_c1()
                    } else {
                        g.index_of_s1()
                    };
                    let nexts2 = if g.is_on_curve2() {
                        g.index_of_c2()
                    } else {
                        g.index_of_s2()
                    };
                    (
                        nexts1,
                        nexts2,
                        g.vertex_first_on_s1().clone(),
                        g.vertex_first_on_s2().clone(),
                        g.is_on_curve1(),
                        g.is_on_curve2(),
                    )
                };
                let cur_on_c1 = {
                    let g = cursd.read().expect("surfdata lock");
                    g.is_on_curve1()
                };
                let cur_on_c2 = {
                    let g = cursd.read().expect("surfdata lock");
                    g.is_on_curve2()
                };
                let mut tol1 = curp1.tolerance().max(nextp1.tolerance());
                let mut tol2 = curp2.tolerance().max(nextp2.tolerance());

                if !curp1.is_on_arc() && nextp1.is_on_arc() {
                    // OCCT L3699: curp1 = nextp1 (the whole CommonPoint).
                    curp1 = nextp1.clone();
                    if (curs1 == nexts1) && !next_on_c1 {
                        // Case when it is not possible to pass along the
                        // border without leaving
                        super::chfi3d_builder_2::change_transition(
                            &brep,
                            &nextp1,
                            &mut curp1,
                            nexts1,
                            &dstr,
                        );
                    }
                } else if curp1.is_on_arc() && !nextp1.is_on_arc() {
                    nextp1 = curp1.clone();
                    if (curs1 == nexts1) && !cur_on_c1 {
                        super::chfi3d_builder_2::change_transition(
                            &brep,
                            &curp1,
                            &mut nextp1,
                            curs1,
                            &dstr,
                        );
                    }
                }

                if !curp2.is_on_arc() && nextp2.is_on_arc() {
                    curp2 = nextp2.clone();
                    if (curs2 == nexts2) && !next_on_c2 {
                        super::chfi3d_builder_2::change_transition(
                            &brep,
                            &nextp2,
                            &mut curp2,
                            curs2,
                            &dstr,
                        );
                    }
                } else if curp2.is_on_arc() && !nextp2.is_on_arc() {
                    nextp2 = curp2.clone();
                    if (curs2 == nexts2) && !cur_on_c2 {
                        super::chfi3d_builder_2::change_transition(
                            &brep,
                            &curp2,
                            &mut nextp2,
                            curs2,
                            &dstr,
                        );
                    }
                }

                curp1.set_tolerance(tol1);
                nextp1.set_tolerance(tol1);
                curp2.set_tolerance(tol2);
                nextp2.set_tolerance(tol2);

                let mut b1 = BndBox::default();
                let mut b2 = BndBox::default();
                if curp1.is_on_arc() {
                    chfi3d_enlarge_box_edge_faces(
                        &brep,
                        curp1.arc(),
                        self.my_ef_map.find(curp1.arc()),
                        curp1.parameter_on_arc(),
                        &mut b1,
                    );
                }
                if curp2.is_on_arc() {
                    chfi3d_enlarge_box_edge_faces(
                        &brep,
                        curp2.arc(),
                        self.my_ef_map.find(curp2.arc()),
                        curp2.parameter_on_arc(),
                        &mut b2,
                    );
                }
                // OCCT L3746: bidst is a null stripe handle.
                chfi3d_enlarge_box_dstr(
                    &brep,
                    &mut dstr,
                    None,
                    &cursd.read().expect("surfdata lock"),
                    &mut b1,
                    &mut b2,
                    false,
                );
                chfi3d_enlarge_box_dstr(
                    &brep,
                    &mut dstr,
                    None,
                    &nextsd.read().expect("surfdata lock"),
                    &mut b1,
                    &mut b2,
                    true,
                );
                tol1 = super::chfi3d_builder_2b::chfi3d_box_diag(&b1);
                tol2 = super::chfi3d_builder_2b::chfi3d_box_diag(&b2);
                curp1.set_tolerance(tol1);
                nextp1.set_tolerance(tol1);
                curp2.set_tolerance(tol2);
                nextp2.set_tolerance(tol2);

                // write the mutated CommonPoints back.
                {
                    let mut g = cursd.write().expect("surfdata lock");
                    *g.change_vertex_last_on_s1() = curp1;
                    *g.change_vertex_last_on_s2() = curp2;
                }
                {
                    let mut g = nextsd.write().expect("surfdata lock");
                    *g.change_vertex_first_on_s1() = nextp1;
                    *g.change_vertex_first_on_s2() = nextp2;
                }
                i += 1;
            }
            // The connections edge/new faces are updated.
            for curhels in &ll {
                let wf2 = curhels.firstparam;
                let wl2 = curhels.lastparam;
                let if_;
                let mut il_;
                let mut nwf = wf2;
                let mut nwl = wl2;
                let mut period = 0.0f64;
                let nbed = spine.base().nb_edges() as i64;
                if periodic {
                    period = spine.base().period();
                    nwf = elclib_in_period(wf2, -self.tolesp, period - self.tolesp);
                    if_ = spine.base().index_of_param(nwf, true) as i64;
                    nwl = elclib_in_period(wl2, self.tolesp, period + self.tolesp);
                    il_ = spine.base().index_of_param(nwl, false) as i64;
                    if nwl < nwf + self.tolesp {
                        il_ += nbed;
                    }
                } else {
                    if_ = spine.base().index_of_param(wf2, true) as i64;
                    il_ = spine.base().index_of_param(wl2, false) as i64;
                }
                if if_ == il_ {
                    // fast processing
                    let mut ifloc = if_;
                    if periodic {
                        ifloc = (if_ - 1) % nbed + 1;
                    }
                    let ej = spine.base().edges(ifloc as usize).clone();
                    let mut i = 1usize; // OCCT 1-based
                    while i <= len {
                        let cursd = seqsurf[i - 1].clone();
                        let (fp, lp, surf) = {
                            let g = cursd.read().expect("surfdata lock");
                            (g.first_spine_param(), g.last_spine_param(), g.surf())
                        };
                        if lp < wf2 + self.tolesp || fp > wl2 - self.tolesp {
                            i += 1;
                            continue;
                        }
                        let li = self.my_evi_map.entry(ej.ptr_id()).or_default();
                        li.push(surf);
                        i += 1;
                    }
                } else if if_ < il_ {
                    let mut wv = vec![0.0f64; (il_ - if_) as usize];
                    let mut i = if_;
                    while i < il_ {
                        let mut iloc = i;
                        if periodic {
                            iloc = (i - 1) % nbed + 1;
                        }
                        let mut wi = spine.base().last_parameter_of(iloc as usize);
                        if periodic {
                            wi = elclib_in_period(wi, wf2, wf2 + period);
                        }
                        let pv = spine.base_mut().value_at(wi);
                        // OCCT L3825: Extrema_LocateExtPC(pv, *curhels, wi,
                        // 1.e-8) — the ElSpine curve machinery is pending
                        // (adaptor_curve() == None).
                        wv[(i - if_) as usize] = wi;
                        if let Some(c) = curhels.adaptor_curve() {
                            if let Some(poc) = extrema_locate_ext_pc(
                                pv,
                                &c,
                                wi,
                                f64::NEG_INFINITY,
                                f64::INFINITY,
                                1.0e-8,
                            ) {
                                wv[(i - if_) as usize] = poc.param;
                            }
                        }
                        i += 1;
                    }
                    let mut i = 1usize; // OCCT 1-based
                    while i <= len {
                        let cursd = seqsurf[i - 1].clone();
                        let (fp, lp, surf) = {
                            let g = cursd.read().expect("surfdata lock");
                            (g.first_spine_param(), g.last_spine_param(), g.surf())
                        };
                        let mut jf = 0i64;
                        let mut jl = 0i64;
                        if lp < wf2 + self.tolesp || fp > wl2 - self.tolesp {
                            i += 1;
                            continue;
                        }
                        let mut j = if_;
                        while j < il_ {
                            jf = j;
                            if fp < wv[(j - if_) as usize] - self.tolesp {
                                break;
                            }
                            j += 1;
                        }
                        let mut j = if_;
                        while j < il_ {
                            jl = j;
                            if lp < wv[(j - if_) as usize] + self.tolesp {
                                break;
                            }
                            j += 1;
                        }
                        let mut j = jf;
                        while j <= jl {
                            let mut jloc = j;
                            if periodic {
                                jloc = (j - 1) % nbed + 1;
                            }
                            let ej = spine.base().edges(jloc as usize).clone();
                            let li = self.my_evi_map.entry(ej.ptr_id()).or_default();
                            li.push(surf);
                            j += 1;
                        }
                        i += 1;
                    }
                }
            }
            self.my_ds = Some(dstr);
        }
        {
            let mut st = stripe.write().expect("stripe lock");
            st.my_spine = Some(spine);
        }
        // restore the ElSpine list taken at entry.
        let mut st = stripe.write().expect("stripe lock");
        if let Some(sp) = st.my_spine.as_mut() {
            sp.base_mut().elspines = ll;
        }
    }
}

/// OCCT ChFi3d_Builder_2.cxx L3282-3294 — ChFi3d_BoxDiag.
pub(crate) fn chfi3d_box_diag(bx: &BndBox) -> f64 {
    let (a, b, c, mut d, mut e, mut f) = bx.get();
    d -= a;
    e -= b;
    f -= c;
    d *= d;
    e *= e;
    f *= f;
    (d + e + f).sqrt()
}
