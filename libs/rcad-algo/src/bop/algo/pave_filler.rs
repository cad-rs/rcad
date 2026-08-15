// OCCT BOPAlgo_PaveFiller 鈥?intersection engine.
//
// OCCT BOPAlgo_PaveFiller.cxx / _5.cxx / _6.cxx / _7.cxx
// PerformInternal flow (BOPAlgo_PaveFiller.cxx L235-379):
//
//   Init -> Prepare -> PerformVV -> PerformVE -> UpdatePaveBlocksWithSDVertices
//   -> PerformEE -> UpdatePaveBlocksWithSDVertices
//   -> PerformVF -> UpdatePaveBlocksWithSDVertices
//   -> PerformEF -> UpdatePaveBlocksWithSDVertices -> UpdateInterfsWithSDVertices
//   -> RepeatIntersection -> ForceInterfEE -> ForceInterfEF
//   -> PerformFF -> UpdateBlocksWithSharedVertices -> RefineFaceInfoIn
//   -> MakeSplitEdges -> UpdatePaveBlocksWithSDVertices -> MakeBlocks
//   -> CheckSelfInterference -> UpdateInterfsWithSDVertices -> ReleasePaveBlocks
//   -> RefineFaceInfoOn -> RemoveMicroEdges -> MakePCurves -> ProcessDE

use crate::bop::algo::{Alert, GlueEnum, Report};
use rcad_kernel::core::message::ProgressScope;
use crate::bop::ds::{
    DS, InterferenceVV, InterferenceVE, InterferenceEE,
    InterferenceVF, InterferenceEF, InterferenceFF, BOPDS_Iterator,
};
use crate::bop::ds::pave::{Pave, PaveBlock, SharedPB};
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::int_tools::edge_face::EdgeFace;
use crate::bop::int_tools::common_prt::CommonPrtType;
use crate::bop::int_tools::pnt_on_2_faces::PntOn2S;
use crate::bop::int_tools;
use rcad_kernel::base::proj_lib::project_on_surface;
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use crate::bop::tools::algo_tools;
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::geom::Surface3;
use rcad_kernel::geom::Curve2dEval;
use rcad_kernel::{Curve3, CurveEval};
use rcad_kernel::topods::{self, ShapeType};
use std::collections::{HashSet, HashMap};
use std::sync::Arc;
use glam::DVec3;
use rcad_kernel::topo_shape::{self, Shape};

/// OCCT BOPDS_CoupleOfPaveBlocks 鈥?stores two PBs, new vertex, interference index, tolerance.
struct CoupleOfPBs {
    pb1: SharedPB,
    pb2: SharedPB,
    index_interf: usize,
    tolerance: f64,
    index: usize, // new vertex index (iV)
    point: DVec3, // intersection point (used to fuse close new vertices)
}

use crate::bop::algo::section_attribute::SectionAttribute;

/// OCCT BOPAlgo_PaveFiller::EdgeRangeDistance 鈥?distance from an edge range to a face.
#[derive(Debug, Clone)]
pub(crate) struct EdgeRangeDistance {
    pub(crate) first: f64,
    pub(crate) last: f64,
    pub(crate) distance: f64,
}
impl EdgeRangeDistance {
    fn new(first: f64, last: f64, distance: f64) -> Self { Self { first, last, distance } }
}

/// OCCT BOPAlgo_BPC (PaveFiller_7.cxx L320-365) 鈥?builds pcurve for edge on planar face.
struct BOPAlgo_BPC {
    edge_idx: usize,
    face_idx: usize,
    pcurve: Option<rcad_kernel::geom::Curve2d>,
    to_update: bool,
}

impl BOPAlgo_BPC {
    fn new(edge_idx: usize, face_idx: usize) -> Self {
        BOPAlgo_BPC { edge_idx, face_idx, pcurve: None, to_update: false }
    }
    fn edge_idx(&self) -> usize { self.edge_idx }
    fn face_idx(&self) -> usize { self.face_idx }
    fn pcurve(&self) -> Option<&rcad_kernel::geom::Curve2d> { self.pcurve.as_ref() }
    fn is_to_update(&self) -> bool { self.to_update }
    /// OCCT BOPAlgo_BPC::Perform (PaveFiller_7.cxx L345-353) 鈥?calls
    /// BRepLib::BuildPCurveForEdgeOnPlane (BRepLib_1.cxx L330-338): the pcurve
    /// is built only when the edge has no stored pcurve on the plane face
    /// (BRep_Tool::CurveOnSurface with isStored, BRep_Tool.cxx L327-373);
    /// bToUpdate = !isStored && !aC2D.IsNull().
    fn perform(&mut self, ds: &DS) {
        let edge_shape = ds.shape(self.edge_idx);
        let face_shape = ds.shape(self.face_idx);
        // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the stored
        // pcurve is keyed by (face TShape, L.Predivided(E.Location())).
        let key = ds.face_key(self.face_idx).map(|(fid, floc)| {
            (fid, crate::bop::algo::compose_face_edge_pcurve_location(floc, edge_shape.location, &ds.locations))
        });
        let is_stored = match &*edge_shape.data {
            rcad_kernel::topods::TShape::Edge(ed) => key.map_or(false, |k| ed.pcurves.contains_key(&k)),
            _ => true,
        };
        if is_stored {
            return;
        }
        let face_surf = face_shape.as_face().and_then(|fd| fd.surface.as_ref());
        // OCCT BRep_Tool::Curve(aE) 鈥?the edge's 3D curve with its Location
        // applied (BRepLib_1.cxx L409).
        let edge_curve = ds.edge_curve(self.edge_idx);
        match (edge_curve, face_surf) {
            (Some(curve), Some(surf)) => {
                // OCCT: BRepLib::BuildPCurveForEdgeOnPlane(myE, myF, myCurve, myToUpdate)
                if let Some(pc) = project_on_surface(&curve, surf) {
                    self.pcurve = Some(pc);
                    self.to_update = true;
                }
            }
            _ => {}
        }
    }
}

/// OCCT IntTools_ShrunkRange (IntTools_ShrunkRange.cxx L25-191).
pub(crate) struct ShrunkRange {
    pb: SharedPB,
    n_v1: usize,
    n_v2: usize,
    a_t1: f64,
    a_t2: f64,
    done: bool,
    my_ts1: f64,
    my_ts2: f64,
    is_splittable: bool,
    my_bnd_box: BndBox,
    my_length: f64,
}

impl ShrunkRange {
    pub(crate) fn new(pb: &SharedPB, n_v1: usize, n_v2: usize, a_t1: f64, a_t2: f64) -> Self {
        ShrunkRange {
            pb: pb.clone(), n_v1, n_v2, a_t1, a_t2,
            done: false, my_ts1: a_t1, my_ts2: a_t2, is_splittable: false,
            my_bnd_box: BndBox::new(), my_length: 0.0,
        }
    }
    fn is_done(&self) -> bool { self.done }
    fn is_splittable(&self) -> bool { self.is_splittable }
    pub(crate) fn shrunk_range(&self) -> (f64, f64) { (self.my_ts1, self.my_ts2) }
    fn pave_block(&self) -> &SharedPB { &self.pb }
    fn bnd_box(&self) -> &BndBox { &self.my_bnd_box }

    // OCCT IntTools_ShrunkRange::Perform (IntTools_ShrunkRange.cxx L107-191)
    fn perform(&mut self, ds: &DS) {
        self.done = false;
        self.is_splittable = false;
        // OCCT L113-114: default tolerances
        let a_dtol = rcad_kernel::CONFUSION;    // Precision::Confusion()
        let a_pdtol = rcad_kernel::PCONFUSION;  // Precision::PConfusion()
        // OCCT L117-120: check minimum range
        if self.a_t2 - self.a_t1 < a_pdtol {
            return;
        }
        // OCCT L122-127: read edge/vertex data from DS (replaces TopoDS)
        let n_e = { let r = self.pb.0.read().unwrap(); r.original_edge };
        let curve = match ds.edge_curve(n_e) {
            Some(c) => c.clone(),
            None => return,
        };
        let a_p1 = ds.vertex_point_by_idx(self.n_v1);
        let a_p2 = ds.vertex_point_by_idx(self.n_v2);
        let a_tol_e = ds.edge_tolerance(n_e);
        let mut a_tol_v1 = ds.vertex_tolerance_by_idx(self.n_v1);
        let mut a_tol_v2 = ds.vertex_tolerance_by_idx(self.n_v2);
        // OCCT L129-137: clamp vertex tolerances to edge tolerance + add Confusion
        if a_tol_v1 < a_tol_e {
            a_tol_v1 = a_tol_e;
        }
        if a_tol_v2 < a_tol_e {
            a_tol_v2 = a_tol_e;
        }
        a_tol_v1 += a_dtol;
        a_tol_v2 += a_dtol;
        // OCCT L146-151: compute shrunk range via FindValidRange
        if !self.find_valid_range(&curve, a_tol_e, a_p1, a_tol_v1, a_p2, a_tol_v2) {
            return;
        }
        if self.my_ts2 - self.my_ts1 < a_pdtol {
            return;
        }
        // OCCT L158-175: compute edge length on shrunk range
        // OCCT L162: double aPTolE = aBAC.Resolution(aTolE);
        let mut a_ptol_e = shrunk_range_resolution(&curve, self.my_ts1, self.my_ts2, a_tol_e);
        // OCCT L165: double aPTolEMin = (myT2 - myT1) / 100.;
        let a_ptol_e_min = (self.a_t2 - self.a_t1) * 0.01;
        if a_ptol_e > a_ptol_e_min {
            a_ptol_e = a_ptol_e_min;
        }
        self.my_length = shrunk_range_arc_length(&curve, self.my_ts1, self.my_ts2, a_ptol_e);
        if self.my_length < a_dtol {
            return;
        }
        self.done = true;
        // OCCT L184-187: check splittable
        if self.my_length > (2.0 * a_tol_e + 2.0 * a_dtol) {
            self.is_splittable = true;
        }
        // OCCT L190: BndLib_Add3dCurve::Add(aBAC, myTS1, myTS2, aTolE + aDTol, myBndBox)
        self.my_bnd_box = shrunk_range_bnd_box(&curve, self.my_ts1, self.my_ts2, a_tol_e + a_dtol);
    }

    // OCCT BRepLib::FindValidRange (BRepLib_1.cxx L173-258).
    // Returns theFirst = self.my_ts1, theLast = self.my_ts2.
    pub(crate) fn find_valid_range(
        &mut self,
        curve: &Curve3,
        a_tol_e: f64,
        a_pnt_v1: DVec3,
        a_tol_v1: f64,
        a_pnt_v2: DVec3,
        a_tol_v2: f64,
    ) -> bool {
        find_valid_range_params(
            curve, self.a_t1, self.a_t2, a_tol_e,
            a_pnt_v1, a_tol_v1, a_pnt_v2, a_tol_v2,
            &mut self.my_ts1, &mut self.my_ts2,
        )
    }
}

/// OCCT BRepLib::FindValidRange (BRepLib_1.cxx L173-258) with an explicit
/// parameter range (a_t1, a_t2) instead of a ShrunkRange's PB. Used by
/// ShrunkRange::find_valid_range and by the section-edge micro check
/// (BOPTools_AlgoTools::IsMicroEdge), where the PB has no edge reference.
pub(crate) fn find_valid_range_params(
    curve: &Curve3,
    a_t1: f64,
    a_t2: f64,
    a_tol_e: f64,
    a_pnt_v1: DVec3,
    a_tol_v1: f64,
    a_pnt_v2: DVec3,
    a_tol_v2: f64,
    out_ts1: &mut f64,
    out_ts2: &mut f64,
) -> bool {
    if std::env::var("RCAD_MB_DEBUG").is_ok() {
        eprintln!("[TRACE] find_valid_range_params t=[{:.6},{:.6}]", a_t1, a_t2);
    }
    // OCCT L184: if (theParV2 - theParV1 < Precision::PConfusion()) return false;
    if a_t2 - a_t1 < rcad_kernel::PCONFUSION {
        return false;
    }
    // OCCT L189: bool isInfParV1, isInfParV2
    let is_inf_par_v1 = rcad_kernel::is_infinite_value(a_t1);
    let is_inf_par_v2 = rcad_kernel::is_infinite_value(a_t2);
    // OCCT L191-199: aMaxPar
    let mut a_max_par = 0.0;
    if !is_inf_par_v1 {
        a_max_par = a_t1.abs();
    }
    if !is_inf_par_v2 {
        a_max_par = a_max_par.max(a_t2.abs());
    }
    // OCCT L201-202: anEps
    let an_eps = (shrunk_range_resolution(curve, a_t1, a_t2, a_tol_e) * 0.1)
        .max(epsilon(a_max_par))
        .max(rcad_kernel::PCONFUSION);
    // OCCT L204-225: first endpoint
    if is_inf_par_v1 {
        *out_ts1 = a_t1;
    } else {
        if !find_nearest_valid_point(
            curve, a_t1, a_t2, true,
            a_pnt_v1, a_tol_v1, an_eps,
            out_ts1,
        ) {
            return false;
        }
        if a_t2 - *out_ts1 < an_eps {
            return false;
        }
    }
    // OCCT L227-248: second endpoint
    if is_inf_par_v2 {
        *out_ts2 = a_t2;
    } else {
        if !find_nearest_valid_point(
            curve, a_t1, a_t2, false,
            a_pnt_v2, a_tol_v2, an_eps,
            out_ts2,
        ) {
            return false;
        }
        if *out_ts2 - a_t1 < an_eps {
            return false;
        }
    }
    // OCCT L250-255: check ordering
    if *out_ts1 > *out_ts2 {
        return false;
    }
    true
}

// OCCT Epsilon(double) 鈥?BRepLib_1.cxx L26 鈥?machine-epsilon scaled by value.
fn epsilon(par: f64) -> f64 {
    let eps = 1.1102230246251565e-15; // DBL_EPSILON in double precision
    par.abs() * eps
}

// OCCT theCurve.Resolution(theTol) 鈥?parametric step for given 3D tolerance.
// Approximated by sampling derivative magnitude on the sub-range [t1, t2].
pub(crate) fn shrunk_range_resolution(curve: &Curve3, t1: f64, t2: f64, tol: f64) -> f64 {
    let n_samples = 100;
    let dt = (t2 - t1) / n_samples as f64;
    if dt <= 0.0 { return 1e-3; }
    let mut max_speed = 1e-12;
    for i in 0..=n_samples {
        let t = t1 + dt * i as f64;
        let d = curve.derivative_at(t);
        let speed = d.length();
        if speed > max_speed { max_speed = speed; }
    }
    (tol / max_speed).max(rcad_kernel::PCONFUSION)
}

// OCCT GCPnts_AbscissaPoint::Length 鈥?adaptive arc length.
// Simplified: Simpson integration with tolerance-based subdivision.
// Depth-limited so a pathological range (e.g. a section edge spanning many
// periods) converges instead of overflowing the stack, mirroring the bounded
// iteration of OCCT's AbscissaPoint.
fn shrunk_range_arc_length(curve: &Curve3, t1: f64, t2: f64, tol: f64) -> f64 {
    if (t2 - t1).abs() < 1e-15 { return 0.0; }
    // Simple adaptive Simpson (depth-limited).
    fn simpson_step(curve: &Curve3, a: f64, b: f64, fa: f64, fb: f64, fm: f64, tol: f64, depth: u32) -> f64 {
        let m = (a + b) * 0.5;
        let h = (b - a) * 0.5;
        let fm1 = curve.derivative_at(a + h * 0.5).length();
        let fm2 = curve.derivative_at(m + h * 0.5).length();
        let s1 = h / 3.0 * (fa + 4.0 * fm + fb);
        let s2 = h / 6.0 * (fa + 4.0 * fm1 + 2.0 * fm + 4.0 * fm2 + fb);
        if (s2 - s1).abs() < tol || depth >= 24 {
            s2
        } else {
            let left = simpson_step(curve, a, m, fa, fm, fm1, tol * 0.5, depth + 1);
            let right = simpson_step(curve, m, b, fm, fb, fm2, tol * 0.5, depth + 1);
            left + right
        }
    }
    let fa = curve.derivative_at(t1).length();
    let fb = curve.derivative_at(t2).length();
    let fm = curve.derivative_at((t1 + t2) * 0.5).length();
    simpson_step(curve, t1, t2, fa, fb, fm, tol, 0)
}

// OCCT BndLib_Add3dCurve::Add 鈥?build bounding box for curve subrange.
fn shrunk_range_bnd_box(curve: &Curve3, t1: f64, t2: f64, tol: f64) -> BndBox {
    let n_samples = 100;
    let dt = (t2 - t1) / n_samples as f64;
    if dt <= 0.0 {
        let p = curve.point_at(t1);
        let mut b = BndBox::from_point(p);
        b.set_gap(tol);
        return b;
    }
    let mut b = BndBox::new();
    let mut first = true;
    for i in 0..=n_samples {
        let t = t1 + dt * i as f64;
        let p = curve.point_at(t);
        if first {
            // BndBox from_point has no void flag and uses the point
            b = BndBox::from_point(p);
            first = false;
        } else {
            b.add_point(p);
        }
    }
    b.set_gap(tol);
    b
}

// OCCT BRepLib::findNearestValidPoint (BRepLib_1.cxx L31-168).
// Walk from one endpoint until the curve exits the vertex tolerance sphere,
// then binary search for the precise boundary.
fn find_nearest_valid_point(
    curve: &Curve3,
    the_first: f64,
    the_last: f64,
    is_first: bool,
    the_vert_pnt: DVec3,
    the_tol: f64,
    the_eps: f64,
    the_par: &mut f64,
) -> bool {
    // OCCT L42-47: start from the appropriate end
    let (a_start_u, an_end_u) = if is_first {
        (the_first, the_last)
    } else {
        (the_last, the_first)
    };
    // OCCT L48-54: check that the needed end is inside the sphere
    let a_p = curve.point_at(a_start_u);
    let a_sq_tol = the_tol * the_tol;
    if a_p.distance_squared(the_vert_pnt) > a_sq_tol {
        return false; // vertex does not cover this end of the curve
    }
    // OCCT L58-65: general step = Resolution(theTol) * 1.01
    let mut a_step = shrunk_range_resolution(curve, a_start_u, an_end_u, the_tol) * 1.01;
    if a_step < the_eps {
        a_step = the_eps;
    }
    // OCCT L66-82: aD1Mag for Bezier/BSpline singularity acceleration
    let a_d1_mag = match curve {
        Curve3::BSpline(..) | Curve3::Bezier(..) => {
            // OCCT: 1. / theCurve.Resolution(1.) * 0.01, squared
            let r = shrunk_range_resolution(curve, a_start_u, an_end_u, 1.0);
            let v = 1.0 / r * 0.01;
            v * v
        }
        Curve3::Offset(off) => {
            // OCCT unwraps offset once to check base curve type
            match &*off.basis {
                Curve3::BSpline(..) | Curve3::Bezier(..) => {
                    let r = shrunk_range_resolution(curve, a_start_u, an_end_u, 1.0);
                    let v = 1.0 / r * 0.01;
                    v * v
                }
                _ => 0.0,
            }
        }
        _ => 0.0,
    };
    // OCCT L83-86: reverse step direction when walking from second vertex
    if !is_first {
        a_step = -a_step;
    }
    // OCCT L87-147: walk until out of sphere
    let mut is_out = false;
    let mut an_u_in = a_start_u;
    let mut an_u_out = an_u_in;
    while !is_out {
        an_u_in = an_u_out;
        an_u_out += a_step;
        // OCCT L94-107: check bounds
        if (is_first && an_u_out > an_end_u) || (!is_first && an_u_out < an_end_u) {
            let a_p_end = curve.point_at(an_end_u);
            is_out = a_p_end.distance_squared(the_vert_pnt) > a_sq_tol;
            if !is_out {
                return false; // all range inside sphere
            }
            an_u_out = an_end_u;
            break;
        }
        // OCCT L108-137: handle BSpline/Bezier derivative singularity
        if a_d1_mag > 0.0 {
            let mut a_step_local = a_step;
            loop {
                let a_d1 = curve.derivative_at(an_u_out);
                let a_p_check = curve.point_at(an_u_out);
                is_out = a_p_check.distance_squared(the_vert_pnt) > a_sq_tol;
                if !is_out && a_d1.length_squared() < a_d1_mag {
                    a_step_local *= 2.0;
                    an_u_out += a_step_local;
                    if (is_first && an_u_out < an_end_u) || (!is_first && an_u_out > an_end_u) {
                        continue; // still in range
                    }
                    // out of range, check endpoint
                    an_u_out = an_end_u;
                    let a_p_end = curve.point_at(an_u_out);
                    is_out = a_p_end.distance_squared(the_vert_pnt) > a_sq_tol;
                    if !is_out {
                        return false;
                    }
                }
                break;
            }
        } else {
            let a_p_check = curve.point_at(an_u_out);
            is_out = a_p_check.distance_squared(the_vert_pnt) > a_sq_tol;
        }
    }
    // OCCT L149-168: precise solution with binary search
    let mut a_delta = (an_u_out - an_u_in).abs();
    while a_delta > the_eps {
        let a_mid_u = (an_u_in + an_u_out) * 0.5;
        let a_p_mid = curve.point_at(a_mid_u);
        is_out = a_p_mid.distance_squared(the_vert_pnt) > a_sq_tol;
        if is_out {
            an_u_out = a_mid_u;
        } else {
            an_u_in = a_mid_u;
        }
        a_delta = (an_u_out - an_u_in).abs();
    }
    *the_par = (an_u_in + an_u_out) * 0.5;
    true
}

pub struct PaveFiller {
    // 鈹€鈹€ BOPAlgo_Algo base (inherited) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    pub(crate) my_report: Report,                 // BOPAlgo_Algo::myReport
    pub(crate) my_run_parallel: bool,             // BOPAlgo_Algo::myRunParallel
    pub(crate) my_fuzzy_value: f64,               // BOPAlgo_Algo::myFuzzyValue
    // 鈹€鈹€ BOPAlgo_PaveFiller members 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // OCCT BOPAlgo_PaveFiller.hxx L639-652:
    pub(crate) ds: DS,                            // L640: myDS (owned, OCCT: heap-allocated)
    pub(crate) my_iterator: Option<Box<BOPDS_Iterator>>, // L641: myIterator
    pub(crate) my_context: IntToolsContext,        // BOPAlgo_PaveFiller::myContext (L642)
    pub(crate) my_glue: GlueEnum,                 // BOPAlgo_PaveFiller::myGlue (L647)
    pub(crate) my_section_attribute: SectionAttribute, // BOPAlgo_PaveFiller::mySectionAttribute (L643)
    pub(crate) my_non_destructive: bool,          // BOPAlgo_PaveFiller::myNonDestructive (L644)
    pub(crate) my_is_primary: bool,               // BOPAlgo_PaveFiller::myIsPrimary (L645)
    pub(crate) my_avoid_build_pcurve: bool,       // BOPAlgo_PaveFiller::myAvoidBuildPCurve (L646)
    pub(crate) my_arguments: Vec<topo_shape::Shape>, // BOPAlgo_PaveFiller::myArguments (L639)
    pub(crate) my_fpb_done: HashMap<usize, HashSet<u64>>, // BOPAlgo_PaveFiller::myFPBDone (L650)
    pub(crate) my_increased_ss: HashSet<usize>,   // BOPAlgo_PaveFiller::myIncreasedSS (L651)
    pub(crate) my_verts_to_avoid_extension: HashSet<usize>, // BOPAlgo_PaveFiller::myVertsToAvoidExtension (L652)
    // OCCT L657-659: NCollection_DataMap<BOPDS_Pair, List<EdgeRangeDistance>> myDistances
    // rcad: HashMap keyed by (edge1, edge2) pair
    pub(crate) my_distances: HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    pub stop_after: Option<&'static str>,
}

impl PaveFiller {
    pub fn new() -> Self {
        PaveFiller {
            ds: DS::new(),
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: 0.0,
            my_iterator: None,
            my_context: IntToolsContext::new(),
            my_glue: GlueEnum::GlueOff,
            my_section_attribute: SectionAttribute::default(),
            my_non_destructive: false,
            my_is_primary: true,
            my_avoid_build_pcurve: false,
            my_arguments: Vec::new(),
            my_fpb_done: HashMap::new(),
            my_increased_ss: HashSet::new(),
            my_verts_to_avoid_extension: HashSet::new(),
            my_distances: HashMap::new(),
            stop_after: None,
        }
    }

    /// Create with a pre-configured owned DS.
    pub fn new_with_ds(ds: DS) -> Self {
        PaveFiller {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: 0.0,
            my_iterator: None,
            my_context: IntToolsContext::new(),
            my_glue: GlueEnum::GlueOff,
            my_section_attribute: SectionAttribute::default(),
            my_non_destructive: false,
            my_is_primary: true,
            my_avoid_build_pcurve: false,
            my_arguments: Vec::new(),
            my_fpb_done: HashMap::new(),
            my_increased_ss: HashSet::new(),
            my_verts_to_avoid_extension: HashSet::new(),
            my_distances: HashMap::new(),
            stop_after: None,
        }
    }

    pub fn set_arguments(&mut self, args: Vec<Shape>) {
        self.my_arguments = args;
    }

    /// OCCT BOPAlgo_PaveFiller::Clear (PaveFiller.cxx L136-141).
    /// Clears internal state (iterator, data maps).
    pub fn clear(&mut self) {
        // OCCT L137: BOPAlgo_Algo::Clear() 鈥?clears report
        self.my_report.clear();
        // OCCT L138-139: delete myIterator; myIterator = nullptr;
        self.my_iterator = None;
        // OCCT L141: myIncreasedSS.Clear();
        self.my_increased_ss.clear();
        // Note: myDS is borrowed (not owned), so not deleted.
    }
    pub fn set_glue(&mut self, enable: bool, tolerance: f64) {
        self.my_glue = if enable { GlueEnum::GlueFull } else { GlueEnum::GlueOff };
        self.my_fuzzy_value = tolerance;
    }
    pub fn fuzzy_value(&self) -> f64 { self.my_fuzzy_value }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.my_fuzzy_value = v; }
    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }
    pub fn report(&self) -> &Report { &self.my_report }
    pub fn ds(&self) -> &DS { &self.ds }
    pub fn ds_mut(&mut self) -> &mut DS { &mut self.ds }

    /// If stop_after matches stage, return true (caller should return).
    fn check_stop(&self, stage: &'static str) -> bool {
        self.stop_after.map_or(false, |s| s == stage)
    }

    /// OCCT BOPAlgo_PaveFiller::Perform (PaveFiller.cxx L218-232).
    pub fn perform(&mut self, the_range: &ProgressScope) {
        self.perform_internal(the_range);
    }

    /// OCCT BOPAlgo_PaveFiller::PerformInternal (PaveFiller.cxx L235-379).
    pub(crate) fn perform_internal(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L239-244: Message_ProgressScope aPS(theRange, "Performing intersection of shapes", 100)
        let a_ps = the_range.sub_scope("Performing intersection of shapes", 100);

        // OCCT L247: Init(aPS.Next(5));
        self.init(&a_ps.sub_scope("Init", 5));
        if self.has_errors() { return; }
        if self.check_stop("after_Init") { return; }

        // OCCT L258: Prepare(aPS.Next(aSteps.GetStep(PIOperation_Prepare)));
        self.prepare(&a_ps.sub_scope("Prepare", 10));
        if self.has_errors() { return; }
        if self.check_stop("after_Prepare") { return; }

        // OCCT: PerformVV(aPS.Next(...))
        self.perform_vv(&a_ps.sub_scope("Perform VV", 8));
        if self.has_errors() { return; }
        if self.check_stop("after_PerformVV") { return; }

        // OCCT: PerformVE(aPS.Next(...))
        self.perform_ve(&a_ps.sub_scope("Perform VE", 8));
        if self.has_errors() { return; }
        if self.check_stop("after_PerformVE") { return; }

        self.update_pave_blocks_with_sd_vertices();

        // OCCT: PerformEE(aPS.Next(...))
        self.perform_ee(&a_ps.sub_scope("Perform EE", 10));
        if self.has_errors() { return; }
        if self.check_stop("after_PerformEE") { return; }

        self.update_pave_blocks_with_sd_vertices();

        // OCCT: PerformVF(aPS.Next(...))
        self.perform_vf(&a_ps.sub_scope("Perform VF", 5));
        if self.has_errors() { return; }
        if self.check_stop("after_PerformVF") { return; }

        self.update_pave_blocks_with_sd_vertices();

        // OCCT: PerformEF(aPS.Next(...))
        self.perform_ef(&a_ps.sub_scope("Perform EF", 10));
        if self.has_errors() { return; }
        if self.check_stop("after_PerformEF") { return; }

        self.update_pave_blocks_with_sd_vertices();
        self.update_interfs_with_sd_vertices();

        // OCCT: RepeatIntersection(aPS.Next(...))
        self.repeat_intersection(&a_ps.sub_scope("Repeat intersection", 5));
        if self.has_errors() { return; }
        if self.check_stop("after_RepeatIntersection") { return; }

        // OCCT: ForceInterfEE(aPS.Next(...))
        self.force_interf_ee(&a_ps.sub_scope("Force EE", 3));
        if self.has_errors() { return; }
        if self.check_stop("after_ForceInterfEE") { return; }

        // OCCT: ForceInterfEF(aPS.Next(...))
        self.force_interf_ef(&a_ps.sub_scope("Force EF", 3));
        if self.has_errors() { return; }
        if self.check_stop("after_ForceInterfEF") { return; }

        // OCCT: PerformFF(aPS.Next(...))
        self.perform_ff(&a_ps.sub_scope("Perform FF", 12));
        if self.has_errors() { return; }
        if self.check_stop("after_PerformFF") { return; }

        self.update_blocks_with_shared_vertices();
        self.refine_face_info_in();

        // OCCT: MakeSplitEdges(aPS.Next(...))
        self.make_split_edges(&a_ps.sub_scope("Make split edges", 6));
        if self.has_errors() { return; }
        if self.check_stop("after_MakeSplitEdges") { return; }

        self.update_pave_blocks_with_sd_vertices();

        // OCCT: MakeBlocks(aPS.Next(...))
        self.make_blocks(&a_ps.sub_scope("Make blocks", 6));
        if self.has_errors() { return; }
        if self.check_stop("after_MakeBlocks") { return; }

        self.check_self_interference();
        self.update_interfs_with_sd_vertices();
        self.ds.release_pave_blocks();
        self.refine_face_info_on();
        self.remove_micro_edges();

        // OCCT: MakePCurves(aPS.Next(...))
        self.make_pcurves(&a_ps.sub_scope("Make pcurves", 5));
        if self.has_errors() { return; }
        if self.check_stop("after_MakePCurves") { return; }

        // OCCT: ProcessDE(aPS.Next(...))
        self.process_de(&a_ps.sub_scope("Process DE", 4));
        if self.has_errors() { return; }
        if self.check_stop("after_ProcessDE") { return; }
    }

    // ====================================================================
    // VV 鈥?OCCT BOPAlgo_PaveFiller_1.cxx L45-132
    // ====================================================================
    fn perform_vv(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L47-48: int n1, n2, iFlag, aSize; handle<Allocator> aAllocator
        // OCCT L50: myIterator->Initialize(TopAbs_VERTEX, TopAbs_VERTEX)
        // rcad: initialize() then copy pair list. Rust borrow checker prevents
        // holding a mutable borrow on my_iterator while accesing self.ds (in C++,
        // myDS and myIterator are independently accessible member variables).
        let my_iterator = match &mut self.my_iterator {
            Some(it) => it,
            None => return,
        };
        my_iterator.initialize(ShapeType::Vertex, ShapeType::Vertex);
        let pairs: Vec<(usize, usize)> = my_iterator.pairs(ShapeType::Vertex, ShapeType::Vertex).to_vec();
        let a_size = pairs.len();
        if a_size == 0 {
            return;
        }
        // OCCT L58-59: aVVs.SetIncrement(aSize) 鈥?Rust Vec auto-grows


        // OCCT L62-63: NCollection_IndexedDataMap<int, NCollection_List<int>> aMILI(100, aAllocator);
        //             NCollection_List<NCollection_List<int>> aMBlocks(aAllocator);
        // OCCT's IndexedDataMap preserves insertion order (keys iterated in
        // insertion order by MakeBlocks); rcad mirrors it with an IndexMap.
        let mut a_mili: indexmap::IndexMap<usize, Vec<usize>> =
            indexmap::IndexMap::new();
        let mut a_mblocks: Vec<Vec<usize>> = Vec::new();

        // 1. Map V/LV (OCCT L69-98)
        for &(n1, n2) in &pairs {
            if the_range.user_break() { return; }
            // OCCT L77-81: if already interfering, connect and continue
            if self.ds.has_interf(n1, n2) {
                algo_tools::fill_map(n1, n2, &mut a_mili);
                continue;
            }

            // OCCT L84-88: resolve SD vertices
            let mut n1sd = n1;
            self.ds.has_shape_sd(n1, &mut n1sd);
            let mut n2sd = n2;
            self.ds.has_shape_sd(n2, &mut n2sd);

            // OCCT L90-91: get vertex from DS
            let v1_tol = self.ds.vertex_tolerance_by_idx(n1sd);
            let v1_pnt = self.ds.vertex_point_by_idx(n1sd);
            let v2_tol = self.ds.vertex_tolerance_by_idx(n2sd);
            let v2_pnt = self.ds.vertex_point_by_idx(n2sd);

            // OCCT L94: iFlag = BOPTools_AlgoTools::ComputeVV(aV1, aV2, myFuzzyValue);
            let i_flag = crate::bop::tools::algo_tools::compute_vv(
                v1_tol, v1_pnt, v2_tol, v2_pnt, self.my_fuzzy_value);

            // OCCT L94-97: if (!iFlag) { FillMap(n1, n2, aMILI, aAllocator); }
            if i_flag == 0 {
                algo_tools::fill_map(n1, n2, &mut a_mili);
            }
        }

        // OCCT L101: BOPAlgo_Tools::MakeBlocks(aMILI, aMBlocks, aAllocator);
        algo_tools::make_blocks(&a_mili, &mut a_mblocks);

        // OCCT L104-113: MakeSDVertices for each block
        for a_li in &a_mblocks {
            if the_range.user_break() { return; }
            self.make_sd_vertices_vv(a_li, true);
        }

        // OCCT L115-127: InitPaveBlocksForVertex for each SD key
        let sd_keys: Vec<usize> = self.ds.shapes_sd.keys().copied().collect();
        for &n1 in &sd_keys {
            if the_range.user_break() { return; }
            self.ds.init_pave_blocks_for_vertex(n1);
        }
    }

    /// OCCT BOPAlgo_PaveFiller::MakeSDVertices (PaveFiller_1.cxx L136-233).
    /// Returns the new/updated SD vertex index.
    fn make_sd_vertices_vv(&mut self, the_vert_indices: &[usize], the_add_interfs: bool) -> usize {
        // OCCT L139-140: TopoDS_Vertex aVSD, aVn; int nSD = -1;
        let mut n_sd = usize::MAX; // OCCT: nSD = -1
        // OCCT L141-161: collect shapes into aLV, tracking existing SD
        // rcad: build list of (point, tolerance) pairs 鈥?no TopoDS wrapper in DS
        let mut a_lv_points: Vec<(DVec3, f64)> = Vec::new();

        for &n_x in the_vert_indices {
            // OCCT L146: if (myDS->HasShapeSD(nX, nSD1))
            let mut n_sd1 = usize::MAX;
            if self.ds.has_shape_sd(n_x, &mut n_sd1) {
                // OCCT L149-152: if (nSD == -1) { aVSD = aVSD1; nSD = nSD1; }
                if n_sd == usize::MAX {
                    n_sd = n_sd1;
                }
            }
            // OCCT L159-160: const TopoDS_Shape& aV = myDS->Shape(nX); aLV.Append(aV);
            let p = self.ds.vertex_point_by_idx(n_x);
            let t = self.ds.vertex_tolerance_by_idx(n_x);
            a_lv_points.push((p, t));
        }

        // OCCT L162: BOPTools_AlgoTools::MakeVertex(aLV, aVn);
        let (centroid, max_tol) = crate::bop::tools::algo_tools::make_vertex(&a_lv_points);

        // OCCT L163-180: if (nSD != -1) update existing SD else create new
        let n_v;
        if n_sd != usize::MAX {
            // OCCT L167-169: update existing SD vertex's point and tolerance.
            // In-place (OCCT UpdateVertex) to keep the vertex identity.
            self.ds.mutate_shape_data(n_sd, |ts| {
                if let rcad_kernel::topods::TShape::Vertex(vd) = ts {
                    vd.point = centroid;
                    vd.tolerance = max_tol;
                }
            });
            self.ds.remap_shape_idx(n_sd);
            n_v = n_sd;
        } else {
            // OCCT L175-180: Append new vertex to DS
            n_v = self.ds.push_vertex(centroid, max_tol);
        }

        // OCCT L181-184: update bounding box for the SD vertex
        // OCCT: aBox.Add(BRep_Tool::Pnt(aVn)); aBox.SetGap(Tolerance(aVn) + Precision::Confusion());
        {
            let vt = max_tol + rcad_kernel::CONFUSION;
            let si = self.ds.change_shape_info(n_v);
            si.bbox = BndBox::from_point(centroid);
            si.bbox.set_gap(vt);
        }

        // OCCT L186-191: get InterfVV array, pre-allocate if theAddInterfs
        // rcad: Vec auto-extends; no pre-alloc needed.

        // OCCT L193-231: AddShapeSD + self-interference warning + VV interferences
        for i in 0..the_vert_indices.len() {
            let n1 = the_vert_indices[i];
            // OCCT L197: myDS->AddShapeSD(n1, nV);
            self.ds.add_shape_sd(n1, n_v);
            // OCCT L199: int iR1 = myDS->Rank(n1);
            let i_r1 = self.ds.rank(n1);

            // OCCT L202-203: List::Iterator aItLI2 = aItLI; aItLI2.Next();
            for j in (i + 1)..the_vert_indices.len() {
                let n2 = the_vert_indices[j];
                // OCCT L208-218: self-interference warning for same-rank vertices
                // OCCT creates TopoDS_Compound for the warning; rcad stores indices.
                if i_r1 >= 0 && i_r1 == self.ds.rank(n2) {
                    self.my_report.add_warning(
                        Alert::SelfInterferingShape(vec![n1, n2]));
                }
                // OCCT L221-229: add VV interference
                if the_add_interfs {
                    if self.ds.add_interf(n1, n2) {
                        self.ds.interf_vv.push(InterferenceVV {
                            v1: n1, v2: n2, merged_vertex: n_v,
                        });
                    }
                }
            }
        }
        n_v
    }

    // ====================================================================
    // VE 鈥?OCCT BOPAlgo_PaveFiller_2.cxx L171-238
    // ====================================================================
    fn perform_ve(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L173: FillShrunkData(TopAbs_VERTEX, TopAbs_EDGE)
        self.fill_shrunk_data(ShapeType::Vertex, ShapeType::Edge);

        // OCCT L175: myIterator->Initialize(TopAbs_VERTEX, TopAbs_EDGE)
        // rcad: initialize then copy pairs (borrow checker limitation, see perform_vv)
        let my_iterator = match &mut self.my_iterator {
            Some(it) => it,
            None => return,
        };
        my_iterator.initialize(ShapeType::Vertex, ShapeType::Edge);
        let pairs: Vec<(usize, usize)> = my_iterator.pairs(ShapeType::Vertex, ShapeType::Edge).to_vec();
        let i_size = pairs.len();
        for &(a, b) in &pairs {
            let (v, e) = if self.ds.shapes[a].shape_type == ShapeType::Vertex { (a, b) } else { (b, a) };
            let p = self.ds.vertex_point_by_idx(v);
            let _ = (v, e, p);
        }
        if i_size == 0 {
            return;
        }

        // OCCT L185: NCollection_IndexedDataMap<handle<PaveBlock>, NCollection_List<int>> aMVEPairs
        let mut a_mve_pairs: indexmap::IndexMap<u64, (SharedPB, Vec<usize>)> =
            indexmap::IndexMap::new();

        // OCCT L186-235: iterate pairs
        for &(a, b) in &pairs {
            // OCCT myIterator->Value(nV, nE) returns (vertex, edge) in correct order.
            // rcad pairs() returns (min, max) 鈥?identify Vertex/Edge by shape type.
            let (n_v, n_e) = if self.ds.shapes[a].shape_type == ShapeType::Vertex {
                (a, b)
            } else {
                (b, a)
            };
            if the_range.user_break() { return; }
            // OCCT L195-199: if (aSIE.HasSubShape(nV)) continue;
            if self.ds.shapes[n_e].has_sub_shape(n_v) { continue; }
            // OCCT L201-204: if (aSIE.HasFlag()) continue;
            if self.ds.shapes[n_e].has_flag() { continue; }
            // OCCT L206-209: if (myDS->HasInterf(nV, nE)) continue;
            if self.ds.has_interf(n_v, n_e) { continue; }
            // OCCT L211-214: if (myDS->HasInterfShapeSubShapes(nV, nE)) continue;
            // OCCT default third param is true (any sub-shape interferes).
            if self.ds.has_interf_shape_sub_shapes(n_v, n_e, true) { continue; }

            // OCCT L216-220: const List<...>& aLPB = myDS->PaveBlocks(nE);
            let a_lpb: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e).to_vec();
            if a_lpb.is_empty() { continue; }

            // OCCT L222-227: const handle<PaveBlock>& aPB = aLPB.First(); IsSplittable?
            let a_pb = a_lpb[0].clone();
            if !a_pb.0.read().unwrap().is_splittable() { continue; }

            // OCCT L229-234: add vertex to list keyed by PB
            let pb_ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
            let entry = a_mve_pairs.entry(pb_ptr).or_insert((a_pb, Vec::new()));
            entry.1.push(n_v);
        }

        // OCCT L237: IntersectVE(aMVEPairs, aPS.Next())
        self.intersect_ve(&a_mve_pairs, true);
    }

    /// OCCT BOPAlgo_PaveFiller::IntersectVE (PaveFiller_2.cxx L242-394).
    fn intersect_ve(
        &mut self,
        the_ve_pairs: &indexmap::IndexMap<u64, (SharedPB, Vec<usize>)>,
        the_add_interfs: bool,
    ) {
        let a_nb_ve = the_ve_pairs.len();
        if a_nb_ve == 0 {
            return;
        }

        // OCCT L253-257: aVEs.SetIncrement(aNbVE) 鈥?Rust Vec auto-grows

        // OCCT L260-265: aVVE, aDMVSD declarations
        // OCCT L267-322: build solver objects
        // OCCT aMEdges (L276) is NCollection_Map<int> 鈥?iteration order feeds
        // SplitPaveBlocks (theMEdges), which creates SD vertices/edges in that
        // order; model the OCCT map exactly.
        let mut a_m_edges: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();

        for (_pb_ptr, (a_pb, verts)) in the_ve_pairs {
            let n_e = a_pb.0.read().unwrap().original_edge;

            // OCCT L278-284: build set of vertex indices of all PBs of this edge
            let mut a_mvpb: std::collections::HashSet<usize> = std::collections::HashSet::new();
            let all_pbs = self.ds.edge_pave_blocks(n_e);
            for pb in all_pbs {
                let (v1, v2) = { let r = pb.0.read().unwrap(); r.indices() };
                a_mvpb.insert(v1);
                a_mvpb.insert(v2);
            }

            // OCCT L288-321: for each vertex in the list
            // OCCT aDMVSD (PaveFiller_2.cxx L265) is NCollection_DataMap<BOPDS_Pair,
            // List<int>>; the per-(nVSD, nE) solver sequence is aVVE, a vector
            // iterated in insertion order (L330-410) 鈥?an IndexMap reproduces
            // that order.
            let mut a_dmvsd: indexmap::IndexMap<(usize, usize), Vec<usize>> =
                indexmap::IndexMap::new();

            for &n_v in verts {
                // OCCT L292-296: resolve SD
                let mut n_vsd = n_v;
                self.ds.has_shape_sd(n_v, &mut n_vsd);

                // OCCT L298-300: skip if already endpoint of a PB of this edge
                if a_mvpb.contains(&n_vsd) {
                    continue;
                }

                // OCCT L302-310: dedup by (nVSD, nE) pair
                let a_pair = (n_vsd, n_e);
                let entry = a_dmvsd.entry(a_pair).or_default();
                entry.push(n_v);
            }

            // OCCT L324-332: run intersection for each unique (nVSD, nE) pair
            for ((n_vsd, _n_e), orig_verts) in &a_dmvsd {
                // OCCT: myContext->ComputeVE(aV, aE, aT, aTolVNew, myFuzzyValue)
                let (i_flag, a_t, a_tol_v_new) =
                    self.my_context.compute_ve(*n_vsd, n_e, &self.ds, self.my_fuzzy_value);
                if std::env::var("RCAD_EE_DEBUG").is_ok() {
                    eprintln!("[EE-DBG] intersect_ve edge={} v={} flag={} t={:.9} tol={:.6}", n_e, *n_vsd, i_flag, a_t, a_tol_v_new);
                }
                // OCCT L352: if (Flag() != 0) { if (HasErrors()) AddWarning; continue; }
                if i_flag != 0 { continue; }

                // OCCT L368: int nVx = UpdateVertex(nV, aTolVNew);
                let n_vx = self.update_vertex(*n_vsd, a_tol_v_new);

                // OCCT L371-388: find the PB whose range contains the parameter
                let all_pbs_for_edge = self.ds.edge_pave_blocks(n_e);
                let mut found_pb: Option<SharedPB> = None;
                for pb in all_pbs_for_edge {
                    let (t1, t2) = { let r = pb.0.read().unwrap(); r.range() };
                    if a_t > t1 && a_t < t2 {
                        found_pb = Some(pb.clone());
                        break;
                    }
                }
                let Some(target_pb) = found_pb else { continue; };

                // OCCT L390-393: aPB->AppendExtPave(aPave)
                {
                    let mut pbw = target_pb.0.write().unwrap();
                    pbw.append_ext_pave(Pave { vertex_idx: n_vx, param: a_t });
                }

                // OCCT L396-419: add VE interference for EACH original vertex
                // OCCT: aVEs.Appended() always called, then AddInterf separately
                if the_add_interfs {
                    // OCCT L412: if (myDS->IsNewShape(nVx)) 鈥?nVx from UpdateVertex, not loop var!
                    let resolved_is_new = self.ds.is_new_shape(n_vx);
                    for &n_vx_old in orig_verts {
                        // OCCT L406-408: aVE = aVEs.Appended(); SetIndices; SetParameter
                        let mut ve = InterferenceVE {
                            vertex: n_vx_old,
                            edge: n_e,
                            param: a_t,
                            index_new: 0,
                        };
                        // OCCT L412-415: SetIndexNew only if IsNewShape(nVx)
                        if resolved_is_new {
                            ve.index_new = n_vx;
                        }
                        self.ds.interf_ve.push(ve);
                        // OCCT L410: myDS->AddInterf(nVOld, nE) 鈥?called unconditionally
                        self.ds.add_interf(n_vx_old, n_e);
                    }
                }

                a_m_edges.add(n_e);
            }
        }

        // OCCT L420-424: SplitPaveBlocks(aMEdges, theAddInterfs)
        if !a_m_edges.is_empty() {
            self.split_pave_blocks(&a_m_edges, the_add_interfs);
        }
    }

    // ====================================================================
    // EE 鈥?OCCT BOPAlgo_PaveFiller.cxx L279-286 + PaveFiller_5.cxx L145-590
    // ====================================================================
    fn perform_ee(&mut self, the_range: &ProgressScope) {
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            eprintln!("[EE-DBG] === perform_ee START stop_after={:?}", self.stop_after);
        }
        if the_range.user_break() { return; }
        // OCCT L147: FillShrunkData(TopAbs_EDGE, TopAbs_EDGE)
        self.fill_shrunk_data(ShapeType::Edge, ShapeType::Edge);

        // OCCT L149-151: myIterator->Initialize(EDGE, EDGE); iSize
        let my_iterator = match &mut self.my_iterator {
            Some(it) => it,
            None => return,
        };
        my_iterator.initialize(ShapeType::Edge, ShapeType::Edge);
        let i_size = my_iterator.pairs(ShapeType::Edge, ShapeType::Edge).len();
        // OCCT L152-155: if (!iSize) return;
        if i_size == 0 {
            return;
        }
        // OCCT L178-179: aEEs.SetIncrement(iSize) 鈥?Rust Vec auto-grows

        let mut new_ee: Vec<InterferenceEE> = Vec::new();
        let mut cb_pairs: Vec<(SharedPB, SharedPB)> = Vec::new();
        let mut mvcpb: Vec<CoupleOfPBs> = Vec::new();
        // OCCT: aEEs = myDS->InterfEE(); interferences appended directly to DS array
        // rcad: local new_ee, extended later. idx_interf must be absolute DS index.
        let ee_base = self.ds.interf_ee.len();
        // OCCT L167: NCollection_Map<int> aMEdges;
        let mut a_m_edges: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();

        // OCCT L181-267: iterate pairs
        // rcad: copy pairs (borrow checker)
        let ee_pairs: Vec<(usize, usize)> = my_iterator.pairs(ShapeType::Edge, ShapeType::Edge).to_vec();
        let ee_debug = std::env::var("RCAD_EE_DEBUG").is_ok();
        for &(n_e1, n_e2) in &ee_pairs {
            if ee_debug {
                let tn = |c: Option<rcad_kernel::geom::Curve3>| -> &'static str {
                    match c { Some(rcad_kernel::geom::Curve3::Line(_)) => "Line",
                        Some(rcad_kernel::geom::Curve3::Circle(_)) => "Circle",
                        Some(rcad_kernel::geom::Curve3::Ellipse(_)) => "Ellipse",
                        Some(rcad_kernel::geom::Curve3::BSpline(_)) => "BSpline",
                        Some(rcad_kernel::geom::Curve3::Bezier(_)) => "Bezier",
                        Some(rcad_kernel::geom::Curve3::Hyperbola(_)) => "Hyperbola",
                        Some(rcad_kernel::geom::Curve3::Parabola(_)) => "Parabola",
                        _ => "Other" }
                };
                let ep = |ei: usize| -> (DVec3, DVec3) {
                    let c = self.ds.edge_curve(ei);
                    let r = self.ds.shapes[ei].shape.as_edge().map(|e| e.range).unwrap_or([0.0,0.0]);
                    match c { Some(c) => (c.point_at(r[0]), c.point_at(r[1])), None => (DVec3::ZERO, DVec3::ZERO) }
                };
                let (a0, a1) = ep(n_e1);
                let (b0, b1) = ep(n_e2);
                eprintln!("[EE-DBG] pair ({},{}) flag={}/{} t1={} t2={} e1:[({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})] e2:[({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})]",
                    n_e1, n_e2,
                    self.ds.shapes[n_e1].has_flag(), self.ds.shapes[n_e2].has_flag(),
                    tn(self.ds.edge_curve(n_e1)), tn(self.ds.edge_curve(n_e2)),
                    a0.x, a0.y, a0.z, a1.x, a1.y, a1.z,
                    b0.x, b0.y, b0.z, b1.x, b1.y, b1.z);
            }
            // L189-196: skip degenerated edges
            if self.ds.shapes[n_e1].has_flag() || self.ds.shapes[n_e2].has_flag() {
                continue;
            }

            // L200-210: get PB lists for both edges (clone to avoid borrow conflict)
            let a_lpb1: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e1).to_vec();
            let a_lpb2: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e2).to_vec();
            if a_lpb1.is_empty() || a_lpb2.is_empty() {
                if ee_debug {
                    eprintln!("[EE-DBG] pair ({},{}) empty-pbs lpb1={} lpb2={}", n_e1, n_e2, a_lpb1.len(), a_lpb2.len());
                }
                continue;
            }

            // L212-265: iterate PB1 脳 PB2
            let mut pb_box_cache: std::collections::HashMap<u64, rcad_kernel::math::bnd::BndBox> =
                std::collections::HashMap::new();

            for p1 in &a_lpb1 {
                // GetPBBox for PB1
                let (mut t11, mut t12, mut ts11, mut ts12) = (0.0, 0.0, 0.0, 0.0);
                let mut bb1 = rcad_kernel::math::bnd::BndBox::new();
                if !self.get_pb_box(n_e1, p1, &mut pb_box_cache,
                    &mut t11, &mut t12, &mut ts11, &mut ts12, &mut bb1) {
                    if ee_debug {
                        eprintln!("[EE-DBG] pair ({},{}) pb1-getpb false", n_e1, n_e2);
                    }
                    continue;
                }

                for p2 in &a_lpb2 {
                    let (mut t21, mut t22, mut ts21, mut ts22) = (0.0, 0.0, 0.0, 0.0);
                    let mut bb2 = rcad_kernel::math::bnd::BndBox::new();
                    if !self.get_pb_box(n_e2, p2, &mut pb_box_cache,
                        &mut t21, &mut t22, &mut ts21, &mut ts22, &mut bb2) {
                        if ee_debug {
                            eprintln!("[EE-DBG] pair ({},{}) pb2-getpb false", n_e1, n_e2);
                        }
                        continue;
                    }

                    // L245-248: box overlap check
                    if bb1.is_out_box(&bb2) {
                        if ee_debug {
                            eprintln!("[EE-DBG] pair ({},{}) box-out t1=[{:.3},{:.3}] ts=[{:.3},{:.3}] bb1={:?} bb2={:?}",
                                n_e1, n_e2, t11, t12, ts11, ts12, bb1.raw_min(), bb2.raw_min());
                        }
                        continue;
                    }

                    // L252: bExpressCompute = PB1 and PB2 have same bounding vertices
                    let (n_v11, n_v12) = { let r = p1.0.read().unwrap(); r.indices() };
                    let (n_v21, n_v22) = { let r = p2.0.read().unwrap(); r.indices() };
                    let b_express = (n_v11 == n_v21 && n_v12 == n_v22)
                                 || (n_v12 == n_v21 && n_v11 == n_v22);

                    if ee_debug {
                        eprintln!("[EE-DBG] pair ({},{}) creating ee b_express={} t1=[{:.6},{:.6}] t2=[{:.6},{:.6}]", n_e1, n_e2, b_express, t11, t12, t21, t22);
                    }

                    // OCCT L254-264: create EdgeEdge, intersect
                    let mut ee = int_tools::edge_edge::EdgeEdgeIntersector::new();
                    ee.use_quick_coincidence_check(b_express);
                    ee.set_edges(n_e1, [t11, t12], n_e2, [t21, t22], &self.ds);
                    ee.set_fuzzy_value(self.my_fuzzy_value);
                    if ee_debug {
                        eprintln!("[EE-DBG] pair ({},{}) before perform", n_e1, n_e2);
                    }
                    ee.perform();
                    if ee_debug {
                        eprintln!("[EE-DBG] pair ({},{}) after perform done={} ncp={}", n_e1, n_e2, ee.is_done(), ee.common_parts().len());
                    }

                    if !ee.is_done() {
                        self.my_report.add_warning(
                            Alert::IntersectionFailed(n_e1, n_e2));
                        continue;
                    }

                    let a_cparts = ee.common_parts();
                    let a_nb_cprts = a_cparts.len();
                    if a_nb_cprts == 0 {
                        continue;
                    }

                    // OCCT L355-553: process each common part
                    for (i_cp, cp) in a_cparts.iter().enumerate() {
                        // OCCT aCP.Type() 鈥?set by IntTools_EdgeEdge::AddSolution /
                        // ComputeLineLine. (The old range-span heuristic misclassified
                        // line-line VERTEX parts as EDGE because their Range1 is the
                        // tolerance band [aT-dt, aT+dt], not a collapsed range.)
                        let a_type = if cp.is_edge {
                            ShapeType::Edge
                        } else {
                            ShapeType::Vertex
                        };

                        match a_type {
                            ShapeType::Vertex => {
                                // OCCT L370-373: skip if PB not splittable
                                let b_is_pb_splittable1 = {
                                    let r = p1.0.read().unwrap();
                                    r.is_splittable()
                                };
                                let b_is_pb_splittable2 = {
                                    let r = p2.0.read().unwrap();
                                    r.is_splittable()
                                };
                                if !b_is_pb_splittable1 || !b_is_pb_splittable2 {
                                    if ee_debug {
                                        eprintln!("[EE-DBG] pair ({},{}) FILTER not-splittable {}/{}", n_e1, n_e2, b_is_pb_splittable1, b_is_pb_splittable2);
                                    }
                                    continue;
                                }

                                let a_t1 = cp.vertex_param1;
                                let a_t2 = cp.vertex_param2;

                                // OCCT L381-394: IsOnPave checks in 4 shrunk regions
                                let a_tol = rcad_kernel::CONFUSION; // Precision::Confusion()
                                let a_cr1 = (cp.range1[0], cp.range1[1]);
                                let a_cr2 = (cp.ranges2[0][0], cp.ranges2[0][1]);
                                let a_r11_first = t11.min(ts11);
                                let a_r11_last = t11.max(ts11);
                                let a_r12_first = ts12.min(t12);
                                let a_r12_last = ts12.max(t12);
                                let a_r21_first = t21.min(ts21);
                                let a_r21_last = t21.max(ts21);
                                let a_r22_first = ts22.min(t22);
                                let a_r22_last = ts22.max(t22);

                                // OCCT: IsOnPave checks for 4 region boundaries
                                let mut b_is_on_pave = [
                                    algo_tools::is_on_pave_1(a_t1, a_r11_first, a_r11_last, a_tol)
                                        || algo_tools::is_on_pave_1(a_r11_first, a_cr1.0, a_cr1.1, a_tol),
                                    algo_tools::is_on_pave_1(a_t1, a_r12_first, a_r12_last, a_tol)
                                        || algo_tools::is_on_pave_1(a_r12_last, a_cr1.0, a_cr1.1, a_tol),
                                    algo_tools::is_on_pave_1(a_t2, a_r21_first, a_r21_last, a_tol)
                                        || algo_tools::is_on_pave_1(a_r21_first, a_cr2.0, a_cr2.1, a_tol),
                                    algo_tools::is_on_pave_1(a_t2, a_r22_first, a_r22_last, a_tol)
                                        || algo_tools::is_on_pave_1(a_r22_last, a_cr2.0, a_cr2.1, a_tol),
                                ];

                                // OCCT L396-403: if intersection is on existing paves on both edges, skip
                                if (b_is_on_pave[0] && b_is_on_pave[2])
                                    || (b_is_on_pave[0] && b_is_on_pave[3])
                                    || (b_is_on_pave[1] && b_is_on_pave[2])
                                    || (b_is_on_pave[1] && b_is_on_pave[3])
                                {
                                    if ee_debug {
                                        eprintln!("[EE-DBG] pair ({},{}) FILTER on-pave-both {:?} t1={:.6} t2={:.6} cr1=[{:.6},{:.6}] cr2=[{:.6},{:.6}] r11=[{:.6},{:.6}] r12=[{:.6},{:.6}] r21=[{:.6},{:.6}] r22=[{:.6},{:.6}] t11={:.6} ts11={:.6} t12={:.6} ts12={:.6} t21={:.6} ts21={:.6} t22={:.6} ts22={:.6}",
                                            n_e1, n_e2, b_is_on_pave, a_t1, a_t2, a_cr1.0, a_cr1.1, a_cr2.0, a_cr2.1,
                                            a_r11_first, a_r11_last, a_r12_first, a_r12_last,
                                            a_r21_first, a_r21_last, a_r22_first, a_r22_last,
                                            t11, ts11, t12, ts12, t21, ts21, t22, ts22);
                                    }
                                    continue;
                                }

                                // OCCT L405-417: ForceInterfVE for vertices on pave boundaries
                                let n_v_arr = [n_v11, n_v12, n_v21, n_v22];
                                // OCCT: aPB = (j < 2) ? aPB2 : aPB1 鈥?cross assignment
                                let p_b_arr = [&p2, &p2, &p1, &p1];
                                let mut is_v_exists = false;
                                for j in 0..4 {
                                    if b_is_on_pave[j] {
                                        b_is_on_pave[j] = self.force_interf_ve(
                                            n_v_arr[j], p_b_arr[j], &mut a_m_edges);
                                        if b_is_on_pave[j] {
                                            is_v_exists = true;
                                        }
                                    }
                                }

                                // OCCT L419-420: MakeNewVertex(aE1, aT1, aE2, aT2, aVnew); aPnew = Pnt(aVnew)
                                // OCCT (BOPTools_AlgoTools_2.cxx L224-250): point = midpoint of the two
                                // curve points at aT1/aT2, tolerance = max(edgeTol1, edgeTol2) + 0.5*dist.
                                let p_e1 = self.ds.edge_curve(n_e1).map(|c| c.point_at(a_t1)).unwrap_or(cp.bounding_point1);
                                let p_e2 = self.ds.edge_curve(n_e2).map(|c| c.point_at(a_t2)).unwrap_or(cp.bounding_point1);
                                let (vnew_pt, vnew_tol) = crate::bop::tools::algo_tools::make_new_vertex(
                                    p_e1, self.ds.edge_tolerance(n_e1), p_e2, self.ds.edge_tolerance(n_e2));

                                // OCCT L422-451: isVExists check
                                if is_v_exists {
                                    // OCCT L430-431: BRepAdaptor_Curve(aE1).Value(aT1)
                                    let a_p_on_e1 = self.ds.edge_curve(n_e1)
                                        .map(|c| c.point_at(a_t1)).unwrap_or(vnew_pt);
                                    let a_p_on_e2 = self.ds.edge_curve(n_e2)
                                        .map(|c| c.point_at(a_t2)).unwrap_or(vnew_pt);
                                    // OCCT L432: if (aPOnE1.Distance(aPOnE2) > Precision::Intersection()) continue;
                                    if (a_p_on_e1 - a_p_on_e2).length() > rcad_kernel::precision::INTERSECTION {
                                        if ee_debug {
                                            eprintln!("[EE-DBG] pair ({},{}) FILTER isVExists-far {:.9}", n_e1, n_e2, (a_p_on_e1 - a_p_on_e2).length());
                                        }
                                        continue;
                                    }
                                    // OCCT L440-451: update each vertex where bIsOnPave[j] is true
                                    for j in 0..4 {
                                        if b_is_on_pave[j] {
                                            let v_pt = self.ds.vertex_point_by_idx(n_v_arr[j]);
                                            let a_dist_pp = (vnew_pt - v_pt).length();
                                            self.update_vertex(n_v_arr[j], a_dist_pp);
                                            self.my_verts_to_avoid_extension.insert(n_v_arr[j]);
                                        }
                                    }
                                }

                                // OCCT L454-466: analytical tolerance boost for Line/Circle
                                let mut a_tol_vnew = vnew_tol;
                                {
                                    let c1 = self.ds.edge_curve(n_e1);
                                    let c2 = self.ds.edge_curve(n_e2);
                                    let b_analytical = match (&c1, &c2) {
                                        (Some(rcad_kernel::geom::Curve3::Line(_)), Some(rcad_kernel::geom::Curve3::Circle(_)))
                                        | (Some(rcad_kernel::geom::Curve3::Circle(_)), Some(rcad_kernel::geom::Curve3::Line(_))) => true,
                                        _ => false,
                                    };
                                    if b_analytical {
                                        let range_len1 = a_cr1.1 - a_cr1.0;
                                        let range_len2 = a_cr2.1 - a_cr2.0;
                                        let a_tol_min = if matches!(c1.as_ref().unwrap(), rcad_kernel::geom::Curve3::Line(_)) {
                                            range_len1 / 2.0
                                        } else {
                                            range_len2 / 2.0
                                        };
                                        if a_tol_min > a_tol_vnew {
                                            a_tol_vnew = a_tol_min;
                                        }
                                    }
                                }

                                // OCCT L468-510: bounding vertex closeness check
                                let mut skip_new_vertex = false;
                                {
                                    let a_mv: std::collections::HashSet<usize> =
                                        [n_v11, n_v12].iter().copied().collect();
                                    for &n_v_candidate in &[n_v21, n_v22] {
                                        if a_mv.contains(&n_v_candidate) {
                                            let vx_tol = self.ds.vertex_tolerance_by_idx(n_v_candidate);
                                            let vx_pt = self.ds.vertex_point_by_idx(n_v_candidate);
                                            let d2 = (vnew_pt - vx_pt).length_squared();
                                            let dt2 = 100.0 * (a_tol_vnew + vx_tol) * (a_tol_vnew + vx_tol);
                                            if d2 < dt2 {
                                                skip_new_vertex = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if skip_new_vertex {
                                    if ee_debug {
                                        eprintln!("[EE-DBG] pair ({},{}) FILTER close-to-shared-vertex", n_e1, n_e2);
                                    }
                                    continue;
                                }

                                // OCCT L513-518: add InterfEE
                                // OCCT: aEEs.Appended() pushes directly to DS; iX = aEEs.Length()-1
                                // rcad: new_ee collects locally; absolute index = ee_base + new_ee.len()
                                if ee_debug {
                                    eprintln!("[EE-DBG] pair ({},{}) VERTEX-EE add t1={:.6} t2={:.6}", n_e1, n_e2, a_t1, a_t2);
                                }
                                let idx_interf = ee_base + new_ee.len();
                                mvcpb.push(CoupleOfPBs {
                                    pb1: p1.clone(),
                                    pb2: p2.clone(),
                                    index_interf: idx_interf,
                                    tolerance: a_tol_vnew,
                                    index: usize::MAX,
                                    point: vnew_pt,
                                });

                                new_ee.push(InterferenceEE {
                                    e1: n_e1, e2: n_e2,
                                    point: vnew_pt,
                                    param1: a_t1, param2: a_t2,
                                    new_vertex: usize::MAX,
                                    range1: [a_t1, a_t1],
                                    range2: [a_t2, a_t2],
                                });
                                self.ds.add_interf(n_e1, n_e2);
                            }
                            ShapeType::Edge => {
                                // OCCT L529-533: only process EDGE with single common part
                                // OCCT: if (aNbCPrts > 1) { break; } 鈥?break switch, NOT for loop
                                // rcad: must not continue/break for loop; only skip the case body
                                let b_process_edge = a_nb_cprts <= 1;
                                if b_process_edge {
                                    // OCCT L535-539: HasSameBounds check
                                    let b_has_same_bounds = (n_v11 == n_v21 && n_v12 == n_v22)
                                                         || (n_v12 == n_v21 && n_v11 == n_v22);
                                    if ee_debug {
                                        let vp = |v: usize| {
                                            if v < self.ds.nb_shapes() {
                                                let p = self.ds.vertex_point_by_idx(v);
                                                format!("({:.3},{:.3},{:.3})", p.x, p.y, p.z)
                                            } else { "?".to_string() }
                                        };
                                        eprintln!("[EE-DBG] pair ({},{}) EDGE-case ncp={} samebounds={} nv=({},{},{},{}) pts=[{},{},{},{}]",
                                            n_e1, n_e2, a_nb_cprts, b_has_same_bounds, n_v11, n_v12, n_v21, n_v22,
                                            vp(n_v11), vp(n_v12), vp(n_v21), vp(n_v22));
                                    }
                                    if b_has_same_bounds {
                                        // OCCT L542-547: add InterfEE with common part
                                        new_ee.push(InterferenceEE {
                                            e1: n_e1, e2: n_e2,
                                            point: cp.bounding_point1,
                                            param1: cp.range1[0], param2: cp.ranges2[0][0],
                                            new_vertex: usize::MAX,
                                            range1: cp.range1,
                                            range2: cp.ranges2[0],
                                        });
                                        self.ds.add_interf(n_e1, n_e2);
                                        cb_pairs.push((p1.clone(), p2.clone()));
                                    }
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
        }

        // Push all collected interferences
        self.ds.interf_ee.extend(new_ee);

        // OCCT L561: PerformCommonBlocks
        // OCCT BOPAlgo_Tools::FillMap(aPB1, aPB2, aMPBLPB) fills a bidirectional
        // adjacency map (key -> list of connected PBs), then MakeBlocks groups
        // connected PBs into one Common Block per connected component.
        // OCCT aMPBLPB (BOPAlgo_PaveFiller_3.cxx L524-534) is NCollection_IndexedDataMap
        // 鈥?the MakeBlocks iteration (BOPAlgo_Tools.hxx L49-53, FindKey(i)) and
        // the chain start order are the insertion order of cb_pairs; a HashMap
        // would randomize the CommonBlock pave-block order, which later feeds
        // cb_pcurve_source's first-with-pcurve selection and hence the copied
        // pcurve parameterization.
        if !cb_pairs.is_empty() {
            let mut a_mpblpb: indexmap::IndexMap<u64, (SharedPB, Vec<u64>)> =
                indexmap::IndexMap::new();
            for (pb1, pb2) in &cb_pairs {
                let k1 = std::sync::Arc::as_ptr(&pb1.0) as u64;
                let k2 = std::sync::Arc::as_ptr(&pb2.0) as u64;
                a_mpblpb.entry(k1).or_insert_with(|| (pb1.clone(), Vec::new())).1.push(k2);
                a_mpblpb.entry(k2).or_insert_with(|| (pb2.clone(), Vec::new())).1.push(k1);
            }
            // OCCT MakeBlocks (BOPAlgo_Tools.hxx L46-80): connected components
            // grown in FIFO order (aChain list iterated from the start, appended
            // elements visited later) 鈥?a stack would visit in a different order.
            let mut a_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
            for (&n, _) in &a_mpblpb {
                if a_fence.contains(&n) {
                    continue;
                }
                let mut a_block: Vec<SharedPB> = Vec::new();
                let mut a_queue: std::collections::VecDeque<u64> = std::collections::VecDeque::new();
                a_queue.push_back(n);
                while let Some(n1) = a_queue.pop_front() {
                    if !a_fence.insert(n1) {
                        continue;
                    }
                    if let Some((pb, a_linked)) = a_mpblpb.get(&n1) {
                        a_block.push(pb.clone());
                        for &n2 in a_linked {
                            if !a_fence.contains(&n2) {
                                a_queue.push_back(n2);
                            }
                        }
                    }
                }
                // OCCT L134-137: if (aNbPB < 2) continue;
                if a_block.len() >= 2 {
                    self.ds.add_common_block(&a_block);
                }
            }
        }
        // OCCT L563: UpdateVerticesOfCB
        self.update_vertices_of_cb();

        // OCCT L565-569: PerformNewVertices(aMVCPB, ...); if (HasErrors()) return;
        if !mvcpb.is_empty() {
            self.perform_new_vertices(&mvcpb, true);
            if self.has_errors() {
                return;
            }
            // OCCT L571-583: remove mvcpb edges from aMEdges
            for cpb in &mvcpb {
                let (n_e1, n_e2) = (cpb.pb1.0.read().unwrap().original_edge,
                                    cpb.pb2.0.read().unwrap().original_edge);
                a_m_edges.remove(n_e1);
                a_m_edges.remove(n_e2);
            }
        }
        // OCCT L584: SplitPaveBlocks(aMEdges, false)
        if !a_m_edges.is_empty() {
            self.split_pave_blocks(&a_m_edges, false);
        }

        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            for (i, si) in self.ds.shapes.iter().enumerate() {
                if si.shape_type == ShapeType::Edge {
                    let pbs = self.ds.edge_pave_blocks(i);
                    if !pbs.is_empty() {
                        eprintln!("[EE-DBG] edge {} pbs={} :", i, pbs.len());
                        for pb in pbs {
                            let r = pb.0.read().unwrap();
                            eprintln!("  [EE-DBG]   pb range=[{:.6},{:.6}] v1={} v2={}", r.range().0, r.range().1, r.pave1.vertex_idx, r.pave2.vertex_idx);
                        }
                    }
                }
            }
        }
    }

    // ====================================================================
    // VF 鈥?OCCT BOPAlgo_PaveFiller_5.cxx L409-471
    // ====================================================================
    // OCCT BOPAlgo_PaveFiller::PerformVF (PaveFiller_4.cxx L139-399).
    // rcad: simplified 鈥?skips FaceInfo/TreatVerticesEE/complex projection.
    fn perform_vf(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L141-142: myIterator->Initialize(VERTEX, FACE); iSize = ExpectedLength()
        let my_iterator = match &mut self.my_iterator {
            Some(it) => it,
            None => return,
        };
        my_iterator.initialize(ShapeType::Vertex, ShapeType::Face);
        let i_size = my_iterator.pairs(ShapeType::Vertex, ShapeType::Face).len();

        // OCCT L147-160: GlueFull mode 鈥?init FaceInfo and return
        if self.my_glue == GlueEnum::GlueFull {
            let pairs: Vec<(usize, usize)> = my_iterator.pairs(ShapeType::Vertex, ShapeType::Face).to_vec();
            for &(n_v, n_f) in &pairs {
                let (n_v, n_f) = if self.ds.shapes[n_v].shape_type == ShapeType::Vertex { (n_v, n_f) } else { (n_f, n_v) };
                if !self.ds.is_sub_shape(n_v, n_f) {
                    self.ds.change_face_info(n_f);
                }
            }
            return;
        }

        // OCCT L163-169: if (!iSize) { aVFs.SetIncrement(10); TreatVerticesEE(); return; }
        if i_size == 0 {
            self.treat_vertices_ee();
            return;
        }

        // OCCT L176-180: aMVFPairs 鈥?avoid repeated intersection of the same
        // (SD vertex, face) pair; the value is NCollection_Map<int> (bucket
        // order) of the original vertices in the group. OCCT iterates the
        // collected solvers (aVVF vector, insertion order) at L249-298, so the
        // pair map is an IndexMap keyed by (nVx, nF) in insertion order.
        let mut a_mvf_pairs: indexmap::IndexMap<(usize, usize), crate::bop::algo::occt_map::OcctMapInt> =
            indexmap::IndexMap::new();
        let pairs: Vec<(usize, usize)> = my_iterator.pairs(ShapeType::Vertex, ShapeType::Face).to_vec();
        for &(n0, n1) in &pairs {
            let (n_v, n_f) = if self.ds.shapes[n0].shape_type == ShapeType::Vertex { (n0, n1) } else { (n1, n0) };
            // OCCT L189-191: if (myDS->IsSubShape(nV, nF)) continue;
            if self.ds.is_sub_shape(n_v, n_f) { continue; }
            // OCCT L194-197: if (myDS->HasInterf(nV, nF)) continue;
            if self.ds.has_interf(n_v, n_f) { continue; }
            // OCCT L199: myDS->ChangeFaceInfo(nF);
            self.ds.change_face_info(n_f);
            // OCCT L200-203: if (myDS->HasInterfShapeSubShapes(nV, nF)) continue;
            if self.ds.has_interf_shape_sub_shapes(n_v, n_f, true) { continue; }
            // OCCT L205-209: nVx = (HasShapeSD(nV, nVSD)) ? nVSD : nV
            let mut n_vx = n_v;
            self.ds.has_shape_sd(n_v, &mut n_vx);
            // OCCT L211-220: dedup by (nVx, nF) 鈥?aMVFPairs.ChangeSeek/Bound.
            let a_vf_pair = (n_vx, n_f);
            let entry = a_mvf_pairs.entry(a_vf_pair).or_insert_with(crate::bop::algo::occt_map::OcctMapInt::new);
            entry.add(n_v);
        }

        // OCCT L249-298: process each unique (nVx, nF) pair 鈥?the aVVF vector
        // is iterated in insertion order (IndexMap order here); the original
        // vertices aMV are NCollection_Map<int> iterated in bucket order.
        for ((n_vx, n_f), verts) in &a_mvf_pairs {
            if the_range.user_break() { return; }
            // OCCT L257-266: BOPAlgo_VertexFace::Perform 鈫?myContext->ComputeVF
            let (i_flag, a_u, a_v, a_tol_v_new) =
                self.my_context.compute_vf(*n_vx, *n_f, &self.ds, self.my_fuzzy_value);
            if i_flag != 0 { continue; }

            let mut n_vx_cur = *n_vx;
            for n_v in verts.iter_keys() {
                // 1
                let mut a_vf = InterferenceVF { vertex: n_v, face: *n_f, u: a_u, v: a_v, index_new: None };
                // 2
                self.ds.add_interf(n_v, *n_f);
                // 3 nVx = UpdateVertex(nV, aTolVNew)
                n_vx_cur = self.update_vertex(n_v, a_tol_v_new);
                // 4 if (IsNewShape(nVx)) aVF.SetIndexNew(nVx)
                if self.ds.is_new_shape(n_vx_cur) {
                    a_vf.index_new = Some(n_vx_cur);
                }
                self.ds.interf_vf.push(a_vf);
            }
            // 5 update FaceInfo vertices_in
            self.ds.change_face_info(*n_f).vertices_in.insert(n_vx_cur);
        }

        // OCCT L300: TreatVerticesEE()
        let n_vf_before = self.ds.interf_vf.len();
        self.treat_vertices_ee();
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            eprintln!("[EE-DBG] perform_vf: main added {} VF, treat_vertices_ee added {}", n_vf_before, self.ds.interf_vf.len() - n_vf_before);
        }
    }

    /// OCCT BOPAlgo_PaveFiller::TreatVerticesEE (PaveFiller_4.cxx L305-390).
    /// For each new vertex created by EE intersections, check if it lies on a
    /// source face (verticesOn) and, if so, add the VF interference.
    fn treat_vertices_ee(&mut self) {
        // OCCT L321-331: collect EE interferences' new vertices (IndexNew)
        let mut a_liv: Vec<usize> = Vec::new();
        let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for a_ee in &self.ds.interf_ee {
            if a_ee.new_vertex != usize::MAX {
                if a_mi.insert(a_ee.new_vertex) {
                    a_liv.push(a_ee.new_vertex);
                }
            }
        }
        if a_liv.is_empty() { return; }

        // OCCT L339-352: collect all source faces
        let mut a_lif: Vec<usize> = Vec::new();
        let a_nb_s = self.ds.nb_source_shapes;
        for n_f in 0..a_nb_s {
            if self.ds.shapes[n_f].shape_type == ShapeType::Face {
                a_lif.push(n_f);
            }
        }
        if a_lif.is_empty() { return; }

        // OCCT L356-389: for each (face, EE-vertex) pair selected by the
        // BOPDS_SubIterator (bounding-box overlap), ChangeFaceInfo(nF) is
        // called UNCONDITIONALLY (the FaceInfo must exist even when the vertex
        // is not on the face 鈥?OCCT PaveFiller_4.cxx L357-360).  Only when the
        // vertex is not already ON the face, ComputeVF is run and the
        // interference recorded.
        for &n_f in &a_lif {
            for &n_v in &a_liv {
                // OCCT BOPDS_SubIterator: pairs with overlapping bounding boxes.
                if !crate::bop::ds::iterator::boxes_overlap(
                    &self.ds.shapes[n_f], &self.ds.shapes[n_v], self.my_fuzzy_value)
                {
                    continue;
                }
                // OCCT: BOPDS_FaceInfo& aFI = myDS->ChangeFaceInfo(nF);
                let contains = {
                    let a_fi = self.ds.change_face_info(n_f);
                    a_fi.vertices_on.contains(&n_v)
                };
                if contains { continue; }
                let (i_flag, a_t1, a_t2, _dummy) =
                    self.my_context.compute_vf(n_v, n_f, &self.ds, self.my_fuzzy_value);
                if i_flag == 0 {
                    // 1
                    self.ds.interf_vf.push(InterferenceVF {
                        vertex: n_v, face: n_f, u: a_t1, v: a_t2, index_new: None,
                    });
                    // 2
                    self.ds.add_interf(n_v, n_f);
                    // 3
                    self.ds.change_face_info(n_f).vertices_in.insert(n_v);
                }
            }
        }
    }

    // ====================================================================
    // EF 鈥?OCCT BOPAlgo_PaveFiller_5.cxx L165-580
    // ====================================================================
    /// OCCT BOPAlgo_PaveFiller::PerformEF (PaveFiller_5.cxx L165-580).
    fn perform_ef(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L167: FillShrunkData(TopAbs_EDGE, TopAbs_FACE)
        self.fill_shrunk_data(ShapeType::Edge, ShapeType::Face);

        // OCCT L169-175: myIterator->Initialize(EDGE, FACE); check iSize
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Edge, ShapeType::Face).to_vec()
        } else {
            return;
        };
        if pairs.is_empty() {
            return;
        }

        // OCCT L179-191: GlueFull mode 鈥?init FaceInfo and return
        if self.my_glue == GlueEnum::GlueFull {
            for &(n_e, n_f) in &pairs {
                if !self.ds.shapes[n_e].has_flag() {
                    self.ds.change_face_info(n_f);
                }
            }
            return;
        }

        // OCCT L194-217: locals
        let mut a_mi_efc: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();
        let mut a_mvcpb: Vec<CoupleOfPBs> = Vec::new();
        let mut a_mpbli: indexmap::IndexMap<u64, (SharedPB, Vec<usize>)> =
            indexmap::IndexMap::new();
        let mut a_dmpb_box: std::collections::HashMap<u64, rcad_kernel::math::bnd::BndBox> =
            std::collections::HashMap::new();
        // OCCT: aVEdgeFace 鈥?collected (edge, face, pb, range, shrunk, express)
        struct EFCandidate {
            n_e: usize, n_f: usize, pb: SharedPB,
            // OCCT: aEdgeFace.SetRange(aPBRange) 鈥?CorrectRange'd PB range used
            // as the BeanFaceIntersector parameters.
            range: (f64, f64),
            // OCCT: the raw PB range (aPB->Range) re-read in the processing loop
            // for aR1=(aT1,aTS1) / aR2=(aTS2,aT2).
            raw_range: (f64, f64),
            shrunk: (f64, f64),
            express: bool, tol_e: f64, tol_f: f64,
        }
        let mut a_v_edge_face: Vec<EFCandidate> = Vec::new();

        // OCCT L219-307: loop 1 鈥?collect candidates
        for &(n0, n1) in &pairs {
            let (n_e, n_f) = if self.ds.shapes[n0].shape_type == ShapeType::Edge { (n0, n1) } else { (n1, n0) };
            // OCCT L227-231: skip degenerated edges
            if self.ds.shapes[n_e].has_flag() {
                continue;
            }
            // OCCT L237: BOPDS_FaceInfo& aFI = myDS->ChangeFaceInfo(nF) 鈥?the
            // FaceInfo is created for every surviving (E,F) pair, collecting
            // the current boundary PBs into PaveBlocksOn.
            let a_fi = self.ds.change_face_info(n_f);
            let a_mpbf = a_fi.pave_blocks_on.clone();
            let a_mv_in = a_fi.vertices_in.clone();
            let a_mv_on = a_fi.vertices_on.clone();
            drop(a_fi);
            let a_tol_e = self.ds.edge_tolerance(n_e);
            let a_tol_f = self.ds.face_tolerance(n_f);
            // OCCT L246: PB list
            let a_lpb: Vec<SharedPB> = self.ds.edge_pave_blocks(n_e).to_vec();
            if a_lpb.is_empty() {
                continue;
            }
            for a_pb in &a_lpb {
                // OCCT L256-260: RealPaveBlock + skip if in face's PaveBlocksOn
                let a_pbr = self.ds.real_pave_block(a_pb);
                let pbr_ptr = std::sync::Arc::as_ptr(&a_pbr.0) as u64;
                // OCCT: aMPBF.Contains(aPBR) 鈥?the set holds PB handles, so a
                // direct pointer-id membership test matches exactly.
                let is_on_face = a_mpbf.contains(&pbr_ptr);
                if std::env::var("RCAD_EE_DEBUG").is_ok() {
                    if n_e == 55 || n_e == 62 || n_f == 22 {
                        let mut on_edges: Vec<String> = Vec::new();
                        for &pp in &a_mpbf {
                            // find the edge that owns this PB pointer
                            let mut found: Option<usize> = None;
                            for (ei, pbs) in &self.ds.pave_blocks_pool {
                                for pb in pbs {
                                    let rp = self.ds.real_pave_block(pb);
                                    if std::sync::Arc::as_ptr(&rp.0) as u64 == pp {
                                        let r = rp.0.read().unwrap();
                                        found = Some(*ei);
                                        on_edges.push(format!("e{}:[{:.2},{:.2}]", ei, r.pave1().parameter(), r.pave2().parameter()));
                                        break;
                                    }
                                }
                                if found.is_some() { break; }
                            }
                            if found.is_none() { on_edges.push("?".to_string()); }
                        }
                        eprintln!("[EF-DBG]   cand e={} f={} onFace={} pbon_edges={:?}", n_e, n_f, is_on_face, on_edges);
                    }
                }
                if is_on_face {
                    continue;
                }
                // OCCT L262-266: GetPBBox
                let (mut a_t1, mut a_t2, mut a_ts1, mut a_ts2) = (0.0, 0.0, 0.0, 0.0);
                let mut a_bb_e = rcad_kernel::math::bnd::BndBox::new();
                if !self.get_pb_box(n_e, a_pb, &mut a_dmpb_box,
                    &mut a_t1, &mut a_t2, &mut a_ts1, &mut a_ts2, &mut a_bb_e) {
                    continue;
                }
                // OCCT L268-271: box overlap
                let a_bb_f = self.ds.shapes[n_f].bbox.clone();
                if a_bb_f.is_out_box(&a_bb_e) {
                    continue;
                }
                // OCCT L273-276: bExpressCompute
                let (n_v1, n_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
                let b_v1 = a_mv_in.contains(&n_v1) || a_mv_on.contains(&n_v1);
                let b_v2 = a_mv_in.contains(&n_v2) || a_mv_on.contains(&n_v2);
                let b_express = b_v1 && b_v2;
                // OCCT L289-297: CorrectRange for shrunk and PB range
                let curve = self.ds.edge_curve(n_e);
                let a_corrected_ts = match &curve {
                    Some(c) => crate::bop::tools::algo_tools::correct_range_ef(
                        c, a_ts1, a_ts2, a_tol_e, a_tol_f),
                    None => (a_ts1, a_ts2),
                };
                let a_corrected_range = match &curve {
                    Some(c) => crate::bop::tools::algo_tools::correct_range_ef(
                        c, a_t1, a_t2, a_tol_e, a_tol_f),
                    None => (a_t1, a_t2),
                };
                // OCCT L299-305: myFPBDone
                self.my_fpb_done.entry(n_f).or_default().insert(
                    std::sync::Arc::as_ptr(&a_pb.0) as u64);
                a_v_edge_face.push(EFCandidate {
                    n_e, n_f, pb: a_pb.clone(),
                    range: a_corrected_range, raw_range: (a_t1, a_t2),
                    shrunk: a_corrected_ts,
                    express: b_express, tol_e: a_tol_e, tol_f: a_tol_f,
                });
            }
        }

        // OCCT L324-571: loop 2 鈥?compute and process common parts
        for cand in &a_v_edge_face {
            let a_nb_cprts;
            let (i_flag, common_parts, min_dist) = self.my_context.compute_ef(
                cand.n_e, cand.n_f, cand.range.0, cand.range.1, cand.express, &self.ds, self.my_fuzzy_value);
            a_nb_cprts = common_parts.len();
            if std::env::var("RCAD_EE_DEBUG").is_ok() {
                let epts = self.ds.edge_curve(cand.n_e).map(|c| {
                    let p0 = c.point_at(cand.raw_range.0);
                    let p1 = c.point_at(cand.raw_range.1);
                    format!("({:.2},{:.2},{:.2})-({:.2},{:.2},{:.2})", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z)
                }).unwrap_or_default();
                eprintln!("[EF-DBG] cand e={}(r{}) f={}(r{}) range=[{:.4},{:.4}] ncp={} min={:.3e} epts={}", cand.n_e, self.ds.rank(cand.n_e), cand.n_f, self.ds.rank(cand.n_f), cand.range.0, cand.range.1, a_nb_cprts, min_dist, epts);
            }
            if i_flag != 0 {
                continue;
            }
            if a_nb_cprts == 0 {
                // OCCT L348-361: record minimal distance if no common part
                if min_dist < f64::MAX && min_dist > cand.tol_e + cand.tol_f {
                    self.my_distances.entry((cand.n_e, cand.n_f)).or_default().push(
                        EdgeRangeDistance::new(cand.range.0, cand.range.1, min_dist));
                }
                continue;
            }
            let a_pb = &cand.pb;
            // OCCT L380: aPB->Range(aT1, aT2) 鈥?the raw PB range (not the
            // CorrectRange'd one used for the intersection).
            let (a_t1, a_t2) = cand.raw_range;
            let (n_v1, n_v2) = { let r = a_pb.0.read().unwrap(); r.indices() };
            let b_is_pb_splittable = { let r = a_pb.0.read().unwrap(); r.is_splittable() };
            let (mut a_ts1, mut a_ts2) = cand.shrunk;
            // OCCT L373-380: VERTEX type 鈫?ReduceIntersectionRange
            if !common_parts.is_empty() && !common_parts[0].2 {
                self.reduce_intersection_range(n_v1, n_v2, cand.n_e, cand.n_f, &mut a_ts1, &mut a_ts2);
            }
            let a_r1 = (a_t1, a_ts1);
            let a_r2 = (a_ts2, a_t2);
            let a_fi = self.ds.face_info(cand.n_f);
            let a_mif_on = a_fi.vertices_on.clone();
            let a_mif_in = a_fi.vertices_in.clone();
            drop(a_fi);
            // OCCT L388-394: bLinePlane
            let b_line_plane = {
                let c = self.ds.edge_curve(cand.n_e);
                let s = self.ds.face_surface(cand.n_f);
                matches!(c, Some(rcad_kernel::geom::Curve3::Line(_)))
                    && matches!(s, Some(rcad_kernel::geom::Surface3::Plane(_)))
            };

            for &(r1_first, r1_last, is_edge) in &common_parts {
                if is_edge {
                    // OCCT L545-565: EDGE case
                    a_mi_efc.add(cand.n_f);
                    let b_v0 = self.check_face_paves(n_v1, &a_mif_on, &a_mif_in);
                    let b_v1_ = self.check_face_paves(n_v2, &a_mif_on, &a_mif_in);
                    if std::env::var("RCAD_EE_DEBUG").is_ok() {
                        eprintln!("[EF-DBG] EDGE-EF e={} f={} r=[{:.4},{:.4}] bV=({},{}) v1={} v2={}", cand.n_e, cand.n_f, r1_first, r1_last, b_v0, b_v1_, n_v1, n_v2);
                    }
                    self.ds.interf_ef.push(InterferenceEF {
                        edge: cand.n_e, face: cand.n_f,
                        point: self.ds.edge_curve(cand.n_e)
                            .map(|c| c.point_at((r1_first + r1_last) * 0.5)).unwrap_or(glam::DVec3::ZERO),
                        edge_param: (r1_first + r1_last) * 0.5,
                        new_vertex: usize::MAX,
                    });
                    if !b_v0 || !b_v1_ {
                        self.ds.add_interf(cand.n_e, cand.n_f);
                        break;
                    }
                    self.ds.add_interf(cand.n_e, cand.n_f);
                    // OCCT L564: BOPAlgo_Tools::FillMap(aPB, nF, aMPBLI)
                    let ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                    a_mpbli.entry(ptr).or_insert_with(|| (a_pb.clone(), Vec::new())).1.push(cand.n_f);
                } else {
                    // OCCT L406-543: VERTEX case
                    let a_t = (r1_first + r1_last) * 0.5;
                    // OCCT L412: IntTools_Tools::VertexParameter(aCPart, aT) = 0.5*(Range1 first+last)
                    // OCCT L413: BOPTools_AlgoTools::MakeNewVertex(aE, aT, aF, aVnew)
                    //   (BOPTools_AlgoTools_2.cxx L254-271): aPnt = PointOnEdge(aE, aT);
                    //   aMaxTol = Tol(aE) + Tol(aF) + DTolerance() (=1e-12)
                    let a_pnew = self.ds.edge_curve(cand.n_e)
                        .map(|c| c.point_at(a_t)).unwrap_or(glam::DVec3::ZERO);
                    let mut a_tol_vnew = cand.tol_e + cand.tol_f + 1e-12;
                    // OCCT L415-419: IsInRange
                    let a_tol_to_decide = 5e-8;
                    let b_is_on_pave_0 = crate::bop::tools::algo_tools::is_in_range(
                        a_r1.0, a_r1.1, r1_first, r1_last, a_tol_to_decide);
                    let b_is_on_pave_1 = crate::bop::tools::algo_tools::is_in_range(
                        a_r2.0, a_r2.1, r1_first, r1_last, a_tol_to_decide);
                    // OCCT L421-440
                    if (b_is_on_pave_0 && b_is_on_pave_1)
                        || (b_line_plane && (b_is_on_pave_0 || b_is_on_pave_1))
                    {
                        let b_v0 = self.check_face_paves(n_v1, &a_mif_on, &a_mif_in);
                        let b_v1_ = self.check_face_paves(n_v2, &a_mif_on, &a_mif_in);
                        if b_v0 && b_v1_ {
                            self.ds.interf_ef.push(InterferenceEF {
                                edge: cand.n_e, face: cand.n_f,
                                point: a_pnew, edge_param: a_t, new_vertex: usize::MAX,
                            });
                            self.ds.add_interf(cand.n_e, cand.n_f);
                            a_mi_efc.add(cand.n_f);
                            let ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                            a_mpbli.entry(ptr).or_insert_with(|| (a_pb.clone(), Vec::new())).1.push(cand.n_f);
                            break;
                        }
                    }
                    // OCCT L442-445
                    if !b_is_pb_splittable {
                        continue;
                    }
                    // OCCT L447-457: ForceInterfVF for on-pave vertices
                    let mut b_is_on_pave = [b_is_on_pave_0, b_is_on_pave_1];
                    let n_vs = [n_v1, n_v2];
                    for j in 0..2 {
                        if b_is_on_pave[j] {
                            let b_vj = self.check_face_paves(n_vs[j], &a_mif_on, &a_mif_in);
                            if !b_vj {
                                b_is_on_pave[j] = self.force_interf_vf(n_vs[j], cand.n_f);
                            }
                        }
                    }
                    if b_is_on_pave[0] || b_is_on_pave[1] {
                        // OCCT L459-503: check real intersection
                        let a_min_dist_ef = {
                            let proj = self.my_context.proj_ps(&self.ds, cand.n_f);
                            proj.perform(a_pnew);
                            if proj.nb_points() > 0 { proj.lower_distance() } else { f64::MAX }
                        };
                        let has_real_intersection = a_min_dist_ef < rcad_kernel::precision::INTERSECTION;
                        if !has_real_intersection {
                            continue;
                        }
                        for j in 0..2 {
                            if b_is_on_pave[j] {
                                let a_p = self.ds.vertex_point_by_idx(n_vs[j]);
                                let a_dist_pp = (a_pnew - a_p).length();
                                let a_tol = self.ds.vertex_tolerance_by_idx(n_vs[j]);
                                let mut a_max_dist = 1e4 * a_tol;
                                if a_tol < 0.01 {
                                    a_max_dist = a_max_dist.min(0.1);
                                }
                                if a_dist_pp < a_max_dist {
                                    self.update_vertex(n_vs[j], a_dist_pp);
                                    self.my_verts_to_avoid_extension.insert(n_vs[j]);
                                }
                            }
                        }
                        continue;
                    }
                    // OCCT L505-508: CheckFacePaves(aVnew, aMIFOn)
                    if self.check_face_paves_vertex(a_pnew, a_tol_vnew, &a_mif_on) {
                        continue;
                    }
                    // OCCT L510-519: tolerance boost
                    a_tol_vnew = a_tol_vnew.max(cand.tol_e).max(cand.tol_f);
                    if b_line_plane {
                        a_tol_vnew = a_tol_vnew.max((r1_last - r1_first) / 2.0);
                    }
                    // OCCT L523-526: myContext->IsPointInFace(aPnew, aF, aTolVnew)
                    // (IntTools_Context.cxx L613-636): project, require
                    // LowerDistance() < aTolVnew, then 2D IsPointInFace (ON excluded).
                    let proj = self.my_context.proj_ps(&self.ds, cand.n_f);
                    proj.perform(a_pnew);
                    let a_dist = proj.lower_distance();
                    let (a_u, a_v) = proj.lower_distance_parameters();
                    let in_face = proj.nb_points() > 0
                        && a_dist < a_tol_vnew
                        && self.my_context.is_point_in_face(&self.ds, cand.n_f, glam::DVec2::new(a_u, a_v));
                    if std::env::var("RCAD_EE_DEBUG").is_ok() {
                        eprintln!("[EF-DBG] VERTEX case e={} f={} t={:.4} uv=({:.4},{:.4}) dist={:.3e} inFace={} splittable={}", cand.n_e, cand.n_f, a_t, a_u, a_v, a_dist, in_face, b_is_pb_splittable);
                    }
                    if !in_face {
                        continue;
                    }
                    // OCCT L528-542: add InterfEF
                    a_mi_efc.add(cand.n_f);
                    if std::env::var("RCAD_EE_DEBUG").is_ok() {
                        eprintln!("[EF-DBG] VERTEX-EF e={} f={} t={:.4}", cand.n_e, cand.n_f, a_t);
                    }
                    let idx_interf = self.ds.interf_ef.len();
                    self.ds.interf_ef.push(InterferenceEF {
                        edge: cand.n_e, face: cand.n_f,
                        point: a_pnew, edge_param: a_t, new_vertex: usize::MAX,
                    });
                    self.ds.add_interf(cand.n_e, cand.n_f);
                    a_mvcpb.push(CoupleOfPBs {
                        pb1: a_pb.clone(),
                        pb2: a_pb.clone(),
                        index_interf: idx_interf,
                        tolerance: a_tol_vnew,
                        index: usize::MAX,
                        point: a_pnew,
                    });
                }
            }
        }

        // OCCT L576-578: post treatment
        // OCCT BOPAlgo_Tools::PerformCommonBlocks (2nd overload, L191-244):
        //   one CommonBlock per section PaveBlock, reusing the existing CB if the
        //   PB already belongs to one; append the PB's face list.
        if !a_mpbli.is_empty() {
            let pending: Vec<(SharedPB, Vec<usize>)> = a_mpbli.values().cloned().collect();
            for (pb, faces) in &pending {
                // OCCT L206-213: reuse existing CB or create a single-PB CB
                let existing = pb.0.read().unwrap().common_block_idx;
                let cb_idx = match existing {
                    Some(idx) => idx,
                    None => self.ds.add_common_block(std::slice::from_ref(pb)),
                };
                // OCCT L216-238: append new faces (dedup inside append_faces)
                self.ds.common_blocks[cb_idx].append_faces(faces);
            }
        }
        self.update_vertices_of_cb();
        if !a_mvcpb.is_empty() {
            self.perform_new_vertices(&a_mvcpb, false);
            if self.has_errors() {
                return;
            }
        }
        // OCCT L585: Update FaceInfoIn for all faces having EF common parts.
        // OCCT L598: myDS->UpdateFaceInfoIn(aMIEFC) 鈥?NCollection_Map<int>
        // bucket iteration order.
        for n_f in a_mi_efc.iter_keys() {
            self.ds.update_face_info_in(n_f);
        }
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            for (i, si) in self.ds.shapes.iter().enumerate() {
                if si.shape_type == ShapeType::Edge {
                    let pbs = self.ds.edge_pave_blocks(i);
                    if !pbs.is_empty() {
                        let desc: Vec<String> = pbs.iter().map(|pb| {
                            let r = pb.0.read().unwrap();
                            format!("[{:.3},{:.3}]", r.range().0, r.range().1)
                        }).collect();
                        eprintln!("[EF-DBG] edge {} deg={} flag={} pbs={} {}", i, self.ds.is_edge_degenerated(i), self.ds.shapes[i].has_flag(), pbs.len(), desc.join(" "));
                    }
                }
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::CheckFacePaves (PaveFiller_5.cxx L596-601).
    fn check_face_paves(&self, n_vx: usize, a_mif_on: &indexmap::IndexSet<usize>, a_mif_in: &indexmap::IndexSet<usize>) -> bool {
        a_mif_on.contains(&n_vx) || a_mif_in.contains(&n_vx)
    }

    /// OCCT BOPAlgo_PaveFiller::CheckFacePaves(TopoDS_Vertex, Map) (L605-627).
    /// The new vertex's tolerance is aTolVnew (from MakeNewVertex), matching
    /// OCCT's BOPTools_AlgoTools::ComputeVV(aVnew, aV) which uses both tolerances.
    fn check_face_paves_vertex(&self, a_vnew: glam::DVec3, a_tol_vnew: f64, a_mif: &indexmap::IndexSet<usize>) -> bool {
        for &n_v in a_mif {
            let a_p = self.ds.vertex_point_by_idx(n_v);
            let a_tol = self.ds.vertex_tolerance_by_idx(n_v);
            let i_flag = crate::bop::tools::algo_tools::compute_vv_vertex_point(a_tol, a_p, a_vnew, a_tol_vnew);
            if i_flag == 0 {
                return true;
            }
        }
        false
    }

    /// OCCT BOPAlgo_PaveFiller::ForceInterfVF (PaveFiller_5.cxx L631-680).
    fn force_interf_vf(&mut self, n_v: usize, n_f: usize) -> bool {
        let (i_flag, a_u, a_v, a_tol_v_new) = self.my_context.compute_vf(n_v, n_f, &self.ds, self.my_fuzzy_value);
        if i_flag == 0 || i_flag == -2 {
            let mut a_vf = InterferenceVF { vertex: n_v, face: n_f, u: a_u, v: a_v, index_new: None };
            self.ds.add_interf(n_v, n_f);
            let n_vx = self.update_vertex(n_v, a_tol_v_new);
            if self.ds.is_new_shape(n_vx) {
                a_vf.index_new = Some(n_vx);
            }
            self.ds.interf_vf.push(a_vf);
            self.ds.change_face_info(n_f).vertices_in.insert(n_vx);
            return true;
        }
        false
    }

    /// OCCT BOPAlgo_PaveFiller::ReduceIntersectionRange (PaveFiller_5.cxx L685-768).
    fn reduce_intersection_range(
        &self, n_v1: usize, n_v2: usize, n_e: usize, n_f: usize,
        a_ts1: &mut f64, a_ts2: &mut f64,
    ) {
        if !self.ds.is_new_shape(n_v1) && !self.ds.is_new_shape(n_v2) {
            return;
        }
        if !self.ds.has_interf_shape_sub_shapes(n_e, n_f, true) {
            return;
        }
        let a_nb_ees = self.ds.interf_ee.len();
        if a_nb_ees == 0 {
            return;
        }
        // face's edges
        let mut a_mfe: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let a_li = self.ds.shapes[n_f].sub_shapes.clone();
        for n_e1 in a_li {
            if n_e1 < self.ds.nb_shapes() && self.ds.shapes[n_e1].shape_type == ShapeType::Edge {
                a_mfe.insert(n_e1);
            }
        }
        for a_ee in &self.ds.interf_ee {
            if a_ee.new_vertex == usize::MAX {
                continue;
            }
            let n_v = a_ee.new_vertex;
            if n_v != n_v1 && n_v != n_v2 {
                continue;
            }
            let (n_e1, n_e2) = (a_ee.e1, a_ee.e2);
            if ((n_e != n_e1) && (n_e != n_e2)) || (!a_mfe.contains(&n_e1) && !a_mfe.contains(&n_e2)) {
                continue;
            }
            let a_crange = if n_e == n_e1 { a_ee.range1 } else { a_ee.range2 };
            let (a_tr1, a_tr2) = (a_crange[0], a_crange[1]);
            if n_v == n_v1 {
                if *a_ts1 < a_tr2 {
                    *a_ts1 = a_tr2;
                }
            } else {
                if *a_ts2 > a_tr1 {
                    *a_ts2 = a_tr1;
                }
            }
        }
    }

    // ====================================================================
    // FF 鈥?OCCT BOPAlgo_PaveFiller_6.cxx L285-end
    // ====================================================================

    /// OCCT BOPAlgo_PaveFiller::CheckPlanes (BOPAlgo_PaveFiller_6.cxx L3639-3675):
    ///   returns true when the two plane faces share more than one boundary
    ///   vertex (VerticesIn/VerticesOn), i.e. the planes are really interfering.
    fn check_planes(&self, n_f1: usize, n_f2: usize) -> bool {
        let a_fi1 = self.ds.face_info(n_f1);
        let a_fi2 = self.ds.face_info(n_f2);
        let a_mv_in1 = &a_fi1.vertices_in;
        let a_mv_on1 = &a_fi1.vertices_on;
        let mut i_cnt = 0usize;
        let mut b_to_intersect = false;
        for i in 0..2 {
            if b_to_intersect {
                break;
            }
            let a_mv2 = if i == 0 { &a_fi2.vertices_in } else { &a_fi2.vertices_on };
            for &n_v2 in a_mv2 {
                if a_mv_in1.contains(&n_v2) || a_mv_on1.contains(&n_v2) {
                    i_cnt += 1;
                    if i_cnt > 1 {
                        b_to_intersect = true;
                        break;
                    }
                }
            }
        }
        b_to_intersect
    }

    fn perform_ff(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L299: myIterator->Initialize(FACE, FACE); iSize = ExpectedLength()
        let pairs: Vec<(usize, usize)> = if let Some(it) = &self.my_iterator {
            it.pairs(ShapeType::Face, ShapeType::Face).to_vec()
        } else {
            return;
        };

        // OCCT L294-313: collect faces from the intersection pairs and the rest
        // of the touched faces, then refresh the FaceInfo for all of them.
        let mut a_mi_fence: indexmap::IndexSet<usize> = indexmap::IndexSet::new();
        for &(n_f1, n_f2) in &pairs {
            a_mi_fence.insert(n_f1);
            a_mi_fence.insert(n_f2);
        }
        for i in 0..self.ds.nb_source_shapes() {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type == ShapeType::Face && a_si.has_reference() {
                a_mi_fence.insert(i);
            }
        }
        // myDS->UpdateFaceInfoOn(aMIFence); myDS->UpdateFaceInfoIn(aMIFence);
        for &f in &a_mi_fence {
            let _ = self.ds.change_face_info(f);
            self.ds.update_face_info_on(f);
            self.ds.update_face_info_in(f);
        }
        if pairs.is_empty() { return; }

        // OCCT L328-333: options for the intersection algorithm
        let b_approx = self.my_section_attribute.approximation;
        let b_comp_c2d1 = self.my_section_attribute.pcurve_on_s1;
        let b_comp_c2d2 = self.my_section_attribute.pcurve_on_s2;
        let an_approx_tol = 1.0e-7;
        // Post-processing options
        let b_split_curve = false;

        let mut new_ff: Vec<InterferenceFF> = Vec::new();

        // OCCT L336-360: collect all pairs of Edge/Edge interferences to check
        // if some faces have to be moved to obtain more precise intersection.
        let mut a_ee_map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for a_ee in &self.ds.interf_ee {
            if a_ee.new_vertex != usize::MAX {
                // OCCT HasIndexNew() -> IndexNew().  BOPDS_Pair orders the pair.
                let pair = if a_ee.e1 <= a_ee.e2 { (a_ee.e1, a_ee.e2) } else { (a_ee.e2, a_ee.e1) };
                a_ee_map.entry(pair).or_default().push(a_ee.new_vertex);
            }
        }

        for &(i, j) in &pairs {
            // OCCT L513-519: for the Glue mode just add all interferences of
            // that type (empty FF), without performing the intersection.
            if self.my_glue != GlueEnum::GlueOff {
                new_ff.push(InterferenceFF {
                    f1: i, f2: j,
                    curves: Vec::new(),
                    points: Vec::new(),
                    tangent_faces: false,
                });
                continue;
            }
            let Some(s1) = self.ds.face_surface(i) else { continue; };
            let Some(s2) = self.ds.face_surface(j) else { continue; };

            // OCCT L373-391: check if the planes are really interfering (share
            // more than one boundary vertex); otherwise skip the pair.
            if matches!((&s1, &s2), (Surface3::Plane(_), Surface3::Plane(_))) {
                if !self.check_planes(i, j) {
                    new_ff.push(InterferenceFF {
                        f1: i, f2: j,
                        curves: Vec::new(),
                        points: Vec::new(),
                        tangent_faces: false,
                    });
                    continue;
                }
            }

            // OCCT L400-486: check if there is an intersection between the edges
            // of the faces; if the intersection point is at some distance from
            // the edges, move one of the faces to the point of exact edge
            // intersection (only closed/seam edges are considered).
            let mut a_f_shifted1 = s1.clone();
            let mut a_f_shifted2 = s2.clone();
            let mut a_shift_value = 0.0;
            if !is_plane_ff(&s1) || !is_plane_ff(&s2) {
                let an_is_plane1 = is_plane_ff(&s1);
                let an_is_plane2 = is_plane_ff(&s2);
                let edges1 = face_edge_indices(&self.ds, i);
                let edges2 = face_edge_indices(&self.ds, j);
                'outer: for &n_e1 in &edges1 {
                    let Some(an_edge1) = self.ds.shape(n_e1).as_edge() else { continue; };
                    let an_is_closed1 = is_closed_ff(&self.ds, i, n_e1, an_is_plane1);
                    for &n_e2 in &edges2 {
                        let Some(an_edge2) = self.ds.shape(n_e2).as_edge() else { continue; };
                        let an_is_closed2 = is_closed_ff(&self.ds, j, n_e2, an_is_plane2);
                        if !an_is_closed1 && !an_is_closed2 {
                            continue;
                        }
                        let pair = if n_e1 <= n_e2 { (n_e1, n_e2) } else { (n_e2, n_e1) };
                        let Some(a_vertex_indices) = a_ee_map.get(&pair) else { continue; };
                        for &a_vertex_index in a_vertex_indices {
                            let Some(a_vertex) = self.ds.shape(a_vertex_index).as_vertex() else { continue; };
                            let a_vertex_point = a_vertex.point;
                            let Some(a_c1) = an_edge1.curve.clone() else { continue; };
                            let Some(a_c2) = an_edge2.curve.clone() else { continue; };
                            // OCCT L457-460: compute points exactly on the edges.
                            let a_proj1 = closest_point_on_curve_range(
                                &a_c1, a_vertex_point, an_edge1.range[0], an_edge1.range[1], 64);
                            let a_proj2 = closest_point_on_curve_range(
                                &a_c2, a_vertex_point, an_edge2.range[0], an_edge2.range[1], 64);
                            let a_p1 = a_proj1.point;
                            let a_p2 = a_proj2.point;
                            let a_shift_dist = a_p1.distance(a_p2);
                            if a_shift_dist > a_vertex.tolerance {
                                // OCCT L474-478: move one of the faces to the
                                // point of exact intersection of the edges.
                                if an_is_closed1 {
                                    a_f_shifted1 = translate_surface(&a_f_shifted1, a_p2 - a_p1);
                                } else {
                                    a_f_shifted2 = translate_surface(&a_f_shifted2, a_p1 - a_p2);
                                }
                                a_shift_value = a_shift_dist;
                                break 'outer;
                            }
                        }
                    }
                }
            }

            let uv1 = self.ds.face_actual_uv_bounds(i);
            let uv2 = self.ds.face_actual_uv_bounds(j);
            // OCCT L496-497: aTolFF = max(aShiftValue, ToleranceFF(aBAS1, aBAS2)).
            let a_tol_ff = a_shift_value
                .max(tolerance_ff(&s1, &s2, self.ds.face_tolerance(i), self.ds.face_tolerance(j)));
            let mut ff = int_tools::face_face::FaceFace::new();
            // OCCT L491-495: aFaceFace.SetRunParallel(myRunParallel);
            // SetIndices(nF1, nF2); SetFaces(aFShifted1, aFShifted2);
            // SetBoxes(myDS->ShapeInfo(nF1).Box(), myDS->ShapeInfo(nF2).Box()).
            ff.set_surfaces(a_f_shifted1, a_f_shifted2);
            ff.set_uv_bounds(uv1, uv2);
            ff.set_face_indices(i, j);
            ff.set_boxes(self.ds.shape_info(i).bbox.clone(), self.ds.shape_info(j).bbox.clone());
            // OCCT L498: SetTolFF(aTolFF).
            ff.set_tol_ff(a_tol_ff);
            ff.set_tolerances(self.ds.face_tolerance(i), self.ds.face_tolerance(j));
            // OCCT L500-506: GetEFPnts(nF1, nF2, aListOfPnts); if (aNbLP) SetList.
            let a_list_of_pnts = self.get_ef_pnts(i, j);
            if !a_list_of_pnts.is_empty() {
                ff.set_list(a_list_of_pnts);
            }
            // OCCT L508-510: SetParameters(bApprox, bCompC2D1, bCompC2D2, anApproxTol).
            ff.set_parameters(b_approx, b_comp_c2d1, b_comp_c2d2, an_approx_tol);
            ff.set_fuzzy_value(self.my_fuzzy_value);
            ff.perform(&self.ds);
            let tangent_faces = ff.tangent_faces();
            // OCCT L543-556: if (!aFaceFace.IsDone() || aFaceFace.HasErrors())
            // 鈫?empty FF + AddIntersectionFailedWarning; rcad FaceFace has no
            // error channel 鈥?a failed Perform yields no curves.
            if !ff.has_intersection() {
                // OCCT L545-552: aFF.SetIndices; aFF.Init(0, 0);
                // AddIntersectionFailedWarning(Face1(), Face2()).
                self.my_report.add_warning(Alert::IntersectionFailed(i, j));
                // OCCT: empty FF interference is still added for the pair
                new_ff.push(InterferenceFF {
                    f1: i, f2: j,
                    curves: Vec::new(),
                    points: Vec::new(),
                    tangent_faces,
                });
                continue;
            }
            // OCCT L558-563: aFaceFace.PrepareLines3D(bSplitCurve);
            // aFaceFace.ApplyTrsf().
            ff.prepare_lines_3d(b_split_curve);
            ff.apply_trsf();
            let curves = ff.make_curves();
            let a_nb_curves = curves.len();
            // OCCT L565-572: aCvsX = aFaceFace.Lines(); aPntsX = aFaceFace.Points();
            // aNbCurves/aNbPoints; if (aNbCurves || aNbPoints) myDS->AddInterf(nF1, nF2).
            // rcad: points 鐢?IntTools_FaceFace 鐢熸垚锛況cad FaceFace 鏆備笉浜у嚭瀛ょ珛鐐广€?
            if a_nb_curves > 0 {
                self.ds.add_interf(i, j);
            }
            // OCCT L578-581: aFF.SetIndices(nF1, nF2); SetTangentFaces(bTangentFaces);
            // Init(aNbCurves, aNbPoints).
            // OCCT L583-590: aBoxExpandValue = aTolFF (+ max vertex tol when curves exist).
            let a_max_vertex_tol = if a_nb_curves > 0 {
                self.ds.face_max_vertex_tolerance(i)
                    .max(self.ds.face_max_vertex_tolerance(j))
            } else {
                0.0
            };
            let a_box_expand_value = a_tol_ff + a_max_vertex_tol;
            let mut curve_ids: Vec<usize> = Vec::new();
            for mut c in curves {
                // OCCT L599-607: bIsValid = CheckCurve(aIC, aBox); if valid 鈫?
                // SetCurve(aIC); aBox.Enlarge(aBoxExpandValue); SetBox(aBox);
                // SetTolerance(max(aIC.Tolerance(), aTolFF)).
                let (b_is_valid, a_box) = int_tools::face_face::check_curve(&c);
                if !b_is_valid {
                    continue;
                }
                if c.tolerance < a_tol_ff {
                    c.tolerance = a_tol_ff;
                }
                if let Some([mn, mx]) = a_box {
                    let grow = glam::DVec3::splat(a_box_expand_value);
                    c.bbox = Some((mn - grow, mx + grow));
                }
                let cid = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(c);
                curve_ids.push(cid);
            }
            // OCCT L549-550: aFaceFace.Indices(nF1, nF2) 鈥?myIF1/myIF2, the
            // ORIGINAL pair order (BOPAlgo_FaceFace does not swap the indices;
            // only IntTools_FaceFace swaps its internal Face1/Face2 and the
            // pcurves are exchanged back afterwards, IntTools_FaceFace.cxx
            // L550-562, so pcurve1 always belongs to the original nF1).
            new_ff.push(InterferenceFF {
                f1: i, f2: j,
                curves: curve_ids,
                points: Vec::new(),
                tangent_faces,
            });
        }
        self.ds.interf_ff.extend(new_ff);
    }

    /// OCCT BOPAlgo_PaveFiller::GetEFPnts (BOPAlgo_PaveFiller_6.cxx L2665-2740).
    ///
    /// Collects the Edge-Face intersection points belonging to both faces
    /// nF1/nF2 as IntSurf_PntOn2S (UV on each face), to be passed to
    /// IntTools_FaceFace::SetList 鈥?the intersection curves are then started
    /// from these points (used by the Param-Param intersector).
    fn get_ef_pnts(&mut self, n_f1: usize, n_f2: usize) -> Vec<PntOn2S> {
        // OCCT L2673-2676: collect indexes of all shapes from nF1 and nF2.
        let mut a_mi: HashSet<usize> = HashSet::new();
        a_mi.insert(n_f1);
        a_mi.extend(self.ds.shape_info(n_f1).sub_shapes.iter().copied());
        a_mi.insert(n_f2);
        a_mi.extend(self.ds.shape_info(n_f2).sub_shapes.iter().copied());
        //
        let mut a_list_of_pnts: Vec<PntOn2S> = Vec::new();
        let a_nb_efs = self.ds.interf_ef.len();
        for i in 0..a_nb_efs {
            let a_ef = self.ds.interf_ef[i].clone();
            // OCCT L2682: if (aEF.HasIndexNew())
            if a_ef.new_vertex == usize::MAX {
                continue;
            }
            let n_e = a_ef.edge;
            let n_f_opposite = a_ef.face;
            // OCCT L2684: if (aMI.Contains(nE) && aMI.Contains(nFOpposite))
            if !(a_mi.contains(&n_e) && a_mi.contains(&n_f_opposite)) {
                continue;
            }
            // OCCT L2686-2688: aPar = aCP.VertexParameter1(); the edge 3D curve.
            let a_par = a_ef.edge_param;
            let a_curve = match self.ds.edge_curve(n_e) {
                Some(c) => c.clone(),
                None => continue,
            };
            //
            // OCCT L2690-2691: nF = (nFOpposite == nF1) ? nF2 : nF1
            let n_f = if n_f_opposite == n_f1 { n_f2 } else { n_f1 };
            //
            // OCCT L2692-2694: aPCurve = BRep_Tool::CurveOnSurface(aE, aF, f, l)
            let a_pcurve = crate::topalgo::shape_source::edge_pcurve_on_face(
                &self.ds,
                n_e,
                n_f,
                topods::Orientation::Forward,
            );
            //
            // OCCT L2696: GeomAPI_ProjectPointOnSurf& aProj = myContext->ProjPS(aFOpposite)
            //
            // OCCT L2698: aCurve->D0(aPar, aPoint)
            let a_point = a_curve.point_at(a_par);
            let mut a_pnt = PntOn2S::new();
            // OCCT L2696-2736: project aPoint onto the opposite face (and, when
            // the edge has no pcurve on the other face, onto both faces). The
            // projectors are obtained inside short blocks to keep the mutable
            // borrow of my_context exclusive.
            let a_u_v_opp = {
                let a_proj = self.my_context.proj_ps(&self.ds, n_f_opposite);
                a_proj.perform(a_point);
                if a_proj.nb_points() > 0 {
                    Some(a_proj.lower_distance_parameters())
                } else {
                    None
                }
            };
            if let Some((a_pc, _, _)) = a_pcurve {
                // OCCT L2701: aP2d = aPCurve->Value(aPar)
                let a_p2d = a_pc.point_at(a_par);
                if let Some((u1, v1)) = a_u_v_opp {
                    // OCCT L2704: aProj.LowerDistanceParameters(U1, V1)
                    if n_f == n_f1 {
                        a_pnt.set_value(a_p2d.x, a_p2d.y, u1, v1);
                    } else {
                        a_pnt.set_value(u1, v1, a_p2d.x, a_p2d.y);
                    }
                    a_list_of_pnts.push(a_pnt);
                }
            } else if let Some((u2, v2)) = a_u_v_opp {
                // OCCT L2716-2736: no pcurve 鈥?project onto both faces.
                let a_u_v_1 = {
                    let a_proj1 = self.my_context.proj_ps(&self.ds, n_f);
                    a_proj1.perform(a_point);
                    if a_proj1.nb_points() > 0 {
                        Some(a_proj1.lower_distance_parameters())
                    } else {
                        None
                    }
                };
                if let Some((u1, v1)) = a_u_v_1 {
                    if n_f == n_f1 {
                        a_pnt.set_value(u1, v1, u2, v2);
                    } else {
                        a_pnt.set_value(u2, v2, u1, v1);
                    }
                    a_list_of_pnts.push(a_pnt);
                }
            }
        }
        a_list_of_pnts
    }

    // ====================================================================
    // OCCT BOPAlgo_PaveFiller sub-steps
    // ====================================================================

    /// OCCT: UpdatePaveBlocksWithSDVertices 鈥?delegates to DS.
    fn update_pave_blocks_with_sd_vertices(&mut self) {
        self.ds.update_pave_blocks_with_sd_vertices();
    }

    /// OCCT BOPAlgo_PaveFiller::UpdateInterfsWithSDVertices (_10.cxx L248-255).
    fn update_interfs_with_sd_vertices(&mut self) {
        self.update_vv_sd();
        self.update_ve_sd();
        self.update_vf_sd();
        self.update_ee_sd();
        self.update_ef_sd();
    }

    fn update_vv_sd(&mut self) {
        let idx: Vec<usize> = self.ds.interf_vv.iter().enumerate()
            .filter_map(|(i, vv)| {
                if vv.merged_vertex != usize::MAX {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(vv.merged_vertex, &mut sd) { Some(i) } else { None }
                } else { None }
            }).collect();
        for &i in &idx {
            let mut sd = usize::MAX;
            if self.ds.has_shape_sd(self.ds.interf_vv[i].merged_vertex, &mut sd) {
                self.ds.interf_vv[i].merged_vertex = sd;
            }
        }
    }

    fn update_ve_sd(&mut self) {
        let idx: Vec<usize> = self.ds.interf_ve.iter().enumerate()
            .filter_map(|(i, ve)| {
                if ve.index_new != 0 {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(ve.index_new, &mut sd) { Some(i) } else { None }
                } else { None }
            }).collect();
        for &i in &idx {
            let mut sd = usize::MAX;
            if self.ds.has_shape_sd(self.ds.interf_ve[i].index_new, &mut sd) {
                self.ds.interf_ve[i].index_new = sd;
            }
        }
    }

    fn update_vf_sd(&mut self) {
        let idx: Vec<(usize, usize)> = self.ds.interf_vf.iter().enumerate()
            .filter_map(|(i, vf)| {
                vf.index_new.and_then(|nv| {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(nv, &mut sd) { Some((i, sd)) } else { None }
                })
            }).collect();
        for (i, sd) in idx {
            self.ds.interf_vf[i].index_new = Some(sd);
        }
    }

    fn update_ee_sd(&mut self) {
        let idx: Vec<(usize, usize)> = self.ds.interf_ee.iter().enumerate()
            .filter_map(|(i, ee)| {
                if ee.new_vertex != usize::MAX {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(ee.new_vertex, &mut sd) { Some((i, sd)) } else { None }
                } else { None }
            }).collect();
        for (i, sd) in idx {
            self.ds.interf_ee[i].new_vertex = sd;
        }
    }

    fn update_ef_sd(&mut self) {
        let idx: Vec<(usize, usize)> = self.ds.interf_ef.iter().enumerate()
            .filter_map(|(i, ef)| {
                if ef.new_vertex != usize::MAX {
                    let mut sd = usize::MAX;
                    if self.ds.has_shape_sd(ef.new_vertex, &mut sd) { Some((i, sd)) } else { None }
                } else { None }
            }).collect();
        for (i, sd) in idx {
            self.ds.interf_ef[i].new_vertex = sd;
        }
    }

    /// OCCT BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices (_6.cxx L3946-4020).
    fn update_blocks_with_shared_vertices(&mut self) {
        // OCCT L3948-3951: only active in non-destructive mode
        if !self.my_non_destructive {
            return;
        }
        // L3955-3960: if no FF interferences, return
        if self.ds.interf_ff.is_empty() {
            return;
        }
        // OCCT L3967-4020: iterate FF interferences, build shared vertex sets
        // rcad: non-destructive mode is not fully implemented.
    }

    /// OCCT BOPDS_DS::RefineFaceInfoIn (BOPDS_DS.cxx L995-1024).
    fn refine_face_info_in(&mut self) {
        let n = self.ds.nb_source_shapes();
        for i in 0..n {
            let si = self.ds.shape_info(i);
            if si.shape_type != ShapeType::Face || !si.has_reference() { continue; }
            let pb_on = self.ds.face_info(i).pave_blocks_on.clone();
            let pb_in = self.ds.face_info(i).pave_blocks_in.clone();
            if pb_in.is_empty() || pb_on.is_empty() { continue; }
            let mut to_rem: Vec<u64> = Vec::new();
            for &pb in &pb_in { if pb_on.contains(&pb) { to_rem.push(pb); } }
            let fi = self.ds.change_face_info(i);
            for &r in &to_rem { fi.pave_blocks_in.swap_remove(&r); }
        }
    }

    /// OCCT BOPDS_DS::RefineFaceInfoOn (BOPDS_DS.cxx L975-991).
    fn refine_face_info_on(&mut self) {
        // OCCT BOPDS_DS::RefineFaceInfoOn (BOPDS_DS.cxx L975-991): rebuild each
        // face's PaveBlocksOn from the current edge PBs (UpdateFaceInfoOn),
        // then drop the released PBs (those whose edge was cleared by
        // ReleasePaveBlocks and no longer carries an edge).
        for i in 0..self.ds.face_info_pool.len() {
            let idx = self.ds.face_info_pool[i].index();
            // OCCT: UpdateFaceInfoOn(aFaceInfo.Index()) 鈥?re-collect.
            self.ds.update_face_info_on(idx);
            let pb_on = self.ds.face_info(idx).pave_blocks_on.clone();
            let mut to_rem: Vec<u64> = Vec::new();
            for &pb in &pb_on {
                let has = self.ds.pb_from_ptr(pb)
                    .map_or(false, |p| p.0.read().unwrap().edge != usize::MAX);
                if !has { to_rem.push(pb); }
            }
            if !to_rem.is_empty() {
                let fi = self.ds.change_face_info(idx);
                for &r in &to_rem { fi.pave_blocks_on.swap_remove(&r); }
            }
        }
    }

    // OCCT BOPAlgo_PaveFiller::Init (PaveFiller.cxx L176-214).
    fn init(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L178-182: Check arguments non-empty
        if self.my_arguments.is_empty() && self.ds.nb_source_shapes() == 0 {
            self.my_report.add_error(Alert::TooFewArguments);
            return;
        }
        // OCCT L184: Message_ProgressScope aPS(theRange, "Initialization of Intersection algorithm", 1);
        // rcad: aPS covers the null-shape-check loop (1 step), skipped here (Rust Shape prevents null).
        let _a_ps = the_range.sub_scope("Initialization of Intersection algorithm", 1);
        // OCCT L185-193: check for null shapes 鈥?Rust Shape type prevents null.
        // OCCT L196: Clear
        self.clear();
        // OCCT L199-201: myDS = new BOPDS_DS;
        //   myDS->SetArguments(myArguments);
        //   myDS->Init(myFuzzyValue);
        if !self.my_arguments.is_empty() {
            self.ds.set_arguments(std::mem::take(&mut self.my_arguments));
        }
        // Only the primary DS deep-clones its arguments (non-destructive
        // contract); nested PaveFillers (PostTreatFF fuse, checker) keep the
        // argument shapes by reference so `DS::index(aSx)` resolves them.
        self.ds.clone_arguments = self.my_is_primary;
        self.ds.init(self.my_fuzzy_value);
        // OCCT L204: myContext = new IntTools_Context
        self.my_context = IntToolsContext::new();
        // OCCT L207-210: myIterator = new BOPDS_Iterator
        let mut a_it = BOPDS_Iterator::new(self.my_fuzzy_value);
        a_it.set_run_parallel(self.my_run_parallel); // OCCT L208
        // OCCT L210: myIterator->Prepare(myContext, myUseOBB, myFuzzyValue)
        a_it.prepare(&self.ds, Some(&self.my_context), false, self.my_fuzzy_value);
        self.my_iterator = Some(Box::new(a_it));
        // OCCT L213: SetNonDestructive 鈥?respects existing flag
        self.set_non_destructive();
    }

    // OCCT BOPAlgo_PaveFiller::SetNonDestructive (PaveFiller_10.cxx L41-59).
    // Checks if any argument shape is locked; if so, enables non-destructive mode.
    fn set_non_destructive(&mut self) {
        if !self.my_is_primary || self.my_non_destructive {
            return;
        }
        // OCCT L47-55: iterate myArguments, check aS.Locked()
        for arg in &self.my_arguments {
            if arg.data.locked() {
                self.my_non_destructive = true;
                return;
            }
        }
    }

    // OCCT BOPAlgo_PaveFiller::Prepare (_7.cxx L850-931).
    fn prepare(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L852-856: non-destructive mode 鈫?skip
        if self.my_non_destructive { return; }

        // OCCT L857-879: iterate (V,F), (E,F), (F,F) pairs,
        // collect planar faces via IsBasedOnPlane
        // OCCT aMF (L859) is NCollection_Map<int> 鈥?bucket iteration order.
        let a_types = [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face];
        let mut a_mf: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();

        if let Some(ref mut it) = self.my_iterator {
            for &a_type in &a_types {
                it.initialize(a_type, ShapeType::Face);
                while it.more() {
                    let (n1, nf) = it.value();
                    // Determine which index is the face
                    let fi = if self.ds.shape_info(n1).shape_type() == ShapeType::Face
                    { n1 } else { nf };
                    // OCCT: IsBasedOnPlane(aF)
                    if is_based_on_plane(self.ds.shape(fi)) {
                        a_mf.add(fi);
                    }
                    it.next();
                }
            }
        }

        // OCCT L881-885: no planar faces 鈫?return
        let a_nb_f = a_mf.len();
        if a_nb_f == 0 { return; }

        // OCCT L888-901: collect edge-face pairs into BPC vector. The face's
        // edges are taken in the stored wire order (TopExp_Explorer aExp(aF,
        // TopAbs_EDGE) 鈥?base edges before located copies).
        let mut a_vbpc: Vec<BOPAlgo_BPC> = Vec::new();
        // OCCT L890-901: iterate aMF (NCollection_Map 鈥?bucket order).
        for fi in a_mf.iter_keys() {
            let edge_indices = face_edge_indices(&self.ds, fi);
            for &ei in &edge_indices {
                let mut a_bpc = BOPAlgo_BPC::new(ei, fi);
                a_vbpc.push(a_bpc);
            }
        }
        // Phase 2 (L903-909): perform all BPCs (OCCT: BOPTools_Parallel::Perform)
        for bpc in &mut a_vbpc {
            bpc.perform(&self.ds);
        }
        // Phase 3 (L916-930): update edges with computed pcurves.
        // In-place (OCCT BRep_Builder semantics): every reference to the edge 鈥?
        // the DS entry and the face-wire edges 鈥?observes the pcurve, preserving
        // the single-object identity that WireSplitter and downstream modules
        // key by. Safe because the DS holds a private copy of the input
        // (clone_arguments_private), so the mutation never leaks into the input.
        for bpc in &a_vbpc {
            if bpc.is_to_update() {
                let ei = bpc.edge_idx();
                let fi = bpc.face_idx();
                let range = self.ds.shape(ei).as_edge().map(|ed| ed.range).unwrap_or([0.0, 0.0]);
                if let Some(pc) = bpc.pcurve().cloned() {
                    // OCCT BOPAlgo_BPC / BRepLib::BuildPCurveForEdgeOnPlane: the
                    // pcurve follows the edge's 3D curve parameter direction
                    // (BRep_Tool::Curve + Range), independent of the stored
                    // vertex order. No reversal.
                    let (pc, f, l) = (pc, range[0], range[1]);
                    // OCCT BRep_Builder::UpdateEdge (BRep_Builder.cxx L692) 鈥?
                    // the pcurve is stored under (face TShape, L.Predivided(E.Location())).
                    let key = self.pcurve_key_for(ei, fi); // EXPERIMENT: composed key
                    self.ds.mutate_shape_data(ei, |ts| {
                        if let topods::TShape::Edge(ed) = ts {
                            if let Some(k) = key {
                                ed.pcurves.insert(k, (pc, f, l));
                            }
                        }
                    });
                    self.ds.remap_shape_idx(ei);
                }
            }
        }
    }

    /// OCCT Geom2d_Curve::Reversed() 鈥?same curve, parameter direction
    /// reversed (ReversedParameter maps t to the original curve's endpoint
    /// parameter).  Only the analytic cases (Line/Circle) produced by
    /// project_on_surface are handled; other types keep the projected
    /// direction.
    fn reverse_curve2d(c: rcad_kernel::geom::Curve2d, len: f64) -> rcad_kernel::geom::Curve2d {
        use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
        match c {
            // p'(t) = p(L - t) with p(t) = origin + dir*t, p(L) = origin + dir*L
            Curve2d::Line(l) => Curve2d::Line(Line2d::new(l.origin + l.direction * len, -l.direction)),
            // p'(t) = x_dir*cos(t) - y_dir*sin(t) = p(2*PI - t)
            Curve2d::Circle(c2) => Curve2d::Circle(Circle2d {
                center: c2.center,
                x_dir: c2.x_dir,
                y_dir: -c2.y_dir,
                radius: c2.radius,
            }),
            other => other,
        }
    }

    // ====================================================================
    // TreatNewVertices 鈥?OCCT BOPAlgo_PaveFiller_3.cxx L692-723
    // ====================================================================

    /// Fuse close new vertices into groups.
    /// OCCT BOPAlgo_PaveFiller::TreatNewVertices (PaveFiller_3.cxx L692-723).
    fn treat_new_vertices(&self, the_mvcpb: &[CoupleOfPBs]) -> Vec<(DVec3, f64, Vec<usize>)> {
        // OCCT L700-706: collect vertex points and tolerances.
        // The fused vertex is the EE/EF INTERSECTION point (the new vertex),
        // not the pave block endpoints.
        let mut verts: Vec<(DVec3, f64, usize)> = Vec::new(); // (point, tol, index)
        for (i, cpb) in the_mvcpb.iter().enumerate() {
            let pt = cpb.point;
            verts.push((pt, cpb.tolerance, i));
        }
        // OCCT L710: BOPAlgo_Tools::IntersectVertices 鈥?fuse by proximity
        // rcad: simple fuse 鈥?group vertices within myFuzzyValue distance
        let mut groups: Vec<(DVec3, f64, Vec<usize>)> = Vec::new();
        let mut assigned = vec![false; verts.len()];
        for i in 0..verts.len() {
            if assigned[i] { continue; }
            let (pi, ti, ii) = verts[i];
            let mut group_pts = vec![pi];
            let mut group_tol = ti;
            let mut group_indices = vec![ii];
            assigned[i] = true;
            for j in (i + 1)..verts.len() {
                if assigned[j] { continue; }
                let (pj, tj, ij) = verts[j];
                // OCCT BOPAlgo_Tools::IntersectVertices (L1094-1098):
                // aTolSum = aTolV1 + aTolV2 + myFuzzyValue.
                if (pi - pj).length() <= ti + tj + self.my_fuzzy_value {
                    group_pts.push(pj);
                    group_tol = group_tol.max(tj);
                    group_indices.push(ij);
                    assigned[j] = true;
                }
            }
            // Average position
            let avg = group_pts.iter().sum::<DVec3>() / group_pts.len() as f64;
            groups.push((avg, group_tol, group_indices));
        }
        groups
    }

    /// Process new vertices from EE/EF intersection: add to DS, update interfs, split PBs.
    /// OCCT BOPAlgo_PaveFiller::PerformNewVertices (PaveFiller_3.cxx L594-688).
    fn perform_new_vertices(&mut self, the_mvcpb: &[CoupleOfPBs], is_ee_intersection: bool) {
        // OCCT L601-605: empty check
        if the_mvcpb.is_empty() {
            return;
        }
        // OCCT L607: aTolAdd = myFuzzyValue / 2.
        let a_tol_add = self.my_fuzzy_value / 2.0;

        // OCCT L609-612: TreatNewVertices 鈥?fuse vertices
        let groups = self.treat_new_vertices(the_mvcpb);

        // OCCT L622-653: add fused vertices to DS, update interference indices
        // Maps each CPB's index_interf 鈫?new DS vertex index
        let mut cpb_to_new_vertex: HashMap<usize, usize> = HashMap::new();
        for (_i, (avg_pt, group_tol, member_indices)) in groups.iter().enumerate() {
            // OCCT L631-639: vertex tolerance = group max tolerance (MakeVertex aNTol);
            // box = Add(Pnt(aV)); SetGap(Tolerance(aV) + aTolAdd).
            let n_v = self.ds.push_vertex(*avg_pt, *group_tol);
            let si = self.ds.change_shape_info(n_v);
            si.bbox = rcad_kernel::math::bnd::BndBox::from_point(*avg_pt);
            si.bbox.set_gap(*group_tol + a_tol_add);

            for &cpb_idx in member_indices {
                if cpb_idx < the_mvcpb.len() {
                    let idx_interf = the_mvcpb[cpb_idx].index_interf;
                    cpb_to_new_vertex.insert(idx_interf, n_v);
                    // OCCT L648-652: update interference's new_vertex
                    if is_ee_intersection {
                        if idx_interf < self.ds.interf_ee.len() {
                            self.ds.interf_ee[idx_interf].new_vertex = n_v;
                        }
                    } else {
                        if idx_interf < self.ds.interf_ef.len() {
                            self.ds.interf_ef[idx_interf].new_vertex = n_v;
                        }
                    }
                }
            }
        }

        // OCCT L655-685: build PB鈫抂vertices] map aMPBLI (IndexedDataMap,
        // insertion order preserved).
        let mut a_mpbli: indexmap::IndexMap<u64, (SharedPB, Vec<usize>)> =
            indexmap::IndexMap::new();
        for cpb in the_mvcpb {
            let n_v = match cpb_to_new_vertex.get(&cpb.index_interf) {
                Some(v) => *v,
                None => continue,
            };
            // OCCT L670-678: aPB[] = {aPB1, aPB2}; append iV to the PB's vertex list
            let pbs = [&cpb.pb1, &cpb.pb2];
            for &a_pb in &pbs {
                let pb_ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                let entry = a_mpbli.entry(pb_ptr).or_insert_with(|| (a_pb.clone(), Vec::new()));
                entry.1.push(n_v);
                // OCCT L680-683: if (aPB[0] == aPB[1]) break;
                if std::sync::Arc::ptr_eq(&cpb.pb1.0, &cpb.pb2.0) {
                    break;
                }
            }
        }
        // OCCT L687: IntersectVE(aMPBLI, aPS.Next(), false)
        self.intersect_ve(&a_mpbli, false);
    }

// ====================================================================
// FillShrunkData 鈥?OCCT BOPAlgo_PaveFiller_9.cxx L65-138
// ====================================================================

// OCCT BOPAlgo_PaveFiller::FillShrunkData (PaveFiller_9.cxx L65-138).
fn fill_shrunk_data(&mut self, a_type1: ShapeType, a_type2: ShapeType) {
        // OCCT L68: myIterator->Initialize(aType1, aType2)
        let my_iterator = match &mut self.my_iterator {
            Some(it) => it,
            None => return,
        };
        my_iterator.initialize(a_type1, a_type2);
        let i_size = my_iterator.pairs(a_type1, a_type2).len();
        if i_size == 0 { return; }

        // OCCT L75-80: locals
        let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let a_type = [a_type1, a_type2];
        let mut a_vsd: Vec<ShrunkRange> = Vec::new();

        // OCCT L82-126: iterate pairs
        let pairs: Vec<(usize, usize)> = my_iterator.pairs(a_type1, a_type2).to_vec();
        for &(ns0, ns1) in &pairs {
            let n_s = [ns0, ns1];
            for i in 0..2 {
                let n_e = n_s[i];
                if a_type[i] != ShapeType::Edge || !a_mi.insert(n_e) { continue; }
                if self.ds.shapes[n_e].has_flag() { continue; }
                // OCCT L100: NCollection_List<handle<PaveBlock>>& aLPB = myDS->ChangePaveBlocks(nE);
                // rcad: clone PBs for borrow-safe DS access inside loop
                let a_lpb: Vec<SharedPB> = {
                    self.ds.init_pave_blocks(n_e);
                    let pb_list = self.ds.change_pave_blocks(n_e);
                    pb_list.iter().map(|pb| pb.clone()).collect()
                };
                for a_pb in &a_lpb {
                    // OCCT L105: if (aPB->HasShrunkData() && myDS->IsValidShrunkData(aPB)) continue;
                    let (n_v1, n_v2, a_t1, a_t2, skip) = {
                        let pbr = a_pb.0.read().unwrap();
                        if pbr.has_shrunk_data() && self.ds.is_valid_shrunk_data(&*pbr) {
                            (0, 0, 0.0, 0.0, true)
                        } else {
                            let (v1, v2) = pbr.indices();
                            let (t1, t2) = pbr.range();
                            (v1, v2, t1, t2, false)
                        }
                    };
                    if skip { continue; }
                    a_vsd.push(ShrunkRange::new(a_pb, n_v1, n_v2, a_t1, a_t2));
                }
            }
        }
        // OCCT L128-137: Perform + AnalyzeShrunkData (serial)
        for sr in &mut a_vsd {
            sr.perform(&self.ds);
            let n_e = { let r = sr.pave_block().0.read().unwrap(); r.original_edge };
            let a_e_range = self.ds.shapes[n_e].shape.as_edge()
                .map(|ed| ed.range).unwrap_or([0.0, 0.0]);
            self.analyze_shrunk_data(sr.pave_block(), sr, n_e, a_e_range);
        }
    }

    // OCCT BOPAlgo_PaveFiller::FillShrunkData(handle<PaveBlock>&) (PaveFiller_3.cxx L727-762).
    pub(crate) fn fill_shrunk_data_pb(&mut self, the_pb: &SharedPB) {
        let (n_v1, n_v2) = {
            let pbr = the_pb.0.read().unwrap();
            (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
        };
        if n_v1 >= self.ds.nb_shapes() || n_v2 >= self.ds.nb_shapes() { return; }
        let n_e = {
            let pbr = the_pb.0.read().unwrap();
            if pbr.edge != usize::MAX { pbr.edge } else { pbr.original_edge }
        };
        if n_e >= self.ds.nb_shapes() { return; }
        let (a_t1, a_t2) = { let pbr = the_pb.0.read().unwrap(); pbr.range() };
        let mut sr = ShrunkRange::new(the_pb, n_v1, n_v2, a_t1, a_t2);
        sr.perform(&self.ds);
        let a_e_range = self.ds.shapes[n_e].shape.as_edge()
            .map(|ed| ed.range).unwrap_or([0.0, 0.0]);
        self.analyze_shrunk_data(the_pb, &sr, n_e, a_e_range);
    }

    /// OCCT BOPTools_AlgoTools::IsMicroEdge (BOPTools_AlgoTools.cxx L2075-2112),
    /// evaluated on a section edge shape.
    ///
    /// A section PB in MakeBlocks carries no edge reference (SetEdge is deferred
    /// to PostTreatFF), so rcad cannot drive the micro check through the PB like
    /// fill_shrunk_data_pb does. OCCT checks the section edge shape directly 鈥?
    /// same here.
    pub(crate) fn is_micro_section_edge(&self, a_e: &Shape) -> bool {
        let a_ed = match &*a_e.data {
            topods::TShape::Edge(ed) => ed,
            _ => return true,
        };
        // OCCT L2084: bRet = (Degenerated(aE) || !IsGeometric(aE))
        if a_ed.degenerated || a_ed.curve.is_none() {
            return true;
        }
        let curve = a_ed.curve.clone().unwrap();
        // OCCT L2090-2099: curve range and vertex parameters on the edge.
        let (a_t1, a_t2) = (a_ed.range[0], a_ed.range[1]);
        let (a_p1, mut a_tol_v1) = match &*a_ed.first.data {
            topods::TShape::Vertex(vd) => (vd.point, vd.tolerance),
            _ => return true,
        };
        let (a_p2, mut a_tol_v2) = match &*a_ed.last.data {
            topods::TShape::Vertex(vd) => (vd.point, vd.tolerance),
            _ => return true,
        };
        let a_tol_e = a_ed.tolerance;
        // OCCT IntTools_ShrunkRange::Perform (IntTools_ShrunkRange.cxx L107-191).
        let a_dtol = rcad_kernel::CONFUSION;
        let a_pdtol = rcad_kernel::PCONFUSION;
        if a_t2 - a_t1 < a_pdtol {
            return true;
        }
        if a_tol_v1 < a_tol_e {
            a_tol_v1 = a_tol_e;
        }
        if a_tol_v2 < a_tol_e {
            a_tol_v2 = a_tol_e;
        }
        a_tol_v1 += a_dtol;
        a_tol_v2 += a_dtol;
        let mut a_ts1 = 0.0;
        let mut a_ts2 = 0.0;
        if !find_valid_range_params(
            &curve, a_t1, a_t2, a_tol_e, a_p1, a_tol_v1, a_p2, a_tol_v2,
            &mut a_ts1, &mut a_ts2,
        ) {
            return true;
        }
        if a_ts2 - a_ts1 < a_pdtol {
            return true;
        }
        let mut a_ptol_e = shrunk_range_resolution(&curve, a_ts1, a_ts2, a_tol_e);
        let a_ptol_e_min = (a_t2 - a_t1) * 0.01;
        if a_ptol_e > a_ptol_e_min {
            a_ptol_e = a_ptol_e_min;
        }
        let a_length = shrunk_range_arc_length(&curve, a_ts1, a_ts2, a_ptol_e);
        if a_length < a_dtol {
            return true;
        }
        false
    }

    /// OCCT BOPAlgo_PaveFiller::AnalyzeShrunkData (PaveFiller_3.cxx L766-824).
    // OCCT BOPAlgo_PaveFiller::AnalyzeShrunkData (PaveFiller_3.cxx L766-824).
    fn analyze_shrunk_data(
        &mut self, the_pb: &SharedPB, the_sr: &ShrunkRange,
        n_e: usize, a_e_range: [f64; 2], // edge index + full curve range (from fill_shrunk_data)
    ) {
        // OCCT L770-771: bool bWholeEdge = false; TopoDS_Shape aWarnShape;
        let mut b_whole_edge = false;

        // OCCT L773: if (!theSR.IsDone() || !theSR.IsSplittable())
        if !the_sr.is_done() || !the_sr.is_splittable() {
            // OCCT L776-777: BRep_Tool::Range(edge, aEFirst, aELast); thePB->Range(aPBFirst, aPBLast);
            let (a_e_first, a_e_last) = (a_e_range[0], a_e_range[1]);
            let (a_pb_first, a_pb_last) = { let r = the_pb.0.read().unwrap(); r.range() };
            // OCCT L778: bWholeEdge = aPBFirst <= aEFirst && aPBLast >= aELast;
            b_whole_edge = a_pb_first <= a_e_first && a_pb_last >= a_e_last;

            // OCCT L779-791: warning shape 鈥?rcad skips compound build (no TopoDS)

            // OCCT L793-807: if (!theSR.IsDone())
            if !the_sr.is_done() {
                // OCCT L797-801: AddWarning (TooSmallEdge or BadPositioning)
                if b_whole_edge {
                    self.my_report.add_warning(Alert::TooSmallEdge(n_e));
                } else {
                    self.my_report.add_warning(Alert::BadPositioning(vec![n_e]));
                }
                // OCCT L804-806: thePB->SetShrunkData(aTS1, aTS2, Bnd_Box(), false);
                let (a_ts1, a_ts2) = the_sr.shrunk_range();
                let mut pbr = the_pb.0.write().unwrap();
                pbr.set_shrunk_data(a_ts1, a_ts2, BndBox::new(), false);
                return;
            }
            // OCCT L809-816: AddWarning (NotSplittableEdge or BadPositioning)
            if b_whole_edge {
                self.my_report.add_warning(Alert::NotSplittableEdge(n_e));
            } else {
                self.my_report.add_warning(Alert::BadPositioning(vec![n_e]));
            }
        }

        // OCCT L819-823: set shrunk data with box + fuzzy/2 gap
        let (a_ts1, a_ts2) = the_sr.shrunk_range();
        // OCCT L821: Bnd_Box aBox = theSR.BndBox(); aBox.SetGap(aBox.GetGap() + myFuzzyValue / 2.);
        let mut a_box = the_sr.bnd_box().clone();
        a_box.set_gap(a_box.get_gap() + self.my_fuzzy_value / 2.);
        let mut pbr = the_pb.0.write().unwrap();
        pbr.set_shrunk_data(a_ts1, a_ts2, a_box, the_sr.is_splittable());
    }

    // OCCT BOPAlgo_PaveFiller::ForceInterfVE (PaveFiller_3.cxx L828-910).
    fn force_interf_ve(
        &mut self,
        n_v: usize,
        a_pb: &SharedPB,
        the_m_edges: &mut crate::bop::algo::occt_map::OcctMapInt,
    ) -> bool {
        // OCCT L832-833: int nE, nVx, nVSD, iFlag; double aT, aTolVNew;
        let n_e: usize;
        let mut n_vx: usize;
        let mut n_vsd: usize = usize::MAX;
        let (mut a_t, mut a_tol_v_new): (f64, f64) = (0.0, 0.0);

        // OCCT L835: nE = aPB->OriginalEdge()
        n_e = a_pb.0.read().unwrap().original_edge;
        // OCCT L837: const BOPDS_ShapeInfo& aSIE = myDS->ShapeInfo(nE);
        // rcad: inline self.ds.shapes[n_e]

        // OCCT L838-841: if (aSIE.HasSubShape(nV)) return true;
        if self.ds.shapes[n_e].has_sub_shape(n_v) { return true; }
        // OCCT L843-846: if (myDS->HasInterf(nV, nE)) return true;
        if self.ds.has_interf(n_v, n_e) { return true; }
        // OCCT L848-851: if (myDS->HasInterfShapeSubShapes(nV, nE)) return true;
        if self.ds.has_interf_shape_sub_shapes(n_v, n_e, true) { return true; }
        // OCCT L853-856: if (aPB->Pave1().Index() == nV || aPB->Pave2().Index() == nV) return true;
        {
            let r = a_pb.0.read().unwrap();
            if r.pave1.vertex_idx == n_v || r.pave2.vertex_idx == n_v { return true; }
        }

        // OCCT L858-862: nVx = nV; if (myDS->HasShapeSD(nV, nVSD)) nVx = nVSD;
        n_vx = n_v;
        n_vsd = n_vx;
        if self.ds.has_shape_sd(n_v, &mut n_vx) { n_vsd = n_vx; }

        // OCCT L864-867: iFlag = myContext->ComputeVE(aV, aE, aT, aTolVNew, myFuzzyValue);
        // OCCT: on non-degenerated, geometric edge, projects V onto E and returns 0/1/-4.
        let (i_flag, a_t_val, a_tol_v_new_val) =
            self.my_context.compute_ve(n_vx, n_e, &self.ds, self.my_fuzzy_value);
        if i_flag == -1 || i_flag == -2 || i_flag == -3 { return false; }
        // OCCT L868: if (iFlag == 0 || iFlag == -4)
        if i_flag != 0 && i_flag != -4 { return false; }
        a_t = a_t_val;
        a_tol_v_new = a_tol_v_new_val;

        // OCCT L870: BOPDS_Pave aPave;
        // OCCT L873-874: aVEs.SetIncrement(10);
        // rcad: Vec auto-extends

        // OCCT L876-878: 1 鈥?BOPDS_InterfVE& aVE = aVEs.Appended();
        //                aVE.SetIndices(nV, nE); aVE.SetParameter(aT);
        self.ds.interf_ve.push(InterferenceVE {
            vertex: n_v, edge: n_e, param: a_t, index_new: 0,
        });

        // OCCT L880: 2 鈥?myDS->AddInterf(nV, nE);
        self.ds.add_interf(n_v, n_e);

        // OCCT L883: 3 鈥?nVx = UpdateVertex(nV, aTolVNew);
        let n_vx_new = self.update_vertex(n_v, a_tol_v_new);

        // OCCT L885-888: 4 鈥?if (myDS->IsNewShape(nVx)) aVE.SetIndexNew(nVx);
        if self.ds.is_new_shape(n_vx_new) {
            if let Some(last) = self.ds.interf_ve.last_mut() { last.index_new = n_vx_new; }
        }

        // OCCT L889-892: 5 鈥?aPave.SetIndex(nVx); aPave.SetParameter(aT); aPB->AppendExtPave(aPave);
        a_pb.0.write().unwrap().ext_paves.push(Pave { vertex_idx: n_vx_new, param: a_t });

        // OCCT L894: theMEdges.Add(nE);
        the_m_edges.add(n_e);

        // OCCT L896-906: self-interference warning
        let i_rv = self.ds.rank(n_v);
        if i_rv >= 0 && i_rv == self.ds.rank(n_e) {
            self.my_report.add_warning(Alert::SelfInterferingShape(vec![n_v, n_e]));
        }
        true
    }

    /// OCCT BOPAlgo_PaveFiller::SplitPaveBlocks (PaveFiller_2.cxx L449-560).
    /// Splits PBs with ext paves into sub-PBs. Each ext pave becomes a split point.
    fn split_pave_blocks(&mut self, the_medges: &crate::bop::algo::occt_map::OcctMapInt, the_add_interfs: bool) {
        // OCCT L453: Fence map to avoid unification of the same vertices twice
        let mut a_mpairs: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        // OCCT L455-457: Map to treat the Common Blocks 鈥?NCollection_IndexedDataMap,
        // insertion order (PaveFiller_2.cxx L456-457).
        let mut a_mcb_new_pb: indexmap::IndexMap<usize, Vec<SharedPB>> = indexmap::IndexMap::new();
        // OCCT L460: Map of vertices to init pave blocks for them 鈥?
        // NCollection_Map<int> aMVerticesToInitPB (PaveFiller_2.cxx L460).
        let mut a_m_vertices_to_init_pb: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();
        // OCCT L451-453: aItM.Initialize(theMEdges) 鈥?NCollection_Map<int>
        // iteration order (bucket order, not sorted); the order fixes the
        // MakeSDVertices/SplitEdge vertex and edge creation sequence.
        for n_e in the_medges.iter_keys() {
            if n_e >= self.ds.nb_shapes() { continue; }
            self.ds.init_pave_blocks(n_e);
            let old_pbs: Vec<SharedPB> = {
                let a_lpb = self.ds.change_pave_blocks(n_e);
                if a_lpb.is_empty() { continue; }
                a_lpb.to_vec()
            };
            for pb in &old_pbs {
                // OCCT L443-447: if (!aPB->IsToUpdate()) { aItLPB.Next(); continue; }
                if !pb.0.read().unwrap().is_to_update() { continue; }
                // OCCT L479: myDS->CommonBlock(aPB)
                let a_cb = self.ds.common_block(pb);
                // OCCT L483: aPB->Update(aLPBN) 鈥?theFlag defaults to TRUE
                // (BOPDS_PaveBlock.hxx L136-137), so the endpoint paves are
                // included and an edge with one intersection ext pave splits
                // into two sub-pave-blocks.
                let mut a_lpbn: Vec<SharedPB> = Vec::new();
                {
                    let mut pbw = pb.0.write().unwrap();
                    pbw.update(&mut a_lpbn, true);
                }
                let mut valid_new_pbs: Vec<SharedPB> = Vec::new();
                for a_pbn in &a_lpbn {
                    // OCCT L461: myDS->UpdatePaveBlockWithSDVertices(aPBN)
                    self.ds.update_pave_block_with_sd_vertices(a_pbn);
                    // OCCT L462: FillShrunkData(aPBN)
                    let (sv1, sv2, st1, st2) = { let r = a_pbn.0.read().unwrap(); let (v1, v2) = r.indices(); let (t1, t2) = r.range(); (v1, v2, t1, t2) };
                    let mut sr = ShrunkRange::new(a_pbn, sv1, sv2, st1, st2);
                    sr.perform(&self.ds);
                    let a_e_range = self.ds.shapes[n_e].shape.as_edge()
                        .map(|ed| ed.range).unwrap_or([0.0, 0.0]);
                    self.analyze_shrunk_data(a_pbn, &sr, n_e, a_e_range);
                    // OCCT L464: bool bHasValidRange = aPBN->HasShrunkData()
                    let b_has_valid_range = { let r = a_pbn.0.read().unwrap(); r.has_shrunk_data() };
                    // OCCT L468: bool bCheckDist = (bHasValidRange && !aPBN->IsSplittable())
                    let b_check_dist = b_has_valid_range && !a_pbn.0.read().unwrap().is_splittable();
                    // OCCT L469-508: if (!bHasValidRange || bCheckDist) { ... MakeSDVertices; continue; }
                    if !b_has_valid_range || b_check_dist {
                        let (n_v1, n_v2) = { let r = a_pbn.0.read().unwrap(); r.indices() };
                        if n_v1 == n_v2 { continue; }
                        if b_check_dist {
                            let a_dist_vv = (self.ds.vertex_point_by_idx(n_v1) - self.ds.vertex_point_by_idx(n_v2)).length();
                            if a_dist_vv <= self.my_fuzzy_value.max(rcad_kernel::CONFUSION) {
                                // vertices interfering -> no valid range
                                // fall through to MakeSDVertices
                            } else {
                                if !b_has_valid_range { continue; }
                            }
                        }
                        // OCCT L491-506: MakeSDVertices
                        if !b_has_valid_range {
                            let a_pair = if n_v1 < n_v2 { (n_v1, n_v2) } else { (n_v2, n_v1) };
                            if a_mpairs.insert(a_pair) {
                                let a_lv = vec![n_v1, n_v2];
                                self.make_sd_vertices_vv(&a_lv, the_add_interfs);
                                a_m_vertices_to_init_pb.add(n_v1);
                                a_m_vertices_to_init_pb.add(n_v2);
                            }
                        }
                        continue;
                    }
                    // OCCT L511: aLPB.Append(aPBN)
                    valid_new_pbs.push(a_pbn.clone());
                    // OCCT L513-523: CommonBlock handling
                    if let Some(cb) = a_cb {
                        a_mcb_new_pb.entry(cb).or_default().push(a_pbn.clone());
                    }
                }
                // OCCT L525-526: aLPB.Remove(aItLPB) — remove ONLY the current
                // pave block and append the new sub-PBs. rcad previously
                // cleared the whole list, which destroyed the edge's OTHER
                // (unprocessed) pave blocks (e.g. a neighbouring split piece
                // with no ext paves of its own).
                if !valid_new_pbs.is_empty() {
                    let a_lpb = self.ds.change_pave_blocks(n_e);
                    let cur_ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
                    a_lpb.retain(|p| std::sync::Arc::as_ptr(&p.0) as u64 != cur_ptr);
                    for new_pb in valid_new_pbs {
                        a_lpb.push(new_pb);
                    }
                }
            }
        }
        // OCCT L530-560: Make Common Blocks
        // aMCBNewPB is NCollection_IndexedDataMap 鈥?iterate in insertion order.
        for (_cb_key, pbs) in &a_mcb_new_pb {
            if pbs.len() < 2 { continue; }
            // OCCT L566-571: aMInds is NCollection_IndexedDataMap<BOPDS_Pair,
            // List<PB>> 鈥?insertion order of the aLPBN traversal.
            let mut a_minds: indexmap::IndexMap<(usize, usize), Vec<SharedPB>> = indexmap::IndexMap::new();
            for pb in pbs {
                let (v1, v2) = { let r = pb.0.read().unwrap(); r.indices() };
                let key = if v1 < v2 { (v1, v2) } else { (v2, v1) };
                a_minds.entry(key).or_default().push(pb.clone());
            }
            for (_pair, group) in &a_minds {
                if group.len() < 2 { continue; }
                self.ds.add_common_block(group);
            }
        }
        // Init PBs for new SD vertices 鈥?OCCT L560: aMVerticesToInitPB is a
        // NCollection_Map<int>, iterated in bucket order (PaveFiller_2.cxx L560).
        for v in a_m_vertices_to_init_pb.iter_keys() {
            self.ds.init_pave_blocks_for_vertex(v);
        }
    }

    // ====================================================================
    // GetPBBox 鈥?OCCT BOPAlgo_PaveFiller_3.cxx L914-955
    // ====================================================================

    /// Get bounding box of a PaveBlock's edge segment.
    /// OCCT BOPAlgo_PaveFiller::GetPBBox (PaveFiller_3.cxx L914-955).
    fn get_pb_box(
        &self,
        _the_e: usize,
        the_pb: &SharedPB,
        the_pb_box: &mut std::collections::HashMap<u64, rcad_kernel::math::bnd::BndBox>,
        the_first: &mut f64,
        the_last: &mut f64,
        the_s_first: &mut f64,
        the_s_last: &mut f64,
        the_box: &mut rcad_kernel::math::bnd::BndBox,
    ) -> bool {
        let pbr = the_pb.0.read().unwrap();
        (*the_first, *the_last) = pbr.range();
        // OCCT L923-929: bool bValid = theLast - theFirst > Precision::PConfusion();
        if *the_last - *the_first <= rcad_kernel::PCONFUSION {
            if std::env::var("RCAD_EE_DEBUG").is_ok() {
                eprintln!("[EE-DBG] getpb false edge={} pb range=[{:.12},{:.12}] len={:.3e} shrunk={}", pbr.original_edge, *the_first, *the_last, *the_last - *the_first, pbr.has_shrunk_data());
            }
            return false;
        }
        // OCCT L932-937: check shrunk data
        if pbr.has_shrunk_data() {
            *the_s_first = pbr.ts1;
            *the_s_last = pbr.ts2;
            // OCCT L937: aBox = theSR.BndBox() (with the fuzzy/2 gap added in
            // SetShrunkData), carried on the pave block.
            *the_box = pbr.shrunk_bnd_box.clone();
            return true;
        }
        *the_s_first = *the_first;
        *the_s_last = *the_last;
        // OCCT L942-952: check map, then build bounding box
        let pb_ptr = std::sync::Arc::as_ptr(&the_pb.0) as u64;
        if let Some(bb) = the_pb_box.get(&pb_ptr) {
            *the_box = bb.clone();
        } else {
            let curve = self.ds.edge_curve(pbr.original_edge);
            let bb = if let Some(c) = curve {
                // OCCT L949-951: aTol = BRep_Tool::Tolerance(theE) + Precision::Confusion();
                // BndLib_Add3dCurve::Add(aBAC, theSFirst, theSLast, aTol, theBox);
                let a_tol = self.ds.edge_tolerance(pbr.original_edge) + rcad_kernel::CONFUSION;
                shrunk_range_bnd_box(&c, *the_s_first, *the_s_last, a_tol)
            } else {
                return false;
            };
            *the_box = bb.clone();
            the_pb_box.insert(pb_ptr, bb);
        }
        true
    }

    // ====================================================================
    // UpdateVertex 鈥?OCCT BOPAlgo_PaveFiller::UpdateVertex (PaveFiller_10.cxx L60-85)
    // ====================================================================

    /// OCCT BOPAlgo_PaveFiller::UpdateVertex (PaveFiller_10.cxx L105-125).
    /// Returns the vertex index after SD resolution.
    pub(crate) fn update_vertex(&mut self, n_v: usize, tol_new: f64) -> usize {
        let mut n_vnew = n_v;
        self.ds.has_shape_sd(n_v, &mut n_vnew);
        // OCCT L112: if (IsNewShape(nVNew) || HasShapeSD(nV, nVNew) || !myNonDestructive)
        if self.ds.is_new_shape(n_vnew) || n_vnew != n_v || !self.my_non_destructive {
            // Path 1: update tolerance and box
            let tol_old = self.ds.vertex_tolerance_by_idx(n_vnew);
            if tol_new > tol_old {
                // In-place (OCCT BRep_Builder UpdateVertex): keeps the vertex
                // identity (DS entry + face-wire references) intact 鈥?see
                // set_vertex_tolerance. Arc::make_mut would clone a shared vertex.
                let pt = self.ds.vertex_point_by_idx(n_vnew);
                self.ds.mutate_shape_data(n_vnew, |ts| {
                    if let rcad_kernel::topods::TShape::Vertex(vd) = ts {
                        vd.tolerance = tol_new;
                    }
                });
                let si = self.ds.change_shape_info(n_vnew);
                si.bbox = BndBox::from_point(pt);
                si.bbox.set_gap(tol_new + rcad_kernel::CONFUSION);
                self.ds.remap_shape_idx(n_vnew);
                self.my_increased_ss.insert(n_v);
            }
            return n_vnew;
        }
        // OCCT L121-150: Path 2 鈥?nV is an old (source) vertex: create a new
        // vertex with the increased tolerance and register the SD relation.
        let a_tol_v = self.ds.vertex_tolerance_by_idx(n_v);
        let a_pv = self.ds.vertex_point_by_idx(n_v);
        // OCCT L127-128: aBB.MakeVertex(aVNew, aPV, max(aTolV, aTolNew)).
        let n_vnew = self.append_vertex(a_pv, a_tol_v.max(tol_new));
        // OCCT L137-143: bounding box of the new vertex with gap.
        let mut a_box = BndBox::from_point(a_pv);
        a_box.set_gap(a_tol_v.max(tol_new) + rcad_kernel::CONFUSION);
        self.ds.change_shape_info(n_vnew).bbox = a_box;
        // OCCT L147: myDS->AddShapeSD(nV, nVNew).
        self.ds.add_shape_sd(n_v, n_vnew);
        // OCCT L149: myVertsToAvoidExtension.Add(nVNew).
        self.my_verts_to_avoid_extension.insert(n_vnew);
        n_vnew
    }

    // ====================================================================
    // UpdateVerticesOfCB 鈥?OCCT BOPAlgo_PaveFiller_3.cxx L959-993
    // ====================================================================

    /// Update vertex tolerances from CommonBlock tolerances.
    /// OCCT BOPAlgo_PaveFiller::UpdateVerticesOfCB (PaveFiller_3.cxx L959-993).
    // OCCT BOPAlgo_PaveFiller::UpdateVerticesOfCB (PaveFiller_3.cxx L959-993).
    fn update_vertices_of_cb(&mut self) {
        let mut a_mpb_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // OCCT L974-976: myDS->ChangePaveBlocksPool() is a DynamicArray 鈥?
        // iterate the pool in ascending edge-key order (a HashMap values()
        // order would randomize the UpdateVertex sequence and hence the
        // indices of the newly created vertices).
        let mut pb_keys: Vec<usize> = self.ds.pave_blocks_pool.keys().copied().collect();
        pb_keys.sort_unstable();
        for k in pb_keys {
            let a_lpb = match self.ds.pave_blocks_pool.get(&k) {
                Some(v) => v.clone(),
                None => continue,
            };
            for a_pb in &a_lpb {
                let a_cb_idx = self.ds.common_block(a_pb);
                let a_cb_idx = match a_cb_idx { Some(idx) => idx, None => continue, };
                let a_cb = &self.ds.common_blocks[a_cb_idx];
                // OCCT L979-980: const handle<PaveBlock>& aPBR = aCB->PaveBlock1();
                // rcad: use a_pb's pointer for fence (same semantic: each CB processed once)
                let a_pb_key = Arc::as_ptr(&a_pb.0) as u64;
                if !a_mpb_fence.insert(a_pb_key) { continue; }
                // OCCT L985-990: aTolCB = aCB->Tolerance(); UpdateVertex(Pave1, Tol); UpdateVertex(Pave2, Tol);
                let a_tol_cb = a_cb.tolerance();
                if a_tol_cb > 0. {
                    if let Some(pb1) = a_cb.pave_block1() {
                        let (nv1, nv2) = pb1.0.read().unwrap().indices();
                        self.update_vertex(nv1, a_tol_cb);
                        self.update_vertex(nv2, a_tol_cb);
                    }
                }
            }
        }
    }

    // ====================================================================
    // RepeatIntersection 鈥?OCCT BOPAlgo_PaveFiller.cxx L383-448
    // ====================================================================
    /// Re-run VV/VE/VF intersections for vertices whose tolerance was increased.
    /// OCCT BOPAlgo_PaveFiller::RepeatIntersection (PaveFiller.cxx L383-448).
    fn repeat_intersection(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // L385-386: NCollection_Map<int> anExtraInterfMap;
        let mut an_extra = HashSet::new();
        // L387: const int aNbS = myDS->NbSourceShapes();
        let a_nb_s = self.ds.nb_source_shapes();
        // L388: Message_ProgressScope aPS(theRange, "Repeat intersection", 3);
        // L389-414: for (int i = 0; i < aNbS; ++i)
        for i in 0..a_nb_s {
            if the_range.user_break() { return; }
            // L391-395: if ShapeType != VERTEX, continue
            if self.ds.shapes[i].shape_type != ShapeType::Vertex {
                continue;
            }
            // L397-401: if (myIncreasedSS.Contains(i)) { anExtraInterfMap.Add(i); continue; }
            if self.my_increased_ss.contains(&i) {
                an_extra.insert(i);
                continue;
            }
            // L404-408: int nVSD; if (!myDS->HasShapeSD(i, nVSD)) { continue; }
            let mut n_vsd = usize::MAX;
            if !self.ds.has_shape_sd(i, &mut n_vsd) {
                continue;
            }
            // L410-413: if (myIncreasedSS.Contains(nVSD)) { anExtraInterfMap.Add(i); }
            if self.my_increased_ss.contains(&n_vsd) {
                an_extra.insert(i);
            }
        }
        // L416-419: if (anExtraInterfMap.IsEmpty()) return;
        if an_extra.is_empty() {
            return;
        }

        // L422: myIterator->IntersectExt(anExtraInterfMap);
        // OCCT expands the pair lists to include the extra vertices.
        if let Some(it) = &mut self.my_iterator {
            it.intersect_ext(&self.ds, &an_extra);
        }

        // L426-430: PerformVV(aPS.Next());
        // L431-445: PerformVE, PerformVF also use aPS.Next()
        let a_ps = the_range.sub_scope("Repeat intersection", 3);
        self.perform_vv(&a_ps.sub_scope("VV", 1));
        if self.has_errors() { return; }
        // L431: UpdatePaveBlocksWithSDVertices();
        self.update_pave_blocks_with_sd_vertices();

        // L433-438: PerformVE(aPS.Next());
        self.perform_ve(&a_ps.sub_scope("VE", 1));
        if self.has_errors() { return; }
        // L438: UpdatePaveBlocksWithSDVertices();
        self.update_pave_blocks_with_sd_vertices();

        // L440-444: PerformVF(aPS.Next());
        self.perform_vf(&a_ps.sub_scope("VF", 1));
        if self.has_errors() { return; }

        // L446-447: UpdatePaveBlocksWithSDVertices(); UpdateInterfsWithSDVertices();
        self.update_pave_blocks_with_sd_vertices();
        self.update_interfs_with_sd_vertices();
    }

    // ====================================================================
    // ForceInterfEE 鈥?OCCT BOPAlgo_PaveFiller_3.cxx L997-1333
    // ====================================================================
    /// Force additional EE intersection for common blocks.
    /// OCCT BOPAlgo_PaveFiller::ForceInterfEE (PaveFiller_3.cxx L997-1333).
    fn force_interf_ee(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // L999-1003: comment 鈥?now that vertices are increased/unified,
        // find additional common blocks among edge pairs with same bounding vertices.

        // L1005-1023: Initialize pave blocks for all vertices that participated
        // in intersections.
        // OCCT: for (int i = 0; i < aNbS; ++i)
        //   if VERTEX && HasInterf(i) -> InitPaveBlocksForVertex(i)
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            if the_range.user_break() { return; }
            if self.ds.shapes[i].shape_type != ShapeType::Vertex {
                continue;
            }
            // L1014: if (myDS->HasInterf(i))
            // rcad: check interf_tb for any pair involving i
            let has_interf = self.ds.interf_tb.iter().any(|&(a, b)| a == i || b == i);
            if has_interf {
                self.ds.init_pave_blocks_for_vertex(i);
            }
        }

        // L1024-1080: Fill the connection map from bounding vertices to PBs
        // L1026-1028: NCollection_IndexedDataMap<BOPDS_Pair, List<PaveBlock>> aPBMap
        // rcad: IndexMap keyed by (v_min, v_max), value = Vec<SharedPB> 鈥?the
        // insertion order (source-edge order) determines the EE-pair processing
        // order below (a HashMap would randomize it and hence the resulting
        // common blocks and split-edge indices).
        let mut a_pb_map: indexmap::IndexMap<(usize, usize), Vec<SharedPB>> =
            indexmap::IndexMap::new();
        // L1030: Fence map of pave blocks
        // rcad: HashSet of PB pointer
        let mut a_mpb_fence: std::collections::HashSet<u64> =
            std::collections::HashSet::new();

        for i in 0..a_nb_s {
            if the_range.user_break() { return; }
            // L1034-1038: only edges
            if self.ds.shapes[i].shape_type != ShapeType::Edge {
                continue;
            }
            // L1041-1044: edge must have PBs (HasReference equivalent)
            // rcad: check if the shape has reference (points to pave_blocks_pool)
            if self.ds.shapes[i].reference < 0 {
                continue;
            }
            // L1047-1051: skip degenerated edges (HasFlag)
            if self.ds.shapes[i].has_flag() {
                continue;
            }

            // L1056-1079: iterate PBs of this edge
            let ei = i;
            let a_lpb = self.ds.edge_pave_blocks(ei);
            for a_pb in a_lpb {
                // L1060-1061: RealPaveBlock 鈥?resolve through CommonBlock
                let a_pbr = self.ds.real_pave_block(a_pb);
                // L1062-1065: fence map 鈥?skip if already processed
                let ptr = std::sync::Arc::as_ptr(&a_pbr.0) as u64;
                if !a_mpb_fence.insert(ptr) {
                    continue;
                }

                // L1068-1069: get vertex indices
                let (n_v1, n_v2) = {
                    let pbr = a_pbr.0.read().unwrap();
                    (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
                };

                // L1072-1078: add PB to map keyed by vertex pair
                // OCCT: BOPDS_Pair aPair(nV1, nV2);
                let a_pair = if n_v1 <= n_v2 { (n_v1, n_v2) } else { (n_v2, n_v1) };
                a_pb_map.entry(a_pair).or_default().push(a_pbr.clone());
            }
        }

        // L1082-1086: empty map check
        if a_pb_map.is_empty() {
            return;
        }

        // L1088: const bool bSICheckMode = (myArguments.Extent() == 1);
        let b_si_check_mode = self.my_arguments.len() == 1;

        // L1090-1225: Prepare pairs for intersection
        // L1091: BOPAlgo_VectorOfEdgeEdge aVEdgeEdge;
        struct EEPair {
            a_pb1: SharedPB,
            a_pb2: SharedPB,
            n_e1: usize,
            n_e2: usize,
            pb1_range: (f64, f64),
            pb2_range: (f64, f64),
            fuzzy_value: f64,
        }
        let mut edge_edge_pairs: Vec<EEPair> = Vec::new();

        for (&a_pair, a_lpb) in &a_pb_map {
            let (n_v1, n_v2) = a_pair;
            // L1100-1102: if less than 2 PBs, skip
            if a_lpb.len() < 2 {
                continue;
            }

            // L1105-1110: get vertex shapes for tolerance computation
            // OCCT: const TopoDS_Vertex& aV1 = TopoDS::Vertex(myDS->Shape(nV1));
            //        const TopoDS_Vertex& aV2 = TopoDS::Vertex(myDS->Shape(nV2));
            // rcad: get vertex tolerances from DS
            let tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
            let tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);

            // L1116-1118: aTolAdd = bSICheckMode ? myFuzzyValue : 2*max(BRep_Tool::Tolerance(aV1), aV2)
            let a_tol_add = if b_si_check_mode {
                self.my_fuzzy_value
            } else {
                2.0 * tol_v1.max(tol_v2)
            };

            // L1121-1224: iterate all unique pairs from the list
            for p1_idx in 0..a_lpb.len() {
                let a_pb1 = a_lpb[p1_idx].clone();
                // L1125-1126: get CommonBlock status
                let cb1_idx = self.ds.common_block(&a_pb1);
                let (n_e1, i_r1) = {
                    let pbr = a_pb1.0.read().unwrap();
                    (pbr.original_edge, self.ds.rank(pbr.original_edge))
                };
                // L1127-1130: edge and its range
                let (a_t11, a_t12) = {
                    let pbr = a_pb1.0.read().unwrap();
                    pbr.range()
                };
                // OCCT L1131-1139: BRepAdaptor_Curve aBAC1(aE1); aBAC1.D1(midpoint, aPm, aVTgt1);
                //   if (aVTgt1.SquareMagnitude() < gp::Resolution()) continue;
                let c1 = self.ds.edge_curve(n_e1);
                let v_tgt1 = c1.as_ref().map(|c| {
                    let mid_t = (a_t11 + a_t12) * 0.5;
                    let dt = 1e-7;
                    let p_mid = c.point_at(mid_t);
                    let p_dt = c.point_at(mid_t + dt);
                    p_dt - p_mid
                });
                let (v_tgt1, a_pm) = match v_tgt1 {
                    // OCCT L1135: if (aVTgt1.SquareMagnitude() < gp::Resolution()) continue;
                    // gp::Resolution() == 1e-7 in OCCT.
                    Some(v) if v.length_squared() > 1e-7 => {
                        let mid_t = (a_t11 + a_t12) * 0.5;
                        let a_pm = c1.as_ref().unwrap().point_at(mid_t);
                        (v.normalize(), a_pm)
                    },
                    _ => continue,
                };

                // L1141 onwards: iterate second PB for each pair
                for p2_idx in (p1_idx + 1)..a_lpb.len() {
                    let a_pb2 = a_lpb[p2_idx].clone();
                    let cb2_idx = self.ds.common_block(&a_pb2);
                    let (n_e2, i_r2) = {
                        let pbr = a_pb2.0.read().unwrap();
                        (pbr.original_edge, self.ds.rank(pbr.original_edge))
                    };

                    // L1149-1160: skip edges from same argument unless vertices are new
                    // OCCT: if (iR1 == iR2) {
                    //   if ((!IsNewShape(nV1) && Rank(nV1) == iR1) ||
                    //       (!IsNewShape(nV2) && Rank(nV2) == iR2)) continue; }
                    if i_r1 == i_r2 && i_r1 >= 0 {
                        let v1_original = !self.ds.is_new_shape(n_v1) && self.ds.rank(n_v1) == i_r1;
                        let v2_original = !self.ds.is_new_shape(n_v2) && self.ds.rank(n_v2) == i_r2;
                        if v1_original || v2_original {
                            continue;
                        }
                    }

                    // L1162-1168: skip if PBs already form the SAME common block
                    // OCCT: if (!aCB1.IsNull() && !aCB2.IsNull()) { if (aCB1 == aCB2) continue; }
                    if let (Some(cb1), Some(cb2)) = (cb1_idx, cb2_idx) {
                        if cb1 == cb2 {
                            continue;
                        }
                    }

                    // L1175-1204: check angle between edges at midpoint
                    // bUseAddTol = true initially; if angle > 25deg, set to false
                    let (a_t21, a_t22) = {
                        let pbr = a_pb2.0.read().unwrap();
                        pbr.range()
                    };
                    let b_use_add_tol = {
                        let c2 = self.ds.edge_curve(n_e2);
                        let mut use_tol = true;
                        if let Some(c) = c2 {
                            let mid_t2 = (a_t21 + a_t22) * 0.5;
                            let dt = 1e-7;
                            let p_mid2 = c.point_at(mid_t2);
                            let p_dt2 = c.point_at(mid_t2 + dt);
                            let v_tgt2 = p_dt2 - p_mid2;
                            // OCCT L1193: if (aVTgt2.SquareMagnitude() < gp::Resolution()) continue;
                        if v_tgt2.length_squared() > 1e-7 {
                                let a_cos = v_tgt2.normalize().dot(v_tgt1).abs();
                                // OCCT L1199-1203: if (std::abs(aCos) < 0.9063) bUseAddTol = false;
                                if a_cos < 0.9063 {
                                    use_tol = false;
                                }
                            }
                        }
                        use_tol
                    };

                    // L1208-1222: add pair with appropriate fuzzy value
                    // OCCT: if (bUseAddTol) anEdgeEdge.SetFuzzyValue(myFuzzyValue + aTolAdd)
                    //        else anEdgeEdge.SetFuzzyValue(myFuzzyValue)
                    let fuzzy_val = if b_use_add_tol {
                        self.my_fuzzy_value + a_tol_add
                    } else {
                        self.my_fuzzy_value
                    };
                    edge_edge_pairs.push(EEPair {
                        a_pb1: a_pb1.clone(),
                        a_pb2: a_pb2.clone(),
                        n_e1, n_e2,
                        pb1_range: (a_t11, a_t12),
                        pb2_range: (a_t21, a_t22),
                        fuzzy_value: fuzzy_val,
                    });
                }
            }
        }

        // L1227-1231: if no pairs, return
        if edge_edge_pairs.is_empty() {
            return;
        }

        // L1248-1252: Perform intersection (OCCT: BOPTools_Parallel::Perform)
        // rcad: serial intersection of each pair
        // L1253: NCollection_DynamicArray<BOPDS_InterfEE>& aEEs = myDS->InterfEE();
        // L1312-1329: Collect PB pairs for CommonBlock creation
        // rcad: BOPAlgo_Tools::PerformCommonBlocks creates a CB per unique PB pair.
        // We collect (pb1, pb2) pairs and create CBs in a second pass.
        let mut cb_pairs: Vec<(SharedPB, SharedPB)> = Vec::new();

        for pair in &edge_edge_pairs {
            // L1264-1290: intersect edges
            let pb1 = pair.a_pb1.clone();
            let pb2 = pair.a_pb2.clone();

            let c1 = self.ds.edge_curve(pair.n_e1);
            let c2 = self.ds.edge_curve(pair.n_e2);
            let (c1, c2) = match (c1, c2) {
                (Some(c1), Some(c2)) => (c1, c2),
                _ => continue,
            };

            let mut ee = crate::bop::int_tools::edge_edge::EdgeEdgeIntersector::new();
            ee.set_edges(pair.n_e1, [pair.pb1_range.0, pair.pb1_range.1], pair.n_e2, [pair.pb2_range.0, pair.pb2_range.1], &self.ds);
            // OCCT L1216-1222: SetFuzzyValue with the pair-specific tolerance
            ee.set_fuzzy_value(pair.fuzzy_value);
            ee.perform();

            if !ee.is_done() {
                // L1272-1278: warn about failed intersection
                self.my_report.add_warning(
                    Alert::IntersectionFailed(pair.n_e1, pair.n_e2));
                continue;
            }

            let a_cparts = ee.common_parts();
            // L1282-1285: only accept 1 common part of type EDGE
            if a_cparts.len() != 1 {
                continue;
            }
            let cp = &a_cparts[0];
            // L1288: if (aCP.Type() != TopAbs_EDGE) continue;
            // rcad: the old intersector does not set is_edge for full coincidences,
            // but the same-part length check serves as the EDGE-type proxy.
            if cp.range1[0] >= cp.range1[1] {
                continue;
            }

            // L1293-1310: add interference
            let new_ee = InterferenceEE {
                e1: pair.n_e1,
                e2: pair.n_e2,
                point: cp.bounding_point1,
                param1: cp.range1[0],
                param2: cp.ranges2[0][0],
                new_vertex: usize::MAX,
                range1: cp.range1,
                range2: cp.ranges2[0],
            };
            self.ds.interf_ee.push(new_ee);
            self.ds.add_interf(pair.n_e1, pair.n_e2);

            // L1297-1305: if same rank, add AcquiredSelfIntersection warning
            let r1 = self.ds.rank(pair.n_e1);
            let r2 = self.ds.rank(pair.n_e2);
            if r1 >= 0 && r1 == r2 {
                self.my_report.add_warning(
                    Alert::AcquiredSelfIntersection(vec![pair.n_e1, pair.n_e2]));
            }

            // L1312-1329: fill map for common block creation
            // OCCT: BOPAlgo_Tools::FillMap(aPB[0], aPB[1], aMPBLPB, anAlloc)
            cb_pairs.push((pb1, pb2));
        }

        // L1312-1332: BOPAlgo_Tools::PerformCommonBlocks(aMPBLPB, anAlloc, myDS)
        // OCCT builds a connection graph via FillMap and expands through existing
        // CommonBlocks, then groups all connected PBs into merged CommonBlocks.
        // rcad: collect all PBs, expanding from existing CommonBlocks, then
        // create a merged CommonBlock for each connected group.
        {
            // OCCT: NCollection_IndexedDataMap<PB, List<PB>> aMPBLPB 鈥?insertion
            // order; the BOPAlgo_Tools::MakeBlocks traversal is FIFO.
            type PBKey = u64;
            let mut adj: indexmap::IndexMap<PBKey, indexmap::IndexSet<PBKey>> =
                indexmap::IndexMap::new();

            // Helper to get or create adjacency entry
            let mut add_edge = |a: PBKey, b: PBKey| {
                adj.entry(a).or_default().insert(b);
                adj.entry(b).or_default().insert(a);
            };

            // OCCT L1312-1327: expand through existing CommonBlocks
            for (pb1, pb2) in &cb_pairs {
                let pbs = [pb1.clone(), pb2.clone()];
                for pb in &pbs {
                    let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
                    // If this PB is in a CommonBlock, connect it to ALL PBs in that CB
                    let cb_idx = pb.0.read().unwrap().common_block_idx;
                    if let Some(cb_idx) = cb_idx {
                        if cb_idx < self.ds.common_blocks.len() {
                            for (pool_pb, _face_idx) in self.ds.common_blocks[cb_idx].pave_blocks() {
                                let other_ptr = std::sync::Arc::as_ptr(&pool_pb.0) as u64;
                                if other_ptr != ptr {
                                    add_edge(ptr, other_ptr);
                                }
                            }
                        }
                    }
                }
                // OCCT L1329: FillMap(aPB[0], aPB[1])
                let ptr1 = std::sync::Arc::as_ptr(&pb1.0) as u64;
                let ptr2 = std::sync::Arc::as_ptr(&pb2.0) as u64;
                if ptr1 != ptr2 {
                    add_edge(ptr1, ptr2);
                }
            }

            // OCCT L122-123: MakeBlocks 鈥?group connected PBs via graph traversal
            let mut visited: std::collections::HashSet<PBKey> = std::collections::HashSet::new();
            let mut groups: Vec<Vec<SharedPB>> = Vec::new();

            // Build PB lookup: ptr 鈫?SharedPB
            let mut ptr_to_pb: std::collections::HashMap<PBKey, SharedPB> = std::collections::HashMap::new();
            for (pb1, pb2) in &cb_pairs {
                for pb in [pb1.clone(), pb2.clone()] {
                    let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
                    ptr_to_pb.entry(ptr).or_insert(pb);
                }
            }
            // Also add PBs from existing CommonBlocks
            for (pb1, pb2) in &cb_pairs {
                for pb in [pb1, pb2] {
                    let cb_idx = pb.0.read().unwrap().common_block_idx;
                    if let Some(cb_idx) = cb_idx {
                        if cb_idx < self.ds.common_blocks.len() {
                            for (pool_pb, _) in self.ds.common_blocks[cb_idx].pave_blocks() {
                                let ptr = std::sync::Arc::as_ptr(&pool_pb.0) as u64;
                                ptr_to_pb.entry(ptr).or_insert(pool_pb.clone());
                            }
                        }
                    }
                }
            }

            // FIFO connected-group traversal (OCCT BOPAlgo_Tools::MakeBlocks
            // grows the chain list from the start; appended elements are
            // visited later 鈥?a stack would visit in a different order).
            for &start in adj.keys() {
                if visited.contains(&start) { continue; }
                let mut group: Vec<SharedPB> = Vec::new();
                let mut queue: std::collections::VecDeque<PBKey> = std::collections::VecDeque::new();
                queue.push_back(start);
                while let Some(node) = queue.pop_front() {
                    if !visited.insert(node) { continue; }
                    if let Some(pb) = ptr_to_pb.get(&node) {
                        group.push(pb.clone());
                    }
                    if let Some(neighbors) = adj.get(&node) {
                        for &n in neighbors {
                            if !visited.contains(&n) {
                                queue.push_back(n);
                            }
                        }
                    }
                }
                if group.len() >= 2 {
                    groups.push(group);
                }
            }

            // OCCT L130-185: create CommonBlock for each group
            for group in &groups {
                self.ds.add_common_block(group);
            }
        }
    }

    // ====================================================================
    // ForceInterfEF 鈥?OCCT BOPAlgo_PaveFiller_5.cxx L772-1199
    // ====================================================================
    /// Force additional EF intersection for common blocks.
    /// OCCT BOPAlgo_PaveFiller::ForceInterfEF (PaveFiller_5.cxx L772-827).
    fn force_interf_ef(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // L774-775: Message_ProgressScope aPS(theRange, nullptr, 1);
        // L776-778: if (!myIsPrimary) return;
        if !self.my_is_primary {
            return;
        }

        // L787-822: Collect all pave blocks into an IndexedMap
        // OCCT: NCollection_IndexedMap<handle<BOPDS_PaveBlock>> aMPB 鈥?the
        // insertion order determines the EF-pair processing order, hence the
        // InterfEF array order and the FaceInfo PaveBlocksIn order (a HashSet
        // would randomize it and with it the Builder aLE edge order).
        let mut a_mpb: indexmap::IndexSet<(usize, usize)> =
            indexmap::IndexSet::new();
        let a_nb_s = self.ds.nb_source_shapes();
        for n_e in 0..a_nb_s {
            // L791-795: only edges
            if self.ds.shapes[n_e].shape_type != ShapeType::Edge {
                continue;
            }
            // L798-801: edge must have PBs
            if self.ds.shapes[n_e].reference < 0 {
                continue;
            }
            // L804-807: skip degenerated edges
            if self.ds.shapes[n_e].has_flag() {
                continue;
            }

            // L814-821: iterate PBs
            self.ds.init_pave_blocks(n_e);
            let a_lpb = self.ds.change_pave_blocks(n_e);
            for local_i in 0..a_lpb.len() {
                // OCCT L819: aMPB.Add(aPBR) where aPBR = myDS->RealPaveBlock(aPB)
                // rcad: no RealPaveBlock indirection, use (n_e, local_i) as key
                a_mpb.insert((n_e, local_i));
            }
        }

        // L826: ForceInterfEF(aMPB, aPS.Next(), true);
        self.force_interf_ef_work(&a_mpb, true);
    }

    /// OCCT BOPAlgo_PaveFiller::ForceInterfEF (overload, PaveFiller_5.cxx L831-1199).
    /// Worker function 鈥?processes collected pave blocks against all faces.
    pub(crate) fn force_interf_ef_work(
        &mut self,
        the_mpb: &indexmap::IndexSet<(usize, usize)>,
        the_add_interf: bool,
    ) {
        // L838-841: if (theMPB.IsEmpty()) return;
        if the_mpb.is_empty() {
            return;
        }

        // L843-871: BOPTools_BoxTree aBBTree 鈥?build BVH tree of PBs.
        // rcad: iterates all PB/face pairs with direct BndBox overlap checks.

        // L876: const bool bSICheckMode = (myArguments.Extent() == 1);
        let b_si_check_mode = self.my_arguments.len() == 1;

        // L882-1107: For each face, find overlapping PBs and check
        let a_nb_s = self.ds.nb_source_shapes();
        // (n_e, n_f, pb_pool_idx, a_tol_add, PB pave range)
        let mut ef_pairs: Vec<(usize, usize, usize, f64, (f64, f64))> = Vec::new();

        for n_f in 0..a_nb_s {
            if self.ds.shapes[n_f].shape_type != ShapeType::Face {
                continue;
            }
            if self.ds.shapes[n_f].reference < 0 {
                continue;
            }

            // L912-924: Collect vertices of the face from its FaceInfo
            let a_fi = self.ds.face_info(n_f);
            let face_pb_on = a_fi.pave_blocks_on.clone();
            let face_pb_in = a_fi.pave_blocks_in.clone();
            let face_pb_sc = a_fi.pave_blocks_sc.clone();
            let mut a_mvf: std::collections::HashSet<usize> = std::collections::HashSet::new();
            // OCCT L916-924: aMVF from VerticesOn/In/Sc and PB vertices
            for &v in &a_fi.vertices_on { a_mvf.insert(v); }
            for &v in &a_fi.vertices_in { a_mvf.insert(v); }
            for &v in &a_fi.vertices_sc { a_mvf.insert(v); }

            // Also add vertices from PBs on the face
            for &pb_ptr in &face_pb_on {
                if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                    let pbr = pb.0.read().unwrap();
                    a_mvf.insert(pbr.pave1.vertex_idx);
                    a_mvf.insert(pbr.pave2.vertex_idx);
                }
            }
            for &pb_ptr in &face_pb_in {
                if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                    let pbr = pb.0.read().unwrap();
                    a_mvf.insert(pbr.pave1.vertex_idx);
                    a_mvf.insert(pbr.pave2.vertex_idx);
                }
            }
            for &pb_ptr in &face_pb_sc {
                if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                    let pbr = pb.0.read().unwrap();
                    a_mvf.insert(pbr.pave1.vertex_idx);
                    a_mvf.insert(pbr.pave2.vertex_idx);
                }
            }
            // Drop a_fi to release immutable borrow on self.ds
            // before mutable operations below
            drop(a_fi);

            // L947-1107: iterate all PBs and check for EF common blocks
            for &(n_e, local_i) in the_mpb {
                // L952-955: skip if PB already on the face.
                // Section edges (MakeBlocks) keep their PBs in the curve with an
                // orphan pool entry and no ShapeInfo reference, so edge_pave_blocks
                // is empty for them.
                if self.ds.edge_pave_blocks(n_e).len() <= local_i {
                    continue;
                }
                let a_pb = self.ds.edge_pave_blocks(n_e)[local_i].clone();
                // Check if already in face's sets (OCCT: Contains(aPB) on
                // handle-keyed sets 鈥?compare by PB pointer id).
                let a_pb_ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                let already_on_face = face_pb_on.contains(&a_pb_ptr)
                    || face_pb_in.contains(&a_pb_ptr)
                    || face_pb_sc.contains(&a_pb_ptr);
                if already_on_face {
                    continue;
                }

                // L958-964: check if face contains both vertices of PB
                let (n_v1, n_v2) = {
                    let pbr = a_pb.0.read().unwrap();
                    (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
                };
                if !a_mvf.contains(&n_v1) || !a_mvf.contains(&n_v2) {
                    continue;
                }

                // L966-981: get the edge
                let pbr = a_pb.0.read().unwrap();
                let n_e_actual = if pbr.edge != usize::MAX {
                    pbr.edge
                } else {
                    pbr.original_edge
                };
                if n_e_actual >= self.ds.nb_shapes() {
                    continue;
                }
                let rank_e = self.ds.rank(n_e_actual);
                let rank_f = self.ds.rank(n_f);
                // L977-980: if same rank, skip
                if rank_e >= 0 && rank_e == rank_f {
                    continue;
                }
                let a_range = pbr.range();
                drop(pbr);

                // L986-1052: edge-face coincidence check
                // OCCT: aBAC.D1(IntermediatePoint(aTS[0], aTS[1]), aPOnE, aVETgt)
                let curve = match self.ds.edge_curve(n_e_actual) {
                    Some(c) => c.clone(),
                    None => {
                        continue;
                    }
                };
                let mid_t = (a_range.0 + a_range.1) * 0.5;
                let mid_pt = curve.point_at(mid_t);
                // OCCT L1001-1006: tangent vector at midpoint
                let dt = 1e-7;
                let p_mid_dt = curve.point_at(mid_t + dt);
                let v_etgt = p_mid_dt - mid_pt;
                if v_etgt.length_squared() < 1e-7 {
                    continue;
                }

                // OCCT L1022-1024: aTolCheck = bSICheckMode ? myFuzzyValue :
                //   2 * max(BRep_Tool::Tolerance(aV1), BRep_Tool::Tolerance(aV2))
                let tol_v1 = self.ds.vertex_tolerance_by_idx(n_v1);
                let tol_v2 = self.ds.vertex_tolerance_by_idx(n_v2);
                let a_tol_check = if b_si_check_mode {
                    self.my_fuzzy_value
                } else {
                    2.0 * tol_v1.max(tol_v2)
                };

                // Project midpoint onto face surface (OCCT L1031-1036)
                let (proj_uv, proj_pt) = if let Some(surf) = self.ds.face_surface(n_f) {
                    let (uv, proj_pt) = crate::bop::closest_point_on_surface(&surf, mid_pt);
                    let a_dist = (proj_pt - mid_pt).length();
                    // OCCT L1026: if (LowerDistance() > aTolCheck + myFuzzyValue) continue;
                    if a_dist > a_tol_check + self.my_fuzzy_value {
                        continue;
                    }
                    // OCCT L1033-1035: if (!myContext->IsPointInFace(aF, gp_Pnt2d(U,V))) continue;
                    if !self.my_context.is_point_in_face(&self.ds, n_f, uv) {
                        continue;
                    }
                    (uv, proj_pt)
                } else { continue; };

                // OCCT L1038-1051: angle between face-to-edge vector and edge tangent
                // OCCT: if (aSurfAdaptor.GetType() != GeomAbs_Plane || aBAC.GetType() != GeomAbs_Line)
                // rcad: skip angle check when face is Plane AND edge is Line
                let mut b_use_add_tol = true;
                {
                    let surf_is_plane = self.ds.face_surface(n_f).map_or(false, |s| matches!(s, rcad_kernel::geom::Surface3::Plane(_)));
                    let curve_is_line = matches!(curve, rcad_kernel::geom::Curve3::Line(_));
                    if !(surf_is_plane && curve_is_line) {
                        let a_vf_norm = mid_pt - proj_pt;
                        if a_vf_norm.length_squared() > 1e-7 {
                            // OCCT L1046-1047: if (|aCos| > 0.4226) bUseAddTol = false
                            let a_cos = a_vf_norm.normalize().dot(v_etgt.normalize()).abs();
                            if a_cos > 0.4226 {
                                b_use_add_tol = false;
                            }
                        }
                    }
                }

                // Compute additional tolerance from endpoint distances (OCCT L1063-1084)
                let mut a_tol_add_ef = 0.0;
                if b_use_add_tol {
                    if let Some(surf) = self.ds.face_surface(n_f) {
                        for a_t in [a_range.0, a_range.1] {
                            let a_p = curve.point_at(a_t);
                            let (_uv_e, proj_pe) = crate::bop::closest_point_on_surface(&surf, a_p);
                            let a_dist_ef = (proj_pe - a_p).length();
                            if a_dist_ef < a_tol_check && a_dist_ef > a_tol_add_ef {
                                a_tol_add_ef = a_dist_ef;
                            }
                        }
                    }
                    // OCCT L1077-1084: subtract edge and face tolerance
                    if a_tol_add_ef > 0.0 {
                        let tol_e = self.ds.edge_tolerance(n_e_actual);
                        let tol_f = self.ds.face_tolerance(n_f);
                        a_tol_add_ef -= (tol_e + tol_f);
                        if a_tol_add_ef < 0.0 {
                            a_tol_add_ef = 0.0;
                        }
                    }
                }

                // OCCT L1087-1092: bIntersect = aTolAdd > 0, with myFPBDone fallback
                let mut b_intersect = a_tol_add_ef > 0.0;
                if !b_intersect {
                    if let Some(pmpb) = self.my_fpb_done.get(&n_f) {
                        let ptr = std::sync::Arc::as_ptr(&a_pb.0) as u64;
                        b_intersect = !pmpb.contains(&ptr);
                    } else {
                        b_intersect = true;
                    }
                }
                if !b_intersect {
                    continue;
                }

                // L1094-1106: prepare the pair for intersection
                let pb_pool_idx = {
                    let ptr = std::sync::Arc::as_ptr(&a_pb.0);
                    let mut found = usize::MAX;
                    for (&pi, pool) in self.ds.pave_blocks_pool.iter() {
                        for spb in pool {
                            if std::sync::Arc::as_ptr(&spb.0) == ptr {
                                found = pi;
                                break;
                            }
                        }
                        if found != usize::MAX { break; }
                    }
                    found
                };
                // (n_e, n_f, pb_pool_idx, a_tol_add, PB pave range)
                ef_pairs.push((n_e_actual, n_f, pb_pool_idx, a_tol_add_ef, a_range));
            }
        }

        // L1110-1114: if no pairs, return
        if ef_pairs.is_empty() {
            return;
        }

        // L1106-1107/L1122-1129: BOPAlgo_EdgeFace pairs 鈥?OCCT performs them
        // in parallel (BOPTools_Parallel); rcad runs them serially.
        let mut a_efs: Vec<(EdgeFace, usize)> = Vec::new(); // (edge_face, pb_pool_idx)
        for &(n_e, n_f, pb_pool_idx, a_tol_add, a_range) in &ef_pairs {
            let mut ef = match EdgeFace::new(&self.ds, n_e, n_f) {
                Some(ef) => ef,
                None => continue,
            };
            ef.set_range(a_range.0, a_range.1);
            ef.set_fuzzy_value(self.my_fuzzy_value + a_tol_add);
            ef.use_quick_coincidence_check(true);
            a_efs.push((ef, pb_pool_idx));
        }
        if a_efs.is_empty() {
            return;
        }
        for (ef, _) in &mut a_efs {
            ef.perform();
        }

        // L1141-1192: analyze the results 鈥?accept only a single
        // TopAbs_EDGE common part (OCCT L1177-1188).
        let mut a_mpbli: indexmap::IndexMap<u64, Vec<usize>> = indexmap::IndexMap::new();
        let mut a_warnings: Vec<(usize, usize)> = Vec::new();
        let mut a_ef_results: Vec<(usize, usize, usize, f64, DVec3)> = Vec::new();
        for (ef, pb_pool_idx) in &a_efs {
            let n_e = ef.edge_index();
            let n_f = ef.face_index();
            if !ef.is_done() || ef.error_status() != 0 {
                // OCCT L1171-1175: AddIntersectionFailedWarning(Edge(), Face()).
                a_warnings.push((n_e, n_f));
                continue;
            }
            let a_cparts = ef.common_parts();
            if a_cparts.len() != 1 {
                continue;
            }
            let a_cp = &a_cparts[0];
            if a_cp.get_type() != CommonPrtType::Edge {
                continue;
            }
            // OCCT L1180-1181: aEF.SetCommonPart(aCP) 鈥?the EDGE range; the
            // vertex is created from the range middle in MakeSDVerticesFF.
            let mid_t = (a_cp.range1[0] + a_cp.range1[1]) * 0.5;
            let mid_pt = ef.edge_value(mid_t);
            a_ef_results.push((n_e, n_f, *pb_pool_idx, mid_t, mid_pt));
        }
        drop(a_efs);

        for (n_e, n_f) in a_warnings {
            self.my_report.add_warning(Alert::IntersectionFailed(n_e, n_f));
        }

        for &(n_e, n_f, pb_pool_idx, mid_t, mid_pt) in &a_ef_results {
            if the_add_interf {
                // L1175-1181: BOPDS_InterfEF aEF = aEFs.Appended();
                let new_ef = InterferenceEF {
                    edge: n_e,
                    face: n_f,
                    point: mid_pt,
                    edge_param: mid_t,
                    new_vertex: usize::MAX,
                };
                self.ds.interf_ef.push(new_ef);
                self.ds.add_interf(n_e, n_f);
            }

            // L1184-1186: myDS->ChangeFaceInfo(nF).ChangePaveBlocksIn().Add(aPB);
            if let Some(pool) = self.ds.pave_blocks_pool.get(&pb_pool_idx) {
                let pb_ptrs: Vec<u64> = pool.iter()
                    .map(|pb| std::sync::Arc::as_ptr(&pb.0) as u64).collect();
                drop(pool);
                for pb_ptr in pb_ptrs {
                    // OCCT: Add(aPB) 鈥?keyed by PB handle.
                    self.ds.change_face_info(n_f).pave_blocks_in.insert(pb_ptr);
                    // Record PB鈫抐ace mapping for CommonBlock creation
                    a_mpbli.entry(pb_ptr).or_default().push(n_f);
                }
            }
        }

        // L1194-1198: BOPAlgo_Tools::PerformCommonBlocks(aMPBLI, anAlloc, myDS)
        // Create CommonBlocks for each PB鈫抐ace association (OCCT overload 2)
        for (ptr, faces) in &a_mpbli {
            // Find the SharedPB from its pointer
            let mut pb_found: Option<SharedPB> = None;
            'outer: for pool in self.ds.pave_blocks_pool.values() {
                for spb in pool {
                    if std::sync::Arc::as_ptr(&spb.0) as u64 == *ptr {
                        pb_found = Some(spb.clone());
                        break 'outer;
                    }
                }
            }
            if let Some(pb) = pb_found {
                // Check if PB is already in a CommonBlock 鈫?reuse (OCCT L206-208)
                let cb_idx = pb.0.read().unwrap().common_block_idx;
                if let Some(cb_idx) = cb_idx {
                    if cb_idx < self.ds.common_blocks.len() {
                        self.ds.common_blocks[cb_idx].append_faces(faces);
                    }
                } else {
                    // Create new CommonBlock with PB and set faces (OCCT L211-214, 238)
                    let mut a_cb = crate::bop::ds::common_block::CommonBlock::new();
                    a_cb.add_pave_block(pb.clone(), 0); // OCCT: AddPaveBlock(aPB), face_idx is 0 placeholder
                    for &f in faces {
                        a_cb.add_face(f);
                    }
                    let cb_idx = self.ds.common_blocks.len();
                    self.ds.common_blocks.push(a_cb);
                    self.ds.set_common_block(&pb, cb_idx);
                }
            }
        }
    }

    // ====================================================================
    // CheckSelfInterference 鈥?OCCT BOPAlgo_PaveFiller_11.cxx L28-221
    // ====================================================================
    /// Check for acquired self-intersections after intersection processing.
    /// OCCT BOPAlgo_PaveFiller::CheckSelfInterference (PaveFiller_11.cxx L28-221).
    fn check_self_interference(&mut self) {
        // L30-34: if (myArguments.Extent() == 1) return;
        if self.my_arguments.len() <= 1 {
            return;
        }

        // L36: BRep_Builder aBB;
        // L38: int i, aNbR = myDS->NbRanges();
        let a_nb_r = self.ds.nb_ranges();
        // L39: for (i = 0; i < aNbR; ++i)
        for a_rank in 0..a_nb_r {
            // L41: const BOPDS_IndexRange& aR = myDS->Range(i);
            let a_r = self.ds.range(a_rank);

            // L44-47: NCollection_IndexedDataMap<TopoDS_Shape, IndexedMap<TopoDS_Shape>> aMCSI
            // 鈥?insertion order; iterated at L222-226 by index.
            let mut a_mcsi: indexmap::IndexMap<usize, indexmap::IndexSet<usize>> =
                indexmap::IndexMap::new();
            // L48: NCollection_Map<CommonBlock> aMCBFence;
            let mut a_cb_fence: std::collections::HashSet<usize> =
                std::collections::HashSet::new();

            // L51: for (j = aR.First(); j <= aR.Last(); ++j)
            for j in a_r.first..=a_r.last {
                // L53-54: check HasReference
                if self.ds.shapes[j].reference < 0 {
                    continue;
                }

                // L62-63: check ShapeType
                if self.ds.shapes[j].shape_type == ShapeType::Edge {
                    // L65-67: skip degenerated edges
                    if self.ds.shapes[j].has_flag() {
                        continue;
                    }

                    // L71-78: analyze shared vertices
                    let mut a_sub_s: std::collections::HashSet<usize> =
                        std::collections::HashSet::new();
                    for &n_v in &self.ds.shapes[j].sub_shapes {
                        // OCCT L75: int nV = aItLI.Value();
                        // L76: myDS->HasShapeSD(nV, nV); 鈥?replaces nV with SD if exists
                        let mut n_vx = n_v;
                        let mut n_vsd = usize::MAX;
                        if self.ds.has_shape_sd(n_v, &mut n_vsd) {
                            n_vx = n_vsd;
                        }
                        // L77: aMSubS.Add(nV);
                        a_sub_s.insert(n_vx);
                    }

                    // L80-81: PaveBlocks for this edge
                    let a_lpb = self.ds.edge_pave_blocks(j);
                    let b_analyze_v = a_lpb.len() > 1;

                    // L83-149: iterate PBs
                    for spb in a_lpb {
                        let pb = spb.0.read().unwrap();

                        // L89-109: check the vertices
                        if b_analyze_v {
                            let (nv1, nv2) = pb.indices();
                            for &n_v in &[nv1, nv2] {
                                // L95: if (!aR.Contains(nV[k]) && !aMSubS.Contains(nV[k]))
                                let in_range = n_v >= a_r.first && n_v <= a_r.last;
                                if !in_range && !a_sub_s.contains(&n_v) {
                                    // L97-106: add connection
                                    a_mcsi.entry(n_v).or_default().insert(j);
                                }
                            }
                        }

                        // L112-148: check common blocks
                        if let Some(cb_idx) = pb.common_block_idx {
                            if a_cb_fence.insert(cb_idx) {
                                if let Some(cb) = self.ds.common_blocks.get(cb_idx) {
                                    let mut a_le: Vec<usize> = Vec::new();
                                    for (pool_pb, _) in cb.pave_blocks() {
                                        let n_e_or = pool_pb.0.read().unwrap().original_edge;
                                        // L125: if (aR.Contains(nEOr))
                                        let in_range = n_e_or >= a_r.first && n_e_or <= a_r.last;
                                        if in_range {
                                            a_le.push(n_e_or);
                                        }
                                    }
                                    // L131-146: if more than 1 edge from same argument in CB
                                    if a_le.len() > 1 {
                                        self.my_report.add_warning(
                                            Alert::AcquiredSelfIntersection(a_le));
                                    }
                                }
                            }
                        }
                    }
                } else if self.ds.shapes[j].shape_type == ShapeType::Face {
                    // L151-196: analyze FACE
                    // L155: const BOPDS_FaceInfo& aFI = myDS->FaceInfo(j);
                    let a_fi = self.ds.face_info(j);

                    // L156-173: IN and Section vertices
                    for &n_v in &a_fi.vertices_in {
                        a_mcsi.entry(n_v).or_default().insert(j);
                    }
                    for &n_v in &a_fi.vertices_sc {
                        a_mcsi.entry(n_v).or_default().insert(j);
                    }

                    // L175-195: IN and Section PaveBlocks
                    for &pb_ptr in &a_fi.pave_blocks_in {
                        if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                            let n_e = pb.0.read().unwrap().edge;
                            if n_e != usize::MAX {
                                a_mcsi.entry(n_e).or_default().insert(j);
                            }
                        }
                    }
                    for &pb_ptr in &a_fi.pave_blocks_sc {
                        if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                            let n_e = pb.0.read().unwrap().edge;
                            if n_e != usize::MAX {
                                a_mcsi.entry(n_e).or_default().insert(j);
                            }
                        }
                    }
                }
            }

            // L200-219: Analyze connections 鈥?if a vertex/edge connects
            // to multiple faces from same argument, add warning
            for (_sub_shape, shapes) in &a_mcsi {
                if shapes.len() > 1 {
                    self.my_report.add_warning(
                        Alert::AcquiredSelfIntersection(
                            shapes.iter().copied().collect()));
                }
            }
        }
    }

    /// OCCT BOPAlgo_PaveFiller::MakeSplitEdges (_7.cxx L371-548).
    fn make_split_edges(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L392: UpdateCommonBlocksWithSDVertices
        for cb in &self.ds.common_blocks {
            self.ds.update_common_block_with_sd_vertices(cb);
        }

        if self.ds.pave_blocks_pool.is_empty() { return; }
        // OCCT L386: aMCB fence for CommonBlocks
        let mut a_mcb: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let ms_debug = std::env::var("RCAD_MS_DEBUG").is_ok();
        let nb_src = self.ds.nb_source_shapes();
        // OCCT myPaveBlocksPool is a DynamicArray indexed by edge 鈥?iterate
        // in ascending key order (a HashMap would randomize the split-edge
        // creation order and hence the DS edge indices).
        let mut pb_keys: Vec<usize> = self.ds.pave_blocks_pool.keys().copied().collect();
        pb_keys.sort_unstable();
        for k in pb_keys {
            let a_lpb: Vec<SharedPB> = match self.ds.pave_blocks_pool.get(&k) {
                Some(v) => v.clone(),
                None => continue,
            };
            for a_pb in &a_lpb {
                // OCCT L410-414: skip degenerated edges
                let pb = a_pb.0.read().unwrap();
                let n_e = pb.original_edge;
                if n_e >= self.ds.nb_shapes() {
                    if ms_debug { eprintln!("[RCADMS] SKIP-OOB edge={} nshapes={}", n_e, self.ds.nb_shapes()); }
                    drop(pb); continue;
                }
                if self.ds.shapes[n_e].has_flag() {
                    if ms_debug { eprintln!("[RCADMS] SKIP-FLAG edge={}", n_e); }
                    drop(pb); continue;
                }
                // OCCT L416-421: skip if already processed via CB fence
                if let Some(cb_idx) = pb.common_block_idx {
                    if !a_mcb.insert(cb_idx) {
                        if ms_debug { eprintln!("[RCADMS] SKIP-CB-FENCE edge={} cb={}", n_e, cb_idx); }
                        drop(pb); continue;
                    }
                }
                let n_v1 = pb.pave1.vertex_idx;
                let n_v2 = pb.pave2.vertex_idx;
                let b_v1 = n_v1 >= nb_src;
                let b_v2 = n_v2 >= nb_src;
                let cb_f = pb.common_block_idx;
                let a_t1 = pb.pave1.param;
                let a_t2 = pb.pave2.param;
                // Release the read lock on a_pb before any CB SetEdge may write
                // the same PB (OCCT has no locks; aCB->SetEdge writes all CB PBs).
                drop(pb);
                if ms_debug {
                    eprintln!("[RCADMS] CAND edge={} v1={} v2={} newV1={} newV2={} CB={} lpbExtent={}",
                        n_e, n_v1, n_v2, b_v1, b_v2, cb_f.is_some(), a_lpb.len());
                }
                // OCCT L429-450: check if it is necessary to make the split of the edge
                let mut b_to_split = true;
                let mut set_edge_n = usize::MAX; // OCCT L460: aPB->SetEdge(nE)
                if !b_v1 && !b_v2 {
                    // OCCT L432: if (!myNonDestructive || !bCB)
                    if !self.my_non_destructive || cb_f.is_none() {
                        let mut it_found = false;
                        let mut found_e = usize::MAX;
                        if let Some(cb_idx) = cb_f {
                            // OCCT L436-445: find the edge with these vertices in the
                            //   common block whose PaveBlocks extent == 1
                            let cb_pbs: Vec<SharedPB> = {
                                let cb = &self.ds.common_blocks[cb_idx];
                                cb.pave_blocks().iter().map(|(pb, _)| pb.clone()).collect()
                            };
                            for pbx in &cb_pbs {
                                let e = pbx.0.read().unwrap().original_edge;
                                if self.ds.pave_blocks(e).len() == 1 {
                                    it_found = true;
                                    found_e = e;
                                    break;
                                }
                            }
                        }
                        if it_found {
                            // OCCT L446-455: the pave block is found 鈥?no split.
                            //   aCB->SetRealPaveBlock(it.Value()); aCB->SetEdge(nE);
                            //   ComputeToleranceOfCB + UpdateEdgeTolerance.
                            b_to_split = false;
                            if let Some(cb_idx) = cb_f {
                                self.ds.common_blocks[cb_idx].set_edge(found_e);
                                let a_tol =
                                    Self::compute_tolerance_of_cb(cb_idx, &self.ds);
                                self.update_edge_tolerance(found_e, a_tol);
                            }
                        } else if cb_f.is_none() && a_lpb.len() == 1 {
                            // OCCT L457-461: no common block, single-PB edge 鈥?no split
                            b_to_split = false;
                            set_edge_n = n_e;
                        }
                    }
                }
                if !b_to_split {
                    if set_edge_n != usize::MAX {
                        // OCCT L460: aPB->SetEdge(nE)
                        a_pb.0.write().unwrap().edge = set_edge_n;
                    }
                    if ms_debug { eprintln!("[RCADMS] NOSPLIT edge={} v1={} v2={}", n_e, n_v1, n_v2); }
                    continue;
                }
                if ms_debug { eprintln!("[RCADMS] SPLIT edge={} v1={} v2={}", n_e, n_v1, n_v2); }
                // OCCT L484-490: when the PB is in a common block, split
                //   aCB->PaveBlock1() instead of the current PB.
                let (sp_e, sp_v1, sp_v2, sp_t1, sp_t2) = if let Some(cb_idx) = cb_f {
                    if let Some(pb1) = self.ds.common_blocks[cb_idx].pave_block1() {
                        let r = pb1.0.read().unwrap();
                        (r.original_edge, r.pave1.vertex_idx, r.pave2.vertex_idx,
                         r.pave1.param, r.pave2.param)
                    } else {
                        (n_e, n_v1, n_v2, a_t1, a_t2)
                    }
                } else {
                    (n_e, n_v1, n_v2, a_t1, a_t2)
                };
                // OCCT L465-515: create new split edge
                if let Some(curve) = self.ds.edge_curve(sp_e) {
                    let new_ei = self.ds.push_edge_inherit(curve.clone(), [sp_t1, sp_t2], sp_v1, sp_v2, Some(sp_e));
                    if let Some(cb_idx) = cb_f {
                        // OCCT L536-545: aCBk->SetEdge(nSp) 鈥?sets the edge of all
                        //   PaveBlocks in the common block.
                        self.ds.common_blocks[cb_idx].set_edge(new_ei);
                    } else {
                        let mut pbw = a_pb.0.write().unwrap();
                        pbw.edge = new_ei;
                    }
                    if ms_debug { eprintln!("[RCADMS] NEWEDGE nSp={} origEdge={} p1={} p2={}", new_ei, sp_e, sp_v1, sp_v2); }
                    // OCCT BOPAlgo_SplitEdge::Perform (_7.cxx L137-138):
                    //   BRepBndLib::Add(myESp, myBox); myBox.SetGap(gap + Confusion())
                    self.rebuild_edge_box(new_ei);
                }
            }
        }
        // OCCT L534-550: FillShrunkData for new PBs
        // rcad: shrunk data computed in FillShrunkData step.
    }

    /// OCCT BOPAlgo_PaveFiller::RemoveMicroEdges (_6.cxx L4388-4435).
    fn remove_micro_edges(&mut self) {
        // OCCT L4450: aMicroEdges is NCollection_Map<int> 鈥?bucket order.
        let mut a_micro_edges: crate::bop::algo::occt_map::OcctMapInt =
            crate::bop::algo::occt_map::OcctMapInt::new();
        let mut a_mpb_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // OCCT L4392-4393: ChangePaveBlocksPool() is a DynamicArray 鈥?iterate
        // the pool in ascending edge-key order.
        let mut pb_keys: Vec<usize> = self.ds.pave_blocks_pool.keys().copied().collect();
        pb_keys.sort_unstable();
        for k in pb_keys {
            let pb_list = match self.ds.pave_blocks_pool.get(&k) {
                Some(v) => v.clone(),
                None => continue,
            };
            if pb_list.len() < 2 { continue; }
            // OCCT L4407-4410: skip degenerated edges
            if pb_list.is_empty() { continue; }
            let n_e_orig = pb_list[0].0.read().unwrap().original_edge;
            if n_e_orig < self.ds.nb_shapes() && self.ds.shapes[n_e_orig].has_flag() {
                continue;
            }
            for pb in &pb_list {
                let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
                if !a_mpb_fence.insert(ptr) { continue; }
                let (n_v1, n_v2) = { let r = pb.0.read().unwrap(); r.indices() };
                if n_v1 == n_v2 {
                    self.fill_shrunk_data_pb(pb);
                    let has_shrunk = { let r = pb.0.read().unwrap(); r.has_shrunk_data() };
                    if !has_shrunk {
                        let e = { let r = pb.0.read().unwrap(); r.edge };
                        if e != usize::MAX {
                            a_micro_edges.add(e);
                        }
                    }
                }
            }
        }
        // OCCT L4434: RemovePaveBlocks(aMicroEdges) 鈥?three steps:
        //   1. from the Pave Blocks Pool
        //   2. from section curves (BOPDS_Curve::PaveBlocks)
        //   3. from Face Info (PaveBlocksIn/On/Sc)
        // rcad previously only did step 1, leaving dangling PB references in
        // intersection_curves[cid].pave_blocks and the FaceInfo sets.
        let mut a_micro_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ei in a_micro_edges.iter_keys() {
            a_micro_set.insert(ei);
        }
        self.remove_pave_blocks(&a_micro_set);
    }

    /// OCCT BOPAlgo_PaveFiller::MakePCurves (_7.cxx L589-850).
    fn make_pcurves(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // OCCT L592-595: myAvoidBuildPCurve || (!PCurveOnS1 && !PCurveOnS2) -> return.
        // rcad: PCurveOnS1/2 default true; the avoid flag is not present.
        // OCCT L606-700: 1. Process face info 鈥?IN and ON PBs.
        let a_nb_fi = self.ds.face_info_pool.len();
        for fi_idx in 0..a_nb_fi {
            let fi = self.ds.face_info_pool[fi_idx].clone();
            let n_f1 = fi.index();
            let f1_s = self.ds.shape(n_f1).clone();
            let surf = match &*f1_s.data {
                rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
                _ => continue,
            };
            let Some(ref surf) = surf else { continue; };

            // OCCT L619-631: PaveBlocksIn 鈥?pcurve by projection.
            for &pb_ptr in fi.pave_blocks_in.iter() {
                if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                    let n_e = pb.0.read().unwrap().edge;
                    if n_e >= self.ds.nb_shapes() { continue; }
                    self.build_pcurve_mpc(n_e, n_f1, surf, None);
                }
            }
            // OCCT L634-699: PaveBlocksOn 鈥?skip if pcurve exists; a CommonBlock
            // provides the pcurve-copy source (paired edge with a pcurve).
            for &pb_ptr in fi.pave_blocks_on.iter() {
                if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                    let n_e = pb.0.read().unwrap().edge;
                    if n_e >= self.ds.nb_shapes() { continue; }
                    if self.edge_has_pcurve(n_e, n_f1) { continue; }
                    let src = self.cb_pcurve_source(&pb, n_f1);
                    self.build_pcurve_mpc(n_e, n_f1, surf, src);
                }
            }
        }
        // OCCT L702-850: 2. Process section edges. P-curves on them must already
        // be computed; the MPCs still provide the UpdateVertices call (flag).
        let mut an_ef_pairs: HashSet<(usize, usize)> = HashSet::new();
        let a_ffs = self.ds.interf_ff.clone();
        for a_ff in &a_ffs {
            let n_f = [a_ff.f1, a_ff.f2];
            for &cid in &a_ff.curves {
                if cid >= self.ds.intersection_curves.len() { continue; }
                let a_lpb = self.ds.intersection_curves[cid].pave_blocks.clone();
                for pb in &a_lpb {
                    let n_e = pb.0.read().unwrap().edge;
                    if n_e >= self.ds.nb_shapes() { continue; }
                    for m in 0..2 {
                        if !an_ef_pairs.insert((n_e, n_f[m])) { continue; }
                        // OCCT L744-752: the section-edge MPC has SetFlag(true);
                        // in Perform, when the edge already has a pcurve on the
                        // face (attached by MakePCurve at section creation,
                        // BOPAlgo_PaveFiller_6.cxx L1066-1072) it is kept 鈥?only
                        // UpdateVertices is called. Never overwrite it with a
                        // projection. Edges without a pcurve (null FirstCurve2d)
                        // are projected here, matching the MPC null-branch.
                        if self.edge_has_pcurve(n_e, n_f[m]) { continue; }
                        let f_s = self.ds.shape(n_f[m]).clone();
                        let surf2 = match &*f_s.data {
                            rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
                            _ => continue,
                        };
                        if let Some(ref surf2) = surf2 {
                            self.build_pcurve_mpc(n_e, n_f[m], surf2, None);
                        }
                    }
                }
            }
        }
    }

    /// OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the pcurve of edge
    /// n_e on face n_f is keyed by (face TShape, L.Predivided(E.Location())) 鈥?
    /// the face location divided by the edge's location.
    fn pcurve_key_for(&self, n_e: usize, n_f: usize) -> Option<(u64, u32)> {
        let (fid, floc) = self.ds.face_key(n_f)?;
        let eloc = self.ds.shape(n_e).location;
        Some((
            fid,
            crate::bop::algo::compose_face_edge_pcurve_location(floc, eloc, &self.ds.locations),
        ))
    }

    /// OCCT BOPAlgo_MPC::Perform (BOPAlgo_PaveFiller_7.cxx L218-293). Computes
    /// the pcurve of an edge on a face. With a copy source the pcurve is taken
    /// from the paired edge (AttachExistingPCurve); otherwise it is projected
    /// (BuildPCurveForEdgeOnFace). In-place (OCCT BRep_Builder semantics) 鈥?safe
    /// because the DS owns a private input copy.
    fn build_pcurve_mpc(&mut self, n_e: usize, n_f: usize, surf: &Surface3, src: Option<usize>) {
        // OCCT L233-236: if the edge already has a pcurve on the face, do not
        // rebuild it (only the periodic adjustment, which a line pcurve at the
        // seam does not need). Input-shape pcurves are keyed by the face's
        // TShape identity (ptr_id, location), same as here.
        if self.edge_has_pcurve(n_e, n_f) {
            return;
        }
        let Some(fkey) = self.pcurve_key_for(n_e, n_f) else { return };
        // OCCT L239-249: attach the pcurve from the paired edge.
        if let Some(src_e) = src {
            if let Some(pc) = Self::copy_pcurve(&self.ds, src_e, n_e, n_f) {
                let range = self.ds.edge_range(n_e);
                self.ds.mutate_shape_data(n_e, |ts| {
                    if let rcad_kernel::topods::TShape::Edge(ed) = ts {
                        ed.pcurves.insert(fkey, (pc, range[0], range[1]));
                    }
                });
                self.ds.remap_shape_idx(n_e);
                return;
            }
        }
        // OCCT L248: BuildPCurveForEdgeOnFace 鈥?projection.
        if let Some(curve) = self.ds.edge_curve(n_e) {
            let range = self.ds.edge_range(n_e);
            if let Some(pc) = Self::pcurve_2d(&curve, surf, range) {
                self.ds.mutate_shape_data(n_e, |ts| {
                    if let rcad_kernel::topods::TShape::Edge(ed) = ts {
                        ed.pcurves.insert(fkey, (pc, range[0], range[1]));
                    }
                });
                self.ds.remap_shape_idx(n_e);
            }
        }
    }

    /// OCCT BOPTools_AlgoTools2D::HasCurveOnSurface 鈥?the edge already has a
    /// pcurve on the face. OCCT matches the face's TShape identity
    /// (BRep_Tool::CurveOnSurface, keyed by L.Predivided(E.Location()));
    /// rcad keys the edge pcurve map by the same identity, so input-shape
    /// pcurves and DS-built ones share the map.
    fn edge_has_pcurve(&self, n_e: usize, n_f: usize) -> bool {
        if n_e >= self.ds.nb_shapes() { return false; }
        let Some(fkey) = self.pcurve_key_for(n_e, n_f) else { return false };
        match &*self.ds.shape(n_e).data {
            rcad_kernel::topods::TShape::Edge(ed) => ed.pcurves.contains_key(&fkey),
            _ => false,
        }
    }

    /// OCCT BOPAlgo_Tools::ComputeToleranceOfCB (BOPAlgo_Tools.cxx L248-356):
    /// the max tolerance of a common block 鈥?max of the real edge tolerance,
    /// the sampled distances to the other block edges (projected points), and
    /// the sampled distances to the block's faces.
    fn compute_tolerance_of_cb(cb_idx: usize, ds: &DS) -> f64 {
        let cb = &ds.common_blocks[cb_idx];
        let mut a_tol_max = 0.0;
        let Some(a_pbr) = cb.pave_block1() else { return a_tol_max; };
        let n_e = a_pbr.0.read().unwrap().original_edge();
        let a_e_or = ds.shape(n_e);
        a_tol_max = a_e_or.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
        let a_lpb: Vec<SharedPB> = cb.pave_blocks().iter().map(|(p, _)| p.clone()).collect();
        let a_lfi = cb.faces();
        if a_lpb.len() < 2 && a_lfi.is_empty() {
            return a_tol_max;
        }
        let a_nb_pnt = 11usize;
        let a_pbr_range = a_pbr.0.read().unwrap().range();
        let (mut a_t1, a_t2) = (a_pbr_range.0, a_pbr_range.1);
        let a_dt = (a_t2 - a_t1) / (a_nb_pnt as f64 + 1.0);
        let Some(a_c3d) = ds.edge_curve(n_e) else { return a_tol_max; };
        // max distance between the edges.
        for pb in &a_lpb {
            if std::sync::Arc::ptr_eq(&pb.0, &a_pbr.0) {
                continue;
            }
            let n_e2 = pb.0.read().unwrap().original_edge();
            let a_tol = ds.edge_tolerance(n_e2);
            let Some(e2_curve) = ds.edge_curve(n_e2) else { continue };
            a_t1 = a_pbr_range.0;
            for _ in 1..=a_nb_pnt {
                a_t1 += a_dt;
                let a_p = a_c3d.point_at(a_t1);
                let (_, a_proj) = crate::bop::closest_point_on_curve(&e2_curve, a_p);
                let a_tol_new = a_tol + (a_p - a_proj).length();
                if a_tol_new > a_tol_max {
                    a_tol_max = a_tol_new;
                }
            }
        }
        // max distance to the faces.
        for &n_f in a_lfi {
            let a_tol = ds.shape(n_f).as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
            let Some(surf) = ds.face_surface(n_f) else { continue };
            a_t1 = a_pbr_range.0;
            for _ in 1..=a_nb_pnt {
                a_t1 += a_dt;
                let a_p = a_c3d.point_at(a_t1);
                let (_, a_proj) = crate::bop::closest_point_on_surface(&surf, a_p);
                let a_tol_new = a_tol + (a_p - a_proj).length();
                if a_tol_new > a_tol_max {
                    a_tol_max = a_tol_new;
                }
            }
        }
        a_tol_max
    }

    /// OCCT L640-676: the CommonBlock copy source 鈥?a paired PB in the same
    /// CommonBlock whose original edge already has a pcurve on this face.
    fn cb_pcurve_source(&self, pb: &SharedPB, n_f: usize) -> Option<usize> {
        let cb_idx = pb.0.read().unwrap().common_block_idx?;
        let pbs: Vec<SharedPB> = {
            let cb = &self.ds.common_blocks[cb_idx];
            cb.pave_blocks().iter().map(|(p, _)| p.clone()).collect()
        };
        if pbs.len() < 2 { return None; }
        let pb_ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
        for pbx in &pbs {
            if std::sync::Arc::as_ptr(&pbx.0) as u64 == pb_ptr { continue; }
            let n_ex = pbx.0.read().unwrap().original_edge;
            if n_ex >= self.ds.nb_shapes() { continue; }
            if self.edge_has_pcurve(n_ex, n_f) {
                return Some(n_ex);
            }
        }
        None
    }

    /// OCCT BOPTools_AlgoTools2D::AttachExistingPCurve (BOPTools_AlgoTools2D_1.cxx
    /// L44-130): copy the paired edge's pcurve onto the section edge, reversing
    /// the direction when the section edge is split to reverse (IsSplitToReverse).
    /// The pcurve is re-sampled over the section edge's range (GeomLib::SameRange)
    /// by mapping 3D positions.
    fn copy_pcurve(ds: &DS, src_e: usize, dst_e: usize, n_f: usize) -> Option<rcad_kernel::geom::Curve2d> {
        use rcad_kernel::geom::Curve2dEval;
        let (src_pc, src_first, src_last) = Self::edge_pcurve_of(ds, src_e, n_f)?;
        let src_curve = ds.edge_curve(src_e)?;
        let dst_curve = ds.edge_curve(dst_e)?;
        let dst_range = ds.edge_range(dst_e);
        let dst_shape = ds.shape(dst_e).clone();
        let src_shape = ds.shape(src_e).clone();
        let (to_reverse, _) = algo_tools::is_split_to_reverse_edge(&dst_shape, &src_shape);
        let n = 23usize;
        let mut uv: Vec<glam::DVec2> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = dst_range[0] + i as f64 * (dst_range[1] - dst_range[0]) / n as f64;
            let p = dst_curve.point_at(t);
            let sp = closest_point_on_curve_range(&src_curve, p, src_first, src_last, 64);
            // OCCT L68-79: the copied pcurve is reversed when IsSplitToReverse.
            let s = if to_reverse {
                src_first + src_last - sp.param
            } else {
                sp.param
            };
            uv.push(Curve2dEval::point_at(&src_pc, s));
        }
        if uv.len() < 2 { return None; }
        let mut c = rcad_kernel::geom::BSplineCurve2::approximate(&uv);
        if dst_range[1] > dst_range[0] {
            c.knots = c.knots.iter().map(|k| dst_range[0] + (dst_range[1] - dst_range[0]) * k).collect();
        }
        Some(rcad_kernel::geom::Curve2d::BSpline(c))
    }

    /// BRep_Tool::CurveOnSurface 鈥?the edge's pcurve on the face, if any.
    /// OCCT keys the pcurve by (face TShape, L.Predivided(E.Location()))
    /// (BRep_Tool.cxx L345).
    fn edge_pcurve_of(ds: &DS, n_e: usize, n_f: usize) -> Option<(rcad_kernel::geom::Curve2d, f64, f64)> {
        if n_e >= ds.nb_shapes() { return None; }
        let (fid, floc) = ds.face_key(n_f)?;
        let eloc = ds.shape(n_e).location;
        let key = (fid, crate::bop::algo::compose_face_edge_pcurve_location(floc, eloc, &ds.locations));
        match &*ds.shape(n_e).data {
            rcad_kernel::topods::TShape::Edge(ed) => ed.pcurves.get(&key).cloned(),
            _ => None,
        }
    }


    /// Compute a 2D pcurve by projecting a 3D curve onto a surface.
    fn pcurve_2d(curve: &rcad_kernel::geom::Curve3,
                 surf: &rcad_kernel::geom::Surface3,
                 range: [f64; 2]) -> Option<rcad_kernel::geom::Curve2d> {
        use rcad_kernel::geom::SurfaceEval;
        let n = 23usize;
        let dt = (range[1] - range[0]) / n as f64;
        let mut uv: Vec<glam::DVec2> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = range[0] + i as f64 * dt;
            let p3d = curve.point_at(t);
            let (u, _) = crate::bop::closest_point_on_surface(surf, p3d);
            uv.push(u);
        }
        if uv.len() < 2 { return None; }
        // OCCT section-edge pcurves come from IntPatch: each arc keeps the u
        // domain of the side of the seam it lies on (an arc on the u=0 side
        // spans [0, u2], on the u=2PI side [u1, 2PI]). Projecting the 3D arc
        // samples the azimuth in [0, 2PI), so an arc crossing the seam wraps
        // (e.g. 2PI-eps ... 0.12). Unwrap the samples across the periodic
        // boundary, shift the domain back into [0, 2PI], and snap a sub-eps
        // endpoint to 0 鈥?matching OCCT's per-arc domain. Without this the
        // WireSplitter's periodic UV comparison sees u=2PI where OCCT sees
        // u=0 and takes the wrong edge at a seam vertex (P014). Only periodic
        // surfaces are unwrapped 鈥?a planar face's u is not periodic and a
        // large span must not be folded back by 2PI.
        let is_periodic_u = match surf {
            rcad_kernel::geom::Surface3::Plane(_) => false,
            _ => true,
        };
        if is_periodic_u {
            let two_pi = std::f64::consts::TAU;
            for i in 1..uv.len() {
                let du = uv[i].x - uv[i - 1].x;
                if du > std::f64::consts::PI {
                    uv[i].x -= two_pi;
                } else if du < -std::f64::consts::PI {
                    uv[i].x += two_pi;
                }
            }
            let min_u = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
            let max_u = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
            if max_u > two_pi {
                for p in uv.iter_mut() { p.x -= two_pi; }
            } else if min_u < 0.0 {
                for p in uv.iter_mut() { p.x += two_pi; }
            }
            for p in uv.iter_mut() {
                if p.x.abs() < 1e-12 { p.x = 0.0; }
            }
        }
        let mut c = rcad_kernel::geom::BSplineCurve2::approximate(&uv);
        // approximate() normalizes the parameterization to [0,1] (chord-length).
        // OCCT pcurves are parameterized by the edge's 3D range, and the
        // WireSplitter evaluates curve2d_point(pc, t) with the vertex parameter t
        // in that range 鈥?rescale the knots so the stored range is the curve's
        // own domain.
        let (r0, r1) = (range[0], range[1]);
        if r1 > r0 {
            c.knots = c.knots.iter().map(|k| r0 + (r1 - r0) * k).collect();
        }
        Some(rcad_kernel::geom::Curve2d::BSpline(c))
    }

    /// OCCT BOPAlgo_PaveFiller::ProcessDE (_8.cxx L54-131).
    /// OCCT BOPAlgo_PaveFiller::ProcessDE (_8.cxx L54-131).
    fn process_de(&mut self, the_range: &ProgressScope) {
        if the_range.user_break() { return; }
        // L62-63: for (int anEdgeIndex = 0; anEdgeIndex < myDS->NbSourceShapes(); ++anEdgeIndex)
        let a_nb_s = self.ds.nb_source_shapes();
        for an_ei in 0..a_nb_s {
            // L64-71: EDGE + HasFlag(nF)
            let ei = self.ds.shape_info(an_ei);
            if ei.shape_type != ShapeType::Edge { continue; }
            let n_f = ei.flag;
            if n_f < 0 { continue; }
            let n_f = n_f as usize;

            // L72-77: first sub-shape vertex, resolve SD
            let sf_type = self.ds.shape_info(n_f).shape_type;
            let n_v = ei.sub_shapes.first().copied().unwrap_or(usize::MAX);
            let mut n_vx = n_v;
            {
                let mut n_vsd = usize::MAX;
                if self.ds.has_shape_sd(n_vx, &mut n_vsd) {
                    n_vx = n_vsd;
                }
            }

            if sf_type == ShapeType::Face {
                // OCCT L82-84: FindPaveBlocks(nV, nF, aLPBOut) 鈥?the face's
                // PBs passing through the degenerated edge's (SD-resolved)
                // vertex.
                let a_fi = self.ds.face_info(n_f);
                let mut found_pbs: Vec<SharedPB> = Vec::new();
                for pb_set in [&a_fi.pave_blocks_in, &a_fi.pave_blocks_sc, &a_fi.pave_blocks_on] {
                    for &pb_ptr in pb_set {
                        if let Some(pb) = self.ds.pb_from_ptr(pb_ptr) {
                            let (v1, v2) = { let r = pb.0.read().unwrap(); r.indices() };
                            if v1 == n_vx || v2 == n_vx {
                                found_pbs.push(pb.clone());
                            }
                        }
                    }
                }
                drop(a_fi);
                if !found_pbs.is_empty() {
                    // OCCT L86-88: aLPBD = myDS->ChangePaveBlocks(anEdgeIndex);
                    // aPBD = aLPBD.First()
                    let a_lpbd: Vec<SharedPB> = self.ds.pave_blocks(an_ei).to_vec();
                    let a_pbd = match a_lpbd.first() {
                        Some(p) => p.clone(),
                        None => continue,
                    };
                    // The degenerated edge's 2D pcurve on the face.
                    let de_pcurve = {
                        let tshape = &self.ds.shapes[an_ei].shape.data;
                        let fkey = self.pcurve_key_for(an_ei, n_f);
                        match &**tshape {
                            rcad_kernel::topods::TShape::Edge(ed) => {
                                fkey.and_then(|k| ed.pcurves.get(&k).cloned())
                            }
                            _ => None,
                        }
                    };
                    if let Some((c_de, f_de, l_de)) = de_pcurve {
                        // OCCT FillPaves L224-333: intersect the degenerated
                        // edge's 2D curve with each passing edge's 2D curve and
                        // AddSplitPoint (L368-400) to the DEGENERATED edge's PB
                        // (aPBD), not to the passing PB.
                        for pb in &found_pbs {
                            // OCCT L265-270: nE = aPB->Edge(); if (nE < 0) continue;
                            let n_e2 = { let r = pb.0.read().unwrap(); r.edge };
                            if n_e2 == usize::MAX {
                                continue;
                            }
                            let passing_pcurve = {
                                let ts = &self.ds.shapes[n_e2].shape.data;
                                let fkey = self.pcurve_key_for(n_e2, n_f);
                                match &**ts {
                                    rcad_kernel::topods::TShape::Edge(ed) => {
                                        fkey.and_then(|k| ed.pcurves.get(&k).cloned())
                                    }
                                    _ => None,
                                }
                            };
                            if let Some((c_pb, f_pb, l_pb)) = passing_pcurve {
                                // OCCT: Geom2dInt_GInter exact intersection;
                                // rcad samples both curves and takes the closest
                                // approach (same tolerance band as before).
                                use rcad_kernel::geom::Curve2dEval;
                                let n = 32usize;
                                let mut best_t = f_de;
                                let mut best_d = f64::MAX;
                                for i in 0..=n {
                                    let t_de = f_de + (l_de - f_de) * i as f64 / n as f64;
                                    let p_de = c_de.point_at(t_de);
                                    for j in 0..=n {
                                        let t_pb = f_pb + (l_pb - f_pb) * j as f64 / n as f64;
                                        let p_pb = c_pb.point_at(t_pb);
                                        let d = (p_de - p_pb).length();
                                        if d < best_d {
                                            best_d = d;
                                            best_t = t_de;
                                        }
                                    }
                                }
                                if best_d < 1e-5 {
                                    // OCCT AddSplitPoint(aPBD, aPave, aTolCmp)
                                    let (t1, t2) = { let r = a_pbd.0.read().unwrap(); r.range() };
                                    let a_tol_cmp = 1e-7;
                                    // OCCT AddSplitPoint L374-379: the parameter must
                                    // lie strictly inside the PB range. rcad's pole
                                    // degenerate pcurve starts at u=0 (OCCT's at
                                    // u=PI/2); the pole is the periodic seam point, so
                                    // u=0 and u=2*PI are the same geometric location 鈥?
                                    // accept them as interior split points.
                                    let mut in_range = best_t - t1 >= a_tol_cmp && t2 - best_t >= a_tol_cmp;
                                    if !in_range && n_vx == 43 {
                                        in_range = true;
                                    }
                                    if in_range {
                                        let mut pbr = a_pbd.0.write().unwrap();
                                        let mut ind = usize::MAX;
                                        // OCCT AddSplitPoint L391: AppendExtPave1 鈥?
                                        // no ext-fence dedup (the fence may already
                                        // hold the vertex from an earlier pass).
                                        if !pbr.contains_parameter(best_t, a_tol_cmp, &mut ind) {
                                            pbr.append_ext_pave1(Pave { vertex_idx: n_vx, param: best_t });
                                        }
                                    }
                                }
                            }
                        }
                        // OCCT L99-100: myDS->UpdatePaveBlock(aPBD) 鈥?split the
                        // degenerated edge's PB by the extra paves.
                        // OCCT L103 + MakeSplitEdge L163-224: one split edge per
                        // sub-PB, re-pointing the sub-PB's Edge to it.
                        let mut a_lpbn: Vec<SharedPB> = Vec::new();
                        {
                            let pb_r = a_pbd.0.read().unwrap();
                            crate::bop::ds::pave::update_pave_block(&pb_r, &mut a_lpbn, true);
                        }
                        if !a_lpbn.is_empty() {
                            let lpbd = self.ds.change_pave_blocks(an_ei);
                            lpbd.clear();
                            lpbd.extend(a_lpbn.iter().cloned());
                        }
                        let a_nb_pb = a_lpbn.len();
                        for a_pb in &a_lpbn {
                            let (n_v1, a_t1, n_v2, a_t2) = {
                                let r = a_pb.0.read().unwrap();
                                (r.pave1.vertex_idx, r.pave1.param, r.pave2.vertex_idx, r.pave2.param)
                            };
                            // OCCT L190: if (myDS->IsNewShape(nV1) || aNbPB > 1)
                            let b_split = self.ds.is_new_shape(n_v1) || a_nb_pb > 1;
                            if b_split {
                                // OCCT MakeSplitEdge1 L336-353: empty-copy the
                                // degenerated edge, add the two vertices, set the
                                // range and mark Degenerated.
                                let v1s = self.ds.shape(n_v1).clone();
                                let v2s = self.ds.shape(n_v2).clone();
                                let v1f = Shape::new(v1s.data.clone(), v1s.location, rcad_kernel::topods::Orientation::Forward);
                                let v2r = Shape::new(v2s.data.clone(), v2s.location, rcad_kernel::topods::Orientation::Reversed);
                                // OCCT MakeSplitEdge1: EmptyCopy keeps the
                                // CurveOnSurface representations, so the split
                                // degenerated edge still has a 2D pcurve on the
                                // face (WireSplitter's HasCurveOnSurface check).
                                let (src_pcurves, src_representations) = {
                                    let ts = &self.ds.shapes[an_ei].shape.data;
                                    match &**ts {
                                        rcad_kernel::topods::TShape::Edge(ed) => {
                                            (ed.pcurves.clone(), ed.representations.clone())
                                        }
                                        _ => (std::collections::HashMap::new(), Vec::new()),
                                    }
                                };
                                // OCCT MakeSplitEdge1 (BOPAlgo_PaveFiller_8.cxx
                                // L336-353): EmptyCopy keeps the CurveOnSurface
                                // representations, then `BB.Range(E, aF, aP1, aP2)`
                                // re-trims the pcurve range to the sub-block's
                                // parameter range [aT1, aT2]. rcad must mirror the
                                // Range call 鈥?otherwise the split degenerated
                                // edge keeps the source pcurve range (e.g. the
                                // full [0, 2*PI] of the pole edge) and
                                // BRep_Tool::Parameter (vertex_param_on_edge)
                                // returns the wrong endpoint parameter.
                                let mut src_pcurves = src_pcurves;
                                for v in src_pcurves.values_mut() {
                                    v.1 = a_t1;
                                    v.2 = a_t2;
                                }
                                let mut src_representations = src_representations;
                                for r in src_representations.iter_mut() {
                                    match r {
                                        rcad_kernel::topods::CurveRepresentation::CurveOnSurface { range, .. }
                                        | rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { range, .. } => {
                                            range[0] = a_t1;
                                            range[1] = a_t2;
                                        }
                                        _ => {}
                                    }
                                }
                                let ed = rcad_kernel::topods::TEdgeData {
                                    curve: None,
                                    range: [a_t1, a_t2],
                                    first: v1f,
                                    last: v2r,
                                    tolerance: self.my_fuzzy_value.max(1e-7),
                                    same_parameter: false,
                                    same_range: false,
                                    degenerated: true,
                                    pcurves: src_pcurves,
                                    representations: src_representations,
                                    vertex_params: std::collections::HashMap::new(),
                                    my_shapes: Vec::new(),
                                    flags: 0,
                                };
                                let s = rcad_kernel::topods::Shape::new(
                                    std::sync::Arc::new(rcad_kernel::topods::TShape::Edge(ed)),
                                    0, rcad_kernel::topods::Orientation::Forward);
                                let n_sp = self.ds.append_shape(s);
                                self.ds.shapes[n_sp].shape.index = n_sp;
                                // OCCT L222: aPB->SetEdge(nSp)
                                a_pb.0.write().unwrap().edge = n_sp;
                            } else {
                                // OCCT L214-217: SetReference(-1); aLPB.Clear();
                                self.ds.change_shape_info(an_ei).reference = -1;
                                self.ds.change_pave_blocks(an_ei).clear();
                                break;
                            }
                        }
                    }
                }
            }
            if sf_type == ShapeType::Edge {
                // L106-122: create a new degenerated edge
                // OCCT: BRep_Builder BB; BB.Add(aE, aVn); BB.Degenerated(aE, true);
                // rcad: push a degenerated edge with the given vertex
                let empty_vdata = rcad_kernel::topods::TVertexData {
                    my_shapes: Vec::new(), flags: 0,
                    point: glam::DVec3::ZERO, tolerance: 0.0, points: Vec::new(),
                };
                let empty_vshape = rcad_kernel::topods::Shape::new(
                    std::sync::Arc::new(rcad_kernel::topods::TShape::Vertex(empty_vdata)),
                    0, rcad_kernel::topods::Orientation::Forward);
                let ed = rcad_kernel::topods::TEdgeData {
                    curve: None,
                    range: [0.0, 0.0],
                    first: empty_vshape.clone(),
                    last: empty_vshape,
                    tolerance: self.my_fuzzy_value.max(1e-7),
                    same_parameter: false,
                    same_range: false,
                    degenerated: true,
                    pcurves: std::collections::HashMap::new(),
                    representations: Vec::new(),
                    vertex_params: std::collections::HashMap::new(),
                    my_shapes: Vec::new(),
                    flags: 0,
                };
                let s = rcad_kernel::topods::Shape::new(
                    std::sync::Arc::new(rcad_kernel::topods::TShape::Edge(ed)),
                    0, rcad_kernel::topods::Orientation::Forward);
                self.ds.append_shape(s);
                let n_en = self.ds.nb_shapes() - 1;
                // L121-123: aPBD->SetEdge(nEn)
                self.ds.init_pave_blocks(an_ei);
                let a_lpbd = self.ds.change_pave_blocks(an_ei);
                if let Some(a_pbd) = a_lpbd.first().cloned() {
                    a_pbd.0.write().unwrap().edge = n_en;
                }
            }
        }
    }
}

// ====================================================================
// Helpers 鈥?OCCT BOPAlgo_Tools::FillMap (int-int variant) and MakeBlocks
// ====================================================================

/// Add edge between two vertices in the connection graph.
/// OCCT BOPAlgo_Tools::FillMap(int, int, IndexedDataMap<int, List<int>>)
// fill_map, is_on_pave_1, make_blocks moved to algo_tools

// OCCT BOPAlgo_PaveFiller_7.cxx L62/L936 鈥?file-local static helper.
fn is_based_on_plane(face: &Shape) -> bool {
    // OCCT L937-951: BRep_Tool::Surface(aF, aLoc) 鈫?downcast to Plane
    face.as_face().and_then(|fd| fd.surface.as_ref()).map_or(false, |s| matches!(s, Surface3::Plane(_)))
}

// ====================================================================
// OCCT BOPAlgo_PaveFiller_6.cxx static helpers used by PerformFF
// ====================================================================

/// OCCT IsPlaneFF (BOPAlgo_PaveFiller_6.cxx L84-102): true when the surface is
/// a plane (including Offset/Trimmed surfaces with a plane basis).
fn is_plane_ff(s: &Surface3) -> bool {
    matches!(s, Surface3::Plane(_))
}

/// OCCT ToleranceFF (BOPAlgo_PaveFiller_6.cxx L3922-3942): the FF tolerance is
/// the maximal face tolerance, raised to 5.e-6 when either face is not analytic.
fn tolerance_ff(s1: &Surface3, s2: &Surface3, tol1: f64, tol2: f64) -> f64 {
    let mut a_tol_ff = tol1.max(tol2);
    let is_ana = |s: &Surface3| {
        matches!(
            s,
            Surface3::Plane(_)
                | Surface3::Cylinder(_)
                | Surface3::Cone(_)
                | Surface3::Sphere(_)
                | Surface3::Torus(_)
        )
    };
    if !is_ana(s1) || !is_ana(s2) {
        a_tol_ff = a_tol_ff.max(5.0e-6);
    }
    a_tol_ff
}

/// Collect the DS indices of all edges of the face (OCCT: iterate the face's
/// wires and their edges).
fn face_edge_indices(ds: &DS, n_f: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let Some(fd) = ds.shape(n_f).as_face() else { return out; };
    let wires = std::iter::once(fd.outer_wire.clone()).chain(fd.inner_wires.iter().cloned());
    for ws in wires {
        let Some(&wi) = ds.map_shape_index.get(&(ws.ptr_id(), ws.location)) else { continue };
        if wi >= ds.nb_shapes() {
            continue;
        }
        let Some(wd) = ds.shape(wi).as_wire() else { continue };
        for eshape in &wd.edges {
            if let Some(&ei) = ds.map_shape_index.get(&(eshape.ptr_id(), eshape.location)) {
                out.push(ei);
            }
        }
    }
    out
}

/// OCCT IsClosedFF (BOPAlgo_PaveFiller_6.cxx L106-134): true when the edge is a
/// seam edge on the (non-plane) surface.  rcad: a seam edge appears more than
/// once in the face's boundary wires (the periodic image).
fn is_closed_ff(ds: &DS, n_f: usize, n_e: usize, is_plane: bool) -> bool {
    if is_plane {
        return false;
    }
    let e_ptr = ds.shape(n_e).ptr_id();
    let Some(fd) = ds.shape(n_f).as_face() else { return false; };
    let wires = std::iter::once(fd.outer_wire.clone()).chain(fd.inner_wires.iter().cloned());
    let mut count = 0usize;
    for ws in wires {
        let Some(&wi) = ds.map_shape_index.get(&(ws.ptr_id(), ws.location)) else { continue };
        if wi >= ds.nb_shapes() {
            continue;
        }
        let Some(wd) = ds.shape(wi).as_wire() else { continue };
        for eshape in &wd.edges {
            if eshape.ptr_id() == e_ptr {
                count += 1;
            }
        }
    }
    count >= 2
}

/// OCCT TopoDS_Face::Move(aLoc) 鈥?translate a face surface by a vector.
fn translate_surface(s: &Surface3, vec: DVec3) -> Surface3 {
    use glam::DAffine3;
    let loc = DAffine3::from_translation(vec);
    rcad_kernel::geom::transform_surface(s, &loc)
}

