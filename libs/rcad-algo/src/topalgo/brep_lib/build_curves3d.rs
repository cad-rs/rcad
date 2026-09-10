// OCCT BRepLib 3d-curve build family (TKTopAlgo/BRepLib/BRepLib.cxx) —
// 1:1 translation:
// - BRepLib::CheckSameRange (cxx L149-183)
// - BRepLib::SameRange (cxx L187-266)
// - static evaluateMaxSegment (cxx L273-297)
// - BRepLib::BuildCurve3d (cxx L301-456; BRepLib.hxx L90-94 defaults
//   Tolerance=1.0e-5, Continuity=GeomAbs_C1, MaxDegree=14, MaxSegment=0)
// - BRepLib::BuildCurves3d(S) (cxx L460-464)
// - BRepLib::BuildCurves3d(S, Tolerance, Continuity, MaxDegree, MaxSegment)
//   (cxx L468-489; BRepLib.hxx L99-103 defaults Continuity=GeomAbs_C1,
//   MaxDegree=14, MaxSegment=0)
//
// Architecture difference (the rcad pool model, topods::BRep): the OCCT
// BRep_Builder edits the shared edge TShape in place through the handle.
// rcad mirrors this through the BRep pool (BRep::edge_mut_inplace) — the
// functions below take `&mut BRep` where OCCT mutates, and `&BRep` where
// OCCT only reads.  TopLoc_Location maps to the u32 location index of the
// rcad Shape.

use std::collections::HashSet;

use rcad_kernel::geom::{BSplineCurve2, BSplineSurface, Curve2d, Curve3, Plane, Surface3};
use rcad_kernel::topo::topods::{
    surface_same, BRep, BRepTool, CurveRepresentation, GeomAbsShape, ShapeType, TShape,
};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::geomalgo::int_patch::{classify_surface_type, GeomAbsSurfaceType};

use super::brep_lib::BRepLib;

impl BRepLib {
    // OCCT BRepLib.cxx L149-183
    /// OCCT BRepLib::CheckSameRange(AnEdge, Tolerance) — true when every
    /// curve representation of the edge carries the same (first, last)
    /// range within the tolerance.
    pub fn check_same_range(the_brep: &BRep, an_edge: &Shape, tolerance: f64) -> bool {
        // OCCT L151: bool IsSameRange = true, first_time_in = true.
        let mut is_same_range = true;
        let mut first_time_in = true;
        // OCCT L156-157: double first, last; current_first = 0., current_last = 0.
        let mut first: f64;
        let mut last: f64;
        let mut current_first: f64 = 0.0;
        let mut current_last: f64 = 0.0;
        // OCCT L153-154: the iterator over the edge curve representations
        // (the rcad representations list in insertion order).
        let a_ed = the_brep.edge(an_edge.clone());
        for a_cr in a_ed.representations.iter() {
            // OCCT L160: while (IsSameRange && an_Iterator.More()).
            if !is_same_range {
                break;
            }
            // OCCT L162-164: down_cast<BRep_GCurve> — the null down-cast
            // (the BRep_CurveOn2Surfaces regularity) is skipped.
            if let Some((a_f, a_l)) = gcurve_range(the_brep, an_edge, a_cr) {
                // OCCT L166-167.
                first = a_f;
                last = a_l;
                if first_time_in {
                    // OCCT L169-172.
                    current_first = first;
                    current_last = last;
                    first_time_in = false;
                } else {
                    // OCCT L176-177.
                    is_same_range = ((current_first - first).abs() <= tolerance)
                        && ((current_last - last).abs() <= tolerance);
                }
            }
            // OCCT L180: an_Iterator.Next().
        }
        // OCCT L182.
        is_same_range
    }

    // OCCT BRepLib.cxx L187-266
    /// OCCT BRepLib::SameRange(AnEdge, Tolerance) — makes every pcurve of
    /// the edge span the same parameter range as the reference range (the
    /// 3d curve range when present, otherwise the first geometric
    /// representation range), by rebuilding each deviating pcurve over the
    /// reference range.
    pub fn same_range(the_brep: &mut BRep, an_edge: &Shape, tolerance: f64) {
        // OCCT L192: Curve2dPtr, Curve2dPtr2, NewCurve2dPtr, NewCurve2dPtr2.
        let mut curve2d_ptr: Option<Curve2d> = None;
        let mut curve2d_ptr2: Option<Curve2d> = None;
        let mut new_curve2d_ptr: Option<Curve2d>;
        let mut new_curve2d_ptr2: Option<Curve2d>;
        // OCCT L193: TopLoc_Location LocalLoc.
        let mut local_loc: u32 = 0;
        // OCCT L195: first_time_in = true, has_curve, has_closed_curve.
        let mut first_time_in = true;
        let mut has_curve = false;
        let mut has_closed_curve = false;
        // OCCT L197: first, current_first, last, current_last.
        let mut first: f64;
        let mut current_first: f64 = 0.0;
        let mut last: f64;
        let mut current_last: f64 = 0.0;
        // OCCT L199-203: C = BRep_Tool::Curve(AnEdge, LocalLoc,
        // current_first, current_last); a present 3d curve fixes the
        // reference range.
        if let Some((_, a_l, a_f, a_ll)) = brep_tool_curve(the_brep, an_edge) {
            local_loc = a_l;
            current_first = a_f;
            current_last = a_ll;
            first_time_in = false;
        }
        // OCCT L205: while (an_Iterator.More()).
        let a_len = the_brep.edge(an_edge.clone()).representations.len();
        let mut an_iterator = 0usize;
        while an_iterator < a_len {
            // OCCT L207-208: down_cast<BRep_GCurve>(an_Iterator.Value())
            // with the null check — the BRep_CurveOn2Surfaces regularity is
            // skipped.
            let a_cr = the_brep.edge(an_edge.clone()).representations[an_iterator].clone();
            if gcurve_is_geometric(&a_cr) {
                // OCCT L210: has_closed_curve = has_curve = false.
                has_closed_curve = false;
                has_curve = false;
                // OCCT L211-212: first = First(); last = Last().
                let (a_f, a_l) = match gcurve_range(the_brep, an_edge, &a_cr) {
                    Some((a_f, a_l)) => (a_f, a_l),
                    None => (0.0, 0.0),
                };
                first = a_f;
                last = a_l;
                // OCCT L213-217: if (IsCurveOnSurface()) — true for the
                // closed-surface representation too.
                if gcurve_is_curve_on_surface(&a_cr) {
                    curve2d_ptr = gcurve_pcurve(&a_cr);
                    has_curve = true;
                }
                // OCCT L218-222: if (IsCurveOnClosedSurface()).
                if gcurve_is_curve_on_closed_surface(&a_cr) {
                    curve2d_ptr2 = gcurve_pcurve2(&a_cr);
                    has_closed_curve = true;
                }
                // OCCT L223: if (has_curve || has_closed_curve).
                if has_curve || has_closed_curve {
                    if first_time_in {
                        // OCCT L227-229.
                        current_first = first;
                        current_last = last;
                        first_time_in = false;
                    }
                    // OCCT L232-233.
                    if ((first - current_first).abs() > rcad_kernel::precision::CONFUSION)
                        || ((last - current_last).abs() > rcad_kernel::precision::CONFUSION)
                    {
                        // OCCT L235-243: if (has_curve) — GeomLib::SameRange
                        // rebuilds the pcurve over the reference range and
                        // the representation takes it (the PCurve setter).
                        if has_curve {
                            new_curve2d_ptr = Some(crate::geomalgo::geom_lib_same_range::same_range(
                                tolerance,
                                curve2d_ptr.as_ref().unwrap(),
                                first,
                                last,
                                current_first,
                                current_last,
                            ));
                            gcurve_set_pcurve(
                                the_brep,
                                an_edge,
                                an_iterator,
                                new_curve2d_ptr.clone(),
                            );
                        }
                        // OCCT L246-255: if (has_closed_curve) — the same
                        // for the second pcurve (the PCurve2 setter).
                        if has_closed_curve {
                            new_curve2d_ptr2 = Some(crate::geomalgo::geom_lib_same_range::same_range(
                                tolerance,
                                curve2d_ptr2.as_ref().unwrap(),
                                first,
                                last,
                                current_first,
                                current_last,
                            ));
                            gcurve_set_pcurve2(
                                the_brep,
                                an_edge,
                                an_iterator,
                                new_curve2d_ptr2.clone(),
                            );
                        }
                    }
                }
            }
            // OCCT L260: an_Iterator.Next().
            an_iterator += 1;
        }
        let _ = local_loc;
        // OCCT L262-263: B.Range(TopoDS::Edge(AnEdge), current_first, current_last).
        builder_range(the_brep, an_edge, current_first, current_last);
        // OCCT L265: B.SameRange(AnEdge, true).
        builder_same_range(the_brep, an_edge, true);
    }

    // OCCT BRepLib.cxx L301-456; BRepLib.hxx L90-94 (defaults
    // Tolerance=1.0e-5, Continuity=GeomAbs_C1, MaxDegree=14, MaxSegment=0)
    /// OCCT BRepLib::BuildCurve3d(AnEdge, Tolerance, Continuity, MaxDegree,
    /// MaxSegment) — computes the 3d curve of the edge from its pcurves.
    /// Returns true when the curve existed or was computed; false when
    /// there is no planar pcurve or the computation failed.
    pub fn build_curve3d(
        the_brep: &mut BRep,
        an_edge: &Shape,
        tolerance: f64,
        continuity: GeomAbsShape,
        max_degree: i32,
        max_segment: i32,
    ) -> bool {
        // OCCT L307-311: int ii, jj.
        let mut ii: i32;
        let mut jj: usize = 0;
        // OCCT L313: TopLoc_Location LocalLoc, L[2], LC.
        let mut local_loc: u32 = 0;
        let mut l: [u32; 2] = [0, 0];
        let mut lc: u32 = 0;
        // OCCT L314: double f, l, fc, lc, first[2], last[2], tolerance,
        // max_deviation, average_deviation (the rcad name l_par stands for
        // the OCCT local `l` — `l` collides with the location array).
        let mut f: f64 = 0.0;
        let mut l_par: f64 = 0.0;
        let mut fc: f64 = 0.0;
        let mut lc_par: f64 = 0.0;
        let mut first: [f64; 2] = [0.0, 0.0];
        let mut last: [f64; 2] = [0.0, 0.0];
        let mut tolerance_local: f64 = 0.0;
        let mut max_deviation: f64 = 0.0;
        let mut average_deviation: f64 = 0.0;
        // OCCT L315: Curve2dPtr, Curve2dArray[2].
        let mut curve2d_ptr: Option<Curve2d> = None;
        let mut curve2d_array: [Option<Curve2d>; 2] = [None, None];
        // OCCT L316: SurfacePtr, SurfaceArray[2].
        let mut surface_ptr: Option<Surface3> = None;
        let mut surface_array: [Option<Surface3>; 2] = [None, None];

        // OCCT L318-325: if the edge has a 3d curve returns true.
        if let Some((a_c, a_l, a_f, a_ll)) = brep_tool_curve(the_brep, an_edge) {
            let _ = a_c;
            local_loc = a_l;
            f = a_f;
            l_par = a_ll;
            return true;
        }
        // OCCT L327-333: this should not exist but UpdateEdge makes funny
        // things if the edge is not same range.
        if !Self::check_same_range(the_brep, an_edge, rcad_kernel::precision::CONFUSION) {
            Self::same_range(the_brep, an_edge, tolerance);
        }

        // OCCT L335-340: search a curve on a plane.
        let mut s: Option<Surface3> = None;
        let mut pc: Option<Curve2d> = None;
        let mut i: i32 = 0;
        let mut p: Option<Plane> = None;
        let mut not_done = true;

        while not_done {
            // OCCT L344: i++.
            i += 1;
            // OCCT L345: BRep_Tool::CurveOnSurface(AnEdge, PC, S, LocalLoc, f, l, i).
            let a_cos = brep_tool_curve_on_surface_index(the_brep, an_edge, i);
            pc = a_cos.pc.clone();
            s = a_cos.s.clone();
            local_loc = a_cos.l;
            f = a_cos.first;
            l_par = a_cos.last;
            // OCCT L346-355: RT = down_cast<Geom_RectangularTrimmedSurface>(S);
            // P = RT.IsNull() ? down_cast<Geom_Plane>(S)
            //                : down_cast<Geom_Plane>(RT->BasisSurface()).
            p = match &s {
                Some(Surface3::Trimmed(a_rt)) => match a_rt.basis.as_ref() {
                    Surface3::Plane(a_pl) => Some(*a_pl),
                    _ => None,
                },
                Some(Surface3::Plane(a_pl)) => Some(*a_pl),
                _ => None,
            };
            // OCCT L356: not_done = !S.IsNull() && P.IsNull().
            not_done = s.is_some() && p.is_none();
        }
        // OCCT L358: if (!P.IsNull()).
        if let Some(a_pl) = p {
            // OCCT L361: gp_Ax2 axes = P->Position().Ax2().
            let axes = plane_position_ax2(&a_pl);
            // OCCT L362: C3d = GeomLib::To3d(axes, PC) — GAP carrier (the
            // TKGeomBase/GeomLib To3d body is not translated; the OCCT
            // failure path is preserved: the null handle makes
            // BuildCurve3d return false).
            let a_c3d_opt = geom_lib_to_3d(&axes, pc.as_ref().unwrap());
            // OCCT L363-366: if (C3d.IsNull()) return false.
            let a_c3d = match a_c3d_opt {
                Some(a_c3d) => a_c3d,
                None => return false,
            };
            // OCCT L368: double First, Last.
            let mut first_out: f64 = 0.0;
            let mut last_out: f64 = 0.0;
            // OCCT L370-371: BRep_Builder B; B.UpdateEdge(AnEdge, C3d, LocalLoc, 0.0e0).
            builder_update_edge_curve3d(the_brep, an_edge, &a_c3d, 0.0e0);
            // OCCT L372: BRep_Tool::Range(AnEdge, S, LC, First, Last).
            let (a_lc, a_f_out, a_l_out) =
                brep_tool_range_on_surface(the_brep, an_edge, s.as_ref().unwrap());
            lc = a_lc;
            first_out = a_f_out;
            last_out = a_l_out;
            // OCCT L373: B.Range(AnEdge, First, Last) — do not forget 3D
            // range (PRO6412).
            builder_range(the_brep, an_edge, first_out, last_out);
        } else {
            // OCCT L377-381: compute the 3d curve using existing surface.
            fc = f;
            lc_par = l_par;
            // OCCT L382: if (!BRep_Tool::Degenerated(AnEdge)).
            if !brep_tool_degenerated(the_brep, an_edge) {
                // OCCT L384: jj = 0.
                jj = 0;
                // OCCT L385-405: for (ii = 1; ii < 3; ii++) — collect the
                // first two pcurves of the edge.
                for a_ii in 1..3 {
                    ii = a_ii;
                    let a_cos = brep_tool_curve_on_surface_index(the_brep, an_edge, ii);
                    curve2d_ptr = a_cos.pc.clone();
                    surface_ptr = a_cos.s.clone();
                    local_loc = a_cos.l;
                    fc = a_cos.first;
                    lc_par = a_cos.last;
                    // OCCT L396-404.
                    if curve2d_ptr.is_some() && jj < 2 {
                        curve2d_array[jj] = curve2d_ptr.clone();
                        surface_array[jj] = surface_ptr.clone();
                        l[jj] = local_loc;
                        first[jj] = fc;
                        last[jj] = lc_par;
                        jj += 1;
                    }
                }
                // OCCT L406-409.
                f = first[0];
                l_par = last[0];
                curve2d_ptr = curve2d_array[0].clone();
                surface_ptr = surface_array[0].clone();

                // OCCT L411-417: the pcurve/surface adaptors and the
                // Adaptor3d_CurveOnSurface over their handles (the rcad
                // value clones stand for the OCCT handle constructions).
                let an_adaptor3d_curve2d = Geom2dAdaptorCurve::new(curve2d_ptr.clone(), f, l_par);
                let an_adaptor3d_surface = GeomAdaptorSurface::new(surface_ptr.clone());
                let an_adaptor3d_curve2d_ptr = an_adaptor3d_curve2d.clone();
                let an_adaptor3d_surface_ptr = an_adaptor3d_surface.clone();
                let curve_on_surface = Adaptor3dCurveOnSurface::new(
                    &an_adaptor3d_curve2d_ptr,
                    &an_adaptor3d_surface_ptr,
                );

                // OCCT L419: handle<Geom_Curve> NewCurvePtr.
                let mut new_curve_ptr: Option<Curve3> = None;

                // OCCT L421-430: GeomLib::BuildCurve3d — GAP carrier (the
                // AdvApprox approximation of GeomLib.cxx is not translated;
                // the OCCT failure path is preserved: the null NewCurvePtr
                // makes BuildCurve3d return false).
                geom_lib_build_curve3d(
                    tolerance,
                    &curve_on_surface,
                    f,
                    l_par,
                    &mut new_curve_ptr,
                    &mut max_deviation,
                    &mut average_deviation,
                    continuity,
                    max_degree,
                    evaluate_max_segment(max_segment, &curve_on_surface),
                );
                // OCCT L431: BRep_Builder B (the pool is the builder arena).
                // OCCT L432: tolerance = BRep_Tool::Tolerance(AnEdge).
                tolerance_local = brep_tool_tolerance(the_brep, an_edge);
                // OCCT L433-435: max_deviation = std::max(tolerance, Tolerance).
                max_deviation = tolerance_local.max(tolerance);
                // OCCT L436-439: if (NewCurvePtr.IsNull()) return false.
                let a_nc = match new_curve_ptr {
                    Some(a_nc) => a_nc,
                    None => return false,
                };
                // OCCT L440: B.UpdateEdge(TopoDS::Edge(AnEdge), NewCurvePtr, L[0], max_deviation).
                builder_update_edge_curve3d(the_brep, an_edge, &a_nc, max_deviation);
                let _ = l[0];
                // OCCT L441-448: if (jj == 1) — if there is only one curve
                // on surface attached to the edge than it can be qualified
                // sameparameter.
                if jj == 1 {
                    builder_same_parameter(the_brep, an_edge, true);
                }
            } else {
                // OCCT L450-453.
                return false;
            }
        }
        // OCCT L455.
        true
    }

    // OCCT BRepLib.cxx L460-464
    /// OCCT BRepLib::BuildCurves3d(S) — the tolerance form with the OCCT
    /// default tolerance 1.0e-5.
    pub fn build_curves3d(the_brep: &mut BRep, the_s: &Shape) -> bool {
        Self::build_curves3d_tol(the_brep, the_s, 1.0e-5)
    }

    // OCCT BRepLib.cxx L468-489; BRepLib.hxx L99-103 (defaults
    // Continuity=GeomAbs_C1, MaxDegree=14, MaxSegment=0)
    /// OCCT BRepLib::BuildCurves3d(S, Tolerance, Continuity, MaxDegree,
    /// MaxSegment) — computes the 3d curves for all the edges of S.
    /// Returns false if one of the computations failed.
    pub fn build_curves3d_tol(the_brep: &mut BRep, the_s: &Shape, the_tolerance: f64) -> bool {
        // The OCCT hxx defaults for the trailing parameters.
        let the_continuity = GeomAbsShape::C1;
        let the_max_degree: i32 = 14;
        let the_max_segment: i32 = 0;
        // OCCT L474: boolean_value, ok = true.
        let mut boolean_value: bool;
        let mut ok = true;
        // OCCT L475: NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>
        // a_counter — the shape set keyed by (TShape, Location).
        let mut a_counter: HashSet<(u64, u32)> = HashSet::new();
        // OCCT L476: TopExp_Explorer ex(S, TopAbs_EDGE).
        let a_ex = bat::explorer(the_s, ShapeType::Edge, ShapeType::Shape);
        // OCCT L478: while (ex.More()).
        for an_edge in a_ex {
            // OCCT L480: if (a_counter.Add(ex.Current())).
            if a_counter.insert((an_edge.ptr_id(), an_edge.location)) {
                // OCCT L482-484: boolean_value = BuildCurve3d(
                //   TopoDS::Edge(ex.Current()), Tolerance, Continuity,
                //   MaxDegree, MaxSegment).
                boolean_value = Self::build_curve3d(
                    the_brep,
                    &an_edge,
                    the_tolerance,
                    the_continuity,
                    the_max_degree,
                    the_max_segment,
                );
                ok = ok && boolean_value;
            }
            // OCCT L486: ex.Next().
        }
        // OCCT L488.
        ok
    }
}

// ---------------------------------------------------------------------------
// OCCT file statics.
// ---------------------------------------------------------------------------

// OCCT BRepLib.cxx L273-297
/// OCCT static evaluateMaxSegment(aMaxSegment, aCurveOnSurface) — returns
/// MaxSegment to pass in approximation, if MaxSegment==0 provided: 30 plus
/// the largest knot count of the B-spline surface/curve under the adaptor.
fn evaluate_max_segment(a_max_segment: i32, a_curve_on_surface: &Adaptor3dCurveOnSurface) -> i32 {
    // OCCT L276-279.
    if a_max_segment != 0 {
        return a_max_segment;
    }
    // OCCT L281-282: the surface and 2d-curve adaptors under the
    // Adaptor3d_CurveOnSurface.
    let a_surf = a_curve_on_surface.get_surface();
    let a_curv2d = a_curve_on_surface.get_curve();
    // OCCT L284: double aNbSKnots = 0, aNbC2dKnots = 0.
    let mut a_nb_s_knots: f64 = 0.0;
    let mut a_nb_c2d_knots: f64 = 0.0;
    // OCCT L286-290: BSpline surface -> max(NbUKnots, NbVKnots).
    if a_surf.get_type() == GeomAbsSurfaceType::BSplineSurface {
        if let Some(a_bspline) = a_surf.bspline() {
            // The rcad knot vectors carry the multiplicities expanded; the
            // OCCT NbUKnots()/NbVKnots() count the distinct knots.
            a_nb_s_knots =
                distinct_knots(&a_bspline.knots_u).max(distinct_knots(&a_bspline.knots_v)) as f64;
        }
    }
    // OCCT L291-294: BSpline 2d curve -> NbKnots.
    if a_curv2d.get_type() == GeomAbsCurveType::BSplineCurve {
        if let Some(a_bspline) = a_curv2d.bspline() {
            a_nb_c2d_knots = distinct_knots(&a_bspline.knots) as f64;
        }
    }
    // OCCT L295-296: (int)(30 + max(aNbSKnots, aNbC2dKnots)).
    (30.0 + a_nb_s_knots.max(a_nb_c2d_knots)) as i32
}

/// The number of distinct (sorted) knot values of an expanded rcad knot
/// vector — the OCCT Geom_BSplineKnots NbKnots() counting (arch.
/// difference: the rcad knots store the multiplicities expanded).
fn distinct_knots(knots: &[f64]) -> usize {
    let mut a_nb = 0usize;
    let mut a_prev = f64::NAN;
    for a_k in knots {
        if *a_k != a_prev {
            a_nb += 1;
            a_prev = *a_k;
        }
    }
    a_nb
}

// ---------------------------------------------------------------------------
// Adaptor re-hosts (TKG3d).  The OCCT adaptors wrap handles and carry
// First/Last; the rcad re-hosts wrap the value geometry (None plays the
// null handle) and read the type/knots from the variants.
// ---------------------------------------------------------------------------

/// OCCT GeomAbs_CurveType (TKG3d/GeomAbs/GeomAbs_CurveType.hxx L23-31).
/// The canonical nine-variant enum lives in rcad_kernel::math; the former
/// local copy was deleted (Rule 4).
use rcad_kernel::math::GeomAbsCurveType;

/// OCCT Geom2dAdaptor_Curve(Curve2d, First, Last) (TKG3d/Geom2dAdaptor) —
/// the pcurve adaptor (a null handle plays None; the trimmed curve is
/// flattened to its basis at the Load, Geom2dAdaptor.cxx).
#[derive(Clone)]
struct Geom2dAdaptorCurve {
    my_curve: Option<Curve2d>,
    my_first: f64,
    my_last: f64,
}

impl Geom2dAdaptorCurve {
    // OCCT Geom2dAdaptor_Curve::Load(C, U1, U2).
    fn new(the_curve: Option<Curve2d>, the_first: f64, the_last: f64) -> Self {
        Geom2dAdaptorCurve {
            my_curve: the_curve,
            my_first: the_first,
            my_last: the_last,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::GetType() — the basis-curve type (the
    /// trimmed curve is stripped at the Load).
    fn get_type(&self) -> GeomAbsCurveType {
        match &self.my_curve {
            Some(Curve2d::Line(_)) => GeomAbsCurveType::Line,
            Some(Curve2d::Circle(_)) => GeomAbsCurveType::Circle,
            Some(Curve2d::Ellipse(_)) => GeomAbsCurveType::Ellipse,
            Some(Curve2d::Hyperbola(_)) => GeomAbsCurveType::Hyperbola,
            Some(Curve2d::Parabola(_)) => GeomAbsCurveType::Parabola,
            Some(Curve2d::Bezier(_)) => GeomAbsCurveType::BezierCurve,
            Some(Curve2d::BSpline(_)) => GeomAbsCurveType::BSplineCurve,
            Some(Curve2d::Offset(_)) => GeomAbsCurveType::OffsetCurve,
            Some(Curve2d::Trimmed(a_tc)) => match a_tc.curve.as_ref() {
                Curve2d::BSpline(_) => GeomAbsCurveType::BSplineCurve,
                _ => GeomAbsCurveType::OtherCurve,
            },
            _ => GeomAbsCurveType::OtherCurve,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::BSpline() — the B-spline handle (the
    /// trimmed curve is stripped to its basis).
    fn bspline(&self) -> Option<&BSplineCurve2> {
        match &self.my_curve {
            Some(Curve2d::BSpline(a_bs)) => Some(a_bs),
            Some(Curve2d::Trimmed(a_tc)) => match a_tc.curve.as_ref() {
                Curve2d::BSpline(a_bs) => Some(a_bs),
                _ => None,
            },
            _ => None,
        }
    }
}

/// OCCT GeomAdaptor_Surface(Surface) (TKG3d/GeomAdaptor) — the surface
/// adaptor (a null handle plays None; the rectangular trimmed surface is
/// stripped to its basis at the Load, GeomAdaptor_Surface.cxx L423-425 —
/// the same rule as geom_adaptor_surface_get_type of brep_lib.rs).
#[derive(Clone)]
struct GeomAdaptorSurface {
    my_surface: Option<Surface3>,
}

impl GeomAdaptorSurface {
    // OCCT GeomAdaptor_Surface::Load(S).
    fn new(the_surface: Option<Surface3>) -> Self {
        GeomAdaptorSurface {
            my_surface: the_surface,
        }
    }

    /// OCCT GeomAdaptor_Surface::GetType().
    fn get_type(&self) -> GeomAbsSurfaceType {
        match &self.my_surface {
            Some(Surface3::Trimmed(a_ts)) => classify_surface_type(&a_ts.basis),
            Some(a_s) => classify_surface_type(a_s),
            None => GeomAbsSurfaceType::OtherSurface,
        }
    }

    /// OCCT Adaptor3d_Surface::BSpline() — the B-spline handle (the
    /// trimmed surface is stripped to its basis).
    fn bspline(&self) -> Option<&BSplineSurface> {
        match &self.my_surface {
            Some(Surface3::BSpline(a_bs)) => Some(a_bs),
            Some(Surface3::Trimmed(a_ts)) => match a_ts.basis.as_ref() {
                Surface3::BSpline(a_bs) => Some(a_bs),
                _ => None,
            },
            _ => None,
        }
    }
}

/// OCCT Adaptor3d_CurveOnSurface(Curve2d, Surf) (TKG3d) — the pcurve-on-
/// surface adaptor over the two handles.
#[derive(Clone)]
struct Adaptor3dCurveOnSurface {
    my_curve2d: Geom2dAdaptorCurve,
    my_surface: GeomAdaptorSurface,
}

impl Adaptor3dCurveOnSurface {
    // OCCT Adaptor3d_CurveOnSurface::Load(Curve2d, Surf).
    fn new(the_curve2d: &Geom2dAdaptorCurve, the_surface: &GeomAdaptorSurface) -> Self {
        Adaptor3dCurveOnSurface {
            my_curve2d: the_curve2d.clone(),
            my_surface: the_surface.clone(),
        }
    }

    /// OCCT Adaptor3d_CurveOnSurface::GetCurve().
    fn get_curve(&self) -> &Geom2dAdaptorCurve {
        &self.my_curve2d
    }

    /// OCCT Adaptor3d_CurveOnSurface::GetSurface().
    fn get_surface(&self) -> &GeomAdaptorSurface {
        &self.my_surface
    }
}

/// OCCT gp_Ax2 (TKMath/gp) — the right-handed axis system (location + main
/// direction + X direction) produced by gp_Ax3::Ax2().
#[derive(Debug, Clone, Copy)]
struct GpAx2 {
    /// OCCT gp_Ax2::Location().
    location: glam::DVec3,
    /// OCCT gp_Ax2::Direction() (the main direction, the plane normal).
    direction: glam::DVec3,
    /// OCCT gp_Ax2::XDirection().
    x_direction: glam::DVec3,
}

/// OCCT Geom_Plane::Position().Ax2() — the plane placement as an axis
/// system (gp_Ax3::Ax2() keeps Location/Direction/XDirection and drops the
/// Y direction; the rcad Plane carries u_dir = the Ax3 XDirection and
/// normal = the Ax3 Direction).
fn plane_position_ax2(the_plane: &Plane) -> GpAx2 {
    GpAx2 {
        location: the_plane.origin,
        direction: the_plane.normal,
        x_direction: the_plane.u_dir,
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (TKBRep/BRep/BRep_Tool.cxx).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Curve(E, L, f, l) (BRep_Tool.cxx L410-452) — the 3d
/// curve with the location transformation applied (the rcad
/// BRepTool::edge_curve_world re-host) and its range; the location out is
/// the edge's own (the rcad curve slot carries no location — arch.
/// difference).
fn brep_tool_curve(the_brep: &BRep, the_e: &Shape) -> Option<(Curve3, u32, f64, f64)> {
    let (a_c, a_range) = the_brep.edge_curve_world(the_e)?;
    Some((a_c, the_e.location, a_range[0], a_range[1]))
}

/// The out-parameter group of the indexed OCCT
/// BRep_Tool::CurveOnSurface(E, C, S, L, f, l, Index).
struct CurveOnSurfaceData {
    /// OCCT: the pcurve handle C.
    pc: Option<Curve2d>,
    /// OCCT: the support surface handle S.
    s: Option<Surface3>,
    /// OCCT: the location L (the rcad u32 index).
    l: u32,
    /// OCCT: First.
    first: f64,
    /// OCCT: Last.
    last: f64,
}

/// OCCT BRep_Tool::CurveOnSurface(E, C, S, L, f, l, Index)
/// (BRep_Tool.cxx L476-534) — the indexed pcurve: the representation list
/// is walked counting one slot per pcurve (two for a closed-surface
/// representation); on a miss every output stays null/zero.
fn brep_tool_curve_on_surface_index(
    the_brep: &BRep,
    the_e: &Shape,
    the_index: i32,
) -> CurveOnSurfaceData {
    // OCCT L478-479: the out parameters (C/S null, L identity, f = l = 0).
    let mut an_out = CurveOnSurfaceData {
        pc: None,
        s: None,
        l: 0,
        first: 0.0,
        last: 0.0,
    };
    // OCCT L480-482: if (Index < 1) return.
    if the_index < 1 {
        return an_out;
    }
    // OCCT L484: int i = 0.
    let mut a_i: i32 = 0;
    // OCCT L486-488: the representation list walk.
    let a_ed = the_brep.edge(the_e.clone());
    for a_cr in a_ed.representations.iter() {
        if gcurve_is_curve_on_surface(a_cr) {
            // OCCT L492: ++i.
            a_i += 1;
            // OCCT L494-500: the index comparison counting two pcurves for
            // the closed-surface representation.
            if a_i == the_index {
                an_out.pc = gcurve_pcurve(a_cr);
            } else if gcurve_is_curve_on_closed_surface(a_cr) {
                a_i += 1;
                if a_i == the_index {
                    an_out.pc = gcurve_pcurve2(a_cr);
                } else {
                    continue;
                }
            } else {
                continue;
            }
            // OCCT L502-505: S = GC->Surface(); L = E.Location() *
            // GC->Location(); GC->Range(First, Last); return.
            an_out.s = gcurve_surface(the_brep, a_cr);
            an_out.l = the_e.location;
            if let Some((a_f, a_l)) = gcurve_range(the_brep, the_e, a_cr) {
                an_out.first = a_f;
                an_out.last = a_l;
            }
            return an_out;
        }
    }
    // OCCT L512-516: not found — nulls and zeros.
    an_out.pc = None;
    an_out.s = None;
    an_out.l = 0;
    an_out.first = 0.0;
    an_out.last = 0.0;
    an_out
}

/// OCCT BRep_Tool::Range(E, S, L, First, Last) (BRep_Tool.cxx L964-992) —
/// the range of the pcurve on the surface; falls back to the 3d range when
/// no representation matches.  Architecture difference: the OCCT match on
/// (surface handle, L.Predivided(E.Location())) maps to the surface value
/// (surface_same) — the rcad representation stores the owning face key and
/// the surface is resolved from the face TShape; a surface carries one
/// representation except for the seam pair stored inside one
/// closed-surface representation, so the surface value identifies it.
fn brep_tool_range_on_surface(the_brep: &BRep, the_e: &Shape, the_s: &Surface3) -> (u32, f64, f64) {
    let a_ed = the_brep.edge(the_e.clone());
    for a_cr in a_ed.representations.iter() {
        if gcurve_is_curve_on_surface(a_cr) {
            if let Some(a_rep_surf) = gcurve_surface(the_brep, a_cr) {
                if surface_same(&a_rep_surf, the_s) {
                    if let Some((a_f, a_l)) = gcurve_range(the_brep, the_e, a_cr) {
                        return (the_e.location, a_f, a_l);
                    }
                }
            }
        }
    }
    // OCCT L987-990: if (!itcr.More()) Range(E, First, Last).
    let (a_f, a_l) = brep_tool_range_3d(the_brep, the_e);
    (the_e.location, a_f, a_l)
}

/// OCCT BRep_Tool::Range(E, First, Last) (BRep_Tool.cxx L931-960) — the
/// range of the first geometric representation (the non-null 3d curve, or
/// the first curve on surface); zeros when the edge carries none.
fn brep_tool_range_3d(the_brep: &BRep, the_e: &Shape) -> (f64, f64) {
    let a_ed = the_brep.edge(the_e.clone());
    for a_cr in a_ed.representations.iter() {
        match a_cr {
            // OCCT L940-949: IsCurve3D with a non-null curve — the rcad
            // Curve3D representation carries no range slot, the edge 3d
            // range is the OCCT BRep_Curve3D First/Last value (arch.
            // difference).
            CurveRepresentation::Curve3D { .. } => {
                return (a_ed.range[0], a_ed.range[1]);
            }
            // OCCT L950-956: IsCurveOnSurface.
            _ if gcurve_is_curve_on_surface(a_cr) => {
                if let Some((a_f, a_l)) = gcurve_range(the_brep, the_e, a_cr) {
                    return (a_f, a_l);
                }
            }
            _ => {}
        }
    }
    // OCCT L959: First = Last = 0.
    (0.0, 0.0)
}

/// OCCT BRep_Tool::Degenerated(E) — the degenerated flag of the edge.
fn brep_tool_degenerated(the_brep: &BRep, the_e: &Shape) -> bool {
    the_brep.edge(the_e.clone()).degenerated
}

/// OCCT BRep_Tool::Tolerance(E) (BRep_Tool.cxx L886-898) — the edge
/// tolerance clamped at Precision::Confusion.
fn brep_tool_tolerance(the_brep: &BRep, the_e: &Shape) -> f64 {
    let a_p = the_brep.edge(the_e.clone()).tolerance;
    if a_p > rcad_kernel::precision::CONFUSION {
        a_p
    } else {
        rcad_kernel::precision::CONFUSION
    }
}

// ---------------------------------------------------------------------------
// BRep_GCurve representation accessors — the OCCT down_cast<BRep_GCurve>
// reads over the curve representation list.  The rcad representations are
// an enum: Curve3D / CurveOnSurface / CurveOnClosedSurface are the
// BRep_GCurve kinds (First/Last carriers), CurveOn2Surfaces (the
// regularity) is the null down-cast.
// ---------------------------------------------------------------------------

/// OCCT BRep_CurveRepresentation::IsCurveOnSurface() — true for the
/// curve-on-surface and closed-surface kinds.
fn gcurve_is_curve_on_surface(a_cr: &CurveRepresentation) -> bool {
    matches!(
        a_cr,
        CurveRepresentation::CurveOnSurface { .. } | CurveRepresentation::CurveOnClosedSurface { .. }
    )
}

/// OCCT BRep_CurveOnSurface::IsCurveOnClosedSurface().
fn gcurve_is_curve_on_closed_surface(a_cr: &CurveRepresentation) -> bool {
    matches!(a_cr, CurveRepresentation::CurveOnClosedSurface { .. })
}

/// OCCT down_cast<BRep_GCurve> null check — a representation carrying
/// First/Last (all but the BRep_CurveOn2Surfaces regularity).
fn gcurve_is_geometric(a_cr: &CurveRepresentation) -> bool {
    !matches!(a_cr, CurveRepresentation::CurveOn2Surfaces { .. })
}

/// OCCT BRep_GCurve::First()/Last() — the representation range; None for
/// the null down-cast.  The rcad Curve3D representation carries no range
/// slot — the edge 3d range is the OCCT BRep_Curve3D First/Last value
/// (arch. difference).
fn gcurve_range(the_brep: &BRep, the_e: &Shape, a_cr: &CurveRepresentation) -> Option<(f64, f64)> {
    match a_cr {
        CurveRepresentation::Curve3D { .. } => {
            let a_ed = the_brep.edge(the_e.clone());
            Some((a_ed.range[0], a_ed.range[1]))
        }
        CurveRepresentation::CurveOnSurface { range, .. } => Some((range[0], range[1])),
        CurveRepresentation::CurveOnClosedSurface { range, .. } => Some((range[0], range[1])),
        CurveRepresentation::CurveOn2Surfaces { .. } => None,
    }
}

/// OCCT BRep_CurveOnSurface::PCurve() — the first pcurve.
fn gcurve_pcurve(a_cr: &CurveRepresentation) -> Option<Curve2d> {
    match a_cr {
        CurveRepresentation::CurveOnSurface { pcurve, .. } => Some(pcurve.clone()),
        CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => Some(pcurve1.clone()),
        _ => None,
    }
}

/// OCCT BRep_CurveOnClosedSurface::PCurve2().
fn gcurve_pcurve2(a_cr: &CurveRepresentation) -> Option<Curve2d> {
    match a_cr {
        CurveRepresentation::CurveOnClosedSurface { pcurve2, .. } => Some(pcurve2.clone()),
        _ => None,
    }
}

/// OCCT BRep_GCurve::Surface() — the support surface.  Architecture
/// difference: the rcad representation stores the owning face key; the
/// surface is resolved from the face TShape payload (the same resolution
/// as the kernel BRepTool::curve_on_surface fallback, topods.rs
/// L2385-2410).
fn gcurve_surface(the_brep: &BRep, a_cr: &CurveRepresentation) -> Option<Surface3> {
    let a_face_key = match a_cr {
        CurveRepresentation::CurveOnSurface { face, .. } => *face,
        CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
        _ => return None,
    };
    let a_ts = the_brep
        .tshapes
        .iter()
        .find(|a_ts| std::sync::Arc::as_ptr(a_ts) as u64 == a_face_key.0)?;
    match a_ts.as_ref() {
        TShape::Face(a_fd) => a_fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_CurveOnSurface::PCurve(C) setter — replaces the first pcurve
/// of the representation at the iterator slot `the_index` (and the
/// pcurves-map entry under its face key).
fn gcurve_set_pcurve(the_brep: &mut BRep, the_e: &Shape, the_index: usize, the_pc: Option<Curve2d>) {
    let the_pc = match the_pc {
        Some(the_pc) => the_pc,
        None => return,
    };
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    let a_face_key = match &a_ed.representations[the_index] {
        CurveRepresentation::CurveOnSurface { face, .. } => *face,
        CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
        _ => return,
    };
    match &mut a_ed.representations[the_index] {
        CurveRepresentation::CurveOnSurface { pcurve, .. } => *pcurve = the_pc.clone(),
        CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => *pcurve1 = the_pc.clone(),
        _ => {}
    }
    if let Some(a_entry) = a_ed.pcurves.get_mut(&a_face_key) {
        a_entry.0 = the_pc;
    }
}

/// OCCT BRep_CurveOnClosedSurface::PCurve2(C) setter — replaces the second
/// pcurve of the representation at the iterator slot `the_index`.
fn gcurve_set_pcurve2(the_brep: &mut BRep, the_e: &Shape, the_index: usize, the_pc: Option<Curve2d>) {
    let the_pc = match the_pc {
        Some(the_pc) => the_pc,
        None => return,
    };
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    match &mut a_ed.representations[the_index] {
        CurveRepresentation::CurveOnClosedSurface { pcurve2, .. } => *pcurve2 = the_pc,
        _ => {}
    }
}

// ---------------------------------------------------------------------------
// BRep_Builder re-hosts (TKBRep/BRep/BRep_Builder.cxx) — in-place pool
// edits (the OCCT builder edits the shared TShape through the handle).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::UpdateEdge(E, C3d, L, Tol) (BRep_Builder.cxx
/// UpdateEdge(E, C, L, Tol) 3d-curve form) — attaches the 3d curve to the
/// edge and raises the tolerance.  Architecture difference: the OCCT
/// BRep_Curve3D representation + location maps to the rcad curve slot of
/// the edge (locations travel on the Shape).
fn builder_update_edge_curve3d(the_brep: &mut BRep, the_e: &Shape, the_c3d: &Curve3, the_tol: f64) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.curve = Some(the_c3d.clone());
    a_ed.tolerance = a_ed.tolerance.max(the_tol);
}

/// OCCT BRep_Builder::Range(E, First, Last) — the 3d range of the edge.
fn builder_range(the_brep: &mut BRep, the_e: &Shape, the_first: f64, the_last: f64) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.range = [the_first, the_last];
}

/// OCCT BRep_Builder::SameParameter(E, B) — the sameparameter flag.
fn builder_same_parameter(the_brep: &mut BRep, the_e: &Shape, the_b: bool) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.same_parameter = the_b;
}

/// OCCT BRep_Builder::SameRange(E, B) — the samerange flag.
fn builder_same_range(the_brep: &mut BRep, the_e: &Shape, the_b: bool) {
    let a_ed = the_brep.edge_mut_inplace(the_e.clone());
    a_ed.same_range = the_b;
}

// ---------------------------------------------------------------------------
// GeomLib GAP carriers (TKGeomBase/GeomLib) — the translated bodies do not
// exist yet; the OCCT failure paths are preserved (the null handles make
// BuildCurve3d return false at the OCCT null checks L363 and L436).
// ---------------------------------------------------------------------------

/// OCCT GeomLib::To3d(Axes, Ptr2d) (GeomLib.cxx To3d) — GAP carrier.
fn geom_lib_to_3d(_the_axes: &GpAx2, _the_ptr2d: &Curve2d) -> Option<Curve3> {
    panic!("GAP: GeomLib::To3d (TKGeomBase/GeomLib not translated)")
}

/// OCCT GeomLib::BuildCurve3d(Tolerance, CurveOnSurface, First, Last,
/// Curve3d, Maxdev, AvDev, Continuity, MaxDegree, MaxSegment)
/// (GeomLib.cxx BuildCurve3d) — GAP carrier.
#[allow(clippy::too_many_arguments)]
fn geom_lib_build_curve3d(
    _the_tolerance: f64,
    _the_curve_on_surface: &Adaptor3dCurveOnSurface,
    _the_first: f64,
    _the_last: f64,
    _the_curve3d: &mut Option<Curve3>,
    _the_maxdev: &mut f64,
    _the_avdev: &mut f64,
    _the_continuity: GeomAbsShape,
    _the_max_degree: i32,
    _the_max_segment: i32,
) {
    panic!("GAP: GeomLib::BuildCurve3d (TKGeomBase/GeomLib not translated)")
}
