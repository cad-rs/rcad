//! OCCT FilletSurf package (TKFillet) — 1:1 translation.
//!
//! Sources:
//!   - FilletSurf_StatusDone.hxx (enum, 3 values)
//!   - FilletSurf_StatusType.hxx (enum, 3 values)
//!   - FilletSurf_ErrorTypeStatus.hxx (enum, 6 values)
//!   - FilletSurf_InternalBuilder.hxx L42-238 / .cxx L59-879
//!   - FilletSurf_Builder.hxx L47-153 / .cxx L32-388
//!
//! OCCT C++ inheritance (FilletSurf_InternalBuilder derives from
//! ChFi3d_FilBuilder) is modeled by composition: the derived builder
//! embeds `ChFi3dFilBuilder` as `base` (the established rcad convention,
//! chfi3d.rs).

use rcad_kernel::geom::{Curve2d, Curve3, CurveEval as _, Plane, Surface3, TrimmedCurve3};
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::{
    face_surface_value, surface_same, CurveRepresentation, Orientation, Shape, TShape,
};
use rcad_kernel::topods;

use super::chfi3d::ChFi3dFilBuilder;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_2::BRepTopAdaptorTopolTool;
use super::brep_blend_line::BRepBlendLine;
use super::chfi3d_builder_6::chfi3d_fil_common_point;
use super::chfi3d_builder_6b::ChFiDSElSpineHandle;
use super::chfi_ds::{
    ChFi3dFilletShape, ChFiDSElSpine, ChFiDSFilSpine, ChFiDSSpineHandle, SharedSurfData,
    SharedStripe,
};
use crate::geomalgo::gtests_stubs::GeomAbsShape;

// =========================================================================
// OCCT FilletSurf_StatusDone.hxx.
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilletSurfStatusDone {
    /// OCCT FilletSurf_IsOk.
    IsOk,
    /// OCCT FilletSurf_IsNotOk.
    IsNotOk,
    /// OCCT FilletSurf_IsPartial.
    IsPartial,
}

// =========================================================================
// OCCT FilletSurf_StatusType.hxx.
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilletSurfStatusType {
    /// OCCT FilletSurf_TwoExtremityOnEdge.
    TwoExtremityOnEdge,
    /// OCCT FilletSurf_OneExtremityOnEdge.
    OneExtremityOnEdge,
    /// OCCT FilletSurf_NoExtremityOnEdge.
    NoExtremityOnEdge,
}

// =========================================================================
// OCCT FilletSurf_ErrorTypeStatus.hxx.
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilletSurfErrorTypeStatus {
    /// OCCT FilletSurf_EmptyList.
    EmptyList,
    /// OCCT FilletSurf_EdgeNotG1.
    EdgeNotG1,
    /// OCCT FilletSurf_FacesNotG1.
    FacesNotG1,
    /// OCCT FilletSurf_EdgeNotOnShape.
    EdgeNotOnShape,
    /// OCCT FilletSurf_NotSharpEdge.
    NotSharpEdge,
    /// OCCT FilletSurf_PbFilletCompute.
    PbFilletCompute,
}

// =========================================================================
// OCCT FilletSurf_InternalBuilder.cxx — file-static helpers.
// =========================================================================

/// OCCT FilletSurf_InternalBuilder.cxx L59-70 — isinlist (file static).
fn isinlist(e: &Shape, l: &[Shape]) -> bool {
    for it in l {
        if e.is_same(it) {
            return true;
        }
    }
    false
}

/// OCCT FilletSurf_InternalBuilder.cxx L72-138 — IntPlanEdge (file static):
/// the intersection of the curve of <ed> with the plane <p>; <w> takes the
/// parameter of the intersection point closest to the plane origin (the
/// edge endpoints are checked through the plane projection, the ElSLib
/// plane evaluation of OCCT L107-136).
///
/// GAP-carrier note: the OCCT IntCurveSurface_HInter / BRepAdaptor_Curve /
/// GeomAdaptor_Surface triple is carried by the rcad IntCurveSurface
/// (geomalgo::int_patch::int_cs) over the edge curve value and the plane
/// Surface3; the endpoint projection is the pure-math plane stand-in.
#[allow(dead_code)]
#[allow(unused_assignments)] // OCCT L135 keeps the last `dist = d` write.
fn int_plan_edge(ed: &Shape, p: &Plane, w: &mut f64, tol3d: f64) -> bool {
    let mut done = false;
    // OCCT L78-79: f = Ed->FirstParameter(); l = Ed->LastParameter();
    let edata = ed.as_edge().expect("not an edge");
    let f = edata.range[0];
    let l = edata.range[1];
    let curve = edata.curve.clone().expect("edge curve");
    // OCCT L80: gp_Pnt Or = P.Location();
    let or = p.origin;
    // OCCT L81-82: Geom_Plane / GeomAdaptor_Surface over P.
    // OCCT L84-87: IntCurveSurface_HInter Intersection; dist = RealLast().
    let mut intersection = crate::geomalgo::int_patch::int_cs::IntCurveSurface::new();
    let mut dist = f64::INFINITY;
    // OCCT L89: Intersection.Perform(Ed, Plan);
    intersection.perform(&curve, &Surface3::Plane(*p), [f, l]);

    // OCCT L91-106: the nearest intersection point wins.
    if intersection.is_done() {
        let nbp = intersection.nb_points();
        for iip in 1..=nbp {
            let ip = intersection.point(iip);
            let pint = ip.p;
            let d = pint.distance(or);
            if d < dist {
                done = true;
                *w = ip.w();
                dist = d;
            }
        }
    }
    // OCCT L107-136: check if the extremities are not solution —
    // ElSLib::Parameters + ElSLib::Value on a plane is the orthogonal
    // projection; dproj = |(Pnt - Or) . Normal|.
    let pdeb = curve.point_at(f);
    let dprojdeb = (pdeb - or).dot(p.normal).abs();
    if dprojdeb < tol3d {
        let d = pdeb.distance(or);
        if d < dist {
            done = true;
            *w = f;
            dist = d;
        }
    }
    let pfin = curve.point_at(l);
    let dprojfin = (pfin - or).dot(p.normal).abs();
    if dprojfin < tol3d {
        let d = pfin.distance(or);
        if d < dist {
            done = true;
            *w = l;
            dist = d;
        }
    }
    done
}

/// OCCT FilletSurf_InternalBuilder.cxx L140-153 — ComputeEdgeParameter
/// (file static): the edge parameter <ped> matching the ElSpine parameter
/// <pelsp> through the normal plane at the guide evaluation point.
///
/// GAP: ChFiDS_ElSpine::D1 — the composite approximation curve on the
/// elementary spine is pending (the ChFi3d_PerformElSpine carrier,
/// chfi3d.rs L3335 note), so the guide evaluation point/tangent P/V and
/// the gp_Pln they build are unavailable; the OCCT not-done path
/// (ComputeEdgeParameter returning false) is kept.
#[allow(dead_code)]
fn compute_edge_parameter(
    spine: &ChFiDSSpineHandle,
    ind: usize,
    pelsp: f64,
    ped: &mut f64,
    tol3d: f64,
) -> bool {
    // OCCT L146: occ::handle<ChFiDS_ElSpine> Guide = Spine->ElSpine(ind);
    let Some(guide) = spine.base().el_spine_of_index(ind) else {
        return false;
    };
    // OCCT L148-152: Guide->D1(pelsp, P, V); gp_Pln pln(P, V);
    // ed = BRepAdaptor_Curve(Spine->CurrentElementarySpine(ind));
    // return IntPlanEdge(ed, pln, ped, tol3d);
    let _ = (guide, pelsp, ped, tol3d);
    false
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) stand-in: the regularity record
/// (BRep_CurveOn2Surfaces) on the edge matching the two face surfaces;
/// GeomAbs_C0 when no record matches.
///
/// GAP-carrier note: rcad-kernel exposes only the BRep_Builder::Continuity
/// setter; this read-back over the CurveOn2Surfaces representations is the
/// BRep_Tool static, kept local to the consumer pending the kernel move.
/// The OCCT location composition (L1/L2 vs the record locations) reduces
/// to the surface value match (rcad surfaces travel as values).
fn brep_tool_continuity(brep: &topods::BRep, e: &Shape, f1: &Shape, f2: &Shape) -> GeomAbsShape {
    let (Some(s1), Some(s2)) = (
        face_surface_value(brep, f1),
        face_surface_value(brep, f2),
    ) else {
        return GeomAbsShape::C0;
    };
    let Some(ts) = brep.tshapes.get(e.index) else {
        return GeomAbsShape::C0;
    };
    let TShape::Edge(ed) = ts.as_ref() else {
        return GeomAbsShape::C0;
    };
    for cr in &ed.representations {
        if let CurveRepresentation::CurveOn2Surfaces {
            surface1,
            surface2,
            continuity,
            ..
        } = cr
        {
            if (surface_same(surface1, &s1) && surface_same(surface2, &s2))
                || (surface_same(surface1, &s2) && surface_same(surface2, &s1))
            {
                // The kernel stores the topods GeomAbsShape; the ChFi3d
                // pipeline carries the geomalgo one (same OCCT enum, same
                // order) — the value-preserving conversion.
                return match continuity {
                    topods::GeomAbsShape::C0 => GeomAbsShape::C0,
                    topods::GeomAbsShape::G1 => GeomAbsShape::G1,
                    topods::GeomAbsShape::C1 => GeomAbsShape::C1,
                    topods::GeomAbsShape::G2 => GeomAbsShape::G2,
                    topods::GeomAbsShape::C2 => GeomAbsShape::C2,
                    topods::GeomAbsShape::C3 => GeomAbsShape::C3,
                    topods::GeomAbsShape::CN => GeomAbsShape::CN,
                };
            }
        }
    }
    GeomAbsShape::C0
}

/// rcad architecture splice: OCCT `st->ChangeSpine()` is a live reference
/// into the stripe; the taken handle is stored back where the OCCT
/// reference lifetime ends.
fn restore_spine(stripe: &SharedStripe, sp: ChFiDSSpineHandle) {
    stripe.write().expect("stripe lock").my_spine = Some(sp);
}

// =========================================================================
// OCCT FilletSurf_InternalBuilder (InternalBuilder.hxx L42-238).
// =========================================================================

/// OCCT FilletSurf_InternalBuilder — derives from ChFi3d_FilBuilder
/// (composition: `base`).
#[derive(Debug, Clone)]
pub struct FilletSurfInternalBuilder {
    pub base: ChFi3dFilBuilder,
}

impl FilletSurfInternalBuilder {
    /// OCCT InternalBuilder.cxx L157-166 — ChFi3d_FilBuilder(S, FShape, Ta),
    /// SetParams(Ta, Tapp3d, Tapp2d, Tapp3d, Tapp2d, 1.e-3),
    /// SetContinuity(GeomAbs_C2, Ta).
    pub fn new(
        brep: &topods::BRep,
        s: &Shape,
        fshape: ChFi3dFilletShape,
        ta: f64,
        tapp3d: f64,
        tapp2d: f64,
    ) -> Self {
        let mut base = ChFi3dFilBuilder::new(brep, s.clone(), fshape, ta);
        base.base.set_params(ta, tapp3d, tapp2d, tapp3d, tapp2d, 1.0e-3);
        base.base.set_continuity(GeomAbsShape::C2, ta);
        FilletSurfInternalBuilder { base }
    }

    /// OCCT: occ::handle<ChFiDS_Stripe> Stripe = myListStripe.First();
    fn first_stripe(&self) -> SharedStripe {
        self.base
            .base
            .my_list_stripe
            .first()
            .cloned()
            .expect("stripe")
    }

    /// OCCT: myListStripe.First()->SetOfSurfData()->Value(Index) (1-based).
    fn surf_data(&self, index: usize) -> SharedSurfData {
        self.first_stripe()
            .read()
            .expect("stripe lock")
            .set_of_surf_data()[index - 1]
            .clone()
    }

    /// OCCT InternalBuilder.cxx L180-317 — Add: creation of spine on a set
    /// of edges.  Return codes: 0 no problem, 1 empty list, 2 non G1 edges,
    /// 3 non G1 adjacent faces, 4 edge is not on the shape, 5 edge is not
    /// alive (not sharp).
    pub fn add(&mut self, e: &[Shape], r: f64) -> i32 {
        if e.is_empty() {
            return 1;
        }
        for cured in e {
            // OCCT L189-193: TopoDS::Edge(It.Value()) null check, then
            // myEFMap.Contains(cured).
            if cured.is_null() {
                return 4;
            }
            if !self.base.base.my_ef_map.contains(cured) {
                return 4;
            }
            // OCCT L198-214: check if the edge is a fracture edge — the two
            // first distinct faces of myEFMap(cured).
            let mut ff1 = Shape::null();
            let mut ff2 = Shape::null();
            for f in self.base.base.my_ef_map.find(cured) {
                if ff1.is_null() {
                    ff1 = f.clone();
                } else {
                    ff2 = f.clone();
                    if !ff2.is_same(&ff1) {
                        break;
                    }
                }
            }
            if ff1.is_null() || ff2.is_null() {
                return 5;
            }
            if ff1.is_same(&ff2) {
                return 5;
            }
            // OCCT L223-226: BRep_Tool::Continuity(cured, ff1, ff2) must be
            // GeomAbs_C0 (the edge is sharp).
            if brep_tool_continuity(&self.base.base.my_brep, cured, &ff1, &ff2) != GeomAbsShape::C0
            {
                return 5;
            }
        }
        // OCCT L228-230: ed = TopoDS::Edge(E.First()); ed.Orientation(FORWARD);
        // ChFi3d_FilBuilder::Add(R, ed);
        let mut ed = e[0].clone();
        ed.orientation = Orientation::Forward;
        self.base.add_radius(r, &ed);
        // OCCT L231-233: st = myListStripe.First(); sp = st->ChangeSpine();
        // periodic = sp->IsPeriodic();
        let stripe = self.first_stripe();
        let mut sp = stripe
            .write()
            .expect("stripe lock")
            .my_spine
            .take()
            .expect("null spine");
        let periodic = sp.base().is_periodic();

        // OCCT L239-246: it is checked if edges of list E are in the contour
        // (the edges that are not in the list stop with the non-G1 code).
        for cured in e {
            if self.base.base.contains(cured) == 0 {
                restore_spine(&stripe, sp);
                return 2;
            }
        }

        // OCCT L248-269: a new FilSpine is filled with the edges of E in the
        // contour order; a gap (a contour edge that is not in E after the
        // list started matching) stops with the non-G1 code.
        let mut newsp = ChFiDSSpineHandle::Fil(ChFiDSFilSpine::new());
        let mut debut = false;
        let mut premierquinyestpas: usize = 0;
        let mut yatrou: i32 = 0;
        let nb_edges = sp.base().nb_edges();
        for i in 1..=nb_edges {
            let cured = sp.base().edges(i).clone();
            if isinlist(&cured, e) {
                debut = true;
                if premierquinyestpas != 0 {
                    yatrou = 1;
                    break;
                }
                newsp.base_mut().set_edges(cured);
            } else if debut && premierquinyestpas == 0 {
                premierquinyestpas = i;
            }
        }
        if !periodic && yatrou != 0 {
            restore_spine(&stripe, sp);
            return 2;
        }
        if periodic && yatrou != 0 {
            // OCCT L274-294: the periodic tail lying before the gap is moved
            // in front of the new spine.
            let mut vraitrou = false;
            let mut a_local_debut = false;
            let mut i = sp.base().nb_edges() as i32;
            while i > yatrou {
                let cured = sp.base().edges(i as usize).clone();
                if isinlist(&cured, e) {
                    if vraitrou {
                        restore_spine(&stripe, sp);
                        return 2;
                    }
                    newsp.base_mut().put_in_first(cured);
                } else if a_local_debut {
                    vraitrou = true;
                }
                a_local_debut = true;
                i -= 1;
            }
        }

        // OCCT L296-301: if the new spine differs from the contour spine it
        // is loaded, the radius is set, and it becomes the stripe spine.
        if newsp.base().nb_edges() != sp.base().nb_edges() {
            newsp.base_mut().load();
            if let Some(fs) = newsp.down_cast_fil_mut() {
                fs.set_radius(r);
            }
            sp = newsp;
        }

        // OCCT L303-315: the ElSpine is immediately constructed.
        let mut hels = ChFiDSElSpine::new();
        let firstparam = sp.base().first_parameter();
        let lastparam = sp.base().last_parameter();
        let (pfirst, tfirst) = sp.base_mut().d1(firstparam);
        let (plast, tlast) = sp.base_mut().d1(lastparam);
        hels.firstparam = firstparam; // L309: hels->FirstParameter(...)
        hels.set_first_point_and_tgt(pfirst, tfirst); // L310
        hels.lastparam = lastparam; // L311
        hels.set_last_point_and_tgt(plast, tlast); // L312
        // OCCT L313: ChFi3d_PerformElSpine(hels, sp, myConti, tolesp);
        // GAP: the composite ElSpine approximation (ChFi3d_Builder_0.cxx
        // ChFi3d_PerformElSpine) is pending in rcad (chfi3d.rs L3335 note);
        // the handle keeps the constructed endpoints/tangents above.
        // OCCT L314: sp->AppendElSpine(hels) — the ChFiDS_FilSpine override
        // (ChFiDS_FilSpine.cxx L369-373) also appends the ComputeLaw
        // composite; the base append is the non-fillet fallback.
        match sp.down_cast_fil_mut() {
            Some(fsp) => fsp.append_el_spine_fil(&hels),
            None => sp.base_mut().append_el_spine(hels),
        } // L314
        sp.base_mut().set_split_done(true); // L315
        restore_spine(&stripe, sp);
        0
    }

    /// OCCT InternalBuilder.cxx L321-337 — Perform.
    pub fn perform(&mut self) {
        // PerformSetOfSurfOnElSpine is enough.
        let stripe = self.first_stripe();
        {
            // OCCT L326-328: HData = Stripe->ChangeSetOfSurfData();
            // HData = new NCollection_HSequence<...>();
            stripe.write().expect("stripe lock").my_hdata = Vec::new();
        }
        // OCCT L329: Spine = Stripe->ChangeSpine();
        let sp = stripe
            .write()
            .expect("stripe lock")
            .my_spine
            .take()
            .expect("null spine");
        // OCCT L330-332: StripeOrientations(Spine, RefOr1, RefOr2, RefChoix);
        let (ref_or1, ref_or2, ref_choix) = self.base.base.stripe_orientations(&sp);
        {
            let mut st = stripe.write().expect("stripe lock");
            st.set_orientation_on_face1(ref_or1); // L333
            st.set_orientation_on_face2(ref_or2); // L334
            st.set_choix(ref_choix); // L335
        }
        restore_spine(&stripe, sp);
        // OCCT L336: PerformSetOfKGen(Stripe, false);
        self.base.base.perform_set_of_k_gen(&stripe, false);
    }

    /// OCCT InternalBuilder.cxx L341-504 — PerformSurf (the fillet/fillet
    /// overload; S1 and S2 both bounded by a TopolTool).  Over the rcad
    /// composition model the OCCT virtual dispatch of PerformSetOfKGen to
    /// this override is not live; the method keeps the OCCT body form.
    #[allow(clippy::too_many_arguments, unreachable_code, dead_code)]
    pub fn perform_surf(
        &mut self,
        seq_data: &mut Vec<SharedSurfData>,
        guide: &ChFiDSElSpineHandle,
        spine: &mut ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &Vector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        // OCCT L364: occ::handle<ChFiDS_SurfData> Data = SeqData(1);
        let data = seq_data.first().cloned().expect("surfdata");
        // OCCT L365-369: fsp = down_cast<ChFiDS_FilSpine>(Spine); null ->
        // Standard_ConstructionError.
        let Some(fsp) = spine.down_cast_fil() else {
            panic!("PerformSurf : this is not the spine of a fillet");
        };
        // OCCT L370: occ::handle<BRepBlend_Line> lin;
        #[allow(unused_mut)] // OCCT SimulData fills the handle (GAP above).
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT L371: TopAbs_Orientation Or = S1->Face().Orientation();
        let _or = s1.face.orientation;
        // OCCT L372-375: !fsp->IsConstant() -> Standard_ConstructionError.
        if !fsp.is_constant() {
            panic!("PerformSurf : no variable radiuses");
        }
        // OCCT L376: (comment) bool maybesingular;
        // OCCT L378: occ::handle<ChFiDS_ElSpine> EmptyGuide;
        let _empty_guide: Option<ChFiDSElSpineHandle> = None;
        // OCCT L395: double PFirst = First;
        let _p_first = *first;

        // OCCT L380-394: BRepBlend_ConstRad Func(S1, S2, Guide);
        // BRepBlend_ConstRadInv FInv(S1, S2, Guide);
        // Func.Set(fsp->Radius(), Choix); FInv.Set(fsp->Radius(), Choix);
        // switch (GetFilletShape()) { Rational -> Func.Set(BlendFunc_Rational);
        // QuasiAngular -> BlendFunc_QuasiAngular; Polynomial ->
        // BlendFunc_Polynomial }.
        // GAP: the rcad BlendFuncConstRad (brep_blend_func_consrad.rs)
        // carries neither the BlendFunction trait wiring nor the fillet
        // shape setter, and the ElSpine guide curve is pending (see the
        // Add GAP); the construction stays a structural pending boundary.
        let _ = (
            guide, choix, i1, i2, max_step, fleche, tol_guide, inside, appro, rec_on_s1,
            rec_on_s2, soldep,
        );

        // OCCT L396-418: done = SimulData(Data, Guide, EmptyGuide, lin, S1,
        // I1, S2, I2, Func, FInv, PFirst, MaxStep, Fleche, TolGuide, First,
        // Last, Inside, Appro, Forward, Soldep, 20, RecOnS1, RecOnS2);
        // GAP: the SimulData overload over Blend_Function
        // (ChFi3d_Builder_6.cxx L4390) has no rcad carrier yet; the OCCT
        // not-done path (L419-422) is kept.
        let done = false;
        if !done {
            return false;
        }
        let lin_ref = lin.as_ref().expect("BRepBlend_Line");
        // OCCT L423-430.
        if lin_ref.start_point_on_first().nb_point_on_rst() != 0 {
            let mut st = data.write().expect("surfdata lock");
            chfi3d_fil_common_point(
                &self.base.base.my_brep,
                lin_ref.start_point_on_first(),
                lin_ref.transition_on_s1(),
                true,
                st.change_vertex(true, 1),
                self.base.base.tolapp3d,
            );
        }
        // OCCT L431-438.
        if lin_ref.end_point_on_first().nb_point_on_rst() != 0 {
            let mut st = data.write().expect("surfdata lock");
            chfi3d_fil_common_point(
                &self.base.base.my_brep,
                lin_ref.end_point_on_first(),
                lin_ref.transition_on_s1(),
                false,
                st.change_vertex(false, 1),
                self.base.base.tolapp3d,
            );
        }
        // OCCT L439-446.
        if lin_ref.start_point_on_second().nb_point_on_rst() != 0 {
            let mut st = data.write().expect("surfdata lock");
            chfi3d_fil_common_point(
                &self.base.base.my_brep,
                lin_ref.start_point_on_second(),
                lin_ref.transition_on_s2(),
                true,
                st.change_vertex(true, 2),
                self.base.base.tolapp3d,
            );
        }
        // OCCT L447-454.
        if lin_ref.end_point_on_second().nb_point_on_rst() != 0 {
            let mut st = data.write().expect("surfdata lock");
            chfi3d_fil_common_point(
                &self.base.base.my_brep,
                lin_ref.end_point_on_second(),
                lin_ref.transition_on_s2(),
                false,
                st.change_vertex(false, 2),
                self.base.base.tolapp3d,
            );
        }
        // OCCT L455: done = CompleteData(Data, Func, lin, S1, S2, Or,
        // false, false, false, false);
        // GAP: the CompleteData overload over Blend_Function
        // (ChFi3d_Builder_6.cxx, rcad complete_data_function) has no
        // BlendFunction carrier for BlendFuncConstRad; the Standard_Failure
        // throw (L458-459) stays as the structural pending boundary.
        let done = true;
        if !done {
            panic!("PerformSurf : Failed approximation!");
        }
        // OCCT L461-480: the Intf walk over the first vertices.
        let mut ok = false;
        if !forward {
            *intf = 0;
            let cpf1 = data
                .read()
                .expect("surfdata lock")
                .vertex_first_on_s1()
                .clone();
            if cpf1.is_on_arc() {
                let f1 = s1.face.clone();
                let mut bid = Shape::null();
                *intf = i32::from(!self.base.base.search_face(spine, &cpf1, &f1, &mut bid));
                ok = *intf != 0;
            }
            let cpf2 = data
                .read()
                .expect("surfdata lock")
                .vertex_first_on_s2()
                .clone();
            if cpf2.is_on_arc() && !ok {
                let f2 = s2.face.clone();
                let mut bid = Shape::null();
                *intf = i32::from(!self.base.base.search_face(spine, &cpf2, &f2, &mut bid));
            }
        }
        // OCCT L481-497: the Intl walk over the last vertices.
        *intl = 0;
        ok = false;
        let cpl1 = data
            .read()
            .expect("surfdata lock")
            .vertex_last_on_s1()
            .clone();
        if cpl1.is_on_arc() {
            let f1 = s1.face.clone();
            let mut bid = Shape::null();
            *intl = i32::from(!self.base.base.search_face(spine, &cpl1, &f1, &mut bid));
            ok = *intl != 0;
        }
        let cpl2 = data
            .read()
            .expect("surfdata lock")
            .vertex_last_on_s2()
            .clone();
        if cpl2.is_on_arc() && !ok {
            let f2 = s2.face.clone();
            let mut bid = Shape::null();
            *intl = i32::from(!self.base.base.search_face(spine, &cpl2, &f2, &mut bid));
        }
        // OCCT L498-499: Data->FirstSpineParam(First);
        // Data->LastSpineParam(Last);
        {
            let mut st = data.write().expect("surfdata lock");
            st.set_first_spine_param(*first);
            st.set_last_spine_param(*last);
        }

        // OCCT L500-503: (comments) maybesingular / SplitSurf / trim.
        true
    }

    /// OCCT InternalBuilder.cxx L506-533 — PerformSurf (the Surf/Rst
    /// overload): throws Standard_DomainError.
    #[allow(clippy::too_many_arguments, dead_code, unused_variables)]
    pub fn perform_surf_surf_rst(
        &mut self,
        seq_data: &mut Vec<SharedSurfData>,
        guide: &ChFiDSElSpineHandle,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        pc1: &Curve2d,
        sref1: &BRepAdaptorSurface,
        pcref1: &Curve2d,
        decroch1: &mut bool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        or2: Orientation,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &Vector,
    ) {
        panic!("BlendFunc_CSConstRad::Section : Not implemented");
    }

    /// OCCT InternalBuilder.cxx L535-562 — PerformSurf (the Rst/Surf
    /// overload): throws Standard_DomainError.
    #[allow(clippy::too_many_arguments, dead_code, unused_variables)]
    pub fn perform_surf_rst_surf(
        &mut self,
        seq_data: &mut Vec<SharedSurfData>,
        guide: &ChFiDSElSpineHandle,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        or1: Orientation,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        pc2: &Curve2d,
        sref2: &BRepAdaptorSurface,
        pcref2: &Curve2d,
        decroch2: &mut bool,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &Vector,
    ) {
        panic!("BlendFunc_CSConstRad::Section : Not implemented");
    }

    /// OCCT InternalBuilder.cxx L564-597 — PerformSurf (the Rst/Rst
    /// overload): throws Standard_DomainError.
    #[allow(clippy::too_many_arguments, dead_code, unused_variables)]
    pub fn perform_surf_rst_rst(
        &mut self,
        seq_data: &mut Vec<SharedSurfData>,
        guide: &ChFiDSElSpineHandle,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        pc1: &Curve2d,
        sref1: &BRepAdaptorSurface,
        pcref1: &Curve2d,
        decroch1: &mut bool,
        or1: Orientation,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        pc2: &Curve2d,
        sref2: &BRepAdaptorSurface,
        pcref2: &Curve2d,
        decroch2: &mut bool,
        or2: Orientation,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p1: bool,
        rec_rst1: bool,
        rec_p2: bool,
        rec_rst2: bool,
        soldep: &Vector,
    ) {
        panic!("BlendFunc_CSConstRad::Section : Not implemented");
    }

    /// OCCT InternalBuilder.cxx L599-602 — Done: the ChFi3d_Builder
    /// protected `done` field.
    pub fn done(&self) -> bool {
        self.base.base.done
    }

    /// OCCT InternalBuilder.cxx L609-612 — NbSurface.
    pub fn nb_surface(&self) -> usize {
        // myListStripe.First()->SetOfSurfData()->Length();
        self.first_stripe()
            .read()
            .expect("stripe lock")
            .set_of_surf_data()
            .len()
    }

    /// OCCT InternalBuilder.cxx L619-624 — SurfaceFillet.
    pub fn surface_fillet(&self, index: usize) -> Surface3 {
        // OCCT L621: isurf = ...->Surf();
        let isurf = self
            .surf_data(index)
            .read()
            .expect("surfdata lock")
            .surf();
        // OCCT L623: myDS->Surface(isurf).Surface();
        let dstr = self.base.base.my_ds.as_ref().expect("DS");
        dstr.surface(isurf).surface().clone()
    }

    /// OCCT InternalBuilder.cxx L631-636 — TolApp3d.
    pub fn tol_app3d(&self, index: usize) -> f64 {
        let isurf = self
            .surf_data(index)
            .read()
            .expect("surfdata lock")
            .surf();
        // OCCT L635: myDS->Surface(isurf).Tolerance();
        let dstr = self.base.base.my_ds.as_ref().expect("DS");
        dstr.surface(isurf).tolerance()
    }

    /// OCCT InternalBuilder.cxx L642-647 — SupportFace1.
    pub fn support_face1(&self, index: usize) -> Shape {
        // OCCT L644: isurf = ...->IndexOfS1();
        let isurf = self
            .surf_data(index)
            .read()
            .expect("surfdata lock")
            .index_of(1);
        // OCCT L646: TopoDS::Face(myDS->Shape(isurf)); the rcad Shape is the
        // untyped handle, the cast is the value copy.
        let dstr = self.base.base.my_ds.as_ref().expect("DS");
        dstr.shape(isurf).clone()
    }

    /// OCCT InternalBuilder.cxx L653-658 — SupportFace2.
    pub fn support_face2(&self, index: usize) -> Shape {
        let isurf = self
            .surf_data(index)
            .read()
            .expect("surfdata lock")
            .index_of(2);
        let dstr = self.base.base.my_ds.as_ref().expect("DS");
        dstr.shape(isurf).clone()
    }

    /// OCCT InternalBuilder.cxx L665-669 — CurveOnFace1.
    pub fn curve_on_face1(&self, index: usize) -> Option<Curve3> {
        // OCCT L667: icurv = ...->InterferenceOnS1().LineIndex();
        let icurv = self
            .surf_data(index)
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .line_index();
        // OCCT L668: myDS->Curve(icurv).Curve();
        let dstr = self.base.base.my_ds.as_ref().expect("DS");
        dstr.curve(icurv).curve().cloned()
    }

    /// OCCT InternalBuilder.cxx L676-680 — CurveOnFace2.
    pub fn curve_on_face2(&self, index: usize) -> Option<Curve3> {
        let icurv = self
            .surf_data(index)
            .read()
            .expect("surfdata lock")
            .interference_on_s2()
            .line_index();
        let dstr = self.base.base.my_ds.as_ref().expect("DS");
        dstr.curve(icurv).curve().cloned()
    }

    /// OCCT InternalBuilder.cxx L686-689 — PCurveOnFace1.
    pub fn pcurve_on_face1(&self, index: usize) -> Option<Curve2d> {
        self.surf_data(index)
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .pcurve_on_face()
            .cloned()
    }

    /// OCCT InternalBuilder.cxx L696-699 — PCurve1OnFillet.
    pub fn pcurve1_on_fillet(&self, index: usize) -> Option<Curve2d> {
        self.surf_data(index)
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .pcurve_on_surf()
            .cloned()
    }

    /// OCCT InternalBuilder.cxx L705-708 — PCurveOnFace2.
    pub fn pcurve_on_face2(&self, index: usize) -> Option<Curve2d> {
        self.surf_data(index)
            .read()
            .expect("surfdata lock")
            .interference_on_s2()
            .pcurve_on_face()
            .cloned()
    }

    /// OCCT InternalBuilder.cxx L714-717 — PCurve2OnFillet.
    pub fn pcurve2_on_fillet(&self, index: usize) -> Option<Curve2d> {
        self.surf_data(index)
            .read()
            .expect("surfdata lock")
            .interference_on_s2()
            .pcurve_on_surf()
            .cloned()
    }

    /// OCCT InternalBuilder.cxx L724-741 — FirstParameter.
    pub fn first_parameter(&self) -> f64 {
        let st = self.first_stripe();
        // OCCT L726-729: sp = st->Spine(); sd = st->SetOfSurfData()->Value(1);
        // p = sd->FirstSpineParam();
        let (mut sp, p) = {
            let guard = st.read().expect("stripe lock");
            let sp = guard.my_spine.clone().expect("null spine");
            let sd = guard.set_of_surf_data()[0].clone();
            let p = sd.read().expect("surfdata lock").first_spine_param();
            (sp, p)
        };
        // OCCT L730-734: ind = 1; if (sp->IsPeriodic()) ind = sp->Index(p);
        let mut ind = 1usize;
        if sp.base().is_periodic() {
            ind = sp.base_mut().index_of_param(p, true);
        }
        // OCCT L735-741: if (ComputeEdgeParameter(...)) return ep;
        let mut ep = 0.0;
        if compute_edge_parameter(&sp, ind, p, &mut ep, self.base.base.tolapp3d) {
            return ep;
        }
        0.0
    }

    /// OCCT InternalBuilder.cxx L747-764 — LastParameter.
    pub fn last_parameter(&self) -> f64 {
        let st = self.first_stripe();
        // OCCT L749-752: sp = st->Spine(); sd = Value(NbSurface());
        // p = sd->LastSpineParam();
        let (mut sp, p) = {
            let guard = st.read().expect("stripe lock");
            let sp = guard.my_spine.clone().expect("null spine");
            let sd = guard.set_of_surf_data()[self.nb_surface() - 1].clone();
            let p = sd.read().expect("surfdata lock").last_spine_param();
            (sp, p)
        };
        // OCCT L753-757: ind = sp->NbEdges(); if periodic ind = sp->Index(p);
        let mut ind = sp.base().nb_edges();
        if sp.base().is_periodic() {
            ind = sp.base_mut().index_of_param(p, true);
        }
        // OCCT L758-763.
        let mut ep = 0.0;
        if compute_edge_parameter(&sp, ind, p, &mut ep, self.base.base.tolapp3d) {
            return ep;
        }
        0.0
    }

    /// OCCT InternalBuilder.cxx L777-795 — StartSectionStatus.
    pub fn start_section_status(&self) -> FilletSurfStatusType {
        let sd = self.surf_data(1);
        let guard = sd.read().expect("surfdata lock");
        let isonedge1 = guard.vertex_first_on_s1().is_on_arc();
        let isonedge2 = guard.vertex_first_on_s2().is_on_arc();

        if isonedge1 && isonedge2 {
            FilletSurfStatusType::TwoExtremityOnEdge
        } else if !isonedge1 && !isonedge2 {
            FilletSurfStatusType::NoExtremityOnEdge
        } else {
            FilletSurfStatusType::OneExtremityOnEdge
        }
    }

    /// OCCT InternalBuilder.cxx L807-826 — EndSectionStatus.
    pub fn end_section_status(&self) -> FilletSurfStatusType {
        let sd = self.surf_data(self.nb_surface());
        let guard = sd.read().expect("surfdata lock");
        let isonedge1 = guard.vertex_last_on_s1().is_on_arc();
        let isonedge2 = guard.vertex_last_on_s2().is_on_arc();

        if isonedge1 && isonedge2 {
            FilletSurfStatusType::TwoExtremityOnEdge
        } else if !isonedge1 && !isonedge2 {
            FilletSurfStatusType::NoExtremityOnEdge
        } else {
            FilletSurfStatusType::OneExtremityOnEdge
        }
    }

    /// OCCT InternalBuilder.cxx L832-847 — Simulate.
    pub fn simulate(&mut self) {
        // ChFi3d_FilBuilder::Simulate(1);
        let stripe = self.first_stripe();
        {
            // OCCT L835-838: HData = Stripe->ChangeSetOfSurfData();
            // HData = new NCollection_HSequence<...>();
            stripe.write().expect("stripe lock").my_hdata = Vec::new();
        }
        // OCCT L839: Spine = Stripe->ChangeSpine();
        let sp = stripe
            .write()
            .expect("stripe lock")
            .my_spine
            .take()
            .expect("null spine");
        // OCCT L840-842: StripeOrientations(Spine, RefOr1, RefOr2, RefChoix);
        let (ref_or1, ref_or2, ref_choix) = self.base.base.stripe_orientations(&sp);
        {
            let mut st = stripe.write().expect("stripe lock");
            st.set_orientation_on_face1(ref_or1);
            st.set_orientation_on_face2(ref_or2);
            st.set_choix(ref_choix);
        }
        restore_spine(&stripe, sp);
        // OCCT L846: PerformSetOfKGen(Stripe, true);
        self.base.base.perform_set_of_k_gen(&stripe, true);
    }

    /// OCCT InternalBuilder.cxx L853-856 — NbSection.
    pub fn nb_section(&self, index_surf: usize) -> usize {
        // OCCT L855: Sect(1, IndexSurf)->Length();
        self.base.sect(1, index_surf).expect("Sect").len()
    }

    /// OCCT InternalBuilder.cxx L864-872 — Section: the trimmed circular
    /// arc of section <index_sec> of surface <index_surf>.
    pub fn section(&self, index_surf: usize, index_sec: usize) -> TrimmedCurve3 {
        // OCCT L869: Sect(1, IndexSurf)->Value(IndexSec).Get(c, deb, fin);
        let sections = self.base.sect(1, index_surf).expect("Sect");
        let (c, deb, fin) = sections[index_sec - 1].get_circ();
        // OCCT L870-871: Gc = new Geom_Circle(c);
        // return new Geom_TrimmedCurve(Gc, deb, fin);
        TrimmedCurve3::new(Curve3::Circle(c), deb, fin)
    }

    // OCCT InternalBuilder.cxx L874-879 — the deprecated Section overload
    // taking `handle<Geom_TrimmedCurve>& Circ` as an out-parameter collapses
    // into the by-value `section` above (no Rust out-parameters).
}

// =========================================================================
// OCCT FilletSurf_Builder (Builder.hxx L47-153 / .cxx L32-388).
// =========================================================================

/// OCCT FilletSurf_Builder — the API facade over the InternalBuilder.
#[derive(Debug, Clone)]
pub struct FilletSurfBuilder {
    /// OCCT: FilletSurf_InternalBuilder myIntBuild.
    pub my_int_build: FilletSurfInternalBuilder,
    /// OCCT: FilletSurf_StatusDone myisdone.
    pub myisdone: FilletSurfStatusDone,
    /// OCCT: FilletSurf_ErrorTypeStatus myerrorstatus.
    pub myerrorstatus: FilletSurfErrorTypeStatus,
}

impl FilletSurfBuilder {
    /// OCCT Builder.cxx L32-67 — the constructor with the OCCT default
    /// tolerances (Ta = 1.0e-2, Tapp3d = 1.0e-4, Tapp2d = 1.0e-5).
    pub fn new(brep: &topods::BRep, s: &Shape, e: &[Shape], r: f64) -> Self {
        FilletSurfBuilder::new_with_tolerances(brep, s, e, r, 1.0e-2, 1.0e-4, 1.0e-5)
    }

    /// OCCT Builder.cxx L32-67 — the constructor with explicit tolerances.
    pub fn new_with_tolerances(
        brep: &topods::BRep,
        s: &Shape,
        e: &[Shape],
        r: f64,
        ta: f64,
        tapp3d: f64,
        tapp2d: f64,
    ) -> Self {
        // OCCT L38: myIntBuild(S, ChFi3d_Polynomial, Ta, Tapp3d, Tapp2d);
        let my_int_build = FilletSurfInternalBuilder::new(
            brep,
            s,
            ChFi3dFilletShape::Polynomial,
            ta,
            tapp3d,
            tapp2d,
        );
        let mut res = FilletSurfBuilder {
            my_int_build,
            myisdone: FilletSurfStatusDone::IsOk, // L40
            myerrorstatus: FilletSurfErrorTypeStatus::EmptyList, // L41
        };
        // OCCT L42-66: add = myIntBuild.Add(E, R); the add code selects the
        // error status.
        let add = res.my_int_build.add(e, r);
        if add != 0 {
            res.myisdone = FilletSurfStatusDone::IsNotOk;
            if add == 1 {
                res.myerrorstatus = FilletSurfErrorTypeStatus::EmptyList;
            } else if add == 2 {
                res.myerrorstatus = FilletSurfErrorTypeStatus::EdgeNotG1;
            } else if add == 3 {
                res.myerrorstatus = FilletSurfErrorTypeStatus::FacesNotG1;
            } else if add == 4 {
                res.myerrorstatus = FilletSurfErrorTypeStatus::EdgeNotOnShape;
            } else if add == 5 {
                res.myerrorstatus = FilletSurfErrorTypeStatus::NotSharpEdge;
            }
        }
        res
    }

    /// OCCT Builder.cxx L73-93 — Perform: computation of the fillet.
    pub fn perform(&mut self) {
        if self.myisdone == FilletSurfStatusDone::IsOk {
            self.my_int_build.perform();
            if self.my_int_build.done() {
                self.myisdone = FilletSurfStatusDone::IsOk;
            } else if self.my_int_build.nb_surface() != 0 {
                self.myisdone = FilletSurfStatusDone::IsPartial;
                self.myerrorstatus = FilletSurfErrorTypeStatus::PbFilletCompute;
            } else {
                self.myisdone = FilletSurfStatusDone::IsNotOk;
                self.myerrorstatus = FilletSurfErrorTypeStatus::PbFilletCompute;
            }
        }
    }

    /// OCCT Builder.cxx L99-102 — IsDone.
    pub fn is_done(&self) -> FilletSurfStatusDone {
        self.myisdone
    }

    /// OCCT Builder.cxx L108-111 — StatusError.
    pub fn status_error(&self) -> FilletSurfErrorTypeStatus {
        self.myerrorstatus
    }

    /// OCCT Builder.cxx L118-125 — NbSurface (raises StdFail_NotDone).
    pub fn nb_surface(&self) -> usize {
        if self.is_done() != FilletSurfStatusDone::IsNotOk {
            return self.my_int_build.nb_surface();
        }
        panic!("StdFail_NotDone: FilletSurf_Builder::NbSurface");
    }

    /// OCCT Builder.cxx L132-139 — SurfaceFillet (raises
    /// Standard_OutOfRange).
    pub fn surface_fillet(&self, index: usize) -> Surface3 {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::SurfaceFillet");
        }
        self.my_int_build.surface_fillet(index)
    }

    /// OCCT Builder.cxx L145-152 — TolApp3d (raises Standard_OutOfRange).
    pub fn tol_app3d(&self, index: usize) -> f64 {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::TolApp3d");
        }
        self.my_int_build.tol_app3d(index)
    }

    /// OCCT Builder.cxx L158-165 — SupportFace1 (raises
    /// Standard_OutOfRange).
    pub fn support_face1(&self, index: usize) -> Shape {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::SupportFace1");
        }
        self.my_int_build.support_face1(index)
    }

    /// OCCT Builder.cxx L171-178 — SupportFace2 (raises
    /// Standard_OutOfRange).
    pub fn support_face2(&self, index: usize) -> Shape {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::SupportFace2");
        }
        self.my_int_build.support_face2(index)
    }

    /// OCCT Builder.cxx L184-191 — CurveOnFace1 (raises
    /// Standard_OutOfRange).
    pub fn curve_on_face1(&self, index: usize) -> Option<Curve3> {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::CurveOnFace1");
        }
        self.my_int_build.curve_on_face1(index)
    }

    /// OCCT Builder.cxx L197-204 — CurveOnFace2 (raises
    /// Standard_OutOfRange).
    pub fn curve_on_face2(&self, index: usize) -> Option<Curve3> {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::CurveOnFace2");
        }
        self.my_int_build.curve_on_face2(index)
    }

    /// OCCT Builder.cxx L210-217 — PCurveOnFace1 (raises
    /// Standard_OutOfRange).
    pub fn pcurve_on_face1(&self, index: usize) -> Option<Curve2d> {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::PCurveOnFace1");
        }
        self.my_int_build.pcurve_on_face1(index)
    }

    /// OCCT Builder.cxx L223-230 — PCurve1OnFillet (raises
    /// Standard_OutOfRange).
    pub fn pcurve1_on_fillet(&self, index: usize) -> Option<Curve2d> {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::PCurve1OnFillet");
        }
        self.my_int_build.pcurve1_on_fillet(index)
    }

    /// OCCT Builder.cxx L236-243 — PCurveOnFace2 (raises
    /// Standard_OutOfRange).
    pub fn pcurve_on_face2(&self, index: usize) -> Option<Curve2d> {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::PCurveOnFace2");
        }
        self.my_int_build.pcurve_on_face2(index)
    }

    /// OCCT Builder.cxx L249-256 — PCurve2OnFillet (raises
    /// Standard_OutOfRange).
    pub fn pcurve2_on_fillet(&self, index: usize) -> Option<Curve2d> {
        if index < 1 || index > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::PCurve2OnFillet");
        }
        self.my_int_build.pcurve2_on_fillet(index)
    }

    /// OCCT Builder.cxx L262-269 — FirstParameter (raises
    /// StdFail_NotDone).
    pub fn first_parameter(&self) -> f64 {
        if self.is_done() == FilletSurfStatusDone::IsNotOk {
            panic!("StdFail_NotDone: FilletSurf_Builder::FirstParameter");
        }
        self.my_int_build.first_parameter()
    }

    /// OCCT Builder.cxx L275-282 — LastParameter (raises StdFail_NotDone).
    pub fn last_parameter(&self) -> f64 {
        if self.is_done() == FilletSurfStatusDone::IsNotOk {
            panic!("StdFail_NotDone: FilletSurf_Builder::LastParameter");
        }
        self.my_int_build.last_parameter()
    }

    /// OCCT Builder.cxx L294-301 — StartSectionStatus (raises
    /// StdFail_NotDone).
    pub fn start_section_status(&self) -> FilletSurfStatusType {
        if self.is_done() == FilletSurfStatusDone::IsNotOk {
            panic!("StdFail_NotDone: FilletSurf_Builder::StartSectionStatus");
        }
        self.my_int_build.start_section_status()
    }

    /// OCCT Builder.cxx L313-320 — EndSectionStatus (raises
    /// StdFail_NotDone).
    pub fn end_section_status(&self) -> FilletSurfStatusType {
        if self.is_done() == FilletSurfStatusDone::IsNotOk {
            panic!("StdFail_NotDone: FilletSurf_Builder::StartSectionStatus");
        }
        self.my_int_build.end_section_status()
    }

    /// OCCT Builder.cxx L326-342 — Simulate.
    pub fn simulate(&mut self) {
        if self.myisdone == FilletSurfStatusDone::IsOk {
            self.my_int_build.simulate();

            if self.my_int_build.done() {
                self.myisdone = FilletSurfStatusDone::IsOk;
            } else {
                self.myisdone = FilletSurfStatusDone::IsNotOk;
                self.myerrorstatus = FilletSurfErrorTypeStatus::PbFilletCompute;
            }
        }
    }

    /// OCCT Builder.cxx L348-359 — NbSection (raises StdFail_NotDone /
    /// Standard_OutOfRange).
    pub fn nb_section(&self, index_surf: usize) -> usize {
        if self.is_done() == FilletSurfStatusDone::IsNotOk {
            panic!("StdFail_NotDone: FilletSurf_Builder::NbSection)");
        } else if index_surf < 1 || index_surf > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::NbSection");
        }
        self.my_int_build.nb_section(index_surf)
    }

    /// OCCT Builder.cxx L367-380 — Section (raises Standard_OutOfRange).
    pub fn section(&self, index_surf: usize, index_sec: usize) -> TrimmedCurve3 {
        if index_surf < 1 || index_surf > self.nb_surface() {
            panic!("Standard_OutOfRange: FilletSurf_Builder::Section NbSurface");
        } else if index_sec < 1 || index_sec > self.nb_section(index_surf) {
            panic!("Standard_OutOfRange: FilletSurf_Builder::Section NbSection");
        }

        self.my_int_build.section(index_surf, index_sec)
    }

    // OCCT Builder.cxx L382-387 — the deprecated Section overload taking
    // `handle<Geom_TrimmedCurve>& Circ` as an out-parameter collapses into
    // the by-value `section` above (no Rust out-parameters).
}
