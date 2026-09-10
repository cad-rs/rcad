//! OCCT ChFi3d_FilBuilder.cxx — the SimulSurf / PerformSurf overload family
//! (E3-S queue 3, the pending stand-ins of chfi3d_builder_2 / _2b).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_FilBuilder.cxx.
//!
//! Coverage of this file:
//!   - SimulParams            (OCCT L112-143, file static)
//!   - SimulSurf  face/face   (OCCT L558-784,  bool return)
//!   - SimulSurf  face/rst    (OCCT L788-1016, void)
//!   - SimulSurf  rst/face    (OCCT L1020-1248, void)
//!   - SimulSurf  rst/rst     (OCCT L1252-1496, void)
//!   - PerformSurf face/face  (OCCT L1538-1717, bool return)
//!   - PerformSurf face/rst   (OCCT L1721-1890, void)
//!   - PerformSurf rst/face   (OCCT L1894-2064, void)
//!   - PerformSurf rst/rst    (OCCT L2068-2270, void)
//!   - SplitSurf              (OCCT L2274-2436)
//!
//! rcad architecture notes: the OCCT methods are ChFi3d_FilBuilder members
//! reached through the ChFi3d_Builder virtual dispatch of CallPerformSurf /
//! PerformSetOfSurfOnElSpine.  Over the rcad composition model
//! (ChFi3dFilBuilder embeds ChFi3dBuilder) the base-level callers cannot
//! name the derived type, so the bodies are modeled on the shared
//! ChFi3dBuilder — they read only state the base carries (my_blend_shape is
//! the OCCT ChFi3d_FilBuilder.hxx L377 `BlendFunc_SectionShape myShape`
//! member, mirrored onto the rcad base by SetFilletShape) and the spine
//! down-cast is performed inside each body exactly as in OCCT.
//!
//! GAPs (each preserves the OCCT failure path, named in place):
//!   - BlendFunc_EvolRadInv (TKFillet/BlendFunc/BlendFunc_EvolRadInv.cxx
//!     L22-470) — the variable-radius SimulData/ComputeData statements take
//!     the OCCT !done route.
//!   - BRepBlend_SurfRstConstRad / SurfRstEvolRad (BRepBlend_SurfRstConstRad.cxx
//!     L1-1048), RstRstConstRad / RstRstEvolRad (BRepBlend_RstRstConstRad.cxx
//!     L1-981), SurfCurvConstRadInv / SurfCurvEvolRadInv
//!     (BRepBlend_SurfCurvConstRadInv.cxx L1-284), SurfPointConstRadInv /
//!     SurfPointEvolRadInv (BRepBlend_SurfPointConstRadInv.cxx L1-276) — the
//!     rst-carrying SimulSurf / PerformSurf overloads take the OCCT
//!     Standard_Failure route.
//!   - SplitSurf needs Geom_Surface::UIso over the stored blend surface (the
//!     kernel has no iso extraction; the existing re-host surface_uiso is
//!     pub(super) to brep_fill) and the bounded math_FunctionRoot
//!     constructor; the call site is kept, the body is the recorded blocker.

use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::{Orientation, Shape};
use rcad_kernel::geom::Circle3;
use glam::DVec2;

use crate::geomalgo::law::law_function::LawFunction;

use super::brep_blend_func_consrad::{BlendFuncConstRad, BlendFuncConstRadInv};
use super::brep_blend_func_evolrad::BlendFuncEvolRad;
use super::brep_blend_line::BRepBlendLine;
use super::chfi3d::ChFi3dBuilder;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_2::{BRepAdaptorCurve2d, BRepTopAdaptorTopolTool};
use super::chfi3d_builder_6::chfi3d_fil_common_point;
use super::chfi3d_builder_6b::{elspine_guide_curve, ChFiDSElSpineHandle};
use super::chfi_ds::{
    ChFiDSCircSection, ChFiDSCircSectionArray, ChFiDSElSpine, ChFiDS_ErrorStatus,
    ChFiDSFilSpine, ChFiDSSpineHandle, SharedSurfData,
};

/// OCCT ChFi3d_FilBuilder.cxx L112-143 — SimulParams (file static): the
/// flexible walking parameters (MaxStep / Fleche) for a simulation.
fn simul_params(
    hguide: &ChFiDSElSpine,
    fsp: &ChFiDSFilSpine,
    max_step: &mut f64,
    fleche: &mut f64,
) {
    let la = hguide.last_parameter();
    let fi = hguide.first_parameter();
    let longueur = la - fi;
    *max_step = longueur * 0.05;
    let mut w;
    let radiussect;
    if fsp.is_constant() {
        radiussect = fsp.radius();
    } else {
        let mut radiussect_max = 0.0f64;
        let lc = fsp
            .law_of(hguide)
            .expect("ChFiDS_FilSpine::Law: no law for this elspine");
        for i in 0..=5i32 {
            w = fi + i as f64 * longueur * 0.2;
            let temp = lc.borrow_mut().value(w);
            if temp > radiussect_max {
                radiussect_max = temp;
            }
        }
        radiussect = radiussect_max;
    }
    *fleche = radiussect * 0.05;
}

impl ChFi3dBuilder {
    /// OCCT ChFi3d_FilBuilder.cxx L558-784 — SimulSurf (the face/face
    /// simulation entry; bool return, CallPerformSurf reprises on false).
    #[allow(clippy::too_many_arguments, unreachable_code)]
    pub fn simul_surf(
        &mut self,
        data: &SharedSurfData,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &[f64; 4],
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        // OCCT L578-582: fsp = down_cast<ChFiDS_FilSpine>(Spine); null ->
        // Standard_ConstructionError.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: SimulSurf : this is not the spine of the fillet"
            ),
        };
        // OCCT L583: occ::handle<BRepBlend_Line> lin;
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT L588-589: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L590-591: sec (null until a branch fills it); pf1..pl2.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        let (mut pf1, mut pl1, mut pf2, mut pl2) =
            (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        // OCCT L593: occ::handle<ChFiDS_ElSpine> EmptyHGuide.  The rcad
        // ElSpine flows by value (no handle identity) — the null additional
        // guide is passed as None.
        let empty_hguide: Option<ChFiDSElSpineHandle> = None;
        // OCCT L595: double PFirst = First;
        let p_first = *first;
        // OCCT L596-603: the intf/intl re-framing.
        if *intf != 0 {
            *first = fsp.base.first_parameter_of(1);
        }
        if *intl != 0 {
            *last = fsp.base.last_parameter_of(fsp.base.nb_edges());
        }
        // OCCT L604: if (fsp->IsConstant()).
        if fsp.is_constant() {
            // OCCT L606-610: BRepBlend_ConstRad Func(S1, S2, HGuide);
            // BRepBlend_ConstRadInv FInv(S1, S2, HGuide);
            // Func.Set(fsp->Radius(), Choix); FInv.Set(fsp->Radius(), Choix);
            // Func.Set(myShape).
            //
            // The rcad ElSpine is carried by value (architecture note); the
            // OCCT handle passed to the functions is materialized here.
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let mut func = BlendFuncConstRad::new(&s1.surface, &s2.surface, &guide);
            let mut finv = BlendFuncConstRadInv::new(&s1.surface, &s2.surface, &guide);
            func.set(fsp.radius(), choix);
            finv.set(fsp.radius(), choix);
            func.set_section_shape(self.my_blend_shape);
            // OCCT L631: Soldep — the rcad call site carries the 4-slot
            // array; the math_Vector form is materialized for the call.
            let mut soldep_v = Vector::new(1, 4);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);
            soldep_v.set(4, soldep[3]);
            // OCCT L611-633: done = SimulData(Data, HGuide, EmptyHGuide, lin,
            // S1, I1, S2, I2, Func, FInv, PFirst, MaxStep, locfleche,
            // TolGuide, First, Last, Inside, Appro, Forward, Soldep, 4,
            // RecOnS1, RecOnS2).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.simul_data_walking(
                &mut dw,
                &hguide_handle,
                empty_hguide.as_ref(),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                &mut func,
                &mut finv,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                &soldep_v,
                4,
                rec_on_s1,
                rec_on_s2,
            );
            // OCCT L634-637.
            if !self.done {
                return false;
            }
            let lin_ref = lin.as_ref().expect("Lin");
            // OCCT L638-661: the section sampling loop.
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u1, mut v1, mut u2, mut v2, mut w, mut p1, mut p2) =
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                (u1, v1) = p.parameters_on_s1();
                (u2, v2) = p.parameters_on_s2();
                w = p.parameter();
                func.section(w, u1, v1, u2, v2, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                if i == 1 {
                    pf1 = DVec2::new(u1, v1);
                    pf2 = DVec2::new(u2, v2);
                }
                if i == nbp {
                    pl1 = DVec2::new(u1, v1);
                    pl2 = DVec2::new(u2, v2);
                }
            }
            sec = Some(arr);
        } else {
            // OCCT L665-669: BRepBlend_EvolRad Func(S1, S2, HGuide,
            // fsp->Law(HGuide)); BRepBlend_EvolRadInv FInv(S1, S2, HGuide,
            // fsp->Law(HGuide)); Func.Set(Choix); FInv.Set(Choix);
            // Func.Set(myShape).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            let mut func = BlendFuncEvolRad::new(&s1.surface, &s2.surface, &guide, law);
            func.set(choix);
            func.set_section_shape(self.my_blend_shape);
            // OCCT L666: BRepBlend_EvolRadInv FInv(S1, S2, HGuide,
            // fsp->Law(HGuide)) — GAP: BlendFunc_EvolRadInv is untranslated
            // (TKFillet/BlendFunc/BlendFunc_EvolRadInv.cxx L22-470); the
            // SimulData statement (OCCT L670-692) cannot be carried.  The
            // OCCT !done route (L693-696) is taken.
            let _ = &func;
            self.done = false;
            if !self.done {
                return false;
            }
            // OCCT L697-720: the section sampling loop (unreachable while
            // the EvolRadInv GAP stands; kept in OCCT statement order).
            let lin_ref = lin.as_ref().expect("Lin");
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u1, mut v1, mut u2, mut v2, mut w, mut p1, mut p2) =
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                (u1, v1) = p.parameters_on_s1();
                (u2, v2) = p.parameters_on_s2();
                w = p.parameter();
                func.section(w, u1, v1, u2, v2, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                if i == 1 {
                    pf1 = DVec2::new(u1, v1);
                    pf2 = DVec2::new(u2, v2);
                }
                if i == nbp {
                    pl1 = DVec2::new(u1, v1);
                    pl2 = DVec2::new(u2, v2);
                }
            }
            sec = Some(arr);
        }
        // OCCT L722-743: Data->SetSimul / Set2dPoints / the four
        // ChFi3d_FilCommonPoint loads.
        {
            let mut dw = data.write().expect("surfdata lock");
            dw.set_simul(sec);
            dw.set_2d_points(pf1, pl1, pf2, pl2);
            let lin_ref = lin.as_ref().expect("Lin");
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.start_point_on_first(),
                lin_ref.transition_on_s1(),
                true,
                dw.change_vertex_first_on_s1(),
                self.tolapp3d,
            );
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.end_point_on_first(),
                lin_ref.transition_on_s1(),
                false,
                dw.change_vertex_last_on_s1(),
                self.tolapp3d,
            );
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.start_point_on_second(),
                lin_ref.transition_on_s2(),
                true,
                dw.change_vertex_first_on_s2(),
                self.tolapp3d,
            );
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.end_point_on_second(),
                lin_ref.transition_on_s2(),
                false,
                dw.change_vertex_last_on_s2(),
                self.tolapp3d,
            );
            // OCCT L744-763: the intf SearchFace block.
            let reverse = !forward || inside;
            if *intf != 0 && reverse {
                let mut ok = false;
                let cp1 = dw.vertex_first_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intf != 0;
                }
                let cp2 = dw.vertex_first_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
            // OCCT L764-782: the intl SearchFace block.
            if *intl != 0 {
                let mut ok = false;
                let cp1 = dw.vertex_last_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intl != 0;
                }
                let cp2 = dw.vertex_last_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
        }
        // OCCT L783.
        true
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1538-1717 — PerformSurf (the face/face
    /// entry; bool return, CallPerformSurf reprises on false).
    #[allow(clippy::too_many_arguments, unreachable_code)]
    pub fn perform_surf(
        &mut self,
        seqsd: &mut Vec<SharedSurfData>,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
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
        soldep: &[f64; 4],
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        // OCCT L1563: Data = SeqData(1).
        let data = seqsd.first().cloned().expect("surfdata");
        // OCCT L1564-1568: fsp down-cast; null -> Standard_ConstructionError.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : this is not the spine of a fillet"
            ),
        };
        // OCCT L1569: bool gd1, gd2, gf1, gf2, maybesingular (uninitialized
        // declarations in OCCT).
        let (mut gd1, mut gd2, mut gf1, mut gf2) = (false, false, false, false);
        let mut maybesingular = false;
        // OCCT L1570: occ::handle<BRepBlend_Line> lin;
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT L1571: TopAbs_Orientation Or = S1->Face().Orientation();
        let or = s1.face.orientation;
        // OCCT L1572: double PFirst = First;
        let p_first = *first;
        // OCCT L1573-1580: the intf/intl re-framing.
        if *intf != 0 {
            *first = fsp.base.first_parameter_of(1);
        }
        if *intl != 0 {
            *last = fsp.base.last_parameter_of(fsp.base.nb_edges());
        }
        // OCCT L1581: if (fsp->IsConstant()).
        if fsp.is_constant() {
            // OCCT L1583-1587: BRepBlend_ConstRad Func(S1, S2, HGuide);
            // BRepBlend_ConstRadInv FInv(S1, S2, HGuide);
            // Func.Set(fsp->Radius(), Choix); FInv.Set(fsp->Radius(), Choix);
            // Func.Set(myShape).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let mut func = BlendFuncConstRad::new(&s1.surface, &s2.surface, &guide);
            let mut finv = BlendFuncConstRadInv::new(&s1.surface, &s2.surface, &guide);
            func.set(fsp.radius(), choix);
            finv.set(fsp.radius(), choix);
            func.set_section_shape(self.my_blend_shape);
            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 4);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);
            soldep_v.set(4, soldep[3]);
            // OCCT L1593-1620: done = ComputeData(Data, HGuide, Spine, lin,
            // S1, I1, S2, I2, Func, FInv, PFirst, MaxStep, Fleche, TolGuide,
            // First, Last, Inside, Appro, Forward, Soldep, intf, intl,
            // gd1, gd2, gf1, gf2, RecOnS1, RecOnS2).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.compute_data(
                &mut dw,
                &hguide_handle,
                Some(spine),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                &mut func,
                &mut finv,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                &soldep_v,
                intf,
                intl,
                &mut gd1,
                &mut gd2,
                &mut gf1,
                &mut gf2,
                rec_on_s1,
                rec_on_s2,
            );
            // OCCT L1626-1629: recovery is possible PMN 14/05/1998.
            if !self.done {
                return false;
            }
            // OCCT L1635: done = CompleteData(Data, Func, lin, S1, S2, Or,
            // gd1, gd2, gf1, gf2) — the L516 overload carries a trailing
            // Reversed flag (false at this call site).
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_function(
                &mut dw, &mut func, lin_ref, s1, Some(s2), or, gd1, gd2, gf1, gf2, false,
            );
            // OCCT L1641-1644.
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            // OCCT L1645.
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
        } else {
            // OCCT L1649-1653: BRepBlend_EvolRad Func(S1, S2, HGuide,
            // fsp->Law(HGuide)); Func.Set(Choix); Func.Set(myShape).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            let mut func = BlendFuncEvolRad::new(&s1.surface, &s2.surface, &guide, law);
            func.set(choix);
            func.set_section_shape(self.my_blend_shape);
            // OCCT L1650: BRepBlend_EvolRadInv FInv(S1, S2, HGuide,
            // fsp->Law(HGuide)) — GAP: BlendFunc_EvolRadInv is untranslated
            // (TKFillet/BlendFunc/BlendFunc_EvolRadInv.cxx L22-470); the
            // ComputeData statement (OCCT L1659-1686) cannot be carried.
            // The OCCT !done route (L1691-1694) is taken.
            let _ = &func;
            self.done = false;
            if !self.done {
                return false;
            }
            // OCCT L1700-1710 (unreachable while the EvolRadInv GAP stands;
            // kept in OCCT statement order).
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_function(
                &mut data.write().expect("surfdata lock"),
                &mut func,
                lin_ref,
                s1,
                Some(s2),
                or,
                gd1,
                gd2,
                gf1,
                gf2,
                false,
            );
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
        }
        // OCCT L1712-1715.
        if maybesingular {
            self.split_surf(seqsd, lin.as_ref().expect("Lin"));
        }
        // OCCT L1716.
        true
    }
}

// =========================================================================
// The rst-carrying SimulSurf / PerformSurf overloads (curve-on-surface
// obstacle on S1 / on S2 / curve-curve).  The blend function families they
// construct (BRepBlend_SurfRstConstRad / SurfRstEvolRad /
// RstRstConstRad / RstRstEvolRad / SurfCurv*Inv / SurfPoint*Inv) are
// untranslated (see the module GAP list); every body carries the OCCT
// Standard_Failure route at the point the missing construction blocks the
// SimulData / ComputeData statement.
// =========================================================================

impl ChFi3dBuilder {
    /// OCCT ChFi3d_FilBuilder.cxx L788-1016 — SimulSurf (face/rst: the
    /// obstacle curve lies on S1).
    #[allow(clippy::too_many_arguments)]
    pub fn simul_surf_face_rst(
        &mut self,
        data: &SharedSurfData,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        hs1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        pc1: &Option<BRepAdaptorCurve2d>,
        hsref1: &BRepAdaptorSurface,
        pcref1: &Option<BRepAdaptorCurve2d>,
        decroch1: &mut bool,
        hs2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        or2: Orientation,
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
        soldep: &[f64; 3],
    ) {
        // OCCT L813-817: fsp down-cast (the OCCT message text is kept
        // verbatim, including its "PerformSurf" prefix).
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : this is not the spine of the fillet"
            ),
        };
        // OCCT L818: occ::handle<BRepBlend_Line> lin;
        let lin: Option<BRepBlendLine> = None;
        // OCCT L821-822: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L823-824: sec; gp_Pnt2d pf, pl, ppcf, ppcl.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        let (mut pf, mut pl, mut ppcf, mut ppcl) =
            (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        let _ = (
            &mut sec, pf, pl, ppcf, ppcl, &mut max_step, locfleche, data, hs1, i1, pc1, hsref1,
            pcref1, decroch1, hs2, i2, or2, fleche, tol_guide, inside, appro, forward, rec_p,
            rec_s, rec_rst, soldep,
        );
        // OCCT L826: double PFirst = First;
        let _p_first = *first;
        if fsp.is_constant() {
            // OCCT L829-852: BRepBlend_SurfRstConstRad func(HS2, HS1, PC1,
            // HGuide); func.Set(HSref1, PCref1); HC->Load(PC1, HS1);
            // BRepBlend_SurfCurvConstRadInv finvc(HS2, HC, HGuide);
            // BRepBlend_SurfPointConstRadInv finvp(HS2, HGuide);
            // BRepBlend_ConstRadInv finv(HS2, HSref1, HGuide);
            // finv.Set(false, PCref1); rad/petitchoix; the Set calls;
            // func.Set(myShape) — GAP: BRepBlend_SurfRstConstRad
            // (BRepBlend_SurfRstConstRad.cxx L1-1048),
            // BRepBlend_SurfCurvConstRadInv (…SurfCurvConstRadInv.cxx L1-284)
            // and BRepBlend_SurfPointConstRadInv (…SurfPointConstRadInv.cxx
            // L1-276) are untranslated; the SimulData statement
            // (OCCT L854-880) cannot be carried.  The OCCT !done route
            // (L881-884) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Failed process!");
        } else {
            // OCCT L914-934: the SurfRstEvolRad branch — same GAP
            // (BRepBlend_SurfRstEvolRad / SurfCurvEvolRadInv /
            // SurfPointEvolRadInv); the SimulData statement
            // (OCCT L935-961) cannot be carried.  The OCCT !done route
            // (L962-965) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Fail !");
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1020-1248 — SimulSurf (rst/face: the
    /// obstacle curve lies on S2).
    #[allow(clippy::too_many_arguments)]
    pub fn simul_surf_rst_face(
        &mut self,
        data: &SharedSurfData,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        hs1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        or1: Orientation,
        hs2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        pc2: &Option<BRepAdaptorCurve2d>,
        hsref2: &BRepAdaptorSurface,
        pcref2: &Option<BRepAdaptorCurve2d>,
        decroch2: &mut bool,
        arrow: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &[f64; 3],
    ) {
        // OCCT L1045-1049: fsp down-cast.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : it is not the spine of a fillet"
            ),
        };
        // OCCT L1050: occ::handle<BRepBlend_Line> lin;
        let lin: Option<BRepBlendLine> = None;
        // OCCT L1053-1054: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L1055-1056: sec; gp_Pnt2d pf, pl, ppcf, ppcl.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        let (mut pf, mut pl, mut ppcf, mut ppcl) =
            (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        let _ = (
            &mut sec, pf, pl, ppcf, ppcl, &mut max_step, locfleche, data, hs1, i1, or1, hs2, i2,
            pc2, hsref2, pcref2, decroch2, arrow, tol_guide, inside, appro, forward, rec_p, rec_s,
            rec_rst, soldep,
        );
        // OCCT L1058: double PFirst = First;
        let _p_first = *first;
        if fsp.is_constant() {
            // OCCT L1061-1084: the SurfRstConstRad branch — GAP (see
            // simul_surf_face_rst for the missing-class anchors); the
            // SimulData statement (OCCT L1086-1112) cannot be carried.
            // The OCCT !done route (L1113-1116) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Failed Processing!");
        } else {
            // OCCT L1146-1166: the SurfRstEvolRad branch — same GAP; the
            // SimulData statement (OCCT L1167-1193) cannot be carried.
            // The OCCT !done route (L1194-1197) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Fail !");
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1252-1496 — SimulSurf (rst/rst: the
    /// curve-curve entry).
    #[allow(clippy::too_many_arguments)]
    pub fn simul_surf_rst_rst(
        &mut self,
        data: &SharedSurfData,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        hs1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        pc1: &Option<BRepAdaptorCurve2d>,
        hsref1: &BRepAdaptorSurface,
        pcref1: &Option<BRepAdaptorCurve2d>,
        decroch1: &mut bool,
        or1: Orientation,
        hs2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        pc2: &Option<BRepAdaptorCurve2d>,
        hsref2: &BRepAdaptorSurface,
        pcref2: &Option<BRepAdaptorCurve2d>,
        decroch2: &mut bool,
        or2: Orientation,
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
        soldep: &[f64; 2],
    ) {
        // OCCT L1283-1287: fsp down-cast.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : it is not the spine of a fillet"
            ),
        };
        // OCCT L1288: occ::handle<BRepBlend_Line> lin;
        let lin: Option<BRepBlendLine> = None;
        // OCCT L1291-1292: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L1293: sec.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        let _ = (
            &mut sec, &mut max_step, locfleche, data, hs1, i1, pc1, hsref1, pcref1, decroch1, or1,
            hs2, i2, pc2, hsref2, pcref2, decroch2, or2, fleche, tol_guide, inside, appro,
            forward, rec_p1, rec_rst1, rec_p2, rec_rst2, soldep,
        );
        // OCCT L1296-1297: int ch1 = 1, ch2 = 2; double PFirst = First;
        let _p_first = *first;
        if fsp.is_constant() {
            // OCCT L1301-1330: BRepBlend_RstRstConstRad func(HS1, PC1, HS2,
            // PC2, HGuide); func.Set(HSref1, PCref1, HSref2, PCref2);
            // HC1/HC2 -> Load; BRepBlend_SurfCurvConstRadInv finv1/finv2;
            // BRepBlend_CurvPointRadInv finvp1/finvp2; the Set calls;
            // func.Set(myShape) — GAP: BRepBlend_RstRstConstRad
            // (BRepBlend_RstRstConstRad.cxx L1-981) and
            // BRepBlend_SurfCurvConstRadInv (…SurfCurvConstRadInv.cxx
            // L1-284) are untranslated; the SimulData statement
            // (OCCT L1332-1362) cannot be carried.  The OCCT !done route
            // (L1363-1366) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Failed processing!");
        } else {
            // OCCT L1386-1417: the RstRstEvolRad branch — same GAP
            // (BRepBlend_RstRstEvolRad); the SimulData statement
            // (OCCT L1419-1449) cannot be carried.  The OCCT !done route
            // (L1451-1454) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Fail !");
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1721-1890 — PerformSurf (face/rst: the
    /// obstacle curve lies on S1).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf_face_rst(
        &mut self,
        seqsd: &mut Vec<SharedSurfData>,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        hs1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        pc1: &Option<BRepAdaptorCurve2d>,
        hsref1: &BRepAdaptorSurface,
        pcref1: &Option<BRepAdaptorCurve2d>,
        decroch1: &mut bool,
        hs2: &BRepAdaptorSurface,
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
        soldep: &[f64; 3],
    ) {
        // OCCT L1747: Data = SeqData(1).
        let data = seqsd.first().cloned().expect("surfdata");
        // OCCT L1748-1752: fsp down-cast.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : this is not the spine of a fillet"
            ),
        };
        // OCCT L1753-1755: lin; PFirst; maybesingular.
        let lin: Option<BRepBlendLine> = None;
        let mut maybesingular = false;
        // OCCT L1757: if (fsp->IsConstant()).
        let _ = (
            data, hs1, i1, pc1, hsref1, pcref1, decroch1, hs2, i2, or2, max_step, fleche,
            tol_guide, inside, appro, forward, rec_p, rec_s, rec_rst, soldep, &mut maybesingular,
            &lin,
        );
        let _p_first = *first; // OCCT L1754.
        if fsp.is_constant() {
            // OCCT L1759-1782: the SurfRstConstRad branch — GAP (see
            // simul_surf_face_rst for the missing-class anchors); the
            // ComputeData statement (OCCT L1784-1809) cannot be carried.
            // The OCCT !done route (L1810-1814) is preserved.
            self.done = false;
            {
                // OCCT L1812: Spine->SetErrorStatus(ChFiDS_WalkingFailure).
                // Boundary note: OCCT mutates through the const handle; the
                // rcad handle is an enum without interior mutability — the
                // status is set on the clone (lost at the throw below, as
                // the OCCT exception also unwinds the stripe state).
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
            }
            panic!("Standard_Failure: PerformSurf : Failed processing!");
        } else {
            // OCCT L1825-1846: the SurfRstEvolRad branch — same GAP; the
            // ComputeData statement (OCCT L1847-1872) cannot be carried.
            // The OCCT !done route (L1873-1877) is preserved.
            self.done = false;
            {
                // OCCT L1875.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
            }
            panic!("Standard_Failure: PerformSurf : Failed processing!");
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1894-2064 — PerformSurf (rst/face: the
    /// obstacle curve lies on S2).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf_rst_face(
        &mut self,
        seqsd: &mut Vec<SharedSurfData>,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        hs1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        or1: Orientation,
        hs2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        pc2: &Option<BRepAdaptorCurve2d>,
        hsref2: &BRepAdaptorSurface,
        pcref2: &Option<BRepAdaptorCurve2d>,
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
        soldep: &[f64; 3],
    ) {
        // OCCT L1920: Data = SeqData(1).
        let data = seqsd.first().cloned().expect("surfdata");
        // OCCT L1921-1925: fsp down-cast.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : this is not the spine of a fillet"
            ),
        };
        // OCCT L1926-1928: lin; PFirst; maybesingular.
        let lin: Option<BRepBlendLine> = None;
        let mut maybesingular = false;
        let _ = (
            data, hs1, i1, or1, hs2, i2, pc2, hsref2, pcref2, decroch2, max_step, fleche,
            tol_guide, inside, appro, forward, rec_p, rec_s, rec_rst, soldep, &mut maybesingular,
            &lin,
        );
        let _p_first = *first; // OCCT L1927.
        if fsp.is_constant() {
            // OCCT L1932-1955: the SurfRstConstRad branch — GAP (see
            // simul_surf_face_rst for the missing-class anchors); the
            // ComputeData statement (OCCT L1957-1982) cannot be carried.
            // The OCCT !done route (L1983-1987) is preserved.
            self.done = false;
            {
                // OCCT L1985.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
            }
            panic!("Standard_Failure: PerformSurf : Failed processing!");
        } else {
            // OCCT L1998-2019: the SurfRstEvolRad branch — same GAP; the
            // ComputeData statement (OCCT L2021-2046) cannot be carried.
            // The OCCT !done route (L2047-2051) is preserved.
            self.done = false;
            {
                // OCCT L2049.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
            }
            panic!("Standard_Failure: PerformSurf : Failed processing!");
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L2068-2270 — PerformSurf (rst/rst: the
    /// curve-curve entry).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf_rst_rst(
        &mut self,
        seqsd: &mut Vec<SharedSurfData>,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        hs1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        pc1: &Option<BRepAdaptorCurve2d>,
        hsref1: &BRepAdaptorSurface,
        pcref1: &Option<BRepAdaptorCurve2d>,
        decroch1: &mut bool,
        or1: Orientation,
        hs2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        pc2: &Option<BRepAdaptorCurve2d>,
        hsref2: &BRepAdaptorSurface,
        pcref2: &Option<BRepAdaptorCurve2d>,
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
        soldep: &[f64; 2],
    ) {
        // OCCT L2100: Data = SeqData(1).
        let data = seqsd.first().cloned().expect("surfdata");
        // OCCT L2101-2105: fsp down-cast.
        let fsp = match spine.down_cast_fil() {
            Some(fsp) => fsp,
            None => panic!(
                "Standard_ConstructionError: PerformSurf : this is not the spine of a fillet"
            ),
        };
        // OCCT L2106-2108: lin; PFirst; maybesingular.
        let lin: Option<BRepBlendLine> = None;
        let mut maybesingular = false;
        let _ = (
            data, hs1, i1, pc1, hsref1, pcref1, decroch1, or1, hs2, i2, pc2, hsref2, pcref2,
            decroch2, or2, max_step, fleche, tol_guide, inside, appro, forward, rec_p1,
            rec_rst1, rec_p2, rec_rst2, soldep, &mut maybesingular, &lin,
        );
        let _p_first = *first; // OCCT L2107.
        if fsp.is_constant() {
            // OCCT L2112-2142: the RstRstConstRad branch — GAP (see
            // simul_surf_rst_rst for the missing-class anchors); the
            // ComputeData statement (OCCT L2144-2173) cannot be carried.
            // The OCCT !done route (L2174-2178) is preserved.
            self.done = false;
            {
                // OCCT L2176.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
            }
            panic!("Standard_Failure: PerformSurf : Failed processing!");
        } else {
            // OCCT L2189-2220: the RstRstEvolRad branch — same GAP; the
            // ComputeData statement (OCCT L2222-2251) cannot be carried.
            // The OCCT !done route (L2253-2257) is preserved.
            self.done = false;
            {
                // OCCT L2255.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
            }
            panic!("Standard_Failure: PerformSurf : Failed processing!");
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L2274-2436 — SplitSurf (the near-singular
    /// sequence split after ComputeData).  GAP carrier / BLOCKER: the body
    /// needs (1) Geom_Surface::UIso over the stored blend surface — the
    /// kernel has no iso extraction and the existing re-host surface_uiso
    /// is pub(super) to brep_fill — and (2) the bounded math_FunctionRoot
    /// constructor (OCCT L2336-2341; the rcad NewtonFunctionRoot carries
    /// only the full-range form).  The OCCT call site (maybesingular) is
    /// kept; until translated the split is skipped.  Recorded as a blocker
    /// in the E3-S session report.
    pub fn split_surf(&mut self, seqsd: &mut Vec<SharedSurfData>, line: &BRepBlendLine) {
        // OCCT L2277-2281: Nbpnt guard (kept: the early return is reachable
        // and preserves the OCCT no-op path for short lines).
        let nbpnt = line.nb_points();
        let _ = seqsd; // OCCT L2285: ref = SeqData(1) — inside the GAP.
        if nbpnt < 3 {
            return;
        }
        // OCCT L2282-2434: blocked by the UIso / bounded-FunctionRoot GAPs
        // (see above).
    }
}
