//! OCCT ChFi3d_FilBuilder.cxx — the SimulSurf / PerformSurf overload family
//! (E3-S queue 3, the pending stand-ins of chfi3d_builder_2 / _2b).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_FilBuilder.cxx.
//!
//! Coverage of this file:
//!   - SimulParams            (OCCT L112-143, file static)
//!   - SimulSurf  face/face   (OCCT L558-784,  bool return)
//!   - PerformSurf face/face  (OCCT L1538-1717, bool return)
//!   - SplitSurf              (OCCT L2274-2436)
//!
//! The rst-carrying SimulSurf / PerformSurf overloads (face/rst, rst/face,
//! rst/rst — OCCT L788-1016 / L1020-1248 / L1252-1496 / L1721-1890 /
//! L1894-2064 / L2068-2270) live in [`super::chfi3d_builder_2d_b`] (split
//! off by the 2000-line rule).
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
//!   - (closed) BRepBlend_SurfRstEvolRad
//!     (brep_blend_surf_rst_evol_rad.rs), RstRstConstRad
//!     (brep_blend_rst_rst_const_rad.rs) and RstRstEvolRad
//!     (brep_blend_rst_rst_evol_rad.rs) — translated 1:1; every rst/rst and
//!     rst/face SimulSurf / PerformSurf branch (constant and variable) calls
//!     the real SimulData / ComputeData / CompleteData statements.
//!   - (closed) BRepBlend_SurfCurvEvolRadInv
//!     (brep_blend_surf_curv_evol_rad_inv.rs) and SurfPointEvolRadInv
//!     (brep_blend_surf_point_evol_rad_inv.rs) — translated 1:1.
//!   - BRepBlend_CurvPointRadInv bound to an Adaptor3d_CurveOnSurface: the
//!     plain-curve port (brep_blend_curv_point_rad_inv.rs) models curv2 as
//!     &Curve3; the HC-bound call sites use the HC-payload port
//!     brep_blend_curv_point_rad_inv_hc.rs whose GetTolerance keeps the
//!     OCCT-pending Adaptor3d_Curve::Resolution / Adaptor2d_Curve2d::
//!     Resolution failure path.
//!   - (closed) BlendFunc_EvolRadInv — translated 1:1 in
//!     blend_func_evol_rad_inv.rs; both 2-face EvolRad arms call the real
//!     SimulData / ComputeData statements.
//!   - (closed) BRepBlend_SurfRstConstRad (brep_blend_surf_rst_const_rad.rs),
//!     BRepBlend_SurfCurvConstRadInv (brep_blend_surf_curv_const_rad_inv.rs)
//!     and BRepBlend_SurfPointConstRadInv
//!     (brep_blend_surf_point_const_rad_inv.rs) — the face/rst constant
//!     arms call the real SimulData / ComputeData / CompleteData statements.
//!   - (closed) SplitSurf — the Geom_Surface::UIso re-host
//!     (crate::brep_fill::brep_fill_sweep::surface_uiso, now pub(crate)) and
//!     the bounded math_FunctionRoot constructor
//!     (rcad_kernel::math::newton_function_root::NewtonFunctionRoot::new_bounded)
//!     unblocked the 1:1 body; ChFi3d_SearchSing lives in
//!     chfi3d_search_sing.rs.

use rcad_kernel::core::precision::p_confusion;
use rcad_kernel::math::newton_function_root::NewtonFunctionRoot;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::Shape;
use rcad_kernel::geom::{Circle3, CurveEval as _, SurfaceEval as _};
use rcad_kernel::math::root::FunctionValue as _;
use glam::DVec2;

use crate::brep_fill::brep_fill_sweep::surface_uiso;
use crate::geomalgo::law::law_function::LawFunction;

use super::brep_blend_func_consrad::{BlendFuncConstRad, BlendFuncConstRadInv};
use super::brep_blend_func_evolrad::BlendFuncEvolRad;
use super::blend_func_evol_rad_inv::BlendFuncEvolRadInv;
use super::brep_blend_line::BRepBlendLine;
use super::chfi3d::ChFi3dBuilder;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_2::BRepTopAdaptorTopolTool;
use super::chfi3d_builder_6::chfi3d_fil_common_point;
use super::chfi3d_builder_6b::{elspine_guide_curve, ChFiDSElSpineHandle};
use super::chfi3d_search_sing::ChFi3dSearchSing;
use super::chfi_ds::{
    ChFiDSCircSection, ChFiDSCircSectionArray, ChFiDSElSpine,
    ChFiDSFilSpine, ChFiDSSpineHandle, ChFiDSSurfData, SharedSurfData,
};

/// OCCT ChFi3d_FilBuilder.cxx L112-143 — SimulParams (file static): the
/// flexible walking parameters (MaxStep / Fleche) for a simulation.
pub(crate) fn simul_params(
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
            let mut func = BlendFuncEvolRad::new(&s1.surface, &s2.surface, &guide, law.clone());
            let mut finv = BlendFuncEvolRadInv::new(&s1.surface, &s2.surface, &guide, law);
            func.set(choix);
            finv.set(choix);
            func.set_section_shape(self.my_blend_shape);
            // OCCT L631 form: Soldep — the rcad call site carries the 4-slot
            // array; the math_Vector form is materialized for the call.
            let mut soldep_v = Vector::new(1, 4);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);
            soldep_v.set(4, soldep[3]);
            // OCCT L670-692: done = SimulData(Data, HGuide, EmptyHGuide, lin,
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
            // OCCT L693-696.
            if !self.done {
                return false;
            }
            // OCCT L697-720: the section sampling loop.
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
            // fsp->Law(HGuide)); BRepBlend_EvolRadInv FInv(S1, S2, HGuide,
            // fsp->Law(HGuide)); Func.Set(Choix); FInv.Set(Choix);
            // Func.Set(myShape).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            let mut func = BlendFuncEvolRad::new(&s1.surface, &s2.surface, &guide, law.clone());
            let mut finv = BlendFuncEvolRadInv::new(&s1.surface, &s2.surface, &guide, law);
            func.set(choix);
            finv.set(choix);
            func.set_section_shape(self.my_blend_shape);
            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 4);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);
            soldep_v.set(4, soldep[3]);
            // OCCT L1659-1686: done = ComputeData(Data, HGuide, Spine, lin,
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
            // OCCT L1691-1694: recovery is possible PMN 14/05/1998.
            if !self.done {
                return false;
            }
            // OCCT L1700-1710.
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

    /// OCCT ChFi3d_FilBuilder.cxx L2274-2436 — SplitSurf (the near-singular
    /// sequence split after ComputeData).
    pub fn split_surf(&mut self, seqsd: &mut Vec<SharedSurfData>, line: &BRepBlendLine) {
        // OCCT L2277: int ii, Nbpnt = Line->NbPoints();
        let nbpnt = line.nb_points();
        // OCCT L2278-2281.
        if nbpnt < 3 {
            return;
        }
        // OCCT L2283: TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
        // OCCT L2285: occ::handle<ChFiDS_SurfData> ref = SeqData(1);
        let ref_data = seqsd[0].clone();
        // OCCT L2288-2292: ISurf = ref->Surf(); Surf = DStr.Surface(ISurf).Surface();
        // Surf->Bounds(UFirst, ULast, VFirst, VLast); Courbe1 = Surf->UIso(UFirst);
        // Courbe2 = Surf->UIso(ULast).  (The kernel Geom_Surface::UIso
        // dispatch is re-hosted by brep_fill_sweep::surface_uiso.)
        // OCCT keeps UFirst / ULast live for the two UIso calls only; the
        // rcad underscore prefix records that they have no later read.
        let (courbe1, courbe2, _u_first, _u_last, _v_first_bound, _v_last_bound) = {
            let dstr = self.my_ds.as_ref().expect("DS");
            let i_surf = ref_data.read().expect("surfdata lock").surf();
            let surf = dstr.surface(i_surf).surface();
            let d = surf.default_domain(); // OCCT Surf->Bounds(...)
            (
                surface_uiso(surf, d[0]),
                surface_uiso(surf, d[1]),
                d[0],
                d[1],
                d[2],
                d[3],
            )
        };
        // OCCT L2293: ChFi3d_SearchSing Fonc(Courbe1, Courbe2);
        let mut fonc = ChFi3dSearchSing::new(&courbe1, &courbe2);

        // OCCT L2295-2297: NCollection_Sequence<double> LesVi;
        // double precedant, suivant, courant; double a, b, c;
        let mut les_vi: Vec<f64> = Vec::new();
        let mut precedant;
        let mut suivant;
        let mut courant;
        let mut a;
        let mut b;
        let mut c;

        // (1) Finds vi so that iso v=vi is punctual
        // OCCT L2300-2303: VFirst / VLast from the interferences.
        let v_first = ref_data
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .parameter_first()
            .min(
                ref_data
                    .read()
                    .expect("surfdata lock")
                    .interference_on_s2()
                    .parameter_first(),
            );
        let v_last = ref_data
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .parameter_last()
            .max(
                ref_data
                    .read()
                    .expect("surfdata lock")
                    .interference_on_s2()
                    .parameter_last(),
            );

        // (1.1) Finds the first point inside
        // OCCT L2306-2308: for (ii = 1; ii <= Nbpnt &&
        // Line->Point(ii).Parameter() < VFirst; ii++) {}
        let mut ii = 1i32;
        while ii <= nbpnt && line.point(ii).parameter() < v_first {
            ii += 1;
        }
        // OCCT L2309-2312.
        if ii == 1 {
            ii += 1;
        }
        // OCCT L2313-2315: P = Line->Point(ii); b = P.Parameter();
        // courant = P.PointOnS1().Distance(P.PointOnS2());
        let p = line.point(ii);
        b = p.parameter();
        courant = p.point_on_s1().distance(p.point_on_s2());
        // OCCT L2316-2318: P = Line->Point(ii - 1); a = P.Parameter();
        // precedant = P.PointOnS1().Distance(P.PointOnS2());
        let p = line.point(ii - 1);
        a = p.parameter();
        precedant = p.point_on_s1().distance(p.point_on_s2());
        // OCCT L2319.
        ii += 1;

        // (1.2) Find a minimum by "points"
        // OCCT L2322-2378.
        while ii <= nbpnt && line.point(ii).parameter() <= v_last {
            // OCCT L2324-2328: the duplication skip.
            while ii <= nbpnt
                && line.point(ii).parameter() < v_last
                && line.point(ii).parameter() - b < p_confusion()
            {
                ii += 1;
            }

            // OCCT L2330-2332.
            let pnt = line.point(ii);
            c = pnt.parameter();
            suivant = pnt.point_on_s1().distance(pnt.point_on_s2());
            // OCCT L2333.
            if (courant < precedant) && (courant < suivant) {
                // (1.3) Find the exact minimum
                // OCCT L2336-2341: math_FunctionRoot Resol(Fonc, (a + c) / 2,
                // tol2d, a, c, 50) — the bounded math_FunctionRoot
                // constructor.
                let resol = NewtonFunctionRoot::new_bounded(
                    &mut fonc,
                    (a + c) / 2.0,
                    self.tol2d,
                    a,
                    c,
                    50,
                );
                // OCCT L2342.
                if resol.is_done() {
                    // OCCT L2344: double Val, racine = Resol.Root();
                    let racine = resol.root();
                    // OCCT L2346: Fonc.Value(Resol.Root(), Val);
                    let mut val = 0.0;
                    if let Some(v) = fonc.value(racine) {
                        val = v;
                    }
                    // OCCT L2347.
                    if val < self.tolapp3d {
                        // the solution (avoiding the risks of confusion)
                        // OCCT L2350-2363.
                        if les_vi.is_empty() {
                            if (racine > v_first + self.tol2d) && (racine < v_last - self.tol2d) {
                                les_vi.push(racine);
                            }
                        } else if (racine > les_vi[les_vi.len() - 1] + self.tol2d)
                            && (racine < v_last - self.tol2d)
                        {
                            les_vi.push(racine);
                        }
                    }
                } else {
                    // OCCT L2366-2371: the CHFI3D_DEB trace is a no-op here.
                }
            }
            // OCCT L2373-2377: update if non duplication.
            a = b;
            precedant = courant;
            b = c;
            courant = suivant;
            // OCCT L2322: the for-loop increment.
            ii += 1;
        }

        // (2) Update of the sequence of SurfData
        // OCCT L2381-2435.
        if !les_vi.is_empty() {
            // OCCT L2383: TopOpeBRepDS_DataStructure& DStru = myDS->ChangeDS();
            let dstru = self.my_ds.as_mut().expect("DS");
            // OCCT L2390: for (ii = 1; ii <= LesVi.Length(); ii++)
            for idx in 0..les_vi.len() {
                let ii = idx as i32 + 1;
                // OCCT L2393: T = LesVi(ii);
                let t = les_vi[idx];
                // (2.0) copy and insertion
                // OCCT L2395-2396: SD = new (ChFiDS_SurfData); SD->Copy(ref);
                let mut sd = ChFiDSSurfData::default();
                sd.copy(&ref_data.read().expect("surfdata lock"));
                // OCCT L2397: SeqData.InsertBefore(ii, SD);
                seqsd.insert(
                    (ii - 1) as usize,
                    std::sync::Arc::new(std::sync::RwLock::new(sd)),
                );
                // OCCT L2398-2399: S = DStru.Surface(ref->Surf());
                // SD->ChangeSurf(DStru.AddSurface(S));
                let s = dstru
                    .surface(ref_data.read().expect("surfdata lock").surf())
                    .clone();
                let new_surf_index = dstru.add_surface(s);
                seqsd[(ii - 1) as usize]
                    .write()
                    .expect("surfdata lock")
                    .change_surf(new_surf_index);
                // OCCT L2400-2401: C1 = DStru.Curve(SD->InterferenceOnS1().LineIndex());
                // SD->ChangeInterferenceOnS1().SetLineIndex(DStru.AddCurve(C1));
                // (C1 is a handle copy in OCCT; the rcad port is an owned
                // record, so the tolerance read of L2414 is captured here.)
                let c1 = dstru
                    .curve(
                        seqsd[(ii - 1) as usize]
                            .read()
                            .expect("surfdata lock")
                            .interference_on_s1()
                            .line_index(),
                    )
                    .clone();
                let c1_tolerance = c1.tolerance();
                let new_c1_index = dstru.add_curve(c1);
                seqsd[(ii - 1) as usize]
                    .write()
                    .expect("surfdata lock")
                    .change_interference_on_s1()
                    .set_line_index(new_c1_index);
                // OCCT L2402-2403: C2 = DStru.Curve(SD->InterferenceOnS2().LineIndex());
                // SD->ChangeInterferenceOnS2().SetLineIndex(DStru.AddCurve(C2));
                let c2 = dstru
                    .curve(
                        seqsd[(ii - 1) as usize]
                            .read()
                            .expect("surfdata lock")
                            .interference_on_s2()
                            .line_index(),
                    )
                    .clone();
                let c2_tolerance = c2.tolerance();
                let new_c2_index = dstru.add_curve(c2);
                seqsd[(ii - 1) as usize]
                    .write()
                    .expect("surfdata lock")
                    .change_interference_on_s2()
                    .set_line_index(new_c2_index);

                // (2.1) Modification of common Point
                // OCCT L2406-2409: the resets.
                {
                    let mut sdw = seqsd[(ii - 1) as usize].write().expect("surfdata lock");
                    sdw.change_vertex_last_on_s1().reset();
                    sdw.change_vertex_last_on_s2().reset();
                }
                {
                    let mut rw = ref_data.write().expect("surfdata lock");
                    rw.change_vertex_first_on_s1().reset();
                    rw.change_vertex_first_on_s2().reset();
                }
                // OCCT L2410-2411: Courbe1->D0(T, P1); Courbe2->D0(T, P2);
                let p1 = courbe1.point_at(t);
                let p2 = courbe2.point_at(t);
                // OCCT L2412: P3d.SetXYZ((P1.XYZ() + P2.XYZ()) / 2);
                let p3d = (p1 + p2) / 2.0;
                // OCCT L2413-2414: VertexTol = P1.Distance(P2);
                // VertexTol += std::max(C1.Tolerance(), C2.Tolerance());
                let mut vertex_tol = p1.distance(p2);
                vertex_tol += c1_tolerance.max(c2_tolerance);

                // OCCT L2416-2423: the point / tolerance loads.
                {
                    let mut sdw = seqsd[(ii - 1) as usize].write().expect("surfdata lock");
                    sdw.change_vertex_last_on_s1().set_point(p3d);
                    sdw.change_vertex_last_on_s2().set_point(p3d);
                    sdw.change_vertex_last_on_s1().set_tolerance(vertex_tol);
                    sdw.change_vertex_last_on_s2().set_tolerance(vertex_tol);
                }
                {
                    let mut rw = ref_data.write().expect("surfdata lock");
                    rw.change_vertex_first_on_s1().set_point(p3d);
                    rw.change_vertex_first_on_s2().set_point(p3d);
                    rw.change_vertex_first_on_s1().set_tolerance(vertex_tol);
                    rw.change_vertex_first_on_s2().set_tolerance(vertex_tol);
                }

                // (2.2) Modification of interferences
                // OCCT L2426-2429.
                {
                    let mut sdw = seqsd[(ii - 1) as usize].write().expect("surfdata lock");
                    sdw.change_interference_on_s1().set_last_parameter(t);
                    sdw.change_interference_on_s2().set_last_parameter(t);
                }
                {
                    let mut rw = ref_data.write().expect("surfdata lock");
                    rw.change_interference_on_s1().set_first_parameter(t);
                    rw.change_interference_on_s2().set_first_parameter(t);
                }

                // Parameters on ElSpine
                // OCCT L2432-2433: SD->LastSpineParam(T); ref->FirstSpineParam(T);
                seqsd[(ii - 1) as usize]
                    .write()
                    .expect("surfdata lock")
                    .set_last_spine_param(t);
                ref_data
                    .write()
                    .expect("surfdata lock")
                    .set_first_spine_param(t);
            }
        }
    }
}

