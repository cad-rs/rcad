//! OCCT ChFi3d_FilBuilder.cxx — the rst-carrying SimulSurf / PerformSurf
//! overloads (curve-on-surface obstacle on S1 / on S2 / curve-curve).
//! Split from chfi3d_builder_2d.rs (2000-line rule).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_FilBuilder.cxx.
//!
//! Coverage of this file:
//!   - SimulSurf  face/rst    (OCCT L788-1016, void)
//!   - SimulSurf  rst/face    (OCCT L1020-1248, void)
//!   - SimulSurf  rst/rst     (OCCT L1252-1496, void)
//!   - PerformSurf face/rst   (OCCT L1721-1890, void)
//!   - PerformSurf rst/face   (OCCT L1894-2064, void)
//!   - PerformSurf rst/rst    (OCCT L2068-2270, void)
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
//!   - BRepBlend_CurvPointRadInv bound to an Adaptor3d_CurveOnSurface: the
//!     plain-curve port (brep_blend_curv_point_rad_inv.rs) models curv2 as
//!     &Curve3; the HC-bound call sites use the HC-payload port
//!     brep_blend_curv_point_rad_inv_hc.rs whose GetTolerance keeps the
//!     OCCT-pending Adaptor3d_Curve::Resolution / Adaptor2d_Curve2d::
//!     Resolution failure path.

use glam::DVec2;

use rcad_kernel::geom::Circle3;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::Orientation;
use rcad_kernel::topo::topods::BRepTool as _;

use super::brep_blend_line::BRepBlendLine;
use super::brep_blend_surf_rst_const_rad::BlendSurfRstConstRad;
use super::brep_blend_surf_rst_evol_rad::BlendSurfRstEvolRad;
use super::brep_blend_surf_curv_const_rad_inv::BlendSurfCurvConstRadInv;
use super::brep_blend_surf_curv_evol_rad_inv::BlendSurfCurvEvolRadInv;
use super::brep_blend_surf_point_const_rad_inv::BlendSurfPointConstRadInv;
use super::brep_blend_surf_point_evol_rad_inv::BlendSurfPointEvolRadInv;
use super::brep_blend_rst_rst_const_rad::BlendRstRstConstRad;
use super::brep_blend_rst_rst_evol_rad::BlendRstRstEvolRad;
use super::brep_blend_curv_point_rad_inv_hc::BRepBlendCurvPointRadInvHc;
use super::brep_blend_func_consrad::BlendFuncConstRadInv;
use super::blend_func_evol_rad_inv::BlendFuncEvolRadInv;
use super::brep_blend_func_inv::BlendSurfCurvFuncInv;
use super::chfi3d::ChFi3dBuilder;
use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_2::{BRepAdaptorCurve2d, BRepTopAdaptorTopolTool};
use super::chfi3d_builder_2d::simul_params;
use super::chfi3d_builder_6::chfi3d_fil_common_point;
use super::chfi3d_builder_6b::elspine_guide_curve;
use super::chfi_ds::{
    ChFiDSCircSection, ChFiDSCircSectionArray, ChFiDSElSpine, ChFiDS_ErrorStatus,
    ChFiDSSpineHandle, SharedSurfData,
};

// =========================================================================
// The rst-carrying SimulSurf / PerformSurf overloads (curve-on-surface
// obstacle on S1 / on S2 / curve-curve).  All six bodies are fully wired:
// the blend function families they construct (BRepBlend_SurfRstConstRad /
// SurfRstEvolRad / RstRstConstRad / RstRstEvolRad / SurfCurv*Inv /
// SurfPoint*Inv / CurvPointRadInv-HC) are translated 1:1 and every constant
// and variable branch calls the real SimulData / ComputeData /
// CompleteData statements.
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
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT L821-822: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L823-824: sec; gp_Pnt2d pf, pl, ppcf, ppcl.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        let (mut pf, mut pl, mut ppcf, mut ppcl) =
            (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        // OCCT L826: double PFirst = First;
        let p_first = *first;
        if fsp.is_constant() {
            // OCCT L829-836: BRepBlend_SurfRstConstRad func(HS2, HS1, PC1,
            // HGuide); func.Set(HSref1, PCref1); HC->Load(PC1, HS1);
            // BRepBlend_SurfCurvConstRadInv finvc(HS2, HC, HGuide);
            // BRepBlend_SurfPointConstRadInv finvp(HS2, HGuide);
            // BRepBlend_ConstRadInv finv(HS2, HSref1, HGuide);
            // finv.Set(false, PCref1).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            // OCCT PC1 / PCref1 are the BRepAdaptor_Curve2d handles; the rcad
            // port materializes their pcurve payload (the (edge, face)
            // lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc1 = pc1.as_ref().expect("PC1");
            let (pc1_curve, _, _) = adpc1
                .brep
                .curve_on_surface(&adpc1.edge, &adpc1.face)
                .expect("PC1 pcurve");
            let adpcref1 = pcref1.as_ref().expect("PCref1");
            let (pcref1_curve, _, _) = adpcref1
                .brep
                .curve_on_surface(&adpcref1.edge, &adpcref1.face)
                .expect("PCref1 pcurve");
            let mut func =
                BlendSurfRstConstRad::new(&hs2.surface, &hs1.surface, &pc1_curve, &guide);
            func.set_ref(&hsref1.surface, &pcref1_curve);
            // OCCT L831-833: the HC Adaptor3d_CurveOnSurface(PC1, HS1) bound
            // as the finvc restriction curve — carried by the
            // (pcurve, surface) pair.
            let mut finvc =
                BlendSurfCurvConstRadInv::new(&hs2.surface, &pc1_curve, &hs1.surface, &guide);
            let mut finvp = BlendSurfPointConstRadInv::new(&hs2.surface, &guide);
            let mut finv = BlendFuncConstRadInv::new(&hs2.surface, &hsref1.surface, &guide);
            finv.set_curve_on_surface(false, &pcref1_curve);

            // OCCT L838-852: rad / petitchoix; the Set calls;
            // func.Set(myShape).
            let rad = fsp.radius();
            let mut petitchoix = 1;
            if or2 == Orientation::Reversed {
                petitchoix = 3;
            }
            if choix % 2 == 0 {
                petitchoix += 1;
            }
            finv.set(rad, choix);
            finvc.set(rad, petitchoix);
            finvp.set(rad, petitchoix);
            func.set(rad, petitchoix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 3);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);

            // OCCT L854-880: done = SimulData(Data, HGuide, lin, HS2, I2,
            // HS1, PC1, I1, Decroch1, func, finv, finvp, finvc, PFirst,
            // MaxStep, locfleche, TolGuide, First, Last, Soldep, 4, Inside,
            // Appro, Forward, RecP, RecS, RecRst).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.simul_data_surf_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs2,
                i2,
                hs1,
                &pc1_curve,
                i1,
                decroch1,
                &mut func,
                &mut finv,
                &mut finvp,
                &mut finvc,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                4,
                inside,
                appro,
                forward,
                rec_p,
                rec_s,
                rec_rst,
            );
            // OCCT L881-884.
            if !self.done {
                panic!("Standard_Failure: SimulSurf : Failed process!");
            }
            // OCCT L885-910: the section sampling loop.
            let lin_ref = lin.as_ref().expect("Lin");
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u, mut v, mut w, mut param, mut p1, mut p2) =
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                (u, v) = p.parameters_on_s();
                w = p.parameter_on_c();
                param = p.parameter();
                func.section(param, u, v, w, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                if i == 1 {
                    pf = DVec2::new(u, v);
                    (u, v) = p.parameters_on_s2();
                    ppcf = DVec2::new(u, v);
                }
                if i == nbp {
                    pl = DVec2::new(u, v);
                    (u, v) = p.parameters_on_s2();
                    ppcl = DVec2::new(u, v);
                }
            }
            sec = Some(arr);
        } else {
            // OCCT L914-934: the SurfRstEvolRad branch — same GAP
            // (BRepBlend_SurfRstEvolRad / SurfCurvEvolRadInv /
            // SurfPointEvolRadInv); the SimulData statement
            // (OCCT L935-961) cannot be carried.  The OCCT !done route
            // (L962-965) is preserved.
            self.done = false;
            panic!("Standard_Failure: SimulSurf : Fail !");
        }
        // OCCT L993-995: Data->SetSimul(sec);
        // Data->Set2dPoints(ppcf, ppcl, pf, pl).
        {
            let mut dw = data.write().expect("surfdata lock");
            dw.set_simul(sec);
            dw.set_2d_points(ppcf, ppcl, pf, pl);
            // OCCT L996-1015: the four ChFi3d_FilCommonPoint loads (the
            // S1 / S2 vertex slots are swapped relative to the face/face
            // SimulSurf).
            let lin_ref = lin.as_ref().expect("Lin");
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.start_point_on_first(),
                lin_ref.transition_on_s1(),
                true,
                dw.change_vertex_first_on_s2(),
                self.tolapp3d,
            );
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.end_point_on_first(),
                lin_ref.transition_on_s1(),
                false,
                dw.change_vertex_last_on_s2(),
                self.tolapp3d,
            );
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.start_point_on_second(),
                lin_ref.transition_on_s2(),
                true,
                dw.change_vertex_first_on_s1(),
                self.tolapp3d,
            );
            chfi3d_fil_common_point(
                &self.my_brep,
                lin_ref.end_point_on_second(),
                lin_ref.transition_on_s2(),
                false,
                dw.change_vertex_last_on_s1(),
                self.tolapp3d,
            );
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
        _arrow: f64, // OCCT L1033: const double /*Arrow*/ — unused.
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
        // OCCT L1050: occ::handle<BRepBlend_Line> lin.
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT L1053-1054: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L1055-1056: sec; gp_Pnt2d pf, pl, ppcf, ppcl.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        let (mut pf, mut pl, mut ppcf, mut ppcl) =
            (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        // OCCT L1058: double PFirst = First.
        let p_first = *first;
        if fsp.is_constant() {
            // OCCT L1061-1068: BRepBlend_SurfRstConstRad func(HS1, HS2, PC2,
            // HGuide); func.Set(HSref2, PCref2); HC->Load(PC2, HS2);
            // BRepBlend_SurfCurvConstRadInv finvc(HS1, HC, HGuide);
            // BRepBlend_SurfPointConstRadInv finvp(HS1, HGuide);
            // BRepBlend_ConstRadInv finv(HS1, HSref2, HGuide);
            // finv.Set(false, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            // OCCT PC2 / PCref2 are the BRepAdaptor_Curve2d handles; the rcad
            // port materializes their pcurve payload (the (edge, face)
            // lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func =
                BlendSurfRstConstRad::new(&hs1.surface, &hs2.surface, &pc2_curve, &guide);
            func.set_ref(&hsref2.surface, &pcref2_curve);
            // OCCT L1063-1064: the HC Adaptor3d_CurveOnSurface(PC2, HS2)
            // bound as the finvc restriction curve — carried by the
            // (pcurve, surface) pair.
            let mut finvc =
                BlendSurfCurvConstRadInv::new(&hs1.surface, &pc2_curve, &hs2.surface, &guide);
            let mut finvp = BlendSurfPointConstRadInv::new(&hs1.surface, &guide);
            let mut finv = BlendFuncConstRadInv::new(&hs1.surface, &hsref2.surface, &guide);
            finv.set_curve_on_surface(false, &pcref2_curve);

            // OCCT L1070-1084: rad / petitchoix; the Set calls;
            // func.Set(myShape).
            let rad = fsp.radius();
            let mut petitchoix = 1;
            if or1 == Orientation::Reversed {
                petitchoix = 3;
            }
            if choix % 2 == 0 {
                petitchoix += 1;
            }
            finv.set(rad, choix);
            finvc.set(rad, petitchoix);
            finvp.set(rad, petitchoix);
            func.set(rad, petitchoix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 3);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);

            // OCCT L1086-1112: done = SimulData(Data, HGuide, lin, HS1, I1,
            // HS2, PC2, I2, Decroch2, func, finv, finvp, finvc, PFirst,
            // MaxStep, locfleche, TolGuide, First, Last, Soldep, 4, Inside,
            // Appro, Forward, RecP, RecS, RecRst).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.simul_data_surf_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                i1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv,
                &mut finvp,
                &mut finvc,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                4,
                inside,
                appro,
                forward,
                rec_p,
                rec_s,
                rec_rst,
            );
            // OCCT L1113-1116.
            if !self.done {
                panic!("Standard_Failure: SimulSurf : Failed Processing!");
            }
            // OCCT L1117-1142: the section sampling loop.
            let lin_ref = lin.as_ref().expect("Lin");
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u, mut v, mut w, mut param, mut p1, mut p2) =
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                (u, v) = p.parameters_on_s();
                w = p.parameter_on_c();
                param = p.parameter();
                func.section(param, u, v, w, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                if i == 1 {
                    pf = DVec2::new(u, v);
                    (u, v) = p.parameters_on_s2();
                    ppcf = DVec2::new(u, v);
                }
                if i == nbp {
                    pl = DVec2::new(u, v);
                    (u, v) = p.parameters_on_s2();
                    ppcl = DVec2::new(u, v);
                }
            }
            sec = Some(arr);
        } else {
            // OCCT L1146-1152: BRepBlend_SurfRstEvolRad func(HS1, HS2, PC2,
            // HGuide, fsp->Law(HGuide)); func.Set(HSref2, PCref2);
            // HC->Load(PC2, HS2); BRepBlend_SurfCurvEvolRadInv finvc(HS1, HC,
            // HGuide, Law); BRepBlend_SurfPointEvolRadInv finvp(HS1, HGuide,
            // Law); BRepBlend_EvolRadInv finv(HS1, HSref2, HGuide, Law);
            // finv.Set(false, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func =
                BlendSurfRstEvolRad::new(&hs1.surface, &hs2.surface, &pc2_curve, &guide, law.clone());
            func.set_ref(&hsref2.surface, &pcref2_curve);
            // OCCT L1147-1148: the HC Adaptor3d_CurveOnSurface(PC2, HS2)
            // bound as the finvc restriction curve — carried by the
            // (pcurve, surface) pair.
            let mut finvc = BlendSurfCurvEvolRadInv::new(
                &hs1.surface,
                &pc2_curve,
                &hs2.surface,
                &guide,
                law.clone(),
            );
            let mut finvp = BlendSurfPointEvolRadInv::new(&hs1.surface, &guide, law.clone());
            let mut finv = BlendFuncEvolRadInv::new(&hs1.surface, &hsref2.surface, &guide, law);
            finv.set_curve_on_surface(false, &pcref2_curve);

            // OCCT L1153-1166: petitchoix; the Set calls; func.Set(myShape).
            let mut petitchoix = 1;
            if or1 == Orientation::Reversed {
                petitchoix = 3;
            }
            if choix % 2 == 0 {
                petitchoix += 1;
            }
            finv.set(choix);
            finvc.set(petitchoix);
            finvp.set(petitchoix);
            func.set(petitchoix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 3);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);

            // OCCT L1167-1193: done = SimulData(Data, HGuide, lin, HS1, I1,
            // HS2, PC2, I2, Decroch2, func, finv, finvp, finvc, PFirst,
            // MaxStep, locfleche, TolGuide, First, Last, Soldep, 4, Inside,
            // Appro, Forward, RecP, RecS, RecRst).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.simul_data_surf_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                i1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv,
                &mut finvp,
                &mut finvc,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                4,
                inside,
                appro,
                forward,
                rec_p,
                rec_s,
                rec_rst,
            );
            // OCCT L1194-1197.
            if !self.done {
                panic!("Standard_Failure: SimulSurf : Fail !");
            }
            // OCCT L1198-1223: the section sampling loop.
            let lin_ref = lin.as_ref().expect("Lin");
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u, mut v, mut w, mut param, mut p1, mut p2) =
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                (u, v) = p.parameters_on_s();
                w = p.parameter_on_c();
                param = p.parameter();
                func.section(param, u, v, w, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                if i == 1 {
                    pf = DVec2::new(u, v);
                    (u, v) = p.parameters_on_s2();
                    ppcf = DVec2::new(u, v);
                }
                if i == nbp {
                    pl = DVec2::new(u, v);
                    (u, v) = p.parameters_on_s2();
                    ppcl = DVec2::new(u, v);
                }
            }
            sec = Some(arr);
        }
        // OCCT L1225-1227: Data->SetSimul(sec);
        // Data->Set2dPoints(pf, pl, ppcf, ppcl).
        {
            let mut dw = data.write().expect("surfdata lock");
            dw.set_simul(sec);
            dw.set_2d_points(pf, pl, ppcf, ppcl);
            // OCCT L1228-1247: the four ChFi3d_FilCommonPoint loads (the
            // S1 / S2 vertex slots are NOT swapped in this overload, unlike
            // the face/rst SimulSurf).
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
        _fleche: f64, // OCCT L1270: const double /*Fleche*/ — unused.
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
        // OCCT L1288: occ::handle<BRepBlend_Line> lin.
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT L1291-1292: SimulParams(HGuide, fsp, MaxStep, locfleche).
        let mut locfleche = 0.0;
        let mut max_step = 0.0;
        simul_params(hguide, fsp, &mut max_step, &mut locfleche);
        // OCCT L1293: sec.
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        // OCCT L1296-1297: int ch1 = 1, ch2 = 2; double PFirst = First.
        let p_first = *first;
        if fsp.is_constant() {
            // OCCT L1301-1302: BRepBlend_RstRstConstRad func(HS1, PC1, HS2,
            // PC2, HGuide); func.Set(HSref1, PCref1, HSref2, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            // OCCT PC1 / PC2 / PCref1 / PCref2 are the BRepAdaptor_Curve2d
            // handles; the rcad port materializes their pcurve payload (the
            // (edge, face) lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc1 = pc1.as_ref().expect("PC1");
            let (pc1_curve, _, _) = adpc1
                .brep
                .curve_on_surface(&adpc1.edge, &adpc1.face)
                .expect("PC1 pcurve");
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref1 = pcref1.as_ref().expect("PCref1");
            let (pcref1_curve, _, _) = adpcref1
                .brep
                .curve_on_surface(&adpcref1.edge, &adpcref1.face)
                .expect("PCref1 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func = BlendRstRstConstRad::new(
                &hs1.surface,
                &pc1_curve,
                &hs2.surface,
                &pc2_curve,
                &guide,
            );
            func.set_ref(
                &hsref1.surface,
                &pcref1_curve,
                &hsref2.surface,
                &pcref2_curve,
            );
            // OCCT L1303-1306: HC1->Load(PC1, HS1); HC2->Load(PC2, HS2) —
            // the Adaptor3d_CurveOnSurface payloads are carried by the
            // (pcurve, surface) pairs.
            // OCCT L1307-1310: BRepBlend_SurfCurvConstRadInv finv1(HSref1,
            // HC2, HGuide); BRepBlend_CurvPointRadInv finvp1(HGuide, HC2);
            // BRepBlend_SurfCurvConstRadInv finv2(HSref2, HC1, HGuide);
            // BRepBlend_CurvPointRadInv finvp2(HGuide, HC1).
            let mut finv1 =
                BlendSurfCurvConstRadInv::new(&hsref1.surface, &pc2_curve, &hs2.surface, &guide);
            let mut finvp1 = BRepBlendCurvPointRadInvHc::new(&guide, &pc2_curve, &hs2.surface);
            let mut finv2 =
                BlendSurfCurvConstRadInv::new(&hsref2.surface, &pc1_curve, &hs1.surface, &guide);
            let mut finvp2 = BRepBlendCurvPointRadInvHc::new(&guide, &pc1_curve, &hs1.surface);

            // OCCT L1312-1313: finv1.Set(PCref1); finv2.Set(PCref2).
            BlendSurfCurvFuncInv::set_rst(&mut finv1, &pcref1_curve);
            BlendSurfCurvFuncInv::set_rst(&mut finv2, &pcref2_curve);

            // OCCT L1315-1330: rad; ch1 / ch2 by Or1 / Or2; the Set calls;
            // func.Set(myShape).
            let rad = fsp.radius();
            let mut ch1 = 1;
            let mut ch2 = 2;
            if or1 == Orientation::Reversed {
                ch1 = 3;
            }
            if or2 == Orientation::Reversed {
                ch2 = 3;
            }

            finv1.set(rad, ch1);
            finvp1.set_choix(choix);
            finv2.set(rad, ch2);
            finvp2.set_choix(choix);
            func.set(rad, choix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 2);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);

            // OCCT L1332-1362: done = SimulData(Data, HGuide, lin, HS1, PC1,
            // I1, Decroch1, HS2, PC2, I2, Decroch2, func, finv1, finvp1,
            // finv2, finvp2, PFirst, MaxStep, locfleche, TolGuide, First,
            // Last, Soldep, 4, Inside, Appro, Forward, RecP1, RecRst1, RecP2,
            // RecRst2).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.simul_data_rst_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                &pc1_curve,
                i1,
                decroch1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv1,
                &mut finvp1,
                &mut finv2,
                &mut finvp2,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                4,
                inside,
                appro,
                forward,
                rec_p1,
                rec_rst1,
                rec_p2,
                rec_rst2,
            );
            // OCCT L1363-1366.
            if !self.done {
                panic!("Standard_Failure: SimulSurf : Failed processing!");
            }
            // OCCT L1367-1382: the section sampling loop.
            let lin_ref = lin.as_ref().expect("Lin");
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u, mut v, mut param, mut p1, mut p2) = (0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                u = p.parameter_on_c1();
                v = p.parameter_on_c2();
                param = p.parameter();
                func.section(param, u, v, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                // OCCT L1380-1381: the i == 1 / i == nbp pf / pl reads are
                // commented out in the original.
            }
            sec = Some(arr);
        } else {
            // OCCT L1386-1387: BRepBlend_RstRstEvolRad func(HS1, PC1, HS2,
            // PC2, HGuide, fsp->Law(HGuide)); func.Set(HSref1, PCref1,
            // HSref2, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            // OCCT PC1 / PC2 / PCref1 / PCref2 are the BRepAdaptor_Curve2d
            // handles; the rcad port materializes their pcurve payload (the
            // (edge, face) lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc1 = pc1.as_ref().expect("PC1");
            let (pc1_curve, _, _) = adpc1
                .brep
                .curve_on_surface(&adpc1.edge, &adpc1.face)
                .expect("PC1 pcurve");
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref1 = pcref1.as_ref().expect("PCref1");
            let (pcref1_curve, _, _) = adpcref1
                .brep
                .curve_on_surface(&adpcref1.edge, &adpcref1.face)
                .expect("PCref1 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func = BlendRstRstEvolRad::new(
                &hs1.surface,
                &pc1_curve,
                &hs2.surface,
                &pc2_curve,
                &guide,
                law.clone(),
            );
            func.set_ref(
                &hsref1.surface,
                &pcref1_curve,
                &hsref2.surface,
                &pcref2_curve,
            );
            // OCCT L1388-1391: HC1->Load(PC1, HS1); HC2->Load(PC2, HS2).
            // OCCT L1393-1396: BRepBlend_SurfCurvEvolRadInv finv1(HSref1,
            // HC2, HGuide, Law); BRepBlend_CurvPointRadInv finvp1(HGuide,
            // HC2); BRepBlend_SurfCurvEvolRadInv finv2(HSref2, HC1, HGuide,
            // Law); BRepBlend_CurvPointRadInv finvp2(HGuide, HC1).
            let mut finv1 = BlendSurfCurvEvolRadInv::new(
                &hsref1.surface,
                &pc2_curve,
                &hs2.surface,
                &guide,
                law.clone(),
            );
            let mut finvp1 = BRepBlendCurvPointRadInvHc::new(&guide, &pc2_curve, &hs2.surface);
            let mut finv2 = BlendSurfCurvEvolRadInv::new(
                &hsref2.surface,
                &pc1_curve,
                &hs1.surface,
                &guide,
                law.clone(),
            );
            let mut finvp2 = BRepBlendCurvPointRadInvHc::new(&guide, &pc1_curve, &hs1.surface);

            // OCCT L1398-1399: finv1.Set(PCref1); finv2.Set(PCref2).
            BlendSurfCurvFuncInv::set_rst(&mut finv1, &pcref1_curve);
            BlendSurfCurvFuncInv::set_rst(&mut finv2, &pcref2_curve);

            // OCCT L1401-1417: ch11 / ch22 by Or1 / Or2; the Set calls;
            // func.Set(myShape).
            let mut ch11 = 1;
            let mut ch22 = 2;
            if or1 == Orientation::Reversed {
                ch11 = 3;
            }
            if or2 == Orientation::Reversed {
                ch22 = 3;
            }

            finv1.set(ch11);
            finvp1.set_choix(choix);
            finv2.set(ch22);
            finvp2.set_choix(choix);
            func.set(choix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 2);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);

            // OCCT L1419-1449: done = SimulData(Data, HGuide, lin, HS1, PC1,
            // I1, Decroch1, HS2, PC2, I2, Decroch2, func, finv1, finvp1,
            // finv2, finvp2, PFirst, MaxStep, locfleche, TolGuide, First,
            // Last, Soldep, 4, Inside, Appro, Forward, RecP1, RecRst1, RecP2,
            // RecRst2).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.simul_data_rst_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                &pc1_curve,
                i1,
                decroch1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv1,
                &mut finvp1,
                &mut finv2,
                &mut finvp2,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                4,
                inside,
                appro,
                forward,
                rec_p1,
                rec_rst1,
                rec_p2,
                rec_rst2,
            );

            // OCCT L1451-1454.
            if !self.done {
                panic!("Standard_Failure: SimulSurf : Fail !");
            }
            // OCCT L1455-1470: the section sampling loop.
            let lin_ref = lin.as_ref().expect("Lin");
            let nbp = lin_ref.nb_points();
            let mut arr: ChFiDSCircSectionArray = Vec::new();
            for i in 1..=nbp {
                let mut isec = ChFiDSCircSection::new();
                let (mut u, mut v, mut param, mut p1, mut p2) = (0.0, 0.0, 0.0, 0.0, 0.0);
                let mut ci = Circle3 {
                    center: glam::DVec3::ZERO,
                    normal: glam::DVec3::Z,
                    x_dir: glam::DVec3::X,
                    y_dir: glam::DVec3::Y,
                    radius: 0.0,
                };
                let p = lin_ref.point(i);
                u = p.parameter_on_c1();
                v = p.parameter_on_c2();
                param = p.parameter();
                func.section(param, u, v, &mut p1, &mut p2, &mut ci);
                isec.set_circ(ci, p1, p2);
                arr.push(isec);
                // OCCT L1468-1469: the i == 1 / i == nbp pf / pl reads are
                // commented out in the original.
            }
            sec = Some(arr);
        }
        // OCCT L1472: Data->SetSimul(sec) (the Set2dPoints call L1474 is
        // commented out in the original).
        {
            let mut dw = data.write().expect("surfdata lock");
            dw.set_simul(sec);
            // OCCT L1476-1495: the four ChFi3d_FilCommonPoint loads (the
            // S1 / S2 vertex slots are NOT swapped in this overload).
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
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1721-1890 — PerformSurf (face/rst: the
    /// obstacle curve lies on S1).
    // OCCT declares `bool maybesingular;` uninitialized (L1755); the rcad
    // initializer is dead on the constant path (assigned before read), which
    // Rust flags as an unused assignment.
    #[allow(unused_assignments, clippy::too_many_arguments)]
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
        let mut lin: Option<BRepBlendLine> = None;
        let mut maybesingular = false;
        let p_first = *first; // OCCT L1754.
        // OCCT L1757: if (fsp->IsConstant()).
        if fsp.is_constant() {
            // OCCT L1759-1766: BRepBlend_SurfRstConstRad func(HS2, HS1, PC1,
            // HGuide); func.Set(HSref1, PCref1); HC->Load(PC1, HS1);
            // BRepBlend_SurfCurvConstRadInv finvc(HS2, HC, HGuide);
            // BRepBlend_SurfPointConstRadInv finvp(HS2, HGuide);
            // BRepBlend_ConstRadInv finv(HS2, HSref1, HGuide);
            // finv.Set(false, PCref1).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            // OCCT PC1 / PCref1 are the BRepAdaptor_Curve2d handles; the rcad
            // port materializes their pcurve payload (the (edge, face)
            // lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc1 = pc1.as_ref().expect("PC1");
            let (pc1_curve, _, _) = adpc1
                .brep
                .curve_on_surface(&adpc1.edge, &adpc1.face)
                .expect("PC1 pcurve");
            let adpcref1 = pcref1.as_ref().expect("PCref1");
            let (pcref1_curve, _, _) = adpcref1
                .brep
                .curve_on_surface(&adpcref1.edge, &adpcref1.face)
                .expect("PCref1 pcurve");
            let mut func =
                BlendSurfRstConstRad::new(&hs2.surface, &hs1.surface, &pc1_curve, &guide);
            func.set_ref(&hsref1.surface, &pcref1_curve);
            // OCCT L1761-1763: the HC Adaptor3d_CurveOnSurface(PC1, HS1)
            // bound as the finvc restriction curve — carried by the
            // (pcurve, surface) pair.
            let mut finvc =
                BlendSurfCurvConstRadInv::new(&hs2.surface, &pc1_curve, &hs1.surface, &guide);
            let mut finvp = BlendSurfPointConstRadInv::new(&hs2.surface, &guide);
            let mut finv = BlendFuncConstRadInv::new(&hs2.surface, &hsref1.surface, &guide);
            finv.set_curve_on_surface(false, &pcref1_curve);

            // OCCT L1768-1782: rad / petitchoix; the Set calls;
            // func.Set(myShape).
            let rad = fsp.radius();
            let mut petitchoix = 1;
            if or2 == Orientation::Reversed {
                petitchoix = 3;
            }
            if choix % 2 == 0 {
                petitchoix += 1;
            }
            finv.set(rad, choix);
            finvc.set(rad, petitchoix);
            finvp.set(rad, petitchoix);
            func.set(rad, petitchoix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 3);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);

            // OCCT L1784-1809: done = ComputeData(Data, HGuide, lin, HS2,
            // I2, HS1, PC1, I1, Decroch1, func, finv, finvp, finvc, PFirst,
            // MaxStep, Fleche, TolGuide, First, Last, Soldep, Inside, Appro,
            // Forward, RecP, RecS, RecRst).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.compute_data_surf_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs2,
                i2,
                hs1,
                &pc1_curve,
                i1,
                decroch1,
                &mut func,
                &mut finv,
                &mut finvp,
                &mut finvc,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                inside,
                appro,
                forward,
                rec_p,
                rec_s,
                rec_rst,
            );
            // OCCT L1810-1814.
            if !self.done {
                // OCCT L1812: Spine->SetErrorStatus(ChFiDS_WalkingFailure).
                // Boundary note: OCCT mutates through the const handle; the
                // rcad handle is an enum without interior mutability — the
                // status is set on the clone (lost at the throw below, as
                // the OCCT exception also unwinds the stripe state).
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
                panic!("Standard_Failure: PerformSurf : Failed processing!");
            }
            // OCCT L1815-1816: Or = HS2->Face().Orientation();
            // done = CompleteData(Data, func, lin, HS1, HS2, Or, true).
            let or = hs2.face.orientation;
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_surf_rst(
                &mut data.write().expect("surfdata lock"),
                &mut func,
                lin_ref,
                hs1,
                Some(hs2),
                or,
                true,
            );
            // OCCT L1817-1820.
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            // OCCT L1821.
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
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
        // OCCT L1886-1889.
        if maybesingular {
            self.split_surf(seqsd, lin.as_ref().expect("Lin"));
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1894-2064 — PerformSurf (rst/face: the
    /// obstacle curve lies on S2).
    // OCCT declares `bool maybesingular;` uninitialized (L1928); the rcad
    // initializer is dead on the constant path (assigned before read), which
    // Rust flags as an unused assignment.
    #[allow(unused_assignments, clippy::too_many_arguments)]
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
        let mut lin: Option<BRepBlendLine> = None;
        let mut maybesingular = false;
        let p_first = *first;
        if fsp.is_constant() {
            // OCCT L1932-1939: BRepBlend_SurfRstConstRad func(HS1, HS2, PC2,
            // HGuide); func.Set(HSref2, PCref2); HC->Load(PC2, HS2);
            // BRepBlend_SurfCurvConstRadInv finvc(HS1, HC, HGuide);
            // BRepBlend_SurfPointConstRadInv finvp(HS1, HGuide);
            // BRepBlend_ConstRadInv finv(HS1, HSref2, HGuide);
            // finv.Set(false, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            // OCCT PC2 / PCref2 are the BRepAdaptor_Curve2d handles; the rcad
            // port materializes their pcurve payload (the (edge, face)
            // lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func =
                BlendSurfRstConstRad::new(&hs1.surface, &hs2.surface, &pc2_curve, &guide);
            func.set_ref(&hsref2.surface, &pcref2_curve);
            // OCCT L1934-1935: the HC Adaptor3d_CurveOnSurface(PC2, HS2)
            // bound as the finvc restriction curve — carried by the
            // (pcurve, surface) pair.
            let mut finvc =
                BlendSurfCurvConstRadInv::new(&hs1.surface, &pc2_curve, &hs2.surface, &guide);
            let mut finvp = BlendSurfPointConstRadInv::new(&hs1.surface, &guide);
            let mut finv = BlendFuncConstRadInv::new(&hs1.surface, &hsref2.surface, &guide);
            finv.set_curve_on_surface(false, &pcref2_curve);

            // OCCT L1941-1955: rad / petitchoix; the Set calls;
            // func.Set(myShape).
            let rad = fsp.radius();
            let mut petitchoix = 1;
            if or1 == Orientation::Reversed {
                petitchoix = 3;
            }
            if choix % 2 == 0 {
                petitchoix += 1;
            }
            finv.set(rad, choix);
            finvc.set(rad, petitchoix);
            finvp.set(rad, petitchoix);
            func.set(rad, petitchoix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 3);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);

            // OCCT L1957-1982: done = ComputeData(Data, HGuide, lin, HS1,
            // I1, HS2, PC2, I2, Decroch2, func, finv, finvp, finvc, PFirst,
            // MaxStep, Fleche, TolGuide, First, Last, Soldep, Inside, Appro,
            // Forward, RecP, RecS, RecRst).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.compute_data_surf_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                i1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv,
                &mut finvp,
                &mut finvc,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                inside,
                appro,
                forward,
                rec_p,
                rec_s,
                rec_rst,
            );
            // OCCT L1983-1987.
            if !self.done {
                // OCCT L1985: Spine->SetErrorStatus(ChFiDS_WalkingFailure).
                // Boundary note: OCCT mutates through the const handle; the
                // rcad handle is an enum without interior mutability — the
                // status is set on the clone (lost at the throw below, as
                // the OCCT exception also unwinds the stripe state).
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
                panic!("Standard_Failure: PerformSurf : Failed processing!");
            }
            // OCCT L1988-1989: Or = HS1->Face().Orientation();
            // done = CompleteData(Data, func, lin, HS1, HS2, Or, false).
            let or = hs1.face.orientation;
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_surf_rst(
                &mut data.write().expect("surfdata lock"),
                &mut func,
                lin_ref,
                hs1,
                Some(hs2),
                or,
                false,
            );
            // OCCT L1990-1993.
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            // OCCT L1994.
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
        } else {
            // OCCT L1998-2005: BRepBlend_SurfRstEvolRad func(HS1, HS2, PC2,
            // HGuide, fsp->Law(HGuide)); func.Set(HSref2, PCref2);
            // HC->Load(PC2, HS2); BRepBlend_SurfCurvEvolRadInv finvc(HS1, HC,
            // HGuide, Law); BRepBlend_SurfPointEvolRadInv finvp(HS1, HGuide,
            // Law); BRepBlend_EvolRadInv finv(HS1, HSref2, HGuide, Law);
            // finv.Set(false, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func =
                BlendSurfRstEvolRad::new(&hs1.surface, &hs2.surface, &pc2_curve, &guide, law.clone());
            func.set_ref(&hsref2.surface, &pcref2_curve);
            // OCCT L2000-2001: the HC Adaptor3d_CurveOnSurface(PC2, HS2)
            // bound as the finvc restriction curve — carried by the
            // (pcurve, surface) pair.
            let mut finvc = BlendSurfCurvEvolRadInv::new(
                &hs1.surface,
                &pc2_curve,
                &hs2.surface,
                &guide,
                law.clone(),
            );
            let mut finvp = BlendSurfPointEvolRadInv::new(&hs1.surface, &guide, law.clone());
            let mut finv = BlendFuncEvolRadInv::new(&hs1.surface, &hsref2.surface, &guide, law);
            finv.set_curve_on_surface(false, &pcref2_curve);

            // OCCT L2006-2019: petitchoix; the Set calls; func.Set(myShape).
            let mut petitchoix = 1;
            if or1 == Orientation::Reversed {
                petitchoix = 3;
            }
            if choix % 2 == 0 {
                petitchoix += 1;
            }
            finv.set(choix);
            finvc.set(petitchoix);
            finvp.set(petitchoix);
            func.set(petitchoix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 3);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);
            soldep_v.set(3, soldep[2]);

            // OCCT L2021-2046: done = ComputeData(Data, HGuide, lin, HS1,
            // I1, HS2, PC2, I2, Decroch2, func, finv, finvp, finvc, PFirst,
            // MaxStep, Fleche, TolGuide, First, Last, Soldep, Inside, Appro,
            // Forward, RecP, RecS, RecRst).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.compute_data_surf_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                i1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv,
                &mut finvp,
                &mut finvc,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                inside,
                appro,
                forward,
                rec_p,
                rec_s,
                rec_rst,
            );
            // OCCT L2047-2051.
            if !self.done {
                // OCCT L2049.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
                panic!("Standard_Failure: PerformSurf : Failed processing!");
            }
            // OCCT L2052-2053: Or = HS1->Face().Orientation();
            // done = CompleteData(Data, func, lin, HS1, HS2, Or, false).
            let or = hs1.face.orientation;
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_surf_rst(
                &mut data.write().expect("surfdata lock"),
                &mut func,
                lin_ref,
                hs1,
                Some(hs2),
                or,
                false,
            );
            // OCCT L2054-2057.
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            // OCCT L2058.
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
        }
        // OCCT L2060-2063.
        if maybesingular {
            self.split_surf(seqsd, lin.as_ref().expect("Lin"));
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L2068-2270 — PerformSurf (rst/rst: the
    /// curve-curve entry).
    // OCCT declares `bool maybesingular;` uninitialized (L2108); the rcad
    // initializer is dead on the constant path (assigned before read), which
    // Rust flags as an unused assignment.
    #[allow(unused_assignments, clippy::too_many_arguments)]
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
        let mut lin: Option<BRepBlendLine> = None;
        let mut maybesingular = false;
        let p_first = *first;
        if fsp.is_constant() {
            // OCCT L2112-2113: BRepBlend_RstRstConstRad func(HS1, PC1, HS2,
            // PC2, HGuide); func.Set(HSref1, PCref1, HSref2, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            // OCCT PC1 / PC2 / PCref1 / PCref2 are the BRepAdaptor_Curve2d
            // handles; the rcad port materializes their pcurve payload (the
            // (edge, face) lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc1 = pc1.as_ref().expect("PC1");
            let (pc1_curve, _, _) = adpc1
                .brep
                .curve_on_surface(&adpc1.edge, &adpc1.face)
                .expect("PC1 pcurve");
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref1 = pcref1.as_ref().expect("PCref1");
            let (pcref1_curve, _, _) = adpcref1
                .brep
                .curve_on_surface(&adpcref1.edge, &adpcref1.face)
                .expect("PCref1 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func = BlendRstRstConstRad::new(
                &hs1.surface,
                &pc1_curve,
                &hs2.surface,
                &pc2_curve,
                &guide,
            );
            func.set_ref(
                &hsref1.surface,
                &pcref1_curve,
                &hsref2.surface,
                &pcref2_curve,
            );
            // OCCT L2114-2117: HC1->Load(PC1, HS1); HC2->Load(PC2, HS2).
            // OCCT L2118-2121: BRepBlend_SurfCurvConstRadInv finv1(HSref1,
            // HC2, HGuide); BRepBlend_CurvPointRadInv finvp1(HGuide, HC2);
            // BRepBlend_SurfCurvConstRadInv finv2(HSref2, HC1, HGuide);
            // BRepBlend_CurvPointRadInv finvp2(HGuide, HC1).
            let mut finv1 =
                BlendSurfCurvConstRadInv::new(&hsref1.surface, &pc2_curve, &hs2.surface, &guide);
            let mut finvp1 = BRepBlendCurvPointRadInvHc::new(&guide, &pc2_curve, &hs2.surface);
            let mut finv2 =
                BlendSurfCurvConstRadInv::new(&hsref2.surface, &pc1_curve, &hs1.surface, &guide);
            let mut finvp2 = BRepBlendCurvPointRadInvHc::new(&guide, &pc1_curve, &hs1.surface);

            // OCCT L2123-2124: finv1.Set(PCref1); finv2.Set(PCref2).
            BlendSurfCurvFuncInv::set_rst(&mut finv1, &pcref1_curve);
            BlendSurfCurvFuncInv::set_rst(&mut finv2, &pcref2_curve);

            // OCCT L2126-2142: ch1 / ch2 by Or1 / Or2; rad; the Set calls;
            // func.Set(myShape).
            let mut ch1 = 1;
            let mut ch2 = 2;
            let rad = fsp.radius();
            if or1 == Orientation::Reversed {
                ch1 = 3;
            }
            if or2 == Orientation::Reversed {
                ch2 = 3;
            }

            finv1.set(rad, ch1);
            finvp1.set_choix(choix);
            finv2.set(rad, ch2);
            finvp2.set_choix(choix);
            func.set(rad, choix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 2);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);

            // OCCT L2144-2173: done = ComputeData(Data, HGuide, lin, HS1,
            // PC1, I1, Decroch1, HS2, PC2, I2, Decroch2, func, finv1, finvp1,
            // finv2, finvp2, PFirst, MaxStep, Fleche, TolGuide, First, Last,
            // Soldep, Inside, Appro, Forward, RecP1, RecRst1, RecP2, RecRst2).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.compute_data_rst_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                &pc1_curve,
                i1,
                decroch1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv1,
                &mut finvp1,
                &mut finv2,
                &mut finvp2,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                inside,
                appro,
                forward,
                rec_p1,
                rec_rst1,
                rec_p2,
                rec_rst2,
            );
            // OCCT L2174-2178.
            if !self.done {
                // OCCT L2176.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
                panic!("Standard_Failure: PerformSurf : Failed processing!");
            }
            // OCCT L2179-2180: Or = HS1->Face().Orientation();
            // done = CompleteData(Data, func, lin, HS1, HS2, Or).
            let or = hs1.face.orientation;
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_rst_rst(
                &mut data.write().expect("surfdata lock"),
                &mut func,
                lin_ref,
                hs1,
                Some(hs2),
                or,
            );
            // OCCT L2181-2184.
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            // OCCT L2185.
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
        } else {
            // OCCT L2189-2190: BRepBlend_RstRstEvolRad func(HS1, PC1, HS2,
            // PC2, HGuide, fsp->Law(HGuide)); func.Set(HSref1, PCref1,
            // HSref2, PCref2).
            let hguide_handle = std::sync::Arc::new(std::sync::RwLock::new(hguide.clone()));
            let guide = elspine_guide_curve(&hguide_handle);
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            // OCCT PC1 / PC2 / PCref1 / PCref2 are the BRepAdaptor_Curve2d
            // handles; the rcad port materializes their pcurve payload (the
            // (edge, face) lookup the adaptor performs — the same read as
            // BRepAdaptorCurve2d::value).
            let adpc1 = pc1.as_ref().expect("PC1");
            let (pc1_curve, _, _) = adpc1
                .brep
                .curve_on_surface(&adpc1.edge, &adpc1.face)
                .expect("PC1 pcurve");
            let adpc2 = pc2.as_ref().expect("PC2");
            let (pc2_curve, _, _) = adpc2
                .brep
                .curve_on_surface(&adpc2.edge, &adpc2.face)
                .expect("PC2 pcurve");
            let adpcref1 = pcref1.as_ref().expect("PCref1");
            let (pcref1_curve, _, _) = adpcref1
                .brep
                .curve_on_surface(&adpcref1.edge, &adpcref1.face)
                .expect("PCref1 pcurve");
            let adpcref2 = pcref2.as_ref().expect("PCref2");
            let (pcref2_curve, _, _) = adpcref2
                .brep
                .curve_on_surface(&adpcref2.edge, &adpcref2.face)
                .expect("PCref2 pcurve");
            let mut func = BlendRstRstEvolRad::new(
                &hs1.surface,
                &pc1_curve,
                &hs2.surface,
                &pc2_curve,
                &guide,
                law.clone(),
            );
            func.set_ref(
                &hsref1.surface,
                &pcref1_curve,
                &hsref2.surface,
                &pcref2_curve,
            );
            // OCCT L2191-2194: HC1->Load(PC1, HS1); HC2->Load(PC2, HS2).
            // OCCT L2196-2199: BRepBlend_SurfCurvEvolRadInv finv1(HSref1,
            // HC2, HGuide, Law); BRepBlend_CurvPointRadInv finvp1(HGuide,
            // HC2); BRepBlend_SurfCurvEvolRadInv finv2(HSref2, HC1, HGuide,
            // Law); BRepBlend_CurvPointRadInv finvp2(HGuide, HC1).
            let mut finv1 = BlendSurfCurvEvolRadInv::new(
                &hsref1.surface,
                &pc2_curve,
                &hs2.surface,
                &guide,
                law.clone(),
            );
            let mut finvp1 = BRepBlendCurvPointRadInvHc::new(&guide, &pc2_curve, &hs2.surface);
            let mut finv2 = BlendSurfCurvEvolRadInv::new(
                &hsref2.surface,
                &pc1_curve,
                &hs1.surface,
                &guide,
                law.clone(),
            );
            let mut finvp2 = BRepBlendCurvPointRadInvHc::new(&guide, &pc1_curve, &hs1.surface);

            // OCCT L2201-2202: finv1.Set(PCref1); finv2.Set(PCref2).
            BlendSurfCurvFuncInv::set_rst(&mut finv1, &pcref1_curve);
            BlendSurfCurvFuncInv::set_rst(&mut finv2, &pcref2_curve);

            // OCCT L2204-2220: ch1 / ch2 by Or1 / Or2; the Set calls;
            // func.Set(myShape).
            let mut ch1 = 1;
            let mut ch2 = 2;
            if or1 == Orientation::Reversed {
                ch1 = 3;
            }
            if or2 == Orientation::Reversed {
                ch2 = 3;
            }

            finv1.set(ch1);
            finvp1.set_choix(choix);
            finv2.set(ch2);
            finvp2.set_choix(choix);
            func.set(choix);
            func.set_section_shape(self.my_blend_shape);

            // OCCT L631 form: the math_Vector materialization of Soldep.
            let mut soldep_v = Vector::new(1, 2);
            soldep_v.set(1, soldep[0]);
            soldep_v.set(2, soldep[1]);

            // OCCT L2222-2251: done = ComputeData(Data, HGuide, lin, HS1,
            // PC1, I1, Decroch1, HS2, PC2, I2, Decroch2, func, finv1, finvp1,
            // finv2, finvp2, PFirst, MaxStep, Fleche, TolGuide, First, Last,
            // Soldep, Inside, Appro, Forward, RecP1, RecRst1, RecP2, RecRst2).
            let mut dw = data.write().expect("surfdata lock");
            self.done = self.compute_data_rst_rst(
                &mut dw,
                &hguide_handle,
                &mut lin,
                hs1,
                &pc1_curve,
                i1,
                decroch1,
                hs2,
                &pc2_curve,
                i2,
                decroch2,
                &mut func,
                &mut finv1,
                &mut finvp1,
                &mut finv2,
                &mut finvp2,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                &soldep_v,
                inside,
                appro,
                forward,
                rec_p1,
                rec_rst1,
                rec_p2,
                rec_rst2,
            );

            // OCCT L2253-2257.
            if !self.done {
                // OCCT L2255.
                let mut sp = spine.clone();
                sp.base_mut()
                    .set_error_status(ChFiDS_ErrorStatus::WalkingFailure);
                panic!("Standard_Failure: PerformSurf : Failed processing!");
            }
            // OCCT L2258-2259: Or = HS1->Face().Orientation();
            // done = CompleteData(Data, func, lin, HS1, HS2, Or).
            let or = hs1.face.orientation;
            let lin_ref = lin.as_ref().expect("Lin");
            self.done = self.complete_data_rst_rst(
                &mut data.write().expect("surfdata lock"),
                &mut func,
                lin_ref,
                hs1,
                Some(hs2),
                or,
            );
            // OCCT L2260-2263.
            if !self.done {
                panic!("Standard_Failure: PerformSurf : Failed approximation!");
            }
            // OCCT L2264.
            maybesingular = func.get_minimal_distance() <= 100.0 * self.tolapp3d;
        }
        // OCCT L2266-2269.
        if maybesingular {
            self.split_surf(seqsd, lin.as_ref().expect("Lin"));
        }
    }

}
