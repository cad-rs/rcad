//! OCCT ChFi3d_ChBuilder.cxx (continued) — the SimulSurf and PerformSurf
//! rst overloads that OCCT implements as Standard_Failure throws
//! (ChFi3d_ChBuilder.cxx L1223-1311, L1854-1945).
//!
//! Split from `chfi3d_builder_chbuilder` per the 2000-line guideline.
//! These are virtual-slot overloads of the ChFi3d_FilBuilder surface:
//! the chamfer builder only realizes the face/face entry; the rst paths
//! deliberately throw, and the translation keeps the exact parameter
//! lists as the OCCT virtual signatures.

use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::Orientation;

use super::chfi3d_builder_0::BRepAdaptorSurface;
use super::chfi3d_builder_6b::ChFiDSElSpineHandle;
use super::chfi3d_builder_2::BRepTopAdaptorTopolTool;
use super::chfi_ds::{ChFiDSSpineHandle, ChFiDSSurfData};
use super::chfi3d::ChFi3dChBuilder;

impl ChFi3dChBuilder {
    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1223-1249 — SimulSurf (the face/rst
    // overload) — throws Standard_Failure.
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn simul_surf_face_rst(
        &mut self,
        _data: &mut ChFiDSSurfData,
        _hguide: &ChFiDSElSpineHandle,
        _spine: &ChFiDSSpineHandle,
        _choix: i32,
        _s1: &BRepAdaptorSurface,
        _i1: &BRepTopAdaptorTopolTool,
        _pc2: &rcad_kernel::geom::Curve2d,
        _s2: &BRepAdaptorSurface,
        _pc2_bis: &rcad_kernel::geom::Curve2d,
        _decroch: &mut bool,
        _i2: &BRepTopAdaptorTopolTool,
        _or: Orientation,
        _tol_guide: f64,
        _fleche: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _appro: bool,
        _forward: bool,
        _rec_on_s1: bool,
        _rec_on_s2: bool,
        _rec_rst: bool,
        _soldep: &Vector,
    ) {
        panic!("Standard_Failure: SimulSurf Not Implemented");
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1251-1277 — SimulSurf (the rst/face
    // overload) — throws Standard_Failure.
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn simul_surf_rst_face(
        &mut self,
        _data: &mut ChFiDSSurfData,
        _hguide: &ChFiDSElSpineHandle,
        _spine: &ChFiDSSpineHandle,
        _choix: i32,
        _s1: &BRepAdaptorSurface,
        _i1: &BRepTopAdaptorTopolTool,
        _or: Orientation,
        _s2: &BRepAdaptorSurface,
        _i2: &BRepTopAdaptorTopolTool,
        _pc2: &rcad_kernel::geom::Curve2d,
        _s2_bis: &BRepAdaptorSurface,
        _pc2_bis: &rcad_kernel::geom::Curve2d,
        _decroch: &mut bool,
        _tol_guide: f64,
        _fleche: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _appro: bool,
        _forward: bool,
        _rec_on_s1: bool,
        _rec_on_s2: bool,
        _rec_rst: bool,
        _soldep: &Vector,
    ) {
        panic!("Standard_Failure: SimulSurf Not Implemented");
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1279-1311 — SimulSurf (the rst/rst
    // overload) — throws Standard_Failure.
    // =====================================================================
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn simul_surf_rst_rst(
        &mut self,
        _data: &mut ChFiDSSurfData,
        _hguide: &ChFiDSElSpineHandle,
        _spine: &ChFiDSSpineHandle,
        _choix: i32,
        _s1: &BRepAdaptorSurface,
        _i1: &BRepTopAdaptorTopolTool,
        _pc1: &rcad_kernel::geom::Curve2d,
        _s2: &BRepAdaptorSurface,
        _pc2: &rcad_kernel::geom::Curve2d,
        _decroch1: &mut bool,
        _i1_bis: &BRepTopAdaptorTopolTool,
        _s1_bis: &BRepAdaptorSurface,
        _i2: &BRepTopAdaptorTopolTool,
        _s2_bis: &BRepAdaptorSurface,
        _pc2_bis: &rcad_kernel::geom::Curve2d,
        _decroch2: &mut bool,
        _or1: Orientation,
        _or2: Orientation,
        _tol_guide: f64,
        _fleche: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _appro: bool,
        _forward: bool,
        _rec_on_s1: bool,
        _rec_on_s2: bool,
        _rec_rst1: bool,
        _rec_rst2: bool,
        _soldep: &Vector,
    ) {
        panic!("Standard_Failure: SimulSurf Not Implemented");
    }

    // OCCT ChFi3d_ChBuilder.cxx L1854-1881 / L1883-1910 / L1912-1945 — the
    // PerformSurf rst overloads — throw Standard_Failure.
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf_face_rst(
        &mut self,
        _seq_data: &mut Vec<std::sync::Arc<std::sync::RwLock<ChFiDSSurfData>>>,
        _hguide: &ChFiDSElSpineHandle,
        _spine: &ChFiDSSpineHandle,
        _choix: i32,
        _s1: &BRepAdaptorSurface,
        _i1: &BRepTopAdaptorTopolTool,
        _pc2: &rcad_kernel::geom::Curve2d,
        _s2: &BRepAdaptorSurface,
        _pc2_bis: &rcad_kernel::geom::Curve2d,
        _decroch: &mut bool,
        _i2: &BRepTopAdaptorTopolTool,
        _or: Orientation,
        _max_step: f64,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _appro: bool,
        _forward: bool,
        _rec_on_s1: bool,
        _rec_on_s2: bool,
        _rec_rst: bool,
        _soldep: &Vector,
    ) {
        panic!("Standard_Failure: PerformSurf Not Implemented");
    }

    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf_rst_face(
        &mut self,
        _seq_data: &mut Vec<std::sync::Arc<std::sync::RwLock<ChFiDSSurfData>>>,
        _hguide: &ChFiDSElSpineHandle,
        _spine: &ChFiDSSpineHandle,
        _choix: i32,
        _s1: &BRepAdaptorSurface,
        _i1: &BRepTopAdaptorTopolTool,
        _or: Orientation,
        _s2: &BRepAdaptorSurface,
        _i2: &BRepTopAdaptorTopolTool,
        _pc2: &rcad_kernel::geom::Curve2d,
        _s2_bis: &BRepAdaptorSurface,
        _pc2_bis: &rcad_kernel::geom::Curve2d,
        _decroch: &mut bool,
        _max_step: f64,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _appro: bool,
        _forward: bool,
        _rec_on_s1: bool,
        _rec_on_s2: bool,
        _rec_rst: bool,
        _soldep: &Vector,
    ) {
        panic!("Standard_Failure: PerformSurf Not Implemented");
    }

    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf_rst_rst(
        &mut self,
        _seq_data: &mut Vec<std::sync::Arc<std::sync::RwLock<ChFiDSSurfData>>>,
        _hguide: &ChFiDSElSpineHandle,
        _spine: &ChFiDSSpineHandle,
        _choix: i32,
        _s1: &BRepAdaptorSurface,
        _i1: &BRepTopAdaptorTopolTool,
        _pc1: &rcad_kernel::geom::Curve2d,
        _s2: &BRepAdaptorSurface,
        _pc2: &rcad_kernel::geom::Curve2d,
        _decroch1: &mut bool,
        _or1: Orientation,
        _i1_bis: &BRepTopAdaptorTopolTool,
        _s1_bis: &BRepAdaptorSurface,
        _i2: &BRepTopAdaptorTopolTool,
        _s2_bis: &BRepAdaptorSurface,
        _pc2_bis: &rcad_kernel::geom::Curve2d,
        _decroch2: &mut bool,
        _or2: Orientation,
        _max_step: f64,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _appro: bool,
        _forward: bool,
        _rec_on_s1: bool,
        _rec_on_s2: bool,
        _rec_rst1: bool,
        _rec_rst2: bool,
        _soldep: &Vector,
    ) {
        panic!("Standard_Failure: PerformSurf Not Implemented");
    }
}
