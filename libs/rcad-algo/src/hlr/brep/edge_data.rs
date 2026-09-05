// OCCT HLRBRep_EdgeData (TKHLR) — the HLR per-edge data record.
//
// HLRBRep_EdgeData.hxx L17-173 + .cxx L24-59 (Set) + .lxx (the flag
// accessors over the EMaskFlags bit system).  The OCCT `myGeometry`
// (HLRBRep_Curve) keeps the loaded edge adaptor through a handle; rcad
// holds the [`Curve`] with the borrowed view lifetime of the HLRBRep_Data
// object graph.
//
// Documented deferral: `Set` reads the TopoDS edge tolerance through
// BRep_Tool::Tolerance and loads the BRepAdaptor_Curve; the rcad caller
// (the Stage 3f Data) passes the loaded [`Curve`] and the tolerance — the
// BRep_Tool read happens at the kernel boundary.

use crate::hlr::algo::edge_status::EdgeStatus;
use crate::hlr::algo::edges_block::MinMaxIndices;

use super::curve::Curve;

/// OCCT HLRBRep_EdgeData.
pub struct EdgeData<'a> {
    my_flags: i32,
    my_hide_count: i32,
    my_v_sta: i32,
    my_v_end: i32,
    my_min_max: MinMaxIndices,
    my_status: EdgeStatus,
    my_geometry: Curve<'a>,
    my_tolerance: f32,
}

/// OCCT `enum EMaskFlags` (hxx L141-158).
mod emask_flags {
    pub const EMASK_SELECTED: i32 = 1;
    pub const EMASK_USED: i32 = 2;
    pub const EMASK_RG1_LINE: i32 = 4;
    pub const EMASK_VERTICAL: i32 = 8;
    pub const EMASK_SIMPLE: i32 = 16;
    pub const EMASK_OUT_LV_STA: i32 = 32;
    pub const EMASK_OUT_LV_END: i32 = 64;
    pub const EMASK_INT_DONE: i32 = 128;
    pub const EMASK_CUT_AT_STA: i32 = 256;
    pub const EMASK_CUT_AT_END: i32 = 512;
    pub const EMASK_VER_AT_STA: i32 = 1024;
    pub const EMASK_VER_AT_END: i32 = 2048;
    pub const EMASK_RG_N_LINE: i32 = 4096;
}

use emask_flags::*;

impl<'a> EdgeData<'a> {
    /// OCCT HLRBRep_EdgeData() (hxx L30-35) — Selected(true).
    pub fn new() -> Self {
        let mut r = EdgeData {
            my_flags: 0,
            my_hide_count: 0,
            my_v_sta: 0,
            my_v_end: 0,
            my_min_max: MinMaxIndices::default(),
            my_status: EdgeStatus::new(),
            my_geometry: Curve::new(),
            my_tolerance: 0.0,
        };
        r.set_selected(true);
        r
    }

    /// OCCT Set (cxx L33-59) — the BRep_Tool tolerance read and the edge
    /// adaptor load are the caller's (the Data) kernel boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn set(
        &mut self,
        rg1_l: bool,
        rg_nl: bool,
        geometry: Curve<'a>,
        tolerance: f32,
        v1: i32,
        v2: i32,
        out1: bool,
        out2: bool,
        cut1: bool,
        cut2: bool,
        start: f64,
        tol_start: f32,
        end: f64,
        tol_end: f32,
    ) {
        self.set_rg1_line(rg1_l);
        self.set_rg_n_line(rg_nl);
        self.set_used(false);
        self.my_geometry = geometry;
        self.my_tolerance = tolerance;
        self.set_v_sta(v1);
        self.set_v_end(v2);
        self.set_out_lv_sta(out1);
        self.set_out_lv_end(out2);
        self.set_cut_at_sta(cut1);
        self.set_cut_at_end(cut2);
        self.my_status.initialize(
            start,
            self.my_geometry.resolution(tol_start as f64) as f32,
            end,
            self.my_geometry.resolution(tol_end as f64) as f32,
        );
    }

    /// OCCT Selected (lxx).
    pub fn selected(&self) -> bool {
        (self.my_flags & EMASK_SELECTED) != 0
    }
    pub fn set_selected(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_SELECTED;
        } else {
            self.my_flags &= !EMASK_SELECTED;
        }
    }

    /// OCCT Rg1Line (lxx).
    pub fn rg1_line(&self) -> bool {
        (self.my_flags & EMASK_RG1_LINE) != 0
    }
    pub fn set_rg1_line(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_RG1_LINE;
        } else {
            self.my_flags &= !EMASK_RG1_LINE;
        }
    }

    /// OCCT RgNLine (lxx).
    pub fn rg_n_line(&self) -> bool {
        (self.my_flags & EMASK_RG_N_LINE) != 0
    }
    pub fn set_rg_n_line(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_RG_N_LINE;
        } else {
            self.my_flags &= !EMASK_RG_N_LINE;
        }
    }

    /// OCCT Vertical (lxx).
    pub fn vertical(&self) -> bool {
        (self.my_flags & EMASK_VERTICAL) != 0
    }
    pub fn set_vertical(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_VERTICAL;
        } else {
            self.my_flags &= !EMASK_VERTICAL;
        }
    }

    /// OCCT Simple (lxx).
    pub fn simple(&self) -> bool {
        (self.my_flags & EMASK_SIMPLE) != 0
    }
    pub fn set_simple(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_SIMPLE;
        } else {
            self.my_flags &= !EMASK_SIMPLE;
        }
    }

    /// OCCT OutLVSta (lxx).
    pub fn out_lv_sta(&self) -> bool {
        (self.my_flags & EMASK_OUT_LV_STA) != 0
    }
    pub fn set_out_lv_sta(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_OUT_LV_STA;
        } else {
            self.my_flags &= !EMASK_OUT_LV_STA;
        }
    }

    /// OCCT OutLVEnd (lxx).
    pub fn out_lv_end(&self) -> bool {
        (self.my_flags & EMASK_OUT_LV_END) != 0
    }
    pub fn set_out_lv_end(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_OUT_LV_END;
        } else {
            self.my_flags &= !EMASK_OUT_LV_END;
        }
    }

    /// OCCT CutAtSta (lxx).
    pub fn cut_at_sta(&self) -> bool {
        (self.my_flags & EMASK_CUT_AT_STA) != 0
    }
    pub fn set_cut_at_sta(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_CUT_AT_STA;
        } else {
            self.my_flags &= !EMASK_CUT_AT_STA;
        }
    }

    /// OCCT CutAtEnd (lxx).
    pub fn cut_at_end(&self) -> bool {
        (self.my_flags & EMASK_CUT_AT_END) != 0
    }
    pub fn set_cut_at_end(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_CUT_AT_END;
        } else {
            self.my_flags &= !EMASK_CUT_AT_END;
        }
    }

    /// OCCT VerAtSta (lxx).
    pub fn ver_at_sta(&self) -> bool {
        (self.my_flags & EMASK_VER_AT_STA) != 0
    }
    pub fn set_ver_at_sta(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_VER_AT_STA;
        } else {
            self.my_flags &= !EMASK_VER_AT_STA;
        }
    }

    /// OCCT VerAtEnd (lxx).
    pub fn ver_at_end(&self) -> bool {
        (self.my_flags & EMASK_VER_AT_END) != 0
    }
    pub fn set_ver_at_end(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_VER_AT_END;
        } else {
            self.my_flags &= !EMASK_VER_AT_END;
        }
    }

    /// OCCT AutoIntersectionDone (lxx).
    pub fn auto_intersection_done(&self) -> bool {
        (self.my_flags & EMASK_INT_DONE) != 0
    }
    pub fn set_auto_intersection_done(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_INT_DONE;
        } else {
            self.my_flags &= !EMASK_INT_DONE;
        }
    }

    /// OCCT Used (lxx).
    pub fn used(&self) -> bool {
        (self.my_flags & EMASK_USED) != 0
    }
    pub fn set_used(&mut self, b: bool) {
        if b {
            self.my_flags |= EMASK_USED;
        } else {
            self.my_flags &= !EMASK_USED;
        }
    }

    /// OCCT HideCount (lxx).
    pub fn hide_count(&self) -> i32 {
        self.my_hide_count
    }
    pub fn set_hide_count(&mut self, i: i32) {
        self.my_hide_count = i;
    }

    /// OCCT VSta (lxx).
    pub fn v_sta(&self) -> i32 {
        self.my_v_sta
    }
    pub fn set_v_sta(&mut self, i: i32) {
        self.my_v_sta = i;
    }

    /// OCCT VEnd (lxx).
    pub fn v_end(&self) -> i32 {
        self.my_v_end
    }
    pub fn set_v_end(&mut self, i: i32) {
        self.my_v_end = i;
    }

    /// OCCT UpdateMinMax (hxx L121-124).
    pub fn update_min_max(&mut self, the_tot_min_max: &MinMaxIndices) {
        self.my_min_max = the_tot_min_max.clone();
    }

    /// OCCT MinMax (hxx L126-128).
    pub fn min_max(&mut self) -> &mut MinMaxIndices {
        &mut self.my_min_max
    }

    /// OCCT Status (lxx).
    pub fn status(&mut self) -> &mut EdgeStatus {
        &mut self.my_status
    }

    /// OCCT Status (const) — the Bounds read path of the Intersector.
    pub fn status_ref(&self) -> &EdgeStatus {
        &self.my_status
    }

    /// OCCT ChangeGeometry (lxx).
    pub fn change_geometry(&mut self) -> &mut Curve<'a> {
        &mut self.my_geometry
    }

    /// OCCT Geometry (lxx).
    pub fn geometry(&self) -> &Curve<'a> {
        &self.my_geometry
    }

    /// OCCT Curve (hxx L163-165).
    pub fn curve(&mut self) -> &mut Curve<'a> {
        &mut self.my_geometry
    }

    /// OCCT Tolerance (lxx).
    pub fn tolerance(&self) -> f32 {
        self.my_tolerance
    }
}

impl Default for EdgeData<'_> {
    fn default() -> Self {
        EdgeData::new()
    }
}
