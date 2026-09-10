//! OCCT HLRAlgo_BiPoint (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_BiPoint.hxx` (L29-260) +
//! `HLRAlgo_BiPoint.cxx` (L25-240).

use glam::{DVec2, DVec3};

/// OCCT HLRAlgo_BiPoint::IndicesT (hxx L32-58).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndicesT {
    pub shape_index: i32,
    pub face_conex1: i32,
    pub face1_pt1: i32,
    pub face1_pt2: i32,
    pub face_conex2: i32,
    pub face2_pt1: i32,
    pub face2_pt2: i32,
    pub min_seg: i32,
    pub max_seg: i32,
    pub seg_flags: i32,
}

impl Default for IndicesT {
    fn default() -> Self {
        IndicesT {
            shape_index: -1,
            face_conex1: 0,
            face1_pt1: 0,
            face1_pt2: 0,
            face_conex2: 0,
            face2_pt1: 0,
            face2_pt2: 0,
            min_seg: 0,
            max_seg: 0,
            seg_flags: 0,
        }
    }
}

/// OCCT HLRAlgo_BiPoint::PointsT (hxx L60-70).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PointsT {
    pub pnt1: DVec3,
    pub pnt2: DVec3,
    pub pnt_p1: DVec3,
    pub pnt_p2: DVec3,
}

impl PointsT {
    /// OCCT PntP12D() — hxx L67.
    pub fn pnt_p1_2d(&self) -> DVec2 {
        DVec2::new(self.pnt_p1.x, self.pnt_p1.y)
    }

    /// OCCT PntP22D() — hxx L69.
    pub fn pnt_p2_2d(&self) -> DVec2 {
        DVec2::new(self.pnt_p2.x, self.pnt_p2.y)
    }
}

/// OCCT enum EMskFlags (hxx L248-255) — private flag bit masks.
pub const E_MSK_RG1_LINE: i32 = 1;
pub const E_MSK_RGN_LINE: i32 = 2;
pub const E_MSK_OUT_LINE: i32 = 4;
pub const E_MSK_INT_LINE: i32 = 8;
pub const E_MSK_HIDDEN: i32 = 16;

/// OCCT HLRAlgo_BiPoint — a segment with its 3D and projected end points,
/// plus face-connection indices and segment flags.
#[derive(Debug, Clone, Copy, Default)]
pub struct BiPoint {
    my_indices: IndicesT,
    my_points: PointsT,
}

impl BiPoint {
    /// OCCT HLRAlgo_BiPoint(X1..ZT2, Index, reg1, regn, outl, intl) —
    /// cxx L25-56.
    #[allow(clippy::too_many_arguments)]
    pub fn from_flags_bool(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        let mut bp = BiPoint {
            my_indices: IndicesT::default(),
            my_points: PointsT {
                pnt1: DVec3::new(x1, y1, z1),
                pnt2: DVec3::new(x2, y2, z2),
                pnt_p1: DVec3::new(xt1, yt1, zt1),
                pnt_p2: DVec3::new(xt2, yt2, zt2),
            },
        };
        bp.my_indices.shape_index = index;
        bp.my_indices.seg_flags = 0;
        bp.set_rg1_line(reg1);
        bp.set_rgn_line(regn);
        bp.set_out_line(outl);
        bp.set_int_line(intl);
        bp.set_hidden(false);
        bp
    }

    /// OCCT HLRAlgo_BiPoint(X1..ZT2, Index, flag) — cxx L60-84.
    #[allow(clippy::too_many_arguments)]
    pub fn from_flags_mask(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        flag: i32,
    ) -> Self {
        let mut bp = BiPoint {
            my_indices: IndicesT::default(),
            my_points: PointsT {
                pnt1: DVec3::new(x1, y1, z1),
                pnt2: DVec3::new(x2, y2, z2),
                pnt_p1: DVec3::new(xt1, yt1, zt1),
                pnt_p2: DVec3::new(xt2, yt2, zt2),
            },
        };
        bp.my_indices.shape_index = index;
        bp.my_indices.seg_flags = flag;
        bp.set_hidden(false);
        bp
    }

    /// OCCT HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, reg1, regn,
    /// outl, intl) — cxx L88-124.
    #[allow(clippy::too_many_arguments)]
    pub fn from_one_face_bool(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        let mut bp = BiPoint {
            my_indices: IndicesT::default(),
            my_points: PointsT {
                pnt1: DVec3::new(x1, y1, z1),
                pnt2: DVec3::new(x2, y2, z2),
                pnt_p1: DVec3::new(xt1, yt1, zt1),
                pnt_p2: DVec3::new(xt2, yt2, zt2),
            },
        };
        bp.my_indices.shape_index = index;
        bp.my_indices.face_conex1 = i1;
        bp.my_indices.face1_pt1 = i1p1;
        bp.my_indices.face1_pt2 = i1p2;
        bp.my_indices.seg_flags = 0;
        bp.set_rg1_line(reg1);
        bp.set_rgn_line(regn);
        bp.set_out_line(outl);
        bp.set_int_line(intl);
        bp.set_hidden(false);
        bp
    }

    /// OCCT HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, flag) —
    /// cxx L128-157.
    #[allow(clippy::too_many_arguments)]
    pub fn from_one_face_mask(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        flag: i32,
    ) -> Self {
        let mut bp = BiPoint {
            my_indices: IndicesT::default(),
            my_points: PointsT {
                pnt1: DVec3::new(x1, y1, z1),
                pnt2: DVec3::new(x2, y2, z2),
                pnt_p1: DVec3::new(xt1, yt1, zt1),
                pnt_p2: DVec3::new(xt2, yt2, zt2),
            },
        };
        bp.my_indices.shape_index = index;
        bp.my_indices.face_conex1 = i1;
        bp.my_indices.face1_pt1 = i1p1;
        bp.my_indices.face1_pt2 = i1p2;
        bp.my_indices.seg_flags = flag;
        bp.set_hidden(false);
        bp
    }

    /// OCCT HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, i2, i2p1, i2p2,
    /// reg1, regn, outl, intl) — cxx L161-202.
    #[allow(clippy::too_many_arguments)]
    pub fn from_two_faces_bool(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        i2: i32,
        i2p1: i32,
        i2p2: i32,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        let mut bp = BiPoint {
            my_indices: IndicesT::default(),
            my_points: PointsT {
                pnt1: DVec3::new(x1, y1, z1),
                pnt2: DVec3::new(x2, y2, z2),
                pnt_p1: DVec3::new(xt1, yt1, zt1),
                pnt_p2: DVec3::new(xt2, yt2, zt2),
            },
        };
        bp.my_indices.shape_index = index;
        bp.my_indices.face_conex1 = i1;
        bp.my_indices.face1_pt1 = i1p1;
        bp.my_indices.face1_pt2 = i1p2;
        bp.my_indices.face_conex2 = i2;
        bp.my_indices.face2_pt1 = i2p1;
        bp.my_indices.face2_pt2 = i2p2;
        bp.my_indices.seg_flags = 0;
        bp.set_rg1_line(reg1);
        bp.set_rgn_line(regn);
        bp.set_out_line(outl);
        bp.set_int_line(intl);
        bp.set_hidden(false);
        bp
    }

    /// OCCT HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, i2, i2p1, i2p2,
    /// flag) — cxx L206-240.
    #[allow(clippy::too_many_arguments)]
    pub fn from_two_faces_mask(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        i2: i32,
        i2p1: i32,
        i2p2: i32,
        flag: i32,
    ) -> Self {
        let mut bp = BiPoint {
            my_indices: IndicesT::default(),
            my_points: PointsT {
                pnt1: DVec3::new(x1, y1, z1),
                pnt2: DVec3::new(x2, y2, z2),
                pnt_p1: DVec3::new(xt1, yt1, zt1),
                pnt_p2: DVec3::new(xt2, yt2, zt2),
            },
        };
        bp.my_indices.shape_index = index;
        bp.my_indices.face_conex1 = i1;
        bp.my_indices.face1_pt1 = i1p1;
        bp.my_indices.face1_pt2 = i1p2;
        bp.my_indices.face_conex2 = i2;
        bp.my_indices.face2_pt1 = i2p1;
        bp.my_indices.face2_pt2 = i2p2;
        bp.my_indices.seg_flags = flag;
        bp.set_hidden(false);
        bp
    }

    /// OCCT Rg1Line() — hxx L193.
    pub fn rg1_line(&self) -> bool {
        (self.my_indices.seg_flags & E_MSK_RG1_LINE) != 0
    }

    /// OCCT Rg1Line(B) — hxx L195-201.
    pub fn set_rg1_line(&mut self, b: bool) {
        if b {
            self.my_indices.seg_flags |= E_MSK_RG1_LINE;
        } else {
            self.my_indices.seg_flags &= !E_MSK_RG1_LINE;
        }
    }

    /// OCCT RgNLine() — hxx L203.
    pub fn rgn_line(&self) -> bool {
        (self.my_indices.seg_flags & E_MSK_RGN_LINE) != 0
    }

    /// OCCT RgNLine(B) — hxx L205-211.
    pub fn set_rgn_line(&mut self, b: bool) {
        if b {
            self.my_indices.seg_flags |= E_MSK_RGN_LINE;
        } else {
            self.my_indices.seg_flags &= !E_MSK_RGN_LINE;
        }
    }

    /// OCCT OutLine() — hxx L213.
    pub fn out_line(&self) -> bool {
        (self.my_indices.seg_flags & E_MSK_OUT_LINE) != 0
    }

    /// OCCT OutLine(B) — hxx L215-221.
    pub fn set_out_line(&mut self, b: bool) {
        if b {
            self.my_indices.seg_flags |= E_MSK_OUT_LINE;
        } else {
            self.my_indices.seg_flags &= !E_MSK_OUT_LINE;
        }
    }

    /// OCCT IntLine() — hxx L223.
    pub fn int_line(&self) -> bool {
        (self.my_indices.seg_flags & E_MSK_INT_LINE) != 0
    }

    /// OCCT IntLine(B) — hxx L225-231.
    pub fn set_int_line(&mut self, b: bool) {
        if b {
            self.my_indices.seg_flags |= E_MSK_INT_LINE;
        } else {
            self.my_indices.seg_flags &= !E_MSK_INT_LINE;
        }
    }

    /// OCCT Hidden() — hxx L233.
    pub fn hidden(&self) -> bool {
        (self.my_indices.seg_flags & E_MSK_HIDDEN) != 0
    }

    /// OCCT Hidden(B) — hxx L235-241.
    pub fn set_hidden(&mut self, b: bool) {
        if b {
            self.my_indices.seg_flags |= E_MSK_HIDDEN;
        } else {
            self.my_indices.seg_flags &= !E_MSK_HIDDEN;
        }
    }

    /// OCCT Indices() — hxx L243.
    pub fn indices(&mut self) -> &mut IndicesT {
        &mut self.my_indices
    }

    /// OCCT Points() — hxx L245.
    pub fn points(&mut self) -> &mut PointsT {
        &mut self.my_points
    }

    /// Immutable accessors (OCCT returns references through non-const
    /// methods; the const reads in the algorithms map to these).
    pub fn indices_ref(&self) -> &IndicesT {
        &self.my_indices
    }

    /// Immutable points access (see `indices_ref`).
    pub fn points_ref(&self) -> &PointsT {
        &self.my_points
    }
}
