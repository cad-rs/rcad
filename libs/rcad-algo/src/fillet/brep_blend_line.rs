//! OCCT BRepBlend_Line (TKFillet/BRepBlend) — 1:1 port of BRepBlend_Line.hxx
//! (L29-121) + BRepBlend_Line.cxx (L19-50) + BRepBlend_Line.lxx (L8-110).
//!
//! Architecture mapping: `NCollection_Sequence<Blend_Point>` maps to
//! `Vec<BlendPoint>` (1-based OCCT access -> 0-based storage); the module
//! re-exports [`IntSurfTypeTrans`] for the IntSurf_TypeTrans fields.

use super::brep_blend_extremity::BRepBlendExtremity;
use super::brep_blend_point::BlendPoint;
pub use crate::geomalgo::int_patch::transitions::TypeTrans as IntSurfTypeTrans;

/// OCCT BRepBlend_Line — instantiation of the form for a single line of
/// walking (BRepBlend_Line.hxx L23).
#[derive(Debug, Clone)]
pub struct BRepBlendLine {
    seqpt: Vec<BlendPoint>,
    tras1: IntSurfTypeTrans,
    tras2: IntSurfTypeTrans,
    stp1: BRepBlendExtremity,
    stp2: BRepBlendExtremity,
    endp1: BRepBlendExtremity,
    endp2: BRepBlendExtremity,
    hass1: bool,
    hass2: bool,
}

impl Default for BRepBlendLine {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepBlendLine {
    /// OCCT BRepBlend_Line() (BRepBlend_Line.cxx L19-30) — empty constructor.
    pub fn new() -> Self {
        BRepBlendLine {
            seqpt: Vec::new(),
            tras1: IntSurfTypeTrans::In,
            tras2: IntSurfTypeTrans::In,
            stp1: BRepBlendExtremity::new(),
            stp2: BRepBlendExtremity::new(),
            endp1: BRepBlendExtremity::new(),
            endp2: BRepBlendExtremity::new(),
            hass1: false,
            hass2: false,
        }
    }

    /// OCCT Clear() (BRepBlend_Line.cxx L32-37).
    pub fn clear(&mut self) {
        self.seqpt.clear();
        self.hass1 = false;
        self.hass2 = false;
    }

    /// OCCT Append(P) (BRepBlend_Line.lxx L9-12).
    pub fn append(&mut self, p: BlendPoint) {
        self.seqpt.push(p);
    }

    /// OCCT Prepend(P) (lxx L14-17).
    pub fn prepend(&mut self, p: BlendPoint) {
        self.seqpt.insert(0, p);
    }

    /// OCCT InsertBefore(Index, P) (lxx L19-22) — 1-based Index.
    pub fn insert_before(&mut self, index: i32, p: BlendPoint) {
        self.seqpt.insert((index - 1) as usize, p);
    }

    /// OCCT Remove(FromIndex, ToIndex) (lxx L24-27) — 1-based inclusive
    /// range.
    pub fn remove(&mut self, from_index: i32, to_index: i32) {
        self.seqpt.drain((from_index - 1) as usize..to_index as usize);
    }

    /// OCCT Set(TranS1, TranS2) (BRepBlend_Line.cxx L39-45).
    pub fn set_transitions(&mut self, tran_s1: IntSurfTypeTrans, tran_s2: IntSurfTypeTrans) {
        self.tras1 = tran_s1;
        self.tras2 = tran_s2;
        self.hass1 = true;
        self.hass2 = true;
    }

    /// OCCT Set(Trans) (BRepBlend_Line.cxx L47-50).
    pub fn set_transition(&mut self, trans: IntSurfTypeTrans) {
        self.tras1 = trans;
        self.hass1 = true;
    }

    /// OCCT SetStartPoints(StartPtOnS1, StartPtOnS2) (lxx L29-35).
    pub fn set_start_points(
        &mut self,
        start_pt1: &BRepBlendExtremity,
        start_pt2: &BRepBlendExtremity,
    ) {
        self.stp1 = start_pt1.clone();
        self.stp2 = start_pt2.clone();
    }

    /// OCCT SetEndPoints(EndPtOnS1, EndPtOnS2) (lxx L37-43).
    pub fn set_end_points(&mut self, end_pt1: &BRepBlendExtremity, end_pt2: &BRepBlendExtremity) {
        self.endp1 = end_pt1.clone();
        self.endp2 = end_pt2.clone();
    }

    /// OCCT NbPoints() (lxx L45-48).
    pub fn nb_points(&self) -> i32 {
        self.seqpt.len() as i32
    }

    /// OCCT Point(Index) (lxx L50-53) — 1-based Index.
    pub fn point(&self, index: i32) -> &BlendPoint {
        &self.seqpt[(index - 1) as usize]
    }

    /// Mutable access for OCCT `line->Point(Index)` mutation patterns.
    pub fn point_mut(&mut self, index: i32) -> &mut BlendPoint {
        &mut self.seqpt[(index - 1) as usize]
    }

    /// OCCT TransitionOnS1() (lxx L55-61).
    pub fn transition_on_s1(&self) -> IntSurfTypeTrans {
        if !self.hass1 {
            panic!("Standard_DomainError: BRepBlend_Line::TransitionOnS1");
        }
        self.tras1
    }

    /// OCCT TransitionOnS2() (lxx L63-69).
    pub fn transition_on_s2(&self) -> IntSurfTypeTrans {
        if !self.hass2 {
            panic!("Standard_DomainError: BRepBlend_Line::TransitionOnS2");
        }
        self.tras2
    }

    /// OCCT StartPointOnFirst() (lxx L71-74).
    pub fn start_point_on_first(&self) -> &BRepBlendExtremity {
        &self.stp1
    }

    /// OCCT StartPointOnSecond() (lxx L76-79).
    pub fn start_point_on_second(&self) -> &BRepBlendExtremity {
        &self.stp2
    }

    /// OCCT EndPointOnFirst() (lxx L81-84).
    pub fn end_point_on_first(&self) -> &BRepBlendExtremity {
        &self.endp1
    }

    /// OCCT EndPointOnSecond() (lxx L86-89).
    pub fn end_point_on_second(&self) -> &BRepBlendExtremity {
        &self.endp2
    }

    /// OCCT TransitionOnS() (lxx L91-94).
    pub fn transition_on_s(&self) -> IntSurfTypeTrans {
        self.tras1
    }
}
