//! OCCT HLRBRep_EdgeIList (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_EdgeIList.hxx` (L29-42) + `.cxx` (L29-182).
//! Static tools over the interference list: `NCollection_List<
//! HLRAlgo_Interference>` maps to `Vec<Interference>` (the landed HLRAlgo
//! precedent), and the list iterators map to index positions.  The
//! `OCCT_DEBUG_SI` blocks (cxx L49-68, L85-87, L103-106) are compiled out
//! in OCCT when the macro is not defined and are not ported; the trailing
//! commented-out block (cxx L130-181) is a comment in OCCT.

use glam::DVec3;

use crate::hlr::algo::interference::Interference;
use crate::hlr::top_cnx::EdgeFaceTransition;

use super::edge_interference_tool::{Data, EdgeInterferenceTool};

/// OCCT HLRBRep_EdgeIList.
pub struct EdgeIList;

impl EdgeIList {
    /// OCCT AddInterference(IL, I, T) — cxx L29-45: insert an interference
    /// in a sorted list.
    pub fn add_interference<D: Data + ?Sized>(
        il: &mut Vec<Interference>,
        i: &Interference,
        t: &EdgeInterferenceTool<'_, D>,
    ) {
        // NCollection_List<HLRAlgo_Interference>::Iterator It(IL) — the
        // iterator position is the index.
        let p = t.parameter_of_interference(i);
        let mut it = 0usize;
        while it < il.len() {
            if p < t.parameter_of_interference(&il[it]) {
                // IL.InsertBefore(I, It)
                il.insert(it, *i);
                return;
            }
            it += 1;
        }
        il.push(*i);
    }

    /// OCCT ProcessComplex(IL, T) — cxx L69-128: process complex
    /// transitions on the list IL.
    pub fn process_complex<D: Data + ?Sized>(
        il: &mut Vec<Interference>,
        t: &EdgeInterferenceTool<'_, D>,
    ) {
        let mut trans_tool = EdgeFaceTransition::new();
        // OCCT gp_Dir locals TgtE, NormE, TgtI, NormI (uninitialized); the
        // neutral default is 0.
        let mut tgt_e = DVec3::ZERO;
        let mut norm_e = DVec3::ZERO;
        let mut tgt_i = DVec3::ZERO;
        let mut norm_i = DVec3::ZERO;
        // OCCT double CurvE, CurvI (uninitialized); the neutral default is 0.
        let mut curv_e: f64 = 0.0;
        let mut curv_i: f64 = 0.0;
        const TOL_ANG: f64 = 0.0001;
        let mut it1 = 0usize;

        while it1 < il.len() {
            // NCollection_List<HLRAlgo_Interference>::Iterator It2(It1);
            // It2.Next();
            let it2 = it1 + 1;
            if it2 < il.len() {
                if t.same_interferences(&il[it1], &il[it2]) {
                    t.edge_geometry(
                        t.parameter_of_interference(&il[it1]),
                        &mut tgt_e,
                        &mut norm_e,
                        &mut curv_e,
                    );
                    trans_tool.reset(tgt_e, norm_e, curv_e);
                    t.interference_boundary_geometry(&il[it1], &mut tgt_i, &mut norm_i, &mut curv_i);
                    trans_tool.add_interference(
                        TOL_ANG,
                        tgt_i,
                        norm_i,
                        curv_i,
                        il[it1].orientation(),
                        il[it1].transition(),
                        il[it1].boundary_transition(),
                    );

                    while it2 < il.len() {
                        if !t.same_interferences(&il[it1], &il[it2]) {
                            break;
                        }

                        t.interference_boundary_geometry(&il[it2], &mut tgt_i, &mut norm_i, &mut curv_i);
                        trans_tool.add_interference(
                            TOL_ANG,
                            tgt_i,
                            norm_i,
                            curv_i,
                            il[it2].orientation(),
                            il[it2].transition(),
                            il[it2].boundary_transition(),
                        );
                        // IL.Remove(It2) — removes the current item and
                        // advances the iterator; the Vec index stays on the
                        // shifted-in successor.
                        il.remove(it2);
                    }
                    // get the cumulated results
                    // It1.ChangeValue().Transition(transTool.Transition());
                    il[it1].set_transition(trans_tool.transition());
                    // It1.ChangeValue().BoundaryTransition(transTool.BoundaryTransition());
                    il[it1].set_boundary_transition(trans_tool.boundary_transition());
                }
            }
            it1 += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec2;
    use crate::hlr::algo::coincidence::Coincidence;
    use crate::hlr::algo::intersection::Intersection;
    use crate::hlr::brep::edge_data::EdgeData;
    use crate::hlr::brep::edge_interference_tool::EdgeInterferenceTool;
    use rcad_kernel::topods::Orientation;

    /// Minimal HLRBRep_Data stand-in (one loaded edge, fixed geometries).
    struct TestData {
        edges: Vec<EdgeData<'static>>,
        edge: i32,
    }

    impl TestData {
        fn new() -> Self {
            let mut ed = EdgeData::new();
            ed.set_v_sta(11);
            ed.set_v_end(22);
            ed.status().initialize(0.0, 1.0e-4, 1.0, 1.0e-4);
            TestData {
                edges: vec![ed],
                edge: 1,
            }
        }
    }

    impl Data for TestData {
        fn e_data_array(&self) -> &[EdgeData<'static>] {
            &self.edges
        }
        fn edge(&self) -> i32 {
            self.edge
        }
        fn local_le_geometry_2d(
            &self,
            _param: f64,
            tg_le: &mut DVec2,
            nm_le: &mut DVec2,
            cr_le: &mut f64,
        ) {
            *tg_le = DVec2::new(1.0, 0.0);
            *nm_le = DVec2::new(0.0, 1.0);
            *cr_le = 0.0;
        }
        fn local_fe_geometry_2d(
            &self,
            _fe: i32,
            _param: f64,
            tg_fe: &mut DVec2,
            nm_fe: &mut DVec2,
            cr_fe: &mut f64,
        ) {
            *tg_fe = DVec2::new(1.0, 0.0);
            *nm_fe = DVec2::new(0.0, 1.0);
            *cr_fe = 0.0;
        }
    }

    fn interference_at(param: f64, index: i32, trans: Orientation, b_trans: Orientation) -> Interference {
        let mut inters = Intersection::new();
        inters.set_parameter(param);
        inters.set_index(index);
        let mut bound = Coincidence::new();
        bound.set_2d(1, param);
        Interference::from_parts(inters, bound, Orientation::Forward, trans, b_trans)
    }

    /// OCCT AddInterference (cxx L29-45): the list stays sorted by the
    /// interference parameter.
    #[test]
    fn edge_ilist_add_interference_keeps_sorted_order() {
        let data = TestData::new();
        let tool = EdgeInterferenceTool::new(&data);
        let mut il: Vec<Interference> = Vec::new();
        EdgeIList::add_interference(&mut il, &interference_at(0.5, 1, Orientation::Forward, Orientation::Forward), &tool);
        EdgeIList::add_interference(&mut il, &interference_at(0.2, 2, Orientation::Forward, Orientation::Forward), &tool);
        EdgeIList::add_interference(&mut il, &interference_at(0.8, 3, Orientation::Forward, Orientation::Forward), &tool);
        EdgeIList::add_interference(&mut il, &interference_at(0.35, 4, Orientation::Forward, Orientation::Forward), &tool);
        let params: Vec<f64> = il
            .iter()
            .map(|i| i.intersection().parameter())
            .collect();
        assert_eq!(params, vec![0.2, 0.35, 0.5, 0.8]);
    }

    /// OCCT ProcessComplex (cxx L69-128): two interferences on the same
    /// geometric locus (equal non-zero intersection index) are merged into
    /// one carrying the cumulated transitions.
    #[test]
    fn edge_ilist_process_complex_merges_same_locus() {
        let data = TestData::new();
        let tool = EdgeInterferenceTool::new(&data);
        // Same locus (index 7).  The INTERNAL interference transitions are
        // re-derived inside TopTrans_CurveTransition::Compare from the
        // geometry: parallel tangents + FORWARD orientation promote the
        // INTERNAL transition to FORWARD (top_trans/mod.rs L128-134), so
        // the cumulated states are OUT before / IN after, i.e. a REVERSED
        // transition; two Forward boundary transitions vote Forward.
        let i1 = interference_at(0.4, 7, Orientation::Internal, Orientation::Forward);
        let i2 = interference_at(0.45, 7, Orientation::Internal, Orientation::Forward);
        // Different locus (index 9), untouched.
        let i3 = interference_at(0.7, 9, Orientation::Reversed, Orientation::Reversed);
        let mut il = vec![i1, i2, i3];
        EdgeIList::process_complex(&mut il, &tool);
        assert_eq!(il.len(), 2);
        // The first pair merged into It1 with the cumulated transitions:
        // tran_first = tran_last = FORWARD gives OUT before / IN after,
        // i.e. the FORWARD transition (the curve enters the face).
        assert_eq!(il[0].intersection().parameter(), 0.4);
        assert_eq!(il[0].transition(), Orientation::Forward);
        assert_eq!(il[0].boundary_transition(), Orientation::Forward);
        // The second locus keeps its own transitions.
        assert_eq!(il[1].intersection().parameter(), 0.7);
        assert_eq!(il[1].transition(), Orientation::Reversed);
        assert_eq!(il[1].boundary_transition(), Orientation::Reversed);
    }
}
