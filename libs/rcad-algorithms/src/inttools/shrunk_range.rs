use crate::bopds::ds::DS;
use crate::bopds::pave::PaveBlock;
use rcad_kernel::geom::Curve3;
use super::curve_range::shrunk_range;

/// ✅ OCCT-aligned: IntTools_ShrunkRange — compute the working (shrunk) range
///   for a 3D curve of an edge.  OCCT IntTools_ShrunkRange.hxx.
///   rcad: operates on DS index types instead of TopoDS handles.
pub struct ShrunkRange {
    edge_idx: usize,
    t1: f64,
    t2: f64,
    v1_tol: f64,
    v2_tol: f64,
    edge_tol: f64,
    ts1: f64,
    ts2: f64,
    is_done: bool,
    is_splittable: bool,
}

impl ShrunkRange {
    /// Empty constructor.
    pub fn new() -> Self {
        ShrunkRange {
            edge_idx: usize::MAX,
            t1: 0.0, t2: 0.0,
            v1_tol: 0.0, v2_tol: 0.0, edge_tol: 0.0,
            ts1: 0.0, ts2: 0.0,
            is_done: false, is_splittable: false,
        }
    }

    /// OCCT: SetData(theEdge, t1, t2, theV1, theV2) — set edge and vertex tolerances.
    ///   rcad: takes edge DS index and vertex tolerances directly.
    pub fn set_data(&mut self, edge_idx: usize, t_range: [f64; 2], v1_tol: f64, v2_tol: f64, edge_tol: f64) {
        self.edge_idx = edge_idx;
        self.t1 = t_range[0];
        self.t2 = t_range[1];
        self.v1_tol = v1_tol;
        self.v2_tol = v2_tol;
        self.edge_tol = edge_tol;
    }

    /// Convenience: set data from a DS edge and its PaveBlock.
    pub fn set_data_from_pb(&mut self, ds: &DS, ei: usize, pb: &PaveBlock) {
        let v1_tol = ds.vertices[pb.pave1.vertex_idx].geom_tol;
        let v2_tol = ds.vertices[pb.pave2.vertex_idx].geom_tol;
        let edge_tol = ds.edges[ei].geom_tol;
        self.set_data(ei, [pb.pave1.param, pb.pave2.param], v1_tol, v2_tol, edge_tol);
    }

    /// OCCT: Perform() — compute the shrunk range.
    pub fn perform(&mut self, curve: &Curve3) {
        match shrunk_range(curve, [self.t1, self.t2], self.v1_tol, self.v2_tol, self.edge_tol) {
            Some([ts1, ts2]) => {
                self.ts1 = ts1;
                self.ts2 = ts2;
                self.is_done = true;
                self.is_splittable = (ts2 - ts1) > 2.0 * self.edge_tol + 2.0 * crate::tolerance::TOLERANCE_ABS;
            }
            None => {
                self.is_done = false;
                self.is_splittable = false;
            }
        }
    }

    /// OCCT: IsDone() — true if shrunk range was computed.
    pub fn is_done(&self) -> bool { self.is_done }

    /// OCCT: IsSplittable() — true if edge can be split (shrunk range large enough).
    pub fn is_splittable(&self) -> bool { self.is_splittable }

    /// OCCT: ShrunkRange(ts1, ts2) — get the computed shrunk range.
    pub fn shrunk_range(&self) -> Option<[f64; 2]> {
        if self.is_done { Some([self.ts1, self.ts2]) } else { None }
    }

    /// Set splittable flag (OCCT SetShrunkData equivalent).
    pub fn set_splittable(&mut self, flag: bool) {
        self.is_splittable = flag;
    }
}

impl Default for ShrunkRange {
    fn default() -> Self { Self::new() }
}
