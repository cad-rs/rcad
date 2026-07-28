//! GeomInt_LineConstructor — splits an IntPatch_Line into
//! valid parameter intervals by classifying vertex-interval midpoints
//! on both face domains.
//!
//! OCCT GeomInt_LineConstructor.hxx / .cxx
//!
//! rcad: domain classification uses IntTools_Context (FClass2d) instead of
//! Adaptor3d_TopolTool.

/// GeomInt_LineConstructor — splits lines by domain-classified
/// vertex intervals.  Stores face indices for domain classification.
pub struct GeomIntLineConstructor {
    f1: usize,
    f2: usize,
}

impl GeomIntLineConstructor {
    /// OCCT L475-476: Load(dom1, dom2, myHS1, myHS2) — store domain tools.
    /// rcad: stores face indices for domain classification via Context.
    pub fn new() -> Self {
        Self { f1: 0, f2: 0 }
    }

    /// Load — initializes with face indices for domain checks.
    pub fn load(&mut self, f1: usize, f2: usize) {
        self.f1 = f1;
        self.f2 = f2;
    }

    pub fn f1(&self) -> usize {
        self.f1
    }
    pub fn f2(&self) -> usize {
        self.f2
    }
}
