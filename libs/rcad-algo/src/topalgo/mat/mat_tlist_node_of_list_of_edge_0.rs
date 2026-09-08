//! OCCT MAT_TListNodeOfListOfEdge — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_TListNodeOfListOfEdge.hxx L27-112 (class + inline impls),
//!         MAT_TListNodeOfListOfEdge_0.cxx L26 (Dummy).
//! The hxx inline implementations are merged here with the _0.cxx body.

use std::sync::{Arc, RwLock};

use super::mat_edge::HandleMatEdge;

/// OCCT `occ::handle<MAT_TListNodeOfListOfEdge>`
pub type HandleMatTListNodeOfListOfEdge = Arc<RwLock<MatTListNodeOfListOfEdge>>;

/// OCCT MAT_TListNodeOfListOfEdge.hxx L27-55
pub struct MatTListNodeOfListOfEdge {
    thenext: Option<HandleMatTListNodeOfListOfEdge>,
    theprevious: Option<HandleMatTListNodeOfListOfEdge>,
    theitem: Option<HandleMatEdge>,
}

impl MatTListNodeOfListOfEdge {
    /// OCCT hxx L61 (inline default constructor)
    pub fn new() -> Self {
        MatTListNodeOfListOfEdge {
            thenext: None,
            theprevious: None,
            theitem: None,
        }
    }

    /// OCCT hxx L65-68 (inline constructor from anitem)
    pub fn new_item(anitem: &HandleMatEdge) -> Self {
        MatTListNodeOfListOfEdge {
            thenext: None,
            theprevious: None,
            theitem: Some(anitem.clone()),
        }
    }

    /// OCCT hxx L72-75 (inline GetItem)
    pub fn get_item(&self) -> Option<HandleMatEdge> {
        self.theitem.clone()
    }

    /// OCCT hxx L79-82 (inline Next getter)
    pub fn next(&self) -> Option<HandleMatTListNodeOfListOfEdge> {
        self.thenext.clone()
    }

    /// OCCT hxx L86-89 (inline Previous getter)
    pub fn previous(&self) -> Option<HandleMatTListNodeOfListOfEdge> {
        self.theprevious.clone()
    }

    /// OCCT hxx L93-96 (inline SetItem)
    pub fn set_item(&mut self, anitem: &HandleMatEdge) {
        self.theitem = Some(anitem.clone());
    }

    /// OCCT hxx L100-104 (inline Next setter; the OCCT handle argument may
    /// be null — MAT_TList Unlink/Permute pass null — hence Option)
    pub fn set_next(&mut self, atlistnode: Option<HandleMatTListNodeOfListOfEdge>) {
        self.thenext = atlistnode;
    }

    /// OCCT hxx L108-112 (inline Previous setter; see set_next)
    pub fn set_previous(&mut self, atlistnode: Option<HandleMatTListNodeOfListOfEdge>) {
        self.theprevious = atlistnode;
    }

    /// OCCT MAT_TListNodeOfListOfEdge_0.cxx L26
    pub fn dummy(&self) {}
}

impl Default for MatTListNodeOfListOfEdge {
    fn default() -> Self {
        Self::new()
    }
}
