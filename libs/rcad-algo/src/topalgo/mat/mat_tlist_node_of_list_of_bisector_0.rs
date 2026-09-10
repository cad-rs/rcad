//! OCCT MAT_TListNodeOfListOfBisector — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_TListNodeOfListOfBisector.hxx L27-113 (class + inline impls),
//!         MAT_TListNodeOfListOfBisector_0.cxx L26 (Dummy).
//! The hxx inline implementations are merged here with the _0.cxx body.

use std::sync::{Arc, RwLock};

use super::mat_bisector::HandleMatBisector;

/// OCCT `occ::handle<MAT_TListNodeOfListOfBisector>`
pub type HandleMatTListNodeOfListOfBisector = Arc<RwLock<MatTListNodeOfListOfBisector>>;

/// OCCT MAT_TListNodeOfListOfBisector.hxx L27-55
pub struct MatTListNodeOfListOfBisector {
    thenext: Option<HandleMatTListNodeOfListOfBisector>,
    theprevious: Option<HandleMatTListNodeOfListOfBisector>,
    theitem: Option<HandleMatBisector>,
}

impl MatTListNodeOfListOfBisector {
    /// OCCT hxx L61 (inline default constructor)
    pub fn new() -> Self {
        MatTListNodeOfListOfBisector {
            thenext: None,
            theprevious: None,
            theitem: None,
        }
    }

    /// OCCT hxx L65-69 (inline constructor from anitem)
    pub fn new_item(anitem: &HandleMatBisector) -> Self {
        MatTListNodeOfListOfBisector {
            thenext: None,
            theprevious: None,
            theitem: Some(anitem.clone()),
        }
    }

    /// OCCT hxx L73-76 (inline GetItem)
    pub fn get_item(&self) -> Option<HandleMatBisector> {
        self.theitem.clone()
    }

    /// OCCT hxx L80-83 (inline Next getter)
    pub fn next(&self) -> Option<HandleMatTListNodeOfListOfBisector> {
        self.thenext.clone()
    }

    /// OCCT hxx L87-90 (inline Previous getter)
    pub fn previous(&self) -> Option<HandleMatTListNodeOfListOfBisector> {
        self.theprevious.clone()
    }

    /// OCCT hxx L94-97 (inline SetItem)
    pub fn set_item(&mut self, anitem: &HandleMatBisector) {
        self.theitem = Some(anitem.clone());
    }

    /// OCCT hxx L101-105 (inline Next setter; the OCCT handle argument may
    /// be null — MAT_TList Unlink/Permute pass null — hence Option)
    pub fn set_next(&mut self, atlistnode: Option<HandleMatTListNodeOfListOfBisector>) {
        self.thenext = atlistnode;
    }

    /// OCCT hxx L109-113 (inline Previous setter; see set_next)
    pub fn set_previous(&mut self, atlistnode: Option<HandleMatTListNodeOfListOfBisector>) {
        self.theprevious = atlistnode;
    }

    /// OCCT MAT_TListNodeOfListOfBisector_0.cxx L26
    pub fn dummy(&self) {}
}

impl Default for MatTListNodeOfListOfBisector {
    fn default() -> Self {
        Self::new()
    }
}
