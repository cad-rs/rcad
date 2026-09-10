//! OCCT MAT_ListOfEdge — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_ListOfEdge.hxx L28-94 (class declaration),
//!         MAT_TList.lxx L19-36 (inline Number/Index/IsEmpty),
//!         MAT_ListOfEdge_0.cxx -> MAT_TList.gxx L19-348 expansion
//!         (all method bodies, including the ~MAT_TList destructor).
//! The _0.cxx `#include <MAT_TList.gxx>` instantiation is expanded inline
//! here for the Item = occ::handle<MAT_Edge> instantiation.
//! OCCT `operator()(anindex)` is C++ operator sugar delegating to
//! Brackets(); the Rust entry point is `brackets`.

use std::sync::{Arc, RwLock};

use super::mat_edge::HandleMatEdge;
use super::mat_tlist_node_of_list_of_edge_0::{
    HandleMatTListNodeOfListOfEdge, MatTListNodeOfListOfEdge,
};

/// OCCT `occ::handle<MAT_ListOfEdge>`
pub type HandleMatListOfEdge = Arc<RwLock<MatListOfEdge>>;

/// OCCT MAT_ListOfEdge.hxx L28-94 (MAT_TList for Item = MAT_Edge)
pub struct MatListOfEdge {
    thefirstnode: Option<HandleMatTListNodeOfListOfEdge>,
    thelastnode: Option<HandleMatTListNodeOfListOfEdge>,
    thecurrentnode: Option<HandleMatTListNodeOfListOfEdge>,
    thecurrentindex: i32,
    thenumberofitems: i32,
}

impl MatListOfEdge {
    /// OCCT MAT_TList.gxx L19-23
    pub fn new() -> Self {
        MatListOfEdge {
            thefirstnode: None,
            thelastnode: None,
            thecurrentnode: None,
            thecurrentindex: 0,
            thenumberofitems: 0,
        }
    }

    /// OCCT MAT_TList.gxx L27-31
    pub fn first(&mut self) {
        self.thecurrentnode = self.thefirstnode.clone();
        self.thecurrentindex = 1;
    }

    /// OCCT MAT_TList.gxx L35-39
    pub fn last(&mut self) {
        self.thecurrentnode = self.thelastnode.clone();
        self.thecurrentindex = self.thenumberofitems;
    }

    /// OCCT MAT_TList.gxx L43-52
    pub fn init(&mut self, anitem: &HandleMatEdge) {
        self.first();
        while self.more() {
            if super::same_handle(
                &self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Init")
                    .read()
                    .unwrap()
                    .get_item(),
                &Some(anitem.clone()),
            ) {
                break;
            }
            self.next();
        }
    }

    /// OCCT MAT_TList.gxx L56-63
    pub fn next(&mut self) {
        if !self.is_empty() {
            let n = self
                .thecurrentnode
                .as_ref()
                .expect("MAT_TList::Next")
                .read()
                .unwrap()
                .next();
            self.thecurrentnode = n;
            self.thecurrentindex = (self.thecurrentindex % self.thenumberofitems) + 1;
        }
    }

    /// OCCT MAT_TList.gxx L67-74
    pub fn previous(&mut self) {
        if !self.is_empty() {
            let p = self
                .thecurrentnode
                .as_ref()
                .expect("MAT_TList::Previous")
                .read()
                .unwrap()
                .previous();
            self.thecurrentnode = p;
            self.thecurrentindex =
                ((self.thecurrentindex + self.thenumberofitems - 2) % self.thenumberofitems) + 1;
        }
    }

    /// OCCT MAT_TList.gxx L78-81
    pub fn more(&self) -> bool {
        !self.thecurrentnode.is_none()
    }

    /// OCCT MAT_TList.gxx L85-88 (Current getter)
    pub fn current(&self) -> Option<HandleMatEdge> {
        self.thecurrentnode
            .as_ref()
            .expect("MAT_TList::Current")
            .read()
            .unwrap()
            .get_item()
    }

    /// OCCT MAT_TList.gxx L92-95 (Current setter)
    pub fn set_current(&self, anitem: &HandleMatEdge) {
        self.thecurrentnode
            .as_ref()
            .expect("MAT_TList::Current")
            .write()
            .unwrap()
            .set_item(anitem);
    }

    /// OCCT MAT_TList.gxx L99-102
    pub fn first_item(&self) -> Option<HandleMatEdge> {
        self.thefirstnode
            .as_ref()
            .expect("MAT_TList::FirstItem")
            .read()
            .unwrap()
            .get_item()
    }

    /// OCCT MAT_TList.gxx L106-109
    pub fn last_item(&self) -> Option<HandleMatEdge> {
        self.thelastnode
            .as_ref()
            .expect("MAT_TList::LastItem")
            .read()
            .unwrap()
            .get_item()
    }

    /// OCCT MAT_TList.gxx L113-116
    pub fn previous_item(&self) -> Option<HandleMatEdge> {
        let prev = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::PreviousItem")
            .read()
            .unwrap()
            .previous();
        let item = prev
            .as_ref()
            .expect("MAT_TList::PreviousItem")
            .read()
            .unwrap()
            .get_item();
        item
    }

    /// OCCT MAT_TList.gxx L120-123
    pub fn next_item(&self) -> Option<HandleMatEdge> {
        let next = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::NextItem")
            .read()
            .unwrap()
            .next();
        let item = next
            .as_ref()
            .expect("MAT_TList::NextItem")
            .read()
            .unwrap()
            .get_item();
        item
    }

    /// OCCT MAT_TList.gxx L127-146
    pub fn brackets(&mut self, anindex: i32) -> Option<HandleMatEdge> {
        if self.thecurrentindex > anindex {
            while self.thecurrentindex != anindex {
                self.thecurrentindex -= 1;
                let prev = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Brackets")
                    .read()
                    .unwrap()
                    .previous();
                self.thecurrentnode = prev;
            }
        } else if self.thecurrentindex < anindex {
            while self.thecurrentindex != anindex {
                self.thecurrentindex += 1;
                let next = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Brackets")
                    .read()
                    .unwrap()
                    .next();
                self.thecurrentnode = next;
            }
        }
        self.thecurrentnode
            .as_ref()
            .expect("MAT_TList::Brackets")
            .read()
            .unwrap()
            .get_item()
    }

    /// OCCT MAT_TList.gxx L150-177
    pub fn unlink(&mut self) {
        let previousisnull = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::Unlink")
            .read()
            .unwrap()
            .previous()
            .is_none();
        let nextisnull = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::Unlink")
            .read()
            .unwrap()
            .next()
            .is_none();

        if self.thecurrentindex != 0 {
            if !nextisnull {
                let next = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Unlink")
                    .read()
                    .unwrap()
                    .next();
                let prev = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Unlink")
                    .read()
                    .unwrap()
                    .previous();
                next.as_ref().expect("MAT_TList::Unlink")
                    .write()
                    .unwrap()
                    .set_previous(prev);
            }
            if !previousisnull {
                let prev = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Unlink")
                    .read()
                    .unwrap()
                    .previous();
                let next = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Unlink")
                    .read()
                    .unwrap()
                    .next();
                prev.as_ref().expect("MAT_TList::Unlink")
                    .write()
                    .unwrap()
                    .set_next(next);
            }

            if self.thecurrentindex == 1 {
                let n = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Unlink")
                    .read()
                    .unwrap()
                    .next();
                self.thefirstnode = n;
            } else if self.thecurrentindex == self.thenumberofitems {
                let p = self
                    .thecurrentnode
                    .as_ref()
                    .expect("MAT_TList::Unlink")
                    .read()
                    .unwrap()
                    .previous();
                self.thelastnode = p;
            }
        }
        self.thenumberofitems -= 1;
        self.thecurrentindex -= 1;
    }

    /// OCCT MAT_TList.gxx L181-205
    pub fn link_before(&mut self, anitem: &HandleMatEdge) {
        self.thenumberofitems += 1;
        if self.thecurrentindex != 0 {
            self.thecurrentindex += 1;
        }

        let mut previous: Option<HandleMatTListNodeOfListOfEdge> = None;

        let node = Arc::new(RwLock::new(MatTListNodeOfListOfEdge::new_item(anitem)));

        let cur_prev = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::LinkBefore")
            .read()
            .unwrap()
            .previous();
        if cur_prev.is_some() {
            previous = cur_prev;
            previous
                .as_ref()
                .expect("MAT_TList::LinkBefore")
                .write()
                .unwrap()
                .set_next(Some(node.clone()));
            node.write().unwrap().set_previous(previous.clone());
        }

        if self.thecurrentindex == 2 {
            self.thefirstnode = Some(node.clone());
        }

        self.thecurrentnode
            .as_ref()
            .expect("MAT_TList::LinkBefore")
            .write()
            .unwrap()
            .set_previous(Some(node.clone()));
        node.write().unwrap().set_next(self.thecurrentnode.clone());
    }

    /// OCCT MAT_TList.gxx L209-230
    pub fn link_after(&mut self, anitem: &HandleMatEdge) {
        self.thenumberofitems += 1;
        let mut next: Option<HandleMatTListNodeOfListOfEdge> = None;

        let node = Arc::new(RwLock::new(MatTListNodeOfListOfEdge::new_item(anitem)));

        let cur_next = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::LinkAfter")
            .read()
            .unwrap()
            .next();
        if cur_next.is_some() {
            next = cur_next;
            next.as_ref()
                .expect("MAT_TList::LinkAfter")
                .write()
                .unwrap()
                .set_previous(Some(node.clone()));
            node.write().unwrap().set_next(next.clone());
        }

        if self.thecurrentindex + 1 == self.thenumberofitems {
            self.thelastnode = Some(node.clone());
        }

        self.thecurrentnode
            .as_ref()
            .expect("MAT_TList::LinkAfter")
            .write()
            .unwrap()
            .set_next(Some(node.clone()));
        node.write().unwrap().set_previous(self.thecurrentnode.clone());
    }

    /// OCCT MAT_TList.gxx L234-253
    pub fn front_add(&mut self, anitem: &HandleMatEdge) {
        self.thenumberofitems += 1;
        if self.thecurrentindex != 0 {
            self.thecurrentindex += 1;
        }

        let node = Arc::new(RwLock::new(MatTListNodeOfListOfEdge::new_item(anitem)));

        if !self.thefirstnode.is_none() {
            self.thefirstnode
                .as_ref()
                .expect("MAT_TList::FrontAdd")
                .write()
                .unwrap()
                .set_previous(Some(node.clone()));
            node.write().unwrap().set_next(self.thefirstnode.clone());
        } else {
            self.thelastnode = Some(node.clone());
        }

        self.thefirstnode = Some(node);
    }

    /// OCCT MAT_TList.gxx L257-273
    pub fn back_add(&mut self, anitem: &HandleMatEdge) {
        self.thenumberofitems += 1;
        let node = Arc::new(RwLock::new(MatTListNodeOfListOfEdge::new_item(anitem)));

        if !self.thelastnode.is_none() {
            self.thelastnode
                .as_ref()
                .expect("MAT_TList::BackAdd")
                .write()
                .unwrap()
                .set_next(Some(node.clone()));
            node.write().unwrap().set_previous(self.thelastnode.clone());
        } else {
            self.thefirstnode = Some(node.clone());
        }

        self.thelastnode = Some(node);
    }

    /// OCCT MAT_TList.gxx L277-310
    pub fn permute(&mut self) {
        let previous = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::Permute")
            .read()
            .unwrap()
            .previous();
        let current = self.thecurrentnode.clone().expect("MAT_TList::Permute");
        let next = self
            .thecurrentnode
            .as_ref()
            .expect("MAT_TList::Permute")
            .read()
            .unwrap()
            .next()
            .expect("MAT_TList::Permute");
        let nextnext = next.read().unwrap().next();
        let null: Option<HandleMatTListNodeOfListOfEdge> = None;

        if let Some(prev) = previous.as_ref() {
            prev.write().unwrap().set_next(Some(next.clone()));
            next.write().unwrap().set_previous(previous.clone());
        } else {
            next.write().unwrap().set_previous(null.clone());
        }
        next.write().unwrap().set_next(Some(current.clone()));
        current.write().unwrap().set_previous(Some(next.clone()));
        if let Some(nn) = nextnext.as_ref() {
            current.write().unwrap().set_next(nextnext.clone());
            nn.write().unwrap().set_previous(Some(current.clone()));
        } else {
            current.write().unwrap().set_next(null.clone());
        }
        if super::same_handle(&self.thefirstnode, &Some(current.clone())) {
            self.thefirstnode = Some(next.clone());
        }
        if super::same_handle(&self.thelastnode, &Some(next.clone())) {
            self.thelastnode = Some(current.clone());
        }
        self.thecurrentindex += 1;
    }

    /// OCCT MAT_TList.gxx L314-318 (OCCT name: Loop)
    pub fn r#loop(&self) {
        self.thelastnode
            .as_ref()
            .expect("MAT_TList::Loop")
            .write()
            .unwrap()
            .set_next(self.thefirstnode.clone());
        self.thefirstnode
            .as_ref()
            .expect("MAT_TList::Loop")
            .write()
            .unwrap()
            .set_previous(self.thelastnode.clone());
    }

    /// OCCT MAT_TList.gxx L322-326
    pub fn dump(&mut self, ashift: i32, alevel: i32) {
        self.first();
        while self.more() {
            self.current()
                .expect("MAT_TList::Dump")
                .read()
                .unwrap()
                .dump(ashift, alevel);
            self.next();
        }
    }

    /// OCCT MAT_TList.lxx L19-22 (Number)
    pub fn number(&self) -> i32 {
        self.thenumberofitems
    }

    /// OCCT MAT_TList.lxx L26-29 (Index)
    pub fn index(&self) -> i32 {
        self.thecurrentindex
    }

    /// OCCT MAT_TList.lxx L33-36 (IsEmpty)
    pub fn is_empty(&self) -> bool {
        self.thenumberofitems == 0
    }
}

/// OCCT MAT_TList.gxx L328-348 (~MAT_TList destructor)
impl Drop for MatListOfEdge {
    fn drop(&mut self) {
        let mut a_node = self.thefirstnode.take();
        while let Some(node) = a_node {
            let a_next = node.read().unwrap().next();
            // aNode->Next(NULL); aNode->Previous(NULL);
            node.write().unwrap().set_next(None);
            node.write().unwrap().set_previous(None);
            a_node = a_next;
        }
        self.thecurrentnode = None; // thecurrentnode.Nullify()
        self.thefirstnode = None; // thefirstnode.Nullify()
        self.thelastnode = None; // thelastnode.Nullify()
        self.thecurrentindex = 0;
        self.thenumberofitems = 0;
    }
}
