//! Typed storage containers for entities, references, and representations.
//!
//! OCCT BRepGraphInc: BRepGraphInc_Storage.hxx
//!
//! Provides DefStore<T>, RefStore<T>, and RepStore<T> with Append/Get/Change,
//! UID allocation, soft-removal flags, and active-count tracking.

use std::marker::PhantomData;
use bitflags::bitflags;

use crate::graph::inc::id::*;

// ── UID (unique ID) allocation ──────────────────────────────────────────────

/// Monotonically increasing unique-ID allocator per entity kind.
#[derive(Debug, Clone)]
pub struct UidAllocator(u64);

impl UidAllocator {
    pub fn new() -> Self { UidAllocator(1) }
    pub fn allocate(&mut self) -> u64 { let v = self.0; self.0 += 1; v }
    pub fn current(&self) -> u64 { self.0 }
    pub fn set_current(&mut self, v: u64) { self.0 = v; }
}

// ── Bit flags (soft removal) ────────────────────────────────────────────────

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct EntityFlags: u8 {
        const NONE   = 0b0000_0000;
        const REMOVED = 0b0000_0001;
    }
}

// ── DefStore ────────────────────────────────────────────────────────────────

/// Typed store for entity definitions.
///
/// Manages a vector of T + per-entry flags + UID allocation.
/// Analogous to OCCT's `DefStore<T>` in BRepGraphInc_Storage.
#[derive(Debug, Clone)]
pub struct DefStore<T> {
    data: Vec<T>,
    flags: Vec<EntityFlags>,
    uids: Vec<u64>,
    uid_alloc: UidAllocator,
    active_count: usize,
}

impl<T> DefStore<T> {
    pub fn new() -> Self {
        DefStore { data: Vec::new(), flags: Vec::new(), uids: Vec::new(), uid_alloc: UidAllocator::new(), active_count: 0 }
    }

    pub fn with_capacity(cap: usize) -> Self {
        DefStore { data: Vec::with_capacity(cap), flags: Vec::with_capacity(cap), uids: Vec::with_capacity(cap), uid_alloc: UidAllocator::new(), active_count: 0 }
    }

    /// Append a definition, allocate its UID and index. Returns the linear index.
    pub fn append(&mut self, def: T) -> usize {
        let idx = self.data.len();
        self.data.push(def);
        self.flags.push(EntityFlags::NONE);
        self.uids.push(self.uid_alloc.allocate());
        self.active_count += 1;
        idx
    }

    pub fn get(&self, idx: usize) -> &T { &self.data[idx] }
    pub fn get_mut(&mut self, idx: usize) -> &mut T { &mut self.data[idx] }
    pub fn uid(&self, idx: usize) -> u64 { self.uids[idx] }
    pub fn flags(&self, idx: usize) -> EntityFlags { self.flags[idx] }
    pub fn set_flags(&mut self, idx: usize, f: EntityFlags) { self.flags[idx] = f; }
    pub fn is_removed(&self, idx: usize) -> bool { self.flags[idx].contains(EntityFlags::REMOVED) }
    pub fn set_removed(&mut self, idx: usize, removed: bool) {
        if removed { self.flags[idx] |= EntityFlags::REMOVED; self.active_count -= 1; }
        else { self.flags[idx] &= !EntityFlags::REMOVED; self.active_count += 1; }
    }
    pub fn len(&self) -> usize { self.data.len() }
    pub fn is_empty(&self) -> bool { self.data.is_empty() }
    pub fn active_count(&self) -> usize { self.active_count }
    pub fn iter(&self) -> impl Iterator<Item = &T> { self.data.iter() }
    pub fn iter_enumerated(&self) -> impl Iterator<Item = (usize, &T)> { self.data.iter().enumerate() }

    /// Reserve space for N definitions with pre-allocated UIDs.
    /// Used for indexed-load preparation.
    pub fn prepare(&mut self, count: usize) {
        self.data.reserve(count);
        self.flags.reserve(count);
        self.uids.reserve(count);
    }
}

// ── RefStore ────────────────────────────────────────────────────────────────

/// Typed store for reference entries.
///
/// Same pattern as DefStore but without per-entry UIDs (refs share the
/// owning entity's generation counter).
#[derive(Debug, Clone)]
pub struct RefStore<T> {
    data: Vec<T>,
    flags: Vec<EntityFlags>,
    active_count: usize,
}

impl<T> RefStore<T> {
    pub fn new() -> Self { RefStore { data: Vec::new(), flags: Vec::new(), active_count: 0 } }
    pub fn with_capacity(cap: usize) -> Self { RefStore { data: Vec::with_capacity(cap), flags: Vec::with_capacity(cap), active_count: 0 } }

    pub fn append(&mut self, entry: T) -> usize {
        let idx = self.data.len();
        self.data.push(entry);
        self.flags.push(EntityFlags::NONE);
        self.active_count += 1;
        idx
    }

    pub fn get(&self, idx: usize) -> &T { &self.data[idx] }
    pub fn get_mut(&mut self, idx: usize) -> &mut T { &mut self.data[idx] }
    pub fn is_removed(&self, idx: usize) -> bool { self.flags[idx].contains(EntityFlags::REMOVED) }
    pub fn set_removed(&mut self, idx: usize, removed: bool) {
        if removed { self.flags[idx] |= EntityFlags::REMOVED; self.active_count -= 1; }
        else { self.flags[idx] &= !EntityFlags::REMOVED; self.active_count += 1; }
    }
    pub fn len(&self) -> usize { self.data.len() }
    pub fn is_empty(&self) -> bool { self.data.is_empty() }
    pub fn active_count(&self) -> usize { self.active_count }
    pub fn iter(&self) -> impl Iterator<Item = &T> { self.data.iter() }
}

// ── RepStore ────────────────────────────────────────────────────────────────

/// Typed store for geometry representations (surfaces, curves, triangulations).
#[derive(Debug, Clone)]
pub struct RepStore<T> {
    data: Vec<T>,
    active_count: usize,
}

impl<T> RepStore<T> {
    pub fn new() -> Self { RepStore { data: Vec::new(), active_count: 0 } }
    pub fn append(&mut self, rep: T) -> usize { let i = self.data.len(); self.data.push(rep); self.active_count += 1; i }
    pub fn get(&self, idx: usize) -> &T { &self.data[idx] }
    pub fn get_mut(&mut self, idx: usize) -> &mut T { &mut self.data[idx] }
    pub fn len(&self) -> usize { self.data.len() }
    pub fn is_empty(&self) -> bool { self.data.is_empty() }
}

// ── Strongly-typed ID wrapper for DefStore access ──────────────────────────

/// Wraps a DefStore with a typed ID for type-safe indexing.
#[derive(Debug, Clone)]
pub struct TypedStore<T, Id> {
    store: DefStore<T>,
    _phantom: PhantomData<Id>,
}

impl<T, Id: Into<u32> + From<u32> + Copy> TypedStore<T, Id> {
    pub fn new() -> Self { TypedStore { store: DefStore::new(), _phantom: PhantomData } }
    pub fn with_capacity(cap: usize) -> Self { TypedStore { store: DefStore::with_capacity(cap), _phantom: PhantomData } }

    pub fn append(&mut self, def: T) -> Id { Id::from(self.store.append(def) as u32) }
    pub fn get(&self, id: Id) -> &T { self.store.get(Into::<u32>::into(id) as usize) }
    pub fn get_mut(&mut self, id: Id) -> &mut T { self.store.get_mut(Into::<u32>::into(id) as usize) }
    pub fn uid(&self, id: Id) -> u64 { self.store.uid(Into::<u32>::into(id) as usize) }
    pub fn flags(&self, id: Id) -> EntityFlags { self.store.flags(Into::<u32>::into(id) as usize) }
    pub fn set_flags(&mut self, id: Id, f: EntityFlags) { self.store.set_flags(Into::<u32>::into(id) as usize, f); }
    pub fn is_removed(&self, id: Id) -> bool { self.store.is_removed(Into::<u32>::into(id) as usize) }
    pub fn set_removed(&mut self, id: Id, removed: bool) { self.store.set_removed(Into::<u32>::into(id) as usize, removed); }
    pub fn len(&self) -> usize { self.store.len() }
    pub fn is_empty(&self) -> bool { self.store.is_empty() }
    pub fn active_count(&self) -> usize { self.store.active_count() }
    pub fn iter(&self) -> impl Iterator<Item = &T> { self.store.iter() }
}
