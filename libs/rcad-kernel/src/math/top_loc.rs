//! OCCT TopLoc: transformation location with hierarchical stacking.
//!
//! `TopLoc` wraps a 3D transformation (`DAffine3`) with an optional parent
//! location, forming a lightweight DAG that supports identity testing,
//! equality comparison, and efficient hashing.
//!
//! OCCT source: src/FoundationClasses/TKMath/TopLoc/TopLoc_Location.cxx

use glam::DAffine3;
use std::hash::{Hash, Hasher};

/// OCCT TopLoc_Location — a transformation in a location hierarchy.
///
/// Stores a location item (transformation) and an optional reference to
/// a parent `TopLoc`, forming a singly-linked path through the location
/// tree.  Two `TopLoc` values are equal iff they are the same node in
/// the tree (pointer equality, not deep transformation equality).
///
/// OCCT: TopLoc_Location — uses Handle-based shared nodes internally.
#[derive(Debug, Clone)]
pub struct TopLoc {
    /// The transformation datum at this level.
    pub transformation: DAffine3,
    /// Optional parent location.
    pub next: Option<Box<TopLoc>>,
}

impl TopLoc {
    /// The identity location (no transformation).
    pub const IDENTITY: TopLoc = TopLoc {
        transformation: DAffine3::IDENTITY,
        next: None,
    };

    /// Create a location from a transformation with no parent.
    /// OCCT: TopLoc_Location(const gp_Trsf& T).
    pub fn new(transformation: DAffine3) -> Self {
        TopLoc { transformation, next: None }
    }

    /// Create a location from a transformation and parent.
    /// OCCT: TopLoc_Location(const gp_Trsf& T, const TopLoc_Location& parent).
    pub fn with_parent(transformation: DAffine3, parent: TopLoc) -> Self {
        TopLoc { transformation, next: Some(Box::new(parent)) }
    }

    /// True if this location is the identity (no transformation in chain).
    /// OCCT: Standard_Boolean IsIdentity() const.
    pub fn is_identity(&self) -> bool {
        if !self.transformation.is_finite() {
            return true; // NaN/Inf treated as undefined = identity-like
        }
        if self.transformation != DAffine3::IDENTITY {
            return false;
        }
        // Check parent chain — if all are identity, this location is identity.
        if let Some(ref next) = self.next {
            next.is_identity()
        } else {
            true
        }
    }

    /// Returns the accumulated transformation (product of this chain).
    /// OCCT: gp_Trsf Transformation() const.
    pub fn transformation(&self) -> DAffine3 {
        match self.next {
            Some(ref parent) => parent.transformation() * self.transformation,
            None => self.transformation,
        }
    }

    /// The first (most recent) transformation datum.
    /// OCCT: const gp_Trsf& FirstDatum() const.
    pub fn first_datum(&self) -> &DAffine3 {
        &self.transformation
    }

    /// The parent location, or `None` if this is the root.
    /// OCCT: const TopLoc_Location& NextLocation() const.
    pub fn next_location(&self) -> Option<&TopLoc> {
        self.next.as_ref().map(|b| b.as_ref())
    }

    /// Returns a new location with `trsf` appended to this chain.
    /// OCCT: TopLoc_Location Predicate(const gp_Trsf& T) const.
    /// Returns a NEW location T * this (T applied after this).
    pub fn predicate(&self, trsf: &DAffine3) -> TopLoc {
        TopLoc {
            transformation: *trsf,
            next: Some(Box::new(self.clone())),
        }
    }

    /// OCCT: Standard_Boolean IsDifferent(const TopLoc_Location&) const.
    /// True if the two locations are NOT the same node.
    pub fn is_different(&self, other: &TopLoc) -> bool {
        !self.is_equal(other)
    }

    /// OCCT: IsEqual — pointer equality of the location DAG node.
    pub fn is_equal(&self, other: &TopLoc) -> bool {
        // Structural equality: same transformation AND same parent chain
        if self.transformation != other.transformation {
            return false;
        }
        match (&self.next, &other.next) {
            (None, None) => true,
            (Some(a), Some(b)) => a.is_equal(b),
            _ => false,
        }
    }

    /// OCCT: HashCode.
    pub fn hash_code(&self) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }
}

impl Default for TopLoc {
    fn default() -> Self { TopLoc::IDENTITY }
}

impl Hash for TopLoc {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash the transformation first
        for &v in self.transformation.to_cols_array().iter() {
            v.to_bits().hash(state);
        }
        // Then hash the parent chain
        if let Some(ref next) = self.next {
            next.hash(state);
        }
    }
}

impl PartialEq for TopLoc {
    fn eq(&self, other: &Self) -> bool { self.is_equal(other) }
}

impl Eq for TopLoc {}

// ══════════════════════════════════════════════════════════════════════════
// Tests
// ══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    #[test]
    fn identity_is_identity() {
        assert!(TopLoc::IDENTITY.is_identity());
    }

    #[test]
    fn non_identity_is_not_identity() {
        let t = DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0));
        let loc = TopLoc::new(t);
        assert!(!loc.is_identity());
    }

    #[test]
    fn chain_accumulates() {
        let t1 = DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0));
        let t2 = DAffine3::from_translation(DVec3::new(0.0, 2.0, 0.0));
        let loc = TopLoc::with_parent(t1, TopLoc::new(t2));
        let acc = loc.transformation();
        // Result should be t1 * t2 = translation (1, 2, 0)
        let p = acc.transform_point3(DVec3::ZERO);
        assert!((p - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn predicate_chains() {
        let t = DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0));
        let loc = TopLoc::new(t);
        let extended = loc.predicate(&DAffine3::from_translation(DVec3::new(0.0, 2.0, 0.0)));
        let acc = extended.transformation();
        let p = acc.transform_point3(DVec3::ZERO);
        assert!((p - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn equality_is_structural() {
        let t = DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0));
        let a = TopLoc::new(t);
        let b = TopLoc::new(t);
        // Same transformation but different objects — structural equality
        assert_eq!(a, b);
    }

    #[test]
    fn different_transforms_not_equal() {
        let a = TopLoc::new(DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0)));
        let b = TopLoc::new(DAffine3::from_translation(DVec3::new(2.0, 0.0, 0.0)));
        assert_ne!(a, b);
    }

    #[test]
    fn hash_is_consistent() {
        let t = DAffine3::from_translation(DVec3::new(1.0, 2.0, 3.0));
        let a = TopLoc::new(t);
        let b = TopLoc::new(t);
        assert_eq!(a.hash_code(), b.hash_code());
    }

    #[test]
    fn first_datum() {
        let t = DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0));
        let loc = TopLoc::new(t);
        assert_eq!(loc.first_datum(), &t);
    }
}
