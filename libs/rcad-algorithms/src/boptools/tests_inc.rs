#[cfg(test)]
mod tests {
    use crate::boptools::*;

    #[test]
    fn test_empty_set() {
        let s = BOPToolsSet::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
        assert_eq!(s.sum, 0);
    }

    #[test]
    fn test_add_single() {
        let mut s = BOPToolsSet::new();
        s.add(5);
        assert!(!s.is_empty());
        assert_eq!(s.len(), 1);
        assert_eq!(s.faces(), &[5]);
    }

    #[test]
    fn test_add_sorted_dedup() {
        let mut s = BOPToolsSet::new();
        s.add(3); s.add(1); s.add(2); s.add(1);
        assert_eq!(s.len(), 3);
        assert_eq!(s.faces(), &[1, 2, 3]);
    }

    #[test]
    fn test_equality() {
        let mut a = BOPToolsSet::new();
        a.add(1); a.add(2); a.add(3);
        let mut b = BOPToolsSet::new();
        b.add(3); b.add(2); b.add(1);
        assert_eq!(a, b);
        b.add(4);
        assert_ne!(a, b);
    }

    #[test]
    fn test_hash_set_dedup() {
        use std::collections::HashSet;
        let mut a = BOPToolsSet::new();
        a.add(1); a.add(2);
        let mut b = BOPToolsSet::new();
        b.add(2); b.add(1);
        let mut c = BOPToolsSet::new();
        c.add(1); c.add(3);

        let mut set = HashSet::new();
        assert!(set.insert(a.clone()));
        // Same content 鈫?no insert (duplicate)
        assert!(!set.insert(b));
        // Different content 鈫?insert
        assert!(set.insert(c));
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_from_slice() {
        let s = BOPToolsSet::from(&[2, 1, 3, 1][..]);
        assert_eq!(s.nb_shapes(), 3);
        assert_eq!(s.faces(), &[1, 2, 3]);
    }
}
