// 1:1 semantic equivalents of OCCT NCollection_Map<int> and
// NCollection_DataMap<K, V> (NCollection_Map.hxx + NCollection_DataMap.hxx +
// NCollection_BaseMap.cxx + NCollection_Primes.cxx +
// NCollection_DefaultHasher.hxx).
//
// The iteration order of these containers is deterministic (bucket array
// order 0..myNbBuckets with head-insert chains), unlike
// std::collections::HashMap whose RandomState seed differs per process. Code
// that iterates them must reproduce this exact order, so a plain HashMap is
// not a valid replacement. This module models the OCCT containers: same bucket
// count growth (NextPrimeForMap), same hash (key % NbBuckets + 1), same
// head-insert chains and same rehash (old buckets 0..Nb, chain head-to-tail,
// pushed to the head of the new bucket).
//
// Key hashing matches OCCT NCollection_DefaultHasher:
//  - int / unsigned int keys: static_cast<size_t> (identity), so
//    key % NbBuckets + 1;
//  - opencascade::handle<T> keys: std::hash<handle<T>> =
//    static_cast<size_t>(reinterpret_cast<uintptr_t>(get())) (Standard_Handle.hxx
//    L455-461) — identity of the raw pointer, which rcad stores as u64
//    (Arc pointer id). For handle keys the OCCT hash is the pointer value
//    modulo the bucket count.

// OCCT NCollection_Primes::THE_PRIME_VECTOR (NCollection_Primes.cxx).
const THE_PRIME_VECTOR: [usize; 24] = [
    101, 1009, 2003, 5003, 10007, 20011, 37003, 57037, 65003, 100019, 209953, 472393, 995329,
    2359297, 4478977, 9437185, 17915905, 35831809, 71663617, 150994945, 301989889, 573308929,
    1019215873, 2038431745,
];

/// OCCT NCollection_Primes::NextPrimeForMap (NCollection_Primes.cxx): the first
/// prime >= theN + 1, or theN + 1 when it exceeds the largest available prime.
fn next_prime_for_map(n: usize) -> usize {
    match THE_PRIME_VECTOR.iter().find(|&&p| p >= n + 1) {
        Some(&p) => p,
        None => n + 1,
    }
}

/// Keys usable in the OCCT map models. `hash_code` returns the raw hash value
/// (before the `% NbBuckets + 1` reduction of OCCT NCollection_Map::HashCode,
/// NCollection_Map.hxx L593-596).
pub trait OcctHashKey: Copy + PartialEq {
    fn hash_code(self) -> usize;
}

impl OcctHashKey for usize {
    fn hash_code(self) -> usize {
        self
    }
}

impl OcctHashKey for u64 {
    fn hash_code(self) -> usize {
        self as usize
    }
}

/// OCCT NCollection_Map::HashCode (NCollection_Map.hxx L593-596):
/// myHasher(theKey) % theUpperBound + 1.
#[inline]
fn occt_hash<K: OcctHashKey>(key: K, upper: usize) -> usize {
    key.hash_code() % upper + 1
}

/// OCCT NCollection_DataMap<K, V> semantic equivalent.
///
/// Each bucket holds a chain of (key, value) pairs; index 0 is the chain head
/// (the most recently inserted element, as OCCT inserts at the head).
pub struct OcctDataMapInt<K, V> {
    // Element count (mySize).
    my_size: usize,
    // Bucket count (myNbBuckets). OCCT default constructor sets it to 1 with
    // a null array, allocated at the first insert.
    my_nb_buckets: usize,
    // Buckets 0..=myNbBuckets, each a chain; index 0 is the chain head.
    my_data: Vec<Vec<(K, V)>>,
}

impl<K: OcctHashKey, V: Default + Clone> OcctDataMapInt<K, V> {
    /// OCCT NCollection_Map() / NCollection_DataMap() empty constructors:
    /// NCollection_BaseMap(1, ...).
    pub fn new() -> Self {
        Self {
            my_size: 0,
            my_nb_buckets: 1,
            my_data: Vec::new(),
        }
    }

    /// OCCT NCollection_BaseMap::Resizable (NCollection_BaseMap.hxx L217):
    /// IsEmpty() || (mySize > myNbBuckets).
    fn resizable(&self) -> bool {
        self.my_size == 0 || self.my_size > self.my_nb_buckets
    }

    /// OCCT NCollection_Map::ReSize (NCollection_Map.hxx L283-312) +
    /// NCollection_BaseMap::BeginResize/EndResize (NCollection_BaseMap.cxx):
    /// rehash into NextPrimeForMap(theExtent) buckets. Old buckets 0..=Nb are
    /// walked in order, each chain head-to-tail, and every node is pushed to
    /// the head of its new bucket.
    fn resize(&mut self, extent: usize) {
        let mut new_buckets = next_prime_for_map(extent);
        if new_buckets <= self.my_nb_buckets {
            if self.my_data.is_empty() {
                // BeginResize: (!myData1) -> theNewBuckets = myNbBuckets.
                new_buckets = self.my_nb_buckets;
            } else {
                // BeginResize returns false: no rehash.
                return;
            }
        }
        let mut new_data = vec![Vec::new(); new_buckets + 1];
        for bucket in &self.my_data {
            for entry in bucket {
                let k = occt_hash(entry.0, new_buckets);
                new_data[k].insert(0, entry.clone());
            }
        }
        self.my_data = new_data;
        self.my_nb_buckets = new_buckets;
    }

    /// OCCT NCollection_DataMap::Bound (NCollection_DataMap.hxx): returns a
    /// reference to the value for theKey, inserting a default-constructed one
    /// when absent.
    pub fn bound(&mut self, key: K) -> &mut V {
        if self.resizable() {
            self.resize(self.my_size);
        }
        let k = occt_hash(key, self.my_nb_buckets);
        if let Some(chain) = self.my_data.get_mut(k) {
            if let Some(pos) = chain.iter().position(|(k2, _)| *k2 == key) {
                return &mut chain[pos].1;
            }
            self.my_size += 1;
            chain.insert(0, (key, V::default()));
            return &mut chain[0].1;
        }
        // Unreachable after the first resize (my_data has myNbBuckets + 1
        // buckets); the empty map resized above.
        unreachable!("OcctDataMapInt: bucket array not allocated")
    }

    /// OCCT NCollection_DataMap::Bind/Add — insert or replace the value.
    /// Returns the previous value when the key existed (HashMap::insert
    /// semantics used by the rcad call sites).
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if self.resizable() {
            self.resize(self.my_size);
        }
        let k = occt_hash(key, self.my_nb_buckets);
        if let Some(chain) = self.my_data.get_mut(k) {
            if let Some(pos) = chain.iter().position(|(k2, _)| *k2 == key) {
                let old = std::mem::replace(&mut chain[pos].1, value);
                return Some(old);
            }
            self.my_size += 1;
            chain.insert(0, (key, value));
            return None;
        }
        unreachable!("OcctDataMapInt: bucket array not allocated")
    }

    /// OCCT NCollection_DataMap::ChangeSeek — the value when present.
    pub fn get_mut(&mut self, key: K) -> Option<&mut V> {
        if self.my_data.is_empty() {
            return None;
        }
        let k = occt_hash(key, self.my_nb_buckets);
        self.my_data[k]
            .iter_mut()
            .find(|(k2, _)| *k2 == key)
            .map(|(_, v)| v)
    }

    /// OCCT NCollection_DataMap::Seek — the value when present.
    pub fn get(&self, key: K) -> Option<&V> {
        if self.my_data.is_empty() {
            return None;
        }
        let k = occt_hash(key, self.my_nb_buckets);
        self.my_data[k]
            .iter()
            .find(|(k2, _)| *k2 == key)
            .map(|(_, v)| v)
    }

    /// OCCT NCollection_DataMap::Contains.
    pub fn contains(&self, key: K) -> bool {
        self.get(key).is_some()
    }

    /// OCCT NCollection_DataMap::Remove — remove the key, returning the value
    /// when present.
    pub fn remove(&mut self, key: K) -> Option<V> {
        if self.my_size == 0 {
            return None;
        }
        let k = occt_hash(key, self.my_nb_buckets);
        if let Some(pos) = self.my_data[k].iter().position(|(k2, _)| *k2 == key) {
            let (_, v) = self.my_data[k].remove(pos);
            self.my_size -= 1;
            return Some(v);
        }
        None
    }

    /// OCCT NCollection_DataMap::Iterator (NCollection_BaseMap.hxx L47-135):
    /// bucket 0..=myNbBuckets, chain head-to-tail.
    pub fn iter(&self) -> impl Iterator<Item = (K, &V)> + '_ {
        self.my_data.iter().flatten().map(|(k, v)| (*k, v))
    }

    /// Keys in the OCCT iteration order.
    pub fn iter_keys(&self) -> impl Iterator<Item = K> + '_ {
        self.my_data.iter().flatten().map(|(k, _)| *k)
    }

    /// OCCT NCollection_DataMap::Clear (release the data; the bucket array is
    /// kept for reuse when doReleaseMemory is false).
    pub fn clear(&mut self) {
        self.my_size = 0;
        for chain in &mut self.my_data {
            chain.clear();
        }
    }

    /// OCCT NCollection_DataMap::Extent (element count).
    pub fn len(&self) -> usize {
        self.my_size
    }

    /// OCCT NCollection_DataMap::IsEmpty.
    pub fn is_empty(&self) -> bool {
        self.my_size == 0
    }
}

/// OCCT NCollection_Map<int> semantic equivalent (NCollection_Map.hxx) —
/// the same bucket machinery as NCollection_DataMap with no values.
pub type OcctMapInt = OcctDataMapInt<usize, ()>;

impl OcctMapInt {
    /// OCCT NCollection_Map::Add (NCollection_Map.hxx L321, addImpl): resize
    /// when Resizable(), then insert at the head of the bucket when the key is
    /// absent. Returns true when the key was newly added.
    pub fn add(&mut self, key: usize) -> bool {
        if self.my_data.is_empty() || self.resizable() {
            self.resize(self.my_size);
        }
        let k = occt_hash(key, self.my_nb_buckets);
        if let Some(chain) = self.my_data.get_mut(k) {
            if chain.iter().any(|(k2, _)| *k2 == key) {
                return false;
            }
            self.my_size += 1;
            chain.insert(0, (key, ()));
            return true;
        }
        unreachable!("OcctMapInt: bucket array not allocated")
    }
}
