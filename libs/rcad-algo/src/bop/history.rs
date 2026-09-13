// BRepTools_History semantic equivalent (BRepTools_History.cxx / .hxx) — the
// history of modifications, generations and removals of the boolean operation.

use rcad_kernel::topods::{ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use std::collections::HashMap;
use std::collections::HashSet;

/// OCCT BRepTools_History::IsSupportedType (BRepTools_History.hxx L145-153):
/// only VERTEX, EDGE, FACE and SOLID are supported.
pub fn is_supported_type(s: &Shape) -> bool {
    match s.shape_type() {
        ShapeType::Vertex | ShapeType::Edge | ShapeType::Face | ShapeType::Solid => true,
        _ => false,
    }
}

/// OCCT history read surface of an algorithm — IsDeleted(theS) / Modified(theS)
/// / Generated(theS) — as consumed by the BRepTools_History template
/// constructor and the Merge template (BRepTools_History.hxx L99-132,
/// L212-217).
///
/// OCCT carries it as the `template <class TheAlgo>` type parameter
/// (theAlgo.IsDeleted / theAlgo.Modified / theAlgo.Generated). Rust has no
/// class templates, so the same three methods are carried by this trait and
/// the constructor takes `&dyn HistoryAlgo` — the call sites name exactly the
/// algorithms OCCT instantiates the template with (BOPAlgo_Builder /
/// BOPAlgo_MakerVolume through their BOPAlgo_BuilderShape base,
/// BRepAlgoAPI_BuilderAlgo).
pub trait HistoryAlgo {
    /// OCCT theAlgo.IsDeleted(aS) — BOPAlgo_BuilderShape::IsDeleted
    /// (BOPAlgo_BuilderShape.hxx L72-75).
    fn is_deleted(&self, the_s: &Shape) -> bool;

    /// OCCT theAlgo.Modified(aS) — BOPAlgo_BuilderShape::Modified
    /// (BOPAlgo_BuilderShape.hxx L52-58).
    fn modified(&self, the_s: &Shape) -> Vec<Shape>;

    /// OCCT theAlgo.Generated(aS) — BOPAlgo_BuilderShape::Generated
    /// (BOPAlgo_BuilderShape.hxx L61-67).
    fn generated(&self, the_s: &Shape) -> Vec<Shape>;
}

/// OCCT BRepTools_History::TRelationType (BRepTools_History.hxx L136-141):
/// the types of the historical relations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TRelationType {
    /// OCCT TRelationType_Removed.
    Removed,
    /// OCCT TRelationType_Generated.
    Generated,
    /// OCCT TRelationType_Modified.
    Modified,
}

/// OCCT BRepTools_History (BRepTools_History.cxx) — myShapeToModified /
/// myShapeToGenerated DataMaps and the myRemoved map.
pub struct BRepToolsHistory {
    // myShapeToModified: DataMap<Shape, List<Shape>>
    my_shape_to_modified: HashMap<(u64, u32), Vec<Shape>>,
    // myShapeToGenerated: DataMap<Shape, List<Shape>>
    my_shape_to_generated: HashMap<(u64, u32), Vec<Shape>>,
    // myRemoved: NCollection_Map<Shape> (TShape + Location)
    my_removed: HashSet<(u64, u32)>,
}

impl BRepToolsHistory {
    pub fn new() -> Self {
        Self {
            my_shape_to_modified: HashMap::new(),
            my_shape_to_generated: HashMap::new(),
            my_removed: HashSet::new(),
        }
    }

    fn key(s: &Shape) -> (u64, u32) {
        (s.ptr_id(), s.location)
    }

    /// OCCT BRepTools_History::BRepTools_History(const NCollection_List
    /// <TopoDS_Shape>& theArguments, TheAlgo& theAlgo) (BRepTools_History.hxx
    /// L99-132) — the history of one algorithm over the given arguments.
    ///
    /// Rust has no class template, so the algorithm is taken through the
    /// [`HistoryAlgo`] trait carrying OCCT's IsDeleted / Modified / Generated.
    pub fn from_algorithm(the_arguments: &[Shape], the_algo: &dyn HistoryAlgo) -> Self {
        // OCCT L103-109: map all argument shapes to save them in history —
        // TopExp::MapShapes of every non-null argument into an
        // NCollection_IndexedMap (pre-order depth-first, deduplicated by
        // TShape + Location).
        let mut an_args_map: Vec<Shape> = Vec::new();
        let mut a_seen: HashSet<(u64, u32)> = HashSet::new();
        for a_arg in the_arguments {
            // OCCT L107: if (!aIt.Value().IsNull()) TopExp::MapShapes(...)
            if a_arg.is_null() {
                continue;
            }
            map_shapes(a_arg, &mut an_args_map, &mut a_seen);
        }

        // OCCT L112-131: copy the history for all supported shapes.
        let mut a_history = Self::new();
        for a_s in &an_args_map {
            // OCCT L116-117: if (!IsSupportedType(aS)) continue;
            if !is_supported_type(a_s) {
                continue;
            }

            // OCCT L119-120: if (theAlgo.IsDeleted(aS)) Remove(aS);
            if the_algo.is_deleted(a_s) {
                a_history.remove(a_s);
            }

            // OCCT L122-125: Modified — for (aIt.Initialize(aModified) ...)
            // AddModified(aS, aIt.Value());
            for a_m in the_algo.modified(a_s) {
                a_history.add_modified(a_s, &a_m);
            }

            // OCCT L127-130: Generated — AddGenerated(aS, aIt.Value());
            for a_g in the_algo.generated(a_s) {
                a_history.add_generated(a_s, &a_g);
            }
        }
        a_history
    }

    /// OCCT BRepTools_History::Merge(const NCollection_List<TopoDS_Shape>&
    /// theArguments, TheAlgo& theAlgo) (BRepTools_History.hxx L212-217):
    ///   Merge(BRepTools_History(theArguments, theAlgo));
    pub fn merge_algorithm(&mut self, the_arguments: &[Shape], the_algo: &dyn HistoryAlgo) {
        self.merge(&Self::from_algorithm(the_arguments, the_algo));
    }

    /// OCCT BRepTools_History::AddGenerated (BRepTools_History.cxx L48-67).
    pub fn add_generated(&mut self, the_initial: &Shape, the_generated: &Shape) {
        if !is_supported_type(the_initial) || !is_supported_type(the_generated) {
            return;
        }
        let list = self
            .my_shape_to_generated
            .entry(Self::key(the_initial))
            .or_default();
        if !list.iter().any(|g| Self::key(g) == Self::key(the_generated)) {
            list.push(the_generated.clone());
        }
    }

    /// OCCT BRepTools_History::AddModified (BRepTools_History.cxx L69-88).
    pub fn add_modified(&mut self, the_initial: &Shape, the_modified: &Shape) {
        if !is_supported_type(the_initial) || !is_supported_type(the_modified) {
            return;
        }
        let list = self
            .my_shape_to_modified
            .entry(Self::key(the_initial))
            .or_default();
        if !list.iter().any(|m| Self::key(m) == Self::key(the_modified)) {
            list.push(the_modified.clone());
        }
    }

    /// OCCT BRepTools_History::Remove (BRepTools_History.cxx L91-108) — unbind
    /// the modifications and add the shape to myRemoved.
    pub fn remove(&mut self, the_removed: &Shape) {
        if !is_supported_type(the_removed) {
            return;
        }
        self.my_shape_to_modified.remove(&Self::key(the_removed));
        self.my_removed.insert(Self::key(the_removed));
    }

    /// OCCT BRepTools_History::Modified.
    pub fn modified(&self, the_initial: &Shape) -> Vec<Shape> {
        self.my_shape_to_modified
            .get(&Self::key(the_initial))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT BRepTools_History::Generated.
    pub fn generated(&self, the_initial: &Shape) -> Vec<Shape> {
        self.my_shape_to_generated
            .get(&Self::key(the_initial))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT BRepTools_History::IsRemoved.
    pub fn is_removed(&self, the_initial: &Shape) -> bool {
        self.my_removed.contains(&Self::key(the_initial))
    }

    /// OCCT BRepTools_History::HasGenerated (BRepTools_History.hxx L192).
    pub fn has_generated(&self) -> bool {
        !self.my_shape_to_generated.is_empty()
    }

    /// OCCT BRepTools_History::HasModified (BRepTools_History.hxx L195).
    pub fn has_modified(&self) -> bool {
        !self.my_shape_to_modified.is_empty()
    }

    /// OCCT BRepTools_History::HasRemoved (BRepTools_History.hxx L198).
    pub fn has_removed(&self) -> bool {
        !self.my_removed.is_empty()
    }

    /// OCCT BRepTools_History::Clear (BRepTools_History.hxx L171-177).
    pub fn clear(&mut self) {
        self.my_shape_to_modified.clear();
        self.my_shape_to_generated.clear();
        self.my_removed.clear();
    }

    /// OCCT BRepTools_History::Merge (BRepTools_History.cxx L203-318):
    /// merges the next history (theHistory23) into this one, propagating the
    /// removals, modifications and generations of the 2-3 range through the
    /// 1-2 range so the result covers the 1-3 span.
    pub fn merge(&mut self, the_history23: &BRepToolsHistory) {
        // OCCT L206-213: nothing to merge when the 2-3 history is empty.
        if !(the_history23.has_modified()
            || the_history23.has_generated()
            || the_history23.has_removed())
        {
            return;
        }

        // Propagate R23 directly and M23 and G23 fully to M12 and G12.
        // Remember the propagated shapes.
        let mut a_r_propagated: HashSet<(u64, u32)> = HashSet::new();
        // OCCT L214: remember the propagated shapes.
        let mut a_m_and_g_propagated: HashSet<(u64, u32)> = HashSet::new();
        // OCCT L210-211: aS1ToGAndM = {&myShapeToGenerated,
        // &myShapeToModified} — aI == 0 walks the generated map, aI == 1 the
        // modified map.
        for a_i in 0..2 {
            let keys: Vec<(u64, u32)> = if a_i == 0 {
                self.my_shape_to_generated.keys().copied().collect()
            } else {
                self.my_shape_to_modified.keys().copied().collect()
            };
            // OCCT L242-252 (the aI != 0 G-additions cross-write into
            // myShapeToGenerated) — collected during the walk and applied
            // afterwards (the walk reads only its own map).
            let mut a_gen_cross_writes: Vec<((u64, u32), Vec<Shape>)> = Vec::new();
            for a_s1_key in keys {
                // The list is taken out of the map so it can be rebuilt in
                // place (OCCT mutates aL12 through ChangeValue + Remove).
                let mut a_l12 = if a_i == 0 {
                    self.my_shape_to_generated.remove(&a_s1_key).unwrap_or_default()
                } else {
                    self.my_shape_to_modified.remove(&a_s1_key).unwrap_or_default()
                };
                // OCCT L219: the G and M additions.
                let mut a_additions: [Vec<Shape>; 2] = [Vec::new(), Vec::new()];
                let mut a_kept: Vec<Shape> = Vec::with_capacity(a_l12.len());
                for a_s2 in a_l12.drain(..) {
                    if the_history23.is_removed(&a_s2) {
                        // OCCT L222-225: R23 propagates directly and aS2 is
                        // dropped from the list.
                        a_r_propagated.insert(Self::key(&a_s2));
                    } else {
                        // OCCT L226-231: the G23 addition (aS2 stays in the
                        // list unless the M branch below removes it).
                        if let Some(a_g) = the_history23
                            .my_shape_to_generated
                            .get(&Self::key(&a_s2))
                        {
                            a_additions[0].extend(a_g.iter().cloned());
                            a_m_and_g_propagated.insert(Self::key(&a_s2));
                        }
                        // OCCT L232-240: the M23 addition; aS2 is replaced by
                        // its modifications when present, otherwise kept.
                        if let Some(a_m) = the_history23
                            .my_shape_to_modified
                            .get(&Self::key(&a_s2))
                        {
                            a_additions[a_i].extend(a_m.iter().cloned());
                            a_m_and_g_propagated.insert(Self::key(&a_s2));
                        } else {
                            a_kept.push(a_s2);
                        }
                    }
                }
                // OCCT L241: add(aL12, aAdditions[aI]).
                a_kept.extend(a_additions[a_i].drain(..));
                match a_i {
                    0 => {
                        self.my_shape_to_generated.insert(a_s1_key, a_kept);
                    }
                    _ => {
                        self.my_shape_to_modified.insert(a_s1_key, a_kept);
                        // OCCT L242-252: when a modified list gained G
                        // additions, they propagate into the generated entry
                        // of aS1 as well.
                        if !a_additions[0].is_empty() {
                            a_gen_cross_writes.push((a_s1_key, std::mem::take(&mut a_additions[0])));
                        }
                    }
                }
            }
            for (a_s1_key, a_adds) in a_gen_cross_writes {
                let entry = self.my_shape_to_generated.entry(a_s1_key).or_default();
                entry.extend(a_adds);
            }
        }

        // OCCT L256-291: propagate M23 and G23 to M12 and G12 sequentially.
        for a_i in 0..2 {
            let keys: Vec<(u64, u32)> = if a_i == 0 {
                the_history23.my_shape_to_generated.keys().copied().collect()
            } else {
                the_history23.my_shape_to_modified.keys().copied().collect()
            };
            for a_s2_key in keys {
                if !a_m_and_g_propagated.contains(&a_s2_key) {
                    let a_m2 = if a_i == 0 {
                        the_history23
                            .my_shape_to_generated
                            .get(&a_s2_key)
                            .cloned()
                            .unwrap_or_default()
                    } else {
                        the_history23
                            .my_shape_to_modified
                            .get(&a_s2_key)
                            .cloned()
                            .unwrap_or_default()
                    };
                    let entry = if a_i == 0 {
                        self.my_shape_to_generated.entry(a_s2_key).or_default()
                    } else {
                        self.my_shape_to_modified.entry(a_s2_key).or_default()
                    };
                    entry.extend(a_m2);
                    // OCCT L286: myRemoved.Remove(aS2).
                    self.my_removed.remove(&a_s2_key);
                }
            }
        }

        // OCCT L294-309: unbind the empty M12 and G12 (an emptied entry
        // becomes a removal).
        for a_i in 0..2 {
            let empty_keys: Vec<(u64, u32)> = if a_i == 0 {
                self.my_shape_to_generated
                    .iter()
                    .filter(|(_, l)| l.is_empty())
                    .map(|(k, _)| *k)
                    .collect()
            } else {
                self.my_shape_to_modified
                    .iter()
                    .filter(|(_, l)| l.is_empty())
                    .map(|(k, _)| *k)
                    .collect()
            };
            for a_s1_key in empty_keys {
                // OCCT L304: myRemoved.Add(aS1).
                self.my_removed.insert(a_s1_key);
                if a_i == 0 {
                    self.my_shape_to_generated.remove(&a_s1_key);
                } else {
                    self.my_shape_to_modified.remove(&a_s1_key);
                }
            }
        }

        // OCCT L312-318: propagate R23 to R12 sequentially.
        for a_s2_key in &the_history23.my_removed {
            if !a_r_propagated.contains(a_s2_key)
                && !self.my_shape_to_modified.contains_key(a_s2_key)
                && !self.my_shape_to_generated.contains_key(a_s2_key)
            {
                self.my_removed.insert(*a_s2_key);
            }
        }
    }

    /// Sanity: TShape referenced for dead-code elimination.
    #[allow(dead_code)]
    fn _shape_type(s: &Shape) -> Option<&'static str> {
        match &*s.data {
            TShape::Vertex(_) => Some("vertex"),
            TShape::Edge(_) => Some("edge"),
            TShape::Face(_) => Some("face"),
            TShape::Solid(_) => Some("solid"),
            _ => None,
        }
    }
}

/// OCCT TopExp::MapShapes(const TopoDS_Shape& S, NCollection_IndexedMap<...>&
/// M) (TopExp.cxx L49-61):
///   M.Add(S);
///   TopoDS_Iterator It(S, cumOri, cumLoc);
///   for (; It.More(); It.Next()) MapShapes(It.Value(), M);
/// A pre-order depth-first walk; the IndexedMap deduplicates by TShape +
/// Location (TopTools_ShapeMapHasher), and an already-mapped shape is still
/// traversed (the Add result is not tested), so `a_seen` gates only the push.
///
/// The direct sub-shapes come from the shared Builder walk
/// (Builder::shape_sub_shapes_static), which mirrors TopoDS_Iterator of the
/// rcad TShape: an edge's vertices carry the composed edge Location, the
/// other children carry the Location indices of the rcad data model.
fn map_shapes(the_shape: &Shape, the_map: &mut Vec<Shape>, a_seen: &mut HashSet<(u64, u32)>) {
    if a_seen.insert((the_shape.ptr_id(), the_shape.location)) {
        the_map.push(the_shape.clone());
    }
    for a_sub in crate::bop::algo::builder::Builder::shape_sub_shapes_static(the_shape) {
        map_shapes(&a_sub, the_map, a_seen);
    }
}
