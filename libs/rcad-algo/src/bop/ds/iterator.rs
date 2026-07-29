// OCCT BOPDS_Iterator — BVH-based pair enumeration.
// OCCT ref: BOPDS_Iterator.hxx / BOPDS_Iterator.cxx
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use rcad_kernel::topods::ShapeType;
use std::collections::HashSet;

/// BOPDS_Iterator — iterates pairs of interfering shapes.
/// OCCT ref: BOPDS_Iterator.hxx
pub struct BOPDS_Iterator<'a> {
    ds: &'a DS,
    fuzzy_tol: f64,
    my_lists: Vec<Vec<(usize, usize)>>,
    my_length: usize, // OCCT: myLength — expected number of interfering pairs
    my_run_parallel: bool, // OCCT: myRunParallel
    current_list: Vec<(usize, usize)>,
    current_pos: usize,
}

/// Maps (ShapeType, ShapeType) to list index (OCCT: BOPDS_Tools::TypeToInteger).
/// Handles both orderings: (Vertex, Edge) and (Edge, Vertex) both map to VE (1).
fn type_to_integer(t1: ShapeType, t2: ShapeType) -> usize {
    match (t1, t2) {
        (ShapeType::Vertex, ShapeType::Vertex) => 0, // VV
        (ShapeType::Vertex, ShapeType::Edge)
        | (ShapeType::Edge, ShapeType::Vertex) => 1, // VE + EV
        (ShapeType::Edge, ShapeType::Edge) => 2,     // EE
        (ShapeType::Vertex, ShapeType::Face)
        | (ShapeType::Face, ShapeType::Vertex) => 3, // VF + FV
        (ShapeType::Edge, ShapeType::Face)
        | (ShapeType::Face, ShapeType::Edge) => 4,   // EF + FE
        (ShapeType::Face, ShapeType::Face) => 5,     // FF
        _ => 0,
    }
}

/// Returns true if bounding boxes of two shapes overlap (with gap tolerance).
fn boxes_overlap(si1: &crate::bop::ds::ShapeInfo, si2: &crate::bop::ds::ShapeInfo, gap: f64) -> bool {
    use rcad_kernel::math::bnd::BndBox;
    let b1 = match (si1.box_min, si1.box_max) {
        (Some(min), Some(max)) => {
            let mut b = BndBox::from_corners(min.x, min.y, min.z, max.x, max.y, max.z);
            b.set_gap(si1.box_gap + gap);
            b
        }
        _ => return true,
    };
    let b2 = match (si2.box_min, si2.box_max) {
        (Some(min), Some(max)) => {
            let mut b = BndBox::from_corners(min.x, min.y, min.z, max.x, max.y, max.z);
            b.set_gap(si2.box_gap + gap);
            b
        }
        _ => return true,
    };
    !b1.is_out_box(&b2)
}

impl<'a> BOPDS_Iterator<'a> {
    pub fn new(ds: &'a DS, fuzzy_tol: f64) -> Self {
        let n = 10;
        let mut my_lists = Vec::with_capacity(n);
        for _ in 0..n { my_lists.push(Vec::new()); }
        BOPDS_Iterator {
            ds, fuzzy_tol,
            my_lists,
            my_length: 0,
            my_run_parallel: false,
            current_list: Vec::new(),
            current_pos: 0,
        }
    }

    // OCCT BOPDS_Iterator.cxx L167-188: ExpectedLength, BlockLength
    pub fn expected_length(&self) -> usize { self.my_length }
    pub fn block_length(&self) -> usize {
        let a_nb_iis = self.my_length;
        let a_cf_predict = 0.5;
        if a_nb_iis <= 1 { return 1; }
        (a_cf_predict * a_nb_iis as f64) as usize
    }

    // OCCT: SetRunParallel / RunParallel
    pub fn set_run_parallel(&mut self, flag: bool) { self.my_run_parallel = flag; }
    pub fn run_parallel(&self) -> bool { self.my_run_parallel }

    // OCCT BOPDS_Iterator.hxx L107: NbExtInterfs
    pub fn nb_ext_interfs() -> usize { 4 }

    // OCCT BOPDS_Iterator::Prepare (BOPDS_Iterator.cxx L247-265).
    // Calls Intersect to populate pair lists by checking AABB overlap.
    pub fn prepare(&mut self, _ctx: Option<&IntToolsContext>, _check_obb: bool, _fuzzy: f64) {
        let a_nb_interf_types = DS::nb_interf_types();
        self.my_length = 0;
        for i in 0..a_nb_interf_types {
            self.my_lists[i].clear();
        }
        if self.ds.nb_shapes() == 0 { return; }
        self.intersect();
    }

    // OCCT BOPDS_Iterator::Intersect (BOPDS_Iterator.cxx L270-359).
    // OCCT: BVH-based pair selection + range-based cross-argument pairing.
    // rcad: brute-force AABB + same range-based filtering (equivalent for correctness).
    fn intersect(&mut self) {
        let a_nb = self.ds.nb_source_shapes();
        // OCCT L279-284: Add shapes with BRep (VERTEX|EDGE|FACE) to BVH
        let is_brep = |st: ShapeType| -> bool {
            matches!(st, ShapeType::Vertex | ShapeType::Edge | ShapeType::Face)
        };

        let a_nb_r = self.ds.nb_ranges();
        // OCCT L306-358: for each range, pair with shapes from OTHER ranges
        for i_r in 0..a_nb_r {
            let a_range = self.ds.range(i_r).clone();
            for i in a_range.first..=a_range.last {
                if i >= a_nb { break; }
                let si1 = &self.ds.shapes[i];
                if !is_brep(si1.shape_type) { continue; } // OCCT L282: HasBRep
                let a_type1 = si1.shape_type;

                for j in (i + 1)..a_nb {
                    if a_range.contains(j) { continue; } // OCCT L320-324: skip same-range
                    let si2 = &self.ds.shapes[j];
                    if !is_brep(si2.shape_type) { continue; } // OCCT L282: HasBRep
                    let a_type2 = si2.shape_type;

                    // OCCT L336-340: avoid interfering shape with its sub-shapes
                    // OCCT uses TypeToInteger(aType) ordering: VERTEX(7)>EDGE(6)>FACE(4)
                    // rcad ShapeType discriminant: Vertex(0)<Edge(1)<Face(3) — same relative order
                    if (a_type1 as isize > a_type2 as isize && si1.has_sub_shape(j))
                        || (a_type2 as isize > a_type1 as isize && si2.has_sub_shape(i))
                    {
                        continue;
                    }

                    if boxes_overlap(si1, si2, self.fuzzy_tol) {
                        let i_x = type_to_integer(a_type1, a_type2);
                        self.my_lists[i_x].push((i.min(j), i.max(j)));
                    }
                }
            }
        }
    }

    /// OCCT BOPDS_Iterator::IntersectExt (BOPDS_Iterator.cxx L182-230).
    /// Expands the pair lists to include all pairs involving shapes from `theMap`.
    /// Used by RepeatIntersection to add pairs with increased-tolerance vertices.
    pub fn intersect_ext(&mut self, the_map: &HashSet<usize>) {
        let n = self.ds.nb_shapes();
        let t_combos = [
            (ShapeType::Vertex, ShapeType::Vertex),
            (ShapeType::Vertex, ShapeType::Edge),
            (ShapeType::Vertex, ShapeType::Face),
        ];
        for &(t1, t2) in &t_combos {
            let idx = type_to_integer(t1, t2);
            let list = &mut self.my_lists[idx];
            for &mi in the_map {
                let si1 = &self.ds.shapes[mi];
                if si1.shape_type != ShapeType::Vertex { continue; }
                for j in 0..n {
                    if j == mi || self.ds.shapes[j].shape_type != t2 { continue; }
                    let si2 = &self.ds.shapes[j];
                    if boxes_overlap(si1, si2, self.fuzzy_tol) {
                        let pair = if mi < j { (mi, j) } else { (j, mi) };
                        if !list.contains(&pair) {
                            list.push(pair);
                        }
                    }
                }
            }
        }
    }

    // OCCT BOPDS_Iterator::Initialize (BOPDS_Iterator.cxx L192-208).
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        self.my_length = 0;
        let i_x = type_to_integer(t1, t2);
        if i_x < self.my_lists.len() {
            // OCCT L202: stable_sort for deterministic iteration order
            self.my_lists[i_x].sort();
            self.current_list = self.my_lists[i_x].clone();
            self.my_length = self.current_list.len();
        } else {
            self.current_list = Vec::new();
        }
        self.current_pos = 0;
    }
    pub fn more(&self) -> bool { self.current_pos < self.current_list.len() }
    pub fn next(&mut self) { self.current_pos += 1; }
    pub fn value(&self) -> (usize, usize) { self.current_list[self.current_pos] }
    pub fn pairs(&self, t1: ShapeType, t2: ShapeType) -> &[(usize, usize)] {
        let idx = type_to_integer(t1, t2);
        if idx < self.my_lists.len() { &self.my_lists[idx] } else { &[] }
    }
}