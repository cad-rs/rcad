use crate::bvh::{Aabb, DsBvh};
use super::DS;
use rcad_kernel::topods::ShapeType;

/// OCCT-aligned: BOPDS_Iterator — BVH-based pair enumeration with type bucketing.
///
/// Builds a single BVH tree over all DS sub-shapes (vertices, edges, faces),
/// finds overlapping AABB pairs, buckets them by (type1, type2) combination,
/// applies stable_sort within each bucket, and provides iteration via
/// `Initialize(T1, T2) → More/Next/Value`.
///
/// OCCT BOPDS_Iterator.hxx / .cxx
pub struct BOPDS_Iterator<'a> {
    ds: &'a DS,
    // Per-type-combo pair buckets, indexed by TypeToInteger(t1, t2) result:
    //   0=VV, 1=VE, 2=EE, 3=VF, 4=EF, 5=FF, 6=VZ, 7=EZ, 8=FZ, 9=ZZ
    my_lists: Vec<Vec<(usize, usize)>>,
    // Current iteration state
    current_list: Vec<(usize, usize)>,  // pairs being iterated (cloned from bucket)
    current_pos: usize,                // index into current_list
    my_run_parallel: bool,
}

impl<'a> BOPDS_Iterator<'a> {
    pub fn new(ds: &'a DS) -> Self {
        let n = 10; // NbInterfTypes = 10 (VV..ZZ)
        let mut my_lists = Vec::with_capacity(n);
        for _ in 0..n {
            my_lists.push(Vec::new());
        }
        BOPDS_Iterator {
            ds,
            my_lists,
            current_list: Vec::new(),
            current_pos: 0,
            my_run_parallel: false,
        }
    }

    pub fn set_run_parallel(&mut self, flag: bool) {
        self.my_run_parallel = flag;
    }

    pub fn run_parallel(&self) -> bool {
        self.my_run_parallel
    }

    /// OCCT BOPDS_Tools::TypeToInteger(ShapeType)
    fn type_to_int(t: ShapeType) -> i32 {
        match t {
            ShapeType::Vertex => 7,
            ShapeType::Edge => 6,
            ShapeType::Face => 4,
            ShapeType::Shell => 3,
            ShapeType::Solid => 2,
            ShapeType::Wire => 5,
            ShapeType::Compound => 0,
            ShapeType::CompSolid => 1,
            _ => 9,
        }
    }

    /// OCCT BOPDS_Tools::TypeToInteger(t1, t2) → bucket index
    fn type_to_bucket(t1: ShapeType, t2: ShapeType) -> i32 {
        let i1 = Self::type_to_int(t1);
        let i2 = Self::type_to_int(t2);
        let ix = i2 * 10 + i1;
        match ix {
            77 => 0,  // VV
            76 | 67 => 1,  // VE
            66 => 2,  // EE
            74 | 47 => 3,  // VF
            64 | 46 => 4,  // EF
            44 => 5,  // FF
            72 | 27 => 6,  // VZ
            62 | 26 => 7,  // EZ
            42 | 24 => 8,  // FZ
            22 => 9,  // ZZ
            _ => -1,
        }
    }

    /// OCCT BOPDS_Iterator::Prepare — build BVH, find all overlapping pairs, bucket by type.
    ///
    /// Builds a single BVH over all shapes (vertices + edges + faces), runs candidate_pairs,
    /// filters by cross-operand, skips shape-subshape pairs, and buckets into my_lists.
    pub fn prepare(&mut self) {
        // Clear all lists (OCCT L254-258)
        for list in &mut self.my_lists {
            list.clear();
        }

        let nv = self.ds.vertices.len();
        let ne = self.ds.edges.len();
        let nf = self.ds.faces.len();
        let total = nv + ne + nf;
        if total < 2 {
            return;
        }

        // Build AABBs for all shapes (verts, edges, faces)
        let mut indices = Vec::with_capacity(total);
        let mut aabbs = Vec::with_capacity(total);

        // Vertices (flat index 0..nv)
        for vi in 0..nv {
            indices.push(vi);
            let pt = self.ds.vertices[vi].point;
            let tol = self.ds.vertices[vi].geom_tol.max(1e-7) + self.ds.fuzzy_tol * 0.5;
            aabbs.push(Aabb { min: pt - glam::DVec3::splat(tol), max: pt + glam::DVec3::splat(tol) });
        }

        // Edges (flat index nv..nv+ne)
        for ei in 0..ne {
            let flat = nv + ei;
            indices.push(flat);
            let e = &self.ds.edges[ei];
            let pts = [self.ds.vertices[e.start_vertex].point, self.ds.vertices[e.end_vertex].point];
            let mut a = Aabb::empty();
            for &p in &pts { a.expand_point(p); }
            let tol = e.geom_tol.max(1e-7) + self.ds.fuzzy_tol * 0.5;
            a.min -= glam::DVec3::splat(tol);
            a.max += glam::DVec3::splat(tol);
            aabbs.push(a);
        }

        // Faces (flat index nv+ne..nv+ne+nf)
        for fi in 0..nf {
            let flat = nv + ne + fi;
            indices.push(flat);
            let f = &self.ds.faces[fi];
            let mut aabb = Aabb::empty();
            for &vi in &f.boundary_verts {
                if vi < nv {
                    aabb.expand_point(self.ds.vertices[vi].point);
                }
            }
            // Sphere: expand to full sphere volume
            if let rcad_kernel::geom::Surface3::Sphere(s) = &f.surface {
                let r = s.radius.abs();
                aabb.expand_point(s.center + glam::DVec3::splat(r));
                aabb.expand_point(s.center - glam::DVec3::splat(r));
            }
            let tol = f.geom_tol.max(1e-7) + self.ds.fuzzy_tol * 0.5;
            aabb.min -= glam::DVec3::splat(tol);
            aabb.max += glam::DVec3::splat(tol);
            aabbs.push(aabb);
        }

        // Build single BVH over all shapes (OCCT L276-291)
        let bvh = DsBvh::build(indices, aabbs);
        let pairs = DsBvh::candidate_pairs(&bvh, &bvh);

        let a_vc = self.ds.a_vertex_count;
        let a_ec = self.ds.a_edge_count;
        let a_fc = self.ds.a_face_count;

        // OCCT L300-357: iterate pairs, determine types, filter, bucket
        for &(ia, ib) in &pairs {
            if ia == ib { continue; }

            // Determine shape type and operand from flat index
            let (t1, t2, op1, op2, s1, s2) = if ia < nv && ib < nv {
                (ShapeType::Vertex, ShapeType::Vertex, ia < a_vc, ib < a_vc, ia, ib)
            } else if ia < nv && ib < nv + ne {
                if ib >= nv {
                    (ShapeType::Vertex, ShapeType::Edge, ia < a_vc, (ib - nv) < a_ec, ia, ib - nv)
                } else {
                    continue;
                }
            } else if ia >= nv && ia < nv + ne && ib < nv {
                (ShapeType::Edge, ShapeType::Vertex, (ia - nv) < a_ec, ib < a_vc, ia - nv, ib)
            } else if ia >= nv && ia < nv + ne && ib >= nv && ib < nv + ne {
                (ShapeType::Edge, ShapeType::Edge, (ia - nv) < a_ec, (ib - nv) < a_ec, ia - nv, ib - nv)
            } else if ia < nv && ib >= nv + ne {
                (ShapeType::Vertex, ShapeType::Face, ia < a_vc, (ib - nv - ne) < a_fc, ia, ib - nv - ne)
            } else if ia >= nv + ne && ib < nv {
                (ShapeType::Face, ShapeType::Vertex, (ia - nv - ne) < a_fc, ib < a_vc, ia - nv - ne, ib)
            } else if ia >= nv && ia < nv + ne && ib >= nv + ne {
                (ShapeType::Edge, ShapeType::Face, (ia - nv) < a_ec, (ib - nv - ne) < a_fc, ia - nv, ib - nv - ne)
            } else if ia >= nv + ne && ib >= nv && ib < nv + ne {
                (ShapeType::Face, ShapeType::Edge, (ia - nv - ne) < a_fc, (ib - nv) < a_ec, ia - nv - ne, ib - nv)
            } else if ia >= nv + ne && ib >= nv + ne {
                (ShapeType::Face, ShapeType::Face, (ia - nv - ne) < a_fc, (ib - nv - ne) < a_fc, ia - nv - ne, ib - nv - ne)
            } else {
                continue;
            };

            // Cross-operand filter: skip same-operand pairs (OCCT uses Ranges)
            if op1 == op2 { continue; }

            // OCCT L335-340: avoid interfering shape with its sub-shapes
            // (not applicable for the 3 shape type combos we handle)
            // Vertex-Edge: skip if vertex is sub-shape of edge
            if t1 == ShapeType::Vertex && t2 == ShapeType::Edge {
                if self.ds.edge_has_vertex(s1, s2) { continue; }
            }
            if t1 == ShapeType::Edge && t2 == ShapeType::Vertex {
                if self.ds.edge_has_vertex(s2, s1) { continue; }
            }
            // Edge-Face: skip if edge is sub-shape of face
            // Face has boundary_edges; only check for edges that are on the face boundary
            // (the HasSubShape check in OCCT prevents processing edge-face pairs where
            //  the edge is a boundary edge of the face, since boundary edges don't need
            //  EF intersection — they already share the face as a parent shape)

            // OCCT L342-352: optional OBB check (skipped for now)

            // Bucket by type combination (OCCT L354-356)
            let bucket = Self::type_to_bucket(t1, t2);
            if bucket >= 0 && (bucket as usize) < self.my_lists.len() {
                let key = if s1 <= s2 { (s1, s2) } else { (s2, s1) };
                self.my_lists[bucket as usize].push(key);
            }
        }

        // stable_sort each bucket (OCCT Initialize L203: std::stable_sort)
        for list in &mut self.my_lists {
            list.sort();
        }
    }

    /// OCCT BOPDS_Iterator::Initialize — select pairs of given type combination.
    ///
    /// Applies stable_sort (already done in Prepare) and sets up iteration.
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        let bucket = Self::type_to_bucket(t1, t2);
        if bucket >= 0 && (bucket as usize) < self.my_lists.len() {
            self.current_list = self.my_lists[bucket as usize].clone();
        } else {
            self.current_list.clear();
        }
        self.current_pos = 0;
    }

    /// OCCT BOPDS_Iterator::More
    pub fn more(&self) -> bool {
        self.current_pos < self.current_list.len()
    }

    /// OCCT BOPDS_Iterator::Next
    pub fn next(&mut self) {
        self.current_pos += 1;
    }

    /// OCCT BOPDS_Iterator::Value
    ///
    /// Returns (index1, index2) where index1 <= index2 (matching OCCT's
    /// BOPDS_Pair which stores min/max indices).
    pub fn value(&self) -> (usize, usize) {
        self.current_list[self.current_pos]
    }

    /// Returns the expected length (OCCT ExpectedLength)
    pub fn expected_length(&self) -> usize {
        self.current_list.len()
    }

    /// Returns a reference to the pre-computed pair list for (t1, t2).
    /// Must be called after `prepare()`.
    pub fn pairs(&self, t1: ShapeType, t2: ShapeType) -> &[(usize, usize)] {
        let bucket = Self::type_to_bucket(t1, t2);
        if bucket >= 0 && (bucket as usize) < self.my_lists.len() {
            &self.my_lists[bucket as usize]
        } else {
            &[]
        }
    }
}
