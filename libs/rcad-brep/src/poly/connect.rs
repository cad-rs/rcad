//! OCCT Poly_Connect — the triangulation adjacency walker (1:1 translation).
//!
//! Source: `OCCT src/FoundationClasses/TKMath/Poly/Poly_Connect.hxx` (L63-161)
//! and `Poly_Connect.cxx` (L34-296). Provides the algorithm to explore, inside
//! a triangulation, the adjacency data for a node or a triangle. For a
//! triangle T with nodes n1, n2, n3 the adjacent triangles AT1/AT2/AT3 share
//! the node pairs (n2,n3), (n3,n1), (n1,n2) respectively, and the adjacent
//! nodes an1/an2/an3 are the third nodes of those triangles (hxx L34-48).
//! This is the primitive the mesh checker's free-link analysis is built on
//! (`MeshTest_CheckTopology.cxx`).

use std::collections::HashSet;

use rcad_kernel::poly::PolyTriangulation;

/// One of the edges starting from a node. OCCT file-scope struct `polyedge`
/// (Poly_Connect.cxx L23-30): `next` (the next edge in the list), `nt[2]`
/// (the two adjacent triangles), `nn[2]` (the two adjacent nodes), `nd` (the
/// second node of the edge). The OCCT intrusive linked list (head per node,
/// prepended on insertion) is carried by a per-node `Vec` where new entries
/// are inserted at position 0 — the same record order the `->next` chain
/// produces, so the linear lookups below see the records in the same order.
#[derive(Debug, Clone, Copy)]
struct PolyEdge {
    /// OCCT `int nd` — the second node of the edge.
    nd: i32,
    /// OCCT `int nt[2]` — the two adjacent triangles.
    nt: [i32; 2],
    /// OCCT `int nn[2]` — the two adjacent nodes.
    nn: [i32; 2],
}

/// Provides an algorithm to explore, inside a triangulation, the adjacency
/// data for a node or a triangle. OCCT Poly_Connect.hxx L63-161.
pub struct PolyConnect<'a> {
    /// OCCT `occ::handle<Poly_Triangulation> myTriangulation`.
    my_triangulation: &'a PolyTriangulation,
    /// OCCT `NCollection_Array1<int> myTriangles` (1-based) — for each node,
    /// the index of a triangle containing it.
    my_triangles: Vec<i32>,
    /// OCCT `NCollection_Array1<int> myAdjacents` (1-based) — six entries per
    /// triangle: t1,t2,t3 (adjacent triangles) then n1,n2,n3 (adjacent nodes).
    my_adjacents: Vec<i32>,
    /// OCCT `int mytr`.
    mytr: i32,
    /// OCCT `int myfirst`.
    myfirst: i32,
    /// OCCT `int mynode`.
    mynode: i32,
    /// OCCT `int myothernode`.
    myothernode: i32,
    /// OCCT `bool mysense`.
    mysense: bool,
    /// OCCT `bool mymore`.
    mymore: bool,
    /// OCCT `TColStd_PackedMapOfInteger myPassedTr`.
    my_passed_tr: HashSet<i32>,
}

impl<'a> PolyConnect<'a> {
    /// Constructs an algorithm to explore the adjacency data of nodes or
    /// triangles for the triangulation theTriangulation.
    /// OCCT Poly_Connect.cxx L46-58 (the uninitialized default constructor
    /// L34-42 is a delayed-init artifact; `load` re-binds the same state).
    pub fn new(the_triangulation: &'a PolyTriangulation) -> Self {
        let mut c = PolyConnect {
            my_triangulation: the_triangulation,
            my_triangles: Vec::new(),
            my_adjacents: Vec::new(),
            mytr: 0,
            myfirst: 0,
            mynode: 0,
            myothernode: 0,
            mysense: false,
            mymore: false,
            my_passed_tr: HashSet::new(),
        };
        c.load(the_triangulation);
        c
    }

    /// Initialize the algorithm to explore the adjacency data of nodes or
    /// triangles for the triangulation theTriangulation.
    /// OCCT Poly_Connect.cxx L62-204 (Load).
    pub fn load(&mut self, the_triangulation: &'a PolyTriangulation) {
        self.my_triangulation = the_triangulation;
        self.mytr = 0;
        self.myfirst = 0;
        self.mynode = 0;
        self.myothernode = 0;
        self.mysense = false;
        self.mymore = false;

        let a_nb_nodes = self.my_triangulation.nb_nodes();
        let a_nb_tris = self.my_triangulation.nb_triangles();
        {
            let a_nb_adjs = 6 * a_nb_tris;
            if self.my_triangles.len() != a_nb_nodes {
                self.my_triangles.resize(a_nb_nodes, 0);
            }
            if self.my_adjacents.len() != a_nb_adjs {
                self.my_adjacents.resize(a_nb_adjs, 0);
            }
        }

        // L86-87: myTriangles.Init(0); myAdjacents.Init(0);
        self.my_triangles.iter_mut().for_each(|v| *v = 0);
        self.my_adjacents.iter_mut().for_each(|v| *v = 0);

        // L89-94: an array of the lists of edges connected to the nodes
        // (NCollection_Array1<polyedge*> anEdges(1, aNbNodes), Init(nullptr);
        // the incremental allocator is not needed — the Vec owns the records).
        let mut an_edges: Vec<Vec<PolyEdge>> = vec![Vec::new(); a_nb_nodes];

        // L96-98: loop on the triangles
        for a_tri_iter in 1..=a_nb_tris as i32 {
            // L102: get the nodes
            let a_tri_nodes = self.triangle_nodes(a_tri_iter);

            // L104-107: update the myTriangles array
            self.my_triangles[(a_tri_nodes[0] - 1) as usize] = a_tri_iter;
            self.my_triangles[(a_tri_nodes[1] - 1) as usize] = a_tri_iter;
            self.my_triangles[(a_tri_nodes[2] - 1) as usize] = a_tri_iter;

            // L109-151: update the edge lists
            for a_node_in_tri in 0..3usize {
                let a_node_next = (a_node_in_tri + 1) % 3; // the following node of the edge
                let an_edge_nodes = if a_tri_nodes[a_node_in_tri] < a_tri_nodes[a_node_next] {
                    [a_tri_nodes[a_node_in_tri], a_tri_nodes[a_node_next]]
                } else {
                    [a_tri_nodes[a_node_next], a_tri_nodes[a_node_in_tri]]
                };

                // L124-137: edge from node 0 to node 1 with node 0 < node 1;
                // insert in the list of node 0; scan for an existing record.
                let list = &mut an_edges[(an_edge_nodes[0] - 1) as usize];
                // the third node of the triangle (L134/L147: 3 - aNodeInTri - aNodeNext)
                let a_third_node = a_tri_nodes[3 - a_node_in_tri - a_node_next];
                let mut found = false;
                for ced in list.iter_mut() {
                    // the edge already exists
                    if ced.nd == an_edge_nodes[1] {
                        // just mark the adjacency if found
                        ced.nt[1] = a_tri_iter;
                        ced.nn[1] = a_third_node;
                        found = true;
                        break;
                    }
                }

                if !found {
                    // L139-150: create the edge if not found (prepended to
                    // the list: ced->next = anEdges[n]; anEdges[n] = ced)
                    list.insert(
                        0,
                        PolyEdge {
                            nd: an_edge_nodes[1],
                            nt: [a_tri_iter, 0],
                            nn: [a_third_node, 0],
                        },
                    );
                }
            }
        }

        // L154-192: now complete the myAdjacents array
        let mut an_adj_index = 0usize; // OCCT 1-based counter starting at 1
        for a_tri_iter in 1..=a_nb_tris as i32 {
            // L159: get the nodes
            let a_tri_nodes = self.triangle_nodes(a_tri_iter);

            // L162-190: for each edge in triangle
            for a_node_in_tri in 0..3usize {
                let a_node_next = (a_node_in_tri + 1) % 3; // the following node of the edge
                let an_edge_nodes = if a_tri_nodes[a_node_in_tri] < a_tri_nodes[a_node_next] {
                    [a_tri_nodes[a_node_in_tri], a_tri_nodes[a_node_next]]
                } else {
                    [a_tri_nodes[a_node_next], a_tri_nodes[a_node_in_tri]]
                };

                // L178-182: find in the list of node 0 (the OCCT `while
                // (ced->nd != anEdgeNodes[1]) ced = ced->next;` walk assumes
                // the record exists — same linear scan here).
                let list = &an_edges[(an_edge_nodes[0] - 1) as usize];
                let ced = list
                    .iter()
                    .find(|ced| ced.nd == an_edge_nodes[1])
                    .expect("Poly_Connect::Load: edge record must exist");

                // L185: Find the adjacent triangle
                let l = if ced.nt[0] == a_tri_iter { 1usize } else { 0usize };

                // OCCT 1-based anAdjIndex maps to 0-based Vec index.
                self.my_adjacents[an_adj_index] = ced.nt[l];
                self.my_adjacents[an_adj_index + 3] = ced.nn[l];
                an_adj_index += 1;
            }
            an_adj_index += 3; // L191: skip the n1/n2/n3 slots
        }
    }

    /// Returns the triangulation analyzed by this tool. OCCT hxx L80.
    pub fn triangulation(&self) -> &PolyTriangulation {
        self.my_triangulation
    }

    /// Returns the index of a triangle containing the node at index N (the
    /// null index 0 when the node is unused). OCCT hxx L84.
    pub fn triangle(&self, n: usize) -> i32 {
        self.my_triangles[n - 1]
    }

    /// Returns in t1, t2 and t3, the indices of the 3 triangles adjacent to
    /// the triangle at index T. Null indices are returned when there are
    /// fewer than 3 adjacent triangles. OCCT hxx L92-98.
    pub fn triangles(&self, t: usize) -> (i32, i32, i32) {
        let index = 6 * (t - 1);
        (
            self.my_adjacents[index],
            self.my_adjacents[index + 1],
            self.my_adjacents[index + 2],
        )
    }

    /// Returns, in n1, n2 and n3, the indices of the 3 nodes adjacent to the
    /// triangle referenced at index T. Null indices are returned when there
    /// are fewer than 3 adjacent nodes. OCCT hxx L105-111.
    pub fn nodes(&self, t: usize) -> (i32, i32, i32) {
        let index = 6 * (t - 1);
        (
            self.my_adjacents[index + 3],
            self.my_adjacents[index + 4],
            self.my_adjacents[index + 5],
        )
    }

    /// Initializes an iterator to search for all the triangles containing
    /// the node referenced at index N. OCCT Poly_Connect.cxx L208-230.
    pub fn initialize(&mut self, n: usize) {
        self.mynode = n as i32;
        self.myfirst = self.triangle(n);
        self.mytr = self.myfirst;
        self.mysense = true;
        self.mymore = self.myfirst != 0;
        self.my_passed_tr.clear();
        self.my_passed_tr.insert(self.mytr);
        if self.more() {
            let no = self.triangle_nodes(self.myfirst);
            let mut i = 0usize;
            while i < 3 {
                if no[i] == self.mynode {
                    break;
                }
                i += 1;
            }
            self.myothernode = no[(i + 2) % 3];
        }
    }

    /// Returns true if there is another element in the iterator. OCCT hxx L136.
    pub fn more(&self) -> bool {
        self.mymore
    }

    /// Advances the iterator defined with the function Initialize to access
    /// the next triangle. Note: There is no action if the iterator is empty.
    /// OCCT Poly_Connect.cxx L234-296 (Next).
    pub fn next(&mut self) {
        let mut t = [0i32; 3];
        (t[0], t[1], t[2]) = self.triangles(self.mytr as usize);
        if self.mysense {
            for &ti in &t {
                if ti != 0 {
                    let n = self.triangle_nodes(ti);
                    for j in 0..3usize {
                        if n[j] == self.mynode && n[(j + 1) % 3] == self.myothernode {
                            self.mytr = ti;
                            self.myothernode = n[(j + 2) % 3];
                            self.my_more_advance(ti);
                            return;
                        }
                    }
                }
            }
            // L260: sinon, depart vers la gauche. ("otherwise, leave to the left")
            let n = self.triangle_nodes(self.myfirst);
            let mut i = 0usize;
            while i < 3 {
                if n[i] == self.mynode {
                    break;
                }
                i += 1;
            }
            self.myothernode = n[(i + 1) % 3];
            self.mysense = false;
            self.mytr = self.myfirst;
            (t[0], t[1], t[2]) = self.triangles(self.mytr as usize);
        }
        if !self.mysense {
            for &ti in &t {
                if ti != 0 {
                    let n = self.triangle_nodes(ti);
                    for j in 0..3usize {
                        if n[j] == self.mynode && n[(j + 2) % 3] == self.myothernode {
                            self.mytr = ti;
                            self.myothernode = n[(j + 1) % 3];
                            self.my_more_advance(ti);
                            return;
                        }
                    }
                }
            }
        }
        self.mymore = false;
    }

    /// Returns the index of the current triangle to which the iterator
    /// points. OCCT hxx L148.
    pub fn value(&self) -> i32 {
        self.mytr
    }

    /// The triangle node triple of triangle T (OCCT
    /// `myTriangulation->Triangle(T).Get(n0, n1, n2)` at cxx L102/L159/L220
    /// etc.).
    fn triangle_nodes(&self, t: i32) -> [i32; 3] {
        let (a, b, c) = self
            .my_triangulation
            .triangle(t as usize)
            .get();
        [a, b, c]
    }

    /// cxx L253-254 / L287-288:
    /// `mymore = !myPassedTr.Contains(mytr); myPassedTr.Add(mytr);`
    fn my_more_advance(&mut self, ti: i32) {
        self.mymore = !self.my_passed_tr.contains(&ti);
        self.my_passed_tr.insert(ti);
    }
}
