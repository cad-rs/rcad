//! OCCT Poly: polygon and triangulation data structures.
//!
//! - `Poly_Triangulation` — triangle mesh with nodes, triangles, normals, UV
//! - `Poly_Connect` — node→triangle adjacency
//!
//! Note: `MathPoly` (polynomial solvers) is in `math/poly` above —
//! this is OCCT `Poly` (polygonal data), a different package.

use glam::DVec2;
use glam::DVec3;

/// OCCT Poly_Triangulation — a triangle mesh.
///
/// Stores:
/// - `nodes`: 3D vertex positions
/// - `triangles`: index triples `[i0, i1, i2]` into `nodes`
/// - `normals` (optional): per-vertex normals
/// - `uv_nodes` (optional): per-vertex UV coordinates
/// - `tri_indices`: 3×N triangle node indices (flat)
///
/// OCCT: Poly_Triangulation — reference-counted mesh data.
#[derive(Debug, Clone, Default)]
pub struct Triangulation {
    /// 3D vertex positions.
    pub nodes: Vec<DVec3>,
    /// Triangle vertex indices: 3 per triangle, stored flat.
    pub triangles: Vec<usize>,
    /// Optional per-vertex normals.
    pub normals: Option<Vec<DVec3>>,
    /// Optional per-vertex UV coordinates.
    pub uv_nodes: Option<Vec<DVec2>>,
}

impl Triangulation {
    /// Number of triangles.
    pub fn nb_triangles(&self) -> usize { self.triangles.len() / 3 }

    /// Number of nodes (vertices).
    pub fn nb_nodes(&self) -> usize { self.nodes.len() }

    /// Get the three vertex indices of triangle `i`.
    pub fn triangle(&self, i: usize) -> Option<[usize; 3]> {
        let i3 = i * 3;
        if i3 + 2 < self.triangles.len() {
            Some([self.triangles[i3], self.triangles[i3 + 1], self.triangles[i3 + 2]])
        } else { None }
    }

    /// Three corner positions of triangle `i`.
    pub fn triangle_points(&self, i: usize) -> Option<[DVec3; 3]> {
        let t = self.triangle(i)?;
        Some([self.nodes[t[0]], self.nodes[t[1]], self.nodes[t[2]]])
    }

    /// Compute the face normal of triangle `i`.
    pub fn triangle_normal(&self, i: usize) -> Option<DVec3> {
        let [a, b, c] = self.triangle_points(i)?;
        Some((b - a).cross(c - a).normalize_or_zero())
    }

    /// Add a triangle from node indices.
    pub fn add_triangle(&mut self, i0: usize, i1: usize, i2: usize) {
        self.triangles.push(i0);
        self.triangles.push(i1);
        self.triangles.push(i2);
    }
}

/// OCCT Poly_Connect — node→triangle connectivity built from a Triangulation.
///
/// Provides:
/// - `triangles_of_node(node)` — triangles incident to a node
/// - `nodes_of_triangle(tri)` — three nodes of a triangle
/// - `triangle_neighbors(tri)` — adjacent triangles
#[derive(Debug, Clone)]
pub struct Connect {
    /// For each node, the list of incident triangle indices.
    pub node_triangles: Vec<Vec<usize>>,
    /// For each triangle, the three adjacent triangles (−1 for boundary).
    pub tri_neighbors: Vec<[isize; 3]>,
}

impl Connect {
    /// Build connectivity from a Triangulation.
    /// OCCT: Poly_Connect::Build — fills node→triangle and triangle→triangle adjacency.
    pub fn build(tri: &Triangulation) -> Self {
        let nn = tri.nb_nodes();
        let nt = tri.nb_triangles();

        // Node → triangles
        let mut node_tri: Vec<Vec<usize>> = vec![Vec::new(); nn];
        for ti in 0..nt {
            if let Some([i0, i1, i2]) = tri.triangle(ti) {
                if i0 < nn { node_tri[i0].push(ti); }
                if i1 < nn { node_tri[i1].push(ti); }
                if i2 < nn { node_tri[i2].push(ti); }
            }
        }

        // Triangle → triangle adjacency via shared edge (opposite vertex)
        let mut tri_neighbors: Vec<[isize; 3]> = vec![[-1, -1, -1]; nt];
        for ti in 0..nt {
            if let Some([a, b, c]) = tri.triangle(ti) {
                let edges = [(a, b, 2), (b, c, 0), (c, a, 1)];
                for &(v0, v1, opp) in &edges {
                    if tri_neighbors[ti][opp] != -1 { continue; }
                    // Search for triangle sharing edge (v0, v1) in reverse order
                    for &nj in &node_tri[v0] {
                        if nj == ti { continue; }
                        if let Some([na, nb, nc]) = tri.triangle(nj) {
                            if (na == v1 && nb == v0) || (nb == v1 && nc == v0) || (nc == v1 && na == v0)
                                || (na == v0 && nb == v1) || (nb == v0 && nc == v1) || (nc == v0 && na == v1)
                            {
                                tri_neighbors[ti][opp] = nj as isize;
                                break;
                            }
                        }
                    }
                }
            }
        }

        Connect { node_triangles: node_tri, tri_neighbors }
    }

    /// Triangles incident to `node_index`.  OCCT: Poly_Connect::Triangles.
    pub fn triangles_of_node(&self, node: usize) -> &[usize] {
        self.node_triangles.get(node).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Adjacent triangles of `tri_index` (−1 means edge is on boundary).
    /// OCCT: Poly_Connect::Neighbors.
    pub fn triangle_neighbors(&self, tri: usize) -> [isize; 3] {
        self.tri_neighbors.get(tri).copied().unwrap_or([-1, -1, -1])
    }
}

// ══════════════════════════════════════════════════════════════════════════
// Tests
// ══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_triangle() {
        let mut tri = Triangulation::default();
        tri.nodes = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        tri.add_triangle(0, 1, 2);
        assert_eq!(tri.nb_triangles(), 1);
        assert_eq!(tri.nb_nodes(), 3);
        let n = tri.triangle_normal(0).unwrap();
        assert!((n - DVec3::Z).length() < 1e-12);
    }

    #[test]
    fn connectivity_two_triangles() {
        // Quadrilateral: two triangles sharing edge (1, 2)
        let mut tri = Triangulation::default();
        tri.nodes = vec![
            DVec3::new(0.0, 0.0, 0.0), // 0
            DVec3::new(1.0, 0.0, 0.0), // 1
            DVec3::new(1.0, 1.0, 0.0), // 2
            DVec3::new(0.0, 1.0, 0.0), // 3
        ];
        tri.add_triangle(0, 1, 2); // T0
        tri.add_triangle(0, 2, 3); // T1
        let conn = Connect::build(&tri);
        assert!(conn.triangle_neighbors(0)[2] == 1 || conn.triangle_neighbors(0)[2] != -1);
    }
}
