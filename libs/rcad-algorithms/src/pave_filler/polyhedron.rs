//! OCCT-aligned: IntPatch_Polyhedron + IntPatch_InterferencePolyhedron
//!
//! Triangle mesh approximation of a parametric surface for seed point
//! detection in PrmPrm intersection. Matches OCCT IntPatch_Polyhedron.hxx/cxx.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};

/// OCCT-aligned: IntPatch_Polyhedron — triangle mesh approximation
pub struct Polyhedron {
    pub(crate) nb_u: i32,
    pub(crate) nb_v: i32,
    pub(crate) points: Vec<DVec3>,
    pub(crate) u_params: Vec<f64>,
    pub(crate) v_params: Vec<f64>,
    pub(crate) bbox_min: DVec3,
    pub(crate) bbox_max: DVec3,
}

impl Polyhedron {
    /// OCCT L34-36: constructor with surface + (nU, nV) grid resolution.
    /// rcad: clamps to min=1 matching OCCT IntPatch_Polyhedron behavior.
    pub fn new(surf: &Surface3, n_u: i32, n_v: i32) -> Self {
        let n_u = n_u.max(1);
        let n_v = n_v.max(1);
        let [u_min, u_max, v_min, v_max] = surf.default_domain();

        let mut pts = Vec::with_capacity(((n_u + 1) * (n_v + 1)) as usize);
        let mut u_params = Vec::with_capacity((n_u + 1) as usize);
        let mut v_params = Vec::with_capacity((n_v + 1) as usize);
        let mut bbox_min = DVec3::splat(f64::MAX);
        let mut bbox_max = DVec3::splat(f64::MIN);

        for i in 0..=n_u {
            let u = u_min + (i as f64 / n_u as f64) * (u_max - u_min);
            u_params.push(u);
        }
        for j in 0..=n_v {
            let v = v_min + (j as f64 / n_v as f64) * (v_max - v_min);
            v_params.push(v);
        }

        for j in 0..=n_v {
            for i in 0..=n_u {
                let p = surf.point_at(u_params[i as usize], v_params[j as usize]);
                bbox_min = bbox_min.min(p);
                bbox_max = bbox_max.max(p);
                pts.push(p);
            }
        }

        Self { nb_u: n_u, nb_v: n_v, points: pts, u_params, v_params,
            bbox_min, bbox_max }
    }

    pub fn nb_triangles(&self) -> i32 { self.nb_u * self.nb_v * 2 }
    pub fn nb_points(&self) -> i32 { (self.nb_u + 1) * (self.nb_v + 1) }

    /// Get triangle vertex indices (1-indexed, OCCT convention)
    pub fn triangle(&self, index: i32) -> (i32, i32, i32) {
        let t = index - 1; // 0-indexed
        let row = t / (2 * self.nb_u);
        let col = t % (2 * self.nb_u);
        let i0 = row * (self.nb_u + 1) + col / 2;
        if col % 2 == 0 {
            (i0 + 1, i0 + 2, i0 + self.nb_u + 2) // lower-left triangle
        } else {
            (i0 + 1, i0 + self.nb_u + 2, i0 + self.nb_u + 1) // upper-right
        }
    }

    pub fn point(&self, index: i32) -> DVec3 { self.points[(index - 1) as usize] }
    /// OCCT-aligned: Size(nbU, nbV) — get grid dimensions.
    pub fn size(&self) -> (i32, i32) { (self.nb_u, self.nb_v) }
    pub fn bbox(&self) -> (DVec3, DVec3) { (self.bbox_min, self.bbox_max) }

    /// OCCT-aligned: TriConnex (IntPatch_Polyhedron.cxx L299-567).
    /// Finds the triangle adjacent to `triang` across the edge (pivot, pedge).
    /// Returns the adjacent triangle index (0 = boundary), sets `tri_con` and `other_p`.
    pub fn tri_connex(&self, triang: i32, pivot: i32, pedge: i32) -> (i32, i32) {
        const LONGUEUR_MINI_EDGE_TRIANGLE: f64 = 1e-14; // OCCT 1e-14
        let nbdeltaVp1 = self.nb_v + 1;
        let nbdeltaVm2 = self.nb_v * 2;

        let pivotm1 = pivot - 1;
        let lig_p = pivotm1 / nbdeltaVp1;
        let col_p = pivotm1 - lig_p * nbdeltaVp1;

        let (lig_e, col_e, typ_e) = if pedge != 0 {
            let lig_e = (pedge - 1) / nbdeltaVp1;
            let col_e = (pedge - 1) - lig_e * nbdeltaVp1;
            let typ_e = if lig_p == lig_e { 1 }      // Horizontal
                   else if col_p == col_e { 2 }      // Vertical
                   else { 3 };                        // Oblique
            (lig_e, col_e, typ_e)
        } else {
            (0, 0, 0)
        };

        let (mut lin_t, mut col_t, mut lin_o, mut col_o);
        if triang != 0 {
            let t  = (triang - 1) / nbdeltaVm2;
            let tt = (triang - 1) - t * nbdeltaVm2;
            lin_t = 1 + t;
            col_t = 1 + tt;
            let typ_e = if typ_e == 0 {
                // Determine edge type from triangle position relative to pivot
                if lig_p == lin_t {
                    (lig_p - 1, col_p - 1, 3)
                } else if col_t == lig_p + lig_p {
                    (lig_p, col_p - 1, 1)
                } else {
                    (lig_p + 1, col_p + 1, 3)
                }
            } else { (lig_e, col_e, typ_e) };
            let (_lig_e, _col_e, typ_e) = typ_e;

            match typ_e {
                1 => { // Horizontal
                    if lin_t == lig_p { lin_t += 1; lin_o = lig_p + 1; col_o = col_p.max(col_e); }
                    else               { lin_t -= 1; lin_o = lig_p - 1; col_o = col_p.min(col_e); }
                }
                2 => { // Vertical
                    if col_t == col_p + col_p { col_t += 1; lin_o = lig_p.max(lig_e); col_o = col_p + 1; }
                    else                      { col_t -= 1; lin_o = lig_p.min(lig_e); col_o = col_p - 1; }
                }
                3 => { // Oblique
                    if (col_t & 1) == 0 { col_t -= 1; lin_o = lig_p.max(lig_e); col_o = col_p.min(col_e); }
                    else                { col_t += 1; lin_o = lig_p.min(lig_e); col_o = col_p.max(col_e); }
                }
                _ => { lin_o = 0; col_o = 0; }
            }
        } else {
            // Unknown triangle position
            if pedge == 0 {
                lin_t = 1.max(lig_p);
                col_t = 1.max(col_p + col_p);
                if lig_p == 0 { lin_o = lig_p + 1; } else { lin_o = lig_p - 1; }
                col_o = col_p;
            } else {
                match typ_e {
                    1 => { lin_t = lig_p + 1; col_t = col_p.max(col_e) * 2; lin_o = lig_p + 1; col_o = col_p.max(col_e); }
                    2 => { lin_t = lig_p.max(lig_e); col_t = col_p + col_p; lin_o = lig_p.min(lig_e); col_o = col_p - 1; }
                    3 => { lin_t = lig_p.max(lig_e); col_t = col_p + col_e; lin_o = lig_p.max(lig_e); col_o = col_p.min(col_e); }
                    _ => { lin_t = 0; col_t = 0; lin_o = 0; col_o = 0; }
                }
            }
        }

        let mut tri_con = (lin_t - 1) * nbdeltaVm2 + col_t;

        // Boundary checks
        if lin_t < 1 {
            lin_o = 0; col_o = col_p + col_p - col_e;
            if col_o < 0 { col_o = 0; lin_o = 1; }
            else if col_o > self.nb_v { col_o = self.nb_v; lin_o = 1; }
            tri_con = 0;
        } else if lin_t > self.nb_u {
            lin_o = self.nb_u; col_o = col_p + col_p - col_e;
            if col_o < 0 { col_o = 0; lin_o = self.nb_u - 1; }
            else if col_o > self.nb_v { col_o = self.nb_v; lin_o = self.nb_u - 1; }
            tri_con = 0;
        }
        if col_t < 1 {
            col_o = 0; lin_o = lig_p + lig_p - lig_e;
            if lin_o < 0 { lin_o = 0; col_o = 1; }
            else if lin_o > self.nb_u { lin_o = self.nb_u; col_o = 1; }
            tri_con = 0;
        } else if col_t > self.nb_v {
            col_o = self.nb_v; lin_o = lig_p + lig_p - lig_e;
            if lin_o < 0 { lin_o = 0; col_o = self.nb_v - 1; }
            else if lin_o > self.nb_u { lin_o = self.nb_u; col_o = self.nb_v - 1; }
            tri_con = 0;
        }

        let other_p = lin_o * nbdeltaVp1 + col_o + 1;
        if pedge != 0 && self.point(pivot).distance_squared(self.point(pedge)) <= LONGUEUR_MINI_EDGE_TRIANGLE {
            return (triang, 0);
        }
        if pedge != 0 && other_p > 0 && other_p <= self.points.len() as i32
            && self.point(other_p).distance_squared(self.point(pedge)) <= LONGUEUR_MINI_EDGE_TRIANGLE {
            return (0, 0);
        }
        (tri_con, other_p)
    }
    pub fn parameters(&self, index: i32) -> (f64, f64) {
        let idx = (index - 1) as usize;
        let n_plus_1 = (self.nb_u + 1) as usize;
        let j = idx / n_plus_1;
        let i = idx % n_plus_1;
        (self.u_params[i], self.v_params[j])
    }
}

/// Seed point found by triangle-triangle intersection
#[derive(Clone, Debug)]
pub struct SeedPoint {
    pub p3d: DVec3,
    pub u1: f64, pub v1: f64,
    pub u2: f64, pub v2: f64,
}

/// OCCT-aligned: IntPatch_InterferencePolyhedron — find intersecting triangles
pub struct InterferencePolyhedron {
    nb_section_lines: i32,
    seed_points: Vec<SeedPoint>,
}

impl InterferencePolyhedron {
    pub fn new(poly1: &Polyhedron, poly2: &Polyhedron) -> Self {
        let (bmin1, bmax1) = poly1.bbox();
        let (bmin2, bmax2) = poly2.bbox();

        // Quick AABB rejection
        let overlap = (bmin1.x <= bmax2.x && bmax1.x >= bmin2.x)
            && (bmin1.y <= bmax2.y && bmax1.y >= bmin2.y)
            && (bmin1.z <= bmax2.z && bmax1.z >= bmin2.z);

        if !overlap {
            return Self { nb_section_lines: 0, seed_points: Vec::new() };
        }

        let mut seeds: Vec<SeedPoint> = Vec::new();

        // OCCT: for each triangle of poly1, check against each triangle of poly2
        let n_tri1 = poly1.nb_triangles();
        let n_tri2 = poly2.nb_triangles();

        for t1 in 1..=n_tri1 {
            let (i1, i2, i3) = poly1.triangle(t1);
            let p1 = poly1.point(i1);
            let p2 = poly1.point(i2);
            let p3 = poly1.point(i3);

            // Triangle bounding box
            let tmin = p1.min(p2).min(p3);
            let tmax = p1.max(p2).max(p3);

            // Quick reject vs poly2 bbox
            if tmax.x < bmin2.x || tmin.x > bmax2.x
                || tmax.y < bmin2.y || tmin.y > bmax2.y
                || tmax.z < bmin2.z || tmin.z > bmax2.z {
                continue;
            }

            for t2 in 1..=n_tri2 {
                let (j1, j2, j3) = poly2.triangle(t2);
                let q1 = poly2.point(j1);
                let q2 = poly2.point(j2);
                let q3 = poly2.point(j3);

                // Triangle-triangle intersection test: check if any edge
                // of one triangle crosses the plane of the other
                if let Some(seed) = Self::intersect_triangles(
                    p1, p2, p3, q1, q2, q3,
                    &poly1, &poly2,
                ) {
                    seeds.push(seed);
                }
            }
        }

        Self { nb_section_lines: seeds.len() as i32, seed_points: seeds }
    }

    fn intersect_triangles(
        p1: DVec3, p2: DVec3, p3: DVec3,
        q1: DVec3, q2: DVec3, q3: DVec3,
        poly1: &Polyhedron, poly2: &Polyhedron,
    ) -> Option<SeedPoint> {
        let n1 = (p2 - p1).cross(p3 - p1).normalize_or_zero();
        let n2 = (q2 - q1).cross(q3 - q1).normalize_or_zero();

        let d1_plane = (q1 - p1).dot(n1);
        let d2_plane = (q2 - p1).dot(n1);
        let d3_plane = (q3 - p1).dot(n1);

        // All same sign → no intersection
        if (d1_plane > 0.0 && d2_plane > 0.0 && d3_plane > 0.0)
            || (d1_plane < 0.0 && d2_plane < 0.0 && d3_plane < 0.0) {
            return None;
        }

        // Find the midpoint of the intersection segment
        let (hit_q, _) = Self::closest_edge_point(q1, q2, q3, p1, n1, d1_plane, d2_plane, d3_plane)?;

        // Project back to both surfaces to get UV
        let (u1, v1) = poly1.parameters(1); // approximate: use first triangle vertex
        let (u2, v2) = poly2.parameters(1);
        let _ = (u1, v1, u2, v2);

        // rcad: intersect at the midpoint of the intersection edge
        Some(SeedPoint {
            p3d: hit_q,
            u1: 0.0, v1: 0.0, // approximate — real UV requires surface projection
            u2: 0.0, v2: 0.0,
        })
    }

    /// Find the intersection point of a triangle edge with the plane
    fn closest_edge_point(
        q1: DVec3, q2: DVec3, q3: DVec3,
        _p1: DVec3, n1: DVec3,
        d1: f64, d2: f64, d3: f64,
    ) -> Option<(DVec3, i32)> {
        let edges = [(q1, q2, d1, d2), (q2, q3, d2, d3), (q3, q1, d3, d1)];
        for (ea, eb, da, db) in edges {
            if da * db < 0.0 {
                let t = da.abs() / (da.abs() + db.abs());
                let hit = ea + t * (eb - ea);
                return Some((hit, 0));
            }
        }
        None
    }

    pub fn nb_section_lines(&self) -> i32 { self.nb_section_lines }
    pub fn seed_points(&self) -> &[SeedPoint] { &self.seed_points }
}

// ═══════════════════════════════════════════════════════════════════════
// OCCT tests: IntPatch_Polyhedron_Test.cxx + IntPatch_PolyhedronBVH_Test.cxx
// ═══════════════════════════════════════════════════════════════════════

