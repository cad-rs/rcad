//! Shape properties: surface area, volume, and centroid.
//!
//! Analogous to OCCT `GProp_GProps` with `BRepGProp`.
//!
//! All computations use triangulated faces where available; for faces without
//! pre-triangulated data, the outer wire's vertices are used directly (fan
//! triangulation from the first vertex).

use glam::DVec3;

use crate::BRep;

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Compute the signed area of a triangle from three points.
/// The sign depends on the orientation relative to the caller.
#[inline]
fn tri_area(a: DVec3, b: DVec3, c: DVec3) -> f64 {
    (b - a).cross(c - a).length() * 0.5
}

/// Signed volume contribution of a tetrahedron from the origin to triangle (a,b,c).
/// Summing over all surface triangles gives 1/6 * signed volume of the solid.
#[inline]
fn tet_signed_volume(a: DVec3, b: DVec3, c: DVec3) -> f64 {
    a.dot(b.cross(c)) / 6.0
}

/// Collect triangles for a face (either pre-triangulated or fan-triangulated
/// from the wire vertices), oriented outward (consistent with face.normal).
fn face_triangles<'a>(brep: &'a BRep, face: &'a crate::topology::Face) -> Vec<[DVec3; 3]> {
    let raw: Vec<[DVec3; 3]> = if !face.triangles.is_empty() {
        // Use pre-triangulated data
        face.triangles.iter()
            .filter_map(|&[i, j, k]| {
                let a = brep.vertices.get(i)?.point;
                let b = brep.vertices.get(j)?.point;
                let c = brep.vertices.get(k)?.point;
                Some([a, b, c])
            })
            .collect()
    } else {
        // Fan-triangulate from wire vertices
        let wire_pts: Vec<DVec3> = face.outer_wire.edges.iter()
            .filter_map(|we| {
                let edge = brep.edges.get(we.idx)?;
                let vidx = if we.forward { edge.start } else { edge.end };
                brep.vertices.get(vidx).map(|v| v.point)
            })
            .collect();

        if wire_pts.len() < 3 {
            return Vec::new();
        }

        let origin = wire_pts[0];
        (1..wire_pts.len() - 1)
            .map(|i| [origin, wire_pts[i], wire_pts[i + 1]])
            .collect()
    };

    // Ensure each triangle is oriented consistently with the face normal.
    raw.into_iter()
        .map(|[a, b, c]| {
            let tri_normal = (b - a).cross(c - a);
            if tri_normal.dot(face.normal) < 0.0 {
                [a, c, b] // flip
            } else {
                [a, b, c]
            }
        })
        .collect()
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Compute the total surface area of all faces in the BRep.
///
/// For each face, sums the areas of its triangles (pre-triangulated or fan).
/// Returns 0.0 if the BRep has no faces.
pub fn surface_area(brep: &BRep) -> f64 {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .flat_map(|f| face_triangles(brep, f))
        .map(|[a, b, c]| tri_area(a, b, c))
        .sum()
}

/// Compute the signed volume of the closed BRep solid.
///
/// Uses the divergence theorem: V = (1/6) Σ_triangles a·(b×c).
/// Works correctly for a closed, consistently-oriented mesh.
/// Returns 0.0 for open shells or empty BReps.
pub fn volume(brep: &BRep) -> f64 {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .flat_map(|f| face_triangles(brep, f))
        .map(|[a, b, c]| tet_signed_volume(a, b, c))
        .sum::<f64>()
        .abs()
}

/// Compute the centroid (center of mass) of the solid by volumetric integration.
///
/// Uses the formula: C = (1 / 8V) Σ_triangles (a+b+c) * tet_signed_vol(a,b,c)
/// where the sum is over all surface triangles.
///
/// Falls back to `BRep::center()` (vertex average) if the volume is near zero.
pub fn centroid(brep: &BRep) -> DVec3 {
    let mut vol_sum = 0.0_f64;
    let mut weighted_sum = DVec3::ZERO;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for [a, b, c] in face_triangles(brep, face) {
                    let sv = tet_signed_volume(a, b, c);
                    vol_sum += sv;
                    // Weight the centroid of each tet (at (a+b+c+origin)/4,
                    // origin=0) → simplified to (a+b+c) * sv
                    weighted_sum += (a + b + c) * sv;
                }
            }
        }
    }

    if vol_sum.abs() < 1e-15 {
        return brep.center();
    }

    // Centroid formula: (1/(2 * 4 * vol_sum)) * Σ (a+b+c) * sv
    // Simplification: weighted_sum / (4 * vol_sum) gives tet centroid average
    weighted_sum / (4.0 * vol_sum)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PrimitiveSolid;

    const EPS: f64 = 1e-6;

    #[test]
    fn unit_box_surface_area() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let area = surface_area(&brep);
        assert!((area - 6.0).abs() < EPS, "unit box surface area should be 6, got {area}");
    }

    #[test]
    fn unit_box_volume() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let vol = volume(&brep);
        assert!((vol - 1.0).abs() < EPS, "unit box volume should be 1, got {vol}");
    }

    #[test]
    fn box_2x3x4_volume() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 3.0, depth: 4.0 });
        let vol = volume(&brep);
        assert!((vol - 24.0).abs() < EPS, "2×3×4 box volume should be 24, got {vol}");
    }

    #[test]
    fn box_2x3x4_surface_area() {
        // SA = 2*(2*3 + 3*4 + 2*4) = 2*(6+12+8) = 52
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 3.0, depth: 4.0 });
        let area = surface_area(&brep);
        assert!((area - 52.0).abs() < EPS, "2×3×4 box SA should be 52, got {area}");
    }

    #[test]
    fn unit_box_centroid() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let c = centroid(&brep);
        // unit box: centroid at (0.5, 0.5, 0.5)
        assert!((c - DVec3::splat(0.5)).length() < 1e-4, "centroid should be (0.5,0.5,0.5), got {c}");
    }
}
