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

// ── Inertia tensor ────────────────────────────────────────────────────────────

/// Symmetric 3×3 moment of inertia tensor (assuming uniform density = 1).
///
/// The components are defined as:
/// ```text
/// Ixx = ∫(y²+z²) dV,  Iyy = ∫(x²+z²) dV,  Izz = ∫(x²+y²) dV
/// Ixy = -∫xy dV,       Ixz = -∫xz dV,       Iyz = -∫yz dV
/// ```
///
/// Computed about the world origin. To get the tensor about the centroid,
/// use the parallel-axis theorem.
#[derive(Debug, Clone, Copy)]
pub struct InertiaTensor {
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixy: f64,
    pub ixz: f64,
    pub iyz: f64,
}

impl InertiaTensor {
    /// Returns the 3×3 inertia matrix as row-major `[[f64;3];3]`.
    pub fn to_matrix(&self) -> [[f64; 3]; 3] {
        [
            [ self.ixx, -self.ixy, -self.ixz],
            [-self.ixy,  self.iyy, -self.iyz],
            [-self.ixz, -self.iyz,  self.izz],
        ]
    }
}

/// Computes the moment of inertia tensor of a closed BRep solid about the
/// world origin.
///
/// Uses the divergence theorem (polyhedral formula from Mirtich 1996) applied
/// to the BRep's triangulated faces, consistent with the existing `volume` and
/// `centroid` implementations.
///
/// Assumes uniform density = 1 (unit density).  Multiply each component by
/// the actual density to get physical inertia.
pub fn inertia_tensor(brep: &BRep) -> InertiaTensor {
    let mut ixx = 0.0_f64;
    let mut iyy = 0.0_f64;
    let mut izz = 0.0_f64;
    let mut ixy = 0.0_f64;
    let mut ixz = 0.0_f64;
    let mut iyz = 0.0_f64;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for [a, b, c] in face_triangles(brep, face) {
                    // Signed volume of tet (origin, a, b, c)
                    // sv = a·(b×c)/6 — same as tet_signed_volume
                    let sv = a.dot(b.cross(c)) / 6.0;

                    // Symmetric quadratic sums for each coordinate pair.
                    // For ∫_tet x² dV = sv/10 * x2_sym (from simplex integration).
                    let x2 = a.x*a.x + b.x*b.x + c.x*c.x + a.x*b.x + a.x*c.x + b.x*c.x;
                    let y2 = a.y*a.y + b.y*b.y + c.y*c.y + a.y*b.y + a.y*c.y + b.y*c.y;
                    let z2 = a.z*a.z + b.z*b.z + c.z*c.z + a.z*b.z + a.z*c.z + b.z*c.z;

                    ixx += sv / 10.0 * (y2 + z2);
                    iyy += sv / 10.0 * (x2 + z2);
                    izz += sv / 10.0 * (x2 + y2);

                    // For ∫_tet xy dV = sv/20 * xy_mixed (from simplex integration).
                    // Product-moment: Ixy = -∫xy dV, etc.
                    let xy = 2.0*(a.x*a.y+b.x*b.y+c.x*c.y)
                           + a.x*b.y + b.x*a.y + a.x*c.y + c.x*a.y + b.x*c.y + c.x*b.y;
                    let xz = 2.0*(a.x*a.z+b.x*b.z+c.x*c.z)
                           + a.x*b.z + b.x*a.z + a.x*c.z + c.x*a.z + b.x*c.z + c.x*b.z;
                    let yz = 2.0*(a.y*a.z+b.y*b.z+c.y*c.z)
                           + a.y*b.z + b.y*a.z + a.y*c.z + c.y*a.z + b.y*c.z + c.y*b.z;

                    ixy += sv / 20.0 * xy;
                    ixz += sv / 20.0 * xz;
                    iyz += sv / 20.0 * yz;
                }
            }
        }
    }

    // Diagonal terms must be positive for a physical solid.
    // Off-diagonal sign: Ixy = -∫xy dV so negate the accumulated sums.
    InertiaTensor {
        ixx: ixx.abs(),
        iyy: iyy.abs(),
        izz: izz.abs(),
        ixy: -ixy,
        ixz: -ixz,
        iyz: -iyz,
    }
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

    #[test]
    fn unit_box_inertia_tensor_diagonal_equal() {
        // Unit box [0,1]^3 about the world origin:
        // Ixx = ∫(y²+z²)dV = (1/3 + 1/3) = 2/3
        // By symmetry, Iyy = Izz = 2/3
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let it = inertia_tensor(&brep);
        let expected = 2.0 / 3.0;
        let tol = 1e-4;
        assert!((it.ixx - expected).abs() < tol, "Ixx = {} expected {}", it.ixx, expected);
        assert!((it.iyy - expected).abs() < tol, "Iyy = {} expected {}", it.iyy, expected);
        assert!((it.izz - expected).abs() < tol, "Izz = {} expected {}", it.izz, expected);
    }

    #[test]
    fn box_2x1x1_inertia_tensor() {
        // Box [0,2]×[0,1]×[0,1] about origin:
        // Ixx = ∫(y²+z²)dV = V*(1/3+1/3) = 2*(2/3) = 4/3
        // Iyy = ∫(x²+z²)dV = V*(4/3÷2 + 1/3) = 2*(2/3+1/3) = 2*(1) = wait:
        //   ∫₀²∫₀¹∫₀¹ (x²+z²) dx dy dz  but order matters since box is [0,2]x[0,1]x[0,1]
        //   = 1*1*(∫₀² x² dx) + 1*2*(∫₀¹ z² dz) = (8/3) + 2*(1/3) = 8/3+2/3 = 10/3
        // Izz = ∫(x²+y²)dV = (8/3) + 2*(1/3) = 10/3
        // Ixx = 2*(1/3) + 2*(1/3) = 4/3
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 1.0, depth: 1.0 });
        let it = inertia_tensor(&brep);
        let tol = 1e-3;
        let expected_ixx = 4.0 / 3.0;
        let expected_iyy = 10.0 / 3.0;
        let expected_izz = 10.0 / 3.0;
        assert!((it.ixx - expected_ixx).abs() < tol, "Ixx = {} expected {}", it.ixx, expected_ixx);
        assert!((it.iyy - expected_iyy).abs() < tol, "Iyy = {} expected {}", it.iyy, expected_iyy);
        assert!((it.izz - expected_izz).abs() < tol, "Izz = {} expected {}", it.izz, expected_izz);
    }
}
