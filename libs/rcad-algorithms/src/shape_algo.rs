//! ShapeAlgo-style additional shape algorithms.
//!
//! Provides utilities for shape analysis and geometry extraction, analogous to
//! OCCT `ShapeAlgo` package. This module includes:
//!
//! - `AlgoContainer`: Container for pluggable shape algorithms
//! - `GetBoxGeometry`: Extract box dimensions from a rcad_kernel::BRep
//! - `GetCylinderGeometry`: Extract cylinder parameters from a rcad_kernel::BRep
//! - `GetSphereGeometry`: Extract sphere parameters from a rcad_kernel::BRep
//! - `GetConeGeometry`: Extract cone parameters from a rcad_kernel::BRep
//! - `GetTorusGeometry`: Extract torus parameters from a rcad_kernel::BRep
//! - `IsPrimitive`: Check if a shape matches a primitive type

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::geom::{
    ConicalSurface, CylindricalSurface, Plane, SphericalSurface, Surface3, ToroidalSurface,
};
use rcad_kernel::topods::{BRep, ShapeRef, TShape};
use std::collections::HashMap;

// =============================================================================
// Geometry Extraction Structures
// =============================================================================

/// Extracted box geometry parameters.
///
/// Represents an axis-aligned box with origin at the minimum corner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoxGeometry {
    /// Origin (minimum corner) of the box.
    pub origin: DVec3,
    /// Dimension along the X axis.
    pub dx: f64,
    /// Dimension along the Y axis.
    pub dy: f64,
    /// Dimension along the Z axis.
    pub dz: f64,
}

/// Extracted cylinder geometry parameters.
///
/// Represents a cylinder defined by origin (center of bottom), axis, radius, and height.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderGeometry {
    /// Origin at the center of the bottom face.
    pub origin: DVec3,
    /// Cylinder axis direction (normalized).
    pub axis: DVec3,
    /// Cylinder radius.
    pub radius: f64,
    /// Cylinder height along the axis.
    pub height: f64,
}

/// Extracted sphere geometry parameters.
///
/// Represents a sphere defined by center and radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereGeometry {
    /// Center of the sphere.
    pub center: DVec3,
    /// Sphere radius.
    pub radius: f64,
}

/// Extracted cone geometry parameters.
///
/// Represents a cone defined by apex, axis, and half-angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConeGeometry {
    /// Apex point of the cone.
    pub apex: DVec3,
    /// Cone axis direction (normalized).
    pub axis: DVec3,
    /// Half-angle of the cone (radians).
    pub angle: f64,
}

/// Extracted torus geometry parameters.
///
/// Represents a torus defined by center, axis, and two radii.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorusGeometry {
    /// Center of the torus.
    pub center: DVec3,
    /// Torus axis direction (normalized).
    pub axis: DVec3,
    /// Distance from center to the center of the tube.
    pub major_radius: f64,
    /// Radius of the tube.
    pub minor_radius: f64,
}

// =============================================================================
// ShapeAlgorithm Trait
// =============================================================================

/// Trait for pluggable shape algorithms.
///
/// Algorithms implementing this trait can be registered with an `AlgoContainer`
/// and executed on rcad_kernel::BRep shapes.
pub trait ShapeAlgorithm: Send + Sync {
    /// Get the name of this algorithm.
    fn name(&self) -> &str;

    /// Execute the algorithm on the given rcad_kernel::BRep.
    ///
    /// Returns `true` if the algorithm succeeded, `false` otherwise.
    fn execute(&self, brep: &rcad_kernel::BRep) -> bool;
}

// =============================================================================
// AlgoContainer
// =============================================================================

/// Container for pluggable shape algorithms.
///
/// Provides a registry for algorithms that can be looked up by name and
/// executed on rcad_kernel::BRep shapes. Analogous to OCCT `ShapeAlgo_AlgoContainer`.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::shape_algo::{AlgoContainer, ShapeAlgorithm};
///
/// let mut container = AlgoContainer::new();
/// // Algorithms can be added via add_algorithm
/// ```
pub struct AlgoContainer {
    /// Registered algorithms indexed by name.
    algorithms: HashMap<String, Box<dyn ShapeAlgorithm>>,
}

impl AlgoContainer {
    /// Create a new empty algorithm container.
    pub fn new() -> Self {
        Self {
            algorithms: HashMap::new(),
        }
    }

    /// Add an algorithm to the container.
    ///
    /// If an algorithm with the same name already exists, it will be replaced.
    pub fn add_algorithm(&mut self, name: &str, algorithm: Box<dyn ShapeAlgorithm>) {
        self.algorithms.insert(name.to_string(), algorithm);
    }

    /// Get an algorithm by name.
    pub fn get_algorithm(&self, name: &str) -> Option<&dyn ShapeAlgorithm> {
        self.algorithms.get(name).map(|b| b.as_ref())
    }

    /// Check if an algorithm exists.
    pub fn has_algorithm(&self, name: &str) -> bool {
        self.algorithms.contains_key(name)
    }

    /// Remove an algorithm by name.
    pub fn remove_algorithm(&mut self, name: &str) -> bool {
        self.algorithms.remove(name).is_some()
    }

    /// Get the number of registered algorithms.
    pub fn len(&self) -> usize {
        self.algorithms.len()
    }

    /// Check if the container is empty.
    pub fn is_empty(&self) -> bool {
        self.algorithms.is_empty()
    }

    /// Get all algorithm names.
    pub fn algorithm_names(&self) -> Vec<&str> {
        self.algorithms.keys().map(|s| s.as_str()).collect()
    }
}

impl Default for AlgoContainer {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// Geometry Extraction Functions
// =============================================================================

/// Collect all face surfaces from the BRep by iterating TShape hierarchy.
///
/// In the new topods::BRep API, surfaces live directly on TShape::Face data.
fn collect_face_surfaces(brep: &BRep) -> Vec<&Surface3> {
    let mut surfaces = Vec::new();
    for ts in &brep.tshapes {
        if let TShape::Face(fd) = ts.as_ref() {
            if let Some(ref surf) = fd.surface {
                surfaces.push(surf);
            }
        }
    }
    surfaces
}

/// Iterate face ShapeRefs for all shells in the first solid of the BRep.
fn first_solid_face_refs(brep: &BRep) -> Vec<ShapeRef> {
    let mut face_refs = Vec::new();
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            for sh_ref in &sd.shells {
                if let TShape::Shell(ref sh_data) = *brep.tshapes[sh_ref.index] {
                    face_refs.extend(sh_data.faces.iter().copied());
                }
            }
            break; // Only first solid
        }
    }
    face_refs
}

/// Extract box geometry from a rcad_kernel::BRep.
///
/// A box is recognized as a solid with 6 planar faces arranged in 3 pairs
/// of parallel faces with appropriate normals.
///
/// Returns `None` if the shape is not a valid box.
pub fn get_box_geometry(brep: &rcad_kernel::BRep) -> Option<BoxGeometry> {
    // Must have exactly one solid
    if brep.solid_count() != 1 {
        return None;
    }
    // Must have exactly one shell
    let shell_count = brep
        .tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Shell(_)))
        .count();
    if shell_count != 1 {
        return None;
    }

    // A box has exactly 6 faces
    if brep.face_count() != 6 {
        return None;
    }

    // Get all face surfaces
    let surfaces = collect_face_surfaces(brep);
    if surfaces.len() != 6 {
        return None;
    }

    // All surfaces must be planes
    let planes: Vec<&Plane> = surfaces
        .iter()
        .filter_map(|s| match s {
            Surface3::Plane(p) => Some(p),
            _ => None,
        })
        .collect();

    if planes.len() != 6 {
        return None;
    }

    // Group planes by normal direction (allowing for opposite directions)
    let mut normal_groups: Vec<(DVec3, Vec<&Plane>)> = Vec::new();

    for plane in &planes {
        let normal = plane.normal.normalize_or_zero();
        let found = normal_groups.iter_mut().find(|(n, _)| {
            let dot = n.dot(normal).abs();
            dot > 0.999
        });

        if let Some((_, group)) = found {
            group.push(*plane);
        } else {
            normal_groups.push((normal, vec![*plane]));
        }
    }

    // A box has 3 pairs of parallel faces
    if normal_groups.len() != 3 {
        return None;
    }

    for (_, group) in &normal_groups {
        if group.len() != 2 {
            return None;
        }
    }

    // Compute box dimensions and origin
    let bbox = brep.bounding_box()?;
    let min_pt = bbox[0];
    let max_pt = bbox[1];

    let dx = (max_pt.x - min_pt.x).abs();
    let dy = (max_pt.y - min_pt.y).abs();
    let dz = (max_pt.z - min_pt.z).abs();

    // Verify the planes correspond to the bounding box
    let tolerance = TOLERANCE_MESH_LEGACY;
    for plane in &planes {
        let d = plane.normal.dot(plane.origin);
        let d_min = plane.normal.dot(min_pt);
        let d_max = plane.normal.dot(max_pt);

        let near_min = (d - d_min).abs() < tolerance;
        let near_max = (d - d_max).abs() < tolerance;

        if !near_min && !near_max {
            return None;
        }
    }

    Some(BoxGeometry {
        origin: min_pt,
        dx,
        dy,
        dz,
    })
}

/// Extract cylinder geometry from a rcad_kernel::BRep.
///
/// A cylinder is recognized as a solid with:
/// - One cylindrical lateral face
/// - Two planar end caps (optional for partial cylinders)
///
/// Returns `None` if the shape is not a valid cylinder.
pub fn get_cylinder_geometry(brep: &rcad_kernel::BRep) -> Option<CylinderGeometry> {
    // Must have exactly one solid
    if brep.solid_count() != 1 {
        return None;
    }
    // Must have exactly one shell
    let shell_count = brep
        .tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Shell(_)))
        .count();
    if shell_count != 1 {
        return None;
    }

    // Get all face surfaces
    let surfaces = collect_face_surfaces(brep);

    // Find the cylindrical surface
    let mut cyl_surf: Option<&CylindricalSurface> = None;
    for surf in &surfaces {
        if let Surface3::Cylinder(c) = surf {
            if cyl_surf.is_some() {
                return None;
            }
            cyl_surf = Some(c);
        }
    }

    let cyl = cyl_surf?;

    // Check that other surfaces are planes (caps)
    for surf in &surfaces {
        match surf {
            Surface3::Cylinder(_) => {}
            Surface3::Plane(_) => {}
            _ => return None,
        }
    }

    // Compute height from the bounding box along the cylinder axis
    let axis = cyl.axis.normalize_or_zero();
    let bbox = brep.bounding_box()?;
    let min_pt = bbox[0];
    let max_pt = bbox[1];

    let mut min_proj = f64::INFINITY;
    let mut max_proj = f64::NEG_INFINITY;

    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let corner = DVec3::new(
                    if i == 0 { min_pt.x } else { max_pt.x },
                    if j == 0 { min_pt.y } else { max_pt.y },
                    if k == 0 { min_pt.z } else { max_pt.z },
                );
                let proj = corner.dot(axis);
                min_proj = min_proj.min(proj);
                max_proj = max_proj.max(proj);
            }
        }
    }

    let height = (max_proj - min_proj).abs();

    // Determine the actual bottom of the cylinder
    let origin_proj = cyl.origin.dot(axis);
    let origin = if (origin_proj - min_proj).abs() < TOLERANCE_MESH_LEGACY {
        cyl.origin
    } else {
        cyl.origin + axis * (min_proj - origin_proj)
    };

    Some(CylinderGeometry {
        origin,
        axis,
        radius: cyl.radius,
        height,
    })
}

/// Extract sphere geometry from a rcad_kernel::BRep.
///
/// A sphere is recognized as a solid with a single spherical face.
///
/// Returns `None` if the shape is not a valid sphere.
pub fn get_sphere_geometry(brep: &rcad_kernel::BRep) -> Option<SphereGeometry> {
    // Must have exactly one solid
    if brep.solid_count() != 1 {
        return None;
    }

    // Get all face surfaces
    let surfaces = collect_face_surfaces(brep);

    // For a sphere, we expect exactly one spherical surface
    let mut sphere_surf: Option<&SphericalSurface> = None;
    for surf in &surfaces {
        match surf {
            Surface3::Sphere(s) => {
                if sphere_surf.is_some() {
                    return None;
                }
                sphere_surf = Some(s);
            }
            _ => return None,
        }
    }

    let sphere = sphere_surf?;

    Some(SphereGeometry {
        center: sphere.center,
        radius: sphere.radius,
    })
}

/// Extract cone geometry from a rcad_kernel::BRep.
///
/// A cone is recognized as a solid with:
/// - One conical lateral face
/// - Optionally one planar base cap
///
/// Returns `None` if the shape is not a valid cone.
pub fn get_cone_geometry(brep: &rcad_kernel::BRep) -> Option<ConeGeometry> {
    // Must have exactly one solid
    if brep.solid_count() != 1 {
        return None;
    }

    // Get all face surfaces
    let surfaces = collect_face_surfaces(brep);

    // Find the conical surface
    let mut cone_surf: Option<&ConicalSurface> = None;
    for surf in &surfaces {
        if let Surface3::Cone(c) = surf {
            if cone_surf.is_some() {
                return None;
            }
            cone_surf = Some(c);
        }
    }

    let cone = cone_surf?;

    // Check that other surfaces are planes (caps)
    for surf in &surfaces {
        match surf {
            Surface3::Cone(_) => {}
            Surface3::Plane(_) => {}
            _ => return None,
        }
    }

    // Get the apex and angle from the cone
    let apex = cone.apex_point();
    let axis = cone.axis.normalize_or_zero();
    let angle = cone.half_angle_rad;

    Some(ConeGeometry { apex, axis, angle })
}

/// Extract torus geometry from a rcad_kernel::BRep.
///
/// A torus is recognized as a solid with a single toroidal face.
///
/// Returns `None` if the shape is not a valid torus.
pub fn get_torus_geometry(brep: &rcad_kernel::BRep) -> Option<TorusGeometry> {
    // Must have exactly one solid
    if brep.solid_count() != 1 {
        return None;
    }

    // Get all face surfaces
    let surfaces = collect_face_surfaces(brep);

    // For a torus, we expect exactly one toroidal surface
    let mut torus_surf: Option<&ToroidalSurface> = None;
    for surf in &surfaces {
        match surf {
            Surface3::Torus(t) => {
                if torus_surf.is_some() {
                    return None;
                }
                torus_surf = Some(t);
            }
            _ => return None,
        }
    }

    let torus = torus_surf?;

    Some(TorusGeometry {
        center: torus.center,
        axis: torus.axis.normalize_or_zero(),
        major_radius: torus.major_radius,
        minor_radius: torus.minor_radius,
    })
}

// =============================================================================
// Primitive Detection Functions
// =============================================================================

/// Check if a rcad_kernel::BRep represents a box.
pub fn is_box(brep: &rcad_kernel::BRep) -> bool {
    get_box_geometry(brep).is_some()
}

/// Check if a rcad_kernel::BRep represents a cylinder.
pub fn is_cylinder(brep: &rcad_kernel::BRep) -> bool {
    get_cylinder_geometry(brep).is_some()
}

/// Check if a rcad_kernel::BRep represents a sphere.
pub fn is_sphere(brep: &rcad_kernel::BRep) -> bool {
    get_sphere_geometry(brep).is_some()
}

/// Check if a rcad_kernel::BRep represents a cone.
pub fn is_cone(brep: &rcad_kernel::BRep) -> bool {
    get_cone_geometry(brep).is_some()
}

/// Check if a rcad_kernel::BRep represents a torus.
pub fn is_torus(brep: &rcad_kernel::BRep) -> bool {
    get_torus_geometry(brep).is_some()
}

// =============================================================================
// Tests
// =============================================================================
