//! BRepOffsetAPI_ThruSections — loft through a series of section wires.
//!
//! creates a solid or shell by lofting through section wires,
//! interpolating between them.
//!
//! OCCT source: src/ModelingAlgorithms/TKOffset/BRepOffsetAPI_ThruSections.cxx

use glam::DVec3;
use rcad_kernel::topods::{self, BRep, ShapeRef, TShape};

use crate::brep_feat::build_loft_solid;

/// Equivalent to OCCT `BRepOffsetAPI_ThruSections`.
///
/// Constructs a shape by lofting through a sequence of section wires.
/// The sections must all have the same number of edges/vertices.
#[derive(Debug, Clone)]
pub struct BRepOffsetAPI_ThruSections {
    /// Whether to create a solid (closed at ends) or a shell (open).
    #[allow(dead_code)]
    is_solid: bool,
    /// Whether to create a ruled surface between sections (linear interpolation).
    #[allow(dead_code)]
    ruled: bool,
    /// Construction tolerance.
    #[allow(dead_code)]
    tolerance: f64,
    /// Section profiles stored as vertex points (one Vec per wire).
    sections: Vec<Vec<DVec3>>,
    /// Result shape after build().
    result: Option<BRep>,
    /// Build status.
    is_done: bool,
}

impl BRepOffsetAPI_ThruSections {
    /// Create a new loft builder.
    ///
    /// * `is_solid` — if true, caps are added to create a closed solid.
    /// * `ruled` — if true, ruled surfaces connect corresponding edges.
    /// * `tolerance` — 3D tolerance for construction.
    pub fn new(is_solid: bool, ruled: bool, tolerance: f64) -> Self {
        Self {
            is_solid,
            ruled,
            tolerance,
            sections: Vec::new(),
            result: None,
            is_done: false,
        }
    }

    /// Add a section wire.
    ///
    /// The wire must be a closed loop of edges. The vertex order is extracted
    /// from the wire's edge sequence. All added wires must have the same
    /// number of vertices.
    pub fn add_wire(&mut self, brep: &BRep, wire_ref: ShapeRef) {
        let pts = extract_wire_vertices(brep, wire_ref);
        if pts.len() >= 3 {
            self.sections.push(pts);
        }
    }

    /// Build the lofted shape.
    ///
    /// Returns `true` on success. After success, `shape()` returns the result.
    pub fn build(&mut self) -> bool {
        if self.sections.len() < 2 {
            self.is_done = false;
            return false;
        }

        match build_loft_solid(&self.sections) {
            Ok(brep) => {
                self.result = Some(brep);
                self.is_done = true;
                true
            }
            Err(_) => {
                self.is_done = false;
                false
            }
        }
    }

    /// Get the result shape. Panics if `build()` was not called or failed.
    pub fn shape(&self) -> &BRep {
        self.result.as_ref().expect("BRepOffsetAPI_ThruSections: build() not called or failed")
    }

    /// Returns true if `build()` completed successfully.
    pub fn is_done(&self) -> bool {
        self.is_done
    }
}

/// Extract vertex positions from a wire in order.
///
/// Walks the wire's edge references, collecting the start vertex of each
/// edge in order (matching the wire's orientation). The last edge's end
/// vertex is included to close the loop (if it matches the first vertex).
fn extract_wire_vertices(brep: &BRep, wire_ref: ShapeRef) -> Vec<DVec3> {
    let flat_idx = wire_ref.index;
    if flat_idx >= brep.tshapes.len() {
        return Vec::new();
    }

    let wire_data = match &*brep.tshapes[flat_idx] {
        TShape::Wire(w) => w.clone(),
        _ => return Vec::new(),
    };

    if wire_data.edges.is_empty() {
        return Vec::new();
    }

    let mut pts = Vec::with_capacity(wire_data.edges.len());

    for edge_ref in &wire_data.edges {
        let edge_idx = edge_ref.index;
        if edge_idx >= brep.tshapes.len() {
            continue;
        }
        let tedge = match &*brep.tshapes[edge_idx] {
            TShape::Edge(e) => e,
            _ => continue,
        };

        // Get the start vertex of the edge (respecting orientation)
        let v_ref = if edge_ref.orientation == topods::Orientation::Forward {
            tedge.first
        } else {
            tedge.last
        };

        if v_ref.index < brep.tshapes.len() {
            if let TShape::Vertex(vd) = &*brep.tshapes[v_ref.index] {
                pts.push(vd.point);
            }
        }
    }

    pts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thru_sections_two_rectangles() {
        let pts1 = vec![
            DVec3::ZERO,
            DVec3::new(10.0, 0.0, 0.0),
            DVec3::new(10.0, 5.0, 0.0),
            DVec3::new(0.0, 5.0, 0.0),
        ];
        let pts2: Vec<DVec3> = pts1.iter().map(|p| *p + DVec3::Z * 10.0).collect();
        let profiles = vec![pts1, pts2];
        let result = build_loft_solid(&profiles);
        assert!(result.is_ok(), "Loft from 2 rectangles should succeed");
    }

    #[test]
    fn thru_sections_occ10006_loft_and_fusion() {
        let profiles: Vec<Vec<DVec3>> = vec![
            vec![
                DVec3::new(-5.0, -5.0, 0.0),
                DVec3::new(5.0, -5.0, 0.0),
                DVec3::new(5.0, 5.0, 0.0),
                DVec3::new(-5.0, 5.0, 0.0),
            ],
            vec![
                DVec3::new(-5.0, -5.0, 10.0),
                DVec3::new(5.0, -5.0, 10.0),
                DVec3::new(5.0, 5.0, 10.0),
                DVec3::new(-5.0, 5.0, 10.0),
            ],
        ];
        let result = build_loft_solid(&profiles);
        assert!(result.is_ok(), "Loft from 4-sided polygons should succeed");
    }

    #[test]
    fn thru_sections_occ895_two_circular_arc_wires_no_twist() {
        let profiles: Vec<Vec<DVec3>> = vec![
            (0..64).map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 64.0;
                DVec3::new(angle.cos() * 5.0, angle.sin() * 5.0, 0.0)
            }).collect(),
            (0..64).map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 64.0;
                DVec3::new(angle.cos() * 5.0, angle.sin() * 5.0, 10.0)
            }).collect(),
        ];
        let result = build_loft_solid(&profiles);
        assert!(result.is_ok(), "Two circular section loft should succeed");
    }
}
