//! Convert a solved [`Sketch`] to a wire [`BRep`] on the XY plane (z = 0).
//!
//! Each entity becomes one or more edges in the BRep.  No faces are created;
//! the result is suitable for visualisation or as input to [`extrude`] after
//! the caller builds a proper closed wire.

use glam::DVec3;
use rcad_kernel::{BRep, Edge, Vertex};
use rcad_kernel::geom::{Circle3, Curve3, Line3};

use crate::entity::EntityKind;
use crate::sketch::Sketch;

impl Sketch {
    /// Build a wire [`BRep`] from all entities in the sketch.
    ///
    /// Each entity is converted to edges on the XY plane (z = 0):
    ///
    /// | Entity | Result |
    /// |--------|--------|
    /// | Point  | A single vertex (no edge) |
    /// | Line   | One edge with a `Line3` curve |
    /// | Circle | One closed edge with a `Circle3` curve |
    /// | Arc    | One edge with a `Circle3` curve, parameter range [start, end] |
    pub fn to_wire_brep(&self) -> BRep {
        let mut brep = BRep::new();

        for (_id, entity) in self.entities.iter().enumerate() {
            let p = &self.params[entity.param_start..entity.param_start + entity.kind.param_count()];

            match entity.kind {
                EntityKind::Point => {
                    // Just a vertex, no edge.
                    brep.vertices.push(Vertex { point: DVec3::new(p[0], p[1], 0.0) });
                    // Pad geom store vectors to keep indices consistent.
                    brep.geom.edge_curve.push(None);
                    brep.geom.edge_curve_range.push(None);
                    brep.geom.edge_degenerated.push(false);
                    brep.geom.edge_tolerance.push(1e-7);
                    brep.geom.edge_same_parameter.push(true);
                    brep.geom.edge_same_range.push(true);
                    brep.geom.edge_pcurves.push(Vec::new());
                }

                EntityKind::Line => {
                    let (x1, y1, x2, y2) = (p[0], p[1], p[2], p[3]);
                    let start = DVec3::new(x1, y1, 0.0);
                    let end = DVec3::new(x2, y2, 0.0);
                    let dir = (end - start).normalize_or_zero();
                    let length = (end - start).length();

                    let v_start = brep.vertices.len();
                    brep.vertices.push(Vertex { point: start });
                    brep.vertices.push(Vertex { point: end });

                    let curve_idx = brep.geom.curves.len();
                    brep.geom.curves.push(Curve3::Line(Line3 { origin: start, direction: dir }));

                    let edge_idx = brep.edges.len();
                    brep.edges.push(Edge { start: v_start, end: v_start + 1 });

                    // Extend geom store parallel arrays.
                    pad_geom_to(&mut brep, edge_idx);
                    brep.geom.edge_curve[edge_idx] = Some(curve_idx);
                    brep.geom.edge_curve_range[edge_idx] = Some([0.0, length]);
                }

                EntityKind::Circle => {
                    let (cx, cy, r) = (p[0], p[1], p[2]);
                    let center = DVec3::new(cx, cy, 0.0);

                    // A closed circle: start == end vertex.
                    let v_idx = brep.vertices.len();
                    brep.vertices.push(Vertex { point: center + DVec3::new(r, 0.0, 0.0) });

                    let curve_idx = brep.geom.curves.len();
                    brep.geom.curves.push(Curve3::Circle(Circle3 {
                        center,
                        normal: DVec3::Z,
                        radius: r,
                    }));

                    let edge_idx = brep.edges.len();
                    brep.edges.push(Edge { start: v_idx, end: v_idx }); // closed

                    pad_geom_to(&mut brep, edge_idx);
                    brep.geom.edge_curve[edge_idx] = Some(curve_idx);
                    brep.geom.edge_curve_range[edge_idx] =
                        Some([0.0, 2.0 * std::f64::consts::PI]);
                }

                EntityKind::Arc => {
                    let (cx, cy, r, t0, t1) = (p[0], p[1], p[2], p[3], p[4]);
                    let center = DVec3::new(cx, cy, 0.0);
                    let start_pt = center + DVec3::new(r * t0.cos(), r * t0.sin(), 0.0);
                    let end_pt = center + DVec3::new(r * t1.cos(), r * t1.sin(), 0.0);

                    let v_start = brep.vertices.len();
                    brep.vertices.push(Vertex { point: start_pt });
                    brep.vertices.push(Vertex { point: end_pt });

                    let curve_idx = brep.geom.curves.len();
                    brep.geom.curves.push(Curve3::Circle(Circle3 {
                        center,
                        normal: DVec3::Z,
                        radius: r,
                    }));

                    let edge_idx = brep.edges.len();
                    brep.edges.push(Edge { start: v_start, end: v_start + 1 });

                    pad_geom_to(&mut brep, edge_idx);
                    brep.geom.edge_curve[edge_idx] = Some(curve_idx);
                    brep.geom.edge_curve_range[edge_idx] = Some([t0, t1]);
                }
            }
        }

        brep
    }
}

/// Ensure all parallel geom-store arrays are at least `edge_idx + 1` long,
/// filling with defaults.
fn pad_geom_to(brep: &mut BRep, edge_idx: usize) {
    let target = edge_idx + 1;
    while brep.geom.edge_curve.len() < target {
        brep.geom.edge_curve.push(None);
    }
    while brep.geom.edge_curve_range.len() < target {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.edge_degenerated.len() < target {
        brep.geom.edge_degenerated.push(false);
    }
    while brep.geom.edge_tolerance.len() < target {
        brep.geom.edge_tolerance.push(1e-7);
    }
    while brep.geom.edge_same_parameter.len() < target {
        brep.geom.edge_same_parameter.push(true);
    }
    while brep.geom.edge_same_range.len() < target {
        brep.geom.edge_same_range.push(true);
    }
    while brep.geom.edge_pcurves.len() < target {
        brep.geom.edge_pcurves.push(Vec::new());
    }
}
