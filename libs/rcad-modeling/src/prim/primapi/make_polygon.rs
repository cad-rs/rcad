//! OCCT BRepBuilderAPI_MakePolygon (TKTopAlgo) — a polygon wire builder.
//!
//! OCCT BRepBuilderAPI_MakePolygon accumulates points via Add() and closes
//! the loop with Close() (BRepBuilderAPI_MakePolygon.cxx); the result is a
//! wire of straight edges.  The rcad port keeps the method surface
//! (new / add / close); each edge follows the BRepLib_MakeEdge::Init vertex
//! orientation contract (first FORWARD, second REVERSED — see
//! make_planar_rect_brep for the rationale).

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{BRep, Orientation};

use crate::BuildError;

/// OCCT BRepBuilderAPI_MakePolygon — point accumulator producing a wire.
#[derive(Debug, Default)]
pub struct MakePolygon {
    points: Vec<DVec3>,
}

impl MakePolygon {
    /// OCCT BRepBuilderAPI_MakePolygon().
    pub fn new() -> Self {
        MakePolygon { points: Vec::new() }
    }

    /// OCCT BRepBuilderAPI_MakePolygon::Add(gp_Pnt).
    pub fn add(&mut self, p: DVec3) {
        self.points.push(p);
    }

    /// OCCT BRepBuilderAPI_MakePolygon::Close() — materialize the closed
    /// polygon wire into `brep` (the flat BRep pool owns the shapes).
    pub fn close(&mut self, brep: &mut BRep) -> Result<Shape, BuildError> {
        let n = self.points.len();
        if n < 3 {
            return Err(BuildError::DegenerateGeometry(
                "polygon needs at least 3 points",
            ));
        }
        // BRepLib_MakeEdge::Init vertex orientation contract: the first
        // endpoint is stored FORWARD, the second REVERSED.
        let rev = |sr: Shape| Shape {
            orientation: Orientation::Reversed,
            ..sr
        };
        let mut vertices = Vec::with_capacity(n);
        for &p in &self.points {
            vertices.push(brep.add_tvertex(p));
        }
        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            let j = (i + 1) % n;
            let delta = self.points[j] - self.points[i];
            let len = delta.length();
            if len < 1e-12 {
                return Err(BuildError::DegenerateGeometry(
                    "zero-length polygon edge",
                ));
            }
            let e = brep.add_tedge(
                Some(Curve3::Line(Line3::new(self.points[i], delta / len))),
                vertices[i].clone(),
                rev(vertices[j].clone()),
                [0.0, len],
            );
            edges.push(e);
        }
        Ok(brep.add_twire(edges))
    }
}
