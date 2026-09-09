//! OCCT ShapeBuild_Vertex (TKShHealing ShapeBuild package).
//!
//! 1:1 translation of `ShapeBuild_Vertex.hxx` L27-53 and
//! `ShapeBuild_Vertex.cxx` L25-72. Provides low-level functions used for
//! constructing vertices.
//!
//! Architecture note: OCCT's methods are `const` member functions reading the
//! vertices through `BRep_Tool` and building the result through a local
//! `BRep_Builder`; rcad TShapes live in the [`BRep`] pool, so the builder
//! target pool is threaded as the leading `&mut BRep` argument (the
//! shape_build module precedents).

use glam::DVec3;
use rcad_kernel::topo::topods::{BRep, BRepTool, Shape};

/// OCCT ShapeBuild_Vertex (`ShapeBuild_Vertex.hxx` L27-53).
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeBuildVertex;

impl ShapeBuildVertex {
    /// OCCT ShapeBuild_Vertex::CombineVertex(V1, V2, tolFactor = 1.0001)
    /// (ShapeBuild_Vertex.cxx L25-34): combines a new vertex from two others.
    /// This new one is the smallest vertex which comprises both of the source
    /// vertices. The function takes into account the positions and tolerances
    /// of the source vertices. The tolerance of the new vertex will be equal
    /// to the minimal tolerance that is required to comprise source vertices
    /// multiplied by tolFactor (in order to avoid errors because of
    /// discreteness of calculations).
    pub fn combine_vertex(&self, brep: &mut BRep, v1: &Shape, v2: &Shape, tol_factor: f64) -> Shape {
        // OCCT delegates with BRep_Tool::Pnt / BRep_Tool::Tolerance of both
        // vertices (world coordinates).
        let (pnt1, tol1) = {
            let p = brep.vertex_position(v1);
            let t = brep.vertex_tolerance(v1);
            (p, t)
        };
        let (pnt2, tol2) = {
            let p = brep.vertex_position(v2);
            let t = brep.vertex_tolerance(v2);
            (p, t)
        };
        self.combine_vertex_points(brep, pnt1, pnt2, tol1, tol2, tol_factor)
    }

    /// OCCT ShapeBuild_Vertex::CombineVertex(pnt1, pnt2, tol1, tol2,
    /// tolFactor = 1.0001) (ShapeBuild_Vertex.cxx L38-72): the same function
    /// as above, except that it accepts two points and two tolerances instead
    /// of vertices.
    pub fn combine_vertex_points(
        &self,
        brep: &mut BRep,
        pnt1: DVec3,
        pnt2: DVec3,
        tol1: f64,
        tol2: f64,
        tol_factor: f64,
    ) -> Shape {
        // gp_Vec v = pnt2.XYZ() - pnt1.XYZ(); dist = v.Magnitude().
        let v = pnt2 - pnt1;
        let dist = v.length();

        let pos;
        let tol;

        // #47 rln 09.12.98 S4054 PRO14323 entity 2844.
        if dist + tol2 <= tol1 {
            pos = pnt1;
            tol = tol1;
        } else if dist + tol1 <= tol2 {
            pos = pnt2;
            tol = tol2;
        } else {
            tol = 0.5 * (dist + tol1 + tol2);
            // szv#4:S4163:12Mar99 anti-exception.
            let s = if dist > 0.0 { (tol2 - tol1) / dist } else { 0.0 };
            pos = 0.5 * ((1.0 - s) * pnt1 + (1.0 + s) * pnt2);
        }

        // OCCT BRep_Builder::MakeVertex(V, pos, tolFactor * tol): a fresh
        // TVertex TShape; rcad's add_tvertex_unique is that fresh-TShape
        // constructor (add_tvertex would share through the position identity
        // cache).
        let vv = brep.add_tvertex_unique(pos);
        brep.vertex_mut(vv.clone()).tolerance = tol_factor * tol;
        vv
    }
}
