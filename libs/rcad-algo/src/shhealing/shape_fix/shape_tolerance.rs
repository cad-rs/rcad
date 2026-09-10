//! 1:1 translation of OCCT `ShapeFix_ShapeTolerance`
//! (`TKShHealing/ShapeFix/ShapeFix_ShapeTolerance.hxx` L17-63 +
//! `ShapeFix_ShapeTolerance.cxx` L1-213, docket row `ShapeFix_ShapeTolerance`).
//!
//! Function-count equation (OCCT ShapeFix_ShapeTolerance.cxx = rcad
//! shape_tolerance.rs): `ShapeFix_ShapeTolerance()` / `LimitTolerance` /
//! `SetTolerance` — 3 OCCT functions = 3 rcad functions.
//!
//! Architecture bridges:
//! - `TopoDS_Shape.TShape()->Tolerance(...)` direct writes — the rcad
//!   `BRep::vertex_mut/edge_mut/face_mut` pool accessors (the same
//!   TShape-handle write; `BRepBuilder::update_*` variants apply Max
//!   semantics, which OCCT's direct `TV->Tolerance()` writes do not).
//! - `BRep_Tool::Tolerance` — the `BRepTool::tolerance` accessor.
//! - `TopExp_Explorer` / `TopExp::Vertices` — the shared re-hosts in
//!   `shape_build/brep_tool.rs` and the `TopExp.cxx` re-host in the package
//!   statics (`shape_fix.rs`).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, ShapeType};

use crate::shhealing::shape_build::brep_tool::topexp_explorer;
use crate::shhealing::shape_fix::shape_fix::top_exp_vertices;

/// OCCT ShapeFix_ShapeTolerance (ShapeFix_ShapeTolerance.hxx L34-63): tool
/// for setting or limiting the tolerances of the sub-shapes of a shape.
///
/// OCCT declares no members — the class is a stateless tool (cxx L31: the
/// defaulted constructor).
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeFixShapeTolerance;

impl ShapeFixShapeTolerance {
    /// OCCT ShapeFix_ShapeTolerance::ShapeFix_ShapeTolerance() (cxx L31):
    /// the defaulted empty constructor.
    pub fn new() -> Self {
        ShapeFixShapeTolerance
    }

    /// OCCT ShapeFix_ShapeTolerance::LimitTolerance (cxx L35-138): limits
    /// the tolerances of the sub-shapes of the given type to the range
    /// [tmin, tmax] (tmin < 0 or a null shape does nothing).  WIRE limits
    /// its edges and their vertices; every other type recurses over
    /// VERTEX + EDGE + FACE.  Returns true when at least one tolerance was
    /// changed.
    pub fn limit_tolerance(
        &self,
        brep: &mut BRep,
        shape: &Shape,
        tmin: f64,
        tmax: f64,
        styp: ShapeType,
    ) -> bool {
        // L40-43: if (shape.IsNull() || tmin < 0) return false.
        if shape.is_null() || tmin < 0.0 {
            return false;
        }
        // L44-46.
        let iamax = tmax >= tmin;
        let mut prec = 0.0f64;
        let mut fait = false;
        if styp == ShapeType::Vertex || styp == ShapeType::Edge || styp == ShapeType::Face {
            // L49-110: the direct-type walk.
            for ex in topexp_explorer(brep, shape, styp) {
                let sh = ex;
                let mut newtol = 0i32;
                if styp == ShapeType::Vertex {
                    // L55-70.
                    let v = sh;
                    prec = brep.tolerance(&v);
                    if iamax && prec > tmax {
                        newtol = 1;
                    } else if prec < tmin {
                        newtol = -1;
                    }
                    if newtol != 0 {
                        // L67-68: TV->Tolerance(newtol > 0 ? tmax : tmin).
                        brep.vertex_mut(v).tolerance = if newtol > 0 { tmax } else { tmin };
                        fait = true;
                    }
                } else if styp == ShapeType::Edge {
                    // L74-89.
                    let e = sh;
                    prec = brep.tolerance(&e);
                    if iamax && prec > tmax {
                        newtol = 1;
                    } else if prec < tmin {
                        newtol = -1;
                    }
                    if newtol != 0 {
                        // L86-87: TE->Tolerance(newtol > 0 ? tmax : tmin).
                        brep.edge_mut(e).tolerance = if newtol > 0 { tmax } else { tmin };
                        fait = true;
                    }
                } else if styp == ShapeType::Face {
                    // L93-108.
                    let f = sh;
                    prec = brep.tolerance(&f);
                    if iamax && prec > tmax {
                        newtol = 1;
                    } else if prec < tmin {
                        newtol = -1;
                    }
                    if newtol != 0 {
                        // L105-106: TF->Tolerance(newtol > 0 ? tmax : tmin).
                        brep.face_mut(f).tolerance = if newtol > 0 { tmax } else { tmin };
                        fait = true;
                    }
                }
            }
        } else if styp == ShapeType::Wire {
            // L112-130: the WIRE walk — edges plus their extremity vertices.
            for ex in topexp_explorer(brep, shape, ShapeType::Edge) {
                let sh = ex;
                let e = sh;
                self.limit_tolerance(brep, &e, tmin, tmax, ShapeType::Edge);
                let (v1, v2) = top_exp_vertices(brep, &e);
                if !v1.is_null() {
                    fait |= self.limit_tolerance(brep, &v1, tmin, tmax, ShapeType::Vertex);
                }
                if !v2.is_null() {
                    fait |= self.limit_tolerance(brep, &v2, tmin, tmax, ShapeType::Vertex);
                }
            }
        } else {
            // L131-136: every other type — all three direct types.
            fait |= self.limit_tolerance(brep, shape, tmin, tmax, ShapeType::Vertex);
            fait |= self.limit_tolerance(brep, shape, tmin, tmax, ShapeType::Edge);
            fait |= self.limit_tolerance(brep, shape, tmin, tmax, ShapeType::Face);
        }
        fait
    }

    /// OCCT ShapeFix_ShapeTolerance::SetTolerance (cxx L142-212): sets the
    /// tolerances of the sub-shapes of the given type to `preci` (a null
    /// shape or a non-positive value does nothing).  WIRE sets its edges and
    /// their vertices; every other type recurses over VERTEX + EDGE + FACE.
    pub fn set_tolerance(&self, brep: &mut BRep, shape: &Shape, preci: f64, styp: ShapeType) {
        //   VERTEX ou EDGE ou FACE : ces types seulement
        //   WIRE : EDGE + VERTEX
        //   Autres : TOUT (donc == WIRE + FACE)
        // L149-152: if (shape.IsNull() || preci <= 0) return.
        if shape.is_null() || preci <= 0.0 {
            return;
        }
        if styp == ShapeType::Vertex || styp == ShapeType::Edge || styp == ShapeType::Face {
            // L155-179: the direct-type walk.
            for ex in topexp_explorer(brep, shape, styp) {
                let sh = ex;
                if styp == ShapeType::Vertex {
                    // L162-163: TV->Tolerance(preci).
                    brep.vertex_mut(sh).tolerance = preci;
                } else if styp == ShapeType::Edge {
                    // L169-170: TE->Tolerance(preci).
                    brep.edge_mut(sh).tolerance = preci;
                } else if styp == ShapeType::Face {
                    // L176-177: TF->Tolerance(preci).
                    brep.face_mut(sh).tolerance = preci;
                }
            }
        } else if styp == ShapeType::Wire {
            // L181-204: the WIRE walk — edges plus their extremity vertices.
            for ex in topexp_explorer(brep, shape, ShapeType::Edge) {
                let sh = ex;
                let e = sh;
                brep.edge_mut(e.clone()).tolerance = preci;
                let (v1, v2) = top_exp_vertices(brep, &e);
                if !v1.is_null() {
                    brep.vertex_mut(v1).tolerance = preci;
                }
                if !v2.is_null() {
                    brep.vertex_mut(v2).tolerance = preci;
                }
            }
        } else {
            // L206-211: every other type — all three direct types.
            self.set_tolerance(brep, shape, preci, ShapeType::Vertex);
            self.set_tolerance(brep, shape, preci, ShapeType::Edge);
            self.set_tolerance(brep, shape, preci, ShapeType::Face);
        }
    }
}
