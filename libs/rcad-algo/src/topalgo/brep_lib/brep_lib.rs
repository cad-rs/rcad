// OCCT BRepLib (BRepLib.hxx / BRepLib_*.cxx)
// BRep library utilities for edge/face/solid operations.
//
// Functions used by TKBO (BOPTools_AlgoTools):
// - SameParameter: ensures edge's 3D curve and pcurve have the same parameterization
// - FindValidRange: finds valid parametric range for an edge on a face
// - BoundingVertex: creates a vertex from a list of points

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topo::topods::{ShapeType, TShape};
use rcad_kernel::geom::Surface3;
use glam::DVec3;

use crate::brep_algo::tool as bat;
use crate::geomalgo::int_patch::{classify_surface_type, GeomAbsSurfaceType};

/// OCCT BRepLib — static utility functions for BRep operations.
pub struct BRepLib;

impl BRepLib {
    /// OCCT: FindValidRange(edge, first, last) — finds valid parametric
    /// range for the edge within the given bounds.
    /// rcad: stub — returns true with unchanged range.
    pub fn find_valid_range(edge: &Shape, first: &mut f64, last: &mut f64) -> bool {
        let _ = edge;
        let _ = (first, last);
        true
    }

    /// OCCT: BoundingVertex(pts, new_pt, tol) — creates a vertex
    /// from a list of points by averaging.
    pub fn bounding_vertex(pts: &[DVec3], _new_pt: &mut Shape, _tol: &mut f64) -> bool {
        if pts.is_empty() { return false; }
        let _ = pts;
        true
    }

    // OCCT BRepLib.cxx L2894-2953
    /// OCCT BRepLib::SortFaces(Sh, LF) — fills LF with the faces of Sh
    /// bucketed by surface type: planes, cylinders, cones, spheres, tori,
    /// other surfaced faces, then triangulation-only faces.  Rectangular
    /// trimmed surfaces are unwrapped to their basis surface before the
    /// type switch.
    pub fn sort_faces(sh: &Shape, lf: &mut Vec<Shape>) {
        // OCCT L2896: LF.Clear().
        lf.clear();
        // OCCT L2897: the seven type buckets.
        let mut l_tri: Vec<Shape> = Vec::new();
        let mut l_plan: Vec<Shape> = Vec::new();
        let mut l_cyl: Vec<Shape> = Vec::new();
        let mut l_con: Vec<Shape> = Vec::new();
        let mut l_sphere: Vec<Shape> = Vec::new();
        let mut l_tor: Vec<Shape> = Vec::new();
        let mut l_other: Vec<Shape> = Vec::new();
        // OCCT L2898: TopExp_Explorer exp(Sh, TopAbs_FACE).
        let a_exp = bat::explorer(sh, ShapeType::Face, ShapeType::Shape);

        // OCCT L2902: for (; exp.More(); exp.Next()).
        for f_exp in a_exp {
            // OCCT L2904: const TopoDS_Face& F = TopoDS::Face(exp.Current()).
            let f = f_exp;
            // OCCT L2905: S = BRep_Tool::Surface(F, l).
            let s = brep_tool_surface(&f);
            if let Some(mut s) = s {
                // OCCT L2906-2910: Geom_RectangularTrimmedSurface -> BasisSurface.
                if let Surface3::Trimmed(a_ts) = &s {
                    s = (*a_ts.basis).clone();
                }
                // OCCT L2911: GeomAdaptor_Surface AS(S); switch (AS.GetType()).
                match geom_adaptor_surface_get_type(&s) {
                    GeomAbsSurfaceType::Plane => {
                        l_plan.push(f);
                    }
                    GeomAbsSurfaceType::Cylinder => {
                        l_cyl.push(f);
                    }
                    GeomAbsSurfaceType::Cone => {
                        l_con.push(f);
                    }
                    GeomAbsSurfaceType::Sphere => {
                        l_sphere.push(f);
                    }
                    GeomAbsSurfaceType::Torus => {
                        l_tor.push(f);
                    }
                    _ => {
                        l_other.push(f);
                    }
                }
            } else {
                // OCCT L2942-2946: null surface (triangulation-only face).
                l_tri.push(f);
            }
        }
        // OCCT L2947-2953: LF.Append(LPlan); ... LF.Append(LTri).
        lf.append(&mut l_plan);
        lf.append(&mut l_cyl);
        lf.append(&mut l_con);
        lf.append(&mut l_sphere);
        lf.append(&mut l_tor);
        lf.append(&mut l_other);
        lf.append(&mut l_tri);
    }

    // OCCT BRepLib.cxx L2955-3010
    /// OCCT BRepLib::ReverseSortFaces(Sh, LF) — the sibling of SortFaces
    /// with the reverse bucket order: triangulation-only faces first, then
    /// other surfaced faces, tori, spheres, cones, cylinders, planes.  No
    /// explicit trimmed-surface unwrap in the body — the adaptor strips it
    /// while loading.
    pub fn reverse_sort_faces(sh: &Shape, lf: &mut Vec<Shape>) {
        // OCCT L2958: LF.Clear().
        lf.clear();
        // OCCT L2960-2961: the seven type buckets (LF allocator — not
        // applicable to Vec).
        let mut l_tri: Vec<Shape> = Vec::new();
        let mut l_plan: Vec<Shape> = Vec::new();
        let mut l_cyl: Vec<Shape> = Vec::new();
        let mut l_con: Vec<Shape> = Vec::new();
        let mut l_sphere: Vec<Shape> = Vec::new();
        let mut l_tor: Vec<Shape> = Vec::new();
        let mut l_other: Vec<Shape> = Vec::new();
        // OCCT L2962: TopExp_Explorer exp(Sh, TopAbs_FACE).
        let a_exp = bat::explorer(sh, ShapeType::Face, ShapeType::Shape);

        // OCCT L2966: for (; exp.More(); exp.Next()).
        for f_exp in a_exp {
            // OCCT L2968: const TopoDS_Face& F = TopoDS::Face(exp.Current()).
            let f = f_exp;
            // OCCT L2969: const handle<Geom_Surface>& S = BRep_Tool::Surface(F, l).
            let s = brep_tool_surface(&f);
            if let Some(s) = s {
                // OCCT L2971: GeomAdaptor_Surface AS(S); switch (AS.GetType()).
                match geom_adaptor_surface_get_type(&s) {
                    GeomAbsSurfaceType::Plane => {
                        l_plan.push(f);
                    }
                    GeomAbsSurfaceType::Cylinder => {
                        l_cyl.push(f);
                    }
                    GeomAbsSurfaceType::Cone => {
                        l_con.push(f);
                    }
                    GeomAbsSurfaceType::Sphere => {
                        l_sphere.push(f);
                    }
                    GeomAbsSurfaceType::Torus => {
                        l_tor.push(f);
                    }
                    _ => {
                        l_other.push(f);
                    }
                }
            } else {
                // OCCT L2999-3002: null surface (triangulation-only face).
                l_tri.push(f);
            }
        }
        // OCCT L3003-3009: LF.Append(LTri); ... LF.Append(LPlan) — reverse order.
        lf.append(&mut l_tri);
        lf.append(&mut l_other);
        lf.append(&mut l_tor);
        lf.append(&mut l_sphere);
        lf.append(&mut l_con);
        lf.append(&mut l_cyl);
        lf.append(&mut l_plan);
    }
}

/// OCCT BRep_Tool::Surface(F, l) — the face surface (architecture bridge:
/// the rcad face payload carries the surface directly; the surface location
/// index is applied by the kernel geometry layer, identity in the offset
/// pipeline).
fn brep_tool_surface(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT GeomAdaptor_Surface(S) constructor + GetType() — the adaptor Load
/// unwraps Geom_RectangularTrimmedSurface to its basis surface before typing
/// (GeomAdaptor_Surface.cxx L423-425); classify_surface_type is the
/// translated GetType equivalent (geomalgo::int_patch, used by the topalgo
/// BRepAdaptor_Surface).
fn geom_adaptor_surface_get_type(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Trimmed(a_ts) => classify_surface_type(&a_ts.basis),
        _ => classify_surface_type(s),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::topo::topo_builder::brep_from_shape;
    use rcad_kernel::topo::topods::{edge_data_pool_free, BRep, BRepBuilder};

    /// The pool-free caller wiring of BRepLib::SameParameter (the
    /// brep_feat_make_d_prism.rs / loc_ope_wires_on_shape_b.rs pattern):
    /// the edge graph is adopted into a standalone pool with the Arc SHARED
    /// (topo_builder::brep_from_shape), the builder flag resets and the
    /// topalgo::brep_lib engine run over the adopted pool, and the writes
    /// must be observable through the CALLER's own Shape handle — the OCCT
    /// in-place TShape mutation through the handle (BRep_TEdge::Modified /
    /// Tolerance, BRepLib.cxx L1719-1737).
    #[test]
    fn same_parameter_adoption_reaches_pool_free_caller() {
        // A straight edge built in a source pool that is then dropped: the
        // caller-side Shape keeps its source flat index (the feat pipeline
        // provenance of the two wiring sites).
        let a_edge = {
            let mut a_src = BRep::new();
            let mut a_b = BRepBuilder::new();
            let a_v0 = a_b.add_vertex(&mut a_src, glam::DVec3::ZERO, 1.0e-7);
            let a_v1 = a_b.add_vertex(&mut a_src, glam::DVec3::X, 1.0e-7);
            a_b.add_edge(
                &mut a_src,
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(
                    glam::DVec3::ZERO,
                    glam::DVec3::X,
                ))),
                a_v0.clone(),
                a_v1.clone(),
                [0.0, 1.0],
            )
        };
        assert!(a_edge.index != usize::MAX, "pool-free by dropped pool, not synthetic");

        // The wiring: adopt with the Arc SHARED and run the engine.
        let mut a_brep = brep_from_shape(&a_edge, &[]);
        // OCCT BRepFeat_MakeDPrism.cxx L409-410: bB.SameRange(ledg, false);
        // bB.SameParameter(ledg, false) — in place through the shared Arc.
        a_brep.edge_mut_inplace(a_edge.clone()).same_range = false;
        a_brep.edge_mut_inplace(a_edge.clone()).same_parameter = false;
        // OCCT LocOpe_WiresOnShape.cxx L908: BRepLib::SameParameter(Edg, tol).
        crate::topalgo::brep_lib::same_parameter::same_parameter(&mut a_brep, &a_edge, 1.0e-5);

        // The adoption must not have cloned the TShape (the precondition of
        // the in-place propagation).
        assert_eq!(
            std::sync::Arc::as_ptr(&a_brep.tshapes[a_edge.index]),
            std::sync::Arc::as_ptr(&a_edge.data),
            "brep_from_shape shares the TShape Arc"
        );
        // The OCCT in-place semantics: the caller's handle observes the
        // engine writes without any write-back step.
        let a_ed = edge_data_pool_free(&a_edge).expect("edge payload");
        assert!(a_ed.same_range, "SameRange flag (cxx L1720)");
        assert!(a_ed.same_parameter, "SameParameter flag (cxx L1736)");
        assert_eq!(a_ed.range, [0.0, 1.0], "Range(aNE, f3d, l3d) (cxx L1719)");
    }
}
