//! OCCT BRepBuilderAPI_Sewing.cxx — the same-parameter edge group:
//! `SameRange` (L110-144), `SameParameter(edge)` (L329-356), the sequence
//! overload `SameParameterEdge` (L358-515), the file-scope statics
//! `findNMVertices` (L517-574) and `ComputeToleranceVertex` (L576-660), and
//! the pair overload `SameParameterEdge` (L662-1188).
//!
//! Re-hosts used here (the sewing.rs conventions):
//! - `GeomLib::SameRange` -> the kernel `same_range_2d` body; its `None` is
//!   the OCCT null-handle outcome the try/catch of `SameRange` swallows.
//! - `BRepLib::SameParameter(E)` -> the crate `topalgo::brep_lib` carrier
//!   (the OCCT void overload reads the edge tolerance itself).
//! - `BRep_Builder::UpdateEdge(E, C, Tol)` -> the local
//!   `builder_update_edge_curve3d` (Arc in-place edit, pool slot contract).
//! - `BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol)` (the surface-keyed
//!   overloads) -> `update_edge_pcurve` / `update_edge_pcurve_closed` keyed
//!   by the face whose surface is S; the kernel read path resolves the
//!   sibling-face reads through the representation surface-value walk
//!   (rcad-kernel topods.rs curve_on_surface fallback = OCCT
//!   BRep_Tool.cxx L350-367).

use glam::DVec3;
use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use std::sync::Arc;

use crate::brep_algo::tool as bat;

use super::{
    iter_no_cumori, make_vertex_at, set_add_bool, top_abs_reverse,
};
use super::BRepBuilderAPISewing;

/// OCCT `BRep_Builder::MakeEdge(E)` — a pool-registered empty edge TShape
/// (architecture difference #3: rcad registers every product in the caller
/// pool; the pool-free bat::builder_make_edge would read back as null).
pub(crate) fn builder_make_edge(brep: &mut BRep) -> Shape {
    let data = Arc::new(TShape::Edge(rcad_kernel::topo::topods::TEdgeData {
        my_shapes: Vec::new(),
        flags: rcad_kernel::topo::topods::tshape_flags::DEFAULT,
        curve: None,
        first: Shape::null(),
        last: Shape::null(),
        range: [0.0, 0.0],
        degenerated: false,
        pcurves: indexmap::IndexMap::new(),
        representations: Vec::new(),
        vertex_params: std::collections::HashMap::new(),
        tolerance: 0.0,
        same_parameter: false,
        same_range: false,
    }));
    brep.tshapes.push(data.clone());
    Shape::from_parts(data, brep.tshapes.len() - 1, 0, Orientation::Forward)
}

/// OCCT `BRep_Builder::UpdateEdge(E, C, Tol)` (BRep_Builder.cxx L152-211) —
/// sets the 3D curve and raises the tolerance.
fn builder_update_edge_curve3d(brep: &mut BRep, edge: Shape, c3d: &Curve3, tol: f64) {
    let ed = brep.edge_mut_inplace(edge);
    ed.curve = Some(c3d.clone());
    ed.tolerance = ed.tolerance.max(tol);
}

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::SameRange(CurvePtr, FirstOnCurve,
    /// LastOnCurve, RequestedFirst, RequestedLast) (cxx L110-144).
    pub(crate) fn same_range(
        &self,
        curve_ptr: &Curve2d,
        first_on_curve: f64,
        last_on_curve: f64,
        requested_first: f64,
        requested_last: f64,
    ) -> Option<Curve2d> {
        // OCCT L113-135: GeomLib::SameRange(PConfusion, ...) inside
        // try/catch; a thrown Standard_Failure leaves NewCurvePtr null — the
        // kernel same_range_2d None is that same outcome.
        let new_curve_ptr = rcad_kernel::geom::same_range_2d(
            PCONFUSION,
            curve_ptr.clone(),
            first_on_curve,
            last_on_curve,
            requested_first,
            requested_last,
        );
        new_curve_ptr
    }

    /// OCCT BRepBuilderAPI_Sewing::SameParameter(edge) (cxx L329-341).
    pub(crate) fn same_parameter(&self, brep: &mut BRep, edge: &Shape) {
        // OCCT L334: BRepLib::SameParameter(edge) inside try/catch — the
        // 1-arg call binds the BRepLib.hxx L161-162 default
        // Tolerance = 1.0e-5; the topalgo::brep_lib engine is the 2-arg
        // edge overload (BRepLib.cxx L1237-1247).  The sewing edges are
        // pool-resident in `brep` (built/updated through the brep-bound
        // builder calls), so the engine's in-place TShape writes reach the
        // edge the way the OCCT engine mutates the shared TShape through
        // the handle.
        crate::topalgo::brep_lib::same_parameter::same_parameter(brep, edge, 1.0e-5);
    }

    /// OCCT BRepBuilderAPI_Sewing::SameParameterEdge(edge, seqEdges,
    /// seqForward, mapMerged, locReShape) (cxx L358-515) — merges the
    /// sequence of sections on one edge.
    pub(crate) fn same_parameter_edge_seq(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        seq_edges: &[Shape],
        seq_forward: &[bool],
        map_merged: &mut super::ShapeSet,
        loc_re_shape: &mut crate::shhealing::shape_build::reshape::ShapeBuildReShape,
    ) -> Shape {
        // Retrieve reference section
        // OCCT L365: TopoDS_Shape aTmpShape = myReShape->Apply(edge);
        let a_tmp_shape = self.my_re_shape.apply(brep, edge, ShapeType::Shape);
        // OCCT L366: TopoDS_Edge Edge1 = TopoDS::Edge(aTmpShape);
        let mut edge1 = a_tmp_shape;
        // OCCT L367-371: aTmpShape = locReShape->Apply(Edge1);
        // if (locReShape != myReShape) Edge1 = TopoDS::Edge(aTmpShape);
        let a_tmp_shape = loc_re_shape.apply(brep, &edge1, ShapeType::Shape);
        if !std::ptr::eq(loc_re_shape as *const _, &self.my_re_shape as *const _) {
            edge1 = a_tmp_shape;
        }
        // OCCT L372: bool isDone = false;
        let mut is_done = false;

        // Create data structures for temporary merged edges
        // OCCT L375-376.
        let mut list_faces1: Vec<Shape> = Vec::new();
        let mut merged_faces: super::ShapeSet = super::ShapeSet::new();

        if self.my_sewing {
            // Fill MergedFaces with faces of Edge1
            // OCCT L381: TopoDS_Shape bnd1 = edge;
            let mut bnd1 = edge.clone();
            // OCCT L382-385.
            if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd1)) {
                bnd1 = v.clone();
            }
            // OCCT L386-396.
            if let Some((_, faces)) = self.my_bound_faces.get(&bat::shape_key(&bnd1)) {
                for itf in faces.clone() {
                    if set_add_bool(&mut merged_faces, &itf) {
                        list_faces1.push(itf.clone());
                    }
                }
            }
        } else {
            // Create presentation edge
            // OCCT L401-403: TopExp::Vertices(Edge1, V1, V2).
            let (mut v1, mut v2) = bat::top_exp_vertices_raw(&edge1);
            let mut v1 = v1.take().unwrap_or_else(Shape::null);
            let mut v2 = v2.take().unwrap_or_else(Shape::null);
            // OCCT L404-409.
            if let Some(v) = self.my_vertex_node.get(&bat::shape_key(&v1)) {
                v1 = v.1.clone();
            }
            if let Some(v) = self.my_vertex_node.get(&bat::shape_key(&v2)) {
                v2 = v.1.clone();
            }

            // OCCT L412-413: TopoDS_Edge NewEdge = Edge1; NewEdge.EmptyCopy();
            let mut new_edge = brep.empty_copied(&edge1);

            // Add the vertices
            // OCCT L416-419.
            let an_edge = bat::oriented(&new_edge, Orientation::Forward);
            let mut a_builder = BRepBuilder::new();
            a_builder.add_to_edge(brep, an_edge.clone(), bat::oriented(&v1, Orientation::Forward));
            a_builder.add_to_edge(brep, an_edge, bat::oriented(&v2, Orientation::Reversed));

            edge1 = new_edge;
        }

        // OCCT L422: bool isForward = true;
        let mut is_forward = true;

        // Merge candidate sections
        // OCCT L425: for (int i = 1; i <= seqEdges.Length(); i++).
        for i in 1..=seq_edges.len() {
            // Retrieve candidate section
            // OCCT L427: const TopoDS_Shape& oedge2 = seqEdges(i);
            let oedge2 = seq_edges[i - 1].clone();

            if self.my_sewing {
                // OCCT L432-439: reshape Edge2 through both reshapers.
                let a_tmp_shape = self.my_re_shape.apply(brep, &oedge2, ShapeType::Shape);
                let mut edge2 = a_tmp_shape;
                let a_tmp_shape = loc_re_shape.apply(brep, &edge2, ShapeType::Shape);
                if !std::ptr::eq(loc_re_shape as *const _, &self.my_re_shape as *const _) {
                    edge2 = a_tmp_shape;
                }

                // Calculate relative orientation
                // OCCT L442-446.
                let mut orientation = seq_forward[i - 1];
                if !is_forward {
                    orientation = !orientation;
                }

                // Retrieve faces information for the second edge
                // OCCT L449-456.
                let mut bnd2 = oedge2.clone();
                if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd2)) {
                    bnd2 = v.clone();
                }
                let list_faces2 = match self.my_bound_faces.get(&bat::shape_key(&bnd2)) {
                    Some((_, l)) => l.clone(),
                    None => continue, // Skip floating edge
                };

                // OCCT L459: int whichSec = 1;
                let mut which_sec: i32 = 1; // Indicates on which edge the pCurve has been reported
                let new_edge = self.same_parameter_edge_pair(
                    brep,
                    &edge1,
                    &edge2,
                    &list_faces1,
                    &list_faces2,
                    orientation,
                    &mut which_sec,
                    true,
                );
                // OCCT L463-466: if (NewEdge.IsNull()) continue;
                if new_edge.is_null() {
                    continue;
                }

                // Record faces information for the temporary merged edge
                // OCCT L469-475.
                for itf in &list_faces2 {
                    if set_add_bool(&mut merged_faces, itf) {
                        list_faces1.push(itf.clone());
                    }
                }

                // Record merged section orientation
                // OCCT L478-481.
                if !orientation && which_sec != 1 {
                    is_forward = !is_forward;
                }
                edge1 = new_edge;
            }

            // Append actually merged edge
            // OCCT L484: mapMerged.Add(oedge2);
            set_add_bool(map_merged, &oedge2);
            // OCCT L485: isDone = true;
            is_done = true;

            // OCCT L487-490: if (!myNonmanifold) break;
            if !self.my_nonmanifold {
                break;
            }
        }

        // OCCT L493-505.
        if is_done {
            // Change result orientation
            edge1.orientation = if is_forward { Orientation::Forward } else { Orientation::Reversed };
        } else {
            edge1 = Shape::null();
        }

        edge1
    }

    /// OCCT BRepBuilderAPI_Sewing::SameParameterEdge(edgeFirst, edgeLast,
    /// listFacesFirst, listFacesLast, secForward, whichSec, firstCall)
    /// (cxx L662-1188) — internal use; merges two section edges keeping the
    /// curve3d/curve2d/range/parametrization of the first section.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn same_parameter_edge_pair(
        &mut self,
        brep: &mut BRep,
        edge_first: &Shape,
        edge_last: &Shape,
        list_faces_first: &[Shape],
        list_faces_last: &[Shape],
        sec_forward: bool,
        which_sec: &mut i32,
        first_call: bool,
    ) -> Shape {
        // Do not process floating edges
        // OCCT L667-671.
        if list_faces_first.is_empty() || list_faces_last.is_empty() {
            return Shape::null();
        }

        // Sort input edges
        // OCCT L674: TopoDS_Edge edge1, edge2;
        let mut edge1: Shape;
        let mut edge2: Shape;
        if first_call {
            // Take the longest edge as first
            // OCCT L678-689.
            let (c3d1, f, l) = match bat::brep_tool_curve(edge_first) {
                Some(v) => v,
                None => return Shape::null(),
            };
            let len1 = super::gcpnts_abscissa_length(&c3d1, f, l);
            let (c3d2, f, l) = match bat::brep_tool_curve(edge_last) {
                Some(v) => v,
                None => return Shape::null(),
            };
            let len2 = super::gcpnts_abscissa_length(&c3d2, f, l);
            if len1 < len2 {
                edge1 = edge_last.clone();
                edge2 = edge_first.clone();
                *which_sec = 2;
            } else {
                edge1 = edge_first.clone();
                edge2 = edge_last.clone();
                *which_sec = 1;
            }
        } else {
            // OCCT L695-705.
            if *which_sec == 1 {
                edge1 = edge_last.clone();
                edge2 = edge_first.clone();
                *which_sec = 2;
            } else {
                edge1 = edge_first.clone();
                edge2 = edge_last.clone();
                *which_sec = 1;
            }
        }

        // OCCT L708-709: BRep_Tool::Range(edge1, first, last);
        let (mut first, mut last) = bat::brep_tool_range(&edge1);
        let _an_builder_holder = (); // OCCT L710: BRep_Builder aBuilder — carried by the pool builder below.

        // To keep NM vertices on edge
        // OCCT L713-714.
        let mut a_seq_nm_vert: Vec<Shape> = Vec::new();
        let mut a_seq_nm_pars: Vec<f64> = Vec::new();
        find_nm_vertices(brep, &edge1, &mut a_seq_nm_vert, &mut a_seq_nm_pars);
        find_nm_vertices(brep, &edge2, &mut a_seq_nm_vert, &mut a_seq_nm_pars);

        // Create new edge
        // OCCT L717-719.
        let mut edge = builder_make_edge(brep);
        edge.orientation = edge1.orientation;

        // Retrieve edge curve
        // OCCT L722-729.
        let c3d = match super::brep_tool_curve_world(brep, &edge1) {
            Some((c, f, l)) => {
                first = f;
                last = l;
                let _ = &mut first;
                c
            }
            None => return Shape::null(),
        };
        // OCCT L730-733: aBuilder.UpdateEdge(edge, c3d, Tol(edge1));
        // aBuilder.Range(edge, first, last); SameRange/SameParameter false.
        builder_update_edge_curve3d(brep, edge.clone(), &c3d, bat::brep_tool_tolerance(&edge1));
        let mut a_builder = BRepBuilder::new();
        a_builder.set_edge_range(brep, edge.clone(), first, last);
        a_builder.set_edge_same_range(brep, edge.clone(), false); // true
        a_builder.set_edge_same_parameter(brep, edge.clone(), false);
        // Create and add new vertices
        // OCCT L735-737.
        {
            // Retrieve original vertices from edges
            // OCCT L739-743.
            let (v11, v12) = bat::top_exp_vertices_raw(&edge1);
            let (v21, v22) = bat::top_exp_vertices_raw(&edge2);
            let v11 = v11.unwrap_or_else(Shape::null);
            let mut v12 = v12.unwrap_or_else(Shape::null);
            let mut v21 = v21.unwrap_or_else(Shape::null);
            let mut v22 = v22.unwrap_or_else(Shape::null);

            // check that edges merged valid way (for edges having length
            // less than specified tolerance; check if edges are closed)
            // OCCT L746-748.
            let is_closed1 = v11.is_same(&v12);
            let is_closed2 = v21.is_same(&v22);
            // OCCT L749-763.
            if !is_closed1 && !is_closed2 {
                if sec_forward {
                    if v11.is_same(&v22) || v12.is_same(&v21) {
                        return Shape::null();
                    }
                } else if v11.is_same(&v21) || v12.is_same(&v22) {
                    return Shape::null();
                }
            }

            // szv: do not reshape here!!!
            // V11 = TopoDS::Vertex(myReShape->Apply(V11));
            // V12 = TopoDS::Vertex(myReShape->Apply(V12));
            // V21 = TopoDS::Vertex(myReShape->Apply(V21));
            // V22 = TopoDS::Vertex(myReShape->Apply(V22));

            // OCCT L772-807.
            let v1_new: Shape;
            let v2_new: Shape;
            if is_closed1 || is_closed2 {
                // at least one of the edges is closed
                if is_closed1 && is_closed2 {
                    // both edges are closed
                    v1_new = compute_tolerance_vertex2(brep, &v11, &v21);
                } else if is_closed1 {
                    // only first edge is closed
                    v1_new = compute_tolerance_vertex3(brep, &v22, &v21, &v11);
                } else {
                    // only second edge is closed
                    v1_new = compute_tolerance_vertex3(brep, &v11, &v12, &v21);
                }
                v2_new = v1_new.clone();
            } else {
                // both edges are open
                // OCCT L813-817.
                let is_old_first = if sec_forward { v11.is_same(&v21) } else { v11.is_same(&v22) };
                let is_old_last = if sec_forward { v12.is_same(&v22) } else { v12.is_same(&v21) };
                let mut v1n = Shape::null();
                let mut v2n = Shape::null();
                if sec_forward {
                    // case if vertices already sewed
                    if !is_old_first {
                        v1n = compute_tolerance_vertex2(brep, &v11, &v21);
                    }
                    if !is_old_last {
                        v2n = compute_tolerance_vertex2(brep, &v12, &v22);
                    }
                } else {
                    if !is_old_first {
                        v1n = compute_tolerance_vertex2(brep, &v11, &v22);
                    }
                    if !is_old_last {
                        v2n = compute_tolerance_vertex2(brep, &v12, &v21);
                    }
                }
                // OCCT L832-838.
                if is_old_first {
                    v1n = v11.clone();
                }
                if is_old_last {
                    v2n = v12.clone();
                }
                v1_new = v1n;
                v2_new = v2n;
            }
            let _ = &mut v12;
            let _ = &mut v21;
            let _ = &mut v22;
            // Add the vertices in the good sense
            // OCCT L840-851.
            let an_edge = bat::oriented(&edge, Orientation::Forward);
            a_builder.add_to_edge(brep, an_edge.clone(), bat::oriented(&v1_new, Orientation::Forward));
            a_builder.add_to_edge(brep, an_edge.clone(), bat::oriented(&v2_new, Orientation::Reversed));

            for k in 1..=a_seq_nm_vert.len() {
                a_builder.add_to_edge(brep, an_edge.clone(), a_seq_nm_vert[k - 1].clone());
            }
        }

        // Retrieve second PCurves
        // OCCT L854-856.
        let mut loc2: u32;

        // OCCT L865-871.
        let itf2: Vec<Shape> = if *which_sec == 1 {
            list_faces_last.to_vec()
        } else {
            list_faces_first.to_vec()
        };
        // OCCT L872: bool isResEdge = false;
        let mut is_res_edge = false;
        for itf2_val in &itf2 {
            // OCCT L874-877.
            let mut c2d21: Option<Curve2d> = None;
            let mut first_old: f64 = 0.0;
            let mut last_old: f64 = 0.0;
            let fac2 = itf2_val.clone();

            // OCCT L879: surf2 = BRep_Tool::Surface(fac2, loc2);
            let (surf2, loc2v) = super::brep_tool_surface_loc(&fac2);
            let surf2 = match surf2 {
                Some(s) => s,
                None => continue,
            };
            loc2 = loc2v;
            // OCCT L880-882.
            let is_seam2 = (self.is_u_closed_surface(brep, &surf2, &edge2, &fac2, loc2)
                || self.is_v_closed_surface(brep, &surf2, &edge2, &fac2, loc2))
                && bat::brep_tool_is_closed_on_surface(&edge2, &fac2);
            if is_seam2 {
                if !self.my_nonmanifold {
                    return Shape::null();
                }
                // OCCT L886-888: c2d21 = CurveOnSurface(edge2.Reversed(), fac2, ...).
                let a_tmp_shape = bat::reversed(&edge2);
                if let Some((c, f, l)) = bat::brep_tool_curve_on_surface(&a_tmp_shape, &fac2) {
                    c2d21 = Some(c);
                    first_old = f;
                    last_old = l;
                }
            }
            // OCCT L889-890: c2d2 = CurveOnSurface(edge2, fac2, firstOld, lastOld).
            let c2d2 = bat::brep_tool_curve_on_surface(&edge2, &fac2);
            let mut c2d2 = match c2d2 {
                Some((c, f, l)) => {
                    first_old = f;
                    last_old = l;
                    c
                }
                None => {
                    if c2d21.is_none() {
                        continue;
                    } else {
                        // OCCT leaves c2d2 null and guards the copy below —
                        // the null-copy crash maps to skipping the face.
                        continue;
                    }
                }
            };

            if let Some(c21) = c2d21.clone() {
                let mut c2d21v = c21;
                // OCCT L897-910.
                if !sec_forward {
                    if matches!(c2d21v, Curve2d::Line(_)) {
                        c2d21v = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                            curve: Box::new(c2d21v),
                            t_min: first_old,
                            t_max: last_old,
                        });
                    }
                    let first2d = first_old; // BUG USA60321
                    let last2d = last_old;
                    first_old = c2d21v.reversed_parameter(last2d);
                    last_old = c2d21v.reversed_parameter(first2d);
                    c2d21v = curve2d_reversed(&c2d21v);
                }
                // OCCT L912: c2d21 = SameRange(c2d21, firstOld, lastOld, first, last).
                c2d21 = self.same_range(&c2d21v, first_old, last_old, first, last);
            }

            // Make second PCurve sameRange with the 3d curve
            // OCCT L916: c2d2 = Copy().
            let mut c2d2v = c2d2.clone();

            // OCCT L918-929.
            if !sec_forward {
                if matches!(c2d2v, Curve2d::Line(_)) {
                    c2d2v = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                        curve: Box::new(c2d2v),
                        t_min: first_old,
                        t_max: last_old,
                    });
                }
                let first2d = first_old;
                let last2d = last_old;
                first_old = c2d2v.reversed_parameter(last2d);
                last_old = c2d2v.reversed_parameter(first2d);
                c2d2v = curve2d_reversed(&c2d2v);
            }

            // OCCT L931: c2d2 = SameRange(c2d2, firstOld, lastOld, first, last).
            let c2d2 = self.same_range(&c2d2v, first_old, last_old, first, last);
            let mut c2d2 = match c2d2 {
                Some(c) => c,
                None => continue,
            };

            // Add second PCurve
            // OCCT L935-938.
            let mut is_seam = false;
            let mut ori = Orientation::Forward;

            // OCCT L941-947.
            let itf1: Vec<Shape> = if *which_sec == 1 {
                list_faces_first.to_vec()
            } else {
                list_faces_last.to_vec()
            };
            for itf1_val in &itf1 {
                if is_seam {
                    break;
                }
                let mut c2d11: Option<Curve2d> = None;
                let fac1 = itf1_val.clone();

                // OCCT L950-951: surf1 = BRep_Tool::Surface(fac1, loc1).
                let (surf1, loc1) = super::brep_tool_surface_loc(&fac1);
                let surf1 = match surf1 {
                    Some(s) => s,
                    None => continue,
                };

                // OCCT L953-957.
                let is_seam1 = (self.is_u_closed_surface(brep, &surf1, &edge1, &fac1, loc1)
                    || self.is_v_closed_surface(brep, &surf1, &edge1, &fac1, loc1))
                    && bat::brep_tool_is_closed_on_surface(&edge1, &fac1);
                let mut c2d1: Option<Curve2d> =
                    bat::brep_tool_curve_on_surface(&edge1, &fac1).map(|(c, _f, _l)| c);
                ori = edge1.orientation;
                if fac1.orientation == Orientation::Reversed {
                    ori = top_abs_reverse(ori);
                }

                if is_seam1 {
                    if !self.my_nonmanifold {
                        return Shape::null();
                    }
                    // OCCT L967-970.
                    let a_tmp_shape = bat::reversed(&edge1);
                    c2d11 = bat::brep_tool_curve_on_surface(&a_tmp_shape, &fac1)
                        .map(|(c, _f, _l)| c);
                    // OCCT L972-978.
                    // INFO: unlike the UpdateEdge sites below, c2d1 here may
                    // be null in OCCT — L966 has no null-check continue
                    // before the L984/L988/L993 calls.  A null pcurve in
                    // BRep_Builder::UpdateEdge does NOT raise: UpdateCurves
                    // treats it as a REMOVAL ("remove the pcurves on <S>
                    // from <lcr> if <C1> or <C2> is null", BRep_Builder.cxx
                    // L248/L288 and L101/L149) plus a no-op
                    // UpdateTolerance(0), and on the freshly created merged
                    // edge (no prior representation on fac1's surface) that
                    // removal is a no-op.  rcad cannot store a null-pcurve
                    // representation (update_edge_pcurve takes Curve2d by
                    // value), so the skip below is the equivalent outcome;
                    // the L985-988 continue (both null) matches OCCT
                    // L996-999.
                    if let Some(c1) = c2d1.clone() {
                        if let Some(c11) = c2d11.clone() {
                            if ori == Orientation::Forward {
                                a_builder.update_edge_pcurve_closed(
                                    brep, edge.clone(), c1, c11, fac1.clone(), 0.0,
                                );
                            } else {
                                a_builder.update_edge_pcurve_closed(
                                    brep, edge.clone(), c11, c1, fac1.clone(), 0.0,
                                );
                            }
                        }
                    }
                } else if let Some(c1) = c2d1.clone() {
                    // OCCT L981-982.
                    a_builder.update_edge_pcurve(brep, edge.clone(), c1, fac1.clone(), 0.0);
                }

                // OCCT L985-988.
                if c2d1.is_none() && c2d11.is_none() {
                    continue;
                }

                // OCCT L991: if (surf2 == surf1).
                if rcad_kernel::topo::topods::surface_same(&surf2, &surf1) {
                    // Merge sections which are on the same face
                    // OCCT L993-1010.
                    if loc2 == loc1 {
                        let uclosed = self.is_u_closed_surface(brep, &surf2, &edge2, &fac2, loc2);
                        let vclosed = self.is_v_closed_surface(brep, &surf2, &edge2, &fac2, loc2);
                        if uclosed || vclosed {
                            // OCCT L1000-1004: c2d1 is non-null here — the
                            // L985-988 continue (cxx L996-999) filtered the
                            // both-null state and the rcad lookups are
                            // orientation-independent, so c2d1/c2d11 are
                            // Some/None together; the unwrap below never
                            // fires.
                            let c1 = c2d1
                                .clone()
                                .expect("sewing SameParameterEdge: null c2d1 past the both-null continue (unreachable; cxx L1010 dereferences it)");
                            let pf = c1.default_domain()[0];
                            let p1n = c1.point_at(first.max(pf));
                            let [c2f, _] = c2d2.default_domain();
                            let p21n = c2d2.point_at(first.max(c2f));
                            let [_, c2l] = c2d2.default_domain();
                            let p22n = c2d2.point_at(last.min(c2l));
                            let a_dist = (p1n.distance(p21n)).min(p1n.distance(p22n));
                            // OCCT L1005-1006: surf2->Bounds(U1, U2, V1, V2).
                            let [su1, su2, sv1, sv2] = surf2.default_domain();
                            let is_seam_v = (uclosed && a_dist > 0.75 * (su2 - su1).abs())
                                || (vclosed && a_dist > 0.75 * (sv2 - sv1).abs());
                            is_seam = is_seam_v;
                            // OCCT L1007-1010.
                            if !is_seam && bat::brep_tool_is_closed_on_surface(&edge, &fac1) {
                                continue;
                            }
                        }
                    }
                }

                // OCCT L1013: isResEdge = true;
                is_res_edge = true;
                // OCCT L1014-1046.
                // INFO: c2d1 is non-null in every state reaching this point
                // (the L985-988 continue filters the both-null state; see
                // the surf2==surf1 block above) — OCCT L1010 dereferences
                // c2d1 before isSeam can become true, so the unwraps below
                // never fire.
                if is_seam {
                    let c1 = c2d1
                        .clone()
                        .expect("sewing SameParameterEdge: null c2d1 in the isSeam UpdateEdge (unreachable; cxx L1010 dereferences it before L1018)");
                    if ori == Orientation::Forward {
                        // OCCT L1016-1018: UpdateEdge(edge, c2d1, c2d2,
                        // surf2, loc2, Confusion) — stored on fac2 (the
                        // face whose surface is surf2); the sibling-face
                        // reads resolve through the surface walk.
                        a_builder.update_edge_pcurve_closed(
                            brep, edge.clone(), c1, c2d2.clone(), fac2.clone(), CONFUSION,
                        );
                    } else {
                        a_builder.update_edge_pcurve_closed(
                            brep, edge.clone(), c2d2.clone(), c1, fac2.clone(), CONFUSION,
                        );
                    }
                } else if is_seam2 {
                    // OCCT L1030-1044.
                    // INFO: c2d21 is non-null whenever is_seam2 holds — the
                    // rcad pcurve lookup is orientation-independent, so the
                    // L498-503 reversed-edge read succeeds exactly when the
                    // L506 edge2/fac2 read does, and its None arm continued
                    // (cxx L896-899); on the OCCT side IsClosed(E, F) (cxx
                    // L884-885) already implies a matching closed
                    // representation whose PCurve/PCurve2 the reversed read
                    // returns non-null (BRep_Tool.cxx L354-357).  The
                    // unwraps below never fire.
                    let mut init_ori = edge2.orientation;
                    let mut sec_ori = edge.orientation;
                    if fac2.orientation == Orientation::Reversed {
                        init_ori = top_abs_reverse(init_ori);
                        sec_ori = top_abs_reverse(sec_ori);
                    }
                    if !sec_forward {
                        init_ori = top_abs_reverse(init_ori);
                    }

                    let c21 = c2d21
                        .clone()
                        .expect("sewing SameParameterEdge: null c2d21 in the isSeam2 UpdateEdge (unreachable; cxx L893 fills it whenever isSeam2)");
                    if init_ori == Orientation::Forward {
                        a_builder.update_edge_pcurve_closed(
                            brep, edge.clone(), c2d2.clone(), c21, fac2.clone(), CONFUSION,
                        );
                    } else {
                        a_builder.update_edge_pcurve_closed(
                            brep, edge.clone(), c21, c2d2.clone(), fac2.clone(), CONFUSION,
                        );
                    }
                } else {
                    // OCCT L1045-1046.
                    a_builder.update_edge_pcurve(
                        brep, edge.clone(), c2d2.clone(), fac2.clone(), CONFUSION,
                    );
                }
            }
        }
        // OCCT L1048-1064.
        let mut tol_reached = INFINITE_VALUE;
        let mut is_same_par = false;
        // OCCT: try { if (isResEdge) SameParameter(edge); if
        // (BRep_Tool::SameParameter(edge)) {...} } catch { isSamePar =
        // false; } — the rcad SameParameter carrier does not throw.
        if is_res_edge {
            self.same_parameter(brep, &edge);
        }
        if edge_same_parameter(&edge) {
            is_same_par = true;
            tol_reached = bat::brep_tool_tolerance(&edge);
        }
        let _ = is_same_par;

        // OCCT L1066-1147.
        if first_call && (!is_res_edge || !is_same_par || tol_reached > self.my_tolerance) {
            let mut which_secn = *which_sec;
            // Try to merge on the second section
            let mut second_ok = false;
            let s_edge = self.same_parameter_edge_pair(
                brep,
                edge_first,
                edge_last,
                list_faces_first,
                list_faces_last,
                sec_forward,
                &mut which_secn,
                false,
            );
            if !s_edge.is_null() {
                let tol_reached_2 = bat::brep_tool_tolerance(&s_edge);
                second_ok =
                    edge_same_parameter(&s_edge) && tol_reached_2 < tol_reached;
                if second_ok {
                    edge = s_edge;
                    *which_sec = which_secn;
                    tol_reached = tol_reached_2;
                }
            }

            if !second_ok && !edge.is_null() {
                // OCCT L1150-1152: GeomAdaptor_Curve c3dAdapt(c3d).
                // Discretize edge curve
                let (nbp, delta_t0) = (23usize, (last - first) / (23 - 1) as f64);
                let _ = delta_t0;
                // OCCT L1153-1157: the 3d sample points of the curve.
                let first3d = first;
                let last3d = last;
                let delta_t = (last3d - first3d) / (nbp as f64 - 1.0);
                let mut c3dpnt: Vec<DVec3> = vec![DVec3::ZERO; nbp];
                use rcad_kernel::geom::CurveEval;
                for i in 1..=nbp {
                    c3dpnt[i - 1] = c3d.point_at(first3d + (i as f64 - 1.0) * delta_t);
                }

                let mut dist = 0.0f64;
                let mut max_tol = -1.0f64;
                let mut more = true;

                let mut j = 1usize;
                while more {
                    // OCCT L1160: BRep_Tool::CurveOnSurface(edge, c2d2,
                    // surf2, loc2, first, last, j) — the j-th representation.
                    let rep = edge_pcurve_at(brep, &edge, j as i32);
                    let (c2d2j, surf2j, loc2j, fr, lr) = match rep {
                        Some(v) => v,
                        None => break,
                    };

                    more = true;
                    {
                        // OCCT L1165-1170.
                        let a_s = if loc2j != 0 {
                            rcad_kernel::geom::transform_surface(
                                &surf2j,
                                &brep.get_location(loc2j),
                            )
                        } else {
                            surf2j.clone()
                        };

                        let mut dist2 = 0.0f64;
                        let delta_t = (lr - fr) / (nbp as f64 - 1.0);
                        for i in 1..=nbp {
                            let a_p2d = c2d2j.point_at(fr + (i as f64 - 1.0) * delta_t);
                            let a_p2 = a_s.point_at(a_p2d.x, a_p2d.y);
                            let a_p1 = c3dpnt[i - 1];
                            dist = a_p2.distance_squared(a_p1);
                            if dist > dist2 {
                                dist2 = dist;
                            }
                        }
                        max_tol = (dist2.sqrt() * (1.0 + 1e-7)).max(CONFUSION);
                    }
                    j += 1;
                }
                // OCCT L1195-1206.
                if max_tol >= 0.0 && max_tol < tol_reached {
                    if tol_reached > self.max_tolerance() {
                        // Set tolerance directly to overwrite too large
                        // tolerance
                        if let TShape::Edge(ed) = Arc::make_mut(&mut edge.data) {
                            ed.tolerance = max_tol;
                        }
                    } else {
                        // just update tolerance with computed distance
                        let mut b = BRepBuilder::new();
                        b.update_edge_tolerance(brep, edge.clone(), max_tol);
                    }
                }
                let mut b = BRepBuilder::new();
                b.set_edge_same_parameter(brep, edge.clone(), true);
            }
        }

        // OCCT L1210-1214.
        let tol_edge1 = bat::brep_tool_tolerance(&edge);
        if tol_edge1 > self.max_tolerance() {
            edge = Shape::null();
        }
        edge
    }
}

/// OCCT `BRep_Tool::SameParameter(E)` (BRep_Tool.cxx L228-244) — the edge
/// TShape flag.
fn edge_same_parameter(e: &Shape) -> bool {
    matches!(&*e.data, TShape::Edge(ed) if ed.same_parameter)
}

/// OCCT static findNMVertices(theEdge, theSeqNMVert, theSeqPars)
/// (cxx L517-574).
pub(crate) fn find_nm_vertices(
    brep: &BRep,
    the_edge: &Shape,
    the_seq_nm_vert: &mut Vec<Shape>,
    the_seq_pars: &mut Vec<f64>,
) -> bool {
    // OCCT L518: TopoDS_Iterator aItV(theEdge, false).
    for sub in iter_no_cumori(the_edge) {
        if sub.orientation == Orientation::Internal || sub.orientation == Orientation::External {
            the_seq_nm_vert.push(sub);
        }
    }
    let nb_v = the_seq_nm_vert.len();
    if nb_v == 0 {
        return false;
    }
    // OCCT L528-532.
    let (c3d, first, last) = match bat::brep_tool_curve(the_edge) {
        Some(v) => v,
        None => return false,
    };
    let _ = brep;
    // OCCT L534-536: Extrema_ExtPC locProj; locProj.Initialize(GAC, first, last).
    let pfirst = c3d.point_at(first);
    let plast = c3d.point_at(last);

    for i in 1..=nb_v {
        let a_v = the_seq_nm_vert[i - 1].clone();
        let pt = bat::brep_tool_pnt(&a_v).unwrap_or(DVec3::ZERO);

        let dist_f2 = pfirst.distance_squared(pt);
        let dist_l2 = plast.distance_squared(pt);
        let mut apar = if dist_f2 > dist_l2 { last } else { first };

        // Project current point on curve
        // OCCT L543-544: locProj.Perform(pt).
        let adaptor = rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor::new(c3d.clone());
        let gac = rcad_kernel::base::extrema_curve_tool::CurveToolHandle::for_curve3(
            &c3d, &adaptor, &adaptor,
        );
        let mut loc_proj = rcad_kernel::base::extrema_ext_pc::ExtremaExtPC::new_point_curve_ranged(
            pt, &gac, first, last, 1.0e-10,
        );
        if loc_proj.is_done() && loc_proj.nb_ext() > 0 {
            let mut dist2_min = dist_f2.min(dist_l2);
            let mut ind_min: usize = 0;
            for ind in 1..=loc_proj.nb_ext() {
                let d_proj2 = loc_proj.square_distance(ind);
                if d_proj2 < dist2_min {
                    ind_min = ind;
                    dist2_min = d_proj2;
                }
            }
            if ind_min != 0 {
                apar = loc_proj.point(ind_min).param;
            }

            the_seq_pars.push(apar);
        }
    }
    true
}

/// OCCT static ComputeToleranceVertex(theV1, theV2, theNewV) (cxx L576-621).
pub(crate) fn compute_tolerance_vertex2(brep: &mut BRep, the_v1: &Shape, the_v2: &Shape) -> Shape {
    // OCCT L577-587.
    let a_eps = f64::EPSILON; // RealEpsilon()
    let a_v = [the_v1.clone(), the_v2.clone()];
    let mut a_p = [DVec3::ZERO; 2];
    let mut a_r = [0.0f64; 2];
    for m in 0..2 {
        a_p[m] = bat::brep_tool_pnt(&a_v[m]).unwrap_or(DVec3::ZERO);
        a_r[m] = bat::brep_tool_tolerance(&a_v[m]);
    }
    // OCCT L590-595: m = 0 (max R), n = 1 (min R).
    let (m, n) = if a_r[0] < a_r[1] { (1usize, 0usize) } else { (0usize, 1usize) };
    // OCCT L597-599.
    let d_r = a_r[m] - a_r[n]; // dR >= 0.
    let a_vd = a_p[n] - a_p[m]; // gp_Vec aVD(aP[m], aP[n]).
    let a_d = a_vd.length();
    // OCCT L601-614.
    if a_d <= d_r || a_d < a_eps {
        make_vertex_at(brep, a_p[m], a_r[m])
    } else {
        let a_rr = 0.5 * (a_r[m] + a_r[n] + a_d);
        let a_xyzr = (a_p[m] + a_p[n] - a_vd * (d_r / a_d)) * 0.5;
        make_vertex_at(brep, a_xyzr, a_rr)
    }
}

/// OCCT static ComputeToleranceVertex(theV1, theV2, theV3, theNewV)
/// (cxx L623-660).
pub(crate) fn compute_tolerance_vertex3(
    brep: &mut BRep,
    the_v1: &Shape,
    the_v2: &Shape,
    the_v3: &Shape,
) -> Shape {
    // OCCT L624-637.
    let a_v = [the_v1.clone(), the_v2.clone(), the_v3.clone()];
    let mut a_p = [DVec3::ZERO; 3];
    let mut a_r = [0.0f64; 3];
    let mut a_xyz = DVec3::ZERO;
    for i in 0..3 {
        a_p[i] = bat::brep_tool_pnt(&a_v[i]).unwrap_or(DVec3::ZERO);
        a_r[i] = bat::brep_tool_tolerance(&a_v[i]);
        a_xyz += a_p[i];
    }
    // OCCT L639-642.
    a_xyz /= 3.0;
    let a_center = a_xyz;
    // OCCT L644-657.
    let mut a_dmax = -1.0f64;
    for i in 0..3 {
        let mut a_di = a_center.distance(a_p[i]);
        a_di += a_r[i];
        if a_di > a_dmax {
            a_dmax = a_di;
        }
    }
    make_vertex_at(brep, a_center, a_dmax)
}

/// OCCT `Geom2d_Curve::Reverse()` — the canonical kernel body (the
/// BOPTools_AlgoTools2D precedent).
pub(crate) fn curve2d_reversed(c: &Curve2d) -> Curve2d {
    rcad_kernel::geom::reverse_curve2d(c)
}

/// OCCT `BRep_Tool::CurveOnSurface(E, C, S, L, First, Last, i)` — the j-th
/// curve-on-surface representation of an edge (BRep_Tool.cxx L272-306
/// iteration form); returns (pcurve, surface, location, first, last).
fn edge_pcurve_at(
    brep: &BRep,
    edge: &Shape,
    index: i32,
) -> Option<(Curve2d, Surface3, u32, f64, f64)> {
    let ed = match edge.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    // The curve-on-surface representations in storage order.
    let mut seen = 0i32;
    for r in &ed.representations {
        if let rcad_kernel::topo::topods::CurveRepresentation::CurveOnSurface {
            face,
            pcurve,
            range,
        } = r
        {
            seen += 1;
            if seen != index {
                continue;
            }
            // Resolve the stored face pointer to the surface + location.
            let ts = brep
                .tshapes
                .iter()
                .find(|ts| Arc::as_ptr(ts) as u64 == face.0)?;
            if let TShape::Face(fd) = ts.as_ref() {
                let surf = fd.surface.clone()?;
                let loc = if face.1 != 0 { face.1 } else { fd.surface_location };
                return Some((pcurve.clone(), surf, loc, range[0], range[1]));
            }
            return None;
        }
    }
    None
}
