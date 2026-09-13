// OCCT LocOpe_SplitShape.cxx L17-1776 — local re-hosts (module b).
//
// This module carries the BRep_Tool / TopExp / BRep_Builder / BRepTools /
// BRepAdaptor_Surface / BRepTopAdaptor_FClass2d / BRepLib_MakeWire re-hosts
// consumed by the LocOpeSplitShape bodies of loc_ope_split_shape.rs. The
// architecture differences are numbered in the main module header.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitShape.cxx

use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_surface,
    builder_add_edge_vertex, builder_add_face_wire, builder_add_wire_edge, oriented,
    shape_is_closed, top_exp_vertices_raw,
};
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_spliter::ShapeSet;
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_degenerated;
use crate::topalgo::brep_lib_make_wire::{builder_update_vertex_parameter, MakeWire};
use crate::topalgo::shape_source::FaceShapeSource;
use glam::{DAffine3, DVec2};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use rcad_kernel::precision::PCONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Local re-hosts (arch. diffs. #2/#3/#4/#5/#10).
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
pub(crate) fn with_orientation(s: &Shape, the_or: Orientation) -> Shape {
    oriented(s, the_or)
}

/// OCCT TopAbs::Reverse(theOr) (TopAbs.hxx L80-91).
pub(crate) fn reverse_orientation(o: Orientation) -> Orientation {
    crate::brep_algo::tool::top_abs_reverse(o)
}

/// OCCT TopoDS_Iterator(S) with CumOri = false — the stored child
/// orientations, NOT composed with the parent (TopoDS_Iterator.cxx L72-80).
/// Used where the OCCT source reads `it.Value().Orientation()` as the raw
/// child orientation (Rebuild cxx L1437) and where Put walks the children.
pub(crate) fn sub_shapes_no_cum_ori(the_s: &Shape) -> Vec<Shape> {
    match the_s.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => vec![ed.first.clone(), ed.last.clone()],
        TShape::Wire(wd) => wd.edges.clone(),
        TShape::Face(fd) => {
            let mut out: Vec<Shape> = Vec::new();
            if matches!(fd.outer_wire.data.as_ref(), TShape::Wire(_)) {
                out.push(fd.outer_wire.clone());
            }
            out.extend(fd.inner_wires.iter().cloned());
            out.extend(fd.internal_vertices.iter().cloned());
            out
        }
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Solid(sd) => {
            let mut out: Vec<Shape> = Vec::new();
            out.extend(sd.shells.iter().cloned());
            out.extend(sd.internal_vertices.iter().cloned());
            out.extend(sd.internal_edges.iter().cloned());
            out
        }
        TShape::CompSolid(cs) => cs.clone(),
        TShape::Compound(cd) => cd.clone(),
    }
}

/// OCCT BRep_Builder::Add(parent, child) for the container kinds Rebuild and
/// RebuildImage build (BRep_Builder.cxx Add overloads; the typed helpers of
/// crate::brep_algo::tool cover Wire/Face/Edge).
pub(crate) fn builder_add_shape(the_parent: &mut Shape, the_child: &Shape) {
    match the_parent.shape_type() {
        ShapeType::Wire => builder_add_wire_edge(the_parent, the_child),
        ShapeType::Face => builder_add_face_wire(the_parent, the_child),
        ShapeType::Edge => builder_add_edge_vertex(the_parent, the_child),
        _ => {
            if let TShape::Shell(sd) = Arc::make_mut(&mut the_parent.data) {
                sd.faces.push(the_child.clone());
                return;
            }
            if let TShape::Solid(sd) = Arc::make_mut(&mut the_parent.data) {
                sd.shells.push(the_child.clone());
                return;
            }
            if let TShape::CompSolid(cs) = Arc::make_mut(&mut the_parent.data) {
                cs.push(the_child.clone());
                return;
            }
            if let TShape::Compound(cd) = Arc::make_mut(&mut the_parent.data) {
                cd.push(the_child.clone());
            }
        }
    }
}

/// OCCT BRep_Builder::Add(E, V) + BRep_Builder::UpdateVertex(V, Par, E, Tol):
/// the Add must precede the Update (the OCCT ori walk reads the edge
/// children) and the mutated vertex is re-bound afterwards (the rcad repair
/// of the Arc fork; topalgo::brep_lib_make_wire precedent).
pub(crate) fn add_vertex_and_update(the_e: &mut Shape, the_v: &mut Shape, the_par: f64, the_tol: f64) {
    builder_add_edge_vertex(the_e, the_v);
    builder_update_vertex_parameter(the_v, the_par, the_e, the_tol);
    builder_add_edge_vertex(the_e, the_v);
}

/// OCCT BRep_Builder::UpdateEdge(E, Tol) — BRep_TEdge::UpdateTolerance keeps
/// the max.
pub(crate) fn builder_update_edge_tol(the_e: &mut Shape, the_tol: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) — the regularity stored on the
/// (E, F1/F2) curve representations (offset/draft_modification.rs re-host).
pub(crate) fn has_continuity(the_e: &Shape, the_f1: &Shape, the_f2: &Shape) -> bool {
    let (Some(s1), Some(s2)) = (brep_tool_surface(the_f1), brep_tool_surface(the_f2)) else {
        return false;
    };
    if let TShape::Edge(ed) = the_e.data.as_ref() {
        return ed.representations.iter().any(|cr| cr.is_regularity_on(&s1, &s2, 0, 0));
    }
    false
}

/// OCCT BRep_Builder::Continuity(E, F1, F2, GeomAbs_CN) (BRep_Builder.cxx
/// L1012-1043) — an existing matching CurveOn2Surfaces takes the continuity,
/// otherwise a new one is appended.
pub(crate) fn builder_continuity(the_e: &mut Shape, the_f1: &Shape, the_f2: &Shape) {
    let (Some(s1), Some(s2)) = (brep_tool_surface(the_f1), brep_tool_surface(the_f2)) else {
        return;
    };
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        // OCCT BRep_Builder::UpdateCurves: an existing matching
        // BRep_CurveOn2Surfaces takes the continuity.
        for i in 0..ed.representations.len() {
            if ed.representations[i].is_regularity_on(&s1, &s2, 0, 0) {
                if let rcad_kernel::topods::CurveRepresentation::CurveOn2Surfaces {
                    continuity,
                    ..
                } = &mut ed.representations[i]
                {
                    *continuity = rcad_kernel::topods::GeomAbsShape::CN;
                    return;
                }
            }
        }
        ed.representations
            .push(rcad_kernel::topods::CurveRepresentation::CurveOn2Surfaces {
                surface1: s1,
                surface2: s2,
                location1: 0,
                location2: 0,
                continuity: rcad_kernel::topods::GeomAbsShape::CN,
            });
    }
}

/// OCCT BRep_Tool::IsClosed(E, F) — the pcurve on F is a closed-surface
/// representation (BRep_Tool.cxx L814-841).
pub(crate) fn brep_tool_is_closed_on_surface(the_e: &Shape, the_f: &Shape) -> bool {
    crate::brep_algo::tool::brep_tool_is_closed_on_surface(the_e, the_f)
}

/// OCCT BRep_Tool::IsClosed(theShape) (BRep_Tool.cxx L1707-1757).
pub(crate) fn brep_tool_is_closed_shape(the_shape: &Shape) -> bool {
    let the_shape = &oriented(the_shape, Orientation::Forward);
    if the_shape.shape_type() == ShapeType::Shell {
        let mut a_map = ShapeSet::new();
        let mut has_bound = false;
        for e in explorer(the_shape, ShapeType::Edge, ShapeType::Shape) {
            if brep_tool_degenerated(&e)
                || e.orientation == Orientation::Internal
                || e.orientation == Orientation::External
            {
                continue;
            }
            has_bound = true;
            if !a_map.add(&e) {
                a_map.remove(&e);
            }
        }
        return has_bound && a_map.is_empty();
    } else if the_shape.shape_type() == ShapeType::Wire {
        let mut a_map = ShapeSet::new();
        let mut has_bound = false;
        for v in explorer(the_shape, ShapeType::Vertex, ShapeType::Shape) {
            if v.orientation == Orientation::Internal || v.orientation == Orientation::External {
                continue;
            }
            has_bound = true;
            if !a_map.add(&v) {
                a_map.remove(&v);
            }
        }
        return has_bound && a_map.is_empty();
    } else if the_shape.shape_type() == ShapeType::Edge {
        let (a_v_first, a_v_last) = top_exp_vertices_raw(the_shape);
        return match (a_v_first, a_v_last) {
            (Some(a), Some(b)) => a.is_same(&b),
            _ => false,
        };
    }
    shape_is_closed(the_shape)
}

/// OCCT BRepTools::OuterWire(F) (BRepTools.cxx L556-590).
pub(crate) fn brep_tools_outer_wire(the_f: &Shape) -> Shape {
    let wires = explorer(the_f, ShapeType::Wire, ShapeType::Shape);
    let Some(first) = wires.first() else {
        return Shape::null();
    };
    let mut wres = first.clone();
    if wires.len() > 1 {
        let mut bounds = brep_tools_uv_bounds_of_wire(the_f, &wres);
        for w in wires.iter().skip(1) {
            let b = brep_tools_uv_bounds_of_wire(the_f, w);
            if (b[0] - bounds[0]) <= PCONFUSION
                && (b[1] - bounds[1]) >= -PCONFUSION
                && (b[2] - bounds[2]) <= PCONFUSION
                && (b[3] - bounds[3]) >= -PCONFUSION
            {
                wres = w.clone();
                bounds = b;
            }
        }
    }
    wres
}

/// OCCT BRepTools::UVBounds(F, W, UMin, UMax, VMin, VMax) (BRepTools.cxx
/// L137-181 through AddUVBounds) — the union of the 2D boxes of the wire's
/// edge pcurves on the face, as [umin, umax, vmin, vmax].
pub(crate) fn brep_tools_uv_bounds_of_wire(the_f: &Shape, the_w: &Shape) -> [f64; 4] {
    let mut umin = rcad_kernel::precision::REAL_LAST;
    let mut umax = -rcad_kernel::precision::REAL_LAST;
    let mut vmin = rcad_kernel::precision::REAL_LAST;
    let mut vmax = -rcad_kernel::precision::REAL_LAST;
    for e in explorer(the_w, ShapeType::Edge, ShapeType::Shape) {
        if let Some((c2d, f, l)) = brep_tool_curve_on_surface(&e, the_f) {
            let b = rcad_kernel::base::bnd_lib::curve2d_bounding_box(&c2d, f, l, 0.0);
            umin = umin.min(b[0]);
            umax = umax.max(b[1]);
            vmin = vmin.min(b[2]);
            vmax = vmax.max(b[3]);
        }
    }
    [umin, umax, vmin, vmax]
}

/// OCCT BRepTools::Update(F) (BRepTools.cxx L383-390) — architecture
/// difference #4: the rcad TEdgeData carries no UV-point cache.
pub(crate) fn brep_tools_update(_the_f: &Shape) {}

/// OCCT BRepAdaptor_Surface (BRepAdaptor_Surface.cxx) — the surface readings
/// the OpenWire/ClosedWire paths consume (architecture difference #6).
pub(crate) struct BRepAdaptorSurface {
    my_surface: Option<Surface3>,
}

impl BRepAdaptorSurface {
    /// OCCT BRepAdaptor_Surface(F, restriction = false).
    pub(crate) fn new(the_f: &Shape) -> Self {
        BRepAdaptorSurface {
            my_surface: brep_tool_surface(the_f),
        }
    }

    /// OCCT BRepAdaptor_Surface::IsUPeriodic().
    pub(crate) fn is_u_periodic(&self) -> bool {
        matches!(
            self.my_surface,
            Some(Surface3::Cylinder(_))
                | Some(Surface3::Cone(_))
                | Some(Surface3::Sphere(_))
                | Some(Surface3::Torus(_))
        )
    }

    /// OCCT BRepAdaptor_Surface::IsVPeriodic().
    pub(crate) fn is_v_periodic(&self) -> bool {
        matches!(
            self.my_surface,
            Some(Surface3::Sphere(_)) | Some(Surface3::Torus(_))
        )
    }

    /// OCCT BRepAdaptor_Surface::UPeriod() (GeomAdaptor_Surface::UPeriod).
    pub(crate) fn u_period(&self) -> f64 {
        use std::f64::consts::TAU;
        match self.my_surface {
            Some(Surface3::Cylinder(_))
            | Some(Surface3::Cone(_))
            | Some(Surface3::Sphere(_))
            | Some(Surface3::Torus(_)) => TAU,
            _ => 0.0,
        }
    }

    /// OCCT BRepAdaptor_Surface::VPeriod() (GeomAdaptor_Surface::VPeriod).
    pub(crate) fn v_period(&self) -> f64 {
        use std::f64::consts::TAU;
        match self.my_surface {
            Some(Surface3::Sphere(_)) => TAU,
            Some(Surface3::Torus(t)) => 2.0 * t.minor_radius,
            _ => 0.0,
        }
    }

    /// OCCT BRepAdaptor_Surface::UResolution(R3d).
    pub(crate) fn u_resolution(&self, r3d: f64) -> f64 {
        match &self.my_surface {
            Some(s) => rcad_kernel::topo::topods::u_resolution_for_surface(s, r3d),
            None => 0.0,
        }
    }

    /// OCCT BRepAdaptor_Surface::VResolution(R3d).
    pub(crate) fn v_resolution(&self, r3d: f64) -> f64 {
        match &self.my_surface {
            Some(s) => rcad_kernel::topo::topods::v_resolution_for_surface(s, r3d),
            None => 0.0,
        }
    }

    /// OCCT BRepAdaptor_Surface::D0(U, V, P).
    pub(crate) fn d0(&self, u: f64, v: f64) -> glam::DVec3 {
        match &self.my_surface {
            Some(s) => s.point_at(u, v),
            None => glam::DVec3::ZERO,
        }
    }
}

/// OCCT BRepTopAdaptor_FClass2d(newFace, Precision::PConfusion())
/// (architecture difference #5).
pub(crate) fn face_classifier(the_face: &Shape) -> crate::topalgo::brep_top_adaptor::fclass2d::FClass2d {
    let surf = brep_tool_surface(the_face).expect("FClass2d face surface");
    let src = FaceShapeSource::new(the_face, surf, &[DAffine3::IDENTITY]);
    let _ = src;
    crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(&src, 0, PCONFUSION)
}

/// OCCT BRepTopAdaptor_FClass2d::PerformInfinitePoint().
pub(crate) fn classif_perform_infinite_point(
    the_face: &Shape,
    the_classif: &crate::topalgo::brep_top_adaptor::fclass2d::FClass2d,
) -> rcad_kernel::topods::State {
    let surf = brep_tool_surface(the_face).expect("FClass2d face surface");
    let src = FaceShapeSource::new(the_face, surf, &[DAffine3::IDENTITY]);
    the_classif.perform_infinite_point(&src)
}

/// OCCT BRepTopAdaptor_FClass2d::Perform(P2d).
pub(crate) fn classif_perform_p(
    the_face: &Shape,
    the_p2d: DVec2,
) -> rcad_kernel::topods::State {
    let surf = brep_tool_surface(the_face).expect("FClass2d face surface");
    let src = FaceShapeSource::new(the_face, surf, &[DAffine3::IDENTITY]);
    let classif = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(&src, 0, PCONFUSION);
    classif.perform(&src, the_p2d, false)
}

/// OCCT gp_Pnt2d::Distance / Precision::IsNegativeInfinite.
pub(crate) fn is_negative_infinite(v: f64) -> bool {
    v <= -rcad_kernel::precision::INFINITE_VALUE
}

/// OCCT Precision::IsPositiveInfinite.
pub(crate) fn is_positive_infinite(v: f64) -> bool {
    v >= rcad_kernel::precision::INFINITE_VALUE
}

/// OCCT BRepLib_MakeWire::Wire().Closed() — the Closed flag of the wire under
/// construction (BRepBuilderAPI_MakeShape::Shape()).
pub(crate) fn make_wire_closed(the_mw: &MakeWire) -> bool {
    shape_is_closed(the_mw.wire())
}
