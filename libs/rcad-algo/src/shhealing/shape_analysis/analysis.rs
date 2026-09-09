//! OCCT ShapeAnalysis package statics (TKShHealing — `ShapeAnalysis.hxx` +
//! `ShapeAnalysis.cxx` L44-299): parameter-period adjustment, shape end
//! bounds, 2D/3D contour cross products, outer-bound queries and face UV
//! bounds.
//!
//! Architecture bridges (the W1 edge.rs style):
//! 1. OCCT `ShapeExtend_WireData` -> `shape_extend::wire_data::WireData`
//!    (the W1-3 translation).
//! 2. OCCT `ShapeAnalysis_Edge` / `ShapeAnalysis_Curve` -> the W1-1 / W2
//!    1:1 translations in this module tree.
//! 3. OCCT `TopoDS_Shape` graph reads -> `rcad_kernel::topods::BRep` pool
//!    accessors; `TopoDS_Iterator(W, false)` -> the wire's raw edge list
//!    (brep_tool::raw_subshapes); `TopExp_Explorer` ->
//!    brep_tool::topexp_explorer; `TopExp::Vertices(W, V1, V2)` -> the
//!    local re-host below (TopExp.cxx L144-180 semantics).
//! 4. OCCT `Bnd_Box2d` -> the kernel `BndBox2d`; `BRep_Tool::Surface(F, L)`
//!    -> the TFace surface value + `Bounds()` -> `SurfaceEval::default_domain`.
//! 5. OCCT `BRepTopAdaptor_FClass2d(F, toluv)` -> the `FClass2dTopol`
//!    re-host (brep_top_adaptor/fclass2d_topol.rs).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::topo::topods::{BRep, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::math::bnd::BndBox2d;

use super::curve::ShapeAnalysisCurve;
use super::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_extend::wire_data::WireData;

// ---------------------------------------------------------------------------
// Free functions (the OCCT statics)
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis::AdjustByPeriod(Val, ToVal, Period)
/// (ShapeAnalysis.cxx L48-62).
pub fn adjust_by_period(the_val: f64, to_val: f64, period: f64) -> f64 {
    let diff = the_val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

/// OCCT ShapeAnalysis::AdjustToPeriod(Val, ValMin, ValMax)
/// (ShapeAnalysis.cxx L66-69).
pub fn adjust_to_period(the_val: f64, val_min: f64, val_max: f64) -> f64 {
    adjust_by_period(the_val, 0.5 * (val_min + val_max), val_max - val_min)
}

/// OCCT TopExp::Vertices(W, V1, V2, CumOri = true) over a wire (TopExp.cxx
/// L144-180): the traversal start/end vertices of the wire (bridge #3).
fn top_exp_wire_vertices(brep: &mut BRep, the_wire: &Shape) -> (Option<Shape>, Option<Shape>) {
    // TopoDS_Iterator(W, CumOri = true) — the orientation-composed edges.
    let edges = crate::shhealing::shape_build::brep_tool::iter_subshapes(
        brep,
        the_wire,
        true,
        false,
    );
    let Some(first_edge) = edges.first() else {
        return (None, None);
    };
    let Some(last_edge) = edges.last() else {
        return (None, None);
    };
    // The edge-level FirstVertex/LastVertex in traversal order.
    let v1 = top_exp_first_vertex(first_edge);
    let v2 = top_exp_last_vertex(last_edge);
    (v1, v2)
}

/// OCCT TopExp::FirstVertex(E, CumOri = true) (the loc_ope_wires_on_shape.rs
/// re-host).
fn top_exp_first_vertex(the_edge: &Shape) -> Option<Shape> {
    let ed = match the_edge.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    if the_edge.orientation.compose(ed.first.orientation) == Orientation::Forward {
        return Some(ed.first.clone());
    }
    if the_edge.orientation.compose(ed.last.orientation) == Orientation::Forward {
        return Some(ed.last.clone());
    }
    None
}

/// OCCT TopExp::LastVertex(E, CumOri = true).
fn top_exp_last_vertex(the_edge: &Shape) -> Option<Shape> {
    let ed = match the_edge.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    if the_edge.orientation.compose(ed.last.orientation) == Orientation::Reversed {
        return Some(ed.last.clone());
    }
    if the_edge.orientation.compose(ed.first.orientation) == Orientation::Reversed {
        return Some(ed.first.clone());
    }
    None
}

/// OCCT ShapeAnalysis::FindBounds(shape, V1, V2) (ShapeAnalysis.cxx L73-102).
pub fn find_bounds(brep: &mut BRep, shape: &Shape, v1: &mut Option<Shape>, v2: &mut Option<Shape>) {
    *v1 = None;
    *v2 = None;
    let ea = ShapeAnalysisEdge::new();
    if shape.shape_type() == ShapeType::Wire {
        let w = shape.clone();
        // invalid work with reversed wires replaced on TopExp
        let (a, b) = top_exp_wire_vertices(brep, &w);
        *v1 = a;
        *v2 = b;
    } else if shape.shape_type() == ShapeType::Edge {
        *v1 = Some(ea.first_vertex(brep, shape));
        *v2 = Some(ea.last_vertex(brep, shape));
    } else if shape.shape_type() == ShapeType::Vertex {
        *v1 = Some(shape.clone());
        *v2 = Some(shape.clone());
    }
}

/// OCCT ReverseSeq(Seq) (cxx L106-110) — the sequence reversal helper.
fn reverse_seq<T>(seq: &mut [T]) {
    seq.reverse();
}

/// OCCT BRep_Tool::CurveOnSurface(edge, face, f, l) — the stored (face)
/// overload matched by the face TShape pointer (the pcurve-key registry).
fn brep_tool_curve_on_surface_face(
    the_edge: &Shape,
    the_face: &Shape,
) -> Option<(Curve2d, f64, f64)> {
    let fptr = std::sync::Arc::as_ptr(&the_face.data) as u64;
    let reps = match the_edge.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return None,
    };
    for r in reps {
        let (key, pc, range) = match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        if key.0 == fptr {
            return Some((pc.clone(), range[0], range[1]));
        }
    }
    None
}

/// OCCT ShapeAnalysis::TotCross2D(sewd, aFace) (cxx L114-150) — the total
/// cross product of the wire's pcurve points.
pub fn tot_cross_2d(sewd: &WireData, a_face: &Shape) -> f64 {
    let mut nbc = 0i32;
    let mut fuv = DVec2::ZERO;
    let mut luv;
    let mut uv0 = DVec2::ZERO;
    let mut totcross = 0.0f64;
    for i in 1..=sewd.nb_edges() {
        let edge = sewd.edge(i);
        let Some((c2d, f2d, l2d)) = brep_tool_curve_on_surface_face(&edge, a_face) else {
            continue;
        };
        nbc += 1;
        let mut seq_pnt: Vec<DVec2> = Vec::new();
        ShapeAnalysisCurve::get_sample_points_2d(&c2d, f2d, l2d, &mut seq_pnt);
        if edge.orientation == Orientation::Reversed {
            // OCCT: edge.Orientation() == 1 (TopAbs_REVERSED).
            reverse_seq(&mut seq_pnt);
        }
        if nbc == 1 {
            fuv = seq_pnt[0];
            uv0 = fuv;
        }
        for pnt in &seq_pnt {
            luv = *pnt;
            totcross += (fuv.x - luv.x) * (fuv.y + luv.y) / 2.0;
            fuv = luv;
        }
    }
    totcross += (fuv.x - uv0.x) * (fuv.y + uv0.y) / 2.0;
    totcross
}

/// OCCT ShapeAnalysis::ContourArea(theWire) (cxx L154-199) — the contour
/// area from the 3D polygon of sample points.
pub fn contour_area(the_wire: &Shape) -> f64 {
    let mut nbc = 0i32;
    let mut fuv = DVec3::ZERO;
    let mut luv;
    let mut uv0 = DVec3::ZERO;
    let mut a_total = DVec3::ZERO;
    // TopoDS_Iterator aIte(theWire, false) — the raw (non-composed) edges.
    let edges = match the_wire.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    };
    for edge in &edges {
        // BRep_Tool::Curve(edge, first, last).
        let Some((c3d, first, last)) = (match edge.data.as_ref() {
            TShape::Edge(ed) => ed
                .curve
                .clone()
                .map(|c| (c, ed.range[0], ed.range[1])),
            _ => None,
        }) else {
            continue;
        };

        let mut a_seq_pnt: Vec<DVec3> = Vec::new();
        if !ShapeAnalysisCurve::get_sample_points(&c3d, first, last, &mut a_seq_pnt) {
            continue;
        }
        nbc += 1;
        if edge.orientation == Orientation::Reversed {
            reverse_seq(&mut a_seq_pnt);
        }
        if nbc == 1 {
            fuv = a_seq_pnt[0];
            uv0 = fuv;
        }
        for pnt in &a_seq_pnt {
            luv = *pnt;
            a_total += luv.cross(fuv); //
            fuv = luv;
        }
    }
    a_total += uv0.cross(fuv); //
    let an_area = a_total.length() * 0.5;
    an_area
}

/// OCCT ShapeAnalysis::IsOuterBound(face) (cxx L203-230).
pub fn is_outer_bound(brep: &mut BRep, face: &Shape) -> bool {
    // TopoDS_Face F = face; F.Orientation(TopAbs_FORWARD) — the FORWARD
    // orientation is the identity composition, so the wire children read
    // as stored.
    let mut w = Shape::null();
    let mut nbw = 0i32;
    let f_wires: Vec<Shape> = match face.data.as_ref() {
        TShape::Face(fd) => {
            let mut v = Vec::new();
            if !fd.outer_wire.is_null() {
                v.push(fd.outer_wire.clone());
            }
            v.extend(fd.inner_wires.iter().cloned());
            v
        }
        _ => Vec::new(),
    };
    for wire in &f_wires {
        w = wire.clone();
        nbw += 1;
    }
    // skl 08.04.2002
    if nbw == 1 {
        // ShapeExtend_WireData sewd = new ShapeExtend_WireData(W).
        let sewd = WireData::new_from_wire(brep, &w, false, true);
        let totcross = tot_cross_2d(&sewd, face);
        totcross >= 0.0
    } else {
        // BRepAdaptor_Surface Ads(F, false): the U/V resolutions of the face
        // surface (the kernel Surface3 re-hosts).
        let surf = match face.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        let Some(surf) = surf else {
            return false;
        };
        let tol = match face.data.as_ref() {
            TShape::Face(fd) => fd.tolerance,
            _ => 0.0,
        };
        let toluv = surf.u_resolution(tol).min(surf.v_resolution(tol));
        let f_class = crate::topalgo::brep_top_adaptor::fclass2d_topol::FClass2dTopol::new(
            std::sync::Arc::new(brep.clone()),
            face,
            toluv,
        );
        let rescl = f_class.perform_infinite_point() == rcad_kernel::topo::topods::State::Out;
        rescl
    }
}

/// OCCT ShapeAnalysis::OuterWire(theFace) (cxx L238-264) — returns the
/// positively oriented wire in the face; if there is no one, returns the
/// last wire of the face.
pub fn outer_wire(brep: &mut BRep, the_face: &Shape) -> Option<Shape> {
    // aF.Orientation(TopAbs_FORWARD) — the FORWARD orientation is the
    // identity composition, so the wire children read as stored.
    let an_it: Vec<Shape> = match the_face.data.as_ref() {
        TShape::Face(fd) => {
            let mut v = Vec::new();
            if !fd.outer_wire.is_null() {
                v.push(fd.outer_wire.clone());
            }
            v.extend(fd.inner_wires.iter().cloned());
            v
        }
        _ => Vec::new(),
    };
    let mut it = an_it.into_iter().peekable();
    while let Some(a_wire) = it.next() {
        // if current wire is the last one, return it without analysis
        if it.peek().is_none() {
            return Some(a_wire);
        }

        // Check if the wire has positive area
        // (ShapeExtend_WireData aSEWD = new ShapeExtend_WireData(aWire)).
        let a_sewd = WireData::new_from_wire(brep, &a_wire, false, true);
        let an_area2d = tot_cross_2d(&a_sewd, the_face);
        if an_area2d >= 0.0 {
            return Some(a_wire);
        }
    }
    None
}

/// OCCT ShapeAnalysis::GetFaceUVBounds(F, UMin, UMax, VMin, VMax)
/// (cxx L268-299).
pub fn get_face_uv_bounds(
    brep: &mut BRep,
    f: &Shape,
    u_min: &mut f64,
    u_max: &mut f64,
    v_min: &mut f64,
    v_max: &mut f64,
) {
    let ex = crate::shhealing::shape_build::brep_tool::topexp_explorer(brep, f, ShapeType::Edge);
    if ex.is_empty() {
        // BRep_Tool::Surface(F, L)->Bounds(UMin, UMax, VMin, VMax).
        if let TShape::Face(fd) = f.data.as_ref() {
            if let Some(surf) = fd.surface.as_ref() {
                let dom = SurfaceEval::default_domain(surf);
                *u_min = dom[0];
                *u_max = dom[1];
                *v_min = dom[2];
                *v_max = dom[3];
            }
        }
        return;
    }

    let mut b = BndBox2d::new();
    let mut sae = ShapeAnalysisEdge::new();
    let sac = ShapeAnalysisCurve;
    for edge in &ex {
        let mut c2d: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if !sae.pcurve_face(brep, edge, f, &mut c2d, &mut cf, &mut cl, false) {
            continue;
        }
        if let Some(c2d) = c2d.as_ref() {
            sac.fill_bnd_box(c2d, cf, cl, 20, true, &mut b);
        }
    }
    // B.Get(UMin, VMin, UMax, VMax).
    if let Some((xmin, ymin, xmax, ymax)) = b.get() {
        *u_min = xmin;
        *v_min = ymin;
        *u_max = xmax;
        *v_max = ymax;
    }
}
