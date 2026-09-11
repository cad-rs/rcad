//! OCCT BRepCheck_Vertex (TKTopAlgo/BRepCheck).
//!
//! Source: `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Vertex.cxx`
//! (L44-339) and `BRepCheck_Vertex.hxx` (L28-56).
//!
//! `Tolerance()` (Vertex.cxx L343-383) is not ported: it is not part of the
//! BRepCheck_Analyzer pass structure and its body is a walk over the vertex
//! point representations whose rcad identity mapping is pending (see the GAP
//! note in `in_context`).

use glam::DVec3;
use rcad_kernel::geom::CurveEval;
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape};

use super::brep_check_result::{
    brep_check_add, brep_tool_tolerance_edge, brep_tool_tolerance_face,
    brep_tool_tolerance_vertex, location_matrix, BRepCheckResultBase, BRepCheckStatus,
};

/// OCCT BRepCheck_Vertex (Vertex.hxx L28-56).
#[derive(Debug)]
pub struct BRepCheckVertex {
    /// OCCT protected base (Result.hxx L82-90).
    pub base: BRepCheckResultBase,
}

impl BRepCheckVertex {
    /// OCCT BRepCheck_Vertex::BRepCheck_Vertex(const TopoDS_Vertex& V)
    /// (Vertex.cxx L44-47).
    pub fn new(brep: &BRep, v: &Shape) -> Self {
        let mut r = BRepCheckVertex { base: BRepCheckResultBase::new() };
        r.base.init(v);
        r.minimum(brep);
        r
    }

    /// OCCT BRepCheck_Vertex::Minimum (Vertex.cxx L51-62) — checks the
    /// existence of a point 3D.
    pub fn minimum(&mut self, _brep: &BRep) {
        if !self.base.my_min {
            // OCCT L56-58: bind a fresh list for myShape.
            self.base.my_map.bound(&self.base.my_shape);
            let lst = self
                .base
                .my_map
                .find_mut(&self.base.my_shape)
                .expect("Minimum: myShape must be bound");
            // OCCT L59.
            lst.push(BRepCheckStatus::NoError);
            // OCCT L60.
            self.base.my_min = true;
        }
    }

    /// OCCT BRepCheck_Vertex::InContext (Vertex.cxx L66-285).
    pub fn in_context(&mut self, brep: &BRep, s: &Shape) {
        // OCCT L68-83: bound check under the (parallel) lock.
        if self.base.my_map.is_bound(s) {
            return;
        }
        self.base.my_map.bind(s.clone(), Vec::new());
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");

        // OCCT L86-98: try to find the vertex among the sub-vertices of S.
        {
            let mut found = false;
            let mut more = false;
            for cur in super::brep_check_result::explorer(brep, s, ShapeType::Vertex) {
                more = true;
                if cur.is_same(&self.base.my_shape) {
                    found = true;
                    break;
                }
            }
            if !more || !found {
                // OCCT L94-97: BRepCheck::Add(lst, BRepCheck_SubshapeNotInShape); return.
                brep_check_add(lst, BRepCheckStatus::SubshapeNotInShape);
                return; // leaves
            }
        }

        // OCCT L100-101: prep = TV->Pnt() (in the vertex frame; locations are
        // applied when comparing).
        let prep: DVec3 = match &*self.base.my_shape.data {
            TShape::Vertex(vd) => vd.point,
            _ => DVec3::ZERO,
        };
        // OCCT L102: Controlp (default-constructed gp_Pnt).
        #[allow(unused_assignments)]
        let mut controlp = DVec3::ZERO;

        let styp = s.shape_type();
        match styp {
            ShapeType::Edge => {
                // OCCT L107-140: try to find the vertex on the edge.
                let itv_children = {
                    // TopoDS_Iterator itv(E.Oriented(TopAbs_FORWARD)) — the
                    // forward-oriented edge enumerates its vertices with
                    // their own orientation.
                    let e_fwd = super::brep_check_result::oriented(s, Orientation::Forward);
                    super::brep_check_result::iterator_subshapes(brep, &e_fwd)
                };
                let mut vfind: Option<Shape> = None;
                let mut multiple = false;
                for itv_value in &itv_children {
                    if itv_value.is_same(&self.base.my_shape) {
                        match &vfind {
                            None => vfind = Some(itv_value.clone()),
                            Some(vf) => {
                                if (vf.orientation == Orientation::Forward
                                    && itv_value.orientation == Orientation::Reversed)
                                    || (vf.orientation == Orientation::Reversed
                                        && itv_value.orientation == Orientation::Forward)
                                {
                                    // the vertex on the edge is at once F and R
                                    multiple = true;
                                }
                                if vf.orientation != Orientation::Forward
                                    && vf.orientation != Orientation::Reversed
                                {
                                    if itv_value.orientation == Orientation::Forward
                                        || itv_value.orientation == Orientation::Reversed
                                    {
                                        vfind = Some(itv_value.clone());
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L142-147: VFind is not null for sure.
                let vfind = vfind.unwrap_or_else(|| self.base.my_shape.clone());
                let orv = vfind.orientation;

                let mut tol = brep_tool_tolerance_vertex(brep, &self.base.my_shape);
                tol = tol.max(brep_tool_tolerance_edge(brep, s)); // to check
                tol *= tol;

                // OCCT L149-242: the walk over the edge curve representations.
                let eloc = s.location;
                let reps = super::brep_check_result::edge_curve_reps(brep, s);
                for cr in &reps {
                    // OCCT L158-159: loc = cr->Location();
                    // L = (Eloc * loc).Predivided(myShape.Location()).
                    let cr_loc = match cr {
                        super::brep_check_result::EdgeCurveRep::Curve3D { location, .. } => *location,
                        super::brep_check_result::EdgeCurveRep::CurveOnSurface { location, .. } => *location,
                        super::brep_check_result::EdgeCurveRep::CurveOnClosedSurface { location, .. } => *location,
                        _ => 0,
                    };
                    let l_mat = location_matrix(brep, eloc, cr_loc)
                        * brep.get_location(self.base.my_shape.location).inverse();

                    match cr {
                        super::brep_check_result::EdgeCurveRep::Curve3D { curve, .. } => {
                            // OCCT L163-204: edge non degenerated.
                            // The TV->Points() walk (L166-180) matches point
                            // representations against the curve handle —
                            // GAP: rcad point representations carry pool
                            // indices without a handle-identity mapping, the
                            // pr loop is neutral (no status set).
                            //
                            // OCCT L181-202: the First/Last control.
                            if orv == Orientation::Forward || orv == Orientation::Reversed {
                                let ed = s.as_edge();
                                let first = ed.map(|e| e.range[0]).unwrap_or(0.0);
                                let last = ed.map(|e| e.range[1]).unwrap_or(0.0);
                                if orv == Orientation::Forward || multiple {
                                    // OCCT L186-192.
                                    controlp = curve.point_at(first);
                                    controlp = l_mat.transform_point3(controlp);
                                    if prep.distance_squared(controlp) > tol {
                                        brep_check_add(lst, BRepCheckStatus::InvalidPointOnCurve);
                                    }
                                }
                                if orv == Orientation::Reversed || multiple {
                                    // OCCT L193-201.
                                    controlp = curve.point_at(last);
                                    controlp = l_mat.transform_point3(controlp);
                                    if prep.distance_squared(controlp) > tol {
                                        brep_check_add(lst, BRepCheckStatus::InvalidPointOnCurve);
                                    }
                                }
                            }
                        }
                        super::brep_check_result::EdgeCurveRep::CurveOnSurface { .. } => {
                            // OCCT L205-240: the pcurve/surface branch —
                            // Controlp is evaluated from the vertex point
                            // representations on the pcurve (L214-239);
                            // GAP: the rcad point-representation identity
                            // mapping is pending — neutral (no status set).
                        }
                        _ => {
                            // Regularity representations are neither a 3D
                            // curve nor a curve on surface (OCCT L161/L205
                            // both false).
                        }
                    }
                }
                // OCCT L243-246.
                if lst.is_empty() {
                    lst.push(BRepCheckStatus::NoError);
                }
            }
            ShapeType::Face => {
                // OCCT L249-258: L = (Floc * TFloc).Predivided(myShape.Location()).
                let floc = s.location;
                let tfloc = s
                    .as_face()
                    .map(|fd| fd.surface_location)
                    .unwrap_or(0);
                let l_mat = location_matrix(brep, floc, tfloc)
                    * brep.get_location(self.base.my_shape.location).inverse();

                #[allow(unused_assignments)]
                let mut tol = brep_tool_tolerance_vertex(brep, &self.base.my_shape);
                tol = tol.max(brep_tool_tolerance_face(brep, s)); // to check
                tol *= tol;
                let _ = tol;

                // OCCT L260-274: the pr->IsPointOnSurface(Su, L) walk —
                // GAP: the rcad point-representation identity mapping is
                // pending — neutral (no status set).
                let _ = l_mat;
                // OCCT L275-278.
                if lst.is_empty() {
                    lst.push(BRepCheckStatus::NoError);
                }
            }
            _ => {
                // OCCT L281-283: default — nothing.
            }
        }
    }

    /// OCCT BRepCheck_Vertex::Blind (Vertex.cxx L289-339) — the body was
    /// removed upstream because of its uselessness.
    pub fn blind(&mut self) {
        if self.base.my_blind {
            return;
        }
        self.base.my_blind = true;
    }
}
