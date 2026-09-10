//! 1:1 translation of OCCT `ShapeFix_FreeBounds`
//! (`TKShHealing/ShapeFix/ShapeFix_FreeBounds.hxx` L17-101 +
//! `ShapeFix_FreeBounds.cxx` L1-124 + `ShapeFix_FreeBounds.lxx` L1-36,
//! docket row `ShapeFix_FreeBounds`).
//!
//! Function-count equation (OCCT ShapeFix_FreeBounds.cxx + lxx = rcad
//! free_bounds.rs): `ShapeFix_FreeBounds()` / `ShapeFix_FreeBounds(shape,
//! sewtoler, closetoler, splitclosed, splitopen)` /
//! `ShapeFix_FreeBounds(shape, closetoler, splitclosed, splitopen)` /
//! `Perform` (cxx L36-124) + `GetClosedWires` / `GetOpenWires` / `GetShape`
//! (lxx) — 7 OCCT functions = 7 rcad functions.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — the `BRep_Builder`/`TopExp_Explorer` re-hosts
//!    over `rcad_kernel::BRep` (the package precedents).
//! 2. `TopoDS_Compound` fields — the `Shape` wrappers over the compounds.
//! 3. `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>` — the
//!    `HashMap<(u64, u32), Shape>` keyed by (TShape ptr, location index)
//!    (the TopTools_ShapeMapHasher IsSame identity).
//! 4. The `TopoDS_Iterator` walk over the edge vertices inside the
//!    replacement loop is snapshotted (the OCCT iterator walks the live
//!    TShape list that the loop mutates; the appended new vertices are
//!    never bound in the map, so the net replacement set is identical).

use std::collections::HashMap;

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, ShapeType};

use crate::shhealing::shape_analysis::free_bounds::ShapeAnalysisFreeBounds;
use crate::shhealing::shape_build::brep_tool::{iter_subshapes, topexp_explorer};
use crate::shhealing::shape_extend::explorer::ShapeExtendExplorer;

/// The map key of a shape (TShape pointer + location index).
#[inline]
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT ShapeFix_FreeBounds (hxx L49-97): outputs the free bounds of the
/// shape (wires of edges referenced by the only face), complementing
/// ShapeAnalysis_FreeBounds with the reduction of the number of open wires.
pub struct ShapeFixFreeBounds {
    /// OCCT myWires (hxx L89): compound of closed wires out of free edges.
    my_wires: Shape,
    /// OCCT myEdges (hxx L90): compound of open wires out of free edges.
    my_edges: Shape,
    /// OCCT myShape (hxx L91): the source shape (possibly modified).
    my_shape: Shape,
    /// OCCT myShared (hxx L92).
    my_shared: bool,
    /// OCCT mySewToler (hxx L93).
    my_sew_toler: f64,
    /// OCCT myCloseToler (hxx L94).
    my_close_toler: f64,
    /// OCCT mySplitClosed (hxx L95).
    my_split_closed: bool,
    /// OCCT mySplitOpen (hxx L96).
    my_split_open: bool,
}

impl Default for ShapeFixFreeBounds {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixFreeBounds {
    /// OCCT ShapeFix_FreeBounds::ShapeFix_FreeBounds() (cxx L36-43): the
    /// empty constructor.
    pub fn new() -> Self {
        ShapeFixFreeBounds {
            my_wires: Shape::null(),
            my_edges: Shape::null(),
            my_shape: Shape::null(),
            my_shared: false,
            my_sew_toler: 0.0,
            my_close_toler: 0.0,
            my_split_closed: false,
            my_split_open: false,
        }
    }

    /// OCCT ShapeFix_FreeBounds::ShapeFix_FreeBounds(shape, sewtoler,
    /// closetoler, splitclosed, splitopen) (cxx L47-60): builds the
    /// forecasting free bounds of the shape (a compound of faces) and
    /// connects the open wires with the tolerance `closetoler`.
    pub fn new_forecast(
        brep: &mut BRep,
        shape: &Shape,
        sewtoler: f64,
        closetoler: f64,
        splitclosed: bool,
        splitopen: bool,
    ) -> Self {
        let mut this = ShapeFixFreeBounds {
            my_wires: Shape::null(),
            my_edges: Shape::null(),
            my_shape: Shape::null(),
            my_shared: false,
            my_sew_toler: sewtoler,
            my_close_toler: closetoler,
            my_split_closed: splitclosed,
            my_split_open: splitopen,
        };
        // L58-59.
        this.my_shape = shape.clone();
        this.perform(brep);
        this
    }

    /// OCCT ShapeFix_FreeBounds::ShapeFix_FreeBounds(shape, closetoler,
    /// splitclosed, splitopen) (cxx L64-76): builds the actual free bounds
    /// of the shape (a compound of shells).
    pub fn new_actual(
        brep: &mut BRep,
        shape: &Shape,
        closetoler: f64,
        splitclosed: bool,
        splitopen: bool,
    ) -> Self {
        let mut this = ShapeFixFreeBounds {
            my_wires: Shape::null(),
            my_edges: Shape::null(),
            my_shape: Shape::null(),
            my_shared: true,
            my_sew_toler: 0.0,
            my_close_toler: closetoler,
            my_split_closed: splitclosed,
            my_split_open: splitopen,
        };
        // L74-75.
        this.my_shape = shape.clone();
        this.perform(brep);
        this
    }

    /// OCCT ShapeFix_FreeBounds::Perform (cxx L80-124).
    pub fn perform(&mut self, brep: &mut BRep) -> bool {
        // L82-90: the analyzer over the shape.
        let safb = if self.my_shared {
            // L85: ShapeAnalysis_FreeBounds(myShape, mySplitClosed,
            // mySplitOpen) — the shared form (CheckInternaledges defaults
            // to False).
            ShapeAnalysisFreeBounds::new_actual(
                brep,
                &self.my_shape,
                self.my_split_closed,
                self.my_split_open,
                false,
            )
        } else {
            // L89: ShapeAnalysis_FreeBounds(myShape, mySewToler,
            // mySplitClosed, mySplitOpen) — the sewing form.
            ShapeAnalysisFreeBounds::new_forecast(
                brep,
                &self.my_shape,
                self.my_sew_toler,
                self.my_split_closed,
                self.my_split_open,
            )
        };

        // L92-93.
        self.my_wires = safb.get_closed_wires().clone();
        self.my_edges = safb.get_open_wires().clone();

        // L95-122.
        if self.my_close_toler > self.my_sew_toler {
            // L97-98: see.SeqFromCompound(myEdges, false).
            let see = ShapeExtendExplorer::new();
            let mut open = see.seq_from_compound(brep, &self.my_edges, false);
            // L99: the map of the original to new connecting vertices.
            let mut vertices: HashMap<(u64, u32), Shape> = HashMap::new();
            // L100-101: newwires = ConnectWiresToWires(open, myCloseToler,
            // myShared, vertices).
            let newwires = ShapeAnalysisFreeBounds::connect_wires_to_wires_vertices(
                brep,
                &mut open,
                self.my_close_toler,
                self.my_shared,
                &mut vertices,
            );
            // L102: myEdges.Nullify().
            self.my_edges = Shape::null();
            // L103: DispatchWires(newwires, myWires, myEdges).
            let mut my_wires = self.my_wires.clone();
            let mut my_edges = self.my_edges.clone();
            ShapeAnalysisFreeBounds::dispatch_wires(
                brep,
                Some(&newwires),
                &mut my_wires,
                &mut my_edges,
            );
            self.my_wires = my_wires;
            self.my_edges = my_edges;

            // L105-121: the shape edges are updated with the new connecting
            // vertices.
            for exp in topexp_explorer(brep, &self.my_shape, ShapeType::Edge) {
                let edge = exp;
                // L108-120: the child vertices (the snapshot bridge #4).
                let children = iter_subshapes(brep, &edge, false, true);
                for child in children {
                    let v = child;
                    // L113-115: vertices.IsBound(V) / Find(V).
                    let Some(new_v) = vertices.get(&shape_key(&v)) else {
                        continue;
                    };
                    let mut b = BRepBuilder::new();
                    let mut new_v = new_v.clone();
                    // L116.
                    new_v.orientation = v.orientation;
                    // L117-118.
                    b.remove_from_edge(brep, edge.clone(), v);
                    b.add_to_edge(brep, edge.clone(), new_v);
                }
            }
        }
        // L123.
        true
    }

    /// OCCT ShapeFix_FreeBounds::GetClosedWires (lxx L15-18): the compound
    /// of closed wires out of free edges.
    pub fn get_closed_wires(&self) -> &Shape {
        &self.my_wires
    }

    /// OCCT ShapeFix_FreeBounds::GetOpenWires (lxx L23-26): the compound of
    /// open wires out of free edges.
    pub fn get_open_wires(&self) -> &Shape {
        &self.my_edges
    }

    /// OCCT ShapeFix_FreeBounds::GetShape (lxx L31-34): the (possibly
    /// modified) source shape.
    pub fn get_shape(&self) -> &Shape {
        &self.my_shape
    }
}
