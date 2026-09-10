//! OCCT ChFi2d_FilletAPI (FilletAPI.hxx L20-88 + FilletAPI.cxx L28-141) —
//! the interface class choosing between the two 2D-fillet algorithms.
//!
//! OCCT provides two algorithms for 2D fillets: `ChFi2d_Builder` (fillet or
//! chamfer for linear and circular edges of a face) and this `ChFi2d_FilletAPI`
//! encapsulation of `ChFi2d_AnaFilletAlgo` (analytic, linear/circular edges
//! with a common point) and `ChFi2d_FilletAlgo` (iteration-recursive, any
//! edge type, common point not required).  This class chooses an algorithm
//! by analyzing the arguments (`IsAnalytical`).
//!
//! Architecture notes (Rust vs OCCT):
//!   - OCCT `Init` overloads -> `init_wire` / `init_edges` (no overloading).
//!   - OCCT default argument `iSolution = -1` -> explicit parameter; pass -1
//!     for the OCCT default.
//!   - `gp_Pln` -> [`rcad_kernel::geom::Plane`]; `gp_Pnt` -> `DVec3`;
//!     `TopoDS_Edge`/`TopoDS_Wire` -> `Shape`.
//!   - `BRepAdaptor_Curve::GetType` -> the edge's `Curve3` variant match
//!     (`Line`/`Circle` cover `GeomAbs_Line`/`GeomAbs_Circle`).
//!
//! Joint-compile note: the callee method signatures follow the OCCT hxx +
//! the rcad naming convention agreed with the parallel Stage 1a agents
//! (ChFi2d_AnaFilletAlgo / ChFi2d_FilletAlgo); reconcile at the joint
//! compile once those files land.

use glam::DVec3;
use rcad_kernel::geom::{CurveEval as _, Curve3, Plane};
use rcad_kernel::topo::topods::Shape;

use super::chfi2d_ana_fillet_algo::ChFi2dAnaFilletAlgo;
use super::chfi2d_fillet_algo::ChFi2dFilletAlgo;

/// OCCT ChFi2d_FilletAPI (FilletAPI.hxx L46-86).
#[derive(Debug, Clone)]
pub struct ChFi2dFilletAPI {
    /// OCCT: ChFi2d_FilletAlgo myFilletAlgo.
    my_fillet_algo: ChFi2dFilletAlgo,
    /// OCCT: ChFi2d_AnaFilletAlgo myAnaFilletAlgo.
    my_ana_fillet_algo: ChFi2dAnaFilletAlgo,
    /// OCCT: bool myIsAnalytical.
    my_is_analytical: bool,
}

impl Default for ChFi2dFilletAPI {
    /// OCCT FilletAPI.cxx L30-34 — the empty constructor; call `init_*`
    /// before `perform`.
    fn default() -> Self {
        ChFi2dFilletAPI {
            my_fillet_algo: ChFi2dFilletAlgo::new(),
            my_ana_fillet_algo: ChFi2dAnaFilletAlgo::new(),
            my_is_analytical: false,
        }
    }
}

impl ChFi2dFilletAPI {
    /// OCCT FilletAPI.cxx L30-34 — ChFi2d_FilletAPI().
    pub fn new() -> Self {
        Default::default()
    }

    /// OCCT FilletAPI.cxx L36-43 — the wire constructor: a wire consisting
    /// of two edges in a plane.
    pub fn new_wire(the_wire: &Shape, the_plane: &Plane) -> Self {
        let mut api = ChFi2dFilletAPI::new();
        api.init_wire(the_wire, the_plane);
        api
    }

    /// OCCT FilletAPI.cxx L45-53 — the two-edges constructor.
    pub fn new_edges(the_edge1: &Shape, the_edge2: &Shape, the_plane: &Plane) -> Self {
        let mut api = ChFi2dFilletAPI::new();
        api.init_edges(the_edge1, the_edge2, the_plane);
        api
    }

    /// OCCT FilletAPI.cxx L55-80 — Init(theWire, thePlane): take the first
    /// two edges of the wire, decide analytical vs iteration-recursive,
    /// forward to the chosen algorithm.
    pub fn init_wire(&mut self, the_wire: &Shape, the_plane: &Plane) {
        // Decide whether we may apply an analytical solution.
        let mut e1: Option<Shape> = None;
        let mut e2: Option<Shape> = None;
        // TopoDS_Iterator walk over the wire's sub-edges (TWireData.edges):
        // first edge -> E1, second -> E2, else break.
        if let Some(wire) = the_wire.as_wire() {
            for etr in &wire.edges {
                if e1.is_none() {
                    e1 = Some(etr.clone());
                } else if e2.is_none() {
                    e2 = Some(etr.clone());
                } else {
                    break;
                }
            }
        }
        if let (Some(e1), Some(e2)) = (&e1, &e2) {
            self.my_is_analytical = self.is_analytical(e1, e2);
        }

        // Initialize the algorithm.
        if self.my_is_analytical {
            self.my_ana_fillet_algo.init_wire(the_wire, the_plane);
        } else {
            self.my_fillet_algo.init_wire(the_wire, the_plane);
        }
    }

    /// OCCT FilletAPI.cxx L82-92 — Init(theEdge1, theEdge2, thePlane).
    pub fn init_edges(&mut self, the_edge1: &Shape, the_edge2: &Shape, the_plane: &Plane) {
        // Decide whether we may apply an analytical solution.
        self.my_is_analytical = self.is_analytical(the_edge1, the_edge2);

        // Initialize the algorithm.
        if self.my_is_analytical {
            self.my_ana_fillet_algo.init_edges(the_edge1, the_edge2, the_plane);
        } else {
            self.my_fillet_algo.init_edges(the_edge1, the_edge2, the_plane);
        }
    }

    /// OCCT FilletAPI.cxx L94-99 — Perform(theRadius): constructs a fillet
    /// edge; returns true if at least one result was found.
    pub fn perform(&mut self, the_radius: f64) -> bool {
        if self.my_is_analytical {
            self.my_ana_fillet_algo.perform(the_radius)
        } else {
            self.my_fillet_algo.perform(the_radius)
        }
    }

    /// OCCT FilletAPI.cxx L101-105 — NbResults(thePoint): number of possible
    /// solutions; thePoint chooses a particular fillet in case of several
    /// (put the intersecting / common point of the edges).  OCCT declares
    /// NbResults non-const on the member algorithms -> `&mut self`.
    pub fn nb_results(&mut self, the_point: DVec3) -> i32 {
        if self.my_is_analytical {
            1
        } else {
            self.my_fillet_algo.nb_results(the_point)
        }
    }

    /// OCCT FilletAPI.cxx L107-116 — Result(thePoint, theEdge1, theEdge2,
    /// iSolution = -1): returns the fillet edge (with the modified edge1 /
    /// edge2 written back) nearest to thePoint when iSolution == -1.
    /// OCCT declares Result non-const -> `&mut self`.
    pub fn result(
        &mut self,
        the_point: DVec3,
        the_edge1: &mut Shape,
        the_edge2: &mut Shape,
        i_solution: i32,
    ) -> Shape {
        if self.my_is_analytical {
            self.my_ana_fillet_algo.result(the_edge1, the_edge2)
        } else {
            self.my_fillet_algo
                .result(the_point, the_edge1, the_edge2, i_solution)
        }
    }

    /// OCCT FilletAPI.cxx L118-141 — IsAnalytical: the analytical solution
    /// is applicable for linear and circular edges having a common point.
    fn is_analytical(&self, the_edge1: &Shape, the_edge2: &Shape) -> bool {
        let mut ret = false;
        let ed1 = the_edge1.as_edge();
        let ed2 = the_edge2.as_edge();
        if let (Some(ed1), Some(ed2)) = (ed1, ed2) {
            let c1 = ed1.curve.as_ref();
            let c2 = ed2.curve.as_ref();
            if let (Some(c1), Some(c2)) = (c1, c2) {
                let is_line_or_circle = |c: &Curve3| {
                    matches!(c, Curve3::Line(_) | Curve3::Circle(_))
                };
                if is_line_or_circle(c1) && is_line_or_circle(c2) {
                    // The edges are lines or arcs of circle.
                    // Now check whether they have a common point.
                    let p11 = c1.point_at(ed1.range[0]);
                    let p12 = c1.point_at(ed1.range[1]);
                    let p21 = c2.point_at(ed2.range[0]);
                    let p22 = c2.point_at(ed2.range[1]);
                    let sq_conf = rcad_kernel::core::precision::CONFUSION
                        * rcad_kernel::core::precision::CONFUSION;
                    if p11.distance_squared(p21) < sq_conf
                        || p11.distance_squared(p22) < sq_conf
                        || p12.distance_squared(p21) < sq_conf
                        || p12.distance_squared(p22) < sq_conf
                    {
                        ret = true;
                    }
                }
            }
        }
        ret
    }
}
