//! ChFi2d_ChamferAPI — OCCT TKFillet 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi2d/
//!   - ChFi2d_ChamferAPI.hxx (L24-57)
//!   - ChFi2d_ChamferAPI.cxx (L24-160)
//!
//! Architecture differences (Rust vs C++ data model):
//!   - TopoDS_Edge / TopoDS_Wire -> `rcad_kernel::topo::topo_shape::Shape`
//!     handles.  OCCT BRep_Tool::Curve reads the 3D curve + range from the
//!     edge TShape (with the shape Location applied); rcad reads the
//!     `Arc<TShape>` payload directly (curve + range), which is identical
//!     for the identity-location shapes this sketcher-level API operates
//!     on.  Location resolution would require the owning BRep table and is
//!     left to the caller (recorded limitation).
//!   - `occ::handle<Geom_Curve>` -> `Curve3` value (Option for the null
//!     pre-Perform state).
//!   - BRepBuilderAPI_MakeEdge / GC_MakeLine have no standalone rcad
//!     translation yet; their behavior is inlined here: the created edges
//!     are new TShapes appended to `my_brep` (an owned empty BRep table),
//!     and the line is built by `rcad_kernel::base::gc::make_line_2p`.

use rcad_kernel::base::gc::make_line_2p;
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve3, CurveEval as _, Line3, Point3};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods;

/// OCCT ChFi2d_ChamferAPI (ChFi2d_ChamferAPI.hxx L24-57) — a class making
/// a chamfer between two linear edges.
#[derive(Debug, Clone)]
pub struct ChFi2dChamferAPI {
    /// OCCT: TopoDS_Edge myEdge1, myEdge2.
    pub my_edge1: Shape,
    pub my_edge2: Shape,
    /// OCCT: occ::handle<Geom_Curve> myCurve1, myCurve2.
    pub my_curve1: Option<Curve3>,
    pub my_curve2: Option<Curve3>,
    /// OCCT: double myStart1, myEnd1, myStart2, myEnd2.
    pub my_start1: f64,
    pub my_end1: f64,
    pub my_start2: f64,
    pub my_end2: f64,
    /// OCCT: bool myCommonStart1, myCommonStart2.
    pub my_common_start1: bool,
    pub my_common_start2: bool,
    /// rcad architecture: the chamfer edge and the trimmed result edges are
    /// created into a BRep TShape table (OCCT: BRepBuilderAPI_MakeEdge
    /// creates them in the global handle graph).  The algorithm owns one.
    pub my_brep: topods::BRep,
}

impl ChFi2dChamferAPI {
    /// OCCT ChFi2d_ChamferAPI.cxx L24-32 — an empty constructor.
    pub fn new() -> Self {
        ChFi2dChamferAPI {
            my_edge1: Shape::null(),
            my_edge2: Shape::null(),
            my_curve1: None,
            my_curve2: None,
            my_start1: 0.0,
            my_end1: 0.0,
            my_start2: 0.0,
            my_end2: 0.0,
            my_common_start1: false,
            my_common_start2: false,
            my_brep: topods::BRep::new(),
        }
    }

    /// OCCT ChFi2d_ChamferAPI.cxx L35-44 — a constructor accepting a wire
    /// consisting of two linear edges.
    pub fn new_wire(the_wire: &Shape) -> Self {
        let mut an_api = ChFi2dChamferAPI::new();
        an_api.init_wire(the_wire);
        an_api
    }

    /// OCCT ChFi2d_ChamferAPI.cxx L47-57 — a constructor accepting two
    /// linear edges (myEdge1 / myEdge2 set in the initializer list).
    pub fn new_edges(the_edge1: &Shape, the_edge2: &Shape) -> Self {
        let mut an_api = ChFi2dChamferAPI::new();
        an_api.my_edge1 = the_edge1.clone();
        an_api.my_edge2 = the_edge2.clone();
        an_api
    }

    /// OCCT ChFi2d_ChamferAPI.cxx L60-80 — initializes the class by a wire
    /// consisting of two linear edges.
    pub fn init_wire(&mut self, the_wire: &Shape) {
        let mut e1 = Shape::null();
        let mut e2 = Shape::null();
        // OCCT TopoDS_Iterator itr(theWire); for (; itr.More(); itr.Next())
        // Architecture bridge: the wire children are the TWireData::edges
        // list read through the Arc<TShape> payload.
        let a_children = the_wire.as_wire().expect("TopoDS::Wire").edges.clone();
        for a_shape in &a_children {
            if e1.is_null() {
                // OCCT: E1 = TopoDS::Edge(itr.Value());
                e1 = a_shape.clone();
            } else if e2.is_null() {
                // OCCT: E2 = TopoDS::Edge(itr.Value());
                e2 = a_shape.clone();
            } else {
                break;
            }
        }
        self.init(&e1, &e2);
    }

    /// OCCT ChFi2d_ChamferAPI.cxx L83-87 — initializes the class by two
    /// linear edges.
    pub fn init(&mut self, the_edge1: &Shape, the_edge2: &Shape) {
        self.my_edge1 = the_edge1.clone();
        self.my_edge2 = the_edge2.clone();
    }

    /// OCCT ChFi2d_ChamferAPI.cxx L91-123 — constructs a chamfer edge.
    /// Returns true if the edge is constructed.
    pub fn perform(&mut self) -> bool {
        // OCCT L93: myCurve1 = BRep_Tool::Curve(myEdge1, myStart1, myEnd1);
        let a_ed1 = self
            .my_edge1
            .as_edge()
            .expect("ChFi2d_ChamferAPI::Perform - myEdge1 is not an edge");
        let a_curve1 = a_ed1
            .curve
            .clone()
            .expect("ChFi2d_ChamferAPI::Perform - myEdge1 has no 3D curve");
        self.my_start1 = a_ed1.range[0];
        self.my_end1 = a_ed1.range[1];
        self.my_curve1 = Some(a_curve1);
        // OCCT L94: myCurve2 = BRep_Tool::Curve(myEdge2, myStart2, myEnd2);
        let a_ed2 = self
            .my_edge2
            .as_edge()
            .expect("ChFi2d_ChamferAPI::Perform - myEdge2 is not an edge");
        let a_curve2 = a_ed2
            .curve
            .clone()
            .expect("ChFi2d_ChamferAPI::Perform - myEdge2 has no 3D curve");
        self.my_start2 = a_ed2.range[0];
        self.my_end2 = a_ed2.range[1];
        self.my_curve2 = Some(a_curve2);
        // OCCT reads the curve handles below; rcad re-reads the stored values.
        let a_curve1 = self.my_curve1.clone().expect("myCurve1");
        let a_curve2 = self.my_curve2.clone().expect("myCurve2");

        // searching for common points
        // OCCT gp_Pnt::IsEqual(theOther, Prec) == (Distance(theOther) <= Prec).
        if a_curve1.point_at(self.my_start1).distance(a_curve2.point_at(self.my_end2)) <= CONFUSION
        {
            self.my_common_start1 = true;
            self.my_common_start2 = false;
        } else {
            if a_curve1.point_at(self.my_end1).distance(a_curve2.point_at(self.my_start2))
                <= CONFUSION
            {
                self.my_common_start1 = false;
                self.my_common_start2 = true;
            } else {
                if a_curve1.point_at(self.my_end1).distance(a_curve2.point_at(self.my_end2))
                    <= CONFUSION
                {
                    self.my_common_start1 = false;
                    self.my_common_start2 = false;
                } else {
                    self.my_common_start1 = true;
                    self.my_common_start2 = true;
                }
            }
        }
        true
    }

    /// OCCT ChFi2d_ChamferAPI.cxx L126-160 — returns the result (chamfer
    /// edge, modified edge1, modified edge2).
    pub fn result(
        &mut self,
        the_edge1: &mut Shape,
        the_edge2: &mut Shape,
        the_length1: f64,
        the_length2: f64,
    ) -> Shape {
        let mut a_result = Shape::null();
        if (self.my_end1 - self.my_start1).abs() < the_length1 {
            return a_result;
        }
        if (self.my_end2 - self.my_start2).abs() < the_length2 {
            return a_result;
        }

        // OCCT L141-144
        let a_common1 = (if self.my_common_start1 {
            self.my_start1
        } else {
            self.my_end1
        }) + (if (self.my_start1 > self.my_end1) ^ self.my_common_start1 {
            the_length1
        } else {
            -the_length1
        });
        let a_common2 = (if self.my_common_start2 {
            self.my_start2
        } else {
            self.my_end2
        }) + (if (self.my_start2 > self.my_end2) ^ self.my_common_start2 {
            the_length2
        } else {
            -the_length2
        });

        let a_curve1 = self.my_curve1.clone().expect("myCurve1");
        let a_curve2 = self.my_curve2.clone().expect("myCurve2");

        // make chamfer edge
        // OCCT L147-151: GC_MakeLine aML(myCurve1->Value(aCommon1),
        // myCurve2->Value(aCommon2)); BRepBuilderAPI_MakeEdge aBuilder(aML.Value(), P1, P2);
        let a_p1: Point3 = a_curve1.point_at(a_common1);
        let a_p2: Point3 = a_curve2.point_at(a_common2);
        let a_line: Line3 = make_line_2p(a_p1, a_p2).expect("GC_MakeLine - aML.Value()");
        a_result = self.make_edge_curve_p1p2(Curve3::Line(a_line), a_p1, a_p2);

        // divide first edge
        // OCCT L153-154: BRepBuilderAPI_MakeEdge aDivider1(myCurve1, aCommon1,
        // (myCommonStart1 ? myEnd1 : myStart1)); theEdge1 = aDivider1.Edge();
        let a_divider1 = self.make_edge_curve_range(
            a_curve1.clone(),
            a_common1,
            if self.my_common_start1 {
                self.my_end1
            } else {
                self.my_start1
            },
        );
        *the_edge1 = a_divider1;
        // divide second edge
        // OCCT L156-157: BRepBuilderAPI_MakeEdge aDivider2(myCurve2, aCommon2,
        // (myCommonStart2 ? myEnd2 : myStart2)); theEdge2 = aDivider2.Edge();
        let a_divider2 = self.make_edge_curve_range(
            a_curve2.clone(),
            a_common2,
            if self.my_common_start2 {
                self.my_end2
            } else {
                self.my_start2
            },
        );
        *the_edge2 = a_divider2;

        a_result
    }

    /// OCCT BRepBuilderAPI_MakeEdge(Curve, P1, P2) (BRepLib_MakeEdge.cxx) —
    /// builds the edge on theCurve between the projections of P1 and P2.
    /// Architecture bridge: rcad has no standalone BRepBuilderAPI_MakeEdge
    /// translation; the edge TShape is created in `my_brep` with vertices at
    /// the two points (parameters obtained by projection onto the curve).
    fn make_edge_curve_p1p2(&mut self, the_curve: Curve3, the_p1: Point3, the_p2: Point3) -> Shape {
        let a_param1 = project_parameter(&the_curve, the_p1);
        let a_param2 = project_parameter(&the_curve, the_p2);
        let a_v1 = self.my_brep.add_tvertex(the_p1);
        let a_v2 = self.my_brep.add_tvertex(the_p2);
        self.my_brep
            .add_tedge(Some(the_curve), a_v1, a_v2, [a_param1, a_param2])
    }

    /// OCCT BRepBuilderAPI_MakeEdge(Curve, p1, p2) — builds the edge on the
    /// curve restricted to the parameter range [p1, p2]; vertices are created
    /// at the curve points of the two parameters.
    fn make_edge_curve_range(&mut self, the_curve: Curve3, the_p1: f64, the_p2: f64) -> Shape {
        let a_v1 = self.my_brep.add_tvertex(the_curve.point_at(the_p1));
        let a_v2 = self.my_brep.add_tvertex(the_curve.point_at(the_p2));
        self.my_brep
            .add_tedge(Some(the_curve), a_v1, a_v2, [the_p1, the_p2])
    }
}

/// OCCT BRepLib_MakeEdge parameter projection — the vertex parameter of a
/// point on the curve.  For a line this is ElCLib::LineParameter
/// ((P - Location) . Direction); other curves fall back to the kernel
/// closest-point projection.
fn project_parameter(the_curve: &Curve3, the_point: Point3) -> f64 {
    match the_curve {
        Curve3::Line(a_line) => (the_point - a_line.origin).dot(a_line.direction),
        _ => {
            let a_domain = the_curve.default_domain();
            rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
                the_curve,
                the_point,
                a_domain[0],
                a_domain[1],
                64,
            )
            .param
        }
    }
}
