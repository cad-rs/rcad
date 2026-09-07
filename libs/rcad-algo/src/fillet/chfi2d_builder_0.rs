//! OCCT ChFi2d_Builder_0 — 1:1 translation.
//!
//! Source: ChFi2d_Builder_0.cxx L60-890 (TKFillet/ChFi2d) — the chamfer
//! half of ChFi2d_Builder (AddChamfer / ComputeChamfer / ModifyChamfer /
//! RemoveChamfer / BuildChamferEdge) plus the file-level helpers
//! (ComputePoint x2, OrientChamfer, IsLineOrCircle).
//!
//! Rust overload disambiguation (OCCT C++ overloads):
//!   - AddChamfer(E1,E2,D1,D2)          -> add_chamfer_edges
//!   - AddChamfer(E,V,D,Ang)            -> add_chamfer_edge_vertex
//!   - ComputeChamfer(V,E1,E2,D1,D2,..) -> compute_chamfer_edges
//!   - ComputeChamfer(V,E1,D,Ang,E2,..) -> compute_chamfer_angle
//!   - ModifyChamfer(Ch,E1,E2,D1,D2)    -> modify_chamfer_edges
//!   - ModifyChamfer(Ch,E,D,Ang)        -> modify_chamfer_angle
//!   - BuildChamferEdge(V,E1,E2,D1,D2,..)     -> build_chamfer_edge_edges
//!   - BuildChamferEdge(V,E1,D,Ang,E2,..)     -> build_chamfer_edge_angle

use glam::{DVec2, DVec3, DQuat};
use rcad_kernel::base::geom_proj_lib::project_on_plane::curve_on_plane;
use rcad_kernel::base::gcpnts::abscissa_point::abscissa_point_parameter;
use rcad_kernel::core::precision::{CONFUSION, p_intersection};
use rcad_kernel::geom::{
    Curve2d, Curve3, CurveEval as _, Line3, Plane, Surface3, SurfaceEval as _,
};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, BRepTool as _, Orientation};

use crate::geomalgo::geom2d_int::GInter;

// OCCT ChFi2d.cxx package functions (translated in fillet/chfi2d.rs).
use super::chfi2d::{
    chfi2d_common_vertex, chfi2d_find_connected_edges, ChFi2dConstructionError,
};
use super::chfi2d_builder::{brep_builder_make_vertex, topexp_explore_face_edges};

// ===========================================================================
// File-level helpers (OCCT ChFi2d_Builder_0.cxx)
// ===========================================================================

/// OCCT ElCLib::Parameter(gp_Lin, gp_Pnt) (ElCLib.cxx L1267-1273): the
/// parameter of the point in the line parameterization.
fn elclib_line_parameter(loc: DVec3, dir: DVec3, p: DVec3) -> f64 {
    (p - loc).dot(dir)
}

/// OCCT GeomAPI::To2d(C, P) (GeomAPI.cxx): projects the curve on the
/// plane (GeomProjLib::Curve2d). rcad: geom_proj_lib::curve_on_plane.
pub(crate) fn geom_api_to_2d(l: &Line3, pl: &Plane) -> Option<Curve2d> {
    curve_on_plane(&Curve3::Line(*l), [f64::NEG_INFINITY, f64::INFINITY], pl)
}

/// OCCT ChFi2d_Builder_0.cxx L702-797 — ComputePoint(V, E, D, Param).
/// Is internally used by <build_chamfer_edge>.
pub(crate) fn chfi2d_compute_point_v_e_d(
    brep: &mut BRep,
    v: &Shape,
    e: &Shape,
    d: f64,
    param: &mut f64,
) -> DVec3 {
    // geometric support
    // OCCT L705: BRepAdaptor_Curve c(E);
    let (c, c_range) = match brep.edge_curve_world(e) {
        Some(x) => x,
        None => return DVec3::ZERO,
    };
    let first = c_range[0];
    let last = c_range[1];

    let mut the_point;
    // OCCT L711: if (c.GetType() == GeomAbs_Line)
    if matches!(c, Curve3::Line(_)) {
        let p1;
        let p2;
        let v1 = brep.first_vertex(e);
        let v2 = brep.last_vertex(e);
        p1 = brep.vertex_position(&v1);
        p2 = brep.vertex_position(&v2);
        let mut my_vec = p2 - p1;
        my_vec = my_vec.normalize(); // OCCT L719: myVec.Normalize();
        my_vec *= d; // OCCT L720: myVec *= D;
        if v2.is_same(v) {
            my_vec *= -1.0; // change the sense of myVec
            the_point = p2 + my_vec; // OCCT L724: thePoint = p2.Translated(myVec);
        } else {
            the_point = p1 + my_vec; // OCCT L728: thePoint = p1.Translated(myVec);
        }

        // OCCT L731: Param = ElCLib::Parameter(c.Line(), thePoint);
        if let Curve3::Line(line) = &c {
            *param = elclib_line_parameter(line.origin, line.direction, the_point);
        }
        // szv:OCC20823-begin
        // OCCT L733: c.D0(Param, thePoint);
        the_point = c.point_at(*param);
        return the_point;
        // szv:OCC20823-end
    } // if (C->IsKind(TYPE ...

    // OCCT L738: if (c.GetType() == GeomAbs_Circle)
    if let Curve3::Circle(cir) = &c {
        let radius = cir.radius; // OCCT L741: double radius = cir.Radius();
        let v1 = brep.first_vertex(e);
        let v2 = brep.last_vertex(e);
        let param1;
        let param2;
        if v.is_same(&v1) {
            param1 = brep_tool_parameter(brep, &v1, e); // OCCT L747: BRep_Tool::Parameter(v1, E)
            param2 = brep_tool_parameter(brep, &v2, e); // OCCT L748
        } else {
            param1 = brep_tool_parameter(brep, &v2, e); // OCCT L751
            param2 = brep_tool_parameter(brep, &v1, e); // OCCT L752
        }
        let delta_alpha = d / radius; // OCCT L755: double deltaAlpha = D / radius;
        if param1 > param2 {
            *param = param1 - delta_alpha;
        } else {
            *param = param1 + delta_alpha;
        }
        // OCCT L764: c.D0(Param, thePoint);
        the_point = c.point_at(*param);
        return the_point;
    } // if (C->IsKind(TYPE ...
    // OCCT L768-796: else — in all other case than lines and circles.
    let v1 = brep.first_vertex(e);
    let v2 = brep.last_vertex(e);
    let p = if v.is_same(&v1) {
        brep.vertex_position(&v1)
    } else {
        brep.vertex_position(&v2)
    };

    // OCCT L783: const GeomAdaptor_Curve& cc = c.Curve();
    let cc = &c;
    // OCCT L784: if (p.Distance(c.Value(first)) <= Precision::Confusion())
    if p.distance(cc.point_at(first)) <= CONFUSION {
        // OCCT L786-787: GCPnts_AbscissaPoint computePoint(cc, D, first);
        //                Param = computePoint.Parameter();
        *param = abscissa_point_parameter(cc, first, last, d, first);
    } else {
        // OCCT L790-792: GCPnts_AbscissaPoint computePoint(cc, D, last);
        //                Param = computePoint.Parameter();
        *param = abscissa_point_parameter(cc, first, last, d, last);
    }
    // OCCT L794: thePoint = cc.Value(Param);
    the_point = cc.point_at(*param);
    the_point
} // ComputePoint

/// OCCT ChFi2d_Builder_0.cxx L801-838 — ComputePoint(F, L, E, Param).
/// Is internally used by <build_chamfer_edge>.
pub(crate) fn chfi2d_compute_point_face_line(
    brep: &mut BRep,
    f: &Shape,
    l: &Curve3,
    e: &Shape,
    param: &mut f64,
) -> DVec3 {
    // OCCT L806-807:
    // BRepAdaptor_Surface Adaptor3dSurface(F);
    // occ::handle<Geom_Plane> refSurf = new Geom_Plane(Adaptor3dSurface.Plane());
    let ref_plane = match brep.face_surface_world(f) {
        Some(Surface3::Plane(pl)) => pl,
        _ => panic!("Standard_TypeMismatch: F is not a plane"),
    };
    // OCCT L808: occ::handle<Geom2d_Curve> lin2d = GeomAPI::To2d(L, refSurf->Pln());
    let lin2d = match l {
        Curve3::Line(l3) => match geom_api_to_2d(l3, &ref_plane) {
            Some(c) => c,
            None => return DVec3::ZERO,
        },
        _ => panic!("Standard_TypeMismatch: L is not a line"),
    };
    // OCCT L809-811: c2d = BRep_Tool::CurveOnSurface(E, F, first, last);
    let (c2d, _first, _last) = match brep.curve_on_surface(e, f) {
        Some((c, fi, la)) => (c, fi, la),
        None => return DVec3::ZERO,
    };
    // OCCT L812-817:
    // Geom2dAdaptor_Curve adaptorL(lin2d);
    // Geom2dAdaptor_Curve adaptorC(c2d);
    // Geom2dInt_GInter    Intersection(adaptorL, adaptorC,
    //                                  Precision::PIntersection(),
    //                                  Precision::PIntersection());
    let intersection = GInter::new_cc(&lin2d, &c2d, p_intersection(), p_intersection());
    let mut param_on_line = 1.0e300f64; // OCCT L818: double paramOnLine = 1E300;
    let mut p2d = DVec2::ZERO; // OCCT L819: gp_Pnt2d p2d;
    if intersection.is_done() {
        let mut i = 1usize;
        while i <= intersection.nb_points() {
            // OCCT L825: IntRes2d_IntersectionPoint iP = Intersection.Point(i);
            let i_p = intersection.point(i);
            if i_p.param_on_first() < param_on_line {
                p2d = i_p.value();
                param_on_line = i_p.param_on_first();
                *param = i_p.param_on_second();
            } // if (iP.ParamOnFirst ...
            i += 1;
        } // while ( i <= ...
    } // if (Intersection.IsDone ...

    // OCCT L836: gp_Pnt thePoint = Adaptor3dSurface.Value(p2d.X(), p2d.Y());
    let the_point = ref_plane.point_at(p2d.x, p2d.y);
    the_point
} // ComputePoint

/// OCCT ChFi2d_Builder_0.cxx L842-865 — OrientChamfer(chamfer, E, V).
/// Is internally used by <build_chamfer_edge>.
pub(crate) fn orient_chamfer(brep: &BRep, chamfer: &mut Shape, e: &Shape, v: &Shape) {
    let v_orient;
    let orient = e.orientation; // OCCT L844: TopAbs_Orientation vOrient, orient = E.Orientation();
    // OCCT L845-854: TopExp::Vertices(E, v1, v2);
    let v1 = brep.first_vertex(e);
    let v2 = brep.last_vertex(e);
    if v1.is_same(v) {
        v_orient = v2.orientation;
    } else {
        v_orient = v1.orientation;
    }

    // OCCT L856-864:
    if (orient == Orientation::Forward && v_orient == Orientation::Forward)
        || (orient == Orientation::Reversed && v_orient == Orientation::Reversed)
    {
        chamfer.orientation = Orientation::Forward;
    } else {
        chamfer.orientation = Orientation::Reversed;
    }
} // OrientChamfer

/// OCCT ChFi2d_Builder_0.cxx L869-890 — IsLineOrCircle(E, F) (external
/// linkage version; ChFi2d_Builder.cxx carries its own file-static
/// copy).
pub(crate) fn is_line_or_circle(brep: &BRep, e: &Shape, f: &Shape) -> bool {
    // OCCT L876: occ::handle<Geom2d_Curve> C = BRep_Tool::CurveOnSurface(E, F, first, last);
    let c = match brep.curve_on_surface(e, f) {
        Some((c, _first, _last)) => c,
        None => return false,
    };
    // OCCT L877-885: down_cast<Geom2d_TrimmedCurve> -> BasisCurve.
    let basis_c = match &c {
        Curve2d::Trimmed(tc) => (*tc.curve).clone(),
        _ => c,
    };

    matches!(basis_c, Curve2d::Circle(_) | Curve2d::Line(_)) // else ...
} // IsLineOrCircle

/// OCCT BRep_Tool::Parameter(V, E) two-argument overload stand-in.
/// rcad stores the vertex parameter on the edge (vertex_params) and the
/// trait lookup ignores the face argument; the ChFi2d edges always lie
/// on newFace, so the stored parameter is the BRep_Tool::Parameter
/// value. Missing storage maps to the OCCT Standard_Failure of the
/// 2-arg overload.
pub(crate) fn brep_tool_parameter(brep: &BRep, v: &Shape, e: &Shape) -> f64 {
    match brep.parameter_on_edge(v, e, &Shape::null()) {
        Some(p) => p,
        None => panic!("Standard_Failure: BRep_Tool::Parameter(V, E) not found"),
    }
}

impl super::chfi2d_builder::ChFi2dBuilder {
    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L71-123 — AddChamfer(E1, E2, D1, D2)
    // =======================================================================

    /// OCCT ChFi2d_Builder::AddChamfer(const TopoDS_Edge& E1,
    /// const TopoDS_Edge& E2, const double D1, const double D2)
    /// (L71-123). Rust overload name: `add_chamfer_edges`.
    pub fn add_chamfer_edges(&mut self, e1: &Shape, e2: &Shape, d1: f64, d2: f64) -> Shape {
        let mut common_vertex_shape = Shape::null();
        let mut e1_mod = Shape::null();
        let mut e2_mod = Shape::null();
        let mut chamfer = Shape::null();

        // OCCT L80: bool hasConnection = ChFi2d::CommonVertex(E1, E2, commonVertex);
        let (cv, has_connection) = chfi2d_common_vertex(e1, e2);
        common_vertex_shape = cv;
        if !has_connection {
            return chamfer;
        }

        if self.is_a_fillet(e1)
            || self.is_a_chamfer(e1)
            || self.is_a_fillet(e2)
            || self.is_a_chamfer(e2)
        {
            self.status = ChFi2dConstructionError::NotAuthorized;
            return chamfer;
        } // if (IsAChamfer ...

        if !is_line_or_circle(&self.my_brep, e1, &self.new_face)
            || !is_line_or_circle(&self.my_brep, e2, &self.new_face)
        {
            self.status = ChFi2dConstructionError::NotAuthorized;
            return chamfer;
        } // if (!IsLineOrCircle ...

        // EE1 and EE2 are copies of E1 and E2 with the good orientation
        // on <new_face>
        let mut ee1 = Shape::null();
        let mut ee2 = Shape::null();
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &common_vertex_shape,
            &mut ee1,
            &mut ee2,
        );
        // OCCT L102-108: if (EE1.IsSame(E2)) { ... }
        if ee1.is_same(e2) {
            let orient = ee1.orientation;
            ee1 = ee2.clone();
            ee2 = e2.clone();
            ee2.orientation = orient;
        }

        self.compute_chamfer_edges(
            &common_vertex_shape,
            &ee1,
            &ee2,
            d1,
            d2,
            &mut e1_mod,
            &mut e2_mod,
            &mut chamfer,
        );
        if self.status == ChFi2dConstructionError::IsDone
            || self.status == ChFi2dConstructionError::FirstEdgeDegenerated
            || self.status == ChFi2dConstructionError::LastEdgeDegenerated
            || self.status == ChFi2dConstructionError::BothEdgesDegenerated
        {
            //  if (status == ChFi2d_IsDone) {
            self.build_new_wire(&ee1, &ee2, &e1_mod, &chamfer, &e2_mod);
            let basis_edge1 = self.basis_edge(&ee1);
            let basis_edge2 = self.basis_edge(&ee2);
            self.up_date_history_new_edge(&basis_edge1, &basis_edge2, &e1_mod, &e2_mod, &chamfer, 2);
            self.status = ChFi2dConstructionError::IsDone;
            return self.chamfers[self.chamfers.len() - 1].clone();
        }
        chamfer
    } // AddChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L127-176 — AddChamfer(E, V, D, Ang)
    // =======================================================================

    /// OCCT ChFi2d_Builder::AddChamfer(const TopoDS_Edge& E,
    /// const TopoDS_Vertex& V, const double D, const double Ang)
    /// (L127-176). Rust overload name: `add_chamfer_edge_vertex`.
    pub fn add_chamfer_edge_vertex(
        &mut self,
        e: &Shape,
        v: &Shape,
        d: f64,
        ang: f64,
    ) -> Shape {
        let mut a_chamfer = Shape::null();
        let mut adj_edge1 = Shape::null();
        let mut adj_edge2 = Shape::null();
        self.status =
            chfi2d_find_connected_edges(&self.new_face, v, &mut adj_edge1, &mut adj_edge2);
        if self.status == ChFi2dConstructionError::ConnexionError {
            return a_chamfer;
        }

        // adjEdge1 is a copy of E  with the good orientation
        // on <new_face>
        // OCCT L141-147: if (adjEdge2.IsSame(E)) { ... }
        if adj_edge2.is_same(e) {
            let orient = adj_edge2.orientation;
            adj_edge2 = adj_edge1.clone();
            adj_edge1 = e.clone();
            adj_edge1.orientation = orient;
        }

        if self.is_a_fillet(&adj_edge1)
            || self.is_a_chamfer(&adj_edge1)
            || self.is_a_fillet(&adj_edge2)
            || self.is_a_chamfer(&adj_edge2)
        {
            self.status = ChFi2dConstructionError::NotAuthorized;
            return a_chamfer;
        } // if (IsAChamfer ...

        if !is_line_or_circle(&self.my_brep, &adj_edge1, &self.new_face)
            || !is_line_or_circle(&self.my_brep, &adj_edge2, &self.new_face)
        {
            self.status = ChFi2dConstructionError::NotAuthorized;
            return a_chamfer;
        } // if (!IsLineOrCircle ...

        let mut e1 = Shape::null();
        let mut e2 = Shape::null();
        self.compute_chamfer_angle(v, &adj_edge1, d, ang, &adj_edge2, &mut e1, &mut e2, &mut a_chamfer);
        if self.status == ChFi2dConstructionError::IsDone
            || self.status == ChFi2dConstructionError::FirstEdgeDegenerated
            || self.status == ChFi2dConstructionError::LastEdgeDegenerated
            || self.status == ChFi2dConstructionError::BothEdgesDegenerated
        {
            //  if (status == ChFi2d_IsDone) {
            self.build_new_wire(&adj_edge1, &adj_edge2, &e1, &a_chamfer, &e2);
            let basis_edge1 = self.basis_edge(&adj_edge1);
            let basis_edge2 = self.basis_edge(&adj_edge2);
            self.up_date_history_new_edge(&basis_edge1, &basis_edge2, &e1, &e2, &a_chamfer, 2);
            self.status = ChFi2dConstructionError::IsDone;
            return self.chamfers[self.chamfers.len() - 1].clone();
        }
        a_chamfer
    }

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L180-212 — ComputeChamfer(V,E1,E2,D1,D2,..)
    // =======================================================================

    /// OCCT ChFi2d_Builder::ComputeChamfer(V, E1, E2, D1, D2, TrimE1,
    /// TrimE2, Chamfer) (L180-212). Rust overload name:
    /// `compute_chamfer_edges`.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_chamfer_edges(
        &mut self,
        v: &Shape,
        e1: &Shape,
        e2: &Shape,
        d1: f64,
        d2: f64,
        trim_e1: &mut Shape,
        trim_e2: &mut Shape,
        chamfer: &mut Shape,
    ) {
        let mut new_extr1 = Shape::null();
        let mut new_extr2 = Shape::null();
        let mut degen1 = false;
        let mut degen2 = false;
        *chamfer =
            self.build_chamfer_edge_edges(v, e1, e2, d1, d2, &mut new_extr1, &mut new_extr2);
        if self.status != ChFi2dConstructionError::IsDone {
            return;
        }
        *trim_e1 = self.build_new_edge_degenerated(e1, v, &new_extr1, &mut degen1);
        *trim_e2 = self.build_new_edge_degenerated(e2, v, &new_extr2, &mut degen2);
        if degen1 && degen2 {
            self.status = ChFi2dConstructionError::BothEdgesDegenerated;
        }
        if degen1 && !degen2 {
            self.status = ChFi2dConstructionError::FirstEdgeDegenerated;
        }
        if !degen1 && degen2 {
            self.status = ChFi2dConstructionError::LastEdgeDegenerated;
        }
        //   TrimE1 = BuildNewEdge(E1, V, newExtr1);
        //  TrimE2 = BuildNewEdge(E2, V, newExtr2);
    } // ComputeChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L216-248 — ComputeChamfer(V,E1,D,Ang,E2,..)
    // =======================================================================

    /// OCCT ChFi2d_Builder::ComputeChamfer(V, E1, D, Ang, E2, TrimE1,
    /// TrimE2, Chamfer) (L216-248). Rust overload name:
    /// `compute_chamfer_angle`.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_chamfer_angle(
        &mut self,
        v: &Shape,
        e1: &Shape,
        d: f64,
        ang: f64,
        e2: &Shape,
        trim_e1: &mut Shape,
        trim_e2: &mut Shape,
        chamfer: &mut Shape,
    ) {
        let mut new_extr1 = Shape::null();
        let mut new_extr2 = Shape::null();
        let mut degen1 = false;
        let mut degen2 = false;
        *chamfer =
            self.build_chamfer_edge_angle(v, e1, d, ang, e2, &mut new_extr1, &mut new_extr2);
        if self.status != ChFi2dConstructionError::IsDone {
            return;
        }
        *trim_e1 = self.build_new_edge_degenerated(e1, v, &new_extr1, &mut degen1);
        *trim_e2 = self.build_new_edge_degenerated(e2, v, &new_extr2, &mut degen2);
        if degen1 && degen2 {
            self.status = ChFi2dConstructionError::BothEdgesDegenerated;
        }
        if degen1 && !degen2 {
            self.status = ChFi2dConstructionError::FirstEdgeDegenerated;
        }
        if !degen1 && degen2 {
            self.status = ChFi2dConstructionError::LastEdgeDegenerated;
        }
        //   TrimE1 = BuildNewEdge(E1, V, newExtr1);
        //   TrimE2 = BuildNewEdge(E2, V, newExtr2);
    } // ComputeChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L252-279 — ModifyChamfer(Ch,E1,E2,D1,D2)
    // =======================================================================

    /// OCCT ChFi2d_Builder::ModifyChamfer(Chamfer, E1, E2, D1, D2)
    /// (L252-279). Rust overload name: `modify_chamfer_edges`.
    pub fn modify_chamfer_edges(
        &mut self,
        chamfer: &Shape,
        _e1: &Shape,
        e2: &Shape,
        d1: f64,
        d2: f64,
    ) -> Shape {
        let a_vertex = self.remove_chamfer(chamfer);
        let mut adj_edge1 = Shape::null();
        let mut adj_edge2 = Shape::null();
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &a_vertex,
            &mut adj_edge1,
            &mut adj_edge2,
        );
        let mut a_chamfer = Shape::null();
        if self.status == ChFi2dConstructionError::ConnexionError {
            return a_chamfer;
        }

        // adjEdge1 and adjEdge2 are copies of E1 and E2 with the good
        // orientation on <new_face>
        // OCCT L269-275: if (adjEdge1.IsSame(E2)) { ... }
        if adj_edge1.is_same(e2) {
            let orient = adj_edge1.orientation;
            adj_edge1 = adj_edge2.clone();
            adj_edge2 = e2.clone();
            adj_edge2.orientation = orient;
        }

        a_chamfer = self.add_chamfer_edges(&adj_edge1, &adj_edge2, d1, d2);
        a_chamfer
    } // ModifyChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L283-306 — ModifyChamfer(Ch,E,D,Ang)
    // =======================================================================

    /// OCCT ChFi2d_Builder::ModifyChamfer(Chamfer, E, D, Ang) (L283-306).
    /// Rust overload name: `modify_chamfer_angle`.
    pub fn modify_chamfer_angle(&mut self, chamfer: &Shape, e: &Shape, d: f64, ang: f64) -> Shape {
        let a_vertex = self.remove_chamfer(chamfer);
        let mut adj_edge1 = Shape::null();
        let mut adj_edge2 = Shape::null();
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &a_vertex,
            &mut adj_edge1,
            &mut adj_edge2,
        );
        let mut a_chamfer = Shape::null();
        if self.status == ChFi2dConstructionError::ConnexionError {
            return a_chamfer;
        }

        // OCCT L297-304:
        if adj_edge1.is_same(e) {
            a_chamfer = self.add_chamfer_edge_vertex(&adj_edge1, &a_vertex, d, ang);
        } else {
            a_chamfer = self.add_chamfer_edge_vertex(&adj_edge2, &a_vertex, d, ang);
        }
        a_chamfer
    } // ModifyChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L310-519 — RemoveChamfer
    // =======================================================================

    /// OCCT ChFi2d_Builder::RemoveChamfer(const TopoDS_Edge& Chamfer)
    /// (L310-519).
    pub fn remove_chamfer(&mut self, chamfer: &Shape) -> Shape {
        let mut common_vertex_shape = Shape::null();

        let mut i = 1usize; // OCCT: int i = 1;
        let mut is_find = false; // OCCT: int IsFind = false;
        while i <= self.chamfers.len() {
            let a_chamfer = self.chamfers[i - 1].clone(); // TopoDS::Edge(chamfers.Value(i))
            if a_chamfer.is_same(chamfer) {
                self.chamfers.remove(i - 1);
                is_find = true;
                break;
            }
            i += 1;
        }
        if !is_find {
            return common_vertex_shape;
        }

        let mut first_vertex;
        let mut last_vertex;
        // OCCT L333: TopExp::Vertices(Chamfer, firstVertex, lastVertex);
        first_vertex = self.my_brep.first_vertex(chamfer);
        last_vertex = self.my_brep.last_vertex(chamfer);

        let mut adj_edge1 = Shape::null();
        let mut adj_edge2 = Shape::null();
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &first_vertex,
            &mut adj_edge1,
            &mut adj_edge2,
        );
        if self.status == ChFi2dConstructionError::ConnexionError {
            return common_vertex_shape;
        }

        let basis_edge1;
        let basis_edge2;
        let e1;
        let mut e2;
        // E1 and E2 are the adjacentes edges to Chamfer

        if adj_edge1.is_same(chamfer) {
            e1 = adj_edge2.clone();
        } else {
            e1 = adj_edge1.clone();
        }
        basis_edge1 = self.basis_edge(&e1);
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &last_vertex,
            &mut adj_edge1,
            &mut adj_edge2,
        );
        if self.status == ChFi2dConstructionError::ConnexionError {
            return common_vertex_shape;
        }
        if adj_edge1.is_same(chamfer) {
            e2 = adj_edge2.clone();
        } else {
            e2 = adj_edge1.clone();
        }
        basis_edge2 = self.basis_edge(&e2);
        let mut connection_e1_chamfer = Shape::null();
        let mut connection_e2_chamfer = Shape::null();
        // OCCT: bool hasConnection = ChFi2d::CommonVertex(basisEdge1, basisEdge2, commonVertex);
        let (cv, mut has_connection) = chfi2d_common_vertex(&basis_edge1, &basis_edge2);
        common_vertex_shape = cv;
        if !has_connection {
            self.status = ChFi2dConstructionError::ConnexionError;
            return common_vertex_shape;
        }
        // OCCT: hasConnection = ChFi2d::CommonVertex(E1, Chamfer, connectionE1Chamfer);
        let (cv1, hc1) = chfi2d_common_vertex(&e1, chamfer);
        connection_e1_chamfer = cv1;
        has_connection = hc1;
        if !has_connection {
            self.status = ChFi2dConstructionError::ConnexionError;
            return common_vertex_shape;
        }
        // OCCT: hasConnection = ChFi2d::CommonVertex(E2, Chamfer, connectionE2Chamfer);
        let (cv2, hc2) = chfi2d_common_vertex(&e2, chamfer);
        connection_e2_chamfer = cv2;
        has_connection = hc2;
        if !has_connection {
            self.status = ChFi2dConstructionError::ConnexionError;
            return common_vertex_shape;
        }

        // rebuild edges on wire
        let mut new_edge1 = Shape::null();
        let mut new_edge2 = Shape::null();
        let mut v;
        let mut v1;
        let mut v2;

        // OCCT L395: TopExp::Vertices(E1, firstVertex, lastVertex);
        first_vertex = self.my_brep.first_vertex(&e1);
        last_vertex = self.my_brep.last_vertex(&e1);
        // OCCT L396: TopExp::Vertices(basisEdge1, v1, v2);
        v1 = self.my_brep.first_vertex(&basis_edge1);
        v2 = self.my_brep.last_vertex(&basis_edge1);
        if v1.is_same(&common_vertex_shape) {
            v = v2.clone();
        } else {
            v = v1.clone();
        }

        if first_vertex.is_same(&v) || last_vertex.is_same(&v) {
            // It means the edge support only one fillet. In this case
            // the new edge must be the basis edge.
            new_edge1 = basis_edge1.clone();
        } else {
            // It means the edge support one fillet on each end.
            if first_vertex.is_same(&connection_e1_chamfer) {
                // OCCT L420: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E1, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e1) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge1 = super::chfi2d_builder::brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &common_vertex_shape,
                    &last_vertex,
                );
                // OCCT L423-424:
                new_edge1.orientation = basis_edge1.orientation;
                new_edge1.location = basis_edge1.location;
            } // if (firstVertex ...
            else if last_vertex.is_same(&connection_e1_chamfer) {
                // OCCT L431: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E1, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e1) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge1 = super::chfi2d_builder::brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &first_vertex,
                    &common_vertex_shape,
                );
                // OCCT L434-435:
                new_edge1.orientation = basis_edge1.orientation;
                new_edge1.location = basis_edge1.location;
            } // else if (lastVertex ...
        } // else ...

        // OCCT L439: TopExp::Vertices(basisEdge2, v1, v2);
        v1 = self.my_brep.first_vertex(&basis_edge2);
        v2 = self.my_brep.last_vertex(&basis_edge2);
        if v1.is_same(&common_vertex_shape) {
            v = v2.clone();
        } else {
            v = v1.clone();
        }

        // OCCT L449: TopExp::Vertices(E2, firstVertex, lastVertex);
        first_vertex = self.my_brep.first_vertex(&e2);
        last_vertex = self.my_brep.last_vertex(&e2);
        if first_vertex.is_same(&v) || last_vertex.is_same(&v) {
            // It means the edge support only one fillet. In this case
            // the new edge must be the basis edge.
            new_edge2 = basis_edge2.clone();
        } else {
            // It means the edge support one fillet on each end.
            if first_vertex.is_same(&connection_e2_chamfer) {
                // OCCT L464: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E2, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e2) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge2 = super::chfi2d_builder::brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &common_vertex_shape,
                    &last_vertex,
                );
                // OCCT L467-468:
                new_edge2.orientation = basis_edge2.orientation;
                new_edge2.location = basis_edge2.location;
            } // if (firstVertex ...
            else if last_vertex.is_same(&connection_e2_chamfer) {
                // OCCT L475: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E2, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e2) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge2 = super::chfi2d_builder::brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &first_vertex,
                    &common_vertex_shape,
                );
                // OCCT L478-479:
                new_edge2.orientation = basis_edge2.orientation;
                new_edge2.location = basis_edge2.location;
            } // else if (lastVertex ...
        } // else ...

        // rebuild the newFace
        // OCCT L484: TopExp_Explorer Ex(newFace, TopAbs_EDGE);
        let ex = topexp_explore_face_edges(&self.my_brep, &self.new_face);
        // OCCT L487-488: BRep_Builder B; B.MakeWire(newWire);
        let mut b = BRepBuilder::new();
        let new_wire = b.make_wire(&mut self.my_brep);

        for the_edge in ex {
            // OCCT L492: const TopoDS_Edge& theEdge = TopoDS::Edge(Ex.Current());
            if !the_edge.is_same(&e1) && !the_edge.is_same(&e2) && !the_edge.is_same(chamfer) {
                b.add_to_wire(&mut self.my_brep, new_wire.clone(), the_edge);
            } else {
                // OCCT L499: if (theEdge == E1) — operator== is IsEqual.
                if the_edge.is_equal(&e1) {
                    b.add_to_wire(&mut self.my_brep, new_wire.clone(), new_edge1.clone());
                } else if the_edge.is_equal(&e2) {
                    b.add_to_wire(&mut self.my_brep, new_wire.clone(), new_edge2.clone());
                }
            } // else
        } // while ...

        // OCCT L511-514:
        // BRepAdaptor_Surface Adaptor3dSurface(refFace);
        // BRepLib_MakeFace mFace(Adaptor3dSurface.Plane(), newWire);
        // newFace.Nullify(); newFace = mFace;
        let ref_plane = match self.my_brep.face_surface_world(&self.ref_face) {
            Some(Surface3::Plane(pl)) => pl,
            _ => panic!("Standard_TypeMismatch: refFace is not a plane"),
        };
        let m_face = b.make_face(
            &mut self.my_brep,
            Some(Surface3::Plane(ref_plane)),
            new_wire.clone(),
        );
        self.new_face = m_face;

        self.up_date_history(&basis_edge1, &basis_edge2, &new_edge1, &new_edge2);

        common_vertex_shape
    } // RemoveChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L523-591 — BuildChamferEdge(V,E1,E2,D1,D2,..)
    // =======================================================================

    /// OCCT ChFi2d_Builder::BuildChamferEdge(V, AdjEdge1, AdjEdge2, D1,
    /// D2, NewExtr1, NewExtr2) (L523-591). The chamfer is computed from a
    /// vertex, two edges and two distances. Rust overload name:
    /// `build_chamfer_edge_edges`.
    #[allow(clippy::too_many_arguments)]
    pub fn build_chamfer_edge_edges(
        &mut self,
        v: &Shape,
        adj_edge1: &Shape,
        adj_edge2: &Shape,
        d1: f64,
        d2: f64,
        new_extr1: &mut Shape,
        new_extr2: &mut Shape,
    ) -> Shape {
        let mut chamfer = Shape::null(); // OCCT L531: TopoDS_Edge chamfer;
        if d1 <= 0.0 || d2 <= 0.0 {
            self.status = ChFi2dConstructionError::ParametersError;
            return chamfer;
        } // if ( D1 <=0 ...

        let mut param1 = 0.0f64;
        let mut param2 = 0.0f64;
        // OCCT L539-540:
        // gp_Pnt p1 = ComputePoint(V, AdjEdge1, D1, param1);
        // gp_Pnt p2 = ComputePoint(V, AdjEdge2, D2, param2);
        let p1 = chfi2d_compute_point_v_e_d(&mut self.my_brep, v, adj_edge1, d1, &mut param1);
        let p2 = chfi2d_compute_point_v_e_d(&mut self.my_brep, v, adj_edge2, d2, &mut param2);

        let tol = CONFUSION; // OCCT L542: double tol = Precision::Confusion();
        // OCCT L543-545: BRep_Builder B; B.MakeVertex(NewExtr1, p1, tol);
        //                B.MakeVertex(NewExtr2, p2, tol);
        *new_extr1 = brep_builder_make_vertex(&mut self.my_brep, p1, tol);
        *new_extr2 = brep_builder_make_vertex(&mut self.my_brep, p2, tol);
        // OCCT L546-547:
        (*new_extr1).orientation = Orientation::Forward;
        (*new_extr2).orientation = Orientation::Reversed;

        // chamfer edge construction
        // OCCT L551: const occ::handle<Geom_Surface> refSurf =
        //   BRep_Tool::Surface(refFace, loc);  (LOCAL surface, no location)
        let _ref_surf = self.my_brep.face_surface(&self.ref_face);
        // OCCT L552-554: gp_Vec myVec(p1, p2); gp_Dir myDir(myVec);
        //                occ::handle<Geom_Line> newLine = new Geom_Line(p1, myDir);
        let my_vec = p2 - p1;
        let my_dir = my_vec.normalize_or_zero();
        let new_line = Curve3::Line(Line3 {
            origin: p1,
            direction: my_dir,
        });
        // OCCT L555: double param = ElCLib::Parameter(newLine->Lin(), p2);
        let param = elclib_line_parameter(p1, my_dir, p2);
        // OCCT L556-561:
        // B.MakeEdge(chamfer, newLine, tol);
        // B.Range(chamfer, 0., param);
        // B.Add(chamfer, NewExtr1); B.UpdateVertex(NewExtr1, 0., chamfer, tol);
        // B.Add(chamfer, NewExtr2); B.UpdateVertex(NewExtr2, param, chamfer, tol);
        // rcad architecture: BRepBuilder::add_edge carries the curve, the
        // two vertices and the range in one step; set_vertex_param carries
        // B.UpdateVertex.
        {
            let mut b = BRepBuilder::new();
            chamfer = b.add_edge(
                &mut self.my_brep,
                Some(new_line.clone()),
                new_extr1.clone(),
                new_extr2.clone(),
                [0.0, param],
            );
            b.set_vertex_param(&mut self.my_brep, chamfer.clone(), new_extr1.clone(), 0.0);
            b.set_vertex_param(&mut self.my_brep, chamfer.clone(), new_extr2.clone(), param);
        }
        orient_chamfer(&self.my_brep, &mut chamfer, adj_edge1, v); // OCCT L562

        // set the orientation of NewExtr1 and NewExtr2 for the adjacent edges
        // OCCT L565-574:
        let v1 = self.my_brep.first_vertex(adj_edge1);
        let v2 = self.my_brep.last_vertex(adj_edge1);
        if v1.is_same(v) {
            (*new_extr1).orientation = v1.orientation;
        } else {
            (*new_extr1).orientation = v2.orientation;
        }

        // OCCT L576-585:
        let v1 = self.my_brep.first_vertex(adj_edge2);
        let v2 = self.my_brep.last_vertex(adj_edge2);
        if v1.is_same(v) {
            (*new_extr2).orientation = v1.orientation;
        } else {
            (*new_extr2).orientation = v2.orientation;
        }
        // OCCT L586-587: B.UpdateVertex(NewExtr1, param1, AdjEdge1, tol);
        //                B.UpdateVertex(NewExtr2, param2, AdjEdge2, tol);
        {
            let mut b = BRepBuilder::new();
            b.set_vertex_param(&mut self.my_brep, adj_edge1.clone(), new_extr1.clone(), param1);
            b.set_vertex_param(&mut self.my_brep, adj_edge2.clone(), new_extr2.clone(), param2);
        }

        self.status = ChFi2dConstructionError::IsDone;
        chamfer
    } // BuildChamferEdge

    // =======================================================================
    // OCCT ChFi2d_Builder_0.cxx L595-698 — BuildChamferEdge(V,E1,D,Ang,E2,..)
    // =======================================================================

    /// OCCT ChFi2d_Builder::BuildChamferEdge(V, AdjEdge1, D, Ang, AdjEdge2,
    /// NewExtr1, NewExtr2) (L595-698). The chamfer is computed from an
    /// edge, a vertex, a distance and an angle. Rust overload name:
    /// `build_chamfer_edge_angle`.
    #[allow(clippy::too_many_arguments)]
    pub fn build_chamfer_edge_angle(
        &mut self,
        v: &Shape,
        adj_edge1: &Shape,
        d: f64,
        ang: f64,
        adj_edge2: &Shape,
        new_extr1: &mut Shape,
        new_extr2: &mut Shape,
    ) -> Shape {
        let mut chamfer = Shape::null(); // OCCT L603: TopoDS_Edge chamfer;
        if d <= 0.0 || ang <= 0.0 {
            self.status = ChFi2dConstructionError::ParametersError;
            return chamfer;
        } // if ( D <= 0 ...

        let mut param1 = 0.0f64;
        let mut param2 = 0.0f64;
        // OCCT L611: gp_Pnt p1 = ComputePoint(V, AdjEdge1, D, param1);
        let p1 = chfi2d_compute_point_v_e_d(&mut self.my_brep, v, adj_edge1, d, &mut param1);
        // OCCT L612-613: gp_Pnt p = BRep_Tool::Pnt(V); gp_Vec myVec(p1, p);
        let p = self.my_brep.vertex_position(v);
        let my_vec = p - p1;

        // compute the tangent vector on AdjEdge2 at the vertex V.
        // OCCT L616: BRepAdaptor_Curve c(AdjEdge2, refFace);
        // rcad architecture: BRepAdaptor_Curve(E, F) is the
        // curve-on-surface adaptor; on the planar refFace the edge's
        // world 3D curve carries the same geometry (the constructor ran
        // BRepLib::BuildCurves3d).
        let (c, c_range) = match self.my_brep.edge_curve_world(adj_edge2) {
            Some(x) => x,
            None => return Shape::null(),
        };
        let first = c_range[0]; // OCCT L618: first = c.FirstParameter();
        let last = c_range[1]; // OCCT L619: last = c.LastParameter();

        // OCCT L622-623: gp_Pnt aPoint; gp_Vec tan; c.D1(first, aPoint, tan);
        let mut a_point = c.point_at(first);
        let mut tan = c.derivative_at(first);
        // OCCT L624: if (aPoint.Distance(p) > Precision::Confusion())
        //   c.D1(last, aPoint, tan);
        if a_point.distance(p) > CONFUSION {
            a_point = c.point_at(last);
            tan = c.derivative_at(last);
        }
        // tangent orientation
        // OCCT L629-639: TopExp::Vertices(AdjEdge2, v1, v2); orient = ...;
        let v1 = self.my_brep.first_vertex(adj_edge2);
        let v2 = self.my_brep.last_vertex(adj_edge2);
        let orient = if v1.is_same(v) {
            v1.orientation
        } else {
            v2.orientation
        };
        // OCCT L640-643: if (orient == TopAbs_REVERSED) tan *= -1;
        if orient == Orientation::Reversed {
            tan = -tan;
        }

        // compute the chamfer geometric support
        // OCCT L646-649:
        // gp_Ax1 RotAxe(p1, tan ^ myVec);          (cross product axis)
        // gp_Vec vecLin = myVec.Rotated(RotAxe, -Ang);
        // gp_Dir myDir(vecLin);
        // occ::handle<Geom_Line> newLine = new Geom_Line(p1, myDir);
        // rcad architecture: the gp_Trsf rotation is expressed by a
        // quaternion.
        let rot_axe = tan.cross(my_vec);
        let vec_lin = DQuat::from_axis_angle(rot_axe.normalize_or_zero(), -ang) * my_vec;
        let my_dir = vec_lin.normalize_or_zero();
        let new_line = Curve3::Line(Line3 {
            origin: p1,
            direction: my_dir,
        });
        // OCCT L650-651: BRep_Builder B1;
        //                B1.MakeEdge(chamfer, newLine, Precision::Confusion());
        // (the OCCT statement builds a placeholder edge that the final
        //  B.MakeEdge below overwrites; the Geom_Line support is what
        //  survives for ComputePoint)
        {
            let mut b1 = BRepBuilder::new();
            chamfer = b1.add_edge(
                &mut self.my_brep,
                Some(new_line.clone()),
                Shape::null(),
                Shape::null(),
                [0.0, 0.0],
            );
        }
        // OCCT L652: gp_Pnt p2 = ComputePoint(refFace, newLine, AdjEdge2, param2);
        let p2 = chfi2d_compute_point_face_line(
            &mut self.my_brep,
            &self.ref_face,
            &new_line,
            adj_edge2,
            &mut param2,
        );

        let tol = CONFUSION; // OCCT L654: double tol = Precision::Confusion();
        // OCCT L655-659:
        *new_extr1 = brep_builder_make_vertex(&mut self.my_brep, p1, tol);
        *new_extr2 = brep_builder_make_vertex(&mut self.my_brep, p2, tol);
        (*new_extr1).orientation = Orientation::Forward;
        (*new_extr2).orientation = Orientation::Reversed;

        // chamfer edge construction
        // OCCT L662: double param = ElCLib::Parameter(newLine->Lin(), p2);
        let param = elclib_line_parameter(p1, my_dir, p2);
        // OCCT L663-668: (same MakeEdge/Range/Add/UpdateVertex sequence as
        // the D1/D2 overload)
        {
            let mut b = BRepBuilder::new();
            chamfer = b.add_edge(
                &mut self.my_brep,
                Some(new_line.clone()),
                new_extr1.clone(),
                new_extr2.clone(),
                [0.0, param],
            );
            b.set_vertex_param(&mut self.my_brep, chamfer.clone(), new_extr1.clone(), 0.0);
            b.set_vertex_param(&mut self.my_brep, chamfer.clone(), new_extr2.clone(), param);
        }
        orient_chamfer(&self.my_brep, &mut chamfer, adj_edge1, v); // OCCT L669

        // set the orientation of NewExtr1 and NewExtr2 for the adjacent edges
        // OCCT L672-692:
        let v1 = self.my_brep.first_vertex(adj_edge1);
        let v2 = self.my_brep.last_vertex(adj_edge1);
        if v1.is_same(v) {
            (*new_extr1).orientation = v1.orientation;
        } else {
            (*new_extr1).orientation = v2.orientation;
        }

        let v1 = self.my_brep.first_vertex(adj_edge2);
        let v2 = self.my_brep.last_vertex(adj_edge2);
        if v1.is_same(v) {
            (*new_extr2).orientation = v1.orientation;
        } else {
            (*new_extr2).orientation = v2.orientation;
        }
        // OCCT L693-694: B.UpdateVertex(NewExtr1, param1, AdjEdge1, tol);
        //                B.UpdateVertex(NewExtr2, param2, AdjEdge2, tol);
        {
            let mut b = BRepBuilder::new();
            b.set_vertex_param(&mut self.my_brep, adj_edge1.clone(), new_extr1.clone(), param1);
            b.set_vertex_param(&mut self.my_brep, adj_edge2.clone(), new_extr2.clone(), param2);
        }

        self.status = ChFi2dConstructionError::IsDone;
        chamfer
    }
}
