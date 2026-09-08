//! OCCT BRepSweep_Translation (TKPrim/BRepSweep) — builds a topology by
//! translation sweep.
//!
//! Sources:
//! - BRepSweep_Translation.hxx L34-184
//! - BRepSweep_Translation.cxx L52-532
//!
//! Architecture difference (Rust has no inheritance): the C++
//! `Translation : BRepSweep_Trsf : NumLinearRegularSweep` chain maps to the
//! composition root owning the [`super::num_linear_regular_sweep::
//! NumLinearRegularSweepCore`] + [`super::trsf::TrsfData`] members and
//! implementing the two slot traits.  The inheritance chains are recorded in
//! the struct docs.

use rcad_kernel::geom::{Curve2d, LinearExtrusionSurface, Surface3};
use rcad_kernel::precision::{p_confusion, CONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};

use super::num_linear_regular_sweep::{NumLinearRegularSweepCore, NumLinearRegularSweepSlots};
use super::sweep_num_shape::SweepNumShape;
use super::tool_rehost::{
    brep_tool_curve_on_surface_seam, brep_tool_is_closed_edge_face, brep_tools_is_really_closed,
    geom_adaptor_surface_get_type, geom_curve_transformed, geom_line,
    GeomAdaptorSurfaceOfLinearExtrusion, GeomAbsSurfaceType,
};
use super::trsf::{self, BRepSweepTrsfSlots, TrsfData};
use super::BRepSweepBuilder;

use crate::brep_algo::tool::{brep_tool_curve, brep_tool_parameter, brep_tool_pnt, brep_tool_surface, brep_tool_tolerance};

/// OCCT BRepSweep_Translation (BRepSweep_Translation.hxx L34-184).
///
/// OCCT inheritance chain: `BRepSweep_Translation : BRepSweep_Trsf :
/// BRepSweep_NumLinearRegularSweep`.
pub struct BRepSweepTranslation {
    /// The base-class members.
    pub core: NumLinearRegularSweepCore,
    /// The BRepSweep_Trsf members (myLocation / myCopy).
    pub trsf: TrsfData,
    /// OCCT: myVec (gp_Vec).
    pub my_vec: glam::DVec3,
    /// OCCT: myCanonize.
    pub my_canonize: bool,
}

impl NumLinearRegularSweepSlots for BRepSweepTranslation {
    fn core(&mut self) -> &mut NumLinearRegularSweepCore {
        &mut self.core
    }
    fn core_ref(&self) -> &NumLinearRegularSweepCore {
        &self.core
    }

    /// OCCT BRepSweep_Translation::MakeEmptyVertex (cxx L104-120) — only
    /// called when the option of construction is with copy.
    fn make_empty_vertex(&mut self, a_gen_v: &Shape, a_dir_v: &SweepNumShape) -> Shape {
        if !self.trsf.my_copy {
            panic!("Standard_ConstructionError: BRepSweep_Translation::MakeEmptyVertex");
        }
        // OCCT L109: P = BRep_Tool::Pnt(Vertex(aGenV)).
        let mut p = brep_tool_pnt(a_gen_v).expect("BRep_Tool::Pnt");
        if a_dir_v.index() == 2 {
            // OCCT L110-113: P.Transform(myLocation.Transformation()).
            p = self.trsf.my_location.transform_point3(p);
        }
        // OCCT L114-118 (modified by jgv, 5.10.01, for buc61008): the vertex
        // tolerance is the generatrix vertex tolerance.
        let tol = brep_tool_tolerance(a_gen_v);
        self.core.my_builder.my_builder.make_vertex(p, tol)
    }

    /// OCCT BRepSweep_Translation::MakeEmptyDirectingEdge (cxx L124-133).
    fn make_empty_directing_edge(&mut self, a_gen_v: &Shape, _a_dir_e: &SweepNumShape) -> Shape {
        let p = brep_tool_pnt(a_gen_v).expect("BRep_Tool::Pnt");
        // OCCT L128-129: gp_Lin L(P, myVec); Geom_Line GL = new Geom_Line(L).
        let gl = geom_line(p, self.my_vec);
        // OCCT L131: MakeEdge(E, GL, BRep_Tool::Tolerance(Vertex(aGenV))).
        let tol = brep_tool_tolerance(a_gen_v);
        self.core.my_builder.my_builder.make_edge_curve(&gl, tol)
    }

    /// OCCT BRepSweep_Translation::MakeEmptyGeneratingEdge (cxx L137-166) —
    /// call only in case of construction with copy.
    fn make_empty_generating_edge(&mut self, a_gen_e: &Shape, a_dir_v: &SweepNumShape) -> Shape {
        if !self.trsf.my_copy {
            panic!("Standard_ConstructionError: BRepSweep_Translation::MakeEmptyVertex");
        }
        let tol = brep_tool_tolerance(a_gen_e);
        if self.core_my_degenerated(a_gen_e) {
            // OCCT L144-148: the degenerated edge — MakeEdge(E);
            // UpdateEdge(E, Tolerance); Degenerated(E, true).
            let new_e = self.core.my_builder.my_builder.make_edge();
            self.core.my_builder.my_builder.update_edge_tol(&new_e, tol);
            self.core.my_builder.my_builder.degenerated(&new_e, true);
            new_e
        } else {
            // OCCT L151-153: C = BRep_Tool::Curve(Edge(aGenE), L, First,
            // Last).  Architecture difference: the rcad 3D curve carries no
            // own location (the OCCT C->Transform(L.Transformation()) with
            // the curve location is the identity no-op here — arch. diff.).
            let c_opt = brep_tool_curve(a_gen_e).map(|(c, _f, _l)| c);
            let mut c = c_opt;
            if let Some(hit) = &mut c {
                // OCCT L156-157: C = C->Copy(); C->Transform(L.Transformation()).
                // OCCT L158-161: if (aDirV.Index() == 2)
                //   C->Transform(myLocation.Transformation()).
                if a_dir_v.index() == 2 {
                    *hit = geom_curve_transformed(hit, &self.trsf.my_location);
                }
            }
            // OCCT L163: MakeEdge(newE, C, BRep_Tool::Tolerance(Edge(aGenE))).
            match &c {
                Some(hit) => self
                    .core
                    .my_builder
                    .my_builder
                    .make_edge_curve(hit, tol),
                None => self.core.my_builder.my_builder.make_edge(),
            }
        }
    }

    /// OCCT BRepSweep_Translation::SetParameters (cxx L170-183) — glues the
    /// parameter of vertices directly included in cap faces.
    fn set_parameters(
        &mut self,
        a_new_face: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_f: &Shape,
        a_gen_v: &Shape,
        _a_dir_v: &SweepNumShape,
    ) {
        // OCCT L177: pnt2d = BRep_Tool::Parameters(Vertex(aGenV), Face(aGenF)).
        let pnt2d = super::tool_rehost::brep_tool_parameters(a_gen_v, a_gen_f);
        // OCCT L178-182: UpdateVertex(Vertex(aNewVertex), pnt2d.X(),
        // pnt2d.Y(), Face(aNewFace), Precision::PConfusion()).
        self.core
            .my_builder
            .my_builder
            .update_vertex_uv_on_face(a_new_vertex, pnt2d.x, pnt2d.y, a_new_face, p_confusion());
    }

    /// OCCT BRepSweep_Translation::SetDirectingParameter (cxx L187-202).
    fn set_directing_parameter(
        &mut self,
        a_new_edge: &Shape,
        a_new_vertex: &mut Shape,
        _a_gen_v: &Shape,
        _a_dir_e: &SweepNumShape,
        a_dir_v: &SweepNumShape,
    ) {
        let mut param = 0.0;
        if a_dir_v.index() == 2 {
            param = self.my_vec.length();
        }
        self.core
            .my_builder
            .my_builder
            .update_vertex_param_on_edge(a_new_vertex, param, a_new_edge, p_confusion());
    }

    /// OCCT BRepSweep_Translation::SetGeneratingParameter (cxx L206-218).
    fn set_generating_parameter(
        &mut self,
        a_new_edge: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_e: &Shape,
        a_gen_v: &Shape,
        _a_dir_v: &SweepNumShape,
    ) {
        // OCCT L212-213: vbid = Vertex(aNewVertex);
        // vbid.Orientation(aGenV.Orientation()).
        a_new_vertex.orientation = a_gen_v.orientation;
        // OCCT L214-217: UpdateVertex(vbid, Parameter(aGenV, aGenE), aNewEdge,
        // PConfusion()).
        let param = brep_tool_parameter(a_gen_v, a_gen_e);
        self.core
            .my_builder
            .my_builder
            .update_vertex_param_on_edge(a_new_vertex, param, a_new_edge, p_confusion());
    }

    /// OCCT BRepSweep_Translation::MakeEmptyFace (cxx L222-279).
    fn make_empty_face(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> Shape {
        let toler;
        let s: Surface3;
        if self.core.my_dir_shape_tool.type_of(a_dir_s) == ShapeType::Edge {
            // OCCT L230-232: C = BRep_Tool::Curve(Edge(aGenS), L, First,
            // Last); toler = Tolerance(Edge(aGenS)).
            let (c0, _first, _last) =
                brep_tool_curve(a_gen_s).expect("BRep_Tool::Curve");
            toler = brep_tool_tolerance(a_gen_s);
            // OCCT L234-237: Tr = L.Transformation(); C = C->Copy();
            // C->Transform(Tr) (extruded surfaces are inverted
            // correspondingly to the topology, so reverse).
            let mut c = c0;
            // OCCT L238-239: gp_Dir D(myVec); D.Reverse().
            let d = -self.my_vec;
            if self.my_canonize {
                // OCCT L243-244: HC = GeomAdaptor_Curve(C, First, Last);
                // AS = GeomAdaptor_SurfaceOfLinearExtrusion(HC, D).
                let as_extr = GeomAdaptorSurfaceOfLinearExtrusion::load(&c, d);
                match as_extr.get_type() {
                    GeomAbsSurfaceType::Plane => {
                        // OCCT L248-249: S = new Geom_Plane(AS.Plane()).
                        s = Surface3::Plane(as_extr.plane());
                    }
                    GeomAbsSurfaceType::Cylinder => {
                        // OCCT L251-252: S = new Geom_CylindricalSurface(
                        // AS.Cylinder()).
                        s = Surface3::Cylinder(as_extr.cylinder());
                    }
                    _ => {
                        // OCCT L254-255: S = new
                        // Geom_SurfaceOfLinearExtrusion(C, D).
                        s = Surface3::LinearExtrusion(LinearExtrusionSurface {
                            profile: Box::new(c.clone()),
                            direction: d.normalize_or_zero(),
                        });
                    }
                }
            } else {
                // OCCT L260-261: S = new Geom_SurfaceOfLinearExtrusion(C, D).
                s = Surface3::LinearExtrusion(LinearExtrusionSurface {
                    profile: Box::new(c.clone()),
                    direction: d.normalize_or_zero(),
                });
            }
            let _ = &mut c;
        } else {
            // OCCT L266-268: S = BRep_Tool::Surface(Face(aGenS), L);
            // toler = Tolerance(Face(aGenS)).
            let mut s0 = brep_tool_surface(a_gen_s).expect("BRep_Tool::Surface");
            toler = brep_tool_tolerance(a_gen_s);
            // OCCT L269-271: Tr = L.Transformation(); S = S->Copy();
            // S->Transform(Tr) — the identity location no-op (arch. diff.).
            // OCCT L272-275: if (aDirS.Index() == 2) S->Translate(myVec).
            if a_dir_s.index() == 2 {
                s0 = rcad_kernel::geom::transform_surface(
                    &s0,
                    &glam::DAffine3::from_translation(self.my_vec),
                );
            }
            s = s0;
        }
        // OCCT L277: MakeFace(F, S, toler).
        self.core.my_builder.my_builder.make_face(&s, toler)
    }

    /// OCCT BRepSweep_Translation::SetPCurve (cxx L283-317) — sets on the
    /// edges of cap faces the same pcurves as the edges of the generating
    /// face.
    fn set_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_f: &Shape,
        a_gen_e: &Shape,
        _a_dir_v: &SweepNumShape,
        _orien: Orientation,
    ) {
        // OCCT L292: isclosed = BRep_Tool::IsClosed(Edge(aGenE), Face(aGenF)).
        let isclosed = brep_tool_is_closed_edge_face(a_gen_e, a_gen_f);
        if isclosed {
            // OCCT L295-301: anE = Edge(aGenE.Oriented(FORWARD));
            // aC1 = CurveOnSurface(anE, Face(aGenF)); anE.Reverse();
            // aC2 = CurveOnSurface(anE, Face(aGenF)).
            let mut an_e = a_gen_e.clone();
            an_e.orientation = Orientation::Forward;
            let a_c1 = brep_tool_curve_on_surface_seam(&an_e, a_gen_f);
            an_e.orientation = Orientation::Reversed;
            let a_c2 = brep_tool_curve_on_surface_seam(&an_e, a_gen_f);
            // OCCT L302-306: UpdateEdge(Edge(aNewEdge), aC1, aC2,
            // Face(aNewFace), PConfusion()).
            if let (Some(c1), Some(c2)) = (&a_c1, &a_c2) {
                self.core
                    .my_builder
                    .my_builder
                    .update_edge_two_pcurves(a_new_edge, &c1.0, &c2.0, a_new_face, p_confusion());
            }
        } else {
            // OCCT L310-316: UpdateEdge(Edge(aNewEdge), CurveOnSurface(
            // Edge(aGenE), Face(aGenF), First, Last), Face(aNewFace),
            // PConfusion()).
            if let Some((c, _f0, _l0)) = brep_tool_curve_on_surface_seam(a_gen_e, a_gen_f) {
                self.core
                    .my_builder
                    .my_builder
                    .update_edge_pcurve(a_new_edge, &c, a_new_face, p_confusion());
            }
        }
    }

    /// OCCT BRepSweep_Translation::SetGeneratingPCurve (cxx L321-378).
    fn set_generating_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        _a_gen_e: &Shape,
        _a_dir_e: &SweepNumShape,
        a_dir_v: &SweepNumShape,
        orien: Orientation,
    ) {
        // OCCT L329: AS = GeomAdaptor_Surface(BRep_Tool::Surface(
        // Face(aNewFace), Loc)).
        let surf = brep_tool_surface(a_new_face).expect("BRep_Tool::Surface");
        if geom_adaptor_surface_get_type(&surf) == GeomAbsSurfaceType::Plane {
            // OCCT L335-364: nothing is done JAG (the commented body is
            // skipped — the OCCT literal empty branch).
        } else {
            // OCCT L367-373: v = 0; if (aDirV.Index() == 2) v =
            // -myVec.Magnitude(); L.SetLocation(gp_Pnt2d(0, v));
            // L.SetDirection(gp_Dir2d(gp_Dir2d::D::X)).
            let mut v = 0.0;
            if a_dir_v.index() == 2 {
                v = -self.my_vec.length();
            }
            let gl = super::tool_rehost::geom2d_line(glam::DVec2::new(0.0, v), glam::DVec2::X);
            // OCCT L375-376: GL = new Geom2d_Line(L); SetThePCurve(B,
            // Edge(aNewEdge), Face(aNewFace), orien, GL).
            set_the_pcurve(
                self,
                a_new_edge,
                a_new_face,
                orien,
                &gl,
                CONFUSION,
            );
        }
    }

    /// OCCT BRepSweep_Translation::SetDirectingPCurve (cxx L382-416).
    fn set_directing_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_e: &Shape,
        a_gen_v: &Shape,
        _a_dir_e: &SweepNumShape,
        orien: Orientation,
    ) {
        let surf = brep_tool_surface(a_new_face).expect("BRep_Tool::Surface");
        if geom_adaptor_surface_get_type(&surf) != GeomAbsSurfaceType::Plane {
            // OCCT L394-395: L.SetLocation(gp_Pnt2d(BRep_Tool::Parameter(
            // Vertex(aGenV), Edge(aGenE)), 0));
            // L.SetDirection(gp_Dir2d(gp_Dir2d::D::NY)).
            let par = brep_tool_parameter(a_gen_v, a_gen_e);
            let gl = super::tool_rehost::geom2d_line(glam::DVec2::new(par, 0.0), glam::DVec2::Y);
            // OCCT L413-414: GL = new Geom2d_Line(L); SetThePCurve(...).
            set_the_pcurve(self, a_new_edge, a_new_face, orien, &gl, CONFUSION);
        }
    }

    /// OCCT BRepSweep_Translation::DirectSolid (cxx L420-436) — compares the
    /// face normal and the direction.
    fn direct_solid(&mut self, a_gen_s: &Shape, _a_dir_s: &SweepNumShape) -> Orientation {
        // OCCT L424-431: BRepAdaptor_Surface surf(Face(aGenS)); surf.D1(
        // midU, midV, P, du, dv).
        let surf = brep_tool_surface(a_gen_s).expect("BRep_Tool::Surface");
        use rcad_kernel::geom::SurfaceEval;
        let dom = surf.default_domain();
        let mid_u = (dom[0] + dom[1]) / 2.0;
        let mid_v = (dom[2] + dom[3]) / 2.0;
        let (_p, du, dv) = surf.derivatives(mid_u, mid_v);
        // OCCT L433: x = myVec.DotCross(du, dv).
        let x = self.my_vec.dot(du.cross(dv));
        // OCCT L434: orient = (x > 0) ? TopAbs_REVERSED : TopAbs_FORWARD.
        if x > 0.0 {
            Orientation::Reversed
        } else {
            Orientation::Forward
        }
    }

    /// OCCT BRepSweep_Translation::GGDShapeIsToAdd (cxx L440-447).
    fn ggd_shape_is_to_add(
        &self,
        _a_new_shape: &Shape,
        _a_new_sub_shape: &Shape,
        _a_gen_s: &Shape,
        _a_sub_gen_s: &Shape,
        _a_dir_s: &SweepNumShape,
    ) -> bool {
        true
    }

    /// OCCT BRepSweep_Translation::GDDShapeIsToAdd (cxx L451-458).
    fn gdd_shape_is_to_add(
        &self,
        _a_new_shape: &Shape,
        _a_new_sub_shape: &Shape,
        _a_gen_s: &Shape,
        _a_dir_s: &SweepNumShape,
        _a_sub_dir_s: &SweepNumShape,
    ) -> bool {
        true
    }

    /// OCCT BRepSweep_Translation::SeparatedWires (cxx L462-469) — here it
    /// always returns false.
    fn separated_wires(
        &self,
        _a_new_shape: &Shape,
        _a_new_sub_shape: &Shape,
        _a_gen_s: &Shape,
        _a_sub_gen_s: &Shape,
        _a_dir_s: &SweepNumShape,
    ) -> bool {
        false
    }

    /// OCCT BRepSweep_Translation::HasShape (cxx L473-518).
    fn has_shape(&self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> bool {
        if self.core.my_dir_shape_tool.type_of(a_dir_s) == ShapeType::Edge {
            if self.core.my_gen_shape_tool.type_of(a_gen_s) == ShapeType::Edge {
                // OCCT L480-485: E = Edge(aGenS); check if the edge is
                // degenerated.
                if self.core_my_degenerated(a_gen_s) {
                    return false;
                }
                // OCCT L504-512: check if the edge is a sewing edge
                // (TopExp_Explorer FaceExp(myGenShape, FACE)).
                for fac in crate::brep_algo::tool::explorer(
                    &self.core.my_gen_shape,
                    ShapeType::Face,
                    ShapeType::Shape,
                ) {
                    if brep_tools_is_really_closed(a_gen_s, &fac) {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// OCCT BRepSweep_Translation::IsInvariant (cxx L522-525) — always false
    /// because here the transformation is a translation.
    fn is_invariant(&self, _a_gen_s: &Shape) -> bool {
        false
    }

    /// OCCT BRepSweep_Trsf::SetContinuity (cxx L97-201) — Translation does
    /// not override; the Trsf implementation serves.
    fn set_continuity(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) {
        trsf::set_continuity(self, a_gen_s, a_dir_s);
    }
}

impl BRepSweepTrsfSlots for BRepSweepTranslation {
    fn trsf(&mut self) -> &mut TrsfData {
        &mut self.trsf
    }
    fn trsf_ref(&self) -> &TrsfData {
        &self.trsf
    }
}

impl BRepSweepTranslation {
    /// OCCT BRepSweep_Translation::BRepSweep_Translation(S, N, L, V, C,
    /// Canonize) (cxx L86-100) — creates a topology by translating <S> with
    /// the vector <V>.  If C is true S subcomponents are copied; if Canonize
    /// is true then generated surfaces are attempted to be canonized in
    /// simple types.
    pub fn new(
        s: &Shape,
        n: &SweepNumShape,
        l: glam::DAffine3,
        v: glam::DVec3,
        c: bool,
        canonize: bool,
    ) -> Self {
        // OCCT L92: BRepSweep_Trsf(BRep_Builder(), S, N, L, C).
        let core = NumLinearRegularSweepCore::new(BRepSweepBuilder::new(), s, n);
        let trsf = TrsfData {
            my_move_memo: std::cell::RefCell::new(indexmap::IndexMap::new()),
            my_location: l,
            my_copy: c,
        };
        let mut me = BRepSweepTranslation {
            core,
            trsf,
            my_vec: v,
            my_canonize: canonize,
        };
        // OCCT L97-98: Standard_ConstructionError_Raise_if(V.Magnitude() <
        // Precision::Confusion(), "BRepSweep_Translation::Constructor").
        if v.length() < CONFUSION {
            panic!("Standard_ConstructionError: BRepSweep_Translation::Constructor");
        }
        // OCCT L99: Init().
        trsf::init(&mut me);
        me
    }

    /// OCCT BRepSweep_Translation::Vec() (cxx L529-532) — the vector of the
    /// Prism; if it is an infinite prism the Vec is unitar.
    pub fn vec(&self) -> glam::DVec3 {
        self.my_vec
    }

    /// OCCT BRep_Tool::Degenerated(E) (the shared read).
    fn core_my_degenerated(&self, e: &Shape) -> bool {
        super::tool_rehost::brep_tool_degenerated(e)
    }
}

/// OCCT SetThePCurve (BRepSweep_Translation.cxx L52-82) — checks if there is
/// already a pcurve on non planar faces; the tolerance argument is
/// Precision::Confusion for Translation (Rotation uses ComputeTolerance).
fn set_the_pcurve(
    me: &mut BRepSweepTranslation,
    e: &Shape,
    f: &Shape,
    o: Orientation,
    c: &Curve2d,
    tol: f64,
) {
    // OCCT L62-65: GP = down_cast<Geom_Plane>(Surface(F, SL)); if the face is
    // not planar, OC = CurveOnSurface(E, F, f, l).
    let oc = match brep_tool_surface(f) {
        Some(Surface3::Plane(_)) => None,
        _ => brep_tool_curve_on_surface_seam(e, f),
    };
    match oc {
        None => {
            // OCCT L67-70: B.UpdateEdge(E, C, F, Precision::Confusion()).
            me.core
                .my_builder
                .my_builder
                .update_edge_pcurve(e, c, f, tol);
        }
        Some((oc_c, _f0, _l0)) => {
            // OCCT L71-81: if (O == REVERSED) UpdateEdge(E, OC, C, F, tol);
            // else UpdateEdge(E, C, OC, F, tol).
            if o == Orientation::Reversed {
                me.core
                    .my_builder
                    .my_builder
                    .update_edge_two_pcurves(e, &oc_c, c, f, tol);
            } else {
                me.core
                    .my_builder
                    .my_builder
                    .update_edge_two_pcurves(e, c, &oc_c, f, tol);
            }
        }
    }
}
