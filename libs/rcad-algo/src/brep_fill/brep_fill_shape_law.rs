//! OCCT BRepFill_ShapeLaw (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_ShapeLaw.hxx (L36-81) + BRepFill_ShapeLaw.lxx (Edge) +
//! BRepFill_ShapeLaw.cxx (whole file L47-489).
//!
//! Architecture differences:
//! - The BRepFill_SectionLaw inheritance maps to the base struct +
//!   ops-trait form of brep_fill_section_law.rs (see its header).
//! - `handle(Law_Function) TheLaw` maps to the rcad
//!   `geomalgo::law::LawFunctionHandle` (Rc<RefCell<dyn LawFunction>>).
//! - `BRep_Tool::Curve(E, f, l)` maps to `TEdgeData::curve + range`; the
//!   owning `BRep` pool is passed as the leading `brep` argument.
//! - `gp_Trsf T; T.SetScale(gp_Pnt(0,0,0), f)` + BRepBuilderAPI_Transform
//!   (cxx L277-281 / L482-486) is a GAP carrier
//!   ([`brep_builder_api_transform`]).
//! - `GeomFill_EvolvedSection` (TKGeomAlgo) is a GAP carrier
//!   ([`evolved_section_new`]).
//! - `BRepLProp::Continuity` is a GAP carrier ([`brep_lprop_continuity`],
//!   the compatible_wires precedent).

use std::cell::RefCell;
use std::rc::Rc;

use glam::DVec3;

use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::geom::{CurveEval, BSplineCurve3, Curve3, Line3, TrimmedCurve3};
use rcad_kernel::topo::topods::GeomAbsShape;
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape};

use crate::brep_fill::brep_fill_section_law::{
    brep_tool_curve, brep_tool_is_closed_edge, expand_knots, reversed_curve,
    reversed_parameter_of, BRepFillSectionLawBase, BRepFillSectionLawOps, CompCurveToBSpline,
    WireExplorerState,
};
use crate::brep_fill::compatible_wires::{brep_tool_parameter, shape_tolerance, wire_edges};
use crate::brep_fill::generator::top_exp_vertices;
use crate::fillet::chfi3d_builder_0::topexp_common_vertex;
use crate::geomalgo::geomfill::section_law::SectionLaw;
use crate::geomalgo::geomfill::trihedron_law::{curve_first_parameter, curve_last_parameter};
use crate::geomalgo::geomfill::uniform_section::UniformSection;
use crate::geomalgo::law::law_function::LawFunctionHandle;

// ---------------------------------------------------------------------------
// GAP carriers
// ---------------------------------------------------------------------------

/// GAP carrier: OCCT GeomFill_EvolvedSection (TKGeomAlgo/GeomFill,
/// GeomFill_EvolvedSection.cxx) — the Law_Function-driven section law is not
/// translated; the construction sites (cxx L176 / L410) keep the OCCT
/// failure path.
pub fn evolved_section_new(_c: &Curve3, _law: &LawFunctionHandle) -> Rc<RefCell<dyn SectionLaw>> {
    panic!(
        "GAP: GeomFill_EvolvedSection (TKGeomAlgo/GeomFill) is not translated — \
         BRepFill_ShapeLaw::Init / ConcatenedLaw (BRepFill_ShapeLaw.cxx L176/L410)"
    )
}

/// GAP carrier: OCCT BRepBuilderAPI_Transform(S, T) (TKTopAlgo) — the shape
/// scale-transform engine is not translated; the Vertex / D0 call sites
/// (cxx L281 / L486) keep the OCCT failure path.
fn brep_builder_api_transform(_s: &Shape, _scale_factor: f64) -> Shape {
    panic!(
        "GAP: BRepBuilderAPI_Transform (TKTopAlgo) is not translated — \
         BRepFill_ShapeLaw::Vertex / D0 (BRepFill_ShapeLaw.cxx L281/L486)"
    )
}

/// GAP carrier: OCCT BRepLProp::Continuity(C1, C2, U1, U2, Eps, TolAngular)
/// (TKTopAlgo/BRepLProp) — the edge-continuity evaluator is not translated
/// (the compatible_wires precedent); the call site (cxx L470) keeps the OCCT
/// failure path.
fn brep_lprop_continuity(
    _brep: &BRep,
    _edge1: &Shape,
    _edge2: &Shape,
    _u1: f64,
    _u2: f64,
    _eps: f64,
    _tol_angular: f64,
) -> GeomAbsShape {
    panic!(
        "GAP: BRepLProp::Continuity (TKTopAlgo) is not translated — \
         BRepFill_ShapeLaw::Continuity (BRepFill_ShapeLaw.cxx L470)"
    )
}

/// The empty wire-explorer state (the default-constructed OCCT member).
fn empty_iterator() -> WireExplorerState {
    WireExplorerState {
        wire: Shape::null(),
        edges: Vec::new(),
        cursor: 0,
    }
}

// ---------------------------------------------------------------------------
// BRepFill_ShapeLaw
// ---------------------------------------------------------------------------

/// OCCT BRepFill_ShapeLaw (hxx L36-81): Build Section Law, with an Vertex,
/// or an Wire.
pub struct BRepFillShapeLaw {
    /// OCCT BRepFill_SectionLaw base sub-object.
    pub base: BRepFillSectionLawBase,
    /// OCCT bool vertex (protected).
    pub vertex: bool,
    /// OCCT TopoDS_Shape myShape (private).
    pub my_shape: Shape,
    /// OCCT handle(NCollection_HArray1<TopoDS_Shape>) myEdges (private).
    pub my_edges: Vec<Shape>,
    /// OCCT handle(Law_Function) TheLaw (private).
    pub the_law: Option<LawFunctionHandle>,
}

impl BRepFillShapeLaw {
    /// OCCT BRepFill_ShapeLaw(V, Build) (cxx L52-74) — Process the case of
    /// Vertex by constructing a line with the vertex in the origin.
    pub fn new_vertex(brep: &BRep, v: &Shape, build: bool) -> Self {
        let mut law = BRepFillShapeLaw {
            base: BRepFillSectionLawBase {
                my_laws: Vec::new(),
                uclosed: false,
                vclosed: true, // constant law
                my_done: false,
                my_indices: std::collections::HashMap::new(),
                my_iterator: empty_iterator(),
            },
            vertex: true,
            my_shape: v.clone(),
            my_edges: vec![v.clone()],
            the_law: None,
        };
        // TheLaw.Nullify();
        if build {
            // gp_Dir D(gp_Dir::D::X); // Following the normal
            // occ::handle<Geom_Line> L = new Geom_Line(BRep_Tool::Pnt(V), D);
            // double Last = 2 * BRep_Tool::Tolerance(V) + PConfusion();
            // occ::handle<Geom_TrimmedCurve> TC = new Geom_TrimmedCurve(L, 0, Last);
            // myLaws->ChangeValue(1) = new GeomFill_UniformSection(TC);
            let last = 2.0 * shape_tolerance(brep, v) + PCONFUSION;
            let line = Curve3::Line(Line3 {
                origin: brep.vertex(v.clone()).point,
                direction: DVec3::X,
            });
            let uniform = UniformSection::new(&line, 0.0, last);
            law.base.my_laws = vec![Rc::new(RefCell::new(uniform))];
        }
        law.base.my_done = true;
        law
    }

    /// OCCT BRepFill_ShapeLaw(W, Build) (cxx L78-86) — Construct an constant
    /// Law.
    pub fn new(brep: &BRep, w: &Shape, build: bool) -> Self {
        let mut law = BRepFillShapeLaw {
            base: BRepFillSectionLawBase {
                my_laws: Vec::new(),
                uclosed: false,
                vclosed: false,
                my_done: false,
                my_indices: std::collections::HashMap::new(),
                my_iterator: empty_iterator(),
            },
            vertex: false,
            my_shape: w.clone(),
            my_edges: Vec::new(),
            the_law: None,
        };
        // TheLaw.Nullify();
        law.init(brep, build);
        law.base.my_done = true;
        law
    }

    /// OCCT BRepFill_ShapeLaw(W, L, Build) (cxx L90-100) — Construct an
    /// evolutive Law.
    pub fn new_with_law(brep: &BRep, w: &Shape, l: LawFunctionHandle, build: bool) -> Self {
        let mut law = BRepFillShapeLaw {
            base: BRepFillSectionLawBase {
                my_laws: Vec::new(),
                uclosed: false,
                vclosed: false,
                my_done: false,
                my_indices: std::collections::HashMap::new(),
                my_iterator: empty_iterator(),
            },
            vertex: false,
            my_shape: w.clone(),
            my_edges: Vec::new(),
            the_law: Some(l),
        };
        law.init(brep, build);
        law.base.my_done = true;
        law
    }

    /// OCCT Init (cxx L106-228) — Case of the wire : Create a table of
    /// GeomFill_SectionLaw.
    pub fn init(&mut self, brep: &BRep, build: bool) {
        self.base.vclosed = true;
        let w = self.my_shape.clone();

        // L116-127: count the usable edges.
        let mut nb_edge = 0usize;
        for e in wire_edges(brep, &w) {
            if !e.is_null() && !brep.edge(e.clone()).degenerated {
                if brep_tool_curve(brep, &e).is_some() {
                    nb_edge += 1;
                }
            }
        }

        // myLaws = new HArray1(1, NbEdge); myEdges = new HArray1(1, NbEdge);
        self.base.my_laws = Vec::new();
        self.my_edges = vec![Shape::null(); nb_edge];
        let mut laws: Vec<Option<Rc<RefCell<dyn SectionLaw>>>> =
            (0..nb_edge).map(|_| None).collect();

        let mut ii = 0usize; // OCCT ii starts at 1; 0-based in the Vec

        // L134-182: fill the tables.
        for e in wire_edges(brep, &w) {
            if !e.is_null() && !brep.edge(e.clone()).degenerated {
                if let Some((mut c, mut first, mut last)) = brep_tool_curve(brep, &e) {
                    self.my_edges[ii] = e.clone();
                    ii += 1;
                    // myIndices.Bind(E, ii)
                    self.base.my_indices.insert(e.ptr_id(), ii as i32);
                    if build {
                        if e.orientation == Orientation::Reversed {
                            // CBis = C->Reversed(); // To avoid the
                            // deterioration of the topology
                            let cbis = reversed_curve(&c);
                            let aux = reversed_parameter_of(&c, first);
                            first = reversed_parameter_of(&c, last);
                            last = aux;
                            c = cbis;
                        }

                        let mut is_closed = brep_tool_is_closed_edge(brep, &e);
                        if is_closed && (curve_first_parameter(&c) - first).abs() > PCONFUSION {
                            is_closed = false; // trimmed curve differs
                        }

                        if (ii > 1) || !is_closed {
                            // Trim C
                            c = Curve3::Trimmed(TrimmedCurve3 {
                                curve: Box::new(c),
                                first,
                                last,
                            });
                        }
                        // otherwise preserve the integrity of the curve
                        let built: Rc<RefCell<dyn SectionLaw>> = match &self.the_law {
                            None => Rc::new(RefCell::new(UniformSection::new(
                                &c,
                                curve_first_parameter(&c),
                                curve_last_parameter(&c),
                            ))),
                            Some(the_law) => evolved_section_new(&c, the_law),
                        };
                        laws[ii - 1] = Some(built);
                    }
                }
            }
        }
        self.base.my_laws = laws.into_iter().flatten().collect();

        // L186-227: Is the law closed by U ?
        let mut uclosed = brep.has_flag(w.clone(), tshape_flags::CLOSED);
        if !uclosed {
            // if not sure about the flag, make check
            let edge1 = self.my_edges[self.my_edges.len() - 1].clone();
            let edge2 = self.my_edges[0].clone();

            let v1;
            let v2;
            if edge1.orientation == Orientation::Reversed {
                v1 = top_exp_vertices(brep, &edge1).0;
            } else {
                v1 = top_exp_vertices(brep, &edge1).1;
            }

            if edge2.orientation == Orientation::Reversed {
                v2 = top_exp_vertices(brep, &edge2).1;
            } else {
                v2 = top_exp_vertices(brep, &edge2).0;
            }
            if v1.is_same(&v2) {
                uclosed = true;
            } else {
                // BRepAdaptor_Curve Curve1(Edge1); Curve2(Edge2);
                // U1 = BRep_Tool::Parameter(V1, Edge1); U2 = ...;
                let u1 = brep_tool_parameter(brep, &v1, &edge1);
                let u2 = brep_tool_parameter(brep, &v2, &edge2);
                let eps = shape_tolerance(brep, &v2) + shape_tolerance(brep, &v1);

                let p1 = brep.edge(edge1.clone()).curve.clone().map(|c| c.point_at(u1));
                let p2 = brep.edge(edge2.clone()).curve.clone().map(|c| c.point_at(u2));
                uclosed = match (p1, p2) {
                    (Some(p1), Some(p2)) => p1.distance(p2) <= eps,
                    _ => false,
                };
            }
        }
        self.base.uclosed = uclosed;
    }

    /// OCCT Edge (BRepFill_ShapeLaw.lxx L25-30).
    pub fn edge(&self, index: i32) -> Shape {
        self.my_edges[(index - 1) as usize].clone()
    }

    /// OCCT IsVertex (cxx L232-235).
    pub fn is_vertex(&self) -> bool {
        self.vertex
    }

    /// OCCT IsConstant (cxx L239-242).
    pub fn is_constant_law(&self) -> bool {
        self.the_law.is_none()
    }

    /// OCCT Vertex (cxx L246-284).
    pub fn vertex(&self, brep: &mut BRep, index: i32, param: f64) -> Shape {
        let mut v = Shape::null();
        if index <= self.my_edges.len() as i32 {
            let e = self.my_edges[(index - 1) as usize].clone();
            if e.orientation == Orientation::Reversed {
                v = top_exp_vertices(brep, &e).1;
            } else {
                v = top_exp_vertices(brep, &e).0;
            }
        } else if index == self.my_edges.len() as i32 + 1 {
            let e = self.my_edges[(index - 2) as usize].clone();
            if e.orientation == Orientation::Reversed {
                v = top_exp_vertices(brep, &e).0;
            } else {
                v = top_exp_vertices(brep, &e).1;
            }
        }

        if let Some(the_law) = &self.the_law {
            // gp_Trsf T; T.SetScale(gp_Pnt(0, 0, 0), TheLaw->Value(Param));
            // V = TopoDS::Vertex(BRepBuilderAPI_Transform(V, T));
            let _scale = the_law.borrow_mut().value(param);
            v = brep_builder_api_transform(&v, _scale);
        }
        v
    }

    /// OCCT VertexTol (cxx L290-351) — Evaluate the hole between 2 edges of
    /// the section.
    pub fn vertex_tol(&self, index: i32, param: f64) -> f64 {
        let mut tol = CONFUSION;
        let i1;
        let i2;
        if (index == 0) || (index == self.my_edges.len() as i32) {
            if !self.base.uclosed {
                return tol; // The least possible error
            }
            i1 = self.my_edges.len() as i32;
            i2 = 1;
        } else {
            i1 = index;
            i2 = i1 + 1;
        }

        // First law: PFirst = BS->Value(Knots(Knots->Length()))
        let p_first = self.section_law_endpoint(i1, param, true);
        // Second law: Tol += PFirst.Distance(BS->Value(Knots->Value(1)))
        let p_last = self.section_law_endpoint(i2, param, false);
        tol += p_first.distance(p_last);
        tol
    }

    /// The shared Geom_BSplineCurve construction of VertexTol (cxx
    /// L317-332 / L334-349): SectionShape / D0(poles) / Knots / Mults ->
    /// Geom_BSplineCurve -> Value(first or last knot).
    fn section_law_endpoint(&self, index: i32, param: f64, first_knot: bool) -> DVec3 {
        let loi = self.base.my_laws[(index - 1) as usize].clone();
        let mut nb_poles = 0usize;
        let mut nb_knots = 0usize;
        let mut degree = 0usize;
        loi.borrow().section_shape(&mut nb_poles, &mut nb_knots, &mut degree);

        let mut poles = vec![DVec3::ZERO; nb_poles];
        let mut weight = vec![0.0f64; nb_poles];
        loi.borrow().d0(param, &mut poles, &mut weight);
        let mut knots = vec![0.0f64; nb_knots];
        loi.borrow().knots(&mut knots);
        let mut mults = vec![0i32; nb_knots];
        loi.borrow().mults(&mut mults);

        // BS = new Geom_BSplineCurve(Poles, Weigth, Knots, Mults, Degree,
        //                             Loi->IsUPeriodic());
        let bs = BSplineCurve3 {
            degree,
            knots: expand_knots(&knots, &mults),
            control_points: poles,
            weights: weight,
            is_periodic: loi.borrow().is_u_periodic(),
        };
        let bs_curve = Curve3::BSpline(bs);
        if first_knot {
            // BS->Value(Knots->Value(Knots->Length()))
            bs_curve.point_at(knots[knots.len() - 1])
        } else {
            // BS->Value(Knots->Value(1))
            bs_curve.point_at(knots[0])
        }
    }

    /// OCCT ConcatenedLaw (cxx L355-415).
    pub fn concatened_law(&self, brep: &BRep) -> Option<Rc<RefCell<dyn SectionLaw>>> {
        if self.base.my_laws.len() == 1 {
            return Some(self.base.my_laws[0].clone());
        }

        if !self.my_shape.is_null() {
            // Concatenation of edges
            let (composite, f, l) = brep_tool_curve(brep, &self.edge(1)).expect("edge curve");
            let tc = Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(composite),
                first: f,
                last: l,
            });
            let mut concat = CompCurveToBSpline::new(&tc);

            let mut bof = true;
            let mut eps_v = 0.0;
            let mut ii = 2i32;
            while ii <= self.my_edges.len() as i32 && bof {
                let (composite, f, l) =
                    brep_tool_curve(brep, &self.edge(ii)).expect("edge curve");
                let tc = Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(composite),
                    first: f,
                    last: l,
                });
                // Bof = TopExp::CommonVertex(Edge(ii-1), Edge(ii), V);
                match topexp_common_vertex(&self.edge(ii - 1), &self.edge(ii)) {
                    Some(v) => {
                        eps_v = shape_tolerance(brep, &v);
                        bof = true;
                    }
                    None => {
                        eps_v = 10.0 * PCONFUSION;
                        bof = false;
                    }
                }
                // Bof = Concat.Add(TC, epsV, true, false, 20);
                bof = concat.add(&tc, eps_v, true, false, 20);
                if !bof {
                    // Bof = Concat.Add(TC, 200 * epsV, true, false, 20);
                    bof = concat.add(&tc, 200.0 * eps_v, true, false, 20);
                }
                ii += 1;
            }
            let composite = concat.bspline_curve();

            return match &self.the_law {
                None => {
                    // Law = new GeomFill_UniformSection(Composite);
                    Some(Rc::new(RefCell::new(UniformSection::new(
                        &composite,
                        curve_first_parameter(&composite),
                        curve_last_parameter(&composite),
                    ))))
                }
                Some(the_law) => {
                    // Law = new GeomFill_EvolvedSection(Composite, TheLaw);
                    Some(evolved_section_new(&composite, the_law))
                }
            };
        }
        None
    }

    /// OCCT Continuity (cxx L419-473).
    pub fn continuity(&self, brep: &BRep, index: i32, tol_angular: f64) -> GeomAbsShape {
        let edge1;
        let edge2;
        if (index == 0) || (index == self.my_edges.len() as i32) {
            if !self.base.uclosed {
                return GeomAbsShape::C0; // The least possible error
            }

            edge1 = self.my_edges[self.my_edges.len() - 1].clone();
            edge2 = self.my_edges[0].clone();
        } else {
            edge1 = self.my_edges[(index - 1) as usize].clone();
            edge2 = self.my_edges[index as usize].clone();
        }

        // TopExp::Vertices(Edge1, vv1, vv2); TopExp::Vertices(Edge2, vv3, vv4);
        let (vv1, vv2) = top_exp_vertices(brep, &edge1);
        let (vv3, vv4) = top_exp_vertices(brep, &edge2);
        let (v1, v2) = if vv1.is_same(&vv3) {
            (vv1, vv3)
        } else if vv1.is_same(&vv4) {
            (vv1, vv4)
        } else if vv2.is_same(&vv3) {
            (vv2, vv3)
        } else {
            (vv2, vv4)
        };

        let u1 = brep_tool_parameter(brep, &v1, &edge1);
        let u2 = brep_tool_parameter(brep, &v2, &edge2);
        let eps = shape_tolerance(brep, &v2) + shape_tolerance(brep, &v1);
        brep_lprop_continuity(brep, &edge1, &edge2, u1, u2, eps, tol_angular)
    }

    /// OCCT D0 (cxx L477-488).
    pub fn d0(&mut self, u: f64, s: &mut Shape) {
        *s = self.my_shape.clone();
        if let Some(the_law) = &self.the_law {
            // gp_Trsf T; T.SetScale(gp_Pnt(0, 0, 0), TheLaw->Value(U));
            // S = BRepBuilderAPI_Transform(S, T);
            let scale = the_law.borrow_mut().value(u);
            *s = brep_builder_api_transform(s, scale);
        }
    }
}

/// OCCT BRepFill_SectionLaw inheritance: BRepFill_ShapeLaw overrides.
impl BRepFillSectionLawOps for BRepFillShapeLaw {
    fn base(&self) -> &BRepFillSectionLawBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut BRepFillSectionLawBase {
        &mut self.base
    }
    fn is_constant(&self) -> bool {
        self.is_constant_law()
    }
    fn is_vertex(&self) -> bool {
        self.is_vertex()
    }
    fn concatened_law(&self, brep: &BRep) -> Option<Rc<RefCell<dyn SectionLaw>>> {
        self.concatened_law(brep)
    }
    fn continuity(&self, brep: &BRep, index: i32, tol_angular: f64) -> GeomAbsShape {
        self.continuity(brep, index, tol_angular)
    }
    fn vertex_tol(&self, _brep: &BRep, index: i32, param: f64) -> f64 {
        self.vertex_tol(index, param)
    }
    fn vertex(&self, brep: &mut BRep, index: i32, param: f64) -> Shape {
        self.vertex(brep, index, param)
    }
    fn d0(&mut self, _brep: &mut BRep, u: f64, s: &mut Shape) {
        self.d0(u, s)
    }
}
