//! OCCT BRepFill_SectionLaw (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_SectionLaw.hxx (L33-95) + BRepFill_SectionLaw.cxx (whole file
//! L28-112).
//!
//! Architecture differences:
//! - The abstract base class maps to the [`BRepFillSectionLawBase`] struct
//!   (the protected members myLaws / uclosed / vclosed / myDone / myIndices)
//!   plus the [`BRepFillSectionLawOps`] trait (the virtual part); the derived
//!   classes (BRepFill_ShapeLaw / BRepFill_NSections) embed the base and
//!   implement the trait (Rust has no inheritance).
//! - The private `BRepTools_WireExplorer myIterator` maps to the
//!   [`WireExplorerState`] cursor over the wire-ordered `TWireData::edges`
//!   list.
//! - `NCollection_DataMap<TopoDS_Shape, int>` maps to `HashMap` keyed by the
//!   `Shape::ptr_id()` identity (TopTools_ShapeMapHasher).
//! - `handle(GeomFill_SectionLaw)` maps to
//!   `Rc<RefCell<dyn geomalgo::geomfill SectionLaw>>`.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use rcad_kernel::geom::Curve3;
use rcad_kernel::topo::topods::GeomAbsShape;
use rcad_kernel::topo::topods::{BRep, Shape, TShape};

use crate::brep_fill::compatible_wires::wire_edges;
use crate::geomalgo::geomfill::section_law::SectionLaw;

/// OCCT BRepFill_SectionLaw protected members (hxx L85-90).
pub struct BRepFillSectionLawBase {
    /// OCCT myLaws (HArray1 of handle(GeomFill_SectionLaw)) — 1-based array
    /// mapped to a Vec (`Value(i)` -> `my_laws[i-1]`).
    pub my_laws: Vec<Rc<RefCell<dyn SectionLaw>>>,
    /// OCCT uclosed.
    pub uclosed: bool,
    /// OCCT vclosed.
    pub vclosed: bool,
    /// OCCT myDone.
    pub my_done: bool,
    /// OCCT myIndices (NCollection_DataMap<TopoDS_Shape, int>).
    pub my_indices: HashMap<u64, i32>,
    /// OCCT myIterator (BRepTools_WireExplorer, private).
    pub my_iterator: WireExplorerState,
}

/// OCCT BRepTools_WireExplorer — the wire-ordered edge cursor.
pub struct WireExplorerState {
    /// The explored wire (OCCT Init argument).
    pub wire: Shape,
    /// The wire-ordered edges (OCCT BRepTools_WireExplorer traversal).
    pub edges: Vec<Shape>,
    /// The cursor (0 = before Init / exhausted; 1-based like OCCT).
    pub cursor: usize,
}

impl WireExplorerState {
    /// OCCT BRepTools_WireExplorer::Init (BRepTools_WireExplorer.cxx
    /// L94-135, the wire-ordered reduction).
    pub fn init(&mut self, brep: &BRep, w: &Shape) {
        self.wire = w.clone();
        self.edges = wire_edges(brep, w);
        self.cursor = 0;
    }

    /// OCCT BRepTools_WireExplorer::More.
    pub fn more(&self) -> bool {
        self.cursor < self.edges.len()
    }

    /// OCCT BRepTools_WireExplorer::Next.
    pub fn next(&mut self) {
        self.cursor += 1;
    }

    /// OCCT BRepTools_WireExplorer::Current.
    pub fn current(&self) -> Shape {
        if self.cursor < self.edges.len() {
            self.edges[self.cursor].clone()
        } else {
            Shape::null()
        }
    }
}

impl BRepFillSectionLawBase {
    /// OCCT NbLaw (cxx L34-37) — Gives the number of elementary law.
    pub fn nb_law(&self) -> i32 {
        self.my_laws.len() as i32
    }

    /// OCCT Law (cxx L41-44).
    pub fn law(&self, index: i32) -> Rc<RefCell<dyn SectionLaw>> {
        self.my_laws[(index - 1) as usize].clone()
    }

    /// OCCT IndexOfEdge (cxx L48-51).
    pub fn index_of_edge(&self, an_edge: &Shape) -> i32 {
        self.my_indices.get(&an_edge.ptr_id()).copied().unwrap_or(0)
    }

    /// OCCT IsUClosed (cxx L55-58).
    pub fn is_uclosed(&self) -> bool {
        self.uclosed
    }

    /// OCCT IsVClosed (cxx L62-65).
    pub fn is_vclosed(&self) -> bool {
        self.vclosed
    }

    /// OCCT IsDone (cxx L69-72).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT Init (cxx L76-79).
    pub fn init(&mut self, brep: &BRep, w: &Shape) {
        self.my_iterator.init(brep, w);
    }

    /// OCCT CurrentEdge (cxx L85-111) — Parses the wire omitting the
    /// degenerated Edges.
    pub fn current_edge(&mut self, brep: &BRep) -> Shape {
        let mut e = Shape::null();
        let mut suivant = false;
        if self.my_iterator.more() {
            e = self.my_iterator.current();
            suivant = !e.is_null() && brep.edge(e.clone()).degenerated;
        }

        while suivant {
            self.my_iterator.next();
            e = self.my_iterator.current();
            suivant =
                !e.is_null() && brep.edge(e.clone()).degenerated && self.my_iterator.more();
        }

        if self.my_iterator.more() {
            self.my_iterator.next();
        }
        e
    }
}

/// OCCT BRepFill_SectionLaw — the virtual part (hxx L58-82) and the
/// `handle<BRepFill_SectionLaw>` slot trait (see the file-header note).
pub trait BRepFillSectionLawOps {
    /// The embedded base sub-object.
    fn base(&self) -> &BRepFillSectionLawBase;
    /// The embedded base sub-object (mutable).
    fn base_mut(&mut self) -> &mut BRepFillSectionLawBase;

    /// OCCT IsConstant — pure virtual: Say if the Law is Constant.
    fn is_constant(&self) -> bool;

    /// OCCT IsVertex — pure virtual: Say if the input shape is a vertex.
    fn is_vertex(&self) -> bool;

    /// OCCT ConcatenedLaw — pure virtual: the law built on a concatenated
    /// section (the OCCT null handle maps to None).  The `brep` pool carries
    /// the BRep_Tool::Curve reads of the derived implementations.
    fn concatened_law(&self, brep: &BRep) -> Option<Rc<RefCell<dyn SectionLaw>>>;

    /// OCCT Continuity — pure virtual.
    fn continuity(&self, brep: &BRep, index: i32, tol_angular: f64) -> GeomAbsShape;

    /// OCCT VertexTol — pure virtual.
    fn vertex_tol(&self, brep: &BRep, index: i32, param: f64) -> f64;

    /// OCCT Vertex — pure virtual.
    fn vertex(&self, brep: &mut BRep, index: i32, param: f64) -> Shape;

    /// OCCT D0 — pure virtual: Apply the Law to a shape, for a given
    /// parameter.
    fn d0(&mut self, brep: &mut BRep, u: f64, s: &mut Shape);
}

/// The section law evaluated from an edge curve — the shared
/// `BRep_Tool::Curve(E, f, l)` mapping used by the derived classes.
pub(super) fn brep_tool_curve(brep: &BRep, e: &Shape) -> Option<(Curve3, f64, f64)> {
    let ed = brep.edge(e.clone());
    let range = ed.range;
    ed.curve.clone().map(|c| (c, range[0], range[1]))
}

/// OCCT BRep_Tool::IsClosed(E) — the single-argument form resolves to
/// IsClosed(TopoDS_Shape); the EDGE branch (BRep_Tool.cxx L1751-1756) is
/// `TopExp::Vertices(E, V1, V2); return !V1.IsNull() && V1.IsSame(V2)`.
pub(super) fn brep_tool_is_closed_edge(brep: &BRep, e: &Shape) -> bool {
    let ed = brep.edge(e.clone());
    !ed.first.is_null() && ed.first.is_same(&ed.last)
}

/// OCCT Geom_TrimmedCurve knots/mults expansion for the rcad flat-knot
/// BSplineCurve3 (the Geom_BSplineCurve(Poles, Weights, Knots, Mults,
/// Degree, Periodic) construction mapping).
pub(super) fn expand_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut flat = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            flat.push(*k);
        }
    }
    flat
}

/// OCCT Geom_Curve::Reversed() over the rcad Curve3 — the basis-curve
/// reversal consumed by BRepFill_ShapeLaw::Init (cxx L151) and
/// BRepFill_Edge3DLaw / BRepFill_ACRLaw (the CBis forms).
///
/// Kernel GAP: only the Line / Circle / Ellipse / BSpline / Bezier / Trimmed
/// reversals are re-hosted (OCCT Geom_Line::Reversed = direction negation,
/// Geom_Conic::Reversed = axis reversal, Geom_BSplineCurve::Reverse =
/// BSplCLib::Reverse family, Geom_TrimmedCurve::Reverse = basis reversal +
/// bound swap per Geom_TrimmedCurve.cxx); the remaining geometric types keep
/// the OCCT failure path.
pub(super) fn reversed_curve(c: &Curve3) -> Curve3 {
    use rcad_kernel::geom::{
        BezierCurve3, BSplineCurve3, Circle3, Ellipse3, Line3, TrimmedCurve3,
    };
    match c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin,
            direction: -l.direction,
        }),
        Curve3::Circle(ci) => Curve3::Circle(Circle3 {
            center: ci.center,
            // gp_Ax2 reversal: the main direction flips, X keeps (the
            // parameterization runs the other way: u -> period - u).
            normal: -ci.normal,
            x_dir: ci.x_dir,
            y_dir: -ci.y_dir,
            radius: ci.radius,
        }),
        Curve3::Ellipse(el) => Curve3::Ellipse(Ellipse3 {
            center: el.center,
            normal: -el.normal,
            major_dir: el.major_dir,
            major_radius: el.major_radius,
            minor_radius: el.minor_radius,
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(bs.reversed()),
        Curve3::Bezier(bz) => Curve3::Bezier(BezierCurve3 {
            control_points: bz.control_points.iter().rev().copied().collect(),
            weights: {
                let mut w = bz.weights.clone();
                w.reverse();
                w
            },
        }),
        Curve3::Trimmed(t) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(reversed_curve(&t.curve)),
            first: reversed_parameter_of(&t.curve, t.last),
            last: reversed_parameter_of(&t.curve, t.first),
        }),
        _ => panic!(
            "GAP: Geom_*::Reversed (TKMath/gp kernel re-host) is not translated \
             for this curve type — BRepFill law reversal"
        ),
    }
}

/// OCCT Geom_Curve::ReversedParameter(U) — the rcad CurveEval mapping.
pub(super) fn reversed_parameter_of(c: &Curve3, u: f64) -> f64 {
    use rcad_kernel::geom::CurveEval;
    c.reversed_parameter(u)
}

/// GAP carrier: OCCT GeomConvert_CompCurveToBSplineCurve (TKGeomBase /
/// GeomConvert) — the 5-argument Add(NewCurve, Tol, After, IgnoreWarning,
/// NumRound) form consumed by BRepFill_ShapeLaw::ConcatenedLaw (cxx L392/395)
/// and BRepFill_NSections::totalsurf.  Not translated; all sites keep the
/// OCCT failure path.
pub(super) struct CompCurveToBSpline;

impl CompCurveToBSpline {
    /// OCCT GeomConvert_CompCurveToBSplineCurve(BasisCurve).
    pub fn new(_basis_curve: &Curve3) -> Self {
        panic!(
            "GAP: GeomConvert_CompCurveToBSplineCurve (TKGeomBase/GeomConvert) \
             is not translated — BRepFill ShapeLaw/NSections concatenation"
        )
    }

    /// OCCT Add(NewCurve, Tol, After, IgnoreWarning, NumRound).
    #[allow(clippy::too_many_arguments)]
    pub fn add(
        &mut self,
        _new_curve: &Curve3,
        _tol: f64,
        _after: bool,
        _ignore_warning: bool,
        _num_round: i32,
    ) -> bool {
        panic!(
            "GAP: GeomConvert_CompCurveToBSplineCurve::Add (TKGeomBase) \
             is not translated — BRepFill ShapeLaw/NSections concatenation"
        )
    }

    /// OCCT BSplineCurve().
    pub fn bspline_curve(&self) -> Curve3 {
        panic!(
            "GAP: GeomConvert_CompCurveToBSplineCurve::BSplineCurve (TKGeomBase) \
             is not translated — BRepFill ShapeLaw/NSections concatenation"
        )
    }
}
