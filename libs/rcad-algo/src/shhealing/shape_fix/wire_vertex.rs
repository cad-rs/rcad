//! 1:1 translation of OCCT `ShapeFix_WireVertex`
//! (`TKShHealing/ShapeFix/ShapeFix_WireVertex.hxx` L17-77 +
//! `ShapeFix_WireVertex.cxx` L1-293, docket row `ShapeFix_WireVertex`).
//!
//! Function-count equation (OCCT ShapeFix_WireVertex.cxx = rcad
//! wire_vertex.rs): `ShapeFix_WireVertex()` / `Init(wire, preci)` /
//! `Init(sbwd, preci)` / `Init(sawv)` / `Analyzer` / `WireData` / `Wire` /
//! `FixSame` / `Fix` — 9 OCCT functions = 9 rcad functions.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `BRep_Builder`/`ShapeAnalysis_Edge` read and
//!    write the TShape graph through `rcad_kernel::BRep`.
//! 2. `occ::handle<ShapeExtend_WireData>` — the
//!    `shape_extend::wire_data::WireData` value.
//! 3. `NCollection_HArray1<TopoDS_Shape>/<double>` — `Vec<Shape>` /
//!    `Vec<f64>` (1-based OCCT index -> index - 1).
//! 4. `V1 == V2` (TopoDS_Shape operator== / IsEqual) — same TShape +
//!    location + orientation (the shape_analysis/wire_vertex.rs bridge #5
//!    precedent).
//! 5. OCCT TopExp::Vertices — the TopExp.cxx re-host in the package
//!    statics (`shape_fix.rs`).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{tshape_flags, BRepBuilder, Orientation};
use rcad_kernel::BRep;

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::wire_vertex::ShapeAnalysisWireVertex;
use crate::shhealing::shape_build::brep_tool::{builder_add, set_flag_inplace};
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::shape_fix::top_exp_vertices;

/// OCCT TopoDS_Shape::operator== (IsEqual): same TShape + location, and the
/// same orientation.
fn shape_is_equal(a: &Shape, b: &Shape) -> bool {
    a.is_same(b) && a.orientation == b.orientation
}

/// OCCT TopoDS_Shape::Reversed: FORWARD <-> REVERSED (other orientations
/// unchanged).
fn reversed_orientation(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// OCCT ShapeFix_WireVertex (hxx L33-75): fixes disconnected edges in the
/// wire on the basis of the pre-analysis made by ShapeAnalysis_WireVertex.
pub struct ShapeFixWireVertex {
    /// OCCT myAnalyzer (hxx L74).
    my_analyzer: ShapeAnalysisWireVertex,
}

impl Default for ShapeFixWireVertex {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixWireVertex {
    /// OCCT ShapeFix_WireVertex::ShapeFix_WireVertex() (cxx L33).
    pub fn new() -> Self {
        ShapeFixWireVertex {
            my_analyzer: ShapeAnalysisWireVertex::new(),
        }
    }

    /// OCCT ShapeFix_WireVertex::Init(wire, preci) (cxx L37-41): loads the
    /// wire, initializes the internal analyzer with the given precision and
    /// performs the analysis.
    pub fn init(&mut self, brep: &mut BRep, wire: &Shape, preci: f64) {
        // L39: sbwd = new ShapeExtend_WireData(wire) (the hxx defaults
        // Chained = theManifold = Standard_True).
        let sbwd = WireData::new_from_wire(brep, wire, true, true);
        // L40.
        self.init_wire_data(brep, sbwd, preci);
    }

    /// OCCT ShapeFix_WireVertex::Init(sbwd, preci) (cxx L45-50).
    pub fn init_wire_data(&mut self, brep: &mut BRep, sbwd: WireData, preci: f64) {
        // L47-49.
        self.my_analyzer.load_wire_data(sbwd);
        self.my_analyzer.set_precision(preci);
        self.my_analyzer.analyze(brep);
    }

    /// OCCT ShapeFix_WireVertex::Init(sawv) (cxx L54-57): loads all the data
    /// on a wire already analysed by ShapeAnalysis_WireVertex.
    pub fn init_analyzer(&mut self, sawv: ShapeAnalysisWireVertex) {
        // L56.
        self.my_analyzer = sawv;
    }

    /// OCCT ShapeFix_WireVertex::Analyzer (cxx L61-64): the internal
    /// analyzer.
    pub fn analyzer(&self) -> &ShapeAnalysisWireVertex {
        &self.my_analyzer
    }

    /// OCCT ShapeFix_WireVertex::WireData (cxx L68-71): the data on the wire
    /// (fixed).
    pub fn wire_data(&self) -> Option<&WireData> {
        // L70.
        self.my_analyzer.wire_data()
    }

    /// OCCT ShapeFix_WireVertex::Wire (cxx L75-78): the resulting wire
    /// (fixed).
    pub fn wire(&self, brep: &mut BRep) -> Shape {
        match self.my_analyzer.wire_data() {
            Some(sbwd) => sbwd.wire(brep),
            None => Shape::null(),
        }
    }

    /// OCCT ShapeFix_WireVertex::FixSame (cxx L82-137): fixes the "Same" or
    /// "Close" statuses (the same vertex may be set without changing the
    /// parameters); returns the count of fixed vertices.
    pub fn fix_same(&mut self, brep: &mut BRep) -> i32 {
        //  FixSame : prend les status "SameCoord" et "Close" et les force a
        //  "Same"  (the OCCT comments stay untranslated)
        // L86-89.
        if !self.my_analyzer.is_done() {
            return 0;
        }

        // L91-92.
        let mut nbfix = 0;
        let mut b = BRepBuilder::new();

        // L94.
        let nb = self
            .my_analyzer
            .wire_data()
            .map(|w| w.nb_edges())
            .unwrap_or(0);

        // L97-135.
        for i in 1..=nb {
            // L99.
            let j = if i == nb { 1 } else { i + 1 };
            // L100.
            let stat = self.my_analyzer.status(i);
            // L101-104.
            if stat != 1 && stat != 2 {
                continue;
            }
            // L105-107: Ici on prend un vertex et on le generalise aux deux
            // edges.
            let e1 = self
                .my_analyzer
                .wire_data()
                .map(|w| w.edge(i))
                .unwrap_or_else(Shape::null);
            let e2 = self
                .my_analyzer
                .wire_data()
                .map(|w| w.edge(j))
                .unwrap_or_else(Shape::null);

            // L109-111.
            let sae = ShapeAnalysisEdge::new();
            let v1 = sae.last_vertex(brep, &e1);
            let v2 = sae.first_vertex(brep, &e2);
            // L112-116: deja fait ... (the OCCT comment stays untranslated).
            if shape_is_equal(&v1, &v2) {
                self.my_analyzer.set_same_vertex(i);
                continue;
            }
            // L117-126.
            if stat == 2 {
                // OK mais en reprenant les tolerances (the OCCT comment
                // stays untranslated).
                let mut crv: Option<rcad_kernel::geom::Curve3> = None;
                let mut cf = 0.0;
                let mut cl = 0.0;
                // L122.
                sae.curve3d(brep, &e1, &mut crv, &mut cf, &mut cl, true);
                // L123.
                let preci = self.my_analyzer.precision();
                b.update_vertex_on_edge(brep, v1.clone(), cl, e1.clone(), preci);
                // L124-125.
                sae.curve3d(brep, &e2, &mut crv, &mut cf, &mut cl, true);
                b.update_vertex_on_edge(brep, v1.clone(), cf, e2.clone(), preci);
            }
            // L127-132: Et remettre ce vtx en commun.
            let mut v1m = v1.clone();
            v1m.orientation = e2.orientation;
            builder_add(brep, &e2, &v1m);
            v1m.orientation = reversed_orientation(e1.orientation);
            builder_add(brep, &e1, &v1m);
            // L133: conclusion.
            self.my_analyzer.set_same_vertex(i);
            // L134.
            nbfix += 1;
        }
        // L136.
        nbfix
    }

    /// OCCT ShapeFix_WireVertex::Fix (cxx L141-293): fixes all statuses
    /// except "Disjoined" — the cases in which a common value has been set,
    /// with or without changing parameters; returns the count of fixed
    /// vertices.
    pub fn fix(&mut self, brep: &mut BRep) -> i32 {
        //  Ici le grand jeu : on repasse partout  (the OCCT comments stay
        //  untranslated)
        // L151-154.
        if !self.my_analyzer.is_done() {
            return 0;
        }

        // L156-159.
        let nb = self
            .my_analyzer
            .wire_data()
            .map(|w| w.nb_edges())
            .unwrap_or(0);
        // L160-168: the first pass — count the statuses to fix (the
        // szv#4:S4163:12Mar99 optimization).
        let mut nbfix = 0;
        for i in 1..=nb {
            if self.my_analyzer.status(i) > 0 {
                nbfix += 1;
            }
        }
        // L169-172.
        if nbfix == 0 {
            return 0;
        }

        // L174: BRep_Builder B.
        let mut b = BRepBuilder::new();

        // L174-180: the per-edge arrays (VI / VJ / EF / UI / UJ).
        let mut vi: Vec<Shape> = vec![Shape::null(); nb as usize];
        let mut vj: Vec<Shape> = vec![Shape::null(); nb as usize];
        let mut ef: Vec<Shape> = vec![Shape::null(); nb as usize];
        let mut ui: Vec<f64> = vec![0.0; nb as usize];
        let mut uj: Vec<f64> = vec![0.0; nb as usize];

        // L182-221: the value-collection pass.
        for i in 1..=nb {
            // L185.
            let j = if i == nb { 1 } else { i + 1 };
            // L186.
            let stat = self.my_analyzer.status(i);

            // L188-192.
            let sae = ShapeAnalysisEdge::new();
            let e_i = self
                .my_analyzer
                .wire_data()
                .map(|w| w.edge(i))
                .unwrap_or_else(Shape::null);
            let e_j = self
                .my_analyzer
                .wire_data()
                .map(|w| w.edge(j))
                .unwrap_or_else(Shape::null);
            let v1 = sae.last_vertex(brep, &e_i);
            let v2 = sae.first_vertex(brep, &e_j);
            vi[(i - 1) as usize] = v1;
            vj[(j - 1) as usize] = v2;

            // L194-196: (the OCCT comment "E.EmptyCopy(); trop d ennuis"
            // stays untranslated).
            ef[(i - 1) as usize] = e_i.clone();

            // L201-202.
            let mut upre = self.my_analyzer.uprevious(i);
            let mut ufol = self.my_analyzer.ufollowing(j);

            // L204-216: szv#4:S4163:12Mar99 optimized.
            let mut crv: Option<rcad_kernel::geom::Curve3> = None;
            let mut cf = 0.0;
            let mut cl = 0.0;
            if stat < 4 {
                sae.curve3d(brep, &e_i, &mut crv, &mut cf, &mut cl, true);
                upre = cl;
            }
            if stat < 3 || stat == 4 {
                sae.curve3d(brep, &e_j, &mut crv, &mut cf, &mut cl, true);
                ufol = cf;
            }

            // L218-219: nbfix ++ (the commented-out OCCT line).
            ui[(i - 1) as usize] = upre;
            uj[(j - 1) as usize] = ufol;
        }

        // L223-226.
        if nbfix == 0 {
            return nbfix;
        }

        // L228-242: the edges are stripped of their former vertices (the
        // OCCT comments stay untranslated).
        for i in 1..=nb {
            // L235-237.
            let mut e1 = ef[(i - 1) as usize].clone();
            e1.orientation = Orientation::Forward;
            // L238.
            let (va, vb) = top_exp_vertices(brep, &e1);
            // L239.
            set_flag_inplace(brep, &e1, tshape_flags::FREE, true);
            // L240-241.
            b.remove_from_edge(brep, e1.clone(), va);
            b.remove_from_edge(brep, e1, vb);
        }

        // L244-284: the fixing pass.
        let prec = self.my_analyzer.precision();
        for i in 1..=nb {
            // L249-250.
            let j = if i == nb { 1 } else { i + 1 };
            // L250.
            let stat = self.my_analyzer.status(i);

            // L253-258.
            let mut v1 = vi[(i - 1) as usize].clone();
            let v2 = vj[(j - 1) as usize].clone();
            let mut e1 = ef[(i - 1) as usize].clone();
            let e2 = ef[(j - 1) as usize].clone();
            let upre = ui[(i - 1) as usize];
            let ufol = uj[(j - 1) as usize];

            // L260-263: Changer les coords ? (the OCCT comment stays
            // untranslated).
            if stat > 2 {
                b.update_vertex_point(brep, v1.clone(), self.my_analyzer.position(i), prec);
            }

            // L266-272: ce qui suit : seulement si vertex a reprendre (the
            // OCCT comments stay untranslated).
            if stat > 0 {
                b.update_vertex_on_edge(brep, v1.clone(), upre, e1.clone(), prec);
                b.update_vertex_on_edge(brep, v1.clone(), ufol, e2.clone(), prec);
                // L270.
                v1.orientation = Orientation::Forward;
                // L271: V1.Orientation (E2.Orientation());
            }

            // L274-280: Comme on a deshabille les edges, il faut tout
            // remettre (the OCCT comments stay untranslated).
            // L275: sur place.
            set_flag_inplace(brep, &e2, tshape_flags::FREE, true);
            // L276.
            builder_add(brep, &e2, &v1);
            // L277.
            v1.orientation = Orientation::Reversed;
            // L278: V1.Orientation (E1.Orientation()); V1.Reverse();
            // L279: sur place.
            set_flag_inplace(brep, &e1, tshape_flags::FREE, true);
            // L280.
            builder_add(brep, &e1, &v1);

            // L282: conclusion.
            self.my_analyzer.set_same_vertex(i);
            // The OCCT V2 binding has no later reference (L254).
            let _ = v2;
        }

        // L286-290: pour finir, MAJ du STW.
        for i in 1..=nb {
            let e = ef[(i - 1) as usize].clone();
            if let Some(sbwd) = self.my_analyzer.wire_data_mut() {
                sbwd.set_edge(&e, i);
            }
        }

        // L292.
        nbfix
    }
}
