//! OCCT BRepFill_SectionPlacement (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_SectionPlacement.hxx (L38-66) + BRepFill_SectionPlacement.cxx
//! (whole file L45-280, including the static SearchParam).
//!
//! Architecture differences:
//! - `handle(BRepFill_LocationLaw) myLaw` maps to
//!   `Rc<RefCell<dyn BRepFillLocationLawOps>>` (the derived-law handle
//!   slot; see the trait note in brep_fill_location_law.rs).
//! - `handle(Geom_Geom) aSection` (a Geom_CartesianPoint or a Geom_Curve)
//!   maps to the local [`SectionGeometry`] enum.
//! - `BRep_Tool::Curve(E, f, l)` maps to `TEdgeData::curve + range`; the
//!   owning `BRep` pool is passed as the leading `brep` argument.
//! - `TopExp_Explorer(S, TopAbs_EDGE / TopAbs_VERTEX)` maps to
//!   `feat::brep_feat_builder::explorer`.
//! - `BRepAdaptor_CompCurve` (TKBRep adaptor) is a GAP carrier
//!   ([`BRepAdaptorCompCurve`]) — the comp-curve adaptor is not translated;
//!   the construction site (cxx L168) keeps the OCCT failure path.
//! - `GeomConvert_CompCurveToBSplineCurve` reuses the offset-package carrier
//!   (offset::brep_offset_inter2d::GeomConvertCompCurveToBSplineCurve).
//! - The GeomFill engine is the batch-2
//!   `geomalgo::geomfill::section_placement::SectionPlacement`
//!   (`Perform(Tol)` -> `perform_confusion`, `Perform(Path, Tol)` ->
//!   `perform_with_path`, `Perform(ParamOnPath, Tol)` -> `perform_at`).

use std::cell::RefCell;
use std::rc::Rc;

use glam::{DVec3, DAffine3};

use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::geom::{Curve3, TrimmedCurve3};
use rcad_kernel::topo::topods::{BRep, Shape, TShape};

use crate::feat::brep_feat_builder::explorer;
use crate::geomalgo::geomfill::section_placement::SectionPlacement;

use crate::brep_fill::brep_fill_location_law::BRepFillLocationLawOps;
use crate::brep_fill::compatible_wires::{
    brep_tool_parameter, shape_tolerance,
};
use crate::brep_fill::generator::top_exp_vertices;
use crate::offset::brep_offset_inter2d::GeomConvertCompCurveToBSplineCurve;

/// OCCT static SearchParam (BRepFill_SectionPlacement.cxx L45-62).
fn search_param(
    brep: &BRep,
    law: &Rc<RefCell<dyn BRepFillLocationLawOps>>,
    ind: i32,
    the_v: &Shape,
) -> f64 {
    let mut t;
    let e = law.borrow().base().edge(ind);
    // BRep_Tool::Parameter(TheV, E)
    t = brep_tool_parameter(brep, the_v, &e);
    if e.orientation == rcad_kernel::topo::topods::Orientation::Reversed {
        // Lf = Law->Law(Ind)->GetCurve()->FirstParameter();
        // Ll = Law->Law(Ind)->GetCurve()->LastParameter();
        let curve = law.borrow().base().law(ind).borrow().get_curve();
        if let Some(curve) = curve {
            let f = brep_tool_range_of(brep, &e).0;
            let l = brep_tool_range_of(brep, &e).1;
            let lf = crate::geomalgo::geomfill::trihedron_law::curve_first_parameter(&curve);
            let ll = crate::geomalgo::geomfill::trihedron_law::curve_last_parameter(&curve);
            // t = Ll - (t - f) * (Ll - Lf) / (l - f);
            t = ll - (t - f) * (ll - lf) / (l - f);
        }
    }
    t
}

/// OCCT BRep_Tool::Range(E, f, l) — the stored edge range.
fn brep_tool_range_of(brep: &BRep, e: &Shape) -> (f64, f64) {
    let r = brep.edge(e.clone()).range;
    (r[0], r[1])
}

/// OCCT BRep_Tool::Curve(E, f, l) — the (curve, first, last) triple.
fn brep_tool_curve(brep: &BRep, e: &Shape) -> Option<(Curve3, f64, f64)> {
    let ed = brep.edge(e.clone());
    let range = ed.range;
    ed.curve.clone().map(|c| (c, range[0], range[1]))
}

/// OCCT Geom_TrimmedCurve(C, f, l).
fn trimmed_curve(c: &Curve3, f: f64, l: f64) -> Curve3 {
    Curve3::Trimmed(TrimmedCurve3 {
        curve: Box::new(c.clone()),
        first: f,
        last: l,
    })
}

/// OCCT handle(Geom_Geom) aSection — either a Geom_CartesianPoint (the
/// vertex case) or a Geom_Curve.
#[allow(clippy::large_enum_variant)]
enum SectionGeometry {
    /// OCCT Geom_CartesianPoint(P).
    Point(DVec3),
    /// OCCT handle(Geom_Curve).
    Curve(Curve3),
}

/// OCCT BRepFill_SectionPlacement (hxx L38-66).
pub struct BRepFillSectionPlacement {
    /// OCCT handle(BRepFill_LocationLaw) myLaw.
    pub my_law: Rc<RefCell<dyn BRepFillLocationLawOps>>,
    /// OCCT TopoDS_Shape mySection.
    pub my_section: Shape,
    /// OCCT gp_Trsf myTrsf.
    pub my_trsf: DAffine3,
    /// OCCT int myIndex.
    pub my_index: i32,
    /// OCCT double myParam.
    pub my_param: f64,
}

impl BRepFillSectionPlacement {
    /// OCCT BRepFill_SectionPlacement(Law, Section, WithContact,
    /// WithCorrection) (cxx L64-74).
    pub fn new(
        brep: &mut BRep,
        law: Rc<RefCell<dyn BRepFillLocationLawOps>>,
        section: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) -> Self {
        let mut place = BRepFillSectionPlacement {
            my_law: law,
            my_section: section.clone(),
            my_trsf: DAffine3::IDENTITY,
            my_index: 0,
            my_param: 0.0,
        };
        // TopoDS_Vertex VNull; VNull.Nullify();
        let v_null = Shape::null();
        place.perform(brep, with_contact, with_correction, &v_null);
        place
    }

    /// OCCT BRepFill_SectionPlacement(Law, Section, Vertex, WithContact,
    /// WithCorrection) (cxx L76-85).
    pub fn new_with_vertex(
        brep: &mut BRep,
        law: Rc<RefCell<dyn BRepFillLocationLawOps>>,
        section: &Shape,
        vertex: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) -> Self {
        let mut place = BRepFillSectionPlacement {
            my_law: law,
            my_section: section.clone(),
            my_trsf: DAffine3::IDENTITY,
            my_index: 0,
            my_param: 0.0,
        };
        place.perform(brep, with_contact, with_correction, vertex);
        place
    }

    /// OCCT Perform (cxx L87-247).
    pub fn perform(
        &mut self,
        brep: &mut BRep,
        with_contact: bool,
        with_correction: bool,
        vertex: &Shape,
    ) {
        // Here we are simply looking for the first valid curve in the section.
        // OCCT TopExp_Explorer anEdgeExplorer(mySection, TopAbs_EDGE)
        let section_edges = explorer(&self.my_section, rcad_kernel::topo::topods::ShapeType::Edge, rcad_kernel::topo::topods::ShapeType::Shape);

        // OCCT walks the explorer to the first edge carrying a non-null
        // 3d curve (cxx L95-110); the explorer cursor position is carried
        // by the (start index) below.
        let mut a_curve: Option<(Curve3, f64, f64)> = None;
        let mut start = section_edges.len();
        for (idx, an_edge) in section_edges.iter().enumerate() {
            if an_edge.is_null() || brep.edge(an_edge.clone()).degenerated {
                continue;
            }
            a_curve = brep_tool_curve(brep, an_edge);
            if a_curve.is_none() {
                continue;
            }
            start = idx;
            break;
        }

        let a_section: SectionGeometry;
        if a_curve.is_none() {
            // No edge found : the section is a vertex
            // OCCT TopExp_Explorer aVertexExplorer(mySection, TopAbs_VERTEX)
            let section_vertices = explorer(&self.my_section, rcad_kernel::topo::topods::ShapeType::Vertex, rcad_kernel::topo::topods::ShapeType::Shape);
            let a_first_vertex = &section_vertices[0];
            let a_point = brep.vertex(a_first_vertex.clone()).point;
            a_section = SectionGeometry::Point(a_point);
        } else {
            let (c, f, l) = a_curve.expect("checked");
            // Geom_TrimmedCurve aTrimmedCurve(aCurve, f, l)
            let mut a_trimmed_curve = trimmed_curve(&c, f, l);
            let mut an_edge_start_param;
            let mut an_edge_end_param;

            // OCCT anEdgeExplorer.Next(); if (anEdgeExplorer.More())
            let has_more = start + 1 < section_edges.len();
            if has_more {
                // GeomConvert_CompCurveToBSplineCurve aBSplineConverter(aTrimmedCurve)
                let mut a_bspline_converter = GeomConvertCompCurveToBSplineCurve::new(&a_trimmed_curve);
                for an_edge in &section_edges[start + 1..] {
                    // avoid null, degenerated edges
                    if an_edge.is_null() || brep.edge(an_edge.clone()).degenerated {
                        continue;
                    }

                    let (c, f0, l0) = match brep_tool_curve(brep, an_edge) {
                        Some(v) => v,
                        None => continue,
                    };
                    an_edge_start_param = f0;
                    an_edge_end_param = l0;

                    // TopExp::Vertices(anEdge, aFirstVertex, aLastVertex);
                    let (a_first_vertex, a_last_vertex) = top_exp_vertices(brep, an_edge);
                    let a_vertex_tolerance = shape_tolerance(brep, &a_first_vertex)
                        .max(shape_tolerance(brep, &a_last_vertex));

                    a_trimmed_curve = trimmed_curve(&c, an_edge_start_param, an_edge_end_param);
                    // Add(TC, min(aPrecisionTolerance, aVertexTolerance));
                    // fallback Add(TC, max(...)).
                    if !a_bspline_converter.add(
                        &a_trimmed_curve,
                        CONFUSION.min(a_vertex_tolerance),
                    ) {
                        a_bspline_converter.add(
                            &a_trimmed_curve,
                            CONFUSION.max(a_vertex_tolerance),
                        );
                    }
                }
                a_curve = Some((a_bspline_converter.bspline_curve(), 0.0, 0.0));
            } else {
                a_curve = Some((a_trimmed_curve.clone(), 0.0, 0.0));
            }

            a_section = SectionGeometry::Curve(a_curve.take().expect("checked").0);
        }

        // GeomFill_SectionPlacement aSectionPlacement(myLaw->Law(1), aSection);
        let mut a_section_placement = match &a_section {
            SectionGeometry::Curve(c) => SectionPlacement::new_loc(
                self.my_law.borrow().base().law(1),
                c,
            ),
            SectionGeometry::Point(p) => {
                SectionPlacement::new_point(self.my_law.borrow().base().law(1), *p)
            }
        };

        // occ::handle<BRepAdaptor_CompCurve> aWireAdaptor =
        //   new BRepAdaptor_CompCurve(myLaw->Wire());
        // aSectionPlacement.Perform(aWireAdaptor, Precision::Confusion());
        let a_wire_adaptor = BRepAdaptorCompCurve::new(brep, &self.my_law.borrow().base().wire());
        a_section_placement.perform_with_path(Some(a_wire_adaptor.comp_curve()), CONFUSION);

        let a_section_param = a_section_placement.parameter_on_path();
        let a_param_confusion = PCONFUSION;

        let mut a_law_index1 = 0i32;
        let mut a_law_index2 = 0i32;
        // In the general case : Localisation via concatenation of the spine
        let mut an_is_interval_found = false;
        let mut a_law_index = 1i32;
        while a_law_index <= self.my_law.borrow().base().nb_law() && !an_is_interval_found {
            let a_curr_knot_param = (a_law_index - 1) as f64;
            let a_next_knot_param = a_law_index as f64;

            // Check if the section parameter is in the interval
            an_is_interval_found = (a_curr_knot_param - a_param_confusion <= a_section_param)
                && (a_next_knot_param + a_param_confusion >= a_section_param);
            if !an_is_interval_found {
                a_law_index += 1;
                continue;
            }

            a_law_index1 = a_law_index;
            if (a_section_param - a_curr_knot_param).abs() < a_param_confusion && a_law_index > 1
            {
                a_law_index2 = a_law_index - 1;
            } else if (a_section_param - a_next_knot_param).abs() < a_param_confusion
                && a_law_index < self.my_law.borrow().base().nb_law()
            {
                a_law_index2 = a_law_index + 1;
            }
            a_law_index += 1;
        }

        if !an_is_interval_found {
            // throw Standard_ConstructionError("Interval is not found")
            panic!("Standard_ConstructionError: Interval is not found");
        }

        // Search of the <Ind1> by vertex <TheV>
        let mut an_is_vertex_on_law = false;
        let mut a_vertex = vertex.clone(); // TopoDS::Vertex(Vertex) cast
        if !a_vertex.is_null() {
            for a_current_law_index in 1..=self.my_law.borrow().base().nb_law() {
                let an_edge = self.my_law.borrow().base().edge(a_current_law_index);
                let (v1, v2) = top_exp_vertices(brep, &an_edge);
                if v1.ptr_id() == a_vertex.ptr_id() || v2.ptr_id() == a_vertex.ptr_id() {
                    an_is_vertex_on_law = true;
                    a_law_index1 = a_current_law_index;
                    a_law_index2 = 0;
                    break;
                }
            }
        }

        // Positioning on the localized edge (or 2 Edges)
        a_section_placement.set_location(self.my_law.borrow().base().law(a_law_index1));
        if !an_is_vertex_on_law {
            a_section_placement.perform_confusion(CONFUSION);
        } else {
            a_section_placement.perform_at(
                search_param(brep, &self.my_law, a_law_index1, &a_vertex),
                CONFUSION,
            );
        }

        self.my_trsf = a_section_placement.transformation(with_contact, with_correction);
        self.my_index = a_law_index1;
        self.my_param = a_section_placement.parameter_on_path();
        let angle = a_section_placement.angle();

        if a_law_index2 != 0 {
            a_section_placement.set_location(self.my_law.borrow().base().law(a_law_index2));
            if !an_is_vertex_on_law {
                a_section_placement.perform_confusion(CONFUSION);
            } else {
                if a_law_index1 == a_law_index2 {
                    // aVertex.Reverse()
                    let reversed = crate::brep_fill::generator::shape_reversed(&a_vertex);
                    a_vertex = reversed;
                }
                a_section_placement.perform_at(
                    search_param(brep, &self.my_law, a_law_index2, &a_vertex),
                    CONFUSION,
                );
            }
            if a_section_placement.angle() > angle {
                self.my_trsf = a_section_placement.transformation(with_contact, with_correction);
                self.my_index = a_law_index2;
                self.my_param = a_section_placement.parameter_on_path();
            }
        }
    }

    /// OCCT Transformation (cxx L249-252) — the const accessor.
    pub fn transformation(&self) -> DAffine3 {
        self.my_trsf
    }

    /// OCCT AbscissaOnPath (cxx L254-257).
    pub fn abscissa_on_path(&mut self) -> f64 {
        self.my_law
            .borrow_mut()
            .base_mut()
            .abscissa(self.my_index, self.my_param)
    }
}

// ---------------------------------------------------------------------------
// GAP carrier
// ---------------------------------------------------------------------------

/// GAP carrier: OCCT BRepAdaptor_CompCurve (TKBRep/BRepAdaptor,
/// BRepAdaptor_CompCurve.cxx) — the comp-curve adaptor over a wire is not
/// translated (plan D3); the construction site
/// (BRepFill_SectionPlacement.cxx L168) keeps the OCCT failure path.
pub struct BRepAdaptorCompCurve;

impl BRepAdaptorCompCurve {
    /// OCCT new BRepAdaptor_CompCurve(W).
    pub fn new(_brep: &BRep, _w: &Shape) -> Self {
        panic!(
            "GAP: BRepAdaptor_CompCurve (TKBRep) is not translated — \
             BRepFill_SectionPlacement::Perform (BRepFill_SectionPlacement.cxx L168)"
        )
    }

    /// The comp-curve as the rcad Curve3 view consumed by
    /// GeomFill_SectionPlacement::Perform(Path, Tol).
    pub fn comp_curve(&self) -> Curve3 {
        panic!(
            "GAP: BRepAdaptor_CompCurve (TKBRep) is not translated — \
             BRepFill_SectionPlacement::Perform (BRepFill_SectionPlacement.cxx L168)"
        )
    }

    /// OCCT BRepAdaptor_CompCurve::D1(U, P, V) — the point/derivative form
    /// consumed by BRepFill_PipeShell::Set(AuxiliarySpine, ...)
    /// (BRepFill_PipeShell.cxx L392-393).
    pub fn d1(&self, _u: f64, _p: &mut DVec3, _v: &mut DVec3) {
        panic!(
            "GAP: BRepAdaptor_CompCurve (TKBRep) is not translated — \
             BRepFill_PipeShell::Set (BRepFill_PipeShell.cxx L392)"
        )
    }
}
