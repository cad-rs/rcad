//! OCCT BRepFill_NSections (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_NSections.hxx (L33-96) + BRepFill_NSections.cxx (whole file
//! L58-882, including the statics EdgeToBSpline / totalsurf).
//!
//! Architecture differences:
//! - The BRepFill_SectionLaw inheritance maps to the base struct +
//!   ops-trait form of brep_fill_section_law.rs (see its header).
//! - `NCollection_Sequence<TopoDS_Shape>` maps to `Vec<Shape>`;
//!   `NCollection_Sequence<double>` / `<gp_Trsf>` to `Vec<f64>` / `Vec<Trsf>`.
//! - `BRepTools_WireExplorer` maps to the wire-ordered `TWireData::edges`;
//!   `BRep_Tool::Curve(E, f, l)` maps to `TEdgeData::curve + range`.
//! - `GeomConvert_ApproxCurve` (TKGeomBase/GeomConvert) is a GAP carrier
//!   ([`GeomConvertApproxCurve`]) — HasResult keeps the OCCT failure path
//!   (false) so the literal GeomConvert::CurveToBSplineCurve fallback is
//!   taken (the rcad kernel convert re-host).
//! - `GeomFill_SectionGenerator::AddCurve / Perform` (the producer half of
//!   the batch-1 SectionGenerator carrier) and `GeomFill_AppSurf`
//!   (TKGeomAlgo/AppBlend, the staged AppBlend gap) are GAP carriers in
//!   [`totalsurf`] — the sites keep the OCCT failure path.
//! - `Geom_BSplineSurface::VIso` maps to the local
//!   [`bspline_surface_viso`] re-host (the homogeneous De Boor evaluation of
//!   the opposite-direction basis, the geomfill/nsections.rs precedent).

use std::cell::RefCell;
use std::rc::Rc;

use glam::DVec3;

use rcad_kernel::base::convert::curve_to_bspline;
use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::geom::{CurveEval, BSplineCurve3, BSplineSurface, Curve3, TrimmedCurve3};
use rcad_kernel::topo::topods::GeomAbsShape;
use rcad_kernel::math::bspl_lib::reparametrize;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, Shape};

use crate::brep_fill::brep_fill_section_law::{
    brep_tool_curve, expand_knots, BRepFillSectionLawBase, BRepFillSectionLawOps,
    WireExplorerState,
};
use crate::brep_fill::compatible_wires::wire_edges;
use crate::brep_fill::generator::top_exp_vertices;
use crate::geomalgo::geomfill::line::Line;
use crate::geomalgo::geomfill::nsections::NSections;
use crate::geomalgo::geomfill::section_generator::SectionGenerator;
use crate::geomalgo::geomfill::section_law::SectionLaw;
use crate::geomalgo::geomfill::trihedron_law::{curve_first_parameter, curve_last_parameter};

/// OCCT gp_Trsf (the myTrsfs element type; the rcad math gp carrier).
use rcad_kernel::math::gp::Trsf;

// ---------------------------------------------------------------------------
// GAP carriers
// ---------------------------------------------------------------------------

/// GAP carrier: OCCT GeomConvert_ApproxCurve (TKGeomBase/GeomConvert,
/// GeomConvert_ApproxCurve.cxx) — the curve approximator is not translated;
/// HasResult keeps the OCCT failure path (false) so the literal fallback
/// GeomConvert::CurveToBSplineCurve (the kernel `curve_to_bspline` re-host)
/// is taken at the EdgeToBSpline call site (cxx L100-109).
struct GeomConvertApproxCurve;

impl GeomConvertApproxCurve {
    /// OCCT GeomConvert_ApproxCurve(Curve, Tol3d, Order, MaxSegments,
    /// MaxDegree).
    #[allow(clippy::too_many_arguments)]
    fn new(_curve: &Curve3, _tol3d: f64, _order: GeomAbsShape, _max_segments: i32, _max_degree: i32) -> Self {
        GeomConvertApproxCurve
    }

    /// OCCT HasResult().
    #[allow(dead_code)]
    fn has_result(&self) -> bool {
        // GAP: the approximator is not translated — no result.
        false
    }

    /// OCCT Curve().
    #[allow(dead_code)]
    fn curve(&self) -> Option<BSplineCurve3> {
        None
    }
}

/// GAP carrier: OCCT GeomFill_AppSurf (TKGeomAlgo/GeomFill, deriving
/// AppBlend_AppSurf) — the AppBlend approximation family is the staged gap
/// (the geomfill/sweep.rs ApproxSweepApproximation precedent); the
/// construction / Perform / read-back sites in totalsurf (cxx L258-263,
/// L266-355) keep the OCCT failure path.
struct GeomFillAppSurf;

impl GeomFillAppSurf {
    /// OCCT GeomFill_AppSurf(Degmin, Degmax, Pres3d, Pres3d, NbIt, KnownP).
    #[allow(clippy::too_many_arguments)]
    fn new(_degmin: i32, _degmax: i32, _pres3d: f64, _pres3d2: f64, _nb_it: i32, _known_p: bool) -> Self {
        panic!(
            "GAP: GeomFill_AppSurf (TKGeomAlgo/AppBlend) is not translated — \
             BRepFill_NSections::totalsurf (BRepFill_NSections.cxx L261)"
        )
    }

    /// OCCT Perform(Line, Section, SpApprox).
    fn perform(&mut self, _line: &Line, _section: &mut SectionGenerator, _sp_approx: bool) {
        panic!(
            "GAP: GeomFill_AppSurf::Perform (TKGeomAlgo/AppBlend) is not translated — \
             BRepFill_NSections::totalsurf (BRepFill_NSections.cxx L263)"
        )
    }

    /// OCCT SurfPoles().
    #[allow(dead_code)]
    fn surf_poles(&self) -> Vec<Vec<DVec3>> {
        panic!("GAP: GeomFill_AppSurf::SurfPoles (TKGeomAlgo/AppBlend)")
    }

    /// OCCT SurfWeights().
    #[allow(dead_code)]
    fn surf_weights(&self) -> Vec<Vec<f64>> {
        panic!("GAP: GeomFill_AppSurf::SurfWeights (TKGeomAlgo/AppBlend)")
    }

    /// OCCT SurfUKnots().
    #[allow(dead_code)]
    fn surf_u_knots(&self) -> Vec<f64> {
        panic!("GAP: GeomFill_AppSurf::SurfUKnots (TKGeomAlgo/AppBlend)")
    }

    /// OCCT SurfVKnots().
    #[allow(dead_code)]
    fn surf_v_knots(&self) -> Vec<f64> {
        panic!("GAP: GeomFill_AppSurf::SurfVKnots (TKGeomAlgo/AppBlend)")
    }

    /// OCCT SurfUMults().
    #[allow(dead_code)]
    fn surf_u_mults(&self) -> Vec<i32> {
        panic!("GAP: GeomFill_AppSurf::SurfUMults (TKGeomAlgo/AppBlend)")
    }

    /// OCCT SurfVMults().
    #[allow(dead_code)]
    fn surf_v_mults(&self) -> Vec<i32> {
        panic!("GAP: GeomFill_AppSurf::SurfVMults (TKGeomAlgo/AppBlend)")
    }

    /// OCCT UDegree().
    #[allow(dead_code)]
    fn u_degree(&self) -> i32 {
        panic!("GAP: GeomFill_AppSurf::UDegree (TKGeomAlgo/AppBlend)")
    }

    /// OCCT VDegree().
    #[allow(dead_code)]
    fn v_degree(&self) -> i32 {
        panic!("GAP: GeomFill_AppSurf::VDegree (TKGeomAlgo/AppBlend)")
    }
}

/// GAP: GeomFill_SectionGenerator::AddCurve — the producer half of the
/// batch-1 SectionGenerator carrier (the AppBlend-consumption accessors are
/// translated); the site (cxx L167 / L223 / L247) keeps the OCCT failure
/// path.
fn gap_section_generator_add_curve(_section: &mut SectionGenerator, _curve: &Curve3) {
    panic!(
        "GAP: GeomFill_SectionGenerator::AddCurve (TKGeomAlgo) is not translated — \
         BRepFill_NSections::totalsurf (BRepFill_NSections.cxx L167)"
    )
}

// ---------------------------------------------------------------------------
// Statics
// ---------------------------------------------------------------------------

/// OCCT static EdgeToBSpline (cxx L65-130) — get curve from edge and convert
/// it to bspline parameterized from 0 to 1 (duplicated in
/// BRepOffsetAPI_ThruSections.cxx).
fn edge_to_bspline(brep: &BRep, the_edge: &Shape) -> Curve3 {
    let mut a_bscurve: Option<BSplineCurve3>;
    if brep.edge(the_edge.clone()).degenerated {
        // degenerated edge : construction of a point curve
        // TopExp::Vertices(theEdge, vl, vf) — (V1 = vl <- first, V2 = vf <- last)
        // aPoles(1) = BRep_Tool::Pnt(vf); aPoles(2) = BRep_Tool::Pnt(vl);
        let (vl, vf) = top_exp_vertices(brep, the_edge);
        let p1 = brep.vertex(vf).point;
        let p2 = brep.vertex(vl).point;
        a_bscurve = Some(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p1, p2],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        });
    } else {
        // get the curve of the edge
        // BRep_Tool::Curve(theEdge, aLoc, aFirst, aLast) — the rcad
        // TEdgeData::curve carries no TopLoc (identity form).
        let (a_curve, a_first, a_last) = brep_tool_curve(brep, the_edge).expect("edge curve");

        // Geom_TrimmedCurve aTrimCurve(aCurve, aFirst, aLast);
        let a_trim_curve = Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(a_curve),
            first: a_first,
            last: a_last,
        });

        // GeomConvert_ApproxCurve anAppr(aCurveTemp, Confusion, GeomAbs_C1,
        //                                 16, 14);
        let an_appr =
            GeomConvertApproxCurve::new(&a_trim_curve, CONFUSION, GeomAbsShape::C1, 16, 14);
        if an_appr.has_result() {
            a_bscurve = an_appr.curve();
        } else {
            a_bscurve = None;
        }

        if a_bscurve.is_none() {
            // aBSCurve = GeomConvert::CurveToBSplineCurve(aTrimCurve);
            a_bscurve = Some(curve_to_bspline(&a_trim_curve, 0));
        }

        // apply transformation if needed (aLoc is the identity carrier)
        // if (!aLoc.IsIdentity()) { aBSCurve->Transform(aLoc.Transformation()); }

        // reparameterize to [0, 1] — BSplCLib::Reparametrize(0., 1., aKnots)
        let mut bs = a_bscurve.expect("bspline");
        let (mut knots, mults) = bs.knots_mults();
        reparametrize(0.0, 1.0, &mut knots);
        bs.set_knots(&knots, &mults);
        a_bscurve = Some(bs);
    }

    // reverse curve if edge is reversed
    let mut bs = a_bscurve.expect("bspline");
    if the_edge.orientation == Orientation::Reversed {
        // aBSCurve->Reverse();
        bs = bs.reversed();
    }
    Curve3::BSpline(bs)
}

/// OCCT static totalsurf (cxx L134-359) — the v-direction BSpline surface of
/// the section grid.
#[allow(clippy::too_many_arguments)]
fn totalsurf(
    brep: &BRep,
    shapes: &[Vec<Shape>], // [edge][section] (the OCCT Array2 (1, NbEdges, 1, NbSects))
    nb_sects: usize,
    nb_edges: usize,
    params: &[f64],
    w1_point: bool,
    w2_point: bool,
    u_closed: bool,
    v_closed: bool,
    my_pres3d: f64,
) -> BSplineSurface {
    let mut jdeb = 1usize;
    let mut jfin = nb_sects;

    // GeomFill_SectionGenerator section;
    let mut section = SectionGenerator::new();
    let mut bs1: Option<Curve3> = None;

    if w1_point {
        jdeb += 1;
        let edge = shapes[0][0].clone();
        // TopExp::Vertices(edge, vl, vf); Extremities(1) = Pnt(vf);
        // Extremities(2) = Pnt(vl);
        let (vl, vf) = top_exp_vertices(brep, &edge);
        let p1 = brep.vertex(vf).point;
        let p2 = brep.vertex(vl).point;
        // Geom_BSplineCurve BSPoint(Extremities, Bounds, Mult, 1)
        let bs_point = Curve3::BSpline(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p1, p2],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        });
        // section.AddCurve(BSPoint);
        gap_section_generator_add_curve(&mut section, &bs_point);
    }

    if w2_point {
        jfin -= 1;
    }

    for j in jdeb..=jfin {
        // case of looping sections
        if j == jfin && v_closed {
            // section.AddCurve(BS1);
            gap_section_generator_add_curve(&mut section, bs1.as_ref().expect("BS1"));
        } else {
            // read the first edge to initialise CompBS;
            let a_prev_edge = shapes[0][j - 1].clone();
            let mut curv_bs = edge_to_bspline(brep, &a_prev_edge);

            // initialization: GeomConvert_CompCurveToBSplineCurve CompBS(curvBS);
            let mut comp_bs =
                crate::brep_fill::brep_fill_section_law::CompCurveToBSpline::new(&curv_bs);

            for i in 2..=nb_edges {
                // read the edge
                let a_next_edge = shapes[i - 1][j - 1].clone();
                curv_bs = edge_to_bspline(brep, &a_next_edge);

                // concatenation
                // Bof = TopExp::CommonVertex(aPrevEdge, aNextEdge, ComV);
                let com_v = crate::fillet::chfi3d_builder_0::topexp_common_vertex(
                    &a_prev_edge,
                    &a_next_edge,
                );
                let eps_v = match &com_v {
                    Some(v) => {
                        // epsV = BRep_Tool::Tolerance(ComV);
                        let v_tol = v.clone();
                        let _ = v_tol;
                        brep_vertex_tolerance(brep, com_v.as_ref().expect("checked"))
                    }
                    None => CONFUSION,
                };
                // Bof = CompBS.Add(curvBS, epsV, true, false, 1);
                let mut bof = comp_bs.add(&curv_bs, eps_v, true, false, 1);
                if !bof {
                    // Bof = CompBS.Add(curvBS, 200 * epsV, true, false, 1);
                    bof = comp_bs.add(&curv_bs, 200.0 * eps_v, true, false, 1);
                }
                let _ = bof;

                // remember previous edge
                // aPrevEdge = aNextEdge;
            }

            // return the final section: BS = CompBS.BSplineCurve();
            // section.AddCurve(BS);
            let bs = comp_bs.bspline_curve();
            gap_section_generator_add_curve(&mut section, &bs);

            // case of looping sections
            if j == jdeb && v_closed {
                bs1 = Some(bs.clone());
            }
        }
    }
    let _ = bs1;

    if w2_point {
        let edge = shapes[nb_edges - 1][nb_sects - 1].clone();
        // TopExp::Vertices(edge, vl, vf); Extremities(1) = Pnt(vf);
        // Extremities(2) = Pnt(vl);
        let (vl, vf) = top_exp_vertices(brep, &edge);
        let p1 = brep.vertex(vf).point;
        let p2 = brep.vertex(vl).point;
        let bs_point = Curve3::BSpline(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p1, p2],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        });
        // section.AddCurve(BSPoint);
        gap_section_generator_add_curve(&mut section, &bs_point);
    }

    // section.SetParam(HPar); section.Perform(PConfusion());
    section.set_param(params.to_vec());
    // (Perform is part of the SectionGenerator producer GAP — unreachable
    // after the AddCurve carrier above; the literal call is kept in form)
    let line = Line::with_points(nb_sects as i32);
    let nb_it = 0;
    let degmin = 2;
    let degmax = 6;
    let known_p = true;
    let mut an_approx = GeomFillAppSurf::new(degmin, degmax, my_pres3d, my_pres3d, nb_it, known_p);
    let sp_approx = true;
    an_approx.perform(&line, &mut section, sp_approx);
    let uperiodic = u_closed;
    let vperiodic = v_closed;
    let mut nup = an_approx.surf_poles().len();
    let mut nvp = an_approx.surf_poles().first().map(|r| r.len()).unwrap_or(0);
    let mut umults = an_approx.surf_u_mults();
    let mut vmults = an_approx.surf_v_mults();

    if uperiodic {
        let nbuk = an_approx.surf_u_knots().len();
        umults[0] -= 1;
        umults[nbuk - 1] -= 1;
        nup -= 1;
    }

    if vperiodic {
        let nbvk = an_approx.surf_v_knots().len();
        vmults[0] -= 1;
        vmults[nbvk - 1] -= 1;
        nvp -= 1;
    }

    // poles / weights read-back (cxx L288-297)
    let poles = an_approx.surf_poles();
    let weights = an_approx.surf_weights();

    // To create non-rational surface if possible (cxx L300-346)
    let tol_eps = 1.0e-13;
    let mut v_rational = false;
    let mut u_rational = false;
    for j in 0..nvp {
        if !v_rational {
            for i in 0..nup - 1 {
                let delta = (weights[i][j] - weights[i + 1][j]).abs();
                if delta > tol_eps {
                    v_rational = true;
                    break;
                }
            }
        }
    }
    for i in 0..nup {
        if !u_rational {
            for j in 0..nvp - 1 {
                let delta = (weights[i][j] - weights[i][j + 1]).abs();
                if delta > tol_eps {
                    u_rational = true;
                    break;
                }
            }
        }
    }
    let mut weights = weights;
    if !v_rational && !u_rational {
        let the_weight = weights[0][0];
        for row in weights.iter_mut() {
            for w in row.iter_mut() {
                *w = the_weight;
            }
        }
    }

    // surface = new Geom_BSplineSurface(poles, weights, UKnots, VKnots,
    //                                    Umults, Vmults, UDegree, VDegree,
    //                                    uperiodic, vperiodic);
    BSplineSurface {
        degree_u: an_approx.u_degree() as usize,
        degree_v: an_approx.v_degree() as usize,
        knots_u: expand_knots(&an_approx.surf_u_knots(), &umults),
        knots_v: expand_knots(&an_approx.surf_v_knots(), &vmults),
        control_points: poles,
        weights,
    }
    .with_periodic_flags(uperiodic, vperiodic, &mut nup, &mut nvp)
}

/// The rcad BSplineSurface carries no periodic flags (see the geomfill
/// nsections.rs header note); the OCCT ctor arguments are consumed here.
trait PeriodicFlags {
    fn with_periodic_flags(self, uperiodic: bool, vperiodic: bool, nup: &mut usize, nvp: &mut usize) -> BSplineSurface;
}
impl PeriodicFlags for BSplineSurface {
    fn with_periodic_flags(self, uperiodic: bool, vperiodic: bool, nup: &mut usize, nvp: &mut usize) -> BSplineSurface {
        let _ = (uperiodic, vperiodic, nup, nvp);
        self
    }
}

// ---------------------------------------------------------------------------
// BRepFill_NSections
// ---------------------------------------------------------------------------

/// OCCT BRepFill_NSections (hxx L33-96).
pub struct BRepFillNSections {
    /// OCCT BRepFill_SectionLaw base sub-object.
    pub base: BRepFillSectionLawBase,
    /// OCCT NCollection_Sequence<TopoDS_Shape> myShapes (private).
    pub my_shapes: Vec<Shape>,
    /// OCCT NCollection_Sequence<gp_Trsf> myTrsfs (private).
    pub my_trsfs: Vec<Trsf>,
    /// OCCT NCollection_Sequence<double> myParams (private).
    pub my_params: Vec<f64>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myEdges (private).
    pub my_edges: Vec<Vec<Shape>>,
    /// OCCT handle(Geom_BSplineSurface) mySurface (private).
    pub my_surface: Option<BSplineSurface>,
    /// OCCT double VFirst.
    pub vfirst: f64,
    /// OCCT double VLast.
    pub vlast: f64,
}

impl BRepFillNSections {
    /// OCCT BRepFill_NSections(S, Build) (cxx L363-379).
    pub fn new(brep: &BRep, s: Vec<Shape>, build: bool) -> Self {
        let my_shapes = s.clone();
        let vfirst = 0.0;
        let vlast = 1.0;
        // for (int i = 1; i <= S.Length(); i++) par.Append(i - 1);
        let par: Vec<f64> = (0..s.len()).map(|i| i as f64).collect();
        let mut this = BRepFillNSections {
            base: empty_base(),
            my_shapes,
            my_trsfs: Vec::new(),
            my_params: par.clone(),
            my_edges: Vec::new(),
            my_surface: None,
            vfirst,
            vlast,
        };
        this.init(brep, &par, build);
        this.base.my_done = true;
        this
    }

    /// OCCT BRepFill_NSections(S, Transformations, P, VF, VL, Build)
    /// (cxx L386-418) — WSeq + Param.
    pub fn new_with_params(
        brep: &BRep,
        s: Vec<Shape>,
        transformations: Vec<Trsf>,
        p: Vec<f64>,
        vf: f64,
        vl: f64,
        build: bool,
    ) -> Self {
        // ok = P is strictly increasing
        let mut ok = true;
        for iseq in 0..p.len().saturating_sub(1) {
            ok = ok && (p[iseq] < p[iseq + 1]);
        }
        if ok {
            let mut this = BRepFillNSections {
                base: empty_base(),
                my_shapes: s,
                my_trsfs: transformations,
                my_params: p.clone(),
                my_edges: Vec::new(),
                my_surface: None,
                vfirst: vf,
                vlast: vl,
            };
            this.init(brep, &p, build);
            this.base.my_done = true;
            this
        } else {
            // myDone = false — the law carries no usable state.
            BRepFillNSections {
                base: BRepFillSectionLawBase {
                    my_done: false,
                    ..empty_base()
                },
                my_shapes: Vec::new(),
                my_trsfs: Vec::new(),
                my_params: Vec::new(),
                my_edges: Vec::new(),
                my_surface: None,
                vfirst: vf,
                vlast: vl,
            }
        }
    }

    /// OCCT Init (cxx L424-651) — Create a table of GeomFill_SectionLaw.
    pub fn init(&mut self, brep: &BRep, p: &[f64], build: bool) {
        let nb_sects = p.len();
        let mut ideb = 1usize;
        let mut ifin = nb_sects;
        let mut w1_point = true;
        let mut w2_point = true;
        let mut first = 0.0;
        let mut last = 0.0;

        // Check if the start and end wires are punctual (L437-456)
        let w = self.my_shapes[0].clone();
        for e in wire_edges(brep, &w) {
            w1_point = w1_point && brep.edge(e.clone()).degenerated;
        }
        if w1_point {
            ideb += 1;
        }
        let w = self.my_shapes[nb_sects - 1].clone();
        for e in wire_edges(brep, &w) {
            w2_point = w2_point && brep.edge(e.clone()).degenerated;
        }
        if w2_point {
            ifin -= 1;
        }

        // Check if the start and end wires are identical (L459)
        let vclosed = self.my_shapes[0].is_same(&self.my_shapes[nb_sects - 1]);
        self.base.vclosed = vclosed;

        // Count the number of non-degenerated edges (L462-470)
        let w = self.my_shapes[ideb - 1].clone();
        let mut nb_edge = 0usize;
        for e in wire_edges(brep, &w) {
            if !brep.edge(e.clone()).degenerated {
                nb_edge += 1;
            }
        }

        // myEdges = new HArray2(1, NbEdge, 1, NbSects)
        self.my_edges = vec![vec![Shape::null(); nb_sects]; nb_edge];

        // Fill tables (L475-546)
        let mut uclosed = true;
        for jj in ideb..=ifin {
            let w = self.my_shapes[jj - 1].clone();

            let mut ii = 0usize;
            for e in wire_edges(brep, &w) {
                if ii >= nb_edge {
                    break;
                }
                if !brep.edge(e.clone()).degenerated {
                    ii += 1;
                    self.my_edges[ii - 1][jj - 1] = e.clone();
                    if e.orientation == Orientation::Forward {
                        self.base.my_indices.insert(e.ptr_id(), ii as i32);
                    } else {
                        self.base.my_indices.insert(e.ptr_id(), -(ii as i32));
                    }
                }
            }

            // Is the law closed by U ? (L500-545)
            let mut w_closed = brep.has_flag(w.clone(), rcad_kernel::topo::topods::tshape_flags::CLOSED);
            if !w_closed {
                // if unsure about the flag, make check
                let edge1 = self.my_edges[nb_edge - 1][jj - 1].clone();
                let edge2 = self.my_edges[0][jj - 1].clone();

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
                    w_closed = true;
                } else {
                    let u1 = crate::brep_fill::compatible_wires::brep_tool_parameter(brep, &v1, &edge1);
                    let u2 = crate::brep_fill::compatible_wires::brep_tool_parameter(brep, &v2, &edge2);
                    let eps =
                        crate::brep_fill::compatible_wires::shape_tolerance(brep, &v2)
                            + crate::brep_fill::compatible_wires::shape_tolerance(brep, &v1);
                    let p1 = brep.edge(edge1.clone()).curve.clone().map(|c| c.point_at(u1));
                    let p2 = brep.edge(edge2.clone()).curve.clone().map(|c| c.point_at(u2));
                    w_closed = match (p1, p2) {
                        (Some(p1), Some(p2)) => p1.distance(p2) <= eps,
                        _ => false,
                    };
                }
            }
            if !w_closed {
                uclosed = false;
            }
        }
        self.base.uclosed = uclosed;

        // point sections at end (L548-569)
        if w1_point {
            let w = self.my_shapes[0].clone();
            let e = wire_edges(brep, &w).first().cloned().unwrap_or_else(|| Shape::null());
            for ii in 0..nb_edge {
                self.my_edges[ii][0] = e.clone();
            }
        }

        if w2_point {
            let w = self.my_shapes[nb_sects - 1].clone();
            let e = wire_edges(brep, &w).first().cloned().unwrap_or_else(|| Shape::null());
            for ii in 0..nb_edge {
                self.my_edges[ii][nb_sects - 1] = e.clone();
            }
        }

        // myLaws = new HArray1(1, NbEdge);
        self.base.my_laws = Vec::new();

        // constexpr double tol = Precision::Confusion();
        let tol = CONFUSION;
        // mySurface = totalsurf(myEdges->Array2(), myShapes.Length(), NbEdge,
        //                       myParams, w1Point, w2Point, uclosed, vclosed, tol);
        let surface = totalsurf(
            brep,
            &self.my_edges,
            self.my_shapes.len(),
            nb_edge,
            &self.my_params,
            w1_point,
            w2_point,
            uclosed,
            vclosed,
            tol,
        );
        let mut surface = surface;

        // Increase the degree so that the position D2 on GeomFill_NSections
        // could be correct (L585-590)
        if surface.degree_v < 2 {
            // mySurface->IncreaseDegree(mySurface->UDegree(), 2) — the V
            // direction degree elevation (the geomfill/nsections.rs
            // increase_degree_v precedent; here to degree 2 in V).
            increase_degree_v(&mut surface, 2);
        }
        self.my_surface = Some(surface);

        // Fill tables (L593-650)
        if build {
            let mut laws: Vec<Rc<RefCell<dyn SectionLaw>>> = Vec::new();
            for ii in 1..=nb_edge {
                let mut nc: Vec<Curve3> = Vec::new();
                for jj in 1..=nb_sects {
                    let e = self.my_edges[ii - 1][jj - 1].clone();
                    let c: Curve3;
                    if brep.edge(e.clone()).degenerated {
                        // point curve
                        let (vl, vf) = top_exp_vertices(brep, &e);
                        let p1 = brep.vertex(vf.clone()).point;
                        let p2 = brep.vertex(vl).point;
                        c = Curve3::BSpline(BSplineCurve3 {
                            degree: 1,
                            knots: vec![0.0, 0.0, 1.0, 1.0],
                            control_points: vec![p1, p2],
                            weights: vec![1.0, 1.0],
                            is_periodic: false,
                        });
                    } else {
                        let (mut c0, mut first0, mut last0) =
                            brep_tool_curve(brep, &e).expect("edge curve");
                        if e.orientation == Orientation::Reversed {
                            // CBis = C->Reversed();
                            let cbis = crate::brep_fill::brep_fill_section_law::reversed_curve(&c0);
                            let aux = crate::brep_fill::brep_fill_section_law::reversed_parameter_of(&c0, first0);
                            first0 = crate::brep_fill::brep_fill_section_law::reversed_parameter_of(&c0, last0);
                            last0 = aux;
                            c0 = cbis;
                        }
                        if (ii > 1) || !crate::brep_fill::brep_fill_section_law::brep_tool_is_closed_edge(brep, &e) {
                            // Cut C
                            c0 = Curve3::Trimmed(TrimmedCurve3 {
                                curve: Box::new(c0),
                                first: first0,
                                last: last0,
                            });
                        }
                        c = c0;
                    }
                    nc.push(c);
                }

                let ufirst = (ii - 1) as f64;
                let ulast = ii as f64;
                // myLaws->ChangeValue(ii) = new GeomFill_NSections(NC, myTrsfs,
                //     myParams, Ufirst, Ulast, VFirst, VLast, mySurface);
                let law = NSections::new_with_reference(
                    nc,
                    self.my_trsfs.clone(),
                    self.my_params.clone(),
                    ufirst,
                    ulast,
                    self.vfirst,
                    self.vlast,
                    self.my_surface.clone(),
                );
                laws.push(Rc::new(RefCell::new(law)));
            }
            self.base.my_laws = laws;
        }
    }

    /// OCCT IsVertex (cxx L655-658).
    pub fn is_vertex(&self) -> bool {
        false
    }

    /// OCCT IsConstant (cxx L662-665).
    pub fn is_constant_law(&self) -> bool {
        false
    }

    /// OCCT Vertex (cxx L669-694).
    pub fn vertex(&self, brep: &mut BRep, index: i32, param: f64) -> Shape {
        // BRep_Builder B; TopoDS_Vertex V; B.MakeVertex(V);
        let mut b = BRepBuilder::new();
        let mut v = b.add_vertex(brep, DVec3::ZERO, 0.0);
        let nb_col = self.my_edges[0].len() as i32;

        if index <= nb_col {
            // Curve = myLaws->Value(Index)->BSplineSurface()->VIso(Param);
            // double first = Curve->FirstParameter(); Curve->D0(first, P);
            // B.UpdateVertex(V, P, Precision::Confusion());
            let loi = self.base.my_laws[(index - 1) as usize].clone();
            let surf = loi.borrow().bspline_surface().cloned();
            if let Some(surf) = surf {
                let curve = Curve3::BSpline(bspline_surface_viso(&surf, param));
                let first = curve_first_parameter(&curve);
                let p = curve.point_at(first);
                b.update_vertex_point(brep, v.clone(), p, CONFUSION);
            }
        } else if index == nb_col + 1 {
            // Curve = myLaws->Value(Index - 1)->BSplineSurface()->VIso(Param);
            // double last = Curve->LastParameter(); Curve->D0(last, P);
            let loi = self.base.my_laws[(index - 2) as usize].clone();
            let surf = loi.borrow().bspline_surface().cloned();
            if let Some(surf) = surf {
                let curve = Curve3::BSpline(bspline_surface_viso(&surf, param));
                let last = curve_last_parameter(&curve);
                let p = curve.point_at(last);
                b.update_vertex_point(brep, v.clone(), p, CONFUSION);
            }
        }

        v
    }

    /// OCCT VertexTol (cxx L700-761) — the same BS-endpoint evaluation as
    /// BRepFill_ShapeLaw::VertexTol.
    pub fn vertex_tol(&self, index: i32, param: f64) -> f64 {
        let mut tol = CONFUSION;
        let nb_col = self.my_edges[0].len() as i32;
        let i1;
        let i2;
        if (index == 0) || (index == nb_col) {
            if !self.base.uclosed {
                return tol; // The least possible error
            }
            i1 = nb_col;
            i2 = 1;
        } else {
            i1 = index;
            i2 = i1 + 1;
        }

        let p_first = section_law_endpoint(&self.base, i1, param, true);
        let p_last = section_law_endpoint(&self.base, i2, param, false);
        tol += p_first.distance(p_last);
        tol
    }

    /// OCCT ConcatenedLaw (cxx L765-786).
    pub fn concatened_law(&self) -> Option<Rc<RefCell<dyn SectionLaw>>> {
        if self.base.my_laws.len() == 1 {
            return Some(self.base.my_laws[0].clone());
        }

        let surface = self.my_surface.as_ref()?;
        // mySurface->Bounds(Ufirst, Ulast, Vfirst, Vlast) — the rcad
        // BSplineSurface carries the domain in its knot vector.
        let ufirst = surface.knots_u.first().copied()?;
        let ulast = surface.knots_u.last().copied()?;
        let vfirst = surface.knots_v.first().copied()?;
        let vlast = surface.knots_v.last().copied()?;

        let mut n_compo: Vec<Curve3> = Vec::new();
        for jj in 1..=self.my_shapes.len() {
            // NCompo.Append(mySurface->VIso(myParams(jj)));
            n_compo.push(Curve3::BSpline(bspline_surface_viso(surface, self.my_params[jj - 1])));
        }
        let law = NSections::new_with_reference(
            n_compo,
            self.my_trsfs.clone(),
            self.my_params.clone(),
            ufirst,
            ulast,
            vfirst,
            vlast,
            Some(surface.clone()),
        );
        Some(Rc::new(RefCell::new(law)))
    }

    /// OCCT Continuity (cxx L790-860).
    pub fn continuity(&self, brep: &BRep, index: i32, tol_angular: f64) -> GeomAbsShape {
        let mut cont = GeomAbsShape::C0;
        let nb_col = self.my_edges[0].len() as i32;

        for jj in 1..=self.my_shapes.len() {
            let edge1;
            let edge2;
            if (index == 0) || (index == nb_col) {
                if !self.base.uclosed {
                    return GeomAbsShape::C0; // The least possible error
                }
                edge1 = self.my_edges[(nb_col - 1) as usize][jj - 1].clone();
                edge2 = self.my_edges[0][jj - 1].clone();
            } else {
                edge1 = self.my_edges[(index - 1) as usize][jj - 1].clone();
                edge2 = self.my_edges[index as usize][jj - 1].clone();
            }

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

            let cont_jj;
            if brep.edge(edge1.clone()).degenerated || brep.edge(edge2.clone()).degenerated {
                cont_jj = GeomAbsShape::CN;
            } else {
                // cont_jj = BRepLProp::Continuity(Curve1, Curve2, U1, U2, Eps,
                //                                 TolAngular);
                // GAP: BRepLProp::Continuity (TKTopAlgo) is not translated —
                // the ShapeLaw precedent; the section continuity test keeps
                // the OCCT failure path.
                panic!(
                    "GAP: BRepLProp::Continuity (TKTopAlgo) is not translated — \
                     BRepFill_NSections::Continuity (BRepFill_NSections.cxx L846)"
                );
            }

            if jj == 1 {
                cont = cont_jj;
            }
            if rank_of(cont) > rank_of(cont_jj) {
                cont = cont_jj;
            }
        }

        cont
    }

    /// OCCT D0 (cxx L864-881).
    pub fn d0(&mut self, brep: &mut BRep, v_param: f64, s: &mut Shape) {
        // BRepLib_MakeWire MW;
        let mut mw = BRepBuilder::new();
        // B.MakeWire(W) — the wire root.
        let wire = mw.make_wire(brep);
        let nb_edge = self.base.my_laws.len();
        for ii in 1..=nb_edge {
            // Curve = myLaws->Value(ii)->BSplineSurface()->VIso(V);
            let loi = self.base.my_laws[ii - 1].clone();
            let surf = loi.borrow().bspline_surface().map(|s| s.clone());
            let Some(surf) = surf else { continue };
            let curve = Curve3::BSpline(bspline_surface_viso(&surf, v_param));
            let first = curve_first_parameter(&curve);
            let last = curve_last_parameter(&curve);
            // TopoDS_Edge E = BRepLib_MakeEdge(Curve, first, last);
            let e = mw.add_edge(
                brep,
                Some(curve),
                Shape::null(),
                Shape::null(),
                [first, last],
            );
            // MW.Add(E);
            mw.add_to_wire(brep, wire.clone(), e);
        }
        // TopAbs_Orientation Orien = TopAbs_FORWARD; S = Wire.Oriented(Orien);
        let oriented = wire; // already FORWARD
        *s = oriented;
    }
}

/// The empty base state before Init.
fn empty_base() -> BRepFillSectionLawBase {
    BRepFillSectionLawBase {
        my_laws: Vec::new(),
        uclosed: false,
        vclosed: false,
        my_done: false,
        my_indices: std::collections::HashMap::new(),
        my_iterator: WireExplorerState {
            wire: Shape::null(),
            edges: Vec::new(),
            cursor: 0,
        },
    }
}

/// The shared Geom_BSplineCurve endpoint evaluation of VertexTol (the
/// ShapeLaw form; the cxx body is the same L727-759).
fn section_law_endpoint(
    base: &BRepFillSectionLawBase,
    index: i32,
    param: f64,
    first_knot: bool,
) -> DVec3 {
    let loi = base.my_laws[(index - 1) as usize].clone();
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

    let bs = BSplineCurve3 {
        degree,
        knots: expand_knots(&knots, &mults),
        control_points: poles,
        weights: weight,
        is_periodic: loi.borrow().is_u_periodic(),
    };
    let bs_curve = Curve3::BSpline(bs);
    if first_knot {
        bs_curve.point_at(knots[knots.len() - 1])
    } else {
        bs_curve.point_at(knots[0])
    }
}

/// The GeomAbs rank of the continuity (the `cont > cont_jj` comparison).
fn rank_of(c: GeomAbsShape) -> i32 {
    match c {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::G1 => 1,
        GeomAbsShape::C1 => 2,
        GeomAbsShape::G2 => 3,
        GeomAbsShape::C2 => 4,
        GeomAbsShape::C3 => 5,
        GeomAbsShape::CN => 6,
    }
}

/// The last-vertex tolerance read (BRep_Tool::Tolerance(ComV)).
fn brep_vertex_tolerance(brep: &BRep, v: &Shape) -> f64 {
    brep.vertex(v.clone()).tolerance
}

/// OCCT Geom_BSplineSurface::VIso(V) — the u-varying iso curve at V (the
/// homogeneous De Boor evaluation of the opposite-direction basis, the
/// geomfill/nsections.rs precedent).
fn bspline_surface_viso(surf: &BSplineSurface, v: f64) -> BSplineCurve3 {
    use rcad_kernel::math::bspl::de_boor_homo;
    let rational = surf
        .weights
        .iter()
        .any(|row| row.iter().any(|&w| w != 1.0));
    let nb_u = surf.control_points.len();
    let nb_v = surf.control_points.first().map(|r| r.len()).unwrap_or(0);
    let mut poles = Vec::with_capacity(nb_u);
    let mut weights = Vec::with_capacity(nb_u);
    for ii in 0..nb_u {
        let column: Vec<DVec3> = (0..nb_v).map(|jj| surf.control_points[ii][jj]).collect();
        let column_w: Vec<f64> = (0..nb_v)
            .map(|jj| if rational { surf.weights[ii][jj] } else { 1.0 })
            .collect();
        let h = de_boor_homo(surf.degree_v, &surf.knots_v, &column, &column_w, v);
        weights.push(h[3]);
        poles.push(DVec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]));
    }
    BSplineCurve3 {
        degree: surf.degree_u,
        knots: surf.knots_u.clone(),
        control_points: poles,
        weights,
        is_periodic: false,
    }
}

/// OCCT Geom_BSplineSurface::IncreaseDegree in the V direction — degree
/// elevation per U row through the 1D kernel (the geomfill/nsections.rs
/// re-host form).
fn increase_degree_v(surf: &mut BSplineSurface, new_degree: usize) {
    use rcad_kernel::math::bspl_lib::increase_degree as bspl_increase_degree;
    let nb_u = surf.control_points.len();
    let nb_v = surf.control_points.first().map(|r| r.len()).unwrap_or(0);
    let rational = surf
        .weights
        .iter()
        .any(|col| col.iter().any(|&w| w != 1.0));
    let (vk, vm) = knots_mults_of(&surf.knots_v);
    let dim = if rational { 4 } else { 3 };
    let new_nb_v = {
        // increase_degree_count_knots keeps the knot count; the pole count
        // grows by (new_degree - degree) per span-boundary effect — the
        // batch-1 form sizes the row from the elevated mults sum.
        vm.iter().map(|&m| m as usize).sum::<usize>() + (new_degree - surf.degree_v)
    };
    let mut new_poles = vec![vec![DVec3::ZERO; new_nb_v]; nb_u];
    let mut new_weights = vec![vec![1.0f64; new_nb_v]; nb_u];
    let mut new_knots = vec![0.0f64; vk.len() + (new_degree - surf.degree_v)];
    let mut new_mults = vm.clone();
    new_mults[0] += (new_degree - surf.degree_v) as i32;
    *new_mults.last_mut().expect("nonempty") += (new_degree - surf.degree_v) as i32;
    for ii in 0..nb_u {
        let column: Vec<DVec3> = (0..nb_v).map(|jj| surf.control_points[ii][jj]).collect();
        let column_w: Vec<f64> = (0..nb_v)
            .map(|jj| if rational { surf.weights[ii][jj] } else { 1.0 })
            .collect();
        let flat_len = nb_v * dim;
        let mut flat = vec![0.0f64; flat_len];
        for jj in 0..nb_v {
            if dim == 4 {
                flat[jj * 4] = column[jj].x * column_w[jj];
                flat[jj * 4 + 1] = column[jj].y * column_w[jj];
                flat[jj * 4 + 2] = column[jj].z * column_w[jj];
                flat[jj * 4 + 3] = column_w[jj];
            } else {
                flat[jj * 3] = column[jj].x;
                flat[jj * 3 + 1] = column[jj].y;
                flat[jj * 3 + 2] = column[jj].z;
            }
        }
        let mut new_flat = vec![0.0f64; new_nb_v * dim];
        let mut new_knots_row = vec![0.0f64; new_knots.len()];
        let mut new_mults_row = new_mults.clone();
        bspl_increase_degree(
            surf.degree_v,
            new_degree,
            false,
            dim,
            &flat,
            &vk,
            &vm,
            &mut new_flat,
            &mut new_knots_row,
            &mut new_mults_row,
        );
        for jj in 0..new_nb_v {
            if dim == 4 {
                let w = new_flat[jj * 4 + 3];
                new_poles[ii][jj] = DVec3::new(
                    new_flat[jj * 4] / w,
                    new_flat[jj * 4 + 1] / w,
                    new_flat[jj * 4 + 2] / w,
                );
                new_weights[ii][jj] = w;
            } else {
                new_poles[ii][jj] = DVec3::new(
                    new_flat[jj * 3],
                    new_flat[jj * 3 + 1],
                    new_flat[jj * 3 + 2],
                );
            }
        }
    }
    surf.degree_v = new_degree;
    surf.control_points = new_poles;
    surf.weights = new_weights;
    // rebuild the flat V knot vector from the elevated (knots, mults)
    surf.knots_v = {
        let mut flat = Vec::with_capacity(new_mults.iter().map(|&m| m as usize).sum());
        for (k, m) in new_knots.iter().zip(new_mults.iter()) {
            for _ in 0..*m {
                flat.push(*k);
            }
        }
        flat
    };
}

/// The distinct knots / multiplicities of a flat knot vector.
fn knots_mults_of(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i == 0 || *k != knots[knots.len() - 1] {
            knots.push(*k);
            mults.push(1);
        } else {
            *mults.last_mut().expect("nonempty") += 1;
        }
    }
    (knots, mults)
}

/// OCCT BRepFill_SectionLaw inheritance: BRepFill_NSections overrides.
impl BRepFillSectionLawOps for BRepFillNSections {
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
    fn concatened_law(&self, _brep: &BRep) -> Option<Rc<RefCell<dyn SectionLaw>>> {
        self.concatened_law()
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
    fn d0(&mut self, brep: &mut BRep, u: f64, s: &mut Shape) {
        self.d0(brep, u, s)
    }
}
