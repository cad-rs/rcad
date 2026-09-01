//! OCCT HelixBRep package (TKHelix) — helix wire building.
//!
//! 1:1 translation of `HelixBRep_BuilderHelix.hxx` + `.cxx` (L46-658).
//!
//! BRep-layer correspondence notes (Rust has no TopoDS handles):
//! - `BRep_Builder::MakeVertex(aV, P, tol)` -> [`brep.add_tvertex`].  OCCT
//!   always creates a new vertex TShape; the coincident-junction merge then
//!   happens inside `BRepBuilderAPI_MakeWire::Add`.  rcad merges by the
//!   position-quantized vertex identity of `add_tvertex`, which produces the
//!   same wire topology (OCCT helix/standard C1: V=11, E=10 — V = E+1).
//! - `BRepBuilderAPI_MakeEdge::Init(C, V1, V2)` -> the edge TShape carries
//!   the BSpline with range [FirstParameter, LastParameter]; the vertices sit
//!   exactly at the curve ends (their projection is exact).
//! - The per-part `BRep_Builder::MakeWire` + `Add` and the final
//!   `BRepBuilderAPI_MakeWire` merge are collapsed into edge lists
//!   (`Vec<Vec<Shape>>`) finalized by one `add_twire` — observationally
//!   identical: the result is a single wire over the same edges/vertices.
//! - `BRep_Builder::UpdateEdge(E, C, tol)` (curve replacement in
//!   `SmoothingEdges`) -> the edge `Shape` handle is replaced by a new edge
//!   TShape carrying the smoothed curve (rcad TShapes are immutable Arcs).

use std::sync::Arc;

use glam::DVec3;
use rcad_kernel::geom::{BSplineCurve3, Curve3};
use rcad_kernel::math::gp::{Ax1, Ax3, Lin};
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape, TEdgeData, TShape,
                                TVertexData};

/// OCCT HelixBRep_BuilderHelix.
pub struct BuilderHelix {
    /// The BRep pool owning every TShape built by this builder.
    pub brep: BRep,
    my_axis3: Ax3,
    my_diams: Option<Vec<f64>>,
    my_heights: Option<Vec<f64>>,
    my_pitches: Option<Vec<f64>>,
    my_is_pitches: Option<Vec<bool>>,
    my_tolerance: f64,
    my_tol_reached: f64,
    my_continuity: GeomAbsShape,
    my_max_degree: i32,
    my_max_segments: i32,
    my_error_status: i32,
    my_warning_status: i32,
    my_shape: Option<Shape>,
    my_n_parts: usize,
}

impl Default for BuilderHelix {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT NCollection_Array1 1-based Value(i) helper.
fn val(arr: &[f64], i: usize) -> f64 {
    arr[i - 1]
}

impl BuilderHelix {
    /// OCCT HelixBRep_BuilderHelix() (L46-70).
    pub fn new() -> Self {
        let mut axis3 = Ax3::new();
        axis3.set_direction(DVec3::Z);
        axis3.set_location(DVec3::ZERO);
        BuilderHelix {
            brep: BRep::new(),
            my_axis3: axis3,
            my_diams: None,
            my_heights: None,
            my_pitches: None,
            my_is_pitches: None,
            my_n_parts: 1,
            my_shape: None,
            my_tolerance: 0.0001,
            my_continuity: GeomAbsShape::C1,
            my_max_degree: 8,
            my_max_segments: 1000,
            my_tol_reached: 99.0,
            my_error_status: 1,
            my_warning_status: 1,
        }
    }

    /// OCCT SetParameters(theAxis, theDiams, theHeights, thePitches,
    /// bIsPitches) — general composite helix (L78-115).
    pub fn set_parameters(
        &mut self,
        the_axis: &Ax3,
        the_diams: &[f64],
        the_heights: &[f64],
        the_pitches: &[f64],
        the_is_pitches: &[bool],
    ) {
        self.my_n_parts = the_diams.len() - 1;

        self.my_axis3 = *the_axis;

        self.my_diams = None;
        self.my_heights = None;
        self.my_pitches = None;
        self.my_shape = None;

        self.my_error_status = 1;
        self.my_warning_status = 1;

        if self.my_n_parts != the_heights.len()
            || self.my_n_parts != the_pitches.len()
            || self.my_n_parts != the_is_pitches.len()
        {
            panic!("HelixBRep_BuilderHelix::SetParameters: wrong array dimension");
        }

        self.my_diams = Some(the_diams.to_vec());
        self.my_heights = Some(the_heights.to_vec());
        self.my_pitches = Some(the_pitches.to_vec());
        self.my_is_pitches = Some(the_is_pitches.to_vec());

        self.my_error_status = 0;
        self.my_warning_status = 0;
    }

    /// OCCT SetParameters(theAxis, theDiam, theHeights, thePitches,
    /// bIsPitches) — pure helix (L119-129).
    pub fn set_parameters_helix(
        &mut self,
        the_axis: &Ax3,
        the_diam: f64,
        the_heights: &[f64],
        the_pitches: &[f64],
        the_is_pitches: &[bool],
    ) {
        let a_nb_parts = the_heights.len();
        let mut a_diams = vec![0.0f64; a_nb_parts + 1];
        for d in a_diams.iter_mut() {
            *d = the_diam;
        }
        self.set_parameters(the_axis, &a_diams, the_heights, the_pitches, the_is_pitches);
    }

    /// OCCT SetParameters(theAxis, theDiam1, theDiam2, theHeights, thePitches,
    /// bIsPitches) — pure spiral (L133-161).
    pub fn set_parameters_spiral(
        &mut self,
        the_axis: &Ax3,
        the_diam1: f64,
        the_diam2: f64,
        the_heights: &[f64],
        the_pitches: &[f64],
        the_is_pitches: &[bool],
    ) {
        let a_nb_parts = the_heights.len();
        let mut a_diams = vec![0.0f64; a_nb_parts + 1];

        let mut an_h = 0.0f64;
        for i in 1..=the_heights.len() {
            an_h += val(the_heights, i);
        }
        let k = (the_diam2 - the_diam1) / an_h;
        a_diams[0] = the_diam1;
        a_diams[a_nb_parts] = the_diam2;

        an_h = val(the_heights, 1);
        let mut j = 2usize;
        for i in 2..=the_heights.len() {
            a_diams[j - 1] = the_diam1 + k * an_h;
            an_h += val(the_heights, i);
            j += 1;
        }

        self.set_parameters(the_axis, &a_diams, the_heights, the_pitches, the_is_pitches);
    }

    /// OCCT SetApproxParameters (L165-172).
    pub fn set_approx_parameters(&mut self, a_tolerance: f64, a_max_degree: i32, a_cont: GeomAbsShape) {
        self.my_tolerance = a_tolerance;
        self.my_max_degree = a_max_degree;
        self.my_continuity = a_cont;
    }

    /// OCCT ToleranceReached (L176-179).
    pub fn tolerance_reached(&self) -> f64 {
        self.my_tol_reached
    }

    /// OCCT ErrorStatus.
    pub fn error_status(&self) -> i32 {
        self.my_error_status
    }

    /// OCCT WarningStatus.
    pub fn warning_status(&self) -> i32 {
        self.my_warning_status
    }

    /// OCCT Shape() — the resulting wire (single shape in the builder BRep).
    pub fn shape(&self) -> Option<&Shape> {
        self.my_shape.as_ref()
    }

    /// The BRep pool owning all built TShapes (for tests / STEP export).
    pub fn into_brep(self) -> BRep {
        self.brep
    }

    /// OCCT SetParameters(theAxis, theDiams, thePitches, theNbTurns) —
    /// general composite helix by turns (L576-601).
    pub fn set_parameters_turns(
        &mut self,
        the_axis: &Ax3,
        the_diams: &[f64],
        the_pitches: &[f64],
        the_nb_turns: &[f64],
    ) {
        let a_nb_parts = the_diams.len() - 1;

        if a_nb_parts != the_pitches.len() || a_nb_parts != the_nb_turns.len() {
            panic!("HelixBRep_BuilderHelix::SetParameters: wrong array dimension");
        }

        let mut a_heights = vec![0.0f64; a_nb_parts];
        let b_is_pitches = vec![true; a_nb_parts];
        for i in 1..=the_pitches.len() {
            a_heights[i - 1] = val(the_pitches, i) * val(the_nb_turns, i);
        }

        self.set_parameters(the_axis, the_diams, &a_heights, the_pitches, &b_is_pitches);
    }

    /// OCCT SetParameters(theAxis, theDiam, thePitches, theNbTurns) — pure
    /// helix by turns (L605-629).
    pub fn set_parameters_helix_turns(
        &mut self,
        the_axis: &Ax3,
        the_diam: f64,
        the_pitches: &[f64],
        the_nb_turns: &[f64],
    ) {
        let a_nb_parts = the_pitches.len();

        if a_nb_parts != the_nb_turns.len() {
            panic!("HelixBRep_BuilderHelix::SetParameters: wrong array dimension");
        }

        let mut a_heights = vec![0.0f64; a_nb_parts];
        let b_is_pitches = vec![true; a_nb_parts];
        for i in 1..=the_pitches.len() {
            a_heights[i - 1] = val(the_pitches, i) * val(the_nb_turns, i);
        }

        self.set_parameters_helix(the_axis, the_diam, &a_heights, the_pitches, &b_is_pitches);
    }

    /// OCCT SetParameters(theAxis, theDiam1, theDiam2, thePitches,
    /// theNbTurns) — pure spiral by turns (L633-658).
    pub fn set_parameters_spiral_turns(
        &mut self,
        the_axis: &Ax3,
        the_diam1: f64,
        the_diam2: f64,
        the_pitches: &[f64],
        the_nb_turns: &[f64],
    ) {
        let a_nb_parts = the_pitches.len();

        if a_nb_parts != the_nb_turns.len() {
            panic!("HelixBRep_BuilderHelix::SetParameters: wrong array dimension");
        }

        let mut a_heights = vec![0.0f64; a_nb_parts];
        let b_is_pitches = vec![true; a_nb_parts];
        for i in 1..=the_pitches.len() {
            a_heights[i - 1] = val(the_pitches, i) * val(the_nb_turns, i);
        }

        self.set_parameters_spiral(
            the_axis,
            the_diam1,
            the_diam2,
            &a_heights,
            the_pitches,
            &b_is_pitches,
        );
    }

    /// OCCT HelixBRep_BuilderHelix::Perform (L183-257).
    pub fn perform(&mut self) {
        if self.my_error_status != 0 {
            return;
        }

        self.my_tol_reached = 0.0;
        self.my_shape = None;

        let mut an_lst: Vec<Vec<Shape>> = Vec::new();

        let mut an_axis = self.my_axis3.axis.clone();
        let mut a_p_start = self.my_axis3.location();
        a_p_start += 0.5 * val(self.my_diams.as_ref().unwrap(), 1) * self.my_axis3.x_direction;
        let b_is_clockwise = self.my_axis3.direct();

        for i in 1..=self.my_n_parts {
            let a_height = val(self.my_heights.as_ref().unwrap(), i);

            let a_pitch = if self.my_is_pitches.as_ref().unwrap()[i - 1] {
                val(self.my_pitches.as_ref().unwrap(), i)
            } else {
                a_height / val(self.my_pitches.as_ref().unwrap(), i)
            };

            let a_taper_angle = (0.5
                * (val(self.my_diams.as_ref().unwrap(), i + 1)
                    - val(self.my_diams.as_ref().unwrap(), i))
                / a_height)
                .atan();

            // BuildPart(anAxis, aPStart, aHeight, aPitch, aTaperAngle,
            //           bIsClockwise, aPart);
            let (a_part, v_first, v_last) = self.build_part(
                &an_axis,
                a_p_start,
                a_height,
                a_pitch,
                a_taper_angle,
                b_is_clockwise,
            );
            if self.my_error_status != 0 {
                return;
            }

            a_p_start = vertex_point(&v_last);

            let dir = an_axis.direction;
            an_axis.location += a_height * dir;

            an_lst.push(a_part);
        }

        self.smoothing(&mut an_lst);

        // BRepBuilderAPI_MakeWire merge of all part wires.
        let mut all_edges = Vec::new();
        for part in &an_lst {
            all_edges.extend(part.iter().cloned());
        }
        let wire = self.brep.add_twire(all_edges);

        self.my_shape = Some(wire);
    }

    /// OCCT HelixBRep_BuilderHelix::BuildPart (L261-427).
    /// Returns the part edges (the OCCT part wire) and its first/last
    /// vertices (TopExp::Vertices result).
    #[allow(clippy::too_many_arguments)]
    fn build_part(
        &mut self,
        the_axis: &Ax1,
        the_p_start: DVec3,
        the_height: f64,
        the_pitch: f64,
        the_taper_angle: f64,
        b_is_clockwise: bool,
    ) -> (Vec<Shape>, Shape, Shape) {
        if self.my_error_status != 0 {
            // OCCT returns the empty wire untouched.
            return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
        }

        self.my_error_status = 0;
        self.my_warning_status = 0;

        // 1. check & prepare data.
        let a_tol_prec = self.my_tolerance;
        let mut a_tol_ang = 1.0e-7;
        // Validate input parameters.
        if the_taper_angle > std::f64::consts::FRAC_PI_2 - a_tol_ang {
            self.my_error_status = 13; // invalid TaperAngle value
            return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
        }
        if the_height < a_tol_prec {
            self.my_error_status = 12; // invalid Height value
            return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
        }
        if the_pitch < a_tol_prec {
            self.my_error_status = 11; // invalid Pitch value
            return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
        }

        let a_lin = Lin::from_ax1(the_axis);
        let a_dist = a_lin.distance(the_p_start);
        if a_dist < a_tol_prec {
            self.my_error_status = 10; // myPStart belongs to the myAxis
            return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
        }

        a_tol_ang = a_tol_prec / a_dist;
        let _a_tol_ang = a_tol_ang;
        let a_angle_start = 0.0f64;

        let a_dir = the_axis.direction;
        let a_vec1 = a_dir;
        let a_m0 = the_axis.location;
        let a_vec = the_p_start - a_m0;
        let a_dm = a_vec1.dot(a_vec);
        let a_m1 = a_m0 + a_dm * a_vec1;
        let a_vec_x = the_p_start - a_m1;
        let a_dir_x = a_vec_x.normalize_or_zero();
        let a_ax2 = rcad_kernel::math::gp::Ax2::new(a_m1, a_dir, a_dir_x);

        let a_two_pi = 2.0 * std::f64::consts::PI;
        let a_c1 = the_pitch / a_two_pi;
        let a_t0 = 0.0f64;
        let mut a_t1 = a_angle_start;

        let mut a_t2 = the_height / a_c1;

        // 2. compute.
        let mut a_bh = crate::helix::helix_geom::builder_helix::BuilderHelix::new();
        a_bh.set_position(&a_ax2);
        a_bh.set_curve_parameters(a_t0, a_t2, the_pitch, a_dist, the_taper_angle, b_is_clockwise);
        a_bh.set_tolerance(self.my_tolerance);
        a_bh.set_approx_parameters(self.my_continuity, self.my_max_degree, self.my_max_segments);

        a_bh.perform();
        let i_err = a_bh.error_status();
        if i_err != 0 {
            self.my_error_status = 2;
            return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
        }

        self.my_tol_reached = self.my_tol_reached.max(a_bh.tolerance_reached());
        // aSC.Assign(aBH.Curves()) — a mutable local copy (OCCT may prepend).
        let mut a_sc: Vec<BSplineCurve3> = a_bh.curves().to_vec();
        if a_t1 < 0.0 {
            let mut a_bh1 = crate::helix::helix_geom::builder_helix::BuilderHelix::new();
            a_bh1.set_position(&a_ax2);
            a_bh1.set_curve_parameters(a_t1, a_t0, the_pitch, a_dist, the_taper_angle, b_is_clockwise);
            a_bh1.set_tolerance(self.my_tolerance);
            a_bh1.set_approx_parameters(self.my_continuity, self.my_max_degree, self.my_max_segments);

            a_bh1.perform();
            let i_err = a_bh1.error_status();
            if i_err != 0 {
                self.my_error_status = 2;
                return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
            }

            self.my_tol_reached = self.my_tol_reached.max(a_bh1.tolerance_reached());
            let a_sc1 = a_bh1.curves();
            for c in a_sc1.iter().rev() {
                a_sc.insert(0, c.clone());
            }
        }

        let a_nb_c = a_sc.len();
        let mut the_part: Vec<Shape> = Vec::new();
        let mut a_v1: Option<Shape> = None;
        let mut a_v2: Option<Shape> = None;
        let mut first_shape: Option<Shape> = None;
        let mut last_shape: Option<Shape> = None;

        for (idx, a_c) in a_sc.iter().enumerate() {
            let i = idx + 1;
            if i == 1 {
                if a_t1 > 0.0 {
                    // OCCT trims the first curve to [aT1, LastParameter] via
                    // Geom_TrimmedCurve — the trimmed range affects the edge
                    // parameter range only.
                    a_t2 = a_c.last_parameter();
                    let range = [a_t1, a_t2];
                    // The trimmed curve keeps the same poles; represent it as
                    // a range restriction.
                    let mut c = a_c.clone();
                    c.knots = c.knots.clone();
                    let _ = range;
                    // (aT1 == 0 for every pipeline caller; full port of the
                    // trimmed curve is deferred until a caller needs it.)
                }
                a_t1 = a_c.first_parameter();
                let a_p1 = rcad_kernel::math::bspl::de_boor(
                    a_c.degree,
                    &a_c.knots,
                    &a_c.control_points,
                    &a_c.weights,
                    a_t1,
                );
                // aBB.MakeVertex(aV1, aP1, myTolReached); aV1 FORWARD.
                let mut v = self.make_vertex(a_p1);
                v.orientation = Orientation::Forward;
                a_v1 = Some(v);
            }

            a_t2 = a_c.last_parameter();
            let a_p2 = rcad_kernel::math::bspl::de_boor(
                a_c.degree,
                &a_c.knots,
                &a_c.control_points,
                &a_c.weights,
                a_t2,
            );
            let mut v2 = self.make_vertex(a_p2);
            v2.orientation = Orientation::Reversed;
            a_v2 = Some(v2);

            // aBME.Init(aC, aV1, aV2); bIsDone check.
            let (b_is_done, a_e) = self.make_edge(a_c, a_v1.as_ref().unwrap(), a_v2.as_ref().unwrap());
            if !b_is_done {
                self.my_error_status = 3;
                return (Vec::new(), Shape::new(null_tshape(), 0, Orientation::Forward), Shape::new(null_tshape(), 0, Orientation::Forward));
            }
            // aBB.UpdateEdge(aE, myTolReached); aBB.Add(thePart, aE).
            let _ = &a_e;
            the_part.push(a_e.clone());
            if first_shape.is_none() {
                first_shape = a_v1.clone();
            }
            last_shape = a_v2.clone();

            // aV1 = aV2; aV1.Orientation(TopAbs_FORWARD).
            let mut v = a_v2.clone().unwrap();
            v.orientation = Orientation::Forward;
            a_v1 = Some(v);
        }

        if self.my_tol_reached > self.my_tolerance {
            self.my_warning_status = 1;
        }

        (
            the_part,
            first_shape.unwrap_or_else(|| Shape::new(null_tshape(), 0, Orientation::Forward)),
            last_shape.unwrap_or_else(|| Shape::new(null_tshape(), 0, Orientation::Forward)),
        )
    }

    /// BRep_Builder::MakeVertex equivalent on the builder's pool (vertex
    /// identity by quantized position — see module docs).
    fn make_vertex(&mut self, p: DVec3) -> Shape {
        self.brep.add_tvertex(p)
    }

    /// BRepBuilderAPI_MakeEdge::Init(C, V1, V2) + BRep_Builder::UpdateEdge
    /// equivalent: creates a new edge TShape over the BSpline with the
    /// [First, Last] parameter range and tolerance `my_tol_reached`.
    fn make_edge(&mut self, c: &BSplineCurve3, v1: &Shape, v2: &Shape) -> (bool, Shape) {
        let curve = Curve3::BSpline(c.clone());
        let range = [c.first_parameter(), c.last_parameter()];
        let index = self.brep.tshapes.len();
        let ts = arc_tshape_edge(&curve, v1, v2, range, self.my_tol_reached);
        self.brep.tshapes.push(ts.clone());
        let mut e = Shape::from_parts(ts, index, 0, Orientation::Forward);
        e.orientation = Orientation::Forward;
        (true, e)
    }

    /// OCCT HelixBRep_BuilderHelix::Smoothing (L431-465).
    fn smoothing(&mut self, the_parts: &mut Vec<Vec<Shape>>) {
        if the_parts.len() == 1 {
            return;
        }

        // BRepTools_WireExplorer over the first part wire — the edges are
        // stored in connected order, so "Current at end" is the last edge.
        let mut a_prev_edge = the_parts[0].last().cloned().unwrap();
        let mut prev_wire = 0usize;

        for wi in 1..the_parts.len() {
            let mut a_next_edge = the_parts[wi].first().cloned().unwrap();

            // Smoothing curves (UpdateEdge may replace the edge TShapes).
            let (new_prev, new_next) = self.smoothing_edges(&a_prev_edge, &a_next_edge);
            // Store the updated edges back into their wires.
            let last_idx = the_parts[prev_wire].len() - 1;
            the_parts[prev_wire][last_idx] = new_prev;
            the_parts[wi][0] = new_next;

            // OCCT: re-walk the wire; aPrevEdge = its last edge.
            a_prev_edge = the_parts[wi].last().cloned().unwrap();
            prev_wire = wi;
        }
    }

    /// OCCT HelixBRep_BuilderHelix::SmoothingEdges (L469-551).  Returns the
    /// (possibly curve-replaced) previous and next edges.
    fn smoothing_edges(&mut self, the_prev: &Shape, the_next: &Shape) -> (Shape, Shape) {
        const EPS_ANG: f64 = 1.0e-7;

        let prev_ed = match the_prev.data.as_ref() {
            TShape::Edge(ed) => ed.clone(),
            _ => return (the_prev.clone(), the_next.clone()),
        };
        let next_ed = match the_next.data.as_ref() {
            TShape::Edge(ed) => ed.clone(),
            _ => return (the_prev.clone(), the_next.clone()),
        };

        let (f1, l1) = (prev_ed.range[0], prev_ed.range[1]);
        let (f2, l2) = (next_ed.range[0], next_ed.range[1]);
        let mut a_cprev = match &prev_ed.curve {
            Some(Curve3::BSpline(b)) => b.clone(),
            _ => return (the_prev.clone(), the_next.clone()),
        };
        let mut a_cnext = match &next_ed.curve {
            Some(Curve3::BSpline(b)) => b.clone(),
            _ => return (the_prev.clone(), the_next.clone()),
        };

        let (p1, mut v1) = {
            let p = rcad_kernel::math::bspl::de_boor(
                a_cprev.degree,
                &a_cprev.knots,
                &a_cprev.control_points,
                &a_cprev.weights,
                l1,
            );
            let v = a_cprev.derivative_at(l1);
            (p, v)
        };
        let (p2, mut v2) = {
            let p = rcad_kernel::math::bspl::de_boor(
                a_cnext.degree,
                &a_cnext.knots,
                &a_cnext.control_points,
                &a_cnext.weights,
                f2,
            );
            let v = a_cnext.derivative_at(f2);
            (p, v)
        };

        if angle(v1, v2) < EPS_ANG {
            return (the_prev.clone(), the_next.clone());
        }

        v1 = 0.5 * (v1 + v2);
        v2 = v1;

        let a_deg_max = rcad_kernel::geom::bspline_ops::BSPLINE_MAX_DEGREE as i32;
        let a_deg = a_cprev.degree as i32;
        let mut b_prev_ok = false;
        let mut b_next_ok = false;

        let mut an_error_status =
            a_cprev.move_point_and_tangent(l1, p1, v1, self.my_tolerance, 1, -1);
        // OCCT initializes anErrorStatus = 1 and calls MovePointAndTangent —
        // the call result above matches; keep the retry loop.
        if an_error_status != 0 {
            for i in (a_deg + 1)..=a_deg_max {
                a_cprev.increase_degree(i as usize);
                an_error_status =
                    a_cprev.move_point_and_tangent(l1, p1, v1, self.my_tolerance, 1, -1);
                if an_error_status == 0 {
                    b_prev_ok = true;
                    break;
                }
            }
        } else {
            b_prev_ok = true;
        }

        let mut the_prev_out = the_prev.clone();
        if b_prev_ok {
            the_prev_out = self.replace_edge_curve(the_prev, &a_cprev, prev_ed.tolerance);
        }

        let a_deg = a_cnext.degree as i32;
        let mut an_error_status =
            a_cnext.move_point_and_tangent(f2, p2, v2, self.my_tolerance, -1, 1);
        if an_error_status != 0 {
            for i in (a_deg + 1)..=a_deg_max {
                a_cnext.increase_degree(i as usize);
                an_error_status =
                    a_cnext.move_point_and_tangent(f2, p2, v2, self.my_tolerance, -1, 1);
                if an_error_status == 0 {
                    b_next_ok = true;
                    break;
                }
            }
        } else {
            b_next_ok = true;
        }

        let mut the_next_out = the_next.clone();
        if b_next_ok {
            the_next_out = self.replace_edge_curve(the_next, &a_cnext, next_ed.tolerance);
        }

        (the_prev_out, the_next_out)
    }

    /// BRep_Builder::UpdateEdge(E, Curve, tol) equivalent — replaces the edge
    /// TShape with a copy carrying the new curve (see module docs).
    fn replace_edge_curve(&mut self, e: &Shape, c: &BSplineCurve3, tolerance: f64) -> Shape {
        let old = match e.data.as_ref() {
            TShape::Edge(ed) => ed.clone(),
            _ => return e.clone(),
        };
        let ts = Arc::new(TShape::Edge(TEdgeData {
            tolerance,
            curve: Some(Curve3::BSpline(c.clone())),
            ..old
        }));
        let index = self.brep.tshapes.len();
        self.brep.tshapes.push(ts.clone());
        Shape::from_parts(ts, index, e.location, e.orientation)
    }
}


fn arc_tshape_edge(
    curve: &Curve3,
    v1: &Shape,
    v2: &Shape,
    range: [f64; 2],
    tolerance: f64,
) -> Arc<TShape> {
    // OCCT BRep_Builder vertex params: BRep_Tool::Parameter mapping computed
    // at the range ends (the vertices sit exactly on them).
    let vertex_params = {
        let mut vp = std::collections::HashMap::new();
        vp.insert(v1.ptr_id(), range[0]);
        vp.insert(v2.ptr_id(), range[1]);
        vp
    };
    Arc::new(TShape::Edge(TEdgeData {
        my_shapes: vec![v1.clone(), v2.clone()],
        flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
        curve: Some(curve.clone()),
        first: v1.clone(),
        last: v2.clone(),
        range,
        degenerated: false,
        pcurves: indexmap::IndexMap::new(),
        representations: Vec::new(),
        vertex_params,
        tolerance,
        same_parameter: true,
        same_range: true,
    }))
}

/// Null TShape placeholder (OCCT null TopoDS handle).
fn null_tshape() -> std::sync::Arc<TShape> {
    std::sync::Arc::new(TShape::Vertex(TVertexData {
        my_shapes: Vec::new(),
        flags: 0,
        point: DVec3::ZERO,
        tolerance: 0.0,
        points: Vec::new(),
    }))
}

/// BRep_Tool::Pnt(vertex).
fn vertex_point(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT gp_Vec::Angle — angle in [0, PI] between two vectors.
fn angle(a: DVec3, b: DVec3) -> f64 {
    const RESOLUTION: f64 = 1.0e-12; // gp::Resolution()
    if a.length() <= RESOLUTION || b.length() <= RESOLUTION {
        panic!("gp_Vec::Angle");
    }
    let cos_ang = (a.dot(b) / (a.length() * b.length())).clamp(-1.0, 1.0);
    cos_ang.acos()
}
