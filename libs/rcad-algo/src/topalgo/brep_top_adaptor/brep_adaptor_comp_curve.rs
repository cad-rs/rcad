//! OCCT BRepAdaptor_CompCurve (TKBRep/BRepAdaptor) — 1:1 translation of
//! `BRepAdaptor_CompCurve.hxx` (L17-188) + `BRepAdaptor_CompCurve.cxx`
//! (L39-597).
//!
//! Architecture differences:
//! - `occ::handle<NCollection_HArray1<BRepAdaptor_Curve>> myCurves` and
//!   `NCollection_HArray1<double> myKnots` map to `Vec`s; the OCCT 1-based
//!   `Value(i)` indexing maps to `[i - 1]` (the established HArray1
//!   reduction, cf. BRepFill_LocationLaw myLaws).
//! - The member element type is the kernel BRepAdaptor_Curve re-host
//!   `rcad_kernel::base::proj_lib::brep_adaptor::BRepAdaptorCurve` (the
//!   GeomAdaptor_TransformedCurve-based edge adaptor, including the
//!   CurveOnSurface fallback and the edge location transformation); the
//!   inherited-member reads below go through its `Adaptor3dCurve` /
//!   `Adaptor3dCurveGeom` implementations.
//! - `GCPnts_AbscissaPoint::Length(C)` (cxx L140) maps to the kernel
//!   [`arc_length`] bridge over the member's underlying curve (the
//!   established GCPnts mapping of brep_fill_location_law.rs).  The
//!   curvilinear-abscissa form over a CurveOnSurface member keeps the OCCT
//!   failure path (the kernel length bridge is Curve3-only).
//! - `TopExp::CommonVertex(E1, E2, V)` (TopExp.cxx L316-334) is mirrored as
//!   the local [`top_exp_common_vertex`] over the shape-only vertex
//!   re-hosts; `TopExp::LastVertex` maps to `top_exp_last_vertex`.
//! - `myWire.Closed()` maps to the wire CLOSED flag (the pool flag read,
//!   cf. BRepFill_LocationLaw::myPath.Closed()).
//! - The adaptor implements the kernel [`Adaptor3dCurve`] /
//!   [`Adaptor3dCurveGeom`] traits — the rcad `Adaptor3d_Curve` instance
//!   surface (`occ::handle<Adaptor3d_Curve>` = `Arc<dyn Adaptor3dCurve>`).
//!   In OCCT the GeomFill consumers hold
//!   `handle(BRepAdaptor_CompCurve)` directly; the current rcad geomfill
//!   ports consume `Curve3` — see the GAP note on
//!   [`BRepAdaptorCompCurve::comp_curve`].

use std::sync::Arc;

use glam::DVec3;

use rcad_kernel::base::gcpnts::abscissa_point::arc_length;
use rcad_kernel::base::proj_lib::adaptor::Adaptor3dCurve;
use rcad_kernel::base::proj_lib::brep_adaptor::BRepAdaptorCurve;
use rcad_kernel::base::proj_lib::geom_adaptor_curve::Adaptor3dCurveGeom;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::Curve3;
use rcad_kernel::math::{GeomAbsCurveType, GeomAbsShape};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape};

use crate::feat::loc_ope_wires_on_shape::{top_exp_last_vertex, top_exp_vertices};
use crate::topalgo::brep_tools_wire_explorer::WireExplorer;

/// OCCT TopExp::CommonVertex(E1, E2, V) (TopExp.cxx L316-334) — the null
/// result maps to None.
fn top_exp_common_vertex(e1: &Shape, e2: &Shape) -> Option<Shape> {
    let (v1, v2) = top_exp_vertices(e1);
    let (v3, v4) = top_exp_vertices(e2);
    if shapes_are_same_opt(&v1, &v3) || shapes_are_same_opt(&v1, &v4) {
        return v1;
    }
    if shapes_are_same_opt(&v2, &v3) || shapes_are_same_opt(&v2, &v4) {
        return v2;
    }
    None
}

/// OCCT TopoDS_Shape::IsSame — the TShape identity test (orientation
/// ignored); None models the OCCT null shape (never same).
fn shapes_are_same_opt(a: &Option<Shape>, b: &Option<Shape>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.ptr_id() == b.ptr_id(),
        _ => false,
    }
}

// ===========================================================================
// BRepAdaptor_CompCurve
// ===========================================================================

/// OCCT BRepAdaptor_CompCurve (hxx L56-186) — the Adaptor3d_Curve over a
/// wire: the edges are sorted along a parametrization whose knots are
/// either the edge indices (KnotByCurvilinearAbcissa = false) or the
/// cumulative curvilinear abscissas.  C0/C1 continuities are not assumed;
/// the class cannot be periodic at all.
pub struct BRepAdaptorCompCurve {
    /// The BRep pool shared state (OCCT keeps the TopoDS handles alive;
    /// rcad: the refcounted pool copy, cf. BRepCurve2d).
    brep: Arc<BRep>,
    /// OCCT myWire.
    my_wire: Shape,
    /// OCCT TFirst.
    t_first: f64,
    /// OCCT TLast.
    t_last: f64,
    /// OCCT PTol.
    p_tol: f64,
    /// OCCT myCurves (NCollection_HArray1<BRepAdaptor_Curve>, 1-based) —
    /// `Value(i)` maps to `my_curves[i - 1]`.
    my_curves: Vec<BRepAdaptorCurve>,
    /// OCCT myKnots (NCollection_HArray1<double>, 1-based, length
    /// NbEdge + 1).
    my_knots: Vec<f64>,
    /// OCCT CurIndex.
    cur_index: i32,
    /// OCCT Forward.
    forward: bool,
    /// OCCT IsbyAC.
    is_by_ac: bool,
}

impl Clone for BRepAdaptorCompCurve {
    /// OCCT ShallowCopy (cxx L79-100) — the member-wise copy.
    fn clone(&self) -> Self {
        BRepAdaptorCompCurve {
            brep: self.brep.clone(),
            my_wire: self.my_wire.clone(),
            t_first: self.t_first,
            t_last: self.t_last,
            p_tol: self.p_tol,
            my_curves: self.my_curves.clone(),
            my_knots: self.my_knots.clone(),
            cur_index: self.cur_index,
            forward: self.forward,
            is_by_ac: self.is_by_ac,
        }
    }
}

/// The 1-based OCCT knot reader (`myKnots->Value(i)`).
fn knot_at(knots: &[f64], i: i32) -> f64 {
    knots[(i - 1) as usize]
}

/// The 1-based OCCT knot writer (`myKnots->SetValue(i, v)`).
fn set_knot(knots: &mut [f64], i: i32, v: f64) {
    knots[(i - 1) as usize] = v;
}

impl BRepAdaptorCompCurve {
    /// OCCT BRepAdaptor_CompCurve() (cxx L39-47) — the undefined curve
    /// with no wire loaded.
    pub fn undefined(brep: &BRep) -> Self {
        BRepAdaptorCompCurve {
            brep: Arc::new(brep.clone()),
            my_wire: Shape::null(),
            t_first: 0.0,
            t_last: 0.0,
            p_tol: 0.0,
            my_curves: Vec::new(),
            my_knots: Vec::new(),
            cur_index: -1,
            forward: false,
            is_by_ac: false,
        }
    }

    /// OCCT BRepAdaptor_CompCurve(theWire, theIsAC = false)
    /// (cxx L49-59).
    pub fn new(brep: &BRep, the_wire: &Shape) -> Self {
        let mut c = BRepAdaptorCompCurve::undefined(brep);
        c.initialize(the_wire, false);
        c
    }

    /// OCCT BRepAdaptor_CompCurve(W, KnotByCurvilinearAbcissa) (cxx L49-59)
    /// with the flag spelled out.
    pub fn new_with_abscissa_flag(brep: &BRep, the_wire: &Shape, the_is_ac: bool) -> Self {
        let mut c = BRepAdaptorCompCurve::undefined(brep);
        c.initialize(the_wire, the_is_ac);
        c
    }

    /// OCCT BRepAdaptor_CompCurve(theWire, theIsAC, theFirst, theLast,
    /// theTolerance) (cxx L61-75).
    pub fn new_with_window(
        brep: &BRep,
        the_wire: &Shape,
        the_is_ac: bool,
        the_first: f64,
        the_last: f64,
        the_tolerance: f64,
    ) -> Self {
        let mut c = BRepAdaptorCompCurve {
            brep: Arc::new(brep.clone()),
            my_wire: the_wire.clone(),
            t_first: the_first,
            t_last: the_last,
            p_tol: the_tolerance,
            my_curves: Vec::new(),
            my_knots: Vec::new(),
            cur_index: -1,
            forward: false,
            is_by_ac: the_is_ac,
        };
        c.initialize_with_window(the_wire, the_is_ac, the_first, the_last, the_tolerance);
        c
    }

    /// OCCT Initialize(W, AC) (cxx L102-174).
    pub fn initialize(&mut self, w: &Shape, ac: bool) {
        let mut nb_edge: usize = 0;

        self.my_wire = w.clone();
        self.p_tol = 0.0;
        self.is_by_ac = ac;

        // Count the non-degenerated edges along the wire explorer.
        let mut wexp = WireExplorer::new();
        // OCCT cxx L112: wexp.Init(myWire) — Init(W) delegates to
        // Init(W, Face()) (BRepTools_WireExplorer.cxx L83-89).
        wexp.init_wire_face(&self.my_wire, &Shape::null());
        while wexp.more() {
            if !self.brep.edge(wexp.current().clone()).degenerated {
                nb_edge += 1;
            }
            wexp.next();
        }

        if nb_edge == 0 {
            return;
        }

        self.cur_index = ((nb_edge + 1) / 2) as i32;
        self.my_curves = Vec::with_capacity(nb_edge);
        self.my_knots = vec![0.0; nb_edge + 1];
        set_knot(&mut self.my_knots, 1, 0.0);

        let mut ii: usize = 0;
        let mut wexp = WireExplorer::new();
        wexp.init_wire_face(&self.my_wire, &Shape::null());
        while wexp.more() {
            let e = wexp.current().clone();
            if !self.brep.edge(e.clone()).degenerated {
                ii += 1;
                // myCurves->ChangeValue(ii).Initialize(E).
                self.my_curves.push(BRepAdaptorCurve::with_edge(&self.brep, &e));
                if ac {
                    // myKnots->SetValue(ii + 1, myKnots->Value(ii)).
                    let prev_knot = knot_at(&self.my_knots, ii as i32);
                    set_knot(&mut self.my_knots, ii as i32 + 1, prev_knot);
                    // myKnots->ChangeValue(ii + 1) +=
                    //   GCPnts_AbscissaPoint::Length(myCurves->ChangeValue(ii))
                    //   — the member's full-range length.
                    let c = &self.my_curves[ii - 1];
                    if c.transformed.is_curve_on_surface() {
                        panic!(
                            "GAP: GCPnts_AbscissaPoint::Length over a \
                             CurveOnSurface member (BRepAdaptor_CompCurve.cxx \
                             L140) — the kernel arc_length bridge is Curve3-only"
                        )
                    }
                    let len = arc_length(
                        c.transformed.geom_curve(),
                        c.first_parameter(),
                        c.last_parameter(),
                    );
                    self.my_knots[ii] += len;
                } else {
                    // myKnots->SetValue(ii + 1, (double)ii).
                    set_knot(&mut self.my_knots, ii as i32 + 1, ii as f64);
                }
            }
            wexp.next();
        }

        self.forward = true; // Default ; The Reverse Edges are parsed.
        if (nb_edge > 2) || ((nb_edge == 2) && (!self.wire_is_closed())) {
            let or = self.my_curves[0].edge().orientation;
            // TopExp::CommonVertex(myCurves->Value(1).Edge(),
            //   myCurves->Value(2).Edge(), VI).
            let vi = top_exp_common_vertex(self.my_curves[0].edge(), self.my_curves[1].edge());
            // TopExp::LastVertex(myCurves->Value(1).Edge()).
            let vl = top_exp_last_vertex(self.my_curves[0].edge());
            if shapes_are_same_opt(&vi, &vl) {
                // The direction of parsing is always preserved
                if or == Orientation::Reversed {
                    self.forward = false;
                }
            } else {
                // The direction of parsing is always reversed
                if or != Orientation::Reversed {
                    self.forward = false;
                }
            }
        }

        self.t_first = 0.0;
        self.t_last = knot_at(&self.my_knots, self.my_knots.len() as i32);
    }

    /// OCCT Initialize(W, AC, First, Last, Tol) (cxx L176-235).
    pub fn initialize_with_window(
        &mut self,
        w: &Shape,
        ac: bool,
        first: f64,
        last: f64,
        tol: f64,
    ) {
        self.initialize(w, ac);
        self.t_first = first;
        self.t_last = last;
        self.p_tol = tol;

        // Trim the extremal curves.
        let mut i1 = self.cur_index;
        let mut i2 = self.cur_index;
        let mut f = self.t_first;
        let mut l = self.t_last;
        let mut d = 0.0;
        self.prepare(&mut f, &mut d, &mut i1);
        self.prepare(&mut l, &mut d, &mut i2);
        self.cur_index = (i1 + i2) / 2; // Small optimization
        if i1 == i2 {
            if l > f {
                // HC = down_cast<BRepAdaptor_Curve>(Value(i1).Trim(f, l, PTol)).
                let hc = self.curve_at(i1).trim_of(f, l, self.p_tol);
                self.set_curve_at(i1, hc);
            } else {
                let hc = self.curve_at(i1).trim_of(l, f, self.p_tol);
                self.set_curve_at(i1, hc);
            }
        } else {
            let mut k = self.curve_at(i1).last_parameter();
            if k > f {
                let hc = self.curve_at(i1).trim_of(f, k, self.p_tol);
                self.set_curve_at(i1, hc);
            } else {
                let hc = self.curve_at(i1).trim_of(k, f, self.p_tol);
                self.set_curve_at(i1, hc);
            }

            k = self.curve_at(i2).first_parameter();
            if k <= l {
                let hc = self.curve_at(i2).trim_of(k, l, self.p_tol);
                self.set_curve_at(i2, hc);
            } else {
                let hc = self.curve_at(i2).trim_of(l, k, self.p_tol);
                self.set_curve_at(i2, hc);
            }
        }
    }

    /// OCCT myWire.Closed() — the wire CLOSED flag.
    fn wire_is_closed(&self) -> bool {
        self.brep
            .has_flag(self.my_wire.clone(), tshape_flags::CLOSED)
    }

    /// OCCT myCurves->Value(i) — the 1-based member access.
    fn curve_at(&self, i: i32) -> &BRepAdaptorCurve {
        &self.my_curves[(i - 1) as usize]
    }

    /// OCCT myCurves->SetValue(i, c).
    fn set_curve_at(&mut self, i: i32, c: BRepAdaptorCurve) {
        self.my_curves[(i - 1) as usize] = c;
    }

    /// OCCT Wire() (cxx L237-240).
    pub fn wire(&self) -> &Shape {
        &self.my_wire
    }

    /// OCCT Edge(U, E, UonE) (cxx L242-249) — returns the edge and the
    /// parameter on it corresponding to U.
    pub fn edge(&self, u: f64) -> (Shape, f64) {
        let mut uon_e = u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut uon_e, &mut d, &mut index);
        (self.curve_at(index).edge().clone(), uon_e)
    }

    /// OCCT FirstParameter() (cxx L251-254).
    pub fn first_parameter(&self) -> f64 {
        self.t_first
    }

    /// OCCT LastParameter() (cxx L256-259).
    pub fn last_parameter(&self) -> f64 {
        self.t_last
    }

    /// OCCT Continuity() (cxx L261-268).
    pub fn continuity(&self) -> GeomAbsShape {
        if self.my_curves.len() > 1 {
            return GeomAbsShape::C0;
        }
        self.curve_at(1).continuity()
    }

    /// OCCT NbIntervals(S) (cxx L270-279).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let mut nb_int = 0usize;
        for c in &self.my_curves {
            nb_int += c.nb_intervals(s);
        }
        nb_int
    }

    /// OCCT Intervals(T, S) (cxx L281-335) — T is the caller-provided
    /// array (T.Length() > NbIntervals()); the OCCT 1-based T(kk) indexing
    /// maps to `t[kk - 1]`.
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let mut kk: usize;

        // First curve (direction of parsing of the edge)
        let n = self.curve_at(1).nb_intervals(s);
        // Ti = new HArray1(1, n + 1) — the boundary array (NbIntervals + 1
        // entries), filled by the member intervals read.
        let mut ti = self.curve_at(1).intervals(s);
        debug_assert_eq!(ti.len(), n + 1);
        let mut f = 0.0;
        let mut delta = 0.0;
        self.inv_prepare(1, &mut f, &mut delta);
        let mut big_f = knot_at(&self.my_knots, 1);
        if delta < 0.0 {
            // invert the direction of parsing
            kk = 1;
            let mut jj = ti.len() as i32;
            while jj > 0 {
                t[kk - 1] = big_f + (ti[(jj - 1) as usize] - f) * delta;
                kk += 1;
                jj -= 1;
            }
        } else {
            kk = 1;
            while kk <= ti.len() {
                t[kk - 1] = big_f + (ti[kk - 1] - f) * delta;
                kk += 1;
            }
        }

        // and the next
        for ii in 2..=self.my_curves.len() {
            let n = self.curve_at(ii as i32).nb_intervals(s);
            if n != ti.len() - 1 {
                // Ti = new HArray1(1, n + 1).
                ti = vec![0.0; n + 1];
            }
            ti.copy_from_slice(&self.curve_at(ii as i32).intervals(s));
            self.inv_prepare(ii as i32, &mut f, &mut delta);
            big_f = knot_at(&self.my_knots, ii as i32);
            if delta < 0.0 {
                // invert the direction of parcing
                let mut jj = ti.len() as i32 - 1;
                while jj > 0 {
                    t[kk - 1] = big_f + (ti[(jj - 1) as usize] - f) * delta;
                    kk += 1;
                    jj -= 1;
                }
            } else {
                for jj in 2..=ti.len() {
                    t[kk - 1] = big_f + (ti[jj - 1] - f) * delta;
                    kk += 1;
                }
            }
        }
    }

    /// OCCT Trim(First, Last, Tol) (cxx L337-344).
    pub fn trim(&self, first: f64, last: f64, tol: f64) -> BRepAdaptorCompCurve {
        BRepAdaptorCompCurve::new_with_window(
            &self.brep,
            &self.my_wire,
            self.is_by_ac,
            first,
            last,
            tol,
        )
    }

    /// OCCT IsClosed() (cxx L346-349).
    pub fn is_closed(&self) -> bool {
        self.wire_is_closed()
    }

    /// OCCT IsPeriodic() (cxx L351-354).
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT Period() (cxx L356-359).
    pub fn period(&self) -> f64 {
        self.t_last - self.t_first
    }

    /// OCCT EvalD0(theU) (cxx L361-367).
    pub fn eval_d0(&self, the_u: f64) -> DVec3 {
        let mut u = the_u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut u, &mut d, &mut index);
        self.curve_at(index).d0(u)
    }

    /// OCCT EvalD1(theU) (cxx L369-377) — (Point, D1) with the chain-rule
    /// factor (`aRes.D1 *= d`).
    pub fn eval_d1(&self, the_u: f64) -> (DVec3, DVec3) {
        let mut u = the_u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut u, &mut d, &mut index);
        let (p, d1) = self.curve_at(index).d1(u);
        (p, d1 * d)
    }

    /// OCCT EvalD2(theU) (cxx L379-388) — (Point, D1, D2).
    pub fn eval_d2(&self, the_u: f64) -> (DVec3, DVec3, DVec3) {
        let mut u = the_u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut u, &mut d, &mut index);
        let (p, d1, d2) = self.curve_at(index).d2(u);
        (p, d1 * d, d2 * (d * d))
    }

    /// OCCT EvalD3(theU) (cxx L390-400) — (Point, D1, D2, D3).
    pub fn eval_d3(&self, the_u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        let mut u = the_u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut u, &mut d, &mut index);
        let (p, d1, d2, d3) = self.curve_at(index).d3(u);
        (p, d1 * d, d2 * (d * d), d3 * (d * d * d))
    }

    /// OCCT EvalDN(theU, theN) (cxx L402-408) — `EvalDN(u, theN) *
    /// std::pow(d, theN)`.
    pub fn eval_dn(&self, the_u: f64, the_n: i32) -> DVec3 {
        let mut u = the_u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut u, &mut d, &mut index);
        self.curve_at(index).dn(u, the_n) * d.powi(the_n)
    }

    /// OCCT Resolution(R3d) (cxx L410-423) — the minimum member resolution.
    pub fn resolution(&self, r3d: f64) -> f64 {
        let mut res = 1.0e200;
        for c in &self.my_curves {
            let r = c.resolution(r3d);
            if r < res {
                res = r;
            }
        }
        res
    }

    /// OCCT GetType() (cxx L425-430).
    pub fn get_type(&self) -> GeomAbsCurveType {
        GeomAbsCurveType::OtherCurve // temporary
    }

    /// OCCT Line() (cxx L432-435).
    pub fn line(&self) -> rcad_kernel::geom::Line3 {
        self.curve_at(1).line()
    }

    /// OCCT Circle() (cxx L437-440).
    pub fn circle(&self) -> rcad_kernel::geom::Circle3 {
        self.curve_at(1).circle()
    }

    /// OCCT Ellipse() (cxx L442-445).
    pub fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
        self.curve_at(1).ellipse()
    }

    /// OCCT Hyperbola() (cxx L447-450).
    pub fn hyperbola(&self) -> rcad_kernel::geom::Hyperbola3 {
        self.curve_at(1).hyperbola()
    }

    /// OCCT Parabola() (cxx L452-455).
    pub fn parabola(&self) -> rcad_kernel::geom::Parabola3 {
        self.curve_at(1).parabola()
    }

    /// OCCT Degree() (cxx L457-460).
    pub fn degree(&self) -> usize {
        self.curve_at(1).degree()
    }

    /// OCCT IsRational() (cxx L462-465).
    pub fn is_rational(&self) -> bool {
        self.curve_at(1).is_rational()
    }

    /// OCCT NbPoles() (cxx L467-470).
    pub fn nb_poles(&self) -> usize {
        self.curve_at(1).nb_poles()
    }

    /// OCCT NbKnots() (cxx L472-475).
    pub fn nb_knots(&self) -> usize {
        self.curve_at(1).nb_knots()
    }

    /// OCCT Bezier() (cxx L477-480).
    pub fn bezier(&self) -> rcad_kernel::geom::BezierCurve3 {
        self.curve_at(1).bezier()
    }

    /// OCCT BSpline() (cxx L482-485).
    pub fn bspline(&self) -> rcad_kernel::geom::BSplineCurve3 {
        self.curve_at(1).bspline()
    }

    /// OCCT ShallowCopy() (cxx L79-100) — the member-wise copy.
    pub fn shallow_copy(&self) -> BRepAdaptorCompCurve {
        self.clone()
    }

    /// OCCT Prepare(W, Delta, theCurIndex) (cxx L496-568) — locate the
    /// parameter in the knot array and reparametrize onto the member edge
    /// window.  When the parameter is close to a node the rule is
    /// determined by the sign of the tolerance offset: negative keeps the
    /// rule preceding the node, positive the rule following it.
    fn prepare(&self, w: &mut f64, delta: &mut f64, the_cur_index: &mut i32) {
        let eps = if *w - self.t_first < self.t_last - *w {
            self.p_tol
        } else {
            -self.p_tol
        };

        let wtest = *w + eps; // Offset to discriminate the nodes

        // Find the index
        let mut trouve = false;
        if knot_at(&self.my_knots, *the_cur_index) > wtest {
            let mut ii = *the_cur_index - 1;
            while ii > 0 && !trouve {
                if knot_at(&self.my_knots, ii) <= wtest {
                    *the_cur_index = ii;
                    trouve = true;
                }
                ii -= 1;
            }
            if !trouve {
                *the_cur_index = 1; // Out of limits...
            }
        } else if knot_at(&self.my_knots, *the_cur_index + 1) <= wtest {
            let mut ii = *the_cur_index + 1;
            while ii <= self.my_curves.len() as i32 && !trouve {
                if knot_at(&self.my_knots, ii + 1) > wtest {
                    *the_cur_index = ii;
                    trouve = true;
                }
                ii += 1;
            }
            if !trouve {
                *the_cur_index = self.my_curves.len() as i32; // Out of limits...
            }
        }

        // Invert ?
        let e = self.curve_at(*the_cur_index).edge();
        let or = e.orientation;
        let reverse = (self.forward && (or == Orientation::Reversed))
            || (!self.forward && (or != Orientation::Reversed));

        // Calculate the local parameter
        // BRep_Tool::Range(E, f, l).
        let range = self.brep.edge(e.clone()).range;
        let f = range[0];
        let l = range[1];
        *delta =
            knot_at(&self.my_knots, *the_cur_index + 1) - knot_at(&self.my_knots, *the_cur_index);
        if *delta > self.p_tol * 1.0e-9 {
            *delta = (l - f) / *delta;
        }

        if reverse {
            *delta *= -1.0;
            *w = l + (*w - knot_at(&self.my_knots, *the_cur_index)) * *delta;
        } else {
            *w = f + (*w - knot_at(&self.my_knots, *the_cur_index)) * *delta;
        }
    }

    /// OCCT InvPrepare(index, First, Delta) (cxx L570-597) — the inverse
    /// reparametrization data: T = Ti + (t - First) * Delta.
    fn inv_prepare(&self, index: i32, first: &mut f64, delta: &mut f64) {
        // Invert?
        let e = self.curve_at(index).edge();
        let or = e.orientation;
        let reverse = (self.forward && (or == Orientation::Reversed))
            || (!self.forward && (or != Orientation::Reversed));

        // Calculate the parameters of reparametrisation
        // such as : T = Ti + (t-First)*Delta
        // BRep_Tool::Range(E, f, l).
        let range = self.brep.edge(e.clone()).range;
        let f = range[0];
        let l = range[1];
        *delta = knot_at(&self.my_knots, index + 1) - knot_at(&self.my_knots, index);
        if l - f > self.p_tol * 1.0e-9 {
            *delta /= l - f;
        }

        if reverse {
            *delta *= -1.0;
            *first = l;
        } else {
            *first = f;
        }
    }
}

// ===========================================================================
// The Adaptor3d_Curve instance surface (the kernel trait pair)
// ===========================================================================

impl Adaptor3dCurve for BRepAdaptorCompCurve {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64 {
        self.first_parameter()
    }

    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64 {
        self.last_parameter()
    }

    /// OCCT Value(U) / D0(U, P) — the EvalD0 composition.
    fn value(&self, u: f64) -> DVec3 {
        self.eval_d0(u)
    }

    /// OCCT D1(U, P, V1).
    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.eval_d1(u)
    }

    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.eval_d2(u)
    }

    /// OCCT D3(U, P, V1, V2, V3).
    fn d3(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        self.eval_d3(u)
    }

    /// OCCT DN(U, N).
    fn dn(&self, u: f64, n: i32) -> DVec3 {
        self.eval_dn(u, n)
    }

    /// OCCT Continuity().
    fn continuity(&self) -> GeomAbsShape {
        self.continuity()
    }

    /// OCCT NbIntervals(S).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.nb_intervals(s)
    }

    /// OCCT Intervals(T, S) — the sized trait encoding over the array form.
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        let mut t = vec![0.0; self.nb_intervals(s)];
        self.intervals(&mut t, s);
        t
    }

    /// OCCT Trim(First, Last, Tol).
    fn trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.trim(first, last, tol))
    }

    /// OCCT IsClosed().
    fn is_closed(&self) -> bool {
        self.is_closed()
    }

    /// OCCT IsPeriodic().
    fn is_periodic(&self) -> bool {
        self.is_periodic()
    }

    /// OCCT Period().
    fn period(&self) -> f64 {
        self.period()
    }

    /// OCCT Resolution(R3d).
    fn resolution(&self, r3d: f64) -> f64 {
        self.resolution(r3d)
    }

    /// OCCT Degree().
    fn degree(&self) -> usize {
        self.degree()
    }

    /// OCCT IsRational().
    fn is_rational(&self) -> bool {
        self.is_rational()
    }

    /// OCCT NbPoles().
    fn nb_poles(&self) -> usize {
        self.nb_poles()
    }

    /// OCCT NbKnots().
    fn nb_knots(&self) -> usize {
        self.nb_knots()
    }

    /// OCCT Bezier().
    fn bezier(&self) -> rcad_kernel::geom::BezierCurve3 {
        self.bezier()
    }

    /// OCCT BSpline().
    fn bspline(&self) -> rcad_kernel::geom::BSplineCurve3 {
        self.bspline()
    }

    /// OCCT ShallowCopy().
    fn shallow_copy(&self) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.shallow_copy())
    }
}

impl Adaptor3dCurveGeom for BRepAdaptorCompCurve {
    /// OCCT GetType() — GeomAbs_OtherCurve (cxx L425-430).
    fn get_type(&self) -> CurveType {
        CurveType::Other
    }

    /// OCCT Line() — the curves(1) delegate.
    fn line(&self) -> rcad_kernel::geom::Line3 {
        self.line()
    }

    /// OCCT Circle().
    fn circle(&self) -> rcad_kernel::geom::Circle3 {
        self.circle()
    }

    /// OCCT Ellipse().
    fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
        self.ellipse()
    }

    /// OCCT Parabola().
    fn parabola(&self) -> rcad_kernel::geom::Parabola3 {
        self.parabola()
    }

    /// OCCT Hyperbola().
    fn hyperbola(&self) -> rcad_kernel::geom::Hyperbola3 {
        self.hyperbola()
    }

    /// OCCT Trim().
    fn trim_geom(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurveGeom> {
        Arc::new(self.trim(first, last, tol))
    }
}

impl BRepAdaptorCompCurve {
    /// OCCT Adaptor3d_Curve::D1(U, P, V) — the base-class out-parameter
    /// form over EvalD1 (Adaptor3d_Curve.hxx L106-114); the consumer form
    /// used by BRepFill_PipeShell::Set (BRepFill_PipeShell.cxx L392-393).
    pub fn d1(&self, u: f64, p: &mut DVec3, v: &mut DVec3) {
        let (pu, vu) = self.eval_d1(u);
        *p = pu;
        *v = vu;
    }

    /// The GeomFill-facing view of the comp curve.  In OCCT the consumers
    /// hold `handle(BRepAdaptor_CompCurve)` directly:
    /// - BRepFill_SectionPlacement::Perform (BRepFill_SectionPlacement.cxx
    ///   L167-168) passes it to GeomFill_SectionPlacement::Perform
    ///   (handle(Adaptor3d_Curve), Tol);
    /// - BRepFill_PipeShell::Set (BRepFill_PipeShell.cxx L399/417/426)
    ///   feeds it to GeomFill_GuideTrihedronAC / GeomFill_GuideTrihedronPlan.
    ///
    /// The rcad geomfill ports consume the kernel `Curve3` value instead of
    /// the `Adaptor3d_Curve` handle, and no exact `Curve3` exists for the
    /// wire comp curve (the OCCT knot parametrization 0..NbEdge with the
    /// per-edge reparametrization is adaptor state, not a Geom curve);
    /// materializing one would be an approximation, which the port
    /// discipline forbids.  This keeps the OCCT failure path at the same
    /// statement until the geomfill Perform/GuideTrihedron entries are
    /// re-hosted over `Arc<dyn Adaptor3dCurve>` (the consumer-side
    /// overload change, outside this unit).
    pub fn comp_curve(&self) -> Curve3 {
        panic!(
            "GAP: the GeomFill-facing Curve3 view of BRepAdaptor_CompCurve — \
             GeomFill consumers (GeomFill_SectionPlacement::Perform(Path, Tol) \
             at BRepFill_SectionPlacement.cxx L168; GeomFill_GuideTrihedronAC/Plan \
             at BRepFill_PipeShell.cxx L399) consume handle(Adaptor3d_Curve); \
             the rcad geomfill Perform/GuideTrihedron ports accept only Curve3 — \
             the handle-typed overload re-host is pending"
        )
    }
}
