//! OCCT BRepFill_Sweep.cxx, part C (members, cxx L2043-2077 + L2213-3207) —
//! CorrectApproxParameters and BuildShell, the companion of
//! [`super::brep_fill_sweep`] (part A) and [`super::brep_fill_sweep_b`]
//! (part-B statics).
//!
//! Architecture differences:
//! - `NCollection_Array2<T>` maps to the local 1-based [`Array2`].
//! - `NCollection_Map / NCollection_DataMap` keyed by TopTools_ShapeMapHasher
//!   map to HashSet/HashMap keyed by `ShapeKey`.
//! - `BRepAdaptor_Curve::Continuity()` / `NbIntervals(GeomAbs_CN)` map to the
//!   local GeomAdaptor_Curve re-hosts ([`geom_adaptor_curve_continuity`] /
//!   [`geom_adaptor_nb_intervals_cn`], GeomAdaptor_Curve.cxx L333-411).

use std::collections::HashSet;

use rcad_kernel::topo::topods::{BRep, BRepBuilder, GeomAbsShape, Orientation, Shape};

use crate::brep_fill::brep_fill_pipe_shell_b::{ShapeHArray2, ShapeToArray2Map};
use crate::brep_fill::brep_fill_sweep::BRepFillSweep;
use crate::brep_fill::generator::ShapeKey;

use rcad_kernel::geom::SurfaceEval;

use super::brep_fill_sweep::null_edge;
use super::brep_fill_sweep_b::is_degen;

// ---------------------------------------------------------------------------
// NCollection_Array2 — the 1-based 2d array
// ---------------------------------------------------------------------------

/// OCCT `NCollection_Array2<T>` (1-based row/col indexing mapped to
/// row-major Vec storage, `Value(i, j)` -> `data[i - r1][j - c1]`).
#[derive(Debug, Clone)]
pub(super) struct Array2<T: Clone> {
    r1: usize,
    c1: usize,
    data: Vec<Vec<T>>,
}

impl<T: Clone> Array2<T> {
    /// OCCT NCollection_Array2(LowerRow, UpperRow, LowerCol, UpperCol) with
    /// an initial value (the `Init` call of the OCCT declarations).
    pub fn new(r1: usize, r2: usize, c1: usize, c2: usize, init: T) -> Self {
        Array2 {
            r1,
            c1,
            data: vec![vec![init; c2 - c1 + 1]; r2 - r1 + 1],
        }
    }

    /// OCCT Value(i, j).
    pub fn get(&self, r: usize, c: usize) -> T {
        self.data[r - self.r1][c - self.c1].clone()
    }

    /// OCCT SetValue(i, j, V).
    pub fn set(&mut self, r: usize, c: usize, v: T) {
        self.data[r - self.r1][c - self.c1] = v;
    }

    /// OCCT ChangeValue(i, j).
    pub fn change_value(&mut self, r: usize, c: usize) -> &mut T {
        &mut self.data[r - self.r1][c - self.c1]
    }
}

// ---------------------------------------------------------------------------
// GeomAdaptor_Curve re-hosts (GeomAdaptor_Curve.cxx L333-411)
// ---------------------------------------------------------------------------

/// OCCT GeomAdaptor_Curve::Continuity (cxx L333-367): analytic curves are
/// CN; a BSpline computes LocalContinuity over its interior knots; a
/// trimmed curve dispatches on the basis (the OCCT adaptor holds the basis).
pub(super) fn geom_adaptor_curve_continuity(
    c: &rcad_kernel::geom::Curve3,
) -> rcad_kernel::math::GeomAbsShape {
    use rcad_kernel::geom::Curve3;
    match c {
        Curve3::BSpline(bs) => {
            // LocalContinuity: the min multiplicity-based continuity at the
            // interior knots (GeomAdaptor_Curve.cxx LocalContinuity).
            let degree = bs.degree as i32;
            let mut cont = 3.min(degree);
            // flat knots: count interior multiplicities
            let mut i = 0usize;
            let n = bs.knots.len();
            while i < n {
                let mut j = i;
                while j < n && (bs.knots[j] - bs.knots[i]).abs() < 1.0e-12 {
                    j += 1;
                }
                let mult = (j - i) as i32;
                let interior = i > 0 || j < n;
                if interior {
                    cont = cont.min(degree - mult);
                }
                i = j;
            }
            match cont.max(0) {
                0 => rcad_kernel::math::GeomAbsShape::C0,
                1 => rcad_kernel::math::GeomAbsShape::C1,
                2 => rcad_kernel::math::GeomAbsShape::C2,
                _ => rcad_kernel::math::GeomAbsShape::C3,
            }
        }
        // GeomAbs_OffsetCurve: GetBasisCurveContinuity downgrade (cxx L340-360)
        Curve3::Offset(_) => rcad_kernel::math::GeomAbsShape::C0,
        // OtherCurve throws; the analytic kinds return CN (cxx L366)
        _ => rcad_kernel::math::GeomAbsShape::CN,
    }
}

/// OCCT GeomAdaptor_Curve::NbIntervals(GeomAbs_CN) (cxx L371-411): 1 when
/// the curve is smooth enough; otherwise the interior knots with
/// multiplicity over degree - aCont (aCont = degree for CN) inside the
/// range (the BSplCLib::Intervals reduction).
pub(super) fn geom_adaptor_nb_intervals_cn(
    c: &rcad_kernel::geom::Curve3,
    first: f64,
    last: f64,
) -> i32 {
    use rcad_kernel::geom::Curve3;
    match c {
        Curve3::BSpline(bs) => {
            // (!IsPeriodic && S <= Continuity()) || S == GeomAbs_C0 -> 1
            if !bs.is_periodic && geom_adaptor_curve_continuity(c) == rcad_kernel::math::GeomAbsShape::CN {
                return 1;
            }
            let a_cont = bs.degree as i32; // case GeomAbs_CN: aCont = aDegree
            let mut intervals = 0i32;
            let mut i = 0usize;
            let n = bs.knots.len();
            while i < n {
                let mut j = i;
                while j < n && (bs.knots[j] - bs.knots[i]).abs() < 1.0e-12 {
                    j += 1;
                }
                let mult = (j - i) as i32;
                let k = bs.knots[i];
                if k > first + 1.0e-12 && k < last - 1.0e-12 && mult > a_cont {
                    intervals += 1;
                }
                i = j;
            }
            intervals + 1
        }
        _ => 1,
    }
}

// ---------------------------------------------------------------------------
// BRepFill_Sweep::CorrectApproxParameters (cxx L2043-2074)
// ---------------------------------------------------------------------------

impl BRepFillSweep {
    /// OCCT CorrectApproxParameters (cxx L2043-2074).
    pub fn correct_approx_parameters(&mut self, brep: &BRep) -> bool {
        // TopoDS_Wire thePath = myLoc->Wire();
        let the_path = self.my_loc.borrow().base().wire();
        // GeomAbs_Shape NewCont = myContinuity;
        let mut new_cont = self.my_continuity;
        let mut new_segmax = self.my_segmax;

        // TopoDS_Iterator iter(thePath);
        let iter = super::brep_fill_sweep_b::topods_iterator(&the_path);
        for it in iter {
            // TopoDS_Edge anEdge = TopoDS::Edge(iter.Value());
            // BRepAdaptor_Curve aBAcurve(anEdge);
            let ed = brep.edge(it.clone());
            let a_ba_curve = ed
                .curve
                .clone()
                .expect("CorrectApproxParameters: edge without 3d curve");
            // GeomAbs_Shape aContinuity = aBAcurve.Continuity();
            let a_continuity = geom_adaptor_curve_continuity(&a_ba_curve);
            // int aNbInterv = aBAcurve.NbIntervals(GeomAbs_CN);
            let a_nb_interv =
                geom_adaptor_nb_intervals_cn(&a_ba_curve, ed.range[0], ed.range[1]);
            if a_continuity < new_cont {
                new_cont = a_continuity;
            }
            if a_nb_interv > new_segmax {
                new_segmax = a_nb_interv;
            }
        }

        let mut corrected = false;
        if new_cont != self.my_continuity || new_segmax != self.my_segmax {
            corrected = true;
        }
        self.my_continuity = new_cont;
        self.my_segmax = new_segmax;
        corrected
    }
}

// ---------------------------------------------------------------------------
// BRepFill_Sweep::BuildShell (cxx L2213-3207)
// ---------------------------------------------------------------------------

impl BRepFillSweep {
    /// OCCT BuildShell (cxx L2213-3207).
    #[allow(clippy::too_many_arguments)]
    pub fn build_shell(
        &mut self,
        brep: &mut BRep,
        _transition: crate::brep_fill::brep_fill_pipe_shell_b::BRepFillTransitionStyle,
        ifirst: i32,
        ilast: i32,
        reversed_edges: &mut HashSet<ShapeKey>,
        tapes: &mut ShapeToArray2Map,
        rails: &mut ShapeToArray2Map,
        extend_first: f64,
        extend_last: f64,
    ) -> bool {
        let mut b = BRepBuilder::new();
        // int NbPath = ILast - IFirst;
        let nb_path = ilast - ifirst;
        // int NbLaw = mySec->NbLaw();
        let nb_law = self.my_sec.borrow().base().nb_law();
        let mut hasdegen = false;
        let const_section = self.my_sec.borrow().is_constant();
        let uclose = self.my_sec.borrow().base().is_uclosed();
        let global_vclose = self.my_loc.borrow().base().is_closed(brep)
            && self.my_loc.borrow().base().is_g1(brep, 0, self.my_tol3d, 1.0e-4) >= 0;
        let vclose = global_vclose
            && self.my_sec.borrow().base().is_vclosed()
            && nb_path == self.my_loc.borrow().base().nb_law();
        self.error = 0.0;

        // (1) Construction of all surfaces

        // (1.1) Construction of Tables
        let mut exch_uv: Array2<bool> = Array2::new(1, nb_law as usize, 1, nb_path as usize, false);
        let mut u_reverse: Array2<bool> = Array2::new(1, nb_law as usize, 1, nb_path as usize, false);
        let mut degenerated: Array2<bool> =
            Array2::new(1, nb_law as usize, 1, nb_path as usize, false);
        // Degenerated.Init(false);
        // No VReverse for the moment...
        let mut tab_err: Array2<f64> = Array2::new(1, nb_law as usize, 1, nb_path as usize, 0.0);
        let mut tab_s: Array2<SurfaceShape> =
            Array2::new(1, nb_law as usize, 1, nb_path as usize, SurfaceShape::None);

        let mut u_edge: Array2<Shape> =
            Array2::new(1, (nb_law + 1) as usize, 1, nb_path as usize, Shape::null());
        let mut v_edge: Array2<Shape> =
            Array2::new(1, nb_law as usize, 1, (nb_path + 1) as usize, Shape::null());
        let mut vertex: Array2<Shape> =
            Array2::new(1, (nb_law + 1) as usize, 1, (nb_path + 1) as usize, Shape::null());

        // TopoDS_Vertex VNULL; VNULL.Nullify(); Vertex.Init(VNULL);

        let mut sec_vertex: Vec<Shape> = vec![Shape::null(); (nb_law + 1) as usize];
        let mut v_error: Vec<f64> = vec![0.0; (nb_law + 1) as usize];
        let mut vi: Vec<f64> = vec![0.0; (nb_path + 1) as usize];

        // Initialization of management of parametric intervals
        // (Case of evolutionary sections)
        // myLoc->CurvilinearBounds(myLoc->NbLaw(), SecDom, Length);
        let mut length = 0.0f64;
        let mut sec_dom = 0.0f64;
        {
            let nb_law_loc = self.my_loc.borrow().base().nb_law();
            self.my_loc
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(nb_law_loc, &mut sec_dom, &mut length);
        }
        // mySec->Law(1)->GetDomain(SecDeb, SecDom);
        let mut sec_deb = 0.0f64;
        {
            let mut sd = 0.0f64;
            self.my_sec
                .borrow()
                .base()
                .law(1)
                .borrow()
                .get_domain(&mut sec_deb, &mut sd);
            sec_dom = sd;
        }
        // SecDom -= SecDeb;
        sec_dom -= sec_deb;
        if ifirst > 1 {
            let (mut lf, mut ll) = (0.0f64, 0.0f64);
            self.my_loc
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(ifirst - 1, &mut lf, &mut ll);
            vi[0] = sec_deb + (ll / length) * sec_dom;
        } else {
            vi[0] = sec_deb;
        }

        // Error a priori on vertices
        if const_section {
            for isec in 1..=(nb_law + 1) {
                v_error[isec as usize - 1] =
                    self.my_sec.borrow().vertex_tol(brep, isec - 1, 0.0);
                sec_vertex[isec as usize - 1] = self.my_sec.borrow().vertex(brep, isec, 0.0);
            }
        }

        // (1.2) Calculate surfaces
        let mut ipath = 1i32;
        let mut ipath_l = ifirst;
        while ipath <= nb_path {
            // GeomFill_Sweep Sweep(myLoc->Law(IPath), KPart);
            let loc_law = self.my_loc.borrow().base().law(ipath_l);
            let mut sweep = crate::geomalgo::geomfill::sweep::Sweep::new(loc_law, self.kpart);
            sweep.set_tolerance_sweep(self.my_tol3d, self.my_bound_tol, self.my_tol2d, self.my_tol_angular);
            sweep.set_force_approx_c1(self.my_force_approx_c1);

            // Case of evolutionary section, definition of parametric
            // correspondence
            if !const_section {
                let (mut lf, mut ll) = (0.0f64, 0.0f64);
                let (mut lf2, mut ll2) = (0.0f64, 0.0f64);
                {
                    let law = self.my_loc.borrow().base().law(ipath_l);
                    law.borrow().get_domain(&mut lf, &mut ll);
                }
                self.my_loc
                    .borrow_mut()
                    .base_mut()
                    .curvilinear_bounds(ipath_l, &mut lf2, &mut ll2);
                vi[ipath as usize] = sec_deb + (ll2 / length) * sec_dom;
                sweep.set_domain(lf, ll, vi[(ipath - 1) as usize], vi[ipath as usize]);
            } else {
                // section is constant
                let (mut lf, mut ll) = (0.0f64, 0.0f64);
                let (mut lf2, mut ll2) = (0.0f64, 0.0f64);
                {
                    let law = self.my_loc.borrow().base().law(ipath_l);
                    law.borrow().get_domain(&mut lf, &mut ll);
                }
                self.my_loc
                    .borrow_mut()
                    .base_mut()
                    .curvilinear_bounds(ipath_l, &mut lf2, &mut ll2);
                vi[ipath as usize] = sec_deb + (ll2 / length) * sec_dom;
            }

            for isec in 1..=nb_law {
                // Sweep.Build(mySec->Law(isec), myApproxStyle, myContinuity,
                //             myDegmax, mySegmax);
                let sec_law = self.my_sec.borrow().base().law(isec);
                sweep.build(sec_law, self.my_approx_style, self.my_continuity, self.my_degmax, self.my_segmax);
                if !sweep.is_done() {
                    return false;
                }
                tab_s.set(isec as usize, ipath as usize, SurfaceShape::Some(sweep.surface().expect("null surface").clone()));
                tab_err.set(isec as usize, ipath as usize, sweep.error_on_surface());
                exch_uv.set(isec as usize, ipath as usize, sweep.exchange_uv());
                u_reverse.set(isec as usize, ipath as usize, sweep.u_reversed());
                if sweep.error_on_surface() > self.error {
                    self.error = sweep.error_on_surface();
                }

                if (ipath == 1) && (extend_first > 0.0) {
                    // GeomLib::ExtendSurfByLength(BndS, ExtendFirst, 1,
                    //                             Sweep.ExchangeUV(), false);
                    let extended = crate::geomalgo::geom_lib_same_range::extend_surf_by_length(
                        tab_s.get(isec as usize, ipath as usize).surface(),
                        extend_first,
                        1,
                        sweep.exchange_uv(),
                        false,
                    );
                    tab_s.set(isec as usize, ipath as usize, SurfaceShape::Some(extended));
                }
                if (ipath == nb_path) && (extend_last > 0.0) {
                    // GeomLib::ExtendSurfByLength(BndS, ExtendLast, 1,
                    //                             Sweep.ExchangeUV(), true);
                    let extended = crate::geomalgo::geom_lib_same_range::extend_surf_by_length(
                        tab_s.get(isec as usize, ipath as usize).surface(),
                        extend_last,
                        1,
                        sweep.exchange_uv(),
                        true,
                    );
                    tab_s.set(isec as usize, ipath as usize, SurfaceShape::Some(extended));
                }
            }
            ipath += 1;
            ipath_l += 1;
        }

        // (2) Construction of Edges
        // Correct <FirstShape> and <LastShape>: reverse modified edges
        super::brep_fill_sweep::reverse_modified_edges(
            brep,
            &mut self.first_shape,
            reversed_edges,
        );
        super::brep_fill_sweep::reverse_modified_edges(
            brep,
            &mut self.last_shape,
            reversed_edges,
        );

        // (2.0) return preexisting Edges and vertices
        let mut is_built: Vec<bool> = vec![false; nb_law as usize];
        let mut start_edges: Vec<Shape> = vec![Shape::null(); nb_law as usize];
        if !self.first_shape.is_null() && (ifirst == 1) {
            self.my_sec.borrow_mut().base_mut().init(brep, &self.first_shape.clone());
            for isec in 1..=nb_law {
                let e = self.my_sec.borrow_mut().base_mut().current_edge(brep);
                let (vfirst, vlast) =
                    crate::brep_fill::brep_fill_pipe::top_exp_vertices(&e);
                v_edge.set(isec as usize, 1, e.clone());
                let v_end = if e.orientation == Orientation::Reversed {
                    vfirst.clone() // TopExp::FirstVertex(E);
                } else {
                    vlast.clone() // TopExp::LastVertex(E);
                };
                vertex.set((isec + 1) as usize, 1, v_end.clone());
                let mut v_end_mut = v_end;
                self.update_vertex(
                    brep,
                    ifirst - 1,
                    isec + 1,
                    tab_err.get(isec as usize, 1),
                    vi[0],
                    &mut v_end_mut,
                );
                vertex.set((isec + 1) as usize, 1, v_end_mut);

                start_edges[isec as usize - 1] = e.clone();
                if tapes.contains_key(&ShapeKey(e.ptr_id())) {
                    is_built[isec as usize - 1] = true;

                    // Initialize VEdge, UEdge, Vertex and myFaces
                    let tape = tapes.get(&ShapeKey(e.ptr_id())).expect("Tapes(E)").clone();
                    for j in 1..=(nb_path + 1) {
                        let mut ve = tape[0][j as usize - 1].clone();
                        // VEdge(isec, j).Reverse(); direction of round is reversed
                        ve.orientation = match ve.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            o => o,
                        };
                        v_edge.set(isec as usize, j as usize, ve);
                    }
                    let ifirst_i = isec + 1;
                    let ilast_i = isec; // direction of round is reversed
                    for j in 1..=nb_path {
                        u_edge.set(ifirst_i as usize, j as usize, tape[1][j as usize - 1].clone());
                    }
                    for j in 1..=nb_path {
                        u_edge.set(ilast_i as usize, j as usize, tape[2][j as usize - 1].clone());
                    }
                    for j in 1..=(nb_path + 1) {
                        vertex.set(ifirst_i as usize, j as usize, tape[3][j as usize - 1].clone());
                    }
                    for j in 1..=(nb_path + 1) {
                        vertex.set(ilast_i as usize, j as usize, tape[4][j as usize - 1].clone());
                    }
                    for j in 1..=nb_path {
                        self.my_faces.as_mut().expect("myFaces")[isec as usize - 1]
                            [(j - 1) as usize] = tape[5][j as usize - 1].clone();
                    }

                    if uclose && isec == 1 {
                        for j in 1..=nb_path {
                            let ue = u_edge.get(1, j as usize);
                            u_edge.set((nb_law + 1) as usize, j as usize, ue);
                        }
                        for j in 1..=(nb_path + 1) {
                            let vv = vertex.get(1, j as usize);
                            vertex.set((nb_law + 1) as usize, j as usize, vv);
                        }
                    }
                    if uclose && isec == nb_law {
                        for j in 1..=nb_path {
                            let ue = u_edge.get((nb_law + 1) as usize, j as usize);
                            u_edge.set(1, j as usize, ue);
                        }
                        for j in 1..=(nb_path + 1) {
                            let vv = vertex.get((nb_law + 1) as usize, j as usize);
                            vertex.set(1, j as usize, vv);
                        }
                    }
                } else {
                    // occ::handle<NCollection_HArray2<TopoDS_Shape>> EmptyArray =
                    //   new NCollection_HArray2<TopoDS_Shape>(1, 6, 1, NbPath + 1);
                    tapes.insert(
                        ShapeKey(e.ptr_id()),
                        vec![vec![Shape::null(); (nb_path + 1) as usize]; 6],
                    );
                    if let Some(vfirst_shape) = non_null(&vfirst) {
                        if rails.contains_key(&ShapeKey(vfirst_shape.ptr_id())) {
                            let rail =
                                rails.get(&ShapeKey(vfirst_shape.ptr_id())).expect("Rails").clone();
                            let ind = if e.orientation == Orientation::Reversed {
                                isec + 1
                            } else {
                                isec
                            };
                            for j in 1..=nb_path {
                                u_edge.set(ind as usize, j as usize, rail[0][j as usize - 1].clone());
                            }
                            for j in 1..=(nb_path + 1) {
                                vertex.set(ind as usize, j as usize, rail[1][j as usize - 1].clone());
                            }
                        }
                    }
                    if let Some(vlast_shape) = non_null(&vlast) {
                        if rails.contains_key(&ShapeKey(vlast_shape.ptr_id())) {
                            let rail =
                                rails.get(&ShapeKey(vlast_shape.ptr_id())).expect("Rails").clone();
                            let ind = if e.orientation == Orientation::Forward {
                                isec + 1
                            } else {
                                isec
                            };
                            for j in 1..=nb_path {
                                u_edge.set(ind as usize, j as usize, rail[0][j as usize - 1].clone());
                            }
                            for j in 1..=(nb_path + 1) {
                                vertex.set(ind as usize, j as usize, rail[1][j as usize - 1].clone());
                            }
                        }
                    }
                }
            }

            if v_edge.get(1, 1).orientation == Orientation::Reversed {
                // Vertex(1, 1) = TopExp::LastVertex(TopoDS::Edge(VEdge(1, 1)));
                let ve = v_edge.get(1, 1);
                let ed = brep.edge(ve.clone());
                vertex.set(1, 1, ed.last.clone());
            } else {
                // Vertex(1, 1) = TopExp::FirstVertex(TopoDS::Edge(VEdge(1, 1)));
                let ve = v_edge.get(1, 1);
                let ed = brep.edge(ve.clone());
                vertex.set(1, 1, ed.first.clone());
            }
            let mut vv = vertex.get(1, 1);
            self.update_vertex(brep, ifirst - 1, 1, tab_err.get(1, 1), vi[0], &mut vv);
            vertex.set(1, 1, vv);
        }

        // double u, v, aux; bool ureverse;
        for isec in 1..=(nb_law + 1) {
            // Return data
            let (surf, ureverse, exuv) = if isec > nb_law {
                (
                    tab_s.get(nb_law as usize, 1),
                    u_reverse.get(nb_law as usize, 1),
                    exch_uv.get(nb_law as usize, 1),
                )
            } else {
                (
                    tab_s.get(isec as usize, 1),
                    u_reverse.get(isec as usize, 1),
                    exch_uv.get(isec as usize, 1),
                )
            };
            let surf = surf.surface().clone();
            let (mut u_first, mut u_last, mut v_first, mut v_last) = surface_bounds(&surf);

            // Choice of parameters
            if ureverse {
                if exuv {
                    let aux = v_first;
                    v_first = v_last;
                    v_last = aux;
                } else {
                    let aux = u_first;
                    u_first = u_last;
                    u_last = aux;
                }
            }
            let (u, v) = if isec != nb_law + 1 {
                (u_first, v_first)
            } else if exuv {
                (u_first, v_last)
            } else {
                (u_last, v_first)
            };

            // construction of vertices
            if vertex.get(isec as usize, 1).is_null() {
                // B.MakeVertex(TopoDS::Vertex(Vertex(isec, 1)), S->Value(u, v),
                //              mySec->VertexTol(isec - 1, Vi(1)));
                let v = b.add_vertex(
                    brep,
                    surf.point_at(u, v),
                    self.my_sec.borrow().vertex_tol(brep, isec - 1, vi[0]),
                );
                vertex.set(isec as usize, 1, v);
            } else {
                // TopLoc_Location Identity; Vertex(isec, 1).Location(Identity);
                let mut vv = vertex.get(isec as usize, 1);
                vv.location = 0;
                // B.UpdateVertex(..., S->Value(u, v), mySec->VertexTol(...));
                let tol = self.my_sec.borrow().vertex_tol(brep, isec - 1, vi[0]);
                b.update_vertex_point(brep, vv.clone(), surf.point_at(u, v), tol);
                vertex.set(isec as usize, 1, vv);
            }
        } // end of for (isec=1; isec<=NbLaw+1; isec++)

        if !self.last_shape.is_null() && (ilast == self.my_loc.borrow().base().nb_law() + 1) {
            self.my_sec.borrow_mut().base_mut().init(brep, &self.last_shape.clone());
            for isec in 1..=nb_law {
                let e = self.my_sec.borrow_mut().base_mut().current_edge(brep);
                if v_edge.get(isec as usize, (nb_path + 1) as usize).is_null() {
                    v_edge.set(isec as usize, (nb_path + 1) as usize, e.clone());
                }

                if vertex.get((isec + 1) as usize, (nb_path + 1) as usize).is_null() {
                    let ve = v_edge.get(isec as usize, (nb_path + 1) as usize);
                    let ed = brep.edge(ve.clone());
                    let v_end = if ve.orientation == Orientation::Reversed {
                        ed.first.clone()
                    } else {
                        ed.last.clone()
                    };
                    vertex.set((isec + 1) as usize, (nb_path + 1) as usize, v_end);
                }
                let mut vv = vertex.get((isec + 1) as usize, (nb_path + 1) as usize);
                self.update_vertex(
                    brep,
                    ilast - 1,
                    isec + 1,
                    tab_err.get(isec as usize, nb_path as usize),
                    vi[nb_path as usize],
                    &mut vv,
                );
                vertex.set((isec + 1) as usize, (nb_path + 1) as usize, vv);
            }

            if vertex.get(1, (nb_path + 1) as usize).is_null() {
                let ve = v_edge.get(1, (nb_path + 1) as usize);
                let ed = brep.edge(ve.clone());
                let v_end = if ve.orientation == Orientation::Reversed {
                    ed.last.clone()
                } else {
                    ed.first.clone()
                };
                vertex.set(1, (nb_path + 1) as usize, v_end);
            }
            let mut vv = vertex.get(1, (nb_path + 1) as usize);
            self.update_vertex(
                brep,
                ilast - 1,
                1,
                tab_err.get(1, nb_path as usize),
                vi[nb_path as usize],
                &mut vv,
            );
            vertex.set(1, (nb_path + 1) as usize, vv);
        }

        for isec in 1..=(nb_law + 1) {
            // Return data
            let (surf, ureverse, exuv) = if isec > nb_law {
                (
                    tab_s.get(nb_law as usize, nb_path as usize),
                    u_reverse.get(nb_law as usize, nb_path as usize),
                    exch_uv.get(nb_law as usize, nb_path as usize),
                )
            } else {
                (
                    tab_s.get(isec as usize, nb_path as usize),
                    u_reverse.get(isec as usize, nb_path as usize),
                    exch_uv.get(isec as usize, nb_path as usize),
                )
            };
            let surf = surf.surface().clone();
            let (mut u_first, mut u_last, mut v_first, mut v_last) = surface_bounds(&surf);

            // Choice of parametres
            if ureverse {
                if exuv {
                    let aux = v_first;
                    v_first = v_last;
                    v_last = aux;
                } else {
                    let aux = u_first;
                    u_first = u_last;
                    u_last = aux;
                }
            }
            let (u, v) = if isec == nb_law + 1 {
                (u_last, v_last)
            } else if exuv {
                (u_last, v_first)
            } else {
                (u_first, v_last)
            };

            // construction of vertex
            if vertex.get(isec as usize, (nb_path + 1) as usize).is_null() {
                let v = b.add_vertex(
                    brep,
                    surf.point_at(u, v),
                    self.my_sec
                        .borrow()
                        .vertex_tol(brep, isec - 1, vi[nb_path as usize]),
                );
                vertex.set(isec as usize, (nb_path + 1) as usize, v);
            } else {
                let mut vv = vertex.get(isec as usize, (nb_path + 1) as usize);
                vv.location = 0;
                let tol = self
                    .my_sec
                    .borrow()
                    .vertex_tol(brep, isec - 1, vi[nb_path as usize]);
                b.update_vertex_point(brep, vv.clone(), surf.point_at(u, v), tol);
                vertex.set(isec as usize, (nb_path + 1) as usize, vv);
            }
        } // end of for (isec=1; isec<=NbLaw+1; isec++)

        // ---------- Creation of Vertex and edge ------------
        let mut ipath = 1i32;
        let mut ipath_l = ifirst;
        while ipath <= nb_path {
            for isec in 1..=nb_law {
                if is_built[isec as usize - 1] {
                    continue;
                }

                let surf = tab_s.get(isec as usize, ipath as usize).surface().clone();
                let mut exuv = exch_uv.get(isec as usize, ipath as usize);
                let (mut u_first, mut u_last, mut v_first, mut v_last) = surface_bounds(&surf);
                if u_reverse.get(isec as usize, ipath as usize) {
                    if exuv {
                        let aux2 = v_first;
                        v_first = v_last;
                        v_last = aux2;
                    } else {
                        let aux2 = u_first;
                        u_first = u_last;
                        u_last = aux2;
                    }
                }

                // (2.1) Construction of new vertices
                if isec == 1 {
                    if ipath == 1 && vertex.get(1, 1).is_null() {
                        // All first
                        if const_section {
                            let sec_v = sec_vertex[0].clone();
                            let mut target = vertex.get(1, 1);
                            self.my_loc.borrow().base().perform_vertex(
                                brep,
                                ipath_l - 1,
                                &sec_v,
                                v_error[0],
                                &mut target,
                                0,
                            );
                            vertex.set(1, 1, target);
                        } else {
                            let sec_v = self.my_sec.borrow().vertex(brep, 1, vi[0]);
                            let sec_tol = self.my_sec.borrow().vertex_tol(brep, 0, vi[0]);
                            let mut target = vertex.get(1, 1);
                            self.my_loc.borrow().base().perform_vertex(
                                brep,
                                ipath_l - 1,
                                &sec_v,
                                sec_tol,
                                &mut target,
                                0,
                            );
                            vertex.set(1, 1, target);
                        }
                    }
                    // the first and the next column
                    if vclose && (ipath == nb_path) {
                        let vv = vertex.get(1, 1);
                        vertex.set(1, (ipath + 1) as usize, vv);
                    } else if vertex.get(1, (ipath + 1) as usize).is_null() {
                        if const_section {
                            let sec_v = sec_vertex[0].clone();
                            let mut target = vertex.get(1, (ipath + 1) as usize);
                            self.my_loc.borrow().base().perform_vertex(
                                brep,
                                ipath_l,
                                &sec_v,
                                tab_err.get(1, ipath as usize) + v_error[0],
                                &mut target,
                                0,
                            );
                            vertex.set(1, (ipath + 1) as usize, target);
                        } else {
                            let sec_v = self.my_sec.borrow().vertex(brep, 1, vi[ipath as usize]);
                            let sec_tol =
                                self.my_sec.borrow().vertex_tol(brep, 0, vi[ipath as usize]);
                            let mut target = vertex.get(1, (ipath + 1) as usize);
                            self.my_loc.borrow().base().perform_vertex(
                                brep,
                                ipath_l,
                                &sec_v,
                                tab_err.get(1, ipath as usize) + sec_tol,
                                &mut target,
                                0,
                            );
                            vertex.set(1, (ipath + 1) as usize, target);
                        }
                        let mut a = vertex.get(1, ipath as usize);
                        let mut bb = vertex.get(1, (ipath + 1) as usize);
                        if self.merge_vertex(brep, &a, &mut bb) {
                            // UEdge(1, ipath) = NullEdge(Vertex(1, ipath));
                            let ne = null_edge(brep, &mut a);
                            u_edge.set(1, ipath as usize, ne);
                        }
                        vertex.set(1, ipath as usize, a);
                        vertex.set(1, (ipath + 1) as usize, bb);
                    }
                }

                if ipath == 1 {
                    if uclose && (isec == nb_law) {
                        let vv = vertex.get(1, 1);
                        vertex.set((isec + 1) as usize, 1, vv);
                    } else if vertex.get((isec + 1) as usize, 1).is_null() {
                        if const_section {
                            let sec_v = sec_vertex[isec as usize].clone();
                            let mut target = vertex.get((isec + 1) as usize, 1);
                            self.my_loc.borrow().base().perform_vertex(
                                brep,
                                ipath_l - 1,
                                &sec_v,
                                tab_err.get(isec as usize, 1) + v_error[isec as usize],
                                &mut target,
                                0,
                            );
                            vertex.set((isec + 1) as usize, 1, target);
                        } else {
                            let sec_v =
                                self.my_sec.borrow().vertex(brep, isec + 1, vi[0]);
                            let sec_tol =
                                self.my_sec.borrow().vertex_tol(brep, isec, vi[0]);
                            let mut target = vertex.get((isec + 1) as usize, 1);
                            self.my_loc.borrow().base().perform_vertex(
                                brep,
                                ipath_l - 1,
                                &sec_v,
                                tab_err.get(isec as usize, 1) + sec_tol,
                                &mut target,
                                0,
                            );
                            vertex.set((isec + 1) as usize, 1, target);
                        }

                        let mut a = vertex.get(isec as usize, 1);
                        let mut bbv = vertex.get((isec + 1) as usize, 1);
                        if self.merge_vertex(brep, &a, &mut bbv) {
                            // VEdge(isec, 1) = NullEdge(Vertex(isec, 1));
                            let ne = null_edge(brep, &mut a);
                            v_edge.set(isec as usize, 1, ne);
                        }
                        vertex.set(isec as usize, 1, a);
                        vertex.set((isec + 1) as usize, 1, bbv);
                    }
                }

                if uclose && (isec == nb_law) {
                    let vv = vertex.get(1, (ipath + 1) as usize);
                    vertex.set((isec + 1) as usize, (ipath + 1) as usize, vv);
                } else if vclose && (ipath == nb_path) {
                    let vv = vertex.get((isec + 1) as usize, 1);
                    vertex.set((isec + 1) as usize, (ipath + 1) as usize, vv);
                } else if vertex.get((isec + 1) as usize, (ipath + 1) as usize).is_null() {
                    if const_section {
                        let sec_v = sec_vertex[isec as usize].clone();
                        let mut target = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                        self.my_loc.borrow().base().perform_vertex(
                            brep,
                            ipath_l,
                            &sec_v,
                            tab_err.get(isec as usize, ipath as usize) + v_error[isec as usize],
                            &mut target,
                            0,
                        );
                        vertex.set((isec + 1) as usize, (ipath + 1) as usize, target);
                    } else {
                        let sec_v =
                            self.my_sec.borrow().vertex(brep, isec + 1, vi[ipath as usize]);
                        let sec_tol =
                            self.my_sec.borrow().vertex_tol(brep, isec, vi[ipath as usize]);
                        let mut target = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                        self.my_loc.borrow().base().perform_vertex(
                            brep,
                            ipath_l,
                            &sec_v,
                            tab_err.get(isec as usize, ipath as usize) + sec_tol,
                            &mut target,
                            0,
                        );
                        vertex.set((isec + 1) as usize, (ipath + 1) as usize, target);
                    }
                }

                // Singular cases
                let mut a_singv = vertex.get(isec as usize, (ipath + 1) as usize);
                let mut b_singv = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                let singv = self.merge_vertex(brep, &a_singv, &mut b_singv);
                vertex.set((isec + 1) as usize, (ipath + 1) as usize, b_singv);
                let mut a_singu = vertex.get((isec + 1) as usize, ipath as usize);
                let mut b_singu = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                let singu = self.merge_vertex(brep, &a_singu, &mut b_singu);
                vertex.set((isec + 1) as usize, (ipath + 1) as usize, b_singu);

                if singu || singv {
                    let d = is_degen(
                        tab_s.get(isec as usize, ipath as usize).surface(),
                        self.my_tol3d.max(tab_err.get(isec as usize, ipath as usize)),
                    );
                    degenerated.set(isec as usize, ipath as usize, d);
                }
                if degenerated.get(isec as usize, ipath as usize) {
                    hasdegen = true;
                    // Particular construction of edges
                    if u_edge.get((isec + 1) as usize, ipath as usize).is_null() {
                        if singu {
                            // Degenerated edge
                            let mut vv = vertex.get((isec + 1) as usize, ipath as usize);
                            let ne = null_edge(brep, &mut vv);
                            u_edge.set((isec + 1) as usize, ipath as usize, ne);
                        } else {
                            // Copy the previous edge
                            let ue = u_edge.get(isec as usize, ipath as usize);
                            u_edge.set((isec + 1) as usize, ipath as usize, ue);
                        }
                    }
                    if v_edge.get(isec as usize, (ipath + 1) as usize).is_null() {
                        if singv {
                            // Degenerated Edge
                            let mut vv = vertex.get(isec as usize, (ipath + 1) as usize);
                            let ne = null_edge(brep, &mut vv);
                            v_edge.set(isec as usize, (ipath + 1) as usize, ne);
                        } else {
                            // Copy the previous edge
                            let ve = v_edge.get(isec as usize, ipath as usize);
                            v_edge.set(isec as usize, (ipath + 1) as usize, ve);
                        }
                    }
                } else {
                    // Construction of edges by isos
                    if exuv {
                        let uv = u_first;
                        u_first = v_first;
                        v_first = uv;
                        let uv = u_last;
                        u_last = v_last;
                        v_last = uv;
                    }

                    // (2.2) Iso-u
                    if isec == 1 && u_edge.get(1, ipath as usize).is_null() {
                        let vv1 = vertex.get(1, ipath as usize);
                        let vv2 = vertex.get(1, (ipath + 1) as usize);
                        if !vv1.is_same(&vv2) {
                            let p1 = brep.vertex(vv1.clone()).point;
                            let p2 = brep.vertex(vv2.clone()).point;
                            if p1.distance(p2) <= self.my_tol3d {
                                vertex.set(1, (ipath + 1) as usize, vv1.clone());
                            }
                        }
                        let ue = super::brep_fill_sweep_b::build_edge_iso(
                            brep,
                            &surf,
                            !exuv,
                            u_first,
                            &vv1,
                            &vv2,
                            self.my_tol3d,
                        );
                        u_edge.set(1, ipath as usize, ue);
                    } else {
                        if u_edge.get(isec as usize, ipath as usize).is_null() {
                            // sweep failed
                            return false;
                        }
                        let ue = u_edge.get(isec as usize, ipath as usize);
                        super::brep_fill_sweep_b::update_edge_iso(brep, &ue, &surf, !exuv, u_first);
                    }

                    if uclose && (isec == nb_law) {
                        if u_edge.get(1, ipath as usize).is_null() {
                            // degenerated case
                            let vv1 = vertex.get((isec + 1) as usize, ipath as usize);
                            let vv2 = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                            let ue = super::brep_fill_sweep_b::build_edge_iso(
                                brep,
                                &surf,
                                !exuv,
                                u_last,
                                &vv1,
                                &vv2,
                                self.my_tol3d,
                            );
                            u_edge.set((isec + 1) as usize, ipath as usize, ue);
                        } else {
                            let ue = u_edge.get(1, ipath as usize);
                            super::brep_fill_sweep_b::update_edge_iso(brep, &ue, &surf, !exuv, u_last);
                            u_edge.set((isec + 1) as usize, ipath as usize, ue);
                        }
                    } else {
                        if u_edge.get((isec + 1) as usize, ipath as usize).is_null() {
                            let vv1 = vertex.get((isec + 1) as usize, ipath as usize);
                            let vv2 = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                            let ue = super::brep_fill_sweep_b::build_edge_iso(
                                brep,
                                &surf,
                                !exuv,
                                u_last,
                                &vv1,
                                &vv2,
                                self.my_tol3d,
                            );
                            u_edge.set((isec + 1) as usize, ipath as usize, ue);
                        } else {
                            let ue = u_edge.get((isec + 1) as usize, ipath as usize);
                            super::brep_fill_sweep_b::update_edge_iso(brep, &ue, &surf, !exuv, u_last);
                        }
                    }

                    // (2.3) Iso-v
                    if ipath == 1 {
                        let vv1 = vertex.get(isec as usize, 1);
                        let vv2 = vertex.get((isec + 1) as usize, 1);
                        // TopoDS_Edge aNewFirstEdge = BuildEdge(...);
                        let a_new_first_edge = super::brep_fill_sweep_b::build_edge_iso(
                            brep,
                            &surf,
                            exuv,
                            v_first,
                            &vv1,
                            &vv2,
                            self.my_tol3d,
                        );
                        if v_edge.get(isec as usize, ipath as usize).is_null() {
                            v_edge.set(isec as usize, ipath as usize, a_new_first_edge);
                        } else {
                            // rebuild first edge
                            let mut old_edge = v_edge.get(isec as usize, ipath as usize);
                            self.rebuild_top_or_bottom_edge(
                                brep,
                                &a_new_first_edge,
                                &mut old_edge,
                                reversed_edges,
                            );
                            v_edge.set(isec as usize, ipath as usize, old_edge.clone());
                            if reversed_edges.contains(&ShapeKey(old_edge.ptr_id())) {
                                let mut ve = old_edge.clone();
                                super::brep_fill_sweep::reverse_edge_in_first_or_last_wire(
                                    brep,
                                    &mut self.first_shape,
                                    &ve,
                                );
                                let mut se = start_edges[isec as usize - 1].clone();
                                se.orientation = match se.orientation {
                                    Orientation::Forward => Orientation::Reversed,
                                    Orientation::Reversed => Orientation::Forward,
                                    o => o,
                                };
                                start_edges[isec as usize - 1] = se;
                                let _ = ve;
                            }
                        }
                    } else {
                        let ve = v_edge.get(isec as usize, ipath as usize);
                        super::brep_fill_sweep_b::update_edge_iso(brep, &ve, &surf, exuv, v_first);
                    }

                    if vclose && (ipath == nb_path) {
                        if v_edge.get(isec as usize, 1).is_null() {
                            // degenerated case
                            let vv1 = vertex.get(isec as usize, (ipath + 1) as usize);
                            let vv2 = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                            let ve = super::brep_fill_sweep_b::build_edge_iso(
                                brep,
                                &surf,
                                exuv,
                                v_last,
                                &vv1,
                                &vv2,
                                self.my_tol3d,
                            );
                            v_edge.set(isec as usize, (ipath + 1) as usize, ve);
                        } else {
                            let ve1 = v_edge.get(isec as usize, 1);
                            super::brep_fill_sweep_b::update_edge_iso(brep, &ve1, &surf, exuv, v_last);
                            v_edge.set(isec as usize, (ipath + 1) as usize, ve1);
                        }
                    } else if v_edge.get(isec as usize, (ipath + 1) as usize).is_null() {
                        let vv1 = vertex.get(isec as usize, (ipath + 1) as usize);
                        let vv2 = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                        let ve = super::brep_fill_sweep_b::build_edge_iso(
                            brep,
                            &surf,
                            exuv,
                            v_last,
                            &vv1,
                            &vv2,
                            self.my_tol3d,
                        );
                        v_edge.set(isec as usize, (ipath + 1) as usize, ve);
                    } else {
                        if ipath != nb_path
                            || vclose
                            || (global_vclose && ilast == self.my_loc.borrow().base().nb_law() + 1)
                        {
                            let ve = v_edge.get(isec as usize, (ipath + 1) as usize);
                            super::brep_fill_sweep_b::update_edge_iso(brep, &ve, &surf, exuv, v_last);
                        } else {
                            // ipath == NbPath && !vclose => rebuild last edge
                            let vv1 = vertex.get(isec as usize, (ipath + 1) as usize);
                            let vv2 = vertex.get((isec + 1) as usize, (ipath + 1) as usize);
                            let a_new_last_edge = super::brep_fill_sweep_b::build_edge_iso(
                                brep,
                                &surf,
                                exuv,
                                v_last,
                                &vv1,
                                &vv2,
                                self.my_tol3d,
                            );
                            let mut old_edge = v_edge.get(isec as usize, (ipath + 1) as usize);
                            self.rebuild_top_or_bottom_edge(
                                brep,
                                &a_new_last_edge,
                                &mut old_edge,
                                reversed_edges,
                            );
                            v_edge.set(isec as usize, (ipath + 1) as usize, old_edge.clone());
                            if reversed_edges.contains(&ShapeKey(old_edge.ptr_id())) {
                                let mut ve = old_edge.clone();
                                super::brep_fill_sweep::reverse_edge_in_first_or_last_wire(
                                    brep,
                                    &mut self.last_shape,
                                    &ve,
                                );
                                let _ = ve;
                            }
                        }
                    }
                }
            } // End of construction of edges
            ipath += 1;
            ipath_l += 1;
        }

        // (3) Construction of Faces
        let mut face = Shape::null();

        let mut ipath = 1i32;
        let mut ipath_l = ifirst;
        while ipath <= nb_path {
            for isec in 1..=nb_law {
                if degenerated.get(isec as usize, ipath as usize) {
                    let ue1 = u_edge.get(isec as usize, ipath as usize);
                    let ue2 = u_edge.get((isec + 1) as usize, ipath as usize);
                    if ue1.is_same(&ue2) {
                        self.my_faces.as_mut().expect("myFaces")[isec as usize - 1]
                            [(ipath_l - 1) as usize] = ue1;
                    } else {
                        let ve = v_edge.get(isec as usize, ipath as usize);
                        self.my_faces.as_mut().expect("myFaces")[isec as usize - 1]
                            [(ipath_l - 1) as usize] = ve;
                    }
                } else if self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                    [(ipath_l - 1) as usize]
                    .is_null()
                {
                    super::brep_fill_sweep_b::build_face(
                        brep,
                        tab_s.get(isec as usize, ipath as usize).surface(),
                        &u_edge.get(isec as usize, ipath as usize),
                        &v_edge.get(isec as usize, ipath as usize),
                        &u_edge.get((isec + 1) as usize, ipath as usize),
                        &v_edge.get(isec as usize, (ipath + 1) as usize),
                        &mut self.my_v_edges_modified,
                        exch_uv.get(isec as usize, ipath as usize),
                        u_reverse.get(isec as usize, ipath as usize),
                        &mut face,
                    );
                    self.my_faces.as_mut().expect("myFaces")[isec as usize - 1]
                        [(ipath_l - 1) as usize] = face.clone();
                }
            }
            ipath += 1;
            ipath_l += 1;
        }

        // (3.1) Reverse the faces that have been built earlier
        for ipath in 1..=nb_path {
            for isec in 1..=nb_law {
                if is_built[isec as usize - 1] {
                    // myFaces->ChangeValue(isec, ipath).Reverse();
                    let mut f = self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                        [(ipath - 1) as usize]
                        .clone();
                    f.orientation = match f.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        o => o,
                    };
                    self.my_faces.as_mut().expect("myFaces")[isec as usize - 1]
                        [(ipath - 1) as usize] = f;
                }
            }
        }

        // (4) History and Continuity
        if hasdegen {
            // (4.1) Degenerated case => Sledgehammer
            let mut comp = b.make_compound(brep, vec![]);
            for isec in 1..=(nb_law + 1) {
                let mut ipath = 1i32;
                let mut ipath_l = ifirst;
                while ipath <= nb_path + 1 {
                    if ipath <= nb_path {
                        self.my_u_edges.as_mut().expect("myUEdges")[isec as usize - 1]
                            [(ipath_l - 1) as usize] =
                            u_edge.get(isec as usize, ipath as usize);
                    }
                    if isec <= nb_law {
                        self.my_v_edges.as_mut().expect("myVEdges")[isec as usize - 1]
                            [(ipath_l - 1) as usize] =
                            v_edge.get(isec as usize, ipath as usize);
                    }
                    if (ipath <= nb_path)
                        && (isec <= nb_law)
                        && !self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                            [(ipath_l - 1) as usize]
                            .is_null()
                        && matches!(
                            self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                                [(ipath_l - 1) as usize]
                                .data
                                .as_ref(),
                            rcad_kernel::topo::topods::TShape::Face(_)
                        )
                    {
                        // B.Add(Comp, myFaces->Value(isec, IPath));
                        b.add_to_compound(
                            brep,
                            comp.clone(),
                            self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                                [(ipath_l - 1) as usize]
                                .clone(),
                        );
                    }
                    ipath += 1;
                    ipath_l += 1;
                }
            }
            // BRepLib::EncodeRegularity(Comp, myTolAngular);
            crate::topalgo::brep_lib_encode_regularity::encode_regularity(brep, &comp, self.my_tol_angular);
        } else {
            // (4.2) General case => Tweezers
            let mut is_g1;
            let mut ff;

            for isec in 1..=(nb_law + 1) {
                if isec > 1 {
                    is_g1 = self.my_sec.borrow().continuity(brep, isec - 1, self.my_tol_angular)
                        >= GeomAbsShape::G1;
                } else {
                    is_g1 = false;
                }
                let mut ipath = 1i32;
                let mut ipath_l = ifirst;
                while ipath <= nb_path {
                    self.my_u_edges.as_mut().expect("myUEdges")[isec as usize - 1]
                        [(ipath_l - 1) as usize] = u_edge.get(isec as usize, ipath as usize);
                    if is_g1 {
                        if isec == nb_law + 1 {
                            ff = self.my_faces.as_ref().expect("myFaces")[0][(ipath_l - 1) as usize]
                                .clone();
                        } else {
                            ff = self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                                [(ipath_l - 1) as usize]
                                .clone();
                        }
                        let an_edge = self.my_u_edges.as_ref().expect("myUEdges")
                            [isec as usize - 1][(ipath_l - 1) as usize]
                            .clone();
                        let f_prev = self.my_faces.as_ref().expect("myFaces")
                            [isec as usize - 2][(ipath_l - 1) as usize]
                            .clone();
                        b.continuity(brep, &an_edge, &f_prev, &ff, GeomAbsShape::G1);
                    }
                    ipath += 1;
                    ipath_l += 1;
                }
            }

            let mut nbpath = nb_path;
            if vclose {
                nbpath += 1; // Another test G1
            }
            let mut ipath = 1i32;
            let mut ipath_l = ifirst;
            while ipath <= nb_path + 1 {
                if (ipath > 1) && (ipath <= nbpath) {
                    is_g1 = self.my_loc.borrow().base().is_g1(
                        brep,
                        ipath_l - 1,
                        self.my_tol3d,
                        self.my_tol_angular,
                    ) >= 0;
                } else {
                    is_g1 = false;
                }
                for isec in 1..=nb_law {
                    self.my_v_edges.as_mut().expect("myVEdges")[isec as usize - 1]
                        [(ipath_l - 1) as usize] = v_edge.get(isec as usize, ipath as usize);
                    if is_g1 {
                        if ipath == nb_path + 1 {
                            ff = self.my_faces.as_ref().expect("myFaces")[isec as usize - 1][0]
                                .clone();
                        } else {
                            ff = self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                                [(ipath_l - 1) as usize]
                                .clone();
                        }
                        let an_edge = self.my_v_edges.as_ref().expect("myVEdges")
                            [isec as usize - 1][(ipath_l - 1) as usize]
                            .clone();
                        // BRepLib::EncodeRegularity(anEdge, FF,
                        //     TopoDS::Face(myFaces->Value(isec, IPath - 1)),
                        //     myTolAngular);
                        let f_prev = self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                            [(ipath_l - 2) as usize]
                            .clone();
                        crate::topalgo::brep_lib_encode_regularity::encode_regularity_edge(
                            brep, &an_edge, &ff, &f_prev, self.my_tol_angular,
                        );
                    }
                }
                ipath += 1;
                ipath_l += 1;
            }
        }

        // (5) Update Tapes and Rails
        if ifirst == 1 && !tapes.is_empty() {
            // works only in case of single shell
            for isec in 1..=nb_law {
                let start = start_edges[isec as usize - 1].clone();
                let tape = tapes.get_mut(&ShapeKey(start.ptr_id())).expect("Tapes(StartEdges)");
                for j in 1..=(nb_path + 1) {
                    tape[0][j as usize - 1] = self.my_v_edges.as_ref().expect("myVEdges")
                        [isec as usize - 1][j as usize - 1]
                        .clone();
                }
                for j in 1..=nb_path {
                    tape[1][j as usize - 1] = self.my_u_edges.as_ref().expect("myUEdges")
                        [isec as usize - 1][j as usize - 1]
                        .clone();
                }
                for j in 1..=nb_path {
                    tape[2][j as usize - 1] = self.my_u_edges.as_ref().expect("myUEdges")
                        [isec as usize][j as usize - 1]
                        .clone();
                }
                for j in 1..=(nb_path + 1) {
                    // (4) Vertex(isec, j)
                    tape[3][j as usize - 1] = vertex.get(isec as usize, j as usize);
                }
                for j in 1..=(nb_path + 1) {
                    // (5) Vertex(isec + 1, j)
                    tape[4][j as usize - 1] = vertex.get((isec + 1) as usize, j as usize);
                }
                for j in 1..=nb_path {
                    tape[5][j as usize - 1] = self.my_faces.as_ref().expect("myFaces")
                        [isec as usize - 1][j as usize - 1]
                        .clone();
                }

                // TopExp::Vertices(TopoDS::Edge(StartEdges(isec)), Vfirst,
                //                  Vlast, true); //with orientation
                let (vfirst, vlast) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(&start);
                if let Some(vs) = non_null(&vfirst) {
                    if !rails.contains_key(&ShapeKey(vs.ptr_id())) {
                        let mut an_array: ShapeHArray2 =
                            vec![vec![Shape::null(); (nb_path + 1) as usize]; 2];
                        for j in 1..=nb_path {
                            an_array[0][j as usize - 1] = self.my_u_edges
                                .as_ref()
                                .expect("myUEdges")[isec as usize - 1][j as usize - 1]
                                .clone();
                        }
                        for j in 1..=(nb_path + 1) {
                            an_array[1][j as usize - 1] = vertex.get(isec as usize, j as usize);
                        }
                        rails.insert(ShapeKey(vs.ptr_id()), an_array);
                    }
                }
                if let Some(vs) = non_null(&vlast) {
                    if !rails.contains_key(&ShapeKey(vs.ptr_id())) {
                        let mut an_array: ShapeHArray2 =
                            vec![vec![Shape::null(); (nb_path + 1) as usize]; 2];
                        for j in 1..=nb_path {
                            an_array[0][j as usize - 1] = self.my_u_edges
                                .as_ref()
                                .expect("myUEdges")[isec as usize][j as usize - 1]
                                .clone();
                        }
                        for j in 1..=(nb_path + 1) {
                            an_array[1][j as usize - 1] =
                                vertex.get((isec + 1) as usize, j as usize);
                        }
                        rails.insert(ShapeKey(vs.ptr_id()), an_array);
                    }
                }
            }
        }

        true
    }
}

/// The rcad `handle(Geom_Surface)` carrier for the TabS array (the OCCT
/// null handle maps to [`SurfaceShape::None`]).
#[derive(Debug, Clone)]
pub(super) enum SurfaceShape {
    None,
    Some(Surface3Value),
}

pub(super) type Surface3Value = rcad_kernel::geom::Surface3;

impl SurfaceShape {
    pub fn surface(&self) -> &rcad_kernel::geom::Surface3 {
        match self {
            SurfaceShape::Some(s) => s,
            SurfaceShape::None => panic!("null surface handle"),
        }
    }
}

/// The surface bounds (S->Bounds(UFirst, ULast, VFirst, VLast)).
pub(super) fn surface_bounds(surf: &rcad_kernel::geom::Surface3) -> (f64, f64, f64, f64) {
    let (_basis, bounds) =
        rcad_kernel::topo::topods::surface_adaptor_basis_and_bounds(surf);
    (bounds[0], bounds[1], bounds[2], bounds[3])
}

/// The OCCT null-handle test over an optional shape.
pub(super) fn non_null(s: &Shape) -> Option<Shape> {
    if s.is_null() {
        None
    } else {
        Some(s.clone())
    }
}
