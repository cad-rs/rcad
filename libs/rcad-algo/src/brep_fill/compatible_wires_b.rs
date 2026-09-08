//! OCCT BRepFill_CompatibleWires::Perform, SameNumberByPolarMethod,
//! SameNumberByACR, ComputeOrigin and SearchOrigin — split from
//! compatible_wires.rs (file-size rule): the inherent impl of the class
//! continues here; the struct and the static helpers live in
//! compatible_wires.rs.

use glam::DVec3;
use std::collections::HashMap;

use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{CurveEval, Plane};
use rcad_kernel::topo::topods::{BRep, Orientation, Shape, ShapeType, TShape};

use super::compatible_wires::{
    find_vertex_by_key, add_new_edge, brep_tool_parameter, build_connected_edges, compute_acr, continuity_rank,
    edge_intersect_on_wire, insert_acr, is_degenerated, plane_of_wire, sample_distance,
    sample_distance_backward, search_root, search_vertex, seq_of_vertices, set_wire_closed,
    shape_tolerance, top_exp_vertices_cumori, top_exp_vertices_stored, transform, vertex_point,
    wire_continuity, wire_edges, GeomAbsShapeAlias, WireFlags, GP_RESOLUTION, REAL_LAST,
};
use super::compatible_wires::CompatibleWires;
use crate::brep_fill::generator::{
    shape_key, shape_oriented, shape_reversed, top_exp_wire_vertices, BRepFillThruSectionErrorStatus,
    ShapeKey,
};

/// OCCT Precision::Confusion().
pub(super) const TOL_CONFUSION: f64 = CONFUSION;

impl CompatibleWires {
    /// OCCT BRepFill_CompatibleWires::Perform (L855-996).
    pub fn perform(&mut self, brep: &mut BRep, with_rotation: bool) {
        self.my_status = BRepFillThruSectionErrorStatus::Done;
        // compute origin and orientation on wires to avoid twisted results
        // and update wires to have same number of edges

        // determination of report:
        // if the number of elements is the same and if the wires have discontinuities
        // by tangency, the report is not carried out by curvilinear abscissa
        let nb_sects = self.my_work.len();
        let mut nbmax = 0usize;
        let mut nbmin = 0usize;
        let mut nb_edges: Vec<usize> = vec![0; nb_sects];
        let mut cont_s = GeomAbsShapeAlias::CN;
        for i in 0..nb_sects {
            let oriented = shape_oriented(&self.my_work[i], Orientation::Forward);
            self.my_work[i] = oriented;
            let w = self.my_work[i].clone();
            let cont = wire_continuity(brep, &w);
            if continuity_rank(cont) < continuity_rank(cont_s) {
                cont_s = cont;
            }
            nb_edges[i] = wire_edges(brep, &w).len();
            if i == 0 {
                nbmin = nb_edges[i];
            }
            nbmax = nbmax.max(nb_edges[i]);
            nbmin = nbmin.min(nb_edges[i]);
        }
        // if the number of elements is not the same or if all wires are at least
        // C1, the report is carried out by curvilinear abscissa of cuts, otherwise
        // a report vertex / Vertex is done
        let report = nbmax != nbmin || continuity_rank(cont_s) >= continuity_rank(GeomAbsShapeAlias::C1);

        // initialization of the map
        for i in 0..nb_sects {
            let w = self.my_work[i].clone();
            for e in wire_edges(brep, &w) {
                self.my_map.entry(shape_key(&e)).or_default().push(e.clone());
            }
        }

        // open/closed sections
        // initialisation of myDegen1, myDegen2
        let mut ideb = 0usize;
        let mut ifin = self.my_work.len().saturating_sub(1);
        // check if the first wire is punctual
        self.my_degen1 = true;
        for e in wire_edges(brep, &self.my_work[ideb].clone()) {
            self.my_degen1 = self.my_degen1 && is_degenerated(brep, &e);
        }
        if self.my_degen1 {
            ideb += 1;
        }
        // check if the last wire is punctual
        self.my_degen2 = true;
        for e in wire_edges(brep, &self.my_work[ifin].clone()) {
            self.my_degen2 = self.my_degen2 && is_degenerated(brep, &e);
        }
        if self.my_degen2 {
            ifin = ifin.saturating_sub(1);
        }

        let mut all_closed = true;
        let mut all_open = true;
        for i in ideb..=ifin {
            let mut wclosed = self.my_work[i].flags_closed(brep);
            if !wclosed {
                // check if the vertices are the same.
                let (v1, v2) = top_exp_wire_vertices(brep, &self.my_work[i]);
                if v1.is_same(&v2) {
                    wclosed = true;
                }
            }
            all_closed = all_closed && wclosed;
            all_open = all_open && !wclosed;
        }

        if all_closed {
            // All sections are closed
            if report {
                // same number of elements
                self.same_number_by_polar_method(brep, with_rotation);
            } else {
                // origin
                self.compute_origin(brep, false);
            }
        } else if all_open {
            // All sections are open
            // origin
            self.search_origin(brep);
            if self.my_status != BRepFillThruSectionErrorStatus::Done {
                return;
            }
            // same number of elements
            if report {
                self.same_number_by_acr(brep, report);
            }
        } else {
            // There are open and closed sections :
            // not processed
            self.my_status = BRepFillThruSectionErrorStatus::NotSameTopology;
        }
    }

    /// OCCT BRepFill_CompatibleWires::SameNumberByPolarMethod (L1008-1507).
    fn same_number_by_polar_method(&mut self, brep: &mut BRep, with_rotation: bool) {
        // initialisation
        let nb_sects = self.my_work.len();
        let mut edge_new_edges: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();

        let mut all_closed = true;
        let mut ideb = 0usize;
        let mut ifin = nb_sects.saturating_sub(1);

        for i in 0..nb_sects {
            // OCCT L1022-1023: BRepCheck_Wire::Closed() == BRepCheck_NoError.
            // GAP: BRepCheck_Wire is not ported (plan §0.6, reported); the
            // closure test uses the wire flag plus the equal-endpoints rule.
            let w = self.my_work[i].clone();
            let mut closed = w.flags_closed(brep);
            if !closed {
                let (v1, v2) = top_exp_wire_vertices(brep, &w);
                closed = v1.is_same(&v2);
            }
            all_closed = all_closed && closed;
        }
        if !all_closed {
            self.my_status = BRepFillThruSectionErrorStatus::NotSameTopology;
            return;
        }

        // sections ponctuelles, sections bouclantes ?
        if self.my_degen1 {
            ideb += 1;
        }
        if self.my_degen2 {
            ifin = ifin.saturating_sub(1);
        }
        let vclosed = !self.my_degen1 && !self.my_degen2 && self.my_work[ideb].is_same(&self.my_work[ifin]);

        // Removing degenerated edges
        for i in ideb..=ifin {
            let mut has_deg_edge = false;
            let edges = wire_edges(brep, &self.my_work[i]);
            for an_edge in &edges {
                if is_degenerated(brep, an_edge) {
                    has_deg_edge = true;
                    break;
                }
            }
            if has_deg_edge {
                let mut new_wire_edges: Vec<Shape> = Vec::new();
                for an_edge in &edges {
                    if !is_degenerated(brep, an_edge) {
                        new_wire_edges.push(an_edge.clone());
                    }
                }
                self.my_work[i] = brep.add_twire(new_wire_edges);
            }
        }

        // Nombre max de decoupes possibles
        let mut nb_max_v = 0usize;
        for i in 0..nb_sects {
            nb_max_v += wire_edges(brep, &self.my_work[i]).len();
        }

        // construction of tables of planes of wires
        let mut pos: Vec<DVec3> = vec![DVec3::ZERO; nb_sects];
        let mut axe: Vec<DVec3> = vec![DVec3::ZERO; nb_sects];
        let mut p = Plane::new(DVec3::ZERO, DVec3::Z);
        for i in ideb..=ifin {
            let w = self.my_work[i].clone();
            if plane_of_wire(brep, &w, &mut p) {
                pos[i] = p.origin;
                axe[i] = p.normal;
            }
        }
        if self.my_degen1 {
            let seq_v = seq_of_vertices(brep, &self.my_work[0]);
            pos[0] = vertex_point(brep, &seq_v[0]);
            axe[0] = axe[ideb];
        }
        if self.my_degen2 {
            let seq_v = seq_of_vertices(brep, &self.my_work[nb_sects - 1]);
            pos[nb_sects - 1] = vertex_point(brep, &seq_v[0]);
            axe[nb_sects - 1] = axe[ifin];
        }

        // construction of RMap, map of reports of wire i to wire i-1
        let mut rmap: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();

        // loop on i
        for i in (ideb + 1..=ifin).rev() {
            let wire1 = self.my_work[i].clone();

            // sequence of vertices of the first wire
            let seq_v = seq_of_vertices(brep, &wire1);
            if seq_v.len() > nb_max_v {
                self.my_status = BRepFillThruSectionErrorStatus::Failed;
                return;
            }

            // loop on vertices of wire1
            for vi in &seq_v {
                // init of RMap for Vi
                rmap.entry(shape_key(vi)).or_default();

                // it is required to find intersection Vi - wire2
                let pi = vertex_point(brep, vi);

                // return Pi in the current plane
                let mut pnew = DVec3::ZERO;
                transform(
                    with_rotation,
                    pi,
                    pos[i],
                    axe[i],
                    pos[i - 1],
                    axe[i - 1],
                    &mut pnew,
                );

                // calculate the intersection
                if pnew.distance(pos[i - 1]) > TOL_CONFUSION {
                    let percent = self.my_percent;
                    let mut vsol = Shape::null();
                    let mut newwire = Shape::null();
                    let wire_prev = self.my_work[i - 1].clone();
                    let new_vertex = edge_intersect_on_wire(
                        brep,
                        pos[i - 1],
                        pnew,
                        percent,
                        &rmap,
                        &wire_prev,
                        &mut vsol,
                        &mut newwire,
                        &mut edge_new_edges,
                    );
                    if new_vertex {
                        self.my_work[i - 1] = newwire;
                    }
                    rmap.entry(shape_key(vi)).or_default().push(vsol);
                }
            } // loop on  ii
        } // loop on  i

        // initialisation of MapVLV, map of correspondences vertex - list of vertices
        let mut map_vlv: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();
        let seq_v = seq_of_vertices(brep, &self.my_work[ideb]);
        let size_map = seq_v.len();
        for vi in &seq_v {
            map_vlv.insert(shape_key(vi), vec![vi.clone()]);
            let mut nb_v = 1usize;
            let mut v0 = vi.clone();
            let mut v1 = Shape::null();
            let mut tantque = search_root(&v0, &rmap, brep, &mut v1);
            while tantque {
                map_vlv.entry(shape_key(vi)).or_default().push(v1.clone());
                nb_v += 1;
                // test on NbV required for looping sections
                if v1.is_same(vi) || nb_v >= self.my_work.len() {
                    tantque = false;
                } else {
                    v0 = v1.clone();
                    tantque = search_root(&v0, &rmap, brep, &mut v1);
                }
            }
        }

        // loop on i
        for i in ideb..ifin {
            let wire1 = self.my_work[i].clone();

            // sequence of vertices of the first wire
            let seq_v = seq_of_vertices(brep, &wire1);
            if seq_v.len() > nb_max_v || seq_v.len() > size_map {
                self.my_status = BRepFillThruSectionErrorStatus::Failed;
                return;
            }

            // next wire
            let wire2 = self.my_work[i + 1].clone();

            // loop on vertices of wire1
            for vi in &seq_v {
                let mut vroot = Shape::null();
                let mut intersect = true;
                if search_root(vi, &map_vlv, brep, &mut vroot) {
                    let lvi = map_vlv.get(&shape_key(&vroot)).cloned().unwrap_or_default();
                    let mut von_w = Shape::null();
                    intersect = !search_vertex(brep, &lvi, &wire2, &mut von_w);
                }

                if intersect {
                    // it is necessary to find intersection Vi - wire2
                    let pi = vertex_point(brep, vi);

                    // return Pi in the current plane
                    let mut pnew = DVec3::ZERO;
                    transform(
                        with_rotation,
                        pi,
                        pos[i],
                        axe[i],
                        pos[i + 1],
                        axe[i + 1],
                        &mut pnew,
                    );

                    // calculate the intersection
                    if pnew.distance(pos[i + 1]) > TOL_CONFUSION {
                        let percent = self.my_percent;
                        let mut vsol = Shape::null();
                        let mut newwire = Shape::null();
                        let new_vertex = edge_intersect_on_wire(
                            brep,
                            pos[i + 1],
                            pnew,
                            percent,
                            &map_vlv,
                            &wire2,
                            &mut vsol,
                            &mut newwire,
                            &mut edge_new_edges,
                        );
                        map_vlv.entry(shape_key(&vroot)).or_default().push(vsol);
                        if new_vertex {
                            self.my_work[i + 1] = newwire;
                        }
                    }
                }
            } // loop on ii
        } // loop on i

        // regularize wires following MapVLV
        let wire = self.my_work[ideb].clone();
        let wire_edges_base = wire_edges(brep, &wire);

        // except for the last if the sections loop
        let mut ibout = ifin;
        if vclosed {
            ibout -= 1;
        }

        for i in (ideb + 1)..=ibout {
            let mut mw_edges: Vec<Shape> = Vec::new(); // BRepLib_MakeWire MW

            let mut an_exp_ix = 0usize; // anExp.Init(wire)
            let ecur = wire_edges_base[an_exp_ix].clone();
            let (mut vf, mut vl) = top_exp_vertices_cumori(brep, &ecur);
            let mut u1 = brep_tool_parameter(brep, &vf, &ecur);
            let mut u2 = brep_tool_parameter(brep, &vl, &ecur);
            let ecur_curve = brep.edge(ecur.clone()).curve.clone();
            let mut pps = ecur_curve
                .as_ref()
                .map(|c| c.point_at(0.1 * (u1 + 9.0 * u2)))
                .unwrap_or(DVec3::ZERO);
            let mut it_f_ix = 0usize;
            let mut it_l_ix = 0usize;
            let vf_list = map_vlv.get(&shape_key(&vf)).cloned().unwrap_or_default();
            let vl_list = map_vlv.get(&shape_key(&vl)).cloned().unwrap_or_default();
            let mut rang = ideb;
            while rang < i && it_f_ix < vf_list.len() && it_l_ix < vl_list.len() {
                it_f_ix += 1;
                it_l_ix += 1;
                rang += 1;
            }
            if it_f_ix >= vf_list.len() || it_l_ix >= vl_list.len() {
                // Correspondence chain shorter than the section index; fail gracefully.
                self.my_status = BRepFillThruSectionErrorStatus::Failed;
                return;
            }
            let v1 = vf_list[it_f_ix].clone();
            let v2 = vl_list[it_l_ix].clone();
            let mut esol: Option<Shape> = None;
            let mut scalmax = 0.0f64;
            for e in wire_edges(brep, &self.my_work[i]) {
                let (vvf, vvl) = top_exp_vertices_cumori(brep, &e);

                // parse candidate edges
                if (v1.is_same(&vvf) && v2.is_same(&vvl)) || (v2.is_same(&vvf) && v1.is_same(&vvl)) {
                    let u1param = brep_tool_parameter(brep, &vvf, &e);
                    let u2param = brep_tool_parameter(brep, &vvl, &e);
                    let curve_e = brep.edge(e.clone()).curve.clone();
                    let pp1 = curve_e
                        .as_ref()
                        .map(|c| c.point_at(0.1 * (u1param + 9.0 * u2param)))
                        .unwrap_or(DVec3::ZERO);
                    let pp2 = curve_e
                        .as_ref()
                        .map(|c| c.point_at(0.1 * (9.0 * u1param + u2param)))
                        .unwrap_or(DVec3::ZERO);

                    let mut pp1 = pp1;
                    let mut pp2 = pp2;
                    for rang in ((ideb + 1)..=i).rev() {
                        transform(with_rotation, pp1, pos[rang], axe[rang], pos[rang - 1], axe[rang - 1], &mut pp1);
                        transform(with_rotation, pp2, pos[rang], axe[rang], pos[rang - 1], axe[rang - 1], &mut pp2);
                    }
                    let mut ns = pps - pos[ideb];
                    ns = ns.normalize_or_zero();
                    let mut n1 = pp1 - pos[ideb];
                    n1 = n1.normalize_or_zero();
                    let mut n2 = pp2 - pos[ideb];
                    n2 = n2.normalize_or_zero();
                    let scal1 = n1.dot(ns);
                    if scal1 > scalmax {
                        scalmax = scal1;
                        esol = Some(e.clone());
                    }
                    let scal2 = n2.dot(ns);
                    if scal2 > scalmax {
                        scalmax = scal2;
                        esol = Some(shape_reversed(&e));
                    }
                }
            } // end of for(; itW.More(); itW.Next())
            match esol {
                Some(es) => {
                    mw_edges.push(es.clone());
                    let esol = es;
                    let mut connected_edges: Vec<Shape> = Vec::new();
                    build_connected_edges(brep, &self.my_work[i], &esol, &v2, &mut connected_edges);

                    let mut an_exp_ix = 0usize;
                    let mut it_ce_ix = 0usize;
                    while an_exp_ix < wire_edges_base.len() && it_ce_ix < connected_edges.len() {
                        let ecur = wire_edges_base[an_exp_ix].clone();
                        let (evf, evl) = top_exp_vertices_cumori(brep, &ecur);
                        let eu1 = brep_tool_parameter(brep, &evf, &ecur);
                        let eu2 = brep_tool_parameter(brep, &evl, &ecur);
                        let ecurve = brep.edge(ecur.clone()).curve.clone();
                        pps = ecurve
                            .as_ref()
                            .map(|c| c.point_at(0.1 * (eu1 + 9.0 * eu2)))
                            .unwrap_or(DVec3::ZERO);

                        let e = connected_edges[it_ce_ix].clone();
                        let (vvf, vvl) = top_exp_vertices_cumori(brep, &e);

                        // parse candidate edges
                        let u1 = brep_tool_parameter(brep, &vvf, &e);
                        let u2 = brep_tool_parameter(brep, &vvl, &e);
                        let ecurve = brep.edge(e.clone()).curve.clone();
                        let mut pp1 = ecurve
                            .as_ref()
                            .map(|c| c.point_at(0.1 * (u1 + 9.0 * u2)))
                            .unwrap_or(DVec3::ZERO);
                        let mut pp2 = ecurve
                            .as_ref()
                            .map(|c| c.point_at(0.1 * (9.0 * u1 + u2)))
                            .unwrap_or(DVec3::ZERO);

                        for rang in ((ideb + 1)..=i).rev() {
                            transform(with_rotation, pp1, pos[rang], axe[rang], pos[rang - 1], axe[rang - 1], &mut pp1);
                            transform(with_rotation, pp2, pos[rang], axe[rang], pos[rang - 1], axe[rang - 1], &mut pp2);
                        }
                        let mut ns = pps - pos[ideb];
                        ns = ns.normalize_or_zero();
                        let mut n1 = pp1 - pos[ideb];
                        n1 = n1.normalize_or_zero();
                        let mut n2 = pp2 - pos[ideb];
                        n2 = n2.normalize_or_zero();
                        let scal1 = n1.dot(ns);
                        let scal2 = n2.dot(ns);
                        let mut e = e;
                        if scal2 > scal1 {
                            e = shape_reversed(&e);
                        }
                        mw_edges.push(e);

                        an_exp_ix += 1;
                        it_ce_ix += 1;
                    }
                    self.my_work[i] = brep.add_twire(mw_edges);
                }
                None => {
                    self.my_status = BRepFillThruSectionErrorStatus::ProfilesInconsistent;
                    return;
                }
            }
            let _ = (&mut vf, &mut vl, &mut u1, &mut u2, &mut pps);
        }

        // blocking sections?
        if vclosed {
            let last_edges = wire_edges(brep, &self.my_work[self.my_work.len() - 1]);
            let first_edges = wire_edges(brep, &self.my_work[0]);
            for (an_edge, a_new_edge) in last_edges.iter().zip(first_edges.iter()) {
                if !an_edge.is_same(a_new_edge) {
                    edge_new_edges.insert(shape_key(an_edge), vec![a_new_edge.clone()]);
                }
            }
            let last = self.my_work.len() - 1;
            self.my_work[last] = self.my_work[0].clone();
        }

        // check the number of edges for debug
        let mut nbmax = 0usize;
        let mut nbmin = 0usize;
        for i in ideb..=ifin {
            let nb_edges = wire_edges(brep, &self.my_work[i]).len();
            if i == ideb {
                nbmin = nb_edges;
            }
            nbmax = nbmax.max(nb_edges);
            nbmin = nbmin.min(nb_edges);
        }
        if nbmin != nbmax {
            self.my_status = BRepFillThruSectionErrorStatus::Failed;
            return;
        }

        // Fill <myMap>
        let keys: Vec<Shape> = {
            let mut ks: Vec<Shape> = Vec::new();
            for (k, _) in self.my_map.iter() {
                ks.push(find_vertex_by_key(brep, *k));
            }
            ks
        };
        for an_edge in &keys {
            let mut list_of_new_edges: Vec<Shape> = Vec::new();

            // for each edge of <myMap> find all newest edges
            // in <EdgeNewEdges> recursively
            add_new_edge(an_edge, &edge_new_edges, &mut list_of_new_edges);

            self.my_map.insert(shape_key(an_edge), list_of_new_edges);
        }
        let _ = INFINITE_VALUE;
    }

    /// OCCT BRepFill_CompatibleWires::SameNumberByACR (L1511-1789).
    fn same_number_by_acr(&mut self, brep: &mut BRep, report: bool) {
        // find the dimension
        let mut ideb = 0usize;
        let mut ifin = self.my_work.len().saturating_sub(1);

        // point sections, blocking  sections?
        if self.my_degen1 {
            ideb += 1;
        }
        if self.my_degen2 {
            ifin = ifin.saturating_sub(1);
        }
        let vclosed = !self.my_degen1 && !self.my_degen2 && self.my_work[ideb].is_same(&self.my_work[ifin]);

        let nb_sects = self.my_work.len();
        let mut nbmax = 0usize;
        let mut nbmin = 0usize;
        let mut nb_edges: Vec<usize> = vec![0; nb_sects];
        for i in 0..nb_sects {
            nb_edges[i] = wire_edges(brep, &self.my_work[i]).len();
            if i == 0 {
                nbmin = nb_edges[i];
            }
            nbmax = nbmax.max(nb_edges[i]);
            nbmin = nbmin.min(nb_edges[i]);
        }

        if nbmax > 1 {
            // several edges

            if report || nbmin < nbmax {
                // insertion of cuts
                let nbdec_full = (nbmax - 1) * nb_sects + 1;
                let mut dec: Vec<f64> = vec![0.0; nbdec_full];
                dec[1] = 1.0;

                let mut wire_len: Vec<f64> = vec![0.0; nb_sects];

                // calculate the table of cuts
                for i in 0..nb_sects {
                    // current wire
                    let wire1 = self.my_work[i].clone();
                    let nb_e = wire_edges(brep, &wire1).len();
                    // length and ACR of the wire
                    let acr = compute_acr(brep, &wire1); // ACR(0..nbE)
                    wire_len[i] = acr[0];
                    // insertion of ACR of the wire in the table of cuts
                    for j in 1..nb_e {
                        // OCCT: for (j = 1; j < ACR.Length() - 1; j++) — the
                        // interior entries (1 .. nbE-1).
                        let acrj = acr[j];
                        let mut k = 1usize;
                        while dec[k] < acrj {
                            k += 1;
                            if k > nbdec_full - 1 {
                                break;
                            }
                        }
                        if dec[k - 1] < acrj && acrj < dec[k] {
                            // shift right from k-1 .. nbdec-2 and insert
                            for l in ((k - 1)..(nbdec_full - 1)).rev() {
                                dec[l + 1] = dec[l];
                            }
                            dec[k - 1] = acrj;
                        }
                    }
                }

                // table of cuts
                let mut k = 1usize;
                while dec[k] < 1.0 {
                    k += 1;
                    if k > nbdec_full - 1 {
                        break;
                    }
                }
                let nbdec = k - 1;
                let dec2: Vec<f64> = dec[1..=nbdec].to_vec();

                // Check of cuts: are all the new edges long enough or not
                let mut cuts_to_remove: Vec<usize> = Vec::new();
                for k in 0..nbdec {
                    let knot1 = dec2[k];
                    let knot2 = if k == nbdec - 1 { 1.0 } else { dec2[k + 1] };
                    let mut all_lengths_null = true;
                    for i in 0..nb_sects {
                        let edge_len = (knot2 - knot1) * wire_len[i];
                        if edge_len > TOL_CONFUSION {
                            all_lengths_null = false;
                            break;
                        }
                    }
                    if all_lengths_null {
                        cuts_to_remove.push(k);
                    }
                }
                let new_nb_dec = nbdec - cuts_to_remove.len();
                let mut dec3: Vec<f64> = Vec::with_capacity(new_nb_dec);
                for k in 0..nbdec {
                    if !cuts_to_remove.contains(&k) {
                        dec3.push(dec2[k]);
                    }
                }
                ///////////////////

                // insertion of cuts in each wire
                for i in 0..nb_sects {
                    let oldwire = self.my_work[i].clone();
                    let mut tol = TOL_CONFUSION;
                    if wire_len[i] > GP_RESOLUTION {
                        tol /= wire_len[i];
                    }
                    let newwire = insert_acr(brep, &oldwire, &dec3, tol);
                    let old_edges = wire_edges(brep, &oldwire);
                    let new_edges = wire_edges(brep, &newwire);
                    let mut an_exp1_ix = 0usize;
                    let mut an_exp2_ix = 0usize;
                    while an_exp1_ix < old_edges.len() {
                        let ecur = old_edges[an_exp1_ix].clone();
                        if !ecur.is_same(&new_edges[an_exp2_ix]) {
                            let mut le: Vec<Shape> = Vec::new();
                            let v1 = top_exp_vertices_cumori(brep, &old_edges[an_exp1_ix]).0;
                            let (mut vf, mut vr) = top_exp_vertices_cumori(brep, &ecur);
                            let mut p1 = DVec3::ZERO;
                            if v1.is_same(&vf) {
                                p1 = vertex_point(brep, &vr);
                            }
                            if v1.is_same(&vr) {
                                p1 = vertex_point(brep, &vf);
                            }
                            let mut v2 = top_exp_vertices_cumori(brep, &new_edges[an_exp2_ix]).0;
                            let (nvf, nvr) = top_exp_vertices_cumori(brep, &new_edges[an_exp2_ix]);
                            vf = nvf;
                            vr = nvr;
                            let mut p2 = DVec3::ZERO;
                            if v2.is_same(&vf) {
                                p2 = vertex_point(brep, &vr);
                            }
                            if v2.is_same(&vr) {
                                p2 = vertex_point(brep, &vf);
                            }
                            while p1.distance(p2) > 1.0e-3 {
                                le.push(new_edges[an_exp2_ix].clone());
                                an_exp2_ix += 1;
                                v2 = top_exp_vertices_cumori(brep, &new_edges[an_exp2_ix]).0;
                                let (nvf, nvr) = top_exp_vertices_cumori(brep, &new_edges[an_exp2_ix]);
                                if v2.is_same(&nvf) {
                                    p2 = vertex_point(brep, &nvr);
                                }
                                if v2.is_same(&nvr) {
                                    p2 = vertex_point(brep, &nvf);
                                }
                                if p1.distance(p2) <= 1.0e-3 {
                                    le.push(new_edges[an_exp2_ix].clone());
                                    an_exp2_ix += 1;
                                }
                            }

                            // find the ancestor whose list contains Ecur and
                            // splice LE in its place
                            let mut found = false;
                            let map_keys: Vec<ShapeKey> = self.my_map.keys().copied().collect();
                            for k in map_keys {
                                if found {
                                    break;
                                }
                                let mut itlist_ix = 0usize;
                                while itlist_ix < self.my_map[&k].len() && !found {
                                    if ecur.is_same(&self.my_map[&k][itlist_ix]) {
                                        let ancestor = k;
                                        // splice: InsertBefore(LE, itlist); Remove(itlist)
                                        let list = self.my_map.get_mut(&ancestor).unwrap();
                                        let mut new_list: Vec<Shape> = Vec::new();
                                        new_list.extend_from_slice(&list[..itlist_ix]);
                                        new_list.extend_from_slice(&le);
                                        new_list.extend_from_slice(&list[itlist_ix + 1..]);
                                        *list = new_list;
                                        found = true;
                                    }
                                    if !found {
                                        itlist_ix += 1;
                                    }
                                }
                            }
                        } else {
                            an_exp2_ix += 1;
                        }
                        an_exp1_ix += 1;
                    }
                    self.my_work[i] = newwire;
                }
            }
        }

        // blocking sections ?
        if vclosed {
            let last = self.my_work.len() - 1;
            self.my_work[last] = self.my_work[0].clone();
        }

        // check the number of edges for debug
        nbmax = 0;
        for i in ideb..=ifin {
            nb_edges[i] = wire_edges(brep, &self.my_work[i]).len();
            if i == ideb {
                nbmin = nb_edges[i];
            }
            nbmax = nbmax.max(nb_edges[i]);
            nbmin = nbmin.min(nb_edges[i]);
        }
        if nbmax != nbmin {
            self.my_status = BRepFillThruSectionErrorStatus::Failed;
        }
    }

    /// OCCT BRepFill_CompatibleWires::ComputeOrigin (L1793-2455).
    fn compute_origin(&mut self, brep: &mut BRep, _polar: bool) {
        // reorganize the wires respecting orientation and origin

        let mut all_closed = true;
        let nb_sects = self.my_work.len();
        let mut ideb = 0usize;
        let mut ifin = nb_sects.saturating_sub(1);

        // point sections, blocking sections
        if self.my_degen1 {
            ideb += 1;
        }
        if self.my_degen2 {
            ifin = ifin.saturating_sub(1);
        }
        let vclosed = !self.my_degen1 && !self.my_degen2 && self.my_work[ideb].is_same(&self.my_work[ifin]);

        for i in ideb..=ifin {
            let mut wclosed = self.my_work[i].flags_closed(brep);
            if !wclosed {
                // check if the vertices are the same.
                let (v1, v2) = top_exp_wire_vertices(brep, &self.my_work[i]);
                if v1.is_same(&v2) {
                    wclosed = true;
                }
            }
            all_closed = all_closed && wclosed;
        }
        if !all_closed {
            self.my_status = BRepFillThruSectionErrorStatus::NotSameTopology;
            return;
        }

        // Consider that all wires have same number of edges (polar==false)
        let mut prev_seq: Vec<Shape> = Vec::new();
        let mut prev_eseq: Vec<Shape> = Vec::new();
        let wire = self.my_work[ideb].clone();
        for e in wire_edges(brep, &wire) {
            // PrevSeq.Append(anExp.CurrentVertex())
            prev_seq.push(top_exp_vertices_cumori(brep, &e).0);
            prev_eseq.push(e);
        }
        let the_length = prev_seq.len();

        let mut nb_samples = 0usize;
        if the_length <= 2 {
            nb_samples = 4;
        }
        let mut first_plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let first_wire = self.my_work[ideb].clone();
        let _is_plane = plane_of_wire(brep, &first_wire, &mut first_plane);
        let mut prev_bary = first_plane.origin;
        let normal_of_first_plane = first_plane.normal;
        for i in (ideb + 1)..=ifin {
            let a_wire = self.my_work[i].clone();

            // Compute offset vector as current bary center projected on first plane
            // to first bary center
            let mut cur_plane = Plane::new(DVec3::ZERO, DVec3::Z);
            let _is_plane = plane_of_wire(brep, &a_wire, &mut cur_plane);
            let mut cur_bary = cur_plane.origin;
            let a_vec = cur_bary - prev_bary;
            let an_offset_proj = normal_of_first_plane * (a_vec.dot(normal_of_first_plane));
            cur_bary -= an_offset_proj; // projected current bary center
            let offset = cur_bary - prev_bary;

            // BB.MakeWire(newwire)
            let mut newwire_edges: Vec<Shape> = Vec::new();

            let wire_edges_list = wire_edges(brep, &a_wire);
            let mut seq_vertices: Vec<Shape> = Vec::new();
            let mut seq_edges: Vec<Shape> = Vec::new();
            for e in &wire_edges_list {
                seq_vertices.push(top_exp_vertices_cumori(brep, e).0);
                seq_edges.push(e.clone());
            }

            let mut min_sum_dist = INFINITE_VALUE;
            let mut jmin = 1usize;
            let mut forward = false;
            if i == self.my_work.len() - 1 && self.my_degen2 {
                // last point section
                jmin = 1;
                forward = true;
            } else {
                for j in 1..=the_length {
                    // Forward
                    let mut sum_dist = 0.0f64;
                    let mut n = 1usize;
                    for k in j..=the_length {
                        let vprev = &prev_seq[n - 1];
                        let pprev = vertex_point(brep, vprev);
                        let v = &seq_vertices[k - 1];
                        let p = vertex_point(brep, v) + offset;
                        sum_dist += pprev.distance(p);
                        if nb_samples > 0 {
                            let prev_edge = &prev_eseq[n - 1];
                            let cur_edge = &seq_edges[k - 1];
                            sum_dist += sample_distance(brep, prev_edge, cur_edge, nb_samples, offset);
                        }
                        n += 1;
                    }
                    for k in 1..j {
                        let vprev = &prev_seq[n - 1];
                        let pprev = vertex_point(brep, vprev);
                        let v = &seq_vertices[k - 1];
                        let p = vertex_point(brep, v) + offset;
                        sum_dist += pprev.distance(p);
                        if nb_samples > 0 {
                            let prev_edge = &prev_eseq[n - 1];
                            let cur_edge = &seq_edges[k - 1];
                            sum_dist += sample_distance(brep, prev_edge, cur_edge, nb_samples, offset);
                        }
                        n += 1;
                    }
                    if sum_dist < min_sum_dist {
                        min_sum_dist = sum_dist;
                        jmin = j;
                        forward = true;
                    }

                    // Backward
                    sum_dist = 0.0;
                    let mut n = 1usize;
                    let mut k = j as isize;
                    while k >= 1 {
                        let vprev = &prev_seq[n - 1];
                        let pprev = vertex_point(brep, vprev);
                        let v = &seq_vertices[(k - 1) as usize];
                        let p = vertex_point(brep, v) + offset;
                        sum_dist += pprev.distance(p);
                        if nb_samples > 0 {
                            // int k_cur = k - 1; if (k_cur == 0) k_cur = theLength;
                            let k_cur = if k - 1 == 0 { the_length } else { (k - 1) as usize };
                            let prev_edge = &prev_eseq[n - 1];
                            let cur_edge = &seq_edges[k_cur - 1];
                            sum_dist += sample_distance_backward(brep, prev_edge, cur_edge, nb_samples, offset);
                        }
                        n += 1;
                        k -= 1;
                    }
                    for k in ((j + 1)..=the_length).rev() {
                        let vprev = &prev_seq[n - 1];
                        let pprev = vertex_point(brep, vprev);
                        let v = &seq_vertices[k - 1];
                        let p = vertex_point(brep, v) + offset;
                        sum_dist += pprev.distance(p);
                        if nb_samples > 0 {
                            let prev_edge = &prev_eseq[n - 1];
                            let cur_edge = &seq_edges[k - 1 - 1];
                            sum_dist += sample_distance_backward(brep, prev_edge, cur_edge, nb_samples, offset);
                        }
                        n += 1;
                    }
                    if sum_dist < min_sum_dist {
                        min_sum_dist = sum_dist;
                        jmin = j;
                        forward = false;
                    }
                }
            }

            prev_seq.clear();
            prev_eseq.clear();
            if forward {
                for j in jmin..=the_length {
                    newwire_edges.push(seq_edges[j - 1].clone());
                    prev_seq.push(seq_vertices[j - 1].clone());
                    prev_eseq.push(seq_edges[j - 1].clone());
                }
                for j in 1..jmin {
                    newwire_edges.push(seq_edges[j - 1].clone());
                    prev_seq.push(seq_vertices[j - 1].clone());
                    prev_eseq.push(seq_edges[j - 1].clone());
                }
            } else {
                for j in (1..=(jmin - 1)).rev() {
                    newwire_edges.push(shape_reversed(&seq_edges[j - 1]));
                    // PrevSeq.Append( SeqVertices(j) );
                    prev_eseq.push(shape_reversed(&seq_edges[j - 1]));
                }
                for j in ((jmin)..=the_length).rev() {
                    newwire_edges.push(shape_reversed(&seq_edges[j - 1]));
                    // PrevSeq.Append( SeqVertices(j) );
                    prev_eseq.push(shape_reversed(&seq_edges[j - 1]));
                }
                for j in (1..=jmin).rev() {
                    prev_seq.push(seq_vertices[j - 1].clone());
                }
                for j in ((jmin + 1)..=the_length).rev() {
                    prev_seq.push(seq_vertices[j - 1].clone());
                }
            }

            let mut newwire = brep.add_twire(newwire_edges);
            // newwire.Closed(true); newwire.Orientation(TopAbs_FORWARD);
            set_wire_closed(brep, &newwire, true);
            newwire = shape_oriented(&newwire, Orientation::Forward);
            self.my_work[i] = newwire;

            prev_bary = cur_bary;
        }

        // blocking sections ?
        if vclosed {
            let last = self.my_work.len() - 1;
            self.my_work[last] = self.my_work[0].clone();
        }
    }

    /// OCCT BRepFill_CompatibleWires::SearchOrigin (L2459-2635).
    fn search_origin(&mut self, brep: &mut BRep) {
        // reorganize the open wires respecting orientation and origin

        let mut all_open = true;
        let mut ideb = 0usize;
        let mut ifin = self.my_work.len().saturating_sub(1);
        if self.my_degen1 {
            ideb += 1;
        }
        if self.my_degen2 {
            ifin = ifin.saturating_sub(1);
        }
        let vclosed = !self.my_degen1 && !self.my_degen2 && self.my_work[ideb].is_same(&self.my_work[ifin]);

        for i in ideb..=ifin {
            all_open = all_open && !self.my_work[i].flags_closed(brep);
        }
        if !all_open {
            self.my_status = BRepFillThruSectionErrorStatus::NotSameTopology;
            return;
        }

        // init
        let wire1 = self.my_work[ideb].clone();
        let wire1 = shape_oriented(&wire1, Orientation::Forward);
        let (vdeb, vfin) = top_exp_wire_vertices(brep, &wire1);
        let mut pdeb = vertex_point(brep, &vdeb);
        let mut pfin = vertex_point(brep, &vfin);
        let mut p0 = Plane::new(DVec3::ZERO, DVec3::Z);
        let mut p = Plane::new(DVec3::ZERO, DVec3::Z);
        let isline0 = !plane_of_wire(brep, &wire1, &mut p0);
        self.my_work[ideb] = wire1.clone();
        // OCC86
        let mut e0 = wire_edges(brep, &wire1)[0].clone();

        for i in (ideb + 1)..=ifin {
            let wire = shape_oriented(&self.my_work[i].clone(), Orientation::Forward);

            let mut seq_edges: Vec<Shape> = Vec::new();
            let mut nb_edges = 0usize;
            let first_e = wire_edges(brep, &wire)[0].clone();
            for e in wire_edges(brep, &wire) {
                seq_edges.push(e);
                nb_edges += 1;
            }
            let (vdeb, vfin) = top_exp_wire_vertices(brep, &wire);
            let isline = !plane_of_wire(brep, &wire, &mut p);

            let parcours;

            if isline0 || isline {
                // particular case of straight segments
                let p1 = vertex_point(brep, &vdeb);
                let p2 = vertex_point(brep, &vfin);
                let dist1 = pdeb.distance(p1) + pfin.distance(p2);
                let dist2 = pdeb.distance(p2) + pfin.distance(p1);
                parcours = dist2 >= dist1;
            } else {
                // OCC86
                let p1 = vertex_point(brep, &vdeb);
                let mut p1o = pdeb;
                let mut p2 = vertex_point(brep, &vfin);
                let mut p2o = pfin;
                if p1.distance(p2) < TOL_CONFUSION || p1o.distance(p2o) < TOL_CONFUSION {
                    // BRepAdaptor_Curve Curve0(E0), Curve(E);
                    // Curve0.D0(Curve0.FirstParameter() + Precision::Confusion(), P2o);
                    // Curve.D0(Curve.FirstParameter() + Precision::Confusion(), P2);
                    let c0 = brep.edge(e0.clone()).curve.clone();
                    let ce = brep.edge(first_e.clone()).curve.clone();
                    if let Some(c0) = c0 {
                        let r0 = brep.edge(e0.clone()).range;
                        p2o = c0.point_at(r0[0] + TOL_CONFUSION);
                    }
                    if let Some(ce) = ce {
                        let re = brep.edge(first_e.clone()).range;
                        p2 = ce.point_at(re[0] + TOL_CONFUSION);
                    }
                }
                let v_deb_fin0 = p2o - p1o;
                let v_deb_fin = p2 - p1;
                let la = v_deb_fin0.length();
                let lb = v_deb_fin.length();
                let a_straight = if la < 1e-300 || lb < 1e-300 {
                    0.0
                } else {
                    (v_deb_fin0.dot(v_deb_fin) / (la * lb)).clamp(-1.0, 1.0).acos()
                };
                parcours = a_straight < std::f64::consts::PI / 2.0;
            }

            // reconstruction of the wire
            let mut newwire_edges: Vec<Shape> = Vec::new();
            if parcours {
                for rang in 0..nb_edges {
                    newwire_edges.push(seq_edges[rang].clone());
                }
            } else {
                for rang in (0..nb_edges).rev() {
                    newwire_edges.push(shape_reversed(&seq_edges[rang]));
                }
            }

            // orientation of the wire
            let mut newwire = brep.add_twire(newwire_edges);
            // OCCT L2610: newwire.Oriented(TopAbs_FORWARD) — the ORIENTED
            // copy is stored (the call result is used, unlike a mutation).
            newwire = shape_oriented(&newwire, Orientation::Forward);
            self.my_work[i] = newwire;

            // passe to the next wire
            if parcours {
                pdeb = vertex_point(brep, &vdeb);
                pfin = vertex_point(brep, &vfin);
            } else {
                pfin = vertex_point(brep, &vdeb);
                pdeb = vertex_point(brep, &vfin);
            }
            p0 = p;
            // isline0 = isline;
            // OCC86
            e0 = first_e;
        }

        // blocking sections ?
        if vclosed {
            let last = self.my_work.len() - 1;
            self.my_work[last] = self.my_work[0].clone();
        }
        let _ = PCONFUSION;
    }
}
