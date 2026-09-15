//! OCCT BRepBuilderAPI_Sewing.cxx — the candidate evaluation group:
//! `EvaluateAngulars` (L1190-1286), `EvaluateDistances` (L1288-1524),
//! `IsMergedClosed` (L1526-1617), `AnalysisNearestEdges` (L1619-1784),
//! `FindCandidates` (L1786-2088) and `MergedNearestEdges` (L4243-4439).
//!
//! Re-hosts:
//! - `BndLib_Add2dCurve::Add(aC2d, f, l, tol, B)` (IsMergedClosed) — the 2d
//!   bounding box of the pcurve subrange; the rcad re-host samples the
//!   pcurve over the range into a `BndBox2d` (the
//!   bop/algo/pave_filler.rs `shrunk_range_bnd_box` precedent; the sampled
//!   box is a documented approximation of the OCCT analytic envelope).
//! - `GCPnts_UniformAbscissa` / `GCPnts_AbscissaPoint::Length` — the
//!   sewing.rs re-hosts.

use glam::DVec3;
use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use std::collections::HashSet;

use crate::brep_algo::tool as bat;

use super::BRepBuilderAPISewing;
use super::{gcpnts_abscissa_length, gcpnts_uniform_abscissa_parameters, set_add, set_add_bool, set_contains};

impl BRepBuilderAPISewing {
    // OCCT note: EvaluateAngulars' call sites are commented out in this
    // OCCT revision; the method surface is kept (dead-code allowed).
    /// OCCT BRepBuilderAPI_Sewing::EvaluateAngulars(sequenceSec, secForward,
    /// tabAng, indRef) (cxx L1190-1286) — called from MergingOfSections only.
    #[allow(dead_code)]
    pub(crate) fn evaluate_angulars(
        &self,
        brep: &BRep,
        sequence_sec: &[Shape],
        sec_forward: &[bool],
        tab_ang: &mut [f64],
        ind_ref: usize,
    ) {
        // OCCT L1191: tabAng.Init(-1.0).
        for v in tab_ang.iter_mut() {
            *v = -1.0;
        }

        // OCCT L1194: int i, j, npt = 4, lengSec = sequenceSec.Length().
        let npt = 4usize;
        let leng_sec = sequence_sec.len();

        // OCCT L1205: NCollection_Array1<gp_Vec> normRef(1, npt).
        let mut norm_ref = vec![DVec3::ZERO; npt];

        // OCCT L1207: for (i = indRef; i <= lengSec; i++).
        for i in ind_ref..=leng_sec {
            // OCCT L1209: edge = TopoDS::Edge(sequenceSec(i)).
            let edge = sequence_sec[i - 1].clone();

            // OCCT L1211-1226.
            let (face, surf, c2d) = {
                let mut bnd = edge.clone();
                if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd)) {
                    bnd = v.clone();
                }
                match self.my_bound_faces.get(&bat::shape_key(&bnd)) {
                    Some((_, faces)) => {
                        let face = faces.first().cloned().unwrap_or_else(Shape::null);
                        // OCCT L1225: surf = BRep_Tool::Surface(face, loc) +
                        // the location transform.
                        let surf = super::brep_tool_surface_world(brep, &face);
                        let c2d = bat::brep_tool_curve_on_surface(&edge, &face)
                            .map(|(c, _f, _l)| c);
                        (face, surf, c2d)
                    }
                    None => {
                        if i == ind_ref {
                            // OCCT L1229: else if (i == indRef) return;
                            return;
                        }
                        // OCCT L1231: else continue;
                        continue;
                    }
                }
            };
            let _ = face;
            let surf = match surf {
                Some(s) => s,
                None => continue,
            };
            let c2d = match c2d {
                Some(c) => c,
                None => continue,
            };

            // OCCT L1234-1240: c3d = BRep_Tool::Curve(edge, loc, first, last)
            // + the location transform.
            let (c3d, first, last) = match super::brep_tool_curve_world(brep, &edge) {
                Some(v) => v,
                None => continue,
            };

            // OCCT L1245-1246: GCPnts_UniformAbscissa uniAbs(adapt, npt,
            // first, last).
            let uni_abs = gcpnts_uniform_abscissa_parameters(&c3d, npt, first, last);

            // OCCT L1248-1250.
            let mut cumulate_angular = 0.0f64;
            let mut nb_computed_angle = 0i32;

            // OCCT L1252-1274.
            for j in 1..=npt {
                // OCCT L1254-1256.
                let par = uni_abs[(if sec_forward[i - 1] || i == ind_ref { j } else { npt - j + 1 }) - 1];
                let p = c2d.point_at(par);
                // OCCT L1259: surf->D1(P.X(), P.Y(), unused, w1, w2).
                let (_unused, w1, w2) = surf.derivatives(p.x, p.y);
                // OCCT L1260: gp_Vec n = w1 ^ w2.
                let n: DVec3 = w1.cross(w2);
                if i == ind_ref {
                    norm_ref[j - 1] = n;
                } else if n.length() > rcad_kernel::core::precision::CONFUSION
                    && norm_ref[j - 1].length() > rcad_kernel::core::precision::CONFUSION
                {
                    // OCCT L1265-1271: angular = n.Angle(normRef(j)).
                    nb_computed_angle += 1;
                    let mut angular =
                        (n.dot(norm_ref[j - 1]) / (n.length() * norm_ref[j - 1].length()))
                            .clamp(-1.0, 1.0)
                            .acos();
                    if angular > std::f64::consts::PI / 2.0 {
                        angular = std::f64::consts::PI - angular;
                    }
                    cumulate_angular += angular;
                }
            }

            // OCCT L1277-1281.
            if nb_computed_angle != 0 {
                tab_ang[i - 1] = cumulate_angular / (nb_computed_angle as f64);
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::EvaluateDistances(sequenceSec, secForward,
    /// tabDst, arrLen, tabMinDist, indRef) (cxx L1288-1524) — evaluate the
    /// distance between the edge with index indRef and the following edges.
    pub(crate) fn evaluate_distances(
        &self,
        brep: &BRep,
        sequence_sec: &[Shape],
        sec_forward: &mut [bool],
        tab_dst: &mut [f64],
        arr_len: &mut [f64],
        tab_min_dist: &mut [f64],
        ind_ref: usize,
    ) {
        // OCCT L1297-1300: the Init() calls.
        for v in sec_forward.iter_mut() {
            *v = true;
        }
        for v in tab_dst.iter_mut() {
            *v = -1.0;
        }
        for v in arr_len.iter_mut() {
            *v = 0.0;
        }
        for v in tab_min_dist.iter_mut() {
            *v = INFINITE_VALUE;
        }
        // OCCT L1301: const int npt = 8.
        let npt = 8usize;
        // OCCT L1302-1303: ptsRef / ptsSec.
        let mut pts_ref = vec![DVec3::ZERO; npt];
        let mut pts_sec = vec![DVec3::ZERO; npt];

        let leng_sec = sequence_sec.len();

        // OCCT L1310: double firstRef = 0., lastRef = 0.
        let mut c3d_ref: Option<rcad_kernel::geom::Curve3> = None;
        let mut first_ref = 0.0f64;
        let mut last_ref = 0.0f64;

        // OCCT L1312: for (i = indRef; i <= lengSec; i++).
        for i in ind_ref..=leng_sec {
            // reading of the edge (attention for the first one: reference)
            // OCCT L1315-1326.
            let sec = sequence_sec[i - 1].clone();
            let (c3d, first, last) = match super::brep_tool_curve_world(brep, &sec) {
                Some(v) => v,
                None => continue,
            };

            if i == ind_ref {
                c3d_ref = Some(c3d.clone());
                first_ref = first;
                last_ref = last;
            }

            // OCCT L1330-1334.
            let mut dist = INFINITE_VALUE;
            let mut dist_for = -1.0f64;
            let mut dist_rev = -1.0f64;
            let mut a_min_dist = INFINITE_VALUE;

            // OCCT L1337-1338.
            let delta_t = (last - first) / (npt as f64 - 1.0);
            let mut a_len_sec2 = 0.0f64;

            let mut nb_found = 0i32;
            for j in 1..=npt {
                // Uniform parameter on curve
                // OCCT L1342-1350.
                let t = if j == 1 {
                    first
                } else if j == npt {
                    last
                } else {
                    first + (j as f64 - 1.0) * delta_t
                };

                // Take point on curve
                // OCCT L1353.
                let pt = c3d.point_at(t);

                if i == ind_ref {
                    // OCCT L1356-1361.
                    pts_ref[j - 1] = pt;
                    if j > 1 {
                        a_len_sec2 += pt.distance_squared(pts_ref[j - 2]);
                    }
                } else {
                    // OCCT L1363-1387.
                    pts_sec[j - 1] = pt;
                    // protection to avoid merging with small sections
                    if j > 1 {
                        a_len_sec2 += pt.distance_squared(pts_sec[j - 2]);
                    }
                    // To evaluate mutual orientation and distance
                    dist = pt.distance(pts_ref[j - 1]);
                    if a_min_dist > dist {
                        a_min_dist = dist;
                    }
                    if dist_for < dist {
                        dist_for = dist;
                    }
                    dist = pt.distance(pts_ref[npt - j]);
                    if a_min_dist > dist {
                        a_min_dist = dist;
                    }
                    if dist_rev < dist {
                        dist_rev = dist;
                    }

                    // Check that point lays between vertices of reference curve
                    // OCCT L1390-1397.
                    let p11 = pts_ref[0];
                    let p12 = pts_ref[npt - 1];
                    let a_vec1 = p11 - pt; // gp_Vec aVec1(pt, p11)
                    let a_vec2 = p12 - pt; // gp_Vec aVec2(pt, p12)
                    let a_vec_ref = p12 - p11; // gp_Vec aVecRef(p11, p12)
                    if (a_vec_ref.dot(a_vec1)) * (a_vec_ref.dot(a_vec2)) < 0.0 {
                        nb_found += 1;
                    }
                }
            }

            // OCCT L1400-1404.
            let a_len_sec = a_len_sec2.sqrt();
            arr_len[i - 1] = a_len_sec;
            // Record mutual orientation
            // OCCT L1405-1406.
            let is_forward = dist_for < dist_rev; // szv debug: <=
            sec_forward[i - 1] = is_forward;

            // OCCT L1408-1416.
            dist = if is_forward { dist_for } else { dist_rev };
            if i == ind_ref || (dist < self.my_tolerance && (nb_found as f64) >= npt as f64 * 0.5) {
                tab_dst[i - 1] = dist;
                tab_min_dist[i - 1] = a_min_dist;
            } else {
                // OCCT L1417-1444.
                nb_found = 0;
                a_min_dist = INFINITE_VALUE;
                dist = -1.0;
                let mut arr_proj = vec![DVec3::ZERO; npt];
                let mut arr_dist = vec![-1.0f64; npt];
                let mut arr_para = vec![0.0f64; npt];
                let c3d_ref_v = c3d_ref.clone().unwrap_or_else(|| c3d.clone());
                if arr_len[ind_ref - 1] >= arr_len[i - 1] {
                    self.project_points_on_curve(
                        &pts_sec, &c3d_ref_v, first_ref, last_ref, &mut arr_dist, &mut arr_para,
                        &mut arr_proj, false,
                    );
                } else {
                    self.project_points_on_curve(
                        &pts_ref, &c3d, first, last, &mut arr_dist, &mut arr_para, &mut arr_proj,
                        false,
                    );
                }
                for j in 1..=npt {
                    if arr_dist[j - 1] < 0.0 {
                        continue;
                    }
                    if dist < arr_dist[j - 1] {
                        dist = arr_dist[j - 1];
                    }
                    if a_min_dist > arr_dist[j - 1] {
                        a_min_dist = arr_dist[j - 1];
                    }
                    nb_found += 1;
                }
                if nb_found > 1 {
                    tab_dst[i - 1] = dist;
                    tab_min_dist[i - 1] = a_min_dist;
                }
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::IsMergedClosed(Edge1, Edge2, face)
    /// (cxx L1526-1617).
    pub(crate) fn is_merged_closed(
        &self,
        brep: &BRep,
        edge1: &Shape,
        edge2: &Shape,
        face: &Shape,
    ) -> bool {
        // Check for closed surface
        // OCCT L1529-1533.
        let (surf_v, loc) = super::brep_tool_surface_loc(face);
        let surf = match surf_v {
            Some(s) => s,
            None => return false,
        };
        let is_u_closed = self.is_u_closed_surface(brep, &surf, edge1, face, loc);
        let is_v_closed = self.is_v_closed_surface(brep, &surf, edge1, face, loc);
        if !is_u_closed && !is_v_closed {
            return false;
        }
        // OCCT L1545-1550.
        let c2d1 = bat::brep_tool_curve_on_surface(edge1, face);
        let c2d2 = bat::brep_tool_curve_on_surface(edge2, face);
        let (c2d1, first2d1, last2d1) = match c2d1 {
            Some(v) => v,
            None => return false,
        };
        let (c2d2, first2d2, last2d2) = match c2d2 {
            Some(v) => v,
            None => return false,
        };

        // OCCT L1565-1586 (the bracketed locals): the 2d bounding boxes.
        let (c1_umin, c1_vmin, c1_umax, c1_vmax) = {
            let mut b1 = BndBox2d::new();
            bnd_lib_add2d_curve(&c2d1, first2d1, last2d1, PCONFUSION, &mut b1);
            b1.get().unwrap_or((0.0, 0.0, 0.0, 0.0))
        };
        let (c2_umin, c2_vmin, c2_umax, c2_vmax) = {
            let mut b2 = BndBox2d::new();
            bnd_lib_add2d_curve(&c2d2, first2d2, last2d2, PCONFUSION, &mut b2);
            b2.get().unwrap_or((0.0, 0.0, 0.0, 0.0))
        };
        let is_u_long_c1 = (c1_vmax - c1_vmin) <= (c1_umax - c1_umin);
        let is_v_long_c1 = (c1_umax - c1_umin) <= (c1_vmax - c1_vmin);
        let is_u_long_c2 = (c2_vmax - c2_vmin) <= (c2_umax - c2_umin);
        let is_v_long_c2 = (c2_umax - c2_umin) <= (c2_vmax - c2_vmin);
        // OCCT L1586: surf->Bounds(SUmin, SUmax, SVmin, SVmax).
        let [s_umin, s_umax, s_vmin, s_vmax] = surf.default_domain();

        // OCCT L1588-1600.
        if is_u_closed && is_v_long_c1 && is_v_long_c2 {
            // Do not merge if not overlapped by V
            let dist = (c2_vmin - c1_vmax).max(c1_vmin - c2_vmax);
            if dist < 0.0 {
                let dist_inner = (c2_umin - c1_umax).max(c1_umin - c2_umax);
                let dist_outer = (s_umax - s_umin) - (c2_umax - c1_umin).max(c1_umax - c2_umin);
                if dist_outer <= dist_inner {
                    return true;
                }
            }
        }
        // OCCT L1601-1613.
        if is_v_closed && is_u_long_c1 && is_u_long_c2 {
            // Do not merge if not overlapped by U
            let dist = (c2_umin - c1_umax).max(c1_umin - c2_umax);
            if dist < 0.0 {
                let dist_inner = (c2_vmin - c1_vmax).max(c1_vmin - c2_vmax);
                let dist_outer = (s_vmax - s_vmin) - (c2_vmax - c1_vmin).max(c1_vmax - c2_vmin);
                if dist_outer <= dist_inner {
                    return true;
                }
            }
        }
        // OCCT L1616.
        false
    }

    /// OCCT BRepBuilderAPI_Sewing::AnalysisNearestEdges(sequenceSec,
    /// seqIndCandidate, seqOrientations, evalDist) (cxx L1619-1784).
    pub(crate) fn analysis_nearest_edges(
        &mut self,
        brep: &BRep,
        sequence_sec: &[Shape],
        seq_ind_candidate: &mut Vec<i32>,
        seq_orientations: &mut Vec<bool>,
        eval_dist: bool,
    ) {
        // OCCT L1621: int workIndex = seqIndCandidate.First().
        let work_index = seq_ind_candidate[0];
        // OCCT L1622-1631.
        let workedge = sequence_sec[(work_index - 1) as usize].clone();
        let mut bnd = workedge.clone();
        let mut workfaces: Vec<Shape> = Vec::new();
        if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd)) {
            bnd = v.clone();
        }
        if let Some((_, faces)) = self.my_bound_faces.get(&bat::shape_key(&bnd)) {
            workfaces = faces.clone();
        }
        if workfaces.is_empty() {
            return;
        }
        // OCCT L1633-1639.
        let mut map_faces: super::ShapeSet = super::ShapeSet::new();
        for l_it in &workfaces {
            set_add_bool(&mut map_faces, l_it);
        }
        let mut seq_not_candidate: Vec<i32> = Vec::new();
        // Separates edges belonging the same face as work edge
        // for exception of edges belonging closed faces
        // OCCT L1642-1643.
        seq_not_candidate.push(work_index);
        // OCCT L1644: for (int i = 1; i <= seqIndCandidate.Length();).
        let mut i = 1usize;
        while i <= seq_ind_candidate.len() {
            let index = seq_ind_candidate[i - 1];
            let mut is_remove = false;
            // OCCT L1647-1652.
            if index == work_index {
                seq_ind_candidate.remove(i - 1);
                seq_orientations.remove(i - 1);
                is_remove = true;
            }
            if !is_remove {
                // OCCT L1655-1662.
                let mut bnd2 = sequence_sec[(index - 1) as usize].clone();
                if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd2)) {
                    bnd2 = v.clone();
                }

                if let Some((_, listfaces)) = self.my_bound_faces.get(&bat::shape_key(&bnd2)) {
                    let listfaces = listfaces.clone();
                    // OCCT L1665-1679.
                    let mut is_merged = true;
                    for l_it in &listfaces {
                        if !is_merged {
                            break;
                        }
                        if set_contains(&map_faces, l_it) {
                            let (surf_v, loc) = super::brep_tool_surface_loc(l_it);
                            if let Some(surf) = surf_v {
                                is_merged = (self
                                    .is_u_closed_surface(brep, &surf, &bnd2, l_it, loc)
                                    || self.is_v_closed_surface(brep, &surf, &bnd2, l_it, loc))
                                    && self.is_merged_closed(
                                        brep,
                                        &sequence_sec[(index - 1) as usize],
                                        &workedge,
                                        l_it,
                                    );
                            }
                        }
                    }
                    // OCCT L1680-1687.
                    if !is_merged {
                        seq_not_candidate.push(index);
                        seq_ind_candidate.remove(i - 1);
                        seq_orientations.remove(i - 1);
                        is_remove = true;
                    }
                } else {
                    // OCCT L1688-1693.
                    seq_ind_candidate.remove(i - 1);
                    seq_orientations.remove(i - 1);
                    is_remove = true;
                }
            }
            // OCCT L1695-1698: if (!isRemove) i++.
            if !is_remove {
                i += 1;
            }
        }
        // OCCT L1699-1707.
        if seq_ind_candidate.is_empty() || seq_not_candidate.len() == 1 {
            return;
        }
        if !eval_dist {
            return;
        }
        // OCCT L1708: NCollection_Array2<double> TotTabDist(1,
        // seqNotCandidate.Length(), 1, seqIndCandidate.Length()).
        let mut tot_tab_dist =
            vec![vec![0.0f64; seq_ind_candidate.len()]; seq_not_candidate.len()];
        let mut map_index: HashSet<i32> = HashSet::new();

        // Definition and removing edges which are not candidate for work edge
        // (they have other nearest edges belonging to the work face)
        // OCCT L1712-1767.
        for k in 1..=seq_not_candidate.len() {
            let index1 = seq_not_candidate[k - 1];
            let edge = sequence_sec[(index1 - 1) as usize].clone();
            let mut tmp_seq: Vec<Shape> = Vec::new();
            tmp_seq.push(edge);
            for kk in 1..=seq_ind_candidate.len() {
                tmp_seq.push(sequence_sec[(seq_ind_candidate[kk - 1] - 1) as usize].clone());
            }

            let leng_sec = tmp_seq.len();
            let mut tab_forward = vec![false; leng_sec];
            let mut tab_dist = vec![-1.0f64; leng_sec];
            let mut arr_len = vec![0.0f64; leng_sec];
            let mut tab_min_dist = vec![INFINITE_VALUE; leng_sec];
            // OCCT L1729-1731: for i1: tabDist(i1) = -1 (the Init form).
            for v in tab_dist.iter_mut() {
                *v = -1.0;
            }

            self.evaluate_distances(
                brep, &tmp_seq, &mut tab_forward, &mut tab_dist, &mut arr_len, &mut tab_min_dist,
                1,
            );
            if k == 1 {
                // OCCT L1734-1742.
                for n in 1..leng_sec {
                    if tab_dist[n] == -1.0 || tab_dist[n] > self.my_tolerance {
                        map_index.insert(n as i32);
                        continue;
                    }
                    tot_tab_dist[k - 1][n - 1] = tab_dist[n];
                    seq_forward_int(&mut tab_forward, n + 1);
                }
                // NOTE: the OCCT body appends tabForward(n + 1) into the
                // seqForward sequence; that sequence is local to this block
                // and unused after — the array form keeps the same reads.
            } else {
                // OCCT L1744-1756.
                for n in 1..leng_sec {
                    if tab_dist[n - 1] == -1.0 || tab_dist[n - 1] > self.my_tolerance {
                        continue;
                    }
                    if tab_dist[n] < tot_tab_dist[0][n - 1] {
                        map_index.insert(n as i32);
                    }
                }
            }
        }

        // OCCT L1769-1778.
        let mut i2 = seq_ind_candidate.len() as i32;
        while i2 >= 1 {
            if map_index.contains(&i2) {
                seq_ind_candidate.remove((i2 - 1) as usize);
                seq_orientations.remove((i2 - 1) as usize);
            }
            i2 -= 1;
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::FindCandidates(seqSections, mapReference,
    /// seqCandidates, seqOrientations) (cxx L1786-2088).
    pub(crate) fn find_candidates(
        &mut self,
        brep: &BRep,
        seq_sections: &mut Vec<Shape>,
        map_reference: &mut indexmap::IndexSet<i32>,
        seq_candidates: &mut Vec<i32>,
        seq_orientations: &mut Vec<bool>,
    ) -> bool {
        // OCCT L1788-1792.
        let nb_sections = seq_sections.len();
        if nb_sections <= 1 {
            return false;
        }
        // Retrieve last reference index
        // OCCT L1794.
        let ind_reference = *map_reference.last().unwrap();
        // OCCT L1795-1799.
        let mut nb_candidates = 0usize;
        let mut faces1: super::ShapeSet = super::ShapeSet::new();

        // OCCT L1801: TopoDS_Edge Edge1 = TopoDS::Edge(seqSections(indReference)).
        let edge1 = seq_sections[(ind_reference - 1) as usize].clone();

        // Retrieve faces for reference section
        // OCCT L1804-1817.
        {
            let mut bnd = edge1.clone();
            if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd)) {
                bnd = v.clone();
            }
            if let Some((_, faces)) = self.my_bound_faces.get(&bat::shape_key(&bnd)) {
                for itf1 in faces.clone() {
                    set_add_bool(&mut faces1, &itf1);
                }
            }
        }

        // Check merging conditions for candidates and remove unsatisfactory
        // OCCT L1821-1838.
        let mut seq_sections_new: Vec<Shape> = Vec::new();
        let mut seq_candidates_new: Vec<i32> = Vec::new();
        for i in 1..=nb_sections {
            if i == ind_reference as usize {
                seq_sections_new.insert(0, edge1.clone());
                seq_candidates_new.insert(0, i as i32);
            } else {
                let edge2 = seq_sections[i - 1].clone();
                seq_sections_new.push(edge2);
                seq_candidates_new.push(i as i32);
            }
        }

        // OCCT L1841-1929.
        let nb_sections_new = seq_sections_new.len();
        if nb_sections_new > 1 {
            // Evaluate distances between reference and other sections
            // OCCT L1845-1850.
            let mut arr_forward = vec![false; nb_sections_new];
            let mut arr_distance = vec![-1.0f64; nb_sections_new];
            let mut arr_len = vec![0.0f64; nb_sections_new];
            let mut arr_min_dist = vec![INFINITE_VALUE; nb_sections_new];
            self.evaluate_distances(
                brep,
                &seq_sections_new,
                &mut arr_forward,
                &mut arr_distance,
                &mut arr_len,
                &mut arr_min_dist,
                1,
            );

            // Fill sequence of candidate indices sorted by distance
            // OCCT L1853-1878.
            for i in 2..=nb_sections_new {
                let a_max_dist = arr_distance[i - 1];
                if a_max_dist >= 0.0 && a_max_dist <= self.my_tolerance && arr_len[i - 1] > self.my_min_tolerance
                {
                    // Reference section is connected to section #i
                    let mut is_inserted = false;
                    let ori = arr_forward[i - 1];
                    let mut j = 1usize;
                    while j <= seq_candidates.len() && !is_inserted {
                        // OCCT L1862: aDelta = arrDistance(i) -
                        // arrDistance(seqCandidates.Value(j)).
                        let a_delta =
                            arr_distance[i - 1] - arr_distance[(seq_candidates[j - 1] - 1) as usize];

                        if a_delta < CONFUSION {
                            if (a_delta).abs() > REAL_SMALL
                                || arr_min_dist[i - 1]
                                    < arr_min_dist[(seq_candidates[j - 1] - 1) as usize]
                            {
                                seq_candidates.insert(j - 1, i as i32);
                                seq_orientations.insert(j - 1, ori);
                                is_inserted = true;
                            }
                        }
                        j += 1;
                    }
                    if !is_inserted {
                        seq_candidates.push(i as i32);
                        seq_orientations.push(ori);
                    }
                }
            }

            // OCCT L1881-1886.
            nb_candidates = seq_candidates.len();
            if nb_candidates == 0 {
                return false; // Section has no candidates to merge
            }

            // Replace candidate indices
            // OCCT L1889-1893.
            for i in 1..=nb_candidates {
                seq_candidates[i - 1] = seq_candidates_new[(seq_candidates[i - 1] - 1) as usize];
            }
        }

        // OCCT L1896-1900.
        if nb_candidates == 0 {
            return false; // Section has no candidates to merge
        }

        // OCCT L1902-1955.
        if self.my_nonmanifold && nb_candidates > 1 {
            let mut seq_new_candidates: Vec<i32> = Vec::new();
            let mut seq_orientations_new: Vec<bool> = Vec::new();
            seq_candidates.insert(0, 1);
            seq_orientations.insert(0, true);
            let mut k = 1usize;
            while k <= seq_sections.len() && seq_candidates.len() > 1 {
                self.analysis_nearest_edges(
                    brep,
                    seq_sections,
                    seq_candidates,
                    seq_orientations,
                    k == 1,
                );
                if k == 1 && seq_candidates.is_empty() {
                    return false;
                }
                if !seq_candidates.is_empty() {
                    seq_new_candidates.push(seq_candidates[0]);
                    seq_orientations_new.push(seq_orientations[0]);
                }
                k += 1;
            }
            // OCCT L1924-1926: Prepend the new candidates.
            let mut new_candidates = seq_new_candidates;
            new_candidates.extend(seq_candidates.iter().cloned());
            *seq_candidates = new_candidates;
            let mut new_orientations = seq_orientations_new;
            new_orientations.extend(seq_orientations.iter().cloned());
            *seq_orientations = new_orientations;
            return true;
        } else {
            // Find the best approved candidate
            // OCCT L2021-2055.
            while nb_candidates != 0 {
                // Retrieve first candidate
                let ind_candidate = seq_candidates[0];
                // Candidate is approved if it is in the map
                if map_reference.contains(&ind_candidate) {
                    break;
                }
                // Find candidates for candidate #indCandidate
                map_reference.insert(ind_candidate); // Push candidate in the map
                let mut seq_candidates1: Vec<i32> = Vec::new();
                let mut seq_orientations1: Vec<bool> = Vec::new();
                let mut sections_copy = seq_sections.clone();
                let is_found = self.find_candidates(
                    brep,
                    &mut sections_copy,
                    map_reference,
                    &mut seq_candidates1,
                    &mut seq_orientations1,
                );
                map_reference.shift_remove(&ind_candidate); // Pop candidate from the map
                let mut is_found = is_found && !seq_candidates1.is_empty();
                if is_found {
                    let ind_candidate1 = seq_candidates1[0];
                    // If indReference is the best candidate for indCandidate
                    // then indCandidate is the best candidate for indReference
                    if ind_candidate1 == ind_reference {
                        break;
                    }
                    // If some other reference in the map is the best
                    // candidate for indCandidate then assume that reference
                    // is the best candidate for indReference
                    if map_reference.contains(&ind_candidate1) {
                        seq_candidates.insert(0, ind_candidate1);
                        nb_candidates += 1;
                        break;
                    }
                    is_found = false;
                }
                if !is_found {
                    // Remove candidate #1
                    seq_candidates.remove(0);
                    seq_orientations.remove(0);
                    nb_candidates -= 1;
                }
            }
        }
        // gka
        // OCCT L2058-2085.
        if nb_candidates > 0 {
            let an_ind = seq_candidates[0];
            let edge2 = seq_sections[(an_ind - 1) as usize].clone();
            let mut bnd = edge2.clone();
            if let Some(v) = self.my_section_bound.get(&bat::shape_key(&bnd)) {
                bnd = v.clone();
            }
            // gka
            if let Some((_, faces)) = self.my_bound_faces.get(&bat::shape_key(&bnd)) {
                let mut is_ok = true;
                for itf2 in faces {
                    if !is_ok {
                        break;
                    }
                    let face2 = itf2.clone();
                    // Check whether condition is satisfied
                    is_ok = !set_contains(&faces1, &face2);
                    if !is_ok {
                        is_ok = self.is_merged_closed(brep, &edge1, &edge2, &face2);
                    }
                }
                if !is_ok {
                    return false;
                }
            }
        }
        // OCCT L2087.
        nb_candidates > 0
    }

    /// OCCT BRepBuilderAPI_Sewing::MergedNearestEdges(edge, SeqMergedEdge,
    /// SeqMergedOri) (cxx L4243-4439).
    pub(crate) fn merged_nearest_edges(
        &mut self,
        brep: &BRep,
        edge: &Shape,
        seq_merged_edge: &mut Vec<Shape>,
        seq_merged_ori: &mut Vec<bool>,
    ) -> bool {
        // Retrieve edge nodes
        // OCCT L4246-4259.
        let (no1, no2) = bat::top_exp_vertices_raw(edge);
        let no1 = no1.unwrap_or_else(Shape::null);
        let no2 = no2.unwrap_or_else(Shape::null);
        let is_node1 = self.my_vertex_node.contains_key(&bat::shape_key(&no1));
        let is_node2 = self.my_vertex_node.contains_key(&bat::shape_key(&no2));
        let nno1 = node_of(brep, &self.my_vertex_node, &no1);
        let nno2 = node_of(brep, &self.my_vertex_node, &no2);
        let _ = (is_node1, is_node2);

        // Fill map of nodes connected to the node #1
        // OCCT L4262-4281.
        let mut map_vert1: super::ShapeSet = super::ShapeSet::new();
        set_add(&mut map_vert1, &nno1);
        if let Some(list) = self.my_cutting_node.get(&bat::shape_key(&nno1)).cloned() {
            for ilv in &list {
                let v1 = ilv.clone();
                set_add(&mut map_vert1, &v1);
                if !is_node1 {
                    if let Some(list2) = self.my_cutting_node.get(&bat::shape_key(&v1)).cloned() {
                        for ilvn in &list2 {
                            set_add(&mut map_vert1, ilvn);
                        }
                    }
                }
            }
        }

        // Fill map of nodes connected to the node #2
        // OCCT L4284-4303.
        let mut map_vert2: super::ShapeSet = super::ShapeSet::new();
        set_add(&mut map_vert2, &nno2);
        if let Some(list) = self.my_cutting_node.get(&bat::shape_key(&nno2)).cloned() {
            for ilv in &list {
                let v1 = ilv.clone();
                set_add(&mut map_vert2, &v1);
                if !is_node2 {
                    if let Some(list2) = self.my_cutting_node.get(&bat::shape_key(&v1)).cloned() {
                        for ilvn in &list2 {
                            set_add(&mut map_vert2, ilvn);
                        }
                    }
                }
            }
        }

        // Find all possible contiguous edges
        // OCCT L4306-4360.
        let mut seq_edges: Vec<Shape> = Vec::new();
        seq_edges.push(edge.clone());
        let mut map_edges: super::ShapeSet = super::ShapeSet::new();
        set_add(&mut map_edges, edge);
        for i in 0..map_vert1.len() {
            let (_, node1) = map_vert1.get_index(i).unwrap().clone();
            let node_sections = match self.my_node_sections.get(&bat::shape_key(&node1)) {
                Some(l) => l.clone(),
                None => continue,
            };
            for ilsec in &node_sections {
                let sec = ilsec.clone();
                if sec.is_same(&edge) {
                    continue;
                }
                // Retrieve section nodes
                // OCCT L4321-4329.
                let (vs1, vs2) = bat::top_exp_vertices_raw(&sec);
                let vs1 = vs1.unwrap_or_else(Shape::null);
                let vs2 = vs2.unwrap_or_else(Shape::null);
                let vs1n = node_of(brep, &self.my_vertex_node, &vs1);
                let vs2n = node_of(brep, &self.my_vertex_node, &vs2);
                if (set_contains(&map_vert1, &vs1n) && set_contains(&map_vert2, &vs2n))
                    || (set_contains(&map_vert1, &vs2n) && set_contains(&map_vert2, &vs1n))
                {
                    if set_add_bool(&mut map_edges, &sec) {
                        // Check for rejected cutting
                        // OCCT L4335-4352.
                        let mut is_rejected = set_contains(&self.my_merged_edges, &sec);
                        if !is_rejected {
                            if let Some(sections) = self.my_bound_sections.get(&bat::shape_key(&sec)).cloned() {
                                for section in &sections {
                                    if !is_rejected && set_contains(&self.my_merged_edges, section) {
                                        is_rejected = true;
                                    }
                                    if is_rejected {
                                        break;
                                    }
                                }
                            }
                        }
                        if !is_rejected {
                            if let Some(bnd_v) = self.my_section_bound.get(&bat::shape_key(&sec)) {
                                let bnd = bnd_v.clone();
                                is_rejected = !self.my_bound_sections.contains_key(&bat::shape_key(&bnd))
                                    || set_contains(&self.my_merged_edges, &bnd);
                            }
                        }

                        if !is_rejected {
                            seq_edges.push(sec);
                        }
                    }
                }
            }
        }

        // OCCT L4362: mapEdges.Clear().
        map_edges.clear();

        // OCCT L4364-4364.
        let mut success = false;

        // OCCT L4366: int nbSection = seqEdges.Length().
        let nb_section = seq_edges.len();
        if nb_section > 1 {
            // Find the longest edge CCI60011
            // OCCT L4368-4386.
            let mut ind_ref = 1usize;
            if self.my_nonmanifold {
                let mut len_ref = 0.0f64;
                for i in 1..=nb_section {
                    let len = match bat::brep_tool_curve(&seq_edges[i - 1]) {
                        Some((c, f, l)) => gcpnts_abscissa_length(&c, f, l),
                        None => 0.0,
                    };
                    if len > len_ref {
                        ind_ref = i;
                        len_ref = len;
                    }
                }
                if ind_ref != 1 {
                    let long_edge = seq_edges[ind_ref - 1].clone();
                    seq_edges[ind_ref - 1] = seq_edges[0].clone();
                    seq_edges[0] = long_edge;
                }
            }

            // Find merging candidates
            // OCCT L4389-4392.
            let mut seq_forward: Vec<bool> = Vec::new();
            let mut seq_candidates: Vec<i32> = Vec::new();
            let mut map_reference = indexmap::IndexSet::new();
            map_reference.insert(ind_ref as i32); // Add index of reference section
            if self.find_candidates(
                brep,
                &mut seq_edges,
                &mut map_reference,
                &mut seq_candidates,
                &mut seq_forward,
            ) {
                let nb_candidates = seq_candidates.len();
                // Record candidate sections
                // OCCT L4396-4409.
                for i in 1..=nb_candidates {
                    // Retrieve merged edge
                    let iedge = seq_edges[(seq_candidates[i - 1] - 1) as usize].clone();
                    let ori = seq_forward[i - 1];
                    seq_merged_edge.push(iedge);
                    seq_merged_ori.push(ori);
                    if !self.my_nonmanifold {
                        break;
                    }
                }
                success = nb_candidates != 0;
            }
        }

        // OCCT L4438.
        success
    }
}

/// OCCT `gp::Resolution()`-independent helper: the RealSmall() constant
/// (Standard_Real.hxx) — the smallest positive representable real.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

/// OCCT `myVertexNode.FindFromKey(v)` — the node of a vertex (the identity
/// when the vertex is not bound).
fn node_of(_brep: &BRep, map: &super::IdxShapeMap, v: &Shape) -> Shape {
    match map.get(&bat::shape_key(v)) {
        Some((_, node)) => node.clone(),
        None => v.clone(),
    }
}

/// OCCT `seqForward.Append(tabForward(n + 1) ? 1 : 0)` — the OCCT
/// NCollection_Sequence<int> append of the boolean flag (the sequence value
/// is dead after the loop; kept for statement parity).
fn seq_forward_int(_tab_forward: &[bool], _n: usize) {}

/// OCCT `BndLib_Add2dCurve::Add(aC2d, f, l, tol, B)` (BndLib_Add2dCurve.cxx
/// L29-36) — the 2d bounding box of the pcurve subrange; the rcad re-host
/// samples the pcurve over the range into `b` (the same documented
/// sampling reduction as the 3d form).
pub(crate) fn bnd_lib_add2d_curve(
    c2d: &rcad_kernel::geom::Curve2d,
    f: f64,
    l: f64,
    tol: f64,
    b: &mut BndBox2d,
) {
    let n = 16usize;
    for i in 0..=n {
        let t = f + (l - f) * (i as f64) / (n as f64);
        let p = c2d.point_at(t);
        b.update_xy(p.x, p.y);
    }
    b.enlarge(tol);
}
