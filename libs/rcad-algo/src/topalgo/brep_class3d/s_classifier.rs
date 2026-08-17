// OCCT BRepClass3d_SClassifier (BRepClass3d_SClassifier.hxx / .cxx)
// Base class for solid classification.
//
// OCCT BRepClass3d_SClassifier::Perform (SClassifier.cxx L203-490) is the
// core point-in-solid algorithm: it builds a line toward a point on a face of
// the solid (SolidExplorer::Segment/OtherSegment), detects line-edge/vertex
// interferences (BndBoxTreeSelectorLine) that make the line unusable
// (faulty-line retry), intersects the line with every face
// (IntCurvesFace_Intersector) and decides IN/OUT from the transition of the
// closest intersection (Trans). The point near a vertex/edge is ON
// (BndBoxTreeSelectorPoint). See temp/sclassifier_alignment.md.

use crate::topalgo::brep_class3d::bnd_box_tree::{BndBoxTreeSelectorLine, BndBoxTreeSelectorPoint};
use crate::topalgo::brep_class3d::solid_explorer::SolidExplorer;
use glam::DVec3;

/// OCCT BRepClass3d_SClassifier — base classification algorithm.
pub struct SClassifier {
    pub my_state: u8,  // 0=unknown, 1=faulty, 2=ON, 3=IN, 4=OUT
}

impl SClassifier {
    pub fn new() -> Self {
        SClassifier { my_state: 0 }
    }

    /// OCCT L203-490: Perform(SolidExplorer, P, Tol) — 1:1 translation.
    pub fn perform(&mut self, explorer: &mut SolidExplorer, p: DVec3, tol: f64) {
        // OCCT L207-212: Reject — a solid without faces is the whole space.
        if explorer.reject(p) {
            self.my_state = 3; // IN
            return;
        }

        // OCCT L214-230: point near a vertex/edge -> ON.
        {
            let (edges, e_tols, e_ranges, verts, v_tols) = explorer.map_ev();
            let mut sel = BndBoxTreeSelectorPoint::new(edges, e_tols, e_ranges, verts, v_tols);
            sel.set_current_point(p);
            if sel.select() > 0 {
                self.my_state = 2; // ON
                return;
            }
        }

        // OCCT L232: mapEF = TopExp::MapShapesAndAncestors(Solid, EDGE, FACE).
        let map_ef = explorer.map_ef();

        // OCCT L240-490: the faulty-line retry loop.
        let mut is_faulty_line = true;
        let mut an_ind_face = 0usize;
        let mut parmin = 0.0f64;
        while is_faulty_line {
            // OCCT L257-260: Segment / OtherSegment.
            let (i_flag, l_origin, l_dir, par) = if an_ind_face == 0 {
                explorer.segment(p, false)
            } else {
                explorer.segment(p, true)
            };
            let a_cur_ind = explorer.face_segment_index();
            if a_cur_ind > an_ind_face {
                an_ind_face = a_cur_ind;
            } else {
                self.my_state = 1; // Faulty
                return;
            }
            // OCCT L274-284: iFlag handling.
            if i_flag == 1 {
                self.my_state = 2; // ON
                return;
            }
            if i_flag == 2 {
                self.my_state = 4; // OUT
                return;
            }
            if i_flag == 3 {
                continue; // point on surface but outside the face
            }

            is_faulty_line = false;
            parmin = f64::MAX;

            // OCCT L288-372: line vs vertices/edges interference.
            let mut near_fault_par = f64::MAX;
            let (edges, e_tols, e_ranges, verts, v_tols) = explorer.map_ev();
            let mut a_sel_line = BndBoxTreeSelectorLine::new(edges, e_tols, e_ranges, verts, v_tols);
            a_sel_line.clear_results();
            a_sel_line.set_current_line(l_origin, l_dir, par);
            let sels_evl = a_sel_line.select();
            let mut lv_ints: Vec<(u64, u32)> = Vec::new();
            if sels_evl > 0 {
                // Line and vertices.
                for &(v_idx, lp) in a_sel_line.vert_params() {
                    // OCCT L299-305: LVInts.Add(V); NearFaultPar = min |LP|.
                    let key = explorer.vertex_key(v_idx);
                    lv_ints.push(key);
                    if lp.abs() < near_fault_par.abs() {
                        near_fault_par = lp;
                    }
                }
                // Line and edges.
                for &(e_idx, param, lpar) in a_sel_line.edge_params() {
                    // OCCT L307-319: ffs = mapEF.FindFromKey(EE); must be 2.
                    let ekey = explorer.edge_key(e_idx);
                    let ffs = match map_ef.get(&ekey) {
                        Some(v) => v.clone(),
                        None => continue,
                    };
                    if ffs.len() != 2 {
                        continue;
                    }
                    let (v1, v2) = explorer.edge_vertices(e_idx);
                    if lv_ints.contains(&v1) || lv_ints.contains(&v2) {
                        continue;
                    }
                    // OCCT L320-327: GetTransi(f1, f2, EE, param, L, tran).
                    let mut tran = 0u8; // 0=Tangent, 1=In, 2=Out
                    let tst = get_transi(explorer, ffs[0], ffs[1], e_idx, param, l_dir, &mut tran);
                    if tst == 1 && lpar.abs() < parmin.abs() {
                        parmin = lpar;
                        trans(parmin, &mut tran, &mut self.my_state);
                    } else if lpar.abs() < near_fault_par.abs() {
                        near_fault_par = lpar;
                    }
                }
            }

            // OCCT L372-460: intersect the line with every face.
            // OCCT iterates shells (InitShell/MoreShell/NextShell) and faces
            // per shell (InitFace/MoreFace/NextFace); RejectShell/RejectFace
            // are always false (SolidExplorer.cxx L1037-1093), so the nested
            // loops reduce to one pass over all faces of the solid. rcad's
            // explorer stores the faces flat, so a single loop over them.
            let n_faces = explorer.nb_faces();
            for fi in 0..n_faces {
                // OCCT L377-397: prolong the segment — the intersector may not
                // find intersection points with the original range due to rough
                // triangulation of a parameterized surface. addW is extended by
                // the face's bounding box (GetAddToParam); minW uses the
                // UNEXTENDED AddW, maxW the extended one.
                let add_w0 = (10.0 * tol).max(0.01 * par);
                let mut add_w = add_w0;
                // OCCT L383-393: the box must be finite (not void, not whole).
                if let Some(bb) = explorer.face_bounding_box(fi) {
                    if bb.0.is_finite() && bb.1.is_finite() {
                        let box_add_w = explorer.get_add_to_param(l_origin, l_dir, par, bb);
                        add_w = add_w.max(box_add_w);
                    }
                }
                let min_w = -add_w0;
                let max_w = (par * 10.0).min(par + add_w);
                let pts = explorer.face_line_intersections(fi, l_origin, l_dir, min_w, max_w);
                if let Some((is_parallel, points)) = pts {
                    if points.is_empty() {
                        // OCCT L380-404: the parallel case — only when the
                        // intersector reports IsParallel().
                        if is_parallel {
                            // OCCT L405-442: check the distance between the
                            // surface and the point (Extrema_ExtPS with
                            // PConfusion and the MIN flag); the nearest
                            // extremum within Tol gives ON when the UV is IN/ON
                            // the face domain (ClassifyUVPoint, tol 1e-7).
                            let dist2 = explorer.point_face_distance(p, fi);
                            if dist2.is_some() {
                                let (d2, uv) = dist2.unwrap();
                                if d2 <= tol * tol {
                                    let st = explorer.classify_uv_point_at(fi, uv);
                                    if st == 1 || st == 2 {
                                        self.my_state = 2; // ON
                                        parmin = 0.0;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    for (w, state, tran, u, v) in points {
                        if w.abs() < parmin.abs() - CONFUSION_EPS {
                            parmin = w;
                            if parmin.abs() <= tol {
                                self.my_state = 2; // ON
                                break;
                            } else if state == 1 {
                                // IN — use the transition.
                                if tran == 0 {
                                    continue; // tangent — ignore
                                }
                                trans(parmin, &mut tran.clone(), &mut self.my_state);
                            } else if state == 2 {
                                // ON — the line is faulty.
                                is_faulty_line = true;
                                break;
                            }
                        }
                    }
                    if self.my_state == 2 {
                        break;
                    }
                } else {
                    self.my_state = 1; // Faulty
                }
            }
            if self.my_state == 2 {
                break;
            }
            // OCCT L462-469: NearFaultPar check.
            if near_fault_par != f64::MAX && parmin.abs() >= near_fault_par.abs() - CONFUSION_EPS {
                is_faulty_line = true;
            }
        }
    }
}

/// Precision::PConfusion() — the parameter comparison tolerance.
const CONFUSION_EPS: f64 = 1e-7;

/// OCCT Trans (SClassifier.cxx L728-765) — state from the closest
/// intersection transition: a line going Out of the solid means the point is
/// IN; In means OUT. A negative parmin reverses the transition.
fn trans(parmin: f64, tran: &mut u8, state: &mut u8) {
    if parmin < 0.0 {
        // 1 = In, 2 = Out.
        *tran = if *tran == 2 { 1 } else { 2 };
    }
    if *tran == 2 {
        *state = 3; // IN — the line is going from inside to outside.
    } else {
        *state = 4; // OUT
    }
}

/// OCCT GetTransi (SClassifier.cxx L654-728) — the transition of the line L
/// crossing the edge shared by faces f1/f2 at the given edge parameter.
/// Returns 1 = OK, 0 = skip, -1 = probably a faulty line.
fn get_transi(
    explorer: &SolidExplorer,
    f1: usize,
    f2: usize,
    _e_idx: usize,
    param: f64,
    l_dir: DVec3,
    trans: &mut u8,
) -> i32 {
    let nf1 = match explorer.face_bound_normal(f1, param) {
        Some(n) => n,
        None => return -1,
    };
    let nf2 = match explorer.face_bound_normal(f2, param) {
        Some(n) => n,
        None => return -1,
    };
    let l_dir = l_dir.normalize_or_zero();
    if l_dir.dot(nf1).abs() < ANGULAR || l_dir.dot(nf2).abs() < ANGULAR {
        // The line is orthogonal to the normals -> tangent.
        return -1;
    }
    if nf1.dot(nf2).abs() > 1.0 - ANGULAR && nf1.cross(nf2).length() < ANGULAR {
        // Parallel normals.
        let ang_d = nf1.dot(l_dir);
        if ang_d.abs() < ANGULAR {
            return -1;
        } else if ang_d > 0.0 {
            *trans = 2; // Out
        } else {
            *trans = 1; // In
        }
        return 1;
    }
    // OCCT L688-707: project LDir on the plane of nf1/nf2.
    let n = nf1.cross(nf2);
    let proj_l = n.cross(l_dir).cross(n);
    let f_ad = nf1.dot(proj_l);
    let s_ad = nf2.dot(proj_l);
    if f_ad < -ANGULAR && s_ad < -ANGULAR {
        *trans = 1; // In
    } else if f_ad > ANGULAR && s_ad > ANGULAR {
        *trans = 2; // Out
    } else {
        return 0;
    }
    1
}

/// Precision::Angular().
const ANGULAR: f64 = 1e-12;
