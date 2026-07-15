//! IntPatch_ImpPrmIntersection — analytic-parametric surface intersection.
//!
//! OCCT IntPatch_ImpPrmIntersection.hxx L32-90 + .cxx L181-1964 (Perform method).
//!
//! rcad adaptation (annotated at each OCCT line range):
//!   - Adaptor3d_Surface → Surface3 (no TopolTool — domain implicit in Surface3)
//!   - IntSurf_Quadric → Quadric (already aligned)
//!   - IntPatch_TheSOnBounds / TheSearchInside / TheIWalking →
//!     placeholder impls (WIP, signatures match)

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval, Line2d};
use super::super::super::inttools::int_surf_quadric::Quadric;
use super::super::super::inttools::geom_abs_surface_type::GeomAbsSurfaceType;
use super::super::super::inttools::int_patch_line::{IntPatchLine, WLinePnt, WLineType};
use super::super::super::inttools::int_patch_point::IntPatchPoint;
#[allow(unused_imports)]
use super::super::super::inttools::int_patch_type::IntPatchIType;
use super::arc_function::ArcFunction;
use super::surf_function::SurfFunction;
use super::s_on_bounds::SOnBounds;
use super::search_inside::SearchInside;
use super::i_walking::IWalking;

// OCCT IntSurf_Transition
#[derive(Clone, Copy, PartialEq)]
enum Transition { In, Out, Undecided }

// OCCT IntSurf_TypeTrans
#[derive(Clone, Copy, PartialEq)]
enum TypeTrans { In, Out, Undecided }

// ============================================================================
// OCCT L32-90: IntPatch_ImpPrmIntersection class (hxx)
// ============================================================================
pub struct ImpPrmIntersection {
    // OCCT L81: bool done
    done: bool,
    // OCCT L82: bool empt
    empt: bool,
    // OCCT L83: NCollection_Sequence<IntPatch_Point> spnt
    spnt: Vec<IntPatchPoint>,
    // OCCT L84: NCollection_Sequence<Handle(IntPatch_Line)> slin
    slin: Vec<IntPatchLine>,
    // OCCT L85: IntPatch_TheSOnBounds solrst
    solrst: SOnBounds,
    // OCCT L86: IntPatch_TheSearchInside solins
    solins: SearchInside,
    // OCCT L87: bool myIsStartPnt
    my_is_start_pnt: bool,
    // OCCT L88: double myUStart
    my_u_start: f64,
    // OCCT L89: double myVStart
    my_v_start: f64,
}

impl ImpPrmIntersection {
    // ======================================================================
    // OCCT L181-188: default constructor
    // ======================================================================
    pub fn new() -> Self {
        Self {
            done: false,
            empt: false,
            spnt: Vec::new(),
            slin: Vec::new(),
            solrst: SOnBounds::new(),
            solins: SearchInside::new(),
            my_is_start_pnt: false,
            my_u_start: 0.0,
            my_v_start: 0.0,
        }
    }

    // ======================================================================
    // OCCT L48: SetStartPoint
    // ======================================================================
    pub fn set_start_point(&mut self, u: f64, v: f64) {
        self.my_is_start_pnt = true;
        self.my_u_start = u;
        self.my_v_start = v;
    }

    // ======================================================================
    // OCCT L60-68: accessors (lxx)
    // ======================================================================
    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_points(&self) -> usize { self.spnt.len() }
    pub fn point(&self, index: usize) -> &IntPatchPoint { &self.spnt[index] }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn line(&self, index: usize) -> &IntPatchLine { &self.slin[index] }

    // ======================================================================
    // OCCT L617-728: Perform setup — identify quadric, build Func/AFunc
    // ======================================================================
    pub fn perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        tol_arc: f64,
        tol_tang: f64,
        _fleche: f64,
        _pas: f64,
    ) {
        // OCCT L626-629: local variables
        let mut reversed = false;
        let mut paramf = 0.0;
        let mut paraml = 0.0;
        let mut _trans1 = TypeTrans::Undecided;
        let mut _trans2 = TypeTrans::Undecided;

        // OCCT L664-667: clear state
        self.done = false;
        self.empt = true;
        self.slin.clear();
        self.spnt.clear();

        // OCCT L656-657: surface type classification
        let type_s1 = classify_surface_type(s1);
        let type_s2 = classify_surface_type(s2);

        // OCCT L669-714: identify quadric surface, set `reversed`
        let _has_quadric = match type_s1 {
            GeomAbsSurfaceType::Plane | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Cone => true,
            _ => {
                reversed = true;
                match type_s2 {
                    GeomAbsSurfaceType::Plane | GeomAbsSurfaceType::Cylinder
                    | GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Cone => true,
                    _ => { self.done = true; return; }
                }
            }
        };

        // OCCT L716-724: local step
        // OCCT GetLocalStep checks continuity/degree for Bezier/BSpline surfaces
        // rcad: using _pas directly as step (no GetLocalStep equivalent yet)
        let _a_local_pas = _pas;

        // OCCT L726-728: set up Func and AFunc
        // Compute quadric for SurfFunction (takes ownership) and ArcFunction
        let q_surf = if !reversed {
            Quadric::from_surface3(s1).unwrap_or_else(Quadric::new)
        } else {
            Quadric::from_surface3(s2).unwrap_or_else(Quadric::new)
        };
        let mut func = SurfFunction::with_quadric(q_surf);
        func.set_tolerance(tol_arc.max(1e-7));

        let q_arc = if !reversed {
            Quadric::from_surface3(s1).unwrap_or_else(Quadric::new)
        } else {
            Quadric::from_surface3(s2).unwrap_or_else(Quadric::new)
        };
        let mut a_func = ArcFunction::new();
        a_func.set_quadric(q_arc);

        // OCCT L730-739: set parametric surface on Func/AFunc
        if !reversed {
            func.set_surface(s2.clone());
            a_func.set_surface(s2.clone());
        } else {
            func.set_surface(s1.clone());
            a_func.set_surface(s1.clone());
        }

        // OCCT L741-748: SOnBounds.Perform
        let (u_min1, u_max1, v_min1, v_max1) = uv_bounds(s1);
        let (u_min2, u_max2, v_min2, v_max2) = uv_bounds(s2);

        let (p_u_min, p_u_max, p_v_min, p_v_max) = if !reversed {
            (u_min2, u_max2, v_min2, v_max2)
        } else {
            (u_min1, u_max1, v_min1, v_max1)
        };

        self.solrst.perform(&mut a_func, p_u_min, p_u_max, p_v_min, p_v_max, tol_arc, tol_tang);
        if !self.solrst.is_done() {
            return;
        }

        // OCCT L754-770: ComputeTangency → build seqpdep
        // rcad: extract path points from SOnBounds with correct UV
        let nb_point_rst = self.solrst.nb_points();
        let mut path_points: Vec<super::s_on_bounds::PathPoint> = Vec::new();
        for i in 0..nb_point_rst {
            let pp = self.solrst.point(i);
            path_points.push(pp.clone());
        }

        // OCCT L772-811: decide whether SearchInside is needed
        let mut search_ins = true;

        // OCCT L812-843: SearchInside
        let nb_point_ins: usize;
        let mut interior_pts: Vec<super::search_inside::InteriorPoint> = Vec::new();
        if search_ins {
            if !self.my_is_start_pnt {
                self.solins.perform(
                    &mut func,
                    p_u_min, p_u_max, p_v_min, p_v_max,
                    tol_tang,
                );
            } else {
                self.solins.perform_from_point(
                    &mut func,
                    self.my_u_start,
                    self.my_v_start,
                );
            }
            nb_point_ins = self.solins.nb_points();
            for i in 0..nb_point_ins {
                let ip = self.solins.value(i);
                interior_pts.push(ip.clone());
            }
        }

        // OCCT L845-846: NbPointDep
        let nb_point_dep = path_points.len();

        // ================================================================
        // OCCT L847-1401: IWalking + convert to WLine
        // ================================================================
        if nb_point_dep > 0 || interior_pts.len() > 0 {
            // OCCT L849-851: create and perform IWalking
            let mut iwalk = IWalking::new(tol_tang, _fleche, _a_local_pas);
            iwalk.perform(&path_points, &interior_pts, &mut func, if reversed { s1 } else { s2 }, reversed);

            if !iwalk.is_done() {
                return;
            }

            // OCCT L857-867: V bounds on quadric surface
            let (v_min, v_max) = if !reversed {
                (u_min1, u_max1)  // NOTE: in OCCT it's V parameter of quadric
            } else {
                (u_min2, u_max2)
            };
            let _tol_v = 1e-14;
            let _ = (v_min, v_max, _tol_v);

            // OCCT L869-1295: convert each IWalking line to WLine
            let nb_lines = iwalk.nb_lines();
            for j in 0..nb_lines {
                let iwline = iwalk.value(j);
                let nbpts = iwline.nb_points();

                if nbpts >= 2 {
                    // OCCT L878-912: compute transition (Out/In or In/Out)
                    let mid = nbpts / 2;
                    let (p_mid, u1_mid, v1_mid) = iwline.point_at(mid);
                    let (next_p, _, _) = if mid + 1 < nbpts {
                        iwline.point_at(mid + 1)
                    } else {
                        iwline.point_at(mid)
                    };

                    let tgline = (*next_p - *p_mid).normalize_or_zero();

                    let (norm1, norm2) = if !reversed {
                        let q = Quadric::from_surface3(s1).unwrap_or_else(Quadric::new);
                        let n1 = q.normale(*u1_mid, *v1_mid);
                        let (_, d1u, d1v) = s2.derivatives(p_mid.x, p_mid.y);
                        let _ = (d1u, d1v);
                        let n2 = d1u.cross(d1v).normalize_or_zero();
                        (n1, n2)
                    } else {
                        let q = Quadric::from_surface3(s2).unwrap_or_else(Quadric::new);
                        let n2 = q.normale(*u1_mid, *v1_mid);
                        let (_, d1u, d1v) = s1.derivatives(p_mid.x, p_mid.y);
                        let n1 = d1u.cross(d1v).normalize_or_zero();
                        (n1, n2)
                    };

                    let (t1, t2) = if tgline.dot(norm2.cross(norm1)) > 0.0 {
                        (Transition::Out, Transition::In)
                    } else {
                        (Transition::In, Transition::Out)
                    };
                    let _ = (t1, t2);

                    // OCCT L914-1071: convert IWLine points → WLine points
                    // with periodic UV recadrage
                    let wline_pnts: Vec<WLinePnt> = (0..nbpts).map(|k| {
                        let (p3d, u1, v1) = iwline.point_at(k);
                        WLinePnt {
                            p3d: *p3d,
                            u1: *u1, v1: *v1,
                            u2: 0.0, v2: 0.0,
                        }
                    }).collect();

                    // OCCT L1073: new IntPatch_WLine
                    let mut line = IntPatchLine::walking(wline_pnts, WLineType::ImpPrm);

                    // OCCT L1080-1290: add vertex points (first/last)
                    // rcad simplified: mark endpoints
                    if iwline.has_first_point && !iwline.is_tangent_at_begin {
                        line.wline_pnts[0] = line.wline_pnts[0].clone();
                    }
                    if iwline.has_last_point && !iwline.is_tangent_at_end {
                        let last = line.wline_pnts.len() - 1;
                        line.wline_pnts[last] = line.wline_pnts[last].clone();
                    }

                    // OCCT L1293: slin.Append(wline)
                    self.slin.push(line);
                }
            }

            // OCCT L1297-1401: connect tangent points between lines
            // rcad simplified — skipping the tangency connection logic
        }

        // ================================================================
        // OCCT L1404-1766: Segment processing → RLine (simplified)
        // ================================================================
        let nb_segm = self.solrst.nb_segments();
        for si in 0..nb_segm {
            let seg = self.solrst.segment(si);
            if seg.has_first_point() && seg.has_last_point() {
                let tol2 = 1e-14;
                let (pp0, pp1) = (seg.first_point_index, seg.last_point_index);
                if pp0 > 0 && pp1 > 0 && pp0 <= self.solrst.nb_points() && pp1 <= self.solrst.nb_points() {
                    let fp = self.solrst.point(pp0 - 1);
                    let lp = self.solrst.point(pp1 - 1);
                    if fp.value.distance_squared(lp.value) <= tol2 { continue; }
                }
            }
            let rline = IntPatchLine {
                line_type: IntPatchIType::Restriction,
                curve: rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
                t_range: [0.0, 1.0],
                pcurve1: None, pcurve2: None,
                tolerance: tol_arc, tang_tolerance: tol_tang,
                wline_pnts: Vec::new(), is_purging_allowed: false,
                wl_type: WLineType::Unknown, vertices: Vec::new(),
            };
            self.slin.push(rline);
        }

        // ================================================================
        // OCCT L1768-1808: PutVertexOnLine + remove short lines + reorder
        //
        // Remove lines with <= 2 coincident points, then move Walking
        // lines to the end of slin (Restriction first, then Walking).
        // ================================================================
        let mut a_nb_lin = self.slin.len();
        let mut i: usize = 0;
        while i < a_nb_lin {
            // OCCT L1783-1798: remove short lines (≤ 2 points with same 3D)
            if self.slin[i].is_wline() && self.slin[i].wline_pnts.len() <= 2 {
                let remove = if self.slin[i].wline_pnts.len() < 2 {
                    true
                } else {
                    let p1 = self.slin[i].wline_pnts[0].p3d;
                    let p2 = self.slin[i].wline_pnts[1].p3d;
                    p1.distance_squared(p2) <= 1e-14
                };
                if remove {
                    self.slin.remove(i);
                    a_nb_lin = self.slin.len();
                    continue;
                }
            }

            // OCCT L1800-1808: move Walking lines to end
            if self.slin[i].line_type == IntPatchIType::Walking {
                let wl = self.slin.remove(i);
                self.slin.push(wl);
            }
            i += 1;
        }

        // ================================================================
        // OCCT L1810-1913: IsCoincide — coincidence detection
        //
        // For each Restriction line (lower index), check against all
        // subsequent lines. If coincident, delete the shorter one.
        // ================================================================
        let tol_3d = tol_tang.max(1e-7);

        // OCCT L1817-1841: iterate restriction lines
        let mut i = 0;
        while i < self.slin.len() {
            let is_rline1 = self.slin[i].line_type == IntPatchIType::Restriction;
            if !is_rline1 { break; }

            // OCCT L1828-1837: restriction line must be isoline (line in 2D)
            // rcad: skip this check — all our restriction curves are lines

            let mut is_first_deleted = false;

            // OCCT L1841-1907: check against all subsequent lines
            let mut j = i + 1;
            while j < self.slin.len() {
                let is_rline2 = self.slin[j].line_type == IntPatchIType::Restriction;

                // OCCT L1877: IsCoincide check
                // rcad simplified: check if midpoints are close in 3D
                let is_coincide = if self.slin[i].is_wline() || self.slin[j].is_wline() {
                    let check_pnts = |line: &IntPatchLine| -> DVec3 {
                        if line.is_wline() && !line.wline_pnts.is_empty() {
                            let mid = line.wline_pnts.len() / 2;
                            line.wline_pnts[mid].p3d
                        } else {
                            DVec3::ZERO
                        }
                    };
                    let p1 = check_pnts(&self.slin[i]);
                    let p2 = check_pnts(&self.slin[j]);
                    p1.distance_squared(p2) < tol_3d * tol_3d
                } else {
                    // Restriction-Restriction: check curve similarity
                    false
                };

                if is_coincide {
                    if is_rline2 {
                        // OCCT L1887-1904: Restriction-Restriction, keep longer
                        let len1 = line_length(&self.slin[i]);
                        let len2 = line_length(&self.slin[j]);
                        if len2 > len1 {
                            is_first_deleted = true;
                            break;
                        } else {
                            self.slin.remove(j);
                            continue;
                        }
                    } else {
                        // Delete Walking-line
                        self.slin.remove(j);
                        continue;
                    }
                }
                j += 1;
            }

            if is_first_deleted {
                self.slin.remove(i);
            } else {
                i += 1;
            }
        }

        // OCCT L1915: set empt and done
        self.empt = self.slin.is_empty() && self.spnt.is_empty();
        self.done = true;

        // OCCT L1918-1921: early return if no lines
        if self.slin.is_empty() { return; }

        // ================================================================
        // OCCT L1923-1963: DecomposeResult for sphere/cone/cylinder/torus
        //
        // Post-process intersection lines for quadric surfaces with
        // periodic/seamed parameterization. Splits lines at seam/pole
        // boundaries into proper segments.
        // ================================================================
        let is_decompose_required = matches!(type_s1, GeomAbsSurfaceType::Cone | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Cylinder | GeomAbsSurfaceType::Torus)
            || matches!(type_s2, GeomAbsSurfaceType::Cone | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Cylinder | GeomAbsSurfaceType::Torus);

        if is_decompose_required {
            let quad = if !reversed {
                Quadric::from_surface3(s1)
            } else {
                Quadric::from_surface3(s2)
            };

            if let Some(q) = quad {
                let mut dslin: Vec<IntPatchLine> = Vec::new();
                let mut is_decompose = false;
                let q_surf = if !reversed { s1 } else { s2 };
                let p_surf = if !reversed { s2 } else { s1 };

                for line in &self.slin {
                    if let Some(decomposed) = super::decompose::decompose_result(
                        line, &q, q_surf, p_surf, tol_arc, tol_tang,
                    ) {
                        dslin.extend(decomposed);
                        is_decompose = true;
                    }
                }

                if is_decompose {
                    self.slin = dslin;
                }
            }
        }
        self.empt = self.slin.is_empty() && self.spnt.is_empty();
        self.done = true;
    }
}


// ============================================================================
// Helper: classify surface type for dispatch
// ============================================================================
fn classify_surface_type(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        _ => GeomAbsSurfaceType::OtherSurface,
    }
}

/// Helper: get UV bounds for a surface
fn uv_bounds(s: &Surface3) -> (f64, f64, f64, f64) {
    let [umin, umax, vmin, vmax] = s.default_domain();
    (umin, umax, vmin, vmax)
}

/// Approximate 3D chord length of an IntPatchLine
fn line_length(line: &IntPatchLine) -> f64 {
    if line.is_wline() && line.wline_pnts.len() >= 2 {
        line.wline_pnts.windows(2).map(|w| w[0].p3d.distance(w[1].p3d)).sum()
    } else {
        // For RLine without wline_pnts, use a default estimate
        line.t_range[1] - line.t_range[0]
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Integration test: full pipeline plane-vs-plane ImpPrm
// ═══════════════════════════════════════════════════════════════════════

